using CUDA
using LinearAlgebra
using Random
using Printf
using Base.Threads

# =====================================================================
# 1. STRICT FEASIBILITY VALIDATOR
# =====================================================================
function is_valid_permutation(p::Vector{Int32}, N::Int)::Bool
    if length(p) != N; return false; end
    seen = falses(N)
    @inbounds for i in 1:N
        val = p[i]
        if val < 1 || val > N || seen[val]; return false; end
        seen[val] = true
    end
    return true
end

# =====================================================================
# 2. EXACT INT64 MINSPM SPARSE GALOIS EVALUATION KERNEL
# =====================================================================
function minspm_tai256c_int64_kernel!(F_sparse_i, F_sparse_j, F_vals, D, P_batch, energies, nnz_terms)
    tid = threadIdx().x
    warp_id = div(tid - 1, 32) + 1
    lane_id = rem(tid - 1, 32) + 1
    
    batch_id = (blockIdx().x - 1) * div(blockDim().x, 32) + warp_id
    shared_vals = @cuDynamicSharedMem(Int64, blockDim().x)
    val = Int64(0)
    
    if batch_id <= size(P_batch, 2)
        idx = lane_id
        while idx <= nnz_terms
            i = F_sparse_i[idx]
            j = F_sparse_j[idx]
            f_val = F_vals[idx]
            p_i = P_batch[i, batch_id]
            p_j = P_batch[j, batch_id]
            
            if p_i >= 1 && p_i <= 256 && p_j >= 1 && p_j <= 256
                val += f_val * D[p_i, p_j]
            end
            idx += 32
        end
    end
    
    shared_vals[tid] = val
    sync_threads()
    
    if lane_id <= 16; shared_vals[tid] += shared_vals[tid + 16]; end; sync_threads()
    if lane_id <= 8;  shared_vals[tid] += shared_vals[tid + 8];  end; sync_threads()
    if lane_id <= 4;  shared_vals[tid] += shared_vals[tid + 4];  end; sync_threads()
    if lane_id <= 2;  shared_vals[tid] += shared_vals[tid + 2];  end; sync_threads()
    if lane_id == 1;  shared_vals[tid] += shared_vals[tid + 1];  end; sync_threads()
    
    if lane_id == 1 && batch_id <= size(P_batch, 2)
        energies[batch_id] = shared_vals[tid]
    end
    return nothing
end

# =====================================================================
# 3. DATA PARSING & GALOIS COMPRESSION
# =====================================================================
function load_qaplib_int64(filepath::String)
    open(filepath, "r") do io
        tokens = split(read(io, String))
        N = parse(Int, tokens[1])
        idx = 2
        F = zeros(Int64, N, N)
        for i in 1:N, j in 1:N; F[i, j] = parse(Int64, tokens[idx]); idx += 1; end
        D = zeros(Int64, N, N)
        for i in 1:N, j in 1:N; D[i, j] = parse(Int64, tokens[idx]); idx += 1; end
        return N, F, D
    end
end

function extract_galois_sparse_tensor(F::Matrix{Int64}, N::Int)
    is_i, is_j, vals = Int32[], Int32[], Int64[]
    for i in 1:N, j in 1:N
        if F[i, j] != 0
            push!(is_i, i); push!(is_j, j); push!(vals, F[i, j])
        end
    end
    return is_i, is_j, vals, length(vals)
end

# =====================================================================
# 4. INTENSIFIED ALNS BATCH GENERATOR WITH MEGA-KICK
# =====================================================================
function generate_intensified_batch!(P_batch, current_p, elite_list, batch_size, N, stagnation_counter)
    n_elites = length(elite_list)
    is_mega_kick = stagnation_counter > 1500 # Trigger escape kick if stagnant
    
    @threads for b in 1:batch_size
        rng = Xoshiro(hash(threadid() + b + time_ns())) 
        P_batch[:, b] .= current_p
        
        if is_mega_kick && (b % 4 == 0)
            # MEGA-KICK: Massively perturb 64 elements to escape deep local basins
            k_kick = rand(rng, 32:64)
            kick_idx = shuffle(rng, 1:N)[1:k_kick]
            P_batch[kick_idx, b] = shuffle(rng, P_batch[kick_idx, b])
        else
            mode = rand(rng)
            if mode < 0.35
                i1, i2 = rand(rng, 1:N, 2)
                P_batch[i1, b], P_batch[i2, b] = P_batch[i2, b], P_batch[i1, b]
            elseif mode < 0.70
                k = rand(rng, 12:28) # Wider LNS span
                idx = sort(shuffle(rng, 1:N)[1:k])
                P_batch[idx, b] .= circshift(P_batch[idx, b], 1)
            else
                if n_elites >= 2
                    e1, e2 = rand(rng, 1:n_elites, 2)
                    P_batch[:, b] .= elite_list[e1]
                    for _ in 1:rand(rng, 4:8)
                        target_i = rand(rng, 1:N)
                        val_needed = elite_list[e2][target_i]
                        pos = findfirst(==(val_needed), P_batch[:, b])
                        if !isnothing(pos)
                            P_batch[target_i, b], P_batch[pos, b] = P_batch[pos, b], P_batch[target_i, b]
                        end
                    end
                end
            end
        end
    end
end

# =====================================================================
# 5. MAIN SOLVER EXECUTION
# =====================================================================
function run_intensified_solver()
    filepath = "tai256c.dat"
    N, F_host, D_host = load_qaplib_int64(filepath)
    is_i, is_j, f_vals, nnz = extract_galois_sparse_tensor(F_host, N)

    println("=========================================================")
    println(" MINSPM Intensified GPU Solver (Targeting < 44,759,292)")
    println("=========================================================")

    d_F_i = CuArray(is_i)
    d_F_j = CuArray(is_j)
    d_F_vals = CuArray(f_vals)
    d_D = CuArray(D_host)

    batch_size = 65_536
    P_host = zeros(Int32, N, batch_size)
    P_gpu = CUDA.zeros(Int32, N, batch_size)
    energies_gpu = CUDA.zeros(Int64, batch_size)
    energies_host = zeros(Int64, batch_size)

    threads = 256
    warps = threads ÷ 32
    blocks = ceil(Int, batch_size / warps)

    current_p = Int32.(randperm(N))
    current_energy = sum(F_host[i, j] * D_host[current_p[i], current_p[j]] for i in 1:N, j in 1:N)
    
    best_global_energy = current_energy
    best_permutation = copy(current_p)

    temperature = 50_000.0   
    cooling_rate = 0.9995     
    max_epochs = 150_000      

    elite_basins = Set{Vector{Int32}}()
    visited_orbits = Set{UInt64}()
    stagnation_counter = 0

    println("\nStarting Intensified Optimization Run...")
    start_time = time()

    for epoch in 1:max_epochs
        elite_list = collect(elite_basins)
        generate_intensified_batch!(P_host, current_p, elite_list, batch_size, N, stagnation_counter)
        copyto!(P_gpu, P_host)

        @cuda threads=threads blocks=blocks shmem=threads*sizeof(Int64) minspm_tai256c_int64_kernel!(
            d_F_i, d_F_j, d_F_vals, d_D, P_gpu, energies_gpu, nnz
        )
        CUDA.synchronize()
        copyto!(energies_host, energies_gpu)

        best_idx = argmin(energies_host)
        best_batch_energy = energies_host[best_idx]
        best_batch_perm = P_host[:, best_idx]

        if !is_valid_permutation(best_batch_perm, N) || best_batch_energy < 44_095_032
            stagnation_counter += 1
            continue
        end

        orbit_sig = hash(best_batch_perm)
        if !(orbit_sig in visited_orbits)
            push!(visited_orbits, orbit_sig)
            delta = Float64(best_batch_energy - current_energy)

            if delta < 0.0 || rand() < exp(-delta / temperature)
                current_p .= best_batch_perm
                current_energy = best_batch_energy

                if current_energy < best_global_energy
                    best_global_energy = current_energy
                    best_permutation .= current_p
                    stagnation_counter = 0 # Reset stagnation on improvement
                    
                    push!(elite_basins, copy(best_permutation))
                    if length(elite_basins) > 150; pop!(elite_basins); end

                    @printf(" [Epoch %5d] 🌟 New Feasible Best Objective: %d (Gap vs Record 44,759,292: %+.3f%%)\n", 
                            epoch, best_global_energy, ((best_global_energy - 44759292.0)/44759292.0)*100)
                    
                    if best_global_energy <= 44_759_292
                        println("\n🎯 SUCCESS: Established new verified feasible upper-bound record!")
                        break
                    end
                else
                    stagnation_counter += 1
                end
            else
                stagnation_counter += 1
            end
        else
            stagnation_counter += 1
        end

        temperature *= cooling_rate
    end
    
    println("\n=========================================================")
    @printf(" FINAL VALIDATED RESULT: %d\n", best_global_energy)
    println(" Best Permutation Vector (Single Line):")
    println(join(best_permutation, " "))
    println("=========================================================")
    CUDA.reclaim()
end

run_intensified_solver()
