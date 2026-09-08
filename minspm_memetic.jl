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
# 2. EXACT INT64 MINSPM SPARSE GALOIS KERNEL
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
# 4. MEMETIC MULTI-START GENERATOR
# =====================================================================
function generate_production_batch!(P_batch, current_p, all_swaps, elite_list, batch_size, N, stagnation_counter)
    n_elites = length(elite_list)
    
    @threads for b in 1:batch_size
        rng = Xoshiro(hash(threadid() + b + time_ns())) 
        P_batch[:, b] .= current_p
        
        if b <= 32640
            # 1. EXACT FNSD 2-OPT (First half of GPU blocks)
            i1, i2 = all_swaps[b]
            P_batch[i1, b], P_batch[i2, b] = P_batch[i2, b], P_batch[i1, b]
        else
            # 2. MEMETIC CROSSOVER & REPAIR (Second half of GPU blocks)
            if n_elites >= 2
                e1, e2 = rand(rng, 1:n_elites, 2)
                # Copy Parent 1
                P_batch[:, b] .= elite_list[e1]
                
                # Identify differing positions
                diff_idx = findall(elite_list[e1] .!= elite_list[e2])
                
                # If they are distinct, randomly resolve a subset to Parent 2
                if length(diff_idx) > 2
                    for _ in 1:rand(rng, 4:min(16, length(diff_idx)))
                        target_i = rand(rng, diff_idx)
                        val_needed = elite_list[e2][target_i]
                        pos = findfirst(==(val_needed), P_batch[:, b])
                        if !isnothing(pos)
                            P_batch[target_i, b], P_batch[pos, b] = P_batch[pos, b], P_batch[target_i, b]
                        end
                    end
                end
                
                # Apply a minor 3-opt mutation to prevent incestuous trapping
                i1, i2, i3 = rand(rng, 1:N, 3)
                P_batch[i1, b], P_batch[i2, b], P_batch[i3, b] = P_batch[i2, b], P_batch[i3, b], P_batch[i1, b]
            else
                # Fallback to Deep LNS if pool is small
                k = rand(rng, 16:32)
                idx = sort(shuffle(rng, 1:N)[1:k])
                P_batch[idx, b] .= circshift(P_batch[idx, b], rand(rng, 1:(k-1)))
            end
        end
    end
end

# =====================================================================
# 5. MAIN SOLVER EXECUTION
# =====================================================================
function run_production_solver()
    filepath = "tai256c.dat"
    N, F_host, D_host = load_qaplib_int64(filepath)
    is_i, is_j, f_vals, nnz = extract_galois_sparse_tensor(F_host, N)

    println("=========================================================")
    println(" MINSPM PRODUCTION LONG-HAUL SOLVER")
    println(" Architecture: Memetic FNSD + Multi-Start Reset")
    println("=========================================================")

    all_swaps = Vector{Tuple{Int, Int}}(undef, 32640)
    idx = 1
    for i in 1:(N-1), j in (i+1):N; all_swaps[idx] = (i, j); idx += 1; end

    d_F_i, d_F_j = CuArray(is_i), CuArray(is_j)
    d_F_vals, d_D = CuArray(f_vals), CuArray(D_host)

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

    # LONG-HAUL PARAMETERS
    max_epochs = 1_000_000      
    restarts = 0

    elite_basins = Set{Vector{Int32}}()
    stagnation_counter = 0

    println("\nStarting Production Optimization Run (1 Million Epochs)...")
    start_time = time()

    for epoch in 1:max_epochs
        # MULTI-START RESET TRIGGER
        if stagnation_counter > 150
            restarts += 1
            # 20% of the time, restart from a Memetic Crossover of elites
            # 80% of the time, nuke the board to explore totally new space
            if rand() < 0.20 && length(elite_basins) >= 2
                elite_list = collect(elite_basins)
                current_p .= elite_list[rand(1:length(elite_list))]
                for _ in 1:16
                    i1, i2 = rand(1:N, 2)
                    current_p[i1], current_p[i2] = current_p[i2], current_p[i1]
                end
            else
                current_p .= Int32.(randperm(N))
            end
            current_energy = sum(F_host[i, j] * D_host[current_p[i], current_p[j]] for i in 1:N, j in 1:N)
            stagnation_counter = 0
        end

        elite_list = collect(elite_basins)
        generate_production_batch!(P_host, current_p, all_swaps, elite_list, batch_size, N, stagnation_counter)
        copyto!(P_gpu, P_host)

        @cuda threads=threads blocks=blocks shmem=threads*sizeof(Int64) minspm_tai256c_int64_kernel!(
            d_F_i, d_F_j, d_F_vals, d_D, P_gpu, energies_gpu, nnz
        )
        CUDA.synchronize()
        copyto!(energies_host, energies_gpu)

        # Evaluate deterministic descent vs memetic crossover
        best_det_idx = argmin(energies_host[1:32640])
        best_det_energy = energies_host[best_det_idx]
        
        best_rand_idx = argmin(energies_host[32641:end]) + 32640
        best_rand_energy = energies_host[best_rand_idx]

        improvement_found = false

        if best_det_energy < current_energy && is_valid_permutation(P_host[:, best_det_idx], N)
            current_p .= P_host[:, best_det_idx]
            current_energy = best_det_energy
            improvement_found = true
        elseif best_rand_energy < current_energy && is_valid_permutation(P_host[:, best_rand_idx], N)
            current_p .= P_host[:, best_rand_idx]
            current_energy = best_rand_energy
            improvement_found = true
        end

        if improvement_found
            stagnation_counter = 0
            if current_energy < best_global_energy
                best_global_energy = current_energy
                best_permutation .= current_p
                
                push!(elite_basins, copy(best_permutation))
                if length(elite_basins) > 500; pop!(elite_basins); end
                
                @printf(" [Epoch %7d | Restarts: %3d] 🌟 NEW BEST: %d (Gap: %+.3f%%)\n", 
                        epoch, restarts, best_global_energy, ((best_global_energy - 44759292.0)/44759292.0)*100)
            end
        else
            stagnation_counter += 1
        end
        
        if best_global_energy <= 44_759_292
            println("\n🎯 SUCCESS: Established new verified feasible upper-bound record!")
            break
        end

        if epoch % 10000 == 0
            elapsed = time() - start_time
            evals = epoch * batch_size
            @printf(" -> Status: %.2f Billion Evals | %.2f Million Evals/sec | Active Elites: %d\n", 
                    evals / 1e9, (evals / elapsed) / 1e6, length(elite_basins))
        end
    end
    
    println("\n=========================================================")
    @printf(" FINAL VALIDATED RESULT: %d\n", best_global_energy)
    println(" Best Permutation Vector (Single Line):")
    println(join(best_permutation, " "))
    println("=========================================================")
    CUDA.reclaim()
end

run_production_solver()
