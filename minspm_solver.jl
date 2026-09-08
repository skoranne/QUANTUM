using CUDA
using LinearAlgebra
using Random
using Printf
using Base.Threads

# =====================================================================
# 1. EXACT QAP / MINSPM SPARSE GALOIS EVALUATION KERNEL
# =====================================================================
function minspm_tai256c_kernel!(F_sparse_i, F_sparse_j, F_vals, D, P_batch, energies, nnz_terms)
    tid = threadIdx().x
    warp_id = div(tid - 1, 32) + 1
    lane_id = rem(tid - 1, 32) + 1
    
    batch_id = (blockIdx().x - 1) * div(blockDim().x, 32) + warp_id
    
    shared_vals = @cuDynamicSharedMem(Float32, blockDim().x)
    val = 0.0f0
    
    if batch_id <= size(P_batch, 2)
        idx = lane_id
        while idx <= nnz_terms
            i = F_sparse_i[idx]
            j = F_sparse_j[idx]
            f_val = F_vals[idx]
            
            p_i = P_batch[i, batch_id]
            p_j = P_batch[j, batch_id]
            
            val += f_val * D[p_i, p_j]
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
# 2. DATA PARSING & GALOIS SPARSITY COMPRESSION
# =====================================================================
function load_qaplib(filepath::String)
    open(filepath, "r") do io
        tokens = split(read(io, String))
        N = parse(Int, tokens[1])
        idx = 2
        F = zeros(Float32, N, N)
        for i in 1:N, j in 1:N
            F[i, j] = parse(Float32, tokens[idx]); idx += 1
        end
        D = zeros(Float32, N, N)
        for i in 1:N, j in 1:N
            D[i, j] = parse(Float32, tokens[idx]); idx += 1
        end
        return N, F, D
    end
end

function extract_galois_sparse_tensor(F::Matrix{Float32}, N::Int)
    # Minimal perfect hash extraction for non-zero flow coefficients
    is_i = Int32[]
    is_j = Int32[]
    vals = Float32[]
    
    for i in 1:N, j in 1:N
        if F[i, j] != 0.0f0
            push!(is_i, i)
            push!(is_j, j)
            push!(vals, F[i, j])
        end
    end
    return is_i, is_j, vals, length(vals)
end

# =====================================================================
# 3. ALNS / SIMULATED ANNEALING BATCH GENERATOR
# =====================================================================
function generate_neighborhood_batch!(P_batch, current_p, elite_basins, batch_size, N)
    @threads for b in 1:batch_size
        P_batch[:, b] .= current_p
        
        mode = rand()
        if mode < 0.50
            # 2-Opt Swap (VNS)
            i1, i2 = rand(1:N, 2)
            P_batch[i1, b], P_batch[i2, b] = P_batch[i2, b], P_batch[i1, b]
        elseif mode < 0.85
            # 4-Opt Block Shift (LNS)
            k = rand(3:6)
            idx = sort(randperm(N)[1:k])
            P_batch[idx, b] .= circshift(P_batch[idx, b], 1)
        else
            # Elite Basin Crossover
            if !isempty(elite_basins)
                elite = rand(elite_basins)
                P_batch[:, b] .= elite
                i1, i2 = rand(1:N, 2)
                P_batch[i1, b], P_batch[i2, b] = P_batch[i2, b], P_batch[i1, b]
            end
        end
    end
end

# =====================================================================
# 4. MAIN SOLVER EXECUTION
# =====================================================================
function run_tai256c_solver()
    filepath = "tai256c.dat"
    if !isfile(filepath)
        error("Place 'tai256c.dat' in the working directory before executing.")
    end

    N, F_host, D_host = load_qaplib(filepath)
    is_i, is_j, f_vals, nnz = extract_galois_sparse_tensor(F_host, N)

    println("=========================================================")
    println(" MINSPM High-Throughput GPU Solver for tai256c")
    println("=========================================================")
    @printf(" Problem Size (N) : %d\n", N)
    @printf(" Galois NNZ Terms : %d (Sparsity Compression Active)\n", nnz)

    # GPU Allocations
    d_F_i = CuArray(is_i)
    d_F_j = CuArray(is_j)
    d_F_vals = CuArray(f_vals)
    d_D = CuArray(D_host)

    batch_size = 65_536
    P_host = zeros(Int32, N, batch_size)
    P_gpu = CUDA.zeros(Int32, N, batch_size)
    energies_gpu = CUDA.zeros(Float32, batch_size)
    energies_host = zeros(Float32, batch_size)

    threads = 256
    warps = threads ÷ 32
    blocks = ceil(Int, batch_size / warps)

    # Initial State (Random Permutation)
    current_p = Int32.(randperm(N))
    
    # Evaluate initial state via CPU baseline reference
    current_energy = Float64(sum(F_host[i, j] * D_host[current_p[i], current_p[j]] for i in 1:N, j in 1:N))
    best_global_energy = current_energy
    best_permutation = copy(current_p)

    temperature = 50000.0
    cooling_rate = 0.997
    max_epochs = 10000

    elite_basins = Set{Vector{Int32}}()
    visited_orbits = Set{UInt64}()

    println("\nStarting Hardware-Accelerated Optimization Run...")
    println("Current Record to Beat: 44,759,292 | Initial Energy: $(round(current_energy))\n")

    start_time = time()

    for epoch in 1:max_epochs
        generate_neighborhood_batch!(P_host, current_p, elite_basins, batch_size, N)
        copyto!(P_gpu, P_host)

        @cuda threads=threads blocks=blocks shmem=threads*sizeof(Float32) minspm_tai256c_kernel!(
            d_F_i, d_F_j, d_F_vals, d_D, P_gpu, energies_gpu, nnz
        )
        CUDA.synchronize()
        copyto!(energies_host, energies_gpu)

        best_idx = argmin(energies_host)
        best_batch_energy = Float64(energies_host[best_idx])
        best_batch_perm = P_host[:, best_idx]

        orbit_sig = hash(best_batch_perm)
        if !(orbit_sig in visited_orbits)
            push!(visited_orbits, orbit_sig)

            delta = best_batch_energy - current_energy
            if delta < 0.0 || rand() < exp(-delta / temperature)
                current_p .= best_batch_perm
                current_energy = best_batch_energy

                if current_energy < best_global_energy
                    best_global_energy = current_energy
                    best_permutation .= current_p
                    push!(elite_basins, copy(best_permutation))
                    if length(elite_basins) > 30; pop!(elite_basins); end

                    @printf(" [Epoch %4d] 🌟 New Best Objective: %.2f (Gap vs 44,759,292: %.3f%%)\n", 
                            epoch, best_global_energy, ((best_global_energy - 44759292.0)/44759292.0)*100)
                    
                    if best_global_energy <= 44759292.0
                        println("\n🎯 SUCCESS: Established new feasible upper-bound record!")
                        println("Permutation Vector: $(best_permutation)")
                    end
                end
            end
        end

        temperature *= cooling_rate

        if epoch % 1000 == 0
            elapsed = time() - start_time
            evals = epoch * batch_size
            @printf(" -> Progress: %.2f Billion Evals | Speed: %.2f Million Evals/sec | Current Temp: %.1f\n", 
                    evals / 1e9, (evals / elapsed) / 1e6, temperature)
        end
    end

    println("\n=========================================================")
    @printf(" FINAL OPTIMIZATION RESULT: %.2f\n", best_global_energy)
    println(" Best Permutation Vector found:")
    display(best_permutation)
    println("\n=========================================================")
    CUDA.reclaim()
end

run_tai256c_solver()
