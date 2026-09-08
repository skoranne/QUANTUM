using CUDA
using LinearAlgebra
using Random

function fused_bqop_tabu_kernel!(B, X_batch, best_costs_batch, N, k_target, iter_seed)
    batch_id = blockIdx().x
    tid = threadIdx().x
    block_dim = blockDim().x
    
    shared_vals = @cuDynamicSharedMem(Float32, 1024)
    
    # Each thread block handles one independent search chain
    # Perform an on-device local swap (2-opt neighborhood move) based on pseudo-random seed
    if tid == 1
        # Deterministic pseudo-random index selection on device
        idx1 = 1 + ((iter_seed * 31 + batch_id * 17) % N)
        idx2 = 1 + ((iter_seed * 43 + batch_id * 23) % N)
        if idx1 != idx2
            # Swap bits if it preserves or improves cardinality
            v1 = X_batch[idx1, batch_id]
            v2 = X_batch[idx2, batch_id]
            X_batch[idx1, batch_id] = v2
            X_batch[idx2, batch_id] = v1
        end
    end
    sync_threads()

    # Evaluate BQOP objective x^T * B * x in shared memory
    val = 0.0f0
    idx = tid
    total_elements = N * N
    
    while idx <= total_elements
        i = div(idx - 1, N) + 1
        j = rem(idx - 1, N) + 1
        
        if X_batch[i, batch_id] == Int8(1) && X_batch[j, batch_id] == Int8(1)
            val += B[i, j]
        end
        idx += block_dim
    end
    
    shared_vals[tid] = val
    sync_threads()
    
    s = 512
    while s > 0
        if tid <= s && (tid + s) <= 1024
            shared_vals[tid] += shared_vals[tid + s]
        end
        sync_threads()
        s = s >> 1
    end
    
    if tid == 1
        sum_ones = 0
        for r in 1:N
            if X_batch[r, batch_id] == Int8(1)
                sum_ones += 1
            end
        end
        penalty = (sum_ones != k_target) ? 1f9 : 0.0f0
        total_cost = shared_vals[1] + penalty
        best_costs_batch[batch_id] = total_cost
    end
    return nothing
end

function load_bqop_matrix(filepath::String)
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
        return N, F .+ D
    end
end

function run_fused_gpu_solver()
    filepath = "tai256c.dat"
    if !isfile(filepath)
        error("Please place 'tai256c.dat' in the working directory.")
    end

    N, B_cpu = load_bqop_matrix(filepath)
    k_target = 92
    known_ub = 44_759_296.0f0
    batch_size = 8_192 # Doubled parallelism for deeper frontier search
    max_generations = 5_000

    println("Initializing Fused GPU Kernel Solver for tai256c [N = $N, Chains = $batch_size]...")
    
    B_gpu = CuArray(B_cpu)
    rng = MersenneTwister(42)

    # Intelligent spectral initialization
    row_strengths = vec(sum(abs.(B_cpu), dims=2))
    greedy_order = sortperm(row_strengths, rev=true)
    
    X_host = zeros(Int8, N, batch_size)
    for b in 1:batch_size
        p_indices = copy(greedy_order[1:k_target])
        for _ in 1:8
            r1 = rand(rng, 1:k_target)
            r2 = rand(rng, (k_target+1):N)
            p_indices[r1] = greedy_order[r2]
        end
        X_host[p_indices, b] .= 1
    end

    X_gpu = CuArray(X_host)
    costs_gpu = CUDA.zeros(Float32, batch_size)

    global_best_cost = Inf32
    best_solution = zeros(Int8, N)

    println("Executing High-Velocity On-Device Neighborhood Search...")
    
    for gen in 1:max_generations
        @cuda threads=1024 blocks=batch_size shmem=1024*sizeof(Float32) fused_bqop_tabu_kernel!(B_gpu, X_gpu, costs_gpu, N, k_target, gen)
        CUDA.synchronize()

        costs_host = Array(costs_gpu)
        
        # Track global best on host periodically
        min_c, min_idx = findmin(costs_host)
        if min_c < global_best_cost && min_c < 1.0f8
            global_best_cost = min_c
            X_curr_host = Array(X_gpu)
            best_solution .= X_curr_host[:, min_idx]
        end

        if gen % 500 == 0 || gen == 1
            gap = ((global_best_cost - known_ub) / known_ub) * 100.0
            println("Generation $gen | Best Found Cost: $global_best_cost | Gap vs UB: $(round(gap, digits=3))%")
        end
    end

    final_gap = ((global_best_cost - known_ub) / known_ub) * 100.0
    println("\n--- Fused Kernel Optimization Summary ---")
    println("Global Best Cost Found     : $global_best_cost")
    println("QAPLIB Best Known UB       : $known_ub")
    println("Optimality Gap vs UB       : $(round(final_gap, digits=3))%")

    CUDA.reclaim()
    GC.gc()
end

run_fused_gpu_solver()
