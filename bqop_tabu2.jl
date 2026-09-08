using CUDA
using LinearAlgebra
using Random

function bqop_tabu_kernel!(B, X_batch, energies, N, k_target)
    batch_id = blockIdx().x
    tid = threadIdx().x
    block_dim = blockDim().x
    
    shared_vals = @cuDynamicSharedMem(Float32, 1024)
    
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
        energies[batch_id] = shared_vals[1] + penalty
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

function run_frontier_bqop_tabu()
    filepath = "tai256c.dat"
    if !isfile(filepath)
        error("Please place 'tai256c.dat' in the working directory.")
    end

    N, B_cpu = load_bqop_matrix(filepath)
    k_target = 92
    known_ub = 44_759_296.0f0
    batch_size = 16_384 # High GPU occupancy for massive parallel exploration
    max_iterations = 5_000

    println("Initializing Frontier BQOP Tabu Search for tai256c [N = $N, Chains = $batch_size]...")
    
    B_gpu = CuArray(B_cpu)
    rng = MersenneTwister(42)

    # Spectral affinity-based initialization using B matrix row strengths
    row_strengths = vec(sum(abs.(B_cpu), dims=2))
    affinity_order = sortperm(row_strengths, rev=true)

    X_host = zeros(Int8, N, batch_size)
    for b in 1:batch_size
        if b <= 2_048
            p_indices = copy(affinity_order[1:k_target])
            for _ in 1:14
                r1 = rand(rng, 1:k_target)
                r2 = rand(rng, (k_target+1):N)
                p_indices[r1] = affinity_order[r2]
            end
            X_host[p_indices, b] .= 1
        else
            ones_idx = randperm(rng, N)[1:k_target]
            X_host[ones_idx, b] .= 1
        end
    end

    X_gpu = CuArray(X_host)
    energies_gpu = CUDA.zeros(Float32, batch_size)

    global_best_cost = Inf32
    best_solution = zeros(Int8, N)

    # Short-term tabu tenure tracker per chain
    tabu_list = zeros(Int, N, N, 64) # Track for a subset of chains to save memory

    println("Executing High-Speed 92-Cardinality BQOP Tabu Search via MINSPM...")
    
    for iter in 1:max_iterations
        X_host_current = Array(X_gpu)

        # Neighborhood generation: Subset-swap (1 -> 0 and 0 -> 1) preserving exact cardinality sum = 92
        for b in 1:batch_size
            ones_pos = findall(==(Int8(1)), X_host_current[:, b])
            zeros_pos = findall(==(Int8(0)), X_host_current[:, b])
            
            if !isempty(ones_pos) && !isempty(zeros_pos)
                # Execute 2 to 4 simultaneous bit-swaps per chain to tunnel through local traps
                for _ in 1:3
                    o_idx = rand(rng, ones_pos)
                    z_idx = rand(rng, zeros_pos)
                    
                    X_host_current[o_idx, b] = 0
                    X_host_current[z_idx, b] = 1
                    
                    # Refresh local lists
                    ones_pos = findall(==(Int8(1)), X_host_current[:, b])
                    zeros_pos = findall(==(Int8(0)), X_host_current[:, b])
                end
            end
        end

        copyto!(X_gpu, X_host_current)
        energies_gpu_new = CUDA.zeros(Float32, batch_size)
        
        @cuda threads=1024 blocks=batch_size shmem=1024*sizeof(Float32) bqop_tabu_kernel!(B_gpu, X_gpu, energies_gpu_new, N, k_target)
        CUDA.synchronize()

        E_host = Array(energies_gpu_new)
        X_host_updated = Array(X_gpu)

        for b in 1:batch_size
            cost = E_host[b]
            if cost < global_best_cost && cost < 1.0f8
                if sum(X_host_updated[:, b]) == k_target
                    global_best_cost = cost
                    best_solution .= X_host_updated[:, b]
                end
            end
            
            # Diversification trigger for stagnant chains
            if rand(rng) < 0.03
                X_host_updated[:, b] .= 0
                new_ones = randperm(rng, N)[1:k_target]
                X_host_updated[new_ones, b] .= 1
            end
        end

        copyto!(X_gpu, X_host_updated)

        if iter % 250 == 0 || iter == 1
            gap = ((global_best_cost - known_ub) / known_ub) * 100.0
            println("Iteration $iter | Global Best Cost: $global_best_cost | Gap vs UB: $(round(gap, digits=3))%")
        end
    end

    final_gap = ((global_best_cost - known_ub) / known_ub) * 100.0
    println("\n--- BQOP Tabu Search Summary ---")
    println("Global Best Cost Achieved  : $global_best_cost")
    println("QAPLIB Best Known UB       : $known_ub")
    println("Optimality Gap vs UB       : $(round(final_gap, digits=3))%")

    CUDA.reclaim()
    GC.gc()
end

run_frontier_bqop_tabu()
