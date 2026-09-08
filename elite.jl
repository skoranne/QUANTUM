using CUDA
using LinearAlgebra
using Random

function minspm_bqop_kernel!(B, X_batch, energies, N, k_target)
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

function run_frontier_bqop_solver()
    filepath = "tai256c.dat"
    if !isfile(filepath)
        error("Please place 'tai256c.dat' in the working directory.")
    end

    N, B_cpu = load_bqop_matrix(filepath)
    k_target = 92
    known_ub = 44_759_296.0f0
    batch_size = 4_096
    max_generations = 1_500

    println("Initializing Frontier-Class BQOP Solver for tai256c [N = $N, Target = $k_target]...")
    
    B_gpu = CuArray(B_cpu)
    rng = MersenneTwister(1337)

    # Intelligent spectral/row-sum initialization
    row_strengths = vec(sum(abs.(B_cpu), dims=2))
    greedy_order = sortperm(row_strengths, rev=true)
    
    X_host = zeros(Int8, N, batch_size)
    for b in 1:batch_size
        if b <= 512
            # Perturbed variants of the greedy optimal structure
            p_indices = copy(greedy_order[1:k_target])
            for _ in 1:12
                r1 = rand(rng, 1:k_target)
                r2 = rand(rng, (k_target+1):N)
                p_indices[r1] = greedy_order[r2]
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
    
    # Elite Pool Archive
    elite_pool = [zeros(Int8, N) for _ in 1:32]
    elite_costs = fill(Inf32, 32)

    println("Executing High-Order Block-Swap & Elite Relinking via MINSPM...")
    
    for gen in 1:max_generations
        energies_gpu_new = CUDA.zeros(Float32, batch_size)
        @cuda threads=1024 blocks=batch_size shmem=1024*sizeof(Float32) minspm_bqop_kernel!(B_gpu, X_gpu, energies_gpu_new, N, k_target)
        CUDA.synchronize()

        E_host = Array(energies_gpu_new)
        X_host_current = Array(X_gpu)

        # Update Global Best & Elite Pool
        for b in 1:batch_size
            cost = E_host[b]
            if cost < global_best_cost && cost < 1.0f8
                if sum(X_host_current[:, b]) == k_target
                    global_best_cost = cost
                    best_solution .= X_host_current[:, b]
                    
                    # Insert into elite pool
                    worst_elite_idx = argmax(elite_costs)
                    if cost < elite_costs[worst_elite_idx]
                        elite_costs[worst_elite_idx] = cost
                        elite_pool[worst_elite_idx] .= X_host_current[:, b]
                    end
                end
            end
        end

        # Multi-Bit Block Swapping & Path Relinking
        for b in 1:batch_size
            ones_pos = findall(==(Int8(1)), X_host_current[:, b])
            zeros_pos = findall(==(Int8(0)), X_host_current[:, b])
            
            # Execute 4-bit block swaps to tunnel through barriers
            if length(ones_pos) >= 4 && length(zeros_pos) >= 4
                for _ in 1:4
                    o_rem = rand(rng, ones_pos)
                    z_add = rand(rng, zeros_pos)
                    X_host_current[o_rem, b] = 0
                    X_host_current[z_add, b] = 1
                    ones_pos = findall(==(Int8(1)), X_host_current[:, b])
                    zeros_pos = findall(==(Int8(0)), X_host_current[:, b])
                end
            end

            # Elite Path Relinking: Pull 15% of stagnant chains toward a random elite solution
            if gen % 15 == 0 && rand(rng) < 0.15 && !isinf(global_best_cost)
                target_elite = elite_pool[rand(rng, 1:length(elite_pool))]
                if sum(target_elite) == k_target
                    diffs = findall(i -> X_host_current[i, b] != target_elite[i], 1:N)
                    if length(diffs) >= 2
                        d_idx = diffs[1]
                        val = target_elite[d_idx]
                        # Swap balancing bit to preserve 92 cardinality
                        if val == 1
                            z_to_clear = findall(i -> X_host_current[i, b] == 1, 1:N)
                            if !isempty(z_to_clear)
                                X_host_current[rand(rng, z_to_clear), b] = 0
                                X_host_current[d_idx, b] = 1
                            end
                        else
                            o_to_clear = findall(i -> X_host_current[i, b] == 0, 1:N)
                            if !isempty(o_to_clear)
                                X_host_current[rand(rng, o_to_clear), b] = 1
                                X_host_current[d_idx, b] = 0
                            end
                        end
                    end
                end
            end
        end

        copyto!(X_gpu, X_host_current)

        if gen % 100 == 0 || gen == 1
            gap = ((global_best_cost - known_ub) / known_ub) * 100.0
            println("Generation $gen | Best Found Cost: $global_best_cost | Gap vs UB: $(round(gap, digits=3))%")
        end
    end

    final_gap = ((global_best_cost - known_ub) / known_ub) * 100.0
    println("\n--- Final Frontier Optimization Summary ---")
    println("Global Best Cost Found     : $global_best_cost")
    println("QAPLIB Best Known UB       : $known_ub")
    println("Optimality Gap vs UB       : $(round(final_gap, digits=3))%")

    CUDA.reclaim()
    GC.gc()
end

run_frontier_bqop_solver()
