using CUDA
using LinearAlgebra
using Random

function minspm_alns_kernel!(B, X_batch, energies, N, k_target)
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

function run_expensive_alns_solver()
    filepath = "tai256c.dat"
    if !isfile(filepath)
        error("Please place 'tai256c.dat' in the working directory.")
    end

    N, B_cpu = load_bqop_matrix(filepath)
    k_target = 92
    known_ub = 44_759_296.0f0
    batch_size = 8_192
    max_generations = 3_000

    println("Initializing Expensive ALNS Ruin-and-Recreate Solver for tai256c [N = $N, Chains = $batch_size]...")
    
    B_gpu = CuArray(B_cpu)
    rng = MersenneTwister(1337)

    # Spectral affinity initialization using B matrix row weights
    row_weights = vec(sum(abs.(B_cpu), dims=2))
    affinity_order = sortperm(row_weights, rev=true)

    X_host = zeros(Int8, N, batch_size)
    for b in 1:batch_size
        # Diverse high-affinity seeding
        p_indices = copy(affinity_order[1:k_target])
        for _ in 1:16
            r1 = rand(rng, 1:k_target)
            r2 = rand(rng, (k_target+1):N)
            p_indices[r1] = affinity_order[r2]
        end
        X_host[p_indices, b] .= 1
    end

    X_gpu = CuArray(X_host)
    energies_gpu = CUDA.zeros(Float32, batch_size)

    global_best_cost = Inf32
    best_solution = zeros(Int8, N)

    println("Executing Heavy-Duty Large Neighborhood Search via MINSPM Oracle...")
    
    for gen in 1:max_generations
        @cuda threads=1024 blocks=batch_size shmem=1024*sizeof(Float32) minspm_alns_kernel!(B_gpu, X_gpu, energies_gpu, N, k_target)
        CUDA.synchronize()

        E_host = Array(energies_gpu)
        X_host_current = Array(X_gpu)

        # Global Best Tracking
        for b in 1:batch_size
            cost = E_host[b]
            if cost < global_best_cost && cost < 1.0f8
                if sum(X_host_current[:, b]) == k_target
                    global_best_cost = cost
                    best_solution .= X_host_current[:, b]
                end
            end
        end

        # Expensive Ruin-and-Recreate (ALNS) Mutation Phase
        for b in 1:batch_size
            ones_pos = findall(==(Int8(1)), X_host_current[:, b])
            zeros_pos = findall(==(Int8(0)), X_host_current[:, b])
            
            # Ruin: Destroy 32 active bits and 32 inactive bits simultaneously
            ruin_size = min(32, length(ones_pos), length(zeros_pos))
            if ruin_size > 0
                destroyed_ones = randperm(rng, length(ones_pos))[1:ruin_size]
                destroyed_zeros = randperm(rng, length(zeros_pos))[1:ruin_size]
                
                for r in destroyed_ones
                    X_host_current[ones_pos[r], b] = 0
                end
                for r in destroyed_zeros
                    X_host_current[zeros_pos[r], b] = 1
                end
                
                # Recreate: Intelligently rebuild using affinity matrix correlation weights
                current_zeros = findall(==(Int8(0)), X_host_current[:, b])
                rebuild_candidates = sort(current_zeros, by=idx -> row_weights[idx], rev=true)
                
                current_sum = sum(X_host_current[:, b])
                needed = k_target - current_sum
                
                if needed > 0 && length(rebuild_candidates) >= needed
                    # Pick top affinity slots to restore exact 92-cardinality sum
                    for c in 1:needed
                        X_host_current[rebuild_candidates[c], b] = 1
                    end
                elseif needed < 0
                    current_ones = findall(==(Int8(1)), X_host_current[:, b])
                    to_clear = sort(current_ones, by=idx -> row_weights[idx])
                    for c in 1:abs(needed)
                        X_host_current[to_clear[c], b] = 0
                    end
                end
            end
        end

        copyto!(X_gpu, X_host_current)

        if gen % 200 == 0 || gen == 1
            gap = ((global_best_cost - known_ub) / known_ub) * 100.0
            println("Generation $gen | Best Found Cost: $global_best_cost | Gap vs UB: $(round(gap, digits=3))%")
        end
    end

    final_gap = ((global_best_cost - known_ub) / known_ub) * 100.0
    println("\n--- ALNS Optimization Summary ---")
    println("Global Best Cost Found     : $global_best_cost")
    println("QAPLIB Best Known UB       : $known_ub")
    println("Optimality Gap vs UB       : $(round(final_gap, digits=3))%")

    CUDA.reclaim()
    GC.gc()
end

run_expensive_alns_solver()
