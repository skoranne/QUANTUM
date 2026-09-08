using CUDA
using LinearAlgebra
using Random

function minspm_bqop_exact_kernel!(B_dict, X_batch, energies, N, k_target)
    batch_id = blockIdx().x
    tid = threadIdx().x
    block_dim = blockDim().x
    
    shared_vals = @cuDynamicSharedMem(Int64, 1024)
    
    val = Int64(0)
    idx = tid
    total_elements = N * N
    
    while idx <= total_elements
        i = div(idx - 1, N) + 1
        j = rem(idx - 1, N) + 1
        
        if X_batch[i, batch_id] == Int8(1) && X_batch[j, batch_id] == Int8(1)
            val += B_dict[i, j]
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
        
        penalty = (sum_ones != k_target) ? Int64(1000000000) : Int64(0)
        energies[batch_id] = shared_vals[1] + penalty
    end
    return nothing
end

function load_bqop_matrix(filepath::String)
    open(filepath, "r") do io
        tokens = split(read(io, String))
        N = parse(Int, tokens[1])
        idx = 2
        
        F = zeros(Int64, N, N)
        for i in 1:N, j in 1:N
            F[i, j] = parse(Int64, tokens[idx]); idx += 1
        end
        
        D = zeros(Int64, N, N)
        for i in 1:N, j in 1:N
            D[i, j] = parse(Int64, tokens[idx]); idx += 1
        end
        return N, F .+ D
    end
end

function get_record_binary_seed(N::Int, k_target::Int)
    record_perm = Int32[
        219, 109, 227, 174, 172, 85, 161, 193, 28, 123, 204, 70, 114, 197, 140, 240, 41, 14, 49, 153, 45, 160, 48, 19, 201, 102, 178, 148, 11, 2, 216, 100, 133, 126, 253, 170, 72, 92, 150, 246, 104, 131, 167, 212, 119, 17, 180, 187, 53, 143, 244, 214, 24, 157, 236, 55, 79, 66, 136, 75, 97, 36, 43, 138, 191, 206, 129, 31, 233, 9, 223, 231, 255, 89, 62, 58, 250, 5, 38, 7, 182, 210, 83, 221, 163, 51, 106, 111, 77, 96, 184, 241, 203, 159, 151, 249, 166, 177, 84, 228, 47, 192, 63, 61, 8, 205, 71, 80, 156, 220, 103, 176, 200, 16, 124, 25, 194, 190, 26, 110, 78, 95, 195, 226, 30, 65, 247, 230, 122, 168, 121, 141, 155, 32, 202, 29, 162, 88, 6, 235, 73, 21, 232, 208, 50, 142, 82, 15, 69, 225, 52, 99, 164, 224, 81, 234, 237, 183, 27, 222, 152, 188, 108, 169, 181, 127,
        118, 171, 154, 56, 39, 251, 207, 165, 239, 105, 242, 1, 149, 139, 217, 93, 245, 254, 40, 116, 87, 135, 3, 4, 76, 10, 54, 86, 189, 144, 101, 68, 64, 218,
        199, 213, 146, 179, 13, 20, 60, 42, 175, 256, 198, 158, 243, 34, 37, 186, 229, 35, 209, 211, 113, 57, 98, 33, 173, 125, 18, 132, 22, 44, 107, 23, 67, 91, 128, 74, 90, 137, 117, 238, 94, 112, 12, 185, 115, 145, 120, 130, 147, 134, 46, 196, 252, 215, 248, 59
    ]
    
    base_x = zeros(Int8, N)
    base_x[record_perm[1:k_target]] .= 1
    return base_x
end

function binary_cardinality_crossover(x1::AbstractVector{Int8}, x2::AbstractVector{Int8}, rng::AbstractRNG, k_target::Int)
    child = copy(x1)
    diff_idx = findall(x1 .!= x2)
    shuffle!(rng, diff_idx)
    
    swap_count = div(length(diff_idx), 2)
    for i in 1:swap_count
        child[diff_idx[i]] = x2[diff_idx[i]]
    end
    
    current_ones = sum(child)
    while current_ones > k_target
        ones_idx = findall(==(Int8(1)), child)
        child[rand(rng, ones_idx)] = 0
        current_ones -= 1
    end
    while current_ones < k_target
        zeros_idx = findall(==(Int8(0)), child)
        child[rand(rng, zeros_idx)] = 1
        current_ones += 1
    end
    return child
end

function run_minspm_bqop_solver()
    filepath = "tai256c.dat"
    if !isfile(filepath)
        error("Please place 'tai256c.dat' in the working directory.")
    end

    N, B_cpu = load_bqop_matrix(filepath)
    k_target = 92
    batch_size = 16_384
    max_generations = 5_000

    println("Initializing 64-Bit Exact MINSPM BQOP Solver [N = $N, Cardinality = $k_target]...")

    base_x = get_record_binary_seed(N, k_target)
    
    baseline_cost = Int64.(base_x)' * B_cpu * Int64.(base_x)
    println("Baseline 92-cardinality BQOP Cost: $baseline_cost")

    B_gpu = CuArray(B_cpu)
    rng = MersenneTwister(31415)

    X_host = zeros(Int8, N, batch_size)
    for b in 1:batch_size
        x = copy(base_x)
        for _ in 1:rand(rng, 1:3)
            ones_in = findall(==(Int8(1)), x)
            zeros_in = findall(==(Int8(0)), x)
            x[rand(rng, ones_in)] = 0
            x[rand(rng, zeros_in)] = 1
        end
        X_host[:, b] .= x
    end

    X_gpu = CuArray(X_host)
    energies_gpu = CUDA.zeros(Int64, batch_size)

    @cuda threads=1024 blocks=batch_size shmem=1024*sizeof(Int64) minspm_bqop_exact_kernel!(B_gpu, X_gpu, energies_gpu, N, k_target)
    CUDA.synchronize()

    E_host = Array(energies_gpu)
    global_best_cost = minimum(E_host)
    best_idx = argmin(E_host)
    global_best_x = copy(X_host[:, best_idx])

    elite_pool = [copy(global_best_x) for _ in 1:64]

    println("Executing Exact BQOP Search Loop...")

    for gen in 1:max_generations
        X_host_new = copy(X_host)
        
        for b in 1:batch_size
            if rand(rng) < 0.40
                parent2 = elite_pool[rand(rng, 1:length(elite_pool))]
                @views X_host_new[:, b] .= binary_cardinality_crossover(X_host[:, b], parent2, rng, k_target)
            else
                ones_in = findall(==(Int8(1)), X_host_new[:, b])
                zeros_in = findall(==(Int8(0)), X_host_new[:, b])
                X_host_new[rand(rng, ones_in), b] = 0
                X_host_new[rand(rng, zeros_in), b] = 1
            end
        end

        copyto!(X_gpu, X_host_new)
        energies_gpu_new = CUDA.zeros(Int64, batch_size)

        @cuda threads=1024 blocks=batch_size shmem=1024*sizeof(Int64) minspm_bqop_exact_kernel!(B_gpu, X_gpu, energies_gpu_new, N, k_target)
        CUDA.synchronize()

        E_host_new = Array(energies_gpu_new)

        for b in 1:batch_size
            if E_host_new[b] < E_host[b] && E_host_new[b] < 1000000000
                @views X_host[:, b] .= X_host_new[:, b]
                E_host[b] = E_host_new[b]
                
                if E_host[b] < global_best_cost
                    global_best_cost = E_host[b]
                    @views global_best_x .= X_host[:, b]
                    elite_pool[rand(rng, 1:length(elite_pool))] .= copy(global_best_x)
                    println("Generation $gen | NEW BQOP RECORD: $global_best_cost")
                end
            end
        end

        if gen % 300 == 0 || gen == 1
            println("Generation $gen | Current Best Cost: $global_best_cost")
        end
    end

    println("\n--- Final MINSPM BQOP Summary ---")
    println("Record Cost Achieved: $global_best_cost")

    open("solution_bqop.txt", "w") do io
        println(io, N)
        println(io, global_best_cost)
        for i in 1:N
            print(io, global_best_x[i], (i == N ? "" : " "))
        end
        println(io)
    end
    println("Saved updated 256-dim binary solution to solution_bqop.txt")
    CUDA.reclaim()
end

run_minspm_bqop_solver()
