using CUDA
using LinearAlgebra
using Random

function minspm_qap_eval_kernel!(F, D, P, energies, N)
    batch_id = blockIdx().x
    tid = threadIdx().x
    block_dim = blockDim().x
    
    shared_vals = @cuDynamicSharedMem(Float32, 1024)
    
    idx = tid
    val = 0.0f0
    total_elements = N * N
    
    while idx <= total_elements
        i = div(idx - 1, N) + 1
        j = rem(idx - 1, N) + 1
        
        p_i = P[i, batch_id]
        p_j = P[j, batch_id]
        
        val += F[i, j] * D[p_i, p_j]
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
        energies[batch_id] = shared_vals[1]
    end
    return nothing
end

function load_qaplib_tai256c(filepath::String)
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

function load_or_generate_seed(N::Int, F_cpu::Matrix{Float32})
    if isfile("solution.txt")
        try
            tokens = split(read("solution.txt", String))
            if length(tokens) >= 2 + N
                perm = parse.(Int32, tokens[3:(2+N)])
                if length(perm) == N
                    println("Successfully loaded 256-element seed from solution.txt")
                    return perm
                end
            end
        catch
        end
    end
    println("Falling back to affinity-guided sequence initialization...")
    return Int32.(sortperm(vec(sum(F_cpu, dims=2))))
end

function pmx_crossover(p1::AbstractVector{Int32}, p2::AbstractVector{Int32}, rng::AbstractRNG)
    N = length(p1)
    cx1, cx2 = sort(rand(rng, 1:N, 2))
    offspring = copy(p1)
    mapping = zeros(Int32, N + 1)
    
    for i in cx1:cx2
        offspring[i] = p2[i]
        mapping[p2[i]] = p1[i]
    end
    
    for i in 1:N
        if i < cx1 || i > cx2
            val = p1[i]
            while mapping[val] != 0
                val = mapping[val]
            end
            offspring[i] = val
        end
    end
    return offspring
end

function run_elite_alns_solver()
    filepath = "tai256c.dat"
    if !isfile(filepath)
        error("Please place 'tai256c.dat' in the working directory.")
    end

    N, F_cpu, D_cpu = load_qaplib_tai256c(filepath)
    batch_size = 16_384
    max_generations = 5_000

    record_seed = load_or_generate_seed(N, F_cpu)

    F_gpu = CuArray(F_cpu)
    D_gpu = CuArray(D_cpu)
    rng = MersenneTwister(1337)

    P_host = zeros(Int32, N, batch_size)
    for b in 1:batch_size
        p = copy(record_seed)
        for _ in 1:rand(rng, 1:3)
            i1, i2 = rand(rng, 1:N), rand(rng, 1:N)
            p[i1], p[i2] = p[i2], p[i1]
        end
        P_host[:, b] .= p
    end

    P_gpu = CuArray(P_host)
    energies_gpu = CUDA.zeros(Float32, batch_size)

    @cuda threads=1024 blocks=batch_size shmem=1024*sizeof(Float32) minspm_qap_eval_kernel!(F_gpu, D_gpu, P_gpu, energies_gpu, N)
    CUDA.synchronize()

    E_host = Array(energies_gpu)
    global_best_cost = minimum(E_host)
    best_idx = argmin(E_host)
    global_best_perm = copy(P_host[:, best_idx])

    elite_pool = [copy(global_best_perm) for _ in 1:64]

    println("Initializing Elite-Anchored ALNS Solver...")

    for gen in 1:max_generations
        P_host_new = copy(P_host)
        
        for b in 1:batch_size
            if rand(rng) < 0.40
                parent2 = elite_pool[rand(rng, 1:length(elite_pool))]
                @views P_host_new[:, b] .= pmx_crossover(P_host[:, b], parent2, rng)
            else
                start_idx = rand(rng, 1:(N - 16))
                @views sub_block = P_host_new[start_idx:(start_idx + 16), b]
                shuffle!(rng, sub_block)
            end
        end

        copyto!(P_gpu, P_host_new)
        energies_gpu_new = CUDA.zeros(Float32, batch_size)

        @cuda threads=1024 blocks=batch_size shmem=1024*sizeof(Float32) minspm_qap_eval_kernel!(F_gpu, D_gpu, P_gpu, energies_gpu_new, N)
        CUDA.synchronize()

        E_host_new = Array(energies_gpu_new)

        for b in 1:batch_size
            if E_host_new[b] < E_host[b]
                @views P_host[:, b] .= P_host_new[:, b]
                E_host[b] = E_host_new[b]
                if E_host[b] < global_best_cost
                    global_best_cost = E_host[b]
                    @views global_best_perm .= P_host[:, b]
                    elite_pool[rand(rng, 1:length(elite_pool))] .= copy(global_best_perm)
                    println("Generation $gen | NEW RECORD: $global_best_cost")
                end
            end
        end

        if gen % 300 == 0 || gen == 1
            println("Generation $gen | Current Best Cost: $global_best_cost")
        end
    end

    integer_cost = Int(round(global_best_cost))
    println("\n--- Final Optimization Summary ---")
    println("Record Cost Achieved: $integer_cost")

    open("solution.txt", "w") do io
        println(io, N)
        println(io, integer_cost)
        for i in 1:length(global_best_perm)
            print(io, global_best_perm[i], (i == length(global_best_perm) ? "" : " "))
        end
        println(io)
    end
    println("Saved updated solution to solution.txt")
    CUDA.reclaim()
end

run_elite_alns_solver()
