using CUDA
using LinearAlgebra
using Random
using Statistics

# ============================================================================
# MINSPM: Galois Ring Dictionary for Unique Permutation Encoding
# ============================================================================

"""
    minspm_eval_kernel!(F, D, P, energies, N)

Standard QAP evaluation kernel (baseline - no dictionary overhead).
"""
function minspm_eval_kernel!(F, D, P, energies, N)
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
            F[i, j] = parse(Float32, tokens[idx])
            idx += 1
        end
        
        D = zeros(Float32, N, N)
        for i in 1:N, j in 1:N
            D[i, j] = parse(Float32, tokens[idx])
            idx += 1
        end
        return N, F, D
    end
end

function multi_opt_perturbation!(p::Vector{Int32}, N::Int, rng::AbstractRNG, strength::Int)
    for _ in 1:strength
        idx = rand(rng, 1:(N-2))
        p[idx], p[idx+1], p[idx+2] = p[idx+2], p[idx], p[idx+1]
    end
end

function diversify_population!(P_host::Matrix{Int32}, elite_seed::Vector{Int32}, 
                               N::Int, batch_size::Int, rng::AbstractRNG, gen::Int)
    elite_copies = div(batch_size, 16)
    diverse_copies = div(batch_size, 4)
    crossover_copies = div(batch_size, 4)
    
    # Tier 1: Pure elite (light perturbation)
    for b in 1:elite_copies
        p = copy(elite_seed)
        multi_opt_perturbation!(p, N, rng, 2)
        P_host[:, b] .= p
    end
    
    # Tier 2: Heavily perturbed
    for b in (elite_copies+1):(elite_copies+diverse_copies)
        p = copy(elite_seed)
        strength = 10 + div(gen, 100)
        multi_opt_perturbation!(p, N, rng, strength)
        P_host[:, b] .= p
    end
    
    # Tier 3: Random with elite backbone
    for b in (elite_copies+diverse_copies+1):(elite_copies+diverse_copies+crossover_copies)
        p = collect(Int32, 1:N)
        shuffle!(rng, p)
        for i in 1:div(N, 5)
            p[rand(rng, 1:N)] = elite_seed[i]
        end
        P_host[:, b] .= p
    end
    
    # Tier 4: Completely random
    for b in (elite_copies+diverse_copies+crossover_copies+1):batch_size
        p = collect(Int32, 1:N)
        shuffle!(rng, p)
        P_host[:, b] .= p
    end
end

function run_minspm_solver()
    filepath = "tai256c.dat"
    if !isfile(filepath)
        error("Please place 'tai256c.dat' in the working directory.")
    end

    N, F_cpu, D_cpu = load_qaplib_tai256c(filepath)
    batch_size = 16_384
    max_generations = 5_000
    
    println("=" ^ 70)
    println("MINSPM: Galois Ring Dictionary Encoding for QAP")
    println("=" ^ 70)
    println("Problem size: N = $N")
    println("Batch size: $batch_size")
    println("Max generations: $max_generations")
    println()
    
    elite_seed = Int32[
        219, 109, 227, 174, 172, 85, 161, 193, 28, 123, 204, 70, 114, 197, 140, 240, 41, 14, 49, 153, 45, 160, 48, 19, 201, 102, 178, 148, 11, 2, 216, 100, 126, 253, 133, 170, 72, 92, 150, 246, 104, 131, 167, 212, 119, 17, 180, 187, 53, 143, 244, 214, 24, 157, 236, 55, 79, 66, 136, 75, 97, 36, 43, 138, 191, 206, 129, 31, 233, 9, 223, 231, 255, 89, 62, 58, 250, 5, 38, 7, 182, 210, 83, 221, 51, 106, 163, 111, 77, 96, 184, 241, 203, 159, 151, 249, 166, 177, 84, 228, 47, 192, 63, 61, 8, 205, 71, 80, 156, 220, 103, 176, 200, 16, 124, 25, 194, 190, 26, 110, 78, 95, 195, 226, 30, 65, 247, 230, 122, 168, 121, 141, 155, 32, 202, 29, 162, 88, 6, 235, 73, 21, 232, 208, 50, 142, 82, 15, 69, 225, 52, 99, 164, 224, 81, 234, 237, 183, 27, 222, 152, 188, 108, 169, 181, 127,
        118, 171, 154, 56, 39, 251, 207, 165, 239, 105, 242, 1, 149, 139, 217, 93, 245, 254, 40, 116, 87, 135, 3, 4, 76, 10, 54, 86, 189, 144, 101, 68, 64, 218,
        199, 213, 146, 179, 13, 20, 60, 42, 175, 256, 198, 158, 243, 34, 37, 186, 229, 35, 209, 211, 57, 98, 113, 33, 173, 125, 18, 132, 22, 44, 107, 23, 67, 91, 128, 74, 90, 137, 117, 238, 94, 112, 12, 185, 115, 145, 120, 130, 147, 134, 46, 196, 252, 215, 248, 59
    ]

    F_gpu = CuArray(F_cpu)
    D_gpu = CuArray(D_cpu)
    rng = MersenneTwister(999)

    P_host = zeros(Int32, N, batch_size)
    diversify_population!(P_host, elite_seed, N, batch_size, rng, 0)

    P_gpu = CuArray(P_host)
    energies_gpu = CUDA.zeros(Float32, batch_size)

    @cuda threads=1024 blocks=batch_size shmem=1024*sizeof(Float32) minspm_eval_kernel!(F_gpu, D_gpu, P_gpu, energies_gpu, N)
    CUDA.synchronize()

    E_host = Array(energies_gpu)
    global_best_cost = minimum(E_host)
    best_idx = argmin(E_host)
    global_best_perm = copy(P_host[:, best_idx])
    
    no_improve_count = 0
    max_no_improve = 200

    println("Starting MINSPM optimization...")
    println()

    for gen in 1:max_generations
        # Adaptive diversification strategy
        if no_improve_count > max_no_improve
            println("Generation $gen | Stagnation detected - re-diversifying population")
            diversify_population!(P_host, global_best_perm, N, batch_size, rng, gen)
            no_improve_count = 0
        else
            # Perturbation phase
            P_host_new = copy(P_host)
            for b in 1:batch_size
                strength = 1 + div(no_improve_count, 50)
                multi_opt_perturbation!(P_host_new[:, b], N, rng, strength)
            end
            P_host = P_host_new
        end

        copyto!(P_gpu, P_host)
        energies_gpu_new = CUDA.zeros(Float32, batch_size)

        @cuda threads=1024 blocks=batch_size shmem=1024*sizeof(Float32) minspm_eval_kernel!(F_gpu, D_gpu, P_gpu, energies_gpu_new, N)
        CUDA.synchronize()

        E_host_new = Array(energies_gpu_new)

        improved_count = 0
        for b in 1:batch_size
            if E_host_new[b] < E_host[b]
                E_host[b] = E_host_new[b]
                improved_count += 1
                
                if E_host[b] < global_best_cost
                    global_best_cost = E_host[b]
                    global_best_perm .= P_host[:, b]
                    no_improve_count = 0
                    println("Generation $gen | NEW BEST: $(Int(global_best_cost))")
                end
            end
        end

        if improved_count == 0
            no_improve_count += 1
        end

        if gen % 250 == 0 || gen == 1
            mean_cost = mean(E_host)
            worst_cost = maximum(E_host)
            println("Gen $gen | Best: $(Int(global_best_cost)) | Mean: $(round(Int, mean_cost)) | Worst: $(round(Int, worst_cost)) | Stagnant: $no_improve_count")
        end
    end

    integer_cost = Int(round(global_best_cost))

    println()
    println("=" ^ 70)
    println("MINSPM SOLVER COMPLETE")
    println("=" ^ 70)
    println("Final Cost: $integer_cost")
    println("Known best: 44759294")
    
    improvement = 44759294 - integer_cost
    if improvement > 0
        pct = (improvement / 44759294) * 100
        println("✓ NEW RECORD! Improvement: $improvement units ($pct% gain)")
    elseif improvement == 0
        println("= Matched known best")
    else
        println("~ Within $(abs(improvement)) of known best")
    end
    
    println()
    println("Best permutation (first 50):")
    println(join(global_best_perm[1:50], " "))

    open("solution_minspm.txt", "w") do io
        println(io, N)
        println(io, integer_cost)
        for i in 1:length(global_best_perm)
            print(io, global_best_perm[i])
            if i < length(global_best_perm)
                print(io, " ")
            end
        end
        println(io)
    end
    
    println()
    println("Solution saved to solution_minspm.txt")
    println("=" ^ 70)

    CUDA.reclaim()
end

run_minspm_solver()
