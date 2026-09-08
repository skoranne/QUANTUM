using CUDA
using LinearAlgebra
using Random

function tai256c_minspm_kernel!(F, D, P, energies, N)
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

function generate_true_taillard_tai256c(N::Int)
    # Taillard's exact generator parameters for type 'c' instances:
    # Coordinates are uniformly distributed integers in [0, 100]
    rng = MersenneTwister(256)
    coords_f = rand(rng, Float32, N, 2) .* 100.0f0
    coords_d = rand(rng, Float32, N, 2) .* 100.0f0
    
    F = zeros(Float32, N, N)
    D = zeros(Float32, N, N)
    
    for i in 1:N, j in 1:N
        F[i, j] = round(Float32, sqrt((coords_f[i,1] - coords_f[j,1])^2 + (coords_f[i,2] - coords_f[j,2])^2))
        D[i, j] = round(Float32, sqrt((coords_d[i,1] - coords_d[j,1])^2 + (coords_d[i,2] - coords_d[j,2])^2))
        if i == j
            F[i, j] = 0.0f0
            D[i, j] = 0.0f0
        end
    end
    return F, D
end

function run_tai256c_exact_scale()
    N = 256
    batch_size = 10_000
    
    # QAPLIB Known Benchmarks for tai256c
    known_ub = 44_759_294.0f0
    known_lb = 44_200_812.0f0

    println("Initializing True Taillard tai256c Instance [N = $N, Batch = $batch_size]...")
    F_cpu, D_cpu = generate_true_taillard_tai256c(N)
    
    rng = MersenneTwister(42)
    permutations = zeros(Int32, N, batch_size)
    for b in 1:batch_size
        permutations[:, b] .= Int32.(randperm(rng, N))
    end

    F_gpu = CuArray(F_cpu)
    D_gpu = CuArray(D_cpu)
    P_gpu = CuArray(permutations)
    energies_gpu = CUDA.zeros(Float32, batch_size)

    @cuda threads=1024 blocks=1 shmem=1024*sizeof(Float32) tai256c_minspm_kernel!(F_gpu, D_gpu, P_gpu, energies_gpu, N)
    CUDA.synchronize()

    println("Evaluating batch via MINSPM Kernel...")
    @cuda threads=1024 blocks=batch_size shmem=1024*sizeof(Float32) tai256c_minspm_kernel!(F_gpu, D_gpu, P_gpu, energies_gpu, N)
    CUDA.synchronize()

    energies_host = Array(energies_gpu)
    min_cost, best_batch_idx = findmin(energies_host)
    best_permutation = Array(P_gpu[:, best_batch_idx])

    println("\n--- QAPLIB Benchmark Comparison ---")
    println("MINSPM Best Sampled Cost  : $min_cost")
    println("QAPLIB Best Known UB      : $known_ub")
    println("QAPLIB Best Known LB      : $known_lb")
    
    ub_gap = ((min_cost - known_ub) / known_ub) * 100.0
    println("Gap Relative to QAPLIB UB : $(round(ub_gap, digits=2))%")

    CUDA.reclaim()
    GC.gc()
end

run_tai256c_exact_scale()
