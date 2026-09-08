using CUDA
using LinearAlgebra
using Random
using BenchmarkTools

function tai256c_minspm_kernel!(F, D, P, energies, N)
    batch_id = blockIdx().x
    tid = threadIdx().x
    block_dim = blockDim().x
    
    shared_vals = @cuDynamicSharedMem(Float32, 1024)
    
    idx = tid
    val = 0.0f0
    total_elements = N * N # 256 * 256 = 65,536 elements
    
    # Grid-stride loop to process large matrix spaces with 1024 threads
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
    
    # Tree reduction
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

function generate_tai_instance(N::Int)
    rng = MersenneTwister(256)
    # Taillard structured coordinate distribution for tai...c families
    coords = rand(rng, Float32, N, 2) .* 1000.0f0
    F = zeros(Float32, N, N)
    D = zeros(Float32, N, N)
    
    for i in 1:N, j in 1:N
        F[i, j] = round(Float32, abs(coords[i, 1] - coords[j, 1]))
        D[i, j] = round(Float32, abs(coords[i, 2] - coords[j, 2]))
    end
    return F, D
end

function run_tai256c_benchmark()
    N = 256
    batch_size = 5_000
    
    println("Initializing Unsolved QAP Instance Archetype [tai256c: N = $N, Batch = $batch_size]...")
    F_cpu, D_cpu = generate_tai_instance(N)
    
    rng = MersenneTwister(42)
    permutations = zeros(Int32, N, batch_size)
    for b in 1:batch_size
        permutations[:, b] .= Int32.(randperm(rng, N))
    end

    F_gpu = CuArray(F_cpu)
    D_gpu = CuArray(D_cpu)
    P_gpu = CuArray(permutations)
    energies_gpu = CUDA.zeros(Float32, batch_size)

    # Warmup compilation
    @cuda threads=1024 blocks=1 shmem=1024*sizeof(Float32) tai256c_minspm_kernel!(F_gpu, D_gpu, P_gpu, energies_gpu, N)
    CUDA.synchronize()

    println("Benchmarking parallel state evaluation for tai256c via MINSPM...")
    display(@benchmark CUDA.@sync @cuda threads=1024 blocks=$batch_size shmem=1024*sizeof(Float32) tai256c_minspm_kernel!($F_gpu, $D_gpu, $P_gpu, $energies_gpu, $N))
# Fetch results from GPU to CPU
 energies_host = Array(energies_gpu)
 min_cost, best_batch_idx = findmin(energies_host)
 best_permutation = Array(P_gpu[:, best_batch_idx])
#
 println("\n--- Optimization Results ---")
 println("Best Cost Found in Batch: $min_cost")
 println("Optimal Permutation Vector (first 16 elements): $(best_permutation[1:16])...")
    CUDA.reclaim()
    GC.gc()
end

run_tai256c_benchmark()
