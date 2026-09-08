using CUDA
using LinearAlgebra
using Random

function esc128_minspm_kernel!(F, D, P, energies, N)
    batch_id = blockIdx().x
    tid = threadIdx().x
    block_dim = blockDim().x
    
    shared_vals = @cuDynamicSharedMem(Float32, 1024)
    
    idx = tid
    val = 0.0f0
    while idx <= (N * N)
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

function generate_esc_instance(N::Int)
    F = zeros(Float32, N, N)
    D = zeros(Float32, N, N)
    
    for i in 1:N, j in 1:N
        # Escallon (esc) instances are characterized by Hamming distance graph structures
        F[i, j] = Float32(count_ones((i - 1) ⊻ (j - 1)))
        D[i, j] = Float32((i == j) ? 0.0f0 : 1.0f0 / (abs(i - j) + 1.0f0))
    end
    return F, D
end

function run_esc128_benchmark()
    N = 128
    batch_size = 10_000
    
    println("Generating 'esc128' QAP Benchmark Structure [N = $N]...")
    F_cpu, D_cpu = generate_esc_instance(N)
    
    rng = MersenneTwister(1337)
    permutations = zeros(Int32, N, batch_size)
    for b in 1:batch_size
        permutations[:, b] .= Int32.(randperm(rng, N))
    end

    F_gpu = CuArray(F_cpu)
    D_gpu = CuArray(D_cpu)
    P_gpu = CuArray(permutations)
    energies_gpu = CUDA.zeros(Float32, batch_size)

    # Warmup compilation
    @cuda threads=1024 blocks=1 shmem=1024*sizeof(Float32) esc128_minspm_kernel!(F_gpu, D_gpu, P_gpu, energies_gpu, N)
    CUDA.synchronize()

    println("Evaluating batch of $batch_size permutations for esc128 via MINSPM...")
    t_start = time()
    @cuda threads=1024 blocks=batch_size shmem=1024*sizeof(Float32) esc128_minspm_kernel!(F_gpu, D_gpu, P_gpu, energies_gpu, N)
    CUDA.synchronize()
    t_elapsed = time() - t_start

    energies_host = Array(energies_gpu)
    min_cost, best_idx = findmin(energies_host)

    println("\nesc128 MINSPM Evaluation Complete:")
    println("Execution Time for Batch: $(round(t_elapsed * 1000, digits=2)) ms")
    println("Best Sampled Cost: $min_cost")

    CUDA.reclaim()
    GC.gc()
end

run_esc128_benchmark()
