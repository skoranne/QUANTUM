using CUDA
using LinearAlgebra
using BenchmarkTools
using Random

# Top-level declarations for global scope visibility
function eval_dense_qubo_n64(Q, X)
    return sum(X .* (Q * X), dims=1)
end

function minspm_kernel_n64!(F, D, P, energies, N)
    batch_id = blockIdx().x
    tid = threadIdx().x
    block_dim = blockDim().x
    
    # Dynamic shared memory for parallel reduction of 4096 elements
    shared_vals = @cuDynamicSharedMem(Float32, 4096)
    
    # Grid-stride loop to handle N^2 = 4096 elements with 1024 threads
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
    
    # Tree reduction in shared memory (reducing 1024 down to 1)
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

function eval_minspm_n64(F, D, P, energies, N, batch_size)
    # 1024 threads per block, 1 block per batch instance
    @cuda threads=1024 blocks=batch_size shmem=4096*sizeof(Float32) minspm_kernel_n64!(F, D, P, energies, N)
    return energies
end

function run_large_qap_benchmark()
    N = 64
    batch_size = 20_000 # Large scale evaluation batch
    
    println("Initializing Large-Scale QAP Benchmark [N = $N, Batch = $batch_size]...")
    
    rng = MersenneTwister(42)
    F_cpu = abs.(rand(rng, Float32, N, N) .* 20.0f0)
    D_cpu = abs.(rand(rng, Float32, N, N) .* 150.0f0)
    
    for i in 1:N
        F_cpu[i,i] = 0.0f0
        D_cpu[i,i] = 0.0f0
    end

    permutations = zeros(Int32, N, batch_size)
    for b in 1:batch_size
        permutations[:, b] .= Int32.(randperm(rng, N))
    end

    # 1. cuBLAS Dense QUBO Setup (The Memory Wall)
    Q_cpu = kron(F_cpu, D_cpu)
    Q_gpu = CuArray(Q_cpu)
    
    num_vars = N^2
    X_dense_cpu = zeros(Float32, num_vars, batch_size)
    for b in 1:batch_size
        for i in 1:N
            X_dense_cpu[(i-1)*N + permutations[i, b], b] = 1.0f0
        end
    end
    X_dense_gpu = CuArray(X_dense_cpu)

    # 2. MINSPM Galois Kernel Setup
    F_gpu = CuArray(F_cpu)
    D_gpu = CuArray(D_cpu)
    P_gpu = CuArray(permutations)
    energies_gpu = CUDA.zeros(Float32, batch_size)

    # Warmup compilations
    eval_dense_qubo_n64(Q_gpu, X_dense_gpu)
    CUDA.@sync eval_minspm_n64(F_gpu, D_gpu, P_gpu, energies_gpu, N, batch_size)

    println("\nBenchmarking execution speeds...")
    
    println("cuBLAS Dense QUBO Backend (N=64):")
    display(@benchmark eval_dense_qubo_n64($Q_gpu, $X_dense_gpu))

    println("\nMINSPM Galois Custom Kernel (N=64):")
    display(@benchmark CUDA.@sync eval_minspm_n64($F_gpu, $D_gpu, $P_gpu, $energies_gpu, $N, $batch_size))

    CUDA.reclaim()
    GC.gc()
end

run_large_qap_benchmark()
