using CUDA
using LinearAlgebra
using BenchmarkTools
using Random

# Top-level declarations for global scope visibility
function eval_dense_qubo(Q, X)
    return sum(X .* (Q * X), dims=1)
end

function minspm_kernel!(F, D, P, energies, N)
    batch_id = blockIdx().x
    tid = threadIdx().x
    
    shared_vals = @cuDynamicSharedMem(Float32, 1024)
    
    i = div(tid - 1, N) + 1
    j = rem(tid - 1, N) + 1
    
    p_i = P[i, batch_id]
    p_j = P[j, batch_id]
    
    shared_vals[tid] = F[i, j] * D[p_i, p_j]
    sync_threads()
    
    step = 512
    while step > 0
        if tid <= step
            shared_vals[tid] += shared_vals[tid + step]
        end
        sync_threads()
        step = step >> 1
    end
    
    if tid == 1
        energies[batch_id] = shared_vals[1]
    end
    return nothing
end

function eval_minspm(F, D, P, energies, N, batch_size)
    @cuda threads=1024 blocks=batch_size shmem=1024*sizeof(Float32) minspm_kernel!(F, D, P, energies, N)
    return energies
end

function run_qaplib_nug32_benchmark()
    N = 32
    batch_size = 50_000 
    
    rng = MersenneTwister(1337)
    F_cpu = abs.(rand(rng, Float32, N, N) .* 15.0f0)
    D_cpu = abs.(rand(rng, Float32, N, N) .* 100.0f0)
    
    for i in 1:N
        F_cpu[i,i] = 0.0f0
        D_cpu[i,i] = 0.0f0
    end

    permutations = zeros(Int32, N, batch_size)
    for b in 1:batch_size
        permutations[:, b] .= Int32.(randperm(rng, N))
    end

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

    F_gpu = CuArray(F_cpu)
    D_gpu = CuArray(D_cpu)
    P_gpu = CuArray(permutations)
    energies_gpu = CUDA.zeros(Float32, batch_size)

    # Warmup compilations
    eval_dense_qubo(Q_gpu, X_dense_gpu)
    CUDA.@sync eval_minspm(F_gpu, D_gpu, P_gpu, energies_gpu, N, batch_size)

    println("Running QAPLIB (Nug32-archetype) Benchmark [Batch = 50,000]:\n")
    
    println("cuBLAS Dense QUBO Backend:")
    display(@benchmark eval_dense_qubo($Q_gpu, $X_dense_gpu))

    println("\nMINSPM Galois Custom Kernel:")
    display(@benchmark CUDA.@sync eval_minspm($F_gpu, $D_gpu, $P_gpu, $energies_gpu, $N, $batch_size))
end

run_qaplib_nug32_benchmark()
