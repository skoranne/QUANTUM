using CUDA
using LinearAlgebra
using JuMP
using BenchmarkTools
using Random

function run_qap_benchmark(N::Int=32, batch_size::Int=10_000)
    println("Initializing MINSPM Galois Dictionary vs Standard QUBO Benchmark")
    println("Problem Size: N = $N, Batch Size = $batch_size permutations\n")

    # 1. Generate Problem Data (Random Flow and Distance matrices)
    F_cpu = rand(Float32, N, N)
    D_cpu = rand(Float32, N, N)
    
    # Generate random permutation states for the batch evaluation
    permutations = zeros(Int32, N, batch_size)
    for b in 1:batch_size
        permutations[:, b] .= Int32.(randperm(N))
    end

    # ==============================================================
    # JUMP FORMULATION (The cuOpt / Classical Solver Logical Model)
    # ==============================================================
    println("--- JuMP Formulation (Logical Overhead) ---")
    model = Model()
    @variable(model, x[1:N, 1:N], Bin)
    
    # Standard QUBO constraints (One-hot rows and columns)
    @constraint(model, [i=1:N], sum(x[i, :]) == 1)
    @constraint(model, [j=1:N], sum(x[:, j]) == 1)
    
    num_vars = N^2
    println("Variables generated: $num_vars binary variables.")
    println("QUBO Matrix size required for backend: $(num_vars) x $(num_vars)\n")

    # ==============================================================
    # DENSE QUBO (cuBLAS / cuTENSOR BACKEND)
    # ==============================================================
    println("--- Standard Dense QUBO (cuBLAS Backend) ---")
    
    # Build Dense Q = kron(F, D)
    Q_cpu = kron(F_cpu, D_cpu)
    Q_gpu = CuArray(Q_cpu)
    
    # Convert permutation dictionary to dense one-hot flattened vectors
    X_dense_cpu = zeros(Float32, num_vars, batch_size)
    for b in 1:batch_size
        for i in 1:N
            X_dense_cpu[(i-1)*N + permutations[i, b], b] = 1.0f0
        end
    end
    X_dense_gpu = CuArray(X_dense_cpu)
    
    # Evaluate batch energy: diag(X^T * Q * X) -> sum(X .* (Q * X))
    function eval_dense_qubo(Q, X)
        QX = Q * X  # cuBLAS SGEMM matrix multiplication
        return sum(X .* QX, dims=1)
    end
    
    qubo_mem_mb = sizeof(Q_cpu) / (1024^2)
    println("Dense QUBO Matrix Memory: $(round(qubo_mem_mb, digits=2)) MB\n")

    # ==============================================================
    # MINSPM GALOIS DICTIONARY ENCODING (Custom CUDA Kernel)
    # ==============================================================
    println("--- MINSPM Galois Dictionary Encoding ---")
    F_gpu = CuArray(F_cpu)
    D_gpu = CuArray(D_cpu)
    P_gpu = CuArray(permutations)
    energies_gpu = CUDA.zeros(Float32, batch_size)
    
    # The dictionary encoding bypassing dense matrices using cache residency
    function minspm_kernel!(F, D, P, energies, N)
        batch_id = blockIdx().x
        tid = threadIdx().x
        
        # Dynamic shared memory for block-level reduction
        shared_vals = @cuDynamicSharedMem(Float32, 1024)
        
        # Map 1D thread (1 to 1024) to 2D tensor space (32x32)
        i = div(tid - 1, N) + 1
        j = rem(tid - 1, N) + 1
        
        # Dictionary Encoding Lookup
        p_i = P[i, batch_id]
        p_j = P[j, batch_id]
        
        # Tensor contraction directly into shared memory
        shared_vals[tid] = F[i, j] * D[p_i, p_j]
        sync_threads()
        
        # Parallel Tree Reduction
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
        threads = 1024 # Fits N=32 perfectly (32^2 = 1024 threads per block)
        blocks = batch_size
        shmem = 1024 * sizeof(Float32)
        @cuda threads=threads blocks=blocks shmem=shmem minspm_kernel!(F, D, P, energies, N)
        return energies
    end

    dict_mem_kb = (sizeof(F_cpu) + sizeof(D_cpu)) / 1024
    println("MINSPM Dictionary Memory: $(round(dict_mem_kb, digits=2)) KB\n")

    # ==============================================================
    # BENCHMARKING (10,000 parallel state evaluations)
    # ==============================================================
    println("Running Benchmarks...")
    
    # Warmup
    eval_dense_qubo(Q_gpu, X_dense_gpu)
    eval_minspm(F_gpu, D_gpu, P_gpu, energies_gpu, N, batch_size)
    
    println("cuBLAS Dense QUBO Backend:")
    display(@benchmark $eval_dense_qubo($Q_gpu, $X_dense_gpu))
    
    println("\nMINSPM Galois Custom Kernel:")
    display(@benchmark CUDA.@sync $eval_minspm($F_gpu, $D_gpu, $P_gpu, $energies_gpu, $N, $batch_size))
end

run_qap_benchmark(32, 10_000)
