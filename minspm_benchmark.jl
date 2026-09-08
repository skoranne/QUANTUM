using CUDA
using LinearAlgebra
using Random
using Printf

# =====================================================================
# 1. BASELINE: DENSE NON-MINSPM QAP EVALUATOR (Float32)
# =====================================================================
function dense_qap_eval_kernel!(F, D, P, energies, N)
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
        
        # Heavy memory indirection and floating point math
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

# =====================================================================
# 2. PROPOSED: EXACT MINSPM BQOP ORACLE (Int64)
# =====================================================================
function minspm_bqop_exact_kernel!(B_dict, X_batch, energies, N)
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
        
        # O(1) Galois Dictionary Lookup - Zero indirection, exact integer math
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
        energies[batch_id] = shared_vals[1]
    end
    return nothing
end

# =====================================================================
# 3. BENCHMARK EXECUTION SUITE
# =====================================================================
function run_ablation_benchmark()
    N = 256
    batch_size = 65_536  # Massive parallel batch for saturation
    runs = 20
    rng = MersenneTwister(1234)

    println("=========================================================")
    println(" MINSPM Oracle vs Dense QAP Hardware Ablation Benchmark")
    println("=========================================================")
    @printf(" Problem Size (N) : %d\n", N)
    @printf(" Batch Size       : %d\n", batch_size)
    @printf(" Benchmark Runs   : %d\n\n", runs)

    # --- Setup Dense Baseline Data ---
    F_host = rand(Float32, N, N)
    D_host = rand(Float32, N, N)
    P_host = zeros(Int32, N, batch_size)
    for b in 1:batch_size
        P_host[:, b] .= randperm(rng, N)
    end
    F_gpu = CuArray(F_host)
    D_gpu = CuArray(D_host)
    P_gpu = CuArray(P_host)
    energies_dense_gpu = CUDA.zeros(Float32, batch_size)

    # --- Setup MINSPM Oracle Data ---
    B_host = rand(Int64, N, N)
    X_host = zeros(Int8, N, batch_size)
    for b in 1:batch_size
        # Generate 92-cardinality vectors
        idx = randperm(rng, N)[1:92]
        X_host[idx, b] .= 1
    end
    B_gpu = CuArray(B_host)
    X_gpu = CuArray(X_host)
    energies_minspm_gpu = CUDA.zeros(Int64, batch_size)

    # --- WARMUP ---
    @cuda threads=1024 blocks=batch_size shmem=1024*sizeof(Float32) dense_qap_eval_kernel!(F_gpu, D_gpu, P_gpu, energies_dense_gpu, N)
    @cuda threads=1024 blocks=batch_size shmem=1024*sizeof(Int64) minspm_bqop_exact_kernel!(B_gpu, X_gpu, energies_minspm_gpu, N)
    CUDA.synchronize()

    # --- BENCHMARK DENSE ---
    println("Benchmarking Dense Non-MINSPM QAP Kernel (Float32)...")
    dense_times = zeros(Float64, runs)
    for r in 1:runs
        dense_times[r] = CUDA.@elapsed begin
            @cuda threads=1024 blocks=batch_size shmem=1024*sizeof(Float32) dense_qap_eval_kernel!(F_gpu, D_gpu, P_gpu, energies_dense_gpu, N)
        end
    end
    dense_avg_time = sum(dense_times) / runs

    # --- BENCHMARK MINSPM ---
    println("Benchmarking Exact MINSPM BQOP Kernel (Int64)...")
    minspm_times = zeros(Float64, runs)
    for r in 1:runs
        minspm_times[r] = CUDA.@elapsed begin
            @cuda threads=1024 blocks=batch_size shmem=1024*sizeof(Int64) minspm_bqop_exact_kernel!(B_gpu, X_gpu, energies_minspm_gpu, N)
        end
    end
    minspm_avg_time = sum(minspm_times) / runs

    # --- REPORTING ---
    dense_throughput = batch_size / dense_avg_time / 1e6
    minspm_throughput = batch_size / minspm_avg_time / 1e6
    speedup = dense_avg_time / minspm_avg_time

    println("\n=========================================================")
    println(" ABLATION RESULTS")
    println("=========================================================")
    @printf(" Baseline (Dense) Time : %.4f ms per batch\n", dense_avg_time * 1000)
    @printf(" MINSPM (Exact) Time   : %.4f ms per batch\n", minspm_avg_time * 1000)
    println("---------------------------------------------------------")
    @printf(" Baseline Throughput   : %.2f Million Evals / sec\n", dense_throughput)
    @printf(" MINSPM Throughput     : %.2f Million Evals / sec\n", minspm_throughput)
    println("---------------------------------------------------------")
    @printf(" MINSPM Speedup        : %.2fx Faster\n", speedup)
    println("=========================================================")

    CUDA.reclaim()
end

run_ablation_benchmark()
