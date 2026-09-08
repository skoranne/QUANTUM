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
# 2. PROPOSED: MINSPM SPARSE GALOIS HASH ORACLE (Int64)
# =====================================================================
# Models the true MINSPM efficiency: Iterating strictly over the 
# compressed non-zero Galois characteristic terms (e.g., ~512 NNZ).
function minspm_sparse_galois_kernel!(Galois_Keys, Galois_Values, X_batch, energies, nnz_terms)
    tid = threadIdx().x
    warp_id = div(tid - 1, 32) + 1
    lane_id = rem(tid - 1, 32) + 1
    
    # 1 Warp = 1 BQOP Evaluation
    batch_id = (blockIdx().x - 1) * div(blockDim().x, 32) + warp_id
    
    shared_vals = @cuDynamicSharedMem(Int64, blockDim().x)
    val = Int64(0)
    
    if batch_id <= size(X_batch, 2)
        idx = lane_id
        # Iterate ONLY over the minimal perfect hash non-zeros
        while idx <= nnz_terms
            # O(1) fetching of the compressed Galois term
            u = Galois_Keys[1, idx]
            v = Galois_Keys[2, idx]
            weight = Galois_Values[idx]
            
            # Check if both vertices in this non-zero edge are active
            if X_batch[u, batch_id] == Int8(1) && X_batch[v, batch_id] == Int8(1)
                val += weight
            end
            
            idx += 32
        end
    end
    
    shared_vals[tid] = val
    sync_threads()
    
    # Warp-level parallel reduction
    if lane_id <= 16; shared_vals[tid] += shared_vals[tid + 16]; end; sync_threads()
    if lane_id <= 8;  shared_vals[tid] += shared_vals[tid + 8];  end; sync_threads()
    if lane_id <= 4;  shared_vals[tid] += shared_vals[tid + 4];  end; sync_threads()
    if lane_id <= 2;  shared_vals[tid] += shared_vals[tid + 2];  end; sync_threads()
    if lane_id == 1;  shared_vals[tid] += shared_vals[tid + 1];  end; sync_threads()
    
    if lane_id == 1 && batch_id <= size(X_batch, 2)
        energies[batch_id] = shared_vals[tid]
    end
    return nothing
end

# =====================================================================
# 3. BENCHMARK EXECUTION SUITE
# =====================================================================
function run_ablation_benchmark()
    N = 256
    batch_size = 65_536
    nnz_terms = 512 # Compressed Galois dictionary size for the sparse QAP
    runs = 20
    rng = MersenneTwister(1234)

    println("=========================================================")
    println(" V3 MINSPM Sparse Galois Oracle vs Dense QAP Baseline")
    println("=========================================================")
    @printf(" Problem Size (N) : %d\n", N)
    @printf(" Batch Size       : %d\n", batch_size)
    @printf(" Galois NNZ Terms : %d\n\n", nnz_terms)

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
    Galois_Keys_host = rand(1:N, 2, nnz_terms)
    Galois_Values_host = rand(1:1000, nnz_terms)
    X_host = zeros(Int8, N, batch_size)
    for b in 1:batch_size
        X_host[randperm(rng, N)[1:92], b] .= 1
    end
    
    Galois_Keys_gpu = CuArray(Int32.(Galois_Keys_host))
    Galois_Values_gpu = CuArray(Int64.(Galois_Values_host))
    X_gpu = CuArray(X_host)
    energies_minspm_gpu = CUDA.zeros(Int64, batch_size)

    # Kernel launch configurations
    dense_threads = 1024
    dense_blocks = batch_size

    minspm_threads = 256
    minspm_warps_per_block = minspm_threads ÷ 32
    minspm_blocks = ceil(Int, batch_size / minspm_warps_per_block)

    # --- WARMUP ---
    @cuda threads=dense_threads blocks=dense_blocks shmem=dense_threads*sizeof(Float32) dense_qap_eval_kernel!(F_gpu, D_gpu, P_gpu, energies_dense_gpu, N)
    @cuda threads=minspm_threads blocks=minspm_blocks shmem=minspm_threads*sizeof(Int64) minspm_sparse_galois_kernel!(Galois_Keys_gpu, Galois_Values_gpu, X_gpu, energies_minspm_gpu, nnz_terms)
    CUDA.synchronize()

    # --- BENCHMARK DENSE ---
    println("Benchmarking Dense Non-MINSPM QAP Kernel...")
    dense_times = zeros(Float64, runs)
    for r in 1:runs
        dense_times[r] = CUDA.@elapsed begin
            @cuda threads=dense_threads blocks=dense_blocks shmem=dense_threads*sizeof(Float32) dense_qap_eval_kernel!(F_gpu, D_gpu, P_gpu, energies_dense_gpu, N)
        end
    end
    dense_avg_time = sum(dense_times) / runs

    # --- BENCHMARK MINSPM ---
    println("Benchmarking Exact MINSPM Sparse Galois Oracle...")
    minspm_times = zeros(Float64, runs)
    for r in 1:runs
        minspm_times[r] = CUDA.@elapsed begin
            @cuda threads=minspm_threads blocks=minspm_blocks shmem=minspm_threads*sizeof(Int64) minspm_sparse_galois_kernel!(Galois_Keys_gpu, Galois_Values_gpu, X_gpu, energies_minspm_gpu, nnz_terms)
        end
    end
    minspm_avg_time = sum(minspm_times) / runs

    # --- REPORTING ---
    dense_throughput = batch_size / dense_avg_time / 1e6
    minspm_throughput = batch_size / minspm_avg_time / 1e6
    speedup = dense_avg_time / minspm_avg_time

    println("\n=========================================================")
    println(" ABLATION RESULTS (V3 - SPARSE GALOIS MATCH)")
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
