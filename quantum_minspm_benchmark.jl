# File   : quantum_minspm_benchmark.jl
# Author : Sandeep Koranne (C) Siemens EDA 2026
# Requires: Pkg.add(["CUDA", "cuTENSOR", "SparseArrays", "BenchmarkTools"])

using CUDA
using SparseArrays
using BenchmarkTools
using LinearAlgebra

# ==============================================================================
# 1. MINSPM (Unique Non-Zero) Custom GPU Kernel & Data Structures
# ==============================================================================

struct CuMINSPM{Tv,Ti,Tidx}
    m::Int
    n::Int
    colptr::CuVector{Ti}
    rowval::CuVector{Ti}
    nz_idx::CuVector{Tidx}   # Low-bandwidth indices (e.g., UInt8)
    nz_dict::CuVector{Tv}    # High-precision codebook fits in L1 cache
end

function minspm_spmv_kernel!(y, x, colptr, rowval, nz_idx, nz_dict, n)
    col = (blockIdx().x - 1) * blockDim().x + threadIdx().x

    if col <= n
        @inbounds xj = x[col]
        if !iszero(xj)
            for k in colptr[col]:(colptr[col+1]-1)
                @inbounds row = rowval[k]
                @inbounds val = nz_dict[nz_idx[k]]
                CUDA.@atomic y[row] += val * xj
            end
        end
    end
    return nothing
end

function mul_minspm_gpu!(y::CuVector, A::CuMINSPM, x::CuVector)
    fill!(y, zero(eltype(y)))
    threads = 256
    blocks = cld(A.n, threads)

    @cuda threads=threads blocks=blocks minspm_spmv_kernel!(
        y, x, A.colptr, A.rowval, A.nz_idx, A.nz_dict, A.n
    )
    return y
end

function CuMINSPM(A::SparseMatrixCSC{Tv,Ti}; index_type=UInt8) where {Tv,Ti}
    unique_vals = unique(A.nzval)
    if length(unique_vals) > typemax(index_type)
        error("Unique values exceed capacity of $index_type codebook.")
    end

    val_to_idx = Dict(v => index_type(i) for (i, v) in enumerate(unique_vals))
    nz_idx_host = index_type[val_to_idx[v] for v in A.nzval]

    return CuMINSPM{Tv,Ti,index_type}(
        A.m, A.n,
        CuArray(A.colptr),
        CuArray(A.rowval),
        CuArray(nz_idx_host),
        CuArray(unique_vals)
    )
end

# ==============================================================================
# 2. Benchmark Execution Suite
# ==============================================================================

function run_benchmark()
    N = 14
    dim = 2^N  # 16,384 x 16,384 matrix (~268 Million elements)

    println("Generating Quantum Operator Tensor (Dim: $dim x $dim)...")

    unique_amplitudes = Float64[1.0, -1.0, 0.5, -0.5, 1/sqrt(2)]
    dense_matrix = zeros(Float64, dim, dim)

    nnz_target = dim * 12
    for _ in 1:nnz_target
        r, c = rand(1:dim), rand(1:dim)
        dense_matrix[r, c] = rand(unique_amplitudes)
    end

    sparse_host = sparse(dense_matrix)
    actual_nnz = nnz(sparse_host)
    println("Sparsity: $(round(actual_nnz / (dim*dim) * 100, digits=4))% ($actual_nnz non-zeros)")

    println("\nLoading Dense Tensor to GPU (for cuBLAS / Dense GEMV)...")
    d_dense_matrix = CuArray(dense_matrix)
    d_x = CUDA.rand(Float64, dim)
    d_y_dense = CUDA.zeros(Float64, dim)

    println("Compressing to MINSPM and loading to GPU...")
    A_gpu_minspm = CuMINSPM(sparse_host, index_type=UInt8)
    d_y_minspm = CUDA.zeros(Float64, dim)

    # Correctness check using native CuArray multiplication vs MINSPM kernel
    mul!(d_y_dense, d_dense_matrix, d_x)
    mul_minspm_gpu!(d_y_minspm, A_gpu_minspm, d_x)

    res_dense = collect(d_y_dense)
    res_minspm = collect(d_y_minspm)

    if isapprox(res_dense, res_minspm, rtol=1e-5)
        println("Correctness verified: Dense GPU and MINSPM outputs match.")
    else
        @warn "Output mismatch!"
    end

    # ==========================================================================
    # 3. Profiling & Benchmarking
    # ==========================================================================
    println("\n--- Benchmarking Dense GPU Matrix-Vector Multiplication ---")
    @btime CUDA.@sync mul!($d_y_dense, $d_dense_matrix, $d_x)

    println("\n--- Benchmarking MINSPM (Dictionary ALUs via CUDA.jl) ---")
    @btime CUDA.@sync mul_minspm_gpu!($d_y_minspm, $A_gpu_minspm, $d_x)

    dense_mb = (dim * dim * sizeof(Float64)) / 1024^2
    minspm_mb = (actual_nnz * sizeof(UInt8) + actual_nnz * sizeof(Int)) / 1024^2

    println("\nMemory Bandwidth Footprint:")
    println("  Dense Matrix : $(round(dense_mb, digits=2)) MB")
    println("  MINSPM Sparse: $(round(minspm_mb, digits=2)) MB")
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_benchmark()
end
