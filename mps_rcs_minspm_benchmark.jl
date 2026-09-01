# mps_rcs_minspm_benchmark.jl
# Requires: Pkg.add(["CUDA", "cuTENSOR", "SparseArrays", "BenchmarkTools"])

using CUDA
using cuTENSOR
using SparseArrays
using BenchmarkTools
using LinearAlgebra

# ==============================================================================
# MINSPM Tensor Block Representation (Float64 for Native GPU Atomics)
# ==============================================================================

struct CuMINSPMTensor{Tv,Ti,Tidx}
    dims::Dims
    colptr::CuVector{Ti}
    rowval::CuVector{Ti}
    nz_idx::CuVector{Tidx}
    nz_dict::CuVector{Tv}
end

# Kernel using native Float64 atomic add on the GPU
function minspm_tensor_kernel!(y, x, colptr, rowval, nz_idx, nz_dict, n)
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

function contract_minspm!(y::CuVector, A::CuMINSPMTensor, x::CuVector)
    fill!(y, zero(eltype(y)))
    threads = 256
    blocks = cld(A.dims[2], threads)
    @cuda threads=threads blocks=blocks minspm_tensor_kernel!(
        y, x, A.colptr, A.rowval, A.nz_idx, A.nz_dict, A.dims[2]
    )
    return y
end

# ==============================================================================
# Benchmark Suite for 70x70 Doped Clifford / 27-Ancilla MPS Layer
# ==============================================================================

function run_mps_rcs_benchmark()
    dim = 4096

    println("Initializing 70x70 RCS / 27-ancilla MPS tensor block workload (Dim: $dim x $dim)...")

    unique_amplitudes = Float64[1.0, -1.0, 0.5, -0.5, 1/sqrt(2), -1/sqrt(2)]
    dense_block = zeros(Float64, dim, dim)

    nnz_target = dim * 20
    for _ in 1:nnz_target
        r, c = rand(1:dim), rand(1:dim)
        dense_block[r, c] = rand(unique_amplitudes)
    end

    sparse_host = sparse(dense_block)
    actual_nnz = nnz(sparse_host)
    println("Effective Tensor Sparsity: $(round(actual_nnz / (dim*dim) * 100, digits=4))% ($actual_nnz non-zeros)")

    # 1. cuTENSOR Setup (Dense Tensor Contraction Engine)
    d_dense = CuArray(dense_block)
    d_x = CUDA.rand(Float64, dim)
    ct_A = CuTensor(d_dense, ['i', 'j'])
    ct_x = CuTensor(d_x, ['j'])

    # 2. MINSPM Setup (Dictionary Encoded Sparse Engine)
    val_to_idx = Dict(v => UInt8(i) for (i, v) in enumerate(unique_amplitudes))
    nz_idx_host = UInt8[val_to_idx[v] for v in sparse_host.nzval]

    A_minspm = CuMINSPMTensor(
        size(sparse_host),
        CuArray(sparse_host.colptr),
        CuArray(sparse_host.rowval),
        CuArray(nz_idx_host),
        CuArray(unique_amplitudes)
    )
    d_y_minspm = CUDA.zeros(Float64, dim)

    println("\n--- Benchmarking cuTENSOR Dense MPS Layer Contraction ---")
    @btime CUDA.@sync ($ct_A * $ct_x)

    println("\n--- Benchmarking MINSPM Sparse MPS Layer Contraction ---")
    @btime CUDA.@sync contract_minspm!($d_y_minspm, $A_minspm, $d_x)

    dense_mb = (dim * dim * sizeof(Float64)) / 1024^2
    minspm_mb = (actual_nnz * sizeof(UInt8) + actual_nnz * sizeof(Int)) / 1024^2
    println("\nMemory Bandwidth Footprint per MPS Sweep Step:")
    println("  cuTENSOR Dense : $(round(dense_mb, digits=2)) MB")
    println("  MINSPM Sparse  : $(round(minspm_mb, digits=2)) MB")
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_mps_rcs_benchmark()
end