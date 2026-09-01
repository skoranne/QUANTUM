# rcs_large_bond_benchmark.jl
# Requires: Pkg.add(["CUDA", "cuTENSOR", "SparseArrays", "BenchmarkTools"])

using CUDA
using cuTENSOR
using SparseArrays
using BenchmarkTools
using LinearAlgebra
using Random

# ==============================================================================
# 1. MINSPM Dictionary-Encoded Sparse Engine
# ==============================================================================
struct CuMINSPMTensor{Tv,Ti,Tidx}
    dims::Dims
    colptr::CuVector{Ti}
    rowval::CuVector{Ti}
    nz_idx::CuVector{Tidx}
    nz_dict::CuVector{Tv}
end

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
# 2. Large Bond Dimension RCS Tensor Block Generation (Chi = 64 -> Dim = 4,096)
# Note: Scaled to 4096 to safely fit the dense CuTensor comparison inside standard VRAM limits.
# ==============================================================================
function generate_large_scale_rcs_tensor()
    chi = 64
    dim = chi * chi

    println("--- Large-Scale RCS Configuration ---")
    println("Target MPS Bond Dimension (chi): $chi")
    println("Active Tensor Block Dimensions : $(dim) × $(dim)")

    unique_amplitudes = Float64[1.0, -1.0, 0.5, -0.5, 1/sqrt(2), -1/sqrt(2)]

    Random.seed!(42)
    fill_rate = 0.002 # 0.2% sparsity
    target_nnz = round(Int, (dim * dim) * fill_rate)

    rows = rand(1:dim, target_nnz)
    cols = rand(1:dim, target_nnz)
    vals = rand(unique_amplitudes, target_nnz)

    sparse_host = sparse(rows, cols, vals, dim, dim)
    return sparse_host, unique_amplitudes
end

# ==============================================================================
# 3. Execution, Metrics, cuTENSOR Comparison, and Accuracy Check
# ==============================================================================
function run_large_bond_comparison()
    sparse_host, codebook = generate_large_scale_rcs_tensor()
    dim = size(sparse_host, 2)

    nnz_count = nnz(sparse_host)
    unique_vals_found = unique(sparse_host.nzval)
    unique_count = length(unique_vals_found)

    println("\n--- Tensor Structure Metrics ---")
    println("Total Non-Zeros (NONZERO)        : $nnz_count")
    println("Number of Unique Non-Zeros (UNIQ): $unique_count")
    println("Unique Codebook Values Found     : $unique_vals_found")

    # Persistent GPU allocation
    d_x = CUDA.rand(Float64, dim)
    d_y_minspm = CUDA.zeros(Float64, dim)
    d_y_dense = CUDA.zeros(Float64, dim)

    val_to_idx = Dict(v => UInt8(i) for (i, v) in enumerate(codebook))
    nz_idx_host = UInt8[get(val_to_idx, v, UInt8(1)) for v in sparse_host.nzval]

    A_minspm = CuMINSPMTensor(
        size(sparse_host),
        CuArray(sparse_host.colptr),
        CuArray(sparse_host.rowval),
        CuArray(nz_idx_host),
        CuArray(codebook)
    )

    # Setup cuTENSOR dense structures
    d_dense = CuArray(Matrix(sparse_host))
    ct_A = CuTensor(d_dense, ['i', 'j'])
    ct_x = CuTensor(d_x, ['j'])

    # Accuracy Validation Check
    println("\n--- Numerical Accuracy Validation ---")
    res_dense_ct = ct_A * ct_x
    copyto!(d_y_dense, res_dense_ct.data)
    contract_minspm!(d_y_minspm, A_minspm, d_x)

    dense_host_out = collect(d_y_dense)
    minspm_host_out = collect(d_y_minspm)

    max_err = maximum(abs.(dense_host_out .- minspm_host_out))
    rel_err = norm(dense_host_out .- minspm_host_out) / norm(dense_host_out)
    println("Maximum Absolute Error : $max_err")
    println("Relative L2 Norm Error : $rel_err")

    # Benchmarking Execution Loops
    sweeps = 50
    println("\n--- Benchmarking Performance ($sweeps Sweeps) ---")

    println("Benchmarking cuTENSOR (Dense Contraction Engine)...")
    time_dense = @elapsed begin
        for _ in 1:sweeps
            res = ct_A * ct_x
            copyto!(d_y_dense, res.data)
        end
        CUDA.synchronize()
    end
    println("cuTENSOR Total Time: $(round(time_dense, digits=4)) seconds")

    println("Benchmarking MINSPM (Dictionary Sparse Engine)...")
    time_sparse = @elapsed begin
        for _ in 1:sweeps
            contract_minspm!(d_y_minspm, A_minspm, d_x)
        end
        CUDA.synchronize()
    end
    println("MINSPM Total Time  : $(round(time_sparse, digits=4)) seconds")

    speedup = time_dense / time_sparse
    println("\nPerformance Speedup Factor: $(round(speedup, digits=2))x")
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_large_bond_comparison()
end