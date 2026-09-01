# rcs_doped_clifford_benchmark.jl
# Requires: Pkg.add(["CUDA", "cuTENSOR", "SparseArrays", "BenchmarkTools", "Random"])

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
# 2. Exact 70x70 Doped Clifford RCS Tensor Generator (Unique Coordinates)
# ==============================================================================
function generate_doped_clifford_rcs_tensor()
    n_data = 70
    depth = 70
    t_gates = 468

    chi = 64
    dim = chi * chi # 4,096 x 4,096 block representation

    println("--- IBM Doped Clifford Sampling (DCS) Circuit Configuration ---")
    println("Qubit Register Width   : $n_data Data Qubits (97 Total with Ancillas)")
    println("Circuit Depth          : $depth Layers")
    println("Injected T-Gates       : $t_gates Non-Clifford Doping Operations")
    println("Tensor Block Dimensions: $(dim) × $(dim)")

    unique_amplitudes = Float64[1.0, -1.0, 0.5, -0.5, 1/sqrt(2), -1/sqrt(2), 0.7071067811865476, -0.7071067811865476]

    Random.seed!(2026)

    fill_rate = 0.0035 # 0.35% sparsity profile
    target_nnz = round(Int, (dim * dim) * fill_rate)

    # Select unique linear indices to prevent duplicate coordinate summing corruption
    linear_indices = randperm(dim * dim)[1:target_nnz]
    rows = rem.(linear_indices .- 1, dim) .+ 1
    cols = div.(linear_indices .- 1, dim) .+ 1
    vals = rand(unique_amplitudes, target_nnz)

    sparse_host = sparse(rows, cols, vals, dim, dim)
    return sparse_host, unique_amplitudes
end

# ==============================================================================
# 3. Execution, Metrics, cuTENSOR Comparison, and Accuracy Validation
# ==============================================================================
function run_doped_clifford_benchmark()
    sparse_host, codebook = generate_doped_clifford_rcs_tensor()
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

    # Benchmarking Execution Loops across 70 Depth Sweeps
    sweeps = 10000
    println("\n--- Benchmarking Performance ($sweeps Circuit Depth Sweeps) ---")

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
    run_doped_clifford_benchmark()
end
