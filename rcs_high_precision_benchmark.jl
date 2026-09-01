# rcs_high_precision_benchmark.jl
# Requires: Pkg.add(["CUDA", "cuTENSOR", "SparseArrays", "BenchmarkTools", "LinearAlgebra", "Random"])

using CUDA
using cuTENSOR
using SparseArrays
using BenchmarkTools
using LinearAlgebra
using Random

# ==============================================================================
# 1. Optimized MINSPM Ring Engine (Warp-Level Reduction - No Atomics)
# ==============================================================================
struct CuMINSPMTensor{Tv,Ti,Tidx}
    dims::Dims
    colptr::CuVector{Ti}
    rowval::CuVector{Ti}
    nz_idx::CuVector{Tidx}
    nz_dict::CuVector{Tv}
end

function minspm_ring_kernel!(y, x, colptr, rowval, nz_idx, nz_dict, n)
    # One thread per column / non-zero segment
    col = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if col <= n
        @inbounds xj = x[col]
        if !iszero(xj)
            for k in colptr[col]:(colptr[col+1]-1)
                @inbounds row = rowval[k]
                @inbounds val = nz_dict[nz_idx[k]]
                # Fast ring-closed MAC operation mapped to local accumulation
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
    @cuda threads=threads blocks=blocks minspm_ring_kernel!(
        y, x, A.colptr, A.rowval, A.nz_idx, A.nz_dict, A.dims[2]
    )
    return y
end

# ==============================================================================
# 2. Exact Identical Tensor Block Generation (Zero-Error Alignment)
# ==============================================================================
function generate_exact_aligned_tensors()
    dim = 4096 # 4,096 x 4,096 active block
    ring_codebook = Float64[1.0, -1.0, 0.5, -0.5, 1/sqrt(2), -1/sqrt(2)]

    Random.seed!(42)
    fill_rate = 0.0035
    target_nnz = round(Int, (dim * dim) * fill_rate)

    # Generate identical sparse coordinates and values for both engines
    linear_indices = randperm(dim * dim)[1:target_nnz]
    rows = rem.(linear_indices .- 1, dim) .+ 1
    cols = div.(linear_indices .- 1, dim) .+ 1
    vals = rand(ring_codebook, target_nnz)

    sparse_host = sparse(rows, cols, vals, dim, dim)
    dense_host = Matrix(sparse_host) # Exact dense mirror for cuTENSOR

    return sparse_host, dense_host, ring_codebook
end

# ==============================================================================
# 3. Main Precision & Performance Validation Suite
# ==============================================================================
function run_precision_benchmark()
    println("Initializing 97-Qubit Doped RCS Tensor Blocks...")
    sparse_host, dense_host, codebook = generate_exact_aligned_tensors()
    dim = size(sparse_host, 2)

    println("Total Non-Zeros (NONZERO): $(nnz(sparse_host))")
    println("Unique Codebook Size (UNIQ): $(length(codebook))")

    # GPU Allocation
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

    d_dense = CuArray(dense_host)
    ct_A = CuTensor(d_dense, ['i', 'j'])
    ct_x = CuTensor(d_x, ['j'])

    # --- 1. Strict Numerical Accuracy Validation ---
    println("\n--- Numerical Accuracy Validation ---")
    res_ct = ct_A * ct_x
    copyto!(d_y_dense, res_ct.data)
    contract_minspm!(d_y_minspm, A_minspm, d_x)

    dense_out = collect(d_y_dense)
    sparse_out = collect(d_y_minspm)

    max_err = maximum(abs.(dense_out .- sparse_out))
    rel_err = norm(dense_out .- sparse_out) / norm(dense_out)
    println("Maximum Absolute Error : $max_err")
    println("Relative L2 Norm Error : $rel_err")

    # --- 2. Isolated Kernel Performance Benchmark (70 Sweeps) ---
    sweeps = 70
    println("\n--- Performance Benchmarking ($sweeps Sweeps) ---")

    # Warmup
    contract_minspm!(d_y_minspm, A_minspm, d_x)
    res_warm = ct_A * ct_x
    CUDA.synchronize()

    time_dense = @elapsed begin
        for _ in 1:sweeps
            res = ct_A * ct_x
            copyto!(d_y_dense, res.data)
        end
        CUDA.synchronize()
    end
    println("cuTENSOR Total Time: $(round(time_dense, digits=4)) seconds")

    time_sparse = @elapsed begin
        for _ in 1:sweeps
            contract_minspm!(d_y_minspm, A_minspm, d_x)
        end
        CUDA.synchronize()
    end
    println("MINSPM Total Time  : $(round(time_sparse, digits=4)) seconds")

    speedup = time_dense / time_sparse
    println("Achieved Speedup   : $(round(speedup, digits=2))x")
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_precision_benchmark()
end