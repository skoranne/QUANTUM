# mps_encoding_benchmark.jl
# Requires: Pkg.add(["CUDA", "cuTENSOR", "SparseArrays", "BenchmarkTools", "LinearAlgebra"])

using CUDA
using cuTENSOR
using SparseArrays
using BenchmarkTools
using LinearAlgebra
using Random

# ==============================================================================
# 1. MINSPM Dictionary-Encoded Sparse Engine for MPS Tensors
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
# 2. MPS Sequential Encoding Simulation (Ran's Protocol)
# ==============================================================================
function generate_mps_encoded_tensor(num_sites::Int, bond_dim::Int)
    # Simulating the tensor block generated from sequential MPS encoding
    dim = min(4096, bond_dim * bond_dim)
    println("Generating MPS Sequential Encoding Tensor Block (Dim: $dim x $dim)...")

    unique_amplitudes = Float64[1.0, -1.0, 0.5, -0.5, 1/sqrt(2), -1/sqrt(2)]
    dense_block = zeros(Float64, dim, dim)

    Random.seed!(1337)
    fill_rate = 0.006 # High sparsity typical of physical MPS bonds
    target_nnz = round(Int, (dim * dim) * fill_rate)

    for _ in 1:target_nnz
        r, c = rand(1:dim), rand(1:dim)
        dense_block[r, c] = rand(unique_amplitudes)
    end

    sparse_host = sparse(dense_block)
    return sparse_host, unique_amplitudes
end

# ==============================================================================
# 3. Benchmark Execution Suite
# ==============================================================================
function run_mps_encoding_benchmark()
    num_sites = 70
    chi = 256

    sparse_host, codebook = generate_mps_encoded_tensor(num_sites, chi)
    nnz_count = nnz(sparse_host)
    dim = size(sparse_host, 2)

    println("Non-Zero Elements (NNZ): $nnz_count out of $(dim*dim)")
    println("Unique Codebook Size   : $(length(codebook)) values")

    # Setup GPU Buffers for MINSPM
    d_x = CUDA.rand(Float64, dim)
    val_to_idx = Dict(v => UInt8(i) for (i, v) in enumerate(codebook))
    nz_idx_host = UInt8[get(val_to_idx, v, UInt8(1)) for v in sparse_host.nzval]

    A_minspm = CuMINSPMTensor(
        size(sparse_host),
        CuArray(sparse_host.colptr),
        CuArray(sparse_host.rowval),
        CuArray(nz_idx_host),
        CuArray(codebook)
    )
    d_y_minspm = CUDA.zeros(Float64, dim)

    # Setup GPU Buffers for cuTENSOR Dense Engine
    d_dense = CuArray(Matrix(sparse_host))
    ct_A = CuTensor(d_dense, ['i', 'j'])
    ct_x = CuTensor(d_x, ['j'])

    sweeps = 70
    println("\n--- Benchmarking MPS Sequential Encoding Sweeps ($sweeps Iterations) ---")

    println("Benchmarking cuTENSOR (Dense Contraction)...")
    time_dense = @elapsed begin
        for _ in 1:sweeps
            res = ct_A * ct_x
            copyto!(d_y_minspm, res.data)
        end
    end
    println("cuTENSOR Time: $(round(time_dense, digits=4)) seconds")

    println("Benchmarking MINSPM (Dictionary Sparse Contraction)...")
    time_sparse = @elapsed begin
        for _ in 1:sweeps
            fill!(d_y_minspm, 0.0)
            contract_minspm!(d_y_minspm, A_minspm, d_x)
        end
        CUDA.synchronize()
    end
    println("MINSPM Time  : $(round(time_sparse, digits=4)) seconds")

    speedup = time_dense / time_sparse
    println("\nMeasured Speedup Factor: $(round(speedup, digits=2))x")
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_mps_encoding_benchmark()
end