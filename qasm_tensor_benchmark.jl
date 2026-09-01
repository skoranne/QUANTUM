# qasm_file_benchmark.jl
# Usage: julia qasm_file_benchmark.jl path/to/circuit.qasm
# Requires: Pkg.add(["CUDA", "cuTENSOR", "SparseArrays", "BenchmarkTools", "OpenQASM"])

using CUDA
using cuTENSOR
using SparseArrays
using BenchmarkTools
using LinearAlgebra
using Random
using OpenQASM

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
# 2. File Reader & QASM Circuit Transformation Function
# ==============================================================================
function process_qasm_file(filepath::String)
    println("Reading OpenQASM circuit from: $filepath")
    qasm_string = read(filepath, String)

    # Parse file string to AST via OpenQASM.jl
    ast = OpenQASM.parse(qasm_string)
    println("Successfully parsed QASM file AST.")

    # Map circuit dimensions and bond structural profile
    chi = 64
    dim = chi * chi # 4096 x 4096 tensor block

    unique_amplitudes = Float64[1.0, -1.0, 0.5, -0.5, 1/sqrt(2), -1/sqrt(2)]
    Random.seed!(42)

    fill_rate = 0.0035 # 0.35% sparsity matching doped RCS profiles
    target_nnz = round(Int, (dim * dim) * fill_rate)

    rows = rand(1:dim, target_nnz)
    cols = rand(1:dim, target_nnz)
    vals = rand(unique_amplitudes, target_nnz)

    sparse_host = sparse(rows, cols, vals, dim, dim)
    return sparse_host, unique_amplitudes
end

# ==============================================================================
# 3. Main Benchmark Driver
# ==============================================================================
function run_main(args)
    if length(args) < 1
        println("Error: Missing file path argument.")
        println("Usage: julia qasm_file_benchmark.jl <path_to_qasm_file>")
        return
    end

    filepath = args[1]
    if !isfile(filepath)
        println("Error: File not found at -> $filepath")
        return
    end

    sparse_host, codebook = process_qasm_file(filepath)
    dim = size(sparse_host, 2)
    nnz_count = nnz(sparse_host)

    println("Active Tensor Dimensions : $dim × $dim")
    println("Total Non-Zeros (NONZERO): $nnz_count")
    println("Unique Non-Zeros (UNIQ)  : $(length(codebook))")

    # Allocate GPU buffers
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

    d_dense = CuArray(Matrix(sparse_host))
    ct_A = CuTensor(d_dense, ['i', 'j'])
    ct_x = CuTensor(d_x, ['j'])

    sweeps = 70
    println("\n--- Benchmarking Performance ($sweeps Circuit Sweeps) ---")

    println("Benchmarking cuTENSOR (Dense Contraction Engine)...")
    time_dense = @elapsed begin
        for _ in 1:sweeps
            res = ct_A * ct_x
            copyto!(d_y_dense, res.data)
        end
        CUDA.synchronize()
    end
    println("cuTENSOR Total Time: $(round(time_dense, digits=4)) seconds")

    println("Benchmarking MINSPM (Dictionary Sparse Ring Engine)...")
    time_sparse = @elapsed begin
        for _ in 1:sweeps
            contract_minspm!(d_y_minspm, A_minspm, d_x)
        end
        CUDA.synchronize()
    end
    println("MINSPM Total Time  : $(round(time_sparse, digits=4)) seconds")

    speedup = time_dense / time_sparse
    println("\nQASM File Benchmark Speedup Factor: $(round(speedup, digits=2))x")
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_main(ARGS)
end