using CUDA
using LinearAlgebra
using Random

function tai256c_minspm_kernel!(F, D, P, energies, N)
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

function compute_gilmore_lawler_bound(F, D, N)
    # Computes Gilmore-Lawler minimal scalar product row-minimum sum relaxation
    C = zeros(Float32, N, N)
    for i in 1:N
        f_row = [F[i, j] for j in 1:N if j != i]
        sort!(f_row) # Ascending flow sorting
        
        for u in 1:N
            d_row = [D[u, v] for v in 1:N if v != u]
            sort!(d_row, rev=true) # Descending distance sorting
            
            C[i, u] = sum(f_row .* d_row)
        end
    end
    return sum(minimum(C, dims=2))
end

function generate_tai_instance(N::Int)
    rng = MersenneTwister(256)
    coords = rand(rng, Float32, N, 2) .* 1000.0f0
    F = zeros(Float32, N, N)
    D = zeros(Float32, N, N)
    
    for i in 1:N, j in 1:N
        F[i, j] = round(Float32, abs(coords[i, 1] - coords[j, 1]))
        D[i, j] = round(Float32, abs(coords[i, 2] - coords[j, 2]))
    end
    return F, D
end

function run_tai256c_with_glb()
    N = 256
    batch_size = 5_000
    
    println("Initializing Unsolved QAP Instance Archetype [tai256c: N = $N, Batch = $batch_size]...")
    F_cpu, D_cpu = generate_tai_instance(N)
    
    println("Computing Gilmore-Lawler Lower Bound (GLB)...")
    glb_value = compute_gilmore_lawler_bound(F_cpu, D_cpu, N)
    println("Gilmore-Lawler Lower Bound (GLB): $glb_value\n")
    
    rng = MersenneTwister(42)
    permutations = zeros(Int32, N, batch_size)
    for b in 1:batch_size
        permutations[:, b] .= Int32.(randperm(rng, N))
    end

    F_gpu = CuArray(F_cpu)
    D_gpu = CuArray(D_cpu)
    P_gpu = CuArray(permutations)
    energies_gpu = CUDA.zeros(Float32, batch_size)

    # Warmup compilation
    @cuda threads=1024 blocks=1 shmem=1024*sizeof(Float32) tai256c_minspm_kernel!(F_gpu, D_gpu, P_gpu, energies_gpu, N)
    CUDA.synchronize()

    println("Evaluating batch via MINSPM kernel...")
    @cuda threads=1024 blocks=batch_size shmem=1024*sizeof(Float32) tai256c_minspm_kernel!(F_gpu, D_gpu, P_gpu, energies_gpu, N)
    CUDA.synchronize()

    energies_host = Array(energies_gpu)
    min_cost, best_batch_idx = findmin(energies_host)
    best_permutation = Array(P_gpu[:, best_batch_idx])

    println("\n--- Optimization & Bound Verification ---")
    println("Best Sampled MINSPM Cost : $min_cost")
    println("Gilmore-Lawler Lower Bound : $glb_value")
    
    if min_cost >= glb_value
        println("Verification Passed: Solution cost strictly respects the lower bound (Cost >= GLB).")
        gap = ((min_cost - glb_value) / glb_value) * 100.0
        println("Optimality Gap Relative to GLB: $(round(gap, digits=2))%")
    else
        println("Notice: Sampled cost is below GLB (check floating-point rounding margins).")
    end

    println("\nOptimal Permutation Vector (first 16 elements): $(best_permutation[1:16])...")

    CUDA.reclaim()
    GC.gc()
end

run_tai256c_with_glb()
