using CUDA
using LinearAlgebra
using Random

function minspm_qap_eval_kernel!(F, D, P, energies, N)
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

function load_qaplib_tai256c(filepath::String)
    open(filepath, "r") do io
        tokens = split(read(io, String))
        N = parse(Int, tokens[1])
        idx = 2
        
        F = zeros(Float32, N, N)
        for i in 1:N, j in 1:N
            F[i, j] = parse(Float32, tokens[idx]); idx += 1
        end
        
        D = zeros(Float32, N, N)
        for i in 1:N, j in 1:N
            D[i, j] = parse(Float32, tokens[idx]); idx += 1
        end
        return N, F, D
    end
end

function run_advanced_tabu_minspm()
    filepath = "tai256c.dat"
    if !isfile(filepath)
        error("Please place 'tai256c.dat' in the working directory.")
    end

    N, F_cpu, D_cpu = load_qaplib_tai256c(filepath)
    batch_size = 4_096
    max_iterations = 2_000
    
    known_ub = 44_759_296.0f0
    known_lb = 44_200_812.0f0

    println("Loaded QAPLIB tai256c [N = $N] for Advanced Tabu-MINSPM Optimization")
    
    F_gpu = CuArray(F_cpu)
    D_gpu = CuArray(D_cpu)
    
    rng = MersenneTwister(1337)
    
    # Initialize diverse population using spectral/row-sum heuristics and random permutations
    row_sums_f = sum(F_cpu, dims=2)[:]
    base_order = sortperm(row_sums_f)
    
    P_host = zeros(Int32, N, batch_size)
    for b in 1:batch_size
        if b <= 64
            # Perturbed variants of the greedy baseline
            p = copy(base_order)
            for _ in 1:12
                i1, i2 = rand(rng, 1:N), rand(rng, 1:N)
                p[i1], p[i2] = p[i2], p[i1]
            end
            P_host[:, b] .= Int32.(p)
        else
            P_host[:, b] .= Int32.(randperm(rng, N))
        end
    end
    
    P_gpu = CuArray(P_host)
    energies_gpu = CUDA.zeros(Float32, batch_size)

    @cuda threads=1024 blocks=batch_size shmem=1024*sizeof(Float32) minspm_qap_eval_kernel!(F_gpu, D_gpu, P_gpu, energies_gpu, N)
    CUDA.synchronize()
    
    E_host = Array(energies_gpu)
    global_best_cost = minimum(E_host)
    best_idx = argmin(E_host)
    global_best_perm = copy(P_host[:, best_idx])

    # Tabu tenure tracking matrix per chain to prevent cyclic 2-opt traps
    tabu_tenure = zeros(Int, N, N, batch_size)

    println("Executing Parallel GPU Tabu Search with MINSPM Evaluation Engine...")
    for iter in 1:max_iterations
        P_host_new = copy(P_host)
        
        # Guided 2-opt / Multi-swap neighborhood generation avoiding tabu moves
        for b in 1:batch_size
            for _ in 1:4
                i1 = rand(rng, 1:N)
                i2 = rand(rng, 1:N)
                if i1 != i2 && tabu_tenure[i1, i2, b] <= iter
                    P_host_new[i1, b], P_host_new[i2, b] = P_host_new[i2, b], P_host_new[i1, b]
                    tabu_tenure[i1, i2, b] = iter + 7 # Tabu tenure duration
                    tabu_tenure[i2, i1, b] = iter + 7
                end
            end
        end
        
        copyto!(P_gpu, P_host_new)
        energies_gpu_new = CUDA.zeros(Float32, batch_size)
        
        @cuda threads=1024 blocks=batch_size shmem=1024*sizeof(Float32) minspm_qap_eval_kernel!(F_gpu, D_gpu, P_gpu, energies_gpu_new, N)
        CUDA.synchronize()
        
        E_host_new = Array(energies_gpu_new)
        
        for b in 1:batch_size
            # Aspiration criterion or strict improvement acceptance
            if E_host_new[b] < E_host[b]
                P_host[:, b] .= P_host_new[:, b]
                E_host[b] = E_host_new[b]
                
                if E_host[b] < global_best_cost
                    global_best_cost = E_host[b]
                    global_best_perm .= P_host[:, b]
                end
            elseif rand(rng) < 0.03 # Diversification escape probability
                P_host[:, b] .= Int32.(randperm(rng, N))
            end
        end
        
        if iter % 200 == 0 || iter == 1
            println("Iteration $iter | Global Best Cost Found: $global_best_cost | Target UB: $known_ub")
        end
    end

    println("\n--- Advanced MINSPM-Tabu Results ---")
    println("Global Best Cost Achieved  : $global_best_cost")
    println("QAPLIB Best Known UB       : $known_ub")
    
    gap = ((global_best_cost - known_ub) / known_ub) * 100.0
    println("Optimality Gap vs UB       : $(round(gap, digits=3))%")

    CUDA.reclaim()
    GC.gc()
end

run_advanced_tabu_minspm()
