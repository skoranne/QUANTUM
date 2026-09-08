using CUDA
using LinearAlgebra
using Random

function minspm_eval_kernel!(F, D, P, energies, N)
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

function load_qaplib_instance(filepath::String)
    open(filepath, "r") do io
        tokens = split(read(io, String))
        if isempty(tokens)
            error("QAPLIB file is empty or missing.")
        end
        N = parse(Int, tokens[1])
        idx = 2
        
        # Parse Flow Matrix F
        F = zeros(Float32, N, N)
        for i in 1:N, j in 1:N
            F[i, j] = parse(Float32, tokens[idx])
            idx += 1
        end
        
        # Parse Distance Matrix D
        D = zeros(Float32, N, N)
        for i in 1:N, j in 1:N
            D[i, j] = parse(Float32, tokens[idx])
            idx += 1
        end
        return N, F, D
    end
end

function run_true_tai256c_optimization()
    filepath = "tai256c.txt"
    if !isfile(filepath)
        error("Please place the official 'tai256c.dat' file from QAPLIB into your working directory.")
    end

    N, F_cpu, D_cpu = load_qaplib_instance(filepath)
    batch_size = 4_096
    num_iterations = 2_500
    
    known_ub = 44_759_296.0f0
    known_lb = 44_200_812.0f0

    println("Loaded True QAPLIB Instance: tai256c [N = $N]")
    
    F_gpu = CuArray(F_cpu)
    D_gpu = CuArray(D_cpu)
    
    rng = MersenneTwister(42)
    
    row_sums_f = sum(F_cpu, dims=2)[:]
    row_sums_d = sum(D_cpu, dims=2)[:]
    base_order = sortperm(row_sums_f)
    
    P_host = zeros(Int32, N, batch_size)
    for b in 1:batch_size
        if b == 1
            P_host[:, b] .= Int32.(base_order)
        else
            P_host[:, b] .= Int32.(randperm(rng, N))
        end
    end
    
    P_gpu = CuArray(P_host)
    energies_gpu = CUDA.zeros(Float32, batch_size)

    @cuda threads=1024 blocks=batch_size shmem=1024*sizeof(Float32) minspm_eval_kernel!(F_gpu, D_gpu, P_gpu, energies_gpu, N)
    CUDA.synchronize()
    
    E_host = Array(energies_gpu)
    best_ever_cost = minimum(E_host)
    
    T = 500_000.0f0
    cooling_rate = 0.996f0

    println("Running High-Throughput GPU Search on True QAPLIB Matrices...")
    for iter in 1:num_iterations
        P_host_new = copy(P_host)
        for b in 1:batch_size
            for _ in 1:8
                i1, i2 = rand(rng, 1:N), rand(rng, 1:N)
                P_host_new[i1, b], P_host_new[i2, b] = P_host_new[i2, b], P_host_new[i1, b]
            end
        end
        
        copyto!(P_gpu, P_host_new)
        energies_gpu_new = CUDA.zeros(Float32, batch_size)
        
        @cuda threads=1024 blocks=batch_size shmem=1024*sizeof(Float32) minspm_eval_kernel!(F_gpu, D_gpu, P_gpu, energies_gpu_new, N)
        CUDA.synchronize()
        
        E_host_new = Array(energies_gpu_new)
        
        for b in 1:batch_size
            delta = E_host_new[b] - E_host[b]
            if delta < 0 || rand(rng) < exp(-delta / max(T, 1e-5))
                P_host[:, b] .= P_host_new[:, b]
                E_host[b] = E_host_new[b]
            end
            
            if E_host[b] < best_ever_cost
                best_ever_cost = E_host[b]
            end
        end
        
        T *= cooling_rate
        if iter % 250 == 0 || iter == 1
            println("Iteration $iter | Current Best Cost: $best_ever_cost | Temperature: $(round(T, digits=2))")
        end
    end

    println("\n--- Final QAPLIB Results ---")
    println("MINSPM Best Found Cost     : $best_ever_cost")
    println("QAPLIB Best Known UB       : $known_ub")
    
    gap = ((best_ever_cost - known_ub) / known_ub) * 100.0
    println("Optimality Gap vs QAPLIB UB: $(round(gap, digits=2))%")

    CUDA.reclaim()
    GC.gc()
end

run_true_tai256c_optimization()
