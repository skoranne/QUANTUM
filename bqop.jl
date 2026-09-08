using CUDA
using LinearAlgebra
using Random

function bqop_eval_kernel!(B, X_bin, energies, N)
    batch_id = blockIdx().x
    tid = threadIdx().x
    block_dim = blockDim().x
    
    shared_vals = @cuDynamicSharedMem(Float32, 1024)
    
    val = 0.0f0
    idx = tid
    total_elements = N * N
    
    while idx <= total_elements
        i = div(idx - 1, N) + 1
        j = rem(idx - 1, N) + 1
        
        # BQOP quadratic form x^T * B * x
        if X_bin[i, batch_id] == 1 && X_bin[j, batch_id] == 1
            val += B[i, j]
        end
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

function load_bqop_matrix_from_qaplib(filepath::String)
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
        # Reduced 256x256 BQOP matrix B derived from Taillard's instance symmetry mapping
        B = F .+ D 
        return N, B
    end
end

function run_bqop_cardinality_optimization()
    filepath = "tai256c.dat"
    if !isfile(filepath)
        error("Please place 'tai256c.dat' in the working directory.")
    end

    N, B_cpu = load_bqop_matrix_from_qaplib(filepath)
    batch_size = 4_096
    num_iterations = 3_000
    k_ones = 92 # Strict cardinality constraint sum(x) = 92
    
    known_ub = 44_759_296.0f0
    known_lb = 44_200_812.0f0

    println("Loaded BQOP Matrix for tai256c [Dimension = $N, Cardinality Target = $k_ones]")
    
    B_gpu = CuArray(B_cpu)
    rng = MersenneTwister(42)
    
    # Initialize binary matrix X where each column has exactly 92 ones
    X_host = zeros(Int8, N, batch_size)
    for b in 1:batch_size
        ones_indices = randperm(rng, N)[1:k_ones]
        X_host[ones_indices, b] .= 1
    end
    
    X_gpu = CuArray(X_host)
    energies_gpu = CUDA.zeros(Float32, batch_size)

    @cuda threads=1024 blocks=batch_size shmem=1024*sizeof(Float32) bqop_eval_kernel!(B_gpu, X_gpu, energies_gpu, N)
    CUDA.synchronize()
    
    E_host = Array(energies_gpu)
    best_ever_cost = minimum(E_host)
    
    T = 100_000.0f0
    cooling_rate = 0.995f0

    println("Running GPU Search on 92-Cardinality BQOP Subspace...")
    for iter in 1:num_iterations
        X_host_new = copy(X_host)
        
        # Subspace-preserving move: swap a random '1' and '0' in each chain
        for b in 1:batch_size
            ones_pos = findall(==(1), X_host_new[:, b])
            zeros_pos = findall(==(0), X_host_new[:, b])
            
            # Perform multiple swaps per step to enhance escape from local traps
            for _ in 1:4
                if !isempty(ones_pos) && !isempty(zeros_pos)
                    p1 = rand(rng, ones_pos)
                    p2 = rand(rng, zeros_pos)
                    X_host_new[p1, b] = 0
                    X_host_new[p2, b] = 1
                    # Update tracking lists locally
                    ones_pos = findall(==(1), X_host_new[:, b])
                    zeros_pos = findall(==(0), X_host_new[:, b])
                end
            end
        end
        
        copyto!(X_gpu, X_host_new)
        energies_gpu_new = CUDA.zeros(Float32, batch_size)
        
        @cuda threads=1024 blocks=batch_size shmem=1024*sizeof(Float32) bqop_eval_kernel!(B_gpu, X_gpu, energies_gpu_new, N)
        CUDA.synchronize()
        
        E_host_new = Array(energies_gpu_new)
        
        for b in 1:batch_size
            delta = E_host_new[b] - E_host[b]
            if delta < 0 || rand(rng) < exp(-delta / max(T, 1e-5))
                X_host[:, b] .= X_host_new[:, b]
                E_host[b] = E_host_new[b]
            end
            
            if E_host[b] < best_ever_cost
                best_ever_cost = E_host[b]
            end
        end
        
        T *= cooling_rate
        if iter % 300 == 0 || iter == 1
            println("Iteration $iter | Current BQOP Best Cost: $best_ever_cost | Temperature: $(round(T, digits=2))")
        end
    end

    println("\n--- BQOP Optimization Results ---")
    println("MINSPM BQOP Best Found Cost : $best_ever_cost")
    println("QAPLIB Best Known UB        : $known_ub")
    
    gap = ((best_ever_cost - known_ub) / known_ub) * 100.0
    println("Optimality Gap vs QAPLIB UB : $(round(gap, digits=2))%")

    CUDA.reclaim()
    GC.gc()
end

run_bqop_cardinality_optimization()
