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

function run_orbit_guided_gpu_solver()
    filepath = "tai256c.dat"
    if !isfile(filepath)
        error("Please place 'tai256c.dat' in the working directory.")
    end

    N, F_cpu, D_cpu = load_qaplib_tai256c(filepath)
    batch_size = 16_384
    max_generations = 6_000
    known_ub = 44_759_296.0f0

    println("Initializing Orbit-Guided Symmetry-Aware Solver for tai256c [N = $N]...")
    
    F_gpu = CuArray(F_cpu)
    D_gpu = CuArray(D_cpu)
    rng = MersenneTwister(2026)
    
    # Compute symmetry orbits using row-sum equivalence classes of B = F + D
    B_matrix = F_cpu .+ D_cpu
    row_signatures = round.(vec(sum(B_matrix, dims=2)), digits=-2)
    unique_sigs = unique(row_signatures)
    
    orbits = [findall(==(sig), row_signatures) for sig in unique_sigs]
    filter!(o -> length(o) > 1, orbits)

    println("Identified $(length(orbits)) symmetry-induced structural orbits for block moves.")

    base_order = sortperm(vec(sum(F_cpu, dims=2)))
    
    P_host = zeros(Int32, N, batch_size)
    for b in 1:batch_size
        if b <= 4_096
            p = copy(base_order)
            for orbit in orbits
                if length(orbit) >= 2 && rand(rng) < 0.5
                    o1, o2 = rand(rng, orbit), rand(rng, orbit)
                    p[o1], p[o2] = p[o2], p[o1]
                end
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

    elite_pool = [copy(global_best_perm) for _ in 1:64]
    elite_costs = fill(global_best_cost, 64)

    println("Executing Orbit-Structured Block Search & Path Relinking via MINSPM...")
    
    for gen in 1:max_generations
        P_host_new = copy(P_host)
        
        for b in 1:batch_size
            if rand(rng) < 0.7 && !isempty(orbits)
                orbit = orbits[rand(rng, 1:length(orbits))]
                if length(orbit) >= 2
                    i1, i2 = rand(rng, orbit), rand(rng, orbit)
                    P_host_new[i1, b], P_host_new[i2, b] = P_host_new[i2, b], P_host_new[i1, b]
                end
            else
                for _ in 1:4
                    i1, i2 = rand(rng, 1:N), rand(rng, 1:N)
                    P_host_new[i1, b], P_host_new[i2, b] = P_host_new[i2, b], P_host_new[i1, b]
                end
            end
        end
        
        copyto!(P_gpu, P_host_new)
        energies_gpu_new = CUDA.zeros(Float32, batch_size)
        
        @cuda threads=1024 blocks=batch_size shmem=1024*sizeof(Float32) minspm_qap_eval_kernel!(F_gpu, D_gpu, P_gpu, energies_gpu_new, N)
        CUDA.synchronize()
        
        E_host_new = Array(energies_gpu_new)
        
        for b in 1:batch_size
            if E_host_new[b] < E_host[b]
                P_host[:, b] .= P_host_new[:, b]
                E_host[b] = E_host_new[b]
                
                if E_host[b] < global_best_cost
                    global_best_cost = E_host[b]
                    global_best_perm .= P_host[:, b]
                    
                    worst_idx = argmax(elite_costs)
                    if E_host[b] < elite_costs[worst_idx]
                        elite_costs[worst_idx] = E_host[b]
                        elite_pool[worst_idx] .= P_host[:, b]
                    end
                end
            elseif rand(rng) < 0.05
                target_elite = elite_pool[rand(rng, 1:length(elite_pool))]
                mismatches = findall(i -> P_host[i, b] != target_elite[i], 1:N)
                if length(mismatches) >= 2
                    m1 = mismatches[1]
                    curr_pos = findfirst(==(target_elite[m1]), P_host[:, b])
                    if !isnothing(curr_pos)
                        P_host[m1, b], P_host[curr_pos, b] = P_host[curr_pos, b], P_host[m1, b]
                    end
                end
            end
        end

        if gen % 200 == 0 || gen == 1
            gap = ((global_best_cost - known_ub) / known_ub) * 100.0
            println("Generation $gen | Global Best Cost: $global_best_cost | Gap vs UB: $(round(gap, digits=3))%")
        end
    end

    final_gap = ((global_best_cost - known_ub) / known_ub) * 100.0
    println("\n--- Orbit-Guided Optimization Summary ---")
    println("Global Best Cost Achieved  : $global_best_cost")
    println("QAPLIB Best Known UB       : $known_ub")
    println("Optimality Gap vs UB       : $(round(final_gap, digits=3))%")

    CUDA.reclaim()
    GC.gc()
end

run_orbit_guided_gpu_solver()
