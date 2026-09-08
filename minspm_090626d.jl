using CUDA
using LinearAlgebra
using Random
using DataStructures

# =====================================================================
# 1. MINSPM GALOIS-ENCODED BOUNDING KERNEL FOR B&B
# =====================================================================
function minspm_galois_bnb_kernel!(B_dict, Fixed_States, Bounds_Out, N, k_target, batch_size)
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
        
        state_i = Fixed_States[i, batch_id]
        state_j = Fixed_States[j, batch_id]
        
        xi = (state_i == 2) ? 1.0f0 : ((state_i == 1) ? 0.0f0 : (Float32(k_target) / Float32(N)))
        xj = (state_j == 2) ? 1.0f0 : ((state_j == 1) ? 0.0f0 : (Float32(k_target) / Float32(N)))
        
        val += B_dict[i, j] * xi * xj
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
        Bounds_Out[batch_id] = shared_vals[1]
    end
    return nothing
end

# =====================================================================
# 2. AUTOMORPHISM ORBIT PRUNING
# =====================================================================
function extract_orbits(B::Matrix{Float64}, N::Int)
    row_sigs = [sort(abs.(B[i, :]), rev=true) for i in 1:N]
    unique_sigs = unique(row_sigs)
    orbits = [findall(s -> row_sigs[s] == sig, 1:N) for sig in unique_sigs]
    filter!(o -> length(o) > 1, orbits)
    return orbits
end

function get_orbit_signature(fixed_state::Vector{Int8}, orbits::Vector{Vector{Int}})
    sig = UInt64(0)
    for orbit in orbits
        sig = hash(sum(fixed_state[orbit]), sig)
    end
    return sig
end

# =====================================================================
# 3. GPU BRANCH-AND-BOUND ENGINE
# =====================================================================
function load_qaplib(filepath::String)
    open(filepath, "r") do io
        tokens = split(read(io, String))
        N = parse(Int, tokens[1])
        idx = 2
        F = zeros(Float64, N, N)
        for i in 1:N, j in 1:N
            F[i, j] = parse(Float64, tokens[idx]); idx += 1
        end
        D = zeros(Float64, N, N)
        for i in 1:N, j in 1:N
            D[i, j] = parse(Float64, tokens[idx]); idx += 1
        end
        return N, F .+ D
    end
end

struct BnBNode
    fixed_states::Vector{Int8} # 0: free, 1: fixed 0, 2: fixed 1
    depth::Int
    bound::Float64
end

Base.isless(a::BnBNode, b::BnBNode) = a.bound < b.bound

function run_galois_gpu_bnb()
    filepath = "tai256c.dat"
    if !isfile(filepath)
        error("Please place 'tai256c.dat' in the working directory.")
    end

    N, B_cpu = load_qaplib(filepath)
    k_target = 92
    upper_bound = 44_759_292.0 
    lower_bound_ref = Ref(44_095_032.0) # Use a mutable Ref to safely update across scopes

    println("Initializing MINSPM Galois-GPU Branch-and-Bound Engine for tai256c...")
    orbits = extract_orbits(B_cpu, N)
    B_gpu = CuArray(Float32.(B_cpu))

    pq = BinaryMinHeap{BnBNode}()
    push!(pq, BnBNode(zeros(Int8, N), 0, lower_bound_ref[]))

    expanded = 0
    pruned = 0
    visited_orbits = Set{UInt64}()
    
    batch_size = 512
    node_batch = BnBNode[]

    println("Traversing B&B search space with Galois orbit pruning...")

    while !isempty(pq) && expanded < 25_000
        node = pop!(pq)

        if node.bound >= upper_bound
            pruned += 1
            continue
        end

        orbit_sig = get_orbit_signature(node.fixed_states, orbits)
        if orbit_sig in visited_orbits
            pruned += 1
            continue
        end
        push!(visited_orbits, orbit_sig)

        push!(node_batch, node)
        expanded += 1

        if length(node_batch) >= batch_size || isempty(pq)
            b_sz = length(node_batch)
            states_host = zeros(Int8, N, b_sz)
            for (idx, nd) in enumerate(node_batch)
                states_host[:, idx] .= nd.fixed_states
            end

            states_gpu = CuArray(states_host)
            bounds_gpu = CUDA.zeros(Float32, b_sz)

            @cuda threads=256 blocks=b_sz shmem=1024*sizeof(Float32) minspm_galois_bnb_kernel!(
                B_gpu, states_gpu, bounds_gpu, N, k_target, b_sz
            )
            CUDA.synchronize()

            bounds_host = Array(bounds_gpu)

            for (idx, nd) in enumerate(node_batch)
                b_val = Float64(bounds_host[idx])
                
                if b_val > lower_bound_ref[]
                    lower_bound_ref[] = b_val
                    println("Expanded $expanded | New Active Lower Bound: $(lower_bound_ref[])")
                end

                if b_val < upper_bound && nd.depth < N
                    free_var = findfirst(==(0), nd.fixed_states)
                    if !isnothing(free_var)
                        s0 = copy(nd.fixed_states)
                        s0[free_var] = 1
                        push!(pq, BnBNode(s0, nd.depth + 1, b_val))

                        s1 = copy(nd.fixed_states)
                        s1[free_var] = 2
                        push!(pq, BnBNode(s1, nd.depth + 1, b_val))
                    end
                end
            end
            empty!(node_batch)
        end
    end

    println("\n--- Galois GPU Branch-and-Bound Complete ---")
    println("Nodes Expanded : $expanded")
    println("Nodes Pruned   : $pruned")
    println("Final B&B Bound: $(lower_bound_ref[])")
    CUDA.reclaim()
end

run_galois_gpu_bnb()
