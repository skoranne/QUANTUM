using CUDA
using LinearAlgebra
using Random
using Printf

# ============================================================
# TAI256C
# Symmetry-aware GPU BQOP / QAP solver
# ============================================================

const INSTANCE_FILE = "tai256c.dat"

const N_EXPECTED = 256
const K_TARGET = 92

# QAPLIB best-known/reference value.
# This is ONLY used for reporting.
const KNOWN_BEST = Int64(44_759_294)

# ------------------------------------------------------------
# Search parameters
# ------------------------------------------------------------

const POPULATION = 16_384
const MAX_EPOCHS = 2_000

const MAX_1SWAP_MOVES = 500

const KICK_SIZES = (2, 4, 8, 16)

const ARCHIVE_INTERVAL = 25
const ELITE_SIZE = 64

const EXPECTED_GROUP_ORDER = 2048

const RNG_SEED = 20260906


# ============================================================
# QAP FILE LOADER
# ============================================================

function load_tai256c(filepath::String)

    if !isfile(filepath)
        error(
            "Cannot find '$filepath'. " *
            "Place tai256c.dat in the current directory."
        )
    end

    tokens = split(read(filepath, String))

    pos = 1

    n = parse(Int, tokens[pos])
    pos += 1

    F = Matrix{Int64}(undef, n, n)

    @inbounds for i in 1:n
        for j in 1:n
            F[i,j] = parse(Int64, tokens[pos])
            pos += 1
        end
    end

    D = Matrix{Int64}(undef, n, n)

    @inbounds for i in 1:n
        for j in 1:n
            D[i,j] = parse(Int64, tokens[pos])
            pos += 1
        end
    end

    if !issymmetric(F)
        error("F is not symmetric.")
    end

    if !issymmetric(D)
        error("D is not symmetric.")
    end

    B = F + D

    return n, F, D, B
end


# ============================================================
# EXACT OBJECTIVE
# ============================================================

function objective_int(
    B::Matrix{Int64},
    x::Vector{Int8}
)

    n = size(B,1)

    result = Int64(0)

    @inbounds for i in 1:n

        if x[i] != 0

            for j in 1:n

                if x[j] != 0
                    result += B[i,j]
                end

            end
        end
    end

    return result
end


# ============================================================
# WEIGHTED GRAPH COLOR REFINEMENT
# ============================================================

function initial_weight_colors(
    B::Matrix{Int64}
)

    n = size(B,1)

    dictionary = Dict{Any,Int}()

    colors = zeros(Int,n)

    next_color = 0

    for v in 1:n

        row = Vector{Int64}(undef,n)

        @inbounds for u in 1:n
            row[u] = B[v,u]
        end

        sort!(row)

        key = (
            B[v,v],
            Tuple(row)
        )

        c = get(dictionary,key,0)

        if c == 0
            next_color += 1
            c = next_color
            dictionary[key] = c
        end

        colors[v] = c
    end

    return colors
end


function refine_weight_colors(
    B::Matrix{Int64},
    colors::Vector{Int};
    max_rounds::Int=50
)

    n = size(B,1)

    current = copy(colors)

    for round in 1:max_rounds

        dictionary = Dict{Any,Int}()

        newcolors = zeros(Int,n)

        next_color = 0

        for v in 1:n

            items =
                Vector{
                    Tuple{Int64,Int}
                }(undef,n)

            @inbounds for u in 1:n
                items[u] =
                    (
                        B[v,u],
                        current[u]
                    )
            end

            sort!(items)

            key = (
                B[v,v],
                Tuple(items)
            )

            c = get(dictionary,key,0)

            if c == 0
                next_color += 1
                c = next_color
                dictionary[key] = c
            end

            newcolors[v] = c
        end

        if newcolors == current
            return current
        end

        current = newcolors
    end

    return current
end


function make_color_classes(
    colors::Vector{Int}
)

    classes =
        Dict{Int,Vector{Int}}()

    for i in eachindex(colors)

        if !haskey(classes,colors[i])
            classes[colors[i]] = Int[]
        end

        push!(
            classes[colors[i]],
            i
        )
    end

    return collect(values(classes))
end


# ============================================================
# EXACT AUTOMORPHISM ENUMERATION
# ============================================================

function enumerate_automorphisms(
    B::Matrix{Int64}
)

    n = size(B,1)

    println()
    println("============================================================")
    println("COMPUTING EXACT AUTOMORPHISM GROUP OF B")
    println("============================================================")

    initial =
        initial_weight_colors(B)

    colors =
        refine_weight_colors(
            B,
            initial
        )

    classes =
        make_color_classes(colors)

    println(
        "Refined color classes = ",
        length(classes)
    )

    println(
        "Class sizes = ",
        sort(length.(classes))
    )

    candidates_by_color =
        Dict{Int,Vector{Int}}()

    for cls in classes
        c = colors[cls[1]]
        candidates_by_color[c] = cls
    end

    candidates =
        Vector{Vector{Int}}(undef,n)

    for i in 1:n
        candidates[i] =
            candidates_by_color[
                colors[i]
            ]
    end

    mapping =
        zeros(Int,n)

    inverse_mapping =
        zeros(Int,n)

    automorphisms =
        Vector{Vector{Int16}}()

    nodes =
        Ref{Int64}(0)

    function choose_vertex()

        best = 0
        best_count = typemax(Int)

        for v in 1:n

            if mapping[v] != 0
                continue
            end

            count = 0

            for w in candidates[v]

                if inverse_mapping[w] == 0
                    count += 1
                end
            end

            if count < best_count

                best_count = count
                best = v

                if count <= 1
                    break
                end
            end
        end

        return best
    end


    function recurse(assigned::Int)

        nodes[] += 1

        if assigned == n

            p =
                Vector{Int16}(undef,n)

            @inbounds for i in 1:n
                p[i] =
                    Int16(mapping[i])
            end

            push!(
                automorphisms,
                p
            )

            if length(automorphisms) % 256 == 0

                println(
                    "  automorphisms = ",
                    length(automorphisms),
                    " | nodes = ",
                    nodes[]
                )
            end

            return
        end

        v =
            choose_vertex()

        if v == 0
            error(
                "Automorphism search failed to select vertex."
            )
        end

        for w in candidates[v]

            if inverse_mapping[w] != 0
                continue
            end

            compatible = true

            @inbounds for u in 1:n

                pu = mapping[u]

                if pu != 0

                    if B[v,u] != B[w,pu]

                        compatible = false
                        break
                    end
                end
            end

            if !compatible
                continue
            end

            mapping[v] = w
            inverse_mapping[w] = v

            recurse(assigned + 1)

            mapping[v] = 0
            inverse_mapping[w] = 0
        end
    end

    recurse(0)

    println()
    println(
        "Automorphism enumeration complete."
    )

    println(
        "Group order = ",
        length(automorphisms)
    )

    println(
        "Expected order = ",
        EXPECTED_GROUP_ORDER
    )

    return automorphisms
end


# ============================================================
# AUTOMORPHISM VERIFICATION
# ============================================================

function verify_automorphism(
    B::Matrix{Int64},
    p::AbstractVector{<:Integer}
)

    n = size(B,1)

    seen = falses(n)

    @inbounds for i in 1:n

        q = Int(p[i])

        if q < 1 || q > n
            return false
        end

        if seen[q]
            return false
        end

        seen[q] = true
    end

    @inbounds for i in 1:n

        for j in 1:n

            if B[i,j] != B[p[i],p[j]
]
                return false
            end
        end
    end

    return true
end


function verify_automorphism_group(
    B::Matrix{Int64},
    group
)

    println()
    println("Verifying automorphisms exactly...")

    for g in eachindex(group)

        if !verify_automorphism(
            B,
            group[g]
        )

            error(
                "Invalid automorphism at index ",
                g
            )
        end

        if g % 256 == 0

            println(
                "  verified ",
                g,
                " / ",
                length(group)
            )
        end
    end

    println(
        "✓ Every automorphism satisfies P'BP = B."
    )

    return true
end


# ============================================================
# SAVE AUTOMORPHISMS
# ============================================================

function save_automorphisms(
    filename,
    group
)

    open(filename,"w") do io

        println(
            io,
            length(group)
        )

        println(
            io,
            length(group[1])
        )

        for p in group

            for i in eachindex(p)

                print(
                    io,
                    p[i]
                )

                if i != length(p)
                    print(io,' ')
                end
            end

            println(io)
        end
    end
end


# ============================================================
# ORBITS
# ============================================================

function compute_orbits(
    group,
    n
)

    seen = falses(n)

    orbits =
        Vector{Vector{Int}}()

    for v in 1:n

        if seen[v]
            continue
        end

        orbit = Int[]

        for p in group
            push!(
                orbit,
                Int(p[v])
            )
        end

        unique!(sort!(orbit))

        for u in orbit
            seen[u] = true
        end

        push!(
            orbits,
            orbit
        )
    end

    return orbits
end


# ============================================================
# INITIAL POPULATION
# ============================================================

function initialize_population(
    n,
    k,
    group,
    population,
    rng
)

    X =
        zeros(
            Float64,
            n,
            population
        )

    ng =
        length(group)

    # First member.
    X[1:k,1] .= 1.0

    for b in 2:population

        idx =
            randperm(rng,n)[1:k]

        @inbounds for i in idx
            X[i,b] = 1.0
        end

        # Exact symmetry transformation.
        if rand(rng) < 0.5

            p =
                group[
                    rand(rng,1:ng)
                ]

            tmp =
                copy(
                    @view X[:,b]
                )

            @inbounds for i in 1:n
                X[i,b] =
                    tmp[p[i]]
            end
        end
    end

    return X
end


# ============================================================
# GPU ENERGY KERNEL
# ============================================================

function energy_kernel!(
    X,
    H,
    E,
    n
)

    b =
        blockIdx().x

    tid =
        threadIdx().x

    shared =
        @cuDynamicSharedMem(
            Float64,
            256
        )

    value = 0.0

    i = tid

    while i <= n

        @inbounds begin
            value +=
                X[i,b] *
                H[i,b]
        end

        i += blockDim().x
    end

    shared[tid] = value

    sync_threads()

    stride = 128

    while stride >= 1

        if tid <= stride

            @inbounds begin
                shared[tid] +=
                    shared[tid + stride]
            end
        end

        sync_threads()

        stride >>= 1
    end

    if tid == 1
        E[b] = shared[1]
    end

    return nothing
end


# ============================================================
# GPU BEST 1-SWAP
# ============================================================

function best_swap_kernel!(
    B,
    X,
    H,
    best_p,
    best_q,
    best_delta,
    n
)

    b =
        blockIdx().x

    tid =
        threadIdx().x

    shared_delta =
        @cuDynamicSharedMem(
            Float64,
            256
        )

    shared_p =
        @cuDynamicSharedMem(
            Int32,
            256,
            256 * sizeof(Float64)
        )

    shared_q =
        @cuDynamicSharedMem(
            Int32,
            256,
            256 * sizeof(Float64) +
            256 * sizeof(Int32)
        )

    local_delta = Inf
    local_p = Int32(0)
    local_q = Int32(0)

    p = tid

    if p <= n

        @inbounds xp = X[p,b]

        if xp > 0.5

            @inbounds hp = H[p,b]

            q = 1

            while q <= n

                @inbounds xq = X[q,b]

                if xq < 0.5

                    @inbounds begin

                        delta =
                            2.0 *
                            (
                                H[q,b] -
                                hp
                            ) +
                            B[p,p] +
                            B[q,q] -
                            2.0 * B[p,q]

                        if delta < local_delta

                            local_delta = delta
                            local_p = Int32(p)
                            local_q = Int32(q)

                        end
                    end
                end

                q += 1
            end
        end
    end

    shared_delta[tid] = local_delta
    shared_p[tid] = local_p
    shared_q[tid] = local_q

    sync_threads()

    stride = 128

    while stride >= 1

        if tid <= stride

            other = tid + stride

            if shared_delta[other] <
               shared_delta[tid]

                shared_delta[tid] =
                    shared_delta[other]

                shared_p[tid] =
                    shared_p[other]

                shared_q[tid] =
                    shared_q[other]
            end
        end

        sync_threads()

        stride >>= 1
    end

    if tid == 1

        best_delta[b] =
            shared_delta[1]

        best_p[b] =
            shared_p[1]

        best_q[b] =
            shared_q[1]
    end

    return nothing
end


# ============================================================
# GPU APPLY SWAP
# ============================================================

function apply_swap_kernel!(
    B,
    X,
    H,
    E,
    best_p,
    best_q,
    best_delta,
    n
)

    b =
        blockIdx().x

    tid =
        threadIdx().x

    p =
        Int(best_p[b])

    q =
        Int(best_q[b])

    delta =
        best_delta[b]

    if p <= 0 || q <= 0
        return nothing
    end

    if delta >= 0.0
        return nothing
    end

    i = tid

    if i <= n

        @inbounds begin

            H[i,b] +=
                B[i,q] -
                B[i,p]

            if i == p
                X[i,b] = 0.0
            elseif i == q
                X[i,b] = 1.0
            end
        end
    end

    if tid == 1
        E[b] += delta
    end

    return nothing
end


# ============================================================
# GPU ONE-SWAP DESCENT
# ============================================================

function gpu_one_swap_descent!(
    B_gpu,
    X_gpu,
    H_gpu,
    E_gpu,
    best_p_gpu,
    best_q_gpu,
    best_delta_gpu,
    n,
    max_moves
)

    population =
        size(X_gpu,2)

    shared_bytes =
        256 * sizeof(Float64) +
        256 * sizeof(Int32) +
        256 * sizeof(Int32)

    total_moves = 0

    for iteration in 1:max_moves

        @cuda threads=256 blocks=population shmem=shared_bytes best_swap_kernel!(
            B_gpu,
            X_gpu,
            H_gpu,
            best_p_gpu,
            best_q_gpu,
            best_delta_gpu,
            n
        )

        CUDA.synchronize()

        deltas =
            Array(
                best_delta_gpu
            )

        improving = false

        @inbounds for b in 1:population

            if deltas[b] < -1e-9

                improving = true
                break
            end
        end

        if !improving
            break
        end

        @cuda threads=256 blocks=population apply_swap_kernel!(
            B_gpu,
            X_gpu,
            H_gpu,
            E_gpu,
            best_p_gpu,
            best_q_gpu,
            best_delta_gpu,
            n
        )

        CUDA.synchronize()

        total_moves += 1
    end

    return total_moves
end


# ============================================================
# GPU RANDOM STATE
# ============================================================

@inline function xorshift64_device(
    x::UInt64
)

    x ⊻=
        x << 13

    x ⊻=
        x >> 7

    x ⊻=
        x << 17

    return x
end


@inline function random_index_device(
    state::UInt64,
    n::Int
)

    state =
        xorshift64_device(state)

    idx =
        Int(
            mod(
                state,
                UInt64(n)
            )
        ) + 1

    return state,idx
end


# ============================================================
# GPU VNS KICK
# ============================================================

function vns_kick_kernel!(
    B,
    X,
    H,
    states,
    kick_size,
    n
)

    b =
        blockIdx().x

    tid =
        threadIdx().x

    shared_p =
        @cuStaticSharedMem(
            Int32,
            1
        )

    shared_q =
        @cuStaticSharedMem(
            Int32,
            1
        )

    state =
        states[b]

    for k in 1:kick_size

        if tid == 1

            p = 0
            q = 0

            for attempt in 1:256

                state,candidate =
                    random_index_device(
                        state,
                        n
                    )

                if X[candidate,b] > 0.5

                    p = candidate
                    break
                end
            end

            for attempt in 1:256

                state,candidate =
                    random_index_device(
                        state,
                        n
                    )

                if X[candidate,b] < 0.5

                    q = candidate
                    break
                end
            end

            shared_p[1] =
                Int32(p)

            shared_q[1] =
                Int32(q)
        end

        sync_threads()

        p =
            Int(shared_p[1])

        q =
            Int(shared_q[1])

        if p > 0 && q > 0

            i = tid

            if i <= n

                @inbounds begin

                    H[i,b] +=
                        B[i,q] -
                        B[i,p]

                    if i == p
                        X[i,b] = 0.0
                    elseif i == q
                        X[i,b] = 1.0
                    end
                end
            end
        end

        sync_threads()
    end

    states[b] =
        state

    return nothing
end


# ============================================================
# GPU SYMMETRY JUMP
# ============================================================

function symmetry_jump_kernel!(
    X,
    Xtmp,
    perm_table,
    choices,
    n
)

    b =
        blockIdx().x

    i =
        threadIdx().x

    if i <= n

        g =
            Int(
                choices[b]
            )

        @inbounds begin

            source =
                Int(
                    perm_table[i,g]
                )

            Xtmp[i,b] =
                X[source,b]
        end
    end

    return nothing
end


# ============================================================
# CPU GRADIENT
# ============================================================

function gradient_cpu(
    B::Matrix{Int64},
    x::Vector{Float64}
)

    return Float64.(B) * x
end


# ============================================================
# CPU PATH RELINKING
# ============================================================

function hamming_distance(
    x,
    y
)

    d = 0

    @inbounds for i in eachindex(x)

        if x[i] != y[i]
            d += 1
        end
    end

    return d
end


function path_relink(
    B::Matrix{Int64},
    start::Vector{Float64},
    target::Vector{Float64}
)

    x =
        copy(start)

    h =
        gradient_cpu(
            B,
            x
        )

    best_x =
        copy(x)

    best_e =
        dot(x,h)

    n =
        length(x)

    for step in 1:n

        best_delta =
            Inf

        best_p = 0
        best_q = 0

        @inbounds for p in 1:n

            if x[p] > 0.5 &&
               target[p] < 0.5

                for q in 1:n

                    if x[q] < 0.5 &&
                       target[q] > 0.5

                        delta =
                            2.0 *
                            (
                                h[q] -
                                h[p]
                            ) +
                            B[p,p] +
                            B[q,q] -
                            2.0 * B[p,q]

                        if delta < best_delta

                            best_delta =
                                delta

                            best_p = p
                            best_q = q
                        end
                    end
                end
            end
        end

        if best_p == 0
            break
        end

        x[best_p] = 0.0
        x[best_q] = 1.0

        @inbounds for i in 1:n

            h[i] +=
                B[i,best_q] -
                B[i,best_p]
        end

        e =
            dot(x,h)

        if e < best_e

            best_e = e
            best_x .= x
        end
    end

    return best_e,best_x
end


# ============================================================
# ELITE ARCHIVE
# ============================================================

function update_archive!(
    archive,
    X,
    E
)

    order =
        sortperm(E)

    limit =
        min(
            length(order),
            ELITE_SIZE
        )

    for pos in 1:limit

        b =
            order[pos]

        e =
            Int64(
                round(
                    E[b]
                )
            )

        x =
            copy(
                @view X[:,b]
            )

        duplicate = false

        for (ae,ax) in archive

            if ae == e &&
               hamming_distance(ax,x) == 0

                duplicate = true
                break
            end
        end

        if !duplicate

            push!(
                archive,
                (e,x)
            )
        end
    end

    sort!(
        archive,
        by = z -> z[1]
    )

    if length(archive) > ELITE_SIZE

        resize!(
            archive,
            ELITE_SIZE
        )
    end

    return nothing
end


# ============================================================
# EXACT CANDIDATE VALIDATION
# ============================================================

function validate_candidate(
    B,
    x
)

    xi =
        Int8.(
            round.(x)
        )

    cardinality =
        sum(xi)

    if cardinality != K_TARGET

        error(
            "Invalid candidate cardinality ",
            cardinality
        )
    end

    cost =
        objective_int(
            B,
            xi
        )

    return cost,xi
end


# ============================================================
# MAIN
# ============================================================

function run_tai256c()

    println()
    println("============================================================")
    println("TAI256C GPU SYMMETRY-AWARE BQOP SOLVER")
    println("============================================================")

    n,F,D,B_int =
        load_tai256c(
            INSTANCE_FILE
        )

    println()
    println("Loaded ",INSTANCE_FILE)
    println("N = ",n)
    println("Target cardinality = ",K_TARGET)

    if n != N_EXPECTED

        error(
            "Expected N=$N_EXPECTED, got $n."
        )
    end

    if !issymmetric(B_int)

        error(
            "B is not symmetric."
        )
    end

    println("B = F + D constructed exactly with Int64.")

    # ========================================================
    # AUTOMORPHISMS
    # ========================================================

    group =
        enumerate_automorphisms(
            B_int
        )

    verify_automorphism_group(
        B_int,
        group
    )

    if length(group) == 1
	        println("WARNING: Only the identity automorphism was found.")
	elseif length(group) == 2
		    println("Found a nontrivial involution; group order = 2.")
	    else
		        println(
				        "Found automorphism group of order ",
					        length(group),
						        "."
							    )
						    end


    save_automorphisms(
        "tai256c_automorphisms.txt",
        group
    )

    println(
        "✓ Saved exact automorphism group."
    )

    orbits =
        compute_orbits(
            group,
            n
        )

    println()
    println(
        "Number of vertex orbits = ",
        length(orbits)
    )

    println(
        "Orbit sizes = ",
        sort(length.(orbits))
    )

    # ========================================================
    # GPU
    # ========================================================

    println()
    println("Initializing CUDA...")

    CUDA.versioninfo()

    B_gpu =
        CuArray(
            Float64.(B_int)
        )

    # ========================================================
    # POPULATION
    # ========================================================

    rng =
        MersenneTwister(
            RNG_SEED
        )

    println()
    println(
        "Generating ",
        POPULATION,
        " feasible solutions..."
    )

    X_host =
        initialize_population(
            n,
            K_TARGET,
            group,
            POPULATION,
            rng
        )

    # ========================================================
    # GPU STATE
    # ========================================================

    X_gpu =
        CuArray(X_host)

    H_gpu =
        CUDA.zeros(
            Float64,
            n,
            POPULATION
        )

    E_gpu =
        CUDA.zeros(
            Float64,
            POPULATION
        )

    best_p_gpu =
        CUDA.zeros(
            Int32,
            POPULATION
        )

    best_q_gpu =
        CUDA.zeros(
            Int32,
            POPULATION
        )

    best_delta_gpu =
        CUDA.zeros(
            Float64,
            POPULATION
        )

    # ========================================================
    # RNG STATES
    # ========================================================

    states_host =
        Vector{UInt64}(
            undef,
            POPULATION
        )

    for b in 1:POPULATION

        states_host[b] =
            UInt64(RNG_SEED) ⊻
            (
                UInt64(b) *
                UInt64(
                    0x9e3779b97f4a7c15
                )
            )
    end

    states_gpu =
        CuArray(states_host)

    # ========================================================
    # INITIAL H
    # ========================================================

    println()
    println("Computing initial B*X on GPU...")

    H_gpu .=
        B_gpu * X_gpu

    @cuda threads=256 blocks=POPULATION shmem=256*sizeof(Float64) energy_kernel!(
        X_gpu,
        H_gpu,
        E_gpu,
        n
    )

    CUDA.synchronize()

    E_host =
        Array(E_gpu)

    initial_idx =
        argmin(E_host)

    initial_candidate =
        Array(
            @view X_gpu[:,initial_idx]
        )

    initial_cost,
    initial_x =
        validate_candidate(
            B_int,
            initial_candidate
        )

    global_best_cost =
        initial_cost

    global_best_x =
        copy(
            initial_candidate
        )

    println()
    println(
        "Initial exact cost = ",
        initial_cost
    )

    println(
        "Reference best     = ",
        KNOWN_BEST
    )

    println(
        "Initial gap        = ",
        initial_cost - KNOWN_BEST
    )

    # ========================================================
    # SYMMETRY TABLE
    # ========================================================

    ng =
        length(group)

    perm_table_host =
        Matrix{Int32}(
            undef,
            n,
            ng
        )

    for g in 1:ng

        p =
            group[g]

        for i in 1:n

            perm_table_host[i,g] =
                Int32(p[i])
        end
    end

    perm_table_gpu =
        CuArray(
            perm_table_host
        )

    choices_host =
        Vector{Int32}(
            undef,
            POPULATION
        )

    for b in 1:POPULATION

        choices_host[b] =
            Int32(
                rand(
                    rng,
                    1:ng
                )
            )
    end

    choices_gpu =
        CuArray(
            choices_host
        )

    Xtmp_gpu =
        CUDA.zeros(
            Float64,
            n,
            POPULATION
        )

    # ========================================================
    # ARCHIVE
    # ========================================================

    archive =
        Tuple{
            Int64,
            Vector{Float64}
        }[]

    update_archive!(
        archive,
        X_host,
        E_host
    )

    total_moves = 0

    last_improvement = 0

    # ========================================================
    # SEARCH
    # ========================================================

    println()
    println("============================================================")
    println("STARTING SEARCH")
    println("============================================================")

    for epoch in 1:MAX_EPOCHS

        # ----------------------------------------------------
        # Local search.
        # ----------------------------------------------------

        moves =
            gpu_one_swap_descent!(
                B_gpu,
                X_gpu,
                H_gpu,
                E_gpu,
                best_p_gpu,
                best_q_gpu,
                best_delta_gpu,
                n,
                MAX_1SWAP_MOVES
            )

        total_moves += moves

        # ----------------------------------------------------
        # Retrieve best.
        # ----------------------------------------------------

        E_host =
            Array(E_gpu)

        idx =
            argmin(E_host)

        candidate =
            Array(
                @view X_gpu[:,idx]
            )

        candidate_cost,
        candidate_x =
            validate_candidate(
                B_int,
                candidate
            )

        if candidate_cost <
           global_best_cost

            global_best_cost =
                candidate_cost

            global_best_x =
                copy(candidate)

            last_improvement =
                epoch

            println()
            println(
                "Epoch ",
                epoch,
                " | NEW BEST = ",
                global_best_cost
            )

            println(
                "Gap vs reference = ",
                global_best_cost -
                KNOWN_BEST
            )
        end

        # ----------------------------------------------------
        # Archive / path relinking.
        # ----------------------------------------------------

        if epoch % ARCHIVE_INTERVAL == 0

            X_host =
                Array(X_gpu)

            update_archive!(
                archive,
                X_host,
                E_host
            )

            if length(archive) >= 2

                for attempt in 1:4

                    a =
                        rand(
                            rng,
                            1:length(archive)
                        )

                    b =
                        rand(
                            rng,
                            1:length(archive)
                        )

                    while b == a

                        b =
                            rand(
                                rng,
                                1:length(archive)
                            )
                    end

                    start =
                        archive[a][2]

                    target =
                        archive[b][2]

                    # Transform target by an exact automorphism.
                    p =
                        group[
                            rand(
                                rng,
                                1:ng
                            )
                        ]

                    target2 =
                        similar(target)

                    @inbounds for i in 1:n
                        target2[i] =
                            target[p[i]]
                    end

                    path_cost,
                    path_x =
                        path_relink(
                            B_int,
                            start,
                            target2
                        )

                    path_x_i8 =
                        Int8.(
                            round.(
                                path_x
                            )
                        )

                    exact_path_cost =
                        objective_int(
                            B_int,
                            path_x_i8
                        )

                    if exact_path_cost <
                       global_best_cost

                        global_best_cost =
                            exact_path_cost

                        global_best_x =
                            copy(path_x)

                        last_improvement =
                            epoch

                        println()
                        println(
                            "Epoch ",
                            epoch,
                            " | PATH RELINK NEW BEST = ",
                            global_best_cost
                        )

                        println(
                            "Gap vs reference = ",
                            global_best_cost -
                            KNOWN_BEST
                        )
                    end

                    # Inject path solution.
                    slot =
                        rand(
                            rng,
                            1:POPULATION
                        )

                    path_h =
                        gradient_cpu(
                            B_int,
                            path_x
                        )

                    X_gpu[:,slot] =
                        CuArray(path_x)

                    H_gpu[:,slot] =
                        CuArray(path_h)

                    E_gpu[slot] =
                        dot(
                            path_x,
                            path_h
                        )
                end
            end
        end

        # ----------------------------------------------------
        # VNS kick.
        # ----------------------------------------------------

        kick_index =
            1 +
            mod(
                (epoch - 1) ÷ 10,
                length(KICK_SIZES)
            )

        kick =
            KICK_SIZES[
                kick_index
            ]

        @cuda threads=256 blocks=POPULATION vns_kick_kernel!(
            B_gpu,
            X_gpu,
            H_gpu,
            states_gpu,
            kick,
            n
        )

        CUDA.synchronize()

        # ----------------------------------------------------
        # Recompute energies.
        # ----------------------------------------------------

        @cuda threads=256 blocks=POPULATION shmem=256*sizeof(Float64) energy_kernel!(
            X_gpu,
            H_gpu,
            E_gpu,
            n
        )

        CUDA.synchronize()

        # ----------------------------------------------------
        # Symmetry diversification after stagnation.
        # ----------------------------------------------------

        if epoch - last_improvement >= 50

            for b in 1:POPULATION

                choices_host[b] =
                    Int32(
                        rand(
                            rng,
                            1:ng
                        )
                    )
            end

            copyto!(
                choices_gpu,
                choices_host
            )

            @cuda threads=256 blocks=POPULATION symmetry_jump_kernel!(
                X_gpu,
                Xtmp_gpu,
                perm_table_gpu,
                choices_gpu,
                n
            )

            CUDA.synchronize()

            X_gpu,
            Xtmp_gpu =
                Xtmp_gpu,
                X_gpu

            H_gpu .=
                B_gpu * X_gpu

            @cuda threads=256 blocks=POPULATION shmem=256*sizeof(Float64) energy_kernel!(
                X_gpu,
                H_gpu,
                E_gpu,
                n
            )

            CUDA.synchronize()

            println(
                "Epoch ",
                epoch,
                " | symmetry diversification"
            )

            last_improvement =
                epoch
        end

        # ----------------------------------------------------
        # Reporting.
        # ----------------------------------------------------

        if epoch == 1 ||
           epoch % 25 == 0

            population_best =
                Int64(
                    round(
                        minimum(E_host)
                    )
                )

            println(
                "Epoch ",
                lpad(epoch,5),
                " | population=",
                population_best,
                " | global=",
                global_best_cost,
                " | gap=",
                global_best_cost -
                KNOWN_BEST,
                " | kick=",
                kick,
                " | moves=",
                total_moves
            )
        end
    end

    # ========================================================
    # FINAL VALIDATION
    # ========================================================

    println()
    println("============================================================")
    println("FINAL VALIDATION")
    println("============================================================")

    final_x =
        Int8.(
            round.(
                global_best_x
            )
        )

    final_cardinality =
        sum(final_x)

    final_cost =
        objective_int(
            B_int,
            final_x
        )

    println(
        "Cardinality = ",
        final_cardinality
    )

    println(
        "Exact cost  = ",
        final_cost
    )

    println(
        "Known best  = ",
        KNOWN_BEST
    )

    println(
        "Difference  = ",
        final_cost -
        KNOWN_BEST
    )

    println(
        "Improvement  = ",
        KNOWN_BEST -
        final_cost
    )

    if final_cardinality !=
       K_TARGET

        error(
            "FINAL CARDINALITY FAILURE"
        )
    end

    # ========================================================
    # WRITE SOLUTION
    # ========================================================

    open(
        "solution_bqop_vns.txt",
        "w"
    ) do io

        println(
            io,
            n
        )

        println(
            io,
            final_cost
        )

        for i in 1:n

            print(
                io,
                final_x[i]
            )

            if i < n
                print(io,' ')
            end
        end

        println(io)
    end

    # ========================================================
    # WRITE STATISTICS
    # ========================================================

    open(
        "tai256c_vns_stats.txt",
        "w"
    ) do io

        println(
            io,
            "instance=tai256c"
        )

        println(
            io,
            "n=$n"
        )

        println(
            io,
            "k=$K_TARGET"
        )

        println(
            io,
            "population=$POPULATION"
        )

        println(
            io,
            "epochs=$MAX_EPOCHS"
        )

        println(
            io,
            "automorphism_group_order=",
            length(group)
        )

        println(
            io,
            "automorphism_orbits=",
            length(orbits)
        )

        println(
            io,
            "reference_best=",
            KNOWN_BEST
        )

        println(
            io,
            "best_cost=",
            final_cost
        )

        println(
            io,
            "improvement=",
            KNOWN_BEST -
            final_cost
        )

        println(
            io,
            "total_1swap_moves=",
            total_moves
        )
    end

    println()
    println("Files:")
    println("  solution_bqop_vns.txt")
    println("  tai256c_automorphisms.txt")
    println("  tai256c_vns_stats.txt")

    println()
    println("============================================================")
    println("SOLVER COMPLETE")
    println("============================================================")

    CUDA.reclaim()

    return final_cost,final_x,group
end


# ============================================================
# RUN
# ============================================================

run_tai256c()
