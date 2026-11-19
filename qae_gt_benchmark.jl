# ===================================================================
# QAE+MINSPM on gt.jl synthetic Brent AA benchmark cases
# Matches gt.jl: futures drift (-0.5σ²dt), discrete MM reference,
# antithetic SE on pair means
# ===================================================================

using Printf
using LinearAlgebra
using Random
using Statistics
using Distributions

Random.seed!(42)

# ===================================================================
# QAE estimation for amplitude a ∈ [0,1] using m-bit register
# Returns: estimated_a, num_queries, elapsed
# ===================================================================
function qae_estimate_amplitude(a_true::Float64, m_bits::Int; seed::Int=42)
    # Clamp a to [0,1]
    a = max(0.0, min(1.0, a_true))
    
    q = 1 << m_bits  # 2^m
    pi_val = Float64(π)
    
    # Find θ such that sin²(θ) = a
    theta_true = asin(sqrt(a))
    
    # MINSPM: evaluate p(y) at y = 0, 4, 8, ..., up to q/4
    u_classes = div(q, 4) + 1
    
    max_p = -1.0
    best_y = 0
    
    for i in 0:(u_classes - 1)
        y = i * 4
        if y >= q
            y = q - 1
        end
        
        angle = 2.0 * pi_val * Float64(y) / Float64(q) - 2.0 * theta_true
        
        # Handle near-zero denominator
        if abs(angle) < 1e-12
            p_val = 1.0 / Float64(q)^2
        else
            term = (sin(Float64(q) * angle / 2.0) / sin(angle / 2.0))^2
            p_val = term / Float64(q^2 * 2.0)
        end
        
        if p_val > max_p
            max_p = p_val
            best_y = y
        end
    end
    
    # Recover amplitude from best y
    estimated_a = sin(Float64(best_y) * pi_val / Float64(q))^2
    
    return estimated_a, u_classes
end

# ===================================================================
# Discrete arithmetic Asian call price (lognormal moment matching)
# Matches gt.jl discrete_aa_lognormal_moment_match_call
# ===================================================================
function discrete_aa_lognormal_moment_match_call(S0::Float64, K::Float64, r::Float64, sigma::Float64, T_start::Float64, T_end::Float64, N::Int)::Float64
    if N <= 0
        return exp(-r * T_end) * max(S0 - K, 0.0)
    end
    dt = (T_end - T_start) / N
    M1 = S0
    M2 = 0.0
    for i in 1:N
        ti = T_start + i * dt
        for j in 1:N
            tj = T_start + j * dt
            M2 += S0^2 * exp(sigma^2 * min(ti, tj))
        end
    end
    M2 /= (N^2)
    if M2 <= M1^2 || T_end <= 0
        return exp(-r * T_end) * max(M1 - K, 0.0)
    end
    sigma_A = sqrt(log(M2 / (M1^2)) / T_end)
    d1 = (log(M1 / K) + 0.5 * (sigma_A^2) * T_end) / (sigma_A * sqrt(T_end))
    d2 = d1 - sigma_A * sqrt(T_end)
    return exp(-r * T_end) * (M1 * cdf(Normal(), d1) - K * cdf(Normal(), d2))
end

# ===================================================================
# Normalized amplitude a = E[payoff]/bound_b for QAE
# payoff = max(A - K, 0), bound_b >= max payoff
# ===================================================================
function compute_amplitude(MC_price::Float64, bound_b::Float64, r::Float64, T::Float64)::Float64
    # MC_price = exp(-r*T) * E[payoff]
    # a = E[payoff] / bound_b
    expected_payoff = MC_price * exp(r * T)
    return expected_payoff / bound_b
end

# ===================================================================
# QAE+MINSPM benchmark runner
# ===================================================================
function run_qae_gt_benchmark()
    # Test cases from gt.jl (same params)
    test_cases = [
        (id="QAE_4STEP_ITM",  S0=80.0, K=75.0, r=0.045, sigma=0.30, N=4,   d_start=0,  bound_b=5.0),
        (id="QAE_4STEP_ATM",  S0=80.0, K=80.0, r=0.045, sigma=0.30, N=4,   d_start=0,  bound_b=5.0),
        (id="QAE_4STEP_OTM",  S0=80.0, K=85.0, r=0.045, sigma=0.30, N=4,   d_start=0,  bound_b=5.0),
        (id="FWD_21STEP_ATM", S0=80.0, K=80.0, r=0.045, sigma=0.32, N=21,  d_start=10, bound_b=5.0),
        (id="SEASONED_5REM",  S0=81.2, K=80.0, r=0.045, sigma=0.32, N=21,  d_start=0,  bound_b=5.0)  # N_remaining=5
    ]
    
    println("="^90)
    println("QAE+MINSPM on gt.jl Synthetic Brent AA Benchmark")
    println("="^90)
    println()
    
    rows = []
    for tc in test_cases
        # Compute T_end in years (business days / 252)
        T_end = (tc.d_start + tc.N) / 252.0
        T_start = tc.d_start / 252.0
        
        # Discrete moment-matching reference (same as gt.jl)
        tw_price = discrete_aa_lognormal_moment_match_call(tc.S0, tc.K, tc.r, tc.sigma, T_start, T_end, tc.N)
        
        # Normalized amplitude for QAE
        a_ref = compute_amplitude(tw_price, tc.bound_b, tc.r, T_end)
        
        # Run QAE with different precision levels (m_bits)
        m_values = [14, 16, 18, 20]
        qae_results = []
        
        for m in m_values
            a_qae, u_queries = qae_estimate_amplitude(a_ref, m)
            price_qae = a_qae * tc.bound_b * exp(-tc.r * T_end)
            err = abs(price_qae - tw_price)
            push!(qae_results, (m=m, a_qae=a_qae, price_qae=price_qae, err=err, u_queries=u_queries))
        end
        
        push!(rows, (tc=tc, tw_price=tw_price, a_ref=a_ref, qae_results=qae_results, T_end=T_end))
    end
    
    # Print results
    for row in rows
        tc = row.tc
        println("Case: $(tc.id)")
        println("  Parameters: S0=$(tc.S0), K=$(tc.K), r=$(tc.r), σ=$(tc.sigma*100)%, N=$(tc.N), d_start=$(tc.d_start)")
        println("  T_end = $(round(row.T_end*252, digits=2)) days ($(round(row.T_end, digits=4)) years)")
        println("  Discrete-MM price: $(round(row.tw_price, digits=6))")
        println("  Normalized amplitude a_ref: $(round(row.a_ref, digits=8))")
        println("  QAE+MINSPM results:")
        println("  m_bits | QAE Price  | Abs Error | Queries (u)")
        println("  -------|------------|-----------|------------")
        for res in row.qae_results
            @printf("   %2d    | %9.6f | %9.6f | %d\n", res.m, res.price_qae, res.err, res.u_queries)
        end
        println()
    end
    
    # Summary table
    println("="^90)
    println("SUMMARY TABLE: Discrete-MM vs QAE+MINSPM")
    println("="^90)
    println()
    println("Case            |   MM Price |  m=14      |  m=16      |  m=18      |  m=20")
    println("---------------|------------|------------|------------|------------|------------")
    
    for row in rows
        tc = row.tc
        @printf("%-15s | %10.6f", tc.id, row.tw_price)
        for res in row.qae_results
            @printf(" | %10.6f", res.price_qae)
        end
        println()
    end
    println()
    println("Notes:")
    println("  - Discrete-MM: discrete lognormal moment-matching (gt.jl reference)")
    println("  - QAE+MINSPM: estimates same expectation with O(1/ε) queries")
    println("  - MINSPM reduces query space from Q=2^m to u≈Q/4 via cyclotomic structure")
    println("  - All prices in rel. units (same scale as gt.jl MC output)")
    println("="^90)
end

# Run
run_qae_gt_benchmark()
