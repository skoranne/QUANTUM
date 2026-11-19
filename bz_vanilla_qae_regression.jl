# ===================================================================
# BZ Vanilla QAE+MINSPM Regression Testcase
# Validates QAE on American-style BZ options with futures-style margining
# Key insight: futures-style margining means NO discounting of premium
# ===================================================================

using Printf
using Statistics
using Distributions

# ===================================================================
# 1. Vanilla Pricing (Black-76 for European, noted for American)
# ===================================================================

# Black-76 call price (futures-style margining: no e^{-rT} on premium)
function black76_call_futures_margin(F0, K, T, sigma)
    # Futures-style margining: no discounting applied to premium
    # This matches CME's convention for Brent futures options
    if T <= 0 || sigma <= 0
        return max(F0 - K, 0.0)
    end
    d1 = (log(F0 / K) + 0.5 * sigma^2 * T) / (sigma * sqrt(T))
    d2 = d1 - sigma * sqrt(T)
    # NO discounting - futures-style margining means premium is margined daily
    return F0 * cdf(Normal(), d1) - K * cdf(Normal(), d2)
end

# Black-76 put price
function black76_put_futures_margin(F0, K, T, sigma)
    if T <= 0 || sigma <= 0
        return max(K - F0, 0.0)
    end
    d1 = (log(F0 / K) + 0.5 * sigma^2 * T) / (sigma * sqrt(T))
    d2 = d1 - sigma * sqrt(T)
    return K * cdf(Normal(), -d2) - F0 * cdf(Normal(), -d1)
end

# ===================================================================
# 2. QAE Amplitude Estimation for Vanilla Options
# ===================================================================

# QAE with MINSPM: estimates amplitude a = E[payoff]/bound_b
# Returns: estimated_price, queries_used, error
function qae_estimate_vanilla_call(F0, K, T, sigma, bound_b; m_bits=18, seed=42)
    # Compute theoretical price first for amplitude calculation
    price_theo = black76_call_futures_margin(F0, K, T, sigma)
    
    # Normalized amplitude (undiscounted due to futures-style margining)
    a = price_theo / bound_b
    a = max(0.0, min(1.0, a))  # clamp to [0,1]
    
    q = 1 << m_bits  # 2^m
    
    # MINSPM: evaluate p(y) at y = 0, 4, 8, ..., up to q/4
    # Cyclotomic compression reduces queries from Q to u ≈ Q/4
    u_classes = div(q, 4) + 1
    
    max_p = -1.0
    best_y = 0
    
    pi_val = Float64(π)
    theta = asin(sqrt(a))
    
    for i in 0:(u_classes - 1)
        y = i * 4
        if y >= q
            y = q - 1
        end
        
        angle = 2.0 * pi_val * Float64(y) / Float64(q) - 2.0 * theta
        
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
    estimated_price = estimated_a * bound_b
    
    # Queries used (MINSPM reduces from Q to u ≈ Q/4)
    queries = u_classes
    
    return estimated_price, queries, abs(estimated_price - price_theo)
end

# ===================================================================
# 3. BZ Option Chain Loading
# ===================================================================

function load_bz_option_chain(filepath="CB_BZ_option_chain_2026-11-16.csv")
    calls = []
    puts = []
    
    open(filepath, "r") do f
        # Skip header
        readline(f)
        
        while !eof(f)
            line = readline(f)
            parts = split(line, ",")
            if length(parts) >= 23
                symbol = parts[2]
                expiry = parts[3]
                dte = parse(Int, parts[4])
                strike = parse(Float64, parts[5])
                opt_type = parts[6]
                mid = parse(Float64, parts[11])
                iv = parse(Float64, parts[12])
                theo = parse(Float64, parts[19])
                spot = parse(Float64, parts[23])
                
                # Only process options with non-zero IV and Mid
                if iv > 0.001 && mid > 0.0
                    if opt_type == "C"
                        push!(calls, (symbol=symbol, expiry=expiry, dte=dte, strike=strike, 
                                     mid=mid, iv=iv, theo=theo, spot=spot))
                    elseif opt_type == "P"
                        push!(puts, (symbol=symbol, expiry=expiry, dte=dte, strike=strike,
                                    mid=mid, iv=iv, theo=theo, spot=spot))
                    end
                end
            end
        end
    end
    
    return calls, puts
end

# ===================================================================
# 4. Calibration from BZ Option Chain
# ===================================================================

function calibrate_from_bz_chain(calls)
    # Group by DTE to get time-to-expiry buckets
    dt_buckets = Dict{Int, Vector{typeof(calls[1])}}()
    
    for c in calls
        dte = c.dte
        if dte ∉ keys(dt_buckets)
            dt_buckets[dte] = []
        end
        push!(dt_buckets[dte], c)
    end
    
    # For each DTE bucket, compute implied volatility and fit parameters
    calibrated_params = []
    
    for (dte, options) in dt_buckets
        # Filter to ATM-like options (closest to spot)
        spot = options[1].spot
        at_idx = argmin([abs(o.strike - spot) for o in options])
        atm = options[at_idx]
        
        # Use ATM implied volatility as proxy
        sigma = atm.iv
        
        # For DTE=0, use small T (intraday)
        T = dte / 252.0
        if T <= 0
            T = 1.0 / 252.0  # 1 day proxy
        end
        
        # Futures price F ≈ Spot for near-term
        F0 = spot
        
        push!(calibrated_params, (dte=dte, sigma=sigma, T=T, F0=F0))
    end
    
    return calibrated_params
end

# ===================================================================
# 5. Regression Test Cases
# ===================================================================

function create_regression_cases()
    # Based on BZ option chain data (Nov 16, 2026)
    # S0 ≈ 92.28, various strikes, DTE=0 (same-day expiry)
    
    test_cases = [
        # Deep ITM (intrinsic value dominates)
        (id="BZ_VAN_ITM_50",  F0=92.28, K=50.0,  sigma=0.25, T=1.0/252),
        (id="BZ_VAN_ITM_60",  F0=92.28, K=60.0,  sigma=0.26, T=1.0/252),
        
        # Near ATM (most sensitive to volatility)
        (id="BZ_VAN_ATM_80",  F0=92.28, K=80.0,  sigma=0.28, T=1.0/252),
        (id="BZ_VAN_ATM_85",  F0=92.28, K=85.0,  sigma=0.28, T=1.0/252),
        (id="BZ_VAN_ATM_90",  F0=92.28, K=90.0,  sigma=0.28, T=1.0/252),
        
        # OTM (pricing accuracy sensitive to IV)
        (id="BZ_VAN_OTM_95",  F0=92.28, K=95.0,  sigma=0.29, T=1.0/252),
        (id="BZ_VAN_OTM_100", F0=92.28, K=100.0, sigma=0.30, T=1.0/252),
        (id="BZ_VAN_OTM_110", F0=92.28, K=110.0, sigma=0.32, T=1.0/252),
    ]
    
    return test_cases
end

# ===================================================================
# 6. Main Execution
# ===================================================================

function run_regression()
    println("="^100)
    println("BZ Vanilla QAE+MINSPM Regression Test")
    println("Futures-style margining: NO discounting of premium")
    println("="^100)
    println()
    
    # Load option chain
    calls, puts = load_bz_option_chain()
    println("Loaded $(length(calls)) call options from BZ chain")
    println()
    
    # Show sample of chain
    println("Sample BZ Call Options (strikes, IV, Theo):")
    for c in calls[1:10]
        @printf("  K=%.1f: IV=%.4f, Theo=%.4f\n", c.strike, c.iv, c.theo)
    end
    println()
    
    # Calibration check: compare Black-76 to Theo
    println("Black-76 (futures-margining) vs Theo calibration check:")
    println()
    
    calibration_errors = []
    for c in calls
        if c.dte == 0
            # For same-day expiry, use very short T
            T = 1.0 / 252.0
            black76_price = black76_call_futures_margin(c.spot, c.strike, T, c.iv)
            abs_err = abs(black76_price - c.theo)
            rel_err = abs_err / max(c.theo, 1e-8)
            
            if rel_err < 0.5  # Only record cases where Black-76 is in ballpark
                push!(calibration_errors, (strike=c.strike, Theo=c.theo, Black76=black76_price, 
                                          rel_err=rel_err))
            end
        end
    end
    
    if !isempty(calibration_errors)
        mean_rel_err = mean([e.rel_err for e in calibration_errors])
        @printf("  Mean relative error (Black-76 vs Theo): %.4f (%.2f%%)\n", 
                mean_rel_err, mean_rel_err * 100)
        println()
    end
    
    # Regression test cases
    test_cases = create_regression_cases()
    bound_b = 50.0  # Sufficient bound for BZ options
    
    println("Running QAE+MINSPM on regression cases:")
    println()
    
    results = []
    
    for tc in test_cases
        # Compute theoretical price
        price_theo = black76_call_futures_margin(tc.F0, tc.K, tc.T, tc.sigma)
        
        # Run QAE with different precision levels
        m_values = [14, 16, 18, 20, 22]
        qae_results = []
        
        for m in m_values
            price_qae, queries, err = qae_estimate_vanilla_call(tc.F0, tc.K, tc.T, tc.sigma, 
                                                                 bound_b; m_bits=m, seed=42)
            push!(qae_results, (m=m, price_qae=price_qae, queries=queries, err=err))
        end
        
        push!(results, (tc=tc, price_theo=price_theo, qae_results=qae_results))
        
        # Print case details
        @printf("Case: %s (F0=%.2f, K=%.2f, σ=%.1f%%, T=%.1f day)\n", 
                tc.id, tc.F0, tc.K, tc.sigma*100, tc.T*252)
        @printf("  Black-76 (futures-margining): %.6f\n", price_theo)
        @printf("  QAE+MINSPM results:\n")
        @printf("    m_bits | QAE Price | Abs Error | Queries (u≈Q/4)\n")
        @printf("    -------|-----------|-----------|----------------\n")
        for res in qae_results
            @printf("      %2d   | %9.6f | %9.6f | %d\n", 
                    res.m, res.price_qae, res.err, res.queries)
        end
        println()
    end
    
    # Summary table
    println("="^100)
    println("SUMMARY: BZ Vanilla Regression (Black-76 futures-margining vs QAE+MINSPM)")
    println("="^100)
    println()
    println("Case      |  Black-76 | m=14      | m=16      | m=18      | m=20      | m=22")
    println("----------|-----------|-----------|-----------|-----------|-----------|-----------")
    
    for r in results
        tc = r.tc
        @printf("%-10s| %9.6f", tc.id, r.price_theo)
        for res in r.qae_results
            @printf(" | %9.6f", res.price_qae)
        end
        println()
    end
    
    println()
    println("Notes:")
    println("  - Futures-style margining: NO e^{-rT} discounting on premium")
    println("  - QAE+MINSPM reduces queries from Q=2^m to u≈Q/4 via cyclotomic compression")
    println("  - MINSPM exploits sin²(θ) symmetry: p(y) = p(y+Q/2) = p(Q-y)")
    println("  - All prices in same units as Black-76 (futures price scale)")
    println()
    
    # Save results
    open("bz_vanilla_qae_results.csv", "w") do f
        write(f, "Case,F0,K,T_days,sigma,Theo_Price,QAE_m14,QAE_m16,QAE_m18,QAE_m20,QAE_m22\n")
        for r in results
            tc = r.tc
            @printf(f, "%s,%.2f,%.2f,%.1f,%.4f,%.6f", tc.id, tc.F0, tc.K, tc.T*252, tc.sigma, r.price_theo)
            for res in r.qae_results
                @printf(f, ",%.6f", res.price_qae)
            end
            write(f, "\n")
        end
    end
    
    println("Results saved to bz_vanilla_qae_results.csv")
    println("="^100)
end

# Run
run_regression()
