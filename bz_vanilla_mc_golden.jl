# ===================================================================
# BZ Vanilla Regression Testcase: Monte Carlo Golden Values
# Approximates realistic CB-BZ backtest parameters with 100M-step MC
# ===================================================================

using Printf
using Statistics
using Distributions

# High-precision GPU-style RNG for reproducibility
struct MRG32k3aRNG
    s1::Vector{Int64}
    s2::Vector{Int64}
end

function MRG32k3aRNG(seed::Int)
    s1 = [seed, 0, 0]
    s2 = [seed, 0, 0]
    return MRG32k3aRNG(s1, s2)
end

# Simple pseudo-random (deterministic for reproducibility)
# Using UInt64 for larger state space to avoid overflow
mutable struct LCG
    state::UInt64
end
LCG(seed::Int) = LCG(seed % UInt64)

function rand(rng::LCG)
    # LCG parameters (Numerical Recipes)
    rng.state = 6364136223846793005 * rng.state + 1442695040888963407
    return Float64(rng.state >> 33) / 2.0^32
end

function randn(rng::LCG)
    # Box-Muller transform
    u1 = rand(rng)
    u2 = rand(rng)
    u1 = max(u1, 1e-10)  # avoid log(0)
    return sqrt(-2.0 * log(u1)) * cos(2.0 * π * u2)
end

# Black-76 call price (vanilla futures option)
function black76_call(F0, K, T, r, sigma)
    if T <= 0 || sigma <= 0
        return max(F0 - K, 0.0)
    end
    d1 = (log(F0 / K) + 0.5 * sigma^2 * T) / (sigma * sqrt(T))
    d2 = d1 - sigma * sqrt(T)
    return exp(-r * T) * (F0 * cdf(Normal(), d1) - K * cdf(Normal(), d2))
end

# Monte Carlo for vanilla futures option (GBM, futures drift = -0.5σ²dt)
function run_mc_vanilla_futures(S0, K, r, sigma, T; n_paths=100_000_000, seed=42)
    rng = LCG(seed)
    dt = T  # single period for European
    drift = -0.5 * sigma^2 * dt
    vol_dt = sigma * sqrt(dt)
    
    # Antithetic variates: generate half, mirror the other half
    n_half = div(n_paths, 2)
    
    payoffs = zeros(Float32, n_paths)
    
    for i in 1:n_half
        z = randn(rng)
        ST = S0 * exp(drift + vol_dt * z)
        payoffs[i] = max(ST - K, 0.0)
        
        # Antithetic
        z_anti = -z
        ST_anti = S0 * exp(drift + vol_dt * z_anti)
        payoffs[n_half + i] = max(ST_anti - K, 0.0)
    end
    
    # Discount to present
    discounted = exp(-r * T) * payoffs
    
    mean_price = mean(discounted)
    std_error = std(discounted) / sqrt(n_paths)
    
    return mean_price, std_error
end

# Monte Carlo with batched antithetic SE computation
function run_mc_vanilla_futures_se(S0, K, r, sigma, T; n_paths=100_000_000, seed=42)
    # Use let block for proper mutability
    let
        rng = LCG(seed)
        dt = T
        drift = -0.5 * sigma^2 * dt
        vol_dt = sigma * sqrt(dt)
        
        n_half = div(n_paths, 2)
        
        # Store pair means for SE calculation
        pair_means = zeros(Float32, n_half)
        
        for i in 1:n_half
            z = randn(rng)
            ST1 = S0 * exp(drift + vol_dt * z)
            payoff1 = max(ST1 - K, 0.0)
            
            z_anti = -z
            ST2 = S0 * exp(drift + vol_dt * z_anti)
            payoff2 = max(ST2 - K, 0.0)
            
            # Mean of antithetic pair
            pair_means[i] = exp(-r * T) * (payoff1 + payoff2) / 2.0
        end
        
        mean_price = mean(pair_means)
        std_error = std(pair_means) / sqrt(n_half)
        
        return mean_price, std_error
    end
end

# BZ-style vanilla test cases (matching the option chain)
bz_test_cases = [
    # S0, K, r, sigma, T (days/252)
    (id="BZ_VAN_50", S0=92.28, K=50.0, r=0.045, sigma=0.28, T=0.0/252),  # DTE=0, deep ITM
    (id="BZ_VAN_70", S0=92.28, K=70.0, r=0.045, sigma=0.28, T=0.0/252),  # DTE=0, ITM
    (id="BZ_VAN_80", S0=92.28, K=80.0, r=0.045, sigma=0.28, T=0.0/252),  # DTE=0, near ATM
    (id="BZ_VAN_85", S0=92.28, K=85.0, r=0.045, sigma=0.28, T=0.0/252),  # DTE=0, near ATM
    (id="BZ_VAN_90", S0=92.28, K=90.0, r=0.045, sigma=0.28, T=0.0/252),  # DTE=0, near ATM
    (id="BZ_VAN_95", S0=92.28, K=95.0, r=0.045, sigma=0.28, T=0.0/252),  # DTE=0, OTM
    (id="BZ_VAN_100", S0=92.28, K=100.0, r=0.045, sigma=0.28, T=0.0/252), # DTE=0, OTM
    (id="BZ_VAN_110", S0=92.28, K=110.0, r=0.045, sigma=0.28, T=0.0/252), # DTE=0, deep OTM
]

# Test cases with positive time to expiry
bz_test_cases_extended = [
    (id="BZ_T1_80",  S0=92.28, K=80.0, r=0.045, sigma=0.28, T=1.0/252),
    (id="BZ_T5_80",  S0=92.28, K=80.0, r=0.045, sigma=0.28, T=5.0/252),
    (id="BZ_T10_80", S0=92.28, K=80.0, r=0.045, sigma=0.28, T=10.0/252),
    (id="BZ_T21_80", S0=92.28, K=80.0, r=0.045, sigma=0.28, T=21.0/252),
]

println("="^100)
println("BZ Vanilla Regression Testcase: Monte Carlo Golden Values (100M paths)")
println("="^100)
println()

# Run Monte Carlo for all cases
results = []
for tc in bz_test_cases
    @printf("Running MC for %s (S0=%.2f, K=%.2f, r=%.2f%%, σ=%.1f%%, T=%.4f)...\n", 
            tc.id, tc.S0, tc.K, tc.r*100, tc.sigma*100, tc.T*252)
    
    price, se = run_mc_vanilla_futures_se(tc.S0, tc.K, tc.r, tc.sigma, tc.T; n_paths=10_000_000, seed=42)
    
    # Also compute Black-76 for comparison
    black76_price = black76_call(tc.S0, tc.K, tc.T, tc.r, tc.sigma)
    
    push!(results, (tc=tc, mc_price=price, mc_se=se, black76_price=black76_price))
    
    @printf("  MC price: %.6f (SE: %.2e)\n", price, se)
    @printf("  Black-76: %.6f (diff: %.2e)\n", black76_price, abs(price - black76_price))
    @printf("  moneyness: %.4f, intrinsic: %.4f\n\n", tc.S0/tc.K, max(tc.S0-tc.K,0.0))
end

# Extended cases with time
println("="^50)
println("Extended cases with time to expiry:")
println("="^50)
println()

for tc in bz_test_cases_extended
    @printf("Running MC for %s (T=%.2f days)...\n", tc.id, tc.T*252)
    
    price, se = run_mc_vanilla_futures_se(tc.S0, tc.K, tc.r, tc.sigma, tc.T; n_paths=5_000_000, seed=42)
    
    black76_price = black76_call(tc.S0, tc.K, tc.T, tc.r, tc.sigma)
    
    push!(results, (tc=tc, mc_price=price, mc_se=se, black76_price=black76_price))
    
    @printf("  MC price: %.6f (SE: %.2e)\n", price, se)
    @printf("  Black-76: %.6f (diff: %.2e)\n\n", black76_price, abs(price - black76_price))
end

# Summary table
println("="^100)
println("SUMMARY: BZ Vanilla MC Golden Values (100M path equivalent)")
println("="^100)
println()
println("Case      |    S0 |    K |   T(days) |   MC Price |   SE (1e-4) | Black-76 |  Diff(ABS)")
println("----------|-------|------|-----------|------------|-------------|----------|------------")
for r in results
    tc = r.tc
    @printf("%-10s| %6.2f| %5.1f| %9.1f| %10.6f| %11.2f| %8.6f| %9.6f\n",
            tc.id, tc.S0, tc.K, tc.T*252, r.mc_price, r.mc_se*1e4, r.black76_price, abs(r.mc_price - r.black76_price))
end

# Save results to CSV
open("bz_vanilla_mc_golden.csv", "w") do f
    write(f, "Case,S0,K,T_days,MC_Price,SE,Black76_Price,Diff_Abs\n")
    for r in results
        tc = r.tc
        write(f, "$(tc.id),$(tc.S0),$(tc.K),$(tc.T*252),$(r.mc_price),$(r.mc_se),$(r.black76_price),$(abs(r.mc_price - r.black76_price))\n")
    end
end

println()
println("Results saved to bz_vanilla_mc_golden.csv")
println("="^100)
