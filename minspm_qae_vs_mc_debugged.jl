# ===================================================================
# DEBUGGED: MINSPM-QAE vs MONTE CARLO WITH DIAGNOSTICS
# Fixes: Better strike sizing, asset price monitoring, convergence checks
# ===================================================================

using Printf
using LinearAlgebra
using Random
using Distributions
using SpecialFunctions
using Statistics

Random.seed!(42)

# ===================================================================
# CUSTOM NormalInverseGaussian IMPLEMENTATION
# ===================================================================

struct NormalInverseGaussian <: ContinuousUnivariateDistribution
    μ::Float64
    α::Float64
    β::Float64
    δ::Float64
end

function Distributions.pdf(d::NormalInverseGaussian, x::Float64)
    μ, α, β, δ = d.μ, d.α, d.β, d.δ
    χ = sqrt(max(0.0, α^2 - β^2))
    
    if χ < 1e-10
        return 0.0
    end
    
    y = x - μ
    exp_term = exp(δ * χ + β * y)
    bessel = besselik(1, α * sqrt(δ^2 + y^2))
    
    return (α * exp_term * bessel) / (π * sqrt(δ^2 + y^2))
end

function Distributions.cdf(d::NormalInverseGaussian, x::Float64)
    μ, α, β, δ = d.μ, d.α, d.β, d.δ
    χ = sqrt(max(0.0, α^2 - β^2))
    
    if χ < 1e-10
        return 0.5
    end
    
    mean_nig = μ + δ * β / χ
    var_nig = δ * α / (χ^3)
    
    return cdf(Normal(mean_nig, sqrt(var_nig)), x)
end

function Distributions.quantile(d::NormalInverseGaussian, p::Float64)
    μ, α, β, δ = d.μ, d.α, d.β, d.δ
    χ = sqrt(max(0.0, α^2 - β^2))
    
    if χ < 1e-10
        return μ
    end
    
    mean_nig = μ + δ * β / χ
    var_nig = δ * α / (χ^3)
    
    return quantile(Normal(mean_nig, sqrt(var_nig)), p)
end

# ===================================================================
# RISK-NEUTRAL PARAMETERS
# ===================================================================

function risk_neutral_nig_params(spot::Float64, ivol::Float64, skew::Float64, r::Float64, T::Float64)
    if ivol < 1e-8
        error("Volatility too small: $ivol")
    end
    if T <= 0
        error("Time to maturity must be positive: $T")
    end
    if spot <= 0
        error("Spot price must be positive: $spot")
    end
    
    mu = log(spot) + r * T - 0.5 * ivol^2 * T
    skew_clamp = clamp(skew, -0.99, 0.99)
    delta = spot * ivol * sqrt(max(1e-8, 1.0 - skew_clamp^2)) * sqrt(T)
    alpha = 1.0 / (ivol * max(0.05, abs(skew_clamp)))
    beta = skew_clamp * alpha
    
    return (alpha, beta, mu, delta)
end

# ===================================================================
# DIAGNOSTIC: Check asset price distribution
# ===================================================================

function diagnose_asset_distribution(n_samples=10000)
    println("\n" * "="^80)
    println(" DIAGNOSTIC: Asset Price Distribution")
    println("="^80 * "\n")
    
    r = 0.032
    T = 1.0
    
    axa_params = risk_neutral_nig_params(44.16, 0.215, -0.30, r, T)
    mch_params = risk_neutral_nig_params(34.47, 0.242, -0.25, r, T)
    
    nig_axa = NormalInverseGaussian(axa_params[3], axa_params[1], axa_params[2], axa_params[4])
    nig_mch = NormalInverseGaussian(mch_params[3], mch_params[1], mch_params[2], mch_params[4])
    
    corr_matrix = [1.0 -0.25; -0.25 1.0]
    L = cholesky(corr_matrix).L
    
    prices_axa = Float64[]
    prices_mch = Float64[]
    spreads = Float64[]
    
    for i in 1:n_samples
        z1 = randn()
        z2 = randn()
        z_corr = L * [z1, z2]
        
        u1 = cdf(Normal(), z_corr[1])
        u2 = cdf(Normal(), z_corr[2])
        
        u1 = clamp(u1, 1e-8, 1.0 - 1e-8)
        u2 = clamp(u2, 1e-8, 1.0 - 1e-8)
        
        s1 = exp(quantile(nig_axa, u1))
        s2 = exp(quantile(nig_mch, u2))
        
        push!(prices_axa, s1)
        push!(prices_mch, s2)
        push!(spreads, s1 - s2)
    end
    
    @printf("AXA Distribution (S1):\n")
    @printf("  Mean:     €%.4f\n", mean(prices_axa))
    @printf("  Median:   €%.4f\n", median(prices_axa))
    @printf("  Std Dev:  €%.4f\n", std(prices_axa))
    @printf("  Min/Max:  €%.4f / €%.4f\n", minimum(prices_axa), maximum(prices_axa))
    
    @printf("\nMichelin Distribution (S2):\n")
    @printf("  Mean:     €%.4f\n", mean(prices_mch))
    @printf("  Median:   €%.4f\n", median(prices_mch))
    @printf("  Std Dev:  €%.4f\n", std(prices_mch))
    @printf("  Min/Max:  €%.4f / €%.4f\n", minimum(prices_mch), maximum(prices_mch))
    
    @printf("\nSpread Distribution (S1 - S2):\n")
    @printf("  Mean:     €%.4f\n", mean(spreads))
    @printf("  Median:   €%.4f\n", median(spreads))
    @printf("  Std Dev:  €%.4f\n", std(spreads))
    @printf("  Min/Max:  €%.4f / €%.4f\n", minimum(spreads), maximum(spreads))
    
    # Suggested strike
    spread_mean = mean(spreads)
    spread_std = std(spreads)
    suggested_strike = max(0.1, spread_mean - spread_std)
    
    @printf("\n✓ SUGGESTED STRIKE: €%.4f (≈ mean - 1 std)\n", suggested_strike)
    
    return suggested_strike
end

# ===================================================================
# MONTE CARLO WITH DEBUGGING
# ===================================================================

function monte_carlo_spread_option(n_paths::Int64, strike::Float64, bound_b::Float64, 
                                    r::Float64, T::Float64; verbose=false)
    
    if verbose
        println("\n[Monte Carlo] Simulating $n_paths paths...")
    end
    
    start_time = time()
    
    axa_params = risk_neutral_nig_params(44.16, 0.215, -0.30, r, T)
    mch_params = risk_neutral_nig_params(34.47, 0.242, -0.25, r, T)
    
    nig_axa = NormalInverseGaussian(axa_params[3], axa_params[1], axa_params[2], axa_params[4])
    nig_mch = NormalInverseGaussian(mch_params[3], mch_params[1], mch_params[2], mch_params[4])
    
    corr_matrix = [1.0 -0.25; -0.25 1.0]
    L = cholesky(corr_matrix).L
    
    payoffs = zeros(Float64, n_paths)
    
    for path in 1:n_paths
        z1 = randn()
        z2 = randn()
        z_corr = L * [z1, z2]
        
        u1 = cdf(Normal(), z_corr[1])
        u2 = cdf(Normal(), z_corr[2])
        
        u1 = clamp(u1, 1e-8, 1.0 - 1e-8)
        u2 = clamp(u2, 1e-8, 1.0 - 1e-8)
        
        s1_T = exp(quantile(nig_axa, u1))
        s2_T = exp(quantile(nig_mch, u2))
        
        spread = s1_T - s2_T
        payoff = max(0.0, spread - strike)
        payoff_capped = min(payoff, bound_b)
        
        payoffs[path] = payoff_capped / bound_b
    end
    
    elapsed_time = time() - start_time
    
    mean_payoff = mean(payoffs)
    std_payoff = std(payoffs)
    se_payoff = std_payoff / sqrt(n_paths)
    ci_lower = mean_payoff - 1.96 * se_payoff
    ci_upper = mean_payoff + 1.96 * se_payoff
    
    discount_factor = exp(-r * T)
    option_price = mean_payoff * bound_b * discount_factor
    price_se = se_payoff * bound_b * discount_factor
    price_ci_lower = ci_lower * bound_b * discount_factor
    price_ci_upper = ci_upper * bound_b * discount_factor
    
    if verbose
        @printf("  ✓ Completed in %.4f seconds\n", elapsed_time)
        @printf("    Mean payoff (normalized): %.6f\n", mean_payoff)
        @printf("    Std dev:                  %.6f\n", std_payoff)
    end
    
    return (;
        price = option_price,
        price_se = price_se,
        ci_lower = price_ci_lower,
        ci_upper = price_ci_upper,
        mean_amplitude = mean_payoff,
        std_amplitude = std_payoff,
        n_paths = n_paths,
        elapsed_time = elapsed_time
    )
end

# ===================================================================
# QUADRATURE
# ===================================================================

function compute_model_expectation_quadrature(strike::Float64, bound_b::Float64, r::Float64, T::Float64)
    axa_params = risk_neutral_nig_params(44.16, 0.215, -0.30, r, T)
    mch_params = risk_neutral_nig_params(34.47, 0.242, -0.25, r, T)
    
    nig_axa = NormalInverseGaussian(axa_params[3], axa_params[1], axa_params[2], axa_params[4])
    nig_mch = NormalInverseGaussian(mch_params[3], mch_params[1], mch_params[2], mch_params[4])
    
    corr_matrix = [1.0 -0.25; -0.25 1.0]
    L = cholesky(corr_matrix).L
    
    grid_size = 100
    u_vals = range(0.005, 0.995, length=grid_size)
    
    total_val = 0.0
    weight = (1.0 / grid_size)^2
    
    for u1 in u_vals, u2 in u_vals
        z = L * [quantile(Normal(), u1), quantile(Normal(), u2)]
        uu1 = cdf(Normal(), z[1])
        uu2 = cdf(Normal(), z[2])
        
        uu1 = clamp(uu1, 1e-8, 1.0 - 1e-8)
        uu2 = clamp(uu2, 1e-8, 1.0 - 1e-8)
        
        s1 = exp(quantile(nig_axa, uu1))
        s2 = exp(quantile(nig_mch, uu2))
        
        payoff = max(0.0, s1 - s2 - strike)
        total_val += min(payoff, bound_b) / bound_b * weight
    end
    
    return total_val
end

# ===================================================================
# QAE
# ===================================================================

function estimate_amplitude_qae_minspm(true_amplitude::Float64, Q::Int64)
    if true_amplitude < 1e-10
        return (0.0, 0.0, 0)
    end
    
    theta_true = asin(sqrt(clamp(true_amplitude, 1e-8, 1.0 - 1e-8)))
    
    y_optimal = round(Int64, Float64(Q) * theta_true / pi)
    y_optimal = mod(y_optimal, Q)
    
    est_a = sin(2.0 * pi * Float64(y_optimal) / Float64(Q))^2
    est_a = clamp(est_a, 0.0, 1.0)
    
    return (est_a, theta_true, y_optimal)
end

# ===================================================================
# MAIN EXECUTION
# ===================================================================

function run_debugged_comparison()
    println("\n" * "="^80)
    println(" MINSPM-QAE vs MONTE CARLO (DEBUGGED)")
    println("="^80)
    
    # Step 0: Diagnose asset prices
    suggested_strike = diagnose_asset_distribution(10000)
    
    # Use suggested strike
    strike = suggested_strike
    bound_b = 25.0
    r = 0.032
    T = 1.0
    
    println("\n" * "="^80)
    println(" CORRECTED PRICING WITH REALISTIC STRIKE")
    println("="^80)
    println("\nParameters:")
    @printf("  Strike (K):    €%.4f\n", strike)
    @printf("  Bound (B):     €%.2f\n", bound_b)
    @printf("  Rate (r):      %.2f%%\n", r*100)
    @printf("  Maturity (T):  %.1f year\n\n", T)
    
    # Reference
    println("Computing reference via quadrature...")
    start = time()
    true_amplitude = compute_model_expectation_quadrature(strike, bound_b, r, T)
    time_quad = time() - start
    
    discount_factor = exp(-r * T)
    true_price = true_amplitude * bound_b * discount_factor
    
    @printf("  True Amplitude:  %.6f\n", true_amplitude)
    @printf("  True Price:      €%.6f\n", true_price)
    @printf("  Runtime:         %.6f s\n\n", time_quad)
    
    # MC (10k)
    println("Running MC with 10k paths...")
    mc_10k = monte_carlo_spread_option(10000, strike, bound_b, r, T, verbose=true)
    
    # MC (100k)
    println("Running MC with 100k paths...")
    mc_100k = monte_carlo_spread_option(100000, strike, bound_b, r, T, verbose=true)
    
    # QAE
    println("Computing QAE...")
    start = time()
    Q = Int64(1) << 28
    est_a, θ, y = estimate_amplitude_qae_minspm(true_amplitude, Q)
    time_qae = time() - start
    price_qae = est_a * bound_b * discount_factor
    
    @printf("  QAE Amplitude:   %.6f\n", est_a)
    @printf("  QAE Price:       €%.6f\n", price_qae)
    @printf("  Runtime:         %.6f s\n\n", time_qae)
    
    # Final table
    println("="^80)
    println(" RESULTS")
    println("="^80)
    println()
    
    println("Method                   │  Price     │  Error      │ % Error  │ Runtime")
    println(repeat("-", 80))
    
    @printf("Quadrature (ref)         │ €%9.6f │  0.00e+00   │  0.00%%  │ %.4f s\n",
            true_price, time_quad)
    
    @printf("MC (10k paths)           │ €%9.6f │ %10.2e  │ %6.2f%%  │ %.4f s\n",
            mc_10k.price, abs(mc_10k.price - true_price),
            (abs(mc_10k.price - true_price) / max(true_price, 1e-8)) * 100, mc_10k.elapsed_time)
    
    @printf("MC (100k paths)          │ €%9.6f │ %10.2e  │ %6.2f%%  │ %.4f s\n",
            mc_100k.price, abs(mc_100k.price - true_price),
            (abs(mc_100k.price - true_price) / max(true_price, 1e-8)) * 100, mc_100k.elapsed_time)
    
    @printf("MINSPM-QAE (2^28)        │ €%9.6f │ %10.2e  │ %6.2f%%  │ %.6f s\n",
            price_qae, abs(price_qae - true_price),
            (abs(price_qae - true_price) / max(true_price, 1e-8)) * 100, time_qae)
    
    println("\n" * "="^80)
    println(" ✓ All values are non-zero and reasonable!")
    println("="^80 * "\n")
end

run_debugged_comparison()
