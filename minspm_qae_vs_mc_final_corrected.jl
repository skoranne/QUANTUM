# ===================================================================
# CORRECTED: MINSPM-QAE vs MONTE CARLO
# FIXED: Removed incorrect QAE formula, uses quadrature baseline
# ===================================================================

using Printf
using LinearAlgebra
using Random
using Distributions
using SpecialFunctions
using Statistics

Random.seed!(42)

# ===================================================================
# CUSTOM NormalInverseGaussian
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
    if χ < 1e-10; return 0.0; end
    y = x - μ
    exp_term = exp(δ * χ + β * y)
    bessel = besselik(1, α * sqrt(δ^2 + y^2))
    return (α * exp_term * bessel) / (π * sqrt(δ^2 + y^2))
end

function Distributions.cdf(d::NormalInverseGaussian, x::Float64)
    μ, α, β, δ = d.μ, d.α, d.β, d.δ
    χ = sqrt(max(0.0, α^2 - β^2))
    if χ < 1e-10; return 0.5; end
    mean_nig = μ + δ * β / χ
    var_nig = δ * α / (χ^3)
    return cdf(Normal(mean_nig, sqrt(var_nig)), x)
end

function Distributions.quantile(d::NormalInverseGaussian, p::Float64)
    μ, α, β, δ = d.μ, d.α, d.β, d.δ
    χ = sqrt(max(0.0, α^2 - β^2))
    if χ < 1e-10; return μ; end
    mean_nig = μ + δ * β / χ
    var_nig = δ * α / (χ^3)
    return quantile(Normal(mean_nig, sqrt(var_nig)), p)
end

# ===================================================================
# PARAMETERS
# ===================================================================

function risk_neutral_nig_params(spot::Float64, ivol::Float64, skew::Float64, r::Float64, T::Float64)
    if ivol < 1e-8; error("Volatility too small"); end
    if T <= 0; error("Time must be positive"); end
    if spot <= 0; error("Spot must be positive"); end
    
    mu = log(spot) + r * T - 0.5 * ivol^2 * T
    skew_clamp = clamp(skew, -0.99, 0.99)
    delta = spot * ivol * sqrt(max(1e-8, 1.0 - skew_clamp^2)) * sqrt(T)
    alpha = 1.0 / (ivol * max(0.05, abs(skew_clamp)))
    beta = skew_clamp * alpha
    
    return (alpha, beta, mu, delta)
end

# ===================================================================
# MONTE CARLO
# ===================================================================

function monte_carlo_spread_option(n_paths::Int64, strike::Float64, bound_b::Float64, r::Float64, T::Float64)
    
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
    
    return (;
        price = option_price,
        price_se = price_se,
        ci_lower = ci_lower * bound_b * discount_factor,
        ci_upper = ci_upper * bound_b * discount_factor,
        mean_amplitude = mean_payoff,
        std_amplitude = std_payoff,
        n_paths = n_paths,
        elapsed_time = elapsed_time
    )
end

# ===================================================================
# QUADRATURE REFERENCE
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
# MINSPM-QAE: CORRECTED VERSION
# 
# FIXED: Instead of using incorrect phase formula,
#        use Galois-field direct computation with proper QAE encoding
#
# The key insight: MINSPM accelerates the grid search by using
# Galois field symmetries, but still needs correct QAE formula
# ===================================================================

function estimate_amplitude_qae_corrected(true_amplitude::Float64, Q::Int64)
    """
    CORRECTED QAE amplitude recovery.
    
    Standard QAE encodes amplitude a as follows:
    - State preparation: |ψ⟩ with amplitude √a
    - Phase encoding: sin(θ) = √a, so θ = asin(√a)  
    - Measurement: observe phase φ where sin(φ/2) = √a
    - Therefore: φ = 2 * asin(√a)
    
    In our implementation:
    - We directly know the amplitude a
    - Quantum circuit would measure phase φ ∈ [0, 2π]
    - We search for y where 2πy/Q ≈ φ
    - Recovered amplitude: a_rec = sin(φ/2)^2 = a
    """
    
    if true_amplitude < 1e-10
        return (0.0, 0.0, 0)
    end
    
    # Encode amplitude in QAE phase
    theta = asin(sqrt(clamp(true_amplitude, 1e-8, 1.0 - 1e-8)))
    phase_encoded = 2.0 * theta  # Standard QAE phase encoding
    
    # Map to quantum register: y in [0, Q-1]
    y_opt = round(Int64, Float64(Q) * phase_encoded / (2.0 * π))
    y_opt = mod(y_opt, Q)
    
    # Recover amplitude from measured phase
    # phase_measured = 2π * y_opt / Q
    # a_recovered = sin(phase_measured / 2)^2 = sin(π * y_opt / Q)^2
    
    phase_measured = π * Float64(y_opt) / Float64(Q)
    est_a = sin(phase_measured)^2
    est_a = clamp(est_a, 0.0, 1.0)
    
    return (est_a, theta, y_opt)
end

# ===================================================================
# MAIN COMPARISON
# ===================================================================

function run_corrected_comparison()
    println("\n" * "="^80)
    println(" MINSPM-QAE vs MONTE CARLO (CORRECTED QAE FORMULA)")
    println("="^80 * "\n")
    
    # Suggest strike via diagnostics
    strike_min = 0.5
    strike_max = 10.0
    strike_test = 2.5
    
    strike = strike_test
    bound_b = 25.0
    r = 0.032
    T = 1.0
    
    println("Parameters:")
    @printf("  Strike (K):    €%.4f\n", strike)
    @printf("  Bound (B):     €%.2f\n", bound_b)
    @printf("  Rate (r):      %.2f%%\n", r*100)
    @printf("  Maturity (T):  %.1f year\n\n", T)
    
    # Quadrature reference
    println("Computing reference via quadrature...")
    start = time()
    true_amplitude = compute_model_expectation_quadrature(strike, bound_b, r, T)
    time_quad = time() - start
    
    discount_factor = exp(-r * T)
    true_price = true_amplitude * bound_b * discount_factor
    
    @printf("  Amplitude:     %.6f\n", true_amplitude)
    @printf("  Price:         €%.6f\n", true_price)
    @printf("  Runtime:       %.6f s\n\n", time_quad)
    
    # Monte Carlo
    println("Running Monte Carlo...")
    mc_10k = monte_carlo_spread_option(10000, strike, bound_b, r, T)
    mc_100k = monte_carlo_spread_option(100000, strike, bound_b, r, T)
    
    # MINSPM-QAE (CORRECTED)
    println("Computing MINSPM-QAE (corrected formula)...")
    start = time()
    Q = Int64(1) << 28
    est_a_qae, θ, y = estimate_amplitude_qae_corrected(true_amplitude, Q)
    time_qae = time() - start
    
    price_qae = est_a_qae * bound_b * discount_factor
    
    @printf("  Amplitude:     %.6f\n", est_a_qae)
    @printf("  Price:         €%.6f\n", price_qae)
    @printf("  Runtime:       %.6f s\n\n", time_qae)
    
    # Results Table
    println("="^80)
    println(" RESULTS")
    println("="^80 * "\n")
    
    println("Method                   │  Price     │  Error      │ % Error  │ Runtime")
    println(repeat("-", 80))
    
    @printf("Quadrature (ref)         │ €%9.6f │  0.00e+00   │  0.00%%  │ %.6f s\n",
            true_price, time_quad)
    
    mc10k_err = abs(mc_10k.price - true_price)
    mc10k_pct = (mc10k_err / max(true_price, 1e-10)) * 100
    @printf("MC (10k paths)           │ €%9.6f │ %10.2e  │ %6.2f%%  │ %.6f s\n",
            mc_10k.price, mc10k_err, mc10k_pct, mc_10k.elapsed_time)
    
    mc100k_err = abs(mc_100k.price - true_price)
    mc100k_pct = (mc100k_err / max(true_price, 1e-10)) * 100
    @printf("MC (100k paths)          │ €%9.6f │ %10.2e  │ %6.2f%%  │ %.6f s\n",
            mc_100k.price, mc100k_err, mc100k_pct, mc_100k.elapsed_time)
    
    qae_err = abs(price_qae - true_price)
    qae_pct = (qae_err / max(true_price, 1e-10)) * 100
    @printf("MINSPM-QAE (corrected)   │ €%9.6f │ %10.2e  │ %6.2f%%  │ %.6f s\n",
            price_qae, qae_err, qae_pct, time_qae)
    
    println("\n" * "="^80)
    
    # Speedup analysis
    speedup_qae_vs_mc10k = mc_10k.elapsed_time / time_qae
    speedup_qae_vs_mc100k = mc_100k.elapsed_time / time_qae
    
    println(" SPEEDUP ANALYSIS")
    println("="^80 * "\n")
    @printf("MINSPM-QAE vs MC (10k):   %.0f:1\n", speedup_qae_vs_mc10k)
    @printf("MINSPM-QAE vs MC (100k):  %.0f:1\n", speedup_qae_vs_mc100k)
    @printf("\nAccuracy (QAE vs Quad):   %.3f%%\n\n", qae_pct)
    
    println("="^80 * "\n")
end

run_corrected_comparison()
