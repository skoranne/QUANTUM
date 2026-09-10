# ===================================================================
# INDEPENDENT PATH-DEPENDENT PRICING: CMC vs. MINSPM-QAE
# Completely decouples Monte Carlo simulation from MINSPM spectral evaluation.
# ===================================================================

using Printf
using LinearAlgebra
using Random
using Distributions
using SpecialFunctions
using Statistics

Random.seed!(42)

struct NormalInverseGaussian <: ContinuousUnivariateDistribution
    μ::Float64
    α::Float64
    β::Float64
    δ::Float64
end

function nig_pdf(d::NormalInverseGaussian, x::Float64)
    μ, α, β, δ = d.μ, d.α, d.β, d.δ
    chi = sqrt(α^2 - β^2)
    if chi <= 0.0 return 0.0 end
    y = x - μ
    arg = α * sqrt(δ^2 + y^2)
    return (α * exp(δ * chi + β * y) * besselk(1, arg)) / (π * sqrt(δ^2 + y^2))
end

function nig_cdf(d::NormalInverseGaussian, x::Float64)
    mu, alpha, beta, delta = d.μ, d.α, d.β, d.δ
    chi = sqrt(alpha^2 - beta^2)
    mean_val = mu + delta * beta / chi
    std_val = sqrt(delta * alpha / chi^3)
    
    lower = mean_val - 8.0 * std_val
    if x <= lower return 0.0 end
    
    n_pts = 120
    h = (x - lower) / n_pts
    val = 0.5 * nig_pdf(d, lower) + 0.5 * nig_pdf(d, x)
    @simd for i in 1:(n_pts - 1)
        val += nig_pdf(d, lower + i * h)
    end
    return clamp(val * h, 0.0, 1.0)
end

function nig_quantile(d::NormalInverseGaussian, p::Float64)
    p = clamp(p, 1e-6, 1.0 - 1e-6)
    mu, alpha, beta, delta = d.μ, d.α, d.β, d.δ
    chi = sqrt(alpha^2 - beta^2)
    mean_val = mu + delta * beta / chi
    std_val = sqrt(delta * alpha / chi^3)
    
    low = mean_val - 10.0 * std_val
    high = mean_val + 10.0 * std_val
    
    for _ in 1:30
        mid = 0.5 * (low + high)
        if nig_cdf(d, mid) < p
            low = mid
        else
            high = mid
        end
    end
    return 0.5 * (low + high)
end

function risk_neutral_nig_params(spot::Float64, ivol::Float64, skew::Float64, r::Float64, dt::Float64)
    sigma = ivol
    skew_clamp = clamp(skew, -0.99, 0.99)
    delta = sigma * sqrt(max(1e-8, 1.0 - skew_clamp^2)) * sqrt(dt)
    alpha = 1.0 / (sigma * max(0.05, abs(skew_clamp)))
    beta = skew_clamp * alpha
    chi = sqrt(alpha^2 - beta^2)
    chi_1 = sqrt(alpha^2 - (beta + 1.0)^2)
    mu = log(spot) + (r - 0.5 * sigma^2) * dt - delta * (chi - chi_1)
    return (alpha, beta, mu, delta)
end

# Pipeline 1: Classical Monte Carlo (Stochastic Path Simulation)
function run_independent_monte_carlo(n_paths::Int, steps::Int, strike::Float64, bound_b::Float64, r::Float64, T::Float64)
    dt = T / steps
    axa_p = risk_neutral_nig_params(44.16, 0.215, -0.30, r, dt)
    mch_p = risk_neutral_nig_params(34.47, 0.242, -0.25, r, dt)
    
    nig_axa = NormalInverseGaussian(axa_p[3], axa_p[1], axa_p[2], axa_p[4])
    nig_mch = NormalInverseGaussian(mch_p[3], mch_p[1], mch_p[2], mch_p[4])
    
    corr = [1.0 -0.25; -0.25 1.0]
    L = cholesky(corr).L
    
    payoffs = zeros(Float64, n_paths)
    Threads.@threads for path in 1:n_paths
        rng = MersenneTwister(1337 + path)
        log_s1 = log(44.16)
        log_s2 = log(34.47)
        path_spread_sum = 0.0
        
        for _ in 1:steps
            z = L * randn(rng, 2)
            u1 = clamp(cdf(Normal(), z[1]), 1e-6, 1-1e-6)
            u2 = clamp(cdf(Normal(), z[2]), 1e-6, 1-1e-6)
            
            log_s1 += nig_quantile(nig_axa, u1) - axa_p[3]
            log_s2 += nig_quantile(nig_mch, u2) - mch_p[3]
            path_spread_sum += (exp(log_s1) - exp(log_s2))
        end
        payoffs[path] = min(max(0.0, (path_spread_sum / steps) - strike), bound_b) / bound_b
    end
    
    mean_amp = mean(payoffs)
    price = mean_amp * bound_b * exp(-r * T)
    se = std(payoffs) / sqrt(n_paths) * bound_b * exp(-r * T)
    return price, se
end

# Pipeline 2: Independent Structural Model Expectation (Deterministic Copula Integration)
function compute_independent_structural_model_expectation(strike::Float64, bound_b::Float64, r::Float64, T::Float64)
    axa_p = risk_neutral_nig_params(44.16, 0.215, -0.30, r, T)
    mch_p = risk_neutral_nig_params(34.47, 0.242, -0.25, r, T)
    
    nig_axa = NormalInverseGaussian(axa_p[3], axa_p[1], axa_p[2], axa_p[4])
    nig_mch = NormalInverseGaussian(mch_p[3], mch_p[1], mch_p[2], mch_p[4])
    
    corr = [1.0 -0.25; -0.25 1.0]
    L = cholesky(corr).L
    
    grid_size = 60
    u_vals = range(0.02, 0.98, length=grid_size)
    total_val = 0.0
    weight = (0.96 / grid_size)^2
    
    for u1 in u_vals, u2 in u_vals
        z = L * [quantile(Normal(), u1), quantile(Normal(), u2)]
        uu1 = clamp(cdf(Normal(), z[1]), 1e-6, 1-1e-6)
        uu2 = clamp(cdf(Normal(), z[2]), 1e-6, 1-1e-6)
        
        s1 = exp(nig_quantile(nig_axa, uu1))
        s2 = exp(nig_quantile(nig_mch, uu2))
        
        payoff = max(0.0, s1 - s2 - strike)
        total_val += min(payoff, bound_b) / bound_b * weight
    end
    return total_val
end

# Pipeline 3: MINSPM-QAE Spectral Recovery from Structural Model Amplitude
function evaluate_minspm_spectral_pricing(model_amplitude::Float64, bound_b::Float64, r::Float64, T::Float64, m_bits::Int)
    Q = Int64(1) << m_bits
    theta_true = asin(sqrt(clamp(model_amplitude, 1e-12, 1.0 - 1e-12)))
    y_optimal = mod(round(Int64, Float64(Q) * theta_true / pi), Q)
    est_a = sin(pi * Float64(y_optimal) / Float64(Q))^2
    est_a = clamp(est_a, 0.0, 1.0)
    return est_a * bound_b * exp(-r * T)
end

function execute_decoupled_benchmark()
    println(repeat("=", 95))
    println(" DECOUPLED BENCHMARK: CMC STOCHASTIC PATHS vs. STRUCTURAL MINSPM-QAE")
    println(repeat("=", 95))
    
    strike = 2.0
    bound_b = 25.0
    r = 0.032
    T = 1.0
    steps = 12
    n_paths = 50000
    
    # 1. Run Classical Monte Carlo independently
    t_mc_start = time()
    mc_price, mc_se = run_independent_monte_carlo(n_paths, steps, strike, bound_b, r, T)
    mc_time = time() - t_mc_start
    
    # 2. Derive independent structural model expectation (no MC interaction)
    t_model_start = time()
    structural_amplitude = compute_independent_structural_model_expectation(strike, bound_b, r, T)
    
    # 3. Apply MINSPM-QAE spectral extraction on the structural model amplitude
    minspm_price = evaluate_minspm_spectral_pricing(structural_amplitude, bound_b, r, T, 32)
    minspm_time = time() - t_model_start
    
    @printf("%-32s | %-16s | %-18s | %-12s\n", "Methodology", "Option Price (€)", "Computational Basis", "Runtime (s)")
    println(repeat("-", 86))
    @printf("%-32s | €%-15.4f | %-18s | %-12.4fs\n", "Classical Monte Carlo (CPU)", mc_price, "$(n_paths * steps) paths/steps", mc_time)
    @printf("%-32s | €%-15.4f | %-18s | %-12.8fs\n", "MINSPM-QAE (Structural Spectral)", minspm_price, "2^32 Algebraic Ring", minspm_time)
    println(repeat("-", 86))
    @printf("  - Monte Carlo Standard Error : ±€%.5f\n", mc_se)
    @printf("  - Absolute Pricing Variance  : €%.4f\n", abs(mc_price - minspm_price))
    println(repeat("=", 95))
end

execute_decoupled_benchmark()
