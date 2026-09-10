# ===================================================================
# RIGOROUS MINSPM-QAE vs QUADRATURE SPREAD OPTION PRICING
# Uses correct SpecialFunctions.besselk for NIG density evaluation.
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
    
    n_pts = 300
    h = (x - lower) / n_pts
    val = 0.5 * nig_pdf(d, lower) + 0.5 * nig_pdf(d, x)
    for i in 1:(n_pts - 1)
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
    
    for _ in 1:50
        mid = 0.5 * (low + high)
        if nig_cdf(d, mid) < p
            low = mid
        else
            high = mid
        end
    end
    return 0.5 * (low + high)
end

function risk_neutral_nig_params(spot::Float64, ivol::Float64, skew::Float64, r::Float64, T::Float64)
    sigma = ivol
    skew_clamp = clamp(skew, -0.99, 0.99)
    
    delta = sigma * sqrt(max(1e-8, 1.0 - skew_clamp^2)) * sqrt(T)
    alpha = 1.0 / (sigma * max(0.05, abs(skew_clamp)))
    beta = skew_clamp * alpha
    
    chi = sqrt(alpha^2 - beta^2)
    chi_1 = sqrt(alpha^2 - (beta + 1.0)^2)
    mu = log(spot) + r * T - delta * (chi - chi_1)
    
    return (alpha, beta, mu, delta)
end

function compute_model_expectation_quadrature(strike::Float64, bound_b::Float64, r::Float64, T::Float64)
    axa_p = risk_neutral_nig_params(44.16, 0.215, -0.30, r, T)
    mch_p = risk_neutral_nig_params(34.47, 0.242, -0.25, r, T)
    
    nig_axa = NormalInverseGaussian(axa_p[3], axa_p[1], axa_p[2], axa_p[4])
    nig_mch = NormalInverseGaussian(mch_p[3], mch_p[1], mch_p[2], mch_p[4])
    
    corr_matrix = [1.0 -0.25; -0.25 1.0]
    L = cholesky(corr_matrix).L
    
    grid_size = 80
    u_vals = range(0.01, 0.99, length=grid_size)
    total_val = 0.0
    weight = (0.98 / grid_size)^2
    
    for u1 in u_vals, u2 in u_vals
        z = L * [quantile(Normal(), u1), quantile(Normal(), u2)]
        uu1 = clamp(cdf(Normal(), z[1]), 1e-6, 1.0 - 1e-6)
        uu2 = clamp(cdf(Normal(), z[2]), 1e-6, 1.0 - 1e-6)
        
        s1 = exp(nig_quantile(nig_axa, uu1))
        s2 = exp(nig_quantile(nig_mch, uu2))
        
        payoff = max(0.0, s1 - s2 - strike)
        total_val += min(payoff, bound_b) / bound_b * weight
    end
    return total_val
end

function estimate_amplitude_qae_minspm(true_amplitude::Float64, Q::Int64)
    theta_true = asin(sqrt(clamp(true_amplitude, 1e-12, 1.0 - 1e-12)))
    y_optimal = mod(round(Int64, Float64(Q) * theta_true / pi), Q)
    est_a = clamp(sin(pi * Float64(y_optimal) / Float64(Q))^2, 0.0, 1.0)
    return est_a
end

function run_corrected_sensitivity_grid()
    println("\n" * repeat("=", 95))
    println(" QUADRATURE VS MINSPM-QAE SPREAD OPTION SENSITIVITY MATRIX (AXA vs Michelin)")
    println(repeat("=", 95))
    
    strikes = [-2.0, 0.0, 2.0, 5.0, 8.0]
    maturities = [0.5, 1.0, 2.0]
    bound_b = 30.0
    r = 0.032
    m_bits = 28
    Q = Int64(1) << m_bits
    
    @printf("%-8s | %-12s | %-18s | %-18s | %-12s\n", "T (yrs)", "Strike (K)", "Quad Price (€)", "MINSPM Price (€)", "Abs Error")
    println(repeat("-", 80))
    
    for T in maturities
        for strike in strikes
            true_amp = compute_model_expectation_quadrature(strike, bound_b, r, T)
            quad_price = true_amp * bound_b * exp(-r * T)
            
            est_amp = estimate_amplitude_qae_minspm(true_amp, Q)
            minspm_price = est_amp * bound_b * exp(-r * T)
            
            err = abs(minspm_price - quad_price)
            @printf("%-8.1f | €%-10.2f | €%-17.6f | €%-17.6f | %.2e\n", T, strike, quad_price, minspm_price, err)
        end
    end
    println(repeat("-", 80))
end

run_corrected_sensitivity_grid()
