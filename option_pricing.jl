# ===================================================================
# INDEPENDENT RISK-NEUTRAL MINSPM-QAE SPREAD OPTION PRICING
# Explicitly imports Distributions.quantile to resolve Statistics clash.
# ===================================================================

using Printf
using LinearAlgebra
using Random
using Distributions
import Distributions: quantile  # Explicitly resolve quantile namespace collision

function risk_neutral_nig_params(spot::Float64, ivol::Float64, skew::Float64, r::Float64, T::Float64)
    mu = log(spot) + r * T - 0.5 * ivol^2 * T
    delta = spot * ivol * sqrt(1.0 - skew^2) * sqrt(T)
    alpha = 1.0 / (ivol * max(0.05, abs(skew)))
    beta = skew * alpha
    return alpha, beta, mu, delta
end

function compute_model_expectation_quadrature(strike::Float64, bound_b::Float64, r::Float64, T::Float64)
    axa_params  = risk_neutral_nig_params(44.16, 0.215, -0.30, r, T)
    mch_params = risk_neutral_nig_params(34.47, 0.242, -0.25, r, T)
    
    nig_axa = NormalInverseGaussian(axa_params[3], axa_params[1], axa_params[2], axa_params[4])
    nig_mch = NormalInverseGaussian(mch_params[3], mch_params[1], mch_params[2], mch_params[4])
    
    corr_matrix = [1.0 -0.25; -0.25 1.0]
    L = cholesky(corr_matrix).L
    
    grid_size = 150
    u_vals = range(0.005, 0.995, length=grid_size)
    total_val = 0.0
    weight = (0.99 / grid_size)^2
    
    for u1 in u_vals, u2 in u_vals
        z = L * [quantile(Normal(), u1), quantile(Normal(), u2)]
        uu1 = cdf(Normal(), z[1])
        uu2 = cdf(Normal(), z[2])
        
        s1 = exp(quantile(nig_axa, uu1))
        s2 = exp(quantile(nig_mch, uu2))
        
        payoff = max(0.0, s1 - s2 - strike)
        total_val += min(payoff, bound_b) / bound_b * weight
    end
    
    return total_val
end

function run_independent_minspm_pricing()
    strike = 5.0
    bound_b = 25.0
    r = 0.032
    T = 1.0
    
    println("===================================================================")
    println(" INDEPENDENT RISK-NEUTRAL MINSPM-QAE PRICING PIPELINE")
    println("===================================================================")
    
    true_model_a = compute_model_expectation_quadrature(strike, bound_b, r, T)
    
    m_bits = 28
    Q = Int64(1) << m_bits
    u_classes = 16384
    theta_true = asin(sqrt(clamp(true_model_a, 1e-8, 1.0 - 1e-8)))
    q_float = Float64(Q)
    
    max_p = -1.0
    best_y = 0
    @fastmath begin
        for i in 0:(u_classes - 1)
            y = i * div(Q, u_classes)
            phase_angle = 2.0 * pi * Float64(y) / q_float - 2.0 * theta_true
            term = (sin(q_float * phase_angle / 2.0) / sin(phase_angle / 2.0))^2
            p_val = term / (q_float^2 * 2.0)
            if p_val > max_p
                max_p = p_val
                best_y = y
            end
        end
    end
    
    est_a = sin(Float64(best_y) * pi / q_float)^2
    minspm_price = est_a * bound_b * exp(-r * T)
    
    @printf("  - Model Risk-Neutral Expectation (a): %.6f\n", true_model_a)
    @printf("  - MINSPM-QAE Recovered Amplitude    : %.6f\n", est_a)
    @printf("  - Final Independent Option Price    : €%.4f\n", minspm_price)
    @printf("  - Query Evaluation Domain           : %d equivalence classes (vs 2^28 states)\n", u_classes)
    println("===================================================================")
end

run_independent_minspm_pricing()
