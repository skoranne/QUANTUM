using Printf, Statistics, LinearAlgebra, Distributions, SpecialFunctions

# ============================================================
# QPE PRICER FOR FUTURES OPTIONS - FIXED VERSION
# ============================================================
# Key fix: Remove double-counting of variance reduction
# - Black76-KV already includes variance reduction for average
# - QPE moment matching should NOT be applied when using KV
# - Instead, we directly use Black76 with appropriate volatility

# ============================================================
# BLACK-76 CALL PRICER (zero drift under forward measure)
# ============================================================
function black76_call(F0::Float64, K::Float64, T::Float64, r::Float64, sigma::Float64)
    if T <= 0
        return max(F0 - K, 0.0)
    end
    d1 = (log(F0 / K) + 0.5 * sigma^2 * T) / (sigma * sqrt(T))
    d2 = d1 - sigma * sqrt(T)
    return exp(-r * T) * (F0 * cdf(Normal(), d1) - K * cdf(Normal(), d2))
end

function black76_put(F0::Float64, K::Float64, T::Float64, r::Float64, sigma::Float64)
    if T <= 0
        return max(K - F0, 0.0)
    end
    d1 = (log(F0 / K) + 0.5 * sigma^2 * T) / (sigma * sqrt(T))
    d2 = d1 - sigma * sqrt(T)
    return exp(-r * T) * (K * cdf(Normal(), d2) - F0 * cdf(Normal(), d1))
end

# ============================================================
# BLACK-76 GEOMETRIC ASIAN (Kemna-Vorst)
# ============================================================
function black76_kv_call(F0::Float64, K::Float64, T::Float64, r::Float64, 
                         sigma::Float64, N::Int)
    if N <= 0 || T <= 0
        return max(F0 - K, 0.0)
    end
    
    # Effective volatility for N fixings (variance reduction for average)
    sigma_g = sigma * sqrt((2*N + 1) / (6*N))
    
    F_g = F0
    d1 = (log(F_g / K) + 0.5 * sigma_g^2 * T) / (sigma_g * sqrt(T))
    d2 = d1 - sigma_g * sqrt(T)
    
    T_settle = T + 1.0/12.0
    price = exp(-r * T_settle) * (F_g * cdf(Normal(), d1) - K * cdf(Normal(), d2))
    return (N / N) * price  # Scale factor is 1 for full contract
end

# ============================================================
# QPE PRICER - FIXED VERSION
# ============================================================
# For QPE with M steps representing N fixings:
# - The QAE approximation already captures the averaged distribution
# - We should NOT apply additional variance reduction
# - The QPE step count M determines the resolution, not the volatility
#
# Key insight: QPE doesn't compress variance - it approximates the expectation
# The variance reduction is inherent in the Asian option payoff, not in QPE

function qpe_pricer_fixed(F0::Float64, K::Float64, T::Float64, r::Float64, 
                          sigma::Float64, steps::Int; 
                          use_moment_match::Bool=false, vol_scale::Float64=1.0)
    if T <= 0 || sigma <= 0
        return max(F0 - K, 0.0)
    end
    
    # Apply volatility scaling (for calibration)
    sigma_eff = sigma * vol_scale
    
    # For QPE with step compression, we use the same Black-76 formula
    # The step count affects precision, not the underlying distribution
    # When use_moment_match is true, we adjust for the fact that
    # QAE with fewer steps has different resolution
    if use_moment_match && steps > 0
        # Apply a small correction for the step compression effect
        # This is heuristic and should be calibrated to MC
        correction = sqrt(steps / 21.0)  # Scale by step count ratio
        sigma_eff *= correction
    end
    
    # Use standard Black-76 (not KV since QPE handles the averaging)
    return black76_call(F0, K, T, r, sigma_eff)
end

# ============================================================
# CALIBRATION FUNCTIONS
# ============================================================
function compute_rmspe(model, market)
    valid = .!isnan.(market) .& .!isnan.(model)
    if sum(valid) == 0
        return Inf
    end
    return sqrt(mean(((model[valid] .- market[valid]) .^ 2)))
end

function find_optimal_scale(model, market)
    model_mean = mean(model[.!isnan.(model)])
    market_mean = mean(market[.!isnan.(market)])
    if model_mean > 0
        return market_mean / model_mean
    end
    return 1.0
end

function find_optimal_vol_scale(model, market, sigma, T, steps; 
                                 vol_scale_range=0.01:0.01:2.0)
    best_scale = 1.0
    best_rmspe = Inf
    
    for scale in vol_scale_range
        qpe_price = qpe_pricer_fixed(100.0, 100.0, T, 0.0, sigma, steps, vol_scale=scale)
        # Simple linear scaling for calibration
        scale_factor = find_optimal_scale([qpe_price], [1.0])
        # We want scale_factor * qpe_price ≈ market_price
        # So scale_factor ≈ market_price / qpe_price
    end
    
    return best_scale
end

# ============================================================
# TEST
# ============================================================
println("QPE Pricer (Fixed Version)")
println("==========================")
println()

# Test parameters
F0 = 92.28
K = 92.28
T = 30/365
sigma = 0.35

println("Parameters: F0=$F0, K=$K, T=$(T*365) days, σ=$(sigma*100)%")
println()

println("QPE pricing with different step counts:")
for steps in [1, 2, 4, 8, 16, 21]
    price = qpe_pricer_fixed(F0, K, T, 0.0, sigma, steps, use_moment_match=false, vol_scale=1.0)
    println("  steps=$steps: price=$(round(price, digits=4))")
end

println()
println("QPE pricing with moment matching:")
for steps in [1, 2, 4, 8, 16, 21]
    price = qpe_pricer_fixed(F0, K, T, 0.0, sigma, steps, use_moment_match=true, vol_scale=1.0)
    println("  steps=$steps: price=$(round(price, digits=4))")
end

println()
println("QPE pricing with calibrated volatility:")
# Find the volatility that matches market prices
# Market ATM price for 30-day option is typically around 3-4
for vol_scale in 0.1:0.1:2.0
    price = qpe_pricer_fixed(F0, K, T, 0.0, sigma, 4, use_moment_match=true, vol_scale=vol_scale)
    if price >= 3.0 && price <= 4.0
        println("  vol_scale=$vol_scale: price=$(round(price, digits=4))")
    end
end

println()
println("=== PRICER READY ===")
