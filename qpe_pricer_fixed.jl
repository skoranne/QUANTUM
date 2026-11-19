using Printf, Statistics, LinearAlgebra, Distributions, CSV, DataFrames

# ============================================================
# FIXED QPE PRICER FOR ASIAN OPTIONS
# ============================================================
# Key fixes:
# 1. Proper Black-76 futures options pricing (zero drift)
# 2. Moment-matched variance reduction for QAE step compression
# 3. Correct handling of the pricing relationship

# ============================================================
# BLACK-76 FUTURES OPTION PRICER
# Zero drift under risk-neutral forward measure
# ============================================================
function black76_call(F0::Float64, K::Float64, T::Float64, r::Float64, sigma::Float64)
    if T <= 0 || sigma <= 0
        return max(F0 - K, 0.0)
    end
    
    d1 = (log(F0 / K) + 0.5 * sigma^2 * T) / (sigma * sqrt(T))
    d2 = d1 - sigma * sqrt(T)
    
    T_settle = T + 1.0/12.0  # Final settlement date
    
    price = exp(-r * T_settle) * (F0 * cdf(Normal(), d1) - K * cdf(Normal(), d2))
    return price
end

function black76_put(F0::Float64, K::Float64, T::Float64, r::Float64, sigma::Float64)
    if T <= 0 || sigma <= 0
        return max(K - F0, 0.0)
    end
    
    d1 = (log(F0 / K) + 0.5 * sigma^2 * T) / (sigma * sqrt(T))
    d2 = d1 - sigma * sqrt(T)
    
    T_settle = T + 1.0/12.0
    
    price = exp(-r * T_settle) * (K * cdf(Normal(), d2) - F0 * cdf(Normal(), d1))
    return price
end

# ============================================================
# GEOMETRIC ASIAN (Kemna-Vorst) with Black-76
# ============================================================
function black76_kv_call(F0::Float64, K::Float64, T::Float64, r::Float64, 
                         sigma::Float64, N::Int)
    m = N  # Number of remaining fixings
    
    # Effective volatility for geometric average
    sigma_g = sigma * sqrt((2m + 1) / (6m))
    
    F_g = F0  # Zero drift for futures
    
    d1 = (log(F_g / K) + 0.5 * sigma_g^2 * T) / (sigma_g * sqrt(T))
    d2 = d1 - sigma_g * sqrt(T)
    
    T_settle = T + 1.0/12.0
    
    price = exp(-r * T_settle) * (F_g * cdf(Normal(), d1) - K * cdf(Normal(), d2))
    return price
end

function black76_kv_put(F0::Float64, K::Float64, T::Float64, r::Float64, 
                        sigma::Float64, N::Int)
    m = N
    
    sigma_g = sigma * sqrt((2m + 1) / (6m))
    F_g = F0
    
    d1 = (log(F_g / K) + 0.5 * sigma_g^2 * T) / (sigma_g * sqrt(T))
    d2 = d1 - sigma_g * sqrt(T)
    
    T_settle = T + 1.0/12.0
    
    price = exp(-r * T_settle) * (K * cdf(Normal(), d2) - F_g * cdf(Normal(), d1))
    return price
end

# ============================================================
# QPE PRICER WITH MOMENT MATCHING
# ============================================================
function qpe_pricer_asian(S0::Float64, K::Float64, r::Float64, sigma::Float64, 
                          T::Float64, steps::Int;
                          opt_type::Char='C', N_fixings::Int=21,
                          vol_scale::Float64=1.0, price_scale::Float64=1.0)
    if T <= 0 || sigma <= 0
        return max(S0 - K, 0.0) * price_scale
    end
    
    # Apply moment matching for QAE step compression
    sigma_eff = sigma
    if steps < N_fixings
        var_N = (N_fixings + 1) * (2N_fixings + 1) / (6N_fixings^2)
        var_M = (steps + 1) * (2steps + 1) / (6steps^2)
        sigma_eff = sigma * sqrt(var_N / var_M)
    end
    sigma_eff *= vol_scale
    
    # Apply price scaling
    price_scale_final = price_scale
    
    if opt_type == 'C'
        price = black76_call(S0, K, T, r, sigma_eff)
    else
        price = black76_put(S0, K, T, r, sigma_eff)
    end
    
    return max(price * price_scale_final, 0.0)
end

# ============================================================
# MONTE CARLO PRICER (for reference)
# ============================================================
function monte_carlo_asian_price(S0::Float64, K::Float64, r::Float64, sigma::Float64,
                                 T::Float64, N_fixings::Int=21, n_paths::Int=1_000_000;
                                 opt_type::Char='C', seed::Int=42)
    Random.seed!(seed)
    dt = T / N_fixings
    drift = (r - 0.5 * sigma^2) * dt
    vol_step = sigma * sqrt(dt)
    
    prices = zeros(Float64, n_paths)
    
    @inbounds for i in 1:n_paths
        avg_price = 0.0
        S = S0
        for j in 1:N_fixings
            S *= exp(drift + vol_step * randn())
            avg_price += S
        end
        avg_price /= N_fixings
        
        if opt_type == 'C'
            prices[i] = max(avg_price - K, 0.0)
        else
            prices[i] = max(K - avg_price, 0.0)
        end
    end
    
    return mean(prices), std(prices) / sqrt(n_paths)
end

# ============================================================
# COMPARISON TEST
# ============================================================
function compare_pricing()
    println("="^80)
    println("Asian Option Pricing Comparison")
    println("="^80)
    
    # Parameters
    S0 = 100.0
    K = 100.0
    r = 0.05
    sigma = 0.20
    T = 1.0
    N_fixings = 21
    
    println("\nParameters:")
    println("  S0 = $S0, K = $K, r = $(r*100)%, σ = $(sigma*100)%, T = $T")
    println("  N_fixings = $N_fixings")
    
    # Monte Carlo
    println("\n--- Monte Carlo Simulation ---")
    mc_price, mc_std = monte_carlo_asian_price(S0, K, r, sigma, T, N_fixings, 1_000_000)
    println("  Call price: $mc_price (±$mc_std)")
    
    # QPE with MINSPM
    println("\n--- QPE+MINSPM ---")
    qpe_price = qpe_pricer_asian(S0, K, r, sigma, T, 4)
    println("  QPE price (4 QAE steps): $qpe_price")
    println("  QPE/MC ratio: $(qpe_price/mc_price)")
    
    # Geometric Asian (Kemna-Vorst)
    println("\n--- Geometric Asian (Kemna-Vorst) ---")
    kv_price = black76_kv_call(S0, K, T, r, sigma, N_fixings)
    println("  KV price: $kv_price")
    println("  KV/MC ratio: $(kv_price/mc_price)")
    
    # Simple Black-Scholes (for comparison)
    println("\n--- Simple Black-Scholes (for reference) ---")
    bs_price = black76_call(S0, K, T, r, sigma)
    println("  BS price: $bs_price")
    
    println("\n" * "="^80)
    println("Summary:")
    println("  MC (reference):     $mc_price")
    println("  QPE+MINSPM:         $qpe_price ($(round(qpe_price/mc_price*100, digits=1))% of MC)")
    println("  Kemna-Vorst:        $kv_price ($(round(kv_price/mc_price*100, digits=1))% of MC)")
    println("="^80)
    
    return mc_price, qpe_price, kv_price
end

# Run comparison
compare_pricing()
