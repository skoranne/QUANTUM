using Printf, Statistics, LinearAlgebra, Distributions, SpecialFunctions, Random, CSV, DataFrames, Dates, Interpolations

# ============================================================
# IMPROVED QPE PRICER FOR BRENT CRUDE OIL OPTIONS
# ============================================================
# Key fixes:
# 1. Black-76 futures drift (zero drift under risk-neutral forward measure)
# 2. Seasoned averaging for contracts inside delivery month
# 3. Moment-matching variance reduction for QAE step compression
# 4. Proper discounting to final settlement date

# ============================================================
# BLACK-76 FUTURES OPTION PRICER (zero drift under forward measure)
# ============================================================
function black76_call(F0::Float64, K::Float64, T::Float64, r::Float64, sigma::Float64)
    # Black-76: futures have zero drift, discount from settlement date
    # dF_t = sigma * F_t * dW_t (under risk-neutral forward measure)
    if T <= 0
        return max(F0 - K, 0.0)
    end
    
    d1 = (log(F0 / K) + 0.5 * sigma^2 * T) / (sigma * sqrt(T))
    d2 = d1 - sigma * sqrt(T)
    
    # Discount from settlement date (first business day of following month)
    # For simplicity, assume T_settle = T + 1/12 (one month after expiry)
    T_settle = T + 1.0/12.0
    
    price = exp(-r * T_settle) * (F0 * cdf(Normal(), d1) - K * cdf(Normal(), d2))
    return price
end

# ============================================================
# SEASONED AVERAGING: Adjust strike for realized fixings
# ============================================================
function compute_seasoned_strike(K::Float64, A_realized::Float64, m::Int, N::Int)
    # K_eff = (N * K - m * A_realized) / (N - m)
    # For contracts inside the delivery month with m realized fixings
    if m >= N
        return 0.0, true  # All fixings completed
    end
    
    K_eff = (N * K - m * A_realized) / (N - m)
    is_deterministic = K_eff <= 0
    
    return K_eff, is_deterministic
end

# ============================================================
# MOMENT MATCHING: Match variance of discrete arithmetic average
# ============================================================
function match_moments_to_qae_steps(N_MC::Int, N_QAE::Int, sigma::Float64, T::Float64)
    # Compute effective volatility for QAE steps
    var_N = (N_MC + 1) * (2N_MC + 1) / (6N_MC^2)
    var_Q = (N_QAE + 1) * (2N_QAE + 1) / (6N_QAE^2)
    
    # Scale volatility to match variance
    sigma_eff = sigma * sqrt(var_N / var_Q)
    
    return sigma_eff, var_N, var_Q
end

# ============================================================
# BLACK-76 GEOMETRIC ASIAN PRICER (Kemna-Vorst with futures drift)
# ============================================================
function black76_kv_call(F0::Float64, K::Float64, T::Float64, r::Float64, 
                         sigma::Float64, N::Int, k::Int=0, S_bar_k::Float64=0.0)
    # M = remaining fixings
    M = N - k
    M <= 0 && return max(S_bar_k - K, 0.0)
    
    # Seasoned effective strike
    K_eff = (N * K - k * S_bar_k) / M
    K_eff <= 0.0 && return max(F0 - K, 0.0)
    
    # Effective volatility for M fixings
    sigma_g = sigma * sqrt((2M + 1) / (6M))
    
    # Black-76: zero drift for futures
    # F_g = F0 (no drift adjustment under forward measure)
    F_g = F0
    
    d1 = (log(F_g / K_eff) + 0.5 * sigma_g^2 * T) / (sigma_g * sqrt(T))
    d2 = d1 - sigma_g * sqrt(T)
    
    # Discount from settlement date
    T_settle = T + 1.0/12.0
    
    price = exp(-r * T_settle) * (F_g * cdf(Normal(), d1) - K_eff * cdf(Normal(), d2))
    return (M / N) * price
end

# ============================================================
# IMPROVED QPE PRICER WITH ALL FIXES
# ============================================================
struct QPEParams
    mc_factor::Float64          # Monte Carlo calibration factor
    vol_scale::Float64          # Volatility scaling
    price_scale::Float64        # Price scaling factor
    use_seasoned::Bool          # Use seasoned averaging
    use_moment_match::Bool      # Use moment-matched volatility
    N_fixings::Int              # Number of fixings in MC (typically 21)
end

function qpe_pricer_improved(F0::Float64, K::Float64, r::Float64, sigma::Float64, 
                            skew::Float64, T::Float64, steps::Int, params::QPEParams;
                            m_realized::Int=0, A_realized::Float64=0.0)
    # Step 1: Apply moment matching if enabled
    sigma_eff = sigma
    if params.use_moment_match && steps < params.N_fixings
        sigma_eff, _, _ = match_moments_to_qae_steps(params.N_fixings, steps, sigma, T)
    end
    
    # Step 2: Apply volatility scaling
    sigma_eff *= params.vol_scale
    
    # Step 3: Apply seasoned strike adjustment if contract is inside delivery month
    K_eff = K
    if params.use_seasoned && m_realized > 0
        K_eff, is_deterministic = compute_seasoned_strike(K, A_realized, m_realized, params.N_fixings)
        if is_deterministic
            # Payoff is deterministic - return intrinsic value
            return max(F0 - K, 0.0)
        end
    end
    
    # Step 4: Compute price using Black-76 KV formula (zero drift for futures)
    price = black76_kv_call(F0, K_eff, T, r, sigma_eff, steps)
    
    # Step 5: Apply Monte Carlo calibration factor
    price *= params.mc_factor
    
    # Step 6: Apply price scaling
    price *= params.price_scale
    
    return max(price, 0.0)
end

# ============================================================
# QPE PRICER DIRECT (simplified interface)
# ============================================================
function qpe_pricer_direct_improved(F0::Float64, K::Float64, r::Float64, sigma::Float64, 
                                   skew::Float64, T::Float64, steps::Int; 
                                   mc_factor::Float64=1.0,
                                   vol_scale::Float64=1.0,
                                   price_scale::Float64=1.0,
                                   use_seasoned::Bool=false,
                                   m_realized::Int=0,
                                   A_realized::Float64=0.0)
    params = QPEParams(mc_factor, vol_scale, price_scale, use_seasoned, true, 21)
    return qpe_pricer_improved(F0, K, r, sigma, skew, T, steps, params,
                              m_realized=m_realized, A_realized=A_realized)
end

# ============================================================
# CALIBRATION WITH IMPROVED MODELING
# ============================================================
function calibrate_improved(model_prices::Vector{Float64}, market_prices::Vector{Float64}, 
                           moneynesses::Vector{Float64})
    # Compute initial calibration using improved model
    model_mean = mean(model_prices)
    market_mean = mean(market_prices)
    
    println("  Model mean: $(round(model_mean, digits=4))")
    println("  Market mean: $(round(market_mean, digits=4))")
    
    if model_mean > 0
        global_scale = market_mean / model_mean
        println("  Global scale factor: $(round(global_scale, digits=4))")
        return global_scale
    end
    
    return 1.0
end

function calibrate_moneyness_dependent(model_prices::Vector{Float64}, market_prices::Vector{Float64}, 
                                      moneynesses::Vector{Float64})
    # Compute moneyness-dependent scaling
    # Bin by moneyness and compute scale for each bin
    
    bins = [0.85, 0.92, 0.98, 1.05, 1.15, 1.30]
    scales = Float64[]
    
    for i in 1:length(bins)-1
        bin_mask = (moneynesses .>= bins[i]) .& (moneynesses .< bins[i+1])
        if sum(bin_mask) >= 5
            bin_model = model_prices[bin_mask]
            bin_market = market_prices[bin_mask]
            
            model_mean = mean(bin_model)
            market_mean = mean(bin_market)
            
            if model_mean > 0
                scale = market_mean / model_mean
                push!(scales, scale)
            else
                push!(scales, 1.0)
            end
        else
            push!(scales, 1.0)
        end
    end
    
    println("  Bin scales: $(round.(scales, digits=2))")
    
    return scales, bins
end

# ============================================================
# VERIFICATION: Three-way comparison
# ============================================================
function verify_three_way(F0::Float64, K::Float64, r::Float64, sigma::Float64, T::Float64)
    # Black-76 call
    black76_price = black76_call(F0, K, T, r, sigma)
    
    # Spot GBM (incorrect for futures)
    d1_spot = (log(F0/K) + (r + 0.5*sigma^2)*T) / (sigma*sqrt(T))
    d2_spot = d1_spot - sigma*sqrt(T)
    spot_gbm_price = exp(-r*T) * (F0 * cdf(Normal(), d1_spot) - K * cdf(Normal(), d2_spot))
    
    println("Three-way comparison:")
    println("  Black-76 (futures, correct):    $(round(black76_price, digits=6))")
    println("  Spot GBM (equity, incorrect):   $(round(spot_gbm_price, digits=6))")
    println("  Difference:                     $(round(black76_price - spot_gbm_price, digits=6))")
    println("  Error: $(round(abs(black76_price - spot_gbm_price)/spot_gbm_price*100, digits=2))%")
    
    return black76_price, spot_gbm_price
end

# ============================================================
# LOAD DATA AND RUN BACKTEST
# ============================================================
println("Improved QPE Pricer for Brent Crude Oil Options")
println("==============================================")
println()

println("--- Verification Tests ---")
println()

println("1. Black-76 vs Spot GBM:")
verify_three_way(75.0, 75.0, 0.045, 0.35, 1.0/12.0)

println()
println("2. Moment matching:")
sigma_eff, var_N, var_Q = match_moments_to_qae_steps(21, 4, 0.35, 1.0/12.0)
println("  Original sigma: 0.35")
println("  Effective sigma: $(round(sigma_eff, digits=4))")
println("  Variance ratio: $(round(var_N/var_Q, digits=4))")

println()
println("3. Seasoned strike:")
K_eff, is_det = compute_seasoned_strike(75.0, 78.0, 10, 21)
println("  Original K: 75.0")
println("  Realized A: 78.0")
println("  Seasoned K_eff: $(round(K_eff, digits=4))")
println("  Is deterministic: $is_det")

println()
println("=== PRICER READY ===")

# Load and test with real data
println()
println("--- Testing with CB_BZ Data ---")

try
    df = DataFrame(CSV.File("/home/skoranne/PROGRAMS/QUANTUM/CB_BZ_option_chain_2026-11-16.csv"))
    println("Loaded $(size(df, 1)) contracts")
    
    # Filter meaningful trades
    mask = df.DTE .> 0 .& df.IV .> 0.01
    df_filtered = df[mask]
    println("Filtered to $(size(df_filtered, 1)) meaningful trades")
    
    # Build simple volatility surface
    spot = df_filtered.Spot_Price[1]
    tenors = sort(unique(df_filtered.DTE[df_filtered.DTE .> 0] / 365.0))
    strikes = sort(unique(df_filtered.Strike))
    println("Vol surface: $(length(tenors)) tenors x $(length(strikes)) strikes")
    
    # Test pricing on subset
    test_contracts = df_filtered[1:min(100, size(df_filtered, 1)), :]
    market_prices = test_contracts.Mid_Price
    implied_vols = test_contracts.IV
    moneynesses = [c.K / c.F0 for c in eachrow(test_contracts)]
    
    # Create params with improved settings
    params = QPEParams(0.78, 1.0, 1.0, true, true, 21)
    
    # Compute prices with improved model
    qpe_prices = Float64[]
    for c in eachrow(test_contracts)
        p = qpe_pricer_direct_improved(c.F0, c.K, c.r, c.IV, c.skew, c.T, 32,
                                      mc_factor=params.mc_factor,
                                      vol_scale=params.vol_scale,
                                      price_scale=params.price_scale,
                                      use_seasoned=params.use_seasoned,
                                      m_realized=0, A_realized=0.0)
        push!(qpe_prices, p)
    end
    
    # Compute RMSPE
    rmspe = sqrt(mean(abs.((qpe_prices .- market_prices) ./ max.(market_prices, 1e-8)).^2))
    println()
    println("Results:")
    println("  QPE mean price: $(round(mean(qpe_prices), digits=4))")
    println("  Market mean price: $(round(mean(market_prices), digits=4))")
    println("  RMSPE: $(round(rmspe*100, digits=2))%")
    
    # Try with global scale calibration
    global_scale = calibrate_improved(qpe_prices, market_prices, moneynesses)
    calibrated_prices = qpe_prices .* global_scale
    rmspe_calibrated = sqrt(mean(abs.((calibrated_prices .- market_prices) ./ max.(market_prices, 1e-8)).^2))
    println()
    println("After global scale calibration:")
    println("  Scale factor: $(round(global_scale, digits=4))")
    println("  RMSPE: $(round(rmspe_calibrated*100, digits=2))%")
    
catch e
    println("Error: $(e)")
end

println()
println("=== DONE ===")
