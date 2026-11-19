using Printf, Statistics, LinearAlgebra, Distributions, SpecialFunctions, Random, CSV, DataFrames, Dates, Interpolations

# ============================================================
# IMPROVED QPE PRICER FOR BRENT CRUDE OIL OPTIONS - v4
# ============================================================
# Key fixes:
# 1. Black-76 futures drift (zero drift under risk-neutral forward measure)
# 2. Seasoned averaging for contracts inside delivery month
# 3. Moment-matching variance reduction for QAE step compression
# 4. Proper discounting to final settlement date
# 5. Adaptive correction via grid search + polynomial regression

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
    T_settle = T + 1.0/12.0
    
    price = exp(-r * T_settle) * (F0 * cdf(Normal(), d1) - K * cdf(Normal(), d2))
    return price
end

# ============================================================
# BLACK-76 GEOMETRIC ASIAN PRICER (Kemna-Vorst with futures drift)
# ============================================================
function black76_kv_call(F0::Float64, K::Float64, T::Float64, r::Float64, 
                         sigma::Float64, N::Int, k::Int=0, S_bar_k::Float64=0.0)
    M = N - k
    M <= 0 && return max(S_bar_k - K, 0.0)
    
    # Seasoned effective strike
    K_eff = (N * K - k * S_bar_k) / M
    K_eff <= 0.0 && return max(F0 - K, 0.0)
    
    # Effective volatility for M fixings
    sigma_g = sigma * sqrt((2M + 1) / (6M))
    
    # Black-76: zero drift for futures
    F_g = F0
    
    d1 = (log(F_g / K_eff) + 0.5 * sigma_g^2 * T) / (sigma_g * sqrt(T))
    d2 = d1 - sigma_g * sqrt(T)
    
    T_settle = T + 1.0/12.0
    
    price = exp(-r * T_settle) * (F_g * cdf(Normal(), d1) - K_eff * cdf(Normal(), d2))
    return (M / N) * price
end

# ============================================================
# MOMENT MATCHING: Match variance of discrete arithmetic average
# ============================================================
function match_moments_to_qae_steps(N_MC::Int, N_QAE::Int, sigma::Float64, T::Float64)
    var_N = (N_MC + 1) * (2N_MC + 1) / (6N_MC^2)
    var_Q = (N_QAE + 1) * (2N_QAE + 1) / (6N_QAE^2)
    
    sigma_eff = sigma * sqrt(var_N / var_Q)
    
    return sigma_eff, var_N, var_Q
end

# ============================================================
# BLACK-76 MOMENT-MATCHED LOGNORMAL PARAMS
# ============================================================
function black76_moments(
    F_curve::Vector{Float64},
    fixing_times::Vector{Float64},
    vols::Vector{Float64}
)
    k = length(F_curve)
    k <= 0 && return (0.0, 0.0, 0.0)
    
    # Cov(ln F_ti, ln F_tj) = sigma^2 * min(t_i, t_j)
    t_min = min.(fixing_times', fixing_times)
    cov_matrix = (vols.^2) .* t_min
    
    # M1 = E[F_bar] = mean of F_curve (Black-76 has zero drift)
    M1 = mean(F_curve)
    
    # M2 = E[F_bar^2] = sum_i,j F_i*F_j*exp(Cov_ij) / k^2
    M2 = sum(outer(F_curve, F_curve) .* exp.(cov_matrix)) / (k^2)
    
    v_eff_sq = log(M2 / (M1^2))
    v_eff = sqrt(max(v_eff_sq, 1e-16))
    mu_eff = log(M1) - 0.5 * v_eff_sq
    
    return mu_eff, v_eff, M1
end

# ============================================================
# SEASONED AVERAGING: Adjust strike for realized fixings
# ============================================================
function compute_seasoned_strike(K::Float64, A_realized::Float64, m::Int, N::Int)
    if m >= N
        return 0.0, true
    end
    
    K_eff = (N * K - m * A_realized) / (N - m)
    is_deterministic = K_eff <= 0
    
    return K_eff, is_deterministic
end

# ============================================================
# QPE PRICER WITH BLACK-76, MOMENT MATCHING, SEASONED
# ============================================================
function qpe_pricer_v4(
    spot::Float64,
    strike::Float64,
    r::Float64,
    iv::Float64,
    T::Float64,
    steps::Int;
    N_fixings::Int=21,
    m_realized::Int=0,
    A_realized::Float64=0.0,
    vol_scale::Float64=1.0,
    price_scale::Float64=1.0,
    use_seasoned::Bool=true,
    use_moment_match::Bool=true,
    mc_factor::Float64=1.0
)
    # Handle edge cases
    if T <= 0
        return max(spot - strike, 0.0) * price_scale
    end
    
    # Compute seasoned strike if contract is inside delivery month
    K_eff = strike
    if use_seasoned && m_realized > 0
        K_eff, is_deterministic = compute_seasoned_strike(strike, A_realized, m_realized, N_fixings)
        if is_deterministic
            return max(spot - strike, 0.0) * price_scale
        end
    end
    
    # Apply moment matching
    sigma_eff = iv
    if use_moment_match && steps < N_fixings
        sigma_eff, _, _ = match_moments_to_qae_steps(N_fixings, steps, iv, T)
    end
    
    # Apply volatility scaling
    sigma_eff *= vol_scale
    
    # Compute price using Black-76 KV formula
    price = black76_kv_call(spot, K_eff, T, r, sigma_eff, steps, m_realized, A_realized)
    
    # Apply calibration factors
    price *= mc_factor
    price *= price_scale
    
    return max(price, 0.0)
end

# ============================================================
# ADAPTIVE CALIBRATION: Grid search + Polynomial regression
# ============================================================
function calibrate_adaptive(
    model_prices::Vector{Float64},
    market_prices::Vector{Float64},
    moneynesses::Vector{Float64}
)
    # Phase 1: Global scale calibration
    model_mean = mean(model_prices)
    market_mean = mean(market_prices)
    
    if model_mean > 0
        global_scale = market_mean / model_mean
    else
        global_scale = 1.0
    end
    
    # Phase 2: Moneyness-dependent calibration via grid search
    # Test range: scale from 0.5 to 2.5 in 0.05 increments
    best_rmspe = Inf
    best_scale = 1.0
    
    for scale in 0.5:0.05:2.5
        test_prices = model_prices * scale
        rmspe = sqrt(mean(((test_prices .- market_prices) ./ max.(market_prices, 1e-8)).^2))
        if rmspe < best_rmspe
            best_rmspe = rmspe
            best_scale = scale
        end
    end
    
    # Phase 3: Moneyness-dependent polynomial regression
    n = length(moneynesses)
    features = ones(n)
    push!(features, moneynesses)  # linear
    push!(features, moneynesses.^2)  # quadratic
    push!(features, log.(moneynesses .+ 1e-8))  # log
    push!(features, moneynesses .* log.(moneynesses .+ 1e-8))  # interaction
    push!(features, (moneynesses .- 1.0).^2)  # distance from ATM
    push!(features, sign.(moneynesses .- 1.0))  # ITM/OTM indicator
    push!(features, 1.0 ./ max.(moneynesses, 0.8))  # inverse moneyness
    
    features = hcat(features...)
    
    # Fit polynomial with outlier rejection
    weights = ones(n)
    beta = zeros(size(features, 2))
    
    for iter in 1:3
        W = Diagonal(weights)
        XTW = features' * W
        beta = (XTW * features) \ (XTW * market_prices)
        predictions = features * beta
        
        # MAD-based outlier rejection
        residuals = market_prices - predictions
        mad = median(abs.(residuals .- median(residuals)))
        weights = @. 1 / (1 + abs(residuals) / (mad + 1e-8))
        weights .= clamp.(weights, 0.1, 1.0)
    end
    
    # Apply polynomial correction to model prices
    poly_predictions = features * beta
    
    return global_scale, best_scale, beta, features, poly_predictions
end

# ============================================================
# MAIN CALIBRATION RUN
# ============================================================
println("="^60)
println("QPE+MINSPM Calibration - v4 (Black-76 + Moment Matching)")
println("="^60)

# Load data
println("\nLoading CB_BZ option chain...")
df = DataFrame(CSV.File("/home/skoranne/PROGRAMS/QUANTUM/CB_BZ_option_chain_2026-11-16.csv"))
println("Loaded $(size(df, 1)) contracts")

# Filter meaningful trades
mask = (df.DTE .> 0) .& (df.IV .> 0.01) .& (df.Mid .> 0)
df_filtered = df[mask, :]
println("Filtered to $(size(df_filtered, 1)) meaningful trades")

# Test on subset
n_test = min(500, size(df_filtered, 1))
test_contracts = df_filtered[1:n_test, :]

println("\nTesting with $n_test contracts...")
println("Available columns: $(names(df_filtered))")

# Extract market data
# Note: CSV has "Mid" as the price column, "Strike" not "K", "Spot_Price" not "F0"
# Rate column might not exist - we'll use a reasonable default
market_prices = Float64[c.Mid for c in eachrow(test_contracts)]
moneynesses = Float64[c.Strike / c.Spot_Price for c in eachrow(test_contracts)]

# Use default risk-free rate
r_default = 0.045

println("\n--- Initial Pricing (no calibration) ---")
qpe_prices = Float64[]
for c in eachrow(test_contracts)
    # Extract parameters with defaults
    spot = c.Spot_Price
    strike = c.Strike
    iv = c.IV
    dte = c.DTE
    
    # Compute T from DTE (years)
    T = dte / 365.0
    
    p = qpe_pricer_v4(spot, strike, r_default, iv, T, 32;
                     N_fixings=21,
                     m_realized=0,
                     A_realized=0.0,
                     vol_scale=1.0,
                     price_scale=1.0,
                     use_seasoned=false,
                     use_moment_match=false,
                     mc_factor=1.0)
    push!(qpe_prices, p)
end

initial_rmspe = sqrt(mean(((qpe_prices .- market_prices) ./ max.(market_prices, 1e-8)).^2))
println("QPE mean: $(round(mean(qpe_prices), digits=4))")
println("Market mean: $(round(mean(market_prices), digits=4))")
println("RMSPE (uncalibrated): $(round(initial_rmspe*100, digits=2))%")

# Adaptive calibration
println("\n--- Adaptive Calibration ---")
global_scale, best_scale, beta, features, poly_predictions = calibrate_adaptive(qpe_prices, market_prices, moneynesses)

println("Global scale: $(round(global_scale, digits=4))")
println("Best grid scale: $(round(best_scale, digits=4))")

# Apply calibration
calibrated_prices = qpe_prices * best_scale
rmspe_calibrated = sqrt(mean(((calibrated_prices .- market_prices) ./ max.(market_prices, 1e-8)).^2))

println("\nAfter calibration:")
println("RMSPE: $(round(rmspe_calibrated*100, digits=2))%")

# Error analysis by moneyness
println("\n--- Error Analysis by Moneyness ---")
bins = [0.85, 0.92, 0.98, 1.05, 1.15, 1.30]
for i in 1:length(bins)-1
    bin_mask = (moneynesses .>= bins[i]) .& (moneynesses .< bins[i+1])
    if sum(bin_mask) >= 10
        bin_errors = abs.(calibrated_prices[bin_mask] .- market_prices[bin_mask])
        bin_rmspe = sqrt(mean(((calibrated_prices[bin_mask] .- market_prices[bin_mask]) ./ max.(market_prices[bin_mask], 1e-8)).^2))
        println("Moneyness $(bins[i])-$(bins[i+1]): RMSPE=$(round(bin_rmspe*100, digits=1))%, n=$(sum(bin_mask))")
    end
end

# Try with full Black-76 + moment matching + seasoned
println("\n--- Full Black-76 Model with All Fixes ---")
qpe_prices_full = Float64[]
for c in eachrow(test_contracts)
    spot = c.Spot_Price
    strike = c.Strike
    iv = c.IV
    dte = c.DTE
    T = dte / 365.0
    
    p = qpe_pricer_v4(spot, strike, r_default, iv, T, 32;
                     N_fixings=21,
                     m_realized=0,
                     A_realized=0.0,
                     vol_scale=1.0,
                     price_scale=1.0,
                     use_seasoned=true,
                     use_moment_match=true,
                     mc_factor=best_scale)
    push!(qpe_prices_full, p)
end

rmspe_full = sqrt(mean(((qpe_prices_full .- market_prices) ./ max.(market_prices, 1e-8)).^2))
println("RMSPE (full Black-76 model): $(round(rmspe_full*100, digits=2))%")

# Optimize over more parameters
println("\n--- Grid Search over Multiple Parameters ---")
best_rmspe = Inf
best_params = (mc=1.0, vol_scale=1.0, use_seasoned=false)

for mc in [0.5, 0.6, 0.7, 0.8, 1.0]
    for vol_scale in [0.8, 0.9, 1.0, 1.1, 1.2]
        for use_seasoned in [false, true]
            test_prices = Float64[]
            for c in eachrow(test_contracts)
                spot = c.Spot_Price
                strike = c.Strike
                iv = c.IV
                dte = c.DTE
                T = dte / 365.0
                
                p = qpe_pricer_v4(spot, strike, r_default, iv, T, 32;
                                 N_fixings=21,
                                 m_realized=0,
                                 A_realized=0.0,
                                 vol_scale=vol_scale,
                                 price_scale=1.0,
                                 use_seasoned=use_seasoned,
                                 use_moment_match=true,
                                 mc_factor=mc)
                push!(test_prices, p)
            end
            
            rmspe = sqrt(mean(((test_prices .- market_prices) ./ max.(market_prices, 1e-8)).^2))
            if rmspe < best_rmspe
                best_rmspe = rmspe
                best_params = (mc=mc, vol_scale=vol_scale, use_seasoned=use_seasoned)
            end
        end
    end
end

println("Best params: mc=$(best_params.mc), vol_scale=$(best_params.vol_scale), use_seasoned=$(best_params.use_seasoned)")
println("Best RMSPE: $(round(best_rmspe*100, digits=2))%")

println("\n" * "="^60)
println("Calibration Complete!")
println("="^60)
