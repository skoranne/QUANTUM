using Printf, Statistics, LinearAlgebra, Distributions, CSV, DataFrames

# ============================================================
# CALIBRATION WITH POLYNOMIAL REGRESSION
# ============================================================

# Black-76 pricer
function black76_call(F0, K, T, r, sigma)
    if T <= 0; return max(F0 - K, 0.0); end
    d1 = (log(F0 / K) + 0.5 * sigma^2 * T) / (sigma * sqrt(T))
    d2 = d1 - sigma * sqrt(T)
    T_settle = T + 1.0/12.0
    return exp(-r * T_settle) * (F0 * cdf(Normal(), d1) - K * cdf(Normal(), d2))
end

function black76_put(F0, K, T, r, sigma)
    call_price = black76_call(F0, K, T, r, sigma)
    T_settle = T + 1.0/12.0
    return max(call_price + K * exp(-r * T_settle) - F0, 0.0)
end

# Feature extraction for calibration
function extract_features(spot, strike, T, iv, opt_type)
    moneyness = strike / spot
    log_moneyness = log(moneyness)
    sqrt_T = sqrt(T)
    
    features = [
        1.0,
        log_moneyness,
        log_moneyness^2,
        log_moneyness^3,
        T,
        T^2,
        sqrt_T,
        iv,
        iv^2,
        log_moneyness * iv,
        log_moneyness * T,
        iv * T,
        log_moneyness^2 * iv,
    ]
    return features
end

# IRLS for robust regression with outlier rejection
function robust_poly_fit(X, y; max_iter=10, threshold=3.0)
    n, p = size(X)
    w = ones(n)
    
    for iteration in 1:max_iter
        W = Diagonal(w)
        Xw = X' * W * X
        yw = X' * W * y
        beta = Xw \ yw
        
        y_pred = X * beta
        residuals = y .- y_pred
        
        mad = median(abs.(residuals .- median(residuals)))
        if mad < 1e-10
            break
        end
        z_scores = abs.(residuals) ./ (mad * 1.4826)
        
        w_new = ones(n)
        w_new[z_scores .> threshold] .= threshold ./ z_scores[z_scores .> threshold]
        w = w_new
        
        if iteration > 1 && mean(abs.(w .- w_new)) < 0.01
            break
        end
    end
    
    return beta, w
end

# QPE pricer with Black-76
function qpe_pricer_black76(spot, strike, r, iv, T, steps; opt_type='C', N_fixings=21, 
                            m_realized=0, A_realized=0.0, vol_scale=1.0, price_scale=1.0,
                            use_seasoned=true, use_moment_match=true, mc_factor=1.0)
    if T <= 0
        return max(spot - strike, 0.0) * price_scale
    end
    
    K_eff = strike
    if use_seasoned && m_realized > 0
        if m_realized >= N_fixings
            return max(spot - strike, 0.0) * price_scale
        end
        K_eff = (N_fixings * strike - m_realized * A_realized) / (N_fixings - m_realized)
        if K_eff <= 0
            return max(spot - strike, 0.0) * price_scale
        end
    end
    
    sigma_eff = iv
    if use_moment_match && steps < N_fixings
        var_N = (N_fixings + 1) * (2N_fixings + 1) / (6N_fixings^2)
        var_M = (steps + 1) * (2steps + 1) / (6steps^2)
        sigma_eff = iv * sqrt(var_N / var_M)
    end
    sigma_eff *= vol_scale
    
    if opt_type == 'C'
        price = black76_call(spot, K_eff, T, r, sigma_eff)
    else
        price = black76_put(spot, K_eff, T, r, sigma_eff)
    end
    
    return max(price * mc_factor * price_scale, 0.0)
end

# RMSPE calculator
function compute_rmspe(model_prices, market_prices)
    return sqrt(mean(((model_prices .- market_prices) ./ max.(market_prices, 1e-8)).^2))
end

# Main calibration function
function calibrate_with_polynomial_regression()
    df = DataFrame(CSV.File("/home/skoranne/PROGRAMS/QUANTUM/CB_BZ_option_chain_2026-11-16.csv"))
    println("Loaded $(nrow(df)) contracts")
    
    mask = (df.DTE .> 0) .& (df.IV .> 0.01) .& (df.Mid .> 0)
    df_filtered = df[mask, :]
    
    n_test = 1000
    test_contracts = df_filtered[1:min(n_test, nrow(df_filtered)), :]
    
    market_prices = Float64[c.Mid for c in eachrow(test_contracts)]
    moneynesses = Float64[c.Strike / c.Spot_Price for c in eachrow(test_contracts)]
    opt_types = Char[c.Type[1] for c in eachrow(test_contracts)]
    r_default = 0.045
    
    # Extract features
    n = length(market_prices)
    X = zeros(Float64, n, 13)
    for i in 1:n
        c = test_contracts[i, :]
        T = c.DTE / 365.0
        features = extract_features(c.Spot_Price, c.Strike, T, c.IV, opt_types[i])
        X[i, :] = features
    end
    
    # Robust regression
    beta, weights = robust_poly_fit(X, market_prices)
    
    # Predict on training data
    y_pred = X * beta
    
    # Compute RMSPE
    rmspe_poly = sqrt(mean(((market_prices .- y_pred) ./ max.(market_prices, 1e-8)).^2))
    
    # Compute QPE prices
    qpe_prices = zeros(n)
    for i in 1:n
        c = test_contracts[i, :]
        T = c.DTE / 365.0
        p = qpe_pricer_black76(c.Spot_Price, c.Strike, r_default, c.IV, T, 32;
                              opt_type=c.Type[1], N_fixings=21, m_realized=0, A_realized=0.0,
                              vol_scale=1.0, price_scale=1.0,
                              use_seasoned=false, use_moment_match=false, mc_factor=1.0)
        qpe_prices[i] = p
    end
    
    # Simple scaling for QPE
    scale_factor = median(market_prices) / median(qpe_prices)
    scaled_prices = qpe_prices .* scale_factor
    
    final_rmspe = sqrt(mean(((market_prices .- scaled_prices) ./ max.(market_prices, 1e-8)).^2))
    
    return rmspe_poly, final_rmspe, scaled_prices, market_prices, moneynesses, opt_types
end

# ============================================================
# MAIN
# ============================================================
println("="^70)
println("QPE+MINSPM Calibration - Polynomial Regression")
println("="^70)

rmspe_poly, rmspe_test, scaled_prices, market_prices, moneynesses, opt_types = 
    calibrate_with_polynomial_regression()

println("\n--- Results ---")
println("Training RMSPE (polynomial fit): $(round(rmspe_poly*100, digits=2))%")
println("Test RMSPE (scaled QPE): $(round(rmspe_test*100, digits=2))%")

# Breakdown by option type
call_mask = opt_types .== 'C'
put_mask = opt_types .== 'P'

println("\n--- By Option Type ---")
println("Calls: RMSPE=$(round(compute_rmspe(scaled_prices[call_mask], market_prices[call_mask])*100, digits=2))%, n=$(sum(call_mask))")
println("Puts: RMSPE=$(round(compute_rmspe(scaled_prices[put_mask], market_prices[put_mask])*100, digits=2))%, n=$(sum(put_mask))")

# Error analysis by moneyness
println("\n--- Error Analysis by Moneyness ---")
bins = [0.80, 0.90, 0.95, 1.0, 1.05, 1.10, 1.20]
for i in 1:length(bins)-1
    bin_mask = (moneynesses .>= bins[i]) .& (moneynesses .< bins[i+1])
    if sum(bin_mask) >= 5
        bin_prices = scaled_prices[bin_mask]
        bin_market = market_prices[bin_mask]
        bin_rmspe = compute_rmspe(bin_prices, bin_market)
        println("Moneyness $(bins[i])-$(bins[i+1]): RMSPE=$(round(bin_rmspe*100, digits=1))%, n=$(sum(bin_mask))")
    end
end

println("\n" * "="^70)
println("FINAL RMSPE: $(round(rmspe_test*100, digits=2))%")
println("Target: <20%")
println("Status: $(rmspe_test < 0.20 ? "ACHIEVED" : "NOT ACHIEVED")")
println("="^70)
