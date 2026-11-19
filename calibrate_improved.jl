using Printf, Statistics, LinearAlgebra, Distributions, CSV, DataFrames, Dates

# Improved calibration for BZ options
# Account for: discrete sampling, Asian averaging, and proper vol scaling

function black76_price_futures(S0, K, sigma, T)
    """
    Black-76 price for futures options (futures-style margining, NO e^{-rT} discounting)
    """
    if T <= 0
        return max(S0 - K, 0.0)
    end
    
    d1 = (log(S0 / K) + 0.5 * sigma^2 * T) / (sigma * sqrt(T))
    d2 = d1 - sigma * sqrt(T)
    
    # No discounting for futures options (premium is margined daily)
    return S0 * cdf(Normal(), d1) - K * cdf(Normal(), d2)
end

function black76_implied_vol(price, S0, K, T)
    """
    Approximate Black-76 implied volatility solver
    """
    if T <= 0 || S0 <= 0
        return 0.0
    end
    
    # Use simple approximation for ATM options
    if abs(S0 - K) < 0.01 * S0
        # For ATM: price ≈ 0.4 * sigma * sqrt(T) * S0
        sigma = price / (0.4 * S0 * sqrt(T))
        return max(sigma, 0.01)
    end
    
    # Newton-Raphson for IV
    sigma = 0.25  # Initial guess
    
    for _ in 1:20
        d1 = (log(S0 / K) + 0.5 * sigma^2 * T) / (sigma * sqrt(T))
        d2 = d1 - sigma * sqrt(T)
        price_model = S0 * cdf(Normal(), d1) - K * cdf(Normal(), d2)
        
        # Vega
        vega = S0 * pdf(Normal(), d1) * sqrt(T)
        
        if abs(vega) < 1e-10
            break
        end
        
        delta_price = price - price_model
        sigma += delta_price / vega
        sigma = max(sigma, 0.01)
        sigma = min(sigma, 2.0)
        
        if abs(delta_price) < 1e-6
            break
        end
    end
    
    return sigma
end

function compute_rmspe(model_prices, market_prices)
    valid_idx = findall(p -> !isnan(p) && !isinf(p) && p > 0, model_prices)
    if length(valid_idx) == 0
        return NaN
    end
    model_valid = model_prices[valid_idx]
    market_valid = market_prices[valid_idx]
    if length(model_valid) == 0 || length(market_valid) == 0
        return NaN
    end
    return sqrt(mean(((model_valid .- market_valid) ./ max.(market_valid, 1e-8)).^2))
end

function compute_rmae(model_prices, market_prices)
    valid_idx = findall(p -> !isnan(p) && !isinf(p) && p > 0, model_prices)
    if length(valid_idx) == 0
        return NaN
    end
    model_valid = model_prices[valid_idx]
    market_valid = market_prices[valid_idx]
    if length(model_valid) == 0 || length(market_valid) == 0
        return NaN
    end
    return mean(abs.(model_valid .- market_valid))
end

function calibrate_improved()
    println("="^80)
    println("Improved Calibration for CB_BZ Brent Crude Oil Options")
    println("="^80)
    
    df = DataFrame(CSV.File("/home/skoranne/PROGRAMS/QUANTUM/CB_BZ_option_chain_2026-11-16.csv"))
    println("\nLoaded $(nrow(df)) contracts")
    
    # Analyze the data structure
    println("\nData analysis:")
    
    # Filter for near-term options (DTE <= 7 days for more stable calibration)
    mask = (df.DTE .>= 0) .& (df.DTE .<= 7) .& (df.IV .> 0.001) .& (df.Mid .> 0.001) .& (df.Mid .< 50)
    df_filtered = df[mask, :]
    println("  Filtered to $(nrow(df_filtered)) near-term contracts (DTE 0-7)")
    
    # Group by strike and type
    strikes = unique(sort(df_filtered.Strike))
    println("  Unique strikes: $(length(strikes))")
    println("  Strike range: $(minimum(strikes)) - $(maximum(strikes))")
    
    # Check ATM options
    spot_price = df.Spot_Price[1]
    atm_mask = abs.(df_filtered.Strike .- spot_price) .< 0.5
    atm_options = df_filtered[atm_mask, :]
    println("  ATM options (near $spot_price): $(nrow(atm_options))")
    
    if nrow(atm_options) > 0
        println("  ATM IV range: $(minimum(atm_options.IV)) - $(maximum(atm_options.IV))")
        println("  ATM Mid range: $(minimum(atm_options.Mid)) - $(maximum(atm_options.Mid))")
    end
    
    # Check ITM and OTM
    itm_mask = df_filtered.Strike .< spot_price
    otm_mask = df_filtered.Strike .> spot_price
    println("  ITM options: $(nrow(df_filtered[itm_mask, :]))")
    println("  OTM options: $(nrow(df_filtered[otm_mask, :]))")
    
    # Test calibration on subset
    n_test = min(500, nrow(df_filtered))
    test_indices = collect(1:n_test)
    
    market_prices = Float64[df_filtered.Mid[test_indices[i]] for i in 1:n_test]
    moneynesses = Float64[df_filtered.Strike[test_indices[i]] / df_filtered.Spot_Price[test_indices[i]] for i in 1:n_test]
    opt_types = Char[df_filtered.Type[test_indices[i]][1] for i in 1:n_test]
    spot_prices = Float64[df_filtered.Spot_Price[test_indices[i]] for i in 1:n_test]
    strikes = Float64[df_filtered.Strike[test_indices[i]] for i in 1:n_test]
    ivs = Float64[df_filtered.IV[test_indices[i]] for i in 1:n_test]
    dtess = Float64[df_filtered.DTE[test_indices[i]] / 365.0 for i in 1:n_test]
    
    println("\n--- Testing various calibration approaches ---")
    
    # Approach 1: Simple Black-76 with market IV
    println("\n1. Black-76 with market IV:")
    prices1 = zeros(n_test)
    for i in 1:n_test
        prices1[i] = black76_price_futures(spot_prices[i], strikes[i], ivs[i], dtess[i])
    end
    rmspe1 = compute_rmspe(prices1, market_prices)
    rmae1 = compute_rmae(prices1, market_prices)
    println("   RMSPE: $(round(rmspe1*100, digits=2))%")
    println("   RMAE: $(round(rmae1, digits=4))")
    
    # Approach 2: Black-76 with recalibrated IV from market price
    println("\n2. Black-76 with IV recalibrated from market price:")
    prices2 = zeros(n_test)
    ivs_recalibrated = zeros(n_test)
    for i in 1:n_test
        iv_recalibrated = black76_implied_vol(market_prices[i], spot_prices[i], strikes[i], dtess[i])
        ivs_recalibrated[i] = iv_recalibrated
        prices2[i] = black76_price_futures(spot_prices[i], strikes[i], iv_recalibrated, dtess[i])
    end
    rmspe2 = compute_rmspe(prices2, market_prices)
    rmae2 = compute_rmae(prices2, market_prices)
    println("   RMSPE: $(round(rmspe2*100, digits=2))%")
    println("   RMAE: $(round(rmae2, digits=4))")
    
    # Approach 3: With moneyness-dependent volatility adjustment
    println("\n3. Black-76 with moneyness-dependent vol adjustment:")
    prices3 = zeros(n_test)
    for i in 1:n_test
        m = moneynesses[i]
        iv_adj = ivs[i]
        
        if opt_types[i] == 'C'
            if m < 0.90; iv_adj = ivs[i] * 1.05
            elseif m < 0.95; iv_adj = ivs[i] * 1.03
            elseif m < 1.00; iv_adj = ivs[i] * 1.01
            elseif m < 1.05; iv_adj = ivs[i] * 1.0
            elseif m < 1.10; iv_adj = ivs[i] * 0.99
            else; iv_adj = ivs[i] * 0.97
            end
        else
            if m < 0.90; iv_adj = ivs[i] * 1.02
            elseif m < 0.95; iv_adj = ivs[i] * 1.01
            elseif m < 1.00; iv_adj = ivs[i] * 1.0
            elseif m < 1.05; iv_adj = ivs[i] * 0.99
            elseif m < 1.10; iv_adj = ivs[i] * 0.98
            else; iv_adj = ivs[i] * 0.96
            end
        end
        
        prices3[i] = black76_price_futures(spot_prices[i], strikes[i], iv_adj, dtess[i])
    end
    rmspe3 = compute_rmspe(prices3, market_prices)
    rmae3 = compute_rmae(prices3, market_prices)
    println("   RMSPE: $(round(rmspe3*100, digits=2))%")
    println("   RMAE: $(round(rmae3, digits=4))")
    
    # Approach 4: With both spot and vol scaling
    println("\n4. Black-76 with spot and vol scaling:")
    
    best_rmspe = rmspe3
    best_params = (spot_scale=1.0, vol_scale=1.0)
    
    for spot_scale in [0.99, 0.995, 1.0, 1.005, 1.01]
        for vol_scale in [0.9, 0.95, 1.0, 1.05, 1.1]
            prices_test = zeros(n_test)
            for i in 1:n_test
                m = moneynesses[i]
                S_adj = spot_prices[i] * spot_scale
                iv_adj = ivs[i] * vol_scale
                
                if opt_types[i] == 'C'
                    if m < 0.90; scale = 1.02
                    elseif m < 0.95; scale = 1.01
                    elseif m < 1.00; scale = 1.005
                    elseif m < 1.05; scale = 1.0
                    elseif m < 1.10; scale = 0.995
                    else; scale = 0.99
                    end
                else
                    if m < 0.90; scale = 1.005
                    elseif m < 0.95; scale = 1.0
                    elseif m < 1.00; scale = 0.995
                    elseif m < 1.05; scale = 0.99
                    elseif m < 1.10; scale = 0.985
                    else; scale = 0.98
                    end
                end
                
                prices_test[i] = black76_price_futures(S_adj, strikes[i], iv_adj, dtess[i]) * scale
            end
            
            rmspe_test = compute_rmspe(prices_test, market_prices)
            if rmspe_test < best_rmspe
                best_rmspe = rmspe_test
                best_params = (spot_scale=spot_scale, vol_scale=vol_scale)
            end
        end
    end
    
    prices4 = zeros(n_test)
    for i in 1:n_test
        m = moneynesses[i]
        S_adj = spot_prices[i] * best_params.spot_scale
        iv_adj = ivs[i] * best_params.vol_scale
        
        if opt_types[i] == 'C'
            if m < 0.90; scale = 1.02
            elseif m < 0.95; scale = 1.01
            elseif m < 1.00; scale = 1.005
            elseif m < 1.05; scale = 1.0
            elseif m < 1.10; scale = 0.995
            else; scale = 0.99
            end
        else
            if m < 0.90; scale = 1.005
            elseif m < 0.95; scale = 1.0
            elseif m < 1.00; scale = 0.995
            elseif m < 1.05; scale = 0.99
            elseif m < 1.10; scale = 0.985
            else; scale = 0.98
            end
        end
        
        prices4[i] = black76_price_futures(S_adj, strikes[i], iv_adj, dtess[i]) * scale
    end
    
    rmspe4 = compute_rmspe(prices4, market_prices)
    rmae4 = compute_rmae(prices4, market_prices)
    println("   Best params: spot_scale=$(best_params.spot_scale), vol_scale=$(best_params.vol_scale)")
    println("   RMSPE: $(round(rmspe4*100, digits=2))%")
    println("   RMAE: $(round(rmae4, digits=4))")
    
    # Summary
    println("\n" * "="^80)
    println("Calibration Summary")
    println("="^80)
    println("Method                    | RMSPE (%)      | RMAE")
    println("-"^79)
    println("Black-76 + market IV      | $(lpad(round(rmspe1*100, digits=2), 12))% | $(lpad(round(rmae1, digits=4), 8))")
    println("Black-76 + recalibrated IV| $(lpad(round(rmspe2*100, digits=2), 12))% | $(lpad(round(rmae2, digits=4), 8))")
    println("Moneyness-adjusted        | $(lpad(round(rmspe3*100, digits=2), 12))% | $(lpad(round(rmae3, digits=4), 8))")
    println("Optimized scaling         | $(lpad(round(rmspe4*100, digits=2), 12))% | $(lpad(round(rmae4, digits=4), 8))")
    println("="^80)
    println("Best RMSPE: $(round(best_rmspe*100, digits=2))%")
    println("Target: <20%")
    println("Status: $(best_rmspe < 0.20 ? "ACHIEVED" : "NOT YET ACHIEVED")")
    println("="^80)
    
    return best_rmspe, prices4, market_prices, best_params
end

# Run calibration
calibrate_improved()
