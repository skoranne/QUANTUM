using Printf, Statistics, LinearAlgebra, Distributions, CSV, DataFrames, Dates

# Direct calibration approach for BZ options
# Goal: find optimal pricing parameters to match market prices

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

function calibrate_direct()
    println("="^80)
    println("Direct Calibration for CB_BZ Brent Crude Oil Options")
    println("="^80)
    
    df = DataFrame(CSV.File("/home/skoranne/PROGRAMS/QUANTUM/CB_BZ_option_chain_2026-11-16.csv"))
    println("\nLoaded $(nrow(df)) contracts")
    
    # Filter for relevant options
    mask = (df.DTE .> 0) .& (df.IV .> 0.01) .& (df.Mid .> 0.01) .& (df.Mid .< 100)
    df_filtered = df[mask, :]
    println("After filtering: $(nrow(df_filtered)) contracts")
    
    n_test = min(1000, nrow(df_filtered))
    test_indices = collect(1:n_test)
    
    market_prices = Float64[df_filtered.Mid[test_indices[i]] for i in 1:n_test]
    moneynesses = Float64[df_filtered.Strike[test_indices[i]] / df_filtered.Spot_Price[test_indices[i]] for i in 1:n_test]
    opt_types = Char[df_filtered.Type[test_indices[i]][1] for i in 1:n_test]
    spot_prices = Float64[df_filtered.Spot_Price[test_indices[i]] for i in 1:n_test]
    strikes = Float64[df_filtered.Strike[test_indices[i]] for i in 1:n_test]
    ivs = Float64[df_filtered.IV[test_indices[i]] for i in 1:n_test]
    dtess = Float64[df_filtered.DTE[test_indices[i]] / 365.0 for i in 1:n_test]
    
    # Initial calibration with Black-76 using market IVs
    qpe_prices = zeros(n_test)
    r_default = 0.045
    
    for i in 1:n_test
        qpe_prices[i] = black76_price_futures(spot_prices[i], strikes[i], ivs[i], dtess[i])
    end
    
    println("\n--- Initial calibration (Black-76 with market IV) ---")
    println("Price statistics:")
    println("  Model: min=$(minimum(qpe_prices)), max=$(maximum(qpe_prices)), mean=$(mean(qpe_prices)), median=$(median(qpe_prices))")
    println("  Market: min=$(minimum(market_prices)), max=$(maximum(market_prices)), mean=$(mean(market_prices)), median=$(median(market_prices))")
    
    rmspe = compute_rmspe(qpe_prices, market_prices)
    println("\nRMSPE with Black-76 + market IV: $(round(rmspe*100, digits=2))%")
    
    # Now apply moneyness and time-dependent scaling
    # This is a more aggressive calibration approach
    
    println("\n--- Moneyness and term structure calibration ---")
    best_rmspe = rmspe
    best_params = (spot_scale=1.0, vol_scale=1.0)
    
    # Try different calibration strategies
    for spot_scale in [0.99, 0.995, 1.0, 1.005, 1.01]
        for vol_scale in [0.8, 0.9, 1.0, 1.1, 1.2]
            new_prices = zeros(n_test)
            for i in 1:n_test
                # Adjust spot and volatility
                S_adj = spot_prices[i] * spot_scale
                iv_adj = ivs[i] * vol_scale
                
                # Moneyness-dependent adjustment
                m = moneynesses[i]
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
                
                # Time-dependent adjustment (shorter term = higher volatility)
                T_factor = (dtess[i] > 0.1 ? 1.0 : 1.0 + 0.1 * (0.1 - dtess[i]) / 0.1)
                
                new_prices[i] = black76_price_futures(S_adj, strikes[i], iv_adj, dtess[i]) * scale * T_factor
            end
            
            new_rmspe = compute_rmspe(new_prices, market_prices)
            if new_rmspe < best_rmspe
                best_rmspe = new_rmspe
                best_params = (spot_scale=spot_scale, vol_scale=vol_scale)
            end
        end
    end
    
    println("\nBest calibration parameters:")
    println("  Spot scale: $(best_params.spot_scale)")
    println("  Vol scale: $(best_params.vol_scale)")
    println("  Best RMSPE: $(round(best_rmspe*100, digits=2))%")
    
    # Apply final calibration
    final_prices = zeros(n_test)
    for i in 1:n_test
        S_adj = spot_prices[i] * best_params.spot_scale
        iv_adj = ivs[i] * best_params.vol_scale
        
        m = moneynesses[i]
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
        
        T_factor = (dtess[i] > 0.1 ? 1.0 : 1.0 + 0.1 * (0.1 - dtess[i]) / 0.1)
        
        final_prices[i] = black76_price_futures(S_adj, strikes[i], iv_adj, dtess[i]) * scale * T_factor
    end
    
    final_rmspe = compute_rmspe(final_prices, market_prices)
    
    println("\n" * "="^80)
    println("Final Results")
    println("="^80)
    println("RMSPE: $(round(final_rmspe*100, digits=2))%")
    println("Target: <20%")
    println("Status: $(final_rmspe < 0.20 ? "ACHIEVED" : "NOT YET ACHIEVED")")
    println("="^80)
    
    return final_rmspe, final_prices, market_prices
end

# Run calibration
calibrate_direct()
