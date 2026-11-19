using Printf, Statistics, LinearAlgebra, Distributions, CSV, DataFrames, Dates
using Optim

# Improved calibration for BZ options
# Focus: futures-style margining (NO e^{-rT} discounting), proper S0 representation

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

function black76_vega(S0, K, sigma, T)
    """Compute vega for Black-76"""
    if T <= 0 || sigma <= 0
        return 0.0
    end
    
    d1 = (log(S0 / K) + 0.5 * sigma^2 * T) / (sigma * sqrt(T))
    return S0 * pdf(Normal(), d1) * sqrt(T)
end

function black76_implied_vol(price, S0, K, T; tol=1e-8, max_iter=100)
    """Newton-Raphson IV solver"""
    if T <= 0 || S0 <= 0
        return 0.0
    end
    
    sigma = 0.25  # Initial guess
    
    for _ in 1:max_iter
        d1 = (log(S0 / K) + 0.5 * sigma^2 * T) / (sigma * sqrt(T))
        d2 = d1 - sigma * sqrt(T)
        price_model = S0 * cdf(Normal(), d1) - K * cdf(Normal(), d2)
        vega = black76_vega(S0, K, sigma, T)
        
        if abs(vega) < 1e-10
            break
        end
        
        diff = price_model - price
        if abs(diff) < tol
            break
        end
        
        delta_sigma = diff / vega
        sigma -= delta_sigma
        sigma = max(sigma, 0.001)
        sigma = min(sigma, 5.0)
    end
    
    return sigma
end

function compute_rmspe(model_prices, market_prices)
    valid_idx = findall(p -> !isnan(p) && !isinf(p) && p > 1e-8, market_prices)
    if length(valid_idx) == 0
        return NaN
    end
    model_valid = model_prices[valid_idx]
    market_valid = market_prices[valid_idx]
    return sqrt(mean(((model_valid .- market_valid) ./ max.(abs.(market_valid), 1e-8)).^2))
end

function compute_rmae(model_prices, market_prices)
    valid_idx = findall(p -> !isnan(p) && !isinf(p) && p > 1e-8, market_prices)
    if length(valid_idx) == 0
        return NaN
    end
    model_valid = model_prices[valid_idx]
    market_valid = market_prices[valid_idx]
    return mean(abs.(model_valid .- market_valid))
end

function calibrate_bz_futures_style()
    println("="^80)
    println("BZ Futures-Style Options Calibration")
    println("="^80)
    
    df = DataFrame(CSV.File("/home/skoranne/PROGRAMS/QUANTUM/CB_BZ_option_chain_2026-11-16.csv"))
    println("\nLoaded $(nrow(df)) contracts")
    
    # Filter: near-term options, positive prices, reasonable IV
    mask = (df.DTE .>= 0) .& (df.DTE .<= 7) .& (df.IV .> 0.001) .& (df.IV .< 5.0) .& 
           (df.Mid .> 1e-6) .& (df.Mid .< 100) .& (df.Spot_Price .> 0)
    df_filtered = df[mask, :]
    println("  Filtered to $(nrow(df_filtered)) contracts")
    
    # Extract data
    n = nrow(df_filtered)
    market_prices = Float64[df_filtered.Mid[i] for i in 1:n]
    spot_prices = Float64[df_filtered.Spot_Price[i] for i in 1:n]
    strikes = Float64[df_filtered.Strike[i] for i in 1:n]
    ivs = Float64[df_filtered.IV[i] for i in 1:n]
    dtess = Float64[df_filtered.DTE[i] / 365.0 for i in 1:n]
    opt_types = Char[df_filtered.Type[i][1] for i in 1:n]
    
    # Key insight: S0 varies by contract but should be constant
    # Use the futures price (Spot_Price) as S0 for each contract
    println("\n  Spot price range: $(minimum(spot_prices)) - $(maximum(spot_prices))")
    println("  Price range: $(minimum(market_prices)) - $(maximum(market_prices))")
    
    # Approach: Calibrate to reproduce market prices using Black-76
    # Model: price = f(S0, K, sigma, T) where S0 is futures price
    
    println("\n--- Calibration Approaches ---\n")
    
    # Approach 1: Direct Black-76 using market IV
    println("1. Direct Black-76 with market IV:")
    prices1 = zeros(n)
    for i in 1:n
        prices1[i] = black76_price_futures(spot_prices[i], strikes[i], ivs[i], dtess[i])
    end
    rmspe1 = compute_rmspe(prices1, market_prices)
    rmae1 = compute_rmae(prices1, market_prices)
    println("   RMSPE: $(round(rmspe1*100, digits=2))%")
    println("   RMAE: $(round(rmae1, digits=4))")
    
    # Approach 2: Recalibrate IV from market price
    println("\n2. IV recalibration from market price:")
    prices2 = zeros(n)
    ivs_recalibrated = zeros(n)
    for i in 1:n
        iv_rec = black76_implied_vol(market_prices[i], spot_prices[i], strikes[i], dtess[i])
        ivs_recalibrated[i] = iv_rec
        prices2[i] = black76_price_futures(spot_prices[i], strikes[i], iv_rec, dtess[i])
    end
    rmspe2 = compute_rmspe(prices2, market_prices)
    rmae2 = compute_rmae(prices2, market_prices)
    println("   RMSPE: $(round(rmspe2*100, digits=2))%")
    println("   RMAE: $(round(rmae2, digits=4))")
    
    # Approach 3: Use constant spot price (average)
    println("\n3. Black-76 with average spot price:")
    S0_avg = mean(spot_prices)
    prices3 = zeros(n)
    for i in 1:n
        prices3[i] = black76_price_futures(S0_avg, strikes[i], ivs[i], dtess[i])
    end
    rmspe3 = compute_rmspe(prices3, market_prices)
    rmae3 = compute_rmae(prices3, market_prices)
    println("   S0_avg: $(round(S0_avg, digits=2))")
    println("   RMSPE: $(round(rmspe3*100, digits=2))%")
    println("   RMAE: $(round(rmae3, digits=4))")
    
    # Approach 4: Use forward price ( Spot * e^(r*T) approximation)
    println("\n4. Black-76 with forward price adjustment:")
    r = 0.02  # Assumed risk-free rate
    prices4 = zeros(n)
    for i in 1:n
        F = spot_prices[i] * exp(r * dtess[i])  # Forward price
        prices4[i] = black76_price_futures(F, strikes[i], ivs[i], dtess[i])
    end
    rmspe4 = compute_rmspe(prices4, market_prices)
    rmae4 = compute_rmae(prices4, market_prices)
    println("   RMSPE: $(round(rmspe4*100, digits=2))%")
    println("   RMAE: $(round(rmae4, digits=4))")
    
    # Approach 5: Moneyness-dependent scaling
    println("\n5. Moneyness-dependent Black-76:")
    prices5 = zeros(n)
    for i in 1:n
        m = strikes[i] / spot_prices[i]
        S0_adj = spot_prices[i]
        
        # Adjust for moneyness
        if opt_types[i] == 'C'
            if m < 0.85; S0_adj = spot_prices[i] * 1.005
            elseif m < 0.90; S0_adj = spot_prices[i] * 1.002
            elseif m < 0.95; S0_adj = spot_prices[i] * 1.001
            elseif m > 1.15; S0_adj = spot_prices[i] * 0.995
            elseif m > 1.10; S0_adj = spot_prices[i] * 0.998
            end
        else
            if m < 0.85; S0_adj = spot_prices[i] * 1.003
            elseif m < 0.90; S0_adj = spot_prices[i] * 1.001
            elseif m > 1.10; S0_adj = spot_prices[i] * 0.997
            elseif m > 1.15; S0_adj = spot_prices[i] * 0.995
            end
        end
        
        prices5[i] = black76_price_futures(S0_adj, strikes[i], ivs[i], dtess[i])
    end
    rmspe5 = compute_rmspe(prices5, market_prices)
    rmae5 = compute_rmae(prices5, market_prices)
    println("   RMSPE: $(round(rmspe5*100, digits=2))%")
    println("   RMAE: $(round(rmae5, digits=4))")
    
    # Approach 6: Full optimization with multiple parameters
    println("\n6. Optimized Black-76 with calibration:")
    
    function objective(params)
        S0_mult, vol_mult, moneyness_weight = params
        S0_mult = max(S0_mult, 0.9)
        S0_mult = min(S0_mult, 1.1)
        vol_mult = max(vol_mult, 0.5)
        vol_mult = min(vol_mult, 2.0)
        moneyness_weight = max(moneyness_weight, 0.0)
        moneyness_weight = min(moneyness_weight, 1.0)
        
        prices_test = zeros(n)
        for i in 1:n
            m = strikes[i] / spot_prices[i]
            
            # Adjust S0
            S0_adj = spot_prices[i] * S0_mult
            
            # Moneyness adjustment
            if opt_types[i] == 'C'
                if m < 0.90; S0_adj *= (1.0 + moneyness_weight * 0.005)
                elseif m > 1.10; S0_adj *= (1.0 - moneyness_weight * 0.005)
                end
            else
                if m < 0.90; S0_adj *= (1.0 + moneyness_weight * 0.003)
                elseif m > 1.10; S0_adj *= (1.0 - moneyness_weight * 0.003)
                end
            end
            
            # Adjust volatility
            iv_adj = ivs[i] * vol_mult
            
            prices_test[i] = black76_price_futures(S0_adj, strikes[i], iv_adj, dtess[i])
        end
        
        return compute_rmspe(prices_test, market_prices)
    end
    
    # Grid search for initial parameters
    best_rmspe = Inf
    best_params = (1.0, 1.0, 0.0)
    
    println("   Grid search...")
    for S0_mult in [0.98, 0.99, 1.0, 1.01, 1.02]
        for vol_mult in [0.95, 1.0, 1.05]
            for mweight in [0.0, 0.5, 1.0]
                rmspe_test = objective((S0_mult, vol_mult, mweight))
                if rmspe_test < best_rmspe
                    best_rmspe = rmspe_test
                    best_params = (S0_mult, vol_mult, mweight)
                end
            end
        end
    end
    println("   Best from grid: S0_mult=$(best_params[1]), vol_mult=$(best_params[2]), mweight=$(best_params[3])")
    println("   RMSPE: $(round(best_rmspe*100, digits=2))%")
    
    # Refine with optimization
    println("   Refining with Nelder-Mead...")
    result = optimize(objective, [best_params[1], best_params[2], best_params[3]], 
                      NelderMead(), 
                      Optim.Options(show_trace=true, iterations=100))
    
    if !Optim.converged(result)
        println("   Optimization did not converge, using grid result")
    end
    
    final_params = Optim.minimizer(result)
    prices6 = zeros(n)
    for i in 1:n
        m = strikes[i] / spot_prices[i]
        S0_adj = spot_prices[i] * final_params[1]
        
        if opt_types[i] == 'C'
            if m < 0.90; S0_adj *= (1.0 + final_params[3] * 0.005)
            elseif m > 1.10; S0_adj *= (1.0 - final_params[3] * 0.005)
            end
        else
            if m < 0.90; S0_adj *= (1.0 + final_params[3] * 0.003)
            elseif m > 1.10; S0_adj *= (1.0 - final_params[3] * 0.003)
            end
        end
        
        iv_adj = ivs[i] * final_params[2]
        prices6[i] = black76_price_futures(S0_adj, strikes[i], iv_adj, dtess[i])
    end
    
    rmspe6 = compute_rmspe(prices6, market_prices)
    rmae6 = compute_rmae(prices6, market_prices)
    println("   Final RMSPE: $(round(rmspe6*100, digits=2))%")
    println("   Final RMAE: $(round(rmae6, digits=4))")
    
    # Summary
    println("\n" * "="^80)
    println("Calibration Summary")
    println("="^80)
    println("Method                         | RMSPE (%)   | RMAE")
    println("-"^79)
    println("Direct (market IV)             | $(lpad(round(rmspe1*100, digits=2), 11))% | $(lpad(round(rmae1, digits=4), 8))")
    println("IV recalibration               | $(lpad(round(rmspe2*100, digits=2), 11))% | $(lpad(round(rmae2, digits=4), 8))")
    println("Average spot price             | $(lpad(round(rmspe3*100, digits=2), 11))% | $(lpad(round(rmae3, digits=4), 8))")
    println("Forward price adjustment       | $(lpad(round(rmspe4*100, digits=2), 11))% | $(lpad(round(rmae4, digits=4), 8))")
    println("Moneyness adjustment           | $(lpad(round(rmspe5*100, digits=2), 11))% | $(lpad(round(rmae5, digits=4), 8))")
    println("Optimized (grid + Nelder-Mead) | $(lpad(round(rmspe6*100, digits=2), 11))% | $(lpad(round(rmae6, digits=4), 8))")
    println("="^80)
    println("Best RMSPE: $(round(rmspe6*100, digits=2))%")
    println("Target: <20%")
    println("Status: $(rmspe6 < 0.20 ? "ACHIEVED" : "NOT YET ACHIEVED")")
    println("="^80)
    
    return rmspe6, prices6, market_prices, final_params
end

# Run calibration
calibrate_bz_futures_style()
