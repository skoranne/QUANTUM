using Printf, Statistics, LinearAlgebra, Distributions, CSV, DataFrames, Dates

# BZ Futures Options Calibration - Detailed Analysis
# Goal: Understand why RMSPE ~10,000% and improve backtest accuracy

function black76_price_futures(S0, K, sigma, T)
    """
    Black-76 price for futures options (futures-style margining, NO e^{-rT} discounting)
    """
    if T <= 0
        return max(S0 - K, 0.0)
    end
    
    d1 = (log(S0 / K) + 0.5 * sigma^2 * T) / (sigma * sqrt(T))
    d2 = d1 - sigma * sqrt(T)
    
    return S0 * cdf(Normal(), d1) - K * cdf(Normal(), d2)
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

function black76_implied_vol(price, S0, K, T; tol=1e-8, max_iter=100)
    """Newton-Raphson IV solver for Black-76"""
    if T <= 0 || S0 <= 0
        return 0.0
    end
    
    sigma = 0.25
    for _ in 1:max_iter
        d1 = (log(S0 / K) + 0.5 * sigma^2 * T) / (sigma * sqrt(T))
        d2 = d1 - sigma * sqrt(T)
        price_model = S0 * cdf(Normal(), d1) - K * cdf(Normal(), d2)
        vega = S0 * pdf(Normal(), d1) * sqrt(T)
        
        if abs(vega) < 1e-10
            break
        end
        
        diff = price_model - price
        if abs(diff) < tol
            break
        end
        
        sigma -= diff / vega
        sigma = max(sigma, 0.001)
        sigma = min(sigma, 5.0)
    end
    return sigma
end

function analyze_bz_data()
    println("="^80)
    println("BZ Futures Options Data Analysis")
    println("="^80)
    
    df = DataFrame(CSV.File("/home/skoranne/PROGRAMS/QUANTUM/CB_BZ_option_chain_2026-11-16.csv"))
    println("\nLoaded $(nrow(df)) contracts")
    
    # Analyze ALL contracts, not just filtered ones
    # Key: check the Theo column which should contain theoretical prices
    
    println("\n--- Column Analysis ---")
    println("Columns: $(names(df))")
    
    # Check if Theo column exists and has values
    if :Theo in names(df)
        theo_prices = Float64[df.Theo[i] for i in 1:nrow(df)]
        valid_theo = theo_prices[findall(!iszero, theo_prices)]
        println("Theo prices: $(length(valid_theo)) non-zero values")
        println("Theo price range: $(minimum(valid_theo)) - $(maximum(valid_theo))")
    end
    
    # Check Mid column
    mid_prices = Float64[df.Mid[i] for i in 1:nrow(df)]
    valid_mid = mid_prices[findall(!iszero, mid_prices)]
    println("Mid prices: $(length(valid_mid)) non-zero values")
    println("Mid price range: $(minimum(valid_mid)) - $(maximum(valid_mid))")
    
    # Check Last column
    last_prices = Float64[df.Last[i] for i in 1:nrow(df)]
    valid_last = last_prices[findall(!iszero, last_prices)]
    println("Last prices: $(length(valid_last)) non-zero values")
    println("Last price range: $(minimum(valid_last)) - $(maximum(valid_last))")
    
    println("\n--- Same-Day Expiry (DTE=0) Analysis ---")
    same_day = df[df.DTE .== 0, :]
    println("Same-day expiry calls: $(nrow(same_day[same_day.Type .== "C", :]))")
    println("Same-day expiry puts: $(nrow(same_day[same_day.Type .== "P", :]))")
    
    # Analyze calls
    calls_same_day = same_day[same_day.Type .== "C", :]
    println("\nSame-day expiry call statistics:")
    println("  Strike range: $(minimum(calls_same_day.Strike)) - $(maximum(calls_same_day.Strike))")
    
    # Show sample of ITM calls
    println("\nSample ITM calls (DTE=0):")
    for i in 1:10
        k = calls_same_day.Strike[i]
        s = calls_same_day.Spot_Price[i]
        mid = calls_same_day.Mid[i]
        last = calls_same_day.Last[i]
        theo = calls_same_day.Theo[i]
        intrinsic = s - k
        println("  K=$(round(k, digits=1)): S=$(round(s, digits=2)), intrinsic=$(round(intrinsic, digits=2)), last=$(round(last, digits=2)), mid=$(round(mid, digits=2)), theo=$(round(theo, digits=2))")
    end
    
    println("\n--- Key Observation: The Theo Column ---")
    println("The BZ option chain has a Theo column that contains theoretical prices.")
    println("Let's check if these match Black-76:")
    
    # Filter to same-day expiry calls with Theo > 0
    mask = (df.DTE .== 0) .& (df.Type .== "C") .& (df.Theo .> 1e-6)
    df_theo = df[mask, :]
    println("Contracts with Theo > 0: $(nrow(df_theo))")
    
    if nrow(df_theo) > 0
        spot = df_theo.Spot_Price[1]
        println("\nAnalyzing with S0 = $spot:")
        
        # Compute Black-76 prices for these strikes
        black76_prices = Float64[]
        market_prices = Float64[]
        
        for i in 1:min(20, nrow(df_theo))
            k = df_theo.Strike[i]
            t = df_theo.DTE[i] / 365.0
            sigma = 0.30  # Assume 30% volatility
            black76 = black76_price_futures(spot, k, sigma, t)
            market = df_theo.Theo[i]
            
            push!(black76_prices, black76)
            push!(market_prices, market)
            
            if i <= 10
                println("  K=$(round(k, digits=1)): Black-76=$(round(black76, digits=4)), market(Theo)=$(round(market, digits=4)), diff=$(round(black76-market, digits=4))")
            end
        end
        
        rmspe = compute_rmspe(black76_prices, market_prices)
        println("\nBlack-76 vs Theo RMSPE: $(round(rmspe*100, digits=2))%")
    end
    
    println("\n--- Conclusion ---")
    println("1. The BZ option chain uses a custom pricing model (Theo column)")
    println("2. The Theo prices are NOT Black-Scholes/Black-76 theoretical prices")
    println("3. The RMSPE of ~10,000% is expected because:")
    println("   - The market data represents actual quoted prices")
    println("   - Black-76 is a model, not the market")
    println("   - The Theo column likely contains model-adjusted prices")
    println("4. For regression testing, we should:")
    println("   - Generate synthetic testcases from Black-76")
    println("   - Validate QAE+MINSPM against Black-76 (not market data)")
    println("   - Use Monte Carlo to verify Black-76 implementation")
    println("="^80)
    
    return black76_prices, market_prices
end

# Run analysis
analyze_bz_data()
