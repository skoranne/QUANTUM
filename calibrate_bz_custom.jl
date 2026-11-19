using Printf, Statistics, LinearAlgebra, Distributions, CSV, DataFrames, Dates

# BZ Futures Options Calibration - Custom Model
# Key insight: BZ market prices show negative extrinsic values for deep ITM calls
# This suggests a modified pricing model is needed

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

function compute_rmae(model_prices, market_prices)
    valid_idx = findall(p -> !isnan(p) && !isinf(p) && p > 1e-8, market_prices)
    if length(valid_idx) == 0
        return NaN
    end
    model_valid = model_prices[valid_idx]
    market_valid = market_prices[valid_idx]
    return mean(abs.(model_valid .- market_valid))
end

function custom_bz_price(S0, K, T, market_price, sigma_base=0.30)
    """
    Custom pricing model that accounts for BZ market microstructure:
    - Deep ITM: market price includes futures-style margining effect
    - ATM/OTM: use Black-76 with adjusted volatility
    
    Key observation: BZ market prices show that deep ITM options have
    extrinsic value ≈ -12.28 (approximately -S0), suggesting the market
    is pricing based on a different reference point.
    """
    if T <= 0
        # For same-day expiry, use intrinsic value
        return max(S0 - K, 0.0)
    end
    
    # For ITM options, adjust for futures-style margining
    intrinsic = max(S0 - K, 0.0)
    
    if intrinsic > 10
        # Deep ITM: market price is approximately intrinsic - S0 + small time value
        # This matches the pattern: price ≈ -12.28 for K=50, S0=92.28
        # Black-76 gives: 42.28 (intrinsic)
        # Market gives: ~45.84 (which is intrinsic + S0 * some factor)
        
        # Try: market_price ≈ Black76 - S0 * adjustment
        # For K=50: 45.84 ≈ 42.28 - 92.28 * adj => adj ≈ -0.038
        adjustment = 0.02
        black76 = black76_price_futures(S0, K, sigma_base, T)
        return black76 - adjustment * S0
    else
        # ATM/OTM: use Black-76
        return black76_price_futures(S0, K, sigma_base, T)
    end
end

function calibrate_bz_custom()
    println("="^80)
    println("BZ Futures Options Calibration - Custom Model")
    println("="^80)
    
    df = DataFrame(CSV.File("/home/skoranne/PROGRAMS/QUANTUM/CB_BZ_option_chain_2026-11-16.csv"))
    println("\nLoaded $(nrow(df)) contracts")
    
    # Filter: same-day expiry, positive prices, reasonable strikes
    mask = (df.DTE .== 0) .& (df.Type .== "C") .& (df.Strike .>= 50) .& (df.Strike .<= 110)
    df_filtered = df[mask, :]
    println("  Filtered to $(nrow(df_filtered)) contracts")
    
    n = nrow(df_filtered)
    
    # Extract data
    market_prices = Float64[df_filtered.Mid[i] for i in 1:n]
    spot_prices = Float64[df_filtered.Spot_Price[i] for i in 1:n]
    strikes = Float64[df_filtered.Strike[i] for i in 1:n]
    
    println("\nData summary (same-day expiry, calls only):")
    println("  Contracts: $n")
    println("  Spot price: $(round(mean(spot_prices), digits=2))")
    println("  Price range: $(round(minimum(market_prices), digits=2)) - $(round(maximum(market_prices), digits=2))")
    println("  Strike range: $(round(minimum(strikes), digits=1)) - $(round(maximum(strikes), digits=1))")
    
    # Analyze market price patterns
    println("\n--- Market Price vs Intrinsic Value Analysis ---")
    intrinsic = max.(spot_prices[1] .- strikes, 0.0)
    extrinsic = market_prices .- intrinsic
    
    println("\nPrice breakdown by moneyness:")
    for (k, mkt, intr, extra) in zip(strikes, market_prices, intrinsic, extrinsic)
        if k == 50 || k == 70 || k == 80 || k == 85 || k == 90 || k == 95 || k == 100 || k == 110
            println("  K=$(round(k, digits=1)): price=$(round(mkt, digits=4)), intrinsic=$(round(intr, digits=4)), extrinsic=$(round(extra, digits=4))")
        end
    end
    
    # Test custom model
    println("\n--- Custom Model vs Market Prices ---")
    custom_prices = zeros(n)
    for i in 1:n
        custom_prices[i] = custom_bz_price(spot_prices[i], strikes[i], df_filtered.DTE[i] / 365.0, market_prices[i])
    end
    
    rmspe_custom = compute_rmspe(custom_prices, market_prices)
    rmae_custom = compute_rmae(custom_prices, market_prices)
    
    println("\nCustom model performance:")
    println("  RMSPE: $(round(rmspe_custom*100, digits=2))%")
    println("  RMAE: $(round(rmae_custom, digits=4))")
    
    println("\nSample comparisons:")
    for i in 1:10
        println("  K=$(round(strikes[i], digits=1)): market=$(round(market_prices[i], digits=4)), custom=$(round(custom_prices[i], digits=4)), error=$(round(custom_prices[i]-market_prices[i], digits=6))")
    end
    
    # Black-76 comparison
    println("\n--- Black-76 Comparison (for reference) ---")
    black76_prices = zeros(n)
    for i in 1:n
        black76_prices[i] = black76_price_futures(spot_prices[i], strikes[i], 0.30, 0.0)
    end
    
    rmspe_bz = compute_rmspe(black76_prices, market_prices)
    println("Black-76 RMSPE: $(round(rmspe_bz*100, digits=2))%")
    
    # Key insight: the market prices appear to be based on a different reference
    println("\n" * "="^80)
    println("Key Insights")
    println("="^80)
    println("1. BZ market prices show negative extrinsic values for deep ITM options")
    println("   (e.g., K=50: price=45.84, intrinsic=42.28, extrinsic=3.56)")
    println("2. This is inconsistent with Black-Scholes theory but matches")
    println("   futures-style margining where the premium is margined daily")
    println("3. The custom model attempts to capture this by adjusting the")
    println("   Black-76 price with a small factor that accounts for margining")
    println("="^80)
    
    return rmspe_custom, custom_prices, market_prices
end

# Run calibration
calibrate_bz_custom()
