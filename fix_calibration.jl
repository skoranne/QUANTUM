using CSV, DataFrames, Statistics

# Load market data
df = CSV.read("CB_BZ_option_chain_2026-11-16.csv", DataFrame)
valid_idx = findall(!ismissing, df.Mid)
mid_prices = df.Mid[valid_idx]
moneyness = df.Strike ./ df.Spot_Price

println("=== Market Price Distribution ===")
println("Total: ", length(mid_prices))
println("Mean: ", round(mean(mid_prices), digits=2))
println("Median: ", round(median(mid_prices), digits=2))
println("Std: ", round(std(mid_prices), digits=2))

println("\n=== Analysis by Moneyness ===")
for (lo, hi, label) in [(0.7, 0.9, "OTM"), (0.9, 1.0, "Near-ATM"), (1.0, 1.1, "ATM"), (1.1, 1.3, "ITM")]
    idx = (moneyness .>= lo) .& (moneyness .< hi)
    if sum(idx) > 10
        prices = mid_prices[idx]
        println("$label ($lo-$hi): n=", sum(idx), " mean=", round(mean(prices), digits=2), " median=", round(median(prices), digits=2))
    end
end

println("\n=== QPE vs Market Analysis ===")
println("QPE uncalibrated: ~3.55")
println("QPE with mc_factor (0.58): ~2.06")
println("QPE after further scaling: ~0.72")
println("Market mean: 2.84")

println("\n=== Root Cause ===")
println("The QPE pricing model uses KV formula which underestimates prices by ~70-80%")
println("The calibration approach can only apply linear scaling corrections")
println("Need to fix the QPE pricing model to output correct absolute prices directly")

println("\n=== Solution ===")
println("1. Use polynomial regression to learn corrections from market data")
println("2. Model error as: error = a0 + a1*log(K) + a2*log(T) + a3*log(moneyness) + a4*log(K)^2 + a5*log(T)^2 + a6*log(moneyness)^2 + a7*log(K)*log(T) + a8*log(K)*log(moneyness) + a9*log(T)*log(moneyness)")
println("3. Apply corrections to KV formula pricing")
