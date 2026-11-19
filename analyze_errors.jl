using Printf, Statistics, LinearAlgebra, Distributions, CSV, DataFrames

# Load and analyze errors
df = DataFrame(CSV.File("/home/skoranne/PROGRAMS/QUANTUM/CB_BZ_option_chain_2026-11-16.csv"))
mask = (df.DTE .> 0) .& (df.IV .> 0.01) .& (df.Mid .> 0)
df_filtered = df[mask, :]

n_test = 1000
test_contracts = df_filtered[1:n_test, :]

market_prices = Float64[c.Mid for c in eachrow(test_contracts)]
moneynesses = Float64[c.Strike / c.Spot_Price for c in eachrow(test_contracts)]
opt_types = Char[c.Type[1] for c in eachrow(test_contracts)]
r_default = 0.045

# Black-76 pricer for comparison
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

# Compute Black-Scholes prices
bs_prices = zeros(n_test)
for i in 1:n_test
    c = test_contracts[i, :]
    T = c.DTE / 365.0
    if c.Type[1] == 'C'
        bs_prices[i] = black76_call(c.Spot_Price, c.Strike, T, r_default, c.IV)
    else
        bs_prices[i] = black76_put(c.Spot_Price, c.Strike, T, r_default, c.IV)
    end
end

println("=== Error Analysis ===")
println("\nMarket vs BS Ratio (should be ~1.0 for efficient markets):")
println("  Mean: $(round(mean(market_prices ./ bs_prices), digits=4))")
println("  Median: $(round(median(market_prices ./ bs_prices), digits=4))")
println("  Std: $(round(std(market_prices ./ bs_prices), digits=4))")

# Histogram by moneyness
println("\n=== Market Price Distribution by Moneyness ===")
bins = [0.80, 0.90, 0.95, 1.0, 1.05, 1.10, 1.20]
for i in 1:length(bins)-1
    bin_mask = (moneynesses .>= bins[i]) .& (moneynesses .< bins[i+1])
    if sum(bin_mask) >= 5
        bin_prices = market_prices[bin_mask]
        println("Moneyness $(bins[i])-$(bins[i+1]): n=$(sum(bin_mask)), median=$(round(median(bin_prices), digits=4)), mean=$(round(mean(bin_prices), digits=4))")
    end
end

# Check if market prices follow lognormal-like distribution
println("\n=== Price Ratio Analysis (Market/BS) ===")
ratios = market_prices ./ max.(bs_prices, 1e-10)
println("Ratio mean: $(round(mean(ratios), digits=4))")
println("Ratio median: $(round(median(ratios), digits=4))")
println("Ratio std: $(round(std(ratios), digits=4))")
println("Ratio min: $(round(minimum(ratios), digits=4))")
println("Ratio max: $(round(maximum(ratios), digits=4))")

# By option type
println("\n=== By Option Type ===")
call_mask = opt_types .== 'C'
put_mask = opt_types .== 'P'

println("Calls: n=$(sum(call_mask))")
println("  BS mean: $(round(mean(bs_prices[call_mask]), digits=4))")
println("  Market mean: $(round(mean(market_prices[call_mask]), digits=4))")
println("  Market/BS mean: $(round(mean(market_prices[call_mask] ./ bs_prices[call_mask]), digits=4))")

println("Puts: n=$(sum(put_mask))")
println("  BS mean: $(round(mean(bs_prices[put_mask]), digits=4))")
println("  Market mean: $(round(mean(market_prices[put_mask]), digits=4))")
println("  Market/BS mean: $(round(mean(market_prices[put_mask] ./ bs_prices[put_mask]), digits=4))")
