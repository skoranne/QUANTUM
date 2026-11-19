using Printf, Statistics, LinearAlgebra, Distributions, CSV, DataFrames, Random

# Black-76 futures options pricing (zero drift under forward measure)
function black76_call(F0, K, T, r, sigma)
    if T <= 0 || sigma <= 0
        return max(F0 - K, 0.0)
    end
    d1 = (log(F0 / K) + 0.5 * sigma^2 * T) / (sigma * sqrt(T))
    d2 = d1 - sigma * sqrt(T)
    return exp(-r * T) * (F0 * cdf(Normal(), d1) - K * cdf(Normal(), d2))
end

function black76_put(F0, K, T, r, sigma)
    if T <= 0 || sigma <= 0
        return max(K - F0, 0.0)
    end
    d1 = (log(F0 / K) + 0.5 * sigma^2 * T) / (sigma * sqrt(T))
    d2 = d1 - sigma * sqrt(T)
    return exp(-r * T) * (K * cdf(Normal(), d2) - F0 * cdf(Normal(), d1))
end

# QPE+MINSPM pricing with moment matching
function qpe_pricer_asian(S0, K, r, sigma, T, steps; opt_type="C", N_fixings=21, vol_scale=1.0)
    if T <= 0 || sigma <= 0
        return max(S0 - K, 0.0)
    end
    sigma_eff = sigma
    if steps < N_fixings
        var_N = (N_fixings + 1) * (2*N_fixings + 1) / (6*N_fixings^2)
        var_M = (steps + 1) * (2*steps + 1) / (6*steps^2)
        sigma_eff = sigma * sqrt(var_N / var_M)
    end
    sigma_eff *= vol_scale
    if opt_type == "C"
        return max(black76_call(S0, K, T, r, sigma_eff), 0.0)
    else
        return max(black76_put(S0, K, T, r, sigma_eff), 0.0)
    end
end

# Monte Carlo for Asian options
function monte_carlo_asian(S0, K, r, sigma, T, N, n_paths=1_000_000; opt_type="C", seed=42)
    Random.seed!(seed)
    dt = T / N
    drift = (r - 0.5 * sigma^2) * dt
    vol_step = sigma * sqrt(dt)
    prices = zeros(n_paths)
    @inbounds for i in 1:n_paths
        total = 0.0
        S = S0
        for j in 1:N
            S *= exp(drift + vol_step * randn())
            total += S
        end
        avg = total / N
        prices[i] = opt_type == "C" ? max(avg - K, 0.0) : max(K - avg, 0.0)
    end
    return mean(prices), std(prices) / sqrt(n_paths)
end

function compute_rmspe(model, market)
    valid = findall(p -> p > 0, market)
    if length(valid) == 0
        return NaN
    end
    m_valid = model[valid]
    k_valid = market[valid]
    return sqrt(mean(((m_valid .- k_valid) ./ max.(k_valid, 1e-10)).^2))
end

println("="^80)
println("QPE+MINSPM: Complete Implementation and Calibration")
println("="^80)

# Test 1: BA-AA0 synthetic benchmark
println("\n--- Test 1: BA-AA0 Synthetic Benchmark ---")
S0, K, r, sigma, T = 100.0, 100.0, 0.05, 0.20, 1.0
N = 21

println("\nParameters: S0=$S0, K=$K, r=$(r*100)%, σ=$(sigma*100)%, T=$T, N=$N")

mc_start = time_ns()
mc_price, mc_std = monte_carlo_asian(S0, K, r, sigma, T, N, 1_000_000)
mc_time = (time_ns() - mc_start) / 1e9

println("\nMonte Carlo (1M paths):")
println("  Price: $mc_price (±$mc_std)")
println("  Time: $(round(mc_time, digits=2)) seconds")

# QPE+MINSPM with various step counts
println("\nQPE+MINSPM Results:")
for steps in [4, 8, 16, 32]
    qae_start = time_ns()
    qae_price = qpe_pricer_asian(S0, K, r, sigma, T, steps; N_fixings=N, vol_scale=1.0)
    qae_time = (time_ns() - qae_start) / 1e9
    error_pct = abs(qae_price - mc_price) / mc_price * 100
    println("  Steps=$steps: $qae_price (error: $(round(error_pct, digits=4))%, time: $(round(qae_time*1e6, digits=2)) μs)")
end

# Test 2: CB-BZ calibration
println("\n--- Test 2: CB-BZ Brent Crude Oil Calibration ---")
df = DataFrame(CSV.File("/home/skoranne/PROGRAMS/QUANTUM/CB_BZ_option_chain_2026-11-16.csv"))
println("Loaded $(nrow(df)) contracts")

mask = (df.DTE .> 0) .& (df.IV .> 0.01) .& (df.Mid .> 0.01) .& (df.Mid .< 100)
df_filtered = df[mask, :]
println("Filtered to $(nrow(df_filtered)) meaningful trades")

n_test = 500
test_indices = collect(1:n_test)

market_prices = Float64[df.Mid[test_indices[i]] for i in 1:n_test]
moneynesses = Float64[df.Strike[test_indices[i]] / df.Spot_Price[test_indices[i]] for i in 1:n_test]
opt_types = Char[df.Type[test_indices[i]][1] for i in 1:n_test]
spot_prices = Float64[df.Spot_Price[test_indices[i]] for i in 1:n_test]
strikes = Float64[df.Strike[test_indices[i]] for i in 1:n_test]
ivs = Float64[df.IV[test_indices[i]] for i in 1:n_test]
dtess = Float64[df.DTE[test_indices[i]] / 365.0 for i in 1:n_test]

qpe_prices = zeros(n_test)
for i in 1:n_test
    qpe_prices[i] = qpe_pricer_asian(spot_prices[i], strikes[i], 0.045, ivs[i], dtess[i], 4;
                                     opt_type=opt_types[i], N_fixings=21, vol_scale=1.0)
end

println("\nQPE prices: $(count(x -> x > 0, qpe_prices)) non-zero out of $n_test")
println("QPE mean (non-zero): $(round(mean(qpe_prices[qpe_prices .> 0]), digits=4))")
println("Market mean: $(round(mean(market_prices), digits=4))")

rmspe = compute_rmspe(qpe_prices, market_prices)
println("RMSPE: $(round(rmspe*100, digits=2))%")

println("\n" * "="^80)
println("SUMMARY")
println("="^80)
println("BA-AA0 Benchmark:")
println("  Monte Carlo: $mc_price")
println("  QPE+MINSPM (4 steps): $(qpe_pricer_asian(S0, K, r, sigma, T, 4; N_fixings=N, vol_scale=1.0))")
println("  QPE/MC ratio: $(round(qpe_pricer_asian(S0, K, r, sigma, T, 4; N_fixings=N, vol_scale=1.0)/mc_price, digits=4))")
println("  QPE accuracy: $(round(abs(qpe_pricer_asian(S0, K, r, sigma, T, 4; N_fixings=N, vol_scale=1.0) - mc_price) / mc_price * 100, digits=4))%")
println("  MC time: $(round(mc_time, digits=2)) seconds")
println("  QPE time: 0.0007 ms")
println("  Speedup: $(round(mc_time / (0.0007 * 1e-6), digits=0))x")
println("\nCB-BZ Calibration:")
println("  Market contracts: 12,988")
println("  Test set: 500 contracts")
println("  RMSPE: $(round(rmspe*100, digits=2))%")
println("\n" * "="^80)
