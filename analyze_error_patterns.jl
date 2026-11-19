using Printf, Statistics, LinearAlgebra, Distributions, CSV, DataFrames

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

println("="^70)
println("Error Pattern Analysis")
println("="^70)

df = DataFrame(CSV.File("/home/skoranne/PROGRAMS/QUANTUM/CB_BZ_option_chain_2026-11-16.csv"))
println("\nLoaded $(nrow(df)) contracts")

mask = (df.DTE .> 0) .& (df.IV .> 0.01) .& (df.Mid .> 0)
df_filtered = df[mask, :]

println("Filtered to $(nrow(df_filtered)) meaningful trades")

n_test = min(1000, nrow(df_filtered))
test_indices = collect(1:n_test)

market_prices = Float64[df.Mid[test_indices[i]] for i in 1:n_test]
moneynesses = Float64[df.Strike[test_indices[i]] / df.Spot_Price[test_indices[i]] for i in 1:n_test]
opt_types = Char[df.Type[test_indices[i]][1] for i in 1:n_test]
r_default = 0.045

qpe_prices = zeros(n_test)
bs_prices = zeros(n_test)
for i in 1:n_test
    idx = test_indices[i]
    T = df.DTE[idx] / 365.0
    qpe_prices[i] = qpe_pricer_black76(df.Spot_Price[idx], df.Strike[idx], r_default, df.IV[idx], T, 32;
                                      opt_type=df.Type[idx][1], N_fixings=21, m_realized=0, A_realized=0.0,
                                      vol_scale=1.0, price_scale=1.0,
                                      use_seasoned=false, use_moment_match=false, mc_factor=1.0)
    if opt_types[i] == 'C'
        bs_prices[i] = black76_call(df.Spot_Price[idx], df.Strike[idx], T, r_default, df.IV[idx])
    else
        bs_prices[i] = black76_put(df.Spot_Price[idx], df.Strike[idx], T, r_default, df.IV[idx])
    end
end

# Analyze errors
qpe_market_ratio = qpe_prices ./ max.(market_prices, 1e-10)
bs_market_ratio = bs_prices ./ max.(market_prices, 1e-10)

println("\n--- Ratio Analysis ---")
println("QPE/Market - Mean: $(round(mean(qpe_market_ratio), digits=4)), Median: $(round(median(qpe_market_ratio), digits=4))")
println("BS/Market - Mean: $(round(mean(bs_market_ratio), digits=4)), Median: $(round(median(bs_market_ratio), digits=4))")

println("\n--- Absolute Price Statistics ---")
println("QPE: mean=$(round(mean(qpe_prices), digits=4)), median=$(round(median(qpe_prices), digits=4)), std=$(round(std(qpe_prices), digits=4))")
println("Market: mean=$(round(mean(market_prices), digits=4)), median=$(round(median(market_prices), digits=4)), std=$(round(std(market_prices), digits=4))")
println("BS: mean=$(round(mean(bs_prices), digits=4)), median=$(round(median(bs_prices), digits=4)), std=$(round(std(bs_prices), digits=4))")

# Error by moneyness and type
println("\n--- Error by Moneyness Bin and Type ---")
bins = [0.80, 0.90, 0.95, 1.0, 1.05, 1.10, 1.20]
for i in 1:length(bins)-1
    bin_mask = (moneynesses .>= bins[i]) .& (moneynesses .< bins[i+1])
    if sum(bin_mask) >= 5
        call_bin_mask = bin_mask .& (opt_types .== 'C')
        put_bin_mask = bin_mask .& (opt_types .== 'P')
        
        if sum(call_bin_mask) >= 3
            call_error = mean(abs.(qpe_prices[call_bin_mask] .- market_prices[call_bin_mask]) ./ max.(market_prices[call_bin_mask], 1e-8))
            println("$(bins[i])-$(bins[i+1]) Calls: error=$(round(call_error*100, digits=1))%, n=$(sum(call_bin_mask))")
        end
        
        if sum(put_bin_mask) >= 3
            put_error = mean(abs.(qpe_prices[put_bin_mask] .- market_prices[put_bin_mask]) ./ max.(market_prices[put_bin_mask], 1e-8))
            println("$(bins[i])-$(bins[i+1]) Puts: error=$(round(put_error*100, digits=1))%, n=$(sum(put_bin_mask))")
        end
    end
end

println("\n" * "="^70)
