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
println("Error Pattern Analysis v4 - Focus on non-zero prices")
println("="^70)

df = DataFrame(CSV.File("/home/skoranne/PROGRAMS/QUANTUM/CB_BZ_option_chain_2026-11-16.csv"))
println("\nLoaded $(nrow(df)) contracts")

# Filter for positive prices only
mask = (df.DTE .> 0) .& (df.IV .> 0.01) .& (df.Mid .> 0.01) .& (df.Mid .< 50)
df_filtered = df[mask, :]

println("Filtered to $(nrow(df_filtered)) contracts with positive prices")

n_test = min(500, nrow(df_filtered))
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
valid_mask = market_prices .> 0.001
qpe_valid = qpe_prices[valid_mask]
market_valid = market_prices[valid_mask]
bs_valid = bs_prices[valid_mask]

qpe_market_ratio = qpe_valid ./ max.(market_valid, 1e-10)
bs_market_ratio = bs_valid ./ max.(market_valid, 1e-8)

println("\n--- Ratio Analysis (non-zero prices) ---")
println("QPE/Market - Mean: $(round(mean(qpe_market_ratio), digits=4)), Median: $(round(median(qpe_market_ratio), digits=4))")
println("BS/Market - Mean: $(round(mean(bs_market_ratio), digits=4)), Median: $(round(median(bs_market_ratio), digits=4))")

println("\n--- Absolute Price Statistics ---")
println("QPE: mean=$(round(mean(qpe_valid), digits=4)), median=$(round(median(qpe_valid), digits=4)), std=$(round(std(qpe_valid), digits=4))")
println("Market: mean=$(round(mean(market_valid), digits=4)), median=$(round(median(market_valid), digits=4)), std=$(round(std(market_valid), digits=4))")
println("BS: mean=$(round(mean(bs_valid), digits=4)), median=$(round(median(bs_valid), digits=4)), std=$(round(std(bs_valid), digits=4))")

# RMSPE calculation
qpe_rmspe = sqrt(mean(((qpe_valid .- market_valid) ./ max.(market_valid, 1e-8)).^2))
bs_rmspe = sqrt(mean(((bs_valid .- market_valid) ./ max.(market_valid, 1e-8)).^2))

println("\n--- RMSPE ---")
println("QPE RMSPE: $(round(qpe_rmspe*100, digits=2))%")
println("BS RMSPE: $(round(bs_rmspe*100, digits=2))%")

# Calibrate QPE
scale_factor = mean(market_valid) / mean(qpe_valid)
qpe_calibrated = qpe_valid * scale_factor
qpe_calibrated_rmspe = sqrt(mean(((qpe_calibrated .- market_valid) ./ max.(market_valid, 1e-8)).^2))
println("\n--- Calibrated QPE ---")
println("Scale factor: $(round(scale_factor, digits=4))")
println("Calibrated RMSPE: $(round(qpe_calibrated_rmspe*100, digits=2))%")

# Error by moneyness and type
println("\n--- Error by Moneyness Bin and Type ---")
bins = [0.80, 0.90, 0.95, 1.0, 1.05, 1.10, 1.20]
for i in 1:length(bins)-1
    bin_mask = (moneynesses[valid_mask] .>= bins[i]) .& (moneynesses[valid_mask] .< bins[i+1])
    if sum(bin_mask) >= 5
        call_bin_mask = bin_mask .& (opt_types[valid_mask] .== 'C')
        put_bin_mask = bin_mask .& (opt_types[valid_mask] .== 'P')
        
        if sum(call_bin_mask) >= 3
            call_qpe = qpe_calibrated[call_bin_mask]
            call_market = market_valid[call_bin_mask]
            call_error = sqrt(mean(((call_qpe .- call_market) ./ max.(call_market, 1e-8)).^2))
            println("$(bins[i])-$(bins[i+1]) Calls: RMSPE=$(round(call_error*100, digits=1))%, n=$(sum(call_bin_mask))")
        end
        
        if sum(put_bin_mask) >= 3
            put_qpe = qpe_calibrated[put_bin_mask]
            put_market = market_valid[put_bin_mask]
            put_error = sqrt(mean(((put_qpe .- put_market) ./ max.(put_market, 1e-8)).^2))
            println("$(bins[i])-$(bins[i+1]) Puts: RMSPE=$(round(put_error*100, digits=1))%, n=$(sum(put_bin_mask))")
        end
    end
end

println("\n--- QPE vs BS comparison ---")
qpe_bs_ratio = qpe_valid ./ max.(bs_valid, 1e-10)
println("QPE/BS - Mean: $(round(mean(qpe_bs_ratio), digits=4)), Median: $(round(median(qpe_bs_ratio), digits=4))")

println("\n" * "="^70)
