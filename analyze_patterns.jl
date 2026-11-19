using Printf, Statistics, LinearAlgebra, Distributions, CSV, DataFrames

# ============================================================
# ANALYZE ERROR PATTERNS
# ============================================================

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

n_test = min(2000, nrow(df_filtered))
test_contracts = df_filtered[1:n_test, :]

market_prices = Float64[c.Mid for c in eachrow(test_contracts)]
moneynesses = Float64[c.Strike / c.Spot_Price for c in eachrow(test_contracts)]
opt_types = Char[c.Type[1] for c in eachrow(test_contracts)]
r_default = 0.045

# Compute QPE and BS prices
qpe_prices = zeros(n_test)
bs_prices = zeros(n_test)
for i in 1:n_test
    c = test_contracts[i, :]
    T = c.DTE / 365.0
    qpe_prices[i] = qpe_pricer_black76(c.Spot_Price, c.Strike, r_default, c.IV, T, 32;
                                      opt_type=c.Type[1], N_fixings=21, m_realized=0, A_realized=0.0,
                                      vol_scale=1.0, price_scale=1.0,
                                      use_seasoned=false, use_moment_match=false, mc_factor=1.0)
    if c.Type[1] == 'C'
        bs_prices[i] = black76_call(c.Spot_Price, c.Strike, T, r_default, c.IV)
    else
        bs_prices[i] = black76_put(c.Spot_Price, c.Strike, T, r_default, c.IV)
    end
end

println("\n--- Price Statistics ---")
println("QPE: mean=$(round(mean(qpe_prices), digits=4)), median=$(round(median(qpe_prices), digits=4)), std=$(round(std(qpe_prices), digits=4))")
println("Market: mean=$(round(mean(market_prices), digits=4)), median=$(round(median(market_prices), digits=4)), std=$(round(std(market_prices), digits=4))")
println("BS: mean=$(round(mean(bs_prices), digits=4)), median=$(round(median(bs_prices), digits=4)), std=$(round(std(bs_prices), digits=4))")

println("\n--- QPE vs Market ---")
ratio = qpe_prices ./ max.(market_prices, 1e-10)
println("QPE/Market ratio: mean=$(round(mean(ratio), digits=4)), median=$(round(median(ratio), digits=4)), std=$(round(std(ratio), digits=4))")

println("\n--- BS vs Market (market efficiency check) ---")
bs_ratio = market_prices ./ max.(bs_prices, 1e-10)
println("Market/BS ratio: mean=$(round(mean(bs_ratio), digits=4)), median=$(round(median(bs_ratio), digits=4)), std=$(round(std(bs_ratio), digits=4))")

println("\n--- QPE vs BS (model comparison) ---")
qpe_bs_ratio = qpe_prices ./ max.(bs_prices, 1e-10)
println("QPE/BS ratio: mean=$(round(mean(qpe_bs_ratio), digits=4)), median=$(round(median(qpe_bs_ratio), digits=4)), std=$(round(std(qpe_bs_ratio), digits=4))")

# Error analysis by moneyness and type
println("\n--- Error by Moneyness and Type ---")

call_mask = opt_types .== 'C'
put_mask = opt_types .== 'P'

bins = [0.80, 0.90, 0.95, 1.0, 1.05, 1.10, 1.20]

for type_mask = [call_mask, put_mask]
    type_name = type_mask == call_mask ? "Calls" : "Puts"
    println("\n$(type_name):")
    
    for i in 1:length(bins)-1
        bin_mask = type_mask .& (moneynesses .>= bins[i]) .& (moneynesses .< bins[i+1])
        if sum(bin_mask) >= 5
            bin_qpe = qpe_prices[bin_mask]
            bin_market = market_prices[bin_mask]
            
            # Simple scaling
            scale = median(bin_market) / median(bin_qpe)
            scaled_qpe = bin_qpe .* scale
            
            rmspe = sqrt(mean(((bin_market .- scaled_qpe) ./ max.(bin_market, 1e-8)).^2))
            
            println("  $(bins[i])-$(bins[i+1]): scale=$(round(scale, digits=2)), RMSPE=$(round(rmspe*100, digits=1))%, n=$(sum(bin_mask))")
        end
    end
end

println("\n" * "="^70)
println("CONCLUSION:")
println("The QPE model underestimates market prices by a factor of ~2.")
println("Simple scaling achieves RMSPE of ~82%, not <20%.")
println("Need a more sophisticated correction function.")
println("="^70)
