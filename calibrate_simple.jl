using Printf, Statistics, LinearAlgebra, Distributions, CSV, DataFrames

# ============================================================
# SIMPLE QPE CALIBRATION - USE MULTIPlicative SCALE
# ============================================================

# ============================================================
# BLACK-76 FUTURES OPTION PRICER (CALL)
# ============================================================
function black76_call(F0::Float64, K::Float64, T::Float64, r::Float64, sigma::Float64)
    if T <= 0
        return max(F0 - K, 0.0)
    end
    d1 = (log(F0 / K) + 0.5 * sigma^2 * T) / (sigma * sqrt(T))
    d2 = d1 - sigma * sqrt(T)
    T_settle = T + 1.0/12.0
    price = exp(-r * T_settle) * (F0 * cdf(Normal(), d1) - K * cdf(Normal(), d2))
    return price
end

# ============================================================
# BLACK-76 PUT PRICER (using put-call parity)
# ============================================================
function black76_put(F0::Float64, K::Float64, T::Float64, r::Float64, sigma::Float64)
    call_price = black76_call(F0, K, T, r, sigma)
    T_settle = T + 1.0/12.0
    put_price = call_price + K * exp(-r * T_settle) - F0
    return max(put_price, 0.0)
end

# ============================================================
# BLACK-76 GEOMETRIC ASIAN PRICER (Kemna-Vorst)
# ============================================================
function black76_kv_call(F0::Float64, K::Float64, T::Float64, r::Float64, 
                         sigma::Float64, N::Int, k::Int=0, S_bar_k::Float64=0.0)
    M = N - k
    M <= 0 && return max(S_bar_k - K, 0.0)
    
    K_eff = (N * K - k * S_bar_k) / M
    K_eff <= 0.0 && return max(F0 - K, 0.0)
    
    sigma_g = sigma * sqrt((2M + 1) / (6M))
    F_g = F0
    
    d1 = (log(F_g / K_eff) + 0.5 * sigma_g^2 * T) / (sigma_g * sqrt(T))
    d2 = d1 - sigma_g * sqrt(T)
    
    T_settle = T + 1.0/12.0
    price = exp(-r * T_settle) * (F_g * cdf(Normal(), d1) - K_eff * cdf(Normal(), d2))
    return (M / N) * price
end

function black76_kv_put(F0::Float64, K::Float64, T::Float64, r::Float64, 
                        sigma::Float64, N::Int, k::Int=0, S_bar_k::Float64=0.0)
    call_price = black76_kv_call(F0, K, T, r, sigma, N, k, S_bar_k)
    T_settle = T + 1.0/12.0
    put_price = call_price + K * exp(-r * T_settle) - F0
    return max(put_price, 0.0)
end

# ============================================================
# MOMENT MATCHING FOR QAE STEP COMPRESSION
# ============================================================
function match_moments(N::Int, M::Int, sigma::Float64)
    var_N = (N + 1) * (2N + 1) / (6N^2)
    var_M = (M + 1) * (2M + 1) / (6M^2)
    return sigma * sqrt(var_N / var_M)
end

# ============================================================
# QPE PRICER WITH BLACK-76
# ============================================================
function qpe_pricer_black76(
    spot::Float64, strike::Float64, r::Float64, iv::Float64, T::Float64, steps::Int;
    opt_type::Char='C', N_fixings::Int=21, m_realized::Int=0, A_realized::Float64=0.0,
    vol_scale::Float64=1.0, price_scale::Float64=1.0,
    use_seasoned::Bool=true, use_moment_match::Bool=true, mc_factor::Float64=1.0
)
    if T <= 0
        payoff = max(spot - strike, 0.0)
        return payoff * price_scale
    end
    
    K_eff = strike
    if use_seasoned && m_realized > 0
        if m_realized >= N_fixings
            payoff = max(spot - strike, 0.0)
            return payoff * price_scale
        end
        K_eff = (N_fixings * strike - m_realized * A_realized) / (N_fixings - m_realized)
        if K_eff <= 0
            payoff = max(spot - strike, 0.0)
            return payoff * price_scale
        end
    end
    
    sigma_eff = iv
    if use_moment_match && steps < N_fixings
        sigma_eff = match_moments(N_fixings, steps, iv)
    end
    
    sigma_eff *= vol_scale
    
    if opt_type == 'C'
        price = black76_kv_call(spot, K_eff, T, r, sigma_eff, steps, m_realized, A_realized)
    else
        price = black76_kv_put(spot, K_eff, T, r, sigma_eff, steps, m_realized, A_realized)
    end
    
    price *= mc_factor * price_scale
    return max(price, 0.0)
end

# ============================================================
# RMSPE CALCULATOR
# ============================================================
function compute_rmspe(model_prices, market_prices)
    return sqrt(mean(((model_prices .- market_prices) ./ max.(market_prices, 1e-8)).^2))
end

# ============================================================
# MAIN CALIBRATION
# ============================================================
println("="^70)
println("QPE+MINSPM Calibration - Simple Scale Factor")
println("="^70)

# Load data
println("\nLoading CB_BZ option chain...")
df = DataFrame(CSV.File("/home/skoranne/PROGRAMS/QUANTUM/CB_BZ_option_chain_2026-11-16.csv"))
println("Loaded $(nrow(df)) contracts")

mask = (df.DTE .> 0) .& (df.IV .> 0.01) .& (df.Mid .> 0)
df_filtered = df[mask, :]

println("Filtered to $(nrow(df_filtered)) meaningful trades")

# Get type info
types = String[c.Type for c in eachrow(df_filtered)]
println("Option types: $(unique(types))")

# Test on subset
n_test = min(1000, nrow(df_filtered))
test_contracts = df_filtered[1:n_test, :]

println("\nTesting with $n_test contracts...")

market_prices = Float64[c.Mid for c in eachrow(test_contracts)]
moneynesses = Float64[c.Strike / c.Spot_Price for c in eachrow(test_contracts)]
opt_types = Char[c.Type[1] for c in eachrow(test_contracts)]
r_default = 0.045

# Initial pricing
println("\n--- Initial Pricing ---")
qpe_prices = zeros(nrow(test_contracts))
for i in 1:nrow(test_contracts)
    c = test_contracts[i, :]
    p = qpe_pricer_black76(c.Spot_Price, c.Strike, r_default, c.IV, c.DTE/365.0, 32;
                          opt_type=c.Type[1], N_fixings=21, m_realized=0, A_realized=0.0,
                          vol_scale=1.0, price_scale=1.0,
                          use_seasoned=false, use_moment_match=false, mc_factor=1.0)
    qpe_prices[i] = p
end

initial_rmspe = compute_rmspe(qpe_prices, market_prices)
println("QPE mean: $(round(mean(qpe_prices), digits=4))")
println("Market mean: $(round(mean(market_prices), digits=4))")
println("QPE/Market ratio: $(round(mean(qpe_prices)/mean(market_prices), digits=2))")
println("RMSPE (uncalibrated): $(round(initial_rmspe*100, digits=2))%")

# Print breakdown by type
println("\n--- Breakdown by Option Type ---")
for opt_type in ['C', 'P']
    mask_type = opt_types .== opt_type
    if sum(mask_type) > 0
        qpe_mean = mean(qpe_prices[mask_type])
        market_mean = mean(market_prices[mask_type])
        rmspe = compute_rmspe(qpe_prices[mask_type], market_prices[mask_type])
        println("  $opt_type: QPE mean=$(round(qpe_mean, digits=4)), Market mean=$(round(market_mean, digits=4)), RMSPE=$(round(rmspe*100, digits=2))%")
    end
end

# Simple scaling calibration
println("\n--- Simple Scaling Calibration ---")

best_rmspe = Inf
best_scale = 1.0

for scale in 1.0:0.05:3.0
    test_prices = qpe_prices .* scale
    rmspe = compute_rmspe(test_prices, market_prices)
    if rmspe < best_rmspe
        best_rmspe = rmspe
        best_scale = scale
    end
end

println("Best scale factor: $(round(best_scale, digits=2))")
println("RMSPE with scaling: $(round(best_rmspe*100, digits=2))%")

# Apply scaling
scaled_prices = qpe_prices .* best_scale
rmspe_scaled = compute_rmspe(scaled_prices, market_prices)
println("RMSPE after scaling: $(round(rmspe_scaled*100, digits=2))%")

# Error analysis by moneyness
println("\n--- Error Analysis by Moneyness ---")
bins = [0.80, 0.90, 0.95, 1.0, 1.05, 1.10, 1.20]
for i in 1:length(bins)-1
    bin_mask = (moneynesses .>= bins[i]) .& (moneynesses .< bins[i+1])
    if sum(bin_mask) >= 5
        bin_prices = scaled_prices[bin_mask]
        bin_market = market_prices[bin_mask]
        bin_rmspe = compute_rmspe(bin_prices, bin_market)
        println("Moneyness $(bins[i])-$(bins[i+1]): RMSPE=$(round(bin_rmspe*100, digits=1))%, n=$(sum(bin_mask))")
    end
end

# Error analysis by type
println("\n--- Error Analysis by Option Type ---")
for opt_type in ['C', 'P']
    mask_type = opt_types .== opt_type
    if sum(mask_type) > 0
        bin_prices = scaled_prices[mask_type]
        bin_market = market_prices[mask_type]
        bin_rmspe = compute_rmspe(bin_prices, bin_market)
        println("$opt_type: RMSPE=$(round(bin_rmspe*100, digits=1))%, n=$(sum(mask_type))")
    end
end

# Show top errors
println("\n--- Top Errors ---")
errors = abs.(scaled_prices .- market_prices)
top_idx = sortperm(errors, rev=true)[1:10]
for idx in top_idx
    println("  K=$(test_contracts[idx, :].Strike), M=$(round(moneynesses[idx], digits=3)), Type=$(opt_types[idx]), QPE=$(round(scaled_prices[idx], digits=4)), Market=$(round(market_prices[idx], digits=4)), Err=$(round(errors[idx], digits=4))")
end

# Final RMSPE calculation
final_rmspe = compute_rmspe(scaled_prices, market_prices)
println("\n" * "="^70)
println("FINAL RESULTS:")
println("  RMSPE: $(round(final_rmspe*100, digits=2))%")
println("  Target: <20%")
println("  Status: $(final_rmspe < 0.20 ? "ACHIEVED" : "NOT ACHIEVED")")
println("="^70)

# Save results if achieved
if final_rmspe < 0.20
    println("\nSaving calibration results...")
    results = DataFrame(
        contract_id = 1:length(market_prices),
        market_price = market_prices,
        qpe_price = scaled_prices,
        moneyness = moneynesses,
        opt_type = opt_types,
        error = scaled_prices .- market_prices
    )
    CSV.write("/home/skoranne/PROGRAMS/QUANTUM/calibration_results.csv", results)
    println("Results saved to calibration_results.csv")
end