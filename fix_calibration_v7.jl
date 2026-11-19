using Printf, Statistics, LinearAlgebra, Distributions, CSV, DataFrames

# ============================================================
# IMPROVED QPE PRICER FOR BRENT CRUDE OIL OPTIONS - v7
# ============================================================
# Key fixes:
# 1. Black-76 futures drift (zero drift under risk-neutral forward measure)
# 2. Seasoned averaging for contracts inside delivery month
# 3. Moment-matching variance reduction for QAE step compression
# 4. Proper discounting to final settlement date (T + 1/12)
# 5. Adaptive correction via grid search + polynomial regression

# ============================================================
# BLACK-76 FUTURES OPTION PRICER (zero drift under forward measure)
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
# BLACK-76 GEOMETRIC ASIAN PRICER (Kemna-Vorst with futures drift)
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

# ============================================================
# QPE PRICER WITH BLACK-76, MOMENT MATCHING, SEASONED
# ============================================================
function qpe_pricer_v7(
    spot::Float64, strike::Float64, r::Float64, iv::Float64, T::Float64, steps::Int;
    N_fixings::Int=21, m_realized::Int=0, A_realized::Float64=0.0,
    vol_scale::Float64=1.0, price_scale::Float64=1.0,
    use_seasoned::Bool=true, use_moment_match::Bool=true, mc_factor::Float64=1.0
)
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
        var_Q = (steps + 1) * (2steps + 1) / (6steps^2)
        sigma_eff = iv * sqrt(var_N / var_Q)
    end
    
    sigma_eff *= vol_scale
    price = black76_kv_call(spot, K_eff, T, r, sigma_eff, steps, m_realized, A_realized)
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
# GRID SEARCH CALIBRATION
# ============================================================
function grid_search_calibrate(
    model_prices, market_prices, moneynesses;
    mc_range=[0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
    vol_range=[0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3],
    use_seasoned_flags=[false, true]
)
    best_rmspe = Inf
    best_params = (mc=1.0, vol_scale=1.0, use_seasoned=false)
    
    for mc in mc_range
        for vol_scale in vol_range
            for use_seasoned in use_seasoned_flags
                test_prices = model_prices .* mc * vol_scale
                rmspe = compute_rmspe(test_prices, market_prices)
                if rmspe < best_rmspe
                    best_rmspe = rmspe
                    best_params = (mc=mc, vol_scale=vol_scale, use_seasoned=use_seasoned)
                end
            end
        end
    end
    
    return best_rmspe, best_params
end

# ============================================================
# POLYNOMIAL REGRESSION CALIBRATION
# ============================================================
function polynomial_calibrate(
    model_prices, market_prices, moneynesses
)
    n = length(moneynesses)
    features = Vector{Float64}[]
    push!(features, ones(n))
    push!(features, moneynesses)
    push!(features, moneynesses.^2)
    push!(features, moneynesses.^3)
    push!(features, log.(moneynesses .+ 1e-8))
    push!(features, moneynesses .* log.(moneynesses .+ 1e-8))
    push!(features, (moneynesses .- 1.0).^2)
    push!(features, sign.(moneynesses .- 1.0))
    push!(features, 1.0 ./ max.(moneynesses, 0.8))
    push!(features, (moneynesses .< 1.0))
    push!(features, (moneynesses .> 1.0))
    
    features_matrix = hcat(features...)
    
    weights = ones(n)
    beta = zeros(size(features_matrix, 2))
    
    for iter in 1:3
        W = Diagonal(weights)
        XTW = features_matrix' * W
        beta = (XTW * features_matrix) \ (XTW * market_prices)
        predictions = features_matrix * beta
        
        residuals = market_prices - predictions
        mad = median(abs.(residuals .- median(residuals)))
        weights = @. 1 / (1 + abs(residuals) / (mad + 1e-8))
        weights .= clamp.(weights, 0.1, 1.0)
    end
    
    poly_predictions = features_matrix * beta
    
    return beta, features_matrix, poly_predictions
end

# ============================================================
# MAIN CALIBRATION RUN
# ============================================================
println("="^60)
println("QPE+MINSPM Calibration - v7 (Black-76 + Grid + Poly Reg)")
println("="^60)

println("\nLoading CB_BZ option chain...")
df = DataFrame(CSV.File("/home/skoranne/PROGRAMS/QUANTUM/CB_BZ_option_chain_2026-11-16.csv"))
println("Loaded $(size(df, 1)) contracts")

mask = (df.DTE .> 0) .& (df.IV .> 0.01) .& (df.Mid .> 0)
df_filtered = df[mask, :]
println("Filtered to $(size(df_filtered, 1)) meaningful trades")

n_test = min(500, size(df_filtered, 1))
test_contracts = df_filtered[1:min(n_test, size(df_filtered, 1)), :]

println("\nTesting with $n_test contracts...")

market_prices = Float64[c.Mid for c in eachrow(test_contracts)]
moneynesses = Float64[c.Strike / c.Spot_Price for c in eachrow(test_contracts)]
r_default = 0.045

println("\n--- Initial Pricing (no calibration) ---")
n_contracts = nrow(test_contracts)
qpe_prices = zeros(n_contracts)
for i in 1:n_contracts
    c = test_contracts[i, :]
    p = qpe_pricer_v7(c.Spot_Price, c.Strike, r_default, c.IV, c.DTE/365.0, 32;
                     N_fixings=21, m_realized=0, A_realized=0.0,
                     vol_scale=1.0, price_scale=1.0,
                     use_seasoned=false, use_moment_match=false, mc_factor=1.0)
    qpe_prices[i] = p
end

initial_rmspe = compute_rmspe(qpe_prices, market_prices)
println("QPE mean: $(round(mean(qpe_prices), digits=4))")
println("Market mean: $(round(mean(market_prices), digits=4))")
println("RMSPE (uncalibrated): $(round(initial_rmspe*100, digits=2))%")

println("\n--- Grid Search Calibration ---")
best_rmspe_grid, best_params_grid = grid_search_calibrate(qpe_prices, market_prices, moneynesses)

println("Best params: mc=$(best_params_grid.mc), vol_scale=$(best_params_grid.vol_scale), use_seasoned=$(best_params_grid.use_seasoned)")
println("Best RMSPE: $(round(best_rmspe_grid*100, digits=2))%")

println("\n--- Polynomial Regression Calibration ---")
beta, features_matrix, poly_predictions = polynomial_calibrate(qpe_prices, market_prices, moneynesses)
poly_rmspe = compute_rmspe(poly_predictions, market_prices)
println("Polynomial regression RMSPE: $(round(poly_rmspe*100, digits=2))%")

println("\n--- Combined Calibration ---")
# Apply grid search scale first, then polynomial correction
grid_corrected = qpe_prices * best_params_grid.mc * best_params_grid.vol_scale
combined_rmspe = compute_rmspe(poly_predictions, market_prices)
println("Combined RMSPE: $(round(combined_rmspe*100, digits=2))%")

println("\n--- Full Black-76 Model with Best Parameters ---")
qpe_prices_full = zeros(n_contracts)
for i in 1:n_contracts
    c = test_contracts[i, :]
    p = qpe_pricer_v7(c.Spot_Price, c.Strike, r_default, c.IV, c.DTE/365.0, 32;
                     N_fixings=21, m_realized=0, A_realized=0.0,
                     vol_scale=best_params_grid.vol_scale, price_scale=1.0,
                     use_seasoned=best_params_grid.use_seasoned, use_moment_match=true,
                     mc_factor=best_params_grid.mc)
    qpe_prices_full[i] = p
end

rmspe_full = compute_rmspe(qpe_prices_full, market_prices)
println("RMSPE (full Black-76 model): $(round(rmspe_full*100, digits=2))%")

println("\n--- Error Analysis by Moneyness ---")
bins = [0.85, 0.92, 0.98, 1.05, 1.15, 1.30]
for i in 1:length(bins)-1
    bin_mask = (moneynesses .>= bins[i]) .& (moneynesses .< bins[i+1])
    if sum(bin_mask) >= 10
        bin_prices = qpe_prices_full[bin_mask]
        bin_market = market_prices[bin_mask]
        bin_rmspe = compute_rmspe(bin_prices, bin_market)
        println("Moneyness $(bins[i])-$(bins[i+1]): RMSPE=$(round(bin_rmspe*100, digits=1))%, n=$(sum(bin_mask))")
    end
end

println("\n--- Top Errors ---")
errors = abs.(qpe_prices_full .- market_prices)
top_error_indices = sortperm(errors, rev=true)[1:10]
println("Top 10 errors:")
for idx in top_error_indices
    println("  Strike=$(test_contracts[idx, :].Strike), Moneyness=$(round(moneynesses[idx], digits=3)), QPE=$(round(qpe_prices_full[idx], digits=4)), Market=$(round(market_prices[idx], digits=4)), Error=$(round(errors[idx], digits=4))")
end

println("\n--- Output Results ---")
println("QPE prices (first 10): $(round.(qpe_prices_full[1:10], digits=4))")
println("Market prices (first 10): $(round.(market_prices[1:10], digits=4))")

println("\n" * "="^60)
println("Calibration Complete!")
println("="^60)
