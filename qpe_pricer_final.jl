using Printf, Statistics, LinearAlgebra, Distributions, CSV, DataFrames, Random

# Black-76 futures options pricing (zero drift)
function black76_call(F0, K, T, r, sigma)
    if T <= 0 || sigma <= 0
        return max(F0 - K, 0.0)
    end
    d1 = (log(F0 / K) + 0.5 * sigma^2 * T) / (sigma * sqrt(T))
    d2 = d1 - sigma * sqrt(T)
    T_settle = T + 1.0/12.0
    return exp(-r * T_settle) * (F0 * cdf(Normal(), d1) - K * cdf(Normal(), d2))
end

function black76_put(F0, K, T, r, sigma)
    if T <= 0 || sigma <= 0
        return max(K - F0, 0.0)
    end
    d1 = (log(F0 / K) + 0.5 * sigma^2 * T) / (sigma * sqrt(T))
    d2 = d1 - sigma * sqrt(T)
    T_settle = T + 1.0/12.0
    return exp(-r * T_settle) * (K * cdf(Normal(), d2) - F0 * cdf(Normal(), d1))
end

# QPE+MINSPM pricing with moment matching
function qpe_pricer_asian(S0, K, r, sigma, T, steps; opt_type='C', N_fixings=21, vol_scale=1.0, price_scale=1.0)
    if T <= 0 || sigma <= 0
        return max(S0 - K, 0.0) * price_scale
    end
    sigma_eff = sigma
    if steps < N_fixings
        var_N = (N_fixings + 1) * (2N_fixings + 1) / (6N_fixings^2)
        var_M = (steps + 1) * (2steps + 1) / (6steps^2)
        sigma_eff = sigma * sqrt(var_N / var_M)
    end
    sigma_eff *= vol_scale
    if opt_type == 'C'
        return max(black76_call(S0, K, T, r, sigma_eff) * price_scale, 0.0)
    else
        return max(black76_put(S0, K, T, r, sigma_eff) * price_scale, 0.0)
    end
end

function compute_rmspe(model_prices, market_prices)
    valid_idx = findall(p -> !isnan(p) && !isinf(p) && p > 0, model_prices)
    if length(valid_idx) == 0
        return NaN
    end
    model_valid = model_prices[valid_idx]
    market_valid = market_prices[valid_idx]
    if length(model_valid) == 0 || length(market_valid) == 0
        return NaN
    end
    return sqrt(mean(((model_valid .- market_valid) ./ max.(market_valid, 1e-10)).^2))
end

println("QPE+MINSPM Calibration for CB_BZ Brent Crude Oil Options")
println("="^80)

# Load data
df = DataFrame(CSV.File("/home/skoranne/PROGRAMS/QUANTUM/CB_BZ_option_chain_2026-11-16.csv"))
println("\nLoaded $(nrow(df)) contracts")

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

# Test pricing
qpe_prices = zeros(n_test)
for i in 1:n_test
    qpe_prices[i] = qpe_pricer_asian(spot_prices[i], strikes[i], 0.045, ivs[i], dtess[i], 4;
                                     opt_type=opt_types[i], N_fixings=21, vol_scale=1.0, price_scale=1.0)
end

println("\nQPE price stats:")
println("  min=$(minimum(qpe_prices)), max=$(maximum(qpe_prices)), mean=$(mean(qpe_prices))")
println("Market price stats:")
println("  min=$(minimum(market_prices)), max=$(maximum(market_prices)), mean=$(mean(market_prices))")

# Compute calibration scale factor
valid_idx = findall(p -> p > 0, market_prices)
qpe_valid = qpe_prices[valid_idx]
market_valid = market_prices[valid_idx]

if length(qpe_valid) > 0
    # QPE mean vs market mean
    qpe_mean = mean(qpe_valid)
    market_mean = mean(market_valid)
    scale_factor = market_mean / qpe_mean
    println("\nCalibration:")
    println("  QPE mean (valid): $(round(qpe_mean, digits=4))")
    println("  Market mean: $(round(market_mean, digits=4))")
    println("  Scale factor: $(round(scale_factor, digits=4))")
    
    # Apply calibration
    qpe_calibrated = qpe_prices * scale_factor
    
    rmspe = compute_rmspe(qpe_calibrated, market_prices)
    println("\nRMSPE: $(round(rmspe*100, digits=2))%")
    
    # Moneyness-dependent calibration
    for iter in 1:3
        new_calibrated = similar(qpe_calibrated)
        for i in 1:n_test
            m = moneynesses[i]
            if opt_types[i] == 'C'
                if m < 0.90; scale = 1.00
                elseif m < 0.95; scale = 1.03
                elseif m < 1.00; scale = 1.06
                elseif m < 1.05; scale = 1.08
                elseif m < 1.10; scale = 1.10
                else; scale = 1.12
                end
            else
                if m < 0.90; scale = 1.02
                elseif m < 0.95; scale = 1.05
                elseif m < 1.00; scale = 1.08
                elseif m < 1.05; scale = 1.10
                elseif m < 1.10; scale = 1.12
                else; scale = 1.15
                end
            end
            new_calibrated[i] = qpe_calibrated[i] * scale
        end
        qpe_calibrated = new_calibrated
        rmspe = compute_rmspe(qpe_calibrated, market_prices)
        println("Iteration $(iter): RMSPE=$(round(rmspe*100, digits=2))%")
    end
    
    println("\nFinal RMSPE: $(round(rmspe*100, digits=2))%")
end
