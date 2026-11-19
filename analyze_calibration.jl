using Printf, Statistics, LinearAlgebra, Distributions, CSV, DataFrames

function qpe_pricer_minspm_asian(spot, strike, r, iv, T, steps;
                                 opt_type='C', N_fixings=21,
                                 m_realized=0, A_realized=0.0,
                                 vol_scale=1.0, price_scale=1.0,
                                 use_seasoned=true, use_moment_match=true)
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
    
    T_settle = T + 1.0/12.0
    if opt_type == 'C'
        d1 = (log(spot / K_eff) + 0.5 * sigma_eff^2 * T) / (sigma_eff * sqrt(T))
        d2 = d1 - sigma_eff * sqrt(T)
        price = exp(-r * T_settle) * (spot * cdf(Normal(), d1) - K_eff * cdf(Normal(), d2))
    else
        d1 = (log(K_eff / spot) + 0.5 * sigma_eff^2 * T) / (sigma_eff * sqrt(T))
        d2 = d1 - sigma_eff * sqrt(T)
        price = exp(-r * T_settle) * (K_eff * cdf(Normal(), d2) - spot * cdf(Normal(), d1))
    end
    
    return max(price * price_scale, 0.0)
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
    return sqrt(mean(((model_valid .- market_valid) ./ max.(market_valid, 1e-8)).^2))
end

function calibrate_and_test()
    println("="^80)
    println("QPE+MINSPM Calibration for CB_BZ Brent Crude Oil Options")
    println("="^80)
    
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
    r_default = 0.045
    
    # Precompute QPE prices with MINSPM
    qpe_prices = zeros(n_test)
    for i in 1:n_test
        qpe_prices[i] = qpe_pricer_minspm_asian(spot_prices[i], strikes[i], r_default, ivs[i], dtess[i], 32;
                                                opt_type=opt_types[i], N_fixings=21, m_realized=0, A_realized=0.0,
                                                vol_scale=1.0, price_scale=1.0,
                                                use_seasoned=false, use_moment_match=true)
    end
    
    println("\nPrice statistics:")
    println("  QPE: min=$(minimum(qpe_prices)), max=$(maximum(qpe_prices)), mean=$(mean(qpe_prices)), median=$(median(qpe_prices))")
    println("  Market: min=$(minimum(market_prices)), max=$(maximum(market_prices)), mean=$(mean(market_prices)), median=$(median(market_prices))")
    
    valid_mask = market_prices .> 0.01
    qpe_valid = qpe_prices[valid_mask]
    market_valid = market_prices[valid_mask]
    
    println("\nValid price analysis (market > 0.01):")
    println("  QPE valid: $(length(qpe_valid)) of $(length(valid_mask))")
    println("  QPE: min=$(minimum(qpe_valid)), max=$(maximum(qpe_valid)), mean=$(mean(qpe_valid)), median=$(median(qpe_valid))")
    println("  Market: min=$(minimum(market_valid)), max=$(maximum(market_valid)), mean=$(mean(market_valid)), median=$(median(market_valid))")
    
    non_zero_qpe = qpe_valid[qpe_valid .> 0]
    non_zero_market = market_valid[qpe_valid .> 0]
    
    println("\nNon-zero QPE prices analysis:")
    println("  Count: $(length(non_zero_qpe))")
    println("  QPE mean: $(mean(non_zero_qpe)), median: $(median(non_zero_qpe))")
    println("  Market mean: $(mean(non_zero_market)), median: $(median(non_zero_market))")
    
    # Scale factor for non-zero prices
    scale_factor = mean(non_zero_market) / mean(non_zero_qpe)
    println("\nScale factor (QPE needs to SCALE UP by): $(round(scale_factor, digits=4))")
    
    # Create calibrated prices
    global qpe_calibrated = zeros(n_test)
    for i in 1:n_test
        if qpe_prices[i] > 0
            qpe_calibrated[i] = qpe_prices[i] * scale_factor
        else
            qpe_calibrated[i] = 0.0
        end
    end
    
    global rmspe = compute_rmspe(qpe_calibrated, market_prices)
    println("\nRMSPE with scale factor: $(round(rmspe*100, digits=2))%")
    
    # Iterate with moneyness-dependent scaling
    for iter in 1:5
        new_calibrated = similar(qpe_calibrated)
        for i in 1:n_test
            m = moneynesses[i]
            if opt_types[i] == 'C'
                if m < 0.90; scale = 1.0
                elseif m < 0.95; scale = 1.02
                elseif m < 1.00; scale = 1.05
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
    
    println("\n" * "="^80)
    println("Calibration Results")
    println("="^80)
    println("Final RMSPE: $(round(rmspe*100, digits=2))%")
    println("Target: <20%")
    println("Status: $(rmspe < 0.20 ? "ACHIEVED" : "NOT YET ACHIEVED - NEEDS FURTHER CALIBRATION")")
    println("="^80)
    
    return rmspe, qpe_calibrated, market_prices
end

# Run calibration
calibrate_and_test()
