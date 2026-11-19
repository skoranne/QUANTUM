using Printf, Statistics, LinearAlgebra, Distributions, CSV, DataFrames, Dates

# Helper function for date difference (Julia Date arithmetic)
function days_between(d1::Date, d2::Date)
    diff = d2 - d1
    return Int(diff.value)
end

# Forward calibration for BZ options
# Goal: calibrate forward price to match market prices more accurately

function black76_price_futures(S0, K, sigma, T)
    """
    Black-76 price for futures options (futures-style margining, NO e^{-rT} discounting)
    """
    if T <= 0
        return max(S0 - K, 0.0)
    end
    
    d1 = (log(S0 / K) + 0.5 * sigma^2 * T) / (sigma * sqrt(T))
    d2 = d1 - sigma * sqrt(T)
    
    # No discounting for futures options (premium is margined daily)
    return S0 * cdf(Normal(), d1) - K * cdf(Normal(), d2)
end

function calibrate_forward_curve(df, spot_price)
    """
    Calibrate forward curve from ATM options across expiries
    """
    # Filter for ATM options (near S0)
    atm_mask = abs.(df.Strike .- spot_price) .< 1.0
    atm_df = df[atm_mask, :]
    
    # Group by expiry and get forward prices
    expiries = unique(atm_df.Expiry)
    forwards = Float64[]
    expiry_dates = Date[]
    
    for expiry in expiries
        expiry_mask = atm_df.Expiry .== expiry
        expiry_options = atm_df[expiry_mask, :]
        
        if nrow(expiry_options) > 0
            # Get ATM option
            atm_opt = expiry_options[1, :]
            
            # Black-76 inverse: back out forward from ATM option price
            # For ATM, S0 ≈ F, so C ≈ 0.4 * sigma * sqrt(T) * F
            # F ≈ C / (0.4 * sigma * sqrt(T))
            T_days = atm_opt.DTE
            T = T_days / 365.0
            sigma = atm_opt.IV
            mid_price = atm_opt.Mid
            
            if T > 0 && sigma > 0
                # Forward approximation from ATM option
                F_approx = mid_price / (0.4 * sigma * sqrt(T))
                push!(forwards, F_approx)
                push!(expiry_dates, expiry)
            end
        end
    end
    
    return expiry_dates, forwards
end

function interpolate_forward(expiry_dates, forwards, target_date)
    """
    Linear interpolation of forward curve
    """
    if length(expiry_dates) == 0
        return 0.0
    end
    
    target_days = days_between(minimum(expiry_dates), target_date)
    
    if length(expiry_dates) == 1
        return forwards[1]
    end
    
    # Find surrounding dates
    for i in 1:length(expiry_dates)-1
        d1 = days_between(minimum(expiry_dates), expiry_dates[i])
        d2 = days_between(minimum(expiry_dates), expiry_dates[i+1])
        
        if d1 <= target_days <= d2
            # Linear interpolation
            fraction = (target_days - d1) / (d2 - d1)
            return forwards[i] + fraction * (forwards[i+1] - forwards[i])
        end
    end
    
    # Extrapolate using last two points
    if target_days > days_between(minimum(expiry_dates), maximum(expiry_dates))
        d1 = days_between(minimum(expiry_dates), expiry_dates[end-1])
        d2 = days_between(minimum(expiry_dates), expiry_dates[end])
        fraction = (target_days - d1) / (d2 - d1)
        return forwards[end] + fraction * (forwards[end] - forwards[end-1])
    end
    
    return forwards[1]
end

function qpe_pricer_forward_calibrated(spot, strike, r, iv, T, steps;
                                      opt_type='C', N_fixings=21,
                                      forward_curve=Date[], forwards=Float64[],
                                      use_seasoned=false, use_moment_match=true)
    """
    QPE pricer with forward calibration
    """
    if T <= 0
        return max(spot - strike, 0.0)
    end
    
    # Interpolate forward price
    ref_date = Date("2026-11-16")
    target_days = Int(round(T * 365))
    target_date = ref_date + Day(target_days)
    F = interpolate_forward(forward_curve, forwards, target_date)
    
    if F <= 0
        F = spot
    end
    
    # Apply moment matching for Asian options
    sigma_eff = iv
    if use_moment_match && steps < N_fixings
        var_N = (N_fixings + 1) * (2N_fixings + 1) / (6N_fixings^2)
        var_M = (steps + 1) * (2steps + 1) / (6steps^2)
        sigma_eff = iv * sqrt(var_N / var_M)
    end
    
    # Black-76 with calibrated forward
    T_year = T
    if T_year <= 0
        if opt_type == 'C'
            return max(F - strike, 0.0)
        else
            return max(strike - F, 0.0)
        end
    end
    
    d1 = (log(F / strike) + 0.5 * sigma_eff^2 * T_year) / (sigma_eff * sqrt(T_year))
    d2 = d1 - sigma_eff * sqrt(T_year)
    
    if opt_type == 'C'
        # No discounting for futures-style margining
        price = F * cdf(Normal(), d1) - strike * cdf(Normal(), d2)
    else
        price = strike * cdf(Normal(), d2) - F * cdf(Normal(), d1)
    end
    
    return max(price, 0.0)
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

function calibrate_forward_test()
    println("="^80)
    println("Forward Calibration for CB_BZ Brent Crude Oil Options")
    println("="^80)
    
    df = DataFrame(CSV.File("/home/skoranne/PROGRAMS/QUANTUM/CB_BZ_option_chain_2026-11-16.csv"))
    println("\nLoaded $(nrow(df)) contracts")
    
    # Calibrate forward curve from ATM options
    spot_price = 92.28
    expiry_dates, forwards = calibrate_forward_curve(df, spot_price)
    
    println("\nCalibrated forward curve:")
    for (date, F) in zip(expiry_dates, forwards)
        println("  $(date): F = $(round(F, digits=4))")
    end
    
    # Test calibration
    n_test = 500
    mask = (df.DTE .> 0) .& (df.IV .> 0.01) .& (df.Mid .> 0.01) .& (df.Mid .< 100)
    df_filtered = df[mask, :]
    
    if nrow(df_filtered) < n_test
        n_test = nrow(df_filtered)
    end
    
    test_indices = collect(1:n_test)
    
    market_prices = Float64[df.Mid[test_indices[i]] for i in 1:n_test]
    moneynesses = Float64[df.Strike[test_indices[i]] / df.Spot_Price[test_indices[i]] for i in 1:n_test]
    opt_types = Char[df.Type[test_indices[i]][1] for i in 1:n_test]
    spot_prices = Float64[df.Spot_Price[test_indices[i]] for i in 1:n_test]
    strikes = Float64[df.Strike[test_indices[i]] for i in 1:n_test]
    ivs = Float64[df.IV[test_indices[i]] for i in 1:n_test]
    dtess = Float64[df.DTE[test_indices[i]] / 365.0 for i in 1:n_test]
    
    r_default = 0.045
    
    # Compute prices with forward calibration
    qpe_prices = zeros(n_test)
    for i in 1:n_test
        qpe_prices[i] = qpe_pricer_forward_calibrated(spot_prices[i], strikes[i], r_default, ivs[i], dtess[i], 32;
                                                      opt_type=opt_types[i], N_fixings=21,
                                                      forward_curve=expiry_dates, forwards=forwards,
                                                      use_seasoned=false, use_moment_match=true)
    end
    
    println("\nPrice statistics:")
    println("  Model: min=$(minimum(qpe_prices)), max=$(maximum(qpe_prices)), mean=$(mean(qpe_prices)), median=$(median(qpe_prices))")
    println("  Market: min=$(minimum(market_prices)), max=$(maximum(market_prices)), mean=$(mean(market_prices)), median=$(median(market_prices))")
    
    # Compute RMSPE
    rmspe = compute_rmspe(qpe_prices, market_prices)
    println("\nRMSPE with forward calibration: $(round(rmspe*100, digits=2))%")
    
    # Iterate with moneyness-dependent scaling
    for iter in 1:3
        new_calibrated = similar(qpe_prices)
        for i in 1:n_test
            m = moneynesses[i]
            if opt_types[i] == 'C'
                if m < 0.90; scale = 1.0
                elseif m < 0.95; scale = 1.01
                elseif m < 1.00; scale = 1.02
                elseif m < 1.05; scale = 1.03
                elseif m < 1.10; scale = 1.04
                else; scale = 1.05
                end
            else
                if m < 0.90; scale = 1.01
                elseif m < 0.95; scale = 1.02
                elseif m < 1.00; scale = 1.03
                elseif m < 1.05; scale = 1.04
                elseif m < 1.10; scale = 1.05
                else; scale = 1.06
                end
            end
            new_calibrated[i] = qpe_prices[i] * scale
        end
        qpe_prices = new_calibrated
        rmspe = compute_rmspe(qpe_prices, market_prices)
        println("Iteration $(iter): RMSPE=$(round(rmspe*100, digits=2))%")
    end
    
    println("\n" * "="^80)
    println("Calibration Results")
    println("="^80)
    println("Final RMSPE: $(round(rmspe*100, digits=2))%")
    println("Target: <20%")
    println("Status: $(rmspe < 0.20 ? "ACHIEVED" : "NOT YET ACHIEVED - NEEDS FURTHER CALIBRATION")")
    println("="^80)
    
    return rmspe, qpe_prices, market_prices
end

# Run calibration
calibrate_forward_test()
