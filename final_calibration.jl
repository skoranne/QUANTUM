using Printf, Statistics, CSV, DataFrames

# ============================================================
# FINAL CALIBRATION - COMPREHENSIVE APPROACH
# ============================================================

function main()
    # Load and filter data
    println("Loading CB_BZ option chain data...")
    df = DataFrame(CSV.File("/home/skoranne/PROGRAMS/QUANTUM/CB_BZ_option_chain_2026-11-16.csv"))

    # Filter data: DTE > 0, IV > 0.01, Mid > 0
    mask = (df.DTE .> 0) .& (df.IV .> 0.01) .& (df.Mid .> 0)
    df_filtered = df[mask, :]

    println("Valid contracts: ", nrow(df_filtered))

    # Add derived columns
    df_filtered[!, :moneyness] = df_filtered.Strike ./ df_filtered.Spot_Price

    # ============================================================
    # ANALYZE ERROR PATTERNS
    # ============================================================
    println("\n" * "="^80)
    println("ERROR PATTERN ANALYSIS")
    println("="^80)

    # Ratio analysis
    ratios = df_filtered.Theo ./ max.(df_filtered.Mid, 1e-10)
    println("\nRatio distribution (Theo/Mid):")
    println("  Min: ", round(minimum(ratios), digits=4))
    println("  Max: ", round(maximum(ratios), digits=4))
    println("  Mean: ", round(mean(ratios), digits=4))
    println("  Median: ", round(median(ratios), digits=4))

    # DTE analysis
    println("\nRatio by DTE:")
    for (dte_max, dte_min) in [(30, 0), (90, 30), (180, 90), (Inf, 180)]
        if dte_max == Inf
            mask_d = df_filtered.DTE .> dte_min
        else
            mask_d = (df_filtered.DTE .> dte_min) .& (df_filtered.DTE .<= dte_max)
        end
        dte_subset = df_filtered[mask_d, :]
        if nrow(dte_subset) > 0
            ratios_d = dte_subset.Theo ./ max.(dte_subset.Mid, 1e-10)
            println("  ", lpad(dte_max == Inf ? ">180" : string(dte_min)*"-"*string(dte_max), 8), 
                    ": mean=", round(mean(ratios_d), digits=4), " (n=", nrow(dte_subset), ")")
        end
    end

    # IV analysis
    println("\nRatio by IV:")
    for (iv_max, iv_min) in [(0.30, 0), (0.40, 0.30), (0.60, 0.40), (Inf, 0.60)]
        if iv_max == Inf
            mask_i = df_filtered.IV .> iv_min
        else
            mask_i = (df_filtered.IV .> iv_min) .& (df_filtered.IV .<= iv_max)
        end
        iv_subset = df_filtered[mask_i, :]
        if nrow(iv_subset) > 0
            ratios_i = iv_subset.Theo ./ max.(iv_subset.Mid, 1e-10)
            println("  ", lpad(string(iv_min)*"-"*string(iv_max), 10), 
                    ": mean=", round(mean(ratios_i), digits=4), " (n=", nrow(iv_subset), ")")
        end
    end

    # ============================================================
    # CALIBRATION STRATEGY
    # ============================================================
    println("\n" * "="^80)
    println("CALIBRATION STRATEGY")
    println("="^80)

    # Based on ratio analysis:
    # - Short DTE: ratio ~1.43, need scale ~0.70
    # - Medium DTE: ratio ~1.08, need scale ~0.93
    # - Long DTE: ratio ~1.04, need scale ~0.96
    # - Very Long DTE: ratio ~1.02, need scale ~0.98
    # - Low IV: ratio ~0.68, need scale ~1.47
    # - Med IV: ratio ~0.90, need scale ~1.11
    # - High IV: ratio ~1.13, need scale ~0.88
    # - Very High IV: ratio ~1.57, need scale ~0.64

    # DTE scales
    dte_means = [1.4282, 1.0763, 1.0399, 1.020]
    dte_scales = 1.0 ./ dte_means

    # IV scales
    iv_means = [0.6788, 0.9047, 1.1344, 1.5676]
    iv_scales = 1.0 ./ iv_means

    # Apply initial calibration
    function get_dte_scale(dte)
        if dte <= 30
            return dte_scales[1]
        elseif dte <= 90
            return dte_scales[2]
        elseif dte <= 180
            return dte_scales[3]
        else
            return dte_scales[4]
        end
    end

    function get_iv_scale(iv)
        if iv <= 0.30
            return iv_scales[1]
        elseif iv <= 0.40
            return iv_scales[2]
        elseif iv <= 0.60
            return iv_scales[3]
        else
            return iv_scales[4]
        end
    end

    df_filtered[!, :calibrated_price] = [df_filtered.Theo[i] * get_dte_scale(df_filtered.DTE[i]) * 
                                          get_iv_scale(df_filtered.IV[i]) for i in 1:nrow(df_filtered)]

    # Calculate RMSPE
    mask_valid = df_filtered.Mid .> 0
    errors = (df_filtered.calibrated_price[mask_valid] .- df_filtered.Mid[mask_valid]) ./ max.(df_filtered.Mid[mask_valid], 1e-10)
    rmspe = sqrt(mean(errors.^2)) * 100

    println("\nInitial calibration RMSPE: ", round(rmspe, digits=2), "%")

    # ============================================================
    # BONUS SCALING TO MINIMIZE RMSPE
    # ============================================================
    println("\nApplying bonus scaling to minimize RMSPE...")

    # The ratio analysis shows we're consistently underestimating
    # Let's find the optimal bonus scale
    best_rmspe = rmspe
    best_bonus = 1.0

    for bonus in 1.0:0.001:1.2
        test_price = df_filtered.calibrated_price * bonus
        test_errors = (test_price[mask_valid] .- df_filtered.Mid[mask_valid]) ./ max.(df_filtered.Mid[mask_valid], 1e-10)
        test_rmspe = sqrt(mean(test_errors.^2)) * 100
        
        if test_rmspe < best_rmspe
            best_rmspe = test_rmspe
            best_bonus = bonus
        end
    end

    println("Optimal bonus scale: ", round(best_bonus, digits=3))
    println("RMSPE after bonus scaling: ", round(best_rmspe, digits=2), "%")

    # Apply final calibration
    df_filtered[!, :final_price] = df_filtered.calibrated_price * best_bonus

    # ============================================================
    # FINAL RESULTS
    # ============================================================
    println("\n" * "="^80)
    println("FINAL CALIBRATION RESULTS")
    println("="^80)

    final_errors = (df_filtered.final_price[mask_valid] .- df_filtered.Mid[mask_valid]) ./ max.(df_filtered.Mid[mask_valid], 1e-10)
    final_rmspe = sqrt(mean(final_errors.^2)) * 100

    println("Final RMSPE: ", round(final_rmspe, digits=2), "%")
    println("Target: <20%")
    println("Status: ", final_rmspe < 20 ? "✓ ACHIEVED" : "✗ NOT ACHIEVED")

    # Error breakdown by DTE
    println("\nError breakdown by DTE:")
    for (dte_max, dte_min) in [(30, 0), (90, 30), (180, 90), (Inf, 180)]
        if dte_max == Inf
            mask_d = df_filtered.DTE .> dte_min
        else
            mask_d = (df_filtered.DTE .> dte_min) .& (df_filtered.DTE .<= dte_max)
        end
        dte_subset = df_filtered[mask_d, :]
        if nrow(dte_subset) > 0
            dte_errors = (dte_subset.final_price .- dte_subset.Mid) ./ max.(dte_subset.Mid, 1e-10)
            dte_rmspe = sqrt(mean(dte_errors.^2)) * 100
            println("  ", lpad(dte_max == Inf ? ">180" : string(dte_min)*"-"*string(dte_max), 8), 
                    ": RMSPE = ", round(dte_rmspe, digits=2), "% (n=", nrow(dte_subset), ")")
        end
    end

    # Error breakdown by IV
    println("\nError breakdown by IV:")
    for (iv_max, iv_min) in [(0.30, 0), (0.40, 0.30), (0.60, 0.40), (Inf, 0.60)]
        if iv_max == Inf
            mask_i = df_filtered.IV .> iv_min
        else
            mask_i = (df_filtered.IV .> iv_min) .& (df_filtered.IV .<= iv_max)
        end
        iv_subset = df_filtered[mask_i, :]
        if nrow(iv_subset) > 0
            iv_errors = (iv_subset.final_price .- iv_subset.Mid) ./ max.(iv_subset.Mid, 1e-10)
            iv_rmspe = sqrt(mean(iv_errors.^2)) * 100
            println("  ", lpad(string(iv_min)*"-"*string(iv_max), 10), 
                    ": RMSPE = ", round(iv_rmspe, digits=2), "% (n=", nrow(iv_subset), ")")
        end
    end

    # Write calibrated output
    println("\n" * "-"^80)
    println("Writing calibrated output...")
    
    # Select columns for output
    output_df = DataFrame(
        Date = df_filtered.Date,
        Symbol = df_filtered.Symbol,
        Expiry = df_filtered.Expiry,
        DTE = df_filtered.DTE,
        Strike = df_filtered.Strike,
        Type = df_filtered.Type,
        Spot_Price = df_filtered.Spot_Price,
        Mid = df_filtered.Mid,
        Theo = df_filtered.Theo,
        IV = df_filtered.IV,
        calibrated_price = df_filtered.calibrated_price,
        final_price = df_filtered.final_price
    )
    
    CSV.write("/home/skoranne/PROGRAMS/QUANTUM/CB_BZ_calibrated.csv", output_df)

    println("\nDone!")
    return final_rmspe
end

# Run main function
rmspe = main()
