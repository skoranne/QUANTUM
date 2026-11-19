using Printf, Statistics, CSV, DataFrames

# ============================================================
# IMPROVED CALIBRATION FOR CB_BZ OPTIONS
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
    # OPTIMIZED CALIBRATION PARAMETERS
    # Based on ratio analysis:
    # DTE groups: short=1.428, med=1.076, long=1.040, very_long=1.020
    # Target ratio = 1.0, so scales = 1/ratio
    # ============================================================

    # DTE-based scaling factors (inverse of mean ratios)
    dte_scales = [1.0 / 1.4282, 1.0 / 1.0763, 1.0 / 1.0399, 1.0 / 1.020]
    
    # IV-based scaling factors (inverse of mean ratios)
    iv_scales = [1.0 / 0.6788, 1.0 / 0.9047, 1.0 / 1.1344, 1.0 / 1.5676]

    println("\nOptimal scales:")
    println("DTE scales: ", round.(dte_scales, digits=3))
    println("IV scales: ", round.(iv_scales, digits=3))

    # Calibration function
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

    # Apply calibration
    df_filtered[!, :calibrated_price] = [df_filtered.Theo[i] * get_dte_scale(df_filtered.DTE[i]) * 
                                          get_iv_scale(df_filtered.IV[i]) for i in 1:nrow(df_filtered)]

    # Calculate RMSPE
    mask_valid = df_filtered.Mid .> 0
    errors = (df_filtered.calibrated_price[mask_valid] .- df_filtered.Mid[mask_valid]) ./ max.(df_filtered.Mid[mask_valid], 1e-10)
    rmspe = sqrt(mean(errors.^2)) * 100

    println("\n" * "="^80)
    println("CALIBRATION RESULTS")
    println("="^80)
    println("RMSPE: ", round(rmspe, digits=2), "%")
    println("Target: <20%")
    println("Status: ", rmspe < 20 ? "✓ ACHIEVED" : "✗ NOT ACHIEVED")

    # Error analysis by DTE
    println("\nError breakdown by DTE:")
    dte_cutoffs = [30, 90, 180, Inf]
    dte_labels = ["short (<30d)", "med (30-90d)", "long (90-180d)", "very_long (>180d)"]

    for (i, dte_label) in enumerate(dte_labels)
        if i == 1
            mask_d = df_filtered.DTE .<= dte_cutoffs[i]
        else
            mask_d = df_filtered.DTE .> dte_cutoffs[i-1]
        end
        dte_subset = df_filtered[mask_d, :]
        if nrow(dte_subset) > 0
            dte_errors = (dte_subset.calibrated_price .- dte_subset.Mid) ./ max.(dte_subset.Mid, 1e-10)
            dte_rmspe = sqrt(mean(dte_errors.^2)) * 100
            println("  ", lpad(dte_label, 15), ": RMSPE = ", round(dte_rmspe, digits=2), "% (n=", nrow(dte_subset), ")")
        end
    end

    # Error analysis by IV
    println("\nError breakdown by IV:")
    iv_cutoffs = [0.30, 0.40, 0.60, Inf]
    iv_labels = ["low (<30%)", "med (30-40%)", "high (40-60%)", "very_high (>60%)"]

    for (i, iv_label) in enumerate(iv_labels)
        if i == 1
            mask_i = df_filtered.IV .<= iv_cutoffs[i]
        else
            mask_i = df_filtered.IV .> iv_cutoffs[i-1]
        end
        iv_subset = df_filtered[mask_i, :]
        if nrow(iv_subset) > 0
            iv_errors = (iv_subset.calibrated_price .- iv_subset.Mid) ./ max.(iv_subset.Mid, 1e-10)
            iv_rmspe = sqrt(mean(iv_errors.^2)) * 100
            println("  ", lpad(iv_label, 15), ": RMSPE = ", round(iv_rmspe, digits=2), "% (n=", nrow(iv_subset), ")")
        end
    end

    # Error analysis by moneyness
    println("\nError breakdown by moneyness:")
    moneyness_cutoffs = [0.95, 1.05, Inf]
    moneyness_labels = ["otm (<95%)", "atm (95-105%)", "itm (>105%)"]

    for (i, m_label) in enumerate(moneyness_labels)
        if i == 1
            mask_m = df_filtered.moneyness .<= moneyness_cutoffs[i]
        else
            mask_m = df_filtered.moneyness .> moneyness_cutoffs[i-1]
        end
        m_subset = df_filtered[mask_m, :]
        if nrow(m_subset) > 0
            m_errors = (m_subset.calibrated_price .- m_subset.Mid) ./ max.(m_subset.Mid, 1e-10)
            m_rmspe = sqrt(mean(m_errors.^2)) * 100
            println("  ", lpad(m_label, 15), ": RMSPE = ", round(m_rmspe, digits=2), "% (n=", nrow(m_subset), ")")
        end
    end

    # Write calibrated output
    println("\n" * "-"^80)
    println("Writing calibrated output...")
    CSV.write("/home/skoranne/PROGRAMS/QUANTUM/CB_BZ_calibrated.csv", df_filtered)

    println("\nDone!")
    return rmspe
end

# Run main function
rmspe = main()
