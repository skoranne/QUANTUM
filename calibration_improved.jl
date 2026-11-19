using Printf, Statistics, LinearAlgebra, Distributions, CSV, DataFrames

# ============================================================
# IMPROVED CALIBRATION WITH DTE AND IV SCALING
# ============================================================

function main()
    # Load data
    println("Loading CB_BZ option chain data...")
    df = DataFrame(CSV.File("/home/skoranne/PROGRAMS/QUANTUM/CB_BZ_option_chain_2026-11-16.csv"))

    # Filter data: DTE > 0, IV > 0.01, Mid > 0
    mask = (df.DTE .> 0) .& (df.IV .> 0.01) .& (df.Mid .> 0)
    df_filtered = df[mask, :]

    println("Valid contracts: ", nrow(df_filtered))

    # Add derived columns
    df_filtered[!, :moneyness] = df_filtered.Strike ./ df_filtered.Spot_Price

    # ============================================================
    # GRID SEARCH CALIBRATION
    # ============================================================

    function compute_rmspe(dte_scales, iv_scales, df)
        total_sq_error = 0.0
        count = 0
        
        for i in 1:nrow(df)
            if df.Mid[i] > 0
                dte = df.DTE[i]
                iv = df.IV[i]
                
                # Find DTE scale
                if dte <= 30
                    d_scale = dte_scales[1]
                elseif dte <= 90
                    d_scale = dte_scales[2]
                elseif dte <= 180
                    d_scale = dte_scales[3]
                else
                    d_scale = dte_scales[4]
                end
                
                # Find IV scale
                if iv <= 0.30
                    i_scale = iv_scales[1]
                elseif iv <= 0.40
                    i_scale = iv_scales[2]
                elseif iv <= 0.60
                    i_scale = iv_scales[3]
                else
                    i_scale = iv_scales[4]
                end
                
                calibrated = df.Theo[i] * d_scale * i_scale
                rel_error = (calibrated - df.Mid[i]) / df.Mid[i]
                total_sq_error += rel_error^2
                count += 1
            end
        end
        
        return sqrt(total_sq_error / count) * 100
    end

    # Grid search with fine resolution
    println("\nPerforming grid search for optimal scales...")

    best_rmspe = Inf
    best_dte_scales = [1.0, 1.0, 1.0, 1.0]
    best_iv_scales = [1.0, 1.0, 1.0, 1.0]

    # Coarse search
    println("Coarse search...")
    for d1 in 1.0:0.1:2.0
        for d2 in 0.9:0.1:1.5
            for d3 in 0.9:0.05:1.2
                for d4 in 0.9:0.05:1.1
                    for i1 in 0.8:0.1:1.4
                        for i2 in 0.7:0.1:1.3
                            for i3 in 0.8:0.1:1.4
                                for i4 in 0.4:0.1:0.9
                                    dte_scales = [d1, d2, d3, d4]
                                    iv_scales = [i1, i2, i3, i4]
                                    rmspe = compute_rmspe(dte_scales, iv_scales, df_filtered)
                                    
                                    if rmspe < best_rmspe
                                        best_rmspe = rmspe
                                        best_dte_scales = copy(dte_scales)
                                        best_iv_scales = copy(iv_scales)
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    println("Coarse search best RMSPE: ", round(best_rmspe, digits=2), "%")
    println("DTE scales: ", round.(best_dte_scales, digits=3))
    println("IV scales: ", round.(best_iv_scales, digits=3))

    # Fine search around best values
    println("\nFine search around best values...")
    for d1 in best_dte_scales[1]-0.15:0.02:best_dte_scales[1]+0.15
        for d2 in best_dte_scales[2]-0.1:0.02:best_dte_scales[2]+0.1
            for d3 in best_dte_scales[3]-0.05:0.01:best_dte_scales[3]+0.05
                for d4 in best_dte_scales[4]-0.05:0.01:best_dte_scales[4]+0.05
                    for i1 in best_iv_scales[1]-0.15:0.02:best_iv_scales[1]+0.15
                        for i2 in best_iv_scales[2]-0.15:0.02:best_iv_scales[2]+0.15
                            for i3 in best_iv_scales[3]-0.15:0.02:best_iv_scales[3]+0.15
                                for i4 in best_iv_scales[4]-0.15:0.02:best_iv_scales[4]+0.15
                                    dte_scales = [d1, d2, d3, d4]
                                    iv_scales = [i1, i2, i3, i4]
                                    rmspe = compute_rmspe(dte_scales, iv_scales, df_filtered)
                                    
                                    if rmspe < best_rmspe
                                        best_rmspe = rmspe
                                        best_dte_scales = copy(dte_scales)
                                        best_iv_scales = copy(iv_scales)
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    println("Fine search best RMSPE: ", round(best_rmspe, digits=2), "%")
    println("DTE scales: ", round.(best_dte_scales, digits=3))
    println("IV scales: ", round.(best_iv_scales, digits=3))

    # ============================================================
    # FINAL CALIBRATION
    # ============================================================
    # Apply final calibrated prices
    function get_dte_scale(dte, scales)
        if dte <= 30
            return scales[1]
        elseif dte <= 90
            return scales[2]
        elseif dte <= 180
            return scales[3]
        else
            return scales[4]
        end
    end

    function get_iv_scale(iv, scales)
        if iv <= 0.30
            return scales[1]
        elseif iv <= 0.40
            return scales[2]
        elseif iv <= 0.60
            return scales[3]
        else
            return scales[4]
        end
    end

    df_filtered[!, :calibrated_price] = [get_dte_scale(df_filtered.DTE[i], best_dte_scales) * 
                                          get_iv_scale(df_filtered.IV[i], best_iv_scales) * 
                                          df_filtered.Theo[i] for i in 1:nrow(df_filtered)]

    # Calculate RMSPE
    mask_valid = df_filtered.Mid .> 0
    errors = (df_filtered.calibrated_price[mask_valid] .- df_filtered.Mid[mask_valid]) ./ max.(df_filtered.Mid[mask_valid], 1e-10)
    rmspe = sqrt(mean(errors.^2)) * 100

    println("\n" * "="^80)
    println("FINAL CALIBRATION RESULTS")
    println("="^80)
    println("RMSPE: ", round(rmspe, digits=2), "%")
    println("Target: <20%")
    println("Status: ", rmspe < 20 ? "✓ ACHIEVED" : "✗ NOT ACHIEVED")

    # Error analysis by category
    println("\nError breakdown by DTE:")
    dte_bins = [0, 30, 90, 180, Inf]
    dte_labels = ["short (<30d)", "med (30-90d)", "long (90-180d)", "very_long (>180d)"]

    for (i, dte_label) in enumerate(dte_labels)
        if i == 1
            mask_d = df_filtered.DTE .<= dte_bins[i+1]
        elseif i == length(dte_bins) - 1
            mask_d = df_filtered.DTE .> dte_bins[i]
        else
            mask_d = (df_filtered.DTE .> dte_bins[i]) .& (df_filtered.DTE .<= dte_bins[i+1])
        end
        dte_subset = df_filtered[mask_d, :]
        if nrow(dte_subset) > 0
            dte_errors = (dte_subset.calibrated_price .- dte_subset.Mid) ./ max.(dte_subset.Mid, 1e-10)
            dte_rmspe = sqrt(mean(dte_errors.^2)) * 100
            println("  ", lpad(dte_label, 15), ": RMSPE = ", round(dte_rmspe, digits=2), "% (n=", nrow(dte_subset), ")")
        end
    end

    println("\nError breakdown by IV:")
    iv_bins = [0.0, 0.3, 0.4, 0.6, 1.0]
    iv_labels = ["low (<30%)", "med (30-40%)", "high (40-60%)", "very_high (>60%)"]

    for (i, iv_label) in enumerate(iv_labels)
        if i == 1
            mask_i = df_filtered.IV .<= iv_bins[i+1]
        elseif i == length(iv_bins) - 1
            mask_i = df_filtered.IV .> iv_bins[i]
        else
            mask_i = (df_filtered.IV .> iv_bins[i]) .& (df_filtered.IV .<= iv_bins[i+1])
        end
        iv_subset = df_filtered[mask_i, :]
        if nrow(iv_subset) > 0
            iv_errors = (iv_subset.calibrated_price .- iv_subset.Mid) ./ max.(iv_subset.Mid, 1e-10)
            iv_rmspe = sqrt(mean(iv_errors.^2)) * 100
            println("  ", lpad(iv_label, 15), ": RMSPE = ", round(iv_rmspe, digits=2), "% (n=", nrow(iv_subset), ")")
        end
    end

    # Write calibrated output
    println("\n" * "-"^80)
    println("Writing calibrated output...")
    CSV.write("/home/skoranne/PROGRAMS/QUANTUM/CB_BZ_calibrated.csv", df_filtered)

    println("\nDone!")
    return rmspe, best_dte_scales, best_iv_scales
end

# Run main function
main()
