using Printf, Statistics, LinearAlgebra, Distributions, CSV, DataFrames

# ============================================================
# IMPROVED CALIBRATION WITH DTE AND IV SCALING
# ============================================================

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
# CALIBRATION STRATEGY
# ============================================================
# Based on analysis:
# - DTE < 30: RMSPE = 80.51% (need large correction)
# - 30 <= DTE < 90: RMSPE = 22.36% (need moderate correction)
# - 90 <= DTE < 180: RMSPE = 15.83% (good, small correction)
# - DTE >= 180: RMSPE = 13.36% (excellent, minimal correction)
# - IV > 60%: RMSPE = 84.99% (need correction)
# - IV 30-40%: RMSPE = 13.89% (excellent)

# DTE-based scaling function
function dte_scale(dte)
    if dte <= 30
        return 1.8 + 0.02 * (30 - dte)  # Scale up for very short DTE
    elseif dte <= 90
        return 1.2 + 0.005 * (90 - dte)  # Moderate scale for short DTE
    elseif dte <= 180
        return 1.0 + 0.002 * (180 - dte)  # Small scale for medium DTE
    else
        return 1.0  # No scaling for long DTE
    end
end

# IV-based scaling function
function iv_scale(iv)
    if iv <= 0.30
        return 1.1  # Low IV - scale up slightly
    elseif iv <= 0.40
        return 0.95  # Optimal IV range - slight downscale
    elseif iv <= 0.60
        return 1.05  # High IV - scale up
    else
        return 0.7 + 0.005 * (iv - 0.60) * 10  # Very high IV - scale down
    end
end

# Combined scaling
function combined_scale(dte, iv)
    d_s = dte_scale(dte)
    i_s = iv_scale(iv)
    return d_s * i_s
end

# Apply calibration
df_filtered[!, :calibrated_price] = df_filtered.Theo .* combined_scale.(df_filtered.DTE, df_filtered.IV)

# Calculate RMSPE
mask_valid = df_filtered.Mid .> 0
errors = (df_filtered.calibrated_price[mask_valid] .- df_filtered.Mid[mask_valid]) ./ max.(df_filtered.Mid[mask_valid], 1e-10)
rmspe = sqrt(mean(errors.^2)) * 100

println("\n" * "="^80)
println("CALIBRATION RESULTS (Initial)")
println("="^80)
println("RMSPE: ", round(rmspe, digits=2), "%")
println("Target: <20%")

# ============================================================
# OPTIMIZATION: Grid search for optimal scales
# ============================================================
println("\n" * "-"^80)
println("Grid search for optimal DTE-based scales...")

function compute_rmspe_for_scales(dte_scales, iv_scales)
    total_sq_error = 0.0
    count = 0
    
    for i in 1:nrow(df_filtered)
        if df_filtered.Mid[i] > 0
            dte = df_filtered.DTE[i]
            iv = df_filtered.IV[i]
            
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
            
            calibrated = df_filtered.Theo[i] * d_scale * i_scale
            rel_error = (calibrated - df_filtered.Mid[i]) / df_filtered.Mid[i]
            total_sq_error += rel_error^2
            count += 1
        end
    end
    
    return sqrt(total_sq_error / count) * 100
end

# Grid search
best_rmspe = Inf
best_dte_scales = [1.0, 1.0, 1.0, 1.0]
best_iv_scales = [1.0, 1.0, 1.0, 1.0]

# Coarse search
for d1 in 1.4:0.2:2.2
    for d2 in 1.1:0.1:1.6
        for d3 in 0.9:0.05:1.2
            for d4 in 0.9:0.05:1.1
                for i1 in 0.9:0.1:1.3
                    for i2 in 0.8:0.1:1.2
                        for i3 in 0.9:0.1:1.3
                            for i4 in 0.5:0.1:0.9
                                dte_scales = [d1, d2, d3, d4]
                                iv_scales = [i1, i2, i3, i4]
                                rmspe = compute_rmspe_for_scales(dte_scales, iv_scales)
                                
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
for d1 in best_dte_scales[1]-0.3:0.05:best_dte_scales[1]+0.3
    for d2 in best_dte_scales[2]-0.2:0.05:best_dte_scales[2]+0.2
        for d3 in best_dte_scales[3]-0.1:0.02:best_dte_scales[3]+0.1
            for d4 in best_dte_scales[4]-0.1:0.02:best_dte_scales[4]+0.1
                for i1 in best_iv_scales[1]-0.2:0.05:best_iv_scales[1]+0.2
                    for i2 in best_iv_scales[2]-0.2:0.05:best_iv_scales[2]+0.2
                        for i3 in best_iv_scales[3]-0.2:0.05:best_iv_scales[3]+0.2
                            for i4 in best_iv_scales[4]-0.2:0.05:best_iv_scales[4]+0.2
                                dte_scales = [d1, d2, d3, d4]
                                iv_scales = [i1, i2, i3, i4]
                                rmspe = compute_rmspe_for_scales(dte_scales, iv_scales)
                                
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
df_filtered[!, :dte_scale] = ifelse.(df_filtered.DTE .<= 30, best_dte_scales[1],
                         ifelse.(df_filtered.DTE .<= 90, best_dte_scales[2],
                         ifelse.(df_filtered.DTE .<= 180, best_dte_scales[3], best_dte_scales[4])))

df_filtered[!, :iv_scale] = ifelse.(df_filtered.IV .<= 0.30, best_iv_scales[1],
                        ifelse.(df_filtered.IV .<= 0.40, best_iv_scales[2],
                        ifelse.(df_filtered.IV .<= 0.60, best_iv_scales[3], best_iv_scales[4])))

df_filtered[!, :calibrated_price] = df_filtered.Theo .* df_filtered.dte_scale .* df_filtered.iv_scale

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
