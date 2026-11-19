using CSV, DataFrames, Statistics, Printf
using Optim

println("Loading CB_BZ option chain data...")
df = DataFrame(CSV.File("CB_BZ_option_chain_2026-11-16.csv"))

mask = (df.DTE .> 0) .& (df.IV .> 0.01) .& (df.Mid .> 0)
df_f = df[mask, :]

println("Valid contracts: ", nrow(df_f))
df_f[!, :moneyness] = df_f.Strike ./ df_f.Spot_Price

function get_dte_idx(dte)
    if dte <= 14; return 1
    elseif dte <= 30; return 2
    elseif dte <= 60; return 3
    elseif dte <= 180; return 4
    else; return 5
    end
end
function get_iv_idx(iv)
    if iv <= 0.25; return 1
    elseif iv <= 0.35; return 2
    elseif iv <= 0.50; return 3
    elseif iv <= 0.70; return 4
    else; return 5
    end
end
function get_m_idx(m)
    if m <= 0.90; return 1
    elseif m <= 0.97; return 2
    elseif m <= 1.03; return 3
    elseif m <= 1.10; return 4
    else; return 5
    end
end

function rmspe_loss(params)
    p = clamp.(params, 0.1, 5.0)
    g = p[1]  # global scale
    d_scales = p[2:6]
    i_scales = p[7:11]
    m_scales = p[12:16]
    total_sq = 0.0
    n = 0
    @inbounds for i in 1:nrow(df_f)
        mid = df_f.Mid[i]
        if mid <= 0; continue; end
        theo = df_f.Theo[i]
        ds = d_scales[get_dte_idx(df_f.DTE[i])]
        ivs = i_scales[get_iv_idx(df_f.IV[i])]
        ms = m_scales[get_m_idx(df_f.moneyness[i])]
        cal = theo * g * ds * ivs * ms
        rel = (cal - mid) / mid
        total_sq += rel * rel
        n += 1
    end
    return sqrt(total_sq / max(n, 1)) * 100.0
end

x0 = [1.0, 0.8462,0.8544,1.0619,0.9925,0.8373, 1.6402,1.3169,1.082,0.7658,0.6039, 1.2437,1.1709,1.0441,0.8585,0.9245]

println("\nOptimizing with global scale + finer bins...")
res = optimize(rmspe_loss, x0, NelderMead(), Optim.Options(iterations=6000, show_trace=false, time_limit=180.0))

println("\n" * "="^60)
println("OPTIMIZATION RESULT (global + finer bins)")
println("="^60)
println("Converged: ", Optim.converged(res))
println("Minimum RMSPE: ", round(Optim.minimum(res), digits=4), "%")
best = clamp.(Optim.minimizer(res), 0.1, 5.0)
println("Global: ", round(best[1], digits=4))
println("DTE: ", round.(best[2:6], digits=4))
println("IV:  ", round.(best[7:11], digits=4))
println("MNY: ", round.(best[12:16], digits=4))

df_f[!, :cal_optim3] = [df_f.Theo[i] * best[1] * best[1+get_dte_idx(df_f.DTE[i])] * best[6 + get_iv_idx(df_f.IV[i])] * best[11 + get_m_idx(df_f.moneyness[i])] for i in 1:nrow(df_f)]
maskv = df_f.Mid .> 0
errs = (df_f.cal_optim3[maskv] .- df_f.Mid[maskv]) ./ max.(df_f.Mid[maskv], 1e-10)
final_rmspe = sqrt(mean(errs.^2)) * 100
println("\nFinal RMSPE: ", round(final_rmspe, digits=2), "%")
println("Target <20%: ", final_rmspe < 20 ? "ACHIEVED" : "NOT ACHIEVED")

CSV.write("CB_BZ_calibrated_optim_v3.csv", df_f)
println("Wrote CB_BZ_calibrated_optim_v3.csv")
