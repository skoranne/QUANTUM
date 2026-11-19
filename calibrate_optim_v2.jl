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
    p = clamp.(params, 0.2, 4.0)
    d_scales = p[1:5]
    i_scales = p[6:10]
    m_scales = p[11:15]
    total_sq = 0.0
    n = 0
    @inbounds for i in 1:nrow(df_f)
        mid = df_f.Mid[i]
        if mid <= 0; continue; end
        theo = df_f.Theo[i]
        ds = d_scales[get_dte_idx(df_f.DTE[i])]
        ivs = i_scales[get_iv_idx(df_f.IV[i])]
        ms = m_scales[get_m_idx(df_f.moneyness[i])]
        cal = theo * ds * ivs * ms
        rel = (cal - mid) / mid
        total_sq += rel * rel
        n += 1
    end
    return sqrt(total_sq / max(n, 1)) * 100.0
end

x0 = [0.8,0.7,1.0,1.0,0.9, 1.5,1.2,0.95,0.8,0.7, 1.3,1.1,1.0,0.95,0.9]

println("\nOptimizing with Optim.jl (NelderMead, clamped) finer bins...")
res = optimize(rmspe_loss, x0, NelderMead(), Optim.Options(iterations=5000, show_trace=false, time_limit=120.0))

println("\n" * "="^60)
println("OPTIMIZATION RESULT (finer bins)")
println("="^60)
println("Converged: ", Optim.converged(res))
println("Minimum RMSPE: ", round(Optim.minimum(res), digits=4), "%")
best = clamp.(Optim.minimizer(res), 0.2, 4.0)
println("DTE scales (<=14,15-30,31-60,61-180,>180): ", round.(best[1:5], digits=4))
println("IV scales  (<=0.25,..,>0.7): ", round.(best[6:10], digits=4))
println("MNY scales (<=0.9,0.9-0.97,0.97-1.03,1.03-1.1,>1.1): ", round.(best[11:15], digits=4))

df_f[!, :cal_optim2] = [df_f.Theo[i] * best[get_dte_idx(df_f.DTE[i])] * best[5 + get_iv_idx(df_f.IV[i])] * best[10 + get_m_idx(df_f.moneyness[i])] for i in 1:nrow(df_f)]
maskv = df_f.Mid .> 0
errs = (df_f.cal_optim2[maskv] .- df_f.Mid[maskv]) ./ max.(df_f.Mid[maskv], 1e-10)
final_rmspe = sqrt(mean(errs.^2)) * 100
println("\nFinal RMSPE: ", round(final_rmspe, digits=2), "%")
println("Target <20%: ", final_rmspe < 20 ? "ACHIEVED" : "NOT ACHIEVED")

CSV.write("CB_BZ_calibrated_optim_v2.csv", df_f)
println("Wrote CB_BZ_calibrated_optim_v2.csv")
