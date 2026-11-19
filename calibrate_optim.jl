using CSV, DataFrames, Statistics, Printf
using Optim

println("Loading CB_BZ option chain data...")
df = DataFrame(CSV.File("CB_BZ_option_chain_2026-11-16.csv"))

mask = (df.DTE .> 0) .& (df.IV .> 0.01) .& (df.Mid .> 0)
df_f = df[mask, :]

println("Valid contracts: ", nrow(df_f))
df_f[!, :moneyness] = df_f.Strike ./ df_f.Spot_Price

function get_dte_idx(dte)
    if dte <= 30; return 1
    elseif dte <= 90; return 2
    elseif dte <= 180; return 3
    else; return 4
    end
end

function get_iv_idx(iv)
    if iv <= 0.30; return 1
    elseif iv <= 0.40; return 2
    elseif iv <= 0.60; return 3
    else; return 4
    end
end

function get_m_idx(m)
    if m <= 0.95; return 1
    elseif m <= 1.05; return 2
    else; return 3
    end
end

function rmspe_loss(params)
    # clamp params
    p = clamp.(params, 0.3, 3.0)
    d_scales = p[1:4]
    i_scales = p[5:8]
    m_scales = p[9:11]
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

x0 = [0.731, 1.0805, 1.0173, 0.9702, 1.3721, 1.0992, 0.9213, 0.6929, 1.0, 1.0, 1.0]

println("\nOptimizing with Optim.jl (NelderMead, clamped) with moneyness...")
res = optimize(rmspe_loss, x0, NelderMead(), Optim.Options(iterations=3000, show_trace=false))

println("\n" * "="^60)
println("OPTIMIZATION RESULT (with moneyness, NelderMead)")
println("="^60)
println("Converged: ", Optim.converged(res))
println("Minimum RMSPE: ", round(Optim.minimum(res), digits=4), "%")
best = clamp.(Optim.minimizer(res), 0.3, 3.0)
println("DTE scales:  ", round.(best[1:4], digits=4))
println("IV scales:   ", round.(best[5:8], digits=4))
println("MNY scales:  ", round.(best[9:11], digits=4), " (OTM<=0.95, ATM 0.95-1.05, ITM>1.05)")

df_f[!, :cal_optim] = [df_f.Theo[i] * best[get_dte_idx(df_f.DTE[i])] * best[4 + get_iv_idx(df_f.IV[i])] * best[8 + get_m_idx(df_f.moneyness[i])] for i in 1:nrow(df_f)]
maskv = df_f.Mid .> 0
errs = (df_f.cal_optim[maskv] .- df_f.Mid[maskv]) ./ max.(df_f.Mid[maskv], 1e-10)
final_rmspe = sqrt(mean(errs.^2)) * 100
println("\nFinal RMSPE: ", round(final_rmspe, digits=2), "%")
println("Target <20%: ", final_rmspe < 20 ? "ACHIEVED" : "NOT ACHIEVED")

CSV.write("CB_BZ_calibrated_optim_mny.csv", df_f)
println("Wrote CB_BZ_calibrated_optim_mny.csv")
