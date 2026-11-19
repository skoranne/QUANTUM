using CSV, DataFrames, Statistics, Distributions

function black76_call_fut(F::Float64, K::Float64, T::Float64, r::Float64, sigma::Float64)
    if T <= 0; return max(F-K,0.0); end
    if sigma <= 1e-12; return exp(-r*T)*max(F-K,0.0); end
    d1 = (log(F/K)+0.5*sigma^2*T)/(sigma*sqrt(T))
    d2 = d1-sigma*sqrt(T)
    return exp(-r*T)*(F*cdf(Normal(),d1)-K*cdf(Normal(),d2))
end
function black76_put_fut(F::Float64, K::Float64, T::Float64, r::Float64, sigma::Float64)
    if T <= 0; return max(K-F,0.0); end
    if sigma <= 1e-12; return exp(-r*T)*max(K-F,0.0); end
    d1 = (log(F/K)+0.5*sigma^2*T)/(sigma*sqrt(T))
    d2 = d1-sigma*sqrt(T)
    return exp(-r*T)*(K*cdf(Normal(),d2)-F*cdf(Normal(),d1))
end

function run()
    df = CSV.read("CB_BZ_option_chain_2026-11-16.csv", DataFrame)
    mask = (df.DTE .> 0) .& (df.IV .> 0) .& (df.Spot_Price .> 0)
    dfv = df[mask, :]
    best_r, best_mean = 0.0, Inf
    for r in 0.0:0.001:0.10
        diffs = Float64[]
        for row in eachrow(dfv)
            F=Float64(row.Spot_Price); K=Float64(row.Strike); T=Float64(row.DTE)/365.0
            sigma=Float64(row.IV); theo=Float64(row.Theo)
            if sigma <= 0; continue; end
            computed = string(row.Type)=="C" ? black76_call_fut(F,K,T,r,sigma) : black76_put_fut(F,K,T,r,sigma)
            push!(diffs, abs(computed-theo))
        end
        m = mean(diffs)
        if m < best_mean
            best_mean = m; best_r = r
        end
    end
    println("Best r by mean|diff|: ", best_r, " with mean|diff|=", round(best_mean, sigdigits=6))
end
run()
