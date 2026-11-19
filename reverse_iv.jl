using CSV, DataFrames, Statistics, Distributions, Printf

function black76_call_fut(F::Float64, K::Float64, T::Float64, r::Float64, sigma::Float64)
    if T <= 0; return max(F - K, 0.0); end
    if sigma <= 1e-12; return exp(-r * T) * max(F - K, 0.0); end
    d1 = (log(F / K) + 0.5 * sigma^2 * T) / (sigma * sqrt(T))
    d2 = d1 - sigma * sqrt(T)
    return exp(-r * T) * (F * cdf(Normal(), d1) - K * cdf(Normal(), d2))
end
function black76_put_fut(F::Float64, K::Float64, T::Float64, r::Float64, sigma::Float64)
    if T <= 0; return max(K - F, 0.0); end
    if sigma <= 1e-12; return exp(-r * T) * max(K - F, 0.0); end
    d1 = (log(F / K) + 0.5 * sigma^2 * T) / (sigma * sqrt(T))
    d2 = d1 - sigma * sqrt(T)
    return exp(-r * T) * (K * cdf(Normal(), d2) - F * cdf(Normal(), d1))
end

function imply_vol_bisect(F,K,T,r,price,is_call; tol=1e-8, maxiter=200)
    if price <= 1e-12
        return 0.0
    end
    lo, hi = 1e-6, 5.0
    for i in 1:maxiter
        mid = 0.5*(lo+hi)
        pm = is_call ? black76_call_fut(F,K,T,r,mid) : black76_put_fut(F,K,T,r,mid)
        if pm > price
            hi = mid
        else
            lo = mid
        end
        if hi - lo < tol
            break
        end
    end
    return 0.5*(lo+hi)
end

function run()
    df = CSV.read("CB_BZ_option_chain_2026-11-16.csv", DataFrame)
    mask = (df.DTE .> 0) .& (df.Mid .> 0) .& (df.Spot_Price .> 0)
    dfv = df[mask, :]
    
    r = 0.0
    count = 0
    for row in eachrow(dfv)
        if count >= 15; break; end
        F = Float64(row.Spot_Price); K = Float64(row.Strike); T = Float64(row.DTE)/365.0
        theo = Float64(row.Theo); iv_file = Float64(row.IV)
        is_call = string(row.Type) == "C"
        iv_back = imply_vol_bisect(F,K,T,r,theo,is_call)
        @printf("%s %s DTE=%d K=%.1f theo=%.4f iv_file=%.4f iv_from_theo(r=0)=%.4f\n",
            row.Symbol, row.Type, row.DTE, K, theo, iv_file, iv_back)
        count += 1
    end
end
run()
