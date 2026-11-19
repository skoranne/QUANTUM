using CSV, DataFrames, Statistics, Distributions, Printf

function black76_call_fut(F::Float64, K::Float64, T::Float64, r::Float64, sigma::Float64)
    if T <= 0
        return max(F - K, 0.0)
    end
    if sigma <= 0 || isnan(sigma)
        return exp(-r * T) * max(F - K, 0.0)
    end
    d1 = (log(F / K) + 0.5 * sigma^2 * T) / (sigma * sqrt(T))
    d2 = d1 - sigma * sqrt(T)
    return exp(-r * T) * (F * cdf(Normal(), d1) - K * cdf(Normal(), d2))
end

function black76_put_fut(F::Float64, K::Float64, T::Float64, r::Float64, sigma::Float64)
    if T <= 0
        return max(K - F, 0.0)
    end
    if sigma <= 0 || isnan(sigma)
        return exp(-r * T) * max(K - F, 0.0)
    end
    d1 = (log(F / K) + 0.5 * sigma^2 * T) / (sigma * sqrt(T))
    d2 = d1 - sigma * sqrt(T)
    return exp(-r * T) * (K * cdf(Normal(), d2) - F * cdf(Normal(), d1))
end

function run_validation()
    df = CSV.read("CB_BZ_option_chain_2026-11-16.csv", DataFrame)
    mask = (df.DTE .> 0) .& (df.IV .> 0.0) .& (df.Spot_Price .> 0)
    dfv = df[mask, :]

    r = 0.045
    n_match = 0
    n_close = 0
    max_diff = 0.0
    diffs = Float64[]

    for row in eachrow(dfv)
        F = Float64(row.Spot_Price)
        K = Float64(row.Strike)
        T = Float64(row.DTE) / 365.0
        sigma = Float64(row.IV)
        theo = Float64(row.Theo)
        
        if sigma <= 0 || isnan(sigma)
            continue
        end
        
        computed = if string(row.Type) == "C"
            black76_call_fut(F, K, T, r, sigma)
        else
            black76_put_fut(F, K, T, r, sigma)
        end
        
        diff = abs(computed - theo)
        push!(diffs, diff)
        if diff < 1e-2
            n_close += 1
        end
        if diff > max_diff
            max_diff = diff
        end
        n_match += 1
    end

    println("Contracts checked (DTE>0, IV>0): ", n_match)
    if length(diffs) > 0
        println("Mean |diff|: ", round(mean(diffs), sigdigits=6))
        println("Max |diff|:  ", round(max_diff, sigdigits=6))
        println("Fraction |diff| < 1e-2: ", round(n_close / n_match * 100, digits=2), "%")
    end
end

println("="^70)
println("VANILLA BLACK-76 VALIDATION: Compare computed vs Theo column")
println("="^70)
run_validation()
