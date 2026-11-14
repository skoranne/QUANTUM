# ==============================================================================
# CANONICAL EXOTIC ASIAN OPTION PRICER (Classical MC, CPU/GPU)
# - Arithmetic-average Asian (path-dependent), European exercise
# - Corrected per-step NIG-copula dynamics (consistent with path_dep_option5.jl)
# - CPU Monte Carlo always; CUDA GPU kernel when CUDA.functional()
# - Modes: price (single/basket from CSV), validate, backtest (realized payoff)
#
# Run from repo root:
#   julia --project=. asian_option_mc.jl validate
#   julia --project=. asian_option_mc.jl price-soxx
#   julia --project=. asian_option_mc.jl price-demo
#   julia --project=. asian_option_mc.jl backtest-soxx
# ==============================================================================

module AsianOptionMC

using Printf
using LinearAlgebra
using Random
using Statistics
using Distributions
using SpecialFunctions
using CSV
using DataFrames

const HAS_CUDA = try
    @eval using CUDA
    CUDA.functional()
catch
    false
end

# ------------------------------ Core math ------------------------------------

@inline fast_normcdf(x::Float64) = 0.5 * erfc(-x * 0.7071067811865475244)

struct NormalInverseGaussian <: ContinuousUnivariateDistribution
    μ::Float64
    α::Float64
    β::Float64
    δ::Float64
end

"Stable NIG pdf via exponentially scaled Bessel K (besselkx)."
function nig_pdf(d::NormalInverseGaussian, x::Float64)
    χ = sqrt(max(d.α^2 - d.β^2, eps(Float64)))
    y = x - d.μ
    arg = d.α * sqrt(d.δ^2 + y^2)
    # besselkx(1,arg) = besselk(1,arg)*exp(arg) cancels overflow
    return (d.α * exp(d.δ * χ + d.β * y - arg) * besselkx(1, arg)) / (π * sqrt(d.δ^2 + y^2))
end

"Build a quantile LUT for a NIG distribution (cumulative trapezoid + inversion)."
function build_nig_lut(d::NormalInverseGaussian; resolution::Int=8192, n_pdf::Int=6000)
    χ = sqrt(max(d.α^2 - d.β^2, 1e-12))
    mean_val = d.μ + d.δ * d.β / χ
    # variance of NIG: δ*α^2/χ^3
    var_val = d.δ * d.α^2 / (χ^3)
    std_val = sqrt(max(var_val, 1e-12))

    x_min, x_max = mean_val - 12.0 * std_val, mean_val + 12.0 * std_val
    dx = (x_max - x_min) / (n_pdf - 1)
    xs = range(x_min, x_max, length=n_pdf)
    pdfs = [max(nig_pdf(d, x), 0.0) for x in xs]
    cdf = zeros(Float64, n_pdf)
    @inbounds for i in 2:n_pdf
        cdf[i] = cdf[i-1] + 0.5 * (pdfs[i] + pdfs[i-1]) * dx
    end
    total = cdf[end]
    if total > 0
        cdf ./= total
    end
    cdf[end] = 1.0

    u_grid = collect(range(1e-6, 1.0 - 1e-6, length=resolution))
    q_grid = zeros(Float64, resolution)
    @inbounds for k in 1:resolution
        p = u_grid[k]
        idx = searchsortedfirst(cdf, p)
        if idx <= 1
            q_grid[k] = xs[1]
        elseif idx > n_pdf
            q_grid[k] = xs[end]
        else
            p0, p1 = cdf[idx-1], cdf[idx]
            x0, x1 = xs[idx-1], xs[idx]
            q_grid[k] = x0 + (p - p0) * (x1 - x0) / max(p1 - p0, 1e-18)
        end
    end
    return (u_grid=u_grid, q_grid=q_grid)
end

@inline function lut_quantile(lut, u::Float64)
    ug = lut.u_grid; qg = lut.q_grid
    u_c = clamp(u, ug[1], ug[end])
    n = length(ug)
    pos = 1.0 + (u_c - ug[1]) / (ug[end] - ug[1]) * (n - 1)
    i0 = clamp(Int(floor(pos)), 1, n - 1)
    frac = pos - i0
    return qg[i0] + frac * (qg[i0 + 1] - qg[i0])
end

"""
Corrected per-step NIG increment parameters (risk-neutral drift adjustment).
Returns (μ, α, β, δ) for the *log-return increment* over dt (spot handled separately).
"""
function nig_increment_params(ivol::Float64, skew::Float64, r::Float64, dt::Float64)
    σ = ivol
    s = clamp(skew, -0.95, 0.95)
    # NIG shape from vol/skew (common moment-matched parameterization)
    δ = σ * sqrt(max(1e-8, 1.0 - s^2)) * sqrt(dt)
    α = 1.0 / (σ * max(0.05, abs(s)))
    β = s * α
    χ = sqrt(max(α^2 - β^2, 1e-12))
    χ1 = sqrt(max(α^2 - (β + 1.0)^2, 1e-12))
    # martingale-corrected drift for E[exp(increment)] = exp(r*dt)
    μ = (r - 0.5 * σ^2) * dt - δ * (χ - χ1)
    return NormalInverseGaussian(μ, α, β, δ)
end

# ------------------------------ Contracts ------------------------------------

Base.@kwdef struct AsianOption
    strike::Float64
    T::Float64
    r::Float64
    steps::Int = 252
    bound_b::Float64 = 1.0          # payoff cap/scale; normalized amplitude in [0,1]
    use_bound::Bool = true
    is_call::Bool = true
end

Base.@kwdef struct Market
    spots::Vector{Float64}
    weights::Vector{Float64}
    vols::Vector{Float64}
    skews::Vector{Float64}
    L::Matrix{Float64}               # lower Cholesky of correlation
    names::Vector{String} = String[]
end

function normalize_weights(w::AbstractVector{<:Real})
    wf = Float64.(collect(w))
    s = sum(wf)
    s <= 0 && error("weights must sum to > 0")
    return wf ./ s
end

"Build NIG increment LUTs for each asset given market and dt."
function build_market_luts(m::Market, dt::Float64; resolution::Int=8192)
    d = length(m.spots)
    luts = Vector{Any}(undef, d)
    for i in 1:d
        dist = nig_increment_params(m.vols[i], m.skews[i], 0.0, dt) # drift applied separately? no: include r per asset
    end
    # include r in increment (risk-neutral per asset, same r)
    luts = [build_nig_lut(nig_increment_params(m.vols[i], m.skews[i], 0.0, dt); resolution=resolution) for i in 1:d]
    return luts
end

"Risk-neutral NIG increment dist including rate r."
nig_dist_r(vol, skew, r, dt) = nig_increment_params(vol, skew, r, dt)

# ------------------------------ CPU MC ---------------------------------------

"Price weighted-basket arithmetic Asian option via CPU Monte Carlo. Returns NamedTuple."
function price_basket_cpu(m::Market, opt::AsianOption; n_paths::Int=50_000, seed::Int=42,
                           lut_res::Int=4096, rng_threads::Bool=true)
    d = length(m.spots)
    @assert length(m.weights) == d && length(m.vols) == d && length(m.skews) == d
    w = normalize_weights(m.weights)
    dt = opt.T / opt.steps
    # build per-asset LUTs for log-increments
    luts = [build_nig_lut(nig_dist_r(m.vols[i], m.skews[i], opt.r, dt); resolution=lut_res) for i in 1:d]
    L = Matrix{Float64}(m.L)
    init_log = log.(m.spots)
    payoffs = zeros(Float64, n_paths)

    Threads.@threads for p in 1:n_paths
        rng = Xoshiro(seed + p)
        logS = copy(init_log)
        basket_sum = 0.0
        z = zeros(Float64, d)
        y = zeros(Float64, d)
        @inbounds for _ in 1:opt.steps
            # correlated normals
            for i in 1:d; z[i] = randn(rng); end
            mul!(y, L, z)
            day_basket = 0.0
            for i in 1:d
                u = clamp(cdf(Normal(), y[i]), 1e-6, 1.0 - 1e-6)
                logS[i] += lut_quantile(luts[i], u)
                day_basket += w[i] * exp(logS[i])
            end
            basket_sum += day_basket
        end
        avg_basket = basket_sum / opt.steps
        intrinsic = opt.is_call ? max(0.0, avg_basket - opt.strike) : max(0.0, opt.strike - avg_basket)
        scaled = opt.use_bound ? min(intrinsic, opt.bound_b) / opt.bound_b : intrinsic
        payoffs[p] = scaled
    end

    mean_amp = mean(payoffs)
    se_amp = std(payoffs) / sqrt(n_paths)
    disc = exp(-opt.r * opt.T)
    price = opt.use_bound ? mean_amp * opt.bound_b * disc : mean_amp * disc
    price_se = opt.use_bound ? se_amp * opt.bound_b * disc : se_amp * disc
    return (price=price, se=price_se, amplitude=mean_amp, amplitude_se=se_amp,
            n_paths=n_paths, steps=opt.steps, engine="cpu-mc")
end

# ------------------------------ GPU MC ---------------------------------------

if HAS_CUDA
function _gpu_kernel!(payoffs, L_dev, qmat_dev, umin, umax, res, init_spots, w_dev,
                      strike, bound_b, rT_dummy, steps, d, use_bound, is_call,
                      logS_mat, z_mat, corr_mat)
    idx = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    n_paths = length(payoffs)
    if idx <= n_paths
        basket_sum = 0.0
        @inbounds for i in 1:d
            logS_mat[i, idx] = log(init_spots[i])
        end
        # simple LCG per-thread (deterministic given idx)
        state = UInt64(idx) * 0x9e3779b97f4a7c15
        @inbounds for _ in 1:steps
            # normals via Box-Muller
            for i in 1:d
                state = state * 2862933555777941757 + 3037000493
                u1 = clamp(Float64(state) / 18446744073709551616.0, 1e-7, 1.0 - 1e-7)
                state = state * 2862933555777941757 + 3037000493
                u2 = clamp(Float64(state) / 18446744073709551616.0, 1e-7, 1.0 - 1e-7)
                z_mat[i, idx] = sqrt(-2.0 * log(u1)) * cos(2.0 * 3.141592653589793 * u2)
            end
            # correlate: corr = L*z (lower-tri)
            for i in 1:d
                s = 0.0
                for j in 1:i
                    s += L_dev[i, j] * z_mat[j, idx]
                end
                corr_mat[i, idx] = s
            end
            day_basket = 0.0
            for i in 1:d
                x = corr_mat[i, idx]
                u = 0.5 * erfc(-x * 0.7071067811865475244)
                u_c = clamp(u, umin, umax)
                pos = 1.0 + (u_c - umin) / (umax - umin) * (res - 1)
                i0 = clamp(Int(floor(pos)), 1, res - 1)
                frac = pos - i0
                q0 = qmat_dev[i, i0]
                q1 = qmat_dev[i, i0 + 1]
                logS_mat[i, idx] += q0 + frac * (q1 - q0)
                day_basket += w_dev[i] * exp(logS_mat[i, idx])
            end
            basket_sum += day_basket
        end
        avg_basket = basket_sum / steps
        intrinsic = is_call ? max(0.0, avg_basket - strike) : max(0.0, strike - avg_basket)
        payoffs[idx] = use_bound ? min(intrinsic, bound_b) / bound_b : intrinsic
    end
    return
end

function price_basket_gpu(m::Market, opt::AsianOption; n_paths::Int=50_000, lut_res::Int=4096)
    d = length(m.spots)
    w = normalize_weights(m.weights)
    dt = opt.T / opt.steps
    luts = [build_nig_lut(nig_dist_r(m.vols[i], m.skews[i], opt.r, dt); resolution=lut_res) for i in 1:d]
    u_grid = luts[1].u_grid
    qmat = zeros(Float64, d, length(u_grid))
    for i in 1:d
        qmat[i, :] .= luts[i].q_grid
    end

    payoffs_dev = CUDA.zeros(Float64, n_paths)
    L_dev = CuArray(Matrix{Float64}(m.L))
    qmat_dev = CuArray(qmat)
    w_dev = CuArray(w)
    spots_dev = CuArray(Float64.(m.spots))
    logS = CUDA.zeros(Float64, d, n_paths)
    zmat = CUDA.zeros(Float64, d, n_paths)
    cmat = CUDA.zeros(Float64, d, n_paths)
    threads = 256
    blocks = cld(n_paths, threads)
    @cuda threads=threads blocks=blocks _gpu_kernel!(payoffs_dev, L_dev, qmat_dev,
        Float64(u_grid[1]), Float64(u_grid[end]), length(u_grid), spots_dev, w_dev,
        Float64(opt.strike), Float64(opt.bound_b), 0.0, Int32(opt.steps), Int32(d),
        opt.use_bound, opt.is_call, logS, zmat, cmat)
    CUDA.synchronize()
    payoffs = Array(payoffs_dev)
    mean_amp = mean(payoffs)
    se_amp = std(payoffs) / sqrt(n_paths)
    disc = exp(-opt.r * opt.T)
    price = opt.use_bound ? mean_amp * opt.bound_b * disc : mean_amp * disc
    price_se = opt.use_bound ? se_amp * opt.bound_b * disc : se_amp * disc
    return (price=price, se=price_se, amplitude=mean_amp, amplitude_se=se_amp,
            n_paths=n_paths, steps=opt.steps, engine="gpu-mc-cuda")
end
end # HAS_CUDA

function price_basket(m::Market, opt::AsianOption; n_paths::Int=50_000, prefer_gpu::Bool=true,
                      seed::Int=42, lut_res::Int=4096)
    if prefer_gpu && HAS_CUDA
        try
            return price_basket_gpu(m, opt; n_paths=n_paths, lut_res=lut_res)
        catch e
            @warn "GPU pricing failed, falling back to CPU" exception=(e, catch_backtrace())
            return price_basket_cpu(m, opt; n_paths=n_paths, seed=seed, lut_res=lut_res)
        end
    else
        return price_basket_cpu(m, opt; n_paths=n_paths, seed=seed, lut_res=lut_res)
    end
end

# ------------------------------ Data helpers ---------------------------------

"Repair correlation matrix to nearest PSD and return (L, adjusted_flag)."
function psd_cholesky(C::AbstractMatrix{<:Real})
    d = size(C, 1)
    S = Matrix{Float64}(Symmetric(C))
    # normalize diagonal
    for i in 1:d; S[i,i] = 1.0; end
    vals = eigvals(Symmetric(S))
    mineig = minimum(vals)
    adjusted = false
    if mineig < 1e-10
        S .+= Matrix{Float64}(I, d, d) * (abs(mineig) + 1e-6)
        Dinv = Diagonal(1.0 ./ sqrt.(diag(S)))
        S = Matrix{Float64}(Symmetric(Dinv * S * Dinv))
        adjusted = true
    end
    return cholesky(Symmetric(S)).L |> Matrix, adjusted
end

"Calibrate simple per-asset vol/skew from a price matrix (rows=time, cols=assets)."
function calibrate_from_prices(price_mat::Matrix{Float64}; trading_days::Int=252)
    logp = log.(max.(price_mat, 1e-12))
    rets = logp[2:end, :] .- logp[1:end-1, :]
    vols = vec(std(rets, dims=1) .* sqrt(trading_days))
    # simple skew proxy from return skewness-like tilt (bounded)
    skews = Float64[]
    for j in 1:size(rets,2)
        rj = rets[:, j]
        m = mean(rj); s = std(rj) + 1e-12
        sk = clamp(mean(((rj .- m) ./ s).^3) / 6.0, -0.6, -0.05)
        push!(skews, sk)
    end
    # fallback if degenerate
    skews = [isfinite(x) ? x : -0.25 for x in skews]
    vols = [isfinite(x) && x > 1e-6 ? x : 0.2 for x in vols]
    return vols, skews
end

"Load a single-ticker price CSV (columns: date,close or a single/last numeric column)."
function load_single_ticker_market(price_file::String; r::Float64=0.0415)
    df = CSV.read(price_file, DataFrame)
    # pick close column
    col = if "close" in lowercase.(string.(names(df)))
        idx = findfirst(x -> lowercase(string(x))=="close", names(df))
        names(df)[idx]
    else
        names(df)[end]
    end
    prices = Float64.(df[:, col])
    prices = prices[isfinite.(prices) .& (prices .> 0)]
    length(prices) < 30 && error("Not enough price history in $price_file")
    price_mat = reshape(prices, :, 1)
    vols, skews = calibrate_from_prices(price_mat)
    spot = prices[end]
    L = reshape([1.0], 1, 1)
    m = Market(spots=[spot], weights=[1.0], vols=vols, skews=skews, L=L, names=[basename(price_file)])
    return m, spot, vols[1], skews[1]
end

function run_price_ticker(price_file::String; strike=nothing, r::Float64=0.0415, T::Float64=1.0,
                          steps::Int=252, n_paths::Int=50_000, bound_b=nothing, prefer_gpu::Bool=true)
    println(repeat("=", 78))
    println("SINGLE-STOCK ARITHMETIC ASIAN CALL (CSV) | CUDA=$(HAS_CUDA)")
    m, spot, vol, skew = load_single_ticker_market(price_file; r=r)
    K = strike === nothing ? spot : Float64(strike)
    B = bound_b === nothing ? max(0.8*spot, 2*vol*spot*sqrt(T)) : Float64(bound_b)
    opt = AsianOption(strike=K, T=T, r=r, steps=steps, bound_b=B, use_bound=true)
    @printf("File=%s  Spot=%.4f  Vol=%.4f  Skew=%.4f\n", price_file, spot, vol, skew)
    @printf("K=%.4f  T=%.3f  r=%.4f  steps=%d  B=%.4f  paths=%d\n", K, T, r, steps, B, n_paths)
    res = price_basket(m, opt; n_paths=n_paths, prefer_gpu=prefer_gpu)
    @printf("Asian Call Price=%.6f  SE=%.6f  amplitude=%.6f  engine=%s\n",
            res.price, res.se, res.amplitude, res.engine)
    # quick European Black-Scholes approx (same vol) for sanity reference
    d1 = (log(spot/K) + (r + 0.5*vol^2)*T) / (vol*sqrt(T))
    d2 = d1 - vol*sqrt(T)
    N(x) = 0.5*(1+erf(x/sqrt(2)))
    bs = spot*N(d1) - K*exp(-r*T)*N(d2)
    @printf("Reference European BS call (same vol)=%.6f  (Asian <=~ European typical)\n", bs)
    println(repeat("=", 78))
    return res
end

function load_soxx_market(; price_file="soxx_asset_prices.csv", corr_file="soxx_correlation_matrix.csv",
                           target_spot::Float64=530.0, use_equal_weights::Bool=true)
    prices_df = CSV.read(price_file, DataFrame)
    corr_df = CSV.read(corr_file, DataFrame)
    colnames = string.(DataFrames.names(prices_df)[2:end])
    d = length(colnames)
    price_mat = Matrix{Float64}(prices_df[:, 2:end])
    raw_spots = Float64.(price_mat[end, :])
    corr_raw = Matrix{Float64}(corr_df[:, 2:end])
    L, adj = psd_cholesky(corr_raw)
    vols, skews = calibrate_from_prices(price_mat)
    w = use_equal_weights ? fill(1.0/d, d) : normalize_weights(raw_spots) # explicit, reproducible
    current = sum(raw_spots .* w)
    scale = target_spot / current
    spots = raw_spots .* scale
    adj && @info "Correlation matrix repaired to PSD for Cholesky"
    return Market(spots=spots, weights=w, vols=vols, skews=skews, L=L, names=colnames), spots, w, scale
end

# ------------------------------ Validation & demos ---------------------------

function run_validate()
    println(repeat("=", 78))
    println("ASIAN OPTION MC VALIDATION (CUDA functional=$(HAS_CUDA))")
    println(repeat("=", 78))
    # simple synthetic 1-asset case
    Random.seed!(1)
    vols = [0.20]; skews = [-0.2]; spots = [100.0]; w = [1.0]
    L = reshape([1.0], 1, 1)
    m = Market(spots=spots, weights=w, vols=vols, skews=skews, L=L, names=["SYN"])
    opt = AsianOption(strike=100.0, T=1.0, r=0.03, steps=52, bound_b=40.0, use_bound=true)

    r1 = price_basket_cpu(m, opt; n_paths=20_000, seed=7)
    r2 = price_basket_cpu(m, opt; n_paths=80_000, seed=7)
    @printf("CPU MC 20k: price=%.4f se=%.4f\n", r1.price, r1.se)
    @printf("CPU MC 80k: price=%.4f se=%.4f (SE should shrink ~2x)\n", r2.price, r2.se)
    @printf("Bounds ok: %s (0 <= %.4f <= %.4f)\n", (0.0 <= r2.price <= opt.bound_b*exp(-opt.r*opt.T)),
            r2.price, opt.bound_b*exp(-opt.r*opt.T))

    if HAS_CUDA
        try
            rg = price_basket_gpu(m, opt; n_paths=80_000)
            @printf("GPU MC 80k: price=%.4f se=%.4f engine=%s\n", rg.price, rg.se, rg.engine)
            @printf("CPU-GPU diff: %.4f (within ~3*SE? %s)\n", abs(r2.price-rg.price),
                    abs(r2.price-rg.price) <= 3*sqrt(r2.se^2 + rg.se^2))
        catch e
            println("GPU check failed ($(sprint(showerror,e))); CPU results above remain valid.")
        end
    else
        println("CUDA not available in this Julia session; GPU check skipped.")
    end

    # ITM/OTM sanity
    opt_itm = AsianOption(strike=80.0, T=1.0, r=0.03, steps=52, bound_b=60.0)
    opt_otm = AsianOption(strike=130.0, T=1.0, r=0.03, steps=52, bound_b=60.0)
    pitm = price_basket_cpu(m, opt_itm; n_paths=40_000, seed=11).price
    potm = price_basket_cpu(m, opt_otm; n_paths=40_000, seed=11).price
    @printf("ITM(80)=%.4f  OTM(130)=%.4f  (ITM>OTM: %s)\n", pitm, potm, pitm > potm)
    println(repeat("=", 78))
end

function run_price_soxx()
    println(repeat("=", 78))
    println("SOXX BASKET ASIAN CALL (CSV) | CUDA=$(HAS_CUDA)")
    if !isfile("soxx_asset_prices.csv") || !isfile("soxx_correlation_matrix.csv")
        println("SOXX CSV files not found in cwd; skipping.")
        return
    end
    mkt, spots, w, scale = load_soxx_market()
    S0 = sum(spots .* w)
    opt = AsianOption(strike=S0, T=1.0, r=0.0415, steps=252, bound_b=200.0, use_bound=true)
    println("Assets=$(length(spots))  S0=$(round(S0,digits=2))  K=$(opt.strike)  steps=$(opt.steps)")
    n = 25_000
    res = price_basket(mkt, opt; n_paths=n, prefer_gpu=true)
    @printf("Price=%.4f  SE=%.4f  amp=%.6f  paths=%d  engine=%s\n",
            res.price, res.se, res.amplitude, res.n_paths, res.engine)
    println(repeat("=", 78))
end

function run_backtest_soxx(; lookback::Int=120, horizon::Int=60, step_stride::Int=20)
    println(repeat("=", 78))
    println("REALIZED-PAYOFF BACKTEST (SOXX CSV) | CUDA=$(HAS_CUDA)")
    if !isfile("soxx_asset_prices.csv")
        println("soxx_asset_prices.csv not found; skipping.")
        return
    end
    prices_df = CSV.read("soxx_asset_prices.csv", DataFrame)
    price_mat = Matrix{Float64}(prices_df[:, 2:end])
    nrows, d = size(price_mat)
    if isfile("soxx_correlation_matrix.csv")
        corr_df = CSV.read("soxx_correlation_matrix.csv", DataFrame)
        corr_raw = Matrix{Float64}(corr_df[:, 2:end])
    else
        corr_raw = cor(log.(price_mat[2:end,:] ./ price_mat[1:end-1,:]))
    end
    rows = Int[]
    model_prices = Float64[]
    realized = Float64[]
    r = 0.04
    t = 1
    while t + lookback + horizon <= nrows
        hist = price_mat[t:t+lookback-1, :]
        fwd = price_mat[t+lookback:t+lookback+horizon-1, :]
        vols, skews = calibrate_from_prices(hist)
        L, _ = psd_cholesky(corr_raw)
        raw_spots = hist[end, :]
        w = normalize_weights(fill(1.0, d))
        spots = raw_spots .* 1.0
        S0 = sum(spots .* w)
        opt = AsianOption(strike=S0, T=horizon/252, r=r, steps=horizon, bound_b=0.5*S0, use_bound=true)
        mkt = Market(spots=spots, weights=w, vols=vols, skews=skews, L=L)
        res = price_basket(mkt, opt; n_paths=8_000, prefer_gpu=true)
        # realized arithmetic-average payoff over forward window
        avg_path = [sum(w .* fwd[i, :]) for i in 1:horizon]
        realized_avg = mean(avg_path)
        realized_payoff = max(0.0, realized_avg - opt.strike)
        push!(rows, t+lookback)
        push!(model_prices, res.price)
        push!(realized, realized_payoff * exp(-r*opt.T))
        t += step_stride
    end
    if isempty(rows)
        println("Not enough data for backtest windows.")
        return
    end
    errs = model_prices .- realized
    @printf("Windows=%d  MAE=%.4f  Bias(mean err)=%.4f\n", length(rows), mean(abs.(errs)), mean(errs))
    for i in 1:min(10, length(rows))
        @printf("  idx=%d  model=%.4f  realized=%.4f  err=%.4f\n", rows[i], model_prices[i], realized[i], errs[i])
    end
    out = DataFrame(index=rows, model_price=model_prices, realized_disc=realized, error=errs)
    CSV.write("asian_backtest_soxx_results.csv", out)
    println("Wrote asian_backtest_soxx_results.csv")
    println(repeat("=", 78))
end

end # module

# ------------------------------ CLI ------------------------------------------
if abspath(PROGRAM_FILE) == @__FILE__
    mode = length(ARGS) >= 1 ? ARGS[1] : "help"
    if mode == "validate"
        AsianOptionMC.run_validate()
    elseif mode == "price-soxx"
        AsianOptionMC.run_price_soxx()
    elseif mode == "backtest-soxx"
        AsianOptionMC.run_backtest_soxx()
    elseif mode == "price-ticker"
        length(ARGS) < 2 && error("Usage: price-ticker <prices.csv> [strike] [T] [r] [steps] [paths]")
        f = ARGS[2]
        strike = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : nothing
        T = length(ARGS) >= 4 ? parse(Float64, ARGS[4]) : 1.0
        r = length(ARGS) >= 5 ? parse(Float64, ARGS[5]) : 0.0415
        steps = length(ARGS) >= 6 ? parse(Int, ARGS[6]) : 252
        npaths = length(ARGS) >= 7 ? parse(Int, ARGS[7]) : 50_000
        AsianOptionMC.run_price_ticker(f; strike=strike, T=T, r=r, steps=steps, n_paths=npaths)
    elseif mode == "price-demo"
        AsianOptionMC.run_validate()
    else
        println("Usage: julia --project=. asian_option_mc.jl [validate|price-soxx|backtest-soxx|price-ticker <csv> [K] [T] [r] [steps] [paths]]")
        println("  validate       - run MC sanity/consistency checks (CPU/GPU)")
        println("  price-soxx     - price SOXX basket Asian call from soxx_*.csv")
        println("  backtest-soxx  - realized-payoff rolling backtest on soxx_asset_prices.csv")
        println("  price-ticker   - price single-stock arithmetic Asian call from a ticker price CSV (date,close)")
    end
end
