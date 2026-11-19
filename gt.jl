# Ground Truth using Monte Carlo

using CUDA
using cuTile
import cuTile as CT
using DataFrames
using Distributions
using Random
using Dates
using Statistics

# cuOpt is optional (may be unavailable); guard usage
const HAS_CUOPT = try
    @eval using JuMP
    @eval using cuOpt
    true
catch e
    @warn "cuOpt/MathOpt not available; skipping cuOpt calibration." exception=(e, catch_backtrace())
    false
end

# ==============================================================================
# 1. Synthetic ICE Brent 1st-Line Futures Generator (with Jump-Diffusion & Roll)
# ==============================================================================
function generate_synthetic_brent_history(;
    n_days::Int = 252,
    s0::Float64 = 78.50,
    mu::Float64 = 0.02,
    sigma::Float64 = 0.32,
    jump_intensity::Float64 = 0.03,
    jump_mean::Float64 = 0.0,
    jump_std::Float64 = 0.04,
    roll_day_of_month::Int = 20,
    roll_spread_mean::Float64 = -0.35,
    seed::Int = 2026
)::DataFrame
    rng = MersenneTwister(seed)
    dt = 1.0 / 252.0

    # Generate business days
    dates = Date[]
    curr = Date(2026, 1, 2)
    while length(dates) < n_days
        if dayofweek(curr) <= 5
            push!(dates, curr)
        end
        curr += Day(1)
    end

    prices = zeros(Float64, n_days)
    is_roll = falses(n_days)
    prices[1] = s0

    for t in 2:n_days
        z = randn(rng)
        diffusion = (mu - 0.5 * sigma^2) * dt + sigma * sqrt(dt) * z

        n_jumps = rand(rng, Poisson(jump_intensity))
        jump = n_jumps > 0 ? randn(rng) * jump_std * sqrt(n_jumps) + jump_mean * n_jumps : 0.0

        prices[t] = prices[t - 1] * exp(diffusion + jump)

        # Simulate ICE Brent front-month calendar roll
        if day(dates[t]) >= roll_day_of_month && day(dates[t - 1]) < roll_day_of_month
            prices[t] += roll_spread_mean
            is_roll[t] = true
        end
    end

    return DataFrame(
        date = dates,
        brent_1st_line = round.(prices, digits=2),
        is_roll_day = is_roll
    )
end

# ==============================================================================
# 2. Turnbull-Wakeman / Discrete moment-matching (lognormal) baseline
# ==============================================================================
function turnbull_wakeman_call_continuous(S0::Float64, K::Float64, r::Float64, sigma::Float64, T_start::Float64, T_end::Float64)::Float64
    tau = T_end - T_start
    if tau <= 1e-8
        return exp(-r * T_end) * max(S0 - K, 0.0)
    end
    M1 = S0
    term1 = (2.0 * exp((sigma^2) * tau) - 2.0 * (1.0 + (sigma^2) * tau)) / ((sigma^4) * (tau^2))
    M2 = (S0^2) * exp((sigma^2) * T_start) * term1

    sigma_A = sqrt(log(M2 / (M1^2)) / T_end)
    d1 = (log(M1 / K) + 0.5 * (sigma_A^2) * T_end) / (sigma_A * sqrt(T_end))
    d2 = d1 - sigma_A * sqrt(T_end)

    return exp(-r * T_end) * (M1 * cdf(Normal(), d1) - K * cdf(Normal(), d2))
end

function discrete_aa_lognormal_moment_match_call(S0::Float64, K::Float64, r::Float64, sigma::Float64, T_start::Float64, T_end::Float64, N::Int)::Float64
    # Discrete arithmetic average of N fixings at equally spaced times from T_start+dt..T_end
    if N <= 0
        return exp(-r * T_end) * max(S0 - K, 0.0)
    end
    dt = (T_end - T_start) / N
    M1 = S0
    M2 = 0.0
    for i in 1:N
        ti = T_start + i * dt
        for j in 1:N
            tj = T_start + j * dt
            M2 += S0^2 * exp(sigma^2 * min(ti, tj))
        end
    end
    M2 /= (N^2)
    if M2 <= M1^2 || T_end <= 0
        return exp(-r * T_end) * max(M1 - K, 0.0)
    end
    sigma_A = sqrt(log(M2 / (M1^2)) / T_end)
    d1 = (log(M1 / K) + 0.5 * (sigma_A^2) * T_end) / (sigma_A * sqrt(T_end))
    d2 = d1 - sigma_A * sqrt(T_end)
    return exp(-r * T_end) * (M1 * cdf(Normal(), d1) - K * cdf(Normal(), d2))
end

# ==============================================================================
# 3. Method A: CUDA.jl SIMT Antithetic Kernel (Zero Intermediate Path Memory)
# ==============================================================================
function apo_mc_simt_kernel!(
    payoffs_pos::CuDeviceVector{Float32},
    payoffs_neg::CuDeviceVector{Float32},
    z_fwd::CuDeviceVector{Float32},
    z_steps::CuDeviceMatrix{Float32},
    S0::Float32, K::Float32, discount::Float32,
    drift_fwd::Float32, vol_fwd::Float32,
    drift_dt::Float32, vol_dt::Float32,
    realized_sum::Float32, n_total_days::Int32, n_remaining::Int32,
    is_call::Bool
)
    idx = (blockIdx().x - Int32(1)) * blockDim().x + threadIdx().x
    if idx <= length(payoffs_pos)
        # Forward start step (if unseasoned and days_to_avg_start > 0)
        zf = z_fwd[idx]
        s_pos = S0 * exp(drift_fwd + vol_fwd * zf)
        s_neg = S0 * exp(drift_fwd - vol_fwd * zf)

        sum_pos = realized_sum
        sum_neg = realized_sum

        @inbounds for step in Int32(1):n_remaining
            z = z_steps[idx, step]
            s_pos *= exp(drift_dt + vol_dt * z)
            s_neg *= exp(drift_dt - vol_dt * z)
            sum_pos += s_pos
            sum_neg += s_neg
        end

        avg_pos = sum_pos / Float32(n_total_days)
        avg_neg = sum_neg / Float32(n_total_days)

        p_pos = is_call ? max(avg_pos - K, 0.0f0) : max(K - avg_pos, 0.0f0)
        p_neg = is_call ? max(avg_neg - K, 0.0f0) : max(K - avg_neg, 0.0f0)

        payoffs_pos[idx] = discount * p_pos
        payoffs_neg[idx] = discount * p_neg
    end
    return nothing
end

# ==============================================================================
# 4. Method B: cuTile.jl Tile-Based Payoff & Moment Aggregation Kernel
# ==============================================================================
const TILE_PATHS = 256

function tile_payoff_kernel!(
    avg_prices::CuDeviceVector{Float32},
    payoffs::CuDeviceVector{Float32},
    K::Float32,
    discount::Float32
)
    bid = CT.blockid(1)
    avg_tile = CT.load(avg_prices, (bid,), (TILE_PATHS,))
    # Idiomatic Julia 1-based broadcasting compiled to NVIDIA Tile IR
    payoff_tile = discount .* max.(avg_tile .- K, 0.0f0)
    CT.store(payoffs, (bid,), payoff_tile)
    return nothing
end

function run_cutile_payoff_evaluation(avg_prices::CuVector{Float32}, K::Float64, discount::Float64)
    n = length(avg_prices)
    @assert n % TILE_PATHS == 0 "Pad path count to a multiple of TILE_PATHS (256)"
    payoffs = CUDA.zeros(Float32, n)
    grid = (cld(n, TILE_PATHS), 1, 1)
    CT.launch(grid, tile_payoff_kernel!, avg_prices, payoffs, Float32(K), Float32(discount))
    CUDA.synchronize()
    return payoffs
end

# ==============================================================================
# 5. Method C: cuOpt.jl GPU LP Martingale-Calibrated Path Reweighting (optional)
# ==============================================================================
if HAS_CUOPT
function cuopt_martingale_calibrated_price(
    terminal_avg_prices::Vector{Float64},
    final_step_prices::Vector{Float64},
    S0::Float64,
    K::Float64,
    discount::Float64,
    sigma::Float64,
    T_expiry::Float64
)
    M = length(terminal_avg_prices)
    u = 1.0 / M
    payoffs = discount .* max.(terminal_avg_prices .- K, 0.0)

    model = Model(cuOpt.Optimizer)
    set_silent(model)

    @variable(model, w[1:M] >= 0.0)
    @variable(model, dev_pos[1:M] >= 0.0)
    @variable(model, dev_neg[1:M] >= 0.0)

    @objective(model, Min, sum(dev_pos[i] + dev_neg[i] for i in 1:M))
    @constraint(model, [i=1:M], w[i] - u == dev_pos[i] - dev_neg[i])
    @constraint(model, sum(w) == 1.0)
    @constraint(model, sum(w[i] * final_step_prices[i] for i in 1:M) == S0)
    target_m2 = (S0^2) * exp((sigma^2) * T_expiry)
    @constraint(model, sum(w[i] * (final_step_prices[i]^2) for i in 1:M) == target_m2)

    optimize!(model)

    w_opt = value.(w)
    calibrated_price = sum(w_opt .* payoffs)
    return calibrated_price
end
else
    function cuopt_martingale_calibrated_price(args...)
        return NaN
    end
end

# ==============================================================================
# 6. Unified Ground-Truth Benchmark Runner
# ==============================================================================
function price_apo_gpu_benchmark(;
    S0::Float64,
    K::Float64,
    r::Float64,
    sigma::Float64,
    n_averaging_days::Int,
    days_to_avg_start::Int = 0,
    realized_fixings::Vector{Float64} = Float64[],
    n_paths::Int = 1_048_576, # Power of 2 (2^20) divisible by TILE_PATHS
    use_cuopt_calibration::Bool = true
)
    dt = 1.0 / 252.0
    n_realized = length(realized_fixings)
    n_remaining = n_averaging_days - n_realized
    @assert n_remaining > 0 "APO has already fully averaged."
    @assert iseven(n_paths) "n_paths must be even for antithetic variates."

    half_paths = n_paths ÷ 2
    T_expiry = (days_to_avg_start + n_remaining) * dt
    discount = exp(-r * T_expiry)

    # Allocate GPU arrays via CUDA.jl CURAND
    z_fwd = (days_to_avg_start > 0 && n_realized == 0) ? CUDA.randn(Float32, half_paths) : CUDA.zeros(Float32, half_paths)
    z_steps = CUDA.randn(Float32, half_paths, n_remaining)
    payoffs_pos = CUDA.zeros(Float32, half_paths)
    payoffs_neg = CUDA.zeros(Float32, half_paths)

    t_fwd = (n_realized == 0) ? days_to_avg_start * dt : 0.0
    drift_fwd = Float32(-0.5 * sigma^2 * t_fwd)
    vol_fwd = Float32(sigma * sqrt(t_fwd))
    drift_dt = Float32(-0.5 * sigma^2 * dt)
    vol_dt = Float32(sigma * sqrt(dt))
    realized_sum = Float32(sum(realized_fixings))

    threads = 256
    blocks = cld(half_paths, threads)

    @cuda threads=threads blocks=blocks apo_mc_simt_kernel!(
        payoffs_pos, payoffs_neg, z_fwd, z_steps,
        Float32(S0), Float32(K), Float32(discount),
        drift_fwd, vol_fwd, drift_dt, vol_dt,
        realized_sum, Int32(n_averaging_days), Int32(n_remaining),
        true
    )
    CUDA.synchronize()

    # Combine antithetic pairs on GPU
    pair_means = 0.5f0 .* (payoffs_pos .+ payoffs_neg)
    mc_price = Float64(mean(pair_means))
    mc_std_err = Float64(std(pair_means)) / sqrt(half_paths)

    # Optional cuOpt Martingale-Calibrated Benchmark on a 25,000-path stratified subset
    cuopt_price = NaN
    if HAS_CUOPT && use_cuopt_calibration && n_realized == 0 && days_to_avg_start == 0
        sub_M = 25_000
        rng_sub = MersenneTwister(99)
        z_sub = randn(rng_sub, sub_M, n_remaining)
        log_inc = (-0.5 * sigma^2 * dt) .+ (sigma * sqrt(dt)) .* z_sub
        paths_sub = S0 .* exp.(cumsum(log_inc, dims=2))
        avg_sub = vec(mean(paths_sub, dims=2))
        final_sub = paths_sub[:, end]
        cuopt_price = cuopt_martingale_calibrated_price(avg_sub, final_sub, S0, K, discount, sigma, T_expiry)
    end
    return (
        mc_price = mc_price,
        mc_std_err = mc_std_err,
        cuopt_calibrated_price = cuopt_price,
        n_remaining_steps = n_remaining
    )
end
function run_validation_suite()
println("=== 1. Generating Synthetic ICE Brent Futures History ===")
hist_df = generate_synthetic_brent_history(n_days=252)
println(first(hist_df, 5), "\n")
test_cases = [
(id="QAE_4STEP_ITM", S0=80.0, K=75.0, r=0.045, sigma=0.30, N=4,  d_start=0,  realized=Float64[]),
(id="QAE_4STEP_ATM", S0=80.0, K=80.0, r=0.045, sigma=0.30, N=4,  d_start=0,  realized=Float64[]),
(id="QAE_4STEP_OTM", S0=80.0, K=85.0, r=0.045, sigma=0.30, N=4,  d_start=0,  realized=Float64[]),
(id="FWD_21STEP_ATM", S0=80.0, K=80.0, r=0.045, sigma=0.32, N=21, d_start=10, realized=Float64[]),
(
id="SEASONED_5REM_ATM", S0=81.20, K=80.0, r=0.045, sigma=0.32, N=21, d_start=0,
realized=[78.5, 79.0, 78.8, 79.2, 79.5, 80.1, 80.3, 79.9, 80.5, 80.8, 81.0, 80.7, 81.1, 81.4, 81.0, 81.2]
)
]
rows = []
for tc in test_cases
res = price_apo_gpu_benchmark(
S0=tc.S0, K=tc.K, r=tc.r, sigma=tc.sigma,
n_averaging_days=tc.N, days_to_avg_start=tc.d_start,
realized_fixings=tc.realized, n_paths=1_048_576,
use_cuopt_calibration=true
)
tw = isempty(tc.realized) ? discrete_aa_lognormal_moment_match_call(
tc.S0, tc.K, tc.r, tc.sigma, tc.d_start / 252.0, (tc.d_start + tc.N) / 252.0, tc.N
) : NaN
push!(rows, (
case_id = tc.id,
S0 = tc.S0,
K = tc.K,
steps = res.n_remaining_steps,
cuda_mc_price = round(res.mc_price, digits=5),
cuda_mc_stderr = round(res.mc_std_err, digits=6),
cuopt_lp_price = isnan(res.cuopt_calibrated_price) ? missing : round(res.cuopt_calibrated_price, digits=5),
tw_analytical = isnan(tw) ? missing : round(tw, digits=5)
))
end
println("=== 2. GPU Ground-Truth Benchmarks for QAE Validation ===")
println(DataFrame(rows))
end
run_validation_suite()

