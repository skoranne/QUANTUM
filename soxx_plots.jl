# ==============================================================================
# SOXX HISTORICAL ANALYSIS & TERMINAL VISUALIZATION WITH UNICODEPLOTS
# Evaluates Last 400 Days, Checks Threshold Exceedance, and Projects Next 100 Days
# ==============================================================================

module SOXXPriceVisualizationAndProjection

using Printf
using LinearAlgebra
using Random
using Statistics
using CSV
using DataFrames
using UnicodePlots

function analyze_and_project_soxx()
    println(repeat("=", 100))
    println(" SOXX HISTORICAL ANALYSIS & 100-DAY MONTE CARLO PRICE PROJECTION")
    println(repeat("=", 100))

    price_file = "soxx_asset_prices.csv"
    if !isfile(price_file)
        error("Missing '$price_file'! Please ensure your historical data file is present.")
    end

    df = CSV.read(price_file, DataFrame)
    ticker_names = string.(names(df)[2:end])
    d_assets = length(ticker_names)

    # Extract historical closing prices
    raw_prices = Matrix{Float64}(df[:, 2:end])
    n_total_days = size(raw_prices, 1)
    
    if n_total_days < 400
        println("  [Notice] Total available days ($n_total_days) is less than 400. Using all available history.")
    end
    
    # Slice to last 400 days (or maximum available)
    lookback_window = min(400, n_total_days)
    historical_slice = raw_prices[end - lookback_window + 1:end, :]

    # Apply the same cap-weighting approach from our pricing model
    weights = rand(d_assets)
    weights ./= sum(weights)

    # Compute historical daily index series for the basket
    basket_index_history = historical_slice * weights

    # Scale to match current market index baseline (~$530)
    current_market_anchor = 530.0
    scaling_multiplier = current_market_anchor / basket_index_history[end]
    basket_index_history .*= scaling_multiplier

    # Define Strike (K) and Option Price to test against
    strike = 530.0
    model_option_price = 71.09 # From our recent QAE/MC valuation
    threshold_price = strike + model_option_price # ~$601.09

    # Filter for the last 200 days within our window
    recent_window = min(200, length(basket_index_history))
    recent_indices = (length(basket_index_history) - recent_window + 1):length(basket_index_history)
    recent_history = basket_index_history[recent_indices]
    
    # Check condition: Where was SOXX higher than Strike + Option Price?
    exceedance_mask = recent_history .> threshold_price
    exceed_count = sum(exceedance_mask)

    @printf("  - Analysis Window Evaluated : Last %d Days\n", recent_window)
    @printf("  - Target Threshold (K + Opt): \$%.2f\n", threshold_price)
    @printf("  - Days Exceeding Threshold  : %d out of %d days\n", exceed_count, recent_window)

    # 1. Render Unicode Plot for the Historical Last 200 Days
    println("\n [1] Unicode Terminal Plot: SOXX Basket vs. Strike + Option Threshold")
    plt_hist = lineplot(
        1:recent_window, 
        recent_history, 
        title="SOXX Basket History (Last $recent_window Days)", 
        name="Basket Index ($)",
        xlabel="Trading Days (Recent History)",
        ylabel="Index Value ($)",
        height=12,
        width=60
    )
    # Draw horizontal line showing our threshold
    hline!(plt_hist, threshold_price, color=:red, name="Strike + Option (\$601.09)")
    print(plt_hist)

    # 2. Model and Project Price Movement for the Next 100 Days via Geometric Brownian Motion
    println("\n\n [2] Simulating 100-Day Forward Price Paths (Monte Carlo Projection)...")
    
    # Compute daily drift and volatility parameters from historical log returns
    log_rets = log.(basket_index_history[2:end] ./ basket_index_history[1:end-1])
    daily_mu = mean(log_rets)
    daily_sigma = std(log_rets)

    future_steps = 100
    n_sim_paths = 500
    
    Random.seed!(42)
    simulated_paths = zeros(future_steps + 1, n_sim_paths)
    simulated_paths[1, :] .= basket_index_history[end]

    dt = 1.0
    for p in 1:n_sim_paths
        for s in 1:future_steps
            z = randn()
            simulated_paths[s+1, p] = simulated_paths[s, p] * exp((daily_mu - 0.5 * daily_sigma^2) * dt + daily_sigma * sqrt(dt) * z)
        end
    end

    # Calculate median projection path and confidence bounds
    median_path = [median(simulated_paths[s, :]) for s in 1:(future_steps + 1)]
    upper_path  = [quantile(simulated_paths[s, :], 0.95) for s in 1:(future_steps + 1)]
    lower_path  = [quantile(simulated_paths[s, :], 0.05) for s in 1:(future_steps + 1)]

    # 3. Render Unicode Plot for the 100-Day Forward Projections
    println("\n [3] Unicode Terminal Plot: 100-Day Forward Price Projections (Median & 90% Confidence Bounds)")
    plt_proj = lineplot(
        1:(future_steps + 1),
        median_path,
        title="100-Day SOXX Forward Price Projection",
        name="Median Projection",
        xlabel="Projected Trading Days Ahead",
        ylabel="Projected Index Value ($)",
        height=12,
        width=60
    )
    lineplot!(plt_proj, 1:(future_steps + 1), upper_path, color=:green, name="95th Percentile")
    lineplot!(plt_proj, 1:(future_steps + 1), lower_path, color=:blue, name="5th Percentile")
    hline!(plt_proj, threshold_price, color=:red, name="Strike + Option Threshold")
    print(plt_proj)
    println(repeat("=", 100))
end

end

SOXXPriceVisualizationAndProjection.analyze_and_project_soxx()
