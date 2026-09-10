using CSV
using DataFrames
using LinearAlgebra

# Ingest saved historical data
prices_df = CSV.read("asset_prices.csv", DataFrame)
corr_df = CSV.read("correlation_matrix.csv", DataFrame)

# Extract vector of latest spot prices
spots = collect(Float64, prices_df[end, 2:end])

# Extract exact correlation matrix and compute Cholesky factor
corr_matrix = Matrix{Float64}(corr_df[:, 2:end])
L = cholesky(Symmetric(corr_matrix)).L

println("Successfully initialized $(length(spots)) real-world assets and correlation structure.")
