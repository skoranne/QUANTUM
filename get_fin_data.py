from datetime import datetime
import numpy as np
import pandas as pd
import yfinance as yf

# Define a representative multi-asset portfolio basket (expandable up to 128 tickers)
tickers = [
    "AAPL",
    "MSFT",
    "GOOGL",
    "AMZN",
    "NVDA",
    "META",
    "TSLA",
    "JPM",
    "V",
    "JNJ",
    "WMT",
    "PG",
    "MA",
    "UNH",
    "HD",
    "DIS",
    "BAC",
    "XOM",
    "PFE",
    "NFLX",
]

# Fetch daily data with multi_level_index=False to ensure flat, clean columns
data = yf.download(
    tickers, start="2024-01-01", end="2026-01-01", multi_level_index=False
)[  # type: ignore
    "Close"
]

# Drop columns with any missing values to maintain matrix integrity
data = data.dropna(axis=1)

# 1. Save raw historical asset closing prices (Rows = Trading Days, Columns = Tickers)
data.to_csv("asset_prices.csv")

# 2. Compute log-returns and generate the empirical Pearson correlation matrix
log_returns = np.log(data / data.shift(1)).dropna()
corr_matrix = log_returns.corr()

# 3. Save correlation matrix to CSV
corr_matrix.to_csv("correlation_matrix.csv")

print(
    f"Successfully exported price matrix {data.shape} and correlation matrix {corr_matrix.shape}."
)
