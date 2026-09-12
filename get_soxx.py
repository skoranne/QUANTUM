from datetime import datetime
import numpy as np
import pandas as pd
import yfinance as yf

# Complete list of major tickers comprising the SOXX semiconductor universe
soxx_tickers = [
    "NVDA",
    "MU",
    "AMD",
    "AVGO",
    "MRVL",
    "INTC",
    "AMAT",
    "TSM",
    "LRCX",
    "KLAC",
    "ADI",
    "TXN",
    "MPWR",
    "TER",
    "NXPI",
    "QCOM",
    "ASML",
    "ALAB",
    "MCHP",
    "CRDO",
    "ON",
    "ASX",
    "ENTG",
    "MTSI",
    "UMC",
    "NVMI",
    "RMBS",
    "STM",
    "SWKS",
    "ARM",
]

print(
    f"Fetching historical daily price data for {len(soxx_tickers)} SOXX semiconductor components..."
)

# Fetch 2 years of daily historical data (flat column structure using multi_level_index=False)
data = yf.download(
    soxx_tickers, start="2024-01-01", end="2026-09-09", multi_level_index=False
)[  # type: ignore
    "Close"
]

# Drop any assets that contain missing columns or incomplete historical gaps
data = data.dropna(axis=1)

# 1. Export raw daily asset pricing matrix (Rows = Trading Days, Columns = Tickers)
data.to_csv("soxx_asset_prices.csv")

# 2. Compute log-returns and generate the empirical cross-asset Pearson correlation matrix
log_returns = np.log(data / data.shift(1)).dropna()
corr_matrix = log_returns.corr()

# 3. Export correlation matrix
corr_matrix.to_csv("soxx_correlation_matrix.csv")

print(
    f"Success! Exported price matrix with shape {data.shape} and correlation matrix with shape {corr_matrix.shape}."
)
