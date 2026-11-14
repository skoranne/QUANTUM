#!/usr/bin/env python3
"""Fetch historical prices for a US ticker via yfinance and write a simple CSV
(date,close) plus print current spot/annualized vol estimate for use with the
Asian option calculator.

Usage:
  python fetch_ticker.py VZ --out vz_prices.csv --period 1y
"""
import argparse
import csv
import math
import sys

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ticker")
    ap.add_argument("--out", default=None)
    ap.add_argument("--period", default="1y")
    args = ap.parse_args()
    try:
        import yfinance as yf
    except Exception as e:
        print(f"yfinance not available: {e}", file=sys.stderr)
        sys.exit(1)
    t = yf.Ticker(args.ticker)
    hist = t.history(period=args.period, auto_adjust=True)
    if hist is None or hist.empty:
        print(f"No data for {args.ticker}", file=sys.stderr)
        sys.exit(2)
    hist = hist.dropna(subset=["Close"])
    closes = hist["Close"].tolist()
    dates = [d.strftime("%Y-%m-%d") for d in hist.index]
    out = args.out or f"{args.ticker.lower()}_prices.csv"
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["date", "close"])
        for d, c in zip(dates, closes):
            w.writerow([d, f"{float(c):.6f}"])
    # stats
    import numpy as np
    prices = np.array(closes, dtype=float)
    rets = np.diff(np.log(prices))
    vol = float(np.std(rets, ddof=1) * math.sqrt(252)) if len(rets) > 1 else float("nan")
    spot = float(prices[-1])
    print(f"ticker={args.ticker} spot={spot:.4f} ann_vol={vol:.4f} rows={len(prices)} out={out}")

if __name__ == "__main__":
    main()
