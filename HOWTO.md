# Asian Option Pricing Tool (Classical MC + MINSPM QCS QPE/QAE)

This repo includes a practical **arithmetic-average Asian option** pricing tool for US underlyings, with:
- **Classical Monte Carlo** (CPU and **CUDA GPU** via Julia) as the high-fidelity baseline
- **MINSPM Julia QCS QPE/QAE** (bounded-entropy / bounded-codebook algebraic phase estimation) scalable to large path spaces (e.g. 128 assets × 252 steps)
- **Circuit-based QAE (QPE)** module in Python (Qiskit Aer, **prefers GPU/cuQuantum-backed Aer** when available, else CPU Aer)
- **Realized-payoff backtesting** from historical price CSVs

European-exercise arithmetic Asian (path-dependent) is the primary instrument. American early-exercise is out of scope.

## Prerequisites

- **Julia** with the repo project environment (`Project.toml` / `Manifest.toml`, Julia 1.13.x recorded). NVIDIA GPU + CUDA for GPU acceleration (optional but supported; CPU fallback works).
- **Python 3** with `qiskit` and `qiskit-aer` for `qae_asian.py` (GPU Aer/cuQuantum used automatically if the environment supports it). `yfinance` for fetching ticker history (`fetch_ticker.py`).

From the repo root:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
# python deps (example): pip install qiskit qiskit-aer yfinance
```

## 1) Classical Monte Carlo pricer (`asian_option_mc.jl`)

Run from the repo root:

```bash
# Sanity/consistency checks (CPU/GPU, SE scaling, bounds, ITM/OTM, CPU vs GPU parity)
julia --project=. asian_option_mc.jl validate

# Price SOXX-weighted basket Asian call from repo CSVs (soxx_asset_prices.csv, soxx_correlation_matrix.csv)
julia --project=. asian_option_mc.jl price-soxx

# Realized-payoff rolling backtest on SOXX history -> asian_backtest_soxx_results.csv
julia --project=. asian_option_mc.jl backtest-soxx

# Price a single US stock arithmetic Asian call from a ticker price CSV (date,close)
julia --project=. asian_option_mc.jl price-ticker <prices.csv> [K] [T] [r] [steps] [paths]
```

Notes:
- Uses corrected per-step NIG-copula path dynamics, Cholesky correlation (PSD-repaired if needed), and reports price, standard error, amplitude, and engine (`gpu-mc-cuda` or `cpu-mc`).
- `price-ticker` defaults: ATM strike `K=S0`, `T=1.0`, `r=0.0415`, `steps=252`, `paths=50000`, and prints a European Black–Scholes reference (same vol) for sanity (arithmetic Asian is typically below the European).

### Fetch current US stock history (VZ, AMZN, etc.)

```bash
python fetch_ticker.py VZ   --period 1y   # writes vz_prices.csv, prints spot/ann vol
python fetch_ticker.py AMZN --period 1y   # writes amzn_prices.csv

julia --project=. asian_option_mc.jl price-ticker vz_prices.csv
julia --project=. asian_option_mc.jl price-ticker amzn_prices.csv
```

Example (1Y ATM arithmetic Asian call, GPU MC, calibrated from 1y history): VZ ~\$2.53 (BS ref ~\$5.84), AMZN ~\$16.67 (BS ref ~\$40.38) — Asian < European is expected due to averaging.

## 2) MINSPM Julia QCS QPE/QAE (`minspm_qcs_qae.jl`)

Bounded-entropy / bounded-codebook MINSPM QCS: 128-asset × 252-step (continuous path space) oracle via (GPU) MC sampling, **quantized amplitude codebook** (bounded H), and **QPE/QAE phase estimation** over a bounded Galois/cyclotomic phase register (`Z_{2^k}`, truncation `k`) with MINSPM direct (O(1)) and refined (u-class) algebraic recovery.

```bash
julia --project=. minspm_qcs_qae.jl
julia --project=. minspm_qcs_qae.jl --n-paths 20000 --k-trunc 8 --m-bits 28 --u-classes 8192
```

Reports: oracle/codebook amplitudes & prices, explicit bounded-register QAE, MINSPM refined and direct-Galois QPE/QAE prices, errors vs oracle, and the bounded entropy/codebook summary. Requires `soxx_asset_prices.csv` and `soxx_correlation_matrix.csv` in the repo root.

## 3) Circuit-based QAE/QPE module (`qae_asian.py`)

Real circuit-based amplitude estimation on a **tractable discretized** Asian payoff (exact full 128×252 state-prep is infeasible on simulators). Prefers **GPU Aer (cuQuantum/cuStateVec-backed when available)**, falls back to CPU Qiskit Aer.

```bash
python qae_asian.py --S0 100 --K 100 --r 0.03 --sigma 0.2 --T 1 --steps 1 --nq-price 3 --m-rounds 5 --shots 256
```

Outputs JSON: discretized reference, circuit-prepared **statevector** amplitude (exact for the circuit), **direct ancilla sampling** estimate (+SE), and canonical QPE-on-Grover estimate, plus backend/qubit info. Use statevector/sampling as the authoritative circuit results.

## Data files

- `soxx_asset_prices.csv`, `soxx_correlation_matrix.csv`, `GEN_CSV/` — historical/calibration inputs for basket pricing and backtests.
- Fetched ticker CSVs (e.g. `vz_prices.csv`, `amzn_prices.csv`) are generated locally.

## Caveats

- Prices are model-based (NIG-copula / discretized assumptions) and depend on calibration inputs and sampling error (see reported SE). Backtests are **realized-payoff** comparisons, not traded option quotes.
- GPU paths require a working CUDA setup (`CUDA.functional()` in Julia; GPU-capable Qiskit Aer for the Python module). All tools fall back to CPU where supported.
- This is a research/pricing assistant tool, not financial advice.
