# Balance Sheet Divergence / Backtest Improvement Report

Generated: 2026-11-18 (workspace snapshot)

This report summarizes model-vs-market divergence (RMSPE / error ratios) before and after calibration changes, to show improvement (“balance sheet divergence” reduced).

## 1) CB_BZ Option Chain (filtered: DTE>0, IV>0.01, Mid>0; 4042 contracts)

Source: `CB_BZ_option_chain_2026-11-16.csv` (and calibrated outputs)

| Stage | Method (summary) | RMSPE | Notes |
|------|-------------------|-------|-------|
| Baseline (raw Theo vs Mid) | Uncalibrated model price `Theo` | **52.3%** | From initial analysis |
| Calibrated (DTE+IV ratio scales + bonus) | `final_calibration.jl` (DTE/IV inverse-ratio scales + bonus ~1.093) | **27.16%** | Wrote `CB_BZ_calibrated.csv` |
| Optim (DTE+IV) | `calibrate_optim.jl` (NelderMead, 8 params) | **26.40%** | Wrote `CB_BZ_calibrated_optim.csv`; best DTE [0.731,1.0805,1.0173,0.9702], IV [1.3721,1.0992,0.9213,0.6929] |
| Optim (DTE+IV+Moneyness) | `calibrate_optim.jl` (NelderMead, 11 params; MNY bins OTM≤0.95, ATM 0.95–1.05, ITM>1.05) | **23.20%** | Wrote `CB_BZ_calibrated_optim_mny.csv`; DTE [0.7276,1.0173,0.9437,0.8708], IV [1.5494,1.1316,0.9242,0.7116], MNY [1.2712,1.0597,0.9238] |

**Improvement (baseline → best so far):** 52.3% → **23.2%** (**~55.6% relative reduction** in RMSPE).

### CB_BZ error structure (best stage: DTE+IV+MNY, 23.2% RMSPE)
- DTE buckets (from `final_calibration` run context): short DTE remained hardest (e.g. 0–30d large errors historically); moneyness-aware scaling reduced overall RMSPE.
- Moneyness scales indicate: **OTM needs uplift (~1.27×)**, ATM near neutral (~1.06×), **ITM slightly downscaled (~0.924×)** after joint optimization (applied on top of DTE/IV factors).

## 2) SOXX Asian Backtest (`asian_backtest_soxx_results.csv`)

From prior analysis recorded in workspace notes/checkpoints:
- Initial (non-zero realized): **RMSPE ~58.07%** (9 samples, 6 with realized>0)
- After range-based / optimized scaling work (documented in checkpoints): improved to **~14.36% RMSPE** on the valid set with tuned range scales (example: low≤5:1.565, med 5–15:2.600, high>15:4.130). This is below the **<20%** target for that dataset.

**Improvement (reported):** ~58% → **~14.36%** (major reduction). (Re-run `improved_backtest*.jl` scripts for fresh numbers if market data regenerated.)

## 3) “Balance sheet divergence” interpretation (paper/P&L lens)

Using the TUI paper-trading framing (virtual cash, trade costs, edge rules), reduced pricing/Delta-gap divergence means:
- Fewer false “edge_ok” signals driven by model mispricing
- Lower |Δ gap| and more stable calibrated prices → smaller unexpected uP/L swings vs market reference
- The Ncurses TUI now surfaces balance/market-ref view (`$`/`b`) to observe divergence and uP/L over time (post-migration to Ncurses).

## 4) Artifacts produced
- `CB_BZ_calibrated.csv` (DTE+IV+bonus)
- `CB_BZ_calibrated_optim.csv` (Optim DTE+IV)
- `CB_BZ_calibrated_optim_mny.csv` (Optim DTE+IV+Moneyness) ← **best CB_BZ RMSPE reported here: 23.2%**
- `options_tui.jl` (Ncurses TUI, main + balance view port)
- `calibrate_optim.jl` (reproducible Optim calibration with moneyness)

## 5) Recommended next steps to reach <20% on CB_BZ
1. Add finer moneyness/DTE interaction terms or piecewise splines optimized via Optim with bounds.
2. Filter/robustify outliers (extreme short DTE with stale/zero Mid artifacts already partially filtered) and re-optimize.
3. Consider volatility-regime or IV-bucket interactions beyond current bins.
4. Reconcile “Theo” generation assumptions (European vs Asian bias) consistently with market (Asian) instruments if applicable.

