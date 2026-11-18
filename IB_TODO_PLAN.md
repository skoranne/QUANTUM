# IBKR (InteractiveBrokers.jl) Integration — Status & TODO Plan

## Goal
Get fast(ish) market data into the MINSPM Options TUI and evolve toward trading on IBKR using Julia (`InteractiveBrokers.jl`) for:
- Universe: **VZ, AMZN, NVDA, AMD, INTC, QCOM, AVGO, AAPL**
- ATM arithmetic-average Asian calls, **+2w / +4w**
- MINSPM coherent moment oracle + MINSPM QCS QPE pricing, Greeks (Δ/Γ/ν), convexity hedge, arbitrage/opportunity signals

## Current status (implemented)

### Dependency
- **`InteractiveBrokers.jl`** added to the Julia project (`Project.toml` / `Manifest.toml`).
- Package/module: `using InteractiveBrokers` (loaded conditionally via `HAS_IB` so the TUI still runs without it).

### TUI integration (`options_tui.jl`)
- **IBKR feed (streaming underlying spots):**
  - Connects to TWS/IB Gateway and starts the async reader (`InteractiveBrokers.start_reader`).
  - Subscribes `reqMktData` (SMART/STK/USD, streaming, non-snapshot) for all 8 tickers.
  - `tickPrice` handler captures fields `LAST`, `DELAYED_LAST`, `BID`, `DELAYED_BID`, `ASK`, `DELAYED_ASK`, `CLOSE`, `DELAYED_CLOSE` into a thread-safe feed (`IBKRFeed` with lock).
  - Optional market-data type request: realtime vs delayed (`--ib-delayed`).
- **Refresh behavior:** each render cycle fetches yfinance history (for vol/skew calibration) and **overrides spot** with latest IBKR prices when connected; MINSPM pricing/Greeks/hedge/signals recompute on the refreshed spot.
- **UI:** prominent `LAST UPDATE` timestamp and a `Data:` status line showing IBKR connection state, mode (realtime/delayed), and `live_prices=n/8` (or fallback/error).
- **Fallback:** if IBKR is off/unavailable/not connected, pricing continues via yfinance-only path. IBKR connection is disconnected on exit (interactive mode).

### CLI flags
```bash
# IBKR live (TWS/IB Gateway running with API enabled)
julia --project=. options_tui.jl --ibkr --ib-port 7497 --ib-client 7 --refresh 30

# one-shot
julia --project=. options_tui.jl --ibkr --ib-port 7497 --ib-client 7 --once

# delayed market data (no live subscription)
julia --project=. options_tui.jl --ibkr --ib-port 7497 --ib-client 7 --ib-delayed --refresh 30

# custom host
julia --project=. options_tui.jl --ibkr --ib-host 127.0.0.1 --ib-port 4002 --ib-client 7

# original yfinance-only mode
julia --project=. options_tui.jl --once
julia --project=. options_tui.jl --refresh 30
```
- `--ibkr` enable IBKR feed
- `--ib-port <int>` (default `7497`) — common: 7497 TWS paper, 7496 TWS live, 4002 Gateway paper, 4001 Gateway live
- `--ib-client <int>` (default `7`) — must be unique per API connection
- `--ib-host <host>` (default `127.0.0.1`)
- `--ib-delayed` request delayed market data type

## Prerequisites (to run against IBKR)
1. **TWS or IB Gateway** running and logged in (start with **paper**).
2. **API enabled:** TWS *Global Configuration → API → Settings*
   - Enable **ActiveX and Socket Clients** (socket API).
   - Confirm **Socket port** matches `--ib-port`.
   - Allow connections from localhost/trusted IP as needed.
   - (Read-Only API can remain enabled for market data; revisit for trading/orders.)
3. **Market data permissions:** without live subscriptions ticks may be delayed or absent — use `--ib-delayed` or ensure US equity data subscriptions are active.
4. **Verification:** run with `--ibkr` and confirm `Data: IBKR: connected ... live_prices=8/8` and that Spot updates with `LAST UPDATE`.

## Known gaps / limitations
- **Underlying spots only:** no live **options-chain** quotes/greeks from IBKR yet; current “ARB/OPPTY” signals compare MINSPM fair vs Euro-BS/QPE references (screening), not executable market option quotes.
- **History/vol calibration:** still sourced from yfinance daily closes (stable, but not tick/IB-historical driven).
- **Connection resilience:** minimal reconnect/error handling; reader/callback errors are surfaced in status but not auto-retried with backoff.
- **No order/execution path:** `placeOrder` / order management / risk controls not integrated.
- Not validated here against a live TWS/Gateway session (syntax-checked only).

## TODO plan

### Phase 1 — Reliable fast underlying data
- [ ] **Validate end-to-end** against TWS paper/Gateway: confirm ticks for all 8 symbols, correct port/clientId, and realtime vs delayed behavior.
- [ ] Improve **price selection/quality:** prefer LAST, fall back to midpoint(BID,ASK), track bid/ask/spread and timestamp/staleness; ignore clearly bad/zero prints.
- [ ] **Reconnect & health:** detect disconnects/errors, auto-reconnect with backoff, resubscribe; show stale-data warning if no ticks within N seconds.
- [ ] (Optional) Pull **historical bars from IBKR** (`reqHistoricalData`) to calibrate vol/skew, reducing dependence on yfinance and aligning with live feed.
- [ ] Config via env/args/file for host/port/clientId/mode; keep secrets/credentials out of repo.

### Phase 2 — Real options market data (executable quotes)
- [ ] **Resolve option contracts** for ATM **+2w/+4w** expiries per underlying:
  - Discover expiries/strikes (`reqSecDefOptParams` / `reqContractDetails`) and select nearest +2w/+4w expirations and ATM strike.
- [ ] **Subscribe to option market data** (`reqMktData` for OPT contracts): bid/ask/last/mid, IV, and model greeks if available (`tickOptionComputation`), or compute greeks from quotes.
- [ ] **Replace screening signals** with **MINSPM fair vs executable bid/ask/mid** edge (accounting for spread, fees/slippage assumptions). Highlight true mispricing/arbitrage candidates only when edge exceeds thresholds.
- [ ] Display chain/quote context in TUI (expiry, strike, bid/ask/mid, IV, OI/volume if available) alongside MINSPM fair and Greeks.
- [ ] Handle options data permissions, pacing limits, and contract-specific quirks.

### Phase 3 — Paper trading & execution (IBKR)
- [ ] **Order construction & routing (paper first):** build OPT orders from signals (`placeOrder`), defaulting to **limit** orders at justified prices; support what-if/checks where available.
- [ ] **Risk & safety controls (mandatory):**
  - Paper-trading mode by default; explicit confirmation/flag for live.
  - Max position/notional/order size, max daily loss, spread/slippage guard, stale-quote/price-band checks, kill switch.
  - Idempotent order handling and dedupe; handle partial fills/cancels.
- [ ] **Order/position lifecycle:** track `orderStatus`, executions, positions/PnL (`reqPositions`, `reqPnL`/account updates) and surface in TUI/logs.
- [ ] **Audit logging:** timestamped log of quotes, MINSPM fair/Greeks, signal, decision, and order/ack/fill for review and backtesting.

### Phase 4 — Robustness & performance
- [ ] Separate async tasks cleanly: feed ingestor, pricer/signal engine, UI/render, execution — with bounded queues/channels and no blocking on the render path.
- [ ] Profiling/latency measurement (tick-to-signal, signal-to-order) and optimization where it matters (data path > compute; MINSPM oracle already fast).
- [ ] Tests/sanity checks: model-vs-market comparisons, replay from recorded ticks, and guardrails for NaN/missing data.
- [ ] (If needed) Evaluate lower-latency data/execution paths (e.g. broker-specific optimizations, dedicated feed handlers) — only if sub-second streaming proves insufficient.

## Notes
- Treat all signals as **research/screening** until validated against executable quotes and risk controls. This is not financial advice.
- “Very fast” for a live TUI/signal workflow = reliable sub-second streaming; HFT/µs latency requires direct feeds, colocation, and specialized infra beyond this scope.
- Keep API ports/client IDs and any credentials/config local (do not commit secrets).
