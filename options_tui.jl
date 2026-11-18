#!/usr/bin/env julia
# MINSPM Option Pricing TUI — ATM Arithmetic Asian, +2w/+4w, VZ AMZN NVDA AMD INTC QCOM AVGO AAPL
# Engine: self-sufficient coherent MINSPM moment oracle + MINSPM QCS QPE (deterministic, no MC sampling)
# Run: julia --project=. options_tui.jl [--refresh 30] [--once]

using Printf, Dates, Statistics, LinearAlgebra, Distributions, SpecialFunctions
include(joinpath(@__DIR__, "coherent_oracle_minspm.jl"))
const CM = CoherentOracleMINSPM

const HAS_UNICODEPLOTS = try
    @eval using UnicodePlots
    true
catch
    false
end

const TICKERS = begin
    f = joinpath(@__DIR__, "tickers.txt")
    if isfile(f)
        [String(strip(l)) for l in readlines(f) if !isempty(strip(l)) && !startswith(strip(l),"#")]
    else
        ["VZ","AMZN","NVDA","AMD","INTC","QCOM","AVGO","AAPL","SPY","QQQ","PLTR"]
    end
end
const TENOR_WEEKS = [1, 2, 4]
const RISK_FREE = 0.0415
const TRADING_DAYS = 252
const BANNER = "designed and implemented by Sandeep Koranne (C) 2026- All rights reserved."
const FIXED_TRADE_COST = 0.50
const TRADE_EDGE_MIN = 1.0
const WARMUP_ITERS = 5
const GAMMA_MAX = 0.05
const DDELTA_MAX = 0.15
const VIRTUAL_CASH_START = 1000.0
const virtual_cash = Ref(VIRTUAL_CASH_START)
const paper_positions = Ref(Dict{Tuple{String,Int},NamedTuple}())
const warmup_count = Ref(0)
const warmup_done = Ref(false)

function py_cmd()
    for p in ("/home/skoranne/miniconda3/bin/python3","/home/skoranne/miniconda3/bin/python","python3","python")
        try
            run(pipeline(`$p -c "import yfinance"`, stdout=devnull, stderr=devnull))
            return p
        catch; end
    end
    return "python3"
end
const PYTHON = Ref{String}("python3")

function fetch_snapshot_lines(tickers::AbstractVector)
    py = PYTHON[]
    syms = join(["\"" * String(t) * "\"" for t in tickers], ",")
    script = "import json,math\n" *
             "import numpy as np\n" *
             "import yfinance as yf\n" *
             "syms=[" * syms * "]\n" *
             "for sym in syms:\n" *
             "  try:\n" *
             "    tk=yf.Ticker(sym)\n" *
             "    hist=tk.history(period='1y',auto_adjust=True)\n" *
             "    closes=[]\n" *
             "    if hist is not None and (not hist.empty) and 'Close' in hist:\n" *
             "      closes=[float(x) for x in hist['Close'].dropna().tolist()]\n" *
             "    spot=None\n" *
             "    try:\n" *
             "      fi=tk.fast_info\n" *
             "      spot=fi.get('lastPrice') or fi.get('last_price') or fi.get('regularMarketPrice')\n" *
             "    except Exception:\n" *
             "      spot=None\n" *
             "    if spot is None and len(closes)>0: spot=closes[-1]\n" *
             "    if spot is None:\n" *
             "      print(sym+'|NONE|')\n" *
             "    else:\n" *
             "      tail=','.join([f'{c:.6f}' for c in closes[-260:]])\n" *
             "      print(sym+'|'+f'{float(spot):.6f}'+'|'+tail)\n" *
             "  except Exception as e:\n" *
             "    print(sym+'|ERR|'+str(e).replace('|','/'))\n"
    result = Dict{String,NamedTuple}()
    try
        out = read(pipeline(`$py -c $script`, stderr=devnull), String)
        for line in split(out, '\n')
            sline = strip(line); isempty(sline) && continue
            parts = split(sline, '|')
            length(parts) < 2 && continue
            sym = String(parts[1])
            if parts[2] == "NONE" || parts[2] == "ERR"
                result[sym] = (spot=NaN, closes=Float64[], ok=false)
                continue
            end
            spot = parse(Float64, parts[2])
            closes = Float64[]
            if length(parts) >= 3 && !isempty(parts[3])
                for x in split(parts[3], ',')
                    sx = strip(x); isempty(sx) && continue
                    push!(closes, parse(Float64, sx))
                end
            end
            result[sym] = (spot=spot, closes=closes, ok=true)
        end
    catch e
        @warn "fetch failed" exception=e
    end
    return result
end


function fetch_option_market_lines(tickers::AbstractVector)
    py = PYTHON[]
    syms = join(["\"" * String(t) * "\"" for t in tickers], ",")
    script = "import math\n" *
             "from datetime import datetime, timezone\n" *
             "try:\n" *
             "    import yfinance as yf\n" *
             "except Exception as e:\n" *
             "    raise SystemExit\n" *
             "syms=[" * syms * "]\n" *
             "r = float(" * string(RISK_FREE) * ")\n" *
             "for sym in syms:\n" *
             "    try:\n" *
             "        tk = yf.Ticker(sym)\n" *
             "        spot = None\n" *
             "        try:\n" *
             "            fi = tk.fast_info\n" *
             "            spot = fi.get('lastPrice') or fi.get('last_price') or fi.get('regularMarketPrice')\n" *
             "        except Exception:\n" *
             "            spot = None\n" *
             "        if spot is None:\n" *
             "            try:\n" *
             "                h = tk.history(period='5d', auto_adjust=True)\n" *
             "                if h is not None and len(h)>0 and 'Close' in h:\n" *
             "                    spot = float(h['Close'].dropna().iloc[-1])\n" *
             "            except Exception:\n" *
             "                spot = None\n" *
             "        if spot is None or not (float(spot) > 0):\n" *
             "            print(sym + '|NONE')\n" *
             "            continue\n" *
             "        spot = float(spot)\n" *
             "        today = datetime.now(timezone.utc).date()\n" *
             "        opts = list(tk.options or [])\n" *
             "        targets = [(7,1),(14,2),(28,4)]\n" *
             "        for target_days, weeks in targets:\n" *
             "            best = None\n" *
             "            best_diff = 1e9\n" *
             "            for e in opts:\n" *
             "                try:\n" *
             "                    ed = datetime.strptime(e, '%Y-%m-%d').date()\n" *
             "                    dd = (ed - today).days\n" *
             "                    if dd < 2:\n" *
             "                        continue\n" *
             "                    if abs(dd - target_days) < best_diff:\n" *
             "                        best_diff = abs(dd - target_days)\n" *
             "                        best = (e, dd)\n" *
             "                except Exception:\n" *
             "                    continue\n" *
             "            if best is None:\n" *
             "                print(sym + '|' + str(weeks) + '|NOEXP')\n" *
             "                continue\n" *
             "            exp_str, ddays = best\n" *
             "            T = max(ddays / 365.0, 1.0/365.0)\n" *
             "            try:\n" *
             "                oc = tk.option_chain(exp_str)\n" *
             "                calls = oc.calls\n" *
             "            except Exception:\n" *
             "                print(sym + '|' + str(weeks) + '|' + exp_str + '|ERR')\n" *
             "                continue\n" *
             "            if calls is None or len(calls) == 0:\n" *
             "                print(sym + '|' + str(weeks) + '|' + exp_str + '|NONE')\n" *
             "                continue\n" *
             "            try:\n" *
             "                calls = calls.copy()\n" *
             "                calls['dsk'] = (calls['strike'].astype(float) - spot).abs()\n" *
             "                row = calls.sort_values('dsk').iloc[0]\n" *
             "            except Exception:\n" *
             "                print(sym + '|' + str(weeks) + '|' + exp_str + '|NONE')\n" *
             "                continue\n" *
             "            K = float(row['strike'])\n" *
             "            bid = float(row.get('bid', 0) or 0)\n" *
             "            ask = float(row.get('ask', 0) or 0)\n" *
             "            lastp = float(row.get('lastPrice', 0) or 0)\n" *
             "            if bid > 0 and ask > 0 and ask >= bid:\n" *
             "                optp = 0.5*(bid+ask)\n" *
             "            elif lastp > 0:\n" *
             "                optp = lastp\n" *
             "            elif bid > 0:\n" *
             "                optp = bid\n" *
             "            elif ask > 0:\n" *
             "                optp = ask\n" *
             "            else:\n" *
             "                optp = float('nan')\n" *
             "            iv = float(row.get('impliedVolatility', 0) or 0)\n" *
             "            if not (iv > 0 and iv < 5.0):\n" *
             "                iv = 0.30\n" *
             "            from math import log, sqrt, exp, erf\n" *
             "            def norm_cdf(x):\n" *
             "                return 0.5*(1.0+erf(x/1.4142135623730951))\n" *
             "            def norm_pdf(x):\n" *
             "                return math.exp(-0.5*x*x)/2.5066282746310002\n" *
             "            if T <= 0:\n" *
             "                delta = 1.0 if spot>K else 0.0; gamma = 0.0; vega = 0.0\n" *
             "            else:\n" *
             "                try:\n" *
             "                    d1 = (log(max(spot,1e-12)/K) + (r + 0.5*iv*iv)*T)/(iv*sqrt(T))\n" *
             "                    delta = norm_cdf(d1)\n" *
             "                    gamma = norm_pdf(d1)/(spot*iv*sqrt(T)) if spot>0 and iv>0 else 0.0\n" *
             "                    vega = spot*norm_pdf(d1)*sqrt(T)/100.0\n" *
             "                except Exception:\n" *
             "                    delta = float('nan'); gamma = float('nan'); vega = float('nan')\n" *
             "            print(sym + '|' + str(weeks) + '|' + ('%.6f'%optp) + '|' + ('%.6f'%delta) + '|' + ('%.6f'%gamma) + '|' + ('%.6f'%vega) + '|OK')\n" *
             "    except Exception as e:\n" *
             "        print(sym + '|FATAL')\n"
    result = Dict{String,Dict{Int,NamedTuple}}()
    try
        out = read(pipeline(`$py -c $script`, stderr=devnull), String)
        for line in split(out, '\n')
            sline = strip(line); isempty(sline) && continue
            parts = split(sline, '|')
            if length(parts) >= 7 && parts[end] == "OK"
                try
                    sym = parts[1]; weeks = parse(Int, parts[2])
                    optp = parse(Float64, parts[3]); md = parse(Float64, parts[4]); mg = parse(Float64, parts[5]); mv = parse(Float64, parts[6])
                    d = get!(result, sym, Dict{Int,NamedTuple}())
                    d[weeks] = (price=optp, delta=md, gamma=mg, vega=mv)
                catch; continue; end
            end
        end
    catch e
        @warn "fetch_option_market_lines failed" exception=e
    end
    return result
end

struct PriceResult
    sym::String; weeks::Int; T::Float64; steps::Int
    spot::Float64; strike::Float64; vol::Float64; skew::Float64
    price::Float64; a::Float64
    delta::Float64; gamma::Float64; vega::Float64
    bs_euro::Float64; qpe_price::Float64
end

function bs_call(S,K,r,σ,T)
    T <= 0 && return max(S-K,0.0)
    σ = max(σ,1e-8)
    d1 = (log(S/K) + (r + 0.5*σ^2)*T) / (σ*sqrt(T))
    d2 = d1 - σ*sqrt(T)
    N(x)=0.5*(1.0+erf(x/sqrt(2.0)))
    return S*N(d1) - K*exp(-r*T)*N(d2)
end

function mc_asian_price(S,K,r,σ,T,steps; npaths=10000)
    T <= 0 && return max(S-K,0.0)
    dt = T / max(steps,1)
    drift = (r - 0.5*σ^2) * dt
    vol = σ * sqrt(dt)
    pay = 0.0
    for _ in 1:npaths
        St = S
        ssum = 0.0
        for _ in 1:steps
            St *= exp(drift + vol * randn())
            ssum += St
        end
        avg = ssum / steps
        pay += max(avg - K, 0.0)
    end
    return exp(-r*T) * pay / npaths
end

function price_one(sym::String, spot::Float64, closes::Vector{Float64}, weeks::Int; r=RISK_FREE)
    T = weeks / 52.0
    steps = max(10, Int(round(T * TRADING_DAYS)))
    vol, sk = length(closes) >= 30 ? CM.calibrate_1d(closes) : (0.30, -0.10)
    vol = clamp(vol, 0.05, 1.50); sk = clamp(sk, -0.8, -0.01)
    K = spot
    Bnd = max(0.8*spot, 2*vol*spot*sqrt(T))
    base = CM.coherent_oracle_moments(spot,K,r,vol,sk,T,steps; bound_b=Bnd)
    theta = asin(sqrt(clamp(base.a,1e-12,1.0-1e-12)))
    _, a_ref = CM.minspm_refined(theta, 1<<24, 4096)
    qpe_price = a_ref * Bnd * exp(-r*T)
    dS = max(spot*0.005, 0.01)
    up = CM.coherent_oracle_moments(spot+dS,K,r,vol,sk,T,steps; bound_b=Bnd)
    dn = CM.coherent_oracle_moments(spot-dS,K,r,vol,sk,T,steps; bound_b=Bnd)
    delta = (up.price_unbounded - dn.price_unbounded) / (2*dS)
    gamma = (up.price_unbounded - 2*base.price_unbounded + dn.price_unbounded) / (dS^2)
    dσ = 0.01
    vup = CM.coherent_oracle_moments(spot,K,r,vol+dσ,sk,T,steps; bound_b=Bnd)
    vdn = CM.coherent_oracle_moments(spot,K,r,max(vol-dσ,0.01),sk,T,steps; bound_b=Bnd)
    vega = (vup.price_unbounded - vdn.price_unbounded) / (2*dσ) / 100.0
    bse = bs_call(spot,K,r,vol,T)
    return PriceResult(sym,weeks,T,steps,spot,K,vol,sk,base.price_unbounded,base.a,
                       delta,gamma,vega,bse,qpe_price)
end

struct Signal; level::Symbol; text::String; end

function analyze(pr::PriceResult)
    sigs = Signal[]
    discount = pr.bs_euro > 1e-6 ? (pr.bs_euro - pr.price)/pr.bs_euro : 0.0
    if pr.price > pr.bs_euro * 1.05
        push!(sigs, Signal(:warn, "Asian>Euro(check)"))
    end
    if pr.price < pr.bs_euro * 0.35
        push!(sigs, Signal(:arb, @sprintf("Cheap vs Euro-ref d=%.0f%%", discount*100)))
    elseif discount < 0.10
        push!(sigs, Signal(:oppty, @sprintf("Tight disc %.0f%%", discount*100)))
    else
        push!(sigs, Signal(:ok, @sprintf("Disc %.0f%%", discount*100)))
    end
    if abs(pr.gamma * pr.spot) > 0.05
        push!(sigs, Signal(:warn, @sprintf("High gamma %.4f", pr.gamma)))
    end
    return sigs
end

const CR="\e[0m"; const CB="\e[1m"; const CD="\e[2m"
const CRED="\e[31m"; const CGRN="\e[32m"; const CYEL="\e[33m"; const CCYN="\e[36m"
const CBGRN="\e[42m"

function sigstr(s::Signal)
    s.level==:arb   && return CBGRN*CB*" ARB "*CR*" "*CGRN*s.text*CR
    s.level==:oppty && return CYEL*"OPPTY"*CR*" "*s.text
    s.level==:warn  && return CRED*"WARN"*CR*" "*s.text
    return CD*s.text*CR
end

function render(results::Vector{PriceResult}, last_update::String, refresh::Int)
    print("\e[2J\e[H")
    println(CB,CCYN,"╔══════════════════════════════════════════════════════════════════════════════╗",CR)
    println(CB,CCYN,"║  MINSPM TUI (Asian ATM, QPE + Coherent Oracle)                               ║",CR)
    println(CB,CCYN,"╚══════════════════════════════════════════════════════════════════════════════╝",CR)
    println(CD,"  ",last_update,"  refresh:",refresh,"s  r=",RISK_FREE," +1w/+2w/+4w   (keys: q then Enter)",CR)
    println(repeat("-", 90))
    @printf("%s%-6s %2s %9s %6s %7s %7s %7s %7s %8s%s\n", CB,"Sym","Tw","Spot","Vol%","Asian","Delta","Gamma","Vega","EuroBS",CR)
    println(repeat("-", 90))
    bysym = Dict{String,Vector{PriceResult}}()
    for pr in results; push!(get!(bysym, pr.sym, PriceResult[]), pr); end
    for sym in TICKERS
        prs = sort(get(bysym, sym, PriceResult[]), by=x->x.weeks); isempty(prs) && continue
        for (i,pr) in enumerate(prs)
            sigs = analyze(pr)
            hasarb = any(s.level==:arb for s in sigs); hasop = any(s.level==:oppty for s in sigs)
            col = hasarb ? CGRN : (hasop ? CYEL : "")
            symc = i==1 ? sym : "      "
            @printf("%s%-6s %+1dw %9.2f %5.1f %7.3f %7.3f %7.4f %7.4f %8.3f%s\n",
                    col,symc,pr.weeks,pr.spot,pr.vol*100,pr.price,pr.delta,pr.gamma,pr.vega,pr.bs_euro,CR)
        end
    end
    println(repeat("-", 90))
    println(CB,"EXPLANATIONS (right panel):",CR)
    println("  Disc%/Signals: shown in analyze() as discount vs EuroBS (Disc = (EuroBS-Asian)/EuroBS). ")
    println("    e.g. 'Disc 46%' means Asian model price is 46% below EuroBS reference (cheaper).")
    println("  Greeks (model, MINSPM bump-reprice): Delta dP/dS, Gamma d²P/dS², Vega dP/dσ (per 1% vol).")
    println("  Market greeks: in Balance view (\$ + Enter) shown as dΔ/dΓ/dV vs Yahoo BS-IV (RED if |dΔ|>=0.05).")
    println("  Trades (paper, Balance view): start \$", VIRTUAL_CASH_START, ", \$", FIXED_TRADE_COST, " cost per new (sym,tw) when mkt price first seen; uP/L = mkt-now - entry.")
    println()
    println(CD, BANNER, CR)
end

function compute_all(snapshot)
    results = PriceResult[]
    for sym in TICKERS
        data = get(snapshot, sym, nothing); data === nothing && continue
        (!data.ok || isnan(data.spot)) && continue
        for w in TENOR_WEEKS
            try; push!(results, price_one(sym, data.spot, data.closes, w))
            catch e; @warn "pricing failed" sym w exception=e; end
        end
    end
    return results
end

function main()
    refresh = 30; once = false
    i = 1
    while i <= length(ARGS)
        if ARGS[i] == "--refresh" && i < length(ARGS); refresh = parse(Int,ARGS[i+1]); i+=2
        elseif ARGS[i] == "--once"; once=true; i+=1
        else; i+=1; end
    end
    PYTHON[] = py_cmd()
    println("MINSPM TUI starting (python=$(PYTHON[]), refresh=$(refresh)s). Keys: q/\$/m + Enter (cooked).")
    if once
        snapshot = fetch_snapshot_lines(TICKERS)
        results = compute_all(snapshot)
        render(results, Dates.format(now(),"yyyy-mm-dd HH:MM:SS"), refresh)
        println("\n(--once) done.")
        return
    end
    quit = Ref(false)
    view_mode = Ref(:main)
    results_ref = Ref(PriceResult[])
    market_opts_ref = Ref(Dict{String,Dict{Int,Any}}())
    @async begin
        try
            while !quit[]
                try
                    line = readline(stdin)
                    cmd = lowercase(strip(line))
                    if cmd == "q"
                        quit[] = true; break
                    elseif cmd == "\$" || cmd == "b" || cmd == "balance"
                        view_mode[] = :balance
                    elseif cmd == "m" || cmd == "main"
                        view_mode[] = :main
                    end
                catch
                    quit[] = true; break
                end
            end
        catch; quit[] = true; end
    end
    while !quit[]
        snapshot = fetch_snapshot_lines(TICKERS)
        results_ref[] = compute_all(snapshot)
        try
            market_opts_ref[] = fetch_option_market_lines(TICKERS)
        catch; market_opts_ref[] = Dict{String,Dict{Int,Any}}(); end
        if !warmup_done[]
            warmup_count[] += 1
            if warmup_count[] >= WARMUP_ITERS
                warmup_done[] = true
            end
        end
        if view_mode[] == :balance
            print("\e[2J\e[H")
            println(CB,CCYN,"╔══════════════════════════════════════════════════════════════════════════════╗",CR)
            println(CB,CCYN,"║  BALANCE / GREEKS vs Yahoo BS-IV (market-ref)                                ║",CR)
            println(CB,CCYN,"╚══════════════════════════════════════════════════════════════════════════════╝",CR)
            println(CD,"  ",Dates.format(now(),"HH:MM:SS"),CR)
            if !warmup_done[]
                println(CYEL, @sprintf("WARMUP %d/%d: no trades opened (calibrating from recent data)...", warmup_count[], WARMUP_ITERS), CR)
            else
                println(@sprintf("Virtual \$%.2f start | cost \$%.2f/trade | open if edge>=\$%.2f & |γ|<=%.2f & |dΔ|<%.2f", VIRTUAL_CASH_START, FIXED_TRADE_COST, FIXED_TRADE_COST+TRADE_EDGE_MIN, GAMMA_MAX, DDELTA_MAX))
            println(CD, "live calibration each refresh; drift filter: block if |1d ret|>3σ_daily", CR)
            end
            println(repeat("-", 76))
            anym = false; dgaps = Float64[]; dlabels = String[]; unreal = 0.0
            for pr in results_ref[]
                mopt = get(get(market_opts_ref[], pr.sym, Dict{Int,Any}()), pr.weeks, nothing)
                if mopt !== nothing && isfinite(get(mopt,:price,NaN)) && mopt.price > 0
                    anym = true
                    key = (pr.sym, pr.weeks)
                    pos = paper_positions[]
                    edge = Float64(mopt.price) - pr.qpe_price
                    dgap_mkt = pr.delta - get(mopt,:delta,NaN)
                    # drift filter: block if |1d return| > 3 * sigma_daily
                    drift_block = false
                    try
                        closes_sym = get(closes_map, pr.sym, Float64[])
                        if length(closes_sym) >= 2
                            r1d = abs((closes_sym[end] - closes_sym[end-1]) / max(closes_sym[end-1], 1e-12))
                            sigd = pr.vol / sqrt(252.0)
                            if sigd > 1e-8 && r1d > 3.0 * sigd
                                drift_block = true
                            end
                        end
                    catch; end
                    can_trade = warmup_done[] && !drift_block && edge >= (FIXED_TRADE_COST + TRADE_EDGE_MIN) && abs(pr.gamma) <= GAMMA_MAX
                    if isfinite(dgap_mkt) && abs(dgap_mkt) >= DDELTA_MAX
                        can_trade = false
                    end
                    if !haskey(pos, key) && can_trade
                        virtual_cash[] = virtual_cash[] - FIXED_TRADE_COST
                        pos2 = copy(pos)
                        pos2[key] = (entry_mkt = Float64(mopt.price), entry_t = now())
                        paper_positions[] = pos2
                        pos = paper_positions[]
                    end
                    if haskey(pos, key)
                        entry_mkt = pos[key].entry_mkt
                        upl = (Float64(mopt.price) - entry_mkt)
                        unreal += upl
                        dgap = pr.delta - get(mopt,:delta,NaN)
                        col = abs(dgap) < 0.05 ? CGRN : CRED
                        rc = ""
                        if upl < -FIXED_TRADE_COST
                            if abs(pr.gamma) > 0.02
                                rc = " | ROOT: convexity/gamma risk (model curvature/bump sensitivity); fix: tighten gamma calc, rebalance size"
                            elseif abs(dgap) >= 0.10
                                rc = " | ROOT: Δ-gap/calib mismatch (stale vol/skew vs market); fix: re-calibrate more often / warmup longer"
                            else
                                rc = " | ROOT: adverse drift: QPE+MINSPM assumed path/vol move not realized (drift filter active: |1d|>3σ blocks opens); fix: reduce edge or tighten recalib"
                            end
                        end
                        println(col, @sprintf(" %-6s+%dw dΔ=%+.3f mkt=%.2f ent=%.2f uP/L=%+.2f%s", pr.sym, pr.weeks, dgap, Float64(mopt.price), entry_mkt, upl, rc), CR)
                    else
                        dgap = pr.delta - get(mopt,:delta,NaN)
                        col = abs(dgap) < 0.05 ? CGRN : CRED
                        if !warmup_done[]
                            status = "(warmup)"
                        elseif drift_block
                            status = "(blocked: |1d ret|>3σ_daily)"
                        elseif abs(pr.gamma) > GAMMA_MAX
                            status = "(blocked: |γ|>=$(GAMMA_MAX))"
                        elseif isfinite(dgap) && abs(dgap) >= DDELTA_MAX
                            status = "(blocked: |dΔ|>=$(DDELTA_MAX))"
                        elseif edge < (FIXED_TRADE_COST + TRADE_EDGE_MIN)
                            status = "(no trade)"
                        else
                            status = "(edge ok)"
                        end
                        println(col, @sprintf(" %-6s+%dw dΔ=%+.3f mkt=%.2f QPE*=%.2f edge=%.2f %s", pr.sym, pr.weeks, dgap, Float64(mopt.price), pr.qpe_price, edge, status), CR)
                    end
                    push!(dgaps, dgap); push!(dlabels, string(pr.sym,"+",pr.weeks))
                else
                    println(CD, @sprintf(" %-6s+%dw (no mkt)", pr.sym, pr.weeks), CR)
                end
            end
            equity = virtual_cash[] + unreal
            println(repeat("-", 76))
            println(@sprintf("Cash: \$%.2f | Unreal P/L: \$%.2f | Equity: \$%.2f", virtual_cash[], unreal, equity))
            # simple heuristic portfolio rebalance suggestion (inverse |gamma|)
            open_syms = String[]
            gammas = Float64[]
            for pr in results_ref[]
                key = (pr.sym, pr.weeks)
                if haskey(paper_positions[], key)
                    push!(open_syms, string(pr.sym,"+",pr.weeks)); push!(gammas, abs(pr.gamma) + 1e-6)
                end
            end
            if !isempty(open_syms)
                w = 1.0 ./ gammas; w = w ./ sum(w)
                sugg = ["$(open_syms[i]):\$" * string(round(equity*w[i],digits=0)) for i in eachindex(open_syms)]
                println(CD, "Rebalance (heuristic, inv|gamma|, spread \$ across open): ", join(sugg, " "), "  (JuMP optimizer TBD)", CR)
            end
            if !warmup_done[]
                println(CYEL, "ROOT-CAUSE: trades blocked during warmup (insufficient calibration).", CR)
            else
                println(CD, "ROOT-CAUSE on losers (uP/L<-cost): convexity/gamma, Δ-gap/calib drift, or adverse QPE move.", CR)
            end
            if anym && HAS_UNICODEPLOTS
                try
                    labs2 = [string(dlabels[i], dgaps[i]>=0 ? "+" : "-") for i in eachindex(dgaps)]
                    println("\n", UnicodePlots.barplot(labs2, abs.(dgaps), title="|Δ gap| model-mktref (BS-IV)"))
                catch e; end
            end
            # cumulative $ per position (unrealized uP/L) plot
            upl_labels = String[]; upl_vals = Float64[]
            for pr in results_ref[]
                key = (pr.sym, pr.weeks)
                mopt = get(get(market_opts_ref[], pr.sym, Dict{Int,Any}()), pr.weeks, nothing)
                if mopt !== nothing && isfinite(get(mopt,:price,NaN)) && mopt.price > 0 && haskey(paper_positions[], key)
                    entry_mkt = paper_positions[][key].entry_mkt
                    upl = Float64(mopt.price) - entry_mkt
                    push!(upl_labels, string(pr.sym,"+",pr.weeks)); push!(upl_vals, upl)
                end
            end
            if !isempty(upl_vals) && HAS_UNICODEPLOTS
                try
                    println("\n", UnicodePlots.barplot(upl_labels, upl_vals, title="Cumulative \$ uP/L per position (mkt-entry)"))
                catch e; end
            elseif !isempty(upl_vals)
                println("Cumulative \$ uP/L per position: " * join(["$(upl_labels[i])=$(round(upl_vals[i],digits=2))" for i in eachindex(upl_vals)], ", "))
            end
            # MC sanity check on first few
            println(repeat("-", 76))
            println(CB, "MC sanity check (GBM, calib σ) vs QPE/model (first 3):", CR)
            cnt = 0
            for pr in results_ref[]
                cnt >= 3 && break
                try
                    mc = mc_asian_price(pr.spot, pr.strike, RISK_FREE, pr.vol, pr.T, max(pr.steps,10); npaths=8000)
                    rel = mc > 1e-8 ? abs(pr.qpe_price - mc) / mc : abs(pr.qpe_price - mc)
                    flag = rel > 0.05 ? CRED*"FLAG mismatch>5%"*CR : CGRN*"OK"*CR
                    println(@sprintf("  %-6s+%dw QPE=%.4f MC=%.4f rel=%.1f%% %s", pr.sym, pr.weeks, pr.qpe_price, mc, rel*100, flag))
                catch e
                    println("  ", pr.sym, "+", pr.weeks, " MC failed: ", e)
                end
                cnt += 1
            end
            println(CD, "Why QPE can be wrong: past-close calib lags regime shifts; model assumes GBM-like dyn vs jumps/skew; finite-diff Greeks+oracle/QPE error; no live IV/path updates. MC above is GBM sanity check under calib σ.", CR)
            println(repeat("-", 76))
            println(CD, "m+\$+q + Enter", CR)
            println(CB, BANNER, CR)
        else
            render(results_ref[], Dates.format(now(),"yyyy-mm-dd HH:MM:SS"), refresh)
            if !warmup_done[]
                println(CYEL, @sprintf("WARMUP %d/%d: calibrating (no trades) ...", warmup_count[], WARMUP_ITERS), CR)
            end
            println("m main, \$ balance, q quit (Enter)")
        end
        for _ in 1:max(1, refresh*5); quit[] && break; sleep(0.2); end
    end
    print("\e[2J\e[H"); println("MINSPM TUI exited.")
end

if abspath(PROGRAM_FILE) == @__FILE__; main(); end
