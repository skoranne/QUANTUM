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

const HAS_NCURSES = try
    @eval using Ncurses_jll
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

module CursesWrap
    using Ncurses_jll
    const KEY_DOWN = Int32(0o402); const KEY_UP = Int32(0o403)
    const KEY_LEFT = Int32(0o404); const KEY_RIGHT = Int32(0o405)
    const A_BOLD = UInt32(1)<<21; const A_REVERSE = UInt32(1)<<18
    function init()
        if !haskey(ENV,"TERMINFO_DIRS") && !haskey(ENV,"TERMINFO")
            ENV["TERMINFO_DIRS"] = joinpath(Ncurses_jll.artifact_dir,"share","terminfo")
        end
        ccall(:setlocale, Ptr{Cchar}, (Cint,Cstring), 0, "")
        return ccall((:initscr, libncursesw), Ptr{Cvoid}, ())
    end
    endwin() = ccall((:endwin, libncursesw), Cint, ())
    cbreak() = ccall((:cbreak, libncursesw), Cint, ())
    noecho() = ccall((:noecho, libncursesw), Cint, ())
    curs_set(v::Integer) = ccall((:curs_set, libncursesw), Cint, (Cint,), v)
    keypad(w::Ptr{Cvoid}, b::Bool) = ccall((:keypad, libncursesw), Cint, (Ptr{Cvoid},Bool), w, b)
    timeout(d::Integer) = ccall((:timeout, libncursesw), Cvoid, (Cint,), d)
    erase() = ccall((:erase, libncursesw), Cint, ())
    refresh() = ccall((:refresh, libncursesw), Cint, ())
    getch() = ccall((:getch, libncursesw), Cint, ())
    getmaxy(w::Ptr{Cvoid}) = ccall((:getmaxy, libncursesw), Cint, (Ptr{Cvoid},), w)
    getmaxx(w::Ptr{Cvoid}) = ccall((:getmaxx, libncursesw), Cint, (Ptr{Cvoid},), w)
    attron(a::UInt32) = ccall((:attron, libncursesw), Cint, (UInt32,), a)
    attroff(a::UInt32) = ccall((:attroff, libncursesw), Cint, (UInt32,), a)
    mvaddstr(y,x,s::AbstractString) = ccall((:mvaddstr, libncursesw), Cint, (Cint,Cint,Cstring), y,x,s)
    mvhline(y,x,ch::Char,n::Integer) = ccall((:mvhline, libncursesw), Cint, (Cint,Cint,UInt32,Cint), y,x,UInt32(ch),n)
end
using .CursesWrap

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

function ncurses_safe_mv(y,x,s,maxx)
    if y < 0 || x < 0 || x >= maxx; return; end
    s2 = String(s)
    m = max(0, maxx - x - 1)
    if length(s2) > m
        s2 = length(s2) > 0 ? s2[1:min(m, end)] : ""
    end
    try CursesWrap.mvaddstr(y,x,s2) catch; end
end

function main()
    if !HAS_NCURSES
        @error "Ncurses_jll not available. This TUI now requires Ncurses. Install Ncurses_jll or use an older version."
        return
    end
    refresh = 30
    i = 1
    while i <= length(ARGS)
        if ARGS[i] == "--refresh" && i < length(ARGS); refresh = parse(Int,ARGS[i+1]); i+=2
        elseif ARGS[i] == "--once"
            # fallback simple print once
            PYTHON[] = py_cmd()
            snapshot = fetch_snapshot_lines(TICKERS)
            results = compute_all(snapshot)
            println("MINSPM (once): ", Dates.format(now(),"yyyy-mm-dd HH:MM:SS"))
            for pr in results
                @printf("%s %+dw spot=%.2f vol=%.1f%% asian=%.3f delta=%.3f euro=%.3f\n", pr.sym, pr.weeks, pr.spot, pr.vol*100, pr.price, pr.delta, pr.bs_euro)
            end
            return
        else; i+=1; end
    end
    PYTHON[] = py_cmd()
    stdscr = CursesWrap.init()
    CursesWrap.cbreak()
    CursesWrap.noecho()
    CursesWrap.curs_set(0)
    CursesWrap.keypad(stdscr,true)
    CursesWrap.timeout(250)

    view_balance = false
    quit = false
    results_ref = PriceResult[]
    market_opts_ref = Dict{String,Dict{Int,Any}}()
    last_fetch = time() - 1e9

    try
        while !quit
            ch = CursesWrap.getch()
            if ch == Int32('q') || ch == Int32('Q')
                quit = true; break
            elseif ch == Int32('$') || ch == Int32('b') || ch == Int32('B')
                view_balance = !view_balance
            elseif ch == Int32('m') || ch == Int32('M')
                view_balance = false
            elseif ch == Int32('r') || ch == Int32('R')
                # force refresh
                last_fetch = time() - 1e9
            end

            # periodic fetch
            if time() - last_fetch >= refresh
                try
                    snap = fetch_snapshot_lines(TICKERS)
                    results_ref = compute_all(snap)
                catch e
                    @warn "fetch/compute failed" exception=e
                end
                try
                    market_opts_ref = fetch_option_market_lines(TICKERS)
                catch e
                    market_opts_ref = Dict{String,Dict{Int,Any}}()
                end
                if !warmup_done[]
                    global warmup_count
                    warmup_count[] += 1
                    if warmup_count[] >= WARMUP_ITERS
                        warmup_done[] = true
                    end
                end
                last_fetch = time()
            end

            CursesWrap.erase()
            maxy, maxx = CursesWrap.getmaxy(stdscr), CursesWrap.getmaxx(stdscr)
            CursesWrap.attron(CursesWrap.A_BOLD)
            ncurses_safe_mv(0,2,"MINSPM TUI (Ncurses, fixed) | [q] quit [\$] balance [m] main [r] refresh", maxx)
            CursesWrap.attroff(CursesWrap.A_BOLD)
            CursesWrap.mvhline(1,0,'-', max(0,maxx-1))
            ncurses_safe_mv(2,2,"$(Dates.format(now(),"yyyy-mm-dd HH:MM:SS")) refresh=$(refresh)s r=$(RISK_FREE)", maxx)

            if !view_balance
                row = 4
                ncurses_safe_mv(row,2,@sprintf("%-6s %2s %9s %6s %8s %8s %8s","Sym","Tw","Spot","Vol%","Asian","Delta","EuroBS"), maxx)
                row += 2
                for pr in results_ref[1:min(length(results_ref), max(0,maxy-row-2))]
                    ncurses_safe_mv(row,2,@sprintf("%-6s %+1dw %9.2f %5.1f %8.3f %8.3f %8.3f",
                        pr.sym,pr.weeks,pr.spot,pr.vol*100,pr.price,pr.delta,pr.bs_euro), maxx)
                    row += 1
                end
                if warmup_done[] == false
                    ncurses_safe_mv(maxy-2,2,"WARMUP $(warmup_count[])/$(WARMUP_ITERS) ... calibrating", maxx)
                end
            else
                row = 3
                ncurses_safe_mv(row,2,"BALANCE / market-ref view (Ncurses)", maxx); row += 1
                if !warmup_done[]
                    ncurses_safe_mv(row,2,"WARMUP $(warmup_count[])/$(WARMUP_ITERS): no trades opened", maxx); row += 1
                else
                    ncurses_safe_mv(row,2,@sprintf("Virtual \$%.2f | cost \$%.2f | edge>=\$%.2f", virtual_cash[], FIXED_TRADE_COST, FIXED_TRADE_COST+TRADE_EDGE_MIN), maxx); row += 1
                end
                row += 1
                unreal = 0.0
                shown = 0
                for pr in results_ref
                    mopt = get(get(market_opts_ref, pr.sym, Dict{Int,Any}()), pr.weeks, nothing)
                    if mopt === nothing || !isfinite(get(mopt,:price,NaN)) || mopt.price <= 0
                        continue
                    end
                    dgap = pr.delta - get(mopt,:delta,NaN)
                    edge = Float64(mopt.price) - pr.qpe_price
                    key = (pr.sym, pr.weeks)
                    pos = paper_positions[]
                    if haskey(pos, key)
                        upl = Float64(mopt.price) - pos[key].entry_mkt
                        unreal += upl
                        line = @sprintf("%-6s+%dw dΔ=%+.3f mkt=%.2f ent=%.2f uP/L=%+.2f", pr.sym,pr.weeks,dgap,Float64(mopt.price),pos[key].entry_mkt,upl)
                    else
                        status = warmup_done[] ? (edge >= (FIXED_TRADE_COST+TRADE_EDGE_MIN) ? "edge_ok" : "no_trade") : "warmup"
                        line = @sprintf("%-6s+%dw dΔ=%+.3f mkt=%.2f QPE=%.2f edge=%.2f %s", pr.sym,pr.weeks,dgap,Float64(mopt.price),pr.qpe_price,edge,status)
                    end
                    if row < maxy - 3
                        ncurses_safe_mv(row,2,line, maxx)
                        row += 1
                        shown += 1
                    end
                end
                row += 1
                equity = virtual_cash[] + unreal
                ncurses_safe_mv(row,2,@sprintf("Cash \$%.2f | Unreal \$%.2f | Equity \$%.2f (shown %d)", virtual_cash[], unreal, equity, shown), maxx); row += 1
                # try to show a small unicode plot if space
                if HAS_UNICODEPLOTS && row < maxy - 6 && shown > 0
                    try
                        upl_labels = String[]; upl_vals = Float64[]
                        for pr in results_ref
                            key = (pr.sym,pr.weeks)
                            mopt = get(get(market_opts_ref, pr.sym, Dict{Int,Any}()), pr.weeks, nothing)
                            if mopt !== nothing && isfinite(get(mopt,:price,NaN)) && haskey(paper_positions[], key)
                                push!(upl_labels, "$(pr.sym)+$(pr.weeks)")
                                push!(upl_vals, Float64(mopt.price) - paper_positions[][key].entry_mkt)
                            end
                        end
                        if !isempty(upl_vals)
                            plt = UnicodePlots.barplot(upl_labels, upl_vals, title="uP/L \$", width=min(40, maxx-4), height=min(8, maxy-row-2))
                            io = IOBuffer(); show(IOContext(io,:color=>false), plt)
                            splot = replace(String(take!(io)), r"\e\[[0-9;]*[a-zA-Z]" => "")
                            for ln in split(splot, '\n')
                                if row < maxy - 1
                                    ncurses_safe_mv(row,2,ln,maxx); row += 1
                                end
                            end
                        end
                    catch e
                        ncurses_safe_mv(row,2,"(plot unavailable)",maxx)
                    end
                end
            end
            CursesWrap.refresh()
        end
    finally
        CursesWrap.endwin()
    end
    println("MINSPM TUI exited.")
end

if abspath(PROGRAM_FILE) == @__FILE__; main(); end
