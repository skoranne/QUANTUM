#!/usr/bin/env julia
# MINSPM Option Pricing TUI — ATM Arithmetic Asian, +2w/+4w, VZ AMZN NVDA AMD INTC QCOM AVGO AAPL
# Engine: self-sufficient coherent MINSPM moment oracle + MINSPM QCS QPE (deterministic, no MC sampling)
# Run: julia --project=. options_tui.jl [--refresh 30] [--once]

using Printf, Dates, Statistics, LinearAlgebra, Distributions, SpecialFunctions
include(joinpath(@__DIR__, "coherent_oracle_minspm.jl"))
const CM = CoherentOracleMINSPM

const TICKERS = ["VZ","AMZN","NVDA","AMD","INTC","QCOM","AVGO","AAPL"]
const TENOR_WEEKS = [2, 4]
const RISK_FREE = 0.0415
const TRADING_DAYS = 252

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

function fetch_snapshot_lines(tickers::Vector{String})
    py = PYTHON[]
    syms = join(["\"" * t * "\"" for t in tickers], ",")
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
    println(CB,CCYN,"  MINSPM Option Pricing TUI — ATM Arithmetic Asian (coherent, self-sufficient, no MC)",CR,
            CD,"  ",last_update,"  refresh:",refresh,"s  [q]uit",CR)
    println(CD,"  Engine: MINSPM moment oracle+QPE (Z_{2^8}) | Greeks: bump-and-reprice | r=$(RISK_FREE) | +2w/+4w",CR)
    println(repeat("─", 160))
    @printf("%s %-6s │ %4s │ %10s │ %10s │ %7s │ %10s │ %9s │ %10s │ %9s │ %10s │ %s%s\n",
            CB,"Sym","Ten","Spot","Strike","Vol","Asian","Delta","Gamma","Vega/1vp","EuroBS","Signals",CR)
    println(repeat("─", 160))
    bysym = Dict{String,Vector{PriceResult}}()
    for pr in results; push!(get!(bysym, pr.sym, PriceResult[]), pr); end
    for sym in TICKERS
        prs = sort(get(bysym, sym, PriceResult[]), by=x->x.weeks); isempty(prs) && continue
        for (i,pr) in enumerate(prs)
            sigs = analyze(pr)
            hasarb = any(s.level==:arb for s in sigs); hasop = any(s.level==:oppty for s in sigs)
            col = hasarb ? CGRN : (hasop ? CYEL : "")
            symc = i==1 ? sym : ""; spotc = i==1 ? @sprintf("%.2f",pr.spot) : ""; volc = i==1 ? @sprintf("%.1f%%",pr.vol*100) : ""
            @printf("%s %-6s │ +%2dw │ %10s │ %10.2f │ %7s │ %10.4f │ %9.4f │ %10.6f │ %9.5f │ %10.4f │ %s",
                    col,symc,pr.weeks,spotc,pr.strike,volc,pr.price,pr.delta,pr.gamma,pr.vega,pr.bs_euro,CR)
            println(join([sigstr(s) for s in sigs], "; "))
        end
        println(CD, repeat("·", 160), CR)
    end
    println()
    println(CB,"  Convexity & Hedge Monitor",CR, CD,"  (gamma P&L per 1% spot move = 0.5*Γ*S²*0.0001; vega per 1 vol pt)",CR)
    println(repeat("─", 140))
    @printf("%s %-6s  %-5s  %10s  %12s  %12s  %12s  %12s  %s%s\n", CB,"Ticker","Ten","Spot","Delta","Gamma","Γ-PnL/1%%","Vega/1vp","Hedge (delta-neutral shares / gamma note)",CR)
    for sym in TICKERS
        for pr in sort(get(bysym,sym,PriceResult[]), by=x->x.weeks)
            conv = 0.5 * pr.gamma * pr.spot^2 * 0.0001
            hedge = -pr.delta
            note = @sprintf("short %.3f sh/opt for Δ-neutral; convexity \$%.4f per 1%% move", hedge, conv)
            flag = abs(conv) > 0.05 ? CRED*" [HIGH Γ]"*CR : ""
            @printf(" %-6s  +%2dw  %10.2f  %12.5f  %12.6f  %12.6f  %12.6f  %s%s\n",
                    sym, pr.weeks, pr.spot, pr.delta, pr.gamma, conv, pr.vega, note, flag)
        end
    end
    println(repeat("─", 140))
    println(CD,"  Signals: ARB/OPPTY = model vs Euro-BS reference (Asian fair < Euro typical). Verify vs executable market quotes before trading.",CR)
    println(CD,"  Convexity hedge: monitor Γ (0.5*Γ*ΔS²). Vega: sensitivity per 1 vol point. MINSPM oracle is deterministic/self-sufficient.",CR)
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
    println("MINSPM TUI starting (python=$(PYTHON[]), refresh=$(refresh)s)... fetching data...")
    if once
        snapshot = fetch_snapshot_lines(TICKERS)
        results = compute_all(snapshot)
        render(results, Dates.format(now(),"yyyy-mm-dd HH:MM:SS"), refresh)
        return
    end
    quit = Ref(false)
    @async begin
        try
            while !quit[]
                if bytesavailable(stdin) > 0
                    c = lowercase(Char(read(stdin, UInt8)))
                    c == 'q' && (quit[] = true; break)
                end
                sleep(0.2)
            end
        catch; quit[] = true; end
    end
    while !quit[]
        snapshot = fetch_snapshot_lines(TICKERS)
        results = compute_all(snapshot)
        render(results, Dates.format(now(),"yyyy-mm-dd HH:MM:SS"), refresh)
        for _ in 1:max(1, refresh*5); quit[] && break; sleep(0.2); end
    end
    print("\e[2J\e[H"); println("MINSPM TUI exited.")
end

if abspath(PROGRAM_FILE) == @__FILE__; main(); end
