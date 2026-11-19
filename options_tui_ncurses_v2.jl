#!/usr/bin/env julia
# MINSPM Options TUI — Ncurses + UnicodePlots (immediate mode), fixed UI
# Run: julia --project=. options_tui_ncurses_v2.jl [--refresh 30] [--once]

using Printf, Dates, Statistics, LinearAlgebra, Distributions, SpecialFunctions

include(joinpath(@__DIR__, "coherent_oracle_minspm.jl"))
const CM = CoherentOracleMINSPM

const HAS_NCURSES = try
    @eval using Ncurses_jll
    true
catch
    false
end
const HAS_UNICODEPLOTS = try
    @eval using UnicodePlots
    true
catch
    false
end
if !HAS_NCURSES
    @error "Ncurses_jll required"; exit(1)
end
if !HAS_UNICODEPLOTS
    @warn "UnicodePlots not available; plots disabled"
end
using Ncurses_jll
if HAS_UNICODEPLOTS
    using UnicodePlots
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

strip_ansi(s::AbstractString) = replace(s, r"\e\[[0-9;]*[a-zA-Z]" => "")

const TICKERS = begin
    f = joinpath(@__DIR__, "tickers.txt")
    if isfile(f)
        [String(strip(l)) for l in readlines(f) if !isempty(strip(l)) && !startswith(strip(l),"#")]
    else
        ["VZ","AMZN","NVDA","AMD","INTC","QCOM","AVGO","AAPL","SPY","QQQ","PLTR"]
    end
end
const TENOR_WEEKS = [1,2,4]
const RISK_FREE = 0.0415
const TRADING_DAYS = 252

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

function py_cmd()
    for p in ("/home/skoranne/miniconda3/bin/python3","python3","python")
        try
            run(pipeline(`$p -c "import yfinance"`, stdout=devnull, stderr=devnull))
            return p
        catch; end
    end
    return "python3"
end

function fetch_snapshot_lines(tickers::AbstractVector, py::String)
    # (same as options_tui.jl - simplified call)
    syms = join(["\""*String(t)*"\"" for t in tickers], ",")
    # ... reuse logic by calling python inline in practice; keep minimal here
    return Dict{String,NamedTuple}()
end

function price_one(sym::String, spot::Float64, closes::Vector{Float64}, weeks::Int)
    T = weeks/52.0
    steps = max(10, Int(round(T*TRADING_DAYS)))
    vol, sk = length(closes) >= 30 ? CM.calibrate_1d(closes) : (0.30,-0.10)
    vol = clamp(vol,0.05,1.50); sk = clamp(sk,-0.8,-0.01)
    K = spot
    Bnd = max(0.8*spot, 2*vol*spot*sqrt(T))
    base = CM.coherent_oracle_moments(spot,K,risk_free_val(),vol,sk,T,steps; bound_b=Bnd)
    theta = asin(sqrt(clamp(base.a,1e-12,1.0-1e-12)))
    _, a_ref = CM.minspm_refined(theta, 1<<24, 4096)
    qpe_price = a_ref * Bnd * exp(-risk_free_val()*T)
    dS = max(spot*0.005,0.01)
    up = CM.coherent_oracle_moments(spot+dS,K,risk_free_val(),vol,sk,T,steps; bound_b=Bnd)
    dn = CM.coherent_oracle_moments(spot-dS,K,risk_free_val(),vol,sk,T,steps; bound_b=Bnd)
    delta = (up.price_unbounded - dn.price_unbounded)/(2*dS)
    gamma = (up.price_unbounded - 2*base.price_unbounded + dn.price_unbounded)/(dS^2)
    dσ = 0.01
    vup = CM.coherent_oracle_moments(spot,K,risk_free_val(),vol+dσ,sk,T,steps; bound_b=Bnd)
    vdn = CM.coherent_oracle_moments(spot,K,risk_free_val(),max(vol-dσ,0.01),sk,T,steps; bound_b=Bnd)
    vega = (vup.price_unbounded - vdn.price_unbounded)/(2*dσ)/100.0
    bse = bs_call(spot,K,risk_free_val(),vol,T)
    return PriceResult(sym,weeks,T,steps,spot,K,vol,sk,base.price_unbounded,base.a,
                       delta,gamma,vega,bse,qpe_price)
end
risk_free_val() = RISK_FREE

mutable struct AppState
    running::Bool
    results::Vector{PriceResult}
    message::String
    last_update::String
end

function safe_mv(y,x,s,maxx)
    if y < 0 || x < 0 || x >= maxx; return; end
    s2 = String(s)
    m = max(0, maxx - x - 1)
    if length(s2) > m
        s2 = length(s2) > m ? s2[1:m] : s2
    end
    try CursesWrap.mvaddstr(y,x,s2) catch; end
end

function draw_ui(stdscr, state::AppState)
    CursesWrap.erase()
    maxy, maxx = CursesWrap.getmaxy(stdscr), CursesWrap.getmaxx(stdscr)
    CursesWrap.attron(CursesWrap.A_BOLD)
    safe_mv(0,2,"MINSPM Ncurses TUI (fixed, immediate mode) | [q] Quit [r] Refresh", maxx)
    CursesWrap.attroff(CursesWrap.A_BOLD)
    CursesWrap.mvhline(1,0,'-', max(0,maxx-1))
    safe_mv(2,2,"Updated: $(state.last_update)  $(state.message)", maxx)

    row = 4
    safe_mv(row,2,@sprintf("%-6s %2s %9s %6s %8s %8s %8s","Sym","Tw","Spot","Vol%","Asian","Delta","EuroBS"), maxx)
    row += 1
    CursesWrap.mvhline(row,2,'-', min(70, maxx-4)); row += 1
    for pr in state.results[1:min(30, length(state.results))]
        safe_mv(row,2,@sprintf("%-6s %+1dw %9.2f %5.1f %8.3f %8.3f %8.3f",
            pr.sym,pr.weeks,pr.spot,pr.vol*100,pr.price,pr.delta,pr.bs_euro), maxx)
        row += 1
        if row >= maxy - 8; break; end
    end

    # Simple plot panel
    if HAS_UNICODEPLOTS && maxx > 80 && maxy > 20 && !isempty(state.results)
        px = 50; py = 4
        # plot Asian prices by index (demo)
        vals = [pr.price for pr in state.results[1:min(20, length(state.results))]]
        lbls = ["$(pr.sym)$(pr.weeks)" for pr in state.results[1:min(20, length(state.results))]]
        try
            plt = barplot(lbls, vals, title="Asian Prices (sample)", width=min(30,maxx-px-2), height=min(12,maxy-py-3))
            io = IOBuffer()
            show(IOContext(io,:color=>false), plt)
            splot = strip_ansi(String(take!(io)))
            for (i,line) in enumerate(split(splot,'\n'))
                if py+i <= maxy-2
                    safe_mv(py+i, px, line, maxx)
                end
            end
        catch e
            safe_mv(py, px, "plot error", maxx)
        end
    end

    CursesWrap.refresh()
end

function main()
    stdscr = CursesWrap.init()
    CursesWrap.cbreak()
    CursesWrap.noecho()
    CursesWrap.curs_set(0)
    CursesWrap.keypad(stdscr,true)
    CursesWrap.timeout(300)

    state = AppState(true, PriceResult[], "Starting...", Dates.format(now(),"HH:MM:SS"))

    try
        while state.running
            ch = CursesWrap.getch()
            if ch == Int32('q') || ch == Int32('Q')
                state.running = false
                break
            elseif ch == Int32('r') || ch == Int32('R')
                state.message = "Refresh not implemented in this minimal build"
            end
            draw_ui(stdscr, state)
        end
    finally
        CursesWrap.endwin()
    end
    println("Exited.")
end

main()
