#!/usr/bin/env julia
# MINSPM Options TUI (Ncurses + UnicodePlots, immediate mode)
# Run: julia --project=. options_tui_ncurses.jl

using Printf, Dates

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
    @error "Ncurses_jll not available. Install it or run options_tui.jl instead."
    exit(1)
end
if !HAS_UNICODEPLOTS
    @error "UnicodePlots not available."
    exit(1)
end

using Ncurses_jll
using UnicodePlots

# --- Minimal Curses wrapper ---
module CursesWrap
    using Ncurses_jll

    const KEY_DOWN  = Int32(0o402)
    const KEY_UP    = Int32(0o403)
    const KEY_LEFT  = Int32(0o404)
    const KEY_RIGHT = Int32(0o405)

    const A_BOLD    = UInt32(1) << 21
    const A_REVERSE = UInt32(1) << 18
    const A_DIM     = UInt32(1) << 17

    function init()
        if !haskey(ENV, "TERMINFO_DIRS") && !haskey(ENV, "TERMINFO")
            ENV["TERMINFO_DIRS"] = joinpath(Ncurses_jll.artifact_dir, "share", "terminfo")
        end
        ccall(:setlocale, Ptr{Cchar}, (Cint, Cstring), 0, "")
        return ccall((:initscr, libncursesw), Ptr{Cvoid}, ())
    end

    endwin() = ccall((:endwin, libncursesw), Cint, ())
    cbreak() = ccall((:cbreak, libncursesw), Cint, ())
    noecho() = ccall((:noecho, libncursesw), Cint, ())
    curs_set(vis::Integer) = ccall((:curs_set, libncursesw), Cint, (Cint,), vis)
    keypad(win::Ptr{Cvoid}, bf::Bool) = ccall((:keypad, libncursesw), Cint, (Ptr{Cvoid}, Bool), win, bf)
    timeout(delay::Integer) = ccall((:timeout, libncursesw), Cvoid, (Cint,), delay)
    erase() = ccall((:erase, libncursesw), Cint, ())
    refresh() = ccall((:refresh, libncursesw), Cint, ())
    getch() = ccall((:getch, libncursesw), Cint, ())

    getmaxy(win::Ptr{Cvoid}) = ccall((:getmaxy, libncursesw), Cint, (Ptr{Cvoid},), win)
    getmaxx(win::Ptr{Cvoid}) = ccall((:getmaxx, libncursesw), Cint, (Ptr{Cvoid},), win)

    attron(attrs::UInt32) = ccall((:attron, libncursesw), Cint, (UInt32,), attrs)
    attroff(attrs::UInt32) = ccall((:attroff, libncursesw), Cint, (UInt32,), attrs)

    mvaddstr(y::Integer, x::Integer, str::AbstractString) =
        ccall((:mvaddstr, libncursesw), Cint, (Cint, Cint, Cstring), y, x, str)

    mvhline(y::Integer, x::Integer, ch::Char, n::Integer) =
        ccall((:mvhline, libncursesw), Cint, (Cint, Cint, UInt32, Cint), y, x, UInt32(ch), n)
end

using .CursesWrap

strip_ansi(s::AbstractString) = replace(s, r"\e\[[0-9;]*[a-zA-Z]" => "")

mutable struct TuiState
    running::Bool
    selected_ticker::Int
    selected_plot::Int  # 1=line,2=damped,3=bar (demo) or mapped to data views
    freq::Float64
    last_key::Int32
    message::String
end

function safe_mvaddstr(y::Int, x::Int, s::AbstractString, maxx::Int)
    if y < 0 || x < 0 || x >= maxx
        return
    end
    s2 = String(s)
    maxlen = max(0, maxx - x - 1)
    if length(s2) > maxlen
        s2 = length(s2) > 0 ? s2[1:min(maxlen, length(s2))] : ""
    end
    try
        CursesWrap.mvaddstr(y, x, s2)
    catch
    end
end

function draw_header(maxy::Int, maxx::Int)
    CursesWrap.attron(CursesWrap.A_BOLD)
    safe_mvaddstr(0, 2, "MINSPM Options TUI | [↑/↓] Menu  [←/→] Param  [q] Quit", maxx)
    CursesWrap.attroff(CursesWrap.A_BOLD)
    CursesWrap.mvhline(1, 0, '-', max(0, maxx-1))
end

function draw_menu(y::Int, x::Int, items::Vector{String}, selected::Int, maxx::Int)
    safe_mvaddstr(y, x, "View:", maxx)
    for (i, item) in enumerate(items)
        row = y + i
        if i == selected
            CursesWrap.attron(CursesWrap.A_REVERSE)
            safe_mvaddstr(row, x + 2, "> $(item) ", maxx)
            CursesWrap.attroff(CursesWrap.A_REVERSE)
        else
            safe_mvaddstr(row, x + 2, "  $(item) ", maxx)
        end
    end
end

function render_plot_view(view_idx::Int, freq::Float64, width::Int, height::Int)
    x = collect(range(0, 4π, length=100))
    plt = if view_idx == 1
        lineplot(x, sin.(freq .* x), title="Sine (f=$(round(freq,digits=2)))", width=width, height=height, border=:ascii)
    elseif view_idx == 2
        lineplot(x, cos.(freq .* x) .* exp.(-0.08 .* x), title="Damped Cosine", width=width, height=height, border=:ascii)
    else
        barplot(["Call", "Put", "Delta", "Vega"], [10*freq, 8*freq, 12, 6*freq], title="Greeks (demo)", width=width, height=height)
    end
    io = IOBuffer()
    show(IOContext(io, :color => false), plt)
    return strip_ansi(String(take!(io)))
end

function main_loop()
    stdscr = CursesWrap.init()
    CursesWrap.cbreak()
    CursesWrap.noecho()
    CursesWrap.curs_set(0)
    CursesWrap.keypad(stdscr, true)
    CursesWrap.timeout(120)

    state = TuiState(true, 1, 1, 1.0, Int32(-1), "Ready")

    menu_items = ["Price Curve", "Damped View", "Greeks Bars"]

    try
        while state.running
            ch = CursesWrap.getch()
            state.last_key = ch

            if ch == Int32('q') || ch == Int32('Q')
                state.running = false
                continue
            elseif ch == CursesWrap.KEY_UP
                state.selected_plot = max(1, state.selected_plot - 1)
            elseif ch == CursesWrap.KEY_DOWN
                state.selected_plot = min(length(menu_items), state.selected_plot + 1)
            elseif ch == CursesWrap.KEY_LEFT
                state.freq = max(0.2, state.freq - 0.1)
            elseif ch == CursesWrap.KEY_RIGHT
                state.freq = min(5.0, state.freq + 0.1)
            end

            CursesWrap.erase()
            maxy, maxx = CursesWrap.getmaxy(stdscr), CursesWrap.getmaxx(stdscr)

            draw_header(maxy, maxx)

            draw_menu(3, 2, menu_items, state.selected_plot, maxx)
            safe_mvaddstr(3 + length(menu_items) + 1, 4, "Freq/Param: $(round(state.freq,digits=2))", maxx)
            safe_mvaddstr(3 + length(menu_items) + 2, 4, "Status: $(state.message)", maxx)

            plot_x = 30
            plot_y = 3
            plot_w = max(20, maxx - plot_x - 2)
            plot_h = max(8, maxy - plot_y - 2)

            if maxx > 60 && maxy > 14
                plot_str = render_plot_view(state.selected_plot, state.freq, plot_w, plot_h)
                for (i, line) in enumerate(split(plot_str, '\n'))
                    if plot_y + i <= maxy - 2
                        safe_mvaddstr(plot_y + i, plot_x, line, maxx)
                    end
                end
            else
                safe_mvaddstr(6, plot_x, "Terminal too small for plot (need >=60x14).", maxx)
            end

            CursesWrap.refresh()
        end
    finally
        CursesWrap.endwin()
    end
end

main_loop()
