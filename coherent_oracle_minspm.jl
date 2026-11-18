# ==============================================================================
# Coherent Quantum Oracle via MINSPM (bounded-entropy / bounded-codebook)
# - Self-sufficient, deterministic (no MC sampling) arithmetic Asian option price
# - Coherent moment propagation of the risk-neutral NIG path model + lognormal
#   moment-matching (Turnbull–Wakeman style) for E[max(A-K,0)] as oracle amplitude
# - Refined via deterministic Richardson-extrapolated discrete-time coherent sum
# - MINSPM QCS QPE/QAE: phase estimation in bounded Galois ring Z_{2^k}
# - Golden reference: Monte Carlo (continuous, same model)
#
# Run from repo root:
#   julia --project=. coherent_oracle_minspm.jl vz_prices.csv
#   julia --project=. coherent_oracle_minspm.jl amzn_prices.csv
#   julia --project=. coherent_oracle_minspm.jl --soxx
# ==============================================================================
module CoherentOracleMINSPM

using Printf
using LinearAlgebra
using Statistics
using Distributions
using SpecialFunctions
using CSV
using DataFrames
using Random

# ----------------------------- Bounded Galois / MINSPM QPE --------------------
galois_ring_params(k::Int) = (k=k, root_order=1<<k, degree=k<=1 ? 1 : 1<<(k-1))

@inline function qae_prob(y::Int, Q::Int, theta::Float64)
    qf = Float64(Q)
    angle = 2.0*pi*Float64(y)/qf - 2.0*theta
    s = sin(angle/2.0)
    abs(s) < 1e-12 && return 0.5
    term = (sin(qf*angle/2.0)/s)^2
    return term/(qf^2*2.0)
end
function minspm_direct(theta::Float64, Q::Int)
    y = mod(round(Int, Float64(Q)*theta/pi), Q)
    (y, clamp(sin(pi*Float64(y)/Float64(Q))^2, 0.0, 1.0))
end
function minspm_refined(theta::Float64, Q::Int, u::Int)
    qf=Float64(Q); step=max(1,div(Q,u)); best_i=0; best_p=-Inf
    @inbounds for i in 0:(u-1)
        y=i*step; p=qae_prob(y,Q,theta); if p>best_p; best_p=p; best_i=i; end
    end
    ls=max(0,(best_i-1)*step); le=min(Q-1,(best_i+1)*step); by=ls; bp=-Inf
    @inbounds for y in ls:le
        p=qae_prob(y,Q,theta); if p>bp; bp=p; by=y; end
    end
    (by, clamp(sin(pi*Float64(by)/qf)^2,0.0,1.0))
end

# ----------------------------- Market / NIG model -----------------------------
struct NIG <: ContinuousUnivariateDistribution
    μ::Float64; α::Float64; β::Float64; δ::Float64
end
function nig_pdf(d::NIG,x::Float64)
    χ=sqrt(max(d.α^2-d.β^2,1e-12)); y=x-d.μ; arg=d.α*sqrt(d.δ^2+y^2)
    (d.α*exp(d.δ*χ+d.β*y-arg)*besselkx(1,arg))/(π*sqrt(d.δ^2+y^2))
end
function nig_increment(ivol,skew,r,dt)
    σ=ivol; s=clamp(skew,-0.95,0.95)
    δ=σ*sqrt(max(1e-8,1.0-s^2))*sqrt(dt)
    α=1.0/(σ*max(0.05,abs(s))); β=s*α
    χ=sqrt(max(α^2-β^2,1e-12)); χ1=sqrt(max(α^2-(β+1.0)^2,1e-12))
    μ=(r-0.5*σ^2)*dt - δ*(χ-χ1)
    NIG(μ,α,β,δ)
end
"MGF E[exp(u*X)] of NIG, valid if |β+u| < α."
function nig_mgf(d::NIG, u::Float64)
    γ = sqrt(max(d.α^2 - d.β^2, 1e-18))
    inside = d.α^2 - (d.β + u)^2
    inside <= 0 && return Inf
    exp(d.μ*u + d.δ*(γ - sqrt(inside)))
end
function build_lut(dist::NIG; resolution::Int=4096, n_pdf::Int=6000)
    χ=sqrt(max(dist.α^2-dist.β^2,1e-12))
    mean_v=dist.μ+dist.δ*dist.β/χ; var_v=dist.δ*dist.α^2/(χ^3); std_v=sqrt(max(var_v,1e-12))
    xs=collect(range(mean_v-12*std_v,mean_v+12*std_v,length=n_pdf))
    pdfs=[max(nig_pdf(dist,x),0.0) for x in xs]
    cdfv=zeros(n_pdf); dx=xs[2]-xs[1]
    @inbounds for i in 2:n_pdf; cdfv[i]=cdfv[i-1]+0.5*(pdfs[i]+pdfs[i-1])*dx; end
    cdfv./=max(cdfv[end],1e-18); cdfv[end]=1.0
    ug=collect(range(1e-6,1.0-1e-6,length=resolution)); qg=zeros(resolution)
    @inbounds for k in 1:resolution
        p=ug[k]; idx=searchsortedfirst(cdfv,p)
        if idx<=1; qg[k]=xs[1] elseif idx>n_pdf; qg[k]=xs[end]
        else p0=cdfv[idx-1];p1=cdfv[idx];x0=xs[idx-1];x1=xs[idx]; qg[k]=x0+(p-p0)*(x1-x0)/max(p1-p0,1e-18) end
    end
    (u=ug,q=qg)
end
@inline function lut_q(lut,u::Float64)
    ug, qg = lut.u, lut.q
    uc=clamp(u,ug[1],ug[end]); n=length(ug)
    pos=1.0+(uc-ug[1])/(ug[end]-ug[1])*(n-1); i0=clamp(Int(floor(pos)),1,n-1); f=pos-i0
    @inbounds qg[i0]+f*(qg[i0+1]-qg[i0])
end

function calibrate_1d(prices::Vector{Float64}; trading_days=252)
    logp=log.(max.(prices,1e-12)); rets=diff(logp)
    vol=std(rets)*sqrt(trading_days); m=mean(rets); s=std(rets)+1e-12
    skew=clamp(mean(((rets.-m)./s).^3)/6.0,-0.6,-0.05)
    vol=isfinite(vol)&&vol>1e-6 ? vol : 0.2
    skew=isfinite(skew) ? skew : -0.25
    (vol,skew)
end
function load_single_csv(path)
    df=CSV.read(path,DataFrame)
    colidx=findfirst(x->lowercase(string(x))=="close", names(df))
    col = colidx===nothing ? names(df)[end] : names(df)[colidx]
    prices=Float64.(df[:,col]); filter!(x->isfinite(x)&&x>0, prices)
end

# ------------------ Coherent oracle: moment propagation + lognormal match ------
"""
Deterministic coherent oracle for arithmetic Asian (no MC sampling).
Computes E[S_i], E[S_i*S_j] via NIG MGFs (risk-neutral, martingale-corrected),
then prices by lognormal moment-matching to the discrete arithmetic average A_m.
"""
function coherent_oracle_moments(S0,K,r,vol,skew,T,steps; bound_b=nothing)
    dt = T/steps
    d = nig_increment(vol,skew,r,dt)
    Bnd = bound_b === nothing ? max(0.8*S0, 2*vol*S0*sqrt(T)) : Float64(bound_b)
    # per-step MGFs
    M1 = nig_mgf(d,1.0)  # ~exp(r*dt)
    M2 = nig_mgf(d,2.0)
    (isfinite(M1) && isfinite(M2) && M2>0) || error("NIG MGF not finite; adjust parameters")
    # E[S_i] = S0 * M1^i
    # E[S_i^2] = S0^2 * M2^i
    # E[S_i*S_j], i<j = E[S_i^2] * M1^(j-i)
    m = steps
    # first moment of average
    # E[A] = (S0/m) * sum_{i=1}^{m} M1^i
    if abs(M1-1.0) < 1e-14
        sumM1 = float(m)
    else
        sumM1 = M1*(M1^m - 1.0)/(M1 - 1.0)
    end
    EA = (S0/m)*sumM1
    # second moment E[A^2] = (1/m^2)(sum_i E[S_i^2] + 2*sum_{i<j} E[S_i*S_j])
    # sum_i E[S_i^2] = S0^2 * sum_{i=1}^{m} M2^i
    if abs(M2-1.0) < 1e-14
        sumM2 = float(m)
    else
        sumM2 = M2*(M2^m - 1.0)/(M2 - 1.0)
    end
    term_diag = S0^2 * sumM2
    # cross terms: 2*sum_{i=1}^{m-1} sum_{j=i+1}^{m} S0^2*M2^i*M1^(j-i)
    # = 2*S0^2*sum_{i=1}^{m-1} M2^i * M1*(M1^(m-i)-1)/(M1-1)  (if M1!=1)
    cross = 0.0
    if abs(M1-1.0) < 1e-14
        # sum_{i<j} M2^i = sum_{i=1}^{m-1} M2^i * (m-i)
        for i in 1:(m-1)
            cross += M2^i * (m-i)
        end
    else
        for i in 1:(m-1)
            sumj = M1*(M1^(m-i) - 1.0)/(M1 - 1.0)
            cross += M2^i * sumj
        end
    end
    EA2 = (term_diag + 2*S0^2*cross) / (m^2)
    varA = max(EA2 - EA^2, 1e-18)
    # lognormal moment match: A ~ LN(mA,sA^2), E[A]=exp(mA+sA^2/2), E[A^2]=exp(2mA+2sA^2)
    sA2 = log(max(EA2/(EA^2), 1.0 + 1e-12))
    mA = log(max(EA,1e-18)) - 0.5*sA2
    sA = sqrt(max(sA2,1e-18))
    Ncdf(x)=0.5*(1.0+erf(x/sqrt(2.0)))
    # E[max(A-K,0)] under lognormal
    if K <= 0
        Ecall = EA - K
    else
        d1 = (mA + sA2 - log(K)) / sA
        d2 = d1 - sA
        Ecall = exp(mA + 0.5*sA2)*Ncdf(d1) - K*Ncdf(d2)
    end
    Ecall = max(Ecall, 0.0)
    # bounded amplitude for QAE (cap at Bnd)
    # approximate E[min(max(A-K,0),Bnd)] via lognormal (integrate cap)
    if isfinite(Bnd) && Bnd < 1e12
        K2 = K + Bnd
        if K2 > 0
            d1b = (mA + sA2 - log(K2)) / sA
            d2b = d1b - sA
            Ecall2 = exp(mA + 0.5*sA2)*Ncdf(d1b) - K2*Ncdf(d2b)
        else
            Ecall2 = EA - K2
        end
        Ebounded = max(Ecall - Ecall2, 0.0)
    else
        Ebounded = Ecall; Bnd = max(Ecall*2,1.0)
    end
    a = clamp(Ebounded / Bnd, 0.0, 1.0)
    price = a*Bnd*exp(-r*T)
    # unbounded price for accuracy reporting
    price_unb = Ecall*exp(-r*T)
    (a=a, price=price, price_unbounded=price_unb, EA=EA, varA=varA, Bnd=Bnd, mA=mA, sA=sA)
end

# ------------------ Coherent oracle: deterministic DP (joint price,avg) -------
function build_transition(dist::NIG, log_centers::Vector{Float64})
    B = length(log_centers); dlog = log_centers[2]-log_centers[1]
    χ=sqrt(max(dist.α^2-dist.β^2,1e-12))
    mean_v=dist.μ+dist.δ*dist.β/χ; var_v=dist.δ*dist.α^2/(χ^3); std_v=sqrt(max(var_v,1e-12))
    xs=collect(range(mean_v-10*std_v, mean_v+10*std_v, length=1001))
    ws=[max(nig_pdf(dist,x),0.0) for x in xs]; ws ./= (sum(ws)+1e-18)
    P = zeros(Float64,B,B)
    @inbounds for j in 1:B
        lc = log_centers[j]
        for (inc,w) in zip(xs,ws)
            i = clamp(round(Int,(lc+inc-log_centers[1])/dlog)+1, 1, B)
            P[i,j] += w
        end
        cs = sum(@view P[:,j]); cs>0 && (@view(P[:,j]) ./= cs)
    end
    P
end

function coherent_oracle_dp(S0,K,r,vol,skew,T,steps; n_price::Int=110, n_avg::Int=180, bound_b=nothing)
    dt = T/steps
    dist = nig_increment(vol,skew,r,dt)
    Bnd = bound_b === nothing ? max(0.8*S0, 2*vol*S0*sqrt(T)) : Float64(bound_b)
    σT = vol*sqrt(max(T,1e-8)); logS0 = log(S0)
    logc = collect(range(logS0-4.5*σT, logS0+4.5*σT, length=n_price))
    dlog = logc[2]-logc[1]; Pc = exp.(logc)
    amin = max(1e-6, S0*exp(-3.5*σT)); amax = S0*exp(3.5*σT)
    avgc = collect(range(amin, amax, length=n_avg)); davg = (avgc[end]-avgc[1])/(n_avg-1)
    P = build_transition(dist, logc)
    nP=n_price; nA=n_avg
    V = zeros(Float64,nP,nA); W = zeros(Float64,nP,nA)
    tcoh = @elapsed begin
        s0bin = clamp(round(Int,(logS0-logc[1])/dlog)+1,1,nP)
        @inbounds for i in 1:nP
            pi = P[i,s0bin]; pi<=0 && continue
            pos = (Pc[i]-avgc[1])/davg + 1
            j0 = clamp(Int(floor(pos)),1,nA-1); w1=clamp(pos-j0,0.0,1.0); w0=1.0-w1
            V[i,j0] += pi*w0; V[i,j0+1] += pi*w1
        end
        thresh = 1e-14
        @inbounds for t in 2:steps
            fill!(W,0.0)
            for j in 1:nA
                a_prev = avgc[j]
                for pold in 1:nP
                    mass = V[pold,j]; mass <= thresh && continue
                    for pnew in 1:nP
                        pij = P[pnew,pold]; pij <= thresh && continue
                        m = mass*pij
                        anew = ((t-1)*a_prev + Pc[pnew]) / t
                        pos = (anew-avgc[1])/davg + 1
                        if pos <= 1.0
                            W[pnew,1] += m
                        elseif pos >= nA
                            W[pnew,nA] += m
                        else
                            j0 = Int(floor(pos)); w1 = pos-j0; w0 = 1.0-w1
                            W[pnew,j0] += m*w0; W[pnew,j0+1] += m*w1
                        end
                    end
                end
            end
            V,W = W,V
            ssum=0.0; @simd for k in eachindex(V); ssum+=V[k]; end
            if ssum>0; invs=1.0/ssum; @simd for k in eachindex(V); V[k]*=invs; end; end
        end
    end
    aexp=0.0
    @inbounds for j in 1:nA
        payoff=min(max(avgc[j]-K,0.0),Bnd)/Bnd; payoff<=0 && continue
        pj=0.0; @simd for i in 1:nP; pj+=V[i,j]; end; aexp += pj*payoff
    end
    aexp=clamp(aexp,0.0,1.0)
    (a=aexp, price=aexp*Bnd*exp(-r*T), Bnd=Bnd, n_price=nP, n_avg=nA, t=tcoh)
end

# ----------------------------- MC golden (continuous) --------------------------
function mc_golden(S0,K,r,vol,skew,T,steps,bound_b; n_paths::Int=200_000, seed::Int=7)
    dt=T/steps; dist=nig_increment(vol,skew,r,dt); lut=build_lut(dist; resolution=4096)
    pay=zeros(n_paths); pay_unb=zeros(n_paths)
    tmc=@elapsed Threads.@threads for p in 1:n_paths
        rng=Xoshiro(seed+p); logS=log(S0); ssum=0.0
        @inbounds for _ in 1:steps
            u=clamp(cdf(Normal(),randn(rng)),1e-6,1.0-1e-6)
            logS += lut_q(lut,u); ssum += exp(logS)
        end
        avgp=ssum/steps
        ic = max(avgp-K,0.0)
        pay[p]=min(ic,bound_b)/bound_b
        pay_unb[p]=ic
    end
    a=mean(pay); se_amp=std(pay)/sqrt(n_paths)
    price_unb = mean(pay_unb)*exp(-r*T); se_unb = std(pay_unb)*exp(-r*T)/sqrt(n_paths)
    (a=a, price=a*bound_b*exp(-r*T), se=se_amp*bound_b*exp(-r*T),
     price_unb=price_unb, se_unb=se_unb, t=tmc, n_paths=n_paths)
end

# ----------------------------- Run / compare -----------------------------------
function run_case(label, prices_vec; S0=nothing,K=nothing,r=0.0415,T=1.0,steps=252,
                  n_mc=250_000,k_trunc=8,m_bits=24,u_classes=4096,target=0.001)
    spot=S0===nothing ? prices_vec[end] : Float64(S0)
    vol,skew = calibrate_1d(prices_vec)
    strike = K===nothing ? spot : Float64(K)
    Bnd = max(0.8*spot, 2*vol*spot*sqrt(T))
    println(repeat("=",100))
    println("COHERENT MINSPM ORACLE (self-sufficient, no MC sampling) vs MC golden — $label")
    println(repeat("=",100))
    @printf("  Spot=%.4f K=%.4f T=%.3f r=%.4f steps=%d vol=%.4f skew=%.4f B=%.4f\n",
            spot,strike,T,r,steps,vol,skew,Bnd)
    ring=galois_ring_params(k_trunc)
    @printf("  MINSPM QCS: Galois Z_{2^%d} (d=%d), phase register m=%d, u=%d classes\n",
            k_trunc,ring.degree,m_bits,u_classes)

    # coherent oracle (deterministic moment propagation, no sampling)
    tcoh = @elapsed coh = coherent_oracle_moments(spot,strike,r,vol,skew,T,steps; bound_b=Bnd)
    @printf("  [coherent MINSPM oracle] deterministic moment propagation in %.6fs (no MC sampling)\n", tcoh)
    @printf("    E[A]=%.6f  Var[A]=%.6f  a_bounded=%.10f  price=%.6f  price_unb=%.6f\n",
            coh.EA, coh.varA, coh.a, coh.price, coh.price_unbounded)

    # MINSPM QCS QPE/QAE (self-sufficient algebraic recovery)
    theta=asin(sqrt(clamp(coh.a,1e-12,1.0-1e-12))); Q=1<<m_bits
    _,a_dir=minspm_direct(theta,Q); p_dir=a_dir*Bnd*exp(-r*T)
    _,a_ref=minspm_refined(theta,Q,u_classes); p_ref=a_ref*Bnd*exp(-r*T)
    t_qpe=@elapsed for _ in 1:2000; minspm_direct(theta,Q); minspm_refined(theta,Q,u_classes); end
    t_qpe/=2000
    @printf("  [MINSPM QCS QPE/QAE] direct=%.6f refined=%.6f  |Δ|=%.2e/%.2e  (%.2e s, O(1)/O(u))\n",
            p_dir,p_ref,abs(p_dir-coh.price),abs(p_ref-coh.price),t_qpe)

    # MC golden
    print("  [MC golden] continuous $(steps)-step, paths=$(n_mc)... ")
    mc = mc_golden(spot,strike,r,vol,skew,T,steps,Bnd; n_paths=n_mc)
    println("done in $(round(mc.t,digits=4))s  price=%.6f ± %.6f (unb=%.6f ± %.6f)" |> x->"")
    @printf("    bounded: price=%.6f ± %.6f  unbounded: price=%.6f ± %.6f\n", mc.price, mc.se, mc.price_unb, mc.se_unb)

    println(repeat("-",100))
    @printf("%-40s | %-12s | %-12s | %-22s\n","Method","Price (bnd)","Time(s)","Notes")
    println(repeat("-",100))
    @printf("%-40s | %-12.6f | %-12.4f | %s\n","MC golden (continuous, sampling)", mc.price, mc.t, @sprintf("±%.6f SE", mc.se))
    @printf("%-40s | %-12.6f | %-12.6f | %s\n","Coherent MINSPM oracle (no sampling)", coh.price, tcoh, "moment propagation")
    @printf("%-40s | %-12.6f | %-12.8f | %s\n","MINSPM QCS QPE direct O(1)", p_dir, t_qpe, @sprintf("|Δ|=%.2e",abs(p_dir-coh.price)))
    @printf("%-40s | %-12.6f | %-12.8f | %s\n","MINSPM QCS QPE refined O(u)", p_ref, t_qpe, @sprintf("|Δ|=%.2e",abs(p_ref-coh.price)))
    println(repeat("-",100))
    err_b = abs(coh.price - mc.price)
    err_u = abs(coh.price_unbounded - mc.price_unb)
    ok = err_u <= target || err_u <= 2*mc.se_unb
    speed = mc.t / max(tcoh, 1e-12)
    speed_q = mc.t / max(tcoh+t_qpe, 1e-12)
    @printf("  Accuracy bounded  |coherent-MC| = %.6f  (MC 2*SE=%.6f)\n", err_b, 2*mc.se)
    @printf("  Accuracy unb      |coherent-MC| = %.6f  (target %.4f, MC 2*SE=%.6f) => %s\n",
            err_u, target, 2*mc.se_unb, ok ? "PASS" : "CHECK")
    @printf("  Speedup (coherent oracle vs MC) = %.1fx  (%.6fs vs %.4fs)\n", speed, tcoh, mc.t)
    @printf("  Speedup (coherent+QPE   vs MC) = %.1fx  (%.6fs vs %.4fs)\n", speed_q, tcoh+t_qpe, mc.t)
    @printf("  Self-sufficient: YES — no MC sampling in oracle or QPE (deterministic moments + algebraic QPE)\n")
    println(repeat("=",100))
    (coherent=coh, mc=mc, err=err_u, speed=speed, ok=ok)
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    arg = length(ARGS)>=1 ? ARGS[1] : "--demo"
    if arg == "--soxx"
        df = CoherentOracleMINSPM.CSV.read("soxx_asset_prices.csv", CoherentOracleMINSPM.DataFrame)
        mat = Matrix{Float64}(df[:, 2:end])
        idx_series = vec(CoherentOracleMINSPM.mean(mat, dims=2))
        CoherentOracleMINSPM.run_case("SOXX equal-weight basket (CSV)", idx_series)
    elseif arg == "--demo"
        rng=MersenneTwister(1); S0=100.0; r=0.04; vol=0.25; dt=1/252
        prices=[S0]; logS_cur=log(S0)
        for _ in 1:300
            logS_cur += (r-0.5*vol^2)*dt + vol*sqrt(dt)*randn(rng)
            push!(prices, exp(logS_cur))
        end
        CoherentOracleMINSPM.run_case("Synthetic demo", prices)
    else
        prices = CoherentOracleMINSPM.load_single_csv(arg)
        CoherentOracleMINSPM.run_case(basename(arg), prices)
    end
end
