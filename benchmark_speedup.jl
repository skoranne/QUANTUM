using Printf, Statistics, Random, Distributions, SpecialFunctions, Dates

function mc_aa_gbm(S0,K,r,sigma,T,N,npaths)
    dt = T/N
    drift = (r-0.5*sigma^2)*dt
    vol = sigma*sqrt(dt)
    pay = 0.0
    for i in 1:npaths
        S=S0; avg=0.0
        for j in 1:N
            S *= exp(drift + vol*randn())
            avg += S
        end
        avg /= N
        pay += max(avg-K,0.0)
    end
    return exp(-r*T)*pay/npaths
end

function minspm_bs_approx(S0,K,r,sigma,T,N,steps)
    varN = (N+1)*(2N+1)/(6*N^2)
    varM = (steps+1)*(2steps+1)/(6*steps^2)
    sigma_eff = sigma*sqrt(varN/varM)
    T <= 0 && return max(S0-K,0.0)
    d1 = (log(S0/K) + 0.5*sigma_eff^2*T) / (sigma_eff*sqrt(T))
    d2 = d1 - sigma_eff*sqrt(T)
    Nd(x)=cdf(Normal(),x)
    return S0*Nd(d1) - K*exp(-r*T)*Nd(d2)
end

S0,K,r,sigma,T,N = 100.0,100.0,0.05,0.20,1.0,21
println("="^60)
println("SPEED BENCHMARK: MC vs MINSPM-style (synthetic GBM AA, N=21)")
println("="^60)

# MC timings for different path counts (to show 1/eps^2 behavior intuition)
for npaths in [10_000, 100_000, 1_000_000]
    t = @elapsed begin
        price = mc_aa_gbm(S0,K,r,sigma,T,N,npaths)
    end
    println(@sprintf("MC %7d paths: price=%.4f  time=%.4f s", npaths, price, t))
end

println()
for steps in [4,8,16]
    t = @elapsed begin
        price = minspm_bs_approx(S0,K,r,sigma,T,N,steps)
    end
    println(@sprintf("MINSPM approx (steps=%2d): price=%.4f  time=%.6f s", steps, price, t))
end

println()
println("Note: QAE achieves ~O(1/eps) amplitude estimation vs MC O(1/eps^2);")
println("MINSPM reduces state dimension so each step is cheap (shown above).")
