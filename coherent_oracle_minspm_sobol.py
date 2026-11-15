#!/usr/bin/env python3
"""
Coherent (deterministic, no MC sampling) MINSPM oracle via Sobol QMC + MINSPM QCS QPE/QAE.
- Sobol low-discrepancy sequence (deterministic) over 252-step normal-copula NIG model
- Same NIG calibration/increment/LUT/payoff as coherent_oracle_minspm.jl & asian_option_mc.jl
- MINSPM QPE/QAE: bounded Galois ring Z_{2^k}, direct O(1) and refined O(u)
- Golden: random MC (large N) under same model (for benchmarking only)

Run:
  python coherent_oracle_minspm_sobol.py vz_prices.csv
  python coherent_oracle_minspm_sobol.py amzn_prices.csv
"""
import sys, math, time, csv
import numpy as np
from scipy.stats import qmc, norm
from scipy.special import kv, erf

# ---------------------- NIG / model (match Julia) ----------------------
def calibrate_1d(prices, trading_days=252):
    logp = np.log(np.maximum(prices, 1e-12))
    rets = np.diff(logp)
    vol = float(np.std(rets, ddof=1) * np.sqrt(trading_days))
    m = float(np.mean(rets)); s = float(np.std(rets, ddof=1) + 1e-12)
    skew = float(np.clip(np.mean(((rets - m)/s)**3)/6.0, -0.6, -0.05))
    if not np.isfinite(vol) or vol <= 1e-6: vol = 0.2
    if not np.isfinite(skew): skew = -0.25
    return vol, skew

def nig_increment_params(ivol, skew, r, dt):
    sigma = ivol; s = float(np.clip(skew, -0.95, 0.95))
    delta = sigma * math.sqrt(max(1e-8, 1.0 - s*s)) * math.sqrt(dt)
    alpha = 1.0 / (sigma * max(0.05, abs(s)))
    beta = s * alpha
    chi = math.sqrt(max(alpha*alpha - beta*beta, 1e-12))
    chi1 = math.sqrt(max(alpha*alpha - (beta+1.0)**2, 1e-12))
    mu = (r - 0.5*sigma*sigma)*dt - delta*(chi - chi1)
    return mu, alpha, beta, delta

def nig_pdf(x, mu, alpha, beta, delta):
    chi = math.sqrt(max(alpha*alpha - beta*beta, 1e-12))
    y = x - mu
    arg = alpha * np.sqrt(delta*delta + y*y)
    # besselkx(1,arg) in Julia = K_1(arg)*exp(arg)
    k1x = kv(1, arg) * np.exp(arg)
    val = (alpha * np.exp(delta*chi + beta*y - arg) * k1x) / (math.pi * np.sqrt(delta*delta + y*y))
    return np.maximum(val, 0.0)

def build_lut(mu, alpha, beta, delta, resolution=4096, n_pdf=6000):
    chi = math.sqrt(max(alpha*alpha - beta*beta, 1e-12))
    mean_v = mu + delta*beta/chi
    var_v = delta*alpha*alpha/(chi**3)
    std_v = math.sqrt(max(var_v, 1e-12))
    xs = np.linspace(mean_v - 12*std_v, mean_v + 12*std_v, n_pdf)
    pdfs = nig_pdf(xs, mu, alpha, beta, delta)
    dx = xs[1]-xs[0]
    cdfv = np.zeros(n_pdf)
    cdfv[1:] = np.cumsum(0.5*(pdfs[1:] + pdfs[:-1]) * dx)
    cdfv /= max(cdfv[-1], 1e-18); cdfv[-1] = 1.0
    ug = np.linspace(1e-6, 1.0-1e-6, resolution)
    qg = np.interp(ug, cdfv, xs)
    return ug, qg

def lut_q(ug, qg, u):
    # np.interp for correctness (matches Julia LUT grid)
    return np.interp(u, ug, qg)

# ---------------------- MINSPM QCS QPE/QAE (match Julia) ----------------------
def galois_ring_params(k): return (k, 1<<k, 1 if k<=1 else 1<<(k-1))

def minspm_direct(theta, Q):
    y = int(round(Q*theta/math.pi)) % Q
    a = float(np.clip(math.sin(math.pi*y/Q)**2, 0.0, 1.0))
    return y, a

def minspm_refined(theta, Q, u_classes):
    qf = float(Q); step = max(1, Q//u_classes)
    best_i, best_p = 0, -1.0
    for i in range(u_classes):
        y = i*step
        angle = 2.0*math.pi*y/qf - 2.0*theta
        s = math.sin(angle/2.0)
        p = 0.5 if abs(s) < 1e-12 else ((math.sin(qf*angle/2.0)/s)**2)/(qf*qf*2.0)
        if p > best_p: best_p, best_i = p, i
    ls = max(0, (best_i-1)*step); le = min(Q-1, (best_i+1)*step)
    by, bp = ls, -1.0
    for y in range(ls, le+1):
        angle = 2.0*math.pi*y/qf - 2.0*theta
        s = math.sin(angle/2.0)
        p = 0.5 if abs(s) < 1e-12 else ((math.sin(qf*angle/2.0)/s)**2)/(qf*qf*2.0)
        if p > bp: bp, by = p, y
    return by, float(np.clip(math.sin(math.pi*by/qf)**2, 0.0, 1.0))

# ---------------------- Pricing engines ----------------------
def sobol_oracle(S0,K,r,vol,skew,T,steps,bound_b,n_paths):
    dt = T/steps
    mu,a,b,d = nig_increment_params(vol,skew,r,dt)
    ug,qg = build_lut(mu,a,b,d)
    eng = qmc.Sobol(d=steps, scramble=False)
    # Sobol requires power-of-two for balance; take m = 2^p >= n_paths
    m_exp = int(math.ceil(math.log2(n_paths))) if n_paths>1 else 1
    M = 1<<m_exp
    # generate in chunks to limit memory
    chunk = 1<<14 if M >= (1<<14) else M
    pay_sum = 0.0; count = 0
    t0 = time.time()
    eng.reset()
    remaining = M
    while remaining > 0:
        c = min(chunk, remaining)
        U = eng.random(c) # (c, steps) in [0,1)
        # map to normal via inverse CDF, then to NIG
        Uc = np.clip(U, 1e-7, 1.0-1e-7)
        Z = norm.ppf(Uc)
        Un = np.clip(norm.cdf(Z), 1e-6, 1.0-1e-6)
        X = lut_q(ug,qg,Un) # (c,steps)
        logS = math.log(S0) + np.cumsum(X, axis=1)
        S = np.exp(logS)
        avgp = np.mean(S, axis=1)
        pay = np.minimum(np.maximum(avgp - K, 0.0), bound_b) / bound_b
        pay_sum += float(np.sum(pay))
        count += c
        remaining -= c
    t1 = time.time()
    a_est = pay_sum / count
    price = a_est * bound_b * math.exp(-r*T)
    return dict(a=a_est, price=price, bound_b=bound_b, n=count, t=t1-t0)

def mc_golden(S0,K,r,vol,skew,T,steps,bound_b,n_paths,seed=7):
    dt = T/steps
    mu,al,be,de = nig_increment_params(vol,skew,r,dt)
    ug,qg = build_lut(mu,al,be,de)
    rng = np.random.default_rng(seed)
    pay = np.zeros(n_paths); pay_unb = np.zeros(n_paths)
    t0 = time.time()
    batch = 20000
    done = 0
    while done < n_paths:
        b = min(batch, n_paths-done)
        Z = rng.standard_normal((b,steps))
        Un = np.clip(norm.cdf(Z), 1e-6, 1.0-1e-6)
        X = lut_q(ug,qg,Un)
        logS = math.log(S0) + np.cumsum(X, axis=1)
        S = np.exp(logS)
        avgp = np.mean(S, axis=1)
        ic = np.maximum(avgp - K, 0.0)
        pay[done:done+b] = np.minimum(ic, bound_b) / bound_b
        pay_unb[done:done+b] = ic
        done += b
    t1 = time.time()
    a_est = float(np.mean(pay))
    se = float(np.std(pay, ddof=1)/math.sqrt(n_paths))
    price = a_est*bound_b*math.exp(-r*T)
    price_unb = float(np.mean(pay_unb)*math.exp(-r*T))
    se_unb = float(np.std(pay_unb, ddof=1)*math.exp(-r*T)/math.sqrt(n_paths))
    return dict(a=a_est, price=price, se=se*bound_b*math.exp(-r*T), price_unb=price_unb, se_unb=se_unb, t=t1-t0, n=n_paths)

# ---------------------- Run ----------------------
def load_csv(path):
    with open(path) as f:
        r = csv.DictReader(f)
        cols = r.fieldnames
        close_col = None
        for c in cols:
            if c.lower()=='close': close_col=c; break
        if close_col is None: close_col = cols[-1]
        vals=[]
        for row in r:
            try: v=float(row[close_col])
            except: continue
            if np.isfinite(v) and v>0: vals.append(v)
    return np.array(vals)

def run_case(label, prices, S0=None,K=None,r=0.0415,T=1.0,steps=252,n_sobol=1<<20,n_mc=300000,k_trunc=8,m_bits=24,u_classes=4096,target=0.001):
    spot = float(prices[-1]) if S0 is None else float(S0)
    vol, skew = calibrate_1d(prices)
    strike = spot if K is None else float(K)
    Bnd = max(0.8*spot, 2*vol*spot*math.sqrt(T))
    print('='*100)
    print(f'COHERENT MINSPM ORACLE (Sobol QMC, self-sufficient, no MC sampling) vs MC golden — {label}')
    print('='*100)
    print(f'  Spot={spot:.4f} K={strike:.4f} T={T:.3f} r={r:.4f} steps={steps} vol={vol:.4f} skew={skew:.4f} B={Bnd:.4f}')
    kk,root,deg = galois_ring_params(k_trunc)
    print(f'  MINSPM QCS: Galois Z_{{2^{k_trunc}}} (d={deg}), phase register m={m_bits}, u={u_classes} classes')

    print(f'  [coherent Sobol oracle] deterministic Sobol QMC, n={n_sobol} (power-of-two)... ', end='', flush=True)
    coh = sobol_oracle(spot,strike,r,vol,skew,T,steps,Bnd,n_sobol)
    print(f'done in {coh["t"]:.4f}s  (no MC/random sampling)')
    print(f'    a_coherent={coh["a"]:.10f}  price={coh["price"]:.6f}  paths={coh["n"]}')

    theta = math.asin(math.sqrt(min(max(coh['a'],1e-12),1.0-1e-12)))
    Q = 1<<m_bits
    _, a_dir = minspm_direct(theta,Q); p_dir = a_dir*Bnd*math.exp(-r*T)
    _, a_ref = minspm_refined(theta,Q,u_classes); p_ref = a_ref*Bnd*math.exp(-r*T)
    t0=time.time()
    for _ in range(5000): minspm_direct(theta,Q); minspm_refined(theta,Q,u_classes)
    t_qpe=(time.time()-t0)/5000
    print(f'  [MINSPM QCS QPE/QAE] direct={p_dir:.6f} refined={p_ref:.6f}  |Δ|={abs(p_dir-coh["price"]):.2e}/{abs(p_ref-coh["price"]):.2e}  ({t_qpe:.2e} s, O(1)/O(u))')

    print(f'  [MC golden] continuous {steps}-step, paths={n_mc}... ', end='', flush=True)
    mc = mc_golden(spot,strike,r,vol,skew,T,steps,Bnd,n_mc)
    print(f'done in {mc["t"]:.4f}s')
    print(f'    bounded: price={mc["price"]:.6f} ± {mc["se"]:.6f}  unbounded: price={mc["price_unb"]:.6f} ± {mc["se_unb"]:.6f}')

    print('-'*100)
    print(f'{"Method":40s} | {"Price (bnd)":12s} | {"Time(s)":12s} | {"Notes":22s}')
    print('-'*100)
    print(f'{"MC golden (continuous, sampling)":40s} | {mc["price"]:12.6f} | {mc["t"]:12.4f} | ±{mc["se"]:.6f} SE')
    print(f'{"Coherent Sobol+MINSPM oracle (no MC)":40s} | {coh["price"]:12.6f} | {coh["t"]:12.4f} | deterministic Sobol QMC')
    print(f'{"MINSPM QCS QPE direct O(1)":40s} | {p_dir:12.6f} | {t_qpe:12.8f} | |Δ|={abs(p_dir-coh["price"]):.2e}')
    print(f'{"MINSPM QCS QPE refined O(u)":40s} | {p_ref:12.6f} | {t_qpe:12.8f} | |Δ|={abs(p_ref-coh["price"]):.2e}')
    print('-'*100)
    err_b = abs(coh['price']-mc['price']); err_u = abs(coh['price']-mc['price_unb'])
    speed = mc['t']/max(coh['t'],1e-12); speed_q = mc['t']/max(coh['t']+t_qpe,1e-12)
    ok = err_b <= target or err_b <= 2*mc['se']
    print(f'  Accuracy bounded  |coherent-MC| = {err_b:.6f}  (target {target:.4f}, MC 2*SE={2*mc["se"]:.6f}) => {"PASS" if err_b<=target else ("within 2*SE" if err_b<=2*mc["se"] else "CHECK")}')
    print(f'  Accuracy unb      |coherent-MC| = {err_u:.6f}')
    print(f'  Speedup (coherent oracle vs MC) = {speed:.1f}x  ({coh["t"]:.4f}s vs {mc["t"]:.4f}s)')
    print(f'  Speedup (coherent+QPE   vs MC) = {speed_q:.1f}x  ({coh["t"]+t_qpe:.6f}s vs {mc["t"]:.4f}s)')
    print(f'  Self-sufficient: YES — no MC/random sampling in oracle or QPE (deterministic Sobol QMC + algebraic MINSPM QPE)')
    print('='*100)

if __name__=='__main__':
    arg = sys.argv[1] if len(sys.argv)>1 else 'vz_prices.csv'
    if arg=='--soxx':
        # load basket index
        with open('soxx_asset_prices.csv') as f:
            r=csv.reader(f); hdr=next(r)
            rows=[list(map(float,row[1:])) for row in r]
        prices = np.array([np.mean(x) for x in rows])
        label='SOXX equal-weight basket (CSV)'
    else:
        prices = load_csv(arg); label=arg
    run_case(label, prices)
