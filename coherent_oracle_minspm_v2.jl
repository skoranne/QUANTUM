# Coherent MINSPM oracle v2: uses deterministic Sobol QMC (self-sufficient, no MC sampling)
# as high-fidelity amplitude oracle, then MINSPM QCS QPE/QAE for phase recovery.
# Calls Python sobol engine (scipy qmc.Sobol, deterministic) for the oracle.
# Run: julia --project=. coherent_oracle_minspm_v2.jl vz_prices.csv [n_sobol_pow2]
include("coherent_oracle_minspm.jl")
using Printf

function run_sobol_oracle(path, n_sobol)
    # call python to get deterministic Sobol amplitude
    cmd = `python3 -c "
import coherent_oracle_minspm_sobol as s, json, math
prices = s.load_csv('$path')
spot=prices[-1]; vol,skew=s.calibrate_1d(prices); K=spot; r=0.0415; T=1.0; steps=252
Bnd=max(0.8*spot,2*vol*spot*math.sqrt(T))
coh = s.sobol_oracle(spot,K,r,vol,skew,T,steps,Bnd,$n_sobol)
print(coh['a']); print(); print(coh['price']); print(); print(coh['t']); print(); print(vol); print(); print(skew); print(); print(Bnd); print(); print(spot)
"`
    out = read(cmd, String)
    lines = split(strip(out), '\n')
    a = parse(Float64, strip(lines[1]))
    price = parse(Float64, strip(lines[2]))
    t = parse(Float64, strip(lines[3]))
    vol = parse(Float64, strip(lines[4]))
    skew = parse(Float64, strip(lines[5]))
    Bnd = parse(Float64, strip(lines[6]))
    spot = parse(Float64, strip(lines[7]))
    (a=a, price=price, t=t, vol=vol, skew=skew, Bnd=Bnd, spot=spot)
end

if abspath(PROGRAM_FILE) == @__FILE__
    path = length(ARGS)>=1 ? ARGS[1] : "vz_prices.csv"
    n_sobol = length(ARGS)>=2 ? parse(Int,ARGS[2]) : (1<<22)
    println("Running deterministic Sobol+MINSPM oracle for $path with n=$n_sobol")
    coh = run_sobol_oracle(path, n_sobol)
    Bnd=coh.Bnd; spot=coh.spot; r=0.0415; T=1.0
    # MINSPM QPE
    theta=asin(sqrt(clamp(coh.a,1e-12,1.0-1e-12))); Q=1<<24
    _,a_dir=CoherentOracleMINSPM.minspm_direct(theta,Q); p_dir=a_dir*Bnd*exp(-r*T)
    _,a_ref=CoherentOracleMINSPM.minspm_refined(theta,Q,4096); p_ref=a_ref*Bnd*exp(-r*T)
    t_qpe=@elapsed for _ in 1:2000; CoherentOracleMINSPM.minspm_direct(theta,Q); CoherentOracleMINSPM.minspm_refined(theta,Q,4096); end
    t_qpe/=2000
    println("Coherent Sobol price: $(coh.price) (a=$(coh.a), t=$(coh.t)s, deterministic, no MC sampling)")
    println("MINSPM QPE direct=$p_dir refined=$p_ref |Δ|=$(abs(p_dir-coh.price)) t=$t_qpe")
end
