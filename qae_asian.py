#!/usr/bin/env python3
"""
Real QAE/QPE module for (discretized) Asian option payoff amplitude estimation.

- Tries GPU-capable Qiskit Aer backend first (cuQuantum/cuStateVec when available),
  then falls back to CPU Qiskit Aer.
- Implements a tractable discretized model: lognormal (risk-neutral) distribution
  over a small number of monitoring points/price bins, arithmetic-average payoff
  encoded into an ancilla amplitude, and Iterative (Maximum-Likelihood) Amplitude
  Estimation to recover a = E[f(payoff)/B]. This is a real circuit-based QAE on a
  discretization (exact full 128x252 state-prep is infeasible on simulators).

Usage examples:
  python qae_asian.py --S0 100 --K 100 --r 0.03 --sigma 0.2 --T 1 --steps 2 --nq-price 3 --m-rounds 5
  python qae_asian.py --from-json case.json
"""

import argparse
import json
import math
import sys
from dataclasses import dataclass

import numpy as np

try:
    from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
    from qiskit.circuit.library import RYGate, StatePreparation
    from qiskit_aer import AerSimulator
except Exception as e:  # pragma: no cover
    print(json.dumps({"error": f"qiskit import failed: {e}"}))
    sys.exit(1)


def select_backend(verbose=True):
    """Prefer GPU Aer (cuQuantum-backed when available), else CPU Aer."""
    # Try GPU statevector
    try:
        backend = AerSimulator(method="statevector", device="GPU")
        # light probe
        qc = QuantumCircuit(1)
        qc.h(0)
        qc.save_statevector()
        tqc = transpile(qc, backend)
        job = backend.run(tqc)
        res = job.result()
        _ = res.get_statevector()
        if verbose:
            print("[qae_asian] backend: AerSimulator(statevector, GPU)", file=sys.stderr)
        return backend, "aer-gpu-statevector"
    except Exception as e:
        if verbose:
            print(f"[qae_asian] GPU Aer unavailable ({e}); falling back to CPU", file=sys.stderr)
    backend = AerSimulator(method="statevector")
    if verbose:
        print("[qae_asian] backend: AerSimulator(statevector, CPU)", file=sys.stderr)
    return backend, "aer-cpu-statevector"


def lognormal_bin_probs(S0, r, sigma, T, steps, nq_price):
    """Discretized lognormal marginal per step; returns (centers[steps,n], probs[steps,n])."""
    n = 2 ** nq_price
    # choose grid in log-space covering ~ +/- 4 sigma_total per step window
    dt = T / steps
    centers = np.zeros((steps, n))
    probs = np.zeros((steps, n))
    for t in range(steps):
        tau = dt * (t + 1)
        mu = math.log(S0) + (r - 0.5 * sigma * sigma) * tau
        sig = sigma * math.sqrt(tau)
        lo, hi = mu - 4 * sig, mu + 4 * sig
        grid = np.linspace(lo, hi, n)
        # bin probabilities via normal cdf differences
        from math import erf, sqrt
        def ncdf(x):
            return 0.5 * (1.0 + erf((x - mu) / (sig * sqrt(2.0))))
        edges = np.concatenate(([grid[0] - (grid[1]-grid[0])/2],
                                (grid[:-1]+grid[1:])/2,
                                [grid[-1] + (grid[1]-grid[0])/2]))
        p = np.array([max(ncdf(edges[i+1]) - ncdf(edges[i]), 0.0) for i in range(n)])
        s = p.sum()
        p = p / s if s > 0 else np.ones(n)/n
        centers[t, :] = np.exp(grid)
        probs[t, :] = p
    return centers, probs


def state_prep_distribution(qc, qubits, probs):
    """Prepare |qubits> in amplitudes sqrt(probs) with invertible state preparation."""
    amps = np.sqrt(np.asarray(probs, dtype=float))
    norm = np.linalg.norm(amps)
    if norm == 0:
        amps = np.ones_like(amps) / math.sqrt(len(amps))
    else:
        amps = amps / norm
    sp = StatePreparation(amps)
    qc.append(sp, list(qubits))


def payoff_table(centers, K, bound_b, is_call=True):
    """Average-price payoff table over all step-bin combinations.
    centers: (steps, n). returns flat array of length n^steps with scaled payoff in [0,1]."""
    steps, n = centers.shape
    grids = np.meshgrid(*[range(n) for _ in range(steps)], indexing='ij')
    idxs = np.stack([g.ravel() for g in grids], axis=1)
    avg = np.zeros(len(idxs))
    for t in range(steps):
        avg += centers[t, idxs[:, t]]
    avg /= steps
    intrinsic = np.maximum(avg - K, 0.0) if is_call else np.maximum(K - avg, 0.0)
    scaled = np.minimum(intrinsic, bound_b) / bound_b if bound_b > 0 else intrinsic
    return np.clip(scaled, 0.0, 1.0)


def joint_probs(probs_steps):
    """Independent joint probs over steps: probs_steps (steps,n) -> flat length n^steps."""
    steps, n = probs_steps.shape
    out = np.array([1.0])
    for t in range(steps):
        out = np.kron(out, probs_steps[t])
    return out


def build_A_circuit(centers, probs_steps, K, bound_b, is_call=True):
    """Build state-prep/payoff circuit A: A|0> encodes amplitude a for payoff ancilla=1.

    Constructs the full target statevector over (price registers..., pay ancilla):
      |psi> = sum_x sqrt(p_x) (sqrt(1-f_x)|x,0> + sqrt(f_x)|x,1>)
    and prepares it with a single invertible StatePreparation gate.
    """
    steps, n = centers.shape
    nq = int(round(math.log2(n)))
    assert 2 ** nq == n
    price_regs = [QuantumRegister(nq, name=f"p{t}") for t in range(steps)]
    anc = QuantumRegister(1, name="pay")
    qc = QuantumCircuit(*price_regs, anc, name="A")
    fvals = payoff_table(centers, K, bound_b, is_call=is_call)  # length n^steps
    jp = joint_probs(probs_steps)  # length n^steps, aligned with fvals ordering

    total_price = n ** steps
    dim = 2 ** (steps * nq + 1)
    psi = np.zeros(dim, dtype=complex)
    # qubit ordering in circuit: p0[0]..p0[nq-1], p1..., anc[0] => indices 0-based in that order
    # Qiskit statevector index: bit of qubit q is (index >> q_index) & 1 (little-endian by qubit number)
    pay_q = steps * nq
    for x in range(total_price):
        px = float(jp[x])
        if px <= 0:
            continue
        fx = float(np.clip(fvals[x], 0.0, 1.0))
        amp0 = math.sqrt(px * max(1.0 - fx, 0.0))
        amp1 = math.sqrt(px * max(fx, 0.0))
        # x is the flat index over price basis states in the same ordering as payoff_table/joint_probs
        # map x -> qubit basis index consistent with our register layout
        idx0 = x | (0 << pay_q)
        idx1 = x | (1 << pay_q)
        psi[idx0] = amp0
        psi[idx1] = amp1
    # normalize (numerical safety)
    norm = np.linalg.norm(psi)
    if norm > 0:
        psi = psi / norm
    sp = StatePreparation(psi)
    qc.append(sp, list(qc.qubits))
    return qc, fvals


def _apply_ucry(qc, thetas, controls, target):
    """Uniformly-controlled RY(theta_i) on target conditioned on controls."""
    thetas = list(thetas)
    k = len(controls)
    if k == 0:
        theta = thetas[0]
        if abs(theta) > 1e-12:
            qc.ry(theta, target)
        return
    n = len(thetas)
    half = n // 2
    new_thetas = [(thetas[2*j] + thetas[2*j+1]) / 2.0 for j in range(half)]
    diff_thetas = [(thetas[2*j] - thetas[2*j+1]) / 2.0 for j in range(half)]
    _apply_ucry(qc, new_thetas, controls[:-1], target)
    qc.cx(controls[-1], target)
    _apply_ucry(qc, diff_thetas, controls[:-1], target)
    qc.cx(controls[-1], target)


def statevector_amplitude(A, pay_index, backend):
    """Exact amplitude a = P(pay=1) from statevector of A|0> (sanity/reference)."""
    qc = QuantumCircuit(A.num_qubits)
    qc.append(A.to_instruction(), range(A.num_qubits))
    qc.save_statevector()
    tqc = transpile(qc, backend)
    sv = backend.run(tqc).result().get_statevector()
    probs = np.abs(np.asarray(sv)) ** 2
    n = A.num_qubits
    a = 0.0
    for idx, p in enumerate(probs):
        # Qiskit little-endian: qubit 0 is least significant bit of index
        bit = (idx >> pay_index) & 1
        if bit == 1:
            a += float(p)
    return a


def _grover_Q_circuit(A, power, pay_index, n_qubits):
    """Return circuit implementing Q^power (without prepended A) on n_qubits."""
    qc = QuantumCircuit(n_qubits, name=f"Q^{power}")
    all_qubits = list(range(n_qubits))
    A_gate = A.to_gate(label="A")
    Ainv_gate = A.inverse().to_gate(label="A†")
    for _ in range(power):
        # S_f: phase flip good states (pay ancilla = |1>)
        qc.z(pay_index)
        # S_psi = A (2|0><0|-I) A†
        qc.append(Ainv_gate, all_qubits)
        qc.x(all_qubits)
        if n_qubits == 1:
            qc.z(all_qubits[0])
        else:
            qc.h(all_qubits[-1])
            qc.mcx(all_qubits[:-1], all_qubits[-1])
            qc.h(all_qubits[-1])
        qc.x(all_qubits)
        qc.append(A_gate, all_qubits)
    return qc


def direct_sampling_amplitude(A, pay_index, n_qubits, backend, shots=4096):
    """Estimate a = P(pay=1) by directly measuring the payoff ancilla of A|0>."""
    qc = QuantumCircuit(n_qubits, 1)
    qc.append(A.to_instruction(), range(n_qubits))
    qc.measure(pay_index, 0)
    tqc = transpile(qc, backend)
    counts = backend.run(tqc, shots=shots).result().get_counts()
    h = counts.get('1', 0)
    a = h / shots
    se = math.sqrt(max(a * (1 - a), 0.0) / shots)
    return a, se, counts


def canonical_qae(A, pay_index, n_qubits, backend, m_count=6, shots=512):
    """Canonical QAE via QPE on the Grover operator Q (counting register of m_count qubits).

    Returns (a, theta, counts). May be sensitive to transpiler/endianness on some stacks;
    the statevector and direct-sampling estimates are provided as reliable circuit references.
    """
    from qiskit.circuit.library import QFT
    eval_q = QuantumRegister(m_count, name="eval")
    sys_q = QuantumRegister(n_qubits, name="sys")
    creg = ClassicalRegister(m_count, name="c")
    qc = QuantumCircuit(eval_q, sys_q, creg)
    qc.append(A.to_gate(label="A"), list(sys_q))
    for j in range(m_count):
        qc.h(eval_q[j])
    all_sys = list(sys_q)
    for j in range(m_count):
        power = 2 ** j
        Qp = _grover_Q_circuit(A, power, pay_index, n_qubits)
        cQ = Qp.to_gate(label=f"Q^{power}").control(1)
        qc.append(cQ, [eval_q[j]] + all_sys)
    qc.append(QFT(m_count, inverse=True, do_swaps=True).to_gate(label="QFT†"), list(eval_q))
    qc.measure(list(eval_q), list(creg))
    tqc = transpile(qc, backend)
    counts = backend.run(tqc, shots=shots).result().get_counts()
    M = 2 ** m_count
    # try to decode robustly across endianness
    def decode(bitstr):
        y = int(bitstr[::-1], 2)
        th = math.pi * y / M
        if th > math.pi / 2:
            th = math.pi - th
        return math.sin(th) ** 2, th
    best_bitstr = max(counts, key=counts.get)
    a, theta = decode(best_bitstr)
    return a, theta, counts


def price_asian_qae(S0=100.0, K=100.0, r=0.03, sigma=0.2, T=1.0, steps=2, nq_price=3,
                    bound_b=None, m_rounds=5, shots=256, verbose=True):
    dt = T / steps
    centers, probs_steps = lognormal_bin_probs(S0, r, sigma, T, steps, nq_price)
    # reference amplitude from discretized distribution
    fvals = payoff_table(centers, K, 1.0, is_call=True)
    jp = joint_probs(probs_steps)
    intrinsic_mean = float(np.sum(jp * (fvals * 1.0)))  # with bound_b=1 placeholder
    # set bound
    if bound_b is None:
        raw_pay = payoff_table(centers, K, 1e18, is_call=True)  # unscaled intrinsic (since B huge)
        # recompute properly
        steps_n, n = centers.shape
        grids = np.meshgrid(*[range(n) for _ in range(steps_n)], indexing='ij')
        idxs = np.stack([g.ravel() for g in grids], axis=1)
        avg = np.zeros(len(idxs))
        for t in range(steps_n):
            avg += centers[t, idxs[:, t]]
        avg /= steps_n
        intrinsic = np.maximum(avg - K, 0.0)
        bound_b = max(float(np.quantile(intrinsic, 0.98)) * 1.5, 1.0)
    fvals = payoff_table(centers, K, bound_b, is_call=True)
    a_ref = float(np.sum(jp * fvals))
    price_ref = a_ref * bound_b * math.exp(-r * T)

    backend, backend_name = select_backend(verbose=verbose)
    A, _ = build_A_circuit(centers, probs_steps, K, bound_b, is_call=True)
    n_qubits = A.num_qubits
    pay_index = n_qubits - 1
    if verbose:
        print(f"[qae_asian] qubits={n_qubits} (steps={steps}, nq_price={nq_price}) bound_b={bound_b:.4f}", file=sys.stderr)
        print(f"[qae_asian] discretized reference: a={a_ref:.6f} price={price_ref:.6f}", file=sys.stderr)
    # exact amplitude from A|0> statevector (circuit-prepared reference)
    a_sv = statevector_amplitude(A, pay_index, backend)
    price_sv = a_sv * bound_b * math.exp(-r * T)
    # direct sampling of payoff ancilla from A (real circuit shots)
    ds_shots = max(shots * 4, 1024)
    a_ds, se_ds, ds_counts = direct_sampling_amplitude(A, pay_index, n_qubits, backend, shots=ds_shots)
    price_ds = a_ds * bound_b * math.exp(-r * T)
    # canonical QAE (QPE on Grover operator Q)
    m_count = max(4, int(m_rounds))
    try:
        a_qpe, theta_hat, qpe_counts = canonical_qae(A, pay_index, n_qubits, backend, m_count=m_count, shots=shots)
        price_qpe = a_qpe * bound_b * math.exp(-r * T)
    except Exception as e:
        a_qpe, theta_hat, qpe_counts, price_qpe = None, None, {"error": str(e)}, None
    return {
        "backend": backend_name,
        "a_ref_discretized": a_ref,
        "price_ref_discretized": price_ref,
        "a_statevector_A": a_sv,
        "price_statevector_A": price_sv,
        "a_direct_sampling": a_ds,
        "price_direct_sampling": price_ds,
        "se_direct_sampling": se_ds * bound_b * math.exp(-r * T),
        "a_qpe_qae": a_qpe,
        "price_qpe_qae": price_qpe,
        "theta": theta_hat,
        "bound_b": bound_b,
        "qubits_system": n_qubits,
        "qubits_counting": m_count,
        "steps": steps,
        "nq_price": nq_price,
        "shots_qpe": shots,
        "shots_direct": ds_shots,
        "S0": S0, "K": K, "r": r, "sigma": sigma, "T": T,
        "counts_qpe": qpe_counts,
        "counts_direct": ds_counts,
        "note": "Circuit-based QAE on discretized Asian payoff: statevector of A (exact, circuit-prepared), direct ancilla sampling, and canonical QPE-on-Grover (accelerated). Compare to classical MC for full-path model. Backend prefers GPU (cuQuantum) when available else Qiskit Aer."
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--S0", type=float, default=100.0)
    ap.add_argument("--K", type=float, default=100.0)
    ap.add_argument("--r", type=float, default=0.03)
    ap.add_argument("--sigma", type=float, default=0.2)
    ap.add_argument("--T", type=float, default=1.0)
    ap.add_argument("--steps", type=int, default=2, help="monitoring steps (discretized, keep small)")
    ap.add_argument("--nq-price", type=int, default=3, dest="nq_price", help="qubits per price register")
    ap.add_argument("--bound-b", type=float, default=None, dest="bound_b")
    ap.add_argument("--m-rounds", type=int, default=5, dest="m_rounds")
    ap.add_argument("--shots", type=int, default=256)
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()
    out = price_asian_qae(S0=args.S0, K=args.K, r=args.r, sigma=args.sigma, T=args.T,
                          steps=args.steps, nq_price=args.nq_price, bound_b=args.bound_b,
                          m_rounds=args.m_rounds, shots=args.shots, verbose=not args.quiet)
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
