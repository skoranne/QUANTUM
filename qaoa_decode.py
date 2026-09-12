"""
Recover the Turan max-clique solution from QAOA simulation output.

Pipeline:
  1. Build the parameterized QAOA circuit (gamma, beta).
  2. Run it on a simulator with many shots -> counts histogram.
  3. Decode each bitstring into a vertex set, fixing Qiskit's
     c[n-1]...c[0] bit ordering.
  4. Score every sample with the QUBO cost AND independently
     re-verify against the actual edge list.
  5. Sweep (gamma, beta) classically (stand-in for a real
     COBYLA/SPSA optimizer loop) to find good parameters.
  6. Report every bitstring tied for the best observed cost.
"""

from itertools import combinations
from qiskit import QuantumCircuit
from qiskit_aer import AerSimulator

# ---- graph definition: perturbed T(8,3) ----
n = 8
vertices = list(range(1, n + 1))
parts = [{1, 2, 3}, {4, 5, 6}, {7, 8}]

edges = set()
for i in vertices:
    for j in vertices:
        if i < j and not any(i in p and j in p for p in parts):
            edges.add((i, j))          # base cross-part edges

edges.add((1, 2))                      # perturbation: add
edges.add((4, 5))                      # perturbation: add
edges.discard((1, 4))                  # perturbation: remove


def is_clique(vertex_set):
    """Ground-truth check -- always verify this way, independent
    of whatever the QUBO/Ising algebra claims."""
    return all((i, j) in edges for i, j in combinations(sorted(vertex_set), 2))


def cost(vertex_set):
    """QUBO cost: -|clique size| + 2*(# violating non-edge pairs)."""
    violations = sum(
        1 for i, j in combinations(sorted(vertex_set), 2) if (i, j) not in edges
    )
    return -len(vertex_set) + 2 * violations


# ---- build the QAOA circuit for given (gamma, beta) ----
def build_circuit(gamma, beta):
    qc = QuantumCircuit(n, n)
    qc.h(range(n))

    for v in (1, 3, 4, 6):             # linear Z terms
        qc.rz(-gamma, v - 1)

    non_edges = [(1, 3), (2, 3), (4, 6), (5, 6), (7, 8), (1, 4)]
    for i, j in non_edges:             # ZZ terms
        a, b = i - 1, j - 1
        qc.cx(a, b)
        qc.rz(gamma, b)
        qc.cx(a, b)

    for q in range(n):                 # mixer
        qc.rx(2 * beta, q)

    qc.measure(range(n), range(n))
    return qc


# ---- decode a Qiskit bitstring, correcting bit order ----
def decode(bitstring):
    bits = bitstring[::-1]             # bits[k] == c[k]
    return {v for v in vertices if bits[v - 1] == "1"}


backend = AerSimulator()


def run_and_score(gamma, beta, shots=4096):
    qc = build_circuit(gamma, beta)
    counts = backend.run(qc, shots=shots).result().get_counts()
    total = sum(counts.values())

    exp_val = 0.0
    best_cost, best_sets = None, []
    for bitstring, count in counts.items():
        vset = decode(bitstring)
        c = cost(vset)
        exp_val += (count / total) * c
        if best_cost is None or c < best_cost:
            best_cost, best_sets = c, [(vset, count)]
        elif c == best_cost:
            best_sets.append((vset, count))

    return exp_val, best_cost, best_sets


# ---- classical outer loop over (gamma, beta) ----
# (a real implementation would use scipy.optimize / COBYLA / SPSA;
#  grid search here just to make the mechanics visible)
best_overall = None
for gamma in (0.2, 0.4, 0.6, 0.8):
    for beta in (0.1, 0.3, 0.5, 0.7):
        exp_val, best_cost, best_sets = run_and_score(gamma, beta)
        if best_overall is None or exp_val < best_overall[0]:
            best_overall = (exp_val, gamma, beta, best_cost, best_sets)

exp_val, gamma_star, beta_star, best_cost, best_sets = best_overall

print(f"Best (gamma, beta) found: ({gamma_star}, {beta_star})")
print(f"<H_C> at these parameters: {exp_val:.3f}")
print(f"Best observed cost: {best_cost}")
print("Recovered candidate(s):")
for vset, count in best_sets:
    print(f"  vertices={sorted(vset)}  shots={count}  "
          f"is_clique(independently verified)={is_clique(vset)}")
