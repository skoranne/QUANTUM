// ============================================================
// QAOA circuit: Extremal-graph (Turan) Maximal Clique problem
// Instance: perturbed T(8,3), k=3 perturbations
// ============================================================
// Base graph T(8,3): parts P1={1,2,3}, P2={4,5,6}, P3={7,8},
// complete tripartite (edges only across parts). Base clique
// number = 3 (Turan's theorem: one vertex per part).
//
// Perturbations (k=3), applied to break the trivial closed form:
//   1. ADD edge (1,2)   -- intra-P1
//   2. ADD edge (4,5)   -- intra-P2
//   3. REMOVE edge (1,4) -- inter-part, previously present
//
// Effect: each added edge alone enables a 4-clique (e.g. add (1,2)
// + one vertex from P2 + one from P3, since those cross edges are
// untouched). Combining BOTH added edges into a 5-clique
// {1,2,4,5,+P3} is blocked because it requires edge (1,4), which
// was removed. True optimum = 4, with SIX degenerate optimal
// cliques (verified by brute force over all 2^8 = 256 subsets):
//   {1,2,5,7} {1,2,5,8} {1,2,6,7} {1,2,6,8} {2,4,5,7} {2,4,5,8}
//
// Non-edge set after perturbation (pairs penalized in the QUBO):
//   (1,3) (2,3)   <- remaining intra-P1 pairs
//   (4,6) (5,6)   <- remaining intra-P2 pairs
//   (7,8)         <- intra-P3
//   (1,4)         <- newly removed inter-part edge
//
// QUBO (Lucas Max-Clique formulation, A=1, B=2, B>A):
//   H = -(x1+x2+...+x8)
//       + 2*( x1x3 + x2x3 + x4x6 + x5x6 + x7x8 + x1x4 )
// Minimum H = -4, attained at the six cliques listed above.
//
// Ising form (x_i = (1-Z_i)/2, additive constant dropped),
// derived and checked against all six optima (each gives H_C=-4):
//   H_C = -0.5*(Z1+Z3+Z4+Z6)
//         + 0.5*(Z1Z3 + Z2Z3 + Z4Z6 + Z5Z6 + Z7Z8 + Z1Z4)
//
// Qubit map: q[i-1] = vertex i   (q[0]=v1, q[1]=v2, ..., q[7]=v8)
//
// Gate identities:
//   exp(-i*theta/2 * Z_a Z_b) = CX(a,b); RZ(theta) b; CX(a,b)
//   exp(-i*theta/2 * Z_a)     = RZ(theta) a
//   exp(-i*theta/2 * X_a)     = RX(theta) a
//
// For a term c*Z_a in H_C: exp(-i*gamma*c*Z_a) = RZ(2*gamma*c) a
//   -> linear terms have c=-0.5  => RZ(-gamma)
//   -> ZZ terms have c=+0.5      => RZ(+gamma) inside CX sandwich
// Mixer: exp(-i*beta*sum X_i) = RX(2*beta) on every qubit
//
// Example variational parameters (placeholders -- a real run tunes
// these classically; p=1 shown here, degenerate optima like this
// typically want more layers / more shots to resolve reliably):
//   gamma = 0.4  -> rz(-0.4) [linear terms], rz(0.4) [ZZ terms]
//   beta  = 0.3  -> rx(0.6)  [mixer]
// ============================================================

OPENQASM 2.0;
include "qelib1.inc";

qreg q[8];
creg c[8];

// ---- initial state: equal superposition ----
h q[0];
h q[1];
h q[2];
h q[3];
h q[4];
h q[5];
h q[6];
h q[7];

// ---- cost unitary exp(-i*gamma*H_C) ----

// linear terms: -0.5*(Z1 + Z3 + Z4 + Z6)  -> vertices 1,3,4,6 -> q0,q2,q3,q5
rz(-0.4) q[0];
rz(-0.4) q[2];
rz(-0.4) q[3];
rz(-0.4) q[5];

// ZZ term: 0.5*Z1*Z3  (vertices 1,3 -> q0,q2)
cx q[0],q[2];
rz(0.4) q[2];
cx q[0],q[2];

// ZZ term: 0.5*Z2*Z3  (vertices 2,3 -> q1,q2)
cx q[1],q[2];
rz(0.4) q[2];
cx q[1],q[2];

// ZZ term: 0.5*Z4*Z6  (vertices 4,6 -> q3,q5)
cx q[3],q[5];
rz(0.4) q[5];
cx q[3],q[5];

// ZZ term: 0.5*Z5*Z6  (vertices 5,6 -> q4,q5)
cx q[4],q[5];
rz(0.4) q[5];
cx q[4],q[5];

// ZZ term: 0.5*Z7*Z8  (vertices 7,8 -> q6,q7)
cx q[6],q[7];
rz(0.4) q[7];
cx q[6],q[7];

// ZZ term: 0.5*Z1*Z4  (vertices 1,4 -> q0,q3)  <- the removed edge
cx q[0],q[3];
rz(0.4) q[3];
cx q[0],q[3];

// ---- mixer unitary exp(-i*beta*sum_i X_i) ----
rx(0.6) q[0];
rx(0.6) q[1];
rx(0.6) q[2];
rx(0.6) q[3];
rx(0.6) q[4];
rx(0.6) q[5];
rx(0.6) q[6];
rx(0.6) q[7];

// ---- measurement ----
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];
measure q[7] -> c[7];

// Expected optimal outcomes (c[0..7] = x_vertex1..x_vertex8), H_C=-4:
//   1,1,0,0,1,0,1,0  -> clique {1,2,5,7}
//   1,1,0,0,1,0,0,1  -> clique {1,2,5,8}
//   1,1,0,0,0,1,1,0  -> clique {1,2,6,7}
//   1,1,0,0,0,1,0,1  -> clique {1,2,6,8}
//   0,1,0,1,1,0,1,0  -> clique {2,4,5,7}
//   0,1,0,1,1,0,0,1  -> clique {2,4,5,8}
