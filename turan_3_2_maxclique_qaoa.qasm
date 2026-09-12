// ============================================================
// QAOA circuit for Maximum Clique on the Turan graph T(3,2)
// ============================================================
// T(3,2): 3 vertices split into parts {1,2} and {3}.
// Complete multipartite -> edges only across parts:
//   Edges:     (1,3), (2,3)
//   Non-edge:  (1,2)
// Max clique size = 2, achieved by {1,3} or {2,3}.
// (T(3,2) is Turan's extremal triangle-free graph on 3 vertices.)
//
// QUBO (Lucas 2014 Max-Clique formulation):
//   x_i = 1  iff vertex i is in the clique
//   A = 1 (reward per vertex), B = 2 (penalty per missing edge, B>A)
//   H = -A(x1+x2+x3) + B*x1*x2      <- only non-edge pair (1,2)
//     = -(x1+x2+x3) + 2*x1*x2
//
// Ising mapping x_i = (1 - Z_i)/2, dropping the additive constant:
//   H_C = 0.5*Z3 + 0.5*Z1*Z2
//
// Verified minimum of H_C over all 8 assignments occurs at
//   (x1,x2,x3) = (1,0,1) or (0,1,1)   <-- the two maximum cliques
// (1,1,1) is penalized for the missing edge (1,2), as required.
//
// Qubit map: q[0]=vertex1 (x1), q[1]=vertex2 (x2), q[2]=vertex3 (x3)
//
// p=1 QAOA:
//   |+++>  --exp(-i*gamma*H_C)-->  --exp(-i*beta*sum_i X_i)-->  measure
//
// Gate identities used:
//   exp(-i*theta/2 * Z_a Z_b)  =  CX(a,b); RZ(theta) b; CX(a,b)
//   exp(-i*theta/2 * Z_a)      =  RZ(theta) a
//   exp(-i*theta/2 * X_a)      =  RX(theta) a
//
// So with H_C = 0.5*Z3 + 0.5*Z1*Z2 :
//   ZZ term (q0,q1): CX; RZ(gamma); CX      [theta = gamma]
//   Z  term (q2):    RZ(gamma)              [theta = gamma]
//   mixer (each qubit): RX(2*beta)          [theta = 2*beta]
//
// Example variational parameters used below (in practice these are
// tuned by a classical outer-loop optimizer, e.g. COBYLA/SPSA):
//   gamma = 0.5   ->  rz(0.5)
//   beta  = 0.4   ->  rx(0.8)
// ============================================================

OPENQASM 2.0;
include "qelib1.inc";

qreg q[3];
creg c[3];

// ---- initial state: equal superposition over all bitstrings ----
h q[0];
h q[1];
h q[2];

// ---- cost unitary exp(-i*gamma*H_C) ----
// ZZ interaction term: 0.5 * Z1 * Z2  (qubits q[0], q[1])
cx q[0],q[1];
rz(0.5) q[1];        // rz(gamma), gamma = 0.5
cx q[0],q[1];

// single-Z term: 0.5 * Z3  (qubit q[2])
rz(0.5) q[2];        // rz(gamma), gamma = 0.5

// ---- mixer unitary exp(-i*beta*sum_i X_i) ----
rx(0.8) q[0];         // rx(2*beta), beta = 0.4
rx(0.8) q[1];
rx(0.8) q[2];

// ---- measurement ----
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];

// Expected dominant outcomes after optimizing (gamma,beta):
//   c[0],c[1],c[2] = 1,0,1  -> clique {vertex1, vertex3}
//   c[0],c[1],c[2] = 0,1,1  -> clique {vertex2, vertex3}
