// ============================================================
// QAOA circuit: extend Exoo's 42-vertex R(5,5)-critical graph
// by one new (43rd) vertex, searching for a connection pattern
// that creates NO new monochromatic K5 (red or blue).
//
// Base graph: Exoo(42) = Cyclic(43) minus vertex 0, with 16
// edges flipped red->blue (G. Exoo, J. Graph Theory 13 (1989)).
// Independently rebuilt and verified here: Exoo(42) has NO red
// K5 and NO blue K5, matching the published result.
//
// Constraint census: 1148 red K4's (new vertex must NOT
// connect red to all 4) and 1170 blue K4's / red-4-
// independent-sets (new vertex must NOT connect blue to all 4).
// Total: 2318 penalty terms.
//
// Qubit map: q[v-1] = 'new vertex connects (red) to existing
// vertex v', for v = 1..42.  q[42],q[43],q[44] = shared ancilla
// (compute AND of 4 literals, apply phase, uncompute; reused
// across every constraint).
//
// Each red-K4 {a,b,c,d} constraint:
//   CCX q[a],q[b] -> anc1 ; CCX anc1,q[c] -> anc2 ;
//   CCX anc2,q[d] -> anc3 ; RZ(gamma) anc3 ; uncompute (reverse)
// Each blue-K4 {a,b,c,d} constraint: same, but each of the 4
// qubits is X-flipped before/after (penalizing all-0 instead
// of all-1, since 'not connected' = blue = would complete a
// blue K5).
//
// gamma = 0.15 (cost angle), beta = 0.3 (mixer angle) -- 
// placeholders; a real run tunes these via a classical outer loop.
// ============================================================

OPENQASM 2.0;
include "qelib1.inc";

qreg q[45];
creg c[45];

// ---- initial state: equal superposition over the 42 connection qubits ----
h q[0];
h q[1];
h q[2];
h q[3];
h q[4];
h q[5];
h q[6];
h q[7];
h q[8];
h q[9];
h q[10];
h q[11];
h q[12];
h q[13];
h q[14];
h q[15];
h q[16];
h q[17];
h q[18];
h q[19];
h q[20];
h q[21];
h q[22];
h q[23];
h q[24];
h q[25];
h q[26];
h q[27];
h q[28];
h q[29];
h q[30];
h q[31];
h q[32];
h q[33];
h q[34];
h q[35];
h q[36];
h q[37];
h q[38];
h q[39];
h q[40];
h q[41];

// ---- cost unitary: one penalty block per forbidden K4 ----
// red K4 (1, 2, 3, 15) -> avoid connecting to all four
ccx q[0],q[1],q[42];
ccx q[42],q[2],q[43];
ccx q[43],q[14],q[44];
rz(0.15) q[44];
ccx q[43],q[14],q[44];
ccx q[42],q[2],q[43];
ccx q[0],q[1],q[42];

// red K4 (1, 2, 3, 23) -> avoid connecting to all four
ccx q[0],q[1],q[42];
ccx q[42],q[2],q[43];
ccx q[43],q[22],q[44];
rz(0.15) q[44];
ccx q[43],q[22],q[44];
ccx q[42],q[2],q[43];
ccx q[0],q[1],q[42];

// red K4 (1, 2, 3, 24) -> avoid connecting to all four
ccx q[0],q[1],q[42];
ccx q[42],q[2],q[43];
ccx q[43],q[23],q[44];
rz(0.15) q[44];
ccx q[43],q[23],q[44];
ccx q[42],q[2],q[43];
ccx q[0],q[1],q[42];

// red K4 (1, 2, 3, 32) -> avoid connecting to all four
ccx q[0],q[1],q[42];
ccx q[42],q[2],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[2],q[43];
ccx q[0],q[1],q[42];

// red K4 (1, 2, 14, 24) -> avoid connecting to all four
ccx q[0],q[1],q[42];
ccx q[42],q[13],q[43];
ccx q[43],q[23],q[44];
rz(0.15) q[44];
ccx q[43],q[23],q[44];
ccx q[42],q[13],q[43];
ccx q[0],q[1],q[42];

// red K4 (1, 2, 14, 32) -> avoid connecting to all four
ccx q[0],q[1],q[42];
ccx q[42],q[13],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[13],q[43];
ccx q[0],q[1],q[42];

// red K4 (1, 2, 15, 22) -> avoid connecting to all four
ccx q[0],q[1],q[42];
ccx q[42],q[14],q[43];
ccx q[43],q[21],q[44];
rz(0.15) q[44];
ccx q[43],q[21],q[44];
ccx q[42],q[14],q[43];
ccx q[0],q[1],q[42];

// red K4 (1, 2, 15, 31) -> avoid connecting to all four
ccx q[0],q[1],q[42];
ccx q[42],q[14],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[14],q[43];
ccx q[0],q[1],q[42];

// red K4 (1, 2, 22, 23) -> avoid connecting to all four
ccx q[0],q[1],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[22],q[44];
rz(0.15) q[44];
ccx q[43],q[22],q[44];
ccx q[42],q[21],q[43];
ccx q[0],q[1],q[42];

// red K4 (1, 2, 22, 24) -> avoid connecting to all four
ccx q[0],q[1],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[23],q[44];
rz(0.15) q[44];
ccx q[43],q[23],q[44];
ccx q[42],q[21],q[43];
ccx q[0],q[1],q[42];

// red K4 (1, 2, 22, 32) -> avoid connecting to all four
ccx q[0],q[1],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[21],q[43];
ccx q[0],q[1],q[42];

// red K4 (1, 2, 24, 31) -> avoid connecting to all four
ccx q[0],q[1],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[23],q[43];
ccx q[0],q[1],q[42];

// red K4 (1, 2, 31, 32) -> avoid connecting to all four
ccx q[0],q[1],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[30],q[43];
ccx q[0],q[1],q[42];

// red K4 (1, 3, 13, 15) -> avoid connecting to all four
ccx q[0],q[2],q[42];
ccx q[42],q[12],q[43];
ccx q[43],q[14],q[44];
rz(0.15) q[44];
ccx q[43],q[14],q[44];
ccx q[42],q[12],q[43];
ccx q[0],q[2],q[42];

// red K4 (1, 3, 13, 23) -> avoid connecting to all four
ccx q[0],q[2],q[42];
ccx q[42],q[12],q[43];
ccx q[43],q[22],q[44];
rz(0.15) q[44];
ccx q[43],q[22],q[44];
ccx q[42],q[12],q[43];
ccx q[0],q[2],q[42];

// red K4 (1, 3, 13, 26) -> avoid connecting to all four
ccx q[0],q[2],q[42];
ccx q[42],q[12],q[43];
ccx q[43],q[25],q[44];
rz(0.15) q[44];
ccx q[43],q[25],q[44];
ccx q[42],q[12],q[43];
ccx q[0],q[2],q[42];

// red K4 (1, 3, 13, 34) -> avoid connecting to all four
ccx q[0],q[2],q[42];
ccx q[42],q[12],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[12],q[43];
ccx q[0],q[2],q[42];

// red K4 (1, 3, 15, 17) -> avoid connecting to all four
ccx q[0],q[2],q[42];
ccx q[42],q[14],q[43];
ccx q[43],q[16],q[44];
rz(0.15) q[44];
ccx q[43],q[16],q[44];
ccx q[42],q[14],q[43];
ccx q[0],q[2],q[42];

// red K4 (1, 3, 15, 28) -> avoid connecting to all four
ccx q[0],q[2],q[42];
ccx q[42],q[14],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[14],q[43];
ccx q[0],q[2],q[42];

// red K4 (1, 3, 17, 19) -> avoid connecting to all four
ccx q[0],q[2],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[18],q[44];
rz(0.15) q[44];
ccx q[43],q[18],q[44];
ccx q[42],q[16],q[43];
ccx q[0],q[2],q[42];

// red K4 (1, 3, 17, 24) -> avoid connecting to all four
ccx q[0],q[2],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[23],q[44];
rz(0.15) q[44];
ccx q[43],q[23],q[44];
ccx q[42],q[16],q[43];
ccx q[0],q[2],q[42];

// red K4 (1, 3, 17, 30) -> avoid connecting to all four
ccx q[0],q[2],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[16],q[43];
ccx q[0],q[2],q[42];

// red K4 (1, 3, 19, 21) -> avoid connecting to all four
ccx q[0],q[2],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[20],q[44];
rz(0.15) q[44];
ccx q[43],q[20],q[44];
ccx q[42],q[18],q[43];
ccx q[0],q[2],q[42];

// red K4 (1, 3, 19, 26) -> avoid connecting to all four
ccx q[0],q[2],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[25],q[44];
rz(0.15) q[44];
ccx q[43],q[25],q[44];
ccx q[42],q[18],q[43];
ccx q[0],q[2],q[42];

// red K4 (1, 3, 19, 32) -> avoid connecting to all four
ccx q[0],q[2],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[18],q[43];
ccx q[0],q[2],q[42];

// red K4 (1, 3, 21, 23) -> avoid connecting to all four
ccx q[0],q[2],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[22],q[44];
rz(0.15) q[44];
ccx q[43],q[22],q[44];
ccx q[42],q[20],q[43];
ccx q[0],q[2],q[42];

// red K4 (1, 3, 21, 28) -> avoid connecting to all four
ccx q[0],q[2],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[20],q[43];
ccx q[0],q[2],q[42];

// red K4 (1, 3, 21, 34) -> avoid connecting to all four
ccx q[0],q[2],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[20],q[43];
ccx q[0],q[2],q[42];

// red K4 (1, 3, 23, 30) -> avoid connecting to all four
ccx q[0],q[2],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[22],q[43];
ccx q[0],q[2],q[42];

// red K4 (1, 3, 24, 26) -> avoid connecting to all four
ccx q[0],q[2],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[25],q[44];
rz(0.15) q[44];
ccx q[43],q[25],q[44];
ccx q[42],q[23],q[43];
ccx q[0],q[2],q[42];

// red K4 (1, 3, 24, 34) -> avoid connecting to all four
ccx q[0],q[2],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[23],q[43];
ccx q[0],q[2],q[42];

// red K4 (1, 3, 26, 28) -> avoid connecting to all four
ccx q[0],q[2],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[25],q[43];
ccx q[0],q[2],q[42];

// red K4 (1, 3, 28, 30) -> avoid connecting to all four
ccx q[0],q[2],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[27],q[43];
ccx q[0],q[2],q[42];

// red K4 (1, 3, 30, 32) -> avoid connecting to all four
ccx q[0],q[2],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[29],q[43];
ccx q[0],q[2],q[42];

// red K4 (1, 3, 32, 34) -> avoid connecting to all four
ccx q[0],q[2],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[31],q[43];
ccx q[0],q[2],q[42];

// red K4 (1, 8, 15, 22) -> avoid connecting to all four
ccx q[0],q[7],q[42];
ccx q[42],q[14],q[43];
ccx q[43],q[21],q[44];
rz(0.15) q[44];
ccx q[43],q[21],q[44];
ccx q[42],q[14],q[43];
ccx q[0],q[7],q[42];

// red K4 (1, 8, 15, 28) -> avoid connecting to all four
ccx q[0],q[7],q[42];
ccx q[42],q[14],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[14],q[43];
ccx q[0],q[7],q[42];

// red K4 (1, 8, 15, 31) -> avoid connecting to all four
ccx q[0],q[7],q[42];
ccx q[42],q[14],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[14],q[43];
ccx q[0],q[7],q[42];

// red K4 (1, 8, 15, 37) -> avoid connecting to all four
ccx q[0],q[7],q[42];
ccx q[42],q[14],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[14],q[43];
ccx q[0],q[7],q[42];

// red K4 (1, 8, 21, 22) -> avoid connecting to all four
ccx q[0],q[7],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[21],q[44];
rz(0.15) q[44];
ccx q[43],q[21],q[44];
ccx q[42],q[20],q[43];
ccx q[0],q[7],q[42];

// red K4 (1, 8, 21, 28) -> avoid connecting to all four
ccx q[0],q[7],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[20],q[43];
ccx q[0],q[7],q[42];

// red K4 (1, 8, 21, 31) -> avoid connecting to all four
ccx q[0],q[7],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[20],q[43];
ccx q[0],q[7],q[42];

// red K4 (1, 8, 21, 37) -> avoid connecting to all four
ccx q[0],q[7],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[20],q[43];
ccx q[0],q[7],q[42];

// red K4 (1, 8, 22, 24) -> avoid connecting to all four
ccx q[0],q[7],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[23],q[44];
rz(0.15) q[44];
ccx q[43],q[23],q[44];
ccx q[42],q[21],q[43];
ccx q[0],q[7],q[42];

// red K4 (1, 8, 24, 26) -> avoid connecting to all four
ccx q[0],q[7],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[25],q[44];
rz(0.15) q[44];
ccx q[43],q[25],q[44];
ccx q[42],q[23],q[43];
ccx q[0],q[7],q[42];

// red K4 (1, 8, 24, 31) -> avoid connecting to all four
ccx q[0],q[7],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[23],q[43];
ccx q[0],q[7],q[42];

// red K4 (1, 8, 24, 37) -> avoid connecting to all four
ccx q[0],q[7],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[23],q[43];
ccx q[0],q[7],q[42];

// red K4 (1, 8, 26, 28) -> avoid connecting to all four
ccx q[0],q[7],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[25],q[43];
ccx q[0],q[7],q[42];

// red K4 (1, 8, 28, 30) -> avoid connecting to all four
ccx q[0],q[7],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[27],q[43];
ccx q[0],q[7],q[42];

// red K4 (1, 8, 30, 37) -> avoid connecting to all four
ccx q[0],q[7],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[29],q[43];
ccx q[0],q[7],q[42];

// red K4 (1, 11, 13, 23) -> avoid connecting to all four
ccx q[0],q[10],q[42];
ccx q[42],q[12],q[43];
ccx q[43],q[22],q[44];
rz(0.15) q[44];
ccx q[43],q[22],q[44];
ccx q[42],q[12],q[43];
ccx q[0],q[10],q[42];

// red K4 (1, 11, 13, 31) -> avoid connecting to all four
ccx q[0],q[10],q[42];
ccx q[42],q[12],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[12],q[43];
ccx q[0],q[10],q[42];

// red K4 (1, 11, 13, 34) -> avoid connecting to all four
ccx q[0],q[10],q[42];
ccx q[42],q[12],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[12],q[43];
ccx q[0],q[10],q[42];

// red K4 (1, 11, 13, 42) -> avoid connecting to all four
ccx q[0],q[10],q[42];
ccx q[42],q[12],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[12],q[43];
ccx q[0],q[10],q[42];

// red K4 (1, 11, 21, 23) -> avoid connecting to all four
ccx q[0],q[10],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[22],q[44];
rz(0.15) q[44];
ccx q[43],q[22],q[44];
ccx q[42],q[20],q[43];
ccx q[0],q[10],q[42];

// red K4 (1, 11, 21, 31) -> avoid connecting to all four
ccx q[0],q[10],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[20],q[43];
ccx q[0],q[10],q[42];

// red K4 (1, 11, 21, 34) -> avoid connecting to all four
ccx q[0],q[10],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[20],q[43];
ccx q[0],q[10],q[42];

// red K4 (1, 11, 21, 42) -> avoid connecting to all four
ccx q[0],q[10],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[20],q[43];
ccx q[0],q[10],q[42];

// red K4 (1, 11, 24, 31) -> avoid connecting to all four
ccx q[0],q[10],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[23],q[43];
ccx q[0],q[10],q[42];

// red K4 (1, 11, 24, 34) -> avoid connecting to all four
ccx q[0],q[10],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[23],q[43];
ccx q[0],q[10],q[42];

// red K4 (1, 11, 24, 42) -> avoid connecting to all four
ccx q[0],q[10],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[23],q[43];
ccx q[0],q[10],q[42];

// red K4 (1, 13, 15, 31) -> avoid connecting to all four
ccx q[0],q[12],q[42];
ccx q[42],q[14],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[14],q[43];
ccx q[0],q[12],q[42];

// red K4 (1, 13, 15, 42) -> avoid connecting to all four
ccx q[0],q[12],q[42];
ccx q[42],q[14],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[14],q[43];
ccx q[0],q[12],q[42];

// red K4 (1, 13, 26, 42) -> avoid connecting to all four
ccx q[0],q[12],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[25],q[43];
ccx q[0],q[12],q[42];

// red K4 (1, 14, 21, 28) -> avoid connecting to all four
ccx q[0],q[13],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[20],q[43];
ccx q[0],q[13],q[42];

// red K4 (1, 14, 21, 34) -> avoid connecting to all four
ccx q[0],q[13],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[20],q[43];
ccx q[0],q[13],q[42];

// red K4 (1, 14, 21, 37) -> avoid connecting to all four
ccx q[0],q[13],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[20],q[43];
ccx q[0],q[13],q[42];

// red K4 (1, 14, 24, 26) -> avoid connecting to all four
ccx q[0],q[13],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[25],q[44];
rz(0.15) q[44];
ccx q[43],q[25],q[44];
ccx q[42],q[23],q[43];
ccx q[0],q[13],q[42];

// red K4 (1, 14, 24, 34) -> avoid connecting to all four
ccx q[0],q[13],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[23],q[43];
ccx q[0],q[13],q[42];

// red K4 (1, 14, 24, 37) -> avoid connecting to all four
ccx q[0],q[13],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[23],q[43];
ccx q[0],q[13],q[42];

// red K4 (1, 14, 26, 28) -> avoid connecting to all four
ccx q[0],q[13],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[25],q[43];
ccx q[0],q[13],q[42];

// red K4 (1, 14, 28, 30) -> avoid connecting to all four
ccx q[0],q[13],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[27],q[43];
ccx q[0],q[13],q[42];

// red K4 (1, 14, 30, 32) -> avoid connecting to all four
ccx q[0],q[13],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[29],q[43];
ccx q[0],q[13],q[42];

// red K4 (1, 14, 30, 37) -> avoid connecting to all four
ccx q[0],q[13],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[29],q[43];
ccx q[0],q[13],q[42];

// red K4 (1, 14, 32, 34) -> avoid connecting to all four
ccx q[0],q[13],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[31],q[43];
ccx q[0],q[13],q[42];

// red K4 (1, 15, 17, 31) -> avoid connecting to all four
ccx q[0],q[14],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[16],q[43];
ccx q[0],q[14],q[42];

// red K4 (1, 15, 17, 37) -> avoid connecting to all four
ccx q[0],q[14],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[16],q[43];
ccx q[0],q[14],q[42];

// red K4 (1, 15, 17, 42) -> avoid connecting to all four
ccx q[0],q[14],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[16],q[43];
ccx q[0],q[14],q[42];

// red K4 (1, 15, 22, 42) -> avoid connecting to all four
ccx q[0],q[14],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[21],q[43];
ccx q[0],q[14],q[42];

// red K4 (1, 15, 28, 42) -> avoid connecting to all four
ccx q[0],q[14],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[27],q[43];
ccx q[0],q[14],q[42];

// red K4 (1, 17, 19, 31) -> avoid connecting to all four
ccx q[0],q[16],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[18],q[43];
ccx q[0],q[16],q[42];

// red K4 (1, 17, 19, 37) -> avoid connecting to all four
ccx q[0],q[16],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[18],q[43];
ccx q[0],q[16],q[42];

// red K4 (1, 17, 19, 42) -> avoid connecting to all four
ccx q[0],q[16],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[18],q[43];
ccx q[0],q[16],q[42];

// red K4 (1, 17, 24, 31) -> avoid connecting to all four
ccx q[0],q[16],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[23],q[43];
ccx q[0],q[16],q[42];

// red K4 (1, 17, 24, 37) -> avoid connecting to all four
ccx q[0],q[16],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[23],q[43];
ccx q[0],q[16],q[42];

// red K4 (1, 17, 24, 42) -> avoid connecting to all four
ccx q[0],q[16],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[23],q[43];
ccx q[0],q[16],q[42];

// red K4 (1, 17, 30, 37) -> avoid connecting to all four
ccx q[0],q[16],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[29],q[43];
ccx q[0],q[16],q[42];

// red K4 (1, 17, 30, 42) -> avoid connecting to all four
ccx q[0],q[16],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[29],q[43];
ccx q[0],q[16],q[42];

// red K4 (1, 19, 21, 31) -> avoid connecting to all four
ccx q[0],q[18],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[20],q[43];
ccx q[0],q[18],q[42];

// red K4 (1, 19, 21, 37) -> avoid connecting to all four
ccx q[0],q[18],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[20],q[43];
ccx q[0],q[18],q[42];

// red K4 (1, 19, 21, 42) -> avoid connecting to all four
ccx q[0],q[18],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[20],q[43];
ccx q[0],q[18],q[42];

// red K4 (1, 19, 26, 42) -> avoid connecting to all four
ccx q[0],q[18],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[25],q[43];
ccx q[0],q[18],q[42];

// red K4 (1, 19, 31, 32) -> avoid connecting to all four
ccx q[0],q[18],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[30],q[43];
ccx q[0],q[18],q[42];

// red K4 (1, 19, 32, 42) -> avoid connecting to all four
ccx q[0],q[18],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[31],q[43];
ccx q[0],q[18],q[42];

// red K4 (1, 21, 22, 23) -> avoid connecting to all four
ccx q[0],q[20],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[22],q[44];
rz(0.15) q[44];
ccx q[43],q[22],q[44];
ccx q[42],q[21],q[43];
ccx q[0],q[20],q[42];

// red K4 (1, 21, 22, 34) -> avoid connecting to all four
ccx q[0],q[20],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[21],q[43];
ccx q[0],q[20],q[42];

// red K4 (1, 21, 22, 42) -> avoid connecting to all four
ccx q[0],q[20],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[21],q[43];
ccx q[0],q[20],q[42];

// red K4 (1, 21, 23, 37) -> avoid connecting to all four
ccx q[0],q[20],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[22],q[43];
ccx q[0],q[20],q[42];

// red K4 (1, 21, 28, 42) -> avoid connecting to all four
ccx q[0],q[20],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[27],q[43];
ccx q[0],q[20],q[42];

// red K4 (1, 22, 24, 34) -> avoid connecting to all four
ccx q[0],q[21],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[23],q[43];
ccx q[0],q[21],q[42];

// red K4 (1, 22, 24, 42) -> avoid connecting to all four
ccx q[0],q[21],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[23],q[43];
ccx q[0],q[21],q[42];

// red K4 (1, 22, 32, 34) -> avoid connecting to all four
ccx q[0],q[21],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[31],q[43];
ccx q[0],q[21],q[42];

// red K4 (1, 22, 32, 42) -> avoid connecting to all four
ccx q[0],q[21],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[31],q[43];
ccx q[0],q[21],q[42];

// red K4 (1, 23, 30, 37) -> avoid connecting to all four
ccx q[0],q[22],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[29],q[43];
ccx q[0],q[22],q[42];

// red K4 (1, 24, 26, 42) -> avoid connecting to all four
ccx q[0],q[23],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[25],q[43];
ccx q[0],q[23],q[42];

// red K4 (1, 26, 28, 42) -> avoid connecting to all four
ccx q[0],q[25],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[27],q[43];
ccx q[0],q[25],q[42];

// red K4 (1, 28, 30, 42) -> avoid connecting to all four
ccx q[0],q[27],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[29],q[43];
ccx q[0],q[27],q[42];

// red K4 (1, 30, 32, 42) -> avoid connecting to all four
ccx q[0],q[29],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[31],q[43];
ccx q[0],q[29],q[42];

// red K4 (2, 3, 4, 16) -> avoid connecting to all four
ccx q[1],q[2],q[42];
ccx q[42],q[3],q[43];
ccx q[43],q[15],q[44];
rz(0.15) q[44];
ccx q[43],q[15],q[44];
ccx q[42],q[3],q[43];
ccx q[1],q[2],q[42];

// red K4 (2, 3, 4, 24) -> avoid connecting to all four
ccx q[1],q[2],q[42];
ccx q[42],q[3],q[43];
ccx q[43],q[23],q[44];
rz(0.15) q[44];
ccx q[43],q[23],q[44];
ccx q[42],q[3],q[43];
ccx q[1],q[2],q[42];

// red K4 (2, 3, 4, 25) -> avoid connecting to all four
ccx q[1],q[2],q[42];
ccx q[42],q[3],q[43];
ccx q[43],q[24],q[44];
rz(0.15) q[44];
ccx q[43],q[24],q[44];
ccx q[42],q[3],q[43];
ccx q[1],q[2],q[42];

// red K4 (2, 3, 4, 33) -> avoid connecting to all four
ccx q[1],q[2],q[42];
ccx q[42],q[3],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[3],q[43];
ccx q[1],q[2],q[42];

// red K4 (2, 3, 15, 25) -> avoid connecting to all four
ccx q[1],q[2],q[42];
ccx q[42],q[14],q[43];
ccx q[43],q[24],q[44];
rz(0.15) q[44];
ccx q[43],q[24],q[44];
ccx q[42],q[14],q[43];
ccx q[1],q[2],q[42];

// red K4 (2, 3, 15, 33) -> avoid connecting to all four
ccx q[1],q[2],q[42];
ccx q[42],q[14],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[14],q[43];
ccx q[1],q[2],q[42];

// red K4 (2, 3, 16, 23) -> avoid connecting to all four
ccx q[1],q[2],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[22],q[44];
rz(0.15) q[44];
ccx q[43],q[22],q[44];
ccx q[42],q[15],q[43];
ccx q[1],q[2],q[42];

// red K4 (2, 3, 16, 32) -> avoid connecting to all four
ccx q[1],q[2],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[15],q[43];
ccx q[1],q[2],q[42];

// red K4 (2, 3, 23, 25) -> avoid connecting to all four
ccx q[1],q[2],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[24],q[44];
rz(0.15) q[44];
ccx q[43],q[24],q[44];
ccx q[42],q[22],q[43];
ccx q[1],q[2],q[42];

// red K4 (2, 3, 23, 33) -> avoid connecting to all four
ccx q[1],q[2],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[22],q[43];
ccx q[1],q[2],q[42];

// red K4 (2, 3, 25, 32) -> avoid connecting to all four
ccx q[1],q[2],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[24],q[43];
ccx q[1],q[2],q[42];

// red K4 (2, 3, 32, 33) -> avoid connecting to all four
ccx q[1],q[2],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[31],q[43];
ccx q[1],q[2],q[42];

// red K4 (2, 4, 14, 16) -> avoid connecting to all four
ccx q[1],q[3],q[42];
ccx q[42],q[13],q[43];
ccx q[43],q[15],q[44];
rz(0.15) q[44];
ccx q[43],q[15],q[44];
ccx q[42],q[13],q[43];
ccx q[1],q[3],q[42];

// red K4 (2, 4, 14, 24) -> avoid connecting to all four
ccx q[1],q[3],q[42];
ccx q[42],q[13],q[43];
ccx q[43],q[23],q[44];
rz(0.15) q[44];
ccx q[43],q[23],q[44];
ccx q[42],q[13],q[43];
ccx q[1],q[3],q[42];

// red K4 (2, 4, 14, 27) -> avoid connecting to all four
ccx q[1],q[3],q[42];
ccx q[42],q[13],q[43];
ccx q[43],q[26],q[44];
rz(0.15) q[44];
ccx q[43],q[26],q[44];
ccx q[42],q[13],q[43];
ccx q[1],q[3],q[42];

// red K4 (2, 4, 14, 35) -> avoid connecting to all four
ccx q[1],q[3],q[42];
ccx q[42],q[13],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[13],q[43];
ccx q[1],q[3],q[42];

// red K4 (2, 4, 16, 18) -> avoid connecting to all four
ccx q[1],q[3],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[17],q[44];
rz(0.15) q[44];
ccx q[43],q[17],q[44];
ccx q[42],q[15],q[43];
ccx q[1],q[3],q[42];

// red K4 (2, 4, 16, 29) -> avoid connecting to all four
ccx q[1],q[3],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[15],q[43];
ccx q[1],q[3],q[42];

// red K4 (2, 4, 18, 20) -> avoid connecting to all four
ccx q[1],q[3],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[19],q[44];
rz(0.15) q[44];
ccx q[43],q[19],q[44];
ccx q[42],q[17],q[43];
ccx q[1],q[3],q[42];

// red K4 (2, 4, 18, 25) -> avoid connecting to all four
ccx q[1],q[3],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[24],q[44];
rz(0.15) q[44];
ccx q[43],q[24],q[44];
ccx q[42],q[17],q[43];
ccx q[1],q[3],q[42];

// red K4 (2, 4, 18, 31) -> avoid connecting to all four
ccx q[1],q[3],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[17],q[43];
ccx q[1],q[3],q[42];

// red K4 (2, 4, 20, 22) -> avoid connecting to all four
ccx q[1],q[3],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[21],q[44];
rz(0.15) q[44];
ccx q[43],q[21],q[44];
ccx q[42],q[19],q[43];
ccx q[1],q[3],q[42];

// red K4 (2, 4, 20, 27) -> avoid connecting to all four
ccx q[1],q[3],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[26],q[44];
rz(0.15) q[44];
ccx q[43],q[26],q[44];
ccx q[42],q[19],q[43];
ccx q[1],q[3],q[42];

// red K4 (2, 4, 20, 33) -> avoid connecting to all four
ccx q[1],q[3],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[19],q[43];
ccx q[1],q[3],q[42];

// red K4 (2, 4, 22, 24) -> avoid connecting to all four
ccx q[1],q[3],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[23],q[44];
rz(0.15) q[44];
ccx q[43],q[23],q[44];
ccx q[42],q[21],q[43];
ccx q[1],q[3],q[42];

// red K4 (2, 4, 22, 29) -> avoid connecting to all four
ccx q[1],q[3],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[21],q[43];
ccx q[1],q[3],q[42];

// red K4 (2, 4, 22, 35) -> avoid connecting to all four
ccx q[1],q[3],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[21],q[43];
ccx q[1],q[3],q[42];

// red K4 (2, 4, 24, 31) -> avoid connecting to all four
ccx q[1],q[3],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[23],q[43];
ccx q[1],q[3],q[42];

// red K4 (2, 4, 25, 27) -> avoid connecting to all four
ccx q[1],q[3],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[26],q[44];
rz(0.15) q[44];
ccx q[43],q[26],q[44];
ccx q[42],q[24],q[43];
ccx q[1],q[3],q[42];

// red K4 (2, 4, 25, 35) -> avoid connecting to all four
ccx q[1],q[3],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[24],q[43];
ccx q[1],q[3],q[42];

// red K4 (2, 4, 27, 29) -> avoid connecting to all four
ccx q[1],q[3],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[26],q[43];
ccx q[1],q[3],q[42];

// red K4 (2, 4, 29, 31) -> avoid connecting to all four
ccx q[1],q[3],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[28],q[43];
ccx q[1],q[3],q[42];

// red K4 (2, 4, 31, 33) -> avoid connecting to all four
ccx q[1],q[3],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[30],q[43];
ccx q[1],q[3],q[42];

// red K4 (2, 4, 33, 35) -> avoid connecting to all four
ccx q[1],q[3],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[32],q[43];
ccx q[1],q[3],q[42];

// red K4 (2, 9, 16, 23) -> avoid connecting to all four
ccx q[1],q[8],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[22],q[44];
rz(0.15) q[44];
ccx q[43],q[22],q[44];
ccx q[42],q[15],q[43];
ccx q[1],q[8],q[42];

// red K4 (2, 9, 16, 29) -> avoid connecting to all four
ccx q[1],q[8],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[15],q[43];
ccx q[1],q[8],q[42];

// red K4 (2, 9, 16, 32) -> avoid connecting to all four
ccx q[1],q[8],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[15],q[43];
ccx q[1],q[8],q[42];

// red K4 (2, 9, 16, 38) -> avoid connecting to all four
ccx q[1],q[8],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[15],q[43];
ccx q[1],q[8],q[42];

// red K4 (2, 9, 22, 23) -> avoid connecting to all four
ccx q[1],q[8],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[22],q[44];
rz(0.15) q[44];
ccx q[43],q[22],q[44];
ccx q[42],q[21],q[43];
ccx q[1],q[8],q[42];

// red K4 (2, 9, 22, 29) -> avoid connecting to all four
ccx q[1],q[8],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[21],q[43];
ccx q[1],q[8],q[42];

// red K4 (2, 9, 22, 32) -> avoid connecting to all four
ccx q[1],q[8],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[21],q[43];
ccx q[1],q[8],q[42];

// red K4 (2, 9, 22, 38) -> avoid connecting to all four
ccx q[1],q[8],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[21],q[43];
ccx q[1],q[8],q[42];

// red K4 (2, 9, 23, 25) -> avoid connecting to all four
ccx q[1],q[8],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[24],q[44];
rz(0.15) q[44];
ccx q[43],q[24],q[44];
ccx q[42],q[22],q[43];
ccx q[1],q[8],q[42];

// red K4 (2, 9, 25, 27) -> avoid connecting to all four
ccx q[1],q[8],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[26],q[44];
rz(0.15) q[44];
ccx q[43],q[26],q[44];
ccx q[42],q[24],q[43];
ccx q[1],q[8],q[42];

// red K4 (2, 9, 25, 32) -> avoid connecting to all four
ccx q[1],q[8],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[24],q[43];
ccx q[1],q[8],q[42];

// red K4 (2, 9, 25, 38) -> avoid connecting to all four
ccx q[1],q[8],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[24],q[43];
ccx q[1],q[8],q[42];

// red K4 (2, 9, 27, 29) -> avoid connecting to all four
ccx q[1],q[8],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[26],q[43];
ccx q[1],q[8],q[42];

// red K4 (2, 9, 29, 31) -> avoid connecting to all four
ccx q[1],q[8],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[28],q[43];
ccx q[1],q[8],q[42];

// red K4 (2, 9, 31, 32) -> avoid connecting to all four
ccx q[1],q[8],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[30],q[43];
ccx q[1],q[8],q[42];

// red K4 (2, 9, 31, 38) -> avoid connecting to all four
ccx q[1],q[8],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[30],q[43];
ccx q[1],q[8],q[42];

// red K4 (2, 12, 14, 24) -> avoid connecting to all four
ccx q[1],q[11],q[42];
ccx q[42],q[13],q[43];
ccx q[43],q[23],q[44];
rz(0.15) q[44];
ccx q[43],q[23],q[44];
ccx q[42],q[13],q[43];
ccx q[1],q[11],q[42];

// red K4 (2, 12, 14, 32) -> avoid connecting to all four
ccx q[1],q[11],q[42];
ccx q[42],q[13],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[13],q[43];
ccx q[1],q[11],q[42];

// red K4 (2, 12, 14, 35) -> avoid connecting to all four
ccx q[1],q[11],q[42];
ccx q[42],q[13],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[13],q[43];
ccx q[1],q[11],q[42];

// red K4 (2, 12, 22, 24) -> avoid connecting to all four
ccx q[1],q[11],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[23],q[44];
rz(0.15) q[44];
ccx q[43],q[23],q[44];
ccx q[42],q[21],q[43];
ccx q[1],q[11],q[42];

// red K4 (2, 12, 22, 32) -> avoid connecting to all four
ccx q[1],q[11],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[21],q[43];
ccx q[1],q[11],q[42];

// red K4 (2, 12, 22, 35) -> avoid connecting to all four
ccx q[1],q[11],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[21],q[43];
ccx q[1],q[11],q[42];

// red K4 (2, 12, 25, 32) -> avoid connecting to all four
ccx q[1],q[11],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[24],q[43];
ccx q[1],q[11],q[42];

// red K4 (2, 12, 25, 35) -> avoid connecting to all four
ccx q[1],q[11],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[24],q[43];
ccx q[1],q[11],q[42];

// red K4 (2, 12, 32, 33) -> avoid connecting to all four
ccx q[1],q[11],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[31],q[43];
ccx q[1],q[11],q[42];

// red K4 (2, 12, 33, 35) -> avoid connecting to all four
ccx q[1],q[11],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[32],q[43];
ccx q[1],q[11],q[42];

// red K4 (2, 14, 16, 32) -> avoid connecting to all four
ccx q[1],q[13],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[15],q[43];
ccx q[1],q[13],q[42];

// red K4 (2, 15, 22, 29) -> avoid connecting to all four
ccx q[1],q[14],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[21],q[43];
ccx q[1],q[14],q[42];

// red K4 (2, 15, 22, 35) -> avoid connecting to all four
ccx q[1],q[14],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[21],q[43];
ccx q[1],q[14],q[42];

// red K4 (2, 15, 22, 38) -> avoid connecting to all four
ccx q[1],q[14],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[21],q[43];
ccx q[1],q[14],q[42];

// red K4 (2, 15, 25, 27) -> avoid connecting to all four
ccx q[1],q[14],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[26],q[44];
rz(0.15) q[44];
ccx q[43],q[26],q[44];
ccx q[42],q[24],q[43];
ccx q[1],q[14],q[42];

// red K4 (2, 15, 25, 35) -> avoid connecting to all four
ccx q[1],q[14],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[24],q[43];
ccx q[1],q[14],q[42];

// red K4 (2, 15, 25, 38) -> avoid connecting to all four
ccx q[1],q[14],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[24],q[43];
ccx q[1],q[14],q[42];

// red K4 (2, 15, 27, 29) -> avoid connecting to all four
ccx q[1],q[14],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[26],q[43];
ccx q[1],q[14],q[42];

// red K4 (2, 15, 29, 31) -> avoid connecting to all four
ccx q[1],q[14],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[28],q[43];
ccx q[1],q[14],q[42];

// red K4 (2, 15, 31, 33) -> avoid connecting to all four
ccx q[1],q[14],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[30],q[43];
ccx q[1],q[14],q[42];

// red K4 (2, 15, 31, 38) -> avoid connecting to all four
ccx q[1],q[14],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[30],q[43];
ccx q[1],q[14],q[42];

// red K4 (2, 15, 33, 35) -> avoid connecting to all four
ccx q[1],q[14],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[32],q[43];
ccx q[1],q[14],q[42];

// red K4 (2, 16, 18, 32) -> avoid connecting to all four
ccx q[1],q[15],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[17],q[43];
ccx q[1],q[15],q[42];

// red K4 (2, 16, 18, 38) -> avoid connecting to all four
ccx q[1],q[15],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[17],q[43];
ccx q[1],q[15],q[42];

// red K4 (2, 18, 20, 32) -> avoid connecting to all four
ccx q[1],q[17],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[19],q[43];
ccx q[1],q[17],q[42];

// red K4 (2, 18, 20, 38) -> avoid connecting to all four
ccx q[1],q[17],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[19],q[43];
ccx q[1],q[17],q[42];

// red K4 (2, 18, 25, 32) -> avoid connecting to all four
ccx q[1],q[17],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[24],q[43];
ccx q[1],q[17],q[42];

// red K4 (2, 18, 25, 38) -> avoid connecting to all four
ccx q[1],q[17],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[24],q[43];
ccx q[1],q[17],q[42];

// red K4 (2, 18, 31, 32) -> avoid connecting to all four
ccx q[1],q[17],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[30],q[43];
ccx q[1],q[17],q[42];

// red K4 (2, 18, 31, 38) -> avoid connecting to all four
ccx q[1],q[17],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[30],q[43];
ccx q[1],q[17],q[42];

// red K4 (2, 20, 22, 32) -> avoid connecting to all four
ccx q[1],q[19],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[21],q[43];
ccx q[1],q[19],q[42];

// red K4 (2, 20, 22, 38) -> avoid connecting to all four
ccx q[1],q[19],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[21],q[43];
ccx q[1],q[19],q[42];

// red K4 (2, 20, 32, 33) -> avoid connecting to all four
ccx q[1],q[19],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[31],q[43];
ccx q[1],q[19],q[42];

// red K4 (2, 22, 23, 35) -> avoid connecting to all four
ccx q[1],q[21],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[22],q[43];
ccx q[1],q[21],q[42];

// red K4 (2, 22, 24, 38) -> avoid connecting to all four
ccx q[1],q[21],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[23],q[43];
ccx q[1],q[21],q[42];

// red K4 (2, 23, 25, 35) -> avoid connecting to all four
ccx q[1],q[22],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[24],q[43];
ccx q[1],q[22],q[42];

// red K4 (2, 23, 33, 35) -> avoid connecting to all four
ccx q[1],q[22],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[32],q[43];
ccx q[1],q[22],q[42];

// red K4 (2, 24, 31, 38) -> avoid connecting to all four
ccx q[1],q[23],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[30],q[43];
ccx q[1],q[23],q[42];

// red K4 (2, 31, 32, 33) -> avoid connecting to all four
ccx q[1],q[30],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[31],q[43];
ccx q[1],q[30],q[42];

// red K4 (3, 4, 16, 26) -> avoid connecting to all four
ccx q[2],q[3],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[25],q[44];
rz(0.15) q[44];
ccx q[43],q[25],q[44];
ccx q[42],q[15],q[43];
ccx q[2],q[3],q[42];

// red K4 (3, 4, 16, 34) -> avoid connecting to all four
ccx q[2],q[3],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[15],q[43];
ccx q[2],q[3],q[42];

// red K4 (3, 4, 17, 24) -> avoid connecting to all four
ccx q[2],q[3],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[23],q[44];
rz(0.15) q[44];
ccx q[43],q[23],q[44];
ccx q[42],q[16],q[43];
ccx q[2],q[3],q[42];

// red K4 (3, 4, 17, 33) -> avoid connecting to all four
ccx q[2],q[3],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[16],q[43];
ccx q[2],q[3],q[42];

// red K4 (3, 4, 24, 26) -> avoid connecting to all four
ccx q[2],q[3],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[25],q[44];
rz(0.15) q[44];
ccx q[43],q[25],q[44];
ccx q[42],q[23],q[43];
ccx q[2],q[3],q[42];

// red K4 (3, 4, 24, 34) -> avoid connecting to all four
ccx q[2],q[3],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[23],q[43];
ccx q[2],q[3],q[42];

// red K4 (3, 4, 25, 26) -> avoid connecting to all four
ccx q[2],q[3],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[25],q[44];
rz(0.15) q[44];
ccx q[43],q[25],q[44];
ccx q[42],q[24],q[43];
ccx q[2],q[3],q[42];

// red K4 (3, 4, 26, 33) -> avoid connecting to all four
ccx q[2],q[3],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[25],q[43];
ccx q[2],q[3],q[42];

// red K4 (3, 5, 15, 17) -> avoid connecting to all four
ccx q[2],q[4],q[42];
ccx q[42],q[14],q[43];
ccx q[43],q[16],q[44];
rz(0.15) q[44];
ccx q[43],q[16],q[44];
ccx q[42],q[14],q[43];
ccx q[2],q[4],q[42];

// red K4 (3, 5, 15, 25) -> avoid connecting to all four
ccx q[2],q[4],q[42];
ccx q[42],q[14],q[43];
ccx q[43],q[24],q[44];
rz(0.15) q[44];
ccx q[43],q[24],q[44];
ccx q[42],q[14],q[43];
ccx q[2],q[4],q[42];

// red K4 (3, 5, 15, 28) -> avoid connecting to all four
ccx q[2],q[4],q[42];
ccx q[42],q[14],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[14],q[43];
ccx q[2],q[4],q[42];

// red K4 (3, 5, 15, 36) -> avoid connecting to all four
ccx q[2],q[4],q[42];
ccx q[42],q[14],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[14],q[43];
ccx q[2],q[4],q[42];

// red K4 (3, 5, 17, 19) -> avoid connecting to all four
ccx q[2],q[4],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[18],q[44];
rz(0.15) q[44];
ccx q[43],q[18],q[44];
ccx q[42],q[16],q[43];
ccx q[2],q[4],q[42];

// red K4 (3, 5, 17, 30) -> avoid connecting to all four
ccx q[2],q[4],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[16],q[43];
ccx q[2],q[4],q[42];

// red K4 (3, 5, 19, 21) -> avoid connecting to all four
ccx q[2],q[4],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[20],q[44];
rz(0.15) q[44];
ccx q[43],q[20],q[44];
ccx q[42],q[18],q[43];
ccx q[2],q[4],q[42];

// red K4 (3, 5, 19, 26) -> avoid connecting to all four
ccx q[2],q[4],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[25],q[44];
rz(0.15) q[44];
ccx q[43],q[25],q[44];
ccx q[42],q[18],q[43];
ccx q[2],q[4],q[42];

// red K4 (3, 5, 19, 32) -> avoid connecting to all four
ccx q[2],q[4],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[18],q[43];
ccx q[2],q[4],q[42];

// red K4 (3, 5, 21, 23) -> avoid connecting to all four
ccx q[2],q[4],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[22],q[44];
rz(0.15) q[44];
ccx q[43],q[22],q[44];
ccx q[42],q[20],q[43];
ccx q[2],q[4],q[42];

// red K4 (3, 5, 21, 28) -> avoid connecting to all four
ccx q[2],q[4],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[20],q[43];
ccx q[2],q[4],q[42];

// red K4 (3, 5, 21, 34) -> avoid connecting to all four
ccx q[2],q[4],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[20],q[43];
ccx q[2],q[4],q[42];

// red K4 (3, 5, 23, 25) -> avoid connecting to all four
ccx q[2],q[4],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[24],q[44];
rz(0.15) q[44];
ccx q[43],q[24],q[44];
ccx q[42],q[22],q[43];
ccx q[2],q[4],q[42];

// red K4 (3, 5, 23, 30) -> avoid connecting to all four
ccx q[2],q[4],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[22],q[43];
ccx q[2],q[4],q[42];

// red K4 (3, 5, 23, 36) -> avoid connecting to all four
ccx q[2],q[4],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[22],q[43];
ccx q[2],q[4],q[42];

// red K4 (3, 5, 25, 26) -> avoid connecting to all four
ccx q[2],q[4],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[25],q[44];
rz(0.15) q[44];
ccx q[43],q[25],q[44];
ccx q[42],q[24],q[43];
ccx q[2],q[4],q[42];

// red K4 (3, 5, 25, 32) -> avoid connecting to all four
ccx q[2],q[4],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[24],q[43];
ccx q[2],q[4],q[42];

// red K4 (3, 5, 26, 28) -> avoid connecting to all four
ccx q[2],q[4],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[25],q[43];
ccx q[2],q[4],q[42];

// red K4 (3, 5, 26, 36) -> avoid connecting to all four
ccx q[2],q[4],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[25],q[43];
ccx q[2],q[4],q[42];

// red K4 (3, 5, 28, 30) -> avoid connecting to all four
ccx q[2],q[4],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[27],q[43];
ccx q[2],q[4],q[42];

// red K4 (3, 5, 30, 32) -> avoid connecting to all four
ccx q[2],q[4],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[29],q[43];
ccx q[2],q[4],q[42];

// red K4 (3, 5, 32, 34) -> avoid connecting to all four
ccx q[2],q[4],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[31],q[43];
ccx q[2],q[4],q[42];

// red K4 (3, 5, 34, 36) -> avoid connecting to all four
ccx q[2],q[4],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[33],q[43];
ccx q[2],q[4],q[42];

// red K4 (3, 10, 17, 24) -> avoid connecting to all four
ccx q[2],q[9],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[23],q[44];
rz(0.15) q[44];
ccx q[43],q[23],q[44];
ccx q[42],q[16],q[43];
ccx q[2],q[9],q[42];

// red K4 (3, 10, 17, 30) -> avoid connecting to all four
ccx q[2],q[9],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[16],q[43];
ccx q[2],q[9],q[42];

// red K4 (3, 10, 17, 33) -> avoid connecting to all four
ccx q[2],q[9],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[16],q[43];
ccx q[2],q[9],q[42];

// red K4 (3, 10, 17, 39) -> avoid connecting to all four
ccx q[2],q[9],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[16],q[43];
ccx q[2],q[9],q[42];

// red K4 (3, 10, 23, 30) -> avoid connecting to all four
ccx q[2],q[9],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[22],q[43];
ccx q[2],q[9],q[42];

// red K4 (3, 10, 23, 33) -> avoid connecting to all four
ccx q[2],q[9],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[22],q[43];
ccx q[2],q[9],q[42];

// red K4 (3, 10, 23, 39) -> avoid connecting to all four
ccx q[2],q[9],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[22],q[43];
ccx q[2],q[9],q[42];

// red K4 (3, 10, 24, 26) -> avoid connecting to all four
ccx q[2],q[9],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[25],q[44];
rz(0.15) q[44];
ccx q[43],q[25],q[44];
ccx q[42],q[23],q[43];
ccx q[2],q[9],q[42];

// red K4 (3, 10, 26, 28) -> avoid connecting to all four
ccx q[2],q[9],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[25],q[43];
ccx q[2],q[9],q[42];

// red K4 (3, 10, 26, 33) -> avoid connecting to all four
ccx q[2],q[9],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[25],q[43];
ccx q[2],q[9],q[42];

// red K4 (3, 10, 26, 39) -> avoid connecting to all four
ccx q[2],q[9],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[25],q[43];
ccx q[2],q[9],q[42];

// red K4 (3, 10, 28, 30) -> avoid connecting to all four
ccx q[2],q[9],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[27],q[43];
ccx q[2],q[9],q[42];

// red K4 (3, 10, 30, 32) -> avoid connecting to all four
ccx q[2],q[9],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[29],q[43];
ccx q[2],q[9],q[42];

// red K4 (3, 10, 32, 33) -> avoid connecting to all four
ccx q[2],q[9],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[31],q[43];
ccx q[2],q[9],q[42];

// red K4 (3, 10, 32, 39) -> avoid connecting to all four
ccx q[2],q[9],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[31],q[43];
ccx q[2],q[9],q[42];

// red K4 (3, 13, 15, 25) -> avoid connecting to all four
ccx q[2],q[12],q[42];
ccx q[42],q[14],q[43];
ccx q[43],q[24],q[44];
rz(0.15) q[44];
ccx q[43],q[24],q[44];
ccx q[42],q[14],q[43];
ccx q[2],q[12],q[42];

// red K4 (3, 13, 15, 33) -> avoid connecting to all four
ccx q[2],q[12],q[42];
ccx q[42],q[14],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[14],q[43];
ccx q[2],q[12],q[42];

// red K4 (3, 13, 15, 36) -> avoid connecting to all four
ccx q[2],q[12],q[42];
ccx q[42],q[14],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[14],q[43];
ccx q[2],q[12],q[42];

// red K4 (3, 13, 23, 25) -> avoid connecting to all four
ccx q[2],q[12],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[24],q[44];
rz(0.15) q[44];
ccx q[43],q[24],q[44];
ccx q[42],q[22],q[43];
ccx q[2],q[12],q[42];

// red K4 (3, 13, 23, 33) -> avoid connecting to all four
ccx q[2],q[12],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[22],q[43];
ccx q[2],q[12],q[42];

// red K4 (3, 13, 23, 36) -> avoid connecting to all four
ccx q[2],q[12],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[22],q[43];
ccx q[2],q[12],q[42];

// red K4 (3, 13, 25, 26) -> avoid connecting to all four
ccx q[2],q[12],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[25],q[44];
rz(0.15) q[44];
ccx q[43],q[25],q[44];
ccx q[42],q[24],q[43];
ccx q[2],q[12],q[42];

// red K4 (3, 13, 26, 33) -> avoid connecting to all four
ccx q[2],q[12],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[25],q[43];
ccx q[2],q[12],q[42];

// red K4 (3, 13, 26, 36) -> avoid connecting to all four
ccx q[2],q[12],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[25],q[43];
ccx q[2],q[12],q[42];

// red K4 (3, 13, 34, 36) -> avoid connecting to all four
ccx q[2],q[12],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[33],q[43];
ccx q[2],q[12],q[42];

// red K4 (3, 15, 17, 33) -> avoid connecting to all four
ccx q[2],q[14],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[16],q[43];
ccx q[2],q[14],q[42];

// red K4 (3, 16, 23, 30) -> avoid connecting to all four
ccx q[2],q[15],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[22],q[43];
ccx q[2],q[15],q[42];

// red K4 (3, 16, 23, 36) -> avoid connecting to all four
ccx q[2],q[15],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[22],q[43];
ccx q[2],q[15],q[42];

// red K4 (3, 16, 23, 39) -> avoid connecting to all four
ccx q[2],q[15],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[22],q[43];
ccx q[2],q[15],q[42];

// red K4 (3, 16, 26, 28) -> avoid connecting to all four
ccx q[2],q[15],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[25],q[43];
ccx q[2],q[15],q[42];

// red K4 (3, 16, 26, 36) -> avoid connecting to all four
ccx q[2],q[15],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[25],q[43];
ccx q[2],q[15],q[42];

// red K4 (3, 16, 26, 39) -> avoid connecting to all four
ccx q[2],q[15],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[25],q[43];
ccx q[2],q[15],q[42];

// red K4 (3, 16, 28, 30) -> avoid connecting to all four
ccx q[2],q[15],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[27],q[43];
ccx q[2],q[15],q[42];

// red K4 (3, 16, 30, 32) -> avoid connecting to all four
ccx q[2],q[15],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[29],q[43];
ccx q[2],q[15],q[42];

// red K4 (3, 16, 32, 34) -> avoid connecting to all four
ccx q[2],q[15],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[31],q[43];
ccx q[2],q[15],q[42];

// red K4 (3, 16, 32, 39) -> avoid connecting to all four
ccx q[2],q[15],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[31],q[43];
ccx q[2],q[15],q[42];

// red K4 (3, 16, 34, 36) -> avoid connecting to all four
ccx q[2],q[15],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[33],q[43];
ccx q[2],q[15],q[42];

// red K4 (3, 17, 19, 33) -> avoid connecting to all four
ccx q[2],q[16],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[18],q[43];
ccx q[2],q[16],q[42];

// red K4 (3, 17, 19, 39) -> avoid connecting to all four
ccx q[2],q[16],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[18],q[43];
ccx q[2],q[16],q[42];

// red K4 (3, 19, 21, 33) -> avoid connecting to all four
ccx q[2],q[18],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[20],q[43];
ccx q[2],q[18],q[42];

// red K4 (3, 19, 21, 39) -> avoid connecting to all four
ccx q[2],q[18],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[20],q[43];
ccx q[2],q[18],q[42];

// red K4 (3, 19, 26, 33) -> avoid connecting to all four
ccx q[2],q[18],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[25],q[43];
ccx q[2],q[18],q[42];

// red K4 (3, 19, 26, 39) -> avoid connecting to all four
ccx q[2],q[18],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[25],q[43];
ccx q[2],q[18],q[42];

// red K4 (3, 19, 32, 33) -> avoid connecting to all four
ccx q[2],q[18],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[31],q[43];
ccx q[2],q[18],q[42];

// red K4 (3, 19, 32, 39) -> avoid connecting to all four
ccx q[2],q[18],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[31],q[43];
ccx q[2],q[18],q[42];

// red K4 (3, 21, 23, 33) -> avoid connecting to all four
ccx q[2],q[20],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[22],q[43];
ccx q[2],q[20],q[42];

// red K4 (3, 21, 23, 39) -> avoid connecting to all four
ccx q[2],q[20],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[22],q[43];
ccx q[2],q[20],q[42];

// red K4 (3, 23, 25, 39) -> avoid connecting to all four
ccx q[2],q[22],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[24],q[43];
ccx q[2],q[22],q[42];

// red K4 (3, 24, 26, 36) -> avoid connecting to all four
ccx q[2],q[23],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[25],q[43];
ccx q[2],q[23],q[42];

// red K4 (3, 24, 34, 36) -> avoid connecting to all four
ccx q[2],q[23],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[33],q[43];
ccx q[2],q[23],q[42];

// red K4 (3, 25, 26, 39) -> avoid connecting to all four
ccx q[2],q[24],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[25],q[43];
ccx q[2],q[24],q[42];

// red K4 (3, 25, 32, 39) -> avoid connecting to all four
ccx q[2],q[24],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[31],q[43];
ccx q[2],q[24],q[42];

// red K4 (4, 6, 16, 18) -> avoid connecting to all four
ccx q[3],q[5],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[17],q[44];
rz(0.15) q[44];
ccx q[43],q[17],q[44];
ccx q[42],q[15],q[43];
ccx q[3],q[5],q[42];

// red K4 (4, 6, 16, 26) -> avoid connecting to all four
ccx q[3],q[5],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[25],q[44];
rz(0.15) q[44];
ccx q[43],q[25],q[44];
ccx q[42],q[15],q[43];
ccx q[3],q[5],q[42];

// red K4 (4, 6, 16, 29) -> avoid connecting to all four
ccx q[3],q[5],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[15],q[43];
ccx q[3],q[5],q[42];

// red K4 (4, 6, 16, 37) -> avoid connecting to all four
ccx q[3],q[5],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[15],q[43];
ccx q[3],q[5],q[42];

// red K4 (4, 6, 18, 20) -> avoid connecting to all four
ccx q[3],q[5],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[19],q[44];
rz(0.15) q[44];
ccx q[43],q[19],q[44];
ccx q[42],q[17],q[43];
ccx q[3],q[5],q[42];

// red K4 (4, 6, 18, 31) -> avoid connecting to all four
ccx q[3],q[5],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[17],q[43];
ccx q[3],q[5],q[42];

// red K4 (4, 6, 20, 22) -> avoid connecting to all four
ccx q[3],q[5],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[21],q[44];
rz(0.15) q[44];
ccx q[43],q[21],q[44];
ccx q[42],q[19],q[43];
ccx q[3],q[5],q[42];

// red K4 (4, 6, 20, 27) -> avoid connecting to all four
ccx q[3],q[5],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[26],q[44];
rz(0.15) q[44];
ccx q[43],q[26],q[44];
ccx q[42],q[19],q[43];
ccx q[3],q[5],q[42];

// red K4 (4, 6, 20, 33) -> avoid connecting to all four
ccx q[3],q[5],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[19],q[43];
ccx q[3],q[5],q[42];

// red K4 (4, 6, 22, 24) -> avoid connecting to all four
ccx q[3],q[5],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[23],q[44];
rz(0.15) q[44];
ccx q[43],q[23],q[44];
ccx q[42],q[21],q[43];
ccx q[3],q[5],q[42];

// red K4 (4, 6, 22, 29) -> avoid connecting to all four
ccx q[3],q[5],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[21],q[43];
ccx q[3],q[5],q[42];

// red K4 (4, 6, 22, 35) -> avoid connecting to all four
ccx q[3],q[5],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[21],q[43];
ccx q[3],q[5],q[42];

// red K4 (4, 6, 24, 26) -> avoid connecting to all four
ccx q[3],q[5],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[25],q[44];
rz(0.15) q[44];
ccx q[43],q[25],q[44];
ccx q[42],q[23],q[43];
ccx q[3],q[5],q[42];

// red K4 (4, 6, 24, 31) -> avoid connecting to all four
ccx q[3],q[5],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[23],q[43];
ccx q[3],q[5],q[42];

// red K4 (4, 6, 24, 37) -> avoid connecting to all four
ccx q[3],q[5],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[23],q[43];
ccx q[3],q[5],q[42];

// red K4 (4, 6, 26, 27) -> avoid connecting to all four
ccx q[3],q[5],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[26],q[44];
rz(0.15) q[44];
ccx q[43],q[26],q[44];
ccx q[42],q[25],q[43];
ccx q[3],q[5],q[42];

// red K4 (4, 6, 26, 33) -> avoid connecting to all four
ccx q[3],q[5],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[25],q[43];
ccx q[3],q[5],q[42];

// red K4 (4, 6, 27, 29) -> avoid connecting to all four
ccx q[3],q[5],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[26],q[43];
ccx q[3],q[5],q[42];

// red K4 (4, 6, 27, 37) -> avoid connecting to all four
ccx q[3],q[5],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[26],q[43];
ccx q[3],q[5],q[42];

// red K4 (4, 6, 29, 31) -> avoid connecting to all four
ccx q[3],q[5],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[28],q[43];
ccx q[3],q[5],q[42];

// red K4 (4, 6, 31, 33) -> avoid connecting to all four
ccx q[3],q[5],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[30],q[43];
ccx q[3],q[5],q[42];

// red K4 (4, 6, 33, 35) -> avoid connecting to all four
ccx q[3],q[5],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[32],q[43];
ccx q[3],q[5],q[42];

// red K4 (4, 6, 35, 37) -> avoid connecting to all four
ccx q[3],q[5],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[34],q[43];
ccx q[3],q[5],q[42];

// red K4 (4, 11, 18, 25) -> avoid connecting to all four
ccx q[3],q[10],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[24],q[44];
rz(0.15) q[44];
ccx q[43],q[24],q[44];
ccx q[42],q[17],q[43];
ccx q[3],q[10],q[42];

// red K4 (4, 11, 18, 31) -> avoid connecting to all four
ccx q[3],q[10],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[17],q[43];
ccx q[3],q[10],q[42];

// red K4 (4, 11, 18, 34) -> avoid connecting to all four
ccx q[3],q[10],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[17],q[43];
ccx q[3],q[10],q[42];

// red K4 (4, 11, 18, 40) -> avoid connecting to all four
ccx q[3],q[10],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[17],q[43];
ccx q[3],q[10],q[42];

// red K4 (4, 11, 24, 31) -> avoid connecting to all four
ccx q[3],q[10],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[23],q[43];
ccx q[3],q[10],q[42];

// red K4 (4, 11, 24, 34) -> avoid connecting to all four
ccx q[3],q[10],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[23],q[43];
ccx q[3],q[10],q[42];

// red K4 (4, 11, 24, 40) -> avoid connecting to all four
ccx q[3],q[10],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[23],q[43];
ccx q[3],q[10],q[42];

// red K4 (4, 11, 25, 27) -> avoid connecting to all four
ccx q[3],q[10],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[26],q[44];
rz(0.15) q[44];
ccx q[43],q[26],q[44];
ccx q[42],q[24],q[43];
ccx q[3],q[10],q[42];

// red K4 (4, 11, 27, 29) -> avoid connecting to all four
ccx q[3],q[10],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[26],q[43];
ccx q[3],q[10],q[42];

// red K4 (4, 11, 27, 34) -> avoid connecting to all four
ccx q[3],q[10],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[26],q[43];
ccx q[3],q[10],q[42];

// red K4 (4, 11, 27, 40) -> avoid connecting to all four
ccx q[3],q[10],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[26],q[43];
ccx q[3],q[10],q[42];

// red K4 (4, 11, 29, 31) -> avoid connecting to all four
ccx q[3],q[10],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[28],q[43];
ccx q[3],q[10],q[42];

// red K4 (4, 11, 31, 33) -> avoid connecting to all four
ccx q[3],q[10],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[30],q[43];
ccx q[3],q[10],q[42];

// red K4 (4, 11, 33, 40) -> avoid connecting to all four
ccx q[3],q[10],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[32],q[43];
ccx q[3],q[10],q[42];

// red K4 (4, 14, 16, 26) -> avoid connecting to all four
ccx q[3],q[13],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[25],q[44];
rz(0.15) q[44];
ccx q[43],q[25],q[44];
ccx q[42],q[15],q[43];
ccx q[3],q[13],q[42];

// red K4 (4, 14, 16, 34) -> avoid connecting to all four
ccx q[3],q[13],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[15],q[43];
ccx q[3],q[13],q[42];

// red K4 (4, 14, 16, 37) -> avoid connecting to all four
ccx q[3],q[13],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[15],q[43];
ccx q[3],q[13],q[42];

// red K4 (4, 14, 24, 26) -> avoid connecting to all four
ccx q[3],q[13],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[25],q[44];
rz(0.15) q[44];
ccx q[43],q[25],q[44];
ccx q[42],q[23],q[43];
ccx q[3],q[13],q[42];

// red K4 (4, 14, 24, 34) -> avoid connecting to all four
ccx q[3],q[13],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[23],q[43];
ccx q[3],q[13],q[42];

// red K4 (4, 14, 24, 37) -> avoid connecting to all four
ccx q[3],q[13],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[23],q[43];
ccx q[3],q[13],q[42];

// red K4 (4, 14, 26, 27) -> avoid connecting to all four
ccx q[3],q[13],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[26],q[44];
rz(0.15) q[44];
ccx q[43],q[26],q[44];
ccx q[42],q[25],q[43];
ccx q[3],q[13],q[42];

// red K4 (4, 14, 27, 34) -> avoid connecting to all four
ccx q[3],q[13],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[26],q[43];
ccx q[3],q[13],q[42];

// red K4 (4, 14, 27, 37) -> avoid connecting to all four
ccx q[3],q[13],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[26],q[43];
ccx q[3],q[13],q[42];

// red K4 (4, 14, 34, 35) -> avoid connecting to all four
ccx q[3],q[13],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[33],q[43];
ccx q[3],q[13],q[42];

// red K4 (4, 14, 35, 37) -> avoid connecting to all four
ccx q[3],q[13],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[34],q[43];
ccx q[3],q[13],q[42];

// red K4 (4, 16, 18, 34) -> avoid connecting to all four
ccx q[3],q[15],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[17],q[43];
ccx q[3],q[15],q[42];

// red K4 (4, 17, 18, 31) -> avoid connecting to all four
ccx q[3],q[16],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[17],q[43];
ccx q[3],q[16],q[42];

// red K4 (4, 17, 18, 40) -> avoid connecting to all four
ccx q[3],q[16],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[17],q[43];
ccx q[3],q[16],q[42];

// red K4 (4, 17, 24, 31) -> avoid connecting to all four
ccx q[3],q[16],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[23],q[43];
ccx q[3],q[16],q[42];

// red K4 (4, 17, 24, 37) -> avoid connecting to all four
ccx q[3],q[16],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[23],q[43];
ccx q[3],q[16],q[42];

// red K4 (4, 17, 24, 40) -> avoid connecting to all four
ccx q[3],q[16],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[23],q[43];
ccx q[3],q[16],q[42];

// red K4 (4, 17, 27, 29) -> avoid connecting to all four
ccx q[3],q[16],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[26],q[43];
ccx q[3],q[16],q[42];

// red K4 (4, 17, 27, 37) -> avoid connecting to all four
ccx q[3],q[16],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[26],q[43];
ccx q[3],q[16],q[42];

// red K4 (4, 17, 27, 40) -> avoid connecting to all four
ccx q[3],q[16],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[26],q[43];
ccx q[3],q[16],q[42];

// red K4 (4, 17, 29, 31) -> avoid connecting to all four
ccx q[3],q[16],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[28],q[43];
ccx q[3],q[16],q[42];

// red K4 (4, 17, 31, 33) -> avoid connecting to all four
ccx q[3],q[16],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[30],q[43];
ccx q[3],q[16],q[42];

// red K4 (4, 17, 33, 35) -> avoid connecting to all four
ccx q[3],q[16],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[32],q[43];
ccx q[3],q[16],q[42];

// red K4 (4, 17, 33, 40) -> avoid connecting to all four
ccx q[3],q[16],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[32],q[43];
ccx q[3],q[16],q[42];

// red K4 (4, 17, 35, 37) -> avoid connecting to all four
ccx q[3],q[16],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[34],q[43];
ccx q[3],q[16],q[42];

// red K4 (4, 18, 20, 34) -> avoid connecting to all four
ccx q[3],q[17],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[19],q[43];
ccx q[3],q[17],q[42];

// red K4 (4, 18, 20, 40) -> avoid connecting to all four
ccx q[3],q[17],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[19],q[43];
ccx q[3],q[17],q[42];

// red K4 (4, 20, 22, 34) -> avoid connecting to all four
ccx q[3],q[19],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[21],q[43];
ccx q[3],q[19],q[42];

// red K4 (4, 20, 22, 40) -> avoid connecting to all four
ccx q[3],q[19],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[21],q[43];
ccx q[3],q[19],q[42];

// red K4 (4, 20, 27, 34) -> avoid connecting to all four
ccx q[3],q[19],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[26],q[43];
ccx q[3],q[19],q[42];

// red K4 (4, 20, 27, 40) -> avoid connecting to all four
ccx q[3],q[19],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[26],q[43];
ccx q[3],q[19],q[42];

// red K4 (4, 20, 33, 40) -> avoid connecting to all four
ccx q[3],q[19],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[32],q[43];
ccx q[3],q[19],q[42];

// red K4 (4, 22, 24, 34) -> avoid connecting to all four
ccx q[3],q[21],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[23],q[43];
ccx q[3],q[21],q[42];

// red K4 (4, 22, 24, 40) -> avoid connecting to all four
ccx q[3],q[21],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[23],q[43];
ccx q[3],q[21],q[42];

// red K4 (4, 22, 34, 35) -> avoid connecting to all four
ccx q[3],q[21],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[33],q[43];
ccx q[3],q[21],q[42];

// red K4 (4, 24, 26, 40) -> avoid connecting to all four
ccx q[3],q[23],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[25],q[43];
ccx q[3],q[23],q[42];

// red K4 (4, 25, 26, 27) -> avoid connecting to all four
ccx q[3],q[24],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[26],q[44];
rz(0.15) q[44];
ccx q[43],q[26],q[44];
ccx q[42],q[25],q[43];
ccx q[3],q[24],q[42];

// red K4 (4, 25, 27, 37) -> avoid connecting to all four
ccx q[3],q[24],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[26],q[43];
ccx q[3],q[24],q[42];

// red K4 (4, 25, 35, 37) -> avoid connecting to all four
ccx q[3],q[24],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[34],q[43];
ccx q[3],q[24],q[42];

// red K4 (4, 26, 27, 40) -> avoid connecting to all four
ccx q[3],q[25],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[26],q[43];
ccx q[3],q[25],q[42];

// red K4 (4, 26, 33, 40) -> avoid connecting to all four
ccx q[3],q[25],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[32],q[43];
ccx q[3],q[25],q[42];

// red K4 (5, 7, 17, 19) -> avoid connecting to all four
ccx q[4],q[6],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[18],q[44];
rz(0.15) q[44];
ccx q[43],q[18],q[44];
ccx q[42],q[16],q[43];
ccx q[4],q[6],q[42];

// red K4 (5, 7, 17, 27) -> avoid connecting to all four
ccx q[4],q[6],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[26],q[44];
rz(0.15) q[44];
ccx q[43],q[26],q[44];
ccx q[42],q[16],q[43];
ccx q[4],q[6],q[42];

// red K4 (5, 7, 17, 30) -> avoid connecting to all four
ccx q[4],q[6],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[16],q[43];
ccx q[4],q[6],q[42];

// red K4 (5, 7, 17, 38) -> avoid connecting to all four
ccx q[4],q[6],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[16],q[43];
ccx q[4],q[6],q[42];

// red K4 (5, 7, 19, 21) -> avoid connecting to all four
ccx q[4],q[6],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[20],q[44];
rz(0.15) q[44];
ccx q[43],q[20],q[44];
ccx q[42],q[18],q[43];
ccx q[4],q[6],q[42];

// red K4 (5, 7, 19, 32) -> avoid connecting to all four
ccx q[4],q[6],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[18],q[43];
ccx q[4],q[6],q[42];

// red K4 (5, 7, 21, 23) -> avoid connecting to all four
ccx q[4],q[6],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[22],q[44];
rz(0.15) q[44];
ccx q[43],q[22],q[44];
ccx q[42],q[20],q[43];
ccx q[4],q[6],q[42];

// red K4 (5, 7, 21, 28) -> avoid connecting to all four
ccx q[4],q[6],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[20],q[43];
ccx q[4],q[6],q[42];

// red K4 (5, 7, 21, 34) -> avoid connecting to all four
ccx q[4],q[6],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[20],q[43];
ccx q[4],q[6],q[42];

// red K4 (5, 7, 23, 25) -> avoid connecting to all four
ccx q[4],q[6],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[24],q[44];
rz(0.15) q[44];
ccx q[43],q[24],q[44];
ccx q[42],q[22],q[43];
ccx q[4],q[6],q[42];

// red K4 (5, 7, 23, 30) -> avoid connecting to all four
ccx q[4],q[6],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[22],q[43];
ccx q[4],q[6],q[42];

// red K4 (5, 7, 23, 36) -> avoid connecting to all four
ccx q[4],q[6],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[22],q[43];
ccx q[4],q[6],q[42];

// red K4 (5, 7, 25, 27) -> avoid connecting to all four
ccx q[4],q[6],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[26],q[44];
rz(0.15) q[44];
ccx q[43],q[26],q[44];
ccx q[42],q[24],q[43];
ccx q[4],q[6],q[42];

// red K4 (5, 7, 25, 32) -> avoid connecting to all four
ccx q[4],q[6],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[24],q[43];
ccx q[4],q[6],q[42];

// red K4 (5, 7, 25, 38) -> avoid connecting to all four
ccx q[4],q[6],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[24],q[43];
ccx q[4],q[6],q[42];

// red K4 (5, 7, 27, 28) -> avoid connecting to all four
ccx q[4],q[6],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[26],q[43];
ccx q[4],q[6],q[42];

// red K4 (5, 7, 27, 34) -> avoid connecting to all four
ccx q[4],q[6],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[26],q[43];
ccx q[4],q[6],q[42];

// red K4 (5, 7, 28, 30) -> avoid connecting to all four
ccx q[4],q[6],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[27],q[43];
ccx q[4],q[6],q[42];

// red K4 (5, 7, 28, 38) -> avoid connecting to all four
ccx q[4],q[6],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[27],q[43];
ccx q[4],q[6],q[42];

// red K4 (5, 7, 30, 32) -> avoid connecting to all four
ccx q[4],q[6],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[29],q[43];
ccx q[4],q[6],q[42];

// red K4 (5, 7, 32, 34) -> avoid connecting to all four
ccx q[4],q[6],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[31],q[43];
ccx q[4],q[6],q[42];

// red K4 (5, 7, 34, 36) -> avoid connecting to all four
ccx q[4],q[6],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[33],q[43];
ccx q[4],q[6],q[42];

// red K4 (5, 7, 36, 38) -> avoid connecting to all four
ccx q[4],q[6],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[35],q[43];
ccx q[4],q[6],q[42];

// red K4 (5, 12, 19, 26) -> avoid connecting to all four
ccx q[4],q[11],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[25],q[44];
rz(0.15) q[44];
ccx q[43],q[25],q[44];
ccx q[42],q[18],q[43];
ccx q[4],q[11],q[42];

// red K4 (5, 12, 19, 32) -> avoid connecting to all four
ccx q[4],q[11],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[18],q[43];
ccx q[4],q[11],q[42];

// red K4 (5, 12, 19, 35) -> avoid connecting to all four
ccx q[4],q[11],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[18],q[43];
ccx q[4],q[11],q[42];

// red K4 (5, 12, 19, 41) -> avoid connecting to all four
ccx q[4],q[11],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[18],q[43];
ccx q[4],q[11],q[42];

// red K4 (5, 12, 25, 26) -> avoid connecting to all four
ccx q[4],q[11],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[25],q[44];
rz(0.15) q[44];
ccx q[43],q[25],q[44];
ccx q[42],q[24],q[43];
ccx q[4],q[11],q[42];

// red K4 (5, 12, 25, 32) -> avoid connecting to all four
ccx q[4],q[11],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[24],q[43];
ccx q[4],q[11],q[42];

// red K4 (5, 12, 25, 35) -> avoid connecting to all four
ccx q[4],q[11],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[24],q[43];
ccx q[4],q[11],q[42];

// red K4 (5, 12, 25, 41) -> avoid connecting to all four
ccx q[4],q[11],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[24],q[43];
ccx q[4],q[11],q[42];

// red K4 (5, 12, 26, 28) -> avoid connecting to all four
ccx q[4],q[11],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[25],q[43];
ccx q[4],q[11],q[42];

// red K4 (5, 12, 28, 30) -> avoid connecting to all four
ccx q[4],q[11],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[27],q[43];
ccx q[4],q[11],q[42];

// red K4 (5, 12, 28, 35) -> avoid connecting to all four
ccx q[4],q[11],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[27],q[43];
ccx q[4],q[11],q[42];

// red K4 (5, 12, 28, 41) -> avoid connecting to all four
ccx q[4],q[11],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[27],q[43];
ccx q[4],q[11],q[42];

// red K4 (5, 12, 30, 32) -> avoid connecting to all four
ccx q[4],q[11],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[29],q[43];
ccx q[4],q[11],q[42];

// red K4 (5, 12, 32, 34) -> avoid connecting to all four
ccx q[4],q[11],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[31],q[43];
ccx q[4],q[11],q[42];

// red K4 (5, 12, 34, 35) -> avoid connecting to all four
ccx q[4],q[11],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[33],q[43];
ccx q[4],q[11],q[42];

// red K4 (5, 12, 34, 41) -> avoid connecting to all four
ccx q[4],q[11],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[33],q[43];
ccx q[4],q[11],q[42];

// red K4 (5, 15, 17, 27) -> avoid connecting to all four
ccx q[4],q[14],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[26],q[44];
rz(0.15) q[44];
ccx q[43],q[26],q[44];
ccx q[42],q[16],q[43];
ccx q[4],q[14],q[42];

// red K4 (5, 15, 17, 35) -> avoid connecting to all four
ccx q[4],q[14],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[16],q[43];
ccx q[4],q[14],q[42];

// red K4 (5, 15, 17, 38) -> avoid connecting to all four
ccx q[4],q[14],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[16],q[43];
ccx q[4],q[14],q[42];

// red K4 (5, 15, 25, 27) -> avoid connecting to all four
ccx q[4],q[14],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[26],q[44];
rz(0.15) q[44];
ccx q[43],q[26],q[44];
ccx q[42],q[24],q[43];
ccx q[4],q[14],q[42];

// red K4 (5, 15, 25, 35) -> avoid connecting to all four
ccx q[4],q[14],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[24],q[43];
ccx q[4],q[14],q[42];

// red K4 (5, 15, 25, 38) -> avoid connecting to all four
ccx q[4],q[14],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[24],q[43];
ccx q[4],q[14],q[42];

// red K4 (5, 15, 27, 28) -> avoid connecting to all four
ccx q[4],q[14],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[26],q[43];
ccx q[4],q[14],q[42];

// red K4 (5, 15, 28, 35) -> avoid connecting to all four
ccx q[4],q[14],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[27],q[43];
ccx q[4],q[14],q[42];

// red K4 (5, 15, 28, 38) -> avoid connecting to all four
ccx q[4],q[14],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[27],q[43];
ccx q[4],q[14],q[42];

// red K4 (5, 15, 35, 36) -> avoid connecting to all four
ccx q[4],q[14],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[34],q[43];
ccx q[4],q[14],q[42];

// red K4 (5, 15, 36, 38) -> avoid connecting to all four
ccx q[4],q[14],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[35],q[43];
ccx q[4],q[14],q[42];

// red K4 (5, 17, 18, 19) -> avoid connecting to all four
ccx q[4],q[16],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[18],q[44];
rz(0.15) q[44];
ccx q[43],q[18],q[44];
ccx q[42],q[17],q[43];
ccx q[4],q[16],q[42];

// red K4 (5, 17, 18, 30) -> avoid connecting to all four
ccx q[4],q[16],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[17],q[43];
ccx q[4],q[16],q[42];

// red K4 (5, 17, 18, 38) -> avoid connecting to all four
ccx q[4],q[16],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[17],q[43];
ccx q[4],q[16],q[42];

// red K4 (5, 17, 19, 35) -> avoid connecting to all four
ccx q[4],q[16],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[18],q[43];
ccx q[4],q[16],q[42];

// red K4 (5, 18, 19, 32) -> avoid connecting to all four
ccx q[4],q[17],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[18],q[43];
ccx q[4],q[17],q[42];

// red K4 (5, 18, 19, 41) -> avoid connecting to all four
ccx q[4],q[17],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[18],q[43];
ccx q[4],q[17],q[42];

// red K4 (5, 18, 25, 32) -> avoid connecting to all four
ccx q[4],q[17],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[24],q[43];
ccx q[4],q[17],q[42];

// red K4 (5, 18, 25, 38) -> avoid connecting to all four
ccx q[4],q[17],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[24],q[43];
ccx q[4],q[17],q[42];

// red K4 (5, 18, 25, 41) -> avoid connecting to all four
ccx q[4],q[17],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[24],q[43];
ccx q[4],q[17],q[42];

// red K4 (5, 18, 28, 30) -> avoid connecting to all four
ccx q[4],q[17],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[27],q[43];
ccx q[4],q[17],q[42];

// red K4 (5, 18, 28, 38) -> avoid connecting to all four
ccx q[4],q[17],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[27],q[43];
ccx q[4],q[17],q[42];

// red K4 (5, 18, 28, 41) -> avoid connecting to all four
ccx q[4],q[17],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[27],q[43];
ccx q[4],q[17],q[42];

// red K4 (5, 18, 30, 32) -> avoid connecting to all four
ccx q[4],q[17],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[29],q[43];
ccx q[4],q[17],q[42];

// red K4 (5, 18, 32, 34) -> avoid connecting to all four
ccx q[4],q[17],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[31],q[43];
ccx q[4],q[17],q[42];

// red K4 (5, 18, 34, 36) -> avoid connecting to all four
ccx q[4],q[17],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[33],q[43];
ccx q[4],q[17],q[42];

// red K4 (5, 18, 34, 41) -> avoid connecting to all four
ccx q[4],q[17],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[33],q[43];
ccx q[4],q[17],q[42];

// red K4 (5, 18, 36, 38) -> avoid connecting to all four
ccx q[4],q[17],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[35],q[43];
ccx q[4],q[17],q[42];

// red K4 (5, 19, 21, 35) -> avoid connecting to all four
ccx q[4],q[18],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[20],q[43];
ccx q[4],q[18],q[42];

// red K4 (5, 19, 21, 41) -> avoid connecting to all four
ccx q[4],q[18],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[20],q[43];
ccx q[4],q[18],q[42];

// red K4 (5, 21, 23, 35) -> avoid connecting to all four
ccx q[4],q[20],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[22],q[43];
ccx q[4],q[20],q[42];

// red K4 (5, 21, 23, 41) -> avoid connecting to all four
ccx q[4],q[20],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[22],q[43];
ccx q[4],q[20],q[42];

// red K4 (5, 21, 28, 35) -> avoid connecting to all four
ccx q[4],q[20],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[27],q[43];
ccx q[4],q[20],q[42];

// red K4 (5, 21, 28, 41) -> avoid connecting to all four
ccx q[4],q[20],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[27],q[43];
ccx q[4],q[20],q[42];

// red K4 (5, 21, 34, 35) -> avoid connecting to all four
ccx q[4],q[20],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[33],q[43];
ccx q[4],q[20],q[42];

// red K4 (5, 21, 34, 41) -> avoid connecting to all four
ccx q[4],q[20],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[33],q[43];
ccx q[4],q[20],q[42];

// red K4 (5, 23, 25, 35) -> avoid connecting to all four
ccx q[4],q[22],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[24],q[43];
ccx q[4],q[22],q[42];

// red K4 (5, 23, 25, 41) -> avoid connecting to all four
ccx q[4],q[22],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[24],q[43];
ccx q[4],q[22],q[42];

// red K4 (5, 23, 35, 36) -> avoid connecting to all four
ccx q[4],q[22],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[34],q[43];
ccx q[4],q[22],q[42];

// red K4 (5, 25, 26, 27) -> avoid connecting to all four
ccx q[4],q[24],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[26],q[44];
rz(0.15) q[44];
ccx q[43],q[26],q[44];
ccx q[42],q[25],q[43];
ccx q[4],q[24],q[42];

// red K4 (5, 25, 26, 38) -> avoid connecting to all four
ccx q[4],q[24],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[25],q[43];
ccx q[4],q[24],q[42];

// red K4 (5, 25, 27, 41) -> avoid connecting to all four
ccx q[4],q[24],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[26],q[43];
ccx q[4],q[24],q[42];

// red K4 (5, 26, 27, 28) -> avoid connecting to all four
ccx q[4],q[25],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[26],q[43];
ccx q[4],q[25],q[42];

// red K4 (5, 26, 28, 38) -> avoid connecting to all four
ccx q[4],q[25],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[27],q[43];
ccx q[4],q[25],q[42];

// red K4 (5, 26, 36, 38) -> avoid connecting to all four
ccx q[4],q[25],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[35],q[43];
ccx q[4],q[25],q[42];

// red K4 (5, 27, 28, 41) -> avoid connecting to all four
ccx q[4],q[26],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[27],q[43];
ccx q[4],q[26],q[42];

// red K4 (5, 27, 34, 41) -> avoid connecting to all four
ccx q[4],q[26],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[33],q[43];
ccx q[4],q[26],q[42];

// red K4 (5, 34, 35, 36) -> avoid connecting to all four
ccx q[4],q[33],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[34],q[43];
ccx q[4],q[33],q[42];

// red K4 (6, 8, 18, 20) -> avoid connecting to all four
ccx q[5],q[7],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[19],q[44];
rz(0.15) q[44];
ccx q[43],q[19],q[44];
ccx q[42],q[17],q[43];
ccx q[5],q[7],q[42];

// red K4 (6, 8, 18, 28) -> avoid connecting to all four
ccx q[5],q[7],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[17],q[43];
ccx q[5],q[7],q[42];

// red K4 (6, 8, 18, 31) -> avoid connecting to all four
ccx q[5],q[7],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[17],q[43];
ccx q[5],q[7],q[42];

// red K4 (6, 8, 18, 39) -> avoid connecting to all four
ccx q[5],q[7],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[17],q[43];
ccx q[5],q[7],q[42];

// red K4 (6, 8, 20, 22) -> avoid connecting to all four
ccx q[5],q[7],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[21],q[44];
rz(0.15) q[44];
ccx q[43],q[21],q[44];
ccx q[42],q[19],q[43];
ccx q[5],q[7],q[42];

// red K4 (6, 8, 20, 33) -> avoid connecting to all four
ccx q[5],q[7],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[19],q[43];
ccx q[5],q[7],q[42];

// red K4 (6, 8, 22, 24) -> avoid connecting to all four
ccx q[5],q[7],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[23],q[44];
rz(0.15) q[44];
ccx q[43],q[23],q[44];
ccx q[42],q[21],q[43];
ccx q[5],q[7],q[42];

// red K4 (6, 8, 22, 29) -> avoid connecting to all four
ccx q[5],q[7],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[21],q[43];
ccx q[5],q[7],q[42];

// red K4 (6, 8, 22, 35) -> avoid connecting to all four
ccx q[5],q[7],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[21],q[43];
ccx q[5],q[7],q[42];

// red K4 (6, 8, 24, 26) -> avoid connecting to all four
ccx q[5],q[7],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[25],q[44];
rz(0.15) q[44];
ccx q[43],q[25],q[44];
ccx q[42],q[23],q[43];
ccx q[5],q[7],q[42];

// red K4 (6, 8, 24, 31) -> avoid connecting to all four
ccx q[5],q[7],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[23],q[43];
ccx q[5],q[7],q[42];

// red K4 (6, 8, 24, 37) -> avoid connecting to all four
ccx q[5],q[7],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[23],q[43];
ccx q[5],q[7],q[42];

// red K4 (6, 8, 26, 28) -> avoid connecting to all four
ccx q[5],q[7],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[25],q[43];
ccx q[5],q[7],q[42];

// red K4 (6, 8, 26, 33) -> avoid connecting to all four
ccx q[5],q[7],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[25],q[43];
ccx q[5],q[7],q[42];

// red K4 (6, 8, 26, 39) -> avoid connecting to all four
ccx q[5],q[7],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[25],q[43];
ccx q[5],q[7],q[42];

// red K4 (6, 8, 28, 29) -> avoid connecting to all four
ccx q[5],q[7],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[27],q[43];
ccx q[5],q[7],q[42];

// red K4 (6, 8, 28, 35) -> avoid connecting to all four
ccx q[5],q[7],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[27],q[43];
ccx q[5],q[7],q[42];

// red K4 (6, 8, 29, 31) -> avoid connecting to all four
ccx q[5],q[7],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[28],q[43];
ccx q[5],q[7],q[42];

// red K4 (6, 8, 29, 39) -> avoid connecting to all four
ccx q[5],q[7],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[28],q[43];
ccx q[5],q[7],q[42];

// red K4 (6, 8, 31, 33) -> avoid connecting to all four
ccx q[5],q[7],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[30],q[43];
ccx q[5],q[7],q[42];

// red K4 (6, 8, 33, 35) -> avoid connecting to all four
ccx q[5],q[7],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[32],q[43];
ccx q[5],q[7],q[42];

// red K4 (6, 8, 35, 37) -> avoid connecting to all four
ccx q[5],q[7],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[34],q[43];
ccx q[5],q[7],q[42];

// red K4 (6, 8, 37, 39) -> avoid connecting to all four
ccx q[5],q[7],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[36],q[43];
ccx q[5],q[7],q[42];

// red K4 (6, 13, 20, 27) -> avoid connecting to all four
ccx q[5],q[12],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[26],q[44];
rz(0.15) q[44];
ccx q[43],q[26],q[44];
ccx q[42],q[19],q[43];
ccx q[5],q[12],q[42];

// red K4 (6, 13, 20, 33) -> avoid connecting to all four
ccx q[5],q[12],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[19],q[43];
ccx q[5],q[12],q[42];

// red K4 (6, 13, 20, 36) -> avoid connecting to all four
ccx q[5],q[12],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[19],q[43];
ccx q[5],q[12],q[42];

// red K4 (6, 13, 20, 42) -> avoid connecting to all four
ccx q[5],q[12],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[19],q[43];
ccx q[5],q[12],q[42];

// red K4 (6, 13, 26, 27) -> avoid connecting to all four
ccx q[5],q[12],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[26],q[44];
rz(0.15) q[44];
ccx q[43],q[26],q[44];
ccx q[42],q[25],q[43];
ccx q[5],q[12],q[42];

// red K4 (6, 13, 26, 33) -> avoid connecting to all four
ccx q[5],q[12],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[25],q[43];
ccx q[5],q[12],q[42];

// red K4 (6, 13, 26, 36) -> avoid connecting to all four
ccx q[5],q[12],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[25],q[43];
ccx q[5],q[12],q[42];

// red K4 (6, 13, 26, 42) -> avoid connecting to all four
ccx q[5],q[12],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[25],q[43];
ccx q[5],q[12],q[42];

// red K4 (6, 13, 27, 29) -> avoid connecting to all four
ccx q[5],q[12],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[26],q[43];
ccx q[5],q[12],q[42];

// red K4 (6, 13, 29, 31) -> avoid connecting to all four
ccx q[5],q[12],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[28],q[43];
ccx q[5],q[12],q[42];

// red K4 (6, 13, 29, 36) -> avoid connecting to all four
ccx q[5],q[12],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[28],q[43];
ccx q[5],q[12],q[42];

// red K4 (6, 13, 29, 42) -> avoid connecting to all four
ccx q[5],q[12],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[28],q[43];
ccx q[5],q[12],q[42];

// red K4 (6, 13, 31, 33) -> avoid connecting to all four
ccx q[5],q[12],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[30],q[43];
ccx q[5],q[12],q[42];

// red K4 (6, 13, 33, 35) -> avoid connecting to all four
ccx q[5],q[12],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[32],q[43];
ccx q[5],q[12],q[42];

// red K4 (6, 13, 35, 36) -> avoid connecting to all four
ccx q[5],q[12],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[34],q[43];
ccx q[5],q[12],q[42];

// red K4 (6, 13, 35, 42) -> avoid connecting to all four
ccx q[5],q[12],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[34],q[43];
ccx q[5],q[12],q[42];

// red K4 (6, 16, 18, 28) -> avoid connecting to all four
ccx q[5],q[15],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[17],q[43];
ccx q[5],q[15],q[42];

// red K4 (6, 16, 18, 36) -> avoid connecting to all four
ccx q[5],q[15],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[17],q[43];
ccx q[5],q[15],q[42];

// red K4 (6, 16, 18, 39) -> avoid connecting to all four
ccx q[5],q[15],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[17],q[43];
ccx q[5],q[15],q[42];

// red K4 (6, 16, 26, 28) -> avoid connecting to all four
ccx q[5],q[15],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[25],q[43];
ccx q[5],q[15],q[42];

// red K4 (6, 16, 26, 36) -> avoid connecting to all four
ccx q[5],q[15],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[25],q[43];
ccx q[5],q[15],q[42];

// red K4 (6, 16, 26, 39) -> avoid connecting to all four
ccx q[5],q[15],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[25],q[43];
ccx q[5],q[15],q[42];

// red K4 (6, 16, 28, 29) -> avoid connecting to all four
ccx q[5],q[15],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[27],q[43];
ccx q[5],q[15],q[42];

// red K4 (6, 16, 29, 36) -> avoid connecting to all four
ccx q[5],q[15],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[28],q[43];
ccx q[5],q[15],q[42];

// red K4 (6, 16, 29, 39) -> avoid connecting to all four
ccx q[5],q[15],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[28],q[43];
ccx q[5],q[15],q[42];

// red K4 (6, 16, 36, 37) -> avoid connecting to all four
ccx q[5],q[15],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[35],q[43];
ccx q[5],q[15],q[42];

// red K4 (6, 16, 37, 39) -> avoid connecting to all four
ccx q[5],q[15],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[36],q[43];
ccx q[5],q[15],q[42];

// red K4 (6, 18, 19, 20) -> avoid connecting to all four
ccx q[5],q[17],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[19],q[44];
rz(0.15) q[44];
ccx q[43],q[19],q[44];
ccx q[42],q[18],q[43];
ccx q[5],q[17],q[42];

// red K4 (6, 18, 19, 31) -> avoid connecting to all four
ccx q[5],q[17],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[18],q[43];
ccx q[5],q[17],q[42];

// red K4 (6, 18, 19, 39) -> avoid connecting to all four
ccx q[5],q[17],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[18],q[43];
ccx q[5],q[17],q[42];

// red K4 (6, 18, 20, 36) -> avoid connecting to all four
ccx q[5],q[17],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[19],q[43];
ccx q[5],q[17],q[42];

// red K4 (6, 19, 20, 33) -> avoid connecting to all four
ccx q[5],q[18],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[19],q[43];
ccx q[5],q[18],q[42];

// red K4 (6, 19, 20, 42) -> avoid connecting to all four
ccx q[5],q[18],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[19],q[43];
ccx q[5],q[18],q[42];

// red K4 (6, 19, 26, 33) -> avoid connecting to all four
ccx q[5],q[18],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[25],q[43];
ccx q[5],q[18],q[42];

// red K4 (6, 19, 26, 39) -> avoid connecting to all four
ccx q[5],q[18],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[25],q[43];
ccx q[5],q[18],q[42];

// red K4 (6, 19, 26, 42) -> avoid connecting to all four
ccx q[5],q[18],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[25],q[43];
ccx q[5],q[18],q[42];

// red K4 (6, 19, 29, 31) -> avoid connecting to all four
ccx q[5],q[18],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[28],q[43];
ccx q[5],q[18],q[42];

// red K4 (6, 19, 29, 39) -> avoid connecting to all four
ccx q[5],q[18],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[28],q[43];
ccx q[5],q[18],q[42];

// red K4 (6, 19, 29, 42) -> avoid connecting to all four
ccx q[5],q[18],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[28],q[43];
ccx q[5],q[18],q[42];

// red K4 (6, 19, 31, 33) -> avoid connecting to all four
ccx q[5],q[18],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[30],q[43];
ccx q[5],q[18],q[42];

// red K4 (6, 19, 33, 35) -> avoid connecting to all four
ccx q[5],q[18],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[32],q[43];
ccx q[5],q[18],q[42];

// red K4 (6, 19, 35, 37) -> avoid connecting to all four
ccx q[5],q[18],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[34],q[43];
ccx q[5],q[18],q[42];

// red K4 (6, 19, 35, 42) -> avoid connecting to all four
ccx q[5],q[18],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[34],q[43];
ccx q[5],q[18],q[42];

// red K4 (6, 19, 37, 39) -> avoid connecting to all four
ccx q[5],q[18],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[36],q[43];
ccx q[5],q[18],q[42];

// red K4 (6, 20, 22, 36) -> avoid connecting to all four
ccx q[5],q[19],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[21],q[43];
ccx q[5],q[19],q[42];

// red K4 (6, 20, 22, 42) -> avoid connecting to all four
ccx q[5],q[19],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[21],q[43];
ccx q[5],q[19],q[42];

// red K4 (6, 22, 24, 36) -> avoid connecting to all four
ccx q[5],q[21],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[23],q[43];
ccx q[5],q[21],q[42];

// red K4 (6, 22, 24, 42) -> avoid connecting to all four
ccx q[5],q[21],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[23],q[43];
ccx q[5],q[21],q[42];

// red K4 (6, 22, 29, 36) -> avoid connecting to all four
ccx q[5],q[21],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[28],q[43];
ccx q[5],q[21],q[42];

// red K4 (6, 22, 29, 42) -> avoid connecting to all four
ccx q[5],q[21],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[28],q[43];
ccx q[5],q[21],q[42];

// red K4 (6, 22, 35, 36) -> avoid connecting to all four
ccx q[5],q[21],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[34],q[43];
ccx q[5],q[21],q[42];

// red K4 (6, 22, 35, 42) -> avoid connecting to all four
ccx q[5],q[21],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[34],q[43];
ccx q[5],q[21],q[42];

// red K4 (6, 24, 26, 36) -> avoid connecting to all four
ccx q[5],q[23],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[25],q[43];
ccx q[5],q[23],q[42];

// red K4 (6, 24, 26, 42) -> avoid connecting to all four
ccx q[5],q[23],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[25],q[43];
ccx q[5],q[23],q[42];

// red K4 (6, 24, 36, 37) -> avoid connecting to all four
ccx q[5],q[23],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[35],q[43];
ccx q[5],q[23],q[42];

// red K4 (6, 26, 27, 28) -> avoid connecting to all four
ccx q[5],q[25],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[26],q[43];
ccx q[5],q[25],q[42];

// red K4 (6, 26, 27, 39) -> avoid connecting to all four
ccx q[5],q[25],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[26],q[43];
ccx q[5],q[25],q[42];

// red K4 (6, 26, 28, 42) -> avoid connecting to all four
ccx q[5],q[25],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[27],q[43];
ccx q[5],q[25],q[42];

// red K4 (6, 27, 28, 29) -> avoid connecting to all four
ccx q[5],q[26],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[27],q[43];
ccx q[5],q[26],q[42];

// red K4 (6, 27, 29, 39) -> avoid connecting to all four
ccx q[5],q[26],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[28],q[43];
ccx q[5],q[26],q[42];

// red K4 (6, 27, 37, 39) -> avoid connecting to all four
ccx q[5],q[26],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[36],q[43];
ccx q[5],q[26],q[42];

// red K4 (6, 28, 29, 42) -> avoid connecting to all four
ccx q[5],q[27],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[28],q[43];
ccx q[5],q[27],q[42];

// red K4 (6, 28, 35, 42) -> avoid connecting to all four
ccx q[5],q[27],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[34],q[43];
ccx q[5],q[27],q[42];

// red K4 (6, 35, 36, 37) -> avoid connecting to all four
ccx q[5],q[34],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[35],q[43];
ccx q[5],q[34],q[42];

// red K4 (7, 9, 19, 21) -> avoid connecting to all four
ccx q[6],q[8],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[20],q[44];
rz(0.15) q[44];
ccx q[43],q[20],q[44];
ccx q[42],q[18],q[43];
ccx q[6],q[8],q[42];

// red K4 (7, 9, 19, 29) -> avoid connecting to all four
ccx q[6],q[8],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[18],q[43];
ccx q[6],q[8],q[42];

// red K4 (7, 9, 19, 32) -> avoid connecting to all four
ccx q[6],q[8],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[18],q[43];
ccx q[6],q[8],q[42];

// red K4 (7, 9, 19, 40) -> avoid connecting to all four
ccx q[6],q[8],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[18],q[43];
ccx q[6],q[8],q[42];

// red K4 (7, 9, 21, 23) -> avoid connecting to all four
ccx q[6],q[8],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[22],q[44];
rz(0.15) q[44];
ccx q[43],q[22],q[44];
ccx q[42],q[20],q[43];
ccx q[6],q[8],q[42];

// red K4 (7, 9, 21, 34) -> avoid connecting to all four
ccx q[6],q[8],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[20],q[43];
ccx q[6],q[8],q[42];

// red K4 (7, 9, 23, 25) -> avoid connecting to all four
ccx q[6],q[8],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[24],q[44];
rz(0.15) q[44];
ccx q[43],q[24],q[44];
ccx q[42],q[22],q[43];
ccx q[6],q[8],q[42];

// red K4 (7, 9, 23, 30) -> avoid connecting to all four
ccx q[6],q[8],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[22],q[43];
ccx q[6],q[8],q[42];

// red K4 (7, 9, 23, 36) -> avoid connecting to all four
ccx q[6],q[8],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[22],q[43];
ccx q[6],q[8],q[42];

// red K4 (7, 9, 25, 27) -> avoid connecting to all four
ccx q[6],q[8],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[26],q[44];
rz(0.15) q[44];
ccx q[43],q[26],q[44];
ccx q[42],q[24],q[43];
ccx q[6],q[8],q[42];

// red K4 (7, 9, 25, 32) -> avoid connecting to all four
ccx q[6],q[8],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[24],q[43];
ccx q[6],q[8],q[42];

// red K4 (7, 9, 25, 38) -> avoid connecting to all four
ccx q[6],q[8],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[24],q[43];
ccx q[6],q[8],q[42];

// red K4 (7, 9, 27, 29) -> avoid connecting to all four
ccx q[6],q[8],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[26],q[43];
ccx q[6],q[8],q[42];

// red K4 (7, 9, 27, 34) -> avoid connecting to all four
ccx q[6],q[8],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[26],q[43];
ccx q[6],q[8],q[42];

// red K4 (7, 9, 27, 40) -> avoid connecting to all four
ccx q[6],q[8],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[26],q[43];
ccx q[6],q[8],q[42];

// red K4 (7, 9, 29, 30) -> avoid connecting to all four
ccx q[6],q[8],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[28],q[43];
ccx q[6],q[8],q[42];

// red K4 (7, 9, 29, 36) -> avoid connecting to all four
ccx q[6],q[8],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[28],q[43];
ccx q[6],q[8],q[42];

// red K4 (7, 9, 30, 32) -> avoid connecting to all four
ccx q[6],q[8],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[29],q[43];
ccx q[6],q[8],q[42];

// red K4 (7, 9, 30, 40) -> avoid connecting to all four
ccx q[6],q[8],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[29],q[43];
ccx q[6],q[8],q[42];

// red K4 (7, 9, 32, 34) -> avoid connecting to all four
ccx q[6],q[8],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[31],q[43];
ccx q[6],q[8],q[42];

// red K4 (7, 9, 34, 36) -> avoid connecting to all four
ccx q[6],q[8],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[33],q[43];
ccx q[6],q[8],q[42];

// red K4 (7, 9, 36, 38) -> avoid connecting to all four
ccx q[6],q[8],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[35],q[43];
ccx q[6],q[8],q[42];

// red K4 (7, 9, 38, 40) -> avoid connecting to all four
ccx q[6],q[8],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[37],q[43];
ccx q[6],q[8],q[42];

// red K4 (7, 14, 21, 28) -> avoid connecting to all four
ccx q[6],q[13],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[20],q[43];
ccx q[6],q[13],q[42];

// red K4 (7, 14, 21, 34) -> avoid connecting to all four
ccx q[6],q[13],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[20],q[43];
ccx q[6],q[13],q[42];

// red K4 (7, 14, 21, 37) -> avoid connecting to all four
ccx q[6],q[13],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[20],q[43];
ccx q[6],q[13],q[42];

// red K4 (7, 14, 27, 28) -> avoid connecting to all four
ccx q[6],q[13],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[26],q[43];
ccx q[6],q[13],q[42];

// red K4 (7, 14, 27, 34) -> avoid connecting to all four
ccx q[6],q[13],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[26],q[43];
ccx q[6],q[13],q[42];

// red K4 (7, 14, 27, 37) -> avoid connecting to all four
ccx q[6],q[13],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[26],q[43];
ccx q[6],q[13],q[42];

// red K4 (7, 14, 28, 30) -> avoid connecting to all four
ccx q[6],q[13],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[27],q[43];
ccx q[6],q[13],q[42];

// red K4 (7, 14, 30, 32) -> avoid connecting to all four
ccx q[6],q[13],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[29],q[43];
ccx q[6],q[13],q[42];

// red K4 (7, 14, 30, 37) -> avoid connecting to all four
ccx q[6],q[13],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[29],q[43];
ccx q[6],q[13],q[42];

// red K4 (7, 14, 32, 34) -> avoid connecting to all four
ccx q[6],q[13],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[31],q[43];
ccx q[6],q[13],q[42];

// red K4 (7, 14, 34, 36) -> avoid connecting to all four
ccx q[6],q[13],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[33],q[43];
ccx q[6],q[13],q[42];

// red K4 (7, 14, 36, 37) -> avoid connecting to all four
ccx q[6],q[13],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[35],q[43];
ccx q[6],q[13],q[42];

// red K4 (7, 17, 19, 29) -> avoid connecting to all four
ccx q[6],q[16],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[18],q[43];
ccx q[6],q[16],q[42];

// red K4 (7, 17, 19, 37) -> avoid connecting to all four
ccx q[6],q[16],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[18],q[43];
ccx q[6],q[16],q[42];

// red K4 (7, 17, 19, 40) -> avoid connecting to all four
ccx q[6],q[16],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[18],q[43];
ccx q[6],q[16],q[42];

// red K4 (7, 17, 27, 29) -> avoid connecting to all four
ccx q[6],q[16],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[26],q[43];
ccx q[6],q[16],q[42];

// red K4 (7, 17, 27, 37) -> avoid connecting to all four
ccx q[6],q[16],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[26],q[43];
ccx q[6],q[16],q[42];

// red K4 (7, 17, 27, 40) -> avoid connecting to all four
ccx q[6],q[16],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[26],q[43];
ccx q[6],q[16],q[42];

// red K4 (7, 17, 29, 30) -> avoid connecting to all four
ccx q[6],q[16],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[28],q[43];
ccx q[6],q[16],q[42];

// red K4 (7, 17, 30, 37) -> avoid connecting to all four
ccx q[6],q[16],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[29],q[43];
ccx q[6],q[16],q[42];

// red K4 (7, 17, 30, 40) -> avoid connecting to all four
ccx q[6],q[16],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[29],q[43];
ccx q[6],q[16],q[42];

// red K4 (7, 17, 37, 38) -> avoid connecting to all four
ccx q[6],q[16],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[36],q[43];
ccx q[6],q[16],q[42];

// red K4 (7, 17, 38, 40) -> avoid connecting to all four
ccx q[6],q[16],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[37],q[43];
ccx q[6],q[16],q[42];

// red K4 (7, 19, 20, 21) -> avoid connecting to all four
ccx q[6],q[18],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[20],q[44];
rz(0.15) q[44];
ccx q[43],q[20],q[44];
ccx q[42],q[19],q[43];
ccx q[6],q[18],q[42];

// red K4 (7, 19, 20, 32) -> avoid connecting to all four
ccx q[6],q[18],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[19],q[43];
ccx q[6],q[18],q[42];

// red K4 (7, 19, 20, 40) -> avoid connecting to all four
ccx q[6],q[18],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[19],q[43];
ccx q[6],q[18],q[42];

// red K4 (7, 19, 21, 37) -> avoid connecting to all four
ccx q[6],q[18],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[20],q[43];
ccx q[6],q[18],q[42];

// red K4 (7, 20, 21, 34) -> avoid connecting to all four
ccx q[6],q[19],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[20],q[43];
ccx q[6],q[19],q[42];

// red K4 (7, 20, 27, 34) -> avoid connecting to all four
ccx q[6],q[19],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[26],q[43];
ccx q[6],q[19],q[42];

// red K4 (7, 20, 27, 40) -> avoid connecting to all four
ccx q[6],q[19],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[26],q[43];
ccx q[6],q[19],q[42];

// red K4 (7, 20, 30, 32) -> avoid connecting to all four
ccx q[6],q[19],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[29],q[43];
ccx q[6],q[19],q[42];

// red K4 (7, 20, 30, 40) -> avoid connecting to all four
ccx q[6],q[19],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[29],q[43];
ccx q[6],q[19],q[42];

// red K4 (7, 20, 32, 34) -> avoid connecting to all four
ccx q[6],q[19],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[31],q[43];
ccx q[6],q[19],q[42];

// red K4 (7, 20, 34, 36) -> avoid connecting to all four
ccx q[6],q[19],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[33],q[43];
ccx q[6],q[19],q[42];

// red K4 (7, 20, 36, 38) -> avoid connecting to all four
ccx q[6],q[19],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[35],q[43];
ccx q[6],q[19],q[42];

// red K4 (7, 20, 38, 40) -> avoid connecting to all four
ccx q[6],q[19],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[37],q[43];
ccx q[6],q[19],q[42];

// red K4 (7, 21, 23, 37) -> avoid connecting to all four
ccx q[6],q[20],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[22],q[43];
ccx q[6],q[20],q[42];

// red K4 (7, 23, 25, 37) -> avoid connecting to all four
ccx q[6],q[22],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[24],q[43];
ccx q[6],q[22],q[42];

// red K4 (7, 23, 30, 37) -> avoid connecting to all four
ccx q[6],q[22],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[29],q[43];
ccx q[6],q[22],q[42];

// red K4 (7, 23, 36, 37) -> avoid connecting to all four
ccx q[6],q[22],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[35],q[43];
ccx q[6],q[22],q[42];

// red K4 (7, 25, 27, 37) -> avoid connecting to all four
ccx q[6],q[24],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[26],q[43];
ccx q[6],q[24],q[42];

// red K4 (7, 25, 37, 38) -> avoid connecting to all four
ccx q[6],q[24],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[36],q[43];
ccx q[6],q[24],q[42];

// red K4 (7, 27, 28, 29) -> avoid connecting to all four
ccx q[6],q[26],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[27],q[43];
ccx q[6],q[26],q[42];

// red K4 (7, 27, 28, 40) -> avoid connecting to all four
ccx q[6],q[26],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[27],q[43];
ccx q[6],q[26],q[42];

// red K4 (7, 28, 29, 30) -> avoid connecting to all four
ccx q[6],q[27],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[28],q[43];
ccx q[6],q[27],q[42];

// red K4 (7, 28, 30, 40) -> avoid connecting to all four
ccx q[6],q[27],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[29],q[43];
ccx q[6],q[27],q[42];

// red K4 (7, 28, 38, 40) -> avoid connecting to all four
ccx q[6],q[27],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[37],q[43];
ccx q[6],q[27],q[42];

// red K4 (7, 36, 37, 38) -> avoid connecting to all four
ccx q[6],q[35],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[36],q[43];
ccx q[6],q[35],q[42];

// red K4 (8, 9, 10, 22) -> avoid connecting to all four
ccx q[7],q[8],q[42];
ccx q[42],q[9],q[43];
ccx q[43],q[21],q[44];
rz(0.15) q[44];
ccx q[43],q[21],q[44];
ccx q[42],q[9],q[43];
ccx q[7],q[8],q[42];

// red K4 (8, 9, 10, 30) -> avoid connecting to all four
ccx q[7],q[8],q[42];
ccx q[42],q[9],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[9],q[43];
ccx q[7],q[8],q[42];

// red K4 (8, 9, 10, 31) -> avoid connecting to all four
ccx q[7],q[8],q[42];
ccx q[42],q[9],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[9],q[43];
ccx q[7],q[8],q[42];

// red K4 (8, 9, 10, 39) -> avoid connecting to all four
ccx q[7],q[8],q[42];
ccx q[42],q[9],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[9],q[43];
ccx q[7],q[8],q[42];

// red K4 (8, 9, 21, 22) -> avoid connecting to all four
ccx q[7],q[8],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[21],q[44];
rz(0.15) q[44];
ccx q[43],q[21],q[44];
ccx q[42],q[20],q[43];
ccx q[7],q[8],q[42];

// red K4 (8, 9, 21, 31) -> avoid connecting to all four
ccx q[7],q[8],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[20],q[43];
ccx q[7],q[8],q[42];

// red K4 (8, 9, 21, 39) -> avoid connecting to all four
ccx q[7],q[8],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[20],q[43];
ccx q[7],q[8],q[42];

// red K4 (8, 9, 22, 29) -> avoid connecting to all four
ccx q[7],q[8],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[21],q[43];
ccx q[7],q[8],q[42];

// red K4 (8, 9, 22, 38) -> avoid connecting to all four
ccx q[7],q[8],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[21],q[43];
ccx q[7],q[8],q[42];

// red K4 (8, 9, 29, 30) -> avoid connecting to all four
ccx q[7],q[8],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[28],q[43];
ccx q[7],q[8],q[42];

// red K4 (8, 9, 29, 31) -> avoid connecting to all four
ccx q[7],q[8],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[28],q[43];
ccx q[7],q[8],q[42];

// red K4 (8, 9, 29, 39) -> avoid connecting to all four
ccx q[7],q[8],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[28],q[43];
ccx q[7],q[8],q[42];

// red K4 (8, 9, 31, 38) -> avoid connecting to all four
ccx q[7],q[8],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[30],q[43];
ccx q[7],q[8],q[42];

// red K4 (8, 9, 38, 39) -> avoid connecting to all four
ccx q[7],q[8],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[37],q[43];
ccx q[7],q[8],q[42];

// red K4 (8, 10, 20, 22) -> avoid connecting to all four
ccx q[7],q[9],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[21],q[44];
rz(0.15) q[44];
ccx q[43],q[21],q[44];
ccx q[42],q[19],q[43];
ccx q[7],q[9],q[42];

// red K4 (8, 10, 20, 30) -> avoid connecting to all four
ccx q[7],q[9],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[19],q[43];
ccx q[7],q[9],q[42];

// red K4 (8, 10, 20, 33) -> avoid connecting to all four
ccx q[7],q[9],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[19],q[43];
ccx q[7],q[9],q[42];

// red K4 (8, 10, 20, 41) -> avoid connecting to all four
ccx q[7],q[9],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[19],q[43];
ccx q[7],q[9],q[42];

// red K4 (8, 10, 22, 24) -> avoid connecting to all four
ccx q[7],q[9],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[23],q[44];
rz(0.15) q[44];
ccx q[43],q[23],q[44];
ccx q[42],q[21],q[43];
ccx q[7],q[9],q[42];

// red K4 (8, 10, 22, 35) -> avoid connecting to all four
ccx q[7],q[9],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[21],q[43];
ccx q[7],q[9],q[42];

// red K4 (8, 10, 24, 26) -> avoid connecting to all four
ccx q[7],q[9],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[25],q[44];
rz(0.15) q[44];
ccx q[43],q[25],q[44];
ccx q[42],q[23],q[43];
ccx q[7],q[9],q[42];

// red K4 (8, 10, 24, 31) -> avoid connecting to all four
ccx q[7],q[9],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[23],q[43];
ccx q[7],q[9],q[42];

// red K4 (8, 10, 24, 37) -> avoid connecting to all four
ccx q[7],q[9],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[23],q[43];
ccx q[7],q[9],q[42];

// red K4 (8, 10, 26, 28) -> avoid connecting to all four
ccx q[7],q[9],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[25],q[43];
ccx q[7],q[9],q[42];

// red K4 (8, 10, 26, 33) -> avoid connecting to all four
ccx q[7],q[9],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[25],q[43];
ccx q[7],q[9],q[42];

// red K4 (8, 10, 26, 39) -> avoid connecting to all four
ccx q[7],q[9],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[25],q[43];
ccx q[7],q[9],q[42];

// red K4 (8, 10, 28, 30) -> avoid connecting to all four
ccx q[7],q[9],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[27],q[43];
ccx q[7],q[9],q[42];

// red K4 (8, 10, 28, 35) -> avoid connecting to all four
ccx q[7],q[9],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[27],q[43];
ccx q[7],q[9],q[42];

// red K4 (8, 10, 28, 41) -> avoid connecting to all four
ccx q[7],q[9],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[27],q[43];
ccx q[7],q[9],q[42];

// red K4 (8, 10, 30, 37) -> avoid connecting to all four
ccx q[7],q[9],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[29],q[43];
ccx q[7],q[9],q[42];

// red K4 (8, 10, 31, 33) -> avoid connecting to all four
ccx q[7],q[9],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[30],q[43];
ccx q[7],q[9],q[42];

// red K4 (8, 10, 31, 41) -> avoid connecting to all four
ccx q[7],q[9],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[30],q[43];
ccx q[7],q[9],q[42];

// red K4 (8, 10, 33, 35) -> avoid connecting to all four
ccx q[7],q[9],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[32],q[43];
ccx q[7],q[9],q[42];

// red K4 (8, 10, 35, 37) -> avoid connecting to all four
ccx q[7],q[9],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[34],q[43];
ccx q[7],q[9],q[42];

// red K4 (8, 10, 37, 39) -> avoid connecting to all four
ccx q[7],q[9],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[36],q[43];
ccx q[7],q[9],q[42];

// red K4 (8, 10, 39, 41) -> avoid connecting to all four
ccx q[7],q[9],q[42];
ccx q[42],q[38],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[38],q[43];
ccx q[7],q[9],q[42];

// red K4 (8, 15, 22, 29) -> avoid connecting to all four
ccx q[7],q[14],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[21],q[43];
ccx q[7],q[14],q[42];

// red K4 (8, 15, 22, 35) -> avoid connecting to all four
ccx q[7],q[14],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[21],q[43];
ccx q[7],q[14],q[42];

// red K4 (8, 15, 22, 38) -> avoid connecting to all four
ccx q[7],q[14],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[21],q[43];
ccx q[7],q[14],q[42];

// red K4 (8, 15, 28, 29) -> avoid connecting to all four
ccx q[7],q[14],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[27],q[43];
ccx q[7],q[14],q[42];

// red K4 (8, 15, 28, 35) -> avoid connecting to all four
ccx q[7],q[14],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[27],q[43];
ccx q[7],q[14],q[42];

// red K4 (8, 15, 28, 38) -> avoid connecting to all four
ccx q[7],q[14],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[27],q[43];
ccx q[7],q[14],q[42];

// red K4 (8, 15, 29, 31) -> avoid connecting to all four
ccx q[7],q[14],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[28],q[43];
ccx q[7],q[14],q[42];

// red K4 (8, 15, 31, 33) -> avoid connecting to all four
ccx q[7],q[14],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[30],q[43];
ccx q[7],q[14],q[42];

// red K4 (8, 15, 31, 38) -> avoid connecting to all four
ccx q[7],q[14],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[30],q[43];
ccx q[7],q[14],q[42];

// red K4 (8, 15, 33, 35) -> avoid connecting to all four
ccx q[7],q[14],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[32],q[43];
ccx q[7],q[14],q[42];

// red K4 (8, 15, 35, 37) -> avoid connecting to all four
ccx q[7],q[14],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[34],q[43];
ccx q[7],q[14],q[42];

// red K4 (8, 15, 37, 38) -> avoid connecting to all four
ccx q[7],q[14],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[36],q[43];
ccx q[7],q[14],q[42];

// red K4 (8, 18, 20, 30) -> avoid connecting to all four
ccx q[7],q[17],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[19],q[43];
ccx q[7],q[17],q[42];

// red K4 (8, 18, 20, 38) -> avoid connecting to all four
ccx q[7],q[17],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[19],q[43];
ccx q[7],q[17],q[42];

// red K4 (8, 18, 20, 41) -> avoid connecting to all four
ccx q[7],q[17],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[19],q[43];
ccx q[7],q[17],q[42];

// red K4 (8, 18, 28, 30) -> avoid connecting to all four
ccx q[7],q[17],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[27],q[43];
ccx q[7],q[17],q[42];

// red K4 (8, 18, 28, 38) -> avoid connecting to all four
ccx q[7],q[17],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[27],q[43];
ccx q[7],q[17],q[42];

// red K4 (8, 18, 28, 41) -> avoid connecting to all four
ccx q[7],q[17],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[27],q[43];
ccx q[7],q[17],q[42];

// red K4 (8, 18, 31, 38) -> avoid connecting to all four
ccx q[7],q[17],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[30],q[43];
ccx q[7],q[17],q[42];

// red K4 (8, 18, 31, 41) -> avoid connecting to all four
ccx q[7],q[17],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[30],q[43];
ccx q[7],q[17],q[42];

// red K4 (8, 18, 38, 39) -> avoid connecting to all four
ccx q[7],q[17],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[37],q[43];
ccx q[7],q[17],q[42];

// red K4 (8, 18, 39, 41) -> avoid connecting to all four
ccx q[7],q[17],q[42];
ccx q[42],q[38],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[38],q[43];
ccx q[7],q[17],q[42];

// red K4 (8, 20, 21, 22) -> avoid connecting to all four
ccx q[7],q[19],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[21],q[44];
rz(0.15) q[44];
ccx q[43],q[21],q[44];
ccx q[42],q[20],q[43];
ccx q[7],q[19],q[42];

// red K4 (8, 20, 21, 33) -> avoid connecting to all four
ccx q[7],q[19],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[20],q[43];
ccx q[7],q[19],q[42];

// red K4 (8, 20, 21, 41) -> avoid connecting to all four
ccx q[7],q[19],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[20],q[43];
ccx q[7],q[19],q[42];

// red K4 (8, 20, 22, 38) -> avoid connecting to all four
ccx q[7],q[19],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[21],q[43];
ccx q[7],q[19],q[42];

// red K4 (8, 21, 22, 35) -> avoid connecting to all four
ccx q[7],q[20],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[21],q[43];
ccx q[7],q[20],q[42];

// red K4 (8, 21, 28, 35) -> avoid connecting to all four
ccx q[7],q[20],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[27],q[43];
ccx q[7],q[20],q[42];

// red K4 (8, 21, 28, 41) -> avoid connecting to all four
ccx q[7],q[20],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[27],q[43];
ccx q[7],q[20],q[42];

// red K4 (8, 21, 31, 33) -> avoid connecting to all four
ccx q[7],q[20],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[30],q[43];
ccx q[7],q[20],q[42];

// red K4 (8, 21, 31, 41) -> avoid connecting to all four
ccx q[7],q[20],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[30],q[43];
ccx q[7],q[20],q[42];

// red K4 (8, 21, 33, 35) -> avoid connecting to all four
ccx q[7],q[20],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[32],q[43];
ccx q[7],q[20],q[42];

// red K4 (8, 21, 35, 37) -> avoid connecting to all four
ccx q[7],q[20],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[34],q[43];
ccx q[7],q[20],q[42];

// red K4 (8, 21, 37, 39) -> avoid connecting to all four
ccx q[7],q[20],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[36],q[43];
ccx q[7],q[20],q[42];

// red K4 (8, 21, 39, 41) -> avoid connecting to all four
ccx q[7],q[20],q[42];
ccx q[42],q[38],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[38],q[43];
ccx q[7],q[20],q[42];

// red K4 (8, 22, 24, 38) -> avoid connecting to all four
ccx q[7],q[21],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[23],q[43];
ccx q[7],q[21],q[42];

// red K4 (8, 24, 26, 38) -> avoid connecting to all four
ccx q[7],q[23],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[25],q[43];
ccx q[7],q[23],q[42];

// red K4 (8, 24, 31, 38) -> avoid connecting to all four
ccx q[7],q[23],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[30],q[43];
ccx q[7],q[23],q[42];

// red K4 (8, 24, 37, 38) -> avoid connecting to all four
ccx q[7],q[23],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[36],q[43];
ccx q[7],q[23],q[42];

// red K4 (8, 26, 28, 38) -> avoid connecting to all four
ccx q[7],q[25],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[27],q[43];
ccx q[7],q[25],q[42];

// red K4 (8, 26, 38, 39) -> avoid connecting to all four
ccx q[7],q[25],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[37],q[43];
ccx q[7],q[25],q[42];

// red K4 (8, 28, 29, 30) -> avoid connecting to all four
ccx q[7],q[27],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[28],q[43];
ccx q[7],q[27],q[42];

// red K4 (8, 28, 29, 41) -> avoid connecting to all four
ccx q[7],q[27],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[28],q[43];
ccx q[7],q[27],q[42];

// red K4 (8, 29, 31, 41) -> avoid connecting to all four
ccx q[7],q[28],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[30],q[43];
ccx q[7],q[28],q[42];

// red K4 (8, 29, 39, 41) -> avoid connecting to all four
ccx q[7],q[28],q[42];
ccx q[42],q[38],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[38],q[43];
ccx q[7],q[28],q[42];

// red K4 (8, 37, 38, 39) -> avoid connecting to all four
ccx q[7],q[36],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[37],q[43];
ccx q[7],q[36],q[42];

// red K4 (9, 10, 11, 23) -> avoid connecting to all four
ccx q[8],q[9],q[42];
ccx q[42],q[10],q[43];
ccx q[43],q[22],q[44];
rz(0.15) q[44];
ccx q[43],q[22],q[44];
ccx q[42],q[10],q[43];
ccx q[8],q[9],q[42];

// red K4 (9, 10, 11, 31) -> avoid connecting to all four
ccx q[8],q[9],q[42];
ccx q[42],q[10],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[10],q[43];
ccx q[8],q[9],q[42];

// red K4 (9, 10, 11, 40) -> avoid connecting to all four
ccx q[8],q[9],q[42];
ccx q[42],q[10],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[10],q[43];
ccx q[8],q[9],q[42];

// red K4 (9, 10, 22, 23) -> avoid connecting to all four
ccx q[8],q[9],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[22],q[44];
rz(0.15) q[44];
ccx q[43],q[22],q[44];
ccx q[42],q[21],q[43];
ccx q[8],q[9],q[42];

// red K4 (9, 10, 22, 32) -> avoid connecting to all four
ccx q[8],q[9],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[21],q[43];
ccx q[8],q[9],q[42];

// red K4 (9, 10, 22, 40) -> avoid connecting to all four
ccx q[8],q[9],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[21],q[43];
ccx q[8],q[9],q[42];

// red K4 (9, 10, 23, 30) -> avoid connecting to all four
ccx q[8],q[9],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[22],q[43];
ccx q[8],q[9],q[42];

// red K4 (9, 10, 23, 39) -> avoid connecting to all four
ccx q[8],q[9],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[22],q[43];
ccx q[8],q[9],q[42];

// red K4 (9, 10, 30, 32) -> avoid connecting to all four
ccx q[8],q[9],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[29],q[43];
ccx q[8],q[9],q[42];

// red K4 (9, 10, 30, 40) -> avoid connecting to all four
ccx q[8],q[9],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[29],q[43];
ccx q[8],q[9],q[42];

// red K4 (9, 10, 31, 32) -> avoid connecting to all four
ccx q[8],q[9],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[30],q[43];
ccx q[8],q[9],q[42];

// red K4 (9, 10, 32, 39) -> avoid connecting to all four
ccx q[8],q[9],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[31],q[43];
ccx q[8],q[9],q[42];

// red K4 (9, 11, 21, 23) -> avoid connecting to all four
ccx q[8],q[10],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[22],q[44];
rz(0.15) q[44];
ccx q[43],q[22],q[44];
ccx q[42],q[20],q[43];
ccx q[8],q[10],q[42];

// red K4 (9, 11, 21, 31) -> avoid connecting to all four
ccx q[8],q[10],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[20],q[43];
ccx q[8],q[10],q[42];

// red K4 (9, 11, 21, 34) -> avoid connecting to all four
ccx q[8],q[10],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[20],q[43];
ccx q[8],q[10],q[42];

// red K4 (9, 11, 21, 42) -> avoid connecting to all four
ccx q[8],q[10],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[20],q[43];
ccx q[8],q[10],q[42];

// red K4 (9, 11, 23, 25) -> avoid connecting to all four
ccx q[8],q[10],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[24],q[44];
rz(0.15) q[44];
ccx q[43],q[24],q[44];
ccx q[42],q[22],q[43];
ccx q[8],q[10],q[42];

// red K4 (9, 11, 23, 36) -> avoid connecting to all four
ccx q[8],q[10],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[22],q[43];
ccx q[8],q[10],q[42];

// red K4 (9, 11, 25, 27) -> avoid connecting to all four
ccx q[8],q[10],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[26],q[44];
rz(0.15) q[44];
ccx q[43],q[26],q[44];
ccx q[42],q[24],q[43];
ccx q[8],q[10],q[42];

// red K4 (9, 11, 25, 38) -> avoid connecting to all four
ccx q[8],q[10],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[24],q[43];
ccx q[8],q[10],q[42];

// red K4 (9, 11, 27, 29) -> avoid connecting to all four
ccx q[8],q[10],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[26],q[43];
ccx q[8],q[10],q[42];

// red K4 (9, 11, 27, 34) -> avoid connecting to all four
ccx q[8],q[10],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[26],q[43];
ccx q[8],q[10],q[42];

// red K4 (9, 11, 27, 40) -> avoid connecting to all four
ccx q[8],q[10],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[26],q[43];
ccx q[8],q[10],q[42];

// red K4 (9, 11, 29, 31) -> avoid connecting to all four
ccx q[8],q[10],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[28],q[43];
ccx q[8],q[10],q[42];

// red K4 (9, 11, 29, 36) -> avoid connecting to all four
ccx q[8],q[10],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[28],q[43];
ccx q[8],q[10],q[42];

// red K4 (9, 11, 29, 42) -> avoid connecting to all four
ccx q[8],q[10],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[28],q[43];
ccx q[8],q[10],q[42];

// red K4 (9, 11, 31, 38) -> avoid connecting to all four
ccx q[8],q[10],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[30],q[43];
ccx q[8],q[10],q[42];

// red K4 (9, 11, 34, 36) -> avoid connecting to all four
ccx q[8],q[10],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[33],q[43];
ccx q[8],q[10],q[42];

// red K4 (9, 11, 36, 38) -> avoid connecting to all four
ccx q[8],q[10],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[35],q[43];
ccx q[8],q[10],q[42];

// red K4 (9, 11, 38, 40) -> avoid connecting to all four
ccx q[8],q[10],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[37],q[43];
ccx q[8],q[10],q[42];

// red K4 (9, 11, 40, 42) -> avoid connecting to all four
ccx q[8],q[10],q[42];
ccx q[42],q[39],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[39],q[43];
ccx q[8],q[10],q[42];

// red K4 (9, 16, 23, 30) -> avoid connecting to all four
ccx q[8],q[15],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[22],q[43];
ccx q[8],q[15],q[42];

// red K4 (9, 16, 23, 36) -> avoid connecting to all four
ccx q[8],q[15],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[22],q[43];
ccx q[8],q[15],q[42];

// red K4 (9, 16, 23, 39) -> avoid connecting to all four
ccx q[8],q[15],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[22],q[43];
ccx q[8],q[15],q[42];

// red K4 (9, 16, 29, 30) -> avoid connecting to all four
ccx q[8],q[15],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[28],q[43];
ccx q[8],q[15],q[42];

// red K4 (9, 16, 29, 36) -> avoid connecting to all four
ccx q[8],q[15],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[28],q[43];
ccx q[8],q[15],q[42];

// red K4 (9, 16, 29, 39) -> avoid connecting to all four
ccx q[8],q[15],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[28],q[43];
ccx q[8],q[15],q[42];

// red K4 (9, 16, 30, 32) -> avoid connecting to all four
ccx q[8],q[15],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[29],q[43];
ccx q[8],q[15],q[42];

// red K4 (9, 16, 32, 34) -> avoid connecting to all four
ccx q[8],q[15],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[31],q[43];
ccx q[8],q[15],q[42];

// red K4 (9, 16, 32, 39) -> avoid connecting to all four
ccx q[8],q[15],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[31],q[43];
ccx q[8],q[15],q[42];

// red K4 (9, 16, 34, 36) -> avoid connecting to all four
ccx q[8],q[15],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[33],q[43];
ccx q[8],q[15],q[42];

// red K4 (9, 16, 36, 38) -> avoid connecting to all four
ccx q[8],q[15],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[35],q[43];
ccx q[8],q[15],q[42];

// red K4 (9, 16, 38, 39) -> avoid connecting to all four
ccx q[8],q[15],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[37],q[43];
ccx q[8],q[15],q[42];

// red K4 (9, 19, 21, 31) -> avoid connecting to all four
ccx q[8],q[18],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[20],q[43];
ccx q[8],q[18],q[42];

// red K4 (9, 19, 21, 39) -> avoid connecting to all four
ccx q[8],q[18],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[20],q[43];
ccx q[8],q[18],q[42];

// red K4 (9, 19, 21, 42) -> avoid connecting to all four
ccx q[8],q[18],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[20],q[43];
ccx q[8],q[18],q[42];

// red K4 (9, 19, 29, 31) -> avoid connecting to all four
ccx q[8],q[18],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[28],q[43];
ccx q[8],q[18],q[42];

// red K4 (9, 19, 29, 39) -> avoid connecting to all four
ccx q[8],q[18],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[28],q[43];
ccx q[8],q[18],q[42];

// red K4 (9, 19, 29, 42) -> avoid connecting to all four
ccx q[8],q[18],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[28],q[43];
ccx q[8],q[18],q[42];

// red K4 (9, 19, 31, 32) -> avoid connecting to all four
ccx q[8],q[18],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[30],q[43];
ccx q[8],q[18],q[42];

// red K4 (9, 19, 32, 39) -> avoid connecting to all four
ccx q[8],q[18],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[31],q[43];
ccx q[8],q[18],q[42];

// red K4 (9, 19, 32, 42) -> avoid connecting to all four
ccx q[8],q[18],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[31],q[43];
ccx q[8],q[18],q[42];

// red K4 (9, 19, 40, 42) -> avoid connecting to all four
ccx q[8],q[18],q[42];
ccx q[42],q[39],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[39],q[43];
ccx q[8],q[18],q[42];

// red K4 (9, 21, 22, 23) -> avoid connecting to all four
ccx q[8],q[20],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[22],q[44];
rz(0.15) q[44];
ccx q[43],q[22],q[44];
ccx q[42],q[21],q[43];
ccx q[8],q[20],q[42];

// red K4 (9, 21, 22, 34) -> avoid connecting to all four
ccx q[8],q[20],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[21],q[43];
ccx q[8],q[20],q[42];

// red K4 (9, 21, 22, 42) -> avoid connecting to all four
ccx q[8],q[20],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[21],q[43];
ccx q[8],q[20],q[42];

// red K4 (9, 21, 23, 39) -> avoid connecting to all four
ccx q[8],q[20],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[22],q[43];
ccx q[8],q[20],q[42];

// red K4 (9, 22, 23, 36) -> avoid connecting to all four
ccx q[8],q[21],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[22],q[43];
ccx q[8],q[21],q[42];

// red K4 (9, 22, 29, 36) -> avoid connecting to all four
ccx q[8],q[21],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[28],q[43];
ccx q[8],q[21],q[42];

// red K4 (9, 22, 29, 42) -> avoid connecting to all four
ccx q[8],q[21],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[28],q[43];
ccx q[8],q[21],q[42];

// red K4 (9, 22, 32, 34) -> avoid connecting to all four
ccx q[8],q[21],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[31],q[43];
ccx q[8],q[21],q[42];

// red K4 (9, 22, 32, 42) -> avoid connecting to all four
ccx q[8],q[21],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[31],q[43];
ccx q[8],q[21],q[42];

// red K4 (9, 22, 34, 36) -> avoid connecting to all four
ccx q[8],q[21],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[33],q[43];
ccx q[8],q[21],q[42];

// red K4 (9, 22, 36, 38) -> avoid connecting to all four
ccx q[8],q[21],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[35],q[43];
ccx q[8],q[21],q[42];

// red K4 (9, 22, 38, 40) -> avoid connecting to all four
ccx q[8],q[21],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[37],q[43];
ccx q[8],q[21],q[42];

// red K4 (9, 22, 40, 42) -> avoid connecting to all four
ccx q[8],q[21],q[42];
ccx q[42],q[39],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[39],q[43];
ccx q[8],q[21],q[42];

// red K4 (9, 23, 25, 39) -> avoid connecting to all four
ccx q[8],q[22],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[24],q[43];
ccx q[8],q[22],q[42];

// red K4 (9, 25, 27, 39) -> avoid connecting to all four
ccx q[8],q[24],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[26],q[43];
ccx q[8],q[24],q[42];

// red K4 (9, 25, 32, 39) -> avoid connecting to all four
ccx q[8],q[24],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[31],q[43];
ccx q[8],q[24],q[42];

// red K4 (9, 25, 38, 39) -> avoid connecting to all four
ccx q[8],q[24],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[37],q[43];
ccx q[8],q[24],q[42];

// red K4 (9, 27, 29, 39) -> avoid connecting to all four
ccx q[8],q[26],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[28],q[43];
ccx q[8],q[26],q[42];

// red K4 (9, 29, 30, 42) -> avoid connecting to all four
ccx q[8],q[28],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[29],q[43];
ccx q[8],q[28],q[42];

// red K4 (9, 30, 32, 42) -> avoid connecting to all four
ccx q[8],q[29],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[31],q[43];
ccx q[8],q[29],q[42];

// red K4 (9, 30, 40, 42) -> avoid connecting to all four
ccx q[8],q[29],q[42];
ccx q[42],q[39],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[39],q[43];
ccx q[8],q[29],q[42];

// red K4 (10, 11, 12, 24) -> avoid connecting to all four
ccx q[9],q[10],q[42];
ccx q[42],q[11],q[43];
ccx q[43],q[23],q[44];
rz(0.15) q[44];
ccx q[43],q[23],q[44];
ccx q[42],q[11],q[43];
ccx q[9],q[10],q[42];

// red K4 (10, 11, 12, 33) -> avoid connecting to all four
ccx q[9],q[10],q[42];
ccx q[42],q[11],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[11],q[43];
ccx q[9],q[10],q[42];

// red K4 (10, 11, 12, 41) -> avoid connecting to all four
ccx q[9],q[10],q[42];
ccx q[42],q[11],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[11],q[43];
ccx q[9],q[10],q[42];

// red K4 (10, 11, 23, 33) -> avoid connecting to all four
ccx q[9],q[10],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[22],q[43];
ccx q[9],q[10],q[42];

// red K4 (10, 11, 23, 41) -> avoid connecting to all four
ccx q[9],q[10],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[22],q[43];
ccx q[9],q[10],q[42];

// red K4 (10, 11, 24, 31) -> avoid connecting to all four
ccx q[9],q[10],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[23],q[43];
ccx q[9],q[10],q[42];

// red K4 (10, 11, 24, 40) -> avoid connecting to all four
ccx q[9],q[10],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[23],q[43];
ccx q[9],q[10],q[42];

// red K4 (10, 11, 31, 33) -> avoid connecting to all four
ccx q[9],q[10],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[30],q[43];
ccx q[9],q[10],q[42];

// red K4 (10, 11, 31, 41) -> avoid connecting to all four
ccx q[9],q[10],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[30],q[43];
ccx q[9],q[10],q[42];

// red K4 (10, 11, 33, 40) -> avoid connecting to all four
ccx q[9],q[10],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[32],q[43];
ccx q[9],q[10],q[42];

// red K4 (10, 12, 22, 24) -> avoid connecting to all four
ccx q[9],q[11],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[23],q[44];
rz(0.15) q[44];
ccx q[43],q[23],q[44];
ccx q[42],q[21],q[43];
ccx q[9],q[11],q[42];

// red K4 (10, 12, 22, 32) -> avoid connecting to all four
ccx q[9],q[11],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[21],q[43];
ccx q[9],q[11],q[42];

// red K4 (10, 12, 22, 35) -> avoid connecting to all four
ccx q[9],q[11],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[21],q[43];
ccx q[9],q[11],q[42];

// red K4 (10, 12, 24, 26) -> avoid connecting to all four
ccx q[9],q[11],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[25],q[44];
rz(0.15) q[44];
ccx q[43],q[25],q[44];
ccx q[42],q[23],q[43];
ccx q[9],q[11],q[42];

// red K4 (10, 12, 24, 37) -> avoid connecting to all four
ccx q[9],q[11],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[23],q[43];
ccx q[9],q[11],q[42];

// red K4 (10, 12, 26, 28) -> avoid connecting to all four
ccx q[9],q[11],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[25],q[43];
ccx q[9],q[11],q[42];

// red K4 (10, 12, 26, 33) -> avoid connecting to all four
ccx q[9],q[11],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[25],q[43];
ccx q[9],q[11],q[42];

// red K4 (10, 12, 26, 39) -> avoid connecting to all four
ccx q[9],q[11],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[25],q[43];
ccx q[9],q[11],q[42];

// red K4 (10, 12, 28, 30) -> avoid connecting to all four
ccx q[9],q[11],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[27],q[43];
ccx q[9],q[11],q[42];

// red K4 (10, 12, 28, 35) -> avoid connecting to all four
ccx q[9],q[11],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[27],q[43];
ccx q[9],q[11],q[42];

// red K4 (10, 12, 28, 41) -> avoid connecting to all four
ccx q[9],q[11],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[27],q[43];
ccx q[9],q[11],q[42];

// red K4 (10, 12, 30, 32) -> avoid connecting to all four
ccx q[9],q[11],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[29],q[43];
ccx q[9],q[11],q[42];

// red K4 (10, 12, 30, 37) -> avoid connecting to all four
ccx q[9],q[11],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[29],q[43];
ccx q[9],q[11],q[42];

// red K4 (10, 12, 32, 33) -> avoid connecting to all four
ccx q[9],q[11],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[31],q[43];
ccx q[9],q[11],q[42];

// red K4 (10, 12, 32, 39) -> avoid connecting to all four
ccx q[9],q[11],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[31],q[43];
ccx q[9],q[11],q[42];

// red K4 (10, 12, 33, 35) -> avoid connecting to all four
ccx q[9],q[11],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[32],q[43];
ccx q[9],q[11],q[42];

// red K4 (10, 12, 35, 37) -> avoid connecting to all four
ccx q[9],q[11],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[34],q[43];
ccx q[9],q[11],q[42];

// red K4 (10, 12, 37, 39) -> avoid connecting to all four
ccx q[9],q[11],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[36],q[43];
ccx q[9],q[11],q[42];

// red K4 (10, 12, 39, 41) -> avoid connecting to all four
ccx q[9],q[11],q[42];
ccx q[42],q[38],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[38],q[43];
ccx q[9],q[11],q[42];

// red K4 (10, 17, 24, 31) -> avoid connecting to all four
ccx q[9],q[16],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[23],q[43];
ccx q[9],q[16],q[42];

// red K4 (10, 17, 24, 37) -> avoid connecting to all four
ccx q[9],q[16],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[23],q[43];
ccx q[9],q[16],q[42];

// red K4 (10, 17, 24, 40) -> avoid connecting to all four
ccx q[9],q[16],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[23],q[43];
ccx q[9],q[16],q[42];

// red K4 (10, 17, 30, 37) -> avoid connecting to all four
ccx q[9],q[16],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[29],q[43];
ccx q[9],q[16],q[42];

// red K4 (10, 17, 30, 40) -> avoid connecting to all four
ccx q[9],q[16],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[29],q[43];
ccx q[9],q[16],q[42];

// red K4 (10, 17, 31, 33) -> avoid connecting to all four
ccx q[9],q[16],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[30],q[43];
ccx q[9],q[16],q[42];

// red K4 (10, 17, 33, 35) -> avoid connecting to all four
ccx q[9],q[16],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[32],q[43];
ccx q[9],q[16],q[42];

// red K4 (10, 17, 33, 40) -> avoid connecting to all four
ccx q[9],q[16],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[32],q[43];
ccx q[9],q[16],q[42];

// red K4 (10, 17, 35, 37) -> avoid connecting to all four
ccx q[9],q[16],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[34],q[43];
ccx q[9],q[16],q[42];

// red K4 (10, 17, 37, 39) -> avoid connecting to all four
ccx q[9],q[16],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[36],q[43];
ccx q[9],q[16],q[42];

// red K4 (10, 20, 22, 32) -> avoid connecting to all four
ccx q[9],q[19],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[21],q[43];
ccx q[9],q[19],q[42];

// red K4 (10, 20, 22, 40) -> avoid connecting to all four
ccx q[9],q[19],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[21],q[43];
ccx q[9],q[19],q[42];

// red K4 (10, 20, 30, 32) -> avoid connecting to all four
ccx q[9],q[19],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[29],q[43];
ccx q[9],q[19],q[42];

// red K4 (10, 20, 30, 40) -> avoid connecting to all four
ccx q[9],q[19],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[29],q[43];
ccx q[9],q[19],q[42];

// red K4 (10, 20, 32, 33) -> avoid connecting to all four
ccx q[9],q[19],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[31],q[43];
ccx q[9],q[19],q[42];

// red K4 (10, 20, 33, 40) -> avoid connecting to all four
ccx q[9],q[19],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[32],q[43];
ccx q[9],q[19],q[42];

// red K4 (10, 22, 23, 35) -> avoid connecting to all four
ccx q[9],q[21],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[22],q[43];
ccx q[9],q[21],q[42];

// red K4 (10, 22, 24, 40) -> avoid connecting to all four
ccx q[9],q[21],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[23],q[43];
ccx q[9],q[21],q[42];

// red K4 (10, 23, 30, 37) -> avoid connecting to all four
ccx q[9],q[22],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[29],q[43];
ccx q[9],q[22],q[42];

// red K4 (10, 23, 33, 35) -> avoid connecting to all four
ccx q[9],q[22],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[32],q[43];
ccx q[9],q[22],q[42];

// red K4 (10, 23, 35, 37) -> avoid connecting to all four
ccx q[9],q[22],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[34],q[43];
ccx q[9],q[22],q[42];

// red K4 (10, 23, 37, 39) -> avoid connecting to all four
ccx q[9],q[22],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[36],q[43];
ccx q[9],q[22],q[42];

// red K4 (10, 23, 39, 41) -> avoid connecting to all four
ccx q[9],q[22],q[42];
ccx q[42],q[38],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[38],q[43];
ccx q[9],q[22],q[42];

// red K4 (10, 24, 26, 40) -> avoid connecting to all four
ccx q[9],q[23],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[25],q[43];
ccx q[9],q[23],q[42];

// red K4 (10, 26, 28, 40) -> avoid connecting to all four
ccx q[9],q[25],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[27],q[43];
ccx q[9],q[25],q[42];

// red K4 (10, 26, 33, 40) -> avoid connecting to all four
ccx q[9],q[25],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[32],q[43];
ccx q[9],q[25],q[42];

// red K4 (10, 28, 30, 40) -> avoid connecting to all four
ccx q[9],q[27],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[29],q[43];
ccx q[9],q[27],q[42];

// red K4 (10, 31, 32, 33) -> avoid connecting to all four
ccx q[9],q[30],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[31],q[43];
ccx q[9],q[30],q[42];

// red K4 (11, 12, 13, 25) -> avoid connecting to all four
ccx q[10],q[11],q[42];
ccx q[42],q[12],q[43];
ccx q[43],q[24],q[44];
rz(0.15) q[44];
ccx q[43],q[24],q[44];
ccx q[42],q[12],q[43];
ccx q[10],q[11],q[42];

// red K4 (11, 12, 13, 33) -> avoid connecting to all four
ccx q[10],q[11],q[42];
ccx q[42],q[12],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[12],q[43];
ccx q[10],q[11],q[42];

// red K4 (11, 12, 13, 34) -> avoid connecting to all four
ccx q[10],q[11],q[42];
ccx q[42],q[12],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[12],q[43];
ccx q[10],q[11],q[42];

// red K4 (11, 12, 13, 42) -> avoid connecting to all four
ccx q[10],q[11],q[42];
ccx q[42],q[12],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[12],q[43];
ccx q[10],q[11],q[42];

// red K4 (11, 12, 24, 34) -> avoid connecting to all four
ccx q[10],q[11],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[23],q[43];
ccx q[10],q[11],q[42];

// red K4 (11, 12, 24, 42) -> avoid connecting to all four
ccx q[10],q[11],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[23],q[43];
ccx q[10],q[11],q[42];

// red K4 (11, 12, 25, 41) -> avoid connecting to all four
ccx q[10],q[11],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[24],q[43];
ccx q[10],q[11],q[42];

// red K4 (11, 12, 34, 41) -> avoid connecting to all four
ccx q[10],q[11],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[33],q[43];
ccx q[10],q[11],q[42];

// red K4 (11, 13, 23, 25) -> avoid connecting to all four
ccx q[10],q[12],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[24],q[44];
rz(0.15) q[44];
ccx q[43],q[24],q[44];
ccx q[42],q[22],q[43];
ccx q[10],q[12],q[42];

// red K4 (11, 13, 23, 33) -> avoid connecting to all four
ccx q[10],q[12],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[22],q[43];
ccx q[10],q[12],q[42];

// red K4 (11, 13, 23, 36) -> avoid connecting to all four
ccx q[10],q[12],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[22],q[43];
ccx q[10],q[12],q[42];

// red K4 (11, 13, 25, 27) -> avoid connecting to all four
ccx q[10],q[12],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[26],q[44];
rz(0.15) q[44];
ccx q[43],q[26],q[44];
ccx q[42],q[24],q[43];
ccx q[10],q[12],q[42];

// red K4 (11, 13, 25, 38) -> avoid connecting to all four
ccx q[10],q[12],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[24],q[43];
ccx q[10],q[12],q[42];

// red K4 (11, 13, 27, 29) -> avoid connecting to all four
ccx q[10],q[12],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[26],q[43];
ccx q[10],q[12],q[42];

// red K4 (11, 13, 27, 34) -> avoid connecting to all four
ccx q[10],q[12],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[26],q[43];
ccx q[10],q[12],q[42];

// red K4 (11, 13, 27, 40) -> avoid connecting to all four
ccx q[10],q[12],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[26],q[43];
ccx q[10],q[12],q[42];

// red K4 (11, 13, 29, 31) -> avoid connecting to all four
ccx q[10],q[12],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[28],q[43];
ccx q[10],q[12],q[42];

// red K4 (11, 13, 29, 36) -> avoid connecting to all four
ccx q[10],q[12],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[28],q[43];
ccx q[10],q[12],q[42];

// red K4 (11, 13, 29, 42) -> avoid connecting to all four
ccx q[10],q[12],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[28],q[43];
ccx q[10],q[12],q[42];

// red K4 (11, 13, 31, 33) -> avoid connecting to all four
ccx q[10],q[12],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[30],q[43];
ccx q[10],q[12],q[42];

// red K4 (11, 13, 31, 38) -> avoid connecting to all four
ccx q[10],q[12],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[30],q[43];
ccx q[10],q[12],q[42];

// red K4 (11, 13, 33, 40) -> avoid connecting to all four
ccx q[10],q[12],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[32],q[43];
ccx q[10],q[12],q[42];

// red K4 (11, 13, 34, 36) -> avoid connecting to all four
ccx q[10],q[12],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[33],q[43];
ccx q[10],q[12],q[42];

// red K4 (11, 13, 36, 38) -> avoid connecting to all four
ccx q[10],q[12],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[35],q[43];
ccx q[10],q[12],q[42];

// red K4 (11, 13, 38, 40) -> avoid connecting to all four
ccx q[10],q[12],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[37],q[43];
ccx q[10],q[12],q[42];

// red K4 (11, 13, 40, 42) -> avoid connecting to all four
ccx q[10],q[12],q[42];
ccx q[42],q[39],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[39],q[43];
ccx q[10],q[12],q[42];

// red K4 (11, 18, 25, 38) -> avoid connecting to all four
ccx q[10],q[17],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[24],q[43];
ccx q[10],q[17],q[42];

// red K4 (11, 18, 25, 41) -> avoid connecting to all four
ccx q[10],q[17],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[24],q[43];
ccx q[10],q[17],q[42];

// red K4 (11, 18, 31, 38) -> avoid connecting to all four
ccx q[10],q[17],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[30],q[43];
ccx q[10],q[17],q[42];

// red K4 (11, 18, 31, 41) -> avoid connecting to all four
ccx q[10],q[17],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[30],q[43];
ccx q[10],q[17],q[42];

// red K4 (11, 18, 34, 36) -> avoid connecting to all four
ccx q[10],q[17],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[33],q[43];
ccx q[10],q[17],q[42];

// red K4 (11, 18, 34, 41) -> avoid connecting to all four
ccx q[10],q[17],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[33],q[43];
ccx q[10],q[17],q[42];

// red K4 (11, 18, 36, 38) -> avoid connecting to all four
ccx q[10],q[17],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[35],q[43];
ccx q[10],q[17],q[42];

// red K4 (11, 18, 38, 40) -> avoid connecting to all four
ccx q[10],q[17],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[37],q[43];
ccx q[10],q[17],q[42];

// red K4 (11, 21, 23, 33) -> avoid connecting to all four
ccx q[10],q[20],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[22],q[43];
ccx q[10],q[20],q[42];

// red K4 (11, 21, 23, 41) -> avoid connecting to all four
ccx q[10],q[20],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[22],q[43];
ccx q[10],q[20],q[42];

// red K4 (11, 21, 31, 33) -> avoid connecting to all four
ccx q[10],q[20],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[30],q[43];
ccx q[10],q[20],q[42];

// red K4 (11, 21, 31, 41) -> avoid connecting to all four
ccx q[10],q[20],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[30],q[43];
ccx q[10],q[20],q[42];

// red K4 (11, 21, 34, 41) -> avoid connecting to all four
ccx q[10],q[20],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[33],q[43];
ccx q[10],q[20],q[42];

// red K4 (11, 23, 25, 41) -> avoid connecting to all four
ccx q[10],q[22],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[24],q[43];
ccx q[10],q[22],q[42];

// red K4 (11, 24, 31, 38) -> avoid connecting to all four
ccx q[10],q[23],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[30],q[43];
ccx q[10],q[23],q[42];

// red K4 (11, 24, 34, 36) -> avoid connecting to all four
ccx q[10],q[23],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[33],q[43];
ccx q[10],q[23],q[42];

// red K4 (11, 24, 36, 38) -> avoid connecting to all four
ccx q[10],q[23],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[35],q[43];
ccx q[10],q[23],q[42];

// red K4 (11, 24, 38, 40) -> avoid connecting to all four
ccx q[10],q[23],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[37],q[43];
ccx q[10],q[23],q[42];

// red K4 (11, 24, 40, 42) -> avoid connecting to all four
ccx q[10],q[23],q[42];
ccx q[42],q[39],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[39],q[43];
ccx q[10],q[23],q[42];

// red K4 (11, 25, 27, 41) -> avoid connecting to all four
ccx q[10],q[24],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[26],q[43];
ccx q[10],q[24],q[42];

// red K4 (11, 27, 29, 41) -> avoid connecting to all four
ccx q[10],q[26],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[28],q[43];
ccx q[10],q[26],q[42];

// red K4 (11, 27, 34, 41) -> avoid connecting to all four
ccx q[10],q[26],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[33],q[43];
ccx q[10],q[26],q[42];

// red K4 (11, 29, 31, 41) -> avoid connecting to all four
ccx q[10],q[28],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[30],q[43];
ccx q[10],q[28],q[42];

// red K4 (12, 13, 25, 26) -> avoid connecting to all four
ccx q[11],q[12],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[25],q[44];
rz(0.15) q[44];
ccx q[43],q[25],q[44];
ccx q[42],q[24],q[43];
ccx q[11],q[12],q[42];

// red K4 (12, 13, 25, 35) -> avoid connecting to all four
ccx q[11],q[12],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[24],q[43];
ccx q[11],q[12],q[42];

// red K4 (12, 13, 26, 33) -> avoid connecting to all four
ccx q[11],q[12],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[25],q[43];
ccx q[11],q[12],q[42];

// red K4 (12, 13, 26, 42) -> avoid connecting to all four
ccx q[11],q[12],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[25],q[43];
ccx q[11],q[12],q[42];

// red K4 (12, 13, 33, 35) -> avoid connecting to all four
ccx q[11],q[12],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[32],q[43];
ccx q[11],q[12],q[42];

// red K4 (12, 13, 34, 35) -> avoid connecting to all four
ccx q[11],q[12],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[33],q[43];
ccx q[11],q[12],q[42];

// red K4 (12, 13, 35, 42) -> avoid connecting to all four
ccx q[11],q[12],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[34],q[43];
ccx q[11],q[12],q[42];

// red K4 (12, 14, 24, 26) -> avoid connecting to all four
ccx q[11],q[13],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[25],q[44];
rz(0.15) q[44];
ccx q[43],q[25],q[44];
ccx q[42],q[23],q[43];
ccx q[11],q[13],q[42];

// red K4 (12, 14, 24, 34) -> avoid connecting to all four
ccx q[11],q[13],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[23],q[43];
ccx q[11],q[13],q[42];

// red K4 (12, 14, 24, 37) -> avoid connecting to all four
ccx q[11],q[13],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[23],q[43];
ccx q[11],q[13],q[42];

// red K4 (12, 14, 26, 28) -> avoid connecting to all four
ccx q[11],q[13],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[25],q[43];
ccx q[11],q[13],q[42];

// red K4 (12, 14, 26, 39) -> avoid connecting to all four
ccx q[11],q[13],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[25],q[43];
ccx q[11],q[13],q[42];

// red K4 (12, 14, 28, 30) -> avoid connecting to all four
ccx q[11],q[13],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[27],q[43];
ccx q[11],q[13],q[42];

// red K4 (12, 14, 28, 35) -> avoid connecting to all four
ccx q[11],q[13],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[27],q[43];
ccx q[11],q[13],q[42];

// red K4 (12, 14, 28, 41) -> avoid connecting to all four
ccx q[11],q[13],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[27],q[43];
ccx q[11],q[13],q[42];

// red K4 (12, 14, 30, 32) -> avoid connecting to all four
ccx q[11],q[13],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[29],q[43];
ccx q[11],q[13],q[42];

// red K4 (12, 14, 30, 37) -> avoid connecting to all four
ccx q[11],q[13],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[29],q[43];
ccx q[11],q[13],q[42];

// red K4 (12, 14, 32, 34) -> avoid connecting to all four
ccx q[11],q[13],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[31],q[43];
ccx q[11],q[13],q[42];

// red K4 (12, 14, 32, 39) -> avoid connecting to all four
ccx q[11],q[13],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[31],q[43];
ccx q[11],q[13],q[42];

// red K4 (12, 14, 34, 35) -> avoid connecting to all four
ccx q[11],q[13],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[33],q[43];
ccx q[11],q[13],q[42];

// red K4 (12, 14, 34, 41) -> avoid connecting to all four
ccx q[11],q[13],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[33],q[43];
ccx q[11],q[13],q[42];

// red K4 (12, 14, 35, 37) -> avoid connecting to all four
ccx q[11],q[13],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[34],q[43];
ccx q[11],q[13],q[42];

// red K4 (12, 14, 37, 39) -> avoid connecting to all four
ccx q[11],q[13],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[36],q[43];
ccx q[11],q[13],q[42];

// red K4 (12, 14, 39, 41) -> avoid connecting to all four
ccx q[11],q[13],q[42];
ccx q[42],q[38],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[38],q[43];
ccx q[11],q[13],q[42];

// red K4 (12, 19, 26, 33) -> avoid connecting to all four
ccx q[11],q[18],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[25],q[43];
ccx q[11],q[18],q[42];

// red K4 (12, 19, 26, 39) -> avoid connecting to all four
ccx q[11],q[18],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[25],q[43];
ccx q[11],q[18],q[42];

// red K4 (12, 19, 26, 42) -> avoid connecting to all four
ccx q[11],q[18],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[25],q[43];
ccx q[11],q[18],q[42];

// red K4 (12, 19, 32, 33) -> avoid connecting to all four
ccx q[11],q[18],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[31],q[43];
ccx q[11],q[18],q[42];

// red K4 (12, 19, 32, 39) -> avoid connecting to all four
ccx q[11],q[18],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[31],q[43];
ccx q[11],q[18],q[42];

// red K4 (12, 19, 32, 42) -> avoid connecting to all four
ccx q[11],q[18],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[31],q[43];
ccx q[11],q[18],q[42];

// red K4 (12, 19, 33, 35) -> avoid connecting to all four
ccx q[11],q[18],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[32],q[43];
ccx q[11],q[18],q[42];

// red K4 (12, 19, 35, 37) -> avoid connecting to all four
ccx q[11],q[18],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[34],q[43];
ccx q[11],q[18],q[42];

// red K4 (12, 19, 35, 42) -> avoid connecting to all four
ccx q[11],q[18],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[34],q[43];
ccx q[11],q[18],q[42];

// red K4 (12, 19, 37, 39) -> avoid connecting to all four
ccx q[11],q[18],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[36],q[43];
ccx q[11],q[18],q[42];

// red K4 (12, 19, 39, 41) -> avoid connecting to all four
ccx q[11],q[18],q[42];
ccx q[42],q[38],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[38],q[43];
ccx q[11],q[18],q[42];

// red K4 (12, 22, 24, 34) -> avoid connecting to all four
ccx q[11],q[21],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[23],q[43];
ccx q[11],q[21],q[42];

// red K4 (12, 22, 24, 42) -> avoid connecting to all four
ccx q[11],q[21],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[23],q[43];
ccx q[11],q[21],q[42];

// red K4 (12, 22, 32, 34) -> avoid connecting to all four
ccx q[11],q[21],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[31],q[43];
ccx q[11],q[21],q[42];

// red K4 (12, 22, 32, 42) -> avoid connecting to all four
ccx q[11],q[21],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[31],q[43];
ccx q[11],q[21],q[42];

// red K4 (12, 22, 34, 35) -> avoid connecting to all four
ccx q[11],q[21],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[33],q[43];
ccx q[11],q[21],q[42];

// red K4 (12, 22, 35, 42) -> avoid connecting to all four
ccx q[11],q[21],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[34],q[43];
ccx q[11],q[21],q[42];

// red K4 (12, 24, 26, 42) -> avoid connecting to all four
ccx q[11],q[23],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[25],q[43];
ccx q[11],q[23],q[42];

// red K4 (12, 25, 26, 39) -> avoid connecting to all four
ccx q[11],q[24],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[25],q[43];
ccx q[11],q[24],q[42];

// red K4 (12, 25, 32, 39) -> avoid connecting to all four
ccx q[11],q[24],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[31],q[43];
ccx q[11],q[24],q[42];

// red K4 (12, 25, 35, 37) -> avoid connecting to all four
ccx q[11],q[24],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[34],q[43];
ccx q[11],q[24],q[42];

// red K4 (12, 25, 37, 39) -> avoid connecting to all four
ccx q[11],q[24],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[36],q[43];
ccx q[11],q[24],q[42];

// red K4 (12, 25, 39, 41) -> avoid connecting to all four
ccx q[11],q[24],q[42];
ccx q[42],q[38],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[38],q[43];
ccx q[11],q[24],q[42];

// red K4 (12, 26, 28, 42) -> avoid connecting to all four
ccx q[11],q[25],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[27],q[43];
ccx q[11],q[25],q[42];

// red K4 (12, 28, 30, 42) -> avoid connecting to all four
ccx q[11],q[27],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[29],q[43];
ccx q[11],q[27],q[42];

// red K4 (12, 28, 35, 42) -> avoid connecting to all four
ccx q[11],q[27],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[34],q[43];
ccx q[11],q[27],q[42];

// red K4 (12, 30, 32, 42) -> avoid connecting to all four
ccx q[11],q[29],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[31],q[43];
ccx q[11],q[29],q[42];

// red K4 (13, 15, 25, 27) -> avoid connecting to all four
ccx q[12],q[14],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[26],q[44];
rz(0.15) q[44];
ccx q[43],q[26],q[44];
ccx q[42],q[24],q[43];
ccx q[12],q[14],q[42];

// red K4 (13, 15, 25, 35) -> avoid connecting to all four
ccx q[12],q[14],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[24],q[43];
ccx q[12],q[14],q[42];

// red K4 (13, 15, 25, 38) -> avoid connecting to all four
ccx q[12],q[14],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[24],q[43];
ccx q[12],q[14],q[42];

// red K4 (13, 15, 27, 29) -> avoid connecting to all four
ccx q[12],q[14],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[26],q[43];
ccx q[12],q[14],q[42];

// red K4 (13, 15, 27, 40) -> avoid connecting to all four
ccx q[12],q[14],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[26],q[43];
ccx q[12],q[14],q[42];

// red K4 (13, 15, 29, 31) -> avoid connecting to all four
ccx q[12],q[14],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[28],q[43];
ccx q[12],q[14],q[42];

// red K4 (13, 15, 29, 36) -> avoid connecting to all four
ccx q[12],q[14],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[28],q[43];
ccx q[12],q[14],q[42];

// red K4 (13, 15, 29, 42) -> avoid connecting to all four
ccx q[12],q[14],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[28],q[43];
ccx q[12],q[14],q[42];

// red K4 (13, 15, 31, 33) -> avoid connecting to all four
ccx q[12],q[14],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[30],q[43];
ccx q[12],q[14],q[42];

// red K4 (13, 15, 31, 38) -> avoid connecting to all four
ccx q[12],q[14],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[30],q[43];
ccx q[12],q[14],q[42];

// red K4 (13, 15, 33, 35) -> avoid connecting to all four
ccx q[12],q[14],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[32],q[43];
ccx q[12],q[14],q[42];

// red K4 (13, 15, 33, 40) -> avoid connecting to all four
ccx q[12],q[14],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[32],q[43];
ccx q[12],q[14],q[42];

// red K4 (13, 15, 35, 36) -> avoid connecting to all four
ccx q[12],q[14],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[34],q[43];
ccx q[12],q[14],q[42];

// red K4 (13, 15, 35, 42) -> avoid connecting to all four
ccx q[12],q[14],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[34],q[43];
ccx q[12],q[14],q[42];

// red K4 (13, 15, 36, 38) -> avoid connecting to all four
ccx q[12],q[14],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[35],q[43];
ccx q[12],q[14],q[42];

// red K4 (13, 15, 38, 40) -> avoid connecting to all four
ccx q[12],q[14],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[37],q[43];
ccx q[12],q[14],q[42];

// red K4 (13, 15, 40, 42) -> avoid connecting to all four
ccx q[12],q[14],q[42];
ccx q[42],q[39],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[39],q[43];
ccx q[12],q[14],q[42];

// red K4 (13, 20, 27, 34) -> avoid connecting to all four
ccx q[12],q[19],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[26],q[43];
ccx q[12],q[19],q[42];

// red K4 (13, 20, 27, 40) -> avoid connecting to all four
ccx q[12],q[19],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[26],q[43];
ccx q[12],q[19],q[42];

// red K4 (13, 20, 33, 40) -> avoid connecting to all four
ccx q[12],q[19],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[32],q[43];
ccx q[12],q[19],q[42];

// red K4 (13, 20, 34, 36) -> avoid connecting to all four
ccx q[12],q[19],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[33],q[43];
ccx q[12],q[19],q[42];

// red K4 (13, 20, 36, 38) -> avoid connecting to all four
ccx q[12],q[19],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[35],q[43];
ccx q[12],q[19],q[42];

// red K4 (13, 20, 38, 40) -> avoid connecting to all four
ccx q[12],q[19],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[37],q[43];
ccx q[12],q[19],q[42];

// red K4 (13, 20, 40, 42) -> avoid connecting to all four
ccx q[12],q[19],q[42];
ccx q[42],q[39],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[39],q[43];
ccx q[12],q[19],q[42];

// red K4 (13, 23, 25, 35) -> avoid connecting to all four
ccx q[12],q[22],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[24],q[43];
ccx q[12],q[22],q[42];

// red K4 (13, 23, 33, 35) -> avoid connecting to all four
ccx q[12],q[22],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[32],q[43];
ccx q[12],q[22],q[42];

// red K4 (13, 23, 35, 36) -> avoid connecting to all four
ccx q[12],q[22],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[34],q[43];
ccx q[12],q[22],q[42];

// red K4 (13, 25, 26, 27) -> avoid connecting to all four
ccx q[12],q[24],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[26],q[44];
rz(0.15) q[44];
ccx q[43],q[26],q[44];
ccx q[42],q[25],q[43];
ccx q[12],q[24],q[42];

// red K4 (13, 25, 26, 38) -> avoid connecting to all four
ccx q[12],q[24],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[25],q[43];
ccx q[12],q[24],q[42];

// red K4 (13, 26, 27, 40) -> avoid connecting to all four
ccx q[12],q[25],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[26],q[43];
ccx q[12],q[25],q[42];

// red K4 (13, 26, 33, 40) -> avoid connecting to all four
ccx q[12],q[25],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[32],q[43];
ccx q[12],q[25],q[42];

// red K4 (13, 26, 36, 38) -> avoid connecting to all four
ccx q[12],q[25],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[35],q[43];
ccx q[12],q[25],q[42];

// red K4 (13, 26, 38, 40) -> avoid connecting to all four
ccx q[12],q[25],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[37],q[43];
ccx q[12],q[25],q[42];

// red K4 (13, 26, 40, 42) -> avoid connecting to all four
ccx q[12],q[25],q[42];
ccx q[42],q[39],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[39],q[43];
ccx q[12],q[25],q[42];

// red K4 (13, 34, 35, 36) -> avoid connecting to all four
ccx q[12],q[33],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[34],q[43];
ccx q[12],q[33],q[42];

// red K4 (14, 16, 26, 28) -> avoid connecting to all four
ccx q[13],q[15],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[25],q[43];
ccx q[13],q[15],q[42];

// red K4 (14, 16, 26, 36) -> avoid connecting to all four
ccx q[13],q[15],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[25],q[43];
ccx q[13],q[15],q[42];

// red K4 (14, 16, 26, 39) -> avoid connecting to all four
ccx q[13],q[15],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[25],q[43];
ccx q[13],q[15],q[42];

// red K4 (14, 16, 28, 30) -> avoid connecting to all four
ccx q[13],q[15],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[27],q[43];
ccx q[13],q[15],q[42];

// red K4 (14, 16, 28, 41) -> avoid connecting to all four
ccx q[13],q[15],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[27],q[43];
ccx q[13],q[15],q[42];

// red K4 (14, 16, 30, 32) -> avoid connecting to all four
ccx q[13],q[15],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[29],q[43];
ccx q[13],q[15],q[42];

// red K4 (14, 16, 30, 37) -> avoid connecting to all four
ccx q[13],q[15],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[29],q[43];
ccx q[13],q[15],q[42];

// red K4 (14, 16, 32, 34) -> avoid connecting to all four
ccx q[13],q[15],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[31],q[43];
ccx q[13],q[15],q[42];

// red K4 (14, 16, 32, 39) -> avoid connecting to all four
ccx q[13],q[15],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[31],q[43];
ccx q[13],q[15],q[42];

// red K4 (14, 16, 34, 36) -> avoid connecting to all four
ccx q[13],q[15],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[33],q[43];
ccx q[13],q[15],q[42];

// red K4 (14, 16, 34, 41) -> avoid connecting to all four
ccx q[13],q[15],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[33],q[43];
ccx q[13],q[15],q[42];

// red K4 (14, 16, 36, 37) -> avoid connecting to all four
ccx q[13],q[15],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[35],q[43];
ccx q[13],q[15],q[42];

// red K4 (14, 16, 37, 39) -> avoid connecting to all four
ccx q[13],q[15],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[36],q[43];
ccx q[13],q[15],q[42];

// red K4 (14, 16, 39, 41) -> avoid connecting to all four
ccx q[13],q[15],q[42];
ccx q[42],q[38],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[38],q[43];
ccx q[13],q[15],q[42];

// red K4 (14, 21, 28, 35) -> avoid connecting to all four
ccx q[13],q[20],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[27],q[43];
ccx q[13],q[20],q[42];

// red K4 (14, 21, 28, 41) -> avoid connecting to all four
ccx q[13],q[20],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[27],q[43];
ccx q[13],q[20],q[42];

// red K4 (14, 21, 34, 35) -> avoid connecting to all four
ccx q[13],q[20],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[33],q[43];
ccx q[13],q[20],q[42];

// red K4 (14, 21, 34, 41) -> avoid connecting to all four
ccx q[13],q[20],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[33],q[43];
ccx q[13],q[20],q[42];

// red K4 (14, 21, 35, 37) -> avoid connecting to all four
ccx q[13],q[20],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[34],q[43];
ccx q[13],q[20],q[42];

// red K4 (14, 21, 37, 39) -> avoid connecting to all four
ccx q[13],q[20],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[36],q[43];
ccx q[13],q[20],q[42];

// red K4 (14, 21, 39, 41) -> avoid connecting to all four
ccx q[13],q[20],q[42];
ccx q[42],q[38],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[38],q[43];
ccx q[13],q[20],q[42];

// red K4 (14, 24, 26, 36) -> avoid connecting to all four
ccx q[13],q[23],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[25],q[43];
ccx q[13],q[23],q[42];

// red K4 (14, 24, 34, 36) -> avoid connecting to all four
ccx q[13],q[23],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[33],q[43];
ccx q[13],q[23],q[42];

// red K4 (14, 24, 36, 37) -> avoid connecting to all four
ccx q[13],q[23],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[35],q[43];
ccx q[13],q[23],q[42];

// red K4 (14, 26, 27, 28) -> avoid connecting to all four
ccx q[13],q[25],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[26],q[43];
ccx q[13],q[25],q[42];

// red K4 (14, 26, 27, 39) -> avoid connecting to all four
ccx q[13],q[25],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[26],q[43];
ccx q[13],q[25],q[42];

// red K4 (14, 27, 28, 41) -> avoid connecting to all four
ccx q[13],q[26],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[27],q[43];
ccx q[13],q[26],q[42];

// red K4 (14, 27, 34, 41) -> avoid connecting to all four
ccx q[13],q[26],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[33],q[43];
ccx q[13],q[26],q[42];

// red K4 (14, 27, 37, 39) -> avoid connecting to all four
ccx q[13],q[26],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[36],q[43];
ccx q[13],q[26],q[42];

// red K4 (14, 27, 39, 41) -> avoid connecting to all four
ccx q[13],q[26],q[42];
ccx q[42],q[38],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[38],q[43];
ccx q[13],q[26],q[42];

// red K4 (14, 34, 35, 36) -> avoid connecting to all four
ccx q[13],q[33],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[34],q[43];
ccx q[13],q[33],q[42];

// red K4 (14, 35, 36, 37) -> avoid connecting to all four
ccx q[13],q[34],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[35],q[43];
ccx q[13],q[34],q[42];

// red K4 (15, 17, 27, 29) -> avoid connecting to all four
ccx q[14],q[16],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[26],q[43];
ccx q[14],q[16],q[42];

// red K4 (15, 17, 27, 37) -> avoid connecting to all four
ccx q[14],q[16],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[26],q[43];
ccx q[14],q[16],q[42];

// red K4 (15, 17, 27, 40) -> avoid connecting to all four
ccx q[14],q[16],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[26],q[43];
ccx q[14],q[16],q[42];

// red K4 (15, 17, 29, 31) -> avoid connecting to all four
ccx q[14],q[16],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[28],q[43];
ccx q[14],q[16],q[42];

// red K4 (15, 17, 29, 42) -> avoid connecting to all four
ccx q[14],q[16],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[28],q[43];
ccx q[14],q[16],q[42];

// red K4 (15, 17, 31, 33) -> avoid connecting to all four
ccx q[14],q[16],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[30],q[43];
ccx q[14],q[16],q[42];

// red K4 (15, 17, 31, 38) -> avoid connecting to all four
ccx q[14],q[16],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[30],q[43];
ccx q[14],q[16],q[42];

// red K4 (15, 17, 33, 35) -> avoid connecting to all four
ccx q[14],q[16],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[32],q[43];
ccx q[14],q[16],q[42];

// red K4 (15, 17, 33, 40) -> avoid connecting to all four
ccx q[14],q[16],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[32],q[43];
ccx q[14],q[16],q[42];

// red K4 (15, 17, 35, 37) -> avoid connecting to all four
ccx q[14],q[16],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[34],q[43];
ccx q[14],q[16],q[42];

// red K4 (15, 17, 35, 42) -> avoid connecting to all four
ccx q[14],q[16],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[34],q[43];
ccx q[14],q[16],q[42];

// red K4 (15, 17, 37, 38) -> avoid connecting to all four
ccx q[14],q[16],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[36],q[43];
ccx q[14],q[16],q[42];

// red K4 (15, 17, 38, 40) -> avoid connecting to all four
ccx q[14],q[16],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[37],q[43];
ccx q[14],q[16],q[42];

// red K4 (15, 17, 40, 42) -> avoid connecting to all four
ccx q[14],q[16],q[42];
ccx q[42],q[39],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[39],q[43];
ccx q[14],q[16],q[42];

// red K4 (15, 22, 29, 36) -> avoid connecting to all four
ccx q[14],q[21],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[28],q[43];
ccx q[14],q[21],q[42];

// red K4 (15, 22, 29, 42) -> avoid connecting to all four
ccx q[14],q[21],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[28],q[43];
ccx q[14],q[21],q[42];

// red K4 (15, 22, 35, 36) -> avoid connecting to all four
ccx q[14],q[21],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[34],q[43];
ccx q[14],q[21],q[42];

// red K4 (15, 22, 35, 42) -> avoid connecting to all four
ccx q[14],q[21],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[34],q[43];
ccx q[14],q[21],q[42];

// red K4 (15, 22, 36, 38) -> avoid connecting to all four
ccx q[14],q[21],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[35],q[43];
ccx q[14],q[21],q[42];

// red K4 (15, 22, 38, 40) -> avoid connecting to all four
ccx q[14],q[21],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[37],q[43];
ccx q[14],q[21],q[42];

// red K4 (15, 22, 40, 42) -> avoid connecting to all four
ccx q[14],q[21],q[42];
ccx q[42],q[39],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[39],q[43];
ccx q[14],q[21],q[42];

// red K4 (15, 25, 27, 37) -> avoid connecting to all four
ccx q[14],q[24],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[26],q[43];
ccx q[14],q[24],q[42];

// red K4 (15, 25, 35, 37) -> avoid connecting to all four
ccx q[14],q[24],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[34],q[43];
ccx q[14],q[24],q[42];

// red K4 (15, 25, 37, 38) -> avoid connecting to all four
ccx q[14],q[24],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[36],q[43];
ccx q[14],q[24],q[42];

// red K4 (15, 27, 28, 29) -> avoid connecting to all four
ccx q[14],q[26],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[27],q[43];
ccx q[14],q[26],q[42];

// red K4 (15, 27, 28, 40) -> avoid connecting to all four
ccx q[14],q[26],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[27],q[43];
ccx q[14],q[26],q[42];

// red K4 (15, 28, 29, 42) -> avoid connecting to all four
ccx q[14],q[27],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[28],q[43];
ccx q[14],q[27],q[42];

// red K4 (15, 28, 35, 42) -> avoid connecting to all four
ccx q[14],q[27],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[34],q[43];
ccx q[14],q[27],q[42];

// red K4 (15, 28, 38, 40) -> avoid connecting to all four
ccx q[14],q[27],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[37],q[43];
ccx q[14],q[27],q[42];

// red K4 (15, 28, 40, 42) -> avoid connecting to all four
ccx q[14],q[27],q[42];
ccx q[42],q[39],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[39],q[43];
ccx q[14],q[27],q[42];

// red K4 (15, 35, 36, 37) -> avoid connecting to all four
ccx q[14],q[34],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[35],q[43];
ccx q[14],q[34],q[42];

// red K4 (15, 36, 37, 38) -> avoid connecting to all four
ccx q[14],q[35],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[36],q[43];
ccx q[14],q[35],q[42];

// red K4 (16, 18, 28, 30) -> avoid connecting to all four
ccx q[15],q[17],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[27],q[43];
ccx q[15],q[17],q[42];

// red K4 (16, 18, 28, 38) -> avoid connecting to all four
ccx q[15],q[17],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[27],q[43];
ccx q[15],q[17],q[42];

// red K4 (16, 18, 28, 41) -> avoid connecting to all four
ccx q[15],q[17],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[27],q[43];
ccx q[15],q[17],q[42];

// red K4 (16, 18, 30, 32) -> avoid connecting to all four
ccx q[15],q[17],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[29],q[43];
ccx q[15],q[17],q[42];

// red K4 (16, 18, 32, 34) -> avoid connecting to all four
ccx q[15],q[17],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[31],q[43];
ccx q[15],q[17],q[42];

// red K4 (16, 18, 32, 39) -> avoid connecting to all four
ccx q[15],q[17],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[31],q[43];
ccx q[15],q[17],q[42];

// red K4 (16, 18, 34, 36) -> avoid connecting to all four
ccx q[15],q[17],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[33],q[43];
ccx q[15],q[17],q[42];

// red K4 (16, 18, 34, 41) -> avoid connecting to all four
ccx q[15],q[17],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[33],q[43];
ccx q[15],q[17],q[42];

// red K4 (16, 18, 36, 38) -> avoid connecting to all four
ccx q[15],q[17],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[35],q[43];
ccx q[15],q[17],q[42];

// red K4 (16, 18, 38, 39) -> avoid connecting to all four
ccx q[15],q[17],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[37],q[43];
ccx q[15],q[17],q[42];

// red K4 (16, 18, 39, 41) -> avoid connecting to all four
ccx q[15],q[17],q[42];
ccx q[42],q[38],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[38],q[43];
ccx q[15],q[17],q[42];

// red K4 (16, 23, 30, 37) -> avoid connecting to all four
ccx q[15],q[22],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[29],q[43];
ccx q[15],q[22],q[42];

// red K4 (16, 23, 36, 37) -> avoid connecting to all four
ccx q[15],q[22],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[35],q[43];
ccx q[15],q[22],q[42];

// red K4 (16, 23, 37, 39) -> avoid connecting to all four
ccx q[15],q[22],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[36],q[43];
ccx q[15],q[22],q[42];

// red K4 (16, 23, 39, 41) -> avoid connecting to all four
ccx q[15],q[22],q[42];
ccx q[42],q[38],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[38],q[43];
ccx q[15],q[22],q[42];

// red K4 (16, 26, 28, 38) -> avoid connecting to all four
ccx q[15],q[25],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[27],q[43];
ccx q[15],q[25],q[42];

// red K4 (16, 26, 36, 38) -> avoid connecting to all four
ccx q[15],q[25],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[35],q[43];
ccx q[15],q[25],q[42];

// red K4 (16, 26, 38, 39) -> avoid connecting to all four
ccx q[15],q[25],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[37],q[43];
ccx q[15],q[25],q[42];

// red K4 (16, 28, 29, 30) -> avoid connecting to all four
ccx q[15],q[27],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[28],q[43];
ccx q[15],q[27],q[42];

// red K4 (16, 28, 29, 41) -> avoid connecting to all four
ccx q[15],q[27],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[28],q[43];
ccx q[15],q[27],q[42];

// red K4 (16, 29, 39, 41) -> avoid connecting to all four
ccx q[15],q[28],q[42];
ccx q[42],q[38],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[38],q[43];
ccx q[15],q[28],q[42];

// red K4 (16, 36, 37, 38) -> avoid connecting to all four
ccx q[15],q[35],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[36],q[43];
ccx q[15],q[35],q[42];

// red K4 (16, 37, 38, 39) -> avoid connecting to all four
ccx q[15],q[36],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[37],q[43];
ccx q[15],q[36],q[42];

// red K4 (17, 18, 19, 31) -> avoid connecting to all four
ccx q[16],q[17],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[18],q[43];
ccx q[16],q[17],q[42];

// red K4 (17, 18, 19, 39) -> avoid connecting to all four
ccx q[16],q[17],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[18],q[43];
ccx q[16],q[17],q[42];

// red K4 (17, 18, 19, 40) -> avoid connecting to all four
ccx q[16],q[17],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[18],q[43];
ccx q[16],q[17],q[42];

// red K4 (17, 18, 30, 40) -> avoid connecting to all four
ccx q[16],q[17],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[29],q[43];
ccx q[16],q[17],q[42];

// red K4 (17, 18, 31, 38) -> avoid connecting to all four
ccx q[16],q[17],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[30],q[43];
ccx q[16],q[17],q[42];

// red K4 (17, 18, 38, 39) -> avoid connecting to all four
ccx q[16],q[17],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[37],q[43];
ccx q[16],q[17],q[42];

// red K4 (17, 18, 38, 40) -> avoid connecting to all four
ccx q[16],q[17],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[37],q[43];
ccx q[16],q[17],q[42];

// red K4 (17, 19, 29, 31) -> avoid connecting to all four
ccx q[16],q[18],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[28],q[43];
ccx q[16],q[18],q[42];

// red K4 (17, 19, 29, 39) -> avoid connecting to all four
ccx q[16],q[18],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[28],q[43];
ccx q[16],q[18],q[42];

// red K4 (17, 19, 29, 42) -> avoid connecting to all four
ccx q[16],q[18],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[28],q[43];
ccx q[16],q[18],q[42];

// red K4 (17, 19, 31, 33) -> avoid connecting to all four
ccx q[16],q[18],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[30],q[43];
ccx q[16],q[18],q[42];

// red K4 (17, 19, 33, 35) -> avoid connecting to all four
ccx q[16],q[18],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[32],q[43];
ccx q[16],q[18],q[42];

// red K4 (17, 19, 33, 40) -> avoid connecting to all four
ccx q[16],q[18],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[32],q[43];
ccx q[16],q[18],q[42];

// red K4 (17, 19, 35, 37) -> avoid connecting to all four
ccx q[16],q[18],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[34],q[43];
ccx q[16],q[18],q[42];

// red K4 (17, 19, 35, 42) -> avoid connecting to all four
ccx q[16],q[18],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[34],q[43];
ccx q[16],q[18],q[42];

// red K4 (17, 19, 37, 39) -> avoid connecting to all four
ccx q[16],q[18],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[36],q[43];
ccx q[16],q[18],q[42];

// red K4 (17, 19, 40, 42) -> avoid connecting to all four
ccx q[16],q[18],q[42];
ccx q[42],q[39],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[39],q[43];
ccx q[16],q[18],q[42];

// red K4 (17, 24, 31, 38) -> avoid connecting to all four
ccx q[16],q[23],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[30],q[43];
ccx q[16],q[23],q[42];

// red K4 (17, 24, 37, 38) -> avoid connecting to all four
ccx q[16],q[23],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[36],q[43];
ccx q[16],q[23],q[42];

// red K4 (17, 24, 38, 40) -> avoid connecting to all four
ccx q[16],q[23],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[37],q[43];
ccx q[16],q[23],q[42];

// red K4 (17, 24, 40, 42) -> avoid connecting to all four
ccx q[16],q[23],q[42];
ccx q[42],q[39],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[39],q[43];
ccx q[16],q[23],q[42];

// red K4 (17, 27, 29, 39) -> avoid connecting to all four
ccx q[16],q[26],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[28],q[43];
ccx q[16],q[26],q[42];

// red K4 (17, 27, 37, 39) -> avoid connecting to all four
ccx q[16],q[26],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[36],q[43];
ccx q[16],q[26],q[42];

// red K4 (17, 29, 30, 42) -> avoid connecting to all four
ccx q[16],q[28],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[29],q[43];
ccx q[16],q[28],q[42];

// red K4 (17, 30, 40, 42) -> avoid connecting to all four
ccx q[16],q[29],q[42];
ccx q[42],q[39],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[39],q[43];
ccx q[16],q[29],q[42];

// red K4 (17, 37, 38, 39) -> avoid connecting to all four
ccx q[16],q[36],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[37],q[43];
ccx q[16],q[36],q[42];

// red K4 (18, 19, 20, 32) -> avoid connecting to all four
ccx q[17],q[18],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[19],q[43];
ccx q[17],q[18],q[42];

// red K4 (18, 19, 20, 40) -> avoid connecting to all four
ccx q[17],q[18],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[19],q[43];
ccx q[17],q[18],q[42];

// red K4 (18, 19, 20, 41) -> avoid connecting to all four
ccx q[17],q[18],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[19],q[43];
ccx q[17],q[18],q[42];

// red K4 (18, 19, 31, 32) -> avoid connecting to all four
ccx q[17],q[18],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[30],q[43];
ccx q[17],q[18],q[42];

// red K4 (18, 19, 31, 41) -> avoid connecting to all four
ccx q[17],q[18],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[30],q[43];
ccx q[17],q[18],q[42];

// red K4 (18, 19, 32, 39) -> avoid connecting to all four
ccx q[17],q[18],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[31],q[43];
ccx q[17],q[18],q[42];

// red K4 (18, 19, 39, 41) -> avoid connecting to all four
ccx q[17],q[18],q[42];
ccx q[42],q[38],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[38],q[43];
ccx q[17],q[18],q[42];

// red K4 (18, 20, 30, 32) -> avoid connecting to all four
ccx q[17],q[19],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[29],q[43];
ccx q[17],q[19],q[42];

// red K4 (18, 20, 30, 40) -> avoid connecting to all four
ccx q[17],q[19],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[29],q[43];
ccx q[17],q[19],q[42];

// red K4 (18, 20, 32, 34) -> avoid connecting to all four
ccx q[17],q[19],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[31],q[43];
ccx q[17],q[19],q[42];

// red K4 (18, 20, 34, 36) -> avoid connecting to all four
ccx q[17],q[19],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[33],q[43];
ccx q[17],q[19],q[42];

// red K4 (18, 20, 34, 41) -> avoid connecting to all four
ccx q[17],q[19],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[33],q[43];
ccx q[17],q[19],q[42];

// red K4 (18, 20, 36, 38) -> avoid connecting to all four
ccx q[17],q[19],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[35],q[43];
ccx q[17],q[19],q[42];

// red K4 (18, 20, 38, 40) -> avoid connecting to all four
ccx q[17],q[19],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[37],q[43];
ccx q[17],q[19],q[42];

// red K4 (18, 25, 32, 39) -> avoid connecting to all four
ccx q[17],q[24],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[31],q[43];
ccx q[17],q[24],q[42];

// red K4 (18, 25, 38, 39) -> avoid connecting to all four
ccx q[17],q[24],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[37],q[43];
ccx q[17],q[24],q[42];

// red K4 (18, 25, 39, 41) -> avoid connecting to all four
ccx q[17],q[24],q[42];
ccx q[42],q[38],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[38],q[43];
ccx q[17],q[24],q[42];

// red K4 (18, 28, 30, 40) -> avoid connecting to all four
ccx q[17],q[27],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[29],q[43];
ccx q[17],q[27],q[42];

// red K4 (18, 28, 38, 40) -> avoid connecting to all four
ccx q[17],q[27],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[37],q[43];
ccx q[17],q[27],q[42];

// red K4 (19, 20, 21, 33) -> avoid connecting to all four
ccx q[18],q[19],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[20],q[43];
ccx q[18],q[19],q[42];

// red K4 (19, 20, 21, 41) -> avoid connecting to all four
ccx q[18],q[19],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[20],q[43];
ccx q[18],q[19],q[42];

// red K4 (19, 20, 21, 42) -> avoid connecting to all four
ccx q[18],q[19],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[20],q[43];
ccx q[18],q[19],q[42];

// red K4 (19, 20, 32, 33) -> avoid connecting to all four
ccx q[18],q[19],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[31],q[43];
ccx q[18],q[19],q[42];

// red K4 (19, 20, 32, 42) -> avoid connecting to all four
ccx q[18],q[19],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[31],q[43];
ccx q[18],q[19],q[42];

// red K4 (19, 20, 33, 40) -> avoid connecting to all four
ccx q[18],q[19],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[32],q[43];
ccx q[18],q[19],q[42];

// red K4 (19, 20, 40, 42) -> avoid connecting to all four
ccx q[18],q[19],q[42];
ccx q[42],q[39],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[39],q[43];
ccx q[18],q[19],q[42];

// red K4 (19, 21, 31, 33) -> avoid connecting to all four
ccx q[18],q[20],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[30],q[43];
ccx q[18],q[20],q[42];

// red K4 (19, 21, 31, 41) -> avoid connecting to all four
ccx q[18],q[20],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[30],q[43];
ccx q[18],q[20],q[42];

// red K4 (19, 21, 33, 35) -> avoid connecting to all four
ccx q[18],q[20],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[32],q[43];
ccx q[18],q[20],q[42];

// red K4 (19, 21, 35, 37) -> avoid connecting to all four
ccx q[18],q[20],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[34],q[43];
ccx q[18],q[20],q[42];

// red K4 (19, 21, 35, 42) -> avoid connecting to all four
ccx q[18],q[20],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[34],q[43];
ccx q[18],q[20],q[42];

// red K4 (19, 21, 37, 39) -> avoid connecting to all four
ccx q[18],q[20],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[36],q[43];
ccx q[18],q[20],q[42];

// red K4 (19, 21, 39, 41) -> avoid connecting to all four
ccx q[18],q[20],q[42];
ccx q[42],q[38],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[38],q[43];
ccx q[18],q[20],q[42];

// red K4 (19, 26, 33, 40) -> avoid connecting to all four
ccx q[18],q[25],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[32],q[43];
ccx q[18],q[25],q[42];

// red K4 (19, 26, 40, 42) -> avoid connecting to all four
ccx q[18],q[25],q[42];
ccx q[42],q[39],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[39],q[43];
ccx q[18],q[25],q[42];

// red K4 (19, 29, 31, 41) -> avoid connecting to all four
ccx q[18],q[28],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[30],q[43];
ccx q[18],q[28],q[42];

// red K4 (19, 29, 39, 41) -> avoid connecting to all four
ccx q[18],q[28],q[42];
ccx q[42],q[38],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[38],q[43];
ccx q[18],q[28],q[42];

// red K4 (19, 31, 32, 33) -> avoid connecting to all four
ccx q[18],q[30],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[31],q[43];
ccx q[18],q[30],q[42];

// red K4 (20, 21, 22, 34) -> avoid connecting to all four
ccx q[19],q[20],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[21],q[43];
ccx q[19],q[20],q[42];

// red K4 (20, 21, 22, 42) -> avoid connecting to all four
ccx q[19],q[20],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[21],q[43];
ccx q[19],q[20],q[42];

// red K4 (20, 21, 34, 41) -> avoid connecting to all four
ccx q[19],q[20],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[33],q[43];
ccx q[19],q[20],q[42];

// red K4 (20, 22, 32, 34) -> avoid connecting to all four
ccx q[19],q[21],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[31],q[43];
ccx q[19],q[21],q[42];

// red K4 (20, 22, 32, 42) -> avoid connecting to all four
ccx q[19],q[21],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[31],q[43];
ccx q[19],q[21],q[42];

// red K4 (20, 22, 34, 36) -> avoid connecting to all four
ccx q[19],q[21],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[33],q[43];
ccx q[19],q[21],q[42];

// red K4 (20, 22, 36, 38) -> avoid connecting to all four
ccx q[19],q[21],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[35],q[43];
ccx q[19],q[21],q[42];

// red K4 (20, 22, 38, 40) -> avoid connecting to all four
ccx q[19],q[21],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[37],q[43];
ccx q[19],q[21],q[42];

// red K4 (20, 22, 40, 42) -> avoid connecting to all four
ccx q[19],q[21],q[42];
ccx q[42],q[39],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[39],q[43];
ccx q[19],q[21],q[42];

// red K4 (20, 27, 34, 41) -> avoid connecting to all four
ccx q[19],q[26],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[33],q[43];
ccx q[19],q[26],q[42];

// red K4 (20, 30, 32, 42) -> avoid connecting to all four
ccx q[19],q[29],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[31],q[43];
ccx q[19],q[29],q[42];

// red K4 (20, 30, 40, 42) -> avoid connecting to all four
ccx q[19],q[29],q[42];
ccx q[42],q[39],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[39],q[43];
ccx q[19],q[29],q[42];

// red K4 (21, 22, 23, 35) -> avoid connecting to all four
ccx q[20],q[21],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[22],q[43];
ccx q[20],q[21],q[42];

// red K4 (21, 22, 34, 35) -> avoid connecting to all four
ccx q[20],q[21],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[33],q[43];
ccx q[20],q[21],q[42];

// red K4 (21, 22, 35, 42) -> avoid connecting to all four
ccx q[20],q[21],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[34],q[43];
ccx q[20],q[21],q[42];

// red K4 (21, 23, 33, 35) -> avoid connecting to all four
ccx q[20],q[22],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[32],q[43];
ccx q[20],q[22],q[42];

// red K4 (21, 23, 35, 37) -> avoid connecting to all four
ccx q[20],q[22],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[34],q[43];
ccx q[20],q[22],q[42];

// red K4 (21, 23, 37, 39) -> avoid connecting to all four
ccx q[20],q[22],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[36],q[43];
ccx q[20],q[22],q[42];

// red K4 (21, 23, 39, 41) -> avoid connecting to all four
ccx q[20],q[22],q[42];
ccx q[42],q[38],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[38],q[43];
ccx q[20],q[22],q[42];

// red K4 (21, 28, 35, 42) -> avoid connecting to all four
ccx q[20],q[27],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[34],q[43];
ccx q[20],q[27],q[42];

// red K4 (22, 23, 35, 36) -> avoid connecting to all four
ccx q[21],q[22],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[34],q[43];
ccx q[21],q[22],q[42];

// red K4 (22, 24, 34, 36) -> avoid connecting to all four
ccx q[21],q[23],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[33],q[43];
ccx q[21],q[23],q[42];

// red K4 (22, 24, 36, 38) -> avoid connecting to all four
ccx q[21],q[23],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[35],q[43];
ccx q[21],q[23],q[42];

// red K4 (22, 24, 38, 40) -> avoid connecting to all four
ccx q[21],q[23],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[37],q[43];
ccx q[21],q[23],q[42];

// red K4 (22, 24, 40, 42) -> avoid connecting to all four
ccx q[21],q[23],q[42];
ccx q[42],q[39],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[39],q[43];
ccx q[21],q[23],q[42];

// red K4 (22, 34, 35, 36) -> avoid connecting to all four
ccx q[21],q[33],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[34],q[43];
ccx q[21],q[33],q[42];

// red K4 (23, 25, 35, 37) -> avoid connecting to all four
ccx q[22],q[24],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[34],q[43];
ccx q[22],q[24],q[42];

// red K4 (23, 25, 37, 39) -> avoid connecting to all four
ccx q[22],q[24],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[36],q[43];
ccx q[22],q[24],q[42];

// red K4 (23, 25, 39, 41) -> avoid connecting to all four
ccx q[22],q[24],q[42];
ccx q[42],q[38],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[38],q[43];
ccx q[22],q[24],q[42];

// red K4 (23, 35, 36, 37) -> avoid connecting to all four
ccx q[22],q[34],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[35],q[43];
ccx q[22],q[34],q[42];

// red K4 (24, 26, 36, 38) -> avoid connecting to all four
ccx q[23],q[25],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[35],q[43];
ccx q[23],q[25],q[42];

// red K4 (24, 26, 38, 40) -> avoid connecting to all four
ccx q[23],q[25],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[37],q[43];
ccx q[23],q[25],q[42];

// red K4 (24, 26, 40, 42) -> avoid connecting to all four
ccx q[23],q[25],q[42];
ccx q[42],q[39],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[39],q[43];
ccx q[23],q[25],q[42];

// red K4 (24, 36, 37, 38) -> avoid connecting to all four
ccx q[23],q[35],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[36],q[43];
ccx q[23],q[35],q[42];

// red K4 (25, 26, 27, 39) -> avoid connecting to all four
ccx q[24],q[25],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[26],q[43];
ccx q[24],q[25],q[42];

// red K4 (25, 26, 38, 39) -> avoid connecting to all four
ccx q[24],q[25],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[37],q[43];
ccx q[24],q[25],q[42];

// red K4 (25, 27, 37, 39) -> avoid connecting to all four
ccx q[24],q[26],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[36],q[43];
ccx q[24],q[26],q[42];

// red K4 (25, 27, 39, 41) -> avoid connecting to all four
ccx q[24],q[26],q[42];
ccx q[42],q[38],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[38],q[43];
ccx q[24],q[26],q[42];

// red K4 (25, 37, 38, 39) -> avoid connecting to all four
ccx q[24],q[36],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[37],q[43];
ccx q[24],q[36],q[42];

// red K4 (26, 27, 28, 40) -> avoid connecting to all four
ccx q[25],q[26],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[27],q[43];
ccx q[25],q[26],q[42];

// red K4 (26, 28, 38, 40) -> avoid connecting to all four
ccx q[25],q[27],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[37],q[43];
ccx q[25],q[27],q[42];

// red K4 (26, 28, 40, 42) -> avoid connecting to all four
ccx q[25],q[27],q[42];
ccx q[42],q[39],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[39],q[43];
ccx q[25],q[27],q[42];

// red K4 (27, 28, 29, 41) -> avoid connecting to all four
ccx q[26],q[27],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[28],q[43];
ccx q[26],q[27],q[42];

// red K4 (27, 29, 39, 41) -> avoid connecting to all four
ccx q[26],q[28],q[42];
ccx q[42],q[38],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[38],q[43];
ccx q[26],q[28],q[42];

// red K4 (28, 29, 30, 42) -> avoid connecting to all four
ccx q[27],q[28],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[29],q[43];
ccx q[27],q[28],q[42];

// red K4 (28, 30, 40, 42) -> avoid connecting to all four
ccx q[27],q[29],q[42];
ccx q[42],q[39],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[39],q[43];
ccx q[27],q[29],q[42];

// blue K4 (1, 4, 5, 9) -> avoid connecting to none of the four
x q[0];
x q[3];
x q[4];
x q[8];
ccx q[0],q[3],q[42];
ccx q[42],q[4],q[43];
ccx q[43],q[8],q[44];
rz(0.15) q[44];
ccx q[43],q[8],q[44];
ccx q[42],q[4],q[43];
ccx q[0],q[3],q[42];
x q[0];
x q[3];
x q[4];
x q[8];

// blue K4 (1, 4, 5, 10) -> avoid connecting to none of the four
x q[0];
x q[3];
x q[4];
x q[9];
ccx q[0],q[3],q[42];
ccx q[42],q[4],q[43];
ccx q[43],q[9],q[44];
rz(0.15) q[44];
ccx q[43],q[9],q[44];
ccx q[42],q[4],q[43];
ccx q[0],q[3],q[42];
x q[0];
x q[3];
x q[4];
x q[9];

// blue K4 (1, 4, 5, 39) -> avoid connecting to none of the four
x q[0];
x q[3];
x q[4];
x q[38];
ccx q[0],q[3],q[42];
ccx q[42],q[4],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[4],q[43];
ccx q[0],q[3],q[42];
x q[0];
x q[3];
x q[4];
x q[38];

// blue K4 (1, 4, 7, 10) -> avoid connecting to none of the four
x q[0];
x q[3];
x q[6];
x q[9];
ccx q[0],q[3],q[42];
ccx q[42],q[6],q[43];
ccx q[43],q[9],q[44];
rz(0.15) q[44];
ccx q[43],q[9],q[44];
ccx q[42],q[6],q[43];
ccx q[0],q[3],q[42];
x q[0];
x q[3];
x q[6];
x q[9];

// blue K4 (1, 4, 7, 12) -> avoid connecting to none of the four
x q[0];
x q[3];
x q[6];
x q[11];
ccx q[0],q[3],q[42];
ccx q[42],q[6],q[43];
ccx q[43],q[11],q[44];
rz(0.15) q[44];
ccx q[43],q[11],q[44];
ccx q[42],q[6],q[43];
ccx q[0],q[3],q[42];
x q[0];
x q[3];
x q[6];
x q[11];

// blue K4 (1, 4, 7, 39) -> avoid connecting to none of the four
x q[0];
x q[3];
x q[6];
x q[38];
ccx q[0],q[3],q[42];
ccx q[42],q[6],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[6],q[43];
ccx q[0],q[3],q[42];
x q[0];
x q[3];
x q[6];
x q[38];

// blue K4 (1, 4, 7, 41) -> avoid connecting to none of the four
x q[0];
x q[3];
x q[6];
x q[40];
ccx q[0],q[3],q[42];
ccx q[42],q[6],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[6],q[43];
ccx q[0],q[3],q[42];
x q[0];
x q[3];
x q[6];
x q[40];

// blue K4 (1, 4, 9, 12) -> avoid connecting to none of the four
x q[0];
x q[3];
x q[8];
x q[11];
ccx q[0],q[3],q[42];
ccx q[42],q[8],q[43];
ccx q[43],q[11],q[44];
rz(0.15) q[44];
ccx q[43],q[11],q[44];
ccx q[42],q[8],q[43];
ccx q[0],q[3],q[42];
x q[0];
x q[3];
x q[8];
x q[11];

// blue K4 (1, 4, 9, 41) -> avoid connecting to none of the four
x q[0];
x q[3];
x q[8];
x q[40];
ccx q[0],q[3],q[42];
ccx q[42],q[8],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[8],q[43];
ccx q[0],q[3],q[42];
x q[0];
x q[3];
x q[8];
x q[40];

// blue K4 (1, 4, 10, 36) -> avoid connecting to none of the four
x q[0];
x q[3];
x q[9];
x q[35];
ccx q[0],q[3],q[42];
ccx q[42],q[9],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[9],q[43];
ccx q[0],q[3],q[42];
x q[0];
x q[3];
x q[9];
x q[35];

// blue K4 (1, 4, 10, 38) -> avoid connecting to none of the four
x q[0];
x q[3];
x q[9];
x q[37];
ccx q[0],q[3],q[42];
ccx q[42],q[9],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[9],q[43];
ccx q[0],q[3],q[42];
x q[0];
x q[3];
x q[9];
x q[37];

// blue K4 (1, 4, 12, 36) -> avoid connecting to none of the four
x q[0];
x q[3];
x q[11];
x q[35];
ccx q[0],q[3],q[42];
ccx q[42],q[11],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[11],q[43];
ccx q[0],q[3],q[42];
x q[0];
x q[3];
x q[11];
x q[35];

// blue K4 (1, 4, 12, 38) -> avoid connecting to none of the four
x q[0];
x q[3];
x q[11];
x q[37];
ccx q[0],q[3],q[42];
ccx q[42],q[11],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[11],q[43];
ccx q[0],q[3],q[42];
x q[0];
x q[3];
x q[11];
x q[37];

// blue K4 (1, 4, 36, 39) -> avoid connecting to none of the four
x q[0];
x q[3];
x q[35];
x q[38];
ccx q[0],q[3],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[35],q[43];
ccx q[0],q[3],q[42];
x q[0];
x q[3];
x q[35];
x q[38];

// blue K4 (1, 4, 36, 41) -> avoid connecting to none of the four
x q[0];
x q[3];
x q[35];
x q[40];
ccx q[0],q[3],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[35],q[43];
ccx q[0],q[3],q[42];
x q[0];
x q[3];
x q[35];
x q[40];

// blue K4 (1, 4, 38, 41) -> avoid connecting to none of the four
x q[0];
x q[3];
x q[37];
x q[40];
ccx q[0],q[3],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[37],q[43];
ccx q[0],q[3],q[42];
x q[0];
x q[3];
x q[37];
x q[40];

// blue K4 (1, 5, 6, 9) -> avoid connecting to none of the four
x q[0];
x q[4];
x q[5];
x q[8];
ccx q[0],q[4],q[42];
ccx q[42],q[5],q[43];
ccx q[43],q[8],q[44];
rz(0.15) q[44];
ccx q[43],q[8],q[44];
ccx q[42],q[5],q[43];
ccx q[0],q[4],q[42];
x q[0];
x q[4];
x q[5];
x q[8];

// blue K4 (1, 5, 6, 10) -> avoid connecting to none of the four
x q[0];
x q[4];
x q[5];
x q[9];
ccx q[0],q[4],q[42];
ccx q[42],q[5],q[43];
ccx q[43],q[9],q[44];
rz(0.15) q[44];
ccx q[43],q[9],q[44];
ccx q[42],q[5],q[43];
ccx q[0],q[4],q[42];
x q[0];
x q[4];
x q[5];
x q[9];

// blue K4 (1, 5, 6, 40) -> avoid connecting to none of the four
x q[0];
x q[4];
x q[5];
x q[39];
ccx q[0],q[4],q[42];
ccx q[42],q[5],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[5],q[43];
ccx q[0],q[4],q[42];
x q[0];
x q[4];
x q[5];
x q[39];

// blue K4 (1, 5, 9, 20) -> avoid connecting to none of the four
x q[0];
x q[4];
x q[8];
x q[19];
ccx q[0],q[4],q[42];
ccx q[42],q[8],q[43];
ccx q[43],q[19],q[44];
rz(0.15) q[44];
ccx q[43],q[19],q[44];
ccx q[42],q[8],q[43];
ccx q[0],q[4],q[42];
x q[0];
x q[4];
x q[8];
x q[19];

// blue K4 (1, 5, 9, 33) -> avoid connecting to none of the four
x q[0];
x q[4];
x q[8];
x q[32];
ccx q[0],q[4],q[42];
ccx q[42],q[8],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[8],q[43];
ccx q[0],q[4],q[42];
x q[0];
x q[4];
x q[8];
x q[32];

// blue K4 (1, 5, 10, 16) -> avoid connecting to none of the four
x q[0];
x q[4];
x q[9];
x q[15];
ccx q[0],q[4],q[42];
ccx q[42],q[9],q[43];
ccx q[43],q[15],q[44];
rz(0.15) q[44];
ccx q[43],q[15],q[44];
ccx q[42],q[9],q[43];
ccx q[0],q[4],q[42];
x q[0];
x q[4];
x q[9];
x q[15];

// blue K4 (1, 5, 10, 29) -> avoid connecting to none of the four
x q[0];
x q[4];
x q[9];
x q[28];
ccx q[0],q[4],q[42];
ccx q[42],q[9],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[9],q[43];
ccx q[0],q[4],q[42];
x q[0];
x q[4];
x q[9];
x q[28];

// blue K4 (1, 5, 16, 20) -> avoid connecting to none of the four
x q[0];
x q[4];
x q[15];
x q[19];
ccx q[0],q[4],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[19],q[44];
rz(0.15) q[44];
ccx q[43],q[19],q[44];
ccx q[42],q[15],q[43];
ccx q[0],q[4],q[42];
x q[0];
x q[4];
x q[15];
x q[19];

// blue K4 (1, 5, 16, 33) -> avoid connecting to none of the four
x q[0];
x q[4];
x q[15];
x q[32];
ccx q[0],q[4],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[15],q[43];
ccx q[0],q[4],q[42];
x q[0];
x q[4];
x q[15];
x q[32];

// blue K4 (1, 5, 16, 40) -> avoid connecting to none of the four
x q[0];
x q[4];
x q[15];
x q[39];
ccx q[0],q[4],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[15],q[43];
ccx q[0],q[4],q[42];
x q[0];
x q[4];
x q[15];
x q[39];

// blue K4 (1, 5, 20, 29) -> avoid connecting to none of the four
x q[0];
x q[4];
x q[19];
x q[28];
ccx q[0],q[4],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[19],q[43];
ccx q[0],q[4],q[42];
x q[0];
x q[4];
x q[19];
x q[28];

// blue K4 (1, 5, 20, 39) -> avoid connecting to none of the four
x q[0];
x q[4];
x q[19];
x q[38];
ccx q[0],q[4],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[19],q[43];
ccx q[0],q[4],q[42];
x q[0];
x q[4];
x q[19];
x q[38];

// blue K4 (1, 5, 29, 33) -> avoid connecting to none of the four
x q[0];
x q[4];
x q[28];
x q[32];
ccx q[0],q[4],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[28],q[43];
ccx q[0],q[4],q[42];
x q[0];
x q[4];
x q[28];
x q[32];

// blue K4 (1, 5, 29, 40) -> avoid connecting to none of the four
x q[0];
x q[4];
x q[28];
x q[39];
ccx q[0],q[4],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[28],q[43];
ccx q[0],q[4],q[42];
x q[0];
x q[4];
x q[28];
x q[39];

// blue K4 (1, 5, 33, 39) -> avoid connecting to none of the four
x q[0];
x q[4];
x q[32];
x q[38];
ccx q[0],q[4],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[32],q[43];
ccx q[0],q[4],q[42];
x q[0];
x q[4];
x q[32];
x q[38];

// blue K4 (1, 5, 39, 40) -> avoid connecting to none of the four
x q[0];
x q[4];
x q[38];
x q[39];
ccx q[0],q[4],q[42];
ccx q[42],q[38],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[38],q[43];
ccx q[0],q[4],q[42];
x q[0];
x q[4];
x q[38];
x q[39];

// blue K4 (1, 6, 7, 10) -> avoid connecting to none of the four
x q[0];
x q[5];
x q[6];
x q[9];
ccx q[0],q[5],q[42];
ccx q[42],q[6],q[43];
ccx q[43],q[9],q[44];
rz(0.15) q[44];
ccx q[43],q[9],q[44];
ccx q[42],q[6],q[43];
ccx q[0],q[5],q[42];
x q[0];
x q[5];
x q[6];
x q[9];

// blue K4 (1, 6, 7, 12) -> avoid connecting to none of the four
x q[0];
x q[5];
x q[6];
x q[11];
ccx q[0],q[5],q[42];
ccx q[42],q[6],q[43];
ccx q[43],q[11],q[44];
rz(0.15) q[44];
ccx q[43],q[11],q[44];
ccx q[42],q[6],q[43];
ccx q[0],q[5],q[42];
x q[0];
x q[5];
x q[6];
x q[11];

// blue K4 (1, 6, 7, 41) -> avoid connecting to none of the four
x q[0];
x q[5];
x q[6];
x q[40];
ccx q[0],q[5],q[42];
ccx q[42],q[6],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[6],q[43];
ccx q[0],q[5],q[42];
x q[0];
x q[5];
x q[6];
x q[40];

// blue K4 (1, 6, 9, 12) -> avoid connecting to none of the four
x q[0];
x q[5];
x q[8];
x q[11];
ccx q[0],q[5],q[42];
ccx q[42],q[8],q[43];
ccx q[43],q[11],q[44];
rz(0.15) q[44];
ccx q[43],q[11],q[44];
ccx q[42],q[8],q[43];
ccx q[0],q[5],q[42];
x q[0];
x q[5];
x q[8];
x q[11];

// blue K4 (1, 6, 9, 41) -> avoid connecting to none of the four
x q[0];
x q[5];
x q[8];
x q[40];
ccx q[0],q[5],q[42];
ccx q[42],q[8],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[8],q[43];
ccx q[0],q[5],q[42];
x q[0];
x q[5];
x q[8];
x q[40];

// blue K4 (1, 6, 10, 25) -> avoid connecting to none of the four
x q[0];
x q[5];
x q[9];
x q[24];
ccx q[0],q[5],q[42];
ccx q[42],q[9],q[43];
ccx q[43],q[24],q[44];
rz(0.15) q[44];
ccx q[43],q[24],q[44];
ccx q[42],q[9],q[43];
ccx q[0],q[5],q[42];
x q[0];
x q[5];
x q[9];
x q[24];

// blue K4 (1, 6, 10, 38) -> avoid connecting to none of the four
x q[0];
x q[5];
x q[9];
x q[37];
ccx q[0],q[5],q[42];
ccx q[42],q[9],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[9],q[43];
ccx q[0],q[5],q[42];
x q[0];
x q[5];
x q[9];
x q[37];

// blue K4 (1, 6, 12, 38) -> avoid connecting to none of the four
x q[0];
x q[5];
x q[11];
x q[37];
ccx q[0],q[5],q[42];
ccx q[42],q[11],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[11],q[43];
ccx q[0],q[5],q[42];
x q[0];
x q[5];
x q[11];
x q[37];

// blue K4 (1, 6, 12, 40) -> avoid connecting to none of the four
x q[0];
x q[5];
x q[11];
x q[39];
ccx q[0],q[5],q[42];
ccx q[42],q[11],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[11],q[43];
ccx q[0],q[5],q[42];
x q[0];
x q[5];
x q[11];
x q[39];

// blue K4 (1, 6, 25, 40) -> avoid connecting to none of the four
x q[0];
x q[5];
x q[24];
x q[39];
ccx q[0],q[5],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[24],q[43];
ccx q[0],q[5],q[42];
x q[0];
x q[5];
x q[24];
x q[39];

// blue K4 (1, 6, 38, 41) -> avoid connecting to none of the four
x q[0];
x q[5];
x q[37];
x q[40];
ccx q[0],q[5],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[37],q[43];
ccx q[0],q[5],q[42];
x q[0];
x q[5];
x q[37];
x q[40];

// blue K4 (1, 6, 40, 41) -> avoid connecting to none of the four
x q[0];
x q[5];
x q[39];
x q[40];
ccx q[0],q[5],q[42];
ccx q[42],q[39],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[39],q[43];
ccx q[0],q[5],q[42];
x q[0];
x q[5];
x q[39];
x q[40];

// blue K4 (1, 7, 10, 16) -> avoid connecting to none of the four
x q[0];
x q[6];
x q[9];
x q[15];
ccx q[0],q[6],q[42];
ccx q[42],q[9],q[43];
ccx q[43],q[15],q[44];
rz(0.15) q[44];
ccx q[43],q[15],q[44];
ccx q[42],q[9],q[43];
ccx q[0],q[6],q[42];
x q[0];
x q[6];
x q[9];
x q[15];

// blue K4 (1, 7, 10, 18) -> avoid connecting to none of the four
x q[0];
x q[6];
x q[9];
x q[17];
ccx q[0],q[6],q[42];
ccx q[42],q[9],q[43];
ccx q[43],q[17],q[44];
rz(0.15) q[44];
ccx q[43],q[17],q[44];
ccx q[42],q[9],q[43];
ccx q[0],q[6],q[42];
x q[0];
x q[6];
x q[9];
x q[17];

// blue K4 (1, 7, 12, 16) -> avoid connecting to none of the four
x q[0];
x q[6];
x q[11];
x q[15];
ccx q[0],q[6],q[42];
ccx q[42],q[11],q[43];
ccx q[43],q[15],q[44];
rz(0.15) q[44];
ccx q[43],q[15],q[44];
ccx q[42],q[11],q[43];
ccx q[0],q[6],q[42];
x q[0];
x q[6];
x q[11];
x q[15];

// blue K4 (1, 7, 12, 18) -> avoid connecting to none of the four
x q[0];
x q[6];
x q[11];
x q[17];
ccx q[0],q[6],q[42];
ccx q[42],q[11],q[43];
ccx q[43],q[17],q[44];
rz(0.15) q[44];
ccx q[43],q[17],q[44];
ccx q[42],q[11],q[43];
ccx q[0],q[6],q[42];
x q[0];
x q[6];
x q[11];
x q[17];

// blue K4 (1, 7, 16, 33) -> avoid connecting to none of the four
x q[0];
x q[6];
x q[15];
x q[32];
ccx q[0],q[6],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[15],q[43];
ccx q[0],q[6],q[42];
x q[0];
x q[6];
x q[15];
x q[32];

// blue K4 (1, 7, 16, 35) -> avoid connecting to none of the four
x q[0];
x q[6];
x q[15];
x q[34];
ccx q[0],q[6],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[15],q[43];
ccx q[0],q[6],q[42];
x q[0];
x q[6];
x q[15];
x q[34];

// blue K4 (1, 7, 18, 33) -> avoid connecting to none of the four
x q[0];
x q[6];
x q[17];
x q[32];
ccx q[0],q[6],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[17],q[43];
ccx q[0],q[6],q[42];
x q[0];
x q[6];
x q[17];
x q[32];

// blue K4 (1, 7, 18, 35) -> avoid connecting to none of the four
x q[0];
x q[6];
x q[17];
x q[34];
ccx q[0],q[6],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[17],q[43];
ccx q[0],q[6],q[42];
x q[0];
x q[6];
x q[17];
x q[34];

// blue K4 (1, 7, 33, 39) -> avoid connecting to none of the four
x q[0];
x q[6];
x q[32];
x q[38];
ccx q[0],q[6],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[32],q[43];
ccx q[0],q[6],q[42];
x q[0];
x q[6];
x q[32];
x q[38];

// blue K4 (1, 7, 33, 41) -> avoid connecting to none of the four
x q[0];
x q[6];
x q[32];
x q[40];
ccx q[0],q[6],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[32],q[43];
ccx q[0],q[6],q[42];
x q[0];
x q[6];
x q[32];
x q[40];

// blue K4 (1, 7, 35, 39) -> avoid connecting to none of the four
x q[0];
x q[6];
x q[34];
x q[38];
ccx q[0],q[6],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[34],q[43];
ccx q[0],q[6],q[42];
x q[0];
x q[6];
x q[34];
x q[38];

// blue K4 (1, 7, 35, 41) -> avoid connecting to none of the four
x q[0];
x q[6];
x q[34];
x q[40];
ccx q[0],q[6],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[34],q[43];
ccx q[0],q[6],q[42];
x q[0];
x q[6];
x q[34];
x q[40];

// blue K4 (1, 9, 12, 18) -> avoid connecting to none of the four
x q[0];
x q[8];
x q[11];
x q[17];
ccx q[0],q[8],q[42];
ccx q[42],q[11],q[43];
ccx q[43],q[17],q[44];
rz(0.15) q[44];
ccx q[43],q[17],q[44];
ccx q[42],q[11],q[43];
ccx q[0],q[8],q[42];
x q[0];
x q[8];
x q[11];
x q[17];

// blue K4 (1, 9, 12, 20) -> avoid connecting to none of the four
x q[0];
x q[8];
x q[11];
x q[19];
ccx q[0],q[8],q[42];
ccx q[42],q[11],q[43];
ccx q[43],q[19],q[44];
rz(0.15) q[44];
ccx q[43],q[19],q[44];
ccx q[42],q[11],q[43];
ccx q[0],q[8],q[42];
x q[0];
x q[8];
x q[11];
x q[19];

// blue K4 (1, 9, 18, 33) -> avoid connecting to none of the four
x q[0];
x q[8];
x q[17];
x q[32];
ccx q[0],q[8],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[17],q[43];
ccx q[0],q[8],q[42];
x q[0];
x q[8];
x q[17];
x q[32];

// blue K4 (1, 9, 18, 35) -> avoid connecting to none of the four
x q[0];
x q[8];
x q[17];
x q[34];
ccx q[0],q[8],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[17],q[43];
ccx q[0],q[8],q[42];
x q[0];
x q[8];
x q[17];
x q[34];

// blue K4 (1, 9, 20, 35) -> avoid connecting to none of the four
x q[0];
x q[8];
x q[19];
x q[34];
ccx q[0],q[8],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[19],q[43];
ccx q[0],q[8],q[42];
x q[0];
x q[8];
x q[19];
x q[34];

// blue K4 (1, 9, 33, 41) -> avoid connecting to none of the four
x q[0];
x q[8];
x q[32];
x q[40];
ccx q[0],q[8],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[32],q[43];
ccx q[0],q[8],q[42];
x q[0];
x q[8];
x q[32];
x q[40];

// blue K4 (1, 9, 35, 41) -> avoid connecting to none of the four
x q[0];
x q[8];
x q[34];
x q[40];
ccx q[0],q[8],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[34],q[43];
ccx q[0],q[8],q[42];
x q[0];
x q[8];
x q[34];
x q[40];

// blue K4 (1, 10, 16, 25) -> avoid connecting to none of the four
x q[0];
x q[9];
x q[15];
x q[24];
ccx q[0],q[9],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[24],q[44];
rz(0.15) q[44];
ccx q[43],q[24],q[44];
ccx q[42],q[15],q[43];
ccx q[0],q[9],q[42];
x q[0];
x q[9];
x q[15];
x q[24];

// blue K4 (1, 10, 16, 27) -> avoid connecting to none of the four
x q[0];
x q[9];
x q[15];
x q[26];
ccx q[0],q[9],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[26],q[44];
rz(0.15) q[44];
ccx q[43],q[26],q[44];
ccx q[42],q[15],q[43];
ccx q[0],q[9],q[42];
x q[0];
x q[9];
x q[15];
x q[26];

// blue K4 (1, 10, 18, 27) -> avoid connecting to none of the four
x q[0];
x q[9];
x q[17];
x q[26];
ccx q[0],q[9],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[26],q[44];
rz(0.15) q[44];
ccx q[43],q[26],q[44];
ccx q[42],q[17],q[43];
ccx q[0],q[9],q[42];
x q[0];
x q[9];
x q[17];
x q[26];

// blue K4 (1, 10, 18, 29) -> avoid connecting to none of the four
x q[0];
x q[9];
x q[17];
x q[28];
ccx q[0],q[9],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[17],q[43];
ccx q[0],q[9],q[42];
x q[0];
x q[9];
x q[17];
x q[28];

// blue K4 (1, 10, 25, 29) -> avoid connecting to none of the four
x q[0];
x q[9];
x q[24];
x q[28];
ccx q[0],q[9],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[24],q[43];
ccx q[0],q[9],q[42];
x q[0];
x q[9];
x q[24];
x q[28];

// blue K4 (1, 10, 25, 36) -> avoid connecting to none of the four
x q[0];
x q[9];
x q[24];
x q[35];
ccx q[0],q[9],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[24],q[43];
ccx q[0],q[9],q[42];
x q[0];
x q[9];
x q[24];
x q[35];

// blue K4 (1, 10, 27, 36) -> avoid connecting to none of the four
x q[0];
x q[9];
x q[26];
x q[35];
ccx q[0],q[9],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[26],q[43];
ccx q[0],q[9],q[42];
x q[0];
x q[9];
x q[26];
x q[35];

// blue K4 (1, 10, 27, 38) -> avoid connecting to none of the four
x q[0];
x q[9];
x q[26];
x q[37];
ccx q[0],q[9],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[26],q[43];
ccx q[0],q[9],q[42];
x q[0];
x q[9];
x q[26];
x q[37];

// blue K4 (1, 10, 29, 38) -> avoid connecting to none of the four
x q[0];
x q[9];
x q[28];
x q[37];
ccx q[0],q[9],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[28],q[43];
ccx q[0],q[9],q[42];
x q[0];
x q[9];
x q[28];
x q[37];

// blue K4 (1, 12, 16, 20) -> avoid connecting to none of the four
x q[0];
x q[11];
x q[15];
x q[19];
ccx q[0],q[11],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[19],q[44];
rz(0.15) q[44];
ccx q[43],q[19],q[44];
ccx q[42],q[15],q[43];
ccx q[0],q[11],q[42];
x q[0];
x q[11];
x q[15];
x q[19];

// blue K4 (1, 12, 16, 27) -> avoid connecting to none of the four
x q[0];
x q[11];
x q[15];
x q[26];
ccx q[0],q[11],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[26],q[44];
rz(0.15) q[44];
ccx q[43],q[26],q[44];
ccx q[42],q[15],q[43];
ccx q[0],q[11],q[42];
x q[0];
x q[11];
x q[15];
x q[26];

// blue K4 (1, 12, 16, 40) -> avoid connecting to none of the four
x q[0];
x q[11];
x q[15];
x q[39];
ccx q[0],q[11],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[15],q[43];
ccx q[0],q[11],q[42];
x q[0];
x q[11];
x q[15];
x q[39];

// blue K4 (1, 12, 18, 27) -> avoid connecting to none of the four
x q[0];
x q[11];
x q[17];
x q[26];
ccx q[0],q[11],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[26],q[44];
rz(0.15) q[44];
ccx q[43],q[26],q[44];
ccx q[42],q[17],q[43];
ccx q[0],q[11],q[42];
x q[0];
x q[11];
x q[17];
x q[26];

// blue K4 (1, 12, 18, 29) -> avoid connecting to none of the four
x q[0];
x q[11];
x q[17];
x q[28];
ccx q[0],q[11],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[17],q[43];
ccx q[0],q[11],q[42];
x q[0];
x q[11];
x q[17];
x q[28];

// blue K4 (1, 12, 20, 29) -> avoid connecting to none of the four
x q[0];
x q[11];
x q[19];
x q[28];
ccx q[0],q[11],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[19],q[43];
ccx q[0],q[11],q[42];
x q[0];
x q[11];
x q[19];
x q[28];

// blue K4 (1, 12, 27, 36) -> avoid connecting to none of the four
x q[0];
x q[11];
x q[26];
x q[35];
ccx q[0],q[11],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[26],q[43];
ccx q[0],q[11],q[42];
x q[0];
x q[11];
x q[26];
x q[35];

// blue K4 (1, 12, 27, 38) -> avoid connecting to none of the four
x q[0];
x q[11];
x q[26];
x q[37];
ccx q[0],q[11],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[26],q[43];
ccx q[0],q[11],q[42];
x q[0];
x q[11];
x q[26];
x q[37];

// blue K4 (1, 12, 29, 38) -> avoid connecting to none of the four
x q[0];
x q[11];
x q[28];
x q[37];
ccx q[0],q[11],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[28],q[43];
ccx q[0],q[11],q[42];
x q[0];
x q[11];
x q[28];
x q[37];

// blue K4 (1, 12, 29, 40) -> avoid connecting to none of the four
x q[0];
x q[11];
x q[28];
x q[39];
ccx q[0],q[11],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[28],q[43];
ccx q[0],q[11],q[42];
x q[0];
x q[11];
x q[28];
x q[39];

// blue K4 (1, 12, 36, 40) -> avoid connecting to none of the four
x q[0];
x q[11];
x q[35];
x q[39];
ccx q[0],q[11],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[35],q[43];
ccx q[0],q[11],q[42];
x q[0];
x q[11];
x q[35];
x q[39];

// blue K4 (1, 16, 20, 25) -> avoid connecting to none of the four
x q[0];
x q[15];
x q[19];
x q[24];
ccx q[0],q[15],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[24],q[44];
rz(0.15) q[44];
ccx q[43],q[24],q[44];
ccx q[42],q[19],q[43];
ccx q[0],q[15],q[42];
x q[0];
x q[15];
x q[19];
x q[24];

// blue K4 (1, 16, 20, 35) -> avoid connecting to none of the four
x q[0];
x q[15];
x q[19];
x q[34];
ccx q[0],q[15],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[19],q[43];
ccx q[0],q[15],q[42];
x q[0];
x q[15];
x q[19];
x q[34];

// blue K4 (1, 16, 25, 33) -> avoid connecting to none of the four
x q[0];
x q[15];
x q[24];
x q[32];
ccx q[0],q[15],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[24],q[43];
ccx q[0],q[15],q[42];
x q[0];
x q[15];
x q[24];
x q[32];

// blue K4 (1, 16, 25, 40) -> avoid connecting to none of the four
x q[0];
x q[15];
x q[24];
x q[39];
ccx q[0],q[15],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[24],q[43];
ccx q[0],q[15],q[42];
x q[0];
x q[15];
x q[24];
x q[39];

// blue K4 (1, 16, 27, 33) -> avoid connecting to none of the four
x q[0];
x q[15];
x q[26];
x q[32];
ccx q[0],q[15],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[26],q[43];
ccx q[0],q[15],q[42];
x q[0];
x q[15];
x q[26];
x q[32];

// blue K4 (1, 16, 27, 35) -> avoid connecting to none of the four
x q[0];
x q[15];
x q[26];
x q[34];
ccx q[0],q[15],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[26],q[43];
ccx q[0],q[15],q[42];
x q[0];
x q[15];
x q[26];
x q[34];

// blue K4 (1, 16, 35, 40) -> avoid connecting to none of the four
x q[0];
x q[15];
x q[34];
x q[39];
ccx q[0],q[15],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[34],q[43];
ccx q[0],q[15],q[42];
x q[0];
x q[15];
x q[34];
x q[39];

// blue K4 (1, 18, 27, 33) -> avoid connecting to none of the four
x q[0];
x q[17];
x q[26];
x q[32];
ccx q[0],q[17],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[26],q[43];
ccx q[0],q[17],q[42];
x q[0];
x q[17];
x q[26];
x q[32];

// blue K4 (1, 18, 27, 35) -> avoid connecting to none of the four
x q[0];
x q[17];
x q[26];
x q[34];
ccx q[0],q[17],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[26],q[43];
ccx q[0],q[17],q[42];
x q[0];
x q[17];
x q[26];
x q[34];

// blue K4 (1, 18, 29, 33) -> avoid connecting to none of the four
x q[0];
x q[17];
x q[28];
x q[32];
ccx q[0],q[17],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[28],q[43];
ccx q[0],q[17],q[42];
x q[0];
x q[17];
x q[28];
x q[32];

// blue K4 (1, 18, 29, 35) -> avoid connecting to none of the four
x q[0];
x q[17];
x q[28];
x q[34];
ccx q[0],q[17],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[28],q[43];
ccx q[0],q[17],q[42];
x q[0];
x q[17];
x q[28];
x q[34];

// blue K4 (1, 20, 25, 29) -> avoid connecting to none of the four
x q[0];
x q[19];
x q[24];
x q[28];
ccx q[0],q[19],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[24],q[43];
ccx q[0],q[19],q[42];
x q[0];
x q[19];
x q[24];
x q[28];

// blue K4 (1, 20, 29, 35) -> avoid connecting to none of the four
x q[0];
x q[19];
x q[28];
x q[34];
ccx q[0],q[19],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[28],q[43];
ccx q[0],q[19],q[42];
x q[0];
x q[19];
x q[28];
x q[34];

// blue K4 (1, 20, 35, 39) -> avoid connecting to none of the four
x q[0];
x q[19];
x q[34];
x q[38];
ccx q[0],q[19],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[34],q[43];
ccx q[0],q[19],q[42];
x q[0];
x q[19];
x q[34];
x q[38];

// blue K4 (1, 25, 29, 33) -> avoid connecting to none of the four
x q[0];
x q[24];
x q[28];
x q[32];
ccx q[0],q[24],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[28],q[43];
ccx q[0],q[24],q[42];
x q[0];
x q[24];
x q[28];
x q[32];

// blue K4 (1, 25, 29, 40) -> avoid connecting to none of the four
x q[0];
x q[24];
x q[28];
x q[39];
ccx q[0],q[24],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[28],q[43];
ccx q[0],q[24],q[42];
x q[0];
x q[24];
x q[28];
x q[39];

// blue K4 (1, 25, 33, 36) -> avoid connecting to none of the four
x q[0];
x q[24];
x q[32];
x q[35];
ccx q[0],q[24],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[32],q[43];
ccx q[0],q[24],q[42];
x q[0];
x q[24];
x q[32];
x q[35];

// blue K4 (1, 25, 36, 40) -> avoid connecting to none of the four
x q[0];
x q[24];
x q[35];
x q[39];
ccx q[0],q[24],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[35],q[43];
ccx q[0],q[24],q[42];
x q[0];
x q[24];
x q[35];
x q[39];

// blue K4 (1, 27, 33, 36) -> avoid connecting to none of the four
x q[0];
x q[26];
x q[32];
x q[35];
ccx q[0],q[26],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[32],q[43];
ccx q[0],q[26],q[42];
x q[0];
x q[26];
x q[32];
x q[35];

// blue K4 (1, 27, 33, 38) -> avoid connecting to none of the four
x q[0];
x q[26];
x q[32];
x q[37];
ccx q[0],q[26],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[32],q[43];
ccx q[0],q[26],q[42];
x q[0];
x q[26];
x q[32];
x q[37];

// blue K4 (1, 27, 35, 38) -> avoid connecting to none of the four
x q[0];
x q[26];
x q[34];
x q[37];
ccx q[0],q[26],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[34],q[43];
ccx q[0],q[26],q[42];
x q[0];
x q[26];
x q[34];
x q[37];

// blue K4 (1, 29, 33, 38) -> avoid connecting to none of the four
x q[0];
x q[28];
x q[32];
x q[37];
ccx q[0],q[28],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[32],q[43];
ccx q[0],q[28],q[42];
x q[0];
x q[28];
x q[32];
x q[37];

// blue K4 (1, 29, 35, 38) -> avoid connecting to none of the four
x q[0];
x q[28];
x q[34];
x q[37];
ccx q[0],q[28],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[34],q[43];
ccx q[0],q[28],q[42];
x q[0];
x q[28];
x q[34];
x q[37];

// blue K4 (1, 29, 35, 40) -> avoid connecting to none of the four
x q[0];
x q[28];
x q[34];
x q[39];
ccx q[0],q[28],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[34],q[43];
ccx q[0],q[28],q[42];
x q[0];
x q[28];
x q[34];
x q[39];

// blue K4 (1, 33, 36, 39) -> avoid connecting to none of the four
x q[0];
x q[32];
x q[35];
x q[38];
ccx q[0],q[32],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[35],q[43];
ccx q[0],q[32],q[42];
x q[0];
x q[32];
x q[35];
x q[38];

// blue K4 (1, 33, 36, 41) -> avoid connecting to none of the four
x q[0];
x q[32];
x q[35];
x q[40];
ccx q[0],q[32],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[35],q[43];
ccx q[0],q[32],q[42];
x q[0];
x q[32];
x q[35];
x q[40];

// blue K4 (1, 33, 38, 41) -> avoid connecting to none of the four
x q[0];
x q[32];
x q[37];
x q[40];
ccx q[0],q[32],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[37],q[43];
ccx q[0],q[32],q[42];
x q[0];
x q[32];
x q[37];
x q[40];

// blue K4 (1, 35, 38, 41) -> avoid connecting to none of the four
x q[0];
x q[34];
x q[37];
x q[40];
ccx q[0],q[34],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[37],q[43];
ccx q[0],q[34],q[42];
x q[0];
x q[34];
x q[37];
x q[40];

// blue K4 (1, 35, 39, 40) -> avoid connecting to none of the four
x q[0];
x q[34];
x q[38];
x q[39];
ccx q[0],q[34],q[42];
ccx q[42],q[38],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[38],q[43];
ccx q[0],q[34],q[42];
x q[0];
x q[34];
x q[38];
x q[39];

// blue K4 (1, 35, 40, 41) -> avoid connecting to none of the four
x q[0];
x q[34];
x q[39];
x q[40];
ccx q[0],q[34],q[42];
ccx q[42],q[39],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[39],q[43];
ccx q[0],q[34],q[42];
x q[0];
x q[34];
x q[39];
x q[40];

// blue K4 (1, 36, 39, 40) -> avoid connecting to none of the four
x q[0];
x q[35];
x q[38];
x q[39];
ccx q[0],q[35],q[42];
ccx q[42],q[38],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[38],q[43];
ccx q[0],q[35],q[42];
x q[0];
x q[35];
x q[38];
x q[39];

// blue K4 (1, 36, 40, 41) -> avoid connecting to none of the four
x q[0];
x q[35];
x q[39];
x q[40];
ccx q[0],q[35],q[42];
ccx q[42],q[39],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[39],q[43];
ccx q[0],q[35],q[42];
x q[0];
x q[35];
x q[39];
x q[40];

// blue K4 (2, 5, 6, 10) -> avoid connecting to none of the four
x q[1];
x q[4];
x q[5];
x q[9];
ccx q[1],q[4],q[42];
ccx q[42],q[5],q[43];
ccx q[43],q[9],q[44];
rz(0.15) q[44];
ccx q[43],q[9],q[44];
ccx q[42],q[5],q[43];
ccx q[1],q[4],q[42];
x q[1];
x q[4];
x q[5];
x q[9];

// blue K4 (2, 5, 6, 11) -> avoid connecting to none of the four
x q[1];
x q[4];
x q[5];
x q[10];
ccx q[1],q[4],q[42];
ccx q[42],q[5],q[43];
ccx q[43],q[10],q[44];
rz(0.15) q[44];
ccx q[43],q[10],q[44];
ccx q[42],q[5],q[43];
ccx q[1],q[4],q[42];
x q[1];
x q[4];
x q[5];
x q[10];

// blue K4 (2, 5, 6, 40) -> avoid connecting to none of the four
x q[1];
x q[4];
x q[5];
x q[39];
ccx q[1],q[4],q[42];
ccx q[42],q[5],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[5],q[43];
ccx q[1],q[4],q[42];
x q[1];
x q[4];
x q[5];
x q[39];

// blue K4 (2, 5, 8, 11) -> avoid connecting to none of the four
x q[1];
x q[4];
x q[7];
x q[10];
ccx q[1],q[4],q[42];
ccx q[42],q[7],q[43];
ccx q[43],q[10],q[44];
rz(0.15) q[44];
ccx q[43],q[10],q[44];
ccx q[42],q[7],q[43];
ccx q[1],q[4],q[42];
x q[1];
x q[4];
x q[7];
x q[10];

// blue K4 (2, 5, 8, 13) -> avoid connecting to none of the four
x q[1];
x q[4];
x q[7];
x q[12];
ccx q[1],q[4],q[42];
ccx q[42],q[7],q[43];
ccx q[43],q[12],q[44];
rz(0.15) q[44];
ccx q[43],q[12],q[44];
ccx q[42],q[7],q[43];
ccx q[1],q[4],q[42];
x q[1];
x q[4];
x q[7];
x q[12];

// blue K4 (2, 5, 8, 40) -> avoid connecting to none of the four
x q[1];
x q[4];
x q[7];
x q[39];
ccx q[1],q[4],q[42];
ccx q[42],q[7],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[7],q[43];
ccx q[1],q[4],q[42];
x q[1];
x q[4];
x q[7];
x q[39];

// blue K4 (2, 5, 8, 42) -> avoid connecting to none of the four
x q[1];
x q[4];
x q[7];
x q[41];
ccx q[1],q[4],q[42];
ccx q[42],q[7],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[7],q[43];
ccx q[1],q[4],q[42];
x q[1];
x q[4];
x q[7];
x q[41];

// blue K4 (2, 5, 10, 13) -> avoid connecting to none of the four
x q[1];
x q[4];
x q[9];
x q[12];
ccx q[1],q[4],q[42];
ccx q[42],q[9],q[43];
ccx q[43],q[12],q[44];
rz(0.15) q[44];
ccx q[43],q[12],q[44];
ccx q[42],q[9],q[43];
ccx q[1],q[4],q[42];
x q[1];
x q[4];
x q[9];
x q[12];

// blue K4 (2, 5, 10, 42) -> avoid connecting to none of the four
x q[1];
x q[4];
x q[9];
x q[41];
ccx q[1],q[4],q[42];
ccx q[42],q[9],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[9],q[43];
ccx q[1],q[4],q[42];
x q[1];
x q[4];
x q[9];
x q[41];

// blue K4 (2, 5, 11, 37) -> avoid connecting to none of the four
x q[1];
x q[4];
x q[10];
x q[36];
ccx q[1],q[4],q[42];
ccx q[42],q[10],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[10],q[43];
ccx q[1],q[4],q[42];
x q[1];
x q[4];
x q[10];
x q[36];

// blue K4 (2, 5, 11, 39) -> avoid connecting to none of the four
x q[1];
x q[4];
x q[10];
x q[38];
ccx q[1],q[4],q[42];
ccx q[42],q[10],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[10],q[43];
ccx q[1],q[4],q[42];
x q[1];
x q[4];
x q[10];
x q[38];

// blue K4 (2, 5, 13, 37) -> avoid connecting to none of the four
x q[1];
x q[4];
x q[12];
x q[36];
ccx q[1],q[4],q[42];
ccx q[42],q[12],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[12],q[43];
ccx q[1],q[4],q[42];
x q[1];
x q[4];
x q[12];
x q[36];

// blue K4 (2, 5, 13, 39) -> avoid connecting to none of the four
x q[1];
x q[4];
x q[12];
x q[38];
ccx q[1],q[4],q[42];
ccx q[42],q[12],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[12],q[43];
ccx q[1],q[4],q[42];
x q[1];
x q[4];
x q[12];
x q[38];

// blue K4 (2, 5, 37, 40) -> avoid connecting to none of the four
x q[1];
x q[4];
x q[36];
x q[39];
ccx q[1],q[4],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[36],q[43];
ccx q[1],q[4],q[42];
x q[1];
x q[4];
x q[36];
x q[39];

// blue K4 (2, 5, 37, 42) -> avoid connecting to none of the four
x q[1];
x q[4];
x q[36];
x q[41];
ccx q[1],q[4],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[36],q[43];
ccx q[1],q[4],q[42];
x q[1];
x q[4];
x q[36];
x q[41];

// blue K4 (2, 5, 39, 40) -> avoid connecting to none of the four
x q[1];
x q[4];
x q[38];
x q[39];
ccx q[1],q[4],q[42];
ccx q[42],q[38],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[38],q[43];
ccx q[1],q[4],q[42];
x q[1];
x q[4];
x q[38];
x q[39];

// blue K4 (2, 5, 39, 42) -> avoid connecting to none of the four
x q[1];
x q[4];
x q[38];
x q[41];
ccx q[1],q[4],q[42];
ccx q[42],q[38],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[38],q[43];
ccx q[1],q[4],q[42];
x q[1];
x q[4];
x q[38];
x q[41];

// blue K4 (2, 6, 7, 10) -> avoid connecting to none of the four
x q[1];
x q[5];
x q[6];
x q[9];
ccx q[1],q[5],q[42];
ccx q[42],q[6],q[43];
ccx q[43],q[9],q[44];
rz(0.15) q[44];
ccx q[43],q[9],q[44];
ccx q[42],q[6],q[43];
ccx q[1],q[5],q[42];
x q[1];
x q[5];
x q[6];
x q[9];

// blue K4 (2, 6, 7, 11) -> avoid connecting to none of the four
x q[1];
x q[5];
x q[6];
x q[10];
ccx q[1],q[5],q[42];
ccx q[42],q[6],q[43];
ccx q[43],q[10],q[44];
rz(0.15) q[44];
ccx q[43],q[10],q[44];
ccx q[42],q[6],q[43];
ccx q[1],q[5],q[42];
x q[1];
x q[5];
x q[6];
x q[10];

// blue K4 (2, 6, 7, 41) -> avoid connecting to none of the four
x q[1];
x q[5];
x q[6];
x q[40];
ccx q[1],q[5],q[42];
ccx q[42],q[6],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[6],q[43];
ccx q[1],q[5],q[42];
x q[1];
x q[5];
x q[6];
x q[40];

// blue K4 (2, 6, 10, 21) -> avoid connecting to none of the four
x q[1];
x q[5];
x q[9];
x q[20];
ccx q[1],q[5],q[42];
ccx q[42],q[9],q[43];
ccx q[43],q[20],q[44];
rz(0.15) q[44];
ccx q[43],q[20],q[44];
ccx q[42],q[9],q[43];
ccx q[1],q[5],q[42];
x q[1];
x q[5];
x q[9];
x q[20];

// blue K4 (2, 6, 10, 34) -> avoid connecting to none of the four
x q[1];
x q[5];
x q[9];
x q[33];
ccx q[1],q[5],q[42];
ccx q[42],q[9],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[9],q[43];
ccx q[1],q[5],q[42];
x q[1];
x q[5];
x q[9];
x q[33];

// blue K4 (2, 6, 11, 17) -> avoid connecting to none of the four
x q[1];
x q[5];
x q[10];
x q[16];
ccx q[1],q[5],q[42];
ccx q[42],q[10],q[43];
ccx q[43],q[16],q[44];
rz(0.15) q[44];
ccx q[43],q[16],q[44];
ccx q[42],q[10],q[43];
ccx q[1],q[5],q[42];
x q[1];
x q[5];
x q[10];
x q[16];

// blue K4 (2, 6, 11, 30) -> avoid connecting to none of the four
x q[1];
x q[5];
x q[10];
x q[29];
ccx q[1],q[5],q[42];
ccx q[42],q[10],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[10],q[43];
ccx q[1],q[5],q[42];
x q[1];
x q[5];
x q[10];
x q[29];

// blue K4 (2, 6, 17, 21) -> avoid connecting to none of the four
x q[1];
x q[5];
x q[16];
x q[20];
ccx q[1],q[5],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[20],q[44];
rz(0.15) q[44];
ccx q[43],q[20],q[44];
ccx q[42],q[16],q[43];
ccx q[1],q[5],q[42];
x q[1];
x q[5];
x q[16];
x q[20];

// blue K4 (2, 6, 17, 34) -> avoid connecting to none of the four
x q[1];
x q[5];
x q[16];
x q[33];
ccx q[1],q[5],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[16],q[43];
ccx q[1],q[5],q[42];
x q[1];
x q[5];
x q[16];
x q[33];

// blue K4 (2, 6, 17, 41) -> avoid connecting to none of the four
x q[1];
x q[5];
x q[16];
x q[40];
ccx q[1],q[5],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[16],q[43];
ccx q[1],q[5],q[42];
x q[1];
x q[5];
x q[16];
x q[40];

// blue K4 (2, 6, 21, 30) -> avoid connecting to none of the four
x q[1];
x q[5];
x q[20];
x q[29];
ccx q[1],q[5],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[20],q[43];
ccx q[1],q[5],q[42];
x q[1];
x q[5];
x q[20];
x q[29];

// blue K4 (2, 6, 21, 40) -> avoid connecting to none of the four
x q[1];
x q[5];
x q[20];
x q[39];
ccx q[1],q[5],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[20],q[43];
ccx q[1],q[5],q[42];
x q[1];
x q[5];
x q[20];
x q[39];

// blue K4 (2, 6, 30, 34) -> avoid connecting to none of the four
x q[1];
x q[5];
x q[29];
x q[33];
ccx q[1],q[5],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[29],q[43];
ccx q[1],q[5],q[42];
x q[1];
x q[5];
x q[29];
x q[33];

// blue K4 (2, 6, 30, 41) -> avoid connecting to none of the four
x q[1];
x q[5];
x q[29];
x q[40];
ccx q[1],q[5],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[29],q[43];
ccx q[1],q[5],q[42];
x q[1];
x q[5];
x q[29];
x q[40];

// blue K4 (2, 6, 34, 40) -> avoid connecting to none of the four
x q[1];
x q[5];
x q[33];
x q[39];
ccx q[1],q[5],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[33],q[43];
ccx q[1],q[5],q[42];
x q[1];
x q[5];
x q[33];
x q[39];

// blue K4 (2, 6, 40, 41) -> avoid connecting to none of the four
x q[1];
x q[5];
x q[39];
x q[40];
ccx q[1],q[5],q[42];
ccx q[42],q[39],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[39],q[43];
ccx q[1],q[5],q[42];
x q[1];
x q[5];
x q[39];
x q[40];

// blue K4 (2, 7, 8, 11) -> avoid connecting to none of the four
x q[1];
x q[6];
x q[7];
x q[10];
ccx q[1],q[6],q[42];
ccx q[42],q[7],q[43];
ccx q[43],q[10],q[44];
rz(0.15) q[44];
ccx q[43],q[10],q[44];
ccx q[42],q[7],q[43];
ccx q[1],q[6],q[42];
x q[1];
x q[6];
x q[7];
x q[10];

// blue K4 (2, 7, 8, 13) -> avoid connecting to none of the four
x q[1];
x q[6];
x q[7];
x q[12];
ccx q[1],q[6],q[42];
ccx q[42],q[7],q[43];
ccx q[43],q[12],q[44];
rz(0.15) q[44];
ccx q[43],q[12],q[44];
ccx q[42],q[7],q[43];
ccx q[1],q[6],q[42];
x q[1];
x q[6];
x q[7];
x q[12];

// blue K4 (2, 7, 8, 42) -> avoid connecting to none of the four
x q[1];
x q[6];
x q[7];
x q[41];
ccx q[1],q[6],q[42];
ccx q[42],q[7],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[7],q[43];
ccx q[1],q[6],q[42];
x q[1];
x q[6];
x q[7];
x q[41];

// blue K4 (2, 7, 10, 13) -> avoid connecting to none of the four
x q[1];
x q[6];
x q[9];
x q[12];
ccx q[1],q[6],q[42];
ccx q[42],q[9],q[43];
ccx q[43],q[12],q[44];
rz(0.15) q[44];
ccx q[43],q[12],q[44];
ccx q[42],q[9],q[43];
ccx q[1],q[6],q[42];
x q[1];
x q[6];
x q[9];
x q[12];

// blue K4 (2, 7, 10, 42) -> avoid connecting to none of the four
x q[1];
x q[6];
x q[9];
x q[41];
ccx q[1],q[6],q[42];
ccx q[42],q[9],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[9],q[43];
ccx q[1],q[6],q[42];
x q[1];
x q[6];
x q[9];
x q[41];

// blue K4 (2, 7, 11, 26) -> avoid connecting to none of the four
x q[1];
x q[6];
x q[10];
x q[25];
ccx q[1],q[6],q[42];
ccx q[42],q[10],q[43];
ccx q[43],q[25],q[44];
rz(0.15) q[44];
ccx q[43],q[25],q[44];
ccx q[42],q[10],q[43];
ccx q[1],q[6],q[42];
x q[1];
x q[6];
x q[10];
x q[25];

// blue K4 (2, 7, 11, 39) -> avoid connecting to none of the four
x q[1];
x q[6];
x q[10];
x q[38];
ccx q[1],q[6],q[42];
ccx q[42],q[10],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[10],q[43];
ccx q[1],q[6],q[42];
x q[1];
x q[6];
x q[10];
x q[38];

// blue K4 (2, 7, 13, 39) -> avoid connecting to none of the four
x q[1];
x q[6];
x q[12];
x q[38];
ccx q[1],q[6],q[42];
ccx q[42],q[12],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[12],q[43];
ccx q[1],q[6],q[42];
x q[1];
x q[6];
x q[12];
x q[38];

// blue K4 (2, 7, 13, 41) -> avoid connecting to none of the four
x q[1];
x q[6];
x q[12];
x q[40];
ccx q[1],q[6],q[42];
ccx q[42],q[12],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[12],q[43];
ccx q[1],q[6],q[42];
x q[1];
x q[6];
x q[12];
x q[40];

// blue K4 (2, 7, 26, 41) -> avoid connecting to none of the four
x q[1];
x q[6];
x q[25];
x q[40];
ccx q[1],q[6],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[25],q[43];
ccx q[1],q[6],q[42];
x q[1];
x q[6];
x q[25];
x q[40];

// blue K4 (2, 7, 39, 42) -> avoid connecting to none of the four
x q[1];
x q[6];
x q[38];
x q[41];
ccx q[1],q[6],q[42];
ccx q[42],q[38],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[38],q[43];
ccx q[1],q[6],q[42];
x q[1];
x q[6];
x q[38];
x q[41];

// blue K4 (2, 7, 41, 42) -> avoid connecting to none of the four
x q[1];
x q[6];
x q[40];
x q[41];
ccx q[1],q[6],q[42];
ccx q[42],q[40],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[40],q[43];
ccx q[1],q[6],q[42];
x q[1];
x q[6];
x q[40];
x q[41];

// blue K4 (2, 8, 11, 17) -> avoid connecting to none of the four
x q[1];
x q[7];
x q[10];
x q[16];
ccx q[1],q[7],q[42];
ccx q[42],q[10],q[43];
ccx q[43],q[16],q[44];
rz(0.15) q[44];
ccx q[43],q[16],q[44];
ccx q[42],q[10],q[43];
ccx q[1],q[7],q[42];
x q[1];
x q[7];
x q[10];
x q[16];

// blue K4 (2, 8, 11, 19) -> avoid connecting to none of the four
x q[1];
x q[7];
x q[10];
x q[18];
ccx q[1],q[7],q[42];
ccx q[42],q[10],q[43];
ccx q[43],q[18],q[44];
rz(0.15) q[44];
ccx q[43],q[18],q[44];
ccx q[42],q[10],q[43];
ccx q[1],q[7],q[42];
x q[1];
x q[7];
x q[10];
x q[18];

// blue K4 (2, 8, 13, 17) -> avoid connecting to none of the four
x q[1];
x q[7];
x q[12];
x q[16];
ccx q[1],q[7],q[42];
ccx q[42],q[12],q[43];
ccx q[43],q[16],q[44];
rz(0.15) q[44];
ccx q[43],q[16],q[44];
ccx q[42],q[12],q[43];
ccx q[1],q[7],q[42];
x q[1];
x q[7];
x q[12];
x q[16];

// blue K4 (2, 8, 13, 19) -> avoid connecting to none of the four
x q[1];
x q[7];
x q[12];
x q[18];
ccx q[1],q[7],q[42];
ccx q[42],q[12],q[43];
ccx q[43],q[18],q[44];
rz(0.15) q[44];
ccx q[43],q[18],q[44];
ccx q[42],q[12],q[43];
ccx q[1],q[7],q[42];
x q[1];
x q[7];
x q[12];
x q[18];

// blue K4 (2, 8, 17, 34) -> avoid connecting to none of the four
x q[1];
x q[7];
x q[16];
x q[33];
ccx q[1],q[7],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[16],q[43];
ccx q[1],q[7],q[42];
x q[1];
x q[7];
x q[16];
x q[33];

// blue K4 (2, 8, 17, 36) -> avoid connecting to none of the four
x q[1];
x q[7];
x q[16];
x q[35];
ccx q[1],q[7],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[16],q[43];
ccx q[1],q[7],q[42];
x q[1];
x q[7];
x q[16];
x q[35];

// blue K4 (2, 8, 19, 34) -> avoid connecting to none of the four
x q[1];
x q[7];
x q[18];
x q[33];
ccx q[1],q[7],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[18],q[43];
ccx q[1],q[7],q[42];
x q[1];
x q[7];
x q[18];
x q[33];

// blue K4 (2, 8, 19, 36) -> avoid connecting to none of the four
x q[1];
x q[7];
x q[18];
x q[35];
ccx q[1],q[7],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[18],q[43];
ccx q[1],q[7],q[42];
x q[1];
x q[7];
x q[18];
x q[35];

// blue K4 (2, 8, 34, 40) -> avoid connecting to none of the four
x q[1];
x q[7];
x q[33];
x q[39];
ccx q[1],q[7],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[33],q[43];
ccx q[1],q[7],q[42];
x q[1];
x q[7];
x q[33];
x q[39];

// blue K4 (2, 8, 34, 42) -> avoid connecting to none of the four
x q[1];
x q[7];
x q[33];
x q[41];
ccx q[1],q[7],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[33],q[43];
ccx q[1],q[7],q[42];
x q[1];
x q[7];
x q[33];
x q[41];

// blue K4 (2, 8, 36, 40) -> avoid connecting to none of the four
x q[1];
x q[7];
x q[35];
x q[39];
ccx q[1],q[7],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[35],q[43];
ccx q[1],q[7],q[42];
x q[1];
x q[7];
x q[35];
x q[39];

// blue K4 (2, 8, 36, 42) -> avoid connecting to none of the four
x q[1];
x q[7];
x q[35];
x q[41];
ccx q[1],q[7],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[35],q[43];
ccx q[1],q[7],q[42];
x q[1];
x q[7];
x q[35];
x q[41];

// blue K4 (2, 10, 13, 19) -> avoid connecting to none of the four
x q[1];
x q[9];
x q[12];
x q[18];
ccx q[1],q[9],q[42];
ccx q[42],q[12],q[43];
ccx q[43],q[18],q[44];
rz(0.15) q[44];
ccx q[43],q[18],q[44];
ccx q[42],q[12],q[43];
ccx q[1],q[9],q[42];
x q[1];
x q[9];
x q[12];
x q[18];

// blue K4 (2, 10, 13, 21) -> avoid connecting to none of the four
x q[1];
x q[9];
x q[12];
x q[20];
ccx q[1],q[9],q[42];
ccx q[42],q[12],q[43];
ccx q[43],q[20],q[44];
rz(0.15) q[44];
ccx q[43],q[20],q[44];
ccx q[42],q[12],q[43];
ccx q[1],q[9],q[42];
x q[1];
x q[9];
x q[12];
x q[20];

// blue K4 (2, 10, 19, 34) -> avoid connecting to none of the four
x q[1];
x q[9];
x q[18];
x q[33];
ccx q[1],q[9],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[18],q[43];
ccx q[1],q[9],q[42];
x q[1];
x q[9];
x q[18];
x q[33];

// blue K4 (2, 10, 19, 36) -> avoid connecting to none of the four
x q[1];
x q[9];
x q[18];
x q[35];
ccx q[1],q[9],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[18],q[43];
ccx q[1],q[9],q[42];
x q[1];
x q[9];
x q[18];
x q[35];

// blue K4 (2, 10, 21, 36) -> avoid connecting to none of the four
x q[1];
x q[9];
x q[20];
x q[35];
ccx q[1],q[9],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[20],q[43];
ccx q[1],q[9],q[42];
x q[1];
x q[9];
x q[20];
x q[35];

// blue K4 (2, 10, 34, 42) -> avoid connecting to none of the four
x q[1];
x q[9];
x q[33];
x q[41];
ccx q[1],q[9],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[33],q[43];
ccx q[1],q[9],q[42];
x q[1];
x q[9];
x q[33];
x q[41];

// blue K4 (2, 10, 36, 42) -> avoid connecting to none of the four
x q[1];
x q[9];
x q[35];
x q[41];
ccx q[1],q[9],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[35],q[43];
ccx q[1],q[9],q[42];
x q[1];
x q[9];
x q[35];
x q[41];

// blue K4 (2, 11, 17, 26) -> avoid connecting to none of the four
x q[1];
x q[10];
x q[16];
x q[25];
ccx q[1],q[10],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[25],q[44];
rz(0.15) q[44];
ccx q[43],q[25],q[44];
ccx q[42],q[16],q[43];
ccx q[1],q[10],q[42];
x q[1];
x q[10];
x q[16];
x q[25];

// blue K4 (2, 11, 17, 28) -> avoid connecting to none of the four
x q[1];
x q[10];
x q[16];
x q[27];
ccx q[1],q[10],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[16],q[43];
ccx q[1],q[10],q[42];
x q[1];
x q[10];
x q[16];
x q[27];

// blue K4 (2, 11, 19, 28) -> avoid connecting to none of the four
x q[1];
x q[10];
x q[18];
x q[27];
ccx q[1],q[10],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[18],q[43];
ccx q[1],q[10],q[42];
x q[1];
x q[10];
x q[18];
x q[27];

// blue K4 (2, 11, 19, 30) -> avoid connecting to none of the four
x q[1];
x q[10];
x q[18];
x q[29];
ccx q[1],q[10],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[18],q[43];
ccx q[1],q[10],q[42];
x q[1];
x q[10];
x q[18];
x q[29];

// blue K4 (2, 11, 26, 30) -> avoid connecting to none of the four
x q[1];
x q[10];
x q[25];
x q[29];
ccx q[1],q[10],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[25],q[43];
ccx q[1],q[10],q[42];
x q[1];
x q[10];
x q[25];
x q[29];

// blue K4 (2, 11, 26, 37) -> avoid connecting to none of the four
x q[1];
x q[10];
x q[25];
x q[36];
ccx q[1],q[10],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[25],q[43];
ccx q[1],q[10],q[42];
x q[1];
x q[10];
x q[25];
x q[36];

// blue K4 (2, 11, 28, 37) -> avoid connecting to none of the four
x q[1];
x q[10];
x q[27];
x q[36];
ccx q[1],q[10],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[27],q[43];
ccx q[1],q[10],q[42];
x q[1];
x q[10];
x q[27];
x q[36];

// blue K4 (2, 11, 28, 39) -> avoid connecting to none of the four
x q[1];
x q[10];
x q[27];
x q[38];
ccx q[1],q[10],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[27],q[43];
ccx q[1],q[10],q[42];
x q[1];
x q[10];
x q[27];
x q[38];

// blue K4 (2, 11, 30, 39) -> avoid connecting to none of the four
x q[1];
x q[10];
x q[29];
x q[38];
ccx q[1],q[10],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[29],q[43];
ccx q[1],q[10],q[42];
x q[1];
x q[10];
x q[29];
x q[38];

// blue K4 (2, 13, 17, 21) -> avoid connecting to none of the four
x q[1];
x q[12];
x q[16];
x q[20];
ccx q[1],q[12],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[20],q[44];
rz(0.15) q[44];
ccx q[43],q[20],q[44];
ccx q[42],q[16],q[43];
ccx q[1],q[12],q[42];
x q[1];
x q[12];
x q[16];
x q[20];

// blue K4 (2, 13, 17, 28) -> avoid connecting to none of the four
x q[1];
x q[12];
x q[16];
x q[27];
ccx q[1],q[12],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[16],q[43];
ccx q[1],q[12],q[42];
x q[1];
x q[12];
x q[16];
x q[27];

// blue K4 (2, 13, 17, 41) -> avoid connecting to none of the four
x q[1];
x q[12];
x q[16];
x q[40];
ccx q[1],q[12],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[16],q[43];
ccx q[1],q[12],q[42];
x q[1];
x q[12];
x q[16];
x q[40];

// blue K4 (2, 13, 19, 28) -> avoid connecting to none of the four
x q[1];
x q[12];
x q[18];
x q[27];
ccx q[1],q[12],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[18],q[43];
ccx q[1],q[12],q[42];
x q[1];
x q[12];
x q[18];
x q[27];

// blue K4 (2, 13, 19, 30) -> avoid connecting to none of the four
x q[1];
x q[12];
x q[18];
x q[29];
ccx q[1],q[12],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[18],q[43];
ccx q[1],q[12],q[42];
x q[1];
x q[12];
x q[18];
x q[29];

// blue K4 (2, 13, 21, 30) -> avoid connecting to none of the four
x q[1];
x q[12];
x q[20];
x q[29];
ccx q[1],q[12],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[20],q[43];
ccx q[1],q[12],q[42];
x q[1];
x q[12];
x q[20];
x q[29];

// blue K4 (2, 13, 28, 37) -> avoid connecting to none of the four
x q[1];
x q[12];
x q[27];
x q[36];
ccx q[1],q[12],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[27],q[43];
ccx q[1],q[12],q[42];
x q[1];
x q[12];
x q[27];
x q[36];

// blue K4 (2, 13, 28, 39) -> avoid connecting to none of the four
x q[1];
x q[12];
x q[27];
x q[38];
ccx q[1],q[12],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[27],q[43];
ccx q[1],q[12],q[42];
x q[1];
x q[12];
x q[27];
x q[38];

// blue K4 (2, 13, 30, 39) -> avoid connecting to none of the four
x q[1];
x q[12];
x q[29];
x q[38];
ccx q[1],q[12],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[29],q[43];
ccx q[1],q[12],q[42];
x q[1];
x q[12];
x q[29];
x q[38];

// blue K4 (2, 13, 30, 41) -> avoid connecting to none of the four
x q[1];
x q[12];
x q[29];
x q[40];
ccx q[1],q[12],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[29],q[43];
ccx q[1],q[12],q[42];
x q[1];
x q[12];
x q[29];
x q[40];

// blue K4 (2, 13, 37, 41) -> avoid connecting to none of the four
x q[1];
x q[12];
x q[36];
x q[40];
ccx q[1],q[12],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[36],q[43];
ccx q[1],q[12],q[42];
x q[1];
x q[12];
x q[36];
x q[40];

// blue K4 (2, 17, 21, 26) -> avoid connecting to none of the four
x q[1];
x q[16];
x q[20];
x q[25];
ccx q[1],q[16],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[25],q[44];
rz(0.15) q[44];
ccx q[43],q[25],q[44];
ccx q[42],q[20],q[43];
ccx q[1],q[16],q[42];
x q[1];
x q[16];
x q[20];
x q[25];

// blue K4 (2, 17, 21, 36) -> avoid connecting to none of the four
x q[1];
x q[16];
x q[20];
x q[35];
ccx q[1],q[16],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[20],q[43];
ccx q[1],q[16],q[42];
x q[1];
x q[16];
x q[20];
x q[35];

// blue K4 (2, 17, 26, 34) -> avoid connecting to none of the four
x q[1];
x q[16];
x q[25];
x q[33];
ccx q[1],q[16],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[25],q[43];
ccx q[1],q[16],q[42];
x q[1];
x q[16];
x q[25];
x q[33];

// blue K4 (2, 17, 26, 41) -> avoid connecting to none of the four
x q[1];
x q[16];
x q[25];
x q[40];
ccx q[1],q[16],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[25],q[43];
ccx q[1],q[16],q[42];
x q[1];
x q[16];
x q[25];
x q[40];

// blue K4 (2, 17, 28, 34) -> avoid connecting to none of the four
x q[1];
x q[16];
x q[27];
x q[33];
ccx q[1],q[16],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[27],q[43];
ccx q[1],q[16],q[42];
x q[1];
x q[16];
x q[27];
x q[33];

// blue K4 (2, 17, 28, 36) -> avoid connecting to none of the four
x q[1];
x q[16];
x q[27];
x q[35];
ccx q[1],q[16],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[27],q[43];
ccx q[1],q[16],q[42];
x q[1];
x q[16];
x q[27];
x q[35];

// blue K4 (2, 17, 36, 41) -> avoid connecting to none of the four
x q[1];
x q[16];
x q[35];
x q[40];
ccx q[1],q[16],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[35],q[43];
ccx q[1],q[16],q[42];
x q[1];
x q[16];
x q[35];
x q[40];

// blue K4 (2, 19, 28, 34) -> avoid connecting to none of the four
x q[1];
x q[18];
x q[27];
x q[33];
ccx q[1],q[18],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[27],q[43];
ccx q[1],q[18],q[42];
x q[1];
x q[18];
x q[27];
x q[33];

// blue K4 (2, 19, 28, 36) -> avoid connecting to none of the four
x q[1];
x q[18];
x q[27];
x q[35];
ccx q[1],q[18],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[27],q[43];
ccx q[1],q[18],q[42];
x q[1];
x q[18];
x q[27];
x q[35];

// blue K4 (2, 19, 30, 34) -> avoid connecting to none of the four
x q[1];
x q[18];
x q[29];
x q[33];
ccx q[1],q[18],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[29],q[43];
ccx q[1],q[18],q[42];
x q[1];
x q[18];
x q[29];
x q[33];

// blue K4 (2, 19, 30, 36) -> avoid connecting to none of the four
x q[1];
x q[18];
x q[29];
x q[35];
ccx q[1],q[18],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[29],q[43];
ccx q[1],q[18],q[42];
x q[1];
x q[18];
x q[29];
x q[35];

// blue K4 (2, 21, 26, 30) -> avoid connecting to none of the four
x q[1];
x q[20];
x q[25];
x q[29];
ccx q[1],q[20],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[25],q[43];
ccx q[1],q[20],q[42];
x q[1];
x q[20];
x q[25];
x q[29];

// blue K4 (2, 21, 30, 36) -> avoid connecting to none of the four
x q[1];
x q[20];
x q[29];
x q[35];
ccx q[1],q[20],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[29],q[43];
ccx q[1],q[20],q[42];
x q[1];
x q[20];
x q[29];
x q[35];

// blue K4 (2, 21, 36, 40) -> avoid connecting to none of the four
x q[1];
x q[20];
x q[35];
x q[39];
ccx q[1],q[20],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[35],q[43];
ccx q[1],q[20],q[42];
x q[1];
x q[20];
x q[35];
x q[39];

// blue K4 (2, 26, 30, 34) -> avoid connecting to none of the four
x q[1];
x q[25];
x q[29];
x q[33];
ccx q[1],q[25],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[29],q[43];
ccx q[1],q[25],q[42];
x q[1];
x q[25];
x q[29];
x q[33];

// blue K4 (2, 26, 30, 41) -> avoid connecting to none of the four
x q[1];
x q[25];
x q[29];
x q[40];
ccx q[1],q[25],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[29],q[43];
ccx q[1],q[25],q[42];
x q[1];
x q[25];
x q[29];
x q[40];

// blue K4 (2, 26, 34, 37) -> avoid connecting to none of the four
x q[1];
x q[25];
x q[33];
x q[36];
ccx q[1],q[25],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[33],q[43];
ccx q[1],q[25],q[42];
x q[1];
x q[25];
x q[33];
x q[36];

// blue K4 (2, 26, 37, 41) -> avoid connecting to none of the four
x q[1];
x q[25];
x q[36];
x q[40];
ccx q[1],q[25],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[36],q[43];
ccx q[1],q[25],q[42];
x q[1];
x q[25];
x q[36];
x q[40];

// blue K4 (2, 28, 34, 37) -> avoid connecting to none of the four
x q[1];
x q[27];
x q[33];
x q[36];
ccx q[1],q[27],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[33],q[43];
ccx q[1],q[27],q[42];
x q[1];
x q[27];
x q[33];
x q[36];

// blue K4 (2, 28, 34, 39) -> avoid connecting to none of the four
x q[1];
x q[27];
x q[33];
x q[38];
ccx q[1],q[27],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[33],q[43];
ccx q[1],q[27],q[42];
x q[1];
x q[27];
x q[33];
x q[38];

// blue K4 (2, 28, 36, 39) -> avoid connecting to none of the four
x q[1];
x q[27];
x q[35];
x q[38];
ccx q[1],q[27],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[35],q[43];
ccx q[1],q[27],q[42];
x q[1];
x q[27];
x q[35];
x q[38];

// blue K4 (2, 30, 34, 39) -> avoid connecting to none of the four
x q[1];
x q[29];
x q[33];
x q[38];
ccx q[1],q[29],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[33],q[43];
ccx q[1],q[29],q[42];
x q[1];
x q[29];
x q[33];
x q[38];

// blue K4 (2, 30, 36, 39) -> avoid connecting to none of the four
x q[1];
x q[29];
x q[35];
x q[38];
ccx q[1],q[29],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[35],q[43];
ccx q[1],q[29],q[42];
x q[1];
x q[29];
x q[35];
x q[38];

// blue K4 (2, 30, 36, 41) -> avoid connecting to none of the four
x q[1];
x q[29];
x q[35];
x q[40];
ccx q[1],q[29],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[35],q[43];
ccx q[1],q[29],q[42];
x q[1];
x q[29];
x q[35];
x q[40];

// blue K4 (2, 34, 37, 40) -> avoid connecting to none of the four
x q[1];
x q[33];
x q[36];
x q[39];
ccx q[1],q[33],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[36],q[43];
ccx q[1],q[33],q[42];
x q[1];
x q[33];
x q[36];
x q[39];

// blue K4 (2, 34, 37, 42) -> avoid connecting to none of the four
x q[1];
x q[33];
x q[36];
x q[41];
ccx q[1],q[33],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[36],q[43];
ccx q[1],q[33],q[42];
x q[1];
x q[33];
x q[36];
x q[41];

// blue K4 (2, 34, 39, 40) -> avoid connecting to none of the four
x q[1];
x q[33];
x q[38];
x q[39];
ccx q[1],q[33],q[42];
ccx q[42],q[38],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[38],q[43];
ccx q[1],q[33],q[42];
x q[1];
x q[33];
x q[38];
x q[39];

// blue K4 (2, 34, 39, 42) -> avoid connecting to none of the four
x q[1];
x q[33];
x q[38];
x q[41];
ccx q[1],q[33],q[42];
ccx q[42],q[38],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[38],q[43];
ccx q[1],q[33],q[42];
x q[1];
x q[33];
x q[38];
x q[41];

// blue K4 (2, 36, 39, 40) -> avoid connecting to none of the four
x q[1];
x q[35];
x q[38];
x q[39];
ccx q[1],q[35],q[42];
ccx q[42],q[38],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[38],q[43];
ccx q[1],q[35],q[42];
x q[1];
x q[35];
x q[38];
x q[39];

// blue K4 (2, 36, 39, 42) -> avoid connecting to none of the four
x q[1];
x q[35];
x q[38];
x q[41];
ccx q[1],q[35],q[42];
ccx q[42],q[38],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[38],q[43];
ccx q[1],q[35],q[42];
x q[1];
x q[35];
x q[38];
x q[41];

// blue K4 (2, 36, 40, 41) -> avoid connecting to none of the four
x q[1];
x q[35];
x q[39];
x q[40];
ccx q[1],q[35],q[42];
ccx q[42],q[39],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[39],q[43];
ccx q[1],q[35],q[42];
x q[1];
x q[35];
x q[39];
x q[40];

// blue K4 (2, 36, 41, 42) -> avoid connecting to none of the four
x q[1];
x q[35];
x q[40];
x q[41];
ccx q[1],q[35],q[42];
ccx q[42],q[40],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[40],q[43];
ccx q[1],q[35],q[42];
x q[1];
x q[35];
x q[40];
x q[41];

// blue K4 (2, 37, 40, 41) -> avoid connecting to none of the four
x q[1];
x q[36];
x q[39];
x q[40];
ccx q[1],q[36],q[42];
ccx q[42],q[39],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[39],q[43];
ccx q[1],q[36],q[42];
x q[1];
x q[36];
x q[39];
x q[40];

// blue K4 (2, 37, 41, 42) -> avoid connecting to none of the four
x q[1];
x q[36];
x q[40];
x q[41];
ccx q[1],q[36],q[42];
ccx q[42],q[40],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[40],q[43];
ccx q[1],q[36],q[42];
x q[1];
x q[36];
x q[40];
x q[41];

// blue K4 (3, 6, 7, 11) -> avoid connecting to none of the four
x q[2];
x q[5];
x q[6];
x q[10];
ccx q[2],q[5],q[42];
ccx q[42],q[6],q[43];
ccx q[43],q[10],q[44];
rz(0.15) q[44];
ccx q[43],q[10],q[44];
ccx q[42],q[6],q[43];
ccx q[2],q[5],q[42];
x q[2];
x q[5];
x q[6];
x q[10];

// blue K4 (3, 6, 7, 12) -> avoid connecting to none of the four
x q[2];
x q[5];
x q[6];
x q[11];
ccx q[2],q[5],q[42];
ccx q[42],q[6],q[43];
ccx q[43],q[11],q[44];
rz(0.15) q[44];
ccx q[43],q[11],q[44];
ccx q[42],q[6],q[43];
ccx q[2],q[5],q[42];
x q[2];
x q[5];
x q[6];
x q[11];

// blue K4 (3, 6, 7, 41) -> avoid connecting to none of the four
x q[2];
x q[5];
x q[6];
x q[40];
ccx q[2],q[5],q[42];
ccx q[42],q[6],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[6],q[43];
ccx q[2],q[5],q[42];
x q[2];
x q[5];
x q[6];
x q[40];

// blue K4 (3, 6, 9, 12) -> avoid connecting to none of the four
x q[2];
x q[5];
x q[8];
x q[11];
ccx q[2],q[5],q[42];
ccx q[42],q[8],q[43];
ccx q[43],q[11],q[44];
rz(0.15) q[44];
ccx q[43],q[11],q[44];
ccx q[42],q[8],q[43];
ccx q[2],q[5],q[42];
x q[2];
x q[5];
x q[8];
x q[11];

// blue K4 (3, 6, 9, 14) -> avoid connecting to none of the four
x q[2];
x q[5];
x q[8];
x q[13];
ccx q[2],q[5],q[42];
ccx q[42],q[8],q[43];
ccx q[43],q[13],q[44];
rz(0.15) q[44];
ccx q[43],q[13],q[44];
ccx q[42],q[8],q[43];
ccx q[2],q[5],q[42];
x q[2];
x q[5];
x q[8];
x q[13];

// blue K4 (3, 6, 9, 41) -> avoid connecting to none of the four
x q[2];
x q[5];
x q[8];
x q[40];
ccx q[2],q[5],q[42];
ccx q[42],q[8],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[8],q[43];
ccx q[2],q[5],q[42];
x q[2];
x q[5];
x q[8];
x q[40];

// blue K4 (3, 6, 11, 14) -> avoid connecting to none of the four
x q[2];
x q[5];
x q[10];
x q[13];
ccx q[2],q[5],q[42];
ccx q[42],q[10],q[43];
ccx q[43],q[13],q[44];
rz(0.15) q[44];
ccx q[43],q[13],q[44];
ccx q[42],q[10],q[43];
ccx q[2],q[5],q[42];
x q[2];
x q[5];
x q[10];
x q[13];

// blue K4 (3, 6, 12, 38) -> avoid connecting to none of the four
x q[2];
x q[5];
x q[11];
x q[37];
ccx q[2],q[5],q[42];
ccx q[42],q[11],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[11],q[43];
ccx q[2],q[5],q[42];
x q[2];
x q[5];
x q[11];
x q[37];

// blue K4 (3, 6, 12, 40) -> avoid connecting to none of the four
x q[2];
x q[5];
x q[11];
x q[39];
ccx q[2],q[5],q[42];
ccx q[42],q[11],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[11],q[43];
ccx q[2],q[5],q[42];
x q[2];
x q[5];
x q[11];
x q[39];

// blue K4 (3, 6, 14, 38) -> avoid connecting to none of the four
x q[2];
x q[5];
x q[13];
x q[37];
ccx q[2],q[5],q[42];
ccx q[42],q[13],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[13],q[43];
ccx q[2],q[5],q[42];
x q[2];
x q[5];
x q[13];
x q[37];

// blue K4 (3, 6, 14, 40) -> avoid connecting to none of the four
x q[2];
x q[5];
x q[13];
x q[39];
ccx q[2],q[5],q[42];
ccx q[42],q[13],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[13],q[43];
ccx q[2],q[5],q[42];
x q[2];
x q[5];
x q[13];
x q[39];

// blue K4 (3, 6, 38, 41) -> avoid connecting to none of the four
x q[2];
x q[5];
x q[37];
x q[40];
ccx q[2],q[5],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[37],q[43];
ccx q[2],q[5],q[42];
x q[2];
x q[5];
x q[37];
x q[40];

// blue K4 (3, 6, 40, 41) -> avoid connecting to none of the four
x q[2];
x q[5];
x q[39];
x q[40];
ccx q[2],q[5],q[42];
ccx q[42],q[39],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[39],q[43];
ccx q[2],q[5],q[42];
x q[2];
x q[5];
x q[39];
x q[40];

// blue K4 (3, 7, 8, 11) -> avoid connecting to none of the four
x q[2];
x q[6];
x q[7];
x q[10];
ccx q[2],q[6],q[42];
ccx q[42],q[7],q[43];
ccx q[43],q[10],q[44];
rz(0.15) q[44];
ccx q[43],q[10],q[44];
ccx q[42],q[7],q[43];
ccx q[2],q[6],q[42];
x q[2];
x q[6];
x q[7];
x q[10];

// blue K4 (3, 7, 8, 12) -> avoid connecting to none of the four
x q[2];
x q[6];
x q[7];
x q[11];
ccx q[2],q[6],q[42];
ccx q[42],q[7],q[43];
ccx q[43],q[11],q[44];
rz(0.15) q[44];
ccx q[43],q[11],q[44];
ccx q[42],q[7],q[43];
ccx q[2],q[6],q[42];
x q[2];
x q[6];
x q[7];
x q[11];

// blue K4 (3, 7, 8, 42) -> avoid connecting to none of the four
x q[2];
x q[6];
x q[7];
x q[41];
ccx q[2],q[6],q[42];
ccx q[42],q[7],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[7],q[43];
ccx q[2],q[6],q[42];
x q[2];
x q[6];
x q[7];
x q[41];

// blue K4 (3, 7, 11, 22) -> avoid connecting to none of the four
x q[2];
x q[6];
x q[10];
x q[21];
ccx q[2],q[6],q[42];
ccx q[42],q[10],q[43];
ccx q[43],q[21],q[44];
rz(0.15) q[44];
ccx q[43],q[21],q[44];
ccx q[42],q[10],q[43];
ccx q[2],q[6],q[42];
x q[2];
x q[6];
x q[10];
x q[21];

// blue K4 (3, 7, 11, 35) -> avoid connecting to none of the four
x q[2];
x q[6];
x q[10];
x q[34];
ccx q[2],q[6],q[42];
ccx q[42],q[10],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[10],q[43];
ccx q[2],q[6],q[42];
x q[2];
x q[6];
x q[10];
x q[34];

// blue K4 (3, 7, 12, 18) -> avoid connecting to none of the four
x q[2];
x q[6];
x q[11];
x q[17];
ccx q[2],q[6],q[42];
ccx q[42],q[11],q[43];
ccx q[43],q[17],q[44];
rz(0.15) q[44];
ccx q[43],q[17],q[44];
ccx q[42],q[11],q[43];
ccx q[2],q[6],q[42];
x q[2];
x q[6];
x q[11];
x q[17];

// blue K4 (3, 7, 12, 31) -> avoid connecting to none of the four
x q[2];
x q[6];
x q[11];
x q[30];
ccx q[2],q[6],q[42];
ccx q[42],q[11],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[11],q[43];
ccx q[2],q[6],q[42];
x q[2];
x q[6];
x q[11];
x q[30];

// blue K4 (3, 7, 18, 22) -> avoid connecting to none of the four
x q[2];
x q[6];
x q[17];
x q[21];
ccx q[2],q[6],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[21],q[44];
rz(0.15) q[44];
ccx q[43],q[21],q[44];
ccx q[42],q[17],q[43];
ccx q[2],q[6],q[42];
x q[2];
x q[6];
x q[17];
x q[21];

// blue K4 (3, 7, 18, 35) -> avoid connecting to none of the four
x q[2];
x q[6];
x q[17];
x q[34];
ccx q[2],q[6],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[17],q[43];
ccx q[2],q[6],q[42];
x q[2];
x q[6];
x q[17];
x q[34];

// blue K4 (3, 7, 18, 42) -> avoid connecting to none of the four
x q[2];
x q[6];
x q[17];
x q[41];
ccx q[2],q[6],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[17],q[43];
ccx q[2],q[6],q[42];
x q[2];
x q[6];
x q[17];
x q[41];

// blue K4 (3, 7, 22, 31) -> avoid connecting to none of the four
x q[2];
x q[6];
x q[21];
x q[30];
ccx q[2],q[6],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[21],q[43];
ccx q[2],q[6],q[42];
x q[2];
x q[6];
x q[21];
x q[30];

// blue K4 (3, 7, 22, 41) -> avoid connecting to none of the four
x q[2];
x q[6];
x q[21];
x q[40];
ccx q[2],q[6],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[21],q[43];
ccx q[2],q[6],q[42];
x q[2];
x q[6];
x q[21];
x q[40];

// blue K4 (3, 7, 31, 35) -> avoid connecting to none of the four
x q[2];
x q[6];
x q[30];
x q[34];
ccx q[2],q[6],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[30],q[43];
ccx q[2],q[6],q[42];
x q[2];
x q[6];
x q[30];
x q[34];

// blue K4 (3, 7, 31, 42) -> avoid connecting to none of the four
x q[2];
x q[6];
x q[30];
x q[41];
ccx q[2],q[6],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[30],q[43];
ccx q[2],q[6],q[42];
x q[2];
x q[6];
x q[30];
x q[41];

// blue K4 (3, 7, 35, 41) -> avoid connecting to none of the four
x q[2];
x q[6];
x q[34];
x q[40];
ccx q[2],q[6],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[34],q[43];
ccx q[2],q[6],q[42];
x q[2];
x q[6];
x q[34];
x q[40];

// blue K4 (3, 7, 41, 42) -> avoid connecting to none of the four
x q[2];
x q[6];
x q[40];
x q[41];
ccx q[2],q[6],q[42];
ccx q[42],q[40],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[40],q[43];
ccx q[2],q[6],q[42];
x q[2];
x q[6];
x q[40];
x q[41];

// blue K4 (3, 8, 11, 14) -> avoid connecting to none of the four
x q[2];
x q[7];
x q[10];
x q[13];
ccx q[2],q[7],q[42];
ccx q[42],q[10],q[43];
ccx q[43],q[13],q[44];
rz(0.15) q[44];
ccx q[43],q[13],q[44];
ccx q[42],q[10],q[43];
ccx q[2],q[7],q[42];
x q[2];
x q[7];
x q[10];
x q[13];

// blue K4 (3, 8, 12, 27) -> avoid connecting to none of the four
x q[2];
x q[7];
x q[11];
x q[26];
ccx q[2],q[7],q[42];
ccx q[42],q[11],q[43];
ccx q[43],q[26],q[44];
rz(0.15) q[44];
ccx q[43],q[26],q[44];
ccx q[42],q[11],q[43];
ccx q[2],q[7],q[42];
x q[2];
x q[7];
x q[11];
x q[26];

// blue K4 (3, 8, 12, 40) -> avoid connecting to none of the four
x q[2];
x q[7];
x q[11];
x q[39];
ccx q[2],q[7],q[42];
ccx q[42],q[11],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[11],q[43];
ccx q[2],q[7],q[42];
x q[2];
x q[7];
x q[11];
x q[39];

// blue K4 (3, 8, 14, 40) -> avoid connecting to none of the four
x q[2];
x q[7];
x q[13];
x q[39];
ccx q[2],q[7],q[42];
ccx q[42],q[13],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[13],q[43];
ccx q[2],q[7],q[42];
x q[2];
x q[7];
x q[13];
x q[39];

// blue K4 (3, 8, 14, 42) -> avoid connecting to none of the four
x q[2];
x q[7];
x q[13];
x q[41];
ccx q[2],q[7],q[42];
ccx q[42],q[13],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[13],q[43];
ccx q[2],q[7],q[42];
x q[2];
x q[7];
x q[13];
x q[41];

// blue K4 (3, 8, 27, 42) -> avoid connecting to none of the four
x q[2];
x q[7];
x q[26];
x q[41];
ccx q[2],q[7],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[26],q[43];
ccx q[2],q[7],q[42];
x q[2];
x q[7];
x q[26];
x q[41];

// blue K4 (3, 9, 12, 18) -> avoid connecting to none of the four
x q[2];
x q[8];
x q[11];
x q[17];
ccx q[2],q[8],q[42];
ccx q[42],q[11],q[43];
ccx q[43],q[17],q[44];
rz(0.15) q[44];
ccx q[43],q[17],q[44];
ccx q[42],q[11],q[43];
ccx q[2],q[8],q[42];
x q[2];
x q[8];
x q[11];
x q[17];

// blue K4 (3, 9, 12, 20) -> avoid connecting to none of the four
x q[2];
x q[8];
x q[11];
x q[19];
ccx q[2],q[8],q[42];
ccx q[42],q[11],q[43];
ccx q[43],q[19],q[44];
rz(0.15) q[44];
ccx q[43],q[19],q[44];
ccx q[42],q[11],q[43];
ccx q[2],q[8],q[42];
x q[2];
x q[8];
x q[11];
x q[19];

// blue K4 (3, 9, 14, 18) -> avoid connecting to none of the four
x q[2];
x q[8];
x q[13];
x q[17];
ccx q[2],q[8],q[42];
ccx q[42],q[13],q[43];
ccx q[43],q[17],q[44];
rz(0.15) q[44];
ccx q[43],q[17],q[44];
ccx q[42],q[13],q[43];
ccx q[2],q[8],q[42];
x q[2];
x q[8];
x q[13];
x q[17];

// blue K4 (3, 9, 14, 20) -> avoid connecting to none of the four
x q[2];
x q[8];
x q[13];
x q[19];
ccx q[2],q[8],q[42];
ccx q[42],q[13],q[43];
ccx q[43],q[19],q[44];
rz(0.15) q[44];
ccx q[43],q[19],q[44];
ccx q[42],q[13],q[43];
ccx q[2],q[8],q[42];
x q[2];
x q[8];
x q[13];
x q[19];

// blue K4 (3, 9, 18, 35) -> avoid connecting to none of the four
x q[2];
x q[8];
x q[17];
x q[34];
ccx q[2],q[8],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[17],q[43];
ccx q[2],q[8],q[42];
x q[2];
x q[8];
x q[17];
x q[34];

// blue K4 (3, 9, 18, 37) -> avoid connecting to none of the four
x q[2];
x q[8];
x q[17];
x q[36];
ccx q[2],q[8],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[17],q[43];
ccx q[2],q[8],q[42];
x q[2];
x q[8];
x q[17];
x q[36];

// blue K4 (3, 9, 20, 35) -> avoid connecting to none of the four
x q[2];
x q[8];
x q[19];
x q[34];
ccx q[2],q[8],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[19],q[43];
ccx q[2],q[8],q[42];
x q[2];
x q[8];
x q[19];
x q[34];

// blue K4 (3, 9, 20, 37) -> avoid connecting to none of the four
x q[2];
x q[8];
x q[19];
x q[36];
ccx q[2],q[8],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[19],q[43];
ccx q[2],q[8],q[42];
x q[2];
x q[8];
x q[19];
x q[36];

// blue K4 (3, 9, 35, 41) -> avoid connecting to none of the four
x q[2];
x q[8];
x q[34];
x q[40];
ccx q[2],q[8],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[34],q[43];
ccx q[2],q[8],q[42];
x q[2];
x q[8];
x q[34];
x q[40];

// blue K4 (3, 9, 37, 41) -> avoid connecting to none of the four
x q[2];
x q[8];
x q[36];
x q[40];
ccx q[2],q[8],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[36],q[43];
ccx q[2],q[8],q[42];
x q[2];
x q[8];
x q[36];
x q[40];

// blue K4 (3, 11, 14, 20) -> avoid connecting to none of the four
x q[2];
x q[10];
x q[13];
x q[19];
ccx q[2],q[10],q[42];
ccx q[42],q[13],q[43];
ccx q[43],q[19],q[44];
rz(0.15) q[44];
ccx q[43],q[19],q[44];
ccx q[42],q[13],q[43];
ccx q[2],q[10],q[42];
x q[2];
x q[10];
x q[13];
x q[19];

// blue K4 (3, 11, 14, 22) -> avoid connecting to none of the four
x q[2];
x q[10];
x q[13];
x q[21];
ccx q[2],q[10],q[42];
ccx q[42],q[13],q[43];
ccx q[43],q[21],q[44];
rz(0.15) q[44];
ccx q[43],q[21],q[44];
ccx q[42],q[13],q[43];
ccx q[2],q[10],q[42];
x q[2];
x q[10];
x q[13];
x q[21];

// blue K4 (3, 11, 20, 35) -> avoid connecting to none of the four
x q[2];
x q[10];
x q[19];
x q[34];
ccx q[2],q[10],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[19],q[43];
ccx q[2],q[10],q[42];
x q[2];
x q[10];
x q[19];
x q[34];

// blue K4 (3, 11, 20, 37) -> avoid connecting to none of the four
x q[2];
x q[10];
x q[19];
x q[36];
ccx q[2],q[10],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[19],q[43];
ccx q[2],q[10],q[42];
x q[2];
x q[10];
x q[19];
x q[36];

// blue K4 (3, 11, 22, 37) -> avoid connecting to none of the four
x q[2];
x q[10];
x q[21];
x q[36];
ccx q[2],q[10],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[21],q[43];
ccx q[2],q[10],q[42];
x q[2];
x q[10];
x q[21];
x q[36];

// blue K4 (3, 12, 18, 27) -> avoid connecting to none of the four
x q[2];
x q[11];
x q[17];
x q[26];
ccx q[2],q[11],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[26],q[44];
rz(0.15) q[44];
ccx q[43],q[26],q[44];
ccx q[42],q[17],q[43];
ccx q[2],q[11],q[42];
x q[2];
x q[11];
x q[17];
x q[26];

// blue K4 (3, 12, 18, 29) -> avoid connecting to none of the four
x q[2];
x q[11];
x q[17];
x q[28];
ccx q[2],q[11],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[17],q[43];
ccx q[2],q[11],q[42];
x q[2];
x q[11];
x q[17];
x q[28];

// blue K4 (3, 12, 20, 29) -> avoid connecting to none of the four
x q[2];
x q[11];
x q[19];
x q[28];
ccx q[2],q[11],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[19],q[43];
ccx q[2],q[11],q[42];
x q[2];
x q[11];
x q[19];
x q[28];

// blue K4 (3, 12, 20, 31) -> avoid connecting to none of the four
x q[2];
x q[11];
x q[19];
x q[30];
ccx q[2],q[11],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[19],q[43];
ccx q[2],q[11],q[42];
x q[2];
x q[11];
x q[19];
x q[30];

// blue K4 (3, 12, 27, 31) -> avoid connecting to none of the four
x q[2];
x q[11];
x q[26];
x q[30];
ccx q[2],q[11],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[26],q[43];
ccx q[2],q[11],q[42];
x q[2];
x q[11];
x q[26];
x q[30];

// blue K4 (3, 12, 27, 38) -> avoid connecting to none of the four
x q[2];
x q[11];
x q[26];
x q[37];
ccx q[2],q[11],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[26],q[43];
ccx q[2],q[11],q[42];
x q[2];
x q[11];
x q[26];
x q[37];

// blue K4 (3, 12, 29, 38) -> avoid connecting to none of the four
x q[2];
x q[11];
x q[28];
x q[37];
ccx q[2],q[11],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[28],q[43];
ccx q[2],q[11],q[42];
x q[2];
x q[11];
x q[28];
x q[37];

// blue K4 (3, 12, 29, 40) -> avoid connecting to none of the four
x q[2];
x q[11];
x q[28];
x q[39];
ccx q[2],q[11],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[28],q[43];
ccx q[2],q[11],q[42];
x q[2];
x q[11];
x q[28];
x q[39];

// blue K4 (3, 12, 31, 40) -> avoid connecting to none of the four
x q[2];
x q[11];
x q[30];
x q[39];
ccx q[2],q[11],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[30],q[43];
ccx q[2],q[11],q[42];
x q[2];
x q[11];
x q[30];
x q[39];

// blue K4 (3, 14, 18, 22) -> avoid connecting to none of the four
x q[2];
x q[13];
x q[17];
x q[21];
ccx q[2],q[13],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[21],q[44];
rz(0.15) q[44];
ccx q[43],q[21],q[44];
ccx q[42],q[17],q[43];
ccx q[2],q[13],q[42];
x q[2];
x q[13];
x q[17];
x q[21];

// blue K4 (3, 14, 18, 29) -> avoid connecting to none of the four
x q[2];
x q[13];
x q[17];
x q[28];
ccx q[2],q[13],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[17],q[43];
ccx q[2],q[13],q[42];
x q[2];
x q[13];
x q[17];
x q[28];

// blue K4 (3, 14, 18, 42) -> avoid connecting to none of the four
x q[2];
x q[13];
x q[17];
x q[41];
ccx q[2],q[13],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[17],q[43];
ccx q[2],q[13],q[42];
x q[2];
x q[13];
x q[17];
x q[41];

// blue K4 (3, 14, 20, 29) -> avoid connecting to none of the four
x q[2];
x q[13];
x q[19];
x q[28];
ccx q[2],q[13],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[19],q[43];
ccx q[2],q[13],q[42];
x q[2];
x q[13];
x q[19];
x q[28];

// blue K4 (3, 14, 20, 31) -> avoid connecting to none of the four
x q[2];
x q[13];
x q[19];
x q[30];
ccx q[2],q[13],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[19],q[43];
ccx q[2],q[13],q[42];
x q[2];
x q[13];
x q[19];
x q[30];

// blue K4 (3, 14, 22, 31) -> avoid connecting to none of the four
x q[2];
x q[13];
x q[21];
x q[30];
ccx q[2],q[13],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[21],q[43];
ccx q[2],q[13],q[42];
x q[2];
x q[13];
x q[21];
x q[30];

// blue K4 (3, 14, 29, 38) -> avoid connecting to none of the four
x q[2];
x q[13];
x q[28];
x q[37];
ccx q[2],q[13],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[28],q[43];
ccx q[2],q[13],q[42];
x q[2];
x q[13];
x q[28];
x q[37];

// blue K4 (3, 14, 29, 40) -> avoid connecting to none of the four
x q[2];
x q[13];
x q[28];
x q[39];
ccx q[2],q[13],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[28],q[43];
ccx q[2],q[13],q[42];
x q[2];
x q[13];
x q[28];
x q[39];

// blue K4 (3, 14, 31, 40) -> avoid connecting to none of the four
x q[2];
x q[13];
x q[30];
x q[39];
ccx q[2],q[13],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[30],q[43];
ccx q[2],q[13],q[42];
x q[2];
x q[13];
x q[30];
x q[39];

// blue K4 (3, 14, 31, 42) -> avoid connecting to none of the four
x q[2];
x q[13];
x q[30];
x q[41];
ccx q[2],q[13],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[30],q[43];
ccx q[2],q[13],q[42];
x q[2];
x q[13];
x q[30];
x q[41];

// blue K4 (3, 14, 38, 42) -> avoid connecting to none of the four
x q[2];
x q[13];
x q[37];
x q[41];
ccx q[2],q[13],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[37],q[43];
ccx q[2],q[13],q[42];
x q[2];
x q[13];
x q[37];
x q[41];

// blue K4 (3, 18, 22, 27) -> avoid connecting to none of the four
x q[2];
x q[17];
x q[21];
x q[26];
ccx q[2],q[17],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[26],q[44];
rz(0.15) q[44];
ccx q[43],q[26],q[44];
ccx q[42],q[21],q[43];
ccx q[2],q[17],q[42];
x q[2];
x q[17];
x q[21];
x q[26];

// blue K4 (3, 18, 22, 37) -> avoid connecting to none of the four
x q[2];
x q[17];
x q[21];
x q[36];
ccx q[2],q[17],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[21],q[43];
ccx q[2],q[17],q[42];
x q[2];
x q[17];
x q[21];
x q[36];

// blue K4 (3, 18, 27, 35) -> avoid connecting to none of the four
x q[2];
x q[17];
x q[26];
x q[34];
ccx q[2],q[17],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[26],q[43];
ccx q[2],q[17],q[42];
x q[2];
x q[17];
x q[26];
x q[34];

// blue K4 (3, 18, 27, 42) -> avoid connecting to none of the four
x q[2];
x q[17];
x q[26];
x q[41];
ccx q[2],q[17],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[26],q[43];
ccx q[2],q[17],q[42];
x q[2];
x q[17];
x q[26];
x q[41];

// blue K4 (3, 18, 29, 35) -> avoid connecting to none of the four
x q[2];
x q[17];
x q[28];
x q[34];
ccx q[2],q[17],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[28],q[43];
ccx q[2],q[17],q[42];
x q[2];
x q[17];
x q[28];
x q[34];

// blue K4 (3, 18, 29, 37) -> avoid connecting to none of the four
x q[2];
x q[17];
x q[28];
x q[36];
ccx q[2],q[17],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[28],q[43];
ccx q[2],q[17],q[42];
x q[2];
x q[17];
x q[28];
x q[36];

// blue K4 (3, 18, 37, 42) -> avoid connecting to none of the four
x q[2];
x q[17];
x q[36];
x q[41];
ccx q[2],q[17],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[36],q[43];
ccx q[2],q[17],q[42];
x q[2];
x q[17];
x q[36];
x q[41];

// blue K4 (3, 20, 29, 35) -> avoid connecting to none of the four
x q[2];
x q[19];
x q[28];
x q[34];
ccx q[2],q[19],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[28],q[43];
ccx q[2],q[19],q[42];
x q[2];
x q[19];
x q[28];
x q[34];

// blue K4 (3, 20, 29, 37) -> avoid connecting to none of the four
x q[2];
x q[19];
x q[28];
x q[36];
ccx q[2],q[19],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[28],q[43];
ccx q[2],q[19],q[42];
x q[2];
x q[19];
x q[28];
x q[36];

// blue K4 (3, 20, 31, 35) -> avoid connecting to none of the four
x q[2];
x q[19];
x q[30];
x q[34];
ccx q[2],q[19],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[30],q[43];
ccx q[2],q[19],q[42];
x q[2];
x q[19];
x q[30];
x q[34];

// blue K4 (3, 20, 31, 37) -> avoid connecting to none of the four
x q[2];
x q[19];
x q[30];
x q[36];
ccx q[2],q[19],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[30],q[43];
ccx q[2],q[19],q[42];
x q[2];
x q[19];
x q[30];
x q[36];

// blue K4 (3, 22, 27, 31) -> avoid connecting to none of the four
x q[2];
x q[21];
x q[26];
x q[30];
ccx q[2],q[21],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[26],q[43];
ccx q[2],q[21],q[42];
x q[2];
x q[21];
x q[26];
x q[30];

// blue K4 (3, 22, 31, 37) -> avoid connecting to none of the four
x q[2];
x q[21];
x q[30];
x q[36];
ccx q[2],q[21],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[30],q[43];
ccx q[2],q[21],q[42];
x q[2];
x q[21];
x q[30];
x q[36];

// blue K4 (3, 22, 37, 41) -> avoid connecting to none of the four
x q[2];
x q[21];
x q[36];
x q[40];
ccx q[2],q[21],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[36],q[43];
ccx q[2],q[21],q[42];
x q[2];
x q[21];
x q[36];
x q[40];

// blue K4 (3, 27, 31, 35) -> avoid connecting to none of the four
x q[2];
x q[26];
x q[30];
x q[34];
ccx q[2],q[26],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[30],q[43];
ccx q[2],q[26],q[42];
x q[2];
x q[26];
x q[30];
x q[34];

// blue K4 (3, 27, 31, 42) -> avoid connecting to none of the four
x q[2];
x q[26];
x q[30];
x q[41];
ccx q[2],q[26],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[30],q[43];
ccx q[2],q[26],q[42];
x q[2];
x q[26];
x q[30];
x q[41];

// blue K4 (3, 27, 35, 38) -> avoid connecting to none of the four
x q[2];
x q[26];
x q[34];
x q[37];
ccx q[2],q[26],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[34],q[43];
ccx q[2],q[26],q[42];
x q[2];
x q[26];
x q[34];
x q[37];

// blue K4 (3, 27, 38, 42) -> avoid connecting to none of the four
x q[2];
x q[26];
x q[37];
x q[41];
ccx q[2],q[26],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[37],q[43];
ccx q[2],q[26],q[42];
x q[2];
x q[26];
x q[37];
x q[41];

// blue K4 (3, 29, 35, 38) -> avoid connecting to none of the four
x q[2];
x q[28];
x q[34];
x q[37];
ccx q[2],q[28],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[34],q[43];
ccx q[2],q[28],q[42];
x q[2];
x q[28];
x q[34];
x q[37];

// blue K4 (3, 29, 35, 40) -> avoid connecting to none of the four
x q[2];
x q[28];
x q[34];
x q[39];
ccx q[2],q[28],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[34],q[43];
ccx q[2],q[28],q[42];
x q[2];
x q[28];
x q[34];
x q[39];

// blue K4 (3, 29, 37, 40) -> avoid connecting to none of the four
x q[2];
x q[28];
x q[36];
x q[39];
ccx q[2],q[28],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[36],q[43];
ccx q[2],q[28],q[42];
x q[2];
x q[28];
x q[36];
x q[39];

// blue K4 (3, 31, 35, 40) -> avoid connecting to none of the four
x q[2];
x q[30];
x q[34];
x q[39];
ccx q[2],q[30],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[34],q[43];
ccx q[2],q[30],q[42];
x q[2];
x q[30];
x q[34];
x q[39];

// blue K4 (3, 31, 37, 40) -> avoid connecting to none of the four
x q[2];
x q[30];
x q[36];
x q[39];
ccx q[2],q[30],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[36],q[43];
ccx q[2],q[30],q[42];
x q[2];
x q[30];
x q[36];
x q[39];

// blue K4 (3, 31, 37, 42) -> avoid connecting to none of the four
x q[2];
x q[30];
x q[36];
x q[41];
ccx q[2],q[30],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[36],q[43];
ccx q[2],q[30],q[42];
x q[2];
x q[30];
x q[36];
x q[41];

// blue K4 (3, 35, 38, 41) -> avoid connecting to none of the four
x q[2];
x q[34];
x q[37];
x q[40];
ccx q[2],q[34],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[37],q[43];
ccx q[2],q[34],q[42];
x q[2];
x q[34];
x q[37];
x q[40];

// blue K4 (3, 35, 40, 41) -> avoid connecting to none of the four
x q[2];
x q[34];
x q[39];
x q[40];
ccx q[2],q[34],q[42];
ccx q[42],q[39],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[39],q[43];
ccx q[2],q[34],q[42];
x q[2];
x q[34];
x q[39];
x q[40];

// blue K4 (3, 37, 40, 41) -> avoid connecting to none of the four
x q[2];
x q[36];
x q[39];
x q[40];
ccx q[2],q[36],q[42];
ccx q[42],q[39],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[39],q[43];
ccx q[2],q[36],q[42];
x q[2];
x q[36];
x q[39];
x q[40];

// blue K4 (3, 37, 41, 42) -> avoid connecting to none of the four
x q[2];
x q[36];
x q[40];
x q[41];
ccx q[2],q[36],q[42];
ccx q[42],q[40],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[40],q[43];
ccx q[2],q[36],q[42];
x q[2];
x q[36];
x q[40];
x q[41];

// blue K4 (3, 38, 41, 42) -> avoid connecting to none of the four
x q[2];
x q[37];
x q[40];
x q[41];
ccx q[2],q[37],q[42];
ccx q[42],q[40],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[40],q[43];
ccx q[2],q[37],q[42];
x q[2];
x q[37];
x q[40];
x q[41];

// blue K4 (4, 5, 8, 13) -> avoid connecting to none of the four
x q[3];
x q[4];
x q[7];
x q[12];
ccx q[3],q[4],q[42];
ccx q[42],q[7],q[43];
ccx q[43],q[12],q[44];
rz(0.15) q[44];
ccx q[43],q[12],q[44];
ccx q[42],q[7],q[43];
ccx q[3],q[4],q[42];
x q[3];
x q[4];
x q[7];
x q[12];

// blue K4 (4, 5, 8, 42) -> avoid connecting to none of the four
x q[3];
x q[4];
x q[7];
x q[41];
ccx q[3],q[4],q[42];
ccx q[42],q[7],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[7],q[43];
ccx q[3],q[4],q[42];
x q[3];
x q[4];
x q[7];
x q[41];

// blue K4 (4, 5, 9, 13) -> avoid connecting to none of the four
x q[3];
x q[4];
x q[8];
x q[12];
ccx q[3],q[4],q[42];
ccx q[42],q[8],q[43];
ccx q[43],q[12],q[44];
rz(0.15) q[44];
ccx q[43],q[12],q[44];
ccx q[42],q[8],q[43];
ccx q[3],q[4],q[42];
x q[3];
x q[4];
x q[8];
x q[12];

// blue K4 (4, 5, 10, 13) -> avoid connecting to none of the four
x q[3];
x q[4];
x q[9];
x q[12];
ccx q[3],q[4],q[42];
ccx q[42],q[9],q[43];
ccx q[43],q[12],q[44];
rz(0.15) q[44];
ccx q[43],q[12],q[44];
ccx q[42],q[9],q[43];
ccx q[3],q[4],q[42];
x q[3];
x q[4];
x q[9];
x q[12];

// blue K4 (4, 5, 10, 42) -> avoid connecting to none of the four
x q[3];
x q[4];
x q[9];
x q[41];
ccx q[3],q[4],q[42];
ccx q[42],q[9],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[9],q[43];
ccx q[3],q[4],q[42];
x q[3];
x q[4];
x q[9];
x q[41];

// blue K4 (4, 5, 13, 39) -> avoid connecting to none of the four
x q[3];
x q[4];
x q[12];
x q[38];
ccx q[3],q[4],q[42];
ccx q[42],q[12],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[12],q[43];
ccx q[3],q[4],q[42];
x q[3];
x q[4];
x q[12];
x q[38];

// blue K4 (4, 5, 39, 42) -> avoid connecting to none of the four
x q[3];
x q[4];
x q[38];
x q[41];
ccx q[3],q[4],q[42];
ccx q[42],q[38],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[38],q[43];
ccx q[3],q[4],q[42];
x q[3];
x q[4];
x q[38];
x q[41];

// blue K4 (4, 7, 8, 12) -> avoid connecting to none of the four
x q[3];
x q[6];
x q[7];
x q[11];
ccx q[3],q[6],q[42];
ccx q[42],q[7],q[43];
ccx q[43],q[11],q[44];
rz(0.15) q[44];
ccx q[43],q[11],q[44];
ccx q[42],q[7],q[43];
ccx q[3],q[6],q[42];
x q[3];
x q[6];
x q[7];
x q[11];

// blue K4 (4, 7, 8, 13) -> avoid connecting to none of the four
x q[3];
x q[6];
x q[7];
x q[12];
ccx q[3],q[6],q[42];
ccx q[42],q[7],q[43];
ccx q[43],q[12],q[44];
rz(0.15) q[44];
ccx q[43],q[12],q[44];
ccx q[42],q[7],q[43];
ccx q[3],q[6],q[42];
x q[3];
x q[6];
x q[7];
x q[12];

// blue K4 (4, 7, 8, 42) -> avoid connecting to none of the four
x q[3];
x q[6];
x q[7];
x q[41];
ccx q[3],q[6],q[42];
ccx q[42],q[7],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[7],q[43];
ccx q[3],q[6],q[42];
x q[3];
x q[6];
x q[7];
x q[41];

// blue K4 (4, 7, 10, 13) -> avoid connecting to none of the four
x q[3];
x q[6];
x q[9];
x q[12];
ccx q[3],q[6],q[42];
ccx q[42],q[9],q[43];
ccx q[43],q[12],q[44];
rz(0.15) q[44];
ccx q[43],q[12],q[44];
ccx q[42],q[9],q[43];
ccx q[3],q[6],q[42];
x q[3];
x q[6];
x q[9];
x q[12];

// blue K4 (4, 7, 10, 15) -> avoid connecting to none of the four
x q[3];
x q[6];
x q[9];
x q[14];
ccx q[3],q[6],q[42];
ccx q[42],q[9],q[43];
ccx q[43],q[14],q[44];
rz(0.15) q[44];
ccx q[43],q[14],q[44];
ccx q[42],q[9],q[43];
ccx q[3],q[6],q[42];
x q[3];
x q[6];
x q[9];
x q[14];

// blue K4 (4, 7, 10, 42) -> avoid connecting to none of the four
x q[3];
x q[6];
x q[9];
x q[41];
ccx q[3],q[6],q[42];
ccx q[42],q[9],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[9],q[43];
ccx q[3],q[6],q[42];
x q[3];
x q[6];
x q[9];
x q[41];

// blue K4 (4, 7, 12, 15) -> avoid connecting to none of the four
x q[3];
x q[6];
x q[11];
x q[14];
ccx q[3],q[6],q[42];
ccx q[42],q[11],q[43];
ccx q[43],q[14],q[44];
rz(0.15) q[44];
ccx q[43],q[14],q[44];
ccx q[42],q[11],q[43];
ccx q[3],q[6],q[42];
x q[3];
x q[6];
x q[11];
x q[14];

// blue K4 (4, 7, 13, 39) -> avoid connecting to none of the four
x q[3];
x q[6];
x q[12];
x q[38];
ccx q[3],q[6],q[42];
ccx q[42],q[12],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[12],q[43];
ccx q[3],q[6],q[42];
x q[3];
x q[6];
x q[12];
x q[38];

// blue K4 (4, 7, 13, 41) -> avoid connecting to none of the four
x q[3];
x q[6];
x q[12];
x q[40];
ccx q[3],q[6],q[42];
ccx q[42],q[12],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[12],q[43];
ccx q[3],q[6],q[42];
x q[3];
x q[6];
x q[12];
x q[40];

// blue K4 (4, 7, 15, 39) -> avoid connecting to none of the four
x q[3];
x q[6];
x q[14];
x q[38];
ccx q[3],q[6],q[42];
ccx q[42],q[14],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[14],q[43];
ccx q[3],q[6],q[42];
x q[3];
x q[6];
x q[14];
x q[38];

// blue K4 (4, 7, 15, 41) -> avoid connecting to none of the four
x q[3];
x q[6];
x q[14];
x q[40];
ccx q[3],q[6],q[42];
ccx q[42],q[14],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[14],q[43];
ccx q[3],q[6],q[42];
x q[3];
x q[6];
x q[14];
x q[40];

// blue K4 (4, 7, 39, 42) -> avoid connecting to none of the four
x q[3];
x q[6];
x q[38];
x q[41];
ccx q[3],q[6],q[42];
ccx q[42],q[38],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[38],q[43];
ccx q[3],q[6],q[42];
x q[3];
x q[6];
x q[38];
x q[41];

// blue K4 (4, 7, 41, 42) -> avoid connecting to none of the four
x q[3];
x q[6];
x q[40];
x q[41];
ccx q[3],q[6],q[42];
ccx q[42],q[40],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[40],q[43];
ccx q[3],q[6],q[42];
x q[3];
x q[6];
x q[40];
x q[41];

// blue K4 (4, 8, 12, 23) -> avoid connecting to none of the four
x q[3];
x q[7];
x q[11];
x q[22];
ccx q[3],q[7],q[42];
ccx q[42],q[11],q[43];
ccx q[43],q[22],q[44];
rz(0.15) q[44];
ccx q[43],q[22],q[44];
ccx q[42],q[11],q[43];
ccx q[3],q[7],q[42];
x q[3];
x q[7];
x q[11];
x q[22];

// blue K4 (4, 8, 12, 36) -> avoid connecting to none of the four
x q[3];
x q[7];
x q[11];
x q[35];
ccx q[3],q[7],q[42];
ccx q[42],q[11],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[11],q[43];
ccx q[3],q[7],q[42];
x q[3];
x q[7];
x q[11];
x q[35];

// blue K4 (4, 8, 13, 19) -> avoid connecting to none of the four
x q[3];
x q[7];
x q[12];
x q[18];
ccx q[3],q[7],q[42];
ccx q[42],q[12],q[43];
ccx q[43],q[18],q[44];
rz(0.15) q[44];
ccx q[43],q[18],q[44];
ccx q[42],q[12],q[43];
ccx q[3],q[7],q[42];
x q[3];
x q[7];
x q[12];
x q[18];

// blue K4 (4, 8, 13, 32) -> avoid connecting to none of the four
x q[3];
x q[7];
x q[12];
x q[31];
ccx q[3],q[7],q[42];
ccx q[42],q[12],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[12],q[43];
ccx q[3],q[7],q[42];
x q[3];
x q[7];
x q[12];
x q[31];

// blue K4 (4, 8, 19, 23) -> avoid connecting to none of the four
x q[3];
x q[7];
x q[18];
x q[22];
ccx q[3],q[7],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[22],q[44];
rz(0.15) q[44];
ccx q[43],q[22],q[44];
ccx q[42],q[18],q[43];
ccx q[3],q[7],q[42];
x q[3];
x q[7];
x q[18];
x q[22];

// blue K4 (4, 8, 19, 36) -> avoid connecting to none of the four
x q[3];
x q[7];
x q[18];
x q[35];
ccx q[3],q[7],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[18],q[43];
ccx q[3],q[7],q[42];
x q[3];
x q[7];
x q[18];
x q[35];

// blue K4 (4, 8, 23, 32) -> avoid connecting to none of the four
x q[3];
x q[7];
x q[22];
x q[31];
ccx q[3],q[7],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[22],q[43];
ccx q[3],q[7],q[42];
x q[3];
x q[7];
x q[22];
x q[31];

// blue K4 (4, 8, 23, 42) -> avoid connecting to none of the four
x q[3];
x q[7];
x q[22];
x q[41];
ccx q[3],q[7],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[22],q[43];
ccx q[3],q[7],q[42];
x q[3];
x q[7];
x q[22];
x q[41];

// blue K4 (4, 8, 32, 36) -> avoid connecting to none of the four
x q[3];
x q[7];
x q[31];
x q[35];
ccx q[3],q[7],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[31],q[43];
ccx q[3],q[7],q[42];
x q[3];
x q[7];
x q[31];
x q[35];

// blue K4 (4, 8, 36, 42) -> avoid connecting to none of the four
x q[3];
x q[7];
x q[35];
x q[41];
ccx q[3],q[7],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[35],q[43];
ccx q[3],q[7],q[42];
x q[3];
x q[7];
x q[35];
x q[41];

// blue K4 (4, 9, 12, 15) -> avoid connecting to none of the four
x q[3];
x q[8];
x q[11];
x q[14];
ccx q[3],q[8],q[42];
ccx q[42],q[11],q[43];
ccx q[43],q[14],q[44];
rz(0.15) q[44];
ccx q[43],q[14],q[44];
ccx q[42],q[11],q[43];
ccx q[3],q[8],q[42];
x q[3];
x q[8];
x q[11];
x q[14];

// blue K4 (4, 9, 13, 28) -> avoid connecting to none of the four
x q[3];
x q[8];
x q[12];
x q[27];
ccx q[3],q[8],q[42];
ccx q[42],q[12],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[12],q[43];
ccx q[3],q[8],q[42];
x q[3];
x q[8];
x q[12];
x q[27];

// blue K4 (4, 9, 13, 41) -> avoid connecting to none of the four
x q[3];
x q[8];
x q[12];
x q[40];
ccx q[3],q[8],q[42];
ccx q[42],q[12],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[12],q[43];
ccx q[3],q[8],q[42];
x q[3];
x q[8];
x q[12];
x q[40];

// blue K4 (4, 9, 15, 41) -> avoid connecting to none of the four
x q[3];
x q[8];
x q[14];
x q[40];
ccx q[3],q[8],q[42];
ccx q[42],q[14],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[14],q[43];
ccx q[3],q[8],q[42];
x q[3];
x q[8];
x q[14];
x q[40];

// blue K4 (4, 10, 13, 19) -> avoid connecting to none of the four
x q[3];
x q[9];
x q[12];
x q[18];
ccx q[3],q[9],q[42];
ccx q[42],q[12],q[43];
ccx q[43],q[18],q[44];
rz(0.15) q[44];
ccx q[43],q[18],q[44];
ccx q[42],q[12],q[43];
ccx q[3],q[9],q[42];
x q[3];
x q[9];
x q[12];
x q[18];

// blue K4 (4, 10, 13, 21) -> avoid connecting to none of the four
x q[3];
x q[9];
x q[12];
x q[20];
ccx q[3],q[9],q[42];
ccx q[42],q[12],q[43];
ccx q[43],q[20],q[44];
rz(0.15) q[44];
ccx q[43],q[20],q[44];
ccx q[42],q[12],q[43];
ccx q[3],q[9],q[42];
x q[3];
x q[9];
x q[12];
x q[20];

// blue K4 (4, 10, 15, 19) -> avoid connecting to none of the four
x q[3];
x q[9];
x q[14];
x q[18];
ccx q[3],q[9],q[42];
ccx q[42],q[14],q[43];
ccx q[43],q[18],q[44];
rz(0.15) q[44];
ccx q[43],q[18],q[44];
ccx q[42],q[14],q[43];
ccx q[3],q[9],q[42];
x q[3];
x q[9];
x q[14];
x q[18];

// blue K4 (4, 10, 15, 21) -> avoid connecting to none of the four
x q[3];
x q[9];
x q[14];
x q[20];
ccx q[3],q[9],q[42];
ccx q[42],q[14],q[43];
ccx q[43],q[20],q[44];
rz(0.15) q[44];
ccx q[43],q[20],q[44];
ccx q[42],q[14],q[43];
ccx q[3],q[9],q[42];
x q[3];
x q[9];
x q[14];
x q[20];

// blue K4 (4, 10, 19, 36) -> avoid connecting to none of the four
x q[3];
x q[9];
x q[18];
x q[35];
ccx q[3],q[9],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[18],q[43];
ccx q[3],q[9],q[42];
x q[3];
x q[9];
x q[18];
x q[35];

// blue K4 (4, 10, 19, 38) -> avoid connecting to none of the four
x q[3];
x q[9];
x q[18];
x q[37];
ccx q[3],q[9],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[18],q[43];
ccx q[3],q[9],q[42];
x q[3];
x q[9];
x q[18];
x q[37];

// blue K4 (4, 10, 21, 36) -> avoid connecting to none of the four
x q[3];
x q[9];
x q[20];
x q[35];
ccx q[3],q[9],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[20],q[43];
ccx q[3],q[9],q[42];
x q[3];
x q[9];
x q[20];
x q[35];

// blue K4 (4, 10, 21, 38) -> avoid connecting to none of the four
x q[3];
x q[9];
x q[20];
x q[37];
ccx q[3],q[9],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[20],q[43];
ccx q[3],q[9],q[42];
x q[3];
x q[9];
x q[20];
x q[37];

// blue K4 (4, 10, 36, 42) -> avoid connecting to none of the four
x q[3];
x q[9];
x q[35];
x q[41];
ccx q[3],q[9],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[35],q[43];
ccx q[3],q[9],q[42];
x q[3];
x q[9];
x q[35];
x q[41];

// blue K4 (4, 10, 38, 42) -> avoid connecting to none of the four
x q[3];
x q[9];
x q[37];
x q[41];
ccx q[3],q[9],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[37],q[43];
ccx q[3],q[9],q[42];
x q[3];
x q[9];
x q[37];
x q[41];

// blue K4 (4, 12, 15, 21) -> avoid connecting to none of the four
x q[3];
x q[11];
x q[14];
x q[20];
ccx q[3],q[11],q[42];
ccx q[42],q[14],q[43];
ccx q[43],q[20],q[44];
rz(0.15) q[44];
ccx q[43],q[20],q[44];
ccx q[42],q[14],q[43];
ccx q[3],q[11],q[42];
x q[3];
x q[11];
x q[14];
x q[20];

// blue K4 (4, 12, 15, 23) -> avoid connecting to none of the four
x q[3];
x q[11];
x q[14];
x q[22];
ccx q[3],q[11],q[42];
ccx q[42],q[14],q[43];
ccx q[43],q[22],q[44];
rz(0.15) q[44];
ccx q[43],q[22],q[44];
ccx q[42],q[14],q[43];
ccx q[3],q[11],q[42];
x q[3];
x q[11];
x q[14];
x q[22];

// blue K4 (4, 12, 21, 36) -> avoid connecting to none of the four
x q[3];
x q[11];
x q[20];
x q[35];
ccx q[3],q[11],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[20],q[43];
ccx q[3],q[11],q[42];
x q[3];
x q[11];
x q[20];
x q[35];

// blue K4 (4, 12, 21, 38) -> avoid connecting to none of the four
x q[3];
x q[11];
x q[20];
x q[37];
ccx q[3],q[11],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[20],q[43];
ccx q[3],q[11],q[42];
x q[3];
x q[11];
x q[20];
x q[37];

// blue K4 (4, 12, 23, 38) -> avoid connecting to none of the four
x q[3];
x q[11];
x q[22];
x q[37];
ccx q[3],q[11],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[22],q[43];
ccx q[3],q[11],q[42];
x q[3];
x q[11];
x q[22];
x q[37];

// blue K4 (4, 13, 19, 28) -> avoid connecting to none of the four
x q[3];
x q[12];
x q[18];
x q[27];
ccx q[3],q[12],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[18],q[43];
ccx q[3],q[12],q[42];
x q[3];
x q[12];
x q[18];
x q[27];

// blue K4 (4, 13, 19, 30) -> avoid connecting to none of the four
x q[3];
x q[12];
x q[18];
x q[29];
ccx q[3],q[12],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[18],q[43];
ccx q[3],q[12],q[42];
x q[3];
x q[12];
x q[18];
x q[29];

// blue K4 (4, 13, 21, 30) -> avoid connecting to none of the four
x q[3];
x q[12];
x q[20];
x q[29];
ccx q[3],q[12],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[20],q[43];
ccx q[3],q[12],q[42];
x q[3];
x q[12];
x q[20];
x q[29];

// blue K4 (4, 13, 21, 32) -> avoid connecting to none of the four
x q[3];
x q[12];
x q[20];
x q[31];
ccx q[3],q[12],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[20],q[43];
ccx q[3],q[12],q[42];
x q[3];
x q[12];
x q[20];
x q[31];

// blue K4 (4, 13, 28, 32) -> avoid connecting to none of the four
x q[3];
x q[12];
x q[27];
x q[31];
ccx q[3],q[12],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[27],q[43];
ccx q[3],q[12],q[42];
x q[3];
x q[12];
x q[27];
x q[31];

// blue K4 (4, 13, 28, 39) -> avoid connecting to none of the four
x q[3];
x q[12];
x q[27];
x q[38];
ccx q[3],q[12],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[27],q[43];
ccx q[3],q[12],q[42];
x q[3];
x q[12];
x q[27];
x q[38];

// blue K4 (4, 13, 30, 39) -> avoid connecting to none of the four
x q[3];
x q[12];
x q[29];
x q[38];
ccx q[3],q[12],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[29],q[43];
ccx q[3],q[12],q[42];
x q[3];
x q[12];
x q[29];
x q[38];

// blue K4 (4, 13, 30, 41) -> avoid connecting to none of the four
x q[3];
x q[12];
x q[29];
x q[40];
ccx q[3],q[12],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[29],q[43];
ccx q[3],q[12],q[42];
x q[3];
x q[12];
x q[29];
x q[40];

// blue K4 (4, 13, 32, 41) -> avoid connecting to none of the four
x q[3];
x q[12];
x q[31];
x q[40];
ccx q[3],q[12],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[31],q[43];
ccx q[3],q[12],q[42];
x q[3];
x q[12];
x q[31];
x q[40];

// blue K4 (4, 15, 19, 23) -> avoid connecting to none of the four
x q[3];
x q[14];
x q[18];
x q[22];
ccx q[3],q[14],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[22],q[44];
rz(0.15) q[44];
ccx q[43],q[22],q[44];
ccx q[42],q[18],q[43];
ccx q[3],q[14],q[42];
x q[3];
x q[14];
x q[18];
x q[22];

// blue K4 (4, 15, 19, 30) -> avoid connecting to none of the four
x q[3];
x q[14];
x q[18];
x q[29];
ccx q[3],q[14],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[18],q[43];
ccx q[3],q[14],q[42];
x q[3];
x q[14];
x q[18];
x q[29];

// blue K4 (4, 15, 21, 30) -> avoid connecting to none of the four
x q[3];
x q[14];
x q[20];
x q[29];
ccx q[3],q[14],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[20],q[43];
ccx q[3],q[14],q[42];
x q[3];
x q[14];
x q[20];
x q[29];

// blue K4 (4, 15, 21, 32) -> avoid connecting to none of the four
x q[3];
x q[14];
x q[20];
x q[31];
ccx q[3],q[14],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[20],q[43];
ccx q[3],q[14],q[42];
x q[3];
x q[14];
x q[20];
x q[31];

// blue K4 (4, 15, 23, 32) -> avoid connecting to none of the four
x q[3];
x q[14];
x q[22];
x q[31];
ccx q[3],q[14],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[22],q[43];
ccx q[3],q[14],q[42];
x q[3];
x q[14];
x q[22];
x q[31];

// blue K4 (4, 15, 30, 39) -> avoid connecting to none of the four
x q[3];
x q[14];
x q[29];
x q[38];
ccx q[3],q[14],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[29],q[43];
ccx q[3],q[14],q[42];
x q[3];
x q[14];
x q[29];
x q[38];

// blue K4 (4, 15, 30, 41) -> avoid connecting to none of the four
x q[3];
x q[14];
x q[29];
x q[40];
ccx q[3],q[14],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[29],q[43];
ccx q[3],q[14],q[42];
x q[3];
x q[14];
x q[29];
x q[40];

// blue K4 (4, 15, 32, 41) -> avoid connecting to none of the four
x q[3];
x q[14];
x q[31];
x q[40];
ccx q[3],q[14],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[31],q[43];
ccx q[3],q[14],q[42];
x q[3];
x q[14];
x q[31];
x q[40];

// blue K4 (4, 19, 23, 28) -> avoid connecting to none of the four
x q[3];
x q[18];
x q[22];
x q[27];
ccx q[3],q[18],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[22],q[43];
ccx q[3],q[18],q[42];
x q[3];
x q[18];
x q[22];
x q[27];

// blue K4 (4, 19, 23, 38) -> avoid connecting to none of the four
x q[3];
x q[18];
x q[22];
x q[37];
ccx q[3],q[18],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[22],q[43];
ccx q[3],q[18],q[42];
x q[3];
x q[18];
x q[22];
x q[37];

// blue K4 (4, 19, 28, 36) -> avoid connecting to none of the four
x q[3];
x q[18];
x q[27];
x q[35];
ccx q[3],q[18],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[27],q[43];
ccx q[3],q[18],q[42];
x q[3];
x q[18];
x q[27];
x q[35];

// blue K4 (4, 19, 30, 36) -> avoid connecting to none of the four
x q[3];
x q[18];
x q[29];
x q[35];
ccx q[3],q[18],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[29],q[43];
ccx q[3],q[18],q[42];
x q[3];
x q[18];
x q[29];
x q[35];

// blue K4 (4, 19, 30, 38) -> avoid connecting to none of the four
x q[3];
x q[18];
x q[29];
x q[37];
ccx q[3],q[18],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[29],q[43];
ccx q[3],q[18],q[42];
x q[3];
x q[18];
x q[29];
x q[37];

// blue K4 (4, 21, 30, 36) -> avoid connecting to none of the four
x q[3];
x q[20];
x q[29];
x q[35];
ccx q[3],q[20],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[29],q[43];
ccx q[3],q[20],q[42];
x q[3];
x q[20];
x q[29];
x q[35];

// blue K4 (4, 21, 30, 38) -> avoid connecting to none of the four
x q[3];
x q[20];
x q[29];
x q[37];
ccx q[3],q[20],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[29],q[43];
ccx q[3],q[20],q[42];
x q[3];
x q[20];
x q[29];
x q[37];

// blue K4 (4, 21, 32, 36) -> avoid connecting to none of the four
x q[3];
x q[20];
x q[31];
x q[35];
ccx q[3],q[20],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[31],q[43];
ccx q[3],q[20],q[42];
x q[3];
x q[20];
x q[31];
x q[35];

// blue K4 (4, 21, 32, 38) -> avoid connecting to none of the four
x q[3];
x q[20];
x q[31];
x q[37];
ccx q[3],q[20],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[31],q[43];
ccx q[3],q[20],q[42];
x q[3];
x q[20];
x q[31];
x q[37];

// blue K4 (4, 23, 28, 32) -> avoid connecting to none of the four
x q[3];
x q[22];
x q[27];
x q[31];
ccx q[3],q[22],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[27],q[43];
ccx q[3],q[22],q[42];
x q[3];
x q[22];
x q[27];
x q[31];

// blue K4 (4, 23, 32, 38) -> avoid connecting to none of the four
x q[3];
x q[22];
x q[31];
x q[37];
ccx q[3],q[22],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[31],q[43];
ccx q[3],q[22],q[42];
x q[3];
x q[22];
x q[31];
x q[37];

// blue K4 (4, 23, 38, 42) -> avoid connecting to none of the four
x q[3];
x q[22];
x q[37];
x q[41];
ccx q[3],q[22],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[37],q[43];
ccx q[3],q[22],q[42];
x q[3];
x q[22];
x q[37];
x q[41];

// blue K4 (4, 28, 32, 36) -> avoid connecting to none of the four
x q[3];
x q[27];
x q[31];
x q[35];
ccx q[3],q[27],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[31],q[43];
ccx q[3],q[27],q[42];
x q[3];
x q[27];
x q[31];
x q[35];

// blue K4 (4, 28, 36, 39) -> avoid connecting to none of the four
x q[3];
x q[27];
x q[35];
x q[38];
ccx q[3],q[27],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[35],q[43];
ccx q[3],q[27],q[42];
x q[3];
x q[27];
x q[35];
x q[38];

// blue K4 (4, 30, 36, 39) -> avoid connecting to none of the four
x q[3];
x q[29];
x q[35];
x q[38];
ccx q[3],q[29],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[35],q[43];
ccx q[3],q[29],q[42];
x q[3];
x q[29];
x q[35];
x q[38];

// blue K4 (4, 30, 36, 41) -> avoid connecting to none of the four
x q[3];
x q[29];
x q[35];
x q[40];
ccx q[3],q[29],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[35],q[43];
ccx q[3],q[29],q[42];
x q[3];
x q[29];
x q[35];
x q[40];

// blue K4 (4, 30, 38, 41) -> avoid connecting to none of the four
x q[3];
x q[29];
x q[37];
x q[40];
ccx q[3],q[29],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[37],q[43];
ccx q[3],q[29],q[42];
x q[3];
x q[29];
x q[37];
x q[40];

// blue K4 (4, 32, 36, 41) -> avoid connecting to none of the four
x q[3];
x q[31];
x q[35];
x q[40];
ccx q[3],q[31],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[35],q[43];
ccx q[3],q[31],q[42];
x q[3];
x q[31];
x q[35];
x q[40];

// blue K4 (4, 32, 38, 41) -> avoid connecting to none of the four
x q[3];
x q[31];
x q[37];
x q[40];
ccx q[3],q[31],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[37],q[43];
ccx q[3],q[31],q[42];
x q[3];
x q[31];
x q[37];
x q[40];

// blue K4 (4, 36, 39, 42) -> avoid connecting to none of the four
x q[3];
x q[35];
x q[38];
x q[41];
ccx q[3],q[35],q[42];
ccx q[42],q[38],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[38],q[43];
ccx q[3],q[35],q[42];
x q[3];
x q[35];
x q[38];
x q[41];

// blue K4 (4, 36, 41, 42) -> avoid connecting to none of the four
x q[3];
x q[35];
x q[40];
x q[41];
ccx q[3],q[35],q[42];
ccx q[42],q[40],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[40],q[43];
ccx q[3],q[35],q[42];
x q[3];
x q[35];
x q[40];
x q[41];

// blue K4 (4, 38, 41, 42) -> avoid connecting to none of the four
x q[3];
x q[37];
x q[40];
x q[41];
ccx q[3],q[37],q[42];
ccx q[42],q[40],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[40],q[43];
ccx q[3],q[37],q[42];
x q[3];
x q[37];
x q[40];
x q[41];

// blue K4 (5, 6, 9, 14) -> avoid connecting to none of the four
x q[4];
x q[5];
x q[8];
x q[13];
ccx q[4],q[5],q[42];
ccx q[42],q[8],q[43];
ccx q[43],q[13],q[44];
rz(0.15) q[44];
ccx q[43],q[13],q[44];
ccx q[42],q[8],q[43];
ccx q[4],q[5],q[42];
x q[4];
x q[5];
x q[8];
x q[13];

// blue K4 (5, 6, 10, 14) -> avoid connecting to none of the four
x q[4];
x q[5];
x q[9];
x q[13];
ccx q[4],q[5],q[42];
ccx q[42],q[9],q[43];
ccx q[43],q[13],q[44];
rz(0.15) q[44];
ccx q[43],q[13],q[44];
ccx q[42],q[9],q[43];
ccx q[4],q[5],q[42];
x q[4];
x q[5];
x q[9];
x q[13];

// blue K4 (5, 6, 11, 14) -> avoid connecting to none of the four
x q[4];
x q[5];
x q[10];
x q[13];
ccx q[4],q[5],q[42];
ccx q[42],q[10],q[43];
ccx q[43],q[13],q[44];
rz(0.15) q[44];
ccx q[43],q[13],q[44];
ccx q[42],q[10],q[43];
ccx q[4],q[5],q[42];
x q[4];
x q[5];
x q[10];
x q[13];

// blue K4 (5, 6, 14, 40) -> avoid connecting to none of the four
x q[4];
x q[5];
x q[13];
x q[39];
ccx q[4],q[5],q[42];
ccx q[42],q[13],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[13],q[43];
ccx q[4],q[5],q[42];
x q[4];
x q[5];
x q[13];
x q[39];

// blue K4 (5, 8, 11, 14) -> avoid connecting to none of the four
x q[4];
x q[7];
x q[10];
x q[13];
ccx q[4],q[7],q[42];
ccx q[42],q[10],q[43];
ccx q[43],q[13],q[44];
rz(0.15) q[44];
ccx q[43],q[13],q[44];
ccx q[42],q[10],q[43];
ccx q[4],q[7],q[42];
x q[4];
x q[7];
x q[10];
x q[13];

// blue K4 (5, 8, 11, 16) -> avoid connecting to none of the four
x q[4];
x q[7];
x q[10];
x q[15];
ccx q[4],q[7],q[42];
ccx q[42],q[10],q[43];
ccx q[43],q[15],q[44];
rz(0.15) q[44];
ccx q[43],q[15],q[44];
ccx q[42],q[10],q[43];
ccx q[4],q[7],q[42];
x q[4];
x q[7];
x q[10];
x q[15];

// blue K4 (5, 8, 13, 14) -> avoid connecting to none of the four
x q[4];
x q[7];
x q[12];
x q[13];
ccx q[4],q[7],q[42];
ccx q[42],q[12],q[43];
ccx q[43],q[13],q[44];
rz(0.15) q[44];
ccx q[43],q[13],q[44];
ccx q[42],q[12],q[43];
ccx q[4],q[7],q[42];
x q[4];
x q[7];
x q[12];
x q[13];

// blue K4 (5, 8, 13, 16) -> avoid connecting to none of the four
x q[4];
x q[7];
x q[12];
x q[15];
ccx q[4],q[7],q[42];
ccx q[42],q[12],q[43];
ccx q[43],q[15],q[44];
rz(0.15) q[44];
ccx q[43],q[15],q[44];
ccx q[42],q[12],q[43];
ccx q[4],q[7],q[42];
x q[4];
x q[7];
x q[12];
x q[15];

// blue K4 (5, 8, 14, 40) -> avoid connecting to none of the four
x q[4];
x q[7];
x q[13];
x q[39];
ccx q[4],q[7],q[42];
ccx q[42],q[13],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[13],q[43];
ccx q[4],q[7],q[42];
x q[4];
x q[7];
x q[13];
x q[39];

// blue K4 (5, 8, 14, 42) -> avoid connecting to none of the four
x q[4];
x q[7];
x q[13];
x q[41];
ccx q[4],q[7],q[42];
ccx q[42],q[13],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[13],q[43];
ccx q[4],q[7],q[42];
x q[4];
x q[7];
x q[13];
x q[41];

// blue K4 (5, 8, 16, 40) -> avoid connecting to none of the four
x q[4];
x q[7];
x q[15];
x q[39];
ccx q[4],q[7],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[15],q[43];
ccx q[4],q[7],q[42];
x q[4];
x q[7];
x q[15];
x q[39];

// blue K4 (5, 8, 16, 42) -> avoid connecting to none of the four
x q[4];
x q[7];
x q[15];
x q[41];
ccx q[4],q[7],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[15],q[43];
ccx q[4],q[7],q[42];
x q[4];
x q[7];
x q[15];
x q[41];

// blue K4 (5, 9, 13, 14) -> avoid connecting to none of the four
x q[4];
x q[8];
x q[12];
x q[13];
ccx q[4],q[8],q[42];
ccx q[42],q[12],q[43];
ccx q[43],q[13],q[44];
rz(0.15) q[44];
ccx q[43],q[13],q[44];
ccx q[42],q[12],q[43];
ccx q[4],q[8],q[42];
x q[4];
x q[8];
x q[12];
x q[13];

// blue K4 (5, 9, 13, 24) -> avoid connecting to none of the four
x q[4];
x q[8];
x q[12];
x q[23];
ccx q[4],q[8],q[42];
ccx q[42],q[12],q[43];
ccx q[43],q[23],q[44];
rz(0.15) q[44];
ccx q[43],q[23],q[44];
ccx q[42],q[12],q[43];
ccx q[4],q[8],q[42];
x q[4];
x q[8];
x q[12];
x q[23];

// blue K4 (5, 9, 13, 37) -> avoid connecting to none of the four
x q[4];
x q[8];
x q[12];
x q[36];
ccx q[4],q[8],q[42];
ccx q[42],q[12],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[12],q[43];
ccx q[4],q[8],q[42];
x q[4];
x q[8];
x q[12];
x q[36];

// blue K4 (5, 9, 14, 20) -> avoid connecting to none of the four
x q[4];
x q[8];
x q[13];
x q[19];
ccx q[4],q[8],q[42];
ccx q[42],q[13],q[43];
ccx q[43],q[19],q[44];
rz(0.15) q[44];
ccx q[43],q[19],q[44];
ccx q[42],q[13],q[43];
ccx q[4],q[8],q[42];
x q[4];
x q[8];
x q[13];
x q[19];

// blue K4 (5, 9, 14, 33) -> avoid connecting to none of the four
x q[4];
x q[8];
x q[13];
x q[32];
ccx q[4],q[8],q[42];
ccx q[42],q[13],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[13],q[43];
ccx q[4],q[8],q[42];
x q[4];
x q[8];
x q[13];
x q[32];

// blue K4 (5, 9, 20, 24) -> avoid connecting to none of the four
x q[4];
x q[8];
x q[19];
x q[23];
ccx q[4],q[8],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[23],q[44];
rz(0.15) q[44];
ccx q[43],q[23],q[44];
ccx q[42],q[19],q[43];
ccx q[4],q[8],q[42];
x q[4];
x q[8];
x q[19];
x q[23];

// blue K4 (5, 9, 20, 37) -> avoid connecting to none of the four
x q[4];
x q[8];
x q[19];
x q[36];
ccx q[4],q[8],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[19],q[43];
ccx q[4],q[8],q[42];
x q[4];
x q[8];
x q[19];
x q[36];

// blue K4 (5, 9, 24, 33) -> avoid connecting to none of the four
x q[4];
x q[8];
x q[23];
x q[32];
ccx q[4],q[8],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[23],q[43];
ccx q[4],q[8],q[42];
x q[4];
x q[8];
x q[23];
x q[32];

// blue K4 (5, 9, 33, 37) -> avoid connecting to none of the four
x q[4];
x q[8];
x q[32];
x q[36];
ccx q[4],q[8],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[32],q[43];
ccx q[4],q[8],q[42];
x q[4];
x q[8];
x q[32];
x q[36];

// blue K4 (5, 10, 13, 14) -> avoid connecting to none of the four
x q[4];
x q[9];
x q[12];
x q[13];
ccx q[4],q[9],q[42];
ccx q[42],q[12],q[43];
ccx q[43],q[13],q[44];
rz(0.15) q[44];
ccx q[43],q[13],q[44];
ccx q[42],q[12],q[43];
ccx q[4],q[9],q[42];
x q[4];
x q[9];
x q[12];
x q[13];

// blue K4 (5, 10, 13, 16) -> avoid connecting to none of the four
x q[4];
x q[9];
x q[12];
x q[15];
ccx q[4],q[9],q[42];
ccx q[42],q[12],q[43];
ccx q[43],q[15],q[44];
rz(0.15) q[44];
ccx q[43],q[15],q[44];
ccx q[42],q[12],q[43];
ccx q[4],q[9],q[42];
x q[4];
x q[9];
x q[12];
x q[15];

// blue K4 (5, 10, 14, 29) -> avoid connecting to none of the four
x q[4];
x q[9];
x q[13];
x q[28];
ccx q[4],q[9],q[42];
ccx q[42],q[13],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[13],q[43];
ccx q[4],q[9],q[42];
x q[4];
x q[9];
x q[13];
x q[28];

// blue K4 (5, 10, 14, 42) -> avoid connecting to none of the four
x q[4];
x q[9];
x q[13];
x q[41];
ccx q[4],q[9],q[42];
ccx q[42],q[13],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[13],q[43];
ccx q[4],q[9],q[42];
x q[4];
x q[9];
x q[13];
x q[41];

// blue K4 (5, 10, 16, 42) -> avoid connecting to none of the four
x q[4];
x q[9];
x q[15];
x q[41];
ccx q[4],q[9],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[15],q[43];
ccx q[4],q[9],q[42];
x q[4];
x q[9];
x q[15];
x q[41];

// blue K4 (5, 11, 14, 20) -> avoid connecting to none of the four
x q[4];
x q[10];
x q[13];
x q[19];
ccx q[4],q[10],q[42];
ccx q[42],q[13],q[43];
ccx q[43],q[19],q[44];
rz(0.15) q[44];
ccx q[43],q[19],q[44];
ccx q[42],q[13],q[43];
ccx q[4],q[10],q[42];
x q[4];
x q[10];
x q[13];
x q[19];

// blue K4 (5, 11, 14, 22) -> avoid connecting to none of the four
x q[4];
x q[10];
x q[13];
x q[21];
ccx q[4],q[10],q[42];
ccx q[42],q[13],q[43];
ccx q[43],q[21],q[44];
rz(0.15) q[44];
ccx q[43],q[21],q[44];
ccx q[42],q[13],q[43];
ccx q[4],q[10],q[42];
x q[4];
x q[10];
x q[13];
x q[21];

// blue K4 (5, 11, 16, 20) -> avoid connecting to none of the four
x q[4];
x q[10];
x q[15];
x q[19];
ccx q[4],q[10],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[19],q[44];
rz(0.15) q[44];
ccx q[43],q[19],q[44];
ccx q[42],q[15],q[43];
ccx q[4],q[10],q[42];
x q[4];
x q[10];
x q[15];
x q[19];

// blue K4 (5, 11, 16, 22) -> avoid connecting to none of the four
x q[4];
x q[10];
x q[15];
x q[21];
ccx q[4],q[10],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[21],q[44];
rz(0.15) q[44];
ccx q[43],q[21],q[44];
ccx q[42],q[15],q[43];
ccx q[4],q[10],q[42];
x q[4];
x q[10];
x q[15];
x q[21];

// blue K4 (5, 11, 20, 37) -> avoid connecting to none of the four
x q[4];
x q[10];
x q[19];
x q[36];
ccx q[4],q[10],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[19],q[43];
ccx q[4],q[10],q[42];
x q[4];
x q[10];
x q[19];
x q[36];

// blue K4 (5, 11, 20, 39) -> avoid connecting to none of the four
x q[4];
x q[10];
x q[19];
x q[38];
ccx q[4],q[10],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[19],q[43];
ccx q[4],q[10],q[42];
x q[4];
x q[10];
x q[19];
x q[38];

// blue K4 (5, 11, 22, 37) -> avoid connecting to none of the four
x q[4];
x q[10];
x q[21];
x q[36];
ccx q[4],q[10],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[21],q[43];
ccx q[4],q[10],q[42];
x q[4];
x q[10];
x q[21];
x q[36];

// blue K4 (5, 11, 22, 39) -> avoid connecting to none of the four
x q[4];
x q[10];
x q[21];
x q[38];
ccx q[4],q[10],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[21],q[43];
ccx q[4],q[10],q[42];
x q[4];
x q[10];
x q[21];
x q[38];

// blue K4 (5, 13, 14, 22) -> avoid connecting to none of the four
x q[4];
x q[12];
x q[13];
x q[21];
ccx q[4],q[12],q[42];
ccx q[42],q[13],q[43];
ccx q[43],q[21],q[44];
rz(0.15) q[44];
ccx q[43],q[21],q[44];
ccx q[42],q[13],q[43];
ccx q[4],q[12],q[42];
x q[4];
x q[12];
x q[13];
x q[21];

// blue K4 (5, 13, 16, 22) -> avoid connecting to none of the four
x q[4];
x q[12];
x q[15];
x q[21];
ccx q[4],q[12],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[21],q[44];
rz(0.15) q[44];
ccx q[43],q[21],q[44];
ccx q[42],q[15],q[43];
ccx q[4],q[12],q[42];
x q[4];
x q[12];
x q[15];
x q[21];

// blue K4 (5, 13, 16, 24) -> avoid connecting to none of the four
x q[4];
x q[12];
x q[15];
x q[23];
ccx q[4],q[12],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[23],q[44];
rz(0.15) q[44];
ccx q[43],q[23],q[44];
ccx q[42],q[15],q[43];
ccx q[4],q[12],q[42];
x q[4];
x q[12];
x q[15];
x q[23];

// blue K4 (5, 13, 22, 37) -> avoid connecting to none of the four
x q[4];
x q[12];
x q[21];
x q[36];
ccx q[4],q[12],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[21],q[43];
ccx q[4],q[12],q[42];
x q[4];
x q[12];
x q[21];
x q[36];

// blue K4 (5, 13, 22, 39) -> avoid connecting to none of the four
x q[4];
x q[12];
x q[21];
x q[38];
ccx q[4],q[12],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[21],q[43];
ccx q[4],q[12],q[42];
x q[4];
x q[12];
x q[21];
x q[38];

// blue K4 (5, 13, 24, 39) -> avoid connecting to none of the four
x q[4];
x q[12];
x q[23];
x q[38];
ccx q[4],q[12],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[23],q[43];
ccx q[4],q[12],q[42];
x q[4];
x q[12];
x q[23];
x q[38];

// blue K4 (5, 14, 20, 29) -> avoid connecting to none of the four
x q[4];
x q[13];
x q[19];
x q[28];
ccx q[4],q[13],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[19],q[43];
ccx q[4],q[13],q[42];
x q[4];
x q[13];
x q[19];
x q[28];

// blue K4 (5, 14, 20, 31) -> avoid connecting to none of the four
x q[4];
x q[13];
x q[19];
x q[30];
ccx q[4],q[13],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[19],q[43];
ccx q[4],q[13],q[42];
x q[4];
x q[13];
x q[19];
x q[30];

// blue K4 (5, 14, 22, 31) -> avoid connecting to none of the four
x q[4];
x q[13];
x q[21];
x q[30];
ccx q[4],q[13],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[21],q[43];
ccx q[4],q[13],q[42];
x q[4];
x q[13];
x q[21];
x q[30];

// blue K4 (5, 14, 22, 33) -> avoid connecting to none of the four
x q[4];
x q[13];
x q[21];
x q[32];
ccx q[4],q[13],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[21],q[43];
ccx q[4],q[13],q[42];
x q[4];
x q[13];
x q[21];
x q[32];

// blue K4 (5, 14, 29, 33) -> avoid connecting to none of the four
x q[4];
x q[13];
x q[28];
x q[32];
ccx q[4],q[13],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[28],q[43];
ccx q[4],q[13],q[42];
x q[4];
x q[13];
x q[28];
x q[32];

// blue K4 (5, 14, 29, 40) -> avoid connecting to none of the four
x q[4];
x q[13];
x q[28];
x q[39];
ccx q[4],q[13],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[28],q[43];
ccx q[4],q[13],q[42];
x q[4];
x q[13];
x q[28];
x q[39];

// blue K4 (5, 14, 31, 40) -> avoid connecting to none of the four
x q[4];
x q[13];
x q[30];
x q[39];
ccx q[4],q[13],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[30],q[43];
ccx q[4],q[13],q[42];
x q[4];
x q[13];
x q[30];
x q[39];

// blue K4 (5, 14, 31, 42) -> avoid connecting to none of the four
x q[4];
x q[13];
x q[30];
x q[41];
ccx q[4],q[13],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[30],q[43];
ccx q[4],q[13],q[42];
x q[4];
x q[13];
x q[30];
x q[41];

// blue K4 (5, 14, 33, 42) -> avoid connecting to none of the four
x q[4];
x q[13];
x q[32];
x q[41];
ccx q[4],q[13],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[32],q[43];
ccx q[4],q[13],q[42];
x q[4];
x q[13];
x q[32];
x q[41];

// blue K4 (5, 16, 20, 24) -> avoid connecting to none of the four
x q[4];
x q[15];
x q[19];
x q[23];
ccx q[4],q[15],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[23],q[44];
rz(0.15) q[44];
ccx q[43],q[23],q[44];
ccx q[42],q[19],q[43];
ccx q[4],q[15],q[42];
x q[4];
x q[15];
x q[19];
x q[23];

// blue K4 (5, 16, 20, 31) -> avoid connecting to none of the four
x q[4];
x q[15];
x q[19];
x q[30];
ccx q[4],q[15],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[19],q[43];
ccx q[4],q[15],q[42];
x q[4];
x q[15];
x q[19];
x q[30];

// blue K4 (5, 16, 22, 31) -> avoid connecting to none of the four
x q[4];
x q[15];
x q[21];
x q[30];
ccx q[4],q[15],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[21],q[43];
ccx q[4],q[15],q[42];
x q[4];
x q[15];
x q[21];
x q[30];

// blue K4 (5, 16, 22, 33) -> avoid connecting to none of the four
x q[4];
x q[15];
x q[21];
x q[32];
ccx q[4],q[15],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[21],q[43];
ccx q[4],q[15],q[42];
x q[4];
x q[15];
x q[21];
x q[32];

// blue K4 (5, 16, 24, 33) -> avoid connecting to none of the four
x q[4];
x q[15];
x q[23];
x q[32];
ccx q[4],q[15],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[23],q[43];
ccx q[4],q[15],q[42];
x q[4];
x q[15];
x q[23];
x q[32];

// blue K4 (5, 16, 31, 40) -> avoid connecting to none of the four
x q[4];
x q[15];
x q[30];
x q[39];
ccx q[4],q[15],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[30],q[43];
ccx q[4],q[15],q[42];
x q[4];
x q[15];
x q[30];
x q[39];

// blue K4 (5, 16, 31, 42) -> avoid connecting to none of the four
x q[4];
x q[15];
x q[30];
x q[41];
ccx q[4],q[15],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[30],q[43];
ccx q[4],q[15],q[42];
x q[4];
x q[15];
x q[30];
x q[41];

// blue K4 (5, 16, 33, 42) -> avoid connecting to none of the four
x q[4];
x q[15];
x q[32];
x q[41];
ccx q[4],q[15],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[32],q[43];
ccx q[4],q[15],q[42];
x q[4];
x q[15];
x q[32];
x q[41];

// blue K4 (5, 20, 24, 29) -> avoid connecting to none of the four
x q[4];
x q[19];
x q[23];
x q[28];
ccx q[4],q[19],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[23],q[43];
ccx q[4],q[19],q[42];
x q[4];
x q[19];
x q[23];
x q[28];

// blue K4 (5, 20, 24, 39) -> avoid connecting to none of the four
x q[4];
x q[19];
x q[23];
x q[38];
ccx q[4],q[19],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[23],q[43];
ccx q[4],q[19],q[42];
x q[4];
x q[19];
x q[23];
x q[38];

// blue K4 (5, 20, 29, 37) -> avoid connecting to none of the four
x q[4];
x q[19];
x q[28];
x q[36];
ccx q[4],q[19],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[28],q[43];
ccx q[4],q[19],q[42];
x q[4];
x q[19];
x q[28];
x q[36];

// blue K4 (5, 20, 31, 37) -> avoid connecting to none of the four
x q[4];
x q[19];
x q[30];
x q[36];
ccx q[4],q[19],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[30],q[43];
ccx q[4],q[19],q[42];
x q[4];
x q[19];
x q[30];
x q[36];

// blue K4 (5, 20, 31, 39) -> avoid connecting to none of the four
x q[4];
x q[19];
x q[30];
x q[38];
ccx q[4],q[19],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[30],q[43];
ccx q[4],q[19],q[42];
x q[4];
x q[19];
x q[30];
x q[38];

// blue K4 (5, 22, 31, 37) -> avoid connecting to none of the four
x q[4];
x q[21];
x q[30];
x q[36];
ccx q[4],q[21],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[30],q[43];
ccx q[4],q[21],q[42];
x q[4];
x q[21];
x q[30];
x q[36];

// blue K4 (5, 22, 31, 39) -> avoid connecting to none of the four
x q[4];
x q[21];
x q[30];
x q[38];
ccx q[4],q[21],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[30],q[43];
ccx q[4],q[21],q[42];
x q[4];
x q[21];
x q[30];
x q[38];

// blue K4 (5, 22, 33, 37) -> avoid connecting to none of the four
x q[4];
x q[21];
x q[32];
x q[36];
ccx q[4],q[21],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[32],q[43];
ccx q[4],q[21],q[42];
x q[4];
x q[21];
x q[32];
x q[36];

// blue K4 (5, 22, 33, 39) -> avoid connecting to none of the four
x q[4];
x q[21];
x q[32];
x q[38];
ccx q[4],q[21],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[32],q[43];
ccx q[4],q[21],q[42];
x q[4];
x q[21];
x q[32];
x q[38];

// blue K4 (5, 24, 29, 33) -> avoid connecting to none of the four
x q[4];
x q[23];
x q[28];
x q[32];
ccx q[4],q[23],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[28],q[43];
ccx q[4],q[23],q[42];
x q[4];
x q[23];
x q[28];
x q[32];

// blue K4 (5, 24, 33, 39) -> avoid connecting to none of the four
x q[4];
x q[23];
x q[32];
x q[38];
ccx q[4],q[23],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[32],q[43];
ccx q[4],q[23],q[42];
x q[4];
x q[23];
x q[32];
x q[38];

// blue K4 (5, 29, 33, 37) -> avoid connecting to none of the four
x q[4];
x q[28];
x q[32];
x q[36];
ccx q[4],q[28],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[32],q[43];
ccx q[4],q[28],q[42];
x q[4];
x q[28];
x q[32];
x q[36];

// blue K4 (5, 29, 37, 40) -> avoid connecting to none of the four
x q[4];
x q[28];
x q[36];
x q[39];
ccx q[4],q[28],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[36],q[43];
ccx q[4],q[28],q[42];
x q[4];
x q[28];
x q[36];
x q[39];

// blue K4 (5, 31, 37, 40) -> avoid connecting to none of the four
x q[4];
x q[30];
x q[36];
x q[39];
ccx q[4],q[30],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[36],q[43];
ccx q[4],q[30],q[42];
x q[4];
x q[30];
x q[36];
x q[39];

// blue K4 (5, 31, 37, 42) -> avoid connecting to none of the four
x q[4];
x q[30];
x q[36];
x q[41];
ccx q[4],q[30],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[36],q[43];
ccx q[4],q[30],q[42];
x q[4];
x q[30];
x q[36];
x q[41];

// blue K4 (5, 31, 39, 40) -> avoid connecting to none of the four
x q[4];
x q[30];
x q[38];
x q[39];
ccx q[4],q[30],q[42];
ccx q[42],q[38],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[38],q[43];
ccx q[4],q[30],q[42];
x q[4];
x q[30];
x q[38];
x q[39];

// blue K4 (5, 31, 39, 42) -> avoid connecting to none of the four
x q[4];
x q[30];
x q[38];
x q[41];
ccx q[4],q[30],q[42];
ccx q[42],q[38],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[38],q[43];
ccx q[4],q[30],q[42];
x q[4];
x q[30];
x q[38];
x q[41];

// blue K4 (5, 33, 37, 42) -> avoid connecting to none of the four
x q[4];
x q[32];
x q[36];
x q[41];
ccx q[4],q[32],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[36],q[43];
ccx q[4],q[32],q[42];
x q[4];
x q[32];
x q[36];
x q[41];

// blue K4 (5, 33, 39, 42) -> avoid connecting to none of the four
x q[4];
x q[32];
x q[38];
x q[41];
ccx q[4],q[32],q[42];
ccx q[42],q[38],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[38],q[43];
ccx q[4],q[32],q[42];
x q[4];
x q[32];
x q[38];
x q[41];

// blue K4 (6, 7, 10, 15) -> avoid connecting to none of the four
x q[5];
x q[6];
x q[9];
x q[14];
ccx q[5],q[6],q[42];
ccx q[42],q[9],q[43];
ccx q[43],q[14],q[44];
rz(0.15) q[44];
ccx q[43],q[14],q[44];
ccx q[42],q[9],q[43];
ccx q[5],q[6],q[42];
x q[5];
x q[6];
x q[9];
x q[14];

// blue K4 (6, 7, 11, 15) -> avoid connecting to none of the four
x q[5];
x q[6];
x q[10];
x q[14];
ccx q[5],q[6],q[42];
ccx q[42],q[10],q[43];
ccx q[43],q[14],q[44];
rz(0.15) q[44];
ccx q[43],q[14],q[44];
ccx q[42],q[10],q[43];
ccx q[5],q[6],q[42];
x q[5];
x q[6];
x q[10];
x q[14];

// blue K4 (6, 7, 12, 15) -> avoid connecting to none of the four
x q[5];
x q[6];
x q[11];
x q[14];
ccx q[5],q[6],q[42];
ccx q[42],q[11],q[43];
ccx q[43],q[14],q[44];
rz(0.15) q[44];
ccx q[43],q[14],q[44];
ccx q[42],q[11],q[43];
ccx q[5],q[6],q[42];
x q[5];
x q[6];
x q[11];
x q[14];

// blue K4 (6, 7, 15, 41) -> avoid connecting to none of the four
x q[5];
x q[6];
x q[14];
x q[40];
ccx q[5],q[6],q[42];
ccx q[42],q[14],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[14],q[43];
ccx q[5],q[6],q[42];
x q[5];
x q[6];
x q[14];
x q[40];

// blue K4 (6, 9, 12, 15) -> avoid connecting to none of the four
x q[5];
x q[8];
x q[11];
x q[14];
ccx q[5],q[8],q[42];
ccx q[42],q[11],q[43];
ccx q[43],q[14],q[44];
rz(0.15) q[44];
ccx q[43],q[14],q[44];
ccx q[42],q[11],q[43];
ccx q[5],q[8],q[42];
x q[5];
x q[8];
x q[11];
x q[14];

// blue K4 (6, 9, 12, 17) -> avoid connecting to none of the four
x q[5];
x q[8];
x q[11];
x q[16];
ccx q[5],q[8],q[42];
ccx q[42],q[11],q[43];
ccx q[43],q[16],q[44];
rz(0.15) q[44];
ccx q[43],q[16],q[44];
ccx q[42],q[11],q[43];
ccx q[5],q[8],q[42];
x q[5];
x q[8];
x q[11];
x q[16];

// blue K4 (6, 9, 14, 15) -> avoid connecting to none of the four
x q[5];
x q[8];
x q[13];
x q[14];
ccx q[5],q[8],q[42];
ccx q[42],q[13],q[43];
ccx q[43],q[14],q[44];
rz(0.15) q[44];
ccx q[43],q[14],q[44];
ccx q[42],q[13],q[43];
ccx q[5],q[8],q[42];
x q[5];
x q[8];
x q[13];
x q[14];

// blue K4 (6, 9, 14, 17) -> avoid connecting to none of the four
x q[5];
x q[8];
x q[13];
x q[16];
ccx q[5],q[8],q[42];
ccx q[42],q[13],q[43];
ccx q[43],q[16],q[44];
rz(0.15) q[44];
ccx q[43],q[16],q[44];
ccx q[42],q[13],q[43];
ccx q[5],q[8],q[42];
x q[5];
x q[8];
x q[13];
x q[16];

// blue K4 (6, 9, 15, 41) -> avoid connecting to none of the four
x q[5];
x q[8];
x q[14];
x q[40];
ccx q[5],q[8],q[42];
ccx q[42],q[14],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[14],q[43];
ccx q[5],q[8],q[42];
x q[5];
x q[8];
x q[14];
x q[40];

// blue K4 (6, 9, 17, 41) -> avoid connecting to none of the four
x q[5];
x q[8];
x q[16];
x q[40];
ccx q[5],q[8],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[16],q[43];
ccx q[5],q[8],q[42];
x q[5];
x q[8];
x q[16];
x q[40];

// blue K4 (6, 10, 14, 15) -> avoid connecting to none of the four
x q[5];
x q[9];
x q[13];
x q[14];
ccx q[5],q[9],q[42];
ccx q[42],q[13],q[43];
ccx q[43],q[14],q[44];
rz(0.15) q[44];
ccx q[43],q[14],q[44];
ccx q[42],q[13],q[43];
ccx q[5],q[9],q[42];
x q[5];
x q[9];
x q[13];
x q[14];

// blue K4 (6, 10, 14, 25) -> avoid connecting to none of the four
x q[5];
x q[9];
x q[13];
x q[24];
ccx q[5],q[9],q[42];
ccx q[42],q[13],q[43];
ccx q[43],q[24],q[44];
rz(0.15) q[44];
ccx q[43],q[24],q[44];
ccx q[42],q[13],q[43];
ccx q[5],q[9],q[42];
x q[5];
x q[9];
x q[13];
x q[24];

// blue K4 (6, 10, 14, 38) -> avoid connecting to none of the four
x q[5];
x q[9];
x q[13];
x q[37];
ccx q[5],q[9],q[42];
ccx q[42],q[13],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[13],q[43];
ccx q[5],q[9],q[42];
x q[5];
x q[9];
x q[13];
x q[37];

// blue K4 (6, 10, 15, 21) -> avoid connecting to none of the four
x q[5];
x q[9];
x q[14];
x q[20];
ccx q[5],q[9],q[42];
ccx q[42],q[14],q[43];
ccx q[43],q[20],q[44];
rz(0.15) q[44];
ccx q[43],q[20],q[44];
ccx q[42],q[14],q[43];
ccx q[5],q[9],q[42];
x q[5];
x q[9];
x q[14];
x q[20];

// blue K4 (6, 10, 15, 34) -> avoid connecting to none of the four
x q[5];
x q[9];
x q[14];
x q[33];
ccx q[5],q[9],q[42];
ccx q[42],q[14],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[14],q[43];
ccx q[5],q[9],q[42];
x q[5];
x q[9];
x q[14];
x q[33];

// blue K4 (6, 10, 21, 25) -> avoid connecting to none of the four
x q[5];
x q[9];
x q[20];
x q[24];
ccx q[5],q[9],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[24],q[44];
rz(0.15) q[44];
ccx q[43],q[24],q[44];
ccx q[42],q[20],q[43];
ccx q[5],q[9],q[42];
x q[5];
x q[9];
x q[20];
x q[24];

// blue K4 (6, 10, 21, 38) -> avoid connecting to none of the four
x q[5];
x q[9];
x q[20];
x q[37];
ccx q[5],q[9],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[20],q[43];
ccx q[5],q[9],q[42];
x q[5];
x q[9];
x q[20];
x q[37];

// blue K4 (6, 10, 25, 34) -> avoid connecting to none of the four
x q[5];
x q[9];
x q[24];
x q[33];
ccx q[5],q[9],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[24],q[43];
ccx q[5],q[9],q[42];
x q[5];
x q[9];
x q[24];
x q[33];

// blue K4 (6, 10, 34, 38) -> avoid connecting to none of the four
x q[5];
x q[9];
x q[33];
x q[37];
ccx q[5],q[9],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[33],q[43];
ccx q[5],q[9],q[42];
x q[5];
x q[9];
x q[33];
x q[37];

// blue K4 (6, 11, 14, 15) -> avoid connecting to none of the four
x q[5];
x q[10];
x q[13];
x q[14];
ccx q[5],q[10],q[42];
ccx q[42],q[13],q[43];
ccx q[43],q[14],q[44];
rz(0.15) q[44];
ccx q[43],q[14],q[44];
ccx q[42],q[13],q[43];
ccx q[5],q[10],q[42];
x q[5];
x q[10];
x q[13];
x q[14];

// blue K4 (6, 11, 14, 17) -> avoid connecting to none of the four
x q[5];
x q[10];
x q[13];
x q[16];
ccx q[5],q[10],q[42];
ccx q[42],q[13],q[43];
ccx q[43],q[16],q[44];
rz(0.15) q[44];
ccx q[43],q[16],q[44];
ccx q[42],q[13],q[43];
ccx q[5],q[10],q[42];
x q[5];
x q[10];
x q[13];
x q[16];

// blue K4 (6, 11, 15, 30) -> avoid connecting to none of the four
x q[5];
x q[10];
x q[14];
x q[29];
ccx q[5],q[10],q[42];
ccx q[42],q[14],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[14],q[43];
ccx q[5],q[10],q[42];
x q[5];
x q[10];
x q[14];
x q[29];

// blue K4 (6, 11, 15, 32) -> avoid connecting to none of the four
x q[5];
x q[10];
x q[14];
x q[31];
ccx q[5],q[10],q[42];
ccx q[42],q[14],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[14],q[43];
ccx q[5],q[10],q[42];
x q[5];
x q[10];
x q[14];
x q[31];

// blue K4 (6, 11, 17, 32) -> avoid connecting to none of the four
x q[5];
x q[10];
x q[16];
x q[31];
ccx q[5],q[10],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[16],q[43];
ccx q[5],q[10],q[42];
x q[5];
x q[10];
x q[16];
x q[31];

// blue K4 (6, 12, 15, 21) -> avoid connecting to none of the four
x q[5];
x q[11];
x q[14];
x q[20];
ccx q[5],q[11],q[42];
ccx q[42],q[14],q[43];
ccx q[43],q[20],q[44];
rz(0.15) q[44];
ccx q[43],q[20],q[44];
ccx q[42],q[14],q[43];
ccx q[5],q[11],q[42];
x q[5];
x q[11];
x q[14];
x q[20];

// blue K4 (6, 12, 15, 23) -> avoid connecting to none of the four
x q[5];
x q[11];
x q[14];
x q[22];
ccx q[5],q[11],q[42];
ccx q[42],q[14],q[43];
ccx q[43],q[22],q[44];
rz(0.15) q[44];
ccx q[43],q[22],q[44];
ccx q[42],q[14],q[43];
ccx q[5],q[11],q[42];
x q[5];
x q[11];
x q[14];
x q[22];

// blue K4 (6, 12, 17, 21) -> avoid connecting to none of the four
x q[5];
x q[11];
x q[16];
x q[20];
ccx q[5],q[11],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[20],q[44];
rz(0.15) q[44];
ccx q[43],q[20],q[44];
ccx q[42],q[16],q[43];
ccx q[5],q[11],q[42];
x q[5];
x q[11];
x q[16];
x q[20];

// blue K4 (6, 12, 17, 23) -> avoid connecting to none of the four
x q[5];
x q[11];
x q[16];
x q[22];
ccx q[5],q[11],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[22],q[44];
rz(0.15) q[44];
ccx q[43],q[22],q[44];
ccx q[42],q[16],q[43];
ccx q[5],q[11],q[42];
x q[5];
x q[11];
x q[16];
x q[22];

// blue K4 (6, 12, 21, 38) -> avoid connecting to none of the four
x q[5];
x q[11];
x q[20];
x q[37];
ccx q[5],q[11],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[20],q[43];
ccx q[5],q[11],q[42];
x q[5];
x q[11];
x q[20];
x q[37];

// blue K4 (6, 12, 21, 40) -> avoid connecting to none of the four
x q[5];
x q[11];
x q[20];
x q[39];
ccx q[5],q[11],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[20],q[43];
ccx q[5],q[11],q[42];
x q[5];
x q[11];
x q[20];
x q[39];

// blue K4 (6, 12, 23, 38) -> avoid connecting to none of the four
x q[5];
x q[11];
x q[22];
x q[37];
ccx q[5],q[11],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[22],q[43];
ccx q[5],q[11],q[42];
x q[5];
x q[11];
x q[22];
x q[37];

// blue K4 (6, 12, 23, 40) -> avoid connecting to none of the four
x q[5];
x q[11];
x q[22];
x q[39];
ccx q[5],q[11],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[22],q[43];
ccx q[5],q[11],q[42];
x q[5];
x q[11];
x q[22];
x q[39];

// blue K4 (6, 14, 15, 23) -> avoid connecting to none of the four
x q[5];
x q[13];
x q[14];
x q[22];
ccx q[5],q[13],q[42];
ccx q[42],q[14],q[43];
ccx q[43],q[22],q[44];
rz(0.15) q[44];
ccx q[43],q[22],q[44];
ccx q[42],q[14],q[43];
ccx q[5],q[13],q[42];
x q[5];
x q[13];
x q[14];
x q[22];

// blue K4 (6, 14, 17, 23) -> avoid connecting to none of the four
x q[5];
x q[13];
x q[16];
x q[22];
ccx q[5],q[13],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[22],q[44];
rz(0.15) q[44];
ccx q[43],q[22],q[44];
ccx q[42],q[16],q[43];
ccx q[5],q[13],q[42];
x q[5];
x q[13];
x q[16];
x q[22];

// blue K4 (6, 14, 17, 25) -> avoid connecting to none of the four
x q[5];
x q[13];
x q[16];
x q[24];
ccx q[5],q[13],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[24],q[44];
rz(0.15) q[44];
ccx q[43],q[24],q[44];
ccx q[42],q[16],q[43];
ccx q[5],q[13],q[42];
x q[5];
x q[13];
x q[16];
x q[24];

// blue K4 (6, 14, 23, 38) -> avoid connecting to none of the four
x q[5];
x q[13];
x q[22];
x q[37];
ccx q[5],q[13],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[22],q[43];
ccx q[5],q[13],q[42];
x q[5];
x q[13];
x q[22];
x q[37];

// blue K4 (6, 14, 23, 40) -> avoid connecting to none of the four
x q[5];
x q[13];
x q[22];
x q[39];
ccx q[5],q[13],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[22],q[43];
ccx q[5],q[13],q[42];
x q[5];
x q[13];
x q[22];
x q[39];

// blue K4 (6, 14, 25, 40) -> avoid connecting to none of the four
x q[5];
x q[13];
x q[24];
x q[39];
ccx q[5],q[13],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[24],q[43];
ccx q[5],q[13],q[42];
x q[5];
x q[13];
x q[24];
x q[39];

// blue K4 (6, 15, 21, 30) -> avoid connecting to none of the four
x q[5];
x q[14];
x q[20];
x q[29];
ccx q[5],q[14],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[20],q[43];
ccx q[5],q[14],q[42];
x q[5];
x q[14];
x q[20];
x q[29];

// blue K4 (6, 15, 21, 32) -> avoid connecting to none of the four
x q[5];
x q[14];
x q[20];
x q[31];
ccx q[5],q[14],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[20],q[43];
ccx q[5],q[14],q[42];
x q[5];
x q[14];
x q[20];
x q[31];

// blue K4 (6, 15, 23, 32) -> avoid connecting to none of the four
x q[5];
x q[14];
x q[22];
x q[31];
ccx q[5],q[14],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[22],q[43];
ccx q[5],q[14],q[42];
x q[5];
x q[14];
x q[22];
x q[31];

// blue K4 (6, 15, 23, 34) -> avoid connecting to none of the four
x q[5];
x q[14];
x q[22];
x q[33];
ccx q[5],q[14],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[22],q[43];
ccx q[5],q[14],q[42];
x q[5];
x q[14];
x q[22];
x q[33];

// blue K4 (6, 15, 30, 34) -> avoid connecting to none of the four
x q[5];
x q[14];
x q[29];
x q[33];
ccx q[5],q[14],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[29],q[43];
ccx q[5],q[14],q[42];
x q[5];
x q[14];
x q[29];
x q[33];

// blue K4 (6, 15, 30, 41) -> avoid connecting to none of the four
x q[5];
x q[14];
x q[29];
x q[40];
ccx q[5],q[14],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[29],q[43];
ccx q[5],q[14],q[42];
x q[5];
x q[14];
x q[29];
x q[40];

// blue K4 (6, 15, 32, 41) -> avoid connecting to none of the four
x q[5];
x q[14];
x q[31];
x q[40];
ccx q[5],q[14],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[31],q[43];
ccx q[5],q[14],q[42];
x q[5];
x q[14];
x q[31];
x q[40];

// blue K4 (6, 17, 21, 25) -> avoid connecting to none of the four
x q[5];
x q[16];
x q[20];
x q[24];
ccx q[5],q[16],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[24],q[44];
rz(0.15) q[44];
ccx q[43],q[24],q[44];
ccx q[42],q[20],q[43];
ccx q[5],q[16],q[42];
x q[5];
x q[16];
x q[20];
x q[24];

// blue K4 (6, 17, 21, 32) -> avoid connecting to none of the four
x q[5];
x q[16];
x q[20];
x q[31];
ccx q[5],q[16],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[20],q[43];
ccx q[5],q[16],q[42];
x q[5];
x q[16];
x q[20];
x q[31];

// blue K4 (6, 17, 23, 32) -> avoid connecting to none of the four
x q[5];
x q[16];
x q[22];
x q[31];
ccx q[5],q[16],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[22],q[43];
ccx q[5],q[16],q[42];
x q[5];
x q[16];
x q[22];
x q[31];

// blue K4 (6, 17, 23, 34) -> avoid connecting to none of the four
x q[5];
x q[16];
x q[22];
x q[33];
ccx q[5],q[16],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[22],q[43];
ccx q[5],q[16],q[42];
x q[5];
x q[16];
x q[22];
x q[33];

// blue K4 (6, 17, 25, 34) -> avoid connecting to none of the four
x q[5];
x q[16];
x q[24];
x q[33];
ccx q[5],q[16],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[24],q[43];
ccx q[5],q[16],q[42];
x q[5];
x q[16];
x q[24];
x q[33];

// blue K4 (6, 17, 32, 41) -> avoid connecting to none of the four
x q[5];
x q[16];
x q[31];
x q[40];
ccx q[5],q[16],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[31],q[43];
ccx q[5],q[16],q[42];
x q[5];
x q[16];
x q[31];
x q[40];

// blue K4 (6, 21, 25, 30) -> avoid connecting to none of the four
x q[5];
x q[20];
x q[24];
x q[29];
ccx q[5],q[20],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[24],q[43];
ccx q[5],q[20],q[42];
x q[5];
x q[20];
x q[24];
x q[29];

// blue K4 (6, 21, 25, 40) -> avoid connecting to none of the four
x q[5];
x q[20];
x q[24];
x q[39];
ccx q[5],q[20],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[24],q[43];
ccx q[5],q[20],q[42];
x q[5];
x q[20];
x q[24];
x q[39];

// blue K4 (6, 21, 30, 38) -> avoid connecting to none of the four
x q[5];
x q[20];
x q[29];
x q[37];
ccx q[5],q[20],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[29],q[43];
ccx q[5],q[20],q[42];
x q[5];
x q[20];
x q[29];
x q[37];

// blue K4 (6, 21, 32, 38) -> avoid connecting to none of the four
x q[5];
x q[20];
x q[31];
x q[37];
ccx q[5],q[20],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[31],q[43];
ccx q[5],q[20],q[42];
x q[5];
x q[20];
x q[31];
x q[37];

// blue K4 (6, 21, 32, 40) -> avoid connecting to none of the four
x q[5];
x q[20];
x q[31];
x q[39];
ccx q[5],q[20],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[31],q[43];
ccx q[5],q[20],q[42];
x q[5];
x q[20];
x q[31];
x q[39];

// blue K4 (6, 23, 32, 38) -> avoid connecting to none of the four
x q[5];
x q[22];
x q[31];
x q[37];
ccx q[5],q[22],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[31],q[43];
ccx q[5],q[22],q[42];
x q[5];
x q[22];
x q[31];
x q[37];

// blue K4 (6, 23, 32, 40) -> avoid connecting to none of the four
x q[5];
x q[22];
x q[31];
x q[39];
ccx q[5],q[22],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[31],q[43];
ccx q[5],q[22],q[42];
x q[5];
x q[22];
x q[31];
x q[39];

// blue K4 (6, 23, 34, 38) -> avoid connecting to none of the four
x q[5];
x q[22];
x q[33];
x q[37];
ccx q[5],q[22],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[33],q[43];
ccx q[5],q[22],q[42];
x q[5];
x q[22];
x q[33];
x q[37];

// blue K4 (6, 23, 34, 40) -> avoid connecting to none of the four
x q[5];
x q[22];
x q[33];
x q[39];
ccx q[5],q[22],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[33],q[43];
ccx q[5],q[22],q[42];
x q[5];
x q[22];
x q[33];
x q[39];

// blue K4 (6, 25, 30, 34) -> avoid connecting to none of the four
x q[5];
x q[24];
x q[29];
x q[33];
ccx q[5],q[24],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[29],q[43];
ccx q[5],q[24],q[42];
x q[5];
x q[24];
x q[29];
x q[33];

// blue K4 (6, 25, 34, 40) -> avoid connecting to none of the four
x q[5];
x q[24];
x q[33];
x q[39];
ccx q[5],q[24],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[33],q[43];
ccx q[5],q[24],q[42];
x q[5];
x q[24];
x q[33];
x q[39];

// blue K4 (6, 30, 34, 38) -> avoid connecting to none of the four
x q[5];
x q[29];
x q[33];
x q[37];
ccx q[5],q[29],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[33],q[43];
ccx q[5],q[29],q[42];
x q[5];
x q[29];
x q[33];
x q[37];

// blue K4 (6, 30, 38, 41) -> avoid connecting to none of the four
x q[5];
x q[29];
x q[37];
x q[40];
ccx q[5],q[29],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[37],q[43];
ccx q[5],q[29],q[42];
x q[5];
x q[29];
x q[37];
x q[40];

// blue K4 (6, 32, 38, 41) -> avoid connecting to none of the four
x q[5];
x q[31];
x q[37];
x q[40];
ccx q[5],q[31],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[37],q[43];
ccx q[5],q[31],q[42];
x q[5];
x q[31];
x q[37];
x q[40];

// blue K4 (6, 32, 40, 41) -> avoid connecting to none of the four
x q[5];
x q[31];
x q[39];
x q[40];
ccx q[5],q[31],q[42];
ccx q[42],q[39],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[39],q[43];
ccx q[5],q[31],q[42];
x q[5];
x q[31];
x q[39];
x q[40];

// blue K4 (7, 8, 11, 16) -> avoid connecting to none of the four
x q[6];
x q[7];
x q[10];
x q[15];
ccx q[6],q[7],q[42];
ccx q[42],q[10],q[43];
ccx q[43],q[15],q[44];
rz(0.15) q[44];
ccx q[43],q[15],q[44];
ccx q[42],q[10],q[43];
ccx q[6],q[7],q[42];
x q[6];
x q[7];
x q[10];
x q[15];

// blue K4 (7, 8, 12, 16) -> avoid connecting to none of the four
x q[6];
x q[7];
x q[11];
x q[15];
ccx q[6],q[7],q[42];
ccx q[42],q[11],q[43];
ccx q[43],q[15],q[44];
rz(0.15) q[44];
ccx q[43],q[15],q[44];
ccx q[42],q[11],q[43];
ccx q[6],q[7],q[42];
x q[6];
x q[7];
x q[11];
x q[15];

// blue K4 (7, 8, 13, 16) -> avoid connecting to none of the four
x q[6];
x q[7];
x q[12];
x q[15];
ccx q[6],q[7],q[42];
ccx q[42],q[12],q[43];
ccx q[43],q[15],q[44];
rz(0.15) q[44];
ccx q[43],q[15],q[44];
ccx q[42],q[12],q[43];
ccx q[6],q[7],q[42];
x q[6];
x q[7];
x q[12];
x q[15];

// blue K4 (7, 8, 16, 42) -> avoid connecting to none of the four
x q[6];
x q[7];
x q[15];
x q[41];
ccx q[6],q[7],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[15],q[43];
ccx q[6],q[7],q[42];
x q[6];
x q[7];
x q[15];
x q[41];

// blue K4 (7, 10, 13, 16) -> avoid connecting to none of the four
x q[6];
x q[9];
x q[12];
x q[15];
ccx q[6],q[9],q[42];
ccx q[42],q[12],q[43];
ccx q[43],q[15],q[44];
rz(0.15) q[44];
ccx q[43],q[15],q[44];
ccx q[42],q[12],q[43];
ccx q[6],q[9],q[42];
x q[6];
x q[9];
x q[12];
x q[15];

// blue K4 (7, 10, 13, 18) -> avoid connecting to none of the four
x q[6];
x q[9];
x q[12];
x q[17];
ccx q[6],q[9],q[42];
ccx q[42],q[12],q[43];
ccx q[43],q[17],q[44];
rz(0.15) q[44];
ccx q[43],q[17],q[44];
ccx q[42],q[12],q[43];
ccx q[6],q[9],q[42];
x q[6];
x q[9];
x q[12];
x q[17];

// blue K4 (7, 10, 15, 16) -> avoid connecting to none of the four
x q[6];
x q[9];
x q[14];
x q[15];
ccx q[6],q[9],q[42];
ccx q[42],q[14],q[43];
ccx q[43],q[15],q[44];
rz(0.15) q[44];
ccx q[43],q[15],q[44];
ccx q[42],q[14],q[43];
ccx q[6],q[9],q[42];
x q[6];
x q[9];
x q[14];
x q[15];

// blue K4 (7, 10, 15, 18) -> avoid connecting to none of the four
x q[6];
x q[9];
x q[14];
x q[17];
ccx q[6],q[9],q[42];
ccx q[42],q[14],q[43];
ccx q[43],q[17],q[44];
rz(0.15) q[44];
ccx q[43],q[17],q[44];
ccx q[42],q[14],q[43];
ccx q[6],q[9],q[42];
x q[6];
x q[9];
x q[14];
x q[17];

// blue K4 (7, 10, 16, 42) -> avoid connecting to none of the four
x q[6];
x q[9];
x q[15];
x q[41];
ccx q[6],q[9],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[15],q[43];
ccx q[6],q[9],q[42];
x q[6];
x q[9];
x q[15];
x q[41];

// blue K4 (7, 10, 18, 42) -> avoid connecting to none of the four
x q[6];
x q[9];
x q[17];
x q[41];
ccx q[6],q[9],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[17],q[43];
ccx q[6],q[9],q[42];
x q[6];
x q[9];
x q[17];
x q[41];

// blue K4 (7, 11, 15, 16) -> avoid connecting to none of the four
x q[6];
x q[10];
x q[14];
x q[15];
ccx q[6],q[10],q[42];
ccx q[42],q[14],q[43];
ccx q[43],q[15],q[44];
rz(0.15) q[44];
ccx q[43],q[15],q[44];
ccx q[42],q[14],q[43];
ccx q[6],q[10],q[42];
x q[6];
x q[10];
x q[14];
x q[15];

// blue K4 (7, 11, 15, 26) -> avoid connecting to none of the four
x q[6];
x q[10];
x q[14];
x q[25];
ccx q[6],q[10],q[42];
ccx q[42],q[14],q[43];
ccx q[43],q[25],q[44];
rz(0.15) q[44];
ccx q[43],q[25],q[44];
ccx q[42],q[14],q[43];
ccx q[6],q[10],q[42];
x q[6];
x q[10];
x q[14];
x q[25];

// blue K4 (7, 11, 15, 39) -> avoid connecting to none of the four
x q[6];
x q[10];
x q[14];
x q[38];
ccx q[6],q[10],q[42];
ccx q[42],q[14],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[14],q[43];
ccx q[6],q[10],q[42];
x q[6];
x q[10];
x q[14];
x q[38];

// blue K4 (7, 11, 16, 22) -> avoid connecting to none of the four
x q[6];
x q[10];
x q[15];
x q[21];
ccx q[6],q[10],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[21],q[44];
rz(0.15) q[44];
ccx q[43],q[21],q[44];
ccx q[42],q[15],q[43];
ccx q[6],q[10],q[42];
x q[6];
x q[10];
x q[15];
x q[21];

// blue K4 (7, 11, 16, 35) -> avoid connecting to none of the four
x q[6];
x q[10];
x q[15];
x q[34];
ccx q[6],q[10],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[15],q[43];
ccx q[6],q[10],q[42];
x q[6];
x q[10];
x q[15];
x q[34];

// blue K4 (7, 11, 22, 26) -> avoid connecting to none of the four
x q[6];
x q[10];
x q[21];
x q[25];
ccx q[6],q[10],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[25],q[44];
rz(0.15) q[44];
ccx q[43],q[25],q[44];
ccx q[42],q[21],q[43];
ccx q[6],q[10],q[42];
x q[6];
x q[10];
x q[21];
x q[25];

// blue K4 (7, 11, 22, 39) -> avoid connecting to none of the four
x q[6];
x q[10];
x q[21];
x q[38];
ccx q[6],q[10],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[21],q[43];
ccx q[6],q[10],q[42];
x q[6];
x q[10];
x q[21];
x q[38];

// blue K4 (7, 11, 26, 35) -> avoid connecting to none of the four
x q[6];
x q[10];
x q[25];
x q[34];
ccx q[6],q[10],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[25],q[43];
ccx q[6],q[10],q[42];
x q[6];
x q[10];
x q[25];
x q[34];

// blue K4 (7, 11, 35, 39) -> avoid connecting to none of the four
x q[6];
x q[10];
x q[34];
x q[38];
ccx q[6],q[10],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[34],q[43];
ccx q[6],q[10],q[42];
x q[6];
x q[10];
x q[34];
x q[38];

// blue K4 (7, 12, 15, 16) -> avoid connecting to none of the four
x q[6];
x q[11];
x q[14];
x q[15];
ccx q[6],q[11],q[42];
ccx q[42],q[14],q[43];
ccx q[43],q[15],q[44];
rz(0.15) q[44];
ccx q[43],q[15],q[44];
ccx q[42],q[14],q[43];
ccx q[6],q[11],q[42];
x q[6];
x q[11];
x q[14];
x q[15];

// blue K4 (7, 12, 15, 18) -> avoid connecting to none of the four
x q[6];
x q[11];
x q[14];
x q[17];
ccx q[6],q[11],q[42];
ccx q[42],q[14],q[43];
ccx q[43],q[17],q[44];
rz(0.15) q[44];
ccx q[43],q[17],q[44];
ccx q[42],q[14],q[43];
ccx q[6],q[11],q[42];
x q[6];
x q[11];
x q[14];
x q[17];

// blue K4 (7, 12, 16, 31) -> avoid connecting to none of the four
x q[6];
x q[11];
x q[15];
x q[30];
ccx q[6],q[11],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[15],q[43];
ccx q[6],q[11],q[42];
x q[6];
x q[11];
x q[15];
x q[30];

// blue K4 (7, 13, 16, 22) -> avoid connecting to none of the four
x q[6];
x q[12];
x q[15];
x q[21];
ccx q[6],q[12],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[21],q[44];
rz(0.15) q[44];
ccx q[43],q[21],q[44];
ccx q[42],q[15],q[43];
ccx q[6],q[12],q[42];
x q[6];
x q[12];
x q[15];
x q[21];

// blue K4 (7, 13, 16, 24) -> avoid connecting to none of the four
x q[6];
x q[12];
x q[15];
x q[23];
ccx q[6],q[12],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[23],q[44];
rz(0.15) q[44];
ccx q[43],q[23],q[44];
ccx q[42],q[15],q[43];
ccx q[6],q[12],q[42];
x q[6];
x q[12];
x q[15];
x q[23];

// blue K4 (7, 13, 18, 22) -> avoid connecting to none of the four
x q[6];
x q[12];
x q[17];
x q[21];
ccx q[6],q[12],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[21],q[44];
rz(0.15) q[44];
ccx q[43],q[21],q[44];
ccx q[42],q[17],q[43];
ccx q[6],q[12],q[42];
x q[6];
x q[12];
x q[17];
x q[21];

// blue K4 (7, 13, 18, 24) -> avoid connecting to none of the four
x q[6];
x q[12];
x q[17];
x q[23];
ccx q[6],q[12],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[23],q[44];
rz(0.15) q[44];
ccx q[43],q[23],q[44];
ccx q[42],q[17],q[43];
ccx q[6],q[12],q[42];
x q[6];
x q[12];
x q[17];
x q[23];

// blue K4 (7, 13, 22, 39) -> avoid connecting to none of the four
x q[6];
x q[12];
x q[21];
x q[38];
ccx q[6],q[12],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[21],q[43];
ccx q[6],q[12],q[42];
x q[6];
x q[12];
x q[21];
x q[38];

// blue K4 (7, 13, 22, 41) -> avoid connecting to none of the four
x q[6];
x q[12];
x q[21];
x q[40];
ccx q[6],q[12],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[21],q[43];
ccx q[6],q[12],q[42];
x q[6];
x q[12];
x q[21];
x q[40];

// blue K4 (7, 13, 24, 39) -> avoid connecting to none of the four
x q[6];
x q[12];
x q[23];
x q[38];
ccx q[6],q[12],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[23],q[43];
ccx q[6],q[12],q[42];
x q[6];
x q[12];
x q[23];
x q[38];

// blue K4 (7, 13, 24, 41) -> avoid connecting to none of the four
x q[6];
x q[12];
x q[23];
x q[40];
ccx q[6],q[12],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[23],q[43];
ccx q[6],q[12],q[42];
x q[6];
x q[12];
x q[23];
x q[40];

// blue K4 (7, 15, 16, 24) -> avoid connecting to none of the four
x q[6];
x q[14];
x q[15];
x q[23];
ccx q[6],q[14],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[23],q[44];
rz(0.15) q[44];
ccx q[43],q[23],q[44];
ccx q[42],q[15],q[43];
ccx q[6],q[14],q[42];
x q[6];
x q[14];
x q[15];
x q[23];

// blue K4 (7, 15, 18, 24) -> avoid connecting to none of the four
x q[6];
x q[14];
x q[17];
x q[23];
ccx q[6],q[14],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[23],q[44];
rz(0.15) q[44];
ccx q[43],q[23],q[44];
ccx q[42],q[17],q[43];
ccx q[6],q[14],q[42];
x q[6];
x q[14];
x q[17];
x q[23];

// blue K4 (7, 15, 18, 26) -> avoid connecting to none of the four
x q[6];
x q[14];
x q[17];
x q[25];
ccx q[6],q[14],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[25],q[44];
rz(0.15) q[44];
ccx q[43],q[25],q[44];
ccx q[42],q[17],q[43];
ccx q[6],q[14],q[42];
x q[6];
x q[14];
x q[17];
x q[25];

// blue K4 (7, 15, 24, 39) -> avoid connecting to none of the four
x q[6];
x q[14];
x q[23];
x q[38];
ccx q[6],q[14],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[23],q[43];
ccx q[6],q[14],q[42];
x q[6];
x q[14];
x q[23];
x q[38];

// blue K4 (7, 15, 24, 41) -> avoid connecting to none of the four
x q[6];
x q[14];
x q[23];
x q[40];
ccx q[6],q[14],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[23],q[43];
ccx q[6],q[14],q[42];
x q[6];
x q[14];
x q[23];
x q[40];

// blue K4 (7, 15, 26, 41) -> avoid connecting to none of the four
x q[6];
x q[14];
x q[25];
x q[40];
ccx q[6],q[14],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[25],q[43];
ccx q[6],q[14],q[42];
x q[6];
x q[14];
x q[25];
x q[40];

// blue K4 (7, 16, 22, 31) -> avoid connecting to none of the four
x q[6];
x q[15];
x q[21];
x q[30];
ccx q[6],q[15],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[21],q[43];
ccx q[6],q[15],q[42];
x q[6];
x q[15];
x q[21];
x q[30];

// blue K4 (7, 16, 22, 33) -> avoid connecting to none of the four
x q[6];
x q[15];
x q[21];
x q[32];
ccx q[6],q[15],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[21],q[43];
ccx q[6],q[15],q[42];
x q[6];
x q[15];
x q[21];
x q[32];

// blue K4 (7, 16, 24, 33) -> avoid connecting to none of the four
x q[6];
x q[15];
x q[23];
x q[32];
ccx q[6],q[15],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[23],q[43];
ccx q[6],q[15],q[42];
x q[6];
x q[15];
x q[23];
x q[32];

// blue K4 (7, 16, 24, 35) -> avoid connecting to none of the four
x q[6];
x q[15];
x q[23];
x q[34];
ccx q[6],q[15],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[23],q[43];
ccx q[6],q[15],q[42];
x q[6];
x q[15];
x q[23];
x q[34];

// blue K4 (7, 16, 31, 35) -> avoid connecting to none of the four
x q[6];
x q[15];
x q[30];
x q[34];
ccx q[6],q[15],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[30],q[43];
ccx q[6],q[15],q[42];
x q[6];
x q[15];
x q[30];
x q[34];

// blue K4 (7, 16, 31, 42) -> avoid connecting to none of the four
x q[6];
x q[15];
x q[30];
x q[41];
ccx q[6],q[15],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[30],q[43];
ccx q[6],q[15],q[42];
x q[6];
x q[15];
x q[30];
x q[41];

// blue K4 (7, 16, 33, 42) -> avoid connecting to none of the four
x q[6];
x q[15];
x q[32];
x q[41];
ccx q[6],q[15],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[32],q[43];
ccx q[6],q[15],q[42];
x q[6];
x q[15];
x q[32];
x q[41];

// blue K4 (7, 18, 22, 26) -> avoid connecting to none of the four
x q[6];
x q[17];
x q[21];
x q[25];
ccx q[6],q[17],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[25],q[44];
rz(0.15) q[44];
ccx q[43],q[25],q[44];
ccx q[42],q[21],q[43];
ccx q[6],q[17],q[42];
x q[6];
x q[17];
x q[21];
x q[25];

// blue K4 (7, 18, 22, 33) -> avoid connecting to none of the four
x q[6];
x q[17];
x q[21];
x q[32];
ccx q[6],q[17],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[21],q[43];
ccx q[6],q[17],q[42];
x q[6];
x q[17];
x q[21];
x q[32];

// blue K4 (7, 18, 24, 33) -> avoid connecting to none of the four
x q[6];
x q[17];
x q[23];
x q[32];
ccx q[6],q[17],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[23],q[43];
ccx q[6],q[17],q[42];
x q[6];
x q[17];
x q[23];
x q[32];

// blue K4 (7, 18, 24, 35) -> avoid connecting to none of the four
x q[6];
x q[17];
x q[23];
x q[34];
ccx q[6],q[17],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[23],q[43];
ccx q[6],q[17],q[42];
x q[6];
x q[17];
x q[23];
x q[34];

// blue K4 (7, 18, 26, 35) -> avoid connecting to none of the four
x q[6];
x q[17];
x q[25];
x q[34];
ccx q[6],q[17],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[25],q[43];
ccx q[6],q[17],q[42];
x q[6];
x q[17];
x q[25];
x q[34];

// blue K4 (7, 18, 33, 42) -> avoid connecting to none of the four
x q[6];
x q[17];
x q[32];
x q[41];
ccx q[6],q[17],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[32],q[43];
ccx q[6],q[17],q[42];
x q[6];
x q[17];
x q[32];
x q[41];

// blue K4 (7, 22, 26, 31) -> avoid connecting to none of the four
x q[6];
x q[21];
x q[25];
x q[30];
ccx q[6],q[21],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[25],q[43];
ccx q[6],q[21],q[42];
x q[6];
x q[21];
x q[25];
x q[30];

// blue K4 (7, 22, 26, 41) -> avoid connecting to none of the four
x q[6];
x q[21];
x q[25];
x q[40];
ccx q[6],q[21],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[25],q[43];
ccx q[6],q[21],q[42];
x q[6];
x q[21];
x q[25];
x q[40];

// blue K4 (7, 22, 31, 39) -> avoid connecting to none of the four
x q[6];
x q[21];
x q[30];
x q[38];
ccx q[6],q[21],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[30],q[43];
ccx q[6],q[21],q[42];
x q[6];
x q[21];
x q[30];
x q[38];

// blue K4 (7, 22, 33, 39) -> avoid connecting to none of the four
x q[6];
x q[21];
x q[32];
x q[38];
ccx q[6],q[21],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[32],q[43];
ccx q[6],q[21],q[42];
x q[6];
x q[21];
x q[32];
x q[38];

// blue K4 (7, 22, 33, 41) -> avoid connecting to none of the four
x q[6];
x q[21];
x q[32];
x q[40];
ccx q[6],q[21],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[32],q[43];
ccx q[6],q[21],q[42];
x q[6];
x q[21];
x q[32];
x q[40];

// blue K4 (7, 24, 33, 39) -> avoid connecting to none of the four
x q[6];
x q[23];
x q[32];
x q[38];
ccx q[6],q[23],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[32],q[43];
ccx q[6],q[23],q[42];
x q[6];
x q[23];
x q[32];
x q[38];

// blue K4 (7, 24, 33, 41) -> avoid connecting to none of the four
x q[6];
x q[23];
x q[32];
x q[40];
ccx q[6],q[23],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[32],q[43];
ccx q[6],q[23],q[42];
x q[6];
x q[23];
x q[32];
x q[40];

// blue K4 (7, 24, 35, 39) -> avoid connecting to none of the four
x q[6];
x q[23];
x q[34];
x q[38];
ccx q[6],q[23],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[34],q[43];
ccx q[6],q[23],q[42];
x q[6];
x q[23];
x q[34];
x q[38];

// blue K4 (7, 24, 35, 41) -> avoid connecting to none of the four
x q[6];
x q[23];
x q[34];
x q[40];
ccx q[6],q[23],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[34],q[43];
ccx q[6],q[23],q[42];
x q[6];
x q[23];
x q[34];
x q[40];

// blue K4 (7, 26, 31, 35) -> avoid connecting to none of the four
x q[6];
x q[25];
x q[30];
x q[34];
ccx q[6],q[25],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[30],q[43];
ccx q[6],q[25],q[42];
x q[6];
x q[25];
x q[30];
x q[34];

// blue K4 (7, 26, 35, 41) -> avoid connecting to none of the four
x q[6];
x q[25];
x q[34];
x q[40];
ccx q[6],q[25],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[34],q[43];
ccx q[6],q[25],q[42];
x q[6];
x q[25];
x q[34];
x q[40];

// blue K4 (7, 31, 35, 39) -> avoid connecting to none of the four
x q[6];
x q[30];
x q[34];
x q[38];
ccx q[6],q[30],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[34],q[43];
ccx q[6],q[30],q[42];
x q[6];
x q[30];
x q[34];
x q[38];

// blue K4 (7, 31, 39, 42) -> avoid connecting to none of the four
x q[6];
x q[30];
x q[38];
x q[41];
ccx q[6],q[30],q[42];
ccx q[42],q[38],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[38],q[43];
ccx q[6],q[30],q[42];
x q[6];
x q[30];
x q[38];
x q[41];

// blue K4 (7, 33, 39, 42) -> avoid connecting to none of the four
x q[6];
x q[32];
x q[38];
x q[41];
ccx q[6],q[32],q[42];
ccx q[42],q[38],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[38],q[43];
ccx q[6],q[32],q[42];
x q[6];
x q[32];
x q[38];
x q[41];

// blue K4 (7, 33, 41, 42) -> avoid connecting to none of the four
x q[6];
x q[32];
x q[40];
x q[41];
ccx q[6],q[32],q[42];
ccx q[42],q[40],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[40],q[43];
ccx q[6],q[32],q[42];
x q[6];
x q[32];
x q[40];
x q[41];

// blue K4 (8, 11, 14, 17) -> avoid connecting to none of the four
x q[7];
x q[10];
x q[13];
x q[16];
ccx q[7],q[10],q[42];
ccx q[42],q[13],q[43];
ccx q[43],q[16],q[44];
rz(0.15) q[44];
ccx q[43],q[16],q[44];
ccx q[42],q[13],q[43];
ccx q[7],q[10],q[42];
x q[7];
x q[10];
x q[13];
x q[16];

// blue K4 (8, 11, 14, 19) -> avoid connecting to none of the four
x q[7];
x q[10];
x q[13];
x q[18];
ccx q[7],q[10],q[42];
ccx q[42],q[13],q[43];
ccx q[43],q[18],q[44];
rz(0.15) q[44];
ccx q[43],q[18],q[44];
ccx q[42],q[13],q[43];
ccx q[7],q[10],q[42];
x q[7];
x q[10];
x q[13];
x q[18];

// blue K4 (8, 11, 16, 17) -> avoid connecting to none of the four
x q[7];
x q[10];
x q[15];
x q[16];
ccx q[7],q[10],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[16],q[44];
rz(0.15) q[44];
ccx q[43],q[16],q[44];
ccx q[42],q[15],q[43];
ccx q[7],q[10],q[42];
x q[7];
x q[10];
x q[15];
x q[16];

// blue K4 (8, 11, 16, 19) -> avoid connecting to none of the four
x q[7];
x q[10];
x q[15];
x q[18];
ccx q[7],q[10],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[18],q[44];
rz(0.15) q[44];
ccx q[43],q[18],q[44];
ccx q[42],q[15],q[43];
ccx q[7],q[10],q[42];
x q[7];
x q[10];
x q[15];
x q[18];

// blue K4 (8, 11, 17, 32) -> avoid connecting to none of the four
x q[7];
x q[10];
x q[16];
x q[31];
ccx q[7],q[10],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[16],q[43];
ccx q[7],q[10],q[42];
x q[7];
x q[10];
x q[16];
x q[31];

// blue K4 (8, 12, 16, 17) -> avoid connecting to none of the four
x q[7];
x q[11];
x q[15];
x q[16];
ccx q[7],q[11],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[16],q[44];
rz(0.15) q[44];
ccx q[43],q[16],q[44];
ccx q[42],q[15],q[43];
ccx q[7],q[11],q[42];
x q[7];
x q[11];
x q[15];
x q[16];

// blue K4 (8, 12, 16, 27) -> avoid connecting to none of the four
x q[7];
x q[11];
x q[15];
x q[26];
ccx q[7],q[11],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[26],q[44];
rz(0.15) q[44];
ccx q[43],q[26],q[44];
ccx q[42],q[15],q[43];
ccx q[7],q[11],q[42];
x q[7];
x q[11];
x q[15];
x q[26];

// blue K4 (8, 12, 16, 40) -> avoid connecting to none of the four
x q[7];
x q[11];
x q[15];
x q[39];
ccx q[7],q[11],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[15],q[43];
ccx q[7],q[11],q[42];
x q[7];
x q[11];
x q[15];
x q[39];

// blue K4 (8, 12, 17, 23) -> avoid connecting to none of the four
x q[7];
x q[11];
x q[16];
x q[22];
ccx q[7],q[11],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[22],q[44];
rz(0.15) q[44];
ccx q[43],q[22],q[44];
ccx q[42],q[16],q[43];
ccx q[7],q[11],q[42];
x q[7];
x q[11];
x q[16];
x q[22];

// blue K4 (8, 12, 17, 36) -> avoid connecting to none of the four
x q[7];
x q[11];
x q[16];
x q[35];
ccx q[7],q[11],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[16],q[43];
ccx q[7],q[11],q[42];
x q[7];
x q[11];
x q[16];
x q[35];

// blue K4 (8, 12, 23, 27) -> avoid connecting to none of the four
x q[7];
x q[11];
x q[22];
x q[26];
ccx q[7],q[11],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[26],q[44];
rz(0.15) q[44];
ccx q[43],q[26],q[44];
ccx q[42],q[22],q[43];
ccx q[7],q[11],q[42];
x q[7];
x q[11];
x q[22];
x q[26];

// blue K4 (8, 12, 23, 40) -> avoid connecting to none of the four
x q[7];
x q[11];
x q[22];
x q[39];
ccx q[7],q[11],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[22],q[43];
ccx q[7],q[11],q[42];
x q[7];
x q[11];
x q[22];
x q[39];

// blue K4 (8, 12, 27, 36) -> avoid connecting to none of the four
x q[7];
x q[11];
x q[26];
x q[35];
ccx q[7],q[11],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[26],q[43];
ccx q[7],q[11],q[42];
x q[7];
x q[11];
x q[26];
x q[35];

// blue K4 (8, 12, 36, 40) -> avoid connecting to none of the four
x q[7];
x q[11];
x q[35];
x q[39];
ccx q[7],q[11],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[35],q[43];
ccx q[7],q[11],q[42];
x q[7];
x q[11];
x q[35];
x q[39];

// blue K4 (8, 13, 14, 17) -> avoid connecting to none of the four
x q[7];
x q[12];
x q[13];
x q[16];
ccx q[7],q[12],q[42];
ccx q[42],q[13],q[43];
ccx q[43],q[16],q[44];
rz(0.15) q[44];
ccx q[43],q[16],q[44];
ccx q[42],q[13],q[43];
ccx q[7],q[12],q[42];
x q[7];
x q[12];
x q[13];
x q[16];

// blue K4 (8, 13, 14, 19) -> avoid connecting to none of the four
x q[7];
x q[12];
x q[13];
x q[18];
ccx q[7],q[12],q[42];
ccx q[42],q[13],q[43];
ccx q[43],q[18],q[44];
rz(0.15) q[44];
ccx q[43],q[18],q[44];
ccx q[42],q[13],q[43];
ccx q[7],q[12],q[42];
x q[7];
x q[12];
x q[13];
x q[18];

// blue K4 (8, 13, 16, 17) -> avoid connecting to none of the four
x q[7];
x q[12];
x q[15];
x q[16];
ccx q[7],q[12],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[16],q[44];
rz(0.15) q[44];
ccx q[43],q[16],q[44];
ccx q[42],q[15],q[43];
ccx q[7],q[12],q[42];
x q[7];
x q[12];
x q[15];
x q[16];

// blue K4 (8, 13, 16, 19) -> avoid connecting to none of the four
x q[7];
x q[12];
x q[15];
x q[18];
ccx q[7],q[12],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[18],q[44];
rz(0.15) q[44];
ccx q[43],q[18],q[44];
ccx q[42],q[15],q[43];
ccx q[7],q[12],q[42];
x q[7];
x q[12];
x q[15];
x q[18];

// blue K4 (8, 13, 17, 32) -> avoid connecting to none of the four
x q[7];
x q[12];
x q[16];
x q[31];
ccx q[7],q[12],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[16],q[43];
ccx q[7],q[12],q[42];
x q[7];
x q[12];
x q[16];
x q[31];

// blue K4 (8, 14, 17, 23) -> avoid connecting to none of the four
x q[7];
x q[13];
x q[16];
x q[22];
ccx q[7],q[13],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[22],q[44];
rz(0.15) q[44];
ccx q[43],q[22],q[44];
ccx q[42],q[16],q[43];
ccx q[7],q[13],q[42];
x q[7];
x q[13];
x q[16];
x q[22];

// blue K4 (8, 14, 17, 25) -> avoid connecting to none of the four
x q[7];
x q[13];
x q[16];
x q[24];
ccx q[7],q[13],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[24],q[44];
rz(0.15) q[44];
ccx q[43],q[24],q[44];
ccx q[42],q[16],q[43];
ccx q[7],q[13],q[42];
x q[7];
x q[13];
x q[16];
x q[24];

// blue K4 (8, 14, 19, 23) -> avoid connecting to none of the four
x q[7];
x q[13];
x q[18];
x q[22];
ccx q[7],q[13],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[22],q[44];
rz(0.15) q[44];
ccx q[43],q[22],q[44];
ccx q[42],q[18],q[43];
ccx q[7],q[13],q[42];
x q[7];
x q[13];
x q[18];
x q[22];

// blue K4 (8, 14, 19, 25) -> avoid connecting to none of the four
x q[7];
x q[13];
x q[18];
x q[24];
ccx q[7],q[13],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[24],q[44];
rz(0.15) q[44];
ccx q[43],q[24],q[44];
ccx q[42],q[18],q[43];
ccx q[7],q[13],q[42];
x q[7];
x q[13];
x q[18];
x q[24];

// blue K4 (8, 14, 23, 40) -> avoid connecting to none of the four
x q[7];
x q[13];
x q[22];
x q[39];
ccx q[7],q[13],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[22],q[43];
ccx q[7],q[13],q[42];
x q[7];
x q[13];
x q[22];
x q[39];

// blue K4 (8, 14, 23, 42) -> avoid connecting to none of the four
x q[7];
x q[13];
x q[22];
x q[41];
ccx q[7],q[13],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[22],q[43];
ccx q[7],q[13],q[42];
x q[7];
x q[13];
x q[22];
x q[41];

// blue K4 (8, 14, 25, 40) -> avoid connecting to none of the four
x q[7];
x q[13];
x q[24];
x q[39];
ccx q[7],q[13],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[24],q[43];
ccx q[7],q[13],q[42];
x q[7];
x q[13];
x q[24];
x q[39];

// blue K4 (8, 14, 25, 42) -> avoid connecting to none of the four
x q[7];
x q[13];
x q[24];
x q[41];
ccx q[7],q[13],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[24],q[43];
ccx q[7],q[13],q[42];
x q[7];
x q[13];
x q[24];
x q[41];

// blue K4 (8, 16, 17, 25) -> avoid connecting to none of the four
x q[7];
x q[15];
x q[16];
x q[24];
ccx q[7],q[15],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[24],q[44];
rz(0.15) q[44];
ccx q[43],q[24],q[44];
ccx q[42],q[16],q[43];
ccx q[7],q[15],q[42];
x q[7];
x q[15];
x q[16];
x q[24];

// blue K4 (8, 16, 19, 25) -> avoid connecting to none of the four
x q[7];
x q[15];
x q[18];
x q[24];
ccx q[7],q[15],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[24],q[44];
rz(0.15) q[44];
ccx q[43],q[24],q[44];
ccx q[42],q[18],q[43];
ccx q[7],q[15],q[42];
x q[7];
x q[15];
x q[18];
x q[24];

// blue K4 (8, 16, 19, 27) -> avoid connecting to none of the four
x q[7];
x q[15];
x q[18];
x q[26];
ccx q[7],q[15],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[26],q[44];
rz(0.15) q[44];
ccx q[43],q[26],q[44];
ccx q[42],q[18],q[43];
ccx q[7],q[15],q[42];
x q[7];
x q[15];
x q[18];
x q[26];

// blue K4 (8, 16, 25, 40) -> avoid connecting to none of the four
x q[7];
x q[15];
x q[24];
x q[39];
ccx q[7],q[15],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[24],q[43];
ccx q[7],q[15],q[42];
x q[7];
x q[15];
x q[24];
x q[39];

// blue K4 (8, 16, 25, 42) -> avoid connecting to none of the four
x q[7];
x q[15];
x q[24];
x q[41];
ccx q[7],q[15],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[24],q[43];
ccx q[7],q[15],q[42];
x q[7];
x q[15];
x q[24];
x q[41];

// blue K4 (8, 16, 27, 42) -> avoid connecting to none of the four
x q[7];
x q[15];
x q[26];
x q[41];
ccx q[7],q[15],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[26],q[43];
ccx q[7],q[15],q[42];
x q[7];
x q[15];
x q[26];
x q[41];

// blue K4 (8, 17, 23, 32) -> avoid connecting to none of the four
x q[7];
x q[16];
x q[22];
x q[31];
ccx q[7],q[16],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[22],q[43];
ccx q[7],q[16],q[42];
x q[7];
x q[16];
x q[22];
x q[31];

// blue K4 (8, 17, 23, 34) -> avoid connecting to none of the four
x q[7];
x q[16];
x q[22];
x q[33];
ccx q[7],q[16],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[22],q[43];
ccx q[7],q[16],q[42];
x q[7];
x q[16];
x q[22];
x q[33];

// blue K4 (8, 17, 25, 34) -> avoid connecting to none of the four
x q[7];
x q[16];
x q[24];
x q[33];
ccx q[7],q[16],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[24],q[43];
ccx q[7],q[16],q[42];
x q[7];
x q[16];
x q[24];
x q[33];

// blue K4 (8, 17, 25, 36) -> avoid connecting to none of the four
x q[7];
x q[16];
x q[24];
x q[35];
ccx q[7],q[16],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[24],q[43];
ccx q[7],q[16],q[42];
x q[7];
x q[16];
x q[24];
x q[35];

// blue K4 (8, 17, 32, 36) -> avoid connecting to none of the four
x q[7];
x q[16];
x q[31];
x q[35];
ccx q[7],q[16],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[31],q[43];
ccx q[7],q[16],q[42];
x q[7];
x q[16];
x q[31];
x q[35];

// blue K4 (8, 19, 23, 27) -> avoid connecting to none of the four
x q[7];
x q[18];
x q[22];
x q[26];
ccx q[7],q[18],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[26],q[44];
rz(0.15) q[44];
ccx q[43],q[26],q[44];
ccx q[42],q[22],q[43];
ccx q[7],q[18],q[42];
x q[7];
x q[18];
x q[22];
x q[26];

// blue K4 (8, 19, 23, 34) -> avoid connecting to none of the four
x q[7];
x q[18];
x q[22];
x q[33];
ccx q[7],q[18],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[22],q[43];
ccx q[7],q[18],q[42];
x q[7];
x q[18];
x q[22];
x q[33];

// blue K4 (8, 19, 25, 34) -> avoid connecting to none of the four
x q[7];
x q[18];
x q[24];
x q[33];
ccx q[7],q[18],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[24],q[43];
ccx q[7],q[18],q[42];
x q[7];
x q[18];
x q[24];
x q[33];

// blue K4 (8, 19, 25, 36) -> avoid connecting to none of the four
x q[7];
x q[18];
x q[24];
x q[35];
ccx q[7],q[18],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[24],q[43];
ccx q[7],q[18],q[42];
x q[7];
x q[18];
x q[24];
x q[35];

// blue K4 (8, 19, 27, 36) -> avoid connecting to none of the four
x q[7];
x q[18];
x q[26];
x q[35];
ccx q[7],q[18],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[26],q[43];
ccx q[7],q[18],q[42];
x q[7];
x q[18];
x q[26];
x q[35];

// blue K4 (8, 23, 27, 32) -> avoid connecting to none of the four
x q[7];
x q[22];
x q[26];
x q[31];
ccx q[7],q[22],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[26],q[43];
ccx q[7],q[22],q[42];
x q[7];
x q[22];
x q[26];
x q[31];

// blue K4 (8, 23, 27, 42) -> avoid connecting to none of the four
x q[7];
x q[22];
x q[26];
x q[41];
ccx q[7],q[22],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[26],q[43];
ccx q[7],q[22],q[42];
x q[7];
x q[22];
x q[26];
x q[41];

// blue K4 (8, 23, 32, 40) -> avoid connecting to none of the four
x q[7];
x q[22];
x q[31];
x q[39];
ccx q[7],q[22],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[31],q[43];
ccx q[7],q[22],q[42];
x q[7];
x q[22];
x q[31];
x q[39];

// blue K4 (8, 23, 34, 40) -> avoid connecting to none of the four
x q[7];
x q[22];
x q[33];
x q[39];
ccx q[7],q[22],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[33],q[43];
ccx q[7],q[22],q[42];
x q[7];
x q[22];
x q[33];
x q[39];

// blue K4 (8, 23, 34, 42) -> avoid connecting to none of the four
x q[7];
x q[22];
x q[33];
x q[41];
ccx q[7],q[22],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[33],q[43];
ccx q[7],q[22],q[42];
x q[7];
x q[22];
x q[33];
x q[41];

// blue K4 (8, 25, 34, 40) -> avoid connecting to none of the four
x q[7];
x q[24];
x q[33];
x q[39];
ccx q[7],q[24],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[33],q[43];
ccx q[7],q[24],q[42];
x q[7];
x q[24];
x q[33];
x q[39];

// blue K4 (8, 25, 34, 42) -> avoid connecting to none of the four
x q[7];
x q[24];
x q[33];
x q[41];
ccx q[7],q[24],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[33],q[43];
ccx q[7],q[24],q[42];
x q[7];
x q[24];
x q[33];
x q[41];

// blue K4 (8, 25, 36, 40) -> avoid connecting to none of the four
x q[7];
x q[24];
x q[35];
x q[39];
ccx q[7],q[24],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[35],q[43];
ccx q[7],q[24],q[42];
x q[7];
x q[24];
x q[35];
x q[39];

// blue K4 (8, 25, 36, 42) -> avoid connecting to none of the four
x q[7];
x q[24];
x q[35];
x q[41];
ccx q[7],q[24],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[35],q[43];
ccx q[7],q[24],q[42];
x q[7];
x q[24];
x q[35];
x q[41];

// blue K4 (8, 27, 32, 36) -> avoid connecting to none of the four
x q[7];
x q[26];
x q[31];
x q[35];
ccx q[7],q[26],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[31],q[43];
ccx q[7],q[26],q[42];
x q[7];
x q[26];
x q[31];
x q[35];

// blue K4 (8, 27, 36, 42) -> avoid connecting to none of the four
x q[7];
x q[26];
x q[35];
x q[41];
ccx q[7],q[26],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[35],q[43];
ccx q[7],q[26],q[42];
x q[7];
x q[26];
x q[35];
x q[41];

// blue K4 (8, 32, 36, 40) -> avoid connecting to none of the four
x q[7];
x q[31];
x q[35];
x q[39];
ccx q[7],q[31],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[35],q[43];
ccx q[7],q[31],q[42];
x q[7];
x q[31];
x q[35];
x q[39];

// blue K4 (9, 12, 15, 18) -> avoid connecting to none of the four
x q[8];
x q[11];
x q[14];
x q[17];
ccx q[8],q[11],q[42];
ccx q[42],q[14],q[43];
ccx q[43],q[17],q[44];
rz(0.15) q[44];
ccx q[43],q[17],q[44];
ccx q[42],q[14],q[43];
ccx q[8],q[11],q[42];
x q[8];
x q[11];
x q[14];
x q[17];

// blue K4 (9, 12, 15, 20) -> avoid connecting to none of the four
x q[8];
x q[11];
x q[14];
x q[19];
ccx q[8],q[11],q[42];
ccx q[42],q[14],q[43];
ccx q[43],q[19],q[44];
rz(0.15) q[44];
ccx q[43],q[19],q[44];
ccx q[42],q[14],q[43];
ccx q[8],q[11],q[42];
x q[8];
x q[11];
x q[14];
x q[19];

// blue K4 (9, 12, 17, 20) -> avoid connecting to none of the four
x q[8];
x q[11];
x q[16];
x q[19];
ccx q[8],q[11],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[19],q[44];
rz(0.15) q[44];
ccx q[43],q[19],q[44];
ccx q[42],q[16],q[43];
ccx q[8],q[11],q[42];
x q[8];
x q[11];
x q[16];
x q[19];

// blue K4 (9, 13, 14, 17) -> avoid connecting to none of the four
x q[8];
x q[12];
x q[13];
x q[16];
ccx q[8],q[12],q[42];
ccx q[42],q[13],q[43];
ccx q[43],q[16],q[44];
rz(0.15) q[44];
ccx q[43],q[16],q[44];
ccx q[42],q[13],q[43];
ccx q[8],q[12],q[42];
x q[8];
x q[12];
x q[13];
x q[16];

// blue K4 (9, 13, 14, 18) -> avoid connecting to none of the four
x q[8];
x q[12];
x q[13];
x q[17];
ccx q[8],q[12],q[42];
ccx q[42],q[13],q[43];
ccx q[43],q[17],q[44];
rz(0.15) q[44];
ccx q[43],q[17],q[44];
ccx q[42],q[13],q[43];
ccx q[8],q[12],q[42];
x q[8];
x q[12];
x q[13];
x q[17];

// blue K4 (9, 13, 17, 28) -> avoid connecting to none of the four
x q[8];
x q[12];
x q[16];
x q[27];
ccx q[8],q[12],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[16],q[43];
ccx q[8],q[12],q[42];
x q[8];
x q[12];
x q[16];
x q[27];

// blue K4 (9, 13, 17, 41) -> avoid connecting to none of the four
x q[8];
x q[12];
x q[16];
x q[40];
ccx q[8],q[12],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[16],q[43];
ccx q[8],q[12],q[42];
x q[8];
x q[12];
x q[16];
x q[40];

// blue K4 (9, 13, 18, 24) -> avoid connecting to none of the four
x q[8];
x q[12];
x q[17];
x q[23];
ccx q[8],q[12],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[23],q[44];
rz(0.15) q[44];
ccx q[43],q[23],q[44];
ccx q[42],q[17],q[43];
ccx q[8],q[12],q[42];
x q[8];
x q[12];
x q[17];
x q[23];

// blue K4 (9, 13, 18, 37) -> avoid connecting to none of the four
x q[8];
x q[12];
x q[17];
x q[36];
ccx q[8],q[12],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[17],q[43];
ccx q[8],q[12],q[42];
x q[8];
x q[12];
x q[17];
x q[36];

// blue K4 (9, 13, 24, 28) -> avoid connecting to none of the four
x q[8];
x q[12];
x q[23];
x q[27];
ccx q[8],q[12],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[23],q[43];
ccx q[8],q[12],q[42];
x q[8];
x q[12];
x q[23];
x q[27];

// blue K4 (9, 13, 24, 41) -> avoid connecting to none of the four
x q[8];
x q[12];
x q[23];
x q[40];
ccx q[8],q[12],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[23],q[43];
ccx q[8],q[12],q[42];
x q[8];
x q[12];
x q[23];
x q[40];

// blue K4 (9, 13, 28, 37) -> avoid connecting to none of the four
x q[8];
x q[12];
x q[27];
x q[36];
ccx q[8],q[12],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[27],q[43];
ccx q[8],q[12],q[42];
x q[8];
x q[12];
x q[27];
x q[36];

// blue K4 (9, 13, 37, 41) -> avoid connecting to none of the four
x q[8];
x q[12];
x q[36];
x q[40];
ccx q[8],q[12],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[36],q[43];
ccx q[8],q[12],q[42];
x q[8];
x q[12];
x q[36];
x q[40];

// blue K4 (9, 14, 15, 18) -> avoid connecting to none of the four
x q[8];
x q[13];
x q[14];
x q[17];
ccx q[8],q[13],q[42];
ccx q[42],q[14],q[43];
ccx q[43],q[17],q[44];
rz(0.15) q[44];
ccx q[43],q[17],q[44];
ccx q[42],q[14],q[43];
ccx q[8],q[13],q[42];
x q[8];
x q[13];
x q[14];
x q[17];

// blue K4 (9, 14, 15, 20) -> avoid connecting to none of the four
x q[8];
x q[13];
x q[14];
x q[19];
ccx q[8],q[13],q[42];
ccx q[42],q[14],q[43];
ccx q[43],q[19],q[44];
rz(0.15) q[44];
ccx q[43],q[19],q[44];
ccx q[42],q[14],q[43];
ccx q[8],q[13],q[42];
x q[8];
x q[13];
x q[14];
x q[19];

// blue K4 (9, 14, 17, 20) -> avoid connecting to none of the four
x q[8];
x q[13];
x q[16];
x q[19];
ccx q[8],q[13],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[19],q[44];
rz(0.15) q[44];
ccx q[43],q[19],q[44];
ccx q[42],q[16],q[43];
ccx q[8],q[13],q[42];
x q[8];
x q[13];
x q[16];
x q[19];

// blue K4 (9, 14, 18, 33) -> avoid connecting to none of the four
x q[8];
x q[13];
x q[17];
x q[32];
ccx q[8],q[13],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[17],q[43];
ccx q[8],q[13],q[42];
x q[8];
x q[13];
x q[17];
x q[32];

// blue K4 (9, 15, 18, 24) -> avoid connecting to none of the four
x q[8];
x q[14];
x q[17];
x q[23];
ccx q[8],q[14],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[23],q[44];
rz(0.15) q[44];
ccx q[43],q[23],q[44];
ccx q[42],q[17],q[43];
ccx q[8],q[14],q[42];
x q[8];
x q[14];
x q[17];
x q[23];

// blue K4 (9, 15, 18, 26) -> avoid connecting to none of the four
x q[8];
x q[14];
x q[17];
x q[25];
ccx q[8],q[14],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[25],q[44];
rz(0.15) q[44];
ccx q[43],q[25],q[44];
ccx q[42],q[17],q[43];
ccx q[8],q[14],q[42];
x q[8];
x q[14];
x q[17];
x q[25];

// blue K4 (9, 15, 20, 24) -> avoid connecting to none of the four
x q[8];
x q[14];
x q[19];
x q[23];
ccx q[8],q[14],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[23],q[44];
rz(0.15) q[44];
ccx q[43],q[23],q[44];
ccx q[42],q[19],q[43];
ccx q[8],q[14],q[42];
x q[8];
x q[14];
x q[19];
x q[23];

// blue K4 (9, 15, 20, 26) -> avoid connecting to none of the four
x q[8];
x q[14];
x q[19];
x q[25];
ccx q[8],q[14],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[25],q[44];
rz(0.15) q[44];
ccx q[43],q[25],q[44];
ccx q[42],q[19],q[43];
ccx q[8],q[14],q[42];
x q[8];
x q[14];
x q[19];
x q[25];

// blue K4 (9, 15, 24, 41) -> avoid connecting to none of the four
x q[8];
x q[14];
x q[23];
x q[40];
ccx q[8],q[14],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[23],q[43];
ccx q[8],q[14],q[42];
x q[8];
x q[14];
x q[23];
x q[40];

// blue K4 (9, 15, 26, 41) -> avoid connecting to none of the four
x q[8];
x q[14];
x q[25];
x q[40];
ccx q[8],q[14],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[25],q[43];
ccx q[8],q[14],q[42];
x q[8];
x q[14];
x q[25];
x q[40];

// blue K4 (9, 17, 20, 26) -> avoid connecting to none of the four
x q[8];
x q[16];
x q[19];
x q[25];
ccx q[8],q[16],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[25],q[44];
rz(0.15) q[44];
ccx q[43],q[25],q[44];
ccx q[42],q[19],q[43];
ccx q[8],q[16],q[42];
x q[8];
x q[16];
x q[19];
x q[25];

// blue K4 (9, 17, 20, 28) -> avoid connecting to none of the four
x q[8];
x q[16];
x q[19];
x q[27];
ccx q[8],q[16],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[19],q[43];
ccx q[8],q[16],q[42];
x q[8];
x q[16];
x q[19];
x q[27];

// blue K4 (9, 17, 26, 41) -> avoid connecting to none of the four
x q[8];
x q[16];
x q[25];
x q[40];
ccx q[8],q[16],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[25],q[43];
ccx q[8],q[16],q[42];
x q[8];
x q[16];
x q[25];
x q[40];

// blue K4 (9, 18, 24, 33) -> avoid connecting to none of the four
x q[8];
x q[17];
x q[23];
x q[32];
ccx q[8],q[17],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[23],q[43];
ccx q[8],q[17],q[42];
x q[8];
x q[17];
x q[23];
x q[32];

// blue K4 (9, 18, 24, 35) -> avoid connecting to none of the four
x q[8];
x q[17];
x q[23];
x q[34];
ccx q[8],q[17],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[23],q[43];
ccx q[8],q[17],q[42];
x q[8];
x q[17];
x q[23];
x q[34];

// blue K4 (9, 18, 26, 35) -> avoid connecting to none of the four
x q[8];
x q[17];
x q[25];
x q[34];
ccx q[8],q[17],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[25],q[43];
ccx q[8],q[17],q[42];
x q[8];
x q[17];
x q[25];
x q[34];

// blue K4 (9, 18, 26, 37) -> avoid connecting to none of the four
x q[8];
x q[17];
x q[25];
x q[36];
ccx q[8],q[17],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[25],q[43];
ccx q[8],q[17],q[42];
x q[8];
x q[17];
x q[25];
x q[36];

// blue K4 (9, 18, 33, 37) -> avoid connecting to none of the four
x q[8];
x q[17];
x q[32];
x q[36];
ccx q[8],q[17],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[32],q[43];
ccx q[8],q[17],q[42];
x q[8];
x q[17];
x q[32];
x q[36];

// blue K4 (9, 20, 24, 28) -> avoid connecting to none of the four
x q[8];
x q[19];
x q[23];
x q[27];
ccx q[8],q[19],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[23],q[43];
ccx q[8],q[19],q[42];
x q[8];
x q[19];
x q[23];
x q[27];

// blue K4 (9, 20, 24, 35) -> avoid connecting to none of the four
x q[8];
x q[19];
x q[23];
x q[34];
ccx q[8],q[19],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[23],q[43];
ccx q[8],q[19],q[42];
x q[8];
x q[19];
x q[23];
x q[34];

// blue K4 (9, 20, 26, 35) -> avoid connecting to none of the four
x q[8];
x q[19];
x q[25];
x q[34];
ccx q[8],q[19],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[25],q[43];
ccx q[8],q[19],q[42];
x q[8];
x q[19];
x q[25];
x q[34];

// blue K4 (9, 20, 26, 37) -> avoid connecting to none of the four
x q[8];
x q[19];
x q[25];
x q[36];
ccx q[8],q[19],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[25],q[43];
ccx q[8],q[19],q[42];
x q[8];
x q[19];
x q[25];
x q[36];

// blue K4 (9, 20, 28, 37) -> avoid connecting to none of the four
x q[8];
x q[19];
x q[27];
x q[36];
ccx q[8],q[19],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[27],q[43];
ccx q[8],q[19],q[42];
x q[8];
x q[19];
x q[27];
x q[36];

// blue K4 (9, 24, 28, 33) -> avoid connecting to none of the four
x q[8];
x q[23];
x q[27];
x q[32];
ccx q[8],q[23],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[27],q[43];
ccx q[8],q[23],q[42];
x q[8];
x q[23];
x q[27];
x q[32];

// blue K4 (9, 24, 33, 41) -> avoid connecting to none of the four
x q[8];
x q[23];
x q[32];
x q[40];
ccx q[8],q[23],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[32],q[43];
ccx q[8],q[23],q[42];
x q[8];
x q[23];
x q[32];
x q[40];

// blue K4 (9, 24, 35, 41) -> avoid connecting to none of the four
x q[8];
x q[23];
x q[34];
x q[40];
ccx q[8],q[23],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[34],q[43];
ccx q[8],q[23],q[42];
x q[8];
x q[23];
x q[34];
x q[40];

// blue K4 (9, 26, 35, 41) -> avoid connecting to none of the four
x q[8];
x q[25];
x q[34];
x q[40];
ccx q[8],q[25],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[34],q[43];
ccx q[8],q[25],q[42];
x q[8];
x q[25];
x q[34];
x q[40];

// blue K4 (9, 26, 37, 41) -> avoid connecting to none of the four
x q[8];
x q[25];
x q[36];
x q[40];
ccx q[8],q[25],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[36],q[43];
ccx q[8],q[25],q[42];
x q[8];
x q[25];
x q[36];
x q[40];

// blue K4 (9, 28, 33, 37) -> avoid connecting to none of the four
x q[8];
x q[27];
x q[32];
x q[36];
ccx q[8],q[27],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[32],q[43];
ccx q[8],q[27],q[42];
x q[8];
x q[27];
x q[32];
x q[36];

// blue K4 (9, 33, 37, 41) -> avoid connecting to none of the four
x q[8];
x q[32];
x q[36];
x q[40];
ccx q[8],q[32],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[36],q[43];
ccx q[8],q[32],q[42];
x q[8];
x q[32];
x q[36];
x q[40];

// blue K4 (10, 13, 14, 18) -> avoid connecting to none of the four
x q[9];
x q[12];
x q[13];
x q[17];
ccx q[9],q[12],q[42];
ccx q[42],q[13],q[43];
ccx q[43],q[17],q[44];
rz(0.15) q[44];
ccx q[43],q[17],q[44];
ccx q[42],q[13],q[43];
ccx q[9],q[12],q[42];
x q[9];
x q[12];
x q[13];
x q[17];

// blue K4 (10, 13, 14, 19) -> avoid connecting to none of the four
x q[9];
x q[12];
x q[13];
x q[18];
ccx q[9],q[12],q[42];
ccx q[42],q[13],q[43];
ccx q[43],q[18],q[44];
rz(0.15) q[44];
ccx q[43],q[18],q[44];
ccx q[42],q[13],q[43];
ccx q[9],q[12],q[42];
x q[9];
x q[12];
x q[13];
x q[18];

// blue K4 (10, 13, 16, 19) -> avoid connecting to none of the four
x q[9];
x q[12];
x q[15];
x q[18];
ccx q[9],q[12],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[18],q[44];
rz(0.15) q[44];
ccx q[43],q[18],q[44];
ccx q[42],q[15],q[43];
ccx q[9],q[12],q[42];
x q[9];
x q[12];
x q[15];
x q[18];

// blue K4 (10, 13, 16, 21) -> avoid connecting to none of the four
x q[9];
x q[12];
x q[15];
x q[20];
ccx q[9],q[12],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[20],q[44];
rz(0.15) q[44];
ccx q[43],q[20],q[44];
ccx q[42],q[15],q[43];
ccx q[9],q[12],q[42];
x q[9];
x q[12];
x q[15];
x q[20];

// blue K4 (10, 13, 18, 21) -> avoid connecting to none of the four
x q[9];
x q[12];
x q[17];
x q[20];
ccx q[9],q[12],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[20],q[44];
rz(0.15) q[44];
ccx q[43],q[20],q[44];
ccx q[42],q[17],q[43];
ccx q[9],q[12],q[42];
x q[9];
x q[12];
x q[17];
x q[20];

// blue K4 (10, 14, 15, 18) -> avoid connecting to none of the four
x q[9];
x q[13];
x q[14];
x q[17];
ccx q[9],q[13],q[42];
ccx q[42],q[14],q[43];
ccx q[43],q[17],q[44];
rz(0.15) q[44];
ccx q[43],q[17],q[44];
ccx q[42],q[14],q[43];
ccx q[9],q[13],q[42];
x q[9];
x q[13];
x q[14];
x q[17];

// blue K4 (10, 14, 15, 19) -> avoid connecting to none of the four
x q[9];
x q[13];
x q[14];
x q[18];
ccx q[9],q[13],q[42];
ccx q[42],q[14],q[43];
ccx q[43],q[18],q[44];
rz(0.15) q[44];
ccx q[43],q[18],q[44];
ccx q[42],q[14],q[43];
ccx q[9],q[13],q[42];
x q[9];
x q[13];
x q[14];
x q[18];

// blue K4 (10, 14, 18, 29) -> avoid connecting to none of the four
x q[9];
x q[13];
x q[17];
x q[28];
ccx q[9],q[13],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[17],q[43];
ccx q[9],q[13],q[42];
x q[9];
x q[13];
x q[17];
x q[28];

// blue K4 (10, 14, 18, 42) -> avoid connecting to none of the four
x q[9];
x q[13];
x q[17];
x q[41];
ccx q[9],q[13],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[17],q[43];
ccx q[9],q[13],q[42];
x q[9];
x q[13];
x q[17];
x q[41];

// blue K4 (10, 14, 19, 25) -> avoid connecting to none of the four
x q[9];
x q[13];
x q[18];
x q[24];
ccx q[9],q[13],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[24],q[44];
rz(0.15) q[44];
ccx q[43],q[24],q[44];
ccx q[42],q[18],q[43];
ccx q[9],q[13],q[42];
x q[9];
x q[13];
x q[18];
x q[24];

// blue K4 (10, 14, 19, 38) -> avoid connecting to none of the four
x q[9];
x q[13];
x q[18];
x q[37];
ccx q[9],q[13],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[18],q[43];
ccx q[9],q[13],q[42];
x q[9];
x q[13];
x q[18];
x q[37];

// blue K4 (10, 14, 25, 29) -> avoid connecting to none of the four
x q[9];
x q[13];
x q[24];
x q[28];
ccx q[9],q[13],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[24],q[43];
ccx q[9],q[13],q[42];
x q[9];
x q[13];
x q[24];
x q[28];

// blue K4 (10, 14, 25, 42) -> avoid connecting to none of the four
x q[9];
x q[13];
x q[24];
x q[41];
ccx q[9],q[13],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[24],q[43];
ccx q[9],q[13],q[42];
x q[9];
x q[13];
x q[24];
x q[41];

// blue K4 (10, 14, 29, 38) -> avoid connecting to none of the four
x q[9];
x q[13];
x q[28];
x q[37];
ccx q[9],q[13],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[28],q[43];
ccx q[9],q[13],q[42];
x q[9];
x q[13];
x q[28];
x q[37];

// blue K4 (10, 14, 38, 42) -> avoid connecting to none of the four
x q[9];
x q[13];
x q[37];
x q[41];
ccx q[9],q[13],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[37],q[43];
ccx q[9],q[13],q[42];
x q[9];
x q[13];
x q[37];
x q[41];

// blue K4 (10, 15, 16, 19) -> avoid connecting to none of the four
x q[9];
x q[14];
x q[15];
x q[18];
ccx q[9],q[14],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[18],q[44];
rz(0.15) q[44];
ccx q[43],q[18],q[44];
ccx q[42],q[15],q[43];
ccx q[9],q[14],q[42];
x q[9];
x q[14];
x q[15];
x q[18];

// blue K4 (10, 15, 16, 21) -> avoid connecting to none of the four
x q[9];
x q[14];
x q[15];
x q[20];
ccx q[9],q[14],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[20],q[44];
rz(0.15) q[44];
ccx q[43],q[20],q[44];
ccx q[42],q[15],q[43];
ccx q[9],q[14],q[42];
x q[9];
x q[14];
x q[15];
x q[20];

// blue K4 (10, 15, 18, 21) -> avoid connecting to none of the four
x q[9];
x q[14];
x q[17];
x q[20];
ccx q[9],q[14],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[20],q[44];
rz(0.15) q[44];
ccx q[43],q[20],q[44];
ccx q[42],q[17],q[43];
ccx q[9],q[14],q[42];
x q[9];
x q[14];
x q[17];
x q[20];

// blue K4 (10, 15, 19, 34) -> avoid connecting to none of the four
x q[9];
x q[14];
x q[18];
x q[33];
ccx q[9],q[14],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[18],q[43];
ccx q[9],q[14],q[42];
x q[9];
x q[14];
x q[18];
x q[33];

// blue K4 (10, 16, 19, 25) -> avoid connecting to none of the four
x q[9];
x q[15];
x q[18];
x q[24];
ccx q[9],q[15],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[24],q[44];
rz(0.15) q[44];
ccx q[43],q[24],q[44];
ccx q[42],q[18],q[43];
ccx q[9],q[15],q[42];
x q[9];
x q[15];
x q[18];
x q[24];

// blue K4 (10, 16, 19, 27) -> avoid connecting to none of the four
x q[9];
x q[15];
x q[18];
x q[26];
ccx q[9],q[15],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[26],q[44];
rz(0.15) q[44];
ccx q[43],q[26],q[44];
ccx q[42],q[18],q[43];
ccx q[9],q[15],q[42];
x q[9];
x q[15];
x q[18];
x q[26];

// blue K4 (10, 16, 21, 25) -> avoid connecting to none of the four
x q[9];
x q[15];
x q[20];
x q[24];
ccx q[9],q[15],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[24],q[44];
rz(0.15) q[44];
ccx q[43],q[24],q[44];
ccx q[42],q[20],q[43];
ccx q[9],q[15],q[42];
x q[9];
x q[15];
x q[20];
x q[24];

// blue K4 (10, 16, 21, 27) -> avoid connecting to none of the four
x q[9];
x q[15];
x q[20];
x q[26];
ccx q[9],q[15],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[26],q[44];
rz(0.15) q[44];
ccx q[43],q[26],q[44];
ccx q[42],q[20],q[43];
ccx q[9],q[15],q[42];
x q[9];
x q[15];
x q[20];
x q[26];

// blue K4 (10, 16, 25, 42) -> avoid connecting to none of the four
x q[9];
x q[15];
x q[24];
x q[41];
ccx q[9],q[15],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[24],q[43];
ccx q[9],q[15],q[42];
x q[9];
x q[15];
x q[24];
x q[41];

// blue K4 (10, 16, 27, 42) -> avoid connecting to none of the four
x q[9];
x q[15];
x q[26];
x q[41];
ccx q[9],q[15],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[26],q[43];
ccx q[9],q[15],q[42];
x q[9];
x q[15];
x q[26];
x q[41];

// blue K4 (10, 18, 21, 27) -> avoid connecting to none of the four
x q[9];
x q[17];
x q[20];
x q[26];
ccx q[9],q[17],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[26],q[44];
rz(0.15) q[44];
ccx q[43],q[26],q[44];
ccx q[42],q[20],q[43];
ccx q[9],q[17],q[42];
x q[9];
x q[17];
x q[20];
x q[26];

// blue K4 (10, 18, 21, 29) -> avoid connecting to none of the four
x q[9];
x q[17];
x q[20];
x q[28];
ccx q[9],q[17],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[20],q[43];
ccx q[9],q[17],q[42];
x q[9];
x q[17];
x q[20];
x q[28];

// blue K4 (10, 18, 27, 42) -> avoid connecting to none of the four
x q[9];
x q[17];
x q[26];
x q[41];
ccx q[9],q[17],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[26],q[43];
ccx q[9],q[17],q[42];
x q[9];
x q[17];
x q[26];
x q[41];

// blue K4 (10, 19, 25, 34) -> avoid connecting to none of the four
x q[9];
x q[18];
x q[24];
x q[33];
ccx q[9],q[18],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[24],q[43];
ccx q[9],q[18],q[42];
x q[9];
x q[18];
x q[24];
x q[33];

// blue K4 (10, 19, 25, 36) -> avoid connecting to none of the four
x q[9];
x q[18];
x q[24];
x q[35];
ccx q[9],q[18],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[24],q[43];
ccx q[9],q[18],q[42];
x q[9];
x q[18];
x q[24];
x q[35];

// blue K4 (10, 19, 27, 36) -> avoid connecting to none of the four
x q[9];
x q[18];
x q[26];
x q[35];
ccx q[9],q[18],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[26],q[43];
ccx q[9],q[18],q[42];
x q[9];
x q[18];
x q[26];
x q[35];

// blue K4 (10, 19, 27, 38) -> avoid connecting to none of the four
x q[9];
x q[18];
x q[26];
x q[37];
ccx q[9],q[18],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[26],q[43];
ccx q[9],q[18],q[42];
x q[9];
x q[18];
x q[26];
x q[37];

// blue K4 (10, 19, 34, 38) -> avoid connecting to none of the four
x q[9];
x q[18];
x q[33];
x q[37];
ccx q[9],q[18],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[33],q[43];
ccx q[9],q[18],q[42];
x q[9];
x q[18];
x q[33];
x q[37];

// blue K4 (10, 21, 25, 29) -> avoid connecting to none of the four
x q[9];
x q[20];
x q[24];
x q[28];
ccx q[9],q[20],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[24],q[43];
ccx q[9],q[20],q[42];
x q[9];
x q[20];
x q[24];
x q[28];

// blue K4 (10, 21, 25, 36) -> avoid connecting to none of the four
x q[9];
x q[20];
x q[24];
x q[35];
ccx q[9],q[20],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[24],q[43];
ccx q[9],q[20],q[42];
x q[9];
x q[20];
x q[24];
x q[35];

// blue K4 (10, 21, 27, 36) -> avoid connecting to none of the four
x q[9];
x q[20];
x q[26];
x q[35];
ccx q[9],q[20],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[26],q[43];
ccx q[9],q[20],q[42];
x q[9];
x q[20];
x q[26];
x q[35];

// blue K4 (10, 21, 27, 38) -> avoid connecting to none of the four
x q[9];
x q[20];
x q[26];
x q[37];
ccx q[9],q[20],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[26],q[43];
ccx q[9],q[20],q[42];
x q[9];
x q[20];
x q[26];
x q[37];

// blue K4 (10, 21, 29, 38) -> avoid connecting to none of the four
x q[9];
x q[20];
x q[28];
x q[37];
ccx q[9],q[20],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[28],q[43];
ccx q[9],q[20],q[42];
x q[9];
x q[20];
x q[28];
x q[37];

// blue K4 (10, 25, 29, 34) -> avoid connecting to none of the four
x q[9];
x q[24];
x q[28];
x q[33];
ccx q[9],q[24],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[28],q[43];
ccx q[9],q[24],q[42];
x q[9];
x q[24];
x q[28];
x q[33];

// blue K4 (10, 25, 34, 42) -> avoid connecting to none of the four
x q[9];
x q[24];
x q[33];
x q[41];
ccx q[9],q[24],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[33],q[43];
ccx q[9],q[24],q[42];
x q[9];
x q[24];
x q[33];
x q[41];

// blue K4 (10, 25, 36, 42) -> avoid connecting to none of the four
x q[9];
x q[24];
x q[35];
x q[41];
ccx q[9],q[24],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[35],q[43];
ccx q[9],q[24],q[42];
x q[9];
x q[24];
x q[35];
x q[41];

// blue K4 (10, 27, 36, 42) -> avoid connecting to none of the four
x q[9];
x q[26];
x q[35];
x q[41];
ccx q[9],q[26],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[35],q[43];
ccx q[9],q[26],q[42];
x q[9];
x q[26];
x q[35];
x q[41];

// blue K4 (10, 27, 38, 42) -> avoid connecting to none of the four
x q[9];
x q[26];
x q[37];
x q[41];
ccx q[9],q[26],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[37],q[43];
ccx q[9],q[26],q[42];
x q[9];
x q[26];
x q[37];
x q[41];

// blue K4 (10, 29, 34, 38) -> avoid connecting to none of the four
x q[9];
x q[28];
x q[33];
x q[37];
ccx q[9],q[28],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[33],q[43];
ccx q[9],q[28],q[42];
x q[9];
x q[28];
x q[33];
x q[37];

// blue K4 (10, 34, 38, 42) -> avoid connecting to none of the four
x q[9];
x q[33];
x q[37];
x q[41];
ccx q[9],q[33],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[37],q[43];
ccx q[9],q[33],q[42];
x q[9];
x q[33];
x q[37];
x q[41];

// blue K4 (11, 14, 15, 19) -> avoid connecting to none of the four
x q[10];
x q[13];
x q[14];
x q[18];
ccx q[10],q[13],q[42];
ccx q[42],q[14],q[43];
ccx q[43],q[18],q[44];
rz(0.15) q[44];
ccx q[43],q[18],q[44];
ccx q[42],q[14],q[43];
ccx q[10],q[13],q[42];
x q[10];
x q[13];
x q[14];
x q[18];

// blue K4 (11, 14, 15, 20) -> avoid connecting to none of the four
x q[10];
x q[13];
x q[14];
x q[19];
ccx q[10],q[13],q[42];
ccx q[42],q[14],q[43];
ccx q[43],q[19],q[44];
rz(0.15) q[44];
ccx q[43],q[19],q[44];
ccx q[42],q[14],q[43];
ccx q[10],q[13],q[42];
x q[10];
x q[13];
x q[14];
x q[19];

// blue K4 (11, 14, 17, 20) -> avoid connecting to none of the four
x q[10];
x q[13];
x q[16];
x q[19];
ccx q[10],q[13],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[19],q[44];
rz(0.15) q[44];
ccx q[43],q[19],q[44];
ccx q[42],q[16],q[43];
ccx q[10],q[13],q[42];
x q[10];
x q[13];
x q[16];
x q[19];

// blue K4 (11, 14, 17, 22) -> avoid connecting to none of the four
x q[10];
x q[13];
x q[16];
x q[21];
ccx q[10],q[13],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[21],q[44];
rz(0.15) q[44];
ccx q[43],q[21],q[44];
ccx q[42],q[16],q[43];
ccx q[10],q[13],q[42];
x q[10];
x q[13];
x q[16];
x q[21];

// blue K4 (11, 14, 19, 22) -> avoid connecting to none of the four
x q[10];
x q[13];
x q[18];
x q[21];
ccx q[10],q[13],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[21],q[44];
rz(0.15) q[44];
ccx q[43],q[21],q[44];
ccx q[42],q[18],q[43];
ccx q[10],q[13],q[42];
x q[10];
x q[13];
x q[18];
x q[21];

// blue K4 (11, 15, 16, 19) -> avoid connecting to none of the four
x q[10];
x q[14];
x q[15];
x q[18];
ccx q[10],q[14],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[18],q[44];
rz(0.15) q[44];
ccx q[43],q[18],q[44];
ccx q[42],q[15],q[43];
ccx q[10],q[14],q[42];
x q[10];
x q[14];
x q[15];
x q[18];

// blue K4 (11, 15, 16, 20) -> avoid connecting to none of the four
x q[10];
x q[14];
x q[15];
x q[19];
ccx q[10],q[14],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[19],q[44];
rz(0.15) q[44];
ccx q[43],q[19],q[44];
ccx q[42],q[15],q[43];
ccx q[10],q[14],q[42];
x q[10];
x q[14];
x q[15];
x q[19];

// blue K4 (11, 15, 19, 30) -> avoid connecting to none of the four
x q[10];
x q[14];
x q[18];
x q[29];
ccx q[10],q[14],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[18],q[43];
ccx q[10],q[14],q[42];
x q[10];
x q[14];
x q[18];
x q[29];

// blue K4 (11, 15, 20, 26) -> avoid connecting to none of the four
x q[10];
x q[14];
x q[19];
x q[25];
ccx q[10],q[14],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[25],q[44];
rz(0.15) q[44];
ccx q[43],q[25],q[44];
ccx q[42],q[19],q[43];
ccx q[10],q[14],q[42];
x q[10];
x q[14];
x q[19];
x q[25];

// blue K4 (11, 15, 20, 39) -> avoid connecting to none of the four
x q[10];
x q[14];
x q[19];
x q[38];
ccx q[10],q[14],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[19],q[43];
ccx q[10],q[14],q[42];
x q[10];
x q[14];
x q[19];
x q[38];

// blue K4 (11, 15, 26, 30) -> avoid connecting to none of the four
x q[10];
x q[14];
x q[25];
x q[29];
ccx q[10],q[14],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[25],q[43];
ccx q[10],q[14],q[42];
x q[10];
x q[14];
x q[25];
x q[29];

// blue K4 (11, 15, 26, 32) -> avoid connecting to none of the four
x q[10];
x q[14];
x q[25];
x q[31];
ccx q[10],q[14],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[25],q[43];
ccx q[10],q[14],q[42];
x q[10];
x q[14];
x q[25];
x q[31];

// blue K4 (11, 15, 30, 39) -> avoid connecting to none of the four
x q[10];
x q[14];
x q[29];
x q[38];
ccx q[10],q[14],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[29],q[43];
ccx q[10],q[14],q[42];
x q[10];
x q[14];
x q[29];
x q[38];

// blue K4 (11, 16, 17, 20) -> avoid connecting to none of the four
x q[10];
x q[15];
x q[16];
x q[19];
ccx q[10],q[15],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[19],q[44];
rz(0.15) q[44];
ccx q[43],q[19],q[44];
ccx q[42],q[16],q[43];
ccx q[10],q[15],q[42];
x q[10];
x q[15];
x q[16];
x q[19];

// blue K4 (11, 16, 17, 22) -> avoid connecting to none of the four
x q[10];
x q[15];
x q[16];
x q[21];
ccx q[10],q[15],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[21],q[44];
rz(0.15) q[44];
ccx q[43],q[21],q[44];
ccx q[42],q[16],q[43];
ccx q[10],q[15],q[42];
x q[10];
x q[15];
x q[16];
x q[21];

// blue K4 (11, 16, 19, 22) -> avoid connecting to none of the four
x q[10];
x q[15];
x q[18];
x q[21];
ccx q[10],q[15],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[21],q[44];
rz(0.15) q[44];
ccx q[43],q[21],q[44];
ccx q[42],q[18],q[43];
ccx q[10],q[15],q[42];
x q[10];
x q[15];
x q[18];
x q[21];

// blue K4 (11, 16, 20, 35) -> avoid connecting to none of the four
x q[10];
x q[15];
x q[19];
x q[34];
ccx q[10],q[15],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[19],q[43];
ccx q[10],q[15],q[42];
x q[10];
x q[15];
x q[19];
x q[34];

// blue K4 (11, 17, 20, 26) -> avoid connecting to none of the four
x q[10];
x q[16];
x q[19];
x q[25];
ccx q[10],q[16],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[25],q[44];
rz(0.15) q[44];
ccx q[43],q[25],q[44];
ccx q[42],q[19],q[43];
ccx q[10],q[16],q[42];
x q[10];
x q[16];
x q[19];
x q[25];

// blue K4 (11, 17, 20, 28) -> avoid connecting to none of the four
x q[10];
x q[16];
x q[19];
x q[27];
ccx q[10],q[16],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[19],q[43];
ccx q[10],q[16],q[42];
x q[10];
x q[16];
x q[19];
x q[27];

// blue K4 (11, 17, 22, 26) -> avoid connecting to none of the four
x q[10];
x q[16];
x q[21];
x q[25];
ccx q[10],q[16],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[25],q[44];
rz(0.15) q[44];
ccx q[43],q[25],q[44];
ccx q[42],q[21],q[43];
ccx q[10],q[16],q[42];
x q[10];
x q[16];
x q[21];
x q[25];

// blue K4 (11, 17, 22, 28) -> avoid connecting to none of the four
x q[10];
x q[16];
x q[21];
x q[27];
ccx q[10],q[16],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[21],q[43];
ccx q[10],q[16],q[42];
x q[10];
x q[16];
x q[21];
x q[27];

// blue K4 (11, 17, 26, 32) -> avoid connecting to none of the four
x q[10];
x q[16];
x q[25];
x q[31];
ccx q[10],q[16],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[25],q[43];
ccx q[10],q[16],q[42];
x q[10];
x q[16];
x q[25];
x q[31];

// blue K4 (11, 17, 28, 32) -> avoid connecting to none of the four
x q[10];
x q[16];
x q[27];
x q[31];
ccx q[10],q[16],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[27],q[43];
ccx q[10],q[16],q[42];
x q[10];
x q[16];
x q[27];
x q[31];

// blue K4 (11, 19, 22, 28) -> avoid connecting to none of the four
x q[10];
x q[18];
x q[21];
x q[27];
ccx q[10],q[18],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[21],q[43];
ccx q[10],q[18],q[42];
x q[10];
x q[18];
x q[21];
x q[27];

// blue K4 (11, 19, 22, 30) -> avoid connecting to none of the four
x q[10];
x q[18];
x q[21];
x q[29];
ccx q[10],q[18],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[21],q[43];
ccx q[10],q[18],q[42];
x q[10];
x q[18];
x q[21];
x q[29];

// blue K4 (11, 20, 26, 35) -> avoid connecting to none of the four
x q[10];
x q[19];
x q[25];
x q[34];
ccx q[10],q[19],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[25],q[43];
ccx q[10],q[19],q[42];
x q[10];
x q[19];
x q[25];
x q[34];

// blue K4 (11, 20, 26, 37) -> avoid connecting to none of the four
x q[10];
x q[19];
x q[25];
x q[36];
ccx q[10],q[19],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[25],q[43];
ccx q[10],q[19],q[42];
x q[10];
x q[19];
x q[25];
x q[36];

// blue K4 (11, 20, 28, 37) -> avoid connecting to none of the four
x q[10];
x q[19];
x q[27];
x q[36];
ccx q[10],q[19],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[27],q[43];
ccx q[10],q[19],q[42];
x q[10];
x q[19];
x q[27];
x q[36];

// blue K4 (11, 20, 28, 39) -> avoid connecting to none of the four
x q[10];
x q[19];
x q[27];
x q[38];
ccx q[10],q[19],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[27],q[43];
ccx q[10],q[19],q[42];
x q[10];
x q[19];
x q[27];
x q[38];

// blue K4 (11, 20, 35, 39) -> avoid connecting to none of the four
x q[10];
x q[19];
x q[34];
x q[38];
ccx q[10],q[19],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[34],q[43];
ccx q[10],q[19],q[42];
x q[10];
x q[19];
x q[34];
x q[38];

// blue K4 (11, 22, 26, 30) -> avoid connecting to none of the four
x q[10];
x q[21];
x q[25];
x q[29];
ccx q[10],q[21],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[25],q[43];
ccx q[10],q[21],q[42];
x q[10];
x q[21];
x q[25];
x q[29];

// blue K4 (11, 22, 26, 37) -> avoid connecting to none of the four
x q[10];
x q[21];
x q[25];
x q[36];
ccx q[10],q[21],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[25],q[43];
ccx q[10],q[21],q[42];
x q[10];
x q[21];
x q[25];
x q[36];

// blue K4 (11, 22, 28, 37) -> avoid connecting to none of the four
x q[10];
x q[21];
x q[27];
x q[36];
ccx q[10],q[21],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[27],q[43];
ccx q[10],q[21],q[42];
x q[10];
x q[21];
x q[27];
x q[36];

// blue K4 (11, 22, 28, 39) -> avoid connecting to none of the four
x q[10];
x q[21];
x q[27];
x q[38];
ccx q[10],q[21],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[27],q[43];
ccx q[10],q[21],q[42];
x q[10];
x q[21];
x q[27];
x q[38];

// blue K4 (11, 22, 30, 39) -> avoid connecting to none of the four
x q[10];
x q[21];
x q[29];
x q[38];
ccx q[10],q[21],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[29],q[43];
ccx q[10],q[21],q[42];
x q[10];
x q[21];
x q[29];
x q[38];

// blue K4 (11, 26, 30, 35) -> avoid connecting to none of the four
x q[10];
x q[25];
x q[29];
x q[34];
ccx q[10],q[25],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[29],q[43];
ccx q[10],q[25],q[42];
x q[10];
x q[25];
x q[29];
x q[34];

// blue K4 (11, 26, 32, 35) -> avoid connecting to none of the four
x q[10];
x q[25];
x q[31];
x q[34];
ccx q[10],q[25],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[31],q[43];
ccx q[10],q[25],q[42];
x q[10];
x q[25];
x q[31];
x q[34];

// blue K4 (11, 26, 32, 37) -> avoid connecting to none of the four
x q[10];
x q[25];
x q[31];
x q[36];
ccx q[10],q[25],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[31],q[43];
ccx q[10],q[25],q[42];
x q[10];
x q[25];
x q[31];
x q[36];

// blue K4 (11, 28, 32, 37) -> avoid connecting to none of the four
x q[10];
x q[27];
x q[31];
x q[36];
ccx q[10],q[27],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[31],q[43];
ccx q[10],q[27],q[42];
x q[10];
x q[27];
x q[31];
x q[36];

// blue K4 (11, 30, 35, 39) -> avoid connecting to none of the four
x q[10];
x q[29];
x q[34];
x q[38];
ccx q[10],q[29],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[34],q[43];
ccx q[10],q[29],q[42];
x q[10];
x q[29];
x q[34];
x q[38];

// blue K4 (12, 15, 16, 20) -> avoid connecting to none of the four
x q[11];
x q[14];
x q[15];
x q[19];
ccx q[11],q[14],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[19],q[44];
rz(0.15) q[44];
ccx q[43],q[19],q[44];
ccx q[42],q[15],q[43];
ccx q[11],q[14],q[42];
x q[11];
x q[14];
x q[15];
x q[19];

// blue K4 (12, 15, 16, 21) -> avoid connecting to none of the four
x q[11];
x q[14];
x q[15];
x q[20];
ccx q[11],q[14],q[42];
ccx q[42],q[15],q[43];
ccx q[43],q[20],q[44];
rz(0.15) q[44];
ccx q[43],q[20],q[44];
ccx q[42],q[15],q[43];
ccx q[11],q[14],q[42];
x q[11];
x q[14];
x q[15];
x q[20];

// blue K4 (12, 15, 18, 21) -> avoid connecting to none of the four
x q[11];
x q[14];
x q[17];
x q[20];
ccx q[11],q[14],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[20],q[44];
rz(0.15) q[44];
ccx q[43],q[20],q[44];
ccx q[42],q[17],q[43];
ccx q[11],q[14],q[42];
x q[11];
x q[14];
x q[17];
x q[20];

// blue K4 (12, 15, 18, 23) -> avoid connecting to none of the four
x q[11];
x q[14];
x q[17];
x q[22];
ccx q[11],q[14],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[22],q[44];
rz(0.15) q[44];
ccx q[43],q[22],q[44];
ccx q[42],q[17],q[43];
ccx q[11],q[14],q[42];
x q[11];
x q[14];
x q[17];
x q[22];

// blue K4 (12, 15, 20, 23) -> avoid connecting to none of the four
x q[11];
x q[14];
x q[19];
x q[22];
ccx q[11],q[14],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[22],q[44];
rz(0.15) q[44];
ccx q[43],q[22],q[44];
ccx q[42],q[19],q[43];
ccx q[11],q[14],q[42];
x q[11];
x q[14];
x q[19];
x q[22];

// blue K4 (12, 16, 17, 20) -> avoid connecting to none of the four
x q[11];
x q[15];
x q[16];
x q[19];
ccx q[11],q[15],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[19],q[44];
rz(0.15) q[44];
ccx q[43],q[19],q[44];
ccx q[42],q[16],q[43];
ccx q[11],q[15],q[42];
x q[11];
x q[15];
x q[16];
x q[19];

// blue K4 (12, 16, 17, 21) -> avoid connecting to none of the four
x q[11];
x q[15];
x q[16];
x q[20];
ccx q[11],q[15],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[20],q[44];
rz(0.15) q[44];
ccx q[43],q[20],q[44];
ccx q[42],q[16],q[43];
ccx q[11],q[15],q[42];
x q[11];
x q[15];
x q[16];
x q[20];

// blue K4 (12, 16, 20, 31) -> avoid connecting to none of the four
x q[11];
x q[15];
x q[19];
x q[30];
ccx q[11],q[15],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[19],q[43];
ccx q[11],q[15],q[42];
x q[11];
x q[15];
x q[19];
x q[30];

// blue K4 (12, 16, 21, 27) -> avoid connecting to none of the four
x q[11];
x q[15];
x q[20];
x q[26];
ccx q[11],q[15],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[26],q[44];
rz(0.15) q[44];
ccx q[43],q[26],q[44];
ccx q[42],q[20],q[43];
ccx q[11],q[15],q[42];
x q[11];
x q[15];
x q[20];
x q[26];

// blue K4 (12, 16, 21, 40) -> avoid connecting to none of the four
x q[11];
x q[15];
x q[20];
x q[39];
ccx q[11],q[15],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[20],q[43];
ccx q[11],q[15],q[42];
x q[11];
x q[15];
x q[20];
x q[39];

// blue K4 (12, 16, 27, 31) -> avoid connecting to none of the four
x q[11];
x q[15];
x q[26];
x q[30];
ccx q[11],q[15],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[26],q[43];
ccx q[11],q[15],q[42];
x q[11];
x q[15];
x q[26];
x q[30];

// blue K4 (12, 16, 31, 40) -> avoid connecting to none of the four
x q[11];
x q[15];
x q[30];
x q[39];
ccx q[11],q[15],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[30],q[43];
ccx q[11],q[15],q[42];
x q[11];
x q[15];
x q[30];
x q[39];

// blue K4 (12, 17, 20, 23) -> avoid connecting to none of the four
x q[11];
x q[16];
x q[19];
x q[22];
ccx q[11],q[16],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[22],q[44];
rz(0.15) q[44];
ccx q[43],q[22],q[44];
ccx q[42],q[19],q[43];
ccx q[11],q[16],q[42];
x q[11];
x q[16];
x q[19];
x q[22];

// blue K4 (12, 17, 21, 36) -> avoid connecting to none of the four
x q[11];
x q[16];
x q[20];
x q[35];
ccx q[11],q[16],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[20],q[43];
ccx q[11],q[16],q[42];
x q[11];
x q[16];
x q[20];
x q[35];

// blue K4 (12, 18, 21, 27) -> avoid connecting to none of the four
x q[11];
x q[17];
x q[20];
x q[26];
ccx q[11],q[17],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[26],q[44];
rz(0.15) q[44];
ccx q[43],q[26],q[44];
ccx q[42],q[20],q[43];
ccx q[11],q[17],q[42];
x q[11];
x q[17];
x q[20];
x q[26];

// blue K4 (12, 18, 21, 29) -> avoid connecting to none of the four
x q[11];
x q[17];
x q[20];
x q[28];
ccx q[11],q[17],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[20],q[43];
ccx q[11],q[17],q[42];
x q[11];
x q[17];
x q[20];
x q[28];

// blue K4 (12, 18, 23, 27) -> avoid connecting to none of the four
x q[11];
x q[17];
x q[22];
x q[26];
ccx q[11],q[17],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[26],q[44];
rz(0.15) q[44];
ccx q[43],q[26],q[44];
ccx q[42],q[22],q[43];
ccx q[11],q[17],q[42];
x q[11];
x q[17];
x q[22];
x q[26];

// blue K4 (12, 18, 23, 29) -> avoid connecting to none of the four
x q[11];
x q[17];
x q[22];
x q[28];
ccx q[11],q[17],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[22],q[43];
ccx q[11],q[17],q[42];
x q[11];
x q[17];
x q[22];
x q[28];

// blue K4 (12, 20, 23, 29) -> avoid connecting to none of the four
x q[11];
x q[19];
x q[22];
x q[28];
ccx q[11],q[19],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[22],q[43];
ccx q[11],q[19],q[42];
x q[11];
x q[19];
x q[22];
x q[28];

// blue K4 (12, 20, 23, 31) -> avoid connecting to none of the four
x q[11];
x q[19];
x q[22];
x q[30];
ccx q[11],q[19],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[22],q[43];
ccx q[11],q[19],q[42];
x q[11];
x q[19];
x q[22];
x q[30];

// blue K4 (12, 21, 27, 36) -> avoid connecting to none of the four
x q[11];
x q[20];
x q[26];
x q[35];
ccx q[11],q[20],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[26],q[43];
ccx q[11],q[20],q[42];
x q[11];
x q[20];
x q[26];
x q[35];

// blue K4 (12, 21, 27, 38) -> avoid connecting to none of the four
x q[11];
x q[20];
x q[26];
x q[37];
ccx q[11],q[20],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[26],q[43];
ccx q[11],q[20],q[42];
x q[11];
x q[20];
x q[26];
x q[37];

// blue K4 (12, 21, 29, 38) -> avoid connecting to none of the four
x q[11];
x q[20];
x q[28];
x q[37];
ccx q[11],q[20],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[28],q[43];
ccx q[11],q[20],q[42];
x q[11];
x q[20];
x q[28];
x q[37];

// blue K4 (12, 21, 29, 40) -> avoid connecting to none of the four
x q[11];
x q[20];
x q[28];
x q[39];
ccx q[11],q[20],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[28],q[43];
ccx q[11],q[20],q[42];
x q[11];
x q[20];
x q[28];
x q[39];

// blue K4 (12, 21, 36, 40) -> avoid connecting to none of the four
x q[11];
x q[20];
x q[35];
x q[39];
ccx q[11],q[20],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[35],q[43];
ccx q[11],q[20],q[42];
x q[11];
x q[20];
x q[35];
x q[39];

// blue K4 (12, 23, 27, 31) -> avoid connecting to none of the four
x q[11];
x q[22];
x q[26];
x q[30];
ccx q[11],q[22],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[26],q[43];
ccx q[11],q[22],q[42];
x q[11];
x q[22];
x q[26];
x q[30];

// blue K4 (12, 23, 27, 38) -> avoid connecting to none of the four
x q[11];
x q[22];
x q[26];
x q[37];
ccx q[11],q[22],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[26],q[43];
ccx q[11],q[22],q[42];
x q[11];
x q[22];
x q[26];
x q[37];

// blue K4 (12, 23, 29, 38) -> avoid connecting to none of the four
x q[11];
x q[22];
x q[28];
x q[37];
ccx q[11],q[22],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[28],q[43];
ccx q[11],q[22],q[42];
x q[11];
x q[22];
x q[28];
x q[37];

// blue K4 (12, 23, 29, 40) -> avoid connecting to none of the four
x q[11];
x q[22];
x q[28];
x q[39];
ccx q[11],q[22],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[28],q[43];
ccx q[11],q[22],q[42];
x q[11];
x q[22];
x q[28];
x q[39];

// blue K4 (12, 23, 31, 40) -> avoid connecting to none of the four
x q[11];
x q[22];
x q[30];
x q[39];
ccx q[11],q[22],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[30],q[43];
ccx q[11],q[22],q[42];
x q[11];
x q[22];
x q[30];
x q[39];

// blue K4 (12, 27, 31, 36) -> avoid connecting to none of the four
x q[11];
x q[26];
x q[30];
x q[35];
ccx q[11],q[26],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[30],q[43];
ccx q[11],q[26],q[42];
x q[11];
x q[26];
x q[30];
x q[35];

// blue K4 (12, 31, 36, 40) -> avoid connecting to none of the four
x q[11];
x q[30];
x q[35];
x q[39];
ccx q[11],q[30],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[35],q[43];
ccx q[11],q[30],q[42];
x q[11];
x q[30];
x q[35];
x q[39];

// blue K4 (13, 14, 17, 22) -> avoid connecting to none of the four
x q[12];
x q[13];
x q[16];
x q[21];
ccx q[12],q[13],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[21],q[44];
rz(0.15) q[44];
ccx q[43],q[21],q[44];
ccx q[42],q[16],q[43];
ccx q[12],q[13],q[42];
x q[12];
x q[13];
x q[16];
x q[21];

// blue K4 (13, 14, 18, 22) -> avoid connecting to none of the four
x q[12];
x q[13];
x q[17];
x q[21];
ccx q[12],q[13],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[21],q[44];
rz(0.15) q[44];
ccx q[43],q[21],q[44];
ccx q[42],q[17],q[43];
ccx q[12],q[13],q[42];
x q[12];
x q[13];
x q[17];
x q[21];

// blue K4 (13, 14, 19, 22) -> avoid connecting to none of the four
x q[12];
x q[13];
x q[18];
x q[21];
ccx q[12],q[13],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[21],q[44];
rz(0.15) q[44];
ccx q[43],q[21],q[44];
ccx q[42],q[18],q[43];
ccx q[12],q[13],q[42];
x q[12];
x q[13];
x q[18];
x q[21];

// blue K4 (13, 16, 17, 21) -> avoid connecting to none of the four
x q[12];
x q[15];
x q[16];
x q[20];
ccx q[12],q[15],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[20],q[44];
rz(0.15) q[44];
ccx q[43],q[20],q[44];
ccx q[42],q[16],q[43];
ccx q[12],q[15],q[42];
x q[12];
x q[15];
x q[16];
x q[20];

// blue K4 (13, 16, 17, 22) -> avoid connecting to none of the four
x q[12];
x q[15];
x q[16];
x q[21];
ccx q[12],q[15],q[42];
ccx q[42],q[16],q[43];
ccx q[43],q[21],q[44];
rz(0.15) q[44];
ccx q[43],q[21],q[44];
ccx q[42],q[16],q[43];
ccx q[12],q[15],q[42];
x q[12];
x q[15];
x q[16];
x q[21];

// blue K4 (13, 16, 19, 22) -> avoid connecting to none of the four
x q[12];
x q[15];
x q[18];
x q[21];
ccx q[12],q[15],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[21],q[44];
rz(0.15) q[44];
ccx q[43],q[21],q[44];
ccx q[42],q[18],q[43];
ccx q[12],q[15],q[42];
x q[12];
x q[15];
x q[18];
x q[21];

// blue K4 (13, 16, 19, 24) -> avoid connecting to none of the four
x q[12];
x q[15];
x q[18];
x q[23];
ccx q[12],q[15],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[23],q[44];
rz(0.15) q[44];
ccx q[43],q[23],q[44];
ccx q[42],q[18],q[43];
ccx q[12],q[15],q[42];
x q[12];
x q[15];
x q[18];
x q[23];

// blue K4 (13, 16, 21, 24) -> avoid connecting to none of the four
x q[12];
x q[15];
x q[20];
x q[23];
ccx q[12],q[15],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[23],q[44];
rz(0.15) q[44];
ccx q[43],q[23],q[44];
ccx q[42],q[20],q[43];
ccx q[12],q[15],q[42];
x q[12];
x q[15];
x q[20];
x q[23];

// blue K4 (13, 17, 21, 32) -> avoid connecting to none of the four
x q[12];
x q[16];
x q[20];
x q[31];
ccx q[12],q[16],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[20],q[43];
ccx q[12],q[16],q[42];
x q[12];
x q[16];
x q[20];
x q[31];

// blue K4 (13, 17, 22, 28) -> avoid connecting to none of the four
x q[12];
x q[16];
x q[21];
x q[27];
ccx q[12],q[16],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[21],q[43];
ccx q[12],q[16],q[42];
x q[12];
x q[16];
x q[21];
x q[27];

// blue K4 (13, 17, 22, 41) -> avoid connecting to none of the four
x q[12];
x q[16];
x q[21];
x q[40];
ccx q[12],q[16],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[21],q[43];
ccx q[12],q[16],q[42];
x q[12];
x q[16];
x q[21];
x q[40];

// blue K4 (13, 17, 28, 32) -> avoid connecting to none of the four
x q[12];
x q[16];
x q[27];
x q[31];
ccx q[12],q[16],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[27],q[43];
ccx q[12],q[16],q[42];
x q[12];
x q[16];
x q[27];
x q[31];

// blue K4 (13, 17, 32, 41) -> avoid connecting to none of the four
x q[12];
x q[16];
x q[31];
x q[40];
ccx q[12],q[16],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[31],q[43];
ccx q[12],q[16],q[42];
x q[12];
x q[16];
x q[31];
x q[40];

// blue K4 (13, 18, 21, 24) -> avoid connecting to none of the four
x q[12];
x q[17];
x q[20];
x q[23];
ccx q[12],q[17],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[23],q[44];
rz(0.15) q[44];
ccx q[43],q[23],q[44];
ccx q[42],q[20],q[43];
ccx q[12],q[17],q[42];
x q[12];
x q[17];
x q[20];
x q[23];

// blue K4 (13, 18, 22, 37) -> avoid connecting to none of the four
x q[12];
x q[17];
x q[21];
x q[36];
ccx q[12],q[17],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[21],q[43];
ccx q[12],q[17],q[42];
x q[12];
x q[17];
x q[21];
x q[36];

// blue K4 (13, 19, 22, 28) -> avoid connecting to none of the four
x q[12];
x q[18];
x q[21];
x q[27];
ccx q[12],q[18],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[21],q[43];
ccx q[12],q[18],q[42];
x q[12];
x q[18];
x q[21];
x q[27];

// blue K4 (13, 19, 22, 30) -> avoid connecting to none of the four
x q[12];
x q[18];
x q[21];
x q[29];
ccx q[12],q[18],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[21],q[43];
ccx q[12],q[18],q[42];
x q[12];
x q[18];
x q[21];
x q[29];

// blue K4 (13, 19, 24, 28) -> avoid connecting to none of the four
x q[12];
x q[18];
x q[23];
x q[27];
ccx q[12],q[18],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[23],q[43];
ccx q[12],q[18],q[42];
x q[12];
x q[18];
x q[23];
x q[27];

// blue K4 (13, 19, 24, 30) -> avoid connecting to none of the four
x q[12];
x q[18];
x q[23];
x q[29];
ccx q[12],q[18],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[23],q[43];
ccx q[12],q[18],q[42];
x q[12];
x q[18];
x q[23];
x q[29];

// blue K4 (13, 21, 24, 30) -> avoid connecting to none of the four
x q[12];
x q[20];
x q[23];
x q[29];
ccx q[12],q[20],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[23],q[43];
ccx q[12],q[20],q[42];
x q[12];
x q[20];
x q[23];
x q[29];

// blue K4 (13, 21, 24, 32) -> avoid connecting to none of the four
x q[12];
x q[20];
x q[23];
x q[31];
ccx q[12],q[20],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[23],q[43];
ccx q[12],q[20],q[42];
x q[12];
x q[20];
x q[23];
x q[31];

// blue K4 (13, 22, 28, 37) -> avoid connecting to none of the four
x q[12];
x q[21];
x q[27];
x q[36];
ccx q[12],q[21],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[27],q[43];
ccx q[12],q[21],q[42];
x q[12];
x q[21];
x q[27];
x q[36];

// blue K4 (13, 22, 28, 39) -> avoid connecting to none of the four
x q[12];
x q[21];
x q[27];
x q[38];
ccx q[12],q[21],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[27],q[43];
ccx q[12],q[21],q[42];
x q[12];
x q[21];
x q[27];
x q[38];

// blue K4 (13, 22, 30, 39) -> avoid connecting to none of the four
x q[12];
x q[21];
x q[29];
x q[38];
ccx q[12],q[21],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[29],q[43];
ccx q[12],q[21],q[42];
x q[12];
x q[21];
x q[29];
x q[38];

// blue K4 (13, 22, 30, 41) -> avoid connecting to none of the four
x q[12];
x q[21];
x q[29];
x q[40];
ccx q[12],q[21],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[29],q[43];
ccx q[12],q[21],q[42];
x q[12];
x q[21];
x q[29];
x q[40];

// blue K4 (13, 22, 37, 41) -> avoid connecting to none of the four
x q[12];
x q[21];
x q[36];
x q[40];
ccx q[12],q[21],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[36],q[43];
ccx q[12],q[21],q[42];
x q[12];
x q[21];
x q[36];
x q[40];

// blue K4 (13, 24, 28, 32) -> avoid connecting to none of the four
x q[12];
x q[23];
x q[27];
x q[31];
ccx q[12],q[23],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[27],q[43];
ccx q[12],q[23],q[42];
x q[12];
x q[23];
x q[27];
x q[31];

// blue K4 (13, 24, 28, 39) -> avoid connecting to none of the four
x q[12];
x q[23];
x q[27];
x q[38];
ccx q[12],q[23],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[27],q[43];
ccx q[12],q[23],q[42];
x q[12];
x q[23];
x q[27];
x q[38];

// blue K4 (13, 24, 30, 39) -> avoid connecting to none of the four
x q[12];
x q[23];
x q[29];
x q[38];
ccx q[12],q[23],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[29],q[43];
ccx q[12],q[23],q[42];
x q[12];
x q[23];
x q[29];
x q[38];

// blue K4 (13, 24, 30, 41) -> avoid connecting to none of the four
x q[12];
x q[23];
x q[29];
x q[40];
ccx q[12],q[23],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[29],q[43];
ccx q[12],q[23],q[42];
x q[12];
x q[23];
x q[29];
x q[40];

// blue K4 (13, 24, 32, 41) -> avoid connecting to none of the four
x q[12];
x q[23];
x q[31];
x q[40];
ccx q[12],q[23],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[31],q[43];
ccx q[12],q[23],q[42];
x q[12];
x q[23];
x q[31];
x q[40];

// blue K4 (13, 28, 32, 37) -> avoid connecting to none of the four
x q[12];
x q[27];
x q[31];
x q[36];
ccx q[12],q[27],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[31],q[43];
ccx q[12],q[27],q[42];
x q[12];
x q[27];
x q[31];
x q[36];

// blue K4 (13, 32, 37, 41) -> avoid connecting to none of the four
x q[12];
x q[31];
x q[36];
x q[40];
ccx q[12],q[31],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[36],q[43];
ccx q[12],q[31],q[42];
x q[12];
x q[31];
x q[36];
x q[40];

// blue K4 (14, 15, 18, 23) -> avoid connecting to none of the four
x q[13];
x q[14];
x q[17];
x q[22];
ccx q[13],q[14],q[42];
ccx q[42],q[17],q[43];
ccx q[43],q[22],q[44];
rz(0.15) q[44];
ccx q[43],q[22],q[44];
ccx q[42],q[17],q[43];
ccx q[13],q[14],q[42];
x q[13];
x q[14];
x q[17];
x q[22];

// blue K4 (14, 15, 19, 23) -> avoid connecting to none of the four
x q[13];
x q[14];
x q[18];
x q[22];
ccx q[13],q[14],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[22],q[44];
rz(0.15) q[44];
ccx q[43],q[22],q[44];
ccx q[42],q[18],q[43];
ccx q[13],q[14],q[42];
x q[13];
x q[14];
x q[18];
x q[22];

// blue K4 (14, 15, 20, 23) -> avoid connecting to none of the four
x q[13];
x q[14];
x q[19];
x q[22];
ccx q[13],q[14],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[22],q[44];
rz(0.15) q[44];
ccx q[43],q[22],q[44];
ccx q[42],q[19],q[43];
ccx q[13],q[14],q[42];
x q[13];
x q[14];
x q[19];
x q[22];

// blue K4 (14, 17, 20, 23) -> avoid connecting to none of the four
x q[13];
x q[16];
x q[19];
x q[22];
ccx q[13],q[16],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[22],q[44];
rz(0.15) q[44];
ccx q[43],q[22],q[44];
ccx q[42],q[19],q[43];
ccx q[13],q[16],q[42];
x q[13];
x q[16];
x q[19];
x q[22];

// blue K4 (14, 17, 20, 25) -> avoid connecting to none of the four
x q[13];
x q[16];
x q[19];
x q[24];
ccx q[13],q[16],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[24],q[44];
rz(0.15) q[44];
ccx q[43],q[24],q[44];
ccx q[42],q[19],q[43];
ccx q[13],q[16],q[42];
x q[13];
x q[16];
x q[19];
x q[24];

// blue K4 (14, 17, 22, 25) -> avoid connecting to none of the four
x q[13];
x q[16];
x q[21];
x q[24];
ccx q[13],q[16],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[24],q[44];
rz(0.15) q[44];
ccx q[43],q[24],q[44];
ccx q[42],q[21],q[43];
ccx q[13],q[16],q[42];
x q[13];
x q[16];
x q[21];
x q[24];

// blue K4 (14, 18, 22, 33) -> avoid connecting to none of the four
x q[13];
x q[17];
x q[21];
x q[32];
ccx q[13],q[17],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[21],q[43];
ccx q[13],q[17],q[42];
x q[13];
x q[17];
x q[21];
x q[32];

// blue K4 (14, 18, 23, 29) -> avoid connecting to none of the four
x q[13];
x q[17];
x q[22];
x q[28];
ccx q[13],q[17],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[22],q[43];
ccx q[13],q[17],q[42];
x q[13];
x q[17];
x q[22];
x q[28];

// blue K4 (14, 18, 23, 42) -> avoid connecting to none of the four
x q[13];
x q[17];
x q[22];
x q[41];
ccx q[13],q[17],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[22],q[43];
ccx q[13],q[17],q[42];
x q[13];
x q[17];
x q[22];
x q[41];

// blue K4 (14, 18, 29, 33) -> avoid connecting to none of the four
x q[13];
x q[17];
x q[28];
x q[32];
ccx q[13],q[17],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[28],q[43];
ccx q[13],q[17],q[42];
x q[13];
x q[17];
x q[28];
x q[32];

// blue K4 (14, 18, 33, 42) -> avoid connecting to none of the four
x q[13];
x q[17];
x q[32];
x q[41];
ccx q[13],q[17],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[32],q[43];
ccx q[13],q[17],q[42];
x q[13];
x q[17];
x q[32];
x q[41];

// blue K4 (14, 19, 22, 25) -> avoid connecting to none of the four
x q[13];
x q[18];
x q[21];
x q[24];
ccx q[13],q[18],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[24],q[44];
rz(0.15) q[44];
ccx q[43],q[24],q[44];
ccx q[42],q[21],q[43];
ccx q[13],q[18],q[42];
x q[13];
x q[18];
x q[21];
x q[24];

// blue K4 (14, 19, 23, 38) -> avoid connecting to none of the four
x q[13];
x q[18];
x q[22];
x q[37];
ccx q[13],q[18],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[22],q[43];
ccx q[13],q[18],q[42];
x q[13];
x q[18];
x q[22];
x q[37];

// blue K4 (14, 20, 23, 29) -> avoid connecting to none of the four
x q[13];
x q[19];
x q[22];
x q[28];
ccx q[13],q[19],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[22],q[43];
ccx q[13],q[19],q[42];
x q[13];
x q[19];
x q[22];
x q[28];

// blue K4 (14, 20, 23, 31) -> avoid connecting to none of the four
x q[13];
x q[19];
x q[22];
x q[30];
ccx q[13],q[19],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[22],q[43];
ccx q[13],q[19],q[42];
x q[13];
x q[19];
x q[22];
x q[30];

// blue K4 (14, 20, 25, 29) -> avoid connecting to none of the four
x q[13];
x q[19];
x q[24];
x q[28];
ccx q[13],q[19],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[24],q[43];
ccx q[13],q[19],q[42];
x q[13];
x q[19];
x q[24];
x q[28];

// blue K4 (14, 20, 25, 31) -> avoid connecting to none of the four
x q[13];
x q[19];
x q[24];
x q[30];
ccx q[13],q[19],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[24],q[43];
ccx q[13],q[19],q[42];
x q[13];
x q[19];
x q[24];
x q[30];

// blue K4 (14, 22, 25, 31) -> avoid connecting to none of the four
x q[13];
x q[21];
x q[24];
x q[30];
ccx q[13],q[21],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[24],q[43];
ccx q[13],q[21],q[42];
x q[13];
x q[21];
x q[24];
x q[30];

// blue K4 (14, 22, 25, 33) -> avoid connecting to none of the four
x q[13];
x q[21];
x q[24];
x q[32];
ccx q[13],q[21],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[24],q[43];
ccx q[13],q[21],q[42];
x q[13];
x q[21];
x q[24];
x q[32];

// blue K4 (14, 23, 29, 38) -> avoid connecting to none of the four
x q[13];
x q[22];
x q[28];
x q[37];
ccx q[13],q[22],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[28],q[43];
ccx q[13],q[22],q[42];
x q[13];
x q[22];
x q[28];
x q[37];

// blue K4 (14, 23, 29, 40) -> avoid connecting to none of the four
x q[13];
x q[22];
x q[28];
x q[39];
ccx q[13],q[22],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[28],q[43];
ccx q[13],q[22],q[42];
x q[13];
x q[22];
x q[28];
x q[39];

// blue K4 (14, 23, 31, 40) -> avoid connecting to none of the four
x q[13];
x q[22];
x q[30];
x q[39];
ccx q[13],q[22],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[30],q[43];
ccx q[13],q[22],q[42];
x q[13];
x q[22];
x q[30];
x q[39];

// blue K4 (14, 23, 31, 42) -> avoid connecting to none of the four
x q[13];
x q[22];
x q[30];
x q[41];
ccx q[13],q[22],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[30],q[43];
ccx q[13],q[22],q[42];
x q[13];
x q[22];
x q[30];
x q[41];

// blue K4 (14, 23, 38, 42) -> avoid connecting to none of the four
x q[13];
x q[22];
x q[37];
x q[41];
ccx q[13],q[22],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[37],q[43];
ccx q[13],q[22],q[42];
x q[13];
x q[22];
x q[37];
x q[41];

// blue K4 (14, 25, 29, 33) -> avoid connecting to none of the four
x q[13];
x q[24];
x q[28];
x q[32];
ccx q[13],q[24],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[28],q[43];
ccx q[13],q[24],q[42];
x q[13];
x q[24];
x q[28];
x q[32];

// blue K4 (14, 25, 29, 40) -> avoid connecting to none of the four
x q[13];
x q[24];
x q[28];
x q[39];
ccx q[13],q[24],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[28],q[43];
ccx q[13],q[24],q[42];
x q[13];
x q[24];
x q[28];
x q[39];

// blue K4 (14, 25, 31, 40) -> avoid connecting to none of the four
x q[13];
x q[24];
x q[30];
x q[39];
ccx q[13],q[24],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[30],q[43];
ccx q[13],q[24],q[42];
x q[13];
x q[24];
x q[30];
x q[39];

// blue K4 (14, 25, 31, 42) -> avoid connecting to none of the four
x q[13];
x q[24];
x q[30];
x q[41];
ccx q[13],q[24],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[30],q[43];
ccx q[13],q[24],q[42];
x q[13];
x q[24];
x q[30];
x q[41];

// blue K4 (14, 25, 33, 42) -> avoid connecting to none of the four
x q[13];
x q[24];
x q[32];
x q[41];
ccx q[13],q[24],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[32],q[43];
ccx q[13],q[24],q[42];
x q[13];
x q[24];
x q[32];
x q[41];

// blue K4 (14, 29, 33, 38) -> avoid connecting to none of the four
x q[13];
x q[28];
x q[32];
x q[37];
ccx q[13],q[28],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[32],q[43];
ccx q[13],q[28],q[42];
x q[13];
x q[28];
x q[32];
x q[37];

// blue K4 (14, 33, 38, 42) -> avoid connecting to none of the four
x q[13];
x q[32];
x q[37];
x q[41];
ccx q[13],q[32],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[37],q[43];
ccx q[13],q[32],q[42];
x q[13];
x q[32];
x q[37];
x q[41];

// blue K4 (15, 16, 19, 24) -> avoid connecting to none of the four
x q[14];
x q[15];
x q[18];
x q[23];
ccx q[14],q[15],q[42];
ccx q[42],q[18],q[43];
ccx q[43],q[23],q[44];
rz(0.15) q[44];
ccx q[43],q[23],q[44];
ccx q[42],q[18],q[43];
ccx q[14],q[15],q[42];
x q[14];
x q[15];
x q[18];
x q[23];

// blue K4 (15, 16, 20, 24) -> avoid connecting to none of the four
x q[14];
x q[15];
x q[19];
x q[23];
ccx q[14],q[15],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[23],q[44];
rz(0.15) q[44];
ccx q[43],q[23],q[44];
ccx q[42],q[19],q[43];
ccx q[14],q[15],q[42];
x q[14];
x q[15];
x q[19];
x q[23];

// blue K4 (15, 16, 21, 24) -> avoid connecting to none of the four
x q[14];
x q[15];
x q[20];
x q[23];
ccx q[14],q[15],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[23],q[44];
rz(0.15) q[44];
ccx q[43],q[23],q[44];
ccx q[42],q[20],q[43];
ccx q[14],q[15],q[42];
x q[14];
x q[15];
x q[20];
x q[23];

// blue K4 (15, 18, 21, 24) -> avoid connecting to none of the four
x q[14];
x q[17];
x q[20];
x q[23];
ccx q[14],q[17],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[23],q[44];
rz(0.15) q[44];
ccx q[43],q[23],q[44];
ccx q[42],q[20],q[43];
ccx q[14],q[17],q[42];
x q[14];
x q[17];
x q[20];
x q[23];

// blue K4 (15, 18, 21, 26) -> avoid connecting to none of the four
x q[14];
x q[17];
x q[20];
x q[25];
ccx q[14],q[17],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[25],q[44];
rz(0.15) q[44];
ccx q[43],q[25],q[44];
ccx q[42],q[20],q[43];
ccx q[14],q[17],q[42];
x q[14];
x q[17];
x q[20];
x q[25];

// blue K4 (15, 18, 23, 24) -> avoid connecting to none of the four
x q[14];
x q[17];
x q[22];
x q[23];
ccx q[14],q[17],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[23],q[44];
rz(0.15) q[44];
ccx q[43],q[23],q[44];
ccx q[42],q[22],q[43];
ccx q[14],q[17],q[42];
x q[14];
x q[17];
x q[22];
x q[23];

// blue K4 (15, 18, 23, 26) -> avoid connecting to none of the four
x q[14];
x q[17];
x q[22];
x q[25];
ccx q[14],q[17],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[25],q[44];
rz(0.15) q[44];
ccx q[43],q[25],q[44];
ccx q[42],q[22],q[43];
ccx q[14],q[17],q[42];
x q[14];
x q[17];
x q[22];
x q[25];

// blue K4 (15, 19, 23, 24) -> avoid connecting to none of the four
x q[14];
x q[18];
x q[22];
x q[23];
ccx q[14],q[18],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[23],q[44];
rz(0.15) q[44];
ccx q[43],q[23],q[44];
ccx q[42],q[22],q[43];
ccx q[14],q[18],q[42];
x q[14];
x q[18];
x q[22];
x q[23];

// blue K4 (15, 19, 23, 34) -> avoid connecting to none of the four
x q[14];
x q[18];
x q[22];
x q[33];
ccx q[14],q[18],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[22],q[43];
ccx q[14],q[18],q[42];
x q[14];
x q[18];
x q[22];
x q[33];

// blue K4 (15, 19, 24, 30) -> avoid connecting to none of the four
x q[14];
x q[18];
x q[23];
x q[29];
ccx q[14],q[18],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[23],q[43];
ccx q[14],q[18],q[42];
x q[14];
x q[18];
x q[23];
x q[29];

// blue K4 (15, 19, 30, 34) -> avoid connecting to none of the four
x q[14];
x q[18];
x q[29];
x q[33];
ccx q[14],q[18],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[29],q[43];
ccx q[14],q[18],q[42];
x q[14];
x q[18];
x q[29];
x q[33];

// blue K4 (15, 20, 23, 24) -> avoid connecting to none of the four
x q[14];
x q[19];
x q[22];
x q[23];
ccx q[14],q[19],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[23],q[44];
rz(0.15) q[44];
ccx q[43],q[23],q[44];
ccx q[42],q[22],q[43];
ccx q[14],q[19],q[42];
x q[14];
x q[19];
x q[22];
x q[23];

// blue K4 (15, 20, 23, 26) -> avoid connecting to none of the four
x q[14];
x q[19];
x q[22];
x q[25];
ccx q[14],q[19],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[25],q[44];
rz(0.15) q[44];
ccx q[43],q[25],q[44];
ccx q[42],q[22],q[43];
ccx q[14],q[19],q[42];
x q[14];
x q[19];
x q[22];
x q[25];

// blue K4 (15, 20, 24, 39) -> avoid connecting to none of the four
x q[14];
x q[19];
x q[23];
x q[38];
ccx q[14],q[19],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[23],q[43];
ccx q[14],q[19],q[42];
x q[14];
x q[19];
x q[23];
x q[38];

// blue K4 (15, 21, 24, 30) -> avoid connecting to none of the four
x q[14];
x q[20];
x q[23];
x q[29];
ccx q[14],q[20],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[23],q[43];
ccx q[14],q[20],q[42];
x q[14];
x q[20];
x q[23];
x q[29];

// blue K4 (15, 21, 24, 32) -> avoid connecting to none of the four
x q[14];
x q[20];
x q[23];
x q[31];
ccx q[14],q[20],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[23],q[43];
ccx q[14],q[20],q[42];
x q[14];
x q[20];
x q[23];
x q[31];

// blue K4 (15, 21, 26, 30) -> avoid connecting to none of the four
x q[14];
x q[20];
x q[25];
x q[29];
ccx q[14],q[20],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[25],q[43];
ccx q[14],q[20],q[42];
x q[14];
x q[20];
x q[25];
x q[29];

// blue K4 (15, 21, 26, 32) -> avoid connecting to none of the four
x q[14];
x q[20];
x q[25];
x q[31];
ccx q[14],q[20],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[25],q[43];
ccx q[14],q[20],q[42];
x q[14];
x q[20];
x q[25];
x q[31];

// blue K4 (15, 23, 24, 32) -> avoid connecting to none of the four
x q[14];
x q[22];
x q[23];
x q[31];
ccx q[14],q[22],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[23],q[43];
ccx q[14],q[22],q[42];
x q[14];
x q[22];
x q[23];
x q[31];

// blue K4 (15, 23, 26, 32) -> avoid connecting to none of the four
x q[14];
x q[22];
x q[25];
x q[31];
ccx q[14],q[22],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[25],q[43];
ccx q[14],q[22],q[42];
x q[14];
x q[22];
x q[25];
x q[31];

// blue K4 (15, 23, 26, 34) -> avoid connecting to none of the four
x q[14];
x q[22];
x q[25];
x q[33];
ccx q[14],q[22],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[25],q[43];
ccx q[14],q[22],q[42];
x q[14];
x q[22];
x q[25];
x q[33];

// blue K4 (15, 24, 30, 39) -> avoid connecting to none of the four
x q[14];
x q[23];
x q[29];
x q[38];
ccx q[14],q[23],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[29],q[43];
ccx q[14],q[23],q[42];
x q[14];
x q[23];
x q[29];
x q[38];

// blue K4 (15, 24, 30, 41) -> avoid connecting to none of the four
x q[14];
x q[23];
x q[29];
x q[40];
ccx q[14],q[23],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[29],q[43];
ccx q[14],q[23],q[42];
x q[14];
x q[23];
x q[29];
x q[40];

// blue K4 (15, 24, 32, 41) -> avoid connecting to none of the four
x q[14];
x q[23];
x q[31];
x q[40];
ccx q[14],q[23],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[31],q[43];
ccx q[14],q[23],q[42];
x q[14];
x q[23];
x q[31];
x q[40];

// blue K4 (15, 26, 30, 34) -> avoid connecting to none of the four
x q[14];
x q[25];
x q[29];
x q[33];
ccx q[14],q[25],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[29],q[43];
ccx q[14],q[25],q[42];
x q[14];
x q[25];
x q[29];
x q[33];

// blue K4 (15, 26, 30, 41) -> avoid connecting to none of the four
x q[14];
x q[25];
x q[29];
x q[40];
ccx q[14],q[25],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[29],q[43];
ccx q[14],q[25],q[42];
x q[14];
x q[25];
x q[29];
x q[40];

// blue K4 (15, 26, 32, 41) -> avoid connecting to none of the four
x q[14];
x q[25];
x q[31];
x q[40];
ccx q[14],q[25],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[31],q[43];
ccx q[14],q[25],q[42];
x q[14];
x q[25];
x q[31];
x q[40];

// blue K4 (15, 30, 34, 39) -> avoid connecting to none of the four
x q[14];
x q[29];
x q[33];
x q[38];
ccx q[14],q[29],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[33],q[43];
ccx q[14],q[29],q[42];
x q[14];
x q[29];
x q[33];
x q[38];

// blue K4 (16, 17, 20, 25) -> avoid connecting to none of the four
x q[15];
x q[16];
x q[19];
x q[24];
ccx q[15],q[16],q[42];
ccx q[42],q[19],q[43];
ccx q[43],q[24],q[44];
rz(0.15) q[44];
ccx q[43],q[24],q[44];
ccx q[42],q[19],q[43];
ccx q[15],q[16],q[42];
x q[15];
x q[16];
x q[19];
x q[24];

// blue K4 (16, 17, 21, 25) -> avoid connecting to none of the four
x q[15];
x q[16];
x q[20];
x q[24];
ccx q[15],q[16],q[42];
ccx q[42],q[20],q[43];
ccx q[43],q[24],q[44];
rz(0.15) q[44];
ccx q[43],q[24],q[44];
ccx q[42],q[20],q[43];
ccx q[15],q[16],q[42];
x q[15];
x q[16];
x q[20];
x q[24];

// blue K4 (16, 17, 22, 25) -> avoid connecting to none of the four
x q[15];
x q[16];
x q[21];
x q[24];
ccx q[15],q[16],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[24],q[44];
rz(0.15) q[44];
ccx q[43],q[24],q[44];
ccx q[42],q[21],q[43];
ccx q[15],q[16],q[42];
x q[15];
x q[16];
x q[21];
x q[24];

// blue K4 (16, 19, 22, 25) -> avoid connecting to none of the four
x q[15];
x q[18];
x q[21];
x q[24];
ccx q[15],q[18],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[24],q[44];
rz(0.15) q[44];
ccx q[43],q[24],q[44];
ccx q[42],q[21],q[43];
ccx q[15],q[18],q[42];
x q[15];
x q[18];
x q[21];
x q[24];

// blue K4 (16, 19, 22, 27) -> avoid connecting to none of the four
x q[15];
x q[18];
x q[21];
x q[26];
ccx q[15],q[18],q[42];
ccx q[42],q[21],q[43];
ccx q[43],q[26],q[44];
rz(0.15) q[44];
ccx q[43],q[26],q[44];
ccx q[42],q[21],q[43];
ccx q[15],q[18],q[42];
x q[15];
x q[18];
x q[21];
x q[26];

// blue K4 (16, 19, 24, 25) -> avoid connecting to none of the four
x q[15];
x q[18];
x q[23];
x q[24];
ccx q[15],q[18],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[24],q[44];
rz(0.15) q[44];
ccx q[43],q[24],q[44];
ccx q[42],q[23],q[43];
ccx q[15],q[18],q[42];
x q[15];
x q[18];
x q[23];
x q[24];

// blue K4 (16, 19, 24, 27) -> avoid connecting to none of the four
x q[15];
x q[18];
x q[23];
x q[26];
ccx q[15],q[18],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[26],q[44];
rz(0.15) q[44];
ccx q[43],q[26],q[44];
ccx q[42],q[23],q[43];
ccx q[15],q[18],q[42];
x q[15];
x q[18];
x q[23];
x q[26];

// blue K4 (16, 20, 24, 25) -> avoid connecting to none of the four
x q[15];
x q[19];
x q[23];
x q[24];
ccx q[15],q[19],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[24],q[44];
rz(0.15) q[44];
ccx q[43],q[24],q[44];
ccx q[42],q[23],q[43];
ccx q[15],q[19],q[42];
x q[15];
x q[19];
x q[23];
x q[24];

// blue K4 (16, 20, 24, 35) -> avoid connecting to none of the four
x q[15];
x q[19];
x q[23];
x q[34];
ccx q[15],q[19],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[23],q[43];
ccx q[15],q[19],q[42];
x q[15];
x q[19];
x q[23];
x q[34];

// blue K4 (16, 20, 25, 31) -> avoid connecting to none of the four
x q[15];
x q[19];
x q[24];
x q[30];
ccx q[15],q[19],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[24],q[43];
ccx q[15],q[19],q[42];
x q[15];
x q[19];
x q[24];
x q[30];

// blue K4 (16, 20, 31, 35) -> avoid connecting to none of the four
x q[15];
x q[19];
x q[30];
x q[34];
ccx q[15],q[19],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[30],q[43];
ccx q[15],q[19],q[42];
x q[15];
x q[19];
x q[30];
x q[34];

// blue K4 (16, 21, 24, 25) -> avoid connecting to none of the four
x q[15];
x q[20];
x q[23];
x q[24];
ccx q[15],q[20],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[24],q[44];
rz(0.15) q[44];
ccx q[43],q[24],q[44];
ccx q[42],q[23],q[43];
ccx q[15],q[20],q[42];
x q[15];
x q[20];
x q[23];
x q[24];

// blue K4 (16, 21, 24, 27) -> avoid connecting to none of the four
x q[15];
x q[20];
x q[23];
x q[26];
ccx q[15],q[20],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[26],q[44];
rz(0.15) q[44];
ccx q[43],q[26],q[44];
ccx q[42],q[23],q[43];
ccx q[15],q[20],q[42];
x q[15];
x q[20];
x q[23];
x q[26];

// blue K4 (16, 21, 25, 40) -> avoid connecting to none of the four
x q[15];
x q[20];
x q[24];
x q[39];
ccx q[15],q[20],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[24],q[43];
ccx q[15],q[20],q[42];
x q[15];
x q[20];
x q[24];
x q[39];

// blue K4 (16, 22, 25, 31) -> avoid connecting to none of the four
x q[15];
x q[21];
x q[24];
x q[30];
ccx q[15],q[21],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[24],q[43];
ccx q[15],q[21],q[42];
x q[15];
x q[21];
x q[24];
x q[30];

// blue K4 (16, 22, 25, 33) -> avoid connecting to none of the four
x q[15];
x q[21];
x q[24];
x q[32];
ccx q[15],q[21],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[24],q[43];
ccx q[15],q[21],q[42];
x q[15];
x q[21];
x q[24];
x q[32];

// blue K4 (16, 22, 27, 31) -> avoid connecting to none of the four
x q[15];
x q[21];
x q[26];
x q[30];
ccx q[15],q[21],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[26],q[43];
ccx q[15],q[21],q[42];
x q[15];
x q[21];
x q[26];
x q[30];

// blue K4 (16, 22, 27, 33) -> avoid connecting to none of the four
x q[15];
x q[21];
x q[26];
x q[32];
ccx q[15],q[21],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[26],q[43];
ccx q[15],q[21],q[42];
x q[15];
x q[21];
x q[26];
x q[32];

// blue K4 (16, 24, 25, 33) -> avoid connecting to none of the four
x q[15];
x q[23];
x q[24];
x q[32];
ccx q[15],q[23],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[24],q[43];
ccx q[15],q[23],q[42];
x q[15];
x q[23];
x q[24];
x q[32];

// blue K4 (16, 24, 27, 33) -> avoid connecting to none of the four
x q[15];
x q[23];
x q[26];
x q[32];
ccx q[15],q[23],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[26],q[43];
ccx q[15],q[23],q[42];
x q[15];
x q[23];
x q[26];
x q[32];

// blue K4 (16, 24, 27, 35) -> avoid connecting to none of the four
x q[15];
x q[23];
x q[26];
x q[34];
ccx q[15],q[23],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[26],q[43];
ccx q[15],q[23],q[42];
x q[15];
x q[23];
x q[26];
x q[34];

// blue K4 (16, 25, 31, 40) -> avoid connecting to none of the four
x q[15];
x q[24];
x q[30];
x q[39];
ccx q[15],q[24],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[30],q[43];
ccx q[15],q[24],q[42];
x q[15];
x q[24];
x q[30];
x q[39];

// blue K4 (16, 25, 31, 42) -> avoid connecting to none of the four
x q[15];
x q[24];
x q[30];
x q[41];
ccx q[15],q[24],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[30],q[43];
ccx q[15],q[24],q[42];
x q[15];
x q[24];
x q[30];
x q[41];

// blue K4 (16, 25, 33, 42) -> avoid connecting to none of the four
x q[15];
x q[24];
x q[32];
x q[41];
ccx q[15],q[24],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[32],q[43];
ccx q[15],q[24],q[42];
x q[15];
x q[24];
x q[32];
x q[41];

// blue K4 (16, 27, 31, 35) -> avoid connecting to none of the four
x q[15];
x q[26];
x q[30];
x q[34];
ccx q[15],q[26],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[30],q[43];
ccx q[15],q[26],q[42];
x q[15];
x q[26];
x q[30];
x q[34];

// blue K4 (16, 27, 31, 42) -> avoid connecting to none of the four
x q[15];
x q[26];
x q[30];
x q[41];
ccx q[15],q[26],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[30],q[43];
ccx q[15],q[26],q[42];
x q[15];
x q[26];
x q[30];
x q[41];

// blue K4 (16, 27, 33, 42) -> avoid connecting to none of the four
x q[15];
x q[26];
x q[32];
x q[41];
ccx q[15],q[26],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[32],q[43];
ccx q[15],q[26],q[42];
x q[15];
x q[26];
x q[32];
x q[41];

// blue K4 (16, 31, 35, 40) -> avoid connecting to none of the four
x q[15];
x q[30];
x q[34];
x q[39];
ccx q[15],q[30],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[34],q[43];
ccx q[15],q[30],q[42];
x q[15];
x q[30];
x q[34];
x q[39];

// blue K4 (17, 20, 23, 26) -> avoid connecting to none of the four
x q[16];
x q[19];
x q[22];
x q[25];
ccx q[16],q[19],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[25],q[44];
rz(0.15) q[44];
ccx q[43],q[25],q[44];
ccx q[42],q[22],q[43];
ccx q[16],q[19],q[42];
x q[16];
x q[19];
x q[22];
x q[25];

// blue K4 (17, 20, 23, 28) -> avoid connecting to none of the four
x q[16];
x q[19];
x q[22];
x q[27];
ccx q[16],q[19],q[42];
ccx q[42],q[22],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[22],q[43];
ccx q[16],q[19],q[42];
x q[16];
x q[19];
x q[22];
x q[27];

// blue K4 (17, 20, 25, 28) -> avoid connecting to none of the four
x q[16];
x q[19];
x q[24];
x q[27];
ccx q[16],q[19],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[24],q[43];
ccx q[16],q[19],q[42];
x q[16];
x q[19];
x q[24];
x q[27];

// blue K4 (17, 21, 25, 36) -> avoid connecting to none of the four
x q[16];
x q[20];
x q[24];
x q[35];
ccx q[16],q[20],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[24],q[43];
ccx q[16],q[20],q[42];
x q[16];
x q[20];
x q[24];
x q[35];

// blue K4 (17, 21, 26, 32) -> avoid connecting to none of the four
x q[16];
x q[20];
x q[25];
x q[31];
ccx q[16],q[20],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[25],q[43];
ccx q[16],q[20],q[42];
x q[16];
x q[20];
x q[25];
x q[31];

// blue K4 (17, 21, 32, 36) -> avoid connecting to none of the four
x q[16];
x q[20];
x q[31];
x q[35];
ccx q[16],q[20],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[31],q[43];
ccx q[16],q[20],q[42];
x q[16];
x q[20];
x q[31];
x q[35];

// blue K4 (17, 22, 25, 28) -> avoid connecting to none of the four
x q[16];
x q[21];
x q[24];
x q[27];
ccx q[16],q[21],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[24],q[43];
ccx q[16],q[21],q[42];
x q[16];
x q[21];
x q[24];
x q[27];

// blue K4 (17, 22, 26, 41) -> avoid connecting to none of the four
x q[16];
x q[21];
x q[25];
x q[40];
ccx q[16],q[21],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[25],q[43];
ccx q[16],q[21],q[42];
x q[16];
x q[21];
x q[25];
x q[40];

// blue K4 (17, 23, 26, 32) -> avoid connecting to none of the four
x q[16];
x q[22];
x q[25];
x q[31];
ccx q[16],q[22],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[25],q[43];
ccx q[16],q[22],q[42];
x q[16];
x q[22];
x q[25];
x q[31];

// blue K4 (17, 23, 26, 34) -> avoid connecting to none of the four
x q[16];
x q[22];
x q[25];
x q[33];
ccx q[16],q[22],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[25],q[43];
ccx q[16],q[22],q[42];
x q[16];
x q[22];
x q[25];
x q[33];

// blue K4 (17, 23, 28, 32) -> avoid connecting to none of the four
x q[16];
x q[22];
x q[27];
x q[31];
ccx q[16],q[22],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[27],q[43];
ccx q[16],q[22],q[42];
x q[16];
x q[22];
x q[27];
x q[31];

// blue K4 (17, 23, 28, 34) -> avoid connecting to none of the four
x q[16];
x q[22];
x q[27];
x q[33];
ccx q[16],q[22],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[27],q[43];
ccx q[16],q[22],q[42];
x q[16];
x q[22];
x q[27];
x q[33];

// blue K4 (17, 25, 28, 34) -> avoid connecting to none of the four
x q[16];
x q[24];
x q[27];
x q[33];
ccx q[16],q[24],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[27],q[43];
ccx q[16],q[24],q[42];
x q[16];
x q[24];
x q[27];
x q[33];

// blue K4 (17, 25, 28, 36) -> avoid connecting to none of the four
x q[16];
x q[24];
x q[27];
x q[35];
ccx q[16],q[24],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[27],q[43];
ccx q[16],q[24],q[42];
x q[16];
x q[24];
x q[27];
x q[35];

// blue K4 (17, 26, 32, 41) -> avoid connecting to none of the four
x q[16];
x q[25];
x q[31];
x q[40];
ccx q[16],q[25],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[31],q[43];
ccx q[16],q[25],q[42];
x q[16];
x q[25];
x q[31];
x q[40];

// blue K4 (17, 28, 32, 36) -> avoid connecting to none of the four
x q[16];
x q[27];
x q[31];
x q[35];
ccx q[16],q[27],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[31],q[43];
ccx q[16],q[27],q[42];
x q[16];
x q[27];
x q[31];
x q[35];

// blue K4 (17, 32, 36, 41) -> avoid connecting to none of the four
x q[16];
x q[31];
x q[35];
x q[40];
ccx q[16],q[31],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[35],q[43];
ccx q[16],q[31],q[42];
x q[16];
x q[31];
x q[35];
x q[40];

// blue K4 (18, 21, 24, 27) -> avoid connecting to none of the four
x q[17];
x q[20];
x q[23];
x q[26];
ccx q[17],q[20],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[26],q[44];
rz(0.15) q[44];
ccx q[43],q[26],q[44];
ccx q[42],q[23],q[43];
ccx q[17],q[20],q[42];
x q[17];
x q[20];
x q[23];
x q[26];

// blue K4 (18, 21, 24, 29) -> avoid connecting to none of the four
x q[17];
x q[20];
x q[23];
x q[28];
ccx q[17],q[20],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[23],q[43];
ccx q[17],q[20],q[42];
x q[17];
x q[20];
x q[23];
x q[28];

// blue K4 (18, 21, 26, 29) -> avoid connecting to none of the four
x q[17];
x q[20];
x q[25];
x q[28];
ccx q[17],q[20],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[25],q[43];
ccx q[17],q[20],q[42];
x q[17];
x q[20];
x q[25];
x q[28];

// blue K4 (18, 22, 26, 37) -> avoid connecting to none of the four
x q[17];
x q[21];
x q[25];
x q[36];
ccx q[17],q[21],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[25],q[43];
ccx q[17],q[21],q[42];
x q[17];
x q[21];
x q[25];
x q[36];

// blue K4 (18, 22, 27, 33) -> avoid connecting to none of the four
x q[17];
x q[21];
x q[26];
x q[32];
ccx q[17],q[21],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[26],q[43];
ccx q[17],q[21],q[42];
x q[17];
x q[21];
x q[26];
x q[32];

// blue K4 (18, 22, 33, 37) -> avoid connecting to none of the four
x q[17];
x q[21];
x q[32];
x q[36];
ccx q[17],q[21],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[32],q[43];
ccx q[17],q[21],q[42];
x q[17];
x q[21];
x q[32];
x q[36];

// blue K4 (18, 23, 24, 27) -> avoid connecting to none of the four
x q[17];
x q[22];
x q[23];
x q[26];
ccx q[17],q[22],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[26],q[44];
rz(0.15) q[44];
ccx q[43],q[26],q[44];
ccx q[42],q[23],q[43];
ccx q[17],q[22],q[42];
x q[17];
x q[22];
x q[23];
x q[26];

// blue K4 (18, 23, 24, 29) -> avoid connecting to none of the four
x q[17];
x q[22];
x q[23];
x q[28];
ccx q[17],q[22],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[23],q[43];
ccx q[17],q[22],q[42];
x q[17];
x q[22];
x q[23];
x q[28];

// blue K4 (18, 23, 26, 29) -> avoid connecting to none of the four
x q[17];
x q[22];
x q[25];
x q[28];
ccx q[17],q[22],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[25],q[43];
ccx q[17],q[22],q[42];
x q[17];
x q[22];
x q[25];
x q[28];

// blue K4 (18, 23, 27, 42) -> avoid connecting to none of the four
x q[17];
x q[22];
x q[26];
x q[41];
ccx q[17],q[22],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[26],q[43];
ccx q[17],q[22],q[42];
x q[17];
x q[22];
x q[26];
x q[41];

// blue K4 (18, 24, 27, 33) -> avoid connecting to none of the four
x q[17];
x q[23];
x q[26];
x q[32];
ccx q[17],q[23],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[26],q[43];
ccx q[17],q[23],q[42];
x q[17];
x q[23];
x q[26];
x q[32];

// blue K4 (18, 24, 27, 35) -> avoid connecting to none of the four
x q[17];
x q[23];
x q[26];
x q[34];
ccx q[17],q[23],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[26],q[43];
ccx q[17],q[23],q[42];
x q[17];
x q[23];
x q[26];
x q[34];

// blue K4 (18, 24, 29, 33) -> avoid connecting to none of the four
x q[17];
x q[23];
x q[28];
x q[32];
ccx q[17],q[23],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[28],q[43];
ccx q[17],q[23],q[42];
x q[17];
x q[23];
x q[28];
x q[32];

// blue K4 (18, 24, 29, 35) -> avoid connecting to none of the four
x q[17];
x q[23];
x q[28];
x q[34];
ccx q[17],q[23],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[28],q[43];
ccx q[17],q[23],q[42];
x q[17];
x q[23];
x q[28];
x q[34];

// blue K4 (18, 26, 29, 35) -> avoid connecting to none of the four
x q[17];
x q[25];
x q[28];
x q[34];
ccx q[17],q[25],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[28],q[43];
ccx q[17],q[25],q[42];
x q[17];
x q[25];
x q[28];
x q[34];

// blue K4 (18, 26, 29, 37) -> avoid connecting to none of the four
x q[17];
x q[25];
x q[28];
x q[36];
ccx q[17],q[25],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[28],q[43];
ccx q[17],q[25],q[42];
x q[17];
x q[25];
x q[28];
x q[36];

// blue K4 (18, 27, 33, 42) -> avoid connecting to none of the four
x q[17];
x q[26];
x q[32];
x q[41];
ccx q[17],q[26],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[32],q[43];
ccx q[17],q[26],q[42];
x q[17];
x q[26];
x q[32];
x q[41];

// blue K4 (18, 29, 33, 37) -> avoid connecting to none of the four
x q[17];
x q[28];
x q[32];
x q[36];
ccx q[17],q[28],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[32],q[43];
ccx q[17],q[28],q[42];
x q[17];
x q[28];
x q[32];
x q[36];

// blue K4 (18, 33, 37, 42) -> avoid connecting to none of the four
x q[17];
x q[32];
x q[36];
x q[41];
ccx q[17],q[32],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[36],q[43];
ccx q[17],q[32],q[42];
x q[17];
x q[32];
x q[36];
x q[41];

// blue K4 (19, 22, 25, 28) -> avoid connecting to none of the four
x q[18];
x q[21];
x q[24];
x q[27];
ccx q[18],q[21],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[24],q[43];
ccx q[18],q[21],q[42];
x q[18];
x q[21];
x q[24];
x q[27];

// blue K4 (19, 22, 25, 30) -> avoid connecting to none of the four
x q[18];
x q[21];
x q[24];
x q[29];
ccx q[18],q[21],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[24],q[43];
ccx q[18],q[21],q[42];
x q[18];
x q[21];
x q[24];
x q[29];

// blue K4 (19, 22, 27, 30) -> avoid connecting to none of the four
x q[18];
x q[21];
x q[26];
x q[29];
ccx q[18],q[21],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[26],q[43];
ccx q[18],q[21],q[42];
x q[18];
x q[21];
x q[26];
x q[29];

// blue K4 (19, 23, 24, 27) -> avoid connecting to none of the four
x q[18];
x q[22];
x q[23];
x q[26];
ccx q[18],q[22],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[26],q[44];
rz(0.15) q[44];
ccx q[43],q[26],q[44];
ccx q[42],q[23],q[43];
ccx q[18],q[22],q[42];
x q[18];
x q[22];
x q[23];
x q[26];

// blue K4 (19, 23, 24, 28) -> avoid connecting to none of the four
x q[18];
x q[22];
x q[23];
x q[27];
ccx q[18],q[22],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[23],q[43];
ccx q[18],q[22],q[42];
x q[18];
x q[22];
x q[23];
x q[27];

// blue K4 (19, 23, 27, 38) -> avoid connecting to none of the four
x q[18];
x q[22];
x q[26];
x q[37];
ccx q[18],q[22],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[26],q[43];
ccx q[18],q[22],q[42];
x q[18];
x q[22];
x q[26];
x q[37];

// blue K4 (19, 23, 28, 34) -> avoid connecting to none of the four
x q[18];
x q[22];
x q[27];
x q[33];
ccx q[18],q[22],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[27],q[43];
ccx q[18],q[22],q[42];
x q[18];
x q[22];
x q[27];
x q[33];

// blue K4 (19, 23, 34, 38) -> avoid connecting to none of the four
x q[18];
x q[22];
x q[33];
x q[37];
ccx q[18],q[22],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[33],q[43];
ccx q[18],q[22],q[42];
x q[18];
x q[22];
x q[33];
x q[37];

// blue K4 (19, 24, 25, 28) -> avoid connecting to none of the four
x q[18];
x q[23];
x q[24];
x q[27];
ccx q[18],q[23],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[24],q[43];
ccx q[18],q[23],q[42];
x q[18];
x q[23];
x q[24];
x q[27];

// blue K4 (19, 24, 25, 30) -> avoid connecting to none of the four
x q[18];
x q[23];
x q[24];
x q[29];
ccx q[18],q[23],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[24],q[43];
ccx q[18],q[23],q[42];
x q[18];
x q[23];
x q[24];
x q[29];

// blue K4 (19, 24, 27, 30) -> avoid connecting to none of the four
x q[18];
x q[23];
x q[26];
x q[29];
ccx q[18],q[23],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[26],q[43];
ccx q[18],q[23],q[42];
x q[18];
x q[23];
x q[26];
x q[29];

// blue K4 (19, 25, 28, 34) -> avoid connecting to none of the four
x q[18];
x q[24];
x q[27];
x q[33];
ccx q[18],q[24],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[27],q[43];
ccx q[18],q[24],q[42];
x q[18];
x q[24];
x q[27];
x q[33];

// blue K4 (19, 25, 28, 36) -> avoid connecting to none of the four
x q[18];
x q[24];
x q[27];
x q[35];
ccx q[18],q[24],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[27],q[43];
ccx q[18],q[24],q[42];
x q[18];
x q[24];
x q[27];
x q[35];

// blue K4 (19, 25, 30, 34) -> avoid connecting to none of the four
x q[18];
x q[24];
x q[29];
x q[33];
ccx q[18],q[24],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[29],q[43];
ccx q[18],q[24],q[42];
x q[18];
x q[24];
x q[29];
x q[33];

// blue K4 (19, 25, 30, 36) -> avoid connecting to none of the four
x q[18];
x q[24];
x q[29];
x q[35];
ccx q[18],q[24],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[29],q[43];
ccx q[18],q[24],q[42];
x q[18];
x q[24];
x q[29];
x q[35];

// blue K4 (19, 27, 30, 36) -> avoid connecting to none of the four
x q[18];
x q[26];
x q[29];
x q[35];
ccx q[18],q[26],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[29],q[43];
ccx q[18],q[26],q[42];
x q[18];
x q[26];
x q[29];
x q[35];

// blue K4 (19, 27, 30, 38) -> avoid connecting to none of the four
x q[18];
x q[26];
x q[29];
x q[37];
ccx q[18],q[26],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[29],q[43];
ccx q[18],q[26],q[42];
x q[18];
x q[26];
x q[29];
x q[37];

// blue K4 (19, 30, 34, 38) -> avoid connecting to none of the four
x q[18];
x q[29];
x q[33];
x q[37];
ccx q[18],q[29],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[33],q[43];
ccx q[18],q[29],q[42];
x q[18];
x q[29];
x q[33];
x q[37];

// blue K4 (20, 23, 24, 28) -> avoid connecting to none of the four
x q[19];
x q[22];
x q[23];
x q[27];
ccx q[19],q[22],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[23],q[43];
ccx q[19],q[22],q[42];
x q[19];
x q[22];
x q[23];
x q[27];

// blue K4 (20, 23, 24, 29) -> avoid connecting to none of the four
x q[19];
x q[22];
x q[23];
x q[28];
ccx q[19],q[22],q[42];
ccx q[42],q[23],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[23],q[43];
ccx q[19],q[22],q[42];
x q[19];
x q[22];
x q[23];
x q[28];

// blue K4 (20, 23, 26, 29) -> avoid connecting to none of the four
x q[19];
x q[22];
x q[25];
x q[28];
ccx q[19],q[22],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[25],q[43];
ccx q[19],q[22],q[42];
x q[19];
x q[22];
x q[25];
x q[28];

// blue K4 (20, 23, 26, 31) -> avoid connecting to none of the four
x q[19];
x q[22];
x q[25];
x q[30];
ccx q[19],q[22],q[42];
ccx q[42],q[25],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[25],q[43];
ccx q[19],q[22],q[42];
x q[19];
x q[22];
x q[25];
x q[30];

// blue K4 (20, 23, 28, 31) -> avoid connecting to none of the four
x q[19];
x q[22];
x q[27];
x q[30];
ccx q[19],q[22],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[27],q[43];
ccx q[19],q[22],q[42];
x q[19];
x q[22];
x q[27];
x q[30];

// blue K4 (20, 24, 25, 28) -> avoid connecting to none of the four
x q[19];
x q[23];
x q[24];
x q[27];
ccx q[19],q[23],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[27],q[44];
rz(0.15) q[44];
ccx q[43],q[27],q[44];
ccx q[42],q[24],q[43];
ccx q[19],q[23],q[42];
x q[19];
x q[23];
x q[24];
x q[27];

// blue K4 (20, 24, 25, 29) -> avoid connecting to none of the four
x q[19];
x q[23];
x q[24];
x q[28];
ccx q[19],q[23],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[24],q[43];
ccx q[19],q[23],q[42];
x q[19];
x q[23];
x q[24];
x q[28];

// blue K4 (20, 24, 28, 39) -> avoid connecting to none of the four
x q[19];
x q[23];
x q[27];
x q[38];
ccx q[19],q[23],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[27],q[43];
ccx q[19],q[23],q[42];
x q[19];
x q[23];
x q[27];
x q[38];

// blue K4 (20, 24, 29, 35) -> avoid connecting to none of the four
x q[19];
x q[23];
x q[28];
x q[34];
ccx q[19],q[23],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[28],q[43];
ccx q[19],q[23],q[42];
x q[19];
x q[23];
x q[28];
x q[34];

// blue K4 (20, 24, 35, 39) -> avoid connecting to none of the four
x q[19];
x q[23];
x q[34];
x q[38];
ccx q[19],q[23],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[34],q[43];
ccx q[19],q[23],q[42];
x q[19];
x q[23];
x q[34];
x q[38];

// blue K4 (20, 25, 28, 31) -> avoid connecting to none of the four
x q[19];
x q[24];
x q[27];
x q[30];
ccx q[19],q[24],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[27],q[43];
ccx q[19],q[24],q[42];
x q[19];
x q[24];
x q[27];
x q[30];

// blue K4 (20, 26, 29, 35) -> avoid connecting to none of the four
x q[19];
x q[25];
x q[28];
x q[34];
ccx q[19],q[25],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[28],q[43];
ccx q[19],q[25],q[42];
x q[19];
x q[25];
x q[28];
x q[34];

// blue K4 (20, 26, 29, 37) -> avoid connecting to none of the four
x q[19];
x q[25];
x q[28];
x q[36];
ccx q[19],q[25],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[28],q[43];
ccx q[19],q[25],q[42];
x q[19];
x q[25];
x q[28];
x q[36];

// blue K4 (20, 26, 31, 35) -> avoid connecting to none of the four
x q[19];
x q[25];
x q[30];
x q[34];
ccx q[19],q[25],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[30],q[43];
ccx q[19],q[25],q[42];
x q[19];
x q[25];
x q[30];
x q[34];

// blue K4 (20, 26, 31, 37) -> avoid connecting to none of the four
x q[19];
x q[25];
x q[30];
x q[36];
ccx q[19],q[25],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[30],q[43];
ccx q[19],q[25],q[42];
x q[19];
x q[25];
x q[30];
x q[36];

// blue K4 (20, 28, 31, 37) -> avoid connecting to none of the four
x q[19];
x q[27];
x q[30];
x q[36];
ccx q[19],q[27],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[30],q[43];
ccx q[19],q[27],q[42];
x q[19];
x q[27];
x q[30];
x q[36];

// blue K4 (20, 28, 31, 39) -> avoid connecting to none of the four
x q[19];
x q[27];
x q[30];
x q[38];
ccx q[19],q[27],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[30],q[43];
ccx q[19],q[27],q[42];
x q[19];
x q[27];
x q[30];
x q[38];

// blue K4 (20, 31, 35, 39) -> avoid connecting to none of the four
x q[19];
x q[30];
x q[34];
x q[38];
ccx q[19],q[30],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[34],q[43];
ccx q[19],q[30],q[42];
x q[19];
x q[30];
x q[34];
x q[38];

// blue K4 (21, 24, 25, 29) -> avoid connecting to none of the four
x q[20];
x q[23];
x q[24];
x q[28];
ccx q[20],q[23],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[28],q[44];
rz(0.15) q[44];
ccx q[43],q[28],q[44];
ccx q[42],q[24],q[43];
ccx q[20],q[23],q[42];
x q[20];
x q[23];
x q[24];
x q[28];

// blue K4 (21, 24, 25, 30) -> avoid connecting to none of the four
x q[20];
x q[23];
x q[24];
x q[29];
ccx q[20],q[23],q[42];
ccx q[42],q[24],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[24],q[43];
ccx q[20],q[23],q[42];
x q[20];
x q[23];
x q[24];
x q[29];

// blue K4 (21, 24, 27, 30) -> avoid connecting to none of the four
x q[20];
x q[23];
x q[26];
x q[29];
ccx q[20],q[23],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[29],q[44];
rz(0.15) q[44];
ccx q[43],q[29],q[44];
ccx q[42],q[26],q[43];
ccx q[20],q[23],q[42];
x q[20];
x q[23];
x q[26];
x q[29];

// blue K4 (21, 24, 27, 32) -> avoid connecting to none of the four
x q[20];
x q[23];
x q[26];
x q[31];
ccx q[20],q[23],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[26],q[43];
ccx q[20],q[23],q[42];
x q[20];
x q[23];
x q[26];
x q[31];

// blue K4 (21, 24, 29, 32) -> avoid connecting to none of the four
x q[20];
x q[23];
x q[28];
x q[31];
ccx q[20],q[23],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[28],q[43];
ccx q[20],q[23],q[42];
x q[20];
x q[23];
x q[28];
x q[31];

// blue K4 (21, 25, 29, 40) -> avoid connecting to none of the four
x q[20];
x q[24];
x q[28];
x q[39];
ccx q[20],q[24],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[28],q[43];
ccx q[20],q[24],q[42];
x q[20];
x q[24];
x q[28];
x q[39];

// blue K4 (21, 25, 30, 36) -> avoid connecting to none of the four
x q[20];
x q[24];
x q[29];
x q[35];
ccx q[20],q[24],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[29],q[43];
ccx q[20],q[24],q[42];
x q[20];
x q[24];
x q[29];
x q[35];

// blue K4 (21, 25, 36, 40) -> avoid connecting to none of the four
x q[20];
x q[24];
x q[35];
x q[39];
ccx q[20],q[24],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[35],q[43];
ccx q[20],q[24],q[42];
x q[20];
x q[24];
x q[35];
x q[39];

// blue K4 (21, 26, 29, 32) -> avoid connecting to none of the four
x q[20];
x q[25];
x q[28];
x q[31];
ccx q[20],q[25],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[28],q[43];
ccx q[20],q[25],q[42];
x q[20];
x q[25];
x q[28];
x q[31];

// blue K4 (21, 27, 30, 36) -> avoid connecting to none of the four
x q[20];
x q[26];
x q[29];
x q[35];
ccx q[20],q[26],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[29],q[43];
ccx q[20],q[26],q[42];
x q[20];
x q[26];
x q[29];
x q[35];

// blue K4 (21, 27, 30, 38) -> avoid connecting to none of the four
x q[20];
x q[26];
x q[29];
x q[37];
ccx q[20],q[26],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[29],q[43];
ccx q[20],q[26],q[42];
x q[20];
x q[26];
x q[29];
x q[37];

// blue K4 (21, 27, 32, 36) -> avoid connecting to none of the four
x q[20];
x q[26];
x q[31];
x q[35];
ccx q[20],q[26],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[31],q[43];
ccx q[20],q[26],q[42];
x q[20];
x q[26];
x q[31];
x q[35];

// blue K4 (21, 27, 32, 38) -> avoid connecting to none of the four
x q[20];
x q[26];
x q[31];
x q[37];
ccx q[20],q[26],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[31],q[43];
ccx q[20],q[26],q[42];
x q[20];
x q[26];
x q[31];
x q[37];

// blue K4 (21, 29, 32, 38) -> avoid connecting to none of the four
x q[20];
x q[28];
x q[31];
x q[37];
ccx q[20],q[28],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[31],q[43];
ccx q[20],q[28],q[42];
x q[20];
x q[28];
x q[31];
x q[37];

// blue K4 (21, 29, 32, 40) -> avoid connecting to none of the four
x q[20];
x q[28];
x q[31];
x q[39];
ccx q[20],q[28],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[31],q[43];
ccx q[20],q[28],q[42];
x q[20];
x q[28];
x q[31];
x q[39];

// blue K4 (21, 32, 36, 40) -> avoid connecting to none of the four
x q[20];
x q[31];
x q[35];
x q[39];
ccx q[20],q[31],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[35],q[43];
ccx q[20],q[31],q[42];
x q[20];
x q[31];
x q[35];
x q[39];

// blue K4 (22, 25, 28, 31) -> avoid connecting to none of the four
x q[21];
x q[24];
x q[27];
x q[30];
ccx q[21],q[24],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[27],q[43];
ccx q[21],q[24],q[42];
x q[21];
x q[24];
x q[27];
x q[30];

// blue K4 (22, 25, 28, 33) -> avoid connecting to none of the four
x q[21];
x q[24];
x q[27];
x q[32];
ccx q[21],q[24],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[27],q[43];
ccx q[21],q[24],q[42];
x q[21];
x q[24];
x q[27];
x q[32];

// blue K4 (22, 25, 30, 31) -> avoid connecting to none of the four
x q[21];
x q[24];
x q[29];
x q[30];
ccx q[21],q[24],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[29],q[43];
ccx q[21],q[24],q[42];
x q[21];
x q[24];
x q[29];
x q[30];

// blue K4 (22, 25, 30, 33) -> avoid connecting to none of the four
x q[21];
x q[24];
x q[29];
x q[32];
ccx q[21],q[24],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[29],q[43];
ccx q[21],q[24],q[42];
x q[21];
x q[24];
x q[29];
x q[32];

// blue K4 (22, 26, 30, 31) -> avoid connecting to none of the four
x q[21];
x q[25];
x q[29];
x q[30];
ccx q[21],q[25],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[29],q[43];
ccx q[21],q[25],q[42];
x q[21];
x q[25];
x q[29];
x q[30];

// blue K4 (22, 26, 30, 41) -> avoid connecting to none of the four
x q[21];
x q[25];
x q[29];
x q[40];
ccx q[21],q[25],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[29],q[43];
ccx q[21],q[25],q[42];
x q[21];
x q[25];
x q[29];
x q[40];

// blue K4 (22, 26, 31, 37) -> avoid connecting to none of the four
x q[21];
x q[25];
x q[30];
x q[36];
ccx q[21],q[25],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[30],q[43];
ccx q[21],q[25],q[42];
x q[21];
x q[25];
x q[30];
x q[36];

// blue K4 (22, 26, 37, 41) -> avoid connecting to none of the four
x q[21];
x q[25];
x q[36];
x q[40];
ccx q[21],q[25],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[36],q[43];
ccx q[21],q[25],q[42];
x q[21];
x q[25];
x q[36];
x q[40];

// blue K4 (22, 27, 30, 31) -> avoid connecting to none of the four
x q[21];
x q[26];
x q[29];
x q[30];
ccx q[21],q[26],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[30],q[44];
rz(0.15) q[44];
ccx q[43],q[30],q[44];
ccx q[42],q[29],q[43];
ccx q[21],q[26],q[42];
x q[21];
x q[26];
x q[29];
x q[30];

// blue K4 (22, 27, 30, 33) -> avoid connecting to none of the four
x q[21];
x q[26];
x q[29];
x q[32];
ccx q[21],q[26],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[29],q[43];
ccx q[21],q[26],q[42];
x q[21];
x q[26];
x q[29];
x q[32];

// blue K4 (22, 28, 31, 37) -> avoid connecting to none of the four
x q[21];
x q[27];
x q[30];
x q[36];
ccx q[21],q[27],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[30],q[43];
ccx q[21],q[27],q[42];
x q[21];
x q[27];
x q[30];
x q[36];

// blue K4 (22, 28, 31, 39) -> avoid connecting to none of the four
x q[21];
x q[27];
x q[30];
x q[38];
ccx q[21],q[27],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[30],q[43];
ccx q[21],q[27],q[42];
x q[21];
x q[27];
x q[30];
x q[38];

// blue K4 (22, 28, 33, 37) -> avoid connecting to none of the four
x q[21];
x q[27];
x q[32];
x q[36];
ccx q[21],q[27],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[32],q[43];
ccx q[21],q[27],q[42];
x q[21];
x q[27];
x q[32];
x q[36];

// blue K4 (22, 28, 33, 39) -> avoid connecting to none of the four
x q[21];
x q[27];
x q[32];
x q[38];
ccx q[21],q[27],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[32],q[43];
ccx q[21],q[27],q[42];
x q[21];
x q[27];
x q[32];
x q[38];

// blue K4 (22, 30, 31, 39) -> avoid connecting to none of the four
x q[21];
x q[29];
x q[30];
x q[38];
ccx q[21],q[29],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[30],q[43];
ccx q[21],q[29],q[42];
x q[21];
x q[29];
x q[30];
x q[38];

// blue K4 (22, 30, 33, 39) -> avoid connecting to none of the four
x q[21];
x q[29];
x q[32];
x q[38];
ccx q[21],q[29],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[32],q[43];
ccx q[21],q[29],q[42];
x q[21];
x q[29];
x q[32];
x q[38];

// blue K4 (22, 30, 33, 41) -> avoid connecting to none of the four
x q[21];
x q[29];
x q[32];
x q[40];
ccx q[21],q[29],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[32],q[43];
ccx q[21],q[29],q[42];
x q[21];
x q[29];
x q[32];
x q[40];

// blue K4 (22, 33, 37, 41) -> avoid connecting to none of the four
x q[21];
x q[32];
x q[36];
x q[40];
ccx q[21],q[32],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[36],q[43];
ccx q[21],q[32],q[42];
x q[21];
x q[32];
x q[36];
x q[40];

// blue K4 (23, 24, 27, 32) -> avoid connecting to none of the four
x q[22];
x q[23];
x q[26];
x q[31];
ccx q[22],q[23],q[42];
ccx q[42],q[26],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[26],q[43];
ccx q[22],q[23],q[42];
x q[22];
x q[23];
x q[26];
x q[31];

// blue K4 (23, 24, 28, 32) -> avoid connecting to none of the four
x q[22];
x q[23];
x q[27];
x q[31];
ccx q[22],q[23],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[27],q[43];
ccx q[22],q[23],q[42];
x q[22];
x q[23];
x q[27];
x q[31];

// blue K4 (23, 24, 29, 32) -> avoid connecting to none of the four
x q[22];
x q[23];
x q[28];
x q[31];
ccx q[22],q[23],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[28],q[43];
ccx q[22],q[23],q[42];
x q[22];
x q[23];
x q[28];
x q[31];

// blue K4 (23, 26, 29, 32) -> avoid connecting to none of the four
x q[22];
x q[25];
x q[28];
x q[31];
ccx q[22],q[25],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[31],q[44];
rz(0.15) q[44];
ccx q[43],q[31],q[44];
ccx q[42],q[28],q[43];
ccx q[22],q[25],q[42];
x q[22];
x q[25];
x q[28];
x q[31];

// blue K4 (23, 26, 29, 34) -> avoid connecting to none of the four
x q[22];
x q[25];
x q[28];
x q[33];
ccx q[22],q[25],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[28],q[43];
ccx q[22],q[25],q[42];
x q[22];
x q[25];
x q[28];
x q[33];

// blue K4 (23, 26, 31, 34) -> avoid connecting to none of the four
x q[22];
x q[25];
x q[30];
x q[33];
ccx q[22],q[25],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[30],q[43];
ccx q[22],q[25],q[42];
x q[22];
x q[25];
x q[30];
x q[33];

// blue K4 (23, 27, 31, 42) -> avoid connecting to none of the four
x q[22];
x q[26];
x q[30];
x q[41];
ccx q[22],q[26],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[30],q[43];
ccx q[22],q[26],q[42];
x q[22];
x q[26];
x q[30];
x q[41];

// blue K4 (23, 27, 32, 38) -> avoid connecting to none of the four
x q[22];
x q[26];
x q[31];
x q[37];
ccx q[22],q[26],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[31],q[43];
ccx q[22],q[26],q[42];
x q[22];
x q[26];
x q[31];
x q[37];

// blue K4 (23, 27, 38, 42) -> avoid connecting to none of the four
x q[22];
x q[26];
x q[37];
x q[41];
ccx q[22],q[26],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[37],q[43];
ccx q[22],q[26],q[42];
x q[22];
x q[26];
x q[37];
x q[41];

// blue K4 (23, 28, 31, 34) -> avoid connecting to none of the four
x q[22];
x q[27];
x q[30];
x q[33];
ccx q[22],q[27],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[30],q[43];
ccx q[22],q[27],q[42];
x q[22];
x q[27];
x q[30];
x q[33];

// blue K4 (23, 29, 32, 38) -> avoid connecting to none of the four
x q[22];
x q[28];
x q[31];
x q[37];
ccx q[22],q[28],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[31],q[43];
ccx q[22],q[28],q[42];
x q[22];
x q[28];
x q[31];
x q[37];

// blue K4 (23, 29, 32, 40) -> avoid connecting to none of the four
x q[22];
x q[28];
x q[31];
x q[39];
ccx q[22],q[28],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[31],q[43];
ccx q[22],q[28],q[42];
x q[22];
x q[28];
x q[31];
x q[39];

// blue K4 (23, 29, 34, 38) -> avoid connecting to none of the four
x q[22];
x q[28];
x q[33];
x q[37];
ccx q[22],q[28],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[33],q[43];
ccx q[22],q[28],q[42];
x q[22];
x q[28];
x q[33];
x q[37];

// blue K4 (23, 29, 34, 40) -> avoid connecting to none of the four
x q[22];
x q[28];
x q[33];
x q[39];
ccx q[22],q[28],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[33],q[43];
ccx q[22],q[28],q[42];
x q[22];
x q[28];
x q[33];
x q[39];

// blue K4 (23, 31, 34, 40) -> avoid connecting to none of the four
x q[22];
x q[30];
x q[33];
x q[39];
ccx q[22],q[30],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[33],q[43];
ccx q[22],q[30],q[42];
x q[22];
x q[30];
x q[33];
x q[39];

// blue K4 (23, 31, 34, 42) -> avoid connecting to none of the four
x q[22];
x q[30];
x q[33];
x q[41];
ccx q[22],q[30],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[33],q[43];
ccx q[22],q[30],q[42];
x q[22];
x q[30];
x q[33];
x q[41];

// blue K4 (23, 34, 38, 42) -> avoid connecting to none of the four
x q[22];
x q[33];
x q[37];
x q[41];
ccx q[22],q[33],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[37],q[43];
ccx q[22],q[33],q[42];
x q[22];
x q[33];
x q[37];
x q[41];

// blue K4 (24, 25, 28, 33) -> avoid connecting to none of the four
x q[23];
x q[24];
x q[27];
x q[32];
ccx q[23],q[24],q[42];
ccx q[42],q[27],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[27],q[43];
ccx q[23],q[24],q[42];
x q[23];
x q[24];
x q[27];
x q[32];

// blue K4 (24, 25, 29, 33) -> avoid connecting to none of the four
x q[23];
x q[24];
x q[28];
x q[32];
ccx q[23],q[24],q[42];
ccx q[42],q[28],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[28],q[43];
ccx q[23],q[24],q[42];
x q[23];
x q[24];
x q[28];
x q[32];

// blue K4 (24, 25, 30, 33) -> avoid connecting to none of the four
x q[23];
x q[24];
x q[29];
x q[32];
ccx q[23],q[24],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[29],q[43];
ccx q[23],q[24],q[42];
x q[23];
x q[24];
x q[29];
x q[32];

// blue K4 (24, 27, 30, 33) -> avoid connecting to none of the four
x q[23];
x q[26];
x q[29];
x q[32];
ccx q[23],q[26],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[32],q[44];
rz(0.15) q[44];
ccx q[43],q[32],q[44];
ccx q[42],q[29],q[43];
ccx q[23],q[26],q[42];
x q[23];
x q[26];
x q[29];
x q[32];

// blue K4 (24, 27, 30, 35) -> avoid connecting to none of the four
x q[23];
x q[26];
x q[29];
x q[34];
ccx q[23],q[26],q[42];
ccx q[42],q[29],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[29],q[43];
ccx q[23],q[26],q[42];
x q[23];
x q[26];
x q[29];
x q[34];

// blue K4 (24, 27, 32, 35) -> avoid connecting to none of the four
x q[23];
x q[26];
x q[31];
x q[34];
ccx q[23],q[26],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[31],q[43];
ccx q[23],q[26],q[42];
x q[23];
x q[26];
x q[31];
x q[34];

// blue K4 (24, 28, 33, 39) -> avoid connecting to none of the four
x q[23];
x q[27];
x q[32];
x q[38];
ccx q[23],q[27],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[32],q[43];
ccx q[23],q[27],q[42];
x q[23];
x q[27];
x q[32];
x q[38];

// blue K4 (24, 29, 32, 35) -> avoid connecting to none of the four
x q[23];
x q[28];
x q[31];
x q[34];
ccx q[23],q[28],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[31],q[43];
ccx q[23],q[28],q[42];
x q[23];
x q[28];
x q[31];
x q[34];

// blue K4 (24, 30, 33, 39) -> avoid connecting to none of the four
x q[23];
x q[29];
x q[32];
x q[38];
ccx q[23],q[29],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[32],q[43];
ccx q[23],q[29],q[42];
x q[23];
x q[29];
x q[32];
x q[38];

// blue K4 (24, 30, 33, 41) -> avoid connecting to none of the four
x q[23];
x q[29];
x q[32];
x q[40];
ccx q[23],q[29],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[32],q[43];
ccx q[23],q[29],q[42];
x q[23];
x q[29];
x q[32];
x q[40];

// blue K4 (24, 30, 35, 39) -> avoid connecting to none of the four
x q[23];
x q[29];
x q[34];
x q[38];
ccx q[23],q[29],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[34],q[43];
ccx q[23],q[29],q[42];
x q[23];
x q[29];
x q[34];
x q[38];

// blue K4 (24, 30, 35, 41) -> avoid connecting to none of the four
x q[23];
x q[29];
x q[34];
x q[40];
ccx q[23],q[29],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[34],q[43];
ccx q[23],q[29],q[42];
x q[23];
x q[29];
x q[34];
x q[40];

// blue K4 (24, 32, 35, 41) -> avoid connecting to none of the four
x q[23];
x q[31];
x q[34];
x q[40];
ccx q[23],q[31],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[34],q[43];
ccx q[23],q[31],q[42];
x q[23];
x q[31];
x q[34];
x q[40];

// blue K4 (25, 28, 31, 34) -> avoid connecting to none of the four
x q[24];
x q[27];
x q[30];
x q[33];
ccx q[24],q[27],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[30],q[43];
ccx q[24],q[27],q[42];
x q[24];
x q[27];
x q[30];
x q[33];

// blue K4 (25, 28, 31, 36) -> avoid connecting to none of the four
x q[24];
x q[27];
x q[30];
x q[35];
ccx q[24],q[27],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[30],q[43];
ccx q[24],q[27],q[42];
x q[24];
x q[27];
x q[30];
x q[35];

// blue K4 (25, 28, 33, 34) -> avoid connecting to none of the four
x q[24];
x q[27];
x q[32];
x q[33];
ccx q[24],q[27],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[32],q[43];
ccx q[24],q[27],q[42];
x q[24];
x q[27];
x q[32];
x q[33];

// blue K4 (25, 28, 33, 36) -> avoid connecting to none of the four
x q[24];
x q[27];
x q[32];
x q[35];
ccx q[24],q[27],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[32],q[43];
ccx q[24],q[27],q[42];
x q[24];
x q[27];
x q[32];
x q[35];

// blue K4 (25, 29, 33, 34) -> avoid connecting to none of the four
x q[24];
x q[28];
x q[32];
x q[33];
ccx q[24],q[28],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[32],q[43];
ccx q[24],q[28],q[42];
x q[24];
x q[28];
x q[32];
x q[33];

// blue K4 (25, 29, 34, 40) -> avoid connecting to none of the four
x q[24];
x q[28];
x q[33];
x q[39];
ccx q[24],q[28],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[33],q[43];
ccx q[24],q[28],q[42];
x q[24];
x q[28];
x q[33];
x q[39];

// blue K4 (25, 30, 31, 34) -> avoid connecting to none of the four
x q[24];
x q[29];
x q[30];
x q[33];
ccx q[24],q[29],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[30],q[43];
ccx q[24],q[29],q[42];
x q[24];
x q[29];
x q[30];
x q[33];

// blue K4 (25, 30, 31, 36) -> avoid connecting to none of the four
x q[24];
x q[29];
x q[30];
x q[35];
ccx q[24],q[29],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[30],q[43];
ccx q[24],q[29],q[42];
x q[24];
x q[29];
x q[30];
x q[35];

// blue K4 (25, 30, 33, 34) -> avoid connecting to none of the four
x q[24];
x q[29];
x q[32];
x q[33];
ccx q[24],q[29],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[32],q[43];
ccx q[24],q[29],q[42];
x q[24];
x q[29];
x q[32];
x q[33];

// blue K4 (25, 30, 33, 36) -> avoid connecting to none of the four
x q[24];
x q[29];
x q[32];
x q[35];
ccx q[24],q[29],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[32],q[43];
ccx q[24],q[29],q[42];
x q[24];
x q[29];
x q[32];
x q[35];

// blue K4 (25, 31, 34, 40) -> avoid connecting to none of the four
x q[24];
x q[30];
x q[33];
x q[39];
ccx q[24],q[30],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[33],q[43];
ccx q[24],q[30],q[42];
x q[24];
x q[30];
x q[33];
x q[39];

// blue K4 (25, 31, 34, 42) -> avoid connecting to none of the four
x q[24];
x q[30];
x q[33];
x q[41];
ccx q[24],q[30],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[33],q[43];
ccx q[24],q[30],q[42];
x q[24];
x q[30];
x q[33];
x q[41];

// blue K4 (25, 31, 36, 40) -> avoid connecting to none of the four
x q[24];
x q[30];
x q[35];
x q[39];
ccx q[24],q[30],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[35],q[43];
ccx q[24],q[30],q[42];
x q[24];
x q[30];
x q[35];
x q[39];

// blue K4 (25, 31, 36, 42) -> avoid connecting to none of the four
x q[24];
x q[30];
x q[35];
x q[41];
ccx q[24],q[30],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[35],q[43];
ccx q[24],q[30],q[42];
x q[24];
x q[30];
x q[35];
x q[41];

// blue K4 (25, 33, 34, 42) -> avoid connecting to none of the four
x q[24];
x q[32];
x q[33];
x q[41];
ccx q[24],q[32],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[33],q[43];
ccx q[24],q[32],q[42];
x q[24];
x q[32];
x q[33];
x q[41];

// blue K4 (25, 33, 36, 42) -> avoid connecting to none of the four
x q[24];
x q[32];
x q[35];
x q[41];
ccx q[24],q[32],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[35],q[43];
ccx q[24],q[32],q[42];
x q[24];
x q[32];
x q[35];
x q[41];

// blue K4 (26, 29, 32, 35) -> avoid connecting to none of the four
x q[25];
x q[28];
x q[31];
x q[34];
ccx q[25],q[28],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[31],q[43];
ccx q[25],q[28],q[42];
x q[25];
x q[28];
x q[31];
x q[34];

// blue K4 (26, 29, 32, 37) -> avoid connecting to none of the four
x q[25];
x q[28];
x q[31];
x q[36];
ccx q[25],q[28],q[42];
ccx q[42],q[31],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[31],q[43];
ccx q[25],q[28],q[42];
x q[25];
x q[28];
x q[31];
x q[36];

// blue K4 (26, 29, 34, 37) -> avoid connecting to none of the four
x q[25];
x q[28];
x q[33];
x q[36];
ccx q[25],q[28],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[33],q[43];
ccx q[25],q[28],q[42];
x q[25];
x q[28];
x q[33];
x q[36];

// blue K4 (26, 30, 31, 34) -> avoid connecting to none of the four
x q[25];
x q[29];
x q[30];
x q[33];
ccx q[25],q[29],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[33],q[44];
rz(0.15) q[44];
ccx q[43],q[33],q[44];
ccx q[42],q[30],q[43];
ccx q[25],q[29],q[42];
x q[25];
x q[29];
x q[30];
x q[33];

// blue K4 (26, 30, 31, 35) -> avoid connecting to none of the four
x q[25];
x q[29];
x q[30];
x q[34];
ccx q[25],q[29],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[30],q[43];
ccx q[25],q[29],q[42];
x q[25];
x q[29];
x q[30];
x q[34];

// blue K4 (26, 30, 35, 41) -> avoid connecting to none of the four
x q[25];
x q[29];
x q[34];
x q[40];
ccx q[25],q[29],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[34],q[43];
ccx q[25],q[29],q[42];
x q[25];
x q[29];
x q[34];
x q[40];

// blue K4 (26, 31, 34, 37) -> avoid connecting to none of the four
x q[25];
x q[30];
x q[33];
x q[36];
ccx q[25],q[30],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[33],q[43];
ccx q[25],q[30],q[42];
x q[25];
x q[30];
x q[33];
x q[36];

// blue K4 (26, 32, 35, 41) -> avoid connecting to none of the four
x q[25];
x q[31];
x q[34];
x q[40];
ccx q[25],q[31],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[34],q[43];
ccx q[25],q[31],q[42];
x q[25];
x q[31];
x q[34];
x q[40];

// blue K4 (26, 32, 37, 41) -> avoid connecting to none of the four
x q[25];
x q[31];
x q[36];
x q[40];
ccx q[25],q[31],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[36],q[43];
ccx q[25],q[31],q[42];
x q[25];
x q[31];
x q[36];
x q[40];

// blue K4 (27, 30, 31, 35) -> avoid connecting to none of the four
x q[26];
x q[29];
x q[30];
x q[34];
ccx q[26],q[29],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[34],q[44];
rz(0.15) q[44];
ccx q[43],q[34],q[44];
ccx q[42],q[30],q[43];
ccx q[26],q[29],q[42];
x q[26];
x q[29];
x q[30];
x q[34];

// blue K4 (27, 30, 31, 36) -> avoid connecting to none of the four
x q[26];
x q[29];
x q[30];
x q[35];
ccx q[26],q[29],q[42];
ccx q[42],q[30],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[30],q[43];
ccx q[26],q[29],q[42];
x q[26];
x q[29];
x q[30];
x q[35];

// blue K4 (27, 30, 33, 36) -> avoid connecting to none of the four
x q[26];
x q[29];
x q[32];
x q[35];
ccx q[26],q[29],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[35],q[44];
rz(0.15) q[44];
ccx q[43],q[35],q[44];
ccx q[42],q[32],q[43];
ccx q[26],q[29],q[42];
x q[26];
x q[29];
x q[32];
x q[35];

// blue K4 (27, 30, 33, 38) -> avoid connecting to none of the four
x q[26];
x q[29];
x q[32];
x q[37];
ccx q[26],q[29],q[42];
ccx q[42],q[32],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[32],q[43];
ccx q[26],q[29],q[42];
x q[26];
x q[29];
x q[32];
x q[37];

// blue K4 (27, 30, 35, 38) -> avoid connecting to none of the four
x q[26];
x q[29];
x q[34];
x q[37];
ccx q[26],q[29],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[34],q[43];
ccx q[26],q[29],q[42];
x q[26];
x q[29];
x q[34];
x q[37];

// blue K4 (27, 31, 36, 42) -> avoid connecting to none of the four
x q[26];
x q[30];
x q[35];
x q[41];
ccx q[26],q[30],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[35],q[43];
ccx q[26],q[30],q[42];
x q[26];
x q[30];
x q[35];
x q[41];

// blue K4 (27, 32, 35, 38) -> avoid connecting to none of the four
x q[26];
x q[31];
x q[34];
x q[37];
ccx q[26],q[31],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[34],q[43];
ccx q[26],q[31],q[42];
x q[26];
x q[31];
x q[34];
x q[37];

// blue K4 (27, 33, 36, 42) -> avoid connecting to none of the four
x q[26];
x q[32];
x q[35];
x q[41];
ccx q[26],q[32],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[35],q[43];
ccx q[26],q[32],q[42];
x q[26];
x q[32];
x q[35];
x q[41];

// blue K4 (27, 33, 38, 42) -> avoid connecting to none of the four
x q[26];
x q[32];
x q[37];
x q[41];
ccx q[26],q[32],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[37],q[43];
ccx q[26],q[32],q[42];
x q[26];
x q[32];
x q[37];
x q[41];

// blue K4 (28, 31, 34, 37) -> avoid connecting to none of the four
x q[27];
x q[30];
x q[33];
x q[36];
ccx q[27],q[30],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[33],q[43];
ccx q[27],q[30],q[42];
x q[27];
x q[30];
x q[33];
x q[36];

// blue K4 (28, 31, 34, 39) -> avoid connecting to none of the four
x q[27];
x q[30];
x q[33];
x q[38];
ccx q[27],q[30],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[33],q[43];
ccx q[27],q[30],q[42];
x q[27];
x q[30];
x q[33];
x q[38];

// blue K4 (28, 31, 36, 39) -> avoid connecting to none of the four
x q[27];
x q[30];
x q[35];
x q[38];
ccx q[27],q[30],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[35],q[43];
ccx q[27],q[30],q[42];
x q[27];
x q[30];
x q[35];
x q[38];

// blue K4 (28, 33, 34, 37) -> avoid connecting to none of the four
x q[27];
x q[32];
x q[33];
x q[36];
ccx q[27],q[32],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[33],q[43];
ccx q[27],q[32],q[42];
x q[27];
x q[32];
x q[33];
x q[36];

// blue K4 (28, 33, 34, 39) -> avoid connecting to none of the four
x q[27];
x q[32];
x q[33];
x q[38];
ccx q[27],q[32],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[33],q[43];
ccx q[27],q[32],q[42];
x q[27];
x q[32];
x q[33];
x q[38];

// blue K4 (28, 33, 36, 39) -> avoid connecting to none of the four
x q[27];
x q[32];
x q[35];
x q[38];
ccx q[27],q[32],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[35],q[43];
ccx q[27],q[32],q[42];
x q[27];
x q[32];
x q[35];
x q[38];

// blue K4 (29, 32, 35, 38) -> avoid connecting to none of the four
x q[28];
x q[31];
x q[34];
x q[37];
ccx q[28],q[31],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[34],q[43];
ccx q[28],q[31],q[42];
x q[28];
x q[31];
x q[34];
x q[37];

// blue K4 (29, 32, 35, 40) -> avoid connecting to none of the four
x q[28];
x q[31];
x q[34];
x q[39];
ccx q[28],q[31],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[34],q[43];
ccx q[28],q[31],q[42];
x q[28];
x q[31];
x q[34];
x q[39];

// blue K4 (29, 32, 37, 40) -> avoid connecting to none of the four
x q[28];
x q[31];
x q[36];
x q[39];
ccx q[28],q[31],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[36],q[43];
ccx q[28],q[31],q[42];
x q[28];
x q[31];
x q[36];
x q[39];

// blue K4 (29, 33, 34, 37) -> avoid connecting to none of the four
x q[28];
x q[32];
x q[33];
x q[36];
ccx q[28],q[32],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[36],q[44];
rz(0.15) q[44];
ccx q[43],q[36],q[44];
ccx q[42],q[33],q[43];
ccx q[28],q[32],q[42];
x q[28];
x q[32];
x q[33];
x q[36];

// blue K4 (29, 33, 34, 38) -> avoid connecting to none of the four
x q[28];
x q[32];
x q[33];
x q[37];
ccx q[28],q[32],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[33],q[43];
ccx q[28],q[32],q[42];
x q[28];
x q[32];
x q[33];
x q[37];

// blue K4 (29, 34, 37, 40) -> avoid connecting to none of the four
x q[28];
x q[33];
x q[36];
x q[39];
ccx q[28],q[33],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[36],q[43];
ccx q[28],q[33],q[42];
x q[28];
x q[33];
x q[36];
x q[39];

// blue K4 (30, 31, 34, 39) -> avoid connecting to none of the four
x q[29];
x q[30];
x q[33];
x q[38];
ccx q[29],q[30],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[33],q[43];
ccx q[29],q[30],q[42];
x q[29];
x q[30];
x q[33];
x q[38];

// blue K4 (30, 31, 35, 39) -> avoid connecting to none of the four
x q[29];
x q[30];
x q[34];
x q[38];
ccx q[29],q[30],q[42];
ccx q[42],q[34],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[34],q[43];
ccx q[29],q[30],q[42];
x q[29];
x q[30];
x q[34];
x q[38];

// blue K4 (30, 31, 36, 39) -> avoid connecting to none of the four
x q[29];
x q[30];
x q[35];
x q[38];
ccx q[29],q[30],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[35],q[43];
ccx q[29],q[30],q[42];
x q[29];
x q[30];
x q[35];
x q[38];

// blue K4 (30, 33, 34, 38) -> avoid connecting to none of the four
x q[29];
x q[32];
x q[33];
x q[37];
ccx q[29],q[32],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[37],q[44];
rz(0.15) q[44];
ccx q[43],q[37],q[44];
ccx q[42],q[33],q[43];
ccx q[29],q[32],q[42];
x q[29];
x q[32];
x q[33];
x q[37];

// blue K4 (30, 33, 34, 39) -> avoid connecting to none of the four
x q[29];
x q[32];
x q[33];
x q[38];
ccx q[29],q[32],q[42];
ccx q[42],q[33],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[33],q[43];
ccx q[29],q[32],q[42];
x q[29];
x q[32];
x q[33];
x q[38];

// blue K4 (30, 33, 36, 39) -> avoid connecting to none of the four
x q[29];
x q[32];
x q[35];
x q[38];
ccx q[29],q[32],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[38],q[44];
rz(0.15) q[44];
ccx q[43],q[38],q[44];
ccx q[42],q[35],q[43];
ccx q[29],q[32],q[42];
x q[29];
x q[32];
x q[35];
x q[38];

// blue K4 (30, 33, 36, 41) -> avoid connecting to none of the four
x q[29];
x q[32];
x q[35];
x q[40];
ccx q[29],q[32],q[42];
ccx q[42],q[35],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[35],q[43];
ccx q[29],q[32],q[42];
x q[29];
x q[32];
x q[35];
x q[40];

// blue K4 (30, 33, 38, 41) -> avoid connecting to none of the four
x q[29];
x q[32];
x q[37];
x q[40];
ccx q[29],q[32],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[37],q[43];
ccx q[29],q[32],q[42];
x q[29];
x q[32];
x q[37];
x q[40];

// blue K4 (30, 35, 38, 41) -> avoid connecting to none of the four
x q[29];
x q[34];
x q[37];
x q[40];
ccx q[29],q[34],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[37],q[43];
ccx q[29],q[34],q[42];
x q[29];
x q[34];
x q[37];
x q[40];

// blue K4 (31, 34, 37, 40) -> avoid connecting to none of the four
x q[30];
x q[33];
x q[36];
x q[39];
ccx q[30],q[33],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[36],q[43];
ccx q[30],q[33],q[42];
x q[30];
x q[33];
x q[36];
x q[39];

// blue K4 (31, 34, 37, 42) -> avoid connecting to none of the four
x q[30];
x q[33];
x q[36];
x q[41];
ccx q[30],q[33],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[36],q[43];
ccx q[30],q[33],q[42];
x q[30];
x q[33];
x q[36];
x q[41];

// blue K4 (31, 34, 39, 40) -> avoid connecting to none of the four
x q[30];
x q[33];
x q[38];
x q[39];
ccx q[30],q[33],q[42];
ccx q[42],q[38],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[38],q[43];
ccx q[30],q[33],q[42];
x q[30];
x q[33];
x q[38];
x q[39];

// blue K4 (31, 34, 39, 42) -> avoid connecting to none of the four
x q[30];
x q[33];
x q[38];
x q[41];
ccx q[30],q[33],q[42];
ccx q[42],q[38],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[38],q[43];
ccx q[30],q[33],q[42];
x q[30];
x q[33];
x q[38];
x q[41];

// blue K4 (31, 35, 39, 40) -> avoid connecting to none of the four
x q[30];
x q[34];
x q[38];
x q[39];
ccx q[30],q[34],q[42];
ccx q[42],q[38],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[38],q[43];
ccx q[30],q[34],q[42];
x q[30];
x q[34];
x q[38];
x q[39];

// blue K4 (31, 36, 39, 40) -> avoid connecting to none of the four
x q[30];
x q[35];
x q[38];
x q[39];
ccx q[30],q[35],q[42];
ccx q[42],q[38],q[43];
ccx q[43],q[39],q[44];
rz(0.15) q[44];
ccx q[43],q[39],q[44];
ccx q[42],q[38],q[43];
ccx q[30],q[35],q[42];
x q[30];
x q[35];
x q[38];
x q[39];

// blue K4 (31, 36, 39, 42) -> avoid connecting to none of the four
x q[30];
x q[35];
x q[38];
x q[41];
ccx q[30],q[35],q[42];
ccx q[42],q[38],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[38],q[43];
ccx q[30],q[35],q[42];
x q[30];
x q[35];
x q[38];
x q[41];

// blue K4 (32, 35, 38, 41) -> avoid connecting to none of the four
x q[31];
x q[34];
x q[37];
x q[40];
ccx q[31],q[34],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[37],q[43];
ccx q[31],q[34],q[42];
x q[31];
x q[34];
x q[37];
x q[40];

// blue K4 (32, 35, 40, 41) -> avoid connecting to none of the four
x q[31];
x q[34];
x q[39];
x q[40];
ccx q[31],q[34],q[42];
ccx q[42],q[39],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[39],q[43];
ccx q[31],q[34],q[42];
x q[31];
x q[34];
x q[39];
x q[40];

// blue K4 (32, 36, 40, 41) -> avoid connecting to none of the four
x q[31];
x q[35];
x q[39];
x q[40];
ccx q[31],q[35],q[42];
ccx q[42],q[39],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[39],q[43];
ccx q[31],q[35],q[42];
x q[31];
x q[35];
x q[39];
x q[40];

// blue K4 (32, 37, 40, 41) -> avoid connecting to none of the four
x q[31];
x q[36];
x q[39];
x q[40];
ccx q[31],q[36],q[42];
ccx q[42],q[39],q[43];
ccx q[43],q[40],q[44];
rz(0.15) q[44];
ccx q[43],q[40],q[44];
ccx q[42],q[39],q[43];
ccx q[31],q[36],q[42];
x q[31];
x q[36];
x q[39];
x q[40];

// blue K4 (33, 34, 37, 42) -> avoid connecting to none of the four
x q[32];
x q[33];
x q[36];
x q[41];
ccx q[32],q[33],q[42];
ccx q[42],q[36],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[36],q[43];
ccx q[32],q[33],q[42];
x q[32];
x q[33];
x q[36];
x q[41];

// blue K4 (33, 34, 38, 42) -> avoid connecting to none of the four
x q[32];
x q[33];
x q[37];
x q[41];
ccx q[32],q[33],q[42];
ccx q[42],q[37],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[37],q[43];
ccx q[32],q[33],q[42];
x q[32];
x q[33];
x q[37];
x q[41];

// blue K4 (33, 34, 39, 42) -> avoid connecting to none of the four
x q[32];
x q[33];
x q[38];
x q[41];
ccx q[32],q[33],q[42];
ccx q[42],q[38],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[38],q[43];
ccx q[32],q[33],q[42];
x q[32];
x q[33];
x q[38];
x q[41];

// blue K4 (33, 36, 39, 42) -> avoid connecting to none of the four
x q[32];
x q[35];
x q[38];
x q[41];
ccx q[32],q[35],q[42];
ccx q[42],q[38],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[38],q[43];
ccx q[32],q[35],q[42];
x q[32];
x q[35];
x q[38];
x q[41];

// blue K4 (33, 36, 41, 42) -> avoid connecting to none of the four
x q[32];
x q[35];
x q[40];
x q[41];
ccx q[32],q[35],q[42];
ccx q[42],q[40],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[40],q[43];
ccx q[32],q[35],q[42];
x q[32];
x q[35];
x q[40];
x q[41];

// blue K4 (33, 37, 41, 42) -> avoid connecting to none of the four
x q[32];
x q[36];
x q[40];
x q[41];
ccx q[32],q[36],q[42];
ccx q[42],q[40],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[40],q[43];
ccx q[32],q[36],q[42];
x q[32];
x q[36];
x q[40];
x q[41];

// blue K4 (33, 38, 41, 42) -> avoid connecting to none of the four
x q[32];
x q[37];
x q[40];
x q[41];
ccx q[32],q[37],q[42];
ccx q[42],q[40],q[43];
ccx q[43],q[41],q[44];
rz(0.15) q[44];
ccx q[43],q[41],q[44];
ccx q[42],q[40],q[43];
ccx q[32],q[37],q[42];
x q[32];
x q[37];
x q[40];
x q[41];

// ---- mixer unitary exp(-i*beta*sum X_i) on the 42 connection qubits ----
rx(0.6) q[0];
rx(0.6) q[1];
rx(0.6) q[2];
rx(0.6) q[3];
rx(0.6) q[4];
rx(0.6) q[5];
rx(0.6) q[6];
rx(0.6) q[7];
rx(0.6) q[8];
rx(0.6) q[9];
rx(0.6) q[10];
rx(0.6) q[11];
rx(0.6) q[12];
rx(0.6) q[13];
rx(0.6) q[14];
rx(0.6) q[15];
rx(0.6) q[16];
rx(0.6) q[17];
rx(0.6) q[18];
rx(0.6) q[19];
rx(0.6) q[20];
rx(0.6) q[21];
rx(0.6) q[22];
rx(0.6) q[23];
rx(0.6) q[24];
rx(0.6) q[25];
rx(0.6) q[26];
rx(0.6) q[27];
rx(0.6) q[28];
rx(0.6) q[29];
rx(0.6) q[30];
rx(0.6) q[31];
rx(0.6) q[32];
rx(0.6) q[33];
rx(0.6) q[34];
rx(0.6) q[35];
rx(0.6) q[36];
rx(0.6) q[37];
rx(0.6) q[38];
rx(0.6) q[39];
rx(0.6) q[40];
rx(0.6) q[41];

// ---- measurement (connection qubits only; ancillas discarded, always |0>) ----
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];
measure q[7] -> c[7];
measure q[8] -> c[8];
measure q[9] -> c[9];
measure q[10] -> c[10];
measure q[11] -> c[11];
measure q[12] -> c[12];
measure q[13] -> c[13];
measure q[14] -> c[14];
measure q[15] -> c[15];
measure q[16] -> c[16];
measure q[17] -> c[17];
measure q[18] -> c[18];
measure q[19] -> c[19];
measure q[20] -> c[20];
measure q[21] -> c[21];
measure q[22] -> c[22];
measure q[23] -> c[23];
measure q[24] -> c[24];
measure q[25] -> c[25];
measure q[26] -> c[26];
measure q[27] -> c[27];
measure q[28] -> c[28];
measure q[29] -> c[29];
measure q[30] -> c[30];
measure q[31] -> c[31];
measure q[32] -> c[32];
measure q[33] -> c[33];
measure q[34] -> c[34];
measure q[35] -> c[35];
measure q[36] -> c[36];
measure q[37] -> c[37];
measure q[38] -> c[38];
measure q[39] -> c[39];
measure q[40] -> c[40];
measure q[41] -> c[41];