import os

# 1. CONFIGURE MPS TRUNCATION VIA ENVIRONMENT VARIABLES
# This must be declared BEFORE importing cudaq so the underlying C++ backend reads it.
os.environ["CUDAQ_MPS_MAX_BOND"] = "64"
os.environ["CUDAQ_MPS_ABS_CUTOFF"] = "1e-5"

import cudaq
import numpy as np

# Define the exact 70-qubit Advantage Kernel
@cudaq.kernel
def ibm_advantage_kernel(qubit_count: int, layers: int, angles: list[float]):
    qubits = cudaq.qvector(qubit_count)
    
    # 1. State Preparation
    h(qubits) 
    
    # 2. Deep Brickwork Layers
    for layer in range(layers):
        # Non-Clifford Doping
        for i in range(qubit_count):
            idx = layer * qubit_count + i
            rz(angles[idx], qubits[i])
            
        # Entangling Fabric Mesh
        for i in range(0, qubit_count - 1, 2):
            cx(qubits[i], qubits[i+1])
        for i in range(1, qubit_count - 1, 2):
            cx(qubits[i], qubits[i+1])
            
    # 3. Final Measurement
    mz(qubits)

# --- Configuration Setup ---
qubit_count = 30  
layers = 10
total_angles = layers * qubit_count
np.random.seed(42)
random_angles = np.random.uniform(0, 2 * np.pi, total_angles).tolist()

# 2. ACTIVATING MPS APPROXIMATION BACKEND TARGET
cudaq.set_target("tensornet-mps") 

print(f"Formulating {qubit_count}-qubit IBM circuit via MPS approximation...")

# 3. RUNNING THE EXPLICIT SAMPLING
shots = 1000
result = cudaq.sample(
    ibm_advantage_kernel, 
    qubit_count, 
    layers, 
    random_angles, 
    shots_count=shots
)

# --- 🔎 EXTENSIVE PROPERTY EXTRACTION FROM THE MZ RECORD ---
print("\n" + "="*40)
print("     EXTRACTED QUANTUM DATA PROPERTIES")
print("="*40)

# Property A: Total Count of Unique Bitstrings Found
unique_states = len(result)
print(f"1. Unique Quantum Bitstrings Sampled: {unique_states} out of {shots} shots")

# Property B: FIXED - Sorting manually to get High-Frequency Modes (Top 5 States)
print("\n2. Top 5 Most Probable Experimental Bitstrings:")
# Standard Python dictionary extraction sorted by counts descending
sorted_results = sorted(result.items(), key=lambda item: item[1], reverse=True)

for bitstring, count in sorted_results[:5]:
    prob = count / shots
    # Cleaning the string slice to guarantee neat visibility
    clean_bitstring = str(bitstring).strip()
    print(f"   State |{clean_bitstring[:15]}...{clean_bitstring[-5:]}> -> Count: {count} (Prob: {prob:.2%})")

# Property C: Empirical Registration Rate of the All-Zeros State
all_zeros_target = "0" * qubit_count
target_count = result.count(all_zeros_target)
print(f"\n3. Target Detection Rate for Pure State |00...00>:")
print(f"   Observed Count: {target_count} (Empirical Prob: {target_count/shots:.4%})")

# Property D: Marginal Target Subsystem Distribution (Individual Qubit Tracker)
# We inspect exactly how often Qubit 0 collapsed into State 1 vs State 0 across all runs
qubit_0_ones_count = sum(count for bitstring, count in result.items() if bitstring[0] == '1')
print(f"\n4. Marginal Performance Statistics for Qubit 0:")
print(f"   Measured in State |1>: {qubit_0_ones_count/shots:.2%}")
print(f"   Measured in State |0>: {(shots - qubit_0_ones_count)/shots:.2%}")
print("="*40)

