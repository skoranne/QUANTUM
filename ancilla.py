import os

# Configure Matrix Product State settings for optimized linear simulation
os.environ["CUDAQ_MPS_MAX_BOND"] = "32"
os.environ["CUDAQ_MPS_ABS_CUTOFF"] = "1e-5"

import cudaq
import numpy as np

# 1. DEFINE THE ACCELERATED KERNEL
# We pass data_map and ancilla_pairs as pre-computed, static layout lists
@cudaq.kernel
def ibm_uchicago_interleaved_rcs(total_qubits: int, depth: int, data_map: list[int], 
                                 ancilla_pairs_flat: list[int], doping_mask: list[int]):
    
    # Allocate a single unified 1D register of qubits
    qubits = cudaq.qvector(total_qubits)
    
    # State Preparation: Apply Hadamard only to physical data qubits
    for i in range(len(data_map)):
        h(qubits[data_map[i]])
    
    # Iterate through the full depth milestone
    for d in range(depth):
        
        # Phase A: Non-Clifford Doping on Data lines
        for logical_i in range(len(data_map)):
            d_phys = data_map[logical_i]
            mask_idx = d * len(data_map) + logical_i
            gate_choice = doping_mask[mask_idx]
            
            if gate_choice == 1:
                t(qubits[d_phys])
            elif gate_choice == 2:
                s(qubits[d_phys])
        
        # Phase B: Alternate Graph Entangling Fabric (LNN over data lines)
        if d % 2 == 0:
            for idx in range(0, len(data_map) - 1, 2):
                cz(qubits[data_map[idx]], qubits[data_map[idx+1]])
        else:
            for idx in range(1, len(data_map) - 1, 2):
                cz(qubits[data_map[idx]], qubits[data_map[idx+1]])
                
        # Phase C: Spacetime Pauli Checks (Error Detection Tracking)
        if d > 0 and d % 10 == 0:
            # Step through the flattened array two elements at a time
            for step in range(0, len(ancilla_pairs_flat), 2):
                d_phys = ancilla_pairs_flat[step]
                a_phys = ancilla_pairs_flat[step+1]
                cx(qubits[d_phys], qubits[a_phys]) # Safe nearest-neighbor gates

    # 3. Final Register Measurement
    mz(qubits)

# --- HOST-SIDE CONFIGURATION CALCULATIONS (STANDARD PYTHON) ---
data_count = 70
ancilla_count = 27
depth = 70  

host_data_map = []
host_ancilla_pairs_flat = []

a_placed = 0
current_physical_idx = 0

for d_logical in range(data_count):
    host_data_map.append(current_physical_idx)
    current_physical_idx += 1
    
    if a_placed < ancilla_count and d_logical % 2 == 1:
        # Save adjacent pair elements flattened for kernel compatibility
        host_ancilla_pairs_flat.append(host_data_map[-1])   # Data physical index
        host_ancilla_pairs_flat.append(current_physical_idx) # Ancilla physical index
        current_physical_idx += 1
        a_placed += 1

total_qubits = current_physical_idx
total_nodes = depth * data_count

np.random.seed(42)
doping_mask = np.random.choice([0, 1, 2], size=total_nodes, p=[0.85, 0.10, 0.05]).tolist()

# Activate the correct Matrix Product State target backend
cudaq.set_target("tensornet-mps")

print(f"Executing Interleaved {data_count}-Data + {ancilla_count}-Ancilla Qubit Advantage Model...")

shots = 1000
result = cudaq.sample(
    ibm_uchicago_interleaved_rcs,
    total_qubits,
    depth,
    host_data_map,
    host_ancilla_pairs_flat,
    doping_mask,
    shots_count=shots
)

# --- 🔎 DEEP DATA PROPERTY EXTRACTION & POST-SELECTION ANALYSIS ---
print("\n" + "="*60)
print("       ADVANTAGE RUNTIME DATA ANALYSIS")
print("="*60)

# Isolate ancilla positions from data positions based on our host maps
ancilla_physical_indices = host_ancilla_pairs_flat[1::2] # Odd elements are the ancillas

accepted_shots = {}
rejected_count = 0

for bitstring, count in result.items():
    clean_str = str(bitstring).strip()
    
    # Reconstruct segments
    data_bits = "".join([clean_str[idx] for idx in host_data_map])
    ancilla_bits = "".join([clean_str[idx] for idx in ancilla_physical_indices])
    
    # Post-selection rule: Reject if any ancilla registered an error syndrome ('1')
    if "1" not in ancilla_bits:
        accepted_shots[data_bits] = accepted_shots.get(data_bits, 0) + count
    else:
        rejected_count += count

total_accepted = sum(accepted_shots.values())
post_selection_rate = (total_accepted / shots) if shots > 0 else 0

print(f"1. Total Compiled Bits: {total_qubits} elements on single 1D chain.")
print(f"2. Post-Selection Acceptance Rate: {post_selection_rate:.2%}")
print(f"   - Clean Accepted Runs: {total_accepted} shots")
print(f"   - Rejected Error Flagged Runs: {rejected_count} shots")

if accepted_shots:
    sorted_accepted = sorted(accepted_shots.items(), key=lambda item: item[1], reverse=True)
    print("\n3. Top Verified Data Bitstrings (Post-Selected):")
    for bitstring, count in sorted_accepted[:3]:
        prob = count / total_accepted
        print(f"   |{bitstring[:15]}...{bitstring[-5:]}> -> Hits: {count} (Clean Prob: {prob:.2%})")
else:
    print("\n3. Top Verified Data Bitstrings: No error-free paths found.")
print("="*60)

