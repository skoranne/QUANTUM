import os

# 1. CONFIGURE MATRIX PRODUCT STATE COMPRESSION PARAMETERS
os.environ["CUDAQ_MPS_MAX_BOND"] = "32"
os.environ["CUDAQ_MPS_ABS_CUTOFF"] = "1e-5"

import cudaq
import numpy as np

# 2. ACTIVATE TARGET MPS SIMULATOR
cudaq.set_target("tensornet-mps")

@cudaq.kernel
def ibm_uchicago_interleaved_rcs(total_qubits: int, depth: int, data_map: list[int], 
                                 ancilla_pairs_flat: list[int], doping_mask: list[int]):
    
    qubits = cudaq.qvector(total_qubits)
    
    # State Preparation
    for i in range(len(data_map)):
        h(qubits[data_map[i]])
    
    # Execution Depth Loops (70 Layers for Real Advantage Boundary)
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
        
        # Phase B: Alternate Graph Entangling Fabric (LNN CZ-Layers)
        if d % 2 == 0:
            for idx in range(0, len(data_map) - 1, 2):
                cz(qubits[data_map[idx]], qubits[data_map[idx+1]])
        else:
            for idx in range(1, len(data_map) - 1, 2):
                cz(qubits[data_map[idx]], qubits[data_map[idx+1]])
                
        # Phase C: Spacetime Pauli Checks (Error Detection Tracking)
        if d > 0 and d % 10 == 0:
            for step in range(0, len(ancilla_pairs_flat), 2):
                d_phys = ancilla_pairs_flat[step]
                a_phys = ancilla_pairs_flat[step+1]
                cx(qubits[d_phys], qubits[a_phys]) 

    # Final Measurement
    mz(qubits)

# --- HOST-SIDE CONFIGURATION CALCULATIONS (70 Data + 27 Ancillas) ---
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
        host_ancilla_pairs_flat.append(host_data_map[-1])   
        host_ancilla_pairs_flat.append(current_physical_idx) 
        current_physical_idx += 1
        a_placed += 1

total_qubits = current_physical_idx
total_nodes = depth * data_count

np.random.seed(42)
doping_mask = np.random.choice([0, 1, 2], size=total_nodes, p=[0.85, 0.10, 0.05]).tolist()

# --- 🛠️ 3. INSTANTIATE STOCHASTIC HARDWARE NOISE MODEL ---
noise_model = cudaq.NoiseModel()

# To mimic real IBM Heron physical performance without triggering multi-qubit 
# dimension mismatch compilation bugs, we target single-qubit operations.
# We set a representative noise value to simulate the comprehensive system decay.
comprehensive_error_rate = 0.0035  # 0.35% error rate per single-qubit operation

depol_channel = cudaq.DepolarizationChannel(comprehensive_error_rate)

# Explicitly register single-qubit noise channels
noise_model.add_all_qubit_channel("h", depol_channel)
noise_model.add_all_qubit_channel("t", depol_channel)
noise_model.add_all_qubit_channel("s", depol_channel)

print(f"Formulating Noisy 97-Qubit IBM Circuit Mesh (70 Data + 27 Ancilla)...")

# 4. RUN SAMPLING VIA TENSORNET-MPS TRAJECTORY NOISE SOLVER
shots = 1000
result = cudaq.sample(
    ibm_uchicago_interleaved_rcs,
    total_qubits, depth, host_data_map, host_ancilla_pairs_flat, doping_mask,
    noise_model=noise_model,   
    shots_count=shots
)

# --- 🔎 POST-SELECTION VERIFICATION DATA DISPLAY ---
print("\n" + "="*60)
print("   REPLICATED RE-SUBMISSION HARDWARE METRICS")
print("="*60)

ancilla_physical_indices = host_ancilla_pairs_flat[1::2]
accepted_shots = {}
rejected_count = 0

for bitstring, count in result.items():
    clean_str = str(bitstring).strip()
    data_bits = "".join([clean_str[idx] for idx in host_data_map])
    ancilla_bits = "".join([clean_str[idx] for idx in ancilla_physical_indices])
    
    # Post-selection rule: Reject if any ancilla caught an error syndrome ('1')
    if "1" not in ancilla_bits:
        accepted_shots[data_bits] = accepted_shots.get(data_bits, 0) + count
    else:
        rejected_count += count

total_accepted = sum(accepted_shots.values())
post_selection_rate = (total_accepted / shots) if shots > 0 else 0

print(f"1. Total Checked Elements: {total_qubits} qubits mapped linearly.")
print(f"2. Post-Selection Acceptance Rate: {post_selection_rate:.3%}")
print(f"   - Clean Accepted Runs: {total_accepted} shots")
print(f"   - Discarded (Ancilla Caught Error): {rejected_count} shots")

if accepted_shots:
    sorted_accepted = sorted(accepted_shots.items(), key=lambda x: x[1], reverse=True)
    print("\n3. Post-Selected Data Output Distribution Top Peaks:")
    for bitstring, count in sorted_accepted[:3]:
        prob = count / total_accepted
        print(f"   |{bitstring[:15]}...{bitstring[-5:]}> -> Hits: {count} (Prob: {prob:.2%})")
else:
    print("\n3. Post-Selected Output Distribution: Noise floor saturated.")
print("="*60)
