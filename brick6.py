import cudaq
import numpy as np

@cudaq.kernel
def ibm_advantage_kernel(qubit_count: int, layers: int, angles: list[float]):
    qubits = cudaq.qvector(qubit_count)
    
    # 1. State Preparation Layer
    h(qubits) 
    
    # 2. Iterate through structured brickwork layers
    for layer in range(layers):
        # Non-Clifford Doping Layer
        for i in range(qubit_count):
            idx = layer * qubit_count + i
            rz(angles[idx], qubits[i])
            
        # Entangling Fabric (Brickwork Lattice)
        # Even qubit couplings
        for i in range(0, qubit_count - 1, 2):
            cx(qubits[i], qubits[i+1])
            
        # Odd qubit couplings
        for i in range(1, qubit_count - 1, 2):
            cx(qubits[i], qubits[i+1])
            
    # NOTE: We remove the internal mz(qubits) call because 
    # cudaq.observe handles measurement externally via the spin operator.

# --- Execution Framework ---
qubit_count = 6  # Scaled down to 6 qubits for easier observation tracking
layers = 2
total_angles = layers * qubit_count
np.random.seed(42)
random_angles = np.random.uniform(0, 2 * np.pi, total_angles).tolist()

# Tensornet works perfectly fine for small counts too
cudaq.set_target("tensornet") 

print(f"Observing expectation value for a {qubit_count}-qubit system...")

# Define the spin operator: Pauli-Z on Qubit 0 (cudaq.spin.z(0))
# This mathematically defines what observable we want to calculate
hamiltonian = cudaq.spin.z(0)

# Run observation instead of sampling
observation_result = cudaq.observe(
    ibm_advantage_kernel, 
    hamiltonian,
    qubit_count, 
    layers, 
    random_angles
)

# Extract and display the expected value
exp_val = observation_result.expectation()
print(f"Expectation Value <Z0>: {exp_val:.6f}")
