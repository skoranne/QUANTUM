# ===================================================================
# Corrected GAP Script: Galois Entropy Bound for AQFT vs Full QFT
# ===================================================================

# 1. AQFT (Truncated at k=3)
# The finest angle is 2*pi / 2^3 = pi/4, represented by the 8th root of unity E(8).
Print("--- AQFT (k=3 Truncation) ---\n");
R3_Phase := E(8);

# Calculate the multiplicative order (Dictionary Depth / Orbit Size)
orbit_size_aqft := Order(R3_Phase);
Print("Phase Dictionary Size (Orbit): ", orbit_size_aqft, "\n");

# Calculate the Galois Extension Degree (VRAM Memory Multiplier)
field_aqft := CF(8); # Cyclotomic Field of order 8
Print("Galois Extension Degree: ", Dimension(field_aqft), "\n\n");


# 2. FULL QFT (e.g., n = 10 qubits)
# The finest angle requires the 1024th root of unity, E(1024).
Print("--- EXACT QFT (n=10 qubits) ---\n");
R10_Phase := E(1024);

# Calculate the multiplicative order (Dictionary Depth / Orbit Size)
orbit_size_qft := Order(R10_Phase);
Print("Phase Dictionary Size (Orbit): ", orbit_size_qft, "\n");

# Calculate the Galois Extension Degree (VRAM Memory Multiplier)
field_qft := CF(1024); # Cyclotomic Field of order 1024
Print("Galois Extension Degree: ", Dimension(field_qft), "\n");
