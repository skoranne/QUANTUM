# ==============================================================================
# GAP SCRIPT: NAM-SU-MASLOV OPTIMIZED COPPERSMITH VS. MINSPM GALOIS LOOKUP
# Uses Cyclic Boundary Conditions to Guarantee Non-Zero Amplitude Validation
# ==============================================================================

Z8_Dictionary := function(phase_idx)
    return E(8)^(phase_idx mod 8);
end;

# 1. Nam-Su-Maslov Optimized Coppersmith Expansion 
# Evaluates the exponential O(2^n) full state superposition space.
Compute_Optimized_Coppersmith := function(n_qubits, y_vec)
    local num_terms, amp_sum, idx, temp_idx, i, start_time, phase_sum, phi_coeffs, phase_val, y_prev1, y_prev2;
    
    num_terms := 2^n_qubits;
    amp_sum := 0 * E(8);
    
    # Pre-compute Nam-Su-Maslov sliding-window phases with cyclic boundary conditions
    phi_coeffs := List([1..n_qubits], k -> 0);
    for i in [1..n_qubits] do
        if i = 1 then
            y_prev1 := y_vec[n_qubits];
            y_prev2 := y_vec[n_qubits-1];
        elif i = 2 then
            y_prev1 := y_vec[1];
            y_prev2 := y_vec[n_qubits];
        else
            y_prev1 := y_vec[i-1];
            y_prev2 := y_vec[i-2];
        fi;

        phase_val := 0;
        if y_vec[i] = 1 then phase_val := phase_val + 4; fi;
        if y_prev1 = 1 then phase_val := phase_val + 2; fi;
        if y_prev2 = 1 then phase_val := phase_val + 1; fi;
        phi_coeffs[i] := phase_val;
    od;
    
    start_time := Runtime();
    
    for idx in [0..(num_terms-1)] do
        phase_sum := 0;
        temp_idx := idx;
        
        for i in [1..n_qubits] do
            if temp_idx mod 2 = 1 then
                phase_sum := phase_sum + phi_coeffs[i];
            fi;
            temp_idx := Int(temp_idx / 2);
        od;
        
        amp_sum := amp_sum + Z8_Dictionary(phase_sum);
    od;
    
    return [amp_sum, Runtime() - start_time];
end;

# 2. MINSPM Galois Dictionary Lookup (O(n) Factorized Ring Product)
Compute_MINSPM_Lookup := function(n_qubits, y_vec)
    local amp_prod, i, phase_val, start_time, y_prev1, y_prev2;
    
    start_time := Runtime();
    amp_prod := 1 * E(8)^0;
    
    for i in [1..n_qubits] do
        if i = 1 then
            y_prev1 := y_vec[n_qubits];
            y_prev2 := y_vec[n_qubits-1];
        elif i = 2 then
            y_prev1 := y_vec[1];
            y_prev2 := y_vec[n_qubits];
        else
            y_prev1 := y_vec[i-1];
            y_prev2 := y_vec[i-2];
        fi;

        phase_val := 0;
        if y_vec[i] = 1 then phase_val := phase_val + 4; fi;
        if y_prev1 = 1 then phase_val := phase_val + 2; fi;
        if y_prev2 = 1 then phase_val := phase_val + 1; fi;
        
        amp_prod := amp_prod * (1 + Z8_Dictionary(phase_val));
    od;
    
    return [amp_prod, Runtime() - start_time];
end;

# ==============================================================================
# BENCHMARK EXECUTION (Sweep across n = 12, 16, 20, 24, 28)
# ==============================================================================
Run_Comparison := function(n_qubits)
    local y_vec, t1, t2;
    
    # State vector with y_i = 1 for all qubits ensures non-zero constructive interference
    y_vec := List([1..n_qubits], i -> 1);
    
    Print("\n========================================================\n");
    Print(" RUNNING COMPARISON FOR n = ", n_qubits, " QUBITS\n");
    Print("========================================================\n");
    
    t1 := Compute_Optimized_Coppersmith(n_qubits, y_vec);
    t2 := Compute_MINSPM_Lookup(n_qubits, y_vec);

    Print("[Method 1: Classical Coppersmith Expansion]\n");
    Print("  - Computed Amplitude : ", t1[1], "\n");
    Print("  - Execution Runtime  : ", Float(t1[2] / 1000), " seconds\n\n");

    Print("[Method 2: MINSPM Galois Dictionary Lookup]\n");
    Print("  - Computed Amplitude : ", t2[1], "\n");
    Print("  - Execution Runtime  : ", Float(t2[2] / 1000), " seconds\n");
end;

# Execute benchmark for target qubit sizes
for n in [12, 16, 20, 24] do
    Run_Comparison(n);
od;
