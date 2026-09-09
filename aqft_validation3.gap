# ===================================================================
# EXHAUSTIVE THEOREM VERIFICATION SUITE
# Tests all 2^14 (16,384) possible y-vectors for K = 4, t = 10, n = 14
# ===================================================================

Generalized_Z_Dictionary := function(phase_idx, K)
    local P;
    P := 2^K;
    return E(P)^(phase_idx mod P);
end;

Coppersmith_Scalar_Gen := function(x_bits, y_bits, n_qubits, K)
    local phase_sum, i, j, k_diff, dist_factor, P;
    P := 2^K;
    phase_sum := 0;
    
    for i in [1..n_qubits] do
        if x_bits[i] = 1 and y_bits[i] = 1 then
            phase_sum := phase_sum + QuoInt(P, 2);
        fi;
        for j in [1..(i-1)] do
            k_diff := i - j;
            if k_diff <= K then
                if x_bits[j] = 1 and y_bits[i] = 1 then
                    dist_factor := QuoInt(P, 2^k_diff);
                    phase_sum := phase_sum + dist_factor;
                fi;
            fi;
        od;
    od;
    return Generalized_Z_Dictionary(phase_sum, K);
end;

# Brute-force superposition expansion over 2^t terms
Compute_Coppersmith_Superposition_Gen := function(n_qubits, t_gates, y_vec, K)
    local num_terms, amp_sum, idx, x_vec, i;
    num_terms := 2^t_gates;
    amp_sum := 0 * E(2^K);
    
    for idx in [1..num_terms] do
        x_vec := List([1..n_qubits], k -> 0);
        for i in [1..t_gates] do
            x_vec[i] := QuoInt(idx - 1, 2^(i-1)) mod 2;
        od;
        amp_sum := amp_sum + Coppersmith_Scalar_Gen(x_vec, y_vec, n_qubits, K);
    od;
    return amp_sum;
end;

# MINSPM analytical product over t variables
Compute_MINSPM_Algebraic_Gen := function(n_qubits, t_gates, y_vec, K)
    local total_amp, j, i, k_diff, dist_factor, coeff_j, P;
    P := 2^K;
    total_amp := 1 * E(P)^0;
    
    for j in [1..t_gates] do
        coeff_j := 0;
        if y_vec[j] = 1 then
            coeff_j := coeff_j + QuoInt(P, 2);
        fi;
        
        for i in [(j+1)..n_qubits] do
            k_diff := i - j;
            if k_diff <= K then
                if y_vec[i] = 1 then
                    dist_factor := QuoInt(P, 2^k_diff);
                    coeff_j := coeff_j + dist_factor;
                fi;
            fi;
        od;
        
        total_amp := total_amp * (1 + Generalized_Z_Dictionary(coeff_j, K));
    od;
    
    return total_amp;
end;

# Exhaustive validation loop across the entire input space
RunExhaustiveTest := function()
    local n, t, K, num_y, y_idx, y_vec, i, cop_val, min_val, start_t, mismatches;
    
    n := 14;
    t := 10;
    K := 4;
    num_y := 2^n;
    
    Print("\n=================================================================\n");
    Print(" EXHAUSTIVE THEOREM VERIFICATION (K = ", K, ", t = ", t, ", n = ", n, ")\n");
    Print(" Enumerating all ", num_y, " possible y-vectors...\n");
    Print("=================================================================\n");
    
    start_t := Runtime();
    mismatches := 0;
    
    for y_idx in [1..num_y] do
        # Decode index into binary vector of length n
        y_vec := List([1..n], k -> 0);
        for i in [1..n] do
            y_vec[i] := QuoInt(y_idx - 1, 2^(i-1)) mod 2;
        od;
        
        cop_val := Compute_Coppersmith_Superposition_Gen(n, t, y_vec, K);
        min_val := Compute_MINSPM_Algebraic_Gen(n, t, y_vec, K);
        
        if cop_val <> min_val then
            mismatches := mismatches + 1;
            Print("  [FAIL] Mismatch at y_idx = ", y_idx, " vector: ", y_vec, "\n");
        fi;
        
        if (y_idx mod 2048) = 0 then
            Print("  -> Progress: Checked ", y_idx, " / ", num_y, " vectors...\n");
        fi;
    od;
    
    Print("-----------------------------------------------------------------\n");
    Print(" Exhaustive Test Complete in ", (Runtime() - start_t) / 1000.0, " seconds.\n");
    Print(" Total Mismatches Found: ", mismatches, "\n");
    Print("=================================================================\n");
end;

RunExhaustiveTest();
