# ===================================================================
# MINSPM vs COPPERSMITH: COMPREHENSIVE VALIDATION & COEFFICIENT PARITY SUITE
# Self-contained GAP 4 Program (Includes Corrected Closed-Form Parity Check)
# ===================================================================

Generalized_Z_Dictionary := function(phase_idx, K)
    local P;
    P := 2^K;
    return E(P)^(phase_idx mod P);
end;

# 1. Generalized Classical Coppersmith Scalar Phase
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

# Brute-force summation for validation baseline O(2^t)
Compute_Coppersmith_Superposition_Gen := function(n_qubits, t_gates, y_vec, K)
    local num_terms, amp_sum, idx, x_vec, i, start_t, end_t;
    num_terms := 2^t_gates;
    amp_sum := 0 * E(2^K);
    
    start_t := Runtime();
    for idx in [1..num_terms] do
        x_vec := List([1..n_qubits], k -> 0);
        for i in [1..t_gates] do
            x_vec[i] := QuoInt(idx - 1, 2^(i-1)) mod 2;
        od;
        amp_sum := amp_sum + Coppersmith_Scalar_Gen(x_vec, y_vec, n_qubits, K);
    od;
    end_t := Runtime();
    
    return [amp_sum, (end_t - start_t) / 1000.0];
end;

# 2. Generalized MINSPM Finite-Ring Analytical Factorization O(t)
Compute_MINSPM_Algebraic_Gen := function(n_qubits, t_gates, y_vec, K)
    local total_amp, j, i, k_diff, dist_factor, coeff_j, P, start_t, end_t;
    P := 2^K;
    
    start_t := Runtime();
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
    end_t := Runtime();
    
    return [total_amp, (end_t - start_t) / 1000.0];
end;

# 3. Definitive Coefficient Parity Check (Corrected Closed-Form vs Scalar Loop)
VerifyCoefficientsMatchCorrected := function()
    local n, t, y_vec, j, i, k_diff, coeff_alg, coeff_scalar_check, all_match;
    
    n := 64;
    t := 26;
    
    y_vec := [
      1,1,0,1,0,1,1,0,
      1,0,1,1,1,0,0,1,
      0,1,1,0,1,1,0,1,
      1,1,0,0,1,0,1,1,
      0,1,0,1,1,0,1,0,
      1,0,1,0,0,1,1,0,
      1,1,0,1,0,1,0,1,
      0,1,1,0,1,0,1,1
    ];
    
    Print("\n========================================================\n");
    Print(" CORRECTED COEFFICIENT PARITY CHECK (MINSPM vs SCALAR)\n");
    Print("========================================================\n");
    
    all_match := true;
    
    for j in [1..t] do
        # Corrected closed-form algebraic formula for K=3:
        # c_j = 4*y_j + 4*y_{j+1} + 2*y_{j+2} + y_{j+3} (mod 8)
        coeff_alg := 4 * y_vec[j];
        if j + 1 <= n then coeff_alg := coeff_alg + 4 * y_vec[j+1]; fi;
        if j + 2 <= n then coeff_alg := coeff_alg + 2 * y_vec[j+2]; fi;
        if j + 3 <= n then coeff_alg := coeff_alg + y_vec[j+3]; fi;
        coeff_alg := coeff_alg mod 8;
        
        # Explicit scalar loop coefficient extraction matching Coppersmith_Scalar
        coeff_scalar_check := 0;
        if y_vec[j] = 1 then
            coeff_scalar_check := coeff_scalar_check + 4;
        fi;
        for i in [(j+1)..n] do
            k_diff := i - j;
            if k_diff <= 3 then
                if y_vec[i] = 1 then
                    coeff_scalar_check := coeff_scalar_check + Int(4 / (2^(k_diff - 1)));
                fi;
            fi;
        od;
        coeff_scalar_check := coeff_scalar_check mod 8;
        
        if coeff_alg <> coeff_scalar_check then
            all_match := false;
        fi;
        
        Print("j = ", j, 
              " | Alg Formula c_j = ", coeff_alg, 
              " | Scalar Loop c_j = ", coeff_scalar_check, 
              " | Match: ", coeff_alg = coeff_scalar_check, 
              " | Factor: 1 + E(8)^", coeff_alg, "\n");
    od;
    
    Print("--------------------------------------------------------\n");
    if all_match then
        Print(" RESULT: 100% Coefficient Parity Match Confirmed Across All Indices!\n");
    else
        Print(" RESULT: Mismatch detected in coefficient mapping.\n");
    fi;
    Print("========================================================\n");
end;

# 4. Comprehensive High-K and Large-t Validation Suite
Run_Full_Validation_Suite := function()
    local n, t, K, y_vec_explicit, y_rand, cop_res, min_res;
    
    # Run coefficient parity checks first
    VerifyCoefficientsMatchCorrected();
    
    Print("\n=================================================================\n");
    Print(" MINSPM vs COPPERSMITH: HIGH-K TRUNCATION & PERFORMANCE BENCHMARK\n");
    Print("=================================================================\n");
    
    # TEST 1: Explicit 64-qubit vector (K=10, t=20)
    n := 64;
    t := 20; 
    K := 10;
    
    y_vec_explicit := [
      1,0,1,1,0,0,1,0,
      1,1,0,1,0,1,1,0,
      0,1,0,1,1,0,1,1,
      1,0,0,1,0,1,0,1,
      1,1,0,0,1,1,0,1,
      0,1,1,0,1,0,0,1,
      1,0,1,1,0,1,1,0,
      0,0,1,0,1,1,0,1
    ];
    
    Print("\n--> TEST 1: Explicit Test Vector (n = 64, t = 20, K = 10, Ring Z_1024)\n");
    Print("    Evaluating 2^20 (1,048,576) explicit superposition terms...\n");
    
    cop_res := Compute_Coppersmith_Superposition_Gen(n, t, y_vec_explicit, K);
    min_res := Compute_MINSPM_Algebraic_Gen(n, t, y_vec_explicit, K);
    
    if cop_res[1] = min_res[1] then
        Print("    [PASS] Exact Algebraic Match Confirmed.\n");
    else
        Print("    [FAIL] Mismatch detected!\n");
    fi;
    
    Print("    - Classical Coppersmith O(2^t) Runtime : ", cop_res[2], " seconds\n");
    Print("    - MINSPM Algebraic O(t) Runtime        : ", min_res[2], " seconds\n");
    
    # TEST 2: Randomized Vector (K=12, t=20)
    K := 12;
    y_rand := List([1..n], i -> Random([0, 1]));
    
    Print("\n--> TEST 2: Randomized Test Vector (n = 64, t = 20, K = 12, Ring Z_4096)\n");
    Print("    Evaluating 2^20 (1,048,576) explicit superposition terms...\n");
    
    cop_res := Compute_Coppersmith_Superposition_Gen(n, t, y_rand, K);
    min_res := Compute_MINSPM_Algebraic_Gen(n, t, y_rand, K);
    
    if cop_res[1] = min_res[1] then
        Print("    [PASS] Exact Algebraic Match Confirmed.\n");
    else
        Print("    [FAIL] Mismatch detected!\n");
    fi;
    
    Print("    - Classical Coppersmith O(2^t) Runtime : ", cop_res[2], " seconds\n");
    Print("    - MINSPM Algebraic O(t) Runtime        : ", min_res[2], " seconds\n");
    Print("\n=================================================================\n");
end;

# Execute everything
Run_Full_Validation_Suite();
