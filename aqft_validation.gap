# ===================================================================
# Generalized Coppersmith vs. MINSPM AQFT Validation Suite
# Includes High-K (K=10) and Large-t (1,048,576 terms) Benchmarks
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

# High-K Explicit and Randomized Test Harness
Run_High_K_Validation_Suite := function()
    local n, t, K, y_vec_explicit, y_rand, cop_res, min_res;
    
    Print("\n=================================================================\n");
    Print(" MINSPM vs COPPERSMITH: HIGH-K TRUNCATION & PERFORMANCE BENCHMARK\n");
    Print("=================================================================\n");
    
    # ---------------------------------------------------------
    # TEST 1: Explicit 64-qubit vector (K=10, t=20)
    # ---------------------------------------------------------
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
    
    # ---------------------------------------------------------
    # TEST 2: Randomized Vector (K=12, t=20)
    # ---------------------------------------------------------
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

Run_High_K_Validation_Suite();
