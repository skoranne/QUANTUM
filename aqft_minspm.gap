# ===================================================================
# MINSPM DICTIONARY vs. NAIIVE PRODUCT BENCHMARK SUITE
# Quantifying the advantage of sparse coefficient accumulation
# ===================================================================

Generalized_Z_Dictionary := function(phase_idx, K)
    local P;
    P := 2^K;
    return E(P)^(phase_idx mod P);
end;

# 1. Brute-force Coppersmith Superposition Expansion O(2^t)
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

Compute_Coppersmith_BruteForce := function(n_qubits, t_gates, y_vec, K)
    local num_terms, amp_sum, idx, x_vec, i, t_start;
    num_terms := 2^t_gates;
    amp_sum := 0 * E(2^K);
    t_start := Runtime();
    for idx in [1..num_terms] do
        x_vec := List([1..n_qubits], k -> 0);
        for i in [1..t_gates] do
            x_vec[i] := QuoInt(idx - 1, 2^(i-1)) mod 2;
        od;
        amp_sum := amp_sum + Coppersmith_Scalar_Gen(x_vec, y_vec, n_qubits, K);
    od;
    return [amp_sum, (Runtime() - t_start) / 1000.0];
end;

# 2. Naive Sequential Algebraic Factorization O(t)
Compute_Naive_Product := function(n_qubits, t_gates, y_vec, K)
    local total_amp, j, i, k_diff, dist_factor, coeff_j, P, t_start;
    P := 2^K;
    t_start := Runtime();
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
    return [total_amp, (Runtime() - t_start) / 1000.0];
end;

# 3. MINSPM Sparse Dictionary Factorization with Multiplicity Accumulation
Compute_MINSPM_Dictionary := function(n_qubits, t_gates, y_vec, K)
    local P, rec_dict, j, i, k_diff, dist_factor, coeff_j, names, key, u, count, total_amp, t_start;
    P := 2^K;
    t_start := Runtime();
    
    # Phase 1: Build sparse frequency dictionary of unique coefficients
    rec_dict := rec();
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
        
        key := String(coeff_j mod P);
        if IsBound(rec_dict.(key)) then
            rec_dict.(key) := rec_dict.(key) + 1;
        else
            rec_dict.(key) := 1;
        fi;
    od;
    
    # Phase 2: Evaluate product over unique dictionary keys raised to multiplicities
    total_amp := 1 * E(P)^0;
    names := RecNames(rec_dict);
    for key in names do
        u := Int(key);
        count := rec_dict.(key);
        # Compute (1 + zeta^u)^count efficiently or via repeated multiplication
        total_amp := total_amp * ((1 + Generalized_Z_Dictionary(u, K))^count);
    od;
    
    return [total_amp, (Runtime() - t_start) / 1000.0, Length(names)];
end;

# Comparative Benchmark Execution
Run_Comparison_Benchmark := function()
    local n, t, K, y_vec, bf_res, naive_res, dict_res;
    
    n := 64;
    t := 20; # 1,048,576 terms
    K := 6;  # Ring Z_64
    
    y_vec := [
      1,0,1,1,0,0,1,0,
      1,1,0,1,0,1,1,0,
      0,1,0,1,1,0,1,1,
      1,0,0,1,0,1,0,1,
      1,1,0,0,1,1,0,1,
      0,1,1,0,1,0,0,1,
      1,0,1,1,0,1,1,0,
      0,0,1,0,1,1,0,1
    ];
    
    Print("\n=================================================================\n");
    Print(" METHOD COMPARISON: BRUTE-FORCE vs. NAIIVE vs. MINSPM DICTIONARY\n");
    Print("=================================================================\n");
    Print("Parameters: n = ", n, ", t = ", t, " (2^t = ", 2^t, " terms), K = ", K, "\n");
    
    # 1. MINSPM Dictionary
    dict_res := Compute_MINSPM_Dictionary(n, t, y_vec, K);
    Print("\n[Method A: MINSPM Dictionary Factorization]\n");
    Print("  - Amplitude Result  : ", dict_res[1], "\n");
    Print("  - Unique Classes (u): ", dict_res[3], " out of ", t, " variables\n");
    Print("  - Execution Runtime : ", dict_res[2], " seconds\n");
    
    # 2. Naive Product
    naive_res := Compute_Naive_Product(n, t, y_vec, K);
    Print("\n[Method B: Naive Sequential Product]\n");
    Print("  - Amplitude Result  : ", naive_res[1], "\n");
    Print("  - Execution Runtime : ", naive_res[2], " seconds\n");
    
    # 3. Brute-Force Coppersmith (Only run for smaller t if needed, but works for 2^20)
    Print("\n[Method C: Brute-Force Coppersmith Expansion O(2^t)]\n");
    Print("  - Running 1,048,576 state evaluations...\n");
    bf_res := Compute_Coppersmith_BruteForce(n, t, y_vec, K);
    Print("  - Amplitude Result  : ", bf_res[1], "\n");
    Print("  - Execution Runtime : ", bf_res[2], " seconds\n");
    
    Print("\n-----------------------------------------------------------------\n");
    Print(" Exact Equivalence Match (A = B = C): ", 
          (dict_res[1] = naive_res[1]) and (naive_res[1] = bf_res[1]), "\n");
    Print("=================================================================\n");
end;

Run_Comparison_Benchmark();
