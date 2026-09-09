# ===================================================================
# Corrected GAP Script: Coppersmith Loop vs. MINSPM Analytical Product
# ===================================================================

Z8_Dictionary := function(phase_idx)
    return E(8)^(phase_idx mod 8);
end;

Coppersmith_Scalar := function(x_bits, y_bits, n_qubits)
    local phase_sum, i, j, k_diff, dist_factor;
    phase_sum := 0;
    for i in [1..n_qubits] do
        if x_bits[i] = 1 and y_bits[i] = 1 then
            phase_sum := phase_sum + 4;
        fi;
        for j in [1..(i-1)] do
            k_diff := i - j;
            if k_diff <= 3 then
                if x_bits[j] = 1 and y_bits[i] = 1 then
                    dist_factor := 4 / (2^(k_diff - 1));
                    phase_sum := phase_sum + Int(dist_factor);
                fi;
            fi;
        od;
    od;
    return Z8_Dictionary(phase_sum);
end;

# 1. Classical Coppersmith Expansion (O(2^t) Brute-Force Summation)
Compute_Coppersmith_Superposition := function(n_qubits, t_gates, y_vec)
    local num_terms, amp_sum, idx, x_vec, i, start_time, end_time, elapsed;
    num_terms := 2^t_gates;
    amp_sum := 0 * E(8);
    
    start_time := Runtime();
    for idx in [1..num_terms] do
        x_vec := List([1..n_qubits], k -> 0);
        for i in [1..t_gates] do
            x_vec[i] := QuoInt(idx - 1, 2^(i-1)) mod 2;
        od;
        amp_sum := amp_sum + Coppersmith_Scalar(x_vec, y_vec, n_qubits);
    od;
    end_time := Runtime();
    
    elapsed := (end_time - start_time) / 1000.0;
    return [amp_sum, elapsed];
end;

# 2. MINSPM Galois Ring Analytical Factorization (O(t * n) Product Formula)
Compute_MINSPM_Algebraic := function(n_qubits, t_gates, y_vec)
    local start_time, end_time, elapsed, total_amp, j, i, k_diff, dist_factor, coeff_j;
    
    start_time := Runtime();
    
    # MINSPM replaces the 2^t loop by exploiting the linear factorization 
    # of the quadratic phase form over the Galois ring GR(2^3, 1):
    # Sum_{x} E(8)^phase(x) = Prod_{j=1}^t (1 + E(8)^coeff_j)
    total_amp := 1 * E(8)^0;
    
    for j in [1..t_gates] do
        coeff_j := 0;
        if y_vec[j] = 1 then
            coeff_j := coeff_j + 4;
        fi;
        
        # Accumulate cross-coupling terms for variable j across all i > j
        for i in [(j+1)..n_qubits] do
            k_diff := i - j;
            if k_diff <= 3 then
                if y_vec[i] =  1 then
                    dist_factor := Int(4 / (2^(k_diff - 1)));
                    coeff_j := coeff_j + dist_factor;
                fi;
            fi;
        od;
        
        # Multiply independent factor for dimension j
        total_amp := total_amp * (1 + Z8_Dictionary(coeff_j));
    od;
    
    end_time := Runtime();
    elapsed := (end_time - start_time) / 1000.0;
    return [total_amp, elapsed];
end;

# Execution and Runtime/Amplitude Comparison
Execute_Benchmark := function()
    local n, t, y_vec, cop_res, min_res;
    n := 64;
    t := 20; # 1,048,576 terms
    y_vec := List([1..n], i -> (i mod 2));
    
    Print("=== Performance & Amplitude Comparison (n = ", n, ", t = ", t, " -> ", 2^t, " terms) ===\n");
    
    cop_res := Compute_Coppersmith_Superposition(n, t, y_vec);
    min_res := Compute_MINSPM_Algebraic(n, t, y_vec);
    
    Print("\n[Method 1: Classical Coppersmith Expansion]\n");
    Print("  - Computed Amplitude : ", cop_res[1], "\n");
    Print("  - Execution Runtime  : ", cop_res[2], " seconds\n");
    
    Print("\n[Method 2: MINSPM Galois Analytical Product]\n");
    Print("  - Computed Amplitude : ", min_res[1], "\n");
    Print("  - Execution Runtime  : ", min_res[2], " seconds\n");
    
    Print("\nExact Algebraic Equivalence Match : ", cop_res[1] = min_res[1], "\n");
end;

Execute_Benchmark();
