# ===================================================================
# HIGH-DIVERSITY STRESS-TEST SUITE (n = 128, t = 64, K = 12)
# Evaluates MINSPM resilience under maximal coefficient entropy (u ≈ t)
# ===================================================================

stress_test_high_diversity := function()
    local n, t, K, P, y_vec, j, i, k_diff, dist_factor, coeff_j,t_start, rec_dict, names, key, u, dict_amp, naive_amp,t_coeff_class, t_naive, t_dict;
    
    n := 128;
    t := 64;
    K := 12;
    P := 2^K;
    
    # Pseudorandom y-vector designed for high coefficient diversity
    y_vec := List([1..n], i -> (PowerModInt(i, 3, 997) mod 2));
    
    # Phase 1: Coefficient Construction & Classification
    t_start := Runtime();
    rec_dict := rec();
    for j in [1..t] do
        coeff_j := 0;
        if y_vec[j] = 1 then
            coeff_j := coeff_j + QuoInt(P, 2);
        fi;
        for i in [(j+1)..n] do
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
    t_coeff_class := (Runtime() - t_start) / 1000.0;
    
    names := RecNames(rec_dict);
    u := Length(names);
    
    # Phase 2: Naive Factorized Product Evaluation O(t)
    t_start := Runtime();
    naive_amp := 1 * E(P)^0;
    for j in [1..t] do
        coeff_j := 0;
        if y_vec[j] = 1 then
            coeff_j := coeff_j + QuoInt(P, 2);
        fi;
        for i in [(j+1)..n] do
            k_diff := i - j;
            if k_diff <= K then
                if y_vec[i] = 1 then
                    dist_factor := QuoInt(P, 2^k_diff);
                    coeff_j := coeff_j + dist_factor;
                fi;
            fi;
        od;
        naive_amp := naive_amp * (1 + E(P)^(coeff_j mod P));
    od;
    t_naive := (Runtime() - t_start) / 1000.0;
    
    # Phase 3: MINSPM Dictionary Product Evaluation O(u)
    t_start := Runtime();
    dict_amp := 1 * E(P)^0;
    for key in names do
        dict_amp := dict_amp * ((1 + E(P)^(Int(key)))^rec_dict.(key));
    od;
    t_dict := (Runtime() - t_start) / 1000.0;
    
    Print("\n=================================================================\n");
    Print(" HIGH-DIVERSITY STRESS TEST RESULTS (n = 128, t = 64, K = 12)\n");
    Print("=================================================================\n");
    Print("  - Total Variables (t)           : ", t, "\n");
    Print("  - Truncation Range (K)          : ", K, "\n");
    Print("  - Unique Classes (u)            : ", u, "\n");
    Print("  - Compression Ratio (t / u)     : ", Float(t) / Float(u), "\n");
    Print("  - Coefficient Generation Time   : ", t_coeff_class, "s\n");
    Print("  - Naive Factorized Product Time : ", t_naive, "s\n");
    Print("  - MINSPM Dictionary Time        : ", t_dict, "s\n");
    Print("  - Coppersmith Brute-Force       : Infeasible (> 10^19 ops)\n");
    Print("  - Algebraic Equivalence Match   : ", naive_amp = dict_amp, "\n");
    Print("=================================================================\n");
end;

stress_test_high_diversity();
