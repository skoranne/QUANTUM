# ===================================================================
# LARGE-SCALE TOPOLOGY BENCHMARK SUITE (n = 100000, t = 50000, K = 6)
# ===================================================================

Generalized_Z_Dictionary := function(phase_idx, K)
    local P;
    P := 2^K;
    return E(P)^(phase_idx mod P);
end;

RunTopologyBenchmark := function(topo_name, y_vec, n, t, K)
    local P, t_start, t_coeff_class, t_prod_dict, t_prod_naive, j, i, k_diff, dist_factor, coeff_j,rec_dict, names, key, u, count, dict_amp, naive_amp, naive_time;
    
    P := 2^K;
    Print("\n-----------------------------------------------------------------\n");
    Print(" Evaluating Topology: ", topo_name, " (n = ", n, ", t = ", t, ")\n");
    Print("-----------------------------------------------------------------\n");
    
    # -----------------------------------------------------------------
    # Phase 1 & 2: Coefficient Construction & Classification (MINSPM Dictionary)
    # -----------------------------------------------------------------
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
    
    # -----------------------------------------------------------------
    # Phase 3A: Dictionary Product Evaluation (O(u) contractions)
    # -----------------------------------------------------------------
    t_start := Runtime();
    dict_amp := 1 * E(P)^0;
    for key in names do
        dict_amp := dict_amp * ((1 + Generalized_Z_Dictionary(Int(key), K))^rec_dict.(key));
    od;
    t_prod_dict := (Runtime() - t_start) / 1000.0;
    
    # -----------------------------------------------------------------
    # Phase 3B: Naive Product Evaluation (O(t) contractions)
    # -----------------------------------------------------------------
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
        naive_amp := naive_amp * (1 + Generalized_Z_Dictionary(coeff_j mod P, K));
    od;
    naive_time := (Runtime() - t_start) / 1000.0;
    
    Print("  - Unique Phase Classes (u)    : ", u, " (out of t = ", t, ")\n");
    Print("  - Compression Ratio (t / u)   : ", Float(t) / Float(u), "\n");
    Print("  - Coefficient + Class Time    : ", t_coeff_class, " seconds\n");
    Print("  - MINSPM Product Eval (O(u))  : ", t_prod_dict, " seconds\n");
    Print("  - Naive Product Eval (O(t))   : ", naive_time, " seconds\n");
    Print("  - Amplitude Equivalence Match : ", dict_amp = naive_amp, "\n");
end;

ExecuteLargeScaleTests := function()
    local n, t, K, y_periodic, y_sparse, y_block;
    n := 100000;
    t := 50000;
    K := 6;
    
    # Topology 1: Highly Periodic (Low u)
    y_periodic := List([1..n], i -> (i mod 2));
    RunTopologyBenchmark("Periodic Alternating", y_periodic, n, t, K);
    
    # Topology 2: Pseudo-Random Sparse (Moderate u)
    y_sparse := List([1..n], i -> When((i mod 7) = 0, 1, 0));
    RunTopologyBenchmark("Sparse Periodic (mod 7)", y_sparse, n, t, K);
    
    # Topology 3: Block-Structured (High u)
    y_block := List([1..n], i -> When((i mod 3) < 2, 1, 0));
    RunTopologyBenchmark("Block Alternating", y_block, n, t, K);
end;

ExecuteLargeScaleTests();
