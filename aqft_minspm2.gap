# ===================================================================
# LARGE-SCALE MINSPM COMPRESSION BENCHMARK: t = 500, n = 1000
# Demonstrates compression advantage when unique classes u << t
# ===================================================================

Generalized_Z_Dictionary := function(phase_idx, K)
    local P;
    P := 2^K;
    return E(P)^(phase_idx mod P);
end;

RunLargeScaleBenchmark := function()
    local n, t, K, y_vec, P, j, i, k_diff, dist_factor, coeff_j,t_start, naive_amp, dict_amp, rec_dict, names, key, u, count,naive_time, dict_time, unique_classes;

    n := 1000;
    t := 500;
    K := 6;
    P := 2^K;

    # Construct a structured periodic y_vec to force a small number of unique coefficient classes u relative to t
    y_vec := List([1..n], i -> (i mod 2));

    Print("\n=================================================================\n");
    Print(" LARGE-SCALE MINSPM COMPRESSION BENCHMARK (t = ", t, ", n = ", n, ")\n");
    Print("=================================================================\n");

    # ---------------------------------------------------------
    # 1. Naive Sequential Product (t individual factor evaluations)
    # ---------------------------------------------------------
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
        naive_amp := naive_amp * (1 + Generalized_Z_Dictionary(coeff_j, K));
    od;
    naive_time := (Runtime() - t_start) / 1000.0;

    # ---------------------------------------------------------
    # 2. MINSPM Dictionary Factorization (u unique classes + multiplicities)
    # ---------------------------------------------------------
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

    dict_amp := 1 * E(P)^0;
    names := RecNames(rec_dict);
    unique_classes := Length(names);

    for key in names do
        u := Int(key);
        count := rec_dict.(key);
        # Compute unique factor raised to its multiplicity: (1 + zeta^u)^count
        dict_amp := dict_amp * ((1 + Generalized_Z_Dictionary(u, K))^count);
    od;
    dict_time := (Runtime() - t_start) / 1000.0;

    Print("  - Total Superposition Variables (t) : ", t, "\n");
    Print("  - Unique Phase Classes (u)          : ", unique_classes, "\n");
    Print("  - Compression Ratio (t / u)         : ", Float(t) / Float(unique_classes), "\n");
    Print("  - Naive Product Runtime             : ", naive_time, " seconds\n");
    Print("  - MINSPM Dictionary Runtime         : ", dict_time, " seconds\n");
    Print("  - Exact Algebraic Match             : ", naive_amp = dict_amp, "\n");
    Print("=================================================================\n");
end;

RunLargeScaleBenchmark();
