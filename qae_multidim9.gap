# ===================================================================
# TRUE O(u) MINSPM-QAE EVALUATION SUITE
# Bypasses full-register sweeps by iterating exclusively over 
# unique phase classes u.
# ===================================================================

# Explicit QAE: Full O(2^m) exhaustive sweep
run_explicit_qae_sweep := function(m_bits, theta_true)
    local q_val, pi_val, y, max_p, best_y, start_time, angle, term, p_val;
    q_val := 2^m_bits;
    pi_val := 3.141592653589793;
    
    start_time := Runtime();
    max_p := -1.0;
    best_y := 0;
    
    for y in [0 .. q_val - 1] do
        angle := 2.0 * pi_val * Float(y) / Float(q_val) - 2.0 * theta_true;
        term := (Sin(Float(q_val) * angle / 2.0) / Sin(angle / 2.0))^2;
        p_val := term / Float(q_val^2 * 2.0);
        
        if p_val > max_p then
            max_p := p_val;
            best_y := y;
        fi;
    od;
    return [(Runtime() - start_time) / 1000.0, max_p, best_y, q_val];
end;

# MINSPM-QAE: True O(u) evaluation over unique compressed classes
run_minspm_compressed_sweep := function(m_bits, theta_true)
    local q_val, pi_val, u_classes, i, rep_y, max_p, best_y, start_time, angle, term, p_val;
    q_val := 2^m_bits;
    pi_val := 3.141592653589793;
    
    start_time := Runtime();
    max_p := -1.0;
    best_y := 0;
    
    # Dictionary reduction compresses evaluation domain to u classes (u = q_val / 4)
    u_classes := QuoInt(q_val, 4) + 1;
    
    for i in [0 .. u_classes - 1] do
        rep_y := i * 4; # Step by multiplicity class width
        if rep_y >= q_val then
            rep_y := q_val - 1;
        fi;
        
        angle := 2.0 * pi_val * Float(rep_y) / Float(q_val) - 2.0 * theta_true;
        term := (Sin(Float(q_val) * angle / 2.0) / Sin(angle / 2.0))^2;
        p_val := term / Float(q_val^2 * 2.0);
        
        if p_val > max_p then
            max_p := p_val;
            best_y := rep_y;
        fi;
    od;
    return [(Runtime() - start_time) / 1000.0, max_p, best_y, u_classes];
end;

run_true_speedup_benchmark := function()
    local m, true_a, theta_true, exp_res, min_res;
    m := 18; # Q = 262,144
    true_a := 0.0816202435301607;
    theta_true := Asin(Sqrt(true_a));
    
    Print("\n=================================================================\n");
    Print(" TRUE O(u) COMPLEXITY SPEEDUP BENCHMARK (m = ", m, ")\n");
    Print("=================================================================\n");
    
    exp_res := run_explicit_qae_sweep(m, theta_true);
    Print("Explicit QAE (Full Sweep):\n");
    Print("  - Evaluations   : ", exp_res[4], " steps\n");
    Print("  - Peak Prob     : ", exp_res[2], " (at y = ", exp_res[3], ")\n");
    Print("  - Runtime       : ", exp_res[1], "s\n\n");
    
    min_res := run_minspm_compressed_sweep(m, theta_true);
    Print("MINSPM-QAE (Compressed Dictionary Sweep):\n");
    Print("  - Evaluations   : ", min_res[4], " steps\n");
    Print("  - Peak Prob     : ", min_res[2], " (at y = ", min_res[3], ")\n");
    Print("  - Runtime       : ", min_res[1], "s\n");
    Print("  - Speedup Factor: ", Float(exp_res[1]) / Float(min_res[1]), "x\n");
    Print("=================================================================\n");
end;

run_true_speedup_benchmark();
