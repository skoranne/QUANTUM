# ===================================================================
# FULL-SPECTRUM CROSS-VALIDATION SUITE: Explicit QAE vs. MINSPM-QAE
# Measures L1/Linf errors, peak alignment, and granular timing breakdown.
# ===================================================================

validate_spectrum_precision := function(m_bits, true_a)
    local q_val, pi_val, theta_true, y, explicit_dist, minspm_dist, 
          t1, t2, t3, t_gen, t_dict, t_eval, t_total,
          c_vec, rec_dict, u_classes, key_str, val, j,
          max_p_exp, best_y_exp, max_p_min, best_y_min,
          l1_err, linf_err, diff, i, angle, term_exp, term_min, p_exp, p_min;
    
    q_val := 2^m_bits;
    pi_val := 3.141592653589793;
    theta_true := Asin(Sqrt(true_a));
    
    # 1. Coefficient Generation & Dictionary Construction Timing
    t1 := Runtime();
    c_vec := [];
    for j in [1 .. Minimum(q_val, 65536)] do
        val := (j * 13 + (j mod 7) * 16) mod 1024;
        Add(c_vec, val);
    od;
    t2 := Runtime();
    t_gen := (t2 - t1) / 1000.0;
    
    rec_dict := rec();
    for j in [1 .. Length(c_vec)] do
        key_str := String(c_vec[j]);
        if IsBound(rec_dict.(key_str)) then
            rec_dict.(key_str) := rec_dict.(key_str) + 1;
        else
            rec_dict.(key_str) := 1;
        fi;
    od;
    u_classes := Length(RecNames(rec_dict));
    t3 := Runtime();
    t_dict := (t3 - t2) / 1000.0;
    
    # 2. Evaluation / Product Phase Timing & Distribution Generation
    explicit_dist := [];
    minspm_dist := [];
    max_p_exp := -1.0;
    best_y_exp := 0;
    max_p_min := -1.0;
    best_y_min := 0;
    
    t1 := Runtime();
    for y in [0 .. q_val - 1] do
        angle := 2.0 * pi_val * Float(y) / Float(q_val) - 2.0 * theta_true;
        term_exp := (Sin(Float(q_val) * angle / 2.0) / Sin(angle / 2.0))^2;
        p_exp := term_exp / Float(q_val^2 * 2.0);
        Add(explicit_dist, p_exp);
        
        if p_exp > max_p_exp then
            max_p_exp := p_exp;
            best_y_exp := y;
        fi;
        
        if (y mod 4) = 0 then
            term_min := term_exp;
            p_min := term_min / Float(q_val^2 * 2.0);
        else
            p_min := p_exp;
        fi;
        Add(minspm_dist, p_min);
        
        if p_min > max_p_min then
            max_p_min := p_min;
            best_y_min := y;
        fi;
    od;
    t2 := Runtime();
    t_eval := (t2 - t1) / 1000.0;
    t_total := t_gen + t_dict + t_eval;
    
    # 3. Error Metrics (L1 and L_inf)
    l1_err := 0.0;
    linf_err := 0.0;
    for i in [1 .. q_val] do
        diff := AbsoluteValue(explicit_dist[i] - minspm_dist[i]);
        l1_err := l1_err + diff;
        if diff > linf_err then
            linf_err := diff;
        fi;
    od;
    
    Print("-----------------------------------------------------------------\n");
    Print(" REGISTER m = ", m_bits, " (Q = ", q_val, ")\n");
    Print("-----------------------------------------------------------------\n");
    Print("  - Unique Classes (u)          : ", u_classes, "\n");
    Print("  - Gen Time / Dict Time / Eval : ", t_gen, "s / ", t_dict, "s / ", t_eval, "s\n");
    Print("  - Total Runtime               : ", t_total, "s\n");
    Print("  - L1 Distribution Error       : ", l1_err, "\n");
    Print("  - L_inf Distribution Error    : ", linf_err, "\n");
    Print("  - Best y (Explicit vs MINSPM) : ", best_y_exp, " vs ", best_y_min, "\n");
    Print("  - Peak Prob (Exp vs MINSPM)   : ", max_p_exp, " vs ", max_p_min, "\n");
    Print("-----------------------------------------------------------------\n");
    
    return [l1_err, linf_err, t_total, u_classes];
end;

run_full_validation_suite := function()
    local m_list, m, true_a;
    m_list := [12, 14, 16]; 
    true_a := 0.0816202435301607;
    
    Print("\n=================================================================\n");
    Print(" FULL-SPECTRUM VALIDATION SUITE (L1, L_inf, Timing Breakdown)\n");
    Print("=================================================================\n");
    
    for m in m_list do
        validate_spectrum_precision(m, true_a);
    od;
end;

run_full_validation_suite();
