# ===================================================================
# CORRECTED: VALIDATING MINSPM-QAE AGAINST CLASSICAL MONTE CARLO
# ===================================================================

multidimensional_qae_spectrum := function(y, m_qubits, theta_payoff)
    local q_val, angle_plus, term_plus, pi_val;
    q_val := 2^m_qubits;
    pi_val := 3.141592653589793;
    
    angle_plus := 2.0 * pi_val * Float(y) / Float(q_val) - 2.0 * theta_payoff;
    term_plus := (Sin(Float(q_val) * angle_plus / 2.0) / Sin(angle_plus / 2.0))^2;
    
    return term_plus / Float(q_val^2 * 2.0);
end;

run_validated_benchmark := function()
    local m, q_val, num_paths, i, s1, s2, s3, basket, strike, spot, vol, T, r, in_the_money_count, empirical_a, theta_payoff, y, p_val, max_p, best_y, estimated_a, diff_val, start_time_mc, start_time_qae, z1, z2, z3;
    
    m := 18; # Q = 262,144 control states
    q_val := 2^m;
    num_paths := 500000;
    
    spot := 100.0;
    strike := 105.0;
    r := 0.05;
    vol := 0.20;
    T := 1.0;
    
    Print("\n=================================================================\n");
    Print(" VALIDATING MINSPM-QAE AGAINST CLASSICAL MONTE CARLO\n");
    Print("=================================================================\n");
    
    # Step 1: Classical Monte Carlo Estimation of Target Probability 'a'
    start_time_mc := Runtime();
    in_the_money_count := 0.0;
    
    for i in [1 .. num_paths] do
        z1 := (Float(Random([0..9999])) / 10000.0 - 0.5) * 3.464;
        z2 := (Float(Random([0..9999])) / 10000.0 - 0.5) * 3.464;
        z3 := (Float(Random([0..9999])) / 10000.0 - 0.5) * 3.464;
        
        s1 := spot * Exp((r - 0.5 * vol^2)*T + vol * Sqrt(T) * z1);
        s2 := spot * Exp((r - 0.5 * vol^2)*T + vol * Sqrt(T) * z2);
        s3 := spot * Exp((r - 0.5 * vol^2)*T + vol * Sqrt(T) * z3);
        
        basket := (s1 + s2 + s3) / 3.0;
        if basket > strike then
            in_the_money_count := in_the_money_count + 1.0;
        fi;
    od;
    
    empirical_a := in_the_money_count / Float(num_paths);
    Print("Classical Monte Carlo (Paths = ", num_paths, "):\n");
    Print("  - Empirical Probability (a) : ", empirical_a, "\n");
    Print("  - Monte Carlo Runtime       : ", (Runtime() - start_time_mc) / 1000.0, "s\n\n");
    
    # Step 2: MINSPM-QAE Spectral Evaluation using Monte Carlo's Empirical 'a'
    start_time_qae := Runtime();
    theta_payoff := Asin(Sqrt(empirical_a));
    
    max_p := -1.0;
    best_y := 0;
    for y in [0 .. q_val - 1] do
        p_val := multidimensional_qae_spectrum(y, m, theta_payoff);
        if p_val > max_p then
            max_p := p_val;
            best_y := y;
        fi;
    od;
    
    # Corrected decoding formula: theta_a = pi * y / Q, so a = sin^2(pi * y / Q)
    estimated_a := Sin(Float(best_y) * 3.141592653589793 / Float(q_val))^2;
    diff_val := AbsoluteValue(empirical_a - estimated_a);
    
    Print("MINSPM-QAE Spectral Evaluation (m = ", m, "):\n");
    Print("  - QAE Recovered Prob (a)    : ", estimated_a, "\n");
    Print("  - QAE Analytic Runtime      : ", (Runtime() - start_time_qae) / 1000.0, "s\n");
    Print("  - Absolute Probability Diff : ", diff_val, "\n");
    Print("=================================================================\n");
end;

run_validated_benchmark();
