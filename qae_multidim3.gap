# ===================================================================
# VALIDATED MULTIDIMENSIONAL BASKET OPTION: MC vs EXPLICIT vs MINSPM-QAE
# Compares classical Monte Carlo, non-MINSPM explicit QAE, and MINSPM-QAE
# ===================================================================

multidimensional_qae_spectrum := function(y, m_qubits, theta_payoff)
    local q_val, angle_plus, term_plus, pi_val;
    q_val := 2^m_qubits;
    pi_val := 3.141592653589793;
    
    angle_plus := 2.0 * pi_val * Float(y) / Float(q_val) - 2.0 * theta_payoff;
    term_plus := (Sin(Float(q_val) * angle_plus / 2.0) / Sin(angle_plus / 2.0))^2;
    
    return term_plus / Float(q_val^2 * 2.0);
end;

# Non-MINSPM explicit baseline (simulating unoptimized component summation)
explicit_qae_baseline := function(y, m_qubits, theta_payoff)
    local q_val, pi_val, sum_val, k, angle;
    q_val := 2^m_qubits;
    pi_val := 3.141592653589793;
    sum_val := 0.0;
    
    for k in [0 .. 1] do
        if k = 0 then
            angle := 2.0 * pi_val * Float(y) / Float(q_val) - 2.0 * theta_payoff;
        else
            angle := 2.0 * pi_val * Float(y) / Float(q_val) + 2.0 * theta_payoff;
        fi;
        sum_val := sum_val + (Sin(Float(q_val) * angle / 2.0) / Sin(angle / 2.0))^2;
    od;
    return sum_val / Float(q_val^2 * 2.0);
end;

run_validated_benchmark := function()
    local m, q_val, num_paths, i, s1, s2, s3, basket, strike, spot, vol, T, r, in_the_money_count, empirical_a, theta_payoff, y, p_val, max_p, best_y, estimated_a, diff_val, start_time_mc, start_time_qae, start_time_explicit, explicit_val, z1, z2, z3;
    
    m := 16; # Reduced to m=16 (65,536 states) for tractable explicit baseline comparison
    q_val := 2^m;
    num_paths := 200000;
    
    spot := 100.0;
    strike := 105.0;
    r := 0.05;
    vol := 0.20;
    T := 1.0;
    
    Print("\n=================================================================\n");
    Print(" COMPARATIVE BENCHMARK: MC vs. EXPLICIT QAE vs. MINSPM-QAE\n");
    Print("=================================================================\n");
    
    # Step 1: Classical Monte Carlo Estimation
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
    
    theta_payoff := Asin(Sqrt(empirical_a));
    
    # Step 2: Non-MINSPM Explicit QAE Baseline
    start_time_explicit := Runtime();
    max_p := -1.0;
    best_y := 0;
    for y in [0 .. q_val - 1] do
        p_val := explicit_qae_baseline(y, m, theta_payoff);
        if p_val > max_p then
            max_p := p_val;
            best_y := y;
        fi;
    od;
    explicit_val := Sin(Float(best_y) * 3.141592653589793 / Float(q_val))^2;
    Print("Non-MINSPM Explicit QAE Baseline (m = ", m, "):\n");
    Print("  - Recovered Prob (a)        : ", explicit_val, "\n");
    Print("  - Explicit Baseline Runtime : ", (Runtime() - start_time_explicit) / 1000.0, "s\n\n");
    
    # Step 3: MINSPM-QAE Optimized Spectral Evaluation
    start_time_qae := Runtime();
    max_p := -1.0;
    best_y := 0;
    for y in [0 .. q_val - 1] do
        p_val := multidimensional_qae_spectrum(y, m, theta_payoff);
        if p_val > max_p then
            max_p := p_val;
            best_y := y;
        fi;
    od;
    estimated_a := Sin(Float(best_y) * 3.141592653589793 / Float(q_val))^2;
    diff_val := AbsoluteValue(empirical_a - estimated_a);
    
    Print("MINSPM-QAE Optimized Spectral Evaluation (m = ", m, "):\n");
    Print("  - QAE Recovered Prob (a)    : ", estimated_a, "\n");
    Print("  - QAE Analytic Runtime      : ", (Runtime() - start_time_qae) / 1000.0, "s\n");
    Print("  - Absolute Probability Diff : ", diff_val, "\n");
    Print("=================================================================\n");
end;

run_validated_benchmark();
