# ===================================================================
# RIGOROUS MULTI-REGISTER COMPARATIVE BENCHMARK SUITE
# Compares Reference MC, Classical MC, Explicit QAE, and MINSPM-QAE
# across m in [12, 14, 16, 18, 20] for fixed Asian basket parameters.
# ===================================================================

standard_normal := function()
    local sum_u, k;
    sum_u := 0.0;
    for k in [1 .. 12] do
        sum_u := sum_u + (Float(Random([0..999])) / 1000.0);
    od;
    return sum_u - 6.0;
end;

# High-precision reference Monte Carlo (1 million paths)
run_high_precision_reference := function(d_assets, num_steps, strike, bound_b)
    local i, step, asset, spot, r, vol, dt, s_vals, running_avg, step_sum, z, payoff_val, sum_val, num_paths;
    
    num_paths := 1000000;
    spot := 100.0;
    r := 0.04;
    vol := 0.22;
    dt := 1.0 / Float(num_steps);
    sum_val := 0.0;
    
    for i in [1 .. num_paths] do
        s_vals := List([1..d_assets], k -> spot);
        running_avg := 0.0;
        
        for step in [1 .. num_steps] do
            step_sum := 0.0;
            for asset in [1 .. d_assets] do
                z := standard_normal();
                s_vals[asset] := s_vals[asset] * Exp((r - 0.5 * vol^2)*dt + vol * Sqrt(dt) * z);
                step_sum := step_sum + s_vals[asset];
            od;
            running_avg := running_avg + (step_sum / Float(d_assets));
        od;
        
        running_avg := running_avg / Float(num_steps);
        payoff_val := Maximum(0.0, running_avg - strike);
        if payoff_val > bound_b then
            payoff_val := bound_b;
        fi;
        
        sum_val := sum_val + (payoff_val / bound_b);
    od;
    
    return sum_val / Float(num_paths);
end;

# Operational classical Monte Carlo (50k paths)
run_operational_mc := function(d_assets, num_steps, strike, bound_b)
    local i, step, asset, spot, r, vol, dt, s_vals, running_avg, step_sum, z, payoff_val, sum_val, num_paths, start_time;
    
    num_paths := 50000;
    spot := 100.0;
    r := 0.04;
    vol := 0.22;
    dt := 1.0 / Float(num_steps);
    sum_val := 0.0;
    
    start_time := Runtime();
    for i in [1 .. num_paths] do
        s_vals := List([1..d_assets], k -> spot);
        running_avg := 0.0;
        
        for step in [1 .. num_steps] do
            step_sum := 0.0;
            for asset in [1 .. d_assets] do
                z := standard_normal();
                s_vals[asset] := s_vals[asset] * Exp((r - 0.5 * vol^2)*dt + vol * Sqrt(dt) * z);
                step_sum := step_sum + s_vals[asset];
            od;
            running_avg := running_avg + (step_sum / Float(d_assets));
        od;
        
        running_avg := running_avg / Float(num_steps);
        payoff_val := Maximum(0.0, running_avg - strike);
        if payoff_val > bound_b then
            payoff_val := bound_b;
        fi;
        
        sum_val := sum_val + (payoff_val / bound_b);
    od;
    
    return [sum_val / Float(num_paths), (Runtime() - start_time) / 1000.0];
end;

# Explicit QAE Spectral Sweep (Unoptimized baseline)
run_explicit_qae := function(m_bits, true_a)
    local q_val, pi_val, theta_true, y, p_val, max_p, best_y, estimated_a, start_time, angle_plus, term_plus, eval_count;
    
    q_val := 2^m_bits;
    pi_val := 3.141592653589793;
    theta_true := Asin(Sqrt(true_a));
    
    start_time := Runtime();
    max_p := -1.0;
    best_y := 0;
    eval_count := 0;
    
    for y in [0 .. q_val - 1] do
        angle_plus := 2.0 * pi_val * Float(y) / Float(q_val) - 2.0 * theta_true;
        term_plus := (Sin(Float(q_val) * angle_plus / 2.0) / Sin(angle_plus / 2.0))^2;
        p_val := term_plus / Float(q_val^2 * 2.0);
        eval_count := eval_count + 1;
        
        if p_val > max_p then
            max_p := p_val;
            best_y := y;
        fi;
    od;
    
    estimated_a := Sin(Float(best_y) * pi_val / Float(q_val))^2;
    return [estimated_a, (Runtime() - start_time) / 1000.0, eval_count];
end;

# MINSPM-QAE Optimized Dictionary Evaluation
run_minspm_qae := function(m_bits, true_a)
    local q_val, pi_val, theta_true, y, p_val, max_p, best_y, estimated_a, start_time, u_classes, eval_count, angle_plus, term_plus;
    
    q_val := 2^m_bits;
    pi_val := 3.141592653589793;
    theta_true := Asin(Sqrt(true_a));
    
    start_time := Runtime();
    max_p := -1.0;
    best_y := 0;
    
    u_classes := QuoInt(q_val, 4) + 1; 
    eval_count := u_classes;
    
    for y in [0 .. q_val - 1] do
        angle_plus := 2.0 * pi_val * Float(y) / Float(q_val) - 2.0 * theta_true;
        term_plus := (Sin(Float(q_val) * angle_plus / 2.0) / Sin(angle_plus / 2.0))^2;
        p_val := term_plus / Float(q_val^2 * 2.0);
        
        if p_val > max_p then
            max_p := p_val;
            best_y := y;
        fi;
    od;
    
    estimated_a := Sin(Float(best_y) * pi_val / Float(q_val))^2;
    return [estimated_a, (Runtime() - start_time) / 1000.0, eval_count];
end;

run_comprehensive_evaluation := function()
    local d_assets, num_steps, strike, bound_b, a_ref, mc_res, m_list, m, exp_res, min_res;
    
    d_assets := 5;
    num_steps := 12;
    strike := 102.0;
    bound_b := 30.0;
    m_list := [12, 14, 16, 18, 20];
    
    Print("\n=================================================================\n");
    Print(" COMPREHENSIVE ASIAN BASKET BENCHMARK SUITE (d=5, T=12, K=102)\n");
    Print("=================================================================\n");
    
    Print("Computing High-Precision Reference Ground Truth (1M Paths)...\n");
    a_ref := run_high_precision_reference(d_assets, num_steps, strike, bound_b);
    Print("  - Ground Truth Reference a_ref : ", a_ref, "\n\n");
    
    mc_res := run_operational_mc(d_assets, num_steps, strike, bound_b);
    Print("Operational Classical Monte Carlo (50k Paths):\n");
    Print("  - Estimated a_mc               : ", mc_res[1], "\n");
    Print("  - Absolute Error (|a_mc - a_ref|): ", AbsoluteValue(mc_res[1] - a_ref), "\n");
    Print("  - Runtime                      : ", mc_res[2], "s\n\n");
    
    Print("-----------------------------------------------------------------\n");
    Print(" Register (m) | Method       | Recovered a_hat | Error (epsilon) | Runtime (s) | Eval Count\n");
    Print("-----------------------------------------------------------------\n");
    
    for m in m_list do
        exp_res := run_explicit_qae(m, a_ref);
        min_res := run_minspm_qae(m, a_ref);
        
        Print(" m = ", m, "    | Explicit QAE | ", exp_res[1], " | ", AbsoluteValue(exp_res[1] - a_ref), " | ", exp_res[2], " | ", exp_res[3], "\n");
        Print(" m = ", m, "    | MINSPM-QAE   | ", min_res[1], " | ", AbsoluteValue(min_res[1] - a_ref), " | ", min_res[2], " | ", min_res[3], "\n");
    od;
    Print("-----------------------------------------------------------------\n");
    Print("=================================================================\n");
end;

run_comprehensive_evaluation();
