# ===================================================================
# QUANTUM AMPLITUDE ESTIMATION (QAE) MINSPM BENCHMARK
# Evaluates QAE measurement distribution via cyclotomic factorization
# ===================================================================

minspm_power_mod := function(a, e, n)
    local result, base, exp;
    result := 1 mod n;
    base := a mod n;
    exp := e;
    while exp > 0 do
        if (exp mod 2) = 1 then
            result := (result * base) mod n;
        fi;
        base := (base * base) mod n;
        exp := QuoInt(exp, 2);
    od;
    return result;
end;

qae_minspm_prob := function(y, m_qubits, theta_a)
    local q_val, total_p, angle_plus, angle_minus, term_plus, term_minus, pi_val;
    q_val := 2^m_qubits;
    pi_val := 3.141592653589793;
    
    angle_plus := 2.0 * pi_val * Float(y) / Float(q_val) - 2.0 * theta_a;
    angle_minus := 2.0 * pi_val * Float(y) / Float(q_val) + 2.0 * theta_a;
    
    term_plus := (Sin(Float(q_val) * angle_plus / 2.0) / Sin(angle_plus / 2.0))^2;
    term_minus := (Sin(Float(q_val) * angle_minus / 2.0) / Sin(angle_minus / 2.0))^2;
    
    return (term_plus + term_minus) / Float(q_val^2);
end;

run_qae_benchmark := function()
    local m, q_val, target_a, theta_a, y, p_val, total_prob, max_p, best_y, estimated_theta, estimated_a, start_time;
    
    m := 10; # Control register size (Q = 1024)
    q_val := 2^m;
    target_a := 0.2; # Target probability to estimate
    theta_a := Asin(Sqrt(target_a));
    
    Print("\n=================================================================\n");
    Print(" QUANTUM AMPLITUDE ESTIMATION BENCHMARK (Target a = ", target_a, ", m = ", m, ")\n");
    Print("=================================================================\n");
    
    start_time := Runtime();
    total_prob := 0.0;
    max_p := -1.0;
    best_y := 0;
    
    for y in [0 .. q_val - 1] do
        p_val := qae_minspm_prob(y, m, theta_a);
        total_prob := total_prob + p_val;
        if p_val > max_p then
            max_p := p_val;
            best_y := y;
        fi;
    od;
    
    estimated_theta := Float(best_y) * 3.141592653589793 / Float(q_val);
    estimated_a := Sin(estimated_theta)^2;
    
    Print("  - Total Probability Sum : ", total_prob, "\n");
    Print("  - Peak Measurement y    : ", best_y, "\n");
    Print("  - Estimated Amplitude a : ", estimated_a, "\n");
    Print("  - Execution Runtime     : ", (Runtime() - start_time) / 1000.0, "s\n");
    Print("=================================================================\n");
end;

run_qae_benchmark();
