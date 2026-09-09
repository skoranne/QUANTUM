# ===================================================================
# STANDARDIZED EUROPEAN OPTION QAE BENCHMARK (Normalized to 1.0)
# ===================================================================

normalized_qae_spectrum := function(y, m_qubits, theta_payoff)
    local q_val, angle_plus, term_plus, pi_val;
    q_val := 2^m_qubits;
    pi_val := 3.141592653589793;
    
    angle_plus := 2.0 * pi_val * Float(y) / Float(q_val) - 2.0 * theta_payoff;
    term_plus := (Sin(Float(q_val) * angle_plus / 2.0) / Sin(angle_plus / 2.0))^2;
    
    return term_plus / Float(q_val^2 * 2.0);
end;

run_standard_qae_benchmark := function()
    local m, q_val, target_prob, theta_payoff, y, p_val, total_prob, max_p, best_y, start_time;
    
    m := 20; 
    q_val := 2^m;
    target_prob := 0.3142; 
    theta_payoff := Asin(Sqrt(target_prob));
    
    start_time := Runtime();
    total_prob := 0.0;
    max_p := -1.0;
    best_y := 0;
    
    for y in [0 .. q_val - 1] do
        p_val := normalized_qae_spectrum(y, m, theta_payoff);
        total_prob := total_prob + p_val;
        if p_val > max_p then
            max_p := p_val;
            best_y := y;
        fi;
    od;
    
    Print("Standard QAE benchmark completed for m = ", m, "\n");
    Print("  - Total Probability Sum : ", total_prob, "\n");
    Print("  - Peak Measurement y    : ", best_y, "\n");
    Print("  - Execution Runtime     : ", (Runtime() - start_time) / 1000.0, "s\n");
end;

run_standard_qae_benchmark();
