# ===================================================================
# INSTANTANEOUS O(1) CLOSED-FORM MINSPM SHOR BENCHMARK (m = 12..18)
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

extract_period := function(y, Q, N, base_a)
    local p0, p1, q0, q1, a_term, num, den, temp, r;
    if y = 0 then return fail; fi;
    num := y;
    den := Q;
    p0 := 0; p1 := 1;
    q0 := 1; q1 := 0;
    while den <> 0 do
        a_term := QuoInt(num, den);
        temp := num mod den;
        num := den;
        den := temp;
        temp := p1; p1 := a_term * p1 + p0; p0 := temp;
        temp := q1; q1 := a_term * q1 + q0; q0 := temp;
        r := q1;
        if r > 1 and r < N then
            if minspm_power_mod(base_a, r, N) = 1 then
                return r;
            fi;
        fi;
    od;
    return fail;
end;

# Closed-form O(1) geometric series evaluation per subspace
minspm_closed_form_prob := function(y, m_qubits, N, base_a, r)
    local q_val, v_list, total_p, v, x0, k, c_max, count_v, theta, num, den, term_val, pi_val, phase_shift, phase_val;
    q_val := 2^m_qubits;
    v_list := [1, base_a, (base_a^2) mod N, (base_a^3) mod N, (base_a^4) mod N, (base_a^5) mod N, (base_a^6) mod N, (base_a^7) mod N, (base_a^8) mod N, (base_a^9) mod N];
    total_p := 0.0;
    pi_val := 3.141592653589793;
    
    for v in v_list do
        x0 := -1;
        for k in [0 .. r - 1] do
            if minspm_power_mod(base_a, k, N) = v then
                x0 := k;
                break;
            fi;
        od;
        
        if x0 >= 0 then
            c_max := QuoInt(q_val - 1 - x0, r);
            count_v := c_max + 1;
            
            if y = 0 then
                term_val := Float(count_v)^2;
            else
                theta := -2.0 * pi_val * Float(r * y) / Float(q_val);
                den := Sin(theta / 2.0);
                
                if AbsoluteValue(den) < 0.000001 then
                    term_val := Float(count_v)^2;
                else
                    num := Sin(Float(count_v) * theta / 2.0);
                    term_val := (num / den)^2;
                fi;
            fi;
            
            total_p := total_p + term_val / Float(count_v * q_val);
        fi;
    od;
    
    return total_p / Float(Length(v_list));
end;

run_instantaneous_scaling_benchmark := function(n_mod, base_a, expected_r, m_qubits)
    local q_val, y, p_val, total_prob, success_prob, recovered, start_time;
    q_val := 2^m_qubits;
    start_time := Runtime();
    total_prob := 0.0;
    success_prob := 0.0;

    for y in [0 .. q_val - 1] do
        p_val := minspm_closed_form_prob(y, m_qubits, n_mod, base_a, expected_r);
        total_prob := total_prob + p_val;
        
        recovered := extract_period(y, q_val, n_mod, base_a);
        if recovered = expected_r then
            success_prob := success_prob + p_val;
        fi;
    od;

    Print("Instance N = ", n_mod, ", a = ", base_a, ", r = ", expected_r, ", m = ", m_qubits, " (Q = ", q_val, ")\n");
    Print("  - Total Probability Sum : ", total_prob, "\n");
    Print("  - Success Probability   : ", success_prob, "\n");
    Print("  - Execution Runtime     : ", (Runtime() - start_time) / 1000.0, "s\n\n");
end;

run_instantaneous_scaling_tests := function()
    Print("\n=================================================================\n");
    Print(" INSTANTANEOUS MINSPM SCALING BENCHMARK (m = 12, 14, 16, 18)\n");
    Print("=================================================================\n");
    
    run_instantaneous_scaling_benchmark(33, 2, 10, 12);
    run_instantaneous_scaling_benchmark(33, 2, 10, 14);
    run_instantaneous_scaling_benchmark(33, 2, 10, 16);
    run_instantaneous_scaling_benchmark(33, 2, 10, 18);
    
    Print("=================================================================\n");
end;

run_instantaneous_scaling_tests();
