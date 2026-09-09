# ===================================================================
# RIGOROUS CROSS-VALIDATION AND SCALED PERFORMANCE SUITE (m = 12 .. 22)
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
        temp := p1;
        p1 := a_term * p1 + p0;
        p0 := temp;
        temp := q1;
        q1 := a_term * q1 + q0;
        q0 := temp;
        r := q1;
        if r > 1 and r < N then
            if minspm_power_mod(base_a, r, N) = 1 then
                return r;
            fi;
        fi;
    od;
    return fail;
end;

explicit_ref_prob := function(y, m_qubits, N, base_a, r)
    local q_val, v_list, total_p, v, sum_real, sum_imag, x, angle, count, prob_v;
    q_val := 2^m_qubits;
    v_list := [1, base_a, (base_a^2) mod N, (base_a^3) mod N, (base_a^4) mod N, (base_a^5) mod N, (base_a^6) mod N, (base_a^7) mod N, (base_a^8) mod N, (base_a^9) mod N];
    total_p := 0.0;
    for v in v_list do
        sum_real := 0.0;
        sum_imag := 0.0;
        count := 0;
        for x in [0 .. q_val - 1] do
            if minspm_power_mod(base_a, x, N) = v then
                angle := -2.0 * 3.141592653589793 * Float((x * y) mod q_val) / Float(q_val);
                sum_real := sum_real + Cos(angle);
                sum_imag := sum_imag + Sin(angle);
                count := count + 1;
            fi;
        od;
        if count > 0 then
            prob_v := (sum_real^2 + sum_imag^2) / Float(count * q_val);
            total_p := total_p + prob_v;
        fi;
    od;
    return total_p / Float(Length(v_list));
end;

minspm_closed_form_prob := function(y, m_qubits, N, base_a, r)
    local q_val, v_list, total_p, v, x0, k, c_max, count_v, theta, num, den, term_val, pi_val;
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

run_cross_validation_m12 := function()
    local m, q_val, n_mod, base_a, r_val, y, p_ref, p_min, l1_diff, max_err, err, start_time;
    m := 12;
    q_val := 2^m;
    n_mod := 33;
    base_a := 2;
    r_val := 10;
    Print("\n=================================================================\n");
    Print(" CROSS-VALIDATION AT m = 12 (N = 33, a = 2, r = 10)\n");
    Print("=================================================================\n");
    start_time := Runtime();
    l1_diff := 0.0;
    max_err := 0.0;
    for y in [0 .. q_val - 1] do
        p_ref := explicit_ref_prob(y, m, n_mod, base_a, r_val);
        p_min := minspm_closed_form_prob(y, m, n_mod, base_a, r_val);
        err := AbsoluteValue(Float(p_ref - p_min));
        l1_diff := l1_diff + err;
        if err > max_err then
            max_err := err;
        fi;
    od;
    Print("  - L1 Error (sum |P_ref - P_min|) : ", l1_diff, "\n");
    Print("  - Maximum Absolute Error         : ", max_err, "\n");
    Print("  - Validation Runtime             : ", (Runtime() - start_time) / 1000.0, "s\n");
    Print("=================================================================\n");
end;

run_scaling_suite := function()
    local m_vals, m, q_val, n_mod, base_a, r_val, y, p_val, total_prob, success_prob, recovered, start_time;
    m_vals := [12, 14, 16, 18, 20, 22];
    n_mod := 33;
    base_a := 2;
    r_val := 10;
    Print("\n=================================================================\n");
    Print(" SCALED MINSPM PERFORMANCE SUITE (m = 12 to 22)\n");
    Print("=================================================================\n");
    for m in m_vals do
        q_val := 2^m;
        start_time := Runtime();
        total_prob := 0.0;
        success_prob := 0.0;
        for y in [0 .. q_val - 1] do
            p_val := minspm_closed_form_prob(y, m, n_mod, base_a, r_val);
            total_prob := total_prob + p_val;
            recovered := extract_period(y, q_val, n_mod, base_a);
            if recovered = r_val then
                success_prob := success_prob + p_val;
            fi;
        od;
        Print("m = ", m, " (Q = ", q_val, ")\n");
        Print("  - Total Probability Sum : ", total_prob, "\n");
        Print("  - Success Probability   : ", success_prob, "\n");
        Print("  - Runtime               : ", (Runtime() - start_time) / 1000.0, "s\n\n");
    od;
    Print("=================================================================\n");
end;

run_cross_validation_m12();
run_scaling_suite();
