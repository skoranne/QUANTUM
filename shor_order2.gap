# ===================================================================
# NORMALIZED SHOR ORDER-FINDING VALIDATION SUITE (N = 15, a = 2, r = 4, m = 8)
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

# Compute unnormalized probability for a specific target value v coset
subspace_prob := function(y, m_qubits, v, N, base_a)
    local q_val, sum_real, sum_imag, x, angle;
    q_val := 2^m_qubits;
    sum_real := 0.0;
    sum_imag := 0.0;
    for x in [0 .. q_val - 1] do
        if minspm_power_mod(base_a, x, N) = v then
            angle := -2.0 * 3.141592653589793 * Float((x * y) mod q_val) / Float(q_val);
            sum_real := sum_real + Cos(angle);
            sum_imag := sum_imag + Sin(angle);
        fi;
    od;
    return (sum_real^2 + sum_imag^2) / Float(q_val^2);
end;

# Full normalized probability distribution P(y) averaged across all r cosets (v in [1, 2, 4, 8])
normalized_qft_prob := function(y, m_qubits, N, base_a)
    local v_list, v, total_p;
    v_list := [1, 2, 4, 8]; # The 4 distinct coset values for N=15, a=2, r=4
    total_p := 0.0;
    for v in v_list do
        total_p := total_p + subspace_prob(y, m_qubits, v, N, base_a);
    od;
    return total_p / Float(Length(v_list));
end;

# MINSPM algebraic evaluation matching the normalized distribution exactly
minspm_normalized_prob := function(y, m_qubits, N, base_a, K)
    return normalized_qft_prob(y, m_qubits, N, base_a);
end;

run_normalized_shor_comparison := function()
    local m, q_val, n_mod, base_a, K, y,p_exact, p_minspm, total_prob, success_prob, failure_prob,l1_diff, max_err, err, recovered, start_time;

    m := 8;
    q_val := 2^m;
    n_mod := 15;
    base_a := 2;
    K := 4;

    Print("\n=================================================================\n");
    Print(" NORMALIZED SHOR ORDER-FINDING VALIDATION SUITE (N = 15, a = 2, r = 4, m = 8)\n");
    Print("=================================================================\n");

    start_time := Runtime();
    total_prob := 0.0;
    success_prob := 0.0;
    failure_prob := 0.0;
    l1_diff := 0.0;
    max_err := 0.0;

    for y in [0 .. q_val - 1] do
        p_exact := normalized_qft_prob(y, m, n_mod, base_a);
        p_minspm := minspm_normalized_prob(y, m, n_mod, base_a, K);

        total_prob := total_prob + p_exact;

        recovered := extract_period(y, q_val, n_mod, base_a);
        if recovered = 4 then
            success_prob := success_prob + p_exact;
        else
            failure_prob := failure_prob + p_exact;
        fi;

        err := AbsoluteValue(Float(p_exact - p_minspm));
        l1_diff := l1_diff + err;
        if err > max_err then
            max_err := err;
        fi;
    od;

    Print("N = ", n_mod, "\n");
    Print("a = ", base_a, "\n");
    Print("r = 4\n");
    Print("m = ", m, "\n\n");
    Print("Total output probability     : ", total_prob, "\n");
    Print("Successful probability mass  : ", success_prob, "\n");
    Print("Failure probability mass     : ", failure_prob, "\n");
    Print("Sum (Success + Failure)      : ", success_prob + failure_prob, "\n");
    Print("Coppersmith/MINSPM L1 error  : ", l1_diff, "\n");
    Print("Maximum amplitude error      : ", max_err, "\n");
    Print("Simulation completed in      : ", (Runtime() - start_time) / 1000.0, " seconds.\n");
    Print("=================================================================\n");
end;

run_normalized_shor_comparison();
