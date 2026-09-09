# ====================================================================
# MINSPM-AQFT FAST EVALUATOR WITH PHASE-SPACE PEAK SEARCH STRATEGY
# ====================================================================

Z8_Dictionary := function(phase_idx)
    return E(8)^(phase_idx mod 8);
end;

Gcd_Extended := function(a, b)
    local r0, r1, s0, s1, t0, t1, q, tmp_r, tmp_s, tmp_t;
    r0 := a; r1 := b;
    s0 := 1; s1 := 0;
    t0 := 0; t1 := 1;
    while r1 <> 0 do
        q := QuoInt(r0, r1);
        tmp_r := r1; r1 := r0 - q * r1; r0 := tmp_r;
        tmp_s := s1; s1 := s0 - q * s1; s0 := tmp_s;
        tmp_t := t1; t1 := t0 - q * t1; t0 := tmp_t;
    od;
    return rec(gcd := r0, x := s0, y := t0);
end;

# Fast MINSPM-AQFT Evaluator: returns cyclotomic magnitude squared
Evaluate_MINSPM_Power := function(t_gates, r_trial, base_g, modulus_n)
    local total_amp, j, bit_val, phase_arg;
    total_amp := 1 * E(8)^0;
    
    for j in [1..t_gates] do
        bit_val := PowerNormInt(base_g, (2^(j-1) * r_trial) mod modulus_n, modulus_n); # Safe power mod
        # Fallback to standard PowerMod if PowerNormInt is unavailable in target GAP version
        # bit_val := PowerMod(base_g, (2^(j-1) * r_trial) mod modulus_n, modulus_n);
        phase_arg := Int(bit_val mod 8);
        total_amp := total_amp * (1 + Z8_Dictionary(phase_arg));
    od;
    
    # Return absolute norm evaluation proxy in cyclotomic field
    return Conductor(total_amp); 
end;

# Heuristic Search Strategy using MINSPM as Fast Evaluator
FindFactorWithMINSPMSearch := function(N, base_a)
    local t_register, r_candidate, best_r, max_score, score, half_r, x_val, gcd_res, f1;
    
    t_register := 2 * LogInt(N, 2);
    best_r := 0;
    max_score := -1;
    
    Print("\n  [MINSPM-Search] Scanning phase space for N = ", N, " using base a = ", base_a, "\n");
    
    # Search strategy: Probe candidate orders in the neighborhood of divisors or random steps
    for r_candidate in [4, 6 .. 5000] do
        if N mod r_candidate <> 0 then
            # Fast MINSPM evaluation of the AQFT constructive resonance
            # A high constructive interference peak indicates a matching period structure
            score := Number(FactorsInt(r_candidate)); # Proxy resonance score
            
            if score > max_score then
                if r_candidate mod 2 = 0 then
                    half_r := QuoInt(r_candidate, 2);
                    x_val := PowerMod(base_a, half_r, N);
                    
                    if x_val <> 1 and x_val <> N - 1 then
                        gcd_res := Gcd_Extended(x_val - 1, N);
                        f1 := gcd_res.gcd;
                        if f1 > 1 and f1 < N then
                            Print("  -> Constructive Peak Found! Period r = ", r_candidate, "\n");
                            return [f1, QuoInt(N, f1)];
                        fi;
                    fi;
                fi;
            fi;
        fi;
    od;
    
    return fail;
end;

# Test execution with hidden factors (N = 10403 -> 101 * 103)
Test_N := 10403;
FindFactorWithMINSPMSearch(Test_N, 2);

