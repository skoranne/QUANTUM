# ==============================================================================
# GAP SCRIPT: QUADRATIC BOOLEAN QAE CHALLENGE
# Demonstrating Algebraic Ring Collapsing for Genuine Entangling Couplings
# ==============================================================================

# ------------------------------------------------------------------------------
# INSTANCE A: t=16, K=4 (Nearest Neighbor Coupling)
# a_i = 1, J_{i,i+1} = 3. 
# ------------------------------------------------------------------------------

Solve_Instance_A_BruteForce := function()
    local t, mod_val, amp_sum, idx, temp_idx, x_vec, i, Q_x, start_time;
    t := 16;
    mod_val := 16; # 2^K = 16
    amp_sum := 0 * E(mod_val);
    
    start_time := Runtime();
    
    # 2^16 Brute Force Enumeration
    for idx in [0..(2^t - 1)] do
        # Extract binary vector x
        x_vec := List([1..t], k -> 0);
        temp_idx := idx;
        for i in [1..t] do
            x_vec[i] := temp_idx mod 2;
            temp_idx := Int(temp_idx / 2);
        od;
        
        # Compute Q(x) = sum(a_i x_i) + sum(J_{i,i+1} x_i x_{i+1})
        Q_x := 0;
        for i in [1..t] do
            Q_x := Q_x + 1 * x_vec[i]; # a_i = 1
            if i > 1 then
                Q_x := Q_x + 3 * x_vec[i-1] * x_vec[i]; # J_{i-1, i} = 3
            fi;
        od;
        
        amp_sum := amp_sum + E(mod_val)^(Q_x mod mod_val);
    od;
    
    return [amp_sum, Runtime() - start_time];
end;

Solve_Instance_A_MPS := function()
    local t, mod_val, dp, next_dp, k, x_prev, x_curr, phase, start_time;
    t := 16;
    mod_val := 16;
    start_time := Runtime();
    
    # dp[1] maps to x_{k} = 0, dp[2] maps to x_{k} = 1
    dp := [1 * E(mod_val)^0, 1 * E(mod_val)^1]; # Initialize for k=1 (a_1 = 1)
    
    for k in [2..t] do
        next_dp := [0 * E(mod_val), 0 * E(mod_val)];
        for x_prev in [0, 1] do
            for x_curr in [0, 1] do
                # Phase added by x_k is: a_k*x_k + J * x_{prev} * x_{curr}
                phase := (1 * x_curr + 3 * x_prev * x_curr) mod mod_val;
                next_dp[x_curr+1] := next_dp[x_curr+1] + dp[x_prev+1] * E(mod_val)^phase;
            od;
        od;
        dp := next_dp;
    od;
    
    return [dp[1] + dp[2], Runtime() - start_time];
end;


# ------------------------------------------------------------------------------
# INSTANCE B: t=48, K=6 (Dense All-to-All Coupling)
# a_i = pseudorandom, J_ij = (i * j) mod 64
# ------------------------------------------------------------------------------

Solve_Instance_B_Algebraic_DP := function()
    local t, mod_val, a, dp, next_dp, k, S_prev, S_new, phase_add, start_time;
    t := 48;
    mod_val := 64; # 2^K = 64
    start_time := Runtime();
    
    # Pseudorandom pattern for a_i (deterministic for reproducibility)
    a := List([1..t], i -> (i * 7) mod mod_val);
    
    # dp array tracks the amplitude sum for each equivalence class of S = sum(j*x_j) mod 64
    dp := List([1..mod_val], i -> 0 * E(mod_val));
    dp[0+1] := 1 * E(mod_val)^0; # Start with S=0 for empty set
    
    for k in [1..t] do
        next_dp := List([1..mod_val], i -> 0 * E(mod_val));
        for S_prev in [0..(mod_val-1)] do
            if dp[S_prev+1] <> 0 * E(mod_val) then
                
                # Case 1: We set x_k = 0
                # S remains unchanged, no phase is added.
                next_dp[S_prev+1] := next_dp[S_prev+1] + dp[S_prev+1];
                
                # Case 2: We set x_k = 1
                # Because J_ij = i*j, the phase accumulated from ALL prior x_j=1 is exactly k * S_prev.
                phase_add := (a[k] + k * S_prev) mod mod_val;
                S_new := (S_prev + k) mod mod_val;
                next_dp[S_new+1] := next_dp[S_new+1] + dp[S_prev+1] * E(mod_val)^phase_add;
                
            fi;
        od;
        dp := next_dp;
    od;
    
    return [Sum(dp), Runtime() - start_time];
end;

# ==============================================================================
# EXECUTION & VALIDATION
# ==============================================================================

Print("\n========================================================\n");
Print(" INSTANCE A: 16 Qubits, Nearest-Neighbor Entanglement\n");
Print("========================================================\n");

res_A_brute := Solve_Instance_A_BruteForce();
Print("[Brute Force Enumeration (2^16 states)]\n");
Print("  - Amplitude : ", res_A_brute[1], "\n");
Print("  - Runtime   : ", Float(res_A_brute[2]/1000), " seconds\n\n");

res_A_mps := Solve_Instance_A_MPS();
Print("[Transfer Matrix / MPS (O(t) operations)]\n");
Print("  - Amplitude : ", res_A_mps[1], "\n");
Print("  - Runtime   : ", Float(res_A_mps[2]/1000), " seconds\n");
Print("  - Match     : ", res_A_brute[1] = res_A_mps[1], "\n");

Print("\n========================================================\n");
Print(" INSTANCE B: 48 Qubits, Dense All-To-All Entanglement\n");
Print("========================================================\n");

res_B_alg := Solve_Instance_B_Algebraic_DP();
Print("[Algebraic State-Collapse DP (O(t * 2^K) operations)]\n");
Print("  - Amplitude : ", res_B_alg[1], "\n");
Print("  - Runtime   : ", Float(res_B_alg[2]/1000), " seconds\n");
Print("  - Note      : Brute force (2^48) would take millennia.\n");
Print("========================================================\n");
