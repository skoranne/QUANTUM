# ===================================================================
# DENSE RANDOM MATRIX FACTORIZATION VERIFICATION (t = 16, K = 8)
# Exhaustively tests all 2^16 = 65,536 states for dense linear phases
# ===================================================================

test_dense_random_matrix := function()
    local p_mod, t_dim, a_mat, b_vec, y_vec, c_vec, brute_sum, fact_prod, x_vec, idx, phase_val, start_time, match_flag, i, j, k, row;
    
    p_mod := 256; # P = 2^K, K = 8
    t_dim := 16;  # 2^16 = 65,536 states
    
    Print("\n=================================================================\n");
    Print(" DENSE RANDOM MATRIX FACTORIZATION EXHAUSTIVE TEST (t = ", t_dim, ")\n");
    Print("=================================================================\n");
    
    start_time := Runtime();
    
    # Fully dense random matrix A and vector b
    a_mat := List([1..t_dim], i -> List([1..t_dim], j -> Random([0..(p_mod - 1)])));
    b_vec := List([1..t_dim], i -> Random([0..(p_mod - 1)]));
    y_vec := List([1..t_dim], i -> Random([0..1]));
    
    c_vec := List([1..t_dim], i -> 0);
    for i in [1..t_dim] do
        row := 0;
        for j in [1..t_dim] do
            row := row + a_mat[i][j] * y_vec[j];
        od;
        c_vec[i] := (row + b_vec[i]) mod p_mod;
    od;
    
    brute_sum := 0 * E(p_mod);
    for idx in [1 .. 2^t_dim] do
        x_vec := List([1..t_dim], k -> 0);
        for k in [1..t_dim] do
            x_vec[k] := QuoInt(idx - 1, 2^(k-1)) mod 2;
        od;
        
        phase_val := 0;
        for i in [1..t_dim] do
            if x_vec[i] = 1 then
                phase_val := phase_val + c_vec[i];
            fi;
        od;
        brute_sum := brute_sum + E(p_mod)^(phase_val mod p_mod);
    od;
    
    fact_prod := 1 * E(p_mod)^0;
    for i in [1..t_dim] do
        fact_prod := fact_prod * (1 + E(p_mod)^(c_vec[i]));
    od;
    
    match_flag := (brute_sum = fact_prod);
    
    Print("  - Exhaustive States Tested      : ", 2^t_dim, "\n");
    Print("  - Dense Brute-Force == Factorized : ", match_flag, "\n");
    Print("  - Verification Runtime          : ", (Runtime() - start_time) / 1000.0, "s\n");
    Print("=================================================================\n");
end;

test_dense_random_matrix();
