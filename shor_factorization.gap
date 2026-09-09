# ===================================================================
# GENERALIZED SPARSE RANDOM MATRIX FACTORIZATION VERIFICATION SUITE
# Tests phi(x,y) = x^T A y + b^T x mod 2^K across varying densities
# ===================================================================

test_sparse_matrix_generalization := function()
    local p_mod, t_dim, densities, d, row, col_idx, k, i, j, a_mat, b_vec, y_vec, c_vec, brute_sum, fact_prod, x_vec, idx, phase_val, rec_dict, u_classes, match_flag;
    
    p_mod := 64; # P = 2^K, K = 6
    t_dim := 14; # t = 14 terms (2^14 = 16384 brute force evaluations)
    densities := [1, 2, 4, 8];
    
    Print("\n=================================================================\n");
    Print(" GENERALIZED RANDOM SPARSE MATRIX FACTORIZATION TEST\n");
    Print("=================================================================\n");
    
    for d in densities do
        a_mat := List([1..t_dim], i -> List([1..t_dim], j -> 0));
        for i in [1..t_dim] do
            col_idx := [];
            while Length(col_idx) < d do
                k := Random([1..t_dim]);
                if not (k in col_idx) then
                    Add(col_idx, k);
                fi;
            od;
            for k in col_idx do
                a_mat[i][k] := Random([1..4]);
            od;
        od;
        
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
        
        rec_dict := rec();
        for i in [1..t_dim] do
            if IsBound(rec_dict.(String(c_vec[i]))) then
                rec_dict.(String(c_vec[i])) := rec_dict.(String(c_vec[i])) + 1;
            else
                rec_dict.(String(c_vec[i])) := 1;
            fi;
        od;
        u_classes := Length(RecNames(rec_dict));
        match_flag := (brute_sum = fact_prod);
        
        Print("Density (non-zeros/row) = ", d, "\n");
        Print("  - Unique Classes (u)            : ", u_classes, " (out of t = ", t_dim, ")\n");
        Print("  - Brute-Force == Factorized     : ", match_flag, "\n\n");
    od;
    Print("=================================================================\n");
end;

test_sparse_matrix_generalization();
