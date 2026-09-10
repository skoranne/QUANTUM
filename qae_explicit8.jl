# ===================================================================
# DUAL-RESONANCE DIAGNOSTIC & MASS CONCENTRATION AUDIT
# Isolates independent lobe contributions C_1(L) and C_2(L).
# ===================================================================

using Printf

function qae_probability(y::Int64, Q::Int64, theta_true::Float64)
    q_float = Float64(Q)
    angle = 2.0 * pi * Float64(y) / q_float - 2.0 * theta_true
    if abs(sin(angle / 2.0)) < 1e-12
        return 1.0 / (2.0 * q_float^2)
    end
    term = (sin(q_float * angle / 2.0) / sin(angle / 2.0))^2
    return term / (q_float^2 * 2.0)
end

function run_dual_resonance_audit_fixed()
    m_bits = 32
    Q = Int64(1) << m_bits
    true_a = 0.0816202435301607
    theta_true = asin(sqrt(true_a))
    
    y_1 = Int64(round((theta_true / pi) * Float64(Q)))
    y_2 = Q - y_1
    
    println("Resonance 1 center (y_1) : $y_1 (y_1 / Q = $(y_1 / Float64(Q)))")
    println("Resonance 2 center (y_2) : $y_2 (y_2 / Q = $(y_2 / Float64(Q)))")
    println("Absolute separation |y_1 - y_2| = $(abs(y_1 - y_2))\n")
    
    println("=========================================================================================")
    println(" DUAL-RESONANCE MASS CONCENTRATION AUDIT (m = 32, Q = 2^32)")
    println("=========================================================================================")
    @printf("%-8s | %-18s | %-18s | %-18s\n", "Window L", "C_1(L) [Lobe 1]", "C_2(L) [Lobe 2]", "C_union(L)")
    println("-----------------------------------------------------------------------------------------")
    
    l_values = [16, 32, 64, 128, 512, 1024, 4096]
    
    for L in l_values
        c1 = 0.0
        c2 = 0.0
        
        for offset in -L:L
            w1 = y_1 + offset
            w2 = y_2 + offset
            if w1 >= 0 && w1 < Q
                c1 += qae_probability(w1, Q, theta_true)
            end
            if w2 >= 0 && w2 < Q
                c2 += qae_probability(w2, Q, theta_true)
            end
        end
        
        c_union = c1 + c2
        @printf("%-8d | %.12f | %.12f | %.12f\n", L, c1, c2, c_union)
    end
    println("-----------------------------------------------------------------------------------------")
end

run_dual_resonance_audit_fixed()
