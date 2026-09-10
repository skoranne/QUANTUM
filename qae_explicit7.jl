# ===================================================================
# DUAL-RESONANCE CUMULATIVE CAPTURED PROBABILITY AUDIT C_2(L)
# Computes mass concentration across both principal QAE resonance lobes.
# ===================================================================

using Printf
using LinearAlgebra

function qae_probability(y::Int64, Q::Int64, theta_true::Float64)
    q_float = Float64(Q)
    angle = 2.0 * pi * Float64(y) / q_float - 2.0 * theta_true
    if abs(sin(angle / 2.0)) < 1e-12
        return 1.0 / (2.0 * q_float^2)
    end
    term = (sin(q_float * angle / 2.0) / sin(angle / 2.0))^2
    return term / (q_float^2 * 2.0)
end

function run_dual_resonance_audit()
    println("=========================================================================================")
    println(" DUAL-RESONANCE CUMULATIVE PROBABILITY AUDIT C_2(L) (m = 32, Q = 2^32)")
    println("=========================================================================================")
    @printf("%-8s | %-18s | %-25s\n", "Window L", "Total States Searched", "Captured Probability C_2(L)")
    println("-----------------------------------------------------------------------------------------")
    
    m_bits = 32
    Q = Int64(1) << m_bits
    true_a = 0.0816202435301607
    theta_true = asin(sqrt(true_a))
    
    # Identify the two principal resonance centers
    y_1 = Int64(round((theta_true / pi) * Float64(Q)))
    y_2 = Q - y_1
    
    l_values = [16, 32, 64, 128, 512, 1024]
    
    for L in l_values
        # Construct dual windows W_1 and W_2 without double counting
        searched_states = Set{Int64}()
        
        for offset in -L:L
            w1 = y_1 + offset
            w2 = y_2 + offset
            if w1 >= 0 && w1 < Q
                push!(searched_states, w1)
            end
            if w2 >= 0 && w2 < Q
                push!(searched_states, w2)
            end
        end
        
        total_p = 0.0
        for y in searched_states
            total_p += qae_probability(y, Q, theta_true)
        end
        
        @printf("%-8d | %-18d | %.12f\n", L, length(searched_states), total_p)
    end
    println("-----------------------------------------------------------------------------------------\n")
end

run_dual_resonance_audit()
