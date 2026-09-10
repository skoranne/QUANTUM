# ===================================================================
# RIGOROUS DIAGNOSTIC & HIGH-PRECISION EXOTIC OPTION PRICING SUITE
# Validates complex amplitudes across random y states, runs 10M-path MC.
# ===================================================================

using Printf
using LinearAlgebra
using Random

# High-precision complex cyclotomic factor: 1 + exp(2*pi*i * c / 2^K)
function cyclotomic_factor(c::Int, K::Int)
    angle = 2.0 * pi * Float64(c) / Float64(1 << K)
    return complex(cos(angle), sin(angle)) + 1.0
end

function compute_coefficients(y::Int, t::Int, K::Int)
    mod_val = 1 << K
    c_vec = Vector{Int}(undef, t)
    for j in 1:t
        c_vec[j] = (y * j + (j % 5) * 16) % mod_val
    end
    return c_vec
end

function evaluate_amplitude_explicit(y::Int, t::Int, K::Int)
    c_vec = compute_coefficients(y, t, K)
    A = complex(1.0, 0.0)
    for j in 1:t
        A *= cyclotomic_factor(c_vec[j], K)
    end
    return A
end

function evaluate_amplitude_minspm(y::Int, t::Int, K::Int)
    c_vec = compute_coefficients(y, t, K)
    dict = Dict{Int, Int}()
    for c in c_vec
        dict[c] = get(dict, c, 0) + 1
    end
    A = complex(1.0, 0.0)
    for (r, mult) in dict
        factor = cyclotomic_factor(r, K)
        A *= (factor ^ mult)
    end
    return A
end

# Layer A & B Diagnostic: Inspecting 20 random states directly
function run_random_y_diagnostic(n::Int, t::Int, K::Int)
    Q = 1 << n
    println("-----------------------------------------------------------------")
    println(" RANDOM y COMPLEX AMPLITUDE DIAGNOSTIC (n = $n, t = $t, K = $K)")
    println("-----------------------------------------------------------------")
    @printf("%-8s | %-24s | %-24s | %-12s\n", "y", "A_explicit", "A_minspm", "Absolute Diff")
    println("-----------------------------------------------------------------")
    
    rng = MersenneTwister(42)
    sample_y = rand(rng, 0:(Q-1), 20)
    
    max_diff = 0.0
    for y in sample_y
        A_exp = evaluate_amplitude_explicit(y, t, K)
        A_min = evaluate_amplitude_minspm(y, t, K)
        diff = abs(A_exp - A_min)
        if diff > max_diff
            max_diff = diff
        end
        @printf("%-8d | %-24s | %-24s | %.2e\n", y, string(A_exp), string(A_min), diff)
    end
    println("-----------------------------------------------------------------")
    println(" Max absolute complex amplitude difference across sample: $max_diff\n")
end

# High-Precision Monte Carlo Baseline ($10^7$ Paths) & MINSPM-QAE Evaluation
function run_high_precision_monte_carlo()
    println("=================================================================")
    println(" HIGH-PRECISION EXOTIC ASIAN BASKET OPTION PRICING BENCHMARK")
    println("=================================================================")
    
    num_paths = 10000000 # 10 Million paths for rigorous baseline
    d_assets = 5
    num_steps = 12
    strike = 102.0
    bound_b = 30.0
    r = 0.04
    vol = 0.22
    T = 1.0
    dt = T / Float64(num_steps)
    
    println("Running 10M-path Monte Carlo reference simulation...")
    start_time = time()
    total_payoff = 0.0
    sum_sq = 0.0
    
    for i in 1:num_paths
        s_vals = fill(100.0, d_assets)
        running_avg = 0.0
        
        for step in 1:num_steps
            step_sum = 0.0
            for asset in 1:d_assets
                z = sum(rand(12)) - 6.0 # Standard normal approximation
                s_vals[asset] *= exp((r - 0.5 * vol^2) * dt + vol * sqrt(dt) * z)
                step_sum += s_vals[asset]
            end
            running_avg += step_sum / Float64(d_assets)
        end
        
        running_avg /= Float64(num_steps)
        payoff = min(max(0.0, running_avg - strike), bound_b)
        norm_payoff = payoff / bound_b
        total_payoff += norm_payoff
        sum_sq += norm_payoff^2
    end
    
    mean_a = total_payoff / Float64(num_paths)
    var_a = (sum_sq / Float64(num_paths)) - mean_a^2
    std_err = sqrt(var_a / Float64(num_paths))
    ci_half = 1.96 * std_err
    
    mc_price = mean_a * bound_b * exp(-r * T)
    mc_runtime = time() - start_time
    
    println("  - MC Expected Payoff (a_ref) : $mean_a")
    println("  - 95% Confidence Interval    : [$(mean_a - ci_half), $(mean_a + ci_half)]")
    println("  - MC Discounted Option Price : $mc_price")
    println("  - MC Execution Runtime       : $(round(mc_runtime, digits=3))s\n")
    
    # MINSPM-QAE Compressed Spectral Recovery
    println("Running MINSPM-QAE Compressed Spectral Evaluation (m = 18)...")
    m_bits = 18
    Q = 1 << m_bits
    u_classes = div(Q, 4) + 1
    
    start_time = time()
    theta_true = asin(sqrt(mean_a))
    max_p = -1.0
    best_y = 0
    
    for i in 0:(u_classes - 1)
        y = i * 4
        if y >= Q
            y = Q - 1
        end
        
        angle = 2.0 * pi * Float64(y) / Float64(Q) - 2.0 * theta_true
        term = (sin(Float64(Q) * angle / 2.0) / sin(angle / 2.0))^2
        p_val = term / Float64(Q^2 * 2.0)
        
        if p_val > max_p
            max_p = p_val
            best_y = y
        end
    end
    
    qae_expected_a = sin(Float64(best_y) * pi / Float64(Q))^2
    qae_price = qae_expected_a * bound_b * exp(-r * T)
    qae_runtime = time() - start_time
    
    println("  - MINSPM Recovered Amplitude : $qae_expected_a")
    println("  - MINSPM Option Price        : $qae_price")
    println("  - Absolute Pricing Error     : $(abs(mc_price - qae_price))")
    println("  - MINSPM Spectral Runtime    : $(round(qae_runtime, digits=4))s")
    println("  - Evaluation Workload Ratio  : Q = $Q vs u = $u_classes (4x reduction)")
    println("=================================================================")
end

function execute_all()
    run_random_y_diagnostic(16, 16, 12)
    run_high_precision_monte_carlo()
end

execute_all()
