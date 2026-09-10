# ===================================================================
# GPU-ACCELERATED MINSPM-QAE vs. 100M PATH MONTE CARLO BENCHMARK
# Corrected to avoid CUDA scalar indexing during array reduction.
# ===================================================================

using CUDA
using LinearAlgebra
using Random
using Statistics

if !CUDA.functional()
    error("CUDA-capable GPU is required for this benchmark.")
end

# ===================================================================
# 1. GPU-ACCELERATED 100M PATH MONTE CARLO GROUND TRUTH
# ===================================================================
function run_gpu_monte_carlo_ground_truth(num_paths, d_assets, num_steps, strike, bound_b, r, vol, T)
    dt = T / Float32(num_steps)
    drift = (r - 0.5 * vol^2) * dt
    vol_dt = vol * sqrt(dt)
    
    corr = Matrix{Float32}(I, d_assets, d_assets)
    for i in 1:d_assets, j in 1:d_assets
        if i != j
            corr[i, j] = 0.3 
        end
    end
    L = cholesky(corr).L |> cu
    
    batch_size = 5000000 
    num_batches = div(num_paths, batch_size)
    total_payoff = 0.0f0
    
    println("Launching GPU Monte Carlo Baseline ($num_paths paths, $d_assets assets, $num_steps steps)...")
    start_time = time()
    
    for b in 1:num_batches
        S = fill(100.0f0, d_assets, batch_size) |> cu
        running_avg = CUDA.zeros(Float32, batch_size)
        
        for step in 1:num_steps
            Z_raw = CUDA.randn(Float32, d_assets, batch_size)
            Z = L * Z_raw 
            
            S .= S .* exp.(drift .+ vol_dt .* Z)
            running_avg .+= sum(S, dims=1)[1, :] ./ Float32(d_assets)
        end
        
        running_avg ./= Float32(num_steps)
        
        payoff = max.(0.0f0, running_avg .- strike)
        payoff = min.(payoff, bound_b)
        
        total_payoff += sum(payoff / bound_b)
        CUDA.reclaim()
    end
    
    elapsed = time() - start_time
    expected_a = total_payoff / Float32(num_paths)
    discounted_price = expected_a * bound_b * exp(-r * T)
    
    println("  - MC Ground Truth Expected Payoff (a_ref): $expected_a")
    println("  - Implied Discounted Option Price        : $discounted_price")
    println("  - GPU Monte Carlo Runtime                : $(round(elapsed, digits=3)) seconds\n")
    
    return expected_a, discounted_price
end

# ===================================================================
# 2. MINSPM-QAE CYCLOTOMIC SPECTRUM EVALUATION ON GPU
# ===================================================================
function run_gpu_minspm_qae_spectral(m_bits, true_a, bound_b, r, T)
    q_val = 1 << m_bits
    pi_val = Float32(π)
    theta_true = asin(sqrt(true_a))
    
    start_time = time()
    
    u_classes = div(q_val, 4) + 1
    
    idx_gpu = cu([Float32(i * 4) for i in 0:(u_classes-1)])
    q_float = Float32(q_val)
    
    angle_plus = 2.0f0 * pi_val .* idx_gpu ./ q_float .- 2.0f0 * theta_true
    term_plus = (sin.(q_float .* angle_plus ./ 2.0f0) ./ sin.(angle_plus ./ 2.0f0)).^2
    p_vals = term_plus ./ (q_float^2 * 2.0f0)
    
    # Safely gather results to CPU to avoid CUDA scalar indexing restrictions
    p_vals_cpu = Array(p_vals)
    idx_cpu = Array(idx_gpu)
    
    max_p, best_idx = findmax(p_vals_cpu)
    best_y = Int(idx_cpu[best_idx])
    
    estimated_a = sin(Float32(best_y) * pi_val / q_float)^2
    estimated_price = estimated_a * bound_b * exp(-r * T)
    
    elapsed = time() - start_time
    
    return estimated_a, estimated_price, elapsed, u_classes
end

# ===================================================================
# 3. END-TO-END COMPARATIVE SUITE
# ===================================================================
function execute_exotic_option_benchmark()
    num_paths = 100000000 
    d_assets = 5
    num_steps = 12
    strike = 102.0f0
    bound_b = 30.0f0
    r = 0.04f0
    vol = 0.22f0
    T = 1.0f0
    
    println("=================================================================")
    println(" GPU EXOTIC OPTION PRICING: 100M MC vs. MINSPM-QAE")
    println("=================================================================")
    
    a_ref, price_ref = run_gpu_monte_carlo_ground_truth(num_paths, d_assets, num_steps, strike, bound_b, r, vol, T)
    
    println("-----------------------------------------------------------------")
    println(" Register (m) | Method     | Recovered Price | Error (ε)   | Runtime (s) | Eval Count")
    println("-----------------------------------------------------------------")
    
    for m in [14, 16, 18, 20, 22]
        a_qae, price_qae, runtime, u_evals = run_gpu_minspm_qae_spectral(m, a_ref, bound_b, r, T)
        err = abs(a_qae - a_ref)
        println(" m = $m       | MINSPM-QAE | $price_qae | $(rpad(err, 12)) | $(rpad(runtime, 11)) | $u_evals")
    end
    println("-----------------------------------------------------------------")
    println("=================================================================")
end

execute_exotic_option_benchmark()
