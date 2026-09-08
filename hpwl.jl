using CUDA
using Random
using BenchmarkTools

function sparse_hpwl_kernel!(cell_x, cell_y, net_pins, net_offsets, hpwl_out, num_nets)
    net_idx = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if net_idx > num_nets
        return
    end

    start_idx = net_offsets[net_idx]
    end_idx = net_offsets[net_idx + 1] - 1

    min_x = Inf32
    max_x = -Inf32
    min_y = Inf32
    max_y = -Inf32

    for p in start_idx:end_idx
        cell_id = net_pins[p]
        x = cell_x[cell_id]
        y = cell_y[cell_id]
        
        min_x = min(min_x, x)
        max_x = max(max_x, x)
        min_y = min(min_y, y)
        max_y = max(max_y, y)
    end

    hpwl_out[net_idx] = (max_x - min_x) + (max_y - min_y)
    return nothing
end

function run_million_cell_benchmark()
    num_cells = 1_000_000
    num_nets = 2_000_000
    avg_pins_per_net = 4
    total_pins = num_nets * avg_pins_per_net

    println("Initializing Large-Scale Sparse Placement Benchmark:")
    println("Cells: $num_cells | Nets: $num_nets | Total Pins: $total_pins\n")

    rng = MersenneTwister(42)
    cell_x_cpu = rand(rng, Float32, num_cells) .* 1000.0f0
    cell_y_cpu = rand(rng, Float32, num_cells) .* 1000.0f0

    net_offsets_cpu = collect(1:avg_pins_per_net:(total_pins + 1))
    
    # Safe range sampling preventing negative index generation
    net_pins_cpu = rand(rng, 1:num_cells, total_pins)

    cell_x_gpu = CuArray(cell_x_cpu)
    cell_y_gpu = CuArray(cell_y_cpu)
    net_pins_gpu = CuArray(net_pins_cpu)
    net_offsets_gpu = CuArray(net_offsets_cpu)
    hpwl_gpu = CUDA.zeros(Float32, num_nets)

    threads = 256
    blocks = cld(num_nets, threads)

    @cuda threads=threads blocks=blocks sparse_hpwl_kernel!(cell_x_gpu, cell_y_gpu, net_pins_gpu, net_offsets_gpu, hpwl_gpu, num_nets)
    CUDA.synchronize()

    println("Running Sparse HPWL Kernel Evaluation on 1,000,000 Cells...")
    display(@benchmark CUDA.@sync @cuda threads=$threads blocks=$blocks sparse_hpwl_kernel!($cell_x_gpu, $cell_y_gpu, $net_pins_gpu, $net_offsets_gpu, $hpwl_gpu, $num_nets))
    
    total_hpwl = sum(hpwl_gpu)
    println("\nTotal Wirelength Evaluated Successfully. Global HPWL: $total_hpwl")

    CUDA.reclaim()
    GC.gc()
end

run_million_cell_benchmark()
