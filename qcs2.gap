# ===================================================================
# Corrected GAP Script: Massive QCS Qubit Scaling Sweep
# ===================================================================

PredictMINSPMScaling := function(max_qubits, k_trunc)
    local n, root_order, ext_degree, vram_bytes, vram_gb;
    
    root_order := 2^k_trunc;
    ext_degree := Dimension(CF(root_order));
    
    Print("=================================================================\n");
    Print(" MINSPM QCS Scaling Matrix (k-truncation: ", k_trunc, ", d: ", ext_degree, ")\n");
    Print("=================================================================\n");
    Print("Qubits (n) | State Vector Size | VRAM Footprint | Hardware Tier\n");
    Print("-----------------------------------------------------------------\n");
    
    for n in [10, 15, 20, 25, 28, 30, 32] do
        if n <= max_qubits then
            vram_bytes := ext_degree * 8 * (2.0^n);
            vram_gb := vram_bytes / (1024.0^3);
            
            Print(String(n), " | ", String(2.0^n), " | ", 
                  String(vram_gb), " GB  | ");
                  
            if vram_gb <= 0.0001 then
                Print("SM L1 / Shared Memory\n");
            elif vram_gb <= 0.04 then
                Print("GPU L2 Cache\n");
            elif vram_gb <= 192.0 then
                Print("Single Enterprise GPU VRAM\n");
            else
                Print("Distributed / Multi-GPU Cluster Required\n");
            fi;
        fi;
    od;
end;

PredictMINSPMScaling(32, 3);
