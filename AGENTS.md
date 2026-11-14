# AGENTS.md

## Scope
- Flat research sandbox for **MINSPM** (Matrix-Induced Non-linear Spectral Product Mapping) experiments.
- Primary stacks: **Julia** (GPU/CUDA-heavy QAP/BQOP/benchmark solvers), **Python** with **Qiskit**/**Qiskit Aer** (QAOA/decoding scripts), **GAP** (`*.gap` quantum/algebra experiments), and **QASM** circuit artifacts. Not a packaged library; most files are standalone experiment scripts.

## Environment & running
- **Julia**: project env in `Project.toml` / `Manifest.toml` (manifest records Julia 1.13.1). From repo root:
  - `julia --project=. -e 'using Pkg; Pkg.instantiate()'` (first-time setup)
  - Run a script: `julia --project=. <script>.jl`
- Many `.jl` scripts (`minspm*.jl`, `bqop*.jl`, `*_solver*.jl`, benchmarks, `rcs_*`, `mps_*`, finance/option scripts using `CUDA`) assume a working **NVIDIA CUDA** GPU and will fail or be meaningless without one. Check `using CUDA` / `CuArray` before expecting CPU-only runs.
- **Python**: ad-hoc scripts (e.g. `qaoa_decode.py` imports `qiskit`, `qiskit_aer`); there is no pinned Python env in-repo—install required packages yourself. Run with `python <script>.py`.
- **GAP**: `*.gap` files are GAP scripts (Shor/AQFT/MINSPM algebra experiments); run with `gap <script>.gap` if GAP is installed.
- No automated test suite, CI, formatter, or lint config in this repo; do not invent test commands.

## Data, paths, and layout quirks
- Scripts often **hardcode input/output filenames and paths** and expect to run from the repo root (e.g. `tai256c.dat`/`tai256c.txt`, `qaplib/` instances, local CSVs). Inspect a script's `load_*`/path constants before running; do not assume CLI args.
- `qaplib/` contains QAPLIB instance data/helpers (has its own `.git`); `GEN_CSV/`, `*.csv`, `*.log`, `solution*.txt`, and similar files in root are experiment inputs/outputs—treat as user data/artifacts, avoid overwriting/deleting casually.
- `ddsim/` is a vendored/built copy of **MQT DDSIM** (includes its own `.git`, `build/`, docs); avoid editing it unless the task explicitly targets DDSIM.
- `README.md` is primarily research/notes content, not operational setup—trust script/config behavior over prose when they differ.

## Conventions
- Preserve the existing flat, script-per-experiment style and naming; prefer adding a focused script over introducing package/framework structure unless requested.
- When modifying Julia experiments, keep consistency with surrounding patterns (CUDA kernels, `CuArray`/`Int32|Float32` choices, batch/kernel conventions vary by script—match the file you edit).
- Root contains many generated/intermediate artifacts (`circuit.ll`, logs, CSVs, QASM variants); do not "clean up" unrelated files. 
