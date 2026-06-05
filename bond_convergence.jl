using LinearAlgebra, Base.Threads
BLAS.set_num_threads(Threads.nthreads())

include("run_tebd.jl")

# --- CLI entry point ---
# Usage:  julia -t auto bond_convergence.jl <chi>
#   e.g.  julia -t auto bond_convergence.jl 64
# Bond-dimension sweep:
#   julia -t auto bond_convergence.jl 64 && julia -t auto bond_convergence.jl 128 && julia -t auto bond_convergence.jl 256
# -t auto sets Julia threads to all available cores; BLAS threads are matched via BLAS.set_num_threads above.
if length(ARGS) != 1
    println(stderr, "Usage: julia -t auto bond_convergence.jl <chi>")
    println(stderr, "  e.g. julia -t auto bond_convergence.jl 64")
    exit(1)
end

χ = tryparse(Int, ARGS[1])
if isnothing(χ) || χ <= 0
    println(stderr, "Error: chi must be a positive integer, got '$(ARGS[1])'")
    exit(1)
end

run_tebd(16, 1.0, 1.0, 0.3, 0.3, 0.0; ttotal=100.0,maxdim=χ, α=0.0)
