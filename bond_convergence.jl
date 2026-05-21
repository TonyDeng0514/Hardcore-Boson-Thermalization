include("run_tebd.jl")

# --- CLI entry point ---
if length(ARGS) != 1
    println(stderr, "Usage: julia bond_convergence.jl <chi>")
    println(stderr, "  e.g. julia bond_convergence.jl 64")
    exit(1)
end

χ = tryparse(Int, ARGS[1])
if isnothing(χ) || χ <= 0
    println(stderr, "Error: chi must be a positive integer, got '$(ARGS[1])'")
    exit(1)
end

run_tebd(20, 1.0, 1.0, 0.98, 0.98, 0.0; ttotal=40.0,maxdim=χ, α=0.0)
