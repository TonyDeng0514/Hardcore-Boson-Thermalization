# Usage: julia -t auto run_ed.jl
# -t auto sets Julia threads to all available cores; BLAS threads are matched below.
# Dense eigen on the N-sector matrix (dim ~ C(L,N)) is LAPACK/BLAS-heavy.

using LinearAlgebra, Base.Threads
BLAS.set_num_threads(Threads.nthreads())

include("ED.jl")

L=18
N=12
t_hop=1.0
V=1.0
tp_hop=0.3
Vp=0.3
μ=0.0
α=0.0

run_ed(L, N, t_hop, V, tp_hop, Vp, μ, α)