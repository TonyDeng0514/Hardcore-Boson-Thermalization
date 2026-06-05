# Usage: julia -t auto run_participation_entropy.jl
# -t auto sets Julia threads to all available cores; BLAS threads are matched below.
# The mul!(ψ_t, T_mat, C_t) inner loop dominates runtime — BLAS threads matter here.

using LinearAlgebra, Base.Threads
BLAS.set_num_threads(Threads.nthreads())

include("participation_entropy.jl")

# L=16
t_hop=1.0
V=1.0
tp_hop=0.3
Vp=0.3

T_end = 100.0

run_participation_entropy(16, t_hop, V, tp_hop,Vp;n_samples=20, t_max=T_end)
run_participation_entropy(14, t_hop, V, tp_hop,Vp;n_samples=20, t_max=T_end)
run_participation_entropy(12, t_hop, V, tp_hop,Vp;n_samples=20, t_max=T_end)
run_participation_entropy(10, t_hop, V, tp_hop,Vp;n_samples=20, t_max=T_end)

run_participation_entropy(16, t_hop, V, 0.0,0.0;n_samples=20, t_max=T_end)
run_participation_entropy(14, t_hop, V, 0.0,0.0;n_samples=20, t_max=T_end)
run_participation_entropy(12, t_hop, V, 0.0,0.0;n_samples=20, t_max=T_end)
run_participation_entropy(10, t_hop, V, 0.0,0.0;n_samples=20, t_max=T_end)