include("participation_entropy.jl")

L=18
t_hop=1.0
V=1.0
tp_hop=0.3
Vp=0.3

run_participation_entropy(L, t_hop, V, tp_hop,Vp;n_samples=1, t_max=30.0)