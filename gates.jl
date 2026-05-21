using ITensorMPS, ITensors
using LinearAlgebra
using Printf


function make_NN_gate(sites, j, t_hop, V, h_eff_j, h_eff_jp1)
    s1, s2 = sites[j], sites[j+1]
    hj =    -t_hop      * op("S+",s1) * op("S-",s2) +
            -t_hop      * op("S-",s1) * op("S+",s2) +
            V           * op("Sz",s1) * op("Sz",s2) +
            -h_eff_j    * op("Sz",s1) * op("Id",s2) +
            -h_eff_jp1  * op("Id",s1) * op("Sz",s2)
    return hj
end

function make_NNN_gate(sites, j, tp_hop, Vp)
    s1, s2 = sites[j], sites[j+2]
    hj =    -tp_hop * op("S+",s1) * op("S-",s2) +
            -tp_hop * op("S-",s1) * op("S+",s2) +
            Vp      * op("Sz",s1) * op("Sz",s2)
    return hj
end

function build_tebd_gates(sites, L, t_hop, V, tp_hop, Vp, μ, τ; α=0.0)
    half_gates = ITensor[]

    for j in 1:(L-1)
        w1 = (j     == 1) ? 1.0 : 0.5
        w2 = (j + 1 == L) ? 1.0 : 0.5
        h_eff_j = -μ + V + Vp + α * j
        h_eff_jp1 = -μ + V + Vp + α * (j+1)
        h = make_NN_gate(sites, j, t_hop, V, w1 * h_eff_j, w2 * h_eff_jp1)
        Gj = exp(-1im * (τ / 2) * h)
        push!(half_gates, Gj)
    end
    
    for j in 1:(L-2)
        h = make_NNN_gate(sites, j, tp_hop, Vp)
        Gj = exp(-1im * (τ / 2) * h)
        push!(half_gates, Gj)
    end

    return [half_gates; reverse(half_gates)]
end