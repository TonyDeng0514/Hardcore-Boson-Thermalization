function build_hamiltonian_mpo(sites, L, t_hop, V, tp_hop, Vp, μ; α=0.0,)
    os = OpSum()
    h_eff = -μ + V + Vp 
    for j=1:L-1
        os += -t_hop,"S+",j,"S-",j+1
        os += -t_hop,"S-",j,"S+",j+1
        os += V,"Sz",j,"Sz",j+1
    end

    for j=1:L-2
        os += -tp_hop,"S+",j,"S-",j+2
        os += -tp_hop,"S-",j,"S+",j+2
        os += Vp,"Sz",j,"Sz",j+2
    end

    for j=1:L
        os += h_eff,"S+",j,"S-",j
        os += α * j, "S+", j, "S-", j
    end

    H = MPO(os,sites)

    return H
end