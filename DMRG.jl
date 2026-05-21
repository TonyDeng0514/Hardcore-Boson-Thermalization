using ITensorMPS, ITensors

t_hop = 1
V = 1

tp_hop = 0.1
Vp = 0.1

μ = 0

h_eff = -μ + V + Vp  # effective on-site field


let
    N = 100
    sites = siteinds("S=1/2",N; conserve_qns=true)

    os = OpSum()
    for j=1:N-1
        os += -t_hop,"S+",j,"S-",j+1
        os += -t_hop,"S-",j,"S+",j+1
        os += V,"Sz",j,"Sz",j+1
    end

    for j=1:N-2
        os += -tp_hop,"S+",j,"S-",j+2
        os += -tp_hop,"S-",j,"S+",j+2
        os += Vp,"Sz",j,"Sz",j+2
    end

    for j=1:N
        os += h_eff,"S+",j,"S-",j
    end

    H = MPO(os,sites)

    state = [isodd(n) ? "Up" : "Dn" for n=1:N]
    psi0  = MPS(sites,state)
    @show flux(psi0)

    nsweeps = 6
    maxdim = [8, 16, 32, 64, 128, 256]
    cutoff = [1E-10]

    energy, psi = dmrg(H,psi0;nsweeps,maxdim,cutoff)
    
    return
end