using LinearAlgebra

include("gates.jl")
include("observable.jl")

function run_tebd(L, t_hop, V, tp_hop, Vp, μ; 
    τ=0.005, ttotal=30.0, maxdim=256, cutoff=1e-10,
    α=0.0, outdir="results",
    )
    
    sites=siteinds("S=1/2", L; conserve_qns=true)
    
    # build sites
    state = [isodd(n) ? "Up" : "Dn" for n=1:L]
    psi0  = MPS(sites,state)
    @show flux(psi0)

    # build half_gates
    gates = build_tebd_gates(sites, L, t_hop, V, tp_hop, Vp, μ, τ; α=α)

    # Measure initial energy
    H   = build_hamiltonian_mpo(sites, L, t_hop, V, tp_hop, Vp, μ; α=α)
    E0  = inner(psi0', H, psi0)

    nsteps = round(Int, ttotal / τ)

    b = L ÷ 2  # center bond for entanglement entropy

    mkpath(outdir)
    outfile_profile = joinpath(outdir, "n_profile_vs_t_L$(L)_chi$(maxdim)_alpha$(α).csv")
    println("\n=============== χ = $maxdim ================")
    psi = psi0
    open(outfile_profile, "w") do io
        println(io, "# E0 = $E0")
        println(io, "time,entropy,truncation_error,energy," * join(["n_$j" for j in 1:L], ","))
        for k in 0:nsteps
            println("step $k / $nsteps   t = $(@sprintf("%.6f", k*τ))   χ = $(maxlinkdim(psi))")

            orthogonalize!(psi, b)
            Sz = real.(expect(psi, "Sz"))
            n = Sz .+ 0.5

            _, S, _ = svd(psi[b], (linkind(psi, b-1), siteind(psi, b)))
            sv = [S[i,i]^2 for i in 1:dim(S,1)]
            sv ./= sum(sv)
            ent = -sum(p * log(p) for p in sv if p > 0)

            energy = real(inner(psi', H, psi))

            trunc_err = 0.0
            if k < nsteps
                psi_new = apply(gates, psi; cutoff, maxdim)
                trunc_err = 1 - norm(psi_new)^2
                normalize!(psi_new)
                psi = psi_new
            end

            println(io, @sprintf("%.6f,%.10e,%.10e,%.10e", k*τ, ent, trunc_err, energy) *
                "," * join([@sprintf("%.10e", n[j]) for j in 1:L], ","))
            flush(io)
            k == nsteps && break
        end
    end
    println("Wrote $(nsteps+1) rows to $(abspath(outfile_profile))")
end