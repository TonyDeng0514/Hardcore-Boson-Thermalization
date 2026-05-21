using LinearAlgebra
using Base.Threads

include("gates.jl")
include("observable.jl")

BLAS.set_num_threads(Threads.nthreads())

function run_tebd(L, t_hop, V, tp_hop, Vp, μ; 
    τ=0.05, ttotal=20.0, maxdim=256, cutoff=1e-10,
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

    nsteps     = round(Int, ttotal / τ)
    times      = [n * τ for n in 0:nsteps]
    n_profile  = zeros(nsteps+1, L)

    begin
        println("\n=============== χ = $maxdim ================")
        psi = psi0
        for k in 0:nsteps
            println("step $k / $nsteps   t = $(@sprintf("%.6f", k*τ))   χ = $(maxlinkdim(psi))")
            Sz = expect(psi, "Sz")
            n_profile[k+1, :] = Sz .+ 0.5
            k == nsteps && break
            psi = apply(gates, psi; cutoff, maxdim)
            normalize!(psi)
        end

        # save full density n_profile
        mkpath(outdir)
        outfile_profile = joinpath(outdir, "n_profile_vs_t_L$(L)_chi$(maxdim)_alpha$(α).csv")
        open(outfile_profile, "w") do io
            println(io, "# E0 = $E0")
            header = "time," * join(["n_$j" for j in 1:L], ",")
            println(io, header)
            for k in eachindex(times)
                row = @sprintf("%.6f", times[k]) *
                    "," * join([@sprintf("%.10e", n_profile[k, j]) for j in 1:L],
                    ",",
                    )
                println(io, row)
            end
        end
        println("Wrote $(length(times)) rows to $(abspath(outfile_profile))")
    end
end