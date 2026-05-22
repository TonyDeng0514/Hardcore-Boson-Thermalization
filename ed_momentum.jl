# ed_momentum.jl — PBC momentum-sector ED for hardcore bosons
#
# Block-diagonalizes H by translational symmetry into L momentum sectors.
# Requires α=0 (Stark tilt breaks translation invariance).
# μ=0: uniform chemical potential shifts all energies by -μN, omitted.
#
# Output: results/ed_ksectors_L{L}_N{N}_tp{tp}_Vp{Vp}.csv
#         columns: sector (m), energy, O_expval  (O = n_{j1}*n_{j2})
#
# Verification: "Total dim" printed at the end must equal binomial(L, N).

using LinearAlgebra, Printf

const L      = 17
const N      = 6
const t_hop  = 1.0
const V      = 1.0
const tp_hop = 0.0
const Vp     = 0.0

# ── bitstring utilities (0-indexed sites; bit j = site j) ────────────────────

translate(s::Int, L::Int) = ((s << 1) | (s >> (L - 1))) & ((1 << L) - 1)

function get_orbit(s::Int, L::Int)
    orbit = Int[s]
    cur = translate(s, L)
    while cur != s
        push!(orbit, cur)
        cur = translate(cur, L)
    end
    orbit
end

# ── enumerate C(L,N) states via Gosper's hack ────────────────────────────────

function n_particle_states(L::Int, N::Int)
    N == 0 && return Int[0]
    states = Int[]
    s = (1 << N) - 1
    lim = 1 << L
    while s < lim
        push!(states, s)
        c = s & -s
        r = s + c
        s = (((r ⊻ s) >> 2) ÷ c) | r
    end
    states
end

# ── state_info[s] = (rep, l, period): T^l(rep) = s ──────────────────────────

function build_state_info(L::Int, N::Int)
    info = Dict{Int,Tuple{Int,Int,Int}}()
    for s in n_particle_states(L, N)
        haskey(info, s) && continue
        orbit  = get_orbit(s, L)
        rep    = minimum(orbit)
        j_rep  = findfirst(==(rep), orbit) - 1   # 0-indexed position of rep
        period = length(orbit)
        for (r, st) in enumerate(orbit)
            info[st] = (rep, mod(r - 1 - j_rep, period), period)
        end
    end
    info
end

# ── list of (rep, period) for sector m  (compatibility: (m*period) % L == 0) ─

function build_k_basis(L::Int, N::Int, m::Int, info::Dict)
    visited = Set{Int}()
    basis   = Tuple{Int,Int}[]
    for s in n_particle_states(L, N)
        s in visited && continue
        rep, _, period = info[s]
        for st in get_orbit(rep, L); push!(visited, st); end
        (m * period) % L == 0 && push!(basis, (rep, period))
    end
    basis
end

# ── build complex Hermitian H for sector m (O(d×L), PBC) ────────────────────

function build_Hk(L::Int, m::Int, t::Float64, V::Float64,
                  tp::Float64, Vp::Float64,
                  info::Dict, basis::Vector)
    k   = 2π * m / L
    d   = length(basis)
    H   = zeros(ComplexF64, d, d)
    col_of = Dict(rep => col for (col, (rep, _)) in enumerate(basis))

    for (col, (ket, Rα)) in enumerate(basis)
        # diagonal: V n_j n_{j+1} + Vp n_j n_{j+2}  (all L bonds, PBC)
        diag = 0.0
        for j in 0:L-1
            nj  = (ket >> j)         & 1
            nj1 = (ket >> ((j+1)%L)) & 1
            nj2 = (ket >> ((j+2)%L)) & 1
            diag += V  * nj * nj1
            diag += Vp * nj * nj2
        end
        H[col, col] += diag

        # hopping helper: add one off-diagonal element
        function hop!(amp::Float64, new_s::Int)
            bra, l, Rβ = info[new_s]
            haskey(col_of, bra) || return
            H[col_of[bra], col] += amp * exp(-im * k * l) * sqrt(Rα / Rβ)
        end

        # NN hopping -t (PBC: wraps at j=L-1 → jn=0)
        for j in 0:L-1
            jn = (j + 1) % L
            ((ket >> j) & 1) != ((ket >> jn) & 1) || continue
            hop!(-t, ket ⊻ (1 << j) ⊻ (1 << jn))
        end

        # NNN hopping -tp (PBC: wraps at j=L-2,L-1)
        if tp != 0.0
            for j in 0:L-1
                jnn = (j + 2) % L
                ((ket >> j) & 1) != ((ket >> jnn) & 1) || continue
                hop!(-tp, ket ⊻ (1 << j) ⊻ (1 << jnn))
            end
        end
    end

    Hermitian(H)
end

# ── ⟨ψ_n|O|ψ_n⟩ for O = n_{j1}*n_{j2} using orbit-average identity ──────────
# ⟨ψ_n|O|ψ_n⟩ = Σ_α |c_{αn}|² × (1/Rα) × Σ_{s∈orbit(α)} O(s)
# (cross-orbit terms vanish because orbits are disjoint)

function obs_expvals(vecs::Matrix{ComplexF64}, basis::Vector, j1::Int, j2::Int, L::Int)
    n_eig = size(vecs, 2)
    obs   = zeros(Float64, n_eig)
    for (col, (rep, period)) in enumerate(basis)
        orbit_O = 0.0
        s = rep
        for _ in 1:period
            orbit_O += Float64((s >> j1) & 1) * Float64((s >> j2) & 1)
            s = translate(s, L)
        end
        orbit_O /= period
        for n in 1:n_eig
            obs[n] += abs2(vecs[col, n]) * orbit_O
        end
    end
    obs
end

# ── main ─────────────────────────────────────────────────────────────────────

# Center bond (0-indexed): sites L÷2-1 and L÷2  →  sites L÷2 and L÷2+1 (1-indexed display)
# j1 = L ÷ 2 - 1
# j2 = L ÷ 2

j1 = 1
j2 = 2

@printf "L=%d  N=%d  t=%.2f  V=%.2f  tp=%.2f  Vp=%.2f\n" L N t_hop V tp_hop Vp
@printf "Building state_info for %d states...\n" binomial(L, N); flush(stdout)
state_info = build_state_info(L, N)

mkpath("results")
outfile = "results/ed_ksectors_L$(L)_N$(N)_tp$(tp_hop)_Vp$(Vp).csv"

open(outfile, "w") do f
    println(f, "# PBC momentum-sector ED  L=$(L) N=$(N) t=$(t_hop) V=$(V) tp=$(tp_hop) Vp=$(Vp)  O=n_$(j1+1)*n_$(j2+1) (1-indexed)")
    println(f, "sector,energy,O_expval")

    dim_total = 0
    for m in 0:L-1
        basis = build_k_basis(L, N, m, state_info)
        isempty(basis) && continue
        Hk         = build_Hk(L, m, t_hop, V, tp_hop, Vp, state_info, basis)
        vals, vecs = eigen(Hk)
        obs        = obs_expvals(vecs, basis, j1, j2, L)
        for n in eachindex(vals)
            @printf f "%d,%.10f,%.10f\n" m vals[n] obs[n]
        end
        dim_total += length(basis)
        @printf "sector m=%2d  dim=%4d\n" m length(basis); flush(stdout)
    end

    @printf "Total dim = %d  (expected binomial(%d,%d) = %d)\n" dim_total L N binomial(L, N)
end

@printf "Saved → %s\n" outfile
