# participation_entropy.jl — time-dependent participation entropy from momentum-sector ED
#
# χ_ent(t) = -ln(Σ_s |⟨s|ψ(t)⟩|^4)   averaged over random Fock initial states (S_z=0)
#
# The Fock-basis IPR starts at 0 for any single Fock initial state and grows toward ln(D)
# as the state thermalizes.  The thermalization time is where the curve saturates.
#
# Key bridge: T_mat[s, k] = ⟨s|E_k⟩ (D×D unitary, built by orbit-DFT across all k-sectors)
#   forward:  C0 = T_mat' * ψ0   →  energy-eigenstate overlaps
#   evolve:   C_t = C0 .* exp(-iE*t)
#   backward: ψ_t = T_mat  * C_t  →  Fock-basis state at time t
#
# L must be even (N = L÷2 for S_z=0).  Stores a D×D complex matrix: ~188 MB at L=14,
# ~1.3 GB at L=16.  Do not attempt L=18 (D=48620, ~37 GB).
#
# Output: results/participation_entropy_L{L}_tp{tp}_Vp{Vp}.csv

using LinearAlgebra, Printf, Random

# ── parameters ───────────────────────────────────────────────────────────────

const L         = 10
const N         = L ÷ 2       # S_z = 0 sector; requires L even
const t_hop     = 1.0
const V         = 1.0
const tp_hop    = 0.3         # NNN hopping: non-zero breaks integrability; 0.0 = integrable
const Vp        = 0.3

const n_samples = 20          # number of random Fock initial states to average over
const t_max     = 30.0
const dt        = 0.005
const rng_seed  = 42

# ── bitstring utilities (unchanged from ed_momentum.jl) ──────────────────────

translate(s::Int, L::Int) = ((s << 1) | (s >> (L - 1))) & ((1 << L) - 1)

function get_orbit(s::Int, L::Int)
    orbit = Int[s]
    cur   = translate(s, L)
    while cur != s; push!(orbit, cur); cur = translate(cur, L); end
    orbit
end

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

function build_state_info(L::Int, N::Int)
    info = Dict{Int,Tuple{Int,Int,Int}}()
    for s in n_particle_states(L, N)
        haskey(info, s) && continue
        orbit  = get_orbit(s, L)
        rep    = minimum(orbit)
        j_rep  = findfirst(==(rep), orbit) - 1
        period = length(orbit)
        for (r, st) in enumerate(orbit)
            info[st] = (rep, mod(r - 1 - j_rep, period), period)
        end
    end
    info
end

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

function build_Hk(L::Int, m::Int, t::Float64, V::Float64,
                  tp::Float64, Vp::Float64, info::Dict, basis::Vector)
    k      = 2π * m / L
    d      = length(basis)
    H      = zeros(ComplexF64, d, d)
    col_of = Dict(rep => col for (col, (rep, _)) in enumerate(basis))
    for (col, (ket, Rα)) in enumerate(basis)
        diag = 0.0
        for j in 0:L-1
            nj  = (ket >> j)          & 1
            nj1 = (ket >> ((j+1)%L))  & 1
            nj2 = (ket >> ((j+2)%L))  & 1
            diag += V * nj * nj1 + Vp * nj * nj2
        end
        H[col, col] += diag
        function hop!(amp::Float64, new_s::Int)
            bra, l, Rβ = info[new_s]
            haskey(col_of, bra) || return
            H[col_of[bra], col] += amp * exp(-im * k * l) * sqrt(Rα / Rβ)
        end
        for j in 0:L-1
            jn = (j + 1) % L
            ((ket >> j) & 1) != ((ket >> jn) & 1) || continue
            hop!(-t, ket ⊻ (1 << j) ⊻ (1 << jn))
        end
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

# ── build Fock↔energy unitary T_mat[s, k] = ⟨s|E_k⟩ ────────────────────────
#
# Each column of T_mat is an energy eigenstate expressed in the Fock basis.
# Construction: for each k-sector m with orbit rep α (period R_α):
#   |α,m⟩ = (1/√R_α) Σ_{l=0}^{R_α-1} exp(2πi·m·l/L) T^l|α⟩
#   ⟨s|α,m⟩ = exp(2πi·m·l_s/L) / √R_α   where T^{l_s}(α) = s
#   ⟨s|E_n^m⟩ = Σ_α vecs[α,n] · ⟨s|α,m⟩
#
# Columns ordered: sector 0 eigenstates (ascending E), then sector 1, ..., sector L-1.

function build_basis_transform(L, N, t_hop, V, tp_hop, Vp, states, state_idx, info)
    D    = length(states)
    T    = zeros(ComplexF64, D, D)
    eigs = zeros(Float64, D)
    col  = 0                          # next free column in T

    for m in 0:L-1
        basis = build_k_basis(L, N, m, info)
        isempty(basis) && continue
        Hk         = build_Hk(L, m, t_hop, V, tp_hop, Vp, info, basis)
        vals, vecs = eigen(Hk)
        k = 2π * m / L
        d = length(vals)

        for (α, (rep, period)) in enumerate(basis)
            s = rep
            for l in 0:period-1
                s_idx = state_idx[s]
                phase = exp(im * k * l) / sqrt(Float64(period))
                for n in 1:d
                    T[s_idx, col + n] += vecs[α, n] * phase
                end
                s = translate(s, L)
            end
        end

        eigs[col+1:col+d] .= vals
        col += d
        @printf "  sector m=%2d  dim=%d\n" m d
    end
    flush(stdout)
    @assert col == D "Eigenbasis incomplete: filled $col of $D columns"
    T, eigs
end

# ── main ─────────────────────────────────────────────────────────────────────

@assert iseven(L) "L must be even for S_z=0 (N=L÷2)"
@printf "L=%d  N=%d  t_hop=%.2f  V=%.2f  tp_hop=%.2f  Vp=%.2f\n" L N t_hop V tp_hop Vp

states    = n_particle_states(L, N)
D         = length(states)
state_idx = Dict(s => i for (i, s) in enumerate(states))
@printf "D = C(%d,%d) = %d  (T_mat: ~%.0f MB)\n" L N D (D^2 * 16 / 1e6)

info = build_state_info(L, N)

@printf "Building eigenbasis (all %d k-sectors)...\n" L; flush(stdout)
T_mat, all_eigvals = build_basis_transform(L, N, t_hop, V, tp_hop, Vp,
                                            states, state_idx, info)

# T_mat should be unitary: ||T†T - I||_F ≲ √D·ε
err = norm(T_mat' * T_mat - I(D)) / sqrt(Float64(D))
@printf "Unitarity check: ||T†T - I||_F / √D = %.2e  (machine ε = %.2e)\n\n" err eps()

# occupation tables for C(t): precomputed once, O(D·L)
occ    = [Bool((states[i] >> (k-1)) & 1)                                      for i in 1:D, k in 1:L]
nn_occ = [Bool((states[i] >> (k-1)) & 1) * Bool((states[i] >> (k % L)) & 1)  for i in 1:D, k in 1:L]

times   = collect(0.0:dt:t_max)
n_times = length(times)
χ_avg   = zeros(Float64, n_times)
χ2_avg  = zeros(Float64, n_times)
C_avg   = zeros(Float64, n_times)
C2_avg  = zeros(Float64, n_times)
ln_D    = log(Float64(D))

rng = MersenneTwister(rng_seed)
C_t = zeros(ComplexF64, D)   # pre-allocated work vectors to avoid per-step allocation
ψ_t = zeros(ComplexF64, D)

@printf "Sampling %d random Fock initial states (S_z=0)...\n" n_samples; flush(stdout)
for sample in 1:n_samples
    s0_idx = rand(rng, 1:D)

    # C0[k] = ⟨E_k|s0⟩ = conj(T_mat[s0_idx, k])
    # (T_mat row s0_idx = ⟨s0|E_k⟩ for all k; conjugate gives ⟨E_k|s0⟩)
    C0 = conj.(T_mat[s0_idx, :])

    χ_s = zeros(Float64, n_times)
    C_s = zeros(Float64, n_times)
    for (ti, t) in enumerate(times)
        @. C_t = C0 * exp((-im * t) * all_eigvals)
        mul!(ψ_t, T_mat, C_t)
        ipr     = sum(x -> abs2(x)^2, ψ_t)   # Σ_s |ψ_s(t)|^4
        χ_s[ti] = -log(ipr)
        prob    = abs2.(ψ_t)
        n_ev    = vec(occ'    * prob)          # ⟨n_k⟩, length L
        nn_ev   = vec(nn_occ' * prob)          # ⟨n_k n_{k+1}⟩, length L
        C_s[ti] = (4.0/L) * sum(nn_ev[k] - n_ev[k] * n_ev[k % L + 1] for k in 1:L)
    end

    χ_avg  .+= χ_s
    χ2_avg .+= χ_s .^ 2
    C_avg  .+= C_s
    C2_avg .+= C_s .^ 2
    @printf "  sample %2d:  χ(t_max)=%.4f  C(t_max)=%.4f  plateau=-1/(L-1)=%.4f\n" sample χ_s[end] C_s[end] (-1.0/(L-1))
    flush(stdout)
end

χ_avg ./= n_samples
χ_std   = sqrt.(max.(χ2_avg ./ n_samples .- χ_avg .^ 2, 0.0))
C_avg  ./= n_samples
C_std    = sqrt.(max.(C2_avg ./ n_samples .- C_avg .^ 2, 0.0))

# ── save ─────────────────────────────────────────────────────────────────────
mkpath("results")
outfile = "results/participation_entropy_L$(L)_tp$(@sprintf("%.3f",tp_hop))_Vp$(@sprintf("%.3f",Vp)).csv"
open(outfile, "w") do f
    println(f, "# χ_ent(t) = -ln(Σ_s |⟨s|ψ(t)⟩|^4)  Fock-basis IPR averaged over $n_samples random Fock initial states")
    println(f, "# C(t) = (4/L) Σ_k [⟨n_k n_{k+1}⟩ - ⟨n_k⟩⟨n_{k+1}⟩]  connected NN density correlator; T=∞ plateau = -1/(L-1) = $(-1.0/(L-1))")
    println(f, "# L=$L  N=$N  t_hop=$t_hop  V=$V  tp=$tp_hop  Vp=$Vp  rng_seed=$rng_seed")
    println(f, "# ln(D) = $ln_D  (ergodic ceiling: χ_ent → ln(D) - ln(2) ≈ $(ln_D - log(2.0)) for thermalizing system)")
    println(f, "time,chi_ent_mean,chi_ent_std,C_mean,C_std")
    for (ti, t) in enumerate(times)
        @printf f "%.4f,%.8f,%.8f,%.8f,%.8f\n" t χ_avg[ti] χ_std[ti] C_avg[ti] C_std[ti]
    end
end
@printf "\nln(D) = %.4f\nSaved → %s\n" ln_D outfile
