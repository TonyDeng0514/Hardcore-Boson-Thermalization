# ED.jl — exact diagonalization for hardcore boson thermalization model
#
# Hamiltonian: NN/NNN hopping + interactions + Stark tilt (matches run_tebd.jl)
# Observable:  O = n_{j1} * n_{j2},  center bond (j1=L÷2, j2=L÷2+1)
# Saves (energy, <O>) for every eigenstate to results/ed_L{L}_N{N}.csv

using LinearAlgebra, SparseArrays, Printf

const L      = 18
const N      = 12

const t_hop  = 1.0
const V      = 1.0
const tp_hop = 0.0
const Vp     = 0.0
const μ      = 0.0
const α      = 0.0

# local 2×2 operators  (basis: |0⟩=empty, |1⟩=occupied)
const _bdag = sparse([2], [1], [1.0], 2, 2)   # b†
const _b    = sparse([1], [2], [1.0], 2, 2)   # b
const _n    = sparse([2], [2], [1.0], 2, 2)   # n = |1⟩⟨1|

# Embed A at site j, B at site k (j < k) in an L-site chain.
# Basis ordering: |σ_1,...,σ_L⟩ → 1-based index 1 + Σ_j σ_j·2^(L-j), σ_j ∈ {0,1}
function two_site_op(A, j, B, k, L)
    kron(
        kron(
            kron(
                kron(sparse(I, 2^(j-1), 2^(j-1)), A),
                sparse(I, 2^(k-j-1), 2^(k-j-1))
            ),
            B
        ),
        sparse(I, 2^(L-k), 2^(L-k))
    )
end

site_op(A, j, L) =
    kron(kron(sparse(I, 2^(j-1), 2^(j-1)), A), sparse(I, 2^(L-j), 2^(L-j)))

function build_hamiltonian_ed(L, t_hop, V, tp_hop, Vp, μ, α)
    dim = 2^L
    H = spzeros(Float64, dim, dim)

    for j in 1:(L-1)
        H += -t_hop * two_site_op(_bdag, j, _b,    j+1, L)
        H += -t_hop * two_site_op(_b,    j, _bdag, j+1, L)
        H +=      V * two_site_op(_n,    j, _n,    j+1, L)
    end

    for j in 1:(L-2)
        H += -tp_hop * two_site_op(_bdag, j, _b,    j+2, L)
        H += -tp_hop * two_site_op(_b,    j, _bdag, j+2, L)
        H +=      Vp * two_site_op(_n,    j, _n,    j+2, L)
    end

    for j in 1:L
        H += (-μ + α * j) * site_op(_n, j, L)
    end

    return H
end

# 1-based global indices of all basis states in the N-particle sector
function n_sector_indices(L, N)
    filter(1:2^L) do idx
        k = idx - 1
        count(j -> div(k, 2^(L-j)) % 2 == 1, 1:L) == N
    end
end

# Occupation (0.0 or 1.0) at site j for each state in the sector
function site_occupation(global_indices, j, L)
    map(global_indices) do idx
        k = idx - 1
        Float64(div(k, 2^(L-j)) % 2)
    end
end

# ─── sector ──────────────────────────────────────────────────────────────────
idx   = n_sector_indices(L, N)
dim_N = length(idx)
@printf "L=%d  N=%d  t=%.2f  V=%.2f  tp=%.2f  Vp=%.2f  μ=%.2f  α=%.2f\n" L N t_hop V tp_hop Vp μ α
@printf "N-sector dimension: %d\n" dim_N
@printf "Estimated memory (H + eigenvecs): ~%.1f MB\n\n" (2 * dim_N^2 * 8 / 1e6)

@printf "Building sparse Hamiltonian...\n"; flush(stdout)
H_sp = build_hamiltonian_ed(L, t_hop, V, tp_hop, Vp, μ, α)

@printf "Extracting N-sector and converting to dense...\n"; flush(stdout)
H_N = Symmetric(Matrix(H_sp[idx, idx]))

@printf "Diagonalizing (dim=%d)...\n" dim_N; flush(stdout)
vals, vecs = eigen(H_N)

# ─── observable: O = n_{j1} * n_{j2}, center bond ────────────────────────────
# O is diagonal in the occupation basis, so no matrix needed.
# ⟨E_n|O|E_n⟩ = vecs[:,n].² ⋅ diag_O
j1, j2 = L ÷ 2, L ÷ 2 + 1
diag_O = site_occupation(idx, j1, L) .* site_occupation(idx, j2, L)

@printf "Computing expectation values...\n"; flush(stdout)
obs_vals = [dot(vecs[:, n].^2, diag_O) for n in 1:dim_N]

# ─── save ─────────────────────────────────────────────────────────────────────
mkpath("results")
outfile = "results/ed_L$(L)_N$(N)_tp$(tp_hop)_Vp$(Vp).csv"
open(outfile, "w") do f
    println(f, "# ED: O = n_$(j1)*n_$(j2),  L=$(L)  N=$(N)  t=$(t_hop)  V=$(V)  tp=$(tp_hop)  Vp=$(Vp)  mu=$(μ)  alpha=$(α)")
    println(f, "energy,O_expval")
    for n in 1:dim_N
        @printf f "%.10f,%.10f\n" vals[n] obs_vals[n]
    end
end
@printf "Saved %d rows → %s\n" dim_N outfile
