# Hardcore-Boson Thermalization

We study the Equations of State (EoS) for hardcore bosons in this project.

The hardcore-boson free model is exactly solvable, so the partition function and EoS are analytically known. What we are interested in is how one may verify this experimentally.

## Research Plan

1. Identify perturbations that thermalize the model (integrable models do not thermalize on their own).
2. Time-evolve the perturbed model for long times to obtain a thermal state.
3. Measure the EoS from the thermal state via linear response theory or high-temperature expansion.

## Target Hamiltonian

The base model is the hardcore-boson hopping model:

$$H = -t \sum_i \left( b^\dagger_i b_{i+1} + \text{h.c.} \right)$$

The perturbed model adds nearest-neighbor (NN) and next-nearest-neighbor (NNN) interactions, as well as a linear Stark tilt $\alpha$:

$$H = -t \sum_i \left( b^\dagger_i b_{i+1} + \text{h.c.} \right) - t' \sum_i \left( b^\dagger_i b_{i+2} + \text{h.c.} \right) + V \sum_i n_i n_{i+1} + V' \sum_i n_i n_{i+2} + \alpha \sum_i j \, n_j$$

## What Has Been Implemented

### TEBD (tensor network)

| Component | Status | Notes |
|---|---|---|
| TEBD time evolution | Done | `run_tebd.jl` + `gates.jl`; second-order Suzuki-Trotter |
| NN gates | Done | Absorbs on-site field `h_eff` with correct boundary weights |
| NNN gates | Done | Acts on sites `(j, j+2)` as two-site gates |
| Hamiltonian MPO | Done | `observable.jl`; used to measure initial energy `E0` |
| Density profile output | Done | Saves `n_profile_vs_t_L{L}_chi{chi}_alpha{alpha}.csv` to `results/` |
| Bond-dimension convergence | Done | `bond_convergence.jl`; CLI script, run as `julia bond_convergence.jl <chi>` |
| DMRG ground state | Done | `DMRG.jl`; standalone, not yet integrated with TEBD |

### Exact Diagonalization (ED)

| Component | Status | Notes |
|---|---|---|
| Full-space ED (OBC) | Done | `ED.jl`; sparse build + dense LAPACK diagonalization in the N-particle sector; supports Stark tilt `α` |
| Momentum-sector ED (PBC) | Done | `ed_momentum.jl`; block-diagonalizes `H` by translational symmetry into `L` complex Hermitian blocks; requires `α=0` |
| Level statistics | Done | Computed from momentum-sector ED eigenvalues; results consistent with ETH |
| Participation entropy χ_ent(t) | Done | `participation_entropy.jl`; Fock-basis IPR entropy averaged over random initial states; results for L=10–16 |
| Connected NN correlator C(t) | Done | Computed alongside χ_ent; plateau value `-1/(L-1)` serves as T=∞ reference |

### Thermalization diagnostics

| Component | Status | Notes |
|---|---|---|
| Thermalization perturbation (NNN) | Done | `tp`, `Vp` break integrability; level statistics and participation entropy confirm thermalization |
| Stark tilt (`α`) | Implemented | Available in TEBD and full-space ED; breaks translation invariance so not used in momentum-sector ED |
| EoS measurement | Not started | Planned via linear response or high-temperature expansion |

## File Map

### TEBD

| File | Role |
|---|---|
| `gates.jl` | Builds NN and NNN two-site Trotter gates; assembles symmetric gate sequence |
| `observable.jl` | Builds the Hamiltonian as an MPO (`build_hamiltonian_mpo`) |
| `run_tebd.jl` | Library: `run_tebd(L, t_hop, V, tp_hop, Vp, μ; ...)` — TEBD driver |
| `bond_convergence.jl` | Runner: bond-dimension convergence study; calls `run_tebd` |
| `DMRG.jl` | Standalone DMRG ground-state finder (not yet integrated with TEBD) |

### Exact Diagonalization

| File | Role |
|---|---|
| `ED.jl` | Library: `run_ed(L, N, t_hop, V, tp_hop, Vp, μ, α)` — full-space sparse ED in the N-particle sector; OBC |
| `run_ed.jl` | Runner for `ED.jl`; uses `julia -t auto` (BLAS threads matched to Julia threads) |
| `ed_momentum.jl` | Library: `run_ed_momentum(L, N, t_hop, V, tp_hop, Vp; j1, j2)` — PBC momentum-sector ED, block-diagonalizes by translation |
| `run_ed_momentum.jl` | Runner for `ed_momentum.jl` |
| `participation_entropy.jl` | Library: `run_participation_entropy(L, t_hop, V, tp_hop, Vp; n_samples, t_max, dt, rng_seed)` — builds full D×D Fock↔energy unitary, computes χ_ent(t) and C(t); includes `ed_momentum.jl` |
| `run_participation_entropy.jl` | Runner: sweeps L=10,12,14,16 for tp=0.3 and tp=0.0 (8 runs total) |

### Notebooks

| File | Role |
|---|---|
| `01_bond_convergence.ipynb` | Bond-dimension convergence analysis |
| `02_therm_time.ipynb` | Thermalization time analysis |
| `06_eth_plot.ipynb` | ETH diagnostic plots |

### Archive

| File | Role |
|---|---|
| `archive/` | Old monolithic TEBD prototype — superseded, kept for reference |

## How to Run

```bash
# Bond-dimension convergence study
julia bond_convergence.jl 256

# Direct TEBD from the REPL
julia -e 'include("run_tebd.jl"); run_tebd(20, 1.0, 1.0, 0.1, 0.1, 0.0; maxdim=64, ttotal=20.0)'

# Full-space ED (edit run_ed.jl for parameters, then)
julia -t auto run_ed.jl

# Momentum-sector ED (edit run_ed_momentum.jl for parameters, then)
julia -t auto run_ed_momentum.jl

# Participation entropy (edit run_participation_entropy.jl for parameters, then)
julia -t auto run_participation_entropy.jl
```

Results are written to `results/` as CSV files.

### Output CSV naming

| Script | Output | Columns |
|---|---|---|
| `run_tebd.jl` | `n_profile_vs_t_L{L}_chi{chi}_alpha{α}.csv` | First row is `# E0 = ...`; rest is density profile vs time |
| `ED.jl` | `ed_L{L}_N{N}_tp{tp}_Vp{Vp}.csv` | `energy, O_expval` |
| `ed_momentum.jl` | `ed_ksectors_L{L}_N{N}_tp{tp}_Vp{Vp}.csv` | `sector, energy, O_expval` |
| `participation_entropy.jl` | `participation_entropy_L{L}_tp{tp}_Vp{Vp}.csv` | `time, chi_ent_mean, chi_ent_std, C_mean, C_std` |

## Physics Notes

- Hardcore bosons are mapped to spin-1/2 via Holstein-Primakoff: `n_j = Sz_j + 1/2`.
- Simulations use half-filling (`N = L/2`); the initial state is a Néel state (alternating Up/Dn).
- Quantum number conservation (`conserve_qns=true`) is enforced throughout TEBD.
- Momentum-sector ED uses PBC and requires `α=0` (Stark tilt breaks translational invariance).
- Participation entropy ergodic ceiling: `χ_ent → ln(D) − ln(2)` for a thermalizing system, where `D = C(L, N)`.
- The connected NN density correlator `C(t) = (4/L) Σ_k [⟨n_k n_{k+1}⟩ − ⟨n_k⟩⟨n_{k+1}⟩]` has T=∞ plateau `-1/(L-1)`.
