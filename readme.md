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

| Component | Status | Notes |
|---|---|---|
| TEBD time evolution | Done | `run_tebd.jl` + `gates.jl`; second-order Suzuki-Trotter |
| NN gates | Done | Absorbs on-site field `h_eff` with correct boundary weights |
| NNN gates | Done | Acts on sites `(j, j+2)` as two-site gates |
| Hamiltonian MPO | Done | `observable.jl`; used to measure initial energy `E0` |
| Density profile output | Done | Saves `n_profile_vs_t_L{L}_chi{chi}_alpha{alpha}.csv` to `results/` |
| Bond-dimension convergence | Done | `bond_convergence.jl`; CLI script, run as `julia bond_convergence.jl <chi>` |
| DMRG ground state | Done | `DMRG.jl`; standalone, not yet integrated with TEBD |
| Thermalization perturbation | In progress | Stark tilt (`α`) implemented; need to identify which parameters thermalize |
| EoS measurement | Not started | Planned via linear response or high-temperature expansion |

## How to Run

```bash
# Bond-dimension convergence study
julia bond_convergence.jl 256

# Direct TEBD from the REPL
julia -e 'include("run_tebd.jl"); run_tebd(20, 1.0, 1.0, 0.1, 0.1, 0.0; maxdim=64, ttotal=20.0)'
```

Results are written to `results/` as CSV files. The first row of each file contains the initial energy `E0` as a comment.

## Physics Notes

- Hardcore bosons are mapped to spin-1/2 via Holstein-Primakoff: `n_j = Sz_j + 1/2`.
- Simulations use half-filling; the initial state is a Néel state (alternating Up/Dn).
- Quantum number conservation (`conserve_qns=true`) is enforced throughout.
