# HO Two-Fermion Sprint Plan — Wiring Minnesota/Moshinsky onto Bargmann-Segal

**Purpose:** Execute EP-2's killer test. Compute two-fermion entanglement entropy on the Bargmann-Segal lattice (Paper 24) with a short-range two-body interaction, compare to He (two-electron Coulomb on S^3). Supports or refutes Paper 27's prediction that entropy signatures differ between Coulomb and HO.

## 1. Scope — existing modules

- `geovac/nuclear/bargmann_graph.py`: single-particle SU(3) lattice. Nodes `(N, l, m_l)` via `enumerate_nodes`, `shell_degeneracy`, `build_bargmann_graph`, `bargmann_ho_spectrum`, `verify_pi_free` — π-free, exact `Fraction`.
- `geovac/nuclear/moshinsky.py`: `moshinsky_bracket`, `_moshinsky_block`, `lab_to_relative_matrix_element`, `_moshinsky_1d`, `_bracket_cartesian`, `_coupled_cartesian_expansion`. Brackets derived combinatorially (rational, per Paper 23).
- `geovac/nuclear/minnesota.py`: `minnesota_params`, `minnesota_potential`, `ho_length_parameter`, `_ho_radial_wf`, `_gaussian_me_ho` (grid trapezoid — numerical), `minnesota_matrix_element_analytical` (currently wraps grid routine despite the name).
- `geovac/nuclear/nuclear_hamiltonian.py`: `SPState`, `DeuteronSpec`, `enumerate_sp_states`, `build_one_body`, `compute_tbme` (pn), `compute_same_species_tbme`, `_build_fci_matrix`, `_build_he4_fci_matrix`, `_enumerate_slater_dets`, `_sd_phase`, `diagonalize_deuteron`, `build_he4_hamiltonian`, `diagonalize_he4`, `hw_scan`.
- `geovac/casimir_ci.py`: `build_fci_matrix`, `build_graph_native_fci`, `compute_he_spectrum` — He on Bargmann/S^3 graph.
- `debug/entanglement_geometry.py`: `build_1rdm_from_singlet_ci`, `compute_entanglement_measures`, `assign_quantum_numbers_to_natural_orbitals`, `run_graph_native_analysis`.

## 2. Gap analysis

Paper 23's code uses `SPState(n_r,l,m_l,m_s)` — same quantum numbers as Bargmann `(N,l,m_l)` with `N=2n_r+l` plus spin. Mapping is identity (no rotation needed); just a bijection `SPState ↔ (N,l,m_l) × m_s`. No isomorphism machinery required — both are the Cartan/uncoupled SU(3) basis. Gap is plumbing: (a) a small adapter registering Bargmann nodes as the orbital sector of `SPState`; (b) a two-fermion (rather than p-n) driver that reuses `compute_same_species_tbme` with antisymmetry; (c) exposing the resulting 1-RDM to `entanglement_geometry.py`.

## 3. Implementation steps (ordered)

1. **Bargmann↔SPState adapter** (`bargmann_sp_adapter.py`). Emit `SPState` list from `enumerate_nodes(N_max)` × {↑,↓}, single species "fermion". Verify diagonal `build_one_body` matches `bargmann_ho_spectrum`. *Status: algebraic.*

2. **Analytical Talmi-integral Minnesota ME** (replace `_gaussian_me_ho` grid). Closed form: `<n l|e^(-κr²)|n' l>` in HO basis is a finite rational combination of Γ(l+3/2+k) with a single `√π` seed; equivalent to Talmi integrals `I_p = Γ(p+3/2)/((1+κb²)^(p+3/2))·(normalizations)`. Implement `_gaussian_me_ho_talmi` using Laguerre-expansion identity (Moshinsky-Brody). *Status: algebraic-pending → algebraic after this step; minimal transcendental seed = `√π` and `(1+κb²)^(1/2)` — per Paper 18 taxonomy.*

3. **Two-fermion TBME builder** (`two_fermion_hamiltonian.py`). Reuse `compute_same_species_tbme` (already does antisymmetrized Moshinsky + Minnesota). Expose `build_two_fermion_hamiltonian(N_max, hw)` returning sparse antisymmetrized FCI matrix on the `C(2N_sp, 2)` Slater-det basis via `_enumerate_slater_dets` and `_sd_phase_two_body`. *Status: algebraic (MT brackets rational, Minnesota MEs analytical after step 2).*

4. **Stretch — exact-rational path.** For a polynomial/separable test V (e.g. `V=λ(r_1·r_2)` or `λ r_{12}^2`), MEs are in ℚ; diagonalize with `sympy.Matrix` or `mpmath` to emit the ground state in exact rationals. Proof-of-principle algebraic two-fermion ground state. *Status: algebraic.*

5. **Sparse diagonalization.** Use `scipy.sparse.linalg.eigsh` (already pattern in `diagonalize_he4`). Block by conserved `M_L`, `M_S`, parity — mirrors `_build_he4_fci_matrix`. *Status: numerical-required (Lanczos).*

6. **1-RDM + entropy.** Feed ground state to `build_1rdm_from_singlet_ci` (generalize signature if needed to take arbitrary orbital basis list) → `compute_entanglement_measures` → `assign_quantum_numbers_to_natural_orbitals`. *Status: algebraic (1-RDM); numerical-required (log eigenvalues; `log` seed, single transcendental).*

7. **Driver + hw scan** parallel to `hw_scan` / `he4_hw_scan` to optimize `hw` and record S(hw).

## 4. Two-fermion ground state

Reuse `_build_he4_fci_matrix` machinery (already antisymmetrized, sparse, good-quantum-number blocked). Do **not** introduce new CI; the He-4 FCI routines are the right target — restrict to 2 particles, same species, Minnesota V. `casimir_ci.build_graph_native_fci` is the He/Coulomb analogue and defines the interface we should mirror for comparability.

## 5. Entanglement measurement interface

`build_1rdm_from_singlet_ci(ci_vector, slater_dets, n_orbitals) -> rho_1` then `compute_entanglement_measures(rho_1) -> {S_vN, S_Renyi2, occupations}`. Orbital labels via `assign_quantum_numbers_to_natural_orbitals(rho_1, bargmann_nodes)`. No changes required if step 1 adapter emits the same `(N,l,m_l,m_s)` tuples the EP-1 pipeline already consumes.

## 6. Validation targets (pre-entropy)

- SP spectrum: `build_one_body` eigenvalues vs `bargmann_ho_spectrum(N_max, hw)` — exact match required (π-free path).
- Talmi MEs: new analytical `_gaussian_me_ho_talmi` vs existing grid `_gaussian_me_ho` to 1e-8 (sanity).
- Deuteron reproduction: call `diagonalize_deuteron` via the adapter → reproduce Paper 23 Track NE deuteron binding at published `(N_shells, hw)`.
- He reference: `compute_he_spectrum` from `casimir_ci` at matching `N_max`, entropy via same EP-1 pipeline — establishes S_He baseline.

## 7. Risks / dead ends

- **Do NOT attempt S^3 conformal projection for HO.** Paper 23 Track NH Fock rigidity theorem: S^3 projection is unique to −Z/r (l-independence within n-shell). HO's shell rule is `N=2n_r+l` (l-dependent), so no S^3 map exists. Paper 24 dual (HO rigidity) says HO lives on Bargmann-Segal/CP² fiber, not S^3. Sprint computes HO entropy **natively on Bargmann-Segal** and compares to He entropy on S^3 — both are ground-state 1-RDM entropies; the geometries differ by design. This is the point of EP-2.
- Ordering: same-species antisymmetrization differs from p-n case — use `compute_same_species_tbme`, not `compute_tbme`.
- Sign conventions: Moshinsky bracket sign in `_bracket_cartesian` must match Paper 23; add regression test against published deuteron MEs.
- Do not reintroduce quadrature for MT brackets — `moshinsky_bracket` is already rational.

## 8. Exit criteria

Report numeric `S_HO(N_max, hw*)` and `S_He(N_max)` at matched `N_max ∈ {2,3,4}`, compute ratio `S_HO / S_He` and per-orbital entropy spectra. Sprint ends when:
- ratio is converged in `N_max` (change < 5% between 3 and 4), AND
- qualitative occupation spectrum (rate of Schmidt-coefficient decay) is tabulated.

Result supports Paper 27's prediction iff Coulomb shows slower-decaying Schmidt spectrum (higher S) than HO at matched `N_max` — the geometric signature of the 1/r long-range tail vs HO's Gaussian confinement.

## 9. Size estimate

**Medium** (3–5 sub-agents, ~1 day). Adapter + analytical Talmi + driver + validation + entropy sweep. Heaviest cost is step 2 (analytical Talmi) and regression against Paper 23 deuteron.

## 10. Papers updated on success

- **Paper 24** (Bargmann-Segal): add Section on two-fermion interacting Hamiltonian and entropy result.
- **Paper 23** (Nuclear shell): cross-reference the adapter and the exact-rational MT/Talmi pipeline; update algebraic registry entry for Minnesota MEs from "algebraic-pending" to "algebraic".
- **Paper 27** (Entanglement signatures / EP-2): primary host for the S_HO vs S_He killer-test result.
- **CLAUDE.md §12 Algebraic Registry**: flip Minnesota HO ME row to `algebraic`, record `√π` as sole transcendental seed.
- **Paper 18** (exchange-constant taxonomy): add HO-two-fermion entry.

## Algebraic opportunity summary

| Step | Status |
|---|---|
| 1 adapter | algebraic |
| 2 Talmi MEs | algebraic-pending → algebraic (seed: √π) |
| 3 TBME build | algebraic |
| 4 stretch rational V | algebraic (exact ℚ) |
| 5 Lanczos diag | numerical-required |
| 6 entropy | algebraic 1-RDM; `log` seed for S_vN |
| 7 hw scan | numerical-required (outer optimize) |

## Critical files for implementation

- `C:\Users\jlout\Desktop\Project_Geometric\geovac\nuclear\bargmann_graph.py`
- `C:\Users\jlout\Desktop\Project_Geometric\geovac\nuclear\moshinsky.py`
- `C:\Users\jlout\Desktop\Project_Geometric\geovac\nuclear\minnesota.py`
- `C:\Users\jlout\Desktop\Project_Geometric\geovac\nuclear\nuclear_hamiltonian.py`
- `C:\Users\jlout\Desktop\Project_Geometric\debug\entanglement_geometry.py`
