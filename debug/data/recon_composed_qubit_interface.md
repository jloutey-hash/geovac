# Reconnaissance: Composed Solver ↔ Qubit Encoder Interface Gap

**Date:** 2026-03-23
**Version:** v1.9.0
**Scope:** Map what `composed_diatomic.py` produces, what `qubit_encoding.py` consumes, and the gap between them.

---

## 1. What Does `qubit_encoding.py` Consume?

Two entry points, both requiring second-quantized integrals in an orbital basis:

### A. `JordanWignerEncoder(lattice_index)`

Expects a duck-typed object with:
- **`n_sp`** — number of spin-orbitals (= 2 × M spatial orbitals)
- **`_H1_spatial`** — sparse `csr_matrix` of shape (M, M), one-electron integrals in spatial orbital basis. Accessed via `.todense()`.
- **`_eri`** — `Dict[Tuple[int,int,int,int], float]`, two-electron integrals in **physicist notation**: `<ab|cd> = ∫ φ_a(1) φ_b(2) (1/r₁₂) φ_c(1) φ_d(2)`. Keys are spatial orbital indices.
- **`V_NN`** (optional) — nuclear repulsion energy (float).

### B. `build_fermion_op_from_integrals(h1, eri, nuclear_repulsion)`

Expects:
- **`h1`** — dense `np.ndarray` of shape `(M, M)`, one-electron integrals in spatial MO basis.
- **`eri`** — dense `np.ndarray` of shape `(M, M, M, M)`, two-electron integrals in **chemist notation**: `(pq|rs) = ∫ φ_p(1) φ_q(1) (1/r₁₂) φ_r(2) φ_s(2)`.
- **`nuclear_repulsion`** — float.

**Key difference:** The `JordanWignerEncoder` uses physicist notation `<ab|cd>` with a sparse dict; `build_fermion_op_from_integrals` uses chemist notation `(pq|rs)` with a dense 4D tensor. Both map spin indices as `sp = 2*spatial + sigma`.

The Gaussian reference pipeline (`gaussian_reference.py`) uses path B: it provides hardcoded `h1` and `eri` arrays (chemist notation), which `build_qubit_hamiltonian()` feeds to `build_fermion_op_from_integrals()` → `jordan_wigner()`.

---

## 2. What Does `composed_diatomic.py` Produce?

The composed solver is a **grid-based PDE pipeline**, not an orbital-basis CI solver. Its outputs are:

| Step | Method | Output |
|------|--------|--------|
| `solve_core()` | Hyperspherical radial solve (He-like) | `E_core` (scalar, Ha) |
| PK derivation | `AbInitioPK` from core density | `A`, `B` (Gaussian barrier params) |
| `_solve_valence_at_R(R)` | `solve_level4_h2_multichannel()` | `E_elec` (scalar, Ha) |
| `scan_pes()` | Loop over R grid | `E_composed(R)` array (PES curve) |
| `fit_spectroscopic_constants()` | Morse fit | `R_eq`, `D_e`, `ω_e`, `B_e` (spectroscopic constants) |
| `build_nuclear_lattice()` | SU(2) rovibrational | Rovibrational spectrum |

**Crucially:** At no point does the composed solver construct:
- One-electron integrals `h1[p,q]` in an orbital basis
- Two-electron integrals `(pq|rs)` or `<ab|cd>` in an orbital basis
- A second-quantized Hamiltonian
- Slater determinants or CI vectors

The Level 4 solver (`solve_level4_h2_multichannel`) works in **molecule-frame hyperspherical coordinates** `(R_e, α)` with angular channels `(l₁, l₂)`. It solves coupled PDEs on finite-difference grids. Its return dict contains `E_elec`, `E_total`, `D_e`, `wavefunction` (radial F(R_e)), and channel metadata — but no integrals in an orbital basis.

The core solver (`CoreScreening` → `hyperspherical_radial.solve_helium`) similarly works on `(R, α)` grids and returns energies and wavefunctions, not integrals.

---

## 3. What Is the Gap?

**Fundamental representation mismatch.** The two systems live in different mathematical worlds:

| Property | Composed Solver | Qubit Encoder |
|----------|----------------|---------------|
| **Representation** | Continuous wavefunctions on FD grids | Discrete orbital basis (second quantization) |
| **Coordinates** | Hyperspherical (R, α, l₁, l₂) | Spin-orbital indices (p, σ) |
| **Electron correlation** | Built into the PDE (adiabatic + multichannel) | Encoded in ERI tensor, solved by CI/VQE |
| **Core-valence coupling** | Z_eff(r) screening + PK barrier (grid potentials) | Would need explicit h1 + ERI cross-terms |
| **Output** | Energies (scalars), PES curves | Pauli string Hamiltonian (operator) |

There is **no direct path** from `ComposedDiatomicSolver.run_all()` → `JordanWignerEncoder`. The solver never constructs the objects the encoder needs.

---

## 4. What Orbital Basis Would Work?

For composed LiH with 2 core + 2 valence electrons (4 total):

### Option A: Full 4-electron encoding
- **Spatial orbitals:** Hydrogenic eigenstates up to `max_n`. At `max_n=3`: 1s, 2s, 2p₋₁, 2p₀, 2p₊₁, 3s, 3p₋₁, 3p₀, 3p₊₁, 3d₋₂, 3d₋₁, 3d₀, 3d₊₁, 3d₊₂ = 14 spatial orbitals per center.
- For LiH with two centers: could go up to ~28 spatial orbitals (56 qubits) but with many redundancies.
- **Problem:** This is the LCAO approach that v0.9.x already tried and is superseded.

### Option B: Natural-geometry orbitals (new approach)
- **Core:** 2 electrons in He-like hyperspherical states on atom A. These are naturally described by ~5 spatial functions (1s², with virtual excitations to 2s, 2p). The hyperspherical solver at `l_max=2` uses ~9 angular channels × `n_alpha` grid points — but the *effective* orbital space for the core ground state is small.
- **Valence:** 2 electrons in Level 4 molecule-frame hyperspherical channels. At `l_max=2, m_max=0`, there are 5 sigma channels: (0,0), (0,2), (2,0), (2,2), (1,1). Each channel has a radial function on the `R_e` grid.
- **Sensible basis size:** ~5 core orbitals + ~5-10 valence orbitals = 10-15 spatial orbitals → 20-30 qubits. This would be competitive with Gaussian STO-3G (12 qubits) or 6-31G (~20 qubits) for LiH.

### Option C: Effective 2-electron valence encoding (simplest)
- Freeze the core entirely (it's already solved to high accuracy).
- Encode only the 2 valence electrons in the Level 4 channel basis.
- ~5 spatial orbitals at `l_max=2` → 10 qubits.
- The PK pseudopotential and Z_eff screening become one-body effective potentials in the orbital basis.
- **This matches the composed solver's philosophy** and gives the smallest qubit count.

---

## 5. Where Do the Cross-Integrals Live?

The core-valence coupling in the composed solver has three components:

### 5a. Z_eff screening (one-body)
`CoreScreening` produces `Z_eff(r) = Z - N_core(r)` where `N_core(r)` is the integrated core density. This is a **radial potential** on a grid, not a matrix element. To express it as h1 integrals: `h1_screen[p,q] = ∫ φ_p(r) [-Z_eff(r)/r + Z_bare/r] φ_q(r) dr`. This is a straightforward one-electron integral that modifies the nuclear attraction.

**Expressible as one-body integrals?** Yes, if we have explicit orbital functions φ_p(r). The screening potential is smooth and well-behaved.

### 5b. Phillips-Kleinman pseudopotential (one-body)
`V_PK(r) = A * exp(-B*r²) / r²` is a **local potential** applied to the valence electrons. In the Level 4 solver, it enters as an additive term in the angular eigenvalue problem. In an orbital basis: `h1_PK[p,q] = ∫ φ_p(r) V_PK(r) φ_q(r) dr`.

**Expressible as one-body integrals?** Yes. It's a local multiplicative potential, so `<p|V_PK|q>` is well-defined. The Gaussian form makes analytical evaluation possible for Gaussian or hydrogenic orbitals.

### 5c. Core-nucleus cross-attraction (one-body)
`V_cross_nuc(R) = -n_core * Z_B * <1s_Z|1/r_B|1s_Z>` is already computed analytically in `_v_cross_nuc_1s()`. This is a **scalar** (R-dependent constant added to the PES), not an operator on the valence electrons. It would enter as part of the nuclear repulsion constant, not as h1 or ERI.

### 5d. Core-valence electron repulsion (two-body — currently absent!)
The composed solver does **not** compute explicit core-valence V_ee. The Z_eff screening is an *averaged* approximation: it replaces the two-body core-valence Coulomb interaction with an effective one-body potential. This is standard in pseudopotential theory but means:
- There are **no two-body cross-integrals** in the current framework.
- The core-valence correlation is entirely captured by Z_eff + PK.
- If we wanted explicit two-body terms, we'd need `<core_i val_j | 1/r₁₂ | core_k val_l>` integrals, which the hyperspherical solver doesn't produce.

**Bottom line:** All cross-terms are one-body (screening + PK) or scalar (cross-nuclear). No two-body cross-integrals exist in the current framework. This is physically justified by the core-valence separation assumption.

---

## 6. Minimal Code Change Sketch

The smallest path to qubit encoding for the composed LiH system:

### Step 1: Define the valence orbital basis
Extract the Level 4 channel functions as effective spatial orbitals. Each channel `(l₁, l₂)` at the equilibrium `R_e` defines a function of `α` (the hyperangle). Discretize to `n_channel` spatial orbitals. At `l_max=2`: 5 orbitals.

### Step 2: Compute h1 matrix elements in the channel basis
For each pair of channels `(c, c')`:
- **Kinetic:** Already computed by the angular solver's finite-difference Laplacian. Extract the kinetic matrix block.
- **Nuclear attraction:** `compute_nuclear_coupling()` already computes V_nuc coupling between channels. Extract at the equilibrium `R_e`.
- **PK pseudopotential:** The PK contribution to the angular Hamiltonian is already computed in `_build_angular_hamiltonian_mc()`. Extract the PK matrix block.
- **Z_eff correction:** The screened nuclear coupling uses `compute_nuclear_coupling_screened()`. This modifies the nuclear attraction matrix elements.

All of these are already computed internally by the Level 4 solver — they just need to be **exposed** rather than immediately consumed by the eigenvalue solver.

### Step 3: Compute ERI matrix elements in the channel basis
The electron-electron repulsion in the Level 4 solver is handled by the Gaunt integral coupling between channels. The angular Hamiltonian already includes V_ee coupling terms. These need to be separated out into a 4-index tensor `<c₁ c₂ | V_ee | c₃ c₄>`.

**This is the hard part.** The current solver computes V_ee as part of the coupled angular eigenvalue problem, mixing it with kinetic and nuclear terms. Separating them requires:
- Identifying which terms in `_build_angular_hamiltonian_mc()` come from V_ee vs V_nuc.
- Computing V_ee matrix elements between channel pairs independently.
- Storing them in the 4-index format the encoder expects.

### Step 4: Add nuclear repulsion and core energy
- `V_NN = Z_A_bare * Z_B / R` (scalar)
- `V_cross_nuc` (scalar, from `_v_cross_nuc_1s`)
- `E_core` (scalar, from hyperspherical solve)
- These enter as the `nuclear_repulsion` parameter.

### Step 5: Feed to `build_fermion_op_from_integrals()`
Assemble `h1` (5×5) and `eri` (5×5×5×5) from Steps 2-3, set `nuclear_repulsion = V_NN + V_cross_nuc + E_core`, and call the existing function.

### Estimated effort
- Steps 1-2: **Medium** — the matrix elements exist inside the solver, need refactoring to expose them.
- Step 3: **Hard** — V_ee coupling is interleaved with other terms in `_build_angular_hamiltonian_mc()`. Requires careful separation and validation.
- Steps 4-5: **Easy** — just arithmetic and a function call.

### Alternative: grid-to-orbital projection
Instead of extracting matrix elements from the PDE solver, compute orbitals on the grid and then evaluate integrals numerically:
1. Solve the Level 4 angular problem at fixed `R_e` to get channel eigenvectors.
2. Construct spatial orbital functions `φ_c(r₁, r₂)` from the channel wavefunctions.
3. Numerically integrate `<φ_a φ_b | 1/r₁₂ | φ_c φ_d>` on the grid.

This avoids refactoring the solver but requires a numerical integration engine for 6D integrals (expensive).

---

## 7. Gaussian LiH Comparison Data

### Available data

| System | Basis | Spatial Orbs | Qubits | Pauli Terms | Source |
|--------|-------|:---:|:---:|:---:|--------|
| LiH | STO-3G | 6 | 12 | ~631 | OpenFermion benchmark (widely published) |
| LiH | 6-31G | 11 | 22 | ~13,000 (est.) | Scaling from Q^3.56 atomic fit |
| LiH | cc-pVDZ | ~28 | ~56 | ~500,000 (est.) | Scaling from Q^4.60 molecular fit |

### LiH STO-3G details
LiH STO-3G is a standard benchmark in quantum computing literature:
- 6 spatial MOs: Li 1s, Li 2s, Li 2p_x, Li 2p_y, Li 2p_z, H 1s
- 12 spin-orbitals → 12 qubits (Jordan-Wigner)
- **631 Pauli terms** — this is a well-known number from OpenFermion/Qiskit tutorials
- FCI energy at R=3.015: approximately -7.882 Ha
- Widely used in VQE papers (Kandala et al., Nature 2017; etc.)

### Scaling estimates
From the validated Gaussian scaling in Paper 14:
- **Atomic systems (single-center):** Q^3.56 (from cc-pVDZ → cc-pVTZ He data)
- **Molecular systems (multi-center):** Q^4.60 (literature estimate, not yet validated with computed integrals for molecules)

For LiH comparison at matched qubit counts:
- GeoVac composed at ~10 qubits (5 valence channels): would need computation
- Gaussian STO-3G at 12 qubits: 631 Pauli terms
- The advantage ratio depends critically on how many GeoVac channels/orbitals are needed for comparable accuracy

### What's needed to make the comparison
1. **Compute** LiH STO-3G integrals (requires PySCF or manual entry from published data)
2. **Or** hardcode the known 631 Pauli term count as a reference point
3. For a fair comparison: match either qubit count or accuracy, not both
4. The GeoVac advantage claim rests on structural sparsity — fewer ERI nonzeros per orbital due to angular momentum selection rules. This should transfer to the composed basis if the channel functions have good angular momentum quantum numbers.

---

## 8. Summary and Recommendations

### The gap is real but bridgeable
The composed solver and qubit encoder operate in fundamentally different representations (grid PDE vs. orbital second-quantization). However:
- The Level 4 solver internally computes all the matrix elements needed (kinetic, nuclear, V_ee coupling between channels)
- They just need to be **extracted and reformatted**, not recomputed
- The core-valence coupling is entirely one-body (screening + PK), which simplifies the interface

### Recommended approach: Option C (valence-only encoding)
1. Freeze core (2e on Li, already solved)
2. Define ~5 valence channel orbitals from Level 4 at l_max=2
3. Extract h1 (5×5) and eri (5×5×5×5) from the angular Hamiltonian construction
4. Feed to `build_fermion_op_from_integrals()` with `nuclear_repulsion = V_NN + V_cross + E_core`
5. Compare 10-qubit GeoVac encoding against 12-qubit Gaussian STO-3G (631 Pauli terms)

### Key risk
The Level 4 "channels" are not standard spatial orbitals — they are correlated two-electron functions labeled by `(l₁, l₂)`. Mapping them to a single-particle orbital basis for second quantization requires care. Each channel describes a *pair* of electrons, not a single electron. The ERI tensor would encode how different channel assignments interact, which is conceptually different from standard orbital ERIs.

This may require defining *effective single-particle orbitals* from the hyperspherical coordinates (e.g., the `l₁` quantum number labels electron 1's angular momentum, so channel functions can be decomposed into single-particle components). This decomposition is physically natural but needs to be implemented and validated.
