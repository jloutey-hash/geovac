# GeoVac Atomic and Molecular Scope Boundary — v2.1.0

GeoVac's composed architecture handles any first-row molecule (atoms H through Ne) built from 1s² cores + valence blocks, and second-row molecules (Na through Ar) via frozen-core [Ne] treatment. The atomic classifier (`geovac/atomic_classifier.py`) maps Z to block decomposition for Z=1-18, and the general composed builder (`geovac/composed_qubit.py`) constructs qubit Hamiltonians from `MolecularSpec` dataclasses. First-row cores (He-like) are solved by the Level 3 hyperspherical solver; second-row cores ([Ne]-like) use analytical Clementi-Raimondi Slater orbital exponents via `geovac/neon_core.py`. Transition metals are out of scope.

---

## First Row (Z=1-10): Fully Supported

| Z | Element | Structure Type | PK Source | Tested in Molecule? |
|:-:|:--------|:---------------|:----------|:--------------------|
| 1 | H | A (bare) | N/A | Yes (H₂, LiH, BeH₂, H₂O, HF, NH₃, CH₄) |
| 2 | He | B (closed shell) | N/A | Yes (standalone) |
| 3 | Li | C (core + 1 valence) | computed (Paper 17) | Yes (LiH) |
| 4 | Be | D (core + 2 valence) | computed (Paper 17) | Yes (BeH₂) |
| 5 | B | D (core + 3 valence) | computed (Track BI) | No |
| 6 | C | D (core + 4 valence) | computed (Track BI) | Yes (CH₄) |
| 7 | N | D (core + 5 valence) | computed (Track BI) | Yes (NH₃) |
| 8 | O | E (core + valence + lone pairs) | Z² scaled | Yes (H₂O) |
| 9 | F | E (core + valence + lone pairs) | computed (Track BI) | Yes (HF) |
| 10 | Ne | B (closed shell, 10e) | Z² scaled | No (no molecular interest) |

**Core treatment:** All first-row atoms Li-Ne have a He-like 1s² core solved by the Level 3 (hyperspherical) solver. PK parameters are derived ab initio from the He-like core solution (Paper 17).

**Valence treatment:** The composed architecture groups valence electrons into pairs (bond pairs, lone pairs). Each pair is treated as a Level 4 block with Z_eff screening. The number of valence blocks scales linearly with valence electron count.

---

## Known Limitations Within First Row

### PK does not reproduce s-p orbital splitting (Track BK: wrong sign)

The l-dependent PK pseudopotential does not capture orbital filling rules (Madelung rule). PK is a repulsive barrier that pushes s-orbitals UP in energy, while in real atoms s-orbitals are MORE bound than p due to core penetration. The effect has the wrong sign and is orders of magnitude too large:

| Atom | Z_eff | GV gap (E(2p)-E(2s)) | HF gap | Sign |
|:-----|:-----:|:---------------------:|:------:|:----:|
| Li | 1 | -0.716 Ha | +0.067 Ha | WRONG |
| C | 4 | -67.77 Ha | +0.272 Ha | WRONG |
| N | 5 | -148.3 Ha | +0.372 Ha | WRONG |
| O | 6 | -284.2 Ha | +0.612 Ha | WRONG |

This is a structural impossibility, not a parameter tuning issue: PK enforces core-valence orthogonality by RAISING the energy of valence orbitals that overlap with the core. No adjustment of PK parameters can fix this — the sign is wrong by construction.

For first-row atoms, filling order is trivial (1s², 2s²2p^x) and is manually assigned in the atomic classifier. For second-row and beyond, r-dependent Z_eff(r) from `algebraic_zeff.py` would be needed. This is an identified research direction, not a current capability.

### Lone-pair coupling is unphysical at Z_eff >= 6

For H₂O and HF, the Slater F⁰ integrals produce unphysical coupling magnitudes for lone pairs: S·F⁰ gives -28 Ha bond-lone and -15 Ha lone-lone coupling for H₂O, exceeding total electronic energy. Bond-bond coupling (~0.5 Ha) is validated and physical. Lone-pair coupling is disabled in production.

### Z² PK scaling is inaccurate

Z²-scaled PK parameters (extrapolated from Li²⁺) have 5-26% errors compared to ab initio values computed via direct CoreScreening solutions (Track BI, v2.0.30). Ab initio PK parameters from Track BI should be used for Z >= 5. The `atomic_classifier.py` module provides both Z²-scaled and ab initio values; ab initio is the production default where available.

---

## Isostructural Pauli Count Invariance (v2.0.30 observation)

The Pauli term count in the composed qubit Hamiltonian depends only on block topology (number and size of blocks), not on Z_eff or atomic identity. Isostructural molecules produce identical term counts:

| Block topology | Molecules | Pauli terms |
|:---------------|:----------|:-----------:|
| 1 core + 1 bond pair (Q=30) | LiH, HF | 334 |
| 1 core + 2 bond pairs (Q=50) | BeH₂, NH₃ | 556 |
| 1 core + 2 bond pairs + 2 lone pairs (Q=70) | H₂O, CH₄ | 778 |

This means the O(Q^2.5) Pauli scaling is a property of the composed architecture itself, independent of chemistry. The scaling exponent and sparsity advantage are structural invariants of the block decomposition.

---

## Second Row (Z=11-18): Implemented (Frozen-Core)

| Atom | Config | Core | Valence | Status |
|:-----|:-------|:-----|:--------|:-------|
| Na (11) | [Ne] 3s¹ | 1s² 2s² 2p⁶ | 1e | Implemented (NaH) |
| Mg (12) | [Ne] 3s² | 1s² 2s² 2p⁶ | 2e | Implemented (MgH₂) |
| Al (13) | [Ne] 3s² 3p¹ | 10e | 3e | Classifier ready, no molecule tested |
| Si (14) | [Ne] 3s² 3p² | 10e | 4e | Implemented (SiH₄) |
| P (15)  | [Ne] 3s² 3p³ | 10e | 5e | Implemented (PH₃) |
| S (16)  | [Ne] 3s² 3p⁴ | 10e | 6e | Implemented (H₂S) |
| Cl (17) | [Ne] 3s² 3p⁵ | 10e | 7e | Implemented (HCl) |
| Ar (18) | [Ne] 3s² 3p⁶ | 10e | 8e | Classifier ready, no molecular interest |

**Implementation (v2.1.0, Tracks CH-CL):** The frozen-core approach (Approach 1 below) has been implemented using Clementi-Raimondi Slater orbital exponents for the [Ne] core. The `FrozenCore(Z)` class in `geovac/neon_core.py` provides analytical Z_eff(r) screening from tabulated STO exponents, requiring no ab initio core solution. The atomic classifier (`geovac/atomic_classifier.py`) has been extended to Z=11-18 with frozen-core classification. Six second-row molecules have been built: NaH, MgH₂, HCl, H₂S, PH₃, SiH₄. All use the balanced coupled builder with cross-center V_ne from multipole expansion.

**Pauli scaling confirmed:** Q^2.50 across 4 second-row molecules, matching the first-row Q^2.5 exponent. Resource table at n_max=2:

| Molecule | Q | Pauli | 1-norm (Ha) | QWC |
|:---------|:--|:------|:------------|:----|
| NaH | 20 | 239 | 191 | 29 |
| MgH₂ | 40 | 1,501 | 303 | 109 |
| HCl | 50 | 2,936 | 1,169 | 133 |
| SiH₄ | 80 | 7,273 | 1,004 | 167 |

**Known limitations:**
- **Frozen-core approximation:** Core polarization neglected (standard assumption in quantum chemistry, but means second-row cores are less accurate than first-row Level 3-solved cores).
- **NaH accuracy:** 2-electron FCI PES shows overattraction at n_max=2 (no equilibrium found). This is consistent with known balanced coupled accuracy limitations at small basis sets (first-row LiH also has significant error at n_max=2).
- **n_max=3 blocked:** The Wigner D-matrix rotation in `shibuya_wulfman.py` currently handles l=0 and l=1 only. Second-row atoms with n_max=3 require l=2 rotation (Ivanic-Ruedenberg recursion), which has not been implemented.
- **1-norm values** include large E_core constants from the frozen-core energy; electronic-only 1-norms should be separated for fair comparison.

**343 new tests passing** (atomic classifier Z=11-18, FrozenCore screening, balanced coupled builders, ecosystem export).

### Design approaches considered

#### Approach 1: Frozen-core tabulation (IMPLEMENTED)

Uses Clementi-Raimondi Slater orbital exponents for analytical Z_eff(r) screening of the [Ne] core. No ab initio core solution required.

- **Pro:** Reuses existing valence infrastructure unchanged. Analytical Z_eff(r) is fast and well-characterized.
- **Con:** Core polarization neglected. Clementi-Raimondi exponents are empirically fitted (introduces external parameters).
- **Status:** Production implementation in `geovac/neon_core.py`.

#### Approach 2: Recursive composition

Treat the 10-electron core as nested composed blocks: 1s² (Level 3) + 2s² (Level 4 with PK from 1s² core) + 2p⁶ (three Level 4 pairs with PK from 1s²+2s²).

- **Pro:** Fully ab initio, no external data needed.
- **Con:** The 2p⁶ subshell has 3 pairs in 3 spatial orientations (m=-1,0,+1). These pairs interact strongly via exchange, which the current composed architecture handles poorly at high Z_eff (lone pair coupling is unphysical at Z_eff >= 6).
- **Feasibility:** Unlikely to be rigorous. The 2s-2p exchange coupling within the [Ne] core is a major correlation effect that recursive composition would struggle to capture.

#### Approach 3: Effective core potential (ECP) from literature

Use published ECPs (e.g., Stuttgart/Cologne) for the [Ne] core and focus GeoVac on the valence electrons only.

- **Pro:** Standard approach in quantum chemistry; well-validated ECPs exist.
- **Con:** Introduces external parameters, losing the ab initio character.
- **Feasibility:** Straightforward if the goal is practical computation rather than foundational purity.

---

## Transition Metals (Z=21-30): Scoping Demonstrated (v2.4.0, Track CZ/DA)

| Atom | Config | Core | Scoping Status |
|:-----|:-------|:-----|:---------------|
| Sc (21) | [Ar] 3d¹ 4s² | 18e | ScH Hamiltonian built: Q=30, 278 Pauli (composed) |
| Ti (22) | [Ar] 3d² 4s² | 18e | TiH Hamiltonian built: Q=30, 278 Pauli (composed) |
| V-Zn (23-30) | [Ar] 3d^n 4s^m | 18e | Not yet tested; same architecture should apply |

**Track CZ scoping result (POSITIVE):** d-orbital blocks are structurally feasible and actually CHEAPER per qubit than s/p blocks. Key measurements:

- d-only block (5 orbitals, Q=10): 56 Pauli terms, Pauli/Q = 5.60 (vs 11.20 for s+p block at same Q)
- d-d ERI density: 4.0% (vs 8.9% for s+p), due to more restrictive Gaunt selection rules at l=2
- Gaunt coefficients: only 20% nonzero for d-d interactions (k = 0, 2, 4)
- d-block Pauli count is 0.50x that of s+p block at same Q

**Track DA result (POSITIVE):** ScH and TiH transition metal hydrides produce valid qubit Hamiltonians:

| Molecule | Encoded e | Blocks | Q | Pauli (comp) | Pauli (bal) | Pauli/Q |
|:---------|:---------:|:------:|:-:|:------------:|:-----------:|:-------:|
| ScH | 3 | 2 | 30 | 278 | 574 | 9.27 |
| TiH | 4 | 2 | 30 | 278 | 574 | 9.27 |

Implementation: [Ar] frozen core (18e via FrozenCore), d-orbital block with `l_min=2` for isolated d-shell.

**Remaining limitations for production use:**

1. **3d/4s near-degeneracy.** The Madelung rule (4s before 3d) is handled by manual block assignment. PK cannot reproduce the s-d energy ordering. Balanced coupled (PK-free) bypasses this.

2. **Strong correlation.** Partially filled d-shells require FCI or equivalent correlated treatment. The quantum computer handles this — the classical solver cannot.

3. **Spin-orbit coupling.** For heavier transition metals, spin-orbit effects become significant. The current framework is non-relativistic.

4. **Atomic classifier.** Z=21-30 still raises `NotImplementedError` in `atomic_classifier.py`. Transition metal specs are built manually via `sch_spec()` / `tih_spec()` rather than through the general classifier.

---

## Near-Term Reachability Summary

| Category | Atoms | Status | Bottleneck |
|:---------|:------|:-------|:-----------|
| Fully operational (first row) | H, He, Li, Be, C, N, O, F | Production | Validation |
| Architecturally ready, untested | B, Ne | Need testing | PK validation, molecular targets |
| Implemented with frozen-core | Na, Mg, Si, P, S, Cl | Production (v2.1.0) | Accuracy at n_max=2; l=2 rotation for n_max=3 |
| Classifier ready, untested | Al, Ar | Need molecular targets | No molecules built yet |
| Scoping demonstrated | Sc (21), Ti (22) | Track CZ/DA (v2.4.0) | d-block feasible, Pauli/Q lower than s/p |
| Not yet tested | V-Zn (23-30) | Architecturally feasible | Same framework, manual spec required |
