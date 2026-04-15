# GeoVac Atomic and Molecular Scope Boundary — v2.9.0

GeoVac's composed architecture handles any first-row molecule (atoms H through Ne) built from 1s² cores + valence blocks. The atomic classifier (`geovac/atomic_classifier.py`) maps Z to block decomposition for Z=1-10, and the general composed builder (`geovac/composed_qubit.py`) constructs qubit Hamiltonians from `MolecularSpec` dataclasses. Second-row atoms (Na-Ar) are feasible but require external core data that the framework does not currently compute. Transition metals are out of scope.

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

## Second Row (Z=11-18): Partially Feasible

| Atom | Config | Core | Valence | Status |
|:-----|:-------|:-----|:--------|:-------|
| Na (11) | [Ne] 3s¹ | 1s² 2s² 2p⁶ | 1e | Needs 10e core |
| Mg (12) | [Ne] 3s² | 1s² 2s² 2p⁶ | 2e | Needs 10e core |
| Al (13) | [Ne] 3s² 3p¹ | 10e | 3e | Needs 10e core |
| Si (14) | [Ne] 3s² 3p² | 10e | 4e | Needs 10e core |
| P (15)  | [Ne] 3s² 3p³ | 10e | 5e | Needs 10e core |
| S (16)  | [Ne] 3s² 3p⁴ | 10e | 6e | Needs 10e core |
| Cl (17) | [Ne] 3s² 3p⁵ | 10e | 7e | Needs 10e core |
| Ar (18) | [Ne] 3s² 3p⁶ | 10e | 8e | No molecular interest |

**The 10-electron core problem:** The first-row composed architecture uses a 2-electron (He-like) core solved by the Level 3 (hyperspherical) solver. Second-row atoms have a 10-electron [Ne] core: 1s² + 2s² + 2p⁶.

### Approach 1: Frozen-core tabulation (RECOMMENDED)

Solve the 10-electron Ne-like core once (using the full N-electron solver or external data), tabulate its properties (energy, density, Z_eff(r), PK parameters), and use these as inputs to the composed valence solver.

- **Pro:** Reuses existing valence infrastructure unchanged.
- **Con:** The 10-electron core solution itself is expensive. The current N-electron solver (Level 4N) handles 4 electrons at l_max=2 with difficulty (Track AK: 750 spectral dim). A 10-electron version would require SO(30) angular machinery.
- **Feasibility:** Yes, if core properties are tabulated externally (e.g., from NIST atomic data or Hartree-Fock calculations). The GeoVac valence solver only needs E_core, n_core(r), and the resulting Z_eff(r) and PK parameters — it does not need to solve the core quantum mechanically.

### Approach 2: Recursive composition

Treat the 10-electron core as nested composed blocks: 1s² (Level 3) + 2s² (Level 4 with PK from 1s² core) + 2p⁶ (three Level 4 pairs with PK from 1s²+2s²).

- **Pro:** Fully ab initio, no external data needed.
- **Con:** The 2p⁶ subshell has 3 pairs in 3 spatial orientations (m=-1,0,+1). These pairs interact strongly via exchange, which the current composed architecture handles poorly at high Z_eff (lone pair coupling is unphysical at Z_eff >= 6).
- **Feasibility:** Unlikely to be rigorous. The 2s-2p exchange coupling within the [Ne] core is a major correlation effect that recursive composition would struggle to capture.

### Approach 3: Effective core potential (ECP) from literature

Use published ECPs (e.g., Stuttgart/Cologne) for the [Ne] core and focus GeoVac on the valence electrons only.

- **Pro:** Standard approach in quantum chemistry; well-validated ECPs exist.
- **Con:** Introduces external parameters, losing the ab initio character.
- **Feasibility:** Straightforward if the goal is practical computation rather than foundational purity.

**Assessment:** Second-row atoms are feasible with approach (1) or (3) but require external core data. A fully ab initio treatment via recursive composition (approach 2) is unlikely to work without significant advances in the many-electron solver. Near-term recommendation: tabulate [Ne] core properties from NIST/HF data and use the existing composed valence infrastructure.

---

## Transition Metals (Z=21-30): Out of Scope

| Atom | Config | Core | Issue |
|:-----|:-------|:-----|:------|
| Sc (21) | [Ar] 3d¹ 4s² | 18e | 3d/4s near-degeneracy |
| Ti-Zn | [Ar] 3d^n 4s^m | 18e | Multi-reference, strong correlation |

Transition metals are out of scope for the foreseeable future, for four compounding reasons:

1. **18-electron core.** The [Ar] core requires solving (or tabulating) an 18-electron system. The recursive composition challenges from the 10-electron case are compounded.

2. **3d/4s near-degeneracy.** The Madelung rule predicts 4s fills before 3d, but 3d and 4s are nearly degenerate. PK cannot reproduce the s-d energy ordering (wrong sign, wrong magnitude — same structural issue as s-p splitting). The composed architecture would need explicit 3d-4s correlation, which is a multi-reference problem.

3. **Strong correlation.** Transition metal chemistry involves partially filled d-shells with strong electron correlation. This is one of the hardest problems in quantum chemistry and is far beyond the current composed architecture's capabilities.

4. **Spin-orbit coupling.** For heavier transition metals (Z > 30), spin-orbit effects become significant. The current framework is non-relativistic.

---

## Exotic Atoms (v2.9.0)

The hyperspherical framework (Level 3) extends to exotic two-particle systems via sign flips in the charge function. No new angular machinery is required; Gaunt selection rules are preserved.

| System | Particles | Charge function | Status | Error |
|:-------|:----------|:----------------|:-------|:------|
| He | e⁻e⁻ + Z=2 nucleus | -Z/sinα - Z/cosα + 1/r₁₂ | Production | 0.019% |
| H⁻ | e⁻e⁻ + Z=1 nucleus | Same as He at Z=1 | Tested (graph-native CI over-binds, standard CI works) | 2.1% (std FCI) |
| PsH | e⁻e⁺ + Z=1 proton | -1/sinα + 1/cosα - 1/r₁₂ | Prototype | 4.1% |
| Ps (positronium) | e⁻e⁺ (no nucleus) | Reduces to 1-body; graph identical to H | Trivial | exact |

**Key differences from standard He:**
- **PsH:** nuclear term sign-flipped for positron (+1/cosα), V_ee attractive (-1/r₁₂), alpha parity mixing required (distinguishable particles double the angular basis). Shallow well (0.042 Ha depth at l_max=3 vs He's 4.47 Ha).
- **H⁻:** same framework as He at Z=1. Graph-native CI violates variational bound (Z < Z_c ≈ 1.84); standard Casimir FCI is properly variational.

**Graph validity boundary:** The graph Laplacian CI is non-variational below Z_c ≈ 1.84, where the rigid kappa=-1/16 inter-shell coupling overestimates correlation for weakly-bound/asymmetric systems. Standard (non-graph) FCI is always variational. Details in CLAUDE.md Section 2.

---

## Near-Term Reachability Summary

| Category | Atoms | Status | Bottleneck |
|:---------|:------|:-------|:-----------|
| Fully operational | H, He, Li, Be, C, N, O, F | Production | Validation |
| Architecturally ready, untested | B, Ne | Need testing | PK validation, molecular targets |
| Feasible with external data | Na-Ar | Approach (1) or (3) | 10e core tabulation |
| Exotic atoms | PsH, H⁻ | Prototype (v2.9.0) | Asymptotic convergence (PsH), graph validity (H⁻) |
| Out of scope | Z > 20 | Fundamental | Multi-reference, many-electron core |
