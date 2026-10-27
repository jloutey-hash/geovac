# GeoVac: Structurally Sparse Qubit Hamiltonians from Graph Theory

![Status](https://img.shields.io/badge/Status-Research-blue) ![Version](https://img.shields.io/badge/Version-5.21.15-blue) ![License](https://img.shields.io/badge/License-MIT-orange) [![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.18738360-blue)](https://doi.org/10.5281/zenodo.18738360)

GeoVac constructs **structurally sparse qubit Hamiltonians** for molecular quantum simulation. The angular momentum selection rules of the hyperspherical harmonic basis enforce block-diagonal electron repulsion integrals, producing Hamiltonians whose Pauli count is **exactly linear in the qubit count across molecules at fixed basis** (N_Pauli = 27.90 × Q main-group, 30.03 × Q d-block) — **51× to 317× fewer Pauli terms** than Gaussian baselines at equal qubits for H₂O (n_max = 1–3; identity-included on both sides), and 76× fewer than cc-pVDZ for LiH (each at its own qubit count: Q = 30 vs 36) (Papers 14, 20). These comparisons are at matched qubit counts, not matched accuracy; the accuracy ceilings of the underlying basis are stated plainly in *Scope and Limitations* below and in [`docs/claims_register.md`](docs/claims_register.md).

```python
from geovac.ecosystem_export import hamiltonian
h = hamiltonian('LiH', R=3.015)  # 838 Pauli terms, 30 qubits
op_qiskit = h.to_qiskit()    # SparsePauliOp
op_of = h.to_openfermion()   # QubitOperator
op_pl = h.to_pennylane()     # qml.Hamiltonian
```

Install the standalone Hamiltonian package (not on PyPI; installs from this repo):
`pip install "geovac-hamiltonians @ git+https://github.com/jloutey-hash/geovac.git#subdirectory=geovac-hamiltonians"`

> **New here? Five documents, in order:** (1) [`papers/synthesis/geovac_field_guide.tex`](papers/synthesis/geovac_field_guide.tex) — what this project is, as a story; (2) [`docs/claims_register.md`](docs/claims_register.md) — every headline claim with its verification tier and falsifier; (3) [`papers/INDEX.md`](papers/INDEX.md) — status map of all ~60 papers (keystone / active / observation / descoped / historical); (4) Paper 7 (the S³ equivalence) and (5) Paper 14 (the headline computational result). The repository deliberately preserves its negative results and closed diagnostic arcs as institutional memory — `debug/archive/` and `papers/archive/` are the record, not the product.

> **Coming from Avery & Avery, *Hyperspherical Harmonics and Their Physical Applications* (World Scientific 2017)?** [`docs/outreach/avery_hyperspherical_map.pdf`](docs/outreach/avery_hyperspherical_map.pdf) maps what GeoVac builds on in that book, chapter by chapter: from the S³ label set through Fock's projection to generalized Sturmians, the Shibuya–Wulfman metric, and the quantum encoding of the isoenergetic secular equation.

---

## Headline Numbers

| Metric | Value |
|:-------|:------|
| Library | **37 systems** via the `hamiltonian()` API (35 composed molecules + He + H₂), three periodic-table rows |
| Pauli count | **N_Pauli = 27.90 × Q**, exactly linear across molecules at fixed basis (30.03 × Q d-block); within-molecule basis exponent 3.17 (small-basis fit; local slope ~3.8 by n_max = 4) |
| Advantage vs Gaussian | **51×–317×** fewer Pauli terms at equal qubits (H₂O, n_max = 1–3); **76×** vs cc-pVDZ (LiH) |
| LiH vs STO-3G | near parity on raw counts: 838 Pauli @ Q = 30 vs 907 raw (3.0× MORE than STO-3G's 276 reduced); 1-norm 34.0 vs 34.3 Ha |
| He accuracy | **0.022%** raw (2D variational, properly variational); **0.004%** after a non-variational cusp extrapolation; zero free parameters |
| H₂ accuracy | **99.97% D_e** (0.053 mHa; prolate explicit-r₁₂ CI, exact integrals, Paper 12); 96.0% via the Level-4 hyperspherical treatment (Paper 15) |
| Algebraic integrals | Exact hypergeometric R^k evaluator eliminates grid quadrature for Slater integrals |
| Mathematical foundation | Fock 1935 S³ conformal equivalence, **18 symbolic checks of the continuum chain** (ten limiting-case or definitional) |

---

## Why This Matters for Quantum Computing

The Pauli term count and 1-norm are the dominant cost factors for near-term (VQE/NISQ) and fault-tolerant (QPE) quantum simulation respectively. GeoVac's basis-intrinsic sparsity — from Gaunt selection rules, not post-hoc optimization — produces qubit Hamiltonians that are structurally cheaper to simulate than Gaussian-basis alternatives at the same qubit count (cheaper per qubit, not more accurate per qubit — see *Scope and Limitations*). The sparsity is compatible with all downstream optimizations (tapering, grouping, tensor factorization). See Paper 14 for the encoding theory and Paper 20 for resource benchmarks.

**Universal vs Coulomb-specific (Paper 22):** The angular sparsity guarantees are universal across spherical fermion systems — they hold for Coulomb, harmonic oscillator, Woods-Saxon, and any other radial potential, depending only on l_max. ERI density at l_max=3 is verified at 6.06% (global-M_L rule; the earlier 1.44% was the retired pair-diagonal rule) regardless of V(r). The S³ conformal projection and Hopf bundle structure are Coulomb-specific (by the Fock rigidity theorem, Paper 23). The sparsity extends to nuclear shell model Hamiltonians (Paper 23: deuteron 16 qubits / 688 Pauli, He-4 16 qubits / 828 Pauli); the conformal machinery does not.

**The HO has its own discretization (Paper 24):** The 3D harmonic oscillator has a discrete graph encoding on the holomorphic sector of S⁵ via the Bargmann-Segal transform, parallel to Fock's discretization of the Coulomb problem on S³. The Bargmann graph is **bit-exactly π-free** in exact rational arithmetic at every finite N_max. In the Kustaanheimo–Stiefel frame the Coulomb problem is itself a constrained 4D oscillator, and π enters both sides the same ways (a measure and spectral sums), so the two differ in representation theory and arithmetic, not in where π enters.

---

## The Math.OA / NCG / Periods Arc

Beyond the quantum-computing applications above, the framework's mathematical structure has been worked out as a discrete almost-commutative spectral triple in the Marcolli–van Suijlekom gauge-network lineage (with the Perez-Sanchez 2024/2025 correction). The math.OA-facing arc lives in `papers/group1_operator_algebras/` and its synthesis; highlights:

- **Paper 38** State-space Gromov–Hausdorff convergence of truncated Dirac triples on SU(2) (Paper 38, main theorem: Fejér-kernel bound with constant 4/π, a kernel constant; O(1/n) [MEASURED] with a tapered amplitude) — **unconditional for the commutative factor** on the true Camporesi–Higuchi window; the almost-commutative extension is open (2026-06-10: translation-seminorm metrization; kernel condition proved on the truthful substrate; dual direction via an exact-fit spinor lifted state)
- **Paper 45** Degeneracy theorem: the natural K⁺-compression route to a "Lorentzian propinquity" annihilates the spatial Dirac — the seminorm vanishes identically (bit-exact falsifier frozen in `tests/`). An earlier convergence claim through this device is **retracted** (2026-06-09; see the paper's History remark)
- **Papers 46–49** (Lorentzian metric-level claims built on the retracted quantity) are archived in `papers/archive/` (2026-09-14; see `docs/retired_papers.md`); the norm-resolvent arrow (Paper 47) and the cocycle-deficit algebra (Paper 49) survive
- **Paper 50** Re-derives the Klebanov–Pufu–Safdi continuum F-theorem values on S³ + S⁵ in closed Hurwitz form (continuum spectra; the truncation never enters)
- **Paper 56** Tannakian reconstruction closed at finite cutoff (5,864 bit-exact zero residuals); infinite-cutoff motivic-Galois identification open
- **Paper 55** A census of the S³ spectral-zeta and heat-trace transcendentals: they sit in pure-Tate / level-4 cyclotomic mixed-Tate rings (chemistry seeds, spectral determinants and Paper 61's elliptic period lie outside the census)

The "C-arc closure" sprints of June 2026 placed the Forced-Count dimension count, non-selection results and Observations in Paper 32 §VIII (tally there) — formal characterizations of where the structural skeleton ends and calibration data begins. Reader's on-ramp: `papers/synthesis/geovac_field_guide.tex`.

### Claims and verification status

**No result in this repository is externally peer-reviewed.** Every headline claim, its verification tier (symbolic proof / measured / panel-verified / internal theorem / conditional / observation / retracted), and its falsifier is catalogued in **[`docs/claims_register.md`](docs/claims_register.md)** — start there before citing anything. Project-internal vocabulary is mapped to standard terminology in [`docs/vocabulary_translation.md`](docs/vocabulary_translation.md). Two self-contained front-door notes for external readers: [`docs/outreach/`](docs/outreach/) (a 3-page math.OA summary and a 2-page periods note with ten-minute-checkable identities).

**Reading paths.** Quantum computing / chemistry: Papers 14 → 22 → 20, then `Quick Start` below. Operator algebras / NCG: `docs/outreach/note_n1_su2_truncations.pdf` → Papers 38 → 45 → 32. Periods / number theory: `docs/outreach/note_n2_s3_identities.pdf` → Papers 28 → 55. Project identity and honest scope: the field guide, then `docs/claims_register.md`.

---

## Quantum Encoding Summary

| Metric | GeoVac (Composed) | Gaussian Baselines | Advantage |
|--------|:-----------------:|:------------------:|:---------:|
| Pauli count across molecules (fixed basis) | 27.90 × Q, exactly linear | O(Q^3.9–4.3) | Structural |
| H₂O, equal qubits (n_max = 1–3) | — | — | 51×–317× |
| LiH vs cc-pVDZ | 838 terms | 63,519 | 76× |
| Atomic path (He), equal qubits | — | — | 0.56× at Q = 10 (unfavorable) to 1.5× at Q = 28; 1-norm 4.0×/7.3× lower, at a 55–58 mHa worse He energy |

### Market Test: GeoVac LiH vs Minimal-Basis STO-3G

| Metric | GeoVac LiH (Q=30) | Gaussian STO-3G | Comparison |
|--------|:-----------------:|:---------------------:|:---------:|
| Pauli terms | 838 | 907 raw / 276 reduced | near parity raw; 3.0× more than reduced |
| 1-norm (λ) | 34.0 Ha | 34.3 Ha | 0.99× (match) |

The decisive wins are the exact Q-linearity across molecules and the 76× reduction vs cc-pVDZ, not the comparison with minimal-basis STO-3G (all counts under the exact global-M_L rule, 2026-08-29).

---

## Molecule Library (37 systems)

| Row | Molecules |
|:---:|:----------|
| 1st (Z=1-10) | H₂, He, LiH, BeH₂, H₂O, HF, NH₃, CH₄ |
| 2nd (Z=11-18) | NaH, MgH₂, HCl, H₂S, PH₃, SiH₄ |
| 3rd (Z=19-36) | KH, CaH₂, GeH₄, AsH₃, H₂Se, HBr |
| Heavier hydrides (frozen cores) | SrH, BaH |
| Multi-center | LiF, CO, N₂, F₂, NaCl |
| Transition metals | ScH, TiH, VH, CrH, MnH, FeH, CoH, NiH, CuH, ZnH |

Isostructural invariance: molecules with the same block topology produce identical Pauli counts regardless of atomic species. The composed count is N_Pauli = 27.90 × Q (main-group) and 30.03 × Q (d-block), exactly linear across molecules at fixed basis. Current per-molecule counts: `hamiltonian(name).n_terms` and Papers 14/20.

### Nuclear Systems (Paper 23)

| System | Qubits | Pauli terms | 1-norm | Notes |
|:-------|:------:|:-----------:|:------:|:------|
| Deuteron (1p+1n) | 16 | 688 | see Paper 23 | Minnesota potential, two-species JW |
| He-4 (2p+2n) | 16 | 828 | see Paper 23 | Same Q, 12.25× Hilbert space |
| Composed deuterium (1p+1n+1e) | 26 | 710 | see Paper 23 | Nuclear + electronic block architecture |

(Counts corrected 2026-08-22, v5.0.0; the earlier 592 / 712 / 614 carried a Moshinsky N_tot guard.)

---

## Classical Validation Benchmarks

| System | Method | Result | Paper |
|:-------|:-------|:-------|:-----:|
| He (2e) | 2D variational (raw) / + cusp extrapolation | **0.022%** / **0.004%** (non-var.) error | 13 |
| He (2e) | Fixed-k = Z hydrogenic CI, physical Hamiltonian (n_max = 7) | **1.88%** (zero parameters, exact integrals; above the HF limit) | 13 |
| H₂ (2e) | Prolate explicit-r₁₂ CI, exact integrals | **99.97% D_e** (0.053 mHa) | 12 |
| H₂ (2e) | Mol-frame hyperspherical | **96.0% D_e** | 15 |
| HeH⁺ (2e) | Prolate explicit-r₁₂ CI | **R_eq −0.5%** at l = 3; energy ~1 mHa from the Kolos–Peek reference (the hyperspherical 93.1% is an adiabatic, non-variational over-binding; the consistent variational value is ~5% of D_e, Paper 15) | 12 |
| LiH (4e) | Prolate all-electron, analytic core (Route C) | **R_eq +0.2%**; energy −8.0490 Ha (analytic explicit-r₁₂) | 19 |
| LiH (4e) | Composed Level 3+4, ab initio PK | **R_eq 5.3%** | 17 |
| BeH₂ (6e) | Composed + exchange coupling | **R_eq 11.7%** | 17 |
| H₂O (10e) | Composed 5-block | **R_eq 19.4%** | 17 |
| LiH (4e) | Balanced coupled, n_max=3 | **0.20% energy** | 19 |

---

## Development Methodology

GeoVac is developed using an **AI-augmented agentic research workflow**. The principal investigator (J. Loutey) provides scientific direction, physical intuition, and quality control. Implementation, numerical exploration, and documentation drafting are performed collaboratively with large language models (Anthropic Claude). All physics results are validated against known analytical solutions, NIST reference data, and a symbolic suite (18 symbolic checks of Fock's S³ conformal chain, ten of them limiting-case or definitional). No result is accepted on the basis of LLM output alone.

This workflow is itself a research contribution — an experiment in whether agentic AI tools can accelerate independent scientific research. The project has no institutional affiliation. Primary dissemination is via GitHub and Zenodo (DOI-stamped releases).

---

## Quick Start

### Installation
```bash
git clone https://github.com/jloutey-hash/geovac.git
cd geovac
pip install -e .
```

### Qubit Hamiltonians (Main API)
```python
from geovac.ecosystem_export import hamiltonian

# Any of the 37 library systems
h = hamiltonian('LiH', R=3.015)
print(f"Qubits: {h.n_qubits}, Pauli terms: {h.n_terms}")
op = h.to_qiskit()  # or .to_openfermion(), .to_pennylane()

# Transition metal hydrides
h_fe = hamiltonian('FeH')
print(f"FeH: Q={h_fe.n_qubits}, Pauli={h_fe.n_terms}")
```

### Atomic Calculations
```python
from geovac import LatticeIndex

# Helium Full CI
idx = LatticeIndex(n_electrons=2, max_n=5, nuclear_charge=2,
                   vee_method='slater_full')  # physical one-body h1 (the default)
E, psi = idx.compute_ground_state()
print(f"He FCI: {E[0]:.6f} Ha")  # accuracy vs cutoff: see Paper 13
```

### Precision He (2D Variational)
```python
from geovac.level3_variational import solve_he_variational_2d

result = solve_he_variational_2d(Z=2, n_basis_R=25, n_basis_alpha=40, l_max=4)
print(f"He: {result['energies'][0]:.6f} Ha, error: {result['error_pct']:.4f}%")
```

---

## Paper Series

| # | Title | Key Result |
|:-:|-------|------------|
| 0 | Geometric Packing | Packing construction of the (n,l,m) label set; κ = −1/16 matched (an Observation) |
| 1 | Spectral Graph Theory | Eigenvalue methods, O(N) scaling |
| **7** | **Dimensionless Vacuum** | **Fock S³ chain (18 symbolic checks, ten limiting-case or definitional), SO(3N) generalization** |
| **14** | **Qubit Hamiltonians** | **N_Pauli = 27.90 × Q exactly linear; 51×–317× vs Gaussian at equal qubits (H₂O)** |
| **16** | **Chemical Periodicity** | **S_N representation theory, atomic classifier** |
| 6 | Quantum Dynamics (archived) | Rabi, spectroscopy, AIMD at O(V) |
| 11 | Molecular Fock Projection | Prolate spheroidal lattice, H₂⁺ to machine precision |
| **12** | Algebraic V_ee | Neumann expansion; H₂ 99.97% D_e (explicit r₁₂, 0.053 mHa), 99.81% re-based CI, 92.4% σ-only |
| **13** | Hyperspherical Lattice | He 0.022% (variational) / 0.004% (cusp extrapolation, non-variational), fiber bundle, algebraic structure |
| **15** | Level 4 Geometry | H₂ 96.0% D_e (structural N-electron treatment); HeH⁺ adiabatic 93.1% is non-variational (variational ~5%) |
| **17** | Composed Geometries | LiH 5.3%, BeH₂ 11.7%, ab initio PK |
| **18** | Exchange Constants | Weyl-Selberg taxonomy, π-free graph principle |
| **22** | **Angular Sparsity Theorem** | **Potential-independent ERI density (6.06% at l_max=3, global-M_L rule)** |
| **24** | **Bargmann-Segal Lattice** | **π-free HO discretization on S⁵ Hardy space** |
| FCI-A | Full CI (Atoms) | Atomic FCI benchmarks (current values in the paper) |
| **19** | Balanced Coupled | 0.20% energy, PK-free, 3-molecule census |
| **20** | Resource Benchmarks | 37-system library, Gaussian comparison, `pip install` |
| **23** | **Nuclear Shell Hamiltonians** | **Deuteron/He-4 qubit encoding, Fock rigidity theorem** |
| 21 | Synthesis (archived) | S³ proof chain, exchange constants, research program |
| **2** | **Fine Structure Constant** | **K = π(B + F − Δ) matches 1/α to 4.8×10⁻⁷ (its cubic's root to 8.8×10⁻⁸); three-sector spectral coincidence; Marcolli-vS gauge-network lineage; the combination rule is an Observation, not derived** |
| 8-9 | Bond Sphere + Sturmian | Structural theorem (guardrail), SO(4) selection rules |
| 25, 28, 30 | Hopf gauge / QED on S³ / SU(2) Wilson | Discrete lattice-gauge structure on the Hopf graph |
| 29 | Ramanujan Hopf | Graph-RH for GeoVac Hopf graphs; Ihara zeta closed forms |
| **32** | **Spectral Triple** | **Explicit construction; §VIII non-selection results and Observations (tally in the paper)** |
| 33, 36 | QED selection rules / Bound-state QED | 1+6+1 partition; Lamb shift at −0.534% one-loop |
| 34, 35 | Projection taxonomy / Time as projection | Two-layer framework; 28 projections × three-axis tagging |
| **38** | **SU(2) state-space GH convergence** | **Unconditional for the commutative factor; compression/lifted-state proof, translation seminorm (see claims register)** |
| 39, 40 | Tensor / Unified convergence | Same-manifold tensor case; 4/π Fejér-moment constant reported universal over compact simple groups (derived at rank 1; fit-sensitive at rank ≥ 2; general-G convergence conditional) |
| 42 | Modular Hamiltonian | Finite-cutoff modular structure; Krein-(3,1) carrier in §11 (Paper 43 archived into it) |
| **45** | **Lorentzian degeneracy theorem** | **K⁺ compression annihilates the spatial Dirac; earlier convergence claim retracted 2026-06-09** |
| 46–49 | Lorentzian metric arc (archived) | Archived 2026-09-14; norm-resolvent arrow + cocycle-deficit algebra survive |
| 50 | CFT₃ partition function | Closed-form re-derivation of the continuum KPS F-theorem values |
| 51, 53 | Gravity arc / flat-plane truncations | S³ spectral action two-term-exact; flat-plane state-space GH theorem, disk-with-cone open |
| **55** | **Periods of GeoVac** | **Cyclotomic mixed-Tate at level 4 (a census of S³ spectral data); master Mellin engine census** |
| **56** | **Tannakian Substrate** | **Finite-cutoff Tannakian reconstruction (5,864 bit-exact residuals); infinite-cutoff motivic-Galois identification open** |
| 57 | Forced/free seam | 77-entry catalogue (43 forced / 4 Observation / 1 admitted / 26 calibration / 3 walls; released-axis sub-domain J added 2026-10-15, J10 2026-10-16, J11 2026-10-17, J12–J14 2026-10-18, door-indexed in `docs/atlas.md`); P5 packing-reachability is an internal-consistency check, not a discriminator |

---

## Scope and Limitations

### What GeoVac Does Well
- **Structurally sparse qubit Hamiltonians:** Pauli count exactly linear in Q across molecules (27.90 × Q), 51×–317× fewer terms than Gaussian at equal qubits (H₂O)
- **Block-diagonal ERIs:** Gaunt selection rules enforce basis-intrinsic sparsity
- **Classical benchmarks:** H₂ 99.97% D_e (Paper 12), LiH R_eq +0.2% (Route C, Paper 19), He 0.022% raw / 0.004% (non-variational cusp extrapolation)
- **Algebraic Slater integrals:** Exact hypergeometric R^k evaluator eliminates grid quadrature (machine precision, 8x speedup)
- **Zero molecular fitting parameters** — all from nuclear charges and geometry

### Current Limitations
- **Classical PES accuracy:** R_eq errors of 5.3-19.4% for composed molecules (PK structural overcounting). Optimal for fixed-geometry quantum simulation, not geometry optimization.
- **Transition metals:** Full first series implemented as hydrides (v2.8.0); non-hydride TM molecules not yet built
- **Polyatomic accuracy:** H₂O R_eq 19.4% (composed, zero parameters)

### What GeoVac Does NOT Replace
- Production quantum chemistry for general molecules
- Gaussian basis integral technology for large systems
- Density functional theory for extended systems

---

## Citation

```
@software{geovac2026,
  author    = {J. Loutey},
  title     = {GeoVac: Structurally Sparse Qubit Hamiltonians from Graph Theory},
  year      = {2026},
  version   = {5.20.2},
  doi       = {10.5281/zenodo.18738360},
  url       = {https://github.com/jloutey-hash/geovac}
}
```

If adapting the acknowledgment for academic papers:

> Computational implementation and documentation drafting were performed with the assistance of large language models (Anthropic Claude) under the author's scientific direction. All results were validated against analytical benchmarks and symbolic check suites.

---

## License

MIT License - See [LICENSE](LICENSE) for details.

**Contact:** Issues and contributions welcome at [https://github.com/jloutey-hash/geovac/issues](https://github.com/jloutey-hash/geovac/issues)
