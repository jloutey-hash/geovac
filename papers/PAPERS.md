# GeoVac Paper Series

**Last Updated:** March 15, 2026

## Reading Guide

The GeoVac papers form a narrative arc from atoms to molecules:

1. **Start here:** Paper 7 (the theoretical foundation -- graph Laplacian = S3 = Schrodinger)
2. **Atoms:** Papers 0, 1 (graph construction, eigenvalue methods), then FCI-A (multi-electron results)
3. **Multi-electron atoms:** Paper 13 (hyperspherical lattice, He at 0.05%, fiber bundle)
4. **Dynamics:** Paper 6 (time evolution, spectroscopy, AIMD on graph Hamiltonians)
5. **Molecules -- the problem:** Paper 8-9 (bond sphere geometry, why single-S3 fails for binding)
6. **Molecules -- the solution:** Paper 11 (prolate spheroidal lattice, H2+ with zero free params)
7. **Molecules -- two-electron:** Paper 12 (algebraic V_ee, Neumann expansion, H2 92.4% D_e)
8. **Molecules -- practical:** FCI-M (LCAO approach for LiH), Paper 10 (nuclear vibration/rotation)
9. **Ab initio spectroscopy:** Paper 13 Sec. IX (full pipeline: PES -> Morse -> nuclear lattice)

## Paper Inventory

### Core (`papers/core/`)

Defensible foundations. Strictly testable, all claims backed by numerical results.

| Paper | File | Status | Key Result |
|:------|:-----|:------:|:-----------|
| Paper 0 | `Paper_0_Geometric_Packing.tex` | Active | Universal constant K = -1/16 |
| Paper 1 | `paper_1_spectrum.tex` | Active | Spectral graph theory, O(N) eigenvalue methods. **v1.2.0: Berry phase k=2.113 retracted** (arg()=0 for real operators); Section IV rewritten with erratum |
| Paper 6 | `Paper_6_Quantum_Dynamics.tex` | Active | O(V) dynamics: Rabi, spectroscopy, AIMD |
| Paper 7 | `Paper_7_Dimensionless_Vacuum.tex` | Active | S3 proof (18/18 symbolic), Schrodinger recovery |
| Paper 10 | `paper_10_nuclear_lattice.tex` | Draft | Graph structures for molecular vibration/rotation |
| Paper 11 | `paper_11_prolate_spheroidal.tex` | Draft | Prolate spheroidal lattice: H2+ 0.21% R_eq, 0.70% E_min |
| FCI-A | `paper_fci_atoms.tex` | Draft | He 0.35%, Li 1.10%, Be 0.90% |
| Paper 12 | `paper_12_neumann_vee.tex` | Active | Neumann V_ee expansion: H2 92.4% D_e, cusp diagnosis (7.6% gap) |
| Paper 13 | `paper_13_hyperspherical.tex` | Active | Hyperspherical lattice: He 0.05%, fiber bundle, autoionization channels, ab initio spectroscopy |
| FCI-M | `paper_fci_molecules.tex` | Scaffold | LiH D_e 1.0% (CP-corrected), R_eq analysis |

### Methods (`papers/methods/`)

Methodological papers -- explorations that led to both positive and negative results.

| Paper | File | Status | Key Result |
|:------|:-----|:------:|:-----------|
| Paper 8-9 | `Paper_8_Bond_Sphere_Sturmian.tex` | Draft | Bond sphere geometry (positive), Sturmian structural theorem (negative), SO(4) selection rules |

### Conjectures (`papers/conjectures/`)

Theoretical physics explorations beyond computational quantum chemistry.

| Paper | Description |
|:------|:-----------|
| Paper 2 | Fine structure constant derivation (geometric ansatz) |
| Paper 3 | Holographic entropy, spectral dimension, central charge |
| Paper 4 | Mass-independence, universality, muonic hydrogen |
| Paper 5 | Comprehensive geometric vacuum framework (synthesis) |
| FAQ | Frequently asked questions |

### Archive (`papers/archive/`)

Superseded material preserved for reference.

| File | Origin |
|:-----|:-------|
| `papers_8_9_originals/` | Original separate Paper 8 and Paper 9 before merge |
| `diagnostic_arc_tables/` | Diagnostic arc details extracted from Paper 8-9 |

## Notes

- Papers 8 and 9 were merged on March 10, 2026. The merged paper moved to `methods/` on March 12, 2026 (v1.0.0).
- FCI papers renamed: `paper_geovac_fci.tex` -> `paper_fci_atoms.tex`, `paper_geovac_lcao_fci.tex` -> `paper_fci_molecules.tex`.
- Paper 11 created on March 12, 2026. Forward references added to Papers 7, 8-9, FCI-A, FCI-M, and 10.
- FCI-M paper scaffolded on March 11, 2026. Section structure and bibliography seeded; prose to be written.
- Papers 12 and 13 added in v1.1.0 (March 15, 2026).
- Paper 1 Berry phase section corrected in v1.2.0 (March 15, 2026): k=2.113 retracted, Section IV rewritten with erratum. See `debug/qa_sprint/berry_phase_reconciliation.md` for full analysis.
- Paper 13 updated in v1.2.0: new Sections X (autoionization) and XI (adiabatic limits).
