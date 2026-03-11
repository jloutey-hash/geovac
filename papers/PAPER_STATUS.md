# GeoVac Paper Status

**Last Updated:** March 11, 2026

| Paper | File | Status | Description |
|:------|:-----|:------:|:------------|
| Paper 0 | `papers/core/Paper_0_Geometric_Packing.tex` | Active | Geometric packing framework and universal constant K = -1/16 |
| Paper 1 | `papers/core/paper_1_spectrum.tex` | Active | Spectral graph theory foundations and eigenvalue methods |
| Paper 6 | `papers/core/Paper_6_Quantum_Dynamics.tex` | Active | O(V) quantum dynamics: Rabi oscillations, delta-kick spectroscopy, AIMD |
| Paper 7 | `papers/core/Paper_7_Dimensionless_Vacuum.tex` | Active | Dimensionless vacuum: graph Laplacian as unit S³, Schrödinger recovery (18/18 proofs) |
| Paper 8-9 | `papers/core/Paper_8_Bond_Sphere_Sturmian.tex` | Draft | Bond Sphere: diatomic molecules as S³ with two weighted poles, Sturmian negative theorem, SO(4) selection rules |
| FCI Paper | `papers/core/paper_geovac_fci.tex` | Draft | Full CI on sparse graph Hamiltonian: multi-electron atoms and LiH via relational lattice indexing |
| LCAO FCI Paper | `papers/core/paper_geovac_lcao_fci.tex` | Scaffold | Topological LCAO FCI for heteronuclear diatomics: LiH benchmark, diagnostic arc, R_eq analysis |
| Paper 8 (archived) | `old_research_archive/papers_8_9_originals/Paper_8_Bond_Sphere.tex` | Archived | Original Paper 8 (Bond Sphere geometry only), superseded by merged Paper 8-9 |
| Paper 9 (archived) | `old_research_archive/papers_8_9_originals/Paper_9_Sturmian_Bond_Sphere.tex` | Archived | Original Paper 9 (Sturmian basis only), superseded by merged Paper 8-9 |

## Tier Classification

- **Core** (`papers/core/`): Defensible foundations. Strictly testable O(N) sparse graph Laplacian implementations.
- **Conjectures** (`papers/conjectures/`): Theoretical physics explorations (Papers 2-5, FAQ). Not listed here; see `papers/conjectures/` directory.

## Notes

- Papers 8 and 9 were merged on March 10, 2026. The merged paper contains all content from both originals plus nine new sections (Runge-Lenz alignment, heteronuclear egg geometry, kinetic membrane, atom-centered failure theorem, uniform diagonal theorem, and MO-FCI specification).
- Paper 8-9 updated on March 11, 2026 (v0.9.37) with three new results: Sturmian structural theorem (H proportional to S), SO(4) selection rules (D2_10_10=1, D2_00_10=0), and harmonic phase locking negative result. Abstract and conclusion revised.
- LCAO FCI paper scaffolded on March 11, 2026. Section structure and bibliography seeded; prose to be written in a future session.
- The FCI tech spec was updated on March 10, 2026 with the complete heteronuclear diagnostic arc table (v0.9.11--v0.9.30), prolate spheroidal implementation documentation, and v0.9.31 MO-FCI architecture specification.
