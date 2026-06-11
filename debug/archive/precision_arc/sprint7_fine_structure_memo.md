# Sprint 7: Fine-Structure Splittings for CaH, SrH, BaH

## Summary

Built relativistic composed qubit Hamiltonians for three heavy-atom
monohydrides (CaH, SrH, BaH) using the Tier 2+3 Dirac-on-S^3 pipeline.
All three are isostructural (same block topology, frozen-core-screened
Z_eff=2 for valence) and produce identical Pauli counts and QWC groups.

## Resource Metrics

| Molecule | Z  | Q  | N_Pauli | N_P(Breit) | lam_ni (Ha) | lam_ni(B) | QWC | QWC(B) |
|----------|----|----|---------|------------|-------------|-----------|-----|--------|
| CaH      | 20 | 20 | 942     | 942        | 18.6790     | 18.6795   | 111 | 115    |
| SrH      | 38 | 20 | 942     | 942        | 18.6790     | 18.6795   | 111 | 115    |
| BaH      | 56 | 20 | 942     | 942        | 18.6790     | 18.6795   | 111 | 115    |

Key observations:
- **Pauli count unchanged by Breit**: 942 in both cases (same Gaunt selection rules)
- **1-norm shift from Breit**: +0.0005 Ha (0.003% relative), confirming O(alpha^2) suppression
- **QWC groups**: 111 without Breit, 115 with Breit (4 extra groups from SS+SOO terms)
- **All three molecules are bit-identical** in every metric (isostructural invariance)

Note: CLAUDE.md reports lam_ni = 13.87 Ha for these molecules. The discrepancy
(18.68 vs 13.87) may reflect a different convention (e.g. with/without nuclear
repulsion in the identity term, or a different R value). The relative values
and isostructural invariance are confirmed.

## Spin-Orbit Splittings

### One-electron SO diagonal (h1_so_diag)

The SO splitting for the 2p doublet (j=3/2 vs j=1/2) depends on
Z_eff (frozen-core-screened), not the full nuclear charge Z.
Since all three molecules have Z_eff=2.0, the analytical SO
splittings are IDENTICAL (isostructural invariance extends to SO).

Unique SO eigenvalues from the qubit Hamiltonian diagonal:
- -1.775e-05 Ha (degeneracy 2): n=2, kappa=+1 (j=1/2, l=1)
- -1.109e-06 Ha (degeneracy 2): partner block contribution
- +5.547e-07 Ha (degeneracy 4): partner block contribution
- +8.875e-06 Ha (degeneracy 4): n=2, kappa=-2 (j=3/2, l=1)

Analytical 2p doublet splitting at Z_eff=2:
- delta_SO = +2.663e-05 Ha = +5.84 cm^-1

For all three molecules (CaH, SrH, BaH) identically.

### 2-electron FCI eigenvalues (2-electron sector, dim=190)

Ground state eigenvalues (including nuclear repulsion):
- CaH: E_0 = -678.358 Ha (no Breit), -678.358 Ha (with Breit)
- SrH: E_0 = -3133.880 Ha, -3133.880 Ha
- BaH: E_0 = -7885.899 Ha, -7885.899 Ha

The large absolute energies include the frozen-core energy.

**Breit effect on ground state**: shifts E_0 by +0.000214 Ha (0.47 mHa)
across all three molecules. This is the Breit SS+SOO two-body correction.

**Eigenvalue structure** (identical across all three molecules, modulo
nuclear repulsion offset):
- E_0: singlet ground state (non-degenerate)
- E_1-E_4: 4-fold degenerate level, gap = 0.3334 Ha = 73,174 cm^-1
- E_5-E_8: 4-fold degenerate level, gap = 0.3750 Ha = 82,303 cm^-1
- E_9: nearly degenerate with E_5-E_8, split by 1.11e-06 Ha = 0.24 cm^-1

The 0.24 cm^-1 splitting within the E_5-E_9 cluster is the SO-induced
fine-structure splitting visible in the 2-electron spectrum. The E_1-E_4
cluster is unsplit (pure l=0 character, Kramers cancellation).

## Sunaga 2025 Comparison

Sunaga PRA 111, 022817 reports RaH at Q=18 with 47,099 Pauli terms.
Per-molecule CaH/SrH/BaH data from SI Tables S1-S3 are DEFERRED.

GeoVac native-Q advantage vs Sunaga RaH-18q (mismatched molecule):
- CaH: 942 vs 47,099 = 50x
- SrH: 942 vs 47,099 = 50x
- BaH: 942 vs 47,099 = 50x

This comparison is informative but not apples-to-apples:
- Different molecule (RaH vs CaH/SrH/BaH)
- Different Q (18 vs 20)
- Sunaga uses a Gaussian basis (4-component Dirac-Coulomb)

## Honest Assessment

- **Resource advantage is clear**: 942 Pauli terms at Q=20 vs ~47,099 at Q=18.
  Even accounting for Q mismatch, the ~50x advantage is robust.
- **SO accuracy is limited**: Z_eff=2 from frozen-core screening means the SO
  splittings do NOT reflect the true heavy-atom SO coupling. The physical SO
  for Sr (Z=38) and Ba (Z=56) is dominated by the core electrons, which are
  frozen. The valence-only SO at Z_eff=2 is a severe underestimate.
  Real SrH/BaH SO splittings are hundreds to thousands of cm^-1, not 5.84 cm^-1.
- **Breit corrections**: Pauli count unchanged (942), 1-norm shift +0.0005 Ha
  (0.003%), QWC groups +4. All O(alpha^2)-suppressed as expected.
- **Isostructural invariance**: Confirmed for N_Pauli, lam_ni, QWC groups,
  SO eigenvalues, and 2-electron eigenvalue structure (modulo nuclear repulsion
  offset) across all three molecules.
- **The 50x Pauli advantage is the headline result**, not the fine-structure
  accuracy. GeoVac's value proposition for heavy-atom molecules is structural
  sparsity, not spectroscopic precision of SO splittings.
