# Track CB: Coupled Composition Scoping — Investigation Report

**Date:** April 3, 2026
**Version:** v2.0.37
**Status:** NEGATIVE RESULT

## Summary

Investigated replacing the Phillips-Kleinman (PK) pseudopotential with explicit
cross-block two-electron integrals (ERIs) in the composed architecture.  The
hypothesis was that PK is unnecessary in second quantization because antisymmetry
is automatic.  **Result: negative.**  Cross-block ERIs add 2.56x Pauli terms and
2.30x 1-norm without improving accuracy, because the composed orbital basis
(different Z per block) creates an asymmetry: cross-block ERIs add core-valence
repulsion but no balancing cross-center nuclear attraction.

## Motivation

PK is the primary accuracy bottleneck in the composed framework (5.3% R_eq for
LiH, 11.7% BeH2, 26% H2O).  Every attempt to fix PK in the classical solver
has failed (6 entries in CLAUDE.md Section 3).  The hypothesis was that PK was
needed only for classical solvers that can't handle inter-group antisymmetry,
and could be replaced by explicit cross-block ERIs on a quantum computer.

## Phase 1: Cross-Block ERI Sparsity Census

LiH at max_n=2 (Q=30, M=15, 3 sub-blocks):

| Metric | Within-block | Cross-block |
|:-------|:-------------|:------------|
| Nonzero ERIs | 195 | 130 |
| Ratio | 1.00 | 0.667 |

The 2/3 cross/within ratio (confirmed from Track BX-3b) reflects Gaunt selection
rule filtering.  Cross-block ERIs only exist between Li_core (Z=3) and the
LiH_bond sub-blocks (Z_eff=1, Z=1).  130 ERIs is modest.

## Phase 2: Cross-Block ERI Magnitudes

From existing BX-3b data (`debug/data/cross_block_mp2_lih.json`):

| Statistic | Value |
|:----------|:------|
| Max |cross ERI| | 0.891 Ha |
| Mean |cross ERI| | 0.142 Ha |
| Median |cross ERI| | 0.074 Ha |
| Min nonzero | 0.0004 Ha |
| Sum |cross ERI| | 18.4 Ha |

Cross-block ERIs are NOT negligible.  The largest (0.891 Ha) is comparable to
within-block ERIs.  This is expected: core (Z=3) and valence (Z=1) orbitals
share the same nuclear center and have substantial radial overlap.

## Phase 3: Coupled Hamiltonian Results

### Qubit Metrics

| Metric | Composed (PK) | Coupled (no PK + cross) | Ratio | Gaussian STO-3G |
|:-------|:-------------:|:-----------------------:|:-----:|:---------------:|
| Qubits | 30 | 30 | 1.00x | 12 |
| Pauli terms | 334 | 854 | 2.56x | 907 |
| 1-norm (Ha) | 37.33 | 85.69 | 2.30x | 34.3 |
| QWC groups | 21 | 88 | 4.19x | 273 |

**Assessment:** The coupled Hamiltonian has 854 Pauli terms (2.56x composed).
This is slightly below Gaussian STO-3G (907) in Pauli count, but the 1-norm
(85.69 Ha) is 2.5x worse than both composed (37.33 Ha) and STO-3G (34.3 Ha).
QWC groups increase by 4.19x, from 21 to 88.

### Comparison to Bounds

Paper 14 Section V.E bounded cross-block Pauli impact at <=2x.  The measured
2.56x slightly exceeds this theoretical bound.  The excess comes from mixed
block-pair cross-block ERIs that were not considered in the original counting
argument.

## Phase 4: Accuracy Assessment

### 2-Electron Valence-Only FCI (orbitals 5-14)

| Configuration | Energy (Ha) | Note |
|:-------------|:-----------:|:-----|
| Composed (PK, no cross) | -7.616 | PK raises valence energy |
| No PK, no cross | -7.948 | Baseline |
| Coupled (no PK + cross) | -7.948 | Cross-block ERIs don't appear in valence subspace |

Cross-block ERIs have **zero effect** on 2-electron valence FCI because both
electrons are in the valence block — cross-block ERIs only activate when
electrons occupy both core and valence blocks simultaneously.

### 4-Electron Full-Space FCI (nuclear_rep = V_NN only)

| Configuration | Energy (Ha) | Error vs exact |
|:-------------|:-----------:|:--------------:|
| No PK, no cross | -7.195 | 10.9% |
| Composed (PK, no cross) | -6.863 | 15.0% |
| **Coupled (no PK + cross)** | **-5.734** | **29.0%** |
| Exact LiH | -8.071 | — |

The coupled approach gives the **worst** energy because cross-block ERIs add
core-valence electron *repulsion* without adding cross-center nuclear
*attraction* (h1 terms where each orbital feels all nuclei).  The asymmetry is:

- **Added:** V_ee cross-block (repulsive, +)
- **Missing:** V_ne cross-center (attractive, -)

Without cross-center h1, the coupled Hamiltonian is energetically unbalanced.

### 4-Electron FCI (nuclear_rep = V_NN + V_cross)

Including V_cross as a constant partially corrects but doesn't resolve the issue
because V_cross is the *expectation value* of the core-H nuclear interaction for
the frozen core state, not the full operator.

## Root Cause Analysis

The composed architecture defines each orbital block with its own effective
nuclear charge (Z_eff).  This means:

1. Core orbitals at Z=3 feel only the Li nucleus (via -Z^2/2n^2 in h1)
2. Valence orbitals at Z_eff=1 feel a screened Li nucleus
3. H-side orbitals at Z=1 feel only the H nucleus

No orbital feels the *other* nucleus.  This is correct for the composed
framework because Z_eff screening absorbs the mean-field cross-center nuclear
attraction.  But when you add explicit cross-block ERIs (which add repulsion
between core and valence electrons), the corresponding cross-center attraction
is missing from h1.

**This is not fixable within the current infrastructure** because cross-center
h1 matrix elements require two-center integrals:
h1_cross[p,p] = <p_Li | -Z_H/|r - R_H| | p_Li>

These are not available in the single-center hydrogenic basis.  Computing them
would require either:
- Prolate spheroidal coordinates (like Level 2)
- Multipole expansion of the cross-center potential
- Numerical integration on a 3D grid

## Structural Finding

The investigation reveals that **cross-block ERIs are the wrong coupling
mechanism for the composed orbital basis**.  The composed architecture works
*because* it replaces multi-center physics with single-center blocks + PK.
Adding cross-block ERIs partially undoes the single-center approximation
(adding V_ee cross-terms) without fully undoing it (missing V_ne cross-terms).

**PK is not redundant in second quantization for the composed basis.**  PK
provides a one-body effective potential that approximates the combined effect
of core-valence repulsion + orthogonality enforcement.  Replacing it with
explicit two-body cross-block ERIs is incomplete without also replacing the
one-body cross-center nuclear attraction.

## Recommendation

**Do not proceed with full coupled composition implementation.**  The approach
is structurally incomplete within the composed framework's single-center orbital
basis.

The correct path to improve accuracy beyond PK is either:
1. **Keep PK** and accept the structural accuracy ceiling (current approach)
2. **Use a shared coordinate system** (e.g., full N-electron hyperspherical,
   Level 4N) where cross-center physics is built into the basis — but this
   sacrifices the O(Q^2.5) sparsity advantage
3. **TC modification** (Track BX-3) which improves accuracy within the existing
   block structure without adding cross-block ERIs

## Files

| File | Description |
|:-----|:------------|
| `geovac/coupled_composition.py` | Coupled Hamiltonian builder + sector-restricted FCI |
| `tests/test_coupled_composition.py` | 13 tests (Gaunt rules, PK removal, Hermiticity, regression) |
| `debug/data/coupled_composition_lih.json` | Raw data from investigation |
| `docs/coupled_composition_scoping.md` | This report |
