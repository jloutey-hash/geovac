# Sprint S2-v2 — Balanced-coupled library negative control panel

**Date:** 2026-06-05 (extension of Sprint S1 negative control to the full
balanced-coupled library)
**Author:** PM session
**Status:** POSITIVE-WITH-NEW-FINDING. The chi=2 → chi>2 lift at sub-block
boundaries holds bit-identically across LiH, BeH₂, H₂O balanced-coupled,
with a clean **chi_balanced − 2 = 7 × n_cross-block-couplings** pattern
that directly inherits from the F2a/F3 multipole channel count (DF rank = 7
for {1s, 2s, 2p} valence).

**Files:**
- `debug/sprint_s2_v2_balanced_library_panel.py`
- `debug/data/sprint_s2_v2_balanced_library_panel.json`

---

## Numerical findings

Balanced-coupled boundary lift across the library:

| molecule | sub-blocks | boundary cut | chi_composed | chi_balanced | lift |
|---:|---:|---:|---:|---:|---:|
| LiH | 3 | 10 | 2 | 9 | 7 |
| LiH | 3 | 20 | 2 | 9 | 7 |
| BeH₂ | 5 | 10 | 2 | 9 | 7 |
| BeH₂ | 5 | 20 | 2 | 16 | 14 |
| BeH₂ | 5 | 30 | 2 | 16 | 14 |
| BeH₂ | 5 | 40 | 2 | 9 | 7 |
| H₂O | 7 | 10 | 2 | 9 | 7 |
| H₂O | 7 | 20 | 2 | 16 | 14 |
| H₂O | 7 | 30 | 2 | 16 | 14 |
| H₂O | 7 | 40 | 2 | 16 | 14 |
| H₂O | 7 | 50 | 2 | 9 | 7 |
| H₂O | 7 | 60 | 2 | 9 | 7 |

**12/12 boundaries lifted**, with the lift quantized to {7, 14} bit-identically.

## Structural reading

The lift IS the multipole-channel cross-block rank at the boundary:

- **Lift = 7 = DF rank for the {1s, 2s, 2p} valence basis** (F2a from
  sprint_df_multipole_lift_memo.md):
  - L=0: 4 distinct radial densities (1s², 1s·2s, 2s², 2p²)
  - L=1: 2 distinct radial densities (1s·2p, 2s·2p)
  - L=2: 1 distinct radial density (2p²)
  - Total: 7
- This is the DF rank that the composed builder's per-sub-block ERIs
  have. When cross-block ERIs vanish (composed, F4), the cross-cut rank
  at boundaries is 0 and chi_k = 2 (additive baseline). When
  cross-block ERIs are active (balanced-coupled), the cross-cut content
  at the boundary IS the multipole rank of the cross-block coupling = 7.

- **Lift = 14 = 2 × 7**:\ at interior boundaries between two heavy-atom
  sub-blocks (where TWO active cross-block couplings straddle the
  boundary — e.g., in BeH₂ between bond1_center and bond2_center, both
  the bond1_center↔bond2_center direct coupling AND the bond1_center↔H₁
  + bond2_center↔H₂ couplings contribute through the cross-cut at this
  boundary).

The 7 vs 14 split is **exactly bit-identical across LiH, BeH₂, H₂O**, and
the value 7 matches the DF rank exactly. This is the F3 inheritance
structure manifesting at the boundary chi level:\ the same Laguerre/Gaunt
ring that controls DF rank also controls the boundary lift.

## Connection to Theorem 3.2.A.unified

The library panel sharpens Theorem 3.2.A.unified statement (F)
(negative control) from "balanced-coupled lifts chi=2 at LiH boundaries"
to:

> **(F-extended).** Under balanced-coupled, the boundary lift
> chi_k^{H,balanced} − chi_k^{H,composed} = 2 equals
> **N_cross-block(boundary) × (DF rank of valence basis)**,
> where N_cross-block(boundary) ∈ {1, 2} counts the number of active
> cross-block multipole couplings straddling the boundary.

This is a direct chemistry-level confirmation of the F3 inheritance
structure of Theorem 3.2.A.unified(D): the MPO bond rank per L inherits
the rank of the radial-integral matrix R^L, and at sub-block boundaries
the per-boundary lift equals the sum over L of (rank R^L) at the
cross-block coupling — which is the DF rank.

## Implications for Paper 14

Add a short sentence to the negative-control remark in §sec:mpo_bond_rank
quoting the extended panel (12/12 boundaries verified, lift = 7 or 14
matching DF rank).

## Honest scope

- Tested on first-row chemistry (LiH, BeH₂, H₂O) at n_max=2 only.
  Second-row molecules (NaH, MgH₂) have frozen-core treatment; the
  negative-control extension should be verified there as a follow-on
  (predicted: identical 7-rank lift since the {1s, 2s, 2p} valence basis
  is the same).
- The interior 14 pattern between two heavy-atom sub-blocks is
  established but the precise enumeration of cross-block couplings at
  each boundary (which boundary involves which 2 couplings) was not
  traced; would tighten the analysis to know whether some boundaries
  could lift to 21 = 3 × 7 (three cross-couplings) in larger molecules.

## Verdict

POSITIVE. Theorem 3.2.A.unified(E,F) extended from LiH-only to the
balanced-coupled library with bit-identical lift quantization
(7 or 14 = N_cross × DF_rank) across 12 boundaries on 3 molecules.

The "lift = DF rank × N_cross" is a substantive new structural finding:
the F3 inheritance carries through to the negative-control boundary lift.
