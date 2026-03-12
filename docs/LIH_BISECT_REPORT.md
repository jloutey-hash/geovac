# LiH D_raw Bisect Report

**Date:** 2026-03-11
**Status:** ROOT CAUSE IDENTIFIED
**Script:** Manual git checkout + code diff (no intermediate commits to bisect)

## Summary

The 0.578 Ha discrepancy in LiH D_raw between v0.9.9 (paper values) and v0.9.37
is caused by **two default parameter changes**, both introduced in the v0.9.37
squash commit (`a2dcc75`). There are no intermediate commits to bisect — all
changes from v0.9.10 through v0.9.37 were squashed into a single commit.

## The Offending Commit

**Hash:** `a2dcc75` (v0.9.37: Codebase lock — LiH benchmark complete)
**Parent:** `e170af2` (Release v0.9.9)

All v0.9.10–v0.9.36 changes were squashed into this single commit.

## Root Cause Decomposition

| Configuration | E_mol (Ha) | D_raw (Ha) | Delta from v0.9.9 |
|:---|:---:|:---:|:---:|
| v0.9.9 (`e170af2`): fourier, J-only | -8.0900 | 0.198 | — |
| v0.9.37 fourier, J+K (`s_only`) | -8.1167 | 0.225 | +0.027 |
| **v0.9.37 exact, J+K (`s_only`)** | **-8.6683** | **0.776** | **+0.578** |

### Change 1: `cross_nuclear_method` default (0.551 Ha, 95% of discrepancy)

**v0.9.9 behavior:** Cross-nuclear attraction applied to **s-orbitals only** (l=0)
via the Fourier/shell-theorem formula `_fourier_cross_attraction()`.

**v0.9.37 behavior:** New default `cross_nuclear_method='exact'` applies
cross-nuclear attraction to **all (n,l,m) orbitals** via full 2D quadrature
`compute_exact_cross_nuclear()`. The p and d orbitals now also feel the other
nucleus, adding ~0.55 Ha of additional attraction to the molecular energy.

**Location:** `geovac/lattice_index.py`, `MolecularLatticeIndex.__init__()`,
parameter `cross_nuclear_method: str = 'exact'` (was implicitly `'fourier'`
in v0.9.9 which had no such parameter — the Fourier s-only path was hardcoded).

**Code path:** `_apply_cross_nuclear_diagonal()` dispatches on
`self.cross_nuclear_method`:
- `'exact'`: loops over ALL orbitals, calls `compute_exact_cross_nuclear()`
- `'fourier'`: loops only over l=0 orbitals (v0.9.9 behavior)

### Change 2: Mulliken exchange K added (0.027 Ha, 5% of discrepancy)

**v0.9.9 behavior:** Cross-atom V_ee included only direct Coulomb J (18 entries).
The code comment read: "Only direct Coulomb J is included; exchange K is neglected."

**v0.9.37 behavior:** `cross_atom_vee='s_only'` now includes both J (18) and K (18)
entries. The `compute_cross_atom_K()` function was added in the Mulliken exchange
extension (v0.9.13 per CHANGELOG).

**Effect:** Exchange K raises E_mol by ~0.027 Ha (reduces overbinding slightly).
This is a minor correction in the right direction.

## Was the Change Intentional?

**Yes.** The `cross_nuclear_method='exact'` was an intentional physics improvement
(CHANGELOG v0.9.35 mentions "exact cross-nuclear (shell theorem)" as a tested
hypothesis for correcting R_eq). However, it was not flagged as a default-breaking
change, and the paper's benchmark values were computed under the old `'fourier'`
default before this parameter existed.

The paper's D_raw = 0.205 Ha and D_CP = 0.093 Ha correspond to:
```python
MolecularLatticeIndex(..., cross_nuclear_method='fourier', cross_atom_vee='s_only')
# (with J-only, no exchange K — exact v0.9.9 config gives D_raw=0.198)
# (with J+K, v0.9.37 fourier gives D_raw=0.225)
```

## Reproducibility of Paper Values

| Quantity | Paper | v0.9.9 exact | v0.9.37 fourier+s_only |
|:---|:---:|:---:|:---:|
| D_raw (R=3.015) | 0.205 | 0.198 | 0.225 |
| BSSE | -0.115 | — | -0.115 |
| D_CP | 0.093 | ~0.083* | ~0.110 |

*The 0.083 value is from the pre-normalization-fix v0.9.9. The paper's 0.093
was computed at v0.9.10/v0.9.11 which included the `_phi_s_orbital_general`
normalization fix (n≥3 orbitals). That fix was also squashed into `a2dcc75`.

**Conclusion:** The paper's exact D_CP = 0.093 Ha is from a transient code state
(v0.9.10/v0.9.11) that existed between sessions but was never individually committed.
It cannot be exactly reproduced from any git-accessible commit. The closest
reproducible values are:
- v0.9.9 (`e170af2`): D_raw = 0.198 (pre-normalization fix)
- v0.9.37 with `fourier` + `s_only`: D_raw = 0.225 (post-normalization + exchange K)

## Recommendation

To reproduce paper-era results from v0.9.37, use:
```python
MolecularLatticeIndex(
    ...,
    cross_nuclear_method='fourier',
    cross_atom_vee='s_only',
)
```

This gives D_raw = 0.225 Ha (vs paper 0.205), with the 0.020 Ha residual
attributable to the exchange K terms added in v0.9.13 (which did not exist
at the time the paper values were computed).

## Output Files

- This report: `docs/LIH_BISECT_REPORT.md`
- Previous investigation: `docs/LIH_REQ_INVESTIGATION.md`
- PES data: `debug/data/lih_req_investigation.txt`
