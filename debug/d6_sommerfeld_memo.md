# D₆ Sommerfeld Fine-Structure Decomposition

**Date:** 2026-04-26
**Status:** Complete

## Result

The sixth Sommerfeld Dirichlet sum at order (Zα)¹², weight 11:

```
D₆ = -(2289/512)ζ(10) + (1589/16)ζ(11)
     - (1617/64)ζ(3)ζ(8) - (1785/64)ζ(4)ζ(7) - (2065/64)ζ(5)ζ(6)
```

Verified to **62 digits** against direct Euler-sum evaluation at 400 dps.

## Product survival rule at p=6

- **ζ(2)ζ(9): ABSENT** (coefficient exactly 0)
- Surviving products: 3 = max(0, floor((2·6-5)/2)), matching prediction
- All three products ζ(3)ζ(8), ζ(4)ζ(7), ζ(5)ζ(6) are present with nonzero rational coefficients

Rule now verified through D₂–D₆ (p=2..6, weights 3..11).

## Method

### Weight-11 Euler sum decompositions

All 7 contributing Euler sums S_{r,11-r} (r=1,2,4,5,6,7,8; r=3 absent because β₃=0) were decomposed against the depth-1 MZV basis {ζ(11), ζ(2)ζ(9), ζ(3)ζ(8), ζ(4)ζ(7), ζ(5)ζ(6)}.

| Sum | ζ(11) | ζ(2)ζ(9) | ζ(3)ζ(8) | ζ(4)ζ(7) | ζ(5)ζ(6) |
|-----|-------|----------|----------|----------|----------|
| S_{1,10} | 6 | -1 | -1 | -1 | -1 |
| S_{2,9} | -27 | 9 | 2 | 6 | 4 |
| S_{4,7} | -329/2 | 84 | 0 | 21 | 4 |
| S_{5,6} | 463/2 | -126 | 0 | -21 | 0 |
| S_{6,5} | -461/2 | 126 | 0 | 21 | 1 |
| S_{7,4} | 331/2 | -84 | 0 | -20 | -4 |
| S_{8,3} | -82 | 36 | 1 | 15 | 6 |

### PSLQ strategy

Standard 6-element PSLQ fails when a coefficient is exactly 0 in the full basis (the relation has a zero entry, which PSLQ cannot find when the extra dimension inflates the search space). The fix: **sub-basis PSLQ** — try all 5-element sub-bases (dropping each element in turn). S_{5,6}, S_{6,5}, and S_{7,4} all had ζ(3)ζ(8) = 0, causing the 6-element PSLQ to fail; dropping ζ(3)ζ(8) from the basis immediately found the correct 5-element relation.

### Analytical assembly

D₆ was assembled from:
```
D₆ = (-2289/512)ζ(10) + (567/128)ζ(11) + Σ_r β_r·(S_{r,11-r} - ζ(11))
```
using exact sympy Rational arithmetic. No numerical PSLQ on D₆ itself was needed.

### Structural observation: ζ(3)ζ(8) sparsity

The ζ(3)ζ(8) column is sparse: nonzero only for S_{1,10} (Euler's formula, always full), S_{2,9}, and S_{8,3} (stuffle partner of S_{3,8}). The middle sums S_{4,7}..S_{7,4} all have zero ζ(3)ζ(8) coefficient.

## K-Sommerfeld structural separation

PSLQ at 200 dps confirms K/π ∉ ℚ-span of {D₂, D₃, D₄, D₅, D₆, 1}. K inhabits the spectral-zeta π-polynomial ring (T9 theorem); D_p lives in the full Euler-Zagier MZV algebra.

## Complete D₂–D₆ table

| p | Weight | ζ(2p-2) | ζ(2p-1) | Products |
|---|--------|---------|---------|----------|
| 2 | 3 | -5/4 | 1 | — |
| 3 | 5 | 19/8 | -11/4 | — |
| 4 | 7 | -205/64 | 71/8 | -9/2·ζ(3)ζ(4) |
| 5 | 9 | 497/128 | -467/16 | 385/32·ζ(3)ζ(6) + 75/8·ζ(4)ζ(5) |
| 6 | 11 | -2289/512 | 1589/16 | -1617/64·ζ(3)ζ(8) - 1785/64·ζ(4)ζ(7) - 2065/64·ζ(5)ζ(6) |

Sign pattern of the single-zeta terms: alternating for ζ(even), but the ζ(odd) coefficient sign at p=6 is **positive** (+1589/16), breaking the alternating pattern (p=2: +1, p=3: -11/4, p=4: +71/8, p=5: -467/16, p=6: +1589/16). Magnitudes grow rapidly.

## Data

- Assembly script: `debug/d6_analytical_assembly.py`
- K-Sommerfeld check: `debug/d6_k_sommerfeld_check.py`
- JSON: `debug/data/d6_analytical_assembly.json`
- Paper 28 §7 updated with D₆ subsection
