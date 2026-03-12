# LiH Hamiltonian Diagnostic Report

**Date:** 2026-03-11
**System:** LiH, nmax=3, R=3.015 bohr, 4 electrons
**Script:** `debug/diagnose_h1_matrix.py`

## 1. Energy Summary Across Configurations

E_sep = E(Li) + E(H) = -7.392086 + -0.500000 = -7.892086 Ha
V_NN  = Z_A * Z_B / R = 3*1/3.015 = 0.995025 Ha

| Configuration | E_elec | V_NN | E_mol | D_raw | D_CP est |
|:---|:---:|:---:|:---:|:---:|:---:|
| fourier+False | -10.3300 | 0.9950 | -9.3350 | 1.4429 | ~1.328 |
| fourier+s_only | -9.1118 | 0.9950 | -8.1167 | 0.2247 | ~0.110 |
| exact+s_only | -9.6633 | 0.9950 | -8.6683 | 0.7762 | ~0.661 |
| exact+True | -9.1125 | 0.9950 | -8.1175 | 0.2254 | ~0.110 |

**Paper target:** D_raw ~ 0.205, D_CP ~ 0.093, E_mol ~ -8.097

## 2. Cross-Nuclear Attraction: Element-by-Element Comparison

### Li orbitals feeling H nucleus (Z_other=1)

| Orbital | (n,l,m) | Pure -Z^2/2n^2 | V_cross (exact) | V_cross (fourier) | Diff |
|:---|:---:|:---:|:---:|:---:|:---:|
| 1s | (1,0,0) | -4.500000 | -0.331675 | -0.331675 | 0.000000 |
| 2s | (2,0,0) | -1.125000 | -0.326884 | -0.326950 | 0.000066 |
| 2p(-1) | (2,1,-1) | -1.125000 | -0.307546 | 0.000000 | -0.307546 |
| 2p(+0) | (2,1,0) | -1.125000 | -0.372892 | 0.000000 | -0.372892 |
| 2p(+1) | (2,1,1) | -1.125000 | -0.307546 | 0.000000 | -0.307546 |
| 3s | (3,0,0) | -0.500000 | -0.229026 | -0.229001 | -0.000025 |
| 3p(-1) | (3,1,-1) | -0.500000 | -0.219559 | 0.000000 | -0.219559 |
| 3p(+0) | (3,1,0) | -0.500000 | -0.292381 | 0.000000 | -0.292381 |
| 3p(+1) | (3,1,1) | -0.500000 | -0.219559 | 0.000000 | -0.219559 |
| 3d(-2) | (3,2,-2) | -0.500000 | -0.236264 | 0.000000 | -0.236264 |
| 3d(-1) | (3,2,-1) | -0.500000 | -0.280978 | 0.000000 | -0.280978 |
| 3d(+0) | (3,2,0) | -0.500000 | -0.361361 | 0.000000 | -0.361361 |
| 3d(+1) | (3,2,1) | -0.500000 | -0.280978 | 0.000000 | -0.280978 |
| 3d(+2) | (3,2,2) | -0.500000 | -0.236264 | 0.000000 | -0.236264 |
| **Total** | --| --| **-4.002913** | **-0.887626** | **-3.115287** |

### H orbitals feeling Li nucleus (Z_other=3)

| Orbital | (n,l,m) | Pure -Z^2/2n^2 | V_cross (exact) | V_cross (fourier) | Diff |
|:---|:---:|:---:|:---:|:---:|:---:|
| 1s | (1,0,0) | -0.500000 | -0.985415 | -0.985415 | -0.000000 |
| 2s | (2,0,0) | -0.125000 | -0.557570 | -0.557776 | 0.000206 |
| 2p(-1) | (2,1,-1) | -0.125000 | -0.601975 | 0.000000 | -0.601975 |
| 2p(+0) | (2,1,0) | -0.125000 | -0.802397 | 0.000000 | -0.802397 |
| 2p(+1) | (2,1,1) | -0.125000 | -0.601975 | 0.000000 | -0.601975 |
| 3s | (3,0,0) | -0.055556 | -0.278002 | -0.278099 | 0.000097 |
| 3p(-1) | (3,1,-1) | -0.055556 | -0.292896 | 0.000000 | -0.292896 |
| 3p(+0) | (3,1,0) | -0.055556 | -0.337937 | 0.000000 | -0.337937 |
| 3p(+1) | (3,1,1) | -0.055556 | -0.292896 | 0.000000 | -0.292896 |
| 3d(-2) | (3,2,-2) | -0.055556 | -0.315896 | 0.000000 | -0.315896 |
| 3d(-1) | (3,2,-1) | -0.055556 | -0.337042 | 0.000000 | -0.337042 |
| 3d(+0) | (3,2,0) | -0.055556 | -0.355553 | 0.000000 | -0.355553 |
| 3d(+1) | (3,2,1) | -0.055556 | -0.337042 | 0.000000 | -0.337042 |
| 3d(+2) | (3,2,2) | -0.055556 | -0.315896 | 0.000000 | -0.315896 |
| **Total** | --| --| **-6.412494** | **-1.821290** | **-4.591204** |

### Physical reasonableness check

- Point-charge limit (Li orbs <- H nuc): -Z_B/R = -0.3317 Ha
- Point-charge limit (H orbs <- Li nuc): -Z_A/R = -0.9950 Ha
- Li 1s cross-nuclear (exact): -0.3317 Ha (OK: |V| < Z_B/R)
- H 1s cross-nuclear (exact): -0.9854 Ha (OK: |V| < Z_A/R)

- Orbitals with cross-nuclear (exact): 28 / 28
- Orbitals with cross-nuclear (fourier): 6 / 28

## 3. Off-Diagonal H1 (Bridge Coupling)

Bridge matrix elements (inter-atom H1 block) should be identical
between exact and fourier since only the diagonal differs.

Max |H1_AB(exact) - H1_AB(fourier)| = 0.00e+00

### Top 5 bridge matrix elements (|H1_AB|)

| Li orbital | H orbital | Value (Ha) |
|:---|:---|:---:|
| 3d(2) | 3d(2) | 0.016658 |
| 3d(-2) | 3d(-2) | 0.016658 |
| 3d(0) | 3d(0) | 0.016658 |
| 3d(1) | 3d(1) | 0.016658 |
| 3d(-1) | 3d(-1) | 0.016658 |

## 4. Cross-Atom V_ee Comparison

| Config | ERI entries | Cross-atom | E_mol (Ha) | Delta from fourier+False |
|:---|:---:|:---:|:---:|:---:|
| fourier+False | --| --| -9.335011 | +0.000000 |
| fourier+s_only | --| --| -8.116739 | +1.218272 |
| exact+s_only | --| --| -8.668316 | +0.666696 |
| exact+True | --| --| -8.117521 | +1.217490 |

## 5. Discrepancy Decomposition

Starting from fourier+no_cross_vee: E_mol = -9.335011

| Change | dE (Ha) | % of total |
|:---|:---:|:---:|
| Add cross-atom J+K (fourier->fourier+s_only) | +1.2183 | 182.7% |
| Change cross-nuclear fourier->exact | -0.5516 | -82.7% |
| **Total** | **+0.6667** | **100%** |

Extra cross-nuclear on p/d orbitals (exact only):
  Li: 11 orbitals, total V_cross = -3.115329 Ha
  H:  11 orbitals, total V_cross = -4.591506 Ha
  Combined extra diagonal shift = -7.706835 Ha

## 6. Configuration Matching Paper Values

| Quantity | Paper | fourier+s_only | exact+s_only | exact+True |
|:---|:---:|:---:|:---:|:---:|
| D_raw (Ha) | 0.205 | 0.2247 | 0.7762 | 0.2254 |
| D_CP est (Ha) | 0.093 | 0.1097 | 0.6612 | 0.1104 |
| E_mol (Ha) | -8.097 | -8.1167 | -8.6683 | -8.1175 |

**Closest match:** `fourier+s_only` (D_raw=0.225) is closest to paper (0.205).
The 0.020 Ha residual comes from Mulliken exchange K terms added post-v0.9.9.

## 7. Root Cause and Recommendation

### Root Cause

The `cross_nuclear_method` default changed from implicit `'fourier'` (v0.9.9)
to explicit `'exact'` (v0.9.35+). The exact method applies cross-nuclear
attraction to ALL (n,l,m) orbitals via 2D quadrature, while fourier only
treats s-orbitals (l=0). This adds ~0.55 Ha of overbinding.

The p/d cross-nuclear integrals are individually reasonable (each < -Z/R),
but their aggregate effect in FCI is too strong because the CI wavefunction
can exploit ALL lowered diagonal elements variationally. Without compensating
l>0 cross-atom V_ee screening, this violates the variational balance that
the fourier/s-only design maintained.

### Recommendation

To reproduce paper results, use:
```python
MolecularLatticeIndex(
    ...,
    cross_nuclear_method='fourier',
    cross_atom_vee='s_only',
)
```

**Critical finding:** `exact+True` gives D_raw=0.225, nearly identical to
`fourier+s_only` (0.225). The two BALANCED configurations converge to the
same answer. Only the IMBALANCED `exact+s_only` overbinds (D_raw=0.776).

This means the exact cross-nuclear integrals are correct --the issue is
purely a mismatch between cross-nuclear scope (all l) and V_ee scope (s-only).

## 8. Energy Decomposition (decompose_energy)

| Component | fourier+s_only | exact+s_only | exact+True | Delta (exact_s - fourier_s) |
|:---|:---:|:---:|:---:|:---:|
| T | -0.1330 | -0.2030 | -0.1329 | -0.0699 |
| V_nA | -8.8523 | -8.8662 | -8.8526 | -0.0139 |
| V_nB | -0.8458 | -0.2446 | -0.8439 | +0.6012 |
| V_cross_A | -0.6632 | -0.6646 | -0.6633 | -0.0014 |
| V_cross_B | -1.7919 | -1.5102 | -1.7920 | +0.2817 |
| V_bridge | -0.0006 | -0.0027 | -0.0006 | -0.0021 |
| V_ee | 3.1751 | 1.8280 | 3.1727 | -1.3471 |
| V_NN | 0.9950 | 0.9950 | 0.9950 | +0.0000 |
| E_total | -8.1167 | -8.6683 | -8.1175 | -0.5516 |

Key: V_cross_A = Li electrons feeling H nucleus, V_cross_B = H electrons feeling Li nucleus
