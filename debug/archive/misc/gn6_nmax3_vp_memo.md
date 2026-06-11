# GN-6: Vacuum Polarization at n_max=3

## Summary

Extended the graph-native QED vacuum polarization from n_max=2 (3x3, all rational) to n_max=3 (13x13). The computation runs in exact sympy arithmetic in ~0.3s total.

## Graph at n_max=3

- 28 Dirac states, 14 Fock nodes, 13 edges (5 T-type, 8 L-type)
- Betti numbers: beta_0=3 (three l-sectors), beta_1=2
- 260 nonzero vertex tensor entries out of 10,192 possible (97.4% sparse)

## Pi matrix structure

Pi decomposes into three independent l-sector blocks with zero cross-block coupling:

| Block | Edges | Size | Rational? |
|-------|-------|------|-----------|
| l=0   | 0,1   | 2x2  | Yes (diagonal) |
| l=1   | 2-8   | 7x7  | Yes |
| l=2   | 9-12  | 4x4  | **No** -- 4 entries with sqrt(6) |

The l=0 block is diagonal (the two T-edges in different shells do not couple). The l=1 block has 4 nonzero off-diagonal pairs, all rational. The l=2 block has off-diagonal entries Pi[9,10] = Pi[12,11] = 64*sqrt(6)/1225 and Pi[10,11] = 192/1225.

## Pi-free certificate

Pi is **pi-free** (no transcendentals) but **NOT all-rational**. The sqrt(6) entries in the l=2 block arise because the l=2 CG coefficients involve sqrt(6) from <2,m_l,1/2,m_s|5/2,m_j>, and the bilinear VP trace does not cancel these when the photon edge connects states with different m_l values. This falsifies the prediction that the bilinear trace always contracts sqrt(rational) CG content back to Q.

## Eigenvalues

13 eigenvalues: 7 rational, 6 irrational (containing sqrt(2), sqrt(7), sqrt(10)). Two eigenvalues are **negative** (-0.186 from l=1, -0.060 from l=2), so Pi is NOT positive semidefinite at n_max=3 (it was PSD at n_max=2). All three n_max=2 eigenvalues (32/225, 32/45, 32/15) embed exactly in the n_max=3 spectrum. All diagonal entries have numerator factor 32.

## Scaling

- Trace: 3872/735 (vs 224/75 at n_max=2), ratio 605/343
- Per-edge trace **decreases**: 224/225 (n_max=2) to 3872/9555 (n_max=3), ratio 0.407
- Interpretation: VP is distributed over more edges but concentrates less per edge

## Data

- JSON: `debug/data/gn6_nmax3_vp.json`
- Script: `debug/gn6_nmax3_vp.py`
- Tests: `tests/test_graph_qed_nmax3.py` (37 tests)
