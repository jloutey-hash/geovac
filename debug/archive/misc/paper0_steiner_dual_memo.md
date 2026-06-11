# Paper-0 Packing Steiner/Fermat/Voronoi Dual Probe — Memo

**Sprint:** conversational follow-up, April 17 2026.
**Conjecture tested** (PI's proposal):
The Paper-0 packing places physical nodes on integer shells $r_k = k \cdot d_0$. A complementary *Steiner / Fermat lattice* constructed from minimum-action (least-length) joins between packing points should land on half-integer shells $r = (k + \tfrac{1}{2}) \cdot d_0$. If so, its shell counts should match the Camporesi–Higuchi Dirac degeneracies $g_n = 2(n+1)(n+2)$ on $S^3$, giving a 2D geometric derivation of the Dirac spectrum's $n + 3/2$ shift.

**Verdict: FALSIFIED** in 2D. The Steiner / Fermat / Voronoi duals of the Paper-0 2D packing do **not** cluster at half-integer radii, and their counts do **not** match the Dirac degeneracy sequence.

## Methodology

Three independent dual constructions on the Paper-0 packing at $k_{\max} = 6$ (37 total points across 6 shells with $2k-1$ angular positions per shell, plus $N_{\text{init}} = 2$ on shell 1):

1. **Delaunay-triangle Fermat points.** For each Delaunay triangle of the packing, compute its Fermat point (Weiszfeld iteration if all angles < 120°; vertex otherwise). Fermat points are the natural 3-terminal Steiner minima of nearest-neighbor triples.
2. **Voronoi vertices.** The classical geometric dual of the packing — points equidistant from 3+ packing points.
3. **Clustering analysis** at half-integer vs integer target radii, with 0.1·$d_0$ tolerance window.

All computations in `debug/compute_paper0_steiner_dual.py`; raw data in `debug/data/paper0_steiner_dual.json`.

## Results

### Delaunay–Fermat dual

61 Delaunay triangles; 57 with interior Fermat points (all angles < 120°), 4 degenerate.

| Target radius | Tolerance ±0.1·$d_0$ | Count | Dirac ref $g_n$ |
|:--:|:--:|:--:|:--:|
| $1.5 d_0$ | | 1 | 4 |
| $2.5 d_0$ | | 4 | 12 |
| $3.5 d_0$ | | 3 | 24 |
| $4.5 d_0$ | | 5 | 40 |
| $5.5 d_0$ | | 1 | 60 |
| $1.0 d_0$ (int) | | 0 | — |
| $2.0 d_0$ (int) | | 2 | — |
| $3.0 d_0$ (int) | | 2 | — |
| $4.0 d_0$ (int) | | 2 | — |
| $5.0 d_0$ (int) | | 4 | — |

Half-integer total: 14. Integer total: 10. Mild preference for half-integer, but the effect is weak and the *counts are nowhere near the Dirac* $g_n = \{4, 12, 24, 40, 60\}$. Plus the full Fermat-radius distribution has significant weight at non-half-integer radii like 3.6, 4.1, 5.1 — i.e., the dual structure is broadly distributed, not sharply concentrated at either sublattice.

### Voronoi-vertex dual

58 Voronoi vertices within $r < 10 d_0$.

| Target radius | Count | Dirac ref $g_n$ |
|:--:|:--:|:--:|
| $1.5 d_0$ | 2 | 4 |
| $2.5 d_0$ | 2 | 12 |
| $3.5 d_0$ | 0 | 24 |
| $4.5 d_0$ | 0 | 40 |
| $5.5 d_0$ | 0 | 60 |
| $2.0 d_0$ (int) | 4 | — |
| $3.0 d_0$ (int) | 3 | — |
| $5.0 d_0$ (int) | 2 | — |
| $7.0 d_0$ (int) | 2 | — |

Half-integer total: 4. Integer total: 11. **Voronoi vertices actually prefer integer radii, not half-integer.** This is the opposite of what the conjecture predicted.

### Full Fermat-radius distribution (histogram, 0.1 bins, showing bins with count ≥ 2)

```
[2.80, 2.90):  2     [4.10, 4.20):  4
[3.00, 3.10):  2     [4.50, 4.60):  5
[3.40, 3.50):  2     [4.60, 4.70):  2
[3.60, 3.70):  4     [4.70, 4.80):  2
[3.90, 4.00):  2     [5.00, 5.10):  2
                     [5.10, 5.20):  8
```

The dominant clusters are at 3.6, 4.1, 4.5, 5.1, 5.2 — a mixture of positions that don't align cleanly with either integer or half-integer radii.

## Interpretation

The Paper-0 packing's Steiner / Fermat / Voronoi dual **is a real structure** (it has specific radii, specific counts, reproducible geometry), but it is **not the Dirac sublattice in disguise**.

Three reasons why the conjecture fails in this setting:

1. **Dimensional mismatch.** The Paper-0 packing is 2D (concentric circles in a plane). The Dirac spectrum's $n + 3/2$ shift lives on $S^3$ (3-sphere). The Dirac spinor bundle doesn't have a reduction to a 2D planar problem that would admit a straightforward Steiner-dual construction. Asking "what is the 2D packing's Steiner dual?" and "what is the Dirac-on-$S^3$ spectrum?" are questions about different-dimensional objects.

2. **Angular-count mismatch.** The Paper-0 packing has $2k - 1 = 1, 3, 5, 7, \ldots$ points per shell. Dirac degeneracies on $S^3$ are $4, 12, 24, 40, \ldots = 2(n+1)(n+2)$ — a quadratic-growth sequence. Dividing: $g_n / (2 \cdot \text{shell}_n) = 4/2, 12/6, 24/10, 40/14 \ldots = 2, 2, 2.4, 2.86, \ldots$ — not a constant ratio. No simple one-to-one mapping exists.

3. **Operator-order origin of the shift.** The Dirac spectrum's half-integer shift is a consequence of Dirac being a *first-order* operator on a *spinor bundle* over $S^3$. The Camporesi–Higuchi formula $|\lambda_n| = n + d/2$ (here $d = 3$) is specifically the spinor eigenvalue formula on $S^d$, which is derived from the representation theory of $\mathrm{Spin}(d+1)$ on spinor harmonics. This is not a geometric-dual statement about the scalar packing; it is a statement about the *double cover* of the rotation group acting on a different bundle.

## What the probe *does* reveal

- The packing has its own nontrivial dual structure with identifiable clusters (at $r \approx 0.73, 1.35, 2.5, 3.6, 4.1, 4.5, 5.1$ for Fermat points; different radii for Voronoi). This is geometric content worth documenting, but not the Dirac content.
- **Voronoi vertices prefer integer radii** (11 vs 4 at half-integer targets). This is structurally interesting in its own right: Voronoi vertices often land AT or NEAR packing-point positions because the Voronoi cells are asymmetric (different numbers of points on different shells).
- **Fermat points mildly prefer half-integer radii** (14 vs 10). The "mildly" is key — this is a small statistical preference, not a structural identity.

## What the probe does *not* settle

- **3D probe still open.** The full Paper-7 S³ Coulomb graph lives in 3D (quantum numbers $(n, \ell, m)$, embedded via Fock projection). A genuine 3D Steiner-dual of *that* graph would be more apt to the Dirac-on-$S^3$ question. This is a larger computation (N=37 points in 2D vs. potentially 100+ points in 3D, and Steiner-tree construction in 3D is harder than 2D), but structurally closer to the conjecture.
- **Different operational principle.** The intuition "first-order minimum-action derivation of the spinor bundle" may still be correct, but realized by a different construction than Steiner trees — e.g., a mid-sphere construction, a Weyl-asymptotic formula on the spinor bundle, or a representation-theoretic dual that's not purely geometric.
- **The packing-vs-Dirac relationship still needs a name.** The integer scalar lattice vs half-integer spinor lattice interleaving is real (it's Paper-18 taxonomy material, Paper-24 Coulomb/HO asymmetry material) — just not as a naive 2D Steiner dual of Paper 0.

## Conclusion

The specific conjecture "Paper-0 packing's Steiner dual sits at half-integer radii with Dirac degeneracy counts" is falsified by direct computation. The 2D Fermat/Voronoi dual of the Paper-0 packing does not align with the Dirac spectrum's $n + 3/2$ sublattice.

The intuition behind it — that the scalar (integer-shell) and spinor (half-integer) lattices are two complementary aspects of the GeoVac geometry, with the half-integer sublattice carrying first-order / least-action content while the integer sublattice carries second-order / packing-axiom content — **remains valid at the level of operator-order classification** (Paper 18 §IV; Paper 24 Coulomb/HO asymmetry). The two sublattices are real and structurally distinct; they just don't stand in a simple geometric-dual relationship on the 2D packing.

This is one more "no common generator" result in the Paper-18 sense: the scalar-Laplacian spectrum and the Dirac spectrum are two categorically different spectral objects on the same manifold, and no simple geometric construction (2D packing dual) unifies them. The unification would need a more sophisticated construction (representation-theoretic, or in full 3D) or may not exist at all.

## Files

- `debug/compute_paper0_steiner_dual.py` — reproducible probe driver.
- `debug/data/paper0_steiner_dual.json` — full data (packing points, Delaunay-Fermat results, Voronoi vertices, histograms, clustering).
- This memo.
