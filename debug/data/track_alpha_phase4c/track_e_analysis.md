# Track alpha-E: S^1 Fiber Fock Weight and zeta(2)

## Summary

**Overall verdict: NEGATIVE.** The hypothesis that F = zeta(2) = pi^2/6
of Paper 2's formula K = pi*(B + F - Delta) comes from the S^1 m-fiber of
the Hopf bundle is not supported. Four independent tests were performed
and none recovered F (or any clean rational multiple of pi^2) from the
fiber data.

The only "near-miss" is that the continuum-scaled path-graph spectral
zeta `(pi/n)^2 * (n^2 - 1)/6 -> pi^2/6` in the large-n limit, but the
finite-size correction at n = 2l_max + 1 = 5 is `-pi^2/150`, which is
not a clean rational multiple of Delta = 1/40.

**Recommendation:** F and Delta remain INDEPENDENT transcendental
injections in Paper 2's formula. They do not arise from the Fock
projection's fiber structure at the Paper 2 shell cutoff. If F has a
geometric origin, it is not in the discrete Hopf fiber at n_max = 3.

## Subtask 1: m-fiber Fock weight matrices

**Verdict: m-fiber is strictly diagonal.** For each (n,l) sector with
n <= 3, the matrix

    W^{(n,l)}_{m,m'} = <n,l,m | (p^2 + p0^2)^{-2} | n,l,m'>

is a scalar multiple of the identity, as expected: the weight w(p)
depends only on |p|^2, and the basis Y_{lm}(Omega_p) is orthonormal in
angle, so angular integration collapses to delta_{m,m'}. The diagonal
values at p_0 = 1 reproduce alpha-C exactly:

| (n,l) | diagonal value | dim |
|:-----:|:--------------:|:---:|
| (1,0) | 7/16  | 1 |
| (2,0) | 5/8   | 1 |
| (2,1) | 3/8   | 3 |
| (3,0) | 5/8   | 1 |
| (3,1) | 17/32 | 3 |
| (3,2) | 11/32 | 5 |

**Exploratory natural variant:** The operator L_x = (L_+ + L_-)/2 in the
l = 1 basis is tridiagonal with entries sqrt(2)/2, confirming that a
rotation-generator-type perturbation DOES couple m by +/- 1 and has
path-graph P_{2l+1} adjacency structure — consistent with alpha-D's
finding that the fiber is a path, not a cycle. But the rotation-
invariant Fock weight itself cannot access this structure; any
non-trivial m-coupling requires breaking the SO(3) symmetry of w(p).

**Conclusion:** The m-fiber is information-free under the standard
Fock weight. Any derivation of F via the fiber must use a different
operator (spectral zeta of the discrete fiber Laplacian, below) or a
symmetry-broken weight.

## Subtask 2: S^1 fiber spectral zeta (four weighting schemes)

Using the path-graph spectrum lambda_j(P_{2l+1}) = 2 - 2 cos(pi j/(2l+1))
for j = 0, 1, ..., 2l, the fiber spectral zeta at s = 1 is

    Sum_{j=1}^{2l} 1 / lambda_j = (n_fiber^2 - 1)/6

where n_fiber = 2l+1. (Closed form verified symbolically; the cycle
C_n form is (n^2-1)/12, which is NOT what the fibers are — alpha-D
confirmed they are paths.)

Weighted sums over n <= 3:

| Weight scheme                 | sum (exact) | sum (numeric) | rational? |
|:------------------------------|:------------|:-------------:|:---------:|
| uniform                       | 20/3        | 6.6666...     | yes |
| degeneracy (2l+1)             | 28          | 28            | yes |
| Casimir l(l+1)                | 88/3        | 29.333...     | yes |
| Hopf (2l+1) l(l+1) (alpha-C)  | 136         | 136           | yes |

**All four weighted sums are pure rationals. None contain pi.** They
cannot equal any rational multiple of pi^2/6 because they themselves
contain no transcendentals.

**Closest "near-misses":**
- uniform = 20/3 ~ 6.667 vs 6 F = pi^2 ~ 9.870 (32% gap)
- All other ratios are > 50% off

**Conclusion: raw path-graph fiber zeta DOES NOT produce F.** The
discrete fiber structure at n_max = 3 is in a purely rational algebraic
world, disjoint from pi^2/6. This is the same algebraic separation that
alpha-D found between graph Laplacian invariants (small integers,
golden ratios) and the Paper 2 transcendental targets.

## Subtask 3: Continuum limit and finite-size correction

The raw path-graph spectral zeta grows as n^2/6 - 1/6. To get a finite
continuum limit we rescale with (pi/n)^2, which corresponds to mapping
the discrete eigenvalues onto the continuum S^1 (length pi) Dirichlet
spectrum k_j^2 = j^2.

Data table:

| n = 2l+1 | s(n) = (n^2-1)/6 | s(n)/n^2 | (pi/n)^2 s(n) | diff from pi^2/6 |
|:--------:|:----------------:|:--------:|:-------------:|:----------------:|
| 1  | 0        | 0         | 0         | -pi^2/6      |
| 3  | 4/3      | 0.1481    | 1.4622    | -pi^2/18     |
| 5  | 4        | 0.1600    | 1.5791    | -pi^2/150    |
| 7  | 8        | 0.1633    | 1.6114    | ...          |
| 9  | 40/3     | 0.1646    | 1.6246    |              |
| 11 | 20       | 0.1653    | 1.6313    |              |
| 21 | 220/3    | 0.1663    | 1.6412    |              |
| inf | -       | 1/6       | pi^2/6    | 0            |

**Verdict (Subtask 3): the continuum limit IS pi^2/6**, as expected
from the well-known continuum S^1 Dirichlet spectral zeta. This is a
mathematical identity, not a new derivation of F.

**Finite-size correction at n = 5 (Paper 2's l_max = 2 cutoff):**

    (pi/5)^2 * (5^2 - 1)/6 - pi^2/6 = 4 pi^2/25 - pi^2/6
                                    = (24 - 25) pi^2 / 150
                                    = -pi^2/150
                                    ~ -0.0658

Ratio to Delta = 1/40:

    -pi^2/150 / (1/40) = -4 pi^2/15
                       ~ -2.632

This is NOT a clean rational multiple of Delta. It is -4 pi^2/15, which
is a new transcendental ratio, not 1/40 nor any simple modification.

**The finite-size correction does NOT identify with Delta.**

## Subtask 4: Cross-checks (n_max invariance)

Since the path-graph spectra depend only on the node count 2l+1, the
p_0 and Z invariances are automatic: neither enters the fiber zeta. The
only non-trivial invariance is n_max.

Comparing continuum-scaled sums at n_max = 3 vs n_max = 4:

| Scheme               | n_max=3 value | n_max=4 value (exact) | n_max=4 numeric |
|:---------------------|:-------------:|:----------------------|:---------------:|
| uniform              | 4.503         | 10228 pi^2 / 11025    | 9.156           |
| degeneracy           | 16.669        | 428 pi^2 / 105        | 40.230          |
| Casimir l(l+1)       | 15.323        | 52568 pi^2 / 11025    | 47.059          |
| Hopf (2l+1)l(l+1)    | 64.920        | 2728 pi^2 / 105       | 256.422         |

**Verdict: NOT n_max-stable.** The continuum-scaled sums grow with
n_max (adding the (4,0..3) cells adds significant contributions under
every weighting). None of the schemes produce an n_max-invariant
quantity equal to F. A genuine fiber-F would have to be stable under
extending the basis, and it is not.

## Overall verdict: NEGATIVE

The "F from S^1 fiber" hypothesis fails on all four tests:

1. **The Fock weight w(p) is strictly m-diagonal**, so the fiber
   direction carries no information from the momentum-space hydrogenic
   basis directly. Any fiber contribution must come from a separate
   spectral object (the path-graph Laplacian) layered on top.

2. **All four weighted sums of the raw fiber zeta are pure rationals.**
   They cannot equal F or any rational multiple of pi^2. Closest ratio
   gap is 32% (uniform vs pi^2).

3. **The continuum-scaled limit is pi^2/6** (as a mathematical
   identity), but the finite-size correction at n = 5 is -pi^2/150,
   which is -4 pi^2/15 in units of Delta = 1/40 — not a rational
   multiple of Delta.

4. **The fiber sums are not n_max-invariant.** A genuine F would have
   to be independent of the UV cutoff; the fiber sums grow with it.

**Interpretation.** F = zeta(2) in Paper 2 is a continuum zeta value,
living in the "embedding" tier of Paper 18's exchange constant
taxonomy (things that enter only when a continuum is projected onto a
discrete graph, or vice versa). At the Paper 2 cutoff n_max = 3, the
discrete Hopf fiber is too small (max dim = 5 for the l=2 cell) to
resolve the pi^2/6 value to the precision Paper 2 requires.

More broadly: the Fock projection is first-order-rigid (per Paper 24's
HO rigidity theorem and Paper 23's Fock rigidity theorem). Its discrete
graph structure is rational (alpha-C's 7/16, 5/8, ...), and its
fiber structure is also rational (path zeta = (n^2-1)/6). The pi in
Paper 2's K formula must enter at a different layer — most likely the
step that embeds the discrete graph into the continuum S^3 manifold
where Weyl-Selberg spectral identities produce pi's naturally.

**Recommendation: F is not in the discrete fiber. Paper 2's F should
be reframed as an embedding exchange constant (Paper 18 taxonomy),
not as a graph invariant.** The alpha-C identity for B and the
alpha-E negative for F together suggest that Paper 2's K formula is
a sum of ingredients from different algebraic tiers:
  - B = 42: rational, from the SO(3) Casimir trace on the base (alpha-C).
  - F = pi^2/6: transcendental, from continuum S^1 zeta regularization
    (alpha-E negative, not from fiber).
  - Delta = 1/40: rational, possibly from a finite-size correction,
    but -pi^2/150 != any rational multiple of 1/40 (alpha-E negative).
  - pi prefactor: from the Fock conformal embedding of S^3 into R^4.

Paper 2's combination K = pi * (B + F - Delta) thus mixes three
different exchange constant tiers and cannot be derived as a single
spectral trace from either the graph or its fiber.

## Files

- `debug/track_alpha_e.py` — full computation
- `debug/data/track_alpha_phase4c/track_e_fiber_zeta.json` — numerical data
- This analysis
