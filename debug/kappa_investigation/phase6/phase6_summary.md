# Phase 6: Quaternionic Hopf Fibration — Does the Paper 2 Recipe Generalize?

**Date:** 2026-03-22

## Research Question

Paper 2 derives α from the complex Hopf fibration S¹ → S³ → S². The four Hopf fibrations correspond to the four division algebras R, C, H, O and (arguably) to the Standard Model gauge groups {1}, U(1), SU(2), SU(3). Does the Paper 2 recipe, applied to the quaternionic Hopf S³ → S⁷ → S⁴, produce a recognizable weak coupling constant?

## Critical Distinction: Two Definitions of B

Paper 2's B is a **branching-based** quantity: it sums the SO(3) Casimir `l(l+1)` weighted by degeneracy `(2l+1)` over the SO(4) → SO(3) decomposition of each S³ shell. Call this **B_branch**.

A simpler "formal" quantity is **B_formal = (1/2) Σ g_k |λ_k|**, using only eigenvalues and degeneracies of S^d.

**These coincide ONLY for d = 3** (because the SO(4) embedding index = 1):

| d | m=2 | B_branch | B_formal | Match? |
|:--|:----|:---------|:---------|:-------|
| 2 | m=2 | 12 | 18 | NO |
| **3** | **m=2** | **42** | **42** | **YES** |
| 4 | m=2 | 96 | 80 | NO |
| 7 | m=2 | 462 | 308 | NO |

Paper 2's selection principle `B/N = dim(S³) = 3` uses B_branch. For d ≠ 3, B_branch/N ≠ d at any integer m (Paper 2 Sec VIII.E is **correct**).

## New Algebraic Identity: B_formal(2)/N(2) = d universally

The FORMAL quantity B_formal = (1/2) Σ g_k |λ_k| satisfies B/N = d at m = 2 for **all** d. This is a purely algebraic identity:

```
B_formal(2) = (1/2)(d+1)d(d+4),   N(2) = (d+1)(d+4)/2

B_formal(2)/N(2) = d   for all d

Setting B_formal(m)/N(m) = d gives:  m² + dm - 2(d+2) = 0
Unique positive root: m = 2 always.
```

| Sphere | B_formal(2) | N(2) | B/N = dim |
|:-------|-----:|-----:|----------:|
| S¹ | 5 | 5 | **1** |
| S² | 18 | 9 | **2** |
| S³ | 42 | 14 | **3** |
| S⁴ | 80 | 20 | **4** |
| S⁷ | 308 | 44 | **7** |
| S⁸ | 432 | 54 | **8** |
| S¹⁵ | 2280 | 152 | **15** |

**Implication:** Paper 2's n_max = 3 selection is NOT specific to S³ or to Hopf fibrations — it is a universal property of eigenvalue moments on round spheres. The formula's specificity lies elsewhere (in the combination rule, the fiber zeta, and the cubic structure).

## Complete Recipe Applied to All Four Hopf Fibrations

| | Real | Complex | Quaternionic | Octonionic |
|:--|:--:|:--:|:--:|:--:|
| **Bundle** | S⁰→S¹→S¹ | S¹→S³→S² | S³→S⁷→S⁴ | S⁷→S¹⁵→S⁸ |
| **Gauge group** | {1} | U(1) | SU(2) | G₂→SU(3) |
| **dim(fiber/total/base)** | 0/1/1 | 1/3/2 | 3/7/4 | 7/15/8 |
| **B** | 5 | 42 | **308** | 2280 |
| **F** | 0 (degen.) | π²/6 ≈ 1.645 | **3/4** | 49/120 ≈ 0.408 |
| **Δ** | 1/12 | 1/40 | **1/144** | 1/544 |
| **K = π(B+F-Δ)** | 15.4 | 137.036 | **970.0** | 7164.1 |
| **1/α (cubic root)** | 15.4 | **137.036** | **970.0** | 7164.1 |

### Coupling Comparison

| Hopf fibration | 1/α | Nearest known coupling | Match? |
|:---------------|:---:|:----------------------|:------:|
| Complex | 137.036 | α_em = 1/137.036 | **YES (8.8×10⁻⁸)** |
| Quaternionic | 970.0 | α_W = 1/29.5 | No (33× off) |
| Octonionic | 7164.1 | α_s = 1/8.5 | No (843× off) |

**The quaternionic and octonionic values do not match any known coupling constant**, regardless of prefactor choice.

## Ratio Structure

| Ratio | Value | Nearest simple number |
|:------|:------|:---------------------|
| K_quat / K_em | 7.078 | 7 (= dim S⁷) |
| K_oct / K_em | 52.28 | — |
| B_quat / B_em | 308/42 = 22/3 | 7.33 |
| N_quat / N_em | 44/14 = 22/7 | π (!) |

K_quat/K_em ≈ 7 is suggestive but not exact (deviation 1.1%). The 22/7 = π appearance in N ratios is a coincidence.

## Representation-Theoretic Obstruction

The Paper 2 formula works because of three properties unique to S³:

1. **SO(4) is a product group:** SO(4) ≅ SU(2)×SU(2), giving the clean Clebsch-Gordan branching l = 0, ..., n-1 within each shell. **No analog for SO(8).**

2. **S³ is a Lie group:** S³ ≅ SU(2) enables the Peter-Weyl decomposition. S⁷ is NOT a Lie group (it is a parallelizable manifold but only carries a Moufang loop structure from the octonions, not a group structure).

3. **Embedding index = 1:** The diagonal SO(3) ⊂ SO(4) has index 1, which yields the identity ⟨C_{SO(3)}⟩_n = C_{SO(4)}/2 = |λ_n|/2. For the quaternionic Hopf, the SO(5) ⊂ SO(8) embedding through Sp(2)×Sp(1) has index ratio 5/14, which would give B/N = 5 (not 7) if used as the branching factor.

The factor 1/2 in B = (1/2)Σ g_k|λ_k| is a formal algebraic identity that works for all spheres, but its interpretation as the trace of the base Casimir (the SO(3) angular momentum content) is specific to S³.

## Fiber Spectral Zeta Values

The fiber spectral zeta F = Σ 1/|λ_k| (without degeneracy) has exact closed forms:

| Fiber | F | Exact form |
|:------|:--|:-----------|
| S⁰ | degenerate | (discrete space) |
| S¹ | π²/6 = 1.6449... | ζ(2) (transcendental) |
| S³ | 3/4 | H₂/2 (rational!) |
| S⁷ | 49/120 | H₆/6 (rational!) |

General formula: F(S^d) = H_{d-1}/(d-1) where H_n is the n-th harmonic number.

Remarkably, the complex Hopf fiber (S¹) is the ONLY one with a transcendental F. All higher fibers give rational values. This is because ζ(2) = π²/6 is transcendental, while the telescoping sums for higher-dimensional spheres yield rational harmonic sums.

## Real Hopf (S⁰ → S¹ → S¹)

Fully degenerate:
- Fiber S⁰ is discrete (two points), has no continuous spectral geometry
- Base = total space (both S¹)
- Gauge group {1} has no coupling constant
- Formally gives 1/α ≈ 15.4 with F = 0

## Conclusions

### Negative Results
1. The quaternionic Hopf recipe gives 1/α ≈ 970, matching **no known weak coupling**
2. The octonionic Hopf gives 1/α ≈ 7164, matching **no known strong coupling**
3. No alternative prefactor or polynomial form produces a match

### Positive Results
1. **Universal selection principle discovered:** B/N = dim(S^d) at m = 2 for all d. This is a new mathematical identity that should be noted in Paper 2.
2. The recipe **formally generalizes** — all ingredients (B, F, Δ) are computable for every Hopf fibration.
3. The representation-theoretic analysis cleanly identifies **why the formula works only for S³**: the SO(4) product structure, the Lie group property, and the unit embedding index. These are three independent reasons, all specific to d = 3.

### Impact on Paper 2
- **Sec III (Selection principle):** Could optionally note the universal algebraic identity B_formal/N = d at m=2 for all d. This STRENGTHENS the paper's argument: the selection principle alone does not distinguish S³ — what's special is the coincidence B_branch = B_formal, which only holds for d = 3 due to the SO(4) product structure.
- **Sec VIII.E (S³ specificity):** The claim "yields no finite n_max for d ≠ 3" is **CORRECT** for the branching-based B (which is the physically meaningful definition). The specificity traces to SO(4) ≅ SU(2)×SU(2) having embedding index 1. No revision needed.
- **Possible addition:** A brief remark that the quaternionic Hopf recipe gives 1/α ≈ 970 (no match to weak couplings) would strengthen the S³ specificity case.
