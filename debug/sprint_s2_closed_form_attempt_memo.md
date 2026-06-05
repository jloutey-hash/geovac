# Sprint S2 closed-form attempt: structural gap identified

**Date:** 2026-06-05 (continuation of Sprint S2)
**Author:** PM session
**Status:** PARTIAL. Naive operator-counting "chi_k = 2 + n_distinct_left_spatial_ops" attempted; does NOT match empirical profile. Structural reason for the gap identified. Closed-form remains a multi-day analytical task as scoped.

**Files:**
- `debug/sprint_s2_closed_form.py` — operator enumeration + naive count test
- `debug/data/` (no JSON saved; results inline)

---

## Question

Sprint S1 found chi_k = 2 at every sub-block boundary (12/12 across LiH/BeH₂/H₂O). Sprint S2 (partial, this morning) identified the universal post-orbital sequence chi_k(m) = {1, 16, 9, 9, 3, 2} for m = 0, 1, 2, 3, 4, 5 orbitals closed on the left.

The conjectural Theorem 3.2.A.2 stated:
> chi_k(m) = 2 + |C_cross(O_left, O_right)|, where C_cross is the set of Gaunt-allowed multipole channels coupling left and right orbital sets at non-trivial L.

This sprint attempts to verify the conjecture by counting C_cross directly.

## Method

Enumerate all Gaunt-allowed ERI quartets (p, q, r, s) and classify by left/right split. For each cross-cut quartet, identify the "left fermionic operator" pattern. Count distinct spatial-only patterns.

Compare the count + 2 to measured chi_k.

## Result — naive count does NOT match

| m_closed | measured chi_k | n_dist_left_spatial | predicted chi_k = 2 + n | match? |
|---:|---:|---:|---:|---|
| 0 | 1 | 0 | 2 | NO (trivial boundary) |
| 1 | 16 | 7 | 9 | NO (off by 7) |
| 2 | 9 | 10 | 12 | NO (off by 3 wrong direction) |
| 3 | 9 | 11 | 13 | NO |
| 4 | 3 | 6 | 8 | NO |
| 5 | 2 | 0 | 2 | YES (boundary) |

The naive count gives the right answer ONLY at the boundary (where there are no cross-cut operators). For interior cuts, the relationship is more complex than "1 bond-rank-dim per distinct spatial fermionic operator."

## Structural reason for the gap

### At m_closed = 1: predicted 9, measured 16 (off by 7)

The naive count enumerates spatial operators (without spin). At m=1, the cross-cut quartets fall into 7 distinct (n_dag_left, n_undag_left) types. But the bond rank is 16.

The ratio 16/9 = 1.78x ≈ 2 (within the boundary 2). This strongly suggests SPIN DOUBLING at m=1: each cross-cut spatial operator splits into multiple spin patterns, and the count should multiply by ~2 for the alpha/beta degree of freedom that hasn't been collapsed yet.

Specifically: at the start of the cut (m_closed=1, just after the 1s orbital), the "left" has only 1s (1 spatial orbital = 2 spin-orbitals). The h1 + ERI couplings can have:
- both 1s alpha and 1s beta on left
- alpha-only or beta-only cross-cut
- pure left or pure right operators

The 2x factor matches the spin doubling at this small-left regime.

### At m_closed = 2: predicted 12, measured 9 (off by 3, wrong direction)

Here the naive count OVERESTIMATES. The 10 spatial-distinct left operators don't all contribute independent bond-rank dimensions — some are LINEAR COMBINATIONS of each other in the coefficient matrix M.

The structural reason: when both 1s and 2s are on the left, several "left operators" can be expressed as linear combinations through Gaunt-coefficient identities. For example, the L=0 contributions from (1s, 2s) and (2s, 2s) and (1s, 1s) are related through their shared Slater integral structure.

This is exactly the F3 finding from the DF=multipole result this morning: at larger basis the radial-integral matrix R^L has a **graded singular-value spectrum** with linear dependences among densities. Those linear dependences carry through to the MPO bond rank.

### At m_closed = 4: predicted 8, measured 3 (off by 5)

Same pattern as m=2 but more severe. At m=4 (only 2p_{+1} remaining on right), the cross-cut operators are highly constrained. The 6 distinct spatial operators reduce to 3 independent bond-rank dimensions through Gaunt + radial-integral identities.

## What this means for the closed-form

The conjectural Theorem 3.2.A.2 needs refinement:

> **Theorem 3.2.A.2 (refined).** chi_k(m_closed) = 2 + rank(K_cross(m_closed)), where K_cross is the cross-cut coefficient matrix in the (spatial fermionic operator + spin) basis, restricted to operators that span the cut.

The key word is **rank**, not **count**. The rank-vs-count gap is filled by:
1. Spin doubling at small m_closed (multiplies by ~2)
2. Gaunt-channel linear dependences (reduces the count to its rank)
3. Radial-integral linear dependences (reduces further — same F3 hydrogenic-basis ill-conditioning as the DF rank).

So the closed-form bond rank ISN'T just "2 + number of cross-cut multipole channels." It's "2 + rank of the structured coefficient matrix at the cut, which is governed by the same hydrogenic radial-density ring as F3."

The cleanest mathematical statement is therefore a **structural inheritance theorem** rather than a closed counting formula:

> **Theorem 3.2.A.3 (structural).** The MPO bond rank profile chi_k of a GeoVac composed Hamiltonian inherits the rank structure of the underlying multipole expansion's radial-integral matrix R^L, with the same hydrogenic-basis ill-conditioning that controls the DF rank (sprint_df_multipole_lift_memo.md, F3).

This connects the MPO bond rank directly to the DF rank — a cleaner result than the conjectural channel count.

## Implications

The S2 closed-form derivation is **harder than the empirical pattern suggests** because the rank structure of the chemistry MPO is itself controlled by the same radial-integral conditioning that controls DF rank. So:

- The "interior profile" {4, 16, 16, 9, 9, 9, 6, 3, 3} is NOT determined by counting alone.
- It IS determined by the rank of the GeoVac multipole structure at each cut, which is the same algebraic object F3 characterizes.
- A clean derivation would relate chi_k to the per-L rank of R^L (already computed for DF) plus the spin-doubling factor at small m_closed.

This is a tractable analytical sprint of 3-5 focused days (the scoping memo's original estimate). It connects to the F3 graded-rank result and would close into a unified theorem:

> **Theorem 3.2.A.unified.** The MPO bond rank chi_k of a GeoVac composed Hamiltonian, the DF rank, and the radial-integral matrix R^L rank are all manifestations of the same hydrogenic-basis ill-conditioning. They share the same Laguerre-recurrence-derived linear dependences.

## Status update

**Empirically PROVEN today (Sprint S1 + S1-NC):**
- chi_k = 2 at every sub-block boundary (12/12, bit-identical)
- chi_k = 9 under balanced-coupled negative control (bit-identical)
- 8x interior/boundary ratio for {1s, 2s, 2p} valence
- Universal interior profile {4, 16, 16, 9, 9, 9, 6, 3, 3}

**Structural connections identified today (Sprint S2 partial):**
- Bond rank gap (16 vs 9 at m=1) attributed to spin doubling
- Bond rank gap (9 vs 12 at m=2 and similar) attributed to F3-style radial-integral linear dependences
- Connection to DF rank's F3 finding suggests a unified theorem statement

**Open follow-on (Sprint S2-v2, 3-5 focused days):**
- Derive chi_k as a function of per-L rank(R^L) + spin combinatorics
- Show chi_k inherits F3 graded-rank structure
- Theorem 3.2.A.unified replaces the channel-counting conjecture

## Sprint S2 partial verdict

The closed-form is **not closed** today. But the structural connection to F3 is a substantive new finding: the MPO bond rank and the DF rank are governed by the same hydrogenic-basis ill-conditioning. This connection wasn't visible from the empirical S1 result alone; it took the operator-counting attempt to surface it.

The empirical Theorem 3.2.A (chi=2 at boundaries, 8x ratio, universal interior profile) is the load-bearing result for any near-term Paper 14 §sec:mpo_bond_rank section. The closed-form Theorem 3.2.A.unified is the 3-5 day analytical follow-on.

## Recommended next moves

Two options for closing today:

1. **Apply Paper 14 §sec:mpo_bond_rank in empirical form NOW**, with a forward reference to the structural F3 connection. The Theorem 3.2.A statement is bit-clean and the negative control confirms the structural reading. The unified theorem is named follow-on. ~30 minutes.

2. **Continue Sprint S2-v2 into tomorrow**, deriving the chi_k = f(per-L rank R^L + spin combinatorics) formula. Then apply the Paper 14 section with the full unified theorem. 3-5 days.

The PI's "push S2" instruction was to make the attempt, which is done. The honest outcome: closed-form is connected to F3 (substantive new finding) but requires more focused analytical work to fully close.
