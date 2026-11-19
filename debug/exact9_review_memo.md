# Cold review of `debug/exact9_independence_memo.md` (order exactly 9, unequal-exponent exchange)

Date 2026-11-18. Reviewer not involved in the memo. Nothing edited outside `debug/exact9_review_*`.
Scripts: `debug/exact9_review_det.py CASE`, `debug/exact9_review_numeric.py`, `debug/exact9_review_numeric2.py`.
Data: `debug/data/exact9_review_det_<case>.json|.log`, `debug/data/exact9_review_numeric.log` (run hit the cap
after N3 at (1, 5/2); the rest is in `exact9_review_numeric2.json|.log`, `..._numeric2_dps85.json`).

## Verdict: SOUND WITH CAVEATS

I could not break either leg. Lemma (L) is a correct, complete hand proof for every real 0 < a < b. Leg (D) was
reproduced by a different algorithm at different R, all five cases, same root structure. The caveats below are
scope statements, not gaps in the argument; C1 and C2 must be in any paper sentence.

## Attack by attack

| # | attack | result | what I did |
|---|---|---|---|
| 1 | A_i, B_i entire; f0 | held | Re-derived from the helper's keys. FM1 (helper label 60) = c(g7 - g6 - g9 + g8) exactly, key by key (LC signs included). Log coefficient of FM1 is L = c(e^{2bR} - e^{2aR} + e^{-2bR} - e^{-2aR}), L(0) = 0; M(0) = 0 by hand. f0' = -FM1/R is what `Calc.diff` and `PX.Algebra.diff` use (same sign, same FM1). f0 = -Lam ln R + entire, Lam = 2c[Cinh(2bR) - Cinh(2aR)]; all A_i, B_i are entire (Ein, exp, and integrals of entire functions). The integration constant of f0 only moves B_1 and does not enter the proof. Numerically (N1): second difference of f0 - A ln R at h = 1e-4, 1e-5 drifts 1.5e-4 with the memo's A_1, and -9.21, -4.61, -2.30 with the sign flipped, zero, or halved, at t = 5/2, sqrt 2, 3. |
| 2 | same nine functions as the helper | held | `nine_functions` = memo table row by row; LC absorbed identically (psi = e^{..}[Ei - 2Ei - LC]). My own FM1 from the memo's g's vs the helper's Func: 2e-63 (N2). The stored names `-E1(x)` are Ei(-x). |
| 3 | Step 3 | held | Coefficient extraction is F2 with rational coefficients, which is how F2 is stated and proved (dominant rate; a nonzero rational function is not exponentially small). rho' + 2b rho = -c/R: a pole of order m at ANY complex point gives order m + 1 >= 2 on the left; right side has one simple pole; so rho is a polynomial and the left side is a polynomial. Infinity needs no separate case. Uses only c != 0, b != 0. |
| 4 | F3 with poles; Ei of positive argument | held | r_i A_i and r_i B_i are rational times entire, so meromorphic on C; F3 as stated covers that. Ei(x) = gamma + ln x + Ein(x) for x > 0 is the principal value the helper's `mp.ei` returns (N2 confirms on both signs). The relation is continued from an interval to R > 0 by the identity theorem. |
| 5 | implication | held | An ODE of order n <= 8 is a C(R)-relation p among K..K^(8); independence of the g_i gives p M = 0; det M != 0 in C(R) gives p = 0. Real vs complex coefficients is immaterial (rank is stable under field extension; M is real). Covers homogeneous ODEs with rational coefficients only. |
| 6 | leg (D) | held | See below. |
| 7 | numeric collocation | held | See below. |
| 8 | l > 0 | held, with C2 | Stored closed forms lie on the same nine functions (`in_nine` not None, fresh solves). Scaling: any power of a in K(R; a, b) = a^k K(aR; 1, b/a) leaves the ODE order unchanged, so the exponent is irrelevant to the conclusion. |

## Leg (D), recomputed

Differences from `exact9_rank.py`: derivative rows from the helper's calculus over Q(t) (`Calc.diff`, `in_nine`),
no hand connection matrix; determinant numerators by evaluation at integer t modulo p = 2^61 - 1 and Newton
interpolation in pure Python, verified on 25 extra points; R = 3/2, 11/7, 9.

| case | cleared degree | t^p (t-1)^q (t+1)^s | residual degree | gcd degree of residuals, 3 pairs | memo |
|---|---|---|---|---|---|
| 1s | 95 | 38, 24, 16 | 17 | 0, 0, 0 | same |
| 2p_sigma 1s | 151 | 59, 33, 22 | 37 | 0, 0, 0 | same |
| 3d0 1s | 216 | 91, 42, 28 | 55 | 0, 0, 0 | same |
| 2p0 2p0 | 219 | 91, 42, 28 | 58 | 0, 0, 0 | same |
| 2p+1 2p+1 | 210 | 91, 42, 28 | 49 | 0, 0, 0 | same |

A gcd of degree 0 modulo p, with the residual degrees equal to the memo's exact ones (no leading coefficient lost),
means the residuals share no root over C. So det M(t0, R) vanishes identically in R only at t0 in {0, 1, -1}.

Sub-attacks:
* **Hand connection matrix.** "Equal at t = 2, 3, 5/2" does not by itself determine entries with a degree-2
  denominator in t. I compared the memo's table with the helper's calculus over Q(t), identically in t: equal.
* **Row clearing hiding (t - t0).** It multiplies the determinant by a polynomial, so it can only add roots.
  Denominator factors over all 81 x Laurent entries: only t, t - 1, t + 1 (none other, all five cases).
* **gcd logic.** If det(t0, R) = 0 identically, every R-coefficient of the numerator vanishes at t0 (algebraic t0),
  or the numerator is zero (transcendental t0); either way t0 is a root at every R_j. Sound.
* **Stored coefficients at other t.** Fresh solves over Q equal the stored Q(t) coefficients at t = 7/3, 3, 17/10
  (1s, 2p_sigma) and 7/3 (the other three). Exact determinants from those fresh solves equal the stored-route
  determinants at R = 3/2, 11/7 and are nonzero. K(1, 5/2, R = 3/2) closed form minus the independent Neumann
  integral: 3e-34; formal K' minus numerical derivative: 4e-19.
* **A finding the memo does not state.** For the four l > 0 cases det M(t, R0) at fixed R0 changes sign on t > 1
  (exact scan, 371 points, R0 = 3/2 and 11/7; the memo's own json has 1 to 4 real roots per R_j). For 1s it has
  one sign. This does not hurt the proof (the roots move with R), but "det != 0 at R = 7/5" is a spot check that
  can fail at an unlucky t, and a guard must not use a single R.

## Numeric test of the conclusion (attack 7)

Collocation of sum_i r_i(R) g_i(R) = 0, deg r_i <= 3 (36 unknowns), 90 Chebyshev points on [0.15, 8], columns
normalised; positive control = a tenth function R^2 g2 + 3 g6 - R g1 (two null vectors expected).

| (a, b) | dps | smallest singular values | control, smallest four |
|---|---|---|---|
| (1, 5/2), helper f0 | 60 | 5.1e-29, 1.4e-27, 1.0e-25 | 1.7e-62, 1.8e-62, 9.0e-32, 7.6e-29 |
| (1, sqrt 2), own f0 | 50 | 3.3123e-35, 8.7e-32, 1.0e-29 | 2.9e-52, 3.5e-52, 2.2e-35, 4.4e-35 |
| (1, sqrt 2) | 85 | 3.3123e-35, 8.7e-32, 1.0e-29 | 1.6e-87, 2.9e-87, 2.2e-35, 4.4e-35 |

The independent case sits at a fixed value when precision is raised; the control's two null vectors follow the
precision floor. No null vector. (The small fixed values are the natural conditioning of 36 smooth columns; this
is supporting evidence only.) Same test on K, ..., K^(8) at (1, 5/2), degree <= 3: smallest 5.5e-29, no gap.

## Caveats

**C1 (must be in the paper sentence).** The result is conditional on the closed form K = sum c_i(t, R) g_i, i.e. on
the existing [SYMBOLIC PROOF] claim, identically in t. I checked it at five more t values than the memo and against
one independent integral, not identically. Validity at every real t > 1 (not only generic t) follows by continuity
in t, since coefficients have poles only at t = 0, +-1; the memo does not say this and should.

**C2 (must be in the paper sentence).** Scope: homogeneous linear ODE in R with polynomial (equivalently rational)
coefficients, fixed exponents, a != b. For 2p_sigma 1s and 3d0 1s, "0 < a < b" means the l > 0 orbital carries the
SMALLER exponent; the mirrored ordering is not symmetric and is not covered (the gcd has no root in (0, 1) either,
but the closed form is not established there). 1s, 2p0 2p0, 2p+1 2p+1 are symmetric, so a < b is no restriction.
The m != 0 integrals stay uncovered. "Exactly" also needs the existing bound <= 9; for 3d0 this removes the paper's
current "no lower order is not checked".

**C3 (does not block).** The lemma is a hand proof with two reviewers' readings and no machine check. It is 20
lines and elementary; I found no step that needs repair. Wording nits: F2 should say "a nonzero rational function
is not exponentially small" rather than "decays at most like a power"; Step 1 should state M(0) = 0 (true) or
keep -M(0) ln R explicitly, as it does.

**C4 (does not block).** DLMF equation numbers remain unverified by me as well; the series were checked
numerically. The Kolchin reference is not load-bearing.

**C5 (guard design).** See the l > 0 sign-change finding above.

## Wrong answers a regression test should reject

1. A bad ratio: the gcd over >= 2 values of R of the determinant numerators has a root in (1, oo). Plant: multiply
   one stored coefficient vector by (t - 2) -> the gcd must gain the root and the test must fail.
2. Rank 8: K's coordinates replaced by a vector whose f0 (or one psi) coefficient is zeroed -> det identically 0;
   the test must fail. A test that only checks "det != 0 at one (t, R)" does not reject a hidden bad ratio and can
   fail spuriously for l > 0.
3. A wrong connection matrix: one sign flipped in a g_i' row (e.g. -2 g3/R -> +2 g3/R in psi_a') -> must differ
   from `Calc.diff` identically in t, not at three points.
4. A wrong log count in psi: bracket Ei(2x) - Ei(x) instead of Ei(2x) - 2 Ei(x) gives A_8 = 0, so psi_a would be
   log-free and Step 4 fails. Pin A_8, A_9 = -e^{-2aR}, -e^{-2bR} numerically.
5. A wrong A_1: sign flipped, zero, or halved (drifts -9.2, -4.6, -2.3 in N1 against 1.5e-4).
6. A wrong f0 rule: f0' = +FM1/R, or FM1 with the g6/g7 or psi_a/psi_b signs swapped; pin FM1 against the nine
   functions and f0 against an independent integral.
7. A stored coefficient that is right at t = 2 only: compare with fresh solves at >= 2 other ratios, one of them
   t = 3 (the Ei-argument coincidence).
8. A new denominator factor in t other than t, t - 1, t + 1 (would void the gcd argument at its root).
9. Scope creep: a sentence claiming the result for a > b in the mixed cases, for m != 0, or at a = b (order 4).
