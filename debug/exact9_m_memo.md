# Exact order 9 for the m != 0 exchange integrals: leg (D) for 2p_pi and 3d_delta

Date 2026-11-18. Diagnostic memo; no paper, test, doc, CLAUDE.md, CHANGELOG or git state touched; nothing SEALED
opened; `debug/pachucki_m_*` read only.

## Verdict: GO (both cases, both algorithms)

For 2p_{+-1,A} 1s_B (2ppi) and 3d_{+-2,A} 1s_B (3ddel), the 9 x 9 coordinate matrix of (K, K', ..., K^(8)) on the
nine functions has determinant not identically zero in R for every real t = b/a > 1. With the independence of the
nine functions over C(R) (`debug/exact9_independence_memo.md`, reviewed in `debug/exact9_review_memo.md`), which
does not depend on the integral, this gives minimal order exactly 9 for every 0 < a < b in both cases.

Argument (as for the m = 0 cases): if det M(t0, R) vanished identically for some t0 > 1, then the numerator
N_j(t) of det M(t, R_j) (rows cleared by powers of t, t - 1, t + 1 and R only, none of which vanishes at t0 > 1)
would vanish at t0 for every R_j. So t0 would be a real root > 1 of gcd_j N_j. The gcd has none.

## Inputs and consistency checks

Input: `nine_coords` of `debug/data/pachucki_m_case_{2ppi,3ddel}_sym_res.json` (Q(t) closed form, a = 1, b = t,
order NAMES9). Both algorithms build K in the a = 1, b = t normalisation only.

| check | 2ppi | 3ddel |
|---|---|---|
| hand connection matrix = `tests/_exact_order.hand_table` = helper `Calc.diff`, identically over Q(t) | yes | yes |
| Q(t) form at t = 2 = Q run `pachucki_m_case_*_1_2.json`, all nine coordinates, exactly | yes | yes |
| Q(t) form at t = 15/8, rescaled, = Q run `pachucki_m_case_*_4-5_3-2.json` | yes | yes |
| stored coordinates round-trip through the helper's `in_nine` | yes | yes |
| derivatives K..K^(8) stay in the nine-function span (helper calculus) | yes | yes |

Scaling, as found by the (4/5, 3/2) check: coordinates 1..8 obey c_i(R; a, b) = a c_i(aR; 1, b/a) exactly; the f0
coordinate carries a^{1 - s0} with s0 = -2 (i.e. f0(R; a, b) = a^{-2} f0(aR; 1, b/a) in the normalisation the
coordinates are stored in), the same for both cases. This only rescales the f0 column by a constant and the
argument R by a, so the order is unaffected; nothing in the leg (D) computation mixes normalisations.

## Denominator audit

Entries of rows 0..8 (all coordinates of K, ..., K^(8)), irreducible denominator factors with largest multiplicity:

| case | algorithm 1 (sympy, Q(t, R)) | algorithm 2 (helper, Q(t)) | other factors |
|---|---|---|---|
| 2ppi | t^3, (t-1)^7, (t+1)^7, R^11 | row clearing total t^27 (t-1)^63 (t+1)^63 | none |
| 3ddel | t^7, (t-1)^9, (t+1)^9, R^13 | row clearing total t^63 (t-1)^81 (t+1)^81 | none |

Only t, t - 1, t + 1 occur (and R). No other factor, so no further real roots to examine.

## Determinant numerators

Algorithm 1 (`debug/exact9_m_rank.py`; hand connection matrix over Q(t, R); per-R row clearing; sympy
DomainMatrix det over Z[t]; gcd and factor stripping in python-flint; real roots by arb isolation, real roots
returned exact-real). Algorithm 2 (`debug/exact9_m_review_det.py`; rows from the helper calculus over Q(t); det
mod p = 2^61 - 1 at integer t, Newton interpolation, checked on 25 extra nodes; pairwise gcd of residual factors
mod p). R sets disjoint.

| case | alg | R | degree | mult t | mult t-1 | mult t+1 | residual deg | residual real roots (in (1,oo)) |
|---|---|---|---|---|---|---|---|---|
| 2ppi | 1 | 1 | 146 | 59 | 34 | 22 | 31 | 5 (1) |
| 2ppi | 1 | 2 | 146 | 59 | 34 | 22 | 31 | 7 (1) |
| 2ppi | 1 | 7/5 | 146 | 59 | 34 | 22 | 31 | 5 (1) |
| 2ppi | 1 | 1/3 | 146 | 59 | 34 | 22 | 31 | 1 (0) |
| 2ppi | 1 | 5 | 146 | 59 | 34 | 22 | 31 | 7 (1) |
| 2ppi | 2 | 3/2 | 146 | 59 | 34 | 22 | 31 | -- |
| 2ppi | 2 | 11/7 | 146 | 59 | 34 | 22 | 31 | -- |
| 2ppi | 2 | 9 | 146 | 59 | 34 | 22 | 31 | -- |
| 3ddel | 1 | 1 | 206 | 91 | 44 | 28 | 43 | 7 (2) |
| 3ddel | 1 | 2 | 206 | 91 | 44 | 28 | 43 | 11 (2) |
| 3ddel | 1 | 7/5 | 206 | 91 | 44 | 28 | 43 | 11 (2) |
| 3ddel | 1 | 1/3 | 206 | 91 | 44 | 28 | 43 | 7 (2) |
| 3ddel | 1 | 5 | 206 | 91 | 44 | 28 | 43 | 9 (2) |
| 3ddel | 2 | 3/2 | 206 | 91 | 44 | 28 | 43 | -- |
| 3ddel | 2 | 11/7 | 206 | 91 | 44 | 28 | 43 | -- |
| 3ddel | 2 | 9 | 206 | 91 | 44 | 28 | 43 | -- |

gcds:

| case | algorithm 1, exact gcd over five R | algorithm 2, gcd degree of residuals mod p (pairs) |
|---|---|---|
| 2ppi | 83010348331692982272 t^59 (t-1)^34 (t+1)^22; no real root in (1, oo) | 0, 0, 0 |
| 3ddel | 16524504796857462531765043200 t^91 (t-1)^44 (t+1)^28; no real root in (1, oo) | 0, 0, 0 |

So the residual factors (everything apart from t, t - 1, t + 1) share no root at all across R, by both algorithms;
the gcd is carried entirely by t, t - 1, t + 1. Both algorithms agree on degree and on the three multiplicities at
every R, although they clear denominators differently.

The residual factors do have real roots > 1 at each single R (one for 2ppi at four of five R, two for 3ddel), and
the exact sign scans of algorithm 2 (191 rational t in (1, 10^6] at R = 3/2 and 11/7) see both signs of det. These
are the R-dependent zero curves of det(t, R) in the (t, R) plane, not order drops: no such root is common to two
R values. This is exactly why a single (t, R) evaluation could not decide the question. Spot values det != 0 at
R = 7/5, t in {2, 3, 3/2, 5, 4/3} are listed in the JSON as illustration only.

## Not established / caveats

* The gcd of algorithm 2 is mod p; a spurious common factor mod p is possible in principle (it would show as gcd
  degree > 0, which did not occur), while a true common factor over Q could be hidden only if p divided its
  leading coefficient, which would also put p in the leading coefficients of both residuals; algorithm 1's gcd is
  exact over Z and agrees.
* "Order exactly 9" also needs the operator to exist (order <= 9), shown by the K^(9) closure in
  `debug/pachucki_m_nonzero_memo.md` over Q(t); not re-derived here.
* t < 1 (a > b) is covered by the a <-> b symmetry only insofar as the closed forms respect it; not examined.

## Files written

* `debug/exact9_m_rank.py` -- algorithm 1 + input checks (connection matrix, t = 2, (4/5, 3/2), denominators).
* `debug/exact9_m_review_det.py` -- algorithm 2 (helper calculus, mod-p interpolation, pairwise gcd, sign scan).
* `debug/data/exact9_m_rank_{2ppi,3ddel}.json`, `debug/data/exact9_m_review_det_{2ppi,3ddel}.json`, and the run
  logs `debug/data/exact9_m_rank_*.log`, `debug/data/exact9_m_review_det_*.log`.

Run times: algorithm 1 26 s (2ppi), 78 s (3ddel); algorithm 2 27 s, 60 s.
