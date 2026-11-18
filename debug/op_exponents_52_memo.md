# Why the (AA|BC) operator has 3/2, 7/2, 9/2 at x* and no 5/2

Diagnostic, 2026-11-18. Nothing registered, no paper, record or git touched, no SEALED file opened.
Scripts: `debug/op52_pieces.py`, `debug/op52_subsets.py`, `debug/op52_groups.py`, `debug/op52_na.py`.
Data: `debug/data/op52_{pieces_A1_1073741441,subsets_A1_1073741441,subsets_A2_1073741381,groups_A1,na_1073741441}.json`
(+ `.log`). The mod-p work reuses MO's cached sequences (`s1_coeffs_*`); MO's frozen engine is unchanged.

## 0. The counting fact used throughout

Near x* every solution is u + (x* - x)^{1/2} v, u and v analytic. The half-integer exponents are 1/2 + ord(v) over a
valuation-adapted basis of the space H of v's: dim H = their number, and 5/2 means some element of H vanishes to order
exactly 2. One function contributes one exponent however many terms its branch has. The generic sum's H is the sum of
the pieces' H's.

## 1. Which components carry which exponent (Q1)

Five of the 21 constants carry x*, all exponentials (no E1 or log piece does). A1, prime 1073741441; A2 (prime
1073741381) gives the same orders and exponent sets for the subsets.

| constant (A1) | endpoint | fed by kernel | piece order | half-integer exps at x* |
|:--|:--|:--|:--:|:--|
| e^{-zB D1} = e^{-117/50} | D1 | NA only | 5 | 3/2 |
| e^{-zB D1 - mu(D2-D1)} = e^{-189/50} | D1 | Y, dY | 6 | 3/2 |
| e^{-zB D1 - mu(D2+D1)} = e^{-513/50} | D1 | Y, dY | 6 | 3/2 |
| e^{+zC D1 - (zB+zC) D2} = e^{-209/50} | D2 | NA, Y, dY | 10 | 7/2 |
| e^{-zC D1 - (zB+zC) D2} = e^{-389/50} | D2 | NA, Y, dY | 10 | 7/2 |

Generic sums over subsets (half-integer exponents only):

| subset | A1 | A2 |
|:--|:--|:--|
| the two Y/dY D1 pieces | 3/2 | 3/2 |
| D1 trio (NA + Y/dY) | 3/2, 7/2 | 3/2, 7/2 |
| D2 pair | 7/2 | 7/2 |
| NA's D1 piece + D2 pair | 3/2, 7/2 | 3/2, 7/2 |
| Y/dY D1 pieces + D2 pair | 3/2, 7/2 | 3/2, 7/2 |
| all five (= all 21) | 3/2, 7/2, 9/2 | 3/2, 7/2, 9/2 |

So dim H = 3 with valuations {1, 3, 4}: the D1 trio gives a 2-dimensional H (valuations 1 and 3), the D2 pair a
1-dimensional H of valuation 3, and together a valuation-3 collision resolved at 4.

**The hypothesis partly holds.** 3/2 comes from point-charge (1/s) components, but exponential ones, not E1. 7/2 is
the p = 3 term. 9/2 is not a p = 5 component but the second of two valuation-3 branches. Nothing carries 5/2.

**The one-electron integral (fresh check, prediction written in the script docstring before running).** Same
endpoint split for record NA's object at MO's N3 geometry: D1 piece e^{-zB D1}: 3/2; D2 pair: 3/2 (one-dimensional);
D1 + pair: 3/2, 5/2; all six: 3/2, 5/2. Held as predicted. In NA both endpoints carry the point charge, and the two
point-charge branches differ at the next order: that is NA's 5/2.

## 2. Component-level derivation (Q2)

**Where the constants come from (derived from the engine's structure).** With B at the origin, the radial
integral splits at r = D1 (the sphere through C) and r = D2 (the sphere through A). The exponents of the five constants
identify their endpoint: D1-endpoint pieces carry e^{-zB D1} times the kernel's value at D1 (NA: no exponential;
Y and dY: i_L(mu D1) k_L(mu D2), split into e^{+-mu D1} e^{-mu D2}); D2-endpoint pieces carry the density at D2,
e^{-(zB+zC) D2} e^{+-zC D1}, times the kernel's jump at D2. The jump of a Yukawa kernel at D2 is a Wronskian of
i_L and k_L, so its exponential factors cancel (e^{+mu D2} e^{-mu D2}): it lands in the same mu-free constants as NA's
jump. That is AB's "record defect 1" (the D2 pair receives NA, Y and dY), and the reason the e^{-2 mu D2} keys are
unused.

**What each endpoint sees.** Write V_A = 1/s - S with S = (1/s + zA) e^{-2 zA s}. The odd powers of s, which are
the ones that make a branch, are (sympy):

| potential | s^-1 | s^1 | s^3 | s^5 |
|:--|:--|:--|:--|:--|
| 1/s (NA) | 1 | 0 | 0 | 0 |
| S (screened: Y with dY) | 1 | **0** | -2 zA^4/3 | -8 zA^6/45 |
| V_A = 1/s - S | 0 | 0 | 2 zA^4/3 | 8 zA^6/45 |

S's s^1 coefficient is mu^2/2 - zA mu = 0 at mu = 2 zA: the bounded-density condition, one term at a time.
A potential term s^p meeting C's cusp |r - C|^q gives the branch (x* - x)^{(p+q+3)/2} (L52 design memo).

- D1 endpoint: the kernel enters by value, so the piece sees one potential at a time. NA's piece sees 1/s,
  leading p = -1, q = 1: 3/2. The Y/dY pieces see S (Y and dY share these constants with rational weight zA, so they
  never separate): leading p = -1: 3/2. The two branches have the same p = -1 term (coefficient 1 in both), so their
  difference is V_A's branch, leading p = 3, q = 1: valuation 3, exponent 7/2. A valuation-2 element would need a
  potential in span{1/s, S} with a nonzero s^1 term, or a second p = -1 branch with a different C-side profile.
  Neither exists: both have zero s^1 term, and the p = -1 parts are identical.
- D2 endpoint: the jump of the whole kernel is V_A's (all three kernels in the same constants), leading p = 3: 7/2.
  One branch only.
- Union: two valuation-3 branches (D1's V_A branch and D2's V_A branch). One combination of them vanishes to
  order 4: 9/2. In NA the same step at valuation 1 gives 5/2.

So the rule reads: H is spanned by one branch per (endpoint, distinct potential seen there); the valuation of each is
(p_min + 3)/2 (q = 1, exponent (p_min + 4)/2); two branches with equal leading term differ at the next order (+1). With p_min = -1
at D1 and 3 for V_A this gives {1, 3, 4} for (AA|BC) and {1, 2} for NA.

**Status of each step.**
- Derived: the counting fact (sec. 0); the endpoint assignment of the five constants (exponent arithmetic of the
  engine); that the D2 jump is mu-free (Wronskian); the s-expansion table and v_1(S) = 0; the branch exponent
  (p + q + 3)/2 for one term (L52 memo, conditional as stated there).
- Matched, not derived: that each (endpoint, potential) contributes exactly one branch dimension (the two
  e^{+-mu D1} halves of Y/dY at D1 collapse to one, and the D2 pair's two e^{+-zC D1} halves collapse to one;
  measured, not proved); that the D1 and D2 valuation-3 branches differ at exactly the next order (4, not higher); that
  no other piece adds to H (measured: the 16 other pieces have no x* point).
- The table explains why no 5/2 appears; it does not prove that 9/2 must be present (a coincidence of the D1 and D2
  branches through order 4 would push it to 11/2; not excluded by the argument).

## 3. Fresh check (Q3)

**The planned amplitude test is not possible.** Over Q (A1, L <= 200, 1500 digits) the components are factorially
divergent formal series (|q_{L,c} c| ~ 10^679 at L = 180); the trio's sum is cancelled by the other 16 constants. A
per-component 3/2 amplitude exists only at operator level.

**What could be measured: the D2 pair.** Its two halves are each divergent (opposite, equal to shown digits) but their
constant-weighted sum converges. The K6 instrument on it: beta_eff = -4.162, -4.086, -4.059, -4.045, -4.040 at
L = 40, 80, 120, 160, 180, monotone toward -4: the pair's own branch leads with 7/2, no 3/2 or 5/2, as the
mechanism says (the D2 endpoint sees only V_A). Written before running and **missed**: that the pair carries all of
f's branch. Its 7/2 amplitude is about -760 times f's (ratio -565, -682, -728, -753, -761 at L = 40..180, still
drifting); f's small 7/2 coefficient is a near-cancellation between the D2 pair and the remaining constants. f itself:
beta_eff(180) = -4.129 (AB: -4.106 at 220).

## 4. For the PI

**Paper 59.** The current sentence ("x* has no 5/2 exponent") could gain one clause: the operator's half-integer
exponents come from the radial endpoints. At the sphere through C the point charge and the screened potential each
give 3/2, and their difference is the smooth one-centre potential (no s^1 term), which gives 7/2. The sphere through A
sees only that difference. That is the same reason f itself has no 3/2 or 5/2. Tier: [ARGUMENT + MEASURED]
(endpoint bookkeeping derived; one dimension per (endpoint, potential) measured at A1, A2, N3 on one prime each).
It would also explain record NA's 5/2 (two point-charge endpoints) and close AB's record defect 1 at operator level.
Backing would be a subset test on cached sequences (minutes), written and fire-tested separately (§9).

**Record.** Drafted as OP52 on the fresh object f' = NA - Y_mu (rho_A ~ e^{-mu s}/s): `docs/predictions/2026-11-18_aa_bc_operator_exponents.md`.
