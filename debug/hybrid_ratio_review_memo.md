# Cold review of `debug/hybrid_ratio_memo.md` (minimal order of four L >= 1 hybrids at every exponent ratio)

Date 2026-11-19. Adversarial reviewer; did not write the memo under review. Nothing under review edited; no paper,
test, doc, CLAUDE.md, CHANGELOG or geovac edit; no git state change.

Scripts: `debug/hybrid_ratio_review_toy.py` (attack 1), `debug/hybrid_ratio_review_fresh.py` (2, 7b),
`debug/hybrid_ratio_review_zero.py` (4, 5), `debug/hybrid_ratio_review_points.py` (3, 6, 7, quadrature).
Data: `debug/data/hybrid_ratio_review_{toy,fresh,zero,points}.json`.

## Verdict: SOUND WITH CAVEATS

Every attack held. The Proposition is correct, and I re-proved it (below). The parametric forms agree with the
unpatched evaluator at 15 fresh (case, region, u) points, from u = 1e-6 to u = 1000, to 59-92 digits. The memo's
own structural claims (proportional E-coefficients, constant zero-gcds, 3 surviving rates in L f) come out the same
when I recompute them with independent code. At two coincidences and two generic points, at exponents not used
before, the existing certified machinery gives orders 4/4/5/5. The caveats are about scope and wording. None is a
defect in the mathematics.

| # | attack | result |
|---|---|---|
| 1 | Proposition order = d_E + k: re-proof, plus 8 toys against a direct rank | tried, held |
| 2 | parametric form at fresh u vs the unpatched evaluator | tried, held |
| 3 | the +a_d atom is e^{a_d R}E1((a_c + a_d)R) | tried, held (memo's correction right) |
| 4 | q's R-coefficients have a nonzero rational gcd | tried, held |
| 5 | L f's three zero-gcds are nonzero rationals | tried, held |
| 6 | order 4 at a coincidence (fresh exponents, exact machinery) | tried, held |
| 7 | scaling; the dropped normalisation is harmless | tried, held |

## 1. The Proposition (attack 1)

**Re-proof.** Let W be the free C(R)-module on the d/dR-closed atom set, with atoms independent over C(R) by Lemma GI.
For a fixed u, M = C(R)[d/dR] f has dimension equal to the minimal order. W_exp is a submodule. That gives the exact
sequence 0 -> M cap W_exp -> M -> pi(M) -> 0.

(i) In W/W_exp an E-atom e^{cR}y behaves like e^{cR}, because (e^{cR}y)' = c e^{cR}y + (rational)·(exp atom). So
W/W_exp is the direct sum over c of (C(R)e^{cR})^{m_c}. Different c give non-isomorphic rank-one, hence simple,
modules, since e^{(c-c')R} is not in C(R). Every submodule N of a semisimple module is the direct sum of its
intersections with the isotypic components, and each intersection equals N's projection onto that component.
Applied to N = pi(M), this gives pi(M) = sum over c of cyc(v_c). Twisting by e^{-cR} turns the c-component into
C(R)^{m_c} with the trivial connection, and cyc(q_1, ..., q_m) is then spanned by the derivative vectors. Its rank
is the rank of the Wronskian matrix, which is dim_C span{q_j}. Two atoms at one shift cause no problem: I checked
d_E = 2 coming from a single shift with independent coefficients, and d_E = 2 coming from three dependent
coefficients (toys T4, T5).

(ii) M cap W_exp = {P f : P pi(f) = 0} = (D·L) f = cyc(L f). Here Ann(pi f) = D·L0, because D is a left PID.
The operator L has order d_E and annihilates every q_j e^{c_j R} coefficient function, so it annihilates pi f.
That makes L = r·L0. The dimension of cyc(L f) is the number of rates (Theorem A). QED. Step (ii) is fine as
written.

**Toys** (`hybrid_ratio_review_toy.py`). The direct order is the rank over Q(R) of (f, ..., f^(N)) in the closed
atom basis, using the tests' own atom calculus. I built L by my own cofactor expansion. All 8 toys agree:

| toy | direct | d_E + k |
|---|---|---|
| T1 q e^{-R}[Ei(R) - Ei(2R)] + R e^{-R} + e^{-3R} | 5 | 1 + 4 |
| T2 same with r = 5q (exponential coefficient in ker L) | **4** | 1 + 3 (rate -1 drops) |
| T3 exp at a closure rate | 3 | 1 + 2 |
| T4 independent coefficients at one shift | 5 | 2 + 3 |
| T5 coefficients q, R, q + R at one shift | 5 | 2 + 3 |
| T6 mirror pair (two shifts), q e^{-R} in kernel, closure rates merging | 5 | 2 + 3 |
| T7 ln R + Ei at one shift | 3 | 1 + 2 |
| T8 exponential at c + kappa (near-miss cancellation) | 3 | 1 + 2 |

T2 is the interaction the task asked about. An exponential coefficient that lies in L's kernel does drop out, and
the count records the drop because k is read from L f itself.

## 2. Fresh-u check (attack 2)

`hybrid_ratio_review_fresh.py` compares the stored atom form at u values outside the memo's 27 points against the
current `hybrid_closed_form`. Atoms are evaluated at 100 digits, at R = 1/7, 6/5, 7/2, 11. I made the test sharper
than the memo's: the ratio must equal **the product of the four `radial_norm`s exactly**, not just some constant.

| case / region | u | max relative deviation |
|---|---|---|
| fam1 mid | 1e-6, 3/7, 1000 | 9e-88, 7e-92, 1e-91 |
| H2 gt | 1e-5, 11/3 | 7e-59, 1e-91 |
| H4 lt | 1e-5, 2/7, 50 | 2e-69, 4e-88, 3e-81 |
| H7 mid | 1e-4, 9/4, 300 | 1e-88, 6e-92, 8e-92 |

The lower figures near u -> 0 come from cancellation in the 100-digit atom sum close to a coincidence. They are
still far beyond 30 digits. Every LOG constant's argument was positive at each u tested (asserted).

Branch risk: an undecided `is_negative` in `e1_moment`/`finite_power_exp` falls to E1, which the collector maps
to the same real continuation -Ei(-kR); `_integrate_with_e1`'s `alpha > 0` raises rather than guesses; a wrong
`_ein` branch would put a log of a negative factor into the constants (asserted against at every fresh u).

Independent route: `hybrid_quadrature` (2-D numerical quadrature) against the closed form at fresh points gives
H7 (2, 7/3), R = 1.3: 3.5e-12; H4 (2, 3/2), R = 1.7: 6.6e-13.

## 3. The +a_d atom (attack 3)

From `hybrid_closed_form` output, parsed with `tests/_minimal_order.Fn`:

| case, Z | (a_p, a_c, a_d) | +a_d shift arguments | E1((a_c+a_d)R) | E1((a_p+a_d)R) |
|---|---|---|---|---|
| H2 (2, 1) | (2, 1, 1/2) | E1(3R/2), E1(7R/2) | yes | no |
| H4 (2, 3/2) | (2, 1, 1/2) | E1(3R/2), E1(7R/2) | yes | no |
| H7 (3, 1) | (5/2, 3/2, 1) | E1(5R/2), E1(5R) | yes | no (would be 7/2) |

The memo's correction is right: the earlier memo's E1((a_p + a_d)R) is a typo.

## 4-5. Zero tests (attacks 4, 5)

`hybrid_ratio_review_zero.py` reads only the stored form. It does not use `hybrid_ratio_analyse.py`; it builds L by
its own cofactor expansion with sympy Poly gcds.

* **fam1 mid.** At each shift the coefficients are proportional with ratios (1, -1), and q is free of constants.
  The gcd of q's R-coefficients is **72** at both shifts. q_+(R)/q_-(-R) = -1 and the shift sum is 0. No E-atom is
  left in L f. The three rates -u-1, -2u-2, -2u-1 have gcds **497664, 62208, 497664**, with no constants and no
  positive roots. Every u-factor of the denominators is in {u, u+1, 2u+1, 3u+2, 4u+3}. These are the memo's
  numbers. Order 5.
* **H7 lt.** q_+(R)/q_-(-R) = **+1** (memo's eps = +1 for H7 confirmed). k = 3, all gcds nonzero rationals (for
  example 74037208411275264), with denominator factors among {u, u+1, u+2, 8u+5, 8u+11}. Order 5.
* **H4 eq1.** L f has two rates, -3 and -1, with gcds 729 and 243. Order 4.

## 6. Coincidence order at fresh exponents (attack 6)

`tests/_minimal_order.py` machinery: an interval-certified r x r minor nonzero for the lower bound, the generic rank
over Q(R) for the upper bound, and every (r+1)-minor tried vanishing.

| case | Z_A, Z_B | (a_p, a_c, a_d) | where | order |
|---|---|---|---|---|
| H7 | 7, 7/2 | (35/6, 7/2, 7/2) | a_d = a_c | **4** |
| H7 | 7, 13/5 | (35/6, 7/2, 13/5) | lt | 5 |
| H4 | 4/3, 6 | (4/3, 2/3, 2) | a_d = a_p + a_c | **4** |
| H2 | 5, 17/3 | (5, 5/2, 17/6) | mid | 5 |

## 7. Scaling and normalisation (attack 7)

f(R; sZ_A, sZ_B) = s f(sR; Z_A, Z_B) holds exactly (deviation 0 at 45 digits) for fam1 (3, 5), s = 3/2, and for
H7 (2, 7/3), s = 5/4. That is what the change of variables predicts for normalised orbitals.

Dropping the normalisation is harmless. Each radial normalisation is a positive rational times Z^{l+3/2}, and Z is
a positive expression in u, so the dropped factor is a nonzero constant in R for every u > 0. Section 2's
exact-norm-product match confirms it numerically.

## Required caveats (for the PM)

1. **Conditional statement.** The result is a theorem conditional on (a) Lemma GI, a hand proof reviewed once with
   caveats, and (b) `hybrid_closed_form_shell` being the integral. Point (b) is validated against independent
   quadrature at points (here to ~1e-12), not proved. The memo's Sec. 5 "now a theorem" should carry both
   conditions in the paper.
2. **The Proposition is a new hand proof.** It held under this review (re-proved, 8 toys). If it is cited in Paper
   58 it needs a backing test. The toy comparison is a natural basis: the quotient count against the direct rank,
   including the kernel-interaction case T2.
3. **Scope is four hybrids, not the class.** H1, H3, H5, H6, H8 and H10 hit the run cap in the build step, with no
   partial result. H9 is class A. Those six keep "proved only at the exponents computed".
4. **The zero tests are exact but run in unreviewed collector code** (`hybrid_ratio_rep.py`). My independent
   recomputation from the stored form covers L f and q, not the collector. The collector is covered by the
   numerical match to ~1e-90 at 15 fresh points: a collector error in one coefficient would show up as an
   R-dependent ratio.
5. Housekeeping, not an error: `geovac/two_center_eri.py` carries an uncommitted `radial_norm` / nsimplify fix,
   and every check, theirs and mine, ran against the working copy. The order cannot depend on it (it changes a
   constant factor only), but the norm-product match in Sec. 2 is against the fixed version.

## Proposed minimal Paper 58 wording

Replace "least order proved only at the exponents computed" with:

> For four hybrids, (2p0 2p0|1s 1s_B), (2p1 2p1|2p1 2p1_B), (2p0 2p0|2p0 3d0_B) and (2p0 3d0|2s 1s_B), the least
> order is proved at every Z_A, Z_B > 0. It is the number of distinct rates plus one: 5 away from the two
> coincidences a_d = a_c and a_d = a_c + a_p, and 4 at each, so each coincidence lowers the order by exactly one
> and there is no other exceptional ratio. The proof assumes the atoms' independence (Lemma GI) and the closed
> form. For the other six panel hybrids with an L >= 1 pair, the least order is proved only at the exponents
> computed.

Tier: proved conditional on GI. Record `rests on: Lemma GI` in the claim-test matrix row.

## Wrong answers a guard should reject

* order 4 at a generic ratio, or 5 at a_d = a_c or a_d = a_c + a_p, at a fresh exponent (table in Sec. 6);
* d_E = 4, treating the two E-atoms at one shift as independent lines (would give 7);
* the Proposition applied without the kernel term, i.e. k counted from f's exponential rates rather than from L f.
  This would miss the T2 drop: toy T2 gives 5 instead of 4;
* E1((a_p + a_d)R) as the +a_d atom (absent at H2 (2, 1), H4 (2, 3/2), H7 (3, 1));
* a normalisation-free form whose ratio to `hybrid_closed_form` is constant in R but is not the product of the
  four `radial_norm`s (a guard testing only "constant ratio" accepts a form off by any u-dependent factor);
* a log or gamma constant surviving into L f, or a zero-gcd of L f (or of q) with a positive root in u;
* eps for H7 given as -1 (it is +1 here, H7 lt; fam1 mid is -1; H4 not independently checked).
