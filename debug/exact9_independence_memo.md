# Order exactly 9 for the unequal-exponent exchange integral: the independence lemma

Date 2026-11-18. Derivation memo; no paper, test, record or git state touched.
Scripts: `debug/exact9_rank.py`, `debug/exact9_lemma_check.py`. Data: `debug/data/exact9_rank_<case>.json`
(+ `.log`), `debug/data/exact9_lemma_check.json`.

**Verdict: GO.** The nine functions are linearly independent over C(R) for every real 0 < a < b, with no
exceptional exponent ratio; and the coordinate matrix of (K, ..., K^(8)) has rank 9 over C(R) at every real
t = b/a > 1 (exact computation), for the 1s integral and for the four sigma = 0 integrals with l > 0. Hence none of
the five satisfies a linear ODE in R of order < 9 with polynomial (even complex) coefficients, for any 0 < a < b.
The proof is my own and unreviewed; see "where this could be wrong".

## 1. Setting

Nine functions on R > 0 (a < b, c = 1/(2(a^2 - b^2)), LC = ln((a+b)/(b-a))):

| i | g_i | g_i' in the nine |
|---|---|---|
| 1 | f0, defined by f0' = -FM1/R, FM1 = c(g7 - g6 - g9 + g8) (Pachucki arXiv:0905.1635 labels 59, 60) | (c/R)(g6 - g7 - g8 + g9) |
| 2 | Ei(-(a+b)R) | g3/R |
| 3 | e^{-(a+b)R} | -(a+b) g3 |
| 4 | e^{-2aR} | -2a g4 |
| 5 | e^{-2bR} | -2b g5 |
| 6 | e^{2aR} Ei(-2(a+b)R) | 2a g6 + g5/R |
| 7 | e^{2bR} Ei(-2(a+b)R) | 2b g7 + g4/R |
| 8 | psi_a = e^{-2aR}[Ei(2(a-b)R) - 2Ei((a-b)R) - LC] | -2a g8 + g5/R - 2 g3/R |
| 9 | psi_b = e^{-2bR}[Ei(2(b-a)R) - 2Ei((b-a)R) - LC] | -2b g9 + g4/R - 2 g3/R |

The derivative column was written by hand and equals the helper's calculus (`Calc.diff` + `in_nine`) at
t = 2, 3, 5/2 (exact, `connection_matrix_equals_helper...: true`).

## 2. The implication, stated carefully

Let V be the C(R)-span of g_1..g_9 inside the germs of analytic functions at a point R0 > 0. V is closed under d/dR
(table). K = sum_i c_i g_i with c_i in Q(t)(R), specialised at the real t; K^(j) = sum_i M_ji g_i, M_ji in R(R).

* If g_1..g_9 are linearly independent over C(R), then dim_{C(R)} span{K, ..., K^(8)} = rank of the 9 x 9 matrix
  M(t0, R) over C(R), which equals its rank over Q(t0)(R) (rank does not change under field extension).
* An ODE sum_{j <= n} p_j(R) K^(j) = 0, p_j in C[R], p_n != 0, n <= 8, is a nontrivial C(R)-linear relation among
  K, ..., K^(8). So if the rank is 9 no such ODE exists; with the proved bound, the minimal order is exactly 9.

Independence is needed over C(R) only (the coefficient field of the ODE), and that is what Section 3 proves. Two
things are needed, and the second was NOT covered by "formal rank 9 over Q(t)(R)":

(L) independence of g_1..g_9 over C(R) for each real 0 < a < b;
(D) det M(t0, R) not identically zero in R **for each real t0 > 1**. Formal rank 9 only gives this outside the
    real zeros of the R-content of det M; those had to be computed (Section 4).

## 3. Lemma (L) and its proof

**Lemma.** For every real 0 < a < b, g_1, ..., g_9 are linearly independent over C(R).

Standard facts used:

(F1) For real x != 0, Ei(x) = gamma + ln|x| + Ein(x), Ein(x) = sum_{k>=1} x^k/(k k!) entire (principal value for
x > 0; Ei(-x) = -E1(x)). DLMF section 6.6 (power series of Ei and E1; I recall 6.6.1 / 6.6.2 but did not re-open the
page: equation numbers unverified). Checked numerically against the `mp.ei` the helper uses (check C4, 1e-37).

(F2) Exponentials e^{lam R} with distinct real lam are linearly independent over C(R). Proof: in a relation
sum r_lam e^{lam R} = 0 take the largest lam with r_lam != 0; then r_lam = -sum_{lam' < lam} r_lam' e^{(lam'-lam)R}
decays exponentially as R -> +oo, while a nonzero rational function decays at most like a power. So r_lam = 0.

(F3) If P, Q are meromorphic on C and P(R) ln R + Q(R) = 0 on a real interval, then P = Q = 0. Proof: otherwise
ln R = -Q/P on the interval, hence (identity theorem) on the slit plane; the right side is single-valued across the
negative axis, the logarithm jumps by 2 pi i.

**Step 1 (local form at R = 0).** By (F1), for R > 0 each g_i = A_i(R) ln R + B_i(R) with A_i, B_i entire:

| i | A_i |
|---|---|
| 2 | 1 |
| 3, 4, 5 | 0 (B_i = the exponential) |
| 6, 7 | e^{2aR}, e^{2bR} |
| 8, 9 | -e^{-2aR}, -e^{-2bR} (the three logs in each bracket sum to 1 - 2 = -1) |
| 1 | A_1 = -Lam - M(0), Lam(R) = int_0^R L(s)/s ds |

For i = 1: FM1 = L ln R + M with L = c(e^{2bR} - e^{2aR} + e^{-2bR} - e^{-2aR}) and M entire; L(0) = 0, so L/R is
entire and Lam is entire with Lam' = L/R. Integrating f0' = -(L ln R + M)/R by parts,
f0 = -(Lam + M(0)) ln R + (entire). (In fact M(0) = 0 and Lam = 2c[Cinh(2bR) - Cinh(2aR)]; the proof only uses
A_1' = -L/R.)

**Step 2 (split the relation).** Suppose sum r_i g_i = 0 on an interval, r_i in C(R). By (F3)

  (A) sum_i r_i A_i = 0,  (B) sum_i r_i B_i = 0.

**Step 3 (r_1 = 0).** (A) reads r_1 A_1 + r_2 + r_6 e^{2aR} + r_7 e^{2bR} - r_8 e^{-2aR} - r_9 e^{-2bR} = 0. If
r_1 != 0 then A_1 = sum_{lam in S} rho_lam e^{lam R}, S = {0, 2a, 2b, -2a, -2b}, rho_lam in C(R). Differentiate
and use A_1' = -L/R:

  sum_{lam in S} (rho_lam' + lam rho_lam) e^{lam R} + (c/R)(e^{2bR} - e^{2aR} + e^{-2bR} - e^{-2aR}) = 0.

The five rates in S are distinct for every 0 < a < b, so by (F2) the coefficient of e^{2bR} vanishes:

  rho' + 2b rho = -c/R,  rho = rho_{2b} in C(R).

No rational function satisfies this (c != 0): if rho has a pole of order m >= 1 anywhere, rho' + 2b rho has a pole
of order m + 1 >= 2 there, but -c/R has only a simple pole; so rho is a polynomial, and then the left side is a
polynomial while the right side is not. Contradiction; r_1 = 0. (This is Liouville's argument that Ei is not
elementary, in the one-line form needed.)

**Step 4.** With r_1 = 0, (A) is a relation among e^{lam R}, lam in {0, 2a, 2b, -2a, -2b}, distinct; by (F2)
r_2 = r_6 = r_7 = r_8 = r_9 = 0. The original relation is now r_3 e^{-(a+b)R} + r_4 e^{-2aR} + r_5 e^{-2bR} = 0,
rates distinct since a < b; by (F2) r_3 = r_4 = r_5 = 0. QED.

**Coverage.** Every real 0 < a < b; no exceptional ratio. The only inputs are the distinctness of
{0, +-2a, +-2b} and of {a+b, 2a, 2b}, and c != 0.

**Coincidences of rates, listed (none touches the proof).**
* Decay rates at R -> +oo (lam + mu of each key): g2, g3: a+b; g4, g7, the LC and Ei(2(b-a)R) parts of g8/g9: 2a;
  g5, g6: 2b; the Ei((a-b)R), Ei((b-a)R) parts of g8, g9: a+b; f0: 2a. So seven of the nine share a rate with
  another function for EVERY a < b. This is why a comparison of growth classes at infinity does not separate them,
  and why the proof works at R = 0, where the log coefficients A_i carry five distinct rates.
* Ei arguments {-(a+b), -2(a+b), +-2(b-a), +-(b-a)}: pairwise distinct except -(a+b) = 2(a-b) at **b = 3a**, where
  g2 and the first term of psi_a contain the same Ei(-4aR). Harmless here (A_2 = 1, A_8 = -e^{-2aR} remain
  independent); it would need separate treatment in a Kolchin-Ostrowski argument (below).
* b = 2a: no coincidence among any of the sets above. det M != 0 was also evaluated explicitly at
  t = 2, 3, 3/2, 5, 4/3 (R = 7/5), all nonzero.

**A second argument for the key step (r_1 = 0), independent of the pole count.** Lam = 2c[Chi(2bR) - Chi(2aR)
- ln(b/a)], and Chi(x) ~ (e^x/2x) sum k!/x^k as x -> +oo (DLMF section 6.12, numbers unverified). If
A_1 = sum rho_lam e^{lam R}, multiplying by e^{-2bR} and letting R -> +oo gives the asymptotic expansion of the
rational function rho_{2b}, a convergent Laurent series at infinity, equal to a constant times
sum k!/(2bR)^{k+1}, which diverges; asymptotic expansions along a ray are unique. Same conclusion.

**Not done: route (2) / Kolchin-Ostrowski.** The standard alternative (the six Ei's are primitives over
F0 = C(R)(e^{aR}, e^{bR}); Kolchin, Amer. J. Math. 90 (1968) 1151, after Ostrowski 1946: algebraic dependence of
primitives forces a constant-coefficient combination to lie in F0) would give more (algebraic independence) but
needs Liouville's argument inside F0, with the cases a/b rational and b = 3a split off. I did not write it out;
it is not needed. No factorisation of the order-9 operator was attempted.

## 4. Condition (D): rank 9 at every real t > 1 (exact)

`debug/exact9_rank.py CASE`: v_0 = stored closed-form coordinates of K over Q(t)(R) (a = 1, b = t;
K(R; a, b) = a K(aR; 1, b/a)), v_{j+1} = dv_j/dR + v_j Gamma with the hand-written connection matrix; M = rows
v_0..v_8. If t0 > 1 were bad (det M(t0, R) = 0 identically in R) then, since the entries have denominators only
R, t, t - 1, t + 1 (checked), t0 would be a root of the numerator N_j(t) of det M(t, R_j) for every rational R_j.
Computed N_j exactly for R_j = 1, 2, 7/5, 1/3, 5 (sympy DomainMatrix over Z[t], rows cleared of denominators) and
their gcd:

| case | entry denominators | deg N_j | N_j = t^p (t-1)^q (t+1)^s x irreducible of degree | gcd over the five R_j | real roots of gcd in (1, oo) |
|---|---|---|---|---|---|
| 1s | R, t-1, t+1 | 95 | 38, 24, 16; 17 | const t^38 (t-1)^24 (t+1)^16 | 0 |
| 2p_sigma 1s | + t | 151 | 59, 33, 22; 37 | const t^59 (t-1)^33 (t+1)^22 | 0 |
| 3d0 1s | + t | 216 | 91, 42, 28; 55 | const t^91 (t-1)^42 (t+1)^28 | 0 |
| 2p0 2p0 | + t | 219 | 91, 42, 28; 58 | const t^91 (t-1)^42 (t+1)^28 | 0 |
| 2p+1 2p+1 | + t | 210 | 91, 42, 28; 49 | const t^91 (t-1)^42 (t+1)^28 | 0 |

The R-dependent irreducible factor differs between the R_j and drops out of the gcd; the common part vanishes only
at t = 0, +-1. So det M(t0, R) is not identically zero for any t0 outside {0, 1, -1}, in particular for every real
t0 > 1. (The zero at t = 1 is the expected collapse to order 4 at equal exponents.)

Other exact checks in the same script, all five cases: stored coefficients = a fresh solve over Q at (1, 2) by
`tests/_pachucki_extension.py` (`solve_all`, `assemble_K`, `in_nine`); det M != 0 at R = 7/5 for
t = 2, 3, 3/2, 5, 4/3. Run times 7 s (1s) to 220 s (3d0).

## 5. Numerical checks of Step 1 (`debug/exact9_lemma_check.py`, dps 40, t = 2, 3, 3/2, 5/4)

| check | result |
|---|---|
| C1: FM1 (helper's label 60, `mp.ei`) = L ln R + M at R = 0.3, 1.1, 2.7; M(0) | diff <= 8e-29; M(0) = 0 |
| C2: contour integral of f0' once round R = 0 (log continued) vs 2 pi i A_1(0.7), A_1 = -2c[Cinh(2bR) - Cinh(2aR)] | diff <= 5e-40 (size 3.8 to 6.9) |
| C3: B_1 = f0 - A_1 ln R with f0 from the helper's quadrature; B_1(1.3) - B_1(0.05) vs the log-free integral -int (M + A_1)/s | diff <= 2e-40 (size 0.9 to 6.0; the ln term is of the same size, so a wrong A_1 would show) |
| C4: g_2..g_9 = A_i ln R + B_i at R = 0.9 | diff <= 2e-37 |

## 6. Result

**Theorem (given the already-established closed form).** For every real 0 < a < b, the 1s exchange integral
K(R; a, b) satisfies no linear ODE in R with polynomial coefficients of order < 9; with the derived bound its
minimal order is exactly 9. The same holds for (2p_sigma 1s|2p_sigma 1s), (3d0 1s|3d0 1s), (2p0 2p0|2p0 2p0) and
(2p+1 2p+1|2p+1 2p+1) (density x^2 + y^2), given their stored closed forms on the same nine functions.

Tier suggestion (PI's call; nothing edited): "exactly 9" can move from [MEASURED] to a proof resting on (i) the
existing symbolic closed form, (ii) the lemma (hand proof, elementary), (iii) the exact determinant computation of
Section 4. A guard, written separately, should pin: the gcd over >= 2 values of R of the determinant numerator has no
root in (1, oo) (rejects a hidden bad ratio); the connection matrix equals the helper's calculus; and the log
coefficients A_i of Step 1 numerically (rejects a wrong sign in psi's log count, which would make A_8 or A_9 vanish
and break Step 4).

## 7. Where this could be wrong

1. **The lemma is about the nine functions as defined, and f0 enters only through f0' = -FM1/R (label 59/60).**
   If the transcription of label 60 were wrong, the closed form of K would be wrong too; the existing tests compare
   f0 with an independent Neumann evaluation, which I did not re-run.
2. **Step 3 uses c != 0 and the log coefficient L not identically zero.** L = 2c(cosh 2bR - cosh 2aR) vanishes
   identically only at a = b; that is the equal-exponent collapse, outside the claim.
3. **(F3) needs the r_i A_i, r_i B_i meromorphic on C.** They are (entire times rational). The relation is assumed
   on a real interval of R > 0, which is where the ODE lives.
4. **Section 4 trusts the stored closed forms identically in t.** I compared them with a fresh solve only at
   t = 2. If a stored coefficient were wrong at other t, the determinant claim would be about the wrong vector. The
   identity in t is the existing [SYMBOLIC PROOF] claim (residuals identically zero over Q(t)), not re-derived here.
   For l > 0 the scaling K(R; a, b) = a K(aR; 1, b/a) is assumed (normalised orbitals; dimensionally forced).
5. **The gcd argument assumes no entry of M has a pole at the candidate t0.** Entry denominators were factored:
   only R, t, t - 1, t + 1. Clearing row denominators multiplies the determinant by such factors only.
6. **"Exactly 9" is a statement about ODEs with polynomial (rational) coefficients in R at fixed a, b.** It says
   nothing about operators with coefficients in a larger field (e.g. involving e^{R}), nor about m != 0 integrals
   (2p_pi, 3d_delta), whose closed forms are not established.
7. **DLMF equation numbers are from memory** (the facts themselves are textbook and were checked numerically); the
   Kolchin reference is cited only for the route not taken.
8. **Unreviewed.** One author, one pass; the proof is short enough to check by hand and should be before any tier
   change.
