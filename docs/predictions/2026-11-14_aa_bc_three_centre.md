# Pre-registered prediction: the three-centre Coulomb integral (AA|BC) in the bond angle (AB)

**Registered 2026-11-14** (chemistry period inventory, row 6), before any exact constant split, operator search or
class test is run on (AA|BC) or on its screened pieces.
**Commit status:** draft; to be committed on its own, before any computation it scores.
**Rests on:** Paper 59 sec:f12 (the density-topology rule, an OBSERVATION: zero or one distinct two-centre density
gives genus zero; every tested instance is two-centre); sec:scope (the one-electron three-centre integral NA, order-8
bond-angle operator, x* its only finite branch point, exponents 3/2 and 5/2); records NA
(`2026-10-28_nuclear_attraction_angle_operator.md`) and IC (`2026-11-01_integrated_class.md`, where NA at case N1 is
NOT TRIANGULAR); `docs/walls/register.md` (genus zero does not mean "lands in {E1, ln, gamma}": K0 is genus zero and
non-Liouvillian); `geovac/two_center_eri.py` V_L (the potential of a 1s density).
**Guardrail check (Papers 8-9, FCI-M, Track DF):** not in scope; no molecular encoding is proposed.

## Object

(AA|BC) = int int rho_A(r1) phi_B(r2) phi_C(r2) / r12, normalized 1s Slater functions, zeta_A, zeta_B, zeta_C, three
distinct centres. B at the origin, C at (0, 0, D1), A at D2 (sin g, 0, cos g); x = cos g = cos(angle C B A).
f(x) = sum_L F_L P_L(x). In NA's notation X = B (zX = zeta_B), Y = C (zH = zeta_C), Z = A.

**The reduction.** rho_A = (zeta_A^3/pi) e^{-2 zeta_A r} has the potential V_A(s) = 1/s - (zeta_A + 1/s) e^{-mu s},
mu = 2 zeta_A. So, with s = |r - A| and rho_BC = phi_B phi_C,

    (AA|BC) = NA(x) - Y_mu(x) - zeta_A W_mu(x),   NA = int rho_BC / s,  Y_mu = int rho_BC e^{-mu s}/s,
                                                  W_mu = int rho_BC e^{-mu s} = -dY_mu/dmu.

NA is the record-NA integral with the point charge at A. Y_mu has the Legendre kernel
mu sum (2L+1) i_L(mu r<) k_L(mu r>) P_L (r<, r> of (r, D2); i_0 = sinh x / x, k_0 = e^{-x}/x, the na_exact convention).

## Cases

| Case | role | D1 | D2 | zeta_B | zeta_C | zeta_A | mu | x* | rho = 1/(x* + sqrt(x*^2 - 1)) |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| A1 | scored | 9/5 | 13/5 | 13/10 | 1 | 9/10 | 9/5 | 125/117 | 9/13 |
| A2 | scored (P1-P3) | 3/2 | 12/5 | 1 | 6/5 | 4/5 | 8/5 | 89/80 | 5/8 |

A1's Coulomb piece is NA's case N1 and A2's is N2, unchanged. zeta_A is chosen so that no exponential rate below
vanishes or coincides with the Coulomb rate zeta_B + zeta_C (an earlier choice, zeta_A = 1 at A1, made
zeta_B - zeta_C - mu = -(zeta_B + zeta_C) and was dropped). Referee points: x = 1/2 (A1), x = -1/3 (A2).

## Controls measured before registering

Drivers `debug/aabc_controls.py` (referee, decay, candidates) and `debug/aabc_na_posctrl.py`; data
`debug/data/aabc_controls.json`, `debug/data/aabc_referee.log`, `debug/data/aabc_decay.log`,
`debug/data/aabc_na_posctrl.json`, `debug/data/aabc_intclass_ctrl_1073741833.log`.

| Control | what | measured |
|:--|:--|:--|
| K1 fit floor | even-tempered n-Gaussian least-squares fit of e^{-r}; closed forms (AA\|AA) = 5 zeta/8 and the equal-exponent two-centre (AA\|BB) at R = 2 | n = 22: 4.9e-12, 1.4e-11; n = 26: 3.0e-13, 2.6e-12 |
| K2 referee, R2 (independent: Gaussian fits, closed-form s-Gaussian ERIs, no use of V_A) | (AA\|BC) at A1, x = 1/2 | n = 18: 0.21286449875185; n = 22: 0.21286449877733; n = 26: 0.21286449875979 |
| | at A2, x = -1/3 | n = 18: 0.22313684936064; n = 22: 0.22313684950325; n = 26: 0.22313684949202 |
| K3 referee, R1 (prolate quadrature of rho_BC against V_A; uses the reduction's potential) | A1, grids 60^2x72 / 90^2x108 / 120^2x144 | 0.2128644844, 0.2128644959, 0.2128644972 |
| | A2, same grids | 0.2231368466, 0.2231368520, 0.2231368505 |
| K4 R1 vs R2 (n = 26), finest grid | A1; A2 | 1.5e-9; 1.0e-9 (R1 converges slowly: V_A is only C^2 at A) |
| K5 NA positive control (`na_exact`, `na_modguess`, unchanged; N1, L <= 1000, weights as na_run.py) | smallest degree at orders 7, 8, two primes | both primes (2^61 - 1; 4611686018427387847): order 7 none (row limit Q = 106); order 8 at Q0 = 31, nullity 1; six constants as listed under P2; 263 s, 224 s |
| K6 decay instrument, positive control on NA (N1, exact rationals, 2000 digits) | beta_eff(L) = log(G_{L+1}/G_{L-1}) / log((L+1)/(L-1)), G_L = F_L / rho^L | N1: L = 40: -2.142; 100: -2.062; 160: -2.040; 220: -2.029 (-> -2, NA3's 3/2 branch); one sign for L >= 5. N2 (A2's Coulomb piece, rho = 5/8): -2.100, -2.044, -2.028, -2.021, identical at 2000 and 2600 digits (`debug/aabc_decay_n2.py`) |
| K6' precision floor | the same at 400 digits | garbage from L = 120 (beta_eff = +1208, mixed signs): the exact split cancels about 1800 digits by L = 240 |
| K7 IC instrument (`intclass_controls.py`, unchanged) | all 16 controls, prime 1073741833 | ALL PASS, 48 s |
| K8 candidates (rule below, evaluated) | rates, E1 and log arguments; vanishing rates | none vanish at A1 or A2 (listed under P2) |

R2's spread between n = 22 and 26 is 1.8e-11 (A1) and 1.1e-11 (A2), so the referee values are
**A1: 0.2128644988 (+-3e-11); A2: 0.2231368495 (+-3e-11)**.

## Predictions (derived before computing; none uses an exact split or operator of Y_mu, W_mu or (AA|BC))

- **AB0 (controls; run first).** K5 reproduces (R, Q0) = (8, 31) with nullity 1 on both primes; K7 passes on both
  IC primes 1073741789 and 1073741833. Falsifier: either fails. Probability 0.97.

- **P1 (the reduction, exact).** The exact split of NA - Y_mu - zeta_A W_mu (Phase-2 engine `debug/aabc_exact.py`,
  na_exact's integration-by-parts and finite-part rules with the Yukawa kernel), evaluated at 2000 digits and summed
  to L = 240 at the referee point, matches the R2 referee to **<= 1e-10** in both cases, and R1 to <= 3e-9.
  Falsifier: a larger difference in either case. Probability 0.85.
  Reported: the per-L agreement with an independent per-L radial quadrature of the screened pieces (target 1e-40).

- **P2 (finite, l-independent constant set).** For L = 2..1000 the set of constants with nonzero rational
  coefficient in F_L / (2 pi N) is the same at every L, and it is contained in the candidate set C built by the
  negative-power rule: an E1 or log constant can arise only where the integrand's e^{b r} component carries a negative
  power of r (the outer Coulomb kernel r^-(L+1), k_L of the moving radius, and, new here, i_L(mu r) split into
  e^{+-mu r}). This rule reproduces NA's measured six constants (K5): e^{-zeta_C D1}, e^{-zeta_B D1},
  e^{+-zeta_C D1 - (zeta_B + zeta_C) D2}, e^{+-zeta_C D1} E1((zeta_B + zeta_C) D2). C at A1 (A2 in brackets):
  - NA's six, unchanged;
  - E1 at D1, arguments -b D1 over the inner rates b = -zeta_B +- zeta_C +- mu: -27/10, 9/10, 189/50, 369/50
    [-27/10, 9/10, 21/10, 57/10], prefactors e^{-zeta_C D1 - mu D2} and, for the two rates shared with the middle
    region, e^{+-zeta_C D1 - mu D2};
  - E1 at D2, arguments (zeta_B + zeta_C - mu) D2 = 13/10 [36/25] with e^{+-zeta_C D1 - mu D2}, and
    (zeta_B + zeta_C + mu) D2 = 533/50 [228/25] with e^{+-zeta_C D1 +- mu D2};
  - e^{-zeta_C D1 - mu D2} ln|b| for the four inner rates, |b| = 1/2, 3/2, 21/10, 41/10 [3/5, 7/5, 9/5, 19/5];
  - exponentials e^{e0} with e0 the endpoint values of these components.
  Sub-claims: (P2a) the set is L-independent and inside C; (P2b) gamma does not appear (its coefficient is minus the
  sum of the residues, zero because the inner-region integrand is analytic at r = 0); (P2c) at least one log appears
  (new against NA, which has none); (P2d) the only mu-free constants are NA's six, and the two E1 constants
  e^{+-zeta_C D1} E1((zeta_B + zeta_C) D2) receive no contribution from Y_mu or W_mu.
  P2 holds if P2a, P2b and P2d hold. Probabilities: P2a 0.85, P2b 0.8, P2c 0.6, P2d 0.85.

- **P3 (singular points).** (P3a) The minimal operator of the generic sequence (NA's search, two primes, L <= 1000,
  R <= 24) has x* as its only finite point with a non-integer exponent, at both cases, and no finite point whose
  position depends on mu: the screened pieces add no Landau point (e^{-mu s}/s is analytic for s != 0). Its
  exponents at x* include 3/2 and 5/2 (inherited: see "Already known"). Probability 0.7.
  (P3b, the new content) **The function f has no (x* - x)^{3/2} or ^{5/2} branch: its leading branch is
  (x* - x)^{7/2}.** Reason: V_A's Fourier transform is (4 pi/k^2)(1 + k^2/(4 zeta_A^2))^-2 ~ 64 pi zeta_A^4 / k^6,
  against the k^-4 cusp of rho_BC at C; k^-10 is |A - C|^7, and |A - C|^2 = 2 D1 D2 (x* - x). The point-charge
  term k^-2 that gives NA its |A - C|^3 cancels between NA and Y_mu + zeta_A W_mu. Test: K6's instrument on the
  exact F_L at 2000 digits, L <= 240: beta_eff(220) in [-4.3, -3.85] and one sign for L >= 20, in both cases
  (NA alone: -2.03). Falsifier: beta_eff(220) above -3.85 (a surviving 3/2 or 5/2 branch) or below -4.3.
  Probability 0.75.

- **P4 (class; IC's instrument, primes 1073741789 and 1073741833, case A1 only).** Each mu-carrying component
  (constant of C not among NA's six) is run through IC's pipeline. **At least one mu-carrying component is clean
  STUCK, confirmed by both operator sets, on both primes:** the screened pieces leave the triangular (two-centre)
  class on their own, not only through NA. Falsifier: every mu-carrying component TRIANGULAR on both primes.
  INCOMPLETE (flags) leaves it UNDECIDED. Probability 0.7. Reported: each residual's IC reading and its lambda set.

**Reported, not scored:** R_min, Q0 and the trade table of the generic sequence; every exponent at x*; the
operators of Y_mu and W_mu alone; P3b's beta_eff for Y_mu alone (predicted -2, the Yukawa kernel is Coulombic at
short range); the whole-f class verdict (see below).

## Already known from the controls (risk of a void prediction, disclosed)

- **The reduction identity at the value level is checked to 1.5e-9 by K4** (R1 uses V_A, R2 does not). P1 scores
  the exact split, which is not built; the identity itself is textbook.
- **The whole-f class verdict is effectively known.** At A1 the Coulomb piece is N1, whose six components IC found
  clean STUCK on both primes. If P2d holds, the two E1 constants of NA appear in f with NA's rational coefficients,
  so under Qbar-independence f's minimal operator contains those components' operators and f is NOT TRIANGULAR.
  This is why P4 scores only the mu-carrying components.
- **P3a's exponents 3/2 and 5/2 are inherited the same way:** the minimal operator of the generic sequence annihilates
  every component, NA's included. P3a's live content is the absence of a new, mu-dependent branch point; P3b is the
  function-level statement and is not implied by any control.
- **P3b has an analogue:** record B2X found the slowest term (x* - x)^{7/2}, decay L^-4, for (XY|XZ). That is a
  different integral; the mechanism here (cancellation of the point-charge term) is derived above, not borrowed.
- The genus prediction of Paper 59's rule (genus zero, one two-centre density) is not scored: P4's NOT TRIANGULAR
  would not contradict genus zero (the walls register's K0 point), and no instrument here reads genus.

## Decision gate

- **VOID** if AB0 fails, or P1 fails (the engine is then wrong and nothing downstream is scored).
- **GO** if P1, P2, P3b hold: the screened reduction is exact, its constants form a fixed finite set, and the
  combination's leading branch at x* is the 7/2 the Fourier argument predicts.
- **STOP** if P3b is missed by a beta_eff that settles at -2 or -3 (the cancellation argument is wrong).
- **BORDERLINE** otherwise. P3a and P4 are scored and reported but do not move the gate.

## Novelty

Silverstone, PNAS 58, 34 (1967) (PMC335592, read in the prior-art memo) already has (1s_a 1s_a | 1s_b 1s_c) as a
Legendre series in the bond angle, each coefficient in exponentials, E_n, logarithms and modified spherical Bessel
functions, for R > R' only, with no closed sum and no statement about the angle function. So "each channel is in the
two-centre class" is his, and P1/P2's channel form is a check of his. New here, if they hold: the l-independent
finite constant set with rational coefficients (P2), the cancellation of NA's 3/2 and 5/2 branches and the 7/2
leading branch (P3b), and any operator and class statement in the angle (P3a, P4). Not searched for this record: the
Silverstone-Kay JCP series (1968-70) beyond its abstracts.

## Caps

1s functions; D1 < D2; two geometries; two primes per test. Function-level statements in the angle at fixed
rational lengths and exponents; no value-level class statement. P4 inherits IC's caps, including the erratum on
degree >= 2 factors with non-integer exponent classes (to be rechecked on any stuck residual). Qbar-independence of
the constants is assumed, as in Paper 59. P3b's decay window is L <= 240; it cannot exclude a 3/2 branch with a
coefficient below about 1e-2 of the 7/2 term at L = 220 (not fire-tested). NOT FOUND means beyond the budget.

## Phase 2 (after this record is committed)

Build `debug/aabc_exact.py` (na_exact plus the Yukawa kernel and its mu-derivative), then: P1 at both referee
points; P2 over L = 0..1000 mod both NA primes; P3a with na_modguess/na_modanalyze; P3b with the K6 instrument;
P4 with `intclass_score.py`'s pipeline on the mu-carrying components of A1.

## Results (2026-11-14, scored as registered; record text above unchanged)

Engine `debug/aabc_exact.py` (na_exact's Field, Bessel tables, J_n and finite-part rules reused; density factor copied;
kernels Y = mu(2L+1) i_L k_L and dY/dmu = (2L+1)[(2L+1) i_L k_L + x i_{L+1} k_L - y i_L k_{L+1}]; f = NA - Y + zeta_A dY).
Self-check: with NA's kernel w_L the engine reproduces `na_exact.coefficients` key-for-key at both cases, over Q and
mod 2^61 - 1. Drivers `debug/aabc_run.py` (p1p3b, perL, p2 (+ P3a search), p4) and `debug/aabc_p3a_points.py`;
data `debug/data/aabc_p1p3b_{2000,2600}.json`, `aabc_perL.json`, `aabc_p2_<p>.json`, `aabc_p3a_<case>_<p>.json`,
`aabc_p3a_points_<case>_<p>[_shift].json`, `aabc_p4_<p>.json` (+ logs), `aabc_ab0_k{5,7}_*.log`.

| Clause | measured | tolerance / criterion | score |
|:--|:--|:--|:--|
| AB0 | K5 rerun: (R, Q0) = (8, 31), nullity 1, order 7 none (Q = 106), both primes (260 s, 392 s). K7: ALL PASS on 1073741789 (44 s) and 1073741833 (45 s) | both reproduce | **HELD** |
| P1 | exact split, 2000 digits, L <= 240 (last term 6e-51 / 2e-61). A1: 0.212864498761243740935; minus stated R2 -3.9e-11 (minus R2 n = 26: +1.5e-12); minus R1 finest +1.56e-9. A2: 0.223136849487580864384; minus stated R2 -1.2e-11 (n = 26: -4.4e-12); minus R1 -1.01e-9. Values identical at 2600 digits | <= 1e-10 (R2), <= 3e-9 (R1), both cases | **HELD** |
| P1 (reported) | per-L: NA, Y, W separately vs independent quadrature (density multipole by Gauss-Legendre in s = \|r - C\|, mpmath Bessel kernels, W by numerical d/dmu), L = 0, 1, 4, 12, 30, both cases: max relative difference 2.1e-41 | target 1e-40 | met |
| P2a | L = 2..1000, both cases, both NA primes: the nonzero-coefficient set is the same 21 constants at every L, all inside C (C as listed = 29 keys; 8 unused, e.g. the E1 at D1 with prefactor e^{-zC D1 - mu D2} for the two rates shared with the middle region) | L-independent, inside C | **HELD** |
| P2b | gamma coefficient zero at every L, both cases, both primes | no gamma | **HELD** |
| P2c | four logs e^{-zC D1 - mu D2} ln\|b\|, \|b\| = 1/2, 3/2, 21/10, 41/10 (A1) [3/5, 7/5, 9/5, 19/5 (A2)] | >= 1 log | **HELD** |
| P2d | constants unchanged when zeta_A -> zeta_A + 1/997: exactly NA's six; Y and dY have zero coefficient on e^{+-zC D1} E1((zB + zC) D2) at every L <= 1000 | as stated | **HELD** |
| P2 | P2a, P2b, P2d | | **HELD** |
| P3a | generic sequence, both cases, both primes: R_min = 17, Q0 = 48, nullity 1 (trade: no operator for R <= 16). At x* (multiplicity 3): exponents 0..13, 3/2, 7/2, 9/2 -- **no 5/2**. Further finite singular points (linear and degree 2..18 factors) change position when zeta_A -> zeta_A + 1/997, and several carry non-integer exponents (a non-rational root at a moving linear point; non-constant exponent sums at moving higher-degree points), on both primes and both cases. Fixed under the shift: x*, +-1, and one rational point with integer exponents {-2, 0..15} (A1: -179401/179400; A2: 2509/1980) | x* sole non-integer point, no mu-dependent point, exponents include 3/2 and 5/2 | **MISSED** |
| P3b | K6 instrument on exact F_L, 2000 digits (identical at 2600): beta_eff(220) = **-4.106 (A1), -4.083 (A2)**, monotone toward -4 (A1: -4.82, -4.51, -4.28, -4.19, -4.14, -4.11 at L = 20, 40, 80, 120, 160, 220); F_L > 0 for all L >= 5, both cases | in [-4.3, -3.85], one sign for L >= 20 | **HELD** |
| P4 | A1, 15 mu-carrying components, IC primes 1073741789 and 1073741833, L = 1200: every component STUCK from both operator sets at the same residual order (order 4: the eight E1 and four ln components; order 5: e^{-162/25}; order 6: e^{-189/50}, e^{-513/50}), but none clean -- every residual carries the flag "lambda edge has 4 roots outside F_p" -> STUCK-UNCONFIRMED; none TRIANGULAR | >= 1 clean STUCK on both primes | **UNDECIDED** (INCOMPLETE) |

**Reported, not scored.** beta_eff(220) of the pieces: NA -2.029 / -2.021, Y alone -2.024 / -2.017 (predicted -2: held),
W alone -3.057 / -3.043 (A1 / A2). The pieces' 3/2 branch cancels in f; W's leading term is one order weaker than
Y's. Constant sets: Y and dY each carry 17 constants. Operator search for Y, W alone and the whole-f class verdict not
run. P4's flag is an instrument cap (the slope-one edge at infinity has rates not in F_p for these primes), not a
TRIANGULAR finding.

**Gate: GO** (AB0 held, P1 held, P2 held, P3b held; P3a missed and P4 undecided do not move the gate).

**New against Silverstone 1967**: the constants of every Legendre coefficient form one fixed set of 21 (NA's six plus
15 mu-carrying: 8 E1, 4 logs, 3 exponentials) with rational coefficients, L-independent from L = 2, with no gamma and
with logs absent from NA; the point-charge branches (x* - x)^{3/2} (and NA's 5/2) cancel and the coefficient decay is
L^{-4} (leading (x* - x)^{7/2}); the generic sequence has an order-17, degree-48 operator whose only half-integer
exponents at x* are 3/2, 7/2, 9/2 and which has further, mu-dependent singular points with non-integer exponents.

**Record defects (reported, not fixed).**
1. "Already known" / P3a's inheritance argument says NA's components appear in f unchanged. Only NA's two E1 components
   do (P2d); the two exponentials e^{+-zC D1 - (zB+zC) D2} also receive mu-free contributions from Y and W
   (measured: both Y and dY carry them; e.g. the outer-region endpoint e^{(s mu - mu) D2} with s = +1). The inheritance of
   NA's 5/2 exponent therefore did not follow from the controls; it is absent.
2. P3a's "no finite point whose position depends on mu" did not exclude points carried by single components that
   cancel in f itself (the operator annihilates each component). The clause is scored as written (on the operator).
3. P4: the record does not anticipate the IC instrument's lambda-outside-F_p flag on these components; scored
   UNDECIDED per "INCOMPLETE (flags) leaves it UNDECIDED".
4. The record's control-run outputs `debug/data/aabc_na_posctrl.json` and `intclass_controls_<p>.json` were
   rewritten by the AB0 reruns (same verdicts; timings differ).
