# Pre-registered prediction: (AA|BC) with p orbitals in the bond angle (AP; inventory row 8, part 2a)

**Registered 2026-11-14**, before any exact constant split, decay exponent or operator is computed for the two
integrals below.
**Commit status:** committed on its own, before any computation it scores.
**Rests on:** record AB (`2026-11-14_aa_bc_three_centre.md`, GO: the 1s integral (1s_A 1s_A | 1s_B 1s_C), its
reduction NA - Y_mu - zeta_A W_mu, the fixed set of 21 constants, the cancellation of NA's 3/2 branch and the
(x* - x)^{7/2} leading branch, beta_eff(220) = -4.11 / -4.08); record R8 (`2026-11-14_row8_two_centre_l2.md`, GO: in the
two-centre layer p and d orbitals add no function type); Paper 59 `sec:scope`; Paper 63 `sec:dictionary`.
**Guardrail check (Papers 8-9, FCI-M, Track DF):** not in scope; no molecular encoding is proposed.

## Objects

AB's geometry and exponents, unchanged: B at the origin, C at (0, 0, D1), A at D2 (sin g, 0, cos g), x = cos g =
cos(angle C B A), f(x) = sum_L F_L P_L(x). Normalized Slater functions.

    PC = (1s_A 1s_A | 1s_B 2pz_C),   2pz_C = sqrt(zC^5/pi) (z - D1) e^{-zC rC}       (p along the B-C bond, away from B)
    PA = (2ps_A 2ps_A | 1s_B 1s_C),  2ps_A = sqrt(zA^5/pi) (s . n) e^{-zA s},  s = r - A, n = A/|A|  (p along B-A)

**PC** keeps AB's potential V_A and changes only the density. Its Legendre multipoles about B follow from the 1s
ones through z - D1 = r cos(theta) - D1 and cos(theta) P_l = ((l+1) P_{l+1} + l P_{l-1})/(2l+1).

**PA, exact reduction (derived and checked symbolically before registering, control K0).** The potential of
(2ps_A)^2 is V0(s) + V2(s) P2(s^.n), with V0 = 1/s plus exponentials and V2 = 3/(zA^2 s^3) plus exponentials with
powers down to s^-3. Moving A along n is d/dD2 at fixed x, and (n.grad_A)^2 U(s) = (U'' - U'/s)(s^.n)^2 + U'/s. So

    PA = d^2/dD2^2  int rho_BC U  +  int rho_BC W0,   mu = 2 zA, e = e^{-mu s},
    U  = 3/(2 zA^2 s) - 3 e/(2 zA^2 s) - 33 e/(16 zA) - 9 s e/8 - zA s^2 e/4,
    W0 = 1/s - e/s - 5 zA e/4 - zA^2 s e/2.

Every radial kernel above is AB's family: 1/s (NA), e^{-mu s}/s (Y) and its mu-derivatives up to order 3. Nothing
outside the family is needed.

## Cases

AB's A1 (D1 = 9/5, D2 = 13/5, zB = 13/10, zC = 1, zA = 9/10, x* = 125/117, rho = 9/13) and A2 (D1 = 3/2, D2 = 12/5,
zB = 1, zC = 6/5, zA = 4/5, x* = 89/80, rho = 5/8). Both integrals at both cases. Referee points x = 1/2 (A1),
x = -1/3 (A2). Every exponential rate is AB's, so no rate vanishes (AB's control K8).

## Controls measured before registering

Driver `debug/aapc_controls.py`; data `debug/data/aapc_controls.json`, `debug/data/aapc_referee.log`.

| Control | what | measured |
|:--|:--|:--|
| K0 | PA reduction: V0, V2 by sympy from the density; (n.grad)^2 U + W0 = V0 + V2 P2 identically; total charge 1 | holds; charge 1 |
| K1 fit floors (geovac closed forms; hydrogenic 2p0 = Slater 2p at zeta = Z/2) | R2's p construction on (1s_A 1s_A \| 1s_B 2p0_B) and (2p0_A 2p0_A \| 1s_B 1s_B), R = 2 | n = 18 / 22 / 26: PC-type 2.2e-9 / 3.1e-10 / 4.0e-11; PA-type 1.0e-12 / 1.1e-11 / -1.1e-12 |
| K2 referee R2 (Gaussian fits; p functions as centre derivatives of s Gaussians; closed-form ERIs and their analytic derivatives; uses no reduction) | PC at A1; A2 | n = 22: -0.2089629965564, -0.2183979664955; n = 26: -0.2089629964862, -0.2183979664315 |
| | PA at A1; A2 | n = 22: 0.2032312569203, 0.2195256149296; n = 26: 0.2032312569631, 0.2195256150153 |
| K3 referee R1 (prolate quadrature of rho_BC against the exact V_A; PA uses V0 + V2 P2) | PC, grids 60^2x72 / 90^2x108 / 120^2x144 | A1: -0.2089629891, -0.2089629949, -0.2089629958; A2: -0.2183979609, -0.2183979731, -0.2183979693 |
| | PA, same grids | A1: 0.2032312576, 0.2032312570, 0.2032312570; A2: 0.2195256151, 0.2195256149, 0.2195256150 |
| K4 R1 vs R2 (n = 26), finest grid | PC A1; A2. PA A1; A2 | 7.1e-10; 2.8e-9. 3.2e-11; 3.4e-11 |
| K5 decay instrument (inherited, not re-run) | AB's K6 on NA: beta_eff -> -2; AB's f: -4.11 / -4.08; W alone -3.06 / -3.04 | separates the 3/2, 5/2 and 7/2 branches at L = 220 |

R2's spread between n = 22 and 26 is 7.0e-11 / 6.4e-11 (PC) and 4.3e-11 / 8.6e-11 (PA), and PC's fit floor is
4e-11. The referee values are **PC A1: -0.20896299649 (+-1.5e-10); PC A2: -0.21839796643 (+-1.5e-10);
PA A1: 0.20323125696 (+-1e-10); PA A2: 0.21952561502 (+-1e-10).** R1 for PC is noisy (non-monotone on A2) and is a
loose check only.

## Predictions (derived before computing; none uses an exact split or decay exponent of PC or PA)

- **AP0 (engine self-checks; run first).** (i) With the factor (z - D1) removed and A's function set back to 1s,
  the extended engine reproduces AB's coefficients key for key at both cases, over Q and mod 2^61 - 1, L <= 30.
  (ii) The mu-jet of the Y kernel reproduces `aabc_exact`'s dY kernel exactly. (iii) The D2-jet (second derivative,
  including the jump term at r = D2) matches a 2000-digit central difference of the exact F_L (h = 10^-300) to
  10^-500 at L = 0, 1, 5, 20, both cases. Falsifier: any of (i)-(iii) fails. Probability 0.9.

- **P1 (exact values).** The exact split, evaluated at 2000 digits and summed to L = 240 at the referee point,
  matches R2 to **<= 3e-10 (PC) and <= 2e-10 (PA)** in both cases, and R1 to <= 5e-9 (PC) and <= 3e-10 (PA).
  Falsifier: a larger difference anywhere. Probability 0.8.

- **P2 (finite, fixed constant set).** Mod both NA primes (2^61 - 1, 4611686018427387847), L <= 1000:
  (P2a) PC's set of constants with nonzero coefficient in F_L is the same at every L >= 3, and PA's at every L >= 2;
  (P2b) both sets lie inside AB's 29-key candidate set C at the same case, with no new E1 or log argument;
  (P2c) gamma does not appear. P2 holds if P2a, P2b and P2c hold for both integrals. Probabilities: P2a 0.8,
  P2b 0.75, P2c 0.85.
  Reported, not scored: whether PC's set is exactly AB's 21 (guess 0.55: the factor r cos(theta) raises powers of r,
  which can remove an E1 or log, and the D1 term keeps them), and which of AB's 8 unused C keys PA uses (the
  D2-derivative of e^{-zC D1 - mu D2} E1((zB + zC + mu) D2) can bring in e^{-857/50} at A1, an unused key of C).

- **P3 (leading branch at x*; the new content).** Both integrals keep AB's leading branch (x* - x)^{7/2}: in
  AB's instrument (exact F_L, 2000 digits, L <= 240), **beta_eff(220) lies in [-4.3, -3.85] for PC and for PA, in
  both cases**, with one sign of F_L for L >= 20.
  Reasons, derived from Fourier transforms at large k:
  - **PC.** The node of 2pz at C does not remove the cusp. The non-analytic part of the density at C is
    -zC (z - D1) r_C = -zC r_C^2 P1, whose transform is ~ k^-5 P1. Against V_A ~ k^-6 that gives
    k^-11 P1 ~ |A - C|^8 P1 = |A - C|^7 (A - C)_z, and (A - C)_z = D2 x* - D1 is nonzero at x*. The wrong
    alternative, "the node kills the cusp", predicts 9/2, i.e. beta -> -5.
  - **PA.** The density's monopole part is like s^3 at s = 0 (k^-6, so V ~ k^-8). Alone it would give k^-12 against
    the 1s cusp of rho_BC, i.e. |A - C|^9 and a 9/2 branch (beta -> -5; this was the planning estimate, which this
    record replaces). But the quadrupole part has the same k^-8 order with angular factor P2, and
    k^-12 P2 ~ |A - C|^9 P2 = |A - C|^7 h2, where h2 = (3 (n.(A - C))^2 - |A - C|^2)/2 is a harmonic quadratic.
    At x*, |A - C|^2 = 0 and n.(A - C) = (D2^2 - D1^2)/(2 D2) is not 0, so the 7/2 branch survives with a nonzero
    coefficient.
  Falsifier: beta_eff(220) above -3.85 or below -4.3 for either integral in either case. Probability: PC 0.75,
  PA 0.65.

**Reported, not scored:** beta_eff of the pieces (PA: d^2/dD2^2 of the U integral and the W0 integral separately;
the 1/s^3 point-quadrupole term is predicted to give the strongest branch, cancelled in the sum); constant counts per
piece; per-L agreement of each piece with an independent radial quadrature (target 1e-40).

## Already known from the controls (risk of a void prediction, disclosed)

- K4 checks the values to 3.4e-11 (PA) and 3e-9 (PC). Neither referee uses the PA reduction: R1 integrates
  against V0 + V2 P2 directly, and R2 uses no potential at all. The reduction is checked only symbolically (K0).
  P1 scores the exact split, which is not built.
- P2b is partly structural: d/dD2 cannot create a new E1 or log argument, only polynomial prefactors and
  exponentials. The live content is that the exponentials it creates stay inside C, and that the set is fixed in L.
- P3 has no control: no decay exponent of PC, PA or their pieces has been computed.

## Decision gate

- **VOID** if AP0 fails or P1 fails (the engine is wrong; nothing downstream is scored).
- **GO** if P1, P2 and P3 hold: p orbitals at either end of the three-centre integral keep the 1s constant set and
  the 7/2 branch, and row 6's 1s structure extends to row 8.
- **STOP** if P3 settles at beta near -5 for PA (the quadrupole argument is wrong) or near -2 or -3 for either
  integral (the point-charge cancellation fails at l > 0).
- **BORDERLINE** otherwise.

## Novelty

Three-centre integrals over Slater p functions are standard numerical objects (Fourier-transform and B-function
methods). We know of no statement about their bond-angle function: fixed constant set, branch at x*, or decay of the
Legendre coefficients. Not searched beyond AB's prior-art check (Silverstone 1967 PNAS: 1s only).

## Caps

Two integrals (one p function at C; one p_sigma pair at A), Slater exponents and geometries of AB, D1 < D2, two
primes. p_pi functions (2px on A or C) and d functions are not covered. They need angular structure off the B-C
axis: PA's 2pz^2 needs the X/T recurrences, and 2px^2 a P_L^2 basis. Exchange at l > 0 and Paper 63's fibre at
l > 0 are row 8, part 2b. Qbar-independence of the constants is assumed, as in Paper 59. P3 cannot exclude a
competing branch whose coefficient is below about 1e-2 of the 7/2 term at L = 220.

## Phase 2 (after this record is committed)

Extend `debug/aabc_exact.py` (na_exact's Field, Bessel tables, Integrals, value unchanged):
- `density_pz`: PC's multipoles from AB's density factor by the cos(theta) recurrence (L-1 and L+1) and the -D1 term.
- `kernel_jet`: mu-derivatives of the Y kernel to order 3 (Bessel recurrences), and D2-derivatives to order 2 of
  every kernel in each region.
- `boundary_D2`: the jump term of d^2/dD2^2 at r = D2, rho_L(D2) [d_D2 K_near - d_D2 K_far] at r = D2.
Then: AP0; P1 at both referee points (2000 and 2600 digits); P2 over L <= 1000 mod both primes; P3 with AB's
instrument on exact F_L at 2000 digits, L <= 240.

## Results (2026-11-14; appended after scoring, text above unchanged)

**Verdict: GO.** AP0, P1, P2 and P3 held in both cases for both integrals.

| clause | A1 PC | A2 PC | A1 PA | A2 PA |
|:--|:--|:--|:--|:--|
| P1: exact value (2000 digits, L <= 240) | -0.208962996502063 | -0.218397966426964 | 0.203231256960264 | 0.219525615007203 |
| P1: minus R2 stated (tol 3e-10 / 2e-10) | -1.2e-11 | 3.0e-12 | 2.6e-13 | -1.3e-11 |
| P1: minus R1 finest (tol 5e-9 / 3e-10) | -7.2e-10 | 2.8e-9 | -3.5e-11 | 2.6e-11 |
| P2a: set fixed over L <= 1000 | yes, 21 keys | yes, 21 | yes, 21 | yes, 21 |
| P2b: inside AB's 29-key C | yes | yes | yes | yes |
| P2c: no gamma | yes | yes | yes | yes |
| P3: beta_eff(220), band [-4.3, -3.85] | -4.123 | -4.101 | -4.086 | -4.059 |
| P3: one sign for L >= 20 | yes (+) | yes (+) | yes (+) | yes (+) |

- **AP0** held (`debug/data/aapc_ap0.log`): AB reproduced key for key over Q and mod 2^61 - 1; mu-jet equals AB's dY
  kernel; D2-jet against the 2000-digit central difference to 4.7e-604 ... 6.8e-613.
- **P1.** Tail terms at L = 240 are below 3e-51. Both primes give identical constant sets (P2); the 2600-digit run
  reproduces the 2000-digit values to all 30 printed digits and beta_eff to six.
- **P2.** At L = 0 only the 7 pure exponentials appear; from L = 1 on, the set is the full 21 in all four cases,
  earlier than the predicted L >= 3 (PC) and L >= 2 (PA). Identical mod both primes.
- **P3.** beta_eff rises toward -4 from below at every L sampled (PC A1: -5.14 at L = 20, -4.33 at 80, -4.12 at
  220), as for AB. F_L > 0 for L >= 5 in all four.

**Reported, not scored.**
- PC's set is exactly AB's 21 at both cases (guess 0.55). The factor r cos(theta) removes no E1 or log.
- PA's set is also exactly AB's 21. It uses none of C's 8 unused keys: the e^{-857/50} candidate does not appear.
- PA pieces: the d^2/dD2^2 U piece has beta_eff(220) = -4.046 (A1), -4.025 (A2), so it carries the 7/2 branch;
  the W0 piece has -5.166 (A1), -5.131 (A2), drifting toward -5, the 9/2 branch. This fits the P3 reasoning: the 7/2
  branch comes from the quadrupole part, which lives in (n.grad)^2 U, while W0 alone has only the monopole-type
  9/2 branch. The pieces were split as d^2/dD2^2 U versus W0, not with the 1/s^3 term separately, so the predicted
  cancellation of that term inside the sum was not isolated.
- Per-piece constant counts and the per-L 1e-40 quadrature check of each piece were not run.

**Disclosed changes to the Phase-2 plan.**
- Engine speed. The first runs mod the two primes ran 6 h with no case finished (a timing test showed L^3 cost,
  ~85 h per PA case at L = 1000). The cause was keeping each Y kernel's double Bessel sum symbolic in (mu, D2, r).
  A factored path (`kernel_eval`: D2 factor and r factor kept apart, mu-derivatives by Leibniz) replaced it before
  any P2 output existed. It gives identical coefficients to the AP0-checked path for PC, PA and AB, both cases, over
  Q and mod both primes, L <= 30; AP0 passes on it (`debug/data/aapc_ap0_fast.log`); and the 2600-digit P1/P3 run
  (fast path) matches the 2000-digit run (old path) at L <= 240. P2 ran on the fast path (~440 s PC, ~870 s PA per
  case per prime).
- Files: `debug/aapc_exact.py`, `debug/aapc_run.py`; data `debug/data/aapc_p1p3_{2000,2600}.json`,
  `debug/data/aapc_p2_{2305843009213693951,4611686018427387847}.json`.

**Consequence for inventory row 8.** Part 2a is closed: a p function at either end of (AA|BC) keeps the 1s constant
set (21 keys, no gamma, no new E1 or log) and AB's 7/2 branch at x*. Part 2b (p_pi and d functions off the B-C
axis, exchange at l > 0, Paper 63's fibre at l > 0) remains.

**Found while writing the guard (2026-11-14).** The jump term of d^2/dD2^2 at r = D2 is identically zero (every L
checked, both cases, over Q and mod p). U has no 1/s singularity (its NA and Y parts cancel at s = 0), so the
first D2-derivative of its kernel is continuous at r = D2. AP0 (iii), described as "including the jump term",
therefore checked the volume part alone; the result is unaffected. Guard: `tests/test_paper59_aa_bc_p.py`
(fire-tested 17/17).
