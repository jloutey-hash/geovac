# Record AS: (AA|BC) with p_pi and d pairs on A (inventory row 8, part 2b: the A side)

*Registered 2026-11-15, before computing, under the standing autonomous-loop authorization (PI, 2026-11-14; the PI
chose this record as the close of row 8 part 2b on 2026-11-15). Engine and controls were built by a sub-agent that
was barred from the scored quantities; the PM read its report, added the amplitude toys and the total re-expansion,
and wrote the predictions. A fresh reviewer read the draft cold and re-derived the mechanism (holds, including at
every B-gradient order); its findings were applied: P4's band widened and re-centred (MUST-FIX: 9/2 terms present
in a and absent in b put the correct ratio about 7% above x*^2 - 1), PX' scored on its total only (b' = -b
identically), the re-expansion check gated, the frozen-file list completed, P0's c restricted to L = 219..221.*

## Object

As in records AB and AP: B at the origin, C = (0, 0, D1), A = D2 (sin g, 0, cos g), x = cos g, D1 < D2; the BC pair
is 1s_B 1s_C. New A pairs:

- **PY = (2p_y,A 2p_y,A | 1s_B 1s_C)**, p_y perpendicular to the plane of the three centres. As a function of A it
  has an m = 0 and an m = 2 part about the B-C axis: PY(x) = sum_L a_L P_L(x) + sum_L b_L P_L^2(x), with
  P_L^2 = (1 - x^2) P_L'' (no Condon-Shortley phase) and (a, b) the 3-D harmonic split (a choice: a function of x
  alone does not fix it). The **total** as a Legendre series, c_j = a_j - j(j-1) b_j + 2(2j+1) sum_{L > j, L = j mod 2}
  b_L (exact identity, checked to 1e-62 at L <= 12).
- **PX' = (2p_x',A 2p_x',A | 1s_B 1s_C)**, p_pi in the plane, perpendicular to B-A. (a', b') is the engine's split,
  with b' = -b identically and x' moving with A, so it is not a 3-D harmonic split; only the total c' is an
  independent PX' prediction.
- **DA = (3d-sigma_A 3d-sigma_A | 1s_B 1s_C)**, sigma axis along B-A; a plain Legendre series F_L.

Engine `debug/aaa_exact.py` (imports `aapc_exact.py`, unchanged): the A potential is written as
sum_j (e.grad)^{2j} U_{2j}(|w|) with radial U's in the engine's kernel family; DA = G4'''' + G2'' + W (D2-derivatives
at fixed x; every jump term at r = D2 is identically zero up to fourth order, exact at L <= 8); PY, PX' from
d^2/dA_y^2 = (lap - d_z^2)/2 - (d_+^2 + d_-^2)/4 acting on G_L(D2) P_L, finite in L' in {L-2, L, L+2}.

## Controls (measured; `debug/aaa_controls.py`, `debug/aaa_ap0.py`, `debug/data/aaa_controls.json`)

| id | what | result |
|:--|:--|:--|
| AS0 | self-checks at L <= 8: exact channels vs direct quadrature of the multipole potentials (no U kernels); PY + PX' + PA = 3 S exactly (key for key, L <= 30); potentials vs Gauss's law; jump terms and the fourth D2-derivative vs a 5-point difference | worst relative 6e-59 (PY a, b), 4e-58 (PA), 1.4e-50 (DA); identity exact; potentials 1e-31; 4th derivative 1e-80 |
| R2 | Gaussian fits n = 30 / 34 (AB's fits; p and d pairs on A by moment identities) at A1 x = 1/2, -1/3 and A2 x = 0, 4/5 | values in the json; n = 30 -> 34 spread <= 1.7e-12; fit floors at n = 34 against exact one-axis cases 8e-15 (p_y), 1.2e-13 (p_z), 2.6e-13 (d) |
| R1 | prolate quadrature 120^2 x 144 | within 8.5e-12 of R2 (p pairs), 4.3e-13 (d) |
| T | toy beta_eff(220) on c_L / rho^L, rho = 1/(x* + sqrt(x*^2 - 1)), x* = (D1^2 + D2^2)/(2 D1 D2): P_L coefficients of (x* - x)^alpha, and P_L^2 coefficients of (1 - x^2)(x* - x)^alpha | A1: P_L 5/2 -3.070, 7/2 -4.119, 9/2 -5.180; P_L^2 5/2 -5.180, 7/2 -6.255, 9/2 -7.342. A2: P_L -3.055, -4.094, -5.143; P_L^2 -5.143, -6.203, -7.274 |
| TA | amplitude toys at L = 218..222: ta (P_L of (x* - x)^{7/2}), tb (P_L^2 of (1 - x^2)(x* - x)^{7/2}), tc (P_L of (x* - x)^{9/2}) | in the json (`amp_toys`); x*^2 - 1 = 0.141433 (A1), 0.237656 (A2) |
| X | small-s expansions of the A potentials (Gauss's formula, sympy): lowest odd power of V_l / s^l | p pair: l = 0 -> s^5, l = 2 -> s^3. d pair: l = 0 -> s^7, l = 2 -> s^5, l = 4 -> s^3 (1s pair, for comparison: s^3) |

**Already known (disclosed):** a driver smoke test at L <= 8 showed PY's constant set equal to AB's 21. Nothing
else about the new coefficients beyond L = 12 has been computed or seen.

## Mechanism (the reason the predictions differ between PY, PX' and DA)

The non-analytic part at x* comes from the cusp |u| of the BC density at C (u = r - C) against the odd powers of
the A potential: a term h_l(w) s^n (n odd) gives |a|^{n+4} h_l(a), a = A - C (the 1s pair: n = 3, |a|^7, record
AB's 7/2 branch). At x*, |a|^2 = 0 (x* > 1, a complex null vector). From X:

- **p pairs:** the monopole reaches only |a|^9; the quadrupole gives |a|^7 h_2,e(a), e the p axis. For **PY**,
  e = e_y and a_y = 0 for every x, so h_2,y(a) = -|a|^2/2 and **the total PY leads with |a|^9, (x* - x)^{9/2}**.
  Its m = 0 and m = 2 parts do not vanish separately: at the null point they are -(3/4) a_z^2 |a|^7 and
  -(3/4) a_x^2 |a|^7 with a_x^2 = D2^2 (1 - x^2), equal and opposite, so **a and b each lead with 7/2 and cancel in
  the total**, with amplitude ratio A_a / A_b = a_z(x*)^2 / D2^2 = **x*^2 - 1** (using a_x^2 + a_z^2 = 0 at x*).
  The cancellation in the PY total holds at every B-gradient order: on the null cone a harmonic built from
  e_y e_y reduces to terms in a_y^2 = 0 or |a|^2. For **PX'**, e . a = D1 sin g does not vanish: the total keeps 7/2.
- **d pair (DA):** l = 4 gives |a|^7 h_4,n(a) with h_4 = (35/8)(n . a)^4 on the null cone and n . a = (D2^2 - D1^2)/(2 D2)
  at x*, nonzero: **DA leads with 7/2**. (l = 2 reaches |a|^9, l = 0 |a|^11.)

## Predictions (scored at A1 and A2)

- **AS0** passed (above). Engine files `debug/aaa_exact.py`, `aapc_exact.py`, `aabc_exact.py`, `na_exact.py`,
  `aabc_run.py` and the driver `debug/aaa_run.py` unchanged from this commit at scoring.
- **P0 (precision gate).** The scored run uses L <= 400 at 3000 and 3600 digits; every series' coefficients at
  L = 219..221 and 398..400 (a, b, F; c at 219..221 only, since c near L = 400 lacks the b tail) and all summed
  values agree between them to relative 1e-10, and `total_reexpansion_check` is below 1e-30 at all four points.
  Failure is VOID.
- **P1 (values).** The series summed to L = 400 match R2 (n = 34) within **1e-10** at both points of both cases, for
  PY, PX' and DA (12 values; summing sum a P_L + sum b P_L^2 for the p pairs). Probability 0.85.
- **P2 (constants).** Mod 2^61 - 1 and 4611686018427387847, up to L = 1000: a_L and F_L (DA) from L = 2, b_L from
  L = 3 (b_2 is the first nonzero), each use exactly AB's 21 constants at the same case, the set fixed in L; no
  Euler gamma. Probability 0.75 (AP,
  AQ and AD all found the 21; the d kernels need d^kY up to k = 5, which AP's did not).
- **P3 (exponents)**, beta_eff at L = 220 on coefficient / rho^L:
  - PY: a in [-4.65, -3.6] (7/2), b in [-6.8, -5.7] (7/2 in the P_L^2 basis), **total c in [-5.7, -4.65]** (9/2).
  - PX': **total c' in [-4.65, -3.6]** (7/2); a', b' reported, not scored.
  - DA: **F in [-4.65, -3.6]** (7/2).
  Probability 0.7 (the PY total: 0.6).
- **P4 (the cancellation, amplitude form).** R = (a_220 / ta_220) / (b_220 / tb_220) lies in **[0.95, 1.20] x
  (x*^2 - 1)** (x*^2 - 1 = 0.1414 at A1, 0.2377 at A2), positive, for PY. The leading ratio is exactly x*^2 - 1;
  at L = 220 the 9/2 terms present in a and absent in b (the trace part of h_2,z, the variation of a_z^2, the
  monopole's s^5 term) raise it, by the reviewer's estimate (not re-derived by the PM) to about 1.08 (A1) and 1.07
  (A2) times x*^2 - 1, falling like 1/L, with B-gradient terms about +-3%. Reported, not scored: R at 219, 221;
  c_220 / tc_220 for PY. Probability 0.6.

**Reported, not scored:** beta_eff at L = 20, ..., 380 for a, b, F (c to L = 220); signs; timings.

## Decision gate

- **VOID** if AS0, P0 or P1 fails, or an engine file changed.
- **GO** if P1, P2, P3 and P4 hold.
- **STOP** if PY's total c has beta_eff(220) above -4.65 at either case (no cancellation: the mechanism is wrong), or
  PY's a or b, PX' or DA is outside [-5.7, -3.0] (a leading branch other than 7/2 or 9/2), or the P4 ratio is
  negative or outside [0.5, 2] times x*^2 - 1.
- **BORDERLINE** otherwise (including P2 failing with P1, P3 and P4 holding).

## Caps

Two geometries, two primes; one p_pi orientation each (PY out of plane, PX' in plane), one d function (m = 0 along
B-A) on A, against a 1s pair at B and C. The (a, b) split is the 3-D harmonic split; P4 is stated in it. Constants
assumed linearly independent over Qbar. d functions with m != 0 on A and exchange integrals are not covered.

## Phase 2 (after this record is committed)

`python debug/aaa_run.py p1p3 400 3000` and `... 400 3600` (all three integrals; writes `debug/data/aaa_p1p3_<dps>.json`
with values, diffs, beta for a, b, c, F, the re-expansion check, P4 ratios); `python debug/aaa_run.py p2 2305843009213693951`
and `... p2 4611686018427387847` (constant sets to L = 1000).

## Results (scored 2026-11-15)

**Verdict: GO.** Engine files and driver unchanged at scoring (`git diff f675ba7f` on all six frozen files empty).
All clauses at 3600 digits unless stated; data `debug/data/aaa_p1p3_{3000,3600}.json`,
`debug/data/aaa_p2_{2305843009213693951,4611686018427387847}.json`.

- **P0 held.** The coefficients a, b, c, F at L = 219..221 and 398..400 of all six integrals and all twelve P1 sums
  agree between 3000 and 3600 digits in every printed digit (25 significant for coefficients, 30 for the sums).
  `total_reexpansion_check` is at most 4e-3001 at every point where it applies (PY and PX', both points, both
  cases; DA has no split).
- **P1 held (12/12).** Differences from R2 (n = 34): 4.2e-14 to 2.5e-13 (tolerance 1e-10).
- **P2 held.** At both primes, to L = 1000: a and F from L = 2, b from L = 3, and a + b, each use exactly AB's 21
  constants, the set fixed in L; no Euler gamma; nothing outside the set.
- **P3 held (10/10).**

  | beta_eff(220) | A1 | A2 | band |
  |:--|:--|:--|:--|
  | PY a | -4.209 | -4.167 | [-4.65, -3.6] |
  | PY b | -6.241 | -6.192 | [-6.8, -5.7] |
  | PY total c | -5.166 | -5.130 | [-5.7, -4.65] |
  | PX' total c' | -4.206 | -4.162 | [-4.65, -3.6] |
  | DA F | -4.151 | -4.107 | [-4.65, -3.6] |

- **P4 held (2/2).** R = (a/ta)/(b/tb) at L = 220: A1 0.15672 = 1.108 (x*^2 - 1); A2 0.25844 = 1.087 (x*^2 - 1)
  (band [0.95, 1.20]; the reviewer's estimate was about 1.08 and 1.07). R falls with L across 219..221 at both.

**Reported, not scored.** beta_eff at L = 380: PY a -4.121 / -4.097, b -6.141 / -6.112, DA F -4.087 / -4.062 (A1 / A2),
all still drifting toward the 7/2 branch from below. c_220 / tc_220 for PY: -0.0817 (A1), -0.0320 (A2). Signs at
L >= 20: PY a, b negative and c positive; PX' the reverse; DA F positive (one sign change below L = 20). Engine time
(3600 digits): 3 s to 4.2 h per integral.

**Reading.** A p_pi or d_sigma pair on A against a 1s pair at B and C keeps AB's 21 constants and the 7/2 branch, as
the p and d functions at C did (AP, AQ, AD). The out-of-plane p_pi pair (PY) is the one case where the total loses
a half-power: its harmonic parts each carry 7/2, and they cancel at the leading order to leave 9/2, as predicted
from e . a = 0 on the null cone; the in-plane PX' keeps 7/2. The amplitude ratio sits 9-11% above x*^2 - 1, at the
reviewer's estimate of the 9/2 admixture.
