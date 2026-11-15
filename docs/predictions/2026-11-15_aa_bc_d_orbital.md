# Record AD: (AA|BC) with a d function at C (inventory row 8, part 2b-ii)

*Registered 2026-11-15, before computing, under the standing autonomous-loop authorization (PI, 2026-11-14).
A fresh reviewer read the draft cold before commit; its findings were applied: P2a's status made consistent,
P3's limits disclosed (it cannot tell a d factor from none; P1 at four points does), an L_max-invariance clause
added to AD0, the fit floor's coverage disclosed, per-case P3 bands, sign failures gated, the L >= 3 reason
corrected.*

## Object

Slater functions, AB's two geometries: B at the origin, C at (0, 0, D1), A at D2 (sin g, 0, cos g), x = cos g.

    PD = (1s_A 1s_A | 1s_B 3ds_C),   3ds_C = N_d h2(r - C) e^{-zC rC},   h2(v) = (3 v_z^2 - |v|^2)/2,
    N_d = sqrt(5 (2 zC)^7 / (4 pi 720))   (normalized Slater 3d, m = 0 along B-C).

Cases A1 (D1 = 9/5, D2 = 13/5, zB = 13/10, zC = 1, zA = 9/10) and A2 (3/2, 12/5, 1, 6/5, 4/5).

**Reduction (derived, not computed).** About B, with z = r t, h2(r - C) = (3 r^2 t^2 - r^2 - 4 D1 r t + 2 D1^2)/2.
So the density multipoles follow from AB's 1s.1s ones by applying the multiply-by-t operator T (PC's recurrence,
(T d)_L = [L d_{L-1} + (L+1) d_{L+1}]/(2L+1) in AB's normalization) once and twice, with powers of r:

    rho^d_L = (1/2) [ 3 r^2 (T^2 d)_L - r^2 d_L - 4 D1 r (T d)_L + 2 D1^2 d_L ],

against AB's kernel NA - Y + zA dY (the 1s potential of A), unchanged. PD = sum_L F_L P_L(x), normalization
2 pi N_B N_d.

## Controls (measured; `debug/aadz_controls.py`, `debug/data/aadz_controls.json`)

| id | what | result |
|:--|:--|:--|
| F | fit floor: R2 against geovac's exact (1s_A 1s_A\|1s_B 3d0_B), Z_A = 1, Z_B = 3, R = 2 | n = 22 / 26 / 30 / 34: -2.7e-9, -4.2e-10, -6.9e-11, -1.1e-11 |
| R2 | Gaussian fits (AB's even-tempered fits of e^{-r}; the d factor as P-derivatives of the s-type ERI, up to second order), n = 22 / 26 / 30 / 34 | A1 x = 1/2: 0.10475756004825808, 0.10475756005424663, 0.10475756007757979, 0.10475756007475293. A1 x = 0: 0.11923238811191147, 0.11923238808250199, 0.11923238810611149, 0.1192323881040522. A2 x = -1/3: 0.11423709524035043, 0.11423709523509168, 0.11423709523315846, 0.11423709523279987. A2 x = 4/5: 0.07383525671510671, 0.07383525669170064, 0.07383525670173083, 0.07383525670081686 |
| R1 | prolate quadrature against the exact V_A, grids 60^2x72 / 90^2x108 / 120^2x144 | A1 1/2: 0.104757573558, 0.104757562026, 0.104757561942; A1 0: 0.119232385313, 0.119232388374, 0.119232388190; A2 -1/3: 0.114237091332, 0.114237099732, 0.114237097146; A2 4/5: 0.073835260172, 0.073835254172, 0.073835257669 |
| T | toy positive control for P3's instrument (record AQ, `debug/data/aaqx_controls.json`): (x* - x)^{7/2}, plain Legendre series, beta_eff at L = 40 / 160 / 220 | A1 -4.616, -4.162, -4.119; A2 -4.500, -4.129, -4.094 |

The d factor is harder on the fit than 2pz (floor 4.2e-10 at n = 26, against 4e-11 for 2pz); n = 34 brings it to
1.1e-11 and the n = 30 -> 34 spread at the four points is <= 2.8e-12. **Referee values (n = 34): A1 x = 1/2
0.10475756007, A1 x = 0 0.11923238810, A2 x = -1/3 0.11423709523, A2 x = 4/5 0.07383525670; tolerance 1e-10.**
R1 sits 0.2e-9 to 1.9e-9 from R2 at its finest grid, is non-monotone on A2, and is reported only (loose 1e-8).
Coverage of F, disclosed: at F's geometry B = C and A is on the axis, so the delta-linear term, the constant term
and the f1 (1 - u_z^2)/R part of gauss_PD never run there; those terms are confirmed only by R1, to about 2e-9. A
bug in them would be O(1), not subtle. F's exponents and geometry also differ from A1 and A2, and R2 is not
monotone in n at A1 (n = 30 vs 34); the n = 22 -> 34 spread is at most 2.6e-11, so 1e-10 keeps a margin of about 4.

## Predictions (none uses an exact split, a constant set or a decay exponent of PD)

- **AD0 (engine self-checks, run first).** (i) AB reproduced key for key through the unchanged path (Q and mod
  2^61 - 1, L <= 30). (ii) F_L at L = 0, 2, 5 (both cases) equals a quadrature in which the density multipole is
  projected directly (h2(r - C) e^{-zC s} against P_L(t), no T recurrence) against V_A's multipole, relative 1e-18
  at 22 digits, and also at L = 12. (A direct quadrature at L near 220 is not feasible: the multipole is about
  1e-36 of its size and oscillates at degree 220.) (iii) F_L for every L <= 220 is unchanged (exactly, mod
  2^61 - 1) when the internal L_max is raised from 240 to 260: catches a truncated top-of-range recurrence
  (T^2 reaches d_{L+2}). (ii) checks the engine against the reduction formula; P1 checks the reduction.
  Probability 0.9.
- **P1 (exact value).** The exact split at 2000 digits summed to L = 240 matches R2 (n = 34) to <= 1e-10 at both
  points of both cases. Probability 0.85.
- **P2a (constant set fixed in L; scored, near-structural).** Mod both primes, the set of constants with nonzero
  coefficient in F_L is the same at every L >= 2, L <= 1000. Probability 0.85. Near-structural: the density is
  built from AB's d_k with r-powers and the T recurrence, as PC's and PX's were. L = 0, 1 are excluded as in AB
  (AB's set was fixed from L = 2); their sets are reported.
  (P2b, P2c; reported, not scored, because AP and AQ imply them) inside AB's 29-key candidate set; no gamma;
  equal to AB's 21.
- **P3 (decay).** In AB's instrument, beta_eff(220) of F_L lies in **T +- 0.17 per case: A1 [-4.29, -3.95],
  A2 [-4.26, -3.92]**, with one sign of F_L for L >= 20. Reason: the non-analytic part of the d density at C is -zC r_C h2(r_C), whose
  transform is ~ h2(k) k^-8; against V_A (k^-6) it is ~ h2(k) k^-14 ~ |a|^7 h2(a), a = A - C (the same pattern as
  PC's |a|^7 a_z: a degree-l harmonic at C gives |a|^7 h_l(a)). At x*, |a|^2 = 0, so h2(a) = (3/2) a_z^2 with
  a_z = D2 x* - D1 = (D2^2 - D1^2)/(2 D1) != 0: the 7/2 branch survives with a nonzero coefficient, beta -> -4. The
  earlier offsets from T were all small (AB, AP, AQ within 0.04). What P3 tests: that this coefficient does not
  cancel; if it did, the next branch is 9/2 (beta near -5), and a surviving 5/2 or 3/2 would give beta near -3 or
  -2. **Disclosed limits.** The coefficient is proportional to a_z(x*)^2, which vanishes only when D1 = D2, so P3
  is close to implied by the derivation plus AP and AQ. And P3 cannot tell a d factor from none: AB's own density
  and PC carry the same 7/2 branch, so an engine that dropped h2 or reused AB's density would pass P2a and P3.
  That is caught by P1 (four points) and AD0 (ii). A real test of the h2(a) mechanism is a D1 = D2 geometry, where
  x* = 1 and the leading term moves to |a|^9; it needs an algebraic-decay instrument and is left to a later
  record. Probability 0.75.

**Reported, not scored:** beta_eff at L = 20, ..., 220; the first L from which the set is fixed; signs for L >= 5.

## Decision gate

- **VOID** if AD0 or P1 fails.
- **GO** if P1, P2a and P3 hold: a d function at C keeps a constant set fixed in L and AB's 7/2 branch.
- **STOP** if beta_eff(220) is below -4.6 in either case (the 7/2 coefficient cancels) or above -3.6 (a 5/2 or 3/2
  branch).
- **BORDERLINE** otherwise, including beta inside the band with a sign change of F_L for L >= 20.

## Caps

One d function (3d, m = 0 along B-C) at C, two geometries, two primes; constants assumed linearly independent over
Qbar. Not covered: m != 0 d functions, d functions on A or as pairs, p_pi pairs on A, exchange at l > 0, Paper 63's
fibre at l > 0. The degree-l pattern is argued, not shown, beyond l = 2.

## Phase 2 (after this record is committed)

Add `dens_dz2` to `debug/aapc_exact.py` (integral 'PD'; a multiply-by-t helper applied once and twice; AB's kernel;
normalization 2 pi N_B N_d), an AD0 driver, and 'PD' scoring in `debug/aapc_run.py`. Then AD0; P1 at 2000 and
2600 digits; P2a over L <= 1000 mod both primes; P3 as above.

## Results (2026-11-15; appended after scoring, text above unchanged)

**Verdict: GO.** AD0, P1, P2a and P3 held in both cases.

| clause | A1 | A2 |
|:--|:--|:--|
| AD0 (i) AB reproduced, Q and mod p | yes | yes |
| AD0 (ii) F_0, F_2, F_5, F_12 vs direct projection (rel., tol 1e-18) | 2.9e-24, 3.2e-23, 1.6e-23, 3.5e-24 | 9.0e-24, 2.4e-23, 5.5e-23, 1.1e-22 |
| AD0 (iii) F_L, L <= 220, unchanged at L_max 240 -> 260 | yes | yes |
| P1: value, first point (x = 1/2; -1/3) | 0.104757560074952 | 0.114237095233031 |
| P1: minus R2 n = 34 (tol 1e-10) | 2.0e-13 | 2.3e-13 |
| P1: value, second point (x = 0; 4/5) | 0.119232388104144 | 0.073835256700850 |
| P1: minus R2 n = 34 (tol 1e-10) | 9.1e-14 | 3.3e-14 |
| R1 finest, reported (loose 1e-8) | -1.9e-9, -8.6e-11 | -1.9e-9, -9.7e-10 |
| P2a: set fixed L = 2..1000, both primes | yes, 21 keys | yes, 21 keys |
| P3: beta_eff(220) (A1 [-4.29, -3.95]; A2 [-4.26, -3.92]) | -4.170 | -4.139 |
| P3: one sign of F_L for L >= 20 | yes (+) | yes (+) |

- The 2600-digit run reproduces the 2000-digit values to all 30 printed digits and beta_eff(220) to six. Tail
  terms at L = 240 are below 1e-51.
- beta_eff approaches -4 from below (A1: -5.71 at L = 20, -4.47 at 80, -4.17 at 220), 0.05 below the toy T at
  L = 220 in both cases; the earlier records sat within 0.04.
- P2 reported: exactly AB's 21 at both cases and both primes, complete from L = 1 (L = 0: the 7 exponentials);
  inside C; no gamma. F_L > 0 for L >= 5.
- Drivers: `debug/aadz_ap0.py` (log `debug/data/aadz_ap0.log`), `debug/aapc_run.py p1p3d / p2d`, engine path
  `dens_dz2` in `debug/aapc_exact.py`. Data `debug/data/aadz_p1p3_{2000,2600}.json`,
  `debug/data/aadz_p2_{2305843009213693951,4611686018427387847}.json`.

**Consequence for inventory row 8.** Part 2b-ii is closed for one d function (3d, m = 0 along B-C, at C): AB's
21 constants and the 7/2 branch with a nonzero coefficient. The values at four points test the d factor; the decay
clause does not distinguish it (disclosed above). Owed in 2b: m != 0 d functions, d or p_pi on A, exchange at
l > 0, Paper 63's fibre at l > 0, and the D1 = D2 test of the |a|^7 h_l(a) mechanism.
