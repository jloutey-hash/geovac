# Record AE: (AA|BC) at equal arms, D1 = D2 (inventory row 8: the test of the |a|^7 h_l(a) mechanism)

*Registered 2026-11-15, before computing, under the standing autonomous-loop authorization (PI, 2026-11-14). A
fresh reviewer read the draft cold before commit; its findings were applied: PC's (1 - x)^{9/2} order has three
contributions that partly cancel (the draft said they "can only add"); the leading amplitudes are now predicted in
closed form and scored (P4), which catches a small surviving 7/2 term that an exponent test would pass; exponent
bands widened for the (1 - x)^{11/2} correction; an engine continuity clause added to AE0. The PM re-derived the
reviewer's constants independently before using them (K(5,3)/K(3,3) = 5/22; the PC bracket vanishes at
zB D = 9).*

## Why this record

Records AP, AQ and AD found the same structure for p and d functions at C (AB's 21 constants, the 7/2 branch), and
each was close to implied by the one before. The reason given in Paper 59 `sec:scope` (derived, untested) is that
an m = 0 harmonic of degree l along B-C at C makes the non-analytic part |a|^7 h_l(a), a = A - C, and h_l(a) does
not vanish at the branch point unless D1 = D2. This record tests that reason where it makes a prediction that the
1s case does not share.

## Object

B at the origin, C at (0, 0, D), A at D (sin g, 0, cos g), x = cos g; equal arms D1 = D2 = D. At x = 1, A sits on C,
|a|^2 = 2 D^2 (1 - x). The integral is analytic in x on [-1, 1] except at x = 1 (A is on a sphere about B; the only
non-analyticity is A meeting C), so its Legendre coefficients decay algebraically, set entirely by the branch at 1.

    S  = (1s_A 1s_A | 1s_B 1s_C)       (AB's integral)
    PC = (1s_A 1s_A | 1s_B 2pz_C)      (AP's)
    PD = (1s_A 1s_A | 1s_B 3ds_C)      (AD's)

Cases E1 (D = 2, zB = 13/10, zC = 1, zA = 9/10), E2 (D = 12/5, zB = 1, zC = 6/5, zA = 4/5). Engine: AB's split with
the middle radial region (D1 < r < D2) empty; nothing else changes.

## Controls (measured; `debug/aaeq_controls.py`, `debug/data/aaeq_controls.json`)

| id | what | result |
|:--|:--|:--|
| R2 | Gaussian fits n = 30 / 34 (AB's fits; s, p, d constructions of AB, AP, AD) | E1 x = 1/2: S 0.22145582331466643 / 0.22145582331291747; PC -0.25281152601248086 / -0.25281152601442053; PD 0.14800720824866198 / 0.14800720825235306. E1 x = 0: S 0.1983956434095259 / 0.19839564340799995; PC -0.25139132907428374 / -0.2513913290749442; PD 0.16155157102830853 / 0.1615515710312944. E2 x = -1/3: S 0.12792757651835862 / 0.12792757651481212; PC -0.1949386995615777 / -0.19493869954406912; PD 0.1641310043955077 / 0.1641310043725679. E2 x = 4/5: S 0.18454821479275932 / 0.18454821479429642; PC -0.2136417206344011 / -0.21364172062861544; PD 0.14248327907176628 / 0.1424832790593707 |
| R1 | prolate quadrature 120^2 x 144 | E1 1/2: S 0.2214558243, PC -0.2528115270, PD 0.1480072077; E1 0: 0.1983956441, -0.2513913306, 0.1615515716; E2 -1/3: 0.1279275763, -0.1949386989, 0.1641310036; E2 4/5: 0.1845482146, -0.2136417206, 0.1424832794 |
| T | toy for the algebraic-decay instrument: exact Legendre coefficients of (1 - x)^alpha (Gamma-function closed form, magnitude checked against direct quadrature at L = 1, 2, 3, 6), beta_eff(L) = log\|c_{L+1}/c_{L-1}\| / log((L+1)/(L-1)) | alpha = 7/2: -7.938 (L = 40), -7.966 (100), -7.977 (160), -7.983 (220); alpha = 9/2: -9.944, -9.961, -9.973, -9.980. The coefficients alternate in sign |

R2's n = 30 -> 34 spread is at most 2.3e-11 (PD at E2); the d-type floor at n = 34 is 1.1e-11 (record AD).
**Referee values: R2 at n = 34; tolerance 1e-10 for all twelve.** R1 sits within about 1.5e-9 and is reported only.

## Predictions (none uses an exact split, constant set or decay exponent of S, PC or PD at equal arms)

- **AE0 (engine self-checks, run first).** (i) With D1 < D2 the engine's output is unchanged: AB's coefficients key
  for key at A1/A2 (Q and mod 2^61 - 1, L <= 30), and PC, PD at A1 identical to the committed v6.7.30 engine
  (L <= 20, mod p). (ii) At E1 and E2, F_L of S, PC, PD at L = 0, 3, 8 equals a quadrature with the density
  multipole projected directly (densities as written, no recurrence) against V_A's multipole, relative 1e-18 at 22
  digits; the quadrature splits at r = D, where the kinks from A and C coincide. (iii) Continuity of the engine
  change: F_L at D1 = D - 10^-40 through the unchanged three-region path equals the D1 = D2 path to relative
  1e-35, L <= 20, both cases, all three integrals. Probability 0.9.
- **P1 (values).** The exact split at 2000 digits summed to L = 240 matches R2 (n = 34) to <= 1e-10 at both points
  of both cases, for all three integrals (12 values). The algebraic tail at L = 240 is about 240^-8 relative for S
  and is reported. Probability 0.85.
- **P2 (constant set).** (P2a) For each integral and case, the set of constants with nonzero coefficient in F_L is
  the same at every L >= 2, L <= 1000, mod both primes. (P2c) No gamma. Probability 0.8 / 0.85. Reported, not
  scored: the set and its size (the middle region is gone, so it is not AB's 21), and whether S, PC, PD share it.
- **P3 (exponent).** beta_eff(220) of F_L:
  - **S: in [-8.4, -7.6]** (|a|^7 = (2D^2 (1 - x))^{7/2}, the bare 7/2 branch at x = 1; toy -7.983).
  - **PC and PD: in [-10.4, -9.6]**, both cases (toy at 9/2: -9.980). Reason (local analysis at a -> 0: the
    non-analytic part of the BC density at C against V_A's odd powers, V_A = zA - (2/3) zA^3 s^2 + (2/3) zA^4 s^3 -
    ...; every other pairing is a smooth convolution). PD: the single degree-9 term h2(a)|a|^7 with
    h2(a) = (3 a_z^2 - |a|^2)/2 = -|a|^2/2 + O(|a|^4) on the sphere (a_z = -|a|^2/(2D) exactly), so (1 - x)^{9/2}.
    PC: three terms at (1 - x)^{9/2}: a_z |a|^7 from the cusp, and from e^{-zB r} expanded about C (factor
    -zB u_z) an isotropic |a|^9 and an h2(a)|a|^7 part; they partly cancel, their sum is proportional to
    (zB D - 9), and E1 (zB D = 2.6) and E2 (2.4) are far from the zero. No (1 - x)^{7/2} term exists for PC or PD
    (it would need a bare |u| in the density, which the angular factor removes). B-gradient corrections beyond
    these, and V_A's s^5 term, start at (1 - x)^{11/2}; the bands allow a (1 - x)^{11/2} admixture up to about 80
    times the leading coefficient (reviewer's toy mix).
  - **The wrong alternative**: if the extra factor did not vanish (the reason in Paper 59 is wrong), PC and PD would
    carry a (1 - x)^{7/2} term and decay like S, beta near -8. The two readings are 2 apart in beta; the toy puts the
    instrument's finite-L offset at 0.02.
  Probability for P3: 0.75 for PC and PD, 0.85 for S.
- **P4 (amplitude; the sharp clause).** The leading coefficient is predicted in closed form (controls `amplitude`):
  S: N_B N_C e^{-zB D} (-zC) v3 K(1,3) |a|^7; PC: N_B N_C' e^{-zB D} zC v3 [3 K33/(2D) + zB (K33 - (99/15) K53)/3]
  |a|^9; PD: N_B N_d e^{-zB D} zC v3 (99/30) K53 |a|^9; v3 = (2/3) zA^4, K(m,n) = c_m c_n / c_{m+n+3} with c_q the
  3-D Fourier constant of |r|^q; F_L^pred = amplitude (2D^2)^alpha c_L(alpha). **F_220 / F_220^pred lies in
  [0.9, 1.1]** for all three integrals in both cases (predicted F_220: E1 S 5.794e-17, PC -2.747e-20, PD
  -1.821e-20; E2 S 1.685e-16, PC -1.186e-19, PD -1.098e-19). The (1 - x)^{11/2} correction is estimated
  at a few percent at L = 220. A surviving (1 - x)^{7/2} term of relative size 1e-4 in PC or PD would move the
  ratio by about 10%; a dropped B-gradient term moves PC's by about 40%. Probability 0.7.

**Reported, not scored:** beta_eff at L = 20, ..., 220; the signs of F_L (the toy alternates); the ratio F_L(PC) /
F_L(S) and F_L(PD) / F_L(S), which the mechanism predicts to fall like L^-2.

## Already known (risk of a void prediction, disclosed)

- The S prediction is the classical one: the same mechanism as Paper 59's water bond-angle cusp-cusp term
  (|Y - Z|^7, l^-8 at equal lengths), in a different integral. Its live content is that the screened 1s potential
  of A leaves that branch unchanged.
- P3 for PC and PD has no control beyond the toy: no decay exponent of any integral at equal arms has been computed.
- The engine change (empty middle region) is exercised by AE0 (ii) and P1, not independently.

## Decision gate

- **VOID** if AE0 or P1 fails.
- **GO** if P1, P2, P3 and P4 hold for all three integrals in both cases: the local analysis predicts the
  equal-arm decay and its amplitude, and p and d functions differ from 1s there by L^-2.
- **STOP** if PC or PD has beta_eff(220) above -8.5 in either case (a (1 - x)^{7/2} term survives: the reason is
  wrong), or S is outside [-8.5, -7.5].
- **STOP** also if P3 holds but P4's ratio is outside [0.5, 2] for any integral (the exponent is right for a wrong
  reason).
- **BORDERLINE** otherwise (including PC or PD between -9.6 and -8.5, or below -10.4; P4 ratio in [0.5, 0.9) or
  (1.1, 2]).

## Caps

Two equal-arm geometries, two primes, three integrals (1s, 2pz, 3d-sigma at C, each against a 1s pair at A); the
degree-l pattern tested at l = 1 and 2 only. Constants assumed linearly independent over Qbar.

## Phase 2 (after this record is committed)

In `debug/aapc_exact.py`, allow D1 = D2 (drop the empty middle region; assert D1 <= D2). A driver
`debug/aaeq_run.py` registers cases E1, E2 in its own process only, runs AE0, P1 at 2000 and 2600 digits, P2 mod
both primes, and P3.
