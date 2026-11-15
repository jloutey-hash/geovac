# Record AE2: equal arms at fresh geometries, with the PC cancellation (follow-up to record AE)

*Registered 2026-11-15, before computing. Follow-up chosen by the PI after record AE's STOP (AE's amplitude clause
failed on a sign convention: the registered Legendre sign rule (-1)^L, read off small L, is wrong for L > alpha + 1;
magnitudes agreed within 0.2%; `docs/predictions/2026-11-15_aa_bc_equal_arms.md`, Results). AE2 re-tests the
closed-form amplitudes with the corrected sign rule at geometries AE never computed, and adds one prediction AE did
not make. A fresh reviewer read the draft cold before commit and its findings were applied: a precision gate
(VOID) and a drift test for PC at E3 (both MUST-FIX); PC's (1 - x)^{11/2} amplitude at E3 registered from a general
local expansion (the reviewer's code, `debug/aaeq2_local.py`, read and run by the PM; it reproduces the PM's
closed-form AE amplitudes at E1 to 12 digits); the correction estimate fixed (it was about 50 times too large) and
P4 tightened to 2%; the E3 tolerance set on R2 alone. The reviewer confirmed the corrected sign rule and that the
PC bracket is complete at (1 - x)^{9/2}.*

## Object

As in AE: B at the origin, C at (0, 0, D), A at D (sin g, 0, cos g), x = cos g, D1 = D2 = D; S = (1s_A 1s_A|1s_B 1s_C),
PC = (1s_A 1s_A|1s_B 2pz_C), PD = (1s_A 1s_A|1s_B 3ds_C). Engine: AE's (empty middle region), unchanged.

Fresh cases:
- **E3: D = 3, zB = 3, zC = 1, zA = 9/10, so zB D = 9.** AE's closed form for PC's (1 - x)^{9/2} coefficient is
  proportional to (zB D - 9) (three terms: a_z|a|^7 from the cusp, and from e^{-zB r} about C an isotropic |a|^9 and
  an h2(a)|a|^7 part; K(5,3)/K(3,3) = 5/22). At E3 it vanishes exactly (evaluates to -5e-54 at 50 digits).
- **E4: D = 5/2, zB = 7/10, zC = 13/10, zA = 11/10** (zB D = 1.75), generic.

## Controls (measured; `debug/aaeq2_controls.py`, `debug/data/aaeq2_controls.json`)

| id | what | result |
|:--|:--|:--|
| T | toy c_L(alpha) = (2L+1)/2 int (1-x)^alpha P_L, exact, sign rule sign(Gamma(alpha+1)/Gamma(alpha+1-L)) (-1)^L checked against quadrature at L = 1, 2, 5, 6, 7, 10, 21 (all match in sign and magnitude); beta_eff at L = 220 / 400 / 580 | 7/2: -7.983, -7.990, -7.993; 9/2: -9.980, -9.988, -9.992; 11/2: -11.977, -11.986, -11.990 |
| A | closed-form leading amplitudes (AE's formulas, corrected sign), predicted F_580 | E3: S 2.4977e-21, PD 2.5473e-25, PC 0 (bracket -5e-54). E4: S 4.7570e-19, PC 5.9876e-23, PD 5.6938e-23 |
| R2, R1 | Gaussian fits n = 30 / 34 and prolate 120^2 x 144, two angles per case (values in the json) | E3: n = 30 -> 34 moves values by 2.8e-10 to 8.5e-10 (B is compact at zB = 3); R1 agrees with n = 34 to 4e-11 to 1.8e-10. E4: spread <= 1.9e-10; R1 within 7e-9 |
| X | general local expansion (`debug/aaeq2_local.py`: e^{-zB r} to degree 5, odd part of e^{-zC r_C}, V_A's odd terms to s^7), coefficients of t^{7/2}, t^{9/2}, t^{11/2}, t = 1 - x | E3: S 0.033048, 0.069269, 0.082152; PC (9/2) -1.7e-18 (float64 rounding in one factor; the 50-digit closed-form bracket is -5e-54), (11/2) -0.0048413; PD -0.028042, -0.062195. E4: S 6.2941, 13.681, 14.014; PC -6.5913, -14.470; PD -6.2680, -6.1094. Predicted F_L for PC at E3: A_{11/2} c_L(11/2) = -8.767e-25 (L = 220), -6.791e-28 (400), -7.8966e-30 (580) |

**Referees: R2 n = 34. Tolerances: E4 2e-10, E3 1e-9, both against R2 alone** (at E3 the fit moves 8.5e-10
between n = 30 and 34; R1 is reported).

## Predictions (none uses an exact split or decay exponent at E3 or E4)

- **AE2-0 (self-checks, run first).** (i) Engine file unchanged since AE's scoring commit. (ii) Per-channel
  direct-projection quadrature at E3 and E4 (S, PC, PD; L = 0, 3, 8; split at r = D), relative 1e-18 at 22 digits.
  (iii) Continuity at D1 = D - 1e-40, relative 1e-35, L <= 20, at 400 digits. Probability 0.9.
- **P0 (precision gate).** F_579, F_580, F_581 of every integral and case, at 2000 and 2600 digits, agree to
  relative 1e-10. Failure is VOID.
- **P1 (values).** Exact split at 2000 digits summed to L = 600 matches R2 within the tolerances above at both
  points of both cases (12 values). Probability 0.85.
- **P3 (exponents) at L = 580** (toy offsets at 580 are under 0.01):
  - S, both cases: [-8.25, -7.75].
  - PD, both cases, and PC at E4: [-10.25, -9.75].
  - **PC at E3: beta_eff in [-12.3, -11.7] at each of L = 220, 400 and 580, with |beta(580) - beta(220)| < 0.05**
    (the new prediction: with the (1 - x)^{9/2} coefficient zero, PC leads with (1 - x)^{11/2} and decays like
    L^-12, where PD keeps L^-10 and S L^-8; the drift test catches a small surviving (1 - x)^{9/2} term, which would
    pull beta toward -10 as L grows; the toy drifts 0.013 over that range).
  The next-order to leading coefficient ratios from X are about 2 (S, PD at E3) to 2.2 (E4); with
  c_L(alpha+1)/c_L(alpha) ~ -2 (alpha+1)^2 / L^2 the correction at L = 580 is under 0.1% for S and PD (the draft's
  "2 to 4%" estimate was about 50 times too large). Probability 0.75 (PC at E3: 0.65).
- **P4 (amplitudes, corrected sign).** F_580 / F_580^pred in **[0.98, 1.02]** for S, PC, PD at E4 and S, PD at E3
  (five ratios; AE's magnitudes matched to 0.2% at L = 220), and **in [0.9, 1.1] for PC at E3** against the
  registered (1 - x)^{11/2} amplitude (its (1 - x)^{13/2} correction is not computed and may be larger after the
  cancellation). Probability 0.7.

**Reported, not scored:** beta_eff at L = 220, 400, 580 for all six; F_L(PC)/F_L(PD) at E3 (predicted to fall like
L^-2); the constant sets (AE showed one fixed set of 17 per case; not re-run).

## Decision gate

- **VOID** if AE2-0, P0 or P1 fails.
- **GO** if P0, P1, P3 and P4 hold: AE's local analysis predicts the equal-arm exponents and amplitudes at fresh
  geometries, and the predicted cancellation at zB D = 9 is real.
- **STOP** if PC at E3 has beta_eff(580) above -10.5 (the cancellation does not happen) or drifts by more than 0.3
  toward -10 between L = 220 and 580, or any P4 ratio is outside [0.5, 2] or negative, or S, PD or E4's PC has beta outside [-8.6, -7.4] / [-10.6, -9.4].
- **BORDERLINE** otherwise.

## Caps

Two fresh equal-arm geometries, three integrals, l <= 2 at C. The (1 - x)^{11/2} order is predicted in exponent only.
AE's own verdict stands as registered (STOP); AE2 does not rescore AE's data. The (1 - x)^{13/2} order is not
computed.

## Phase 2 (after this record is committed)

`debug/aaeq2_run.py`: registers E3, E4 in its own process, runs AE2-0, then P1/P3/P4 at L <= 600 at 2000 and 2600
digits, reading predictions from `debug/data/aaeq2_controls.json`.
