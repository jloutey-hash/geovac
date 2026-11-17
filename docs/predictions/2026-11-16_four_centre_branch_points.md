# Pre-registered: the four-centre integral (XY|ZW) as a function of one angle has two contact branch points, each (x* - x)^{7/2} (record FC)

**Registered:** 2026-11-16, before any scored computation. The PI chose the route: "numeric-only first", over the exact
matched-exponent slice. The cold review before commit raised one MUST-FIX, now applied: the old FC3 (one sign) could
not fail, so it moved into FC0 as an instrument check. It also raised three SHOULD-FIX items, now applied: the
pre-Landau-point assumption is stated, control K4ext extends the exact AB reduction to L = 220, and the
rate-sensitivity formula is corrected.
**Inventory:** row 7 of `docs/chemistry_period_inventory.md`, column 3 (singular points), which is MISSING.
**Owner:** Paper 63 `sec:four_centre`.

## Object

Four distinct centres and unnormalised 1s Slater products:
- rho_XY = e^{-zX|r-X|} e^{-zY|r-Y|} and rho_ZW = e^{-zZ|r-Z|} e^{-zW|r-W|};
- I(x) = (rho_XY | rho_ZW) = sum_L F_L P_L(x).

**The one-parameter slice:**
- X sits at the origin and Y = D1 z-hat.
- Z = a n-hat and W = b n-hat lie on one ray from X, with 0 < a < b and x = n-hat . z-hat.
- Equivalently, the ray XZW turns rigidly about X against the fixed bond XY.

The exponents are generic; nothing is matched.

Both densities are axial about axes through X, so F_L is a two-dimensional radial integral of their Legendre
multipoles about X: F_L = (4 pi/(2L+1))^2 ∫∫ rho1_L rho2_L r<^L / r>^{L+1} r1^2 r2^2. The multipole rho1_L is exact
(one factor on X; i_L/k_L form). The multipole rho2_L is a high-precision angular quadrature on panels graded toward
the cusp singularities, done in arb at 340 bits (level A) and 440 bits (level B). Engine: `debug/fc_engine.py`,
self-contained.

## Theory (why two points, why 7/2)

- The geometry enters only through the distances |Y - Z| and |Y - W|. X-Z and X-W are fixed; Y-X and Z-W are within a
  density. That gives two contact (Landau) points in x:
  - x_Z* = (D1^2 + a^2)/(2 D1 a), with rho_Z = min(D1, a)/max(D1, a);
  - x_W* = (D1^2 + b^2)/(2 D1 b), with rho_W = min(D1, b)/max(D1, b).
- There are no point charges, so the leading non-analytic term at a contact is cusp against cusp through 1/r12. In
  momentum space that is k^-2 k^-4 k^-4, which gives |Y - W|^7.
- Since |Y - W|^2 = 2 D1 b (x_W* - x), the term is (x* - x)^{7/2}. In decay terms, F_L ~ rho^L L^-4.
- Record B2X measured this cusp-cusp |.|^7 term in the shared-centre case (XY|XZ). Record AB found a leading 7/2 in
  (AA|BC).
- Nothing at 5/2 or 3/2 is expected:
  - the cusp at X never meets a cusp of rho_ZW;
  - a cusp against the other density's smooth part is analytic;
  - the first correction to cusp-cusp is |R|^7 times an analytic factor;
  - the amplitude is a product of positive cusp strengths, so it does not cancel.
- **Positivity (known in advance, not a prediction).** rho1_L and rho2_L are >= 0 for every L, because each density
  is absolutely monotone in the angle cosine, and the radial kernel is positive. So F_L > 0 for every L by
  construction. By a Pringsheim-type argument, positive Legendre coefficients put the dominant singularity on the
  real axis at x > 1, which excludes complex-x competitors.
- **Assumption (stated, not derived):** no real singularity in (1, x_dom*) other than the two contacts. Real
  geometry has singularities only at coincidences of centres, but x > 1 is complexified geometry, so this does not
  follow automatically. The [-7, -2] STOP clause is what catches a violation.

## Cases

| case | a, b, D1 | zX, zY, zZ, zW | x_Z*, rho_Z | x_W*, rho_W | dominant |
|:--|:--|:--|:--|:--|:--|
| G1 | 1, 9/4, 8/5 | 13/10, 1, 6/5, 1 | 89/80, 5/8 | 3049/2880, 32/45 | W |
| G2 | 1, 4, 6/5 | 13/10, 1, 6/5, 7/10 | 61/60, 5/6 | 109/60, 3/10 | Z |

The W contact dominates G1 and the Z contact dominates G2, so the two cases together test both points. The
subdominant point is suppressed by (5/8 / 32/45)^220 ~ e^{-28} at G1 and far more at G2. It is not scored.

## Controls (measured by the engine agent; `debug/data/fc_controls.json`)

The engine agent was barred from coefficient values at L >= 40, from beta_eff of the engine's sequence at L >= 40,
and from any branch-point amplitude. Nothing barred was printed or stored.

| id | what | result |
|:--|:--|:--|
| K1 | rho2_L at L = 0, 1, 4, 12, 30, r on and around a, b, D1, vs 2-D mp.quad at 60 digits; rho1_L likewise | worst relative 5.6e-51 (rho2), 4.9e-40 (rho1) |
| K2 | referee R2: McMurchie-Davidson on even-tempered n-Gaussian fits, n = 22 vs 26, at x = 1/2, -1/3 | spread 3.2e-12, 1.8e-10 (G1); 1.7e-9, 1.2e-10 (G2) |
| K3 | sum to L = 260 minus R2(n = 26) | G1 -1.1e-11, -2.5e-11; G2 -6.2e-11, -4.1e-11 (referee-limited); level A - B on the sums 8e-36 (G1), 2.7e-30 (G2) |
| K4 | reduction b = a, Z = W, mapped onto record AB's case A1, vs AB's exact `aabc_exact` coefficients | L <= 30, worst relative 1.1e-28; value at x = 1/2 0.21286449876 vs AB's referee 0.2128644988 |
| K5 | precision: level A vs B, which differ in bits, theta panels, radial grid, sub-nodes and cutoff | relative self-difference 1.7e-23 .. 3.2e-23 (L = 219..221), 6.0e-21 .. 9.1e-21 (L = 258..260), both cases |
| K4ext | (added after the cold review) the engine's beta_eff on K4A1 = AB's case A1, at L = 20, 40, ..., 220, vs AB's exact-route beta_eff stored in `debug/data/aabc_p1p3b_2600.json` (`debug/fc_k4ext.py`, `debug/data/fc_k4ext.json`) | max abs difference 4.1e-16 at both levels (float precision of AB's stored values); all F_L > 0 |
| K6 | toys: Legendre coefficients of (x* - x)^nu at each case's dominant x*, beta_eff at L = 60, 120, 220 | G1: nu = 5/2 -3.265, -3.136, -3.075; 7/2 -4.450, -4.231, -4.128; 9/2 -5.684, -5.351, -5.194. G2: 5/2 -3.473, -3.251, -3.140; 7/2 -4.790, -4.421, -4.236; 9/2 -6.180, -5.634, -5.356 |

Timing: about 3 minutes per case for both levels to L = 260.

## Predictions

The scored run is `python debug/fc_run.py G1` and `... G2`, giving `debug/data/fc_run_<case>.json`. Beta_eff is
measured on G_L = F_L / rho_dom^L, from level B.

- **FC0 (instrument gate):**
  - the relative self-difference of A and B is at most 1e-15 at L = 219..221 and 258..260 (a reproducibility
    check of K5);
  - |beta_eff A - B| at L = 200, 220 and 240 is at most 1e-6;
  - every F_L > 0 (positivity is forced, so a negative coefficient is a numerical bug);
  - the frozen files are unchanged.

  Failure is VOID. The large-L normalisation was checked independently by K4ext; FC0 does not re-test it.
- **FC1 (the branch, both cases):** beta_eff(220) lies between the toy midpoints around 7/2:
  - G1 in **[-4.661, -3.602]**;
  - G2 in **[-4.796, -3.688]**.

  Probability 0.75.
- **FC2 (settled):** beta_eff at L = 200 and 240 lie in the same bands, both cases. The bands reuse the L = 220
  toy midpoints; the toys drift by about 0.01 over that range, which is harmless. Probability 0.7.

**Reported, not scored:**
- beta_eff at L = 20..240;
- the amplitude fits, with beta fixed at -4 and with beta free;
- the K-series of the fits.

## Decision gate

- **VOID:** FC0 fails, a frozen file changed, or a step has not finished 48 h after launch.
- **GO:** FC0, FC1 and FC2 hold.
- **STOP:** at either case, beta_eff(220) is above the 5/2-7/2 midpoint (-3.602 at G1, -3.688 at G2), meaning a
  stronger singularity than cusp-cusp. STOP also applies if beta_eff(220) lies outside [-7, -2]: the rate rho_dom is
  then wrong, so another singular point dominates or the contact point is absent. Sensitivity: a rate error
  rho_true/rho_dom shifts beta_eff by about L log(rho_true/rho_dom), so at L = 220 this clause catches a rate error of
  about 1% or more. The subdominant point would shift beta_eff by about -28.
- **BORDERLINE:** otherwise. This includes beta_eff below the 7/2-9/2 midpoint, a cancellation of the 7/2 amplitude.

## Caps

- 1s only, unnormalised products, one family of slices (X, Z, W collinear on one ray), two geometries.
- A decay exponent read at L <= 240, not an operator. The value class and an exact recurrence (the row's other
  columns) are not addressed. The exact matched-exponent slice stays open.
- The subdominant point is not tested where it is subdominant. Each point is tested where it dominates.

## Frozen (this commit; unchanged at scoring)

| file | sha256 |
|:--|:--|
| `debug/fc_engine.py` | `9fef5cf60eeec833283f2211050af56e38d163b71b4fcf6a6bbe1744f56707aa` |
| `debug/fc_run.py` | `5a1f63a96e4c99f4a6dacc7bd45af80114eb2660e77d64d3aee77cac7c9cd5ff` |
| `debug/fc_controls.py` | `648bd789cd2112f64094ef49df98c243cd2fa6a0a2c3b53160a0ae7fd8c67942` |
| `debug/data/fc_controls.json` | `0acb0902b80fb2baa6eb7cc437dd5b4f9fd933e987f7d16aa601e07ec36b02b0` |
| `debug/fc_k4ext.py` | `ae8c30e0c6240b55295e3baf0c989da792334df7654128503a68f01aa8f1987c` |
| `debug/data/fc_k4ext.json` | `c990ba6a4a9bc0a22be836959416fcc046b290b656cb7da54d96e26ecadb1a48` |

Environment: Python 3.14.0, python-flint 0.9.0, mpmath 1.3.0. The controls also import `na_exact`, `aabc_exact` and
`aabc_controls`, which are not frozen here; that is acceptable because the controls were already measured. The
scoring driver imports only `fc_engine`.

## Phase 2 (after this record is committed and the commit verified)

`python debug/fc_run.py G1` and `python debug/fc_run.py G2`.

## Results (scored 2026-11-16; record text above unchanged)

The scored runs were launched after commit 378fb8ca was verified. The frozen files show no diff at scoring. Each case
took about 3 minutes.

**Verdict: GO (FC0, FC1, FC2 held).**

- **FC0 held.**
  - Relative A-B self-difference: at most 2.4e-23 (G1) and 3.2e-23 (G2) at L = 219..221; at most 7.8e-21 and
    9.1e-21 at L = 258..260.
  - |beta_eff A - B| at L = 200, 220, 240: at most 2.1e-20.
  - Every F_L > 0 in both cases.
- **FC1 held.** beta_eff(220) is **-4.1226** at G1 (band [-4.661, -3.602]; the 7/2 toy gives -4.128) and **-4.2357**
  at G2 (band [-4.796, -3.688]; toy -4.236).
- **FC2 held.** beta_eff(200, 240) is -4.1346, -4.1126 at G1 and -4.2584, -4.2166 at G2.

**Reported, not scored:**
- beta_eff at L = 20, 40, 60, 80, 100, 120, 160:
  - G1: -5.190, -4.635, -4.429, -4.326, -4.264, -4.221, -4.167;
  - G2: -5.836, -5.106, -4.789, -4.611, -4.498, -4.421, -4.320.
- Free-beta fits (level B, L0 = 130 and 173, K = 3): beta = -3.99998 and -3.99999 at G1; -3.99988 and -3.99993 at
  G2.
- Amplitudes on F/rho^L with beta fixed at -4 (unnormalised densities): 0.049602 (G1, the W contact) and 1.1669e-4
  (G2, the Z contact).

Reading: in both cases the decay follows the toy of a pure (x* - x)^{7/2} branch at the predicted rate, L by L, to
about 0.005 at L = 220. The free fit lands on -4 to 1e-4. Each contact point, Y-W in G1 and Y-Z in G2, is a 7/2
branch point. This fills row 7's singular-point column for the slice, numerically; no operator is claimed.
