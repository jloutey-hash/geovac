# Pre-registered prediction: where the three-centre integral is singular as a function of the bond angle (WA0-WA4)

**Registered:** 2026-10-27, before any computation that scores a prediction in this file. PI-directed ("Yes",
to the bond-angle version of the water question).
**Commit status:** committed on its own, before any computation in this file.
**Rests on:** Paper 59 `sec:momentum` (Eqs. `eq:angular`, `eq:reduced`) and the shared-centre channel sum of
`tests/test_pairwise_three_center.py`.

## The question

Take the shared-centre three-centre integral (XY|XZ) over 1s Slater orbitals, X at the origin, |XY| = D1,
|XZ| = D2, bond angle gamma, x = cos gamma. As an analytic function of x, where are its singular points? In
particular, is the linear geometry (x = -1, CO2-like) singular or special, and where does bent water
(x = cos 104.5 deg = -0.25) sit?

## What I know before computing (disclosed)

- **The earlier pitch was wrong.** I first predicted singular points at the linear angles, because
  |W|^2 = s^2 D1^2 + t^2 D2^2 - 2 s t D1 D2 x is a perfect square at x = +-1. But Paper 59's reduction enters x only
  through j0(k|W|), which is entire in |W|^2. A perfect square alone gives no singularity.
- **Hand reasoning, not checked.** A singularity in x needs the k-integral to diverge for complex W. The decay is
  e^{-k (D1 sqrt(c1) + D2 sqrt(c2))}. Following that to the end points of the Feynman parameters puts the onset at
  (Y - Z)^2 = 0, i.e. x* = (D1^2 + D2^2)/(2 D1 D2): the two outer centres coalescing. For D1 = D2 that is x* = 1
  (gamma = 0). For D1 != D2 it is real and above 1, outside the physical range.
- **Analogy, not a derivation.** 1/|Y - Z| has Legendre coefficients D<^l / D>^(l+1) exactly. Its only singularity in
  x is x*, and its Bernstein-ellipse parameter is x* + sqrt(x*^2 - 1) = D>/D<.
- **Method.** The shared-centre channel sum is exactly a Legendre series:
  (XY|XZ) = sum_l F_l P_l(x), F_l = (2(2l+1)/pi) int dk A_l^{XY}(k) A_l^{XZ}(k) (the test file's formula). So the
  decay of F_l locates the nearest singularity of the function of x. Geometric decay r^l means analytic inside the
  Bernstein ellipse of parameter 1/r. A sign pattern (-1)^l means the singularity is on the x < 0 side.
- Paper 59 already reports that the water channel sum converges to about 1e-9 by l = 12 at the water geometry.
- Prior art: Ozdogan and Ruiz treat three-centre Slater repulsion integrals analytically for LINEAR geometries
  (Paper 59, `ozdogan2012`). Whether that reflects a property of the function or of the coordinates is part of this
  question.

## Setup

- Evaluator: the channel formula of `tests/test_pairwise_three_center.py`, with `G` in closed form, re-implemented
  in a new driver `debug/wa_channels.py` with quadrature tightened as needed. No change to `geovac/` or `tests/`.
- Geometries: symmetric water D1 = D2 = 1.809. Asymmetric pairs (D1, D2) = (1.809, 2.5) and (1.809, 1.0).
- Exponents: set E1, all 1s with zeta = 1. Set E2, zeta = 2.25 on X and zeta = 1.24 on Y and Z (O valence-like and
  H-like; a proxy, not a fitted basis).
- l range: 0 to 40, or until |F_l| reaches the quadrature noise floor, whichever comes first. The noise floor is
  measured by halving the grids. Only l with |F_l| at least 100 times the floor are used.
- Fits: log|F_l| = a - p log l + l log r over the usable l at or above 15.

## Predictions

**WA0 (controls; run FIRST; WA1-WA4 void if any fails).**
- (a) The channel sum at the test geometry (X = 0, Y = (0,0,2), Z = (1.5,0,0.5), zeta = 1) reproduces the test's
  GROUND_TRUTH 0.2049417218 to 1e-8.
- (b) The fit applied to the exact coefficients D<^l/D>^(l+1) of 1/|Y - Z| recovers r = D</D> to 1e-6 at both
  asymmetric pairs, and p = 0 and r = 1 at D1 = D2.
- (c) Applied to 1/|Y + Z| (coefficients (-1)^l D<^l/D>^(l+1)), the alternation test of WA3 fires.

**WA1 (asymmetric: the nearest singularity is the coalescence point).** At both asymmetric pairs and both exponent
sets, the fitted r equals D</D> within 1 percent: 0.7236 for (1.809, 2.5) and 0.5528 for (1.809, 1.0). Equivalently,
the nearest singularity of the function of x is at x* = (D1^2 + D2^2)/(2 D1 D2), and it does not depend on the
exponents.
- Falsifier: r differs from D</D> by more than 1 percent in any of the four cases.
- Probability 0.7.

**WA2 (symmetric water: singular at gamma = 0, algebraic decay).** At D1 = D2 = 1.809, for both exponent sets, the
fitted |log r| is below 0.01 (no geometric decay), and the log-log slope p is stable to +-0.5 between the halves of
the usable l range. p is reported, not predicted.
- Falsifier: |log r| at or above 0.01, or p unstable.
- Probability 0.6. The risk is that a high-order singularity at x = 1 looks geometric over a short l range.

**WA3 (the linear geometry is not singular).** In all six cases, F_l has no alternating component. Fit
F_l = S_l + (-1)^l T_l with S and T smooth in l, using adjacent triples. Then |T_l / S_l| < 1e-3 for usable l at or
above 15. Read with WA1 and WA2, this puts the nearest singularity on the gamma = 0 side, not at or near
gamma = 180 deg.
- Falsifier: |T_l / S_l| at or above 1e-3 in any case.
- Probability 0.75.

**WA4 (reported, not scored).** Where bent water sits: the Bernstein-ellipse distance from x = -0.25 to the nearest
singularity in each case. The p of WA2 and the lowest odd power of |Y - Z| it implies (an x-singularity (1 - x)^alpha
gives coefficients about l^(-2 alpha - 1)).

## Decision gate

- **GO:** WA0 passes and WA1 and WA3 hold. Then the bond angle has one finite singular point, the coalescence of the
  outer centres. The linear geometry is an ordinary point, and bent water is a generic one.
- **BORDERLINE:** WA1 holds but WA3 fails, or WA2 is ambiguous. The singular point is located, but its side or
  type is not established.
- **STOP:** WA0 fails, or WA1 fails. A failure of WA1 with clean controls is itself a finding: a singularity in x
  other than coalescence.

## Caps

- At most 2 h wall, 8 processes. Quadrature refinement at most two halvings per case.
- If the noise floor leaves fewer than 15 usable l in a case, that case is void and reported.

## Novelty status (before computing)

- Legendre-coefficient decay as a locator of singularities (Bernstein ellipse, Darboux): KNOWN, standard.
- Multipole (Barnett-Coulson, Lowdin alpha-function) expansions about the shared centre: KNOWN.
- Linear-geometry analytic treatment: KNOWN (Ozdogan and Ruiz).
- The singularity structure of the three-centre Slater integral in the bond angle: NOT SEARCHED. Treat as OPEN.

## Addendum 1 (2026-10-27, before any computation; predictions above unchanged)

- WA0(a) names the wrong reference. 0.2049417218 is the 8-Gaussian McMurchie-Davidson value, which carries the fit
  error. The channel sum computes the Slater integral, whose referee in the test file is REFEREE_NG16 =
  0.204941652260 (Paper 59: the Slater value is 7e-8 lower). WA0(a) as registered would fail on that 7e-8 for a
  reason unrelated to the method. Found by the PM reading the test file, before running anything.
- WA0(a) is scored against REFEREE_NG16 to 1e-9 (the test's own tolerance). The disposition also reports the
  literal registered comparison.

## Disposition

(Appended after computing. The predictions above are not edited.)

### Disposition (appended 2026-10-27 by the PM; computed in the main session)

Driver `debug/wa_channels.py` (modes `assign`, `control`, `cases`), scorer `debug/wa_analyze.py`; data
`debug/data/wa_channels.json`; logs `debug/data/wa_{assign,control,cases,analyze}.log`. Wall time about 1 min.

**Before scoring.**
- Exponent assignment: the test's G has zA on X. The on-axis density transform at zX = 2.25, zH = 1.24 agrees with a
  direct 2D quadrature (0.25230424 vs 0.25230424 at k = 0.7; 0.06075227 at k = 2), and the swapped assignment does
  not. E2 uses the test's assignment.
- *Deviation, disclosed.* My first WA3 estimator fitted 4 parameters to 4 consecutive l. That is interpolation, and it
  read the curvature of a plain geometric sequence as alternation (4e-3 and 2e-2 on 1/|Y - Z|). It was replaced,
  before any case was computed, by the estimator the record names: detrend by the WA1 fit, then the second
  difference on adjacent triples. On the controls it gives about 1e-15 for 1/|Y - Z| and 3.0 for 1/|Y + Z|.

| Item | Result | Verdict |
|:--|:--|:--|
| WA0(a) | channel sum 0.204941652531; REFEREE_NG16 diff 2.7e-10 (below 1e-9). Literal registered reference: diff -6.9e-8 | PASS as amended (Addendum 1); FAIL literally |
| WA0(b) | 1/|Y-Z|: r = 0.72360000 and 0.55279160 (exact); D1 = D2: r = 1, p = 0 | PASS |
| WA0(c) | 1/|Y+Z|: alternation fires (3.0) in all three | PASS |
| WA1 | E1: r = 0.74243 (+2.6%), 0.56848 (+2.8%); E2: 0.74139 (+2.5%), 0.56440 (+2.1%). Usable l 15..36 and 15..23 | **MISSED** (all four over 1%) |
| WA2 | E1: log r = -0.0083, p = 7.60, halves 7.50 / 7.72. E2: log r = -0.0107, p = 7.50, halves 7.38 / 7.67 | **MISSED** (E2 log r over 0.01; E1 holds) |
| WA3 | alternation at most 3.5e-4 in all six; every F_l at l >= 15 positive | HELD |

**Gate: STOP as registered** (WA1 failed). The registered reading of a WA1 failure was "a singularity other than
coalescence". The evidence below points instead to the registered fit model being too crude, but that evidence was
not registered.

**Unregistered follow-up (does not change any score).** `debug/wa_unregistered.py`, log
`debug/data/wa_unregistered.log`. E1 only, grid level 2, l to 60.
- Adding a 1/l term to the fit: r = 0.72983 (+0.86%) and 0.55589 (+0.56%) for the asymmetric pairs. Symmetric water:
  r = 1.00057, p = 8.05.
- In every registered case the upper half of the l range gave r closer to D</D> than the lower half.
- E1 and E2 rates agree to a few tenths of a percent, so the rate does not depend on the exponents.
- Reading (not established): the nearest singularity is the coalescence of the outer centres, approached slowly by
  the l-range available. For symmetric water p near 8 would be a (1 - x)^(7/2) singularity at gamma = 0, i.e. the
  lowest odd power of |Y - Z| being the seventh.

**WA4 (reported).** Bent water, x = -0.25, sits in the interior of the analyticity region in every case; the
symmetric-water singular point is at gamma = 0 (H atoms merging), not at the linear geometry.

**What stands.** WA3: no alternating component, so the linear geometry gamma = 180 deg is not the nearest
singularity in any case tested, and the nearest one is on the gamma = 0 side. WA1's location is not confirmed at the
registered tolerance. A clean test needs a fit model with subleading terms and more digits at higher l, registered
in advance.
