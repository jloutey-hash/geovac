# Pre-registered prediction: the Landau-pinch size and the Gaussian crossover on fresh cases, controls measured first (WX)

**Registered:** 2026-10-28, before any F_l of a case in this file is computed. PI-directed ("Yes", to one fresh
registration covering both tests, on new cases, with every control measured beforehand on cases that are not scored).
**Commit status:** committed on its own, before any computation it scores.
**Rests on:** WL (`2026-10-28_water_landau_unequal.md`, STOP: control (c) bound set without measuring) and WG
(`2026-10-28_water_gaussian_contact.md`, STOP: control (b) precision not set). The derivations, the extraction rule and
the crossover model are theirs, unchanged. This record re-tests them on fresh cases with measured controls.

## Part L: the size of the singularity at unequal bond lengths

**Claim tested (WL's derivation).** The leading bond-angle singular term at x* = (D1^2 + D2^2)/(2 D1 D2) is
Phi (2 D1 D2)^(7/2) (x* - x)^(7/2), with
Phi = (pi^2/630) 24 zH^2 NH^2 NX^2 int int s u exp(-zX (R1(s) + R2(u))) ds du,
the two cusp functions averaged over the segment between the H atoms (`debug/wl_local.py`, unchanged). The naive size
uses the endpoint s = u = 0.

**Fresh cases.** None has been computed. Both bond lengths and the H exponent now vary, not only D2.

| Case | zX | zH | D1 | D2 | x* | Phi | naive | Phi/naive |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| P | 4 | 1.0 | 1.6 | 2.2 | 1.0511364 | 3.18633105060983e-8 | 2.5442706004598e-8 | 1.2524 |
| Q | 3 | 1.5 | 2.0 | 1.5 | 1.0416667 | 9.54576814901055e-6 | 8.96163909641325e-6 | 1.0652 |
| R | 3 | 1.24 | 1.5 | 2.1 | 1.0571429 | 2.81467987214106e-6 | 2.56301551974454e-6 | 1.0982 |
| S | 4 | 1.5 | 2.1 | 1.6 | 1.0372024 | 3.365406376737e-7 | 2.88228809305434e-7 | 1.1676 |

**Extraction (WL's, unchanged).** Linear least squares onto the exact Legendre coefficients of (x* - x)^(7/2 + j),
j = 0..nb-1, over l = 16, 18, ..., 40, rows divided by the j = 0 coefficient. A = the j = 0 coefficient over
Phi (2 D1 D2)^(7/2). **Change from WL:** the F_l are channel sums at grid level 3 (`F_coeffs(..., level=3)` of
`tests/test_paper59_water_bond_angle.py`), and A is scored at nb = 5.

**Controls measured before registering** (`debug/wx_measure_L.py`, data `debug/data/wx_measure_L_{A,C,D}.json`; the
three old WL cases, already voided, not scored here):

| Old case | level 2 vs 3, max over l = 16..40 | level 3 vs 4 | A5 at level 2, 3, 4 |
|:--|:--|:--|:--|
| A | 3.1e-6 | 1.2e-8 | 0.999989, 1.000011, 1.000011 |
| C | 5.1e-6 | 2.5e-7 | 1.000016, 1.000015, 1.000017 |
| D | 8.1e-6 | 4.1e-7 | 0.999965, 0.999930, 0.999933 |

So level 3 sits within 4.1e-7 of level 4 on these cases, and A at nb = 5 moves by at most 4e-5 between levels. Level 4
takes about 7 min per case on one core. The synthetic control (d) gave errors of -3.1e-3 to -7.9e-3 at x* = 1.033 to
1.066 (WL disposition); the fresh x* lie inside that range.

**What I know before comparing (disclosed).** On the four old WL cases, A at nb = 5 was within 3.5e-5 of 1 and A at
nb = 4 within 3.8e-4 of A at nb = 5 (WL, reported only). The bounds below are set from those with margin. No F_l of
P, Q, R or S has been computed at any level.

**Predictions.**
- **WX-L0 (controls; run FIRST).**
  - (a) The five monomial moments of WL, exact; Phi at D1 = D2 = 1.809, zX = zH = 1, over WA6's C, is 1 to 1e-20.
  - (b) The channel sum reproduces the Slater referee of the test file to 1e-9.
  - (c) Per case, |F_l(level 4)/F_l(level 3) - 1| < 1e-5 for every l in 16..40 (measured worst 4.1e-7; margin 24).
    A case failing (c) is void.
  - (d) WL's synthetic function, (beta, gamma, delta, e) = (5, -6, 3, 0.3), at each live case's x*: amplitude within
    2e-2 of 1 at nb = 4 (measured worst 7.9e-3).
- **WX-La (the size, nothing fitted in it).** In every live case, |A - 1| < 1e-3 at nb = 5, and |A(nb = 4) -
  A(nb = 5)| < 2e-3. Falsifier: either bound broken in a live case. Probability 0.85.
- **WX-Lb (the pinch average, not the endpoint).** In every live case, |A_naive - 1| > 10 |A - 1| at nb = 5, where
  A_naive = A Phi/naive. Falsifier: fails in any live case. Probability 0.85.
- **Gate L.** GO: L0 (a), (b), (d) pass, at least three cases survive (c), La and Lb hold. BORDERLINE: La holds and Lb
  fails in one case only. STOP: otherwise.

## Part G: Gaussian orbitals, no cusp and the crossover

**Claims tested (WG's derivation and model, unchanged).** (i) Single s Gaussians: F_l > 0, the step ratio q_l =
F_(l+1)(2l+3)/(F_l T1) follows 1 + c2/l^2 + c3/l^3 + c4/l^4 + c5/l^5 with WG's coefficients on the symmetric line, and
there is no power-law tail. (ii) STO-3G and STO-6G (WG's transcriptions of the Psi4 files) follow the Slater
coefficients S_l up to a crossover l* (first l >= 6 with F_l/S_l < 1/2) placed by WG's model, and fall far faster
than l^-8 past it.

Code `debug/wx_gauss.py`: WG's `wg_common.py` series form and `wg_derive.py` model with the bond length D as a
parameter. Check before registering: at D = 1.809 it reproduces WG's model l* = 11 (E1-3) and 32 (E1-6) and its
ratios to three digits (`debug/wx_model.py E1`).

**Fresh cases.** Symmetric geometry, D the X-H length; none computed by any route.

| Case | kind | zX or alpha_X | zH or alpha_H | D |
|:--|:--|:--|:--|:--|
| g3 | single | 0.6 | 1.5 | 2.2 |
| g4 | single | 0.5 | 1.2 | 1.6 |
| E3-3, E3-6 | STO-3G, STO-6G | 1.5 | 1.0 | 2.2 |
| E4-3, E4-6 | STO-3G, STO-6G | 3.0 | 1.5 | 1.6 |

Single-Gaussian T = alpha_H^2 D^2/(alpha_X + alpha_H) and q coefficients (c2, c3, c4, c5):
g3: T = 5.185714286; -11.90862, 47.63449, 102.8869, -1729.425. g4: T = 2.168470588; -3.344037, 13.37615, -21.14847,
-56.99906.

Model crossover (`debug/wx_model.py E3 E4`, data `debug/data/wx_model_E3_E4.json`):

| Case | D sqrt(alpha_max) | model l* | model ratio at l*-1, l* | l*/(D sqrt(alpha_max)) |
|:--|:--|:--|:--|:--|
| E3-3 | 3.284 | 13 | 0.533, 0.337 | 3.96 |
| E3-6 | 10.574 | 39 | 0.546, 0.466 | 3.69 |
| E4-3 | 3.582 | 14 | 0.589, 0.384 | 3.91 |
| E4-6 | 11.536 | 43 | 0.503, 0.432 | 3.73 |

E4-6 is close to the edge: its model ratio at l = 42 is 0.503.

**What is computed for scoring.** F_l, l = 0..63, series form at 50 digits, all primitive quartets. S_l =
`F_coeffs(D, D, zX, zH, 63, level=2)`. ratio_l = F_l/S_l. Local power p_l = log(F_l/F_(l+4))/log((l+4)/l).

**Controls measured before registering** (neighbour cases = each fresh case with D + 0.05, not scored; and old sets):

| Control | Measured on | Result | Data |
|:--|:--|:--|:--|
| k-space route vs series, 50 digits, l = 0, 5, 10, 20, 30 | g3n, g4n | worst 1.3e-31, 3.0e-24 | `wx_measure_G_kspace_g3n.json`, `_g4n` |
| same, l = 0, 10, 20 | E3-3n, E4-3n | worst 4.4e-44, 6.7e-44 | `wx_measure_G_kspace_E3-3n.json`, `_E4-3n` |
| WG's old cases at 50 digits (WG rerun) | g1, g2, E1-3, E1-6 | worst 2.9e-19 (g1, l = 30) | `wg_ctrl_b_*.json` |
| Slater S_l level 2 vs 3, l = 0..63 | E1 (old), E3n, E4n | worst 1.3e-6, 5.0e-7, 1.1e-6 | `wx_measure_G_slater_*.json` |
| resummation at 104.5 deg, l = 0..60, vs direct Boys | g3n, E3-3n | 1.5e-48, 2.4e-50 | `wx_measure_G2.json` |
| q_l clause, worst residual / bound, l = 40..60; p_30 | g3n, g4n | 0.14, 0.35; 80.2, 107.3 | `wx_measure_G2.json` |
| Slater local power at l = 21, 23, 42, 50 | E3n, E4n | 7.75-7.91, 7.77-7.91 | `wx_measure_G3.json` |

All data files in `debug/data/`. The working precision of every control is stated below; none is left to the scorer.

**What I know before comparing (disclosed).** On WG's old cases (reported only), the observed l* was model l* or
model l* - 1 in all four, l*(6)/l*(3) was 3.20 and 3.25 against sqrt(alpha_max(6)/alpha_max(3)) = 3.22, and p at
model l* + 10 was 24 to 46. The windows below are set from those, tighter than WG's.

**Predictions.**
- **WX-G0 (controls; run FIRST).**
  - (a) Slater S_l: |F_coeffs level 3 / level 2 - 1| < 2e-5 for every l = 0..63, per set (measured worst 1.3e-6;
    margin 15). A set failing (a) voids its two contracted cases.
  - (b) k-space route vs series form, both at 50 working digits: relative < 1e-12, for g3 and g4 at l = 0, 5, 10, 20,
    30, and for E3-3 and E4-3 at l = 0, 10, 20 (measured worst 3.0e-24). A single case failing (b) is void; a
    contracted STO-3G case failing (b) voids its set's two cases. STO-6G is not run through (b) (cost; WG's E1-6 run
    agreed to 4.2e-48).
  - (c) Resummation, 50 digits: sum_l F_l P_l(cos 104.5 deg), l = 0..60, equals the direct Boys evaluation to a
    relative 1e-14, for g3 and E3-3.
  - (d) WG's toy q_l check (T = 1.2, l = 40: residual after l^-5 below 1e-7), re-run.
  - (e) Negative control: Slater local power p_l at l = 21, 23, 42, 50 lies in 7..9, for E3 and E4.
- **WX-Ga (no cusp, no power).** For g3 and g4: F_l > 0 for l = 0..63; for every l = 40..60,
  |q_l - (1 + c2/l^2 + c3/l^3 + c4/l^4 + c5/l^5)| < 0.5 max(|c5|, 10)/l^5; and p_30 > 50. Probability 0.95.
  Novelty: KNOWN. A check.
- **WX-Gb1 (sign).** F_l > 0 for l = 0..63 in every live contracted case. Probability 0.97. KNOWN. A check.
- **WX-Gb2 (the crossover).** In every live contracted case l*_obs exists and lies in [model l* - 2, model l* + 1]
  (E3-3: 11..14; E3-6: 37..40; E4-3: 12..15; E4-6: 41..44); and in each set with both cases live,
  2.9 <= l*_obs(6)/l*_obs(3) <= 3.5. Probability 0.7. Novelty: PROBABLY KNOWN as a mechanism (basis-set lore), not
  searched for this integral; a CHECK, not a finding.
- **WX-Gb3 (no power past the crossover).** In every live contracted case, at l0 = model l* + 10 (E3-3: 23; E3-6: 49;
  E4-3: 24; E4-6: 53): p_l0 > 16 and p_(l0+5) > p_l0. Probability 0.85. KNOWN in principle. A check.
- **Gate G.** GO: G0 (c), (d), (e) pass; g3 and g4 both live; at least three of the four contracted cases live; Ga,
  Gb1, Gb2 and Gb3 hold. BORDERLINE: as GO, but Gb2 fails. STOP: otherwise.

## Scoring rules for both parts

- The two gates are scored separately. A failure in one part does not void the other.
- Every control runs at the precision written above. A control that fails as run voids what it says it voids, whatever
  the cause. No rerun at another precision is scored.
- No scored F_l is compared with anything until every control of its part has run.

## Caps

- At most 60 min wall, at most 4 processes. No change to `geovac/`, `tests/` or `papers/` during the run.
- Part L: 1s Slater orbitals, zX = 3 and 4, zH = 1.0 to 1.5; subleading terms are fitted nuisances, not derived.
- Part G: symmetric geometry, s functions, the two published contractions. Not a claim about p or d Gaussians or
  unequal bond lengths.

## Novelty status (before comparing)

- Landau equations and Feynman-parameter pinch analysis: KNOWN (standard). The collinear-pinch size of the cusp-cusp
  term for the bond-angle series: NOT SEARCHED, treat as OPEN.
- Part G: see each item. Nothing in Part G can count as a finding without a literature search.

## Disposition

(Appended after comparing. The predictions above are not edited.)
