# Pre-registered prediction: the finite-nucleus defect for other charge profiles

**Registered:** 2026-10-18, before any computation in this file (PI-directed: "we can start on the
owed items"; the first owed item of atlas door 24, "other nuclear profiles beyond the leading term").
**Extends:** `2026-10-18_finite_nucleus.md` (the uniform sphere) and its derivation section (the series
through x^5, script `debug/finite_nucleus_threshold_series.py`); seam row J14; Paper 23
`sec:finite_nucleus`.
**Commit status:** committed on its own, before any computation for these profiles.

## The question

For the uniform sphere the s-wave defect is -(2/5)x^2 + (32/63)x^3 - x^4[(4/25)ln(2x) + (8/25)gamma_E
+ rational] + ..., and the x^3 term is the Friar term. Does that structure belong to the sphere, or to
any round nucleus? If it is general, the price of a finite nucleus is a list of consumed moments, one
per order, and the weight-one constants at fourth order are universal.

## Profiles

- **Thin shell:** all charge on a sphere of radius R; x = ZR; <r^2> = R^2; the self-convolution of
  two uniform shells has <d^3> = 16 R^3/5 (distance density d/(2R^2) on [0, 2R]). Inside, the
  potential is the constant -Z/R. Exactly solvable.
- **Gaussian:** charge density proportional to exp(-r^2/(2b^2)); y = Zb; <r^2> = 3b^2; the
  self-convolution is a Gaussian of variance 2b^2 per axis, so <d^3> = 32 b^3/sqrt(pi). The potential
  is -Z erf(r/(sqrt(2) b))/r. Solved numerically.

## Predictions

- **NP1 (leading):** delta_s = -(2/3) Z^2 <r^2> for every profile: -(2/3)x^2 (shell), -2y^2 (Gaussian).
- **NP2 (next order, the Friar term):** + (1/3) Z^3 <r^3>_(2): + (16/15)x^3 (shell),
  + 32/(3 sqrt(pi)) y^3 (Gaussian). So for the Gaussian the defect carries pi^(-1/2) at third order,
  a number class inherited from the profile, not from the electron's skeleton.
- **NP3 (fourth order, the universal part):** the coefficients of ln(Z) and of gamma_E at fourth order
  are -delta_0^2 and -2 delta_0^2, with delta_0 the leading term, for every profile; the rest of the
  fourth-order term depends on the profile. For the shell the whole fourth-order term is
  -delta_0^2 [ln(2x) + 2 gamma_E] plus a rational multiple of x^4. Confidence moderate: it rests on
  the log coming from the exterior irregular function times the square of the leading defect.

## Tests and decision rules

1. Shell: the same exact series derivation as the uniform sphere's (interior u = sin(sqrt(2x) s)),
   through x^5; then an exact numerical check by matching at E = 0 at 50 digits for x = 1e-2, 3e-3,
   1e-3.
2. Gaussian: integrate the E = 0 radial equation from the origin to r = 12 b at 30 digits, match to
   the exterior zero-energy pair; y = 1e-3 .. 3.2e-2; fit delta/(-2y^2) in powers of y with and
   without y^2 ln y; the y coefficient must equal -16/(3 sqrt(pi)) to 1e-3 relative, and the y^2 ln y
   coefficient must equal 2 (that is, -delta_0^2 / delta_0 with delta_0 = -2y^2) to 2e-2 relative.
3. Bounded fits state their ranges.

## Falsifiers

1. NP1: a leading term not set by <r^2> alone.
2. NP2: a third-order term other than the Friar term for either profile.
3. NP3: an ln or gamma_E coefficient at fourth order other than -delta_0^2, -2 delta_0^2.

## Honest caps

- Textbook in substance: the Friar term is Friar's (1979); the fourth-order logarithm of
  finite-size corrections is known in the precision literature. The contribution, if any, is the
  accounting: which orders consume which moments, and which constants are universal.
- Non-relativistic, s-waves only, infinite nuclear mass.

## Result (2026-10-18, same session, after commit 8292bf51)

Script `debug/nuclear_profiles_door.py`; data `debug/data/nuclear_profiles_door.json`. All three
registered predictions hold, for both profiles.

- **Shell, exact series (sympy, through x^5):**
  mu = -(2/3) x^2 + (16/15) x^3 - x^4 [ (4/9) ln(2x) + (8/9) gamma_E + 127/315 ]
       + x^5 [ (64/45) ln(2x) + (128/45) gamma_E - 4448/2835 ] + ...
  NP1 (-2/3) and NP2 (16/15, the Friar term with <d^3> = 16 R^3/5) hold symbolically; NP3 holds
  symbolically: the fourth-order term minus -delta_0^2 [ln(2x) + 2 gamma_E] is the rational -127/315 x^4.
  Exact matching at 50 digits differs from the series by 6.9e-12, 9.0e-15, 5.4e-18 at x = 1e-2, 3.3e-3,
  1e-3, about 5 to 7 times x^6, slowly drifting as an x^6 ln x term would.
- **Gaussian, numerical (odefun at 30 digits to r = 12 b; y = 1e-3 .. 3.2e-2, 11 points):** with log
  columns the fit residual is 3e-11 (without them 9e-7). NP1: the leading ratio is 1 to 1.2e-9. NP2:
  the y coefficient is -3.009016002 against -16/(3 sqrt(pi)) = -3.009011112, relative error 1.6e-6
  (gate 1e-3). NP3: the y^2 ln y coefficient is 1.99545 against 2, relative error 2.3e-3 (gate 2e-2).
  So for a Gaussian nucleus pi^(-1/2) enters the defect at third order, from the profile.

Not registered, derived for both the uniform sphere and the shell: at fourth and fifth order the
weight-one constants appear only in the combination ln(2ZR) + 2 gamma_E, the expansion of the
exterior irregular function Y_1. A measured observation across two exactly solved profiles, not a
theorem.

What this says, in one line: the structure is not the sphere's: it holds for all three round profiles tested. The
leading order consumes <r^2>, the next the Friar moment <r^3>_(2), and these carry whatever number
class the profile carries (rational for the sphere and the shell, pi^(-1/2) for the Gaussian); from
fourth order on, a universal piece -delta_0^2 [ln(2ZR) + 2 gamma_E] enters for every profile tested (the Gaussian's gamma_E coefficient was not separated), the
weight-one class of row J1, with the rest profile-dependent.

Home (2026-10-18, PI-approved): Paper 23 `sec:finite_nucleus` (eq:defect_s now the derived series; an "Other profiles" paragraph; the price paragraph rewritten as one consumed moment per order); seam J14 and atlas door 24 updated; J14 recorded as the seventh instance of the multi-focal composition wall (Paper 57 `sec:unity`, PI decision); guard `tests/test_paper23_nuclear_profiles.py`, written as its own step and fire-tested 9/9.
