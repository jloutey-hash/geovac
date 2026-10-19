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
