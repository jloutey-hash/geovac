# Pre-registered prediction: the size of the bond-angle singularity at unequal bond lengths, from the Landau pinch (WL)

**Registered:** 2026-10-28, before any F_l in this file is computed. PI-directed ("Yeah, let's do this", to a Landau
map of the integral's singular points, with unequal-length water as the first test).
**Commit status:** committed on its own, before any computation it scores.
**Rests on:** WA6 (`2026-10-27_water_cusp_next_term.md`, GO) and its machinery `debug/wa6_local.py`; the channel formula
`F_coeffs` of `tests/test_paper59_water_bond_angle.py`.

## The question

At unequal bond lengths D1 != D2 the singular point of the bond-angle series is x* = (D1^2 + D2^2)/(2 D1 D2) > 1,
where |Y - Z|^2 = 0 with Y and Z on a complex bond angle. WA5 found that the cusp-cusp size C of the equal-length case
does not carry over, because every gradient term enters at the same order there. That size has been the open cap of
Paper 59's water passage. Can it be derived, with nothing fitted?

## The derivation (done before registering)

- The cusp-cusp singular part is S ~ Phi |w|^7 near the null cone, w = Y - Z. At a null w, a Taylor term of degree
  (d1, d2) in the two cusp functions G1 (about Y) and G2 (about Z) contributes to |w|^7 with coefficient
  (pi^2/630) (-1)^d1 (d1+1)! (d2+1)! 24/(d1+d2+4)!. Derived from the Fourier route of WA6 (most singular part of each
  inverse transform) and checked exactly against `debug/wa6_local.py` for five monomials: (a_z, 1) -2/5; (1, b_z)
  2/5; (a_x, b_x) -2/15; (a_z^2, 1) 1/5; (a_z^2, b_z) 2/35.
- These are the moments of 24 s u on the simplex s, u >= 0, s + u <= 1. So
  Phi = (pi^2/630) 24 int int s u G1(-s w) G2(u w) ds du:
  the cusp functions averaged over the segment between the two H atoms. This is the collinear pinch of the Landau
  equations for the three "propagators" |a|, |b|, 1/|w + a - b|: all three vanish together when a and b lie along
  the null w.
- On the null cone |Y + s(Z - Y)|^2 = (1 - s) D1^2 + s D2^2, real and positive, and sinh(zH |a|)/|a| = zH. So
  Phi = (pi^2/630) 24 zH^2 NH^2 NX^2 int int s u exp(-zX (R1(s) + R2(u))) ds du,
  R1(s) = sqrt((1 - s) D1^2 + s D2^2), R2(u) = sqrt((1 - u) D2^2 + u D1^2).
  At D1 = D2 it equals WA6's C exactly (ratio 1.000...). Script `debug/wl_local.py`.
- Leading term: L_l = the Legendre coefficients of Phi (2 D1 D2)^(7/2) (x* - x)^(7/2), computed by 40-digit
  quadrature.
- "Naive" size, for contrast: C with e^{-zX (D1 + D2)}, the endpoint value s = u = 0. WA5 used it.

## Cases (fresh; none has been computed)

zX = 4 is chosen, with an outlying D2, so that Phi and the naive size differ by 8-52 percent. D1 = 1.809, 1s
Slater orbitals, zH the H exponent.

| Case | zX | zH | D2 | rho = D</D> | Phi | naive | Phi/naive |
|:--|:--|:--|:--|:--|:--|:--|:--|
| A | 4 | 1.24 | 1.4 | 0.774 | 8.70813295609469e-7 | 7.93129092961136e-7 | 1.0979 |
| B | 4 | 1.24 | 2.2 | 0.822 | 3.57535903642545e-8 | 3.23296898684584e-8 | 1.1059 |
| C | 4 | 1.24 | 2.6 | 0.696 | 9.91259923516535e-9 | 6.52725181228884e-9 | 1.5186 |
| D | 2.25 | 1.24 | 2.6 | 0.696 | 2.82650368571769e-6 | 2.60618713474370e-6 | 1.0845 |

## The extraction (fixed before computing)

The subleading terms are not derived. They are (x* - x)^(9/2), (x* - x)^(11/2), ..., from the next contractions and from
Phi's variation off the cone. So the leading size is extracted from F_l by linear least squares onto the exact Legendre
coefficients of (x* - x)^(7/2 + j), j = 0..nb-1, over l = 16, 18, ..., 40, each row divided by the j = 0 coefficient.
A = the j = 0 coefficient found, divided by Phi (2 D1 D2)^(7/2). Reported for nb = 4 and nb = 5; scored on nb = 4.
The subleading coefficients are fitted nuisances. The leading size is not fitted: it is compared with Phi.

## What I know before comparing (disclosed)

- No F_l for cases A-D has been computed. F_l for WA5's unequal cases (D2 = 2.5, 1.0; zX = 1 and 2.25) are stored, and
  their ratios to the naive size at two l are known (WA5 disposition). They are not used here, and are reported only
  after scoring.
- The extraction was tried on synthetic functions only (`debug/wl_extrap.py` and a scratch test):
  (x* - x)^(7/2) (1 + beta (x* - x) + gamma (x* - x)^2 + delta (x* - x)^4) + e (e^{2x}) at x* of D2 = 1.4 and 2.6, with
  (beta, gamma, delta, e) up to (5, -6, 3, 0.3) and (-4, 6, -5, 1): amplitude error at nb = 4 at most 1.3e-2, at nb = 5
  at most 1e-5. A plain extrapolation in 1/l through four points was tried first and rejected (errors up to 7e-2).
- So the extraction's own error is about 1e-2 when the remainder outside the basis is as large as in the synthetic case.
  How large water's remainder is, is not known.

## Predictions

**WL-0 (controls; run FIRST).**
- (a) The five monomial moments above, re-run in the scorer, are exact; Phi at D1 = D2 = 1.809 (E1) over C is 1 to
  1e-20.
- (b) The channel sum reproduces the Slater referee (the `F_coeffs` referee check of the test file) to 1e-9.
- (c) Per case, |F_l(level 2)/F_l(level 1) - 1| < 1e-7 for every l in 16..40. A case failing (c) is void.
- (d) The synthetic extraction with (beta, gamma, delta, e) = (5, -6, 3, 0.3) at each live case's x* gives an
  amplitude within 2e-2 of 1 at nb = 4.

**WLa (the size, nothing fitted in it).** In every non-void case, |A - 1| < 2e-2 at nb = 4, and the nb = 4 and nb = 5
values of A agree within 1e-2.
- Falsifier: either bound broken in any non-void case.
- Probability 0.6.

**WLb (it is the pinch average, not the endpoint).** In every non-void case, |A_naive - 1| > 3 |A - 1|, where A_naive
is the extracted amplitude divided by the naive size times (2 D1 D2)^(7/2) (so A_naive = A Phi/naive).
- Falsifier: fails in any non-void case.
- Probability 0.6.

## Decision gate

- **GO:** WL-0 (a), (b), (d) pass, at least three cases survive (c), and WLa and WLb hold. Then the size of the
  bond-angle singularity at unequal bond lengths is the Landau collinear-pinch average of the two cusp functions along
  the segment between the H atoms, with nothing fitted in it.
- **BORDERLINE:** WLa holds and WLb fails in one case only.
- **STOP:** any of WL-0 (a), (b), (d) fails, fewer than three cases survive (c), or WLa fails.

## Caps

- At most 45 min wall; at most 4 processes. No change to `geovac/` or `tests/`.
- 1s Slater orbitals, zX = 2.25 and 4, one H exponent. The subleading terms are fitted nuisances, not derived.
- "Landau pinch" names the mechanism of the derivation; this record does not test the general Landau map.

## Novelty status (before comparing)

- Landau equations and Feynman-parameter pinch analysis: KNOWN (particle physics, standard).
- The collinear-pinch size of the cusp-cusp term, and its use for the bond-angle series: NOT SEARCHED. Treat as OPEN.

## Disposition

(Appended after comparing. The predictions above are not edited.)
