# Pre-registered prediction: a finite nucleus (n stops being an integer)

**Registered:** 2026-10-18, before any computation in this file (PI-directed: "you can start on that
stuff too", with the PI's intuition that the nucleus carries a skeleton of its own embedded in the
electron's, or a second projection).
**Extends:** seam sub-domain J and `docs/atlas.md`; Paper 23 `sec:rigidity` (S^3 needs exactly -Z/r),
`sec:ne` and `sec:cross-register` (the composed nuclear-electronic register at two focal lengths);
Paper 34 `sec:proj_charge_density` (the Foldy/Friar projection, radius consumed); the multi-focal
composition wall (`memory/multi_focal_wall_pattern.md`; Paper 34 calls the Breit term its cleanest
instance).
**Commit status:** committed on its own, before any computation, with `2026-10-18_zero_energy.md`.

## The door

A point nucleus gives exactly -Z/r, and Paper 23's rigidity theorem says the S^3 projection needs
exactly -Z/r. Give the nucleus a size and SO(4) must break. The question is what breaks, what
survives, and what the price is. Model: a uniformly charged sphere of radius R, potential
-Z(3R^2 - r^2)/(2R^3) inside and -Z/r outside. Atomic units throughout; the nuclear mass is infinite
except in the control.

## Where the PI's intuition sits in the corpus (disclosed, not a prediction)

The nucleus as a second skeleton inside the electron's is already built, as an architecture: Paper 23
`sec:ne` composes a nuclear register (oscillator basis, Paper 24's skeleton) with an electron register
(the (n, l, m) labels), "at categorically different focal lengths", and `sec:cross-register` couples
them by a two-focal-length Coulomb operator. The multi-focal wall says the labels of two such
skeletons couple cleanly but their spatial focal lengths do not compose: the scale ratio is consumed.
Paper 34 catalogues the consumed number for this case, the charge radius. What this door adds is the
accounting on the electron's skeleton: what the second scale releases there.

A formal match noted before computing, recorded and not read: inside a uniform sphere the electron
moves in an isotropic oscillator potential, the operator of Paper 24. The interior solution is a
Kummer function of r^2 and the exterior one a Coulomb function, so the whole problem is one matching
of an oscillator solution to a Coulomb solution at r = R. That oscillator is the electron's, not the
nucleons' (Paper 23's shell model is a different oscillator); any reading gets the curve-fit audit.

## What I know before computing (disclosed)

Textbook: the finite-size shift Delta E = (2 pi / 3) Z |psi(0)|^2 <r^2> for s states (Foldy; Friar
1979); quantum-defect theory and Seaton's theorem (the defect at threshold equals the continuum phase
shift divided by pi; Seaton 1958, 1983). Worked by hand before computing: first-order perturbation
with R_nl(r) ~ A_nl r^l near the origin and the J13 product A_nl^2 dn/dE = 4Z (2Z)^(2l)
prod_{s=1..l}(1 - s^2/n^2) / ((2l+1)!)^2.

## FN0 (control: finite nuclear mass releases nothing)

With nuclear mass M the levels are -mu Z^2/(2 n^2), mu = M/(M+1): one rescaling of the focal length,
SO(4) exact, every defect zero.

## FN1 (what the release leaves standing)

l, m and the Gaunt rules survive (the nucleus is round). n survives as a label: the number of radial
nodes n - l - 1 is unchanged for ZR small. What is released is n being a whole number, and with it
the l-degeneracy within a shell: E_nl = -Z^2/(2 (n - delta_nl)^2).

## FN2 (main: the defect)

To leading order in ZR,

  delta_nl = -c_l (Z R)^(2l+2) prod_{s=1..l}(1 - s^2/n^2),
  c_l = 4 * 2^(2l) k_l / ((2l+1)!)^2,  k_l = 1/(2l+2) - 3/(2(2l+3)) + 1/(2(2l+5)),

so c_0 = 2/5 and c_1 = 1/105: the s-wave defect is the same for every n, -(2/5)(ZR)^2 =
-(2/3) Z^2 <r^2>, and the p-wave defect carries the factor (1 - 1/n^2). The relative correction to
the s-wave defect at the next order is of order ZR, with an n-independent coefficient (from the cusp
of the wavefunction at the origin).

## FN3 (continuation across threshold)

Since 1 - s^2/n^2 = 1 + 2 s^2 E / Z^2, the leading defect is a polynomial in E, and it continues to
positive energy as the short-range phase shift: the exact phase shift Delta_l(E), relative to
Coulomb, satisfies Delta_l(E)/pi = -c_l (ZR)^(2l+2) prod (1 + 2 s^2 E / Z^2) at leading order, and,
exactly (Seaton), lim_{n -> infinity} delta_nl = lim_{E -> 0+} Delta_l(E)/pi.

## FN4 (the price)

The release consumes one number per charge moment of the nucleus, here R. The leading coefficients
(c_l, the product over s) are rational. Whether a log(ZR) enters at higher order is measured, not
predicted. Predicted verdict: FREE for the price, the first released-axis door whose price is
calibration data rather than a number class fixed by the geometry; FORCED for what survives and for
the rational structure of the leading defect. The pi-free structure is predicted: no pi enters the
defect at leading order (the 1/(4 pi) of |psi(0)|^2 cancels against the 4 pi of the volume integral).

## Tests and decision rules

1. Exact bound energies by matching the interior Kummer solution to the exterior Whittaker function
   W_{nu, l+1/2}(2 Z r / nu) at r = R, at 40 digits, for l = 0, 1, 2, n up to at least 40, ZR in
   {1e-3, 3e-3, 1e-2, 3e-2}; delta_nl = n - nu. Fit delta against powers of ZR.
2. Exact phase shifts by matching to F_l and G_l at r = R, at five energies, the same radii.
3. Seaton's limit checked from both sides; FN3's leading form checked at each energy.
4. A second route for FN2: first-order perturbation by direct quadrature of the potential difference
   against the hydrogen radial function.
5. FN0 symbolic. Bounded searches state their bounds.

## Falsifiers

1. FN1: a node count that changes at small ZR; an l missing.
2. FN2: c_0 != 2/5 or c_1 != 1/105; an n-dependent leading s-wave defect; a p-wave defect without
   the factor (1 - 1/n^2).
3. FN3: Seaton's limit failing; a leading phase shift not given by the continued polynomial.
4. FN4: an irrational leading coefficient; a pi in the leading defect.

## Honest caps, stated before the result

- One nuclear model (the uniform sphere); real nuclei have other profiles, which change the
  higher-order coefficients but not the leading term's form (the leading term depends on <r^2> only).
- Non-relativistic, one electron. The relativistic case, where the finite nucleus becomes compulsory
  near Z alpha = 1, is not touched.
- The physics is textbook (Foldy; Friar 1979; Seaton). The contribution, if any, is the accounting
  on the electron's skeleton, and the link to J13's product.

## Home, if confirmed

Proposed, PI decides: an atlas door and a seam J row (n released from integrality by a second
scale), with its home in Paper 23 beside the rigidity theorem, or in Paper 34 beside the Foldy/Friar
projection; a guard in `tests/`.

## Amendment (2026-10-18, after commit bd2435d1, before any computation)

Noticed while writing the computation, not registered: first-order perturbation with the J13 result
gives tan Delta_l ~ pi mu(E) / (1 - e^(-2 pi eta)), where mu(E) is the continued leading defect, because
the energy-normalized continuum density is the continued bound density times that factor. So FN3's
registered leading form, Delta_l/pi = mu(E) with no factor, is expected to fail at energies where
e^(-2 pi eta) is not small (about 4 percent at eta = 1/2), while Seaton's limit at threshold is
unaffected (the factor tends to 1). FN3 stays as registered and is scored as registered; this note
records the expectation before the test.
