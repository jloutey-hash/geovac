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

## Result (2026-10-18, same session, after commits bd2435d1 and 5ba129aa)

Script `debug/finite_nucleus_door.py`; data `debug/data/finite_nucleus_door.json`. Every quantity is one
matching of the interior Kummer solution to the exterior Whittaker or Coulomb function at r = R, at
50-60 digits. FN0, FN1, FN2 and FN4 hold; FN3 holds at threshold and fails as registered away from it,
as the amendment expected; the amended form holds.

- **FN0: CONFIRMED [SYMBOLIC].** The reduced mass only rescales the focal length; every defect is zero.
- **FN1: CONFIRMED.** At ZR = 0.01 the radial node count equals n - l - 1 for all eight states tested
  (1s, 2s, 3s, 5s, 2p, 3p, 4p, 3d). l and m survive by construction (a round nucleus); the defects
  below depend on l, so the l-degeneracy within a shell is gone. The defect depends on Z and R only
  through ZR, exactly (Z = 2 at R/2 reproduces Z = 1 at R to all digits).
- **FN2: CONFIRMED.** Over ZR = 1e-3 .. 3.2e-2 and n up to 40, delta_nl divided by the registered
  leading form tends to 1: to 1e-10 for l = 0 once log terms are allowed (below), and to 1e-11 for
  l = 1, 2. So c_0 = 2/5 and c_1 = 1/105 hold, the s-wave defect is the same for every n, and the
  factor prod(1 - s^2/n^2) is present for l = 1 and 2 (c_2 = 1/28350 from the same formula). The next
  order is n-independent, as registered: for l = 0 the relative coefficient of ZR is -1.2698412695 at
  n = 1 and at n = 5, which is -80/63 to 4e-10 (an identification by agreement, not a PSLQ
  certification); first-order perturbation alone gives -5/6, the rest is second order. Second route:
  first-order perturbation by direct quadrature reproduces the leading form with its -5/6 ZR
  correction for 2s, 5s, 2p, 5p.
- **FN3: Seaton's limit CONFIRMED; the registered continuum form FAILS away from threshold;
  the amended form HOLDS.** At ZR = 0.01 the bound defects extrapolated to n -> infinity and the
  phase shifts extrapolated to k -> 0 agree to 1e-15 (l = 0) and 4e-15 (l = 1). Away from threshold,
  at ZR = 1e-3, the phase shift over pi divided by the continued leading defect is 0.99873, 1.00060
  and 1.04384 at k = 1/2, 1, 2 for l = 0 (a 4.4 percent miss at k = 2, eta = 1/2), so the registered
  form is falsified. With J13's factor included, tan Delta_l / pi = mu(E) / (1 - e^(-2 pi eta)), the
  ratio is 0.9987284 at all three energies, equal to the bound side's ratio at the same ZR (the next
  order, -80/63 ZR, continues too), and likewise 0.99935 for l = 1. So the defect continues across
  threshold as one analytic function of E, and the continuum adds only J13's factor.
- **FN4: CONFIRMED, and the log question answered.** The leading coefficients are rational
  (2/5, 1/105, 1/28350, and the product over s), and no pi enters the leading defect. Measured, not
  registered: a log does enter. For the s-wave, a fit in powers of ZR alone stalls at a residual of
  1.7e-8, while adding (ZR)^k log(ZR) columns for k = 2, 3, 4 drops it to 6e-16; the coefficient of
  (ZR)^2 log(ZR), relative to the leading term, is 0.400002, consistent with 2/5. So
  delta_s = -(2/5) x^2 [1 - (80/63) x + (2/5) x^2 log x + ...], x = ZR. The price is the consumed
  radius R; the coefficients multiplying its powers are rational, with logs of R from second order on.

Defects in my own test code, caught before any verdict: the root finder started at n and failed on
the steep matching function (fixed by starting the secant at the leading-order estimate); the first
node count counted a convention sign flip of the exterior function at r = R as a node, giving 2s and
3p one node too many (fixed by matching the exterior sign to the interior at R).

What the door says, in one line: a finite nucleus keeps l, m, Gaunt and the node label n, releases n
from being a whole number, and the defect it opens is a consumed radius times rational coefficients
(with logs from second order); the defect continues across threshold as one function of energy, and
the continuum adds only J13's factor. The interior is the electron's oscillator matched to its
Coulomb exterior at the nuclear radius; that formal match with Paper 24 is recorded, not read.

Home (2026-10-18, PI-approved): Paper 23 `sec:finite_nucleus`, beside the rigidity theorem, with a pointer from Paper 34 `sec:proj_charge_density`; seam row J14 (C); atlas door 24 (FREE for the radius); guard `tests/test_paper23_finite_nucleus_door.py`, written as its own step and fire-tested 12/12. While fire-testing, the guard's node check was found not to sample the interior (its first grid point lay outside the nucleus); a point inside was added, after which the sign-matching plant fires.

## Owed item: the measured coefficients derived (2026-10-18, after v5.20.26)

Script `debug/finite_nucleus_threshold_series.py`. Seaton's limit lets the s-wave defect be computed at
E = 0: inside, with s = r/R, the radial equation is u'' = -x(3 - s^2)u, a power series in x = ZR;
outside, the energy-normalized zero-energy pair is sqrt(2r) J_1(sqrt(8r)) and -sqrt(2r) Y_1(sqrt(8r))
(Wronskian checked). Matching at r = R gives, exactly to the order kept,

  mu = -(2/5) x^2 + (32/63) x^3 - x^4 [ (4/25) ln(2x) + (8/25) gamma_E + 28477/225225 ]
       + x^5 [ (128/315) ln(2x) + (256/315) gamma_E - 30486592/57432375 ] + O(x^6 ln x).

- The measured -80/63 is derived: (32/63)/(-2/5). The measured 2/5 on x^2 ln x is derived: (4/25)/(2/5).
- The x^3 term is the Friar term (1/3) Z^3 <r^3>_(2): for the uniform sphere the self-convolution
  moment is exactly 32 R^3/21, computed symbolically from the overlap volume of two spheres. So the
  next order consumes a second, distinct number of the nucleus, not a power of R.
- New: at fourth order the weight-one constants gamma_E and ln 2 enter, the same class as the
  two-centre seed set of row J1, and still no pi. The x^4 term has the form
  -delta_0^2 [ln(2x) + 2 gamma_E] + rational, delta_0 = -(2/5) x^2.
- Check against the exact matching: at x = 1e-2, 3.3e-3, 1e-3 the threshold defect minus the series is
  1.6e-12, 2.1e-15, 1.3e-18, that is about 1.3 to 1.6 times x^6, slowly drifting as an x^6 ln x term
  would.
