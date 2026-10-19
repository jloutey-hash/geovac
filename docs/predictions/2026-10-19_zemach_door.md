# Pre-registered prediction: the hyperfine door, a finite magnetization distribution (Zemach)

**Registered:** 2026-10-19, before any computation in this file (PI-directed: "Both please", after the
suggestion to re-price the multi-focal wall's hyperfine instance as a door, the magnetic twin of J14).
**Extends:** `2026-10-18_finite_nucleus.md`, `2026-10-18_nuclear_profiles.md` (seam row J14, atlas door
24); Paper 23 `sec:finite_nucleus`; the multi-focal wall's instance HF-4 (hyperfine coupling point-like in
the proton position, no magnetization distribution).
**Commit status:** committed on its own, before any computation.

## The door

The hyperfine interaction of an s electron with a point nuclear moment is the Fermi contact term, set by
|psi(0)|^2. A real nucleus spreads its magnetization over a second length scale, and its charge over a
third. The door releases the point magnetization. Question: is the price again a list of consumed
moments, one per order, and how do two distributions (charge and magnetization) compose?

## What I know before computing (disclosed)

From memory, not consulted: Zemach (1956), the leading correction -2 Z m alpha r_Z with r_Z the first
moment of the convolution of the charge and magnetization densities; Bohr and Weisskopf (1950); Friar
(1979); higher-order size corrections to hyperfine structure (Friar and Payne, 2005; Karshenboim), where a
logarithm at second order is probably known. Derived by hand while drafting: (a) for a spherical
magnetization density the dipolar field averages to zero over an s electron and a uniformly magnetized
shell has no interior field, so the hyperfine factor is the overlap below; (b) the short-distance Taylor
series of point-Coulomb s densities, |psi_n(r)|^2 / |psi_n(0)|^2 = 1 - 2Zr + Z^2 (5/3 + 1/(3 n^2)) r^2 + ...;
(c) the log coefficient of ZM3, from the irregular zero-energy solution 1/r - 2Z ln r, cross-checked
against J14's verified -delta_0^2 ln term (the same Green-function piece reproduces it). Nothing was
computed.

## Model

Non-relativistic, s states, infinite nuclear mass, spin magnetization only, Z = 1 with x = ZR. The
hyperfine factor of level n is H_n = integral of rho_M |psi_n|^2 d^3r divided by |psi_n^point(0)|^2,
psi_n normalized; delta_n = H_n - 1. Charge: point, or a uniform sphere of radius R. Magnetization
(normalized): point, a uniform ball of radius R, or a thin shell of radius R. Zemach moments
r_Z = mean distance between a charge point and a magnetization point: uniform charge with point
magnetization 3R/4; uniform charge with uniform ball 36R/35; uniform charge with shell 6R/5; point charge
with uniform ball 3R/4.

## Predictions

- **ZM1 (leading, Zemach):** delta_n = -2 Z r_Z + O(x^2 ln x), the same for every n (the Kato cusp fixes
  the linear term of every s density). Coefficients of x: -3/2 (uniform charge, point magnetization);
  -72/35 (uniform ball, uniform ball); -12/5 (uniform charge, shell magnetization). The two
  distributions compose at this order through one convolution moment and nothing else.
- **ZM2 (a point charge brings no logarithm):** with a point charge, delta_n is an entire function of x;
  its Taylor coefficients are polynomials in 1/n^2 with rational coefficients for rational-moment
  profiles. For the uniform-ball magnetization: delta_n = -(3/2) x + (1 + 1/(5 n^2)) x^2 + O(x^3).
  Checked exactly for n = 1..4 through x^6.
- **ZM3 (the composition logarithm):** with a uniform-sphere charge, delta_n contains x^2 ln x with
  coefficient (4/3) Z^2 <r^2>_E / R^2 = 4/5, which is -2 delta_0 with delta_0 = -(2/3) Z^2 <r^2>_E, J14's
  leading quantum defect; the same for every magnetization profile and every n. So the logarithm comes
  only from the charge distribution, and it is J14's composition piece reappearing one order earlier.
- **ZM4 (number class):** through the orders computed, every coefficient for these profiles is rational
  except the x^2 ln x term and its companion constant, which is weight one (ln 2, gamma_E, and possibly
  harmonic numbers of n): the class of door 11. The companion constant's form is not registered.

## The price, predicted

FREE: the magnetization distribution is consumed, through its convolution with the charge distribution,
one moment per order. FORCED: the labels; the n-independence of the leading term; the logarithm's
coefficient, fixed by J14's delta_0. If this holds, the hyperfine instance of the multi-focal wall is a
door with a price list, like J14, rather than a wall. Changing its status is a /walls question for the PI.

## Tests and decision rules

1. ZM2: sympy, point-Coulomb psi_n exactly, n = 1..4, series in x through x^6.
2. ZM1 and ZM3: 40 digits. Interior: exact power series for the uniform-sphere potential. Exterior: the
   Whittaker function W_(nu, 1/2)(2r/nu). Energy from matching log-derivatives at r = R. Normalization by
   quadrature. x = 1e-5 .. 3e-3, 12 log-spaced points; fit columns x, x^2 ln x, x^2, x^3 ln x, x^3, x^4.
   n = 1 and n = 2, each of the three magnetization profiles with the uniform-sphere charge.
3. Gates: the x coefficient within 1e-4 relative of -2 r_Z / R; the x^2 ln x coefficient within 2e-2
   relative of 4/5. Fit stability reported by dropping the last column.

## Falsifiers

1. ZM1: an x coefficient off its Zemach value, or depending on n.
2. ZM2: a logarithm, or a coefficient that is not a polynomial in 1/n^2.
3. ZM3: no x^2 ln x term, or a coefficient other than 4/5, or one depending on n or on the magnetization
   profile.

## Honest caps

- Textbook in substance (Zemach; Bohr and Weisskopf; Friar and Payne). The contribution, if any, is the
  accounting on the door, and the tie between the hyperfine logarithm and J14's quantum defect.
- Non-relativistic; no recoil, no QED, no orbital (Bohr-Weisskopf) currents; model profiles.

## Result (2026-10-19, same session, after commit be16fccd)

Script `debug/zemach_door.py`; data `debug/data/zemach_door.json`. ZM1, ZM2 and ZM3 hold. ZM4 fails in
part: an unregistered x^3 ln x term is present.

Code defect fixed before scoring: the first interior series used r^2/R^2 for the r^2 term of 2V
instead of r^2/R^3. A smoke test caught it: the quantum defect came out near -(2/3) x^3 instead of
J14's -(2/5) x^2. After the fix the defect at x = 1e-3 is -3.99491399e-7, matching J14's derived series
-(2/5) x^2 + (32/63) x^3 to the printed digits; 40 and 50 digits agree to all printed figures.

- **Zemach moments (sympy, exact):** 3/4, 36/35, 6/5 for point, ball and shell magnetization with a
  uniform charge, as registered.
- **ZM2 (exact, n = 1..6, through x^6):** no logarithm at any order; the x coefficient is -3/2 and the x^2
  coefficient 1 + 1/(5 n^2) for every n; each coefficient c_k(n) is a polynomial in 1/n^2 of degree
  floor(k/2), fitted on the first floor(k/2) + 1 values of n and reproducing all six.
- **ZM1 (40 digits, x = 1e-5 .. 3e-3, 12 points, n = 1, 2, 3):** the x coefficient is -1.5, -2.0571428571,
  -2.4 for point, ball and shell, relative error 1e-10 or better with the registered columns, the same
  for every n (gate 1e-4).
- **ZM3:** the x^2 ln x coefficient is 0.800000 for all nine cases (three profiles, three n), relative
  error 4e-6 or better with the registered columns and 5e-7 with x^4 ln x added (gate 2e-2). A fit with
  no logarithm leaves a residual 1.8e-8 against 2e-15. So the hyperfine logarithm is -2 delta_0, set by
  the charge distribution alone, the same for every magnetization profile and every n.
- **ZM4, failed in part:** the fit also requires x^3 ln x: dropping that column raises the residual from
  2e-15 to 6e-12. Its coefficient is the same for every n and depends on the profile; it equals
  (4/5) times the Zemach coefficient plus -1.0147, the last the same for all three profiles to 2e-5
  between fits of one column set and to 1e-2 across column sets, so not identified. The registration
  said that only the x^2 ln x term and its companion constant were non-rational, so it is falsified.
  The companion's class is consistent with weight one: after removing ZM2's magnetization part, the
  x^2 coefficient is (4/5) [psi(n) - ln n - 1/n] plus a constant that depends on the profile and not on
  n (residual spread 1e-6 over n = 1, 2, 3). That form is post-hoc: read off n = 1, 2, written down
  before the n = 3 run, then confirmed by it (prediction 0.8544435, measured 0.85444352). The profile
  constants (1.61973, 1.60095, 1.626 for point, ball, shell) are not identified.

What this says, in one line: the magnetization door is a price list, like J14. Two distributions
compose at leading order through one convolution moment, the Zemach radius, the same for every level.
The logarithm comes only from the charge distribution, with coefficient -2 delta_0 from J14's quantum
defect. A point charge brings no logarithm at any order. The n-dependence enters at second order in
weight-one form. Predicted verdict FREE (the distribution consumed) with a forced shape; whether the
hyperfine instance HF-4 leaves the multi-focal wall is a /walls question for the PI.

**Disposition (2026-10-19, PI-delegated: "You choose").** HF-4 stays an instance of the multi-focal
composition wall, priced like J14 (walls register note of this date; Paper 57 `sec:unity`). Home: Paper 23
`sec:finite_nucleus`, "A finite magnetization". Guard `tests/test_paper23_hyperfine_door.py`, fire-tested 8/8.
Seam J16 (C), atlas door 26.
