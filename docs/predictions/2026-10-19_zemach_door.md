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
