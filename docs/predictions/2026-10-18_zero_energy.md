# Pre-registered prediction: the zero-energy door (n released to infinity)

**Registered:** 2026-10-18, before any computation in this file (PI-directed: "You can start #1", after
the spin and continuum doors, `2026-10-18_large_spin.md` and `2026-10-18_hydrogen_continuum.md`).
**Extends:** seam row J13 and atlas door 23 (the continuum); Paper 7 `sec:continuum`.
**Commit status:** committed on its own, before any computation, with `2026-10-18_finite_nucleus.md`.

## The door

Bander and Itzykson's projection has three shapes: the sphere S^3 for E < 0, the hyperboloid H^3 for
E > 0, and flat space at E = 0. The flat case is the third member. It is reached from the bound side
as n -> infinity, the sphere's radius growing without bound, which is the same move as the large-spin
door: a compact group contracts. Here SO(4) contracts to the Euclidean group E(3) = ISO(3).

## What I know before computing (disclosed)

Textbook, from memory: the quantum Runge-Lenz vector A = (p x L - L x p)/2 - Z r/r and Pauli's
relations [A_i, A_j] = -2i H eps_ijk L_k, A.L = 0, A^2 = Z^2 + 2H(L^2 + 1); the zero-energy Coulomb
functions u(r) proportional to sqrt(r) J_{2l+1}(sqrt(8 Z r)); the Laguerre-to-Bessel limit
(Mehler-Heine type). Worked by hand before computing: n^(3/2) R_nl(r) -> sqrt(2) Z J_{2l+1}(sqrt(8Zr))
/ sqrt(r); and the bound density per unit energy is A_nl^2 (dn/dE) phi_E(r)^2, with phi_E the regular
solution normalized at the origin (entire in E) and A_nl^2 dn/dE = 4Z (2Z)^(2l) prod_{s=1..l}
(1 - s^2/n^2) / ((2l+1)!)^2, a polynomial in 1/n^2 = -2E/Z^2. So most of what follows is nearly a
corollary of J13, as J4 was of Paper 59; the registration exists to test it, and to place the
symmetry and the price.

## ZE1 (what the release leaves standing: the contraction)

On zero-energy states Pauli's commutator [A_i, A_j] vanishes and A^2 = Z^2, A.L = 0: the zero-energy
solutions carry the E(3) representation with momentum radius Z and helicity zero. l, m and the Gaunt
rules survive; the cap l <= n - 1 is released, and each l >= 0 occurs once at E = 0.

## ZE2 (main: the limit, from both sides)

At fixed r, l and Z:

  R_nl(r)^2 n^3 / Z^2  ->  2 J_{2l+1}(sqrt(8 Z r))^2 / r      as n -> infinity,

with relative corrections in powers of 1/n^2 only (no 1/n term, although e^(-Zr/n) alone would
suggest one). From the continuum side the energy-normalized density (2/(pi k)) F_l^2 / r^2 tends to
the same function as k -> 0+, with corrections in powers of k^2 plus the factor 1/(1 - e^(-2 pi eta)),
whose departure from 1 is smaller than any power of k.

## ZE3 (the price)

The zero-energy radial density 2 J_{2l+1}(sqrt(8Zr))^2 / r is r^(2l) times a power series in r with
rational coefficients (times (8Z)^(2l+1)), an entire function with no pi and no exponential of pi. The
only pi at zero energy is the solid angle, Y_00^2 = 1/(4 pi), so the contact density stays Z/pi. The
exponential of pi of J13 is flat at E = 0 and does not reach the threshold function. Tag, to be
checked against Paper 18 and Paper 34 as part of the test: Bessel functions of the square root of a
rational multiple of r, a 0F1 class; predicted to have no rung on the ladder and no entry in the
projection catalogue under that name.

## Tests and decision rules

1. ZE1: the operators A and L built symbolically in Cartesian coordinates (sympy), applied to the
   zero-energy functions for (l, m) = (0, 0), (1, 0), (1, 1), (2, 0); the commutators and A^2 - Z^2
   evaluated at random points to at least 25 digits; a control at E < 0 where [A_x, A_y] must not
   vanish.
2. ZE2: bound side with mpmath Laguerre polynomials at n up to at least 1000, l = 0, 1, 2, Z = 1, 2,
   r in {0.5, 2, 10}; fit the ratio in 1/n with integer powers, and a fit allowing odd powers whose
   odd coefficients must vanish. Continuum side with mpmath's Coulomb function at k down to 1e-3.
3. ZE3: the Taylor coefficients of J_{2l+1}(2 sqrt x)^2 / x^(2l+1) checked rational symbolically to
   order 10; Paper 18 and Paper 34 searched for Bessel / 0F1 entries.

## Falsifiers

1. ZE1: [A_i, A_j] nonzero on a zero-energy state, or A^2 != Z^2 there.
2. ZE2: a different limit function or constant; a nonzero 1/n coefficient.
3. ZE3: a pi or exponential of pi in the radial threshold function.

## Honest caps, stated before the result

- Textbook mathematics (Pauli 1926; Bander and Itzykson 1966; the Bessel limit of Laguerre
  polynomials). The contribution, if any, is the accounting.
- Non-relativistic, one electron, point nucleus.
- ZE2 is close to a corollary of J13; a pass adds little beyond the symmetry and the price.

## Home, if confirmed

Proposed, PI decides: an extension of Paper 7 `sec:continuum` (the third shape) and of seam row J13,
or a row of its own; a guard in `tests/`.

## Result (2026-10-18, same session, after commit bd2435d1)

Script `debug/zero_energy_door.py`; data `debug/data/zero_energy_door.json`. All three registered
predictions hold.

- **ZE1: CONFIRMED.** The quantum Runge-Lenz vector and L, built symbolically in Cartesian
  coordinates, applied to the zero-energy functions for (l, m) = (0,0), (1,0), (1,1), (2,0) at
  Z = 1, 2: [A_x, A_y] psi, (A^2 - Z^2) psi and (A.L) psi vanish to 1e-33 or better relative to psi, at
  three random points. The zero-energy function was entered as its Bessel power series truncated at
  18 terms; at 12 terms the residuals were 1e-18 to 1e-23 and they fell by the truncation's own
  factor, so they are truncation, not physics. Control at E < 0 (2p, m = 1): the commutator is
  nonzero and equals -2i E L_z psi to 7e-32, and A^2 psi = (Z^2/4) psi to 7e-32, which also validates
  the operator code. So at E = 0 the Runge-Lenz algebra is the Euclidean algebra E(3), with
  A^2 = Z^2 and A.L = 0.
- **ZE2: CONFIRMED.** Bound side, n = 200..1000, l = 0, 1, 2, Z = 1, 2, r = 0.5, 2, 10: the limit
  R_nl^2 n^3/Z^2 -> 2 J_{2l+1}(sqrt(8Zr))^2 / r holds to 5e-16 or better, the 1/n coefficient is
  at most 2e-12 and the 1/n^3 coefficient at most 3e-6 against 1/n^2 coefficients of order 0.08 to 170.
  Continuum side, k = 1/100 .. 1/800, l = 0, 1, r = 2: the energy-normalized density tends to the
  same function to 2e-9, with a k^1 coefficient at most 3e-6 against k^2 coefficients of order 0.08
  to 7. Not registered, measured: the bound 1/n^2 coefficient equals minus Z^2 times the continuum
  k^2 coefficient to three or four digits in all four cases (for example -7.3513 and 7.3501 at
  Z = 1, l = 0), so the two sides share one energy coefficient: one function of E across threshold.
- **ZE3: CONFIRMED.** The Taylor coefficients of J_{2l+1}(2 sqrt x)^2 / x^(2l+1) are rational to the
  order checked (l = 0, 1, 2), and the series matches sympy's Bessel function to 2e-23 or better. So
  the zero-energy radial density has no pi and no exponential of pi; the only pi at E = 0 is the
  solid angle 1/(4 pi). Tag: Papers 18 and 34 have no entry for Bessel functions of position (the
  corpus's Bessel content is the K_0 moments of the period ladder, Paper 59, a different object); so
  the threshold function is an entire function with rational Taylor coefficients, a 0F1 class, and no
  rung is assigned.

Defects in my own test code, caught before any verdict: the first ZE1 run applied the operators to
the Bessel functions symbolically and did not finish (stopped; replaced by the truncated series, whose
truncation error was then shown to be the residual).

What the door says, in one line: at E = 0 SO(4) contracts to E(3) with A^2 = Z^2, l, m and Gaunt
survive with every l present, and the threshold function is a pi-free entire function that both sides
approach in whole powers of E; it is nearly a corollary of J13, and adds the symmetry and the price.

Home (2026-10-18, PI-approved): Paper 7 `sec:continuum`, as the third shape (eq:zero_energy); seam row J13 and atlas door 23 extended; guard `tests/test_paper7_zero_energy_door.py`, written as its own step and fire-tested 9/9.
