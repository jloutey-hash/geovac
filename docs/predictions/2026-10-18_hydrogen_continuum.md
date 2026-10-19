# Pre-registered prediction: releasing n into the hydrogen continuum

**Registered:** 2026-10-18, before any computation in this file (PI-directed, with the large-spin
door `2026-10-18_large_spin.md`: "let's do both").
**Extends:** the released-axis rows of `docs/forced_free_seam.md` (sub-domain J) and `docs/atlas.md`,
which index no positive-energy door; Paper 7 (Fock's projection; Bander and Itzykson's hyperboloid
for scattering states is cited there as prior art and never computed in the corpus).
**Commit status:** committed on its own, before any computation.

## The door

Bound hydrogen lives on S^3 with the compact symmetry SO(4); n labels its finite irreducible
pieces. At positive energy Fock's projection lands on the unit hyperboloid H^3 instead, the symmetry
is SO(3,1), and n becomes continuous (Bander and Itzykson 1966). The release is the continuation
n -> i eta, eta = Z/k, E = k^2/2 (atomic units). The compact SO(3) subgroup, and with it l and m,
is untouched. This is the PI's picture for the prolate case, applied to n: release one label and ask
what stays forced.

## What I know before computing (disclosed)

Textbook, from memory: the hydrogen oscillator strengths, f(1s -> np) = 2^8 n^5 (n-1)^(2n-4) /
(3 (n+1)^(2n+4)); Stobbe's photoionization formula (1930), with its factors e^(-4 eta arccot eta)
and 1/(1 - e^(-2 pi eta)); the Sommerfeld factor 2 pi eta/(1 - e^(-2 pi eta)) for the continuum
density at the origin; continuity of the oscillator-strength density across threshold (Seaton's
quantum-defect theory). The continuation of ((n-1)/(n+1))^(2n) to e^(-4 eta arccot eta) and the
threshold value 2^8 e^(-4)/(3 Z^2) were worked by hand before any computation. Everything below is
informed. What is tested is the accounting.

## Conventions

Bound states |n l m>, E_n = -Z^2/(2 n^2), density of states dn/dE = n^3/Z^2. Continuum states
energy-normalized, <E l m|E' l' m'> = delta(E - E') delta_ll' delta_mm'. Oscillator strengths in the
length form, summed over final m. The continuation of powers uses the principal logarithm:
((n-a)/(n+a))^(2n) = exp(2n log((n-a)/(n+a))).

## HC1 (what the release leaves standing)

(a) E_n continues to +k^2/2. (b) The S^3 Casimir n^2 - 1 continues to -(1 + eta^2), minus the
spectrum 1 + rho^2 (rho = eta) of the Laplacian on the unit hyperboloid, the sign being the curvature
flip S^3 -> H^3: the S^3 zonal function sin(n chi)/sin(chi), continued by n -> i rho and chi -> i chi,
is the H^3 zonal function sin(rho chi)/sinh(chi) up to a constant. (c) l and m survive with their
Gaunt selection rules; the only thing released is the cap l <= n - 1: at every k > 0 every l >= 0
occurs.

## HC2 (the contact density across threshold)

The s-wave contact density per unit energy is exactly Z/pi on the bound side, for every n:
|psi_n00(0)|^2 dn/dE = (Z^3/(pi n^3)) (n^3/Z^2). On the continuum side it is
(Z/pi) / (1 - e^(-2 pi eta)). It is continuous at threshold, and the difference beyond threshold
vanishes faster than any power of E.

## HC3 (main: one universal factor)

For each tested transition from a bound state i, let Q_i(n) = f(i -> n) dn/dE, the bound-side
oscillator-strength density, a closed form in n. Prediction: the continuum density is

  df/dE (i -> E) = Q_i(i eta) / (1 - e^(-2 pi eta)),

with the same factor for every initial state and every l. For 1s this gives
df/dE = 2^8 eta^8 e^(-4 eta arccot eta) / (3 Z^2 (1 + eta^2)^4 (1 - e^(-2 pi eta))), Stobbe's
formula, and at threshold both sides tend to 2^8 e^(-4)/(3 Z^2). Tested transitions: 1s -> p,
2s -> p, 2p -> s, 2p -> d, and the contact density of HC2. Confidence: high for 1s and the contact
density; moderate-high that the factor is l-independent (my argument: it is the ratio of the
continuum normalization, through |Gamma(l + 1 + i eta)|^2 = (pi eta/sinh(pi eta)) prod_{s=1..l}
(s^2 + eta^2), to the continued bound normalization, through 1/prod_{s=-l..l}(n + s), and the
l-dependent products cancel).

## HC4 (the price)

On the bound side every tested Q_i(n) is rational at each integer n. The release lets in two kinds
of number. (i) Exponentials of eta times arccotangents, the continuation of the rational factors
((n-a)/(n+a))^(2n); their threshold limits, such as e^(-4), are already the n -> infinity limits of
the bound side's rational sequences, so they are shared across threshold. (ii) The single factor
1/(1 - e^(-2 pi eta)), an exponential of pi times eta, which the bound side never shows and which is
flat at threshold. The Coulomb phase arg Gamma(l + 1 + i eta) lives in the amplitudes and is
predicted to cancel from every tested squared quantity. Tag, to be checked against Paper 18 and
Paper 34 as part of the test: Paper 18's ladder is predicted to have no rung for e^(2 pi eta), so the
door adds a class, not a rung; the pi in Z/pi is the solid-angle measure, Y_00^2 = 1/(4 pi).

Disclosed caution, not a prediction: 1/(1 - e^(-x)) at x = 2 pi eta has the Bose-Einstein form,
and the corpus meets that form at temporal compactification (Paper 35, beta = 2 pi). A shared
functional form is the "X matches Y" case. No thermal reading is registered, and any reading of this
kind gets the curve-fit audit before it is written anywhere.

## Tests and decision rules

1. HC1: symbolic (sympy) for (a), (b) and the zonal functions; (c) from the regular Coulomb
   functions existing at every l.
2. HC2, HC3: two independent routes. Route A: the closed form Q_i(i eta)/(1 - e^(-2 pi eta)).
   Route B: direct radial quadrature of the dipole matrix element against mpmath's regular Coulomb
   function F_l(eta, kr), energy-normalized, at k in {0.1, 0.5, 1, 2, 5} for Z = 1 and Z = 2, to at
   least 25 digits. Bound sides in exact arithmetic for n up to at least 60, including the threshold
   limits.
3. HC3 universality: the ratio route-B / Q_i(i eta) for each transition, which must equal
   1/(1 - e^(-2 pi eta)) for all of them at the stated precision.
4. HC4: the tested squared quantities checked for Gamma or log content by PSLQ at the stated
   precision against a basis with and without arg Gamma(l + 1 + i eta), log k and gamma_E, with a
   positive control.
5. Bounded searches state their bounds.

## Falsifiers

1. HC1: the continued zonal function failing the H^3 eigen-equation; an l missing at some k.
2. HC2: the bound side not exactly Z/pi; the continuum not (Z/pi)/(1 - e^(-2 pi eta)).
3. HC3: a tested transition whose continuum density differs from Q_i(i eta)/(1 - e^(-2 pi eta)) at
   the stated precision, or a factor that depends on l or on the initial state.
4. HC4: a tested squared quantity that needs a transcendental outside the stated set.

## Honest caps, stated before the result

- Non-relativistic hydrogen-like ions only: point nucleus, no fine structure, one electron. The
  door at two centres or for many electrons is not touched.
- Squared quantities only. The amplitudes and their phases are not the object.
- The mathematics is textbook (Sommerfeld; Stobbe 1930; Bethe and Salpeter; Seaton 1958; Bander and
  Itzykson 1966). The contribution, if any, is the accounting on the seam: what the release leaves
  standing (l, m, Gaunt), and that the whole continuum price beyond the continued skeleton is one
  factor.

## Home, if confirmed

Proposed, PI decides: an atlas door and a seam J row (n released at positive energy), with its home
in Paper 7 beside Fock's projection, and a guard in `tests/`.

## Result (2026-10-18, same session, after commit fb4fdd85)

Script `debug/hydrogen_continuum_door.py`; data `debug/data/hydrogen_continuum_door.json`. All four
registered predictions hold.

- **HC1: CONFIRMED [SYMBOLIC].** E_n continues to k^2/2; n^2 - 1 continues to -(1 + rho^2); the S^3
  zonal function satisfies its eigen-equation with n^2 - 1, the H^3 zonal function sin(rho chi)/sinh
  chi satisfies its own with 1 + rho^2, and the continued S^3 function is exactly i times the H^3
  one. l and m survive: the continuum densities below are computed at l = 0, 1, 2 at every k, from
  the same initial states.
- **HC2: CONFIRMED.** The bound side is Z/pi symbolically. The continuum contact density, from
  mpmath's own regular Coulomb function at the origin, equals (Z/pi)/(1 - e^(-2 pi eta)) to 4e-29 at
  k in {0.1, 0.5, 1, 2, 5}, Z = 1, 2.
- **HC3: CONFIRMED, l-independent.** Route A, Q_i(n) as one closed form in n, matches the exact
  rational oscillator strengths at integer n = 2..60 to 2e-39 for all four transitions, and the
  registered 1s closed form is exact there. Route B, direct quadrature against mpmath's Coulomb
  functions, energy-normalized. Over 40 rows (1s -> p, 2s -> p, 2p -> s, 2p -> d; five k; Z = 1, 2)
  the ratio route B / Q_i(i eta) equals 1/(1 - e^(-2 pi eta)) to 4e-30; Q_i(i eta) is real to
  6e-30. Thresholds: at n = 1000 and eta = 1000 the two sides agree to 6e-6 for every transition,
  bracketing the 1s value 2^8 e^(-4)/3 = 1.5629345 from above (bound, 1.5629387) and below
  (continuum, 1.5629304).
- **HC4: CONFIRMED.** At eta = 1 and 2, each continuum density with the factor 1/(1 - e^(-2 pi eta))
  and the exponential e^(-4 eta arccot(eta/n_i)) stripped off is a rational (for 1s at eta = 1,
  16/3, so df/dE = (8/3)/sinh(pi)); PSLQ gives zero weight to the Coulomb phase
  arg Gamma(l + 1 + i eta), log k and gamma_E (bound 1e15, tolerance 1e-35), and a planted
  7/13 + 2 gamma_E is recovered. Tag: neither Paper 18 nor Paper 34 has any entry for an exponential
  of pi (searched: e^pi, e^(-2 pi), sinh pi, coth, Gelfond), so the door adds a class, not a rung.

Not registered, structural: in every tested transition the 2F1 of the radial integral has c - b a
non-positive integer, so Euler's transformation makes it (1 - z)^(c-a-b) times a finite polynomial,
with 1 - z = (n - n_i)/(n + n_i). That is where the factors ((n - n_i)/(n + n_i))^(2n), and after
continuation the exponentials e^(-4 eta arccot(eta/n_i)), come from: one closed form in n covers
the bound and free sides, and the only thing the continuum adds beyond it is the factor.

Defects in my own test code, caught before any verdict: mpmath's 2F1 failed to converge at integer
n (catastrophic cancellation in a terminating sum) and at the 2s initial state with eta = 2; both
were replaced by the Euler-transformed closed form above, validated at integer n against the
independent exact-rational route.

What the door says, in one line: releasing n leaves l, m and the Gaunt rules standing, removes only
the cap l <= n - 1, and the whole continuum is the continued bound skeleton times one factor
1/(1 - e^(-2 pi eta)), an exponential of pi that the bound side never shows and that is flat at
threshold. The shared Bose-Einstein form is recorded, not read. The mathematics is textbook; this is
the accounting.

Home (2026-10-18, PI-approved): Paper 7 `sec:continuum`; seam row J13 (F); atlas door 23; guard `tests/test_paper7_hydrogen_continuum_door.py`, written as its own step and fire-tested 12/12. The guard adds a check the computation did not make: the factor holds at the level of the radial functions, for every r and l <= 3, as the ratio of the two normalizations.
