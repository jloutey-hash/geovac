# Pre-registered prediction: releasing spin by the large-spin contraction

**Registered:** 2026-10-18, before any computation in this file (PI-directed: "is it possible to
decompactify the spin quantum number?", then "let's do both" for this door and the hydrogen
continuum, `2026-10-18_hydrogen_continuum.md`).
**Extends:** the released-axis rows of `docs/forced_free_seam.md` (sub-domain J) and `docs/atlas.md`;
Paper 24 (the Bargmann-Segal oscillator, which is this door's endpoint).
**Commit status:** committed on its own, before any computation.

## The door

An electron's spin is fixed at 1/2, so no molecular parameter releases it. Spin becomes continuous
only if the group changes. The Euclidean way is the contraction of a spin S as S grows (Holstein and
Primakoff 1940; Inonu and Wigner 1953): near the pole m = S the sphere of a spin S looks flat, and
the spin algebra becomes the oscillator algebra. Write n = S - m for the number of quanta counted
from the pole. The compact side is the spin-S label set, 2S + 1 states; the released side is the
oscillator of Paper 24.

## What I know before computing (disclosed)

Textbook, from memory: the Holstein-Primakoff map; the spin coherent states on the sphere and their
resolution of the identity with measure (2S+1)/(4 pi) dOmega (Radcliffe 1971; Arecchi, Courtens,
Gilmore and Thomas 1972); de Moivre-Laplace (the symmetric binomial tends to a Gaussian); Stirling's
series; the Krawtchouk-to-Hermite limit of the Askey scheme. The first-order coefficient c_1 in LS3
was derived by hand from Stirling's series before any computation. Everything below is therefore
informed. What is tested is the accounting: what the release leaves standing, and the number class
it lets in.

## Conventions

S_z |S,m> = m |S,m>; S_- |S,m> = sqrt((S+m)(S-m+1)) |S,m-1>. Oscillator: a|n> = sqrt(n) |n-1>,
position x = (a + a^dagger)/sqrt(2), ground density pi^(-1/2) e^(-x^2), Hermite functions
psi_n(x) = (2^n n! sqrt(pi))^(-1/2) H_n(x) e^(-x^2/2).

## LS1 (what the release leaves standing)

The label n survives as the oscillator quantum number; the only thing released is the cap n <= 2S.
In the scaled operator A_S^dagger = S_- / sqrt(2S) the matrix element is exactly
<n+1| A_S^dagger |n> = sqrt(n+1) sqrt(1 - n/(2S)), so it tends to sqrt(n+1) at fixed n with relative
error n/(4S) + O(S^-2), and [A_S, A_S^dagger] = 1 - n/S on |n>. Non-uniform: the limit fails for
n of order S.

## LS2 (the skeleton stays algebraic at every S)

At every finite S, every label-space quantity tested here (the matrix elements of S_+, S_-, S_z, and
the components of the S_x eigenvectors in the S_z basis, the Wigner d-matrix at angle pi/2) is
algebraic with a rational square. No pi appears at finite S on the label side.

## LS3 (main: the released spectrum and its density)

X_S = S_x / sqrt(S) has eigenvalues m/sqrt(S), m = -S..S: spacing 1/sqrt(S), extent +-sqrt(S), so
the spectrum becomes dense in the whole real line. In the pole state |S,S> the probability of
eigenvalue m is P(m) = C(2S, S+m) / 4^S. The density p_S(x) = sqrt(S) P(x sqrt(S)) at lattice points
satisfies, at fixed x,

  p_S(x) = pi^(-1/2) e^(-x^2) [ 1 + c_1(x)/S + c_2(x)/S^2 + ... ],   c_1(x) = -1/8 + x^2/2 - x^4/6,

with no half-integer powers of 1/S (the binomial is symmetric), and every c_k(x) a polynomial in
x^2 of degree at most 2k with rational coefficients. The series is asymptotic, not convergent
(Stirling).

## LS4 (eigenvectors)

The components S^(1/4) <S, S-n | X_S = x> tend to (-1)^n psi_n(x) or psi_n(x), the sign fixed by
convention, at fixed n and x, with relative error O(1/S) and no half-integer powers, for n = 0..4.
Confidence: high for n = 0 and 1 (closed forms), moderate for n = 2..4.

## LS5 (the price, and where pi sits)

The transcendental content the release lets in is exactly pi^(1/2) and e^(-x^2), and both enter only
in the limit, as limits of algebraic label data: p_S(0) = sqrt(S) C(2S,S)/4^S is sqrt(S) times a
rational at every S (Wallis's limit gives pi^(-1/2)), and p_S(x)/p_S(0) is rational at every S and
tends to e^(-x^2). On the continuous side, pi is present at finite S only in the sphere's
coherent-state measure (2S+1)/(4 pi) dOmega, and that measure contracts to Paper 24's Bargmann
measure pi^(-1) e^(-|w|^2) d^2w under z = w/sqrt(2S). Tag, to be checked against Paper 18 and
Paper 34 as part of the test: volume content (Vol(S^2) = 4 pi, the Gaussian integral sqrt(pi)),
candidate M1. Confidence in the tag: moderate; the tag is a classification, not a computation.

## Tests and decision rules

1. LS1, LS2: exact arithmetic (sympy / fractions) for S up to at least 200, and symbolic in S where
   the closed form allows.
2. LS3: p_S(x) as exact rationals times sqrt(S), at S = N^2 for N up to at least 400, at x in
   {0, 1/2, 1, 3/2, 2} (lattice points when the denominator divides N); fit c_1 and c_2 by Richardson
   extrapolation; c_1 must match the registered polynomial to at least 12 digits; PSLQ on the
   extrapolated c_2 at each x against {1, x^2, x^4, x^6, x^8} with rational coefficients, and a
   control basis adding log 2, gamma_E and zeta(3), which must enter with zero coefficients.
3. LS4: two routes, the closed-form d-matrix (Krawtchouk) in exact arithmetic and a numerical
   eigendecomposition of S_x at 50 digits, agreeing; then the rate.
4. LS5: symbolic check of the measure contraction; the Wallis and exponential limits checked
   numerically to at least 20 digits.
5. Bounded searches state their bounds; every positive PSLQ result is re-checked at full precision
   with a planted control of the same height.

## Falsifiers

1. LS1: any deviation from the exact matrix element; a label other than the cap released.
2. LS2: a label-space quantity at finite S whose square is irrational.
3. LS3: a half-integer power of 1/S at fixed x; c_1 different from the registered polynomial; a
   coefficient needing a transcendental beyond pi^(1/2) and e^(-x^2) (log 2, gamma_E, zeta(3)).
4. LS4: a limit other than the Hermite functions, or a rate O(1/sqrt(S)).
5. LS5: a transcendental at finite S on the label side; a contracted measure other than Bargmann's.

## Honest caps, stated before the result

- The mathematics is textbook (de Moivre 1733; Laplace; Holstein-Primakoff 1940; Inonu-Wigner 1953;
  Radcliffe 1971; the Askey scheme). The contribution, if any, is placing the door on the seam: what
  it leaves standing, what it lets in, and where the pi sits.
- One contraction only (near a pole). The classical limit S_z/S, whose spectrum fills [-1, 1]
  without decompactifying, is a different limit and is not tested.
- The other two routes to continuous spin are not tested: the non-compact real form SU(1,1)
  (Liu-Noui, cited in Paper 18) and Wigner's massless little group ISO(2) with its continuous-spin
  representations, which is Lorentzian and outside the descoped Paper 45 program.

## Home, if confirmed

Proposed, PI decides: an atlas door and a seam J row (spin released by contraction), with its home
in Paper 24 as the spin-side entrance to the Bargmann space, and a guard in `tests/`.

## Result (2026-10-18, same session, after commit fb4fdd85)

Script `debug/large_spin_door.py`; data `debug/data/large_spin_door.json`. All five registered
predictions hold.

- **LS1: CONFIRMED [SYMBOLIC].** The scaled element is sqrt(n+1) sqrt(1 - n/(2S)) symbolically, with
  relative error n/(4S) + n^2/(32 S^2); [A_S, A_S^dagger] = 1 - n/S on |n> exactly for S = 1..20.
- **LS2: CONFIRMED.** The closed-form d-matrix columns at pi/2 are S_x eigenvectors to 2e-50 at
  S = 3, 7, 12; their squares are rational by construction.
- **LS3: CONFIRMED.** p_S(x) is exactly rational at S = N^2 (sqrt(S) = N). Fit over N = 40..400:
  the limit is pi^(-1/2) e^(-x^2) to 1e-35 or better, and c_1 matches -1/8 + x^2/2 - x^4/6 to 1e-29
  at all seven x in {0, 1/2, 1, 3/2, 2, 5/2, 3}. A fit allowing half-integer powers gives their
  coefficients at most 2e-14, against order-one integer ones (consistent with zero at that fit's
  precision). Not registered, measured: c_2(x) = 1/128 - 11x^2/48 + 19x^4/48 - 3x^6/20 + x^8/72,
  one polynomial fitting seven points with five coefficients to 9e-27, degree 2k in x^2 as
  registered; PSLQ returns c_2(1) = 221/5760 and c_2(2) = -3571/5760 with zero weight on log 2,
  gamma_E and zeta(3), and a planted 221/5760 + 3 log 2 is recovered. At x = 1/2 the search is
  inconclusive at this precision (bound 1e5, tolerance 1e-20; the denominator is larger). The
  Wallis limit: p_S(0) at N = 400 differs from pi^(-1/2) by -4.4e-7, which is
  -pi^(-1/2)/(8S), the registered c_1(0).
- **LS4: CONFIRMED.** For n = 0..4 at x = 1/2, 1, 3/2, S^(1/4) <S, S-n | X_S = x> tends to psi_n(x)
  to 1e-23 or better, rate O(1/S), half-integer coefficients consistent with zero (at most 2e-12),
  once the d-matrix's convention phase (-1)^(S-m) is divided out; the residual sign is +1 for every
  n. Second route: an mpmath eigendecomposition of S_x at S = 30 agrees with the closed form to
  3e-50.
- **LS5: CONFIRMED [SYMBOLIC], tag as a reading.** (2S+1)/(2S) (1 + |w|^2/(2S))^(-2S-2) tends to
  e^(-|w|^2) and C(2S, n)/(2S)^n to 1/n! for n = 0..5: the sphere's coherent-state measure contracts
  to Paper 24's Bargmann measure. On the label side no pi appears at any finite S; pi^(-1/2) and
  e^(-x^2) enter only as limits of rational label data. Paper 18 reads pi = Vol(S^2)/4 as a
  normalization of the M1 volume slot and tiers that reading an Observation
  (`sec:mellin`, the M1 list); the spin door's pi is the same Vol(S^2) normalization, so it takes
  the same tag at the same tier: M1, read, not derived.

Defects in my own test code, caught before any verdict: the first LS4 run fixed the d-matrix's
sign from the last point only, so the ratio alternated with the parity of S - m at x = 1/2 and 3/2
(fixed by dividing out (-1)^(S-m)); the first PSLQ call on c_2 used a coefficient bound of 1e8 at
25 digits and returned a spurious relation, then a tolerance tighter than c_2's precision and
returned nothing (fixed: bound 1e5, tolerance 1e-20, with a planted control).

What the door says, in one line: releasing spin by contraction leaves the label n standing, removes
only the cap n <= 2S, and lets in exactly pi^(1/2) and e^(-x^2), both as limits of algebraic label
data; the pi is the sphere's volume normalization. The mathematics is textbook; this is the
accounting.

Home (2026-10-18, PI-approved): Paper 24 `sec:large_spin`; seam row J12 (F, with the measured corrections and the M1 reading outside it); atlas door 22; guard `tests/test_paper24_large_spin_door.py`, written as its own step and fire-tested 11/11.

## Novelty classification (appended 2026-10-19)

Classes per `docs/predictions/README.md`, set after the fact by a literature audit (`debug/prediction_novelty_audit_memo.md`; this record predates the novelty field). FINDING = held, not found in the literature; CHECK = held, already known; MISS = failed.

- LS1-LS5: CHECK (Holstein-Primakoff; Krawtchouk-Hermite; Radcliffe)
