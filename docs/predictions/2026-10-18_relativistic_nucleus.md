# Pre-registered prediction: the finite nucleus near and past Z alpha = 1 (Dirac)

**Registered:** 2026-10-18, before any computation in this file (PI-directed: "let's work those owed
items"; the remaining owed item of atlas door 24, "the relativistic case, where a finite nucleus
becomes compulsory near Z alpha = 1").
**Extends:** `2026-10-18_finite_nucleus.md` and `2026-10-18_nuclear_profiles.md` (the non-relativistic
door, seam row J14); Paper 23 `sec:finite_nucleus`; Paper 36 (bound-state QED on the Dirac S^3).
**Commit status:** committed on its own, before any computation.

## The door

For a point nucleus the Dirac 1s level is W = c^2 gamma (rest energy included), gamma = sqrt(1 - (Z alpha)^2),
alpha = 1/c, atomic units. At Z alpha = 1, gamma = 0 and W = 0; past it gamma is imaginary, gamma = i sigma,
sigma = sqrt((Z alpha)^2 - 1), and the point-Coulomb Dirac operator has no unique self-adjoint
realization: a boundary condition at the origin, one real number, must be supplied. A finite nucleus
supplies it. So the non-relativistic picture of door 24, where the radius is a small correction, turns
into one where the radius is compulsory. The control parameter Z alpha itself contains a consumed
constant, alpha (atlas door 9).

## What I know before computing (disclosed)

Textbook, from memory: the Sommerfeld-Dirac levels; that finite-size shifts scale as R^(2 gamma)
(Shabaev's formula, from memory, not consulted); Case's analysis of singular potentials (1950);
the critical charge near Z = 170 for realistic nuclear radii (Pomeranchuk and Smorodinsky 1945;
Pieper and Greiner 1969; Zeldovich and Popov 1971); log-periodic level structure for supercritical
attraction, the "fall to the centre". Nothing below was computed beforehand.

## Model and method

Uniformly charged sphere of radius R, Dirac radial equations for G = r g and F = r f, kappa = -1 (s 1/2)
and kappa = +1 (p 1/2); c = 137.035999084. Inside, an exact power series in r (the potential is a
polynomial). Outside, the exterior equations integrated inward from large r in the variable ln r,
starting on the decaying solution, and matched at r = R. Check: as R -> 0 at Z alpha < 1 the levels must
reproduce the Sommerfeld formula.

## Predictions

- **RD1 (what survives):** kappa, m_j and parity survive for a round nucleus, as l, m did
  non-relativistically.
- **RD2 (the exponent):** for Z alpha < 1 the 1s finite-size shift Delta W = W(R) - W(0) scales as
  R^(2 gamma): d ln(Delta W)/d ln R -> 2 gamma as R -> 0, checked at Z alpha = 0.3, 0.6, 0.9 to 1e-3
  absolute. At Z = 1 it reduces to the non-relativistic (2/5)(ZR)^2 Z^2 law times 1 + O(alpha^2).
- **RD3 (at Z alpha = 1):** W(R) -> 0 as R -> 0, slower than any power of R: over R = 1e-3 .. 1e-9 the
  local slope d ln W / d ln R falls below 0.05 and keeps falling. Moderate confidence on the form;
  a 1/ln(1/R) law is my guess and is not registered.
- **RD4 (the critical charge):** with R = 1.2 (2.5 Z)^(1/3) fm (bohr = 52917.72109 fm), the 1s level
  reaches W = -c^2 at Z_cr between 165 and 180.
- **RD5 (discrete scale invariance past Z alpha = 1):** at fixed Z alpha = 1.2 (sigma = 0.6633), the radii
  R_k at which successive s 1/2 levels reach W = -c^2 satisfy R_k / R_(k+1) -> e^(pi/sigma) = 114.0 as
  k grows, to 2 percent by the third level.

## The price, predicted

For Z alpha < 1 the consumed radius enters as R^(2 gamma), an algebraic power of a consumed number with
exponent an algebraic function of the consumed alpha. Past Z alpha = 1 the radius enters at leading
order, and only through ln R taken modulo pi/sigma: the continuous scale invariance of the point problem
is broken to a discrete one by one consumed number. Predicted verdict: FREE, the radius compulsory
rather than corrective.

## Falsifiers

1. RD2: a limiting slope differing from 2 gamma by more than 1e-3.
2. RD3: a power law in R with exponent above 0.05 holding over the range.
3. RD4: Z_cr outside [165, 180].
4. RD5: a ratio not tending to e^(pi/sigma).

## Honest caps

- One electron, no vacuum polarization or self-energy, a uniform-sphere nucleus, nuclear radius law
  fixed as above. Real supercritical physics (pair production, the charged vacuum) is not touched.
- The physics is textbook (references above). The contribution, if any, is the accounting on the door:
  what the release keeps, and how the consumed number's role changes at Z alpha = 1.
