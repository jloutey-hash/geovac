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

## Result (2026-10-19)

Script `debug/relativistic_nucleus_door.py`, log `debug/data/relativistic_nucleus_door.log`. The run
was stopped by the operating system for low memory after RD3 finished, before RD4. RD2 and RD3 are
read from that log. RD4 and RD5 were then run as separate stages (logs `_rd4.log`, `_rd5.log`); RD5
now integrates the exterior one grid cell at a time, so its memory stays bounded.

**Code defects found and fixed before the scored run (not physics):**
1. The first matching function was the ratio of interior and exterior F/G. It has a narrow pole in W,
   and the secant step jumped past c^2, where the square root turns complex. Replaced by the
   Wronskian form, which is smooth in W, R and Z, and a secant start at the point-nucleus level.
2. At Z alpha = 1 the shift is not small, so a start at the point level (W = 0) fails the same way.
   Added a bracketing search that walks down from W = 0.6 c^2 to the first sign change.
3. RD4: mpmath's anderson and illinois root finders both stalled in the bracket. Replaced by a
   hand-written Illinois step with a step cap and a progress print.

**RD2: CONFIRMED.** Local slopes d ln(Delta W)/d ln R, smallest R first:

| Z alpha | 2 gamma | slopes (R = 1e-8..1e-5) | miss at smallest R |
|:--|:--|:--|:--|
| 0.3 | 1.907878 | 1.907876, 1.907857, 1.907665 | 2e-6 |
| 0.6 | 1.600000 | 1.599995, 1.599947, 1.599471 | 5e-6 |
| 0.9 | 0.871780 | 0.871717, 0.871360, 0.869146 | 6e-5 |

All inside 1e-3. The clause "at Z = 1 it reduces to the non-relativistic law" was not computed in
this run; it is unscored.

**RD1:** holds by construction (the potential is round, so kappa, m_j and parity label the levels).
The p 1/2 shift at Z alpha = 0.6 was computed but only written to the JSON, which the stop prevented.

**RD3: MISSED on the registered threshold.** W/c^2 at Z alpha = 1:

| R | 1e-3 | 1e-4 | 1e-5 | 1e-6 | 1e-7 | 1e-8 | 1e-9 |
|:--|:--|:--|:--|:--|:--|:--|:--|
| W/c^2 | 0.34484 | 0.20752 | 0.14377 | 0.10900 | 0.08751 | 0.07300 | 0.06259 |

Local slopes per decade: 0.221, 0.159, 0.120, 0.095, 0.079, 0.067. The slope keeps falling, so no
power law holds over the range and falsifier 2 does not fire. But the registered number was "below
0.05" by R = 1e-9, and the last slope is 0.067. Scored as missed. The unregistered guess
W ~ 1/ln(1/R) also fails: W ln(1/R)/c^2 runs 2.38, 1.91, 1.66, 1.51, 1.41, 1.34, 1.30, still falling.
W goes to zero slower than any power but faster than 1/ln(1/R) over this range.

**RD4: CONFIRMED.** With R = 1.2 (2.5 Z)^(1/3) fm, the matching function at W = -c^2 changes sign
once between Z = 160 and 184, between 168 and 170. The root is Z_cr = 169.7065, residual 8e-26, at
R = 1.704e-4 bohr (9.02 fm). Inside [165, 180]. Consistent with the literature values near 170 for
similar radius laws.

**RD5: CONFIRMED.** At Z alpha = 1.2, sigma = 0.663325, e^(pi/sigma) = 113.9921. The radii at which
successive s 1/2 levels reach W = -c^2:

| level | 1 | 2 | 3 | 4 | 5 |
|:--|:--|:--|:--|:--|:--|
| R_k (bohr) | 1.04464e-4 | 9.13329e-7 | 8.01197e-9 | 7.02853e-11 | 6.16581e-13 |

Ratios R_k/R_(k+1): 114.377, 113.995, 113.992, 113.992. Within 0.34 percent at the first ratio and
3e-5 by the second, against the registered 2 percent by the third.

Cross-check between the stages: at Z = 164.4 the deepest level dives at R = 1.04e-4 bohr, below the
radius law's 1.69e-4 bohr at that charge, where RD4's scan has not yet dived (D = +0.16). Consistent.

**Score: RD1 by construction, RD2 confirmed, RD3 missed on its threshold, RD4 and RD5 confirmed.**

**Price, read off the result.** Below Z alpha = 1 the radius enters as R^(2 gamma). At Z alpha = 1
it enters more slowly than any power. Past it, the radius fixes the levels only through ln R modulo
pi/sigma: the continuous scale invariance of the point problem is broken to a discrete one, with
period e^(pi/sigma), by one consumed number. The radius is compulsory, not corrective. Verdict as
predicted: FREE. The transcendental e^(pi/sigma) is a function of the consumed alpha through sigma,
with pi from the log-periodic phase; it is not on the skeleton side.
