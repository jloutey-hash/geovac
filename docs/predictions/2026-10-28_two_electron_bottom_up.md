# Pre-registered prediction: the two-electron bond-angle structure, built bottom-up (B2)

**Registered:** 2026-10-28, before any two-electron bottom-up derivation. PI-directed ("Let's take that natural next
step": build the two-electron equation up from the label rules).
**Commit status:** committed on its own, before the derivation is run.
**Rests on:** AO3 (STOP: no x-operator below about 1970 unknowns; 31 constants; components of orders 6-13), BU (GO:
the one-electron operator derived from the label rules through a first-order system in L), NA (the one-electron
singular structure).

## Why the recurrence in L is enough

The Legendre coefficients' recurrence in L carries the angle equation's singular structure: a characteristic root
lam maps to the point x = (lam + 1/lam)/2, and a term (x0 - x)^alpha there to formal solutions lam^L L^beta with
beta = -alpha - 1/2. **Calibration** (`debug/b2_asym.py` on BU's order-9 recurrence for NA's N1, both primes,
measured before registering): roots 9/13 (twice; beta -2, -3), 13/9 (beta -2), both mapping to x* = 125/117 (NA's
alpha 3/2, 5/2); roots -300/299 and -299/300 (beta 3/2), mapping to -179401/179400 (NA's alpha -2 point); root 0
twice. x = +-1 and the apparent block do not appear. So a branch point (half-integer alpha) shows as integer beta,
and an integer alpha as half-integer beta.

## The task

Case W1 of AO3 (D1 = 17/10, D2 = 5/2, zX = 7/5, zH = 11/10; x* = 457/425, rho = D</D> = 17/25). Derive, from the label
rules and exact radial integration with no fit to coefficient data (as BU), a first-order system in L reproducing the
two-electron coefficients, then the recurrence in L of the generic sequence (AO3's random weights), then its large-L
structure by `debug/b2_asym.py`. The x-operator is optional (likely too large).

## Predictions

- **B2-0 (checks).** The system's output equals the exact coefficients of `debug/ao_exact.py` (generic sequence,
  random.Random(2026) weights as AO3) on L <= 300 on both primes, and the derived recurrence annihilates them.
- **B2-1 (the gate).** A recurrence in L is derived within the agent's budget and passes B2-0. Probability 0.6.
- **B2-2 (only branch point).** Among the rational characteristic roots, integer beta occurs only at lam = 17/25 and
  25/17 (both mapping to x*); every other rational root has only half-integer beta. Non-rational roots are reported.
  Probability 0.65.
- **B2-3 (the contact form).** At lam = 17/25 every beta is an integer <= -4 and -4 occurs: the slowest bond-angle
  branch is (x* - x)^(7/2), the cusp-cusp term. Probability 0.7.

**Reported, not scored:** system dimension, recurrence order and degree, all roots and exponents, and the x-operator
if obtained.

**Gate.** STOP if B2-1 fails. GO if B2-1, B2-2, B2-3 hold. BORDERLINE otherwise.

## Caps

One geometry, 1s functions; the linear independence of the 31 constants over Q is assumed for the generic sequence.

## Disposition (appended 2026-10-28, after the derivation; nothing above edited)

**Outcome: BORDERLINE.** B2-0, B2-1 and B2-2 held. B2-3 missed as worded.

Derivation by a sub-agent (memo `debug/b2_memo.md`; scripts `debug/b2_{engine,derive,data,check,charpoly}.py`).
- Route: the integral split so that each inner integral starts at 0 or ends at infinity, so no constants enter. Bessel
  modules with integration-by-parts identities reduced mod p, closed under the L-shift: a first-order system of
  dimension 180.
- Recurrence in L by the cyclic-vector method, 4200 points, rational reconstruction: order 41, degree 926, the same on
  both primes.

**PM verification** (`debug/b2_pmcheck.py`). Independence: `b2_derive.py` reads only `b2_engine.py` and two helper
functions from `bu_derive.py`, with no coefficient data. Annihilation: the generic sequence was rebuilt here straight from
`debug/ao_exact.py` with AO3's weights, to L = 400, past the 300 the record names. Result: 0 of 360 rows nonzero on both
primes. With one coefficient changed, 360 of 360 rows are nonzero. `debug/b2_asym.py` was re-run on both primes and gives
identical output.

**Large-L structure** (both primes):

| root lam | x | betas |
|:--|:--|:--|
| 17/25 (x5) | 457/425 = x* | -4, -5, -6, 3, 1 |
| 25/17 (x4) | x* | -4, -5, 3, 1 |
| -51/121, -121/51 | -8621/6171 | -1/2 |
| -1875/2057, -2057/1875 | -3873437/3856875 | -1/2 |
| 0 (x10) | none | none |

There are also two irreducible quadratics, each twice, with real negative roots mapping to x < -1.

- **B2-2 held.** Integer beta occurs only at 17/25 and 25/17. The two other rational points lie past x = -1 and carry
  only beta = -1/2, which corresponds to an integer alpha.
- **B2-3 missed.** At 17/25 the contact exponent beta = -4, i.e. (x* - x)^(7/2), is present and is the slowest negative
  one. But beta = 3 and beta = 1 also occur, so not every beta is <= -4.

**Open, not scored.** It is not known whether beta = 3 and 1 belong to the integral's coefficients, or only to the
recurrence's larger solution space. The recurrence's space has dimension 41. In BU the derived recurrence carried one
solution more than the sequence needed, and the 25/17 solutions grow, so they too must be absent from the decaying
coefficients. Settling this needs data: either a large-L numeric test of F_L or the sequence's minimal recurrence.
