# Open-confinement aperture flag — two-sided (electron + gravity)

**Date:** 2026-05-30
**Verdict:** POSITIVE-STRUCTURAL (universal arrow), with system-specific rate and bound geometry.
**Flag:** Opening a confinement aperture is a *symmetry contraction accompanied by an entropy climb*, and the
**arrow is universal across the electron and gravity sectors**; the *rate* and *bound geometry* are system-specific.

## Survey (Task 0, explorer — all citations web-verified)
`debug/explorer_open_confinement_inheritance_memo.md`. INHERIT: Bander-Itzykson 1966 (SO(4)->E(3)->SO(3,1)
contraction, parameter = our p0 = sqrt(-2E)); Seaton 1983 (quantum defect bridges E<0/E>0); Bialynicki-Birula-
Mycielski 1975 (entropic floor S_x+S_p >= d(1+ln pi)); Rennie/Gayral non-unital triples (the R^3/H^3 limit object).
BUILD (out of scope, named): full spatial propinquity S^3 -> R^3/H^3 — no published result runs this direction;
the temporal S^1 -> R machinery (Papers 47/48) does NOT transport directly (curved+double-limit+signature-crossing).
Gift: no published work frames the gravity beta->inf cigar limit as an Inonu-Wigner contraction — that framing is ours.

## A — contraction (backbone, inherited)
- Electron: p0 = sqrt(-2E) IS Bander-Itzykson's contraction parameter. p0 -> 0 (shells accumulate at threshold) =
  SO(4) -> E(3); E>0/H^3/SO(3,1) is the scattering sector (RH-B's dead end, deliberately not crossed).
- Gravity: Paper 47 beta -> inf, S^1 -> R, U(1) -> R is the temporal-axis contraction (proven; here re-read as a
  contraction, an open framing).

## B — entropy climb (stretch, computed)
ELECTRON (`debug/open_confinement_electron_entropy.py`, verified vs exact 1s S_x = 3+ln(pi) = 4.14473 to all digits;
norms at 1e-16):

  n   p0=1/n     E_n       S_x      S_p    S_x+S_p   above BBM floor(6.434)
  1   1.0000  -0.50000   4.145    2.422    6.567    0.132
  2   0.5000  -0.12500   8.111   -0.758    7.353    0.919
  4   0.2500  -0.03125  12.075   -3.146    8.929    2.495
  7   0.1429  -0.01020  15.299   -4.938   10.361    3.926
  10  0.1000  -0.00500  17.365   -6.052   11.313    4.879

  - S_x+S_p climbs monotonically as aperture p0 -> 0; margin above inherited BBM floor grows 0.13 -> 4.88.
  - Split: S_x rises (orbital delocalizes, Dr ~ n^2), S_p falls (momentum collapses to threshold, Dp ~ 1/n).
    The climb lives in the scale-invariant SUM = phase-space-volume growth (~ n per dim).
  - Rate: climb fit n>=5 gives slope 2.63 in ln n, rising toward the 3 ln n (3D phase-space) prediction from below.
  - LAYER 2: S_x+S_p is a log-of-density, transcendental, numerically integrated — NOT bit-exact, NOT pi-free.
    Exactly as the costume-table predicts: entropy is the open-regime currency.

GRAVITY (`debug/open_confinement_gravity.py`, verified vs exact first law dM=TdS to 5e-9, C=-8 pi M^2 to 1e-8):

  - Aperture 1/beta -> 0 (beta=8 pi M -> inf, T -> 0, M -> inf). S_BH = 4 pi M^2 = beta^2/(16 pi) ~ (1/aperture)^2:
    QUADRATIC entropy divergence as the aperture opens. Same arrow as the electron.
  - 2pi in T_H = 1/(8 pi M) is the M1 Hopf-base-measure / Matsubara-circumference signature
    (beta * kappa_g = 2pi four-witness law) — temperature = confinement of imaginary time.

## Stretch #3 — density-of-states at the open end (the direct twin)
- ELECTRON: from our shells g_n = n^2, N(|E|) ~ |E|^{-3/2} (exponent = 1.5 at every decade, exact) => dN/dE ~ |E|^{-5/2}.
- GRAVITY: rho(M) = exp(S) = exp(4 pi M^2) => exponential divergence.
- Both DOS diverge as the aperture opens. Confinement removes states; opening floods them.

## Universality verdict (honest scope)
- ARROW: UNIVERSAL. Both entropy and DOS diverge monotonically as the aperture opens. This is the confidence-giving claim.
- RATE: SYSTEM-SPECIFIC. Electron power-law (|E|^{-5/2}, S~log) from the Coulomb 1/n^2 ladder; gravity exponential
  (exp M^2, S~M^2) from the area law. Not universal — and shouldn't be.
- BOUND: the real contrast. Electron entropy sits ABOVE the BBM floor (open above, ground state nearly saturates it);
  gravity entropy SATURATES the holographic ceiling A/4 (bound achieved). Floor vs saturated ceiling.
  Caveat: electron S_x+S_p is an uncertainty-entropy of a pure state; gravity S_BH is a thermodynamic/entanglement
  entropy. Different entropy SPECIES (cf. TD Track 5: state-side entropy vs spectral-side pi). The floor/ceiling
  contrast may partly reflect the species difference, not only the sector — flagged, not resolved.

## What this plants, and what it does not
- PLANTED: a verified two-sided instance of "open the aperture -> contraction + entropy climb," anchored on both ends
  to inherited exact results (BBM floor; Schwarzschild first law). Confinement = a physical regime with a measurable
  open-end signature, the same arrow for electrons and black holes.
- NOT crossed: the actual scattering sector (E>0, H^3) and the full spatial propinquity (S^3 -> R^3/H^3) remain the
  multi-month BUILD frontier. Explorer staging: do the E=0 / E(3) / flat-R^3 fixed point first.
- The "entropy = thermodynamic face of confinement" reading is PROPOSED-PRINCIPLE strength, now with two verified
  data points, not a theorem.

## Files
debug/open_confinement_electron_entropy.py, debug/check_1s.py, debug/open_confinement_gravity.py,
debug/data/open_confinement_electron_entropy.json, debug/data/open_confinement_gravity.json,
debug/explorer_open_confinement_inheritance_memo.md
