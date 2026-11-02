# Pre-registered prediction: how fast the ring momentum dies on a squashed hydrogen ring (SQ)

**Registered 2026-11-02** (PI: "Yep, let's do that.", in reply to the proposed ellipse test), before any scored
case is computed.
**Commit status:** committed on its own, before any computation it scores.
**Rests on:** record HR (`docs/predictions/2026-11-01_hydrogen_ring.md`: ring momentum k is exact on the round
ring); record TR (`docs/predictions/2026-11-02_hydrogen_torus.md`: the tube winding leaks as (r/R)^|delta m|).
**Guardrail check (Papers 8-9, FCI-M, Track DF):** not in scope. One atom-centred 1s function per hydrogen; no
single-centre or shared-exponent molecular basis; no graph concatenation.

## Why

A review of the corpus's geometries (2026-11-02) read three label-breaking laws side by side. Bending a linear
molecule opens m at order s^(2 Lambda) (atlas door 15, Renner); curving a tube into a torus leaks the tube winding
as (r/R)^|delta m| (TR); a second nucleus turns the forbidden coupling on through an exponential tail (door 12).
The reading: a change of shape breaks a label as a power, and the power counts how much the label changes. Two
shape instances is thin, and they measure different quantities (an energy splitting; a matrix-element share). SQ is
a third shape instance with a sharp prediction, including a case where the power is not the obvious one.

## The order rule (derived before any computation)

Put atom j at parametric angle phi_j = 2 pi j / N on an ellipse: position R((1 + eps) cos phi_j, (1 - eps) sin phi_j, 0),
that is z_j = R(e^{i phi_j} + eps e^{-i phi_j}). Rotating the molecule by 2 pi / N maps the configuration with complex
deformation eps to the one with eps e^{2i(2 pi / N)}. The Lowdin site Hamiltonian is a smooth function of (eps,
conj eps), and Lowdin orthogonalization commutes with the rotation. So a term eps^{p+} conj(eps)^{p-} changes the ring
momentum by 2(p+ - p-) mod N: the squash is a "spurion" of charge 2. The part of h and of (pq|rs) that changes k by
c (mod N, up to sign) is therefore of order |eps|^{p(c)}, with

  p(c) = min |q| over integers q with 2q = +-c (mod N).

- Even N: odd c is exactly zero (the half-turn j -> j + N/2 survives the squash); even c has p = c/2.
- Odd N: even c has p = c/2; odd c has p = (N - c)/2. This wrap-around order is the sharp part: at N = 11 a change
  of k by 1 needs five powers of eps, not one.

The rule gives the order as a lower bound; that the leading coefficient is nonzero is the generic expectation, and
the predictions test it.

## Method (fixed before scoring)

- **Driver:** `debug/ellipse_ring_study.py`, sha256
  `5ffed8b7b47bfd39c8a7ba03fce408f28846975dfdb8b327d11b54c0eb78cf02`. It loads `debug/hring_study.py` (sha256
  unchanged, `afd9102800943aa8c45e4c6718209aad5f4823d9ff45f592e9c7f4cc97534f28`) for the basis and integral helpers.
- **Basis:** one Slater 1s per H, zeta 1.24, fitted by 3 Gaussians; Lowdin site orbitals; complex DFT over j.
- **Measure:** for each class c = |delta k| mod N (c = 1 ... floor(N/2)), the share of the Frobenius norm of h and
  of (pq|rs) carried by entries in that class.
- **Slopes:** local log-log slope of a class share against eps between consecutive eps values.
- **Scored cases** (`python debug/ellipse_ring_study.py case N 1.0 0 0.005 0.01 0.02 0.04 0.08 0.16 --out study`;
  rows to `debug/data/ellipse_ring_study.jsonl`): spacing **s = 1.0 A** (the controls used 1.2), N = 9, 10, 11, 12,
  13. N = 10, 12, 13 have no control at any spacing. Closest H-H distance at eps = 0.16 is about 0.84 A.
- **Floor:** a slope is scored only when both shares of its step are at least 1e-13; steps below the floor are
  reported and not scored.

## What I know before comparing (disclosed)

Five non-scored control rings at s = 1.2 A, N = 6, 7, 8, 9, 11, each at eps = 0, 0.005, 0.01, 0.02, 0.04, 0.08,
0.16: `debug/data/ellipse_ring_controls.jsonl`, 35 rows, sha256
`9a8bbd1ea10c3b820a400440a122955947068fe2f7ebdbeea59f294fdd529acb`. The order rule was derived before the controls.

| Measure | Control result |
|:--|:--|
| Every class share at eps = 0 | at most 1.8e-15 |
| Forbidden classes (odd c at even N), every eps | at most 1.6e-15 |
| Wrap-around orders seen | N = 7, c = 1: 3.00; N = 9, c = 1: 4.00, c = 3: 3.00; N = 11, c = 1: 5.00, c = 3: 4.00, c = 5: 3.00 |
| Largest deviation of a slope from p(c), step 0.005 -> 0.01 | 0.007 |
| Largest deviation, any step ending at eps <= 0.08 | 0.050 (N = 11, h, c = 1, step 0.04 -> 0.08) |
| Largest deviation, step 0.08 -> 0.16 | 0.097 |
| Smallest allowed share at eps = 0.005 | 7.5e-13 (N = 11, g, c = 1) |

## Predictions

- **SQ1 (exact zeros).** Every class share at eps = 0, and every odd-c share at even N (10, 12) at every eps, is at
  most 1e-13, for h and g.
- **SQ2 (order at small squash).** For every allowed class at N = 9-13, h and g, the slope of the step
  0.005 -> 0.01 is within 0.02 of p(c), wherever the step is above the floor.
- **SQ3 (order over the range).** For the same classes, every slope of a step ending at eps <= 0.08 is within 0.10
  of p(c), wherever the step is above the floor.
- **SQ4 (wrap-around).** In particular, at odd N the c = 1 class has slope near (N - 1)/2: 4 at N = 9, 5 at N = 11,
  6 at N = 13, each within the SQ3 bar wherever scored. (A subset of SQ2-SQ3, listed because it is the sharp case;
  N = 13 is the first order never seen in a control.)
- Steps 0.08 -> 0.16 are reported and not scored.

## Novelty status (set before computing)

- The mechanism is standard. Treating a symmetry-breaking parameter as a "spurion" that carries a fixed charge, and
  reading the order of a breaking effect off that charge, is the spurion method of particle physics (for example
  the quark mass in chiral perturbation theory). Momentum conserved only modulo the lattice size is Umklapp. The
  shape-based instances in molecules are Renner (bending, 1934) and Encinosa-Williamson (torus curvature; TR).
- SQ is a check of that counting in a molecule, and a third shape instance for the atlas reading. It is not new
  physics. No scan was run for this specific squashed-ring computation.

## Controls

The five control rings above, at a different spacing, set every bar. No scored case was run before this record was
committed.

## Decision

- **GO:** SQ1, SQ2 and SQ3 held.
- **BORDERLINE:** SQ1 held; some SQ2 or SQ3 slope missed its bar but every scored slope is within 0.5 of p(c).
- **STOP:** SQ1 missed, or any scored slope is 0.5 or more from p(c) (the order rule is wrong for that class).

## Reading if GO, stated in advance

On a squashed ring the ring momentum dies as eps^{p(c)}, with p(c) counted from the squash's charge of 2,
including the wrap-around order at odd N. With the bending and torus cases, that makes three shape deformations
that break a label as a power fixed by charge counting, against one distance door (a second nucleus) that breaks it
through an exponential tail. That supports writing the reading into the atlas as a dated observation. It says
nothing about energies, and nothing about how large the leading coefficients are.

## Disposition (2026-11-02, appended after scoring; predictions above unedited)

**Decision: GO.** SQ1, SQ2 and SQ3 held; SQ4 held. Scored rows `debug/data/ellipse_ring_study.jsonl` (35 rows,
s = 1.0 A, N = 9-13; sha256 `55063b4932000fba75c242346c6e64d39cc1029c63d906f5875e8dcb1b224a7e`), driver hash as
registered.

| Prediction | Result | Verdict |
|:--|:--|:--|
| SQ1 exact zeros <= 1e-13 | max 2.1e-15 (every class at eps = 0; odd c at N = 10, 12, every eps) | HELD |
| SQ2 first-step slope within 0.02 of p(c) | max deviation 0.0056 | HELD |
| SQ3 slopes to eps = 0.08 within 0.10 of p(c) | max deviation 0.045 (N = 13, h, c = 1, step 0.04 -> 0.08) | HELD |
| SQ4 wrap-around, c = 1 at odd N | N = 9: 4.001-4.033; N = 11: 5.001-5.041; N = 13: 5.999-6.045 | HELD |

Below the floor, not scored: N = 13, c = 1, step 0.005 -> 0.01 (shares 6.0e-15 and 7.8e-15 at eps = 0.005, at
the rounding level; slopes 5.86 and 5.98 shown for completeness). Steps 0.08 -> 0.16, reported only: every one
within 0.08 of p(c). N = 13's c = 1 (order 6) and c = 3 (order 5) are orders no control reached at that class.

**Reading.** As stated in advance: the ring momentum dies as eps^{p(c)}, counted from the squash's charge of 2,
with the wrap-around order at odd N. Three shape deformations (bending, torus curvature, squashing) now break a
label as a power fixed by charge counting; one distance door (a second nucleus) breaks it through an exponential
tail. This is the known spurion and Umklapp counting, checked in a molecule.

**Caps.** One minimal Slater 1s per H, no energies; planar rings only; the coefficients of the leading powers are
measured, not predicted; the bending instance (door 15) measures an energy splitting, while TR and SQ measure
matrix-element shares, so the three agree on the counting, not on one quantity.
