# Pre-registered prediction: the integral's own differential equation as an evaluator, timed against quadrature (OE1-OE4)

**Registered:** 2026-10-26, before any computation in this file. PI direction: "Doesn't this research into the
structure of these integrals help us run them optimally", then, on the two proposed tests, "Let's do both."
**Prompted by:** records IB (`2026-10-24_ibp_kite.md`) and PF (`2026-10-23_price_formula_slots.md`). The operators
found there were used only to check values. This record asks whether they are cheaper evaluators.
**Extends:** `2026-10-24_ibp_kite.md`, `2026-10-23_price_formula_slots.md`.
**Commit status:** committed on its own, before any computation in this file.

## The question

Y(t) off the symmetric slice, at (a, b, gamma_1, gamma_2) = (1, 2, 3/2, 3), is computed today in two stages:

1. the kite Phi~(w) at 247 nodes of an exp-sinh grid, w from 1.2e-16 to 8.4e15, each node by adaptive 2-D
   tanh-sinh quadrature at 60 digits (`debug/pf_3_grid.py`);
2. Y(t) and its theta-derivatives as one 247-term weighted sum per t (`ring_ext_QG_3_values.Y_theta`).

Stage 1 is the cost. The kite phi(z), z = w^2, obeys an exact order-11 operator over Q(z) found by integration by
parts (`debug/data/ib_pf_op_Q.json`). Stage 2 has its own exact order-13 operator M13 in t
(`debug/data/pf_34_m13_exact.json`). The question: does solving these equations beat the quadrature, at equal
accuracy, and where?

## What I know before computing (disclosed)

- The grid log records 7893.7 s wall on 15 worker processes for the 247 nodes
  (`debug/data/pf_grid_d60_h0.03125.log`). CPU time was not logged. Tasks are handed out one node at a time, so the
  workers were roughly balanced, but 15 x 7893.7 s is an upper bound on its CPU time, not a measurement.
- The IBP operator came from 80 samples per prime, 685-687 s single-process per prime
  (`debug/data/ib_pf_p0_run.log`, `ib_pf_p1_run.log`), at two primes for the lift and a third for the check, plus
  the lift itself, whose time was not logged. I will time a re-run of the lift step if it is cheap, and disclose
  it either way.
- Exponents at z = 0, computed today by `debug/ib_exponents_Q.py` before this record, to plan the method: -4, -1/2,
  0, 1/2, 1, 3/2, 2, 5/2, 3, 4, 5. phi is analytic at 0, so it lies in the span of the solutions with exponents
  0, 1, 2, 3, 4, 5, and its first six Taylor coefficients fix it if those six solutions are log-free. Whether they
  are is not checked.
- All genuine singular points of the kite operator lie on the negative real z axis (record IB, IB4: -1/4 to
  -81/4). The positive axis, where every node lies, meets only apparent singular points, if any. Where they sit is
  not checked.
- The Taylor coefficients at z = 0 are known to 280 digits from a costly run (`debug/pf_4b_series_direct.py`, 300
  coefficients). Those are not used as input. Only the first six are needed, and they will be recomputed by the same
  script at the needed precision and timed.
- Nothing about the operator route's values or times has been computed.

## Model

**Route O for stage 1, the kite.**
- Input: phi_0 to phi_5, the first six Taylor coefficients at z = 0, by `debug/pf_4b_series_direct.py` at about 80
  digits. Timed.
- The Taylor series of phi at 0 from the operator's recurrence, used for |z| up to about 1/8, inside the radius 1/4
  set by the operator's nearest singular point.
- Analytic continuation from there to every node, through the right half-plane: Taylor steps at rational centres,
  each step at most half the distance to the nearest root of the leading coefficient, with the series of phi at
  each new centre built from the operator. Complex detours only if an apparent singular point lies on or near the
  positive axis.
- Large nodes: continue as far as the stepping keeps the digits needed. If error growth stops it, the remaining
  nodes, where w^2 Phi~ is smallest, are reported separately, and route O is scored on what it reaches.
- Output: Phi~ at the 247 nodes, then stage 2 unchanged.
- Timing: single process, wall time, on this machine.

**Route O for stage 2, Y's own operator, marginal cost.** From theta^j Y, j = 0 to 12, at t0 = 1.7 on the grid,
continue M13 in t to t = 0.8, 1.3, 2.5, 2.6, 3.8, avoiding its singular points. Time it per t, against the time of
one grid sum per t.

## Predictions

**OE1 (accuracy).** The route-O kite, put through stage 2, gives Y at t = 0.8, 1.3, 1.7 and 2.5 that agrees with the
grid's Y to 1e-40 relative at all four t.
- Falsifier: worse than 1e-40 at any of the four.
- Confidence: moderate. The risk is error growth along the path to large z.

**OE2 (cost, the main claim).** Route O for the kite, input included, takes at most 789 s wall on one process. That
is a tenth of the grid's wall time on fifteen.
- Falsifier: more than 789 s.
- Confidence: moderate. Six input coefficients at low degree should be cheap. The unknown is the number of
  continuation steps and their cost at 80 digits with coefficients of degree up to 113.

**OE3 (cost, derivation included).** Route O plus the IBP derivation (the three logged prime runs, about 2060 s,
plus the lift) takes less than a tenth of the grid's CPU bound, 15 x 7893.7 s.
- Falsifier: a tenth or more.
- Confidence: moderate to high, if OE2 holds.

**OE4 (where the structure does not pay).** Y's own operator M13 reproduces the grid's Y at t = 0.8, 1.3, 2.5, 2.6
and 3.8 to 1e-35, but its cost per t is higher than one grid sum per t.
- Falsifier: worse than 1e-35, or route O cheaper per t than the grid sum.
- Confidence: moderate. The grid sum per t is 247 weighted terms. The operator has order 13 and degree 75.

## Novelty status

- Evaluating Feynman integrals by numerically solving their differential equations, with analytic continuation by
  series: KNOWN in particle physics. Recalled, not read here: DiffExp (Hidding, Comput. Phys. Commun. 269, 108125,
  2021), AMFlow (Liu and Ma, Comput. Phys. Commun. 283, 108565, 2023), and the series method of Moriello
  (JHEP 01 (2020) 150). To be checked against the sources before any disposition cites them. OE1 to OE3 are checks
  of the method on our family if they hold.
- For explicitly correlated molecular integrals of this family, with half-integer density lines and squared geminal
  lines: NOT FOUND in a brief search, not to the standard of a novelty claim. OE2 and OE3, as a measured cost ratio
  on such an integral, would be at most a modest finding.
- OE4: no prior statement found or sought. It is a property of this two-stage computation.

## Controls

- Node check: route O's Phi~ against the grid's phi at every node it reaches, relative error reported node by
  node. The grid's own accuracy is about 1e-60.
- Input check: the six recomputed coefficients against the 280-digit series, to the precision used.
- Recurrence check: the Taylor coefficients 6 to 40 generated by the operator, against the 280-digit series.
- Negative control: the same continuation with the operator perturbed, one polynomial coefficient times 201/200,
  must fail the node check by many orders of magnitude.
- Stage 2: M13 perturbed the same way must fail OE4's accuracy band.

## Caps stated in advance

- One parameter set. The operator is for (1, 2; 3/2, 3) only. A new parameter set needs a new IBP derivation, which
  is why OE3 includes it.
- Wall time on one machine, at different process counts. The grid's CPU time is an upper bound.
- Not rigorous. No error bounds are certified, only checked against the grid.
- A production code needs the integrals for many orbital pairs at once. This measures one integral.
- If the input coefficients cannot be had cheaply, or the path cannot reach the nodes, that is a STOP for
  feasibility, reported and not scored.
