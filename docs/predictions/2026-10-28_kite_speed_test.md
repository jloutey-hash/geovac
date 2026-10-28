# Pre-registered prediction: the kite operator against a float64 quadrature at chemistry precision (KS1-KS5)

**Registered:** 2026-10-28, before any timing in this file. Draft by a sub-agent; reviewed and committed by the PM.
**Prompted by:** record OE (`2026-10-26_operator_evaluator.md`), whose speed ratio was measured at 60 digits against
this program's own 60-digit grid, and record MCT (`2026-10-27_master_count_tabulation.md`), which says a timing test
against a production bar is "a separate later step with its own record".
**Extends:** `2026-10-26_operator_evaluator.md`.
**Commit status:** to be committed on its own, before any timing it scores.

## The question

Record OE found that the kite's order-11 operator rebuilds the kite at the 247 nodes of a 60-digit grid far faster
than the 60-digit quadrature that built the grid. Chemistry needs about 1e-9 to 1e-12, not 60 digits. At that
precision the natural rival is plain float64 quadrature, not 60-digit ball arithmetic. Does the operator route still
pay there, on one CPU core, with and without the cost of deriving the operator?

## What the kite route computes (read from the records, not assumed)

- The object is Y(t) = <rho_i rho_j rho_k exp(-gamma_1 r_13) exp(-gamma_2 r_23) / r_12>, the exponential-geminal
  correlated triangle, at (a, b, gamma_1, gamma_2) = (1, 2, 3/2, 3), apex scale t (`debug/pf_common.py` docstring).
- The densities are the "fingerprint fibre" densities of Paper 59: rho~(k) = (c k^2 + 1)^(-3/2), so rho is
  proportional to K0(a r) in position space (`debug/sprint_ring_ext_QG_memo.md`). They are not 1s densities. Y is one
  point of a fibre, an integrand of an outer integral that has not been done. It is not a molecular integral.
- Y is computed in two stages. Stage 1 is the kite phi(z), z = w^2, a 2-D Feynman-parameter integral at each node of
  an exp-sinh w-grid. It does not depend on t. Stage 2 is one weighted sum over the nodes per t. Stage 1 is the cost.
- The operator is exact for this one parameter set. A new (a, b, gamma_1, gamma_2) needs a new derivation. Record PK
  (`2026-10-26_parametric_kite_operator.md`) stopped a one-parameter family on sample count. So the derivation is a
  per-parameter-set cost, amortized only over t and over the nodes, never over parameter sets.

## Why there is no fair comparison with Chen-Silberman-Segal-Greif

Read from the primary source, arXiv:2610.12608v1 (math.NA, 8 Oct 2026), downloaded to a folder outside the repo:

- It times two-electron four-centre Coulomb repulsion integrals over Slater-type orbitals. Accuracy: ten reference
  integrals (Table 1), largest absolute discrepancy 1.7e-9 against a FORTRAN reference code; values are 0.003 to 0.31.
- Two rates. "A preliminary resident-cache benchmark computes approximately 4.2 x 10^5 unique integrals per second on
  an NVIDIA V100 at the contraction stage" (abstract; Sec. 6 adds "not an end-to-end timing including cache
  generation"). And, end to end over whole molecules, "we sustain about 40,000 ERIs per second end-to-end" (Sec. 6,
  Fig. 1: 78 to about 50,000 ERIs/s by workload, rising with molecule size).
- The rate is batch throughput: pair coefficients are cached and reused across up to about 2e7 integrals per molecule.

Y is none of these things. It has three electrons, two geminals, K0 fibre densities, one centre, and no orbital-pair
reuse. No known identity writes a four-centre STO ERI as a sum of Y values, or the reverse. Any ratio of our time per
Y to their time per ERI would compare different integrals on different hardware in different batch regimes. **No
comparison with arXiv:2610.12608 is made in this record.**

What would make one possible: a GeoVac route to the same four-centre STO ERIs (for example Paper 59's two-centre
engine extended to four centres), run on their ten Table-1 integrals and centres at 1e-10 absolute, timed per integral
and per molecule batch on one CPU core, with a GPU-to-core factor stated as a range. The V100's FP64 peak (about 7 to
8 TFLOP/s) against one Zen 3 core (tens of GFLOP/s) gives a factor of order 1e2; this is recalled, not verified, and a
real normalization would need a measured kernel on both.

Note for the PM: Paper 59 (line ~162) quotes "about 4 x 10^5 integrals per second on one GPU", and record MCT quotes
"about 4e4 ... end to end". Both are in the source. They are different rates (contraction stage, end to end). Paper 59
does not say its figure excludes the cache build.

## The baseline chosen, and why

**Route Q (baseline): float64 tanh-sinh quadrature of the same route-B double integral, numpy-vectorized.** In u, a
tanh-sinh rule graded at both ends (u and 1-u kept separately). In v, the same rule on two panels split at the ridge
v = 1 - u, where k2 = (1-u-v)^2 z vanishes and the kernel peaks at large z. This is the standard numerical route a
production code would use for a smooth 2-D integral at 1e-12, and it runs at the hardware's native precision.

Rejected candidates, with reasons:
- Gaussian-fitted geminals with closed-form Gaussian integrals (the F12 route). Records GF and PG show six-term fits
  miss Y by 2e-6 to 13%. They do not reach the target precision, so the comparison would not be at matched accuracy.
  The K0 fibre densities are also not Gaussians.
- The 60-digit arb tanh-sinh grid of record OE. It is the baseline record OE already beat, and it spends 60 digits
  that chemistry does not need. Using it here would repeat OE, not test it.
- Chen et al. See above.

## What I know before timing (disclosed)

- Record OE, at 100 working digits: route O kite 62.7 s wall on one process (248 Taylor steps), input coefficients
  0.96 s; the 60-digit grid 7893.7 s wall on 15 processes. Derivation: one prime re-measured at 3,065 s of processor
  time; three primes about 9,200 s, assuming p10 and p13 cost what p5 did; the lift 12.8 s. These numbers are not
  re-scored here. The derivation total enters KS4 as a known constant.
- Feasibility checks run for this draft, accuracy only, no clock read (scripts below):
  - `debug/ks_check.py`: a tensor Gauss rule without grading (Chebyshev in u, Legendre in v split at the ridge) gives
    phi to about 1e-15 for w up to 6, but fails at w = 40 (3e-3 at 80 x 80 nodes). Grading is needed.
  - `debug/ks_check2.py`: the graded float64 tanh-sinh rule above. Node errors at inner step h = 1/16: 3e-16 at
    w = 4.4, 4e-13 at w = 93, 4e-5 at w = 1.4e5, 0.5 at w = 1.6e13. Weighted into Y over the 247-node grid, the
    worst single-node contribution is 1e-14 of Y at h = 1/8 and 9e-18 at h = 1/16. So the large-w nodes where float64
    fails do not matter at 1e-12, as they did not matter in record OE.
  - `debug/ks_check3.py`: the outer w-grid can be coarser. From the stored 60-digit phi, step 1/8 with a cut at
    1e-16 of Y keeps 40 nodes and gives Y at t = 0.8, 1.3, 1.7, 2.5 to 6.5e-14, 1.2e-13, 1.9e-13, 3.8e-13. Step 1/4
    gives 4e-7 to 1.4e-5, too coarse for 1e-9.
- So route Q's accuracy at 1e-12 is essentially known. KS1 is therefore a control, not a prediction.
- Nothing about route O below 100 digits has been run: neither its accuracy nor its time. Whether the walk survives at
  20 to 40 digits, given that rounding near apparent roots broke it once in record OE, is not known.
- No time of route Q has been read.
- The machine: AMD Ryzen 7 5825U, 8 cores, 16 threads. Another agent is computing at the same time.

## Model

**Common to both routes.**
- One parameter set, (1, 2, 3/2, 3). Targets T12 = 1e-12 and T9 = 1e-9, relative, on Y at t = 0.8, 1.3, 1.7, 2.5.
- Reference: the stored 60-digit grid's Y (`debug/data/pf_grid_d60_h0.03125.json`), good to about 1e-60.
- Outer grid, fixed before timing and shared: exp-sinh step 1/8, nodes kept where |term| > 1e-16 |Y| (T12) or
  1e-13 |Y| (T9) at any of the four t. These are chosen from the stored reference, before any route runs.
- Stage 2, the four weighted sums, is the same code for both routes. It is timed and reported but cannot differ.
- One "evaluation" = phi at every node of the shared outer grid, for one parameter set, then Y at the four t.

**Route Q (baseline).** `debug/ks_route_q.py`, float64 numpy. Inner tanh-sinh step from the fixed ladder
h = 1/4, 1/8, 1/16, 1/32 (kmax = 5/h). The scored setting is the first rung, in that order, that meets the target at
all four t. The rung is chosen on accuracy alone, before any time is read.

**Route O (operator).** `debug/ks_route_o.py`, reusing `debug/oe_engine.py` and the walker of `debug/oe_kite.py`
unchanged, on the shared outer grid only.
- Input: phi_0 to phi_5 by `debug/pf_4b_series_direct.py` on one process, at the working precision plus 10 digits.
- Working precision from the fixed ladder 20, 25, 30, 40, 60, 100 digits. The scored setting is the first rung that
  meets the target at all four t, chosen on accuracy alone.
- Two costs, reported apart: the one-off load (reading `debug/data/ib_pf_op_Q.json` and building the polynomial
  coefficients) and the walk (input plus Taylor steps). The walk plus the load is route O's per-parameter-set cost
  when the operator is given.
- A complex128 numpy port of the walk may be attempted, with a cap of 1 hour of agent effort. If it meets the target
  it is scored alongside the flint form, and route O takes the faster of the two. If it does not, it is reported and
  not scored.
- Derivation: 9,200 s of processor time from record OE's addendum, not re-run here (known constant, disclosed).

**Timing protocol.**
- Processor time, `time.process_time()`, inside one process, with OMP_NUM_THREADS, MKL_NUM_THREADS and
  OPENBLAS_NUM_THREADS set to 1 before numpy loads. Python start-up and imports are excluded and reported.
- Five repeats per route per target, run interleaved (O, Q, O, Q, ...), median scored, minimum and maximum reported.
- Ratios: R_marg = route O (load + walk + stage 2) / route Q (stage 1 + stage 2). R_tot = (route O + 9,200 s) /
  route Q. A ratio above 1 means route Q is faster.

## Controls (run first)

- **C1 (baseline accuracy).** Route Q at its scored rung meets the target at all four t. Largely known (above).
- **C2 (operator accuracy).** Route O at its scored rung meets the target at all four t. If no rung up to 100 digits
  does, route O is void at that target, reported as a feasibility STOP, and KS2 fails.
- **C3 (agreement).** At the scored rungs, routes O and Q agree with each other on Y to within the target at all four
  t, and with the reference.
- **C4 (negative, operator).** Route O with one coefficient of P_5 times 201/200, as in record OE, misses Y by more
  than 1e-6 at every t.
- **C5 (negative, baseline).** Route Q with the slice assignment k2 = (u-v)^2 z (`assign='wrong'` in pf_common)
  misses Y by more than 1e-6 at every t.
- **C6 (timer, positive).** Route Q at inner step h and h/2: the time ratio lies between 2.5 and 6 (node count rises
  four-fold). This checks that the timer resolves the work.
- **C7 (reproducibility).** For each route, max/min of the five repeats is at most 1.3. If not, run up to ten repeats
  and report; if still not, the timings are UNCLEAR and the gate is not scored.

## Predictions (scored at T12; T9 reported with the same tests)

**KS2 (route O survives low precision).** Route O meets T12 at all four t at some rung of 40 digits or fewer.
- Falsifier: it needs 60 digits or more, or no rung meets it.
- Probability 0.7. The risk is the error growth near apparent roots seen in record OE.

**KS3 (the baseline wins when the operator is given).** R_marg > 1 at T12.
- Falsifier: R_marg <= 1.
- Probability 0.85. Route Q needs about 40 nodes times about 2 x 81^2 vectorized float64 kernel calls. Route O needs
  over a hundred Taylor steps, each an order-11 recurrence with coefficients of degree up to 113, in multiprecision.

**KS3b (by a wide margin).** R_marg > 10 at T12.
- Falsifier: R_marg <= 10.
- Probability 0.6.

**KS4 (the derivation decides).** R_tot > 1000 at T12.
- Falsifier: R_tot <= 1000, that is, route Q costs more than about 9 s of processor time.
- Probability 0.95.

**KS5 (the verdict does not depend on the target).** The gate outcome at T9 is the same as at T12.
- Falsifier: different outcomes.
- Probability 0.85.

## Decision gate (scored at T12, controls C1 to C7 passed)

- **GO:** R_tot < 1. The operator route, derivation included, beats float64 quadrature at chemistry precision.
  Probability 0.01.
- **BORDERLINE:** R_tot >= 1 and R_marg < 1. Once known, the operator evaluates faster, but deriving it costs more
  than it saves for one parameter set. It would pay only if one operator served many parameter sets, which record PK
  did not reach. Probability 0.14.
- **STOP:** R_marg >= 1. At chemistry precision the operator route does not pay even with the operator given; record
  OE's ratio was a property of 60-digit arithmetic. Probability 0.85.
- If C2 fails at T12 (no rung reaches it), the outcome is STOP for feasibility, reported, and KS3 to KS5 unscored.

## Novelty status

- Solving a Feynman integral's differential equation by series as an evaluator: KNOWN (Moriello, arXiv:1907.13234;
  Hidding, arXiv:2006.05510; Liu and Ma, arXiv:2201.11669; checked in record OE). Known too, as general practice, that
  such solvers pay at high precision and against costly multi-loop numerics, rather than against cheap low-dimensional
  float64 quadrature. That is recalled, not read from a source here.
- KS2 to KS5 on this family: NOT FOUND, and not searched beyond record OE's brief search. Any result is at most a
  measured cost on one integral, a CHECK of expectations if STOP.

## Caps stated in advance

- One parameter set, one machine, one core, one integral. No production batch, no orbital-pair reuse.
- The derivation cost is from record OE (one prime measured, two assumed equal). It is not re-measured here.
- The baseline is plain numpy. A compiled or batched baseline would be faster still; a compiled route O might narrow
  the gap. Neither is tried beyond the optional complex128 port.
- Another agent shares the machine. Processor time limits, but does not remove, contention effects (cache, clocks).
- Accuracy is checked against the 60-digit reference, not certified.
- Nothing here bears on Chen et al., on physical three- or four-centre ERIs, or on any molecular energy.

## Phase-2 caps and files

- At most 4 of my processes at once; processes I did not start are never touched. Every script has a main guard.
- At most 2 hours of wall time for phase 2. Exceeding it is a STOP, reported with what was reached.
- New files only: `debug/ks_route_q.py`, `debug/ks_route_o.py`, `debug/ks_time.py`, `debug/data/ks_*`. Existing
  scripts are imported, never edited. Phase-1 checks: `debug/ks_check.py`, `debug/ks_check2.py`, `debug/ks_check3.py`.
- Disposition appended below after scoring; the text above is not edited.

## Disposition (appended 2026-10-28; predictions above unchanged)

**Gate at T12: STOP.** R_marg is about 500: with the operator given, float64 quadrature is about 500 times faster.
KS2 FAILED, KS3 HELD, KS3b HELD, KS4 HELD, KS5 UNCLEAR (not scored as registered, because C7 failed at T9).
Scripts `debug/ks_common.py`, `debug/ks_route_q.py`, `debug/ks_route_o.py`, `debug/ks_time.py`,
`debug/ks_time_q_loop.py`; data `debug/data/ks_outer_{T12,T9}.json`, `debug/data/ks_timing.{json,log}`,
`debug/data/ks_timing_qloop.{json,log}`. One child process at a time; phase 2 took well under the 2 h cap.

**Shared outer grid.** T12: 40 nodes, w from 4.9e-6 to 653. T9: 36 nodes, w from 7.5e-5 to 149. With the exact
60-digit phi, stage 2 on these grids gives Y to 6.5e-14 to 3.8e-13 (T12) and 1.2e-13 to 4.6e-13 (T9).

**Ladders (accuracy only, chosen before any clock was read).**
- Route Q, max relative error on Y over the four t: h = 1/4: 3.7e-9 (fails both targets); h = 1/8: 4.6e-13 (T12),
  5.3e-13 (T9). Scored rung h = 1/8 at both targets.
- Route O at T12: 20, 25 and 30 digits give nan or inf; 40 digits gives 1.8e278; 60 digits gives 3.8e-13. At T9:
  20 to 40 digits blow up (inf, inf, 8.7e263, 4.8e138); 60 digits gives 4.6e-13. Scored rung 60 digits at both.
  The failure below 60 digits is a blow-up, not a gradual loss. Its cause was not diagnosed.
- Route O used 76 main-path steps (80 in all, 2 branch walks) at T12, 69 (73) at T9, against 248 in record OE.

**Controls.**
- C1 PASS (above). C2 PASS, at 60 digits only. C3 PASS: O and Q agree on Y to 1.4e-14 to 7.4e-14 at all four t, and
  both meet the target against the reference.
- C4 PASS: perturbed operator, Y wrong by 1.4e-2 to 6.8e-2 (both targets). C5 PASS: slice assignment, Y wrong by
  6.9e-2 to 1.1e-1.
- C6 PASS: Q at h = 1/16 over h = 1/8, median ratio 5.0 (registered run); 4.94 in the looped run below.
- C7: T12 PASS (O 1.14, Q 1.00). T9 FAIL: after ten repeats Q's max/min was 1.5 (O 1.10). As registered, the T9
  timings are UNCLEAR and the T9 gate is not scored.

**Timer resolution (not foreseen; disclosed).** On this machine `time.process_time()` ticks at 15.625 ms. Route Q's
whole evaluation is 2 to 3 ticks: every registered Q repeat read 0.03125 s, one read 0.046875 s. That is the T9 C7
failure. As a supplement, not part of the registered scoring, Q was repeated 50 times inside each of five children
(`debug/ks_time_q_loop.py`): median 33.1 ms (T12, spread 1.08), 163.8 ms at h = 1/16 (spread 1.13), 30.9 ms (T9,
spread 1.07).

**Timings (processor time, one process, medians).**

| | T12 | T9 |
|:--|:--|:--|
| Route O load (operator, roots) | 0.08 s | 0.08 s |
| Route O input (six coefficients, 70 digits) | 0.22 s | 0.23 s |
| Route O walk (60 digits) | 16.08 s | 16.02 s |
| Route O total, spread | 16.39 s, 1.14 | 16.34 s, 1.10 |
| Route Q total, registered (tick-limited) | 0.031 s, 1.00 | 0.031 s, 1.50 |
| Route Q total, looped supplement | 0.0331 s, 1.08 | 0.0309 s, 1.07 |
| R_marg, registered / looped | 525 / 495 | 523 / 528 |
| R_tot (with 9,200 s derivation), registered / looped | 2.9e5 / 2.8e5 | 2.9e5 / 3.0e5 |

Imports before each run: about 0.30 s (O) and 0.23 s (Q), excluded as registered. Stage 2 was below one tick for
both routes.

**Scores.**
- **KS2 FAILED.** Route O needs 60 digits at both targets; nothing at 40 or fewer reaches either. **MISS.**
- **KS3 HELD.** R_marg = 525 at T12. **UNCLEAR** (novelty not searched at scoring).
- **KS3b HELD.** R_marg > 10 by a factor of about 50. **UNCLEAR**, same reason.
- **KS4 HELD.** R_tot = 2.9e5; route Q costs 0.03 s, not 9 s. **UNCLEAR**, same reason.
- **KS5 UNCLEAR.** Not scorable as registered (C7 failed at T9). In the looped supplement the T9 outcome is the same
  STOP (R_marg 528).

**Departures from the model, all disclosed.**
- Route O's detour around a real apparent root is walked only if 0.65 rho is below the largest node's z. The shorter
  grids end before the later roots, and walking past the last node would only add steps. Otherwise oe_engine and the
  oe_kite walker were imported unchanged, including its 1200-bit floor for the local operator coefficients.
- The input coefficients were computed in-process by `pf_4b_series_direct.u_chunk` (80 x 80 'sub', one process),
  not by a subprocess as in record OE, so that process_time counts them.
- The complex128 port of the walk was not attempted. The flint walk fails at 40 digits, so a 16-digit port cannot
  reach either target.
- The looped Q timing is an addition to the plan, forced by the timer tick. The registered scoring stands without it.

**Reading.** The 1/13 of record OE is a property of 60-digit arithmetic, where the rival was 60-digit ball
quadrature. At chemistry precision the kite's integrand is a smooth 2-D integral, and a graded float64 rule gets Y to
5e-13 in about 33 ms. The operator walk cannot drop below 60 digits on this path, and at 60 digits it costs about
16 s, nearly all of it in the walk. Counting the derivation, which is owed afresh for every parameter set, the
quadrature is about 3e5 times cheaper. This says nothing about Chen et al., whose integrals are a different object.
