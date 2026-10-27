# Pre-registered prediction: how much of the LiH explicitly correlated run is the triangle class (LT1-LT7)

**Registered:** 2026-10-26, before any instrumented run. PI direction: find out whether the triangle class is worth
targeting in the real LiH calculation, i.e. whether it is a large share of LiH's cost.
**Prompted by:** `2026-10-26_operator_evaluator.md` (OE1-OE4), where a one-centre model of the triangle class was
evaluated from its exact differential equation about 40x cheaper than by quadrature (OE3, ratio 1/41).
**Extends:** `2026-10-26_operator_evaluator.md`.
**Commit status:** committed on its own, before the measurement (phase 2).

## The question

The corpus-best LiH energy (E_R12 = -8.049045088639625 Ha, exp geminal with gamma = 1.5, pi reference at truncation
1e-3) is assembled from pieces: Fbar and sigma2, h_Vne and g_Vne, h_T and g_T, then g_Vee. Inside g_Vee one term
class is a three-body "triangle": two geminals share a vertex and the Coulomb line sits on the opposite edge
(`geovac/lih_r12ci/triangle.py`). Two questions:

1. What share of the run's wall time goes to the triangle class?
2. How many distinct triangle kernels K_mu(1,2) = INT mu(3) f(r13) f(r23) d3 does the run build?

## Configurations

- **Target (corpus-best):** golden item `pi_marriage_exp1.5_1e-3` in `tests/data/lih_golden/golden.json`.
  nact 196, E_R12 = -8.049045088639625. Stored cost: 10,538 s and 8.44 GB peak (golden run). The original
  2026-10-03 run logged 11,123 s and 15.1 GB peak with older code (`debug/data/lih_exp_geminal_exp1.5_1e-3.log`).
  Single process, one BLAS thread.
- **Reduced proxy:** golden item `pi_marriage_exp1.5_1e-2`. Same code paths, same geminal, same grid (72x44, 28
  azimuthal points), same triangle code. Only the det truncation differs, so nact is 64 instead of 196. Stored cost:
  2,008 s and 2.81 GB. Golden E_R12 = -8.048140609752963.
- Both are driven by `debug/lih_triangle_share.py`, which mirrors the golden call and wraps timers around existing
  functions by monkeypatching. `geovac/` is not edited.

## What I know before computing (disclosed)

I read the code and the 1e-3 log. The log is a record of an earlier run of the same calculation, so the phase times
below are not blind. I disclose them and say which predictions lean on them.

- Cumulative stamps in the 1e-3 log: Fbar and sigma2 done at 905.6 s; h_Vne and g_Vne at 4,468.8 s; h_T and g_T at
  7,855.3 s; g_Vee non-triangle at 8,902.1 s; triangle phase done at 11,122.7 s. So the triangle phase took 2,220.6 s
  of 11,122.8 s, which is 20.0 percent. The h_Vne/g_Vne block took 3,563 s (32 percent) and the h_T/g_T block
  3,387 s (30 percent). These stamps include the time before 905.6 s (primitive load, FCI, natural orbitals).
- The same log reports `tri_solve_rows` = 2,412,864 and `tri_cm` = 196. That is about 0.92 ms per Coulomb row if the
  triangle phase were all solves.
- From the code (`gvee_sparse.tri_entries`): per middle density cm, the middle density has azimuthal number
  mu_mid. Combinations (A, sa, B, sb) with sa*A + sb*B = mu_mid and A, B in 0..4 give, per cm: mu_mid = 0: 9
  combinations, of which 5 distinct symmetric-kernel eigendecompositions (A = 0..4, the cache key ignores the
  signs); |mu_mid| = 1: 8 combinations, all A != B, each a randomized SVD that is not cached; |mu_mid| = 2: 7
  combinations, 1 eigendecomposition (A = B = 1) and 6 SVDs.
- The reference at 1e-3 has 196 ordered pairs, 88 of which carry m != 0 (`npi_pairs` in the stored npz). The
  natural orbitals have four pi orbitals (m = +1, -1, -1, +1). So |mu_mid| = 2 needs a pi pair with opposite signs,
  at most 8 of the ordered pairs. I did not build the reference to count the exact split.
- Before computing, I estimated the matvec cost of one K_mu factorization: about 7 GFlop (grid of 3,168 points, rank
  30 plus 16 oversample, power 2), so about 0.25 s on one thread.
- Phase 2 will run both configurations once. Nothing about the instrumented values is known.

## Model

Wall-clock timers around: `build_primitives`, `build_ref_from_prims`, `fbar`, `sigma2`, `hg_Vne`, `hg_T`,
`hg_Vee_sparse`, every `vee_reduce_sparse` call (the non-triangle g_Vee), `eval_triangles` and `tri_entries` (the
triangle phase), and inside it `coul_mode_stack_any` (the Coulomb mode solves), `Kmu_modes` and `_rand_svd_middle`
(the K_mu builders), and `assemble_c`. Counters: number of cm seen, K_mu eigendecompositions that miss the mode
cache, SVD calls.

**Definitions.**
- Triangle share = time in `eval_triangles` / total wall of the run (primitive load through `assemble_c`).
- The triangle collection step (recording deferred triangle terms) runs inside the non-triangle g_Vee loop and is
  not separable. It is counted as non-triangle. This biases the share low.
- "Distinct triangle kernels": two counts are reported. (a) K_mu middle densities mu, times geminals. The exp run has
  one geminal, so this is the number of cm. (b) Mode kernels actually factorized: eigendecompositions that miss
  the cache plus SVD calls.
- A cProfile run at 1e-2 is a cross-check on the wrapper totals only.

## Predictions

**LT1 (share at the target).** At truncation 1e-3 the triangle share of total wall time is 20 percent, band 15 to
27 percent.
- Falsifier: outside 15 to 27 percent.
- Confidence: high that it is inside the band, because the earlier log gave 20.0 percent. The risk is that the
  earlier run shared the machine with other jobs.

**LT2 (kernel counts).** At 1e-3 the run builds K_mu for 196 middle densities (one geminal). The number of mode
kernels factorized is 1,200 (band 1,000 to 1,400): about 540 symmetric eigendecompositions and about 690 SVDs.
At 1e-2 the run builds 64 K_mu, with 300 to 450 mode kernels.
- Falsifier: `tri_cm` other than 196 (64 at 1e-2), or a mode-kernel total outside the bands.
- Confidence: high for 196 and 64 (the log gives 196). Moderate for the mode-kernel bands, which rest on my
  reading of the combination rule and an unmeasured split of m values.

**LT3 (where the triangle time goes).** Inside the triangle phase at 1e-3, the Coulomb mode solves
(`coul_mode_stack_any`) take 75 to 90 percent of the triangle time. The K_mu factorizations (eigh plus SVD) take 5
to 20 percent.
- Falsifier: solves below 75 percent or factorizations above 20 percent.
- Confidence: moderate. The matvec estimate of 0.25 s per kernel and 0.92 ms per solve row both come from my
  reading, not a profile.

**LT4 (rank among the pieces).** At 1e-3 the triangle phase is not the largest piece. The h_Vne/g_Vne block and the
h_T/g_T block are each larger than it, and the non-triangle g_Vee is smaller than it.
- Falsifier: any other order.
- Confidence: moderate to high, because the log already shows this order. A different order would mean the earlier
  run was distorted.

**LT5 (proxy share).** At the reduced truncation 1e-2 the triangle share is 25 percent, band 10 to 40 percent. It is
not assumed to match the 1e-3 value.
- Falsifier: outside 10 to 40 percent.
- Confidence: low. Triangle cost depends on the number of distinct (cm, left, right) entries, which grows faster than
  nact (1,156 unique index triples at 1e-2 against 56,644 at 1e-3, from the build plan), while the Vne and T blocks
  grow with the number of kept dets. The direction of the difference is not known to me.

**LT6 (the replaceable part).** The part of the run an operator-based evaluator of the K_mu kernels could replace is
the K_mu factorization time, not the Coulomb solves. At 1e-3 that time is at most 5 percent of total wall.
- Falsifier: factorization time above 5 percent of total wall.
- Confidence: moderate. It follows from LT1 and LT3 (20 percent times 5 to 20 percent gives 1 to 4 percent).
- Consequence if it holds: even a free K_mu build saves at most about 1.05x on the whole run. A free triangle phase
  (the solves too) caps the speedup at 1/(1 - 0.20) = 1.25x (band 1.18x to 1.37x from LT1).

**LT7 (control).** The instrumented run at 1e-2 reproduces the golden E_R12 = -8.048140609752963 to every stored
digit, and at 1e-3 reproduces -8.049045088639625 (if the 1e-3 run is done). Timer overhead is below 1 percent of
wall, checked by the sum of the timed top-level pieces against the total.
- Falsifier: any difference in a stored digit of E_R12, c, h, or any stored piece; or overhead of 1 percent or more.
- Confidence: high. The wrappers only call the original and read a clock. A difference would mean the wrapping
  changed the call order or the random streams.

## Novelty status

This is a cost measurement of this repository's own code. No claim of novelty is made, and none sought. The
operator-evaluator method it is set against is the known one credited in `2026-10-26_operator_evaluator.md`
(Moriello 2020, Hidding 2021, Liu and Ma 2023).

## Controls

- Reproduction of the golden energies (LT7). The golden items are bit-identical in a fresh single-thread process, so
  the instrumented run uses the same process setup (one BLAS thread, `configure_geminal(1.5)` before imports,
  `USE_EXACT_NEUMANN = True`).
- Sum check: total wall against the sum of the top-level timed pieces plus the primitive and reference build.
- Cross-check of the wrapper total for `eval_triangles` against the cProfile cumulative time at 1e-2.
- Cross-check of the mode-kernel count: eigendecomposition cache misses plus SVD calls against the count implied by
  the m-value histogram of the kept pairs (written to the output as `mA_counts`) and the combination rule above.
- Cross-check of `tri_cm` and `tri_solve_rows` against the printed stats line of `hg_Vee_sparse`.

## Caps stated in advance

- One machine, one run per configuration, wall time with another heavy job possibly running. Shares are stated to
  the band widths above; differences of a few percent are not claimed.
- A share is not a speedup. LT6 gives the Amdahl cap only; no evaluator is built here.
- The triangle collection step is counted as non-triangle (bias low, stated above).
- Peak memory is read from the process but is not a prediction.
- The 1e-3 run takes about 3 h and 8.4 GB on one process. If it cannot be run, the proxy at 1e-2 is scored
  alone, LT1 and LT4 are reported as not scored, and the log-derived 20.0 percent is stated as an earlier-run
  value, not a measurement here.
- If the instrumented run does not reproduce the golden digits, the timings are void and reported as such.
- Up to 4 worker processes are allowed; the run uses 1.

## Disposition (appended 2026-10-26; predictions above unchanged)

**LT1 HELD, LT2 HELD, LT3 HELD (edge), LT4 HELD, LT5 MISSED, LT6 HELD. LT7: digits HELD; the overhead clause is NOT CONFIRMED, since the sum check leaves 1.1% and the falsifier as written (1 percent or more) fires; the cProfile control was not run.** (Header and the LT7 bullet corrected by the PM before commit.)
Script `debug/lih_triangle_share.py`; data `debug/data/lih_triangle_share_{0.01,0.001}.json` and the two `.log` files.
One process, one BLAS thread, both runs fresh. Wall times here (735 s at 1e-2, 10,187 s at 1e-3) differ from the
golden runtimes (2,008 s, 10,538 s), so shares are from one run each on this machine.

| Quantity | 1e-2 (nact 64) | 1e-3 (nact 196, target) |
|:--|:--|:--|
| Total wall | 735.4 s | 10,187.0 s |
| Triangle phase (`eval_triangles`) | 352.0 s, 47.9% | 2,159.2 s, 21.2% |
| h_Vne + g_Vne | 61.4 s, 8.3% | 3,070.4 s, 30.1% |
| h_T + g_T | 158.1 s, 21.5% | 2,275.5 s, 22.3% |
| g_Vee non-triangle (`vee_reduce_sparse`) | 16.1 s, 2.2% | 1,705.3 s, 16.7% |
| Fbar + sigma2 | 55.7 s, 7.6% | 860.3 s, 8.4% |
| Primitive load + reference build | 10.5 s | 3.7 s |
| Middle densities cm | 64 | 196 |
| Eigendecomposition builds / SVD calls | 192 / 204 (396) | 548 / 688 (1,236) |
| Factorization time (eigh + SVD) | 159.4 s (21.7% of wall) | 333.1 s (3.3% of wall) |
| Coulomb mode solves (all calls) | 189.2 s | 2,000.2 s |
| `tri_solve_rows` | 194,232 | 2,412,864 |
| E_R12 | -8.048140609752963 | -8.049045088639625 |

- **LT1 HELD.** 21.2% at 1e-3, inside 15 to 27 percent. The earlier log gave 20.0%.
- **LT2 HELD.** `tri_cm` is 196 and 64. Mode kernels 1,236 (band 1,000 to 1,400) and 396 (band 300 to 450). The m
  histogram at 1e-3 is {0: 108, +-1: 40 each, +-2: 4 each}, so 8 pairs with |mA| = 2, exactly the bound I used.
  The eigh count 548 and SVD count 688 match my 540 and 690.
- **LT3 HELD at the edge.** The factorizations are 333.1 s of 2,159.2 s, which is 15.4% (band 5 to 20). The
  `coul_mode_stack_any` timer reads 2,000.2 s, which is 92.6% of the triangle time, above the stated 75 to 90 band.
  That timer is not restricted to the triangle: the three sub-timers sum to 2,333 s, which exceeds the 2,159 s
  phase by 174 s, so some solve calls run in the non-triangle code. The triangle's solves plus its small
  matrix products are the remainder 2,159 - 333 = 1,826 s, 84.6% at most. The stated falsifier (solves below 75
  percent or factorization above 20 percent) did not fire. I score it HELD, with the 92.6% reading recorded as
  an over-count. At 1e-2 the sub-timers (348 s) match the phase (352 s).
- **LT4 HELD.** At 1e-3 the order is h_Vne/g_Vne 3,070 s, h_T/g_T 2,276 s, triangle 2,159 s, non-triangle g_Vee 1,705 s.
  The h_T/g_T block is larger than the triangle by only 5 percent.
- **LT5 MISSED.** At 1e-2 the share is 47.9%, outside 10 to 40 percent. The proxy overstates the target share by
  a factor 2.3. The triangle cost at 1e-2 is mostly fixed per-cm work (factorizations 21.7% of wall), while the Vne
  and T blocks shrank much more. The 1e-2 run is not a reliable proxy for the share.
- **LT6 HELD.** Factorization time is 3.3% of total wall at 1e-3 (cap 5%). Amdahl cap for a free triangle phase:
  1/(1 - 0.212) = 1.27x, inside 1.18x to 1.37x. For a free K_mu build: 1/(1 - 0.033) = 1.03x.
- **LT7 digits HELD; overhead clause not confirmed (see header); the cProfile cross-check was not run.** E_R12, c, h and E_trunc
  match the golden repr at both truncations (1e-2: -0.23890174597345007, 0.09609089427285555; 1e-3:
  -0.22689482141009634, 0.0899638869289805), the digits listed in the golden record. The 18-piece-level comparison
  against golden.json was not run; the four printed values were. The top-level timers sum to 653.9 s of 735.4 s
  (1e-2) and 10,074.6 s of 10,187.0 s (1e-3). The unaccounted 81 s and 112 s are untimed work (most likely the
  `GeminalM` construction and the loop glue); I did not isolate it. That is 1.1 percent at 1e-3, so the
  under-1-percent overhead claim is not confirmed as stated, though that gap is not necessarily timer cost.

**Reading.** In the corpus-best LiH run the triangle phase is about a fifth of the wall time, third of the four big
pieces. A free triangle phase would speed the run by about 1.27x. The K_mu kernel build, the part the operator
evaluator of record 2026-10-26_operator_evaluator.md would replace, is about 3% of the run (1.03x at best). The rest
of the triangle time is the prolate-Neumann Coulomb mode solves (about 2.4 million rows) and their products, which that
evaluator does not address. The larger targets are the h_Vne/g_Vne block (30%) and the h_T/g_T block (22%).
