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
