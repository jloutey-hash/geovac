# Avery load-vs-generate, step 1 (atoms): the gain is real but modest

**Date:** 2026-10-11 · **Type:** diagnostic (no geovac/ / papers / CHANGELOG / version edits)
**Driver:** `debug/avery_load_vs_generate_diag.py` · **Data:** `debug/data/avery_load_vs_generate_diag.json`, `debug/data/avery_lvg_run_n6.log`
**Supersedes for this question:** the August I/O ladder (`sprint_io_ladder_*_memo.md`), which counted tables and wrote the cost formula but never priced generation; its LG-3 lambda numbers predate the 2026-08-29 exact-rule fix and the 2026-09-25 He-hop removal and are not used here.

## Question
James Avery's intuition: the angular part is algebraic and can be produced on the device from labels, so only the small 1-D radial part has to be loaded. Does that beat loading every integral in a qubitized block encoding, and by how much?

## Key design point (corrects the August model)
Generating the coefficient on the device does not slot into standard coherent alias sampling: that method loads precomputed (alt, keep) tables, not coefficient values. The working version is **two-stage alias sampling over the multipole expansion**: stage 1 over (k, l's, m's) with weight |angular| * sum_n |R^k|, stage 2 conditioned on (k, l's) over the n-quartet with weight |R^k| / sum_n |R^k|. That is exact, needs **no on-device 3j arithmetic**, and loads two tables: angular (depends on l_max only) and radial (m-independent). The price is that LCU terms are per multipole, so lambda_F = sum_k |term_k| >= lambda_V.

## Result (single centre, hydrogenic, all n <= N, l < n, all m)
| N | M | L (load) | T1 ang | T2 rad | table ratio | Toffoli load | Toffoli fact | lambda_F/lambda_V | **net (lambda x Toffoli)** |
|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 2 | 5 | 27 | 12 | 10 | 1.2x | 27 | 23 | 1.01 | **1.16x** |
| 3 | 14 | 950 | 202 | 148 | 2.7x | 328 | 231 | 1.02 | **1.39x** |
| 4 | 30 | 14,425 | 1,694 | 1,128 | 5.1x | 1,412 | 787 | 1.04 | **1.72x** |
| 5 | 55 | 129,896 | 9,456 | 5,901 | 8.5x | 4,361 | 1,965 | 1.07 | **2.08x** |
| 6 | 91 | 809,790 | 38,155 | 23,832 | 13.1x | 11,407 | 4,139 | 1.09 | **2.52x** |

Toffoli = one QROAM pass per PREPARE (Berry et al. 2019 cost, min_k ceil(L/k)+b(k-1)), word b = 20 keep bits + alt-index bits for both schemes.

**Validation:** four textbook hydrogenic integrals (5/8, 17/81, 16/729, 93/512+(4/25)45/512) to 1e-16; 4-fold symmetry on 200 random quartets; bulk lambda and nonzero count equal a brute-force sum over all M^4 quartets at N=2,3 (lambda to 1e-13, counts exact).

## Reading
- **Direction: James's intuition holds for atoms.** The tables shrink a lot (13x at N=6, growing) and the net advantage grows with basis size.
- **Magnitude: modest.** QROAM costs ~sqrt(table), so a 13x smaller table buys ~2.8x fewer Toffolis, and the per-multipole LCU costs a slowly growing lambda penalty (1% -> 9%). Net ~1.2x -> 2.5x over N = 2..6, rising slowly. Not the orders-of-magnitude regime of first-quantized plane waves.
- **Implementation note worth telling him:** no 3j calculator on the chip is needed. "Generate" is best realized as "load a tiny angular table + a radial table and combine in two stages."

## Caveats (do not overclaim)
1. **Baseline is plain sparse loading.** Double factorization / THC also compress the table; the atomic multipole split is itself a symmetry-exact low-rank factorization, so a DF baseline may capture much of the same gain. This is the fair-baseline question for James.
2. Only the PREPARE data-loading step is costed; SELECT, reflections and precision overheads are common to both and dilute the ratio.
3. Counts use ordered/4 symmetry (slightly undercounts LOAD; conservative toward LOAD). Orbital-index words not charged to LOAD (conservative toward LOAD).
4. **Atoms only.** At two centres l is not conserved (Papers 14/58); the angular factorization around one centre fails for cross-centre terms. Molecules are untested here and are the hard case.
5. Basis = hydrogenic (orthonormal). Shared-k Sturmians after Loewdin keep the same angular factorization (the Sturmian overlap is diagonal in (l,m) and m-independent), so table sizes carry over; their lambda differs (August memo: Loewdin inflates it).

## Next (step 2, only after James weighs in on the baseline)
Add a double-factorization baseline on the same atomic tensors; extend to a diatomic to measure how much of the gain survives two centres.
