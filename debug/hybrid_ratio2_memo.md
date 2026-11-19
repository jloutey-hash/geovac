# Minimal order of the remaining panel hybrids at every exponent ratio (H1, H3, H5, H6, H8, H10)

Date 2026-11-19. Extension of `debug/hybrid_ratio_memo.md` (same mathematics, faster build). Unreviewed. No paper,
test, doc, geovac file or git state touched; the existing `debug/hybrid_ratio_*.py` files were imported, never edited.

## 0. Verdict

**GO for all six cases.** In every one of the 30 runs (6 cases x lt, mid, gt, eq1, eq2) the structure is exactly the
one the first memo proved for fam1, H2, H4, H7:

  order in R = 5 in each open region (lt, mid, gt), for every u > 0;
  order in R = 4 at a_d = a_c (eq1) and at a_d = a_c + a_p (eq2);
  no other exceptional ratio.

Every run: d_E = 2, k = 3 (open) or 2 (coincidence), zero-gcds nonzero rationals, label checks clean, cross-check
against the unpatched evaluator agrees to the working precision. No FINDING.

| case | panel id | lt | mid | gt | eq1 | eq2 | d_E | k (open/eq) | label issues | cross-check max rel. spread (40 dig.) |
|---|---|---|---|---|---|---|---|---|---|---|
| (3d0 3d0 \| 3d0 1s_B) | H1 | 5 | 5 | 5 | 4 | 4 | 2 | 3 / 2 | none | 2.2e-31 |
| (3d0 3d0 \| 1s 1s_B) | H3 | 5 | 5 | 5 | 4 | 4 | 2 | 3 / 2 | none | 5.6e-34 |
| (3d1 3d1 \| 3d1 2p1_B) | H5 | 5 | 5 | 5 | 4 | 4 | 2 | 3 / 2 | none | 1.6e-29 |
| (3d2 3d2 \| 3s 2s_B) | H6 | 5 | 5 | 5 | 4 | 4 | 2 | 3 / 2 | none | 5.1e-31 |
| (3p1 3d1 \| 3d0 2p0_B) | H8 | 5 | 5 | 5 | 4 | 4 | 2 | 3 / 2 | none | 3.4e-28 |
| (3d2 3d0 \| 3d-1 3d1_B) | H10 | 5 | 5 | 5 | 4 | 4 | 2 | 3 / 2 | none | 1.2e-24 |

With the four cases of the first memo, every hybrid of Paper 58's l <= 2 panel except the s-type H9 (covered
elsewhere, not run) now has "order = number of distinct rates + 1" at every exponent choice, both bounds,
conditional on the same items as before (Lemma GI, the Proposition of the first memo's Sec. 2, the evaluator's
algebra on symbolic rates).

## 1. What each run showed (same fields as the first memo)

* **E-part.** Two shifts c = -a_d and c = +a_d; at each shift exactly two E-labels (at eq1/eq2 one of the -a_d pair is
  ln R), with coefficients equal and opposite (`ratio_to_first` = -1 in all 60 shift groups). So d_E = 2 everywhere.
  The gcd of q's R-coefficients is a nonzero rational in every run (e.g. H1 mid 69120, H10 mid 345600, H3 eq1
  885735); q vanishes identically for no u > 0. q carries no constant (gamma, log).
* **L f.** Exactly the exponential rates of f survive: three in each open region, two at each coincidence (the -a_d
  rate merges with -a_c at eq1, with -(a_p + a_c) at eq2). Each surviving coefficient's numerator gcd over Q[u] is a
  nonzero rational (e.g. H10 mid 859963392000000, 20639121408000000, 1981355655168000000), and no constant survives
  into L f (no gamma, no log), so k is the full rate count for every u > 0 with no assumption on the logs.
* **Constants of f.** Logs of rate ratios in the open regions (e.g. gt: ln u, ln(u+2), ln(u+4), ln(u+6)); gamma
  (plus ln 3; ln 2, ln 3, ln 5 for H3) at the coincidences. They sit only on exponentials that L annihilates.
* **Labels** (`hybrid_ratio_summary.check`, imported unchanged): no two rates collide, no two E-labels at one shift
  collide, no Ei argument vanishes, and no u-factor of any coefficient denominator or zero-gcd has a root in u > 0, in
  all 30 runs. So the open regions contain no exceptional ratio.
* **Cross-check.** `hybrid_ratio2_check.py` (fresh process; never imports the engine, so geovac's
  `hybrid_closed_form` runs with its own normalised orbitals and numeric rational exponents) at u = 1 and u = 1/3 in
  each open region and at the single point of each coincidence, at R = 0.9, 1.7, 2.6: the ratio closed form /
  parametric form is one constant per point (the dropped normalisation), max relative spread over the three R as in
  the table. The looser spreads are precision, not mismatch: rerun at 80 digits, H10 gt fell from 1.2e-24 to 2.4e-65,
  H10 lt to 6.1e-69 (H5 gt in Sec. 4). 48 points in all (2 per open region, 1 per coincidence).

## 2. Optimisations (no change to the mathematics)

Profiling (`hybrid_ratio2_profile.py`, `hybrid_ratio2_profile_rep.py`; outputs `debug/data/hybrid_ratio2_profile_*.txt`)
on H6 found the time in three places, all sympy bookkeeping:

1. **One process per piece.** The closed form is a sum over multipole terms (L, L') and the two shell branches
   (x < R, x > R), and within a branch it is linear in the terms of the fixed-x inner double integral. Each piece =
   (term, branch, chunk of the inner double's Add-terms, round-robin) is built, integrated and collected in its own
   process (`hybrid_ratio2_job3.py`); `hybrid_ratio2_combine.py` adds the collected coefficient dicts exactly and runs
   `hybrid_ratio_analyse.analyse` (unchanged). Chunks are taken *before* the weight is multiplied in, so the
   expensive x -> R + w expansion is not repeated per chunk (a first version chunked after it and repeated it; that
   version is `hybrid_ratio2_job.py`). The per-term x-power assertion still runs on every product term.
2. **Integral accumulation.** `_integrate_with_e1` builds `total += ...` over its terms; re-flattening the growing sympy
   Add was 436 of 584 s on the profiled piece. The job calls it on ONE canonical term at a time (it is a term-by-term
   sum of the same calls) and collects each result immediately.
3. **Collector arithmetic.** About 90% of collection time was the pure-Python heuristic gcd of sympy's sparse field
   Q(u, R). `hybrid_ratio2_fastrep.py` rebinds the collector module's field object (in-process only) to the same field
   with numerator/denominator held as python-flint `fmpq_mpoly` (gcd in C). The walk, keys and constants are the
   collector's own code; finished coefficients are serialised as exact rational coefficient lists and rebuilt in the
   sympy field for the analysis.

Net effect on the profiled heavy piece (H6 lt, L = 4, x > R, chunk 0/6): 665 s -> 125 s.

**Validation (each change checked by exact reproduction before use):**
* v2 chunking on fam1 mid (the first memo's run): every coefficient of the collected form identical to
  `debug/data/hybrid_ratio_runs/2p0_2p0_1s_1s__mid.json` (same key set, 0 unequal coefficients), same d_E, k, order,
  rate list, constants, denominator factors. The `shifts`/`L` printouts differ only in which label of a pair is called
  "first" (dict order), i.e. an overall sign of q and of L; mathematically the same.
* v1 (chunk after expansion) vs v2 (chunk before) on H6 lt: identical forms, identical results.
* v3 (term-by-term integration + flint field) vs v2: all 10 pieces of fam1 mid and the heavy H6 lt (t2, gt, 0/6) piece
  identical key by key in Q(u, R) (`hybrid_ratio2_cmp_pieces.py`). Pieces from v2 and v3 were mixed in the final sums
  (v2 pieces already on disk were kept) only after this check.
* Every final form is also checked end to end against the unpatched evaluator (Sec. 1).

## 3. Timings

All 30 runs built, none near the 90-minute per-process cap. Machine: 16 cores, 13-14 worker processes.
Final driver (`hybrid_ratio2_driver3.py`): 1239 piece jobs, ~97 min wall for all six cases. Summed piece CPU
(includes the slower v2 pieces kept for H6 and part of H3):

| case | pieces | CPU (min) | longest piece (s) |
|---|---|---|---|
| H3 | 135 | 53 | 208 |
| H6 | 135 | 527 (mostly v2) | 1461 (v2) |
| H8 | 225 | 247 | 515 |
| H10 | 225 | 349 | 755 |
| H5 | 270 | 249 | 468 |
| H1 | 315 | 170 | 432 |

Combine + analysis: 2 to 190 s per run. Cross-checks: 20 to 50 s per evaluator build.

## 4. Where this could be wrong (in addition to the first memo's Sec. 6)

1. **Chunk linearity.** Correct only if every step after the inner double is term-wise linear. `_integrate_with_e1`
   expands and canonicalises term by term; `expand_fixed` and `.subs` are linear; the x-power assertion is per term.
   The exact v1/v2/v3 reproductions and the end-to-end cross-check would have caught a non-linear step.
2. **Mixed piece versions.** H6 and part of H3 sum v2 and v3 pieces; v2 = v3 exactly on the tested pieces, and the
   summed forms pass the end-to-end cross-check.
3. **Cross-check independence.** The reference is the same shell algebra at numeric exponents, not quadrature (as
   in the first memo); a shared algebra error would not show.
4. **Precision of the looser cross-checks.** H5 gt at 80 digits: 1.6e-29 -> 1.4e-69. The `hybrid_ratio2_checks/`
   files for H10 gt, H10 lt and H5 gt hold the 80-digit reruns; the table above quotes the 40-digit pass.

## 5. Files

New scripts `debug/hybrid_ratio2_*.py` (job = v2 piece, job3 = v3 piece used, fastrep, combine, driver, driver3,
check, summary, cmp_pieces, profile, profile_rep). Data under `debug/data/`: `hybrid_ratio2_runs/` (analysis + atom
form per run, plus the fam1-mid validation run), `hybrid_ratio2_results.json` (summary with label checks and
cross-checks), `hybrid_ratio2_checks/`, `hybrid_ratio2_pieces/`, `hybrid_ratio2_logs/`, `hybrid_ratio2_driver_C.out`,
`hybrid_ratio2_profile_*.txt`; superseded `hybrid_ratio2_pieces_v1/`, `hybrid_ratio2_logs_v1/`,
`hybrid_ratio2_runs_v1check/`.
