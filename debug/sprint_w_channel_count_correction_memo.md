# Corrected SO(4) channel count W carried through Paper 28 (2026-11-07)

Driver: `debug/w_channel_count_correction.py` (sections A to F). Data:
`debug/data/w_channel_count_correction.json`, log alongside it.
Decision gate: **GO.** Every value was computed by two routes and the routes
agree. Tests are green, the Paper 28 loci are updated, and each period claim
has a verdict.

## 1. W and I: closed form and proof

**Corrected rule** (the PM's fix in `geovac/qed_vertex.so4_channel_count`).
Inside the window (|n1-n2| <= q <= n1+n2, n1+n2+q odd, q >= 1), both
spinors carry same-sign content ((n+1)/2, n/2). Each vector component
((q+-1)/2, (q-+1)/2) must pass the triangle and integrality conditions in
both SU(2) factors.

**Theorem.** W(n1,n2,q) = 2*[window]. Consequences:

- W is never 1.
- W is symmetric in (n1,n2).
- I(n1,n2) = sum_q W = 2 min(n1,n2), with I(0,n) = 0.

**Proof.** Take n1 <= n2 and d = n2 - n1. The window is
q in {d+1, d+3, ..., n1+n2-1}, which has n1 elements.

- Component V_A:
  - L triangle ((n1+1)/2, (q+1)/2, (n2+1)/2) holds for d-1 <= q <= n1+n2+1.
  - R triangle (n1/2, (q-1)/2, n2/2) holds for d+1 <= q <= n1+n2+1.
- Component V_B:
  - L triangle holds for d+1 <= q <= n1+n2+3.
  - R triangle holds for d-1 <= q <= n1+n2-1.
- All four integrality conditions reduce to n1+n2+q odd.

Each component is therefore open at every window q.

**Checks.**

- The geovac count matches an independent no-geovac implementation, and
  equals 2*window, for every triple with n1,n2 <= 40 and all q.
- The old count, re-implemented, reproduces 2min - 1 - [n1=n2] for
  n1,n2 <= 12. That confirms we are comparing against the right object.

**Symmetry (P28 sec. three-loop).** P28 said "W is NOT symmetric". That
was true of the old count and false of the corrected one. The text is fixed.

Old vs new W at the triples the tests pinned:

| (n1,n2,q) | old | new |
|---|---|---|
| (1,1,1) | 0 | 2 |
| (1,2,2) | 1 | 2 |
| (2,1,2) | 1 | 2 |
| (2,2,1) | 0 | 2 |
| (2,2,3) | 2 | 2 |
| (3,3,1) | 0 | 2 |

**Consequence.** The CG weight is a uniform factor 2 on the vertex-restricted
sums, so every "CG-weighted" sum equals 2^(#vertices/2) times its
vertex-restricted counterpart.

## 2. Values, old -> new

Abbreviations:
- **Route-geo**: geovac machinery with the corrected count.
- **Route-ind**: independent code, no geovac import (written in the driver).

| Quantity | Old | New | Routes, agreement |
|---|---|---|---|
| I(n1,n2) | 2min-1-[n1=n2] | 2min(n1,n2) | proof; geo = ind for n <= 40 |
| I(n,n) | 2n-2 | 2n | same |
| S_sunset | 2S_min - P^2 - Q = 3.89451620519784... | 2S_min = 4.95987387606844510882715900165876... | exact Fraction identity at N = 30 (geo); S_min closed form vs Hurwitz-tail sum at 70 dps (5e-70) |
| S^(3) (a=4, p=1) | 31.5725612075120... | **40.65538756358542254774433996911705224...** | closed form (200 dps) vs route 1 (old canonical + directly summed correction): 3e-120. Fully independent extrapolation (direct partial sums, ln^k N/N^j fit, no t-values): 6e-15, leading tail coefficient -32.0 (analytic) |
| S^(3) rigorous bracket, N = 4e6 | [31.57063, 31.57300] | [40.65339, 40.65578] | direct summation; L <= 4(ln m+1) checked (max ratio 0.91) |
| S^(3) closed-form ring | Q[pi^2, ln2, zeta(3), zeta(5), zeta(7), zeta(9), zeta(5,3)] | Q[pi^2, ln2, zeta(3), zeta(5), zeta(5,3)] | symbolic |
| S^(3)(N), N = 5/10/15/20/50 | 1.0616/3.5721/5.9083/7.8719/14.7984 | 2.753126/6.835894/10.181385/12.836576/21.580724 | geo = exact Fraction (N <= 20) and mp O(N^2) route, to 1e-20 |
| Increment fit, N in [20,120] | "270 N^-1.31" | 21.6 N^-1.23 (local slope 0.88 -> 1.42); tail at N = 120 ~ 12.0 (was ~7) | ind |
| S^(4) | bracket [~316.5, <0.5 wide] | [409.2860, 409.4247] (N = 4e6) | ind bracket; exact identities S4 = 8 M4 (n-space) = C4 - 6PC - 3S_o^2 + 27S_oP^2 - 27P^4 (o-space), geo, N = 8 |
| Sigma(n_ext), one loop (infinite) | (P28 table withdrawn earlier) | 2 sum_{k<=n} T(k) = 2nP - 2 sum_{m<n}(n-m)phi(m): 0, 1.92411, 3.23383, 4.22367, 5.01843, 5.68204 (n = 0..5), monotone | geo truncated = ind Fraction (n <= 5, N <= 8); closed form = tail route to 1e-25 |
| Two-loop SE raw total, n_max = 2..6 | 3->5 ratios > 2 | 1054.357, 2994.778, 6439.294, 11902.547, 19773.162 | geo = exact Fraction ind route |
| Truncated continuum VP, n_ch = 1 / 2 | 0 / 538361856/3603000625 | 73728/390625 / 1097979568128/2251875390625 | geo = ind exactly for n_ch = 0..5 |
| C_VP (n_max = 2, 3, 4, 5) | 0, 50471424/1779441125 = 0.0284, 0.0474, 0.0599 | 6912/109375 = 0.0632, 102935584512/1112150703125 = 0.0926, 0.1098, 0.1195 | ind; graph traces are W-free |
| C_SE (2..5) | 0, 0.275, 0.453, 0.582 | 0.503, 0.671, 0.793, 0.890 | ind |
| C_VP/C_SE (3, 4, 5) | 0.1030/0.1047/0.1029; "0.1035 +- 0.0009, CV 0.83%, ~3/29" | 0.1380/0.1384/0.1342; 0.137 +- 0.002, CV 1.7%; 0.1257 at n = 2 | ind |
| Fits over n = 3..5 | C_VP ~ n^1.48, C_SE ~ n^1.48 (text also said 0.0070 n^0.93) | C_VP ~ 0.054 n^0.50 (R^2 0.99), C_SE ~ 0.36 n^0.56 (R^2 1.00); C_F2 ~ n^0.57 unchanged | ind |
| C x F2 / (alpha/2pi) (n = 2, 3, 4, 5) | ---, 45.7, 64.9, --- ("grows") | 128.0, 149.3, 150.2, 143.6 (flat near 0.17, ~150x) | ind |
| Graph/continuum VP at n = 3 | ~35x | ~10.8x | ind |
| Per-mode rho range | 0 to 1.66 | 0 to 2.51 (graph 0, cont != 0 and the reverse both still occur) | saved graph shell sums + ind continuum |
| Census rule 5 (tab:selection_rules) | 48 of 48 (100%) have W = 0 | 8 of 48 (17%) have W = 0 | archived census driver with patched count |
| P34 vector-F2 overshoot, n = 3 | 12.7x | 41.5x (147.0x at n = 4); the 3.6x vector/scalar factor is unchanged | ind |

## 3. Three-loop structure (what changed algebraically)

**Old relation** (from `s3_decomp_setup`):

  S3_old = C - 2G - 16 P S_min - 8P^3 + 8PQ + R3

**Corrected relation:**

  S3 = C - 12 P S_min - 9 P^3

Both are bit-exact in Fraction arithmetic at N = 12, against geovac and the
independent count. The difference is

  Delta = 2G + 4 P S_min - P^3 - 8PQ - R3,

also bit-exact. Delta is a combination of three kinds of term:

- odd-weight doubles t(even, odd), all depth <= 1;
- lambda products;
- the depth <= 1 S_min.

All the triple-t and weight-10 content lives in C, and C is unchanged. So:

- all four trailing-1 identifications stand;
- the stuffle witnesses stand;
- the w10 reduction stands;
- the Charlton-Hoffman pair sums stand;
- the a-collapse stands.

The 11 t2 closed forms used by G were re-verified against direct Levin sums
at 60 dps, to about 1e-61. The closed-form G was re-verified against a
direct numerical G (Hurwitz inner tails, no t-value table) to 2e-58.

**New closed form** (25 terms; now printed as P28 eq:s3_closed_form):

    S3 = 128 pi^2 ln^2 2 - 32/3 pi^4 ln^2 2
       + 64 pi^2 ln2 - 304/3 pi^4 ln2 + 16 pi^6 ln2 - 2/3 pi^8 ln2
       - 96 pi^2 ln2 z3 + 16 pi^4 ln2 z3 - 80 pi^2 ln2 z5
       + 20 pi^2 z3^2 - 7 pi^4 z3^2 + 70 pi^2 z3 z5
       - 24 pi^2 z3 + 40 pi^4 z3 - 9 pi^6 z3 + pi^8 z3/2
       - 20 pi^2 z5 + 30 pi^4 z5 - 5/2 pi^6 z5
       + 6 pi^2 zeta(5,3)
       + 8/3 pi^4 + 374/45 pi^6 - 8231/3780 pi^8 + 124571/680400 pi^10 - pi^12/192

**Features of the closed form:**

- zeta(7), zeta(9), zeta(11) are absent entirely (Delta's zeta(7), zeta(9)
  and zeta(11) cancel against the old form's).
- Every rational constant is absent.
- The weight-10 block contribution
  -7pi^4 z3^2 + 70pi^2 z3 z5 + 6pi^2 zeta(5,3) + 179/170100 pi^10 is
  unchanged, so eq:s3_w10_contribution stands verbatim, with "Ident" now
  the corrected identified part.

## 4. Verdicts on the period claims

1. **Level-2 MZV ring.**
   - **HOLDS.** S_sunset = 2 S_min lies in Q[pi^2, ln2, zeta(3), zeta(5)].
   - **HOLDS, sharpened.** S^(3) lies in Q[pi^2, ln2, zeta(3), zeta(5),
     zeta(5,3)], a proper subring of the old statement.
2. **Realized depth <= k-1.**
   - **HOLDS at k = 2** (depth <= 1) **and at k = 3** (depth <= 2, from
     zeta(5,3) only).
   - **HOLDS at k = 4.** The corrected o-space relation keeps C4 at
     coefficient +1. The other cores (C, S_o^2, products) have depth <= 3.
     The stuffle-rank argument concerns C4's word set, which the
     correction does not touch.
   - The tower prose ("loop order bounds depth; realized depth grows more
     slowly") is unchanged.
3. **zeta(5,3) the only constant beyond single zetas:** **HOLDS.** Its
   coefficient is still 6 pi^2, because Delta carries no zeta(5,3). The
   cancellation of t(5,1), t(5,3), t(7,1), Li4(1/2) and the a-doubles is a
   property of C alone, so it is unchanged.
4. **Pure / mixed Tate:** **HOLDS.** The old sunset sentence ("correction
   terms P^2, Q are depth-1 elements of MT(Z), pure-Tate") is moot, since
   there are no correction terms. The sunset is 2 S_min, a mixed-Tate
   period over Z[1/2].

No PSLQ re-run was needed. The new closed form is algebraic (old closed form
plus an exactly decomposed correction), and three numerical routes confirm
it, one of them fully independent of the t-value machinery.

## 5. Paper 28 edits (W-dependent loci only)

- Abstract item (6): S^(3) value and bracket. Abstract "C x F2 diverges
  from alpha/2pi" changed to "stays about 150x alpha/(2pi)".
- sec:setup: W in {0,1,2} changed to W = 2 on the window, 0 elsewhere.
- sec:smin intro: the sunset is 2 S_min.
- subsec:smin_channel: (iii) I(n,n) = 2n; (iv) I = 2min, with the proof
  inline (replaces "proved in dossier round-1c").
- subsec:smin_sunset: S_sunset = 2 S_min and its value. The P and Q
  definitions (eq:P_def, eq:Q_def, cited nowhere else) were removed with the
  old relation.
- sec:three_loop:
  - W symmetric, so L = R.
  - Fit 21.6 N^-1.23, local slope 0.9 to 1.4, tail about 12, 21.58 vs 40.66.
  - Bracket and eq:s3_value.
  - The decomposition sentence now reads S3 = C - 12PS_min - 9P^3, with
    12 triples and 18 doubles.
  - The stage-1 error paragraph is labelled as concerning the earlier
    count's sum (value 31.5726).
  - Full-assembly paragraph: zeta(7), zeta(9), zeta(11) and rationals
    cancel; the identified part lies in pi^2, ln2, zeta(3), zeta(5).
  - eq:s3_ring updated; eq:s3_closed_form added; the verification sentence
    now cites route 1 at 3e-120; falsifier pointers updated.
- Prop. depth_k: k = 2 bullet names S_sunset = 2 S_min; k = 3 ring; k = 4
  now prints the o-space relation.
- tab:three_loop_convergence rows and caption.
- tab:selection_rules: rule 5 metric.
- sec:continuum_bridge: C and the ~11x figure.
- tab:c_times_f2 and the following paragraph.
- tab:three_diagram, the fits paragraph and the VP/SE paragraph. The
  "two-tier calibration structure" claim is replaced: with exponents 0.50,
  0.56 and 0.57 over three points there is no tier separation, and the
  types differ in magnitude only.
- rho range.
- Discussion paragraph on the projection constants.
- Conclusion paragraph on S^(3).
- Open question on the VP/SE ratio (3/29 removed).

The paper compiles cleanly over three passes. Gates: check_latex_escapes
group5 PASS; check_numeric_consistency group5 PASS. No registry key moved
(none of these quantities is registered).

## 6. Tests changed (each pins the corrected value; independent route noted in the test)

| Test file | Change |
|---|---|
| test_qed_vertex | Dipole test (W = 2*[n1 = n2 >= 1]); new W = 2*window test; W_total = 2min for n >= 0; ratio exactly 2 |
| test_graph_qed_continuum_bridge | W(1,1,1) = W(1,2,2) = W(2,1,2) = 2. VP(1) and VP(2) exact, plus an exact match to an inline 2min route for n = 0..4. n_max = 2: C = 6912/109375, CALIBRATION tier. n_max = 3: C = 102935584512/1112150703125. Decomposition has no W = 0 allowed triple |
| test_paper28_gap_closure_pins | New rows; added an inline exact-Fraction N = 5 check; slow N = 50 row |
| test_smin_decomposition | I = 2min, diagonal 2n via CG; sunset group rewritten (exact N = 30 identity with CG, value 2 S_min). The old 3.8945 pin was removed with the old relation |
| test_s3_decomposition | I = 2min; L = 2A; new o-space relation test. Canonical value now 40.6554...; the bracket must exclude 31.5726 as well as the stage-1 figures. New TestCorrectedClosedForm (closed form at 70 dps; independent route at 45 dps) |
| test_s4_stage1 | Chain collapse with the corrected kernel, = 8 M4, o-space relation. The 15-term identity is kept, explicitly for the pre-correction kernel. Bracket recomputed live at N = 1e6 ([409.03, 409.46]), excluding 316.57 and the old json bracket |
| test_two_loop_self_energy | Exact pins of raw totals at 3, 4, 5 (independent Fraction route) plus local slope > 2.5. The old "ratio > 2 per step" is false for the corrected sum (1.85 at 4 -> 5) |
| test_qed_three_loop | true_tail updated to 40.6554 - 12.8366 (the assertion passed before, but pinned the old value) |

All 17 consumer files plus test_fock_projection and test_fock_laplacian
pass: 995 passed, 38 skipped, 0 failed, with -n auto.

**`test_paper55_..._rhJ_closed_form[5]`: flaky, not W-related.** It
relies on the module-level `mp.dps = 60`, with REL_TOL 1e-50. At dps 50 the
s = 5 and s = 6 cases fail, and at dps 30 all three fail, so any earlier
test in the same xdist worker that lowers mp.dps breaks it. The fix is an
autouse dps fixture; that is out of scope and not done.

## 7. Out-of-P28 loci (NOT edited)

| Locus | Old | New |
|---|---|---|
| papers/group3_foundations/paper_55_periods_of_geovac.tex:997 | S^(3) = 31.5726... | 40.6554... |
| paper_55:1744-1745 | bracket [31.57063, 31.57300]; value 31.5726 | [40.65339, 40.65578]; 40.6554 |
| paper_55:1757-1758 | "polynomial in pi^2, ln2, odd zetas, zeta(11) cancelling" | pi^2, ln2, zeta(3), zeta(5); zeta(7), zeta(9), zeta(11) all cancel |
| paper_55:1767, 1769 | ring Q[pi^2, ln2, zeta(odd), zeta(5,3)]; "assembly gate 1.2e-198" | Q[pi^2, ln2, zeta(3), zeta(5), zeta(5,3)]; the gate is now route 1 at 3e-120 plus extrapolation at 6e-15 |
| papers/synthesis/group5_qed_gauge_synthesis.tex:306 | S^(3) = 31.5726... | 40.6554... |
| papers/group6_precision_observations/paper_34_projection_taxonomy.tex:4046-4049 | "C_VP x F2 diverges"; "12.7x overshoot at n_max = 3" | stays ~150x alpha/(2pi); 41.5x (vector F2) |
| papers/group5_qed_gauge/paper_36_bound_state_qed.tex:320 | "W_i in {0,1,2}" | W_i in {0,2} (W = 2 on the window) |
| paper_36 tab:ls8a_native (n = 3) | C = 18.97 | 18.96 (2994.7777/(4pi)^2 = 18.9648). The rest of the table matches my independent route, including n = 6 = 19773.2 |
| tests/test_paper55_witness_panels_computed.py:59, 190-224 | pins S3 = 31.5725612075 via an inline old kernel L = 2A - P - phi | 40.6553875636 with L = 2A (extrapolation reproduces it to 6e-15) |
| tests/test_s3_w10_identification.py | PIN['S3'] and IDENT_FORM are the old-count closed form (true identities of that sum, but labelled S^(3)) | new form is in test_s3_decomposition |
| debug/data/s4_bracket.json, spectral_projection_*.json | old-kernel artifacts | provenance only |
| docs/claim_test_matrix.md:126 | S^(3) = 31.5725612075 | 40.6553875636 |
| docs/validation_benchmarks.md:303, 327 | C = 50471424/1779441125; VP/SE ~ 3/29 | 102935584512/1112150703125; 0.137 +- 0.002 |
| docs/topic_to_paper_lookup.md:279, 291 | "C x F2 diverges"; "two-tier calibration" | flat ~150x; two-tier withdrawn |
| docs/curve_fit_audit_memo.md:73; docs/qa/group5.done.md:69-71, 118 | 0.1035 ~ 3/29 | withdrawn |
| docs/development_frontier_archive.md:68 | S^(3) = 31.5726 bracket | chronicle, dated; leave |
| CLAUDE.md sec. 3 row + docs/failed_approaches_ledger.md:51 | "C x F2 grows with n_max, doesn't converge to alpha/(2pi)" | the negative stands (no convergence to alpha/2pi), but the mechanism is now "stays ~150x", not "grows". PI/PM call; append-only rules apply |

Papers 2, 25, 30, 41 and 51 contain no printed W-dependent value; that was
checked by grep. Paper 33 rows R3/#3 are generic. Paper 18:2307 is generic.

## 8. Further defects / notes

- **Docstring-only edits** (same defect class, now false statements):
  - `geovac/graph_qed_continuum_bridge.py`: module header, the
    compute_projection_constant NOTE, mode_ratio_analysis and
    classify_projection_constant. These said the n_max = 2 continuum VP is
    zero because W(1,1,1) = 0.
  - `geovac/qed_vertex.py`: effective_pair_weight said "2min-1-delta";
    two_loop_sunset_weighted said "1 at boundaries, 2 in interior" and
    called S_min "NOT identifiable".
  - `geovac/qed_three_loop.py`: "W NOT symmetric".
  - No logic was changed.
- The archived census driver (`debug/archive/qed_arc/gn_selection_rule_census.py`)
  and `spectral_projection_mode_resolved.py` carry their own copies of the
  old count. The census copy adds an integrality check, which is why it
  reported 100% W = 0 (every unit empty, as the PM found). These are
  archived files, not edited; they were run with the corrected count by
  monkeypatching.
- P28's "C_VP ~ 0.0070 n^0.93 (R^2 = 0.95)" never matched its own caption
  (n^1.48). That is a pre-existing inconsistency, now superseded.
- The two-tier calibration reading and the 3/29 open question were the
  substantive claims that did not survive. The CxF2 falsification survives
  on magnitude rather than growth.
