# Sprint: `/qa group1` re-run-2 of Batch 1 (rr1b) + dedicated citation sweep — 2026-06-21 (v4.31.0)

Second certifying re-run of Batch 1 (Papers 45–49 + synthesis), then — per PI
direction on the result — a **dedicated exhaustive citation sweep** of 45–49 to
drain the C4 tail in one pass (rather than keep iterating re-runs).

## 1. rr1b verdict: FAIL (calibrated) → remediated

14-agent panel, worktree `../geovac-qa-seed-group1-rr1b` (removed; no seed leaked).
Seed key `debug/qa/group1_rr1b_seed_key.json`.
- **Sensitivity 4/4:** S-claims-C14 (p46 theorem-title "(established)")→claims-46;
  S-citation-C4 (p47 Bunce–Deddens vol/year)→citation-47; S-code-C2 (Datta chain
  `==0`→`>=0`)→code-49; S-synthesis-C9 (Paper-47 composite "full metric level")→synthesis.
- **Specificity 5/5:** controls affirmed; **the v4.30.0 keystone fixes HELD**
  (claims-49 confirmed strictness "generically"; code-47 confirmed the analytic
  proof-by-argument framing) — no keystone defect this round.

**Finding profile shifted:** claims/code converging (only SMALL residuals), but
the C4 citation dimension surfaced 2 new LARGE issues no prior pass caught →
motivated the dedicated sweep.

## 2. Genuine MATERIAL defects remediated

**Claims (SMALL, rr1b):**
- p47 Remark 7.6: "Theorem 7.3 holds on the natural substrate" → "construction
  restricted to the natural substrate; even there the metric-level closure is
  descoped (inherits Paper 45's degenerate seminorm; only norm-resolvent/spectral
  established)".
- p47 L2 "cb-norm" contradiction: "$\cbnorm{S}=2/(\nmax+1)$" → "the constant
  $2/(\nmax+1)$ is the Plancherel-mass maximum (the map is UCP, cb-norm 1)",
  matching the paper's own Status-note correction (iii).
- p49 ×4 summary spots (abstract / honest-scope / aggregate-status / conclusion):
  the v4.30.0 panel-body relabel (Datta→Umegaki) hadn't propagated; the 66–81 nats
  are now labeled "illustrative Umegaki relative-entropy", with the load-bearing
  $D_{\max}$ chain pointed to the 96/96 test.

**Citations LARGE (rr1b):**
- p47 `hekkelman_mcdonald2024` misattribution: the abstract + §"Position vs
  Hekkelman–McDonald" credited it with "Tauberian estimates on flat tori /
  extension to ℝ^d via covering map / compact-to-non-compact machinery" — NOT in
  arXiv:2412.00628 (which is a noncommutative integral on spectrally truncated
  triples + quantum ergodicity). Corrected the characterization; the outer-arrow
  machinery is the (correctly cited) Reed–Simon §XIII.16.
- p48 `muller2022` fabrication: "Lorentzian GH convergence and Cauchy slabs, CMP
  391 (2022) 855–882" — no such paper. → real Müller "Lorentzian GH theory and
  finiteness results," Gen. Relativ. Gravit. 54 (2022) Art. 117, arXiv:1912.00988.

## 3. Dedicated citation sweep (PI direction) — C4 tail drained

Exhaustive per-paper verify (one citation-reviewer per paper, opus, every bibitem
web-checked). Produced the complete remaining tail; applied ~25 corrections across
45–49 (many duplicated across p48/p49's shared Lorentzian-cluster bibliography):
- **Wrong titles (correct ID):** minguzzi_suhr ("Lorentzian metric spaces and
  their GH convergence"), leimbach ("for tori"), latremoliere2018 (+"quantum"),
  sakovich_sormani, mondino_ryborz_samann (+author V.→ Ryborz), ketterer,
  nieuviarts2025_v2 ("Emergence of Time…"), sormani_vega ("Null distance on a
  spacetime"), franco_eckstein ("An algebraic formulation of causality…"),
  latremoliere2026_spectral_c1, martinetti2026_adjacent.
- **Wrong metadata:** kostant1965 year 1965→1969; kubota author T.→H.;
  connes_rovelli end-page 2917→2918; toyota M.→R.; franco label 2014→2013.
- **Fabricated/unverifiable → real work:** farsi_latremoliere2024 (non-existent
  un-IDed preprint → "Collapse in NCG and spectral continuity," arXiv:2404.00240);
  bertozzini (filled in SIGMA 6 (2010) 067, arXiv:1007.4094).
- **Duplicate bibitems removed:** hekkelman_mcdonald2024b in p47 (uncited→deleted)
  and p48 (cited→repointed to hekkelman_mcdonald2024, then deleted).

## 4. Reconcile catches
- **p48 PI-grade caveat resolved:** the p48 sweep flagged that Mondino–Sämann
  Def 3.6/3.8/Thm 6.2 (underpinning Bridge B3/B4) might not exist (the old
  "Latrémolière Thm 5.5" class). The p49 sweep **independently verified** MS Def
  2.3/3.8/4.4 exist in arXiv:2504.10380 and say what's attributed → NOT a LARGE;
  the MS definitions are grounded.
- p48 `latremoliere_metric_propinquity_2015` "The **dual** GH propinquity" is the
  real 2015 JMPA paper (correct) — not the "dual" error (that was only in
  `_metric_st_2017`, fixed v4.30.0).
- p49 bizi `1611.07062` correct (the rr1 seed was worktree-only).

## 5. Honest scope / deferred
- All 5 papers compile errors=0/undef=0 (p45 27pp, p46 26pp, p47 19pp, p48 30pp,
  p49 36pp); the two bibitem deletions/repoints verified non-breaking. C5/C11/C13/
  C14/C15 group1 PASS. Worktree removed; no seed leaked.
- **Residual NIT (cosmetic, optional):** stale cite-KEY strings whose displayed
  bibitems are now correct (e.g. `latremoliere2018`/`_metric_st_2017` key-years,
  `bousso_casini_fisher_maldacena2020` key author names, `_proceedings` suffix) —
  invisible to readers; left as-is.
- **Certification status:** with the C4 tail now drained (exhaustive verify, not
  sampling) and the claims/code dimensions converged, the **next clean re-run
  should certify Batch 1**. Batches 2/3 re-runs still pending — they share much of
  this Lorentzian-cluster bibliography, so the sweep likely reduced their tail too,
  but each still needs its own clean pass.

## 6. Next
Clean re-run of Batch 1 (PI-invoked) to certify; then Batches 2/3 re-runs.
