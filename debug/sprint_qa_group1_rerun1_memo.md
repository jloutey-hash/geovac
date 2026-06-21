# Sprint: `/qa group1` CERTIFYING re-run — Batch 1 (Papers 45–49 + synthesis) — 2026-06-21 (v4.30.0)

PI-invoked `/qa group1` **certifying re-run** (convert the rc1/rc2/rc3
FAIL→remediated cycles into a trustworthy PASS). Lorentzian-first; Batch 1 =
45–49 + synthesis. PI directive on the FAIL: **remediate all now (my judgment)**.

## 1. Verdict: FAIL (calibrated) → remediated

14-agent panel (claims ×5, citation ×5, code ×3 [45/47/49], synthesis ×1),
path-pinned to worktree `../geovac-qa-seed-group1-rr1` (removed; no seed leaked;
real corpus clean on all 4 anchors). Seed key `debug/qa/group1_rr1_seed_key.json`.

- **Sensitivity 4/4:** S-claims-C14 (p48 T6 "(established)") → claims-48;
  S-citation-C4 (p49 bizi 1611.07062→1611.09062) → citation-49; S-code-C2 (p45
  K⁺ annihilation `<TOL`→`<1e12`) → code-45; S-synthesis-C9 (Paper 45 "first
  Lorentzian propinquity convergence instance") → synthesis.
- **Specificity 5/5:** M1–M5 controls all affirmed.
- **The re-run worked as designed:** calibrated, and it caught residuals rc1
  missed **plus two keystone-correctness issues all prior passes missed**.

## 2. Genuine MATERIAL defects remediated (PM-verified vs primary text/code)

**Keystone-level (PI flagged "remediate all, my judgment"; both independently verified):**
- **p49 — strict super-additivity was FALSE as stated.** Claimed "strict
  inequality whenever the three KMS states are pairwise distinct." I reproduced
  code-49's counterexample: distinct *commuting* states give D_max chain deficit
  = **exactly 0** (saturation). Restated to "strict **generically** (substantially
  positive on the off-orbit panel and all 96/96 tested D_max cells; can saturate
  for non-generic configurations whose D_max ratios align; pairwise distinctness
  necessary not sufficient)" at all ~6 sites (theorem `thm:strict_super_additivity`
  l.1732, D_max ingredient l.1772, restatement (ii) l.1967, l.1794, abstract
  l.241, intro l.472).
- **p49 — the numerical panel was mislabeled.** `tab:uhlmann_panel` was captioned
  "Datta max-divergence chain inequality deficits" but the driver
  (`q1prime_phase2b3_panel_compute.py`) computes `_relative_entropy` = **Umegaki**
  (I confirmed the source). The numbers (66.998/68.720/81.256) are *correct
  Umegaki values* → **relabeled** to "Umegaki relative-entropy chain deficits
  (illustrative)", **demoted** to an illustrative cross-check (the Umegaki chain
  is not generic — the paper's own finding), and **pointed** the load-bearing
  D_max positivity to the genuine 96/96 test `tests/test_wh7_b3_phase3_sprint2.py`
  (code-49 confirmed SOUND). Not recomputed — the numbers were right for what
  they are.
- **p47 — three surviving keystones (outer-arrow norm-resolvent, Main theorem,
  three-carrier) had no backing test.** They are analytic results (the T→∞
  de-compactification is not computationally implemented). Added a Status Remark
  marking `thm:main` as analytic **proof-by-argument** (Reed–Simon §XIII.16),
  with R^outer established and R^inner descoped — closing the provenance gap
  honestly without claiming a nonexistent test.

**C14 zombies (rc1-incomplete):**
- p46 l.286 "Paper 45 closed the convergence-theorem leg" → "established the
  K⁺-restricted weak-form **degeneracy** theorem (convergence claim retracted)".
- p48 §6.3 "the synthetic-side closure exists; the bridge transports it back" →
  conditional/descoped ("would exist on the repaired seminorm … rests on Paper
  45's degenerate Λ_prop, so T6 is a proposed/conditional route"); + §6.3 title
  "(proposed/descoped)" + l.2121 "deepest **proposed** content … descoped".
- p47 `thm:main` descope Remark (above) + §6 section title "(descoped)".

**C7:** p46 l.1289 "Paper 38's SU(2) Riemannian propinquity bound" → "state-space GH bound".

**C4 (load-bearing):** `latremoliere_metric_st_2017` spurious "dual" + wrong
vol/year (Adv. Math. 411 (2023) 108790 → **404 (2022) Paper 108393**) fixed in
**p48 + p49** (cited in body). `farsi_latremoliere2024` WRONG-ID (JFA 286/110293
is Matringe–Offen–Yang, not Farsi–Latrémolière) **removed** from p47 (bibitem +
the 4-cite bundle; the claim is carried by the 3 real companion cites).

## 3. Reconcile catches (over-flags / not-the-error)
- p48 `latremoliere_metric_propinquity_2015` ("The **dual** GH propinquity," JMPA
  103 (2015) 303–351) is **correct** — "The Dual GH Propinquity" is a real 2015
  Latrémolière paper; the grep "dual" hit there is not an error. The genuine
  "dual" error was only in `latremoliere_metric_st_2017`.
- p49 bizi `1611.09062` was the **seed** (worktree-only; real corpus has the
  correct `1611.07062`) — no real-corpus fix.

## 4. Honest scope / deferred
- **All MATERIAL + C7 + load-bearing-citation defects fixed.** p46/47/48/49
  compile errors=0/undef=0 (26/19/30/36 pp). C5/C11/C13/C14/C15 group1 PASS.
- **Deferred (NIT cosmetic-citation sweep, rc1/rc2/rc3 precedent):** ~10
  title-paraphrase / key-year slips on **non-load-bearing adjacent-framework**
  refs — p48 muller2022 (wrong title/journal/vol), sakovich/mondino_ryborz/
  minguzzi_suhr/ketterer/nieuviarts (paraphrased titles, correct IDs); p47
  `hekkelman_mcdonald2024b` phantom duplicate (already self-flagged "PI to
  consolidate"); p49 bousso cite-key author mismatch; p46 minguzzi_suhr/leimbach
  titles + latremoliere2018 title-missing-"quantum"; franco 2014→2013 labels.
- **Certification status:** Batch 1 re-run is **FAIL→remediated**, NOT yet a
  clean PASS. A *clean* re-run (zero material defects) is what certifies — the
  re-run finding deep issues is the gate working; certification requires
  convergence (a subsequent clean pass). Batches 2/3 re-run still pending.

## 5. Next
Re-run Batch 1 again (should now be clean) OR proceed to Batch 2/3 re-runs, per
PI. Certification = a clean calibrated pass on all three batches.
