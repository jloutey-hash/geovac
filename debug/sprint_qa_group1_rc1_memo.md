# Sprint: `/qa group1` re-certification Batch 1 (Papers 45–49 + synthesis) — 2026-06-19 (v4.27.0)

PI-invoked `/qa group1` re-certification re-runs, batched ~5 papers (token
budget). PI reordered to **Lorentzian-first**:\ Batch 1 = the C14-heavy descoped
core (45, 46, 47, 48, 49) + synthesis. PI directive: **remediate now, defer
re-run**; and **build the deterministic inline-arXiv-ID cross-check**.

## 1. Verdict: FAIL (calibrated) → remediated

12-agent panel (claims ×5, citation ×5, code ×1 [p45], synthesis ×1),
path-pinned to worktree `../geovac-qa-seed-group1-rc1` (removed; no seed leaked;
real corpus confirmed clean). Seed key `debug/qa/group1_rc1_seed_key.json`.

- **Sensitivity 4/5:** S-claims-C8 (p45 annihilation→"distance established")
  caught by claims-45; S-claims-C14 (p47 G2-metric→"established") by claims-47;
  S-synthesis-C9 (p48 inheritance→"established") by synthesis; S-code-C2 (p45
  test `allclose(U,U)` tautology) by code-45. **S-citation-C4 MISSED** — an
  inline `arXiv:2504.10830` (transposition of the bibitem's `2504.10380`)
  slipped past the bibitem-focused citation reviewers.
- **Specificity 5/5:** all controls (annihilation thm, norm-resolvent, §Q1
  closure, devastato cite, product-carrier signature-agnostic) affirmed sound.

**Headline: the B2 full-sweep (v4.24.0) was incomplete.** The re-cert surfaced
residual "Status-note-only descope" zombies + citation defects B2 missed. All
findings PM-verified against the real corpus (§9 reconcile) before fixing.

## 2. The citation calibration gap → new deterministic gate C15

The citation reviewers verify *bibitems* thoroughly (they caught bizi, Datta,
leimbach, hekkelman) but missed an inline-prose transposed digit. Per PI
direction, built **`debug/qa/check_inline_arxiv.py` (C15)**:\ flags any inline
arXiv ID that is a near-match (same length, Hamming ≤ 2) of a bibitem ID without
equalling it — the transposed/typo'd-ID class, deterministically + cheaply.
Complements (does not replace) the C4 LLM reviewer (which still owns wrong-title
/ wrong-venue — right ID, wrong metadata). Test `tests/test_inline_arxiv_check.py`
(4/4); registered as shared C15 in `docs/qa/criteria.md`; group1 PASS.

## 3. Verified genuine defects remediated (all B2-incomplete)

**C14 residual zombies** (descoped re-asserted as established):
- p46 Appendix B `thm:enlarged_main` "β-L5 closure" + opening prose ("records
  the closure", "answers both questions affirmatively") → proposed route,
  proof-sketch grade, descoped (consistent with §B.5's own "not a closure"); +
  inline descope tag in the theorem body.
- p47 `thm:g2_metric` title "G2-metric closure" → "G2-metric route (descoped)".
- p48 T6 "closes G2 metric-level" cluster ×4 (contribution list, honest-scope
  §1.4, theorem title+body, §9.2) → "proposed/conditional route (descoped)".
- p49 "confirms B4′ convergence transport at machine precision" → "reproduces
  the closed-form rate (consistency check, not a metric convergence; descoped)".
- synthesis l.1593 "all four properties hold theorem-grade" → "(B1′)–(B3′)
  theorem-grade; (B4′) Λ-inheritance descoped"; abstract Paper-46 4/π
  co-citation dropped (Paper 46 is descoped, contributes no live Lorentzian 4/π).

**C7** (Paper 38 = state-space GH, not Latrémolière propinquity): p48 ×4
"Riemannian propinquity hypertopology" → "state-space Gromov–Hausdorff".

**C4 citations** (B2 fixed these in p47 only / new):
- bizi_brouder_besnard2018 confabulated title "Spectral action in Lorentzian
  signature, CQG 35 175004" → "Space and time dimensions…, JMP 59 (2018) 062303"
  in **p48 + p49** (arXiv:1611.07062 was always correct).
- leimbach_vs2024 wrong-author/venue → sole-author M. Leimbach, to appear JNCG
  (not van Suijlekom co-author, not Adv. Math.) in p48 + p49.
- p49 hekkelman wrong-title → "A noncommutative integral on spectrally truncated
  spectral triples…" (arXiv:2412.00628).
- **p49 Datta "Theorem 11" misattribution** — the max-divergence chain
  inequality `D_max(1‖3) ≤ D_max(1‖2)+D_max(2‖3)` is NOT Theorem 11 of Datta
  2009. PM verification:\ it is **elementary** (immediate from
  `D_max(ρ‖σ)=log min{λ:ρ≤λσ}` by transitivity of operator order) and
  framework-verified (96/96, `test_wh7_b3_phase3_sprint2.py`). So the math + the
  p49 headline (strict super-additivity) STAND; only the false theorem-number
  attribution was wrong. Reframed to "elementary from the D_max definition;
  framework-verified; Datta for the max-divergence" — a citation fix, not a
  retraction (the reviewer over-stated severity; reconcile-rule catch).

## 4. Honest scope / deferred
- **Material findings all fixed.** All 5 papers compile errors=0/undefined=0;
  C5/C11/C13/C14/C15 group1 PASS.
- **Deferred (NIT, B2 precedent):** ~12 cosmetic citation title-paraphrases with
  correct arXiv IDs (ketterer/sakovich/nieuviarts/sormani-vega/mondino-ryborz
  paraphrased titles, franco 2014→2013 year, kostant 1965→1969 year, bousso body
  co-author list, latremoliere-dual title) + p49 bousso cite-key. Logged for a
  cosmetic title sweep.
- **Re-run deferred** (PI direction):\ Batch 1 re-run after Batches 2–3.

## 5. Next
Batch 2 (Papers 42, 43, 44, 53 — Lorentzian foundations + disk), then Batch 3
(29, 39, 40, 50, 52 — results backbone), then re-run all to certify.
