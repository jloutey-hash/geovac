# Sprint: `/qa group1` Bite B sub-bite 2 (Papers 45–49) — 2026-06-18 (v4.24.0)

PI-invoked `/qa group1` Bite B, sub-bite 2 = the **descoped/partial Lorentzian
core**: Papers 45 (degeneracy theorem), 46 (strong-form, descoped), 47 (norm-
resolvent, partial), 48 (Krein–MS bridge, partial), 49 (OSLPLS, partial) +
synthesis.

## 1. Verdict: FAIL (calibrated) → remediated (full sweep, PI direction)

13-agent panel (claims ×5, citation ×5, code ×1 [+1 deferred], synthesis ×1),
path-pinned to worktree. **Calibrated: sensitivity 5/5** (S-claims-C8 p45,
S-claims-C14 p46, S-citation-C4 p48, S-code-C2 p45-test, S-synthesis-C9 all
caught), **specificity 5/5** (every honest degeneracy/descope control affirmed
SOUND). Heaviest FAIL of the campaign. Seed key `debug/qa/group1_B2_seed_key.json`.

## 2. The defining pattern: "Status-note-only descope"

Every one of Papers 45–49 has an **honest, correct Status note** recording the
descope — but the descope was **never propagated into the theorem bodies,
abstracts, appendices, or header comments**, which still re-asserted the
descoped/retracted Lorentzian-metric convergence as established. Each paper
contradicted itself. Remediation = propagate each paper's own Status-note
framing into the zombie re-assertions (no judgment calls; the correct wording
already existed in-paper).

**Fixed:**
- **p46** — Appendix B "closes G1$'$" → "proposes a route (descoped; reproduces
  the degenerate seminorm)"; §6.6/6.7/8.5 "convergence closed / structurally
  complete" → descoped; §7 the nonexistent **Latrémolière Thm 5.5 / Def 3.5**
  (the fabrication that drove the descope) re-pointed to van Suijlekom
  state-space GH; §1.2 "closes the strong-form leg" → "addresses (descopes)".
- **p47** — header comment + Honest-scope block (the metric-level Thm 7.3
  "closure") → descoped; §7 "This closes Q1" → "would close … descoped".
- **p48** — §2.2 "Paper 45 establishes the convergence → 0 / Paper 46 free
  upgrade" → degeneracy/rate-formula; T3 "first quantitative pLGH-convergence
  panel" → "rate panel (descoped)"; T6 "G2 closure established at theorem-grade
  rigor" → "proposed route, conditional/descoped"; §6.5 "genuinely new content"
  → "conditional, proposed".
- **p49** — B4$'$ theorem tagged DESCOPED (hypothesis fails — degenerate
  seminorm); aggregate Bridge Theorem "all four properties" → "three (B1$'$–B3$'$)
  established + B4$'$ descoped"; "free upgrade" / "verified bit-exact" /
  "all load-bearing checks pass" → rate-formula/state-level + B4$'$ descoped;
  "closes Q1$'$" → "at the OSLPLS/state level (metric leg descoped)".
- **p45** — Paper-46 cross-ref bibitem "closing the strong-form gap" → descoped;
  `test_paper45_asymptotic_rate` docstring relabeled (it verifies the surviving
  SU(2) spatial rate, not a Lorentzian Λ^L convergence).
- **synthesis** — no edit (the real narrative was already correct; the only
  synthesis hit was the worktree seed).

## 3. Citations (verified vs arXiv + primary PDFs)

- `devastato_lizzi_martinetti2018` wrong-ID in 45/46 (fabricated titles "Time as
  an emergent property…") → **Lorentz signature and twisted spectral triples**,
  Devastato-Farnsworth-Lizzi-Martinetti, JHEP 03 (2018) 089, arXiv:1710.04965.
- `hekkelman_mcdonald2024` fabricated/wrong-title in 46 (T^d/2403.18619=OpenMP,
  removed), 47/48 (wrong title on real 2412.00628 → corrected).
- `bizi_brouder_besnard2018` p47 **confabulated title+venue** ("Spectral action
  in Lorentzian signature, CQG 35 175004") → real "Space and time dimensions…,
  JMP 59 (2018) 062303" (arXiv was correct).
- `latremoliere_metric_st_2017` p47 wrong "dual"/vol/year → Adv. Math. 404
  (2022) 108393.
- `latremoliere2018` (if present) → "quantum" GH propinquity, TAMS 368 (2016).
- **Theorem-number verification (PI-requested):** van den Dungen Prop 4.1,
  Nieuviarts Def 2.2, Mondino–Sämann Def 2.3/3.8/4.4, the Latrémolière
  hypertopology Def/Thm numbers — all **GROUNDED** (no fix).

All 5 papers compile errors=0 / undefined=0; C11/C13/C14 PASS; p45 tests green.

## 4. Honest scope / deferred
- **Material defects all fixed** (the systematic C14 zombies + the LARGE/SMALL
  citations). The surviving math (degeneracy theorems, norm-resolvent §4,
  OSLPLS/cocycle state-level results, B1$'$–B3$'$) is genuinely backed.
- **Deferred NIT-grade title slips** (real arXiv IDs, verbatim-title drift only,
  citation-reviewers rated NIT): p49 ×5 (martinetti/nieuviarts/ketterer/
  latremoliere-C1/hekkelman titles), p48 minor (kostant year, minguzzi/
  mondino-ryborz/nieuviarts/ketterer titles), p47 hekkelman prose attribution +
  reed_simon §-pointer; p45 test docstring 2 residual label lines. Logged for a
  cosmetic citation-title sweep.

## 5. Next
Bite B sub-bite 3 (Papers 39, 52, 53) — the last group1 sub-bite.
