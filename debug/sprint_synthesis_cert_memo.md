# Sprint memo — `/qa synthesis` certification + §2 compaction (2026-07-04)

**Verdict:** the **synthesis layer is CERTIFIED** (7th and FINAL `/qa` target) → the whole-corpus branch-QA sweep is **COMPLETE** (trunk + group1–6 + synthesis). Releases `v4.67.1`–`v4.67.4`. Companion to the group6 cert (`v4.67.0`, `debug/sprint_group4_8thcert_memo.md`-style; group6 notes in `docs/qa/group6.done.md`).

Canonical DoD + per-run chronicle: `docs/qa/synthesis.done.md`. Seed keys (gitignored): `debug/qa/synthesis_{,delta_,redelta_}seed_key.json`.

---

## 1. Scope

`/qa synthesis` = the top-of-tree readiness gate, fired after all six paper groups certified. Target:
- **PRIMARY:** the field guide `papers/synthesis/geovac_field_guide.tex` (~931 lines, the cross-group "what is GeoVac" narrative, 27 cites) — the one *uncertified* synthesis document, full review.
- **SECONDARY:** the six group syntheses (each already the C9 dimension of its own branch cert) — checked for cross-group consistency + faithfulness to the now-frozen certified papers.

**Code/test-backing dimension WAIVED** (PI direction at freeze): the synthesis layer inherits its backing from the already-certified group-synthesis papers; C1/C2 is out of scope, not an unexercised gating dimension. Gating dims = claims (C3/C5/C6/C8) + citations (C4) + synthesis-faithfulness (C9) + deterministic (C10–C17); verdict = AND across those four.

## 2. The arc (four moves)

**(A) 1st cert FULL run = INCONCLUSIVE — PM setup error.** Panel calibrated (6/6 seeds, specificity clean across 3 agents: claims-FG, cross-group, citation). BUT I built the seed worktree from HEAD while `geovac_field_guide.tex` carried **uncommitted working-tree edits** (pre-session PI fixes that had already corrected He 0.019%→0.004%, 28→37 molecules, "two"→"six" syntheses). So the claims reviewers reviewed **stale text** — their "genuine" findings were already fixed. The claims dimension was therefore not validly exercised on the real corpus ⇒ INCONCLUSIVE. The **citation dimension was valid** (bibliography identical HEAD↔working-tree) and found genuine MATERIAL.

**(B) Field-guide re-run on the correct text = CLEAN.** Rebuilt the worktree synced to the working tree, re-ran a combined claims/faithfulness reviewer with fresh seeds (He-0.4% drift, Lorentzian re-zombie) — both caught, controls verified faithful, zero genuine defects. Confirmed the corrected field guide is clean.

**(C) Citation delta (PI-directed, "clean = pass") = DEFECTS → remediated (v4.67.1/.2).** Calibrated 2/2 (re-corrupted the two v4.67.1-fixed cites). Scanning the full group1 bibliography, it surfaced **genuine wrong-arXiv/venue citations in the already-certified group1 synthesis** — a class its branch cert missed:
- `hekkelman2022` arXiv **2206.13744 = a Kerr–Melvin black-hole imaging paper** → real "Truncated geometry on the circle," LMP 112 (2022) 20, **arXiv:2111.13865**.
- `hekkelman_mcdonald2024` arXiv **2403.18619 = an OpenMP shortest-path CS paper** → wrong ID removed, flagged pending author confirmation (real ID unconfirmable via search).
- `latremoliere_metric_st_2017` "Adv. Math. **415** (2023), 108876, 88pp" → **404 (2022), 108393, 56pp** (arXiv:1811.10843 journal-ref).
- Plus title/label SMALLs: `hawkins2000` (title mashup), `latremoliere2025_hypertopology` (title paraphrase), `mondino_samann2025` ("Synthetic" extraneous), `connes1995` (in-text label 1995→1994), `farsi_latremoliere2024` (no ID, unverifiable → flagged).
- Earlier (v4.67.1): `latremoliere2018`→TAMS 368/2016, field-guide `glanois`→JNT 160/2016, `connes_vs2021` p.2067, D6 Paper-49 metric-level scoping, CLAUDE.md §2 κ line→Observation.

All verified against primary sources (arXiv abstract fetches + web).

**(D) Citation re-delta = CLEAN + calibrated → CERTIFIED (v4.67.3).** 2/2 fresh seeds caught, specificity clean, all 7 remediated cites re-verified correct, deterministic 7/7. One genuine non-seed NIT (`minguzzi_suhr2024` wrong title, *correct referent* — fixed on sight). Zero MATERIAL. Per PI direction, a clean calibrated citation delta = the certifying pass.

**(E) §2 compaction round 6 (v4.67.4).** With the sweep complete, CLAUDE.md §2's ~30 accumulated `/qa` bullets collapsed to the 7 CERTIFIED one-liners + a rollup pointer; 26 superseded per-run bullets moved verbatim to `docs/development_frontier_archive.md`. CLAUDE.md 152 KB → **140 KB** (repo-health WARN cleared).

## 3. Calibration

| Run | Sensitivity | Specificity |
|---|---|---|
| 1st FULL (claims-FG, cross-group, citation) | 6/6 | clean |
| Field-guide re-run | 2/2 | clean |
| Citation delta | 2/2 | clean |
| Citation re-delta | 2/2 | clean |

Every dispatched agent caught its plant; no known-good control was ever flagged MATERIAL. The panel demonstrably detects dirt — including a real citation drift a branch cert had passed over.

## 4. Files changed (papers/docs only — no production code)

- `papers/synthesis/geovac_field_guide.tex` — glanois/deligne… no: glanois JNT, connes_vs2021 page, D6 Paper-49 scoping (v4.67.1); (field-guide PI edits were pre-existing, committed in v4.67.1).
- `papers/synthesis/group1_operator_algebras_synthesis.tex` — latremoliere2018, connes1995 (v4.67.1); hekkelman2022, hekkelman_mcdonald2024, latremoliere_metric_st_2017, hawkins2000, latremoliere2025, mondino_samann, farsi_latremoliere flag (v4.67.2); minguzzi_suhr2024 title (v4.67.3). Recompiles ERRORS=0.
- `CLAUDE.md` — §2 κ-line→Observation (v4.67.1); §2 synthesis bullet + version bumps; §2 compaction (v4.67.4).
- `CHANGELOG.md`, `docs/qa/synthesis.done.md` (DoD → CERTIFIED), `docs/development_frontier_archive.md` (compaction).
- `memory/branch_qa_sweep_phase.md` — sweep COMPLETE.

## 5. Standing lessons

1. **Worktree-sync before seeding (the setup-error lesson).** For any `/qa` target with uncommitted in-scope edits, SYNC the working-tree files into the seed worktree before planting — a worktree-from-HEAD silently reviews stale text and de-validates the dimension (the same class as the group3 path-pin lesson, inverted). Recorded in `docs/qa/synthesis.done.md` + `memory/branch_qa_sweep_phase.md`.
2. **The cross-cutting synthesis cert earns its keep.** It caught wrong-arXiv-ID drift (two IDs resolving to a black-hole paper and an OpenMP paper) that the group1 *branch* cert's citation pass had missed. A branch cert is not a permanent guarantee for its bibliography; a full external-citation audit of all group syntheses is a standing candidate (PI-deferred).
3. **A "clean delta" can still surface one more NIT each cycle** on a large drifted bibliography (delta: 3 LARGE; re-delta: 1 NIT) — the convergence is strong but not instantaneous; the NIT carve-out (correct-referent metadata typo = fix-on-sight) is the stop condition, per the group2 thin-residual-asymptote precedent.

## 6. Honest scope

- **Theorem grade / closed:** nothing new proven this sprint — it is a QA-certification + citation-remediation sprint. The certified *content* (WH1 PROVEN, the two-layer split, the Lorentzian degeneracy theorem, etc.) was proven in prior sprints; this sprint verified the field guide + syntheses *report* it at the correct tier.
- **Numerical/verified observation:** all citation corrections are primary-source-verified (arXiv abstract fetches for hekkelman2022/2206.13744, hekkelman_mcdonald/2403.18619, latremoliere_metric_st/1811.10843, minguzzi/2209.14384; web for glanois JNT 160, latremoliere2018 TAMS 368). The re-delta calibration (2/2 + clean specificity) is the measured discrimination.
- **Structural sketch / flagged, NOT closed:** `hekkelman_mcdonald2024` — the real arXiv ID for "Spectral truncations of Tᵈ and quantum metric geometry" could not be located via search; the wrong ID was removed and the entry flagged "pending author confirmation." `farsi_latremoliere2024` — title "Continuity of spectral propinquity for spectrally truncated spheres" does not resolve to any confirmable paper; flagged "pending confirmation," retained as prior-art lineage credit (not fabricated, not verified).
- **Certification ceiling:** "CERTIFIED" = survived the calibrated detectors across the exercised dimensions (claims re-run, citations re-delta, synthesis-faithfulness, deterministic). Groups 2–6 synthesis bibliographies got only *light spot-checks* this cert (they passed their own branch certs); the field guide + group1 got the deep scan. Code/test-backing was WAIVED (inherited).
- **Named open follow-ons:** (1) full external-citation audit of all group syntheses (PI-deferred); (2) `hekkelman_mcdonald2024` + `farsi_latremoliere2024` bibliographic confirmation (author-level); (3) the standing debug/-refs debt (§9 / C14 advisory).

## 7. Process artifacts

- Version arc: `v4.67.1` (remediation checkpoint) → `v4.67.2` (group1-bib citation fixes) → `v4.67.3` (CERTIFIED) → `v4.67.4` (§2 compaction). All committed + tagged + pushed to `origin/main`.
- No `/regression` run: v4.67.1–.4 touched no `geovac/` production code (the last production change was `thermal_tensor_triple.py` in v4.67.0/group6, verified then with the 190-test group6 suite green + 42 thermal/kg tests).
- Hard-prohibitions (§13.5): none touched — the κ-line edit *removed* a stale "derivable" overclaim (→ Observation); K = π(B+F−Δ) stayed Observation throughout (C12 green).
