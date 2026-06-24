# /qa group1 — WHOLE-GROUP CONFIRMATION re-cert (post v4.48.0/.1) — run notes (2026-06-23)

Worktree from cad3394 (v4.48.1). Seed key: debug/qa/group1_recert_seed_key.json.

## Deterministic layer (whole group, first): ALL PASS
C5, C10 (compile), C11, C13, C14, C15, C16 — C16 now also scans the geovac lorentzian module (blind-spot fix).

## Calibration scorecard
Sensitivity **9/10 seeds caught**; specificity clean (all F1–F5 + M1–M6 controls held).
- CODE: CODE-p40 ✓, CODE-p50 ✓ (2/2).
- CLAIMS: CLAIMS-A ✓, CLAIMS-B ✓, **CLAIMS-C ✗ MISSED** (p46 "established at the metric level" — reviewer harmonised it with the "(degenerate)" qualifier in the same sentence), CLAIMS-D ✓ (3/4).
- CITATION: CITE-A ✓, CITE-B ✓, CITE-C ✓ (3/3; each returned the correct value).
- **SYNTHESIS: SYNTH ✓ CAUGHT** — the sharpened claims-reviewer flagged the "established K⁺ bound was →0" lead-in DESPITE the retraction 12 lines below. **Calibration RECOVERED** (the v4.48.1 fix worked).

## Specificity — the v4.48.0/.1 fixes HELD
p45 lorentzian_theorem_statement() = K⁺-annihilation (code-45 ✓); p53 plane docstring = heat-attribution (code-53 ✓); p46 C₃=1 / √ withdrawn (code-46 ✓); synthesis eq:p46_main √-withdrawn (synthesis ✓); p29 int_alg backfill GENUINE (code-29 ✓); p42/p43 tripwires clean (claims-B ✓). No control regressed.

## VERDICT: FAIL → remediated (NOT certified).
3 genuine non-seed MATERIAL defects (all SMALL), all fixed:
1. **code-45**: geovac/lorentzian_propinquity_compact_temporal.py module-docstring inner Theorem block (l.44) + priority claim (l.84) stated the retracted convergence theorem as live prose under (not within) the banner → inline WITHDRAWN tags added.
2. **code-29**: paper_29 Cor int_alg worked example + obs:galois cited a NON-EXISTENT factor 4s²+1 (live S⁵ N_max=2 factor is **(2s²+1)²**, verified via ihara_zeta_bass) → corrected to 2s²+1 / s²+2; test_paper29_int_alg.py updated to the real factor + a live-recompute tie-in.
3. **claims-C**: paper_47 §1 G1 gap-list "closed on the natural substrate in Paper 46" → "descoped (rate-formula only, metric-level open)".
+ synthesis Finding-1 (L5 "Latrémolière propinquity assembly") = OVER-FLAG, reconciled NIT (matches Paper 38's own lemma name; prior-reconciled v4.43.3).

## Process/calibration fixes (the durable wins)
- **claims-reviewer.md**: the descope enumerate-and-quote rule GENERALIZED from synthesis-only to per-paper C14 too (the claims-C miss class) — "a nearby hedge does NOT cure a live-sounding clause."
- **C16**: added geovac/lorentzian_propinquity_compact_temporal.py to the propinquity + withdrawn-c3op entries (the code-docstring blind spot code-45 named — C16 had scanned papers/*.tex only).

## Coverage limitation of THIS run
Wave B (rest of code: p39, p42, p43, p44, p49, analytic-trio 47/48/52) NOT re-run this turn + no separate completeness-critic — the verdict is FAIL regardless (genuine defects found), and those papers had no v4.48.0 changes + were clean in the prior whole-group run. claims-C chunk under-calibrated this run (missed its descope seed; now addressed by the prompt fix, needs re-validation).

## Status
group1 NOT yet certified. Convergence: zero LARGE this run, all residuals SMALL, synthesis calibration recovered, two calibration/blind-spot mechanisms now closed. The recurring reality: each fresh whole-group pass surfaces a few new latent SMALL items + occasionally a calibration gap (the "fresh adversaries" lesson). PI decision pending: keep iterating toward a clean PASS, or set the cert bar at "calibrated panel / zero-LARGE / all-SMALL-fixed / deterministic-green."
