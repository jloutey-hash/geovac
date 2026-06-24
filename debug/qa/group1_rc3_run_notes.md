# /qa group1 WHOLE-GROUP re-cert #2 (hardened detectors) — run notes (2026-06-23, v4.48.3)

Worktree from 5110a4a (v4.48.2). Deterministic layer (whole group, first): ALL PASS (C5/C10/C11/C13/C14/C15/C16).

## Calibration: PERFECT — 10/10 seeds caught, specificity clean.
- CODE-p46 ✓, CODE-p53 ✓; CLAIMS-A ✓, CLAIMS-B ✓, **CLAIMS-C ✓ (the seed it MISSED last run — the generalized per-paper C14 rule WORKED)**, CLAIMS-D ✓; CITE-A ✓, CITE-B ✓, CITE-C ✓; **SYNTH ✓ (recovery held)**.
- Both recovery tests passed. All controls (v4.48.2 fixes + priors) held — no regression.

## Dimension verdicts
- **Claims / Synthesis / Citation: MATERIAL-CLEAN** (zero non-seed defects — first time; the detector hardening converged them).
- **Code: 4 genuine SMALL non-seed defects, all the same code-docstring/cross-paper-consistency "stale-prose one-scope-over" class** (no logic/tolerance/false-positive defects beyond seeds):
  1. p45 lorentzian module SIBLING docstrings (class LorentzianPropinquityBound l.466 + compute_* l.534 + table l.634) — retracted convergence untagged (v4.48.2 fixed only the module docstring).
  2. p29 test_paper29_int_alg.py docstring/comment (l.16,60) — stale "4s²+1 worked example" (v4.48.2 fixed the body).
  3. gh_convergence.py (l.26,495) + gh_convergence_tensor.py (l.1112) — state-space GH mislabeled "Latrémolière propinquity".
  4. p48 T5 (l.2028) — stale "Π_W-odd per input (I2)" (Paper 43 corrected to mixed-parity).
  + test-hygiene: test_paper50_f_theorem.py module-global mp.dps=60 polluted in mixed batches → spurious failures.

## VERDICT: FAIL → remediated (NOT certified).
All 5 fixed comprehensively (swept whole lorentzian module + gh_convergence pair, not just flagged lines): inline WITHDRAWN tags on the 3 sibling docstrings; p29 docstring→2s²+1; "state-space GH" relabels; p48 T5 mixed-parity; test_paper50 autouse dps fixture. C16 PASS; 47 passed (incl. mixed-batch dps check); p48 errors=0. Worktree removed, leak-scan clean.

## TRAJECTORY (3 whole-group runs)
- v4.48.0: 5 MATERIAL (LARGE+SMALL), synthesis under-calibrated.
- v4.48.2: 3 SMALL, synthesis RECOVERED, claims-C under-calibrated.
- v4.48.3: 4 SMALL (all code-docstring class), claims-C RECOVERED, **calibration PERFECT 10/10**, claims/synth/citation MATERIAL-CLEAN.
Monotone convergence: severity LARGE→SMALL→SMALL-docstring-only; calibration 9/10→9/10→10/10; LLM-judgment dimensions → CLEAN.

## CERT-BAR DECISION POINT (PI)
Residual is entirely mechanical code-docstring stale-prose (not science, not judgment, not in papers). Strict /qa PASS = zero verified MATERIAL ⇒ still FAIL. But: calibrated-perfect, zero-LARGE, claims/synth/citation clean, deterministic-green. Options: (a) iterate strictly to zero (each pass finds ~1 more stale docstring one-scope-over — the asymptote); (b) set the bar at "calibrated + zero-LARGE + claims/synth/citation clean + deterministic-green + code-docstring-stale-prose = fix-on-sight NIT" → certify now (frozen-criteria change, PI-only); (c) build a deterministic docstring-scan gate first, then iterate. Recommendation: (c)-then-(b).
