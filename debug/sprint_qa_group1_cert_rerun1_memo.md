# Sprint memo — `/qa group1` CERTIFYING re-run #1, Batch 1 (cert1)

**Date:** 2026-06-21
**Version:** v4.35.0
**Target:** group1 Batch 1 — Papers 45, 46, 47, 48, 49 + group1 synthesis
**Verdict:** FAIL (calibrated 5/5 sensitivity, 5/5 specificity) → remediated. NOT yet a certified PASS.
**Seed key:** `debug/qa/group1_cert1_seed_key.json`

The first *clean* certifying re-run of the Lorentzian batch — the goal being to convert the rr1/rr1b FAIL→remediated into a certified PASS. It did not certify: the calibrated panel surfaced genuine defects beyond the planted seeds. This is the convergence reality (a fresh adversarial pass with different focus finds new genuine defects); certification needs another clean iteration.

## 1. Panel & calibration

16 agents (claims ×5, citation ×5 exhaustive, code ×5, synthesis ×1), path-pinned to worktree `../geovac-qa-seed-group1-cert1` (base HEAD 6aa0f54), forbidden from the real corpus. Worktree removed; real-corpus seed-leak scan CLEAN.

5 planted seeds, all caught (sensitivity 5/5): S-claims-C14 (P45 abstract)→claims-45; S-claims-C7 (P48 l.1364)→claims-48; S-citation-C4 (P49 che_perales title)→citation-49; S-code-C2 (test_p45 T3 `<1e9`)→code-45/46/48; S-synthesis-C9 (synth "re-established")→synthesis. 5 controls, none false-flagged (specificity 5/5). **Panel calibrated ⇒ verdict trustworthy.**

Deterministic dimension (C5/C11/C12/C13/C14/C15) PASS before and after remediation.

## 2. Genuine MATERIAL defects remediated (real corpus)

- **p46 C14 LARGE** — the informal Main theorem (l.320–339) asserted the withdrawn strong-form Lorentzian quantum-GH convergence as a live Latrémolière quantum-metric distance (`Λ_prop ≤ … → 0`), repaired only by the distant formal theorem / Status note. Added inline **(Descoped.)** note (degenerate Krein seminorm ⇒ rate-formula, not a metric; convergence open). + l.1424 "the present convergence theorem" → "the present (descoped) rate-formula result".
- **p48 C14** — §1.2 overview (l.510–514) and §3 opener (`sec:krein_ppqms`) asserted "all nine Latrémolière axioms transport at theorem-grade rigor" with no inline descope, contradicting the Status note. Added degeneracy/descope notes at both sites (categorical/structural lift, not a metric PPQMS; metric-level descoped).
- **p47 C1** — thm:outer and thm:three_carriers (the PARTIAL paper's surviving keystones) had no test + no matrix row; three_carriers wasn't tagged proof-by-argument. Tagged the three_carriers proof analytic; registered both + the rest of the cluster's surviving content in `docs/claim_test_matrix.md` (new Batch-1 section).
- **C4 ×3 (web-verified):** `farsi_latremoliere2024` was a fabricated "crossed products by amenable group actions, JFA 286 (2024) 110293" in p48 **and** p49 (no such F–L paper; the real 2024 F–L work is "Collapse in NCG and spectral continuity," arXiv:2404.00240 — what p46 already had) → repointed both to Collapse; `kubota2026` "T. Kubota" → "H. Kubota" (arXiv:2605.09101); `bertozzini…` title → "Modular Theory, Non-Commutative Geometry and Quantum Gravity" (SIGMA 6 (2010) 067).

## 3. Reconcile catches (NIT / not acted as MATERIAL)

- p46 F1 temporal-Lipschitz-invisibility is a genuine algebraic identity (proof-by-argument); registered as such in the matrix (panel Frob=0 is illustration). No prose retraction needed.
- p45 the cited `test_p45_kplus_degeneracy.py` T3 is the looser companion to the tight `<1e-9` annihilation in `test_lorentzian_toeplitz_kplus.py` — the keystone IS tightly backed (non-cited sibling). (In the real corpus T3 = `< TOL`; the `<1e9` was the planted seed.)
- p47 Reed–Simon §XIII.16 section-pointer imprecision (book correct, result proved self-contained) — NIT.
- Cite-key/year cosmetic mismatches across the cluster — NIT.

## 4. Verification

- p46/47/48/49 compile errors=0/undef=0; p45 + synthesis unchanged in the real corpus (their only findings were the planted seeds).
- C5/C11/C12/C13/C14/C15 group1 PASS.
- Worktree removed; no seed leaked.

## 5. Honest scope

This is **FAIL→remediated**, not a certified PASS. The genuine defects (p46/p48 descope zombies, p47 unregistered keystones, 3 citation fabrications/slips) were residuals the rr1/rr1b passes (different reviewer focus) did not surface. A clean calibrated re-run with zero genuine material defects is still required to certify Batch 1. Batches 2 and 3 also still need their clean cert re-runs.

## 6. cert2 (re-run #2, v4.35.1)

Second clean cert re-run after the §1–5 (v4.35.0) remediation. 16-agent panel, **calibrated 5/5 sens (5 FRESH seeds, different sites) / 5/5 spec**. Seed key `debug/qa/group1_cert2_seed_key.json`. Worktree removed, no leak.

**The v4.35.0 cert1 fixes HELD — none re-flagged.** The claims reviewers explicitly recognized the p46/p48 inline descope notes and the p47 proof-by-argument matrix registrations as the correct state (specificity confirmed on the remediation itself, not just the controls).

**Verdict: FAIL (calibrated) → remediated — converging.** Exactly ONE genuine new defect (vs cert1's four):
- **p49 §11 C4 author misattribution (SMALL):** prose (l.2406–2407) credited arXiv:2007.00230 ("Gravity Dual of Connes Cocycle Flow") to "Bousso, Casini, Fisher, Maldacena" → web-verified **Bousso, Chandrasekaran, Rath, Shahbazi-Moghaddam** (bibitem already correct). Non-load-bearing positioning §. Prose fixed; stale cite key left as NIT.

NITs (not acted): connes_rovelli page 2918→2917; kubota title leading "A"; p47 §7.5 "Λ_prop" labels the SU(2) γ-rate; p47 latremoliere2018 §5-AF inline pointer; stale cite-key years.

**Convergence:** cert1 = 4 genuine (2 LARGE descope zombies + 3 citations) → cert2 = 1 SMALL non-load-bearing citation. One more clean calibrated re-run (zero genuine defects) certifies Batch 1.
