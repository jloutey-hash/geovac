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

## 7. cert3 (re-run #3, v4.36.0)

Third clean cert re-run. 16-agent panel, **calibrated 5/5 sens (5 fresh seeds) / 5/5 spec**; cert1+cert2 fixes HELD (none re-flagged). Seed key `debug/qa/group1_cert3_seed_key.json`. Worktree removed, no leak.

**Verdict: FAIL (calibrated) → remediated — convergence NOT monotone (cert1=4, cert2=1, cert3=3).** Three genuine defects, all reconciled before acting:
- **p48 C14 LARGE (Bridge-Theorem B4/T2 descope zombie):** the abstract Bridge summary ("all four bridge properties hold at theorem-grade rigor") + `thm:convergence_transport` (B4) + `thm:synthetic_compactness` (T2) asserted Mondino–Sämann pLGH *convergence* as established, while the premise $\BigDeth\to0$ rides on the degenerate Krein metric $L^K$ (Status note descopes T3/T6 but not B4/T2). SAME class as cert1's §3 fix — a deeper Bridge-Theorem read found it. Fixed: inline descope on the abstract B4 item (B1/B3 survive; B4 descoped), `thm:convergence_transport`, `thm:synthetic_compactness` (pre-compactness survives, convergence descoped).
- **synthesis C4 MATERIAL (wrong-ID → different paper):** `bousso_etal2020` = arXiv:2008.03319, which **web-resolves to Akers–Penington "Leading order corrections to the quantum extremal surface prescription"** (a different paper), under wrong authors/title, while the prose describes the kink-transform/Connes-cocycle gravity dual (p49's BCRS paper). Repointed bibitem + prose to "Gravity dual of Connes cocycle flow," Bousso–Chandrasekaran–Rath–Shahbazi-Moghaddam, PRD 102 (2020) 066008, arXiv:2007.00230. The cert3 synthesis reviewer flagged it as "same as p49"; the reconcile web-verify showed it was worse (wrong-ID resolving to a different paper).
- **p46 C8 SMALL:** Appendix-B L1' $\prop_{\mathrm{achievable}}=1$ declarative inside the descoped `thm:enlarged_main` proof → proof-sketch reservation added.

Reconcile note: the synthesis Bousso finding was elevated SMALL→MATERIAL by web-verify (the cert3 reviewer's "same as p49" framing was incomplete; it's a distinct wrong-ID confusion).

**Systemic finding (the load-bearing takeaway):** cert1 (§3) and cert3 (B4/T2) are the same defect class — the Lorentzian cluster (P46/P47/P48) carries many theorem statements whose descope lives only in a Status note, not inline on each theorem. Each fresh deep read finds another un-tagged instance, so blind cert re-runs won't converge efficiently. **Recommended: a dedicated inline descope-tagging sweep of P46/P47/P48** (enumerate every theorem/proposition asserting a metric-level Lorentzian/propinquity convergence; ensure each carries an inline descope tag), analogous to the folded-in C4 citation sweeps, then re-run cert. This drains the tail in one pass.

---

## cert4 (v4.38.0, 2026-06-22) — FAIL→remediated; descope sweep HELD, but prose/abstract zombies remain

Base HEAD ee522a6 (the v4.37.0 descope-sweep state). Calibrated panel: **5/5 sensitivity** (all 5 fresh seeds caught — code C2 vacuous tol, synthesis C9 descoped→established, citation C4 kunzinger title, claims C14 p45 signature-agnostic→Lorentzian, claims C7 p48 Paper-38 propinquity mislabel), **5/5 specificity** (controls M1–M5 + ALL prior fixes + the v4.37.0 descope tags NOT re-flagged — the sweep held).

**Verdict: FAIL→remediated.** The theorem-only descope sweep held perfectly (no theorem-statement zombie re-surfaced), but cert4 found the *next* class — descope zombies in **§1 Introduction prose** and **abstract enumerations**, which the theorem-statement sweep did not cover.

Genuine material defects (all remediated):
- **p47 §1 Introduction (l.216–224), C14:** the program restatement "The Lorentzian propinquity program of Papers 45–46 establishes that the finite-cutoff Lorentzian Krein spectral triple converges, in the Latrémolière propinquity, to the continuum…" stated the DESCOPED metric-convergence as established, no inline tag. Fix: "set out to establish" + explicit "that metric/propinquity-level convergence is descoped (Paper 45 degenerate K⁺ seminorm); only the norm-resolvent / spectral leg survives."
- **p48 thm:compact_agreement (l.1371–1386), code/C8 over-reach:** "the Krein hypertopology reduces bit-exactly to the Paper 38 SU(2) state-space GH hypertopology" over-states the JOINT rate. Per `central_fejer_compact_temporal.py` l.500, `gamma_l1(N_t=1) = gamma_su2 + T/4` — only the operator-system substrate and the SU(2) rate factor reduce bit-exact; the joint L1 rate carries an additive U(1) temporal offset T/4. Statement + proof reworded to distinguish substrate (bit-exact) from joint rate (SU(2) factor + T/4).
- **p48 abstract "Main theorem" (l.176–182), C14:** the B4 item ("convergence transport via Plancherel-nested Berezin compatibility") listed under "holds at theorem-grade rigor" with no descope — a DIFFERENT abstract location from the cert3-fixed §1.2 overview (l.382). Fix: inline "(B4 descoped: rests on the degenerate Krein metric of Paper 45; conditional pending substrate repair)."
- **p49 driver ref (l.2131–2133), C13/C14-advisory:** cited a nonexistent driver `q1prime_phase2b3_panel_compute.py`. Repointed to the permanent backing `tests/test_lorentzian_propinquity_foundation.py` (§9 cite-permanent-records policy). C13 now resolves it live.
- **NIT — `geovac/lorentzian_propinquity_compact_temporal.py` docstring:** stated the RETRACTED K⁺-restricted weak-form propinquity theorem as the live "Paper 45 headline." Added a prominent `*** DESCOPED / THEOREM RETRACTED ***` banner (P45 K⁺ degeneracy; falsifier `tests/test_p45_kplus_degeneracy.py`; what survives = signature-agnostic product-carrier convergence).

OVER-FLAG reconciled (NOT a defect, NOT changed): **che_perales bibitem venue** "Differential Geometry and its Applications, Volume 103 (2026)" — citation-47/49 reviewers flagged it as a fabricated venue on an unpublished preprint; web-verify CONFIRMED the paper ("Gromov's Compactness Theorem for the Intrinsic Timed Hausdorff Distance," Che–Perales–Sormani) IS published in DGA Vol. 103 (June 2026). The arXiv page simply hasn't been updated with the journal-ref. Venue is correct.

Verification: p47/p48/p49 compile errors=0, zero undefined cites/refs; C5/C11/C13/C14/C15 group1 PASS; foundation + p45-degeneracy + wh7-b1 tests 57 passed / 2 slow-skipped.

**Convergence status:** cert4 found 4 new genuine defects of NEW classes (§1-prose + abstract-enumeration zombies, a rate-vs-substrate over-reach, a dangling driver). The cert loop is still finding genuine material per pass — NOT yet a clean certified PASS. Next sweep candidate (if the PI wants to drain this tail too): an inline descope scan of the §1 Introduction + abstract *prose* of P45–49 (the theorem-statement sweep is done; the prose/abstract layer is the remaining systematic home). Then cert5.
