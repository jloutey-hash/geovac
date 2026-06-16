# Group 4 (Quantum Computing) — carry-forward findings for QA

Findings surfaced while QA-ing another branch that **land on group4 papers** and must be dispositioned when group4 is certified. Read this before co-writing `docs/qa/group4.done.md`.

## CF-1 — The composed Pauli-count advantage rests partly on the pair-diagonal ERI approximation (from `/qa group3`, 2026-06-16)

**What.** The production composed pipeline (`geovac/composed_qubit.py::_ck_coefficient`, `q = mc - ma`) realizes the **pair-diagonal** ERI selection rule (m_a=m_c AND m_b=m_d), not the physically-correct global-M_L Coulomb rule (m_a+m_b=m_c+m_d). It silently drops genuinely-nonzero m-swap ERIs (e.g. ⟨p₊₁p₋₁|p₋₁p₊₁⟩, angular factor ≈ 0.49). This is a legitimate *sparsifying approximation*, but the group4 headline numbers were computed under it and are currently presented as if they were exact angular-selection-rule sparsity.

**Measured impact (LiH, the one tractable full-FCI case):**
- N_Pauli 333 → **837 (2.51× denser)** under the physical rule.
- The market-test line **"LiH 334 Pauli, 2.7× fewer than STO-3G (907)"** re-prices to **837 vs 907 ≈ 1.08× — parity.**
- Energy effect small (ΔE ≈ 1.2 mHa, sub-kcal/mol) — so this is a *sparsity/advantage* issue, not an accuracy one.
- The O(Q^2.5) *scaling structure* and the potential-independence theorem (Paper 22) are unaffected; what re-prices is the *magnitude of the multiplier* and whether part of it is undisclosed approximation.

**Lands on:** Paper 14 (qubit encoding — 51×–1712× vs Gaussian, LiH-vs-STO-3G market test, O(Q^2.5)); Paper 20 (resource benchmarks).

**Disposition for group4 QA (do NOT do now — group4's bite):**
1. Decide: keep the pair-diagonal approximation (then **disclose** it and re-state the multipliers honestly) or switch production to the global rule (denser, re-prices the headline).
2. Quantify the re-pricing across the 40-molecule library — **cheap: Pauli-count only, build-only, no FCI** (LiH is 2.51×; systems with more same-l valence electrons will re-price harder).
3. The *energy* deltas for BeH₂/H₂O need a tractable correlated method (full number-projected FCI is intractable at the composed active space — that was the 2 h hang in the diagnostic), but the energy axis is a chemistry footnote, not the group4 concern.

**Evidence:** `debug/qa/group3_density_diagnostic_memo.md` (LiH result + sanity), `debug/qa/group3_followup_density_layers_memo.md` (mechanism + full doc mapping), driver `debug/qa/group3_density_diagnostic.py`.

**Group3-side disposition (done / in progress):** the foundations-side fix is just the D-vs-D_pd relabel (universal selection-rule density D = 6.06% at l_max=3; production-realized pair-diagonal D_pd = 1.44%) in Papers 22/31, the group3 synthesis, the claims register, and CLAUDE.md §1.6 (the §1.6 edit is PI-only).
