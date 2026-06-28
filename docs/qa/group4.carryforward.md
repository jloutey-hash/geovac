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

**Disposition for group4 QA:**
1. Decide: keep the pair-diagonal approximation (then **disclose** it and re-state the multipliers honestly) or switch production to the global rule (denser, re-prices the headline). *(Pending PI decision; the sweep below shows option A — disclose — is well-supported.)*
2. ~~Quantify the re-pricing across the library~~ **DONE 2026-06-28** — see below + `debug/qa/group4_cf1_library_sweep_memo.md`.
3. The *energy* deltas for BeH₂/H₂O need a tractable correlated method (full number-projected FCI is intractable at the composed active space — that was the 2 h hang in the diagnostic), but the energy axis is a chemistry footnote, not the group4 concern.

### CF-1 QUANTIFIED — library-wide re-pricing sweep (2026-06-28, group4 pre-work)

Driver `debug/qa/group4_cf1_library_sweep.py` (build-only, no FCI) over all 35 composed library molecules. **The re-pricing factor is constant within valence class** (the dropped m-swaps are a valence-l angular property, not molecule geometry):

| Valence class | N_Pauli/Q pair-diag (production) | N_Pauli/Q global-M_L (physical) | re-pricing |
|:--------------|:-------------------------------:|:-------------------------------:|:----------:|
| main-group s/p (23 mols) | **11.10** | **27.9** | **2.51×** |
| d-block TM hydrides (10 mols) | **9.23** | **30.0** | **3.25×** |

Two material consequences:
1. **Scaling survives, prefactor re-prices.** O(Q^2.5) / "N_Pauli = 11.10·Q exact" / linear-in-Q are ROBUST (constant ratio ⇒ unchanged log-log slope). What re-prices: the multiplier magnitude, and the absolute market-test line — **Paper 20's "LiH 334 vs STO-3G 907" → 837 vs 907 ≈ parity** under the physical rule (does not survive); Paper 14's "51×–1712× vs Gaussian" → ≈20×–680× (one-to-two orders; qualitative advantage survives).
2. **The "d-orbitals are sparser" claim REVERSES.** Papers 14 + 20 both say d-block 9.23 < main-group 11.10 "due to more restrictive Gaunt rules." Under the physical rule d-block is **30.0·Q — denser than main-group 27.9·Q**: the low pair-diagonal 9.23 is an artifact of dropping MORE m-swaps at higher l (3.25× vs 2.51×), not sparser physics. MATERIAL for disposition.

Library-size note (**DECIDED 37, ship as-is — PI 2026-06-28; applied corpus-wide**): `_SYSTEM_REGISTRY` ships **37 systems** (35 composed molecules + `he` + `h2`). Five corpus numbers had disagreed (docstring 28 / Paper 14 30 / registry 37 / Paper 20 38 / CLAUDE.md 40); the 3 organics CH₂O/C₂H₂/C₂H₆ Paper 20 listed are **not buildable** (bond-length entries only, no registry entry). All synced to 37: Paper 14, Paper 20 (organics removed from abstract + multi-center table 8→5), synthesis, docstring, test comment, claims_register, CLAUDE.md §2+§1.1+§1.5. LiH "333" in the sweep is the raw-builder count; the shipping `hamiltonian()` API gives **334** (matches Paper 20 + CLAUDE.md) — 334 is correct, no fix.

**Evidence:** `debug/qa/group3_density_diagnostic_memo.md` (LiH result + sanity), `debug/qa/group3_followup_density_layers_memo.md` (mechanism + full doc mapping), driver `debug/qa/group3_density_diagnostic.py`.

**Group3-side disposition (done / in progress):** the foundations-side fix is just the D-vs-D_pd relabel (universal selection-rule density D = 6.06% at l_max=3; production-realized pair-diagonal D_pd = 1.44%) in Papers 22/31, the group3 synthesis, the claims register, and CLAUDE.md §1.6 (the §1.6 edit is PI-only).
