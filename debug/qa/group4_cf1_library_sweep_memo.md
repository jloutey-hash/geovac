# CF-1 library-wide re-pricing sweep — memo (group4 QA pre-work)

**Date:** 2026-06-28
**Driver:** `debug/qa/group4_cf1_library_sweep.py` (READ-ONLY; runtime monkeypatch of
`composed_qubit._ck_coefficient`, restored after each build; BUILD-ONLY, no FCI).
**Data:** `debug/qa/group4_cf1_library_sweep.json`.
**Predecessor:** `debug/qa/group3_density_diagnostic_memo.md` (LiH-only, +FCI),
`docs/qa/group4.carryforward.md` (CF-1 statement).

## Task

The PI directed (group4 pre-work): *"Quantify CF-1 re-pricing (library sweep)"* /
*"Quantify the full library first, then decide."* CF-1 = the composed Pauli-sparsity
headlines were computed under the **pair-diagonal** ERI rule (`q = mc - ma`, keeps only
m_a=m_c ∧ m_b=m_d), which silently drops genuinely-nonzero m-swap ERIs vs the physical
**global-M_L** Coulomb rule (m_a+m_b=m_c+m_d, via `casimir_ci._gaunt_ck`). Quantify the
N_Pauli re-pricing across the whole composed library, build-only.

## Result — the re-pricing factor is constant within valence class

35 composed library molecules built under both rules. N_Pauli stays **exactly linear in Q
under BOTH rules**; only the coefficient moves:

| Valence class | molecules | N_Pauli / Q (pair-diag, production) | N_Pauli / Q (global-M_L, physical) | re-pricing |
|:--------------|:---------:|:-----------------------------------:|:----------------------------------:|:----------:|
| main-group (s/p) | 23 (LiH…BaH, CO/N2/F2/LiF/NaCl) | **11.10** | **27.9** | **2.51×** |
| d-block (TM hydrides) | 10 (ScH…ZnH) | **9.23** | **30.0** | **3.25×** |

- Sweep stats: n=35 built, ratio min 2.51 / median 2.51 / max 3.25 / mean 2.72.
- **LiH: 333 → 837 (2.51×)** — reproduces the CF-1 memo's LiH measurement bit-for-bit.
- The constancy is structural: the dropped m-swaps are a property of the **valence-l angular
  content**, not the molecule's geometry, so every s/p molecule drops the same fraction.

## Two material consequences for the group4 headlines

1. **The O(Q^2.5)/linear-in-Q *scaling* survives untouched; only the *prefactor* re-prices.**
   Because the ratio is a constant within each class, the log-log slope is unchanged — the
   "N_Pauli = 11.10·Q exact" / "O(Q^2.5) universal" structural headlines are ROBUST. What
   re-prices is the *magnitude of the multiplier* and any *absolute* market-test line:
   - **LiH market test (Paper 20 abstract):** "334 Pauli vs STO-3G 907" → **837 vs 907 ≈
     0.92× (parity)** under the physical rule. This specific line does not survive.
   - **"51×–1712× / two-or-more orders of magnitude vs Gaussian" (Paper 14):** divide the
     main-group multiplier by 2.51 → still ≈ 20×–680× (one-to-two orders). The *qualitative*
     sparsity advantage survives; the headline *magnitude* re-prices by the constant factor.

2. **The "d-orbitals are sparser" claim REVERSES under the physical rule.** Both Paper 14 and
   Paper 20 abstracts state the d-block coefficient (9.23) is *lower* than main-group (11.10)
   "due to more restrictive Gaunt selection rules at higher angular momentum." Under the
   physical global-M_L rule the d-block coefficient is **30.0·Q — DENSER than main-group's
   27.9·Q.** The pair-diagonal 9.23 is low precisely because higher-l blocks have MORE m-swap
   ERIs to drop (3.25× vs 2.51×), not because the physics is sparser. This is a
   pair-diagonal artifact, not a selection-rule property — MATERIAL for the CF-1 disposition.

## Disposition options (for the PI, group4 cert)

CF-1 carryforward step 1 is a strategic choice; the sweep makes it cheap to state either way:
- **(A) Keep pair-diagonal, DISCLOSE it.** Add one sentence + the constant re-pricing factor
  (2.51× main-group / 3.25× d-block) to Papers 14/20, drop/retract the LiH-vs-STO-3G
  parity-losing market-test line, and qualify the "d-block sparser" claim as pair-diagonal-
  specific. Cheapest; preserves the (robust) scaling headlines; honest about the prefactor.
- **(B) Switch production to global-M_L.** Re-prices every headline number (11.10→27.9 etc.),
  denser Hamiltonians, larger but physically-faithful. Heavier; touches `composed_qubit` +
  every downstream count + the 40-molecule library tables.

The sweep shows (A) is well-supported: the re-pricing is a single disclosable constant per
valence class, and the scaling structure (the load-bearing Paper 14/20 result) is unaffected.

## Library-size reconciliation (resolved 2026-06-28)

`ecosystem_export._SYSTEM_REGISTRY` = **37 buildable systems** via `hamiltonian()`: `he`,
`h2`, 18 main-group hydrides, 5 multi-center (LiF/CO/N2/F2/NaCl), 10 TM hydrides (ScH–ZnH),
SrH, BaH. (The 35 composed molecules swept above + `he` + `h2` = 37.) The five corpus
numbers disagreed: ecosystem docstring **28** (stale: said "14 main-group hydrides"),
Paper 14 **30**, code registry **37**, Paper 20 **38**, CLAUDE.md §1 **40**. **The 3 organics
CH2O/C2H2/C2H6 that Paper 20 lists as library members are NOT buildable** — they have only
bond-length entries in `molecular_spec.py:733`, no registry entry; `37 + 3 = 40` is why
CLAUDE.md's "40" is aspirational. Fixed this pass: the ecosystem docstring → 37 (correct
breakdown); the `test_ecosystem_export` comment → 37. **Flagged to PI** (canonical-count
decision touches the non-editable §1.5 "40-molecule library" headline): adopt 37 (drop the
organics overclaim) or register the 3 organics to reach a genuine 40; Paper 14 "30" / Paper 20
"38" + organics sync to whichever the PI picks.

## LiH 333 vs 334 (resolved 2026-06-28 — NOT a defect)

The sweep's "LiH 333" is the **raw `build_composed_hamiltonian` count** (excludes the
identity/constant term). The shipping `ecosystem_export.hamiltonian('LiH')` API gives
**334**, matching Paper 20 + CLAUDE.md. So **334 is canonical and correct**; the re-pricing is
equivalently 334→838 (ratio 2.51× unchanged). No paper fix needed.

## Honest scope

- Build-only Pauli counts; no energy axis (FCI intractable at composed active space — the
  CF-1 memo's 2 h hang; energy effect is sub-kcal/mol per the LiH FCI in the group3 memo,
  a chemistry footnote not the group4 concern).
- The sweep covers the 35 composed multi-center molecules where CF-1 applies; atomic systems
  (O(Q^3.15), Paper 14) use the atomic path and are not re-priced by this rule.
