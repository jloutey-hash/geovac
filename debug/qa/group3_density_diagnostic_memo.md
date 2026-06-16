# group3 density diagnostic — pair-diagonal vs global-M_L ERI rule (LiH)

**Date:** 2026-06-16. **Driver:** `debug/qa/group3_density_diagnostic.py` (read-only; production code monkeypatched at runtime only, restored after). **Purpose:** decide whether `composed_qubit._ck_coefficient`'s pair-diagonal restriction (q=mc−ma) is intended physics or a latent bug, by measuring what the dropped m-swap ERIs do to energy + Pauli count.

## Sanity cross-check (PASSED — diagnostic is trustworthy)

| l_max | pair-diagonal % | global-M_L % |
|:--|:--|:--|
| 1 | 7.81 | 14.84 |
| 2 | 2.76 | 8.52 |
| 3 | **1.44** | **6.06** |

- Baseline reproduces the corpus's quoted **1.44%** at l_max=3; global rule gives **6.06%**.
- m-swap probe c²(p₊₁,p₋₁): pair-diagonal **−0.0000** (dropped) vs global **−0.4899** (retained). The patch genuinely flips the rule.

## LiH result (FCI completed both passes)

| | baseline (pair-diag) | modified (global-M_L) | Δ / ratio |
|:--|:--|:--|:--|
| E_FCI (Ha, composed-internal ref) | −14.143000 | −14.144222 | **ΔE = −1.22 mHa (−0.0086%)** |
| N_Pauli | **333** | **837** | **2.51× denser** |
| nonzero ERIs | 195 | 321 | 1.65× |

Baseline N_Pauli 333 reproduces the published composed-LiH ~334 → sanity ✓.
(The −14.14 absolute is the composed builder's internal energy, not the physical −8.07 Ha total; only the **ΔE between conventions** is the physics-impact measure.)

## Interpretation

- **The dropped m-swap ERIs are real, but their energy effect is small.** Including them lowers the LiH ground state by **~1.2 mHa ≈ 0.77 kcal/mol — sub-chemical-accuracy.** Sign is physical (dropping attractive terms had raised the energy). So the pair-diagonal restriction is **not exact physics** (it discards genuinely nonzero Coulomb terms) but is a **mild, defensible energy approximation**.
- **The cost lands on the headline sparsity/Pauli advantage, not the energy.** Under the physical rule LiH is **2.51× denser** (333→837 Pauli). The market-test claim "LiH 334 Pauli, 2.7× fewer than STO-3G (907)" re-prices to **837 vs 907 ≈ 1.08× fewer (parity)**. The flagship advantage rests materially on this approximation.

**Verdict:** neither "clean intended exact physics" nor "catastrophic bug" — a **legitimate sparsifying approximation that the corpus currently presents as if it were the exact angular density.** Integrity fix = **disclose + re-price**, not necessarily a code change:
1. Label 1.44% as the *realized/approximate pair-diagonal* density `D_pd`; 6.06% is the exact global-M_L density (Paper 22 / 31 / synthesis / CLAUDE.md §1.6 / claims_register).
2. Caveat the Pauli-count advantage as resting on the pair-diagonal approximation (~2.5× on LiH), and re-state the LiH-vs-STO-3G comparison honestly.
3. Whether to *switch* production to the global rule is a cost/benefit call for the PI (≈1 mHa better energy for ~2.5× more Pauli terms) — likely keep the approximation, but disclosed.

## BeH₂ / H₂O

Full number-projected FCI at the composed active space is **intractable** through this driver — BeH₂'s baseline FCI ran > 2 h without completing (6 e⁻; H₂O's 10 e⁻ worse). The **energy** delta for BeH₂/H₂O needs a tractable correlated method, not full FCI. The **Pauli-count re-pricing** (the advantage question) is cheap (build-only, no FCI) and can be obtained in seconds if wanted.
