# Sprint memo — full deterministic variance-optimized Jastrow (He), Stage-2 of the accuracy capstone

**Date:** 2026-08-23/24  **Verdict: the physical variance-optimized Jastrow CONFIRMS (does not resolve) the v5.0.7 accuracy↔cost tension.**

**Question (PI-chosen "full deterministic optimizer" build).** v5.0.7 left the accuracy leg open: the cheap non-Herm TC operator is γ-fragile, and the literature fix is a *flexible, variance-optimized* Jastrow (Haupt et al. JCP 158 224105; López Ríos–Alavi JCP 163 084107). Does building that — rigorously — make the cheap TC operator (κ_V O(1), no overlap matrix) chemically accurate on GeoVac's Coulomb-Sturmian basis?

Driver: `debug/tc_flex_variance_he.py` (+ audits `debug/_audit_flexvar.py`, `debug/_phys_flexvar.py`); data `debug/data/tc_flex_variance_he.json`, `debug/data/_phys_flexvar.log`. Fast engine: `debug/tc_flex_jastrow_he.py` (H̃ is exactly quadratic in the Jastrow coefficients).

## What was built (each piece validated)
- **Fast polynomial engine.** The TC 2-body operator is degree-2 in the Jastrow coefficients c (u′ linear ⇒ K linear; w = 1/r − Σc∇²u_μ − (Σc u_μ′)² linear+quadratic). Precompute per-exponent/per-pair asym tensors once; assemble H̃(c) by cheap combination. Reduces to the base engine's `w_kernel` at M=1 (validated: γ=1.0 → +5.7 mHa, matches Stage-1).
- **Well-posed objective (the Stage-1 correction).** The naive internal σ²=Σ|⟨D|H̃|Φ⟩|² is ILL-POSED (Stage-1: spurious σ²→0 at wild |c|~128, garbage energy). The correct object is the **e^{2J}-weighted local-energy variance** (zero-variance principle): using (H−E)e^J=e^J(H̃−E), Var = wᵀW2 w / (ΦᵀW2Φ) with W2 = FCI matrix of the e^{2u(r12)} operator (the VMC weight the naive version dropped), Ew = (ΦᵀW2 H̃Φ)/(ΦᵀW2Φ), w = H̃Φ − Ew·Φ.
- **RHF reference (required).** The bare h1-aufbau det is a terrible reference (⟨Φ|H|Φ⟩ = −1.845 Ha; weight 0.72 in the FCI ground) → E_weighted came out +610 mHa (nonsense). Fix: RHF (E_HF = −2.8616) carried as a CI-VECTOR reference Φ_HF in the Löwdin det basis (transform_1/2 are hard-coded for the symmetric Löwdin matrix, so an orbital-rotation to the HF basis breaks them — stay in Löwdin, use the vector). After the fix the weighted variance has a genuine interior minimum and E_weighted is sane (+18–33 mHa single-geminal).

## Result
| Jastrow | u(0) | shape | Var | E_weighted (VMC) | **transplant dTC** | κ_V |
|---|---|---|---|---|---|---|
| flexible, unconstrained | +0.95 | UNphysical (anti-hole) | 0.0231 | +12.7 mHa | **−0.5 mHa** | 1.2–1.5 |
| flexible, physical (u≤0), β=1..100 | −0.31 | physical hole ✓ | 0.0288 | +17 mHa | **+7.7..+8.2 mHa** | 1.2–1.4 |
| single geminal (best) | <0 | physical | 0.044 | +18 | ~0 (tuned γ) | O(1) |
| **R12-CI, same basis (v5.0.2)** | — | enlarged space | — | — | **0.80 mHa** | (κ(S)≈200) |

- **The −0.5 mHa is a conditioning artifact.** Near-collinear exponents give the (well-posed) variance a flat plateau; unconstrained, it is *minimized* by an oversized-|c|, unphysical **positive-u** Jastrow (Var 0.023 < physical 0.029) that fits an accurate energy by cancellation — the corpus's numerical-luck-at-high-condition-number pattern. Robust across the CMAX sweep (|c| always pins to the bound; Var flat 0.02314→0.02312 over CMAX 4→20).
- **The physical answer is +7.7 mHa** — robust across β=1/10/100 — a real ~3× improvement over plain (+25) at κ_V≈1.2 and no overlap matrix, but NOT chemical accuracy and no better than a hand-tuned single geminal. The physical zero-variance floor (Var≈0.029, not ~0) says the ansatz (single HF det × pair-Jastrow, folded into the operator) cannot represent exact He to sub-mHa in this basis.

## Verdict — tension CONFIRMED and SHARPENED
The sub-mHa accuracy lives in the geminal-**enlarged** wavefunction space (R12-CI = 0.80 mHa on the *same* basis, because it adds the geminal as a full CI *function*). The cheap TC operator deliberately does not carry that space (it folds the geminal into H̃ and keeps the small orbital space) → caps at ~8 mHa. That enlarged space is exactly the one whose *encoding* Track B showed is conditioning-limited (κ(S)≈200). So: **cheap operator OR enlarged-space accuracy, not both — the accuracy is provably in the part that is expensive to encode.** The physicality audit was the decisive step (without it we would have banked a false −0.5 mHa GO).

## Honest scope / not done
Minimal single-common-k s-only Sturmian basis, single HF-det reference, M=5 pair-Jastrow. A richer basis/reference would lower the ~8 mHa floor, but the *structural* separation (enlarged-space accuracy vs cheap operator) is robust and matches Track B. Li/Be + resource estimate not run — they would re-confirm the He conclusion at more cost (3-body terms, more electrons). The fast polynomial engine + weighted-variance + RHF-vector-reference machinery are reusable if the direction is resumed.
