# Sprint G6-FP — Fierz-Pauli decomposition of the (1,1) graviton subspace

**Date:** 2026-05-31
**Version:** v3.34.0
**Verdict:** POSITIVE-STRUCTURAL. The diagonal-SU(2) decomposition (1,1) = J=0 ⊕ J=1 ⊕ J=2 is clean, but the spectral action S^(2) is J-blind: identical eigenvalues across all three sectors (proven by Schur's lemma). The Fierz-Pauli decomposition requires the metric identification (δD ↔ δg), not just the spectral action.

## 1. Summary

The (1,1) irrep of SO(4) = SU(2)_L × SU(2)_R decomposes under the diagonal SU(2) as:

  (1,1) = J=0 (trace, 1 dim) ⊕ J=1 (antisymmetric, 3 dim) ⊕ J=2 (graviton, 5 dim)

We computed S^(2) analytically within each J-sector at n_max = 1, 2, 3, 4. The headline finding:

**Theorem (J-blindness).** The spectral-action second variation S^(2) has identical eigenvalue spectra across J = 0, 1, 2 within the (1,1) subspace. Cross-sector eigenvalue multiplicities scale as 1:3:5 = dim(J=0):dim(J=1):dim(J=2).

**Proof.** S = Tr(f(D²/Λ²)) is invariant under unitary conjugation, hence under SO(4) ⊃ diag SU(2). S^(2) is therefore equivariant. By Schur's lemma, S^(2) is proportional to the identity on each irreducible J-component within a given sector pair (a,b). The proportionality constant is w(λ_a, λ_b) — the analytical weight depending only on eigenvalues, not angular quantum numbers. □

## 2. Dimensions

| n_max | dim H | J=0 | J=1 | J=2 | Total (1,1) |
|:-----:|:-----:|:---:|:---:|:---:|:-----------:|
| 1 | 16 | 6 | 14 | 26 | 46 |
| 2 | 40 | 16 | 40 | 72 | 128 |
| 3 | 80 | 26 | 66 | 118 | 210 |
| 4 | 140 | 36 | 92 | 164 | 292 |

Within-sector / cross-sector split at n_max=3:
- J=0: 6 within, 20 cross
- J=1: 6 within, 60 cross
- J=2: 18 within, 100 cross

Within-sector multiplicities scale as 1:1:3 (not 1:3:5) because Hermitianization (V+V^T)/2 reduces the antisymmetric content: each within-sector (1,1) block contributes 1 J=0 + 1 J=1 + 3 J=2 = 5 modes (not the full 9).

## 3. S^(2) eigenvalue structure (n_max=3, Λ²=6)

All J-sectors share the same eigenvalue set. Representative (J=2):

| S^(2) eigenvalue | Mult (J=2) | Laplacian ||[D,V]||² | Type |
|:-----------------|:-----------|:---------------------|:-----|
| +0.133 | 6 | 0 | Within-sector n=2 |
| +0.127 | 6 | 0 | Within-sector n=1 |
| +0.121 | 20 | 4 = (2)² | Cross Δn=2 |
| +0.096 | 20 | 16 = (4)² | Cross Δn=1 (n=0↔1) |
| +0.066 | 6 | 0 | Within-sector n=3 |
| -0.025 | 20 | 36 = (6)² | Cross Δn=1 (n=1↔2) |
| -0.074 | 20 | 16 = (4)² | Cross (gauge orbit) |
| -0.159 | 20 | 64 = (8)² | Cross (gauge orbit) |

J=0 and J=1 have identical eigenvalue values with multiplicities scaled by 1/5 and 3/5 respectively (for cross-sector) or 1/3 and 1/3 (for within-sector).

## 4. Analytical S^(2) formula

The weight w(λ_k, λ_l) in S^(2)[V] = Σ_{kl} w_{kl} |V_{kl}|² is:

**Same eigenvalue (λ_k = λ_l = λ):**
  w = e^{-λ²/Λ²} (4λ²/Λ⁴ − 2/Λ²)

**Different eigenvalues (λ_k ≠ λ_l):**
  w = (λ_k+λ_l)²/Λ⁴ · φ₂(λ_k²/Λ², λ_l²/Λ²) − (a_k+a_l)/Λ²
  where φ₂(x,y) = (e^{-x} − e^{-y})/(y−x), a_i = e^{-λ_i²/Λ²}

Validated against 5-point stencil at n_max=2: max relative error 1.24×10⁻⁶.

## 5. Physical interpretation

The J-blindness theorem establishes a clean scope boundary:

**What the spectral action S^(2) provides:**
- Existence of positive kinetic eigenvalues in the graviton sector (G6-Full)
- Bit-exact integer kinetic spectrum (2k)² (discrete Laplacian)
- Correct sign structure: positive for physical, negative for gauge (SD A₁)
- 2% approach to Lichnerowicz ratio 13/6 at Λ²=4 (G6-Full)

**What the spectral action S^(2) cannot provide:**
- Distinction between TT (2 DOF), longitudinal (3 DOF), and trace (1 DOF) within (1,1)
- The Fierz-Pauli mass structure (m² splitting between J sectors)
- The identification of which modes are physical gravitons vs conformal modes

**What provides the missing structure:**
The metric identification δD ↔ δg_μν, which is the inner-fluctuation / Connes-Chamseddine dictionary. This map breaks the SO(4) symmetry to SO(3) (isometry group on the tangent bundle) and distinguishes symmetric (metric) from antisymmetric (non-metric) perturbations. On the discrete substrate, this identification emerges in the propinquity limit (Paper 38).

## 6. Relation to prior diagnostics

The J-blindness theorem unifies three earlier findings:
- **GD-5 (irrep-blindness):** A_λ is the same for all irreps within a sector — special case of J-blindness for within-sector modes.
- **GD-1 (9→2 as continuum phenomenon):** The polarization reduction requires the metric identification, which is a continuum-limit construction — consistent with J-blindness (spectral action alone doesn't reduce DOF within (1,1)).
- **GD-4 (helicity question):** Whether helicity-0 or helicity-2 modes carry the graviton — irrelevant at the S^(2) level since J-blindness makes all helicities equivalent.

## 7. Honest scope

**Reached:**
- Clean J=0,1,2 decomposition with exact dimensions ✓
- J-blindness theorem with Schur's lemma proof ✓
- Analytical S^(2) formula (validated to 1.2e-6) ✓
- Clear scope boundary: spectral action vs metric identification ✓

**Not reached (requires metric identification / propinquity):**
- Isolation of 2 TT graviton polarizations from 5 traceless-symmetric modes
- Fierz-Pauli mass structure (if any)
- Graviton propagator with correct residue
- Linearized Einstein equations from spectral action + metric identification

**Structural closing of G6 at the spectral-action level:** The J-blindness theorem shows that the spectral action S^(2) has extracted ALL the information it can from the (1,1) subspace. Further graviton physics requires the metric identification, which is a propinquity-level construction. The "multi-month G6" is therefore not a deeper spectral-action computation — it's a propinquity / CC-dictionary task.

## 8. Files

- `debug/g6_fierz_pauli.py` — computation driver
- `debug/data/g6_fierz_pauli.json` — eigenvalue data at n_max=1..4
- `debug/g6_fierz_pauli_memo.md` — this memo
