# Kluth-Litim 2020 cross-check of α-X's S⁵ Seeley-DeWitt values

**Scope.** Independently verify the sympy-computed S⁵ Seeley-DeWitt
coefficients in `debug/alpha_sprint_a/compute_s5_sd_coeffs.py` against
an external published reference.

**One-line verdict.** **MATCH.** All six α-X values (b₀, b₂, b₄ for
the scalar Laplacian and squared Dirac operator on unit S⁵) are
reproduced exactly by Kluth-Litim 2020 Table 1 after (i) dividing α-X
by Vol(S⁵) to match KL's volume-normalized convention, and (ii) using
KL Eq. 13 to transform the scalar-Laplacian reference at E=0 to
E_KL = −R/4 (Lichnerowicz) and multiplying by dim_spin = 4.

---

## 1. Reference

**Kluth, Y. and Litim, D.F.**, *"Heat kernel coefficients on the
sphere in any dimension,"* EPJ C **80**:269 (2020),
arXiv:1910.00543v2 [hep-th]. The paper derives all heat-kernel
coefficients for Laplacians on scalars, transverse vectors, and TT
tensors on fully symmetric backgrounds in any dimension, confirmed
via Euler-Maclaurin spectral sums. Relevant: Sec. II.A (definitions),
Eq. 13 (E-transformation), Sec. III.B + Table 1 (scalar coefficients
for d = 2..6).

**KL does NOT explicitly tabulate the Dirac operator.** For D², the
coefficients were obtained from KL's scalar table via Eq. 13 with
E_KL = −R/4 (Lichnerowicz), then multiplied by dim_spin = 4 (5D
spinor dimension). On a constant-curvature manifold the curvature
endomorphism is scalar, so the heat kernel splits as
dim_spin × scalar-result-at-E=−R/4 — the standard reduction.

---

## 2. KL Table 1, d = 5, scalar Laplacian (E = 0)

| Coefficient | KL value | On unit S⁵ (R_scalar = 20) |
|:-----------:|:--------:|:---------------------------:|
| b₀^(0) | 1 | 1 |
| b₂^(0) | R / 6 | 10/3 |
| b₄^(0) | R² / 75 | 16/3 |

**Normalization.** KL Eq. 11: `Tr U_E = Vol/(4πt)^{d/2} · Σ b·tⁿ`,
so KL's b is **volume-normalized**: b = (1/Vol)·Tr[ã_k]. α-X's a_k
includes Vol explicitly; comparison requires α-X / Vol.

---

## 3. Comparison

**Scalar.** α-X / Vol on unit S⁵ (Vol = π³):

| Coeff | α-X (Vol-divided) | KL value | Match |
|:-----:|:-----------------:|:--------:|:-----:|
| b₀ | 1 | 1 | ✓ |
| b₂ | 10/3 | 10/3 | ✓ |
| b₄ | 16/3 | 16/3 | ✓ |

**Squared Dirac D²** (KL Eq. 13 with E_KL = −R/4, then ×4):

| Coeff | α-X (Vol-divided) | KL-derived (R=20) | Match |
|:-----:|:-----------------:|:----------------:|:-----:|
| b₀ | 4 | 4·1 = 4 | ✓ |
| b₂ | −20/3 | 4·(−R/12) = −20/3 | ✓ |
| b₄ | 14/3 | 4·(7R²/2400) = 14/3 | ✓ |

The b₄ Dirac arithmetic is the nontrivial step:

```
b₄^Dirac_scalar(E=−R/4) = (E²/2)·b₀ + E·b₂ + b₄
                        = R²/32 − R²/24 + R²/75
                        = R²·(75 − 100 + 32)/2400 = 7R²/2400
4 × 7R²/2400 = 7R²/600;  at R=20: 7·400/600 = 14/3.  ✓
```

All six values reproduce α-X exactly as rational symbolic equalities
(see `debug/alpha_sprint_a/compare_s5_sd.py` for the runnable check).

---

## 4. Verdict: MATCH

**α-X's sympy computation is independently verified.** No
substantive normalization discrepancies. The only convention
adjustments are: (a) Vol-normalization (KL bare-b, α-X Vol-included
a) — trivial reclassification, no physics change; (b) sign convention
on Lichnerowicz E_KL = −R/4 (KL's operator is P = −∇² − E_KL, which
reproduces D² = −∇² + R/4). Both conventions are standard
(KL follows Avramidi/Vassilevich; α-X follows the GeoVac S³ module
`geovac/qed_vacuum_polarization.py`). The curvature polynomial
content on unit S⁵ is identical.

---

## 5. Implications for α-EB v2 / α-X interpretation

1. **α-X's sympy values are correct.** Unit-S⁵ scalar SD: b₀=1,
   b₂=R/6, b₄=R²/75. Dirac with Lichnerowicz: b₀=4, b₂=−R/3,
   b₄=7R²/600.

2. **α-X's PARTIAL-NEGATIVE verdict stands unchanged.** The match
   confirms π³ = Vol(S⁵) factors out of every SD coefficient
   (KL Eq. 11 makes this manifest). But the absence of any
   α³-order structural mechanism in the Connes-Chamseddine
   spectral-action expansion on a 5-manifold (odd d → odd Λ powers
   only, no tree-level α³) is **not** contradicted by KL — KL extends
   the same expansion to higher orders (b₆, b₈, b₁₀ at d=5: R³/1500,
   R⁴/45000, R⁵/2250000 — all rational × R^k after normalization),
   none of which produces an α³-multiplier. The numerical coincidence
   R_predict ≈ π³·α³ at 0.25% remains shape-consistent (π³ = Vol(S⁵))
   but derivationally absent.

3. **Paper 2 §IV.G recommendation unchanged.** Defer α-EB v2
   integration beyond a structural-hint footnote. Paper 24 π-free
   certificate intact.

4. **Higher-order extension available.** If a later sprint extends
   α-X beyond a₄, KL Table 1 provides the independent reference
   through b₁₀ in d=5 automatically.

---

**Files.**
- `debug/alpha_sprint_a/kluth_litim_s5_verification_memo.md` (this memo)
- `debug/alpha_sprint_a/compare_s5_sd.py` (verification script, prints MATCH)
- `debug/alpha_sprint_a/kluth_litim_extracted.txt` (raw PDF text, 76 KB)

**Cross-reference.**
- α-X memo: `debug/alpha_sprint_a/paper24_pi_cubed_crosscheck_memo.md`
- α-X data: `debug/data/alpha_sprint_a_s5_sd_coeffs.json`
- α-X script: `debug/alpha_sprint_a/compute_s5_sd_coeffs.py`
- KL: arXiv:1910.00543v2, EPJ C 80:269 (2020)
