# Track SM-F: Higher Hopf Fibrations and Standard Model Couplings

**Sprint:** α / SM connection, wildcard speculative track.
**Date:** 2026-04-10
**Status:** FAIL (as a derivation of SM couplings); PARTIAL (as a set of
well-defined spectral invariants with isolated numerical coincidences).

---

## 1. Setup

Paper 2 extracts α from the complex Hopf bundle S¹ → S³ → S² via three
spectral invariants:

- **B = 42**: degeneracy-weighted Casimir trace on the base S²,
  summed over the S³ Fock shells (n,l), n ≤ 3, l ≤ n−1.
- **F = π²/6**: Dirichlet series of the S³ shell degeneracies
  D_{n²}(d_max=4) = ζ_R(2).
- **Δ = 1/40 = 1/(|λ_{n_max}|·N(n_max−1)) = 1/(8·5)**: finite-size
  boundary term.
- **K = π(B + F − Δ) = 137.036064**, α = 1/K.

The selection principle that picks n_max = 3 is
**B/N = dim(S³) = 3**, uniquely solved by the closed-form
B(m)/N(m) = 3(m+2)(m−1)/10.

The higher Hopf fibrations are

| Division algebra | Fibration | dim base | dim total | fiber |
|:---|:---|:---:|:---:|:---:|
| ℝ | S⁰ → S¹ → S¹ | 1 | 1 | S⁰ |
| ℂ (Paper 2) | S¹ → S³ → S² | 2 | 3 | S¹ |
| ℍ | S³ → S⁷ → S⁴ | 4 | 7 | S³ |
| 𝕆 | S⁷ → S¹⁵ → S⁸ | 8 | 15 | S⁷ |

with expected gauge group assignments (Furey/Dixon/Baez literature):
ℂ → U(1)_em, ℍ → SU(2)_L-like, 𝕆 → color/exceptional.

---

## 2. Computed B-analogs: closed forms

Sum **B(L) = Σ_{l=0..L} dim(V_l) · l(l + N−1)** directly on the base
spherical harmonics (no analog of a Fock shell structure exists on
S⁴ or S⁸; Paper 24's HO rigidity theorem forbids the Fock projection
outside the Coulomb/S³ case).

### Base S⁴ (quaternionic)

Degeneracy: dim(V_l) = (l+1)(l+2)(2l+3)/6. Eigenvalue: l(l+3).

Closed forms (sympy, exact):

  B_4(L) / N_4(L) = **2 L(L+4) / 3**

Table (L = 1..10):

| L | B | N | B/N |
|---:|---:|---:|---:|
| 1 | 20 | 6 | 10/3 |
| 2 | 160 | 20 | **8** |
| 3 | 700 | 50 | 14 |
| 4 | 2240 | 105 | 64/3 |
| 5 | 5880 | 196 | **30** |
| 6 | 13440 | 336 | 40 |
| 7 | 27720 | 540 | 154/3 |
| 8 | 52800 | 825 | **64** |
| 9 | 94380 | 1210 | 78 |
| 10 | 160160 | 1716 | 280/3 |

### Base S⁸ (octonionic)

Closed form:

  B_8(L) / N_8(L) = **4 L(L+8) / 5**

Table (L = 1..10):

| L | B | N | B/N |
|---:|---:|---:|---:|
| 1 | 72 | 10 | 36/5 |
| 2 | 864 | 54 | 16 |
| 3 | 5544 | 210 | 132/5 |
| 4 | 25344 | 660 | 192/5 |
| 5 | 92664 | 1782 | **52** |
| 6 | 288288 | 4290 | 336/5 |
| 7 | 792792 | 9438 | 84 |
| 8 | 1976832 | 19305 | 512/5 |
| 9 | 4550832 | 37180 | 612/5 |
| 10 | 9801792 | 68068 | 144 |

---

## 3. Selection principle does NOT transfer

Paper 2's selection equation is B/N = dim(total space) = dim(S³) = 3,
uniquely hit at n_max = 3 via 3(m+2)(m−1)/10 = 3.

For the higher Hopfs we test all natural targets:

| base | target | equation | integer solution? |
|:---|:---|:---|:---|
| S⁴ | B/N = 4 (dim base) | 2L(L+4)/3 = 4 → L² + 4L − 6 = 0 | L = −2+√10 ≈ 1.16, **no** |
| S⁴ | B/N = 7 (dim total S⁷) | 2L(L+4)/3 = 7 → L² + 4L − 21/2 = 0 | L ≈ 1.79, **no** |
| S⁴ | B/N = 3 (dim S³ fiber) | 2L(L+4)/3 = 3 → L² + 4L − 9/2 = 0 | L ≈ 0.92, **no** |
| S⁴ | B/N = 8/3 (SU(5) ratio) | 2L(L+4)/3 = 8/3 → L² + 4L − 4 = 0 | L ≈ 0.83, **no** |
| S⁸ | B/N = 8 (dim base) | 4L(L+8)/5 = 8 → L² + 8L − 10 = 0 | L ≈ 1.12, **no** |
| S⁸ | B/N = 15 (dim total) | 4L(L+8)/5 = 15 → L² + 8L − 75/4 = 0 | L ≈ 1.90, **no** |
| S⁸ | B/N = 7 (dim fiber S⁷) | 4L(L+8)/5 = 7 → L² + 8L − 35/4 = 0 | L ≈ 0.97, **no** |

**Finding:** The closed-form ratios 2L(L+4)/3 and 4L(L+8)/5 have NO
integer roots at any natural target. The selection principle
B/N = d does not extend to higher Hopf fibrations.

(A Fock-shell-style cutoff (n,l), n ≤ n_max, l ≤ n−1, with S⁴ or S⁸
degeneracies also fails: S⁴ gives B/N = 20/7, 180/27, 880/77, …
and S⁸ gives B/N = 36/5·... — no integer hit at any dim target.)

---

## 4. F-analogs: Dirichlet series of the base degeneracies

Paper 2: F = D_{n²}(s=4) = ζ_R(2) = π²/6. Here the exponent s = 4 = d_max
from Paper 0's packing axiom, and the weight is the S³ Fock degeneracy.

For S⁴ with degeneracy polynomial (l³/3 + 3l²/2 + 13l/6 + 1):

  F_{S⁴}(s) = (1/3) ζ(s−3) + (3/2) ζ(s−2) + (13/6) ζ(s−1) + ζ(s)

Evaluated at natural exponents:

| s | F_{S⁴}(s) exact | num |
|:---|:---|:---|
| 5 | ζ(5) + (13/6)ζ(4) + (3/2)ζ(3) + (1/3)ζ(2) = ζ(5) + (3/2)ζ(3) + (13π⁴/540) + π²/18 | 5.7334 |
| 7 | ζ(7) + (13/6)ζ(6) + (3/2)ζ(5) + (1/3)ζ(4) = ζ(7) + (3/2)ζ(5) + 13π⁶/5670 + π⁴/270 | 5.1288 |

For S⁸ the degeneracy is a degree-7 polynomial in l, and
F_{S⁸}(s) is a sum of 7 shifted Riemann zetas. At s = 9, 11, 13 the
result is a rational combination of ζ(2) through ζ(9).

**None of these values resembles a clean Paper-2-style π^k/rational**:
unlike F_{S³}(s=4) which collapses to ζ(2) because the degeneracy is
the single monomial n², the higher bases have polynomial degeneracies
whose Dirichlet series is a *sum* of Riemann zetas that does not
simplify to a single transcendental.

---

## 5. K and α combinations

K_analog(L, s) = π (B(L) + F(s) − Δ(L)), tested over all reasonable
(L, s) combinations (see sm_f_higher_hopf.json for the full grid).

**Quaternionic S⁴:** K values explode with L because B grows as
~L⁴/3. At L = 1, K ≈ 80; at L = 2, K ≈ 521; at L = 5, K ≈ 18490.
No K value lands near α_W⁻¹ ≈ 29.5, α_s⁻¹ ≈ 8.467, or any
SM coupling constant within 10% relative error.

**Octonionic S⁸:** K values are even larger (B ~ L⁸/315). No hits.

Scan of the full (L, s) grid for hits within 10% of any SM target
(α_em⁻¹, α_W⁻¹, α_s⁻¹, 8/3, sin²θ_W, and their inverses, vs. K, 1/K,
K/π):

**Zero hits.**

The *raw* K combination rule fails for the higher Hopf fibrations.

---

## 6. Isolated numerical coincidences (speculative)

Although the full K formula gives garbage, three B/N ratio hits are
worth noting:

1. **L = 5 on S⁴: B/N = 30.000** (exactly).
   Compared with 1/α_W(M_Z) ≈ 29.5, this is a 1.7% relative error.
   B/N = 2·5·9/3 = 30 is a clean integer Diophantine coincidence of
   the ratio function 2L(L+4)/3, not a spectral identity. It is the
   L = 5 point in a sequence that also gives B/N = 8, 14, 30, 40, 52.

2. **L = 2 on S⁴: B/N = 8**, compared with 1/α_s(M_Z) ≈ 8.467
   (5.5% relative error). This is even weaker — α_s runs strongly with
   scale and there is no canonical value to match.

3. **L = 8 on S⁴: B = 52800**, **B/N = 64**.
   The value 64 = 2⁶ = dim Cl(6), which is the Clifford algebra used
   in Furey-style SM fermion constructions. However, the B value 52800
   has no independent meaning, and N = 825 = 52800/64 is arbitrary.
   This is a numerological coincidence, not a structural link.

None of these survive as structural derivations. They are integer hits
of the polynomial 2L(L+4)/3 at small L values.

---

## 7. Why the construction fails for higher Hopfs

The Paper 2 selection principle depends on three structural features
specific to the complex Hopf:

(a) **Fock projection exists:** The S³ base has the Coulomb Fock
projection (Paper 24's HO rigidity theorem: only the Coulomb potential
gives a linear projection from a one-electron spectrum onto the
round-sphere Laplacian; the quaternionic and octonionic Hopfs have no
such canonical quantum-mechanical origin).

(b) **B/N closed form is quadratic:** For S³,
B(m)/N(m) = 3(m+2)(m−1)/10 is a single quadratic in m, with discriminant
matched to an integer target at m = 3. For S⁴ and S⁸, the ratios
2L(L+4)/3 and 4L(L+8)/5 are also quadratics in L, but their integer
roots at natural targets (4, 7, 3, 8/3, 8, 15, 7) are all irrational.

(c) **Fiber ζ collapses to a single monomial:** D_{n²}(4) = ζ(2) because
the S³ degeneracy is n² (single monomial). For S⁴ and S⁸ the degeneracy
is a polynomial, so F becomes a mixed sum of ζ(2)…ζ(7), losing Paper 2's
clean π²/6 structure.

The Paper 24 asymmetry (π-free for the HO, π-calibrated for the Coulomb)
therefore extends: only the complex Hopf has the specific structural
ingredients needed to assemble α, and the higher division-algebra
analogs do NOT produce the weak or strong couplings by the same
mechanism. The Furey program's success in capturing SM fermion content
from ℂ⊗𝕆 does not imply that the *coupling constants* follow from
extending Paper 2's spectral combination rule to S⁴ or S⁸.

---

## 8. Verdict

**FAIL** on the original pass criterion (at least one of K_4, K_8 within
10% of α_W or α_s). No K value across any (L, s) combination produces a
hit within 10% of any SM coupling.

**PARTIAL** as a set of well-defined spectral invariants:
the computation is rigorous, the closed-form ratios are exact, and the
investigation cleanly rules out the B/N = dim selection principle for
the higher Hopf bases.

**Unexpected numerical hits:**
- B/N = 30 at L=5 on S⁴ (1.7% from α_W⁻¹), purely an integer
  coincidence of 2L(L+4)/3.
- B/N = 8 at L=2 on S⁴ (5.5% from α_s⁻¹), same mechanism.
- B = 52800 at L=8 on S⁴ where B/N = 64 = dim Cl(6) —
  numerological match to the Furey construction.

None of these survive as structural derivations. The combination rule
K = π(B + F − Δ) is a Coulomb/complex-Hopf-specific identity with no
clean extension to the quaternionic or octonionic bases.

**Recommendation:** The SM gauge couplings do not arise from the direct
higher-Hopf extension of Paper 2's construction. If there is a
GeoVac/SM link, it must come from a different mechanism —
possibly Furey-style Clifford-ideal constructions driven by the S³/ℂ
base (which Paper 2 already uses for U(1)_em) combined with a
*separate* group-theoretic argument for SU(2) and SU(3), not from a
spectral-invariant combination on S⁴ and S⁸.
