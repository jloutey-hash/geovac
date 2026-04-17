# Sprint 6 DW: C^(k) Convention Check — Technical Memo

**Sprint:** Sprint 6 Track DW (sub-agent).
**Date:** 2026-04-15.
**Status:** **CLEAN NEGATIVE** — the C^(k) ↔ Y^(k) convention factor does NOT
close the Drake derivation. The obstruction is structural: the bipolar expansion
of C^(2)(r̂₁₂)/r₁₂³ gives radial kernels that are fundamentally different from
Drake's M^k integrals.

---

## 1. The hypothesis tested

Sprint 5 DV found a numerical near-hit: the "equal-β" ansatz for the bipolar
prefactors gave β ≈ 1/√(4π) = 0.2821 (0.64% off). The hypothesis was that
DV used Y^(k) convention where the Breit-Pauli operator uses Racah's C^(k),
and this √(4π/(2k+1)) conversion factor would turn the near-hit into an exact
match.

## 2. What was computed

### 2.1 Brink-Satchler expansion coefficients (exact sympy)

For C^(2)_M(r̂₁₂)/r₁₂³, the bipolar expansion in C^(k) convention gives
channels (k₁, k₂) with k₁ + k₂ = 2:

| Channel | Coefficient |
|:--------|:------------|
| (0, 2)  | √5 ≈ 2.236 |
| (1, 1)  | −√6 ≈ −2.449 |
| (2, 0)  | √5 ≈ 2.236 |

### 2.2 Angular selection (confirmed from DV)

For He (1s)(2p) ³P, Gaunt selection restricts:
- **Direct** (0,1 → 0,1): ONLY (0,2) survives
  - ⟨0‖C⁰‖0⟩ = 1, ⟨1‖C²‖1⟩ = −√30/5
- **Exchange** (0,1 → 1,0): ONLY (1,1) survives
  - ⟨0‖C¹‖1⟩ = −1, ⟨1‖C¹‖0⟩ = +1

No higher-order channels (k₁+k₂ = 4, 6, ...) contribute — the angular
selection is complete at the leading term.

### 2.3 The STRUCTURAL OBSTRUCTION: radial kernel mismatch

The bipolar expansion of the irregular solid harmonic I^K_M = C^K/r^{K+1}
gives radial kernels from the addition theorem:

| Channel | Region I (r₁ < r₂) | Region II (r₁ > r₂) |
|:--------|:---------------------|:---------------------|
| (0, 2)  | 1/r₂³               | r₂²/r₁              |
| (1, 1)  | r₁/r₂²              | r₂/r₁²              |
| (2, 0)  | r₁²/r₂              | 1/r₁³               |

Drake's M^k integrals use kernel r_<^k / r_>^{k+3}:

| M^k | Region I | Region II |
|:----|:---------|:----------|
| M⁰  | 1/r₂³   | 1/r₁³    |
| M¹  | r₁/r₂⁴  | r₂/r₁⁴  |
| M²  | r₁²/r₂⁵ | r₂²/r₁⁵ |

**The power laws are completely different.** For example:
- Bipolar (0,2) Region II: r₂²/r₁ (exponent −1)
- Drake M⁰ Region II: 1/r₁³ (exponent −3)
- Drake M² Region II: r₂²/r₁⁵ (exponent −5)

The bipolar expansion kernel has the form r^{l₁}/r^{l₂+1} with l₁+l₂=2,
while Drake's kernel has the form r^k/r^{k+3} — these are different power laws
with different r_> exponents. **No choice of angular convention (C^(k) vs Y^(k))
can bridge a difference in radial power laws.**

### 2.4 Numerical verification

Direct numerical computation of the bipolar radial integrals at Z=1:

| Integral | Value |
|:---------|------:|
| N^(0,2)_dir (kernel 1/r_>³) | +2.155×10⁻² |
| M⁰_dir (Drake, kernel 1/r_>³) | +2.927×10⁻² |

Even N^(0,2)_dir ≠ M⁰_dir despite apparently matching kernels in Region I.
The discrepancy comes from Region II: bipolar (0,2) Region II has kernel
r₂²/r₁ while M⁰ Region II has kernel 1/r₁³.

## 3. Why the near-hit β ≈ 1/√(4π) was a coincidence

The angular coupling coefficients in the C^(k) convention are:
- Direct (0,2): √5 × (−√30/5) = −√6 ≈ −2.449
- Exchange (1,1): (−√6) × (−√5/3) = √30/3 ≈ 1.826

The ratio 1/√(4π) ≈ 0.282 is close to the geometric mean of 1/√6 and 1/√30,
which naturally arises in the angular coupling. This is a numerical
coincidence from the relative scale of the Wigner 3j/9j symbols, not a
meaningful convention factor.

## 4. Where Drake's M^k kernel r_<^k/r_>^{k+3} comes from

Drake's r_<^k/r_>^{k+3} kernel arises from the **distributional
regularization** of 1/r₁₂³ as a two-body operator. The naive Gegenbauer
expansion of 1/r₁₂³ gives r_<^k/r_>^{k+2} (with a (k+1) prefactor), but
the Breit-Pauli operator involves the **second angular derivative** of
1/r₁₂, which converts the Coulomb kernel r_<^k/r_>^{k+1} to
r_<^k/r_>^{k+3} (two extra powers of 1/r_>).

This is NOT the same as the irregular solid harmonic bipolar expansion,
which uses the addition theorem for C^(K)(r₁₂)/r₁₂^{K+1} without the
distributional regularization.

The Drake kernel is the correct one for computing Breit-Pauli matrix
elements. The bipolar expansion gives a different (incorrect for this
purpose) decomposition.

## 5. Conclusion

**The C^(k) convention hypothesis is NEGATIVE.** The structural
obstruction identified in Sprint 5 DV is NOT a convention issue — it is
a fundamental mismatch between:

1. The **irregular solid harmonic** bipolar expansion (kernel r^{l₁}/r^{l₂+1})
2. The **Breit-Pauli distributional** decomposition (kernel r^k/r^{k+3})

These two decompositions produce different radial power laws. Drake's
(3/50, −2/5) coefficients encode the mapping between these two bases —
they are not simple angular coupling coefficients but rather a
**piecewise radial identity** relating the two kernel types.

A first-principles derivation of (3/50, −2/5) therefore requires
establishing the piecewise radial identity that connects the bipolar
N^{k₁,k₂} integrals to Drake's M^k integrals for hydrogenic (1s)(2p)
wavefunctions. This is a radial-integral algebra problem, not an angular
algebra problem.

## 6. Recommendation

**Close this sub-task with clean negative status.** The hypothesis was
well-defined and definitively refuted.

For future work: the path to closing (3/50, −2/5) is through the
**piecewise radial identity** connecting N^{0,2} and N^{1,1} to M^0, M^1, M^2.
This requires evaluating the symbolic integrals in both decompositions
for hydrogenic (1s)(2p) and algebraically transforming one into the other.
Estimated 2-4 sub-agent hours.

## 7. Files

| File | Purpose |
|:-----|:--------|
| `debug/dv_ck_convention.py` | Full derivation script with numerical verification |
| `debug/dv_ck_convention_memo.md` | This memo |

**NOT MODIFIED:** `geovac/breit_integrals.py`, `tests/test_breit_integrals.py`
(no production changes, no new tests — negative result).
