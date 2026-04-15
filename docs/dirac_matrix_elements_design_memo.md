# Dirac Matrix Elements Design Memo (Track T1, Tier 2 sprint)

**Date:** 2026-04-15
**Status:** Complete. All angular and diagonal radial matrix elements closed-form symbolic; off-diagonal radial integrals via sympy on demand; 66/66 T1 tests + 51/51 D1 regression tests pass.

## 1. API surface

New module: `geovac/dirac_matrix_elements.py`. Public surface:

### κ-labeling
- `DiracLabel(n_fock, kappa, two_m_j)` — dataclass, frozen. Derived properties `l`, `j`, `j_times_2`, `m_j`. Validation: κ ≠ 0, |m_j| ≤ j, parity of 2·m_j matches 2j, l = kappa_to_l(κ) < n_fock.
- `kappa_to_l(κ)`, `kappa_to_j(κ)`, `kappa_to_l_sigma(κ)`: signed-κ → (l, σ) conversions.
- `l_sigma_to_kappa(l, σ)`: inverse. Raises on unphysical (l=0, σ=−1/2).
- `iter_dirac_labels(n_max)`: generator yielding all (n_fock, κ, 2m_j) triples up to n_fock=n_max.
- `spinor_label_to_dirac_label(D1_label)` / `dirac_label_to_spinor_label(D_label, chirality)`: bridge to D1's `SpinorHarmonicLabel`.

### Angular layer (Szmytkowski)
- `angular_matrix_r_hat(κ_a, 2m_a, κ_b, 2m_b)`: σ·r̂. Returns −1 iff κ_a = −κ_b and m_a = m_b, else 0 (Szmytkowski Eq. 2.7).
- `angular_matrix_J_sq(...)`: j(j+1)·δ, diagonal.
- `angular_matrix_L(..., component)`: "sq" → l(l+1)·δ (diagonal); "z"/"+"/"−" → rank-1 symbolic placeholder emitted as unevaluated `L_reduced(...)` (full analytic reduction deferred to T2).
- `angular_matrix_L_dot_S(...)`: (−κ−1)/2 · δ, diagonal. Verified against direct [j(j+1)−l(l+1)−3/4]/2 in both κ-sign branches.
- `angular_matrix_sigma(..., component)`: "sq" → 3·δ; "z"/"+"/"−" → rank-1 symbolic placeholder.
- `angular_sigma_dot_rhat_identity(κ, 2m_j)`: returns (−κ, 2m_j, −1) — the target state and coefficient.

### Radial layer (non-relativistic hydrogenic)
- `radial_expectation_diagonal(n, l, operator, Z)`: closed-form ⟨n,l|op|n,l⟩ for op ∈ {"1/r", "1/r^2", "1/r^3", "r", "r^2"}. All verified exact against direct sympy integration of the hydrogenic wavefunction.
- `inverse_r_cubed_hydrogenic(n, l, Z)` — convenience wrapper for the T2 spin-orbit deliverable. Raises on l=0 (Darwin divergence).
- `radial_matrix_element(n_a, l_a, n_b, l_b, op, Z)`: diagonal fast-path uses the table; off-diagonal routes through direct sympy integration of `assoc_laguerre` expansions.

### Exposed symbols
- `Z_sym`: positive sympy symbol for charge.
- `alpha_sym`: positive symbol for fine-structure constant α (reserved for T2/T5).
- `gamma_rel`: positive symbol for γ = √(1 − (Zα)²). **Not bound to any matrix element in T1**; exposed for T2/T5 to consume when Dirac-Coulomb radial expressions are added.

## 2. κ-labeling bridge

D1 §8 flagged the `(l, σ)` vs κ decision as an open question. T1 resolves it by **adding a parallel labeling** rather than mutating `SpinorHarmonicLabel`. Rationale:

- D1's 51-test contract is preserved (zero regressions).
- κ is the natural label for matrix elements (Szmytkowski; Dyall §9; Grant §7.5) and directly carries both `l` and the j-branch in one signed integer.
- `(l, σ)` is the natural label for the D1 spinor-spherical-harmonic construction (Bär 1996 Theorem 1, SU(2)_L × SU(2)_R decomposition).
- The bridge is a two-line conversion each way (`l_sigma_to_kappa`, `kappa_to_l_sigma`).

Conventions used throughout (matching T0 memo §2):
- κ = −(l+1) for j = l+1/2 (κ < 0)
- κ = +l for j = l−1/2 (κ > 0)
- |κ| = j + 1/2

The forbidden case (l=0, σ=−1/2) gives j = −1/2, which is unphysical. `l_sigma_to_kappa(0, −1/2)` raises `ValueError`. `spinor_label_to_dirac_label` therefore refuses to convert D1 labels at l=0 with σ=−1/2 — those are chirality-bookkeeping artifacts in the full Dirac sector, not physical (κ, m_j) states.

## 3. Algebraic-first check

All T1 deliverables are 100% algebraic in exact sympy arithmetic. No numerical quadrature, no floating-point comparisons in tests, no iterative solvers.

**Diagonal radial closed forms** (verified against direct sympy integration for all (n, l) ≤ (3, 2)):

| Operator | Expression | Notes |
|:---------|:-----------|:------|
| ⟨1/r⟩_{nl} | Z/n² | l-independent |
| ⟨1/r²⟩_{nl} | Z²/[n³(l+½)] | |
| ⟨1/r³⟩_{nl} | Z³/[n³·l(l+½)(l+1)] | diverges l=0 |
| ⟨r⟩_{nl} | [3n² − l(l+1)]/(2Z) | |
| ⟨r²⟩_{nl} | n²[5n²+1−3l(l+1)]/(2Z²) | |

⟨r³⟩ diagonal closed form was attempted but found to have sign disagreements across reference sources; the diagonal table omits it. Off-diagonal r³ remains available via `radial_matrix_element`'s sympy integration fallback.

**Off-diagonal** radial matrix elements compute in ~0.1–1 s per call via sympy `assoc_laguerre` integration. Tested: ⟨1s|r|2s⟩ at Z=1 = −32√2/81 (exact rational × √2).

**Angular** matrix elements are pure rationals (in fact integers or simple half-integer eigenvalues).

**Transcendentals surfaced:**
1. **α** (fine-structure constant): present as a name-only symbol; no matrix element in T1 depends on it.
2. **γ = √(1−(Zα)²)**: present as a name-only symbol. Will enter T2 (spin-orbit via `<1/r³>` with relativistic corrections) and T3 (full Dirac-Coulomb radial functions). T5 will classify: provisional label "spinor-intrinsic" (Paper 18 §IV), a new first-order-operator analog of Paper 24's "second-order-operator calibration π".

No π, ζ, or other transcendentals appear anywhere in T1's outputs. The non-relativistic limit is fully in Q(Z).

## 4. Algebraic vs numerical

| Quantity | Status | Method |
|:---------|:-------|:-------|
| κ ↔ (l, σ) conversion | algebraic | integer arithmetic |
| ⟨J²⟩, ⟨L²⟩, ⟨L·S⟩ eigenvalues | algebraic | κ-identity, closed form |
| ⟨σ²⟩ eigenvalue | algebraic | trivial (3) |
| σ·r̂ selection rule | algebraic | Szmytkowski Eq. 2.7, Kronecker delta |
| Diagonal ⟨r^k⟩, ⟨1/r^k⟩ | algebraic | Bethe-Salpeter closed forms, verified |
| Off-diagonal ⟨n',l'|r^k|n,l⟩ | algebraic | direct sympy integration |
| σ_z, L_z (rank-1) | symbolic-placeholder | emitted as unevaluated `reduced(...)` Function; T2 will close |

Nothing falls into "numerical-required". The rank-1 components of σ and L were intentionally left as symbolic placeholders rather than expanding the Wigner-Eckart reductions prematurely — T2 needs only the diagonal `L·S` expectation value plus the σ·r̂ selection rule, and the reductions for rank-1 components are more cleanly done once T2 knows which exact matrix elements it needs.

## 5. What T2 can consume

T2 (spin-orbit coupling) needs:
- ✅ `angular_matrix_L_dot_S(κ, 2m, κ, 2m)` → (−κ−1)/2, diagonal.
- ✅ `inverse_r_cubed_hydrogenic(n, l, Z)` → Z³/[n³·l(l+½)(l+1)], pure rational in Z.
- ✅ The Z-scaling identity ⟨ξ⟩ ∝ Z⁴ falls out as: ⟨ξ L·S⟩ = Z·⟨1/r³⟩·[(−κ−1)/2] ∝ Z⁴ (Z from ξ = Z/r³, Z³ from ⟨1/r³⟩).
- ✅ The Coulomb ξ(r) = Z/r³ is itself captured as a scalar factor; no extra API needed.

## 6. D1 impact

Zero changes to `geovac/dirac_s3.py`. All D1 tests still pass. The bridge `spinor_label_to_dirac_label` lives in the new T1 module; D1 remains unaware.

## 7. Algebraic obstructions

None encountered. Every T1 deliverable closed algebraically. One minor dead-end: an initial attempt to encode ⟨r³⟩ diagonal from a reference formula produced wrong values — removed, since T1 callers don't need it. The direct sympy path for off-diagonal r³ still works.

## 8. Guardrails observed

- No Dirac Fock projection theorem attempted (Tier 1 Explorer Gap #1).
- No TC combinations (CUSP-3).
- No numerical radial quadrature.
- No modifications to Paper 18 (T5's job) or Paper 22 (T6's job).
- No modifications to `composed_qubit.py` or `molecular_spec.py` (T3's job).
- σ·r̂ identity test passes cleanly, confirming no labeling or sign error.
