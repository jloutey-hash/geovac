# Internal Multi-Focal Architecture for Atomic Excited States — Design Memo

**Sprint:** post-Track-4 He-oscillator named follow-on
**Date:** 2026-05-09
**Phase:** A (Architecture design)
**Status:** Pre-implementation memo, design only.

---

## 1. Motivation: the Track-4 diagnosis

Track 4 (CLAUDE.md §1.8 §V.C.5 He oscillator strength autopsy, today)
found that the framework reproduces

    f(2¹P → 1¹S, He) ≈ 0.444  vs Drake 0.276    (+61% residual, plateaued at n_max = 7)

with a clean smoking gun: the 2¹P excited state is **variationally
violated** by 0.08 Ha (E_GeoVac = −2.204 vs Drake −2.124), the
signature of the small-Z graph-validity-boundary mechanism (Z_c ≈ 1.84,
He at Z = 2 is just above) acting on excited Rydberg-class states.

The Lyman α sanity check reproduced 24576/59049 to machine precision,
so the angular projections (Wigner 3j §III.8, vector-photon §III.11)
are exact. **The failure lives in the radial / CI sector.** The single
common Sturmian exponent k = Z = 2 is qualitatively wrong: it cannot
represent both a doubly-screened compact 1s electron and a singly-screened
diffuse 2p electron simultaneously.

This is the structural analog of the cross-register multi-focal wall
that `geovac/cross_register_vne.py` closed for chemistry observables:

| Wall | Mechanism | Closure module |
|:-----|:---------|:---------------|
| W1a (cross-register V_eN) | electron + nucleus on distinct λ | `cross_register_vne.py` (~990 lines) |
| W1c (cross-center screened V_ne) | core orbital ⊗ valence orbital on distinct centers + λ | `cross_center_screened_vne.py` |
| **internal multi-focal** (this sprint) | **two electrons on the same nucleus at distinct λ_orb per orbital** | **`internal_multifocal.py` (this sprint)** |

The architectural pattern is the same: each orbital gets its own focal
length, and the cross-exponent matrix elements (overlap, kinetic, dipole,
two-electron Coulomb) need to be computed against this enlarged basis.

---

## 2. The multi-focal basis

We label each orbital by **four** quantum numbers (n, l, m, λ) instead of
the standard (n, l, m). Concretely:

```
Orbital_p = (n_p, l_p, m_p, λ_p)
```

with the radial wavefunction

```
R_p(r) = N_p (2 λ_p r)^{l_p} e^{-λ_p r} L_{n_p − l_p − 1}^{2 l_p + 1}(2 λ_p r)
```

normalized to ∫|R_p|² r² dr = 1. This is **standard L2 hydrogenic normalization
at independent exponent λ_p** — exactly what `_hydrogenic_poly_coeffs_lam`
in `shibuya_wulfman.py` already produces. (NOT the Coulomb-Sturmian
normalization where ∫R_p R_q r dr = δ_pq; we want orbitals that look
like physical bound states for each focal length, not Coulomb-Sturmians.)

For the He case, the natural per-orbital choice (Slater rules / variational):

| Orbital | n | l | λ | Why |
|:--------|:-:|:-:|:-:|:----|
| 1s | 1 | 0 | 27/16 = 1.6875 | He variational ground state, also Slater rule (Z − 5/16) |
| 2s | 2 | 0 | (Z − 0.85)/2 = 0.575 | Slater rule (n=2 valence, 1s shields by 0.85) |
| 2p | 2 | 1 | (Z − 0.85)/2 = 0.575 | Slater rule, same as 2s |
| 3s, 3p, 3d, ... | 3, ... | 0, 1, 2 | depends on screening | screened valence convention |

The choice is **physical, not fitted**: Slater rules are textbook
empirical screening predictions; variational 1s minimizes the
single-electron Hartree-Fock-like energy. This sprint uses Slater
exponents; alternative conventions are documented and selectable.

---

## 3. The new matrix elements

The multi-focal CI builds three classes of matrix elements that the
single-exponent path doesn't have to deal with:

### 3.1 Cross-exponent overlap S_pq

```
S_pq = ⟨p|q⟩ = δ_{l_p l_q} δ_{m_p m_q} ∫₀^∞ R_p(r) R_q(r) r² dr
```

For l_p = l_q (same angular character but possibly different n and λ),
this is **non-zero** in general and is the key new content. With both
R's of the form e^{−λr} × poly(r) at exponents λ_p and λ_q, the integral
is a finite sum of terms ∫₀^∞ r^k e^{−(λ_p + λ_q) r} dr = k! / (λ_p+λ_q)^{k+1}.
Closed form, machine precision.

For (n_p = n_q, l_p = l_q, λ_p = λ_q), S_pq = 1 by normalization.
For (l_p ≠ l_q), S_pq = 0 by angular orthogonality.

### 3.2 Cross-exponent dipole ⟨p|r|q⟩

```
⟨p|z|q⟩ = c_1(l_p, m_p, l_q, m_q) · D_pq
D_pq = ∫₀^∞ R_p(r) R_q(r) r · r² dr = ∫₀^∞ R_p R_q r³ dr
```

Selection: c_1 nonzero requires l_q = l_p ± 1 and m_p = m_q.

The radial integral D_pq with R's at distinct λ_p, λ_q closes in the
same finite-sum form as S_pq: each polynomial term r^k contributes
(k+3)! / (λ_p+λ_q)^{k+4}.

### 3.3 Cross-exponent two-electron Slater integral R^k_{pqrs}

```
R^k_{p₁q₁; p₂q₂} = ∫∫ R_{p₁}(r₁) R_{q₁}(r₁) (r_<^k/r_>^{k+1}) R_{p₂}(r₂) R_{q₂}(r₂) r₁² r₂² dr₁ dr₂
```

Here the four orbitals can each have their own λ. Crucially the **angular
factors c_k are unchanged** — Gaunt selection rules are angular content,
independent of the radial exponents. **Multipole termination at L_max =
min(l_p+l_q, l_r+l_s) is preserved** (Paper 22 angular sparsity theorem;
the proof depends only on Wigner 3j, which doesn't see λ).

The radial pair-densities ρ_{p₁q₁}(r) = R_{p₁}(r) R_{q₁}(r) r² and similarly
ρ_{p₂q₂} have decay rates λ_{p₁} + λ_{q₁} and λ_{p₂} + λ_{q₂} respectively
(can differ). The double radial integral splits at r₁ = r₂ and reduces to
one-dimensional incomplete-gamma integrals via the same machinery
`shibuya_wulfman._split_integral_analytical` already implements at fixed
R_AB. Here the integration is over r₂ (the inner) and r₁ (the outer)
both from 0 to ∞. The architecture is exactly the cross-register
double-radial integral with R_AB → r₁ floating and one Gauss-Laguerre
quadrature on r₁.

**Implementation:** delegate to a 1D Gauss-Laguerre quadrature on r₁
with the natural weight e^{−(λ_{p₁} + λ_{q₁}) r₁}, evaluating the inner
r₂ integral analytically at each node via `_split_integral_analytical`.
This is bit-identical to the cross-register pattern
(`_cross_register_double_integral_via_gauss`) with both registers being
the same electronic register but at distinct λ's.

**Safety:** at λ_p = λ_q = ... = λ (matched-exponent limit), this must
reduce bit-identical (modulo Gauss-Laguerre quadrature noise) to the
single-exponent path that `casimir_ci.two_electron_integral` uses
(which calls `compute_rk_float` with k_orb scaling). Regression test
non-negotiable.

### 3.4 One-body kinetic + nuclear attraction

The Hamiltonian on the multi-focal basis is

    H = T + V_{Ne} + V_{ee}
    T_{pq}    = ⟨p| (−½ ∇²) |q⟩
    V_{Ne,pq} = ⟨p| (−Z/r) |q⟩
    V_{ee}    = sum over Slater integrals

For l_p = l_q only (kinetic and Coulomb don't change l):
- T_pq is closed-form via -½ ∇² acting on hydrogenic radial functions.
  For orbitals at λ, T |R_λ⟩ = (-½ λ² + λ²/n²) |R_λ⟩ + cross terms ...
  the cleanest approach is to evaluate ⟨p| -½ ∇² |q⟩ as an integral of
  R_p (-½ d²/dr² + l(l+1)/(2r²)) R_q with r² Jacobian.
- ⟨p | 1/r | q⟩ is the same integral form as the dipole with k=−1
  instead of k=+1.

For the **first-pass implementation**, we use a clean approach: the
single-exponent diagonal hydrogenic eigenvalue formula −Z²/(2n²) does
NOT apply when λ ≠ Z/n. Instead we evaluate T and V_Ne explicitly
from the radial integrals.

For ⟨p| T |q⟩ on hydrogenic Sturmian-style basis at λ:
  −½ d²/dr² R_λ(r) acting on R_{n,l}^λ(r) ... we use the identity
  
    (−½ ∇² − λ Z/(n_λ_eff r) + λ²/(2 n_λ_eff²)) |n,l,λ_natural⟩ = 0
  
  only if λ = Z/n. For arbitrary (n, λ), R_{n,l}^λ is NOT an eigenfunction
  of the natural Hamiltonian; it's a basis function.
  
  Numerically simpler: evaluate the matrix elements directly via the
  explicit Laguerre form. This is closed-form (sum of incomplete gammas).

**For Phase B**, we keep the kinetic + V_Ne evaluation simple and
delegate to direct radial integration (sympy or finite-sum closed form).
The CI sector inherits all the angular machinery from `casimir_ci.py`
unchanged.

---

## 4. Module structure

Single new module `geovac/internal_multifocal.py` (~600-800 lines):

```
geovac/internal_multifocal.py
├── Class MultifocalSpec
│     orbitals: list of (n, l, m, lambda) tuples
│     Z_nuc: nuclear charge
│     label: human-readable
│
├── Cross-exponent radial primitives
│   ├── overlap_radial(n_p, l_p, lam_p, n_q, l_q, lam_q) -> float
│   ├── matrix_element_rk(n_p, l_p, lam_p, n_q, l_q, lam_q, k) -> float
│   │     # ⟨R_p | r^k | R_q⟩, closed-form sum of factorial / (lam_p+lam_q)^...
│   ├── kinetic_radial(n_p, l_p, lam_p, n_q, l_q, lam_q) -> float
│   ├── slater_rk_multifocal(n1, l1, lam1, n3, l3, lam3,
│   │                        n2, l2, lam2, n4, l4, lam4, k) -> float
│
├── Single-particle assemblies
│   ├── h1_multifocal(spec) -> (h1_matrix, orbitals_list)
│   │     # T + V_Ne, full matrix in (orbital_p, orbital_q) labels
│   ├── overlap_multifocal(spec) -> overlap_matrix
│   ├── dipole_z_multifocal(spec) -> dipole_z_matrix
│
├── CI on multi-focal basis (mirrors casimir_ci.build_singlet_LM_subblock)
│   ├── build_singlet_LM_subblock_multifocal(spec, L, M_L)
│   │     -> (H, configs, orbitals)
│   │
│   │   IMPORTANT: must handle non-orthogonal basis (S ≠ I).
│   │   Two paths:
│   │     (a) Generalized eigenvalue H c = E S c (preferred, exact).
│   │     (b) Symmetric-orthogonalization S^{-1/2} H S^{-1/2} = H'
│   │         then standard eigh on H' (cheaper).
│   │   First-pass: use scipy.linalg.eigh with the b kwarg
│   │   (generalized symmetric eigenproblem).
│
└── Driver
    └── compute_oscillator_strength_multifocal(spec_init, spec_final, ...)
```

Existing `casimir_ci.py` is **extended** in one minor way:

* Add a `subblock_multifocal=True` flag (or just import the new module
  and use its functions). Phase B opts for "import the new module"
  rather than touching `casimir_ci.py` core to keep the existing
  bit-identical regression unaffected.

---

## 5. Verification plan (Phase B tests, target ≥ 10)

Each test is named `test_{topic}_{specifier}`:

1. **test_overlap_self_normalized**: at λ_p = λ_q = Z/n, S_pp = 1 (basis normalization).
2. **test_overlap_orthogonal_n_at_matched_lambda**: at λ = Z/n_p = Z/n_q, S_{1s,2s} = 0 (standard orthogonality at matched λ).
3. **test_overlap_nonzero_at_mismatched_lambda**: at λ_1s ≠ λ_2s, S_{1s,2s} ≠ 0 (cross-exponent overlap is the new content).
4. **test_dipole_lyman_alpha_match**: with λ_1s = λ_2p = 1, |⟨1s|z|2p_{m=0}⟩| reproduces the standard hydrogenic value 128√2/243 to 1e-12.
5. **test_dipole_at_mismatched_lambda_finite**: at λ_1s = 1.69, λ_2p = 0.575, dipole is well-defined and finite.
6. **test_slater_matched_lambda_regression**: at all four λ's equal to 1.0, slater_rk_multifocal reproduces compute_rk_float bit-identically (machine-precision, modulo GL quadrature noise <1e-10).
7. **test_slater_F0_1s1s_at_lambda**: at λ = 1.0, F^0(1s,1s) reproduces the textbook 5/8.
8. **test_slater_gaunt_selection_preserved**: at distinct λ's, F^k(1s,1s,2p,2p) is zero for forbidden k (via angular factor); cross-exponent ERIs preserve Gaunt.
9. **test_h1_hermitian**: T + V_Ne matrix is Hermitian (symmetric, real).
10. **test_overlap_hermitian**: S matrix is Hermitian (symmetric, real, positive-definite).
11. **test_he_singleton_sanity**: at single-orbital basis (just 1s at λ=27/16), eigenvalue of (T + V_Ne)/S = -27²/(2·16²) - Z·(27/16)·… — verifies hydrogenic energy formula at variational λ.
12. **test_he_2config_correlation**: 1s² + 1s2s configuration, λ_1s=27/16 λ_2s=0.575, gives ground state below pure 1s² (correlation energy nonzero, sign correct).
13. **test_radial_integral_machine_precision**: hand-computed cross-exponent overlap S_{1s(λ=1), 1s(λ=2)} = (2·sqrt(2))³/(1+2)³ = 16√2/27 reproduced to 1e-14.
14. **test_dipole_z_selection_rule**: ⟨2s|z|3s⟩ = 0 (l-selection); ⟨2p_0|z|2p_0⟩ = 0 (parity selection).

Plus 1 slow Phase C sanity:

15. **test_he_oscillator_phase_c_basic** (slow): at the simplest multifocal basis (1s @ 1.69, 2s @ 0.575, 2p @ 0.575), compute oscillator strength and verify it's in the [0.1, 0.6] range — physical bounds, not yet a precision test.

---

## 6. Honest limitations of the Phase B implementation

* **Non-Coulomb-Sturmian basis**: this is a **standard hydrogenic L2 basis at
  arbitrary per-orbital λ**, NOT the Coulomb-Sturmian basis (different
  normalization). The Coulomb-Sturmian closure property identified in
  Sprint Calc-P is intrinsic to the Coulomb-Sturmian normalization and
  is NOT what this architecture targets. We are not pursuing Sturmian
  closure with this architecture — we are addressing the orthogonal
  question of whether using physical screening exponents per orbital
  closes the multi-electron observable.

* **No graph κ adjacency at the multi-focal level**: the graph one-body
  Hamiltonian h1 = h1_diag + κ·(−A) of `casimir_ci._build_graph_h1` is
  built on the (n, l, m) graph, NOT on the (n, l, m, λ) extended labels.
  In Phase B we use the **physical kinetic + V_Ne** matrix (no κ
  adjacency) on the multi-focal basis. This means we are testing whether
  the physical (non-κ-corrected) eigenvalue problem on the multi-focal
  basis closes the residual. The κ adjacency is the framework's
  graph-native ingredient; without it we lose the "graph-native"
  property and become a standard variational hydrogenic CI.
  
  This is the correct experiment: Track 4 named κ as the over-binding
  mechanism; replacing the basis is one fix, removing κ is another. The
  multi-focal basis with physical h1 tests whether the basis fix alone
  closes the residual. If yes: we know κ + single-focal was the problem
  and one fix at a time suffices. If no: the basis fix needs to be
  combined with a κ-fix or the diagnosis needs sharpening.

* **Gauss-Laguerre quadrature for the outer Slater integral**: the
  cross-exponent Slater R^k uses a GL quadrature on r₁ at n_quad ~ 200
  for machine-precision (same as `cross_register_vne`). At λ_p = λ_q
  the inner r₂ integral closes via incomplete gamma; the outer r₁
  integral is then a polynomial × exponential against incomplete gamma.
  GL quadrature for this combination converges at n_quad ~ 100-200.

* **Generalized eigenproblem H c = E S c**: scipy.linalg.eigh with the
  `b=S_matrix` kwarg handles this directly. Numerical care: S must be
  positive-definite. At very small λ separation this can become
  ill-conditioned; we catch this with a condition-number check.

* **First-pass scope**: just enough to test the Phase C oscillator-strength
  question. NOT a full re-architecture of the He CI infrastructure.

---

## 7. Phase C test: what's a "win"

Drake reference: f = 0.27616.
Track 4 single-focal residual: +61% (f = 0.444).

| Phase C f outcome | Verdict | Implication |
|:------------------|:--------|:------------|
| < 0.30 (≤ 9% residual) | **WIN** | Diagnosis validated; multi-focal architecture extends to multi-electron precision. |
| 0.30–0.35 (9-26% residual) | **PARTIAL** | Multi-focal helps significantly but other physics (correlation, κ adjustment) needed. |
| 0.35–0.40 (26-45% residual) | **WEAK** | Multi-focal helps modestly; primary residual is elsewhere. |
| > 0.40 (≥ 45% residual) | **NEGATIVE** | Diagnosis broken; the basis isn't the dominant issue. Sharpen scope. |

The honest pre-experiment expectation: **PARTIAL is most likely.**
Reasoning:
- Multi-focal addresses the basis-architecture inadequacy clearly named
  in Track 4 §6 (single λ can't represent two length scales).
- It does NOT fix the κ-induced inter-shell over-binding directly.
  The graph-validity-boundary at small Z is a separate effect (the κ
  perturbation grows as 1/(8Z²); at He it's ~3% of the Hartree shell
  spacing, enough to over-bind diffuse excited states).
- A clean-win outcome would require both effects to be predominantly
  basis-driven, not κ-driven. If κ contributes ~half the residual,
  even a perfect basis would leave ~30% residual.

A PARTIAL outcome would still be a major positive: it confirms the
diagnosis class (basis architecture matters, not just angular machinery
or vector-photon promotion), and it cleanly separates the two
contributions for Phase D scoping (κ-on-multifocal as the next
follow-on).

---

## 8. Files this Phase A produces

* `debug/internal_multifocal_design_memo.md` (this file)

Phase B produces:
* `geovac/internal_multifocal.py` (~600-800 lines)
* `tests/test_internal_multifocal.py` (≥ 10 tests, target 14)

Phase C produces:
* `debug/calc_track_he_oscillator_v2.py` (driver)
* `debug/data/he_oscillator_v2.json` (results)
* `debug/internal_multifocal_implementation_memo.md` (results + scoping)

No production GeoVac code is modified outside `geovac/internal_multifocal.py`
in Phase B. The existing single-focal CI path remains bit-identical.
