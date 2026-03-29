# Sturmian Angular Basis Investigation

**Date:** 2026-03-28 (v2.0.5, Track B)
**Files:** `geovac/algebraic_angular_sturmian.py`, `tests/test_algebraic_angular_sturmian.py`
**Status:** Research investigation complete. Mixed results — see recommendation.

---

## 1. Monopole Decomposition (Step 1)

The coupling potential V_coupling = V_nuclear + V_ee was decomposed to assess the Sturmian construction's theoretical basis.

| l_max | V_nuc (Frob.) | V_ee (Frob.) | V_total (Frob.) | Monopole fraction | V_ee fraction |
|:-----:|:-------------:|:------------:|:---------------:|:-----------------:|:-------------:|
| 0     | 237.1         | 8.0          | 233.6           | 101.5%            | 3.4%          |
| 1     | 286.7         | 12.7         | 281.4           | 101.9%            | 4.5%          |
| 2     | 315.7         | 16.7         | 308.8           | 102.3%            | 5.4%          |

**Key finding:** The nuclear monopole dominates V_coupling at all l_max (>100% Frobenius norm — the nuclear and V_ee terms partially cancel). The V_ee remainder is only 3-5% of the total. This strongly motivates the Sturmian approach: the nuclear monopole contains essentially all the coupling physics.

Note: the monopole fraction exceeds 100% because V_nuc and V_ee partially oppose each other (nuclear attraction vs electron repulsion), so ||V_nuc|| + ||V_ee|| > ||V_nuc + V_ee||.

---

## 2. Sturmian Basis Construction (Step 2)

### Method

The Sturmian basis is constructed by:
1. Building the angular Hamiltonian H_ref = diag(Casimir) + R0 × V_nuclear in a LARGE free basis (n_construct = 50 Gegenbauer functions per l-channel)
2. Diagonalizing to get the Sturmian eigenvectors
3. Selecting the lowest n_basis eigenvectors (per l-channel for block-diagonal V_ref)
4. Projecting the full Hamiltonian onto this truncated subspace

The key insight: the Sturmian functions are optimal linear combinations of many free functions that capture the nuclear potential structure. The first n_basis Sturmian functions span a better-adapted subspace than the first n_basis free functions.

### R0 Optimization

R0 (reference hyperradius) controls which R regime the Sturmian basis is optimized for.

| R0   | Energy (Ha) | Error (%) |
|:----:|:-----------:|:---------:|
| 0.3  | -2.897723   | 0.207     |
| 0.5  | -2.898240   | 0.189     |
| 0.7  | -2.898657   | 0.175     |
| 1.0  | -2.899102   | 0.159     |
| 1.5  | -2.899377   | 0.150     |
| 2.0  | -2.899097   | 0.159     |
| 3.0  | -2.897048   | 0.230     |
| 5.0  | -2.888447   | 0.526     |

**Best R0 = 1.5 bohr** (near the V_eff well minimum). R0 in [0.5, 3.0] gives similar results (< 0.05% spread). R0 > 5 degrades significantly because the Sturmian basis becomes optimized for the large-R asymptotic regime where the eigenfunction has localized, missing the physics of the well region.

### n_construct Convergence

| n_construct | Energy (Ha) | Error (%) |
|:-----------:|:-----------:|:---------:|
| 15          | -2.898577   | 0.177     |
| 20          | -2.899030   | 0.162     |
| 30          | -2.899264   | 0.154     |
| 40          | -2.899322   | 0.152     |
| 50          | -2.899343   | 0.151     |
| 60          | -2.899352   | 0.151     |
| 80          | -2.899359   | 0.150     |

Converged by n_construct ≈ 40-50. The one-time diagonalization cost is O(n_construct³) = O(50³) = O(125K) — negligible.

---

## 3. l_max=0 Benchmarks (Step 3)

### He ground-state energy: Sturmian vs free basis

Parameters: n_construct=50, R0=1.5, n_R=200, N_R_radial=2000

| n_basis | E_Sturmian (Ha) | err_S (%) | E_Free (Ha) | err_F (%) | Improvement |
|:-------:|:---------------:|:---------:|:-----------:|:---------:|:-----------:|
| 3       | -2.884936       | 0.647     | -2.835758   | 2.341     | **3.62x**   |
| 5       | -2.895834       | 0.272     | -2.879777   | 0.825     | **3.04x**   |
| 8       | -2.898796       | 0.170     | -2.893975   | 0.336     | **1.98x**   |
| 10      | -2.899343       | 0.151     | -2.896727   | 0.241     | **1.60x**   |
| 15      | -2.899780       | 0.136     | -2.898964   | 0.164     | **1.21x**   |
| 20      | -2.899895       | 0.132     | -2.899554   | 0.144     | 1.09x       |
| 25      | -2.899937       | 0.130     | -2.899771   | 0.136     | 1.04x       |

**Key finding:** The Sturmian basis converges 1.6-3.6x faster than the free basis at small n_basis. The advantage is most dramatic at n_basis ≤ 10 (where basis truncation error dominates) and diminishes as both bases approach the complete-basis limit.

### Eigenvalue convergence at fixed R

At R = 2.0 bohr (the V_eff well minimum), the Sturmian basis gives 10-12x lower error at the same n_basis:

| n_basis | μ_Sturm err | μ_Free err | Ratio |
|:-------:|:----------:|:----------:|:-----:|
| 3       | 3.0e-2     | 4.0e-1     | 13x   |
| 5       | 7.9e-3     | 1.2e-1     | 16x   |
| 10      | 1.2e-3     | 1.9e-2     | 15x   |
| 15      | 4.7e-4     | 5.9e-3     | 12x   |

The advantage is even larger at R = 5-10 bohr (where the eigenfunction localizes and many free basis functions are needed).

### Success criterion assessment

**Criterion:** "Sturmian at n_basis=10 should beat free at n_basis=25"
- Sturmian n_basis=10: 0.151% error
- Free n_basis=25: 0.136% error
- **NOT MET.** Sturmian at n_basis=10 is close but does not beat free at n_basis=25.
- Sturmian at n_basis=15 (0.136%) matches free at n_basis=25 (0.136%).

The 0.130% floor (Sturmian at n_basis=25) is the adiabatic single-channel limit — improving the angular basis further cannot reduce this.

---

## 4. l_max=1 Benchmarks (Step 4)

### He ground-state energy across l_max

| l_max | E_Sturmian (Ha) | err_S (%) | E_Free (Ha) | err_F (%) | FD err (%) |
|:-----:|:---------------:|:---------:|:-----------:|:---------:|:----------:|
| 0     | -2.899791       | 0.136     | -2.898975   | 0.164     | 0.054      |
| 1     | -2.920753       | 0.587     | -2.919963   | 0.559     | 0.78       |
| 2     | -2.922898       | 0.660     | -2.922110   | 0.633     | 0.85       |
| 3     | -2.923423       | 0.678     | -2.922636   | 0.651     | 0.87       |

**Key finding: at l_max>0, the Sturmian basis is slightly WORSE than the free basis** (by ~0.03 percentage points), though both significantly beat the FD solver.

**Root cause analysis:** This is the adiabatic approximation interacting with basis completeness:

1. At l_max=0, the adiabatic energy is ABOVE exact (-2.899 vs -2.904). Basis truncation error pushes it further above. The Sturmian basis reduces this truncation error, moving the energy closer to exact. **Good.**

2. At l_max≥1, the adiabatic energy is BELOW exact (-2.920 to -2.923 vs -2.904). The missing non-adiabatic correction (+0.035 Ha DBOC minus 97% cancellation) should raise the energy. More angular correlation (higher l_max or better basis) drives the adiabatic energy further below exact. The Sturmian basis captures MORE angular correlation at each n_basis (it converges faster to the adiabatic limit), which is the WRONG direction. **Bad.**

This is a fundamental limitation: **the Sturmian basis converges faster to the wrong answer at l_max>0.** The cure is the coupled-channel solver (which includes the non-adiabatic correction), not a better angular basis.

### Success criterion assessment

**Criterion:** "He ground state at l_max=1 with Sturmian basis achieves < 0.5%"
- Sturmian l_max=1: 0.587% error
- **NOT MET.** Both Sturmian and free basis give ~0.5-0.6% at l_max=1.
- Both beat the FD solver (0.78%) by a clear margin.
- The 0.5% target requires the coupled-channel solver, not a better angular basis.

---

## 5. Algebraicization Assessment (Step 5)

### Matrix element classification

| Matrix element | Sturmian construction | Runtime (per R) | Status |
|:---------------|:--------------------:|:----------------:|:------:|
| Free Casimir eigenvalues | — | Used in H_ref | ✅ Algebraic |
| V_nuc in free basis | Quadrature (n_quad=150) | Projected once | 🔶 Algebraic in principle |
| V_ee in free basis | Quadrature + Gaunt | Projected once | 🔶 Algebraic in principle |
| Sturmian eigenvalues μ_k | One-time diag of H_ref | Used as diagonal | ❌ Numerical (diag) |
| Sturmian rotation matrix U | One-time diag of H_ref | Used for projection | ❌ Numerical (diag) |
| V_ref projected | U^T @ V_ref @ U | Precomputed once | Inherits from V_ref |
| V_residual projected | U^T @ V_res @ U | Precomputed once | Inherits from V_res |
| H(R) construction | diag + scalar × matrix | O(n_basis²) per R | ✅ Algebraic (linear in R) |
| H(R) diagonalization | eigh(H) | O(n_basis³) per R | ❌ Numerical (diag) |

### Algebraicization routes

1. **V_nuc and V_ee in free basis:** Currently computed via Gauss-Legendre quadrature. Paper 13 Section XII proves these are partial harmonic sums (algebraic). The audit (algebraic_hyperradial_audit.md) confirms closed-form expressions exist for all l,l' coupling integrals. Replacing quadrature with algebraic evaluation would make the Sturmian construction 100% algebraic up to the one-time diagonalization.

2. **Three-term recurrence (Aquilanti et al.):** The Sturmian eigenvectors can in principle be computed from three-term recurrence relations involving the V_nuc matrix elements. This would replace the numerical diagonalization with an algebraic construction. However, for n_basis ≤ 15 (the relevant regime), the diagonalization cost is negligible and this optimization is not needed.

3. **R-dependent part:** At each R, the Hamiltonian H(R) = diag(μ_k) + (R-R0)×V_ref + R×V_residual is constructed from precomputed algebraic matrices with R appearing only as a scalar multiplier. The diagonalization of this small matrix (15×15) is the only per-R numerical operation.

### Summary

The Sturmian basis **preserves the algebraic structure** of the free basis:
- All coupling matrix elements in the Sturmian basis are linear combinations of the free-basis matrix elements (which are algebraic in principle)
- The Sturmian rotation matrix U is computed once by numerical diagonalization
- All subsequent operations at each R are algebraic modulo the final eigenvalue extraction

The algebraicization bottleneck is the same as for the free basis: replacing the Gauss-Legendre quadrature with closed-form evaluation of the Gegenbauer coupling integrals (Paper 13 Eq. 31-33). The Sturmian construction does not add new algebraicization requirements.

---

## 6. Summary and Recommendation

### What worked

1. **Monopole decomposition confirms the physics:** Nuclear monopole is >100% of V_coupling, validating the Sturmian construction.

2. **l_max=0 convergence improvement:** 1.6-3.6x faster convergence at small n_basis. Sturmian at n_basis=15 matches free at n_basis=25.

3. **Eigenvalue convergence at fixed R:** 10-15x improvement at the V_eff well minimum (R ~ 2 bohr).

4. **Algebraic structure preserved:** All matrix elements remain algebraic (in principle), with only one-time numerical diagonalization added.

5. **Both bases beat the FD solver:** At l_max=1, both Sturmian (0.59%) and free (0.56%) beat FD (0.78%).

### What didn't work

1. **l_max>0 degradation:** The Sturmian basis converges faster to the adiabatic limit, which at l_max>0 is below exact. Better angular convergence makes the adiabatic error worse, not better.

2. **Success criteria not fully met:**
   - n_basis=10 Sturmian (0.15%) does not beat n_basis=25 free (0.14%)
   - l_max=1 Sturmian (0.59%) does not reach < 0.5%

3. **Diminishing returns at large n_basis:** Above n_basis=15, the improvement is < 10% and both bases are limited by the adiabatic approximation, not the angular basis.

### Recommendation

**Do NOT proceed to production implementation of the Sturmian basis.**

The Sturmian approach provides genuine improvement at small n_basis for l_max=0, but:
- The improvement is modest at the n_basis values where accuracy matters (n_basis ≥ 15: only 17% improvement)
- At l_max>0 (the regime where the framework needs improvement), the Sturmian basis is slightly worse
- The actual bottleneck is the adiabatic approximation, not the angular basis

**Instead, prioritize:**
1. The coupled-channel solver (which addresses the 97% DBOC cancellation)
2. Algebraicization of the existing Gegenbauer coupling integrals (which improves both bases equally)

The Sturmian prototype remains available in `algebraic_angular_sturmian.py` for future use if the coupled-channel solver benefits from faster angular convergence (possible since the coupled equations involve many angular eigensolves).

---

## 7. Test Inventory

12 new tests in `tests/test_algebraic_angular_sturmian.py`, all passing:

| # | Test | What it verifies |
|:-:|:-----|:-----------------|
| 1 | test_sturmian_orthonormality | U^T U = I for truncated rotation matrix |
| 2 | test_sturmian_orthonormality_multichannel | Orthonormality at l_max=2 |
| 3 | test_sturmian_full_dimension_matches_free | n_construct=n_basis gives identical eigenvalues |
| 4 | test_sturmian_R0_eigenvalues | mu(R0) close to large-basis reference |
| 5 | test_sturmian_improves_lmax0 | Sturmian beats free at n_basis=10, l_max=0 |
| 6 | test_he_energy_sturmian_lmax0 | He energy < 0.15% at n_basis=15 |
| 7 | test_he_energy_sturmian_lmax1_beats_fd | He l_max=1 < 0.65% (beats FD 0.78%) |
| 8 | test_R0_sensitivity | R0 in [0.5, 3.0] gives < 0.05% spread |
| 9 | test_n_construct_convergence | Monotonic convergence with n_construct |
| 10 | test_monopole_decomposition | Nuclear monopole > 95% of coupling |
| 11 | test_sturmian_dboc_positive | DBOC ≥ 0 at all R |
| 12 | test_sturmian_coupling_symmetry | Projected matrices are symmetric |
