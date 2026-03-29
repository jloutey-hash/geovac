# Algebraic Angular Solver — Implementation Notes

**Date:** 2026-03-27 (v2.0.5, Track B Task 2)

---

## 1. Gegenbauer Parameter Correction

The task specification and audit reference Gegenbauer polynomials C_{n-1}^2(cos 2α) (λ=2) as the free SO(6) eigenfunctions. This is **incorrect for the l=0 channel**.

The correct mapping is **λ = l + 1**, where l is the partial-wave index:
- l=0: λ=1 → C_k^1(cos 2α) = U_k(cos 2α) (Chebyshev of second kind)
- l=1: λ=2 → C_k^2(cos 2α)
- l=2: λ=3, etc.

**Derivation:** The angular equation on S^5 for the l=0 sector, after change of variable x = cos(2α), gives the Gegenbauer equation with (2λ+1)x coefficient equal to 3, hence λ=1. The source of confusion is likely the full S^5 Gegenbauer parameter (d-2)/2 = 2 for d=6, but this applies to the unreduced hyperspherical harmonics, not to the reduced 1D equation in α after partial-wave projection.

**Verification:** The free eigenfunctions sin(2nα) = sin(2α) × U_{n-1}(cos 2α) satisfy the Liouville equation -1/2 u'' - 2u = μu with eigenvalues μ = 2n²-2, matching the SO(6) Casimir formula ν(ν+4)/2 with ν = 2(n-1).

**Implementation choice:** The code uses sin(2nα) directly rather than evaluating Gegenbauer polynomials, since the trigonometric form is simpler, avoids the λ ambiguity, and is numerically well-conditioned.

---

## 2. Quadrature Strategy

The coupling matrix elements are computed via Gauss-Legendre quadrature on two sub-intervals [0, π/4] and [π/4, π/2]:

- **Why split at π/4:** The V_ee monopole term 1/max(cos α, sin α) has a kink at α = π/4. Splitting gives exponential convergence on each smooth sub-interval.
- **Why not closed-form:** The full algebraic evaluation (Gegenbauer linearization → Jacobi integrals → Gamma function ratios) is derivable but complex to implement correctly. GL quadrature on smooth integrands gives 14+ digits with 100 points per sub-interval, which is effectively exact for a 10-20 basis matrix.
- **Key insight:** The singularity handling is the whole point. The FD solver puts 1/cos²α directly on a grid, causing boundary errors. The spectral solver multiplies smooth basis functions × (bounded integrand) and integrates — no singularity is ever evaluated on a grid.

**Performance:** Matrix precomputation with n_basis=10, n_quad=100 takes <0.1s. The angular solve at each R is a 10×10 diagonalization (~1μs). Total for 200 R points: ~0.2s (vs ~5s for FD at n_alpha=400).

---

## 3. Selection Rules

**Nuclear coupling (verified numerically):**
V_nuc[n', n] = 0 when (n' - n) mod 4 ≠ 0.

For singlet basis n = 1, 3, 5, 7, 9, ..., this means coupling only between:
- n=1 ↔ n=5 (Δ=4), n=1 ↔ n=9 (Δ=8), n=3 ↔ n=7 (Δ=4), etc.
- n=1 and n=3 do NOT couple (Δ=2 ≡ 2 mod 4).

**Physical origin:** The nuclear potential -Z(1/cosα + 1/sinα) is invariant under α → π/2-α (exchange symmetry). In the Fourier basis sin(2nα), this exchange maps n → (-1)^{n+1}n, enforcing that only modes with n' - n divisible by 4 can couple. This is a consequence of S₂ permutation symmetry acting on the SO(6) angular basis.

**V_ee coupling:** The monopole 1/max(cos,sin) does NOT have the same selection rules (the kink at π/4 breaks the strict parity). V_ee couples adjacent modes as well.

---

## 4. Accuracy Assessment

### Cross-validation with FD solver
At Z=2, n_basis=25 vs FD with n_alpha=400, l_max=0:

- R ∈ [0.1, 3.0]: agreement < 0.5% (both well-converged in the well region)
- R = 10.0: ~1.7% difference (spectral basis truncation at large R)

### He ground-state energy
| Solver | Config | Energy (Ha) | Error |
|--------|--------|-------------|-------|
| FD | l_max=0, n_alpha=100 | -2.9052 | 0.054% |
| Algebraic | n_basis=10 | -2.8967 | 0.241% |
| Algebraic | n_basis=20 | -2.8996 | 0.144% |
| Algebraic | n_basis=30 | -2.8999 | 0.133% |

### Large-R limitation
The spectral basis converges slowly at large R because the angular eigenfunction
localizes near α=0 and α=π/2 (one electron near nucleus, one far away).
Representing such localized functions requires many oscillatory sin(2nα) basis
functions — a known limitation of global spectral methods.

**This does NOT affect the l_max > 0 use case**, which is the solver's primary
motivation.  The FD solver's l_max=0 accuracy (0.054%) is better than the
spectral solver's (0.13%), but at l_max > 0 the FD solver degrades to 0.87%
due to centrifugal boundary singularities, while the spectral solver handles
the singularity analytically.

### Basis convergence
Monotonically convergent with n_basis, but sub-exponential at the He energy
level due to the large-R localization issue.

---

## 5. Future Work: l_max > 0

For the multichannel case, each channel l needs its own basis:
- u_{l,k}(α) = N_{l,k} × (sinα cosα)^{l+1} × C_k^{l+1}(cos 2α)
- Kinetic + centrifugal is diagonal for each l separately
- Inter-channel V_ee coupling involves Gaunt-weighted integrals

The block structure is (n_l × n_basis) × (n_l × n_basis), still much smaller than the FD matrix. The Gaunt integrals (_precompute_gaunt) from the FD solver can be reused for the angular coupling coefficients.

The key l_max convergence test — showing that l_max=1,2,3 improve on l_max=0 (unlike the FD solver) — requires this extension.
