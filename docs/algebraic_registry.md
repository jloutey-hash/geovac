# Algebraic registry (extracted from CLAUDE.md S12, live document)

> Extracted verbatim from CLAUDE.md on 2026-06-10 (v3.110.0 state) during compaction round 2. This file is the canonical archive; the CLAUDE.md section now holds only the compact working form.

## S12 (update here, not in CLAUDE.md)

## 12. Algebraic Registry

Tracks which matrix elements at each level are computed algebraically vs numerically. Status: **algebraic** (closed-form from quantum numbers), **algebraic (implicit)** (defined by polynomial equation P=0 with known coefficient ring; pointwise diag is computational convenience, not mathematical necessity), **algebraic-pending** (algebraic route identified but production code still uses quadrature), **numerical-required** (no known algebraic replacement).

### Level 3 (Hyperspherical — He)

| Matrix Element | Status | Notes |
|:---------------|:------:|:------|
| SO(6) Casimir eigenvalues | algebraic | Exact: ν(ν+4)/2 from representation theory |
| Centrifugal barrier | algebraic | Diagonal in Gegenbauer spectral basis (confirmed Track B, v2.0.6) |
| Nuclear coupling | algebraic | Partial harmonic sums (Eqs. 31-32, Paper 13) |
| V_ee coupling | algebraic-pending | Split-region Legendre structure confirmed; GL quadrature still used in production |
| Centrifugal matrix elements | algebraic | Diagonal in Gegenbauer basis (v2.0.6) |
| Hyperradial overlap (S) | algebraic | Pentadiagonal M2 moment matrix from three-term Laguerre recurrence. S = M2/(8α³). Machine-precision agreement with quadrature (< 1e-14 relative). (Track H, v2.0.10) |
| Hyperradial kinetic (K) | algebraic | Derivative kernel B_n = -n/2 L_{n-1} + 1/2 L_n + (n+1)/2 L_{n+1} (tridiagonal). K = bbᵀ/(4α). Machine-precision agreement (< 1e-14 relative). 11× build speedup. (Track H, v2.0.10) |
| Hyperradial potential (V_eff) | algebraic at l_max=0; numerical-required at l_max≥1 | l_max=0: V_eff algebraic via Stieltjes moment decomposition (single transcendental seed e^a·E₁(a), Track P2 v2.0.13). l_max≥1: μ(R) is algebraic over Q(π,√2) (Track P1), but V_eff integrals are non-elementary (radical obstruction: √Δ(R) at l_max=1, Cardano at l_max=2). Quadrature required. |
| Hellmann-Feynman P-matrix | algebraic | Exact from R-independent dH/dR (v2.0.6) |
| Q-matrix (second derivative coupling) | algebraic | Exact Q = PP + dP/dR computed from Hellmann-Feynman quantities (v2.0.6) |
| Coupled-channel radial solve | numerical-required | Coupled ODE integration on R grid |
| Adiabatic potential curves U(R) | algebraic (implicit) | μ(R) satisfies P(R,μ)=0, polynomial in both R and μ, coefficients in Q(π,√2)[R]. Degree L+1 in both variables at l_max=L. Point-by-point diagonalization is computational convenience. (Track P1, v2.0.12) |
| TC double commutator | algebraic | [[H, R sin(α)], R sin(α)] = -kron(S_R, I + cos²α). Exact closed form, terminates at 2nd order for 2 electrons. The cos²α matrix is block-diagonal in l (Y_l orthogonality). |
| TC single commutator | algebraic | [H, R sin(α)] = (-d/dR)⊗sin(α) + (1/R)⊗[Λ²+15/8, sin(α)] + I⊗[C, sin(α)]. D_R via Laguerre derivative (x weight, anti-symmetric). Angular commutator via IBP + V_extra decomposition. |
| Angular Casimir decomposition | algebraic | V_extra = diag(casimir) - IBP_free, where IBP_free[j,k] = <χ_j'|χ_k'>. Extracts the S⁴ metric cotangent-derivative terms algebraically from the known eigenvalues. Angular derivatives from d/du C_n^λ(u) = 2λ C_{n-1}^{λ+1}(u). |
| sin/sin²/cos² angular matrices | algebraic | Block-diagonal in l by Y_l(θ₁₂) orthogonality. Cross-channel elements are zero because sin(α) preserves θ₁₂ quantum number l. (v2.8.2) |

### Level 2 (Prolate Spheroidal — H₂⁺)

| Matrix Element | Status | Notes |
|:---------------|:------:|:------|
| Angular η-eigenfunctions | algebraic | Legendre spectral basis (n_basis=50), same as Mitnik et al. angular component |
| Azimuthal m | algebraic | Exact separation of variables, integer quantum number |
| Radial ξ-solver (m=0) | algebraic | Ordinary Laguerre basis (n_basis=20, α-adapted). 250× dimension reduction vs FD. Three-term recurrence for ALL matrix elements — zero quadrature. Machine-precision agreement with quadrature (< 1e-14). (v2.0.9) |
| Radial ξ-solver (m≠0) | algebraic | Associated Laguerre basis L_n^{|m|}(x) with weight x^|m|·e^{-x}. Partial-fraction decomposition of centrifugal 1/x singularity into lowered moment M_{-1} (algebraic, DLMF 18.9.13) + Stieltjes integral J (three-term recurrence). Single transcendental seed e^a·E₁(a); all other elements algebraic. Associated basis converges faster than ordinary for m=1 (stable by N=10). (Track J, v2.0.10) |
| Separation parameter c² | numerical-required | Iterative root-finding (Brent method) for self-consistency between angular and radial equations |
| Coupled-channel ceiling | characterized | Error floor 0.19-0.20% from adiabatic approximation; algebraic convergence ~l_max⁻²; 3 channels sufficient; n_basis and R-grid converged. Sub-0.1% requires non-adiabatic (2D variational) solver. (v2.0.8) |

### Level 4 (Mol-Frame Hyperspherical — H₂)

| Matrix Element | Status | Notes |
|:---------------|:------:|:------|
| Hyperradial overlap (S) | algebraic | Same Laguerre three-term recurrence as Level 3. Pentadiagonal M2 moment matrix. (Track I, v2.0.10) |
| Hyperradial kinetic (K) | algebraic | Same derivative kernel as Level 3. Tridiagonal expansion. (Track I, v2.0.10) |
| Hyperradial potential (U_eff) | numerical-required | Adiabatic curve from angular sweep. μ(ρ) is piecewise-smooth in ρ (NOT algebraic — split-region Legendre expansion creates kinks at min/max boundaries). Structurally different from Level 3. (Track S, v2.0.13) |
| Angular eigenvalue sweep (FD) | numerical-required | Point-by-point diagonalization of H_ang at ~130 ρ-values. FD: 1000×1000 matrices, 99% of wall time. Structurally irreducible: no global P(ρ,μ)=0 exists (Track S). |
| Angular eigenvalue sweep (spectral) | numerical-required | Jacobi polynomial basis: 50×50 matrices, 269× speedup, ~50% of wall time. SO(6) Casimir free spectrum + precomputed V_ee coupling. `angular_method='spectral'`. Optimal approach given Track S's structural finding. (Track K, v2.0.11) |
| 2D tensor product assembly | numerical-required | H_ang evaluated at each quadrature point for (R_e, α) tensor product. Dense kronecker assembly. |
| Nuclear coupling (split-region Legendre) | algebraic | Gaunt integrals, exact via 3j triangle inequality (Paper 15). |
| Cusp factor f(α) | negative result | Alpha-only cusp factor f(α) = 1 + (R_e/2)sin(2α) does not improve D_e. Cusp is 2D in (α, θ₁₂), not separable in α alone. Slow convergence dominated by θ₁₂ Gegenbauer expansion. (Track U, v2.0.14) |
| 1/r₁₂ graph absorption | structural obstruction | Cannot absorb 1/r₁₂ into conformally weighted S⁵ Laplacian. Green's function singularity mismatch: 1/d³ (S⁵) vs 1/d¹ (Coulomb). Cusp is an embedding exchange constant. (Track W, v2.0.14) |
| Cusp energy correction | corrective (post-processing) | Schwartz partial-wave extrapolation: ΔE_cusp ~ -A/(l_max+2)⁴. A = (10/π)⟨δ³(r₁₂)⟩. He: 0.10% at l_max=2 (from 0.24%). H₂: R-dependent, ~1.7 mHa differential, ~1.0 pp D_e. Basis-dependent, not exchange constant. (Track X, v2.0.14) |

### Level 4N (Full N-Electron Mol-Frame Hyperspherical — LiH)

| Matrix Element | Status | Notes |
|:---------------|:------:|:------|
| SO(12) Casimir eigenvalues | algebraic | Exact from SO(3N) representation theory. Free spectrum validated. (Track AK, v2.0.20) |
| S₄ Young projector | algebraic | Character-based projection onto [2,2] singlet irrep. Optimized channel-space eigendecompose. (Track AJ, v2.0.20) |
| Cross-pair V_ee (Gaunt direction) | algebraic | Gaunt integral structure for angular coupling between electron pairs. (Track AK, v2.0.20) |
| Cross-pair V_ee (hyperangular) | numerical-required | 3D numerical hyperangular integration. ρ-independent (precomputed once). (Track AK, v2.0.20) |
| Angular eigenvalue sweep | numerical-required | Point-by-point diagonalization at each ρ. Spectral basis: 750 dim at l_max=2 (1000× compression from FD). (Track AK, v2.0.20) |
| Hyperradial solve | numerical-required | Same Laguerre spectral basis as Level 4. (Track AJ, v2.0.20) |

### Level 5 (Composed Geometries — LiH, BeH₂, H₂O)

| Matrix Element | Status | Notes |
|:---------------|:------:|:------|
| Z_eff screening function | algebraic | Laguerre spectral expansion: n(r) projected onto L_k(2βr)·r²·e^{-2βr} basis. N_core(r) = C_inf − e^{-X}·P(X), single transcendental seed e^{-2βr}. Production default `zeff_method='spectral_laguerre'`. Auto-fallback to spline for Z≥4 (polynomial cancellation). (Track N, v2.0.12; Track R wiring, v2.0.13) |
| Slater F^k integrals | algebraic | Two implementations: (1) Laguerre moment decomposition (Track O, v2.0.12): F^k = c_A^T · R^k · c_B, single transcendental seed γ(1,x). (2) Hypergeometric R^k evaluator (v2.8.1, threshold-dispatched in 2026-05-08 cleanup): `geovac/hypergeometric_slater.py` with exact Fraction-arithmetic for arbitrary n_max. `compute_rk_float()` is correct at all n via threshold dispatch — n ≤ 4 uses the fast pure-float path (machine precision, ~240 μs/call, bit-identical to pre-fix), n ≥ 5 delegates to `compute_rk_algebraic` (exact Fraction) and casts to float. **Pre-fix the float path silently degraded for n ≥ 6** (rel err 1e-7 at n=6, 1e-4 at n=8, 0.2 at n=10, fully wrong sign at n ≥ 12) due to catastrophic cancellation between O(10^25) Laguerre-product terms summing to O(10^11) — exceeded float64 mantissa. Post-fix performance for uncached calls: n=8: 0.54 s, n=10: 1.9 s, n=12: 6.5 s, n=20: ~1 min; cached via `_CACHE_FLOAT` so paid once per unique quartet. `k > l_a + l_b` Gaunt-violating quartets (which `two_electron_integral` filters at the `k_min/k_max` step) now raise a clear `ValueError` instead of crashing on factorial-of-negative. 38 regression tests in `tests/test_hypergeometric_slater.py` cover n=6..20 across ssss/pppp/spsp/mixed quartets, verified against `compute_rk_exact` ground truth. Validated 144/145 table entries in `casimir_ci.py`; found+fixed F²(2p,2p) typo (43/512→45/512). Corrected systematic grid bias (0.06-0.44% per integral). 8x speedup over grid quadrature for small n; large-n cached-Fraction-cast for correctness. Production default for graph-native CI. |
| Phillips-Kleinman PK | algebraic | Ab initio from core eigenvector projection (Paper 17 Sec IV). Gaussian PK parameters A, B from core screening. |
| Inter-fiber exchange coupling | algebraic | Full 1-RDM exchange via channel overlap S(R) and Slater F^0 (algebraic Laguerre). Bond-bond coupling validated (~0.5 Ha). Lone pair coupling disabled at Z_eff>4 (physics limitation). |
| Cross-center V_ne (multipole) | algebraic | Exact termination at L = 2*l_max by Gaunt selection rules. 33 terms at n_max=2, 168 at n_max=3. m-diagonal. Analytical evaluation via incomplete gamma functions (machine precision, zero quadrature). Replaces grid-based trapezoid. (Track CD, v2.0.39-40) |

### Level 5 — Spin-ful composed (Tier 2)

| Matrix Element | Status | Notes |
|:---------------|:------:|:------|
| Szmytkowski angular σ·r̂, J², L·S, L², σ² | algebraic | Exact eigenvalues in (κ, m_j); integer / half-integer Kronecker deltas |
| Diagonal hydrogenic ⟨r^k⟩, ⟨1/r^k⟩ | algebraic | Bethe-Salpeter rationals; ⟨1/r⟩=Z/n², ⟨1/r³⟩=Z³/[n³l(l+½)(l+1)] (Z³ diverges l=0, suppressed by Kramers) |
| Off-diagonal ⟨n'l'\|r^k\|n l⟩ | algebraic | Direct sympy integration of assoc_laguerre; ~0.1–1s per call |
| Breit-Pauli H_SO in (κ, m_j) | algebraic | Closed form H_SO = −Z⁴α²(κ+1)/[4n³l(l+½)(l+1)], diagonal; exact Kramers at l=0 |
| jj-coupled angular coefficient X_k(κ_a, m_a, κ_c, m_c) | algebraic | sympy wigner_3j, cached, full-Gaunt selection |
| Spinor composed two-body ERI | algebraic | X_k·X_k·R^k factorization, radial via hypergeometric_slater (exact Fraction or machine-float) |
| Relativistic γ = √(1−(Zα)²) | algebraic | Degree-2 algebraic over ℚ(α²) via γ² + (Zα)² = 1; reserved symbol in T1/T2, not bound at Tier 2 |
| Dirac-Coulomb radial ⟨r^{-1}⟩ (all states) | algebraic | Hellmann-Feynman: Z(γn_r+κ²)/(γN_D³), exact for all n_r. NR limit Z/n². (T7, v2.12.0) |
| Dirac-Coulomb radial ⟨r^{-2}⟩, ⟨r^{-3}⟩ (n_r=0) | algebraic | Pochhammer ratios from single-term wavefunction structure. (T7, v2.12.0) |
| Dirac-Coulomb radial ⟨r^{-2}⟩ (all states) | algebraic | Kramers-Pasternak direct integration via confluent hypergeometric polynomial expansion. Works for all n_r, all κ. (v2.18.1) |
| Dirac-Coulomb radial ⟨r^{-3}⟩ (|κ|≥2) | algebraic | Kramers-Pasternak direct integration. p₁/₂ (|κ|=1, n_r≥1) limitation: sympy integral does not converge for near-singular integrand. (v2.18.1) |
| Darwin + mass-velocity α⁴ corrections | algebraic | Closed-form rationals in (Z, n, l) times α²; verified Dirac formula exact. (T8, v2.12.0) |
| Breit retarded radial R^k_BP(n₁l₁,n₂l₂;n₃l₃,n₄l₄) | algebraic | Exact Fraction/sympy arithmetic via r_<^k / r_>^{k+3} kernel; Z³ scaling verified; `geovac/breit_integrals.py`. (v2.18.0) |
| Seeley-DeWitt coefficients a₀, a₁, a₂ on S³ | algebraic | Exact sympy Rational on unit S³ (a₀=a₁=√π, a₂=√π/8). `geovac/qed_vacuum_polarization.py`. (v2.18.2) |
| QED vertex even/odd Dirac split D_even(s), D_odd(s) | algebraic | Exact Hurwitz zeta at quarter-integer shifts: D_even via ζ(s,3/4), D_odd via ζ(s,5/4). At s=4: D_even = π²/2 − π⁴/24 − 4G + 4β(4), D_odd = π²/2 − π⁴/24 + 4G − 4β(4). PSLQ-verified at 80 digits. `geovac/qed_vertex.py`. (v2.19.1) |
| Screened ⟨1/r³⟩ from FrozenCore Z_eff(r) | numerical-required | Radial Schrodinger equation with Z_eff(r)/r potential solved on uniform grid (FD, eigh_tridiagonal). ⟨1/r³⟩ = ∫\|u\|²/r³ dr from normalized wavefunction. Enhancement 12-144× over hydrogenic at Z_eff=2 for Ca/Sr/Ba 2p. `geovac/neon_core.py` (screened_r3_inverse, screened_xi_so). (v2.19.4) |

---


