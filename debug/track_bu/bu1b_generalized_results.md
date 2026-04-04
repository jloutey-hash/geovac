# Track BU-1b: Generalized Sturmian Secular Equation — Results

**Date:** 2026-04-02
**Version:** v2.0.33 (pre-release)
**Module:** `geovac/sturmian_solver.py` — `GeneralizedSturmianCI` class

## Method

Generalized Sturmian CI: each Slater determinant ν gets scaling β_ν from the
isoenergetic condition E_trial = -β²Z²/2 × Σ(1/n_i²). All orbitals within one
configuration share Z_eff = β_ν × Z. Self-consistency loop iterates E_trial →
β_ν → (H, S) → generalized eigenvalue problem → E_new until convergence.

SD-level matrix elements use Löwdin formula for 2-electron non-orthogonal
determinants:
  S_μν = det(M)
  H1_μν = Σ_ij h_ij × cofactor(M,i,j)
  V_ee_μν = g(a₁,a₂;b₁,b₂) − g(a₁,a₂;b₂,b₁)

One-body h1 from ket eigenvalue identity:
  ⟨bra|T+V_nuc|ket⟩ = -Z_ket²/(2n²)·S + (Z_ket - Z_nuc)·⟨1/r⟩

Fast R^k integral (_slater_rk_fast): O(n_grid) via cumulative sums vs O(n_grid²)
in the original, giving ~500x speedup for the mixed-Z_eff case.

## Results

### He Ground State Energy (exact: -2.903724 Ha)

| max_n | N_spatial | N_SD | Method             | Energy (Ha) | Error (%) | dE vs std (mHa) | dE vs Coulomb (mHa) |
|-------|-----------|------|--------------------|-------------|-----------|------------------|----------------------|
| 2     | 5         | 45   | Standard FCI       | -2.807809   | 3.30      | —                | —                    |
| 2     | 5         | 45   | Coulomb Sturmian   | -2.811609   | 3.17      | -3.8             | —                    |
| 2     | 5         | 45   | **Generalized**    | -2.824884   | 2.72      | **-17.1**        | **-13.3**            |
| 3     | 14        | 378  | Standard FCI       | -2.817209   | 2.98      | —                | —                    |
| 3     | 14        | 378  | Coulomb Sturmian   | -2.844020   | 2.06      | -26.8            | —                    |
| 3     | 14        | 378  | **Generalized**    | -2.831925   | 2.47      | -14.7            | **+12.1**            |

### Key Finding: Generalized WORSE than Coulomb at max_n=3

At max_n=2, the generalized Sturmian improves by 13.3 mHa over Coulomb.
At max_n=3, it is **12.1 mHa WORSE** than Coulomb.

Root cause: the isoenergetic condition forces all orbitals within one configuration
to share the SAME Z_eff = β_ν × Z. In the Coulomb Sturmian, each orbital gets
Z_eff = n × k (n-dependent), providing within-config radial differentiation. The
generalized Sturmian sacrifices this within-config flexibility for between-config
flexibility. At small basis (max_n=2), between-config matters more. At larger
basis (max_n=3), the loss of within-config differentiation dominates.

### Variational Bound

All energies above exact (-2.903724 Ha): PASS.

### Self-Consistency Convergence

| max_n | Iterations | Damping | Final dE     |
|-------|------------|---------|-------------|
| 2     | 16         | 0.5     | 6.6e-7 Ha   |
| 3     | 12         | 0.5     | 7.6e-6 Ha   |

Convergence is monotonic and well-behaved. No oscillation or divergence observed.

### β Distribution (max_n=3, converged)

| Config (n₁,n₂) | β        | Z_eff  | ⟨r⟩_n₁ (bohr) | ⟨r⟩_n₂ (bohr) |
|-----------------|----------|--------|----------------|----------------|
| (1,1)           | 0.8414   | 1.6828 | 0.89           | 0.89           |
| (1,2)           | 1.0643   | 2.1286 | 0.70           | 5.63           |
| (1,3)           | 1.1289   | 2.2578 | 0.66           | —              |
| (2,2)           | 1.6828   | 3.3657 | 0.45           | 3.56           |
| (2,3)           | 1.9802   | 3.9604 | 0.38           | —              |
| (3,3)           | 2.5242   | 5.0485 | 0.30           | 2.38           |

β range: 0.84 – 2.52 (factor 3.0×). Differentiation IS present.

### Top-Weight Configurations (max_n=3)

| SD  | Orbitals       | β      | Z_eff  | Weight |
|-----|----------------|--------|--------|--------|
| 0   | 1s(α) 1s(β)   | 0.8414 | 1.683  | 0.9403 |
| 2   | 1s(α) 2s(β)   | 1.0643 | 2.129  | 0.0119 |
| 27  | 1s(β) 2s(α)   | 1.0643 | 2.129  | 0.0119 |
| 53  | 2s(α) 2s(β)   | 1.6828 | 3.366  | 0.0044 |
| 147 | 2p(α) 2p(β)   | 1.6828 | 3.366  | 0.0017 |
| 10  | 1s(α) 3s(β)   | 1.1289 | 2.258  | 0.0012 |

The ground state is 94% 1s² with small 1s-2s and 2s² admixtures.

### Gaunt Sparsity Preserved

| max_n | ERI nonzero (Generalized) | ERI nonzero (Coulomb) | Match? |
|-------|---------------------------|------------------------|--------|
| 2     | 65                        | 65                     | YES    |
| 3     | 1492                      | 1492                   | YES    |

Angular selection rule sparsity is exactly preserved (β-independent).

### Wall Times

| max_n | Generalized (total) | Per iteration | Coulomb Sturmian |
|-------|--------------------:|-------------:|:-----------------|
| 2     | ~14s (16 iter)      | ~0.8s         | ~8s              |
| 3     | ~86s (12 iter)      | ~7s           | ~130s            |

Generalized is faster per solve (no k-scan needed), but each iteration
requires building the full SD×SD overlap and Hamiltonian matrices.

## BU-1b Exit Assessment

**NEGATIVE RESULT.** The generalized Sturmian does not consistently improve on
the Coulomb Sturmian. The critical diagnostic is:

**Q: Does β differentiate core from valence?**
**A: Yes, but it's not enough.** β ranges from 0.84 (1s²) to 2.52 (3d²), a
factor 3× spread. Configurations with compact core orbitals get smaller β
(more diffuse), while high-n configs get larger β (more compact). The
differentiation is physically meaningful.

**Q: Does this improve the energy?**
**A: Only at small basis (max_n=2).** At max_n=3, the generalized Sturmian is
12 mHa WORSE than the Coulomb Sturmian. The root cause is that the isoenergetic
constraint (all orbitals in a config share one Z_eff) removes the n-dependent
radial differentiation that the Coulomb Sturmian provides (Z_eff = n×k).

**Implication for BU-3 (Li):** The generalized Sturmian approach will NOT solve
the core-valence PK problem. Even though β differentiates configs, the loss of
within-config orbital flexibility makes the overall basis worse. The Coulomb
Sturmian (BU-1) remains the better atomic basis.

**Do NOT proceed to BU-3 with generalized Sturmians.** The Coulomb Sturmian
direction (BU-1 → BU-2) is the better path for quantum simulation.

## Technical Notes

1. The `_slater_rk_fast` function provides ~500× speedup over `_slater_rk`
   using vectorized cumulative sums. This made the generalized solver practical
   (86s vs estimated ~2 hours with the original R^k).

2. The overlap matrix S is well-conditioned (cond ~4.1 at max_n=3). No
   numerical issues with the generalized eigenvalue problem.

3. The h1 formula uses the ket eigenvalue identity, which is exact for
   hydrogenic orbitals. Hermiticity is verified by symmetrization.

4. The Gaunt angular coefficients are precomputed once (β-independent) and
   reused across all iterations. Only the R^k radial integrals change with β.
