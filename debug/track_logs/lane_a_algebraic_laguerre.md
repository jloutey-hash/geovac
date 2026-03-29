# Lane A — Algebraic Laguerre Matrix Elements (Level 2)

**Date:** 2026-03-28
**Track:** Level 2 algebraicization — eliminate quadrature from spectral radial solver
**Status:** COMPLETE (m=0)

## Summary

Replaced all Gauss-Laguerre quadrature in the Level 2 spectral radial solver with
algebraic three-term recurrence relations for m=0 (sigma states). The radial solver
now constructs H and S matrices without any numerical integration.

## Mathematical Approach

All matrix elements reduce to **Laguerre moments**:
```
M_k[i,j] = ∫₀^∞ x^k · L_i(x) · L_j(x) · e^{-x} dx
```

These are computed exactly via the three-term recurrence:
```
x·L_n(x) = -(n+1)L_{n+1}(x) + (2n+1)L_n(x) - n·L_{n-1}(x)
```

which gives closed-form band matrices:
- M0 = identity (orthonormality)
- M1 = tridiagonal: M1[i,i] = 2i+1, M1[i,i±1] = -(i+1) or -i
- M2 = pentadiagonal: M2[j,j] = 6j²+6j+2, etc.

### Matrix element derivation

**Overlap:** S_mn = δ_mn/(2α)

**Potential:** q(ξ) = A + aξ - c²ξ², with ξ = 1 + x/(2α), becomes
q(x) = q₀ + q₁x + q₂x². Then V = (1/(2α))[q₀M₀ + q₁M₁ + q₂M₂].

**Kinetic:** Define f_n = L_n + 2∑_{j<n} L_j (Laguerre expansion of the basis
derivative structure). Then K = -(1/2)·F·M₁·Fᵀ - (1/(8α))·F·M₂·Fᵀ
where F[n,j] = 2 for j<n, 1 for j=n.

## Results

| Metric | Value |
|--------|-------|
| Overlap max diff (algebraic vs quadrature) | 4.4e-15 (machine precision) |
| Hamiltonian max diff | 4.1e-12 |
| Energy diff at R=2.0 | ~1e-14 Ha |
| Max PES diff (20 points) | 1.6e-14 Ha |
| Speed improvement | ~1.6x (eliminates GL root computation) |
| R_eq agreement | identical |
| E_min agreement | identical |

## m≠0 Status

The m²/(ξ²-1) term introduces a 1/x singularity in Laguerre x-space. This means
matrix elements involve ∫ L_i·L_j·e^{-x}/x dx, which diverges at x=0 and has
no finite Laguerre moment expansion. The algebraic method raises
NotImplementedError for m≠0, and the solver falls back to quadrature automatically.

The m=0 case covers the H₂⁺ ground state (1sσ_g) and all σ excited states, which
are the primary use cases.

## Files Modified

- `geovac/prolate_spheroidal_lattice.py`:
  - Added `_laguerre_moment_matrices(N)` — builds M0, M1, M2
  - Added `_build_laguerre_matrices_algebraic(n_basis, alpha, A, a_param, c2, m)`
  - Refactored quadrature code into `_build_laguerre_matrices_quadrature()`
  - Added `matrix_method` parameter ('quadrature' default, 'algebraic')
  - Updated `scan_h2plus_pes()` to accept `matrix_method`

- `tests/test_spectral_radial.py`:
  - Added `TestLaguerreMomentMatrices` (6 tests)
  - Added `TestAlgebraicVsQuadratureMatrices` (5 tests)
  - Added `TestAlgebraicSolverEnergies` (8 tests)

## Test Count

- 37 original tests: all pass (zero regressions)
- 19 new tests: all pass
- Total: 56 tests

## Escalation Items

None. This is a numerical improvement within existing channel structure.
No changes to quantum numbers, selection rules, or discrete structure.
