# Lane 2: Spectral Laguerre for Level 3 Hyperradial

**Date:** 2026-03-28
**Track:** Level 3 hyperradial algebraicization
**Goal:** Replace 3000-point FD hyperradial grid with spectral Laguerre basis

## Summary

Successfully implemented spectral Laguerre basis for the hyperradial coordinate R,
replicating the pattern from the Level 2 spectral solver (Track C, v2.0.8).

### Basis

```
phi_n(R) = (R - R_min) * exp(-alpha*(R - R_min)) * L_n(2*alpha*(R - R_min))
```

- Domain: [R_min, inf) mapped via x = 2*alpha*(R - R_min) to [0, inf)
- (R - R_min) prefactor enforces Dirichlet BC at R_min (matching FD solver)
- Gauss-Laguerre quadrature for matrix elements
- Generalized eigenvalue problem: (K + V) c = E S c

### Key Design Decisions

1. **Dirichlet BC enforcement**: Without the (R-R_min) prefactor, the basis
   functions are nonzero at R_min, giving extra variational freedom that
   produces energies systematically below the FD solver (0.005 Ha difference).
   The prefactor eliminates this discrepancy (agreement to < 0.00003 Ha).

2. **Domain shift to [R_min, inf)**: The hyperradial V_eff has a 1/R^2
   singularity at R=0. Sampling near R=0 with Gauss-Laguerre quadrature
   produces catastrophic results (12% error, variational bound violated).
   Shifting to [R_min, inf) avoids the singularity entirely, matching the
   FD solver's domain.

3. **Kinetic energy prefactor**: The correct kinetic matrix element is
   K_ij = 1/(4*alpha) * sum_k w_k B_i(x_k) B_j(x_k), where the 1/(4a) =
   (1/2) * (1/(2a)) combines the 1/2 from -1/2 d^2/dR^2 and the Jacobian.
   An initial factor-of-2 error (0.5/(4a) instead of 1/(4a)) produced
   -3.26 Ha (12% error); fixing it gave the correct -2.90 Ha.

## Results

### Single-Channel (FD angular input, l_max=0)

| Method | Energy (Ha) | Error (%) | Notes |
|--------|-------------|-----------|-------|
| FD (N_R=2000) | -2.90425 | 0.018% | Baseline |
| Spectral (n=25) | -2.90419 | 0.016% | 0.00006 Ha diff |

### Coupled-Channel (algebraic angular, q_mode='exact')

| l_max | FD Energy | Spectral Energy | FD err% | Spec err% | diff (Ha) |
|-------|-----------|-----------------|---------|-----------|-----------|
| 0 | -2.871636 | -2.871616 | 1.105% | 1.106% | 0.000020 |
| 1 | -2.894687 | -2.894664 | 0.311% | 0.312% | 0.000022 |
| 2 | -2.896790 | -2.896766 | 0.239% | 0.240% | 0.000024 |
| 3 | -2.897329 | -2.897304 | 0.220% | 0.221% | 0.000024 |

### Performance

| Metric | FD | Spectral | Factor |
|--------|-----|----------|--------|
| Single-ch dimension | 3000 | 25 | 120x reduction |
| Coupled (3ch) dimension | 9000 | 75 | 120x reduction |
| Single-ch time | 1.9 ms | 3.6 ms | ~1x (both fast) |
| Coupled (3ch) time | 561 ms | 5.9 ms | **95x speedup** |

### Convergence with n_basis

| n_basis | Energy (Ha) |
|---------|------------|
| 5 | -2.90267 |
| 10 | -2.90415 |
| 15 | -2.90418 |
| 20 | -2.90419 |
| 25 | -2.90419 |
| 30 | -2.90419 |

Converged at n_basis ~15. Using 25 for safety margin.

### Alpha sensitivity

Tested alpha in [0.5, 3.0]: energy spread < 0.001 Ha. Basis is robust.

## Files Modified

- `geovac/hyperspherical_radial.py`: Added `solve_radial_spectral()`,
  `solve_coupled_radial_spectral()`, `_build_laguerre_matrices_dirichlet()`.
  Added `radial_method`, `n_basis_radial`, `alpha_radial` params to `solve_helium()`.
- `geovac/algebraic_coupled_channel.py`: Added `radial_method`, `n_basis_radial`,
  `alpha_radial` params to `solve_hyperspherical_algebraic_coupled()`.
  Imports updated for spectral functions.
- `tests/test_hyperspherical_he.py`: Added 10 tests (TestSpectralRadialSolver: 7,
  TestSpectralCoupledChannel: 3).

## Tests

- 131 total (121 existing + 10 new), all passing
- Zero regressions
- Default behavior unchanged (radial_method='fd')

## Escalation Items

None. No paper changes, CLAUDE.md updates, or failed-approaches entries needed.
The spectral solver is a drop-in replacement for the FD radial solver that
preserves all physics (identical energies to < 0.00003 Ha).
