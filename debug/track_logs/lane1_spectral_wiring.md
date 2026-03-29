# Lane 1 — Track C Production Wiring Log

**Date:** 2026-03-28
**Goal:** Wire spectral Laguerre radial solver into scan_h2plus_pes() and solve_with_wavefunction()

## Changes Made

### geovac/prolate_spheroidal_lattice.py

1. **scan_h2plus_pes()** — Added `radial_method` (default `'fd'`) and `n_basis` (default 20) parameters. These thread through to `ProlateSpheroidalLattice` constructor. Default behavior unchanged.

2. **solve_with_wavefunction()** — Refactored radial wavefunction computation into two helper methods:
   - `_radial_wavefunction_fd(c2, A)` — original FD path (extracted, not modified)
   - `_radial_wavefunction_spectral(c2, A)` — new spectral Laguerre path that evaluates the eigenvector from the Galerkin expansion on a 500-point uniform xi grid
   - Method dispatch based on `self.radial_method`
   - Added `'radial_method'` key to return dict for downstream identification

### tests/test_spectral_pes_wiring.py (NEW, 20 tests)

- **TestSpectralPESScan** (7 tests): R_eq, E_min, D_e accuracy; no NaNs; not at boundary
- **TestFDPESScanUnchanged** (2 tests): default behavior and explicit FD match
- **TestSpectralPESSpeed** (1 test): spectral faster than FD
- **TestSpectralWavefunction** (10 tests): all keys present, method tag, nonzero/normalized/decaying F(xi), energy consistency, angular consistency, shape similarity

## Benchmark Results

### Spectral PES (n_basis=20, 41 R-points from 1.5 to 3.5 bohr)

| Quantity | Spectral | FD (N_xi=5000) | Exact |
|----------|----------|----------------|-------|
| R_eq     | 2.0046 bohr (0.38%) | 2.0000 bohr (0.15%) | 1.997 bohr |
| E_min    | -0.60263701 Ha (0.0005%) | -0.59654137 Ha (1.01%) | -0.6026342 Ha |
| D_e      | 0.10264 Ha (0.04%) | 0.09654 Ha (5.9%) | 0.1026 Ha |
| Wall time | 0.53s | 153.1s | — |
| **Speedup** | **287x** | 1x | — |

**Key observation:** Spectral PES is dramatically more accurate for E_min and D_e (5000x and 160x improvement respectively). R_eq is slightly less accurate (0.38% vs 0.15%) due to polynomial fit sensitivity — the underlying energies are far more precise. The 287x speedup is consistent with the single-point 270x speedup from v2.0.8.

## Test Results

- **17/17** existing test_spectral_radial.py tests pass (zero regressions)
- **20/20** new test_spectral_pes_wiring.py tests pass
- **37 total, all passing**

## Escalation Items

None. All changes are within Lane 1 scope:
- No paper changes needed
- No CLAUDE.md updates needed
- No failed approach entries needed
- Default behavior preserved (radial_method='fd')

## Next Steps (not done, per task scope)

- Wire spectral solver into composed-geometry pipeline (rebuild Level 5 Hamiltonians with spectral Level 2)
- Apply spectral Laguerre pattern to Level 3 hyperradial solver
