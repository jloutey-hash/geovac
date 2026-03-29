# Lane 3: Level 3 Coupled-Channel Convergence Study

**Date:** 2026-03-28
**Goal:** Push He coupled-channel toward sub-0.1% error

## Results: q_mode='exact', n_basis=15, n_channels=3

| l_max | Energy (Ha) | Error % | Wall Time | Angular Dim |
|------:|------------:|--------:|----------:|------------:|
| 0 | -2.871636 | 1.105 | 0.6s | 15 |
| 1 | -2.894687 | 0.311 | 0.8s | 30 |
| 2 | -2.896790 | 0.239 | 1.0s | 45 |
| 3 | -2.897329 | 0.220 | 1.4s | 60 |
| 4 | -2.897522 | 0.214 | 1.6s | 75 |
| 5 | -2.897607 | 0.211 | 2.0s | 90 |

## Convergence Rate

Error decrements (l_max 1->2->3->4->5): 0.072, 0.019, 0.007, 0.003
Fit: error ~ 0.20 + 0.27 / l_max^2 (algebraic, not exponential)

## Sensitivity Checks

### n_channels (l_max=5)
- n_channels=3: 0.2107%
- n_channels=5: 0.2101% (delta: 0.0006%)
- Channels 4-5 have weights < 1e-5. 3 channels is converged.

### n_basis (l_max=3)
- n_basis=10: 0.283%
- n_basis=15: 0.220%
- n_basis=20: 0.204%
- n_basis=25: 0.198%
- n_basis=30: 0.195%
- Converges to ~0.19% as n_basis -> infinity.

### Radial grid (l_max=3, n_basis=15)
- N_R_radial=3000: 0.220%
- N_R_radial=5000: 0.221%
- Radial grid fully converged at 3000.

## Assessment

**Sub-0.1% is NOT reachable** with the current coupled-channel P+Q machinery.
The error converges to a structural floor of ~0.19-0.20%, set by the adiabatic
channel truncation (3 channels capture >99.9% of wavefunction weight, but the
adiabatic approximation itself has irreducible error).

The coupled-channel solver reduces single-channel error by ~60% (from 0.56% to
0.22% at l_max=3), but cannot eliminate the adiabatic approximation error.

**Route to sub-0.1%:** Variational 2D solver (Paper 15 Section VI.D) that avoids
the adiabatic approximation entirely. This is consistent with the Track A
diagnosis (adiabatic approximation is the fundamental bottleneck).

## Profiling

Q-matrix closure computation is 61% of angular+coupling time, but total time
is <2s per l_max value. No optimization needed at current l_max range.

## Tests Added

- Test 14: Extended to l_max=1-5 (was 1-3)
- Test 16: l_max=5 error in [0.15%, 0.25%] (convergence floor characterization)
- Test 17: n_channels=5 matches n_channels=3 to <1 mHa (channel convergence)
- All 17/17 tests pass.
