# LiH FCI Convergence Summary

**Date:** 2026-03-08
**Version:** v0.9.10 (post cross-atom J normalization fix)
**System:** LiH at R=3.015 Bohr, 4 electrons, CP-corrected binding energy
**Experiment:** D_e = 0.0924 Ha (2.515 eV)

## Convergence Table

| nmax | N_SD | D_e^CP (Ha) | Error (%) | Time (s) | Notes |
|------|------|-------------|-----------|----------|-------|
| 2 | 4,845 | 0.270 | 192% | 6 | BSSE > D_e; basis too small |
| 3 | 367,290 | 0.093 | **1.0%** | 82 | First reliable result |
| 4 | 8,214,570 | — | — | >1800 | Infeasible (>5.4 GB, timeout) |
| expt | — | 0.092 | — | — | Huber & Herzberg 1979 |

## Key Finding

The v0.9.10 normalization fix to `_phi_s_orbital_general` (n >= 3 orbitals)
changed the nmax=3 LiH binding energy from:

- **Before (v0.9.9):** D_e^CP = 0.083 Ha (10% error)
- **After (v0.9.10):** D_e^CP = 0.093 Ha (**1.0% error**)

The improvement comes from correct n=3 cross-atom J values. Previously,
Φ_3s(0) = 113 instead of 1.0, producing unphysical cross-atom repulsion
that was either clipped or distorted the CI wavefunction.

## Scaling Wall

LiH FCI basis size grows as C(4·nmax², 4) ~ nmax⁸:

- nmax 2→3: 76× more SDs, 14× longer
- nmax 3→4: 22× more SDs, exceeds 30 min and 5.4 GB

Reaching nmax=4 requires sparse SD representation or symmetry-restricted CI.
