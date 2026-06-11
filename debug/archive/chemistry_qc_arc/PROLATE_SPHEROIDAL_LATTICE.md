# Prolate Spheroidal Lattice for H2+

**Date:** 2026-03-12
**Status:** PASS (E_min < 1% error, zero free parameters)

## Motivation

The S3 graph Laplacian (Paper 7) solves atoms exactly because the
hydrogen Schrodinger equation separates in momentum space via Fock's
stereographic projection. For molecules, the analogous natural
coordinate system is prolate spheroidal coordinates (xi, eta), where
the two-center Coulomb problem separates.

**Key question:** Can we build a "molecular Fock projection" by
discretizing the Laplace-Beltrami operator in prolate spheroidal
coordinates?

## Method

The H2+ Schrodinger equation separates in prolate spheroidal coords:

```
xi equation:  d/dxi[(xi^2-1)dF/dxi] + (A + a*xi - c^2*xi^2)F = 0
eta equation: d/deta[(1-eta^2)dG/deta] + (-A + c^2*eta^2 + b*eta)G = 0
```

where:
- a = R*(Z_A+Z_B), b = R*(Z_B-Z_A)
- c^2 = -R^2*E_elec/2 (positive for bound states)
- A = angular separation constant (matches the two equations)

**Solver:**
1. Angular equation: spectral method (Legendre basis, exact to machine precision)
2. Radial equation: self-adjoint FD with Neumann BC at xi=1, Dirichlet at xi_max
3. Root-finding: brentq in c^2 until the radial top eigenvalue = 0

**Zero free parameters.** All edge weights come from the metric. Only inputs: Z_A, Z_B, R, grid resolution.

## Results

### Single point (R=2.0 bohr, exact E_total = -0.6026 Ha)

| N_xi | c^2 | E_total | Error |
|------|-----|---------|-------|
| 500 | 2.114 | -0.5568 | 7.60% |
| 1000 | 2.156 | -0.5779 | 4.11% |
| 2000 | 2.179 | -0.5893 | 2.21% |
| 5000 | 2.193 | -0.5965 | 1.01% |
| 10000 | 2.198 | -0.5990 | 0.59% |

Convergence: ~O(1/N_xi) due to FD discretization of the radial equation.

### PES and Spectroscopic Constants (N_xi=8000)

| Quantity | Computed | Exact | Error |
|----------|----------|-------|-------|
| R_eq | 2.001 bohr | 1.997 | **0.21%** |
| E_min | -0.5984 Ha | -0.6026 | **0.70%** |
| D_e | 0.0984 Ha | 0.1026 | 4.09% |
| k | 0.1032 | ~0.10 | ~3% |

### PES Near Minimum

```
  R(bohr)  E_total(Ha)
    1.50   -0.578321
    1.60   -0.586876
    1.70   -0.592584
    1.80   -0.596096
    1.90   -0.597909
    2.00   -0.598402  <-- minimum region
    2.10   -0.597870
    2.20   -0.596545
    2.50   -0.589446
    3.00   -0.573047
    4.00   -0.541195
    5.00   -0.518967
    6.00   -0.505785
```

### Tests: 11/11 pass

- Basic solver sanity (4 tests)
- Grid convergence (2 tests)
- PES shape and spectroscopic constants (3 tests)
- Dissociation limit (1 test)
- Heteronuclear HeH2+ (1 test)

## Comparison with 1D Axial Approach (v1)

| Method | R_eq err | E_min err | Free params |
|--------|----------|-----------|-------------|
| 1D axial (rho0=1.0) | 33% | 30% | 1 (rho0) |
| 1D axial (rho0=1.2) | 55% | 10% | 1 (rho0) |
| **Prolate spheroidal** | **0.21%** | **0.70%** | **0** |

The prolate spheroidal approach is 40x more accurate and parameter-free.

## What Failed: 2D Grid Approach

Initially tried solving the UNSEPARATED 2D (xi, eta) eigenvalue
problem on a tensor product grid. This converged extremely slowly
(28% error at 50x25 = 1250 unknowns) due to:

1. **Wrong BCs**: Dirichlet at eta=+/-1 forces psi=0 at nuclei (wrong for sigma_g)
2. **Resolution**: Uniform grid wastes points in the asymptotic region
3. **Coupling**: The 2D operator mixes xi and eta, requiring much finer grids

The separated approach (two 1D equations + root-finding) is vastly
more efficient because it exploits the exact separation of variables.

## Significance

This confirms the hypothesis from the user's proposal:

> "Every separable quantum system has a natural lattice -- the lattice
> that discretizes the Laplace-Beltrami operator in the coordinate
> system where separation occurs."

- **Atoms:** S3 lattice (from SO(4) symmetry, Fock projection)
- **Diatomics:** Prolate spheroidal lattice (from two-center separation)

The GeoVac program isn't about S3 specifically -- it's about finding
the natural conformal geometry for each class of quantum systems.

## Files

- `geovac/prolate_spheroidal_lattice.py` -- ProlateSpheroidalLattice class
- `tests/test_prolate_h2plus.py` -- 11 tests (all passing)
- `debug/data/prolate_h2plus_pes.txt` -- PES data
- `debug/validate_prolate_h2plus.py` -- Full validation script

## Remaining Work

1. **Higher grid order:** Replace 2nd-order FD with 4th-order or spectral methods for the radial equation to achieve faster convergence
2. **Excited states:** sigma_u, pi states (m=1) -- already supported via m parameter
3. **Two-electron extension:** H2 would require a 2-electron prolate spheroidal CI
4. **Connection to S3:** Establish the formal map between the prolate spheroidal manifold and the S3 topology for the united-atom limit
