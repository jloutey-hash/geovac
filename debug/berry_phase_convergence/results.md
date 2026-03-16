# Berry Phase Convergence Study

**Date:** 2026-03-15
**Status:** Complete (7/7 tests passing)
**Tests:** `tests/test_berry_phase.py` -- TestBerryPhaseComputation (5), TestExponentConvergence (2)

---

## 1. Definition

A **plaquette** on the geometric lattice is the closed loop:

    |n,l,m> --T+--> |n+1,l,m> --L+--> |n+1,l,m+1> --T---> |n,l,m+1> --L---> |n,l,m>

The **plaquette phase** (discrete holonomy) is:

    theta = ln(w_T+) + ln(w_L+) - ln(w_T-) - ln(w_L-)

For topological weights w = 1/(n1*n2):
- T+ from (n,l,m) to (n+1,l,m): w = 1/(n(n+1))
- L+ from (n+1,l,m) to (n+1,l,m+1): w = 1/(n+1)^2
- T- from (n+1,l,m+1) to (n,l,m+1): w = 1/(n(n+1))
- L- from (n,l,m+1) to (n,l,m): w = 1/n^2

The T+ and T- weights are identical, so they cancel exactly:

    theta(n) = ln(1/(n+1)^2) - ln(1/n^2) = -2*ln((n+1)/n)

This is **independent of l and m** -- purely a radial (n-dependent) quantity.

---

## 2. Analytical Formula

    theta_per_plaq(n) = -2*ln(1 + 1/n)

Taylor expansion for large n:

    = -2/n + 1/n^2 - 2/(3n^3) + ...

Leading order: ~2/n, so the per-plaquette phase decays as **k = 1 exactly**.

---

## 3. Total Berry Phase at Shell n

The number of valid plaquettes at shell n is:

    N_plaq(n) = sum_{l=0}^{n-1} (2l) = n(n-1)

(l=0 contributes 0 plaquettes since there's only m=0; l=1 contributes 2; etc.)

Total Berry phase at shell n:

    Theta_total(n) = n(n-1) * (-2*ln((n+1)/n))

For large n: ~ n(n-1) * 2/n = 2(n-1) ~ 2n. The total phase **grows** linearly.

---

## 4. Convergence Study

### Per-plaquette exponent k_per (fitted from |theta| = A*n^(-k))

| max_n | k_per | |k-1| | k_asymp (n>10) |
|:-----:|:-----:|:-----:|:--------------:|
| 10 | 0.899 | 0.101 | -- |
| 20 | 0.928 | 0.072 | -- |
| 50 | 0.955 | 0.045 | -- |
| 100 | 0.970 | 0.030 | -- |
| 150 | 0.976 | 0.024 | 0.988 |
| 200 | 0.980 | 0.020 | 0.990 |

**k_per converges monotonically to 1.0 from below.** The sub-unity values at finite nmax are due to the 1/n^2 correction term in the Taylor expansion.

### Numerical vs analytical agreement

    max |numerical - analytical| per plaquette: 6.84e-15 (machine precision)

The plaquette Berry phase computation is **exact** -- it reproduces the analytical formula to machine precision.

---

## 5. Detailed Phase Values (max_n=100)

| n | Total phase | N_plaq | Per-plaquette | Analytical |
|:---:|:---------:|:------:|:------------:|:----------:|
| 2 | -1.6219 | 2 | 0.8109 | 0.8109 |
| 5 | -7.2929 | 20 | 0.3646 | 0.3646 |
| 10 | -17.156 | 90 | 0.1906 | 0.1906 |
| 20 | -37.081 | 380 | 0.0976 | 0.0976 |
| 50 | -97.033 | 2450 | 0.0396 | 0.0396 |
| 99 | -195.02 | 9702 | 0.0201 | 0.0201 |

---

## 6. Discrete Curvature

The second difference of the per-plaquette phase (discrete curvature) scales as:

    d^2(theta)/dn^2 ~ n^(-2.94)

converging to n^(-3) (the third term in the Taylor expansion).

---

## 7. Comparison with Paper 1

Paper 1 reports k ~ 2.113 for a "phase angle" on the lattice. The plaquette Berry phase defined here gives k = 1.0 exactly. These are **different quantities**:

- **Paper 1's phase**: Likely the argument of a complex eigenvalue or a spectral angle of the graph Laplacian, which involves the full Hamiltonian H = kinetic_scale*(D-A+W) and its eigenvalue structure.
- **This study's phase**: The discrete holonomy from circulating around a single plaquette in the (n,m) plane, computed purely from adjacency matrix weights.

The plaquette Berry phase is a **geometric** property of the graph (curvature of the connection), while Paper 1's phase is a **spectral** property (related to eigenvalues). The two are connected through the Chern-Gauss-Bonnet theorem on the discrete bundle, but they are not identical.

### What the plaquette phase tells us

1. **The topological weights 1/(n1*n2) create non-zero curvature** on the lattice. Binary weights (all 1.0) give zero curvature everywhere.

2. **The curvature is purely radial** -- it depends only on n, not on l or m. This is consistent with the SO(4) symmetry of the hydrogen atom: the curvature of the S^3 topology is isotropic in the angular directions.

3. **The total holonomy grows linearly** with shell number n. This means the lattice accumulates phase as you move outward -- consistent with the physical picture of the conformal factor Omega = 2p_0/(p^2 + p_0^2) creating more "stretching" at higher n.

4. **The per-plaquette curvature decays as 1/n** -- the lattice becomes flatter at large n, approaching the flat-space (continuum) limit.

---

## 8. Energy Corrections

Since k = 1 exactly (not 2), the Berry phase does not produce corrections at the same order as the energy levels E_n = -1/(2n^2). The per-plaquette phase 2/n is at a different order than the eigenvalue corrections.

If one defines an "energy correction" from the Berry phase as:

    dE_Berry ~ theta_per_plaq^2 / (2n^2) ~ 4/n^2 * 1/(2n^2) = 2/n^4

this would be a fourth-order correction, much smaller than relativistic corrections (~alpha^2/n^3) or QED corrections (~alpha^3*ln(alpha)/n^3). It is numerically negligible for spectroscopy.

---

## 9. Files Created

| File | Description |
|------|-------------|
| `geovac/berry_phase.py` | Berry phase module (4 functions) |
| `tests/test_berry_phase.py` | Validation tests (7 tests, 7/7 passing) |
| `debug/berry_phase_convergence/results.md` | This analysis |

---

## 10. Conclusion

The discrete plaquette Berry phase on the geometric lattice with topological weights is **analytically solvable**: theta(n) = -2*ln((n+1)/n), giving a per-plaquette decay exponent k = 1.0 exactly. This is distinct from Paper 1's spectral phase angle (k ~ 2.113), which involves the full Hamiltonian eigenvalue structure rather than the pure graph geometry.

The key insight is that the 1/(n1*n2) edge weights create a **non-zero but analytically tractable curvature** on the lattice fiber bundle. This curvature is purely radial (n-dependent), isotropic in angular quantum numbers (l, m), and decays as 1/n -- the lattice becomes flat in the continuum limit, consistent with the Fock projection recovering flat R^3 from curved S^3.
