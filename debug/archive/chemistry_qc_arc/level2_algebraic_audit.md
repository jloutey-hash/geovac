# Algebraic Boundary Audit: Level 2 Prolate Spheroidal Solver

**Date:** 2026-03-28 (v2.0.6)
**Scope:** Read-only audit of the algebraic vs. numerical boundary in the Level 2 (H2+) prolate spheroidal solver.
**Goal:** Map every matrix element type, assess algebraic status, identify what can be replaced by spectral/algebraic methods, and survey relevant literature on Sturmian approaches in prolate spheroidal coordinates.

---

## 1. Literature Review

### 1.1 Mitnik, Ancarani & Gasaneo (2021)

**Paper:** "Generalized Sturmian Functions in prolate spheroidal coordinates," Molecular Physics 119(8), e1881179. arXiv: 2006.06616.

**Approach:** Develops Generalized Sturmian Functions (GSF) as a spectral basis for the separated prolate spheroidal equations. The GSF are solutions to an auxiliary Sturm-Liouville problem with a *generating potential* (typically Coulomb) and an energy-dependent eigenvalue parameter. Two computational schemes are presented:

1. **Iterative 1D:** Solves the angular and radial separated equations alternately in their respective GSF bases. The angular equation is expanded in angular GSFs (which reduce to spheroidal harmonics), producing a separation constant A. The radial equation is then expanded in radial GSFs (generalized Laguerre-type functions adapted to the xi coordinate), and the energy eigenvalue is extracted from the radial diagonalization. The two equations are iterated to self-consistency.

2. **Direct 2D:** Builds the full Hamiltonian matrix in a tensor product basis of angular x radial GSFs and diagonalizes directly. No iteration needed, but larger matrix.

**Key mathematical ingredients:**
- Angular basis: associated Legendre polynomials (same as GeoVac's current angular solver)
- Radial basis: Coulomb Sturmian functions adapted to the xi coordinate, with built-in correct asymptotic decay exp(-alpha*xi) and regularity at xi=1
- Matrix elements: computed without numerical derivatives (the GSF basis construction handles boundary conditions intrinsically). The overlap and Hamiltonian matrix elements between GSFs involve integrals of products of Sturmian functions times polynomial/rational potentials, which can be evaluated via recurrence relations

**Key results:**
- Achieves very accurate results for H2+ ground and excited states with minimal basis sets (typically 5-15 radial and 5-10 angular functions)
- No numerical derivatives required -- all matrix elements computed from basis function properties
- Exponential convergence with basis size (spectral method)
- Robust over wide range of R, including atomic limit R -> 0

**What GeoVac already has:**
- The angular solver is already spectral (Legendre basis, n_basis=50). This is essentially the same as Mitnik et al.'s angular component
- The self-consistency loop (Brent root-finding in c^2) serves the same purpose as their iterative 1D scheme
- The separated equation structure is identical

**What would be new:**
- Replacing the FD radial solver with a spectral radial basis (Coulomb Sturmians or generalized Laguerre functions in xi)
- Computing radial matrix elements algebraically via recurrence relations instead of FD stencils
- The direct 2D tensor product approach (not currently used by GeoVac)

### 1.2 Kereselidze, Chkadua, Defrance & Ogilvie (2016)

**Paper:** "Derivation, properties and application of Coulomb Sturmians defined in spheroidal coordinates," Molecular Physics 114(1), 148-161.

See also: Kereselidze & Ogilvie (2018), "The Hydrogen-Atom Problem and Coulomb Sturmian Functions in Spheroidal Coordinates," Advances in Quantum Chemistry.

**Approach:** Derives Coulomb Sturmian amplitude functions directly in prolate spheroidal coordinates, presented in closed algebraic form. The key insight is that the separated radial and angular equations for the two-center Coulomb problem are instances of Heun's confluent equation, and the Sturmian functions are polynomial solutions of this equation.

**Key mathematical ingredients:**
- Spheroidal Sturmian functions defined as polynomial solutions of the confluent Heun equation
- Closed-form expressions relating spheroidal Sturmians to spherical Sturmians via hybridization coefficients that depend on R
- The hybridization (mixing of spherical l-values into spheroidal orbitals) is quantified analytically -- each spheroidal orbital is a known linear combination of spherical harmonics

**Key results:**
- Ground state energy E(1sigma_g) = -1.102634186 a.u. at R_eq with only 6 Legendre polynomials in the angular expansion
- Demonstrates that spheroidal Sturmians are the "most appropriate" basis for diatomic calculations due to natural adaptation to two-center geometry
- Closed-form algebraic expressions for the basis functions
- Strong convergence properties inherited from completeness of the Sturmian basis

**What GeoVac already has:**
- The angular Legendre expansion (which is the angular part of the Kereselidze construction)
- The separation constant solver

**What would be new:**
- The radial Sturmian functions defined via confluent Heun equation solutions
- Closed-form matrix element evaluation using Heun function properties
- Analytical hybridization coefficients connecting spheroidal and spherical bases

### 1.3 Other Relevant Work

**Singor, Savage, Bray, Schneider & Fursa (2023):** "Continuum solutions to the two-center Coulomb problem in prolate spheroidal coordinates," Computer Physics Communications 282, 108514. Extends the spheroidal harmonic expansion approach to positive-energy continuum states. Uses numerical methods for the quasi-radial equation but validates against known phase shifts.

**Leaver (various):** Solutions to the generalized spheroidal wave equation using continued fraction methods. These provide an alternative to direct diagonalization for computing separation constants and could accelerate the angular solver.

**General observation:** The literature confirms that the angular part of the Level 2 problem is well-handled by spectral methods (Legendre expansion), and the radial part is the primary target for algebraicization. The Mitnik et al. GSF approach and the Kereselidze Heun equation approach are complementary routes to the same goal: replacing the FD radial grid with a spectral radial basis.

### 1.4 Critical Distinction from Failed Sturmian Approaches

**CLAUDE.md Section 3 records that shared-p0 Sturmians fail for molecules** (Papers 8-9, Structural Theorem). This failure is fundamentally different from what the literature papers propose:

- **GeoVac's failed approach (Papers 8-9):** Used Sturmian functions on a *single S3* with a *shared momentum scale p0* for both atomic centers. The Hamiltonian becomes proportional to the overlap matrix (H ~ S), making eigenvalues R-independent. This is an inherent limitation of the single-sphere topology.

- **Mitnik et al. / Kereselidze et al.:** Use Sturmian functions in *prolate spheroidal coordinates*, where the two-center geometry is built into the coordinate system itself. The separation parameter c^2 = -R^2 E/2 provides the R-dependence that the single-S3 approach lacks. There is no shared p0 -- the energy enters through c^2, and the Sturmian scaling parameter is adapted to the radial equation at each energy.

The literature approaches are spectral bases for the *already-separated* equations in prolate spheroidal coordinates -- they replace the FD grid in the radial equation, not the coordinate system. This is directly analogous to how the Level 3 algebraic angular solver (Track B, v2.0.6) replaced the FD alpha grid with a Gegenbauer spectral basis while keeping the hyperspherical coordinate system.

---

## 2. Complete Matrix Element Inventory

The Level 2 solver consists of two coupled components: an angular solver (spectral) and a radial solver (FD), linked by Brent root-finding in c^2.

### Architecture Overview

```
ProlateSpheroidalLattice.solve():
  Input: R, Z_A, Z_B, N_xi, xi_max, m, n_angular, n_radial

  Loop (Brent root-finding in c^2):
    1. Angular solver: A(c^2)          [SPECTRAL - Legendre basis]
    2. Radial solver: lambda_max(c^2,A) [FD - 5000-point grid]
    3. Residual: f(c^2) = lambda_max

  Output: E_elec = -2*c2/R^2
```

### 2.1 Angular Equation: Spectral Solver

**Equation (Paper 11, Eq. 6):**
$$\frac{d}{d\eta}\left[(1-\eta^2)\frac{dG}{d\eta}\right] + \left(-A + c^2\eta^2 + b\eta - \frac{m^2}{1-\eta^2}\right)G = 0$$

**Implementation:** `molecular_sturmian.py:_angular_sep_const()` (lines 34-64)

Expands G(eta) in associated Legendre polynomials P_r^m(eta), r = m, m+1, ..., m+N-1, producing a dense eigenvalue problem of dimension N (default N=50).

#### 2.1.1 Diagonal: -r(r+1)

**Formula:** The free Legendre eigenvalue for P_r^m: -r(r+1).

**Implementation:** `H = np.diag(-r_vals * (r_vals + 1))` (line 59)

**Algebraic status:** :white_check_mark: **Already algebraic.** Pure quantum number expression. These are the eigenvalues of the Legendre operator d/deta[(1-eta^2)d/deta] - m^2/(1-eta^2).

#### 2.1.2 eta^2 coupling: c^2 * nu^2

**Formula:** The matrix elements of eta^2 in the orthonormal Legendre basis. Uses the three-term recurrence eta * P_r^m = alpha_r P_{r+1}^m + beta_r P_{r-1}^m, so eta^2 = (eta-matrix)^2, which is pentadiagonal.

**Implementation:** `nu_mat` is constructed (lines 52-57) as the tridiagonal matrix representation of eta in the orthonormal Legendre basis. Then `H += c**2 * (nu_mat @ nu_mat)` (line 60) gives the pentadiagonal eta^2 contribution.

The nu_mat elements are:
$$\nu_{r,r+1} = \frac{(r - m + 1)}{(2r+1)} \sqrt{\frac{\text{norm}_{r+1}}{\text{norm}_r}}$$
where norm_r = 2/(2r+1) * (r+m)!/(r-m)!.

**Algebraic status:** :white_check_mark: **Already algebraic.** The recurrence coefficients are closed-form expressions of quantum numbers (r, m) and factorials. The matrix-matrix product nu^2 is exact.

#### 2.1.3 Heteronuclear coupling: b * nu

**Formula:** The matrix elements of b*eta in the Legendre basis, where b = R(Z_B - Z_A).

**Implementation:** `H += b * nu_mat` (line 61). Uses the same tridiagonal nu_mat.

**Algebraic status:** :white_check_mark: **Already algebraic.** Same recurrence coefficients as 2.1.2.

#### 2.1.4 Norm coefficients

**Formula:** norm_r = 2/(2r+1) * (r+m)!/(r-m)!

**Implementation:** Lines 47-49, computed via `factorial()`.

**Algebraic status:** :white_check_mark: **Already algebraic.** Factorials of integers.

#### 2.1.5 Diagonalization

**Implementation:** `np.linalg.eigvalsh(H)` (line 63), returning the (N-1-n_sph)-th eigenvalue from the sorted spectrum.

**Algebraic status:** :x: **Numerically required.** The eigenvalue A(c^2) at finite c^2 is transcendental -- it is not a closed-form function of c^2. However, the matrix being diagonalized is small (50x50) and has algebraic entries. The diagonalization is O(N^3) = O(50^3) ~ 10^5 ops, which is negligible.

**Summary: The angular solver is essentially 100% algebraic.** Matrix elements are closed-form; only the final eigenvalue extraction is numerical (small dense diag).

---

### 2.2 Radial Equation: Finite-Difference Solver

**Equation (Paper 11, Eq. 5):**
$$\frac{d}{d\xi}\left[(\xi^2-1)\frac{dF}{d\xi}\right] + \left(A + a\xi - c^2\xi^2 - \frac{m^2}{\xi^2-1}\right)F = 0$$

where a = R(Z_A + Z_B).

**Implementation:** `prolate_spheroidal_lattice.py:_radial_top_eigenvalue()` (lines 89-125)

Uses a uniform FD grid in xi with N=5000 points, xi_min = 1.0005, xi_max = 25.0. Self-adjoint stencil with half-point evaluation of p(xi) = xi^2 - 1.

#### 2.2.1 Kinetic operator: d/dxi[(xi^2-1) dF/dxi]

**Formula (Paper 11, Eq. 10-11):** Self-adjoint FD stencil:
$$\frac{1}{h^2}\left[p_{i+1/2}(F_{i+1} - F_i) - p_{i-1/2}(F_i - F_{i-1})\right]$$
where p_{i+1/2} = (xi_i + h/2)^2 - 1.

**Implementation:** Lines 103-113.
- `p_plus = (xi + h/2)**2 - 1`
- `p_minus = (xi - h/2)**2 - 1`
- `diag = -(p_plus + p_minus)/h^2 + q`
- `off = p_plus[:-1]/h^2`

**Algebraic status:** :x: **Irreducibly numerical in the current formulation.** The FD stencil discretizes a continuous differential operator on a grid. However, if the radial solution were expanded in a spectral basis (e.g., Coulomb Sturmians S_n(xi; alpha) with built-in exp(-alpha*xi) decay and (xi^2-1)^{m/2} behavior at xi=1), the kinetic matrix elements between basis functions would be:

$$K_{nn'} = \int_1^\infty S_n(\xi) \frac{d}{d\xi}\left[(\xi^2-1)\frac{dS_{n'}}{d\xi}\right] d\xi$$

For Sturmian bases derived from the confluent Heun equation (Kereselidze et al.), or for generalized Laguerre functions in a mapped coordinate (Mitnik et al.), these integrals evaluate via recurrence relations or closed-form expressions. **Algebraic-pending.**

#### 2.2.2 Nuclear attraction: a * xi

**Formula:** The linear potential a*xi where a = R(Z_A + Z_B).

**Implementation:** Part of `q = A + self._a * xi - c2 * xi**2` (line 104). Evaluated pointwise on the FD grid.

**Algebraic status:** :large_orange_diamond: **Algebraic in principle.** In a spectral radial basis {phi_n(xi)}, the matrix element is:

$$\langle \phi_n | \xi | \phi_{n'} \rangle = \int_1^\infty \phi_n(\xi) \cdot \xi \cdot \phi_{n'}(\xi) \cdot w(\xi) \, d\xi$$

For polynomial or Sturmian bases, this is a three-term recurrence (xi * phi_n = alpha_n phi_{n+1} + beta_n phi_n + gamma_n phi_{n-1}), making the matrix tridiagonal with algebraic entries. **Known closed form for Laguerre-type bases.**

#### 2.2.3 Quadratic confinement: -c^2 * xi^2

**Formula:** The confining potential -c^2 * xi^2.

**Implementation:** Part of `q` (line 104). Evaluated pointwise.

**Algebraic status:** :large_orange_diamond: **Algebraic in principle.** Same as 2.2.2 but with xi^2. In a spectral basis, xi^2 = (xi-matrix)^2, which is pentadiagonal. Matrix elements are algebraic via the recurrence relation.

#### 2.2.4 Centrifugal barrier: -m^2/(xi^2-1)

**Formula:** The centrifugal term for m != 0 states.

**Implementation:** Lines 105-106: `q -= self.m**2 / xi2_1`.

**Algebraic status:** :large_orange_diamond: **Algebraic in principle.** For Sturmian bases with the built-in factor (xi^2-1)^{m/2}, the centrifugal matrix elements involve integrals of (xi^2-1)^{-1} against the basis functions. These are related to digamma functions or partial fraction integrals. For polynomial bases in a mapped coordinate (e.g., t = xi - 1), these evaluate to finite sums. The Mitnik et al. approach handles this by building the (xi^2-1)^{m/2} behavior into the basis functions themselves, so the centrifugal term becomes part of the basis definition rather than a matrix element.

#### 2.2.5 Separation constant: A (diagonal shift)

**Formula:** The constant A from the angular solver, added to the diagonal.

**Implementation:** Part of `q = A + ...` (line 104).

**Algebraic status:** :white_check_mark: **Already algebraic** (as a matrix element -- it is just a scalar shift). The value of A itself is transcendental (from angular diagonalization), but its contribution to the radial matrix is trivially A * I (identity).

#### 2.2.6 Neumann boundary condition at xi=1

**Formula:** For m=0, the flux p(xi) dF/dxi -> 0 as xi -> 1 (natural BC since p(1) = 0).

**Implementation:** Lines 115-116: `diag[0] = -p_plus[0]/h^2 + q[0]` (sets p_minus = 0 at first grid point).

**Algebraic status:** :large_orange_diamond: **Handled differently in spectral basis.** In a Sturmian/spectral approach, the regularity condition at xi=1 is built into the basis functions (they are chosen to be regular at the singular point). No explicit boundary condition imposition is needed -- the basis automatically satisfies it. This is one of the key advantages of the spectral approach (Mitnik et al.: "no derivatives have to be calculated numerically" because the GSF "obey appropriate physical boundary conditions").

#### 2.2.7 Dirichlet BC at xi_max

**Formula:** F(xi_max) = 0.

**Implementation:** Implicit in the FD grid (F is not defined beyond xi_max).

**Algebraic status:** :large_orange_diamond: **Handled differently in spectral basis.** Sturmian basis functions have built-in exponential decay exp(-alpha*xi), automatically satisfying the asymptotic BC. The decay parameter alpha is determined by c^2 (the energy). No artificial truncation domain needed.

#### 2.2.8 Tridiagonal diagonalization

**Implementation:** `eigh_tridiagonal(diag, off, eigvals_only=True)` (line 118). Returns all N=5000 eigenvalues; the top eigenvalue (or n_radial-th from top) is extracted.

**Algebraic status:** :x: **Numerically required** (in any formulation). Even with algebraic matrix elements, the eigenvalue of the radial matrix at finite c^2 is transcendental. However, the matrix dimension would shrink from N=5000 (FD) to N~10-20 (spectral basis), reducing cost from O(5000) to O(20^3) per evaluation.

---

### 2.3 Root-Finding: Self-Consistency in c^2

**Equation:** Find c^2 such that lambda_max(c^2, A(c^2)) = 0.

**Implementation:** `prolate_spheroidal_lattice.py:solve()` (lines 143-186). Uses `scipy.optimize.brentq` with xtol=1e-12. Typically ~50 iterations.

**Algebraic status:** :x: **Irreducibly numerical.** The self-consistency condition requires evaluating the angular and radial problems at multiple c^2 values. This is a transcendental equation in c^2.

**Note:** Each Brent iteration calls both the angular solver (cheap: 50x50 diag) and the radial solver (expensive: 5000-point tridiag). The radial solver dominates the cost. Replacing the radial FD with a spectral basis would reduce each iteration from O(5000) to O(20^3) ~ O(8000), roughly comparable but with much higher accuracy.

---

### 2.4 Wavefunction Evaluation

**Implementation:** `solve_with_wavefunction()` (lines 193-284). Reconstructs F(xi) and G(eta) after solving.

- **G(eta):** Evaluated from Legendre coefficients using `scipy.special.lpmv`. Algebraic (just basis function evaluation).
- **F(xi):** Extracted as FD eigenvector, then interpolated. This would become algebraic if the radial basis were spectral (just evaluate basis functions at desired points).

---

### 2.5 PES Scan and Spectroscopic Constants

**Implementation:** `scan_h2plus_pes()` (lines 287-330) and `fit_spectroscopic_constants()` (lines 333-364).

- PES scan: loops over R values, calling `solve()` at each. Independent evaluations.
- Spectroscopic fit: quadratic polynomial fit near minimum. Algebraic (polyfit).

---

## 3. Coupling Integral Census

Unlike the Level 3 solver, the Level 2 solver has no inter-channel coupling (it is a single-electron, fully separated problem). The "coupling" is between the angular and radial equations through the shared parameters (c^2, A), enforced by root-finding rather than matrix coupling.

### 3a. Angular matrix elements (all algebraic)

| Element | Formula | Size | Status |
|:--------|:--------|:-----|:------:|
| Legendre eigenvalue | -r(r+1) | N_ang diagonal | :white_check_mark: Algebraic |
| eta recurrence (nu) | Factorials of (r, m) | N_ang tridiagonal | :white_check_mark: Algebraic |
| eta^2 = nu^2 | Matrix product | N_ang pentadiagonal | :white_check_mark: Algebraic |
| Norm coefficients | 2/(2r+1) * (r+m)!/(r-m)! | N_ang | :white_check_mark: Algebraic |

### 3b. Radial matrix elements (all FD-numerical)

| Element | Formula | Size | Status |
|:--------|:--------|:-----|:------:|
| Kinetic d/dxi[p dF/dxi] | SA FD stencil, p = xi^2-1 | N_xi tridiagonal | :large_orange_diamond: Algebraic-pending |
| Nuclear attraction a*xi | Pointwise on grid | N_xi diagonal | :large_orange_diamond: Algebraic-pending |
| Quadratic -c^2*xi^2 | Pointwise on grid | N_xi diagonal | :large_orange_diamond: Algebraic-pending |
| Centrifugal -m^2/(xi^2-1) | Pointwise on grid | N_xi diagonal | :large_orange_diamond: Algebraic-pending |
| Separation constant A | Scalar shift | N_xi diagonal | :white_check_mark: Algebraic |
| Neumann BC at xi=1 | Modified first row | 1 element | :large_orange_diamond: Built into spectral basis |
| Dirichlet BC at xi_max | Implicit truncation | - | :large_orange_diamond: Built into spectral basis |

---

## 4. Feasibility Assessment

### 4.1 Current Algebraic Fraction

| Component | Method | Status | Fraction of total cost |
|:----------|:-------|:------:|:----------------------:|
| Angular matrix construction | Algebraic (Legendre recurrence) | :white_check_mark: | ~1% |
| Angular diagonalization | Numerical (50x50 dense) | :x: Required | ~1% |
| Radial matrix construction | FD stencil (5000 points) | :large_orange_diamond: Pending | ~8% |
| Radial diagonalization | Numerical (5000 tridiag) | :x: Required | ~40% |
| Root-finding iterations (~50x) | Numerical (Brent) | :x: Required | multiplier |

The radial FD solver dominates cost: ~50 Brent iterations x (angular diag + radial diag) per R point. The angular part is cheap (~1ms per eval). The radial tridiag diag at N=5000 is ~10ms per eval, so ~500ms per energy point.

### 4.2 What Could Be Made Algebraic

**Target: Replace the N=5000 FD radial grid with an N~10-20 spectral radial basis.**

Two candidate spectral bases for the xi equation:

#### Option A: Mapped Laguerre / Generalized Sturmian (Mitnik et al.)

Define radial basis functions:
$$\phi_n(\xi) = (\xi^2 - 1)^{m/2} \cdot e^{-\alpha(\xi-1)} \cdot L_n^{(2\alpha)}(2\alpha(\xi-1))$$

where L_n are generalized Laguerre polynomials, alpha = sqrt(c^2) (the asymptotic decay rate), and the prefactor handles the xi=1 singularity.

Matrix elements in this basis:
- **Kinetic:** Tridiagonal (Laguerre three-term recurrence under differentiation)
- **a*xi = a*(1 + (xi-1)):** Tridiagonal (shift + Laguerre recurrence for (xi-1))
- **-c^2*xi^2:** Pentadiagonal (Laguerre recurrence squared)
- **-m^2/(xi^2-1):** Requires evaluation, but the (xi^2-1)^{m/2} prefactor cancels the singularity, leaving a smooth integral

All matrix elements are **finite sums of Gamma function ratios** -- algebraic.

#### Option B: Coulomb Spheroidal Sturmians (Kereselidze et al.)

Define basis functions as polynomial solutions of the confluent Heun equation. These are the "natural" eigenfunctions of the auxiliary Sturmian problem with a Coulomb generating potential in spheroidal coordinates.

Advantages: Basis functions are exact eigenfunctions of a closely related operator, so convergence is exponential. Closed-form expressions exist.

Disadvantages: More complex to implement; Heun function recurrence relations are less standard than Laguerre.

#### Option C: Chebyshev/mapped polynomial (simplest)

Map xi in [1, inf) to t in [-1, 1] via t = 1 - 2/(xi - 1 + 1) or similar algebraic map. Expand F in Chebyshev polynomials T_n(t). Matrix elements computed via Gauss-Chebyshev quadrature (high-order, exponentially convergent).

This is the simplest to implement but matrix elements are not strictly algebraic (they use quadrature). However, with ~30-50 quadrature points, the accuracy would far exceed the current N=5000 FD grid.

### 4.3 The Irreducible Numerical Core

Three components remain irreducibly numerical regardless of algebraicization:

1. **Angular eigenvalue A(c^2):** Transcendental function of c^2. However, computed from a small (50x50) algebraic matrix, so this is already efficient.

2. **Radial eigenvalue lambda(c^2, A):** Transcendental. With a spectral basis, computed from a small (~15x15) matrix with algebraic entries.

3. **Root-finding in c^2:** Brent iteration over the transcendental residual. Cannot be eliminated. But each iteration becomes much cheaper with a spectral radial solver.

### 4.4 Expected Performance Improvements

#### Accuracy

The current solver achieves 0.70% error for H2+ at N_xi=5000 (Paper 11, Table III). The convergence is O(1/N_xi) -- suboptimal for a second-order FD scheme, attributed to:
1. Grid offset delta=5e-4 from the xi=1 singularity
2. First-order Neumann BC implementation
3. Nonuniform behavior of p(xi) = xi^2-1 near xi=1

Paper 11, Appendix A explicitly notes: "A spectral discretization of the radial equation (e.g., a mapped Chebyshev method with coordinate transformation xi = 1 + t^2 to resolve the xi=1 endpoint) would eliminate these issues and achieve exponential convergence."

**With a spectral radial basis:**
- Exponential convergence in basis size (vs O(1/N) for FD)
- N~15 basis functions should exceed N=5000 FD accuracy
- The xi=1 singularity is handled exactly by the basis (built-in (xi^2-1)^{m/2} factor)
- Asymptotic decay handled exactly (built-in exp(-alpha*xi))
- Expected: sub-0.01% accuracy with N_radial ~ 15, vs 0.70% currently

**Comparison with literature:**
- Kereselidze et al. achieve E(1sigma_g) = -1.102634186 a.u. with 6 Legendre polynomials (angular), implying high radial accuracy. The exact value is E = -1.1026342144949 a.u. (Bates et al.), so their result is accurate to ~10^-7 a.u. (~10^-5 %).
- Mitnik et al. report "very accurate results with minimal basis sets" (5-15 radial functions), consistent with exponential convergence.

#### Speed

| Operation | Current (FD) | Spectral (est.) | Speedup |
|:----------|:-------------|:-----------------|:--------|
| Radial matrix build | O(N_xi) = O(5000) | O(N_basis^2) = O(225) | ~20x |
| Radial eigensolve | O(N_xi) = O(5000) tridiag | O(N_basis^3) = O(3375) dense | ~1.5x |
| Per Brent iteration | ~10ms | ~0.5ms | ~20x |
| Full solve (50 iters) | ~500ms | ~25ms | ~20x |
| PES scan (50 R points) | ~25s | ~1.3s | ~20x |

The dominant cost shifts from the radial eigensolve (currently O(N_xi) tridiag) to the root-finding iteration count. With a more accurate solver, the Brent convergence tolerance could be tightened without significant cost increase.

#### Matrix Dimension Reduction

| Component | Current | Spectral | Reduction |
|:----------|:--------|:---------|:---------:|
| Angular basis | 50 | 50 (unchanged) | 1x |
| Radial grid/basis | 5000 | 10-20 | 250-500x |
| Total matrix dimension | 5000 | 10-20 | 250-500x |

---

## 5. Comparison with Level 3 Algebraicization (Track B)

The Level 3 algebraic angular solver (v2.0.6, `geovac/algebraic_angular.py`) provides a direct template for the Level 2 radial algebraicization:

| Aspect | Level 3 (He, Track B) | Level 2 (H2+, this audit) |
|:-------|:---------------------|:--------------------------|
| Target equation | Angular: H_ang(alpha) at fixed R | Radial: L_xi at fixed c^2, A |
| Current method | FD grid (200-800 points) | FD grid (5000 points) |
| Proposed spectral basis | Gegenbauer C_n^2(cos 2alpha) | Laguerre/Sturmian in xi |
| Dimension reduction | 200-800 -> 10-30 | 5000 -> 10-20 |
| Singularity handling | alpha=0, pi/2 (centrifugal) | xi=1 (coordinate singularity) |
| Key matrix elements | Nuclear coupling, V_ee coupling | a*xi, c^2*xi^2, m^2/(xi^2-1) |
| Algebraic status of elements | Proven algebraic (Paper 13 Sec XII) | Known algebraic (Laguerre recurrence) |
| Literature precedent | Krivec & Mandelzweig (variational Gegenbauer) | Mitnik et al. 2021, Kereselidze et al. 2016 |

The Level 2 case is actually **simpler** than Level 3 in several ways:
- No inter-channel coupling (single electron, fully separated)
- The radial potential is a low-degree polynomial in xi (linear + quadratic + inverse), vs the Level 3 angular potential which involves 1/cos(alpha) and 1/sin(alpha)
- Standard spectral basis functions (Laguerre polynomials) exist for the semi-infinite interval [1, inf), vs the finite interval [0, pi/2] for Level 3

---

## 6. Summary of Algebraic Status

| Matrix Element Type | Current Method | Algebraic Status | Replacement |
|:--------------------|:---------------|:-----------------|:------------|
| Angular Legendre eigenvalues -r(r+1) | Closed form | :white_check_mark: Already algebraic | None needed |
| Angular eta recurrence (nu_mat) | Factorial formula | :white_check_mark: Already algebraic | None needed |
| Angular eta^2 = nu^2 | Matrix product | :white_check_mark: Already algebraic | None needed |
| Angular norm coefficients | Factorial formula | :white_check_mark: Already algebraic | None needed |
| Heteronuclear b*eta | Tridiag from nu | :white_check_mark: Already algebraic | None needed |
| Angular eigenvalue A(c^2) | Dense diag (50x50) | :x: Transcendental | Keep (cheap) |
| Radial kinetic d/dxi[p dF/dxi] | FD stencil (N=5000) | :large_orange_diamond: Algebraic-pending | Laguerre/Sturmian recurrence |
| Radial nuclear attraction a*xi | FD pointwise | :large_orange_diamond: Algebraic-pending | Laguerre three-term recurrence |
| Radial confinement -c^2*xi^2 | FD pointwise | :large_orange_diamond: Algebraic-pending | Laguerre five-term (recurrence^2) |
| Radial centrifugal -m^2/(xi^2-1) | FD pointwise | :large_orange_diamond: Algebraic-pending | Built into basis prefactor |
| Separation constant A (shift) | Scalar | :white_check_mark: Already algebraic | None needed |
| Neumann BC at xi=1 | Modified FD row | :large_orange_diamond: Algebraic-pending | Built into basis |
| Asymptotic BC at xi_max | FD truncation | :large_orange_diamond: Algebraic-pending | Built into basis (exp decay) |
| Radial eigenvalue lambda(c^2, A) | Tridiag diag (N=5000) | :x: Transcendental | Dense diag (~15x15) |
| Self-consistency c^2 | Brent root-finding | :x: Irreducibly numerical | Keep (fewer iters needed) |

**Bottom line:** The angular equation is already 100% algebraic (spectral Legendre basis with closed-form matrix elements). The radial equation is 100% numerical (FD grid), but every matrix element has a known algebraic replacement via Laguerre/Sturmian recurrence relations. The algebraicization target is clear and well-supported by literature.

---

## 7. Recommended Implementation Path

### Priority 1: Spectral Radial Solver (High Impact, Medium Difficulty)

**Goal:** Replace the N=5000 FD radial grid with an N~15 spectral basis.

**Recommended basis:** Mapped Laguerre polynomials (Option A from Section 4.2).

**Steps:**

1. **Define the mapped coordinate and basis functions.**
   - Map: t = xi - 1, t in [0, inf)
   - Basis: phi_n(xi) = (xi^2-1)^{m/2} * exp(-alpha*t) * L_n^{(k)}(2*alpha*t)
   - where alpha = sqrt(c^2), k = m+1 (to match singular behavior)

2. **Derive radial matrix elements algebraically.**
   - Overlap: <phi_n|phi_n'> via Laguerre orthogonality (diagonal if basis chosen correctly)
   - Kinetic: d/dxi[p dF/dxi] via integration by parts + Laguerre recurrence
   - a*xi = a*(1+t): tridiagonal (Laguerre recurrence for t)
   - c^2*xi^2 = c^2*(1+t)^2: pentadiagonal
   - Verify against FD results at multiple c^2 values

3. **Build `SpectralRadialSolver` class.**
   - Input: c^2, A, m, n_basis, Z_A, Z_B, R
   - Construct algebraic matrix H_rad (n_basis x n_basis)
   - Diagonalize, return top eigenvalue (and eigenvector for wavefunction)

4. **Integrate into `ProlateSpheroidalLattice`.**
   - Replace `_radial_top_eigenvalue()` with spectral version
   - Keep Brent root-finding unchanged
   - Validate: same energies as FD (to FD accuracy), then demonstrate superior convergence

5. **Benchmark.**
   - Convergence: energy vs n_basis at fixed R
   - Accuracy: compare to exact H2+ at multiple R
   - Speed: wall time comparison

**Expected outcome:** Sub-0.01% accuracy with N_basis~15, ~20x speedup, elimination of all FD discretization artifacts (grid offset, boundary error, sub-optimal convergence rate).

**Estimated effort:** 2-3 sub-agent sessions. The mathematics is well-established (Laguerre recurrences are textbook). The main work is deriving the specific matrix elements for the prolate spheroidal radial operator and validating against the existing FD solver.

### Priority 2: Wavefunction Quality (Medium Impact, Low Difficulty)

Once the spectral radial solver is in place, the radial wavefunction F(xi) is available as a linear combination of known basis functions -- no interpolation needed. This improves:
- V_ee quadrature for two-electron extensions (Paper 11, Section VIII.D)
- Density evaluation for composed geometry (Paper 17)
- Visualization

### Priority 3: Perturbation Series for A(c^2) (Low Impact, Medium Difficulty)

The angular separation constant A(c^2) could be expanded as a power series in c^2:
$$A(c^2) = A_0 + A_1 c^2 + A_2 c^4 + \ldots$$
where A_0 = -l(l+1) and A_k are algebraic (computable from Legendre recurrence). This would eliminate the angular diagonalization entirely, replacing it with a polynomial evaluation. However, the angular solver is already fast (~1ms), so the practical impact is minimal. A Pade approximant in c^2 could extend the useful range.

### Priority 4: Continued Fraction for Root-Finding (Low Impact, High Difficulty)

The separation constant A(c^2) and the radial eigenvalue lambda(c^2, A) satisfy a joint implicit equation. For the spheroidal wave equation (b=0), continued fraction representations exist (Flammer 1957, Leaver). These could potentially replace the Brent root-finding with a direct algebraic solution. However, the heteronuclear case (b != 0) complicates this, and the root-finding is already reliable. Low priority.

---

## 8. Risk Assessment

| Risk | Likelihood | Mitigation |
|:-----|:-----------|:-----------|
| Laguerre recurrence for kinetic term is unstable | Low | Use orthonormal basis; validate against FD at each step |
| Spectral basis converges slowly for large c^2 (deep binding) | Medium | Adapt alpha parameter to c^2; use 2-3 extra basis functions |
| Heteronuclear case (b != 0) breaks simplifications | Low | The b-dependence is entirely in the angular equation, which is already spectral |
| Implementation complexity exceeds benefit | Low | The FD -> spectral replacement is a well-trodden path (cf. Level 3 Track B) |
| Confused with failed Sturmian approach (Papers 8-9) | Medium | Document distinction clearly (Section 1.4 above). The spectral radial basis operates within the prolate spheroidal separation, not on a single S3 |

---

## 9. Connection to Broader GeoVac Algebraicization Program

The Level 2 radial algebraicization would extend the algebraic registry (CLAUDE.md Section 12) as follows:

### Level 2 (Prolate Spheroidal -- H2+)

| Matrix Element | Current Status | Target Status |
|:---------------|:--------------:|:-------------:|
| Angular Legendre eigenvalues | algebraic | algebraic |
| Angular eta/eta^2 recurrence | algebraic | algebraic |
| Angular eigenvalue A(c^2) | numerical (50x50 diag) | numerical (50x50 diag) |
| Radial kinetic operator | numerical (FD, N=5000) | algebraic (Laguerre recurrence) |
| Radial nuclear attraction a*xi | numerical (FD pointwise) | algebraic (Laguerre recurrence) |
| Radial confinement c^2*xi^2 | numerical (FD pointwise) | algebraic (Laguerre recurrence) |
| Radial centrifugal m^2/(xi^2-1) | numerical (FD pointwise) | algebraic (built into basis) |
| Radial BCs (xi=1, xi_max) | numerical (FD boundary) | algebraic (built into basis) |
| Radial eigenvalue lambda | numerical (N=5000 tridiag) | numerical (~15x15 dense diag) |
| Self-consistency c^2 | numerical (Brent) | numerical (Brent) |

This would make Level 2 the second level (after Level 1, which is fully algebraic on S3) where all matrix elements are algebraic, with only eigenvalue extraction and self-consistency remaining numerical. The pattern is consistent across levels: algebraic matrix elements + numerical small-matrix diagonalization.

---

## 10. Literature References

1. Mitnik, D.M., Lopez, F.A. & Ancarani, L.U. (2021). "Generalized Sturmian Functions in prolate spheroidal coordinates." Molecular Physics 119(8), e1881179. [arXiv:2006.06616](https://arxiv.org/abs/2006.06616)

2. Kereselidze, T., Chkadua, G., Defrance, P. & Ogilvie, J.F. (2016). "Derivation, properties and application of Coulomb Sturmians defined in spheroidal coordinates." Molecular Physics 114(1), 148-161. [PDF](http://www.cecm.sfu.ca/personal/ogilvie/jfopub/16Sturmi.pdf)

3. Kereselidze, T. & Ogilvie, J.F. (2018). "The Hydrogen-Atom Problem and Coulomb Sturmian Functions in Spheroidal Coordinates." Advances in Quantum Chemistry. [ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S0065327618300030)

4. Singor, A., Savage, J.S., Bray, I., Schneider, B.I. & Fursa, D.V. (2023). "Continuum solutions to the two-center Coulomb problem in prolate spheroidal coordinates." Computer Physics Communications 282, 108514.

5. Bates, D.R., Ledsham, K. & Stewart, A.L. (1953). "Wave Functions of the Hydrogen Molecular Ion." Phil. Trans. R. Soc. A 246, 215-240.

6. Flammer, C. (1957). Spheroidal Wave Functions. Stanford University Press.
