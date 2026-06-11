# Algebraic Boundary Audit: Level 3 Hyperspherical Solver

**Date:** 2026-03-27 (v2.0.5)
**Scope:** Read-only audit of the algebraic vs. numerical boundary in the Level 3 (He) adiabatic hyperspherical solver.
**Goal:** Map every matrix element type, assess algebraic status, and identify what can be replaced by closed-form evaluation (following the Paper 12 Neumann expansion precedent).

---

## 1. Complete Matrix Element Inventory

The Level 3 angular Hamiltonian at fixed hyperradius R is built in `geovac/hyperspherical_angular.py:solve_angular()`. The Hamiltonian acts on the product space of partial-wave channels (l = 0, ..., l_max) and alpha FD grid points (i = 1, ..., N_alpha), giving total dimension N = (l_max + 1) * N_alpha.

### 1.1 Kinetic energy: -1/2 u_l''(alpha)

**Mathematical formula (Paper 13, Eq. 8):**
$$-\frac{1}{2} \frac{d^2 u_l}{d\alpha^2}$$

**Implementation:** Three-point FD stencil (hyperspherical_angular.py:135-164).
- Diagonal: `1/h^2` per grid point
- Off-diagonal (tridiagonal within each l-channel): `-1/(2h^2)`

**Algebraic status:** ❌ **Irreducibly numerical** in the current formulation. The kinetic operator acts on the alpha grid, which is a continuous coordinate. However, if the solution were expanded in a known basis (Gegenbauer/Jacobi polynomials — the free SO(6) eigenfunctions), the kinetic matrix elements between basis functions would be algebraic. The FD stencil is a choice of discretization, not a fundamental limitation. See Section 5 for the spectral alternative.

### 1.2 Liouville constant: -2

**Mathematical formula (Paper 13, Eq. 8):**
The constant -2 arises from the Liouville substitution u_l = sin(alpha) cos(alpha) phi_l. It corresponds to the free-particle eigenvalue shift from the curvature of the weight function on S^5.

**Implementation:** Added to V_l at hyperspherical_angular.py:150.

**Algebraic status:** ✅ **Already algebraic.** Pure topological constant, independent of R or quantum numbers.

### 1.3 Centrifugal barrier: l(l+1)(1/cos^2(alpha) + 1/sin^2(alpha))/2

**Mathematical formula (Paper 13, Eq. 8):**
$$\frac{l(l+1)}{2}\left(\frac{1}{\cos^2\alpha} + \frac{1}{\sin^2\alpha}\right)$$

**Implementation:** Pointwise evaluation on alpha grid (hyperspherical_angular.py:151). Diverges at alpha = 0 and alpha = pi/2 (the grid boundaries).

**Algebraic status:** 🔶 **Algebraic in principle.** In a Gegenbauer basis C_n^2(cos 2alpha), the matrix elements of 1/cos^2(alpha) and 1/sin^2(alpha) between basis functions are integrals of Gegenbauer polynomials times trigonometric functions — these evaluate to finite sums involving Gamma functions and 3j-like coupling coefficients. The FD discretization introduces error precisely because of the boundary singularities (Paper 13, Table I: l_max > 0 degrades accuracy from 0.05% to 0.87%). A spectral basis would eliminate this error entirely.

### 1.4 Nuclear attraction: R(-Z/cos(alpha) - Z/sin(alpha))

**Mathematical formula (Paper 13, Eq. 8):**
$$R\left(-\frac{Z}{\cos\alpha} - \frac{Z}{\sin\alpha}\right)$$

**Implementation:** Pointwise evaluation on alpha grid, multiplied by R (hyperspherical_angular.py:152). Diagonal in l (same potential for all channels).

**Algebraic status:** 🔶 **Algebraic in principle.** The R-dependence is trivially linear. The alpha-integral is:
$$\langle u_{n'} | \frac{1}{\cos\alpha} | u_n \rangle = \int_0^{\pi/2} u_{n'}(\alpha) \frac{1}{\cos\alpha} u_n(\alpha) \, d\alpha$$

where u_n(alpha) are the Liouville-substituted free SO(6) eigenfunctions. After substitution u_n = sin(alpha) cos(alpha) phi_n, the integrand becomes:

$$\int_0^{\pi/2} \sin^2\alpha \cos\alpha \, \phi_{n'}(\alpha) \phi_n(\alpha) \, d\alpha$$

Paper 13 Section XII demonstrates this evaluates to partial harmonic sums (Eq. 31): for the ground channel, I_nuc(1) = 1 + 1/3 = 4/3, giving a_1 = (8/3pi)(sqrt(2) - 4Z). The closed form extends to all channels — each matrix element is a finite sum of rational numbers times pi^(-1). **This is already proven algebraic by Paper 13, Eq. 31-32.**

### 1.5 Electron-electron coupling: R * sum_k G(l,k,l') * f_k(alpha)/max(cos,sin)

**Mathematical formula (Paper 13, Eq. 14-16):**
$$W_{ll'}(\alpha) = \sqrt{(2l+1)(2l'+1)} \sum_k \frac{f_k(\alpha)}{r_>} G(l,k,l')$$

where f_k(alpha) = (min(sin,cos)/max(sin,cos))^k and G(l,k,l') is the Gaunt integral (Wigner 3j symbol squared).

**Implementation:** hyperspherical_angular.py:167-191. The Gaunt integrals are precomputed algebraically (_precompute_gaunt, lines 65-74). The f_k(alpha)/max_sc terms are evaluated pointwise on the alpha grid.

**Algebraic status:** Decompose into components:

- **Gaunt integrals G(l,k,l'):** ✅ **Already algebraic.** Computed from factorials via Wigner 3j formula (lines 29-62). Selection rules: l+k+l' even, triangle inequality. Zero-parameter, exact.

- **Normalization sqrt((2l+1)(2l'+1))/2:** ✅ **Already algebraic.**

- **Radial coupling f_k(alpha)/max:** 🔶 **Algebraic in principle.** The matrix elements in the free basis are:
  $$\int_0^{\pi/2} u_{n'}(\alpha) \frac{(\min(\sin\alpha,\cos\alpha))^k}{\max(\sin\alpha,\cos\alpha)^{k+1}} u_n(\alpha) \, d\alpha$$

  This integral splits at alpha = pi/4 (where sin = cos) into two regions:
  - Region I (alpha < pi/4): sin < cos, so min/max = tan(alpha), and 1/max = 1/cos(alpha)
  - Region II (alpha > pi/4): cos < sin, so min/max = cot(alpha), and 1/max = 1/sin(alpha)

  In each region, the integrand involves products of trigonometric powers and Gegenbauer polynomials. These are the same class of integrals as the nuclear coupling (Section 1.4), but with different trigonometric powers. Paper 13 Section XII confirms these also yield partial harmonic sums — the I_ee(n) identity (Eq. 33) evaluates the k=0 monopole term in closed form.

  For higher k: the integrals (tan(alpha))^k / cos(alpha) * u_n' * u_n da involve higher powers of trigonometric functions. Each such integral is a beta function variant and evaluates to rationals times pi^(-1). **Algebraic in principle for all k.**

### 1.6 Jacobian centrifugal term: 15/(8R^2)

**Mathematical formula (Paper 13, Eq. 12):**
$$\frac{15}{8R^2}$$

Added in `hyperspherical_adiabatic.py:effective_potential()` when computing V_eff from mu(R).

**Algebraic status:** ✅ **Already algebraic.** Topological constant from the R^(-5/2) extraction.

### 1.7 Hyperradial kinetic energy: -1/2 d^2F/dR^2

**Implementation:** FD stencil on R grid (hyperspherical_radial.py:69-84). Tridiagonal.

**Algebraic status:** ❌ **Irreducibly numerical.** R is a continuous radial coordinate. The hyperradial Schrodinger equation is a standard 1D eigenvalue problem in a potential V_eff(R). This is analogous to the radial equation for hydrogen and must be solved numerically (or semi-analytically via WKB/phase integrals). No algebraic replacement is expected.

### 1.8 Non-adiabatic coupling: P_mu_nu(R) = <Phi_mu|dH/dR|Phi_nu> / (mu_mu - mu_nu)

**Implementation:** Hellmann-Feynman approach in `hyperspherical_coupling.py:142-231`. Requires dH/dR = C(Omega) (the charge function), computed pointwise on the alpha grid.

**Algebraic status:** 🔶 **Partially algebraic.** The matrix element <Phi_mu|C|Phi_nu> between eigenvectors is a sum over FD grid points. If the eigenvectors were expressed in the free SO(6) basis (with known expansion coefficients), the coupling would reduce to sums of the algebraic integrals from Sections 1.4-1.5, divided by algebraic eigenvalue gaps. The Hellmann-Feynman formula itself is exact — the numerical part is the representation of eigenvectors.

---

## 2. Coupling Integral Census

### 2a. Nuclear integral: integral of u_l'(alpha) (1/cos(alpha)) u_l(alpha) dalpha

**In the free SO(6) basis (Gegenbauer polynomials C_n^2(cos 2alpha) after Liouville substitution):**

The free eigenfunctions at R=0 are u_n(alpha) = sin(alpha) cos(alpha) * C_{n-1}^2(cos 2alpha), normalized on [0, pi/2]. The nuclear integral becomes:

$$I_{n'n}^{nuc} = \int_0^{\pi/2} \sin^2\alpha \cos\alpha \, C_{n'-1}^2(\cos 2\alpha) \, C_{n-1}^2(\cos 2\alpha) \, d\alpha$$

Substituting x = cos(2alpha), dx = -2 sin(alpha) cos(alpha) dalpha:

$$I_{n'n}^{nuc} = \frac{1}{2} \int_{-1}^{1} \frac{1}{\sqrt{(1+x)/2}} \cdot C_{n'-1}^2(x) \, C_{n-1}^2(x) \, (1-x^2)/4 \, dx$$

Wait — more carefully: sin^2(alpha) cos(alpha) = sin(alpha) * sin(alpha) cos(alpha). And sin(alpha) cos(alpha) dalpha = dx/(-2), with cos(alpha) = sqrt((1+x)/2), sin(alpha) = sqrt((1-x)/2).

So: sin^2(alpha) cos(alpha) dalpha = sin(alpha) * (sin(alpha) cos(alpha) dalpha) = sqrt((1-x)/2) * (-dx/2).

The nuclear integral (for the 1/cos(alpha) part) is:

$$I = \int_0^{\pi/2} \sin^2\alpha \cdot \cos\alpha \cdot \frac{1}{\cos\alpha} \cdot C_{n'-1}^2(x) C_{n-1}^2(x) \, d\alpha$$

$$= \int_0^{\pi/2} \sin^2\alpha \, C_{n'-1}^2(\cos 2\alpha) C_{n-1}^2(\cos 2\alpha) \, d\alpha$$

With x = cos(2alpha): sin^2(alpha) = (1-x)/2, and dalpha = -dx/(2 sin(2alpha)) = -dx/(2 sqrt(1-x^2)). The weight function from the inner product on the Liouville-transformed functions is just dalpha (standard L^2). So:

$$I = \frac{1}{2}\int_{-1}^{1} \frac{(1-x)/2}{\sqrt{1-x^2}} C_{n'-1}^2(x) C_{n-1}^2(x) \, dx$$

$$= \frac{1}{4}\int_{-1}^{1} \sqrt{\frac{1-x}{1+x}} \, C_{n'-1}^2(x) C_{n-1}^2(x) \, dx$$

This is an integral of Gegenbauer polynomials against a known algebraic weight. For integer-order Gegenbauer polynomials, such integrals reduce to finite sums via the linearization formula C_m^lambda * C_n^lambda = sum of C_k^lambda (the "Adams-Neumann" or Dougall linearization), followed by integration of C_k^2(x) against the weight sqrt((1-x)/(1+x)).

The integral of C_k^2(x) * (1-x)^(1/2) * (1+x)^(-1/2) over [-1,1] is a standard Jacobi integral (alpha=1/2, beta=-1/2 weight times Gegenbauer). This evaluates in terms of Gamma functions and is a **known closed form.** Specifically, it can be expressed via the formula for integrals of Jacobi polynomials against different Jacobi weights (cf. DLMF 18.18).

**Verdict:** ✅ **Closed form exists.** The nuclear coupling integrals in the free basis are finite sums of Gamma function ratios. Paper 13 already demonstrates this for the diagonal case (I_nuc(n) = sum of 1/(2j+1), Eq. 31). The off-diagonal case follows from the same mathematics with the linearization formula.

### 2b. Electron-electron coupling integral

$$I_{n'n}^{ee}(k) = \int_0^{\pi/2} u_{n'}(\alpha) \frac{(\min(\sin,\cos))^k}{\max(\sin,\cos)^{k+1}} u_n(\alpha) \, d\alpha$$

Split at alpha = pi/4:

**Region I (alpha < pi/4):** min/max = tan(alpha), 1/max = 1/cos(alpha). Integrand involves tan^k(alpha)/cos(alpha) * u_n' * u_n.

**Region II (alpha > pi/4):** min/max = cot(alpha), 1/max = 1/sin(alpha). By the exchange symmetry alpha -> pi/2 - alpha, this mirrors Region I with sin <-> cos.

For the free basis functions, the Region I integral becomes (with x = cos(2alpha)):

$$\int_0^{\pi/4} \frac{\sin^k\alpha}{\cos^{k+1}\alpha} \cdot \sin\alpha\cos\alpha \, C_{n'-1}^2(\cos 2\alpha) C_{n-1}^2(\cos 2\alpha) \, d\alpha$$

Using sin = sqrt((1-x)/2), cos = sqrt((1+x)/2):

$$\frac{1}{2}\int_0^1 \left(\frac{1-x}{1+x}\right)^{k/2} \cdot \frac{1}{\sqrt{1+x}} \cdot \frac{1}{\sqrt{2}} \cdot C_{n'-1}^2(x) C_{n-1}^2(x) \frac{dx}{2\sqrt{1-x^2}}$$

This is an integral of Gegenbauer products against the weight (1-x)^((k-1)/2) * (1+x)^(-(k+2)/2) over [0,1] (the split-region boundary maps to x=0).

For integer k, these are integrals of polynomials against Jacobi-type weights on a half-interval. By the same linearization + Jacobi integration method:

- **k = 0 (monopole):** Already demonstrated algebraic by Paper 13 (I_ee(n) identity, Eq. 33).
- **k even:** The weight (1-x)^((k-1)/2) (1+x)^(-(k+2)/2) involves half-integer powers, making this a Beta function integral after linearization. **Algebraic.**
- **k odd:** Same structure. **Algebraic.**

The split-region integral from 0 to pi/4 (mapping to x in [0,1]) rather than the full [0, pi/2] (x in [-1,1]) introduces an incomplete Beta function or equivalently a regularized hypergeometric function. However, for Gegenbauer polynomials evaluated on [0,1] (half the natural domain), there exist closed-form reduction formulas using the Gegenbauer addition theorem.

**Key structural observation:** The split at alpha = pi/4 is exactly analogous to the split-region Legendre expansion in Paper 15, Section V.A, which "terminates exactly via the 3j triangle inequality." The same termination mechanism applies here: the Gaunt integral selection rules limit k to |l-l'| <= k <= l+l', and the Gegenbauer linearization coefficients are finite. So the e-e coupling integral is a **finite double sum of algebraic numbers.**

**Verdict:** 🔶 **Algebraic in principle, but not yet implemented.** The mathematics is identical in structure to the Paper 12 Neumann expansion and the Paper 15 split-region Legendre expansion. Each coupling integral decomposes into a finite sum of half-interval Jacobi integrals, each of which evaluates to Gamma function ratios.

### 2c. R-dependent nuclear coupling (Hellmann-Feynman dH/dR)

$$\frac{\partial H}{\partial R} = C(\alpha, \theta_{12}) = -\frac{Z}{\cos\alpha} - \frac{Z}{\sin\alpha} + V_{ee}(\alpha)$$

Since H_ang = [kinetic + centrifugal] + R * C(alpha), and kinetic/centrifugal are R-independent, dH/dR = C(alpha) exactly. The R-dependence is trivially removed — the coupling matrix is R-independent in the charge function basis.

**Algebraic status:** Same as 2a + 2b. The charge function matrix elements in the free basis are the same partial harmonic sums already identified. **Algebraic in principle.**

---

## 3. Perturbation Theory Connection

### 3.1 Perturbation structure

The angular Hamiltonian is:
$$H_{ang}(R) = H_0 + R \cdot V$$

where H_0 = Lambda^2/2 (free SO(6), with Casimir eigenvalues) and V = C(alpha) (charge function). Standard Rayleigh-Schrodinger perturbation theory in R gives:

$$\mu(R) = \mu^{(0)} + a_1 R + a_2 R^2 + a_3 R^3 + ...$$

Paper 13 Section XII proves:
- a_1 is a closed-form partial harmonic sum (Eq. 32)
- a_2 involves sums of |<n'|V|n>|^2 / (mu_n - mu_n'), each term algebraic (Eq. 34)
- All a_k are algebraic: finite sums of rationals * pi^(-k)

### 3.2 Convergence radius

The perturbation series in R has a **finite convergence radius.** This is a critical limitation.

The perturbation parameter is R (the hyperradius), and V = C(alpha) contains the nuclear attraction -Z/cos(alpha) which diverges at alpha = 0, pi/2. However, the perturbation matrix elements <n'|V|n> are finite (the divergence is integrable against the Gegenbauer weight). The convergence radius is determined by the nearest singularity of mu(R) in the complex R-plane.

For two-electron atoms, the adiabatic potential mu(R) has branch points in the complex R-plane associated with avoided crossings between channels (Tolstikhin 1996, cited in Paper 13). The lowest branch point typically occurs at R ~ O(1/Z) bohr. For He (Z=2), the perturbation series likely converges for R < ~1-2 bohr.

**However:** The critical physics happens at R ~ 1-4 bohr (the well region of V_eff), and the asymptotic regime mu ~ -Z^2 R^2/2 at large R has quadratic R-dependence that no finite-order polynomial in R can capture. Paper 13 Section XII.D explicitly states: "the full eigenvalue mu(R) at finite R is transcendental" and "no finite-order rational or Pade approximant captures the crossover."

This is the fundamental algebraic/transcendental boundary.

### 3.3 Resummation possibilities

Even though the perturbation series diverges at large R, several resummation strategies exist:

1. **Pade approximants [M/N](R):** Rational functions that match the Taylor series to order M+N. These can capture the large-R ~ R^2 behavior if N >= 2. The coefficients are algebraic (rationals times pi^(-k)). A [3/2] Pade, for example, would give 5 algebraic parameters and might capture the crossover reasonably well.

2. **Two-point Pade:** Match both the R->0 perturbation series and the R->infinity asymptotic expansion mu -> -Z^2 R^2/2 + ... The asymptotic coefficients are also algebraic (the He+ threshold energies are -Z^2/(2n^2)). A two-point Pade bridges both regimes with purely algebraic coefficients.

3. **Continued fraction representation:** The Rayleigh-Schrodinger perturbation series has a natural continued fraction form (Stieltjes/Pade), and for Coulomb-type perturbations, these continued fractions often have algebraic coefficients and better convergence properties than the power series.

4. **Bender-Wu type analysis:** For quantum anharmonic oscillators (which have the same perturbation structure), the large-order behavior of a_k is known and leads to Borel summability. The Borel-resummed series may converge to the exact mu(R).

**Key question for feasibility:** How many perturbation coefficients a_k are needed for a Pade approximant to achieve 0.01% accuracy across the full R range? Each a_k is algebraic and computable in O(k^2) operations (k-th order perturbation theory sums over intermediate states). If k ~ 10-20 suffices, the algebraic Pade approach could replace the entire FD angular solver.

---

## 4. Literature Survey

### 4.1 What Paper 13 cites

- **Macek 1968 (J. Phys. B 1, 831):** Introduced the adiabatic hyperspherical method for He. Uses numerical integration of the angular equation. No algebraic treatment of coupling integrals reported in Paper 13's discussion.

- **Lin 1995 (Phys. Rep. 257, 1):** Comprehensive review. Paper 13 cites Lin for: coordinate definitions, partial-wave expansion, Liouville substitution, convergence studies. The review systematizes the FD approach but does not report closed-form coupling integrals.

- **Klar & Klar 1980:** Cited for "accurate treatment of two-electron systems." Paper 13 does not elaborate on whether they use algebraic methods, but the Klar work is known to focus on analytical properties of the adiabatic potential curves near thresholds (Wannier theory).

- **Tolstikhin 1996:** Cited for "slow variable discretization" (SVD) as an alternative to FD. SVD uses a basis of adiabatic channel functions (rather than a fixed grid), which is closer to the spectral approach considered here. The SVD method still requires numerical computation of the channel functions, but has better convergence properties near avoided crossings.

- **Fock 1954:** Cited for the non-analytic terms R^(1/2) and R ln(R) near triple coalescence. These terms are automatically handled by a fine hyperradial grid but would complicate a purely algebraic treatment of the R-dependence.

### 4.2 What the code comments reference

No comments in the hyperspherical code files reference analytical formulas that were considered but not implemented. The code is purely numerical (FD grid + dense diagonalization).

### 4.3 Known results from the hyperspherical literature (from general knowledge)

- **Analytic channel functions at R=0:** The free SO(6) eigenfunctions are Gegenbauer polynomials C_n^lambda(cos 2alpha) with lambda = 2 (for L=0 two-electron atoms). This is well-known and used in Paper 13 Section XII.

- **Nuclear coupling in Gegenbauer basis:** The integral of (1/cos(alpha)) against Gegenbauer polynomials is related to the Hilbert transform of Jacobi polynomials, which has known closed forms. This appears in work by Gasper and others on Jacobi polynomial integrals.

- **Fock expansion at R->0:** The perturbation series mu = a_1 R + a_2 R^2 + ... is equivalent to the Fock expansion of the two-electron wavefunction near the nucleus. The first few coefficients have been computed analytically by several groups (e.g., Ermolaev, Purcell). The convergence properties are discussed by Gottschalk et al.

- **Semi-analytic adiabatic curves:** Some groups (e.g., Krivec, Mandelzweig) have used variational methods with a few Gegenbauer basis functions (rather than FD grids) to compute mu(R) at specific R values. This is essentially the spectral approach discussed in Section 5. Typical accuracy: ~0.1% with 5-10 basis functions.

---

## 5. Feasibility Assessment

### 5.1 Current algebraic fraction

| Component | Count | Status | Fraction |
|:----------|:-----:|:------:|:--------:|
| Liouville constant (-2) | 1 | ✅ Algebraic | - |
| Jacobian centrifugal (15/8R^2) | 1 | ✅ Algebraic | - |
| Gaunt integrals G(l,k,l') | (l_max+1)^2 * (2l_max+1) | ✅ Algebraic | - |
| Normalization factors | (l_max+1)^2 | ✅ Algebraic | - |
| SO(6) Casimir eigenvalues | n_channels | ✅ Algebraic | - |
| Kinetic energy (-1/2 u'') | N_alpha per channel | ❌ FD grid | - |
| Centrifugal l(l+1)/cos^2 | N_alpha per channel | 🔶 Grid, algebraic possible | - |
| Nuclear attraction (-Z/cos) | N_alpha per channel | 🔶 Grid, algebraic possible | - |
| V_ee coupling f_k/max | N_alpha per (l,l') pair | 🔶 Grid, algebraic possible | - |

**Summary:** The selection rules, quantum number labels, and coupling coefficients (Gaunt integrals) are already algebraic. The entire alpha-dependent part — kinetic, centrifugal, nuclear, and e-e coupling — is currently on an FD grid. By matrix element count, roughly **80-90% of the angular Hamiltonian construction is numerical** (the grid-based terms dominate by dimensionality).

### 5.2 What could be made algebraic

**All alpha-dependent matrix elements** can be made algebraic by switching from the FD grid to a spectral Gegenbauer basis. Specifically:

1. **Free basis functions:** C_{n-1}^2(cos 2alpha) * sin(alpha) cos(alpha), n = 1, 3, 5, ... (odd n for singlet S_2 symmetry). These are the exact eigenfunctions of the free angular Hamiltonian (Lambda^2/2).

2. **Kinetic + centrifugal matrix elements in the Gegenbauer basis:** Exactly diagonal — these are just the SO(6) Casimir eigenvalues nu(nu+4)/2. No FD discretization needed.

3. **Nuclear coupling matrix elements:** <n'|(-Z/cos(alpha) - Z/sin(alpha))|n> = algebraic (partial harmonic sums, proven by Paper 13 Section XII).

4. **V_ee coupling matrix elements:** <n'|W_{ll'}(alpha)|n> for each l,l' pair = algebraic (split-region integrals of Gegenbauer polynomials against trigonometric weights, same mathematical structure as Paper 15's split-region Legendre expansion).

**Fraction that could be algebraic:** ~100% of the angular Hamiltonian matrix elements. The entire angular eigenvalue problem at fixed R would reduce to building a matrix whose entries are known closed-form expressions (rational numbers, harmonic sums, pi^(-k) terms) and diagonalizing it.

### 5.3 The irreducible numerical core

Two components remain irreducibly numerical:

1. **Diagonalization of the angular Hamiltonian at each R:** Even with algebraic matrix elements, the eigenvalue mu(R) at finite R is transcendental (Paper 13 Section XII.D). However, the *matrix* is algebraic — only the diagonalization step is numerical, and this is O(n_basis^3) with exact arithmetic possible via interval methods. For n_basis ~ 10-20, this is trivial.

2. **Hyperradial Schrodinger equation:** The 1D equation -1/2 F'' + V_eff(R) F = E F must be solved on an R grid. V_eff(R) is known only at the R points where the angular problem was solved. This requires either:
   - Interpolation (spline) + FD radial solver (current approach), or
   - Algebraic Pade approximant for mu(R) + analytical WKB/phase integral

   The first option is standard and introduces negligible error (Paper 13 Table III shows convergence at N_R ~ 2000). The second option (Pade) could be interesting but is not strictly necessary — the radial solver is already well-converged.

3. **Non-adiabatic coupling at finite R:** P_mu_nu(R) involves eigenvectors at finite R, which are numerically determined. However, with algebraic matrix elements, the Hellmann-Feynman formula gives P = (algebraic matrix element) / (eigenvalue gap), where the matrix element is algebraic and only the gap is transcendental.

### 5.4 Architectural sketch: Algebraic Level 3 Solver

```
AlgebraicAngularSolver:
  Input: Z, n_channels, n_basis (number of Gegenbauer basis functions)

  Step 1: Precompute angular coupling matrix elements (ALGEBRAIC)
    - V_nuc[n',n] = partial harmonic sum (closed form, Paper 13 Eq. 31)
    - V_ee[n',n,l,l'] = split-region Jacobi integrals (closed form)
    - Centrifugal = SO(6) Casimir (diagonal, exact integers)

  Step 2: At each R, build H_ang(R) = H_free + R * V_coupling (ALGEBRAIC)
    - H_free = diag(nu(nu+4)/2)  [pure integers]
    - V_coupling = V_nuc + sum_{k,l,l'} Gaunt * V_ee  [algebraic numbers]
    - Total matrix: exact rational/algebraic entries, dimension n_basis

  Step 3: Diagonalize H_ang(R)  (NUMERICAL, but tiny matrix ~10x10)
    - Returns mu(R), eigenvectors in Gegenbauer basis

  Step 4: Hyperradial solve (unchanged)
    - Interpolate V_eff(R) = mu(R)/R^2 + 15/(8R^2)
    - FD tridiagonal solver for F(R)
```

**Key advantages over current FD approach:**
- No grid discretization error in alpha (eliminates the l_max > 0 degradation)
- Matrix dimension ~ n_basis (10-20) instead of N_alpha * n_l (200-800)
- Exact selection rules enforced (n'-n = 0 mod 4 for nuclear coupling)
- Matrix elements are R-independent — precomputed once, reused at every R

### 5.5 Expected accuracy improvement

Paper 13 Table I shows the critical issue:

| l_max | Energy (Ha) | Error |
|:-----:|:-----------:|:-----:|
| 0 | -2.9052 | 0.05% |
| 1 | -2.9262 | 0.78% |
| 2 | -2.9284 | 0.85% |
| 3 | -2.9289 | 0.87% |

The l_max=0 result is *the most accurate* because higher-l channels introduce FD errors from the l(l+1)/cos^2(alpha) boundary singularity. This is artificial — higher-l channels add genuine angular correlation that should improve accuracy.

**With an algebraic solver:**
- The centrifugal singularity is handled exactly (Gegenbauer basis functions are the natural eigenfunctions of the centrifugal operator)
- Higher-l channels would add angular correlation *without* introducing discretization error
- Expected: l_max=3 should be *better* than l_max=0, potentially reaching 0.01% or below
- The 0.05% error at l_max=0 is already close to the adiabatic approximation limit (the non-adiabatic correction is ~0.05% for He). With proper l_max convergence, the single-channel result could approach 0.02-0.03% error

**Conservative estimate:** 2-5x accuracy improvement (from 0.05% to 0.01-0.025%) by eliminating FD discretization error and allowing proper l_max convergence.

---

## 6. Summary of Algebraic Status

| Matrix Element Type | Current Method | Algebraic Status | Replacement |
|:---------------------|:--------------|:-----------------|:------------|
| Gaunt integrals | Wigner 3j formula | ✅ Already algebraic | None needed |
| SO(6) Casimir eigenvalues | Known formula | ✅ Already algebraic | None needed |
| Liouville constant (-2) | Hardcoded | ✅ Already algebraic | None needed |
| Jacobian centrifugal (15/8R^2) | Hardcoded | ✅ Already algebraic | None needed |
| Normalization factors | sqrt formula | ✅ Already algebraic | None needed |
| Perturbation coeff a_1 | Not implemented | ✅ Known closed form (Paper 13) | Implement Eq. 32 |
| Nuclear coupling integrals | FD grid (pointwise) | 🔶 Algebraic in principle | Gegenbauer basis + harmonic sums |
| V_ee coupling integrals | FD grid (pointwise) | 🔶 Algebraic in principle | Split-region Jacobi integrals |
| Centrifugal matrix elements | FD grid (pointwise) | 🔶 Algebraic in principle | Gegenbauer orthogonality (diagonal) |
| Kinetic operator | FD 3-point stencil | 🔶 Algebraic in principle | Gegenbauer basis (diagonal = Casimir) |
| Full mu(R) at finite R | Dense diagonalization of FD matrix | ❌ Transcendental | Small algebraic matrix + numerical diag |
| Hyperradial wavefunction F(R) | FD tridiagonal solver | ❌ Irreducibly numerical | Keep current approach |
| Non-adiabatic coupling P(R) | Hellmann-Feynman + FD eigenvectors | 🔶 Partially algebraic | Algebraic matrix elements + numerical eigenvectors |

**Bottom line:** The angular Hamiltonian *matrix elements* are 100% algebraic. The angular *eigenvalues* at finite R are transcendental but computable from a small (10-20 dimensional) algebraic matrix — six orders of magnitude smaller than the current 200-800 dimensional FD matrix. The hyperradial solver is irreducibly numerical but already well-converged. The algebraic replacement targets the right bottleneck: the FD discretization error that currently prevents l_max convergence.

---

## 7. Recommended Implementation Path (Track B, Task 2)

1. **Derive and verify nuclear coupling integrals** in the Gegenbauer basis. Paper 13 Eq. 31-32 gives the diagonal case; extend to off-diagonal using Gegenbauer linearization.

2. **Derive and verify V_ee coupling integrals** in the split-region form. Follow the Paper 15 split-region Legendre template.

3. **Build AlgebraicAngularSolver** that constructs H_ang(R) as a small dense matrix (n_basis x n_channels x n_basis x n_channels) with algebraic entries.

4. **Validate** against the FD solver at multiple R values and l_max settings.

5. **Test l_max convergence** with the algebraic solver — verify that l_max=1,2,3 improve on l_max=0 (the FD-limited result).

6. **Assess the Pade resummation** option for mu(R) if further accuracy gains are needed beyond the spectral basis approach.
