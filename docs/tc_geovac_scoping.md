# Transcorrelated GeoVac: Feasibility Scoping

**Track BX-1 | April 2, 2026**

## 1. Background

The electron-electron cusp is the documented accuracy bottleneck for GeoVac's Level 4 (H2) and Level 3 (He) solvers. Paper 12 diagnosed a 7.6% D_e gap arising from the Kato cusp condition: the exact wavefunction has a derivative discontinuity at r_12 = 0 that no smooth partial-wave expansion can represent at finite l_max. The partial-wave energy error converges as 1/(l+1/2)^4 (Schwartz 1962), making brute-force convergence impractical. The current best H2 result is 96.0% D_e at l_max=6 with CBS extrapolation ~97%.

Three direct attacks on the cusp have failed (CLAUDE.md Section 3): (i) an alpha-only Jastrow factor f(alpha) that ignores the theta_12 dependence, (ii) absorption of 1/r_12 into a conformally weighted S^5 graph Laplacian (Green function singularity mismatch: 1/d^3 on S^5 vs 1/d^1 for Coulomb), and (iii) a theta_12-adapted angular basis. The cusp is 2D in (alpha, theta_12) and structurally irreducible. Paper 18 classifies 1/r_12 as an "embedding exchange constant" -- the transcendental price of resolving the interacting two-electron problem within the hyperspherical angular space.

The transcorrelated (TC) method takes a different approach: rather than improving the basis to represent the cusp, it transforms the Hamiltonian to remove the cusp entirely. Given a Jastrow factor J(r_12), the similarity transformation H_TC = e^{-J} H e^{J} produces a non-Hermitian Hamiltonian whose eigenfunctions are smooth (cusp-free), converging exponentially rather than algebraically in the partial-wave expansion. The key question is whether this transformation is compatible with GeoVac's angular structure.

## 2. Angular Symmetry Analysis

**Does J = f(r_12) preserve Gaunt selection rules?**

Yes, with a nuance. The Jastrow factor f(r_12) depends only on the interelectron distance |r_1 - r_2|, which is rotationally invariant: it commutes with the total angular momentum operators L^2 and L_z. The similarity transformation e^{-J} H e^{J} therefore preserves total angular momentum quantum numbers. The TC Hamiltonian H_TC has the same block-diagonal structure in L as the original.

The Gaunt selection rules survive because they derive from rotational symmetry (Wigner 3j symbols enforce the triangle inequality on angular momentum coupling). Since J is a scalar under rotations, e^{-J} preserves which angular channels couple -- no new channels beyond those allowed by the triangle inequality appear.

However, **sparsity within the allowed channels is not preserved**. The Baker-Campbell-Hausdorff (BCH) expansion generates new operator structures:

- [T, J] produces one-body and two-body gradient terms (nabla_i J) that, when expanded in the partial-wave basis, populate more of the allowed radial integrals within each angular block.
- [[T, J], J] produces (nabla J)^2 terms -- a multiplicative potential that couples all radial channels within each angular block.

In GeoVac's factorized ERI structure (angular Gaunt x radial Slater R^k), the TC transformation leaves the angular Gaunt factors unchanged but modifies and densifies the radial integrals. The number of nonzero R^k-type integrals within each angular channel will increase.

**Bottom line:** Angular quantum numbers and selection rules are preserved. The angular sparsity (which channels couple) survives. The radial sparsity (how many integrals per channel are nonzero) degrades.

## 3. TC-GeoVac Implementation Architecture

**BCH expansion.** For a two-electron system with J = f(r_12):

    H_TC = H + [H, J] + (1/2)[[H, J], J] + ...

For a two-electron system with a symmetric Jastrow, the BCH series terminates at second order in the kinetic commutators (since [V, J] = 0 when V and J are both functions of coordinates only). The resulting TC Hamiltonian has three new terms beyond the original H:

1. **One-body gradient term** from [T_1 + T_2, J]: a first-derivative operator -(nabla_i f) . nabla_i. In hyperspherical coordinates, nabla_i J involves partial derivatives of f(r_12) with respect to (alpha, theta_12), which are smooth functions away from the cusp point.

2. **Scalar potential** from (nabla J)^2: a multiplicative (local) potential (1/2)|nabla f|^2. For the exact cusp Jastrow J = -(1/2)r_12, this equals 1/4 -- a constant that shifts all energies uniformly and is trivially absorbed.

3. **Modified V_ee**: The cusp-causing 1/r_12 singularity is cancelled by the -nabla^2 J term. For J = -(1/2)r_12, the Laplacian nabla^2(r_12) = 2/r_12 (in 3D), so (-1/2)(-1/2)(2/r_12) = 1/(2r_12) partially cancels 1/r_12. With the full BCH at second order, the net effect is that the 1/r_12 divergence is removed and replaced by smooth, finite potentials.

**New integrals required.** In GeoVac's Level 4 hyperspherical basis {Phi_nu(alpha) Y_{l,m}(theta, phi)}, the TC terms require:

- <nu, l | nabla_alpha f | nu', l'>: gradient matrix elements in the hyperangle. These involve derivatives of r_12(alpha, theta_12) = R_e sqrt(1 - sin(2alpha) cos(theta_12)), which are smooth trigonometric functions -- computable analytically via the same Gaunt machinery.
- <nu, l | (nabla f)^2 | nu', l'>: products of trigonometric functions -- straightforward angular integrals.
- <nu, l | nabla^2 f | nu', l'>: the Laplacian of r_12 in hyperspherical coordinates. This has a 1/r_12 singularity (from the 2/r_12 term in nabla^2 r_12), but it exactly cancels the original V_ee singularity -- that is the entire point of the TC transformation.

The resulting TC Hamiltonian matrix is non-Hermitian (the similarity transformation is not unitary). This means eigenvalues are real (provably, for the exact Jastrow) but left and right eigenvectors differ. Standard symmetric eigensolvers (Lanczos, Davidson) cannot be used; one needs either:
- A biorthogonal solver (left + right eigenvectors).
- The Hermitized variant H_TC^dag H_TC (doubles the matrix size but restores Hermiticity).
- The xTC (explicitly transcorrelated) approach, which uses a Hermitian approximation to H_TC by symmetrizing the one-body terms.

**Three-body terms for N > 2 electrons.** For systems with more than two electrons (Level 5 composed molecules), the TC transformation generates genuine three-body operators from [V_ee, J] when J depends on a single pair r_ij. These three-body terms couple three electrons simultaneously and cannot be decomposed into one- and two-body operators. Options:
- Use a pair-specific Jastrow J = sum_{i<j} u(r_ij): three-body terms arise from cross-pair kinetic commutators [T_i, u(r_jk)]. These are zero for the isotropic case (i not in {j,k}), so the only three-body terms come from [nabla_i^2, u(r_ij)] acting on pairs sharing electron i.
- Approximate: the NOCI-TC approach truncates three-body terms. Error is typically < 1 mHa for light atoms.
- For GeoVac's composed geometry (Level 5), where electron pairs are in separate blocks, three-body terms between blocks are zero by construction (no shared electrons). Within-block terms are present but tractable (each block has at most 2 electrons for first-row atoms).

## 4. Minimum Viable Experiment: TC He

**System:** He (2 electrons, Level 3 hyperspherical). The simplest two-electron system where BCH terminates exactly.

**Jastrow:** J = -(1/2) r_12. This is the exact Kato cusp condition. In Level 3 hyperspherical coordinates: r_12 = R sqrt(1 - sin(2alpha) cos(theta_12)), where R is the hyperradius and (alpha, theta_12) are hyperangular coordinates.

**Implementation steps:**

1. Compute the three TC correction matrices in the existing Level 3 angular basis {Gegenbauer polynomials C_nu^2(cos 2alpha) x P_l(cos theta_12)}:
   - nabla J matrix elements: derivatives of r_12 with respect to (alpha, theta_12), expanded in the angular basis. These are trigonometric/algebraic functions -- computable from Gaunt-like integrals.
   - (nabla J)^2 matrix: a scalar potential; standard angular integration.
   - nabla^2 J matrix: contains 1/r_12 that cancels V_ee; the net V_ee + nabla^2 J is smooth and integrable.

2. Assemble H_TC = H_original + H_correction as a non-Hermitian matrix.

3. Solve for eigenvalues using a general (non-symmetric) eigensolver (numpy.linalg.eig or scipy sparse equivalent).

4. Compare: standard He at l_max=2 gives 0.24% error (0.10% with Schwartz post-correction). TC-He at l_max=2 should converge exponentially rather than algebraically, potentially reaching < 0.05% without any post-correction.

**Comparison targets:**

| Method | l_max | Error (%) | Notes |
|:-------|:-----:|:---------:|:------|
| Standard Level 3 | 0 | 0.16 | Single-channel, algebraic basis |
| Standard Level 3 | 2 | 0.24 | Uncorrected |
| Standard + Schwartz | 2 | 0.10 | Post-processing correction |
| TC Level 3 (target) | 2 | < 0.05 | Exponential convergence expected |
| Exact | -- | 0.00 | -2.9037 Ha |

If TC at l_max=2 beats standard at l_max=5 (without Schwartz), the basis compression is demonstrated.

**Key technical risk:** The r_12 function in hyperspherical coordinates contains a square root: r_12 = R sqrt(1 - sin(2alpha) cos(theta_12)). Its gradient has a 1/sqrt(...) singularity at the coalescence point (alpha=pi/4, theta_12=0), and its Laplacian has a 1/r_12 singularity. After cancellation with V_ee, the net potential is smooth -- but verifying this cancellation numerically at the matrix element level requires care (analytic cancellation before discretization, not numerical subtraction of large terms).

## 5. Assessment

| Criterion | Status |
|:----------|:-------|
| Angular quantum numbers preserved? | Yes -- J is rotationally invariant |
| Gaunt selection rules preserved? | Yes -- triangle inequality unchanged |
| Angular sparsity preserved? | Yes -- same allowed channels |
| Radial sparsity preserved? | No -- more nonzero integrals per channel |
| BCH terminates (2 electrons)? | Yes -- finite order |
| Three-body terms (N > 2)? | Present but manageable in composed geometry |
| Non-Hermiticity? | Requires modified solver (biorthogonal or xTC) |
| New integral types needed? | Yes -- gradient and Laplacian of r_12 in hypersp. basis |
| Analytic cancellation feasible? | Likely -- r_12 has known analytic structure |

**Feasibility classification: REQUIRES NEW INTEGRALS.**

The TC transformation is theoretically compatible with GeoVac's angular framework. The selection rules and angular quantum numbers are preserved. The main implementation cost is computing new radial/hyperangular matrix elements for the gradient and Laplacian of the Jastrow factor in the existing basis. These integrals are non-trivial but analytically structured (trigonometric functions in hyperspherical coordinates). The non-Hermitian eigenproblem is a secondary complication, solvable with standard numerical tools.

The approach is NOT straightforward (it requires new integral evaluation code and a non-Hermitian solver), but it is also not infeasible. The critical path is verifying that the V_ee + nabla^2 J cancellation produces numerically stable matrix elements.

## 6. Track BX-2 Plan (Conditional)

If the PI approves proceeding, the implementation would be split into three sub-agent tasks:

**Sub-agent 1: TC integral derivation (theory).**
- CONTEXT: Paper 13 (hyperspherical coordinates), Paper 15 Sec V (cusp location in Level 4).
- TASK: Derive analytic expressions for nabla r_12, (nabla r_12)^2, and nabla^2 r_12 in Level 3 hyperspherical coordinates (R, alpha, theta_12). Verify the V_ee cancellation algebraically. Express all TC correction terms as expansions in the Gegenbauer x Legendre angular basis.
- SUCCESS: Closed-form matrix element expressions for all TC correction terms.

**Sub-agent 2: TC He implementation.**
- CONTEXT: `geovac/hyperspherical.py`, `geovac/algebraic_angular.py`, sub-agent 1 output.
- TASK: Implement `tc_correction_matrices(l_max, R_grid)` returning the three TC correction matrices. Assemble H_TC = H + H_corr. Solve with `numpy.linalg.eig`. Compare against standard solver at l_max=0,1,2,3.
- SUCCESS: TC He at l_max=2 achieves < 0.10% error without Schwartz post-correction.
- CONSTRAINTS: Do NOT modify existing solver code. TC solver is a new module `geovac/tc_solver.py`.

**Sub-agent 3: Sparsity impact analysis.**
- CONTEXT: `geovac/composed_qubit.py`, Paper 14, sub-agent 2 output.
- TASK: Count nonzero matrix elements in H_TC vs H at matched l_max. Quantify the radial densification. Estimate Pauli term count for a TC qubit Hamiltonian.
- SUCCESS: Report with sparsity ratios and projected Pauli scaling exponent.
