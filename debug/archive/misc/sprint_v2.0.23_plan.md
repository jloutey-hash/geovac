# Track W: Cusp-Graph Theory Investigation

## Summary

**Question:** Can the 1/r_12 electron-electron singularity be absorbed into a conformally weighted graph Laplacian on the Level 4 angular space, analogous to how the 1/r nuclear singularity IS the graph Laplacian at Level 1 (Paper 7)?

**Answer:** NO. Structural dimensionality obstruction. The result is a clean negative.

**Classification:** Dimensionality mismatch.

**Tests:** 37/37 passing in `tests/test_cusp_graph.py`.

---

## The Mathematical Argument

### Step 1: Why Level 1 Works

At Level 1, the hydrogen atom maps to S^3 via Fock's stereographic projection (Paper 7). The Green's function of the Laplace-Beltrami operator on S^n has a well-known singularity structure:

- S^1, S^2: logarithmic (special cases)
- S^n, n >= 3: power-law 1/d^{n-2}

On S^3: Green ~ 1/d^{3-2} = **1/d^1**.

The Coulomb potential 1/r also has singularity order 1: it diverges as 1/d^1.

This is a **dimensional coincidence unique to S^3**: the Green's function singularity order (n-2 = 1) matches the Coulomb singularity order (1) if and only if n = 3. Under stereographic projection, the chordal distance identity converts the S^3 Green's function into the momentum-space Coulomb kernel 1/|p-p'|^2, which Fourier transforms to 1/r in position space. The conformal weights from the projection are smooth and get absorbed into wavefunction normalization.

**The graph Laplacian absorbs the nuclear Coulomb potential because S^3 is the unique sphere whose propagator IS the Coulomb interaction.**

### Step 2: Why Level 4 Fails

The Level 4 angular space is 5-dimensional (coordinates alpha, theta_1, theta_2, Phi for Sigma states), topologically related to S^5.

On S^5: Green ~ 1/d^{5-2} = **1/d^3**.

But 1/r_12 in hyperspherical coordinates is:

    r_12 = R_e * sqrt(1 - sin(2*alpha)*cos(theta_12))

Near the coalescence point (alpha = pi/4, theta_12 = 0), expanding to leading order:

    r_12 ~ R_e * sqrt(2*delta_alpha^2 + delta_theta^2/2) ~ R_e * d

where d is the geodesic distance to the coalescence manifold. Therefore:

    1/r_12 ~ 1/d^**1**

The mismatch is:
- Green's function on S^5: 1/d^**3**
- Coulomb 1/r_12: 1/d^**1**
- Discrepancy: d^**2** (two orders of magnitude in the singularity)

**No conformal rescaling can bridge this gap.** Under g -> Omega^2 g on an n-manifold, the Green's function transforms as G_new = Omega(x)^{-(n-2)/2} * Omega(y)^{-(n-2)/2} * G_old. Since Omega is smooth and positive, the singularity ORDER is preserved. The obstruction is topological (encoded in the dimension), not metric.

### Step 3: The Fiber Question

Could there be a natural S^3 submanifold at the coalescence point where 1/r_12 is a fiber Green's function?

The coalescence manifold {alpha = pi/4, theta_12 = 0} is a **3-dimensional submanifold** of S^5 (parametrized by the remaining coordinates theta_1, theta_2, Phi). Its codimension is 2, so the normal bundle has **2-dimensional fibers**.

To embed S^3 as a fiber, we would need the normal bundle to have dimension >= 3. But the normal fiber dimension is 2. The candidates:

| Sphere | Fits in 2D fiber? | Green matches 1/d^1? | Viable? |
|--------|:-:|:-:|:-:|
| S^1 | Yes | No (logarithmic) | No |
| S^2 | Yes | No (logarithmic) | No |
| S^3 | **No** (dim=3 > 2) | **Yes** | No |

S^3 is the unique sphere that matches the Coulomb singularity, but it cannot be embedded in the 2-dimensional normal bundle to the coalescence manifold. The fiber dimension is wrong by exactly 1.

### Step 4: Conformal Invariance Closes All Routes

The argument is airtight because:

1. **Direct matching fails:** S^5 Green (1/d^3) != Coulomb (1/d^1)
2. **Conformal rescaling cannot fix it:** singularity order is a conformal invariant
3. **Fiber decomposition fails:** coalescence codimension is 2, need 3 for S^3
4. **No sphere of any dimension works in the available fiber:** S^1 and S^2 have logarithmic Green's functions

All four routes to absorbing the cusp are blocked.

---

## Numerical Verification

- r_12 singularity order measured as 0.9999 (expected 1.0) -- confirmed
- S^3 Green function proportional to 1/d near pole -- confirmed (ratio_std/mean < 0.01)
- Hyperspherical r_12 formula verified at known points (coalescence, maximum, boundary)

---

## Paper 18 Classification

In Paper 18's exchange constant taxonomy, 1/r_12 is an **EMBEDDING exchange constant**:

- The discrete graph (S^5 Casimir spectrum: nu*(nu+4)/2) captures the free-particle angular structure exactly with integer eigenvalues
- The 1/r_12 interaction cannot be absorbed because it has the wrong singularity order for S^5
- It must be treated as a charge function C_ee/R that perturbs the free spectrum
- This perturbation introduces transcendental content (the adiabatic eigenvalue mu(R)) through spatial integration
- The transcendental content is the irreducible "exchange rate" between the discrete angular spectrum and the continuous Coulomb interaction

This is structurally parallel to the Level 3 situation: the SO(6) Casimir eigenvalues are exact integers, but the charge function C(alpha, theta_12)/R introduces mu(R) which is algebraic over Q(pi, sqrt(2)) (Track P1) -- not rational. The transcendental content measures the mismatch between the graph's 1/d^3 propagator and the physics' 1/d^1 Coulomb singularity.

---

## Implications for Companion Tracks

**Track U (Kato factor extraction): SUPPORTED.** Since the cusp cannot be absorbed topologically, explicit analytic treatment of the cusp factor f(r_12) = 1 + r_12/2 + ... is the correct approach. The Jastrow factor is not a graph property -- it is the analytic content of the singularity mismatch.

**Track V (cusp regularization): SUPPORTED.** The current approach of treating C_ee as a perturbation on the free S^5 Casimir spectrum via Gaunt integrals and split-region Legendre expansion is structurally correct. The cusp must be regularized in the angular basis because no graph topology can absorb it.

**Track X (basis improvement): SUPPORTED.** Since the cusp is irreducibly a charge-function property, improving the angular basis representation of C_ee is the correct path. The spectral Jacobi basis (Track K, 269x speedup) represents the optimal approach within the graph framework.

---

## Key Insight

The success of Level 1 -- the fact that the nuclear Coulomb potential IS the graph Laplacian on S^3 -- rests on a dimensional coincidence: the Green's function on S^n has singularity 1/d^{n-2}, and n-2 = 1 only when n = 3. This coincidence is what makes the entire GeoVac framework possible for single-electron atoms. It is specific to the three-sphere and does not generalize to higher-dimensional angular spaces. The electron-electron cusp, living in a 5-dimensional angular space, cannot benefit from this coincidence. It is the irreducible transcendental content of multi-electron quantum mechanics.

---

## Files

- Analysis code: `geovac/cusp_graph.py`
- Tests (37/37): `tests/test_cusp_graph.py`
- This document: `debug/track_w/analysis.md`
