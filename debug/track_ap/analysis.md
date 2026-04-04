# Track AP: 2D Variational Solver Scoping for 4 Electrons

Date: 2026-04-01
Status: SCOPING COMPLETE

## Executive Summary

**Verdict: FEASIBLE WITH OPTIMIZATION at l_max=2; REQUIRES ALGORITHMIC BREAKTHROUGH at l_max>=3.**

A full 2D (hyperradial + all hyperangles) variational solver for the 4-electron
problem is feasible at l_max=2 using iterative eigensolvers on sparse matrices,
with estimated wall time ~10-30 min/point. At l_max=3 the tensor product grows
to 62,500 and becomes borderline. A PARTIAL 2D approach -- treating R_e + alpha2
(the inter-pair hyperangle) while keeping alpha1, alpha3 adiabatic -- is the
recommended path, reducing dimensions by 5x while capturing the dominant
non-adiabatic coupling.

---

## (a) DIMENSIONS: Hyperangular Structure

### Coordinate system (from n_electron_solver.py)

4 electrons in mol-frame hyperspherical coordinates have:
- 1 hyperradius R_e (radial)
- 3 hyperangles: alpha1 (pair 1-2), alpha2 (inter-pair), alpha3 (pair 3-4)
- 8 direction angles: (theta_i, phi_i) for i=1..4

After angular momentum channel expansion (integrating out the 8 direction
angles into (l_i, m_i) labels), 3 continuous hyperangles remain.

### Fast vs slow classification

The hyperangles have different physical roles:

| Angle   | Physical meaning          | R_e coupling | Classification |
|---------|--------------------------|--------------|----------------|
| alpha2  | Inter-pair (core/valence) | STRONG       | FAST           |
| alpha1  | Pair 1-2 (intra-core)     | MODERATE     | SLOW           |
| alpha3  | Pair 3-4 (intra-valence)  | MODERATE     | SLOW           |

**alpha2 is the fast variable.** It controls the ratio rho_34/rho_12 (valence
vs core pair radii). As R_e changes, the core/valence partition shifts, and
alpha2 tracks this directly. This is the 4-electron analog of the single alpha
at Level 4 (which partitions two electrons between inner and outer).

alpha1 and alpha3 are the intra-pair angles. They control how each electron
pair distributes between its two members. These are approximately separable
from R_e -- they describe angular correlation within each pair, which is
relatively R_e-independent (the pair structure is set by angular momentum,
not by the overall scale).

**Evidence:** The nuclear coupling V_nuc depends on rho = R/(2*R_e) and
enters through min(s_i, rho)/max(s_i, rho). Since s_i depends on all three
alphas but rho_A, rho_B depend on R_e, the strongest R_e-coupling is through
the alpha2 dependence (which controls the overall scale splitting between
the two pairs). The V_ee coupling is R_e-independent in the charge function
formulation (it only depends on s_i ratios), confirming that alpha1, alpha3
are slow variables.

---

## (b) TENSOR PRODUCT SIZE

### Spectral angular dimensions (from n_electron_spectral.py)

| l_max | Raw channels | S4 [2,2] channels | Basis/channel (5^3) | Spectral dim |
|-------|-------------|-------------------|---------------------|--------------|
| 0     | 1           | 1                 | 125                 | 125          |
| 1     | 16          | 1                 | 125                 | 125          |
| 2     | 81          | 6                 | 125                 | 750          |
| 3     | 256         | 20                | 125                 | 2,500        |
| 4     | 625         | 50                | 125                 | 6,250        |

The spectral dimension is n_channels_S4 * n_basis_per_angle^3.

### Full 2D tensor product (spectral angular x Laguerre radial)

| l_max | Angular dim | Radial dim (n_basis) | Tensor product | Notes |
|-------|-------------|---------------------|----------------|-------|
| 0     | 125         | 25                  | 3,125          | Trivial |
| 1     | 125         | 25                  | 3,125          | Trivial |
| 2     | 750         | 25                  | 18,750         | Target |
| 3     | 2,500       | 25                  | 62,500         | Borderline |
| 4     | 6,250       | 25                  | 156,250        | Intractable |

### Does structure reduce the tensor product?

**Partially.** The Hamiltonian H = T_Re (x) I_ang + H_ang(R_e)/R_e^2 has
Kronecker structure, but H_ang(R_e) depends on R_e through rho = R/(2*R_e),
so the angular blocks are R_e-dependent. This means the full tensor product
must be built (no block-diagonal reduction). However:

1. **Sparsity:** The kinetic energy T_Re is tridiagonal (Laguerre), and
   H_ang is sparse in the spectral basis (SO(12) Casimir is diagonal,
   V_ee coupling is rho-independent, V_nuc has limited channel coupling
   via Gaunt selection rules). The overall matrix is SPARSE.

2. **Reduced radial basis:** n_basis=25 is conservative. At Level 4,
   n_basis=20 is optimal. For the 4-electron case, n_basis=15-20 may
   suffice, reducing the tensor product by 20-40%.

3. **Reduced angular basis:** n_basis_per_angle=5 may be excessive for
   the slow variables alpha1, alpha3. Using n_basis_per_angle=3 for
   these and 5 for alpha2 would give 3*5*3 = 45 basis functions per
   channel instead of 125, reducing spectral dim by 2.8x.

---

## (c) EIGENSOLVE FEASIBILITY

### Dense eigensolve

| Dimension | Flops (N^3)     | Memory (N^2 doubles) | Time estimate    |
|-----------|-----------------|----------------------|------------------|
| 18,750    | 6.6 x 10^12     | 2.8 GB               | ~30 min          |
| 62,500    | 2.4 x 10^14     | 31 GB                | ~15 hours        |
| 156,250   | 3.8 x 10^15     | 195 GB               | ~10 days         |

Dense eigensolve is INFEASIBLE at l_max>=3 due to memory. At l_max=2
it is borderline (2.8 GB fits in RAM but 30 min per PES point is slow
for a 12-point scan).

### Iterative eigensolve (Lanczos/ARPACK via scipy.sparse.linalg.eigsh)

The Hamiltonian is sparse with structure:
- Radial kinetic: tridiagonal in R_e, identity in angular = O(N_ang) nnz per row
- Angular blocks: at each R_e quadrature point, H_ang is n_ch-blocked with
  coupling only between Gaunt-connected channels. Sparsity ~10-20% for small l_max.

Estimated nnz per row: ~50-200 (channel coupling + radial coupling)
Total nnz: ~3.5M-37.5M at l_max=2 (18,750 x 200)

| Dimension | Matvec cost | Iterations (k~100) | Total time      | Feasible? |
|-----------|-------------|---------------------|-----------------|-----------|
| 18,750    | ~7M flops   | 100-200             | ~10-30 sec      | YES       |
| 62,500    | ~75M flops  | 200-400             | ~5-30 min       | YES       |
| 156,250   | ~500M flops | 300-600             | ~1-5 hours      | BORDERLINE|

**Critical issue:** The angular Hamiltonian H_ang(R_e) must be evaluated at
each Gauss-Laguerre quadrature point for the spectral radial basis. With
n_quad ~ 60-80 quadrature points, building H_ang costs:
- l_max=2: 750 x 750 matrix x 80 points = 45M elements to compute
- Each matrix build: ~0.4s (from Track AK timing at spectral dim=750)
- Total: 80 x 0.4s = 32s per PES point (just for matrix assembly)

The matrix assembly dominates, not the eigensolve. This is the same pattern
as Level 4 (Track K finding: angular sweep dominates).

### Sparse structure assessment

The 2D Hamiltonian IS sparse:
- T_Re (x) I_ang: banded (bandwidth = N_ang, 3 bands)
- sum_q w_q L_n(x_q) L_m(x_q) H_ang(R_q): this is a sum of rank-N_ang
  outer products weighted by quadrature, creating a dense-ish coupling.
  However, H_ang itself is sparse (~20% fill at l_max=2 due to Gaunt
  selection rules).

Net sparsity: ~5-15% fill. Enough for iterative methods to win over dense.

---

## (d) COMPARISON TO LEVEL 4

### Level 4 2D solver reference (from level4_multichannel.py)

At Level 4, the 2D solver uses:
- n_ch channels x n_alpha alpha-grid = N_ang (e.g., 5 channels x 200 = 1,000; or spectral: 5 x 10 = 50)
- n_Re radial points (FD: ~300; Laguerre: ~25)
- Tensor product: 50 x 25 = 1,250 (spectral)

### Scaling comparison

| Property              | Level 4 (2e)  | Level 4N (4e, l_max=2) | Ratio |
|-----------------------|---------------|------------------------|-------|
| Hyperangles           | 1 (alpha)     | 3 (alpha1,2,3)         | 3     |
| S4 channels           | 5 (l_max=2)   | 6                      | 1.2   |
| Basis per channel     | 10            | 125 (5^3)              | 12.5  |
| Angular dim           | 50            | 750                    | 15    |
| Radial dim            | 25            | 25                     | 1     |
| Tensor product        | 1,250         | 18,750                 | 15    |
| Angular build time    | 0.002s        | 0.4s                   | 200   |

### Wall time estimates

Level 4 2D solver at spectral dimensions: ~0.5-2s per PES point (from
composed_diatomic.py benchmarks, where 2D solver is one component).

4-electron 2D solver at l_max=2:
- Matrix assembly: 80 quadrature points x 0.4s = **32s**
- Iterative eigensolve: 100 iterations x 18,750 matvec = **~10s**
- Total per PES point: **~40-60s**
- Full PES scan (12 points): **~8-12 min**

Compare to current adiabatic solver: 69s/point at l_max=2 (Track AM, FD grid).

**The 2D solver would be COMPARABLE in cost to the current adiabatic FD solver.**
This is because the adiabatic solver does a dense eigensolve of dim=2592 at
each of ~130 rho-points (130 x 2592^3/3 operations ~ 2.3T flops), while the
2D solver does one sparse iterative solve of dim=18,750.

At l_max=3:
- Matrix assembly: 80 x 2s (estimated from 750->2500 scaling) = **160s**
- Iterative eigensolve: 200 iterations x 62,500 matvec = **~150s**
- Total per PES point: **~5-10 min**
- Full PES scan: **~1-2 hours**

Compare to adiabatic at l_max=3: ~15 min/point (Track AM).

---

## (e) PARTIAL 2D: Which Hyperangle(s) to Prioritize

### Recommendation: alpha2 (inter-pair angle)

**alpha2 is the clear priority** for a partial 2D approach. Evidence:

1. **Physical role:** alpha2 controls core/valence partition. As R_e changes,
   the core electrons contract and valence electrons expand. alpha2 tracks
   this differential scaling directly. This is the dominant source of
   non-adiabatic coupling (the adiabatic approximation breaks down when
   the core/valence partition changes rapidly with R_e).

2. **Analogy to Level 4:** At Level 4, the single hyperangle alpha controls
   the inner/outer electron partition, and the 2D solver treating (R_e, alpha)
   variationally reduces l_max drift by 4x. The 4-electron alpha2 plays the
   same role.

3. **alpha1, alpha3 are intra-pair:** They control how each pair's two
   electrons distribute relative to each other. This is angular correlation
   within a pair, which is less R_e-sensitive.

4. **Coupling structure:** The nuclear potential V_nuc depends on rho = R/(2*R_e)
   through min(s_i, rho)/max(s_i, rho). Since s_i = f(alpha1, alpha2, alpha3),
   and rho depends on R_e, the R_e-coupling enters all three alphas. But the
   dominant contribution is through alpha2, which sets the overall scale of
   s_1,s_2 (core) vs s_3,s_4 (valence). The alpha1, alpha3 dependence is
   secondary (redistribution within a fixed-scale pair).

### Partial 2D dimensions

Treat (R_e, alpha2) as tensor product; keep (alpha1, alpha3, channels) adiabatic.

At each alpha2 point, solve a 2D adiabatic problem in (alpha1, alpha3, channels)
to get eigenvalues mu(alpha2; R_e). Then solve the (R_e, alpha2) problem with
these eigenvalues as the effective potential.

| l_max | alpha2 basis | Radial basis | Tensor product | Feasible? |
|-------|-------------|--------------|----------------|-----------|
| 2     | 5 x 6 ch = 30 | 25         | 750            | TRIVIAL   |
| 3     | 5 x 20 ch = 100 | 25       | 2,500          | TRIVIAL   |
| 4     | 5 x 50 ch = 250 | 25       | 6,250          | EASY      |

Wait -- this doesn't account for the (alpha1, alpha3) adiabatic solve. The
partial 2D approach would:

1. For each (R_e, alpha2) point: solve the (alpha1, alpha3) eigenvalue
   problem in the channel basis. Dimension: n_ch x n_basis_1 x n_basis_3
   = 6 x 5 x 5 = 150 at l_max=2.
2. Take the lowest eigenvalue mu(R_e, alpha2).
3. Assemble the (R_e, alpha2) Hamiltonian and find the ground state.

Step 1 costs: 150^3/3 ~ 1.1M flops per (R_e, alpha2) point.
Grid: 25 R_e x 30 alpha2 = 750 points. Total: 0.8G flops ~ 0.5s.
Step 3: 750 x 750 eigensolve ~ 0.14G flops ~ 0.1s.
**Total: ~1s per PES point.** 700x faster than current adiabatic solver.

### Alternative partial 2D: (R_e, alpha2) with alpha1, alpha3 channels

A cleaner approach: expand alpha1 and alpha3 in a Jacobi basis (5 functions each),
making them part of the "channel" label. The 2D solver then works in
(R_e, alpha2) with an expanded channel set:

n_effective_channels = n_S4_channels x n_basis_1 x n_basis_3 = 6 x 5 x 5 = 150
alpha2 basis: 5 functions
Radial basis: 25 functions

2D tensor product: 150 x 5 x 25 = 18,750

This is the same as the full 2D! The partial approach only saves if we
use an adiabatic approximation for alpha1, alpha3 (reducing 150 to ~1-3
effective channels by taking the lowest adiabatic eigenvalues).

**Corrected partial 2D dimensions:**

| l_max | Adiabatic alpha1,3 channels | alpha2 basis | Radial | Tensor product |
|-------|----------------------------|-------------|--------|----------------|
| 2     | 3 (lowest eigenstates)      | 30          | 25     | 2,250          |
| 3     | 3                           | 100         | 25     | 7,500          |
| 4     | 3                           | 250         | 25     | 18,750         |

This is a meaningful reduction: 18,750 -> 2,250 (8x) at l_max=2.

---

## DIMENSION TABLE (Summary)

| l_max | Angular dim | Radial dim | Full 2D tensor | Partial 2D tensor | Iterative? | Est. time (full) | Est. time (partial) |
|-------|-------------|------------|----------------|-------------------|------------|------------------|---------------------|
| 0     | 125         | 25         | 3,125          | 375               | YES        | ~1s              | <1s                 |
| 1     | 125         | 25         | 3,125          | 375               | YES        | ~1s              | <1s                 |
| 2     | 750         | 25         | 18,750         | 2,250             | YES        | ~40-60s          | ~5-10s              |
| 3     | 2,500       | 25         | 62,500         | 7,500             | YES        | ~5-10 min        | ~30-60s             |
| 4     | 6,250       | 25         | 156,250        | 18,750            | BORDERLINE | ~1-5 hrs         | ~5-10 min           |

---

## FEASIBILITY VERDICT

### Full 2D solver

| l_max | Verdict                          | Rationale                                   |
|-------|----------------------------------|---------------------------------------------|
| 0-1   | FEASIBLE                         | Tensor product < 3,200; trivial eigensolve  |
| 2     | FEASIBLE WITH OPTIMIZATION       | 18,750 sparse iterative; ~1 min/point       |
| 3     | REQUIRES ALGORITHMIC BREAKTHROUGH| 62,500 dim; ~10 min/point; matrix assembly dominant |
| 4     | INTRACTABLE                      | 156,250 dim; memory prohibitive for dense   |

### Partial 2D solver (recommended)

| l_max | Verdict                    | Rationale                                        |
|-------|----------------------------|--------------------------------------------------|
| 0-2   | FEASIBLE                   | Tensor product < 2,300; fast eigensolve          |
| 3     | FEASIBLE WITH OPTIMIZATION | 7,500 dim; ~1 min/point; practical for PES scans |
| 4     | FEASIBLE WITH OPTIMIZATION | 18,750 dim; ~10 min/point; same as full l_max=2  |

### Key insight

The partial 2D approach at l_max=4 has the SAME dimension as the full 2D
at l_max=2. This makes it the clear recommendation: it extends the feasible
l_max by 2 levels compared to the full 2D approach.

---

## RECOMMENDATIONS

1. **Build partial 2D solver treating (R_e, alpha2) variationally.**
   Keep alpha1, alpha3 adiabatic (take lowest 1-3 eigenstates).

2. **alpha2 is the priority variable** -- it controls the core/valence
   partition and has the strongest R_e-coupling.

3. **Expected impact:** At Level 4, the 2D solver reduced l_max drift
   by 4x (from +0.400 to +0.100 bohr/l_max). A similar improvement
   at Level 4N could shift R_eq from ~1.0 toward ~1.5-2.0 bohr at
   l_max=2, though convergence to 3.015 would still require higher l_max.

4. **Implementation path:**
   - Step 1: Build alpha2-adiabatic solver (alpha1, alpha3 eigenproblem
     at fixed alpha2). This is a 2D eigenproblem of dim ~150 at l_max=2.
   - Step 2: Assemble (R_e, alpha2) tensor product Hamiltonian.
   - Step 3: Solve with iterative eigensolver (ARPACK/Lanczos).
   - Test at l_max=2 first; if R_eq improves, extend to l_max=3.

5. **Risk:** The adiabatic approximation may not be the dominant error
   source at Level 4N. Track AM showed R_eq stuck at ~1.0 for both
   l_max=2 and l_max=3, suggesting the bottleneck may be channel
   truncation rather than adiabatic error. The 2D solver would test
   this hypothesis.
