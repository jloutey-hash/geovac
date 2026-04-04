# Track AG: Full 4-Electron Mol-Frame Hyperspherical Solver for LiH

## Scoping Analysis

**Date:** March 31, 2026
**Track:** AG (v2.0.19 sprint)
**Goal:** Scope what a full 4-electron mol-frame hyperspherical solver for LiH would look like, with no PK, no Z_eff, no composed geometry.

**Verdict: INTRACTABLE at production l_max. POTENTIALLY TRACTABLE at l_max=1 for proof-of-concept.**

---

## (a) Coordinates

For N=4 electrons + 2 fixed nuclei (Born-Oppenheimer), the electronic configuration space is R^12.

**Coordinate decomposition:**

| Quantity | N=2 (Level 4) | N=4 (Full LiH) |
|----------|:---:|:---:|
| Config dim | 6 | 12 |
| Hyperradius R_e | 1 | 1 |
| Angular manifold | S^5 | S^11 |
| Angular dim | 5 | 11 |
| Isometry group | SO(6) | SO(12) |
| Casimir formula | nu(nu+4)/2 | nu(nu+10)/2 |
| Centrifugal term | 15/(8R_e^2) | 99/(8R_e^2) |
| Hyperangles (alpha) | 1 | 3 |
| Direction angles | 4 (theta_1, phi_1, theta_2, phi_2) | 8 (theta_i, phi_i for i=1..4) |

The natural decomposition uses a Jacobi tree:
- R_e = sqrt(r_1^2 + r_2^2 + r_3^2 + r_4^2) -- hyperradius
- alpha_1 = arctan(r_2/r_1) -- pair (1,2) correlation
- alpha_2 = arctan(rho_34/rho_12) -- inter-pair, where rho_ij = sqrt(r_i^2 + r_j^2)
- alpha_3 = arctan(r_4/r_3) -- pair (3,4) correlation
- 8 direction angles (theta_i, phi_i) for 4 electrons

The 3 hyperangles alpha_1, alpha_2, alpha_3 replace the single alpha of Level 4.

---

## (b) Angular Hilbert Space

### Channel count table

Channels are labeled by (l_1, m_1, l_2, m_2, l_3, m_3, l_4, m_4) with M_total = sum(m_i) = 0.

| l_max | L4 homo | L4 hetero | 4e sigma | 4e full (M=0) | 4e L=0 coupled (atomic) |
|:---:|:---:|:---:|:---:|:---:|:---:|
| 0 | 1 | 1 | 1 | 1 | 1 |
| 1 | 2 | 4 | 16 | 70 | 14 |
| 2 | 5 | 9 | 81 | 1,107 | 91 |
| 3 | 8 | 16 | 256 | (est. ~8,000) | (est. ~350) |
| 4 | 13 | 25 | 625 | (est. ~40,000) | (est. ~1,000) |

Key observations:
- **4-electron sigma channels scale as (l_max+1)^4** (vs (l_max+1)^2 for 2-electron). This is the fundamental dimensionality explosion.
- **Full M=0 channels are far larger:** at l_max=2, 1,107 channels vs 81 sigma-only. The m-coupling adds a factor of ~14x.
- **The L=0 atomic count** (relevant for understanding how many coupled angular momenta configurations give total L=0) is 91 at l_max=2.

### Spectral basis dimension

The Level 4 angular eigenproblem has 1 continuous variable (alpha) discretized with n_basis Jacobi polynomials per channel.

For 4 electrons, the angular eigenproblem has 3 continuous variables (alpha_1, alpha_2, alpha_3), requiring a tensor-product spectral basis with n_basis^3 functions per channel.

| l_max | L4 dim (n_b=10) | 4e sigma dim (after S_4, n_b=10) | Ratio |
|:---:|:---:|:---:|:---:|
| 0 | 10 | 1,000 | 100x |
| 1 | 40 | 2,000 | 50x |
| 2 | 90 | 13,000 | 144x |
| 3 | 160 | 42,000 | 262x |
| 4 | 250 | 104,000 | 416x |

The 4-electron angular matrix is 100-400x larger than the Level 4 matrix at the same l_max.

---

## (c) Symmetry Reduction

### SO(3) spatial rotation

Total angular momentum L=0 restriction. For the molecular problem (axial symmetry), only M_total = 0 is enforced (not full L=0). The M=0 constraint already reduces the channel count from (2l+1)^4 per l-set to a constrained sum.

### S_4 electron exchange

| S_4 irrep | Young diagram | Dimension | Role |
|-----------|:---:|:---:|------|
| [4] | fully symmetric | 1 | Paired with fully antisymmetric spin (S=2 quintet) |
| [3,1] | standard | 3 | Paired with [2,1,1] spin (S=1 triplet) |
| [2,2] | two-row | 2 | **LiH ground state** (S=0 singlet) |
| [2,1,1] | standard conjugate | 3 | Paired with [3,1] spin |
| [1,1,1,1] | fully antisymmetric | 1 | Paired with fully symmetric spin |

**LiH ground state (1Sigma+, S=0):** The spatial wavefunction transforms under the [2,2] irrep of S_4.

Burnside's lemma gives dim^2/|S_4| = 4/24 = 1/6 as the approximate fraction of channels surviving antisymmetrization. This means roughly **1/6 of the uncoupled M=0 channels survive the S_4 projection**.

### Combined reduction factors

| Symmetry | Reduction factor |
|----------|:---:|
| M_total = 0 | ~1/(l_max+1)^2 from full product |
| S_4 [2,2] projection | ~1/6 |
| L_total = 0 (atomic, partial) | variable |
| Inversion (homonuclear only) | 1/2 (not applicable for LiH) |

Even with all reductions, the 4-electron channel space is orders of magnitude larger than Level 4.

---

## (d) Separation Structure

**The rho-collapse DOES generalize to 4 electrons.**

At Level 4, the angular eigenvalue depends on (R, R_e) only through rho = R/(2R_e), collapsing a 2D parameter sweep to 1D.

For 4 electrons, each nuclear distance is:
  r_{iA} = |r_i - R_A| = R_e * h_iA(Omega, R/R_e)

where h_iA depends on the angular coordinates Omega and the ratio R/R_e = 2*rho only. Since R enters exclusively as rho, the angular eigenvalue problem at fixed R_e depends on (R, R_e) only through rho = R/(2R_e), exactly as at Level 4.

**V_ee is still rho-independent** (6 electron pairs, all involving only angular coordinates).

**V_nuc is rho-dependent** (8 terms: 4 electrons x 2 nuclei, each with split-region Legendre structure).

The rho-grid (~100-200 points) remains the same as Level 4. This is a significant positive finding: the parameter dimensionality does not increase.

---

## (e) Spectral Basis: SO(12) Casimir

| nu | mu_free | l_eff | Degeneracy |
|:---:|:---:|:---:|:---:|
| 0 | 0.0 | 4.0 | 1 |
| 1 | 5.5 | 5.0 | 12 |
| **2** | **12.0** | **6.0** | **77** |
| 3 | 19.5 | 7.0 | 352 |
| 4 | 28.0 | 8.0 | 1,287 |
| 5 | 37.5 | 9.0 | 4,004 |
| 6 | 48.0 | 10.0 | 11,011 |

**Key result from Paper 16:** For 4-electron singlet ([2,2] spatial irrep), the minimum grand angular momentum is nu_min = N - lambda_1 = 4 - 2 = 2, giving mu_free = 12.0. States with nu = 0 and nu = 1 are Pauli-forbidden for the singlet.

**Comparison:** At Level 4 (SO(6)), the singlet ground state has nu_min = 0 with mu_free = 0. The 4-electron problem starts at a much higher angular momentum baseline.

The SO(12) degeneracy grows explosively: 77 at nu=2, 1,287 at nu=4, 11,011 at nu=6. This is the group-theoretic origin of the Hilbert space explosion.

---

## (f) Coupling Matrices

### V_ee (6 electron pairs)

C(4,2) = 6 pairs: (1,2), (1,3), (1,4), (2,3), (2,4), (3,4).

Each pair 1/r_{ij} has a multipole expansion involving Gaunt integrals. The coupling structure:
- **Pair (1,2):** Couples (l_1, l_2) -> (l_1', l_2') via Gaunt(l_1', k, l_1) * Gaunt(l_2', k, l_2). Leaves l_3, l_4 unchanged. This is the same structure as Level 4 V_ee.
- **Pair (3,4):** Analogous, couples (l_3, l_4). Leaves l_1, l_2 unchanged.
- **Cross pairs (1,3), (1,4), (2,3), (2,4):** Couple across the Jacobi tree. For example, pair (1,3) couples l_1 and l_3 simultaneously while leaving l_2 and l_4 unchanged. These cross-pair couplings involve different hyperangles and create a MUCH denser coupling matrix.

**The cross-pair V_ee couplings are the structurally new element.** At Level 4, V_ee couples only 2 angular momenta via one Gaunt product. At 4 electrons, each pair couples 2 angular momenta, but the 4 "cross" pairs entangle angular momenta across the Jacobi tree, preventing factorization of the coupling matrix.

Gaunt integrals still apply for each individual pair, and the selection rules still enforce sparsity (triangle inequality on each pair). But the COMBINED sparsity is much lower because different pairs constrain different angular momentum indices.

### V_nuc (4 electrons x 2 nuclei = 8 terms)

Each nuclear attraction -Z_X/|r_i - R_X| has the SAME split-region Legendre form as Level 4:

  -Z_X * sum_k (min(s_i, rho_X)/max(s_i, rho_X))^k * P_k(cos theta_iX) / max(s_i, rho_X)

The nuclear coupling is diagonal in the angular momenta of all electrons EXCEPT electron i. This means each nuclear term couples only (l_i, l_i') while leaving the other three l-values unchanged.

**Nuclear coupling structure is the same as Level 4** -- just replicated 4 times (one per electron). This is manageable.

---

## (g) Feasibility Estimate

### Wall-time estimates (sigma-only, n_basis=10 per alpha)

| l_max | 4e angular dim | Time/rho-point | Total PES time | L4 total | Ratio |
|:---:|:---:|:---:|:---:|:---:|:---:|
| 0 | 1,000 | 1.0 s | 130 s | 0.13 ms | 10^6 x |
| 1 | 2,000 | 8.0 s | 17 min | 8.3 ms | 1.2 x 10^5 x |
| 2 | 13,000 | 37 min | 3.3 days | 95 ms | 3.0 x 10^6 x |
| 3 | 42,000 | 21 hr | 111 days | 0.53 s | 1.8 x 10^7 x |
| 4 | 104,000 | 13 days | 4.6 years | 2.0 s | 7.2 x 10^7 x |

**Notes on the estimate:**
- Times assume 1 GFLOP/s effective eigensolve throughput (conservative for modern hardware).
- The 4e angular problem has 3 hyperangles (tensor product n_basis^3 = 1000 per channel), vs 1 hyperangle (n_basis = 10 per channel) for Level 4.
- S_4 antisymmetry reduction (~1/6) is included.
- Full M-coupled channels would be ~14x worse than sigma-only.
- Only the LOWEST eigenvalue is needed (Lanczos/Davidson could reduce cost from O(N^3) to O(N^2) for sparse matrices), but the matrix is dense (potential coupling is all-to-all).

### Assessment by l_max

- **l_max=0:** 130 seconds, trivially tractable. But l_max=0 captures no angular correlation -- this gives only the s-wave result.
- **l_max=1:** 17 minutes, tractable for proof-of-concept. This captures the basic core-valence structure. 16 sigma channels (2 antisymmetrized after S_4 projection).
- **l_max=2:** 3.3 days for a single PES point. A 20-point PES scan would take ~2 months. Tractable on a cluster but not practical for development.
- **l_max=3+:** Intractable.

### Comparison to Level 4

The Level 4 solver at l_max=6 (61 channels, 754 s/point) is well within workstation capability. The 4-electron solver at l_max=2 would require ~3 x 10^6 times more computation.

---

## (h) Antisymmetry

### How it manifests in hyperspherical coordinates

At Level 4 (2 electrons), exchange symmetry (12) acts as alpha -> pi/2 - alpha, producing the gerade constraint l_1 + l_2 = even.

For 4 electrons with S_4 permutation symmetry, the action on the Jacobi tree is:

| Permutation | Action on hyperangles | Action on angular momenta |
|-------------|----------------------|--------------------------|
| (12) | alpha_1 -> pi/2 - alpha_1 | swap (l_1, m_1) <-> (l_2, m_2) |
| (34) | alpha_3 -> pi/2 - alpha_3 | swap (l_3, m_3) <-> (l_4, m_4) |
| (13) | Mixes ALL three alpha's | swap (l_1, m_1) <-> (l_3, m_3) |
| (24) | Mixes ALL three alpha's | swap (l_2, m_2) <-> (l_4, m_4) |
| (1234) | Complex rearrangement | Cyclic permutation |

**Critical observation:** The transpositions (12) and (34) are "within-pair" and act simply on individual hyperangles. But transpositions that cross the Jacobi tree ((13), (14), (23), (24)) mix ALL three hyperangles simultaneously and entangle all angular momenta.

The [2,2] projection operator is:

  P_{[2,2]} = (2/24) * sum_{g in S_4} chi_{[2,2]}(g) * D(g)

where D(g) is the representation matrix of g on the channel basis. This requires computing all 24 representation matrices, but is a one-time setup cost.

### The gerade analog

For 2 electrons: l_1 + l_2 = even (single constraint from (12) exchange).

For 4 electrons with [2,2] symmetry:
- l_1 + l_2 = even (from (12) exchange)
- l_3 + l_4 = even (from (34) exchange)
- A MIXED constraint from cross-pair transpositions that CANNOT be decomposed into independent constraints on each pair.

This mixed constraint is the fundamental difference from composed geometry, which treats each pair independently.

---

## (i) Pauli Exclusion and Core Separation

**Answer: PARTIAL.** Antisymmetry provides the CONSTRAINT but not the SEPARATION.

The [2,2] projector ensures no two electrons can occupy the same single-particle state. The 1s^2 core is correctly excluded from the valence space at the level of the wavefunction.

However, the 4-electron hyperspherical solver does NOT separate core and valence into independent subspaces. All 4 electrons live in the SAME 12-dimensional configuration space. The core electrons contribute to the full coupled-channel problem.

**This is both the advantage and the cost:**

| Aspect | Full 4e solver | Composed geometry (Level 5) |
|--------|:-:|:-:|
| Core-valence correlation | Exact | Approximated by Z_eff |
| Core polarization | Included | Missing |
| PK pseudopotential | Not needed | Required (empirical R_ref) |
| Antisymmetry | Exact (S_4) | Approximate (orthogonality constraint) |
| Angular dimension | ~13,000 at l_max=2 | ~90 at l_max=2 |
| Wall time | ~3 days/point | ~1 s/point |

**The Pauli centrifugal barrier:** From Paper 16, nu_min = 2 gives mu_free = 12.0. This centrifugal term in the hyperradial equation acts as a "Pauli barrier" that keeps the 4 electrons from collapsing to the same spatial region. It is the group-theoretic encoding of the exclusion principle.

**Would the 1s^2 core automatically separate?** Not in the solver, but in the WAVEFUNCTION. At large R_e (electron cloud much larger than nucleus), the wavefunction would naturally develop a nodal structure with two tight (core) and two diffuse (valence) electrons. This separation would be an OUTPUT of the solver, not an INPUT -- unlike composed geometry where it is assumed.

---

## Structural Findings

### What the full 4e solver would gain over composed geometry

1. **Exact core-valence correlation.** No PK pseudopotential, no Z_eff approximation. The l_max divergence problem (CLAUDE.md Section 3) cannot occur because there is no channel-blind PK.

2. **Exact antisymmetry.** The [2,2] projector handles all 24 permutations correctly, including the cross-pair exchanges that composed geometry cannot represent.

3. **No empirical parameters.** No R_ref, no PK width, no Z_eff screening function. The solver is fully ab initio.

4. **Core polarization.** The core electrons respond to the valence electrons, an effect absent in composed geometry.

### What makes it intractable

1. **Dimensional explosion.** Going from 2 to 4 electrons changes the angular space from S^5 (5D) to S^11 (11D). The channel count scales as (l_max+1)^4 instead of (l_max+1)^2.

2. **Three hyperangles instead of one.** The spectral basis has n_basis^3 functions per channel instead of n_basis. This alone accounts for a 1000x increase at n_basis=10.

3. **Cross-pair cusps.** At Level 4, the e-e cusp is at alpha=pi/4, theta_{12}=0 -- a 2D submanifold in the 5D angular space. For 4 electrons, the 6 pairwise cusps include 4 "cross-pair" cusps (e.g., r_{13}=0) that live in complicated submanifolds cutting across the Jacobi tree. These cannot be resolved by simple boundary conditions on individual hyperangles.

4. **Dense V_ee coupling.** The 6 electron pairs create coupling between all 4 angular momentum indices simultaneously, preventing the sparse structure that makes Level 4 efficient.

### The fundamental bottleneck: angular dimensionality

The bottleneck is NOT the rho-sweep (which is still 1D) or the hyperradial solve (which is still a 1D ODE). It is the ANGULAR eigensolve at each rho-point.

At Level 4, the angular eigensolve is a ~50x50 matrix (spectral basis, l_max=2). At 4 electrons, it would be a ~13,000 x 13,000 matrix. The O(N^3) eigensolve dominates the cost.

---

## Feasibility Verdict

### INTRACTABLE at production l_max (l_max >= 2)

The angular dimension of ~13,000 at l_max=2 (sigma-only, after S_4 reduction) requires ~3 days per PES point. A production PES scan with 20 R-points would take ~2 months. With full M-coupling (1,107 raw channels), the cost increases another ~14x.

At l_max=4 (which is needed for 96% D_e at Level 4), the angular dimension exceeds 100,000 and the computation is completely infeasible.

### POTENTIALLY TRACTABLE for proof-of-concept at l_max=1

At l_max=1 (16 sigma channels, ~2 after S_4), the angular dimension is ~2,000 and a single PES point takes ~17 minutes. A coarse 20-point PES scan would take ~6 hours. This is feasible for a proof-of-concept that demonstrates:
- Correct antisymmetry structure
- Automatic core-valence separation in the wavefunction
- PK-free equilibrium (or lack thereof)
- Comparison with composed geometry at the same l_max

### Recommendations

1. **Do NOT implement a production 4-electron solver.** The dimensional scaling makes it non-competitive with composed geometry for routine calculations.

2. **Consider a proof-of-concept at l_max=0-1** to validate the S_4 projection machinery and test whether the 1s^2 core separates automatically. This would be a valuable diagnostic (like Track AD) rather than a production tool.

3. **The composed geometry approach is vindicated.** The PK pseudopotential and Z_eff screening, while approximate, reduce the problem from ~13,000 to ~90 angular dimensions at l_max=2 -- a factor of 144x. The l_max divergence and empirical R_ref are prices worth paying for this compression.

4. **The cross-pair cusps at 4 electrons** represent a qualitatively new challenge that has no analog at Level 4. Any future attempt at a full N-electron solver would need to address these.

5. **Sparse iterative eigensolvers** (Lanczos, Davidson) could reduce the O(N^3) dense eigensolve to O(N^2 * k) for k lowest eigenvalues. This would help by ~1-2 orders of magnitude but does not change the fundamental scaling.

6. **Selected CI / importance truncation** could reduce the channel space by keeping only the most important angular configurations. This is the standard approach in traditional quantum chemistry (CIPSI, heat-bath CI) and could potentially make l_max=2 feasible, but would require significant algorithmic development.

---

## Code Artifacts

- `geovac/n_electron_scope.py` -- Scoping module with dimension counting, symmetry analysis, and feasibility estimates
- `tests/test_n_electron_scope.py` -- 37 tests validating all dimension counts and symmetry properties
- `debug/track_ag/analysis.md` -- This analysis document
