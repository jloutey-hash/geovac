# Feasibility Assessment: Hyperspherical Lattice for Helium

**Date:** 2026-03-14
**Status:** Phase 1 — Mathematical Reconnaissance
**Target:** He ground state, sub-1% error (exact: -2.9037 Ha)
**Current best:** -2.8508 Ha (1.82% error, S³ lattice with variational Z_eff)

---

## 1. Can the Prolate Root-Finding Strategy Be Reused?

### Prolate Strategy (Paper 11)
1. **Angular solver:** Legendre spectral expansion in η → tridiagonal eigenvalue
   problem → separation constant A(c²)
2. **Radial solver:** Self-adjoint FD in ξ → tridiagonal eigenvalue problem →
   λ_max(c², A)
3. **Root-finding:** Brent's method on f(c²) = λ_max = 0

### Hyperspherical Analogue

| Step | Prolate (H₂⁺) | Hyperspherical (He) | Compatible? |
|:---|:---|:---|:---:|
| Angular basis | P_l^m(η) | {φ_l(α)} coupled by V_ee | **YES** (matrix eigenvalue) |
| Angular solve | Tridiagonal EVP | Block-tridiagonal EVP | **YES** (slightly larger) |
| Angular output | A (scalar) | U_μ(R) (curve per channel) | **SIMILAR** |
| Radial grid | Uniform FD in ξ | Uniform FD in R | **YES** (identical) |
| Radial solve | Tridiagonal EVP | Coupled-channel FD | **MODIFIED** |
| Root-finding | Brent on c² | Direct eigenvalue of coupled radial | **DIFFERENT** |

**Verdict:** The overall strategy is **compatible** but requires two modifications:

1. **The angular problem is a coupled matrix eigenvalue problem** (not a single
   tridiagonal). For l_max = 3 and L=0, this is a 4×4 matrix at each R —
   trivially solvable.

2. **The radial problem includes inter-channel coupling.** In the adiabatic
   approximation (single channel), it reduces to exactly the prolate strategy.
   For coupled channels, it becomes a matrix Schrödinger equation — standard
   in nuclear physics, well-understood numerically.

---

## 2. How Many Coupled Channels for 1% Accuracy?

### Literature Benchmarks (Macek 1968, Lin 1995)

| Channels | Method | E (Ha) | Error |
|:---:|:---|:---:|:---:|
| 1 | Adiabatic (lowest ¹S curve only) | -2.879 | 0.85% |
| 2 | Coupled (lowest 2 ¹S curves) | -2.900 | 0.12% |
| 3 | Coupled (3 ¹S curves) | -2.9033 | 0.014% |
| 5 | Coupled (5 ¹S curves) | -2.90370 | 0.001% |
| Exact | Pekeris (1958), Drake (1999) | -2.90372 | 0 |

**Key finding:** **2 coupled channels give 0.12% error** — already sub-1%.
Even the crude single-channel adiabatic approximation (0.85%) beats our
current S³ result (1.82%).

### Why So Few Channels Suffice

The He ground state is dominated by the s-wave (l₁ = l₂ = 0) configuration.
The first channel (ν=1) is approximately the 1s² configuration. The second
channel (ν=2) mixes in the 2p² configuration (l₁ = l₂ = 1), which captures
the angular correlation.

Higher channels (d-wave: l₁ = l₂ = 2, etc.) contribute < 0.01 Ha each.

---

## 3. Expected Graph Structure

### Nodes

```
Node set: V = {(R_i, ν) : i = 1,...,N_R  and  ν = 1,...,N_ch}

Total nodes: N = N_R × N_ch
```

**Estimated sizes:**
- N_R ≈ 200–500 (hyperradial grid points, R ∈ [0, R_max ≈ 20 bohr])
- N_ch = 2 (for 0.12% accuracy) to 5 (for 0.001%)
- **Total: 400–2500 nodes** — extremely small!

For comparison:
- S³ lattice at nmax=5: ~50 nodes
- Prolate lattice: ~8000 nodes (N_ξ)
- Full CI for He at nmax=5: 215,820 Slater determinants

### Edges

```
Type 1: Radial FD (within channel)
  (R_i, ν) ↔ (R_{i±1}, ν)
  Count: ~2 × N_R × N_ch
  Bandwidth: 1 (tridiagonal)

Type 2: Non-adiabatic coupling (between channels, same R)
  (R_i, ν) ↔ (R_i, μ)  for ν ≠ μ
  Count: N_R × N_ch × (N_ch - 1)
  Bandwidth: N_ch (block structure)

Total edges: ~N_R × (2N_ch + N_ch²) ≈ N_R × N_ch²
For N_R=300, N_ch=2: ~1200 edges → very sparse graph
```

### Matrix Structure

The coupled-channel Hamiltonian is a **block tridiagonal matrix**:

```
H = [ H₁₁  H₁₂  0    0    ...  ]    H_νν = -½ d²/dR² + V_ν(R)
    [ H₂₁  H₂₂  H₂₃  0    ...  ]    H_νμ = -P_νμ d/dR - ½Q_νμ
    [ 0     H₃₂  H₃₃  H₃₄  ...  ]
    [ ...                         ]

Block size: N_ch × N_ch
Number of blocks: N_R
Total matrix: (N_R · N_ch) × (N_R · N_ch) ≈ 600 × 600
```

This is a **standard sparse eigenvalue problem** — solvable in milliseconds.

---

## 4. Is There a Natural Kinetic Scale?

### The Question

On S³, the universal kinetic scale κ = -1/16 maps the graph Laplacian
eigenvalues to hydrogen energies. Does a similar constant exist for the
hyperradial problem?

### Analysis

**S³ (hydrogen):**
```
Graph eigenvalue: λ_n = -(n² - 1)
Physical energy:  E_n = -1/(2n²)
Mapping:          E_n = κ · λ_n / n²   with κ = -1/16 ← NOT QUITE

More precisely: E_n = -Z²/(2n²) follows from the Fock projection
p₀² = -2E → p₀ = Z/n, which is the energy-shell constraint.
The -1/16 maps between the GRAPH Laplacian and the S³ Laplace-Beltrami.
```

**S⁵ (helium):**
```
Hyperangular eigenvalue: Λ²_K = -K(K+4) on unit S⁵
Gegenbauer polynomial:   C²_K(cos χ) with eigenvalue K(K+4)
```

The S⁵ eigenvalues follow the pattern K(K+4) = (K+2)² - 4, which is the
5D analogue of the S³ pattern n(n+2) = (n+1)² - 1 (with the shift from
dimension d: l(l+d-2) for the Laplace-Beltrami on S^{d-1}).

**For a graph Laplacian discretizing S⁵:**
```
Candidate kinetic scale: κ₅ = ?

The S³ result κ = -1/16 relates to the conformal dimension:
  In d=3: κ = -1/(4d) = -1/12?  No, it's -1/16.
  Actually κ = -1/(2d(d-1))? No.

From Paper 0: κ = -1/16 = -1/(4·4) = -1/(4·(d+1)) for d=3?
  For d=5: κ₅ = -1/(4·6) = -1/24?

This is SPECULATIVE. The kinetic scale for S⁵ must be derived from a
Fock-type projection in 6D, which has not been done.
```

**Verdict:** There is likely a natural kinetic scale for S⁵, but it has
not been established. This is a **Phase 2 question**. For the initial solver,
the hyperradial approach does not need it — the angular problem is solved
numerically (not via a graph Laplacian), and the radial problem uses standard
FD.

The kinetic scale question becomes relevant only if we want to build a
**discrete S⁵ graph** analogous to the S³ paraboloid lattice. For the
adiabatic hyperspherical solver, it is not needed.

---

## 5. Computational Cost Estimate

### Angular Solve (at each R)

For L=0, l_max channels:
- Build coupling matrix: O(l_max²) Gaunt integrals (precomputed)
- α-grid integration: O(N_α × l_max²) for the matrix elements
- Eigenvalue decomposition: O(l_max³)
- Total per R point: O(N_α × l_max²)
- N_α ≈ 50 (Gauss-Legendre quadrature in α)
- l_max = 3 (for 0.01% accuracy)

**Cost per R: ~O(500) operations — negligible**

### Radial Solve

- Build N_ch coupled FD equations on N_R grid
- Solve (N_R × N_ch) × (N_R × N_ch) sparse eigenvalue problem
- N_R ≈ 300, N_ch ≈ 2–5
- Matrix size: 600–1500

**Cost: O(N_R × N_ch²) ≈ O(1000) — negligible**

### Total

The bottleneck is computing U_μ(R) at each R grid point:
- N_R angular eigenvalue problems, each O(N_α × l_max²)
- Total: O(N_R × N_α × l_max²) ≈ O(300 × 50 × 9) ≈ 135,000 operations

**Expected runtime: < 1 second.** This is vastly faster than Full CI.

---

## 6. Comparison with Existing GeoVac Approaches

| Method | He Energy (Ha) | Error | Runtime | Basis Size |
|:---|:---:|:---:|:---:|:---:|
| S³ lattice (nmax=5, Z_eff) | -2.8508 | 1.82% | ~0.1s | 50 nodes |
| Full CI (nmax=5, Slater F⁰) | ~-2.854 | 1.71% | ~120s | 215k SDs |
| Hylleraas (9 terms, prolate) | ~94.7% D_e | — | ~10s | 9 basis fns |
| **Hyperspherical (2 ch.)** | **-2.900** | **0.12%** | **<1s** | **600 nodes** |
| **Hyperspherical (3 ch.)** | **-2.9033** | **0.014%** | **<1s** | **900 nodes** |

The hyperspherical approach is projected to be both **more accurate** and
**faster** than every existing GeoVac method for He.

---

## 7. Risk Assessment

### Low Risk
- **Mathematical framework:** Well-established (Macek 1968, Lin 1995, thousands of papers)
- **Numerical methods:** Standard (FD + spectral + eigenvalue)
- **Channel convergence:** Literature confirms 2 channels suffice for sub-1%
- **Code reuse:** Angular solver is a Legendre spectral problem (same as prolate η)

### Medium Risk
- **Non-adiabatic coupling computation:** Requires numerical differentiation of
  channel functions with respect to R. Must use smooth interpolation or
  analytic derivatives. Could introduce numerical noise at avoided crossings.
- **Cusp boundary condition:** Enforcing the Kato cusp at (α=π/4, θ₁₂=0)
  in the angular solver may require special basis functions or cusp factors.
- **R → 0 behavior:** The Fock logarithmic terms R^{1/2} and R ln R appear
  at very small R. Standard FD may need logarithmic grid spacing near R=0.

### Low-Medium Risk
- **Connection to S³ lattice:** The hyperspherical approach works in physical
  space, not momentum space. Connecting it to Fock's S³ projection is
  mathematically interesting but not necessary for the solver.

---

## 8. Recommended Phase 2 Plan

### Phase 2A: Single-Channel Adiabatic Solver
1. Implement charge function C(α, θ₁₂) for Z=2
2. L=0 partial-wave reduction: coupled equations in α
3. Solve angular eigenvalue problem at a grid of R values → U₁(R)
4. Solve hyperradial equation with U₁(R) potential → E₁
5. **Validation target:** E ≈ -2.879 Ha (0.85% error)

### Phase 2B: Two-Channel Coupled Solver
1. Extract second adiabatic curve U₂(R)
2. Compute non-adiabatic coupling P₁₂(R) (numerical R-derivative)
3. Solve coupled 2-channel radial equation → E
4. **Validation target:** E ≈ -2.900 Ha (0.12% error)

### Phase 2C: Graph-Theoretic Interpretation
1. Cast the coupled-channel equation as a graph Hamiltonian
2. Identify the natural kinetic scale (if it exists)
3. Compare graph eigenvalue spectrum with S³ and prolate lattices
4. Write up as Paper 13 (or section of Paper 13)

---

## 9. Answer to the Key Question

> *Is the hyperspherical lattice for He structurally the same kind of object
> as the atomic S³ lattice and the molecular prolate spheroidal lattice?*

**YES, with one important generalization.**

All three lattices share the same fundamental structure:
- **Nodes** carry quantum number labels
- **Edges** encode transition amplitudes
- **The Hamiltonian is a sparse matrix** on this graph

The generalization is that the hyperspherical lattice is a **coupled channel
graph** — a fiber bundle where the angular (channel) structure varies with the
radial coordinate R. The S³ lattice is a trivial bundle (no radial coordinate).
The prolate lattice is also trivial (exact separation → no inter-channel
coupling). The hyperspherical lattice is the first non-trivial bundle in GeoVac.

**Structurally, the difference is:**

| Feature | S³ | Prolate | Hyperspherical |
|:---|:---:|:---:|:---:|
| Angular solve | Algebraic (SO(4)) | Spectral (1D ODE) | Spectral (coupled 1D) |
| Radial solve | None (eigenvalue) | FD (1D ODE) | FD (coupled 1D) |
| Inter-channel coupling | None | None | P_μν(R) ≠ 0 |
| Separation | Exact (Fock) | Exact (ξ,η,φ) | Approximate (adiabatic) |
| Graph type | Fixed | Fixed | R-dependent |

The natural geometry principle **does extend** to multi-particle systems. The
geometry is more complex (a bundle, not a product), but the computational
strategy (spectral angular + FD radial) is the same. The path to Paper 13
is clear.

**What is structurally different:**
The S³ lattice exploits an exact symmetry (SO(4)) to avoid any numerical solve.
The prolate lattice exploits exact separation to decouple angular and radial.
The hyperspherical lattice has neither — it requires numerical solution of a
coupled system. But the coupling is weak (2 channels suffice) and the matrices
are tiny (~600 × 600), so this is computationally trivial.

The real insight is that **the cusp — the feature that limits all other
approaches — becomes a boundary condition on the angular problem**, not a
singularity in any coordinate. This is the natural geometry principle at work:
choose coordinates where the hard physics becomes tractable.
