# Helium in Hyperspherical Coordinates: Full Hamiltonian

**Date:** 2026-03-14
**Status:** Phase 1 — Mathematical Setup
**System:** He (Z=2, two electrons, one center)

---

## 1. Coordinate Definition

### 1.1 From Cartesian to Hyperspherical

Two electrons at positions **r₁**, **r₂** ∈ ℝ³ define a 6D configuration space.
Remove the trivial center-of-mass motion (fixed nucleus approximation) to get
6 internal coordinates.

**Hyperspherical coordinates** (Whitten-Smith):

```
R = √(r₁² + r₂²)                 — hyperradius,  R ∈ [0, ∞)
α = arctan(r₂/r₁)                — hyperangle,   α ∈ [0, π/2]
θ₁, φ₁                            — direction of r₁ (spherical angles)
θ₂, φ₂                            — direction of r₂ (spherical angles)
```

**Inverse relations:**

```
r₁ = R cos α
r₂ = R sin α
```

The interelectron angle θ₁₂ is defined by:

```
cos θ₁₂ = r̂₁ · r̂₂ = cos θ₁ cos θ₂ + sin θ₁ sin θ₂ cos(φ₁ - φ₂)
```

---

## 2. Interelectron Distance in Hyperspherical Coordinates

```
r₁₂² = r₁² + r₂² - 2r₁r₂ cos θ₁₂
     = R²cos²α + R²sin²α - 2R²cos α sin α cos θ₁₂
     = R²(1 - sin 2α cos θ₁₂)
```

Therefore:

```
┌─────────────────────────────────────────────┐
│  r₁₂ = R √(1 - sin 2α cos θ₁₂)            │
└─────────────────────────────────────────────┘
```

**Cusp analysis:**
- r₁₂ = 0 requires sin 2α cos θ₁₂ = 1, i.e., **α = π/4 and θ₁₂ = 0**
- This is a *manifold* in the 5D hyperangular space, not a point
- At fixed R, the cusp is a **boundary condition** on the angular problem
- The hyperradius R itself has **no singularity** at the cusp — the
  non-analyticity (R^{1/2}, R ln R) is local in R, not angular

---

## 3. The Full Hamiltonian

### 3.1 Kinetic Energy

The 6D Laplacian in hyperspherical coordinates separates as:

```
∇₁² + ∇₂² = ∂²/∂R² + (5/R) ∂/∂R + (1/R²) Λ²(Ω̂)
```

where **Λ²** is the **grand angular momentum operator** (Laplace-Beltrami on S⁵):

```
Λ² = -1/(sin²α cos²α) ∂/∂α (sin²α cos²α ∂/∂α)
     + l̂₁²/cos²α + l̂₂²/sin²α
```

Here l̂₁² and l̂₂² are the angular momentum operators for electrons 1 and 2.

The kinetic energy operator is:

```
T = -½(∇₁² + ∇₂²)
  = -½[∂²/∂R² + (5/R)∂/∂R] - Λ²/(2R²)
```

It is conventional to extract the R^{-5/2} Jacobian factor. Writing
Ψ = R^{-5/2} F(R) Φ(Ω̂), the hyperradial kinetic energy becomes:

```
T_R = -½ d²F/dR² + (K(K+4) + 15/4)/(2R²) F
```

where K is the grand angular momentum quantum number (eigenvalue of Λ²).

### 3.2 Nuclear Attraction

```
V_nuc = -Z/r₁ - Z/r₂ = -Z/(R cos α) - Z/(R sin α)
      = -(Z/R)(1/cos α + 1/sin α)
```

This is **separable** in R and α: it factors as (1/R) × (angular function).

### 3.3 Electron Repulsion

```
V_ee = 1/r₁₂ = 1/(R√(1 - sin 2α cos θ₁₂))
     = (1/R) · 1/√(1 - sin 2α cos θ₁₂)
```

Also factors as (1/R) × (angular function). The angular part depends on
**both** α and θ₁₂.

### 3.4 The Charge Function

All potential terms scale as 1/R. Define the **charge function** (angular
potential landscape):

```
┌──────────────────────────────────────────────────────────────┐
│  C(α, θ₁₂) = -Z(1/cos α + 1/sin α)                        │
│             + 1/√(1 - sin 2α cos θ₁₂)                      │
│                                                              │
│  V_total = C(α, θ₁₂) / R                                   │
└──────────────────────────────────────────────────────────────┘
```

**Key features of C(α, θ₁₂):**
- Nuclear singularities at α = 0 (electron 1 at nucleus) and α = π/2
  (electron 2 at nucleus) — these are coordinate singularities, removed
  by the centrifugal term from l̂₁², l̂₂²
- Electron-electron coalescence at α = π/4, θ₁₂ = 0 — the **cusp**
- Wannier saddle at α = π/4, θ₁₂ = π: C_saddle = -2√2 Z + 1/√2

For Z = 2: C_saddle ≈ -4.950

---

## 4. Separation Structure

### 4.1 What Separates

The Hamiltonian at fixed R is:

```
H_ang(R) = Λ²/(2R²) + C(α, θ₁₂)/R
         = (1/R)[Λ²/(2R) + C(α, θ₁₂)]
```

At each fixed R, this is a **self-adjoint eigenvalue problem on S⁵**:

```
[Λ²/2 + R·C(α, θ₁₂)] Φ_μ(R; Ω̂) = U_μ(R) Φ_μ(R; Ω̂)
```

The eigenvalues U_μ(R) are the **adiabatic potential curves**.

### 4.2 What Couples

**Within the angular problem** (at fixed R):
- Λ² couples different (l₁, l₂) via the α-dependent terms
- V_ee couples different l₁, l₂ channels via the θ₁₂ dependence
- Nuclear attraction couples different α states

**Between radial and angular:**
- Non-adiabatic couplings: ⟨Φ_μ|∂/∂R|Φ_ν⟩ and ⟨Φ_μ|∂²/∂R²|Φ_ν⟩
- These are **off-diagonal** in channel index ν and drive transitions
  between adiabatic curves

### 4.3 L=0 Partial-Wave Reduction

For He ground state (L=0, M=0), expand in coupled angular momenta:

```
Φ(α, r̂₁, r̂₂) = Σ_l φ_l(α) · [Y_l(r̂₁) ⊗ Y_l(r̂₂)]^{L=0}
               = Σ_l φ_l(α) · (-1)^l/√(2l+1) · P_l(cos θ₁₂)
```

This reduces the 5D angular problem to a set of **coupled 1D equations in α**:

```
[-1/(sin²α cos²α) d/dα(sin²α cos²α d/dα) + l(l+1)/cos²α + l(l+1)/sin²α] φ_l(α)
+ R · Σ_{l'} V_{ll'}(α) φ_{l'}(α) = 2R · U_μ(R) · φ_l(α)
```

where the coupling matrix is:

```
V_{ll'}(α) = ⟨P_l(cos θ₁₂)| C(α, θ₁₂) |P_{l'}(cos θ₁₂)⟩_{θ₁₂}
```

**Nuclear attraction** is diagonal in l (no θ₁₂ dependence):

```
V_{ll'}^{nuc}(α) = -Z(1/cos α + 1/sin α) · δ_{ll'}
```

**Electron repulsion** couples different l channels:

```
⟨P_l| 1/√(1 - sin 2α cos θ) |P_{l'}⟩_θ = Σ_k c_{kll'} (sin 2α)^k
```

The multipole expansion of 1/r₁₂ gives:

```
1/r₁₂ = (1/R) Σ_k (r_</r_>)^k P_k(cos θ₁₂)
       = (1/R) Σ_k f_k(α) P_k(cos θ₁₂)
```

where:

```
f_k(α) = {  (tan α)^k  if α ≤ π/4  (r₂ ≤ r₁)
          {  (cot α)^k  if α > π/4  (r₂ > r₁)
```

The coupling integral becomes:

```
V_{ll'}^{ee}(α) = (2/R) Σ_k f_k(α) ⟨P_l|P_k|P_{l'}⟩

⟨P_l|P_k|P_{l'}⟩ = ∫₋₁¹ P_l(x) P_k(x) P_{l'}(x) dx
                   = 2 · (l  k  l')²    ← Wigner 3j symbol squared
                         (0  0  0 )
```

This is nonzero only when l+k+l' is even and |l-l'| ≤ k ≤ l+l'.

---

## 5. The Hyperradial Equation

After solving the angular problem at each R to get U_μ(R), the
hyperradial equation for channel μ is:

```
┌───────────────────────────────────────────────────────────────┐
│  [-½ d²/dR² + U_μ(R)/R + (μ_eff² - 1/4)/(2R²)] F_μ(R)     │
│  - Σ_{ν≠μ} [P_μν(R) dF_ν/dR + ½ Q_μν(R) F_ν(R)] = E F_μ(R)│
└───────────────────────────────────────────────────────────────┘
```

where:
- μ_eff = K + 2 (effective angular momentum, K = grand angular momentum)
- P_μν(R) = ⟨Φ_μ|∂/∂R|Φ_ν⟩ — first-derivative coupling
- Q_μν(R) = ⟨Φ_μ|∂²/∂R²|Φ_ν⟩ — second-derivative coupling

### 5.1 Asymptotic Behavior

**R → 0:**
```
U_μ(R) → (K(K+4))/2 + R·C_min + O(R²)

where C_min is the minimum of C(α, θ₁₂) ≈ -4.950 for Z=2
```

The centrifugal barrier (K+2)²/(2R²) dominates — wavefunction vanishes.

**R → ∞:**
```
U_μ(R)/R → -Z²/2 + (ionization threshold)
```

The lowest curve tends to the He⁺(1s) + e⁻ threshold at E = -Z²/2 = -2.0 Ha.

### 5.2 Potential Curve Structure

For Z=2 He, the lowest adiabatic curves are:

```
μ=1:  Deep well, minimum at R ≈ 1.5 bohr, supports bound states
      U₁(R)/R → -2.0 Ha as R → ∞ (He⁺ 1s threshold)

μ=2:  Shallower well, avoided crossing with μ=1 near R ≈ 3 bohr
      U₂(R)/R → -0.5 Ha as R → ∞ (He⁺ 2s threshold)

μ≥3:  Weakly bound or unbound Rydberg channels
```

---

## 6. Volume Element and Inner Products

### 6.1 Jacobian

The volume element in hyperspherical coordinates is:

```
dV₆ = R⁵ dR · sin²α cos²α dα · sin θ₁ dθ₁ dφ₁ · sin θ₂ dθ₂ dφ₂
     = R⁵ dR · dΩ₅
```

where dΩ₅ is the volume element on S⁵.

### 6.2 Inner Products

Hyperradial:
```
⟨F_μ|F_ν⟩_R = ∫₀^∞ F_μ(R) F_ν(R) dR    (after extracting R^{-5/2})
```

Hyperangular (L=0 sector):
```
⟨φ_l|φ_{l'}⟩_α = ∫₀^{π/2} φ_l(α) φ_{l'}(α) sin²α cos²α dα
```

The weight function sin²α cos²α ensures proper normalization on S⁵.

---

## 7. Summary of Separation Structure

```
┌─────────────────────────────────────────────────────────┐
│              FULL 6D PROBLEM                            │
│  H = T_R + Λ²/(2R²) + C(α,θ₁₂)/R                     │
│                                                         │
│  Step 1: Fix L=0 → reduce to 3D (R, α, θ₁₂)           │
│                                                         │
│  Step 2: Expand θ₁₂ in Legendre → coupled 1D in α      │
│          (finite set of coupled ODEs, l = 0,...,l_max)   │
│                                                         │
│  Step 3: Solve α-equations at each R → U_μ(R) curves    │
│          (matrix eigenvalue problem, same as prolate η)  │
│                                                         │
│  Step 4: Solve hyperradial equation with U_μ(R)/R       │
│          (1D Schrödinger with effective potential,       │
│           same as prolate ξ equation)                    │
│                                                         │
│  Step 5: Non-adiabatic coupling P_μν, Q_μν for          │
│          multi-channel accuracy                          │
└─────────────────────────────────────────────────────────┘
```

**Structural isomorphism with prolate spheroidal lattice:**

| Prolate (H₂⁺) | Hyperspherical (He) |
|:---|:---|
| ξ (radial) | R (hyperradial) |
| η (angular) | α (hyperangular) |
| φ (azimuthal, separates) | θ₁₂ (interelectron, partial-wave expansion) |
| Separation constant A | Adiabatic eigenvalue U_μ(R) |
| Legendre spectral in η | Legendre spectral in α + coupled channels |
| FD grid in ξ | FD grid in R |
| Brent root-finding on c² | Brent root-finding on E (or direct eigenvalue) |
| 1 uncoupled channel | N_ch coupled channels (N_ch ≈ 2 for 1% accuracy) |
