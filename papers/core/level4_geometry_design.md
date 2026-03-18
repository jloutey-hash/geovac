# Level 4 Natural Geometry: Two-Center, Two-Electron Systems

**Status:** Design Document (Mathematical Specification)
**Date:** March 16, 2026
**Target System:** H₂ (two-center, two-electron)
**Prerequisite Papers:** 7 (S³/Fock), 11 (prolate spheroidal), 12 (Neumann V_ee), 13 (hyperspherical)

---

## 1. The Level 4 Problem

The natural geometry hierarchy resolves singularities by coordinate choice:

| Level | System | Geometry | Singularity Resolved | Paper |
|:---:|:---|:---|:---|:---:|
| 1 | H (1-center, 1e⁻) | S³ (Fock) | 1/r | 7 |
| 2 | H₂⁺ (2-center, 1e⁻) | Prolate spheroid | 1/r_A + 1/r_B | 11 |
| 3 | He (1-center, 2e⁻) | Hyperspherical | 1/r₁₂ (e-e cusp) | 13 |
| **4** | **H₂ (2-center, 2e⁻)** | **?** | **All five cusps** | **This document** |

H₂ has **five Coulomb singularities**: 1/r₁A, 1/r₁B, 1/r₂A, 1/r₂B, and 1/r₁₂.
Paper 12 showed that prolate spheroidal CI (which resolves the four nuclear cusps)
saturates at 92.4% of the dissociation energy D_e. The remaining 7.6% gap is the
electron-electron cusp, which is a coordinate singularity in (ξ₁, η₁, ξ₂, η₂) space.

The hypothesis: Level 4 requires a **fiber bundle** that combines prolate spheroidal
structure (nuclear cusps) with hyperspherical structure (e-e cusp).

---

## 2. Coordinate System

### 2.1 Physical Setup

Two nuclei A, B with charges Z_A, Z_B separated by R along the z-axis, positioned
at z = -R/2 (A) and z = +R/2 (B). Two electrons at positions **r₁**, **r₂** measured
from the molecular midpoint.

Full configuration space: 7D = (R, r₁, r₂) with R the internuclear distance and
each rᵢ ∈ ℝ³. After separating the nuclear center-of-mass and exploiting axial
symmetry (Σ states, total M = 0), the effective space is:

- **R** — internuclear distance (1D)
- **Electronic coordinates** — 5D effective (from 6D after M = 0 projection)

### 2.2 Molecule-Frame Hyperspherical Coordinates

Define electron positions in the body-fixed frame (z-axis along internuclear axis):

**Hyperradial coordinates** (from r₁, r₂ measured from molecular midpoint):

    R_e = √(r₁² + r₂²)          electronic hyperradius
    α   = arctan(r₂/r₁)         correlation angle, α ∈ [0, π/2]

so that r₁ = R_e cos α, r₂ = R_e sin α.

**Angular coordinates** (direction of each electron in the molecular frame):

    θ₁  = polar angle of r̂₁ from z-axis     θ₁ ∈ [0, π]
    θ₂  = polar angle of r̂₂ from z-axis     θ₂ ∈ [0, π]
    Φ   = φ₁ - φ₂  relative azimuthal angle  Φ ∈ [0, 2π)

For Σ states (M = m₁ + m₂ = 0), the wavefunction depends on Φ = φ₁ - φ₂, not
on φ₁ and φ₂ independently. This eliminates one angle.

**Full coordinate set** (7D):

    (R, R_e, α, θ₁, θ₂, Φ, φ_cm)

where φ_cm = φ₁ + φ₂ is cyclic for Σ states and integrates out trivially.

**Effective coordinates** (6D for Σ states at fixed R):

    (R_e, α, θ₁, θ₂, Φ)

### 2.3 Relation to Prolate Spheroidal Coordinates

For each electron, the prolate spheroidal coordinates (ξᵢ, ηᵢ) are:

    ξᵢ = (rᵢA + rᵢB) / R       ξᵢ ∈ [1, ∞)
    ηᵢ = (rᵢA - rᵢB) / R       ηᵢ ∈ [-1, 1]

where rᵢA, rᵢB are distances from electron i to nuclei A, B. These relate to
the body-frame spherical coordinates via:

    rᵢA = √(rᵢ² + R²/4 + rᵢR cos θᵢ)
    rᵢB = √(rᵢ² + R²/4 - rᵢR cos θᵢ)

so ξᵢ and ηᵢ are functions of (rᵢ, θᵢ) = (R_e · trig(α), θᵢ):

    ξᵢ = ξᵢ(R_e, α, θᵢ; R)
    ηᵢ = ηᵢ(R_e, α, θᵢ; R)

The key point: ξᵢ, ηᵢ are **not** independent of the hyperspherical coordinates.
They are derived quantities. The hyperspherical coordinates (R_e, α) are more
fundamental because they directly encode the e-e cusp structure.

### 2.4 Volume Element (Jacobian)

The volume element in molecule-frame hyperspherical coordinates:

    dV = R_e⁵ cos²α sin²α · sin θ₁ sin θ₂ · dR_e dα dθ₁ dθ₂ dΦ

This follows from the standard 6D two-particle volume element:

    dV₆ = r₁² r₂² sin θ₁ sin θ₂ dr₁ dr₂ dθ₁ dθ₂ dφ₁ dφ₂

with the substitution r₁ = R_e cos α, r₂ = R_e sin α giving the Jacobian
∂(r₁,r₂)/∂(R_e,α) = R_e, combined with r₁²r₂² = R_e⁴ cos²α sin²α.

After extracting the wavefunction as:

    Ψ = R_e^(-5/2) (cos α sin α)^(-1) F(R_e) Φ(α, θ₁, θ₂, Φ)

the Jacobian centrifugal potential is 15/(8R_e²), identical to the He case (Paper 13).
This is a consequence of the 6D electronic space — the hyperspherical decomposition
is independent of the external potential.

### 2.5 The Controlling Parameter ρ

The ratio of nuclear separation to electronic hyperradius defines the regime:

    ρ ≡ R / (2R_e)

| ρ | Physical Regime | Effective Geometry |
|:---:|:---|:---|
| ρ → 0 | United atom (R_e ≫ R) | Level 3 (hyperspherical He-like) |
| ρ ~ 1 | Molecular (R_e ~ R) | **Level 4** (full coupling) |
| ρ → ∞ | Separated atoms (R_e ≪ R) | Two independent Level 1 atoms |

This dimensionless ratio is the **natural adiabatic parameter** for the Level 4
geometry. The angular eigenvalue problem (Section 3) depends on ρ parametrically,
just as the hyperspherical charge function depends on Z.

---

## 3. Hamiltonian and Separation Structure

### 3.1 Full Electronic Hamiltonian (at fixed R)

In atomic units, with R fixed:

    H_elec = -½∇₁² - ½∇₂² + V(r₁, r₂; R)

where the potential is:

    V = V_nuc + V_ee + V_NN

    V_nuc = -Z_A/r₁A - Z_B/r₁B - Z_A/r₂A - Z_B/r₂B     (nuclear attraction)
    V_ee  = 1/r₁₂                                          (electron repulsion)
    V_NN  = Z_A Z_B / R                                     (nuclear repulsion, constant)

### 3.2 Kinetic Energy in Hyperspherical Coordinates

The 6D electronic Laplacian decomposes as (standard result, cf. Lin 1995 Eq. 2.6):

    ∇₁² + ∇₂² = ∂²/∂R_e² + (5/R_e) ∂/∂R_e + (1/R_e²) Λ²

where Λ² is the **grand angular momentum** on S⁵. After the substitution
Ψ = R_e^(-5/2) F(R_e) Φ(Ω):

    -½[d²F/dR_e² + 2μ(R_e)/R_e² F + 15/(8R_e²) F] = E F

where μ(R_e) is the eigenvalue of the angular problem. This is structurally
identical to the He case (Paper 13, Eq. 12).

The grand angular operator Λ² has the same form as in He:

    Λ² = -(1/sin²α cos²α) ∂_α[sin²α cos²α ∂_α·]
         + l̂₁²/cos²α + l̂₂²/sin²α

where l̂ᵢ² are the angular momentum operators for electrons 1 and 2. In He,
these act on spherical harmonics. In H₂, the **molecular frame** breaks spherical
symmetry, and l̂ᵢ² must be replaced by their body-frame representations.

### 3.3 The Molecular Charge Function

**This is the central object of Level 4.**

In He, the charge function is (Paper 13, Eq. 8):

    C_He(α, θ₁₂) = -Z(1/cos α + 1/sin α) + 1/√(1 - sin 2α cos θ₁₂)

In H₂, the potential does not factor into a single 1/R_e term. Instead:

    V(R_e, Ω; R) = (1/R_e) C_ee(α, θ₁₂) + V_nuc(R_e, α, θ₁, θ₂; R)

where the e-e part retains the He-like form:

    C_ee(α, θ₁₂) = 1/√(1 - sin 2α cos θ₁₂)

but the nuclear attraction does **not** scale as 1/R_e:

    -Z/rᵢA = -Z / √(rᵢ² + R²/4 + rᵢR cos θᵢ)

For electron 1 (r₁ = R_e cos α):

    -Z/r₁A = -Z / √(R_e² cos²α + R²/4 + R_e R cos α cos θ₁)
            = -Z / (R_e cos α) · 1/√(1 + ρ²/cos²α + 2ρ cos θ₁/cos α)

Using ρ = R/(2R_e), this becomes:

    -Z/r₁A = -(Z/R_e) · g_A(α, θ₁; ρ)

where:

    g_A(α, θ₁; ρ) = 1 / (cos α · √(1 + ρ²/cos²α + 2ρ cos θ₁/cos α))

Similarly for nucleus B (cos θ → -cos θ in g_B).

**The molecular charge function** is therefore:

    C_mol(α, θ₁, θ₂, Φ; ρ) = C_ee(α, θ₁₂)
                               - Z[g_A(α, θ₁; ρ) + g_B(α, θ₁; ρ)]
                               - Z[g_A(α, θ₂; ρ) + g_B(α, θ₂; ρ)]

where θ₁₂ is the interelectron angle:

    cos θ₁₂ = cos θ₁ cos θ₂ + sin θ₁ sin θ₂ cos Φ

**Critical difference from He:** C_mol depends on ρ = R/(2R_e), making the angular
eigenvalue problem **R_e-dependent** through the nuclear terms. The e-e cusp term
C_ee is ρ-independent (purely angular), but the nuclear terms depend on ρ.

### 3.4 Coupling Structure Summary

    ┌──────────────────────────────────────────────────────────┐
    │  Term             │ Couples which coordinates?           │
    │──────────────────────────────────────────────────────────│
    │  T_radial         │ R_e only                             │
    │  Λ²/R_e²          │ α, (θ₁,θ₂,Φ) — angular only         │
    │  V_ee = C_ee/R_e  │ α, θ₁₂(θ₁,θ₂,Φ), R_e (via 1/R_e)  │
    │  V_nuc            │ R_e, α, θ₁, θ₂ (all but Φ for σ)   │
    │  V_NN             │ None (constant at fixed R)           │
    │  15/(8R_e²)       │ R_e only (Jacobian centrifugal)      │
    └──────────────────────────────────────────────────────────┘

The nuclear potential V_nuc is the **sole source of non-separability** between
the hyperradial coordinate R_e and the angular coordinates. In He (ρ = 0), V_nuc
becomes -Z(1/cos α + 1/sin α)/R_e and the system is separable in the adiabatic
sense. For H₂, the ρ-dependence means the angular eigenstates Φ_ν change shape
as R_e varies — exactly the content of the non-adiabatic coupling.

### 3.5 Location of the Electron-Electron Cusp

The e-e cusp 1/r₁₂ diverges when r₁₂ → 0. In these coordinates:

    r₁₂² = r₁² + r₂² - 2r₁r₂ cos θ₁₂
          = R_e²(1 - sin 2α cos θ₁₂)

so r₁₂ = 0 requires:

    sin 2α cos θ₁₂ = 1  →  α = π/4  AND  θ₁₂ = 0

This is identical to the He case (Paper 13). The cusp locus is a codimension-2
surface in the 5D angular space, independent of ρ and R_e. The Kato cusp
condition (∂Ψ/∂r₁₂|_{r₁₂=0} = ½Ψ|_{r₁₂=0}) imposes a boundary condition on
the angular eigenfunction at (α = π/4, θ₁₂ = 0).

**This is the key advantage of the hyperspherical decomposition:** the cusp
is always a boundary condition, never a coordinate singularity, regardless of
the molecular geometry.

---

## 4. Adiabatic Approximation

### 4.1 Two-Stage Born-Oppenheimer (Double Adiabatic)

The Level 4 problem has a natural **double adiabatic** structure:

```
    ┌─────────────────────────────────────────────────────────┐
    │  STAGE 1 (OUTER): Nuclear Born-Oppenheimer              │
    │  ─────────────────────────────────────────               │
    │  Slow variable: R (internuclear distance)                │
    │  Fast variables: all electronic coordinates              │
    │  Eigenvalue: E_elec(R) → molecular PES                   │
    │                                                          │
    │  ┌───────────────────────────────────────────────────┐   │
    │  │  STAGE 2 (INNER): Electronic Hyperspherical       │   │
    │  │  ─────────────────────────────────────────        │   │
    │  │  Slow variable: R_e (electronic hyperradius)      │   │
    │  │  Fast variables: (α, θ₁, θ₂, Φ) — angular        │   │
    │  │  Eigenvalue: μ_ν(R_e; R) → adiabatic channels     │   │
    │  │  Parameter: ρ = R/(2R_e)                           │   │
    │  └───────────────────────────────────────────────────┘   │
    │                                                          │
    │  1D radial equation in R_e with potential V_ν(R_e; R)    │
    │  → electronic energy E_elec(R)                           │
    │                                                          │
    │  Then: 1D nuclear equation in R with potential E_elec(R) │
    │  → rovibrational levels                                  │
    └─────────────────────────────────────────────────────────┘
```

### 4.2 The Angular Eigenvalue Problem (Fixed R_e, R)

At fixed R and R_e, solve for eigenstates of the angular Hamiltonian:

    [½Λ² + R_e · C_mol(α, θ₁, θ₂, Φ; ρ)] Φ_ν(Ω; R_e, R) = μ_ν(R_e; R) Φ_ν

This is the **generalization of Paper 13 Eq. 9** to two centers.

**Partial-wave expansion** (L = 0, Σ_g sector):

For the single-center He case, the angular eigenfunctions expand as:

    Φ(α, θ₁₂) = Σ_l φ_l(α) · [(-1)^l/√(2l+1)] P_l(cos θ₁₂)

where the selection rule l₁ = l₂ = l is forced by L = 0 and parity.

For H₂, the molecular field breaks spherical symmetry. Each electron's angular
part must be expanded in **Legendre polynomials of cos θᵢ** (the angle from the
molecular axis), not spherical harmonics of the full solid angle:

    Φ(α, θ₁, θ₂, Φ) = Σ_{l₁ l₂ l} φ_{l₁ l₂ l}(α) · Y_{l₁}^0(θ₁) Y_{l₂}^0(θ₂) ·
                        [coupling coefficients for θ₁₂ via l]

For σ-only configurations (m₁ = m₂ = 0), this simplifies to:

    Φ(α, θ₁, θ₂) = Σ_{l₁ l₂} φ_{l₁ l₂}(α) · P_{l₁}(cos θ₁) P_{l₂}(cos θ₂)

with the Σ_g symmetry constraint φ_{l₁ l₂}(α) = φ_{l₂ l₁}(π/2 - α).

**Coupling through V_ee:** The 1/r₁₂ term couples different (l₁, l₂) channels
through the multipole expansion:

    1/r₁₂ = Σ_l (4π/(2l+1)) (r_<^l / r_>^(l+1)) Σ_m Y_l^m*(θ₁,φ₁) Y_l^m(θ₂,φ₂)

which in α coordinates (with r₁ = R_e cos α, r₂ = R_e sin α) gives:

    C_ee = (1/R_e) Σ_l (4π/(2l+1)) f_l(α) Σ_m Y_l^m*(θ₁,φ₁) Y_l^m(θ₂,φ₂)

where f_l(α) = (tan α)^l for α < π/4, (cot α)^l for α > π/4. The Gaunt
integral structure is identical to the He case (Paper 13, Sec. III.C).

**Coupling through V_nuc:** The nuclear attraction couples l₁ channels for
electron 1 and l₂ channels for electron 2 independently. In the Legendre basis:

    ⟨P_{l₁'}|g_A(α, θ₁; ρ)|P_{l₁}⟩ = G_{l₁' l₁}(α; ρ)

These integrals involve Legendre functions of the nuclear attraction kernel
and can be computed by expanding g_A in Legendre polynomials of cos θᵢ:

    1/√(1 + a² + 2a cos θ) = Σ_l (-a)^l P_l(cos θ)    for |a| < 1

where a = ρ/cos α. This converges for ρ < cos α, i.e., when the electron is
farther from the midpoint than the nucleus. The convergence fails when an
electron is near a nucleus — precisely where the nuclear cusp lives.

**Open Question 1:** The Legendre expansion of the nuclear attraction converges
conditionally when ρ/cos α > 1 (electron between the nuclei, close to a nucleus).
A re-expansion around each nucleus (cf. Morse & Feshbach, §10.3) or a numerical
quadrature strategy is needed for this regime. See Section 7.

### 4.3 Liouville Substitution (α equation)

Following Paper 13, substitute u_{l₁ l₂}(α) = sin α cos α · φ_{l₁ l₂}(α) to
convert the weighted Sturm-Liouville problem to standard Schrödinger form:

    -½ u'' + V_eff(α; ρ) u = μ u

with Dirichlet boundary conditions u(0) = u(π/2) = 0.

The effective potential has contributions:

    V_eff(α) = V_centrifugal(α) + V_Liouville(α) + R_e · V_nuc(α) + R_e · V_ee(α)

where:
- V_centrifugal = [l₁(l₁+1)/cos²α + l₂(l₂+1)/sin²α] / 2  (angular momentum barriers)
- V_Liouville = -2  (curvature of S⁵ metric, cf. Paper 13 Eq. 14)
- V_nuc = nuclear contribution (ρ-dependent, couples channels)
- V_ee = e-e Gaunt coupling (ρ-independent in α, couples channels via l)

### 4.4 Effective Radial Equation in R_e

Given the angular eigenvalues μ_ν(R_e; R), the hyperradial equation is:

    -½ d²F_ν/dR_e² + [μ_ν(R_e; R)/R_e² + 15/(8R_e²)] F_ν = E_elec(R) F_ν

with F_ν(0) = F_ν(∞) = 0.

**This is structurally identical to the He radial equation** (Paper 13, Eq. 18),
but with μ_ν depending on R as a parameter. The non-adiabatic corrections
(off-diagonal Born-Oppenheimer coupling between channels ν, ν') have the same
form as Paper 13, Section VII.

### 4.5 Outer Equation: Nuclear Motion

The electronic energy E_elec(R) from the inner problem serves as the potential
for nuclear motion:

    -½μ_nuc d²χ/dR² + [E_elec(R) + Z_A Z_B/R + J(J+1)/(2μ_nuc R²)] χ = E_total χ

This is standard Born-Oppenheimer. The nuclear lattice approach (Paper 10)
can be applied directly.

---

## 5. Limiting Cases

### 5.1 Reduction to Level 3: United Atom (R → 0)

When R → 0 (ρ → 0 for all R_e):

    g_A(α, θ; 0) = 1/cos α     (nucleus A at midpoint)
    g_B(α, θ; 0) = 1/cos α     (nucleus B at midpoint, same position)

The nuclear potential becomes:

    V_nuc → -(Z_A + Z_B)(1/cos α + 1/sin α) / R_e

which is exactly the He-like charge function with Z = Z_A + Z_B. All θ-dependence
drops out. The angular problem reduces to the 1D α-equation of Paper 13. ✓

### 5.2 Reduction to Level 2: One-Electron Limit

Remove electron 2 (formally: set α = 0, which sends r₂ → 0, r₁ → R_e).
The electronic problem becomes:

    H = -½∇₁² - Z_A/r₁A - Z_B/r₁B

In body-frame coordinates (r₁ = R_e, θ₁):

    r₁A = √(R_e² + R²/4 + R_e R cos θ₁)
    r₁B = √(R_e² + R²/4 - R_e R cos θ₁)

The separation in prolate spheroidal coordinates (ξ₁, η₁) follows. The
hyperradial coordinate R_e maps to a function of ξ₁: specifically,
R_e² = (R/2)²(ξ₁² + η₁² - 1), and the angular variable θ₁ maps to η₁.

The adiabatic separation (fixing R_e, solving for θ₁) becomes the η-equation
of Paper 11, while the R_e equation becomes the ξ-equation. ✓

### 5.3 Reduction to Level 1: United Atom, One Electron

Taking both R → 0 and removing one electron: standard hydrogen on S³. ✓

### 5.4 Separated Atom Limit (R → ∞)

When R → ∞ (ρ → ∞):
- Each electron localizes on one nucleus: α → 0 (electron 1 on A, electron 2
  far away) or α → π/2 (reversed)
- The angular eigenvalue problem develops two isolated wells in α near 0 and π/2
- Each well reproduces a hydrogen-like atom with Z = Z_A or Z_B
- The dissociation limit is E → E_atom(Z_A) + E_atom(Z_B)

---

## 6. Fiber Bundle Structure

### 6.1 Bundle Diagram

```
    Total space E (7D)
    ─────────────────
    │
    ├── Base B₁: R_nuc ∈ (0, ∞)        ← nuclear geometry (Stage 1)
    │   │
    │   └── Fiber F₁ over each R:       ← electronic problem (6D)
    │       │
    │       ├── Base B₂: R_e ∈ (0, ∞)   ← electronic hyperradius (Stage 2)
    │       │   │
    │       │   └── Fiber F₂ over each R_e:  ← angular channels (5D)
    │       │       │
    │       │       ├── α ∈ [0, π/2]         ← correlation angle
    │       │       ├── θ₁ ∈ [0, π]          ← mol-frame polar (e⁻ 1)
    │       │       ├── θ₂ ∈ [0, π]          ← mol-frame polar (e⁻ 2)
    │       │       └── Φ ∈ [0, 2π)          ← relative azimuth
    │       │
    │       └── Connection: P_νμ(R_e) = ⟨Φ_ν|∂/∂R_e|Φ_μ⟩
    │           (non-adiabatic coupling, Berry connection)
    │
    └── Connection: P_νμ(R) = ⟨ψ_ν^elec|∂/∂R|ψ_μ^elec⟩
        (nuclear non-adiabatic coupling)
```

### 6.2 Structure Group

The He hyperspherical fiber bundle (Paper 13, Sec. VII) has structure group O(N_ch)
where N_ch is the number of adiabatic channels. The Level 4 bundle has:

- **Inner fiber** (angular channels at fixed R_e): O(N_ang) where N_ang is the
  number of angular channels (determined by l₁_max, l₂_max truncation)
- **Outer fiber** (electronic states at fixed R): O(N_elec) where N_elec is the
  number of hyperradial adiabatic states

The full molecular fiber bundle is a **composition** of two adiabatic bundles:

    E → B₂ → B₁

with the outer connection encoding the derivative of the entire electronic
eigenstate with respect to R.

### 6.3 Comparison with Paper 13

| Feature | Level 3 (He, Paper 13) | Level 4 (H₂) |
|:---|:---|:---|
| Base | R_e (hyperradius) | R × R_e (two slow variables) |
| Fiber | (α, θ₁₂) — 2D effective | (α, θ₁, θ₂, Φ) — 4D effective |
| Charge function | C(α, θ₁₂) — ρ-independent | C_mol(α, θ₁, θ₂, Φ; ρ) — ρ-dependent |
| Angular PDE | 1D in α (after PW expansion) | 2D in (α, coupled l₁l₂ channels) |
| Non-adiabatic | P_νμ(R_e) — 1 parameter | P_νμ(R_e; R) — 2 parameters |
| Cusp location | α = π/4, θ₁₂ = 0 | α = π/4, θ₁₂ = 0 (identical) |

---

## 7. Lattice Construction

### 7.1 Quantum Number Labels (Nodes)

Each node in the Level 4 graph is labeled by:

    |R_i, R_e_j, ν⟩

where:
- R_i: nuclear grid point (i = 1, ..., N_R)
- R_e_j: electronic hyperradial grid point (j = 1, ..., N_Re)
- ν: angular channel index (ν = 1, ..., N_ch)

The angular channels are themselves labeled by (l₁, l₂, n_α) where n_α is the
α-quantum number within the (l₁, l₂) channel.

For the σ-only sector (m₁ = m₂ = 0), channel labels reduce to (l₁, l₂) pairs
with the constraint l₁ + l₂ even (for Σ_g symmetry).

**Channel counting** (σ-only, Σ_g):

| l_max | Channels (l₁,l₂) | Pairs with l₁+l₂ even |
|:---:|:---:|:---:|
| 0 | (0,0) | 1 |
| 1 | (0,0),(1,1) | 2 |
| 2 | (0,0),(1,1),(0,2),(2,0),(2,2) | 5 |
| 3 | (0,0),(1,1),(0,2),(2,0),(2,2),(1,3),(3,1),(3,3) | 8 |

Each channel has N_α grid points in the α discretization.

### 7.2 Selection Rules (Edges)

Edges arise from three sources:

**A. Hyperradial kinetic coupling (tridiagonal in R_e):**

    (R_i, R_e_j, ν) ↔ (R_i, R_e_{j±1}, ν)

Same channel, adjacent R_e grid points. Self-adjoint FD stencil as in Paper 13.

**B. Non-adiabatic coupling (off-diagonal in ν at same R_e):**

    (R_i, R_e_j, ν) ↔ (R_i, R_e_j, μ)    for μ ≠ ν

Mediated by the Berry connection P_νμ(R_e; R). Selection rules from angular
symmetry: ΔL = 0, Δπ = 0 (same total symmetry manifold).

**C. Nuclear kinetic coupling (tridiagonal in R):**

    (R_i, R_e_j, ν) ↔ (R_{i±1}, R_e_j, ν)

Same electronic state, adjacent R grid points.

**D. Nuclear non-adiabatic coupling:**

    (R_i, R_e_j, ν) ↔ (R_i, R_e_j, μ)

Mediated by ⟨ψ_ν^elec|∂/∂R|ψ_μ^elec⟩ at the outer BO level.

**Sparsity:** The full Hamiltonian matrix is block-structured:
- Block-tridiagonal in R_e (kinetic, within each (R_i, ν))
- Block-diagonal couplings in ν (non-adiabatic, within each (R_i, R_e_j))
- Block-tridiagonal in R (nuclear kinetic, within each (R_e_j, ν))

### 7.3 Matrix Dimensions

At the angular level (innermost problem, fixed R_i, R_e_j):

    N_ang = N_ch × N_α

For the inner adiabatic radial problem (fixed R_i, single channel):

    N_inner = N_Re

For the full coupled-channel problem at fixed R_i:

    N_inner_coupled = N_ch × N_Re

For the outer nuclear problem:

    N_outer = N_R

**Total matrix** (full coupled problem, single Σ_g sector):

    N_total = N_R × N_ch × N_Re

With typical values N_R = 50, N_Re = 200, N_ch = 5 (l_max = 2):

    N_total = 50,000

This is large but sparse — each node connects to at most 2(tridiagonal) + N_ch
(non-adiabatic) + 2(nuclear tridiagonal) ≈ 9 neighbors.

### 7.4 Reusable GeoVac Modules

| Module | Reuse | Adaptation Needed |
|:---|:---|:---|
| `hyperspherical_angular.py` | α-discretization, Liouville substitution, Gaunt coupling | Add l₁ ≠ l₂ channels; add nuclear coupling G_{l'l}(α; ρ) |
| `hyperspherical_adiabatic.py` | Adiabatic curve computation (sweep R_e) | Add ρ = R/(2R_e) parameter; compute curves at each R |
| `hyperspherical_radial.py` | Self-adjoint FD radial solver | Reuse directly for R_e equation |
| `hyperspherical_coupling.py` | Non-adiabatic coupling P_νμ | Extend to two-parameter (R_e, R) dependence |
| `prolate_spheroidal_lattice.py` | Angular η-expansion, ξ-radial solver | Reference for nuclear attraction integrals |
| `neumann_vee.py` | Neumann expansion of 1/r₁₂ | Use for validating C_ee against Neumann at large R |
| `molecular_sturmian.py` | Legendre spectral expansion | Reuse for θᵢ angular basis |

---

## 8. Implementation Strategy

### 8.1 Phased Approach

**Phase 1: σ-only, single channel (l₁ = l₂ = 0)**

Solve the simplest version: only s-wave electrons, single adiabatic channel.
The angular problem reduces to 1D in α (identical structure to Paper 13 He
ground state). The molecular charge function at (l₁ = l₂ = 0) is:

    C_mol^{00}(α; ρ) = ⟨Y_0^0(θ₁) Y_0^0(θ₂)|C_mol|Y_0^0(θ₁) Y_0^0(θ₂)⟩

which is an integral over θ₁, θ₂ that can be evaluated by quadrature or by
the Legendre expansion of g_A, g_B.

Expected result: capture the dominant electron-electron cusp physics. Compare
with Paper 12 Neumann V_ee results.

**VALIDATED:** Phase 1 recovers 31% D_e — far below Paper 12 (92.4%). The l=0
projection destroys directional bonding information. See Appendix A.

**Phase 2: Multi-channel σ (l_max = 2-4)**

Add l₁, l₂ > 0 angular channels. The Gaunt coupling structure from Paper 13
carries over. The new element is the nuclear coupling matrix G_{l₁'l₁}(α; ρ).

**IMPORTANT:** Both orderings (l₁, l₂) and (l₂, l₁) must be included as
separate channels because they have different centrifugal potentials.
Gerade symmetry requires l₁+l₂ even; the first nuclear coupling to (0,0)
is quadrupolar, requiring l_max ≥ 2.

**VALIDATED:** Phase 2 exceeds Paper 12 at l_max=4 (95.5% D_e). See Appendix A.

**Phase 3: Include π orbitals (m ≠ 0)**

Allow m₁ = ±1, m₂ = ∓1 with M = 0. This brings in Φ-dependence and doubles
the channel count. Required for excited states (Π, Δ symmetry).

**Phase 4: Full coupled-channel with nuclear dynamics**

Sweep over R, compute E_elec(R) as a molecular PES, solve nuclear equation.
Compare with exact H₂ benchmarks (Kolos & Wolniewicz).

### 8.2 Computational Bottleneck

The innermost angular eigenvalue problem must be solved at **every (R_i, R_e_j)
grid point** — a 2D parameter sweep. With N_R × N_Re = 50 × 200 = 10,000 solves
of an N_ang-dimensional eigenvalue problem, the total cost is:

    Cost ~ N_R × N_Re × N_ang³ (if using dense diagonalization)
         ~ 10,000 × 200³ = 8 × 10¹⁰    (for N_ang = 200)

This is prohibitive. Mitigation strategies:

1. **Interpolation in ρ:** Since the angular problem depends on (R_e, R) only
   through ρ = R/(2R_e), compute angular eigenvalues on a 1D ρ-grid (N_ρ ~ 50)
   and interpolate. This reduces the sweep from 2D to 1D.

2. **Sparse eigensolver:** Only the lowest N_ch eigenvalues are needed. Use
   Lanczos/ARPACK (scipy.sparse.linalg.eigsh) to avoid dense O(N³).

3. **Perturbation theory:** At small ρ, the nuclear coupling is a perturbation
   on the He-like angular problem. At large ρ, the separated-atom limit provides
   a starting point.

---

## 9. Open Questions

### 9.1 Nuclear Attraction Integrals (Critical)

The Legendre expansion of g_A(α, θ; ρ) diverges when ρ/cos α > 1 (electron
between midpoint and nucleus A, closer to A than to midpoint). This corresponds
to ξ < (1 + ρ²/cos²α)^{1/2}, which includes the important bonding region.

**Options:**
- (a) Split the integration domain in α at α* where ρ/cos α* = 1, and use
  different expansions in each region (Morse & Feshbach re-expansion).
- (b) Compute the angular matrix elements ⟨P_{l'}|g_A|P_l⟩ by direct numerical
  quadrature in θ, avoiding the Legendre expansion of g_A altogether.
- (c) Use a prolate spheroidal basis for the θ-dependence (η instead of cos θ)
  which may have better convergence properties near the nuclear cusps.

Recommendation: start with (b) for Phase 1, investigate (c) for optimization.

### 9.2 Cusp Boundary Condition (Important)

In He (Paper 13), the Kato cusp condition at α = π/4 is automatically satisfied
by the finite-difference discretization (the cusp is "seen" through the Gaunt
coupling at the grid boundary). Does this carry over to H₂?

The cusp condition is:

    (∂Ψ/∂r₁₂)_{r₁₂=0} = ½ Ψ|_{r₁₂=0}

In the (α, θ₁₂) coordinates, r₁₂ = 0 is at α = π/4, θ₁₂ = 0. The cusp is a
**local** condition and should be the same in H₂ as in He — the nuclear potential
is finite at the cusp point and does not modify the cusp structure. But this needs
verification.

### 9.3 Convergence of the Partial-Wave Expansion (RESOLVED)

**RESOLVED:** The partial-wave expansion converges much more slowly than in He.
l_max=0 gives only 31% D_e (vs 99.95% for He). The key difference: gerade
symmetry requires l₁+l₂ even, so the first coupling to (0,0) is quadrupolar
(Δl = 2). The jump at l_max=2 (31% → 88%) is dramatic. Convergence:

| l_max | N_ch | D_e % |
|:---:|:---:|:---:|
| 0 | 1 | 31% |
| 1 | 2 | 37% |
| 2 | 5 | 88% |
| 3 | 8 | 89% |
| 4 | 13 | 95.5% |

Grid convergence at l_max=3 shows the l-expansion (not FD resolution) is
the limiting factor. Further improvement requires l_max ≥ 5.

### 9.4 Non-Adiabatic Coupling Computation

Computing P_νμ(R_e; R) = ⟨Φ_ν|∂Φ_μ/∂R_e⟩ requires the derivative of the
angular eigenstates with respect to R_e. Since the angular problem depends on
R_e only through ρ = R/(2R_e), this can be computed as:

    P_νμ(R_e; R) = -R/(2R_e²) · ⟨Φ_ν|∂Φ_μ/∂ρ⟩

The ρ-derivative can be obtained by finite differences on the angular
eigenstates computed at nearby ρ values, or analytically via the
Hellmann-Feynman theorem applied to the ρ-dependent nuclear coupling.

### 9.5 Relationship to James-Coolidge Coordinates

James and Coolidge (1933) used the coordinates (ξ₁, η₁, ξ₂, η₂, r₁₂) for
their variational calculation of H₂. Our coordinate system replaces (ξ₁, ξ₂)
with (R_e, α) and handles r₁₂ through the cusp boundary condition rather
than an explicit coordinate. It would be valuable to establish the explicit
coordinate transformation between the two systems and verify that the
James-Coolidge basis functions can be re-expressed in our coordinates.

### 9.6 Symmetry Exploitation

The H₂ wavefunction has several symmetries:
- **Σ_g:** Inversion through midpoint (gerade)
- **Singlet/triplet:** Electron exchange
- **Nuclear exchange:** For homonuclear (Z_A = Z_B), permutation of nuclei

In the (R_e, α, θ₁, θ₂, Φ) coordinates:
- Electron exchange: (α, θ₁, θ₂) → (π/2 - α, θ₂, θ₁)
- Inversion: (θ₁, θ₂) → (π - θ₁, π - θ₂)

These constrain the angular eigenfunctions and halve the effective Hilbert space.
The implementation should exploit these from the start.

### 9.7 Accuracy Target (PARTIALLY RESOLVED)

Paper 12 achieves 92.4% D_e with prolate spheroidal CI (nuclear cusps resolved,
e-e cusp unresolved). Paper 13 achieves 0.05% error for He (e-e cusp resolved
but one-center).

**Status:** Phase 2 at l_max=4 achieves 95.5% D_e, exceeding Paper 12. The
remaining 4.5% gap likely requires:
- Higher partial waves (l_max = 5-6), ~18-25 channels
- Non-adiabatic coupling between hyperradial channels
- Phase 3 (π-orbital channels, m ≠ 0)

Updated targets:
- Phase 2 (multi-channel σ, l_max=4): **95.5% D_e** ✓ (exceeds Paper 12)
- Phase 3 (add π orbitals): > 98% D_e
- Phase 4 (full coupled-channel): < 0.5% error on D_e

### 9.8 Scaling to Larger Systems

The double-adiabatic structure generalizes:
- More electrons: higher-dimensional hyperspherical decomposition (S^{3N-1}
  for N electrons). The angular problem dimension grows, but the structure
  (base = hyperradius, fiber = angular channels) is preserved.
- More nuclei: the nuclear base space becomes multi-dimensional (all R_ij
  internuclear distances). The inner electronic problem is solved at each
  nuclear geometry.

The computational cost scales as O(N_nuc^geom × N_ρ × N_ang³), which may be
tractable for small molecules (3-4 atoms) but will require further innovation
for larger systems.

---

## 10. References

1. **Fock, V. A.** (1935) Z. Phys. 98, 145. Hydrogen atom in momentum space,
   stereographic projection to S³.

2. **Morse, P. M. & Feshbach, H.** (1953) *Methods of Theoretical Physics*,
   Ch. 10. Two-center coordinates, prolate spheroidal harmonics.

3. **James, H. M. & Coolidge, A. S.** (1933) J. Chem. Phys. 1, 825.
   Variational H₂ calculation in prolate spheroidal + r₁₂ coordinates.

4. **Macek, J.** (1968) J. Phys. B 1, 831. Hyperspherical approach to
   two-electron atoms.

5. **Lin, C. D.** (1995) Phys. Rep. 257, 1. Hyperspherical coordinate approach
   to atomic and other Coulombic three-body systems. [Comprehensive review]

6. **Kolos, W. & Wolniewicz, L.** (1968) J. Chem. Phys. 49, 404. Exact
   H₂ ground state (benchmark for accuracy targets).

7. **GeoVac Papers:**
   - Paper 7: Dimensionless vacuum / S³ conformal geometry
   - Paper 11: Prolate spheroidal lattice (Level 2)
   - Paper 12: Neumann V_ee expansion, cusp diagnosis
   - Paper 13: Hyperspherical lattice (Level 3), fiber bundle structure

---

---

## Appendix A: Validated Numerical Results (March 16, 2026)

### A.1 Phase 1 (l_max=0, single channel)

**Implementation:** `geovac/level4_sigma_channel.py`
**Tests:** `tests/test_level4_sigma.py` (12 tests, all pass)

At R = 1.4 bohr, n_alpha=200, n_Re=400, n_quad=32:
- E_total = -1.0546 Ha (exact: -1.1745 Ha)
- D_e = 0.0546 Ha = **31.3%** of exact
- United-atom limit (rho → 0): matches He l_max=0 eigenvalue to machine precision

**Diagnosis:** The l=0 projection averages over electron directions isotropically,
destroying the directional sigma_g bonding. Nuclear attraction is ~30% weaker
than united-atom at asymmetric alpha values.

### A.2 Phase 2 (multichannel, l_max=0-4)

**Implementation:** `geovac/level4_multichannel.py`
**Tests:** `tests/test_level4_multichannel.py` (12 tests, all pass)

**Convergence at R = 1.4 bohr** (n_alpha=150, n_Re=300, n_quad=24):

| l_max | N_ch | E_total (Ha) | D_e (Ha) | D_e % | Time |
|:---:|:---:|:---:|:---:|:---:|:---:|
| 0 | 1 | -1.0538 | 0.0538 | 30.8% | 2s |
| 1 | 2 | -1.0651 | 0.0651 | 37.3% | 9s |
| 2 | 5 | -1.1532 | 0.1532 | 87.8% | 49s |
| 3 | 8 | -1.1543 | 0.1543 | 88.5% | 118s |
| **4** | **13** | **-1.1667** | **0.1667** | **95.5%** | **333s** |

**Key findings:**
1. l_max=2 is the critical threshold — gerade symmetry (l1+l2 even) means
   the first nuclear coupling to (0,0) is quadrupolar (Delta l = 2)
2. Both orderings (l1,l2) and (l2,l1) must be included as separate channels
3. Grid convergence at l_max=3 plateaus at ~89% — l-expansion is the bottleneck
4. **l_max=4 achieves 95.5%, exceeding Paper 12's 92.4%**

### A.3 Bugs Found and Fixed

1. **Missing channels (critical):** Initial l1<=l2 restriction excluded (2,0),
   (3,1), etc. Nuclear coupling requires ONE l-index to match:
   (0,0) → (0,2) via electron-2 AND (0,0) → (2,0) via electron-1.
   Missing half the coupling channels.

2. **Normalization mismatch (critical):** Nuclear coupling was in raw P_l basis
   but e-e coupling in orthonormal hat_P = sqrt((2l+1)/2) * P_l basis.
   Fixed by multiplying nuclear integrals by sqrt((2l'+1)(2l+1)) factors.

---

**End of Design Document**

*Next step: Phase 3 — add pi-orbital channels (m != 0) to approach 98%+ D_e,
or Phase 4 — sweep R to compute the full PES and rovibrational spectrum.*
