# Literature Survey: Hyperspherical Methods for Two-Electron Systems

**Date:** 2026-03-14
**Status:** Phase 1 — Mathematical Reconnaissance
**Focus:** Quantum numbers, channel structure, and mapping to GeoVac (n,l,m)

---

## 1. Fock (1954): Exact Two-Electron Expansion Near R = 0

**Reference:** V.A. Fock, *Izv. Akad. Nauk SSSR, Ser. Fiz.* **18**, 161 (1954)

### Key Results

The exact two-electron wavefunction near the triple coalescence point
(both electrons at the nucleus) has the expansion:

```
Ψ(R, Ω̂) = Σ_{n=0}^∞ Σ_{k=0}^{⌊n/2⌋} R^n (ln R)^k f_{nk}(Ω̂)
```

where **R = √(r₁² + r₂²)** is the hyperradius and **Ω̂** collectively denotes
the five hyperangular coordinates.

**Critical non-analytic terms:**
- `R^(1/2)` terms — half-integer powers of the hyperradius
- `R ln R` terms — logarithmic corrections

These arise from the electron-electron cusp (Kato condition) and **cannot be
represented by any finite polynomial basis** in single-electron coordinates.
This is the fundamental reason Paper 12's prolate spheroidal CI saturates at
92.4% D_e.

### Implication for GeoVac
The hyperradius R is the *natural* radial variable for two-electron systems.
A lattice discretized in R (not r₁, r₂ separately) can represent the cusp
terms as local power-law behavior at small R, analogous to how the S³ lattice
naturally represents 1/r via Fock's 1935 stereographic projection.

---

## 2. Macek (1968): Adiabatic Hyperspherical Method

**Reference:** J. Macek, *J. Phys. B* **1**, 831 (1968)

### The Adiabatic Approximation

Macek introduced the **adiabatic hyperspherical method**: treat R as a slow
(adiabatic) variable and solve for the hyperangular motion at each fixed R.

**Procedure:**
1. At fixed R, solve the **hyperangular eigenvalue problem**:

   ```
   [Λ²(Ω̂) + V(R, Ω̂)] Φ_μ(R; Ω̂) = U_μ(R) Φ_μ(R; Ω̂)
   ```

   where Λ² is the grand angular momentum operator on S⁵ and V contains
   nuclear attraction + electron repulsion (both scale as 1/R at fixed angles).

2. The eigenvalues U_μ(R) are **adiabatic potential curves** — effective
   potentials for the hyperradial motion.

3. Solve the **hyperradial equation**:

   ```
   [-d²/dR² + (μ(μ+4) + 15/4)/R² + U_μ(R)/R] F_μ(R) = E F_μ(R)
   ```

   (after extracting the centrifugal term from U_μ).

### Channel Functions and Quantum Numbers

Each channel μ is labeled by:
- **L** — total orbital angular momentum (conserved, good quantum number)
- **M** — z-projection of L (conserved)
- **S** — total spin (0 = singlet, 1 = triplet)
- **π** — parity under inversion: (-1)^{l₁+l₂} (conserved)
- **ν** — channel index within a given (L, S, π) symmetry block

For He ground state (¹S): L=0, S=0, π=even. The channels are labeled by ν alone.

### Key Finding: Channel Convergence
Macek showed that for He ground state, **a single adiabatic channel** (the
lowest ¹S curve) gives E ≈ -2.879 Ha (0.85% error). Adding the first
excited channel via coupling gives E ≈ -2.900 Ha (0.12% error). The coupled
2-channel result nearly matches the exact -2.9037 Ha.

### Implication for GeoVac
This is very promising: **2 coupled channels suffice for sub-1% accuracy**.
The adiabatic potential curves U_μ(R) play the same structural role as the
separation constant A in the prolate spheroidal lattice — they couple the
angular and radial problems.

---

## 3. Lin (1995): Comprehensive Review of Hyperspherical Methods

**Reference:** C.D. Lin, *Phys. Rep.* **257**, 1 (1995)

### Coordinate Systems

Lin systematizes three common hyperspherical coordinate choices:

**(a) Whitten-Smith coordinates** (most common for atoms):
```
R = √(r₁² + r₂²)           — hyperradius
α = arctan(r₂/r₁)           — hyperangle, α ∈ [0, π/2]
θ₁, φ₁                       — direction of r₁
θ₂, φ₂                       — direction of r₂
```
The angle between r₁ and r₂ is θ₁₂ (not an independent coordinate — derived
from θ₁, φ₁, θ₂, φ₂).

**(b) Smith-Whitten democratic coordinates:**
```
R, Θ, Φ — body-frame hyperangles
α, β, γ — Euler angles for overall rotation
```
More symmetric but more complex.

**(c) Delves coordinates** (for scattering):
```
ρ = √(μ₁₂ r₁₂² + μ₃ r₃²)   — mass-weighted hyperradius
θ = arctan(...)               — arrangement channel angle
```

### For He (Central Field): Partial-Wave Reduction

For a central potential (He, not H₂), total angular momentum L is conserved.
The wavefunction decomposes as:

```
Ψ(r₁, r₂) = Σ_{l₁,l₂} R^{-5/2} F_{l₁l₂}^L(R) Y_{l₁l₂}^{LM}(r̂₁, r̂₂, α)
```

where the **coupled spherical harmonics** are:

```
Y_{l₁l₂}^{LM}(r̂₁, r̂₂) = Σ_{m₁m₂} ⟨l₁m₁ l₂m₂|LM⟩ Y_{l₁}^{m₁}(r̂₁) Y_{l₂}^{m₂}(r̂₂)
```

For L=0 (He ground state), the angular structure simplifies enormously:

```
Y_{ll}^{00}(r̂₁, r̂₂) = (-1)^l/√(2l+1) · P_l(cos θ₁₂)
```

So the L=0 wavefunction depends on only **three** coordinates: R, α, θ₁₂.

### Quantum Number Table for Hyperspherical Harmonics

The **hyperspherical harmonics on S⁵** are labeled by:

| Quantum Number | Range | Physical Meaning |
|:---|:---|:---|
| **K** (grand angular momentum) | 0, 1, 2, ... | Total hyperangular excitation |
| **l₁** | 0, 1, ..., K | Angular momentum of electron 1 |
| **l₂** | 0, 1, ..., K-l₁ | Angular momentum of electron 2 |
| **L** | \|l₁-l₂\|, ..., l₁+l₂ | Total angular momentum |
| **M** | -L, ..., L | z-projection |

Constraint: **K = 2n_α + l₁ + l₂** where n_α = 0, 1, 2, ... is the
hyperangular vibrational quantum number.

Degeneracy of shell K (for all L, ignoring symmetry): **(K+1)(K+2)²(K+3)/12**

### Mapping to Paper 0's (n, l, m) Structure

| S³ Atomic Lattice | S⁵ Hyperspherical |
|:---|:---|
| n (principal) | K (grand angular momentum) |
| l (angular momentum) | (l₁, l₂) pair |
| m (magnetic) | M (total z-projection) |
| 2n² states per shell | (K+1)(K+2)²(K+3)/12 per shell |
| SO(4) symmetry | SO(6) ⊃ SO(3) × SO(3) |

**Critical observation:** The S³ lattice has SO(4) symmetry (hydrogen atom).
The S⁵ hyperangular space has SO(6) symmetry, which decomposes as
SO(6) ⊃ SO(3)₁ × SO(3)₂ (angular momenta of electrons 1 and 2). For L=0
states, this further reduces because the coupled system locks l₁ = l₂.

---

## 4. Klar & Klar (1980): Threshold Laws

**Reference:** H. Klar and M. Klar, *J. Phys. B* **13**, 1057 (1980)

### Wannier Ridge and the Cusp

Klar & Klar analyzed the **Wannier ridge** — the configuration where both
electrons are equidistant from the nucleus (α = π/4) and on opposite sides
(θ₁₂ = π). This is the saddle point of the potential surface.

The effective potential at fixed R is:

```
C(α, θ₁₂) = -Z/cos α - Z/sin α + 1/√(1 - sin 2α cos θ₁₂)
```

The charge function C encodes all angular dependence. At the Wannier ridge
(α = π/4, θ₁₂ = π):

```
C_saddle = -2√2 Z + 1/√2

For Z=2 (He): C_saddle = -4√2 + 1/√2 ≈ -4.950
```

### Cusp Location in Hyperspherical Coordinates

The electron-electron coalescence (r₁₂ = 0) occurs when:

```
r₁₂ = R√(1 - sin 2α cos θ₁₂) = 0
→ sin 2α cos θ₁₂ = 1
→ α = π/4, θ₁₂ = 0
```

This is a **manifold** in the 5D hyperangular space, not a point singularity.
The cusp condition becomes a **boundary condition** on the channel functions
Φ_μ at (α = π/4, θ₁₂ = 0).

### Kato Cusp Condition in Hyperspherical Coordinates

```
lim_{r₁₂→0} (1/Ψ)(∂Ψ/∂r₁₂) = 1/2

In hyperspherical coordinates, this constrains the derivative of the
channel function at the coalescence manifold.
```

### Implication for GeoVac
The cusp is a **boundary condition on the hyperangular problem**, not a
singularity in the hyperradial coordinate. This is structurally analogous to
the boundary condition F(ξ=1) in the prolate spheroidal lattice. The lattice
can enforce it naturally through its edge structure.

---

## 5. Bartlett et al. (1935), Hylleraas (1929): Explicit Correlation

**References:**
- E.A. Hylleraas, *Z. Phys.* **54**, 347 (1929)
- J.H. Bartlett, *Phys. Rev.* **51**, 661 (1937)

### Historical Context

Hylleraas's variational ansatz for He:

```
Ψ(s, t, u) = e^{-ks} Σ c_{lmn} s^l t^{2m} u^n
```

where s = r₁+r₂, t = r₁-r₂, u = r₁₂. This achieves extraordinary accuracy
(6-term: 0.003% error) precisely because u = r₁₂ is an **explicit coordinate**
that naturally represents the cusp.

### Connection to Hyperspherical

In hyperspherical coordinates:
```
s = r₁ + r₂ = R(cos α + sin α)
t = r₁ - r₂ = R(cos α - sin α)
u = r₁₂ = R√(1 - sin 2α cos θ₁₂)
```

The Hylleraas variables (s, t, u) are smooth functions of (R, α, θ₁₂).
The hyperspherical approach **contains** Hylleraas-type correlations
implicitly through the coupled channel structure, without needing explicit
r₁₂ terms.

---

## 6. Summary: Quantum Numbers of the Hyperspherical System

### Complete Quantum Number Set

| Label | Name | Range | Conservation |
|:---|:---|:---|:---|
| **R** | Hyperradius | (0, ∞) | Not conserved (radial coord) |
| **K** | Grand angular momentum | 0, 1, 2, ... | Approximate (broken by 1/r_i, 1/r₁₂) |
| **l₁** | Ang. mom. electron 1 | 0, ..., K | Approximate (broken by 1/r₁₂) |
| **l₂** | Ang. mom. electron 2 | 0, ..., K-l₁ | Approximate (broken by 1/r₁₂) |
| **L** | Total angular momentum | \|l₁-l₂\|, ..., l₁+l₂ | **Exact** |
| **M** | z-projection of L | -L, ..., L | **Exact** |
| **S** | Total spin | 0 or 1 | **Exact** |
| **π** | Parity | ±1 | **Exact** |
| **ν** | Adiabatic channel index | 1, 2, ... | Approximate (non-adiabatic coupling) |

### For He Ground State (¹S)
- L = 0, M = 0, S = 0, π = +1 (all fixed)
- K labels the hyperangular excitation
- l₁ = l₂ (forced by L=0 coupling)
- The problem reduces to 3 effective coordinates: R, α, θ₁₂
- With L=0 partial-wave reduction: only α remains as the hyperangular variable
  (after expanding θ₁₂ dependence in Legendre polynomials P_l(cos θ₁₂))

### Mapping to GeoVac Hierarchy

```
Atomic S³ lattice:     nodes = (n, l, m)      ← SO(4) on S³
Prolate lattice:       nodes = (n_ξ, n_η, m)  ← separation in (ξ, η, φ)
Hyperspherical lattice: nodes = (R_i, ν)       ← adiabatic channels at R grid points

Edges:
  S³:     l ↔ l±1 (angular), n ↔ n±1 (radial)
  Prolate: ξ_i ↔ ξ_{i±1} (radial FD), η via Legendre coupling
  Hyper:   R_i ↔ R_{i±1} (hyperradial FD), ν ↔ ν' (non-adiabatic coupling)
```

The hyperspherical lattice is structurally the **same kind of object** — a
discrete graph where nodes carry quantum labels and edges encode physical
couplings. The key difference is that the "angular" problem is itself a
coupled PDE (in α, θ₁₂) rather than a separable ODE.
