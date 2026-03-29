# Algebraic PK Replacement: Data Flow Audit

**Date:** 2026-03-27
**Scope:** Read-only audit of PK pseudopotential data flow from core solve through Level 4 angular Hamiltonian injection.

---

## 1. Core Channel Data: What Does the Level 3 Solve Produce?

### Data Structure

The Level 3 angular eigenvalue problem is solved in `hyperspherical_angular.py:solve_angular()` (line 77). It returns:

```python
evals: ndarray (n_channels,)           # angular eigenvalues μ(R)
evecs: ndarray (n_channels, N)         # eigenvectors, N = (l_max+1) * n_alpha
```

The eigenvectors come directly from `scipy.linalg.eigh(H)` (line 194), transposed so each row is one eigenstate.

### Channel Expansion Coefficients c_l(α)

The eigenvector is a **discrete array on the α-grid**, indexed as:

```
psi[l * n_alpha : (l+1) * n_alpha]  →  c_l(α_i) for i = 0..n_alpha-1
```

where `l` ∈ [0, l_max] is the partial-wave index. These are **not** stored as continuous functions or expansion coefficients in some basis — they are FD grid values after Liouville substitution (`u_l(α) = sin(α) cos(α) φ_l(α)`).

### Extraction of Channel Weights

`core_screening.py:compute_core_density_algebraic()` (lines 160–316) extracts channel data by:

```python
_, vecs = solve_angular(R, Z, l_max, n_alpha, n_channels=1)
u = vecs[0]  # shape ((l_max+1) * n_alpha,)
ang_density = np.zeros(n_alpha)
for l_idx in range(n_l):
    ang_density += u[l_idx * n_alpha:(l_idx + 1) * n_alpha] ** 2
```

This gives `D(α; R) = Σ_l |c_l(α)|²`, the angular density at hyperangle α. The algebraic density method builds cubic splines of `D(α)` at each R, then evaluates at `α = arccos(r/R)` to map onto the radial coordinate.

### Shape Summary

| Quantity | Shape | Type |
|:---------|:------|:-----|
| Full eigenvector (ground state) | `((l_max+1) * n_alpha,)` | ndarray, float64 |
| Per-channel slice c_l(α) | `(n_alpha,)` | ndarray, float64 |
| Alpha grid | `(n_alpha,)` | Uniform FD: α_i = (i+1)h, h = π/(2(n_alpha+1)) |
| Channel weight ∫|c_l|²dα | scalar | float |

---

## 2. PK Injection Point

### Location

**File:** `level4_multichannel.py`
**Function:** `build_angular_hamiltonian()` (lines 752–987)
**PK block:** Lines 923–959

### How It Works

The PK is added as a **diagonal potential** V(α) on the finite-difference grid, independently per channel:

```python
# Lines 934-936: Compute angle-averaged PK potential on alpha grid
V_pk_e1, V_pk_e2 = compute_pk_pseudopotential(
    alpha_grid, rho_A, R_e, pk_potentials, rho_B=rho_B,
)

# Lines 946-959: Inject into diagonal of H
for ic, (l1, m1, l2, m2) in enumerate(channels_4):
    # ... determine w1, w2 from channel mode ...
    for i in range(n_alpha):
        ii = idx(ic, i)                    # = ic * n_alpha + i
        H[ii, ii] += R_e * w1 * V_pk_e1[i]  # electron 1
        H[ii, ii] += R_e * w2 * V_pk_e2[i]  # electron 2
```

### Matrix Structure at Injection

The angular Hamiltonian `H` has shape `(n_ch * n_alpha, n_ch * n_alpha)` with block structure:

```
       ch=0       ch=1       ch=2
ch=0 [ Diag+Kin | Nuc+Vee  | Nuc+Vee  ]
ch=1 [ Nuc+Vee  | Diag+Kin | Nuc+Vee  ]
ch=2 [ Nuc+Vee  | Nuc+Vee  | Diag+Kin ]
```

- **Diagonal blocks** (ch, ch): Tridiagonal kinetic + diagonal centrifugal + Liouville + core screening + **PK** + nuclear self-coupling
- **Off-diagonal blocks** (ch, ch'): Nuclear coupling + e-e coupling (both via Gaunt integrals)

The current PK adds **only to diagonal elements** `H[ii, ii]`, meaning it acts as a local potential in α-space within each channel's diagonal block. It does **not** couple channels.

### Upstream Data Flow

```
CoreScreening(Z).solve()
    → AbInitioPK(core, n_core)
        → pk_dict() = {'C_core': A, 'beta_core': B, 'atom': 'A'}
            → ComposedDiatomicSolver.pk_potentials = [pk_dict]
                → solve_level4_h2_multichannel(pk_potentials=...)
                    → compute_adiabatic_curve_mc(pk_potentials=...)
                        → solve_angular_multichannel(pk_potentials=...)
                            → build_angular_hamiltonian(pk_potentials=...)
                                → compute_pk_pseudopotential()  [lines 500-606]
                                → H[ii, ii] += R_e * w * V_pk[i]  [lines 956-959]
```

The `compute_pk_pseudopotential()` function (lines 500–606) evaluates the monopole (k=0) angle-average of `V_ECP(r) = C exp(-βr²)/r²`:

```
⟨V_ECP⟩(s, ρ, R_e) = C / (4 R_e² s ρ) × [E₁(β a²) − E₁(β b²)]
```

where `s = cos(α)` or `sin(α)`, `ρ = R/(2R_e)`, and E₁ is the exponential integral. Returns two arrays of shape `(n_alpha,)`: one for electron 1, one for electron 2.

---

## 3. Current l-Dependent Mode

### Implementation

In `build_angular_hamiltonian()` lines 945–959:

```python
pk_ch_mode = pk_potentials[0].get('channel_mode', 'channel_blind')

for ic, (l1, m1, l2, m2) in enumerate(channels_4):
    if pk_ch_mode == 'l_dependent':
        w1 = 1.0 if l1 == 0 else 0.0   # δ_{l1, 0}
        w2 = 1.0 if l2 == 0 else 0.0   # δ_{l2, 0}
    else:
        w1 = 1.0
        w2 = 1.0
    if w1 == 0.0 and w2 == 0.0:
        continue
    for i in range(n_alpha):
        ii = idx(ic, i)
        H[ii, ii] += R_e * w1 * V_pk_e1[i]
        H[ii, ii] += R_e * w2 * V_pk_e2[i]
```

### Behavior

- **channel_blind**: Full PK barrier applied to all channels uniformly.
- **l_dependent**: PK is **completely zeroed** for l > 0 channels. There is no scaling — it is a hard on/off per electron via δ_{l,0}.

For the (l1, l2) = (0, 2) channel: electron 1 (l1=0) gets full PK; electron 2 (l2=2) gets zero PK.
For (1, 1): both electrons get zero PK.
For (0, 0): both electrons get full PK.

### Configuration

Set via `pk_channel_mode='l_dependent'` in `ComposedDiatomicSolver` constructor (propagated to `pk_dict()` as `'channel_mode'` key).

---

## 4. Interface Specification for an Algebraic PK Projector

### What the Current Interface Accepts

The Level 4 angular Hamiltonian is a dense matrix of shape `(N, N)` where `N = n_ch × n_alpha`. The PK currently injects via diagonal elements only. However, the Hamiltonian `H` is just a NumPy array — any modification to `H[ii, jj]` at the right location would work.

### What an Algebraic Projector Would Provide

The exact Phillips-Kleinman projector is:

```
V_PK = (E_val − ε_core) |core⟩⟨core|
```

In the Level 4 channel basis, `|core⟩` must be expressed as a vector in the `(channel, α)` space:

```
|core⟩ → c_{l1,l2}(α_i)    shape: (n_ch * n_alpha,)
```

The projector becomes a **rank-1 matrix**:

```
M_PK[ii, jj] = (E_val − ε_core) × core_vec[ii] × core_vec[jj]
```

with shape `(n_ch * n_alpha, n_ch * n_alpha)`.

### Injection Point

The matrix should be added inside `build_angular_hamiltonian()` at lines ~960 (after current PK block, before return), as:

```python
H += (E_val - E_core) * np.outer(core_vec, core_vec)
```

Or equivalently, the projector could be passed as a new parameter alongside `pk_potentials` — e.g., `pk_projector: Optional[dict]` containing:

```python
{
    'core_vec': ndarray (n_ch * n_alpha,),  # |core⟩ in Level 4 basis
    'energy_shift': float,                   # (E_val − ε_core)
}
```

### Required Normalization

The core vector must be normalized: `∫ Σ_{l1,l2} |c_{l1,l2}(α)|² dα = 1`, which in FD approximation is `Σ_i core_vec[i]² × h_alpha ≈ 1`.

The R_e scaling factor that multiplies other potentials (line 958: `R_e * ...`) may or may not apply — this depends on whether the projector energy is in charge-function units or absolute energy units. The current PK is in charge-function units (divided by R_e), so the `R_e *` factor converts to energy. The algebraic projector `(E_val − ε_core)` is already in Ha, so the injection would be:

```python
H += R_e * (E_val - E_core) * np.outer(core_vec, core_vec)
```

if the Hamiltonian equation is `H u = μ u` where `μ` is in charge-function units (dimensionless × R_e).

---

## 5. Gap Analysis

### Key Question: Level 3 → Level 4 Basis Overlap

**The Level 3 and Level 4 angular bases use different channel labeling schemes.**

| Property | Level 3 (hyperspherical.py) | Level 4 (level4_multichannel.py) |
|:---------|:---------------------------|:---------------------------------|
| System | 2-electron atom (He, Li²⁺) | 2-electron molecule (H₂, LiH valence) |
| Channel label | Single index `l` | Pair `(l1, l2)` [or `(l1, m1, l2, m2)`] |
| Physical meaning | l = shared angular momentum of e-e pair | l1, l2 = individual electron angular momenta about molecular axis |
| Coordinate | Hyperangle α ∈ (0, π/2), single radius-ratio | Hyperangle α ∈ (0, π/2), plus molecular orientation |
| Nuclear potential | Spherical: −Z/cos α − Z/sin α | Molecular: split-region Legendre expansion |
| α-grid coupling | Coupled across channels by V_ee | Coupled across channels by V_nuc + V_ee |

### The Overlap Is NOT Trivially Available

The Level 3 channel `l` is the **total** angular momentum quantum number of the two-electron pair on S⁵ (the coupled angular momentum from the bipolar harmonics `Y_l^m(r̂₁) Y_l^m(r̂₂)` contracted to L=0).

The Level 4 channels `(l1, l2)` are **individual** angular momenta of each electron about the internuclear axis (Legendre polynomials `P_{l1}(cos θ₁) P_{l2}(cos θ₂)`).

These are related by a **Clebsch-Gordan transformation** — the Level 3 `l`-channel contains contributions from multiple Level 4 `(l1, l2)` pairs, and vice versa. Specifically, in the atomic limit (R → 0), the Level 4 `(l1, l2)` channel with l1 = l2 = l maps to the Level 3 channel `l`, but at finite R this mapping is mixed by the molecular field.

### What Is Missing

To construct `|core⟩` in the Level 4 basis, we need one of:

1. **Direct projection (recommended):** Solve the Level 3 angular problem at the same R_e (hyperradius) used in the Level 4 solve, extract the ground-state eigenvector `ψ_core(α)` as a function of α, then compute the overlap:

   ```
   core_vec[(l1,l2), α_i] = ∫ dΩ ψ_core(α, Ω) × Y_{l1}(Ω₁) Y_{l2}(Ω₂)
   ```

   This is a Clebsch-Gordan decomposition and can be done analytically via 3j symbols. The core (1s²) is dominated by `l=0`, which maps to `(l1, l2) = (0, 0)` plus small corrections from `(l, l)` pairs.

2. **Atomic-limit approximation:** For a 1s² core with Z_eff ≈ Z−0.3125 (screened), the core is >98% in the `l=0` channel. In the Level 4 basis, this maps almost entirely to `(0, 0)`. The projection would be approximately:

   ```
   core_vec ≈ δ_{(l1,l2),(0,0)} × u_core(α)
   ```

   where `u_core(α)` is the Level 3 ground-state α-wavefunction restricted to the `l=0` slice. This is the data already available from `solve_angular()`.

3. **Self-consistent Level 4 core solve:** Run the Level 4 angular solver itself with the core Z (no screening, no PK) and extract the ground state. This automatically gives `|core⟩` in the Level 4 `(l1, l2)` basis with no basis transformation needed. The tradeoff: this is a molecular solve (includes internuclear axis), which the atomic core doesn't physically "see" — but at the short distances where the core lives (r < 1 bohr for Li), the molecular field is negligible.

### Recommended Path

**Option 2 (atomic-limit approximation)** is the simplest and most physically justified:

- The 1s² core is 98%+ in the `l=0` channel (Paper 13, Section XII confirms weak channel coupling for ground state)
- In the Level 4 basis, `l=0` maps to `(l1, l2) = (0, 0)` with `m1 = m2 = 0`
- The α-dependence `u_core(α)` is already computed during `CoreScreening.solve()` → `solve_angular(R, Z, l_max, n_alpha, n_channels=1)`
- The α-grids must match (same n_alpha and grid spacing) between Level 3 and Level 4, or interpolation is needed

**What's needed to implement:**
1. Store the Level 3 ground-state eigenvector (currently discarded after density extraction in `compute_core_density_algebraic`)
2. Extract the `l=0` slice: `u_core_l0 = vecs[0][0:n_alpha]`
3. Normalize: `u_core_l0 /= sqrt(sum(u_core_l0**2) * h_alpha)`
4. Build Level 4 core vector: `core_vec = zeros(n_ch * n_alpha)`, then `core_vec[0:n_alpha] = u_core_l0` (channel 0 = (0,0))
5. Inject: `H += R_e * (E_val - E_core_per_electron) * outer(core_vec, core_vec)`

### Risk: α-Grid Mismatch

The Level 3 solve uses `n_alpha=100` (default in `CoreScreening`) while the Level 4 solve uses `n_alpha` set by the caller (typically 200 in `composed_diatomic.py`). If grids differ, cubic spline interpolation of the core wavefunction from Level 3 → Level 4 α-grid is needed. The α-grids are both uniform on (0, π/2) with Dirichlet BCs, so interpolation is straightforward.

### Risk: R_e Dependence

The Level 3 solve at a fixed R gives `ψ_core(α; R)`. In the Level 4 solve, R_e sweeps over a grid to build the adiabatic curve. The core eigenvector `|core⟩` must be evaluated at the **same R_e** as the Level 4 angular solve. Currently, `compute_core_density_algebraic()` already solves the angular problem at multiple R values and caches splines — this infrastructure can be reused.

---

## Summary

| Question | Answer |
|:---------|:-------|
| What does Level 3 produce? | Dense ndarray `(n_channels, (l_max+1)*n_alpha)`, FD grid values of `c_l(α)` |
| Where is PK injected? | `level4_multichannel.py:build_angular_hamiltonian()` lines 923–959, diagonal elements only |
| How does l-dependent PK work? | Hard δ_{l,0}: full PK for l=0 channels, zero for l>0 — no scaling |
| What shape should algebraic PK be? | Rank-1 matrix `(n_ch*n_alpha, n_ch*n_alpha)` or equivalently a vector `(n_ch*n_alpha,)` + energy scalar |
| Is the Level 3 → Level 4 overlap trivial? | **No.** Different channel labeling (single `l` vs pair `(l1,l2)`). But for 1s² core, the atomic-limit approximation maps `l=0` → `(0,0)` at >98% fidelity, making the overlap nearly trivial in practice. |
| What's missing? | (a) Persist Level 3 eigenvector (currently discarded), (b) α-grid interpolation if n_alpha differs, (c) R_e-dependent core vector evaluation |
