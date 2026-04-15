# Spin-Orbit Coupling Design Memo (Track T2)

**Sprint:** Dirac-on-S³ Tier 2, Track T2.
**Deliverable:** `geovac/spin_orbit.py`, `tests/test_spin_orbit.py`.
**Consumes:** T1 closed-form matrix elements (`geovac/dirac_matrix_elements.py`).

## Goal

Assemble the Breit–Pauli leading-order spin-orbit Hamiltonian `H_SO` on
the Dirac-on-S³ basis (D1) using the T1 closed forms for `⟨L·S⟩` and
`⟨1/r³⟩`. All matrix elements are produced algebraically — no
numerical quadrature, no diagonalisation. The result is diagonal in
(n, κ, m_j).

## Derivation

In atomic units the Breit–Pauli spin-orbit operator for one electron
in a central potential `V(r)` is

```
H_SO = ξ(r) · L·S ,    ξ(r) = (1/(2 m² c²)) · (1/r) · dV/dr .
```

For Coulomb `V(r) = −Z/r` and α = 1/c (α is the fine-structure
constant),

```
dV/dr = Z/r² ,    ξ(r) = Z α² / (2 r³) .
```

Hence

```
H_SO = (Z α² / 2) · (1/r³) · L·S .     (*)
```

In the κ-native basis, T1 provides both factors:

| Operator | Closed form | Source |
|:---------|:------------|:-------|
| ⟨κ, m_j\|L·S\|κ, m_j⟩ | −(κ + 1)/2 | T1 Szmytkowski reduction |
| ⟨n, l\|1/r³\|n, l⟩   | Z³ / [n³ · l(l+½)(l+1)] | Bethe–Salpeter §3 |

Both are diagonal. L·S commutes with J_z so H_SO is m_j-independent.
Substituting into (*) gives the closed-form matrix element

```
⟨n, κ, m_j | H_SO | n, κ, m_j⟩
    = (Z α² / 2) · [−(κ+1)/2] · Z³ / [n³ · l(l+½)(l+1)]
    = −Z⁴ α² · (κ+1) / [4 · n³ · l(l+½)(l+1)] ,     (**)
```

where `l = kappa_to_l(κ)` and `(l + ½) = (2l+1)/2` in the
rational-arithmetic implementation. Equation (**) is pure
sympy — a single `sp.Expr`, no numerical evaluation.

## The l = 0 (Kramers) cancellation

For s-states only κ = −1 is allowed. Then `κ + 1 = 0` so the
numerator of (**) vanishes, and `H_SO = 0` identically for every
s-state at every n and Z. This is the Foldy–Wouthuysen / Kramers
cancellation: the apparent divergence of `⟨1/r³⟩` at l = 0 is exactly
matched by the zero of `⟨L·S⟩`, and the Breit–Pauli SO operator is
well-defined on s-states with no regularisation needed.

**Implementation.** `so_diagonal_matrix_element(n, κ, …)` short-circuits
the `l == 0` case and returns `Integer(0)` *before* any attempt to
evaluate `inverse_r_cubed_hydrogenic(n, 0)` (which would raise
ValueError). This is explicit rather than relying on sympy to
simplify `0 · oo`. See `tests/test_spin_orbit.py::TestKramersL0`.

## Z⁴ scaling

For fixed (n, κ) with l ≥ 1, equation (**) is exactly proportional to
Z⁴. For Tier-2 targets:

| Atom | Z | Z⁴ | H_SO relative to hydrogen |
|:-----|:-:|:---|:--------------------------|
| H    | 1  | 1          | 1× |
| Li   | 3  | 81         | 81× |
| Be   | 4  | 256        | 256× |
| Sr   | 38 | 2 085 136  | ~2.1 × 10⁶ × |

This is why SrH needs a relativistic treatment: the Breit–Pauli SO
energy at 2p₃/₂ in hydrogen is α²/96 ≈ 5.6 × 10⁻⁷ Ha ≈ 0.12 cm⁻¹,
whereas in Sr it is ≈ 1.2 × 10⁰ Ha ≈ 2.5 × 10⁵ cm⁻¹. (The actual
Sr SO is smaller — the full Dirac-operator γ = √(1−(Zα)²) correction
and valence screening reduce this — but the order-of-magnitude
message stands.) `verify_z4_scaling(n, κ, Z_values)` evaluates (**)
symbolically at requested Z and returns the map; the corresponding
regression test confirms the ratio is exactly (Z/Z_ref)⁴.

## Hydrogen fine-structure benchmark

The standard benchmark is the 2p doublet:

| State | κ | ⟨L·S⟩ | ⟨1/r³⟩ (Z=1) | H_SO |
|:------|:-:|:------|:--------------|:-----|
| 2p₃/₂ | −2 | +1/2 | 1/24 | +α²/96 |
| 2p₁/₂ | +1 | −1   | 1/24 | −α²/48 |

The spin-orbit splitting is

```
Δ_SO(2p) = H_SO(2p_{3/2}) − H_SO(2p_{1/2})
         = α²/96 − (−α²/48)
         = α²/96 + 2α²/96  =  3α²/96  =  α²/32 .
```

**Note on the α²/16 quote.** The full Dirac-equation fine-structure
splitting in hydrogen at n = 2 is α²·Z⁴·(1/16 − 1/12) = ... (and
reduces to α⁴/32·m_e in cgs/SI; α²/16 Hartree is sometimes quoted
in short-hand for the *total* 2p–2s-summed fine-structure level
energy relative to the bare Schrödinger eigenvalue, which is NOT
the SO splitting alone). Equation (**) is the spin-orbit piece
only — no Darwin term, no mass-velocity (relativistic KE)
correction. The closed form (**) gives α²/32 for Δ_SO(2p) and
this is the number the tests pin.

Bringing the Darwin and mass-velocity terms in (T6/relativistic
extension) would recover the full Dirac fine-structure ladder at
α⁴·Z⁴ order. That is explicitly out of scope for T2 and deferred
to either a future relativistic-correction track or to replacement
of (**) by the full Martínez-y-Romero Dirac-Coulomb radial
expectation ⟨r⁻³⟩_D (which carries γ = √(1−(Zα)²) corrections).

## Paper 18 classification (provisional)

Per Paper 18's exchange-constant taxonomy:

- The α² prefactor is **spinor-intrinsic content**: it appears
  first-order in α², carries no π, no transcendentals beyond α
  itself, and is present already at the one-electron Breit–Pauli
  reduction. It is the same α² that appears in the classical
  Thomas-precession factor of ½ absorbed into the Breit–Pauli
  form.
- The radial factor Z⁴/[n³ · l(l+½)(l+1)] is a **pure rational in
  (Z, n, l)**. No π, no Riemann ζ, no e^a E₁(a). It comes from the
  Bethe–Salpeter hydrogenic ⟨1/r³⟩ which T1 registered as
  algebraic.
- The operator is **diagonal in (n, κ, m_j)**; no basis rotation or
  quadrature is incurred.

Therefore the entire H_SO block is **algebraic** in the Paper 18
taxonomy. T5 will formalise the classification row for Section 12
(Algebraic Registry, Level 5 Composed) with this as input.

## API summary

| Function | Purpose |
|:---------|:--------|
| `so_diagonal_matrix_element(n, κ, Z, α)` | One eigenvalue, symbolic or numeric. |
| `build_so_hamiltonian_block(n_max, Z, α)` | Full diagonal block over `iter_dirac_labels(n_max)`. |
| `verify_z4_scaling(n, κ, Z_values, α)` | Evaluate across Z; regression helper. |
| `SOHamiltonianBlock` (dataclass) | `.diag: {DiracLabel → Expr}`, `.eigenvalues()`, `.as_vector(labels)`. |

## What T3 consumes

T3 (composed-pipeline integration) needs:

1. `so_diagonal_matrix_element(n, κ, Z)` — direct call for a single
   orbital's SO contribution, to be added as a diagonal update to
   the one-body matrix h₁ in the composed block Hamiltonian.
2. The fact that H_SO is diagonal in (n, κ, m_j), so integration
   into the composed pipeline is a one-body diagonal update. No
   off-block couplings, no new selection rules.
3. The α² factor stays symbolic until the user substitutes a
   numerical α, enabling T6's taxonomy classification to reason
   about H_SO's exchange-constant content without premature
   binding.

## Deferred / out of scope

- **Full Dirac-Coulomb radial factor** (γ corrections, Martínez-y-Romero
  recursions): deferred to a later relativistic-extension track. At
  Z α ≪ 1 the Breit–Pauli form is accurate to < 1 %.
- **Darwin + mass-velocity terms**: needed to reproduce the full
  α⁴ fine-structure ladder; out of scope for T2.
- **Zeeman SO-extension** (`so_zeeman_correction`): the task
  description flagged this as an optional stretch goal. On inspection,
  a Zeeman term `g_s B · S` is algebraic but couples to m_j and is
  itself a one-body diagonal update of the same structural type as
  H_SO; combined with H_SO it requires diagonalising in m_j for each
  (n, κ) block. That's still algebraic but an extra implementation
  surface we do not need for T3. **Deferred**; not implemented.
