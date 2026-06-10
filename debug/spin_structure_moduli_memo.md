# Spin-structure moduli on the Layer-1 substrate, and the chirality-reflection obstruction

**Date:** 2026-06-08
**Drivers:** `debug/spin_config_symmetry_decomposition.py` (scalar Fock graph),
`debug/spin_config_dirac_substrate.py` (Dirac-S³ graph).
**Verification:** `tests/test_spin_structure_obstruction.py` (bit-exact).
**Origin:** exploratory conversation — "can we put spin on hemispheres; how do
symmetries (de)compose on the Layer-1 substrate; can different projections carry
different spin configurations?"

## Framing

"Arrangements that preserve coupling" on the bare graph = **flat ℤ₂ connections**
= the cycle space H¹(G; F₂) = ker(B mod 2), dimension β₁. These are exactly the
**inequivalent spin structures** on the substrate (a spin structure is a choice
of sign-transport around loops — an element of H¹(·; ℤ₂)). There are 2^β₁ of them.
A graph symmetry σ (node permutation preserving adjacency) acts on this space as an
F₂-linear map; the configs that *also* preserve σ are the fixed subspace
ker(M_σ − I), so #symmetric = 2^dim(fixed). All arithmetic is exact over F₂ (π-free,
Layer-1).

## Scalar Fock graph (FockGraphHodge: L±/T± ladders, Δl = 0)

The graph factorises into one component per ℓ-sector; each ℓ-component is a 2D grid
P_{(n_max−ℓ)} × P_{(2ℓ+1)} in the (n, m) plane, so the cycles are literally the grid
plaquettes. β₁ = Σ_ℓ (n_max−ℓ−1)(2ℓ).

| n_max | β₁ | coupling-preserving (2^β₁) | Hopf m→−m symmetric | radial-inv. symmetric |
|:-----:|:--:|:--:|:--:|:--:|
| 2 | 0 | 1 | 1 | 1 |
| 3 | 2 | 4 | **2** | 4 |
| 4 | 8 | 256 | **16** | 64 |

- **Hopf reflection m→−m** acts on the two n=3 plaquettes as the pure swap
  `[[0,1],[1,0]]`: 2 symmetric structures survive — `00` (no flux) and `11` (flux
  through both plaquettes); `10`/`01` (single-plaquette flux) form one broken orbit.
  *Single-plaquette flux = "spin on one hemisphere" = the symmetry-breaking sector.*
- **Radial inversion n→reflect** is a hidden combinatorial automorphism that acts as
  the identity on the n=3 cycle space (orientation-only, invisible mod 2) — a physical
  symmetry constrains the spin structure; a hidden one need not. (At n_max=4 it begins
  to act, fixed dim 6 of 8.)
- #Hopf-symmetric = 2^(number of plaquette orbits under m→−m).

**Submission confirmed at the scalar level:** ≥2 inequivalent symmetry-respecting
spin structures from n_max=3 on.

## Dirac-S³ graph (DiracLattice, mode='atomic', E1 dipole Rule B, Δl=±1, |Δm_j|≤1)

| n_max | β₁ | coupling-preserving (2^β₁) | m_j-symmetric (2^fixed) |
|:-----:|:--:|:--:|:--:|
| 2 | 11 | 2,048 | 64 (fixed 6) |
| 3 | 79 | ~6.0×10²³ | ~1.1×10¹² (fixed 40) |
| 4 | 253 | ~1.4×10⁷⁶ | ~1.7×10³⁸ (fixed 127) |

Two reflections behave oppositely:

**Magnetic σ_mj : (n,κ,m_j) → (n,κ,−m_j).** Graph automorphism at every n_max
(adjacency depends only on |Δm_j|). **Fixed-point-free on nodes** — every state pairs
with its m_j-reverse, none self-paired, because half-integer m_j is never 0. This is
**Kramers degeneracy** straight from the combinatorics (contrast: the scalar m=0 nodes
were fixed points). The moduli of symmetric spin structures *grows enormously* vs the
scalar graph (2 → 10¹² at n_max=3): the coupling did not collapse the moduli, it
enlarged it.

**Chirality σ_χ : χ → −χ (genuine spin reflection).** NOT a node permutation at any
n_max — and not because it is "broken," but because the two chirality sectors have
**unequal sizes**:

| n_max | #χ=+1 (κ<0) | #χ=−1 (κ>0) | imbalance |
|:-----:|:--:|:--:|:--:|
| 2 | 8 | 2 | 6 |
| 3 | 20 | 8 | 12 |
| 4 | 40 | 20 | 20 |

Per orbital block (n, ℓ): the χ=+1 multiplet (κ=−(ℓ+1), j=ℓ+½) has **2ℓ+2** states,
the χ=−1 multiplet (κ=+ℓ, j=ℓ−½) has **2ℓ** — a mismatch of exactly **+2, every block,
never cancelling** (ℓ=0 gives 2 vs 0). Global imbalance = **n_max(n_max+1)** =
2 × (number of orbital blocks). There is no bijection between the chirality halves, so a
chirality-reflecting node permutation **does not exist**.

This +2-per-block chiral asymmetry is the discrete **Dirac index / spectral asymmetry**
of the Camporesi–Higuchi spinor bundle (positive-chirality eigenvalue +(n+3/2),
negative −(n+3/2); the count asymmetry is η-invariant-flavoured). It is **robust**: pure
node-counting, independent of adjacency rule (A vs B) and mode (atomic vs s3). The
moduli counts (2^β₁) are Rule-dependent; the chirality obstruction is not.

## Why this matters (corpus connection)

This is the **one-body, graph-level structural root** of the three relativistic-ℤ₂
tapering sprints (κ-parity v3.92.0, m_j-parity direct v3.94.0, rotated v3.95.0), all
NEGATIVE. Those sprints tried to build a ℤ₂ stabiliser from a chirality/spin sign-count;
the memory `gaunt_parity_protects_hopf_z2_not_relativistic.md` gives the *Hamiltonian-level*
reason ("M_J is a sum, not a parity"; jj-coupled ERI sign $(-1)^{j_a+j_b+j_c+j_d}$ not
forced even). The present finding gives the deeper *substrate-level* reason: the chirality
reflection is not even a graph automorphism, because the spinor bundle's two chirality
halves differ in dimension by 2 per orbital block. You cannot taper with a reflection that
isn't a permutation. The magnetic m_j reflection, by contrast, IS a clean Kramers
automorphism and supports an enormous symmetric spin-structure moduli.

Through-line of the conversation: "put spin on hemispheres" → on the scalar graph the
hemisphere flux is the symmetry-breaking sector; on the Dirac graph genuine spin cannot
be hemisphere-reflected at all, obstructed by an exact +2-per-block chiral asymmetry =
the framework's discrete Dirac index. The Hopf-bundle non-triviality I gestured at
conceptually is, concretely, this counting obstruction.

## Scope / honesty

- Flat-ℤ₂ spin *structures* at the one-body / substrate level; the full two-body
  Hamiltonian selection is a layer below (this explains, does not contradict, the
  Hamiltonian-level negatives).
- Dirac numbers are Rule-B, mode='atomic'. The **chirality imbalance is rule- and
  mode-independent** (node counting); the moduli sizes are not.
- π-free, bit-exact F₂ throughout — Layer-1 skeleton, consistent with the
  discrete-for-skeleton discipline.

## Recommended corpus capture (NOT yet applied — for PI greenlight)

Short observation in **Paper 29** (Dirac-S³ graph / Hopf-ℤ₂ block section): the Dirac
graph admits a fixed-point-free Kramers automorphism (m_j→−m_j) but no chirality
reflection, the obstruction being the 2ℓ+2 vs 2ℓ spinor-bundle imbalance (discrete Dirac
index n_max(n_max+1)). Reads well to the NCG audience as a graph-theoretic Dirac-index
statement. Would need its `tests/test_spin_structure_obstruction.py` gate (done).
