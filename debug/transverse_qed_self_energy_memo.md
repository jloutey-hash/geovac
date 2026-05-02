# Transverse QED Self-Energy: Plaquette-Recovered Ground-State Protection

**Date:** 2026-05-01  
**Track:** GN-QED extension (transverse sector)  
**Verdict:** POSITIVE — ground-state structural zero recovered from graph topology

## Summary

The transverse photon propagator G_T = (d₁d₁ᵀ)⁺, built from Paper 30's Wilson
gauge plaquettes, recovers the continuum QED ground-state protection theorem
(Paper 28 Theorem 4: Σ(n_ext=0) = 0) which was broken on the scalar Fock graph
(pendant-edge theorem: Σ(GS) = 2(n_max-1)/n_max).

## Results

| n_max | Σ_scalar(GS) | Σ_transverse(GS) | Plaquettes | Co-exact rank |
|:-----:|:------------:|:-----------------:|:----------:|:-------------:|
| 3     | 4/3          | 0 (exact)         | 4          | 2             |
| 4     | 3/2          | 0 (exact)         | 33         | 8             |

Both verified to machine precision (< 1e-16).

## Mechanism: Topological Protection via Pendant Edge

The ground state |1, κ=-1, m_j=±1/2⟩ maps via CG projection to the Fock node
n=1 (the unique node in the first shell). This node has degree 1 on the Fock
graph — it is a **pendant vertex** connected to the rest of the graph by a
single edge e₀.

A pendant edge cannot participate in any closed walk (plaquette) because you
would need to leave the pendant vertex and return to it, but it has only one
neighbor. Therefore:

    d₁[e₀, f] = 0  for all faces f

This means the e₀ row of d₁ is identically zero, so:

    (d₁d₁ᵀ)[e₀, :] = 0  and  (d₁d₁ᵀ)[:, e₀] = 0

Therefore G_T[e₀, :] = G_T[:, e₀] = 0, and the self-energy sum

    Σ[a, b] = Σ_{e,e'} V[a,:,e] · G_T[e,e'] · V[b,:,e']ᵀ

vanishes whenever V couples through e₀. Since the GS node connects ONLY via
e₀, Σ(GS) = 0 identically.

This is a **theorem**, not a numerical observation: for any graph where the
ground state is a pendant vertex, the co-exact (transverse) self-energy
vanishes on the ground state.

## Contrast with Scalar (Longitudinal) Propagator

The scalar propagator G_scalar = L₁⁺ = (BᵀB)⁺ does NOT have this property.
The pendant edge e₀ participates fully in L₁ (it has L₁[e₀,e₀] = 2 from
its two endpoint degrees), so G_scalar[e₀, e₀] = (n_max-1)/n_max > 0, and
Σ_scalar(GS) = 2(n_max-1)/n_max (the pendant-edge theorem from GN-5).

## Continuum Interpretation

In continuum QED on S³, the structural zero Σ(n_ext=0) = 0 arises because
vertex parity forces n₁ + n₂ + q to be odd, and with n₂ = n_ext = 0 this
requires 2n₁ to be odd — impossible. The root cause is that the vertex
parity selection rule restricts the photon mode sums.

The graph recovers the SAME protection through a different but structurally
equivalent mechanism: the transverse sector (co-exact Hodge) topologically
decouples from pendant vertices. The "vertex parity" of the continuum has its
graph analog in the "plaquette exclusion of pendant edges."

## What Does NOT Improve

The transverse propagator does not improve the other selection rule tests:

- **Diagonal dominance:** 0.56-0.57 (transverse) vs 0.64-0.65 (scalar)
- **Cross-shell suppression:** 41-43% within-shell (transverse) vs 87-88% (scalar)

This is expected and correct: transverse photons carry angular momentum (L≥1)
and should mediate inter-shell transitions. The scalar propagator's high
within-shell fraction is an artifact of the longitudinal sector (gauge modes
that should be removed by gauge fixing).

## Refined Selection Rule Partition

Previous (from VQ sprint):
- 1 always survives (Gaunt/CG sparsity)
- 3 spinor-recoverable (Dirac graph nodes)
- 4 vector-photon-required

New (from transverse sprint):
- 1 always survives (Gaunt/CG sparsity)
- 3 spinor-recoverable (Dirac graph nodes)
- **1 plaquette-recoverable** (GS structural zero from co-exact topology)
- 3 vector-quantum-number-required (SO(4) channel count, Ward identity, charge conjugation)

The boundary between "graph-intrinsic" and "calibration" has shifted: GS
protection is graph-intrinsic via the Hodge decomposition, not requiring
any photon quantum numbers.

## Paper 18 Taxonomic Placement

This result sharpens the three-tier structure:
- **Tier 1 (graph-intrinsic):** Gaunt sparsity + spinor quantum numbers + co-exact topology
- **Tier 2 (calibration):** photon (L, M_L) quantum numbers within the co-exact sector
- **Tier 3 (embedding):** continuum spectral density matching (projection constants C_VP, C_SE, C_F2)

## Structural Significance

Paper 30's SU(2) Wilson gauge plaquettes were originally motivated by lattice
gauge theory formalism. This result gives them a concrete QED role: they define
the face structure whose boundary operator separates longitudinal (gauge) from
transverse (physical) photon sectors, and this separation alone — without any
labeling of modes — recovers one continuum selection rule.

Papers 25 and 30's Hodge decomposition (L₁ = BᵀB is the exact/longitudinal
sector) is now COMPLEMENTED by d₁d₁ᵀ as the co-exact/transverse sector, and
only the latter is the physical photon propagator.

## Key Files

- `debug/transverse_qed_self_energy.py` — computation script
- `debug/data/transverse_qed_self_energy.json` — numerical results
- `debug/fock_photon_projection.py` — prerequisite establishing L₁ = longitudinal
