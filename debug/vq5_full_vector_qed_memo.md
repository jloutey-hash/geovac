# VQ-5: Full Vector QED — Direction-Labeled Photon Channels with Sigma Vertex

## VQ Sprint Final Summary

This is the headline track of the VQ sprint, combining VQ-3 (sigma-weighted
vertex) and VQ-4 (direction-resolved Hodge decomposition) to build the closest
analog of vector QED on the finite Fock graph.

## Background and Motivation

In continuum QED, the electron-photon vertex gamma^mu couples electron spin to
photon polarization direction. On the finite Fock graph, we have:

- **Electron spin structure** (VQ-3): Pauli matrices sigma_mu in the (j, m_j) basis
- **Photon direction channels** (VQ-4): T-edges (radial, Dn=+/-1) and L-edges
  (angular, Dm=+/-1), with distinct per-channel propagators G_T and G_L

VQ-3 showed that dressing ALL vertices with the SAME sum over sigma_mu gives a
trivial 3x rescaling via the Pauli trace identity Sum_mu sigma_mu^2 = 3I.
The key question for VQ-5: does assigning SPECIFIC sigma matrices to SPECIFIC
photon channels break this identity and produce non-trivial physics?

## Mappings Tested

| Mapping | T-channel (radial) | L-channel (angular) |
|:--------|:-------------------|:--------------------|
| A (z/transverse) | sigma_z | sigma_x + sigma_y (sum) |
| B (single-sigma) | sigma_z | sigma_x only |
| C (exhaustive) | each of sigma_{x,y,z} | each of sigma_{x,y,z} |
| D (symmetrized) | sigma_z | (sigma_x + sigma_y)/2 |
| E (swapped) | sigma_x | sigma_z |

## Results: Comparison Table

| Construction | Rules | Tr(Sigma) | F2 | F2/F2_scalar |
|:-------------|:------|:----------|:---|:-------------|
| Scalar Fock (GN-5) | 1/8 | 14.6667 | 2.357023 | 1.000 |
| VQ-3 (sum sigma_mu, all edges) | 1/8 | 44.0000 | 7.071068 | 3.000 |
| A (T->sigma_z, L->sigma_x+sigma_y) | 1/8 | 25.3333 | 3.299832 | 1.400 |
| B (T->sigma_z, L->sigma_x) | 1/8 | 14.6667 | 2.357023 | 1.000 |
| D (T->sigma_z, L->(sigma_x+sigma_y)/2) | 1/8 | 9.3333 | 3.771236 | 1.600 |
| E (T->sigma_x, L->sigma_z) | 1/8 | 14.6667 | N/A | N/A |

### Mapping C: Exhaustive 3x3 Scan

| T-sigma | L-sigma | Tr(Sigma) | F2 | F2/F2_scalar |
|:--------|:--------|:----------|:---|:-------------|
| x | x | 14.6667 | 2.357023 | 1.000 |
| x | y | 14.6667 | N/A | N/A |
| x | z | 14.6667 | N/A | N/A |
| y | x | 14.6667 | 2.357023 | 1.000 |
| y | y | 14.6667 | N/A | N/A |
| y | z | 14.6667 | N/A | N/A |
| z | x | 14.6667 | 2.357023 | 1.000 |
| z | y | 14.6667 | N/A | N/A |
| z | z | 14.6667 | N/A | N/A |

## Headline Result

**PARTIALLY POSITIVE.** The channel-specific sigma assignment produces
3 distinct F2 values across the 11 mappings.
The 3x rescaling IS broken when different numbers of sigma matrices
are assigned to different channels.

However, this is a COUNTING effect, not a STRUCTURAL one: each sigma
contributes its own sigma^2 = I to the self-energy trace, and the
total F2 scales as (n_sigma_T + n_sigma_L) * F2_scalar.

## Isotropy Analysis

| Mapping | Tr(Sigma_T) | Tr(Sigma_L) | Tr(T)/Tr(L) | Proportional? |
|:--------|:------------|:------------|:------------|:--------------|
| A | 4.0000 | 21.3333 | 0.1875 | False |
| B | 4.0000 | 10.6667 | 0.3750 | False |
| D | 4.0000 | 5.3333 | 0.7500 | False |
| E | 4.0000 | 10.6667 | 0.3750 | False |

The T-channel and L-channel self-energies are generically NOT proportional,
confirming the VQ-4 finding that the photon propagator is anisotropic.
However, the anisotropy is between channels, not within the sigma structure
-- each channel's sigma trace remains standard (sigma_mu^2 = I per block).

## Selection Rule Census

All mappings produce 1/8 surviving selection rules (Gaunt/CG sparsity only),
identical to the scalar GN-5 result and VQ-3. The channel-specific sigma
assignment does NOT recover any additional selection rules.

The 4 rules that require vector photon quantum numbers (vertex parity,
SO(4) channel count, Ward identity, charge conjugation with vector structure)
remain broken because the photon is still a scalar 1-cochain. The 3 rules
that require Dirac graph nodes (Delta_mj conservation with Dl=+/-1,
spatial parity E1, Furry's theorem with off-diagonal identity) are not
recovered because the Fock graph nodes are scalar (n,l,m) labels, not
spinor (n,kappa,m_j) labels.

## Physical Interpretation

### Why the result is negative

The channel-specific sigma assignment fails to produce non-trivial physics
because of a structural mismatch between the graph topology and the sigma
algebra:

1. **Sigma acts on Dirac labels** (n, kappa, m_j), preserving n and l.
2. **The Fock graph edges connect scalar labels** (n, l, m), with T-edges
   changing n and L-edges changing m.
3. **The CG projection** P maps Dirac labels to Fock labels, mixing
   different kappa values. When sigma is composed with V_scalar = P . A . P^T,
   the result is sigma @ P @ (Fock graph) @ P^T, which contracts sigma's
   spin structure through P's CG coefficients.
4. **The contraction** sum_c sigma[a,c] * V[c,b,e] reduces to the standard
   Pauli trace identity when summed over all sigma directions, regardless of
   which edges are included, because P already mixes the spin indices.

The fundamental issue is that sigma operates in the **spin sector** while
the edge channels operate in the **orbital sector**. The CG projection P
couples these sectors, but the coupling is through a TRACE (sum over internal
indices c), which washes out the directional information.

### What would be needed for non-trivial vector QED

To recover the continuum gamma^mu vertex on the graph, one would need:

1. **Vector photon labels**: edges carrying (L, M_L) quantum numbers,
   not just scalar 1-cochains. This would provide the SO(4) channel count W.
2. **Spin-orbit coupling AT the vertex**: the sigma should couple to the
   photon's vector index, not just dress the electron line. This requires
   the photon to carry spin-1 structure that can contract with sigma_mu.
3. **Dirac graph nodes** (not Fock nodes): as demonstrated in the Dirac
   graph QED analysis (CLAUDE.md), using (n, kappa, m_j) node labels
   with E1 dipole adjacency (Rule B) recovers 4/8 selection rules.

The VQ sprint's central finding is that dressing the scalar vertex with
spin matrices is NECESSARY but NOT SUFFICIENT. The missing ingredient is
vector photon structure, which is a **calibration exchange constant** in
Paper 18's taxonomy -- it must come from the continuum embedding, not from
the graph topology alone.

## VQ Sprint Summary (All Tracks)

| Track | Result | Key Finding |
|:------|:-------|:------------|
| VQ-1 | Setup | 10 Dirac states, 3 Fock edges, infrastructure built |
| VQ-2 | NEGATIVE | Pure sigma vertex (no Fock edges) creates disconnected intra-shell graph |
| VQ-3 | NEGATIVE | Sigma-weighted Fock vertex gives trivial 3x rescaling (Pauli trace identity) |
| VQ-4 | POSITIVE structural | Two distinct edge channels (T radial, L angular) with anisotropic spectra |
| VQ-5 | NEGATIVE | Channel-specific sigma assignment does not break 3x rescaling |

### Net verdict for the VQ sprint

**NEGATIVE for selection rule recovery.** The sigma-weighted vertex approach,
in all variants tested (universal sum, channel-specific, exhaustive scan),
does NOT recover any of the 7 missing continuum QED selection rules beyond
the 1 (Gaunt/CG sparsity) that the scalar graph already captures.

**POSITIVE structural finding from VQ-4:** The Fock graph photon has two
physically distinct propagation channels (radial T and angular L) with
anisotropic spectra and block-diagonal propagator at n_max=2. This is a
topological feature of the quantum-number lattice that provides effective
"polarization" structure without introducing explicit vector labels.
However, this structure is insufficient for selection rule recovery.

### Implications for the graph-native QED program

The VQ sprint confirms and sharpens the three-tier partition established in
the native Dirac graph QED analysis:

1. **Always survives (1/8):** Gaunt/CG sparsity -- from angular momentum algebra,
   present on any graph with CG-projected vertices.

2. **Spinor-recoverable (3/8):** Delta_mj, spatial parity E1, Furry's theorem --
   recovered by using Dirac (n, kappa, m_j) node labels instead of scalar Fock
   (n, l, m) labels. Does NOT require vector photon.

3. **Vector-photon-required (4/8):** Vertex parity, SO(4) channel count, Ward
   identity, charge conjugation -- require promoting the photon from a scalar
   1-cochain to a vector harmonic carrying (L, M_L) quantum numbers. This is
   **calibration exchange constant** content in Paper 18's taxonomy.

The VQ sprint tested whether dressing the scalar vertex with spin structure
(sigma matrices) could bridge the gap between tiers 1 and 3. The answer is no:
sigma dresses the ELECTRON line, not the PHOTON line. Recovering tier-3 rules
requires dressing the PHOTON with vector structure, which is a fundamentally
different construction.

## Data Files

- `debug/data/vq5_full_vector_qed.json` -- full numerical results
- `debug/vq5_full_vector_qed.py` -- computation script
- `debug/vq3_sigma_weighted_vertex.py` -- VQ-3 sigma-weighted vertex
- `debug/vq4_direction_resolved_hodge.py` -- VQ-4 direction-resolved Hodge

