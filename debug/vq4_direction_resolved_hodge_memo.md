# VQ-4: Direction-Resolved Hodge Decomposition of the Fock Graph

## Summary

The Fock scalar graph edges are classified by their quantum-number
change (Dn, Dm) into two physically distinct channels:

- **T-channel** (radial): Dn = +/-1, Dm = 0 (inter-shell transitions)
- **L-channel** (angular): Dn = 0, Dm = +/-1 (intra-shell transitions)

Per-channel Hodge decompositions (sub-incidence B_mu, sub-Laplacians
L0_mu and L1_mu, photon propagator G_mu = L1_mu^+) are computed at
n_max = 2 and n_max = 3.

## Finding 0: Three Channels Reduce to Two

The naive classification into three channels (L+, L-, T) yields only
two non-empty channels. The L- channel (Dm = -1) is **empty** at both
n_max = 2 and n_max = 3. This is structural: the canonical edge
orientation (lower-index node -> higher-index node) always assigns
Dm = +1 to angular edges because the state ordering places m in
ascending order within each (n, l) sub-shell. Each undirected L-edge
represents BOTH the m -> m+1 and m+1 -> m transitions; the
"direction" is a labeling convention, not a physical distinction.

The physically meaningful decomposition is **T vs L** (radial vs angular).

## Finding 1: Edge Counts and Asymmetry

| n_max | V | E | E_T | E_L | E_T/E_L |
|:------|:--|:--|:----|:----|:--------|
| 2 | 5 | 3 | 1 | 2 | 0.500 |
| 3 | 14 | 13 | 5 | 8 | 0.625 |

At n_max = 2: 1 T-edge and 2 L-edges. The T-channel is the minority.
At n_max = 3: 5 T-edges and 8 L-edges. L-edges dominate because the
l = 2 shell contributes 4 angular edges (a path graph of length 5).

## Finding 2: Per-Channel Connectivity

**n_max = 2:**

- Channel T: 2 nodes touched, 1 connected components, connected = True
- Channel L+: 3 nodes touched, 1 connected components, connected = True

**n_max = 3:**

- Channel T: 9 nodes touched, 4 connected components, connected = False
- Channel L+: 11 nodes touched, 3 connected components, connected = False

The T-channel connects nodes across different n-shells (same l, m).
The L-channel connects nodes within the same (n, l) sub-shell across m.
Neither channel alone connects the full graph -- photon propagation
requires BOTH radial and angular hops.

## Finding 3: Spectral Anisotropy

**n_max = 2:** ANISOTROPIC (T and L+ spectra differ)

- T  nonzero L1 eigenvalues: [2.0]
- L+ nonzero L1 eigenvalues: [1.0, 3.0]

**n_max = 3:** ANISOTROPIC (T and L+ spectra differ)

- T  nonzero L1 eigenvalues: [1.0, 2.0, 2.0, 2.0, 3.0]
- L+ nonzero L1 eigenvalues: [0.381966, 1.0, 1.0, 1.381966, 2.618034, 3.0, 3.0, 3.618034]

The T and L channels have **different spectra**, confirming the
photon propagation is inherently anisotropic on the Fock graph.
The T-channel (radial) carries different eigenvalues than the
L-channel (angular). This is a structural feature of the graph
topology, not a dynamical effect.

## Finding 4: Propagator Additivity

Does the sum of per-channel propagators (embedded in the full edge
space) equal the full propagator?

| n_max | ||G_T + G_L - G_full||_F | Relative diff | Match? |
|:------|:------------------------|:--------------|:-------|
| 2 | 0.00e+00 | 0.00e+00 | True |
| 3 | 1.08e+00 | 3.37e-01 | False |

At n_max = 2: **YES** -- the propagator decomposes additively.
At n_max = 3: **NO** -- the propagator does NOT decompose additively (diff = 1.08e+00).

## Finding 5: Cross-Channel Coupling in G_full

The full propagator G_full can be partitioned into T-T, L-L, and
T-L blocks. If the T-L cross-block is nonzero, the photon
propagation mixes radial and angular channels.

| n_max | ||G_{TT}|| | ||G_{LL}|| | ||G_{TL}|| | Block-diag? |
|:------|:-----------|:-----------|:-----------|:------------|
| 2 | 0.500000 | 1.054093 | 0.000000 | True |
| 3 | 1.190333 | 2.955987 | 0.185179 | False |

At n_max = 2: G_full is **block-diagonal** in T/L channels.
The photon propagates independently in each direction channel.
At n_max = 3: G_full has **nonzero cross-channel coupling**.
The photon propagation mixes radial and angular channels.

## Finding 6: Per-Channel Propagator Traces

| n_max | Tr(G_T) | Tr(G_L) | Tr(G_full) | Tr(G_T)+Tr(G_L) |
|:------|:--------|:--------|:-----------|:----------------|
| 2 | 0.500000 | 1.333333 | 1.833333 | 1.833333 |
| 3 | 2.833333 | 6.666667 | 7.700000 | 9.500000 |

## Structural Interpretation

The direction-resolved Hodge decomposition reveals that the Fock graph
photon has two structurally distinct propagation channels:

1. **Radial (T)**: photon hops between shells (Dn = +/-1, same l and m).
   These are the edges that connect different energy levels.

2. **Angular (L)**: photon hops within a shell (Dm = +/-1, same n and l).
   These are the edges that rotate the magnetic quantum number.

The spectra are generically different (anisotropic photon). The full
propagator's T-L cross-block determines whether radial and angular
photon modes couple or propagate independently.

This gives the scalar graph photon effective "polarization" structure
from pure topology, without introducing explicit vector labels. The
anisotropy is a consequence of the Fock graph having two geometrically
distinct edge types built into its quantum-number lattice.

## Data Files

- `debug/data/vq4_direction_resolved_hodge.json` -- full numerical results
- `debug/vq4_direction_resolved_hodge.py` -- computation script

