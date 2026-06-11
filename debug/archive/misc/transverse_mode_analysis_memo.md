# Transverse Mode Analysis: Can Co-Exact Eigenvectors Carry Photon Quantum Numbers?

**Date:** 2026-05-01  
**Track:** GN-QED extension (transverse mode-resolved)  
**Verdict:** STRUCTURAL NEGATIVE — plaquette topology gives correct eigenvalues but wrong support

## Question

After the transverse propagator G_T = (d₁d₁ᵀ)⁺ recovered the GS structural zero
(pendant-edge theorem), the question is: can the individual co-exact eigenvectors
be labeled by photon angular momentum q, such that the mode-resolved self-energy
enforces the remaining 3 selection rules (SO(4) channel count, Ward identity,
charge conjugation)?

## Method

1. Diagonalize d₁d₁ᵀ to get eigenvectors {v_k} with eigenvalues {μ_k}
2. Assign q_k by nearest continuum eigenvalue n(n+2)
3. Decompose Σ_T = Σ_k Σ_k into per-mode contributions
4. Test whether mode k with assigned q_k only couples electron states
   (n₁, n₂) satisfying the triangle inequality |n₁-n₂| ≤ q ≤ n₁+n₂-2
5. Test Ward identity via ||[Σ_T, H₀]|| / ||Σ_T||

## Results

### Eigenvalue Matching

| n_max | Modes | Best match (residual=0) | Worst residual |
|:-----:|:-----:|:-----------------------:|:--------------:|
| 3     | 2     | q=3 (exact 15.0)        | q=2 (+1.0)     |
| 4     | 8     | q=6 (exact 48.0)        | q=5 (+5.0)     |

At n_max=4, modes split into two quality tiers:
- Low-q modes (0-3, assigned q=3,4): residuals 0.8-2.8 (reasonable finite-size)
- High-q modes (4-7, assigned q=5,6): residuals 3.0-5.0 (poor assignment)

### Triangle Inequality Test

| n_max | Mode | q_assigned | Verdict | Violation rate |
|:-----:|:----:|:----------:|:-------:|:--------------:|
| 3     | 0    | q=2        | PASS    | 0%             |
| 3     | 1    | q=3        | **FAIL**| 30%            |
| 4     | 0-3  | q=3,4      | PASS    | 0%             |
| 4     | 4-5  | q=5        | **FAIL**| 66%            |
| 4     | 6    | q=5        | **FAIL**| 63%            |
| 4     | 7    | q=6        | **FAIL**| 86%            |

**Pattern:** modes PASS the triangle inequality if and only if their assigned q
is less than or equal to the maximum shell separation available in the graph
(max Δn = 1 on Fock graph). Since the graph only has edges with Δn ≤ 1, any
mode with q > 2 can couple electron states (n, n+1) that should be forbidden
by |n₁-n₂| ≤ q ≤ n₁+n₂-2 in the continuum.

### Ward Identity Diagnostic

| n_max | ||[Σ_T, H₀]|| / ||Σ_T|| |
|:-----:|:--------------------------:|
| 3     | 0.587                     |
| 4     | 0.573                     |

Both values ~58% — far from zero. The transverse self-energy does not commute
with the free Hamiltonian, meaning it mediates transitions between shells with
no Ward-identity constraint on the tensor structure.

### Eigenvector Support Analysis

Key finding: co-exact eigenvectors at n_max=4 split into two structural classes:

**Class A (modes 0-3):** localized on upper shell-pairs (3,3), (3,4), (4,4) with
ZERO weight on (1,2), (2,2), (2,3). These modes "don't see" the lower shells.

**Class B (modes 4-7):** delocalized across ALL shell-pairs including (2,2) and
(2,3). These are the modes that violate the triangle inequality.

The mechanism is clear: Class A modes live far from the graph boundary where the
finite-size effects of nearest-neighbor connectivity are minimal. Class B modes
span the full graph and are distorted by the boundary.

## Structural Obstruction

The obstruction is **topological, not spectral**:

1. **Eigenvalues converge** to n(n+2) for lower modes (correct spectrum)
2. **Eigenvectors are confined** to adjacent shell-pairs (wrong support)

In the continuum, a photon of angular momentum q is a 1-form on S³ that couples
ANY two electron shells satisfying the triangle inequality. On the Fock graph,
all edges connect shells with |Δn| ≤ 1, so the eigenvectors of d₁d₁ᵀ are
linear combinations of these Δn=±1 and Δn=0 edges. Even if the eigenvalue is
correct, the eigenvector CANNOT have the correct support because it's built from
edges that are too short-range.

This is the photon analog of the pendant-edge theorem in reverse: the graph's
nearest-neighbor connectivity PROTECTS the GS (pendant edge excluded from
plaquettes → GS zero) but PREVENTS the remaining selection rules (edges too
short-range → triangle inequality violated).

## Conclusion: Refined Selection Rule Partition (Final)

| Count | Category | Mechanism | Status |
|:-----:|:---------|:----------|:------:|
| 1     | Always survives | Gaunt/CG angular sparsity | GRAPH-INTRINSIC |
| 3     | Spinor-recoverable | Δm_j, parity, Furry (Dirac graph nodes) | GRAPH-INTRINSIC |
| 1     | Plaquette-recoverable | GS structural zero (co-exact topology) | GRAPH-INTRINSIC |
| 3     | Vector-quantum-number-required | SO(4) channel count, Ward, charge conj. | CALIBRATION |

The 5/8 ratio (graph-intrinsic / total) is maximal without promoting the photon
from a topological 1-cochain to a continuum vector harmonic carrying (L, M_L).

**What would be needed to recover the final 3:** The photon propagator must
become a matrix in an internal angular-momentum space (not just a scalar on edges),
with entries indexed by (q, m_q). This is equivalent to replacing the graph's edge
Laplacian with a fiber bundle where each edge carries a vector space of photon
modes — essentially the continuum QED structure projected onto the graph skeleton.

## Paper 18 Taxonomic Placement (Sharpened)

- **Tier 1 (graph-intrinsic):** Gaunt sparsity + spinor quantum numbers + co-exact
  topology → 5 selection rules from pure graph structure
- **Tier 2 (calibration):** photon (L, M_L) quantum numbers within the co-exact
  sector → remaining 3 selection rules. This is the vector structure calibration:
  promoting 1-cochains to sections of a vector bundle.
- **Tier 3 (embedding):** continuum spectral density matching (C_VP, C_SE, C_F₂)

The Tier 1 / Tier 2 boundary is now precisely located: it sits at the transition
from topological (co-exact eigenvalues and pendant exclusion) to geometric
(angular momentum labels on the fiber). The graph gives you the correct
number of modes and approximately correct energies, but not the angular structure
within each mode.

## Key Files

- `debug/transverse_mode_analysis.py` — diagnostic script
- `debug/data/transverse_mode_analysis.json` — numerical results
- `debug/transverse_qed_self_energy.py` — prerequisite (GS zero recovery)
- `debug/transverse_qed_self_energy_memo.md` — prerequisite memo
