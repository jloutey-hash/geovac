# VQ-3: Sigma-weighted Fock graph vertex -- Memo

**Date:** 2026-05-01
**Track:** VQ-3 (Vector QED sprint)
**Status:** COMPLETE -- NEGATIVE

## Headline

Inserting Pauli sigma matrices at the scalar CG vertex (V_mu_e = sigma_mu @ V_e) does NOT recover any of the 4 missing continuum QED selection rules. The sigma-weighted construction produces exactly 3x the scalar self-energy trace and 3x the scalar F_2, consistent with the trace identity sum_mu sigma_mu^2 = 3I. The selection rule count remains 1/8 (Gaunt/CG sparsity only), the same as scalar QED. Two previously surviving rules are broken (Delta_m_j conservation, Furry's theorem). This confirms that the missing selection rules require promoting the photon from a scalar 1-cochain to a vector harmonic, not dressing the vertex with spin matrices.

## Construction

The sigma-weighted vertex combines the inter-shell connectivity of the Fock graph (CG projection vertex V_e from `geovac/graph_qed_vertex.py`) with the spin structure of the Pauli matrices:

    V_mu[a,b,e] = sum_c sigma_mu[a,c] * V_scalar[c,b,e]
                = (sigma_mu @ V_e)[a,b]

for mu in {x, y, z} and each Fock graph edge e. The sigma matrices are built in the coupled (j, m_j) DiracLabel basis via CG coefficients, preserving n and l while coupling different kappa values (j changes while l stays fixed).

Self-energy uses Hermitian contraction (V^dagger instead of V^T) since sigma_mu are Hermitian:

    Sigma_vec = sum_mu sum_{e,e'} G_gamma[e,e'] * (sigma_mu @ V_e) @ (sigma_mu @ V_{e'})^dagger

This fixes VQ-2's disconnected-graph problem: VQ-2 used sigma alone (intra-shell only), giving a graph with no inter-shell edges. VQ-3 uses sigma @ V_scalar, preserving the inter-shell connectivity.

## Quantitative Results (n_max=2)

### Self-energy

| Quantity | Scalar | Sigma-weighted |
|:---------|:-------|:---------------|
| Tr(Sigma) | 44/3 = 14.667 | 44.000 |
| Zero eigenvalues | 5 | 0 |
| Positive semidefinite | Yes | Yes |
| GS block | [[1,1],[1,1]] | [[3,-1],[-1,3]] |
| GS structural zero | No | No |
| Eigenvalues | {0(x5), 4/3, 2(x2), 4, 16/3} | {4/3, 2(x2), 8/3, 4(x3), 16/3, 8, 32/3} |

**Trace ratio = 3 exactly.** This is the Pauli trace identity: sum_mu Tr(sigma_mu^2) = 3 * Tr(I).

### Vertex correction and F_2

| Quantity | Scalar | Sigma-weighted |
|:---------|:-------|:---------------|
| F_2 | 5*sqrt(2)/3 = 2.357 | 5*sqrt(2) = 7.071 |
| F_2 number field | Q(sqrt(2)) | Q(sqrt(2)) |
| F_2 ratio | 1 | 3.000 |

**F_2 ratio = 3 exactly.** F_2_sigma = 3 * F_2_scalar. Number field unchanged (both in Q(sqrt(2))).

### Selection rule census

| Rule | Scalar (GN-5) | VQ-3 sigma-weighted |
|:-----|:--------------|:--------------------|
| 1. Delta m_j conservation | PASS (|Dm_j| <= 1) | **FAIL** (|Dm_j| up to 2) |
| 2. Spatial parity E1 | FAIL | FAIL |
| 3. Gaunt/CG sparsity | PASS | PASS |
| 4. Vertex parity | FAIL | FAIL |
| 5. SO(4) channel count | FAIL | FAIL |
| 6. Charge conjugation | FAIL | **FAIL** (was FAIL before too) |
| 7. Furry's theorem | PASS (scalar) | **FAIL** (sigma breaks it) |
| 8. Ward identity | FAIL | FAIL |
| **Total** | **1/8** | **1/8** |

The sigma dressing preserves Gaunt/CG sparsity (1/8) but BREAKS two rules that the scalar vertex satisfied:

- **Delta_m_j**: sigma introduces kappa-changing transitions (Delta_kappa = +/-3) that couple different m_j values with |Delta_m_j| = 2 (beyond the continuum |Delta_m_j| <= 1).
- **Furry's theorem**: the sigma-weighted tadpole sum_e (sigma_mu @ V_e) is nonzero because sigma_z is diagonal (not off-diagonal like the identity on the scalar graph). The scalar tadpole sum_e V_e vanishes by the graph's T+/T- cancellation; sigma_z breaks this cancellation for within-shell edges.

### Vertex transition structure

The sigma-weighted vertex has (Delta_n, Delta_l, Delta_kappa) values:

    (-1, 0, 0), (0, 0, -3), (0, 0, 0), (0, 0, 3), (1, 0, 0)

The scalar vertex has only (Delta_n, Delta_l, Delta_kappa) = {(-1,0,0), (0,0,0), (1,0,0)}. The new Delta_kappa = +/-3 entries come from sigma coupling kappa=-1 (s_{1/2}) to kappa=-2 (p_{3/2}) via the CG transformation -- but this changes j, not l, so Delta_l = 0 throughout.

### Number field

- Self-energy entries are NOT rational (irrational entries present)
- F_2^2 = 50 exactly (rational), so F_2 = 5*sqrt(2) in Q(sqrt(2))
- Same number field as scalar QED

### Nonzero entry count

- Scalar vertex: 48 nonzero entries (across all edges)
- Sigma-weighted: 136 nonzero entries (2.83x inflation)

## Physics interpretation

### Why sigma dressing gives exactly 3x

The factor of 3 is the universal Pauli trace: sum_mu sigma_mu @ sigma_mu = 3I. Since the self-energy quadratic form is linear in the vertex-vertex product, and V^dagger = V^T for the scalar vertex part (real matrices), the sigma-weighted self-energy reduces to:

    Sigma_vec = sum_mu sigma_mu @ Sigma_scalar @ sigma_mu

and the trace identity gives Tr(Sigma_vec) = 3 * Tr(Sigma_scalar). The 3x factor carries through to F_2 because both numerator and denominator of the F_2 extraction formula scale uniformly.

### Why selection rules are not recovered

The sigma matrices preserve n and l but change kappa (and hence j). This means:

1. **E1 spatial parity (Delta_l = +/-1)** -- Not recovered because sigma preserves l. The Fock graph has Delta_l = 0, and sigma @ V still has Delta_l = 0. E1 requires the PHOTON to carry angular momentum, which means vector harmonics on the edge, not spin matrices on the node.

2. **Vertex parity (n1+n2+q odd)** -- Not recovered because sigma preserves n. The Fock graph's n structure (both same-shell and adjacent-shell edges) is unchanged.

3. **SO(4) channel count** -- Not recovered because the photon is still a scalar 1-cochain. W(n1,n2,q) = 0 requires the photon to have vector quantum numbers.

4. **Ward identity** -- Not recovered because [D, Sigma_vec] != 0 (the sigma dressing breaks the Dirac operator's diagonal structure in the kappa-mixed basis).

### Comparison to the three-tier selection rule partition

The selection rule census across all constructions:

| Tier | Rules | Scalar Fock | VQ-2 sigma-only | VQ-3 sigma@V | Dirac Rule B |
|:-----|:------|:-----------|:----------------|:-------------|:-------------|
| Always survives | Gaunt/CG | YES | YES | YES | YES |
| Spinor-recoverable | Dm_j, E1, Furry | 0/3 | 1/3 | -1/3 (breaks 2) | 3/3 |
| Vector-required | VP, SO4, Ward, C | 0/4 | 0/4 | 0/4 | 0/4 |

VQ-3 is WORSE than the scalar vertex for spinor-recoverable rules: it breaks Delta_m_j and Furry. The Dirac graph Rule B (which uses native Dirac nodes with E1 adjacency) recovers all 3 spinor-recoverable rules. This confirms the partition: spinor-recoverable rules require the GRAPH TOPOLOGY to encode parity-flipping (Delta_l = +/-1) transitions, not just spin matrices at the vertex.

## Implications

1. **Sigma vertex dressing is a dead end for selection rule recovery.** The 3x trace factor is a universal algebraic identity (Pauli trace), not new physics. All observables are trivially rescaled versions of the scalar QED results.

2. **The selection rule deficit is in the photon, not the electron vertex.** The 4 vector-required rules (vertex parity, SO(4) channel count, charge conjugation, Ward identity) all require the photon to carry vector quantum numbers. No amount of vertex modification on the electron side can recover them while the photon remains a scalar 1-cochain.

3. **The Dirac graph Rule B approach (GN-sprint, native_dirac_graph_memo.md) remains the correct intermediate construction.** It recovers 4/8 rules by using E1 adjacency (Delta_l = +/-1) as the graph topology, which naturally encodes spatial parity. The sigma-weighted vertex cannot replicate this because it preserves l.

4. **Next step for the VQ sprint:** the vector photon promotion (making the edge carry angular momentum quantum numbers) is the remaining construction to test. This is the only path to the 4 vector-required rules.

## Key files

- Script: `debug/vq3_sigma_weighted_vertex.py`
- Data: `debug/data/vq3_sigma_weighted_vertex.json`
- This memo: `debug/vq3_sigma_weighted_vertex_memo.md`

## Comparison to VQ-2

| Feature | VQ-2 (sigma-only) | VQ-3 (sigma @ V_scalar) |
|:--------|:-------------------|:------------------------|
| Graph connectivity | Intra-shell only (disconnected) | Full Fock graph (inter-shell T+/T-) |
| Photon propagator | N/A (no edges) | Same L_1^+ as scalar QED |
| Selection rules | 1/8 (Dm_j only) | 1/8 (Gaunt only) |
| F_2 | Undefined (no propagator) | 5*sqrt(2) = 3 * scalar |
| Verdict | NEGATIVE (disconnected) | NEGATIVE (trivial rescaling) |
