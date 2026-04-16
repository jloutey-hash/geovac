# Sprint 1: Native Dirac Graph + QED-on-S³ Scoping

**Version target:** v2.13.0
**Date:** April 2026
**Tracks:** Two parallel tracks (independent, no blocking dependencies)

---

## Track D6: Native Dirac Graph — Three-ℤ₂ Unification

**Goal:** Determine whether Paper 0's orientability ℤ₂, Szmytkowski's κ-parity ℤ₂, and Dirac chirality ℤ₂ are the same symmetry under Fock projection. If yes, build a native (n, κ, m_j) graph. If no, document which ℤ₂ breaks and why.

**Principle:** Algebraic Deconstruction + Transcendental Cataloging

### Sub-tracks

**D6-A: Enumerate the three ℤ₂s explicitly**

Write down the action of each symmetry on both labeling schemes:

| ℤ₂ | Name | Action on (n, l, m, s) | Action on (n, κ, m_j) | Source |
|-----|------|------------------------|------------------------|--------|
| σ | Paper 0 orientability | hemisphere parity on S² base | ? | Paper 0 §2 |
| P | Szmytkowski κ-parity | l ↔ l, s ↔ -s (spin flip) | κ → -κ | T1, dirac_matrix_elements.py |
| C | Dirac chirality γ⁵ | chirality ±1 sector swap | ? | D1, dirac_s3.py, Phase 4I D4 |

Deliverable: explicit matrix representation of each ℤ₂ acting on the 14-node n_max=3 graph. Identify fixed points, orbits, and quotient graphs.

**D6-B: Fock projection compatibility**

The scalar Fock projection maps S³ → ℝ³ via stereographic projection with p₀ = √(-2E_n). Check:
1. Does each ℤ₂ commute with the Fock projection? (i.e., does the symmetry exist on both sides of the conformal map?)
2. The Dirac-Coulomb accidental degeneracy E(n,κ) = E(n,-κ-1) (6 pairs confirmed in T8) — is this the orbit of one of the three ℤ₂s on the graph?
3. If the κ→-κ map is a graph automorphism, what are its eigenvalues? Do they reproduce any known spectral invariant?

Deliverable: for each ℤ₂, state whether it commutes with Fock projection (proof or counterexample). Classify the 6 Dirac degeneracy pairs by which ℤ₂ relates them.

**D6-C: Native Dirac graph construction (conditional on D6-A/B positive)**

If the three ℤ₂s unify (or at least two do):
1. Build a `DiracLattice` class with (n, κ, m_j) nodes and edges determined by the Dirac selection rules from T1 (Szmytkowski angular matrix elements)
2. Edge weights from T1's exact rational matrix elements
3. Verify: does the DiracLattice adjacency spectrum reproduce the Camporesi-Higuchi eigenvalues?
4. Verify: does α → 0 recover the scalar GeometricLattice spectrum?
5. π-free certificate: all eigenvalues and edge weights must be exact rationals (extending D1's certificate)

If negative: document which ℤ₂ is incompatible with Fock projection and why. This is a valuable structural result.

### Success criteria
- All three ℤ₂ actions written as explicit permutation matrices on the 14-node graph
- Fock compatibility proved or disproved for each
- Either: DiracLattice class with passing tests, OR: documented obstruction theorem

### Failed approaches to check
- None directly overlap. The Dirac sector is genuinely new territory.
- Track α-D (Phase 4B) built the n_max=3 graph and its S² quotient but did NOT examine graph automorphisms or ℤ₂ actions.

### Papers affected
- Paper 0: if positive, upgrade "structural correspondence" to theorem
- Paper 7: Dirac-Fock projection section
- Paper 18: ℤ₂ classification in exchange-constant taxonomy

### Files to read
- `geovac/lattice.py` (scalar graph construction)
- `geovac/dirac_s3.py` (Tier 1 Dirac infrastructure, SpinorHarmonicLabel)
- `geovac/dirac_matrix_elements.py` (DiracLabel, κ bridge, angular matrix elements)
- `geovac/hopf_bundle.py` (Hopf projection, fiber analysis)
- `debug/data/track_alpha_phase4b/track_d_graph_morphism.json` (14-node graph data)

---

## Track Q1: QED Photon Field on S³ — Edge Laplacian Scoping

**Goal:** Determine whether the photon field on S³ has a natural discretization on the GeoVac graph. Bounded to 1 sprint. Clean positive/negative outcome.

**Principle:** Natural Geometry Search

### Sub-tracks

**Q1-A: Literature survey — lattice QED on compact manifolds**

Targeted search for:
1. Lattice gauge theory on S³ specifically (Wilson, Lüscher, Luscher-Weisz)
2. Chern-Simons theory on S³ (Witten 1989) — this is a *topological* gauge theory, different from dynamical QED, but the S³ geometry is shared
3. Compact QED on curved backgrounds — does conformal coupling of the photon give a natural S³ formulation?
4. The photon's conformal weight in 4D: the Maxwell action is conformally invariant in 4D. Does this give a natural S³ projection analogous to Fock's matter projection?

Deliverable: 1-page summary of what's known. Key question: is there a "Fock projection for photons"?

**Q1-B: Edge Laplacian of the n_max=3 graph**

In standard lattice gauge theory, matter lives on nodes, gauge fields live on edges (links). The GeoVac graph at n_max=3 has 14 nodes and 13 edges.

Compute:
1. The incidence matrix B (14×13) of the n_max=3 graph
2. The edge Laplacian L_edge = B^T B (13×13)
3. The node Laplacian L_node = B B^T (14×14) — should match the existing graph Laplacian
4. Spectrum of L_edge: all 13 eigenvalues
5. Check: do any spectral invariants of L_edge (trace, determinant, zeta function, truncated Casimir) match K = π(B + F - Δ) ingredients?
   - B = 42? F = π²/6? Δ = 1/40? K/π = B + F - Δ?

This is a pure algebraic computation. The edge Laplacian is the combinatorial Hodge Laplacian on 1-forms — it's the natural "photon Laplacian" in graph gauge theory.

**Q1-C: Hodge decomposition and gauge structure (conditional on Q1-B finding anything)**

If Q1-B produces spectral hits:
1. Decompose the edge space into exact (im B^T), coexact (im B), and harmonic (ker L_edge) components
2. The harmonic component counts β₁ = first Betti number = loops in the graph. For a tree (13 edges, 14 nodes, 1 component): β₁ = 0. For a graph with cycles: β₁ > 0.
3. Gauge invariance in lattice QED corresponds to the exact component. How many gauge degrees of freedom?

If Q1-B produces no hits: document as clean negative. The photon field does not naturally live on this graph's edges.

### Success criteria
- Edge Laplacian spectrum computed exactly (rational arithmetic)
- Each K ingredient checked against edge spectral invariants
- Literature summary: is there a "Fock projection for photons"?
- Clean positive (spectral hit → open Q1-C) or clean negative (no hits → shelve)

### Failed approaches to check
- Phase 4B-4H exhausted sphere-spectral approaches for the α combination rule on the *node* Laplacian. The edge Laplacian is a genuinely different object — not in the dead-end list.
- SM-running (Phase 4H) is dead. This is NOT a perturbative QED calculation — it's a graph-combinatorial question.

### Papers affected
- Paper 2: if positive, the combination rule gets a graph-theoretic interpretation (matter on nodes, photon on edges)
- Paper 18: photon exchange constant classification

### Files to read
- `debug/data/track_alpha_phase4b/track_d_graph_morphism.json` (graph topology)
- `geovac/lattice.py` (adjacency matrix construction)
- `geovac/hopf_bundle.py` (Hopf fiber structure)

---

## Sprint 1 PM Prompt

```
Read CLAUDE.md, docs/sprint1_tier4_plan.md, and the files listed in each track.
Dispatch Track D6 and Track Q1 as parallel sub-agents.

Track D6 (Native Dirac Graph):
  Sub-agent 1 (D6-A+B): Enumerate all three ℤ₂ actions on the 14-node
  n_max=3 graph. Build explicit permutation matrices. Test Fock projection
  compatibility. Classify the 6 Dirac degeneracy pairs. Pure algebra —
  no new modules, just analysis + a debug/ script with JSON output.

  Sub-agent 2 (D6-C, conditional): If D6-A/B positive, build DiracLattice
  class. If negative, write obstruction theorem.

Track Q1 (QED Edge Laplacian):
  Sub-agent 3 (Q1-A): Literature survey. WebSearch for lattice QED on S³,
  conformal photon projection, Hodge Laplacian on Hopf graph.

  Sub-agent 4 (Q1-B): Compute edge Laplacian of the n_max=3 graph.
  Build incidence matrix from the adjacency data. Exact rational eigenvalues.
  Check all K ingredients against edge spectral invariants.

Gate: D6-C waits for D6-A/B results. Q1-C waits for Q1-B results.
Everything else runs in parallel.

Exit criteria:
- D6: three ℤ₂ permutation matrices + Fock compatibility verdict + either
  DiracLattice or obstruction theorem
- Q1: edge Laplacian spectrum + K-ingredient check + literature summary +
  positive/negative verdict
```

---

## Sprint 2-3 Dependencies

Sprint 2 (SS/SOO) depends on Sprint 1 only if D6 produces a native Dirac graph:
- If D6 positive: implement Breit interaction on DiracLattice
- If D6 negative: implement Breit on existing bolt-on architecture

Sprint 3 (heavy atoms) is independent. Can start anytime after Sprint 1.
