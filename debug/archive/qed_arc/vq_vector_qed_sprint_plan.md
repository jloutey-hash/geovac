# VQ Sprint: Vector QED on the Fock Graph

**Date:** 2026-05-01
**Motivation:** The graph-native QED construction (GN-1 through GN-7) gives scalar QED: the photon is a scalar 1-cochain on edges with no polarization. This recovers 1/8 continuum selection rules (scalar Fock) or 4/8 (Dirac Rule B graph). The 4 missing rules (vertex parity, SO(4) channel count, Ward identity, charge conjugation) all require vector photon structure.

**Key insight from conversation:** S³ = SU(2) has three intrinsic "directions" — the left-invariant vector fields X₁, X₂, X₃, which are the angular momentum operators L_x, L_y, L_z. The graph Laplacian uses their sum of squares (L²). To get vector structure, we need the individual L_i. The L-edges (Δm = ±1) carry directional information; the T-edges (Δn = ±1) are omnidirectional.

**Prior results:**
- VQ-1: σ·L decomposition → zero on all graph edges (intra-shell only). NEGATIVE.
- VQ-2: pure σ_μ vertex → disconnected graph (2/8 rules, worse than baseline). NEGATIVE.
- Lesson: need to combine graph adjacency (inter-shell) WITH directional structure (σ_μ).

## Track VQ-3: σ-weighted CG vertex (scalar photon + vector vertex)

**Goal:** Test whether inserting σ_μ at the CG vertex — keeping the same graph and same photon propagator — recovers any of the 4 missing selection rules.

**Construction:**
For each Fock graph edge e and each direction μ ∈ {x, y, z}:
```
V_μ[a,b,e] = Σ_c σ_μ[a,c] · V_scalar[c,b,e]
```
where V_scalar is the existing CG vertex from `graph_qed_vertex.py`.

The combined self-energy:
```
Σ_vec[a,b] = Σ_μ Σ_{e,e'} G_γ[e,e'] · (σ_μ · V_e)[a,:] · G_e · (σ_μ · V_{e'})^T[:,b]
```

This is "QED with vector vertex but scalar photon" — the photon still has no polarization, but the electron-photon coupling now includes the Pauli matrix.

**Deliverables:**
1. Three vertex matrices V_x, V_y, V_z (per edge) at n_max=2
2. Combined self-energy Σ_vec and its trace, eigenvalues, Hermiticity
3. Vertex correction Λ_vec and F₂ extraction
4. Full 8-rule selection rule census
5. Comparison table: scalar (1/8), Dirac Rule B (4/8), σ-weighted (?/8)
6. Number field characterization (π-free certificate)
7. Ground state self-energy Σ_vec(GS) — does structural zero recover?

**Key physics question:** Does σ_μ at the vertex restore vertex parity (n₁+n₂+q odd)?
The continuum vertex parity comes from γ^μ flipping parity. If σ_μ in the CG projection provides a similar parity flip, this rule might recover.

**Files to read:** `geovac/graph_qed_vertex.py`, `geovac/graph_qed_self_energy.py`, `debug/vq2_sigma_vertex_qed.py` (for σ_μ construction in Dirac basis)

## Track VQ-4: Direction-resolved edge Hodge decomposition

**Goal:** Classify Fock graph edges by their directional quantum numbers and build per-channel Hodge decompositions. This is the infrastructure track — no QED yet, just characterizing the directional structure of the graph.

**Construction:**
Classify edges by (Δn, Δm) into three channels:
- Channel (+): L_+ edges, Δm = +1 (angular, "east" on S²)
- Channel (−): L_- edges, Δm = -1 (angular, "west" on S²)
- Channel (0): T-edges, Δm = 0 (radial/z-direction, Δn = ±1)

For each channel:
1. Build the sub-incidence matrix B_μ (V × E_μ)
2. Build the sub-edge Laplacian L₁^{(μ)} = B_μ^T · B_μ
3. Compute spectrum, Betti numbers, connected components
4. Compute photon propagator G_γ^{(μ)} = (L₁^{(μ)})⁺

**Deliverables:**
1. Edge counts per channel at n_max = 2, 3
2. Channel-specific Betti numbers and spectra
3. Is each channel connected? (Important: if a channel is disconnected, its photon can't propagate globally)
4. Per-channel photon propagator matrices
5. Spectral comparison: is the Hodge decomposition "isotropic" (all channels identical) or "anisotropic"?

**Key physics question:** Does the z-channel (T-edges) have different topological properties than the ±-channels (L-edges)? If so, the photon propagation is inherently anisotropic, which would give the photon effective "polarization" even without explicit vector labels.

**Files to read:** `geovac/fock_graph_hodge.py` (for Hodge machinery), `geovac/graph_qed_photon.py`

## Track VQ-5: Full vector QED (direction channels × σ_μ vertex)

**Goal:** Build the complete vector QED construction combining direction-labeled photon channels (from VQ-4) with σ_μ vertex coupling (from VQ-3). This is the headline track.

**Construction:**
Three direction-specific QED channels:
```
Σ_vec[a,b] = Σ_μ Σ_{e∈E_μ, e'∈E_μ} G_γ^{(μ)}[e,e'] · (σ_μ · V_e)[a,:] · G_e · (σ_μ · V_{e'})^T[:,b]
```

Each direction μ has:
- Its own edge set E_μ (from VQ-4's classification)
- Its own photon propagator G_γ^{(μ)} (from VQ-4's per-channel Hodge)
- Its own vertex coupling σ_μ (from VQ-3's σ-modified CG vertex)

The total self-energy sums over all three channels.

**Deliverables:**
1. Per-channel self-energy Σ_x, Σ_y, Σ_z and combined Σ_vec = Σ_x + Σ_y + Σ_z
2. Per-channel vertex correction Λ_x, Λ_y, Λ_z and combined
3. Per-channel and combined F₂ extraction
4. Full 8-rule selection rule census on the combined vertex
5. Comparison table across ALL constructions:
   | Construction | Rules | Σ(GS) | F₂ | Number field |
   |-------------|-------|--------|-----|-------------|
   | Scalar Fock (GN-5) | 1/8 | ≠0 (pendant) | 5√2/3 | ℚ(√2,√3,√6) |
   | Dirac Rule B | 4/8 | ≠0 | — | ℚ(√2,√17,...) |
   | VQ-2 pure σ | 2/8 | ≠0 | — | ℚ(√2,...) |
   | VQ-3 σ-weighted | ?/8 | ? | ? | ? |
   | VQ-5 full vector | ?/8 | ? | ? | ? |
6. π-free certificate for the full construction
7. Isotropy check: is Σ_x + Σ_y + Σ_z isotropic (all three equal)?

**Key physics question:** Does the COMBINATION of direction-labeled channels AND σ_μ vertex recover any of the 4 missing rules? If vertex parity recovers, that's a significant result — it means the "vector structure" is accessible from the graph's intrinsic geometry without importing continuum physics.

**Dependency:** Requires VQ-4 edge classification and per-channel Hodge. Can reuse VQ-3's σ_μ vertex construction.

## Execution Plan

**Phase 1 (parallel):** VQ-3 and VQ-4 run simultaneously
**Phase 2 (sequential):** VQ-5 runs after VQ-3 and VQ-4 complete

## Exit Criteria

- **POSITIVE:** Any of the 4 missing rules (vertex parity, SO(4) channel count, Ward identity, charge conjugation) recovered → major result, write up for Paper 28
- **PARTIAL:** No new rules recovered, but structural insight gained (e.g., anisotropic Hodge, cleaner algebraic structure) → document as structural finding
- **NEGATIVE:** No new rules, no structural insight → confirms "graph = scalar, continuum = vector" partition definitively. Still valuable as a classification result.
