# Sprint 5 Track S5: S^5 Bargmann-Segal gauge-structure analysis

**Date:** 2026-04-15
**Status:** Complete
**Verdict:** MIXED — abelian U(1) transfers verbatim; SU(3) does not; the m_l-quotient is not a CP^2 discretization. Sharpens Paper 24's Coulomb/HO asymmetry.

## Task

Sprint 4 Track QG produced `papers/synthesis/paper_25_hopf_gauge_structure.tex` — the framework observation that GeoVac's S^3 Hopf graph at finite n_max is simultaneously a Fock-projection discretization of S^3, a triangle-free simplicial 1-complex with discrete Hodge decomposition, and a Wilson-type U(1) lattice gauge structure. Paper 25 §VII.1 lists an open question: does the S^5 Bargmann-Segal lattice (Paper 24) have an analogous gauge structure?

Three candidate outcomes:
- (a) Abelian U(1) analog — simple, same mechanism as S^3
- (b) Non-abelian SU(3) analog — the S^5 → CP^2 fibration is complex, so SU(3) may be natural
- (c) No gauge structure — sharpens the Paper 24 asymmetry

## Summary of findings

| Component | S^3 Coulomb (Paper 25) | S^5 Bargmann (this sprint) |
|:---|:---|:---|
| Nodes | (n, l, m_l) | (N, l, m_l) |
| Edges | L± angular, T± radial | dipole ΔN=±1, Δl=±1 |
| π-free | yes | yes |
| Orientation | L± direction | creation-operator direction (lower N → higher N) |
| Incidence matrix B, L_0=BB^T, L_1=B^T B | well-defined | well-defined |
| U(1) gauge action covariant under ψ_v → e^(iχ_v) ψ_v | yes | **yes** |
| β_1 (cycles) | 2 at n_max=3 | **110 at N_max=5** |
| Node Laplacian L_0 computes physical spectrum | yes, via κ = −1/16 | **no — HO spectrum is in the diagonal** |
| SU(3)-covariant lattice structure | N/A | **no (see below)** |
| m_l-quotient matches base-space Laplacian | yes (6-sector S^2 base, spectrum {0,0,0,1,3,6}) | **no (12 sectors; spectrum does not match CP^2 Fubini-Study λ_k = 4k(k+2))** |

## Part 1-2: Graph construction

Used `geovac/nuclear/bargmann_graph.build_bargmann_graph(5)` directly (the production module implements Paper 24 with `fractions.Fraction` throughout).

Reproduced Paper 24 certificate:
- Nodes: **56** (target 56) ✓
- Edges: **165** (target 165) ✓
- All diagonal entries rational (Fraction) ✓
- All adjacency entries rational (Fraction) ✓
- `pi_free: True` ✓

Node count by shell N: [1, 3, 6, 10, 15, 21] = Σ (N+1)(N+2)/2, ✓.
Edge count by (N → N+1) pair: {0→1: 3, 1→2: 12, 2→3: 27, 3→4: 48, 4→5: 75}.
Edge class breakdown by (l → l'): selects l → l±1 (dipole), with 3j + radial squared factors giving 57 distinct rational weight values.

See `debug/data/s5_bargmann_graph.json` for full breakdown.

## Part 3: Edge Laplacian, β_1, Hodge structure

Built the signed incidence matrix B ∈ Z^{56 × 165} with the creation-operator orientation (edge from lower-N endpoint to higher-N endpoint). Computed in exact integer / sympy arithmetic.

- **L_0 = B B^T** (56 × 56), node Laplacian: has one zero eigenvalue, i.e. β_0 = 1 (connected graph).
- **L_1 = B^T B** (165 × 165), edge Laplacian: has **β_1 = 110 zero eigenvalues**.
- **SVD theorem check**: nonzero eigenvalues of L_0 and L_1 agree exactly (max diff = 0.0, sympy exact arithmetic).
- **Euler characteristic**: V − E = 56 − 165 = −109; β_0 − β_1 = 1 − 110 = −109 ✓ (consistent).

Interpretation: β_1 = 110 means the S^5 Bargmann lattice at N_max = 5 has 110 independent cycle classes, each carrying a potentially non-trivial Wilson-loop holonomy in the abelian U(1) gauge interpretation.

## Part 4: Complex Hopf quotient (m_l fiber collapse)

Collapsing m_l fibers at fixed (N, l) gives **12 sectors** at N_max = 5:
```
[(0,0), (1,1), (2,0), (2,2), (3,1), (3,3), (4,0), (4,2), (4,4), (5,1), (5,3), (5,5)]
```
with fiber dimensions (2l+1) = [1, 3, 1, 5, 3, 7, 1, 5, 9, 3, 7, 11].

**This is NOT a faithful discretization of CP^2.**

The true CP^2 base of the complex Hopf fibration S^1 → S^5 → CP^2 is a 4-real-dimensional Kähler manifold with scalar Laplacian eigenvalues (Fubini-Study, holomorphic sectional curvature K = 4):
```
λ_k = 4 k (k+2),   k = 0, 1, 2, ...
deg λ_k = (k+1)^3
```
giving the sequence (0, 12, 32, 60, 96, 140, ...).

The quotient graph Laplacian (multiplicity-weighted) has the 12 eigenvalues:
```
{0, 2.22, 4.87, 4.94, 11.65, 12.11, 18.74, 28.09, 33.58, 50.61, 63.77, 99.43}
```

Rescaling comparison to CP^2:
- Best least-squares linear scale: quotient ≈ 0.1391 × CP^2, but max relative residual 24.98%
- Ratios vary from 0.08 to 0.19 (non-uniform)

Under the total-weight and fiber-average schemes the spectrum is also irregular. None of the three tested weight schemes matches the CP^2 spectrum to any consistent rescaling.

**Verdict (c) is confirmed for the spectral question:** the m_l quotient is NOT a discretization of CP^2 in the Fubini-Study Laplacian sense. It is a combinatorial SO(3)-labeled graph on 12 sectors, which has too few nodes to faithfully represent the 4-dimensional CP^2 Laplace-Beltrami spectrum. This is consistent with the structural feature that **the Bargmann graph's edges carry SO(3)-adapted (l, m_l) labels, not full CP^2 coordinates**.

## Part 5: Gauge-structure tests

### (a) Abelian U(1) test — POSITIVE

The incidence matrix B gives a well-defined exterior derivative d_0 = B^T from 0-cochains (node wavefunctions) to 1-cochains (edge phases). A node-local gauge transformation ψ_v → e^(iχ_v) ψ_v induces

 U_{vw} → e^(-iχ_v) U_{vw} e^(iχ_w)

on every edge amplitude U_{vw} = ⟨w | z_q | v⟩ (dipole matrix element). This is the standard Wilson U(1) gauge transformation on a simplicial 1-complex. The fact that the UNGAUGED amplitudes U_{vw} happen to be real positive (Condon-Shortley) does not prevent the U(1) bundle from existing — it just means the "default" trivialization has no holonomy. β_1 = 110 cycles carry the available holonomy degrees of freedom.

**Conclusion:** the abelian U(1) structure of Paper 25 transfers verbatim to the S^5 Bargmann graph. The combinatorial Hodge vocabulary (L_0, L_1, Hodge decomposition, Wilson loops) applies without modification.

### (b) Non-abelian SU(3) test — NOT NATURAL

The SU(3) (N, 0) irreps transform into each other under the dipole operator z_q (which is the fundamental 3 of SU(3)). But the Bargmann graph is built in the **SO(3)-adapted (l, m_l) basis** within each (N, 0) irrep, so the graph edges carry SO(3)-covariant labels, not SU(3)-covariant labels.

For an SU(3) lattice gauge structure, one would need:
1. **SU(3)-covariant node labels**: instead of (l, m_l), use the full SU(3) weight-space labels of the (N, 0) irrep.
2. **SU(3) link variables**: instead of U(1) phases, link variables would be elements of SU(3).
3. **A natural gauge action**: node-local SU(3) rotations would act on the link variables by left/right multiplication.

The obstruction: **transitions between (N, 0) and (N+1, 0) are intertwiners between DIFFERENT irreps, not SU(3) group elements**. The intertwiner
 V_{(N, 0)} ⊗ V_{(1, 0)} → V_{(N+1, 0)}
projects onto the (N+1, 0) irrep in the tensor product. It is a Clebsch-Gordan coupling, not a Wilson link variable.

A proper Wilson SU(3) lattice gauge theory has a single SU(3) irrep on every site (matter) and SU(3) group elements on every link (gauge). The Bargmann graph has **different SU(3) irreps on different shells**, which is structurally a representation tower (an SU(3) Hilbert space) rather than a lattice gauge theory.

**Conclusion:** the natural gauge group of the Bargmann graph as built is U(1), not SU(3). An SU(3)-covariant reformulation of the HO on S^5 is a different construction (possibly interesting, but outside the Paper 25 Wilson-Hodge frame).

### (c) CP^2 base spectrum test — NEGATIVE

Already discussed in Part 4.

## Part 6: Structural verdict

**MIXED: (a) positive + (c) negative on the spectral side.**

The S^5 Bargmann graph:

**DOES** carry the same abelian U(1) Wilson-Hodge structure as the S^3 Coulomb graph:
- Signed incidence B, node Laplacian L_0 = B B^T, edge Laplacian L_1 = B^T B all well-defined.
- β_1 = 110 (vs β_1 = 2 for S^3 at n_max = 3).
- Node-local U(1) gauge transformations act covariantly on edge amplitudes.
- Discrete Hodge decomposition R^{165} = im(d_0) ⊕ ker(L_1) applies.

**DOES NOT** carry the Paper 25 physical content that makes the S^3 gauge structure nontrivial:
- On S^3, the node Laplacian **computes the physical spectrum** via κ = −1/16: (D − A) ψ = (n² − 1) ψ, and kappa · (n² − 1) = Rydberg. **On S^5, the Bargmann node Laplacian is spectrally inert**: the HO spectrum ℏω(N + 3/2) lives in the DIAGONAL, not in D − A. The graph adjacency encodes only dipole transitions (spectroscopic content), not spectral content. This was already stated by Paper 24 §II.D ("the graph adjacency encodes dipole transitions only, NOT the spectrum"), and the gauge analysis confirms it at the level of the abelian U(1) interpretation: the Wilson dictionary applies combinatorially (L_1 = B^T B), but the node-matter propagator content is reduced — on S^3, L_0 propagates matter and generates physical energies; on S^5, L_0 only diffuses on the dipole graph, the physical Hamiltonian is elsewhere (the diagonal).

**DOES NOT** lift to a Wilson SU(3) structure on the (N, 0) symmetric irreps:
- Transitions between distinct (N, 0) irreps are intertwiners, not group elements. The Wilson construction requires a fixed group fiber on every link.

**DOES NOT** give a spectrally faithful CP^2 discretization via m_l-quotient:
- 12 sectors is too few to represent the 4-dimensional CP^2 Laplace-Beltrami spectrum.
- Relative residual after best linear rescaling is 25% (no uniform rescaling reconciles the two sequences).

### Sharpening Paper 24's Coulomb/HO asymmetry

Paper 24 established the asymmetry at the level of **spectrum**: on S^3, the graph Laplacian (D − A) computes the Coulomb spectrum via κ (nonlinear projection, carries calibration π); on S^5, the graph Laplacian is only spectroscopic (first-order Euler operator on the diagonal, linear-affine projection, zero calibration π). The sprint adds a **gauge-structure layer** to the same asymmetry:

> The abelian U(1) Wilson-Hodge structure is UNIVERSAL (it applies to any finite simplicial 1-complex with a chosen orientation). On both S^3 and S^5, the incidence matrix B, the node Laplacian L_0, and the edge Laplacian L_1 are well-defined and obey the discrete Hodge theorem. However, the **physical content** of these objects differs categorically:
>
> **On S^3 Coulomb:** L_0 is a spectrum-computing operator (Laplace-Beltrami after scaling by κ). The U(1) gauge interpretation makes the pair (matter propagator, photon propagator) = (L_0, L_1) a nontrivial double, and Paper 25's framework observation places Paper 2's (B, F, Δ) ingredients in distinct corners of this double.
>
> **On S^5 HO:** L_0 is not a spectrum-computing operator — the HO spectrum is in the diagonal, not in D − A. The U(1) gauge structure is still there combinatorially, but it has **reduced physical content**: the edge Laplacian L_1 generates cycle classes and carries β_1 = 110 Wilson-loop degrees of freedom, but there is no "matter propagator" role for L_0 to fulfill. The HO is already diagonal; nothing needs to be propagated.

This is a new dimension of the Coulomb/HO asymmetry: the **physical interpretability** of the universal Wilson-Hodge structure depends on whether the graph Laplacian carries spectrum-computing content. It does on S^3, it does not on S^5. Calibration π (Paper 24) and spectrum-computing Laplacian (this sprint) go together in the Coulomb case; neither appears in the HO case.

An additional consequence: the Paper 2 analog "K = π(B + F − Δ) as a spectral invariant of the Hopf gauge structure" cannot transfer to S^5 in any natural form, because:
1. There is no `B` analog — the SU(3) → SO(3) Casimir trace on the (N, 0) shells is not naturally weighted by a gauge-bundle measure, because there is no `S^1 → S^5 → CP^2` Wilson-like framework that links the Casimir to the gauge structure.
2. There is no `F` analog — the Fock-degeneracy Dirichlet D_{g_n}(s) = D_{(N+1)(N+2)/2}(s) at s = packing-exponent d_max = 6 (dim C^3 = 6) gives
   D_{g_n}(6) = (1/2)[ζ(4) + 3 ζ(5) + 2 ζ(6)],
   which is a mixed even/odd zeta combination — no clean ζ_R(2) = π²/6 analog emerges.
3. There is no `Δ` analog — the single-level Dirac mode count on S^5 would be computed by the Camporesi-Higuchi spectrum on S^5, but Paper 23 established Fock rigidity is unique to S^3 (not S^5), so the Dirac mode count does not naturally enter an α-like combination on S^5.

All three points are consistent with Paper 24 §V (Coulomb/HO asymmetry) and Paper 23 §III (Fock projection rigidity).

## Part 7: Paper updates

Two targeted updates are appropriate:

1. **Paper 24** (core): add a short subsection after §V.C (HO two-fermion entanglement corollary) titled "Gauge-structure corollary to the HO rigidity theorem", documenting that the abelian U(1) Wilson-Hodge structure transfers combinatorially but the spectrum-computing role of the node Laplacian does not. This sharpens the Coulomb/HO asymmetry of §IV into a two-layer statement: spectrum-computing AND gauge-structure-physical-content are both Coulomb-specific; the universal content is the combinatorial Wilson-Hodge vocabulary, not the physical Wilson interpretation.

2. **Paper 25** (synthesis): update §VII.1 ("Higher Hopf bundles — open question") to reference this result: the abelian U(1) structure transfers, SU(3) does not, CP^2 spectral identification fails, Coulomb/HO asymmetry is sharpened rather than bridged. Frame as "question ANSWERED, with MIXED verdict: abelian U(1) extension is trivial (universal); SU(3) non-abelian extension is not natural; K-analog on S^5 does not exist."

See `debug/data/s5_graph_spectrum.json` and `debug/data/s5_bargmann_graph.json` for full numerical data.

## Files

- `debug/s5_bargmann_segal_graph.py` — graph construction at N_max=5, reproduces Paper 24 certificate
- `debug/s5_edge_laplacian_analysis.py` — incidence matrix, Laplacians, quotient, gauge tests
- `debug/data/s5_bargmann_graph.json` — basic graph stats (56, 165, π-free, β_1 = 110)
- `debug/data/s5_graph_spectrum.json` — full spectra, quotient analyses, verdict
- `debug/s5_gauge_structure_memo.md` — this memo

## Runtime

Graph build: < 1s (pure fractions.Fraction arithmetic).
Node-Laplacian sympy eigenvalues (56×56): ~15 s.
Edge-Laplacian sympy eigenvalues (165×165): ~1 minute.
Quotient sympy eigenvalues (12×12, three schemes): ~3 s.
Total run: ~2 minutes.
