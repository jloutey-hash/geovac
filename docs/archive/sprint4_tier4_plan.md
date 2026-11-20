# Sprint 4: Drake Derivation + T3 Regression Fix + QED-on-S³ Framing

**Version target:** v2.15.0
**Date:** April 2026
**Tracks:** Three parallel tracks, all independent

---

## Track DD: First-Principles Drake Coefficient Derivation

**Goal:** Derive the Drake 1971 combining coefficients A_SS = α²(3/50·M²_dir − 2/5·M²_exch), A_SOO = α²(3/2·M¹_dir − 1·M¹_exch) from first principles using Breit-Pauli tensor operators + ³P angular coupling. Confirm structurally that BF-D's rational-search result matches a Racah 9j identity.

**Principle:** Algebraic Deconstruction + Transcendental Cataloging

### Background

Sprint 3 Track BF-D found the Drake coefficients by brute-force rational search over the production `geovac/breit_integrals.py` output. The coefficients have small denominators (2, 5, 10, 25, 50), strongly suggesting a clean Racah algebra origin. BF-D's `debug/bf_d_racah_derivation.py` was incomplete — the sub-agent made a first-principles attempt using sympy's Wigner 9j but did not close the derivation.

The Breit-Pauli two-body operator for the ³P manifold has a standard tensor decomposition:
- H_SS: rank-2 spin tensor coupled to rank-2 radial tensor via a rank-0 combining structure (SS is a scalar operator in J-space, so rank-0 total)
- H_SOO: rank-1 × rank-1 coupled similarly

For a specific LS-coupled ³P state |(1s, 2p); L=1, S=1; J⟩, the matrix element of H_SS or H_SOO is:

⟨³P_J | H | ³P_J⟩ = Σ_k (angular factor via 9j) × (radial M^k or N^k integral)

The angular factor involves a 9j symbol coupling the two-electron spin and orbital spaces to the ³P total angular momentum. The 3/50, -2/5, etc. coefficients should emerge as products of 3j, 6j, 9j symbols for this coupling.

### Sub-tracks

**DD-A: Compute the 9j angular factors symbolically**

1. Read `debug/bf_d_racah_derivation.py` (BF-D's incomplete first-principles attempt)
2. For the He 2³P states (L=1, S=1, J ∈ {0,1,2}), compute the 9j angular factors:
   - SS: ⟨(l_a=0, l_b=1; L=1) (s=1/2, s=1/2; S=1) J | T^{(2)}_q(σ_1⊗σ_2) · T^{(2)}_{-q}(r̂) | same⟩
   - SOO: analogous with rank-1 coupling
3. Use sympy.physics.wigner for 3j, 6j, 9j. Keep everything symbolic.

**DD-B: Derive the combining coefficients**

1. For each k ∈ {0, 1, 2}, compute the angular factor × radial factor in LS basis.
2. Project onto the direct vs exchange radial integrals using Slater-rule decomposition.
3. The Drake coefficients should appear as rational combinations of 3j/6j/9j values.
4. Verify: the derived A_SS and A_SOO match BF-D's rational-search result (3/50, -2/5, 3/2, -1).

**DD-C: Generalize and document**

1. If the derivation succeeds, document the generalization: what's the formula for arbitrary LS multiplet? (E.g., ³S, ¹D, etc.)
2. Add a test in `tests/test_breit_integrals.py` verifying the Drake coefficients symbolically.
3. Update `geovac/breit_integrals.py` docstrings with the Racah identity references.
4. Update Paper 14 §V paragraph with the derivation reference (currently says "found by rational search").

### Success criteria
- Drake coefficients derived symbolically from Wigner 9j (sympy exact)
- Match to BF-D's rational-search result verified
- Test in test_breit_integrals.py that asserts the identity
- If successful: generalization formula for other LS multiplets

### Failed approaches to check
- BF-D's `bf_d_racah_derivation.py` attempt (incomplete) — understand why it was incomplete before re-attempting

### Files to read first
- `debug/bf_d_racah_derivation.py` (incomplete first-principles attempt)
- `debug/bf_d_coef_search.py` (the brute-force search that found 3/50, -2/5, 3/2, -1)
- `debug/bf_d_benchmark_memo.md` (full BF-D context)
- `geovac/breit_integrals.py` (the production module — don't modify during derivation, only update docstrings at the end)
- Bethe-Salpeter §39 or Johnson "Atomic Structure Theory" Ch. 8 (if accessible)

---

## Track TR: Tier-2 T3 Regression Fix

**Goal:** Fix the jj-vs-LS angular coupling discrepancy in the Tier-2 T3 relativistic builder. DC-B (Sprint 2) found that at α=0, the spinor FCI does not reproduce scalar FCI exactly — gap 0.95 → 1.44 → 1.66 mHa growing with n_max. Must diagonalize angular couplings via Clebsch-Gordan jj↔LS transformation.

**Principle:** Algebraic Deconstruction (the jj↔LS transformation is exact via CG coefficients)

### Background

DC-B §5 observed: at n_max=1, both spinor and scalar FCI give E=-13.500 Ha exactly (single 1s² config). At n_max≥2, the gap emerges and grows. Both bases share the same radial R^k integrals via `hypergeometric_slater`, so the discrepancy is in the angular part.

Suspected causes (from DC-B memo):
1. Missing phase or normalization factor in the √((2j+1)(2j'+1)) prefactor of the spinor X_k coefficients
2. jj → LS projector implicitly applied in the finite basis but not reproducing the LS singlet cleanly
3. Selection rule mismatch between the scalar c^k and spinor X_k Gaunt coefficients

The problem is purely algebraic: the two bases should be related by a unitary transformation (Clebsch-Gordan coefficients coupling (l, s) to j). In a complete basis, scalar and spinor FCI must give identical energies at α=0. The ~1 mHa gap says the transformation is incomplete or has a bug.

### Sub-tracks

**TR-A: Isolate the bug at n_max=2**

1. Read DC-B's memo and script (`debug/dc_b_convergence_memo.md`, `debug/dc_b_dirac_cusp_convergence.py`).
2. Read `geovac/composed_qubit_relativistic.py` — find the X_k (spinor Gaunt) construction, the reduced 3j prefactor, and the jj basis enumeration.
3. At n_max=2, He (Z=2), α=0:
   - Diagonalize the scalar FCI matrix → get the exact eigenvalue
   - Diagonalize the spinor FCI matrix → get the 1 mHa-higher eigenvalue
   - Project the scalar ground state onto the spinor basis via explicit CG coefficients
   - Compare matrix elements between the two FCI Hamiltonians in the projected basis
   - Identify where they differ (which angular channel has the wrong weight)

4. Produce a minimum reproducing test case showing the single matrix element that's off.

**TR-B: Fix the bug**

1. Based on TR-A's diagnosis, locate the root cause in `composed_qubit_relativistic.py`.
2. Fix it — likely a missing normalization factor, a wrong phase, or an incomplete sum.
3. Verify the fix doesn't break existing tests.
4. Re-run DC-B at Z=4, n_max=2,3,4 to confirm the spinor(α=0) == scalar FCI to machine precision.

**TR-C: Add regression test**

1. Write `test_relativistic_alpha0_matches_scalar` in `tests/test_spin_ful_composed.py`:
   - At Z=4 (Be²⁺), n_max=2, α=0: assert |E_spinor − E_scalar| < 1e-10 Ha
   - Repeat for n_max=3 (if time permits — n_max=3 is Q=28, may be slow)
2. Add to the existing Tier 2 test suite.

### Success criteria
- Bug located and fixed
- spinor(α=0) FCI matches scalar FCI to < 1e-10 Ha across n_max=2,3
- Regression test added
- No existing test broken (including 251 Tier 1-3 tests + Sprint 3's 116 new tests)

### Files to read first
- `debug/dc_b_convergence_memo.md` §5 "Secondary finding" (the regression observation)
- `debug/dc_b_dirac_cusp_convergence.py` (Step 2 shows the gap)
- `geovac/composed_qubit_relativistic.py` (where the bug lives)
- `geovac/dirac_matrix_elements.py` (Szmytkowski angular + X_k construction)
- `tests/test_spin_ful_composed.py` (existing regression tests)

---

## Track QG: QED-on-S³ Framing Paper

**Goal:** Draft a positioning paper documenting the framework-level observation that GeoVac's Hopf decomposition IS a concrete instance of graph Hodge theory applied to lattice gauge theory — a connection nobody in either community has articulated. This is positioning prose, not a computation, but it's still delegatable given the rich context in Sprint 1's outputs.

**Principle:** Transcendental Cataloging (the paper catalogs where gauge-theoretic structure enters) + Natural Geometry Search (the Hopf bundle as the natural gauge geometry)

### Background

Sprint 1 Track Q1 surfaced three independent observations that, together, form a paper-worthy framing:

1. **No "Fock projection for photons" exists** (massless dispersion → light cone, not sphere). But the photon has a natural S³ home via conformal compactification (Maxwell conformally invariant in 4D).

2. **Paper 2 (conjectural) interprets the Hopf bundle as the electron-photon geometry**: total space S³ = electron manifold, base S² = photon momentum sphere, fiber S¹ = U(1) gauge phase. α is "the cost of projecting S³ → S²".

3. **Graph Hodge theory** (Lim 2020, Schaub et al. 2020) defines the edge Laplacian L₁ = B^T B as the discrete analog of the Hodge Laplacian on 1-forms. Nobody has connected this to lattice gauge theory, despite the obvious correspondence (matter on nodes = 0-forms, gauge fields on edges = 1-forms).

4. **GeoVac's Hopf decomposition is a concrete instance**: the S² quotient graph (Phase 4B α-D) IS the discrete photon propagation manifold if Paper 2's interpretation holds.

### What's new (and what's NOT)

**New:** The synthesis. Nobody has written "GeoVac's Hopf graph is a lattice gauge theory where matter lives on S³ nodes and gauge fields on S¹ fibers."

**Not new:** Any specific component (Hopf bundle, lattice gauge theory, conformal Maxwell, graph Hodge theory are all well-established separately).

**Conjectural:** Paper 2's physical interpretation of the Hopf components. The framing paper must be careful to distinguish:
- Mathematical observation (the graph decomposition exists and has a Hodge-theoretic structure)
- Physical interpretation (that S² base = photon, S¹ fiber = gauge) — Paper 2 conjecture, not proven

### Sub-tracks

**QG-A: Literature review (depth)**

Starting from Sprint 1 Q1's survey (`debug/q1a_literature_survey.md`):
1. Deepen the lattice-gauge-on-sphere search: Luscher compact-manifold lattice gauge theory, Witten Chern-Simons on S³, modern work (post-2015) on simplicial complexes and higher-form gauge theories.
2. Deepen the graph Hodge theory search: recent (2020+) work on simplicial gauge theory, topological data analysis applied to physics.
3. Identify the specific gap: why hasn't the connection been made? (Usually: the lattice-gauge community uses hypercubic flat lattices; the Hodge-theory community isn't thinking about physics; the GeoVac-like approach is too niche.)

**QG-B: The structural observation, formalized**

Write up the key claim:
- The GeoVac S³ graph at n_max=m has V_m = Σ_{n=1}^m n² nodes labeled by (n,l,m_l) and E_m edges given by the angular/radial transition rules.
- The Hopf quotient under U(1) (collapsing m_l-fibers at fixed (n,l)) gives an (n,l)-labeled quotient graph with V_m^S² = Σ_{n=1}^m n = m(m+1)/2 nodes.
- The edge Laplacian L₁ = B^T B of the full graph is the discrete Hodge-1 Laplacian.
- The Hopf quotient + the node/edge decomposition correspond to a lattice gauge structure: matter on nodes, gauge field on edges (fibers), photon propagation on the S² base.
- Within Paper 2's conjectural interpretation, α is a specific combinatorial invariant of this structure.

**QG-C: Draft the paper**

Target: `papers/synthesis/paper_25_hopf_gauge_structure.tex` (or `papers/observations/` if that fits better — check existing structure).

Sections:
1. **Introduction:** The Hopf bundle as the electron-photon geometry (Paper 2 recap).
2. **Graph Hodge theory primer:** Node Laplacian (0-forms), edge Laplacian (1-forms), Betti numbers, exact/coexact decomposition.
3. **The GeoVac Hopf graph as a lattice gauge theory:** Specific mapping matter ↔ nodes, gauge ↔ edges, photon ↔ S² base.
4. **What this predicts/doesn't predict:** No new α derivation (Paper 2 stays conjectural); but the framework statement makes Paper 2's K = π(B + F - Δ) formula a specific combinatorial invariant.
5. **Related work:** Luscher, Witten CS, Lim 2020, Schaub 2020. Why the gap remained.
6. **Open questions:** Does the S⁵ Bargmann-Segal lattice (Paper 24) have an analogous gauge structure? Does the rank-1 (SOO) / rank-2 (SS) tensor decomposition (Paper 22 Sprint 2 extension) have a gauge-theoretic reading?

Length target: 15-20 pages (synthesis paper format).

Tone: observational, not assertive. The framework observation is sharp; its physical interpretation is conjectural. Frame accordingly.

### Success criteria
- Draft paper in `papers/synthesis/` or `papers/observations/` with clear structural content
- Literature review deeper than Sprint 1 Q1
- No over-claiming: the observation is the framework connection, not a new α derivation
- Cross-reference integration: Paper 2 (conjectural α), Paper 22 (angular sparsity), Paper 24 (HO rigidity), Papers 18 (exchange constants)

### Rhetoric constraints (CRITICAL)

Per CLAUDE.md §1.5 rhetoric rule:
- Do NOT assert physical priority of the graph description over continuous QM.
- Paper 2's interpretation is conjectural; frame all physical claims accordingly.
- The concrete advance is the framework observation (Hopf-as-lattice-gauge), not a physics breakthrough.

Per CLAUDE.md §1 project context:
- This is an AI-augmented independent research project. Target: GitHub + Zenodo DOI.
- Do NOT format for a specific journal or suggest traditional peer review.
- Frame as scientific observation worth documenting, not as a product claim.

### Files to read first
- `papers/conjectures/paper_2_alpha.tex` (especially §IV post-Tier-1 rewrite)
- `papers/core/paper_22_angular_sparsity.tex` (rank-2 Breit extension Sprint 2)
- `papers/core/paper_24_bargmann_segal.tex` (HO rigidity, π-free graph)
- `papers/core/paper_18_exchange_constants.tex` (taxonomy — where does the gauge structure fit?)
- `papers/synthesis/paper_21_geometric_vacuum_synthesis.tex` (synthesis format template)
- `debug/q1a_literature_survey.md` (Sprint 1 QED-on-S³ survey)
- `debug/q1_reframing_memo.md` (Sprint 1 reframing after PI correction)
- `debug/data/track_alpha_phase4b/track_d_graph_morphism.json` (S² quotient graph data)

---

## Sprint 4 PM Prompt

```
Read CLAUDE.md, docs/sprint4_tier4_plan.md, and:
- docs/sprint2_final_summary.md (context for TR)
- docs/sprint3_final_summary.md (context for DD)
- debug/q1a_literature_survey.md + debug/q1_reframing_memo.md (context for QG)

Dispatch three parallel sub-agents:

Track DD (Drake Derivation):
  Sub-agent 1: First-principles derivation of the Drake coefficients
  (3/50, -2/5, 3/2, -1) via Wigner 9j symbolic computation. Verify match
  with BF-D's rational search. Add test + Paper 14 §V update. Starting
  point: debug/bf_d_racah_derivation.py (incomplete).

Track TR (Tier-2 T3 Regression Fix):
  Sub-agent 2: Locate and fix the jj↔LS coupling bug in
  geovac/composed_qubit_relativistic.py that causes spinor FCI at α=0
  to differ from scalar FCI by 0.95-1.66 mHa. Add regression test. No
  existing tests should break.

Track QG (QED-on-S³ Framing Paper):
  Sub-agent 3: Draft papers/synthesis/paper_25_hopf_gauge_structure.tex
  (or similar path) documenting the framework observation that GeoVac's
  Hopf decomposition IS a concrete instance of graph Hodge theory applied
  to lattice gauge theory. Observational framing, NOT a new α derivation.
  Follow CLAUDE.md §1.5 rhetoric rule.

Algebraic-first: DD uses exact sympy; TR is a real bug (no numerical
papering-over); QG is prose but should reference exact math.

Exit criteria:
- DD: coefficients derived symbolically, test passes, Paper 14 §V updated
  (or honest-negative memo if derivation doesn't close)
- TR: regression test passes at < 1e-10 Ha agreement, no existing tests broken
- QG: draft paper in papers/synthesis/, properly cross-referenced, rhetoric-
  consistent with CLAUDE.md §1.5
```

---

## What Sprint 4 Resolves

1. **Drake coefficients structural origin** (DD): from rational-search to Racah identity. Strengthens Sprint 3's He 2³P result.

2. **Tier-2 T3 structural correctness** (TR): spinor FCI at α=0 must equal scalar FCI. Closing this loose end fixes a ~1 mHa discrepancy that affects all Tier 2+ resource claims at n_max≥2.

3. **Framework observation formalized** (QG): Paper 25 (or wherever it lands) documents the graph-Hodge-to-lattice-gauge connection. The observation is framework-level, not a specific physics result — it sharpens the interpretation of existing work (Paper 2's α conjecture, Paper 22's angular sparsity, Paper 24's HO rigidity) by placing them in a gauge-theoretic context.

## Sprint 5 Candidates (preview)

Surfaced by Sprints 1-4 but not addressed:
- Li/Be fine structure via core polarization (Tier 3+ extension)
- [Rn] frozen core for RaH direct Sunaga comparison
- First-principles generalization of Drake coefficients to arbitrary LS multiplets (if DD succeeds)
- Operator-order × bundle grid extensions (Paper 18 taxonomy completeness)
