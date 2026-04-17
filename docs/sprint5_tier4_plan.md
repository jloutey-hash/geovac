# Sprint 5: Drake Varshalovich + S⁵ Gauge Structure + Li/Be Core Polarization

**Version target:** v2.16.0
**Date:** April 2026
**Tracks:** Three parallel tracks, all independent

---

## Track DV: Drake Direct/Exchange Mixing via Varshalovich Bipolar Harmonics

**Goal:** Close Sprint 4 DD's honest negative. Derive the direct/exchange mixing ratios (3/50, −2/5, 3/2, −1) from the full bipolar harmonic expansion of 1/r₁₂³ (Varshalovich §5.17, Brink-Satchler App. 5). The J-pattern is already closed via Racah 6j; DV completes the spatial-sector derivation.

**Principle:** Algebraic Deconstruction + Transcendental Cataloging

### Background

Sprint 4 DD derived the Drake J-coefficients f_SS = (-2, +1, -1/5) and f_SOO = (+2, +1, -1) from Racah 6j symbolically. The J-pattern has a clean structural origin: both normalize to -6 via the identity (-1)^(L+S+J)·6j{L,S,J;S,L,k}.

**The open gap:** The spatial coefficients 3/50, -2/5, 3/2, -1 in front of M²_dir, M²_exch, M¹_dir, M¹_exch were NOT derived from 9j alone. DD's direct-Slater-determinant computation with universal r<^K/r>^(K+3) kernels gave the correct factorization structure f(J)·(c_d·M_d + c_e·M_e) with constant c_d/c_e = -3√5/5 across J, but the J-pattern disagreed with f_SS — showing that the r<^K/r>^(K+3) kernel is only one of several (k_1, k_2)-channel-specific kernels in the bipolar harmonic expansion of Y^(K)(r̂₁₂)/r₁₂^(K+1).

The missing machinery: the bipolar harmonic expansion decomposes Y^K(r̂₁₂)/r₁₂^(K+1) into a sum over (k_1, k_2, K) triples with specific (k_1, k_2)-dependent radial kernels. Varshalovich §5.17 gives the exact formulas; Brink-Satchler App. 5 gives the coupling coefficients. The Drake 1971 combining coefficients should emerge as specific projections onto the ³P two-electron channels.

### Sub-tracks

**DV-A: Derive the bipolar harmonic expansion**

1. Read `debug/dd_drake_derivation.md` and `debug/dd_drake_direct_sd.py` for DD's direct SD computation and where it stopped.
2. Starting from Varshalovich eq. (5.17.?) or Brink-Satchler eq. (App. 5.?):
   $$\frac{Y^K_M(\hat r_{12})}{r_{12}^{K+1}} = \sum_{k_1, k_2} A^{k_1 k_2 K}(r_1, r_2) \cdot [Y^{k_1}(\hat r_1) \otimes Y^{k_2}(\hat r_2)]^K_M$$
   where the radial kernel A^{k_1 k_2 K}(r_1, r_2) has explicit closed form involving r<, r>, and the (k_1, k_2, K)-triangle.
3. For 1/r₁₂³ Breit-Pauli (K=1 for SOO, K=2 for SS), enumerate the allowed (k_1, k_2) triples.

**DV-B: Project onto 1s-2p ³P angular channels**

1. For He 2³P_J with configuration 1s(l=0) · 2p(l=1), the angular wavefunction is a specific LS-coupled product.
2. Apply the bipolar expansion to the Breit-Pauli matrix element ⟨1s, 2p; ³P_J | 1/r₁₂³ × σ·σ | 1s, 2p; ³P_J⟩.
3. The allowed (k_1, k_2) for 1s·2p → 1s·2p direct are constrained by l_a=l_c=0 (k_1=0) and l_b=l_d=1 (k_2=0 or 2); for exchange l_a=0, l_c=1 gives k_1=1, and l_b=1, l_d=0 gives k_2=1.
4. Each channel contributes a specific (radial integral) × (angular coefficient via 9j or 6j).
5. Sum over channels: the coefficients 3/50, -2/5, 3/2, -1 should emerge as the projection result.

**DV-C: Verify and document**

1. Confirm the Sprint 4 BF-D rational-search values match the derived expressions exactly (sympy symbolic equality).
2. Generalize: what's the formula for arbitrary LS multiplet (e.g., Li 2²P which has (1s, 1s, 2p) and needs a different channel count)?
3. Update `debug/dd_drake_derivation.md` with the closed derivation.
4. Update `geovac/breit_integrals.py` docstring to reference the Varshalovich machinery.
5. Update Paper 14 §V paragraph: change the current "J-pattern derived from Racah 6j, direct/exchange mixing ratios found by rational search" to "both J-pattern and mixing ratios derived from Varshalovich bipolar harmonic expansion."
6. Add test in `tests/test_breit_integrals.py` that asserts the full symbolic equivalence.

### Success criteria
- Bipolar harmonic expansion of 1/r₁₂³ written symbolically
- ³P direct/exchange mixing ratios (3/50, -2/5, 3/2, -1) derived symbolically (exact sympy)
- Match with BF-D rational search confirmed to machine precision
- Test added + Paper 14 §V updated
- (Stretch) Generalization formula for other LS multiplets (Li 2²P, Be 2s2p ³P)

### Files to read first
- `debug/dd_drake_derivation.md` (Sprint 4 DD memo; partial derivation)
- `debug/dd_drake_direct_sd.py` (direct SD computation that hit the dead end)
- `debug/dd_drake_derivation.py` (sympy 6j/9j scaffolding)
- `geovac/breit_integrals.py` (production module — do NOT modify except docstring)
- Varshalovich "Quantum Theory of Angular Momentum" §5.17 (if accessible)
- Brink & Satchler "Angular Momentum" Appendix 5 (if accessible)

### Honest-negative contingency

If Varshalovich/Brink-Satchler are inaccessible and reconstructing the bipolar expansion from first principles is taking too long, the minimum acceptable outcome is a clear characterization of where the derivation stops. DO NOT tune coefficients to match BF-D — derive from first principles, then verify. A clean partial result is better than a false positive.

---

## Track S5: S⁵ Bargmann-Segal Gauge Structure

**Goal:** Extend Paper 25's framework observation (Hopf graph as lattice gauge theory) to the S⁵ Bargmann-Segal HO lattice (Paper 24). Does the HO sector have an analogous abelian U(1)-on-edges decomposition, or does SU(3) force a non-abelian gauge reading?

**Principle:** Natural Geometry Search + Transcendental Cataloging

### Background

Paper 25 (Sprint 4 QG) formalized: GeoVac's S³ Hopf graph at n_max is simultaneously a discretization of S³ (Paper 7), a simplicial 1-complex with discrete Hodge decomposition (Lim 2020), and a Wilson-type lattice gauge structure with matter on nodes and U(1) on edges via ladder-operator phases. Under Paper 2's conjectural interpretation, the S² base = photon, S¹ fiber = gauge.

Paper 24 (Bargmann-Segal) established the 3D HO discretization on the Hardy-space sector of S⁵, with nodes (N, l, m_l) and edges generated by SU(3) (N, 0) symmetric rep transitions (ΔN=±1, Δl=±1). The HO rigidity theorem (Paper 24 Theorem 3) says: the 3D isotropic HO is the unique central potential whose spectrum arises from the Euler operator on the Hardy space H²(S⁵) restricted to (N, 0) SU(3) irreps.

**The open question** (Paper 25 §VII, subsection "S⁵ Bargmann-Segal analog"): Does the S⁵ HO lattice have an analogous gauge structure? Three possibilities:

1. **Abelian U(1) on edges.** The SU(3) (N, 0) reps have no nontrivial U(1) subgroup beyond the trivial phase. If edges carry only U(1) phases, we get a lattice QED analog but with different spectrum.

2. **Non-abelian SU(3) on edges.** The SU(3) transitions naturally carry SU(3) matrix elements. This would be a lattice SU(3) gauge theory (like QCD), but with matter in a very specific representation.

3. **No gauge structure.** The S⁵ HO lattice has a different structural role — the Hopf fibration S¹ → S³ → S² is not the natural fibration for S⁵. S⁵ has the complex Hopf fibration S¹ → S⁵ → CP² instead, which has a different gauge structure.

Option 3 is the most likely honest negative; it would sharpen Paper 24's structural asymmetry finding (Coulomb second-order with calibration π vs HO first-order π-free).

### Sub-tracks

**S5-A: Literature + Paper 24 re-reading**

1. Read `papers/core/paper_24_bargmann_segal.tex` carefully — especially §IV (Coulomb/HO asymmetry), §V (HO rigidity theorem).
2. Literature: complex Hopf fibration S¹ → S⁵ → CP² (Hopf 1931, Nash-Sen), its role in CP² gauge theory.
3. SU(3) (N,0) Young-diagram reps: how do they decompose under the U(1) Cartan and under the CP² geometry?
4. Brief WebSearch: "lattice gauge theory on S⁵", "CP² Hopf", "SU(3) symmetric rep edge Laplacian"

**S5-B: Construct the S⁵ HO graph explicitly at small N_max**

1. Following Paper 24's π-free certificate machinery: nodes (N, l, m_l) with N ∈ {0, 1, ..., N_max}, l ∈ {N, N-2, ..., 1 or 0}, m_l ∈ {-l, ..., +l}.
2. Edges from dipole transitions ΔN=±1, Δl=±1, Δm_l ∈ {-1, 0, +1}. Edge weights are exact rational SU(3) Clebsch-Gordan-like matrix elements.
3. At N_max=5 (Paper 24 certified bit-exact π-free): 56 nodes, 165 edges.

**S5-C: Edge Laplacian + complex Hopf quotient**

1. Compute the signed incidence matrix B (56 × 165 at N_max=5).
2. Compute L_edge = B^T B (165 × 165). First Betti number = 165 - 56 + 1 = 110 (assuming connected).
3. Compute the complex Hopf quotient: collapse m_l-fibers at fixed (N, l) to get the CP² base graph. Nodes: Σ_{N≤5} (number of l values for each N) = 6 + 6 + 6 + ... quite a lot.
4. Compare to Paper 25's S² quotient: does the CP² quotient have analogous spectral structure?

**S5-D: Test for gauge structure**

1. Does the edge Laplacian of the S⁵ HO graph contain U(1)-like phases from SU(3) CG coefficients?
2. Or does it carry SU(3) matrix elements structurally?
3. Check: under α → 0, does any gauge-like structure reduce to classical HO dynamics? (α is the fine structure constant; it doesn't appear in the HO, so the question is whether ANY coupling constant plays an analogous role.)

**S5-E: Structural conclusion**

Three possible outcomes:
1. **Positive:** S⁵ HO has an abelian or non-abelian lattice gauge structure analogous to S³ Hopf. Paper 25 gets a §VII.1 subsection.
2. **Clean negative:** The complex Hopf fibration CP² base is structurally different; no direct gauge analog. Paper 24's HO rigidity theorem gains an extension about the absence of gauge structure.
3. **Mixed:** Some aspects carry over (edge Laplacian as discrete Hodge-1), but the gauge interpretation doesn't (no natural coupling constant playing α's role).

### Success criteria

- S⁵ HO graph constructed at N_max=5 with exact-rational edge weights (Paper 24 certification pattern)
- Edge Laplacian computed with spectrum
- CP² quotient graph computed and analyzed
- Clear structural verdict: gauge / no gauge / mixed

### Failed approaches to check

- Not much overlap with existing dead ends (this is new territory)
- Paper 24's HO rigidity theorem bounds what's possible — any proposed gauge structure must respect first-order operator + (N,0) SU(3) rep restriction

### Files to read first

- `papers/core/paper_24_bargmann_segal.tex` (HO discretization + rigidity theorem)
- `papers/synthesis/paper_25_hopf_gauge_structure.tex` (Paper 25 §VII.1 open question)
- `debug/qg_paper25_memo.md` (Sprint 4 Paper 25 context)

### Honest-negative contingency

If the structural analysis shows the HO sector has NO gauge analog, that's a valuable result — it sharpens Paper 24's Coulomb/HO asymmetry thesis by extending it from transcendental content (π vs no π) to structural content (gauge vs no gauge). Document as such.

---

## Track CP: Li/Be Fine Structure via Core Polarization

**Goal:** Close Sprint 3 BF-E honest negatives. Extend Breit-Pauli fine-structure benchmarks from He 2³P (Sprint 3: 0.20% error) to Li 2²P (currently 211% with leading Z_eff α²) and Be 2s2p ³P (currently 78%). Target: <20% error on both.

**Principle:** Algebraic Deconstruction (the core-polarization operator is a known closed-form addition)

### Background

Sprint 3 BF-E found:
- **Li 2²P:** NOT met with physical Z_eff. Best fit Z_eff=0.534 gives 0.5% but is unphysically small. Standard Z_eff=1 (Z_eff⁴ convention) gives +118%, matching the T8 baseline behavior. Root cause: a single valence electron above a closed 1s² core needs explicit core-polarization — the 1s² core responds to the 2p electron's electric field, screening the effective 2p potential. No SS contribution at S=1/2.
- **Be 2s2p ³P:** NOT met. Best result (88.4% max error) has correct sign structure but magnitudes off by 50-90%. Root cause: near-degenerate 2s/2p shell requires multi-configuration treatment — the 2³P state mixes with other ³L states via configuration interaction, which pure perturbation theory on 1s²·2s2p doesn't capture.

### Sub-tracks

**CP-A (Li): Core polarization potential**

Standard atomic physics treatment for an alkali-like single-valence atom (Migdalek & Bylicki 1998, Johnson Ch. 9):

$$V_{\rm cp}(r) = -\frac{\alpha_d}{2 r^4} W(r)$$

where α_d is the dipole polarizability of the closed core (for Li: α_d ≈ 0.192 a.u.) and W(r) is a cutoff function (standard form: W(r) = 1 - exp(-(r/r_c)^6) for a cutoff radius r_c ≈ 0.55 a.u.).

1. Add V_cp as a one-body operator on the 2p valence orbital.
2. Compute its expectation value ⟨V_cp⟩_{2p} symbolically (Laguerre × r^4 moment + cutoff integral).
3. Add to the leading Z_eff·α² fine-structure correction.
4. Target: Li 2²P doublet splitting within 20% of NIST 0.3354 cm⁻¹.

**CP-B (Be): Multi-configuration treatment**

The Be 2s2p ³P mixes with nearby ³D, ³S states. Standard treatment (Fischer multi-configuration):

1. Identify the mixing configurations: 2s2p ³P mixes primarily with 2p²¹D (through two-body ERI) and 2p3p ³P (small; higher shell).
2. Build a small CI matrix in the 2-config basis: {2s2p ³P, 2p² ¹D → ³P sector via projection}.
3. Diagonalize → get the ³P multiplet energies.
4. Apply Breit-Pauli SS + SOO via the Sprint 3 machinery, but evaluated on the mixed wavefunction.
5. Target: Be 2s2p ³P multiplet span within 20% of NIST 3.15 cm⁻¹.

**CP-C: Generalize and document**

1. If CP-A and CP-B meet targets: update Paper 14 §V fine-structure table with working Li and Be rows.
2. Document the extensions: core polarization for alkalis, multi-config for open-shell ns·np configurations.
3. If the targets are met for Li but not Be (or vice versa), document honestly.

### Success criteria

- Li 2²P doublet splitting within 20% of NIST
- Be 2s2p ³P multiplet span within 20% of NIST
- Paper 14 §V updated with working rows (or honest negative documentation)
- Tests added for the new closed-form operators

### Files to read first

- `debug/bf_d_benchmark_memo.md` (Sprint 3 BF-D/E baseline)
- `debug/bf_d_verify.py` (BF-E Li and Be attempts)
- `geovac/breit_integrals.py` (Sprint 3 production module — use unchanged)
- `geovac/spin_orbit.py` (T2 single-particle pattern)
- `papers/core/paper_14_qubit_encoding.tex` §V (fine-structure table)

### Failed approaches to check

- Don't try to tune Z_eff to match NIST (Sprint 3 BF-E shows this is unphysical)
- Don't replicate BR-C's minimal Bethe-Salpeter §39 scoping (Sprint 2 showed it's under-scoped for fine structure)

### Honest-negative contingency

If core polarization alone doesn't close Li to <20% (possible — Migdalek's treatment typically gets 10-20% but can be worse with simple cutoff functions), document the limit and identify the next layer (e.g., quantum-electrodynamic QED Lamb shift corrections at ~1% level). Don't claim the target was met with unphysical parameter tuning.

---

## Sprint 5 PM Prompt

```
Read CLAUDE.md, docs/sprint5_tier4_plan.md, and:
- docs/sprint4_final_summary.md (context for DV)
- papers/synthesis/paper_25_hopf_gauge_structure.tex (context for S5)
- debug/bf_d_benchmark_memo.md (context for CP)

Dispatch three parallel sub-agents:

Track DV (Drake Varshalovich):
  Sub-agent 1: Derive the bipolar harmonic expansion of 1/r₁₂³ symbolically
  and project onto He 2³P channels to reproduce (3/50, -2/5, 3/2, -1) from
  first principles. Starting point: debug/dd_drake_direct_sd.py (Sprint 4
  dead-end). Add test + Paper 14 §V update.

Track S5 (S⁵ gauge structure):
  Sub-agent 2: Construct the S⁵ HO graph at N_max=5 (56 nodes, 165 edges)
  per Paper 24 certification. Compute edge Laplacian + complex Hopf (S¹→S⁵→CP²)
  quotient. Determine whether the HO sector has an abelian/non-abelian gauge
  structure analog to Paper 25's S³ case, or a structural gap.

Track CP (Core polarization):
  Sub-agent 3: Implement core polarization V_cp for Li 2²P (Migdalek-Bylicki)
  and multi-config for Be 2s2p ³P. Target <20% error vs NIST on both.
  Integrate with Sprint 3's geovac/breit_integrals.py machinery.

Algebraic-first: DV uses exact sympy bipolar expansion; S5 uses exact
rational SU(3) CG coefficients; CP adds closed-form V_cp (known form).

Exit criteria:
- DV: mixing ratios derived symbolically, or clean honest-negative on exact
  machinery needed
- S5: structural verdict on S⁵ gauge structure (positive/negative/mixed)
- CP: Li and Be <20% (or honest-negative per atom)
```

---

## Sprint 6 Candidates (preview)

After Sprint 5:
- **Paper 18 taxonomy consolidation** (Sprint 4 carried over): distributional embedding + log-embedding + Racah-intrinsic vs tensor-dep housekeeping.
- **[Rn] frozen core + RaH matched comparison** (Sprint 3 open): engineering extension.
- **Paper 25 §VII open questions** depending on S5 outcome.
- **Tier-3 T7 Kramers-Pasternak** (long-deferred): closed-form Dirac-Coulomb ⟨r^{-2}⟩, ⟨r^{-3}⟩ for n_r ≥ 1 states.

---

## What Sprint 5 Resolves

1. **Drake fine-structure fully derived** (DV): J-pattern + mixing ratios both from Racah/Varshalovich algebra, no rational search. Strengthens He 2³P claim to fully first-principles.

2. **S⁵ gauge structure verdict** (S5): Extends Paper 25's framework observation to the HO sector, either confirming an abelian/non-abelian analog or sharpening Paper 24's Coulomb/HO asymmetry.

3. **Li and Be fine structure** (CP): Closes Sprint 3's honest negatives. Paper 14 §V fine-structure table gains working rows.
