# Unified U(1) × SU(2) × SU(3) Gauge Content of GeoVac: Synthesis and Gap Analysis

**Date:** 2026-05-04
**Status:** Scoping memo (no paper edits)
**Verdict:** GeoVac now contains all three Standard-Model gauge groups as Wilson lattice constructions on three structurally distinct sub-manifolds, sharing the universal SU(N) kinetic coefficient 1/(4·N_c). The construction is a graph-spectral instance of the Marcolli–van Suijlekom 2014 gauge-network framework with the Perez-Sanchez 2024 correction. **It is not a Standard Model**: the three sectors live on three different sub-manifolds rather than a single spectral triple, the matter content is one species per construction rather than three generations, no Higgs-as-inner-fluctuation mechanism is in place, and hypercharge / electroweak unification structure is structurally absent. The cross-manifold gap is the most consequential.

---

## §1. What we have: the unified gauge content

Three Wilson lattice gauge constructions, built independently and confirmed at machine precision, exist as siblings on GeoVac sub-manifolds:

| Construction | Group | Host manifold | Graph | Source | Kinetic coefficient | Matter rep at node |
|:-------------|:-----:|:--------------|:------|:-------|:-------------------:|:-------------------|
| Paper 25 | U(1) | S³ (Hopf base) | Fock-projected S³ Coulomb graph (Paper 7) | `Paper 25` | 1/(4·1) = 1/4 (abelian → trivial) | scalar 0-cochain (electron states) |
| Paper 30 | SU(2) | S³ = SU(2) | Same S³ Coulomb graph | `Paper 30` | 1/8 = 1/(4·2) | spin-1/2 Dirac (Paper 14 Tier 2) |
| Sprint ST-SU3 | SU(3) | S⁵ (Bargmann-Segal Hardy sector) | (N, l, m_l) HO dipole graph (Paper 24) | `geovac/su3_wilson_s5.py` | 1/12 = 1/(4·3) | (N, 0) symmetric SU(3) irreps (CG-coupled) |

The three kinetic coefficients align on the **universal SU(N) Wilson formula**

> 1 / (4 · N_c) per plaquette per generator component,

derived from the Wilson action S_W = β · Σ_P (1 − (1/N_c)·Re Tr U_P) by expanding U_e = exp(i A_e^a T^a) with the standard tr(T^a T^b) = (1/2)δ^{ab} normalization. The N_c-dependence is structural (more colors share fixed (1/N_c) Re Tr normalization while the trace identity is held fixed). Each non-abelian construction has been verified independently against this prediction:

- Paper 30 §5.2: 1/8 (SU(2)) — symbolic + numerical, machine precision.
- ST-SU3 §4: 1/12 (SU(3)) — symbolic + numerical, rel. err ~ O(ε²).
- Both also verify that the eight (resp. three) generator components decouple at quadratic order; non-abelian self-interactions enter at O(A⁴) via [T^a, T^b] = i f^{abc} T^c.

The **maximal-torus reductions** are also clean siblings:

- Paper 30: SU(2) → U(1) reduces to Paper 25 exactly. Paper 25 is the Cartan-subalgebra projection of an SU(2) theory.
- ST-SU3: SU(3) → U(1) × U(1) reduces to a *coupled* abelian theory with character (1/3)(cos a_P + cos b_P + cos(a_P+b_P)). The two U(1) sectors do not decouple — the third phase −(a+b) is fixed by det = 1.

In each case, **the edge Laplacian L₁ = B^⊤B is the weak-coupling kinetic operator** on g-valued 1-cochains. Paper 25's "discrete Hodge-1 Laplacian / gauge propagator" is therefore the same combinatorial object across all three constructions; only its target Lie algebra changes. This is the single mathematical identity unifying the three sectors.

**Plaquette structure:** primitive closed non-backtracking walks (Ihara cycles, Paper 29). The S³ graph is sparse (β₁ = 0, 2, 8 at n_max = 2, 3, 4). The S⁵ Bargmann graph is plaquette-rich (β₁ = 6, 23, 56 at N_max = 2, 3, 4). The Wilson loop measurements show monotonic ⟨W⟩ in β with the expected SU(N)-color suppression (SU(3) lower than SU(2) at small β; SU(2) lower than U(1) by the same logic).

---

## §2. Lineage placement: Marcolli–van Suijlekom 2014, Perez-Sanchez 2024

The published precedent for placing all of (graph + matter on nodes + Wilson connection on edges + finite spectral triple) on a single object is Marcolli–van Suijlekom, "Gauge networks in noncommutative geometry," J. Geom. Phys. **75** (2014), arXiv:1301.3480. In that framework:

- Each vertex carries a finite spectral triple (A_v, H_v, D_v).
- Each oriented edge carries a connection in a chosen Lie group G.
- The composite spectral action of the network is Wilson lattice gauge theory in the leading order, with matter coupling determined by the choice of A_v and H_v.

Each of GeoVac's three Wilson constructions is a specific instance of this framework with a specific choice of vertex algebra, vertex Hilbert space, and gauge group:

| | A_v | H_v | G_e | Continuum content |
|:-:|:----|:-----|:-----|:------------------|
| Paper 25 | ℂ | ℂ (scalar) | U(1) | electromagnetism on S³ |
| Paper 30 | ℂ | ℂ² (Dirac) | SU(2) | non-abelian YM on S³ |
| ST-SU3 | ℂ | ⊕_N ℂ^{dim(N,0)} | SU(3) | non-abelian YM on S⁵ |

The **Perez-Sanchez 2024 correction** (arXiv:2401.03705 with the 2025 follow-up arXiv:2508.17338) clarifies that the continuum limit of Marcolli–van Suijlekom gauge networks is **Wilson Yang–Mills WITHOUT a Higgs field**. This is a structural correction to the original 2014 Higgs claim: the gauge-network spectral action recovers Yang–Mills, but the Higgs mechanism in Connes-Chamseddine spectral triples enters through a different ingredient (the inner fluctuation of D coupled to the off-diagonal block of an almost-commutative algebra), which gauge-network constructions do not automatically inherit.

GeoVac's three constructions sit precisely in the Perez-Sanchez-corrected lineage: each is gauge-network-style Wilson YM, none carries an automatic Higgs sector. This is consistent with the Sprint TS / WH1 reading that GeoVac is an almost-commutative spectral triple of finite cutoff but has not produced a continuum YM-Higgs sector by default.

The three GeoVac constructions therefore represent **three independent Wilson-YM instances of a single published mathematical framework**, distinguished by their host manifold and gauge group. This is a valid synthesis statement under CLAUDE.md §1.5: the Wilson-YM-on-graph machinery is published; what is novel is the assembly of *three* such instances on three GeoVac sub-manifolds with a unified kinetic coefficient.

---

## §3. Gap analysis: what GeoVac does NOT have for any SM-level claim

The honest distinction between "we have the three SM gauge groups" and "we have the Standard Model" lives in five structural features of the SM that GeoVac does not currently realize.

### 3.1 Three generations of matter

The Standard Model has three generations of fermions, related by structurally identical reps of the gauge group with different mass eigenvalues. GeoVac has one species per construction:

- Paper 30 has spin-1/2 Dirac matter (Paper 14 Tier 2) — one generation's worth.
- ST-SU3 has SU(3)-shell matter (one tower, no flavor-doubling).

Sprint 4H Track SM-B (April 2026) tested the naive shell-to-generation map (Σ_f N_c Q_f² = 8 = |λ_3| at three generations). The result was a **clean negative**: the per-generation contribution 8/3 is charge-universal across all three SM generations, while every per-shell S³ invariant varies with n. There is no shell→generation map consistent with charge-universal per-generation structure. The 8 = 8 equality was a numerical coincidence with no representation-theoretic content.

The actual structural obstruction is therefore at the level of the spectral triple: GeoVac's H_GV is a single-species Hilbert space (electronic for S³ Coulomb, nuclear-HO for S⁵ Bargmann). To carry three generations, H_GV would need to be replaced by H_GV ⊗ ℂ³ (or similar) with the generation index acting trivially under the gauge group but nontrivially under the mass operator. There is no GeoVac structure currently selecting this generation-tripling: the Fock projection (Paper 7), the SU(3) HO Bargmann projection (Paper 24), and the operator-system truncation (WH1 R2) all produce single-species Hilbert spaces.

This is a real and structural gap, not just a "we haven't computed it yet" gap.

### 3.2 Higgs as inner fluctuation of D

In Connes' Standard Model construction, the Higgs field appears as the non-commutative part of the inner fluctuation D ↦ D + [a, D] for a in an almost-commutative algebra A = C^∞(M) ⊗ A_F with A_F a finite matrix algebra (the Pati-Salam / ℂ ⊕ ℍ ⊕ M_3(ℂ) algebra). The Higgs is the off-diagonal component of the fluctuation between the two factors of A_F.

GeoVac's D_GV is the Camporesi-Higuchi Dirac on the Fock-projected S³ graph (Sprint TS, R3.5). Whether D_GV admits inner fluctuations depends on what algebra A is taken:

- If A_GV = C(graph nodes) (commutative diagonal multipliers, Paper 32 §III): inner fluctuations vanish — the Higgs sector is empty. This is the C*-envelope reading from WH1 R1.
- If A_GV is upgraded to the operator system O_{n_max} = P_{n_max} C^∞(S³) P_{n_max} (WH1 R2): there is more structure (off-diagonal Gaunt elements), but the operator system is *-closed but not a multiplicative algebra. Inner fluctuations of D depend on a multiplicative structure, so the natural Higgs construction does not directly apply. The ab ∉ O obstruction (witness pair M_{2,1,0}², 14.9% residual) is exactly the place a Higgs-style inner fluctuation would have to live; whether one can be defined in the operator-system setting is open in the NCG literature.
- For an almost-commutative extension A_GV ⊗ M_n(ℂ) (the Marcolli-vS gauge-network case): inner fluctuations do exist and produce gauge bosons — Paper 30 / ST-SU3 Wilson actions are the leading-order content. But the off-diagonal block structure that produces a Higgs in Connes' SM is NOT automatically in GeoVac's almost-commutative algebra. Per Perez-Sanchez 2024, Marcolli-vS gauge networks give YM-without-Higgs in the continuum.

The Higgs gap is therefore: GeoVac would need a specific almost-commutative structure with an off-diagonal-block algebra (analog of ℂ ⊕ ℍ ⊕ M_3(ℂ)) before inner fluctuations of D_GV could yield a Higgs field. None of the three Wilson constructions selects such a structure.

### 3.3 Hypercharge assignments and electroweak unification

The SM has the rep-theoretic relation Q = T_3 + Y_W/2 (Gell-Mann–Nishijima) tying the U(1)_em to T_3 of SU(2)_L and U(1)_Y. The three GeoVac Wilson constructions have:

- Paper 25's U(1) is on the S³ Hopf graph — interpretable as electromagnetic, but not as U(1)_Y.
- Paper 30's SU(2) is on the same S³ — could be SU(2)_L, but no chiral structure is built in (the Camporesi-Higuchi Dirac is non-chiral by default; chirality enters only through the offdiag CH variant of WH1 R3.5, and even there it's chirality-grading not weak-isospin chirality).
- The three constructions are on different manifolds — there is no co-location of U(1) and SU(2) needed for electroweak unification.

The hypercharge structure lives in the SM's representation content (left vs. right fermion fields with different Y_W). GeoVac's Hilbert space carries no left-right asymmetry as a built-in feature. Left-right asymmetry on the truncated operator system would have to come from an additional grading on H_GV, which would need to be argued for from the geometry rather than imposed.

### 3.4 The cross-manifold obstruction (most consequential gap)

Connes' Standard Model uses **a single spectral triple** (A_SM, H_SM, D_SM) on a single manifold (M = spacetime), with the gauge group emerging from inner automorphisms of an almost-commutative algebra A_SM = C^∞(M) ⊗ A_F. All three SM gauge groups (U(1)_Y, SU(2)_L, SU(3)_c) emerge from one A_F = ℂ ⊕ ℍ ⊕ M_3(ℂ) acting on one Hilbert space.

GeoVac has **three separate spectral triples on three separate sub-manifolds**:

- Paper 25 U(1) on S³ Coulomb graph.
- Paper 30 SU(2) on S³ = SU(2) Coulomb graph (same graph, different triple structure).
- ST-SU3 SU(3) on S⁵ Bargmann-Segal graph (different graph, different triple).

The S³ constructions can in principle be unified — both live on the same graph, same Hilbert space, same Dirac operator. Combining U(1) × SU(2) on the same S³ would produce something formally analogous to the electroweak sector. This is a clean target for a future Paper 25 × Paper 30 unification sprint.

The S⁵ SU(3) construction sits on a structurally different graph. There is no single GeoVac spectral triple containing both the S³ Coulomb graph and the S⁵ Bargmann graph as commuting factors; the four-layer Coulomb/HO asymmetry (Paper 24 §V; sharpened by ST-SU3 §5.4) explicitly states that S⁵ has no analog of S³'s spectrum-computing L_0, no calibration π, and no natural matter coupling for Wilson SU(3). The most one can say is that the Bargmann graph and the Coulomb graph could be encoded as orthogonal factors of a tensor-product Hilbert space H_GV = H_S³ ⊗ H_S⁵ — but there is no operator on this product space that geometrically unifies them.

This is the **single most consequential structural gap**. Connes' SM is a unification: one spectral triple, three gauge groups simultaneously. GeoVac is at present a **direct sum** of three Wilson-YM instances. The cross-manifold coupling that would make them sectors of one theory does not exist.

### 3.5 Other gaps worth naming

- **CKM matrix / fermion mixing**: requires off-diagonal mass matrices in the chiral basis. GeoVac has no chiral basis built in.
- **Anomaly cancellation**: SM hypercharge assignments satisfy Σ Y³ = 0, Σ Y = 0 per generation. GeoVac has no analog of this structural constraint.
- **Running couplings to GUT scale**: requires β-functions for each gauge factor on a common spectral cutoff. GeoVac's three Wilson actions have separate cutoffs (n_max for S³, N_max for S⁵), and there is no common renormalization-group flow.
- **Lorentzian signature**: all three Wilson constructions are Euclidean lattice gauge theories on compact manifolds. There is no Wick rotation infrastructure for the Bargmann graph, and Sprint Track RH-B (April 2026) closed the Wick-rotation route on the Coulomb side as a clean negative.

---

## §4. Recommendation

**Recommended:** Section VIII appendix on Paper 32 (spectral triple).

**Rejected alternatives:**

- *New Paper 37 synthesis*: tempting but premature. A standalone synthesis paper would need a load-bearing structural claim beyond "the three SM gauge groups appear on three different GeoVac sub-manifolds." That claim by itself is a catalogue, not a theorem; the four gaps in §3 mean the catalogue does not yet make a unified statement. A paper-length treatment would need to show either (a) progress on the cross-manifold obstruction (a single triple containing all three Wilson constructions, even tentatively), or (b) a clean structural reason why the three sectors must live on different sub-manifolds in GeoVac (and what that says about how GeoVac differs from Connes' SM lineage). Neither (a) nor (b) is in hand.

- *Memo only*: too weak. The unified kinetic coefficient 1/(4·N_c) verified independently across three constructions is a substantive cross-paper structural fact; the lineage placement in Marcolli-vS / Perez-Sanchez is published context that CLAUDE.md §1.5 wants to make explicit. Both deserve a permanent home in the paper record, not just in `debug/`.

**Why Paper 32 §VIII is the right venue:**

Paper 32 is the synthesis paper that explicitly constructs the GeoVac spectral triple (A_GV, H_GV, D_GV) and reads Papers 25 / 28 / 30 / 31 as projections of the same triple. The §VIII case-exhaustion theorem (added v2.27.3 in Sprint TS-E1) already states the three-mechanism structure of π-sources via Paper 34 projections. A natural follow-on appendix would:

1. Tabulate the three Wilson constructions as instances of the gauge-network framework, with the unified 1/(4·N_c) coefficient as the cross-construction structural identity.
2. Place them in the Marcolli-vS / Perez-Sanchez lineage (replacing the existing Paper 25 / Paper 30 / Paper 31 cross-references with a unified citation block).
3. State the four-gap analysis from §3 above as the honest structural distance to Connes' SM, with the cross-manifold gap (§3.4) named as the dominant obstruction.
4. Identify the S³ U(1) × SU(2) co-location as the most concrete near-term unification target (already on the same graph, same triple structure differs only at the gauge-group level).

This treatment lives naturally inside Paper 32's spectral-triple framing, requires no new computational sprint, and avoids the over-claiming risk a standalone Paper 37 would carry. It also closes a thread that is currently open across three separate papers (25, 30, ST-SU3 memo) without committing to a unification narrative the framework does not yet support.

**Estimated length:** ~3–4 pages added to Paper 32 §VIII as a new subsection (call it §VIII.B "The unified U(1) × SU(2) × SU(3) gauge content and its limits"). PI to confirm before any paper edit per CLAUDE.md §13.7.

---

## §5. Files referenced

- `papers/group5_qed_gauge/paper_25_hopf_gauge_structure.tex` (Paper 25, U(1) on S³)
- `papers/group5_qed_gauge/paper_30_su2_wilson.tex` (Paper 30, SU(2) on S³ = SU(2))
- `papers/group1_operator_algebras/paper_32_spectral_triple.tex` (Paper 32, recommended host)
- `geovac/su2_wilson_gauge.py` (Paper 30 implementation)
- `geovac/su3_wilson_s5.py` (ST-SU3 implementation)
- `debug/st_su3_wilson_memo.md` (ST-SU3 sprint memo)
- `debug/s5_gauge_structure_memo.md` (Sprint 5 Track S5 prior negative)
- `debug/wh1_round1_connes_vs_pdf_verification.md` (WH1 Marcolli-vS lineage placement)

External references: Marcolli–van Suijlekom, J. Geom. Phys. 75 (2014), arXiv:1301.3480; Perez-Sanchez, arXiv:2401.03705 and arXiv:2508.17338.
