# Sprint L0 — 28-Projection Lorentzian Transfer Audit

**Date:** 2026-05-16
**Author:** L0 PM (Claude)
**Source files persisted:** `debug/lorentzian_partition_transfer_memo.md` (L3 input verbatim); `debug/data/lorentzian_l0_audit.json` (machine-readable 28-row table)
**Trigger:** Nieuviarts 2502.18105v3 (May 2025) twisted-spectral-triple Riemannian → pseudo-Riemannian morphism. Trigger fired per `debug/lorentzian_literature_update_2026_05_16.md`; revisits the May 2026 NO-GO scoping at the documentation level.
**Status:** First-pass documentation deliverable. Paper-edit application is the load-bearing output.

---

## §1. Scope and discipline

This audit classifies all 28 entries of Paper 34's projection dictionary against the Riemannian → Lorentzian extension. Four classification buckets:

| Bucket | Definition |
|:-------|:-----------|
| **TRANSFERS_FREELY** | Signature-blind: lives on the V (variable) or D (dimension) axis of Paper 34 only, with no Riemannian volume form, no Mellin transform, no Latrémolière propinquity. Holds at (3,0) and (3,1) identically. |
| **WICK_MAP_FREE** | Reduces to the M1 sub-mechanism of the master Mellin engine via the four-witness Wick-rotation theorem (Hawking + Sewell + Bisognano-Wichmann + Unruh, codified Sprint Unruh-pendant 2026-05-10). β = 2π/(surface gravity); 2π = M1 Vol(S¹). Bridge supplied by published QFT machinery at the metric-functional level. |
| **EUCLIDEAN_SPECIFIC** | Requires Krein-space machinery, Lorentzian-native heat-trace formalism, or Lorentzian propinquity. Named blocker B1/B2/B3. |
| **MIXED** | Parts transfer freely, parts require extension. Specify which sub-mechanisms transfer and which are blocked. |

This is a **first-pass classification by an analyst (the L0 PM), not a theorem.** Each row's confidence rating reflects how robust the classification is. Sprint L1 (modular Hamiltonian computation) and Sprint L2 (BBB Krein lift + Connes axiom audit at (3,1)) will refine these classifications with actual computations. Several MIXED rows may shift in either direction.

The classification uses the Nieuviarts 2502.18105v3 trigger-firing literature update; whether the Nieuviarts twist morphism applies cleanly to GeoVac S³ = SU(2) is the parallel L2-Nieuviarts-scoping pass and is NOT assumed here.

---

## §2. Headline distribution

Of the 28 projections (count by bucket; full per-row classification in §3):

- **TRANSFERS_FREELY: 17 projections** (61%) — these inherit the universal sector of Paper 31's A/D partition, projected to the Sig/Op axis of the Riemannian → Lorentzian extension. They are signature-blind because they live on the SU(2) angular algebra (Wigner 3j/D, Gaunt termination, bipolar harmonic, Young tableau symmetry projection), on parameter-only mechanisms (rest-mass, Drake-Swainson asymptotic subtraction, Phillips-Kleinman), on consumption of external scalar inputs from QCD (charge radius, Zemach radius, quadrupole moment), or on the state-side complement of the dictionary (apparatus identity / von Neumann entropy, which is PSLQ-disjoint from the master Mellin engine per Sprint TD Track 5).
- **WICK_MAP_FREE: 4 projections** (14%) — observation/temporal-window §III.15, vector-photon promotion §III.11, gauge choice §III.26, Wick rotation §III.27 itself. All carry the M1 transcendental signature 2π = Vol(S¹) and inherit the Lorentzian reading via the four-witness theorem at no structural cost at the metric-functional level. §III.27 is the bridge projection by construction.
- **EUCLIDEAN_SPECIFIC: 5 projections** (18%) — Fock conformal §III.1, Hopf bundle §III.2, stereographic §III.4, spectral action §III.6, Camporesi-Higuchi spinor lift §III.7. These require a Krein-space spectral triple at (3,1) per BBB 2018 (blocker B2) or Lorentzian-native heat-trace formalism (blocker B3). Sprint L2 territory.
- **MIXED: 2 projections** (7%) — Wilson plaquette §III.10 (combinatorial vocabulary universal, Marcolli-vS YM correspondence Euclidean), Breit retardation §III.16 (angular content TRANSFERS, kinematic content needs Lorentzian re-derivation but result stays in ℚ(α)).

The distribution is informative: **the framework's structural-skeleton scope reading (CLAUDE.md §2, 2026-05-07) survives the (3,1) extension cleanly.** 61% transfers freely + 14% Wick-maps cleanly = 75% of the dictionary has a Lorentzian reading available at no further structural cost. The remaining 25% sits in the EUCLIDEAN_SPECIFIC / MIXED bucket, which is exactly the load-bearing structural extension content that defines Sprint L1, L2, L3.

---

## §3. The 28-row classification table

(Source: `debug/data/lorentzian_l0_audit.json`. Confidence ratings: high / medium-high / medium / medium-low. One-line mechanism per row.)

| # | Projection | Sec | V | D | T | Class | Conf | Mechanism / Blocker |
|:-:|:-----------|:----|:--|:--|:--|:-----:|:----:|:--------------------|
| 1 | Fock conformal | III.1 | Z, E | energy | κ; π via Vol(S³) | E.S. | high | Stereographic R³→S³ uses Riemannian metric; pseudo-Riemannian sphere quotient unpublished. B2. |
| 2 | Hopf bundle | III.2 | α | dimensionless | π = Vol(S²)/4 | E.S. | high | Riemannian fibration; pseudo-Riemannian Hopf base unpublished. |
| 3 | Bargmann-Segal | III.3 | ℏω | energy | none (π-free) | T.F. | medium | First-order complex Euler op, linear-affine projection; π-free certificate is signature-blind at the discrete-graph level. |
| 4 | Stereographic | III.4 | r (length) | length | conformal factor | E.S. | high | Riemannian conformal map; Lorentzian conformal compactification differs. B2. |
| 5 | Sturmian reparam | III.5 | Z, n (re-use) | preserves | rational-preserved | T.F. | high | Algebraic re-labeling on algebra A side; signature-blind. |
| 6 | Spectral action | III.6 | Λ, α | energy via Λ | √π·ℚ SD; π^{2k}·ℚ obs | E.S. | high | Tr f(D²/Λ²) is Euclidean heat-kernel (M2). B3. |
| 7 | C-H spinor lift | III.7 | α via (Zα)² | dimensionless | half-int Hurwitz → odd ζ; G, β(4) | E.S. | med-high | Riemannian spin manifold; Krein spinor at (3,1) unpublished. B2. |
| 8 | Wigner 3j | III.8 | none | dimensionless | ℚ[√(2k+1)]_k | T.F. | high | SU(2) rep theory; signature-blind. Canonical universal sector. |
| 9 | Wigner D rotation | III.9 | R_AB | length | ℚ[√2, √3, √6] | T.F. | high | Spatial rotations between nuclear centers; SU(2) rep theory. |
| 10 | Wilson plaquette | III.10 | β = 1/g² | dimensionless | SU(2) Haar | MIXED | medium | Combinatorial vocabulary universal; Marcolli-vS = YM Euclidean. |
| 11 | Vector-photon | III.11 | none | dimensionless | 1/(4π) per loop | W.M.F. | medium | 1/(4π) is M1 Hopf-base measure → Wick-maps; photon algebra at (3,1) needs B2. |
| 12 | Mol-frame hypersph | III.12 | R | length | piecewise-smooth in R | T.F. | med-high | Angular separation on algebra A side; relativistic molecular framework needs §III.7 lift. |
| 13 | Drake-Swainson | III.13 | K (cancels) | dimensionless | flow tier; rational denom | T.F. | high | Central-tier (parametric) per Track TS-C; signature-blind. |
| 14 | Rest-mass | III.14 | m | mass | trivial (ring-preserving) | T.F. | high | D²→D²+m²𝟙 is central; canonical bridge to relativistic kinematics. |
| 15 | Observation / temporal-window | III.15 | β (or T = 1/β) | time | 2π·ℚ per Matsubara | W.M.F. | high | THE Wick-engine: Matsubara 2πk/β is the four-witness M1 mechanism. |
| 16 | Breit retardation | III.16 | m_l/m_n | energy | α⁴·ℚ | MIXED | medium | Angular content TRANSFERS; kinematic content needs Lorentzian re-derivation but stays in ℚ(α). |
| 17 | Charge density (Foldy/Friar) | III.17 | r_E | length | ℚ(α) | T.F. | high | Layer-2 scalar input from QCD; signature-blind at framework boundary. |
| 18 | Magnetization (Zemach) | III.18 | r_Z | length | ℚ(α) | T.F. | high | Same as Foldy/Friar; Layer-2 scalar. |
| 19 | Tensor multipole (Q_N) | III.19 | Q_N | length² | ℚ via Wigner 3j | T.F. | high | Rank-2 via Wigner-Eckart; pure SU(2) at angular level. |
| 20 | Phillips-Kleinman | III.20 | {φ_c, E_c} | dimensionless | ring-preserving | T.F. | high | Pseudopotential projection consuming core orbitals; signature-blind. |
| 21 | Multipole / Gaunt termination | III.21 | L (integer) | preserves | ring-preserving | T.F. | high | L_max = 2 l_max is Wigner 3j triangle inequality; pure SU(2). |
| 22 | Bipolar / Drake combining | III.22 | (k₁, k₂, K) | preserves | pure rational | T.F. | high | Multi-electron angular content; pure SU(2). |
| 23 | Symmetry / Young | III.23 | irrep tableau λ | dimensionless | none | T.F. | high | S_N permutation rep theory; Pauli exclusion is signature-blind. |
| 24 | Adiabatic / BO | III.24 | R; m_e/M_n | length | none at proj step | T.F. | med-high | Scale separation; relativistic BO is well-defined. |
| 25 | Coupled-channel | III.25 | R, ν | length | algebraic-implicit; piecewise | T.F. | medium | ODE machinery signature-blind; algebraic-implicit π inherits from Fock §III.1 (possible refinement to MIXED). |
| 26 | Gauge choice | III.26 | none | preserves | gauge-dep prop / gauge-inv obs | W.M.F. | med-high | 1/(4π) Coulomb gauge = M1; Lorentzian gauge fixing is standard QFT. |
| 27 | Wick rotation | III.27 | signature | preserves | 2π·ℚ via Vol(S¹) (M1) | W.M.F. | high | THE BRIDGE. By construction. Sprint L1 deepens to operator-system level. |
| 28 | Apparatus identity | III.28 | ρ; β = 1/T | preserves | von Neumann; PSLQ-disjoint | T.F. | med-high | State-side complement of dictionary; PSLQ-disjoint from M1∪M2∪M3 per Sprint TD Track 5. |

**Abbreviations:** T.F. = TRANSFERS_FREELY; W.M.F. = WICK_MAP_FREE; E.S. = EUCLIDEAN_SPECIFIC; MIXED = MIXED.

---

## §4. Sequencing recommendation

The audit identifies four natural sprint tiers:

**Tier 0 (this sprint, L0):** Documentation only. Sixteen TRANSFERS_FREELY projections inherit Paper 31's universal sector cleanly; defer to Paper 31 §8 + Paper 34 §V.E placement (see §6 below). One additional TRANSFERS_FREELY entry (§III.28 apparatus identity, state-side) deserves its own paragraph because it surfaces the state-side sub-dictionary as a structurally distinct piece of the partition. **Status: this memo + paper edits = closure for the 17 leaves.**

**Tier 1 (Sprint L1, ~4-8 weeks, MEDIUM-HIGH priority):** Operator-system literal-identification of the Wick-rotation projection §III.27. The principal falsifier is: compute the modular Hamiltonian K on T_{n_max} for a Rindler-wedge-restricted Camporesi-Higuchi state, verify σ_{i·2π} = identity at finite n_max, then take the GH limit using Paper 38's lemmas. **Does NOT require B2 (Krein lift)** — operator-system extension is internal to the existing Riemannian propinquity machinery, just specialized to the wedge-restricted algebra. One operator-level proof lifts all four Wick-rotation faces simultaneously (Hawking, Sewell, BW, Unruh) because the bridge runs through the same KMS-Tomita-Takesaki algebra. Already named in Paper 32 §VIII rem:bisognano_wichmann_reading and Paper 34 §VIII open question entry on operator-system-level Lorentzian extension.

**Tier 2 (Sprint L2, ~3-6 months):** BBB Krein lift + Connes axiom audit at (3,1) signature. Gates: blocker B2. Five projections in the EUCLIDEAN_SPECIFIC bucket (§III.1 Fock, §III.2 Hopf, §III.4 stereographic, §III.6 spectral action, §III.7 spinor lift) all wait here. The Nieuviarts 2502.18105v3 trigger-firing literature update suggests a possible morphism-based shortcut; the parallel L2-Nieuviarts-scoping pass tests whether the twist morphism applies cleanly to GeoVac S³ = SU(2) (odd-dim caveat to verify; BBB only treats (m,n) ∈ (ℤ/8)² explicitly). Side benefit of Sprint L2: refines two MIXED rows (§III.10 Wilson plaquette, §III.16 Breit retardation) into definitive classifications.

**Tier 3 (Sprint L3, ~6-12 months):** Lorentzian propinquity construction (blocker B1). Becomes load-bearing only if the operator-system Lorentzian extension (deepening of §III.27 + §III.7) requires Lorentzian Latrémolière propinquity to take the GH limit at finite n_max with Krein-space Hilbert space. May be shortened significantly by Nieuviarts morphism if the L2 scoping verifies applicability. **Status: probably skippable for Sprint L1 if the modular Hamiltonian closure works on the existing Riemannian truncated triple; becomes load-bearing for full Sprint L2 closure.**

**Trivially closed at (3,1) — Sprint L2 prediction:** M3 sub-mechanism of the master Mellin engine (vertex-parity Hurwitz / Dirichlet-L content in Camporesi-Higuchi vertex sums, Paper 28 §QED-vertex) becomes structurally empty on chirality-symmetric spectrum because the (3,1) BBB Table 2 entry flips the {J, γ_5} anticommutation table. This is a Sprint L2 prediction worth verifying as a consistency check (it would simplify the master Mellin engine's three-mechanism partition to two on the Lorentzian side, which is a non-trivial structural observation if true).

---

## §5. Where a future standalone paper would draw from

Per PI instruction, **no new paper is drafted at this stage.** The current edits to Papers 31, 34, 32 are the right level. For the record, when Sprint L1 and Sprint L2 land actual computations, a standalone paper at the natural NCG venue (J. Geom. Phys. or Letters in Mathematical Physics) would draw from:

- **This L0 audit memo + JSON** (the 28-row partition, classification mechanism, and sequencing logic).
- **Paper 31 §8 (post-L0 extension)** — the Sig/Op partition framework as the structural backbone.
- **Paper 38 §6.3** — the Lorentzian propinquity open question; Sprint L1 closure or L3 construction would land Theorem-level content here.
- **Paper 34 §III.27** — the Wick-rotation projection with the four-witness theorem as published anchor.
- **Paper 32 §VIII case-exhaustion theorem + GH-convergence theorem** — the structural framework that proves WH1 (Riemannian PROVEN); Lorentzian extension lifts to a (3,1) version under B2.
- **Sprint L1 computational results** (modular Hamiltonian on T_{n_max}, σ_{i·2π} = identity).
- **Sprint L2 computational results** (Connes axiom audit at (3,1) per BBB sign table; sub-mechanism trivialization predictions).
- **Bizi-Brouder-Besnard 2018 (arXiv:1611.07062)** as the published Krein-spectral-triple framework.
- **Nieuviarts 2502.18105 (May 2025)** as the morphism-based emergence prescription (pending L2 scoping verification of applicability).

The natural title for the future paper would be something like *"Lorentzian Extension of the GeoVac Spectral Triple: Wick-Rotation Bridge and Operator-System Closure on T_{S³}"*. The current sprint produces the structural skeleton; that paper would carry the actual L1 + L2 computations.

---

## §6. Paper-update rationale and placement decisions

### §6.1 Paper 31 — Sig/Op partition section

Paper 31 (Universal/Coulomb Partition) is the natural anchor for the (Riemannian → Lorentzian) Sig/Op partition because the structural template (universal sector lives on A; specific sector lives on D) extends verbatim to (signature-blind sector lives on V/D axes of Paper 34; signature-specific sector requires Krein/Lorentzian extension). The natural location is a **new §8 "Signature partition"**, placed after the existing §7 "Implications / What requires per-system rederivation" and before §9 "Open Questions" (current numbering: §7 Implications, §8 Open Questions — so the new section becomes §8, pushing Open Questions to §9). This placement is parallel to how §5 (Three-Layer Sharpening / Sprint 5) extended the partition from two to three layers.

Forward references from new §8: Paper 32 §VIII (case-exhaustion theorem + WH1), Paper 34 §V.E or §VIII open-question table (depending on §6.2 decision), Paper 38 §6.3 (Lorentzian propinquity open question as blocker B1).

### §6.2 Paper 34 — placement decision (§V.E vs. §VIII)

The 28-projection transfer table can sit naturally in either of two places:

**Option A: New §V.E "Lorentzian transfer audit (Sprint L0)".** Parallel to §V.B (off-precision matches), §V.C (Roothaan autopsies), §V.D (literature convention exposures). This treats the transfer audit as a standalone partition table living in the Empirical Matches Catalogue's running-catalogue layer. Pros: parallel structure with the other running catalogues; future updates (Sprint L1 / L2 results) flow into the same section. Cons: §V is "Empirical Matches Catalogue" by name; a partition-table audit isn't strictly an empirical match.

**Option B: Extend §VIII open question 6 (Operator-system-level Lorentzian extension).** The existing open-question entry already names operator-system Lorentzian extension as the principal falsifier for the §III.27 structural-correspondence verdict. Adding the 28-row transfer table as a sub-table of this open-question entry sharpens the open question into a structured roadmap. Pros: naturally extends an existing entry; no re-numbering. Cons: §VIII is "Open Questions" and the audit is now structured / partial-closure content, not just an open question.

**Decision (PM): Option A.** The 28-row table is a structural partition statement, not an open question. Open Questions §6 still names the operator-system extension as the principal Sprint L1 falsifier, but the audit itself is a partition deliverable. §V.E is the cleaner home because (i) the section will accumulate Sprint L1 + L2 results as they land, parallel to how §V.C accumulates new Roothaan autopsies; (ii) the explicit "running catalogue" header in §V.D is already the template for accumulating sprint-by-sprint results; (iii) a forward reference from §V.E to §VIII open question 6 keeps the falsifier visible. §VIII §6 gets a brief two-sentence cross-reference back to §V.E.

### §6.3 Paper 32 §VIII subsection

Paper 32 §VIII already contains the master Mellin engine case-exhaustion theorem (thm:pi_source_case_exhaustion), the GH-convergence theorem (thm:gh_convergence), the BW reading remark (rem:bisognano_wichmann_reading with the Unruh pendant extension), the SM unified gauge appendix §VIII.B, and the Sprint H1 verdict §VIII.C. A new subsection §VIII.E "Lorentzian transfer audit (Sprint L0)" naturally extends this with: (i) cross-reference to the new Paper 34 §V.E for the per-projection partition; (ii) a sharpening of the case-exhaustion theorem's M1 sub-mechanism as the load-bearing free-transfer engine via the four-witness theorem; (iii) explicit statement that M3 becomes trivially closed on chirality-symmetric spectrum under BBB Table 2 at (3,1), as a Sprint L2 prediction.

### §6.4 CLAUDE.md §2 entry

A sprint summary bullet documenting the L0 sprint, its inputs (lit update trigger + L3 partition memo + L4 BH entropy scoping), tracks executed, headline findings (61%/14%/18%/7% bucket split), and the paper-edit list. Standard sprint-bullet template per the recent (May 16) XCWG arc and Paper 40 entries.

### §6.5 Memory file

A new memory file `memory_lorentzian_l0_audit.md` capturing the partition's structural reading for future PM session restoration. CLAUDE.md MEMORY.md gets a one-line bullet pointing to the new memory file.

---

## §7. Honest scope and what L0 does not claim

- This is documentation, not theorem. None of the 28 classifications is a proved structural statement; they are first-pass analyst classifications.
- The 75% (TRANSFERS_FREELY + WICK_MAP_FREE) figure is encouraging but should not be over-interpreted. The remaining 25% (EUCLIDEAN_SPECIFIC + MIXED) contains the structurally consequential projections: spectral action, Fock conformal, Hopf bundle, spinor lift, stereographic. Without those, GeoVac at (3,1) cannot reproduce α-conjecture, calibration κ, or the Hopf-base measure observables that hold the framework together.
- The trigger-firing event (Nieuviarts 2502.18105 morphism) is named as making B1 possibly shorter, NOT as confirmed to apply to GeoVac. The applicability scoping (odd-dim caveat to verify on S³ = SU(2)) is the parallel L2-Nieuviarts-scoping pass.
- The Bekenstein-Hawking scoping memo (`debug/bekenstein_hawking_sprint_scoping_memo.md`) is correctly DEFERRED. BH entropy is Euclidean spectral-action work and does not require Lorentzian extension; it is GO-WITH-PREREQUISITES (8-14 weeks) and outside Sprint L0/L1/L2/L3 scope.
- Sprint L1's operator-system closure does NOT promote the Wick-rotation projection §III.27 from "structural correspondence" to "literal identification" until the falsifier (σ_{i·2π} = identity on T_{n_max}) is computed; the L0 audit honestly notes this gap.

---

## §8. Files modified summary

1. `debug/lorentzian_partition_transfer_memo.md` — persisted L3 input (created)
2. `debug/lorentzian_l0_audit_memo.md` — this memo (created)
3. `debug/data/lorentzian_l0_audit.json` — 28-row machine-readable table (created)
4. `papers/core/paper_31_universal_coulomb_partition.tex` — new §8 Signature partition + §7 forward reference + bibliography additions
5. `papers/observations/paper_34_projection_taxonomy.tex` — new §V.E Lorentzian transfer audit + §VIII open question 6 cross-reference back
6. `papers/synthesis/paper_32_spectral_triple.tex` — new §VIII.E subsection on Lorentzian transfer audit
7. `CLAUDE.md` §2 — new sprint bullet for Sprint L0
8. `MEMORY.md` — new memory file pointer
9. `memory/memory_lorentzian_l0_audit.md` — memory file

---

## §9. Forward plan flags

- **Sprint L1 (modular Hamiltonian, 4-8 weeks):** Compute K on T_{n_max} for a Rindler-wedge-restricted Camporesi-Higuchi state. Falsifier: σ_{i·2π} = identity at finite n_max. Closure promotes Paper 34 §III.27 from structural correspondence to literal identification and clears WICK_MAP_FREE bucket fully.
- **Sprint L2-Nieuviarts-scoping (1-2 weeks, parallel):** Test whether Nieuviarts 2502.18105 twist morphism applies cleanly to GeoVac S³ = SU(2) (odd-dim caveat). If yes, shortens Sprint L2 budget significantly. If no, defaults to BBB 2018 standard prescription.
- **Sprint L2 full (3-6 months):** BBB Krein lift + Connes axiom audit at (3,1). Refines 5 EUCLIDEAN_SPECIFIC + 2 MIXED rows to definitive classifications. Predicts M3 trivialization on chirality-symmetric spectrum (consistency check).
- **Sprint L3 (6-12 months, possibly skippable):** Lorentzian propinquity construction. Becomes load-bearing only if Sprint L1/L2 closure requires Krein-Latrémolière propinquity.
- **Standalone Lorentzian-extension paper:** Drafted after L1 + L2 land. Title placeholder *"Lorentzian Extension of the GeoVac Spectral Triple"*. Outline above §5.

End of memo.
