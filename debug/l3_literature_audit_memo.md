# Sprint L3 — Literature Audit (Track 2)

**Date:** 2026-05-17
**Author:** L3-scoping (Track 2 literature search)
**Trigger:** Sprint L2 closed-at-finite-cutoff (2026-05-16/17). Need to scope whether Sprint L3 (continuum Lorentzian propinquity construction) has any published precedent, partial precedent, or named entry point. Decision target: 1–3 month Sprint L3a with specific scaffold, vs. multi-month solo new mathematics.
**Method:** WebSearch + WebFetch, 14 searches across arXiv, 9 PDF/abstract verifications. Cross-checked against May-2026 reference list in Paper 42 §3.2, `debug/lorentzian_l0_audit_memo.md`, and `debug/sprint_l2_synthesis_memo.md`.
**Honest scope:** abstract-level literature scan, not full PDF reads. Verdicts are based on titles, abstracts, and author research profiles. Findings should be sanity-checked at proof level if any sprint launches from them.

---

## §1. Lorentzian / pseudo-Riemannian propinquity status

### Verdict: NO PUBLISHED LORENTZIAN PROPINQUITY EXISTS AS OF MAY 2026.

Latrémolière's propinquity framework — both the original (arXiv:1302.4058, 2013) and the metric-spectral-triple Gromov–Hausdorff propinquity (arXiv:1811.10843, 2018; Adv. Math. 2022) — is strictly Riemannian. There is no published or pre-printed Lorentzian/pseudo-Riemannian extension authored by Latrémolière himself.

**Latrémolière 2024–2026 publication trajectory** (all Riemannian):
- arXiv:2404.00240 (Farsi–Latrémolière, "Collapse in Noncommutative Geometry and Spectral Continuity", v3 Oct 2025): Riemannian spin manifolds, U(1) principal bundles over Riemannian. **No Lorentzian content.**
- arXiv:2504.11715 (Farsi–Latrémolière, "Continuity for the spectral propinquity of the Dirac operators associated with an analytic path of Riemannian metrics", April 2025): the title says it — analytic path of *Riemannian* metrics.
- Latrémolière 2026 book (existence already cited in CLAUDE.md §2 L2-A scoping): the published-book successor to the 2017 propinquity arXiv writeup, still Riemannian.

**Critical clarification on the "Lorentzian GH" literature.** There is a vibrant 2024–2025 literature on Lorentzian/measured Gromov–Hausdorff convergence — but it is **synthetic Lorentzian metric geometry**, not noncommutative-geometry propinquity. The two strands of research currently run in disjoint communities and have no published bridge:

| Paper | arXiv | Date | Authors | Setting | Relevance to L3 |
|:------|:------|:-----|:--------|:--------|:----------------|
| Lorentzian metric spaces & GH-convergence | 2209.14384 | 2022 → 2025 (v3) | Müller (+collaborators) | Synthetic Lorentzian length spaces | TANGENTIAL |
| Gromov-Hausdorff metrics & dimensions of Lorentzian length spaces | 2209.12736 | 2022 | (multi-author) | Synthetic, dimension theory | TANGENTIAL |
| Lorentzian metric spaces & GH-convergence: unbounded case | 2412.04311 | Dec 2024, rev May 2025 | Bykov, Minguzzi, Suhr | Synthetic Lorentzian, no spectral triples | TANGENTIAL |
| Lorentzian Gromov-Hausdorff convergence and pre-compactness | 2504.10380 | April 2025, rev Dec 2025 | Mondino, Sämann | Synthetic, causal-diamond ε-nets | TANGENTIAL |
| Gromov's reconstruction theorem & measured Lorentzian GH | 2506.10852 | June 2025 | Braun, Sämann | Synthetic measured Lorentzian | TANGENTIAL |
| Results on Lorentzian metric spaces | 2510.24423 | Oct 2025 | Minguzzi | Synthetic, Cauchy time functions | TANGENTIAL |
| Intrinsic Timed Hausdorff Convergence | 2511.18389 | Nov 2025 | (multi-author) | Synthetic | TANGENTIAL |

**Why TANGENTIAL not USEFUL_TOOL:** these papers work with two-point Lorentzian distance functions on synthetic spaces (causal diamonds, time-separation functions, length-space axioms). They do NOT define an operator-system / C*-algebra / spectral-triple analog. There is no Lipschitz seminorm from a Dirac operator, no Lip-norm convergence, no Kantorovich–Rubinstein dual lift to spectral triples. The Mondino–Sämann construction (ε-nets of causal diamonds) is the closest in spirit to Latrémolière propinquity but lives entirely on the metric-space side.

**The natural reading:** Sprint L3 would need to *create* the bridge between the Mondino–Sämann–Braun synthetic Lorentzian GH program (mature, 2022–2025, ~7 papers) and the Latrémolière–Farsi metric-spectral-triple program (mature, 2013–2025, ~30+ papers). Both halves exist; the bridge does not. This is original NCG-mathematics work, but the literature scaffold on each side is much stronger than the May-2026 L2-A audit credited it. **The "no published scaffold" framing of Sprint L0 is correct on the propinquity side but understates the Lorentzian-synthetic-GH side.**

### Closest analog

The Mondino–Sämann 2504.10380 (April 2025) causal-diamond ε-net definition is the most propinquity-analogous Lorentzian construction. It builds GH-convergence from a *family of distances* (one per causal-diamond scale), which is structurally parallel to Latrémolière's tunneling-pair construction. A Sprint L3a entry point would be: lift the causal-diamond ε-net to operator-system data by replacing the time-separation function with a Krein-space modular flow generator.

---

## §2. Krein-space spectral triple convergence status

### Verdict: NO PUBLISHED CONVERGENCE THEOREM FOR FINITE-CUTOFF KREIN TRIPLES. Krein-space spectral triples themselves are a moderately active 2018–2026 area.

**Bizi–Brouder–Besnard 2018** (arXiv:1611.07062) remains the canonical Krein-classification reference (the (m, n) ∈ (ℤ/8)² classification). I found NO direct 2024–2026 follow-up paper extending BBB's framework with operator-system truncations, convergence theorems, or Krein-positive-cone state spaces.

**What does exist** (USEFUL_TOOL tier):
- Pierre Martinetti, "Twisted Standard Model and its Krein structure -- in memoriam Manuele Filaci" (arXiv:2603.03216, March 2026): CIRM Marseille proceedings. Reviews SM contributions in NCG, examines twisting the spectral triple, and shows twisting **induces an inner product transforming the Hilbert space into a Krein space**. Does NOT address truncation convergence. Connects Krein structure to twisted-spectral-triple emergence (Devastato–Lizzi–Martinetti lineage).
- Bär–Strohmaier global index theorem extensions (van den Dungen–Ronge 2021; Damaschke 2021; Shen–Wrochna 2022; Bär–Strohmaier 2024). These are Lorentzian-side index theorems that generalize Strohmaier 2006 — but they do not produce a GH-convergence framework. **USEFUL_TOOL** for index-theoretic content; not a propinquity scaffold.
- "A Lorentzian Equivariant Index Theorem" (arXiv:2602.16547, 2026): further extension of the Bär–Strohmaier program. Does not address convergence.

**What does NOT exist:**
- A "Latrémolière propinquity for Krein-space spectral triples."
- A "Connes–van Suijlekom truncation theorem on a Krein triple."
- A "Berezin reconstruction on a Krein space."
- A spectral truncation Λ on a Krein operator at signature (3, 1) with bounded reach as Λ → ∞.

**Reachable / open question:** can Connes–vS 2021 (CMP, spectral truncations on S¹/Toeplitz) be extended to the (m, n) = (4, 6) Bizi–Brouder–Besnard Krein triple? This is a self-contained mathematical question with all ingredients on the shelf. Nobody has written it. **Open.**

The Krein-quantization side (arXiv:2505.19632, "Krein space quantization and New Quantum Algorithms", May 2025) is about quantum algorithms for non-unitary evolution — IRRELEVANT to spectral-triple convergence, despite the keyword overlap.

---

## §3. Nieuviarts twist morphism status

### Verdict: NO ODD-DIM EXTENSION HAS APPEARED. S³ DISQUALIFICATION HOLDS THROUGH MAY 2026.

**Latest Nieuviarts arXiv versions** (cross-checked against Paper 42 §3.2 footnote bibitems):

| Paper | arXiv | Latest version | Date | Dim restriction |
|:------|:------|:---------------|:-----|:----------------|
| Signature change by a morphism of spectral triples | 2402.05839 | v6 | March 2, 2026 | **Still even-dim only** ("In the case of even-dimensional manifolds, we demonstrate how this construction implements a local signature change via the parity operator induced by the twist") |
| Emergence of Lorentz symmetry from almost-commutative twisted spectral triple | 2502.18105 | v3 | May 6, 2025 | Not updated since May 2025 |
| Emergence of Time from a Twisted Spectral Triple in Almost-Commutative Geometry | 2512.15450 | v2 | May 11, 2026 | Proceedings synthesis. Same dimensional scope as 2502.18105. |

**Critical:** the v6 (March 2026) revision of 2402.05839 still explicitly restricts to even dimensions. The author's own justification (quoted in Paper 42 §3.2 footnote): *"the fact that no such procedure for twisting odd-dimensional manifold's spectral triple has been found will be one of the justifications to focus on the study of even-dimensional manifolds"*. **Two years on, the odd-dim restriction is acknowledged-and-unresolved**, not provisional pending v_next. Sprint L2's NO-GO verdict for GeoVac S³ = SU(2) holds.

**Even-dim applications watch-list** (would shorten any future cross-manifold Lorentzian sprint, e.g. W2b T_S³ ⊗ T_S⁵ where S⁵ × S¹ is even-dim, or any sprint touching T_S^{2n}):
- 2512.15450 v2 (May 2026): proceedings synthesis of the entire Nieuviarts arc. Best entry point for someone wanting to understand the construction.
- 2603.03216 (Martinetti, March 2026): explicitly extends the twist-morphism / Krein-structure dictionary to the Standard Model setting.

**Verdict for GeoVac L3:** the Nieuviarts shortcut named in the Sprint L0 audit (`debug/lorentzian_l0_audit_memo.md` §4 Tier 3) is **structurally unavailable** to GeoVac's S³ = SU(2) sector. The author has had 2 years to find an odd-dim extension; none has appeared. Sprint L3 cannot use this shortcut. Direct BBB Krein lift remains the path.

---

## §4. Connes–Strohmaier 2006 application status

### Verdict: ACTIVE BUT INDEX-THEORETIC, NOT CONVERGENCE-THEORETIC. NOT the natural scaffold for L3.

Strohmaier 2006 ("On Noncommutative and Pseudo-Riemannian Geometry", arXiv:math-ph/0110001 lineage) established the temporal Lorentzian spectral triple framework. The 2024–2026 literature on this lineage has gone in two directions, neither matching L3's needs:

**(a) Bär–Strohmaier index-theoretic extensions** (the productive 2020–2026 thread):
- Braverman 2020
- van den Dungen–Ronge 2021
- Damaschke 2021
- Shen–Wrochna 2022
- Bär–Strohmaier 2024
- A Lorentzian Equivariant Index Theorem (arXiv:2602.16547, 2026)

These are Lorentzian / globally hyperbolic spacetime index theorems. They do not produce a Latrémolière-style propinquity convergence theorem. **USEFUL_TOOL** for any future GeoVac index-theoretic question on the Lorentzian side; not the scaffold for L3.

**(b) Causality-from-algebraic-data thread (Franco–Eckstein lineage):**
- Franco 2011 thesis (arXiv:1108.0592)
- Franco–Eckstein 2014 (arXiv:1409.1480, "Noncommutative geometry, Lorentzian structures and causality")
- Franco–Eckstein 2015 (arXiv:1508.01917, "Two roads to noncommutative causality")
- 2025 follow-ups in kappa-Minkowski space (arXiv:2503.24192)

This thread defines Lorentzian distance formulas on Krein-based spectral triples and recovers causality from algebraic data. Closer to what L3 wants — but still does not provide a propinquity-style convergence framework. **USEFUL_TOOL** for the algebraic-data scaffold; not L3's missing piece.

**The Connes–Strohmaier scaffold IS adequate for the static structure (Lorentzian spectral triple at fixed cutoff), which Sprint L2 already used via the van den Dungen 2016 Proposition 4.1 lift.** What is missing is the convergence-as-cutoff-removed direction.

**Verdict for L3:** Strohmaier 2006 lineage is the right vocabulary, but L3's missing piece is a convergence theorem, which the Strohmaier program has never produced.

---

## §5. van den Dungen 2016 successors

### Verdict: VAN DEN DUNGEN HAS PIVOTED TO INDEX THEORY. NO CONVERGENCE-OF-TRUNCATIONS DIRECTION.

van den Dungen's published 2024–2025 work has moved decisively into Dirac–Schrödinger index theory and the Callias theorem:

| Paper | Venue | Date |
|:------|:------|:-----|
| Dirac-Schrödinger operators, index theory and spectral flow | J. London Math. Soc., Vol. 112, e70301 | 2025 |
| Generalised Dirac-Schrödinger operators and the Callias Theorem | Forum of Math. Sigma, Vol. 13:e11 | 2025 |
| (van den Dungen–Ronge 2021) Bär–Strohmaier index extension | (index-theoretic, Lorentzian) | 2021 |

**No published successor to van den Dungen 2016 (arXiv:1505.01939) directly addresses:**
- Truncation of Krein spectral triples
- Convergence of finite-rank Krein operators
- Coupling Krein spectral triples to matter at the operator-system level

The 2016 paper itself remains the canonical reference for the Krein-spectral-triple lift via Proposition 4.1. Sprint L2-C used this verbatim. There is no published refinement of Prop 4.1 to a convergent-truncation framework.

**Mechanism:** van den Dungen's 2020–2025 trajectory has been into the Bär–Strohmaier–Callias index-theoretic continuation, NOT into the metric / convergence direction. This is fine for GeoVac in the static / fixed-cutoff sense (L2 is closed), but means the natural author who could extend Prop 4.1 to convergence is not currently working on it.

---

## §6. Best-bet mathematical entry points for Sprint L3

If a Sprint L3 launches, three layered entry points are realistic, in increasing order of solo-mathematics cost:

### Entry point A (CRITICAL_PRECEDENT, ~1–3 months): Connes–vS spectral truncation extended to BBB Krein triple.

**Specific question:** define the analog of Connes–vS 2021 truncated operator system $P_n A P_n$ on the (m, n) = (4, 6) Bizi–Brouder–Besnard Krein spectral triple. Prove that the propagation number remains finite, and that the analog of the prop = 2 Toeplitz match (WH1 Round 2 result) extends to the Krein setting.

**Why this is the cleanest first move:**
- All ingredients on the shelf (Connes–vS 2021, BBB 2018, Sprint L2's bit-exact Connes axiom audit at (4, 6))
- No need for Lorentzian propinquity yet (this is structure of the truncated operator system at fixed cutoff, not convergence)
- Direct extension of WH1 Round 2 + Sprint L2-D into a published-paper-quality theorem
- Sprint-scale (1–3 months)

**Who should write:** the Sprint L3a author, drawing on the WH1 Round 2 prop = 2 machinery and Sprint L2-D verified BBB sign signature.

### Entry point B (USEFUL_TOOL, ~3–6 months): Krein-positive cone as state space + Kantorovich–Rubinstein dual.

**Specific question:** the Sprint L2-C Krein-positive cone $\mathcal{K}^+$ is a half-dimensional sub-cone of the Krein space at every finite cutoff. Define a Wasserstein–Kantorovich distance on the Krein-positive states with a Lip-norm coming from $D_L$ (the L2-C Lorentzian Dirac).

**Why this is the second move:**
- Builds on Paper 38 Lemma L4 (Berezin reconstruction) — Berezin maps from Krein-positive states are well-defined, but their Wasserstein lift is open
- Closest published analog: Franco–Eckstein 2014 Lorentzian distance formula on Krein spectral triples + Connes–vS 2021 state-space convergence
- May land "Krein-Latrémolière propinquity at finite cutoff" as a finite-dim theorem before tackling the limit

**Who should write:** an L3b author drawing on Franco–Eckstein causality lineage + Sprint L2-C K^+ construction.

### Entry point C (multi-month, ~6–12 months): Bridge synthetic Lorentzian GH (Mondino–Sämann) to spectral-triple propinquity.

**Specific question:** Mondino–Sämann 2504.10380 defines Lorentzian GH convergence via causal-diamond ε-nets. Lift this to operator-system data: replace the time-separation function with a Krein modular flow generator from Sprint L2-E (the K_L^{α, W} integer-spectrum operator). Build a tunneling-pair construction that descends to the Mondino–Sämann ε-net at the synthetic level and to the Krein spectral triple at the operator level.

**Why this is the third move:**
- This IS the "create the bridge" original-NCG-mathematics task that Sprint L2-A audit budgeted at 6–12 months
- Both halves are mature; bridge is open
- Would produce a "Lorentzian propinquity" as a publishable framework, not just a GeoVac-internal construction
- High-value, high-risk

**Who should write:** a multi-month NCG-mathematics sprint, possibly in collaboration with Mondino, Sämann, or Latrémolière. The PI may want to consider whether outreach to those authors is appropriate before committing.

---

## §7. Risk-of-being-scooped assessment

### Active researchers in adjacent directions (May 2026)

| Researcher(s) | Direction | Could scoop L3? |
|:--------------|:----------|:----------------|
| Latrémolière (Denver) | Riemannian spectral propinquity | No — Lorentzian extension not on stated agenda |
| Farsi (Boulder) | Collapse, continuity, Riemannian | No — Lorentzian extension not on stated agenda |
| Mondino (Oxford) | Synthetic Lorentzian GH | **Possibly** — but operator-side bridge would require collaboration with NCG specialist |
| Sämann (Vienna) | Synthetic Lorentzian GH, measured | **Possibly** — same caveat |
| Bär (Potsdam), Strohmaier (Leeds) | Lorentzian index theory | Unlikely — index-theoretic direction, not convergence |
| van den Dungen (Bonn) | Dirac–Schrödinger index | No — pivoted out of metric direction |
| Nieuviarts | Twist morphism (even-dim) | No — odd-dim restriction acknowledged |
| Martinetti (Genoa) | Twisted SM Krein structure | No — focused on SM particle content, not convergence |
| Franco, Eckstein | Lorentzian causality from algebra | Possibly — but no convergence theorem in their program |

**Net assessment:** the Mondino–Sämann thread is the highest scoop-risk direction. Their work is moving toward measured Lorentzian GH (June 2025: Braun–Sämann) which is one step from "operator-system Lorentzian GH" if they pick up an NCG collaborator. However, their published trajectory has been entirely synthetic, with no spectral-triple content. The risk is real but moderate.

**Latrémolière scoop-risk = LOW.** His 2024–2025 trajectory is firmly Riemannian and his 2026 book consolidates Riemannian propinquity. No signal that a Lorentzian propinquity volume is in preparation.

**Recommendation:** if Sprint L3 launches, pursue Entry Point A (Connes–vS × BBB Krein truncation) first. It is sprint-scale, has lowest scoop risk (no active researcher in the intersection), and would produce a publishable theorem within 1–3 months. Entry Points B and C can be sequenced after.

---

## §8. Net verdict

### Recommendation: SPRINT L3a (Entry Point A, Connes–vS × BBB Krein, 1–3 months) is justified by the literature.

**Why L3a not full L3:**

1. **Entry Point A has a specific scaffold.** The Connes–vS 2021 truncation framework + BBB 2018 Krein classification + Sprint L2-D bit-exact verified axioms at (m, n) = (4, 6) are all on the shelf. Writing the truncated-operator-system theorem on the Krein triple is sprint-scale.

2. **Entry Points B and C have less scaffold.** The full Lorentzian propinquity construction (Entry Point C) is genuinely original NCG-mathematics work. The May-2026 Sprint L2-A audit (6–12 months) estimate was correct on this part. The 2024–2025 synthetic Lorentzian GH literature on the metric-space side is mature, but the bridge to operator systems has not been written.

3. **The Nieuviarts shortcut is structurally unavailable to GeoVac.** Two years on, the odd-dim restriction holds. Sprint L0's "may shorten Sprint L2 budget" note is now confirmed-false for GeoVac.

4. **Scoop risk is moderate-low.** Mondino–Sämann is the only direction at moderate scoop risk, but their work is entirely synthetic.

**Practical sequencing:**

- **Default next step (if L3 is selected):** Sprint L3a — Connes–vS truncated operator system on BBB (4, 6) Krein triple. 1–3 months. Sprint-scale, publishable.
- **Optional (after L3a lands):** Sprint L3b — Krein-positive-cone Wasserstein-Kantorovich + Berezin lift. 3–6 months.
- **Multi-month commitment, defer until L3a + L3b land:** Sprint L3c — bridge synthetic Lorentzian GH (Mondino–Sämann) to operator-system. 6–12 months. This would close Paper 42 §10 O1 in the GH limit.

**If Sprint L3 is NOT selected:** the alternatives (Paper 43 drafting, state-side dictionary expansion, multi-focal precision catalogue) remain valid per Sprint L2 §6 next-direction options. Sprint L0's note that L3 "may be skippable" is honest — the finite-cutoff Lorentzian closure (Sprint L2) is sufficient for the four-witness Wick-rotation theorem.

**Recommended next action:** PI selects between (A) draft Paper 43 first to consolidate Sprint L2's deliverables, or (B) launch Sprint L3a immediately to extend the Lorentzian arc into a new theorem before Paper 43 is finalized. The two are not mutually exclusive — L3a's deliverable would fit naturally into Paper 43 §11 or a Paper 44.

---

## §9. Files cross-referenced

- `CLAUDE.md` §1.7 WH1 entry (Sprint L2 closure paragraph)
- `papers/standalone/paper_42_modular_hamiltonian_four_witness.tex` §3.2 (Nieuviarts NO-GO footnote), §5.2 (Nieuviarts twist morphism subsection)
- `debug/lorentzian_l0_audit_memo.md` §4 (sequencing recommendation, L3 budget)
- `debug/sprint_l2_synthesis_memo.md` §5.1 (Sprint L3 named as open frontier)
- `debug/lorentzian_partition_transfer_memo.md` (cited L3 input)

End of memo.

---

### Sources (arXiv identifiers verified during this audit)

**Propinquity / metric spectral triples (Riemannian only):**
- arXiv:1302.4058 — Latrémolière 2013, Quantum GH Propinquity
- arXiv:1506.04341 — Latrémolière 2015, Propinquity book preprint
- arXiv:1811.10843 — Latrémolière 2018, GH propinquity for metric spectral triples (Adv. Math. 2022)
- arXiv:2112.11000 — Latrémolière 2021, Continuity of Dirac spectrum for spectral propinquity
- arXiv:2301.00274 — Farsi–Latrémolière 2023, Convergence of inductive sequences
- arXiv:2302.09117 — 2023, Isometry groups and GH convergence
- arXiv:2403.16323 — 2024, Spectral triples on noncommutative solenoids
- arXiv:2404.00240 — Farsi–Latrémolière 2024 (v3 Oct 2025), Collapse in NCG
- arXiv:2504.11715 — Farsi–Latrémolière April 2025, Continuity on Riemannian path
- arXiv:2410.15454 — UCP maps GH convergence

**Lorentzian / pseudo-Riemannian spectral triples (Krein-based):**
- arXiv:1505.01939 — van den Dungen 2016, Krein spectral triples and fermionic action
- arXiv:1611.07062 — Bizi–Brouder–Besnard 2018, (m,n) classification
- arXiv:1611.07830 — Besnard 2016, Spacetimes in NCG Part I
- arXiv:1210.6575 — Strohmaier 2012/2014, Temporal Lorentzian spectral triples
- arXiv:1710.04965 — Devastato–Lizzi–Martinetti 2018, Lorentz signature & twisted spectral triples
- arXiv:1409.1480 — Franco–Eckstein 2014, Lorentzian structures & causality
- arXiv:1508.01917 — Franco–Eckstein 2015, Two roads to NC causality
- arXiv:2603.03216 — Martinetti March 2026, Twisted SM and Krein structure

**Nieuviarts twist-morphism (even-dim only):**
- arXiv:2402.05839 — v6 March 2026, Signature change by morphism
- arXiv:2502.18105 — v3 May 2025, Emergence of Lorentz symmetry
- arXiv:2512.15450 — v2 May 2026, Emergence of Time (proceedings synthesis)

**Synthetic Lorentzian GH convergence (no spectral triples):**
- arXiv:2209.14384 — Lorentzian metric spaces & GH convergence
- arXiv:2209.12736 — GH metrics & dimensions of Lorentzian length spaces
- arXiv:2412.04311 — Bykov–Minguzzi–Suhr Dec 2024, unbounded case
- arXiv:2504.10380 — Mondino–Sämann April 2025, GH and pre-compactness
- arXiv:2506.10852 — Braun–Sämann June 2025, Reconstruction and measured GH
- arXiv:2510.24423 — Minguzzi Oct 2025, Results on Lorentzian metric spaces
- arXiv:2511.18389 — Nov 2025, Intrinsic Timed Hausdorff Convergence

**Lorentzian index theory (Bär–Strohmaier lineage):**
- arXiv:2602.16547 — 2026, Lorentzian Equivariant Index Theorem
- J. London Math. Soc. 112:e70301 (2025) — van den Dungen, Dirac-Schrödinger and spectral flow
- Forum of Math. Sigma 13:e11 (2025) — van den Dungen, Callias Theorem
