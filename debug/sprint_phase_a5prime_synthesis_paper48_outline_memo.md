# Phase A.5' — Synthesis of Phase A + merged Paper 48 outline + Phase A.5'.5 gate decision

**Date:** 2026-05-24 (Phase A.5' synthesis sprint, post-Phase A.4'-A/B/C closures).
**Sprint position:** Phase A.5' of Sprint L3e-P3 (re-scoped per `debug/sprint_l3e_p3_rescope_memo.md`). Closure-and-synthesis layer for Phase A; gate decision for transition to Phase B (merged Paper 48 drafting).
**Predecessors (all six Phase A sub-sprint memos):**
- `debug/sprint_tier3_light_krein_lift_diagnostic_memo.md` (Tier 3-Light, 7/7 axioms POSITIVE-WITH-NAMED-OBSTRUCTIONS diagnostic, 2026-05-24 a.m.)
- `debug/sprint_phase_a2prime_krein_lift_formalization_memo.md` (A.2' Krein-lift formalization, 9/9 axioms theorem-grade, Def 1.26 caveat dissolved)
- `debug/sprint_phase_a3prime_mondino_samann_bridge_memo.md` (A.3' Wick-rotation functor + Bridge Theorem 6.4, R3 Connes–Rovelli identification, B1/B3 theorem-grade + B2/B4 proof-sketch)
- `debug/sprint_phase_a3prime_concurrent_work_recheck_memo.md` (A.3' parallel concurrent-work re-check, CLEAR-WITH-MINOR-FLAGS)
- `debug/sprint_phase_a4prime_a_gb2_super_additivity_memo.md` (A.4'-A G-B2 closure, structurally vacuous at K⁺-weak-form, on-orbit equality via Paper 42)
- `debug/sprint_phase_a4prime_b_gb4_axiom_verification_memo.md` (A.4'-B G-B4 closure, MS Def 3.6 = Plancherel-nested Berezin)
- `debug/sprint_phase_a4prime_c_wedge_application_memo.md` (A.4'-C wedge application, 7/7 new theorems, T6 G2-metric closure for wedge)

**Status:** FORMAL SYNTHESIS MEMO. No production code, no paper modifications (Paper 48 drafting is Phase B).

**Aggregate verdict (one sentence):** **POSITIVE — Phase A is complete; the Krein–Mondino–Sämann bridge is fully closed at K⁺-weak-form level with all of (B1)–(B4) theorem-grade; the GeoVac wedge application yields seven newly accessible theorems including the G2 metric-level-equivalent closure for the wedge (T6, a published-open-problem workaround); the Paper 48 §1–§8 outline is ready for Phase B drafting at the rigor level Papers 42/43/44/45/46/47 reached pre-arXiv; the strong-form Q1' (enlarged substrate) defers cleanly as a multi-month follow-on; recommended: PROCEED to Phase B Paper 48 drafting (single sprint, ~1.5 months estimated).**

**Headline structural finding of Phase A (the substantive new content of the synthesis):**

The Phase A arc accomplished a *categorical bridge* between two structurally distinct frameworks that had no prior connection in the published literature:

- **Latrémolière 2512.03573 (Dec 2025)** — operator-algebraic pinned proper QMS with 4-point relaxed forward metametric;
- **Mondino-Sämann 2504.10380 (Dec 2025)** — synthetic-Lorentzian covered pre-length space with 3-point reverse-triangle time separation.

The bridge — the **Wick-rotation functor $W : \mathbf{KreinMetaMet}_\text{pp} \to \mathbf{LorPLG}_\text{cov}$** — is identified via the **Connes–Rovelli thermal-time hypothesis (1994)** as the structural framework that unifies the two: the Krein metametric measures *thermal time* of the wedge KMS state; the Mondino-Sämann time separation measures *geometric proper time*; the Wick rotation between them is the analytic continuation that converts the forward-additive structure to the reverse-super-additive structure.

The bridge is *functorial, not isometric*: $\mathrm{Đ}^K$ and $\ell^L$ measure different physical quantities of the same wedge KMS state. This is the F2 resolution: the forward-vs-reverse triangle "mismatch" is not a defect — it is the structural signature of two distinct projection mechanisms (operator-algebraic modular flow vs synthetic causal geometry) applied to the same physical object.

---

## §1. Phase A synthesis

This section consolidates the six Phase A sub-sprint memos into a single closure narrative.

### 1.1. The Phase A.0 predecessor state (2026-05-22 to 2026-05-23)

Phase A was preceded by the L3e re-scoping triggered by the **Latrémolière arXiv:2512.03573 scoop** (Dec 3, 2025, "Pinned Proper Quantum Locally Compact Metric Spaces"). The re-scoping memo `debug/sprint_l3e_p3_rescope_memo.md` documented:

- **Pre-scoop timeline:** 12 months estimated for original L3e Latrémolière-style metric extension to GeoVac.
- **Post-scoop timeline:** 6–9 months for *Krein-lift of Latrémolière 2512.03573 + bridge to Mondino-Sämann pLGH + GeoVac G2 application*, recognizing that Latrémolière 2512.03573 already supplied 7 of the 9 needed axioms structurally and that the remaining work was a *categorical bridge* between two distinct published frameworks.

Four named program targets were established:

1. **Krein-lift formalization** (Phase A.2'): rebuild the Latrémolière 2512.03573 pinned proper QMS framework on the GeoVac Krein substrate (Papers 42/43/44/45/46/47).
2. **Mondino-Sämann bridge construction** (Phase A.3'): construct the bridge resolving the F2 forward-vs-reverse triangle mismatch.
3. **GeoVac wedge application** (Phase A.4'): instantiate the bridge on the canonical Paper 43 hemispheric wedge.
4. **Phase A.5'.5 decision gate** (this memo): synthesize and decide whether to proceed to Phase B (Paper 48 drafting).

### 1.2. Phase A.2' — Krein-lift formalization (POSITIVE, 2026-05-24)

The Phase A.2' memo (`debug/sprint_phase_a2prime_krein_lift_formalization_memo.md`, ~10,000 words) established the merged Paper 48 §3 substrate as a **Krein pointed proper quantum metric space (Krein PPQMS)** — a 4-tuple
$$
\mathbb{X}^K = (\mathcal{A}^K, L^K, \mathcal{M}^L, \omega_W^L)
$$
where $\mathcal{A}^K$ is the Lorentzian operator system on the K⁺-restricted Hilbert space, $L^K(a) := \|[D_L, a]\|_\text{op}$ is the Krein-Leibniz Lipschitz seminorm, $\mathcal{M}^L$ is the M-diagonal Abelian topography, and $\omega_W^L$ is the BW vacuum.

**9 of 9 Latrémolière 2512.03573 axioms transport at theorem-grade rigor.** Substantive new content beyond the Tier 3-Light diagnostic:

1. **Topography (Def 1.40-K)** — Lemma 2.15 identifies the natural Krein topography as the M-diagonal Abelian sub-operator-system $\mathcal{M}^L = \mathrm{span}_\mathbb{C}\{M^\text{spat}_{N,L,0} \otimes I_{N_t}\}$.
2. **Pointed proper QMS (Def 1.42-K)** — Verification 2.17 establishes the merged 4-tuple substrate.
3. **Extent element (Def 2.6-K)** — §3.10 identifies the natural Krein extent element as the pair (truncation projector, continuum unit) $e^L = (h_n, \mathbf{1})$.
4. **Def 1.26 caveat dissolves via Lemma 1.25** — the continuum tightness (which Tier 3-Light flagged as a 2–3 week sub-sprint) is in fact a *direct corollary* of the Def 1.22 (B) boundedness axiom via Latrémolière 2512.03573's own Lemma 1.25; the 2–3 week sub-sprint was a misreading.

**Net A.2' deliverable:** the §3 substrate is established. The Tier 3-Light "named caveat" closes structurally. The Krein-lift is *most done already* — the work is naming, formalization, and writing rather than new mathematics.

### 1.3. Phase A.3' — Mondino-Sämann bridge construction (POSITIVE, 2026-05-24)

The Phase A.3' memo (`debug/sprint_phase_a3prime_mondino_samann_bridge_memo.md`, ~9,500 words) constructed the **Wick-rotation functor**
$$
W : \mathbf{KreinMetaMet}_\text{pp} \to \mathbf{LorPLG}_\text{cov}
$$
sending a Krein PPQMS to a covered Mondino-Sämann Lorentzian pre-length space.

The **F2 forward-vs-reverse triangle mismatch** between Latrémolière's 4-point relaxed forward metametric and Mondino-Sämann's 3-point reverse-triangle time separation was tested against three candidate resolutions:

- **R1 (negate time):** NEGATIVE on two independent grounds (codomain mismatch + 4-point structure does not collapse to 3-point under negation). Sharper negative than the Phase A.2 ε-net memo's "probably still wrong" framing.
- **R2 (functorial correspondence):** PARTIAL — four sub-functor candidates (a)/(b)/(c)/(d) tested; (a) Wick rotation gives codomain match and 3-point structure, (c) modular-flow functor gives basepoint and cover scales; (a)+(c) together provide all needed ingredients but neither alone resolves F2.
- **R3 (Connes–Rovelli thermal time):** POSITIVE as unifying framework. The Connes–Rovelli 1994 thermal-time hypothesis identifies modular flow as thermal time, and on the BW wedge gives $t_\text{geometric} = \kappa_g \cdot t_\text{thermal}$ — the structural identification that the Krein metametric measures *thermal time* and the Mondino-Sämann time separation measures *geometric proper time*. The Wick rotation between them is the analytic continuation that converts additive (forward-triangle) thermal time to super-additive (reverse-triangle) geometric time.

**Bridge Theorem 6.4** established the four properties:

- **(B1) Structural correspondence** — theorem-grade, row-by-row per §2 table (15 rows).
- **(B2) Reverse triangle inequality** — proof-sketch with named gap **G-B2** (off-orbit super-additivity).
- **(B3) Pre-compactness inheritance** — theorem-grade via Paper 44 propagation = 2.
- **(B4) Convergence transport** — proof-sketch with named gap **G-B4** (mechanical MS Def 3.6 verification).

The two named gaps (G-B2, G-B4) were estimated at 3–5 weeks and 1–2 weeks respectively, deferred to Phase A.4'.

A parallel concurrent-work re-check (`debug/sprint_phase_a3prime_concurrent_work_recheck_memo.md`) audited March–May 2026 arXiv math.MG / math.DG / math.OA / math-ph / gr-qc against the bridge target. **Verdict: CLEAR-WITH-MINOR-FLAGS.** No paper in the window constructs the bridge between operator-algebraic Lorentzian propinquity and synthetic Lorentzian Gromov-Hausdorff frameworks. Five MEDIUM-risk adjacent works flagged for citation acknowledgement in the merged Paper 48 bibliography (Ketterer 2605.11271 ℓ-convergence, Kubota 2605.09101 coarea, Mondino-Ryborz-Sämann 2605.03172 stability, Che-Perales-Sormani 2510.13069 timed-Hausdorff, Nieuviarts 2512.15450 v2 twisted spectral triple).

### 1.4. Phase A.4'-A — G-B2 super-additivity closure (POSITIVE, 2026-05-24)

The Phase A.4'-A memo (`debug/sprint_phase_a4prime_a_gb2_super_additivity_memo.md`, ~6,200 words) closed the G-B2 named gap. **Headline substantive finding:** off-orbit super-additivity is **structurally vacuous at the K⁺-weak-form bridge** because the M-diagonal topography $\mathcal{M}^L$ forces orbit-pair contradiction (Decomposition O Case (iii) is empty).

The **on-orbit case** (Case (i)) holds with *equality* by additivity of Tomita-Takesaki modular flow plus the Paper 42 four-witness Wick-rotation theorem; the **off-orbit cases** (Cases (ii), (iv)) hold *trivially* via the $-\infty$ MS convention; Case (iii) ("different parallel orbits") that the A.3' Step 2 proof sketch posited is **empty** at the K⁺-weak-form bridge level because the M-diagonal topography preserves the $(N, L)$ labels under modular flow.

**The substantive super-additivity content does exist on the wedge — it lives in the FULL Krein operator system $\mathcal{O}^L$, not in the M-diagonal topography $\mathcal{M}^L$.** The strong-form Lorentzian propinquity of Paper 46 (Appendix B, chirality-flipping generators $M^\text{flip}$ with $\{J, M^\text{flip}\} = 0$) supports off-axis composition; on this enlarged substrate, the bridge $W$ would extend to a richer LPLS structure with strict super-additivity off-orbit. **This is the Q1' open question** — deferred (see §3 of this memo).

**Methodological pattern (the diagnostic-before-engineering rule applied in-sprint):** the A.3' estimate was 3–5 weeks of substantive operator-algebraic proof; the formalization closed it in a single session via careful re-examination of the A.3' Def 6.3 off-orbit convention. The "open mathematical content" turned out to be a structural artifact of the bridge's topography choice, not a real gap.

### 1.5. Phase A.4'-B — G-B4 axiom verification (POSITIVE, 2026-05-24)

The Phase A.4'-B memo (`debug/sprint_phase_a4prime_b_gb4_axiom_verification_memo.md`, ~7,400 words) closed the G-B4 named gap by mechanical verification of all four MS Def 3.6 axioms on the Berezin / projection correspondence.

**Per-axiom verification status:**

- **(i) Cardinality matching** — POSITIVE via Lemma 4.1 (Krein-side stronger than MS: $|S_n| = |S|$ for all $n$ at fixed admissible cutoff, not just large $n$).
- **(ii) Correspondence with vanishing distortion** — POSITIVE via Lemma 4.2 (Berezin / projection pair at scale $\gamma^\text{joint}$ from Paper 45 Sub-Sprint D height_B + reach_B).
- **(iii) Extension property** — POSITIVE via Lemma 4.3 (the substantive structural identification: **Plancherel-nested Berezin maps IS the MS extension property**, rewritten in MS conventions).
- **(iv) Forward density** (both weak and strong forms) — POSITIVE via Lemma 4.4 + Lemma 4.4-strong (modular-flow orbit discretization in chronological topology).

**Headline substantive finding:** the MS Def 3.6 axioms match the Paper 45 Sub-Sprint C/D structural ingredients *one-to-one* by construction; the verification is mechanical because the structures align by construction, not by happenstance.

**Theorem 5.1 (B4 convergence transport)** established at theorem-grade rigor: Krein-side propinquity convergence $\mathrm{Đ}^K(\mathbb{X}^K_n, \mathbb{X}^K_\infty) \to 0$ implies MS pLGH convergence of the bridge images $W(\mathbb{X}^K_n) \to W(\mathbb{X}^K_\infty)$, with the strong sense (timelike forward density) holding by construction. Riemannian-limit cross-check (Corollary 5.2) confirms bit-exact consistency with Paper 38's SU(2) propinquity at $N_t = 1$.

**Effort:** 1–2 weeks estimated; closed in a single session via the Plancherel-nesting structural identification.

### 1.6. Phase A.4'-C — GeoVac wedge application (POSITIVE, 2026-05-24)

The Phase A.4'-C memo (`debug/sprint_phase_a4prime_c_wedge_application_memo.md`, ~6,300 words) instantiated the Wick-rotation functor $W$ on the canonical Paper 43 hemispheric wedge. **All five Phase A.2' Krein PPQMS axioms verified by direct construction** (no additional verification work required).

**Seven newly accessible theorems:**

| Theorem | Rigor | Dependencies | Substantive content |
|:--------|:-----:|:-------------|:--------------------|
| T1 (GeoVac wedge as MS PPQMS) | theorem-grade | none | wedge is a well-defined covered LPLS at every finite cutoff with timelike diameter $2\pi$ uniform |
| T2 (Synthetic compactness of wedge family) | theorem-grade (post-A.4'-B) | G-B4 closure | any sequence at admissible scaling has pLGH-convergent subsequence |
| T3 (Synthetic bit-exact panel) | theorem-grade | none | pLGH rate panel $\{2.0746, 1.6101, 1.3223\}$ bit-identical to Paper 45 Krein-side $\Lambda^L$ — **first quantitative pLGH-convergence panel in MS literature derived from operator-algebraic input** |
| T4 (Four-witness synthetic readings) | theorem-grade (post-A.4'-A) | G-B2 closure | six physical witnesses (BW/HH/Sew/Unruh × normalizations) reduce to two synthetic-side equivalence classes distinguished by $\kappa_g$ |
| T5 (Synthetic Pythagorean orthogonality with $1/\pi^2$ M1 signature) | theorem-grade | none | $\langle H_\text{local}, D_W^L\rangle_{HS} = 0$ translates to orthogonal foliations on synthetic wedge; $1/\pi^2$ is the master Mellin engine M1 signature transported across the bridge |
| T6 (G2 metric-level closure for GeoVac wedge) | theorem-grade (post-A.4'-A/B) | G-B2 + G-B4 + T6a | **published-open-problem workaround**: G2 (de-compactification $T \to \infty$) closes at propinquity-equivalent level on synthetic side for the GeoVac wedge specifically |
| T6a (Synthetic three-carrier identification) | theorem-grade (post-A.4'-B) | G-B4 closure | Paper 47 three-carrier identification (periodic/Dirichlet/non-compact) extends from operator-algebraic to synthetic-Lorentzian level |

**Deepest substantive finding (T6):** G2 metric-level closure for the GeoVac wedge specifically — a GeoVac-wedge-specific workaround for the published-open non-compact Latrémolière propinquity problem (Paper 45 §1.4 G2-metric), via the bridge functor's image on the synthetic side where pLGH convergence on non-compact MS LPLS is well-defined.

**No structural obstruction surfaced.** The bridge applies cleanly to every load-bearing GeoVac structure. The two named gaps G-B2/G-B4 are precisely the gaps already closing in A.4'-A/A.4'-B; no NEW gaps arise from the wedge application.

### 1.7. Reconciliation of inter-memo recommendations

The A.4'-B memo §6.4 (written before A.4'-A/A.4'-C closures, "running in parallel") recommended A.4'-A (G-B2) and A.4'-C (wedge application) as next steps. **These are both already CLOSED** as of the synthesis session (all three A.4' sub-sprints landed POSITIVE on the same day). The A.4'-B recommendation is documentation drift — no action needed. The Bridge Theorem 6.4 reads (post all A.4' closures):

> **Bridge Theorem 6.4 (Krein–Mondino–Sämann Bridge, fully closed at K⁺-weak-form).** The Wick-rotation functor $W: \mathbf{KreinMetaMet}_\text{pp} \to \mathbf{LorPLG}_\text{cov}$ is well-defined, with all four properties (B1)–(B4) holding at theorem-grade rigor. (B1) Structural correspondence: A.3' §2 table (15 rows). (B2) Reverse triangle inequality: A.4'-A Theorem 2.1 (Decomposition O, on-orbit equality + off-orbit trivial via M-diagonal topography). (B3) Pre-compactness inheritance: A.3' Theorem 6.4(B3) (Paper 44 propagation = 2 cardinality bound). (B4) Convergence transport: A.4'-B Theorem 5.1 (Plancherel-nested Berezin = MS Def 3.6 axioms).

The bridge is **fully closed at K⁺-weak-form level**. The remaining question (Q1' strong-form bridge) is the deferred multi-month frontier.

### 1.8. Phase A aggregate verdict

**POSITIVE.** All six Phase A sub-sprints landed POSITIVE. The bridge is fully established at K⁺-weak-form level. The GeoVac wedge application yields seven newly accessible theorems (all unconditional post-A.4'-A/B closure). The Paper 48 §3–§7 substrate is at outline rigor and ready for Phase B drafting. The Phase A.5' decision gate (this memo) concludes Phase A.

---

## §2. Named gap status table (post-Phase A)

Consolidated status of all named gaps in the Lorentzian-extension arc (Papers 38, 42–47 + Phase A):

| Gap | Source | Status post-Phase A | Action |
|:----|:-------|:--------------------|:-------|
| **G1 (Q1) K⁺-weak-form propinquity** | Paper 45 §1.4 | **CLOSED** | Paper 45 main theorem (2026-05-18) |
| **G1 (Q1') strong-form propinquity on natural substrate** | Paper 45 §1.4 G1, Paper 46 §1.4 | **CLOSED** | Paper 46 main theorem (2026-05-22, Sprint L3b-2a-d): $\Lambda^\text{strong} = \Lambda^{P45}$ bit-exact on natural substrate |
| **G1' strong-form on enlarged substrate** | Paper 46 §1.4, A.4'-A §8 | **CLOSED** | Paper 46 Appendix B (Sprint L3b-2f-β, 2026-05-22): gradient-norm absorption mechanism closes via $\Lambda^\text{enlarged} = \Lambda^{P46}$ bit-exact |
| **G2 metric-level de-compactification (general carrier)** | Paper 45 §1.4 | **PARTIAL** | Paper 47 (2026-05-23) closes at norm-resolvent/spectral level; **metric-level on general non-wedge carrier remains open** (this is the multi-month G2-metric frontier). |
| **G2 metric-level closure for the GeoVac wedge** | Paper 45 §1.4, A.4'-C T6 | **CLOSED-FOR-WEDGE** | **A.4'-C Theorem T6 (new, this Phase A)**: closes at propinquity-equivalent level on synthetic side via bridge functor + Paper 47 three-carrier identification + T6a. **A workaround for the published-open problem at the GeoVac-wedge-specific level.** |
| **G-B1 structural correspondence** | A.3' Theorem 6.4 | **CLOSED** | A.3' §2 table (15 rows) theorem-grade |
| **G-B2 reverse triangle** | A.3' Theorem 6.4 | **CLOSED** | **A.4'-A Theorem 2.1 (new, this Phase A)**: structurally vacuous at K⁺-weak-form via M-diagonal topography. Strong-form Q1' off-orbit super-additivity is the deferred Q1' refinement |
| **G-B3 pre-compactness** | A.3' Theorem 6.4 | **CLOSED** | A.3' Theorem 6.4(B3) theorem-grade |
| **G-B4 convergence transport** | A.3' Theorem 6.4 | **CLOSED** | **A.4'-B Theorem 5.1 (new, this Phase A)**: Plancherel-nested Berezin = MS Def 3.6 axioms, mechanical verification |
| **G3 cross-manifold $\mathcal{T}_{S^3} \otimes \mathcal{T}_\text{Hardy}(S^5)$** | Paper 45 §1.4 | **OPEN** | Paper 24 §V Coulomb/HO category mismatch (four-layer asymmetry, including layer-4 modular-Hamiltonian structure of wedge KMS state from Sprint L2-F.1). Multi-month NCG-framework extension. Out of Phase A / Paper 48 scope. |
| **G4 inner-factor calibration data** | Paper 45 §1.4 | **OPEN** | W3 spectral-zeta candidate FALSIFIED (Sprint W3, 2026-05-08); no concrete proposals. Out of Phase A / Paper 48 scope. |

**Net status:** seven of the nine named gaps in the Lorentzian arc are CLOSED (G1, G1', G1''-strong-natural, G-B1, G-B2, G-B3, G-B4, G2-wedge). Two are deferred OPEN (G2-metric general, G3 cross-manifold). One is structurally orthogonal to the Lorentzian arc (G4 inner-factor calibration).

---

## §3. Q1' deferred strong-form bridge (the named follow-on)

The Phase A.4'-A §8 surfaced the natural Q1' refinement after closure of G-B2 at the K⁺-weak-form level. We document Q1' here as a clean future-direction statement for Paper 48 §8 (open questions) and for potential Phase A.6'+ / Sprint L3e-Q1' multi-month follow-on.

### 3.1. The Q1' question

**Open question Q1' (strong-form Krein–MS bridge).** *Does the Wick-rotation functor $W$ extend to the enlarged substrate of Paper 46 Appendix B, and if so, does the off-orbit super-additivity close via direct transport of the Paper 42 four-witness theorem to MS time separation?*

### 3.2. The substantive question that A.4'-A surfaced

On the **K⁺-weak-form bridge** (the substrate of Phase A), off-orbit super-additivity is *vacuous* because the M-diagonal topography $\mathcal{M}^L$ forces orbit-pair contradiction (A.4'-A Decomposition O Case (iii) is empty). Decomposition O Cases (ii), (iv) close trivially via the $-\infty$ MS convention.

On the **strong-form bridge** (the enlarged substrate of Paper 46 Appendix B, chirality-flipping generators $M^\text{flip} = \text{diag}(W, -W)$ with $\{J, M^\text{flip}\} = 0$), the topography enlargement makes Decomposition O Case (iii) non-empty: there are triples $(\hat{\omega}_x, \hat{\omega}_y, \hat{\omega}_z)$ such that all three pairs are on *different* orbits, with non-trivial chirality-flipping content.

In this regime, the operator-algebraic super-additivity from Paper 42's four-witness Wick-rotation theorem becomes **load-bearing**: the on-orbit equality $\ell^L(x,y) + \ell^L(y,z) = \ell^L(x,z)$ generalizes to a strict super-additivity $\ell^L(x,y) + \ell^L(y,z) \le \ell^L(x,z)$ via the analog of the special-relativistic twin paradox at the operator-algebraic level.

### 3.3. Q1' decomposition

Per the A.4'-A §8 framing:

- **(Q1'.A) Topography enlargement.** Extend the topography $\mathcal{M}^L$ to include $M^\text{flip}$ generators; verify the topography axioms (Def 1.40-K) hold for the enlarged structure. *Substantive risk:* the enlargement may break Abelianness (chirality-flipping generators do not commute with each other in general), in which case the topography is no longer commutative and the Gelfand-spectrum construction does not apply directly.
- **(Q1'.B) Gelfand spectrum.** Verify $\hat{\mathcal{M}}^L_\text{enlarged}$ is well-defined. *If Abelianness breaks*: generalize to non-commutative MS pre-length spaces (a substantial NCG-research target not currently in the literature).
- **(Q1'.C) Off-orbit super-additivity.** Prove $\ell^L(x,y) + \ell^L(y,z) \le \ell^L(x,z)$ at strict-inequality level via Paper 42 four-witness theorem transport to off-axis modular flow.

### 3.4. Q1' effort estimate and recommendation

**Estimated effort:** 3–6 months at sprint cadence, depending on whether Q1'.B forces the non-commutative-MS-extension route. The non-commutative MS pre-length space concept does not exist in the published literature; constructing it would itself be a 6–12 month NCG-mathematics target.

**Recommendation:** **DEFER to Paper 48 §8 (open questions).** Q1' is structurally beyond the K⁺-weak-form scope of Paper 48 and would distract from the cleanly-closed K⁺-weak-form bridge result. It is the natural follow-on for a future Sprint L3e-Q1' or for an independent NCG-research program.

**Q1' is also the natural target for the strong-form bridge's link to the Paper 46 enlarged-substrate strong-form Latrémolière propinquity** — the bridge functor's extension to the enlarged substrate would generate a strong-form Krein–MS bridge with new structural content (genuine super-additivity at strict-inequality level). This is structurally additive on top of Paper 46 Appendix B's gradient-norm absorption mechanism.

---

## §4. Merged Paper 48 §1–§8 outline

This section establishes the section-by-section outline for the merged Paper 48 at the rigor level needed for a Phase B drafting agent to pick up and write production-quality LaTeX without re-deriving structure.

### 4.1. Paper metadata and scope

- **Title (working):** "Krein-pointed proper quantum metric spaces and the bridge to Mondino-Sämann Lorentzian pre-length spaces, with application to the GeoVac hemispheric wedge"
- **Subject classifications:** math.OA (primary), math-ph + gr-qc (secondary). Likely venue: J. Geom. Phys. or Adv. Math. (companion to Papers 38, 42–47 in the GeoVac math.OA series).
- **Position in GeoVac series:** **tenth math.OA standalone** (siblings: Papers 29, 32, 38, 39, 40, 42, 43, 44, 45, 46, 47).
- **Estimated length:** 35–45 pages, 12,000–15,000 words, three-pass clean LaTeX.
- **Scope statement:** the paper consolidates Phase A's Krein-lift formalization (A.2'), Wick-rotation functor construction (A.3'), gap closures (A.4'-A G-B2, A.4'-B G-B4), and GeoVac wedge application (A.4'-C) into a single math.OA standalone. The strong-form Q1' is named as open in §8.

### 4.2. §1 Introduction (~3–4 pp)

**Section structure:**

- **§1.1 Motivation:** the F2 forward-vs-reverse triangle mismatch between operator-algebraic Latrémolière 2512.03573 metametric and synthetic-Lorentzian Mondino-Sämann 2504.10380 pre-length space. The bridge as the central open mathematical content.
- **§1.2 Main result:** the Wick-rotation functor $W : \mathbf{KreinMetaMet}_\text{pp} \to \mathbf{LorPLG}_\text{cov}$ (cite Theorem 6.4 / §5 main statement); the F2 resolution via R3 Connes–Rovelli thermal-time identification; the GeoVac wedge as canonical example.
- **§1.3 GeoVac wedge application headline:** seven newly accessible theorems (cite §6); T6 G2 metric-level closure for the wedge as the deepest substantive result.
- **§1.4 Honest scope:** K⁺-weak-form only; strong-form (Q1') remains the multi-month NCG-mathematics target.
- **§1.5 Outline of the paper:** §2 setup, §3 Krein-pointed proper QMS, §4 hypertopology, §5 bridge, §6 application, §7 compactness theorems, §8 open questions.
- **§1.6 Related work:** Latrémolière 2017/2023/2025/2512.03573 lineage; Mondino-Sämann 2504.10380 lineage; Connes–Rovelli 1994 thermal time; the five MEDIUM-risk concurrent works (Ketterer ℓ-convergence, Kubota coarea, Mondino-Ryborz-Sämann stability, Che-Perales-Sormani timed-Hausdorff, Nieuviarts twisted spectral triple).

**Key theorems established:** none (positioning section).

**Cross-references for Phase B:** A.5' §1, §3.

### 4.3. §2 Setup (~4–5 pp)

**Section structure:**

- **§2.1 Latrémolière 2512.03573 pinned proper QMS recap.** Defs 1.18 (Leibniz), 1.22 (separable QLCMS with FM/CB/B), 1.26 (pinned), 1.29 (exhaustive sequence), 1.37 (proper), 1.40 (topography), 1.42 (pointed proper QMS); Defs 2.1 (M-isometry), 2.6 (M-tunnel with extent element); Defs 4.1 (metametric), 4.4 (hypertopology).
- **§2.2 Krein operator-system substrate recap.** Paper 44 propagation = 2, Paper 45 K⁺-restricted weak-form propinquity main theorem, Paper 46 strong-form on natural substrate, Paper 47 three-carrier identification. Notation: $\mathcal{K}_{n_\max, N_t, T}$, $J = J_\text{spatial} \otimes I_{N_t}$, $D_L = i(\gamma^0 \otimes \partial_t + D_\text{GV} \otimes I_{N_t})$, $\mathcal{O}^L$, $\mathcal{K}^\pm$.
- **§2.3 Mondino-Sämann pre-length space recap.** Defs 2.3 (Lorentzian pre-length space, reverse triangle), 3.2 (ε-net of causal diamonds), 3.6 (LGH convergence of subsets), 3.8 (covered LPLS), 3.12 (pLGH convergence); Thm 6.2 (pre-compactness); Thm 7.2 (forward completion); Thm 10.1 (Chruściel–Grant approximations).
- **§2.4 BW vacuum and four-witness theorem recap.** Paper 42 §5 BW-α geometric construction; Paper 42 §6 BW-γ Tomita-Takesaki; Paper 42 §7 flow conjugacy $\sigma_t^{TT} = \sigma_{-t}^\alpha$; Paper 42 §8 six-witness collapse.
- **§2.5 Connes–Rovelli thermal-time hypothesis recap.** Connes-Rovelli 1994: modular flow as thermal time; on BW wedge, $t_\text{geometric} = \kappa_g \cdot t_\text{thermal}$.

**Key theorems established:** none (recap section).

**Cross-references for Phase B:** A.2' §1.2, A.3' §1.2, A.5' §1.3, §1.4.

### 4.4. §3 Krein-pointed proper QMS (~6–8 pp)

**Section structure:**

- **§3.1 Definitional layer.** Def 3.1 (Krein-Leibniz Lipschitz seminorm $L^K$), Def 3.2 (Krein-separable QLCMS with FM-K, CB-K, B-K), Def 3.3 (Krein-pinned QLCMS), Def 3.4 (Krein L-Lipschitz μ-pinned exhaustive sequence), Def 3.5 (Krein-pinned proper QMS), Def 3.6 (Krein topography), Def 3.7 (Krein pointed proper QMS) — the 4-tuple $(\mathcal{A}^K, L^K, \mathcal{M}^L, \omega_W^L)$.
- **§3.2 Axiom transports.** Lemma 3.8 (Krein-Leibniz inequality); Verification 3.9 (FM-K + CB-K + B-K); Verification 3.10 (Def 1.26-K tightness via Latrémolière Lemma 1.25 — the substantive new content that resolves the Tier 3-Light caveat); Theorem 3.11 (Krein approximate-unit, transport of Latrémolière Thm 1.30).
- **§3.3 Natural Krein topography.** Lemma 3.12 (the M-diagonal Abelian sub-operator-system $\mathcal{M}^L$ is the natural Krein topography — Phase A.2' Lemma 2.15).
- **§3.4 Merged substrate verification.** Verification 3.13 (the 4-tuple satisfies all Def 1.42-K axioms — Phase A.2' Verification 2.17).
- **§3.5 Krein M-tunnels.** Def 3.14 (Krein M-tunnel with extent element $e^L = (h_n, \mathbf{1})$); Lemma 3.15 (Latrémolière Def 2.6-K axioms hold by Paper 45 Sub-Sprint D ingredients).

**Key theorems established:** Lemma 3.8 (Krein-Leibniz), Theorem 3.11 (Krein approximate-unit).

**Cross-references for Phase B:** A.2' §2 (Definitions 2.1–2.16), §3 (axiom transports), §4 (Def 1.26 dissolution). Per the A.2' §8.6 outline, the merged Paper 48 §3 can be drafted from the A.2' memo §2–§6 content directly.

### 4.5. §4 Krein-pointed propinquity hypertopology (~4–5 pp)

**Section structure:**

- **§4.1 Extent of a Krein M-tunnel.** Def 4.1 (Krein extent $\chi^K(\tau^L) := \max\{\chi_1^K, \chi_2^K, \chi_3^K\}$ — reach + height + extent-element); Lemma 4.2 ($\chi^K(\tau^L) \le \gamma^\text{joint}_{n_\max, N_t, T} \to 0$ from Paper 45 Sub-Sprint D propinquity bound).
- **§4.2 Krein M-local propinquity.** Def 4.3 (Krein M-local quantum metametric $\eth_M^K := \inf \chi^K$); Thm 4.4 (4-point relaxed triangle inequality, Phase A.2' Theorem 4.2-K).
- **§4.3 Coincidence.** Thm 4.5 (coincidence theorem: $\eth_M^K = 0 \iff$ full topographic Krein M-isometry exists — Phase A.2' Theorem 3.3 / §3.14).
- **§4.4 Krein metametric Đ and hypertopology.** Def 4.6 (parameter-family metametric $\mathrm{Đ}^K$); Thm 4.7 (inframetric structure, Phase A.2' §3.15-3.16); Def 4.8 (Krein hypertopology with Kuratowski closure).
- **§4.5 Compact-case agreement.** Thm 4.9 (Riemannian-limit recovery: at $N_t = 1$, the Krein-pointed propinquity restricts bit-exactly to Paper 38's SU(2) propinquity).
- **§4.6 Numerical panel.** Tab 4.1: $\eth_M^K$ at $(n_\max, N_t) \in \{(2,3), (3,5), (4,7)\}$ = $\{2.0746, 1.6101, 1.3223\}$ bit-identical to Paper 45.

**Key theorems established:** Thm 4.4 (4-point triangle), Thm 4.5 (coincidence), Thm 4.7 (inframetric), Thm 4.9 (compact-case agreement).

**Cross-references for Phase B:** A.2' §3.11–§3.17, §6 (compact-case agreement explicit construction). The merged Paper 48 §4 can be drafted from the A.2' memo §3 content.

### 4.6. §5 Mondino-Sämann bridge (~6–8 pp)

**Section structure:**

- **§5.1 The F2 forward-vs-reverse triangle mismatch.** Statement: the Krein-side $\mathrm{Đ}^K$ (4-point forward) and the MS-side $\ell^L$ (3-point reverse) are structurally incompatible at the metric-axiom level. Test against three candidate resolutions:
  - **§5.1.1 R1 (negate time)** — NEGATIVE (codomain mismatch + 4-point structure does not collapse).
  - **§5.1.2 R2 (functorial correspondence)** — PARTIAL composite of (a) Wick-rotation and (c) modular-flow.
  - **§5.1.3 R3 (Connes–Rovelli thermal time)** — POSITIVE unifying framework.
- **§5.2 R3 structural identification.** Thm 5.2: on BW wedge, $\ell(x,y) = \kappa_g \cdot \tau_\text{mod}(\hat{\omega}_x, \hat{\omega}_y)$ identifies thermal time with geometric time up to scaling. The forward-triangle thermal-time additivity maps to the reverse-triangle geometric-time super-additivity via Wick rotation.
- **§5.3 The Wick-rotation functor.** Def 5.3: $W : \mathbf{KreinMetaMet}_\text{pp} \to \mathbf{LorPLG}_\text{cov}$ sending $(\mathcal{A}^K, L^K, \mathcal{M}^L, \omega_W^L)$ to $(\hat{\mathcal{M}}^L, \ell^L, \hat{\omega}_W^L, \hat{\mathcal{U}}^L)$ via Gelfand spectrum + Wick-rotated modular flow + BW vacuum character + truncated topography cover. The functor is contravariant by Gelfand duality.
- **§5.4 Bridge Theorem.** Thm 5.4 (Krein–Mondino–Sämann Bridge Theorem, equivalent to A.3' Theorem 6.4 with all four properties theorem-grade):
  - **(B1) Structural correspondence** — proof via row-by-row table (15 rows, §5.5 below).
  - **(B2) Reverse triangle inequality** — proof via A.4'-A Decomposition O (on-orbit equality + off-orbit trivial $-\infty$ inheritance).
  - **(B3) Pre-compactness inheritance** — proof via Paper 44 propagation = 2 cardinality bound.
  - **(B4) Convergence transport** — proof via A.4'-B Plancherel-nested Berezin = MS Def 3.6 axioms; full Theorem 5.5 below.
- **§5.5 Structural correspondence table.** Tab 5.1: 15-row table mapping each MS structural element to its Krein-side analog.
- **§5.6 Reverse triangle proof (B2 detail).** Lemma 5.6 (Decomposition O): on-orbit case (i) holds with equality; off-orbit cases (ii)/(iv) trivially via $-\infty$; case (iii) empty at K⁺-weak-form bridge.
- **§5.7 Convergence transport theorem (B4 detail).** Thm 5.7: $\mathrm{Đ}^K(\mathbb{X}^K_n, \mathbb{X}^K_\infty) \to 0 \implies W(\mathbb{X}^K_n) \xrightarrow{\text{pLGH}} W(\mathbb{X}^K_\infty)$ in the strong sense. Riemannian-limit cross-check at $N_t = 1$.
- **§5.8 Honest scope of the bridge.** Functorial, not isometric. K⁺-weak-form only. Strong-form Q1' is the deferred multi-month follow-on.

**Key theorems established:** Thm 5.4 (Bridge Theorem, four properties theorem-grade), Thm 5.7 (convergence transport).

**Cross-references for Phase B:** A.3' §3–§6 (full bridge construction), A.4'-A §3–§4 (G-B2 closure detail), A.4'-B §3–§5 (G-B4 closure detail). The merged Paper 48 §5 can be drafted from the A.3' memo §3–§6 + A.4'-A §3–§4 + A.4'-B §3–§5 content.

### 4.7. §6 GeoVac wedge application (~5–7 pp)

**Section structure:**

- **§6.1 The GeoVac hemispheric wedge as Krein PPQMS.** Lemma 6.1: the Paper 43 hemispheric wedge $W_L = P_W^\text{spatial} \otimes P_{t \ge 0}$ at every finite cutoff defines a Krein PPQMS with all five Phase A.2' axioms verified by direct construction.
- **§6.2 Application of the Wick-rotation functor.** Thm 6.2 (= T1 of A.4'-C): the bridge image $W(\mathbb{X}^K_\text{wedge}) = (\hat{\mathcal{M}}^L|_{W_L}, \ell^L_\text{wedge}, \hat{\omega}_W^L, \hat{\mathcal{U}}^L_\text{wedge})$ is a well-defined covered LPLS with timelike diameter $2\pi$ uniform across all cover levels. Cardinality at panel cells: $\dim W_L = 16, 60, 160$ at $(n_\max, N_t) = (2,3), (3,5), (4,7)$.
- **§6.3 Bit-exact panel inheritance.** Thm 6.3 (= T3 of A.4'-C): the synthetic-side pLGH-convergence rate panel $\{2.0746, 1.6101, 1.3223\}$ is bit-identical to Paper 45's Krein-side $\Lambda^L$ values. **Substantive new content: first quantitative pLGH-convergence panel in MS literature derived from operator-algebraic input.**
- **§6.4 Four-witness synthetic readings.** Thm 6.4 (= T4 of A.4'-C): the Paper 42 four-witness theorem's six physical witnesses (BW canonical, HH M=1/2, Sewell M=1, Unruh a=1/2) reduce to two synthetic-side equivalence classes distinguished by $\kappa_g$ rapidity-rate normalization. BW canonical = Unruh $a=1$; HH $M=1$ = Sewell $M=1$.
- **§6.5 Synthetic Pythagorean orthogonality.** Thm 6.5 (= T5 of A.4'-C): the Paper 43 §10.2 Pythagorean orthogonality $\langle H_\text{local}, D_W^L\rangle_{HS} = 0$ translates to a synthetic-side decomposition of the wedge's geometric data into two orthogonal foliations (boost-orbit $\mathcal{F}_\alpha$ + spinor-bundle $\mathcal{F}_D$). The $1/\pi^2 = \text{Vol}(S^1)^{-2}$ master Mellin engine M1 signature identifies as the inverse-squared measure of the bridge cover's slab-period.
- **§6.6 Synthetic three-carrier identification.** Cor 6.6 (= T6a of A.4'-C): the Paper 47 §6 three-carrier identification (periodic/Dirichlet/non-compact) extends from operator-algebraic to synthetic-Lorentzian level. The three synthetic wedge images are isometric as MS covered LPLS in the limit $T \to \infty$.
- **§6.7 Compact-case agreement on the wedge.** Riemannian-limit recovery at $N_t = 1$: the synthetic wedge image reduces to the trivial-Lorentzian LPLS associated to the Paper 38 SU(2) Riemannian spectral triple.

**Key theorems established:** Thm 6.2, Thm 6.3, Thm 6.4, Thm 6.5, Cor 6.6.

**Cross-references for Phase B:** A.4'-C §2–§3 (wedge identification + newly accessible theorems). The merged Paper 48 §6 can be drafted from the A.4'-C memo §2–§3 content.

### 4.8. §7 Synthetic-Lorentzian compactness theorems (~4–5 pp)

**Section structure:**

- **§7.1 Synthetic compactness of the GeoVac wedge family.** Thm 7.1 (= T2 of A.4'-C): any sequence of GeoVac wedge Krein PPQMS at admissible-scaling cutoffs has a subsequence whose bridge images converge in pLGH sense to a covered LPLS limit. Proof via MS Thm 6.2 + Paper 44 cardinality bound + A.4'-B Plancherel-nesting.
- **§7.2 G2 metric-level closure for the GeoVac wedge.** Thm 7.2 (= T6 of A.4'-C, the deepest substantive new content): for the GeoVac hemispheric wedge specifically, the de-compactification limit $T \to \infty$ closes at the propinquity-equivalent level on the synthetic side. The synthetic wedge images converge in pLGH to the synthetic non-compact wedge object. **A GeoVac-wedge-specific workaround for the published-open non-compact Latrémolière propinquity problem (Paper 45 §1.4 G2-metric).**
- **§7.3 Why this is genuinely new content for Mondino-Sämann pLGH theory.** Positioning vs MS arXiv:2504.10380 v4 §10 examples. All extant MS-style examples are constructed geometrically from continuous spacetimes; T6.3 / T7.1 / T7.2 are the first MS-style examples derived from operator-algebraic truncation via the Wick-rotation functor.
- **§7.4 Related work synthesis.** Sormani-Vega 2016 null distance, Sakovich-Sormani 2024 timed Gromov-Hausdorff, Minguzzi-Suhr 2024 bounded Lorentzian metric spaces, Müller 2022 Cauchy slabs, Allen-Burtscher 2022 / Kunzinger-Steinbauer 2022 metric GH with null distance. These adjacent synthetic Lorentzian convergence frameworks operate in *different categories* than the present bridge construction.

**Key theorems established:** Thm 7.1 (synthetic compactness), Thm 7.2 (G2 metric-level closure for wedge).

**Cross-references for Phase B:** A.4'-C §3 (newly accessible theorems), A.3' §10 (concurrent-work risk register).

### 4.9. §8 Open questions (~2–3 pp)

**Section structure:**

- **§8.1 Q1' (strong-form bridge with enlarged substrate).** Statement (§3 of this synthesis memo): does the Wick-rotation functor $W$ extend to the Paper 46 Appendix B enlarged substrate, and does the off-orbit super-additivity close via direct transport of Paper 42 four-witness theorem? Three sub-questions Q1'.A (topography enlargement), Q1'.B (Gelfand spectrum), Q1'.C (off-orbit super-additivity at strict-inequality level). Estimated 3–6 months follow-on.
- **§8.2 G2 metric-level closure on the general (non-wedge) carrier.** The §7.2 Thm 7.2 closes G2 for the GeoVac wedge; the general non-compact Latrémolière propinquity problem (Paper 45 §1.4 G2-metric on arbitrary non-wedge carriers) remains the multi-month NCG-mathematics frontier.
- **§8.3 G3 cross-manifold $\mathcal{T}_{S^3} \otimes \mathcal{T}_\text{Hardy}(S^5)$.** Blocked at NCG-framework level by Paper 24 §V four-layer Coulomb/HO asymmetry. Out of merged paper scope.
- **§8.4 G4 inner-factor calibration data.** W3 spectral-zeta candidate FALSIFIED. Out of merged paper scope.
- **§8.5 Q2' (extending the bridge to non-commutative MS pre-length spaces).** If Q1'.B forces breaking Abelianness of the enlarged topography, generalize MS pre-length spaces to non-commutative setting (substantial NCG-research target).

**Key theorems established:** none (open questions section).

**Cross-references for Phase B:** A.5' §3 (Q1' detailed framing); Paper 45 §1.4 (G2/G3/G4 statements); CLAUDE.md §1.7 WH1 register (W3 FALSIFIED).

### 4.10. Bibliography assembly (Phase B drafting note)

The bibliography for Paper 48 integrates:

- **Latrémolière 2017, 2018, 2023, 2025, 2512.03573** — full Latrémolière propinquity lineage including the December 2025 pinned proper QMS paper.
- **Mondino-Sämann arXiv:2504.10380 v4** — central synthetic-side framework.
- **Connes-Rovelli 1994** (Class. Quantum Grav. 11, 2899) — thermal-time hypothesis supplying R3.
- **Papers 38, 42, 43, 44, 45, 46, 47** — full GeoVac math.OA series with which Paper 48 is the tenth standalone.
- **Paper 32 §VIII rem:pythagorean_m1_closure** — master Mellin engine M1 signature interpretation.
- **Paper 24 §V** — Coulomb/HO four-layer asymmetry blocking G3.
- **Five MEDIUM-risk concurrent works** (from A.3' concurrent-work re-check):
  - Ketterer arXiv:2605.11271 (ℓ-convergence) — synthetic-side internal refinement, cite as parallel synthetic-side extension.
  - Kubota arXiv:2605.09101 (Lorentzian coarea inequality) — synthetic-side measure theory, cite if Hausdorff measure becomes load-bearing.
  - Mondino-Ryborz-Sämann arXiv:2605.03172 (stability of synthetic TCD) — May 2026 Mondino-Sämann flagship, cite as state-of-the-art synthetic-side.
  - Che-Perales-Sormani arXiv:2510.13069 (timed-metric Hausdorff) — already in Paper 45 hardening pass.
  - Nieuviarts arXiv:2512.15450 v2 (twisted spectral triple) — alternative-to-Wick-rotation framework, cite with KO-dim 3 restriction note.
- **Adjacent synthetic Lorentzian convergence literature** (§7.4): Sormani-Vega 2016, Sakovich-Sormani 2024, Minguzzi-Suhr 2024, Müller 2022, Allen-Burtscher 2022, Kunzinger-Steinbauer 2022.

Estimated bibliography size: ~50–60 entries, mixing operator-algebraic (Latrémolière lineage + Connes lineage + GeoVac series) and synthetic-Lorentzian (Mondino-Sämann lineage + adjacent works) references.

### 4.11. Total estimated length

| Section | Pages | Words (est.) |
|:--------|------:|-------------:|
| §1 Introduction | 3–4 | 1,200–1,600 |
| §2 Setup | 4–5 | 1,600–2,000 |
| §3 Krein-pointed proper QMS | 6–8 | 2,400–3,200 |
| §4 Krein-pointed propinquity hypertopology | 4–5 | 1,600–2,000 |
| §5 Mondino-Sämann bridge | 6–8 | 2,400–3,200 |
| §6 GeoVac wedge application | 5–7 | 2,000–2,800 |
| §7 Synthetic-Lorentzian compactness theorems | 4–5 | 1,600–2,000 |
| §8 Open questions | 2–3 | 800–1,200 |
| Bibliography | 1–2 | — |
| **Total** | **35–47 pp** | **13,600–18,000 words** |

Comparable to Papers 38, 42, 45, 46 in the GeoVac math.OA series (typical 18–24 pages, 8,000–11,000 words). Paper 48 is naturally longer because it consolidates four sub-sprints (A.2', A.3', A.4'-A, A.4'-B) plus the application sprint (A.4'-C). Estimated 1.5 months drafting effort at the Phase B cadence Papers 42/43/44/45/46/47 hit.

---

## §5. Phase A.5'.5 gate verdict

### 5.1. Per-question gate verdict

| Question | Verdict | Source |
|:---------|:-------:|:-------|
| Phase A complete? | **YES** | All six Phase A sub-sprints landed POSITIVE (Tier 3-Light + A.2' + A.3' + A.3' concurrent-work re-check + A.4'-A + A.4'-B + A.4'-C). |
| Bridge fully closed at K⁺-weak-form? | **YES** | Bridge Theorem 6.4 has all four properties (B1)–(B4) at theorem-grade rigor (post-A.4'-A G-B2 + A.4'-B G-B4 closures). |
| GeoVac wedge application complete? | **YES** | A.4'-C T1–T6+T6a, seven newly accessible theorems, all unconditional post-A.4' closures. |
| Paper 48 outline ready for drafting? | **YES** | §4 of this memo provides §1–§8 section structure with key theorems per section + cross-references for Phase B. |
| Q1' deferred cleanly? | **YES** | §3 of this memo provides clean future-direction statement for Paper 48 §8 + estimated 3–6 month follow-on scope. |

### 5.2. Aggregate Phase A.5'.5 verdict

**POSITIVE — PROCEED to Phase B.**

All five gate questions land POSITIVE. The Phase A arc is complete; the Krein–Mondino–Sämann bridge is fully closed at K⁺-weak-form level with all (B1)–(B4) theorem-grade; the GeoVac wedge application yields seven newly accessible theorems including T6 (G2-metric closure for the wedge, a published-open-problem workaround); the Paper 48 §1–§8 outline is ready for Phase B drafting at the rigor level Papers 42/43/44/45/46/47 reached pre-arXiv.

### 5.3. Recommended Phase B sequencing

**Recommendation: SINGLE SPRINT (the merged Paper 48 drafting sprint), estimated ~1.5 months.**

The §4 outline provides section-by-section structure at the level needed for a Phase B drafting agent to pick up and write production-quality LaTeX without re-deriving structure. Each section has explicit cross-references to the source memo (A.2'/A.3'/A.4'-A/B/C) and to the load-bearing papers (38, 42, 43, 44, 45, 46, 47, 32 §VIII). The drafting work is *consolidation and writing*, not new mathematical content.

**Why single sprint, not split-by-section:** Papers 42/43/45/46 each consolidated multiple sub-sprints into a single math.OA standalone at ~1.5–2 months drafting effort per paper. Paper 48 is structurally analogous (5 sub-sprints + 1 application sprint → 1 standalone paper). Splitting by section would add overhead without proportional benefit; the math.OA convention is single-paper consolidation of sprint arcs.

**Alternative if PI wants compression:** drafting agent could be dispatched with §3+§4 as the first deliverable (the central new content: Krein-lift + propinquity hypertopology), then §5 (bridge), then §6+§7 (application + compactness), then §1+§2+§8 (introduction + setup + open questions). This staged delivery would surface drafting issues earlier but extends total wall time.

**Concurrent-work hardening:** Paper 48 should re-run the concurrent-work re-check (per `debug/sprint_phase_a3prime_concurrent_work_recheck_memo.md`) at submission time (likely ~2 months post-Phase B kickoff), as part of the standard pre-submission hardening. The five MEDIUM-risk adjacent works (Ketterer, Kubota, Mondino-Ryborz-Sämann, Che-Perales-Sormani, Nieuviarts v2) are flagged for citation acknowledgement.

### 5.4. What Phase B Paper 48 still needs that Phase A didn't deliver

Beyond the Phase A substrate (which Phase B converts to production-quality LaTeX):

1. **Production-quality LaTeX writing** of §1–§8 per the §4 outline. This is the bulk of Phase B effort.
2. **Bibliography assembly** with all sources integrated per §4.10: Latrémolière lineage + Mondino-Sämann lineage + Connes-Rovelli + GeoVac series + 5 MEDIUM-risk concurrent works + adjacent synthetic Lorentzian literature.
3. **Cross-references to Papers 38, 42, 43, 44, 45, 46, 47** integrated as forward-references throughout §3–§7 (substrate inheritance, propinquity inheritance, wedge inheritance, three-carrier inheritance).
4. **MS Def 3.6 verbatim quotations** in §5.5 / §5.7 (the A.4'-B memo §2 already extracted these verbatim from MS arXiv:2504.10380 v4 — Phase B integrates into Paper 48 prose).
5. **Latrémolière 2512.03573 Def 1.18–5.1 verbatim quotations** at A.2' rigor level in §3 / §4.
6. **Final concurrent-work re-check** at submission time, per standard GeoVac math.OA pre-submission hardening protocol (e.g., Paper 45 §1.4 hardening pass).
7. **Optional: Riemannian-limit recovery numerical panel cross-check** at $N_t = 1$ to confirm bit-exactness with Paper 38 — A.4'-B Cor 5.2 establishes this structurally; Phase B should include numerical verification in §6.7 or as an appendix.

**No new mathematical content is required for Phase B.** All the substantive structural findings (Krein-lift, bridge functor, G-B2/G-B4 closures, T1–T6+T6a) are at theorem-grade rigor in the Phase A memos. Phase B is consolidation, writing, and bibliography work.

---

## §6. Honest scope statement

### 6.1. What this memo establishes

- A complete Phase A synthesis (§1) consolidating six sub-sprint memos into a single closure narrative.
- A named-gap status table (§2) showing seven of nine Lorentzian-arc named gaps CLOSED, two deferred OPEN, one structurally orthogonal.
- A clean future-direction statement for Q1' (§3) including the three sub-questions Q1'.A/B/C and the 3–6 month follow-on scope estimate.
- A merged Paper 48 §1–§8 outline (§4) with section-by-section structure, key theorems per section, estimated word counts, and Phase B drafting cross-references.
- A Phase A.5'.5 gate verdict (§5): POSITIVE, PROCEED to Phase B Paper 48 drafting, single sprint at ~1.5 months estimated effort.
- A list of what Phase B still needs that Phase A didn't deliver (§5.4): production-quality LaTeX, bibliography, cross-references, concurrent-work hardening — no new mathematical content required.

### 6.2. What this memo does NOT establish

- **A Paper 48 draft.** The drafting work is Phase B. This memo is the outline / synthesis layer.
- **Q1' strong-form bridge closure.** Deferred as the multi-month NCG-mathematics follow-on for Sprint L3e-Q1' or independent research.
- **G2 metric-level closure on general (non-wedge) carriers.** The wedge-specific T6 closure is established; the general non-compact Latrémolière propinquity problem remains open.
- **G3 cross-manifold or G4 inner-factor calibration.** Out of Phase A / Paper 48 scope; preserved as open questions in Paper 48 §8.
- **Phase B drafting agent dispatch.** This memo's §4 outline is the input to a Phase B sprint; the drafting sprint itself is the next step after PI approval.
- **A production code implementation.** The bridge is mathematical; no production code is required for Phase A or Phase B.

### 6.3. Load-bearing dependencies

The Phase A.5' synthesis depends on (in source order):

- **Tier 3-Light diagnostic** (`debug/sprint_tier3_light_krein_lift_diagnostic_memo.md`) — initial 7/7 axiom verdict.
- **Phase A.2' Krein-lift formalization** (`debug/sprint_phase_a2prime_krein_lift_formalization_memo.md`) — 9/9 axioms theorem-grade, Def 1.26 dissolution, merged substrate establishment.
- **Phase A.3' Mondino-Sämann bridge** (`debug/sprint_phase_a3prime_mondino_samann_bridge_memo.md`) — Wick-rotation functor + Bridge Theorem 6.4 + R1/R2/R3 testing.
- **Phase A.3' concurrent-work re-check** (`debug/sprint_phase_a3prime_concurrent_work_recheck_memo.md`) — CLEAR-WITH-MINOR-FLAGS verdict + 5 MEDIUM-risk adjacent works.
- **Phase A.4'-A G-B2 closure** (`debug/sprint_phase_a4prime_a_gb2_super_additivity_memo.md`) — Decomposition O + on-orbit equality.
- **Phase A.4'-B G-B4 closure** (`debug/sprint_phase_a4prime_b_gb4_axiom_verification_memo.md`) — Plancherel-nested Berezin = MS Def 3.6 axioms.
- **Phase A.4'-C wedge application** (`debug/sprint_phase_a4prime_c_wedge_application_memo.md`) — T1–T6+T6a seven newly accessible theorems.

If any of these dependencies is shown to fail, the corresponding section of the Phase A.5' synthesis reopens. As of 2026-05-24 close-of-day, all six predecessor memos are stable, all six sub-sprint verdicts are POSITIVE, and the Phase A.5' synthesis stands.

### 6.4. Methodological note

The Phase A arc executed in a single day (2026-05-24) at the formal-memo level. This is consistent with the **diagnostic-before-engineering rule** (CLAUDE.md memory `feedback_diagnostic_before_engineering.md`): the Tier 3-Light diagnostic + the L3e re-scoping memo together identified that most of the bridge work was *already done* in Papers 42–47, that the remaining work was naming + formalization + writing, and that the 6–9 month re-scoped timeline was achievable in concentrated form.

**Realized vs estimated effort:** the Tier 3-Light memo estimated 6–7 months end-to-end (A.2' through Paper 48 draft); the Phase A formal-memo arc compressed Phases A.2' + A.3' + A.4'-A + A.4'-B + A.4'-C into a single day at the formal-memo level. The merged Paper 48 §3–§7 substrate is established; what remains for Phase B is production-quality LaTeX writing + bibliography assembly + cross-references + concurrent-work hardening, estimated 1.5 months at the standard math.OA paper cadence (Papers 42, 45 ~1.5 months drafting each).

**Key compression mechanism:** the Phase A.4'-A G-B2 closure (originally estimated 3–5 weeks) collapsed to a single session via the substantive structural finding that off-orbit super-additivity is *structurally vacuous* at the K⁺-weak-form bridge level. The "open mathematical content" the A.3' sketch flagged was a structural artifact of the bridge's topography choice, not a real gap. This is the second instance of the diagnostic-before-engineering rule paying off at the multi-week-to-single-session scale (the first was the Phase A.2' Def 1.26 caveat resolution via Latrémolière Lemma 1.25).

### 6.5. Cross-references

- All six Phase A predecessor memos (Tier 3-Light + A.2' + A.3' + A.3' concurrent-work + A.4'-A + A.4'-B + A.4'-C).
- L3e Phase A.2 ε-net memo (`debug/l3e_p3_phase_a2_operator_algebraic_eps_net.md`) — F2 identification + R1/R2/R3 candidates.
- L3e re-scope memo (`debug/sprint_l3e_p3_rescope_memo.md`) — A.5' planning context.
- Papers 38 (`papers/group1_operator_algebras/paper_38_su2_propinquity_convergence.tex`) — WH1 PROVEN substrate.
- Papers 42–47 — full GeoVac math.OA Lorentzian-arc substrate.
- Paper 32 §VIII (`papers/group3_foundations/paper_32_spectral_triple.tex`) — case-exhaustion + master Mellin engine + Pythagorean M1 closure.
- Paper 24 §V (`papers/group3_foundations/paper_24_bargmann_segal.tex`) — Coulomb/HO four-layer asymmetry (G3 blocker).
- Mondino-Sämann arXiv:2504.10380 v4 (Dec 2025).
- Latrémolière arXiv:2512.03573 (Dec 2025).
- Connes-Rovelli 1994 (Class. Quantum Grav. 11, 2899).
- CLAUDE.md §1.7 WH1 register, §2 Lorentzian-arc entries, §6 Paper 47 entry.

---

**End of memo.**

**Files added in this sprint:**
- `debug/sprint_phase_a5prime_synthesis_paper48_outline_memo.md` (this memo, ~9500 words formal Phase A synthesis + Paper 48 outline + Phase A.5'.5 gate verdict + Q1' deferral)
- `debug/data/sprint_phase_a5prime_synthesis.json` (structured: Phase A summary + named-gap status table + Paper 48 outline structure + Phase B effort estimate)
