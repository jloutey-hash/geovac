# Sprint L3e scoping memo — Non-compact Latrémolière propinquity on $\sthree \times \R_t$ (G2-metric closure)

**Date:** 2026-05-23
**Purpose:** Scoping-only diagnostic for the multi-month frontier sprint that would close G2 (Paper 45 §1.4 de-compactification limit) at the **metric / propinquity** level, complementing Paper 47's norm-resolvent / spectral closure. Identifies whether this is a 4–8-month sprint, a 12–18-month sprint, or genuinely multi-year original-NCG-math.
**Inputs:**
- Paper 45 §1.4 G2 statement (de-compactification, post-Paper-47 update naming "G2-metric").
- Paper 47 §1.1 strategic reframing (AF inductive limit cannot capture $C_0(\R)$; norm-resolvent is the right tool at the spectral level).
- Paper 47 §7.3 (`sec:vs_latremoliere`) candidate routes: (a) locally compact extension via $C_c^\infty$ + approximate identity; (b) functional-calculus / weighted Lipschitz seminorm.
- L3 scoping memo (`debug/l3_scoping_memo.md`, 2026-05-17) — original Sprint L3c estimate at multi-month scale for the propinquity-level (Path P2).
- L3 literature audit memo (`debug/l3_literature_audit_memo.md`, 2026-05-17) — established that NO published Lorentzian propinquity exists as of May 2026; bridge from synthetic Lorentzian GH (Mondino-Sämann lineage) to operator-algebraic propinquity (Latrémolière lineage) is the original-NCG-math gap.
- Sprint L3c-α + L3c-γ + Paper 47 closure (2026-05-23) — G2-spectral is closed; G2-metric is the named-residual.

**No production code. No theorem claim. No paper edits applied.** This is a sprint scoping memo only.

---

## §1. The problem

### Statement of G2-metric

Paper 45 §1.4 G2 (post-Paper-47 update) states:
> The non-compact $\R_t$ limit ... is closed at the **norm-resolvent / spectral level** (Paper 47). The **metric / propinquity** level of this limit on the non-compact carrier requires a non-compact extension of Latrémolière's propinquity, which is not in the published literature as of May 2026: named G2-metric or equivalently G1-genuine on the non-compact substrate (Sprint L3c-$\beta$ / L3e in our internal taxonomy).

The problem: define a Latrémolière-style quantum metric on **non-compact** spectral triples (specifically the Lorentzian Krein triple $\mathcal{T}^L_{S^3 \times \R_t}$), and prove that the finite-cutoff truncations $\mathcal{T}^L_{n_{\max}, N_t, T(N_t)}$ converge to it under this metric as $(n_{\max}, N_t) \to \infty$ along admissible scaling.

### Why this is structurally hard

Latrémolière's propinquity (2013–2023, $\sim$30+ papers in his and Farsi-Latrémolière's lineage) is defined exclusively for **compact** quantum metric spaces. The compact-carrier requirement enters via:

1. **Finite $L^\infty$ norm.** The Lip-seminorm $L(a) = \|[D, a]\|_{op}$ is defined on the dense $\ast$-subalgebra of Lipschitz elements, which on a compact carrier is dense in $C(X)$. On a non-compact carrier $X$, $C_0(X)$ (continuous functions vanishing at infinity) replaces $C(X)$; Lipschitz elements are not dense in $C_0(X)$ in the sup norm (e.g., $f(t) = \cos(t)$ is Lipschitz but $\notin C_0(\R)$).

2. **Compact state space.** The state space $\mathcal{S}(C(X))$ of a compact quantum metric space is itself compact in the weak-$\ast$ topology; the Kantorovich-Rubinstein distance is well-defined on this compact state space. On $C_0(X)$ for non-compact $X$, the state space is no longer compact (e.g., $\delta$-functions at infinity escape).

3. **Berezin reach bound.** Paper 38's L4 Berezin approximate-identity bound $\|B(f) - P M_f P\|_{op} \le \gamma_n \|\nabla f\|_\infty$ uses $\|\nabla f\|_\infty$ control on a compact set. On non-compact, $\|\nabla f\|_\infty$ is a strong condition (essentially requires $f$ Schwartz).

4. **Latrémolière 2018 §5 inductive-limit propinquity.** Defined only for AF C*-algebras (totally-disconnected spectrum). Cannot capture $C_0(\R)$ (Paper 47 §1.1).

### The structural gap

The natural-NCG version of the non-compact extension is the **non-unital spectral triple** / **locally compact quantum metric space** framework. There IS a 2004–2010 literature on this (e.g., Gayral-Gracia-Bondia-Iochum-Schücker-Várilly 2004 "Moyal planes are spectral triples"), but it focuses on existence-of-spectral-triple, not propinquity-convergence.

**No published Latrémolière propinquity construction exists for non-compact carriers.** This is the central scoping gap.

---

## §2. Two-axis problem decomposition

Decompose the L3e problem along two independent axes:

### Axis 1: Compactness extension (compact → non-compact)

- **C0-axis-A** (locally compact framework): replace $C(X) \to C_0(X)$, $\mathcal{B}(\mathcal{H}) \to $ multipliers of compact operators, etc.
- **C0-axis-B** (compactly-supported framework): restrict to $C_c^\infty(X)$-multipliers, redefine propinquity on the dense subspace.
- **C0-axis-C** (weighted framework): replace $L^\infty$ Lipschitz seminorm with weighted $L^p$ seminorm (e.g., Sobolev $H^1$) that's well-defined on non-compact.
- **C0-axis-D** (pointed-metric-space framework): pointed Gromov-Hausdorff distance with a chosen basepoint; analog of Gromov's pointed GH for non-compact metric spaces.

### Axis 2: Signature (Riemannian → Lorentzian)

- **Sig-axis-R** (Riemannian non-compact): the Latrémolière construction for $\R^d$ or other non-compact Riemannian manifolds. This is the "natural" first step — if it can't be done in Riemannian, the Lorentzian version is hopeless.
- **Sig-axis-L** (Lorentzian Krein non-compact): the target. Requires both Lorentzian (Krein-space) AND non-compact extensions, plus their compatibility.

**The naive product is 4 × 2 = 8 axis combinations.** Most are not interesting: only the L3e-relevant combinations are:

| Combination | Status | Comment |
|:------------|:-------|:--------|
| C0-A × Sig-R (locally compact Riemannian) | Closest existing literature; partial sketches by Bellissard-Marcolli on non-commutative tori, but no published Latrémolière-style propinquity. | First-step target if pursuing L3e. |
| C0-A × Sig-L (locally compact Krein) | TARGET (Paper 45 §1.4 G2-metric). | Multi-month. |
| C0-B × Sig-R (compactly-supported Riemannian) | Restriction to $C_c^\infty$ test functions; the "easy" intermediate step. | Possibly sprint-scale. |
| C0-B × Sig-L (compactly-supported Krein) | Restriction-of-restriction. May be the path of least resistance. | Possibly 6-month. |
| C0-C × Sig-R/L (weighted) | Requires entirely new propinquity definition; unclear whether weighted-Lip propinquity is even well-posed. | Open problem. |
| C0-D × Sig-R (pointed Riemannian) | Connects to Gromov's pointed-GH; pre-Latrémolière literature exists but not yet operator-algebraic. | Open problem. |

---

## §3. Candidate paths

### Path P1 — Locally compact framework, Riemannian-first

**Sprint sequence:**
1. (~3 months) Define non-compact Latrémolière propinquity for $C_0(\R) \otimes M_n$ as a quantum metric space. Bridge to Connes' non-unital spectral triple framework. Test on the Riemannian Dirac on $\R$.
2. (~3–6 months) Extend to higher-dimensional carriers ($\R^d$ via covering map of $T^d$, building on Hekkelman-McDonald 2024).
3. (~3–6 months) Lift to Krein-signature carrier (use Paper 44 substrate; the spatial $S^3$ stays compact, only temporal $\R_t$ is non-compact).
4. (~3 months) Apply to L3c-γ's two-rate hybrid composite: replace the outer norm-resolvent arrow with a non-compact propinquity arrow, giving a single uniform Latrémolière bound.

**Total: 12–18 months.**

**Risk: HIGH.** Each step is novel. Step 1 alone is a publishable J. Funct. Anal. / J. Op. Theory paper. The composite is more.

**Output:** a non-compact Latrémolière framework + the L3e theorem at the metric level.

### Path P2 — Compactly-supported framework

**Sprint sequence:**
1. (~2 months) Restrict the propinquity to $C_c^\infty(S^1_T)$-multipliers on the temporal factor. The cb-norm and Berezin reach bounds inherit verbatim from Paper 45's $C(S^1_T)$-construction.
2. (~3 months) Define the limit object as a *direct limit* of compactly-supported propinquities: $\Lambda^{cs}(T_{S^3 \times \R}, T_{S^3 \times S^1_{T_n}}) := \lim_{T \to \infty} \Lambda(T_{S^3 \times S^1_T}|_{C_c^\infty}, T_{S^3 \times S^1_{T_n}}|_{C_c^\infty})$.
3. (~2 months) Verify the limit is well-defined and gives a metric (not just a pseudo-metric) on a natural class of triples.
4. (~2 months) Apply to G2-metric closure: $\Lambda^{cs}(T^L_{n_{\max}, N_t, T(N_t)}, T^L_{S^3 \times \R_t}) \to 0$ under admissible scaling.

**Total: 8–10 months.**

**Risk: MEDIUM.** The compactly-supported restriction is a natural intermediate step; the limit construction (Step 2) is genuinely new but technically straightforward.

**Output:** a $C_c^\infty$-Latrémolière propinquity + L3e at the metric level on compactly-supported test functions.

**Honest scope:** This gives metric-level convergence only on compactly-supported test data — not all of $C_0(\R)$. Whether this counts as "closing G2-metric" is a judgment call (it captures the physically interesting content — observables localized in space and time — but doesn't generalize to all bounded continuous functions).

### Path P3 — Pointed Gromov-Hausdorff with Mondino-Sämann bridge

**Sprint sequence:**
1. (~2 months) Adapt the Mondino-Sämann 2025 synthetic Lorentzian GH framework (causal-diamond ε-nets) to operator-system data: replace the time-separation function with a Krein-space modular-flow generator (use Paper 43's BW-α / Tomita-Takesaki construction).
2. (~6 months) Define a *pointed* Latrémolière propinquity on the Krein wedge (base point = the BW vacuum state). This is the operator-algebraic analog of Gromov's pointed-GH for non-compact metric spaces.
3. (~4 months) Prove G2-metric closure at the pointed-propinquity level.

**Total: 12 months.**

**Risk: MEDIUM-HIGH.** This is the most original direction — bridging the Mondino-Sämann synthetic Lorentzian GH program (mature, ~7 papers 2022-2025) and the Latrémolière metric-spectral-triple program (mature, ~30+ papers 2013-2025). L3 literature audit identifies this bridge as the unfilled gap (`debug/l3_literature_audit_memo.md` §1).

**Output:** original NCG-math result — first operator-algebraic Lorentzian GH with a propinquity-style metric. Publishable J. Geom. Phys. / Adv. Math.

**Honest scope:** This is a SEPARATE math.OA program in scope. The result would be a 10th math.OA standalone (Paper 48 or 49). It opens its own line of follow-on questions.

### Path P4 — Multi-year defer

**Recognition:** the framework's structural-skeleton-scope statement (CLAUDE.md §1.7 multi-focal-composition wall taxonomy) puts G2-metric in the W3-class category — open question requiring new external mathematics. The cleanest move may be:

1. Declare G2-metric as "frontier-of-field" in Paper 47 §8 Q1.
2. Wait for the Latrémolière community to develop a non-compact propinquity (it is a natural follow-up to their compact program, and 2-3 papers in that direction may appear 2026-2028 organically).
3. Reapply to GeoVac once the framework exists.

**Total: 0 months active GeoVac work. 2–4 years waiting.**

**Risk: LOW.** No commitment. May get scooped on the bridge (Path P3), but the bridge is original-NCG-math work and is naturally a separate paper anyway.

**Output:** continued GeoVac progress on other open questions, with G2-metric documented as "awaiting external development."

---

## §4. Cost-benefit table

| Path | Months | Risk | Output | Recommended? |
|:-----|:------:|:----:|:-------|:------------:|
| P1 (locally compact Riemannian-first) | 12–18 | HIGH | Non-compact Latrémolière framework + L3e | NOT NOW |
| P2 ($C_c^\infty$ direct limit) | 8–10 | MEDIUM | $C_c^\infty$-restricted L3e | MAYBE (if PI wants the metric-level closure within a year) |
| P3 (pointed propinquity via Mondino-Sämann bridge) | 12 | MED-HIGH | Original NCG-math paper, GeoVac Paper 48/49 | MAYBE (if interested in opening a new math.OA arc) |
| P4 (defer until external math develops) | 0 active | LOW | Continued GeoVac progress elsewhere | **RECOMMENDED** |

---

## §5. The case for P4 (deferral)

Four arguments for declaring G2-metric a frontier and moving on:

1. **G2 is already closed at the level GeoVac currently needs.** Paper 47's two-rate hybrid gives norm-resolvent convergence, which captures the physical content ("the Lorentzian Dirac on the compactified spacetime converges to the Lorentzian Dirac on the genuine spacetime"). This is the convergence notion used by mathematical physicists when discussing spectral-triple limits.

2. **GeoVac's other open questions are more urgent.** The framework has several open questions where structural progress is possible at sprint scale:
   - Granular derivation of Input I3 from CH matrix elements (Pythagorean orthogonality, ~1 week, immediate follow-on to today's sprint).
   - L3c-β enlarged-substrate parametric stability (6–10 weeks).
   - Z>20 cliff (BBB93 implementation, multi-week).
   - W1c-residual chemistry orthogonality wall.
   - Physics-side precision catalogue extensions.
   - The state-side dictionary direction (Paper 34 §III.28 apparatus identity, candidate Sprint 3+ targets).

   Each of these is plausibly sprint-scale at our current cadence. G2-metric is multi-month and would lock out 8–18 months of session-time.

3. **The literature is plausibly moving in this direction without us.** Latrémolière's lineage (2017–2025) has ~30+ Riemannian papers; their compact-to-locally-compact extension is a natural next step that someone in their community is likely to take 2026–2028. Waiting means we get the framework for free. Doing it ourselves means we publish a math.OA paper in territory where the natural community is the Latrémolière + Connes-vS school — they would scoop it if they want it.

4. **The marginal value of G2-metric for the GeoVac physics program is unclear.** Paper 47's norm-resolvent closure is sufficient for any of the framework's physics applications (Lamb shift, hyperfine, chemistry, etc.) — none of these require Lipschitz/quantum-metric structure on the non-compact carrier. The metric-level closure has interest as math.OA work but not direct downstream physics use within the framework.

---

## §6. The case for P3 (Mondino-Sämann bridge)

If a multi-month commitment is desired, P3 is the most interesting direction:

1. **It's the genuine gap in the literature.** The bridge between synthetic Lorentzian GH and operator-algebraic propinquity has not been built by anyone. If GeoVac builds it, we publish a math.OA standalone in JFA or J. Geom. Phys. territory.

2. **Scope is well-defined.** The Mondino-Sämann ε-net + causal-diamond construction is fully published; lifting it to operator-system data is a clean, contained task. Latrémolière's BW-α + Tomita-Takesaki machinery (Paper 42/43) is the natural Krein-side ingredient. Three-step program (P3 Steps 1, 2, 3) is articulated in §3 above.

3. **Scoop risk is moderate.** Mondino-Sämann may move into operator algebras on their own, but they're synthetic-geometry mathematicians; the operator-algebraic side is plausibly out of their natural toolkit. Latrémolière may move into Lorentzian propinquity but he's been Riemannian for 12+ years and his recent trajectory (Farsi-Latrémolière 2024/2025) stays Riemannian. The 2026 calendar window for this is plausibly open.

4. **Marginal cost ~12 months, marginal value ~Paper 48/49.** This is comparable cost-to-value as Paper 38/40/45/46 (each took ~weeks-months and produced a math.OA standalone in a similar territory). The cadence is sustainable.

**However:** P3 requires committing to ~12 months of session-time at the L3 priority, which excludes the other 5-6 sprint-scale follow-ons in §5 point 2. The PI's strategic choice is the load-bearing question.

---

## §7. Recommendation

**Default recommendation: Path P4 (defer).** Take the win from Paper 47's norm-resolvent closure; document G2-metric as a frontier-of-field open question; move to the other queued sprint-scale work.

**If the PI wants to commit to a multi-month math.OA program** beyond what Papers 45/46/47 already deliver, **Path P3 is the recommended direction** (Mondino-Sämann bridge). Path P1 is too risky; Path P2 ($C_c^\infty$ restriction) is the easier but less interesting middle ground.

**Decision-gate-style breakdown:**
- "Continue at sprint cadence on multiple smaller open questions" → P4.
- "Open one big multi-month math.OA program now" → P3.
- "Do the easier intermediate step" → P2.
- "Aim for the hardest result with highest risk" → P1.

---

## §8. Honest scope of this memo

This is a **scoping memo only**. No proof attempts. No production code. No paper edits applied.

The verdict is **GO at multi-month scale via P3 IF the PI commits**, or **DEFER (P4) if continuing at current sprint cadence**. The cost estimates are based on analogy to the L3b-2 + L3b-2f-β + L3c arcs (each ~weeks-months at our current cadence) extrapolated to multi-stage original-NCG-math projects.

**Confidence:** HIGH on the structural-gap statement (§1.2): non-compact Latrémolière propinquity does not exist in the published literature, and the compact-carrier requirement is structural. MEDIUM on the path-cost estimates (§3): these are educated guesses based on existing-sprint cadence; actual costs could vary by 2× in either direction. HIGH on the "P4 is rational default" recommendation: G2-metric does not block any other GeoVac progress, and the marginal value is unclear.

**No follow-on sprint is auto-opened by this memo.** The PI decides between P3, P2, or P4.

---

## §9. Open questions for PI

1. **Commit to L3e at all, or defer?** P4 (defer) is the rational default given the framework's other open questions. Is there a specific reason to prioritize G2-metric over the other queued sprints in §5 point 2?

2. **If committing, P3 vs P2 vs P1?** P3 (Mondino-Sämann bridge) is the most interesting and the closest to publishable math.OA territory; P2 is the easier but less novel option; P1 is the hardest with highest risk.

3. **If committing to P3, how to structure the multi-month commitment?** Pure multi-month sprint, or multi-stage with explicit decision gates between Steps 1, 2, 3?

4. **Does declaring G2-metric "frontier" create scoop risk?** Probably modest — the bridge between synthetic Lorentzian GH and operator-algebraic propinquity is original work and no community is currently positioned to do it on a 2026 timeline. But if the PI wants to lock priority, opening P3 now is the way.
