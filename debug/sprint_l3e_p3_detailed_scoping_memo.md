# Sprint L3e-P3 detailed scoping — Pointed Latrémolière propinquity via Mondino-Sämann synthetic Lorentzian GH bridge

**Date:** 2026-05-23 (post-master-theorem, post-quick-test of #1).
**Purpose:** Detailed multi-month scoping for Path P3 of the L3e scoping memo. Pointed Latrémolière propinquity on non-compact Lorentzian Krein spectral triples, via a bridge between the Mondino-Sämann synthetic Lorentzian GH program (mature, 2022-2025) and the Latrémolière metric-spectral-triple propinquity program (mature, 2013-2025). Closes G2-metric (Paper 45 §1.4 / Paper 47 §8 Q1) at the metric level.
**Predecessor:** `debug/sprint_l3e_scoping_memo.md` (Path P3 named as "recommended direction IF committing"; 12-month estimate, MEDIUM-HIGH risk).
**Scope:** This memo is the detailed sprint plan for P3, with phase breakdown, decision gates, concurrent-work audit, and explicit Paper 48/49 deliverables. **The PI has committed to pursuing this direction; this scoping memo is the entry point.**

---

## §1. Closure target

### What G2-metric is

Paper 47 §1.1 + Paper 45 §1.4 G2:

> The non-compact $\R_t$ limit is closed at the **norm-resolvent / spectral level** (Paper 47). The **metric / propinquity** level of this limit on the non-compact carrier requires a non-compact extension of Latrémolière's propinquity, which is not in the published literature as of May 2026: named G2-metric.

The target is a theorem of the form:
$$\Lambda^{\mathrm{pointed}}_L\bigl(\Tcal^L_{n_{\max}, N_t, T(N_t)},\ \Tcal^L_{\sthree \times \R_t}\bigr) \to 0$$
as $(n_{\max}, N_t) \to \infty$ along admissible scaling $T(N_t) \to \infty$, at the **pointed Latrémolière propinquity** level (a new construction, the principal output of Phase B).

The base point on both sides is the **BW wedge vacuum state** $\omega^* = \omega_W^L$ (Paper 43 §4.2). The pointed propinquity construction restricts the metric structure to neighborhoods of $\omega^*$ in the state space, making it well-defined on non-compact carriers.

### What success looks like (12-month horizon)

- **Paper 48** (math.OA standalone, ~Phase B output): *Pointed Latrémolière propinquity for non-compact quantum metric spaces*. Defines pointed quantum compact metric space, pointed tunneling pair, pointed propinquity metric. Proves metric axioms, compact-case agreement with Latrémolière, non-compact-case extension via base-state restriction.
- **Paper 49** (math.OA standalone, ~Phase C output): *Pointed propinquity convergence on truncated Lorentzian Krein spectral triples: G2-metric closure for $\sthree \times \R_t$*. Applies Paper 48 machinery to GeoVac's Lorentzian wedge, proves pointed-propinquity convergence to the genuine non-compact Lorentzian Dirac.
- **Possibly merged**: if scope tightens, a single Paper 48 with both pointed-propinquity definition and GeoVac G2-metric closure as application.

### Why this is the right direction (vs P1/P2/P4)

From L3e scoping memo §3-§6:
- P1 (locally compact framework, Riemannian-first): 12-18 months HIGH risk, less interesting
- P2 ($C_c^\infty$ direct limit): 8-10 months MEDIUM risk, easier but less novel
- P3 (Mondino-Sämann bridge): **12 months MED-HIGH risk, most interesting direction**
- P4 (defer): rational default if not committing, but loses 2026 calendar window

P3 fills the genuine published-literature gap: no operator-algebraic Lorentzian propinquity exists, no bridge between synthetic Lorentzian GH and operator-algebraic propinquity exists. Both communities (Mondino-Sämann synthetic, Latrémolière-Farsi operator-algebraic) are mature; P3 builds the missing bridge.

---

## §2. Phase A — Synthetic-to-operator-algebraic bridge (~2 months)

### Goal

Define a precise correspondence between the Mondino-Sämann synthetic Lorentzian GH framework (causal-diamond ε-nets) and operator-algebraic data (truncated Krein wedge spectral triples). Show that GeoVac's truncated wedge IS a synthetic-Lorentzian ε-net in the appropriate sense.

### Sub-deliverables

**A.1 — Literature deep-dive (2–3 weeks).**
Extract precise definitions from:
- Mondino-Sämann 2024 (Adv. Math.) "Synthetic Lorentzian GH convergence" — pre-length spaces with causal diamonds, time-separation function $\tau(x, y)$.
- Mondino-Sämann 2025 (arXiv:2504.10380, April 2025) "Lorentzian GH convergence and pre-compactness" — causal-diamond ε-nets, pre-compactness theorem.
- Bykov-Minguzzi-Suhr 2024 (arXiv:2412.04311) "Lorentzian metric spaces & GH-convergence: unbounded case" — non-compact extension on synthetic side.
- Braun-Sämann 2025 (arXiv:2506.10852) "Gromov's reconstruction theorem & measured Lorentzian GH" — measure-theoretic version.
- Che-Perales-Sormani 2025 (arXiv:2510.13069) "Synthetic Lorentzian Gromov-Hausdorff convergence" — recent extension.
- Minguzzi 2025 (arXiv:2510.24423) "Results on Lorentzian metric spaces" — Cauchy time functions, structure theorems.

Output: `debug/l3e_p3_phase_a1_literature_audit.md`. Specifically extract:
- Definition of causal-diamond ε-net at scale ε on a Lorentzian length space
- Definition of time-separation-function-induced metric structure
- Pre-compactness theorem statement
- Bykov-Minguzzi-Suhr unbounded extension mechanism

**A.2 — Define operator-algebraic ε-net (3–4 weeks).**
Define the operator-algebraic analog of a causal-diamond ε-net:
- Causal diamond at scale ε → wedge region $W_L$ at modular-flow time $\beta = \varepsilon^{-1}$ (Paper 43 §4.1)
- Time-separation function $\tau(x, y)$ on causally-related point pairs → KMS-state-based correlation function $\tau_L(\omega_1, \omega_2)$ on Krein-positive state pairs
- ε-net structure → truncated operator system $\Op^L_{n_{\max}, N_t}$ at finite cutoff (Paper 44 substrate)

Output: explicit definition of "Krein-pointed ε-net" as an operator-algebraic data structure. Verify it satisfies the synthetic-side axioms (triangle inequality, GH-compatible).

**A.3 — Correspondence theorem (2–3 weeks).**
Prove the correspondence: every Krein-pointed ε-net at scale ε in the operator-algebraic sense corresponds to a causal-diamond ε-net (in the Mondino-Sämann sense) on the underlying spacetime. Specifically, given a truncated Lorentzian Krein spectral triple $\Tcal^L_{n_{\max}, N_t, T}$ at cutoff $(n_{\max}, N_t, T)$, exhibit a causal-diamond ε-net on $\sthree \times S^1_T$ at scale $\varepsilon = O(\max(1/n_{\max}, 1/N_t))$.

Output: `debug/l3e_p3_phase_a3_correspondence_theorem.md` + verification at small panel $(n_{\max}, N_t) = (2, 3)$.

**A.4 — GeoVac wedge as ε-net (2–3 weeks).**
Apply the correspondence to GeoVac's specific truncated Krein wedge constructions (Paper 43 §4 hemispheric wedge + Paper 44 substrate). Verify that the BW vacuum state $\omega_W^L$ corresponds to the canonical Bisognano-Wichmann base point in the synthetic causal-diamond structure.

Output: explicit identification of GeoVac wedge ε-net data with Mondino-Sämann causal-diamond ε-net at every panel cell.

**A.5 — Phase A sprint memo + decision gate (1 week).**
Write up Phase A closure with verdict: positive (proceed to Phase B), negative (terminate program with documented obstruction), or mixed (partial bridge + named gap requires Phase B redesign).

Output: `debug/l3e_p3_phase_a5_sprint_memo.md`.

### Phase A risk

**HIGH.** The bridge between synthetic length-space framework and operator-algebraic framework may surface structural mismatches at the conceptual level:
- Synthetic side uses causal-diamond structure as primary; operator-algebraic side uses Hilbert-space inner product as primary
- Time-separation function on synthetic side is a *finite* number per point pair; operator-algebraic correlation function may be unbounded
- Pre-compactness theorem may not have an operator-algebraic analog (Krein-positive state space at infinite cutoff)

If Phase A returns NEGATIVE, the program terminates here with documentation. ~2 months sunk cost. This is the diagnostic-before-engineering checkpoint.

### Phase A decision gate

Continue to Phase B only if Phase A.5 returns POSITIVE or MIXED-WITH-REDESIGN. If NEGATIVE, write up as cleanly-named negative result (potentially publishable as "structural obstruction to operator-algebraic Lorentzian GH bridge") and return to other GeoVac sprint queue.

---

## §3. Phase B — Pointed Latrémolière propinquity framework (~6 months)

### Goal

Define a pointed Latrémolière propinquity on (possibly non-compact) quantum compact metric spaces. Prove the standard metric properties (triangle inequality, metric structure on equivalence classes). Verify the compact-case agreement with the unpointed Latrémolière propinquity.

### Sub-deliverables

**B.1 — Pointed QCMS definition (4–6 weeks).**
Define **pointed quantum compact metric space** as a triple $(A, L, \omega^*)$ where:
- $A$ is a C*-algebra
- $L : A_{sa}^+ \to [0, \infty]$ is a Lipschitz seminorm (lower semicontinuous, lower bounded, satisfies Lip duality)
- $\omega^* \in \mathcal{S}(A)$ is a distinguished base state

Verify the standard axioms (Latrémolière 2018 / Farsi-Latrémolière 2024) extend to the pointed version. Define pointed quantum compact metric isomorphism (with base-state alignment).

Output: `debug/l3e_p3_phase_b1_pointed_qcms.md`.

**B.2 — Pointed tunneling pair (4–6 weeks).**
Define **pointed tunneling pair** between $(A_1, L_1, \omega_1^*)$ and $(A_2, L_2, \omega_2^*)$: a quintuple $(P, B, \pi_1, \pi_2, \omega^*_P)$ such that:
- $P$ is a C*-algebra with Lipschitz seminorm $L_P$
- $\pi_j : P \to A_j$ are quotient morphisms preserving Lipschitz seminorms (within reach/height)
- $\omega^*_P$ is a base state on $P$ with $\pi_j^*(\omega_j^*) = \omega^*_P$ (alignment)

Define reach, height, and base-state alignment-distance:
- $\mathrm{reach}_j(\pi_j) := \sup\{ |L_P(\pi_j(a_j)) - L_j(a_j)| : a_j \in A_j\}$
- $\mathrm{height}_j(\pi_j) := $ usual Latrémolière height
- $\mathrm{align}(P, \omega^*_P) := \max_j \mathrm{dist}(\omega^*_P, \pi_j^*(\omega^*_j))$

Output: `debug/l3e_p3_phase_b2_pointed_tunneling.md`.

**B.3 — Pointed propinquity metric (4–6 weeks).**
Define **pointed propinquity** as infimum over pointed tunneling pairs:
$$\Lambda^{\mathrm{pointed}}\bigl((A_1, L_1, \omega_1^*), (A_2, L_2, \omega_2^*)\bigr) := \inf_{(P, B, \pi_1, \pi_2, \omega^*_P)} \max\bigl(\mathrm{reach}, \mathrm{height}, \mathrm{align}\bigr)$$

Prove:
- $\Lambda^{\mathrm{pointed}} \ge 0$ with equality iff isomorphism (modulo base-state preserving conjugation)
- Triangle inequality
- Compact-case agreement: when both $(A_j, L_j)$ are compact in the Latrémolière sense, $\Lambda^{\mathrm{pointed}}$ reduces to $\Lambda^{\mathrm{Latr}}$ plus a base-state-alignment correction that vanishes when $\omega^*_j = \omega^*_P$ is canonical.

Output: theorem statement + proof + Paper 48 §3 draft.

**B.4 — Non-compact extension via base-state restriction (4–6 weeks).**
For non-compact carriers (e.g., $C_0(\R)$), the standard Latrémolière framework fails because the Lipschitz ball isn't precompact in $L^\infty$. The pointed version restricts to base-state-neighborhoods:
$$\mathcal{B}_R(\omega^*) := \{\omega \in \mathcal{S}(A) : d_{KR}(\omega, \omega^*) \le R\}$$
where $d_{KR}$ is the Kantorovich-Rubinstein distance induced by $L$.

Define **R-bounded pointed propinquity** $\Lambda^{\mathrm{pointed}, R}$ restricting to states in $\mathcal{B}_R(\omega^*)$, prove it's well-defined for non-compact $(A, L)$ when the base-state-neighborhood is compact.

Output: theorem statement + proof + Paper 48 §4 draft.

**B.5 — Paper 48 draft (2–4 weeks).**
Write up Phase B as standalone math.OA paper. Outline:
- §1 Introduction (motivation: G2-metric, Mondino-Sämann bridge)
- §2 Pointed QCMS definition + axioms
- §3 Pointed tunneling pair + pointed propinquity
- §4 Non-compact extension via base-state restriction
- §5 Compact-case agreement theorem
- §6 Examples (Krein wedge from Paper 43, $C_0(\R)$, Hekkelman-McDonald torus, etc.)
- §7 Open questions

Target: 25-30 pages, 50-60 bibitems. arXiv-ready by ~6 months in.

### Phase B risk

**MEDIUM-HIGH.** Significant original NCG math. Specific risk points:
- B.3 triangle inequality may require strengthened axioms beyond standard Latrémolière
- B.4 R-bounded restriction may have unexpected non-locality (R-dependence may be the wrong limit)
- Reviewer pushback on pointed-propinquity definition matching the published Mondino-Sämann ε-net concept

Mitigation: B.5 paper draft should be reviewed by an external NCG expert before Phase C.

### Phase B decision gate

Continue to Phase C only if Phase B.5 paper draft is internally consistent and externally reviewed positive. Phase B output (Paper 48) is itself a publishable math.OA standalone — if it lands cleanly, the framework's math.OA arc gains a 10th standalone paper regardless of whether Phase C succeeds.

---

## §4. Phase C — G2-metric closure via pointed propinquity (~4 months)

### Goal

Apply Paper 48's pointed propinquity machinery to the specific G2 closure question. Prove pointed-propinquity convergence of the truncated Lorentzian Krein wedge spectral triples to the genuine non-compact $\sthree \times \R_t$ limit.

### Sub-deliverables

**C.1 — BW vacuum as canonical base state (4–6 weeks).**
Identify the Bisognano-Wichmann wedge vacuum state $\omega_W^L = e^{-K_\alpha^W}/Z$ (Paper 43 §4.2) as the canonical base state at every truncated cutoff $(n_{\max}, N_t, T)$ and at the continuum limit. Verify:
- $\omega_W^L$ is Krein-positive (lives in $\Kplus$)
- $\omega_W^L$ has well-defined Kantorovich-Rubinstein neighborhoods at every cutoff
- $\omega_W^L$ converges (in some appropriate sense) along the (n_max, N_t, T) sequence

Output: structural identification of BW base state across cutoffs.

**C.2 — Pointed propinquity bound (4–6 weeks).**
Prove:
$$\Lambda^{\mathrm{pointed}, R}_L\bigl(\Tcal^L_{n_{\max}, N_t, T(N_t)},\ \Tcal^L_{\sthree \times \R_t}\bigr) \to 0$$
as $(n_{\max}, N_t) \to \infty$ along admissible scaling $T(N_t) \to \infty$, with base states $\omega_W^L$ at each cutoff and continuum limit, with $R$ taken to capture the physically relevant state-space neighborhood.

Combine Paper 47's two-rate hybrid (inner propinquity + outer norm-resolvent) with Paper 48's pointed-propinquity framework. The inner arrow gets pointed-restricted automatically; the outer arrow's norm-resolvent convergence implies state-space convergence in a restricted neighborhood.

Output: theorem statement + proof + Paper 49 §4 draft.

**C.3 — Physical applications (2–4 weeks).**
Connect to bound-state QFT in curved spacetime:
- Hawking radiation: pointed propinquity at the cigar-geometry base state gives operator-algebraic statement of Hawking thermal emission
- Unruh effect: pointed propinquity at the Rindler wedge base state recovers Unruh temperature
- Bound-state QFT (Lamb shift, hyperfine, etc.) in curved background: pointed propinquity provides convergence framework for cutoff-dependent calculations

Output: at minimum, structural-correspondence statements (analogous to Sprint TD Track 4 / Sprint Unruh-pendant, but at the pointed-propinquity level). Full bound-state QFT applications may need separate sprints.

**C.4 — Paper 49 draft (2–4 weeks).**
Write up Phase C as standalone math.OA paper. Outline:
- §1 Introduction (G2-metric closure motivation, summary of Paper 48 machinery)
- §2 Setup (Krein wedge, BW base state, R-bounded neighborhoods)
- §3 Main theorem (pointed propinquity convergence)
- §4 Proof (Paper 48 machinery + Paper 47 two-rate hybrid)
- §5 Three-carrier identification at the pointed level (Paper 47 §6 analog)
- §6 Physical applications (Hawking, Unruh, bound-state QFT)
- §7 Open questions

Target: 20-25 pages. arXiv-ready by ~12 months in.

### Phase C risk

**MEDIUM.** Phase C is a direct application of Phase B; if B lands, C should follow without surprises. Specific risk:
- C.2 may surface that the inner-propinquity convergence (Paper 45/46) needs strengthening to be pointed-compatible
- C.3 physical applications may require additional sprint-scale work (Hawking calculation at the pointed-propinquity level, etc.)

### Phase C decision gate

Phase C output (Paper 49) closes G2-metric at the pointed-propinquity level. If C.3 physical applications surface obstructions, Paper 49 still publishable as math.OA standalone with G2-metric closure proper; downstream physics applications become separate sprints.

---

## §5. Decision gates summary

| Gate | At | Decision | If POSITIVE | If NEGATIVE |
|:-----|:---|:---------|:------------|:------------|
| 1 | End of Phase A.5 (~2 months) | Bridge correspondence exists? | Proceed to Phase B | Terminate program; document obstruction; return to GeoVac sprint queue |
| 2 | End of Phase B.5 (~7 months) | Pointed propinquity is a valid metric? | Paper 48 published; proceed to Phase C | Paper 48 lands as partial result; revisit Phase C scope |
| 3 | End of Phase C.4 (~12 months) | G2-metric closed at pointed-propinquity level? | Paper 49 published; G2-metric CLOSED | Paper 49 lands as G2-metric partial closure; named gap documented |

---

## §6. Concurrent-work risk audit

L3e scoping memo §6 noted moderate scoop risk. **Refreshed assessment (2026-05-23):**

| Risk | Source | Timeline | Mitigation |
|:-----|:-------|:---------|:-----------|
| Mondino-Sämann move into operator algebras | Mondino-Sämann lineage (Italian synthetic geometry community) | 2026-2028 plausibly | Phase A.1 literature audit catches any post-May-2026 publication |
| Latrémolière moves into Lorentzian | Latrémolière + Farsi-Latrémolière | 2026-2028 unlikely (trajectory stays Riemannian per `debug/l3_literature_audit_memo.md` §1) | Phase A.1 + Phase B.5 re-audit |
| Independent pointed-propinquity definition | Anyone in NCG metric geometry | Always possible | Pre-submit Paper 48 to arXiv ASAP |
| Hekkelman-McDonald non-compact extension | Hekkelman + McDonald | 2026-2027 possibly | Their framework is Tauberian, not pointed; may complement rather than scoop |

**Recommendation:** Run a literature re-check at the start of each Phase (A.1, B.1, C.1) — total ~6 hours of work spread over 12 months.

---

## §7. Output / Paper deliverables

| Output | Phase | Type | Estimated submission |
|:-------|:------|:-----|:---------------------|
| Phase A sprint memo | A.5 | Internal | ~2 months in |
| Paper 48 (pointed propinquity) | B.5 | math.OA standalone | ~7 months in (arXiv preprint) |
| Phase C sprint memo | C.4 | Internal | ~12 months in |
| Paper 49 (G2-metric closure) | C.4 | math.OA standalone | ~12 months in (arXiv preprint) |

**The Phase B output (Paper 48) is the most original-NCG-math contribution of the program.** It's a publishable math.OA paper regardless of whether Phase C succeeds. If Phase B lands and Phase C surfaces obstructions, the framework still ships a 10th math.OA standalone defining pointed propinquity.

---

## §8. Strategic implications

### For the GeoVac framework

This is the framework's largest math.OA commitment yet. After Papers 38/39/40/42/43/44/45/46/47 (today's count: 9 math.OA standalones), Papers 48/49 would bring the math.OA arc to **11 standalones over ~14 months**. This represents the framework's claim to original-NCG-math contributions beyond the spectral-triple-on-S³ infrastructure — specifically, an extension of Latrémolière propinquity to non-compact carriers via Mondino-Sämann bridging.

### For the broader NCG community

The bridge between synthetic Lorentzian GH (Mondino-Sämann lineage) and operator-algebraic propinquity (Latrémolière lineage) is the unfilled gap between two mature communities. Building it would:
- Open Lorentzian QFT applications of propinquity machinery (Hawking, Unruh, bound-state QFT in curved spacetime)
- Open synthetic-geometry interpretations of operator-algebraic Lorentzian constructions
- Position GeoVac at the bridge — a "transmitter" between the two communities

### For physics applications

Phase C.3 connects directly to physics:
- **Hawking radiation:** Sprint TD Track 4 reproduced $T_H = 1/(8\pi M)$ at the structural-correspondence level (M1 Hopf-base measure signature). Pointed propinquity would lift this from "structural correspondence" to "literal pointed-metric convergence" — analogous to how Sprint L1 lifted four-witness Wick rotation from "structural correspondence" to "literal identification at operator-system level."
- **Unruh effect:** Sprint Unruh-pendant similarly. Pointed propinquity would lift the structural identification.
- **Bound-state QFT on curved background:** the framework's Lamb shift, hyperfine, etc., calculations currently use flat-space Dirac-S³. Pointed propinquity opens the door to curved-spacetime versions.

These are downstream physics applications, NOT in the Phase C scope, but the framework would unlock them.

---

## §9. Recommendation

**OPEN Sprint L3e-P3 with Phase A as the first concrete sub-sprint.**

Specifically: open Sprint L3e-P3-A (literature deep-dive + bridge construction) as the next concrete track. This is the LOW-risk diagnostic phase that determines whether Phases B and C are tractable. Estimated 2 months.

Decision at Phase A.5: if the bridge works (POSITIVE or MIXED-WITH-REDESIGN), commit to Phase B. If NEGATIVE, terminate program with documentation.

**Do NOT open the full Phase A-B-C stack as a single multi-month commitment** — the Phase A diagnostic determines whether the full program is reachable.

### Alternative: "low-fly" first sub-sprint

If even 2 months feels too much commitment up front, an alternative entry point is **Sprint L3e-P3-A.1 only** (literature deep-dive, ~3 weeks). This produces a comprehensive literature audit memo and a refined Phase A.2 plan with cost estimates. After A.1, decide whether to continue with A.2-A.5.

### What's at stake

If Phase A goes positive (probability ~60%, my honest assessment): the framework opens a multi-month math.OA program that produces 2 standalone math papers and unlocks Lorentzian QFT on curved spacetime.

If Phase A goes negative (probability ~40%): the framework documents a clean negative result on the operator-algebraic Lorentzian GH bridge, recovers to the GeoVac sprint queue, ~2 months sunk cost.

The diagnostic-before-engineering rule strongly favors opening Phase A as a contained 2-month commitment with explicit decision gate, NOT committing to the full 12-month stack.

---

## §10. Honest scope

This memo:
- IS a scoping memo for a multi-month commitment
- IS NOT a substitute for the actual Phase A work
- Is informed by the May-17 L3 literature audit memo but does NOT re-do the literature review
- Assumes the Mondino-Sämann lineage and Latrémolière lineage are both correctly characterized as "no published bridge exists" — a Phase A.1 deliverable is to re-verify this

**Confidence levels:**
- HIGH on the structural-gap statement (no published bridge between Mondino-Sämann and Latrémolière) per May-17 literature audit
- HIGH on the path P3 being the "most interesting direction if committing" per L3e scoping memo
- MEDIUM on the Phase A success probability (~60%) — based on analogy to other math.OA sprints; could be higher or lower depending on whether specific structural obstructions surface in Phase A.1-A.3
- MEDIUM on the 12-month total cost estimate — based on Papers 38/39/40/42-47 cadence (each ~weeks-months at our cadence); cumulative 14-month estimate could vary by ±4 months

**No follow-on sprint is auto-opened by this memo.** The PI's "yes please" was the trigger to draft this scoping memo; opening Sprint L3e-P3-A is the next PI decision.

---

## §11. Immediate next steps (if PI confirms)

1. **Open Sprint L3e-P3-A.1** (literature deep-dive, ~3 weeks). Three sub-tasks:
   - Mondino-Sämann lineage: extract precise definitions of causal-diamond ε-net, time-separation function, pre-compactness theorem from the 6-paper canonical list (§2 Phase A.1).
   - Latrémolière lineage: review the 2024-2025 updates to Latrémolière propinquity, looking for any pointed / locally compact extension proposals.
   - Hekkelman-McDonald 2024 (arXiv:2412.00628): re-read with non-compact-extension lens; assess whether their Tauberian framework can complement pointed propinquity.

2. **Output**: `debug/l3e_p3_phase_a1_literature_audit.md` (~5000 words, comparable to the L3 literature audit memo).

3. **Decision at end of A.1**: continue with A.2-A.5 (the bridge construction, ~6-8 weeks total) or pause for PI review.

The first concrete step is a 3-week literature deep-dive. Total contractual commitment if opened: 3 weeks. Full Phase A: 8-9 weeks. Full P3 program: ~12 months.
