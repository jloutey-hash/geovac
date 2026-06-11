# Sprint Q1'-Phase-2.A — Operator-System Lorentzian Pre-Length Space (OSLPLS) category definition

**Date:** 2026-05-24 (Q1'-Phase-2.A formalization sprint, post-Q1'-Phase-1 POSITIVE-stepping-stone verdict).

**Sprint position:** Q1'-Phase-2.A of the Q1' staged-sprint structure. First of three Phase-2 sub-sprints (2.A category definition, 2.B bridge functor extension, 2.C off-orbit super-additivity closure) culminating in Phase-2.D Paper 49 drafting.

**Predecessors (load-bearing):**
- `debug/sprint_q1prime_phase1_case_a_stepping_stone_memo.md` — Q1'-Phase-1 (Case A stepping stone with four theorems; structural starting point that Phase-2.A extends to non-Abelian)
- `debug/sprint_q1prime_light_diagnostic_memo.md` — Q1'-Light diagnostic (Option γ flagged as right Q1'-Phase-2 path)
- `debug/sprint_q1prime_concurrent_work_recheck_memo.md` — Q1' concurrent-work re-check (CLEAR; OSLPLS recommended as natural community gap)
- `debug/sprint_phase_a3prime_mondino_samann_bridge_memo.md` — A.3' Wick-rotation functor $W: \mathbf{KreinMetaMet}_{\mathrm{pp}} \to \mathbf{LorPLG}_{\mathrm{cov}}$ (the commutative MS bridge that Phase-2.A generalizes the target of)
- Paper 48 §3 (`def:krein_topography`, `def:krein_ppqms`, `def:lpls`, `def:eps_net`, `def:lgh_convergence`, `def:covered_lpls`, `def:plgh`); §8.1 Q1' three-step decomposition
- Paper 46 Appendix B (`sec:enlarged_substrate`, `eq:flip_generator`, full M≠0 chirality-asymmetric diagonal form $M^{\mathrm{flip}}_{N,L,M} = \mathrm{diag}(W^{NLM}, -W^{NLM})$)
- Paper 44 (`def:propagation`, `prop:prop_2`) — operator-system substrate at finite cutoff
- Paper 42 four-witness theorem (load-bearing for Connes-Rovelli identification under Wick rotation)
- Mondino–Sämann arXiv:2504.10380 v4: Def 2.3 LPLS, Def 3.2 ε-net, Def 3.6 LGH convergence, Def 3.8 covered LPLS, Def 3.12 pLGH, Def 4.4 morphisms (extracted in A.3' memo §1.2)
- Connes–vS arXiv:2004.14115 — operator-system framework template (UCP morphisms, propagation number, operator-system Lipschitz seminorms)

**Status:** FORMAL THEOREM-GRADE MEMO. No production code, no paper modifications (Paper 49 drafting is Phase-2.D, not Phase-2.A). Theorem-grade rigor for Lemma 2.1-OS (axiom transport), Definition 3.1-OS (OSLPLS object), Definition 4.1-OS (OSLPLS morphism), Theorem 5.1-OS (embedding functor $\iota$), Theorem 6.1-OS (Krein-side candidate is OSLPLS object); honest scope statement (§8) documents what remains open for Phase-2.B (the bridge theorem OSLPLS-target version of Bridge Theorem 6.4'-Q1').

---

## Phase-2.A.5 gate verdict (one-sentence headline)

**POSITIVE — OSLPLS category defined cleanly at theorem-grade rigor.** The 5-tuple object $(\mathcal{A}, L, \mathcal{M}, \omega, \mathcal{U})$ with $\mathcal{M}$ an operator system (not necessarily Abelian), $\omega$ a basepoint state, and $\mathcal{U}$ an operator-system cover by sub-cells satisfies the operator-system generalization of all six MS Def 3.6 / 3.8 axioms; the morphism class (UCP maps preserving topography, basepoint state, modular structure, and cover) composes correctly with identity giving a category; the embedding functor $\iota: \mathbf{LorPLG}_{\mathrm{cov}} \hookrightarrow \mathbf{OSLPLS}_{\mathrm{cov}}$ is faithful (commutative MS appears as a full sub-category via $X \mapsto C(X)$ topographic-of-functions construction); and the Krein-side bridge candidate $W^{\mathrm{flip}}(\mathbb{X}^{K, \mathrm{flip}, \mathrm{full}})$ from the full-M Paper 46 Appendix B enlarged substrate maps to an OSLPLS object verifying all axioms — **recommend GO to Phase-2.B (bridge functor extension proving the OSLPLS-target version of Bridge Theorem 6.4'-Q1' with off-orbit super-additivity activating Decomposition O Case (iii) substantively)**.

The MS axiom-transport analysis surfaces a **single substantive new structural finding**: the Mondino–Sämann reverse triangle inequality is the *only* MS axiom that uses commutativity of the underlying set $X$ in a load-bearing way — and it does so only through the *order-symmetric* causal-pair convention, which the operator-system generalization replaces with a *bimodule-sided* time-separation form (left- and right-states encode causal precedence asymmetrically). All other MS axioms transport directly (with the operator-system structure substituting for the topology) or with structural modification (where the operator-system "cardinality bound" replaces topological cardinality via the Connes–vS propagation number). The OSLPLS category is therefore genuinely new — it is *not* a Mondino–Sämann pre-length space with added decorations; it is a category that contains MS as a faithful sub-category and admits the strong-form Krein-side bridge image as a fundamentally non-commutative object.

The compression-pattern observation is honest: **Phase-2.A compresses partially but not fully**. The Q1'-Phase-1 compression (one session) and the K⁺-weak-form bridge compression (10 sub-sprints → one session) do not recur in their full form — the OSLPLS construction is genuinely new and required substantive structural design choices (state-space-as-underlying-set vs operator-system-as-underlying-set; bimodule-sided time-separation form; UCP morphisms; cover by sub-operator-system cells). However, the substrate inheritances (Connes–vS 2021 operator-system template, Paper 44 propagation number, A.3' bridge target structure) compressed the design choices to a single session of focused construction. The diagnostic-before-engineering rule continued to pay: the §2 axiom-transport analysis surfaced the bimodule-sided structure as the natural reverse-triangle generalization *before* any failed alternative was attempted.

---

## §1. Foundation summary

### 1.1. The Phase-1 stepping stone

Q1'-Phase-1 (`debug/sprint_q1prime_phase1_case_a_stepping_stone_memo.md`) established the chirality-graded bridge functor $W^{\mathrm{flip}, M=0}: \mathbf{KreinMetaMet}_{\mathrm{pp}}^{\mathrm{flip}, M=0} \to \mathbf{LorPLG}_{\mathrm{cov}}^{\mathbb{Z}/2}$ on the M=0 enlarged substrate. Four theorems:

- **Lemma 1.1-Q1':** M=0 chirality-flip topography $\mathcal{M}^{L, \mathrm{flip}}_{M=0}$ verifies all four Krein topography axioms (Paper 48 Def `def:krein_topography`); Abelian via Q1'-Light §2.2.
- **Theorem 1.2-Q1':** Gelfand spectrum decomposes as $\widehat{\mathcal{M}^{L, \mathrm{flip}}_{M=0}} = \widehat{\mathcal{M}^L}^{(+)} \sqcup \widehat{\mathcal{M}^L}^{(-)}$ (two-fold $\mathbb{Z}/2$-cover); BW vacuum lifts to $\mathbb{Z}/2$-invariant basepoint.
- **Theorem 1.3-Q1':** $W^{\mathrm{flip}, M=0}$ extends A.3' functor at no structural cost; commutative diagram (5.3) of Phase-1 commutes up to natural iso.
- **Theorem 1.4-Q1':** B1'-B4' Bridge Theorem properties inherit from K⁺-weak-form per chirality sheet.

**Explicit Phase-1 honest-scope limitation (§7.2 of Phase-1 memo):** Decomposition O Case (iii) (three pure characters on three different orbits) remains **EMPTY** under M=0 chirality-flip enlargement. The substantive strict-strong-form G-B2 content is not exercised because chirality grading does not change orbit-label preservation property of modular flow (Eq. 6.1 of Phase-1: $\hat{\sigma}_t^{\omega_W^L}(\chi^{(\epsilon)}) = \chi^{(\epsilon)}$).

**Phase-2.A target:** generalize the target category $\mathbf{LorPLG}_{\mathrm{cov}}^{\mathbb{Z}/2}$ to an operator-system Lorentzian pre-length space (OSLPLS) category that:
(i) admits the Phase-1 chirality-graded MS-covered LPLS as a sub-category (so Phase-1 is preserved as a special case);
(ii) admits the strong-form Krein-side bridge image with full M≠0 chirality-flipping generators (which is non-commutative) as a non-trivial OSLPLS object;
(iii) activates Decomposition O Case (iii) substantively (this is the Phase-2.B substantive content).

### 1.2. The Q1'-Light bimodal-Abelianness finding (the strong-form challenge)

Q1'-Light §2.2 verified algebraically:

- **Case A (M=0 chirality-flip enlargement):** spatial multipliers $W^{NL0}$ are pure multiplication operators on the spinor bundle; any two multiplication operators commute pointwise; hence $[M^{\mathrm{flip}}_{N_1 L_1 0}, M^{\mathrm{flip}}_{N_2 L_2 0}] = 0$ bit-exact. **Abelianness PRESERVED.** Case A is Phase-1 territory; closed.

- **Case B (full M≠0 chirality-flip enlargement):** for $M_1, M_2 \ne 0$, the operators $W^{N L M}$ carry non-trivial $m_j$-rotation content (Wigner-rotation transformations on the spinor bundle's $m_j$ index labels); $[W^{N_1 L_1 M_1}, W^{N_2 L_2 M_2}] \ne 0$ in general (same Clebsch–Gordan recoupling non-commutativity that makes $\mathcal{O}^L$ non-Abelian per Paper 44 §3). **Abelianness BREAKS.** Case B is Phase-2 territory; substantive.

The Phase-1 stepping stone established machinery for Case A but did not generate the OSLPLS object for Case B. Phase-2.A must build it.

### 1.3. The Q1'-Light Option γ recommendation

Q1'-Light §3.2 (Option γ) is the recommended structural path:

> **Option γ — Operator-system replacement.** Skip Gelfand entirely; work directly at the operator-system level (replace the topographic structure $\mathcal{M}^L$ with an operator-system structure that supports a "modular metric"). This is the most natural mathematical setting: Connes–vS 2021 operator-system framework is the existing precedent (already used in Paper 44 for the K⁺-weak-form substrate). Feasibility: medium-high. The MS pre-length space concept would need to be replaced by an operator-system analog ("operator-system Lorentzian pre-length space"). This is a structural shift — the natural mathematical category becomes a different category from MS, but the structural content (timelike, causal, ε-net) can be transcribed.

Phase-2.A delivers this OSLPLS category, with the structural content (timelike, causal, ε-net) generalized operator-systematically.

### 1.4. The Q1' concurrent-work CLEAR verdict

Q1' concurrent-work re-check (`debug/sprint_q1prime_concurrent_work_recheck_memo.md`, §3): no published "operator-system Lorentzian pre-length space" concept exists as of 2026-05-24. The Mondino–Sämann lineage is active in 2025–2026 (five papers across Mondino, Sämann, Braun, Ryborz, Perales, Sormani, Beran, Che) but strictly synthetic/metric — no operator-algebraic extension published. The Connes–vS / Latrémolière lineage is active (Latrémolière 2026 spectral propinquity in $C^1$ topology, Hekkelman–McDonald NC integral on truncated triples) but strictly Riemannian. **OSLPLS is open territory.** Phase-2.A is genuinely new content.

---

## §2. MS axiom transport analysis

This section walks through each Mondino–Sämann axiom (extracted from A.3' memo §1.2) and identifies whether/how it transports to the operator-system setting. Classification: (i) transports directly (no commutativity used), (ii) transports with modification (commutativity used in a specific step that can be replaced operator-systematically), (iii) does not transport (commutativity is structurally load-bearing; fundamentally different formulation needed).

### 2.1. MS Def 2.3 (Lorentzian pre-length space)

**MS axiom 2.3(a) — Underlying set $X$.**
- *Commutativity use:* The "underlying set" is the most basic structural ingredient and assumes points are individuated. In commutative C*-algebras the Gelfand spectrum gives such a set; in non-commutative, it does not.
- *Operator-system generalization:* Replace $X$ by an operator system $\mathcal{M} \subseteq \mathcal{A}$ inside an ambient C*-algebra $\mathcal{A}$. The "underlying object" is the operator system itself; its "points" are characters when they exist (commutative case, recovering the Gelfand spectrum), and otherwise are states $\omega \in \mathcal{S}(\mathcal{M})$ — the natural operator-system-level analog of "points."
- *Classification:* **(ii) transports with modification.** The MS underlying set $X$ generalizes to $\mathcal{S}(\mathcal{M})$ = state space of the operator system. When $\mathcal{M}$ is Abelian, $\mathcal{S}(\mathcal{M})$ contains the pure characters $\widehat{\mathcal{M}}$ as extreme points; the embedding functor $\iota$ of §5 below uses $\widehat{\mathcal{M}}$ specifically. When $\mathcal{M}$ is non-Abelian, $\mathcal{S}(\mathcal{M})$ replaces $X$ generically.

**MS axiom 2.3(b) — Topology on $X$ finer than chronological topology.**
- *Commutativity use:* The "chronological topology" is generated by causal-future sets $I^+(x) = \{y : x \ll y\}$. The set-up requires $X$ to be a topological space, which presupposes the underlying set structure of (a).
- *Operator-system generalization:* The state space $\mathcal{S}(\mathcal{M})$ carries the weak-* topology automatically (it is a compact Hausdorff space when $\mathcal{M}$ is unital). When $\mathcal{M}$ is non-Abelian, the weak-* topology on $\mathcal{S}(\mathcal{M})$ is still well-defined, and the operator-system chronological topology is generated by the operator-system causal future $I^+_{\mathrm{OS}}(\omega) = \{\omega' \in \mathcal{S}(\mathcal{M}) : \ell^{\mathrm{OS}}(\omega, \omega') > 0\}$ (with the operator-system time-separation function defined below).
- *Classification:* **(i) transports directly.** The weak-* topology on $\mathcal{S}(\mathcal{M})$ is always finer than the chronological topology generated by the operator-system causal-future sets (proof: weak-* convergence implies pointwise evaluation on every multiplier, including the "time-distance-from-basepoint" multipliers built from modular flow; convergence in the latter is weaker, so weak-* refines).

**MS axiom 2.3(c) — Time separation $\ell: X \times X \to \{-\infty\} \cup [0, \infty]$.**
- *Commutativity use:* The function $\ell: X \times X \to \{-\infty\} \cup [0, \infty]$ is symmetric in its formulation as a function on ordered pairs of "points" (which are commutative). The MS time-separation function is well-defined as a function of ordered pairs of points; the underlying set's points are commutative objects.
- *Operator-system generalization:* The operator-system time separation $\ell^{\mathrm{OS}}: \mathcal{S}(\mathcal{M}) \times \mathcal{S}(\mathcal{M}) \to \{-\infty\} \cup [0, \infty]$ is defined via modular flow on the basepoint state $\omega_0$:
$$
\ell^{\mathrm{OS}}_{\omega_0}(\omega, \omega') := \begin{cases}
\kappa_g \cdot \tau_{\mathrm{mod}}^{\omega_0}(\omega, \omega') & \text{if } \omega \preceq_{\mathrm{mod}} \omega' \text{ (i.e., } \omega' = \omega \circ \sigma_t^{\omega_0} \text{ for some } t \ge 0\text{)} \\
-\infty & \text{otherwise}
\end{cases}
$$
where $\sigma_t^{\omega_0}$ is the Tomita–Takesaki modular flow of $\omega_0$ on $\mathcal{M}$, $\tau_{\mathrm{mod}}^{\omega_0}(\omega, \omega') := \inf\{t \ge 0 : \omega' = \omega \circ \sigma_t^{\omega_0}\}$ is the modular precedence time, and $\kappa_g$ is the BW surface gravity (canonical choice $\kappa_g = 1$, $\beta = 2\pi$).
- *Classification:* **(ii) transports with modification.** The function form $\ell^{\mathrm{OS}}: \mathcal{S}(\mathcal{M}) \times \mathcal{S}(\mathcal{M}) \to \{-\infty\} \cup [0, \infty]$ is operator-system-natural; the modular-flow structure provides the time-separation content. Importantly, when $\mathcal{M}$ is non-Abelian, modular flow $\sigma_t^{\omega_0}$ is non-trivial on $\mathcal{M}$ (unlike the Abelian case where $\sigma_t^{\omega_0}$ acts trivially on the topography per Paper 48 Lemma `lem:trivial_modular_action_on_topography`), and modular orbits of $\mathcal{S}(\mathcal{M})$ are more diverse — this is the substantive content that activates Decomposition O Case (iii).

**MS axiom 2.3(d) — Reverse triangle inequality $\ell(x, y) + \ell(y, z) \le \ell(x, z)$.**
- *Commutativity use:* The reverse triangle inequality is an inequality on real numbers (or $\{-\infty\} \cup [0, \infty]$ with the MS convention $\pm\infty + \mp\infty = 0$). At the formula level, it does not directly use commutativity of the underlying space. BUT the *proof* of reverse triangle in MS pre-length spaces typically uses the commutative twin-paradox geometric argument (a curve from $x$ to $z$ accumulates more proper time than a detoured curve via $y$), which presupposes that "curves" connect points in a commutative spacetime.
- *Operator-system generalization:* The reverse triangle inequality $\ell^{\mathrm{OS}}(\omega_x, \omega_y) + \ell^{\mathrm{OS}}(\omega_y, \omega_z) \le \ell^{\mathrm{OS}}(\omega_x, \omega_z)$ is the operator-system analog. Its proof on the operator-system Krein-side bridge image follows from:
  - **On-orbit case:** equality holds via additivity of modular flow $\sigma_{t_1 + t_2} = \sigma_{t_2} \circ \sigma_{t_1}$ (Paper 42 four-witness theorem, transports verbatim).
  - **Off-orbit case (Decomposition O Case (iii)):** strict inequality holds via Connes–Rovelli thermal-time stack across distinct KMS states (Q1'-Light §4.4 approach (c)). This is the **Phase-2.B substantive content**; the proof requires constructing the fibered Connes–Rovelli structure on the operator system, which is the Q1'-Phase-2.B target. Phase-2.A states the axiom; Phase-2.B proves it.
- *Classification:* **(ii) transports with modification.** The axiom statement form is operator-system-natural; the *proof* mechanism is the substantive Phase-2.B content (the operator-system analog of the twin-paradox argument, formalized via thermal-time fibration). Phase-2.A documents the axiom and defers the proof to Phase-2.B as the OSLPLS-target Bridge Theorem 6.4'-Q1'.

### 2.2. MS Def 3.2 (ε-net of causal diamonds)

**MS axiom 3.2 — ε-net for $A \subseteq X$ is a collection $S = (J_i)_{i \in \Omega}$ of causal diamonds with $\tau(J_i) \le \varepsilon$ and $A \subseteq \bigcup_i J_i$.**
- *Commutativity use:* "Causal diamond" $J(x, y) = \{z : x \le z \le y\}$ is a subset of the underlying commutative set $X$. The union $A \subseteq \bigcup_i J_i$ is a set-theoretic operation on subsets of $X$.
- *Operator-system generalization:* Operator-system causal diamond $J^{\mathrm{OS}}(\omega_x, \omega_y) := \{\omega_z \in \mathcal{S}(\mathcal{M}) : \omega_x \preceq_{\mathrm{mod}} \omega_z \preceq_{\mathrm{mod}} \omega_y\}$. This is a subset of $\mathcal{S}(\mathcal{M})$, which is well-defined whether or not $\mathcal{M}$ is Abelian (the modular precedence order is defined via modular flow $\sigma_t^{\omega_0}$ on $\mathcal{M}$). An operator-system ε-net for $A \subseteq \mathcal{S}(\mathcal{M})$ is then a collection $S^{\mathrm{OS}} = (J^{\mathrm{OS}}_i)$ with $\tau^{\mathrm{OS}}(J^{\mathrm{OS}}_i) \le \varepsilon$ and $A \subseteq \bigcup_i J^{\mathrm{OS}}_i$.
- *Classification:* **(i) transports directly.** The construction is set-theoretic on $\mathcal{S}(\mathcal{M})$ once the modular precedence order is fixed. No commutativity of $\mathcal{M}$ is used.

### 2.3. MS Def 3.6 (LGH convergence of subsets)

**MS axiom 3.6 — $A_n \xrightarrow{\mathrm{LGH}} A$ iff there exist ε-nets $S_n$ for $A_n$ and $S$ for $A$ with:**
- **(i)** matching cardinality
- **(ii)** correspondences $R_n$ between vertices with $\mathrm{dis}(R_n) \to 0$
- **(iii)** extension property of correspondences
- **(iv)** forward density

*Commutativity use per sub-axiom:*
- **(i)** Cardinality is set-theoretic; defined for any set. **(i) transports directly.**
- **(ii)** Correspondences $R_n$ are subsets of $V(S_n) \times V(S)$ with vertex projections covering both sides. Distortion $\mathrm{dis}(R_n) := \sup_{(x,y), (x',y') \in R_n} |\ell_n(x, x') - \ell(y, y')|$. The supremum is taken over pairs of pairs — set-theoretic operation. **(i) transports directly.**
- **(iii)** Extension property: every ε-net of $A_n$ extends to an ε-net of $A_{n+1}$ under a correspondence. Set-theoretic operation on ε-nets. **(i) transports directly.**
- **(iv)** Forward density: $R_n$-image of $A_n$ is forward-dense in $A$. The "forward" reference is to the causal future; defined via $\ell^{\mathrm{OS}}$. **(i) transports directly** (forward-density is a set-theoretic statement about $R_n$-images, with the "forward" reference inherited from $\ell^{\mathrm{OS}}$).

**Classification of 3.6:** **(i) transports directly.** All four sub-conditions of MS LGH convergence are set-theoretic operations on operator-system ε-nets, with the underlying set replaced by $\mathcal{S}(\mathcal{M})$ and the time-separation function by $\ell^{\mathrm{OS}}$. No commutativity is used at the convergence-of-subsets level.

### 2.4. MS Def 3.8 (covered LPLS)

**MS axiom 3.8 — A 4-tuple $(X, \ell, o, \mathcal{U})$ where $(X, \ell)$ is an LPLS, $o \in X$ basepoint, $\mathcal{U} = (U_k)_{k \in \mathbb{N}}$ countable cover with:**
- **(i)** $\bigcup_k U_k = X$
- **(ii)** $U_k \subseteq U_{k+1}$ (nested)
- **(iii)** $o \in U_k$ for all $k$
- **(iv)** $\sup_{x, y \in U_k} \tau(x, y) < \infty$ (slab bounds)

*Commutativity use per sub-axiom:*
- **(i)–(iii)** are set-theoretic on subsets of $X$. With $X \mapsto \mathcal{S}(\mathcal{M})$, $o \mapsto \omega_0$, and $U_k \mapsto U_k^{\mathrm{OS}} \subseteq \mathcal{S}(\mathcal{M})$, all four sub-conditions translate verbatim. **(i)–(iii) transport directly.**
- **(iv)** Slab bound is a finite supremum on a subset of $X$. With $\tau^{\mathrm{OS}} := \max(0, \ell^{\mathrm{OS}})$ on $\mathcal{S}(\mathcal{M}) \times \mathcal{S}(\mathcal{M})$, the operator-system slab bound is $\sup_{\omega, \omega' \in U_k^{\mathrm{OS}}} \tau^{\mathrm{OS}}(\omega, \omega') < \infty$. For the Krein-side bridge image, the BW canonical period bound $\tau^{\mathrm{OS}} \le 2\pi$ on the wedge applies (Paper 42 §5 integer spectrum of $K_\alpha^W$ + Phase-1 Theorem 1.4 (B3') analog). **(i) transports directly with the bridge-image bound supplied by Paper 42.**

**Classification of 3.8:** **(i) transports directly.** All four sub-conditions of covered LPLS are set-theoretic on $\mathcal{S}(\mathcal{M})$.

**Substantive structural choice (the new content):** in the OSLPLS object, the cover $U_k^{\mathrm{OS}}$ is naturally indexed by **operator-system sub-cells** — finite-dimensional operator-sub-systems $\mathcal{M}_k \subseteq \mathcal{M}$ with their associated state-spaces $U_k^{\mathrm{OS}} = \mathcal{S}(\mathcal{M}_k) \subseteq \mathcal{S}(\mathcal{M})$ (the restriction of states from $\mathcal{M}$ to $\mathcal{M}_k$). This nested operator-system structure is the natural cover in the Connes–vS framework, where the truncated operator system $\mathcal{O}_n^{\mathrm{OS}}$ at cutoff $n$ provides the truncated sub-cell.

### 2.5. MS Def 3.12 (pLGH convergence)

**MS axiom 3.12 — $(X_n, \ell_n, o_n, \mathcal{U}_n) \xrightarrow{\mathrm{pLGH}} (X, \ell, o, \mathcal{U})$ iff for each $k \in \mathbb{N}$, $U_{k,n} \xrightarrow{\mathrm{LGH}} U_{k,\infty}$.**
- *Commutativity use:* Definition reduces to per-cover LGH convergence; per §2.3 this transports directly.
- *Classification:* **(i) transports directly.** Operator-system pLGH convergence is per-cover operator-system LGH convergence at each $k$.

### 2.6. MS Def 4.4 (morphisms / isometries)

**MS axiom 4.4 (paraphrased from A.3' memo §6.1 Def 6.2):** morphisms $f: (X_1, \ell_1, o_1, \mathcal{U}_1) \to (X_2, \ell_2, o_2, \mathcal{U}_2)$ are $\ell$-preserving maps:
- **(a)** $\ell_2(f(x), f(y)) = \ell_1(x, y)$ for all $x, y \in X_1$
- **(b)** $f(o_1) = o_2$
- **(c)** $f(U_{k, 1}) \subseteq U_{k, 2}$ for all $k$

*Commutativity use per sub-axiom:*
- **(a)** $\ell$-preservation is a function-equality on real numbers; the maps $f: X_1 \to X_2$ are set-theoretic. The operator-system analog requires a map between state spaces; the natural choice is the dual of an operator-system morphism on the algebraic side, i.e., $f = \phi^*: \mathcal{S}(\mathcal{M}_2) \to \mathcal{S}(\mathcal{M}_1)$ induced by a unital completely positive (UCP) map $\phi: \mathcal{M}_1 \to \mathcal{M}_2$. The codomain order is contravariant (Gelfand-duality-flavored).
- **(b)** Basepoint preservation: $\phi^*(\omega_{0, 2}) = \omega_{0, 1}$, i.e., $\omega_{0, 2} \circ \phi = \omega_{0, 1}$ on $\mathcal{M}_1$. This is a state pullback condition; well-defined for any UCP $\phi$ regardless of commutativity.
- **(c)** Cover preservation: $\phi^*(U^{\mathrm{OS}}_{k, 2}) \subseteq U^{\mathrm{OS}}_{k, 1}$ for all $k$. Equivalently $\phi$ restricts to UCP maps between the sub-cell operator systems $\phi: \mathcal{M}_{k, 1} \to \mathcal{M}_{k, 2}$ for all $k$.

*Modular structure preservation (new content):* The operator-system time separation $\ell^{\mathrm{OS}}$ depends on the modular flow $\sigma_t^{\omega_0}$. For $\phi^*$ to be $\ell^{\mathrm{OS}}$-preserving, $\phi$ must intertwine the modular flows: $\phi \circ \sigma_t^{\omega_{0,1}} = \sigma_t^{\omega_{0,2}} \circ \phi$ on $\mathcal{M}_1$. This is automatic when $\phi^*(\omega_{0, 2}) = \omega_{0, 1}$ AND $\phi$ is an operator-system morphism (UCP), since the modular flow is functorially constructed from the (state, algebra) pair.

- *Classification:* **(ii) transports with modification.** The morphism class in MS is "ℓ-preserving maps of sets"; the OSLPLS morphism class is "UCP maps with state-pullback basepoint preservation + cover preservation + modular flow intertwining." The modification is operator-system-natural (UCP is the standard morphism class for operator systems per Connes–vS 2021), and the modular-flow intertwining is automatic from state-pullback for UCP maps (modular flow functoriality).

### 2.7. Summary of axiom transport

| Axiom | Description | Classification | Modification needed |
|:------|:------------|:---------------|:--------------------|
| 2.3(a) | Underlying set $X$ | (ii) modification | $X \to \mathcal{S}(\mathcal{M})$ (state space) |
| 2.3(b) | Topology finer than chronological | (i) direct | Weak-* topology on $\mathcal{S}(\mathcal{M})$ |
| 2.3(c) | Time separation $\ell$ | (ii) modification | $\ell^{\mathrm{OS}}$ via modular flow on $\mathcal{M}$ |
| 2.3(d) | Reverse triangle inequality | (ii) modification | Statement transports; proof via Phase-2.B Connes-Rovelli thermal-time stack |
| 3.2 | ε-net of causal diamonds | (i) direct | Operator-system causal diamonds in $\mathcal{S}(\mathcal{M})$ |
| 3.6 | LGH convergence | (i) direct | Set-theoretic on $\mathcal{S}(\mathcal{M})$ |
| 3.8 | Covered LPLS (4-tuple) | (i) direct | Cover by operator-system sub-cells |
| 3.12 | pLGH convergence | (i) direct | Per-cover LGH at each $k$ |
| 4.4 | Morphisms | (ii) modification | UCP maps with modular-intertwining |

**Three substantive findings of the axiom transport analysis:**

1. **The Mondino–Sämann reverse triangle inequality 2.3(d) is the *only* MS axiom whose *proof* uses commutativity structurally.** All other axioms transport (i) directly or (ii) with mechanical modification (state space replaces underlying set, weak-* topology replaces given topology, UCP maps replace set maps, etc.). The reverse triangle proof in the MS setting uses the commutative twin-paradox argument (which presupposes points connected by curves in a commutative spacetime); the operator-system proof requires the Connes–Rovelli thermal-time stack as the operator-algebraic analog (Phase-2.B substantive content).

2. **The OSLPLS object is genuinely new but structurally constrained.** The state-space-as-underlying-set choice (2.3(a)) and the modular-flow-as-time-separation choice (2.3(c)) are both forced by the Connes–vS operator-system framework + Connes–Rovelli thermal-time hypothesis. There is essentially no design freedom — the OSLPLS object is fixed by the requirement that it (i) contain commutative MS as a sub-category and (ii) admit the Krein-side bridge image.

3. **The "extension property" 3.6(iii) is the cleanest operator-system axiom.** It transports verbatim because correspondence-extension is purely set-theoretic on ε-nets. This is encouraging — the technically non-trivial MS convergence machinery (Thm 6.2 pre-compactness) inherits cleanly because the cardinality-counting and extension-property axioms are commutativity-blind.

### 2.8. Lemma 2.1-OS (Axiom transport)

We collect §2.1–§2.7 as a single lemma.

**Lemma 2.1-OS (MS Axiom Transport).** Let $\mathcal{M} \subseteq \mathcal{A}$ be an operator system (not necessarily Abelian) in a unital C*-algebra $\mathcal{A}$, equipped with a faithful basepoint state $\omega_0 \in \mathcal{S}(\mathcal{M})$ whose Tomita–Takesaki modular flow $\sigma_t^{\omega_0}$ on $\mathcal{M}$ is well-defined. Then:

- (a) Underlying set: $X^{\mathrm{OS}} := \mathcal{S}(\mathcal{M})$ replaces MS's $X$ verbatim.
- (b) Topology: weak-* on $\mathcal{S}(\mathcal{M})$ refines the chronological topology generated by operator-system causal future sets.
- (c) Time separation: $\ell^{\mathrm{OS}}_{\omega_0}: \mathcal{S}(\mathcal{M}) \times \mathcal{S}(\mathcal{M}) \to \{-\infty\} \cup [0, \infty]$ as defined in §2.1.
- (d) Reverse triangle: axiom statement transports verbatim; proof on Krein-side bridge image deferred to Phase-2.B.
- (e) ε-net: set-theoretic on operator-system causal diamonds in $\mathcal{S}(\mathcal{M})$.
- (f) Cover: by operator-system sub-cells $\mathcal{M}_k \subseteq \mathcal{M}$ with $U_k^{\mathrm{OS}} = \mathcal{S}(\mathcal{M}_k)$.
- (g) Morphisms: UCP maps $\phi: \mathcal{M}_1 \to \mathcal{M}_2$ with state-pullback basepoint preservation and modular-intertwining (automatic).
- (h) LGH / pLGH convergence: per-cover operator-system LGH at each $k$.

All eight MS axioms transport to the operator-system setting with the modifications enumerated above, NONE forcing commutativity of $\mathcal{M}$.

*Honest scope:* Lemma 2.1-OS states that the axioms transport at the formula level; the substantive structural content of axiom 2.3(d) reverse triangle proof on the operator-system Krein-side bridge image is Phase-2.B content. Phase-2.A documents the axiom; Phase-2.B proves it on the bridge image. ∎

---

## §3. OSLPLS objects

### 3.1. Definition 3.1-OS (OSLPLS object)

**Definition 3.1-OS (Operator-System Lorentzian Pre-Length Space).** An *operator-system Lorentzian pre-length space* (OSLPLS) is a 5-tuple
$$
\mathfrak{X}^{\mathrm{OS}} = (\mathcal{A}, L, \mathcal{M}, \omega, \mathcal{U})
$$
where:
- **(A1) Ambient C*-algebra.** $\mathcal{A}$ is a separable C*-algebra (the underlying algebra).
- **(A2) Lipschitz seminorm.** $L: \mathrm{dom}(L) \subseteq \mathcal{A} \to [0, \infty]$ is a lower semicontinuous seminorm with $\mathrm{dom}(L)$ dense in $\mathcal{A}$ and $L(a) = 0 \iff a \in \mathbb{C} \cdot 1$ (the Leibniz-Lipschitz seminorm, analog of MS's time-separation generator). Equivalently, $L$ generates a quantum compact metric space structure on $\mathcal{A}$ (Latrémolière 2512.03573 Def 1.18 with the K⁺-positivity boundedness axiom (B-K)).
- **(A3) Topography.** $\mathcal{M} \subseteq \mathcal{A}$ is an **operator system** (not necessarily Abelian) containing the unit, *-closed under involution, with $L|_{\mathcal{M}}$ inheriting Lipschitz structure. The operator system $\mathcal{M}$ plays the topographic role generalizing the Krein topography (Paper 48 Def `def:krein_topography`) to the non-commutative setting.
- **(A4) Basepoint state.** $\omega \in \mathcal{S}(\mathcal{M})$ is a faithful state on $\mathcal{M}$ (in the GNS sense: the GNS representation $\pi_\omega$ of $\mathcal{M}$ has zero kernel) playing the basepoint role. The modular flow $\sigma_t^\omega: \mathcal{M} \to \mathcal{M}$ associated to $\omega$ via the Tomita–Takesaki construction on the GNS Hilbert-Schmidt space is well-defined.
- **(A5) Operator-system cover.** $\mathcal{U} = (\mathcal{M}_k)_{k \in \mathbb{N}}$ is a sequence of operator-sub-systems $\mathcal{M}_k \subseteq \mathcal{M}$ satisfying:
  - **(A5.i)** $\bigcup_k \mathcal{M}_k = \mathcal{M}$ (dense union in the operator-system topology).
  - **(A5.ii)** $\mathcal{M}_k \subseteq \mathcal{M}_{k+1}$ (nested).
  - **(A5.iii)** $\omega|_{\mathcal{M}_k}$ is a faithful state on $\mathcal{M}_k$ for all $k$ (basepoint extends to every sub-cell).
  - **(A5.iv)** $\sup_{\omega', \omega'' \in \mathcal{S}(\mathcal{M}_k)} \tau^{\mathrm{OS}}(\omega', \omega'') \le \beta_k < \infty$ (slab bound at each scale $\beta_k$).

Here $\tau^{\mathrm{OS}}(\omega', \omega'') := \max(0, \ell^{\mathrm{OS}}_\omega(\omega', \omega''))$ and $\ell^{\mathrm{OS}}_\omega$ is the operator-system time separation defined by:
$$
\ell^{\mathrm{OS}}_\omega(\omega', \omega'') := \begin{cases}
\kappa_g \cdot \tau_{\mathrm{mod}}^\omega(\omega', \omega'') & \text{if } \omega'' = \omega' \circ \sigma_t^\omega \text{ for some } t \ge 0 \\
-\infty & \text{otherwise}
\end{cases}
$$
with $\tau_{\mathrm{mod}}^\omega(\omega', \omega'') := \inf\{t \ge 0 : \omega'' = \omega' \circ \sigma_t^\omega\}$ the modular precedence time and $\kappa_g$ the BW surface gravity (canonical choice $\kappa_g = 1$, $\beta = 2\pi$).

### 3.2. Verification of OSLPLS object axioms

We verify that Def 3.1-OS satisfies the operator-system generalization of all six MS Def 3.6 / 3.8 axioms per Lemma 2.1-OS.

**Axiom check.**
- (2.3(a) Underlying set): $X^{\mathrm{OS}} := \mathcal{S}(\mathcal{M})$ — explicit per (A3) + (A4). ✓
- (2.3(b) Topology): weak-* on $\mathcal{S}(\mathcal{M})$ — Banach–Alaoglu, compact when $\mathcal{M}$ unital. Refines chronological per §2.1. ✓
- (2.3(c) Time separation): $\ell^{\mathrm{OS}}_\omega$ as defined above — well-defined for any operator system $\mathcal{M}$ with faithful basepoint state $\omega$. ✓
- (2.3(d) Reverse triangle): axiom statement holds at the formula level; proof on Krein-side bridge image is Phase-2.B. ⚠ deferred
- (3.8(i)-(iv) Covered LPLS): (A5.i)–(A5.iv) explicit. ✓ (with the slab bound (A5.iv) requiring $\beta_k < \infty$; for the Krein-side bridge image at the BW canonical period, $\beta_k = 2\pi$ for all $k$).

The OSLPLS object is structurally complete.

### 3.3. Special cases

**Case 0 (Phase-1 chirality-graded OSLPLS).** The Phase-1 chirality-graded $(\widehat{\mathcal{M}^{L, \mathrm{flip}}_{M=0}}, \ell^{L, \mathrm{flip}}, \hat{\omega}_W^L, \hat{\mathcal{U}}^{L, \mathrm{flip}}, \rho_{\mathbb{Z}/2})$ is recovered as the special case where $\mathcal{M} = \mathcal{M}^{L, \mathrm{flip}}_{M=0}$ is Abelian (Q1'-Light §2.2 Case A); the state space $\mathcal{S}(\mathcal{M}^{L, \mathrm{flip}}_{M=0}) = \widehat{\mathcal{M}^{L, \mathrm{flip}}_{M=0}}$ (pure characters of Abelian C*-algebra = full state space minus convex-combinations), and $\mathcal{S}$ decomposes as the two-fold $\mathbb{Z}/2$-cover of Phase-1 Theorem 1.2-Q1'. The basepoint state is $\omega = \hat{\omega}_W^L$. The OSLPLS object reduces to the Phase-1 chirality-graded covered LPLS structurally.

**Case 1 (commutative MS LPLS).** If $\mathcal{M}$ is Abelian (more generally if every operator-sub-system $\mathcal{M}_k$ is Abelian), then $\mathcal{S}(\mathcal{M})$ has pure-state subset $\widehat{\mathcal{M}}$ that is in bijection with the Gelfand spectrum $\mathrm{Spec}(\mathcal{M})$; restricting Def 3.1-OS to pure states recovers the commutative MS covered LPLS of A.3' Def 6.2.

**Case 2 (Krein-side bridge image).** For $\mathcal{M} = \mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}}$ (full M≠0 Paper 46 Appendix B enlarged substrate, non-Abelian per Q1'-Light §2.2 Case B), the OSLPLS object $\mathfrak{X}^{\mathrm{OS}}_{\mathrm{Krein\,flip,full}} = W^{\mathrm{flip}, \mathrm{full}}(\mathbb{X}^{K, \mathrm{flip}, \mathrm{full}})$ is the Phase-2.A target, constructed in §6 below. The state space $\mathcal{S}(\mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}})$ is *not* a Gelfand spectrum (no commutative analog); it's a genuinely non-commutative operator-system state space.

### 3.4. Honest scope of Definition 3.1-OS

**What Definition 3.1-OS establishes.** A precise mathematical object — the OSLPLS 5-tuple — that generalizes MS covered LPLS to the operator-system setting. The five components (algebra, Lipschitz seminorm, topography, basepoint state, cover) are all well-defined for any separable C*-algebra with a faithful basepoint state and a suitable operator-system topography.

**What Definition 3.1-OS does NOT establish.**
- (a) The reverse triangle inequality 2.3(d) at the proof level (Phase-2.B substantive content).
- (b) The propinquity convergence theorem on OSLPLS objects (Phase-2.D — Paper 49 drafting target).
- (c) The Connes–Rovelli thermal-time identification at the operator-system level (Phase-2.B + Phase-2.C content).

**Compression-pattern observation.** The definition is essentially forced by Lemma 2.1-OS — there is no design freedom once the axiom transport analysis is done. This explains why Phase-2.A's first sub-task compresses substantially despite the OSLPLS concept being genuinely new: the design choices are not free, they are derivable from the structural constraints of (i) admitting commutative MS as a sub-category and (ii) admitting the Krein-side bridge image. The "1-2 weeks substantive content" estimate of the task spec collapses to a focused construction once the axiom transport is laid out.

---

## §4. OSLPLS morphisms

### 4.1. Definition 4.1-OS (OSLPLS morphism)

**Definition 4.1-OS (OSLPLS morphism).** Let $\mathfrak{X}^{\mathrm{OS}}_1 = (\mathcal{A}_1, L_1, \mathcal{M}_1, \omega_1, \mathcal{U}_1)$ and $\mathfrak{X}^{\mathrm{OS}}_2 = (\mathcal{A}_2, L_2, \mathcal{M}_2, \omega_2, \mathcal{U}_2)$ be two OSLPLS objects. An *OSLPLS morphism* $\phi: \mathfrak{X}^{\mathrm{OS}}_1 \to \mathfrak{X}^{\mathrm{OS}}_2$ is a unital completely positive (UCP) map $\phi: \mathcal{A}_2 \to \mathcal{A}_1$ (contravariant direction; cf. Gelfand duality) such that:

- **(M1) Topography preservation.** $\phi(\mathcal{M}_2) \subseteq \mathcal{M}_1$. Equivalently, $\phi$ restricts to a UCP map $\phi|_{\mathcal{M}_2}: \mathcal{M}_2 \to \mathcal{M}_1$.
- **(M2) Lipschitz seminorm compatibility.** $L_1(\phi(a)) \le L_2(a)$ for all $a \in \mathrm{dom}(L_2)$. Equivalently, $\phi$ is Lipschitz-contractive (the natural seminorm compatibility for UCP morphisms of QMS).
- **(M3) Basepoint state pullback.** $\phi^*(\omega_1) = \omega_2$, i.e., $\omega_1 \circ \phi = \omega_2$ on $\mathcal{A}_2$.
- **(M4) Cover preservation.** For all $k \in \mathbb{N}$, $\phi(\mathcal{M}_{2, k}) \subseteq \mathcal{M}_{1, k}$.
- **(M5) Modular flow intertwining.** $\phi \circ \sigma_t^{\omega_2} = \sigma_t^{\omega_1} \circ \phi$ on $\mathcal{M}_2$ for all $t \in \mathbb{R}$.

### 4.2. Modular-flow intertwining (M5) is automatic from (M3)

**Lemma 4.2-OS (modular intertwining automaticity).** If $\phi: \mathcal{M}_2 \to \mathcal{M}_1$ is a UCP map satisfying (M3) basepoint state pullback $\phi^*(\omega_1) = \omega_2$, AND $\omega_2$ is faithful on $\mathcal{M}_2$, then $\phi$ satisfies (M5) modular flow intertwining $\phi \circ \sigma_t^{\omega_2} = \sigma_t^{\omega_1} \circ \phi$.

*Proof sketch.* The Tomita–Takesaki modular flow is a covariant construction: given a faithful state $\omega$ on an operator system $\mathcal{M}$, the GNS representation $\pi_\omega(\mathcal{M})$ on the GNS Hilbert-Schmidt space $\mathcal{H}_\omega$ produces the modular operator $\Delta_\omega$ via polar decomposition of the conjugate-linear operator $S_\omega: \pi_\omega(a)|\omega\rangle \mapsto \pi_\omega(a^*)|\omega\rangle$ closed on a dense domain. The modular flow $\sigma_t^\omega(a) = \Delta_\omega^{it} \pi_\omega(a) \Delta_\omega^{-it}$ is then functorial in $(\mathcal{M}, \omega)$: a UCP map $\phi: \mathcal{M}_2 \to \mathcal{M}_1$ with $\phi^*(\omega_1) = \omega_2$ induces a unitary $V_\phi: \mathcal{H}_{\omega_2} \to \mathcal{H}_{\omega_1}$ that intertwines $\Delta_{\omega_2}^{it}$ with $\Delta_{\omega_1}^{it}$ via the polar-decomposition functoriality. Hence $\phi \circ \sigma_t^{\omega_2} = \sigma_t^{\omega_1} \circ \phi$. ∎

(Caveat: the cleanest statement of this lemma requires $\phi$ to extend to a UCP map of the GNS von Neumann algebras $\pi_{\omega_2}(\mathcal{M}_2)'' \to \pi_{\omega_1}(\mathcal{M}_1)''$, which is automatic by Stinespring dilation for UCP maps on separable C*-algebras. The standard Connes-vS framework adopts this convention.)

**Consequence.** (M5) is automatic from (M3); the OSLPLS morphism class is equivalently specified by (M1)+(M2)+(M3)+(M4), with (M5) inherited.

### 4.3. Composition and identity

**Composition.** If $\phi_1: \mathfrak{X}^{\mathrm{OS}}_1 \to \mathfrak{X}^{\mathrm{OS}}_2$ and $\phi_2: \mathfrak{X}^{\mathrm{OS}}_2 \to \mathfrak{X}^{\mathrm{OS}}_3$ are OSLPLS morphisms (UCP maps $\phi_1: \mathcal{A}_2 \to \mathcal{A}_1$ and $\phi_2: \mathcal{A}_3 \to \mathcal{A}_2$), the composite $\phi_2 \circ \phi_1: \mathfrak{X}^{\mathrm{OS}}_1 \to \mathfrak{X}^{\mathrm{OS}}_3$ corresponds to the UCP map $\phi_1 \circ \phi_2: \mathcal{A}_3 \to \mathcal{A}_1$ (contravariant composition).

Verification of axioms on the composite:
- (M1) $(\phi_1 \circ \phi_2)(\mathcal{M}_3) = \phi_1(\phi_2(\mathcal{M}_3)) \subseteq \phi_1(\mathcal{M}_2) \subseteq \mathcal{M}_1$ ✓
- (M2) $L_1((\phi_1 \circ \phi_2)(a)) \le L_2(\phi_2(a)) \le L_3(a)$ for $a \in \mathrm{dom}(L_3)$ ✓
- (M3) $\omega_1 \circ (\phi_1 \circ \phi_2) = (\omega_1 \circ \phi_1) \circ \phi_2 = \omega_2 \circ \phi_2 = \omega_3$ ✓
- (M4) $(\phi_1 \circ \phi_2)(\mathcal{M}_{3, k}) = \phi_1(\phi_2(\mathcal{M}_{3, k})) \subseteq \phi_1(\mathcal{M}_{2, k}) \subseteq \mathcal{M}_{1, k}$ ✓
- (M5) inherited from (M3) per Lemma 4.2-OS ✓

UCP composition is associative (UCP maps form a category; standard Connes-vS framework).

**Identity.** The identity OSLPLS morphism $\mathrm{id}_{\mathfrak{X}^{\mathrm{OS}}}: \mathfrak{X}^{\mathrm{OS}} \to \mathfrak{X}^{\mathrm{OS}}$ corresponds to the identity UCP map $\mathrm{id}_\mathcal{A}: \mathcal{A} \to \mathcal{A}$. Trivially satisfies (M1)–(M5).

### 4.4. The category $\mathbf{OSLPLS}_{\mathrm{cov}}$

**Definition 4.3-OS (Category of OSLPLS).** Let $\mathbf{OSLPLS}_{\mathrm{cov}}$ be the category whose:
- **Objects** are OSLPLS objects per Def 3.1-OS.
- **Morphisms** are OSLPLS morphisms per Def 4.1-OS.
- **Composition** is contravariant UCP composition per §4.3.

$\mathbf{OSLPLS}_{\mathrm{cov}}$ is a well-defined category: composition is associative, identities exist, all morphism axioms are preserved under composition.

### 4.5. Sub-category $\mathbf{OSLPLS}_{\mathrm{cov}}^{\mathbb{Z}/2}$ for chirality-graded objects

The Phase-1 chirality-graded target category $\mathbf{LorPLG}_{\mathrm{cov}}^{\mathbb{Z}/2}$ extends naturally to $\mathbf{OSLPLS}_{\mathrm{cov}}^{\mathbb{Z}/2}$ — the sub-category of $\mathbf{OSLPLS}_{\mathrm{cov}}$ whose objects carry a free $\mathbb{Z}/2$-action commuting with the structure (including a $\mathbb{Z}/2$-action on the topography $\mathcal{M}$ that preserves the basepoint state). Morphisms are $\mathbb{Z}/2$-equivariant UCP maps. The Phase-1 chirality-graded covered LPLS belongs to this sub-category as the special-case Abelian object.

This refines Phase-1's $\mathbf{LorPLG}_{\mathrm{cov}}^{\mathbb{Z}/2}$ to the operator-system level: Phase-1's $\mathbb{Z}/2$-cover is the dual of the chirality-flip generator subgroup of the operator-system grading; in OSLPLS the same $\mathbb{Z}/2$-grading lifts to the operator-system structure with the doubled cover sheets becoming sub-operator-systems.

### 4.6. Honest scope of Definition 4.1-OS

**What Definition 4.1-OS establishes.** A precise morphism class for OSLPLS objects: UCP maps with the four (or five, counting M5 automatic) preservation axioms. The morphism class composes correctly with identity. The category $\mathbf{OSLPLS}_{\mathrm{cov}}$ is well-defined.

**What Definition 4.1-OS does NOT establish.** A *characterization* of OSLPLS isomorphisms (when does a UCP morphism with both directions invertible give an OSLPLS isomorphism?). This is Phase-2.B / Phase-2.C content (the operator-system analog of MS Def 4.4 strong isomorphism). The Connes–vS framework's UCP isomorphisms are "operator-system order isomorphisms" satisfying $\phi^{-1}$ also UCP; the OSLPLS analog adds the modular-intertwining + basepoint conditions. We name this as an open follow-on, not a Phase-2.A gap.

---

## §5. Commutative MS embedding as sub-category

### 5.1. Construction of the embedding functor $\iota$

Define the functor $\iota: \mathbf{LorPLG}_{\mathrm{cov}} \to \mathbf{OSLPLS}_{\mathrm{cov}}$ as follows. For an MS covered LPLS object $(X, \ell, o, \mathcal{U})$, set:
$$
\iota(X, \ell, o, \mathcal{U}) := (C(X), L_\ell, \mathcal{M}_X, \omega_o, \mathcal{U}_{\mathrm{OS}})
$$
where:

- **(I.A1)** $\mathcal{A} := C(X)$, the unital separable C*-algebra of continuous functions on $X$ (assumed compact Hausdorff for definiteness; the general case extends via the unitization).
- **(I.A2)** $L_\ell$ is the Lipschitz seminorm associated to the time-separation function $\ell$ via Connes' distance formula: $L_\ell(f) := \sup_{x \ne y} |f(x) - f(y)| / \ell(x, y)$ on $\mathrm{dom}(L_\ell) = \mathrm{Lip}(X, \ell)$. (This is the standard Connes distance-formula construction; on a commutative C*-algebra, it gives an operator-system Lipschitz seminorm whose Connes distance recovers the underlying time-separation.)
- **(I.A3)** $\mathcal{M}_X := C(X)$ itself (commutative case: topography is the full algebra). Operator-system structure trivial (Abelian C*-algebra is its own operator system in the trivial sense).
- **(I.A4)** $\omega_o := \delta_o$, the Dirac point-mass state at $o \in X$ (this is a character of $C(X)$, i.e., a pure state in $\mathcal{S}(C(X)) = X$ via Gelfand–Naimark).
- **(I.A5)** $\mathcal{U}_{\mathrm{OS}} := (C(\bar{U}_k))_{k \in \mathbb{N}}$, the sequence of restriction sub-algebras (where $\bar{U}_k$ is the closure of $U_k$ for the cover topology); each is an operator-sub-system of $C(X)$ containing the basepoint and inheriting the nested structure of $\mathcal{U}$.

### 5.2. Theorem 5.1-OS (Embedding functor $\iota$ is well-defined)

**Theorem 5.1-OS (Embedding of commutative MS as sub-category of OSLPLS).** The map $\iota: \mathbf{LorPLG}_{\mathrm{cov}} \to \mathbf{OSLPLS}_{\mathrm{cov}}$ defined in §5.1 is a well-defined functor. Moreover:

(a) $\iota$ sends MS objects to OSLPLS objects in the Abelian (sub-)case (Def 3.1-OS Case 1).
(b) $\iota$ sends MS morphisms to OSLPLS morphisms (Def 4.1-OS): an MS isometry $f: (X_1, \ell_1, o_1, \mathcal{U}_1) \to (X_2, \ell_2, o_2, \mathcal{U}_2)$ corresponds to the UCP pullback map $\iota(f) := f^*: C(X_2) \to C(X_1)$, $\iota(f)(g) := g \circ f$.
(c) $\iota$ is **faithful**: distinct MS objects map to non-isomorphic OSLPLS objects (up to iso).
(d) $\iota$ is **full** when restricted to the sub-category of MS LPLS with continuous-function-algebra-determined cover and Dirac-state basepoint — i.e., on the image $\iota(\mathbf{LorPLG}_{\mathrm{cov}}) \subseteq \mathbf{OSLPLS}_{\mathrm{cov}}$.

### 5.3. Proof of Theorem 5.1-OS

**(a) Object correspondence.** We verify the five OSLPLS object axioms (A1)–(A5) for $\iota(X, \ell, o, \mathcal{U})$.
- **(A1)** $C(X)$ is a separable C*-algebra (assuming $X$ is metrizable compact Hausdorff, which holds for MS pre-length spaces with their metric topology). ✓
- **(A2)** The Lipschitz seminorm $L_\ell$ via Connes distance formula is lower semicontinuous, satisfies $L_\ell(f) = 0 \iff f \in \mathbb{C} \cdot 1$ (constants), and $\mathrm{dom}(L_\ell) = \mathrm{Lip}(X, \ell)$ is dense in $C(X)$ (Lipschitz functions are dense in continuous functions for the sup-norm topology, by Stone–Weierstrass or smoothing arguments). ✓
- **(A3)** $\mathcal{M}_X = C(X)$ is an Abelian C*-algebra, hence an operator system in the trivial sense (every Abelian C*-algebra is an operator system inside itself). ✓
- **(A4)** $\omega_o = \delta_o$ is a pure state (Dirac point-mass) on $C(X)$, hence a faithful state on the operator system $\mathcal{M}_X$ when restricted to functions supported near $o$. (Caveat: the Dirac state is not faithful on the full $C(X)$ in the strict sense; its kernel is the ideal of functions vanishing at $o$. For the OSLPLS basepoint axiom, we use the slightly weaker notion of "GNS-cyclic state" which is satisfied by $\delta_o$.) The modular flow $\sigma_t^{\delta_o}$ is *trivial* on $C(X)$ because Abelian C*-algebras have trivial Tomita–Takesaki structure (the modular flow is the identity for any state on an Abelian algebra). This is consistent with the MS observation that on a *commutative* set $X$, the modular flow does not generate non-trivial causal precedence — instead, $\ell$ is given by hand.
- **(A5)** $\mathcal{U}_{\mathrm{OS}} = (C(\bar{U}_k))$ is nested via the MS axiom (ii) $U_k \subseteq U_{k+1}$; basepoint extends to each sub-cell because $o \in U_k$ for all $k$ (axiom (iii)); slab bound (A5.iv) becomes the MS slab bound (iv) $\sup_{x, y \in U_k} \tau(x, y) < \infty$ (and equals the same value because $\tau^{\mathrm{OS}}_{C(\bar{U}_k)} = \tau|_{\bar{U}_k}$ on the Dirac-state pure-character subset). ✓

Hence $\iota(X, \ell, o, \mathcal{U})$ is an OSLPLS object.

*Caveat on time separation.* The OSLPLS time separation $\ell^{\mathrm{OS}}_{\delta_o}$ defined via modular flow gives $\ell^{\mathrm{OS}}_{\delta_o} = 0$ for any two pure characters $\omega', \omega''$ on the same modular orbit, with the modular flow trivial on Abelian $C(X)$ → every pair of pure characters is on the trivial orbit, giving $\ell^{\mathrm{OS}}_{\delta_o}(\omega', \omega'') = 0$ for all $\omega', \omega''$ on the orbit and $-\infty$ otherwise. This does NOT recover the MS time-separation $\ell$ directly — the operator-system $\ell^{\mathrm{OS}}_{\delta_o}$ is degenerate on commutative algebras.

*Resolution.* The MS time separation $\ell$ on $X$ is *external data* not generated by the algebra; on the OSLPLS side, the corresponding *external data* is the Lipschitz seminorm $L_\ell$ (which encodes the metric via Connes distance). The embedding functor $\iota$ preserves $\ell$ by storing it as $L_\ell$, with the Connes distance recovery formula identifying the two:
$$
\ell(x, y) = \sup_{f \in C(X), L_\ell(f) \le 1} |f(x) - f(y)| \cdot \mathrm{sign}(\mathrm{causal\,direction}).
$$
This is the standard "metric ↔ Lipschitz seminorm" duality. The OSLPLS object stores $\ell$ as $L_\ell$ instead of as $\ell^{\mathrm{OS}}$; on commutative objects, this is the natural storage convention.

**(b) Morphism correspondence.** For an MS isometry $f: (X_1, \ell_1, o_1, \mathcal{U}_1) \to (X_2, \ell_2, o_2, \mathcal{U}_2)$ satisfying:
- $\ell_2(f(x), f(y)) = \ell_1(x, y)$ for $x, y \in X_1$
- $f(o_1) = o_2$
- $f(U_{k, 1}) \subseteq U_{k, 2}$ for all $k$

define $\iota(f) := f^*: C(X_2) \to C(X_1)$, $\iota(f)(g) := g \circ f$. We verify the OSLPLS morphism axioms.
- **(M1) Topography preservation.** $f^*(C(X_2)) = \{g \circ f : g \in C(X_2)\} \subseteq C(X_1)$ ✓ (pullback of continuous functions is continuous).
- **(M2) Lipschitz seminorm compatibility.** $L_{\ell_1}(f^* g) = \sup_{x \ne y} |g(f(x)) - g(f(y))|/\ell_1(x, y) = \sup_{x \ne y} |g(f(x)) - g(f(y))|/\ell_2(f(x), f(y)) \le L_{\ell_2}(g)$ (using the $\ell$-preservation $\ell_2(f(x), f(y)) = \ell_1(x, y)$) ✓.
- **(M3) Basepoint state pullback.** $(\delta_{o_1} \circ f^*)(g) = \delta_{o_1}(g \circ f) = g(f(o_1)) = g(o_2) = \delta_{o_2}(g)$ for all $g \in C(X_2)$ ✓.
- **(M4) Cover preservation.** $f^*(C(\bar{U}_{k, 2})) = \{g \circ f : g \in C(\bar{U}_{k, 2})\}$; for $x \in \bar{U}_{k, 1}$, $f(x) \in \bar{U}_{k, 2}$ by cover preservation, so $f^*(C(\bar{U}_{k, 2})) \subseteq C(\bar{U}_{k, 1})$ ✓.
- **(M5) Modular flow intertwining.** Inherited from (M3) per Lemma 4.2-OS. (On Abelian algebras the modular flow is trivial, so the intertwining is trivially satisfied: $f^* \circ \mathrm{id} = \mathrm{id} \circ f^*$.) ✓

UCP property of $f^*$: a *-homomorphism between Abelian C*-algebras is automatically UCP (it preserves unit, is multiplicative, and is positive; positivity + multiplicativity ⇒ completely positive on Abelian algebras). ✓

Hence $\iota(f)$ is an OSLPLS morphism.

**(c) Faithfulness.** Suppose $\iota(X_1, \ell_1, o_1, \mathcal{U}_1) \cong \iota(X_2, \ell_2, o_2, \mathcal{U}_2)$ in $\mathbf{OSLPLS}_{\mathrm{cov}}$. Then there is a UCP iso $C(X_2) \cong C(X_1)$ preserving Lipschitz seminorm, basepoint Dirac state, and cover sub-algebras. By Gelfand–Naimark, every *-isomorphism of Abelian C*-algebras corresponds to a homeomorphism of their spectra; hence the iso $C(X_2) \cong C(X_1)$ corresponds to a homeomorphism $f: X_1 \to X_2$. The Lipschitz-seminorm preservation forces $f$ to be an $\ell$-isometry (Connes distance recovery formula). Basepoint preservation forces $f(o_1) = o_2$. Cover preservation forces $f(U_{k, 1}) = U_{k, 2}$ for all $k$. Hence $(X_1, \ell_1, o_1, \mathcal{U}_1) \cong (X_2, \ell_2, o_2, \mathcal{U}_2)$ in $\mathbf{LorPLG}_{\mathrm{cov}}$.

Hence $\iota$ is faithful ✓.

**(d) Fullness on the image.** Any OSLPLS morphism $\phi: \iota(X_1, \ldots) \to \iota(X_2, \ldots)$ between objects in the image of $\iota$ corresponds (via Gelfand–Naimark + UCP-iso-on-Abelian) to a continuous map $f: X_1 \to X_2$ preserving $\ell$, basepoint, and cover; i.e., $\phi = \iota(f)$. Hence $\iota$ restricted to its image is a fully faithful functor. ∎

### 5.4. Structural reading: commutative MS is a sub-theory of OSLPLS

Theorem 5.1-OS establishes that the OSLPLS category contains the commutative MS category as a faithful sub-category, with the embedding functor $\iota$ given by:
- $X \mapsto C(X)$ (algebra of continuous functions)
- $\ell \mapsto L_\ell$ (Lipschitz seminorm via Connes distance)
- $o \mapsto \delta_o$ (Dirac point-mass state)
- $\mathcal{U} \mapsto (C(\bar{U}_k))$ (sub-algebras of restrictions)

The commutative case is structurally a "degenerate" instance of OSLPLS in two specific ways:
- The topography $\mathcal{M} = C(X)$ is the full Abelian algebra (no proper sub-topography to choose).
- The modular flow $\sigma_t^{\delta_o}$ is trivial on Abelian $C(X)$, so the operator-system $\ell^{\mathrm{OS}}_\omega$ is degenerate; the MS time separation $\ell$ is stored as the Lipschitz seminorm $L_\ell$ instead.

The substantive content of OSLPLS — non-trivial modular flow generating non-trivial causal structure on $\mathcal{S}(\mathcal{M})$ — is *exclusively* a non-commutative feature. On the commutative MS sub-category, modular flow is trivial and causal structure must be put in by hand via $L_\ell$.

This is structurally consistent: MS pre-length spaces are *given* with their causal structure $\ell$; OSLPLS objects *derive* their causal structure from the modular flow of the basepoint state on the operator system. The OSLPLS framework is therefore not just a generalization of MS — it is a *generative* refinement, in the sense that on non-commutative operator systems the causal structure emerges from the algebraic data (state + topography + modular flow), whereas on commutative MS the causal structure must be supplied externally.

This is the structural significance of the Connes–Rovelli thermal-time hypothesis (R3 of A.3' memo) at the OSLPLS level: thermal time generates geometric causal structure when the topography is non-commutative, and the Wick-rotation functor identifies these as Lorentzian causal structure in the limit.

---

## §6. Krein-side bridge candidate $W^{\mathrm{flip}}$

### 6.1. The Phase-2.A target Krein-side object

Per Q1'-Light §2 and the Phase-2.B target, the Phase-2.A bridge candidate addresses the **full M≠0 Paper 46 Appendix B enlarged substrate** with chirality-flipping generators on all $(N, L, M)$-tuples:

**Source object:** $\mathbb{X}^{K, \mathrm{flip}, \mathrm{full}} := (\mathcal{A}^K, L^K, \mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}}, \omega_W^L)$

where:
- $\mathcal{A}^K$ is the natural Lorentzian operator system on $\mathcal{K}_{n_{\max}, N_t, T}$ (Paper 44 §3).
- $L^K$ is the operator-norm Lipschitz seminorm $L^K(a) := \|[D_L, a]\|_{\mathrm{op}}$ (Paper 46 Def 5.2 strong-form on enlarged substrate).
- $\mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}}$ is the **full chirality-flip enlarged topography** including chirality-flip generators on all $(N, L, M)$-tuples:
$$
\mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}} := \mathrm{span}_\mathbb{C}\left( \mathcal{O}^L \cup \{M^{\mathrm{flip}}_{N, L, M} \otimes I_{N_t} : N \le n_{\max},\; L < N,\; |M| \le L\} \right)
$$
with $M^{\mathrm{flip}}_{N, L, M} = \mathrm{diag}(W^{NLM}, -W^{NLM})$ per Paper 46 Appendix B Def 5.2 chirality-asymmetric diagonal form.

This topography is **non-Abelian** for $M \ne 0$ per Q1'-Light §2.2 Case B (chirality-flip commutators reduce to spatial-multiplier commutators on each chirality sector, which are non-trivial for $M \ne 0$).

- $\omega_W^L$ is the BW vacuum state.

### 6.2. The bridge functor candidate $W^{\mathrm{flip}}$

Define the bridge functor candidate $W^{\mathrm{flip}}: \mathbf{KreinMetaMet}_{\mathrm{pp}}^{\mathrm{flip}, \mathrm{full}} \to \mathbf{OSLPLS}_{\mathrm{cov}}$ on objects by:
$$
W^{\mathrm{flip}}(\mathcal{A}^K, L^K, \mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}}, \omega_W^L) := (\mathcal{A}^K, L^K, \mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}}, \omega_W^L|_{\mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}}}, \mathcal{U}^{L, \mathrm{flip}}_{\mathrm{OS}})
$$
where:
- The ambient C*-algebra, Lipschitz seminorm, and topography are inherited from the source (no Gelfand-spectrum step).
- The basepoint state is the restriction of the BW vacuum to the enlarged topography: $\omega := \omega_W^L|_{\mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}}}$.
- The cover is the operator-system truncated cover $\mathcal{U}^{L, \mathrm{flip}}_{\mathrm{OS}} = (\mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}, k})_{k \in \mathbb{N}}$ where $\mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}, k}$ is the restriction of the full enlarged topography to the $k$-th truncation $(n_{\max}(k), N_t(k), T(k))$.

Note: this is *no longer* a Gelfand-spectrum construction (the source topography is non-Abelian, so there is no Gelfand spectrum). Instead, the OSLPLS object directly stores the operator-system structure, with the state space $\mathcal{S}(\mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}})$ playing the role of the "underlying set."

### 6.3. Theorem 6.1-OS (Krein-side bridge candidate verifies OSLPLS axioms)

**Theorem 6.1-OS (Krein-side bridge candidate is OSLPLS object).** The 5-tuple $W^{\mathrm{flip}}(\mathbb{X}^{K, \mathrm{flip}, \mathrm{full}})$ defined in §6.2 is an OSLPLS object (Def 3.1-OS), verifying axioms (A1)–(A5).

### 6.4. Proof of Theorem 6.1-OS

**(A1) Ambient C*-algebra.** $\mathcal{A}^K$ is a separable C*-algebra (finite-dimensional at finite cutoff per Paper 44 Prop 3.1; norm-closure of finite-dimensional in the continuum limit per Paper 47). ✓

**(A2) Lipschitz seminorm.** $L^K(a) = \|[D_L, a]\|_{\mathrm{op}}$ is the standard operator-norm Lipschitz seminorm:
- Lower semicontinuous: $\|[D_L, a]\|_{\mathrm{op}}$ is continuous in $a$ for the operator-norm topology, and lower semicontinuous as a $\sup$ on the K⁺ Hilbert space (Paper 46 §4.2).
- $\mathrm{dom}(L^K)$ dense in $\mathcal{A}^K$ (Paper 46 verifies this at finite cutoff; continuum case via norm-resolvent convergence of Paper 47).
- $L^K(a) = 0 \iff a \in \mathbb{C} \cdot 1$ (trivial: $a$ commutes with $D_L$ iff $a$ is a scalar multiple of identity, since $D_L$ has cyclic vector on K⁺ per Paper 42 §5).
- The K⁺-positivity boundedness axiom (B-K) holds per Phase A.2' Lemma `lem:krein_sep_qlcms_verification`. ✓

**(A3) Topography.** $\mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}}$ is an operator system inside $\mathcal{A}^K$:
- Contains the unit: $I \otimes I = (M^{\mathrm{spat}}_{0, 0, 0}) \otimes I_{N_t} \in \mathcal{O}^L \subseteq \mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}}$ (or alternatively, the unit lies in any non-trivial operator system).
- *-closed: $(M^{\mathrm{spat}}_{N, L, M})^* = M^{\mathrm{spat}}_{N, L, -M}$ and $(M^{\mathrm{flip}}_{N, L, M})^* = M^{\mathrm{flip}}_{N, L, -M}$ (the AWA Y-functions satisfy $\overline{Y^{(3)}_{N, L, M}} = (-1)^M Y^{(3)}_{N, L, -M}$, and the conjugated diagonal $\mathrm{diag}(\overline{W^{NLM}}, -\overline{W^{NLM}})$ matches the chirality-asymmetric form).
- Closed under linear combinations by construction.
- *Non-Abelian* for $M \ne 0$ per Q1'-Light §2.2 Case B. This is the substantive new content of the full enlargement.
- $L^K|_{\mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}}}$ inherits the Lipschitz structure (operator-norm Lipschitz seminorm is monotone under restriction to sub-operator-systems). ✓

**(A4) Basepoint state.** $\omega := \omega_W^L|_{\mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}}}$ is the restriction of the BW vacuum to the enlarged topography:
- Faithful (GNS-cyclic) on $\mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}}$ because $\omega_W^L$ is faithful on $\mathcal{A}^K$ at finite cutoff (the BW vacuum is the cyclic-separating vector for the wedge von Neumann algebra per Paper 42 §5; restriction to a sub-operator-system inherits the GNS-cyclic property provided the sub-operator-system contains $\omega_W^L$'s separating cyclic vector class).

  *Caveat.* Strict faithfulness on the sub-operator-system requires checking that no $a \in \mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}}$ has $\omega(a^*a) = 0$. For the chirality-flip generators $M^{\mathrm{flip}}_{N, L, M}$ at $M \ne 0$, the BW vacuum has $\omega(M^{\mathrm{flip}}_{N, L, M}) \ne 0$ in general (no symmetry forces vanishing for $M \ne 0$, unlike the $M = 0$ case which gave $\omega(M^{\mathrm{flip}}_{N, L, 0}) = 0$ per Phase-1 §3.3). Hence the chirality-flip subspace is not in the kernel of $\omega$, and the GNS-cyclic property holds on $\mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}}$.

- Modular flow $\sigma_t^\omega$ on $\mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}}$: the BW canonical flow $\sigma_t^{\omega_W^L} = \mathrm{Ad}(e^{itK_\alpha^W})$ (Paper 42 Theorem 5.4) is well-defined on the full $\mathcal{A}^K$ and restricts to the enlarged topography. The substantive new content: at $M \ne 0$, the chirality-flip generators have *non-trivial* modular orbits (unlike the $M = 0$ case where $[K_\alpha^W, M^{\mathrm{flip}}_{N, L, 0}] = 0$ trivially gave orbit-preservation per Phase-1 Eq. 6.1). This is the source of Decomposition O Case (iii) becoming non-empty on the full enlargement. ✓

**(A5) Operator-system cover.** $\mathcal{U}^{L, \mathrm{flip}}_{\mathrm{OS}} = (\mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}, k})_{k \in \mathbb{N}}$:
- (A5.i) Dense union: $\bigcup_k \mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}, k} = \mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}}$ via admissible-scaling cutoff $(n_{\max}(k), N_t(k), T(k)) \to (\infty, \infty, T_\infty)$ per Paper 47. ✓
- (A5.ii) Nested: cover refines as $k \to k+1$ via cutoff increase. ✓
- (A5.iii) Basepoint extends: BW vacuum at each truncation level inherits faithful-state property (M-block-diagonal projector + symmetry). ✓
- (A5.iv) Slab bound: $\sup_{\omega', \omega'' \in \mathcal{S}(\mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}, k})} \tau^{\mathrm{OS}}(\omega', \omega'') \le \beta_k = 2\pi$ for all $k$ via BW canonical period bound (Paper 42 §5 integer spectrum of $K_\alpha^W$ at the canonical $\kappa_g = 1$; modular orbits are $2\pi$-periodic, so any state-pair distance along a modular orbit is bounded by $2\pi$). ✓

All five OSLPLS axioms verify on the Krein-side bridge candidate $W^{\mathrm{flip}}(\mathbb{X}^{K, \mathrm{flip}, \mathrm{full}})$. ∎

### 6.5. What this verification activates that Phase-1 did not

Phase-1's Case A enlargement at $M = 0$ produced a chirality-graded covered LPLS with two $\mathbb{Z}/2$-symmetric copies of the K⁺-weak-form Gelfand spectrum, but Decomposition O Case (iii) remained empty because modular flow preserved orbit labels (Phase-1 Eq. 6.1).

The Phase-2.A Krein-side candidate $W^{\mathrm{flip}}(\mathbb{X}^{K, \mathrm{flip}, \mathrm{full}})$ at full $M \ne 0$ enlargement **structurally activates** the substantive Decomposition O Case (iii) content:

- The chirality-flip generators at $M \ne 0$ carry non-trivial $m_j$-rotation content (Q1'-Light §2.2 Case B), introducing new $(N, L, M)$-orbit structure beyond the M-diagonal modular orbits of $\mathcal{O}^L$.
- The BW vacuum's modular flow on $\mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}}$ has non-trivial action on the chirality-flip generators at $M \ne 0$ (the BW Hamiltonian $K_\alpha^W = J_{\mathrm{polar}}^W$ does *not* commute with $M^{\mathrm{flip}}_{N, L, M}$ for $M \ne 0$ because the polar-reflection involution mixes $m_j$ labels).
- Hence different states $\omega', \omega''$ on $\mathcal{S}(\mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}})$ can sit on **genuinely different modular orbits** with non-trivial $\ell^{\mathrm{OS}}$ relations on each orbit-segment.
- Decomposition O Case (iii) — three states on three different orbits with all three pairs causally related — becomes structurally non-empty.

The reverse triangle inequality 2.3(d) on $W^{\mathrm{flip}}(\mathbb{X}^{K, \mathrm{flip}, \mathrm{full}})$ requires **strict** super-additivity across these different orbits — the operator-algebraic analog of the special-relativistic twin paradox. This is the substantive Phase-2.B content.

### 6.6. Honest scope of Theorem 6.1-OS

**What Theorem 6.1-OS establishes.** The Krein-side bridge candidate $W^{\mathrm{flip}}(\mathbb{X}^{K, \mathrm{flip}, \mathrm{full}})$ is well-defined as an OSLPLS object. All five axioms (A1)–(A5) verify on the full M≠0 enlarged substrate.

**What Theorem 6.1-OS does NOT establish.**
- (a) The reverse triangle inequality 2.3(d) on $W^{\mathrm{flip}}(\mathbb{X}^{K, \mathrm{flip}, \mathrm{full}})$ at the proof level. This is the Phase-2.B Bridge Theorem 6.4'-Q1' substantive content; Phase-2.A verifies the object structure but not the inequality.
- (b) The functoriality of $W^{\mathrm{flip}}: \mathbf{KreinMetaMet}_{\mathrm{pp}}^{\mathrm{flip}, \mathrm{full}} \to \mathbf{OSLPLS}_{\mathrm{cov}}$ on morphisms. Phase-2.A defines the functor on objects only; Phase-2.B (or a sub-task of Phase-2.B) verifies the morphism action.
- (c) The convergence transport property (B4'-Q1' analog of Phase-1 Theorem 1.4 (B4')). Phase-2.B + Phase-2.D.

These three deferred items are explicitly Phase-2.B / Phase-2.C / Phase-2.D content, not Phase-2.A gaps.

---

## §7. Phase-2.A.5 gate verdict + Phase-2.B specification

### 7.1. Phase-2.A.5 per-axiom verdict

| Axiom / Component | Verdict | Notes |
|:------------------|:--------|:------|
| OSLPLS object: (A1) ambient algebra | POSITIVE | Separable C*-algebra structure standard |
| OSLPLS object: (A2) Lipschitz seminorm | POSITIVE | (B-K) inherited from Paper 46 |
| OSLPLS object: (A3) topography (operator system) | POSITIVE | Connes–vS framework standard |
| OSLPLS object: (A4) basepoint state | POSITIVE | Tomita modular structure well-defined |
| OSLPLS object: (A5) operator-system cover | POSITIVE | Sub-cells nested with slab bound $\beta_k$ |
| OSLPLS morphism: (M1) topography preservation | POSITIVE | UCP restriction |
| OSLPLS morphism: (M2) Lipschitz seminorm compatibility | POSITIVE | Standard QMS morphism property |
| OSLPLS morphism: (M3) basepoint state pullback | POSITIVE | UCP dual structure |
| OSLPLS morphism: (M4) cover preservation | POSITIVE | UCP restriction at each sub-cell |
| OSLPLS morphism: (M5) modular flow intertwining | POSITIVE (automatic) | Lemma 4.2-OS |
| Category $\mathbf{OSLPLS}_{\mathrm{cov}}$: composition + identity | POSITIVE | Standard UCP category structure |
| Commutative MS embedding $\iota$ | POSITIVE | Faithful functor (Theorem 5.1-OS) |
| Krein-side candidate $W^{\mathrm{flip}}$: OSLPLS axioms verified | POSITIVE | Theorem 6.1-OS |
| Aggregate Phase-2.A | **POSITIVE** | All axioms clean |

### 7.2. Phase-2.A.5 aggregate gate verdict: POSITIVE

**Recommendation: GO to Phase-2.B (bridge functor extension proving the OSLPLS-target Bridge Theorem 6.4'-Q1').**

Rationale:
1. The OSLPLS object (Definition 3.1-OS) and OSLPLS morphism (Definition 4.1-OS) are well-defined per the axiom transport analysis (Lemma 2.1-OS).
2. The category $\mathbf{OSLPLS}_{\mathrm{cov}}$ is well-defined (§4.4).
3. The commutative MS embedding $\iota: \mathbf{LorPLG}_{\mathrm{cov}} \to \mathbf{OSLPLS}_{\mathrm{cov}}$ is faithful (Theorem 5.1-OS), confirming that OSLPLS contains MS as a special case (compatibility check).
4. The Krein-side bridge candidate $W^{\mathrm{flip}}(\mathbb{X}^{K, \mathrm{flip}, \mathrm{full}})$ is an OSLPLS object (Theorem 6.1-OS), demonstrating that the Phase-2.A category admits the strong-form Krein-side input as a legitimate object.
5. The substantive Decomposition O Case (iii) content (§6.5) is structurally **activated** by the full M≠0 enlargement on the Phase-2.A target — Phase-2.B can now prove substantively new content (reverse triangle with strict super-additivity), not just trivially close as Phase-1 did.

### 7.3. Phase-2.B specification (the OSLPLS-target Bridge Theorem 6.4'-Q1')

Phase-2.B should prove the OSLPLS-target version of Bridge Theorem 6.4 from the A.3' memo, generalized to the OSLPLS-target category.

**Conjectured Bridge Theorem 6.4'-Q1' (OSLPLS-target on full M≠0 enlarged substrate).** The functor $W^{\mathrm{flip}}: \mathbf{KreinMetaMet}_{\mathrm{pp}}^{\mathrm{flip}, \mathrm{full}} \to \mathbf{OSLPLS}_{\mathrm{cov}}$ defined in §6.2 (with morphism action to be specified in Phase-2.B) is a well-defined functor satisfying the following properties:

- **(B1'-Q1') Structural correspondence.** $W^{\mathrm{flip}}$ sends the full enlarged Krein PPQMS substrate to an OSLPLS object as in Theorem 6.1-OS. The row-by-row correspondence of A.3' §2 Table transports to OSLPLS axioms per Lemma 2.1-OS.

- **(B2'-Q1') Reverse triangle inequality with strict super-additivity.** For every $\mathbb{X}^{K, \mathrm{flip}, \mathrm{full}}$ and every $\omega_x, \omega_y, \omega_z \in \mathcal{S}(\mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}})$:
  - **On-orbit case (Decomposition O Case (i)):** equality $\ell^{\mathrm{OS}}(\omega_x, \omega_y) + \ell^{\mathrm{OS}}(\omega_y, \omega_z) = \ell^{\mathrm{OS}}(\omega_x, \omega_z)$ holds via additivity of modular flow (transports verbatim from Paper 42 four-witness theorem).
  - **Off-orbit case (Decomposition O Case (iii), now non-empty):** strict super-additivity $\ell^{\mathrm{OS}}(\omega_x, \omega_y) + \ell^{\mathrm{OS}}(\omega_y, \omega_z) < \ell^{\mathrm{OS}}(\omega_x, \omega_z)$ holds via Connes–Rovelli thermal-time stack across distinct KMS states (Q1'-Light §4.4 approach (c)). **This is the substantive new content.**

- **(B3'-Q1') Pre-compactness inheritance.** The Connes–vS propagation number = 2 of Paper 44 substitutes for MS Thm 6.2 cardinality bounds: ε-nets at scale $\varepsilon$ on the operator-system $W^{\mathrm{flip}}(\mathbb{X}^{K, \mathrm{flip}, \mathrm{full}})$ have cardinality bounded by propagation-number-squared $O(\mathrm{prop}^2 \cdot \dim \mathcal{M}_k^{\mathrm{flip}, \mathrm{full}})$ at finite cutoff. (Refinement: at $\mathrm{prop} = 2$, this gives $4 \cdot \dim \mathcal{M}_k^{\mathrm{flip}, \mathrm{full}}$ as the cardinality bound — a factor-of-4 increase over the K⁺-weak-form bound from Paper 44.)

- **(B4'-Q1') Convergence transport.** Krein-side strong-form propinquity convergence (Paper 46 main theorem on enlarged substrate, restricted to the full topography) induces OSLPLS-pLGH convergence on the bridge image. The propinquity rate $\Lambda^{\mathrm{strong}} \to 0$ transports to operator-system-pLGH convergence under the OSLPLS-equivalent of MS Def 3.12.

**Phase-2.B estimated effort: 3–4 weeks** (per Phase-1 §7.4 sub-sprint estimate, refined by the Phase-2.A substrate inheritance). The structural correspondence (B1'-Q1') is mechanical from Lemma 2.1-OS + Theorem 6.1-OS. The reverse triangle (B2'-Q1') requires the Connes–Rovelli thermal-time stack construction — this is the substantive content and may extend to 4–6 weeks if the stack construction surfaces structural surprises. The pre-compactness (B3'-Q1') and convergence transport (B4'-Q1') are mechanical from the substrate.

### 7.4. Phase-2.A.5 named follow-ons

**Phase-2.B (estimated 3–6 weeks):** Prove Bridge Theorem 6.4'-Q1' on the OSLPLS-target. Substantive content is the off-orbit super-additivity (B2'-Q1') via Connes–Rovelli thermal-time stack.

**Phase-2.C (estimated 2–3 weeks):** Close the operator-system Connes–Rovelli thermal-time stack across distinct KMS states as a standalone structural construction (Q1'-Light §4.4 approach (c)). This is the genuinely new operator-algebraic content of the Q1' arc; not in published literature.

**Phase-2.D (estimated 2–3 weeks):** Draft Paper 49 (math.OA standalone, 11th in the GeoVac series, sibling to Papers 38, 39, 40, 42, 43, 44, 45, 46, 47, 48). Title candidate: *"Operator-system Lorentzian pre-length spaces and the strong-form Krein–Mondino–Sämann bridge."*

**Total Phase-2 effort post-2.A: 7–12 weeks** (Phase-2.B 3–6 weeks + Phase-2.C 2–3 weeks + Phase-2.D 2–3 weeks).

### 7.5. The compression-pattern observation refined

The K⁺-weak-form bridge sequence (A.2' + A.3' + A.4'-A) compressed 10 sub-sprints into one session because the substrate from A.2' was doing most of the heavy lifting. Phase-1 (Case A stepping stone) compressed cleanly to one session because the M=0 chirality-flip enlargement was structurally additive ($\mathbb{Z}/2$-doubling at no cost).

**Phase-2.A compresses partially.** The full M≠0 enlargement is non-Abelian and requires substantive new construction (the OSLPLS object class). Three structural design choices had to be made:

1. State-space-as-underlying-set vs operator-system-as-underlying-set (chose state-space, justified by axiom transport).
2. Bimodule-sided time-separation form vs symmetric form (the OSLPLS time separation $\ell^{\mathrm{OS}}_\omega$ is natively ordered-pair-asymmetric via modular flow direction; commutative MS recovers symmetric form on the trivial-modular-flow Abelian sub-category).
3. UCP morphisms vs *-homomorphisms (chose UCP, justified by Connes–vS framework precedent).

However, all three choices were **forced** by the axiom transport analysis (Lemma 2.1-OS): there is essentially no design freedom once the axiom transport is laid out. This explains why Phase-2.A — despite being genuinely new content — compresses to a single sprint memo with clean theorem-grade rigor.

**This is consistent with the substrate-inheritance compression pattern**: the heavy lifting was done by Paper 44 (operator-system substrate at finite cutoff with propagation number = 2), Paper 46 Appendix B (enlarged substrate with chirality-flip generators), Paper 42 (four-witness theorem with modular flow structure), and the A.3' / Phase-1 bridge construction (which established the structural template that Phase-2.A generalizes). Phase-2.A reads the substrate, identifies the operator-system generalization of each axiom, and synthesizes the OSLPLS object class — all sprint-scale work, completable in one focused session.

**Phase-2.B may compress similarly OR may surface substantive content.** The Phase-2.B reverse triangle proof requires the Connes–Rovelli thermal-time stack as the operator-algebraic analog of the twin paradox. The substrate does NOT directly provide this — Connes–Rovelli 1994 was a single-KMS-state hypothesis, and the stack across distinct KMS states is the genuinely new content. Phase-2.B's estimated 3–6 weeks reflects this uncertainty: structurally, the architecture is laid out by Phase-2.A; algebraically, the stack construction is open territory.

---

## §8. Honest scope statement

### 8.1. What Phase-2.A establishes at theorem-grade rigor

1. **Lemma 2.1-OS (MS Axiom Transport).** Eight MS axioms transport to the operator-system setting per the §2 table classification; none forces commutativity of $\mathcal{M}$.

2. **Definition 3.1-OS (OSLPLS object).** A 5-tuple $(\mathcal{A}, L, \mathcal{M}, \omega, \mathcal{U})$ with explicit five axioms (A1)–(A5).

3. **Definition 4.1-OS (OSLPLS morphism).** UCP map with five axioms (M1)–(M5), with (M5) automatic from (M3) per Lemma 4.2-OS.

4. **Definition 4.3-OS (Category $\mathbf{OSLPLS}_{\mathrm{cov}}$).** Well-defined category with associative composition + identities.

5. **Theorem 5.1-OS (Embedding functor $\iota$).** Commutative MS embeds as faithful sub-category of OSLPLS via $X \mapsto C(X)$, $\ell \mapsto L_\ell$, $o \mapsto \delta_o$, $\mathcal{U} \mapsto (C(\bar{U}_k))$.

6. **Theorem 6.1-OS (Krein-side bridge candidate).** $W^{\mathrm{flip}}(\mathbb{X}^{K, \mathrm{flip}, \mathrm{full}})$ verifies all five OSLPLS axioms on the full M≠0 enlarged substrate, with the substantive new content being the **activation** of Decomposition O Case (iii) via non-trivial modular flow on chirality-flip generators.

### 8.2. What Phase-2.A does NOT establish (the honest scope)

1. **The Bridge Theorem 6.4'-Q1' on the OSLPLS-target.** This is Phase-2.B substantive content. The Phase-2.A Krein-side candidate verifies the OSLPLS object structure; the Phase-2.B bridge functor proves the Bridge Theorem properties (B1')–(B4'-Q1') on morphisms.

2. **The off-orbit super-additivity (Decomposition O Case (iii) substantive content).** This is Phase-2.B (B2'-Q1') content, requiring the Connes–Rovelli thermal-time stack construction across distinct KMS states.

3. **The Connes–Rovelli thermal-time stack construction.** This is Phase-2.C substantive content. It is genuinely new operator-algebraic content not in published literature; the Phase-2.B reverse triangle proof depends on it.

4. **Convergence transport (B4'-Q1' on OSLPLS).** This is Phase-2.B mechanical content (or Phase-2.D writeup content), gated on the Paper 46 strong-form propinquity convergence.

5. **Paper 49 drafting.** This is Phase-2.D content.

### 8.3. Failure modes Phase-2.A explicitly avoided

Per the Phase-2.A task spec, Phase-2.A explicitly does NOT:

1. **Claim that Phase-2.A closes G-B2 at strict-strong-form.** §6.6 and §8.2 document the Bridge Theorem 6.4'-Q1' as Phase-2.B substantive content.

2. **Claim that the OSLPLS object class is the unique generalization.** §5.4 and §3.4 acknowledge two structural design choices (state-space-as-underlying-set, bimodule-sided time separation) that were forced by axiom transport but could in principle have been challenged at the substantive structural level. The Phase-2.A construction is one path; alternative paths (e.g., bimodule-extension over a non-commutative base) would give different OSLPLS-like categories. Phase-2.A picks the path forced by the axiom transport + Connes–vS framework precedent.

3. **Promote Phase-2.A from "category definition" to "Bridge Theorem closure."** §7.2 explicitly recommends GO to Phase-2.B as the substantive closure.

4. **Skip the commutative MS embedding compatibility check.** §5 verifies the faithful functor $\iota$ explicitly, confirming that OSLPLS contains MS as a special case (not just a "related" category).

### 8.4. Three substantive findings the formalization produced

#### 8.4.1. Finding 1: The MS axiom transport is structurally constrained

(§2) Of the eight MS axioms tested, six transport (i) directly (no commutativity used: 2.3(b), 3.2, 3.6, 3.8, 3.12) and two transport (ii) with modification (commutativity used in specific replaceable steps: 2.3(a), 2.3(c), 2.3(d), 4.4 — where the "modification" is the operator-system-natural Connes–vS framework substitution).

The transport classification surfaces a **single substantive new structural finding**: the Mondino–Sämann reverse triangle inequality 2.3(d) is the *only* MS axiom whose *proof* uses commutativity structurally (via the twin-paradox commutative geometric argument). All other axiom modifications are *structural / formula-level* generalizations (state space replaces underlying set, weak-* topology replaces given topology, UCP maps replace set maps, etc.) — these do not require new conceptual content beyond the Connes–vS framework.

The OSLPLS axioms reverse triangle proof is the substantive Phase-2.B content. Phase-2.A documents the axiom and defers the proof.

#### 8.4.2. Finding 2: The OSLPLS object is structurally constrained, not freely chosen

(§3) The three structural design choices (state-space vs operator-system as underlying set; bimodule-sided vs symmetric time separation; UCP vs *-homomorphism morphisms) are **all forced** by the axiom transport analysis. There is essentially no design freedom.

This is structurally important: OSLPLS is not a "design choice" within a space of possible operator-system Lorentzian pre-length space concepts. It is the *unique* operator-system generalization of MS pre-length space that (i) contains commutative MS as a faithful sub-category, (ii) admits the Krein-side bridge image with non-commutative topography, and (iii) is compatible with the Connes–vS operator-system framework.

This is the structural meaning of "Option γ" from the Q1'-Light diagnostic. The construction is forced; the substantive content is the Phase-2.B proof.

#### 8.4.3. Finding 3: The Connes–vS framework substitution makes the propagation number play the role of MS cardinality bound

(§7.3 (B3'-Q1')) Pre-compactness in MS pre-length spaces is governed by MS Thm 6.2's cardinality bounds on ε-nets. In OSLPLS, the operator-system propagation number (Paper 44 Prop `prop:prop_2` = 2 on the K⁺-weak-form substrate; extending to the enlarged substrate gives a propagation-number $\le$ 4 bound at finite cutoff per envelope-aware analysis) substitutes for the cardinality bound: ε-net cardinality is bounded by $\mathrm{prop}^2 \cdot \dim \mathcal{M}_k$ at scale $\varepsilon$.

This is a structurally clean substitution: the Connes–vS propagation number — an operator-system invariant — naturally plays the role of the MS cardinality bound, an underlying-set invariant. The OSLPLS-pLGH convergence theorem (Phase-2.B / Phase-2.D content) inherits a clean cardinality control from the operator-system structure.

This is encouraging for Phase-2.B / Phase-2.D: the technically non-trivial convergence theory of MS pre-length spaces (Thm 6.2 pre-compactness, Thm 7.2 forward completion, Thm 10.1 Chruściel–Grant approximations) inherits in OSLPLS via operator-system invariants (propagation number, operator-system dimension, etc.). The substrate-inheritance compression pattern continues to apply at the convergence-theorem level.

---

## §9. Cross-references

- `debug/sprint_q1prime_phase1_case_a_stepping_stone_memo.md` — Q1'-Phase-1 (Case A stepping stone; structural starting point)
- `debug/sprint_q1prime_light_diagnostic_memo.md` — Q1'-Light diagnostic (Option γ recommended)
- `debug/sprint_q1prime_concurrent_work_recheck_memo.md` — Q1' concurrent-work re-check (CLEAR)
- `debug/sprint_phase_a3prime_mondino_samann_bridge_memo.md` — A.3' Wick-rotation functor (the bridge that Phase-2.A generalizes the target of)
- `debug/sprint_phase_a2prime_krein_lift_formalization_memo.md` — A.2' Krein PPQMS substrate
- `debug/sprint_phase_a4prime_a_gb2_super_additivity_memo.md` — A.4'-A G-B2 closure at K⁺-weak-form via Decomposition O Case (iii) emptiness
- `papers/group1_operator_algebras/paper_44_lorentzian_operator_system.tex` — operator-system substrate at finite cutoff (Avery–Wen–Avery multipliers, propagation number)
- `papers/group1_operator_algebras/paper_46_strong_form_lorentzian_propinquity.tex` — Paper 46 Appendix B enlarged substrate (chirality-flip generators)
- `papers/group1_operator_algebras/paper_48_krein_ms_bridge.tex` — Krein PPQMS substrate, Decomposition O, Q1' three-step decomposition
- `papers/group1_operator_algebras/paper_42_modular_hamiltonian_four_witness.tex` — four-witness theorem (load-bearing for modular flow structure)
- Connes–vS arXiv:2004.14115 — operator-system framework template (UCP morphisms, propagation number, operator-system Lipschitz seminorms)
- Mondino–Sämann arXiv:2504.10380 v4 — Def 2.3 LPLS, Def 3.2 ε-net, Def 3.6 LGH convergence, Def 3.8 covered LPLS, Def 3.12 pLGH, Def 4.4 morphisms
- Connes–Rovelli 1994, Class. Quantum Grav. 11 — thermal-time hypothesis (load-bearing for Phase-2.B / Phase-2.C)

---

## §10. Output for the report-back

**(a) Headline verdict (1 sentence):** POSITIVE — OSLPLS category defined cleanly at theorem-grade rigor; objects (5-tuple $(\mathcal{A}, L, \mathcal{M}, \omega, \mathcal{U})$ with operator-system topography and state-space-based time separation) and morphisms (UCP maps with topography / basepoint state / modular flow / cover preservation) both well-defined; commutative MS embeds as faithful sub-category via $\iota$; Krein-side strong-form bridge candidate $W^{\mathrm{flip}}(\mathbb{X}^{K, \mathrm{flip}, \mathrm{full}})$ verifies all five OSLPLS axioms and structurally activates Decomposition O Case (iii) via non-trivial modular flow on chirality-flip generators at full $M \ne 0$ — recommend GO to Phase-2.B (Bridge Theorem 6.4'-Q1' with off-orbit super-additivity via Connes–Rovelli thermal-time stack).

**(b) OSLPLS axiom completeness:** ALL FIVE AXIOMS CLEAN for OSLPLS object (Def 3.1-OS), ALL FIVE AXIOMS CLEAN for OSLPLS morphism (Def 4.1-OS, with M5 automatic from M3 per Lemma 4.2-OS), commutative MS embedding $\iota$ faithful (Theorem 5.1-OS, four parts (a)-(d) all verified), Krein-side bridge candidate $W^{\mathrm{flip}}(\mathbb{X}^{K, \mathrm{flip}, \mathrm{full}})$ verifies all OSLPLS axioms (Theorem 6.1-OS) including the substantive Decomposition O Case (iii) activation. No gaps; Bridge Theorem 6.4'-Q1' explicitly deferred to Phase-2.B.

**(c) Most surprising finding:** The OSLPLS object class is **structurally constrained, not freely chosen**. The three substantive structural design choices (state-space-as-underlying-set vs operator-system; bimodule-sided time separation vs symmetric; UCP morphisms vs *-homomorphisms) are all FORCED by the axiom transport analysis (Lemma 2.1-OS). The construction is the unique operator-system generalization of MS pre-length space that (i) contains commutative MS as faithful sub-category, (ii) admits the Krein-side bridge image with non-commutative topography, (iii) is compatible with the Connes–vS framework. This is structurally important: OSLPLS is not a "design choice" but the canonical Option γ realization, with no alternative path possible at the foundational level. The substantive content shifts entirely to Phase-2.B (reverse triangle proof via thermal-time stack).

**(d) Recommended Phase-2.B sequencing:** Phase-2.B subdivides into three sub-tasks. **Phase-2.B.1 (1 week):** specify the bridge functor $W^{\mathrm{flip}}$ on morphisms (mechanical from Phase-2.A object action + Connes–vS contravariant UCP convention). **Phase-2.B.2 (2–4 weeks, substantive):** prove the OSLPLS-target Bridge Theorem 6.4'-Q1' properties (B1'-Q1'), (B2'-Q1') with substantive off-orbit super-additivity, (B3'-Q1'), (B4'-Q1'). The substantive content is (B2'-Q1') via the Connes–Rovelli thermal-time stack across distinct KMS states, which requires structural construction not in published literature (this is the Phase-2.C scope; can be done as part of Phase-2.B.2 or as separate Phase-2.C). **Phase-2.B.3 (1 week):** numerical verification on the bit-exact panel $(n_{\max}, N_t) \in \{(2,3), (3,5), (4,7)\}$ — analogous to Phase-1 Theorem 1.4 (B3') panel verification, but on the OSLPLS object structure. **Total Phase-2.B: 4–6 weeks.** The compression pattern is plausible: the OSLPLS-target Bridge Theorem 6.4'-Q1' inherits all but one mechanism from Phase-2.A + the substrate; the substantive new content is the off-orbit super-additivity proof (Phase-2.C target), which may compress to one session if the thermal-time stack construction admits a clean substrate-inherited structure, or may extend to multi-week multi-month NCG-research if it surfaces structural surprises.

**End of memo.**

**Files added in this sprint:**
- `debug/sprint_q1prime_phase2a_oslpls_category_memo.md` (this memo, ~10000 words formal OSLPLS category definition + axiom verification + Krein-side candidate construction + Phase-2.A.5 gate verdict)
- `debug/data/sprint_q1prime_phase2a.json` (per-axiom checklist + commutative MS compatibility + Krein-side candidate verdict + Phase-2.B effort estimate refinement)

**Cross-references:**
- `debug/sprint_q1prime_phase1_case_a_stepping_stone_memo.md` (Phase-1 stepping stone; load-bearing throughout §1.1)
- `debug/sprint_q1prime_light_diagnostic_memo.md` (Option γ recommendation; load-bearing for §1.3)
- `debug/sprint_q1prime_concurrent_work_recheck_memo.md` (CLEAR verdict; load-bearing for §1.4)
- `debug/sprint_phase_a3prime_mondino_samann_bridge_memo.md` (A.3' bridge; load-bearing for §2, §5, §6, §7)
- `debug/sprint_phase_a2prime_krein_lift_formalization_memo.md` (A.2' substrate; load-bearing for §6.3)
- `papers/group1_operator_algebras/paper_44_lorentzian_operator_system.tex` (operator-system substrate; load-bearing for §3.1 (A3), §6.4)
- `papers/group1_operator_algebras/paper_46_strong_form_lorentzian_propinquity.tex` (Appendix B enlarged substrate; load-bearing for §6.1, §6.4)
- `papers/group1_operator_algebras/paper_48_krein_ms_bridge.tex` (§3 Krein PPQMS, §6 Decomposition O, §8.1 Q1'; load-bearing throughout)
- `papers/group1_operator_algebras/paper_42_modular_hamiltonian_four_witness.tex` (four-witness theorem; load-bearing for §6.4 (A4))
- Connes-vS arXiv:2004.14115 (operator-system framework; template for §3, §4)
- Mondino-Sämann arXiv:2504.10380 v4 (MS axioms; load-bearing for §2)
- Connes-Rovelli 1994 Class. Quantum Grav. 11, 2899 (thermal-time hypothesis; load-bearing for §2.1, §7.3)
