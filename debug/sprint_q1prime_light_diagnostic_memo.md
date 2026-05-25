# Sprint Q1'-Light — Diagnostic of strong-form Krein–MS bridge feasibility on Paper 46 Appendix B enlarged substrate

**Date:** 2026-05-24

**Sprint type:** DIAGNOSTIC ONLY. No production code, no paper modifications, no Q1' theorem attempted.

**Source materials read:**
- `papers/group1_operator_algebras/paper_46_strong_form_lorentzian_propinquity.tex` — especially Appendix B (Definition 5.2 enlarged substrate, eq.~(1.41) chirality-asymmetric flip form, Theorem 5.4 strict-strong-form bound, Reading 1/2/3 interpretation §5.5)
- `papers/group1_operator_algebras/paper_48_krein_ms_bridge.tex` — §3 Krein-pointed proper QMS substrate, §2 main bridge theorem, §6 (B2) Decomposition O lemma, §8.1 Q1' three-step decomposition
- `papers/group1_operator_algebras/paper_42_modular_hamiltonian_four_witness.tex` — four-witness theorem on the BW wedge (Thm 5.4 BW-α + Thm 6.3 BW-γ + Thm 7.1 flow conjugacy + Cor 8.1 six-witness collapse)
- `papers/group1_operator_algebras/paper_44_lorentzian_operator_system.tex` — Avery–Wen–Avery basis for the spatial 3-Y multipliers + Krein-positive substrate
- `debug/sprint_phase_a4prime_a_gb2_super_additivity_memo.md` — A.4'-A G-B2 closure at K⁺-weak-form (Decomposition O Case (iii) empty at K⁺-weak-form)
- `debug/sprint_phase_a3prime_mondino_samann_bridge_memo.md` — R3 Connes–Rovelli thermal time as winning candidate; bridge functor W definition; F2 mismatch resolution

**Aggregate verdict (1-sentence):** **PARTIAL-WITH-NAMED-OBSTRUCTIONS — Q1'.A topography enlargement compresses ONLY in the M=0 sub-case (where Abelianness is preserved structurally), but in that sub-case the enlargement does NOT activate Case (iii) of Decomposition O and so does NOT exercise the substantive strong-form content; activating the substantive content requires the full M≠0 enlargement, on which chirality-flipping generators do not commute and Q1'.B fragments into Option α/β/γ workarounds vs Option δ non-commutative MS extension; recommend Q1'-staged sprint plan with Q1'-Phase-1 (M=0 Abelian sub-case, sprint-scale 1–3 weeks) as a stepping stone that establishes machinery but does NOT close G-B2 at strict-strong-form, followed by Q1'-Phase-2 (Option γ operator-system replacement, 1–2 months) OR Q1'-Phase-3 (Option δ non-commutative MS, 6–12 months) as the substantive multi-month follow-on.**

---

## §1. Foundation summary

### 1.1. The substrate distinction

The four-tuple architecture across the relevant papers:

| Substrate | Paper | Defining condition | Status |
|:----------|:------|:-------------------|:-------|
| Natural ($\mathcal{O}^L$) | Paper 44 Prop 5.1 | $[J_L, a] = 0$ for every multiplier | Foundation |
| K⁺-weak-form bridge topography ($\mathcal{M}^L$) | Paper 48 Lemma 2.15 | Natural ∩ $\{M = 0\}$ | Bridge target |
| Enlarged ($\mathcal{O}^L_{\mathrm{enlarged}}$) | Paper 46 Def 5.2 | Natural ∪ chirality-flip generators ($\{J_L, M^{\mathrm{flip}}\} = 0$) | Strong-form bound |
| Enlarged topography ($\mathcal{M}^L_{\mathrm{enlarged}}$) | Q1'.A target | Enlarged ∩ "Abelian topography axiom" | Open |

The natural topography $\mathcal{M}^L$ is the *intersection* of the K⁺-preservation constraint (commutes with $J_L$) AND the M=0 restriction (preserves m_j = 0 sector). The K⁺-weak-form bridge of Paper 48 uses this intersection. The strong-form bound of Paper 46 uses the union of $\mathcal{O}^L$ with chirality-flip generators (no M restriction). The Q1'-target enlarged topography is *not yet defined* — that's the substantive content of Q1'.A.

### 1.2. The chirality-asymmetric DIAGONAL form

Paper 46 Appendix B Remark 5.3 carefully distinguishes two candidate enlargement forms:

**(Off-block-diagonal form, NOT used):** $M^{\mathrm{flip}} = \begin{pmatrix} 0 & W \\ W^* & 0 \end{pmatrix}_\chi$ in the chiral basis. This was the "naïve" form contemplated in earlier Paper 46 §5.2. For Hermitian $W$, $\{\gamma^0, M^{\mathrm{flip}}\} = \begin{pmatrix} W^* + W & 0 \\ 0 & W + W^* \end{pmatrix} \neq 0$ in general — *fails* the $\{J_L, M^{\mathrm{flip}}\} = 0$ condition.

**(Chirality-asymmetric diagonal form, USED in Paper 46 Appendix B Definition 5.2):** $M^{\mathrm{flip}}_{N,L,M} = \mathrm{diag}(W^{NLM}, -W^{NLM})$ in the chiral basis. The anti-commutator computation is bit-exact zero:
$\{\gamma^0, \mathrm{diag}(W, -W)\} = \begin{pmatrix} 0 & I \\ I & 0 \end{pmatrix} \begin{pmatrix} W & 0 \\ 0 & -W \end{pmatrix} + \begin{pmatrix} W & 0 \\ 0 & -W \end{pmatrix}\begin{pmatrix} 0 & I \\ I & 0 \end{pmatrix} = \begin{pmatrix} 0 & -W + W \\ W - W & 0 \end{pmatrix} = 0.$

This is the load-bearing form. $W^{NLM}$ is an **Avery–Wen–Avery Weyl-spinor multiplier** of Paper 44 §3 — i.e., a spatial multiplication operator by the 3-Y function $Y^{(3)}_{NLM}$ acting on the appropriate spinor bundle component.

### 1.3. The substantive question for Q1'-Light

Each of Q1'.A, Q1'.B, Q1'.C separates into a sub-case structure that we test below.

---

## §2. Q1'.A topography enlargement test — the substantive structural test

The substantive test is: **do the chirality-flip generators $\{M^{\mathrm{flip}}_{N,L,M}\}$ commute with each other AND with the existing M-diagonal generators?**

### 2.1. Algebraic structure of $[M^{\mathrm{flip}}_{1}, M^{\mathrm{flip}}_{2}]$

For two chirality-flip generators $M^{\mathrm{flip}}_{N_1, L_1, M_1} = \mathrm{diag}(W^1, -W^1)$ and $M^{\mathrm{flip}}_{N_2, L_2, M_2} = \mathrm{diag}(W^2, -W^2)$:
\[
[M^{\mathrm{flip}}_1, M^{\mathrm{flip}}_2]
= \mathrm{diag}(W^1, -W^1) \mathrm{diag}(W^2, -W^2) - \mathrm{diag}(W^2, -W^2) \mathrm{diag}(W^1, -W^1)
\]
The product of two block-diagonal $2 \times 2$ matrices is block-diagonal, with the upper block being $W^1 \cdot W^2$ and the lower block being $(-W^1)(-W^2) = W^1 W^2$:
\[
\mathrm{diag}(W^1, -W^1) \mathrm{diag}(W^2, -W^2) = \mathrm{diag}(W^1 W^2, W^1 W^2).
\]
Symmetrically:
\[
\mathrm{diag}(W^2, -W^2) \mathrm{diag}(W^1, -W^1) = \mathrm{diag}(W^2 W^1, W^2 W^1).
\]
So:
\[
[M^{\mathrm{flip}}_1, M^{\mathrm{flip}}_2] = \mathrm{diag}([W^1, W^2], [W^1, W^2]).
\]

**Key observation:** the chirality-flip-commutator reduces to the spatial-multiplier commutator $[W^{N_1 L_1 M_1}, W^{N_2 L_2 M_2}]$ on each chirality sector.

The chirality-asymmetric doubling **preserves the spatial commutation structure** — it does not introduce new non-commutativity beyond what already exists in the underlying spatial multipliers $\{W^{NLM}\}$.

### 2.2. Case decomposition on $M$ quantum number

**Case A — $M = 0$ restriction (chirality-flip topography $\mathcal{M}^{L, \mathrm{flip}}_{M=0}$):**

$W^{NL0}$ for $M = 0$ is a multiplication operator by the M=0 spatial 3-Y function $Y^{(3)}_{NL0}$ on the spinor bundle. Multiplication operators on the same Hilbert space commute (any two functions $f, g$ on a space satisfy $f \cdot g = g \cdot f$ as multiplication operators).

Hence $[W^{N_1 L_1 0}, W^{N_2 L_2 0}] = 0$ for all $(N_1, L_1), (N_2, L_2)$ — and so:
\[
[M^{\mathrm{flip}}_{N_1 L_1 0}, M^{\mathrm{flip}}_{N_2 L_2 0}] = \mathrm{diag}(0, 0) = 0 \quad \text{bit-exact.}
\]

**Furthermore**, the chirality-flip generators with $M = 0$ commute with the existing M-diagonal topography generators $M^{\mathrm{spat}}_{N', L', 0} \otimes I_{N_t} \in \mathcal{M}^L$, because the existing topography generators are also spatial multiplication operators by $Y^{(3)}_{N' L' 0}$, and any two multiplication operators commute.

**ABELIANNESS PRESERVED in the M=0 sub-case** (Case A).

**Case B — Full $M \neq 0$ enlarged substrate (chirality-flip topography $\mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}}$):**

$W^{NLM}$ for $M \neq 0$ is NOT a multiplication operator on the spinor bundle — it carries non-trivial $m_j$ structure and acts on the spatial part as a Wigner-rotation transformation on the m_j-quantum-number labels. For two such operators with different $(M_1, M_2)$:
\[
[W^{N_1 L_1 M_1}, W^{N_2 L_2 M_2}] \neq 0 \quad \text{in general.}
\]

The non-commutativity is the same Clebsch–Gordan recoupling structure that makes the operator system $\mathcal{O}^L$ non-Abelian in the first place (Paper 44 §3 propagation number > 1).

Hence $[M^{\mathrm{flip}}_{N_1 L_1 M_1}, M^{\mathrm{flip}}_{N_2 L_2 M_2}] \neq 0$ in general.

**ABELIANNESS BREAKS in the full enlarged sub-case** (Case B).

### 2.3. Q1'.A verdict per sub-case

| Sub-case | Abelianness | Topography axioms (Def 1.40-K) | Q1'.A status |
|:---------|:------------|:-------------------------------|:-------------|
| Case A (M=0 chirality-flip enlargement) | PRESERVED | All four hold (Abelian, contains strictly positive $K_\alpha^W$, character $\hat{\omega}_W^L$, approximate unit) | COMPRESSES — bit-exact verification |
| Case B (full M≠0 chirality-flip enlargement) | BREAKS | (a) Abelian FAILS, (c) character may not exist | GENUINELY HARDER — requires structural workaround |

The Paper 48 §8.1 Q1'.A statement ("the enlargement may break Abelianness — chirality-flipping generators do not commute with each other in general") is **confirmed at the operator-algebraic level** as a Case-B statement.

### 2.4. The substantive obstruction at Q1'.A

A Case-A-only enlargement (M=0 chirality-flip) closes Q1'.A trivially, but at the cost of **not exercising the substantive content of the strong-form bridge**. The reason:

The Paper 48 §6 Decomposition O Case (iii) — three pure characters $(\hat{\omega}_x, \hat{\omega}_y, \hat{\omega}_z)$ such that all three pairs lie on different modular orbits — is empty at the K⁺-weak-form (M=0) topography because modular flow preserves the $(N, L)$ labels (Lemma 6.1.4 of Paper 48), and two characters on different orbits have different $(N, L)$, which forces $\ell^L = -\infty$ via Cases (ii) or (iv).

Adding chirality-flip generators with $M = 0$ does NOT introduce new $(N, L)$ labels — it only adds a chirality-asymmetric sign factor on the existing labels. Hence Decomposition O Case (iii) remains empty under the Case-A enlargement. The strong-form bridge structurally requires Case-B enlargement (full $M \neq 0$ chirality-flip content) to activate Case (iii), at which point Abelianness breaks.

**Q1'.A is bimodal**: the Abelian sub-case (Case A) compresses but doesn't exercise the substantive content; the substantive sub-case (Case B) breaks Abelianness and forces Q1'.B into substantial workaround territory.

---

## §3. Q1'.B Gelfand spectrum or non-commutative extension

### 3.1. Case A — Gelfand spectrum well-defined

In Case A (M=0 chirality-flip enlargement), $\mathcal{M}^{L, \mathrm{flip}}_{M=0}$ is an Abelian C*-algebra, and its Gelfand spectrum is well-defined. The characters of $\mathcal{M}^{L, \mathrm{flip}}_{M=0}$ are determined by the eigenvalues of the spatial multiplication operators $W^{NL0}$, which are labeled by points of the underlying spatial S³ wedge plus a binary chirality sign (the "upper" vs "lower" chiral sector that the $\pm W$ acts on).

This gives a Gelfand spectrum $\widehat{\mathcal{M}^{L, \mathrm{flip}}_{M=0}} = \widehat{\mathcal{M}^L} \sqcup \widehat{\mathcal{M}^L}$ (two-fold cover labeled by chirality sign), with the modular flow acting trivially on the chirality label (because the BW geometric Hamiltonian $K_\alpha = J_{\mathrm{polar}}$ has integer spectrum on the spinor bundle and commutes with the chirality grading — Paper 42 Thm 5.4).

**Q1'.B at Case A compresses** to a mechanical verification that the two-fold-cover spectrum inherits the modular-causal structure of the underlying $\widehat{\mathcal{M}^L}$, with each chirality sheet getting its own copy of the modular flow. This is sprint-scale work, on the order of 1–2 weeks.

### 3.2. Case B — Abelianness breaks, structural workaround needed

In Case B (full enlarged substrate), $\mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}}$ is non-Abelian. Gelfand–Naimark does not apply directly. Four workaround options:

**Option α — Joint spectrum of commuting subfamilies + covering structure.** The full enlarged substrate decomposes into a family of commuting subalgebras (one per M-eigenvalue), and one can take the joint spectrum of each subfamily and assemble a covering structure. **Feasibility:** medium. The covering would have non-trivial overlap structure tied to the Clebsch–Gordan recoupling of different M-sectors. Estimated effort: 4–8 weeks of focused work to assemble the cover and verify its compatibility with the MS Def 3.6 axioms. The structural concept is essentially within sprint scale, but the bookkeeping is delicate.

**Option β — Twisted Gelfand spectrum.** Use a twisted-Gelfand construction for self-adjoint operator subsystems (Effros–Ruan, operator-system theory). The twisted spectrum exists but is not a topological space — it's a quantum measurable structure. Feasibility: medium-low. Connecting it to the MS pre-length-space framework (which requires a topological space) requires additional structural work. Estimated effort: 2–3 months.

**Option γ — Operator-system replacement.** Skip Gelfand entirely; work directly at the operator-system level (replace the topographic structure $\mathcal{M}^L$ with an operator-system structure that supports a "modular metric"). This is the most natural mathematical setting: Connes–vS 2021 operator-system framework is the existing precedent (already used in Paper 44 for the K⁺-weak-form substrate). Feasibility: medium-high. The MS pre-length space concept would need to be replaced by an operator-system analog ("operator-system Lorentzian pre-length space"). This is a structural shift — the natural mathematical category becomes a different category from MS, but the structural content (timelike, causal, ε-net) can be transcribed. Estimated effort: 6–10 weeks for the construction and verification. **This is the most natural Q1'-staged path.**

**Option δ — Genuine non-commutative MS pre-length spaces.** Build the non-commutative MS pre-length-space concept from scratch. Mondino–Sämann pre-length spaces are commutative (an underlying set of "points" with timelike-precedence and time-separation). A non-commutative version would generalize "points" to "states on an operator system" and "timelike-precedence" to a modular-flow-induced partial order on states. **This is the 6–12 month NCG-research target named in Paper 48 §8.1 and §8.4.** It is genuinely outside the published literature. Estimated effort: 6–12 months.

### 3.3. Q1'.B verdict per option

| Option | Mathematical setting | Effort estimate | Sprint vs multi-month |
|:-------|:---------------------|:----------------|:----------------------|
| α (joint spectrum + cover) | Existing Gelfand + cover bookkeeping | 4–8 weeks | Sprint-scale boundary |
| β (twisted Gelfand) | Operator-system Effros–Ruan | 2–3 months | Multi-month |
| γ (operator-system replacement) | Connes–vS 2021 + MS-analog | 6–10 weeks | Sprint-scale boundary |
| δ (genuine non-commutative MS) | New NCG concept | 6–12 months | Multi-month |

Options α and γ are at the sprint–multi-month boundary. Options β and δ are clearly multi-month. The Paper 48 §8.4 Q2' open question explicitly names Option δ as a substantial NCG-research target.

**Q1'.B as a standalone task is NOT cleanly sprint-scale.** Even Options α/γ require 1.5–2 months of focused work to land the construction.

---

## §4. Q1'.C off-orbit super-additivity feasibility

### 4.1. The Paper 42 four-witness theorem as load-bearing input

The four-witness theorem of Paper 42 establishes the following at finite cutoff on the BW wedge (Hartle–Hawking + Sewell + Bisognano–Wichmann + Unruh, all collapsed to one operator-system construction per Cor 8.1):

1. **(Thm 5.4 BW-α)** $K_\alpha^W = J_{\mathrm{polar}}^W$ has integer spectrum two_m_j on the spinor bundle, restricted to the wedge.
2. **(Thm 6.3 BW-γ)** $K_{TT}^W = -\log \Delta$ obtained from polar decomposition of $S = J_{TT}\Delta^{1/2}$ on the GNS Hilbert–Schmidt space gives an equivalent modular Hamiltonian.
3. **(Thm 7.1)** Flow conjugacy: $\sigma_t^{TT}(a) = \sigma_{-t}^\alpha(a)$ at bit-exact finite cutoff.
4. **(Cor 8.1)** Six-witness collapse: the four physical witnesses reduce to the same operator-system construction.

The **key structural property** these theorems establish: $\sigma_{2\pi}^{TT}(a) = a$ bit-exact for every multiplier $a \in \mathcal{O}^L$, i.e., **modular flow has period $2\pi$**.

### 4.2. The on-orbit case (already closed)

Within a single modular orbit, additivity of modular flow gives:
\[
\sigma_{t_1 + t_2}^\alpha(a) = \sigma_{t_2}^\alpha(\sigma_{t_1}^\alpha(a)).
\]
The Paper 48 §6.2 proof of (B2) Case (i) uses this directly:
\[
\ell^L(\hat{\omega}_x, \hat{\omega}_y) + \ell^L(\hat{\omega}_y, \hat{\omega}_z) = \kappa_g(t_1 + t_2) = \ell^L(\hat{\omega}_x, \hat{\omega}_z) \quad \text{(equality at on-orbit case).}
\]

This is the **load-bearing transport of Paper 42 four-witness theorem to the bridge** at the K⁺-weak-form level. At the strong-form level (Q1'.C), the on-orbit case continues to apply verbatim — the modular flow on a single orbit is still $2\pi$-periodic and additive, whether the topography is M=0-restricted (Case A) or full M≠0 (Case B). On-orbit at strong-form is **automatic from Paper 42**.

### 4.3. The off-orbit case at strong-form — what's actually different

At the K⁺-weak-form, Case (iii) of Decomposition O is empty (three different orbits is structurally impossible under the M=0 topography). At strong-form (full M≠0), Case (iii) becomes non-empty: there are triples $(\hat{\omega}_x, \hat{\omega}_y, \hat{\omega}_z)$ such that $(x, y)$, $(y, z)$, $(x, z)$ lie on three different orbits, with the orbits not all sharing a single $(N, L)$ label.

For such a triple, the substantive question is: does the proper-time accumulation $\ell^L(x, y) + \ell^L(y, z)$ satisfy a *strict* super-additivity $\le \ell^L(x, z)$?

The physical content here is the **special-relativistic twin paradox**: a worldline that detours through $y$ before reaching $z$ accumulates *less* proper time than a direct worldline from $x$ to $z$. This is well-known in classical Lorentzian geometry. The Mondino–Sämann reverse triangle inequality $\ell(x, y) + \ell(y, z) \le \ell(x, z)$ is exactly this property at the level of synthetic Lorentzian pre-length spaces.

### 4.4. Three candidate approaches to transport Paper 42 to off-axis

**(a) Direct four-witness transport via covering structure.** Lift the four-witness theorem from a single wedge to a cover of the union of boost orbits via a structurally compatible covering. Feasibility: depends on the covering construction in Q1'.B. If Q1'.B is closed via Option α or γ, this approach becomes mechanically derivable from the per-cover-cell modular flow. Estimated effort: 2–4 weeks ON TOP of Q1'.B's resolution.

**(b) Modular-flow composition across orbits.** Use the Tomita–Takesaki modular flow's composition properties to derive super-additivity across different orbits. The Tomita–Takesaki structure (Paper 42 §6) is single-orbit by construction — the polar decomposition $S = J_{TT}\Delta^{1/2}$ on the GNS Hilbert–Schmidt space gives ONE modular flow per pin state. Different orbits would require different KMS states, hence different modular structures. Composition across distinct modular structures is not a standard Tomita–Takesaki construction. Feasibility: medium-low. Estimated effort: 2–4 months.

**(c) Connes–Rovelli thermal-time stack.** The thermal-time hypothesis (R3 winning candidate at the K⁺-weak-form bridge, per A.3' memo §5) identifies modular flow on a single KMS state with geometric time on the corresponding wedge. **Across different wedges (different orbits), thermal-time stacks as a fibered structure** over the orbit-parameter space. Geometric super-additivity (twin paradox) is the dual of additivity-along-each-fiber + concavity-across-fibers. Feasibility: medium-high. This is the most natural physical approach but requires constructing the fibration explicitly. Estimated effort: 1.5–3 months.

### 4.5. Q1'.C verdict

The on-orbit content of Q1'.C is automatic from Paper 42 (transports verbatim from K⁺-weak-form to strong-form). The off-axis content is the substantive new ingredient, with three candidate approaches:

| Approach | Effort | Sprint-scale boundary |
|:---------|:-------|:----------------------|
| (a) Direct four-witness transport via cover | 2–4 weeks (gated on Q1'.B) | Sprint-scale IF Q1'.B closes via α/γ |
| (b) Modular-flow composition | 2–4 months | Multi-month |
| (c) Connes–Rovelli thermal-time stack | 1.5–3 months | Multi-month |

**Q1'.C as a standalone task is NOT cleanly sprint-scale at strong-form (off-orbit).** The on-orbit content is automatic, but the off-orbit content requires substantive new structural work that depends on Q1'.B's resolution.

---

## §5. Aggregate verdict and recommended Q1' sprint structure

### 5.1. Aggregate verdict

| Sub-task | Compression assessment | Bottleneck |
|:---------|:-----------------------|:-----------|
| Q1'.A | Bimodal (Case A compresses, Case B doesn't) | Substantive content requires Case B → Abelianness breaks |
| Q1'.B | NOT cleanly sprint-scale (best options 1.5–2 months) | Non-Abelianness forces structural workaround |
| Q1'.C | On-orbit automatic; off-orbit requires Q1'.B + additional work | Cascades from Q1'.B |

**Verdict: PARTIAL-WITH-NAMED-OBSTRUCTIONS.** The Q1' three-step decomposition has:
- Q1'.A "POSITIVE in Case A only, but Case A doesn't exercise substantive content" — a structural bifurcation, not pure compression
- Q1'.B "PARTIAL with four options, none cleanly sprint-scale" — best path (Option γ operator-system replacement) is 6–10 weeks, but this changes the mathematical setting away from MS
- Q1'.C "On-orbit automatic, off-orbit gated on Q1'.B" — cannot be tackled independently

The original Paper 48 §8.1 effort estimate of 3–6 months at sprint cadence (expanding to 6–12 months if Q1'.B forces Option δ) is **confirmed by this diagnostic** as approximately correct. The diagnostic refines:

- The 3-month minimum is realistic IF Q1'.B closes via Option α or γ (sprint-scale workarounds)
- The 6-month + scenario applies IF Q1'.B forces Option β or δ
- The 6–12 month NCG-research scenario applies IF Option δ is the only resolution

### 5.2. Recommended Q1' sprint structure

**Recommendation: Q1'-STAGED with three sub-sprints, NOT a single 1–3 week burst.**

**Q1'-Phase-1 (sprint-scale, 1–3 weeks): Case A enlargement as stepping stone.** Close the M=0 chirality-flip topography enlargement at theorem-grade rigor. This establishes:
- The two-fold-cover Gelfand spectrum $\widehat{\mathcal{M}^{L, \mathrm{flip}}_{M=0}} = \widehat{\mathcal{M}^L} \sqcup \widehat{\mathcal{M}^L}$
- The chirality-graded bridge functor $W^{\mathrm{flip}}_{M=0}: \mathbf{KreinMetaMet}^{\mathrm{flip}}_{M=0, \mathrm{pp}} \to \mathbf{LorPLG}^{2:1}_{\mathrm{cov}}$
- Bit-exact verification at finite cutoff $(n_{\max}, N_t) \in \{(2, 3), (3, 5)\}$

**Deliverable:** sprint-scale memo + Paper 48 §8.1 update naming Phase-1 as Q1'-Phase-1-Done. **Does NOT close G-B2 at strict-strong-form** — that requires Phase-2 (full M ≠ 0).

**Q1'-Phase-2 (6–10 weeks, IF the agent decides to pursue the non-Abelian extension): Option γ operator-system replacement.** Build an "operator-system Lorentzian pre-length space" (OSLPLS) concept that generalizes Mondino–Sämann to non-commutative topographies. Verify that the strong-form bridge extends to $W^{\mathrm{flip}}_{\mathrm{full}}: \mathbf{KreinMetaMet}^{\mathrm{flip}}_{\mathrm{full}, \mathrm{pp}} \to \mathbf{OSLPLG}_{\mathrm{cov}}$. Close G-B2 at strict-strong-form via the Connes–Rovelli thermal-time stack (Q1'.C approach (c)) on the new substrate.

**Deliverable:** Paper 49 draft (math.OA), 11th in the GeoVac series. Closes Q1' at the operator-system-replacement level. **Does NOT close Q2'** (non-commutative MS pre-length space concept) — that remains open.

**Q1'-Phase-3 (6–12 months, multi-month NCG-research): Option δ non-commutative MS pre-length spaces.** Build the non-commutative MS pre-length-space concept from scratch, generalizing both Mondino–Sämann (commutative side) and Q1'-Phase-2 (operator-system side). This is the substantive NCG-research target the literature has not yet attempted.

**Deliverable:** Multi-paper NCG-research arc, comparable in scope to Latrémolière propinquity development. Closes Q1' + Q2' together.

### 5.3. Per-sub-task summary table

| Sub-task | Sub-case | Verdict | Effort |
|:---------|:---------|:--------|:-------|
| Q1'.A | Case A (M=0 flip) | COMPRESSES | 1–3 weeks |
| Q1'.A | Case B (full enlarged) | BREAKS ABELIANNESS | Cascades to Q1'.B |
| Q1'.B | Option α (joint spectrum + cover) | Sprint-multi-month boundary | 4–8 weeks |
| Q1'.B | Option β (twisted Gelfand) | Multi-month | 2–3 months |
| Q1'.B | Option γ (operator-system replacement) | Best sprint-scale path | 6–10 weeks |
| Q1'.B | Option δ (non-commutative MS) | Multi-month NCG-research | 6–12 months |
| Q1'.C | On-orbit (single orbit) | AUTOMATIC from Paper 42 | 0 weeks |
| Q1'.C | Off-orbit approach (a) | Gated on Q1'.B α/γ | 2–4 weeks ON TOP |
| Q1'.C | Off-orbit approach (b) | Multi-month | 2–4 months |
| Q1'.C | Off-orbit approach (c) | Sprint-multi-month boundary | 1.5–3 months |

### 5.4. R3 Connes–Rovelli winning candidate at strong-form?

The R3 Connes–Rovelli thermal-time hypothesis was the K⁺-weak-form winner (A.3' memo §5, used the F2 mismatch resolution to construct the bridge functor $W$). At strong-form, the question is whether R3 *still* works, or whether the F2 mismatch resolution needs a different candidate.

**R3 status at strong-form:** R3 continues to work for the on-orbit content (Case (i) of Decomposition O, automatic from Paper 42). For the off-orbit content (Case (iii), the substantive new ingredient at strong-form), R3 gives the structural framework but does not by itself produce the strict super-additivity inequality. The off-orbit super-additivity requires the **Connes–Rovelli thermal-time stack across distinct KMS states** (Q1'.C approach (c)), which is a fibration generalization of R3 that does not exist in the literature.

So **R3-extended (with the thermal-time stack)** is the natural candidate at strong-form, but it is NOT a drop-in replacement for K⁺-weak-form R3 — it is a fibered generalization requiring substantive new construction. The Connes–Rovelli 1994 paper itself does not address the fibered case.

---

## §6. Honest scope statement

This is a **diagnostic memo**. The verdict, sub-task effort estimates, and recommended sprint structure are based on:
1. The structural content of Paper 46 Appendix B Definition 5.2 (chirality-asymmetric diagonal form $M^{\mathrm{flip}} = \mathrm{diag}(W, -W)$)
2. The structural content of Paper 48 §3 (Krein-pointed proper QMS substrate), §6 (Decomposition O lemma), §8.1 (Q1' three-step decomposition)
3. The A.4'-A memo §6.5 finding that off-orbit super-additivity at K⁺-weak-form is structurally trivial via Decomposition O Case (iii) emptiness
4. The Paper 42 four-witness theorem as load-bearing input for the on-orbit modular-flow additivity

The diagnostic surfaces **one substantive new finding beyond the Paper 48 §8.1 named decomposition**: the Q1'.A topography enlargement is **structurally bimodal** between Case A (M=0 chirality-flip) and Case B (full M≠0 chirality-flip). Case A compresses cleanly but does NOT activate Decomposition O Case (iii) (the structurally novel off-orbit configurations). Case B activates Case (iii) but breaks Abelianness. This bifurcation is implicit in Paper 48 §8.1 (the "topography enlargement may break Abelianness" remark) but not made structurally explicit — the Q1'-Light diagnostic makes it explicit.

The diagnostic does NOT:
- Write any Q1' theorems or formalization
- Attempt to close Q1'.A Case A even at sprint scale (that's Q1'-Phase-1 of the recommended structure)
- Build the non-commutative MS extension or any of the four Q1'.B options
- Modify any production code or papers

**The diagnostic-before-engineering rule** (CLAUDE.md memory): the diagnostic surfaces the Q1'.A Case A / Case B bifurcation as a structural property that informs the Q1' sprint plan. Without the diagnostic, a "single 1–3 week sprint" attempt at Q1' would likely have hit the Case-A vs Case-B bifurcation mid-sprint and either stalled (if attempted on Case B blindly) or produced a misleading partial result (if attempted on Case A without recognizing it doesn't exercise the substantive content).

The "POSITIVE-COMPRESSES" verdict suggested in the decision gate **is NOT supported** by the diagnostic. The "NEGATIVE-multi-month" verdict is **partially supported** — Q1' is genuinely multi-month at the substantive content level, with best-case 6–10 weeks for the operator-system-replacement path (Phase 2) and 6–12 months for the non-commutative MS path (Phase 3). The **PARTIAL-WITH-NAMED-OBSTRUCTIONS** verdict is the honest reading: Q1' is a staged multi-month follow-on, but the first stage (Q1'-Phase-1, Case A) compresses to sprint scale and establishes a substantial portion of the machinery.

---

## §7. Cross-references

- `papers/group1_operator_algebras/paper_46_strong_form_lorentzian_propinquity.tex` — Appendix B definition of enlarged substrate and chirality-asymmetric flip form
- `papers/group1_operator_algebras/paper_48_krein_ms_bridge.tex` — §3 Krein PPQMS substrate, §6 Decomposition O, §8.1 Q1' three-step decomposition (the source of the three-task structure tested in this diagnostic)
- `papers/group1_operator_algebras/paper_42_modular_hamiltonian_four_witness.tex` — four-witness theorem on the BW wedge (load-bearing for on-orbit content)
- `papers/group1_operator_algebras/paper_44_lorentzian_operator_system.tex` — Avery–Wen–Avery Weyl-spinor multipliers $W^{NLM}$ (load-bearing for the spatial-multiplier structure of the chirality-flip generators)
- `debug/sprint_phase_a4prime_a_gb2_super_additivity_memo.md` — A.4'-A G-B2 closure at K⁺-weak-form via Decomposition O Case (iii) emptiness; §8 Q1' refinement
- `debug/sprint_phase_a3prime_mondino_samann_bridge_memo.md` — R3 Connes–Rovelli thermal-time hypothesis as winning candidate; bridge functor $W$ definition; §6.6 Q1' open question

---

## §8. Output for the report-back

**Headline verdict (1 sentence):** PARTIAL-WITH-NAMED-OBSTRUCTIONS — Q1'.A topography enlargement is structurally bimodal (M=0 Case A compresses but doesn't exercise substantive content; full M≠0 Case B activates the substantive off-orbit content but breaks Abelianness), forcing Q1'.B into four-option structural-workaround territory of which the best sprint-scale path (Option γ operator-system replacement) is 6–10 weeks and the genuine non-commutative MS extension (Option δ) remains a 6–12 month NCG-research target.

**Q1'.A Abelianness result:** BIMODAL. Chirality-flip generators commute trivially in the M=0 restriction (multiplication operators always commute) → ABELIANNESS PRESERVED in Case A. Chirality-flip generators DO NOT commute in the full enlarged substrate (M≠0 carries non-trivial m_j structure) → ABELIANNESS BREAKS in Case B. The substantive Decomposition O Case (iii) content requires Case B, which forces non-Abelian Q1'.B workarounds.

**Most surprising finding:** The Q1'.A Case A vs Case B bifurcation is implicit in Paper 48 §8.1 but not made explicit. Case A would be tempting to attempt as a "Q1' first-pass closure" since it preserves Abelianness — but it does NOT activate the substantive strong-form content (Decomposition O Case (iii) remains empty under M=0 chirality-flip enlargement, just as it is for the K⁺-weak-form bridge). A sprint that closes Case A without recognizing this would land a misleading-positive result. The diagnostic-before-engineering rule (CLAUDE.md memory) caught this in 1.5 hours of structural analysis.

**Recommended Q1' sprint structure:** STAGED, NOT single sprint. Q1'-Phase-1 (sprint-scale, 1–3 weeks): close Case A enlargement at theorem-grade rigor as stepping stone establishing the chirality-graded bridge functor at the two-fold-cover Gelfand spectrum level. Q1'-Phase-2 (6–10 weeks): close Q1' at strong-form via Option γ operator-system replacement, producing Paper 49 (math.OA) as the 11th GeoVac series paper. Q1'-Phase-3 (6–12 months): only if pursued — build the non-commutative MS pre-length-space concept from scratch (NCG-research arc, comparable in scope to Latrémolière propinquity development).

**End of memo.**
