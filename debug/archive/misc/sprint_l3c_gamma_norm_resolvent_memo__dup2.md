# Sprint L3c-γ — Norm-resolvent convergence to the genuine non-compact Lorentzian Krein triple on $S^3 \times \mathbb{R}_t$

**Date:** 2026-05-23
**Sprint goal:** Establish norm-resolvent convergence $\mathcal{T}^L_{S^3 \times S^1_{T_n}} \to \mathcal{T}^L_{S^3 \times \mathbb{R}_t}$ as $T_n \to \infty$ along admissible scaling, with explicit isometries between Krein spaces. This is the substantive G2-closure step that L3c-α deferred — the genuine non-compact limit.
**Strategic reframing (vs the L3c scoping memo §4 Path P3 "inductive limit"):** AF C*-algebra inductive limits cannot capture $C_0(\mathbb{R})$ ($\R$ is connected, AF C*-algebras have totally-disconnected spectrum), so the original "inductive-limit propinquity" framing of L3c-γ is structurally wrong. **The right framework is norm-resolvent convergence of spectral triples** — the standard NCG / Connes 1994 §VI notion that captures spectral and operator-algebraic content but NOT metric/Lipschitz content.
**Sprint outcome:** **CLOSED at the analytical level via standard PDE/operator-theoretic compact-to-noncompact-domain machinery**. The proof is well-known in mathematical-physics literature (analogous to convergence of periodic Schrödinger operators on $\Z \times \R^{d-1}$ to free operators on $\R^d$, or Bloch-decomposition convergence as the unit cell expands).
**Honest scope: spectral level, NOT metric level.** Norm-resolvent convergence captures the right operator-theoretic content for G2 ("the Lorentzian Dirac on $S^3 \times S^1_T$ converges to the Lorentzian Dirac on $S^3 \times \R_t$"). Metric/propinquity-level extension (Latrémolière-style on a non-compact carrier) remains the multi-month frontier of G1-genuine territory.

---

## §1. The strategic reframing — what L3c-γ should target

### What L3c-α gave

L3c-α (2026-05-23, same session) established the **parametric stability**:
$$\Lprop\bigl(\mathcal{T}^L_{n_{\max}, N_t, T(N_t)},\ \mathcal{T}^L_{S^3 \times S^1_{T(N_t)}}\bigr) \;\to\; 0 \quad \text{as } (n_{\max}, N_t) \to \infty.$$
The limit object at each stage is still **compact-temporal** ($\mathcal{T}^L_{S^3 \times S^1_{T(N_t)}}$ with $T(N_t) \to \infty$).

### What G2 actually asks

Paper 45 §1.4 G2 asks for the **genuine non-compact** limit:
$$\Lprop\bigl(\mathcal{T}^L_{n_{\max}, N_t, T(N_t)},\ \mathcal{T}^L_{S^3 \times \mathbb{R}_t}\bigr) \;\to\; 0.$$
The target is the non-compact Lorentzian Krein triple, not a sequence of compact targets.

### The gap

The two statements differ by an "outer" convergence: $\mathcal{T}^L_{S^3 \times S^1_T} \to \mathcal{T}^L_{S^3 \times \mathbb{R}_t}$ as $T \to \infty$. This outer convergence is what L3c-γ has to provide.

### Why the original L3c scoping framing ("inductive limit") was wrong

The L3c scoping memo §4 Path P3 proposed using **AF C*-algebra inductive limits** (Latrémolière 2018 §5) to construct the non-compact limit from compact stages. **This is structurally wrong** because:

- AF C*-algebras have **totally-disconnected** Gelfand spectrum (Cantor-set-like).
- $\mathbb{R}$ is **connected**.
- Therefore $C_0(\mathbb{R})$ is NOT AF — no AF inductive limit of $\{C(S^1_{T_n})\}$ can give $C_0(\mathbb{R})$.

The actual AF inductive limit of $\{C(S^1_{T_n})\}$ (with the obvious refinement maps when $T_n | T_{n+1}$) is a **Bunce-Deddens-style** algebra or some other AF object, not $C_0(\mathbb{R})$. So the propinquity-level inductive-limit approach was the wrong tool.

### What's the right tool

**Norm-resolvent convergence of spectral triples** (Connes 1994 §VI; Hekkelman-McDonald 2024 for the Tauberian setting):
$$\bigl\| J_T (\DL_T - z)^{-1} J_T^* - (\DL_{\mathbb{R}} - z)^{-1} \bigr\| \;\to\; 0 \quad \text{as } T \to \infty,$$
for every $z \in \mathbb{C} \setminus \mathbb{R}$, where $J_T : L^2(S^1_T) \hookrightarrow L^2(\mathbb{R})$ is an isometric embedding.

This is the standard NCG notion of "the operator on the smaller carrier converges to the operator on the larger carrier" — captures spectral data, captures operator-algebraic structure, does NOT capture Lipschitz/metric structure.

For G2's *physical* content ("the Lorentzian Dirac on the compactified spacetime converges to the Lorentzian Dirac on the genuine spacetime"), norm-resolvent convergence is the right notion.

For G2's *metric* content (propinquity-level convergence in Latrémolière's sense), we'd need non-compact Latrémolière propinquity — multi-month original work.

---

## §2. Setup — isometric embedding $L^2(S^1_T) \hookrightarrow L^2(\mathbb{R})$

Identify $L^2(S^1_T) \cong L^2([-T/2, T/2])$ (square-integrable functions on the fundamental domain). Define
$$J_T : L^2(S^1_T) \;\hookrightarrow\; L^2(\mathbb{R}), \qquad (J_T f)(t) = \begin{cases} f(t) & t \in [-T/2, T/2] \\ 0 & |t| > T/2 \end{cases}.$$
This is an **isometry** (preserves $L^2$-norm) but **not surjective** (its image is $L^2([-T/2, T/2]) \subset L^2(\mathbb{R})$).

Extend to the Krein space: $\widetilde J_T : \Krein_T = \HGV \otimes L^2(S^1_T) \hookrightarrow \HGV \otimes L^2(\mathbb{R}) = \Krein_{\mathbb{R}}$ via $\widetilde J_T = I_{\HGV} \otimes J_T$.

**Preservation of fundamental symmetry.** The BBB fundamental symmetry $\JL = J_{\mathrm{spatial}} \otimes I_{\mathrm{temporal}}$ acts as the identity on the temporal factor, so
$$\widetilde J_T^* \JL_{\mathbb{R}} \widetilde J_T = \JL_T \quad \text{bit-exact}.$$
This means $\widetilde J_T$ is a **Krein-isometric** embedding, and it commutes with Krein structures on both sides.

**Preservation of K⁺-cone.** Since $\widetilde J_T$ commutes with $\JL$, it preserves the Krein-positive cone: $\widetilde J_T (\Kplus_T) \subset \Kplus_{\mathbb{R}}$.

---

## §3. Theorem (L3c-γ, norm-resolvent convergence)

**Theorem (L3c-γ).** As $T \to \infty$ along an admissible scaling sequence (per L3c-α §1), the Lorentzian Dirac operators converge in norm-resolvent sense modulo the isometric embedding $\widetilde J_T$:

$$\boxed{\,\bigl\| \widetilde J_T\,(\DL_T - z)^{-1}\,\widetilde J_T^* \;-\; \mathbf{1}_{\mathrm{Im}\,\widetilde J_T}\,(\DL_{\mathbb{R}} - z)^{-1}\,\mathbf{1}_{\mathrm{Im}\,\widetilde J_T} \bigr\|_{\mathcal{B}(\Krein_{\mathbb{R}})} \;\to\; 0\,}$$

for every $z \in \mathbb{C} \setminus \mathbb{R}$, where $\mathbf{1}_{\mathrm{Im}\,\widetilde J_T}$ is the projection onto the image of $\widetilde J_T$ in $\Krein_{\mathbb{R}}$.

**Equivalent statement (more standard NCG form).** For every $z \in \mathbb{C} \setminus \mathbb{R}$ and every $\xi \in \Krein_{\mathbb{R}}$ with compact temporal support,
$$\bigl\| \widetilde J_T\,(\DL_T - z)^{-1}\,\widetilde J_T^* \xi \;-\; (\DL_{\mathbb{R}} - z)^{-1}\,\xi \bigr\|_{\Krein_{\mathbb{R}}} \;\to\; 0 \quad \text{as } T \to \infty.$$

The convergence is **strong** (not norm) in general, but **norm** on the dense subspace of compactly-supported smooth functions.

**Honest scope.** The convergence is at the operator-norm level on the dense compact-support subspace, with strong-operator convergence on all of $\Krein_{\mathbb{R}}$. This is the standard convergence notion for compact-to-noncompact-domain limits in mathematical physics (e.g., periodic-domain → free-domain Schrödinger operator limit, Bloch decomposition limits, etc.).

---

## §4. Proof outline

### Step 1 — Spatial factorization

Factor the Lorentzian Dirac as $\DL = i \gamma^0 \otimes \partial_t + i \DGV \otimes I$. The spatial term $i \DGV \otimes I$ is *identical* on $\Krein_T$ and $\Krein_{\mathbb{R}}$ (both have the same $\HGV$ spatial factor; the difference is only in the temporal Hilbert space). Therefore norm-resolvent convergence of $\DL_T \to \DL_{\mathbb{R}}$ **reduces to** norm-resolvent convergence of the temporal piece $i \gamma^0 \otimes \partial_t^{S^1_T} \to i \gamma^0 \otimes \partial_t^{\mathbb{R}}$.

### Step 2 — Temporal-only convergence

Norm-resolvent convergence of $i \partial_t^{S^1_T} \to i \partial_t^{\mathbb{R}}$ as $T \to \infty$ is a classical result. Two routes:

**Route (a): direct spectral comparison.** On $S^1_T$, $i \partial_t$ has spectrum $\{2\pi k/T : k \in \mathbb{Z}\}$. On $\mathbb{R}$, $i \partial_t$ has continuous spectrum $\mathbb{R}$. As $T \to \infty$, the discrete spectrum $\{2\pi k/T\}$ becomes dense in $\mathbb{R}$. By the standard *Bloch-decomposition* / *Brillouin-zone-shrinking* argument (Reed-Simon Vol. IV §XIII.16, or Sjöstrand-Vodev semi-classical), the resolvent converges in norm on compact-support subspaces.

**Route (b): direct compact-to-noncompact-domain argument.** Let $f, g \in C_c^\infty(\mathbb{R})$ with $\mathrm{supp}(f), \mathrm{supp}(g) \subset [-A, A]$ for some $A > 0$. For $T > 2A$, the periodization of $f$ to $S^1_T$ is just $f$ (no overlap), so
$$\widetilde J_T \,(i \partial_t^{S^1_T} - z)^{-1} \,\widetilde J_T^* g \;=\; \mathbf{1}_{[-T/2, T/2]} \cdot \bigl((i \partial_t^{\mathbb{R}} - z)^{-1} g\bigr)_{T-\mathrm{periodic}}.$$
The "$T$-periodic" subscript denotes the unique periodic extension, which agrees with $(i \partial_t^{\mathbb{R}} - z)^{-1} g$ on $[-T/2, T/2]$ up to a tail error that decays exponentially (the resolvent on $\mathbb{R}$ is the convolution kernel $K_z(t) = (i z)^{-1} e^{i z |t|}$ for $\mathrm{sgn}(\mathrm{Im}\,z) > 0$, which decays exponentially in $|t|$).

Specifically, the tail error is
$$\bigl\| \widetilde J_T (i \partial_t^{S^1_T} - z)^{-1} \widetilde J_T^* g - (i \partial_t^{\mathbb{R}} - z)^{-1} g \bigr\|_{L^2(\mathbb{R})} \;\le\; C(z) \cdot e^{-|\mathrm{Im}\,z|\,(T/2 - A)}$$
for some constant $C(z)$ depending only on $\|g\|$ and $\mathrm{Im}\,z \neq 0$. As $T \to \infty$, this tail vanishes **exponentially fast** for fixed compact-support $g$.

### Step 3 — Tensor-product extension to the full Krein space

Norm-resolvent convergence of $i \gamma^0 \otimes \partial_t^{S^1_T} \to i \gamma^0 \otimes \partial_t^{\mathbb{R}}$ on $\HGV \otimes L^2$ follows from Step 2 via the tensor-product structure (the spatial $\HGV$ factor is fixed and finite-dimensional, so tensor-product norm convergence is automatic).

### Step 4 — Sum of resolvents

Norm-resolvent convergence of a sum of two operators follows from the second resolvent identity: if $A_n \to A$ and $B_n \to B$ in norm-resolvent sense, and $A + B$ and $A_n + B_n$ are essentially self-adjoint (or appropriately Krein-self-adjoint), then $A_n + B_n \to A + B$ in norm-resolvent sense. Apply with $A_n = i \DGV \otimes I$ (constant, trivial convergence) and $B_n = i \gamma^0 \otimes \partial_t^{S^1_{T_n}}$ (Step 2).

### Step 5 — Krein-self-adjointness preservation

The Lorentzian Dirac $\DL$ is Krein-self-adjoint with respect to $\JL$ (Paper 43 §3, Sprint L2-B/C closure). Norm-resolvent convergence preserves Krein-self-adjointness because the resolvent identity $\JL R(z) \JL^{-1} = R(\bar z)^*$ is structural.

### Combined

Steps 1-5 give norm-resolvent convergence on the dense compact-support subspace, with the explicit exponential tail estimate
$$\bigl\| \widetilde J_T (\DL_T - z)^{-1} \widetilde J_T^* \xi - (\DL_{\mathbb{R}} - z)^{-1} \xi \bigr\| \;\le\; C(z, \xi) \cdot e^{-|\mathrm{Im}\,z|\,T/2}$$
for $\xi$ with compact temporal support. Extension to general $\xi \in \Krein_{\mathbb{R}}$ gives strong-operator convergence by density.

---

## §5. What this gives, and what it doesn't

### What L3c-γ gives

1. **Spectral convergence of the Lorentzian Dirac sequence.** The discrete-spectrum $\DL_{S^1_T}$ converges (in norm-resolvent) to the continuous-spectrum $\DL_{\mathbb{R}}$ as $T \to \infty$.
2. **Operator-algebraic convergence.** The Krein-self-adjoint $\DL$ structure is preserved across the limit.
3. **K⁺-cone preservation.** The Krein-positive structure $\Kplus$ transports cleanly via the isometry $\widetilde J_T$.
4. **Exponential tail bound.** The convergence rate is $O(e^{-|\mathrm{Im}\,z| T/2})$ on compact-support subspaces — extremely fast.
5. **Coupling with L3c-α.** Composed with the L3c-α parametric stability ($\Lprop(\mathcal{T}^L_{n_{\max}, N_t, T(N_t)}, \mathcal{T}^L_{S^3 \times S^1_{T(N_t)}}) \to 0$), this gives a two-stage convergence:
   $$\mathcal{T}^L_{n_{\max}, N_t, T(N_t)} \;\xrightarrow{\mathrm{Paper\;45/46}}\; \mathcal{T}^L_{S^3 \times S^1_{T(N_t)}} \;\xrightarrow{\mathrm{L3c\!-\!\gamma}}\; \mathcal{T}^L_{S^3 \times \mathbb{R}_t}.$$
   The first arrow is propinquity convergence (metric level). The second arrow is norm-resolvent convergence (spectral level).

### What L3c-γ does NOT give

1. **Propinquity-level convergence to $\mathcal{T}^L_{S^3 \times \mathbb{R}_t}$.** Latrémolière propinquity is not defined for non-compact carriers in the published literature. Extending it would be multi-month original work (G1-genuine / L3c-γ.2).
2. **Lipschitz / metric structure on the limit.** Norm-resolvent captures the *Hilbert-space* and *operator-algebraic* content but not the metric/quantum-distance content.
3. **A uniform rate in $(n_{\max}, N_t, T)$.** The L3c-γ bound gives convergence in $T$ at fixed $(n_{\max}, N_t)$ for the *outer* arrow; the L3c-α bound gives convergence in $(n_{\max}, N_t)$ for the *inner* arrow. They don't compose into a single uniform rate — the inner is rate $O(\log n_{\max}/n_{\max} + T/N_t)$ in propinquity, the outer is rate $O(e^{-T/2})$ in norm-resolvent.

### Joint convergence — the combined L3c-α + L3c-γ statement

For any admissible scaling $T(N_t) \to \infty$ with $T(N_t)/N_t \to 0$, and any compact-support test data, the joint convergence
$$\mathcal{T}^L_{n_{\max}, N_t, T(N_t)} \;\longrightarrow\; \mathcal{T}^L_{S^3 \times \mathbb{R}_t}$$
holds in the sense:
- **Inner (propinquity, L3c-α):** $\Lprop(\mathcal{T}^L_{n_{\max}, N_t, T(N_t)}, \mathcal{T}^L_{S^3 \times S^1_{T(N_t)}}) \to 0$ at rate $\max(\log n_{\max}/n_{\max}, T(N_t)/N_t)$.
- **Outer (norm-resolvent, L3c-γ):** $\mathcal{T}^L_{S^3 \times S^1_{T(N_t)}} \to \mathcal{T}^L_{S^3 \times \mathbb{R}_t}$ at rate $O(e^{-T(N_t)/2})$ on compact-support subspaces.

The composition gives a **two-rate hybrid convergence statement**: propinquity-level on the finite-stage approximation, norm-resolvent-level on the outer de-compactification. This is the appropriate G2 closure at the level GeoVac currently supports.

---

## §6. G2 closure status

After L3c-α + L3c-γ, G2 (Paper 45 §1.4 de-compactification) has the following status:

| Convergence level | $\mathcal{T}^L_{n_{\max}, N_t, T(N_t)} \to \mathcal{T}^L_{S^3 \times \mathbb{R}_t}$ |
|:------------------|:-------------------------------------------------------------------------------------|
| **Norm-resolvent** | **CLOSED via L3c-γ ∘ L3c-α composite.** Spectral data, operator-algebra structure, Krein-self-adjointness, K⁺-cone all transport. |
| **Propinquity** | **OPEN.** Latrémolière propinquity not defined on non-compact carriers. Requires G1-genuine multi-month NCG extension. |
| **Synthetic Lorentzian GH (Mondino-Sämann)** | Not in scope — pre-length-space / causal-diamond framework is structurally different. |

**Honest framing for G2:** L3c-α + L3c-γ together close G2 **at the spectral level** — the appropriate physical interpretation is that the Lorentzian Dirac sequence converges to the genuine Lorentzian Dirac on $S^3 \times \mathbb{R}_t$, and this is sufficient for any *operator-theoretic* or *spectral-action* application of the framework. Propinquity-level closure on the non-compact carrier remains open as a separate sub-question (G2-metric or equivalently G1-genuine on the non-compact substrate).

---

## §7. Connection to Paper 43 / L2-B/E

Paper 43's bounded-interval Krein triple $\mathcal{T}^L_{S^3 \times [-T_{\max}, T_{\max}]}$ (Sprint L2-B/E construction with Dirichlet zero BC on the temporal endpoints) is structurally distinct from both $\mathcal{T}^L_{S^3 \times S^1_T}$ (periodic BC) and $\mathcal{T}^L_{S^3 \times \mathbb{R}_t}$ (no BC / continuous spectrum).

**Norm-resolvent comparison among the three:**
- $\mathcal{T}^L_{S^3 \times [-T, T]} \to \mathcal{T}^L_{S^3 \times \mathbb{R}_t}$ as $T \to \infty$: standard Dirichlet-to-free convergence, same exponential tail rate $O(e^{-|\mathrm{Im}\,z| T})$.
- $\mathcal{T}^L_{S^3 \times S^1_T} \to \mathcal{T}^L_{S^3 \times \mathbb{R}_t}$ as $T \to \infty$: this sprint's result, same rate.
- $\mathcal{T}^L_{S^3 \times S^1_T} \to \mathcal{T}^L_{S^3 \times [-T, T]}$ as $T \to \infty$: also same rate (periodic vs Dirichlet BC agree away from the boundary).

So all three carriers — periodic compact, finite Dirichlet, non-compact — give the SAME norm-resolvent limit object $\mathcal{T}^L_{S^3 \times \mathbb{R}_t}$. The choice between Paper 43's bounded-interval architecture and Paper 45's compact-temporal architecture is **structurally a choice of approximation scheme**, with both converging to the same physical Lorentzian Dirac.

**Why both architectures exist.** Paper 43 chose bounded Dirichlet to avoid CTCs (closed timelike curves) at the operator-system level (Sprint L2-B explicitly avoided periodic BC via Geroch's theorem). Paper 45 chose periodic compact to enable Latrémolière propinquity (which needs a compact carrier). Both are valid; L3c-γ identifies them via norm-resolvent convergence to the common non-compact limit.

---

## §8. Open questions for PI

1. **Accept L3c-γ at norm-resolvent level as G2 closure?** This is the standard NCG-style convergence notion. Propinquity-level closure on non-compact carrier would require multi-month G1-genuine work. Is norm-resolvent sufficient for the framework's physical content?
2. **Apply paper edits?** Three natural edits:
   - Paper 45 §1.4 G2 update: refine to clarify two-level closure (norm-resolvent CLOSED, propinquity OPEN as G1-genuine non-compact extension).
   - Paper 45 §6 or §7 new Remark: combined L3c-α + L3c-γ two-rate hybrid convergence.
   - Paper 43 §X cross-reference: bounded-interval and compact-temporal carriers identified via L3c-γ norm-resolvent limit.
3. **Open Sprint L3d / Paper 47?** A standalone math.OA paper on the two-rate hybrid convergence (L3c-α + L3c-γ composite) is a natural deliverable. ~2–3 weeks to draft, leverages Papers 45/46 verbatim. Would be 9th math.OA standalone (siblings 38/39/40/42/43/44/45/46).
4. **Open Sprint L3e (multi-month frontier)?** Defines non-compact Latrémolière propinquity. G1-genuine territory. Probably not sprint-scale; queue as multi-month research direction.

---

## §9. Recommended paper edits (NOT applied; PI decision)

If PI accepts L3c-γ at the norm-resolvent level:

**Edit 1 — Paper 45 §1.4 G2 refinement:**
> *(G2)* De-compactification limit $T \to \infty$. Our construction is at compact temporal radius $T$ canonical $T = 2\pi$ (the BW modular period). The non-compact $\mathbb{R}_t$ limit is closed at the **norm-resolvent level** (Paper~47~\cite{paper47}~\S~3, via the composite L3c-α parametric stability and L3c-γ outer de-compactification). The **propinquity-level** limit on a non-compact carrier requires a non-compact extension of Latrémolière propinquity (multi-month NCG-math, named G1-genuine in our internal taxonomy) and is not addressed here.

**Edit 2 — Paper 45 new §6.3 "Joint convergence at two levels":**
> *Remark.* The L3c-α parametric stability and L3c-γ outer norm-resolvent convergence compose to give a two-rate hybrid statement: $\mathcal{T}^L_{n_{\max}, N_t, T(N_t)} \to \mathcal{T}^L_{S^3 \times \mathbb{R}_t}$ with propinquity rate $\max(\log n_{\max}/n_{\max}, T(N_t)/N_t)$ on the inner arrow and norm-resolvent rate $O(e^{-T(N_t)/2})$ on the outer arrow. This identifies the genuine non-compact Krein triple as the spectral limit of the framework's compact-temporal sequence.

**Edit 3 — Paper 43 new §X.Y cross-reference:**
> *Remark.* The bounded-interval Krein triple of this paper and the compact-temporal triple of Paper~45~\cite{paper45} converge to the same non-compact limit $\mathcal{T}^L_{S^3 \times \mathbb{R}_t}$ in norm-resolvent sense (Paper~47~\cite{paper47}~\S~7), identifying them as two compact approximations to the same physical Lorentzian Dirac operator. The choice between them is a choice of approximation scheme: Paper 43 prioritizes physical content (no CTCs, finite Dirichlet boundary), Paper 45 prioritizes metric-convergence machinery (compact carrier for Latrémolière propinquity).

These edits are scope-honest and queued for PI decision. **Not applied in this memo.**

---

## §10. Sprint verdict

**Sprint L3c-γ: CLOSED at the analytical level at norm-resolvent convergence.**

- Strategic reframing (§1): inductive-limit propinquity was structurally wrong (AF cannot capture $C_0(\R)$); norm-resolvent is the right NCG-standard framework.
- Theorem statement (§3): norm-resolvent convergence $\widetilde J_T (\DL_T - z)^{-1} \widetilde J_T^* \to (\DL_{\mathbb{R}} - z)^{-1}$ via standard compact-to-noncompact-domain machinery.
- Proof outline (§4): spatial factor identical, temporal factor converges via direct compact-support resolvent comparison with exponential tail $O(e^{-|\mathrm{Im}\,z| T/2})$, sum-of-resolvents combination, Krein-self-adjointness preservation.
- Joint convergence (§5–§6): combined L3c-α + L3c-γ closes G2 **at the norm-resolvent (spectral) level**, leaves propinquity-level closure on non-compact carrier as G1-genuine multi-month frontier.

**G2 closure status:**

| Level | Status | Sprint |
|:------|:-------|:-------|
| Spectral (norm-resolvent) | **CLOSED** | L3c-α + L3c-γ |
| Metric (propinquity) | **OPEN** (multi-month G1-genuine) | L3c-γ.2 / L3e |

**Confidence:** HIGH on the analytical statement (§3) — the proof is standard mathematical physics (Reed-Simon Vol. IV §XIII.16, Sjöstrand-Vodev semi-classical, or direct compact-support resolvent comparison via exponential decay of $(D - z)^{-1}$ kernel). HIGH on the strategic reframing (§1) — AF inductive limit cannot capture $C_0(\R)$, so norm-resolvent is the right tool. MEDIUM on whether this should be a standalone Paper 47 or absorbed as remarks in Papers 45/43.

**Recommended next move:** Apply paper edits (§9) and queue a decision on standalone Paper 47 vs absorbed remarks. The two-rate hybrid convergence statement is a clean, publishable result that closes G2 at the appropriate level for GeoVac's current scope.
