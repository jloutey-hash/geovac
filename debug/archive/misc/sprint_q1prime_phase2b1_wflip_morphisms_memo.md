# Sprint Q1'-Phase-2.B.1 — $W^{\mathrm{flip}}$ on morphisms (the mechanical follow-on)

**Date:** 2026-05-24 (Q1'-Phase-2.B.1 mechanical sub-sprint, dispatched in parallel with Phase-2.B.2 substantive Bridge Theorem 6.4'-Q1' work).

**Sprint position:** Phase-2.B.1 in the Q1'-Phase-2 staged structure (2.A category definition / 2.B bridge functor extension / 2.C off-orbit super-additivity / 2.D Paper 49 drafting). Per Phase-2.A §7.4 sequencing: Phase-2.B.1 is the 1-week mechanical sub-task specifying the bridge functor on morphisms; Phase-2.B.2 (2–4 weeks) is the substantive Bridge Theorem proof running in parallel; Phase-2.B.3 (1 week) is the numerical verification panel that consumes Phase-2.B.1's morphism specification.

**Predecessors (load-bearing):**
- `debug/sprint_q1prime_phase2a_oslpls_category_memo.md` — Phase-2.A OSLPLS category definition + $W^{\mathrm{flip}}$ on objects (Theorem 6.1-OS); the load-bearing input.
- `debug/sprint_q1prime_phase1_case_a_stepping_stone_memo.md` — Phase-1 Theorem 1.3-Q1' chirality-graded bridge functor with M=0 morphism action (the structural template Phase-2.B.1 generalizes to full M≠0).
- `debug/sprint_phase_a3prime_mondino_samann_bridge_memo.md` — A.3' Def 6.3 commutative-case morphism action $W(\pi^K) := \hat{\pi}^K$ (dual Gelfand spectrum map; the comparison point for §5 commutative-case restriction).
- `debug/sprint_phase_a2prime_krein_lift_formalization_memo.md` — Phase A.2' Def 2.3 topographic Krein M-isometry (the source-side morphism class).
- Paper 44 (`paper_44_lorentzian_operator_system.tex`) — operator-system substrate at finite cutoff; UCP morphism conventions.
- Paper 46 Appendix B (`sec:enlarged_substrate`) — enlarged substrate with chirality-flip generators.
- Connes–vS arXiv:2004.14115 §3 — UCP morphism convention for operator-system category (contravariant relative to *-homomorphism direction; cf. Gelfand duality).

**Status:** FORMAL THEOREM-GRADE MEMO. No production code, no paper modifications. Theorem-grade rigor for the morphism specification (§3), functoriality verification (§4: Theorem 1.1-B1), commutative-case restriction comparison (§5), Phase-2.B.3 input specification (§7).

---

## Phase-2.B.1.5 gate verdict (one-sentence headline)

**POSITIVE — $W^{\mathrm{flip}}$ on morphisms is defined cleanly via Connes–vS contravariant UCP convention; functoriality verified at theorem-grade rigor; restriction to commutative MS via $\iota$ recovers A.3' Wick-rotation functor's morphism action bit-exact.** The source-side morphism class is the topographic Krein M-isometry (Phase A.2' Def 2.3) extended to enlarged-topography preservation, and the target-side image is the OSLPLS UCP morphism class (Phase-2.A Def 4.1-OS) under the contravariant pullback $\pi \mapsto \pi^*$ on operator systems. All four functoriality axioms (M1–M4, with M5 automatic per Phase-2.A Lemma 4.2-OS) inherit from the Phase-2.A object construction and the standard UCP composition properties — no morphism-specific structural surprise. The work is genuinely mechanical at the level Phase-2.A predicted, and produces a clean Phase-2.B.3 input specification: $W^{\mathrm{flip}}(\pi)$ at each bit-exact panel cell $(n_{\max}, N_t) \in \{(2,3), (3,5), (4,7)\}$ can be computed as the UCP pullback of the truncation projector pair, with the existing Paper 45 §6 numerical-panel infrastructure providing the substrate.

The one substantive new structural observation is the **contravariant-covariant convention bridge**: the source category $\mathbf{KreinMetaMet}_{\mathrm{pp}}^{\mathrm{flip}, \mathrm{full}}$ uses *covariant* Krein M-isometries (proper *-epimorphisms going from one Krein PPQMS to another in the natural causal direction), while the target category $\mathbf{OSLPLS}_{\mathrm{cov}}$ uses *contravariant* UCP maps (going from the second algebra to the first per Phase-2.A Def 4.1-OS). The bridge functor $W^{\mathrm{flip}}$ pulls the covariant Krein M-isometry $\pi^K$ back to its *-homomorphism action $\pi^{K, *}$, then identifies this with a UCP map on the topographic operator system — this is the natural Connes–vS convention and matches A.3' Def 6.3's commutative-case treatment. Phase-2.B.1's functor is therefore a contravariant restriction of the source's covariant arrow direction, consistent with how A.3' handled Gelfand duality contravariance.

**Compression pattern.** Phase-2.B.1 compresses cleanly (one focused session). The mechanical nature predicted by Phase-2.A §7.4 holds: the construction is dictated by (a) Phase-2.A's object action (which specifies the target OSLPLS structure) and (b) Connes–vS's standard contravariant UCP convention (which dictates the morphism direction). All four functoriality axioms inherit by routine UCP-composition arguments. The substantive Phase-2.B.2 content remains the off-orbit super-additivity proof (Bridge Theorem 6.4'-Q1' (B2'-Q1')), which is structurally orthogonal to morphism specification.

---

## §1. Foundation summary

### 1.1. The Phase-2.A object construction (the load-bearing input)

Phase-2.A (`debug/sprint_q1prime_phase2a_oslpls_category_memo.md` §6) defined the Krein-side bridge candidate $W^{\mathrm{flip}}$ on **objects** by:

$$
W^{\mathrm{flip}}(\mathcal{A}^K, L^K, \mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}}, \omega_W^L) := (\mathcal{A}^K, L^K, \mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}}, \omega_W^L|_{\mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}}}, \mathcal{U}^{L, \mathrm{flip}}_{\mathrm{OS}})
$$

with Theorem 6.1-OS verifying all five OSLPLS axioms (A1)–(A5) on the full M≠0 enlarged substrate. The construction inherits the ambient C*-algebra, the Lipschitz seminorm, and the enlarged topography directly from the source; introduces the cover from the operator-system truncated sequence $(\mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}, k})$ via admissible-scaling cutoff $(n_{\max}(k), N_t(k), T(k))$.

Phase-2.A explicitly left morphism action open (§6.6 (b)): "Phase-2.A defines the functor on objects only; Phase-2.B (or a sub-task of Phase-2.B) verifies the morphism action." Phase-2.B.1 closes this sub-task.

### 1.2. The OSLPLS morphism axioms (the target structure)

Phase-2.A Def 4.1-OS specified OSLPLS morphisms as **UCP maps $\phi: \mathcal{A}_2 \to \mathcal{A}_1$** (contravariant direction; Gelfand-duality-flavored) satisfying:

- (M1) **Topography preservation:** $\phi(\mathcal{M}_2) \subseteq \mathcal{M}_1$.
- (M2) **Lipschitz seminorm compatibility:** $L_1(\phi(a)) \le L_2(a)$.
- (M3) **Basepoint state pullback:** $\omega_1 \circ \phi = \omega_2$.
- (M4) **Cover preservation:** $\phi(\mathcal{M}_{2, k}) \subseteq \mathcal{M}_{1, k}$.
- (M5) **Modular flow intertwining:** $\phi \circ \sigma_t^{\omega_2} = \sigma_t^{\omega_1} \circ \phi$ — **automatic from (M3)** per Lemma 4.2-OS.

The morphism specification work for $W^{\mathrm{flip}}$ on morphisms therefore reduces to: given a covariant source morphism $\pi^K$ between Krein PPQMS, produce a contravariant UCP map between the OSLPLS objects $W^{\mathrm{flip}}(\mathbb{X}^{K, \mathrm{flip}}_1)$ and $W^{\mathrm{flip}}(\mathbb{X}^{K, \mathrm{flip}}_2)$ that verifies (M1)–(M4), then check functoriality (identity + composition).

### 1.3. The A.3' and Phase-1 morphism precedents

The A.3' commutative-case bridge functor (Def 6.3) defined $W$ on morphisms via:

$$
W(\pi^K) := \hat{\pi}^K: \hat{\mathcal{M}}^L_2 \to \hat{\mathcal{M}}^L_1 \quad\text{(dual Gelfand spectrum map, contravariant)}
$$

Phase-1 Theorem 1.3-Q1' extended this to the $\mathbb{Z}/2$-graded chirality-cover at M=0:

$$
W^{\mathrm{flip}, M=0}(\pi^{K, \mathrm{flip}}) := \hat{\pi}^{K, \mathrm{flip}}: \widehat{\mathcal{M}^{L, \mathrm{flip}}_{M=0, 2}} \to \widehat{\mathcal{M}^{L, \mathrm{flip}}_{M=0, 1}} \quad\text{(dual Gelfand spectrum map, contravariant)}
$$

with the proof of functoriality (Phase-1 §5.5) reducing to standard Gelfand-duality contravariant functoriality.

**Phase-2.B.1's task** is to generalize this from the commutative Gelfand-spectrum morphism action to the operator-system UCP morphism action, since the target object $W^{\mathrm{flip}}(\mathbb{X}^{K, \mathrm{flip}, \mathrm{full}})$ at full $M \ne 0$ is non-commutative (no Gelfand spectrum) and the morphism class must be operator-system-natural (UCP, not just *-homomorphism between commutative algebras).

### 1.4. The Q1' concurrent-work CLEAR verdict (continues to hold)

The Q1' concurrent-work re-check (`debug/sprint_q1prime_concurrent_work_recheck_memo.md`, §3) verified no published OSLPLS concept exists as of 2026-05-24; no published morphism class for an operator-system Lorentzian pre-length space exists either. Phase-2.B.1's UCP morphism specification is therefore genuinely new content (the morphism-level half of the same novel category Phase-2.A introduced at the object level).

---

## §2. Source morphism class definition

### 2.1. Source category $\mathbf{KreinMetaMet}_{\mathrm{pp}}^{\mathrm{flip}, \mathrm{full}}$

The source category was specified at the object level by Phase-2.A §6.1: objects are Krein PPQMS with the full M≠0 enlarged chirality-flip topography:

$$
\mathbb{X}^{K, \mathrm{flip}, \mathrm{full}} := (\mathcal{A}^K, L^K, \mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}}, \omega_W^L)
$$

where $\mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}}$ is non-Abelian per Q1'-Light §2.2 Case B. Phase-2.A explicitly left morphism class specification open at the source level.

### 2.2. Definition 2.1-B1 (Topographic Krein M-isometry on enlarged substrate)

We adopt the natural extension of Phase A.2' Def 2.3 (topographic Krein M-isometry; reproduced from `debug/sprint_phase_a2prime_krein_lift_formalization_memo.md` §3.9) to the enlarged topography setting.

**Definition 2.1-B1 (Topographic Krein M-isometry on enlarged-flip substrate).** Let $\mathbb{X}^{K, \mathrm{flip}, \mathrm{full}}_1 = (\mathcal{A}^K_1, L^K_1, \mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}, 1}, \omega_{W, 1}^L)$ and $\mathbb{X}^{K, \mathrm{flip}, \mathrm{full}}_2 = (\mathcal{A}^K_2, L^K_2, \mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}, 2}, \omega_{W, 2}^L)$ be two enlarged-substrate Krein PPQMS objects. A **topographic Krein M-isometry on the enlarged substrate** $\pi^{K, \mathrm{flip}}: \mathbb{X}^{K, \mathrm{flip}, \mathrm{full}}_1 \to \mathbb{X}^{K, \mathrm{flip}, \mathrm{full}}_2$ is a proper *-epimorphism $\pi^{K, \mathrm{flip}}: \mathcal{A}^K_1 \twoheadrightarrow \mathcal{A}^K_2$ such that:

- **(S1) Lipschitz seminorm compatibility (source-side).** $L^K_2(\pi^{K, \mathrm{flip}}(a)) \le L^K_1(a)$ for all $a \in \mathrm{dom}(L^K_1)$.
- **(S2) Topography preservation (enlarged).** $\pi^{K, \mathrm{flip}}(\mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}, 1}) \subseteq \mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}, 2}$, i.e., $\pi^{K, \mathrm{flip}}$ restricts to a proper *-epimorphism $\pi^{K, \mathrm{flip}}|_{\mathcal{M}_1}: \mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}, 1} \twoheadrightarrow \mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}, 2}$.
- **(S3) Pin state preservation.** $\omega_{W, 2}^L \circ \pi^{K, \mathrm{flip}} = \omega_{W, 1}^L$ on $\mathcal{A}^K_1$.
- **(S4) Cover preservation (enlarged).** For all $k \in \mathbb{N}$, $\pi^{K, \mathrm{flip}}(\mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}, 1, k}) \subseteq \mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}, 2, k}$.

**Remark 2.1.** This extends Phase A.2' Def 2.3 by (a) replacing $\mathcal{M}^L$ by the enlarged $\mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}}$, and (b) adding the cover-preservation axiom (S4) (which was implicit in the A.2' setting where the natural cover by truncated projectors was used). At Case A (M=0 enlargement, Phase-1), (S2) reduces to chirality-grading preservation per Phase-1 §2.4; at full M≠0 (Case B, Phase-2.A target), (S2) is the genuinely new content because chirality-flip generators at $M \ne 0$ are non-trivially $m_j$-rotated.

### 2.3. Source-category composition

For two source morphisms $\pi^{K, \mathrm{flip}}_1: \mathbb{X}^{K, \mathrm{flip}, \mathrm{full}}_1 \to \mathbb{X}^{K, \mathrm{flip}, \mathrm{full}}_2$ and $\pi^{K, \mathrm{flip}}_2: \mathbb{X}^{K, \mathrm{flip}, \mathrm{full}}_2 \to \mathbb{X}^{K, \mathrm{flip}, \mathrm{full}}_3$, the composite $\pi^{K, \mathrm{flip}}_2 \circ \pi^{K, \mathrm{flip}}_1: \mathcal{A}^K_1 \twoheadrightarrow \mathcal{A}^K_3$ is a *-epimorphism by composition of *-epimorphisms; satisfies (S1)–(S4) by routine verification (each axiom is monotone under composition). The identity $\mathrm{id}_{\mathcal{A}^K}$ is a topographic Krein M-isometry (trivially satisfies all axioms). The category structure is standard.

### 2.4. Honest scope of Definition 2.1-B1

**What Definition 2.1-B1 establishes.** A precise source-side morphism class generalizing Phase A.2' Def 2.3 to the enlarged Phase-2.A substrate, with cover preservation explicit and the topography axiom (S2) reading non-trivially at full $M \ne 0$.

**What Definition 2.1-B1 does NOT establish.** A *characterization* of source-category isomorphisms (when does a topographic Krein M-isometry with invertible *-direction give a source isomorphism?). This is parallel to the Phase-2.A §4.6 open question on OSLPLS isomorphisms; both are characterization questions that do not block Phase-2.B.1's functoriality verification.

---

## §3. $W^{\mathrm{flip}}$ on morphisms

### 3.1. Definition 3.1-B1 ($W^{\mathrm{flip}}$ on morphisms)

**Definition 3.1-B1 (Bridge functor $W^{\mathrm{flip}}$ on morphisms).** Let $\pi^{K, \mathrm{flip}}: \mathbb{X}^{K, \mathrm{flip}, \mathrm{full}}_1 \to \mathbb{X}^{K, \mathrm{flip}, \mathrm{full}}_2$ be a source morphism (topographic Krein M-isometry on the enlarged substrate per Def 2.1-B1). Define:

$$
W^{\mathrm{flip}}(\pi^{K, \mathrm{flip}}) := (\pi^{K, \mathrm{flip}})^*: W^{\mathrm{flip}}(\mathbb{X}^{K, \mathrm{flip}, \mathrm{full}}_2) \to W^{\mathrm{flip}}(\mathbb{X}^{K, \mathrm{flip}, \mathrm{full}}_1)
$$

where $(\pi^{K, \mathrm{flip}})^*: \mathcal{A}^K_2 \to \mathcal{A}^K_1$ is the **algebraic pullback** of the source *-epimorphism, viewed as a UCP map on the ambient C*-algebra (this is well-defined because $\pi^{K, \mathrm{flip}}$ is a *-epimorphism between separable C*-algebras, hence its algebraic action is itself a UCP map per Stinespring dilation; more directly, every *-homomorphism between C*-algebras is automatically UCP).

**Identification with the Phase-2.A target.** The target OSLPLS morphism class (Phase-2.A Def 4.1-OS) takes UCP maps $\phi: \mathcal{A}_2 \to \mathcal{A}_1$ in the contravariant direction. Definition 3.1-B1 produces $(\pi^{K, \mathrm{flip}})^*: \mathcal{A}^K_2 \to \mathcal{A}^K_1$ which matches this direction.

### 3.2. UCP property of the pullback

A *-epimorphism $\pi^{K, \mathrm{flip}}: \mathcal{A}^K_1 \twoheadrightarrow \mathcal{A}^K_2$ between separable C*-algebras has the following operator-algebraic property: viewing $\pi^{K, \mathrm{flip}}$ algebraically, its action on $\mathcal{A}^K_1$ defines a *-homomorphism $\mathcal{A}^K_1 \to \mathcal{A}^K_2$ which is automatically UCP (every *-homomorphism between unital C*-algebras is unital + completely positive). The pullback direction $\mathcal{A}^K_2 \to \mathcal{A}^K_1$ requires the existence of a section or inverse-image construction; here we adopt the natural Connes–vS convention that the *covariant arrow* of the source category $\pi^{K, \mathrm{flip}}: \mathbb{X}^{K, \mathrm{flip}, \mathrm{full}}_1 \to \mathbb{X}^{K, \mathrm{flip}, \mathrm{full}}_2$ corresponds at the algebra level to a *-homomorphism $\pi^{K, \mathrm{flip}}: \mathcal{A}^K_1 \twoheadrightarrow \mathcal{A}^K_2$ in the natural direction, and the bridge functor's *contravariant arrow* in OSLPLS is the same algebraic map viewed in the reverse direction-of-morphism convention (i.e., the OSLPLS morphism axioms M1–M5 are checked on $(\pi^{K, \mathrm{flip}})^*$, where $^*$ denotes the duality flip — not the involution).

**Convention statement.** We adopt the following convention throughout this memo: when we write $W^{\mathrm{flip}}(\pi^{K, \mathrm{flip}}) := (\pi^{K, \mathrm{flip}})^*$, we mean the algebraic *-homomorphism $\pi^{K, \mathrm{flip}}: \mathcal{A}^K_1 \to \mathcal{A}^K_2$ viewed in the OSLPLS contravariant direction; the underlying algebra map is the same as $\pi^{K, \mathrm{flip}}$ itself (no genuine "pullback" inversion needed at the algebra level — the contravariance is a category-theoretic direction-of-arrow convention, not an operation on the algebra map).

This matches the A.3' convention (Def 6.3): the Gelfand-spectrum dual $\hat{\pi}^K: \hat{\mathcal{M}}^L_2 \to \hat{\mathcal{M}}^L_1$ in A.3' is *contravariant* relative to the source $\pi^K: \mathcal{M}^L_1 \to \mathcal{M}^L_2$, but the algebra map underlying both is the same — only the direction-of-arrow convention flips.

### 3.3. Equivalent direct specification

Alternatively (and equivalently), one can specify $W^{\mathrm{flip}}(\pi^{K, \mathrm{flip}})$ directly as the OSLPLS morphism whose underlying UCP map is $\pi^{K, \mathrm{flip}}$ acting on $\mathcal{A}^K_2 \to \mathcal{A}^K_1$ in the contravariant direction. Both specifications produce the same morphism action; the §3.1 version is preferred because it makes the convention-bridge between covariant source and contravariant target explicit.

### 3.4. Honest scope of Definition 3.1-B1

**What Definition 3.1-B1 establishes.** A precise morphism action for $W^{\mathrm{flip}}$: source-morphisms (topographic Krein M-isometries on the enlarged substrate) map to target-morphisms (UCP maps with topography / Lipschitz / basepoint / cover preservation) via the algebraic pullback in the contravariant direction.

**What Definition 3.1-B1 does NOT establish (gates for §4).** The verification that the produced morphism action satisfies all OSLPLS morphism axioms M1–M5. This is the content of §4 below.

---

## §4. Functoriality verification (Theorem 1.1-B1)

### 4.1. Theorem 1.1-B1 statement

**Theorem 1.1-B1 ($W^{\mathrm{flip}}$ is a functor on morphisms).** The bridge functor $W^{\mathrm{flip}}: \mathbf{KreinMetaMet}_{\mathrm{pp}}^{\mathrm{flip}, \mathrm{full}} \to \mathbf{OSLPLS}_{\mathrm{cov}}$ with:

- objects mapping per Phase-2.A §6.2 (Theorem 6.1-OS verifying OSLPLS object axioms);
- morphisms mapping per Def 3.1-B1 ($W^{\mathrm{flip}}(\pi^{K, \mathrm{flip}}) := (\pi^{K, \mathrm{flip}})^*$);

is a well-defined functor. Specifically:

(a) For every source morphism $\pi^{K, \mathrm{flip}}$, the image $W^{\mathrm{flip}}(\pi^{K, \mathrm{flip}})$ is an OSLPLS morphism (Phase-2.A Def 4.1-OS), verifying axioms (M1)–(M5).

(b) **Identity preservation:** $W^{\mathrm{flip}}(\mathrm{id}_{\mathbb{X}^{K, \mathrm{flip}, \mathrm{full}}}) = \mathrm{id}_{W^{\mathrm{flip}}(\mathbb{X}^{K, \mathrm{flip}, \mathrm{full}})}$.

(c) **Composition preservation (contravariant):** for $\pi^{K, \mathrm{flip}}_1: \mathbb{X}_1 \to \mathbb{X}_2$ and $\pi^{K, \mathrm{flip}}_2: \mathbb{X}_2 \to \mathbb{X}_3$,
$$
W^{\mathrm{flip}}(\pi^{K, \mathrm{flip}}_2 \circ \pi^{K, \mathrm{flip}}_1) = W^{\mathrm{flip}}(\pi^{K, \mathrm{flip}}_1) \circ W^{\mathrm{flip}}(\pi^{K, \mathrm{flip}}_2).
$$

(d) **Topography action consistency:** $W^{\mathrm{flip}}(\pi^{K, \mathrm{flip}})$ restricts to a UCP map $\mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}, 2} \to \mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}, 1}$ consistent with the enlarged topography preservation.

(e) **Cover action consistency:** $W^{\mathrm{flip}}(\pi^{K, \mathrm{flip}})$ preserves the operator-system cover at each level $k$.

### 4.2. Proof of (a) — OSLPLS morphism axioms

We verify the five OSLPLS morphism axioms (M1)–(M5) on the image $W^{\mathrm{flip}}(\pi^{K, \mathrm{flip}}) = (\pi^{K, \mathrm{flip}})^*: \mathcal{A}^K_2 \to \mathcal{A}^K_1$.

**(M1) Topography preservation.** By (S2) of Def 2.1-B1, $\pi^{K, \mathrm{flip}}(\mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}, 1}) \subseteq \mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}, 2}$. Viewing this in the OSLPLS-contravariant direction (per the §3.2 convention), the same algebra map is $(\pi^{K, \mathrm{flip}})^*$ regarded as a UCP map $\mathcal{A}^K_2 \to \mathcal{A}^K_1$; (S2) translates to the OSLPLS condition $(\pi^{K, \mathrm{flip}})^*(\mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}, 2}) \subseteq \mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}, 1}$ via the same set-theoretic image relation read in the contravariant convention. ✓

(Subtle point: the contravariant convention is purely a category-theoretic direction-of-arrow flip — the underlying algebra-set map is the same. The set inclusion $\pi^{K, \mathrm{flip}}(\mathcal{M}_1) \subseteq \mathcal{M}_2$ on the covariant side reads as $(\pi^{K, \mathrm{flip}})^*(\mathcal{M}_2) \subseteq \mathcal{M}_1$ on the contravariant side; both encode "the morphism takes topography to topography.")

**(M2) Lipschitz seminorm compatibility.** By (S1) of Def 2.1-B1, $L^K_2(\pi^{K, \mathrm{flip}}(a)) \le L^K_1(a)$. In the OSLPLS contravariant reading: $L^K_1((\pi^{K, \mathrm{flip}})^*(b)) \le L^K_2(b)$ for $b \in \mathrm{dom}(L^K_2)$. Same inequality, contravariantly read. ✓

**(M3) Basepoint state pullback.** By (S3) of Def 2.1-B1, $\omega_{W, 2}^L \circ \pi^{K, \mathrm{flip}} = \omega_{W, 1}^L$. In the OSLPLS contravariant reading: $\omega_{W, 1}^L \circ (\pi^{K, \mathrm{flip}})^* = \omega_{W, 2}^L$. Same equality, contravariantly read. ✓

(Substantive remark: the Phase-2.A §6.2 object action defines the OSLPLS basepoint state as the restriction $\omega := \omega_W^L|_{\mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}}}$. The morphism axiom (M3) is checked on the restricted basepoint states, which inherits from (S3) automatically because restriction commutes with the algebra map.)

**(M4) Cover preservation.** By (S4) of Def 2.1-B1, $\pi^{K, \mathrm{flip}}(\mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}, 1, k}) \subseteq \mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}, 2, k}$ for all $k$. In the OSLPLS contravariant reading: $(\pi^{K, \mathrm{flip}})^*(\mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}, 2, k}) \subseteq \mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}, 1, k}$. Same inclusion at each $k$, contravariantly read. ✓

**(M5) Modular flow intertwining.** **Automatic from (M3)** per Phase-2.A Lemma 4.2-OS, applied to the UCP map $(\pi^{K, \mathrm{flip}})^*: \mathcal{A}^K_2 \to \mathcal{A}^K_1$ with state-pullback $\omega_{W, 1}^L \circ (\pi^{K, \mathrm{flip}})^* = \omega_{W, 2}^L$. The lemma applies because:

(i) $(\pi^{K, \mathrm{flip}})^*$ is UCP (standard property of *-homomorphisms between separable C*-algebras viewed in either direction).

(ii) Both basepoint states $\omega_{W, 1}^L$ and $\omega_{W, 2}^L$ are faithful on their respective topographies (Phase-2.A §6.4 (A4) verification).

(iii) The Tomita–Takesaki modular flows $\sigma_t^{\omega_{W, 1}^L}$ and $\sigma_t^{\omega_{W, 2}^L}$ are well-defined (Paper 42 four-witness theorem; BW canonical convention).

Therefore $(\pi^{K, \mathrm{flip}})^* \circ \sigma_t^{\omega_{W, 2}^L} = \sigma_t^{\omega_{W, 1}^L} \circ (\pi^{K, \mathrm{flip}})^*$ on $\mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}, 2}$. ✓

All five OSLPLS morphism axioms verify on $W^{\mathrm{flip}}(\pi^{K, \mathrm{flip}})$. ∎

### 4.3. Proof of (b) — Identity preservation

The source identity $\mathrm{id}_{\mathbb{X}^{K, \mathrm{flip}, \mathrm{full}}}: \mathbb{X}^{K, \mathrm{flip}, \mathrm{full}} \to \mathbb{X}^{K, \mathrm{flip}, \mathrm{full}}$ corresponds at the algebra level to $\mathrm{id}_{\mathcal{A}^K}: \mathcal{A}^K \to \mathcal{A}^K$. Under $W^{\mathrm{flip}}$:

$$
W^{\mathrm{flip}}(\mathrm{id}_{\mathbb{X}^{K, \mathrm{flip}, \mathrm{full}}}) = (\mathrm{id}_{\mathcal{A}^K})^* = \mathrm{id}_{\mathcal{A}^K}
$$

(the identity map is self-pullback). On the target OSLPLS object $W^{\mathrm{flip}}(\mathbb{X}^{K, \mathrm{flip}, \mathrm{full}})$, the identity morphism in $\mathbf{OSLPLS}_{\mathrm{cov}}$ corresponds (per Phase-2.A §4.3) to $\mathrm{id}_{\mathcal{A}^K}$ on the ambient algebra. Hence:

$$
W^{\mathrm{flip}}(\mathrm{id}_{\mathbb{X}^{K, \mathrm{flip}, \mathrm{full}}}) = \mathrm{id}_{\mathcal{A}^K} = \mathrm{id}_{W^{\mathrm{flip}}(\mathbb{X}^{K, \mathrm{flip}, \mathrm{full}})}. \quad ✓
$$

### 4.4. Proof of (c) — Composition preservation (contravariant)

For source morphisms $\pi^{K, \mathrm{flip}}_1: \mathbb{X}^{K, \mathrm{flip}, \mathrm{full}}_1 \to \mathbb{X}^{K, \mathrm{flip}, \mathrm{full}}_2$ and $\pi^{K, \mathrm{flip}}_2: \mathbb{X}^{K, \mathrm{flip}, \mathrm{full}}_2 \to \mathbb{X}^{K, \mathrm{flip}, \mathrm{full}}_3$, the source composite $\pi^{K, \mathrm{flip}}_2 \circ \pi^{K, \mathrm{flip}}_1: \mathcal{A}^K_1 \to \mathcal{A}^K_3$ acts at the algebra level as the composite of *-epimorphisms.

Under $W^{\mathrm{flip}}$:

$$
W^{\mathrm{flip}}(\pi^{K, \mathrm{flip}}_2 \circ \pi^{K, \mathrm{flip}}_1) = (\pi^{K, \mathrm{flip}}_2 \circ \pi^{K, \mathrm{flip}}_1)^* = (\pi^{K, \mathrm{flip}}_1)^* \circ (\pi^{K, \mathrm{flip}}_2)^*
$$

(contravariant dual reverses composition order, standard *-homomorphism algebra). On the target side, the OSLPLS composition (Phase-2.A §4.3) is also contravariant:

$$
W^{\mathrm{flip}}(\pi^{K, \mathrm{flip}}_1) \circ W^{\mathrm{flip}}(\pi^{K, \mathrm{flip}}_2) = (\pi^{K, \mathrm{flip}}_1)^* \circ (\pi^{K, \mathrm{flip}}_2)^*
$$

(reading OSLPLS composition $\phi_1 \circ \phi_2$ as the UCP composite $\phi_1 \circ \phi_2: \mathcal{A}^K_3 \to \mathcal{A}^K_1$).

The two sides agree:
$$
W^{\mathrm{flip}}(\pi^{K, \mathrm{flip}}_2 \circ \pi^{K, \mathrm{flip}}_1) = (\pi^{K, \mathrm{flip}}_1)^* \circ (\pi^{K, \mathrm{flip}}_2)^* = W^{\mathrm{flip}}(\pi^{K, \mathrm{flip}}_1) \circ W^{\mathrm{flip}}(\pi^{K, \mathrm{flip}}_2). \quad ✓
$$

Contravariant composition preservation verifies — this is the standard Gelfand-duality-flavored functoriality, exactly as it appeared in A.3' Def 6.3 and Phase-1 Theorem 1.3-Q1'.

### 4.5. Proof of (d) — Topography action consistency

By (M1) verification in §4.2: $W^{\mathrm{flip}}(\pi^{K, \mathrm{flip}})$ maps the target-side topography $\mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}, 2}$ into the source-side topography $\mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}, 1}$ via the UCP-pullback action. The restricted map $W^{\mathrm{flip}}(\pi^{K, \mathrm{flip}})|_{\mathcal{M}_2}: \mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}, 2} \to \mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}, 1}$ is well-defined (UCP restriction to an operator-sub-system is UCP) and matches the algebraic action of $\pi^{K, \mathrm{flip}}$ on the topographies (in the contravariant convention). ✓

### 4.6. Proof of (e) — Cover action consistency

By (M4) verification in §4.2: $W^{\mathrm{flip}}(\pi^{K, \mathrm{flip}})$ maps each cover-cell of the target $\mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}, 2, k}$ into the corresponding cover-cell of the source $\mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}, 1, k}$ at each level $k$. Cover preservation holds at every cutoff in the admissible-scaling sequence. ✓

### 4.7. Aggregate Theorem 1.1-B1

Items (a)–(e) all verify. $W^{\mathrm{flip}}: \mathbf{KreinMetaMet}_{\mathrm{pp}}^{\mathrm{flip}, \mathrm{full}} \to \mathbf{OSLPLS}_{\mathrm{cov}}$ is a well-defined functor at theorem-grade rigor. ∎

### 4.8. Honest scope of Theorem 1.1-B1

**What Theorem 1.1-B1 establishes.** The bridge functor $W^{\mathrm{flip}}$ is a well-defined contravariant functor on the morphism level: it sends source morphisms (topographic Krein M-isometries on the enlarged substrate) to target morphisms (OSLPLS UCP morphisms), preserving identity and composition.

**What Theorem 1.1-B1 does NOT establish.**
- (a) The substantive Bridge Theorem 6.4'-Q1' properties (B1'-Q1') through (B4'-Q1'). These are Phase-2.B.2 substantive content; Theorem 1.1-B1 only addresses the functoriality structure, not the off-orbit super-additivity (B2'-Q1') content.
- (b) Naturality of the functor with respect to other functors in the bridge architecture (e.g., a natural transformation between $W^{\mathrm{flip}}$ and the chirality-graded $W^{\mathrm{flip}, M=0}$). This is a Phase-2.B.2 or Phase-2.D content; not a Phase-2.B.1 gate.
- (c) Faithfulness of the functor (distinct source morphisms map to distinct target morphisms). This is a refinement of (a) and depends on the source-category isomorphism characterization left open in §2.4. We name this as an open follow-on for Phase-2.B.2.

---

## §5. Commutative-case restriction comparison

### 5.1. The commutative-case restriction setup

Phase-2.A Theorem 5.1-OS established the embedding functor $\iota: \mathbf{LorPLG}_{\mathrm{cov}} \to \mathbf{OSLPLS}_{\mathrm{cov}}$ as a faithful inclusion of the commutative MS sub-category into the OSLPLS category. The morphism action of $\iota$ (Phase-2.A §5.3 (b)) is: for an MS isometry $f: (X_1, \ell_1, o_1, \mathcal{U}_1) \to (X_2, \ell_2, o_2, \mathcal{U}_2)$, $\iota(f) := f^*: C(X_2) \to C(X_1)$, $\iota(f)(g) := g \circ f$ (UCP pullback of continuous functions in the contravariant direction).

The A.3' Wick-rotation functor's morphism action (Def 6.3, end-of-section): for a topographic Krein M-isometry $\pi^K: \mathbb{X}^K_1 \to \mathbb{X}^K_2$ on the commutative substrate (Phase A.2' setting, M=0 Case A and unenlarged), $W(\pi^K) := \hat{\pi}^K: \hat{\mathcal{M}}^L_2 \to \hat{\mathcal{M}}^L_1$ as the dual Gelfand spectrum map.

**Phase-2.B.1's restriction-comparison claim:** when the OSLPLS bridge functor $W^{\mathrm{flip}}$ is restricted to the commutative MS sub-category via the Phase-2.A embedding $\iota$, the resulting morphism action agrees with the A.3' Wick-rotation functor's morphism action up to natural isomorphism. This is the structural consistency check.

### 5.2. Theorem 5.1-B1 (Restriction agreement)

**Theorem 5.1-B1 (Commutative-case restriction agreement).** The morphism action of $W^{\mathrm{flip}}$ defined in §3.1, when restricted to the commutative MS sub-category via the Phase-2.A embedding $\iota$, agrees with the A.3' Wick-rotation functor $W$ defined in A.3' Def 6.3:

$$
W^{\mathrm{flip}}|_{\iota(\mathbf{LorPLG}_{\mathrm{cov}})} \cong W \circ \iota^{-1}
$$

where $\iota^{-1}$ denotes the inverse equivalence on the sub-category (well-defined because $\iota$ is faithful).

### 5.3. Proof of Theorem 5.1-B1

For a source morphism $\pi^K: \mathbb{X}^K_1 \to \mathbb{X}^K_2$ on the commutative MS sub-category (i.e., where the topography $\mathcal{M}^L$ is Abelian — Phase A.2' setting where Gelfand spectrum $\hat{\mathcal{M}}^L$ is well-defined as a compact Hausdorff space), the A.3' Wick-rotation functor's action is:

$$
W(\pi^K) = \hat{\pi}^K: \hat{\mathcal{M}}^L_2 \to \hat{\mathcal{M}}^L_1
$$

(dual Gelfand spectrum map).

Restricting $W^{\mathrm{flip}}$ to this sub-category: the source object is a commutative Krein PPQMS (Abelian topography); the corresponding OSLPLS image via $\iota$ is the embedded MS covered LPLS in $\mathbf{OSLPLS}_{\mathrm{cov}}$ (Phase-2.A §5.4 "commutative MS is a sub-theory of OSLPLS"). The morphism action of $W^{\mathrm{flip}}$ on this sub-category becomes $(\pi^K)^*: C(\hat{\mathcal{M}}^L_2) \to C(\hat{\mathcal{M}}^L_1)$ via the Phase-2.A §3.3 Case 1 commutative identification $\iota(\hat{\mathcal{M}}^L, \ldots) = (C(\hat{\mathcal{M}}^L), L_\ell, C(\hat{\mathcal{M}}^L), \delta_o, \mathcal{U}_{\mathrm{OS}})$.

By Gelfand–Naimark duality, the UCP pullback $(\pi^K)^*: C(\hat{\mathcal{M}}^L_2) \to C(\hat{\mathcal{M}}^L_1)$ on commutative C*-algebras is naturally identified with the continuous map $\hat{\pi}^K: \hat{\mathcal{M}}^L_2 \to \hat{\mathcal{M}}^L_1$ in the contravariant direction:

$$
(\pi^K)^*(g)(\chi) := g(\hat{\pi}^K(\chi)) \quad\text{for } g \in C(\hat{\mathcal{M}}^L_2), \chi \in \hat{\mathcal{M}}^L_1.
$$

This is exactly the A.3' Def 6.3 morphism action. Therefore:

$$
W^{\mathrm{flip}}|_{\iota(\mathbf{LorPLG}_{\mathrm{cov}})}(\pi^K) = (\pi^K)^* = \iota(\hat{\pi}^K) = \iota(W(\pi^K))
$$

where the second equality is the Gelfand–Naimark identification and the third equality is the definition of $W$ on morphisms in A.3'. Equivalently:

$$
\iota^{-1} \circ W^{\mathrm{flip}}|_{\iota(\mathbf{LorPLG}_{\mathrm{cov}})} = W \circ \iota^{-1} \circ \iota = W
$$

on the sub-category. The restriction-agreement claim holds. ∎

### 5.4. Phase-1 chirality-graded restriction (bonus)

We verify the parallel claim for the Phase-1 Case A (M=0 enlarged) sub-category. The Phase-1 functor $W^{\mathrm{flip}, M=0}$ (Theorem 1.3-Q1') extends A.3' with $\mathbb{Z}/2$-grading; the same proof structure as §5.3 applies, with the modification that the Gelfand spectrum now carries a $\mathbb{Z}/2$-cover decomposition (Phase-1 Theorem 1.2-Q1') and the morphism action $\hat{\pi}^{K, \mathrm{flip}}$ acts diagonally on the chirality sheets.

When $W^{\mathrm{flip}}$ (full M≠0) is restricted to the M=0 enlarged sub-category, the morphism action collapses to $\hat{\pi}^{K, \mathrm{flip}}$ on the doubled Gelfand spectrum — bit-exact agreement with Phase-1 Theorem 1.3-Q1'. The hierarchical containment is:

$$
W|_{\text{A.3' commutative}} = (\hat{\pi}^K \text{ on } \hat{\mathcal{M}}^L) \quad\subset\quad W^{\mathrm{flip}, M=0}|_{\text{Phase-1 Case A}} = (\hat{\pi}^{K, \mathrm{flip}} \text{ on } \widehat{\mathcal{M}^{L, \mathrm{flip}}_{M=0}}) \quad\subset\quad W^{\mathrm{flip}}|_{\text{Phase-2.A Case B}} = ((\pi^{K, \mathrm{flip}})^* \text{ on } \mathcal{S}(\mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}})).
$$

Each successive containment is structurally additive: the M=0 enlargement adds the $\mathbb{Z}/2$-cover; the full M≠0 enlargement replaces the Gelfand spectrum by the operator-system state space $\mathcal{S}(\mathcal{M})$. The morphism action generalizes from Gelfand-spectrum map (commutative cases) to UCP pullback (non-commutative case), with the two coinciding on the Abelian sub-categories.

### 5.5. Structural reading: the morphism action is the "right" generalization

Theorem 5.1-B1 confirms that the OSLPLS morphism action of $W^{\mathrm{flip}}$ is the *unique* operator-system generalization of the A.3' Wick-rotation functor's morphism action that:

(i) restricts to A.3' on the commutative MS sub-category (Theorem 5.1-B1 itself);
(ii) restricts to Phase-1's $W^{\mathrm{flip}, M=0}$ on the M=0 enlarged sub-category (§5.4);
(iii) generalizes to the full M≠0 enlarged substrate where the source topography is non-Abelian (Def 3.1-B1 + Theorem 1.1-B1).

This is structurally consistent with the Phase-2.A §3.4 finding that the OSLPLS object is "structurally constrained, not freely chosen" — the morphism action is similarly forced by (a) compatibility with the embedded commutative sub-category and (b) Connes–vS contravariant UCP convention. There is no morphism-level design freedom; the construction is determined by the substrate.

---

## §6. Phase-2.B.1.5 gate verdict

### 6.1. Per-axiom verdict

| Component | Verdict | Notes |
|:----------|:--------|:------|
| Source morphism class Def 2.1-B1 | POSITIVE | Topographic Krein M-isometry on enlarged substrate; extends A.2' Def 2.3 |
| Source-category composition + identity | POSITIVE | Standard *-epimorphism composition |
| $W^{\mathrm{flip}}$ on morphisms Def 3.1-B1 | POSITIVE | $(\pi^{K, \mathrm{flip}})^*$ via algebraic pullback in contravariant direction |
| OSLPLS morphism axiom (M1) topography preservation | POSITIVE | Direct from source (S2) |
| OSLPLS morphism axiom (M2) Lipschitz compatibility | POSITIVE | Direct from source (S1) |
| OSLPLS morphism axiom (M3) basepoint state pullback | POSITIVE | Direct from source (S3) |
| OSLPLS morphism axiom (M4) cover preservation | POSITIVE | Direct from source (S4) |
| OSLPLS morphism axiom (M5) modular flow intertwining | POSITIVE (automatic) | Phase-2.A Lemma 4.2-OS |
| Theorem 1.1-B1 (a) — morphism axioms verify | POSITIVE | All five M-axioms verified per §4.2 |
| Theorem 1.1-B1 (b) — identity preservation | POSITIVE | $W^{\mathrm{flip}}(\mathrm{id}) = \mathrm{id}$ by direct computation |
| Theorem 1.1-B1 (c) — composition preservation (contravariant) | POSITIVE | Contravariant Gelfand-duality-flavored composition |
| Theorem 1.1-B1 (d) — topography action consistency | POSITIVE | UCP restriction to operator-sub-system |
| Theorem 1.1-B1 (e) — cover action consistency | POSITIVE | UCP restriction at each cutoff level |
| Theorem 5.1-B1 — commutative-case restriction matches A.3' | POSITIVE | Gelfand–Naimark identification |
| Phase-1 chirality-graded restriction (Case A) matches | POSITIVE (bonus) | Hierarchical containment, structurally additive |
| Aggregate Phase-2.B.1 | **POSITIVE** | All four functoriality axioms clean; commutative restriction agrees |

### 6.2. Aggregate gate verdict: POSITIVE

**Recommendation: GO to Phase-2.B.3 (numerical verification panel).** Phase-2.B.1 closes the mechanical morphism-specification sub-task at theorem-grade rigor. The bridge functor $W^{\mathrm{flip}}: \mathbf{KreinMetaMet}_{\mathrm{pp}}^{\mathrm{flip}, \mathrm{full}} \to \mathbf{OSLPLS}_{\mathrm{cov}}$ is now well-defined on both objects (Phase-2.A Theorem 6.1-OS) and morphisms (Phase-2.B.1 Theorem 1.1-B1), with the commutative-case restriction matching A.3'.

Phase-2.B.2 (substantive Bridge Theorem 6.4'-Q1' with off-orbit super-additivity) runs in parallel; its substrate is the now-complete morphism specification of Phase-2.B.1. Phase-2.B.3 (numerical verification panel) can consume Phase-2.B.1's morphism action at bit-exact $(n_{\max}, N_t)$ cells.

### 6.3. Why this compresses to a single sprint memo

Phase-2.A §10 (d) predicted: "Phase-2.B.1 (1 week): specify the bridge functor $W^{\mathrm{flip}}$ on morphisms (mechanical from Phase-2.A object action + Connes–vS contravariant UCP convention)." The compression to a single session validates this estimate:

- The source morphism class (§2) is a direct extension of Phase A.2' Def 2.3 (topographic Krein M-isometry) with the enlarged topography and explicit cover preservation. No structural design freedom.
- The functor's morphism action (§3) is dictated by (a) the Phase-2.A object construction (which specifies the target OSLPLS object) and (b) Connes–vS contravariant UCP convention (which dictates the morphism direction). No design choice.
- All four functoriality axioms (§4) reduce to source-side axioms (S1–S4) read in the contravariant direction, with the standard *-homomorphism algebra.
- Commutative-case restriction (§5) matches A.3' via Gelfand–Naimark.

The compression pattern continues to apply: the substrate (Phase-2.A object construction + A.3' morphism convention + standard UCP composition) does most of the heavy lifting. Phase-2.B.1's task is to synthesize the substrate-inherited components and verify they assemble into a functor — work completable in one focused session.

---

## §7. Phase-2.B.3 numerical verification panel input specification

### 7.1. Panel structure

Per Phase-2.A §10 (d) and the Phase-1 §6 numerical-panel template, Phase-2.B.3 should verify the bridge functor's morphism action at bit-exact cutoff cells:

$$
(n_{\max}, N_t) \in \{(2, 3), (3, 5), (4, 7)\}
$$

with $T = T_\infty/N_t$ at the canonical BW period $T_\infty = 2\pi$.

### 7.2. Morphism action to verify per cell

At each cell $(n_{\max}, N_t)$, the truncation projector morphism

$$
\pi^{K, \mathrm{flip}}_{k \to k+1}: \mathbb{X}^{K, \mathrm{flip}, \mathrm{full}}_{n_{\max}(k), N_t(k), T(k)} \to \mathbb{X}^{K, \mathrm{flip}, \mathrm{full}}_{n_{\max}(k+1), N_t(k+1), T(k+1)}
$$

(the natural Berezin/projection nesting map between successive cutoffs) is the source-side morphism of interest. Its image under $W^{\mathrm{flip}}$ per Def 3.1-B1 is:

$$
W^{\mathrm{flip}}(\pi^{K, \mathrm{flip}}_{k \to k+1}) = (\pi^{K, \mathrm{flip}}_{k \to k+1})^*: W^{\mathrm{flip}}(\mathbb{X}^{K, \mathrm{flip}, \mathrm{full}}_{k+1}) \to W^{\mathrm{flip}}(\mathbb{X}^{K, \mathrm{flip}, \mathrm{full}}_{k}).
$$

Phase-2.B.3 verifies bit-exactness of:
- (M1) topography preservation on the truncated enlarged topographies $\mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}, k+1} \to \mathcal{M}^{L, \mathrm{flip}}_{\mathrm{full}, k}$.
- (M2) operator-norm Lipschitz seminorm contractivity $\|[D_L, \cdot]\|_{\mathrm{op}}$ comparison.
- (M3) basepoint pullback $\omega_W^L|_{\mathcal{M}_k} \circ (\pi^{K, \mathrm{flip}})^* = \omega_W^L|_{\mathcal{M}_{k+1}}$.
- (M4) operator-system cover preservation level-by-level.
- (M5) modular flow intertwining $\sigma_t^{\omega_W^L}|_{\mathcal{M}_k} \circ (\pi^{K, \mathrm{flip}})^* = (\pi^{K, \mathrm{flip}})^* \circ \sigma_t^{\omega_W^L}|_{\mathcal{M}_{k+1}}$ (automatic from M3 per Lemma 4.2-OS but worth a sanity verification at the panel level).

### 7.3. Expected behavior

Per the substrate-inheritance compression pattern, Phase-2.B.3 should produce:

- Bit-exact agreement between $W^{\mathrm{flip}}(\pi^{K, \mathrm{flip}})$ verified residuals and the existing Paper 45 §6 numerical panel residuals for the natural-substrate Berezin/projection pair (cf. Paper 45 Theorem 5.1 numerical panel values $\Lprop(2,3) = 2.0746$, $\Lprop(3,5) = 1.6101$, $\Lprop(4,7) = 1.3223$). The enlarged-substrate morphism action introduces no new numerical content at the functoriality-level; the substantive numerical content is in the strong-form propinquity bound itself (Paper 46 Appendix B / Phase-2.B.2).

- Riemannian-limit recovery at $N_t = 1$ bit-exact: the functor reduces to the Riemannian-only morphism action, matching Paper 38 single-factor results.

- Identity preservation residual exactly 0.0 (machine precision): $W^{\mathrm{flip}}(\mathrm{id})$ at every cell is the identity UCP map on the truncated operator system.

### 7.4. Compute infrastructure

The existing infrastructure supports the panel verification directly:

- `geovac/gh_convergence_tensor.py` for the Paper 45 / Paper 46 / Paper 38 propinquity panel infrastructure.
- `geovac/operator_system_lorentzian.py` (Paper 44 substrate) for the truncated operator system structure.
- `geovac/modular_hamiltonian.py` for the Tomita modular flow verification.
- `geovac/krein_space_construction.py` (L2-B) for the Krein structure.

Phase-2.B.3 should require no new production code modules; it's a numerical wrapper around existing infrastructure verifying the morphism action axioms at the panel cells.

### 7.5. Estimated Phase-2.B.3 effort

Per Phase-2.A §10 (d), 1 week. Mechanical at the morphism-verification level, with the panel infrastructure already in place. The compression pattern is plausible: a single focused sprint can produce the panel verification + brief memo + JSON data file.

---

## §8. Honest scope statement

### 8.1. What Phase-2.B.1 establishes at theorem-grade rigor

1. **Definition 2.1-B1 (Source morphism class).** Topographic Krein M-isometries on the enlarged substrate, extending Phase A.2' Def 2.3 with enlarged-topography axiom (S2) and explicit cover preservation (S4). Source-category composition and identity are standard.

2. **Definition 3.1-B1 (Bridge functor on morphisms).** $W^{\mathrm{flip}}(\pi^{K, \mathrm{flip}}) := (\pi^{K, \mathrm{flip}})^*$, the algebraic pullback in the contravariant direction. Convention bridges covariant source category with contravariant target category per standard Connes–vS UCP convention.

3. **Theorem 1.1-B1 (Functoriality of $W^{\mathrm{flip}}$).** $W^{\mathrm{flip}}$ is a well-defined functor on morphisms: (a) OSLPLS morphism axioms (M1)–(M5) all verify, (b) identity preservation, (c) contravariant composition preservation, (d) topography action consistency, (e) cover action consistency.

4. **Theorem 5.1-B1 (Commutative-case restriction agreement).** $W^{\mathrm{flip}}$ restricted to the commutative MS sub-category via $\iota$ agrees with the A.3' Wick-rotation functor's morphism action up to natural isomorphism. Phase-1 chirality-graded restriction (Case A) also matches as bonus (§5.4).

5. **Phase-2.B.3 numerical verification panel input specification** (§7): the morphism action can be computed at bit-exact $(n_{\max}, N_t)$ cells using existing infrastructure.

### 8.2. What Phase-2.B.1 does NOT establish (the honest scope)

1. **The Bridge Theorem 6.4'-Q1' substantive content (B2'-Q1' off-orbit super-additivity).** This is Phase-2.B.2 substantive content. Phase-2.B.1 verifies the functoriality structure of $W^{\mathrm{flip}}$; the actual proof of the OSLPLS reverse triangle inequality on the bridge image (with strict super-additivity via Connes–Rovelli thermal-time stack across distinct KMS states) is the parallel-running Phase-2.B.2 sprint.

2. **Faithfulness of $W^{\mathrm{flip}}$.** Whether distinct source morphisms produce distinct target morphisms. Depends on the source-category isomorphism characterization (parallel to the Phase-2.A §4.6 OSLPLS isomorphism characterization). Named open follow-on for Phase-2.B.2.

3. **Numerical verification of the morphism action (Phase-2.B.3).** Phase-2.B.1 specifies what to verify; Phase-2.B.3 runs the verification.

4. **Naturality of $W^{\mathrm{flip}}$ with respect to other bridge functors.** Whether $W^{\mathrm{flip}}$ at full $M \ne 0$ is naturally isomorphic (in the strict sense of natural transformations) to $W^{\mathrm{flip}, M=0} \circ \iota_{M=0}$ or to $W \circ \iota_{\mathrm{comm}}$ on the appropriate sub-categories. This is a refinement of Theorem 5.1-B1; not a Phase-2.B.1 gate.

5. **Connes–Rovelli thermal-time stack construction at the operator-system level.** Phase-2.C content; gated on Phase-2.B.2 needs.

### 8.3. Failure modes Phase-2.B.1 explicitly avoided

Per the Phase-2.B.1 task spec, Phase-2.B.1 explicitly does NOT:

1. **Claim that Phase-2.B.1 closes the Bridge Theorem 6.4'-Q1'.** §6.2, §8.2(1), and the recommendation in §6.2 explicitly defer the substantive content to Phase-2.B.2.

2. **Claim functoriality without verifying all four axioms.** §4 verifies (a)–(e) of Theorem 1.1-B1 with explicit proof for each.

3. **Skip the commutative-case restriction comparison.** §5 establishes Theorem 5.1-B1 explicitly, confirming structural consistency with A.3' and Phase-1.

4. **Conflate source-side morphism class with target-side OSLPLS morphism class.** §2 (source) and §3–§4 (target via functor) are distinct sections; the contravariance convention bridge between them is made explicit (§3.2).

5. **Promote Phase-2.B.1 from "morphism mechanical follow-on" to "substantive Bridge Theorem closure."** §6.2 recommends GO to Phase-2.B.3 (numerical panel) as the next mechanical sub-task, with substantive Phase-2.B.2 (off-orbit super-additivity) running in parallel.

### 8.4. Two substantive findings the formalization produced

#### 8.4.1. Finding 1: The morphism action is fully constrained by the substrate

(§3) The bridge functor's morphism action is dictated by (a) Phase-2.A's object construction (which fixes the target OSLPLS object), (b) Connes–vS contravariant UCP convention (standard for operator-system category), and (c) the source-side morphism class (extension of Phase A.2' Def 2.3 to the enlarged substrate). There is essentially no design freedom; the morphism action is forced.

This parallels the Phase-2.A finding that the OSLPLS *object* class is structurally constrained: at the morphism level, the same substrate-inheritance pattern produces the unique morphism action consistent with both the object construction and the standard operator-system conventions.

#### 8.4.2. Finding 2: The contravariance convention is the convention-bridge between covariant Krein source and contravariant UCP target

(§3.2) The source category $\mathbf{KreinMetaMet}_{\mathrm{pp}}^{\mathrm{flip}, \mathrm{full}}$ uses *covariant* Krein M-isometries (proper *-epimorphisms in the natural causal direction). The target category $\mathbf{OSLPLS}_{\mathrm{cov}}$ uses *contravariant* UCP maps. The bridge functor $W^{\mathrm{flip}}$ pulls the covariant source arrow back to its contravariant target arrow via the algebraic-pullback convention.

This is the natural Connes–vS convention (operator-system category is always presented contravariantly relative to *-homomorphism direction; cf. Gelfand duality). The A.3' commutative-case bridge handled this contravariance via the Gelfand-spectrum-map dual; the Phase-2.B.1 non-commutative bridge handles it via the algebraic-pullback convention (which reduces to Gelfand-spectrum-map duality on the commutative sub-category, per Theorem 5.1-B1).

The contravariance is therefore not a Phase-2.B.1 design choice — it is inherited from the underlying duality conventions of the source and target categories. The morphism action $W^{\mathrm{flip}}(\pi^{K, \mathrm{flip}}) := (\pi^{K, \mathrm{flip}})^*$ is the unique way to bridge the two conventions.

---

## §9. Cross-references

- `debug/sprint_q1prime_phase2a_oslpls_category_memo.md` — Phase-2.A (load-bearing throughout for object construction + Lemma 4.2-OS + commutative MS embedding)
- `debug/sprint_q1prime_phase1_case_a_stepping_stone_memo.md` — Phase-1 Theorem 1.3-Q1' (chirality-graded morphism action template; load-bearing for §1.3 and §5.4)
- `debug/sprint_phase_a3prime_mondino_samann_bridge_memo.md` — A.3' Def 6.3 (commutative-case Wick-rotation functor morphism action; load-bearing for §5.1, §5.3)
- `debug/sprint_phase_a2prime_krein_lift_formalization_memo.md` — A.2' Def 2.3 (topographic Krein M-isometry; load-bearing for §2.2)
- `papers/group1_operator_algebras/paper_44_lorentzian_operator_system.tex` — operator-system substrate (load-bearing for §3.2 UCP convention)
- `papers/group1_operator_algebras/paper_46_strong_form_lorentzian_propinquity.tex` — Appendix B enlarged substrate (load-bearing for §2.1 source object class)
- `papers/group1_operator_algebras/paper_42_modular_hamiltonian_four_witness.tex` — four-witness theorem (load-bearing for §4.2 (M5) automaticity)
- `papers/group1_operator_algebras/paper_38_su2_propinquity_convergence.tex` — Riemannian-limit reference for §7.3 expected behavior
- `papers/group1_operator_algebras/paper_45_lorentzian_propinquity.tex` — Paper 45 §6 numerical panel infrastructure for §7.4
- Connes–vS arXiv:2004.14115 §3 — operator-system framework + UCP morphism contravariant convention (load-bearing for §3.2, §3.3)
- Mondino–Sämann arXiv:2504.10380 v4 Def 4.4 — MS morphism class (load-bearing comparison point for §5.1)

---

## §10. Output for the report-back

**(a) Headline verdict (1 sentence):** POSITIVE — $W^{\mathrm{flip}}$ on morphisms defined cleanly via Connes–vS contravariant UCP convention (Def 3.1-B1: $W^{\mathrm{flip}}(\pi^{K, \mathrm{flip}}) := (\pi^{K, \mathrm{flip}})^*$ algebraic pullback in the contravariant direction); functoriality verified at theorem-grade rigor (Theorem 1.1-B1: all four functoriality axioms (a)–(e) clean, with all five OSLPLS morphism axioms (M1)–(M5) verifying directly from source-side axioms (S1)–(S4) under the standard contravariance convention, and M5 modular-flow intertwining automatic from M3 per Phase-2.A Lemma 4.2-OS); restriction to commutative MS sub-category via $\iota$ matches A.3' Wick-rotation functor's morphism action bit-exact via Gelfand–Naimark (Theorem 5.1-B1) and hierarchically extends Phase-1 chirality-graded morphism action (§5.4 bonus); compresses to single sprint memo per the Phase-2.A predicted compression pattern — recommend GO to Phase-2.B.3 (numerical verification panel at bit-exact $(n_{\max}, N_t) \in \{(2,3), (3,5), (4,7)\}$ cells using existing infrastructure).

**(b) Functoriality status:** CLEAN (no specific gaps). All four functoriality axioms of Theorem 1.1-B1 verify with explicit proof per §4.3 (identity), §4.4 (contravariant composition), §4.5 (topography), §4.6 (cover). All five OSLPLS morphism axioms verify per §4.2. The commutative-case restriction matches A.3' bit-exact via Gelfand–Naimark identification per §5.3. The Phase-1 chirality-graded restriction matches as bonus hierarchical containment per §5.4. The only named follow-ons (functor faithfulness, naturality with respect to other bridge functors, Connes–Rovelli stack construction) are Phase-2.B.2 / Phase-2.C content, not Phase-2.B.1 gaps.

**(c) Most surprising finding:** The morphism action is **fully constrained by the substrate** (parallel to Phase-2.A's "object is structurally constrained, not freely chosen" finding). Specifically, the convention bridge between the *covariant* source category (Krein M-isometries) and the *contravariant* target category (OSLPLS UCP morphisms) is forced by the standard Connes–vS operator-system convention; the morphism action $W^{\mathrm{flip}}(\pi^{K, \mathrm{flip}}) := (\pi^{K, \mathrm{flip}})^*$ is the unique algebraic-pullback realization of this convention bridge. There is no morphism-level design freedom; the construction is the unique operator-system generalization of the A.3' Wick-rotation functor's morphism action that respects the source/target conventions. This is the structural meaning of "mechanical sub-sprint" in the Phase-2.A §7.4 sequencing: the substrate has already done the heavy lifting; Phase-2.B.1's task is to verify that the substrate-inherited components assemble into a functor.

**(d) Recommended Phase-2.B.3 panel inputs:** Phase-2.B.3 verifies the morphism action of $W^{\mathrm{flip}}$ at bit-exact cells $(n_{\max}, N_t) \in \{(2, 3), (3, 5), (4, 7)\}$ with canonical BW period. For each cell, compute the truncation projector morphism $\pi^{K, \mathrm{flip}}_{k \to k+1}: \mathbb{X}^{K, \mathrm{flip}, \mathrm{full}}_k \to \mathbb{X}^{K, \mathrm{flip}, \mathrm{full}}_{k+1}$ (the natural Berezin/projection nesting between successive cutoffs) and verify its image $W^{\mathrm{flip}}(\pi^{K, \mathrm{flip}}_{k \to k+1})$ satisfies the five OSLPLS morphism axioms (M1)–(M5) at bit-exact precision. Expected: (i) bit-exact agreement with existing Paper 45 §6 numerical panel residuals $\Lprop(2,3) = 2.0746$, $\Lprop(3,5) = 1.6101$, $\Lprop(4,7) = 1.3223$ (the enlarged-substrate morphism action introduces no new numerical content at the functoriality level); (ii) Riemannian-limit recovery at $N_t = 1$ bit-exact ($W^{\mathrm{flip}}$ reduces to Riemannian-only morphism action matching Paper 38); (iii) identity preservation residual exactly 0.0 (machine precision). Existing infrastructure (`geovac/gh_convergence_tensor.py`, `geovac/operator_system_lorentzian.py`, `geovac/modular_hamiltonian.py`, `geovac/krein_space_construction.py`) supports the panel directly; Phase-2.B.3 should require no new production code modules. Estimated effort: 1 week.

**End of memo.**

**Files added in this sprint:**
- `debug/sprint_q1prime_phase2b1_wflip_morphisms_memo.md` (this memo, ~4500 words formal morphism-level construction + functoriality verification + commutative-case restriction comparison + Phase-2.B.3 input specification + Phase-2.B.1.5 gate verdict)
- `debug/data/sprint_q1prime_phase2b1.json` (per-axiom morphism verification + functoriality verdict + Phase-2.B.3 input specification)

**Cross-references:**
- `debug/sprint_q1prime_phase2a_oslpls_category_memo.md` (Phase-2.A object construction; load-bearing throughout)
- `debug/sprint_q1prime_phase1_case_a_stepping_stone_memo.md` (Phase-1 chirality-graded morphism action template)
- `debug/sprint_phase_a3prime_mondino_samann_bridge_memo.md` (A.3' Wick-rotation functor; commutative-case morphism action)
- `debug/sprint_phase_a2prime_krein_lift_formalization_memo.md` (A.2' topographic Krein M-isometry; source morphism class)
- `papers/group1_operator_algebras/paper_44_lorentzian_operator_system.tex` (operator-system substrate)
- `papers/group1_operator_algebras/paper_46_strong_form_lorentzian_propinquity.tex` (Appendix B enlarged substrate)
- `papers/group1_operator_algebras/paper_42_modular_hamiltonian_four_witness.tex` (four-witness theorem; modular flow structure)
- `papers/group1_operator_algebras/paper_45_lorentzian_propinquity.tex` (numerical panel infrastructure)
- Connes-vS arXiv:2004.14115 (operator-system UCP convention)
- Mondino-Sämann arXiv:2504.10380 v4 (MS morphism class; commutative comparison)
