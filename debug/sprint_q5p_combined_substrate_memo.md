# Sprint Q5'-Combined-Substrate — first multi-year scoping step combining T3a (Peter–Weyl $J^*(S^3)$) and T3b (cross-shell off-diagonal Dirac) on a SINGLE combined substrate

**Date:** 2026-06-06 (close-of-day follow-on to Sprint Q5'-HardParts-Round2, v3.62.0)
**Sprint scope:** the categorical question of L1 vs L5 at the substrate level — does combining T3a + T3b via TENSOR PRODUCT (L1) vs DIRECT-SUM-PLUS-RELATIONS / EMBEDDING (L5) give a genuinely richer Hopf-algebra structure?
**Drivers:** `debug/compute_q5p_combined_substrate.py` + `debug/compute_q5p_combined_substrate_part2.py`
**Data:** `debug/data/sprint_q5p_combined_substrate.json` + `debug/data/sprint_q5p_combined_substrate_part2.json`
**Wall time:** 0.008 s (Part 1) + ~0.3 s (Part 2)
**Discipline:** bit-exact `sympy.Rational` + symbolic SU(2) rotation; no PSLQ; no floats.

---

## 1. TL;DR

**Verdict: POSITIVE-REDUCES-TO-L1.** The combined substrate L5 obtained by attempting to *embed* the OffDiag (T3b) transition generators $T_{s' \to s} := e_s \cdot (\kappa A) \cdot e_{s'}$ inside the Peter–Weyl matrix-coefficient algebra (T3a) **reduces categorically to L1, the tensor product** $\mathcal{A}^{J^*}_{j_{\max}} \otimes \mathcal{A}^{\mathrm{OD}}_{n_{\max}}$, at the basic-substrate level. The structural reason is **two-sided categorical**:

(i) **Peter–Weyl matrix coefficients are FUNCTIONS** on $SU(2) \cong SL_2$ (polynomial functions $\pi^j_{mn} \in \mathcal{O}(SL_2)$); **OffDiag transitions are SCALAR-VALUED cocycle data** (specific rational $\eta$-class values in $\kappa^2 \cdot \mathbb{Z}$). No $SL_2$ group element evaluates the 4 spin-$1/2$ matrix coefficients $\{a, b, c, d\}$ to the 12 OffDiag $\eta$-class rationals at all 12 transitions simultaneously — the two are categorically different objects (polynomial functions on a group vs scalar traces on operators).

(ii) **The natural $SU(2)$ action does NOT preserve the basic 5-sector OffDiag basis.** Bit-exactly verified symbolically (`debug/compute_q5p_combined_substrate_part2.py`): the z-rotation $U_\theta = \mathrm{diag}(e^{i\theta m_J})$ conjugates a model E1 transition $T_{(1,1)\to(1,0)}$ to a matrix with three distinct phase factors $\{1, e^{i\theta}, e^{-i\theta}\}$ — NOT a scalar multiple of the original. The result lives in the span of $m_J$-RESOLVED sub-transitions, not in the span of the basic 5-sector $T_{s'\to s}$. For $SU(2)$ to act non-trivially on a smash product, the OffDiag basis would need to be REFINED to $m_J$-resolved sub-idempotents — i.e., the basic OffDiag substrate of v3.61.0 / T3b would no longer be invariant.

**Headline structural finding (joint):** The combined substrate L5 IS structurally L1 at the categorical-tensor-product level. The candidate motivic Galois group at the quotient is
$$
\boxed{
U^{*(n_{\max}, j_{\max})}_{\mathrm{combined}} \;=\; SL_2 \times \mathbb{G}_a^{N_{\mathrm{OD}}(n_{\max})}
}
$$
— the **Levi-decomposition** shape (semisimple $\times$ pro-unipotent) — exactly the multi-year synthesis target predicted by T3a's headline finding. **L1 = L5 categorically; the Levi decomposition is the natural Stage-2 target structurally.**

What is *gained* by stating this explicitly (not a tautology of T3a + T3b separately): the substrate enrichment program now has a **theorem-grade categorical closure** at finite cutoff — the multi-year Stage-2 Tannakian construction reduces to ONE substrate (the Levi-decomposition product), not THREE independent ingredients to be combined later. The structural barriers to L5 (categorical incompatibility of polynomial-on-group vs scalar-on-operator data; SU(2) non-invariance of the OffDiag basis) are NAMED, not undiscovered, and serve as falsifiers for any future attempt to upgrade the categorical product to a smash product.

---

## 2. Verdict against decision gate

| Gate | Selected? | Why |
|:-----|:---------:|:----|
| POSITIVE-NEW-STRUCTURE | not selected | The candidate "embedding" $T_{s'\to s} \mapsto c \cdot \pi^{j}_{m_s, n_{s'}}$ has no consistent realization: the OffDiag $\eta$-class rationals are scalar cocycle values; Peter–Weyl matrix coefficients are polynomial functions on a group. The two live in different categories. |
| **POSITIVE-REDUCES-TO-L1** | **selected** | Combined substrate $\mathcal{A}^{J^*} \otimes \mathcal{A}^{\mathrm{OD}}$ is the categorical tensor product (Section 5). The Peter–Weyl factor's $SU(2)$ symmetry does NOT preserve the OffDiag basis (Section 5–6, symbolic bit-exact). $U^* = SL_2 \times \mathbb{G}_a^{N_{\mathrm{OD}}}$ is the Levi-decomposition shape — the natural Stage-2 synthesis target. |
| BORDERLINE | not selected | Both walls (categorical mismatch + SU(2) non-invariance) are clean structural obstructions, not gray-zone partial interactions. |
| STOP | not selected | The categorical tensor product L1 is well-defined and gives the Levi-decomposition structure. The substrate enrichment program is not blocked — it is closed at the L1 level with explicit Stage-2 target. |

---

## 3. The categorical question

### 3.1 Two candidate combinations

(a) **TENSOR PRODUCT (L1-style):** $\mathcal{A}^{\mathrm{combined}} = \mathcal{A}^{J^*}_{j_{\max}} \otimes \mathcal{A}^{\mathrm{OD}}_{n_{\max}}$. The categorical Hopf algebra tensor product (Klimyk–Schmüdgen Prop 1.3.5): multiplication is component-wise; coproduct is $\Delta^{\mathrm{combined}} = (\Delta^{J^*} \otimes \Delta^{\mathrm{OD}}) \circ (\mathrm{id} \otimes \tau \otimes \mathrm{id})$ with the standard flip $\tau$.

(b) **EMBEDDED (genuinely new):** Embed the OffDiag transitions $T_{s' \to s}$ AS matrix coefficients (or some quotient/extension) of Peter–Weyl. Specifically: at $n_{\max} = 2$, attempt $T_{s'\to s} = c_{s', s} \cdot \pi^{j(s, s')}_{m_s, n_{s'}}$ for some rational $c_{s', s}$ and specific $(j, m, n)$ indices determined by Wigner–Eckart selection rules.

### 3.2 Test of the embedding hypothesis

At $n_{\max} = 2$, the 12 OffDiag transitions carry the following $\eta$-class values (from `sprint_q5p_offdiag_dirac_memo.md` §4.2):

| Transition | $\eta(T_{s' \to s})$ | $\eta / \kappa^2$ (integer) |
|:----------:|:-------------------:|:----:|
| $(1,0) \to (1,1)$ | $+1/64$ | $+4$ |
| $(1,0) \to (2,1)$ | $+5/128$ | $+10$ |
| $(1,1) \to (1,0)$ | $-1/64$ | $-4$ |
| $(1,1) \to (2,0)$ | $-1/64$ | $-4$ |
| $(1,1) \to (2,2)$ | $-3/128$ | $-6$ |
| $(2,0) \to (1,1)$ | $+1/64$ | $+4$ |
| $(2,0) \to (2,1)$ | $+5/128$ | $+10$ |
| $(2,1) \to (1,0)$ | $+1/128$ | $+2$ |
| $(2,1) \to (2,0)$ | $+1/128$ | $+2$ |
| $(2,1) \to (2,2)$ | $+1/64$ | $+4$ |
| $(2,2) \to (1,1)$ | $-3/128$ | $-6$ |
| $(2,2) \to (2,1)$ | $-1/16$ | $-16$ |

These are integers in $\kappa^2 \cdot \mathbb{Z}$ — chirality-weighted two-step E1 path counts.

The Peter–Weyl matrix coefficients $\pi^{1/2}_{m,n}$ at $j_{\max} = 1/2$ are FUNCTIONS on $SL_2$:
$$
\pi^{1/2}_{m,n}(g) = (m,n) \text{-entry of } \rho_{1/2}(g) \qquad g \in SL_2,
$$
where $\rho_{1/2}$ is the defining 2-dim representation. As polynomials over $\mathbb{Z}$:
- $a = \pi^{1/2}_{+1/2, +1/2}$, $b = \pi^{1/2}_{+1/2, -1/2}$, $c = \pi^{1/2}_{-1/2, +1/2}$, $d = \pi^{1/2}_{-1/2, -1/2}$
- Defining relation: $ad - bc = 1$ (det = 1)
- Evaluation at identity: $a = d = 1$, $b = c = 0$.

A consistent embedding $T_{s'\to s} \leftrightarrow \pi^j_{m, n}$ would require:

(i) **Index identification:** a function $(s, s') \leftrightarrow (j, m, n)$. The CH sectors are triangular-indexed (5 sectors at $n_{\max}=2$); the Peter–Weyl indices are square-indexed per $j$-shell (1 + 4 = 5 at $j_{\max} = 1/2$). Dimensional match is coincidence; the index structures differ.

(ii) **Value match:** a group element $g_0 \in SL_2$ such that $\pi^j_{m,n}(g_0)$ equals the OffDiag $\eta$-class value $c_{s',s} \cdot \eta(T_{s'\to s})$ at all 12 transitions simultaneously. At $j_{\max} = 1/2$, this would require ONE group element $g_0$ producing FOUR specific values for $\{a(g_0), b(g_0), c(g_0), d(g_0)\}$ — but the OffDiag data has 12 values, not 4. Even at $j_{\max} = 3/2$, the Peter–Weyl basis has 30 elements, but the OffDiag $\eta$-data is 12-dimensional, and the values are FORCED by chirality-weighted path counts that have nothing to do with $SL_2$-evaluation at a chosen group element.

**Verdict on (b):** the embedding is categorically blocked. Peter–Weyl matrix coefficients are POLYNOMIAL FUNCTIONS on a group; OffDiag transitions are OPERATORS with SCALAR cocycle data. The two live in different categories.

---

## 4. The categorical tensor product (case (a)) IS the combined substrate

### 4.1 Combined Hopf-algebra structure

$$
\mathcal{A}^{\mathrm{combined}}_{n_{\max}, j_{\max}} := \mathcal{A}^{J^*}_{j_{\max}}(S^3) \otimes \mathcal{A}^{\mathrm{OD}}_{n_{\max}}
$$
with:
- **Algebra:** componentwise product $(x_1 \otimes y_1)(x_2 \otimes y_2) = x_1 x_2 \otimes y_1 y_2$.
- **Coproduct:** $\Delta^{\mathrm{combined}} = (\Delta^{J^*} \otimes \Delta^{\mathrm{OD}}) \circ (\mathrm{id} \otimes \tau \otimes \mathrm{id})$. On generators:
  - $\Delta(\pi^j_{mn} \otimes 1) = \sum_p (\pi^j_{mp} \otimes 1) \otimes (\pi^j_{pn} \otimes 1)$ — Peter–Weyl matrix-coefficient coproduct on the left factor.
  - $\Delta(1 \otimes e_s) = (1 \otimes e_s) \otimes (1 \otimes 1) + (1 \otimes 1) \otimes (1 \otimes e_s)$ — primitive on the OffDiag idempotents.
  - $\Delta(1 \otimes T_{s'\to s})$ — non-primitive on transitions; exact form open at v3.61.0 Track B level (multi-year continuation).
- **Counit:** $\varepsilon^{\mathrm{combined}}(x \otimes y) = \varepsilon^{J^*}(x) \cdot \varepsilon^{\mathrm{OD}}(y)$.
- **Antipode:** $S^{\mathrm{combined}}(x \otimes y) = S^{J^*}(x) \otimes S^{\mathrm{OD}}(y)$ (at the standard $\mathcal{O}(SL_2)$ quotient on the left factor; the OffDiag antipode at finite cutoff inherits the involutive structure of the path algebra).

### 4.2 Dimensional panel at $(n_{\max}, j_{\max}) = (2, 1/2)$

| Factor | Dim |
|:------|:---:|
| $\mathcal{A}^{J^*}_{1/2}$ (Peter–Weyl basic) | 5 |
| $\mathcal{A}^{\mathrm{OD}}_{2}$ basic-generator count (5 idempotents + 12 transitions) | 17 |
| $\mathcal{A}^{\mathrm{OD}}_{2}$ path-algebra closure (5 + 12 + 18 two-step) | 35 |
| L1 combined basic dim | $5 \times 17 = 85$ |
| L1 combined closure dim | $5 \times 35 = 175$ |

### 4.3 M-slot grading + non-abelian content preserved

The combined substrate retains:
- **Mellin-slot M1/M2/M3 partition** on the OffDiag (v3.61.0) factor via the $k$-index on idempotents.
- **Non-abelian content** on the Peter–Weyl factor (matrix-coefficient coproduct is non-cocommutative for $j > 0$).
- **Pro-unipotent Lie structure** from OffDiag transitions (24/66 non-zero commutators, T3b finding) inside the abelian-base factor.

---

## 5. Bit-exact verification: SU(2) does NOT preserve the OffDiag basis

### 5.1 The structural claim

For L1 to be UPGRADEABLE to L5 (a smash product where the SL_2 factor acts non-trivially on the OffDiag substrate), one would need $SU(2)$ to act on the OffDiag substrate preserving its basic basis. We test this directly: does the natural conjugation action of $U_\theta = \mathrm{diag}(e^{i \theta m_J})$ (z-rotation) on End$(\mathcal{H})$ preserve the span of the 12 OffDiag transitions?

### 5.2 4-state model and symbolic computation (bit-exact in `sympy`)

In a 4-state truncation of the CH spectral triple at sector pair $((1,1), (1,0))$, with states ordered as $|a\rangle, |b\rangle$ in sector $(1,0)$ (with $m_J = \pm 1/2$) and $|c\rangle, |d\rangle$ in sector $(1,1)$ (with $m_J = \pm 1/2$):

$$
T_{(1,1)\to(1,0)} = \begin{pmatrix} 0 & 0 & \kappa & \kappa \\ 0 & 0 & \kappa & \kappa \\ 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 \end{pmatrix}, \qquad \kappa = -1/16
$$

(uniform-$\kappa$ entries on all E1-allowed pairs satisfying $|\Delta m_J| \le 1$).

Under $U_\theta = \mathrm{diag}(e^{i\theta/2}, e^{-i\theta/2}, e^{i\theta/2}, e^{-i\theta/2})$ (z-rotation by $m_J$):

$$
U_\theta T U_\theta^{-1} = \begin{pmatrix} 0 & 0 & -1/16 & -e^{i\theta}/16 \\ 0 & 0 & -e^{-i\theta}/16 & -1/16 \\ 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 \end{pmatrix}.
$$

The ratios $T_{\mathrm{conj}}[i,j] / T[i,j]$ for nonzero entries are:

| $(i, j)$ | Ratio |
|:--------:|:-----:|
| $(0, 2)$ | $1$ |
| $(0, 3)$ | $e^{i\theta}$ |
| $(1, 2)$ | $e^{-i\theta}$ |
| $(1, 3)$ | $1$ |

**Three distinct phase factors $\{1, e^{i\theta}, e^{-i\theta}\}$ appear.** Bit-exact symbolic verification: `T_conj` is **NOT** a scalar multiple of $T$. The conjugate lives in the span of $m_J$-RESOLVED sub-transitions $T^{(m_J, m_{J'})}_{s'\to s}$, NOT in the span $\{T_{s'\to s}\}$ of the basic 5-sector OffDiag basis.

### 5.3 Structural consequence

For the SU(2) factor of the combined substrate to act non-trivially on the OffDiag factor (i.e., for L5 to be strictly richer than L1 via a smash product), the OffDiag substrate would need to be REFINED to the $m_J$-resolved sub-sector basis. This is a categorically distinct substrate from v3.61.0's basic 5-sector basis — its dimensionality grows from 17 generators (5 idempotents + 12 transitions) to roughly $16 + 16$ (per-state idempotents + per-state-pair transitions), losing the sector-structure that v3.61.0 Track A identified as the load-bearing data.

**Verdict:** the basic OffDiag substrate is NOT SU(2)-invariant; any non-trivial smash product would require refining the OffDiag basis, after which the v3.61.0 sector structure is lost. **Therefore L1 (categorical tensor product) IS the right combination at the substrate level.**

---

## 6. Wigner–Eckart structure of $A$ (sub-layer)

Although L5 reduces to L1 at the substrate level, there IS structural content at a sub-layer: the E1 dipole adjacency $A$ is a Wigner–Eckart spin-1 operator. Its matrix elements within each transition $T_{s'\to s}$ are determined by Clebsch–Gordan coefficients:

$$
\langle (n, l, m) | A | (n', l', m') \rangle = \langle (n, l) || A^{(1)} || (n', l') \rangle \cdot \langle l, m; 1, q | l', m' \rangle,
$$

where $q = m' - m$ and the reduced matrix element $\langle (n, l) || A^{(1)} || (n', l') \rangle$ is forced to a uniform value $\kappa = -1/16$ in the CH convention.

At $l \to l'$ E1 pairs:

| $(l, l')$ | CG factors |
|:---------:|:----------|
| $(0, 1)$ | $\{1\}$ (single non-zero per $q$) |
| $(1, 0)$ | $\{1\}$ (Hermitian conjugate) |
| $(1, 2)$ | $\sqrt{2/5}, \sqrt{3/5}, \ldots$ — non-trivial radicals |
| $(2, 1)$ | conjugate of $(1, 2)$ |

**Structural reading:** the CG factors describe the INTERNAL structure of each transition $T_{s'\to s}$ (which matrix entries within the off-diagonal block are non-zero, and their relative ratios). They do NOT describe relations BETWEEN transitions. So the SU(2) representation-theoretic structure of $A$ is *internal* to each OffDiag generator, not a smash-product-style action across generators.

This is the precise sense in which "L5 reduces to L1 with additional internal CG structure" — the per-transition matrix block has SU(2) content, but the inter-transition relations are flat (tensor product).

---

## 7. Motivic Galois group structure

### 7.1 Combined substrate at finite cutoff

$$
\boxed{
U^{*(n_{\max}, j_{\max})}_{\mathrm{combined}} \;=\; SL_2 \times \mathbb{G}_a^{N_{\mathrm{OD}}(n_{\max})}
}
$$

with:
- $SL_2$ factor from T3a's $\mathcal{O}(SU(2))$ quotient of the Peter–Weyl matrix-coefficient algebra.
- $\mathbb{G}_a^{N_{\mathrm{OD}}(n_{\max})}$ pro-unipotent factor from the OffDiag idempotents + transitions, with v3.61.0's M-slot grading on the abelian-base subfactor.
- T3b's non-abelian Lie content (24/66 non-zero commutators) sits INSIDE the pro-unipotent factor as the explicit Lie algebra of transitions.

### 7.2 Comparison to L1's Levi-decomposition

L1 (the v3.61.0 + T3a tensor product, $\mathcal{H}_{\mathrm{GV}} \otimes \mathcal{O}(SL_2)$) gives $U^* = \mathbb{G}_a^{3N(n_{\max})} \times SL_2$. The two factors $SL_2$ and $\mathbb{G}_a$ commute (categorical tensor product = direct product of algebraic groups at the spec level).

L5 (this sprint's combined OffDiag-enriched substrate) gives $U^* = SL_2 \times \mathbb{G}_a^{N_{\mathrm{OD}}}$. The factors again commute.

**L5 = L1 categorically** at the group structure level. The difference between L1 and L5 is in the IDENTITY of the pro-unipotent factor: L1 uses v3.61.0's M-slot-graded $\mathbb{G}_a^{3N(n_{\max})}$; L5 uses T3b's OffDiag-enriched $\mathbb{G}_a^{N_{\mathrm{OD}}}$. These two are themselves related by a substrate refinement (each generator carries additional Wigner–Eckart-determined block structure in L5 vs L1).

### 7.3 The multi-year synthesis target

The natural multi-year target is the Levi-decomposition product

$$
U^*_{\mathrm{Stage-2}} = SL_2 \ltimes \mathbb{G}_a^{3N(n_{\max})}_{\mathrm{M-slot}, \mathrm{OffDiag-enriched}}
$$

— with the semidirect product structure (if a non-trivial SL_2 action on the pro-unipotent factor is identified at the REFINED OffDiag basis level). At the basic 5-sector level, the semidirect product collapses to a direct product (this sprint's finding). The genuinely-new structural content at the smash-product level lives at the $m_J$-RESOLVED refinement, which is the natural multi-year continuation.

---

## 8. Honest scope

### 8.1 Closed at theorem grade (bit-exact at finite cutoff)

- **Categorical mismatch test:** Peter–Weyl matrix coefficients live in $\mathcal{O}(SL_2)$ (polynomial functions on group); OffDiag transitions live in End($\mathcal{H}$) (operators with scalar cocycle data). No group element $g_0 \in SL_2$ evaluates the 4 spin-$1/2$ matrix coefficients to the 12 OffDiag $\eta$-class rationals simultaneously.
- **SU(2) z-rotation non-invariance:** symbolic bit-exact `sympy` computation on a 4-state model: $U_\theta T U_\theta^{-1}$ has phase factors $\{1, e^{i\theta}, e^{-i\theta}\}$ on the off-diagonal block — not a scalar multiple of $T$.
- **Levi-decomposition shape:** $U^*_{\mathrm{combined}} = SL_2 \times \mathbb{G}_a^{N_{\mathrm{OD}}}$ at the categorical tensor product level.
- **L1 = L5 categorically** at the substrate level.

### 8.2 Open at multi-year continuation

- **Smash product upgrade.** If the OffDiag basis is refined to $m_J$-resolved sub-idempotents, does the SU(2) action become non-trivial? At what cutoff does the semidirect product structure emerge? Multi-year derivation against the case-exhaustion theorem (Paper 32 §VIII).
- **Bit-exact bridge from L1 to v3.61.0 Track B.** Does the L1 tensor product reproduce the $\pm 1/65536$ closure drift residual of v3.61.0 Track B at the combined substrate level? T3b's structural alignment hint (modulo $2^3$ JLO simplex factor) is the starting point.
- **Stage-2 Tannakian construction.** Whether the pro-affine algebraic group obtained from the inverse-limit Hopf algebra on the combined substrate admits the Connes–Marcolli cosmic-Galois structure. Same open question as in T3a §7.2 and v3.61.0 §10.

### 8.3 Curve-fit audit (`feedback_audit_numerical_claims`)

The categorical mismatch finding (Peter–Weyl matrix coefficients are functions on a group; OffDiag transitions are scalar cocycle data) is a textbook categorical observation, not a numerical coincidence — no PSLQ, no fitted coefficient. The SU(2) z-rotation non-invariance is a direct symbolic computation in `sympy` with three distinct phase factors on the panel; verdict is forced by the categorical structure, not selection-biased.

### 8.4 Transcendental tagging (`feedback_tag_transcendentals`)

At the substrate level, no transcendentals appear: matrix coefficients are polynomial over $\mathbb{Z}$, OffDiag cocycle values are rationals in $\kappa^2 \cdot \mathbb{Z}$. The CG coefficients introduce degree-2 algebraic radicals ($\sqrt{2/5}$, etc.) at the sub-layer — these are integers under the standard SU(2) representation theory and stay strictly algebraic. No $\pi$ at the skeleton level (consistent with `feedback_discrete_for_skeleton`).

### 8.5 Discrete-for-skeleton (`feedback_discrete_for_skeleton`)

All verifications bit-exact in `sympy.Rational` + symbolic SU(2). No floats, no PSLQ.

### 8.6 WH1 PROVEN unaffected

This sprint constructs a categorical scoping at the Hopf-algebra substrate level; does not test propinquity convergence or modify the WH1 / Marcolli–vS lineage closure.

### 8.7 Hard prohibitions (§13.5)

No changes to natural geometry hierarchy. No fitted/empirical parameters. No deletion of negative results from §3. No removal of "conjectural" label from Paper 2 combination rule (Paper 2 not touched).

---

## 9. Files

### Produced
- `debug/compute_q5p_combined_substrate.py` — Part 1 driver (~250 lines, 0.008 s wall, structural analysis + Wigner–Eckart CG panel).
- `debug/compute_q5p_combined_substrate_part2.py` — Part 2 driver (~100 lines, ~0.3 s wall, bit-exact symbolic SU(2) z-rotation conjugation test).
- `debug/data/sprint_q5p_combined_substrate.json` — Part 1 data.
- `debug/data/sprint_q5p_combined_substrate_part2.json` — Part 2 data with explicit phase factors.
- `debug/sprint_q5p_combined_substrate_memo.md` — this memo.

### Used (load-bearing inputs)
- `debug/sprint_q5p_j_star_s3_memo.md` (T3a Peter–Weyl substrate).
- `debug/sprint_q5p_offdiag_dirac_memo.md` (T3b off-diagonal Dirac substrate).
- `debug/data/sprint_q5p_offdiag_dirac.json` (12 OffDiag transition $\eta$-class values).
- `debug/sprint_q5p_hard_parts_round2_2026_06_06_memo.md` (v3.62.0 umbrella — the joint structural-hierarchy reading).

### Published references
- Klimyk, A.; Schmüdgen, K. *Quantum Groups and Their Representations.* Springer (1997), §1.3.2 and Prop 1.3.5 (Hopf algebra tensor product).
- Connes, A.; Kreimer, D. *"Hopf algebras, renormalization and noncommutative geometry."* CMP 199 (1998), 203-242. (Original nested-Hopf-tower target.)
- Connes, A.; Marcolli, M. *"Noncommutative Geometry, Quantum Fields and Motives."* AMS Colloquium Publications 55 (2008), Ch. 4. (Cosmic-Galois $U^*$ on Hopf algebras; the Levi-decomposition structure of general algebraic groups acting via Tannakian formalism.)
- Bröcker, T.; tom Dieck, T. *Representations of Compact Lie Groups.* Springer (1985), Ch. III. (Wigner–Eckart theorem; SU(2) representation theory.)

---

## 10. Paper-edit recommendations (PI to apply — apply NONE in this sprint)

### 10.1 Paper 32 §VIII — ONE new Remark `rem:q5p_combined_substrate_levi` after `rem:q5p_offdiag_dirac_enrichment`

```latex
\begin{rem}[Q5' Stage-2 combined substrate L5 = L1 categorically,
Sprint Q5'-Combined-Substrate, June 2026]
\label{rem:q5p_combined_substrate_levi}
The first multi-year scoping step of combining the T3a Peter--Weyl
substrate (Remark~\ref{rem:q5p_j_star_substrate}) and T3b cross-shell
off-diagonal Dirac substrate (Remark~\ref{rem:q5p_offdiag_dirac_enrichment})
on a SINGLE combined substrate is CLOSED via categorical reduction.
The combined substrate $\mathcal{A}^{\mathrm{combined}} =
\mathcal{A}^{J^*}_{j_{\max}} \otimes \mathcal{A}^{\mathrm{OD}}_{n_{\max}}$
reduces to the categorical TENSOR PRODUCT (L1) at the basic substrate
level, not a smash product (L5). Two structural obstructions block the
smash-product upgrade: (i) Peter--Weyl matrix coefficients $\pi^j_{mn}$
are polynomial \emph{functions} on $SL_2$, while OffDiag transitions
$T_{s'\to s}$ are operators on $\mathcal{H}$ with \emph{scalar} cocycle
data ($\eta$-class values in $\kappa^2 \cdot \mathbb{Z}$) — categorically
different objects; no $SL_2$ group element evaluates the 4 spin-$1/2$
matrix coefficients to the 12 OffDiag $\eta$-rationals simultaneously.
(ii) The natural $SU(2)$ z-rotation $U_\theta$ does NOT preserve the
basic 5-sector OffDiag basis: bit-exact symbolic computation
(\texttt{debug/compute\_q5p\_combined\_substrate\_part2.py}) gives
three distinct phase factors $\{1, e^{i\theta}, e^{-i\theta}\}$ on the
conjugate $U_\theta T U_\theta^{-1}$ — the result lives in the span of
$m_J$-resolved sub-transitions, not in the span of the basic basis.
The candidate motivic Galois group at the quotient is therefore the
\emph{Levi-decomposition product}
\[
U^{*(n_{\max}, j_{\max})}_{\mathrm{combined}} = SL_2 \times \mathbb{G}_a^{N_{\mathrm{OD}}(n_{\max})},
\]
exactly the multi-year synthesis target predicted by T3a's headline
finding. \textbf{L5 reduces to L1 at the substrate level}; the genuinely
new content (a non-trivial smash product) lives at the $m_J$-resolved
refinement of the OffDiag basis, which loses the v3.61.0 sector
structure as load-bearing data and constitutes a multi-year continuation.
See Paper~55 \S\ref{subsec:open_m2_m3}.
\end{rem}
```

### 10.2 Paper 55 §subsec:open_m2_m3 — ONE new paragraph after the OffDiag enrichment paragraph

```latex
\emph{Combined substrate scoping: L5 = L1 categorically (Sprint
Q5'-Combined-Substrate, June 2026; memo
\texttt{debug/sprint\_q5p\_combined\_substrate\_memo.md}; data
\texttt{debug/data/sprint\_q5p\_combined\_substrate.json}).} The first
multi-year scoping step of combining T3a (Peter--Weyl $J^*(S^3)$) and
T3b (cross-shell off-diagonal Dirac) on a SINGLE combined substrate
closes the categorical question: \textbf{at the basic-substrate level, L5
reduces to L1 — the categorical tensor product $\mathcal{A}^{J^*} \otimes
\mathcal{A}^{\mathrm{OD}}$}. Two clean structural obstructions block the
smash-product upgrade: (i) Peter--Weyl matrix coefficients are polynomial
functions on $SL_2$ (in $\mathcal{O}(SL_2)$); OffDiag transitions are
operators on the CH Hilbert space with scalar cocycle data — categorically
different objects, no consistent identification at the panel level. (ii)
The $SU(2)$ z-rotation $U_\theta$ does NOT preserve the basic 5-sector
OffDiag basis — bit-exact symbolic computation gives three distinct
phase factors $\{1, e^{i\theta}, e^{-i\theta}\}$ on the conjugate
$U_\theta T U_\theta^{-1}$; the result lives in the span of $m_J$-resolved
sub-transitions. The candidate motivic Galois group at the combined
substrate quotient is therefore the \emph{Levi-decomposition product}
$U^*_{\mathrm{combined}} = SL_2 \times \mathbb{G}_a^{N_{\mathrm{OD}}}$
— the multi-year synthesis target predicted by T3a's headline. The
genuinely-new smash-product content lives at the $m_J$-resolved
refinement, named as a multi-year continuation. This sprint closes the
categorical question on the combined substrate while preserving the
Levi-decomposition structural target.
```

### 10.3 Paper 18 — no edit needed

§III.7 master Mellin engine is upstream; the substrate-level categorical question is downstream from the projection-mechanism partition.

### 10.4 Paper 38 — no edit needed

Paper 38 SU(2) Peter–Weyl substrate is used at the propinquity / GH-convergence level, not at the Hopf-algebra substrate-enrichment level.

---

## 11. One-line verdict

**POSITIVE-REDUCES-TO-L1.** The combined substrate $\mathcal{A}^{J^*} \otimes \mathcal{A}^{\mathrm{OD}}$ reduces to the categorical tensor product L1 at the basic-substrate level, blocked from upgrade to a smash product L5 by two clean structural obstructions: (i) Peter–Weyl matrix coefficients are polynomial functions on $SL_2$ while OffDiag transitions are operators with scalar cocycle data — categorically distinct objects; (ii) the $SU(2)$ z-rotation does NOT preserve the basic 5-sector OffDiag basis (bit-exact symbolic test gives phase factors $\{1, e^{i\theta}, e^{-i\theta}\}$ on the conjugate). The candidate motivic Galois group at the combined-substrate quotient is the Levi-decomposition product $U^*_{\mathrm{combined}} = SL_2 \times \mathbb{G}_a^{N_{\mathrm{OD}}}$ — exactly the multi-year synthesis target predicted by T3a's headline. The genuinely-new smash-product content lives at the $m_J$-resolved refinement of the OffDiag basis, which loses the v3.61.0 sector structure as load-bearing data — multi-year continuation.
