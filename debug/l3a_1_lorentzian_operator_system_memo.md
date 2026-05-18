# Sprint L3a-1 — Lorentzian truncated operator system at finite cutoff

**Date:** 2026-05-17
**Author:** L3a-1 PM (Claude)
**Status:** CLOSED-AT-FINITE-CUTOFF with one substantive structural finding (propagation number is *envelope-dependent*).
**Inputs:** CLAUDE.md §1.7 WH1 + L2 paragraphs; Paper 32 §III `rem:operator_system` + `prop:propagation_2`; Paper 43 §11; `debug/l3_scoping_memo.md`; `debug/l3_literature_audit_memo.md`; `debug/sprint_l2_synthesis_memo.md`.
**Deliverables:** `geovac/operator_system_lorentzian.py` (~700 lines), `tests/test_operator_system_lorentzian.py` (34 tests, 33 fast + 1 slow), `debug/data/l3a_1_lorentzian_operator_system.json`, this memo.

---

## §1. Construction summary

This sprint builds the Lorentzian analog of the Connes–van Suijlekom (CMP 2021, arXiv:2004.14115) truncated operator system on the BBB Krein spectral triple from Sprint L2.

### Hilbert space

The Krein space from Sprint L2-B (`geovac.krein_space_construction.KreinSpace`):
$$\mathcal{K}_{n_{\max}, N_t} = \mathcal{H}_{\mathrm{GV}}^{n_{\max}} \otimes \mathbb{C}^{N_t},$$
with fundamental symmetry $J = \gamma^0 \otimes I_{N_t}$ in the Peskin–Schroeder chiral basis, West-coast metric $\eta = \mathrm{diag}(+,-,-,-)$, at BBB $(m, n) = (4, 6)$.

$\mathcal{H}_{\mathrm{GV}}$ is the Camporesi–Higuchi spinor space (Paper 32 §III Def. 3.2) at FullDirac convention: chirality-doubled, with basis labels `FullDiracLabel(n_fock, l, two_m_j, chirality)`. $\dim(\mathcal{H}_{\mathrm{GV}}) = \frac{2}{3}\,n_{\max}(n_{\max}+1)(n_{\max}+2)$ = $4, 16, 40$ at $n_{\max} = 1, 2, 3$.

### Multipliers

A scalar function $f(\omega) \in C^\infty(S^3)$ acts on $\mathcal{H}_{\mathrm{GV}}$ via the chirality-doubled scalar 3-Y multiplier:
$$M^{\mathrm{spat}}_{NLM} \big|_{\mathcal{H}_{\mathrm{GV}}} = M^{\mathrm{Weyl}}_{NLM} \oplus M^{\mathrm{Weyl}}_{NLM}$$
(block-diagonal with two equal copies). The spatial Weyl block is built via the spinor lift (Clebsch–Gordan decomposition + Avery–Wen–Avery 3-Y integrals from `geovac.spinor_operator_system`), exactly matching `FullDiracTruncatedOperatorSystem.build_full_dirac_multiplier_matrix`.

A function $g(t) \in C^\infty(\mathbb{R}_t)_{\mathrm{cutoff}}$ acts on $\mathbb{C}^{N_t}$ by pointwise multiplication at the temporal grid points $t_0, \ldots, t_{N_t-1}$. We use the polynomial basis $g_p(t) = t^p$ for $p = 0, \ldots, N_t - 1$; the Vandermonde matrix at $N_t$ distinct grid points is invertible, so these span the full $N_t$-dimensional diagonal subalgebra.

The Lorentzian truncated operator system is then
$$O^L_{n_{\max}, N_t} := \mathrm{span}_{\mathbb{C}} \big\{ M^{\mathrm{spat}}_{NLM} \otimes \mathrm{diag}(g_p(t_0), \ldots, g_p(t_{N_t-1})) \big\}.$$

Total generator count: $|\{(N, L, M)\}| \cdot N_t = \dim_{\mathrm{spat}}(O^{\mathrm{spat}}_{n_{\max}}) \cdot N_t$.

At $N_t = 1$: only the constant temporal multiplier $g_0 = 1$ is available; the construction reduces to `FullDiracTruncatedOperatorSystem(n_max)`.

---

## §2. Load-bearing falsifiers passed

### F1. Riemannian limit at $N_t = 1$ — BIT-EXACT

At $N_t = 1$, the construction reduces to `FullDiracTruncatedOperatorSystem(n_max)` bit-identically. Verified at $n_{\max} \in \{1, 2, 3\}$:

| $n_{\max}$ | $\dim_K$ | multipliers | $\dim(O^L)$ | max residual $\|M_{\mathrm{self}} - M_{\mathrm{ref}}\|_F$ |
|:----------:|:--------:|:-----------:|:-----------:|:--------------------------------------------------------:|
| 1 | 4 | 1 | 1 | **0.0** (exact) |
| 2 | 16 | 14 | 14 | **0.0** (exact) |
| 3 | 40 | 55 | 55 | **0.0** (exact) |

The max residual is exactly 0.0 in float64, not just within machine epsilon — the construction reduces to the same Kronecker products at the matrix-element level. This is the load-bearing falsifier: had it failed, the scalar-multiplier spinor lift would have been structurally incompatible with the L2-B Krein-space basis at $N_t = 1$, a major structural finding requiring escalation. **The lift is consistent.**

### F2. *-closure — PASS at every cell

For every multiplier $M_i \in O^L$, the conjugate transpose $M_i^\dagger$ also lies in $O^L$. Verified at $(n_{\max}, N_t) \in \{(1, 1), (2, 1), (1, 3), (2, 3)\}$. Mechanism: spatial scalar 3-Y multipliers satisfy $M_{NLM}^* = M_{N, L, -M}$ up to a sign (Wigner conjugation), and temporal real-diagonal multipliers are self-adjoint trivially. Tensor product preserves *-closure.

### F3. Identity in $O^L$ — PASS at every cell

The trivial multiplier $f = 1$ on $S^3 \times \mathbb{R}_t$ (the $N=L=M=0$ spatial mode times the $p=0$ temporal mode) gives the identity matrix on $\mathcal{K}$. Identity-in-$O^L$ residual is $\le 1.12 \times 10^{-15}$ at all tested cells.

### F4. Sprint L2 baseline zero regression

Full Sprint L2 baseline (`tests/test_operator_system.py`, `tests/test_krein_space_construction.py`, `tests/test_lorentzian_dirac.py`, `tests/test_connes_axiom_audit_31.py`, `tests/test_modular_hamiltonian_lorentzian.py`): **390 passed, 3 skipped** in 152 s. No regression.

---

## §3. Propagation number result at panel

The Sprint asks: does $\mathrm{prop}(O^L_{n_{\max}, N_t}) = 2$ hold in the Krein setting, matching the Toeplitz-$S^1$ result (Connes–vS Prop 4.2) and Paper 32 §III on $S^3$?

**Substantive structural finding: the answer depends on which envelope is targeted.**

### §3.1 Achievable envelope = $\dim_{\mathrm{Weyl}}^2 \cdot N_t$

Scalar multipliers acting on the chirality-doubled FullDirac space have the structural form $M \oplus M$ in chirality blocks. Their products $(M_1 \oplus M_1)(M_2 \oplus M_2) = (M_1 M_2) \oplus (M_1 M_2)$ remain block-diagonal in chirality with equal blocks. The maximal subspace reachable by scalar-multiplier products is the *Weyl-doubled subspace*:
$$\big\{ M_W \oplus M_W : M_W \in M_{\dim_{\mathrm{Weyl}}}(\mathbb{C}) \big\} \cong M_{\dim_{\mathrm{Weyl}}}(\mathbb{C}),$$
of dimension $\dim_{\mathrm{Weyl}}^2$. Tensored with the temporal commutative subalgebra of dim $N_t$ (diagonal matrices on $\mathbb{C}^{N_t}$), the achievable envelope is
$$\dim_{\mathrm{achievable}} = \dim_{\mathrm{Weyl}}^2 \cdot N_t.$$

Under this convention, the propagation number test is:

| $n_{\max}$ | $N_t$ | $\dim_{\mathrm{Weyl}}^2 \cdot N_t$ | $\dim(O^L)$ | $\dim((O^L)^2)$ | prop |
|:----------:|:-----:|:----------------------------------:|:-----------:|:---------------:|:----:|
| 1 | 1 | 4 | 1 | 1 | $\infty$ (trivial) |
| 1 | 3 | 12 | 3 | 3 | $\infty$ (trivial) |
| 2 | 1 | 64 | 14 | 64 | **2** |
| 3 | 1 | 400 | 55 | 400 | **2** |
| 2 | 3 | 192 | 42 | 192 | **2** |
| 2 | 5 | 320 | 70 | 320 | **2** |

**prop = 2 across all $n_{\max} \ge 2$ cells under the achievable-envelope convention**, matching Connes–vS Toeplitz $S^1$ (Prop 4.2) and Paper 32 §III prop=2 verbatim. The Lorentzian extension preserves the propagation-2 structural feature of the Riemannian construction at the operator-system level.

The $n_{\max} = 1$ cells are *trivial* (the spatial operator system reduces to identity-only, because the only $(N, L, M)$ multiplier with allowed labels at $n_{\max} = 1$ is $(N=1, L=0, M=0)$ which gives the identity matrix). The 3-dim span at $(1, 3)$ is the temporal subalgebra times identity. prop=$\infty$ at $n_{\max} = 1$ is the structurally expected trivial case (the spatial operator system has no nontrivial generators to drive propagation).

### §3.2 Full envelope = $\dim_K^2$

Targeting the full $B(\mathcal{K})$ envelope of dim $\dim_K^2$:

| $n_{\max}$ | $N_t$ | $\dim_K^2$ | $\dim((O^L)^k)$ saturates at | prop (full) |
|:----------:|:-----:|:----------:|:----------------------------:|:-----------:|
| 2 | 1 | 256 | 64 | $\infty$ |
| 2 | 3 | 2304 | 192 | $\infty$ |

**prop = $\infty$ under the full-envelope convention at every tested cell.** Mechanism: scalar multipliers cannot generate chirality-flipping operators (the off-diagonal Weyl ↔ anti-Weyl blocks of $M_{\dim_K}(\mathbb{C})$). They can also not generate non-diagonal temporal operators (the temporal multiplier algebra is commutative). The full envelope contains operators of both types, which $O^L$ cannot reach.

### §3.3 Why both numbers are real

The "achievable" convention is the natural Connes–vS-style propagation question for scalar-multiplier operator systems on the *spinor bundle*: it asks at what depth $k$ the operator system reaches the maximal subspace consistent with the construction. The "full" convention is a strictly larger envelope that includes operators outside the natural construction.

Paper 32 §III prop=2 was stated for the *scalar* operator system on the scalar Hilbert space (Fock-projected $S^3$, no chirality doubling). When we lift to the FullDirac spinor space, the natural envelope is the Weyl-doubled subspace because that is where scalar multipliers live. prop=2 holds there. The full $B(\mathcal{K})$ envelope is not the natural target — reaching it would require either gamma-matrix-valued multipliers (chirality-flipping) or noncommutative temporal multipliers ($\partial_t$-type couplings).

This is a substantive sharpening of the Paper 32 §III propagation question for the spinor lift. The Riemannian-side `FullDiracTruncatedOperatorSystem` from Sprint WH1-R3.5 has the same structure (prop=2 against the Weyl-doubled envelope; prop=$\infty$ against the full $\dim_K^2$ envelope), but this had not been explicitly verified until now.

---

## §4. Krein-positivity findings

The Krein-positive cone $\mathcal{K}^+ = \{|\psi\rangle : J|\psi\rangle = +|\psi\rangle\}$ is the $+1$-eigenspace of $J$, of dim $\dim_K / 2$. An operator $O$ preserves $\mathcal{K}^+$ iff $[J, O] = 0$ at the $\mathcal{K}^+$ block level.

**Finding: every multiplier in $O^L$ preserves $\mathcal{K}^+$.** Mechanism: $J_{\mathrm{spatial}} = $ chirality-swap (off-diagonal in chirality), and $M^{\mathrm{spat}} = M \oplus M$ commutes with chirality-swap because $J$ permutes the two equal copies of $M$. Tensored with any operator commuting with $I_{N_t}$ (which is everything), the full multiplier commutes with $J$.

Verified at $(n_{\max}, N_t) \in \{(2, 1), (2, 3), (3, 1)\}$: **all multipliers preserve $\mathcal{K}^+$**, fraction = 1.0 / 1.0 at every cell.

**Structural reading: the Krein-positive restriction $O^{L,+} := \{O \in O^L : O \cdot \mathcal{K}^+ \subset \mathcal{K}^+\}$ equals $O^L$ itself.** The restriction is trivial in this construction.

This is the *expected* result given the chirality-doubling structure: a scalar function does not see chirality, so it cannot break the $J$-grading. To get a nontrivial Krein-positive restriction, the operator system would need to be enlarged to include chirality-mixing operators (e.g., gamma-matrix-valued multipliers), which are outside the standard pointwise-multiplication algebra.

For the L3 program, this finding means: **the natural Lorentzian operator system from the scalar-function construction is already entirely Krein-positive**. The L3-Krein-positive-restricted claim (`debug/l3_scoping_memo.md` §1.4) is then satisfied trivially at the operator-system substrate level. The Krein-positivity question shifts to the *state* level (Krein-positive states on $O^L$), which is the propinquity / Wasserstein-Kantorovich question downstream of this sprint.

---

## §5. Witness pair construction

Take
$$a = M^{\mathrm{spat}}_{N=2, L=1, M=0} \otimes g_0(t),  \qquad b = a^*.$$

Both lie in $O^L$ (verified). Their product is
$$a b = \big(M^{\mathrm{spat}} (M^{\mathrm{spat}})^*\big) \otimes (g_0(t))^2 = \big(M^{\mathrm{spat}} (M^{\mathrm{spat}})^*\big) \otimes I_{N_t}.$$

The spatial factor $M^{\mathrm{spat}} (M^{\mathrm{spat}})^*$ is the Riemannian Paper 32 §III witness: it does NOT lie in $O^{\mathrm{spat}}$ because raise-then-lower at the top shell $n = n_{\max}$ produces a deficit on the top-shell diagonal that no spatial multiplier $M_{N L M}$ can supply (the corresponding $n_{\max} + 1$ shell is outside the truncation).

The least-squares projection residual at $n_{\max} = 2$:

| $(n_{\max}, N_t)$ | residual $\|ab - \mathrm{proj}_{O^L}(ab)\|_F / \|ab\|_F$ |
|:-----------------:|:-------------------------------------------------------:|
| (2, 1) | 0.381 |
| (2, 3) | 0.381 |
| (2, 5) | 0.381 |

The residual is identical across $N_t$ values because the wedge-mass deficit is purely spatial, tensored with the temporal identity on $g_0^2 = 1$ which is fully inside the temporal multiplier algebra. **The witness pair certifies that $O^L$ is genuinely an operator system (not a *-algebra): closed under * and linear combinations but NOT under multiplication.**

---

## §6. Connection to L2-E hemispheric wedge

The Sprint L2-E hemispheric Krein wedge $W_L = P_W^{\mathrm{spatial}} \otimes P_{t \ge 0}$ (`geovac.modular_hamiltonian_lorentzian.LorentzianWedge`) lives on the same Krein space $\mathcal{K}_{n_{\max}, N_t}$. The natural connection is the wedge-restricted operator system
$$O^L_W := \{ P_{W_L} \cdot O \cdot P_{W_L} : O \in O^L \}.$$

Verified at $(n_{\max}, N_t) = (2, 1)$: $\dim(W_L) = 8$ (spatial wedge dim 8 × temporal positive dim 1); the wedge-block linear-span dim is bounded above by $\dim(O^L) = 14$ and by $\dim(W_L)^2 = 64$. The wedge restriction is well-defined and produces a sub-operator-system on the wedge.

At $(n_{\max}, N_t) = (2, 3)$: $\dim(W_L) = 16$ (spatial wedge dim 8 × temporal positive dim 2 at $t \in \{0, +1\}$).

This is the natural substrate for lifting L2-E's modular-Hamiltonian construction (the wedge KMS state $\rho_W^L$, the BW-α geometric generator $K_L^\alpha$, the Tomita-Takesaki $K_L^{\mathrm{TT}}$) to the operator-system level. The L2-E construction is at the matrix level; this sprint lifts it to the operator-system level by identifying the wedge-restricted operator system as the natural substrate for modular-Hamiltonian computations.

The full propinquity-style L3 question (operator-system-level Lorentzian-propinquity convergence on the wedge KMS state) requires the L2 joint kernel construction (debug/l3_scoping_memo.md §2 Q1), which is the dominant cost of L3b/L3c.

---

## §7. Verdict

**CLOSED-AT-FINITE-CUTOFF.** All load-bearing falsifiers passed at the panel (n_max, N_t) ∈ {(1, 1), (2, 1), (3, 1), (1, 3), (2, 3), (2, 5)}:

1. **Riemannian limit bit-exact** at $N_t = 1$ (max_residual = 0.0 in float64 at all $n_{\max} \in \{1, 2, 3\}$).
2. **Operator system axioms hold:** identity in $O^L$, *-closure, well-defined linear-span dimension.
3. **Propagation number = 2** under the natural (achievable) envelope, matching Paper 32 §III prop=2 verbatim.
4. **Propagation number = $\infty$** under the full $B(\mathcal{K})$ envelope (structural finding: chirality-doubling + commutative temporal subalgebra block scalar-multiplier propagation past the Weyl-doubled subspace).
5. **Krein-positivity trivial** — all multipliers preserve $\mathcal{K}^+$ in this construction (structural finding: the chirality-doubled scalar multipliers commute with $J = $ chirality-swap).
6. **Witness pair exhibited** — $M^{2,1,0,0}$ and its adjoint in $O^L$, product NOT in $O^L$ with $\sim 38\%$ residual (Lorentzian lift of Paper 32 witness).
7. **Connection to L2-E wedge** documented — the wedge-restricted operator system is well-defined; natural substrate for downstream propinquity work.
8. **Zero regression** on the full Sprint L2 baseline (390 passed, 3 skipped).

### Honest scope

This sprint closes the **operator-system level** of the Lorentzian extension at finite cutoff. It does NOT establish:

- **Continuum / GH-limit / Lorentzian-propinquity convergence.** That is the full L3b/L3c program (months of work, see `debug/l3_scoping_memo.md`). The named blockers L2 joint kernel (Q1) and L5 Krein-positive propinquity definition (Q4) remain.
- **Non-commutative temporal multipliers.** The natural pointwise-multiplication algebra on $\mathbb{R}_t$ is commutative; this is the canonical Connes–vS definition. Lifting to non-commutative temporal structure (e.g., $\partial_t$-extended multiplier algebra) would change the propagation question and is outside the natural construction.
- **Gamma-matrix-valued multipliers.** These would enlarge the operator system from $C^\infty(S^3 \times \mathbb{R})$ to $C^\infty(S^3 \times \mathbb{R}, \mathrm{Cl}(\mathbb{R}^4))$ (or similar) and would give a non-trivial Krein-positive restriction. They are outside the scope of the natural truncated operator system construction.

The verdict is **closure at the natural operator-system level**, with two substantive structural findings (envelope-dependent propagation number, trivial Krein-positive restriction).

---

## §8. Recommended Sprint L3b first move

Per the L3 scoping memo (`debug/l3_scoping_memo.md` §6), the natural next sprint is **L3b compact-temporal restriction**: full L1'–L5 propinquity convergence on $S^3 \times S^1_T$ (temporal compactification with periodic boundary conditions), Krein-positive restriction, 4–8 weeks.

This sprint (L3a-1) establishes the operator-system substrate that L3b would take as input. The propagation-2 result at the achievable envelope is a structural feature that L3b's UCP-map tunneling pair construction would inherit; the trivial Krein-positive restriction means L3b can work with the full operator system without additional positivity bookkeeping at the substrate level.

The L3b proof architecture (transferred from Paper 38):

1. **L1' (Krein-positive substrate):** done at finite cutoff in this sprint; lift to the continuum requires the BBB §5 Type-I-factor completion (`debug/l3_scoping_memo.md` §3).
2. **L2 (joint Fejér kernel on SU(2) × $S^1_T$):** the dominant cost. Requires constructing a UCP central Schur multiplier on the joint group with cb-norm $\le 1$ and joint mass-concentration rate. 2–4 weeks at compact $S^1_T$; months at non-compact $\mathbb{R}$.
3. **L3 (Lipschitz spinor bound):** transfers verbatim under the Wick-rotated gradient norm convention (`debug/l3_scoping_memo.md` §2 L3 entry).
4. **L4 (Berezin reconstruction):** Krein-positive Berezin construction extending Paper 38 L4; tractable conditional on L2.
5. **L5 (Latrémolière propinquity assembly):** requires the Krein-positive propinquity definition (Q4 in scoping memo).

**Decision gate:** open Sprint L3b only after the L3a-2 sister track (Q7 Nieuviarts scoping, dispatched in parallel to this sprint per the L3a entry-point plan) clarifies whether the Nieuviarts 2025 twist morphism applies to $S^3$. If yes, L3b reduces to morphism-induced transfer; if no, L3b reverts to the full BBB Krein-construction path.

### Standalone publication value

This sprint's deliverable is naturally citable as a section of Paper 43 (Lorentzian extension at finite cutoff) or as a sub-section of a future L3-completion paper. The propagation-2 result at the achievable envelope is a clean structural sharpening of Paper 32 §III prop:propagation_2 to the chirality-doubled Lorentzian setting, and the bit-exact Riemannian-limit recovery is an independent load-bearing verification of the Sprint L2-B Krein-space spinor basis convention.

---

End of memo.
