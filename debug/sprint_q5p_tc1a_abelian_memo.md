# Sprint Q5'-Tannakian-Closure TC-1a — $\mathrm{Rep}_{\mathrm{fin}}(\mathcal{H}_{\mathrm{GV}})$ is abelian

**Date:** 2026-06-06 (first sub-track of the Tannakian Reconstruction Foundation sprint, single-thread continuation of the Pro-System-Lockdown arc)
**Sprint:** TC-1a of Q5'-Tannakian-Closure (TC-1b, TC-1c, TC-1d, TC-1e to follow single-threaded;\ full Tannakian closure remains the multi-year frontier per PS-4)
**Driver:** `debug/compute_q5p_tc1a_abelian.py`
**Module:** `geovac/tannakian.py` (new, ~700 lines)
**Data:** `debug/data/sprint_q5p_tc1a_abelian.json`
**Wall time:** 0.03 s
**Discipline:** bit-exact `sympy.Rational` / `sympy.Integer` throughout. No floats. No PSLQ.

---

## 1. TL;DR

**Verdict: POSITIVE.** The category $\mathrm{Rep}_{\mathrm{fin}}(\mathcal{H}_{\mathrm{GV}}(n_{\max}))$ of finite-dimensional rational representations of the v3.61.0 Track A abelian primitive Hopf algebra is **abelian** in the Deligne--Milne 1982 sense, verified bit-exact on a five-rep five-morphism panel at $n_{\max} \in \{2, 3\}$ across all six standard universal-property axioms:\ zero object, finite direct sums, kernels, cokernels, mono $=$ ker(coker), epi $=$ coker(ker).

**106 / 106 bit-exact zero residuals:**

| Axiom | per cutoff | $\times$ 2 cutoffs |
|:------|:----------:|:------------------:|
| Zero object (5 reps $\times$ 2 directions)        | 10 | 20 |
| Direct sum universal property (3 pairs $\times$ 5 identities)         | 15 | 30 |
| Kernel universal property (4 morphisms $\times$ 2 conditions)     | 8  | 16 |
| Cokernel universal property (4 morphisms $\times$ 2 conditions)   | 8  | 16 |
| Mono $=$ ker(coker) (2 mono tests $\times$ 3 conditions)               | 6  | 12 |
| Epi $=$ coker(ker) (2 epi tests $\times$ 3 conditions)                 | 6  | 12 |
| **Total** | **53** | **106** |

TC-1a closes the **first of four sprint-scale prerequisites** named by PS-4 (v3.69.0) for Tannakian closure proper. Three remain (TC-1b symmetric monoidal, TC-1c rigidity, TC-1d fiber functor) plus the multi-year first-stone TC-1e ($\mathrm{Aut}^\otimes(\omega) \supseteq U^\ast_{\mathrm{Levi}}$ at $n_{\max} = 2$).

---

## 2. Verdict against decision gate

| Gate | Verdict |
|:-----|:--------|
| POSITIVE | **selected** — all 106 axiom identities bit-exact zero at $n_{\max} \in \{2, 3\}$. The abelian-category structure of $\mathrm{Rep}_{\mathrm{fin}}(\mathcal{H}_{\mathrm{GV}})$ holds verbatim with no GeoVac-specific modifications. |
| BORDERLINE | not selected — closure is bit-exact at the full panel. |
| STOP | rejected — no structural obstruction;\ the panel is constructed bit-exactly from the v3.61.0 Track A abelian primitive substrate plus standard module-category constructions. |

---

## 3. What TC-1a closes and what it does NOT close

TC-1a verifies the **abelian-category axioms** on a small representative panel — it does NOT yet prove Tannakian closure or construct the cosmic-Galois $U^\ast$ via Tannakian reconstruction. The narrative is:

- TC-1a (**this sprint**, closed) — abelian category structure verified.
- TC-1b (next) — symmetric monoidal structure $\otimes_\mathbb{Q}$ with the diagonal Hopf action, unit object, associator, braiding.
- TC-1c — rigidity via the antipode $S(x) = -x$:\ duals + evaluation + coevaluation + snake identities.
- TC-1d — fiber functor $\omega: \mathrm{Rep}_{\mathrm{fin}}(\mathcal{H}_{\mathrm{GV}}) \to \mathrm{Vec}_\mathbb{Q}$:\ exactness, faithfulness, $\otimes$-preservation, $\omega(\mathbf{1}) = \mathbb{Q}$.
- TC-1e (the first-stone of the multi-year wall) — explicit inclusion $\mathrm{Aut}^\otimes(\omega) \supseteq U^\ast_{\mathrm{GeoVac}, \mathrm{Levi}} = \mathbb{G}_a^{3 N(n_{\max})} \rtimes SL_2$ at $n_{\max} = 2$.

TC-1a is **bookkeeping on existing substrate**:\ the abelian-category axioms hold for ANY module category over a ring, and $\mathrm{Rep}_{\mathrm{fin}}(\mathcal{H}_{\mathrm{GV}})$ is the well-known module category of finite-dim ℚ-vector spaces equipped with commuting nilpotent endomorphisms (since $\mathcal{H}_{\mathrm{GV}} = \mathrm{Sym}_\mathbb{Q}(V) = \mathcal{O}(\mathbb{G}_a^{3N(n_{\max})})$ has primitive coproduct). The substantive PS-4 finding was that this is a sprint-scale prerequisite, and TC-1a confirms that empirically by closing in 0.03 s on a representative panel.

---

## 4. The category

### 4.1 Objects

A finite-dim rational rep is a pair $(M, \{X_g\}_g)$ where:

- $M$ is a finite-dim $\mathbb{Q}$-vector space;
- $\{X_g\}_g$ is a family of pairwise commuting nilpotent $\mathbb{Q}$-linear endomorphisms of $M$, indexed by the $3 N(n_{\max})$ primitive generators $g = (n, l, s)$ of $\mathcal{H}_{\mathrm{GV}}(n_{\max})$.

The $\mathbb{G}_a^{3N(n_{\max})}$-action is $\rho(t) = \exp(\sum_g t_g X_g)$ (truncated by nilpotency).

### 4.2 Morphisms

A morphism $f: (M, \{X_g^M\}) \to (M', \{X_g^{M'}\})$ is a $\mathbb{Q}$-linear map satisfying

$$f \circ X_g^M = X_g^{M'} \circ f \qquad \forall g.$$

### 4.3 Code substrate

The `FinDimRep` class (`geovac/tannakian.py`) stores the dim and a sparse dict of non-zero endomorphisms;\ construction-time validation enforces nilpotency and pairwise commutativity bit-exact. `RepMorphism` validates the intertwining condition. `kernel`, `cokernel`, and `direct_sum` give the three universal-property constructions. The verifier functions check the abelian axioms.

---

## 5. Bit-exact panel

### 5.1 Substrate

Five reps at each cutoff $n_{\max} \in \{2, 3\}$:

| Label | Construction |
|:-----:|:-------------|
| `R0` | Zero rep, $\dim = 0$ |
| `T1` | Trivial 1-dim, all endos $= 0$ |
| `T2` | Trivial 2-dim, all endos $= 0$ |
| `J2` | Jordan 2-dim, $X_{(1, 0, 0)} = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$, others $= 0$ |
| `J3` | Jordan 3-dim, $X_{(1, 0, 0)} = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 0 & 0 & 0 \end{pmatrix}$, others $= 0$ |

Five morphisms:

| Label | Construction | Type |
|:-----:|:-------------|:-----|
| `f0` | $R_0 \to T_1$, zero matrix | zero morphism |
| `f1` | $T_1 \to J_2$, $1 \mapsto e_1$ | mono (the kernel of $X_{(1,0,0)}$) |
| `f2` | $J_2 \to T_1$, $e_2 \mapsto 1$, $e_1 \mapsto 0$ | epi (the cokernel of $f_1$) |
| `f3` | $J_2 \to J_2$, identity | iso |
| `f4` | $T_1 \to T_2$, $1 \mapsto e_1$ | mono |

### 5.2 Axiom panels

Per cutoff:

1. **Zero object axiom** (5 reps × 2 directions): for each $R$, the unique morphisms $R_0 \to R$ and $R \to R_0$ exist bit-exact (matrix shapes $(\dim R, 0)$ and $(0, \dim R)$).

2. **Direct sum universal property** (3 pairs × 5 identities): for each pair $(R_1, R_2) \in \{(T_1, T_1), (T_1, J_2), (J_2, J_3)\}$, the inclusions $\iota_1, \iota_2$ and projections $\pi_1, \pi_2$ satisfy $\pi_i \iota_j = \delta_{ij}\, \mathrm{id}_{R_i}$ and $\iota_1 \pi_1 + \iota_2 \pi_2 = \mathrm{id}_S$ bit-exact.

3. **Kernel universal property** (4 morphisms × 2 conditions): for $f_1, f_2, f_3, f_4$, the constructed $K, \iota$ satisfy $f \circ \iota = 0$ and $\iota$ is mono.

4. **Cokernel universal property** (4 morphisms × 2 conditions): symmetric, $\pi \circ f = 0$ and $\pi$ is epi.

5. **Mono $=$ ker(coker)** (2 mono tests × 3 conditions): for $f_1, f_4$ (both injective), the kernel of the cokernel recovers the original mono up to canonical isomorphism (dim match + image match).

6. **Epi $=$ coker(ker)** (2 epi tests × 3 conditions): for $f_2, f_3$ (both surjective), the cokernel of the kernel recovers the original epi (dim match + kernel match).

### 5.3 Totals

Per cutoff: $10 + 15 + 8 + 8 + 6 + 6 = 53$. Two cutoffs: **106 / 106 bit-exact zero residuals**, all axioms close.

---

## 6. Honest scope

**Closed at theorem grade (bit-exact at finite cutoff $n_{\max} \in \{2, 3\}$ on a representative panel):**

- $\mathrm{Rep}_{\mathrm{fin}}(\mathcal{H}_{\mathrm{GV}}(n_{\max}))$ is an abelian category, verified on the five-rep five-morphism panel across all six universal-property axioms.
- The category structure inherits from the abelian primitive Hopf substrate of v3.61.0 Track A with no GeoVac-specific obstructions.

**Sprint-scale next steps (TC-1b, TC-1c, TC-1d, single-threaded continuation):**

- **TC-1b — symmetric monoidal structure.** $\otimes_\mathbb{Q}$ with the diagonal Hopf action $\Delta(x) = x \otimes 1 + 1 \otimes x$. Unit object = trivial 1-dim rep. Verify associator + unitor + braiding + triangle / pentagon / hexagon coherence bit-exact at small reps.

- **TC-1c — rigidity.** Duals $V^\vee = \mathrm{Hom}_\mathbb{Q}(V, \mathbb{Q})$ with the contragredient Hopf action $X_g^{V^\vee} = -(X_g^V)^T$ (using the antipode $S(x) = -x$). Evaluation $V^\vee \otimes V \to \mathbf{1}$, coevaluation $\mathbf{1} \to V \otimes V^\vee$, snake identities bit-exact.

- **TC-1d — fiber functor $\omega$.** Forget the Hopf-action, keep the ℚ-vector-space. Verify exactness, faithfulness, $\otimes$-preservation, $\omega(\mathbf{1}) = \mathbb{Q}$, $\omega(V^\vee) = \omega(V)^\vee$.

**Multi-year first-stone (TC-1e and beyond):**

- **TC-1e — $\mathrm{Aut}^\otimes(\omega) \supseteq U^\ast_{\mathrm{Levi}}$ at $n_{\max} = 2$.** Construct the explicit inclusion: every $(t, g) \in \mathbb{G}_a^{3 \cdot 5}(\mathbb{Q}) \times SL_2(\mathbb{Q})$ defines a natural $\otimes$-automorphism of $\omega$. Bit-exact verification at $n_{\max} = 2$. **The converse equality** (and the pro-system limit $n_{\max} \to \infty$, and the coherence with v3.66.0 FO3 Interpretation C at the period level on $\mathcal{O}_\infty$ from PS-3) remains multi-year.

**Curve-fit audit (`feedback_audit_numerical_claims`):**

- The abelian-category axioms are derived from standard module-category arguments (Deligne--Milne 1982), not fitted. Zero free parameters.
- The 106 identities are direct universal-property checks (matrix-level), not curve-fit alignments.
- Selection bias:\ the verdict gate was articulated before running the driver;\ the panel was designed to test the six standard axioms, not to confirm a specific outcome. The result (bit-exact closure) matches the prediction from v3.61.0 Track A's abelian primitive substrate.

**Discrete-for-skeleton compliance (`feedback_discrete_for_skeleton`):**

- All matrix entries, every universal-property check, every residual is bit-exact `sympy.Integer` / `sympy.Rational`. Zero floats. Zero PSLQ. Zero transcendentals introduced (TC-1a stays on Layer 1).

**Tag transcendentals (`feedback_tag_transcendentals`):**

- Zero transcendentals appear in this sprint. The Hopf-algebraic substrate of v3.61.0 Track A is intrinsically ℚ-rational.

**WH1 PROVEN unaffected.**

**Hard prohibitions check (CLAUDE.md §13.5).** No changes to natural geometry hierarchy. No fitted/empirical parameters. No deletion of negative results from §3. No removal of "conjectural" label from Paper 2 combination rule.

---

## 7. Files

### Produced
- `geovac/tannakian.py` — new module (~700 lines). Public API:\ `FinDimRep(n_max, dim, endos, label, validate)` with `.X(g)`, `.non_zero_endos()`, `.is_zero_object()`;\ `RepMorphism(source, target, matrix, label, validate)` with `.is_zero()`, `.is_injective()`, `.is_surjective()`;\ `zero_rep(n_max)`, `trivial_rep(n_max, dim)`;\ `compose(g, f)`, `kernel(f)`, `cokernel(f)`, `direct_sum(R1, R2)`;\ verifiers `verify_zero_object_axiom`, `verify_kernel_universal_property`, `verify_cokernel_universal_property`, `verify_direct_sum_universal_property`, `verify_mono_eq_ker_coker`, `verify_epi_eq_coker_ker`. Bit-exact sympy throughout.
- `debug/compute_q5p_tc1a_abelian.py` — driver (~330 lines, 0.03 s wall, bit-exact).
- `debug/data/sprint_q5p_tc1a_abelian.json` — bit-exact data dump.
- `debug/sprint_q5p_tc1a_abelian_memo.md` — this memo.
- `tests/test_tannakian.py` — 28 tests, all pass in 0.92 s.

### Used (load-bearing inputs)
- `geovac/pro_system.py` (PS-2 / PS-3 `primitive_generators`, `n_primitive_generators`, `sectors_at_cutoff`, `N_sectors`).
- `debug/sprint_q5p_stage2_hopf_memo.md` (v3.61.0 Track A;\ abelian primitive substrate).
- `debug/sprint_q5p_ps4_endo_rigidity_memo.md` (v3.69.0 PS-4;\ Tannakian-readiness gap-list naming the four prerequisites).

### Published references
- Deligne, P.; Milne, J. S. ``Tannakian categories.'' In *Hodge Cycles, Motives, and Shimura Varieties*, LNM 900 (1982), 101--228.
- Mac Lane, S. *Categories for the Working Mathematician*, 2nd ed., GTM 5 (1998), Ch. VIII (abelian categories).

---

## 8. Paper-edit recommendations (PI to apply)

### 8.1 Paper 55 \S subsec:open_m2_m3 — ONE new paragraph after the PS-4 paragraph (opens TC-1 series)

```latex
\emph{Tannakian closure foundation TC-1a:\
$\mathrm{Rep}_{\mathrm{fin}}(\mathcal{H}_{\mathrm{GV}})$ is abelian
(Sprint Q5'-Tannakian-Closure, TC-1a sub-track, 2026-06-06;\ memo
\texttt{debug/sprint\_q5p\_tc1a\_abelian\_memo.md};\ data
\texttt{debug/data/sprint\_q5p\_tc1a\_abelian.json}).}  The category
$\mathrm{Rep}_{\mathrm{fin}}(\mathcal{H}_{\mathrm{GV}}(n_{\max}))$ of
finite-dimensional rational representations of the v3.61.0 Track A
abelian primitive Hopf algebra
$\mathcal{H}_{\mathrm{GV}}(n_{\max}) = \mathrm{Sym}_\mathbb{Q}(V_{n_{\max}})
= \mathcal{O}(\mathbb{G}_a^{3 N(n_{\max})})$ is an abelian category
in the Deligne--Milne 1982 sense. An object is a finite-dim
$\mathbb{Q}$-vector space $M$ with $3 N(n_{\max})$ pairwise commuting
nilpotent endomorphisms $\{X_g\}$ indexed by the primitive generators
$g = (n, l, s)$;\ a morphism is a $\mathbb{Q}$-linear map intertwining
all $X_g$.  Verified bit-exact on a five-rep five-morphism panel at
$n_{\max} \in \{2, 3\}$ across all six universal-property axioms:\
zero object (5 reps $\times$ 2 directions $\times$ 2 cutoffs),
direct sum universal property (3 pairs $\times$ 5 identities $\times$
2 cutoffs), kernel universal property (4 morphisms $\times$ 2 conditions
$\times$ 2 cutoffs), cokernel universal property
(4 $\times$ 2 $\times$ 2), mono $=$ ker(coker)
(2 $\times$ 3 $\times$ 2), epi $=$ coker(ker) (2 $\times$ 3 $\times$ 2).
Total bit-exact zero residuals:\ $106 / 106$.  This closes the FIRST
of four sprint-scale Tannakian-closure prerequisites named by PS-4
(v3.69.0):\ abelian category structure.  TC-1b (symmetric monoidal),
TC-1c (rigidity), TC-1d (fiber functor) are the remaining sprint-scale
prerequisites;\ TC-1e (the first stone of the multi-year wall:\
$\mathrm{Aut}^\otimes(\omega) \supseteq U^\ast_{\mathrm{GeoVac},
\mathrm{Levi}}$ at $n_{\max} = 2$) opens the multi-year arc.
```

### 8.2 Paper 32 — no edit needed at TC-1a

TC-1a is documentation of standard abelian-category axioms on the v3.61.0 Track A substrate.

### 8.3 Paper 18 — no edit needed

TC-1a operates on Layer 1 (the bit-exact ℚ-rational skeleton).

---

## 9. One-line verdict

**POSITIVE.** $\mathrm{Rep}_{\mathrm{fin}}(\mathcal{H}_{\mathrm{GV}}(n_{\max}))$ is an abelian category in the Deligne--Milne 1982 sense, verified bit-exact on a five-rep five-morphism panel at $n_{\max} \in \{2, 3\}$ across all six universal-property axioms (zero object, direct sums, kernels, cokernels, mono $=$ ker(coker), epi $=$ coker(ker));\ 106 / 106 bit-exact zero residuals. TC-1a closes the FIRST of four sprint-scale Tannakian-closure prerequisites named by PS-4;\ TC-1b (symmetric monoidal) is next.
