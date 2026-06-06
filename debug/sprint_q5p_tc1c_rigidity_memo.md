# Sprint Q5'-Tannakian-Closure TC-1c — rigidity on $\mathrm{Rep}_{\mathrm{fin}}(\mathcal{H}_{\mathrm{GV}})$

**Date:** 2026-06-06 (third sub-track of the Tannakian Reconstruction Foundation sprint, single-thread continuation of TC-1a / TC-1b)
**Sprint:** TC-1c of Q5'-Tannakian-Closure (TC-1a closed v3.70.0; TC-1b closed v3.71.0; TC-1d, TC-1e to follow)
**Driver:** `debug/compute_q5p_tc1c_rigidity.py`
**Module:** `geovac/tannakian.py` (TC-1c additions ~360 lines; total ~1700 lines after TC-1c)
**Data:** `debug/data/sprint_q5p_tc1c_rigidity.json`
**Wall time:** 0.02 s
**Discipline:** bit-exact `sympy.Rational` / `sympy.Integer` throughout. No floats. No PSLQ.

---

## 1. TL;DR

**Verdict: POSITIVE.** $\mathrm{Rep}_{\mathrm{fin}}(\mathcal{H}_{\mathrm{GV}}(n_{\max}))$ is a **rigid** symmetric monoidal category. Every finite-dim rep $V$ admits a dual $V^\vee = \mathrm{Hom}_\mathbb{Q}(V, \mathbb{Q})$ with the contragredient Hopf action $X_g^{V^\vee} = -(X_g^V)^T$ (determined by the antipode $S(x_g) = -x_g$ on primitive generators, v3.61.0 Track A), together with evaluation $\mathrm{ev}_V : V^\vee \otimes V \to \mathbf{1}$ and coevaluation $\mathrm{coev}_V : \mathbf{1} \to V \otimes V^\vee$ satisfying both snake identities bit-exact on a four-rep panel at $n_{\max} \in \{2, 3\}$.

**50 / 50 bit-exact zero residuals** (25 per cutoff):

| Axiom | per cutoff | $\times$ 2 cutoffs |
|:------|:----------:|:------------------:|
| Dual action ($X_g^{V^\vee} = -(X_g^V)^T$, nilpotency, commutativity preserved):\ 4 reps | 4 | 8 |
| Evaluation morphism intertwines:\ 4 reps | 4 | 8 |
| Coevaluation morphism intertwines:\ 4 reps | 4 | 8 |
| First snake $(\mathrm{ev}_V \otimes \mathrm{id}_{V^\vee}) \circ (\mathrm{id}_{V^\vee} \otimes \mathrm{coev}_V) = \mathrm{id}_{V^\vee}$:\ 4 reps | 4 | 8 |
| Second snake $(\mathrm{id}_V \otimes \mathrm{ev}_V) \circ (\mathrm{coev}_V \otimes \mathrm{id}_V) = \mathrm{id}_V$:\ 4 reps | 4 | 8 |
| Double dual $(V^\vee)^\vee = V$ (strict equality of endo data):\ 4 reps | 4 | 8 |
| Unit self-dual $\mathbf{1}^\vee = \mathbf{1}$ | 1 | 2 |
| **Total** | **25** | **50** |

TC-1c closes the **third of four sprint-scale Tannakian-closure prerequisites** named by PS-4 (v3.69.0). TC-1d (fiber functor $\omega$) remains as the sprint-scale stone before the multi-year TC-1e first stone (Tannakian closure $\mathrm{Aut}^\otimes(\omega) = U^\ast_{\mathrm{GeoVac}, \mathrm{Levi}}$).

---

## 2. Verdict against decision gate

| Gate | Verdict |
|:-----|:--------|
| POSITIVE | **selected** — bit-exact closure at the full panel:\ dual action, evaluation, coevaluation, both snake identities, double dual and unit self-dual all bit-exact zero residual at $n_{\max} \in \{2, 3\}$. |
| BORDERLINE | not selected. |
| STOP | rejected — no structural obstruction;\ the abelian primitive Hopf substrate (antipode $S(x) = -x$) makes the dual structure canonical, and the negative-transpose preserves nilpotency and pairwise commutativity automatically (anti-involution of the matrix algebra). |

---

## 3. The rigidity structure

### 3.1 Dual rep (contragredient Hopf action)

For each $V = (M, \{X_g^M\})$ the dual rep $V^\vee = (\mathrm{Hom}_\mathbb{Q}(M, \mathbb{Q}), \{X_g^{V^\vee}\})$ has underlying space $\mathbb{Q}^{\dim M}$ in the dual basis $\{\phi_i\}$ with $\phi_i(e_j) = \delta_{ij}$, and endomorphisms
$$X_g^{V^\vee} = -(X_g^V)^T.$$

The minus sign is determined by the antipode of the abelian primitive Hopf algebra $\mathcal{H}_{\mathrm{GV}}(n_{\max})$ (v3.61.0 Track A):\ $S(x_g) = -x_g$ on primitive generators, lifted multiplicatively to the symmetric algebra. The transpose enforces the duality pairing in the dual basis:\ $(X_g \phi)(v) := -\phi(X_g v)$ corresponds to $\phi \mapsto -X^T \phi$ in coordinates.

**Nilpotency** transports because $X^T$ is nilpotent iff $X$ is (transpose is an automorphism of the matrix algebra and preserves rank).

**Pairwise commutativity** transports because $[A^T, B^T] = -[A, B]^T$;\ in particular $A^T B^T = (BA)^T$ and $B^T A^T = (AB)^T$, so $[A, B] = 0$ implies $[A^T, B^T] = 0$, hence $[-A^T, -B^T] = 0$.

Realised as `dual_rep(V)` in `geovac/tannakian.py`.

### 3.2 Evaluation morphism

$\mathrm{ev}_V : V^\vee \otimes V \to \mathbf{1}$ is the canonical pairing $\phi_i \otimes e_j \mapsto \delta_{ij}$. In the lex basis of $V^\vee \otimes V$ (position $i \cdot \dim V + j$ for $\phi_i \otimes e_j$), the matrix is
$$\mathrm{ev}_V[0, \;\; i \cdot \dim V + j] = \delta_{ij}.$$

Shape $1 \times \dim(V)^2$. Realised as `evaluation_morphism(V)`.

**Intertwining** with the diagonal Hopf action follows because $\mathbf{1}$ has zero action, so the requirement reduces to
$$\mathrm{ev}_V \cdot X_g^{V^\vee \otimes V} = 0$$
for every primitive generator $g$. Expanding $X_g^{V^\vee \otimes V} = -(X^T) \otimes I + I \otimes X$ on the lex basis,
$$(\mathrm{ev}_V \cdot X_g^{V^\vee \otimes V})[0, i \cdot d + k]
  = -\sum_a X[a, i] \cdot \delta_{ak} + \sum_a \delta_{ai} \cdot X[k, a] \cdot 0 + \ldots$$
collapses (with the explicit Kronecker indices) to a difference $-X[k, i] + X[k, i] = 0$. Verified bit-exact for all 4 reps at $n_{\max} \in \{2, 3\}$.

### 3.3 Coevaluation morphism

$\mathrm{coev}_V : \mathbf{1} \to V \otimes V^\vee$ sends $1 \mapsto \sum_i e_i \otimes \phi_i$. In the lex basis of $V \otimes V^\vee$ (position $i \cdot \dim V + j$),
$$\mathrm{coev}_V[i \cdot \dim V + j, \;\; 0] = \delta_{ij}.$$

Shape $\dim(V)^2 \times 1$. Realised as `coevaluation_morphism(V)`.

**Intertwining** with the diagonal Hopf action follows by the symmetric computation:\ $X_g^{V \otimes V^\vee} \cdot \mathrm{coev}_V = 0$ because $X[k, i] \otimes \delta_{ij} + \delta_{ki} \otimes (-X[j, k])$ summed on $i = j$ telescopes to zero. Verified bit-exact for all 4 reps at $n_{\max} \in \{2, 3\}$.

### 3.4 Snake (zigzag) identities

In our canonical lex bases the associator is the identity matrix and the unitors are identity matrices, so the snake identities reduce to clean matrix identities (no unitor or associator chasing needed).

**Second snake on $V$:**
$$(\mathrm{id}_V \otimes \mathrm{ev}_V) \circ (\mathrm{coev}_V \otimes \mathrm{id}_V) = \mathrm{id}_V.$$
- $\mathrm{coev}_V \otimes \mathrm{id}_V$ is the Kronecker product of $\mathrm{coev}_V$ (size $d^2 \times 1$) with $I_V$ (size $d \times d$):\ shape $d^3 \times d$.
- $\mathrm{id}_V \otimes \mathrm{ev}_V$ is the Kronecker product of $I_V$ with $\mathrm{ev}_V$ (size $1 \times d^2$):\ shape $d \times d^3$.
- Their composite is $d \times d$, and the identity asserts it equals $I_d$.

**First snake on $V^\vee$:** structurally dual:
$$(\mathrm{ev}_V \otimes \mathrm{id}_{V^\vee}) \circ (\mathrm{id}_{V^\vee} \otimes \mathrm{coev}_V) = \mathrm{id}_{V^\vee}.$$

Both bit-exact on T1, T2, J2, J3 at $n_{\max} \in \{2, 3\}$. Realised as `verify_snake_identity_first(V)` and `verify_snake_identity_second(V)`.

### 3.5 Double dual is strict equality, not just iso

Because the antipode is involutive on the abelian primitive Hopf algebra ($S^2 = \mathrm{id}$) and because $(M^T)^T = M$, the double dual has $X_g^{(V^\vee)^\vee} = -(-X_g^T)^T = X_g$ on the nose. Strict equality of the endomorphism dictionaries holds;\ no intermediate canonical isomorphism is needed.

This is a structural simplification on the abelian Hopf substrate that does **not** in general hold in non-commutative Hopf algebras (where $S^2$ can be a non-trivial inner automorphism). For $\mathcal{H}_{\mathrm{GV}}$ at the current Track A substrate the canonical map $V \to (V^\vee)^\vee$ is the identity on the lex basis.

### 3.6 Unit is self-dual

$\mathbf{1}^\vee = \mathbf{1}$ strictly:\ $\dim = 1$, no non-zero endos, so $-(0)^T = 0$. Realised as `verify_unit_self_dual(n_max)`.

---

## 4. Bit-exact panel

### 4.1 Substrate

| Label | Construction |
|:-----:|:-------------|
| $T_1$ | Tensor unit $\mathbf{1}$, $\dim = 1$, all $X_g = 0$ |
| $T_2$ | Trivial 2-dim, all $X_g = 0$ |
| $J_2$ | Jordan 2-dim, $X_{(1, 0, 0)} = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$ |
| $J_3$ | Jordan 3-dim, $X_{(1, 0, 0)} = $ canonical 3x3 nilpotent |

### 4.2 Per-cutoff panel totals

Per cutoff $n_{\max} \in \{2, 3\}$:

1. **Dual action** on 4 reps:\ $X_g^{V^\vee} = -(X_g^V)^T$ matches `dual_rep(V).X(g)` and preserves nilpotency + commutativity. **4 / 4**.
2. **Evaluation intertwines** on 4 reps:\ $\mathrm{ev}_V \cdot X_g^{V^\vee \otimes V} = 0$. **4 / 4**.
3. **Coevaluation intertwines** on 4 reps:\ $X_g^{V \otimes V^\vee} \cdot \mathrm{coev}_V = 0$. **4 / 4**.
4. **First snake** on 4 reps:\ $\dim V \times \dim V$ matrix identity. **4 / 4**.
5. **Second snake** on 4 reps:\ $\dim V \times \dim V$ matrix identity. **4 / 4**.
6. **Double dual** on 4 reps:\ strict equality of endo data. **4 / 4**.
7. **Unit self-dual:** **1 / 1**.

**Per cutoff: 25 / 25 bit-exact.**

**Two cutoffs: 50 / 50 bit-exact zero residuals.**

---

## 5. Honest scope

**Closed at theorem grade (bit-exact at finite cutoff $n_{\max} \in \{2, 3\}$ on a representative panel):**

- $\mathrm{Rep}_{\mathrm{fin}}(\mathcal{H}_{\mathrm{GV}}(n_{\max}))$ is a **rigid** symmetric monoidal category.
- Every finite-dim rep $V$ admits a dual $V^\vee$ with contragredient Hopf action $X_g^{V^\vee} = -(X_g^V)^T$.
- $\mathrm{ev}_V$ and $\mathrm{coev}_V$ are explicit Kronecker-delta morphisms intertwining the diagonal Hopf action.
- Both snake identities hold bit-exact as $\dim V \times \dim V$ matrix identities (with associator and unitors equal to identity in lex bases).
- $(V^\vee)^\vee = V$ as **strict equality** of endo data (because the antipode is involutive).
- $\mathbf{1}^\vee = \mathbf{1}$ on the nose.

**Structural surprises (clean, non-blocking):**

- *Strict equality, not iso.* The double dual on the abelian primitive Hopf substrate is **strict equality**, not merely a canonical isomorphism. This is a structural simplification specific to the abelian Track A Hopf substrate where $S(x_g) = -x_g$ and $S^2 = \mathrm{id}$. For Track B (non-abelian, projected v3.62.0 Levi-decomp $\mathbb{G}_a^{3N} \rtimes SL_2$), $S^2$ will in general be inner-non-trivial and the double dual will only be canonically isomorphic, not strictly equal. The strict-equality reading should NOT be claimed beyond TC-1c.

- *No sign convention pathology.* The negative transpose convention $X^{V^\vee} = -X^T$ (as opposed to $+X^T$) is required for intertwining $\mathrm{ev}_V$ and $\mathrm{coev}_V$ with the diagonal Hopf action;\ flipping the sign breaks the panel bit-exactly at the Jordan reps. The $-1$ is the antipode signature, not a free convention choice.

- *Snake identities reduce to identity-matrix algebra.* Because the associator and unitors are identity matrices in lex bases (TC-1b finding), the snake identities collapse to plain matrix identities — no associator chasing, no unitor cancellation. The Kronecker products give the right shapes ($d^3 \times d$ and $d \times d^3$) and the composite is $d \times d$.

**Sprint-scale next steps (single-threaded continuation):**

- **TC-1d — fiber functor $\omega$.** $\omega: \mathrm{Rep}_{\mathrm{fin}}(\mathcal{H}_{\mathrm{GV}}) \to \mathrm{Vec}_\mathbb{Q}$ forgets the Hopf action keeping the $\mathbb{Q}$-vector space. Verify exactness, faithfulness, $\otimes$-preservation $\omega(M \otimes N) \cong \omega(M) \otimes \omega(N)$, $\omega(\mathbf{1}) = \mathbb{Q}$, $\omega(V^\vee) = \omega(V)^\vee$. Bit-exact at small reps.

**Multi-year first-stone (TC-1e):**

- $\mathrm{Aut}^\otimes(\omega) \supseteq U^\ast_{\mathrm{GeoVac}, \mathrm{Levi}} = \mathbb{G}_a^{3 N(n_{\max})} \rtimes SL_2$ at $n_{\max} = 2$. The converse equality (full Tannakian closure on the chosen Hopf substrate) remains multi-year.

**Numerical observation:**

- The associator being the identity in the canonical lex basis is the same structural property identified at TC-1b (and re-exploited here):\ the snake identities at TC-1c reduce to pure Kronecker-matrix identities, so the rigidity verification is the cheapest panel of the four sprint-scale prerequisites (0.02 s wall, smaller than TC-1b's 0.12 s wall).
- The most expensive single check at TC-1c is `verify_double_dual_iso` on J3 at $n_{\max} = 3$, which involves two transpose operations and an equality check on a $3 \times 3$ matrix dictionary. Approximate wall: < 1 ms.

**Curve-fit audit (`feedback_audit_numerical_claims`):**

- The rigidity structure is **derived** from the abelian primitive antipode $S(x_g) = -x_g$ of v3.61.0 Track A and the standard Kronecker-product encoding of the tensor product, not fitted. Zero free parameters.
- The 50 identities are direct matrix comparisons (negative transposes and Kronecker products of well-defined matrices), not curve-fit alignments.
- Selection bias:\ the panel was constructed (T1, T2, J2, J3) to test every shape relevant at TC-1c — trivial 1-dim, trivial 2-dim (no action but non-trivial shape), 2-dim Jordan (the smallest non-trivial nilpotent), 3-dim Jordan (the first higher-order Jordan block).

**Discrete-for-skeleton compliance (`feedback_discrete_for_skeleton`):**

- Every matrix, every intertwining check, every snake identity is bit-exact `sympy.Integer` / `sympy.Rational`. Zero floats. Zero PSLQ. Zero transcendentals (TC-1c stays on Layer 1).

**Tag transcendentals (`feedback_tag_transcendentals`):**

- Zero transcendentals appear in this sprint.

**WH1 PROVEN unaffected.**

**Hard prohibitions check (CLAUDE.md §13.5) clean.**

---

## 6. Files

### Produced
- `geovac/tannakian.py` — TC-1c additions ~360 lines (total ~1700 lines). Public API additions:\ `dual_rep(V)`, `evaluation_morphism(V)`, `coevaluation_morphism(V)`. Verifiers:\ `verify_dual_action(V)`, `verify_evaluation_intertwines(V)`, `verify_coevaluation_intertwines(V)`, `verify_snake_identity_first(V)`, `verify_snake_identity_second(V)`, `verify_double_dual_iso(V)`, `verify_unit_self_dual(n_max)`.
- `debug/compute_q5p_tc1c_rigidity.py` — driver (~270 lines, 0.02 s wall, bit-exact).
- `debug/data/sprint_q5p_tc1c_rigidity.json` — bit-exact data dump.
- `debug/sprint_q5p_tc1c_rigidity_memo.md` — this memo.
- `tests/test_tannakian_rigidity.py` — 30 new tests for TC-1c (all pass in 0.92 s). Combined with TC-1a + TC-1b tests in `tests/test_tannakian.py` (49 tests), the full tannakian suite is **79 tests, all passing in 1.11 s**.

### Used (load-bearing inputs)
- TC-1a substrate (`geovac/tannakian.py` `FinDimRep`, `RepMorphism`, `compose`).
- TC-1b substrate (`geovac/tannakian.py` `tensor_rep`, `tensor_morphism`, `unit_object`, associator-as-identity, braiding-as-swap).
- v3.61.0 Track A (abelian primitive Hopf $\mathcal{H}_{\mathrm{GV}}$, antipode $S(x_g) = -x_g$).
- v3.69.0 PS-4 (Tannakian-readiness gap-list naming the four prerequisites).

### Published references
- Deligne, P.; Milne, J. S. ``Tannakian categories'' (1982), §1 (rigidity in tensor categories).
- Mac Lane, S. *Categories for the Working Mathematician* (1998), Ch. VII (monoidal categories), Ch. VII.7 (duality and rigidity).

---

## 7. Paper-edit recommendations (PI to apply)

### 7.1 Paper 55 §subsec:open_m2_m3 — ONE new paragraph after the TC-1b paragraph

```latex
\emph{Tannakian closure foundation TC-1c:\ rigidity (Sprint
Q5'-Tannakian-Closure, TC-1c sub-track, 2026-06-06;\ memo
\texttt{debug/sprint\_q5p\_tc1c\_rigidity\_memo.md};\ data
\texttt{debug/data/sprint\_q5p\_tc1c\_rigidity.json}).}  The category
$\mathrm{Rep}_{\mathrm{fin}}(\mathcal{H}_{\mathrm{GV}}(n_{\max}))$ is
\emph{rigid}.  Every finite-dim rep $V$ admits a dual
$V^\vee = \mathrm{Hom}_\mathbb{Q}(V, \mathbb{Q})$ with the
contragredient Hopf action
$X_g^{V^\vee} = -(X_g^V)^T$ determined by the antipode
$S(x_g) = -x_g$ of the abelian primitive Hopf algebra (v3.61.0
Track A).  The evaluation
$\mathrm{ev}_V : V^\vee \otimes V \to \mathbf{1}$
is the canonical pairing $\phi_i \otimes e_j \mapsto \delta_{ij}$ in the
dual basis, and the coevaluation
$\mathrm{coev}_V : \mathbf{1} \to V \otimes V^\vee$
sends $1 \mapsto \sum_i e_i \otimes \phi_i$;\ both are explicit
Kronecker-delta morphisms that intertwine the diagonal Hopf action.
Both snake (zigzag) identities reduce to matrix identities in the
canonical lex bases (where the associator and unitors are identity
matrices, cf.\ TC-1b) and are verified bit-exact on a four-rep panel
$(T_1, T_2, J_2, J_3)$ at $n_{\max} \in \{2, 3\}$ across seven axiom
panels totalling 50 identities:\ contragredient action (4 reps $\times$
2 cutoffs $= 8$), evaluation intertwining (4 $\times$ 2 $= 8$),
coevaluation intertwining (4 $\times$ 2 $= 8$), first snake
(4 $\times$ 2 $= 8$), second snake (4 $\times$ 2 $= 8$), double dual
$(V^\vee)^\vee = V$ as strict equality of endo data
(4 $\times$ 2 $= 8$), and unit self-dual $\mathbf{1}^\vee = \mathbf{1}$
(1 $\times$ 2 $= 2$).  Total bit-exact zero residuals:\ $50 / 50$.
\emph{Structural surprise (sub-percent open):} the double dual is
\emph{strict equality} of endomorphism data on the Track-A abelian
substrate (because $S^2 = \mathrm{id}$ via $S(x_g) = -x_g$ and
$(M^T)^T = M$);\ this should not be claimed beyond TC-1c, since a
non-abelian extension of the Hopf substrate would introduce a
non-trivial canonical isomorphism in place of strict equality.  This
closes the THIRD of four sprint-scale Tannakian-closure prerequisites
named by PS-4 (v3.69.0).  TC-1d (fiber functor $\omega$) is next;\
TC-1e (multi-year first stone) opens the heart of Tannakian closure.
```

### 7.2 Paper 32 — no edit needed at TC-1c

### 7.3 Paper 18 — no edit needed

---

## 8. One-line verdict

**POSITIVE.** $\mathrm{Rep}_{\mathrm{fin}}(\mathcal{H}_{\mathrm{GV}})$ is rigid:\ contragredient duals exist with explicit evaluation and coevaluation morphisms, and both snake identities hold bit-exact at $n_{\max} \in \{2, 3\}$. 50 / 50 zero residuals. THIRD of four sprint-scale Tannakian-closure prerequisites closed;\ TC-1d (fiber functor) is next. **Structural surprise (clean, non-blocking):** double dual is strict equality, not iso, on the abelian Track-A substrate — a simplification specific to $S(x_g) = -x_g$ and hence $S^2 = \mathrm{id}$.
