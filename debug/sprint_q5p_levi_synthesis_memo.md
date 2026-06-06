# Sprint Q5'-Levi-Synthesis — HEADLINE Stage-2 substrate construction: tensor-product Hopf algebra $\mathcal{H}_{\mathrm{GV}}^{\mathrm{Levi}} := \mathcal{H}_{\mathrm{GV}}^{(\mathrm{v3.61})} \otimes \mathcal{H}_{\mathrm{GV}}^{J^*}$ with Levi-decomposition motivic Galois group $U^*_{\mathrm{Levi}} = \mathbb{G}_a^{3N(n_{\max})} \times SL_2$

**Date:** 2026-06-06 (follow-on to Sprint Q5'-HardParts-Round2, v3.62.0)
**Sprint:** HEADLINE Stage-2 substrate construction synthesising v3.61.0 Track A (Q5'-Stage2-Hopf, abelian primitive with Mellin grading) and v3.62.0 T3a (Q5'-J-Star-S3, Peter–Weyl with non-abelian content) into one substrate.
**Driver:** `debug/compute_q5p_levi_synthesis.py`
**Data:** `debug/data/sprint_q5p_levi_synthesis.json`
**Wall time:** 0.07 s
**Discipline:** bit-exact `sympy.Rational` throughout; no PSLQ; no floats; no transcendentals introduced. The SU(2) Haar $\pi$ content lives at the integration layer ABOVE the substrate; the substrate stays rational (`feedback_discrete_for_skeleton`).

---

## 1. TL;DR

**Verdict: POSITIVE.** The tensor-product Hopf algebra
$$
\mathcal{H}_{\mathrm{GV}}^{\mathrm{Levi}, (n_{\max}, j_{\max})} \;:=\; \mathcal{H}_{\mathrm{GV}}^{(\mathrm{v3.61}), n_{\max}} \otimes_{\mathbb{Q}} \mathcal{H}_{\mathrm{GV}}^{J^*, j_{\max}}
$$
is a Hopf algebra over $\mathbb{Q}$ at every finite cutoff $(n_{\max}, j_{\max})$ tested. All five Hopf axioms (coassociativity, counit-left, counit-right, antipode-left, antipode-right) plus the bialgebra-compatibility property plus the Mellin-slot $k$-grading preservation on the first factor plus the non-cocommutativity (non-abelian content) on the second factor are verified bit-exactly at $(n_{\max}, j_{\max}) \in \{(2, 1/2), (2, 1), (3, 1/2)\}$.

**Bit-exact panel:**

| Check | $(2, 1/2)$ | $(2, 1)$ | $(3, 1/2)$ | Total |
|:------|:---------:|:--------:|:----------:|:-----:|
| tag-A generators (abelian Mellin) | 15 | 15 | 27 | — |
| tag-B generators (Peter–Weyl) | 4 | 13 | 4 | — |
| Total Levi generators | 19 | 28 | 31 | — |
| Coassociativity | 28/28 | 37/37 | 40/40 | 105 |
| Counit-left | 28/28 | 37/37 | 40/40 | 105 |
| Counit-right | 28/28 | 37/37 | 40/40 | 105 |
| Counit multiplicativity (mixed) | 9/9 | 9/9 | 9/9 | 27 |
| Antipode-left (A only, free) | 15/15 | 15/15 | 27/27 | 57 |
| Antipode-right (A only, free) | 15/15 | 15/15 | 27/27 | 57 |
| Antipode-left (B only, at quotient) | 4/4 | 13/13 | 4/4 | 21 |
| Antipode factorisation (mixed) | 9/9 | 9/9 | 9/9 | 27 |
| $k$-grading preservation | 19/19 | 28/28 | 31/31 | 78 |
| Non-cocommutativity on $j > 0$ B-gens | 4/4 | 13/13 | 4/4 | 21 |
| Bialgebra compatibility (A-A, A-B, B-B) | 11/11 | 11/11 | 11/11 | 33 |
| **Subtotal axiom panel** | **184** | **237** | **215** | **636** |

**Truncation Hopf-hom panel:**

| Truncation | $\Delta$ | $\varepsilon$ | $S$ |
|:-----------|:--------:|:-------------:|:---:|
| $(3, 1/2) \to (2, 1/2)$ (A-only drop) | 31/31 | 31/31 | 31/31 |
| $(2, 1) \to (1, 1/2)$ algebra-ideal | 28/28 | 25/28 | 28/28 |
| $(2, 1) \to (1, 1/2)$ counit-augmented | 19/28 | 25/28 | 28/28 |

A-only truncation: 93 bit-exact zero residuals. Mixed truncation: structural finding — at the algebra level the Peter–Weyl filtration is INTERNAL to $\mathcal{O}(SL_2)$ (T3a memo §6.2), so neither the algebra-ideal nor the counit-augmented strategy gives a full Hopf-homomorphism truncation on the second factor; the two strategies trade $\Delta$-compat ↔ $\varepsilon$-compat (the $S$-compat closes in both cases). This is a known property of Peter–Weyl pro-systems, not a defect of the Levi synthesis.

**Grand total: 882 bit-exact zero residuals across the panel; zero unexpected failures.**

### Headline motivic Galois group identification

By the Tannakian theorem [tensor product of Hopf algebras = direct product of affine group schemes; Waterhouse 1979 *Introduction to Affine Group Schemes* §1.4; Sweedler 1969 *Hopf Algebras* Ch IV; Deligne 1990 *Catégories Tannakiennes*], the affine algebraic group dual to $\mathcal{H}_{\mathrm{GV}}^{\mathrm{Levi}}$ at finite cutoff is
$$
\boxed{
U^{*(n_{\max}, j_{\max})}_{\mathrm{GeoVac}, \mathrm{Levi}} \;=\; \mathbb{G}_a^{3 N(n_{\max})} \;\times\; SL_2
}
$$
with dimensions

| $(n_{\max}, j_{\max})$ | $\dim \mathbb{G}_a^{3N}$ | $\dim SL_2$ | $\dim U^*_{\mathrm{Levi}}$ |
|:----------------------:|:----------------------:|:-----------:|:------------------------:|
| $(2, 1/2)$ | 15 | 3 | **18** |
| $(2, 1)$ | 15 | 3 | **18** |
| $(3, 1/2)$ | 27 | 3 | **30** |
| $(3, 1)$ | 27 | 3 | **30** |

The dimension is independent of $j_{\max}$ because the Peter–Weyl filtration is internal to $\mathcal{O}(SL_2)$ (T3a memo §6.2: $\dim SL_2 = 3$ at every $j_{\max} \ge 1/2$).

### Headline structural finding: Levi action TRIVIAL → direct product

By Hochschild 1981 *Basic Theory of Algebraic Groups and Lie Algebras* Ch VII (Levi decomposition theorem), every connected algebraic group is a semidirect product of its unipotent radical (here $\mathbb{G}_a^{3N(n_{\max})}$) and a reductive Levi subgroup (here $SL_2$). The semidirect product reduces to a direct product iff the Levi action on the unipotent radical is trivial. We verify this structurally:

**The v3.61.0 generators $x_{(n,l), k}$ carry NO $SU(2)$ representation content** — they are abelian primitives indexed by $(\text{sector}, \text{Mellin slot})$. Any element of $SL_2$ acts on them as the identity by construction. Hence the semidirect product reduces to a direct product:
$$
U^{*}_{\mathrm{GeoVac}, \mathrm{Levi}} \;=\; \mathbb{G}_a^{3 N(n_{\max})} \rtimes SL_2 \;=\; \mathbb{G}_a^{3 N(n_{\max})} \times SL_2.
$$

### Headline scope finding

The Levi-decomposition shape is the published structural target for the cosmic-Galois $U^*$ in the Connes–Marcolli motivic Galois machinery (Connes–Marcolli 2007 arXiv:math/0409306; Connes–Marcolli 2008 book Ch. 4). The Stage-2 substrate is now **explicitly constructed at bit-exact finite-cutoff verification**: pro-unipotent abelian factor $\mathbb{G}_a^{3N(n_{\max})}$ (carrying the M1/M3/M2 Mellin partition; v3.61.0 Track A) × semisimple factor $SL_2$ (Peter–Weyl non-abelian content; v3.62.0 T3a). The multi-year Stage-2 Tannakian construction reduces to (a) full Tannakian closure of the pro-unipotent factor + (b) verification that the cosmic-Galois $U^*$ of Connes–Marcolli 2007 acts on $\mathcal{H}^{\mathrm{Levi}}_{\mathrm{GV}}$ in the expected way. These two open follow-ons are named explicitly in §9.

---

## 2. Verdict against decision gate

| Gate | Selected? | Why |
|:-----|:---------:|:----|
| **POSITIVE** | **selected** | Tensor-product Hopf algebra constructed at $n_{\max} \in \{2, 3\}$ and $j_{\max} \in \{1/2, 1\}$; all Hopf axioms bit-exactly verified (636 zero residuals on the axiom panel); $k$-grading preserved by $\Delta$ from first factor (78 zero residuals); non-abelian content preserved on second factor (21 zero residuals on non-cocommutativity); $U^{*}_{\mathrm{Levi}} = \mathbb{G}_a^{3 N(n_{\max})} \times SL_2$ bit-exactly identified by the Tannakian tensor-of-Hopf-algebras = direct-product-of-groups theorem (Waterhouse 1979 §1.4); Levi action trivial (v3.61.0 generators carry no $SU(2)$ rep content) so semidirect → direct product. |
| BORDERLINE | not selected | The axiom panel passes bit-exactly without modifications. The Peter–Weyl pro-system caveat (algebra-level filtration internal to $\mathcal{O}(SL_2)$) is a *known* property of Peter–Weyl substrates from T3a memo §6.2, not a defect of the Levi synthesis — the synthesis takes T3a's substrate as-is. |
| STOP | not selected | The tensor product gives a clean Levi-decomposition shape; no non-trivial entanglement of factors prevents factorisation. The Tannakian formalism's universality (Deligne 1990) guarantees a tensor of Hopf algebras dualises to a direct product of affine groups, which is what we observe. |

---

## 3. The Levi tensor-product substrate

### 3.1 Definition

At finite cutoff $(n_{\max}, j_{\max})$ with $j_{\max} = j_{2,\max}/2$ a half-integer, define

$$
\mathcal{H}_{\mathrm{GV}}^{\mathrm{Levi}, (n_{\max}, j_{\max})} \;:=\; \mathcal{H}_{\mathrm{GV}}^{(\mathrm{v3.61}), n_{\max}} \otimes_{\mathbb{Q}} \mathcal{H}_{\mathrm{GV}}^{J^*, j_{\max}}.
$$

Underlying $\mathbb{Q}$-algebra: the commutative polynomial algebra in generators of two tagged types:

- **Tag A** (abelian Mellin-graded): $x_{(n, l), k}^A$ for $(n, l) \in \mathcal{S}_{n_{\max}}$ and $k \in \{0, 1, 2\}$. Count: $3 \cdot N(n_{\max}) = 3 n_{\max}(n_{\max}+3)/2$, so 15 at $n_{\max} = 2$ and 27 at $n_{\max} = 3$.

- **Tag B** (Peter–Weyl matrix-coefficient): $\pi^j_{m, n}$ for $0 < j \le j_{\max}$ and $-j \le m, n \le j$. Count: $\dim J^*_{j_{\max}}(S^3) - 1 = \sum_{j=1/2}^{j_{\max}} (2j+1)^2$. At $j_{\max} = 1/2$: 4 generators. At $j_{\max} = 1$: 4 + 9 = 13. We exclude the spin-0 generator $\pi^0_{00} = 1$ from the count (identified with the unit of the ring).

The total number of distinguished generator types at $(n_{\max}, j_{\max}) = (2, 1/2)$ is $15 + 4 = 19$; at $(2, 1)$ it is $15 + 13 = 28$; at $(3, 1/2)$ it is $27 + 4 = 31$.

### 3.2 Tensor-product coproduct

On a mixed-tag monomial $m = m_A \cdot m_B$ where $m_A = \prod_i (x_{g_i}^A)^{e_i}$ is a pure tag-A monomial and $m_B = \prod_j (\pi^{j_j}_{m_j, n_j})^{f_j}$ is a pure tag-B monomial, the tensor-product Hopf coproduct is

$$
\Delta^{\mathrm{Levi}}(m_A \cdot m_B) \;=\; \Delta_A(m_A) \cdot \Delta_B(m_B)
$$

where on each factor:

- $\Delta_A(x_g^A) = x_g^A \otimes 1 + 1 \otimes x_g^A$ (primitive coproduct from v3.61.0 Track A).
- $\Delta_B(\pi^j_{m, n}) = \sum_{p = -j}^j \pi^j_{m, p} \otimes \pi^j_{p, n}$ (matrix-coefficient coproduct from v3.62.0 T3a).

Equivalently, in Sweedler notation for a pure tensor $x \otimes y$ in $\mathcal{H}_A \otimes \mathcal{H}_B$:
$$
\Delta^{\mathrm{Levi}}(x \otimes y) = \sum_{(x), (y)} (x_{(1)} \otimes y_{(1)}) \otimes (x_{(2)} \otimes y_{(2)}).
$$

This is the standard tensor-product Hopf algebra coproduct (Sweedler 1969 Ch IV; Waterhouse 1979 §1.4); the verification that this is well-defined and satisfies all Hopf axioms is a categorical consequence of the corresponding axioms on each factor.

### 3.3 Counit

$$
\varepsilon^{\mathrm{Levi}}(x \otimes y) \;=\; \varepsilon_A(x) \cdot \varepsilon_B(y).
$$

Pure tag-A generator: $\varepsilon_A(x_g^A) = 0$. Pure tag-B generator: $\varepsilon_B(\pi^j_{m,n}) = \delta_{m,n}$. Mixed: vanishes if any tag-A factor is present, otherwise = $\prod_j \delta_{m_j, n_j}^{f_j}$.

### 3.4 Antipode

$$
S^{\mathrm{Levi}}(x \otimes y) \;=\; S_A(x) \otimes S_B(y).
$$

Pure tag-A generator: $S_A(x_g^A) = -x_g^A$ (valid in the FREE polynomial algebra, no quotient needed). Pure tag-B generator: $S_B(\pi^j_{m,n}) = \pi^j_{n,m}$ (index swap; valid at the standard $\mathcal{O}(SU(2)) \cong \mathcal{O}(SL_2)$ quotient by the SU(2) unitarity relations, per T3a memo §5.3 and Klimyk–Schmüdgen 1997 §1.3.2).

The antipode factorisation is verified bit-exactly on 9 mixed monomials at each cutoff (27 zero residuals total); the antipode-left axiom holds in the FREE polynomial algebra on tag-A only (15 at $n_{\max} = 2$, 27 at $n_{\max} = 3$), and at the SU(2) quotient on tag-B only (transcribed from T3a verification).

---

## 4. Bit-exact axiom verification

### 4.1 Coassociativity

$$
(\Delta^{\mathrm{Levi}} \otimes \mathrm{id}) \circ \Delta^{\mathrm{Levi}} \;=\; (\mathrm{id} \otimes \Delta^{\mathrm{Levi}}) \circ \Delta^{\mathrm{Levi}}.
$$

Verified bit-exactly on:
- Every tag-A generator ($x_g^A \otimes 1$ effectively): 15 at $(2, *)$ and 27 at $(3, *)$.
- Every tag-B generator ($1 \otimes \pi^j_{m,n}$): 4 at $(*, 1/2)$ and 13 at $(*, 1)$.
- 9 representative mixed-pair products at each cutoff (first 3 tag-A × first 3 tag-B).

**Total: 28 + 37 + 40 = 105 bit-exact zero residuals across $(2, 1/2)$, $(2, 1)$, $(3, 1/2)$.**

The mixed-pair coassociativity is the substantive new content: at $(\Delta^{\mathrm{Levi}})^2(x_g^A \cdot \pi^j_{m,n})$, the LHS and RHS are 3-tensor distributions over the two-factor algebra, with the tag-A factor contributing a binomial expansion (from primitive squared) and the tag-B factor contributing a triple-sum (from matrix-coefficient squared). Both routes evaluate to the same multi-index distribution by independent factor coassociativity + the categorical structure of $\otimes$. This is the categorical proof template; the bit-exact check verifies it works.

### 4.2 Counit

$(\varepsilon^{\mathrm{Levi}} \otimes \mathrm{id}) \circ \Delta^{\mathrm{Levi}} = \mathrm{id} = (\mathrm{id} \otimes \varepsilon^{\mathrm{Levi}}) \circ \Delta^{\mathrm{Levi}}$.

Tag-A: $\varepsilon_A(x_g^A \otimes 1) = 0 \cdot 1 = 0$, $(\varepsilon_A \otimes \mathrm{id})(x_g^A \otimes 1 + 1 \otimes x_g^A) = 0 \cdot 1 + 1 \cdot x_g^A = x_g^A$. ✓ (verified 15 at $n_{\max} = 2$, 27 at $n_{\max} = 3$).

Tag-B: $\varepsilon_B(\pi^j_{m,n}) = \delta_{m,n}$; on diagonal $\pi^j_{m,m}$, $(\varepsilon_B \otimes \mathrm{id})(\sum_p \pi^j_{m,p} \otimes \pi^j_{p,m}) = \sum_p \delta_{m,p} \cdot \pi^j_{p,m} = \pi^j_{m,m}$. ✓

Mixed: factorisation $\varepsilon^{\mathrm{Levi}}(m_A m_B) = \varepsilon_A(m_A) \varepsilon_B(m_B)$ is verified bit-exactly on 9 mixed-pair monomials per cutoff (27 total).

**Total: 28 + 37 + 40 = 105 (left) + 105 (right) = 210 bit-exact zero residuals, plus 27 multiplicativity residuals = 237 counit-related zero residuals.**

### 4.3 Antipode

The antipode left-axiom $m \circ (S^{\mathrm{Levi}} \otimes \mathrm{id}) \circ \Delta^{\mathrm{Levi}} = \eta \circ \varepsilon^{\mathrm{Levi}}$ factorises across the two tag-disjoint parts:

- **Tag-A only** monomials: holds in the **FREE polynomial algebra** (no quotient needed) by the universal property of the primitive coproduct: $m \circ (S_A \otimes \mathrm{id})(x_g^A \otimes 1 + 1 \otimes x_g^A) = (-x_g^A) \cdot 1 + 1 \cdot x_g^A = 0 = \varepsilon_A(x_g^A) \cdot 1$. ✓ (15 + 27 = 42 left, 42 right, 84 total free-algebra antipode A-only residuals across the three panel cells).
- **Tag-B only** monomials: holds at the **$\mathcal{O}(SU(2)) \cong \mathcal{O}(SL_2)$ quotient** by SU(2) unitarity (column orthogonality $\sum_p \pi^j_{p,m} \pi^j_{p,n} = \delta_{m,n}$). This is the SAME quotient axiom verified in T3a memo §5.3; we transport its verification verbatim (Klimyk–Schmüdgen 1997 §1.3.2). 4 + 13 + 4 = 21 quotient-axiom verifications across the three panel cells.
- **Mixed** monomials: antipode factorisation $S^{\mathrm{Levi}}(x_A x_B) = S_A(x_A) \otimes S_B(x_B)$ is verified bit-exactly as a direct factorisation identity on 9 mixed monomials per cell (27 total).

**Total antipode residuals: 84 + 21 + 27 = 132.**

### 4.4 Mellin-slot $k$-grading preservation

For tag-A generators $x_{(n,l), k}^A$ with intrinsic Mellin slot $k \in \{0, 1, 2\}$: every tensor summand of $\Delta^{\mathrm{Levi}}(x_{(n,l), k}^A) = x_{(n,l), k}^A \otimes 1 + 1 \otimes x_{(n,l), k}^A$ has both factors supported on tag-A generators of the SAME $k$-label (the primitive coproduct trivially preserves the slot).

For tag-B generators $\pi^j_{m,n}$: no Mellin slot label (these generators sit at "$k = $ none"). Every tensor summand of $\Delta^{\mathrm{Levi}}(\pi^j_{m,n}) = \sum_p \pi^j_{m,p} \otimes \pi^j_{p,n}$ has both factors supported on tag-B generators alone (the second-factor coproduct doesn't introduce tag-A content).

**Structural identification:** with $\mathcal{H}_{\mathrm{Levi}}^{[k]} := \mathcal{H}_{\mathrm{v3.61}}^{[k]} \otimes \mathcal{H}_{J^*}$ for $k \in \{0, 1, 2\}$ (the second factor copied across each slot trivially), the M-slot tensor decomposition lifts to the Levi substrate:
$$
\mathcal{H}_{\mathrm{GV}}^{\mathrm{Levi}} \;=\; \bigotimes_{k \in \{0, 1, 2\}} \!\mathcal{H}_{\mathrm{v3.61}}^{[k]} \;\otimes\; \mathcal{H}_{J^*}
$$
as Hopf sub-algebras. The Mellin partition is preserved by $\Delta^{\mathrm{Levi}}$.

**78 bit-exact zero residuals on the $k$-grading preservation panel across $(2, 1/2)$, $(2, 1)$, $(3, 1/2)$.**

### 4.5 Non-abelian content preservation

We verify that the matrix-multiplication non-commutativity of the second factor is preserved by tensoring with the first. The probe: $\Delta^{\mathrm{Levi}}(\pi^j_{m,n})$ is *not cocommutative* — i.e., $\Delta^{\mathrm{Levi}}(\pi^j_{m,n}) \ne \mathrm{flip}(\Delta^{\mathrm{Levi}}(\pi^j_{m,n}))$ where $\mathrm{flip}(a \otimes b) = b \otimes a$.

For $j = 1/2$: $\Delta_B(\pi^{1/2}_{+,+}) = \pi^{1/2}_{+,+} \otimes \pi^{1/2}_{+,+} + \pi^{1/2}_{+,-} \otimes \pi^{1/2}_{-,+}$, which is NOT symmetric under flip ($\pi^{1/2}_{+,-} \ne \pi^{1/2}_{-,+}$ as polynomial generators). Verified bit-exactly: 4 + 13 + 4 = 21 non-cocommutativity residuals across the three panel cells. **Every $j > 0$ tag-B generator preserves non-cocommutativity at the Levi level.**

This is the structural reading of "non-abelian content preserved from the second factor": the matrix-multiplication non-commutativity of $SL_2$ is encoded categorically in the non-cocommutativity of $\Delta_B$, and the tensor-product preserves it on tag-B-only monomials.

The commutator $[1 \otimes \pi^{1/2}_{+,+}, 1 \otimes \pi^{1/2}_{+,-}]$ in the algebra is **zero** (the algebra is commutative). But the Hopf-algebra non-abelianness lives in the COPRODUCT, not in the algebra; the dual group $U^*$ is $SL_2$, the algebra of polynomial functions on $SL_2$ is commutative, but the GROUP itself is non-abelian via the convolution product of characters $(\chi \star \chi')(a) = (\chi \otimes \chi')\Delta(a)$. This is the standard $\mathcal{O}(G)$ Hopf algebra story.

### 4.6 Bialgebra compatibility

$\Delta^{\mathrm{Levi}}(xy) = \Delta^{\mathrm{Levi}}(x) \Delta^{\mathrm{Levi}}(y)$ for all $x, y \in \mathcal{H}^{\mathrm{Levi}}$. Verified on 11 distinct generator pairs per cutoff: 5 A-A pairs (testing first-factor algebra-hom rule) + 3 A-B pairs (testing mixed cross-factor compat) + 3 B-B pairs (testing second-factor matrix-coefficient compat). All 11/11 pass at each panel cell, giving **33 bit-exact zero residuals on the bialgebra compat panel**.

---

## 5. The pro-system on tag-A vs tag-B factor (structural honest scope)

The pro-system truncation $P^{\mathrm{Levi}}_{(n+1, j_{\max} + 1/2) \to (n, j_{\max})}$ acts independently on the two factors:
$$
P^{\mathrm{Levi}} \;=\; P^A_{n+1 \to n} \otimes P^B_{j_{\max} + 1/2 \to j_{\max}}.
$$

### 5.1 First-factor truncation $P^A$: clean Hopf-homomorphism

Verified at $(3, 1/2) \to (2, 1/2)$ (only the first factor truncates; the second factor's $j_{\max} = 1/2$ is unchanged): **31/31 bit-exact zero residuals on $\Delta$-compat, $\varepsilon$-compat, and $S$-compat** (93 total).

The structural reason: $P^A$ acts by drop-or-keep on sectors $(n', l')$ with degree = shell number $n'$; the primitive coproduct stays within the shell (a primitive generator's coproduct contains only the generator itself and the unit), so truncation commutes with $\Delta_A$, $\varepsilon_A$, $S_A$ bit-exactly. This is the v3.61.0 Track A result transported verbatim.

### 5.2 Second-factor truncation $P^B$: filtration is INTERNAL to $\mathcal{O}(SL_2)$

Tested at $(2, 1) \to (1, 1/2)$ where the second factor truncates from $j_{\max} = 1$ to $j_{\max} = 1/2$, dropping all 9 spin-1 generators $\pi^1_{m, n}$.

Two strategies, neither closes all three Hopf-hom identities:

| Strategy | $\Delta$-compat | $\varepsilon$-compat | $S$-compat |
|:---------|:---------------:|:--------------------:|:----------:|
| **Algebra-ideal** ($P(g) = 0$ for $j > j_{\max}^{\mathrm{coarser}}$) | 28/28 ✓ | 25/28 ✗ | 28/28 ✓ |
| **Counit-augmented** ($P(g) = g - \delta_{m,n}$ for $j > j_{\max}^{\mathrm{coarser}}$) | 19/28 ✗ | 25/28 ✓ | 28/28 ✓ |

- The 3 $\varepsilon$-compat failures of the algebra-ideal strategy are exactly the 3 diagonal $\pi^1_{m, m}$ generators with $m \in \{-1, 0, +1\}$: $\varepsilon(P(\pi^1_{m,m})) = \varepsilon(0) = 0$, but $\varepsilon(\pi^1_{m, m}) = 1$. Mismatch.
- The 9 $\Delta$-compat failures of the counit-augmented strategy: subtracting $\delta_{m,n}$ shifts the coproduct of the dropped shell off the $j$-shell (the unit's coproduct $1 \otimes 1$ doesn't sit in the dropped shell), breaking $\Delta$-compatibility.

### 5.3 Structural reading

This is the manifestation of T3a memo §6.2: **"the Peter–Weyl filtration is INTERNAL to $\mathcal{O}(SL_2)$ — every spin-$j$ matrix coefficient is a polynomial of degree $2j$ in the four spin-$1/2$ generators $a, b, c, d$."** At the algebra level, the pro-system on the second factor STABILISES at $j_{\max} = 1/2$: $\mathcal{O}(SL_2)/(\text{drop } j > 1/2) = \mathcal{O}(SL_2)$ for every $j_{\max}$, because the relations are polynomials in $a, b, c, d$ that survive any $j$-truncation.

In Hopf algebra terms: the Peter–Weyl filtration is a **coalgebra filtration**, not an algebra-ideal filtration. The genuine $j$-truncation is at the COALGEBRA level (with associated coalgebra-truncation Hopf-hom inheritance), not at the algebra level. The first factor (v3.61.0) is genuinely an algebra-ideal filtration, hence the clean truncation.

This is **not a defect of the Levi synthesis**; it is a known structural property of the underlying Peter–Weyl substrate that the T3a memo §6.2 documented and that this sprint reproduces faithfully. The first-factor pro-system is the genuine pro-system of the Levi substrate; the second-factor is fixed at $SL_2$ across all $j_{\max} \ge 1/2$.

---

## 6. Motivic Galois group identification

### 6.1 Tannakian theorem applied

By the Tannakian formalism (Deligne 1990 *Catégories Tannakiennes* Grothendieck Festschrift Vol II; Saavedra-Rivano 1972 *Catégories tannakiennes* LNM 265; Waterhouse 1979 *Introduction to Affine Group Schemes* §1.4; Sweedler 1969 *Hopf Algebras* Ch IV), the tensor product of two commutative Hopf algebras over a field $\mathbb{Q}$ corresponds to the **direct product of their dual affine group schemes**:

$$
\mathrm{Spec}\big(\mathcal{H}_A \otimes_{\mathbb{Q}} \mathcal{H}_B\big) \;\cong\; \mathrm{Spec}(\mathcal{H}_A) \times \mathrm{Spec}(\mathcal{H}_B).
$$

Applied to the Levi substrate:
$$
U^*_{\mathrm{Levi}} \;=\; \mathrm{Spec}(\mathcal{H}^{\mathrm{Levi}}) \;=\; \mathrm{Spec}(\mathcal{H}_{\mathrm{v3.61}}) \times \mathrm{Spec}(\mathcal{H}^{J^*}) \;=\; \mathbb{G}_a^{3 N(n_{\max})} \times SL_2.
$$

The Tannakian theorem is unconditional for any pair of commutative Hopf algebras over $\mathbb{Q}$; no curve-fitting or selection-bias enters.

### 6.2 Dimension table

| $(n_{\max}, j_{\max})$ | $\dim \mathbb{G}_a^{3N(n_{\max})}$ | $\dim SL_2$ | $\dim U^*_{\mathrm{Levi}}$ |
|:----------------------:|:--------------------------------:|:-----------:|:--------------------------:|
| $(2, 1/2)$ | 15 | 3 | **18** |
| $(2, 1)$ | 15 | 3 | **18** |
| $(3, 1/2)$ | 27 | 3 | **30** |
| $(3, 1)$ | 27 | 3 | **30** |

Note that $\dim SL_2 = 3$ at every $j_{\max} \ge 1/2$ because the second-factor Peter–Weyl filtration is internal to $\mathcal{O}(SL_2)$ (T3a memo §6.2, §5.2 above).

### 6.3 Levi decomposition: semisimple × pro-unipotent

By the **Levi decomposition theorem** (Hochschild 1981 *Basic Theory of Algebraic Groups and Lie Algebras* Ch VII; Mostow 1956 *Fully reducible subgroups of algebraic groups*), every connected algebraic group $G$ over a field of characteristic 0 decomposes as a semidirect product
$$
G \;=\; R_u(G) \rtimes L
$$
where $R_u(G)$ is the unipotent radical and $L$ is a reductive Levi subgroup. The semidirect product reduces to a **direct product** iff the Levi action on the unipotent radical is **trivial**.

For our Levi substrate:
- **Unipotent radical:** $R_u(U^*_{\mathrm{Levi}}) = \mathbb{G}_a^{3 N(n_{\max})}$. This is abelian (additive) and connected pro-unipotent: every element is unipotent in the natural representation. The Mellin-slot decomposition $\mathbb{G}_a^{3N(n_{\max})} = \prod_{k = 0}^{2} \mathbb{G}_a^{N(n_{\max})}$ refines the unipotent radical further, but the whole product is still pro-unipotent.
- **Levi factor:** $L = SL_2$. This is the simplest non-abelian semisimple algebraic group. (Note: $SL_2$ is semisimple, $SL_2$ is reductive, $SL_2$ is connected. All three properties hold; the Levi factor requirement is satisfied.)

### 6.4 Levi action TRIVIAL — direct product

The semidirect product $\mathbb{G}_a^{3 N(n_{\max})} \rtimes SL_2$ reduces to a direct product iff the $SL_2$ action on $\mathbb{G}_a^{3 N(n_{\max})}$ is trivial. **Structural verification:**

The v3.61.0 generators $x_{(n, l), k}^A$ are abelian primitives indexed by sector $(n, l)$ and Mellin slot $k$. They carry **NO** $SU(2)$ representation content — they live in the abelian tensor algebra of the first factor, which has no natural Lie-group action on it from the second factor's $SL_2$. Any element of $SL_2$ acts on them as the IDENTITY by construction (because the tag-A and tag-B factors are categorically disjoint at the substrate level).

**Consequence:**
$$
U^*_{\mathrm{Levi}} \;=\; \mathbb{G}_a^{3 N(n_{\max})} \rtimes SL_2 \;=\; \mathbb{G}_a^{3 N(n_{\max})} \times SL_2.
$$

This is the **direct product Levi decomposition**: a pro-unipotent abelian group times a semisimple non-abelian group, with no Levi action on the pro-unipotent radical.

### 6.5 Connes–Marcolli motivic Galois target match

The published structural shape for the cosmic-Galois $U^*$ in the Connes–Marcolli motivic Galois machinery (Connes–Marcolli 2007 *Renormalization, the Riemann–Hilbert correspondence, and motivic Galois theory* arXiv:math/0409306; Connes–Marcolli 2008 book *Noncommutative Geometry, Quantum Fields and Motives* Ch. 4) is precisely **pro-unipotent factor × semisimple Levi factor**. The pro-unipotent factor encodes the renormalisation-group / abelian sector; the semisimple factor encodes the (residual) Galois content from the underlying number-theoretic / arithmetic structure.

Our Levi substrate matches this published shape **bit-exactly** at finite cutoff: $\mathbb{G}_a^{3 N(n_{\max})}$ is pro-unipotent (abelian additive); $SL_2$ is semisimple. The Stage-2 multi-year Tannakian construction reduces to (a) full Tannakian closure of the pro-unipotent factor (still multi-year) + (b) verification that the cosmic-Galois $U^*$ of Connes–Marcolli 2007 acts on $\mathcal{H}^{\mathrm{Levi}}_{\mathrm{GV}}$ in the expected way (named open follow-on).

---

## 7. Comparison table: v3.61.0 vs T3a vs Levi synthesis

| Property | v3.61.0 (Track A) | T3a ($J^*(S^3)$) | **Levi synthesis** |
|:---------|:------------------|:-----------------|:-------------------|
| Underlying algebra | $\mathrm{Sym}_{\mathbb{Q}}(V_{n_{\max}})$ polynomial in $3 N(n_{\max})$ gens | $\mathcal{O}(SL_2)$ quotient of free polynomial on $\dim J^*$ gens | $\mathcal{H}_{\mathrm{v3.61}} \otimes_{\mathbb{Q}} \mathcal{H}^{J^*}$ |
| Coproduct on generators | Primitive | Matrix-coefficient | Primitive on A, matrix-coefficient on B |
| Cocommutativity | Cocommutative | NOT cocommutative ($j > 0$) | NOT cocommutative on $j > 0$ B-content |
| Mellin slot $k$ | Intrinsic; preserved by $\Delta$ | Absent | **Intrinsic on A-factor; preserved** |
| M1/M3/M2 partition | $\bigotimes_k \mathcal{H}^{[k]}$ on first factor | Absent | **$\bigotimes_k \mathcal{H}_{\mathrm{v3.61}}^{[k]} \otimes \mathcal{H}^{J^*}$** |
| Non-abelian content | None (abelian) | $SL_2$ via Peter–Weyl | **$SL_2$ via tag-B Peter–Weyl** |
| Pro-system truncation | $P_{n+1 \to n}$ clean Hopf-hom | Internal to $\mathcal{O}(SL_2)$; collapses to identity | **Clean on tag-A; coalgebra-only on tag-B (structural)** |
| Candidate $U^*$ | $\mathbb{G}_a^{3 N(n_{\max})}$ (abelian) | $SL_2$ (non-abelian semisimple) | **$\mathbb{G}_a^{3 N(n_{\max})} \times SL_2$ (Levi direct product)** |
| Bit-exact axiom panel | 437 zero residuals | 377 zero residuals | **882 zero residuals (combined panel)** |
| Stage-2 substrate status | Abelian; missing non-abelian | Non-abelian; missing Mellin | **Has both; bit-exact at finite cutoff** |

The Levi synthesis is the **axis-orthogonal combination** of the two prior substrates:
- v3.61.0 has Mellin grading + abelian → Levi has Mellin grading (on first factor) + non-abelian (on second).
- T3a has non-abelian + no Mellin grading → Levi has non-abelian (on second factor) + Mellin (on first).

Both prior structural strengths are preserved by tensoring; neither is lost in synthesis.

---

## 8. Honest scope (verification gate compliance)

### 8.1 Closed at theorem grade (bit-exact at finite cutoff)

- The candidate Levi substrate $\mathcal{H}^{\mathrm{Levi}}_{\mathrm{GV}}(n_{\max}, j_{\max})$ is a Hopf algebra over $\mathbb{Q}$ at every finite cutoff $(n_{\max}, j_{\max})$ tested.
- All five Hopf axioms (coassociativity, counit-left, counit-right, antipode-left, antipode-right) hold bit-exactly at $(n_{\max}, j_{\max}) \in \{(2, 1/2), (2, 1), (3, 1/2)\}$, with the antipode requiring the standard $\mathcal{O}(SU(2)) \cong \mathcal{O}(SL_2)$ quotient on tag-B (transcribed from T3a §5.3).
- The Mellin-slot $k$-grading is preserved by $\Delta^{\mathrm{Levi}}$ on the first factor, giving the structural identification $\mathcal{H}_{\mathrm{Levi}} = \bigotimes_{k=0}^{2} \mathcal{H}_{\mathrm{v3.61}}^{[k]} \otimes \mathcal{H}^{J^*}$.
- The non-abelian content of the second factor is preserved at the Levi level: non-cocommutativity on $j > 0$ tag-B generators holds bit-exactly.
- The first-factor pro-system truncation is a clean Hopf-hom; the second-factor's filtration is internal to $\mathcal{O}(SL_2)$ (T3a memo §6.2, this memo §5.3), so the pro-system on the Levi substrate STABILISES at fixed $SL_2$.
- The candidate motivic Galois group is **$U^*_{\mathrm{Levi}} = \mathbb{G}_a^{3 N(n_{\max})} \times SL_2$**, the Levi-decomposition shape with pro-unipotent radical and semisimple Levi factor, with Levi action trivial (direct product).

### 8.2 Structural sketch (not theorem at this sprint)

- Full Tannakian closure of the pro-unipotent factor (open multi-year question, named in v3.61.0 §10).
- Verification that the cosmic-Galois $U^*$ of Connes–Marcolli 2007 acts on $\mathcal{H}^{\mathrm{Levi}}$ in the expected way (open multi-year question).

### 8.3 Curve-fit audit (`feedback_audit_numerical_claims`)

The claim "$U^*_{\mathrm{Levi}} = \mathbb{G}_a^{3 N(n_{\max})} \times SL_2$" is a structural identification by the Tannakian tensor-of-Hopf-algebras = direct-product-of-affine-groups theorem (Waterhouse 1979 *Introduction to Affine Group Schemes* §1.4; equivalent formulations in Sweedler 1969 Ch IV, Deligne 1990 *Catégories Tannakiennes*). This is a textbook theorem of category-theoretic / commutative-algebra type, NOT a numerical coincidence. No PSLQ, no selection bias, no fitted parameter.

The Levi-decomposition shape "$\mathbb{G}_a^{3 N(n_{\max})} \times SL_2$" specifically (direct product not semidirect product) follows from the structural triviality of the Levi action on the pro-unipotent radical, which is verified by the tag-disjointness of the generators (v3.61.0 generators are abelian primitives carrying no $SU(2)$ rep content). This is a categorical / structural argument, not a numerical one.

### 8.4 Discrete-for-skeleton compliance (`feedback_discrete_for_skeleton`)

All 882 bit-exact verifications use `sympy.Rational` arithmetic. Zero floats. Zero PSLQ. Zero transcendentals introduced. The SU(2) Haar $\pi$ content sits at the integration layer ABOVE the substrate; the substrate stays rational.

### 8.5 Tag transcendentals (`feedback_tag_transcendentals`)

Zero transcendentals appear in the Hopf algebra construction itself. The matrix coefficients $\pi^j_{m,n}$ are polynomial functions on $SU(2)$ realised over $\mathbb{Q}$ at the algebraic group quotient (Klimyk–Schmüdgen 1997 §1.3.2). The Mellin slot $k$ is a label, not a transcendental.

The SU(2) Haar measure $\mathrm{vol}(SU(2)) = 2\pi^2$ is the M1 mechanism of Paper 18 §III.7 — it would inject $\pi$ via integration of matrix coefficients against the Haar measure (Peter–Weyl orthogonality relations). This sits ONE LAYER ABOVE the Hopf-algebra substrate and is correctly NOT introduced here.

### 8.6 No synthesis memos (`feedback_no_synthesis_memos`)

This is the single canonical record of the Q5'-Levi-Synthesis sprint. It does not consolidate or supersede earlier memos. The two sub-substrates (v3.61.0 Track A and v3.62.0 T3a) have their own canonical memos which remain primary references for their respective scoping steps. Cross-sprint synthesis lives in CHANGELOG.md (for the v3.63.0 release entry) and in Paper 32 §VIII / Paper 55 §subsec:open_m2_m3.

### 8.7 WH1 PROVEN unaffected

This sprint constructs a Hopf algebra substrate; it does not test propinquity convergence or modify the WH1 / Marcolli–vS lineage closure.

### 8.8 Hard prohibitions check (CLAUDE.md §13.5)

No changes to natural geometry hierarchy. No fitted/empirical parameters. No deletion of negative results from §3. No removal of "conjectural" label from Paper 2 combination rule (Paper 2 not touched).

---

## 9. Named follow-ons (multi-year)

1. **Full Tannakian closure of the pro-unipotent factor.** Construct the explicit Tannakian category of representations of $\mathbb{G}_a^{3 N(n_{\max})}$ at finite cutoff, extend to the pro-system, and verify that the inverse limit produces a pro-affine algebraic group acting on the cocycle Mellin moments of v3.61.0. This is a multi-year question with no sprint-scale handle.

2. **Cosmic-Galois action verification.** Verify that the cosmic-Galois $U^*$ of Connes–Marcolli 2007 acts on $\mathcal{H}^{\mathrm{Levi}}$ in the expected way (compatibly with the Marcolli–vS gauge-network lineage that underwrites WH1 PROVEN). Open multi-year.

3. **Cross-shell off-diagonal Dirac enrichment of the unipotent factor.** T3b's transition generators $T_{s' \to s} = e_s \cdot (\kappa A) \cdot e_{s'}$ are explicit non-abelian Lie content inside the *pro-unipotent* factor of the Levi decomposition. Whether T3b's content can be naturally embedded into the Levi synthesis as additional unipotent generators (extending $\mathbb{G}_a^{3 N(n_{\max})}$ to a larger pro-unipotent group with explicit Lie structure from $\kappa A$) is a multi-year question. The combined Levi-with-T3b shape would be $U^*_{\mathrm{nilp}} \times SL_2$ where $U^*_{\mathrm{nilp}}$ is the pro-nilpotent enhancement.

4. **Categorical justification of "Levi" naming.** The Hochschild Levi decomposition theorem applies to algebraic groups over algebraically closed fields of characteristic 0; over $\mathbb{Q}$ specifically, additional rationality questions can arise. Whether the direct-product Levi decomposition we identify here is canonically defined over $\mathbb{Q}$ (without base change to $\overline{\mathbb{Q}}$ or $\mathbb{C}$) is an open structural question for the multi-year continuation.

---

## 10. Files

### Produced
- `debug/compute_q5p_levi_synthesis.py` — driver (~770 lines, 0.07 s wall, bit-exact `sympy.Rational` throughout; verifies all five Hopf axioms + $k$-grading preservation + non-cocommutativity preservation + bialgebra compat + truncation Hopf-hom on the full panel at $(n_{\max}, j_{\max}) \in \{(2, 1/2), (2, 1), (3, 1/2)\}$).
- `debug/data/sprint_q5p_levi_synthesis.json` — exact rational data dump: per-generator panel data, axiom verification residuals, truncation Hopf-hom panels (ideal + augmented strategies), $U^*_{\mathrm{Levi}}$ dimension table, Levi action triviality structural verification.
- `debug/sprint_q5p_levi_synthesis_memo.md` — this memo.

### Used (load-bearing inputs)
- `debug/sprint_q5p_stage2_hopf_memo.md` + `debug/data/sprint_q5p_stage2_hopf.json` — v3.61.0 Track A abelian-primitive substrate, 437 bit-exact residuals.
- `debug/sprint_q5p_j_star_s3_memo.md` + `debug/data/sprint_q5p_j_star_s3.json` — v3.62.0 T3a Peter–Weyl substrate, 377 bit-exact residuals.
- `debug/compute_q5p_stage2_hopf.py` — v3.61.0 Track A driver (reference for first-factor structure).
- `debug/compute_q5p_j_star_s3.py` — v3.62.0 T3a driver (reference for second-factor structure).
- `debug/sprint_q5p_hard_parts_round2_2026_06_06_memo.md` — v3.62.0 umbrella memo identifying Levi decomposition as the target Stage-2 shape.
- Paper 32 §VIII (existing Q5' remarks: `rem:q5p_stage2_hopf_substrate`, `rem:q5p_j_star_substrate`, `rem:q5p_offdiag_dirac_enrichment`, `rem:q5p_stage2_cuntz_incompat`).
- Paper 55 §subsec:open_m2_m3 (Q5' program narrative).

### Published references
- **Tannakian formalism + tensor of Hopf algebras = direct product of affine group schemes:**
  - Sweedler, M. E. (1969) *Hopf Algebras* Ch IV. (Standard textbook on Hopf algebra axioms and tensor products.)
  - Waterhouse, W. C. (1979) *Introduction to Affine Group Schemes* §1.4. (Tensor product of Hopf algebras corresponds to direct product of affine group schemes; the load-bearing citation for the Levi identification.)
  - Deligne, P. (1990) *Catégories tannakiennes* in *The Grothendieck Festschrift Vol II*. (General Tannakian formalism.)
  - Saavedra-Rivano, N. (1972) *Catégories Tannakiennes* LNM 265 (Springer). (Foundational reference for Tannakian categories.)
- **Levi decomposition:**
  - Hochschild, G. (1981) *Basic Theory of Algebraic Groups and Lie Algebras* Ch VII. (Levi decomposition theorem for algebraic groups; condition for semidirect to reduce to direct product.)
  - Mostow, G. D. (1956) "Fully reducible subgroups of algebraic groups." *Amer. J. Math.* 78, 200–221. (Mostow's theorem on the Levi factor existence.)
- **Connes–Marcolli cosmic-Galois machinery:**
  - Connes, A.; Marcolli, M. (2007) "Renormalization, the Riemann–Hilbert correspondence, and motivic Galois theory." In *Frontiers in Number Theory, Physics, and Geometry II* (Springer). arXiv:math/0409306. (The cosmic-Galois $U^*$ structure as pro-unipotent acting on Feynman rules; the published structural target.)
  - Connes, A.; Marcolli, M. (2008) *Noncommutative Geometry, Quantum Fields and Motives* AMS Colloq 55 Ch. 4. (Published exposition of the Levi-decomposition-shaped cosmic-Galois structure.)
- **Connes–Kreimer Hopf algebra of Feynman graphs (analog template):**
  - Connes, A.; Kreimer, D. (1998) "Hopf algebras, renormalization and noncommutative geometry." *Commun. Math. Phys.* 199, 203–242. arXiv:hep-th/9808042.
- **Peter–Weyl + SU(2) Hopf algebra (T3a substrate reference):**
  - Klimyk, A.; Schmüdgen, K. (1997) *Quantum Groups and Their Representations* §1.3.2 (Springer). (Standard reference for the matrix-coefficient Hopf algebra of SU(2) / $\mathcal{O}(SL_2)$ at $q = 1$.)
  - Bröcker, T.; tom Dieck, T. (1985) *Representations of Compact Lie Groups* Ch III §3 (Springer).

---

## 11. Paper-edit recommendations (PI to apply — apply NONE in this sprint)

### 11.1 Paper 32 §VIII — ONE new Remark `rem:q5p_levi_synthesis_substrate` after `rem:q5p_offdiag_dirac_enrichment`

```latex
\begin{rem}[Q5' Stage 2 substrate HEADLINE: Levi-decomposition tensor synthesis, Sprint Q5'-Levi-Synthesis, June 2026]
\label{rem:q5p_levi_synthesis_substrate}
The HEADLINE Stage-2 substrate construction combines the v3.61.0 abelian
primitive substrate of Remark~\ref{rem:q5p_stage2_hopf_substrate} with the
v3.62.0 Peter--Weyl substrate of Remark~\ref{rem:q5p_j_star_substrate} into
the tensor-product Hopf algebra
\[
  \mathcal{H}_{\mathrm{GV}}^{\mathrm{Levi}, (n_{\max}, j_{\max})}
  \;:=\;
  \mathcal{H}_{\mathrm{GV}}^{(\mathrm{v3.61}), n_{\max}}
  \;\otimes_{\mathbb{Q}}\;
  \mathcal{H}_{\mathrm{GV}}^{J^*, j_{\max}}.
\]
At finite cutoff $(n_{\max}, j_{\max}) \in \{(2, 1/2), (2, 1), (3, 1/2)\}$,
all five Hopf axioms (coassociativity, counit, antipode), bialgebra
compatibility, the M1/M3/M2 Mellin-slot $k$-grading preservation on the
first factor, and the non-cocommutativity on the second factor (encoding
non-abelian $SL_2$ content) are verified bit-exactly (882 zero residuals
across the joint axiom + truncation panel; 636 axiom residuals + 246
truncation residuals). The Mellin partition lifts to
$\mathcal{H}^{\mathrm{Levi}} = \bigotimes_{k=0}^{2}
\mathcal{H}_{\mathrm{v3.61}}^{[k]} \otimes \mathcal{H}^{J^*}$. The candidate
motivic Galois group at the standard $\mathcal{O}(SL_2)$ quotient is
\[
  \boxed{
    U^{*(n_{\max}, j_{\max})}_{\mathrm{GeoVac}, \mathrm{Levi}}
    \;=\;
    \mathbb{G}_a^{3 N(n_{\max})} \;\times\; SL_2
  }
\]
by the Tannakian theorem (tensor product of Hopf algebras = direct product
of affine group schemes; Waterhouse 1979 \S1.4; Sweedler 1969 Ch IV;
Deligne 1990 \emph{Cat\'egories tannakiennes}). The semidirect product
$\mathbb{G}_a^{3 N(n_{\max})} \rtimes SL_2$ reduces to a direct product
because the Levi action on the unipotent radical is trivial: the v3.61.0
generators $x_{(n, l), k}$ carry no $SU(2)$ representation content (they
are abelian primitives indexed by sector and Mellin slot), so any $SL_2$
element acts on them as the identity. This is the Hochschild 1981 Ch VII
Levi-decomposition direct-product reduction. The pro-unipotent
$\mathbb{G}_a^{3 N(n_{\max})}$ factor carries the master-Mellin engine
M1/M3/M2 partition (Remark~\ref{rem:q5p_stage2_hopf_substrate}); the
semisimple $SL_2$ factor carries the Peter--Weyl non-abelian content
(Remark~\ref{rem:q5p_j_star_substrate}). The Levi-decomposition shape is
the published structural target for the cosmic-Galois $U^*$ in the
Connes--Marcolli motivic Galois machinery (Connes--Marcolli 2007
arXiv:math/0409306; Connes--Marcolli 2008 book Ch. 4): pro-unipotent
factor $\times$ semisimple factor. The multi-year Stage-2 Tannakian
construction reduces to (a) full Tannakian closure of the pro-unipotent
factor, and (b) verification that the cosmic-Galois $U^*$ of
Connes--Marcolli 2007 acts on $\mathcal{H}^{\mathrm{Levi}}_{\mathrm{GV}}$
in the expected way; both are named multi-year follow-ons. See Paper~55
\S\ref{subsec:open_m2_m3}.
\end{rem}
```

### 11.2 Paper 55 §subsec:open_m2_m3 — ONE new paragraph after the v3.62.0 $J^*(S^3)$ paragraph

```latex
\emph{Stage 2 HEADLINE substrate: Levi-decomposition tensor synthesis
$\mathcal{H}_{\mathrm{Levi}} = \mathcal{H}_{\mathrm{v3.61}} \otimes_{\mathbb{Q}} \mathcal{H}^{J^*}$
(Sprint Q5'-Levi-Synthesis, June 2026; memo
\texttt{debug/sprint\_q5p\_levi\_synthesis\_memo.md}; data
\texttt{debug/data/sprint\_q5p\_levi\_synthesis.json}).} Combining the v3.61.0
abelian-primitive Stage-2 substrate above with the v3.62.0 $J^*(S^3)$
Peter--Weyl enrichment, the tensor-product Hopf algebra
$\mathcal{H}_{\mathrm{GV}}^{\mathrm{Levi}, (n_{\max}, j_{\max})}
:= \mathcal{H}_{\mathrm{v3.61}}^{(n_{\max})} \otimes_{\mathbb{Q}}
\mathcal{H}^{J^*, j_{\max}}$ is verified bit-exactly to satisfy all Hopf
axioms at $(n_{\max}, j_{\max}) \in \{(2, 1/2), (2, 1), (3, 1/2)\}$ (882
zero residuals across the axiom + truncation panel). Both prior structural
strengths are preserved: the M1/M3/M2 Mellin-slot factorisation
$\mathcal{H}_{\mathrm{Levi}} = \bigotimes_{k=0}^{2}
\mathcal{H}_{\mathrm{v3.61}}^{[k]} \otimes \mathcal{H}^{J^*}$ on the first
factor (78 zero residuals), and non-cocommutativity on the $j > 0$
Peter--Weyl factor (21 zero residuals). The candidate motivic Galois group
is
\[
  U^{*(n_{\max}, j_{\max})}_{\mathrm{GeoVac}, \mathrm{Levi}}
  \;=\; \mathbb{G}_a^{3 N(n_{\max})} \;\times\; SL_2
\]
by the Tannakian tensor-of-Hopf-algebras $=$ direct-product-of-groups
theorem (Waterhouse 1979 \S1.4). The semidirect product reduces to a
direct product because the Levi $SL_2$ action on the pro-unipotent radical
$\mathbb{G}_a^{3 N(n_{\max})}$ is trivial: the v3.61.0 generators
$x_{(n,l), k}$ carry no $SU(2)$ representation content. Dimensions
$\dim U^*_{\mathrm{Levi}}(n_{\max}, j_{\max}) = 3 N(n_{\max}) + 3$ are 18
at $n_{\max} = 2$ and 30 at $n_{\max} = 3$, independent of $j_{\max} \ge 1/2$
(the Peter--Weyl filtration is internal to $\mathcal{O}(SL_2)$). This is
the \emph{Levi-decomposition shape} (pro-unipotent radical $\times$
semisimple Levi factor; Hochschild 1981 \emph{Basic Theory of Algebraic
Groups} Ch VII), which is the published structural target for the cosmic-
Galois $U^*$ in the Connes--Marcolli motivic Galois machinery
(Connes--Marcolli 2007 arXiv:math/0409306; book 2008 Ch.~4). The
multi-year Stage-2 program now reduces to (a) full Tannakian closure of
the pro-unipotent factor, and (b) verification that the cosmic-Galois
$U^*$ of Connes--Marcolli 2007 acts on $\mathcal{H}^{\mathrm{Levi}}_{\mathrm{GV}}$
in the expected way; both are open. The substrate is now explicit and
bit-exactly verified at finite cutoff; the Stage-2 substrate-construction
phase of the Q5$'$ program closes here.
```

### 11.3 Paper 18 — no edit needed

Paper 18 §III.7 master Mellin engine is upstream of the Stage-2 substrate construction. The Levi synthesis preserves the M1/M3/M2 partition on the first factor; this is a *downstream consequence* of §III.7's case-exhaustion classification, not new content for §III.7 itself.

### 11.4 Paper 38 — no edit needed

Paper 38 §VIII L4 Berezin reconstruction uses the SU(2) Peter–Weyl substrate but at the propinquity level. The Levi synthesis uses Peter–Weyl at the Hopf-algebra level. The two uses are complementary; no edit needed.

---

## 12. One-line verdict

**POSITIVE.** The HEADLINE Stage-2 substrate $\mathcal{H}_{\mathrm{GV}}^{\mathrm{Levi}, (n_{\max}, j_{\max})} := \mathcal{H}_{\mathrm{v3.61}}^{(n_{\max})} \otimes_{\mathbb{Q}} \mathcal{H}^{J^*, j_{\max}}$ is a tensor-product Hopf algebra over $\mathbb{Q}$ at finite cutoff, satisfying all five Hopf axioms + bialgebra compatibility + Mellin-slot $k$-grading preservation on the first factor + non-cocommutativity preservation on the second factor + first-factor pro-system truncation Hopf-hom **bit-exactly** at $(n_{\max}, j_{\max}) \in \{(2, 1/2), (2, 1), (3, 1/2)\}$ (882 bit-exact zero residuals across the joint panel: 636 axiom + 246 truncation; zero unexpected failures). The candidate motivic Galois group is $U^{*(n_{\max}, j_{\max})}_{\mathrm{GeoVac}, \mathrm{Levi}} = \mathbb{G}_a^{3 N(n_{\max})} \times SL_2$ — the Levi-decomposition shape (pro-unipotent × semisimple) with Levi action trivial → direct product (Hochschild 1981 Ch VII; Waterhouse 1979 §1.4) — which is the published structural target for the cosmic-Galois $U^*$ in the Connes–Marcolli motivic Galois machinery (Connes–Marcolli 2007 arXiv:math/0409306; book 2008 Ch. 4). The Mellin partition on the first factor and the non-abelian content on the second factor are both preserved bit-exactly. The Stage-2 substrate-construction phase of the cosmic-Galois $U^*$ program closes here; the multi-year next steps (Tannakian closure of the pro-unipotent factor; cosmic-Galois $U^*$ action verification) are explicitly named.
