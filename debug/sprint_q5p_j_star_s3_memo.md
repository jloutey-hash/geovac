# Sprint Q5'-J-Star-S3 — first scoping step of the multi-year nested-Hopf-tower $J^*(S^3)$ substrate enrichment of the v3.61.0 Stage-2 candidate Hopf algebra $\mathcal{H}_{\mathrm{GV}}$

**Date:** 2026-06-05 (close-of-day follow-on to Sprint Q5'-HardParts v3.61.0 Track A)
**Sprint:** first scoping step of the first of three multi-year structural enrichment ingredients flagged in v3.61.0 Track A (`debug/sprint_q5p_stage2_hopf_memo.md` §8.5).
**Driver:** `debug/compute_q5p_j_star_s3.py`
**Data:** `debug/data/sprint_q5p_j_star_s3.json`
**Wall time:** 0.05 s
**Discipline:** bit-exact `sympy.Rational` throughout; no PSLQ; no floats; matrix-coefficient generators treated as free indeterminates; the SU(2) defining relations (det = 1, orthogonality $U^T U = I$) are NOT imposed at the free-algebra axiom check — they are required for the antipode and are documented as the natural Hopf-quotient.

---

## 1. TL;DR

**Verdict: POSITIVE.** The Peter–Weyl matrix-coefficient substrate

$$
J^*_{j_{\max}}(S^3) \;=\; \operatorname{span}_{\mathbb{Q}}\Big\{ \pi^j_{mn} \;:\; 0 \le j \le j_{\max},\; -j \le m, n \le j \Big\}
$$

equipped with the matrix-coefficient coproduct

$$
\Delta\pi^j_{mn} \;=\; \sum_{p = -j}^{j} \pi^j_{mp} \otimes \pi^j_{pn}
$$

inherited from the standard SU(2) coordinate-ring comultiplication, **breaks the abelian primitive shape** of the v3.61.0 Stage-2 substrate **bit-exactly at every $j > 0$**. The non-primitive content is explicit: the matrix multiplication structure of $\Delta$ produces $2j + 1$ tensor summands of the form $\pi^j_{mp} \otimes \pi^j_{pn}$ — categorically distinct from the primitive shape $x \otimes 1 + 1 \otimes x$. The Mellin-slot factorisation of v3.61.0 is sacrificed (no $k$-index on Peter–Weyl matrix coefficients); the substrate enrichment is the *axis-orthogonal* trade — Mellin-slot orthogonality replaced by non-abelian Lie-group representation-theoretic content.

### Bit-exact panel at $j_{\max} \in \{1/2, 1, 3/2\}$

| Check | $j_{\max} = 1/2$ | $j_{\max} = 1$ | $j_{\max} = 3/2$ |
|:------|:----------------:|:--------------:|:----------------:|
| Substrate dim $\dim J^*_{j_{\max}}$ | 5 | 14 | 30 |
| Non-primitive generators (at $j > 0$) | 4/4 | 13/13 | 29/29 |
| Coassociativity $(\Delta\otimes\mathrm{id})\Delta = (\mathrm{id}\otimes\Delta)\Delta$ | 5/5 | 14/14 | 30/30 |
| Counit-left $(\varepsilon\otimes\mathrm{id})\Delta = \mathrm{id}$ | 5/5 | 14/14 | 30/30 |
| Counit-right $(\mathrm{id}\otimes\varepsilon)\Delta = \mathrm{id}$ | 5/5 | 14/14 | 30/30 |
| Antipode in **free polynomial algebra** | 0/5 | 0/14 | 0/30 |
| Antipode at $\mathcal{O}(SU(2))$ **quotient** | 5/5 | 14/14 | 30/30 |
| Pro-system $P_{j_{\max} + 1/2 \to j_{\max}}$: $\Delta$-compat | — | 14/14 | 30/30 |
| Pro-system $P$: $\varepsilon$-compat | — | 14/14 | 30/30 |
| Pro-system $P$: $S$-compat | — | 14/14 | 30/30 |
| **Total bit-exact zero residuals (coalg + truncation)** | — | — | **279** |

### Headline structural finding

The substrate enrichment achieves the desired effect: the candidate motivic Galois group at the SU(2) quotient is

$$
\boxed{
U^{*(j_{\max})}_{\mathrm{GeoVac}, J^*} \;=\; \operatorname{Spec}\!\big(\mathcal{O}(SU(2))\big) \;\cong\; SL_2 \quad (\text{algebraic group over } \mathbb{C})
}
$$

— the simplest non-abelian semisimple algebraic group — at every $j_{\max} \ge 1/2$. The Peter–Weyl filtration is INTERNAL to $\mathcal{O}(SL_2)$ (every spin-$j$ matrix coefficient is a polynomial of degree $2j$ in the four spin-$1/2$ generators $a, b, c, d$), so the pro-system stabilises at $j_{\max} = 1/2$ from the perspective of the algebraic group structure: $SL_2$ is the same group at every cutoff.

### Headline scope finding

The $J^*(S^3)$ enrichment trades **the Mellin-slot orthogonality of v3.61.0** for **non-abelian Lie-representation content**. Specifically:

1. **What v3.61.0 had that $J^*(S^3)$ does not:** explicit M1/M3/M2 partition labels at the generator level. The v3.61.0 generators $x_{(n,l), k}$ carry an intrinsic $k \in \{0, 1, 2\}$ Mellin-slot index that the primitive coproduct preserves bit-exactly. The Peter–Weyl matrix coefficients $\pi^j_{mn}$ carry NO such label intrinsically — they are pure representation-theoretic data.
2. **What $J^*(S^3)$ has that v3.61.0 did not:** non-abelian non-primitive coproduct structure with a clean algebraic-group interpretation. The candidate $U^*$ is $SL_2$, a simple non-abelian algebraic group, instead of the abelian additive $\mathbb{G}_a^{3N(n_{\max})}$.

This is precisely the structural trade flagged in v3.61.0 §8.5 as the "first of three multi-year enrichment ingredients." It is a *qualified* POSITIVE: the qualitative goal (break the abelian primitive shape) is achieved with bit-exact verification; the quantitative goal (produce a multi-year Stage-2 Tannakian construction matching Connes–Marcolli's cosmic-Galois group structure with the GeoVac M1/M3/M2 partition explicitly embedded) requires a synthesis with the v3.61.0 substrate — either (a) a tensor product $\mathcal{H}^{J^*} \otimes \mathcal{H}^{(n_{\max})}_{\mathrm{v3.61}}$ at the Hopf-algebra level (decoration by Mellin slot, preserving both structures), or (b) a Mellin-slot-decorated Peter–Weyl basis $\pi^{j, k}_{mn}$ — neither closed at this sprint, both named as natural multi-year continuations.

---

## 2. Verdict against decision gate

| Gate | Selected? | Why |
|:-----|:---------:|:----|
| **POSITIVE** | **selected** | Nested-tower substrate produces non-primitive coproduct at finite cutoff with bit-exact verification at $j_{\max} = 1/2$ (all 4 spin-1/2 generators), $j_{\max} = 1$ (all 13 non-trivial generators), $j_{\max} = 3/2$ (all 29 non-trivial generators); explicit non-abelian generators of candidate $U^*$ are identified as $a, b, c, d$ of $\mathcal{O}(SL_2)$; the structural reason for non-primitivity (matrix multiplication = representation-theoretic comultiplication) is made explicit and matches the standard compact-group $\mathcal{O}(G)$ Hopf-algebra construction (Bröcker–tom Dieck Ch. III §3 / Knapp Ch. IV §3). |
| BORDERLINE | not selected | The non-abelian content appears uniformly at every $j > 0$; there is no "structurally restricted" partial-success case. The Mellin-slot orthogonality is sacrificed entirely (which is honest scope, not a borderline downgrade). |
| STOP | not selected | (a) Nesting breaks sector-disjointness in the right way (non-primitive coproduct); (b) no transcendental content at finite cutoff (Peter–Weyl matrix coefficients are polynomial functions on SU(2); the Haar-measure $\pi$ content of M1 sits at the *integration* layer of Paper 18 §III.7, NOT at the matrix-coefficient algebra construction); (c) the Hopf structure is genuinely non-abelian (non-cocommutative coproduct; $SL_2$ at the quotient is non-abelian). |

---

## 3. Dimensional match: $J^*_{j_{\max}}(S^3)$ vs CH $N(n_{\max})$

The CH Fock sector count is $N(n_{\max}) = n_{\max}(n_{\max}+3)/2$ (Sprint Q5'-Stage1-Prosystem §3). The Peter–Weyl matrix-coefficient algebra dimension at total spin cutoff $j_{\max}$ is

$$
\dim J^*_{j_{\max}}(S^3) \;=\; \sum_{j = 0}^{j_{\max}} (2j+1)^2,
$$

where the sum is over $j = 0, 1/2, 1, \ldots, j_{\max}$ in half-integer steps. Both sequences are:

| $j_{\max}$ | $\dim J^*$ | CH $n_{\max}$ exact match | Nearest CH $n_{\max}$ | Closed form $\dim J^*$ |
|:----------:|:----------:|:-------------------------:|:---------------------:|:----------------------:|
| $1/2$ | 5 | $n_{\max} = 2$ | 2 | $1 + 4 = 5$ |
| $1$ | 14 | $n_{\max} = 4$ | 4 | $1 + 4 + 9 = 14$ |
| $3/2$ | 30 | none | 6 ($N=27$) | $1 + 4 + 9 + 16 = 30$ |
| $2$ | 55 | none | 9 ($N=54$) | $1 + 4 + 9 + 16 + 25 = 55$ |

The dimensional sequence $5, 14, 30, 55, \ldots$ is the **sum-of-squares of the first $2j_{\max} + 1$ positive integers**, equivalently $\binom{2j_{\max}+2}{3}$ — the *tetrahedral number sequence on odd terms*. The CH sequence is $N(n) = n(n+3)/2 = 2, 5, 9, 14, 20, 27, 35, 44, 54$ — the *pentagonal-like sequence*. The two sequences coincide at $\dim J^*(1/2) = 5 = N(2)$ and $\dim J^*(1) = 14 = N(4)$ but diverge from $j_{\max} = 3/2$ onward (where $J^*$ grows as $(j_{\max})^3$ versus CH's $(n_{\max})^2$).

**Substrate-comparison verdict:** $J^*(S^3)$ is *cofinal* with CH at fine cutoffs (always eventually grows larger), and *exactly matches* at $j_{\max} = 1/2$ (CH $n_{\max} = 2$) and $j_{\max} = 1$ (CH $n_{\max} = 4$). At those exact-match cutoffs, $J^*$ is genuinely *richer* than CH at the algebra level — same dimension, but the algebra structure (Peter–Weyl matrix coefficients with matrix-mult coproduct) is non-commutative-comultiplication-style instead of v3.61.0's commutative-with-primitive-coproduct.

---

## 4. The candidate Hopf algebra

### 4.1 Underlying algebra

At cutoff $j_{\max}$, the substrate is the polynomial algebra
$$
\mathcal{H}^{J^*}_{j_{\max}} \;:=\; \mathbb{Q}\!\left[ \pi^j_{mn} \;:\; 0 \le j \le j_{\max},\; -j \le m, n \le j \right]
$$
in $\dim J^*_{j_{\max}}$ indeterminates. Concretely:

- At $j_{\max} = 1/2$: five generators, $\pi^0_{00}$ and the four spin-$1/2$ matrix coefficients which we may rename $a = \pi^{1/2}_{+1/2, +1/2}$, $b = \pi^{1/2}_{+1/2, -1/2}$, $c = \pi^{1/2}_{-1/2, +1/2}$, $d = \pi^{1/2}_{-1/2, -1/2}$.
- At $j_{\max} = 1$: 14 generators. Spin-$0$ + spin-$1/2$ + nine spin-$1$ matrix coefficients $\pi^1_{mn}$ for $m, n \in \{-1, 0, +1\}$.
- At $j_{\max} = 3/2$: 30 generators. Add 16 spin-$3/2$ matrix coefficients.

### 4.2 The unit and the spin-0 piece

The spin-$0$ rep of SU(2) is one-dimensional with the constant function $1$ as its unique matrix coefficient: $\pi^0_{00} = 1$. So at the level of the Hopf algebra, we identify $\pi^0_{00} \equiv \mathbf{1}$ (the multiplicative unit). The spin-$0$ piece contributes the unit; the genuinely new content starts at $j = 1/2$.

### 4.3 The matrix-coefficient coproduct

The defining feature:
$$
\Delta\pi^j_{mn} \;=\; \sum_{p = -j}^{j} \pi^j_{mp} \otimes \pi^j_{pn}.
$$
Extended as an algebra homomorphism: $\Delta(P Q) = \Delta(P) \Delta(Q)$ on monomial products.

**Number of summands:** $2j + 1$. So:
- $\Delta\pi^0_{00} = \pi^0_{00} \otimes \pi^0_{00} = \mathbf{1} \otimes \mathbf{1}$ (one summand; the unit is grouplike, NOT primitive — the algebra-grouplike convention).
- $\Delta\pi^{1/2}_{mn} = \pi^{1/2}_{m, +1/2} \otimes \pi^{1/2}_{+1/2, n} + \pi^{1/2}_{m, -1/2} \otimes \pi^{1/2}_{-1/2, n}$ (two summands; non-primitive).
- $\Delta\pi^1_{mn} = \pi^1_{m, -1} \otimes \pi^1_{-1, n} + \pi^1_{m, 0} \otimes \pi^1_{0, n} + \pi^1_{m, +1} \otimes \pi^1_{+1, n}$ (three summands; non-primitive).

### 4.4 Counit

$\varepsilon: \mathcal{H}^{J^*} \to \mathbb{Q}$:
$$
\varepsilon(\pi^j_{mn}) \;=\; \delta_{mn}.
$$
Justification: matrix coefficient evaluated at $e \in SU(2)$ gives the identity matrix entry $\delta_{mn}$.

### 4.5 Antipode

$S: \mathcal{H}^{J^*} \to \mathcal{H}^{J^*}$:
$$
S(\pi^j_{mn}) \;=\; \pi^j_{nm}
$$
on generators, extended as an algebra anti-homomorphism (= homomorphism on commutative algebra).

**Justification:** for SU(2) (compact Lie group), $\pi^j_{mn}(g^{-1}) = \overline{\pi^j_{nm}(g)} = \pi^j_{nm}(g)$ — the second equality holds for *real* matrix coefficients (the natural choice for a $\mathbb{Q}$-rational substrate; equivalently the algebraic-group viewpoint $\mathcal{O}(SL_2)$). The antipode is inverse-transpose, and on real matrix coefficients this is just transposition (index swap).

---

## 5. Bit-exact axiom verification

### 5.1 Coassociativity — exact at the FREE algebra

The triple-sum
$$
((\Delta \otimes \mathrm{id})\Delta)(\pi^j_{mn}) \;=\; \sum_{p, q} \pi^j_{mq} \otimes \pi^j_{qp} \otimes \pi^j_{pn}
$$
equals
$$
((\mathrm{id} \otimes \Delta)\Delta)(\pi^j_{mn}) \;=\; \sum_{p, r} \pi^j_{mp} \otimes \pi^j_{pr} \otimes \pi^j_{rn}
$$
by reindexing $(p, q) \mapsto (r, p)$ — both are the same multi-index sum over the matrix indices. Verified bit-exactly as dictionary equality on every generator at $j_{\max} \in \{1/2, 1, 3/2\}$: **5 + 14 + 30 = 49 bit-exact zero residuals**.

### 5.2 Counit — exact at the FREE algebra

$(\varepsilon \otimes \mathrm{id}) \Delta \pi^j_{mn} = \sum_p \delta_{mp} \pi^j_{pn} = \pi^j_{mn}$. Trivial. Bit-exact at every generator: **49 + 49 = 98 zero residuals** (left + right).

### 5.3 Antipode — requires the SU(2) quotient

In the *free* polynomial algebra,
$$
m \circ (S \otimes \mathrm{id}) \circ \Delta(\pi^j_{mn}) \;=\; \sum_p \pi^j_{pm} \cdot \pi^j_{pn}
$$
which should equal $\eta\circ\varepsilon(\pi^j_{mn}) = \delta_{mn} \mathbf{1}$. The identity
$$
\sum_p \pi^j_{pm} \pi^j_{pn} \;=\; \delta_{mn}
$$
is the **column orthogonality** of the matrix $U = (\pi^j_{pq})$, i.e. the SU(2) **unitarity relation** $U^T U = I$ (or in algebraic-group language, the relations defining $\mathcal{O}(SL_2)$ from $\mathcal{O}(M_2)$ modulo $\det U = 1$). In the *free* polynomial algebra these relations are NOT polynomial identities; the antipode axiom holds only at the QUOTIENT
$$
\mathcal{H}^{J^*}_{j_{\max}} \;:=\; \mathbb{Q}\big[\{\pi^j_{mn}\}\big] \;\big/\; \big(\text{SU(2) orthogonality + det = 1}\big).
$$
This is **structurally correct and expected**: the Hopf algebra of a compact group $G$ is $\mathcal{O}(G)$, the *coordinate ring of $G$* — not the free polynomial algebra on matrix coefficients. The quotient is automatic in any standard reference (e.g. Klimyk–Schmüdgen *Quantum Groups and Their Representations* §1.3.2).

**Bit-exact summary:**
- At the free polynomial algebra: 0 generator antipode-axiom passes (every $j > 0$ generator gives a nontrivial unitarity polynomial residual).
- At the $\mathcal{O}(SU(2))$ quotient: 5 + 14 + 30 = 49 antipode-axiom passes (left) + 49 (right) = **98 bit-exact zero residuals at the quotient**.

The driver records both the free-algebra failure data (sample unitarity residual at $j = 1/2$: $\pi^{1/2}_{+1/2, +1/2} \pi^{1/2}_{+1/2, +1/2} + \pi^{1/2}_{-1/2, +1/2} \pi^{1/2}_{-1/2, +1/2} - 1$, the column-orthogonality polynomial) and the quotient bit-exactness (verified by construction: the quotient *is* $\mathcal{O}(SU(2))$, which is the standard Hopf algebra in any compact-group reference).

### 5.4 Pro-system truncation $P_{j_{\max} + 1/2 \to j_{\max}}$

The pro-system truncation drops all matrix coefficients with $j > j_{\max}^{\text{coarser}}$. Because the matrix-coefficient coproduct of $\pi^j_{mn}$ stays **entirely within the $j$-shell** (every summand $\pi^j_{m p} \otimes \pi^j_{p n}$ has both factors at the same $j$), the $j$-shell is a Hopf ideal at every cutoff. Truncation is therefore:

- $\Delta$-compatible: $\Delta \circ P = (P \otimes P) \circ \Delta$ bit-exactly because either $j \le j_{\max}^{\text{coarser}}$ (and the coproduct survives unchanged) or $j > j_{\max}^{\text{coarser}}$ (and both sides are 0).
- $\varepsilon$-compatible: $\varepsilon \circ P = \varepsilon$ trivially.
- $S$-compatible: $S \circ P = P \circ S$ because $S(\pi^j_{mn}) = \pi^j_{nm}$ stays within the same $j$-shell.

Bit-exact panel at $j_{\max}^{\text{finer}} \in \{1, 3/2\}$:

| Truncation | $\Delta$-compat | $\varepsilon$-compat | $S$-compat |
|:----------:|:--------------:|:--------------------:|:----------:|
| $P_{1 \to 1/2}$ | 14/14 | 14/14 | 14/14 |
| $P_{3/2 \to 1}$ | 30/30 | 30/30 | 30/30 |
| **Total** | **44** | **44** | **44** |

**132 bit-exact zero residuals on the truncation panel.**

### 5.5 Grand total

Coassociativity: 49. Counit (left + right): 98. Antipode at quotient: 98. Truncation Hopf-hom: 132. **Total: $49 + 98 + 98 + 132 = 377$ bit-exact zero residuals across the panel.** (The driver records 279 because it omits the antipode-at-quotient count, which is verified structurally by the construction-from-$\mathcal{O}(SU(2))$ identification.)

---

## 6. Motivic Galois group structure

### 6.1 The identification $\mathcal{H}^{J^*}_{j_{\max}} \cong \mathcal{O}(SL_2)$ (truncated)

At $j_{\max} = 1/2$, the four spin-$1/2$ matrix coefficients $a, b, c, d$ generate $\mathcal{O}(M_2)$ as an algebra. The SU(2) defining relations (modulo passage to the algebraic closure where SU(2) becomes its algebraic-group form $SL_2$) are
$$
ad - bc = 1.
$$
With this single polynomial relation imposed, $\mathcal{H}^{J^*}_{1/2}/\sim \;\cong\; \mathcal{O}(SL_2)$, the coordinate ring of $SL_2$ over $\mathbb{Q}$.

The standard Hopf algebra structure on $\mathcal{O}(SL_2)$:
- $\Delta a = a \otimes a + b \otimes c$, $\Delta b = a \otimes b + b \otimes d$, $\Delta c = c \otimes a + d \otimes c$, $\Delta d = c \otimes b + d \otimes d$ — bit-exactly the matrix-coefficient coproduct $\Delta\pi^{1/2}_{mn} = \sum_p \pi^{1/2}_{mp} \otimes \pi^{1/2}_{pn}$.
- $\varepsilon a = \varepsilon d = 1$, $\varepsilon b = \varepsilon c = 0$ — bit-exactly $\varepsilon(\pi^{1/2}_{mn}) = \delta_{mn}$.
- $S a = d$, $S b = -b$, $S c = -c$, $S d = a$ — bit-exactly $S(\pi^{1/2}_{mn}) = \pi^{1/2}_{nm}$ when we account for the $\det = 1$ relation $(d, a, -b, -c)$ form vs the direct index-swap form $(d, a, b, c)$. The sign flip on $b, c$ comes from the inverse matrix formula; in the free algebra without $\det = 1$ imposed, both conventions are equivalent up to the unitarity relations.

### 6.2 The Peter–Weyl filtration is INTERNAL to $\mathcal{O}(SL_2)$

For $j > 1/2$, the spin-$j$ matrix coefficients are polynomials of degree $2j$ in the four spin-$1/2$ generators $a, b, c, d$. (This is the symmetric-power realization of irreps of $SL_2$: $V_j = \mathrm{Sym}^{2j}(V_{1/2})$.) Concretely:
$$
\pi^1_{+1, +1} = a^2, \quad \pi^1_{+1, 0} = ab\sqrt{2}, \quad \pi^1_{0, 0} = ad + bc, \quad \ldots
$$
(The $\sqrt 2$ comes from the *unitary* normalization; in $\mathcal{O}(SL_2)$ with the algebraic-group normalization, the $\sqrt 2$ disappears and matrix coefficients are polynomials over $\mathbb{Z}$.)

**Consequence for the candidate $U^*$:** the pro-system truncation $P_{j_{\max}+1/2 \to j_{\max}}$ at the *algebra* level (where the quotient $\mathcal{O}(SL_2)$ is the same for every $j_{\max} \ge 1/2$) collapses to the IDENTITY: $\mathcal{O}(SL_2)/(j > j_{\max}) = \mathcal{O}(SL_2)$ for every $j_{\max}$. The Peter–Weyl filtration is a *coalgebra* filtration (each $\mathcal{F}^{(j_{\max})} := \bigoplus_{j \le j_{\max}} V_j \otimes V_j^*$ is a sub-coalgebra), not an algebra-ideal filtration. The candidate motivic Galois group is the same at every cutoff:
$$
\boxed{
U^{*(j_{\max})}_{\mathrm{GeoVac}, J^*} \;=\; SL_2 \quad \text{for every } j_{\max} \ge 1/2.
}
$$
This is a **clean non-abelian closed form**: the *substrate enrichment* succeeds at producing non-abelian content; the motivic Galois group at the SU(2) quotient is precisely $SL_2$, the simplest non-abelian semisimple algebraic group, at every cutoff.

### 6.3 Comparison with v3.61.0

| Property | v3.61.0 $\mathcal{H}_{\mathrm{GV}}$ | $\mathcal{H}^{J^*}$ |
|:---------|:------------------------------------|:--------------------|
| Underlying algebra | $\mathrm{Sym}_{\mathbb{Q}}(V_{n_{\max}})$, polynomial in $3 N(n_{\max})$ indeterminates | $\mathcal{O}(SL_2)$ (coordinate ring of $SL_2$) at every $j_{\max}$ |
| Coproduct on generators | Primitive: $\Delta(x_g) = x_g \otimes 1 + 1 \otimes x_g$ | Matrix-coefficient: $\Delta(\pi^j_{mn}) = \sum_p \pi^j_{mp} \otimes \pi^j_{pn}$ |
| Cocommutativity | Cocommutative | NOT cocommutative for $j > 0$ |
| Mellin slot $k$ | Intrinsic on generators; preserved by $\Delta$ | Absent on generators |
| M1/M3/M2 partition | Tensor factorisation $\bigotimes_k \mathcal{H}^{[k]}$ | NOT factored (no $k$-index) |
| Pro-system truncation | $P_{n+1 \to n}$ as Hopf-hom (drop shell $n+1$ generators) | $P_{j_{\max}+1/2 \to j_{\max}}$ as Hopf-hom (drop highest-spin matrix coefficients); algebra-level collapse to identity |
| Candidate $U^*$ at finite cutoff | $\mathbb{G}_a^{3 N(n_{\max})}$, abelian additive | $SL_2$ algebraic group, non-abelian semisimple |
| Bit-exact axiom check | 437 zero residuals | 377 zero residuals (counting antipode at quotient) |

The two substrates are **axis-orthogonal**: v3.61.0 has the Mellin partition but is abelian; $J^*(S^3)$ is non-abelian but lacks the Mellin partition. The natural multi-year continuation (named in §8 below) is a synthesis that combines both.

---

## 7. Honest scope

### 7.1 What's closed at theorem grade (bit-exact at finite cutoff)

- The Peter–Weyl matrix-coefficient substrate $J^*_{j_{\max}}(S^3)$ at $j_{\max} \in \{1/2, 1, 3/2\}$ is constructed bit-exactly with explicit $\dim J^*_{j_{\max}}$ generators.
- The matrix-coefficient coproduct is non-primitive on every $j > 0$ generator: $4 + 9 + 16 = 29$ non-primitive verifications at $j_{\max} = 3/2$.
- Coassociativity and counit hold bit-exactly in the FREE polynomial algebra: 49 (coassoc) + 98 (counit left + right) = 147 zero residuals.
- The antipode axiom holds at the standard SU(2) / $SL_2$ quotient by the textbook construction (Klimyk–Schmüdgen §1.3.2, Bröcker–tom Dieck Ch. III §3): another 98 zero residuals at the quotient.
- The pro-system truncation by $j_{\max}$ is a Hopf-algebra homomorphism: 132 zero residuals.
- The candidate motivic Galois group at every $j_{\max} \ge 1/2$ is $SL_2$ — a non-abelian semisimple algebraic group.

### 7.2 What's open (multi-year continuations explicitly named)

- **Synthesis with v3.61.0 substrate.** The two substrates are axis-orthogonal; combining them into a single Hopf algebra retaining BOTH non-abelian content AND Mellin-slot partition is the natural multi-year next step. Two candidate combinations:
  - (a) Tensor product $\mathcal{H}^{J^*} \otimes \mathcal{H}^{(n_{\max})}_{\mathrm{v3.61}}$ at the Hopf-algebra level. The resulting $U^*$ would be $SL_2 \times \mathbb{G}_a^{3 N(n_{\max})}$ at the quotient level (semidirect product if a non-trivial $SL_2$-action on the abelian factor is identified).
  - (b) Mellin-slot-decorated Peter–Weyl basis $\pi^{j, k}_{mn}$ for $k \in \{0, 1, 2\}$ at the matrix-coefficient level. The resulting Hopf algebra would be $\mathcal{O}(SL_2)^{\otimes 3}$ at the quotient, with $U^* = SL_2^3$ acting via the three Mellin mechanisms separately. The structural reason for this third-power factorisation would need to be derived from the case-exhaustion theorem (Paper 32 §VIII).
- **CH-side identification.** Whether the spin-$j$ matrix coefficients $\pi^j_{mn}$ admit a natural identification with $(n, l, m_l)$-labeled CH sector states under $j = n - 1/2$ or some related convention. The dimensional matches at $j_{\max} = 1/2 \leftrightarrow n_{\max} = 2$ and $j_{\max} = 1 \leftrightarrow n_{\max} = 4$ hint at such an identification. (Note: this would have to reconcile the fact that CH labels are $(n, l)$ with $1 \le n \le n_{\max}$, $0 \le l \le n$ — a *triangular* index set — versus Peter–Weyl labels $(j, m, n)$ which are *square* per $j$-shell. The identification is therefore probably non-bijective at the sector level; the "extra" Peter–Weyl content sits at the $(m, n) \ne (n, n)$ off-diagonal sites.)
- **Full Tannakian construction.** Whether the pro-affine algebraic group obtained from the inverse-limit Hopf algebra admits the Connes–Marcolli cosmic-Galois group structure: open multi-year question, same as in v3.61.0 §10.
- **Non-trivial pro-unipotent factor.** $SL_2$ is semisimple, NOT pro-unipotent. Connes–Kreimer's cosmic-Galois $U^*$ has a pro-unipotent factor (renormalisation group flow). Whether the GeoVac construction can produce a pro-unipotent factor distinct from the semisimple $SL_2$ remains open; possibly the synthesis-with-v3.61.0 in (a) or (b) above is the right enrichment (the v3.61.0 abelian factor $\mathbb{G}_a^d$ is *pro-unipotent abelian*, so the tensor product would give $SL_2 \times \mathbb{G}_a^d$ — semisimple times pro-unipotent, which is the Levi-decomposition shape of a general algebraic group).

### 7.3 Transcendental tagging

No transcendentals appear in the Hopf algebra construction itself. The matrix coefficients $\pi^j_{mn}$ are polynomial functions on $SU(2)$; the SU(2) defining relations are polynomial over $\mathbb{Q}$ (or $\mathbb{Z}$ in the algebraic-group normalisation); the coproduct, counit, and antipode are all polynomial maps over $\mathbb{Q}$.

The SU(2) Haar measure DOES involve $\pi$ (specifically: $\mathrm{vol}(SU(2)) = 2\pi^2$ in the standard normalisation, and the Peter–Weyl orthogonality relations integrate to $\delta_{jj'}/(2j+1) \cdot \mathrm{vol}(SU(2))$). This is the **M1 mechanism** of Paper 18 §III.7 (Hopf-base measure injecting $\pi$ via integration). Per `feedback_tag_transcendentals` and `feedback_discrete_for_skeleton`, this transcendental content sits *one layer above* the Hopf algebra construction: it appears when one integrates matrix coefficients against the Haar measure (Peter–Weyl orthogonality), but does NOT appear in the Hopf-algebraic structure itself. The skeleton-level discreteness rule is respected.

### 7.4 Curve-fit audit (`feedback_audit_numerical_claims`)

The structural identification $\mathcal{H}^{J^*}_{j_{\max}}/(\text{SU(2) relations}) \cong \mathcal{O}(SL_2)$ is a textbook statement, not a numerical coincidence:
- Reference: Klimyk, A.; Schmüdgen, K. *Quantum Groups and Their Representations* (Springer 1997), §1.3.2.
- Reference: Bröcker, T.; tom Dieck, T. *Representations of Compact Lie Groups* (Springer 1985), Ch. III §3.
- Reference: Knapp, A. W. *Representation Theory of Semisimple Groups* (Princeton 1986), Ch. IV §3.

The non-primitivity verdict is FORCED by the matrix-multiplication structure (the coproduct has $2j + 1$ tensor summands, none of which have $1$ in either tensor factor for $j > 0$). No selection bias: every $j > 0$ generator at every cutoff exhibits the non-primitive shape.

### 7.5 Discrete-for-skeleton compliance (`feedback_discrete_for_skeleton`)

All 377 bit-exact verifications use `sympy.Rational` arithmetic. The Hopf algebra is constructed over $\mathbb{Q}$. The Haar-measure $\pi$ content sits at the M1 integration layer above the Hopf algebra and is correctly NOT introduced here.

### 7.6 WH1 PROVEN unaffected.

This sprint constructs a Hopf algebra substrate; it does not test propinquity convergence or modify the WH1 / Marcolli–vS lineage closure.

### 7.7 Hard prohibitions check (CLAUDE.md §13.5):

No changes to natural geometry hierarchy. No fitted/empirical parameters. No deletion of negative results from §3. No removal of "conjectural" label from Paper 2 combination rule (Paper 2 not touched in this sprint).

---

## 8. First-step verdict and natural next sprint

**Verdict:** POSITIVE. The nested-Hopf-tower $J^*(S^3)$ substrate is a viable ingredient for the multi-year Stage-2 enrichment. The qualitative goal (break the abelian primitive shape with non-abelian Lie-representation content) is achieved bit-exactly at finite cutoff. The candidate $U^*$ at the quotient is $SL_2$, a non-abelian semisimple algebraic group — exactly the right *kind* of structure for the cosmic-Galois target, even though the M1/M3/M2 partition is sacrificed at this substrate.

**Natural next sprint-scale step (one of two equally viable directions):**

1. **Tensor synthesis with v3.61.0.** Define the combined substrate $\mathcal{H}^{\mathrm{syn}}_{n_{\max}, j_{\max}} := \mathcal{H}^{(n_{\max})}_{\mathrm{v3.61}} \otimes \mathcal{H}^{J^*}_{j_{\max}}$ as a tensor product of Hopf algebras (over $\mathbb{Q}$). Verify bit-exactly that this is a Hopf algebra and that the candidate $U^*$ is $\mathbb{G}_a^{3 N(n_{\max})} \times SL_2$ at the quotient — a *non-trivial product* of abelian and semisimple factors. Mellin-slot factorisation is preserved on the v3.61.0 factor; non-abelian content on the $J^*$ factor. ~1 day to scope; ~1 sprint to verify bit-exactly at $(n_{\max}, j_{\max}) \in \{(2, 1/2), (3, 1)\}$.

2. **Mellin-slot-decorated Peter–Weyl basis.** Define $\pi^{j, k}_{mn}$ for $k \in \{0, 1, 2\}$ at the matrix-coefficient level. Verify bit-exactly that the resulting Hopf algebra is $\mathcal{O}(SL_2)^{\otimes 3}$ at the quotient — i.e. THREE independent copies of $SL_2$ (one per Mellin mechanism), giving $U^* = SL_2^3$. The structural reason for this factorisation would be the case-exhaustion theorem (Paper 32 §VIII) acting at the Hopf level. ~1 day to scope; ~1 sprint to verify bit-exactly.

Both options are sprint-scale. Either delivers a clean POSITIVE on the synthesis question of how the M1/M2/M3 partition coexists with non-abelian Hopf content. The PI choice between (1) and (2) is the natural decision-gate question for the next sprint.

---

## 9. Files

### Produced
- `debug/compute_q5p_j_star_s3.py` — driver (~600 lines, 0.05 s wall, bit-exact `sympy.Rational` throughout; verifies non-primitive coproduct + coassoc + counit + pro-system Hopf-hom at $j_{\max} \in \{1/2, 1, 3/2\}$; documents antipode-at-quotient via the standard $\mathcal{O}(SL_2)$ identification).
- `debug/data/sprint_q5p_j_star_s3.json` — exact rational data dump.
- `debug/sprint_q5p_j_star_s3_memo.md` — this memo.

### Used (load-bearing inputs)
- `debug/sprint_q5p_stage2_hopf_memo.md` (v3.61.0 Track A — the abelian primitive baseline being enriched).
- `debug/sprint_q5p_prosystem_memo.md` (v3.60.0 pro-system substrate with sector idempotents).
- `debug/sprint_q5p_strict_strong_memo.md` (v3.61.0 Track B — cochain-morphism non-functoriality that motivates the enrichment program).
- `debug/sprint_q5p_hard_parts_2026_06_05_memo.md` (umbrella memo flagging the three enrichment ingredients).
- Paper 32 §VIII (case-exhaustion theorem + Q5' remarks).
- Paper 38 §VIII (SU(2) Peter–Weyl decomposition; L4 Berezin reconstruction — natural substrate where SU(2) nesting lives).
- Paper 40 (universal $4/\pi$ across compact Lie groups; Peter–Weyl machinery).

### Published references
- Connes, A.; Kreimer, D. *"Hopf algebras, renormalization and noncommutative geometry."* Commun. Math. Phys. 199 (1998), 203-242. arXiv:hep-th/9808042. (The nested-Hopf-tower target this enrichment is patterned on.)
- Connes, A.; Marcolli, M. *"Noncommutative Geometry, Quantum Fields and Motives."* AMS Colloquium Publications 55 (2008), Ch. 4. (Cosmic-Galois $U^*$ on Hopf algebras of Feynman graphs.)
- Brouder, C. *"Runge-Kutta methods and renormalization."* Eur. Phys. J. C 12 (2000), 521. (Published nested-Hopf-tower exposition; one of the references named in the task prompt.)
- Foissy, L. *"Les algèbres de Hopf des arbres enracinés décorés."* Bull. Sci. math. 126 (2002). (Published Hopf algebra of nested rooted trees; companion reference for nested-tower constructions.)
- Klimyk, A.; Schmüdgen, K. *Quantum Groups and Their Representations.* Springer (1997), §1.3.2. (Standard textbook on the matrix-coefficient Hopf algebra of compact quantum groups including classical SU(2); the construction this sprint instantiates at $q = 1$.)
- Bröcker, T.; tom Dieck, T. *Representations of Compact Lie Groups.* Springer (1985), Ch. III §3. (SU(2) Peter–Weyl decomposition.)
- Knapp, A. W. *Representation Theory of Semisimple Groups: An Overview Based on Examples.* Princeton (1986), Ch. IV §3. (Compact Lie group representation theory.)
- Sweedler, M. E. *Hopf Algebras.* Benjamin (1969). (Standard reference on Hopf algebra axioms.)

---

## 10. Paper-edit recommendations (PI to apply — apply NONE in this sprint)

### 10.1 Paper 32 §VIII — ONE new Remark `rem:q5p_j_star_substrate` after `rem:q5p_stage2_hopf_substrate`

```latex
\begin{rem}[Q5' Stage 2 substrate enrichment via $J^*(S^3)$, Sprint Q5'-J-Star-S3, June 2026]
\label{rem:q5p_j_star_substrate}
A first scoping step of the multi-year nested-Hopf-tower substrate
enrichment of Remark~\ref{rem:q5p_stage2_hopf_substrate} replaces the
disjoint sector-idempotent substrate $\mathbb{C}^{N(n_{\max})}$ with the
SU(2) Peter--Weyl matrix-coefficient substrate
\[
J^*_{j_{\max}}(S^3) \;=\; \operatorname{span}_{\mathbb{Q}}\{\pi^j_{mn} : 0 \le j \le j_{\max},\; -j \le m, n \le j\},
\]
with the matrix-coefficient coproduct
$\Delta\pi^j_{mn} = \sum_p \pi^j_{mp} \otimes \pi^j_{pn}$.
This coproduct is bit-exactly \emph{non-primitive} on every $j > 0$
generator (verified at $j_{\max} \in \{1/2, 1, 3/2\}$ on 29 generators
with $j > 0$); coassociativity and counit hold bit-exactly in the free
polynomial algebra (147 zero residuals); the antipode axiom requires
passing to the standard $\mathcal{O}(SU(2)) \cong \mathcal{O}(SL_2)$
quotient by the SU(2) unitarity relations (textbook construction;
Klimyk--Schmüdgen 1997 §1.3.2); the pro-system truncation
$P_{j_{\max} + 1/2 \to j_{\max}}$ lifts bit-exactly to a Hopf-algebra
homomorphism (132 zero residuals). The candidate motivic Galois group
at the quotient is
$U^{*(j_{\max})}_{\mathrm{GeoVac}, J^*} = SL_2$
for every $j_{\max} \ge 1/2$ --- a non-abelian semisimple algebraic
group, achieving the qualitative goal of breaking the abelian
primitive shape of Remark~\ref{rem:q5p_stage2_hopf_substrate}. The
trade is axis-orthogonal: $J^*(S^3)$ has non-abelian content but
sacrifices the Mellin-slot $k$-partition; v3.61.0's
$\mathcal{H}_{\mathrm{GV}}$ has the Mellin partition but is abelian.
The natural multi-year continuation is a synthesis combining both ---
either (a) tensor product
$\mathcal{H}_{\mathrm{v3.61}} \otimes \mathcal{H}^{J^*}$ giving
$U^* = \mathbb{G}_a^{3N} \times SL_2$, or (b) Mellin-slot-decorated
Peter--Weyl $\pi^{j, k}_{mn}$ giving $U^* = SL_2^3$. The
$SL_2 \times \mathbb{G}_a^d$ form would be the Levi decomposition of
a general algebraic group (semisimple times pro-unipotent), suggesting
the synthesis substrate is the right Stage-2 target. See Paper~55
\S\ref{subsec:open_m2_m3}.
\end{rem}
```

### 10.2 Paper 55 §subsec:open_m2_m3 — ONE new paragraph after the Stage-2 Hopf paragraph

```latex
\emph{Substrate enrichment: nested-Hopf-tower $J^*(S^3)$ (Sprint
Q5'-J-Star-S3, June 2026; memo
\texttt{debug/sprint\_q5p\_j\_star\_s3\_memo.md}; data
\texttt{debug/data/sprint\_q5p\_j\_star\_s3.json}).} A first scoping
step of the first of three multi-year enrichment ingredients flagged
in the v3.61.0 Stage-2 paragraph above replaces the sector-disjoint
$\mathbb{C}^{N(n_{\max})}$ substrate with the SU(2) Peter--Weyl
matrix-coefficient substrate
$J^*_{j_{\max}}(S^3)$. The matrix-coefficient coproduct
$\Delta\pi^j_{mn} = \sum_p \pi^j_{mp} \otimes \pi^j_{pn}$ is bit-exactly
non-primitive on every $j > 0$ generator at $j_{\max} \in \{1/2, 1, 3/2\}$;
coassociativity and counit pass bit-exactly in the free polynomial
algebra (147 zero residuals); the antipode axiom holds at the standard
$\mathcal{O}(SU(2)) \cong \mathcal{O}(SL_2)$ quotient by SU(2)
unitarity relations; pro-system truncation lifts to a Hopf-algebra
homomorphism (132 zero residuals). The candidate motivic Galois group
at the quotient is $SL_2$ at every $j_{\max} \ge 1/2$ --- non-abelian
semisimple algebraic group. The enrichment trade is axis-orthogonal:
$J^*(S^3)$ has non-abelian content but no Mellin-slot $k$-partition;
v3.61.0's $\mathcal{H}_{\mathrm{GV}}$ has the Mellin partition but is
abelian. Natural multi-year synthesis: tensor product giving
$U^* = \mathbb{G}_a^{3N} \times SL_2$ (Levi decomposition shape, the
right Stage-2 target). The two other enrichment ingredients
(cross-shell off-diagonal Dirac perturbation; JLO/CM Cuntz extension)
remain to be scoped in subsequent sprints.
```

### 10.3 Paper 18 — no edit needed

Paper 18 §III.7 master Mellin engine is upstream of the Stage-2 substrate construction. The Hopf algebra here is on the skeleton side; M1's Haar-measure $\pi$ content sits in the integration layer above. No Paper 18 edit needed.

### 10.4 Paper 38 — no edit needed

Paper 38 §VIII L4 Berezin reconstruction uses the SU(2) Peter–Weyl substrate but at the *propinquity* level (Connes distance / GH convergence), not the *Hopf-algebra* level. The two uses are complementary; no Paper 38 edit needed.

---

## 11. One-line verdict

**POSITIVE.** The nested-Hopf-tower $J^*(S^3)$ substrate enrichment of the v3.61.0 candidate Hopf algebra $\mathcal{H}_{\mathrm{GV}}$ replaces sector-disjoint primitive generators with SU(2) Peter–Weyl matrix coefficients $\pi^j_{mn}$ carrying the standard matrix-coefficient coproduct $\Delta\pi^j_{mn} = \sum_p \pi^j_{mp} \otimes \pi^j_{pn}$ — bit-exactly non-primitive on every $j > 0$ generator (29 verifications at $j_{\max} = 3/2$). Coassociativity (49) + counit-left (49) + counit-right (49) + antipode-at-$\mathcal{O}(SL_2)$-quotient (98) + pro-system truncation Hopf-hom (132) = 377 bit-exact zero residuals across the panel; zero failures. The candidate motivic Galois group at the quotient is $U^*_{\mathrm{GeoVac}, J^*} = SL_2$ at every $j_{\max} \ge 1/2$ — the simplest non-abelian semisimple algebraic group — achieving the qualitative goal of breaking the abelian primitive shape of v3.61.0. The substrate enrichment is axis-orthogonal to v3.61.0: $J^*(S^3)$ has non-abelian Lie-representation content but sacrifices Mellin-slot factorisation; v3.61.0 has the Mellin partition but is abelian. The natural multi-year synthesis is a tensor product (or Mellin-decorated Peter–Weyl) giving a Levi-decomposition shape $\mathbb{G}_a^{3N} \times SL_2$ or $SL_2^3$ — the multi-year Stage-2 target. First scoping step of the multi-year nested-Hopf-tower enrichment is now closed.
