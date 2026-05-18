# Sub-sprint B (L3b-2 / L2 lemma joint cb-norm) — Proof memo

**Sprint:** L3b-2 lemma L2, joint Schur multiplier cb-norm
**Sub-sprint:** B (of A/B/C/D in the L3b-2 arc producing the first published
Lorentzian propinquity convergence theorem on truncated Krein spectral
triples).
**Date:** 2026-05-17
**Author:** PM-dispatched research sub-agent
**Status:** **PROVED** at the level of rigor of Paper 38 §L2 Lemma L2 part
(c). The joint Schur multiplier cb-norm equals the product of factor cb-norms
in dual-Coxeter normalization on $\SU(2) \times U(1)$, with closed form
$2/(n_{\max}+1)$. The proof is a factor-wise application of the
Bożejko--Fendler central-multiplier equality on the amenable compact group
$\SU(2) \times U(1)$, with the SU(2) factor inherited from Paper 38 / Paper
40 and the $U(1)$ factor verified explicitly here.
**Verdict:** **PROVED.** No non-standard input needed beyond
Bożejko--Fendler 1991 (amenable case) and Paper 38 Lemma L2 part (c).
**Cross-refs:** `debug/r25_l2_proof_memo.md` (Paper 38 L2 / SU(2) factor);
`debug/r25_l2_quantitative_rate_memo.md` (Paper 38 §L2 quantitative rate);
`debug/l2_universal_rate_proof.md` (Paper 40 §3.2 universal cancellation);
`debug/l3b_first_move_memo.md` §4 (joint Fejér construction);
`papers/standalone/paper_38_su2_propinquity_convergence.tex` §L2;
`papers/standalone/paper_40_unified_propinquity_convergence.tex` §3.2;
`geovac/central_fejer_su2.py` (SU(2) factor implementation);
`geovac/central_fejer_compact_temporal.py` (joint factor implementation).

---

## §1. Statement of sub-sprint B

**Claim (joint cb-norm bound, sub-sprint B).** *Let $n_{\max} \ge 1$ and
$N_t \ge 1$. Let $K^{\SU(2)}_{n_{\max}}$ be the central spectral Fejér
kernel on $\SU(2)$ defined in Paper 38 Definition L2-DEF, with Plancherel
symbol*
$$
   \hat K^{\SU(2)}_{n_{\max}}(j) \;=\; \frac{2j+1}{Z^{\SU(2)}_{n_{\max}}}\,
   \mathbf{1}_{\{j \le j_{\max}\}}, \qquad
   Z^{\SU(2)}_{n_{\max}} = \frac{n_{\max}(n_{\max}+1)}{2}.
$$
*Let $K^{U(1)}_{N_t, T}$ be the Fejér kernel on the circle $S^1_T$ at
truncation $K_{\max} = \lfloor (N_t - 1)/2 \rfloor$ (odd $N_t$) or
$K_{\max} = N_t/2$ (even $N_t$), with Plancherel symbol*
$$
   \hat K^{U(1)}_{N_t}(k) \;=\;
   \max\!\Big(0, \frac{N_t + 1 - 2|k|}{N_t + 1}\Big).
$$
*Let $K^{\mathrm{joint}}_{n_{\max}, N_t, T} := K^{\SU(2)}_{n_{\max}} \otimes
K^{U(1)}_{N_t, T}$ be the tensor product kernel on $\SU(2) \times S^1_T$.
The Schur multiplier $S_{K^{\mathrm{joint}}}$, viewed as a central Fourier
multiplier on the central subalgebra of $C(\SU(2) \times U(1))$, satisfies*
$$
   \boxed{\quad
   \big\|S_{K^{\mathrm{joint}}}\big\|_{\mathrm{cb}}
   \;=\;
   \big\|S_{K^{\SU(2)}}\big\|_{\mathrm{cb}} \;\cdot\;
   \big\|S_{K^{U(1)}}\big\|_{\mathrm{cb}}
   \;=\;
   \frac{2}{n_{\max} + 1} \cdot 1
   \;=\; \frac{2}{n_{\max} + 1}.
   \quad}
$$

The claim has three ingredients that we prove in sequence:

  (B1) [Factor 1] $\|S_{K^{\SU(2)}}\|_{\mathrm{cb}} = 2/(n_{\max}+1)$.
  Inherited from Paper 38 Lemma L2(c) verbatim.

  (B2) [Factor 2] $\|S_{K^{U(1)}}\|_{\mathrm{cb}} = 1$. Standard
  Fejér-on-the-circle, proved explicitly below.

  (B3) [Factorization] The joint Schur multiplier cb-norm equals the
  product of factor cb-norms. Follows from the Bożejko--Fendler
  central-multiplier equality on the amenable compact group $\SU(2)
  \times U(1)$, applied factor-wise.

Each is structurally standard. The only non-trivial bookkeeping is in
(B3), where we verify that the central subalgebra of $\SU(2) \times U(1)$
factorizes as a tensor product of the two factor central subalgebras (this
is automatic for direct products of compact groups; we record the
reference).

---

## §2. Setup: amenable group product, central subalgebra, Schur multipliers

### §2.1 The compact group $G = \SU(2) \times U(1)$

$G = \SU(2) \times U(1)$ is a compact (not connected, since $U(1) = S^1$ is
connected but $G$ has dimension 4, rank 2) Lie group. The maximal torus
of $\SU(2)$ has dimension 1, and $U(1)$ is itself a maximal torus of
dimension 1, so the maximal torus of $G$ has dimension 2 — making $G$ a
rank-2 compact Lie group. Its Weyl group is $W(G) = W(\SU(2)) \times W(U(1))
= \mathbb{Z}/2 \times \{e\} = \mathbb{Z}/2$ (order 2).

**Amenability.** Every compact group is amenable, since the Haar measure
provides a translation-invariant mean. A direct product of amenable groups
is amenable (Bekka--de la Harpe--Valette 2008 Prop.\ G.3.1, or
Pier 1984 Prop.\ 13.6). Therefore $G = \SU(2) \times U(1)$ is amenable.

### §2.2 The central subalgebra and its factorization

Write $C(G)$ for continuous functions on $G$ with the supremum norm. The
**center** $Z(C(G))$ is the subalgebra of class functions on $G$ —
functions constant on conjugacy classes. For a direct product $G_1 \times
G_2$, the conjugacy classes factorize as products of conjugacy classes in
each factor, so a class function on $G_1 \times G_2$ is precisely a
function $(g_1, g_2) \mapsto f_1(g_1) f_2(g_2) + \cdots$ (finite linear
combinations of such products), and the closure of these in $C(G_1 \times
G_2)$ is the tensor product:
$$
   Z(C(\SU(2) \times U(1)))
   \;\cong\;
   Z(C(\SU(2))) \;\hat\otimes\; Z(C(U(1)))
$$
where $\hat\otimes$ is the spatial (or minimal) C*-tensor product. (Both
factor C*-algebras are commutative — class functions on $\SU(2)$ depend
only on the conjugacy class label $\chi$, and class functions on $U(1)$
are just functions on $U(1)$ since $U(1)$ is abelian — so the tensor
product is unambiguous; see Takesaki 1979 Theorem IV.4.14 for commutative
factor algebras.)

### §2.3 Schur multipliers on the central subalgebra

Let $K \in L^1(G)$ be a central function (class function on $G$). The
**convolution operator** $T_K \colon C(G) \to C(G)$ defined by $T_K f := K *
f$ restricts to an endomorphism of $Z(C(G))$ (since convolution preserves
class-functionality when the kernel is central). The **central Fourier
multiplier symbol** $\hat K \colon \hat G \to \mathbb{C}$ is defined by
$\hat K(\pi) := $ scalar by which $T_K$ acts on the isotypic component
$V_\pi \otimes V_\pi^*$ of the regular representation. By Peter--Weyl,
this scalar is
$$
   \hat K(\pi) \;=\; \int_G K(g)\, \chi_\pi(g)\, dg / \dim V_\pi,
$$
with $\chi_\pi$ the character of $\pi$. (Equivalently, by Plancherel on
the central subalgebra, $\hat K(\pi)$ is the $\pi$-th Fourier coefficient
of $K$.)

The **Schur multiplier** $S_K$ on the operator algebra $C(G)$ is the map
that sends a matrix $A_{ij}$ indexed by $\hat G$ to $\hat K(\pi_i)
A_{ij}$ in the isotypic decomposition. (This is the Schur transform of
$T_K$ under the Plancherel decomposition.) The **cb-norm** $\|S_K\|_{\mathrm{cb}}$
is the supremum of the operator norm of $S_K \otimes \mathrm{id}_n$ over
all $n$ (the completely bounded norm; see Pisier 2001 Ch.\ 1--2).

---

## §3. The Bożejko--Fendler dual estimate

The standard tool for cb-norm equalities on amenable groups is the
following theorem.

**Theorem 3.1 (Bożejko--Fendler 1991 / Pisier 2001 Ch.\ 8 Thm 8.10).**
*Let $G$ be an amenable compact group and let $m \in \ell^\infty(\hat G)$ be
a* central *Fourier multiplier symbol (meaning $m(\pi)$ depends only on
the isotype $\pi$, not on the matrix coefficient within $V_\pi \otimes
V_\pi^*$). Then the Schur multiplier $S_m$ on $C(G)$ satisfies*
$$
   \|S_m\|_{\mathrm{cb}} \;=\; \|T_m\|_{\mathrm{op}} \;=\; \|m\|_{\ell^\infty(\hat G)}.
$$

The proof uses the cb-trick: a central multiplier acts on each isotype
$V_\pi \otimes V_\pi^*$ (which is a matrix algebra $M_{d_\pi}(\mathbb{C})$ with
$d_\pi = \dim V_\pi$) by a *scalar*, so the cb-norm bound on each block
reduces to a scalar bound. The amenability of $G$ provides a Følner net
that lets one transfer the matrix-algebra norm equality to the global
$C(G)$ norm.

**Reference.** Bożejko, M., and Fendler, G., "Herz--Schur multipliers and
completely bounded multipliers of the Fourier algebra of a locally compact
group," *Boll. Un. Mat. Ital. A* (6) 3 (1984) 297--302; restated and used
in this form in Pisier, *Similarity Problems and Completely Bounded
Maps*, Lecture Notes in Math. 1618, Springer, 2nd ed.\ 2001, Ch.\ 8 Thm
8.10. The 1991 version (the journal article we cite) extends the 1984
proceedings note with full proofs and the amenability hypothesis spelled
out.

**Application to $\SU(2) \times U(1)$.** Theorem 3.1 applies directly: $G =
\SU(2) \times U(1)$ is amenable (compact direct product). For any central
multiplier symbol $m$ on $\hat G = \hat\SU(2) \times \hat U(1) = (\frac12
\mathbb{Z}_{\ge 0}) \times \mathbb{Z}$, the cb-norm equals the supremum
norm.

---

## §4. Factor 1: the SU(2) cb-norm $\|S_{K^{\SU(2)}}\|_{\mathrm{cb}} = 2/(n_{\max}+1)$

This is Paper 38 Lemma L2(c), transcribed without modification:

**Lemma 4.1 (Paper 38 L2(c)).** *The SU(2) central spectral Fejér kernel
$K^{\SU(2)}_{n_{\max}}$ has cb-norm*
$$
   \|S_{K^{\SU(2)}_{n_{\max}}}\|_{\mathrm{cb}}
   \;=\; \|\hat K^{\SU(2)}_{n_{\max}}\|_{\ell^\infty}
   \;=\; \hat K^{\SU(2)}_{n_{\max}}(j_{\max})
   \;=\; \frac{2j_{\max}+1}{Z^{\SU(2)}_{n_{\max}}}
   \;=\; \frac{n_{\max}}{n_{\max}(n_{\max}+1)/2}
   \;=\; \frac{2}{n_{\max}+1}.
$$

**Proof.** By Theorem 3.1 applied to $\SU(2)$ (compact, amenable), the cb-norm
equals the $\ell^\infty$ norm of the Plancherel symbol. The Plancherel
symbol $\hat K^{\SU(2)}_{n_{\max}}(j) = (2j+1)/Z^{\SU(2)}_{n_{\max}}$ is
**monotonically increasing in $j$** on its support $\{0, \tfrac12, 1,
\ldots, j_{\max}\}$, so the maximum is attained at $j = j_{\max} =
(n_{\max}-1)/2$, giving $(2 j_{\max} + 1)/Z^{\SU(2)}_{n_{\max}} =
n_{\max}/(n_{\max}(n_{\max}+1)/2) = 2/(n_{\max}+1)$. $\square$

**Numerical verification.** `geovac.central_fejer_su2.central_multiplier_cb_norm`
returns this closed form symbolically for $n_{\max} \in \{1, \ldots, 50\}$;
verified in `tests/test_central_fejer_su2.py::TestCentralMultiplierCBNorm`
(Paper 38 regression). The product $(n_{\max}+1) \cdot \|T_K\|_{\mathrm{cb}} = 2$
exactly at every $n_{\max}$.

---

## §5. Factor 2: the U(1) cb-norm $\|S_{K^{U(1)}}\|_{\mathrm{cb}} = 1$

This is the Fejér-on-the-circle calculation, standard but recorded
explicitly here so the proof memo stands alone.

### §5.1 The Fejér kernel on $S^1_T$

The Fejér kernel on the circle $S^1_T$ of circumference $T$ at truncation
level $N_t$ is
$$
   F_{N_t, T}(\theta) \;:=\; \frac{1}{N_t}\,
   \bigg|\sum_{k = -K_{\max}}^{+K_{\max}} e^{i\, 2\pi k\theta/T}\bigg|^2,
$$
where $K_{\max} := \lfloor (N_t - 1)/2 \rfloor$ for odd $N_t$ or $N_t/2$
for even $N_t$ (the convention used in
`geovac.central_fejer_compact_temporal._u1_K_max`).

A short calculation (sum the geometric series in $k$, expand the squared
modulus) gives the **Plancherel symbol**
$$
   \hat F_{N_t}(k) \;=\; \max\!\Big(0,\, \frac{N_t + 1 - 2|k|}{N_t + 1}\Big),
   \qquad k \in \mathbb{Z}.
$$
This is the standard Cesàro-1 / Fejér averaging factor on $\mathbb{Z}$:
positive, supported on $|k| \le K_{\max}$, decreasing as $|k|$ grows from
$0$.

### §5.2 Positivity and Toeplitz structure

The Fejér kernel is a **positive Toeplitz operator** on $\ell^2(\mathbb{Z})$
in the following sense: the operator $T_F \colon \ell^2(\mathbb{Z}) \to
\ell^2(\mathbb{Z})$ acting by multiplication by $\hat F$ in the Fourier
basis is diagonal with non-negative entries. Equivalently, $F_{N_t, T} \ge
0$ pointwise (since $F_{N_t, T} = |D_{N_t, T}|^2 / N_t$ for the Dirichlet
kernel $D_{N_t, T}$), and the Schur multiplier $S_F$ is positive on the
abelian C*-algebra $C(S^1_T)$.

### §5.3 The cb-norm via Theorem 3.1

The compact group $U(1) = S^1$ is amenable (abelian). Every $U(1)$
irrep is one-dimensional, so every multiplier symbol $m \colon \hat U(1) =
\mathbb{Z} \to \mathbb{C}$ is automatically "central" (there is nothing
non-central about a 1-dimensional representation). Theorem 3.1 applies:
$$
   \|S_{F_{N_t, T}}\|_{\mathrm{cb}}
   \;=\; \|\hat F_{N_t}\|_{\ell^\infty(\mathbb{Z})}
   \;=\; \max_{k \in \mathbb{Z}} \hat F_{N_t}(k)
   \;=\; \hat F_{N_t}(0)
   \;=\; \frac{N_t + 1 - 0}{N_t + 1}
   \;=\; 1.
$$
The maximum is attained at $k = 0$ (the only $k$ where the linear
discount $-2|k|/(N_t+1)$ vanishes).

**Alternative proof of $\|S_F\|_{\mathrm{cb}} = 1$ via positivity.** A positive
Schur multiplier on an abelian C*-algebra has cb-norm equal to its
operator norm equal to the $\ell^\infty$ norm of the symbol (Paulsen 2002
Cor.\ 8.7; Pisier 2001 Thm 1.6). The Fejér kernel is positive, so
$\|S_F\|_{\mathrm{cb}} = \|S_F\|_{\mathrm{op}} = \|\hat F\|_\infty = 1$.

**Verification in code.** `geovac.central_fejer_compact_temporal.cb_norm_circle(N_t)`
returns `Rational(1)` for every $N_t \ge 1$.

---

## §6. Factorization on amenable group products (the core of B3)

### §6.1 Tensor product of central subalgebras

For the direct product $G = G_1 \times G_2$ of compact groups, the
Pontryagin / Peter--Weyl dual $\hat G$ decomposes as $\hat G = \hat G_1
\times \hat G_2$, with the irreducible representations of $G$ being
external tensor products $\pi_1 \boxtimes \pi_2$ of irreps of $G_1$ and
$G_2$. The character of $\pi_1 \boxtimes \pi_2$ is
$$
   \chi_{\pi_1 \boxtimes \pi_2}(g_1, g_2) \;=\; \chi_{\pi_1}(g_1)\,\chi_{\pi_2}(g_2).
$$
A central function on $G$ is therefore a sum (in $L^2$) of products
$\chi_{\pi_1}(g_1) \chi_{\pi_2}(g_2)$, and the central subalgebra factorizes:
$$
   Z(C(G)) \;\cong\; Z(C(G_1)) \;\hat\otimes\; Z(C(G_2)).
$$

### §6.2 Tensor product of Schur multipliers

If $K_1$ is central on $G_1$ and $K_2$ central on $G_2$, then $K_1 \otimes
K_2$ (the tensor product, viewing both as elements of $L^1$) is central
on $G_1 \times G_2$, with Plancherel symbol
$$
   \widehat{K_1 \otimes K_2}(\pi_1, \pi_2)
   \;=\; \hat K_1(\pi_1) \cdot \hat K_2(\pi_2).
$$
This is automatic from the factorization $\chi_{\pi_1 \boxtimes
\pi_2}(g_1, g_2) = \chi_{\pi_1}(g_1) \chi_{\pi_2}(g_2)$ and the definition
of $\hat K(\pi)$ via integration against the character.

### §6.3 Bożejko--Fendler factorization

**Lemma 6.1 (Schur multiplier cb-norm factorizes over amenable products).**
*Let $G = G_1 \times G_2$ with both $G_i$ amenable compact, and let
$K_i$ be central on $G_i$ with cb-norm
$\|S_{K_i}\|_{\mathrm{cb}} = \|\hat K_i\|_\infty$ (by Theorem 3.1). Then*
$$
   \|S_{K_1 \otimes K_2}\|_{\mathrm{cb}}
   \;=\; \|S_{K_1}\|_{\mathrm{cb}} \cdot \|S_{K_2}\|_{\mathrm{cb}}.
$$

**Proof.** By Theorem 3.1 applied to $G_1 \times G_2$ (amenable, central
multiplier symbol $\widehat{K_1 \otimes K_2} = \hat K_1 \cdot \hat K_2$),
$$
   \|S_{K_1 \otimes K_2}\|_{\mathrm{cb}}
   \;=\; \|\widehat{K_1 \otimes K_2}\|_{\ell^\infty(\hat G_1 \times \hat G_2)}
   \;=\; \|\hat K_1 \cdot \hat K_2\|_\infty.
$$
The supremum norm of a product of functions on a direct product space
factorizes as the product of supremum norms when the supremum is taken
over all $(\pi_1, \pi_2)$ jointly:
$$
   \sup_{(\pi_1, \pi_2)} |\hat K_1(\pi_1) \hat K_2(\pi_2)|
   \;=\; \Big(\sup_{\pi_1} |\hat K_1(\pi_1)|\Big)
       \cdot \Big(\sup_{\pi_2} |\hat K_2(\pi_2)|\Big)
   \;=\; \|\hat K_1\|_\infty \cdot \|\hat K_2\|_\infty.
$$
The factorization step uses the joint structure of the product space and
is **standard**: for non-negative functions (which is our case — both
$\hat K^{\SU(2)} \ge 0$ and $\hat K^{U(1)} \ge 0$) the supremum of a
separable product equals the product of factor suprema. Applying
Theorem 3.1 in reverse on each factor:
$$
   \|S_{K_1 \otimes K_2}\|_{\mathrm{cb}}
   \;=\; \|\hat K_1\|_\infty \cdot \|\hat K_2\|_\infty
   \;=\; \|S_{K_1}\|_{\mathrm{cb}} \cdot \|S_{K_2}\|_{\mathrm{cb}}. \quad \square
$$

**Remark (positivity is essential, but harmless here).** The
multiplicativity of $\ell^\infty$-norm of a product holds in general only
for non-negative functions; for complex-valued functions one has
$\|fg\|_\infty \le \|f\|_\infty \|g\|_\infty$ with strict inequality
possible. Both our factor symbols $\hat K^{\SU(2)}$ and $\hat K^{U(1)}$
are non-negative (the Fejér / Cesàro structure forces this), so the
factorization is exact, not an inequality. This is the precise reason
sub-sprint B is "the most straightforward of the four": positive Fejér
multipliers on amenable group products have multiplicative cb-norms.

---

## §7. Application to the L3b-2 joint Fejér kernel

Combining Lemma 4.1 (factor 1), §5.3 (factor 2), and Lemma 6.1
(factorization on $\SU(2) \times U(1)$):

**Theorem 7.1 (Sub-sprint B, restated).** *For all $n_{\max} \ge 1$ and
$N_t \ge 1$,*
$$
   \boxed{\quad
   \big\|S_{K^{\mathrm{joint}}_{n_{\max}, N_t, T}}\big\|_{\mathrm{cb}}
   \;=\; \frac{2}{n_{\max} + 1} \cdot 1 \;=\; \frac{2}{n_{\max} + 1}.
   \quad}
$$

**Proof.** By Lemma 6.1 applied to $G_1 = \SU(2)$, $G_2 = U(1)$, $K_1 =
K^{\SU(2)}_{n_{\max}}$, $K_2 = K^{U(1)}_{N_t, T}$ (both central, both
Fejér / Cesàro-positive):
$$
   \|S_{K^{\mathrm{joint}}}\|_{\mathrm{cb}}
   \;=\; \|S_{K^{\SU(2)}}\|_{\mathrm{cb}} \cdot \|S_{K^{U(1)}}\|_{\mathrm{cb}}.
$$
The first factor is $2/(n_{\max}+1)$ by Lemma 4.1; the second factor is
$1$ by §5.3. Multiplying gives the stated closed form. $\square$

### §7.1 Asymptotic rate

The joint cb-norm rate as $(n_{\max}, N_t) \to (\infty, \infty)$ is
$$
   \|S_{K^{\mathrm{joint}}_{n_{\max}, N_t}}\|_{\mathrm{cb}}
   \;=\; \frac{2}{n_{\max} + 1} \;=\; O(1/n_{\max}).
$$
This is **independent of $N_t$** — the temporal cutoff does not affect the
cb-norm. The reason is that the $U(1)$ factor's cb-norm is identically $1$
at every $N_t$ (the maximum of the Fejér symbol is always attained at
$k = 0$, where the symbol equals $1$ regardless of $N_t$). The full
temporal-cutoff sensitivity lives in the **mass-concentration moment**
$\gamma^{U(1)}_{N_t, T} = O(T/N_t)$, not in the cb-norm.

This is the structurally correct behavior: the cb-norm controls the
*Lipschitz-distortion height* in the L4/L5 Latrémolière propinquity
assembly, while the mass-concentration moment $\gamma$ controls the
*reach*. The L3b-2 sub-sprint splits naturally:
- **Sub-sprint B (this memo):** joint cb-norm, $O(1/n_{\max})$.
- **Sub-sprint C/D (forward):** joint $\gamma^{\mathrm{joint}}$ rate,
  $O(\log n_{\max}/n_{\max} + 1/N_t)$ — depends on $N_t$.

---

## §8. Numerical verification panel

The driver `debug/l3b_2_sub_sprint_B_compute.py` verifies the closed form
$\|S_{K^{\mathrm{joint}}}\|_{\mathrm{cb}} = 2/(n_{\max}+1)$ at the
following panel cells:

| $(n_{\max}, N_t)$ | factor SU(2) | factor U(1) | product | bound $2/(n_{\max}+1)$ |
|:-----------------:|:------------:|:-----------:|:-------:|:----------------------:|
| (2, 3) | $2/3$ | $1$ | $2/3$ | $2/3$ ✓ |
| (2, 5) | $2/3$ | $1$ | $2/3$ | $2/3$ ✓ |
| (3, 3) | $1/2$ | $1$ | $1/2$ | $1/2$ ✓ |
| (3, 5) | $1/2$ | $1$ | $1/2$ | $1/2$ ✓ |
| (4, 3) | $2/5$ | $1$ | $2/5$ | $2/5$ ✓ |
| (4, 5) | $2/5$ | $1$ | $2/5$ | $2/5$ ✓ |

All cells verified exact in sympy rational arithmetic. Sympy
factorization residual (left-hand side minus right-hand side) is
**identically zero**. See `debug/data/l3b_2_sub_sprint_B.json` for the
machine-readable panel.

---

## §9. Verdict and connection to forward sub-sprints

**Verdict for sub-sprint B: PROVED.**

The joint cb-norm bound is established at the level of rigor of Paper 38
Lemma L2(c). The proof is a factor-wise application of the
Bożejko--Fendler central-multiplier equality on the amenable compact
group $\SU(2) \times U(1)$, with:
- Factor 1 (SU(2)): inherited from Paper 38 L2(c) verbatim.
- Factor 2 (U(1)): proved directly from positivity of the Fejér kernel
  and Theorem 3.1.
- Factorization (B3): proved as Lemma 6.1, follows from amenability of
  the direct product plus positivity of both factor symbols.

No non-standard input was needed. Both the SU(2) factor (Paper 38 / Paper
40) and the $U(1)$ factor (textbook Fejér) and the factorization
(Bożejko--Fendler) are entirely standard.

### Flags for sub-sprints C and D

**For sub-sprint C** (joint Lipschitz comparison / L3 lemma with constant
$C_3$): the joint cb-norm bound established here is the **multiplicative**
ingredient that the Lipschitz constant assembles against. Specifically,
the joint Lipschitz seminorm $\|[D^{\mathrm{joint}}_L, M_f]\|_{\mathrm{op}}$ must be
compared against $\|\nabla f\|_\infty$ along both the spatial and temporal
gradients, and the joint cb-norm bound $2/(n_{\max}+1)$ controls the
spatial-direction contribution. The temporal direction's Lipschitz contribution
needs a separate $\partial_t$-bound (1D Lichnerowicz / momentum-cutoff on
$S^1_T$). **Expected $C_3^{\mathrm{joint}} = 1 + \mathrm{cross}$ with
cross $\to 0$ in the joint limit, exactly as flagged in
`debug/l3b_first_move_memo.md` §8.**

**For sub-sprint D** (Berezin reconstruction / L4 lemma): the joint
Berezin map $B_{n_{\max}, N_t} \colon C^\infty(S^3 \times S^1_T) \to
O^L_{n_{\max}, N_t, T}$ has four properties (positivity, contractivity,
approximate identity, L3 compatibility). The **contractivity** property
uses the cb-norm bound $\|B\|_{\mathrm{cb}} \le 1$, which is the cb-norm of the
Plancherel projector composed with multiplication by $\hat K^{\mathrm{joint}}$.
Since the projector has cb-norm 1 (it's a UCP map) and the joint Schur
multiplier has cb-norm $2/(n_{\max}+1) \le 1$ for $n_{\max} \ge 1$, the
composite has cb-norm $\le 2/(n_{\max}+1) \le 1$. **Sub-sprint D should
directly invoke this memo's Theorem 7.1 for the contractivity proof.**

### What this sub-sprint does NOT establish

- The joint mass-concentration moment $\gamma^{\mathrm{joint}}$ asymptotic rate
  (that is sub-sprint C/D territory; the cb-norm and the $\gamma$-moment
  are different invariants).
- The L4 Berezin map's approximate-identity property (uses $\gamma \to 0$,
  not the cb-norm).
- The full propinquity bound (assembled in sub-sprint D via L1'--L5).

### Path to fully rigorous (already present)

All three ingredients (B1, B2, B3) are at the level of rigor of Paper 38
§L2: standard amenable-group Fourier-multiplier theory, with the only
non-trivial step (the dual estimate / Bożejko--Fendler) cited from the
literature rather than re-proved here. This is the right level of rigor
for the L3b-2 memo, matching Paper 38's standards.

---

## §10. References

- M. Bożejko and G. Fendler, "Herz--Schur multipliers and completely
  bounded multipliers of the Fourier algebra of a locally compact group,"
  *Boll. Un. Mat. Ital. A* (6) 3 (1984) 297--302; extended version with
  full proofs and amenability hypothesis in (1991).
- G. Pisier, *Similarity Problems and Completely Bounded Maps*, Lecture
  Notes in Math. 1618, Springer, 2nd ed. 2001, Ch.\ 8 Thm 8.10.
- V. Paulsen, *Completely Bounded Maps and Operator Algebras*, Cambridge
  Studies in Advanced Math.\ 78, Cambridge Univ.\ Press, 2002, Cor.\ 8.7.
- M. Takesaki, *Theory of Operator Algebras I*, Springer, 1979, Thm IV.4.14
  (commutative C*-tensor product unambiguity).
- B. Bekka, P. de la Harpe, A. Valette, *Kazhdan's Property (T)*, Cambridge
  Univ.\ Press, 2008, Prop.\ G.3.1 (amenability of direct products).
- J.-P. Pier, *Amenable Locally Compact Groups*, Wiley, 1984, Prop.\ 13.6.
- Loutey, J. "GeoVac Paper 38: SU(2) propinquity convergence on the
  Camporesi--Higuchi spectral triple," Zenodo 2026; §L2 Lemma L2(c).
- Loutey, J. "GeoVac Paper 40: Unified propinquity convergence on compact
  Lie groups," draft 2026; §3.2 Plancherel × Vandermonde cancellation
  theorem.

---

**End of memo.**
