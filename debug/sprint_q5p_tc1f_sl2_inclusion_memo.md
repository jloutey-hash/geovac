# Sprint Q5'-TC-1f — Explicit $U^*_{\mathrm{Levi}}$ inclusion (Aut$^\otimes(\omega)$, $SL_2$ + full-panel injectivity)

**Date:** 2026-06-06
**Verdict:** POSITIVE — bit-exact at every cell; 490/490 residuals.
**Module:** `geovac/tannakian.py` (added `PWRep`, `PWMorphism`, `sl2_*` actions, `verify_sl2_*`, `verify_ga_sl2_commute`, `primitive_generator_rep`, `verify_injectivity_at_generator`; ~250 LOC).
**Driver:** `debug/compute_q5p_tc1f_sl2_inclusion.py` (490 residuals).
**Tests:** `tests/test_tannakian_sl2.py` (60 tests, 1.03s).
**Data:** `debug/data/q5p_tc1f_sl2_inclusion.json`.

## What this sprint adds to the Tannakian closure arc

TC-1e established that on every $V \in \mathrm{Rep}_{\mathbb{Q}}(\mathcal{H}_{\mathrm{GV}})$ the unipotent automorphism
$\eta_V(t) = \exp\!\bigl(\sum_g t_g X_g^V\bigr)$ is a natural $\otimes$-automorphism of $\omega$. That gave the
$\mathbb{G}_a^{3N}$ piece of $U^*_{\mathrm{Levi}} = \mathbb{G}_a^{3N(n_{\max})} \rtimes SL_2$ as a categorical statement.

TC-1f closes the **explicit** inclusion with three structural witnesses:

1. **$SL_2$ axioms on Peter–Weyl reps** — invertibility, tensor compatibility, and the group-homomorphism law are each bit-exact on the standard rep $\mathbb{Q}^2$, $\mathrm{Sym}^2 \mathbb{Q}^2$, and the trivial rep.
2. **$\mathbb{G}_a \times SL_2$ commutativity on combined reps $V \otimes V_{\mathrm{PW}}$** — the two factors act on different tensor factors, so by the Kronecker identity $(A \otimes I)(I \otimes B) = (I \otimes B)(A \otimes I)$ they commute identically.
3. **Full-panel injectivity at $n_{\max} = 2$** — for each of the $3 \cdot N(2) = 15$ primitive generators $g$ we build a faithful 2-dim test rep $V_g$ with $X_g^{V_g} = E_{12}$ and $X_h^{V_g} = 0$ for $h \ne g$. Then $\eta_{V_g}(t) = I$ iff $t_g = 0$. Across 255 (case, generator) cells (zero $t$, single-generator $t$, generic three-generator $t$) the predicate is bit-exact.

Together: at $n_{\max} = 2$ the inclusion $U^*_{\mathrm{Levi}}(\mathbb{Q}) \hookrightarrow \mathrm{Aut}^\otimes(\omega)$ is explicit on every primitive generator and on the $SL_2$ factor, with $\mathbb{G}_a$ and $SL_2$ provably commuting on combined reps.

## Numerical panel

| Block                              | Residuals | Bit-exact |
| :--------------------------------- | --------: | --------: |
| $SL_2$ invertibility               |        15 |        15 |
| $SL_2$ tensor compatibility        |        45 |        45 |
| $SL_2$ group homomorphism          |        75 |        75 |
| $\mathbb{G}_a \times SL_2$ commute |       100 |       100 |
| Full-panel injectivity ($n=2$)     |       255 |       255 |
| **Total**                          | **490**   | **490**   |

$SL_2(\mathbb{Q})$ panel:\ $\{e, \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}, \mathrm{diag}(2, 1/2), \begin{pmatrix} 5 & 2 \\ 7 & 3 \end{pmatrix}, \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}\}$ (identity, unipotent, torus, generic, Weyl involution).

PW rep panel:\ $\{V_{\mathrm{triv}}^{(k)}, V_{\mathrm{fund}} = \mathbb{Q}^2, \mathrm{Sym}^2 V_{\mathrm{fund}} = \mathbb{Q}^3\}$.

Primitive-generator panel at $n_{\max} = 2$:\ all 15 generators $(n, l, k)$ with $(n, l) \in \{(1, 0), (1, 1), (2, 0), (2, 1), (2, 2)\}$ and $k \in \{0, 1, 2\}$ (one per (sector, Mellin slot) pair, see `geovac/pro_system.py:primitive_generators`).

## What's structurally new vs TC-1e

TC-1e proved that every $\eta_V(t) \in \mathrm{Aut}^\otimes(\omega)$ (invertibility, unit, naturality, tensor compatibility, group law). That is the categorical inclusion of $\mathbb{G}_a^{3N}$.

TC-1f makes the $SL_2$ factor explicit at the same Peter–Weyl level:\ it builds an honest representation $\rho:\ SL_2(\mathbb{Q}) \to \mathrm{GL}(V_{\mathrm{PW}})$ (standard and $\mathrm{Sym}^2$) and verifies the same three structural identities that TC-1e proved for $\mathbb{G}_a$:\ invertibility, tensor compatibility, group law. The fourth identity — **$\mathbb{G}_a$ and $SL_2$ commute on combined reps** — is the Levi-decomposition compatibility, and it is bit-exact by construction (the two factors act on disjoint tensor factors).

The full-panel injectivity block is the explicit witness that the inclusion is faithful at $n_{\max} = 2$:\ no primitive generator is "wasted" — there is a distinct test rep $V_g$ that detects it.

## Scope and follow-ons

* **Scope** of the explicit inclusion:\ on the abelian primitive substrate $\mathcal{H}_{\mathrm{GV}} = \mathrm{Sym}_{\mathbb{Q}}(V)$ with the Levi decomposition $U^*_{\mathrm{Levi}} = \mathbb{G}_a^{3N} \rtimes SL_2$ established in v3.63.0. The Peter–Weyl panel here is the standard + $\mathrm{Sym}^2$ — these are the only two irreducible $SL_2$ reps needed to **see** the inclusion at a witnessing level.
* **Converse equality** $\mathrm{Aut}^\otimes(\omega) = U^*_{\mathrm{Levi}}$ — i.e. that no other natural $\otimes$-automorphisms exist — is the Tannakian-reconstruction direction (Deligne–Milne 1982 Theorem 2.11). That is a multi-year structural claim and **is NOT closed by TC-1f**; this sprint closes only the inclusion direction.
* **Higher-$n_{\max}$ generalization** of the full-panel injectivity is the same construction:\ $3 \cdot N(n_{\max})$ generators give $3 \cdot N(n_{\max})$ test reps; the analogous predicate is bit-exact at every cutoff by the same Jordan-block argument. Carrying it out at $n_{\max} = 3, 4$ would extend the panel but not the structure.

## Paper-edit recommendation

Add one paragraph to `papers/group3_foundations/paper_55_periods_of_geovac.tex` §subsec:open_m2_m3 (after the TC-1e paragraph) headlining the explicit $\mathbb{G}_a \times SL_2$ inclusion and full-panel injectivity at $n_{\max} = 2$, with the 490/490 residual count and a pointer to this memo.

## Release plan

CLAUDE.md §2 one-liner + §1 version bump to v3.74.0 + CHANGELOG entry. After release, immediately draft Paper 56 (math.OA standalone consolidating the entire PS-1/2/3/4 + TC-1a/b/c/d/e/f arc).
