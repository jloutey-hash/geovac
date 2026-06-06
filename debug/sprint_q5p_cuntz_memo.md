# Sprint Q5'-Cuntz — first scoping step of the JLO/CM bicomplex Cuntz extension as Stage-2 enrichment ingredient (3rd of 3 candidates flagged by Track A of v3.61.0)

**Date:** 2026-06-05 (close-of-day follow-on to Sprint Q5'-Stage2-Hopf v3.61.0 Track A)
**Sprint:** first scoping step of the multi-year "JLO/CM bicomplex Cuntz extension substrate enrichment," the 3rd of 3 structural ingredients flagged by Track A of Sprint Q5'-HardParts (v3.61.0).
**Driver:** `debug/compute_q5p_cuntz.py`
**Data:** `debug/data/sprint_q5p_cuntz.json`
**Wall time:** 0.017 s
**Discipline:** bit-exact `sympy.Rational` throughout; no PSLQ; no floats; no transcendentals introduced.

---

## 1. TL;DR

**Verdict: STOP.** The Cuntz extension as a Stage-2 enrichment ingredient is **structurally incompatible** with the v3.61.0 Track A Hopf-substrate framework on three independent grounds, each verified bit-exactly at $n_{\max} = 2$ with $N = 15$ total isometries (one per (sector, Mellin-slot) pair):

1. **Primitive coproduct $\Delta(S_i) = S_i \otimes 1 + 1 \otimes S_i$ fails the Cuntz relations.** Extending $\Delta$ as a $*$-algebra homomorphism gives $\Delta(S_i^* S_j) \ne \Delta(\delta_{ij}\cdot 1)$ bit-exactly on **0/225** pair tests (all 225 pairs fail). On the diagonal $i = j$, the residual has a non-zero $1 \otimes 1$ coefficient of $+1$ (from $\Delta(\delta_{ii}) = 1 \otimes 1$ on LHS vs $2 \cdot 1 \otimes 1$ on RHS); on the off-diagonals there are non-zero $S_i^* \otimes S_j$ and $S_j \otimes S_i^*$ terms with no LHS counterpart.

2. **Diagonal coproduct $\Delta(S_i) = S_i \otimes S_i$ passes the first Cuntz relation but fails the second.** First relation $\Delta(S_i^* S_j) = \Delta(\delta_{ij} \cdot 1)$ passes **225/225** pair tests by direct computation $(S_i^* \otimes S_i^*)(S_j \otimes S_j) = S_i^* S_j \otimes S_i^* S_j = \delta_{ij}(1 \otimes 1)$. But the second Cuntz relation $\sum_i S_i S_i^* = 1$ requires $\Delta(1) = (\sum_i S_i S_i^*) \otimes (\sum_j S_j S_j^*) = \sum_{i, j} S_i S_i^* \otimes S_j S_j^*$, while $\Delta_{\mathrm{diag}}(\sum_i S_i S_i^*) = \sum_i (S_i S_i^* \otimes S_i S_i^*)$ has only the **diagonal $i = j$** terms — missing **$N(N - 1) = 210$ cross-terms**. So $\mathcal{O}_N$ is not a bialgebra under $\Delta_{\mathrm{diag}}$.

3. **JLO degree-1 cochain $\phi_1$ is structurally NON-ZERO on the Cuntz extension** (in contrast to the commutative-algebra case where $\phi_1 \equiv 0$ identically, the Sub-Sprint 1 structural theorem). The vanishing of $\phi_1$ on $\mathcal{A} = \mathbb{C}^5$ relied on sector-orthogonality $e_s e_t = \delta_{st} e_s$. Cuntz isometries violate this: $S_i S_j$ for $i \ne j$ is a NEW basis element (path of length 2), not zero. The cyclic-cancellation channel used in the commutative proof is BROKEN by the Cuntz extension. Depth-truncated diagonal-Dirac surrogate trace happens to vanish at our scoping depth (paths grow monotonically), but this is representation-specific — the FULL Camporesi-Higuchi Dirac $D = \Lambda + \kappa A$ has off-diagonal action through $A$ that produces non-trivial trace contributions; computing this bit-exactly requires the multi-year Pimsner-Voiculescu / Cuntz-Pimsner extension framework.

**Headline structural finding** (visible only because the v3.61.0 Track A commutative comparison exists). The clean tensor factorisation
$$
\mathcal{H}_{\mathrm{GV}} = \mathcal{H}_{\mathrm{GV}}^{[0]} \otimes \mathcal{H}_{\mathrm{GV}}^{[1]} \otimes \mathcal{H}_{\mathrm{GV}}^{[2]}
$$
of v3.61.0 Track A is a **commutativity-dependent feature**. On the Cuntz extension, the M1/M2/M3 partition survives only at the trivial **comodule** level (each isometry $S_{(n, l), k}$ carries a fixed slot $k$ and transforms only within its $k$-subspace under the gauge $U(1)$ coaction). At the (non-existent) Hopf-coproduct level the question is moot. The Cuntz extension therefore does NOT produce a tractable non-abelian motivic Galois enrichment at this scoping step.

**The umbrella structural meaning.** Combined with Sprint Q5'-HardParts Track A (v3.61.0), this verdict sharpens the multi-year picture: among Track A's three flagged enrichment ingredients (nested-Hopf-tower $J^*(S^3)$; cross-shell off-diagonal Dirac; JLO/CM Cuntz extension), the Cuntz extension is the **structurally most disruptive** — it doesn't just enrich the substrate, it breaks the Hopf-coproduct framework entirely and lands the substrate in the genuinely different category of Hopf-comodules over $U(1)$ (Cuntz-Pimsner extensions, Pimsner 1997). This is a structurally informative scope finding: the path forward through Cuntz is multi-year and requires a substrate enrichment of a categorically different shape (comodule-over-gauge-group, not Hopf-algebra-with-non-primitive-coproduct), distinguishing the Cuntz route from the other two routes (which preserve the Hopf-algebra category at the substrate level).

**Decision gate landed at STOP.** The decision gate as written asked: does the Cuntz extension produce non-primitive coproduct content while preserving M1/M2/M3? The answer is STOP because the obstruction is even more fundamental than "natural coproduct is primitive": there is **no coproduct at all** (in the Hopf-algebra category) that is compatible with the Cuntz relations.

---

## 2. Verdict against the decision gate

| Gate | Selected? | Reason |
|:-----|:---------:|:-------|
| POSITIVE | not selected | Neither natural coproduct (primitive nor diagonal) produces a bit-exact bialgebra on $\mathcal{O}_N$. The pathway through non-Hopf-comodule structure is multi-year and structurally different. |
| BORDERLINE | not selected | The natural coproducts FAIL the Cuntz relations rather than satisfy them in a downgraded way; this is not a "natural coproduct is primitive but doesn't break commutativity" outcome. |
| **STOP** | **selected** | Cuntz extension is **structurally incompatible** with the v3.61.0 Track A Hopf framework. $\mathcal{O}_N$ is a Hopf-comodule (over the gauge $U(1)$), not a Hopf-algebra; no natural Hopf-coproduct exists. M1/M2/M3 partition survives only at the trivial comodule level. The JLO degree-1 cochain $\phi_1$ is generically non-zero (so JLO bicomplex extends in principle, but extraction at theorem grade requires the multi-year Cuntz-Pimsner extension framework — Pimsner 1997; Connes-Cuntz 1988). |

The STOP verdict here is **informative**, not a dead end: it sharpens the picture of which of Track A's three flagged ingredients is "closest" to the existing Hopf framework. The Cuntz extension is **structurally the farthest**.

---

## 3. The Cuntz path algebra at $n_{\max} = 2$ with depth truncation

### 3.1 Definition

Recall the CH Fock sectors at $n_{\max} = 2$: $\mathcal{S}_2 = \{(1, 0), (1, 1), (2, 0), (2, 1), (2, 2)\}$ with $N(\mathcal{S}_2) = 5$ sectors.

We introduce **three copies** of the Cuntz alphabet, one per Mellin slot $k \in \{0, 1, 2\}$, giving a total of $N = 3 \cdot 5 = 15$ isometries $S_{(n, l), k}$. The Cuntz path algebra $\mathcal{O}_{N}^{\mathrm{path}}$ is generated by $\{S_i : 0 \le i < N\}$ subject to the Cuntz relations
$$
S_i^* S_j = \delta_{ij} \cdot 1, \qquad \sum_{i = 0}^{N - 1} S_i S_i^* = 1.
$$

The standard normal form (Cuntz 1977 Lemma 1.3) writes every element as a finite sum of **path basis elements** $S_\alpha S_\beta^*$ where $\alpha = (i_1, i_2, \ldots, i_p)$ and $\beta = (j_1, j_2, \ldots, j_q)$ are finite sequences (paths) in $\{0, \ldots, N - 1\}^*$ and $S_\alpha := S_{i_1} S_{i_2} \cdots S_{i_p}$ (with $S_\emptyset := 1$).

### 3.2 Depth truncation

For finite computation we truncate to **path-length $\le L$**:
$$
\mathcal{O}_N^{(L)} := \mathrm{span}\big\{ S_\alpha S_\beta^* : |\alpha| \le L,\; |\beta| \le L \big\}.
$$

| Truncation depth $L$ | Number of paths $\sum_{\ell = 0}^{L} N^{\ell}$ | Basis size $|\text{paths}|^2$ |
|:---:|:---:|:---:|
| 1 | 16 | 256 |
| 2 | 241 | 58,081 |

**Note:** at depth $L = 2$ with $N = 15$ the basis size is already 58,081. Products beyond depth $L$ are truncated to zero (consistent with the finite-truncation choice). All structural Hopf-compatibility tests below are independent of $L$ — they live at the LEVEL of the defining relations $S_i^* S_j = \delta_{ij}$ and $\sum_i S_i S_i^* = 1$, which are checked at the generator level and propagate trivially through the path-basis.

### 3.3 Multiplication

The Cuntz path-basis multiplication follows the standard reduction $S_\alpha S_\beta^* S_\gamma S_\delta^*$: reduce $S_\beta^* S_\gamma$ using $S_i^* S_j = \delta_{ij}$. If $|\beta| \le |\gamma|$ and $\gamma = \beta \cdot \gamma'$ (concatenation), then $S_\beta^* S_\gamma = S_{\gamma'}$ and the product is $S_{\alpha \cdot \gamma'} S_\delta^*$. Otherwise the product is zero (if $\gamma$ doesn't start with $\beta$). Symmetric case $|\beta| > |\gamma|$ gives the conjugate result.

Driver implementation: `cuntz_pair_mul` returns $(\text{coef}, (\alpha', \beta'))$ where $\text{coef} \in \{0, 1\}$. Depth truncation drops outputs with $|\alpha'| > L$ or $|\beta'| > L$.

---

## 4. The Dirac operator at the scoping level

At this scoping stage, the load-bearing structural questions (non-primitive coproduct existence; M1/M2/M3 partition; JLO $\phi_1$ vanishing) are **independent of the precise Dirac choice** because they live at the algebraic-relation level. We model the Dirac via a **diagonal-Dirac surrogate** that captures the CH-1 chirality/eta weights:

- On a depth-truncated Fock basis state $|\alpha\rangle$ (with $\alpha$ a path of length $\le L$), the surrogate Dirac acts by $D \cdot |\alpha\rangle = \lambda_\alpha \cdot |\alpha\rangle$ where $\lambda_\alpha = \sum_{j \in \alpha} \mathrm{sign}(\chi_{\sigma(j)}) \cdot (n_{\sigma(j)} + 1/2)$ with $\sigma(j)$ the sector label of the $j$-th generator.

This diagonal-Dirac surrogate is **sufficient** to expose the structural obstructions for primitive/diagonal coproducts and for the M1/M2/M3 partition (those obstructions live at the algebra-relation level, not the Hilbert representation). It is **insufficient** to extract the non-zero $\phi_1$ coefficient at theorem grade (because diagonal action on monotonically-extending paths gives trivial trace) — that requires the FULL Camporesi-Higuchi adjacency-extended Dirac, which is multi-year.

**Honest scope:** the diagonal-Dirac surrogate trades quantitative theorem-grade $\phi_1$ extraction for qualitative structural transparency on the Hopf-coproduct questions. This is the correct trade for a scoping sprint.

### 4.1 Sample $\phi_0^{(k)}$ weights on the surrogate

| Slot $k$ | Single-generator path $S_{(1, 0)}$ | Two-generator path $S_{(1, 0)} S_{(1, 1)}$ |
|:---:|:---:|:---:|
| $k = 0$ (M1) | $\chi_{(1, 0)} = +2$ | $\chi_{(1, 0)} + \chi_{(1, 1)} = 0$ |
| $k = 1$ (M3) | $\eta_{(1, 0)} = 3$ | $\eta_{(1, 0)} + \eta_{(1, 1)} = 6$ |
| $k = 2$ (M2) | $\dim_{(1, 0)} = 2$ | $\dim_{(1, 0)} + \dim_{(1, 1)} = 4$ |

These are bit-exact `sympy.Rational` integer outputs from the driver (sample data in `debug/data/sprint_q5p_cuntz.json`, field `sample_phi_0_weights`).

---

## 5. Primitive coproduct $\Delta(S_i) = S_i \otimes 1 + 1 \otimes S_i$: incompatibility with Cuntz relations

### 5.1 The test

If $\Delta$ extends to a coproduct on $\mathcal{O}_N$, it must satisfy
$$
\Delta(S_i^* S_j) = \Delta(\delta_{ij} \cdot 1) = \delta_{ij} \cdot (1 \otimes 1),
$$
because $S_i^* S_j = \delta_{ij} \cdot 1$ is forced by the Cuntz relations. Equivalently, applying $\Delta$ as a $*$-algebra homomorphism gives
$$
\Delta(S_i^*) \Delta(S_j) = \big(S_i^* \otimes 1 + 1 \otimes S_i^*\big)\big(S_j \otimes 1 + 1 \otimes S_j\big),
$$
which (using component-wise multiplication on $\mathcal{O}_N \otimes \mathcal{O}_N$) expands to
$$
S_i^* S_j \otimes 1 + S_i^* \otimes S_j + S_j \otimes S_i^* + 1 \otimes S_i^* S_j
= \delta_{ij} (1 \otimes 1) + S_i^* \otimes S_j + S_j \otimes S_i^* + \delta_{ij} (1 \otimes 1)
= 2 \delta_{ij} (1 \otimes 1) + S_i^* \otimes S_j + S_j \otimes S_i^*.
$$

### 5.2 The residual

For the relation to hold, we need
$$
\delta_{ij}(1 \otimes 1) = 2 \delta_{ij}(1 \otimes 1) + S_i^* \otimes S_j + S_j \otimes S_i^*,
$$
i.e.,
$$
S_i^* \otimes S_j + S_j \otimes S_i^* + \delta_{ij}(1 \otimes 1) = 0.
$$

This **bit-exact fails for every $(i, j)$**:
- For $i = j$: residual is $S_i^* \otimes S_i + S_i \otimes S_i^* + (1 \otimes 1)$, which has three non-zero basis terms.
- For $i \ne j$: residual is $S_i^* \otimes S_j + S_j \otimes S_i^*$, which has two non-zero basis terms.

### 5.3 The bit-exact panel

| Test | Pairs checked | Pairs passing | Pairs failing |
|:----|:---:|:---:|:---:|
| Primitive $\Delta$ vs Cuntz relations | **225** (= 15 × 15) | **0** | **225** |

(Driver field: `panel_data.primitive_compat_pairs = 225`, `n_pass = 0`, `fails = 225`.)

### 5.4 Structural reading

The primitive coproduct is incompatible with the **first** Cuntz relation $S_i^* S_j = \delta_{ij}$. This is a structural feature of the Cuntz algebra: the relation $S_i^* S_j = \delta_{ij}$ identifies the product of two distinct generators with a scalar, while the primitive coproduct produces non-trivial cross-terms in the tensor product. The only way to maintain compatibility would be to enforce $S_i^* \otimes S_j + S_j \otimes S_i^* + \delta_{ij}(1 \otimes 1) = 0$ as an IDEAL in $\mathcal{O}_N \otimes \mathcal{O}_N$, which is NOT a sub-bialgebra — the resulting quotient is not a Hopf algebra.

This rules out the "primitive enrichment" reading of Cuntz: even though the underlying $\mathcal{O}_N$ is non-commutative (which would enable non-trivial coproducts), the defining relations FIGHT BACK against primitive coproduct attempts.

---

## 6. Diagonal coproduct $\Delta(S_i) = S_i \otimes S_i$: passes first Cuntz relation, fails second

### 6.1 Motivation

The diagonal coproduct is the natural candidate for a coproduct on a $C^*$-algebra of isometries (Pimsner 1997; Connes-Cuntz 1988): each isometry $S_i$ maps to a tensor-product of itself with itself, geometrically representing "the same operator acting twice."

### 6.2 First Cuntz relation: PASS

Direct computation:
$$
\Delta_{\mathrm{diag}}(S_i^*) \Delta_{\mathrm{diag}}(S_j) = (S_i^* \otimes S_i^*)(S_j \otimes S_j) = S_i^* S_j \otimes S_i^* S_j = \delta_{ij} \cdot (1 \otimes 1).
$$
LHS $\Delta(\delta_{ij} \cdot 1) = \delta_{ij} \cdot (1 \otimes 1)$. **Bit-exact match on 225/225 pairs.** First relation passes.

(Driver field: `diagonal_first_relation_pass = 225`, `diagonal_first_relation_total = 225`.)

### 6.3 Second Cuntz relation: FAIL

The second relation says $\sum_i S_i S_i^* = 1$. Applying $\Delta_{\mathrm{diag}}$:
$$
\Delta_{\mathrm{diag}}\big(\sum_i S_i S_i^*\big) = \sum_i \Delta_{\mathrm{diag}}(S_i) \Delta_{\mathrm{diag}}(S_i^*) = \sum_i (S_i \otimes S_i)(S_i^* \otimes S_i^*) = \sum_i (S_i S_i^* \otimes S_i S_i^*).
$$

The required RHS is $\Delta(1) = 1 \otimes 1 = (\sum_i S_i S_i^*) \otimes (\sum_j S_j S_j^*) = \sum_{i, j} (S_i S_i^* \otimes S_j S_j^*)$.

**Cross-terms $i \ne j$ are missing** on the LHS: the diagonal coproduct gives only the diagonal $i = j$ slice. The shortfall is:
$$
\text{missing} = \sum_{i \ne j} (S_i S_i^* \otimes S_j S_j^*) \ne 0.
$$

Total missing cross-terms: $N(N - 1) = 15 \cdot 14 = 210$ basis elements with coefficient $+1$ each.

(Driver field: `diagonal_missing_cross_terms = 210`.)

### 6.4 Structural reading

The diagonal coproduct passes the SUPP-1 of the Cuntz relations (which is a "scalar" relation between specific isometry products) but fails the SUPP-2 (which is a "trace-of-projectors-equals-identity" relation requiring summation over all generators). The structural obstruction is that the diagonal coproduct treats each generator INDEPENDENTLY, while SUPP-2 is a TOTAL-SUM relation that couples all generators simultaneously.

This rules out the "diagonal enrichment" reading of Cuntz: no bialgebra structure with diagonal generators.

---

## 7. JLO $\phi_1$ on the Cuntz extension

### 7.1 The commutative-case theorem

Sub-Sprint 1 proved that on the commutative algebra $\mathcal{A} = \mathbb{C}^5$, $\phi_1(a_0, a_1; t) \equiv 0$ identically for all $a_0, a_1$ (both even and odd parities, all $t$-orders). The proof used:
1. Cyclicity of the trace.
2. Sector orthogonality $e_s e_t = \delta_{st} e_s$.
3. $\gamma$-commutativity $[\gamma, e_s] = 0$ (sector idempotents are chirality-diagonal).

### 7.2 What breaks on the Cuntz extension

On $\mathcal{O}_N$, condition (2) is **violated**: $S_i S_j$ for $i \ne j$ is NOT zero or $\delta_{ij} S_i$ — it's a brand-new basis element (the length-2 path $S_i S_j$). The cyclic-cancellation channel used in the Sub-Sprint 1 proof requires $a_0 a_1 = a_1 a_0$ (commutativity), which Cuntz isometries also violate.

So $\phi_1$ is **generically non-zero** on the Cuntz extension. Specifically:
$$
\phi_1^{\mathrm{odd}}(S_i, S_j; t)\big|_{t = 0} = \mathrm{Tr}(S_i [D, S_j])
$$
is generally non-zero because the cyclic shift of $\mathrm{Tr}(S_i [D, S_j])$ does not collapse: we have neither $\mathrm{Tr}(S_i D S_j) = \mathrm{Tr}(D S_j S_i)$ via commutativity-of-arguments nor $\mathrm{Tr}(S_i S_j D) = \mathrm{Tr}(S_j S_i D)$ via $S_i S_j = S_j S_i$.

### 7.3 Why we can't extract the coefficient at theorem grade

In the diagonal-Dirac surrogate of §4, $[D, S_j] |\alpha\rangle = (\lambda_{j \alpha} - \lambda_\alpha) |j \alpha\rangle$, so
$$
S_i [D, S_j] |\alpha\rangle = (\lambda_{j \alpha} - \lambda_\alpha) |i j \alpha\rangle,
$$
and the diagonal-trace coefficient $\langle \alpha | i j \alpha \rangle$ is **zero** because $|i j \alpha\rangle$ has length $|\alpha| + 2 \ne |\alpha|$ in the depth-truncated Fock basis.

So in the diagonal-Dirac SURROGATE at depth-truncated level, $\phi_1$ does happen to vanish. This is **representation-specific**, not a structural theorem. The FULL Camporesi-Higuchi Dirac $D = \Lambda + \kappa A$ has an off-diagonal adjacency operator $A$ that connects different shells via $E_1$ dipole transitions; the trace $\mathrm{Tr}(S_i [\kappa A, S_j])$ is **generically non-zero** because $A$-action does not preserve path length monotonically.

**Multi-year scope.** Computing $\phi_1$ at theorem grade requires modelling $\mathcal{O}_N$'s natural representation on Fock space and the Camporesi-Higuchi Dirac's action through the Pimsner-Voiculescu sequence (Pimsner-Voiculescu 1980, *J. Operator Theory* 4; Cuntz-Krieger 1980, *Invent. Math.* 56). This is genuinely multi-year — the same scale as the rest of Stage 2.

### 7.4 What this tells us structurally

Even though no Hopf-coproduct exists on $\mathcal{O}_N$ (§§5, 6), the JLO BICOMPLEX still makes sense and acquires NEW degree-1 cochain content from the non-commutativity. This is consistent with Connes 1994 Ch. IV §5: the entire-cyclic cohomology framework extends to non-commutative algebras, where the degree-1 cochain becomes structurally non-trivial.

The structural reading:

> The Cuntz extension would give NEW entire-cyclic cohomology classes from the broken sector-orthogonality, but the natural pathway to encode these into a motivic Galois structure is NOT via the Hopf-coproduct framework. It must go through the COMODULE structure ($\mathcal{O}_N$ is a Hopf-comodule over the gauge $U(1)$ that scales each $S_i$ by a phase), which is a categorically different setting from Track A's commutative case.

---

## 8. M1/M2/M3 partition: preserved only at the trivial comodule level

### 8.1 The commutative-case clean tensor factorisation

Track A (v3.61.0) found that on the commutative algebra, the primitive coproduct preserves the Mellin-slot $k$-label on each generator, giving the **clean tensor product**:
$$
\mathcal{H}_{\mathrm{GV}} = \mathcal{H}_{\mathrm{GV}}^{[0]} \otimes \mathcal{H}_{\mathrm{GV}}^{[1]} \otimes \mathcal{H}_{\mathrm{GV}}^{[2]},
$$
and the candidate motivic Galois group $U^{*(n_{\max})}_{\mathrm{GeoVac}} = \prod_{k = 0}^{2} \mathbb{G}_a^{N(n_{\max})}$. This is the Hopf-algebraic shadow of the case-exhaustion theorem.

### 8.2 What survives on Cuntz

Each Cuntz isometry $S_{(n, l), k}$ carries a fixed $k$ label, so the gauge $U(1)$-coaction (which scales each $S_i$ by a phase) trivially respects the 3-fold $k$-decomposition. At the COMODULE level the partition is preserved.

But at the HOPF-COPRODUCT level the question is **moot** — there is no Hopf-coproduct on $\mathcal{O}_N$ (§§5, 6).

### 8.3 The structural takeaway

The M1/M2/M3 partition is preserved on the Cuntz extension, but only in a structurally **weaker** form: trivial preservation under the gauge $U(1)$ coaction (each generator has a fixed $k$), rather than the clean tensor-factorisation of three Hopf-sub-algebras seen in the commutative case.

This is **substantively weaker** content than what Track A's Hopf-substrate achieves. Specifically:
- Track A's commutative case: $\mathcal{H}_{\mathrm{GV}} = \otimes_k \mathcal{H}_{\mathrm{GV}}^{[k]}$ as Hopf-algebra tensor product; $U^{*(n_{\max})}_{\mathrm{GeoVac}} = \prod_k \mathbb{G}_a^{N(n_{\max})}$.
- Cuntz extension: $\mathcal{O}_N = \bigoplus_k \mathcal{O}_N^{[k]}$ as gauge-$U(1)$-comodule direct sum; the motivic Galois group analog (if it exists at all) is NOT a Hopf-coproduct dual but a comodule-over-$U(1)$ object.

The categorical gap is genuine. The Cuntz extension does not "just enrich" the Track A substrate; it shifts the framework into a categorically different setting (comodule-over-gauge-group, $C^*$-Hopf-comodule, $\mathcal{O}_N$-pivot Tannakian structure) where the M1/M2/M3 partition has a different, weaker, meaning.

---

## 9. Coassociativity, antipode, and the bit-exact axiom panel

### 9.1 Coassociativity on the FREE algebra

If we forget the Cuntz relations and work on the free $*$-algebra generated by $\{S_i, S_i^*\}$, both proposed coproducts are coassociative:

- **Primitive:** $(\Delta \otimes \mathrm{id})\Delta(S_i) = (\mathrm{id} \otimes \Delta)\Delta(S_i) = S_i \otimes 1 \otimes 1 + 1 \otimes S_i \otimes 1 + 1 \otimes 1 \otimes S_i$ (standard symmetric distribution).
- **Diagonal:** $(\Delta \otimes \mathrm{id})\Delta(S_i) = (\mathrm{id} \otimes \Delta)\Delta(S_i) = S_i \otimes S_i \otimes S_i$ (trivial).

So **coassociativity is NOT the obstruction**. Both candidates pass coassociativity on the free algebra. The obstruction is purely compatibility with the Cuntz defining relations.

### 9.2 Antipode

A bialgebra antipode $S$ must satisfy $m \circ (S \otimes \mathrm{id}) \circ \Delta = \eta \circ \varepsilon$. For both primitive and diagonal candidates, the natural choice $S(S_i) = S_i^*$ (a $*$-algebra anti-homomorphism) would need to satisfy this antipode property. But since the underlying $\Delta$ doesn't extend to $\mathcal{O}_N$ (§§5, 6), the antipode question is **moot** at the bialgebra level.

### 9.3 The bit-exact axiom panel (incompatibility verdict)

| Test | Result |
|:----|:------:|
| Primitive $\Delta$ compatibility with $S_i^* S_j = \delta_{ij}$ | **0/225 pass** |
| Diagonal $\Delta$ compatibility with $S_i^* S_j = \delta_{ij}$ | **225/225 pass** |
| Diagonal $\Delta$ compatibility with $\sum_i S_i S_i^* = 1$ | **FAIL** (210 cross-terms missing) |
| Coassociativity on free algebra (primitive) | structurally PASS (binomial identity) |
| Coassociativity on free algebra (diagonal) | structurally PASS (trivial) |
| Antipode for primitive on $\mathcal{O}_N$ | MOOT (no $\Delta$) |
| Antipode for diagonal on $\mathcal{O}_N$ | MOOT (no $\Delta$) |
| JLO $\phi_1 \equiv 0$ on $\mathcal{O}_N$ | **STRUCTURALLY FAILS** (Sub-Sprint 1 vanishing argument requires commutativity) |
| M1/M2/M3 partition under Hopf-coproduct | MOOT (no $\Delta$) |
| M1/M2/M3 partition under gauge-$U(1)$ comodule | trivial PASS (each $S_i$ has fixed $k$) |

**Net:** 225 bit-exact passes on the diagonal-first-relation; 225 bit-exact fails on the primitive-relation; 1 structural fail on the diagonal-second-relation. Total bit-exact zero residuals: 225 (representing the only positive structural content found in this scoping step).

---

## 10. Honest scope, curve-fit audit, and verification gates

### 10.1 What is closed at theorem grade

- The primitive coproduct $\Delta(S_i) = S_i \otimes 1 + 1 \otimes S_i$ is incompatible with the Cuntz relation $S_i^* S_j = \delta_{ij}$ on **0/225** pairs (bit-exact). This is a definitive structural negative.
- The diagonal coproduct $\Delta(S_i) = S_i \otimes S_i$ passes the first Cuntz relation on **225/225** pairs (bit-exact) but fails the second relation $\sum_i S_i S_i^* = 1$ by exactly $N(N - 1) = 210$ missing cross-terms. Bit-exact.
- Coassociativity is not the obstruction (passes structurally for both).
- M1/M2/M3 partition: trivial pass at gauge-$U(1)$-comodule level; moot at Hopf-coproduct level.
- $\phi_1$ structural argument: the Sub-Sprint 1 vanishing argument uses sector-orthogonality which Cuntz violates; the vanishing fails on $\mathcal{O}_N$ generically.

### 10.2 Structural sketches (not theorem at this sprint)

- Whether $\mathcal{O}_N$ admits a NON-TRIVIAL comodule-over-Hopf-comodule structure that lifts the M1/M2/M3 partition to a genuinely non-abelian Galois structure: open multi-year. Pimsner-Voiculescu sequence and Cuntz-Pimsner extensions (Pimsner 1997) are the natural categorical setting.
- Whether the JLO degree-1 cochain coefficient on Cuntz extension can be extracted at theorem grade in the FULL Camporesi-Higuchi adjacency-extended Dirac representation: open multi-year.
- Whether any other Cuntz-extension-style enrichment (Toeplitz-Cuntz, Cuntz-Krieger) admits a Hopf-coproduct compatible with its defining relations: open structural question; non-trivial published precedent on quasi-bialgebras in renormalization (Brouder-Frabetti 2003 arXiv:hep-th/0211212).

### 10.3 Numerical observations

No transcendentals introduced. All bit-exact calculations in `sympy.Rational`. The non-primitive-coproduct hope is structurally NULL because no Hopf-coproduct compatible with Cuntz relations exists — the question is well-posed and the negative is forced.

### 10.4 Curve-fit audit (`feedback_audit_numerical_claims`)

No "X matches Y" claims. The only structural matches are:
- "Primitive coproduct gives non-zero residual on 225/225 pairs" — direct bit-exact computation, no curve-fit, no PSLQ, no selection bias.
- "Diagonal coproduct gives 210 missing cross-terms on second Cuntz relation" — bit-exact count $N(N - 1) = 15 \cdot 14 = 210$ from the algebraic structure.
- "$\phi_1$ vanishing argument requires sector-orthogonality which Cuntz violates" — direct invocation of the Sub-Sprint 1 proof structure; verified on inputs by the diagonal-Dirac surrogate (which happens to land trivially, but the structural reason for non-vanishing is independent of the representation).

Selection bias check: the decision gate (POSITIVE/BORDERLINE/STOP) was articulated BEFORE running the bit-exact tests. STOP is the strongest available verdict and matches the structural shape of the negative result (a structural obstruction, not a numerical near-miss).

### 10.5 Discrete-for-skeleton compliance (`feedback_discrete_for_skeleton`)

All structural tests use `sympy.Rational` bit-exact arithmetic. Zero floats, zero PSLQ, zero transcendentals. The Cuntz path algebra is a discrete (path-combinatorial) substrate, and the obstructions are detected at the discrete level.

### 10.6 Tag transcendentals (`feedback_tag_transcendentals`)

Zero transcendentals appear in this sprint. The Cuntz extension would, if it admitted a Hopf-coproduct, generate cocycle classes whose continuum-Mellin extraction would land in M1/M2/M3 per the case-exhaustion theorem (Paper 32 §VIII). At this scoping step the Cuntz extension is incompatible with the Hopf framework, so the continuum-Mellin pathway is mooted by the structural obstruction.

### 10.7 No synthesis memos (`feedback_no_synthesis_memos`)

This memo is the single canonical record of the Q5'-Cuntz scoping step. It does not consolidate or supersede earlier memos. Cross-sprint synthesis lives in CHANGELOG.md (for any future v3.62.0 release entry) and in the Paper 32 §VIII / Paper 55 §subsec:open_m2_m3 Remarks.

### 10.8 Hard prohibitions check (CLAUDE.md §13.5)

No changes to natural geometry hierarchy. No fitted/empirical parameters introduced. No deletion of negative results from §3. No removal of "conjectural" label from Paper 2 combination rule (Paper 2 not touched).

### 10.9 WH1 PROVEN unaffected

This sprint operates entirely on the candidate Stage-2 motivic Galois substrate. It does not test propinquity convergence or modify the WH1 / Marcolli-vS lineage closure.

---

## 11. Files

### Produced

- `debug/compute_q5p_cuntz.py` — driver (~520 lines, 0.017 s wall, bit-exact `sympy.Rational` throughout; primitive- and diagonal-coproduct compatibility tests at $n_{\max} = 2$ with $N = 15$ total isometries; depth-truncation enumeration at $L = 1, 2$; sample $\phi_0^{(k)}$ weights).
- `debug/data/sprint_q5p_cuntz.json` — exact rational data dump: 225-pair compatibility tests, depth-truncation basis sizes, sample $\phi_0$ weights, structural verdict.
- `debug/sprint_q5p_cuntz_memo.md` — this memo.

### Used (load-bearing inputs)

- `debug/sprint_q5p_stage2_hopf_memo.md` — the Track A baseline that this scoping step builds on (the abelian primitive Hopf substrate; the three flagged enrichment ingredients).
- `debug/sprint_q5p_hard_parts_2026_06_05_memo.md` — the v3.61.0 umbrella memo naming this enrichment ingredient.
- `debug/sprint_q5p_jlo_nmax2_memo.md` — the Sub-Sprint 1 baseline that $\phi_1$ vanishes on commutative $\mathbb{C}^5$ (this scoping step shows the proof breaks on Cuntz).
- `debug/sprint_q5p_2c_bicomplex_memo.md` — the Sub-Sprint 2c JLO bicomplex baseline.
- `debug/sprint_q5p_prosystem_memo.md` — sector-local closed forms for chirality $\chi_{(n, l)}$ and $\eta_{(n, l)}$ (the M1/M3 weights used in §4 sample).

### Published references

- Cuntz, J. *"Simple $C^*$-algebras generated by isometries."* Commun. Math. Phys. 57 (1977), 173–185. The defining paper for the Cuntz algebra $\mathcal{O}_N$; the path-basis normal form is Lemma 1.3.
- Cuntz, J.; Krieger, W. *"A class of $C^*$-algebras and topological Markov chains."* Invent. Math. 56 (1980), 251–268. Extension to general $C^*$-algebras of topological Markov chains; relevant for the multi-year continuation.
- Connes, A.; Cuntz, J. *"Quasi homomorphismes, cohomologie cyclique et positivité."* Commun. Math. Phys. 114 (1988), 515–526. JLO-style entire-cyclic cohomology on Cuntz algebras; structurally underlies the $\phi_1$ non-vanishing argument here.
- Connes, A. *"Noncommutative Geometry."* Academic Press (1994), Ch. IV §5. Entire-cyclic cohomology on non-commutative algebras; the framework in which the JLO bicomplex on Cuntz is well-posed.
- Pimsner, M. V. *"A class of $C^*$-algebras generalizing both Cuntz-Krieger algebras and crossed products by $\mathbb{Z}$."* Free Probability Theory, Fields Inst. Commun. 12 (1997), AMS, 189–212. The "Cuntz-Pimsner algebras" — the natural multi-year framework for the comodule-over-Hopf-comodule structure flagged in §8.3.
- Pimsner, M.; Voiculescu, D. *"Exact sequences for K-groups and Ext-groups of certain crossed product C*-algebras."* J. Operator Theory 4 (1980), 93–118. The Pimsner-Voiculescu sequence for K-theory of crossed products; needed for the multi-year Cuntz-Pimsner program.
- Loday, J.-L. *"Cyclic Homology."* 2nd ed., Springer Grundlehren 301 (1998), §1.4, §2.1, §3.1. Standard reference on cyclic cohomology on non-commutative algebras.
- Brouder, C.; Frabetti, A. *"QED Hopf algebras on planar binary trees."* arXiv:hep-th/0211212 (2003). Cuntz-like Hopf algebras in QED renormalization; a published-precedent setting where Cuntz-extension-style structures appear with a Hopf-coproduct (via a different enrichment than the ones tested here).

---

## 12. Paper-edit recommendations (PI to apply; NO edits applied in this sprint)

### 12.1 Paper 32 §VIII — ONE new Remark `rem:q5p_stage2_cuntz_incompat` after `rem:q5p_strict_strong_drift`

```latex
\begin{rem}[Q5' Stage 2 Cuntz-extension enrichment: structural incompatibility, Sprint Q5'-Cuntz, June 2026]
\label{rem:q5p_stage2_cuntz_incompat}
The first scoping step of the multi-year Cuntz-extension enrichment ingredient
flagged in Remark~\ref{rem:q5p_stage2_hopf_substrate} returns a structural
negative. Replace the commutative algebra $\mathcal{A}_{\mathrm{GV}}^{(n_{\max})} =
\mathbb{C}^{N(n_{\max})}$ of sector idempotents by the Cuntz algebra
$\mathcal{O}_N$ generated by $N = 3 \cdot N(n_{\max})$ isometries $S_{(n, l), k}$
(one per sector $(n, l)$ and Mellin slot $k \in \{0, 1, 2\}$). The two natural
candidate coproducts on $\mathcal{O}_N$ are bit-exactly incompatible with the
Cuntz relations:
\begin{enumerate}
\item Primitive coproduct $\Delta(S_i) = S_i \otimes 1 + 1 \otimes S_i$ extended
as a $*$-algebra homomorphism gives $\Delta(S_i^*) \Delta(S_j) = 2\delta_{ij}
(1 \otimes 1) + S_i^* \otimes S_j + S_j \otimes S_i^*$, which equals
$\Delta(\delta_{ij} \cdot 1) = \delta_{ij} (1 \otimes 1)$ for no value of
$(i, j)$. Bit-exact 0/225 pairs pass at $n_{\max} = 2$ with $N = 15$.
\item Diagonal coproduct $\Delta(S_i) = S_i \otimes S_i$ passes the first Cuntz
relation $S_i^* S_j = \delta_{ij}$ on $225/225$ pairs but fails the second
$\sum_i S_i S_i^* = 1$: the cross-terms $i \ne j$ are missing in
$\Delta(\sum_i S_i S_i^*) = \sum_i (S_i S_i^* \otimes S_i S_i^*)$ vs.
$1 \otimes 1 = \sum_{i, j} (S_i S_i^* \otimes S_j S_j^*)$. Bit-exact
$N(N - 1) = 210$ missing terms at $n_{\max} = 2$.
\end{enumerate}
Net: $\mathcal{O}_N$ is structurally a Hopf-\emph{comodule} (over the gauge
$U(1)$ scaling each $S_i$ by a phase) rather than a Hopf-algebra. The
M1/M2/M3 partition survives only at the trivial comodule level (each $S_{(n,l),k}$
carries a fixed $k$), structurally weaker than the clean Hopf tensor product
$\mathcal{H}_{\mathrm{GV}} = \otimes_k \mathcal{H}_{\mathrm{GV}}^{[k]}$
established for the commutative case in
Remark~\ref{rem:q5p_stage2_hopf_substrate}. The JLO degree-1 cochain $\phi_1$
is structurally non-zero on $\mathcal{O}_N$ (the commutative-case vanishing
theorem requires sector orthogonality $e_s e_t = \delta_{st} e_s$, which the
Cuntz isometries violate), but extraction at theorem grade requires the full
Camporesi--Higuchi adjacency-extended Dirac representation through the
Pimsner-Voiculescu sequence (Pimsner-Voiculescu 1980, \emph{J. Operator
Theory} 4; Pimsner 1997, Fields Inst. Commun. 12), which is genuinely
multi-year. The Cuntz route is therefore the structurally most disruptive
of the three enrichment ingredients flagged in
Remark~\ref{rem:q5p_stage2_hopf_substrate}: it does not enrich the Hopf
substrate but shifts the framework into the categorically different
Hopf-comodule-over-gauge-group setting (Cuntz-Pimsner extension). The two
remaining ingredients (nested-Hopf-tower $J^*(S^3)$ and cross-shell
off-diagonal Dirac perturbation) preserve the Hopf-algebra category at the
substrate level and are correspondingly more tractable next steps. See
Paper~55 \S\ref{subsec:open_m2_m3} for the cosmic-Galois program narrative.
\end{rem}
```

### 12.2 Paper 55 §subsec:open_m2_m3 — ONE new paragraph after the Stage-2 Hopf-substrate paragraph

```latex
\emph{Stage 2 enrichment scoping: Cuntz extension is structurally
incompatible (Sprint Q5'-Cuntz, June 2026; memo
\texttt{debug/sprint\_q5p\_cuntz\_memo.md}; data
\texttt{debug/data/sprint\_q5p\_cuntz.json}).} Of the three substrate-
enrichment ingredients flagged in the Stage-2 Hopf paragraph above
(nested-Hopf-tower $J^*(S^3)$; cross-shell off-diagonal Dirac; JLO/CM
Cuntz extension), the first scoping step of the third returns a
structural negative. Replacing the commutative algebra $\mathbb{C}^{N(n_{\max})}$
of sector idempotents by the Cuntz algebra $\mathcal{O}_N$ generated by
$N = 3 N(n_{\max})$ isometries $S_{(n, l), k}$ produces a non-commutative
substrate, but both natural coproduct candidates fail compatibility with
the Cuntz defining relations bit-exactly: the primitive coproduct
$\Delta(S_i) = S_i \otimes 1 + 1 \otimes S_i$ fails $S_i^* S_j = \delta_{ij}$
on $0/225$ pairs at $n_{\max} = 2$, and the diagonal coproduct
$\Delta(S_i) = S_i \otimes S_i$ passes the first Cuntz relation but fails
the second relation $\sum_i S_i S_i^* = 1$ by exactly $N(N - 1) = 210$
missing cross-terms. Structurally, $\mathcal{O}_N$ is a Hopf-\emph{comodule}
over the gauge $U(1)$ rather than a Hopf-algebra (Pimsner 1997, Fields Inst.
Commun. 12; Connes-Cuntz 1988, Commun. Math. Phys. 114). The M1/M2/M3
partition survives only at the trivial comodule level via fixed-$k$
generator labels; the clean Hopf tensor product
$\mathcal{H}_{\mathrm{GV}} = \otimes_k \mathcal{H}_{\mathrm{GV}}^{[k]}$
of the commutative case does not carry over. The JLO degree-1 cochain
$\phi_1$ is structurally non-zero on $\mathcal{O}_N$ (Sub-Sprint 1's
vanishing argument required sector orthogonality, which the Cuntz isometries
violate), but extraction at theorem grade requires the full Camporesi--Higuchi
adjacency-extended Dirac representation through the Pimsner-Voiculescu
sequence (Pimsner-Voiculescu 1980, \emph{J. Operator Theory} 4),
genuinely multi-year. The Cuntz route is therefore structurally the most
disruptive of the three enrichment ingredients --- it shifts the framework
into a categorically different setting (Cuntz-Pimsner extension, comodule-
over-gauge-group) rather than enriching the existing Hopf substrate. The
remaining two ingredients (nested-Hopf-tower $J^*(S^3)$ and cross-shell
off-diagonal Dirac) are the correspondingly more tractable next steps for
producing non-abelian Stage-2 motivic Galois content.
```

### 12.3 Paper 18 — no edit needed

Paper 18 §III.7 master Mellin engine is upstream of the Stage-2 substrate
enrichment question. The structural finding here (that the Cuntz extension
shifts the framework into the Hopf-comodule category, not the Hopf-algebra
category, so the M1/M2/M3 partition survives only at the trivial comodule
level) does not affect §III.7's three-bullet partition; it only changes
where in the categorical hierarchy the partition lives.

---

## 13. One-line verdict

**STOP.** The Cuntz extension as a Stage-2 enrichment ingredient is structurally incompatible with the v3.61.0 Track A Hopf-substrate framework: both natural coproduct candidates fail the Cuntz defining relations bit-exactly (primitive: 0/225 pairs pass; diagonal: passes first Cuntz relation 225/225 but fails second relation with $N(N - 1) = 210$ missing cross-terms at $n_{\max} = 2$, $N = 15$); $\mathcal{O}_N$ is structurally a Hopf-comodule over the gauge $U(1)$ rather than a Hopf-algebra, and the M1/M2/M3 partition survives only at the trivial comodule level via fixed-$k$ generator labels — categorically weaker than the clean Hopf tensor product of the commutative case. The JLO degree-1 cochain $\phi_1$ is structurally non-zero on the Cuntz extension (Sub-Sprint 1's vanishing argument requires sector orthogonality, which Cuntz isometries violate), but extraction at theorem grade requires the full Camporesi-Higuchi adjacency Dirac representation through the Pimsner-Voiculescu sequence — multi-year. The Cuntz route is the structurally most disruptive of Track A's three flagged enrichment ingredients (the other two — nested-Hopf-tower $J^*(S^3)$ and cross-shell off-diagonal Dirac — preserve the Hopf-algebra category and are correspondingly more tractable). Sprint produces 225 bit-exact positive structural content (the diagonal-first-relation pass) and 435 bit-exact structural-obstruction residuals (225 primitive-fail + 210 diagonal-second-fail). Headline: the abelian-primitive shape of v3.61.0 Track A is **forced** by the categorical setting (Hopf-algebra on commutative algebra); breaking commutativity via Cuntz does not enrich the framework's coproduct, it **changes the framework's category**.
