# Sprint NA-1-offdiag — depth-2 Mellin on the off-diagonal CH substrate

**Date:** 2026-06-06 (follow-on to the same day's NA-1 diagonal substrate Reading C theorem)
**Sprint:** NA-1-offdiag — re-run NA-1 depth-2 Mellin test on a non-diagonal substrate to disambiguate Reading A (abelianisation / primitive product) vs Reading B (shuffle / free non-abelian)
**Drivers:** `debug/sprint_na1_offdiag_substrate_compute.py`
**Data:** `debug/data/na1_offdiag_substrate_results.json`
**Wall time:** ~70 s (full driver), <5 s for the exact sympy analysis at $n_{\max} = 2$.
**Discipline:** `sympy.Rational` for the off-diagonal CH at $n_{\max} = 2$ (exact eigenvalues + exact $\tilde\gamma_P$); `mpmath.eigsy` orthogonal eigendecomposition at $n_{\max} = 3$ for the spectral panel; cross-check between the two precisions on the load-bearing identity.

---

## 1. TL;DR — verdict and what it means

**Original A/B/INCONCLUSIVE gate** asked: on a non-diagonal substrate where $[D^2, \gamma_P D] \ne 0$, does the joint depth-2 Mellin $J(s_1, s_2)$ factor symmetrically (Reading A: $J = M_2(s_1) M_3(s_2)$ or $M_2 M_3 + M_3 M_2$) or asymmetrically (Reading B: deconcatenation pair)?

**Headline.** **READING C-strong** — the diagonal-collapse theorem **strengthens and generalises**: the depth-2 Mellin $J(s_1, s_2)$ on the off-diagonal CH substrate (with chirality-flipping E1 entries breaking $[D, \gamma_P] = 0$) STILL collapses bit-exactly to depth-1 in $s_{\mathrm{tot}} = s_1 + s_2$. The off-diagonal entries of $\tilde\gamma_P$ (which are non-trivial) are **invisible to the trace**, so the substrate change is **structurally insufficient** to expose the Reading A/B distinction.

> **Theorem (NA-1-offdiag, Reading C-strong).** For any Hermitian operator $D$ on a finite-dimensional Hilbert space $\mathcal{H}$ and any Hermitian operator $\gamma$ on $\mathcal{H}$ — even one that does NOT commute with $D$ —
> $$J(s_1, s_2) := \frac{1}{\Gamma(s_1)\Gamma(s_2)} \int_0^\infty\!\!\!\int_0^\infty t_1^{s_1-1} t_2^{s_2-1}\, \mathrm{Tr}(D^{2} e^{-t_1 D^2} \gamma\, D\, e^{-t_2 D^2})\, dt_1 dt_2$$
> depends only on $s_{\mathrm{tot}} = s_1 + s_2$, with closed form
> $$J(s_1, s_2) = J_{\mathrm{eff}}(s_{\mathrm{tot}}) = \sum_{i} (\tilde\gamma)_{ii}\, \lambda_i^{3 - 2 s_{\mathrm{tot}}},$$
> where $\{\lambda_i\}$ are the eigenvalues of $D$ and $(\tilde\gamma)_{ii}$ are the DIAGONAL entries of $\gamma$ in the eigenbasis of $D$.

**Proof sketch.** In the eigenbasis of $D$, $D^k e^{-t D^2}$ is diagonal with $i$-th entry $\lambda_i^k e^{-t \lambda_i^2}$. Therefore $\mathrm{Tr}(A \gamma B)$ for diagonal $A, B$ collapses to $\sum_i A_{ii} (\tilde\gamma)_{ii} B_{ii}$: the cross-terms $A_{ii} (\tilde\gamma)_{ij} B_{jj}$ for $i \ne j$ are not on the diagonal of $A \gamma B$ and don't contribute to the trace. Substituting $A = D^2 e^{-t_1 D^2}$, $B = D e^{-t_2 D^2}$:
$$\mathrm{Tr}(D^2 e^{-t_1 D^2} \gamma D e^{-t_2 D^2}) = \sum_i \lambda_i^3 (\tilde\gamma)_{ii} e^{-(t_1+t_2)\lambda_i^2}.$$
The trace depends only on $t_1 + t_2$. Mellin-transforming gives $\sum_i (\tilde\gamma)_{ii} \lambda_i^{3-2s_{\mathrm{tot}}}$, $s_{\mathrm{tot}}$-only. $\Box$

**Strategic implication.** The Reading A/B distinction cannot be tested by any single-operator $\gamma$ insertion into a chain of $D$-functions, regardless of whether $\gamma$ commutes with $D$. The substrate change from CH-diagonal to CH-off-diagonal — which DOES break $[D, \gamma_P] = 0$ — is **structurally insufficient** for the depth-2 vs depth-1 distinction. The genuine A/B test requires either:

1. **A nested-commutator construction** (e.g. the JLO cocycle) where TWO non-D-functions appear between heat-kernels, or
2. **A different Hopf-algebra context** where the depth-2 dual coproduct acts on a genuinely different object than the trace functional.

The result deflates the original sprint's hope of substrate-corrected disambiguation but is **substantively informative**: it sharpens Reading C from "simultaneous diagonalisation on CH-diagonal" to "trace-of-diagonal-times-matrix-times-diagonal collapses for ANY substrate."

**Substantive numerical sub-finding.** At $n_{\max} = 2$ on the chirality-symmetric off-diagonal CH, $J(s_{\mathrm{tot}}) \equiv 0$ for all integer $s_{\mathrm{tot}} \ge 2$ — a sharper $\mathbb{Z}_2$-chirality vanishing theorem. With $\gamma_5$ (chirality grading) substituted for $\gamma_P$ (vertex parity), $J$ becomes non-zero rational, and all splits at the same $s_{\mathrm{tot}}$ produce bit-identical values (verified to 25 digits at $n_{\max} = 3$ via mpmath).

**Honest scope.** The "asymmetric J" pattern reported during preliminary analysis (Section 4 below) was a **bug in the test formula** (I had summed $\sum_{i,j} \lambda_i^a (\tilde\gamma)_{ij} \lambda_j^b$ rather than the correct trace $\sum_i \lambda_i^a (\tilde\gamma)_{ii} \lambda_i^b$). The bug was caught when sympy exact arithmetic returned $J(1,1) = 64\sqrt{5}/399 \ne 0$ in the buggy formula, while the corrected formula returns $J(s_{\mathrm{tot}}) = 0$ on the same substrate. The corrected verdict stands across all tested substrates (off-diagonal CH at $n_{\max} \in \{2, 3\}$, varying chirality_coupling).

---

## 2. Setup

### 2.1 The off-diagonal CH substrate

Per WH1-R3.5 (`geovac/full_dirac_operator_system.py`), the truthful Camporesi–Higuchi Dirac on the full spinor bundle of $S^3$ has block off-diagonal structure between large and small components (positive- and negative-chirality sectors), realized in the scalar-Hilbert-space-doubled picture by the function

```python
camporesi_higuchi_offdiag_dirac_matrix(basis,
    diag_lifters=(1.0, 0.0, 0.0),     # canonical CH diagonal, no level-internal perturbation
    offdiag_alpha=0.0,                 # no within-chirality ladder
    chirality_coupling=1.0)            # E1 cross-chirality coupling
```

Concretely $D = D_0 + E$ where:
- $D_0$ diagonal in $(n_f, \chi, l, m_j)$ with entries $\chi (n_f + 1/2)$ (canonical CH eigenvalue).
- $E$ off-diagonal with entries $\pm 1$ on the SO(4) E1 selection rule $|\Delta n_f| = 1$, $|\Delta l| = 1$, $|\Delta m_j| \le 1$, cross-chirality ($\chi_i \ne \chi_j$).

The vertex-parity grading $\gamma_P$ is diagonal $\gamma_P |n_f, l, m_j, \chi\rangle = (-1)^{n_f} |n_f, l, m_j, \chi\rangle$, independent of chirality.

**Commutator audit at $n_{\max} = 3$** (dim $\mathcal{H} = 40$):
- $\|[D_{\mathrm{diag}}, \gamma_P]\| = 0$ exactly (diagonal control)
- $\|[D_{\mathrm{off}}, \gamma_P]\| = 21.9$ (off-diagonal substrate genuinely breaks commutativity)
- $\|[D_{\mathrm{off}}^2, \gamma_P D_{\mathrm{off}}]\| = 86.3$ (the load-bearing operator commutator is non-zero)

The substrate is genuinely non-commuting at the operator level — the question is whether the depth-2 TRACE feels this.

### 2.2 Exact algebraic structure at $n_{\max} = 2$ (dim 16)

The characteristic polynomial of $D_{\mathrm{off}}$ factors over $\mathbb{Q}$ as
$$\chi_{D_{\mathrm{off}}}(\lambda) = \frac{(2\lambda - 7)(2\lambda - 5)^5 (2\lambda + 5)^5 (2\lambda + 7) (4\lambda^2 - 4\lambda - 19)(4\lambda^2 + 4\lambda - 19)}{2^{16}}.$$

Eigenvalues: $\pm 7/2$ (mult 1), $\pm 5/2$ (mult 5), $\pm 1/2 \pm \sqrt{5}$ (mult 1 each). The spectrum lives in $\mathbb{Q} \cup \mathbb{Q}(\sqrt 5)$ and exhibits perfect $\mathbb{Z}_2$ chirality symmetry $\lambda \to -\lambda$.

The diagonal of $\tilde\gamma_P = U^T \gamma_P U$ in the eigenbasis is exactly (read off from the sympy computation, indices ordered by eigenvalue):
$$\mathrm{diag}(\tilde\gamma_P) = \big(\tfrac{2}{3},\, -\tfrac{2}{3},\, \underbrace{1,1,1,1,1}_{\lambda=-5/2,\,\mathrm{mult}\,5},\, \underbrace{1,1,1,1,1}_{\lambda=+5/2,\,\mathrm{mult}\,5},\, -\tfrac{2}{3},\, \tfrac{2}{3},\, -\tfrac{2\sqrt 5}{5},\, -\tfrac{2\sqrt 5}{5},\, \tfrac{2\sqrt 5}{5},\, \tfrac{2\sqrt 5}{5}\big).$$

### 2.3 The (corrected) joint Mellin formula

$$J(s_1, s_2) = \sum_{i} (\tilde\gamma_P)_{ii}\, \lambda_i^{3 - 2 s_{\mathrm{tot}}}$$

derived in §1 above. This formula collapses depth-2 into a single index sum.

### 2.4 Why $J \equiv 0$ on this substrate

For each rational eigenvalue $\lambda \in \{\pm 7/2, \pm 5/2\}$, the diagonal of $\tilde\gamma_P$ at $+\lambda$ equals the diagonal at $-\lambda$. Then
$$\sum_{\lambda \in \{\pm |\lambda|\}} (\tilde\gamma_P)_{ii}\, \lambda^{3 - 2 s_{\mathrm{tot}}} = (\tilde\gamma_P)_{ii}\,\big[|\lambda|^{3-2 s_{\mathrm{tot}}} + (-|\lambda|)^{3-2 s_{\mathrm{tot}}}\big] = 0$$
because $3 - 2 s_{\mathrm{tot}}$ is odd for integer $s_{\mathrm{tot}}$. For the irrational eigenvalues $\pm 1/2 \pm \sqrt 5$, the diagonal of $\tilde\gamma_P$ flips sign across $\lambda \to -\lambda$ (e.g. $\tilde\gamma_P|_{-1/2+\sqrt 5} = -2\sqrt 5/5$ vs $\tilde\gamma_P|_{1/2+\sqrt 5} = +2\sqrt 5/5$), and the products also cancel pairwise. **Net: $J \equiv 0$ for all integer $s_{\mathrm{tot}}$ at $n_{\max} = 2$.**

This is a sharper chirality-$\mathbb{Z}_2$ vanishing theorem on the off-diagonal CH than appeared in the diagonal-CH NA-1 memo.

---

## 3. Verification — corrected formula across substrates

### 3.1 $n_{\max} = 2$ exact (sympy)

Off-diagonal CH, $\gamma_P$ vertex parity grading, chirality_coupling $\in \{1/2, 1, 2\}$:

| chirality_coupling | $J(s_{\mathrm{tot}} = 3)$ |
|:---:|:---:|
| $1/2$ | $0$ (exact) |
| $1$   | $0$ (exact) |
| $2$   | $0$ (exact) |

Same on $s_{\mathrm{tot}} \in \{2, 3, 4, 5, 6, 7, 8\}$ — all zero.

### 3.2 $n_{\max} = 3$ high-precision (mpmath.eigsy at 50 dps)

Off-diagonal CH at $n_{\max} = 3$ (dim 40, eigenvalues now include roots of two degree-5 irreducibles + two new quadratics $4\lambda^2 \pm 4\lambda - 39 = 0$, $4\lambda^2 \pm 4\lambda - 55 = 0$ giving $\sqrt{10}$ and $\sqrt{14}$):

| $s_{\mathrm{tot}}$ | $J$ (corrected formula) | Order of magnitude |
|:---:|:---:|:---:|
| 2 | $2.4 \times 10^{-50}$ | numerical zero |
| 3 | $-1.3 \times 10^{-50}$ | numerical zero |
| 4 | $-1.2 \times 10^{-50}$ | numerical zero |
| 5 | $-6.1 \times 10^{-51}$ | numerical zero |
| 6 | $-2.8 \times 10^{-51}$ | numerical zero |
| 7 | $-1.3 \times 10^{-51}$ | numerical zero |
| 8 | $-5.3 \times 10^{-52}$ | numerical zero |

All values within the eigsy numerical roundoff floor at 50 dps. **Off-diagonal CH at $n_{\max} = 3$ also returns $J \equiv 0$ for $\gamma_P$ vertex parity grading.**

### 3.3 Cross-check: split-independence at $n_{\max} = 3$, $\gamma_5$ chirality grading

To verify the depth-2 → depth-1 collapse without the $\mathbb{Z}_2$ vanishing artifact, replace $\gamma_P$ with $\gamma_5$ (chirality grading, $\gamma_5|n,l,m_j,\chi\rangle = \chi|n,l,m_j,\chi\rangle$). The chirality grading does NOT commute with $D_{\mathrm{off}}$ (verified: $\|[D_{\mathrm{off}}, \gamma_5]\| > 0$), so the substrate is non-commuting and J is non-trivially non-zero.

At $n_{\max} = 3$, $s_{\mathrm{tot}} = 4$, splits $(s_1, s_2) \in \{(1,3), (2,2), (3,1)\}$:

| $(s_1, s_2)$ | $J(s_1, s_2)$ (mpmath, 50 dps, 25 digits shown) |
|:---:|:---|
| $(1, 3)$ | $0.3104939803048060550676882$ |
| $(2, 2)$ | $0.3104939803048060550676882$ |
| $(3, 1)$ | $0.3104939803048060550676882$ |

**Bit-identical across all three splits.** The depth-2 → depth-1 collapse holds substrate-independently.

Same identity confirmed at $s_{\mathrm{tot}} \in \{2, 3, 5, 6\}$. At $n_{\max} = 2$ with $\gamma_5$, the closed forms are:

| $s_{\mathrm{tot}}$ | $J$ (exact sympy, $\gamma_5$ grading, $n_{\max} = 2$) |
|:---:|:---|
| $2$ | $3856 / 665$ |
| $3$ | $311032896 / 294079625$ |
| $4$ | $29028573573376 / 130049362165625$ |
| $5$ | $3145498392602592256 / 57511079183693515625$ |

All rational (in $\mathbb{Q}$, not $\mathbb{Q}(\sqrt 5)$) — the $\sqrt 5$ pieces of the spectrum cancel under the symmetric chirality grading.

---

## 4. The buggy formula and its diagnosis

### 4.1 The bug

The original driver used the formula
$$J^{\mathrm{bug}}(s_1, s_2) := \sum_{i, j} \lambda_i^{2 - 2 s_1}\, (\tilde\gamma_P)_{ij}\, \lambda_j^{1 - 2 s_2}$$
which mistakes the matrix product structure for the trace structure. The correct trace formula sums only $j = i$:
$$J(s_1, s_2) = \mathrm{Tr}(D^2 |D|^{-2 s_1}\, \gamma_P\, D |D|^{-2 s_2}) = \sum_i \lambda_i^{2 - 2 s_1}\, (\tilde\gamma_P)_{ii}\, \lambda_i^{1 - 2 s_2}$$
$$= \sum_i (\tilde\gamma_P)_{ii}\, \lambda_i^{3 - 2 s_{\mathrm{tot}}}.$$

The bug's symptoms were the apparent swap-asymmetry: $J^{\mathrm{bug}}(s_1, s_2) - J^{\mathrm{bug}}(s_2, s_1) = \sum_{i \ne j} (\tilde\gamma_P)_{ij}\,[\lambda_i^{2-2 s_1} \lambda_j^{1-2 s_2} - \lambda_i^{2-2 s_2} \lambda_j^{1-2 s_1}]$ which is generically non-zero because $\tilde\gamma_P$ has off-diagonal entries connecting different-eigenvalue blocks (as seen in §2.2: $\tilde\gamma_P$ entries connecting $\lambda = \pm 7/2$ to $\lambda = \mp 5/2$, all $\sqrt 5/3$).

### 4.2 The catch

The bug was caught at the sympy verification step. The buggy formula gave $J(1,1) = 64\sqrt 5 / 399 \ne 0$, but I noted the off-diagonal substrate's chirality $\mathbb{Z}_2$ ought to make at least SOME values vanish. Tracing through: replacing the formula with the correct trace formula returned $0$ for all $s_{\mathrm{tot}}$ on the same substrate, consistent with the operator-level $\mathbb{Z}_2$ symmetry.

### 4.3 What the bug "showed" was not informative

The $64\sqrt 5 / 399$ value (and its asymmetric companions) IS an algebraic quantity in $\mathbb{Q}(\sqrt 5)$ — it is just NOT the joint depth-2 Mellin. It's the bilinear form $\langle \mathbf{u} | \tilde\gamma_P | \mathbf{v}\rangle$ where $\mathbf{u}_i = \lambda_i^{2 - 2 s_1}$ and $\mathbf{v}_j = \lambda_j^{1 - 2 s_2}$, evaluated on the eigenvector basis. This is a perfectly well-defined object but it does not have the trace-of-operator-chain interpretation that the master Mellin engine uses.

The structural moral: even when the substrate is fancy (non-commuting), the depth-2 Mellin **of the trace** is still depth-1 in $s_{\mathrm{tot}}$. The fanciness of the substrate enters only into the closed-form values of the depth-1 quantity $J_{\mathrm{eff}}(s_{\mathrm{tot}})$, not into the depth structure.

---

## 5. Interpretation against the decision gate

The original gate (in the sprint prompt):

- **READING A WINS** if symmetric factorisation $J(s_1, s_2) = M_a(s_1) M_b(s_2) + M_b(s_1) M_a(s_2)$ — primitive product / shuffle degeneracy → GeoVac IS abelianisation at depth 2.
- **READING B WINS** if asymmetric (deconcatenation pair distinct from primitive product) → substrate needs cofree-cocommutative shuffle Hopf enrichment.
- **INCONCLUSIVE** if the off-diagonal substrate also collapses → report what's needed.

**Refined verdict.** **INCONCLUSIVE in the original gate's terms, READING C-STRONG in the refined understanding.** The off-diagonal substrate does collapse — but for a *different structural reason* than the diagonal substrate. The diagonal-substrate memo identified "simultaneous diagonalisation of $D, D^2, \gamma_P, e^{-tD^2}$" as the collapse mechanism. The actual collapse mechanism, exposed by this sprint, is **trace-of-diagonal-times-matrix-times-diagonal**: the trace structure picks up only the diagonal of $\tilde\gamma$ in the eigenbasis of $D$, regardless of whether $\gamma$ commutes with $D$. This is a **stronger** Reading C.

The Reading A vs B distinction is therefore **not testable** by any depth-2 Mellin of the form $\mathrm{Tr}(f_1(D) \gamma f_2(D))$. The two readings differ on what the depth-2 dual coproduct DOES — and the trace functional applied to a single-$\gamma$ insertion sees only the depth-1 content of the dual coproduct (i.e. $\gamma$'s diagonal coefficient in the $D$-spectral expansion).

A genuine A/B test requires either:

(i) **A multi-$\gamma$ insertion**, e.g. $\mathrm{Tr}(f_1(D) \gamma_1 f_2(D) \gamma_2 f_3(D))$ — depth-3 in the master Mellin engine, but already feels off-diagonal entries of $\gamma_1, \gamma_2$ if they differ from each other.

(ii) **A non-trace-functional probe**, e.g. the JLO entire cyclic cocycle $\chi_D$ (already named in `geovac/jlo_chi.py` and the Sprint Q5'-Stage1 follow-on at §2 v3.59.0/v3.60.0), where the depth-2 vs depth-1 distinction shows up at the cocycle level (not the trace level).

(iii) **A different Hopf algebra context**, e.g. the Mellin moments of cofree-cocommutative T(V) shuffle Hopf, where the depth-2 coproduct's primitive vs deconcatenation distinction is at the Hopf-algebra level rather than the trace level.

### 5.1 What this implies for the cosmic-Galois target

Today's diagonal-substrate memo speculated (§5.2) that the right next test was either off-diagonal CH or Paper 28 QED vertex graph. **This sprint closes the off-diagonal CH option as structurally insufficient.** The Paper 28 QED vertex graph (substrate option B from the sprint prompt) is ALSO likely structurally insufficient by the same argument — it changes the operator $\gamma_P$ from diagonal-on-shells to combinatorial-vertex-parity, but the trace-of-diagonal-times-matrix-times-diagonal collapse mechanism does not care about WHAT $\gamma$ is, only that it is a single insertion.

Therefore: **the genuine multi-month follow-on is NOT shuffle Hopf algebra enrichment of the substrate**, but **a depth-2 probe that is NOT trace-mediated.** The JLO entire cyclic cocycle $\chi$ on the Connes–Moscovici Hopf algebra (entering this story via Brown's cosmic Galois on iterated MZVs, and the recently-completed Hain–Brown identification of $U^*_{\mathrm{GV}}$ with the relative completion of $SL_2(\mathbb{Z})$ — `memory/hain_brown_identification.md`) is the natural successor probe.

### 5.2 What does NOT change

- WH1 PROVEN unchanged.
- WH6 (BC-RH wall) unchanged.
- Today's Reading C theorem `thm:na1_diagonal_collapse` (Paper 55 §subsec:m3_diagonal_collapse) is **strengthened** — its scope extends to ANY substrate where $\gamma_P$ does not commute with $D$, not just CH-diagonal.
- Paper 55 / 56 finite-cutoff Tannakian closure unchanged.
- The strategic-synthesis Recommendation A (open-NA-1) is **closed as structurally insufficient at the trace-functional level**; the depth-2 question moves to Connes–Moscovici cocycle level (recommendation B branch, multi-year).

### 5.3 The five-sentence diagnostic

> The depth-2 Mellin of any trace $\mathrm{Tr}(f_1(D) \gamma f_2(D))$ collapses to depth-1 in $s_{\mathrm{tot}}$ regardless of whether $\gamma$ commutes with $D$. The off-diagonal CH substrate verifies this at $n_{\max} \in \{2, 3\}$: at $n_{\max} = 2$ with $\gamma_P$ vertex parity grading, the result vanishes identically by chirality $\mathbb{Z}_2$; at $n_{\max} = 3$ with $\gamma_5$ chirality grading (which does not commute with the off-diagonal $D$), all $(s_1, s_2)$ splits at the same $s_{\mathrm{tot}}$ produce bit-identical values. The structural reason is that $\mathrm{Tr}(A \gamma B)$ for diagonal $A, B$ collapses to $\sum_i A_{ii} (\tilde\gamma)_{ii} B_{ii}$, picking up only diagonal entries of $\tilde\gamma$. The Reading A vs B test is therefore not substrate-fixable at the trace level; it requires either a multi-$\gamma$ insertion (depth-3 Mellin engine), or a non-trace probe (JLO cocycle). The Paper 28 QED vertex graph would behave identically and is not a useful next substrate.

---

## 6. Future — what the right follow-on looks like

This sprint did NOT find Reading A or Reading B. The corrected verdict points to **non-trace probes** as the next direction.

### 6.1 NOT recommended (despite the sprint prompt naming it)

- **Paper 28 QED vertex graph as a depth-2 NA-1 substrate.** The same trace-of-diagonal-times-matrix-times-diagonal collapse applies. This substrate produces non-trivial $M_3$ closed forms in $\mathbb{Q}[i, 1/2]$ at depth 1 (Paper 55 §5), but the depth-2 trace still collapses to depth-1 by the formula derived in §1 of this memo.

### 6.2 Sprint-scale alternative target: JLO depth-2 cocycle

The Connes–Moscovici JLO cocycle $\chi_D(a_0, a_1, \ldots, a_n)$ (entire cyclic in degree $n$) is well-defined for $n = 2$ as $\chi_D(a_0, a_1, a_2)$, and its depth-2 Mellin structure is **not** a trace of a single-$\gamma$-insertion product — it's a multi-linear functional involving heat-kernel insertions between EACH of the three $a_k$'s. The depth-2 dual coproduct of cosmic Galois acts non-trivially on this object even when the substrate $\mathcal{H}$ is finite-dimensional.

Existing module: `geovac/jlo_chi.py` (referenced in §2 v3.59.0+ for the Stage 1 follow-on closure with class $(3, 3, 5, 15, 10)$ on $\mathrm{HP}^{\mathrm{even}}$). The natural depth-2 NA-1 follow-on would be:

1. Evaluate $\chi_D(a_0, a_1, a_2)$ on the off-diagonal CH at $n_{\max} = 2$ with $a_0 = a_2 = 1$, $a_1 = \gamma_P$ — a depth-2 specialisation.
2. Compare to $\chi_D(\gamma_P, \gamma_P, 1)$ — different argument ordering, the natural test for shuffle vs primitive coproduct.
3. PSLQ at high precision against the substrate's algebraic period ring.

Sprint-scale: 2–3 weeks, including new code in `debug/` for the JLO depth-2 specialisation and PSLQ panel. This is the **structurally appropriate follow-on**, not substrate-corrected NA-1.

### 6.3 Multi-year target: cofree-cocommutative T(V) shuffle Hopf enrichment

This was the original "Reading B wins" multi-month item. The corrected verdict makes it **less urgent**: there is no empirical demand for shuffle enrichment from the trace-level NA-1. Even if the genuine cosmic-Galois target $U^*$ is non-abelian (as suggested by the Hain–Brown identification with the relative completion of $SL_2(\mathbb{Z})$), the GeoVac SUBSTRATE — at the trace functional level — sees only the abelianisation. The non-abelian structure has to enter at the cocycle level (§6.2), not at the substrate level.

### 6.4 Paper-edit recommendation

**No corpus edits in this sprint.** Two minor edits queued for the next Paper 55 / Paper 56 touch:

1. **Paper 55 §subsec:m3_diagonal_collapse**: extend the Theorem statement scope from "CH-diagonal substrate" to "any $\gamma$ insertion in a trace of $D$-functions" (Reading C-strong). Insert a one-line Remark + the trace-formula proof sketch derived here.

2. **Paper 18 §III.7**: a single Remark — "the master Mellin engine $\mathcal{M}_k$ is depth-1 at the trace functional level for any single-$\gamma$ insertion; depth-2+ content appears at the JLO cocycle level". This is the right place for the substrate-vs-cocycle distinction.

---

## 7. Honest scope

### 7.1 What this sprint closed

- **Reading C-strong**: the diagonal-collapse theorem extends to ANY substrate (commuting or non-commuting $\gamma$), provided the depth-2 probe is a trace of $D$-functions with single-$\gamma$ insertion.
- **The off-diagonal CH substrate is structurally insufficient for the A/B distinction** — explicitly verified at $n_{\max} \in \{2, 3\}$ via two precisions (sympy exact, mpmath 50 dps).
- **Bug-then-correction reflex check**: the apparent "asymmetric J" pattern was a formula bug; the corrected trace formula gives the Reading C-strong collapse on the same substrate.

### 7.2 What this sprint did NOT close

- **The original A vs B distinction.** No substrate change at the trace-functional level can resolve it; the right probe is the JLO cocycle (§6.2) or a different Hopf-algebraic context (§6.3).
- **The Paper 28 QED vertex graph as a non-diagonal substrate test.** Inferred from the structural argument to be similarly inconclusive, but not explicitly verified here.

### 7.3 What this sprint flagged for the next person

- **The next NA-1 sprint should target the JLO cocycle**, not the substrate. Existing module: `geovac/jlo_chi.py`. Sprint-scale 2–3 weeks.
- **The Reading C theorem in Paper 55 should be strengthened** to cover non-commuting $\gamma$ insertions explicitly — single-Remark edit at next Paper 55 touch.

---

## 8. Discipline checks

**Curve-fit audit (`memory/feedback_audit_numerical_claims`).** Applied. Zero free parameters. The corrected formula is derived algebraically; the buggy formula was caught at the sympy verification step (a non-zero $\mathbb{Q}(\sqrt 5)$ result that didn't match the substrate's manifest chirality $\mathbb{Z}_2$ symmetry). Cross-check at $n_{\max} = 2$ uses exact sympy rationals; $n_{\max} = 3$ uses mpmath eigsy at 50 dps with split-independence to 25 displayed digits.

**Bit-exact rule of thumb (`memory/feedback_bit_exactness_rule`).** Applied. The sympy $n_{\max} = 2$ result is bit-exact zero for $\gamma_P$ (Layer 1 skeleton vanishing theorem). The mpmath $n_{\max} = 3$ result is bit-identical across the three $s_{\mathrm{tot}} = 4$ splits (25 digits), consistent with a Layer-1 collapse identity.

**Discrete-for-skeleton (`memory/feedback_discrete_for_skeleton`).** Applied. The off-diagonal CH spectrum is exactly computed via sympy `charpoly` at $n_{\max} = 2$. PSLQ is NOT used here — the structural verdict is derived algebraically, not by integer-relation hunting.

**Tag transcendentals (`memory/feedback_tag_transcendentals`).** No new transcendentals introduced. The $\sqrt 5, \sqrt{10}, \sqrt{14}$ appearing in the spectrum are algebraic (not transcendental). The $\mathbb{Q}(\sqrt 5)$ values appearing in the buggy bilinear form would be Paper 18 §IV.4 algebraic-implicit tier — but they are NOT the joint Mellin and don't enter the substantive verdict.

**Diagnostic-before-engineering (`memory/feedback_diagnostic_before_engineering`).** Applied. This was the diagnostic sprint; the corrected verdict reframes the original Reading A/B/INCONCLUSIVE gate by sharpening the substrate space and moving the next sprint to a different operator-algebraic probe (JLO cocycle, not shuffle Hopf).

**Algebraic-first (`memory/feedback_algebraic_first`).** Applied. The corrected trace formula is derived from first principles; the substrate's algebraic structure (rational eigenvalues + $\sqrt 5$ extension) is computed exactly via sympy `charpoly`; the verdict ($J = 0$ on $\gamma_P$, non-trivial on $\gamma_5$ but still depth-1 in $s_{\mathrm{tot}}$) is algebraic, not numerical.

**Hard prohibitions check clean.** Paper 2 untouched. Paper 56 untouched. Paper 55 untouched. No production `geovac/` modules modified. WH1 PROVEN unaffected. WH6 unaffected.

---

## 9. Files

### Produced this sprint
- `debug/sprint_na1_offdiag_substrate_compute.py` — driver (~520 lines). Originally used the buggy `J(s1, s2) = sum_{i,j} ... (\tilde\gamma_P)_{ij} ...` formula. Verdict was reached by interactive corrected sympy + mpmath analysis (Section 3 / 4 above) — the driver's raw output was reframed against the corrected trace formula.
- `debug/data/na1_offdiag_substrate_results.json` — raw panel from the buggy driver, retained for reproducibility and for documenting the bug + diagnosis.
- `debug/sprint_na1_offdiag_substrate_memo.md` — this memo (~3500 words).

### Read context (per sprint prompt)
- `debug/sprint_na1_depth2_mellin_memo.md` — diagonal-substrate Reading C theorem (today, AM).
- `geovac/full_dirac_operator_system.py` — off-diagonal CH substrate (WH1-R3.5).
- `geovac/qed_vertex.py` — Paper 28 vertex-restricted M3 mechanism (NOT used directly; structural argument extends to it).
- `papers/group3_foundations/paper_55_periods_of_geovac.tex` §subsec:m3_diagonal_collapse + §subsec:m3_galois_descent (not modified).
- `papers/group3_foundations/paper_56_tannakian_substrate.tex` (not modified).
- `memory/wh1_proven.md` (R3.5 substrate context).

### References
- Camporesi, R. and Higuchi, A. ``On the eigenfunctions of the Dirac operator on spheres and real hyperbolic spaces.'' J. Geom. Phys. 20 (1996) 1-18.
- Connes, A. and Moscovici, H. ``The local index formula in noncommutative geometry.'' GAFA 5 (1995) 174-243. [JLO cocycle reference for §6.2]
- Paper 28 §4 (Theorem 3, parity discriminant $D_{\mathrm{even}}(s) - D_{\mathrm{odd}}(s) = 2^{s-1}(\beta(s) - \beta(s-2))$).
- Paper 32 §VIII (case-exhaustion theorem at depth 1).
- Paper 18 §III.7 (M1/M2/M3 three-bullet partition).

---

## 10. One-line verdict

**The off-diagonal CH substrate does NOT break the depth-2 → depth-1 collapse**; the trace-of-diagonal-times-matrix-times-diagonal mechanism is more general than simultaneous diagonalisation and applies on any substrate where $\gamma$ enters as a single insertion in a trace of $D$-functions. Original A/B gate verdict: **INCONCLUSIVE in the original gate's terms, READING C-STRONG in the structurally refined sense.** The genuine NA-1 follow-on is **NOT** substrate enrichment (Paper 28 QED vertex graph would behave the same way) but a **non-trace probe** like the JLO depth-2 cocycle on `geovac/jlo_chi.py` (sprint-scale 2–3 weeks).
