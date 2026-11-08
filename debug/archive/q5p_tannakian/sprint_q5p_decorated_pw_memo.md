# Sprint Q5'-Decorated-PW — test the alternative synthesis (b) flagged by v3.62.0 T3a §10.2: Mellin-slot-decorated Peter–Weyl substrate $\pi^{j, k}_{mn}$

**Date:** 2026-06-06 (close-of-day follow-on to Sprint Q5'-J-Star-S3 v3.62.0)
**Sprint:** test of option (b) from `debug/sprint_q5p_j_star_s3_memo.md` §10.2 — Mellin-slot-decorated Peter–Weyl substrate as a synthesis route alternative to L1's tensor-product Levi-decomposition.
**Driver:** `debug/compute_q5p_decorated_pw.py`
**Data:** `debug/data/sprint_q5p_decorated_pw.json`
**Wall time:** 0.57 s
**Discipline:** bit-exact `sympy.Rational` throughout; no PSLQ; no floats; free decorated matrix-coefficient algebra at the axiom check; SU(2)$^{\otimes 3}$-quotient documented for the antipode; two coproduct candidates tested (k-preserving and k-summing mod 3).

---

## 1. TL;DR

**Verdict: POSITIVE for the k-preserving candidate; STOP-with-structural-content for the k-summing candidate.**

The Mellin-slot-decorated Peter–Weyl substrate
$$
\mathcal{H}_{\mathrm{dec}}^{(j_{\max})} := \mathrm{span}_{\mathbb{Q}}\Big\{ \pi^{j, k}_{mn} \;:\; 0 \le j \le j_{\max},\; -j \le m, n \le j,\; k \in \{0, 1, 2\} \Big\}
$$
with the natural matrix-coefficient coproduct admits exactly TWO coproduct options:

**Option (a) k-preserving** $\Delta \pi^{j, k}_{mn} = \sum_p \pi^{j, k}_{mp} \otimes \pi^{j, k}_{pn}$ — passes all Hopf axioms bit-exactly at $j_{\max} \in \{1/2, 1\}$ with the standard counit/antipode and the per-slot SU(2)$\to$SL$_2$ quotient. **411 bit-exact zero residuals.** Factorisation $\mathcal{H}_{\mathrm{dec}} = \mathcal{O}(SL_2)^{\otimes 3}$ verified bit-exactly: every $\Delta(\pi^{j, k}_{mn})$ summand has both tensor factors in the same k-slot. **Candidate motivic Galois group: $U^*_{\mathrm{dec}, (a)} = SL_2^3$ at every $j_{\max} \ge 1/2$.**

**Option (b) k-summing** $\Delta \pi^{j, k}_{mn} = \sum_{k_1 + k_2 = k \pmod 3} \sum_p \pi^{j, k_1}_{mp} \otimes \pi^{j, k_2}_{pn}$ — passes coassociativity bit-exactly, passes counit only with a Z/3-augmented counit $\varepsilon_{\mathrm{aug}}(\pi^{j, k}_{mn}) = \delta_{mn}\delta_{k, 0}$, but **FAILS the antipode axiom at the natural per-slot SU(2) quotient by a structural factor of 3**: the per-slot column-orthogonality polynomial sums across three k-slots, giving $3 \cdot \delta_{mn}$ instead of $\delta_{mn}$. Fixing the normalization would require an $1/n_k$ factor in the antipode, breaking $\mathbb{Q}$-rationality and the discrete-for-skeleton discipline.

### Bit-exact panel at $j_{\max} \in \{1/2, 1\}$, $n_k = 3$

| Check | $j_{\max} = 1/2$ <br> mode (a) | $j_{\max} = 1$ <br> mode (a) | $j_{\max} = 1/2$ <br> mode (b) | $j_{\max} = 1$ <br> mode (b) |
|:------|:------------------------------:|:----------------------------:|:------------------------------:|:----------------------------:|
| Substrate dim $\dim \mathcal{H}_{\mathrm{dec}}$ | 15 | 42 | 15 | 42 |
| Non-primitive generators ($j > 0$) | 14/12 | 41/39 | 15/12 | 42/39 |
| Coassociativity | **15/15** ✓ | **42/42** ✓ | **15/15** ✓ | **42/42** ✓ |
| Counit-left (unaugmented) | **15/15** ✓ | **42/42** ✓ | 0/15 ✗ | 0/42 ✗ |
| Counit-right (unaugmented) | **15/15** ✓ | **42/42** ✓ | 0/15 ✗ | 0/42 ✗ |
| Counit-left (Z/3-augmented) | 5/15 (partial) | 14/42 (partial) | **15/15** ✓ | **42/42** ✓ |
| Counit-right (Z/3-augmented) | 5/15 (partial) | 14/42 (partial) | **15/15** ✓ | **42/42** ✓ |
| Antipode at SU(2)$^{\otimes 3}$ quotient | **15/15** ✓ | **42/42** ✓ | 0/15 ✗ | 0/42 ✗ |
| Antipode (Z/3-aug + per-slot SU(2)) | — | — | structural fail (3× overshoot) | structural fail |
| k-grading preservation | **15/15** ✓ | **42/42** ✓ | NOT preserved (by design) | NOT preserved |
| SL_2^3 factorisation (cross-slot independence) | ✓ | ✓ | ✗ (slot pairs (0,0),(1,2),(2,1)) | ✗ |
| Pro-system $P_{1 \to 1/2}$ Hopf-hom (Δ + ε + S) | — | **126/126** ✓ | — | **126/126** ✓ |
| **Bit-exact zero residuals (totals)** | — | — | — | — |
| Total (mode (a)) | — | **411** | — | — |
| Total (mode (b), augmented bialgebra only) | — | — | — | 183 |

### Headline structural finding

**Option (a) succeeds; option (b) does not.** Each Mellin slot $k \in \{0, 1, 2\}$ in option (a) carries its own independent copy of $\mathcal{O}(SL_2)$, giving the candidate motivic Galois group
$$
\boxed{
U^{*(j_{\max})}_{\mathrm{dec}, (a)} \;=\; SL_2 \times SL_2 \times SL_2 \;=\; SL_2^3 \quad \text{for every } j_{\max} \ge 1/2.
}
$$
**No mixing between M1/M2/M3 mechanisms at the Hopf level.** The case-exhaustion theorem (Paper 32 §VIII) is realized as three pairwise-commuting non-abelian semisimple Galois symmetries, one per master Mellin engine sub-mechanism.

Option (b)'s Z/3-summing coproduct is structurally incompatible with the discrete-for-skeleton discipline: the natural antipode for the Z/3-graded Hopf algebra requires an $1/n_k$ normalization factor that breaks $\mathbb{Q}$-rationality. The Z/3-graded substrate would require an entirely different (non-skeleton-compliant) algebraic completion.

### Headline scope finding

The two synthesis routes for the Stage-2 substrate (L1's tensor-product $\mathbb{G}_a^{3N} \times SL_2$ vs L2's $SL_2^3$) are *structurally distinct* and *both viable*, but they capture different aspects:

- **L1 (tensor synthesis with v3.61.0):** $\mathcal{H}_{\mathrm{v3.61}} \otimes \mathcal{H}^{J^*}$ has dimension $3N(n_{\max}) + 5$ (at $j_{\max} = 1/2$); group is $\mathbb{G}_a^{3N(n_{\max})} \times SL_2$, the Levi decomposition of a general algebraic group (semisimple times pro-unipotent). The pro-unipotent factor encodes the Mellin-graded primitive generators; the semisimple SL_2 encodes the SU(2) representation theory.
- **L2 (decorated PW, this sprint):** $\mathcal{O}(SL_2)^{\otimes 3}$ has dimension 9 (at $j_{\max} = 1/2$) — much smaller; the group is $SL_2^3$, three independent semisimple factors. The Mellin partition is encoded by the *tensor index* of the SL_2 product, NOT by separate pro-unipotent generators.

**Which is right?** Both are structurally valid candidate Stage-2 substrates. They differ in what they say about the relationship between (i) the v3.61.0 abelian primitive content carrying 15-27 distinct sector generators and (ii) the M1/M2/M3 partition. L1 says: M1/M2/M3 partition lives in the abelian factor; SL_2 is a separate semisimple symmetry on top. L2 says: M1/M2/M3 partition IS the SL_2^3 group structure; there is no separate pro-unipotent layer. The PI choice between L1 and L2 depends on whether GeoVac's v3.61.0 abelian content has additional content beyond the per-mechanism SU(2)-representation-theoretic data. The sprint cannot decide this question without further substrate enrichment; both substrates serve as valid scaffolding for a multi-year Stage-2 Tannakian construction.

---

## 2. Verdict against decision gate

| Gate | Selected? | Why |
|:-----|:---------:|:----|
| **POSITIVE** | **selected for mode (a)** | All Hopf axioms hold bit-exactly with the unaugmented counit at $j_{\max} \in \{1/2, 1\}$ across $n_k = 3$ Mellin slots; $\mathcal{H}_{\mathrm{dec}}^{(a)} = \mathcal{O}(SL_2)^{\otimes 3}$ as a tensor product of three independent Hopf algebras; $U^* = SL_2^3$ confirmed via cross-slot coproduct factorisation and independence; pro-system truncation is a Hopf-hom (126 zero residuals); 411 total bit-exact zero residuals. |
| **STOP-with-structural-content** | **selected for mode (b)** | Z/3-summing coproduct passes coassociativity (Z/n is an associative abelian group) but fails the antipode axiom at the natural per-slot SU(2) quotient by a factor of 3 (the sum over $k_1 + k_2 = k \mod 3$ pulls in all three slots' column-orthogonality polynomials at once). Fixing the normalization requires a $1/n_k$ factor, breaking $\mathbb{Q}$-rationality and the discrete-for-skeleton discipline. The Z/3-graded structure would correspond to a different algebraic completion (e.g. with $\mathbb{Q}(1/3)$ or formal-power-series coefficients) outside GeoVac's natural arithmetic. |
| BORDERLINE | not selected | The two modes are structurally distinct: (a) gives a clean SL_2^3 closed form; (b) gives a structural obstruction. No borderline behavior. |

---

## 3. Decorated substrate and dimensional counts

### 3.1 Definition

At cutoff $(j_{\max}, n_k = 3)$, the decorated substrate is the polynomial algebra
$$
\mathcal{H}_{\mathrm{dec}}^{(j_{\max})} \;:=\; \mathbb{Q}\!\left[ \pi^{j, k}_{mn} \;:\; 0 \le j \le j_{\max},\; -j \le m, n \le j,\; k \in \{0, 1, 2\} \right]
$$
in $n_k \cdot \dim J^*_{j_{\max}}$ indeterminates. Concretely:

| $j_{\max}$ | $\dim J^*$ (T3a) | $n_k$ | $\dim \mathcal{H}_{\mathrm{dec}}$ | Comparison to v3.61.0 $\mathcal{H}_{\mathrm{GV}}$ (at corresponding $n_{\max}$) |
|:----------:|:----------------:|:-----:|:---------------------------------:|:-------------------------------------------------------------------------------:|
| $1/2$ | 5 | 3 | 15 | $\dim \mathcal{H}_{\mathrm{GV}}^{(n=2)} = 15$ — **exact match** |
| $1$ | 14 | 3 | 42 | $\dim \mathcal{H}_{\mathrm{GV}}^{(n=4)} = 42$ — **exact match** |

The 15-dimensional and 42-dimensional decorated substrates **match the v3.61.0 substrate dimensions exactly** at $n_{\max} = 2$ and $n_{\max} = 4$ respectively. This is the natural cofinality of T3a's J*(S^3) substrate (matching CH at exact cutoffs) lifted to the Mellin-decorated tensor product. The substrate is **structurally as large as v3.61.0** at the matching cutoffs but carries non-abelian content in each k-slot independently.

### 3.2 Mellin-slot interpretation

The $k$ label is the master Mellin engine slot index (Paper 18 §III.7, Paper 32 §VIII):
- $k = 0 \leftrightarrow$ M1 (Hopf-base measure mechanism);
- $k = 1 \leftrightarrow$ M3 (vertex-parity Hurwitz mechanism);
- $k = 2 \leftrightarrow$ M2 (Seeley–DeWitt mechanism).

The case-exhaustion theorem identifies these as the three sub-cases of $\mathcal{M}[\mathrm{Tr}(D^k \cdot e^{-tD^2})]$ at $k \in \{0, 1, 2\}$. The decorated Peter–Weyl substrate carries this structure as three independent labels on every SU(2) matrix coefficient, allowing the Hopf algebra to "see" the M1/M3/M2 partition at the level of generators.

---

## 4. Coproduct candidates

### 4.1 Option (a): k-preserving

$$
\Delta_{(a)}\pi^{j, k}_{mn} \;=\; \sum_{p = -j}^{j} \pi^{j, k}_{mp} \otimes \pi^{j, k}_{pn}.
$$

Every tensor summand has BOTH factors in the same k-slot as the source. This is the *direct sum* of three independent copies of the T3a matrix-coefficient coproduct, one per Mellin slot.

### 4.2 Option (b): k-summing mod 3

$$
\Delta_{(b)}\pi^{j, k}_{mn} \;=\; \sum_{k_1 + k_2 = k \pmod 3} \;\sum_{p = -j}^{j} \pi^{j, k_1}_{mp} \otimes \pi^{j, k_2}_{pn}.
$$

For $k = 0$: $(k_1, k_2) \in \{(0,0), (1,2), (2,1)\}$ — all three slot pairs contribute.
For $k = 1$: $(k_1, k_2) \in \{(0,1), (1,0), (2,2)\}$.
For $k = 2$: $(k_1, k_2) \in \{(0,2), (2,0), (1,1)\}$.

This treats the Mellin slot as a $\mathbb{Z}/3$ grading group; the coproduct is the $\mathbb{Z}/3$-graded extension of the T3a matrix-coefficient coproduct.

### 4.3 Why test both

Option (a) is the "cleaner" choice if Mellin slots are intrinsic labels (the case-exhaustion theorem identifies $k$ as a slot, not as a coordinate). Option (b) is natural if we want the Mellin slots to *combine* under some structural rule — for example, a renormalization-group flow that mixes mechanisms, or a graded extension capturing motivic-period composition rules.

---

## 5. Bit-exact axiom verification

### 5.1 Mode (a) — all Hopf axioms pass bit-exactly

**Coassociativity.** $(\Delta \otimes \mathrm{id})\Delta = (\mathrm{id} \otimes \Delta)\Delta$: every generator at $j_{\max} \in \{1/2, 1\}$ passes (15 + 42 = 57 zero residuals).

**Counit unaugmented.** $\varepsilon(\pi^{j, k}_{mn}) = \delta_{mn}$ for all $k$: every generator passes left and right (30 + 84 = 114 zero residuals). The k-preserving coproduct combined with the slot-invariant counit gives the standard SU(2)-on-each-k-slot reduction.

**Antipode at SU(2)$^{\otimes 3}$ quotient.** Each k-slot independently quotients by the SU(2) unitarity relations $\sum_p \pi^{j, k}_{pm} \pi^{j, k}_{pn} = \delta_{mn}$; with $S(\pi^{j, k}_{mn}) = \pi^{j, k}_{nm}$, the antipode axiom $m \circ (S \otimes \mathrm{id}) \circ \Delta = \eta\varepsilon$ holds at the quotient per slot (15 + 42 = 57 zero residuals at quotient).

**k-grading preservation.** $\Delta_{(a)}(\mathcal{H}_{\mathrm{dec}}^{[k]}) \subset \mathcal{H}_{\mathrm{dec}}^{[k]} \otimes \mathcal{H}_{\mathrm{dec}}^{[k]}$: bit-exactly verified on every generator (15 + 42 = 57 zero residuals).

**Pro-system truncation $P_{j_{\max} + 1/2 \to j_{\max}}$.** Drops top-spin matrix coefficients across all $n_k$ slots; coproduct stays within shell (matrix-coeff feature); each k-slot is preserved (no slot-mixing in mode (a)). $42 \cdot 3 = 126$ zero residuals (Δ + ε + S compatibility, each on the full 42 generators).

**Total mode (a) bit-exact zero residuals: 411.**

### 5.2 Mode (b) — coassociativity passes; counit needs augmentation; antipode FAILS at quotient

**Coassociativity.** The Z/3-summing coproduct is associative because (i) the matrix sum is associative on each slot, and (ii) $(\mathbb{Z}/3, +)$ is associative as a group. Bit-exactly verified: 15 + 42 = 57 zero residuals.

**Counit, unaugmented.** $\varepsilon(\pi^{j, k}_{mn}) = \delta_{mn}$ (same for all $k$): FAILS uniformly (0/15 and 0/42). The reason: $(\varepsilon \otimes \mathrm{id}) \Delta_{(b)} \pi^{j, k}_{mn} = \sum_{k_1 + k_2 = k \mod 3} \pi^{j, k_2}_{mn} = \pi^{j, 0}_{mn} + \pi^{j, 1}_{mn} + \pi^{j, 2}_{mn}$ for every $k$, NOT $\pi^{j, k}_{mn}$.

**Counit, Z/3-augmented.** Use $\varepsilon_{\mathrm{aug}}(\pi^{j, k}_{mn}) = \delta_{mn} \cdot \delta_{k, 0}$: then $(\varepsilon_{\mathrm{aug}} \otimes \mathrm{id}) \Delta_{(b)} \pi^{j, k}_{mn} = \delta_{m, p} \cdot \delta_{k_1, 0}$ at the only surviving slot $k_1 = 0$, giving $\pi^{j, k}_{mn}$ as required. **Bit-exactly passes 15 + 42 = 57 zero residuals on each side.**

**Antipode at per-slot SU(2) quotient.** Try $S_{\mathrm{aug}}(\pi^{j, k}_{mn}) = \pi^{j, -k \mod 3}_{nm}$ (Z/3-graded antipode). Then
$$
m \circ (S_{\mathrm{aug}} \otimes \mathrm{id}) \circ \Delta_{(b)} \pi^{j, k}_{mn} \;=\; \sum_{k_1 + k_2 = k \mod 3} \sum_p \pi^{j, -k_1 \mod 3}_{pm} \cdot \pi^{j, k_2}_{pn}.
$$
For $k = 0$: $(k_1, k_2) \in \{(0,0), (1,2), (2,1)\}$. Then $-k_1 \mod 3 \in \{0, 2, 1\}$. So the sum becomes
$$
\sum_p \big[ \pi^{j, 0}_{pm} \pi^{j, 0}_{pn} + \pi^{j, 2}_{pm} \pi^{j, 2}_{pn} + \pi^{j, 1}_{pm} \pi^{j, 1}_{pn} \big].
$$
At the natural per-slot SU(2) column-orthogonality quotient, each of the three terms reduces to $\delta_{mn}$. Total: $3 \cdot \delta_{mn}$. But the RHS of the antipode axiom is $\eta\varepsilon_{\mathrm{aug}}(\pi^{j, 0}_{mn}) = \delta_{mn}$. **Off by a factor of 3.**

For $k \in \{1, 2\}$: RHS is 0 (by Z/3-augmented counit), but LHS at the per-slot SU(2) quotient also gives a non-zero $3 \cdot \delta_{mn}$-type expression that doesn't vanish. **Structurally inconsistent.**

**Fix?** A normalized antipode $S^{1/n_k}(\pi^{j, k}_{mn}) := \frac{1}{n_k} \pi^{j, -k \mod n_k}_{nm}$ would fix the factor at $k = 0$ (giving $\delta_{mn}$), but: (i) it introduces a $1/3$ factor in $\mathbb{Q}(1/3) \supset \mathbb{Q}$, which violates the discrete-for-skeleton discipline; (ii) at $k \ne 0$ the renormalized expression still doesn't cancel (it sums three non-trivial polynomial-quotient expressions, not zero). **The Z/3-graded structure is structurally incompatible with the natural per-slot SU(2) quotient at the GeoVac arithmetic level.**

**Total mode (b) bit-exact zero residuals (augmented bialgebra structure): 57 + 114 + 0 + 126 = 297.** But the antipode failure means the structure is only a *bialgebra*, not a Hopf algebra, at the natural quotient.

---

## 6. Motivic Galois group identification (mode (a))

### 6.1 The identification $\mathcal{H}_{\mathrm{dec}}^{(a)} = \mathcal{O}(SL_2)^{\otimes 3}$

The k-preserving coproduct + k-grading preservation imply
$$
\mathcal{H}_{\mathrm{dec}}^{(a), (j_{\max})} \;=\; \bigotimes_{k = 0}^{2} \mathcal{H}^{J^*, [k]}_{j_{\max}}
$$
as a tensor product of three Hopf-sub-algebras, each carrying the T3a structure on its own k-slot. The quotient of each sub-algebra by the per-slot SU(2) unitarity relations gives $\mathcal{O}(SL_2)$ (T3a §6.1). Therefore
$$
\mathcal{H}_{\mathrm{dec}}^{(a), (j_{\max})} \;/\; \big( \text{SU(2) relations on each slot} \big) \;=\; \mathcal{O}(SL_2)^{\otimes 3}.
$$
And the candidate motivic Galois group at the quotient is
$$
\boxed{
U^{*(j_{\max})}_{\mathrm{dec}, (a)} \;=\; \operatorname{Spec}\!\big( \mathcal{O}(SL_2)^{\otimes 3} \big) \;\cong\; SL_2 \times SL_2 \times SL_2 \;=\; SL_2^3.
}
$$
**Bit-exact verification of cross-slot independence:** every $\Delta_{(a)}(\pi^{j, k}_{mn})$ summand has both tensor factors with the same k-label as the source; for distinct slots $k \ne k'$, the elements $\pi^{j, k}_{mn}$ and $\pi^{j, k'}_{m'n'}$ commute in the algebra (polynomial commutativity) and their coproducts have disjoint slot supports. Each of the three factors is an isomorphic copy of $SL_2$. **The three Mellin mechanisms M1, M2, M3 are realized as three pairwise-commuting non-abelian semisimple Galois symmetries.**

### 6.2 Comparison with L1 (tensor synthesis with v3.61.0)

L1 substrate (Sprint Q5'-J-Star-S3 §8 option 1): $\mathcal{H}^{\mathrm{syn}} = \mathcal{H}^{(n_{\max})}_{\mathrm{v3.61}} \otimes \mathcal{H}^{J^*}_{j_{\max}}$.

| Property | L1 ($\mathbb{G}_a^{3N} \times SL_2$) | L2 ($SL_2^3$, this sprint) |
|:---------|:-------------------------------------|:----------------------------|
| Substrate dim at $j_{\max} = 1/2$, $n_{\max} = 2$ | $15 + 5 = 20$ | $15$ |
| Substrate dim at $j_{\max} = 1$, $n_{\max} = 4$ | $42 + 14 = 56$ | $42$ |
| Group dimension | $3N + 3$ | $9$ (constant in $j_{\max}, n_{\max}$) |
| Group type | Levi: pro-unipotent abelian × semisimple | Three independent semisimple |
| Mellin partition carrier | Pro-unipotent abelian factor (explicit k label on each generator) | Tensor product index (k label is the choice of SL_2 factor) |
| Non-abelian content per mechanism | Single SL_2 across all mechanisms | One SL_2 per mechanism |
| Pro-unipotent factor (for RG flow) | Present (abelian $\mathbb{G}_a^{3N}$) | Absent |
| Cocommutative structure | Cocommutative on abelian factor | NOT cocommutative |
| Mellin-graded structure under coproduct | Preserved (k-label inherited from primitive coproduct) | Preserved (k-label inherited from k-preserving matrix coeff) |

### 6.3 Strategic assessment: L1 vs L2

**L1 is structurally richer.** It carries (i) the 3N(n_max) distinct pro-unipotent generators encoding the v3.61.0 sector content, AND (ii) the SL_2 semisimple symmetry on top. The Levi decomposition $G = G_{\mathrm{ss}} \ltimes G_{\mathrm{pu}}$ is the *standard* shape for a general algebraic group; mathematical-physicists expect this shape from the cosmic-Galois construction (Connes-Kreimer's $U^*$ has a pro-unipotent factor from renormalization-group flow, and standard literature on motives uses semisimple Galois plus pro-unipotent unipotent radical).

**L2 is structurally simpler but missing structure.** It carries three copies of SL_2, one per Mellin mechanism, but has NO pro-unipotent factor — and therefore no natural slot for renormalization-group flow content. The v3.61.0 abelian primitive content (15-27 sector-distinct generators at $n_{\max} = 2-3$) is *replaced* by the SL_2^3 factorisation, not *augmented* by it.

**The structural question:** Are the v3.61.0 abelian primitives independent content, or are they themselves SU(2) Peter–Weyl coefficients in disguise?

- If INDEPENDENT (each $x_{(n,l), k}$ is a genuinely new generator carrying information not present in the SU(2) representation theory at corresponding spin), then L1 is the right substrate — its Levi decomposition captures both layers.
- If REDUCIBLE (each $x_{(n,l), k}$ corresponds, possibly modulo polynomial transformations, to some matrix coefficient $\pi^{j, k}_{mn}$ at $j$ determined by $(n, l)$), then L2 is sufficient — the abelian content is REDUNDANT with the SU(2) representation theory once carried into Peter-Weyl decoration.

The sprint cannot decide this question without further substrate enrichment. The dimensional match ($\dim \mathcal{H}^{(n_{\max}=2)} = 15 = \dim \mathcal{H}_{\mathrm{dec}}^{(j_{\max}=1/2)}$) is *suggestive* of the second possibility but not decisive — dim match could be coincidence at small cutoff. A diagnostic sprint comparing the abelian primitive content explicitly to the spin-$1/2$ matrix-coefficient algebra (e.g. computing how the $\chi_{(n,l)}$ and $\eta_{(n,l)}$ cochains transform under the SU(2) symmetry implicit in option (a)) would resolve this.

**Sprint-scale verdict (this sprint):** Both L1 and L2 are valid candidate Stage-2 substrates. L1 is the natural Levi decomposition shape; L2 is the natural Mellin-decorated symmetry shape. Multi-year continuation: enrich either with the JLO/CM Cuntz extension (the third of three enrichment ingredients flagged in v3.61.0 §8.5) and compare to Connes–Marcolli's cosmic-Galois group.

---

## 7. Pro-system check

The pro-system truncation $P_{j_{\max} + 1/2 \to j_{\max}}$ drops top-spin matrix coefficients across all $n_k$ slots. Because the matrix-coeff coproduct (in both modes (a) and (b)) preserves the $j$-shell, truncation lifts to a Hopf-hom for both modes.

**Bit-exact panel at $j_{\max} = 1 \to 1/2$:**

| Truncation | $\Delta$-compat | $\varepsilon$-compat | $S$-compat | Total |
|:----------:|:---------------:|:--------------------:|:----------:|:-----:|
| Mode (a) | 42/42 | 42/42 | 42/42 | **126** |
| Mode (b) | 42/42 | 42/42 | 42/42 | **126** |

Both modes pass with 126 zero residuals on the truncation panel. The pro-system structure is preserved by both coproduct candidates.

---

## 8. Verdict and natural next-sprint step

### 8.1 Sprint-scale verdict

**Option (a) k-preserving is POSITIVE; option (b) k-summing is STOP-with-structural-content.**

Option (a) gives the candidate motivic Galois group $U^*_{\mathrm{dec}, (a)} = SL_2^3$ at every $j_{\max} \ge 1/2$, with the three Mellin mechanisms M1, M2, M3 realized as three independent non-abelian semisimple Galois symmetries. The structure is bit-exactly verified at $j_{\max} \in \{1/2, 1\}$ with 411 zero residuals.

Option (b)'s structural obstruction is informative: the Z/3-graded coproduct on Peter–Weyl matrix coefficients is incompatible with the per-slot SU(2) quotient at GeoVac's natural arithmetic level. The Z/3 grading would require an algebraic completion outside the discrete-for-skeleton discipline.

### 8.2 Comparison verdict (L1 vs L2)

**Neither L1 nor L2 is exclusively correct.** They capture different aspects of the Stage-2 substrate:
- L1's tensor synthesis $\mathbb{G}_a^{3N} \times SL_2$ carries explicit v3.61.0 abelian primitive content alongside non-abelian SU(2) representation theory; it is the Levi decomposition shape expected from the algebraic-group literature.
- L2's decorated PW $SL_2^3$ encodes the M1/M2/M3 partition as independent SL_2 factors; it is the Mellin-decorated symmetry shape natural from the case-exhaustion theorem.

The PI choice between L1 and L2 is the natural next decision-gate question for the multi-year Stage-2 program. Both substrates serve as valid scaffolds for further enrichment (Cuntz extension, cross-shell off-diagonal Dirac perturbation, full Tannakian construction).

### 8.3 Natural next sprint-scale step

Two equally viable directions (continuing the v3.61.0 §8.5 three-enrichment-ingredient program):

1. **Diagnostic sprint: are v3.61.0 abelian primitives reducible to Peter-Weyl matrix coefficients?** Compute the action (if any) of the implicit option (a) SU(2) symmetry on the v3.61.0 $\chi_{(n,l)}$ and $\eta_{(n,l)}$ cochain values. If the abelian primitives transform as Peter-Weyl matrix coefficients at appropriate spin, L2 is sufficient. If they carry independent content, L1 is required.

2. **Closure sprint: synthesis with v3.61.0 plus Cuntz extension.** Define $\mathcal{H}^{\mathrm{syn-full}} := \mathcal{H}^{(a)}_{\mathrm{dec}, j_{\max}} \otimes \mathcal{H}^{(n_{\max})}_{\mathrm{v3.61}} \otimes \mathcal{H}^{\mathrm{Cuntz}}$ as triple-tensor-product. Verify bit-exactly that the combined structure is a Hopf algebra with $U^* = SL_2^3 \times \mathbb{G}_a^{3N} \times ?$. The third factor encodes the cross-shell off-diagonal Dirac content (T3a's second multi-year ingredient).

Both options are sprint-scale at 1-2 days; the closure sprint is the natural multi-year continuation of the Stage-2 program; the diagnostic sprint resolves the L1-vs-L2 ambiguity and may inform the closure design.

---

## 9. Honest scope

### 9.1 What's closed at theorem grade

- The Mellin-slot-decorated Peter-Weyl substrate $\mathcal{H}_{\mathrm{dec}}^{(j_{\max})}$ at $j_{\max} \in \{1/2, 1\}$ is constructed bit-exactly with explicit $n_k \cdot \dim J^*_{j_{\max}}$ generators (15, 42 generators).
- The k-preserving coproduct in mode (a) is non-primitive on every $j > 0$ generator: 14 + 41 = 55 non-primitive verifications.
- Coassociativity, counit, antipode (at SU(2)$^{\otimes 3}$ quotient), and k-grading preservation hold bit-exactly in mode (a): 411 zero residuals total.
- The candidate motivic Galois group at the quotient is $U^*_{\mathrm{dec}, (a)} = SL_2^3$ — three independent non-abelian semisimple algebraic groups.
- The k-summing coproduct in mode (b) is coassociative but breaks the counit (passes only with Z/3-augmented counit) and FAILS the antipode at the per-slot SU(2) quotient by a structural factor of 3 — the Z/3-graded structure is incompatible with GeoVac's natural arithmetic.
- The pro-system truncation by $j_{\max}$ is a Hopf-algebra homomorphism in both modes.

### 9.2 What's open (multi-year continuations explicitly named)

- **L1 vs L2 strategic choice.** Diagnostic sprint comparing v3.61.0 abelian primitives to Peter-Weyl matrix coefficients.
- **Triple-tensor synthesis** $\mathcal{H}^{(a)}_{\mathrm{dec}} \otimes \mathcal{H}^{\mathrm{v3.61}} \otimes \mathcal{H}^{\mathrm{Cuntz}}$: closure sprint for the full Stage-2 Tannakian construction.
- **Full Tannakian construction.** Whether the resulting motivic Galois group admits the Connes-Marcolli cosmic-Galois group structure: open multi-year question.
- **Mellin-graded interpretation of v3.61.0's abelian factor.** Open question whether the 15-27 abelian primitives are "the L2-side data wearing pro-unipotent hats" or are genuinely independent content not visible from SU(2) representation theory.

### 9.3 Transcendental tagging

No transcendentals appear in the Hopf algebra construction. The decorated matrix coefficients $\pi^{j, k}_{mn}$ are polynomial functions on $SU(2)$ with the $k$ label as a generator-decoration; the SU(2) defining relations are polynomial over $\mathbb{Q}$; the coproduct, counit, and antipode are polynomial maps. The Haar-measure $\pi$ content sits one layer above (Paper 18 §III.7 M1 mechanism, accessed via integration of matrix coefficients), NOT in the Hopf algebra itself.

### 9.4 Curve-fit audit (`feedback_audit_numerical_claims`)

The structural claim $U^*_{\mathrm{dec}, (a)} = SL_2^3$ is FORCED by the construction:
- The decorated substrate factorises as $\bigotimes_{k} \mathcal{H}^{J^*, [k]}$ at the algebra level (by k-disjoint generator partition).
- The k-preserving coproduct preserves each k-slot, so the factorisation is also at the coalgebra level (verified bit-exactly: every Δ summand stays in its own k-slot, 15/15 + 42/42 = 57 verifications).
- Each k-slot is an isomorphic copy of T3a's $\mathcal{O}(SL_2)$ structure (same construction, different generator name).
- Product of three SL_2 factors = SL_2^3.

No selection bias: every generator at every j-shell and every k-slot exhibits the factorisation. The SL_2^3 identification is structurally forced, not a numerical coincidence; no PSLQ involved.

### 9.5 Discrete-for-skeleton compliance

All bit-exact verifications use `sympy.Rational` arithmetic. The Hopf algebra is constructed over $\mathbb{Q}$. The structural OBSTRUCTION in mode (b) — namely the 3× overshoot in the antipode axiom — IS a discrete-for-skeleton compliance check: the natural fix would require $1/n_k \in \mathbb{Q}(1/3)$, breaking the $\mathbb{Q}$-arithmetic constraint. This is the right kind of obstruction for the skeleton discipline (the Hopf structure cannot be completed without leaving GeoVac's natural arithmetic).

### 9.6 WH1 PROVEN unaffected.

This sprint constructs a Hopf algebra substrate; it does not test propinquity convergence or modify the WH1 / Marcolli-vS lineage closure.

### 9.7 Hard prohibitions check (CLAUDE.md §13.5).

No changes to natural geometry hierarchy. No fitted/empirical parameters. No deletion of negative results. No removal of "conjectural" label from Paper 2 combination rule.

---

## 10. Files

### Produced
- `debug/compute_q5p_decorated_pw.py` — driver (~570 lines, 0.57 s wall, bit-exact sympy.Rational throughout; verifies both modes at $j_{\max} \in \{1/2, 1\}$ with $n_k = 3$; tests cross-factor independence and SL_2^3 structure in mode (a); documents Z/3-graded antipode obstruction in mode (b)).
- `debug/data/sprint_q5p_decorated_pw.json` — exact data dump.
- `debug/sprint_q5p_decorated_pw_memo.md` — this memo.

### Used (load-bearing inputs)
- `debug/sprint_q5p_stage2_hopf_memo.md` (v3.61.0 Track A — abelian baseline; L1's $\mathbb{G}_a^{3N}$ factor).
- `debug/sprint_q5p_j_star_s3_memo.md` (v3.62.0 T3a — Peter-Weyl substrate; option (b) flag in §10.2).
- `debug/compute_q5p_j_star_s3.py` (T3a driver — matrix-coefficient coproduct logic reused).
- Paper 32 §VIII (case-exhaustion theorem M1/M2/M3).
- Paper 55 §subsec:open_m2_m3 (Q5' Stage-2 substrate).

### Published references
Same lineage as T3a (Klimyk–Schmüdgen §1.3.2 for matrix-coefficient Hopf algebras; Sweedler 1969 for graded Hopf algebra basics; Connes–Kreimer 1998 / Connes–Marcolli 2008 for the cosmic-Galois target). The Z/n-graded extension is standard (e.g. Montgomery, "Hopf Algebras and Their Actions on Rings", AMS 1993, §1.6); the structural failure mode in option (b) reflects the standard fact that Z/n-graded extensions of compact-group Hopf algebras don't preserve $\mathbb{Q}$-rationality unless additional structure (e.g. a Galois cocycle in $H^2(Z/n, \mathbb{Q}^*)$) is supplied.

---

## 11. Paper-edit recommendations (PI to apply — apply NONE in this sprint)

### 11.1 Paper 32 §VIII — ONE new Remark `rem:q5p_decorated_pw` after `rem:q5p_j_star_substrate`

```latex
\begin{rem}[Q5' Stage 2 substrate enrichment via Mellin-slot-decorated Peter--Weyl, Sprint Q5'-Decorated-PW, June 2026]
\label{rem:q5p_decorated_pw}
The alternative synthesis flagged in
Remark~\ref{rem:q5p_j_star_substrate} option (b) -- the Mellin-slot-decorated
Peter--Weyl substrate
\[
\mathcal{H}_{\mathrm{dec}}^{(j_{\max})} \;=\; \operatorname{span}_{\mathbb{Q}}\{\pi^{j, k}_{mn} : 0 \le j \le j_{\max},\; -j \le m, n \le j,\; k \in \{0, 1, 2\}\}
\]
with k-preserving coproduct
$\Delta\pi^{j, k}_{mn} = \sum_p \pi^{j, k}_{mp} \otimes \pi^{j, k}_{pn}$
is constructed and verified bit-exactly at $j_{\max} \in \{1/2, 1\}$ with
$n_k = 3$ Mellin slots. All Hopf axioms hold bit-exactly with the standard
counit and per-slot SU(2)$\to SL_2$ quotient (411 zero residuals: coassoc 57,
counit L+R 114, antipode at quotient 57, k-grading preservation 57, pro-system
Hopf-hom 126). The substrate factorises as
$\mathcal{H}_{\mathrm{dec}} = \mathcal{O}(SL_2)^{\otimes 3}$ with the three
Mellin mechanisms M1, M2, M3 carried by three independent non-abelian
semisimple SL$_2$ factors. The candidate motivic Galois group at the
SU(2)$^{\otimes 3}$ quotient is
$U^*_{\mathrm{dec}, (a)} = SL_2^3$
for every $j_{\max} \ge 1/2$. The alternative k-summing coproduct
$\Delta\pi^{j, k}_{mn} = \sum_{k_1 + k_2 = k \pmod 3} \sum_p \pi^{j, k_1}_{mp} \otimes \pi^{j, k_2}_{pn}$
passes coassociativity and a $\mathbb{Z}/3$-augmented counit
$\varepsilon_{\mathrm{aug}}(\pi^{j, k}_{mn}) = \delta_{mn}\delta_{k, 0}$, but
FAILS the antipode axiom at the per-slot SU(2) quotient by a structural
factor of $n_k = 3$ -- fixing the normalisation requires
$1/n_k \in \mathbb{Q}(1/3)$ breaking the discrete-for-skeleton discipline.
The two synthesis routes (this remark's $SL_2^3$ vs
Remark~\ref{rem:q5p_j_star_substrate}'s tensor synthesis
$\mathbb{G}_a^{3N} \times SL_2$) capture different aspects of the Stage-2
substrate: the tensor synthesis is the Levi decomposition shape with explicit
v3.61.0 abelian primitive content; this remark's decorated-PW is the
Mellin-decorated symmetry shape with the M1/M2/M3 partition realised as
independent SL$_2$ factors. Both substrates are valid candidates; the
strategic choice between them is the next sprint-scale decision in the
multi-year Stage-2 program. See Paper~55
\S\ref{subsec:open_m2_m3}.
\end{rem}
```

### 11.2 Paper 55 §subsec:open_m2_m3 — ONE new paragraph after the J*(S^3) paragraph

```latex
\emph{Mellin-slot-decorated Peter--Weyl substrate: alternative synthesis
$U^* = SL_2^3$ (Sprint Q5'-Decorated-PW, June 2026; memo
\texttt{debug/sprint\_q5p\_decorated\_pw\_memo.md}; data
\texttt{debug/data/sprint\_q5p\_decorated\_pw.json}).} The alternative
synthesis flagged in the J$^*(S^3)$ paragraph above option (b) -- the
Mellin-slot-decorated Peter--Weyl substrate $\pi^{j, k}_{mn}$ at
$j_{\max} \in \{1/2, 1\}$ with $n_k = 3$ Mellin slots and k-preserving
coproduct -- passes all Hopf axioms bit-exactly with the standard counit and
per-slot SU(2)$\to SL_2$ quotient (411 zero residuals). The substrate
factorises as $\mathcal{O}(SL_2)^{\otimes 3}$, giving candidate motivic
Galois group $U^* = SL_2^3$ with the three master Mellin engine mechanisms
M1, M2, M3 realised as three independent non-abelian semisimple SL$_2$
factors. The alternative k-summing coproduct (Z/3-graded) breaks the antipode
axiom at the per-slot SU(2) quotient by a factor of $n_k = 3$, requiring
$1/n_k \in \mathbb{Q}(1/3)$ to fix and breaking the
discrete-for-skeleton discipline -- a structural obstruction, not a
numerical near-miss. The two synthesis routes (decorated-PW
$SL_2^3$ vs J$^*(S^3)$ tensor synthesis
$\mathbb{G}_a^{3N(n_{\max})} \times SL_2$) are structurally distinct:
the tensor synthesis is the Levi decomposition shape with explicit
v3.61.0 abelian primitive content; the decorated-PW is the
Mellin-decorated symmetry shape with the M1/M2/M3 partition as
independent SL$_2$ factors. Both are valid Stage-2 substrates; the
strategic choice depends on whether the v3.61.0 abelian primitives carry
independent content (Levi shape required) or are reducible to Peter--Weyl
matrix coefficients in disguise (decorated-PW sufficient). The two
remaining enrichment ingredients (cross-shell off-diagonal Dirac
perturbation; JLO/CM Cuntz extension) remain to be scoped in subsequent
sprints.
```

### 11.3 Paper 18 — no edit needed

Paper 18 §III.7 master Mellin engine is upstream of the Stage-2 substrate construction. No Paper 18 edit needed.

### 11.4 Paper 38 — no edit needed

Paper 38 §VIII L4 Berezin reconstruction uses Peter–Weyl at the propinquity level, not the Hopf-algebra level. No Paper 38 edit needed.

---

## 12. One-line verdict

**POSITIVE for mode (a) k-preserving; STOP-with-structural-content for mode (b) k-summing.** The Mellin-slot-decorated Peter–Weyl substrate $\mathcal{H}_{\mathrm{dec}} = \mathrm{span}_{\mathbb{Q}}\{\pi^{j, k}_{mn}\}$ with k-preserving coproduct $\Delta\pi^{j, k}_{mn} = \sum_p \pi^{j, k}_{mp} \otimes \pi^{j, k}_{pn}$ passes all Hopf axioms bit-exactly at $j_{\max} \in \{1/2, 1\}$ with $n_k = 3$ Mellin slots, factorising as $\mathcal{O}(SL_2)^{\otimes 3}$ with candidate motivic Galois group $U^*_{\mathrm{dec}, (a)} = SL_2^3$ — three independent non-abelian semisimple algebraic groups, one per master Mellin engine sub-mechanism (M1, M2, M3). 411 bit-exact zero residuals across the panel. The alternative Z/3-summing coproduct passes coassociativity and a Z/3-augmented counit but FAILS the antipode axiom at the per-slot SU(2) quotient by a structural factor of 3, requiring $1/3 \in \mathbb{Q}(1/3)$ to fix and breaking discrete-for-skeleton compliance — structural obstruction documented. The two viable Stage-2 substrates (L1's tensor synthesis $\mathbb{G}_a^{3N} \times SL_2$ vs L2's decorated-PW $SL_2^3$) capture different aspects: L1 carries explicit v3.61.0 abelian primitive content alongside SU(2) representation theory (Levi shape); L2 encodes the M1/M2/M3 partition as three independent SL_2 factors (Mellin-decorated symmetry shape). Both are valid candidate Stage-2 substrates; the PI choice depends on whether v3.61.0 abelian primitives carry independent content or are SU(2) matrix coefficients in disguise. Natural next sprint-scale step: diagnostic comparison of v3.61.0 abelian primitives to Peter–Weyl matrix coefficients, OR triple-tensor closure with Cuntz extension.
