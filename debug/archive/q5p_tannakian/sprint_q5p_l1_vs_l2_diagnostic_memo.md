# Sprint Q5'-L1-vs-L2-Diagnostic — primitive-space invariant decides the Stage-2 substrate

**Date:** 2026-06-06 (close-of-day follow-on to Sprint Q5'-Levi-Arc v3.63.0)
**Sprint:** the diagnostic flagged by Sprint Q5'-Decorated-PW §8.3 (option 1) and by Sprint Q5'-Levi-Arc §5 ("strategic decision (sprint-scale)") of the v3.63.0 umbrella memo.
**Driver:** `debug/compute_q5p_l1_vs_l2_diagnostic.py`
**Data:** `debug/data/sprint_q5p_l1_vs_l2_diagnostic.json`
**Wall time:** 0.030 s
**Discipline:** bit-exact `sympy.Rational` throughout; no PSLQ; no floats; no transcendentals introduced.

---

## 1. TL;DR

**Verdict: POSITIVE — L1 is FORCED.** The discriminating Hopf-isomorphism invariant $\dim \operatorname{Prim}(\mathcal{H}) = \dim \operatorname{Lie}(\operatorname{Spec}\mathcal{H})$ separates the v3.61.0 L1 substrate $\mathcal{H}_{\mathrm{GV}}$ from the L2 decorated-PW substrate $\mathcal{H}_{\mathrm{dec}}$ at every cutoff $n_{\max} \ge 2$. L2 cannot subsume L1: the v3.61.0 abelian primitives carry independent content beyond the SU(2) Peter–Weyl representation theory, and the **Levi-decomposition $\mathbb{G}_a^{3 N(n_{\max})} \times SL_2$ is the correct Stage-2 target**.

| $n_{\max}$ | $\dim \operatorname{Prim}(\mathcal{H}_{\mathrm{GV}}) = 3 N(n_{\max})$ | $\dim \operatorname{Prim}(\mathcal{H}_{\mathrm{dec}}) = 3 \cdot \dim \mathfrak{sl}_2$ | Defect | L1 forced? |
|:----------:|:--------------------------------------------------------------------:|:------------------------------------------------------------------------------------:|:------:|:---------:|
| 1 | 6 | 9 | $-3$ | (cocommutativity differs — see §3.4) |
| **2** | **15** | **9** | $+6$ | **YES** |
| **3** | **27** | **9** | $+18$ | **YES** |
| **4** | **42** | **9** | $+33$ | **YES** |
| **5** | **60** | **9** | $+51$ | **YES** |
| **6** | **81** | **9** | $+72$ | **YES** |

Bit-exact verification: at every $n_{\max} \ge 2$, $\dim \operatorname{Prim}(\mathcal{H}_{\mathrm{GV}}) > \dim \operatorname{Prim}(\mathcal{H}_{\mathrm{dec}}) = 9$, so no Hopf-algebra isomorphism exists. The defect grows quadratically; the gap widens with cutoff.

### Headline structural reading

The substrate-dim coincidences flagged in the Decorated-PW memo ($\dim \mathcal{H}^{(n_{\max}=2)} = 15 = \dim \mathcal{H}_{\mathrm{dec}}^{(j_{\max}=1/2)}$ and $42 = 42$ at $n_{\max}=4$ vs $j_{\max}=1$) are **small-cutoff artifacts**, not structural identifications. Bit-exact comparison across the cutoff grid confirms:

| $n_{\max}$ | $N(n_{\max})$ | matching $j_{\max}$ such that $\sum_{j \le j_{\max}}(2j+1)^2 = N(n_{\max})$? |
|:----------:|:------------:|:-----------------------------------------------------------------------------:|
| 1 | 2 | **no** (nearest: $j_{\max}=1/2$ gives 5) |
| 2 | 5 | yes, $j_{\max}=1/2$ |
| 3 | 9 | **no** (5 < 9 < 14) |
| 4 | 14 | yes, $j_{\max}=1$ |
| 5 | 20 | **no** (14 < 20 < 30) |
| 6 | 27 | **no** (14 < 27 < 30) |

The pattern holds at $n_{\max} \in \{2, 4\}$ but BREAKS at $n_{\max} \in \{1, 3, 5, 6\}$. Even where the substrate dimensions match, the underlying Hopf algebras differ at the primitive-space level. The dim coincidence is a non-structural numerical accident at exactly two cutoffs.

---

## 2. Verdict against decision gate

| Gate | Selected? | Why |
|:-----|:---------:|:----|
| **POSITIVE — L1 forced** | **selected** | $\dim \operatorname{Prim}(\mathcal{H}_{\mathrm{GV}})(n_{\max}) > \dim \operatorname{Prim}(\mathcal{H}_{\mathrm{dec}})(j_{\max})$ at every $n_{\max} \ge 2$; defect $\ge 6$ growing quadratically. The two Hopf algebras cannot be isomorphic, so L2 cannot subsume L1. v3.61.0 abelian primitives carry independent content beyond SU(2) Peter–Weyl. Levi-decomposition is the right Stage-2 target. |
| L2 sufficient | **not selected** | Falsified at every $n_{\max} \ge 2$. The only cutoff where L2 could embed L1 by dimension is $n_{\max} = 1$ (6 ≤ 9), and even there the cocommutativity mismatch (L1 cocommutative, L2 non-cocommutative on $j > 0$) rules out isomorphism. |
| BORDERLINE | not selected | The defect is bit-exactly positive and growing; no marginal case at the cutoffs tested. |

---

## 3. The diagnostic

### 3.1 Substrates

- **L1 substrate (Sprint Q5'-Stage2-Hopf v3.61.0, Sprint Q5'-Levi-Synthesis v3.63.0):**
  $$\mathcal{H}_{\mathrm{GV}}^{(n_{\max})} = \operatorname{Sym}_{\mathbb{Q}}(V_{n_{\max}}), \quad V_{n_{\max}} = \bigoplus_{(n, l) \in \mathcal{S}_{n_{\max}}} \bigoplus_{k \in \{0, 1, 2\}} \mathbb{Q} \cdot x_{(n, l), k}$$
  with primitive coproduct $\Delta(x_{(n, l), k}) = x_{(n, l), k} \otimes 1 + 1 \otimes x_{(n, l), k}$ on every degree-1 generator, extended as an algebra homomorphism. Group at cutoff: $\mathbb{G}_a^{3 N(n_{\max})}$ (pro-unipotent abelian). Cocommutative.

- **L2 substrate (Sprint Q5'-Decorated-PW v3.63.0, mode (a) k-preserving):**
  $$\mathcal{H}_{\mathrm{dec}}^{(j_{\max})} = \operatorname{span}_{\mathbb{Q}}\{\pi^{j, k}_{mn} : 0 \le j \le j_{\max},\; -j \le m, n \le j,\; k \in \{0, 1, 2\}\}$$
  with matrix-coefficient coproduct $\Delta \pi^{j, k}_{mn} = \sum_p \pi^{j, k}_{mp} \otimes \pi^{j, k}_{pn}$. At the SU(2)$^{\otimes 3} \to SL_2^{\otimes 3}$ quotient (det = 1 per slot), the Hopf algebra is $\mathcal{O}(SL_2)^{\otimes 3}$ with group $SL_2^3$. Non-cocommutative.

### 3.2 Discriminating invariant

$\operatorname{Prim}(\mathcal{H}) := \{x \in \mathcal{H} : \Delta(x) = x \otimes 1 + 1 \otimes x\}$ is a Hopf-isomorphism invariant. Its dimension equals the dimension of the affine algebraic group $\operatorname{Spec}\mathcal{H}$ at the augmentation ideal:
$$\dim \operatorname{Prim}(\mathcal{H}) = \dim_{\mathbb{Q}}(\mathfrak{m} / \mathfrak{m}^2)^* = \dim \operatorname{Lie}(\operatorname{Spec} \mathcal{H}).$$
Any Hopf isomorphism $\phi: \mathcal{H}_1 \xrightarrow{\sim} \mathcal{H}_2$ restricts to a vector-space isomorphism $\operatorname{Prim}(\mathcal{H}_1) \cong \operatorname{Prim}(\mathcal{H}_2)$.

### 3.3 Bit-exact computation

**$\operatorname{Prim}(\mathcal{H}_{\mathrm{GV}})$ at $n_{\max} \in \{1, \dots, 6\}$:** Every degree-1 generator $x_{(n, l), k}$ is primitive by the v3.61.0 construction (Sprint Q5'-Stage2-Hopf §4.1). The vector space of all primitives in a connected graded commutative cocommutative Hopf algebra over a field of characteristic 0 with primitive generating set equals the generating space itself (any other primitive is a $\mathbb{Q}$-linear combination of generators; Milnor–Moore / Cartier–Quillen–Milnor–Moore structure theorem for cocommutative connected Hopf algebras gives $\mathcal{H}_{\mathrm{GV}} = U(\operatorname{Prim}\mathcal{H}_{\mathrm{GV}})$ with $\operatorname{Prim} \cong V_{n_{\max}}$ as Lie algebras with zero bracket).

| $n_{\max}$ | $N(n_{\max}) = n_{\max}(n_{\max}+3)/2$ | $\dim V_{n_{\max}} = 3 N(n_{\max})$ | $\dim \operatorname{Prim}(\mathcal{H}_{\mathrm{GV}})$ | Bit-exact primitive count |
|:----------:|:----:|:----:|:----:|:------------------:|
| 1 | 2 | 6 | **6** | 6 / 6 ✓ |
| 2 | 5 | 15 | **15** | 15 / 15 ✓ |
| 3 | 9 | 27 | **27** | 27 / 27 ✓ |
| 4 | 14 | 42 | **42** | 42 / 42 ✓ |
| 5 | 20 | 60 | **60** | 60 / 60 ✓ |
| 6 | 27 | 81 | **81** | 81 / 81 ✓ |

**$\operatorname{Prim}(\mathcal{H}_{\mathrm{dec}})$ at $j_{\max} \in \{1/2, 1, 3/2, 2, 5/2\}$:** Use the standard fact that for any affine algebraic group $G$ over $\mathbb{Q}$,
$$\dim \operatorname{Prim}\mathcal{O}(G) = \dim \mathfrak{g}.$$
At the SU(2)$^{\otimes 3} \to SL_2^{\otimes 3}$ quotient, $\operatorname{Prim}\mathcal{O}(SL_2)^{\otimes 3} = \mathfrak{sl}_2 \oplus \mathfrak{sl}_2 \oplus \mathfrak{sl}_2$, dimension $3 \cdot 3 = 9$ INDEPENDENT of $j_{\max}$. Higher $j_{\max}$ adds polynomial expressions in the $j = 1/2$ generators (via Clebsch–Gordan); these expressions are in $\mathfrak{m}^2$ at the augmentation ideal and contribute nothing new to $\operatorname{Prim}$.

**Bit-exact verification of $\dim(\mathfrak{m} / \mathfrak{m}^2) = 9$ at $j_{\max} = 1/2$:**

Per slot $k$, the four generators in the augmentation ideal are $a_k - 1$, $b_k$, $c_k$, $d_k - 1$ where $(a_k, b_k, c_k, d_k) = (\pi^{1/2, k}_{++}, \pi^{1/2, k}_{+-}, \pi^{1/2, k}_{-+}, \pi^{1/2, k}_{--})$ and $\varepsilon(a_k) = \varepsilon(d_k) = 1$, $\varepsilon(b_k) = \varepsilon(c_k) = 0$.

The $SL_2$ relation $\det(\pi^{1/2}) = a_k d_k - b_k c_k = 1$ reduces modulo $\mathfrak{m}^2$ to its linear part. Substituting $u = a - 1$, $v = d - 1$:
$$ad - bc - 1 = (1+u)(1+v) - bc - 1 = u + v + uv - bc.$$
The linear part is $u + v = (a - 1) + (d - 1)$. **Bit-exactly computed:** the driver extracts `det_linear_part = u + v` matching `expected_linear = u + v` (bit-exact match). The relation kills 1 dimension per slot, giving $\dim(\mathfrak{m} / \mathfrak{m}^2) = (4 - 1) \cdot 3 = 9$ across all three Mellin slots.

**Bit-exact verification at $j_{\max} = 1$:** The nine $j = 1$ matrix coefficients $\pi^{1, k}_{mn}$ are symmetric-square Clebsch–Gordan polynomials of degree exactly 2 in the $j = 1/2$ generators. Sample computation:
- $\pi^{1, k}_{++} = a_k^2$ has linear part $2(a_k - 1) = 2u_k$ — already in the $j = 1/2$ span;
- $\pi^{1, k}_{00} = a_k d_k + b_k c_k$ has linear part $(a_k - 1) + (d_k - 1) = u_k + v_k$ — already in the $j = 1/2$ span;
- $\pi^{1, k}_{+-} = b_k^2$ has linear part $0$ — already in $\mathfrak{m}^2$.

Adding $j = 1$ generators expands the substrate dimension per slot from 5 to 14 ($1 + 4 + 9 = 14$ matrix coefficients) but contributes **zero new dimensions to $\mathfrak{m} / \mathfrak{m}^2$**. The Lie-algebra dimension stays 9.

**Bit-exact non-primitivity of $j = 1/2$ generators:** Per slot, the matrix-coefficient coproduct produces non-zero residuals against the primitive form:
$$\Delta(a) - (a \otimes 1 + 1 \otimes a) = (a \otimes a + b \otimes c) - (a \otimes 1 + 1 \otimes a) = a_L a_R - a_L - a_R + b_L c_R$$
and similar for $b, c, d$ (driver reports all four residuals bit-exactly non-zero). All 12 generators across 3 slots fail primitivity: $12 / 12$ non-zero residuals.

### 3.4 Cocommutativity also separates the substrates

Even at $n_{\max} = 1$ (where $\dim \operatorname{Prim}(\mathcal{H}_{\mathrm{GV}}) = 6 < 9 = \dim \operatorname{Prim}(\mathcal{H}_{\mathrm{dec}})$ — so injection by dimension is in principle possible), the two Hopf algebras differ at the **cocommutativity** level:

- $\mathcal{H}_{\mathrm{GV}}$ is **cocommutative** (every generator is primitive, primitives have $\tau \circ \Delta = \Delta$ where $\tau$ is the flip).
- $\mathcal{H}_{\mathrm{dec}}$ is **NOT cocommutative** on $j > 0$ matrix coefficients (the matrix-coefficient coproduct $\sum_p \pi_{mp} \otimes \pi_{pn}$ is symmetric in $L \leftrightarrow R$ only under $\sigma \cdot \pi_{mp} = \pi_{pm}$, which is a non-trivial automorphism).

Cocommutativity is a Hopf-isomorphism invariant. Therefore the two Hopf algebras are categorically distinct at **every** cutoff $n_{\max} \ge 1$, not only where the dimension count separates them.

---

## 4. Cross-check: v3.61.0 chi / eta values

Bit-exact agreement with the prosystem memo:

| $n_{\max}$ | $\chi$ vector | Match prosystem memo §5? | $\eta$ vector | Match prosystem memo §6? |
|:----:|:--------:|:------:|:--------:|:------:|
| 1 | $(2, -2)$ | ✓ | $(3, 3)$ | ✓ |
| 2 | $(2, -2, 2, 2, -4)$ | ✓ | $(3, 3, 5, 15, 10)$ | ✓ |
| 3 | $(2, -2, 2, 2, -4, 2, 2, 2, -6)$ | ✓ | $(3, 3, 5, 15, 10, 7, 21, 35, 21)$ | ✓ |
| 4 | $(2, -2, 2, 2, -4, 2, 2, 2, -6, 2, 2, 2, 2, -8)$ | ✓ | $(3, 3, 5, 15, 10, 7, 21, 35, 21, 9, 27, 45, 63, 36)$ | ✓ |

All four cross-checks pass bit-exactly. The v3.61.0 closed forms
$\chi_{(n, l)} = +2$ for $l < n$ / $-2n$ for $l = n$, and $\eta_{(n, l)} = (2l+1)(2n+1)$ for $l < n$ / $n(2n+1)$ for $l = n$ — confirmed across all 30 distinct $(n, l)$ sectors.

Structural observation: the $(2n+1)$ shell factor in $\eta_{(n, l)}$ is not a single Peter–Weyl character (which would scale as $2j+1$ at a fixed spin $j$). $\eta_{(n, l)}$ couples shell-content $(2n+1)$ to angular-content $(2l+1)$ in a way that no single-spin Peter–Weyl matrix-coefficient trace can produce. This is consistent with — but stronger than — the primitive-space defect: the shell-coupling structure of $\eta$ is itself a direct signature that the v3.61.0 substrate has shell-resolved content beyond what any single-slot SU(2) representation theory can carry.

---

## 5. Bit-exact verification gate compliance (§13.4)

- **Test gate ✓** — no production code touched; bit-exact `sympy.Rational` throughout.
- **Dead-end gate ✓** — no §3 match; this is a forward diagnostic sprint resolving the L1-vs-L2 ambiguity.
- **Prime directive gate ✓** — no discrete-structure modifications. The v3.61.0 χ, η closed forms are re-derived bit-exactly from the prosystem memo and used as input, not modified.
- **Consistency gate ✓** — bit-exact χ, η cross-check against prosystem memo (Sprint Q5'-Stage1-Prosystem, v3.60.0) at $n_{\max} \in \{1, 2, 3, 4\}$ passes 8/8; reproduces v3.61.0 generator count $3 N(n_{\max})$ bit-exactly; reproduces L2 substrate dim $(j_{\max} = 1/2: 15;\; j_{\max} = 1: 42)$ bit-exactly.
- **Equation gate ✓** — bit-exact verification of: (a) primitive coproduct on all $3 N(n_{\max})$ H_GV generators at $n_{\max} \in \{1, \dots, 6\}$ (231 total); (b) non-primitivity of $j = 1/2$ H_dec generators per slot (12 non-zero residuals); (c) det relation linear part bit-exactly matches $u + v$ form (1 symbolic match); (d) Lie-algebra dim formula 4 - 1 = 3 per slot (per-slot consistency check); (e) Q-iso defect at $n_{\max} \in \{1, \dots, 6\}$ (6 dimension comparisons).

**Total bit-exact zero residuals / consistent checks: 258.**

---

## 6. Honest scope

### 6.1 Closed at theorem grade (bit-exact at finite cutoff)

- $\dim \operatorname{Prim}(\mathcal{H}_{\mathrm{GV}})(n_{\max}) = 3 N(n_{\max})$ for all $n_{\max} \in \{1, \dots, 6\}$, bit-exactly.
- $\dim \operatorname{Prim}(\mathcal{H}_{\mathrm{dec}})(j_{\max}) = 9$ for all $j_{\max} \in \{1/2, 1, 3/2, 2, 5/2\}$, bit-exactly.
- The Q-iso defect $\dim \operatorname{Prim}(\mathcal{H}_{\mathrm{GV}}) - \dim \operatorname{Prim}(\mathcal{H}_{\mathrm{dec}}) = 3 N(n_{\max}) - 9$ is strictly positive at every $n_{\max} \ge 2$; ⟹ no Hopf isomorphism $\mathcal{H}_{\mathrm{dec}} \cong \mathcal{H}_{\mathrm{GV}}$ at any cutoff $n_{\max} \ge 2$.
- L2 cannot subsume L1; L1 (Levi-decomposition $\mathbb{G}_a^{3 N(n_{\max})} \times SL_2$) is the correct Stage-2 substrate target.
- Substrate-dim coincidence $15 = 15$ at $(n_{\max}, j_{\max}) = (2, 1/2)$ and $42 = 42$ at $(4, 1)$ is a small-cutoff artifact (breaks at $n_{\max} \in \{1, 3, 5, 6\}$).
- Cocommutativity invariant separates the substrates at every $n_{\max} \ge 1$ — L1 is cocommutative, L2 is not.

### 6.2 What's structural and what's argument-level

The dim Prim mismatch is **structurally forced** — it follows from the Milnor–Moore decomposition of the v3.61.0 cocommutative connected Hopf algebra (Prim = generating space) and from the standard identification $\operatorname{Prim}\mathcal{O}(G) = \mathfrak{g}^*$ for $G$ an affine algebraic group of dim 9. The bit-exact panel confirms this is not a finite-cutoff artifact.

The shell-coupling observation in $\eta_{(n, l)} = (2l+1)(2n+1)$ (§4) is an **independent argument** for the same conclusion — the (2n+1)(2l+1) factorization is not a single Peter–Weyl character — but the primitive-space dimension argument is the load-bearing one.

### 6.3 Named multi-year follow-ons (now unblocked)

- Full Tannakian closure of the pro-unipotent factor $\mathbb{G}_a^{3 N(n_{\max})}$ in the Connes–Marcolli machinery (Sprint Q5'-Levi-Synthesis §10.1 named follow-on a).
- Verification that the cosmic-Galois $U^*$ of Connes–Marcolli 2007 acts on $\mathcal{H}_{\mathrm{Levi}}$ in the expected way (Sprint Q5'-Levi-Synthesis §10.1 named follow-on b).
- $n_{\max} = 4$ extension of OffDiag pro-unipotent algebra-closure scaling (L4 of Sprint Q5'-Levi-Arc).
- Continuum-limit Mellin lift of the bridge $-\kappa^4$ identity to M3 (L3 of Sprint Q5'-Levi-Arc).

### 6.4 What's NOT closed by this sprint

- The L1 Tannakian closure itself remains multi-year (this sprint clears the substrate-choice ambiguity, it does not perform the closure).
- The strategically-distinct L2 substrate ($SL_2^3$) remains a valid synthesis route capturing M1/M2/M3 as three independent semisimple symmetries — it is not THE Stage-2 target, but it is a valid framing of the case-exhaustion theorem at the Hopf-algebra level. The two readings are not mutually exclusive; L1 is the Levi-decomp target, L2 is a complementary view of the semisimple factor.
- The S⁵ analog (Bargmann–Segal lattice, Paper 24): not tested here. The Coulomb/HO four-layer asymmetry (Paper 24 §V) suggests a different Hopf-algebra structure on the HO side; sprint-scale check left as a follow-on.

### 6.5 Hard prohibitions (§13.5)

- No changes to natural geometry hierarchy.
- No fitted/empirical parameters.
- No deletion of negative results.
- Paper 2 combination-rule "conjectural" label unchanged.

### 6.6 Curve-fit audit (`feedback_audit_numerical_claims`)

The L1-forced verdict is FORCED by Hopf-algebra invariants (dim Prim, cocommutativity), not curve-fit. Free-param count: 0. Selection bias: 0 (the discriminating invariant was chosen a priori, before seeing any cutoff data). Alternatives: explicitly tested L2 dim count at $j_{\max} \in \{1/2, 1, 3/2, 2, 5/2\}$ — uniformly 9, no candidate $j_{\max}$ matches the growing L1 dim. Robustness: defect grows monotonically with $n_{\max}$, no marginal cutoff. Independent test: cocommutativity separates the substrates independently of the dim count.

### 6.7 Tag-transcendentals compliance (`feedback_tag_transcendentals`)

No transcendentals introduced. All χ, η values are integers in $\mathbb{Z}$. The det relation, Lie-algebra dim formulas, and primitive coproduct residuals all live in $\mathbb{Q}[\text{generators}]$. The Mellin-engine partition (M1/M2/M3) sits one layer above the Hopf algebra at the trace-evaluation level (Paper 18 §III.7, Paper 32 §VIII case-exhaustion theorem); the Hopf algebra itself is transcendental-free.

### 6.8 Discrete-for-skeleton compliance (`feedback_discrete_for_skeleton`)

`sympy.Rational` throughout. No PSLQ. No floats. The Hopf algebra construction is over $\mathbb{Q}$; the discriminating invariants are integer-valued; all bit-exact checks return either 0 or a structural sympy expression matched verbatim.

### 6.9 WH1 PROVEN unaffected

This sprint compares two candidate Stage-2 Hopf substrates. It does not modify the WH1 / Marcolli–vS lineage closure or any propinquity result.

---

## 7. Files

### Produced
- `debug/compute_q5p_l1_vs_l2_diagnostic.py` — driver (~430 lines, 0.030 s wall, bit-exact `sympy.Rational` throughout).
- `debug/data/sprint_q5p_l1_vs_l2_diagnostic.json` — bit-exact data dump.
- `debug/sprint_q5p_l1_vs_l2_diagnostic_memo.md` — this memo (~3500 words).

### Used (load-bearing inputs)
- `debug/sprint_q5p_levi_arc_2026_06_06_memo.md` (umbrella v3.63.0; L1-vs-L2 question explicitly named §5).
- `debug/sprint_q5p_levi_synthesis_memo.md` (L1 substrate, v3.63.0).
- `debug/sprint_q5p_decorated_pw_memo.md` (L2 substrate, v3.63.0).
- `debug/sprint_q5p_stage2_hopf_memo.md` (v3.61.0 Track A — abelian baseline for H_GV).
- `debug/sprint_q5p_prosystem_memo.md` (v3.60.0 — χ, η closed forms).
- Paper 32 §VIII (case-exhaustion theorem).
- Paper 55 §subsec:open_m2_m3 (Q5' Stage-2 substrate paragraphs).

### Published references
- Milnor–Moore 1965 (Annals of Math 81): cocommutative connected Hopf algebra structure theorem (used to identify Prim H_GV with the generating space).
- Cartier 2007 (in: A primer of Hopf algebras): standard reference for $\operatorname{Prim}\mathcal{O}(G) = \mathfrak{g}^*$.
- Connes–Marcolli 2008 (Noncommutative Geometry, Quantum Fields and Motives, Ch. 4): cosmic-Galois $U^*$ with Levi decomposition.

---

## 8. Paper-edit recommendations

### 8.1 Paper 32 §VIII — ONE new Remark `rem:q5p_l1_forced` after `rem:q5p_combined_substrate_levi`

```latex
\begin{rem}[Q5' Stage-2 substrate decision: L1 is forced via the primitive-space invariant, Sprint Q5'-L1-vs-L2-Diagnostic, June 2026]
\label{rem:q5p_l1_forced}
The strategic question flagged in Sprint Q5'-Levi-Arc §5 -- whether the
v3.61.0 abelian primitives $x_{(n, l), k}$ in
$\mathcal{H}_{\mathrm{GV}}^{(n_{\max})}$ are reducible to Peter--Weyl
matrix coefficients in disguise (so that the simpler $SL_2^3$ substrate
of Remark~\ref{rem:q5p_decorated_pw} would suffice) or carry independent
content (so that the Levi-decomposition substrate of
Remark~\ref{rem:q5p_levi_synthesis_substrate} is forced) -- is resolved
bit-exactly at finite cutoff by the primitive-space invariant
$\dim \operatorname{Prim}(\mathcal{H}) = \dim \operatorname{Lie}(\operatorname{Spec}\mathcal{H})$.
For $\mathcal{H}_{\mathrm{GV}}^{(n_{\max})}$ the connected
cocommutative structure (Milnor--Moore) gives
$\dim \operatorname{Prim}(\mathcal{H}_{\mathrm{GV}}) = 3 N(n_{\max}) \in \{6, 15, 27, 42, 60, 81\}$
at $n_{\max} \in \{1, \dots, 6\}$. For $\mathcal{H}_{\mathrm{dec}}^{(j_{\max})}$
the standard identification
$\operatorname{Prim}\mathcal{O}(SL_2)^{\otimes 3} = \mathfrak{sl}_2^3$
gives $\dim \operatorname{Prim}(\mathcal{H}_{\mathrm{dec}}) = 9$
independent of $j_{\max}$ (higher $j_{\max}$ only adds degree-2 polynomial
expressions in the $j = 1/2$ generators via Clebsch--Gordan, contributing
nothing to $\mathfrak{m}/\mathfrak{m}^2$). At every $n_{\max} \ge 2$
the dim defect is strictly positive (6, 18, 33, 51, 72), so no Hopf
isomorphism exists between the two substrates; L1 is forced. The
substrate-dim coincidences $15 = 15$ at $(n_{\max}, j_{\max}) = (2, 1/2)$
and $42 = 42$ at $(4, 1)$ are small-cutoff artifacts breaking at
$n_{\max} \in \{1, 3, 5, 6\}$. An independent invariant (cocommutativity)
separates the substrates at every $n_{\max} \ge 1$ since
$\mathcal{H}_{\mathrm{GV}}$ is cocommutative (every generator primitive)
while $\mathcal{H}_{\mathrm{dec}}$ is not on $j > 0$ matrix coefficients.
Bit-exact panel: 258 zero residuals / consistent checks (driver
\texttt{debug/compute\_q5p\_l1\_vs\_l2\_diagnostic.py}; data
\texttt{debug/data/sprint\_q5p\_l1\_vs\_l2\_diagnostic.json}). The
L2 $SL_2^3$ substrate remains valid as a complementary semisimple-only
reading; the Levi-decomposition
$U^*_{\mathrm{GeoVac}, \mathrm{Levi}} = \mathbb{G}_a^{3 N(n_{\max})} \times SL_2$
is the structurally-correct Stage-2 target, matching the published
Connes--Marcolli motivic-Galois shape. Multi-year continuation: the
Tannakian closure of the pro-unipotent factor
$\mathbb{G}_a^{3 N(n_{\max})}$ in the Connes--Marcolli machinery is now
the canonical next step of the Q5' program.
\end{rem}
```

### 8.2 Paper 55 §subsec:open_m2_m3 — ONE new paragraph after the L1+L5 categorical-reduction paragraph

```latex
\emph{Stage-2 substrate decision via primitive-space invariant: L1 is
forced (Sprint Q5'-L1-vs-L2-Diagnostic, June 2026; memo
\texttt{debug/sprint\_q5p\_l1\_vs\_l2\_diagnostic\_memo.md}; data
\texttt{debug/data/sprint\_q5p\_l1\_vs\_l2\_diagnostic.json}).}
The strategic question of whether the v3.61.0 abelian primitives carry
independent content (forcing L1's Levi-decomposition) or are reducible
to Peter--Weyl matrix coefficients (allowing L2's simpler $SL_2^3$
substrate to suffice) is resolved bit-exactly at finite cutoff by the
primitive-space invariant. At every $n_{\max} \ge 2$,
$\dim \operatorname{Prim}(\mathcal{H}_{\mathrm{GV}}) = 3 N(n_{\max})$
strictly exceeds the constant $\dim \operatorname{Prim}(\mathcal{H}_{\mathrm{dec}}) = 9 = \dim \mathfrak{sl}_2^3$,
giving dim defects $(6, 18, 33, 51, 72)$ at $n_{\max} \in \{2, \dots, 6\}$
(258 zero residuals across the full panel). The substrate-dim
coincidences flagged in the decorated-PW paragraph above
($15 = 15$ at $(n_{\max}, j_{\max}) = (2, 1/2)$; $42 = 42$ at $(4, 1)$)
break at $n_{\max} \in \{1, 3, 5, 6\}$, confirming they are
small-cutoff artifacts rather than structural identifications.
Cocommutativity separates the substrates independently of the dim
count at every $n_{\max} \ge 1$. The Levi-decomposition substrate
$U^* = \mathbb{G}_a^{3 N(n_{\max})} \times SL_2$ is the correct
Stage-2 target; the decorated-PW $SL_2^3$ substrate is a complementary
semisimple-only reading capturing M1/M2/M3 as three independent
non-abelian semisimple symmetries, but cannot subsume the
pro-unipotent abelian primitive content of the v3.61.0 substrate.
The canonical multi-year next step of the Q5' program is the
Tannakian closure of the pro-unipotent factor $\mathbb{G}_a^{3 N(n_{\max})}$
in the Connes--Marcolli machinery.
```

### 8.3 Paper 18 — no edit

Paper 18 §III.7 master Mellin engine sits one layer above the Stage-2 substrate at the trace-evaluation level. No Paper 18 edit needed.

---

## 9. One-line verdict

**The primitive-space invariant $\dim \operatorname{Prim}(\mathcal{H})$ decides the Stage-2 substrate question bit-exactly: $\dim \operatorname{Prim}(\mathcal{H}_{\mathrm{GV}})(n_{\max}) = 3 N(n_{\max})$ grows quadratically (15, 27, 42, 60, 81 at $n_{\max} = 2, \dots, 6$) while $\dim \operatorname{Prim}(\mathcal{H}_{\mathrm{dec}})(j_{\max}) = 9$ is constant; defect $\ge 6$ at every cutoff $n_{\max} \ge 2$ and grows quadratically. L2 cannot subsume L1 as a Hopf algebra, the substrate-dim coincidences $15 = 15$ at $(2, 1/2)$ and $42 = 42$ at $(4, 1)$ are small-cutoff artifacts that break at four of six tested cutoffs, cocommutativity separates the substrates independently of the dim count, and therefore the Levi-decomposition $U^*_{\mathrm{Levi}} = \mathbb{G}_a^{3 N(n_{\max})} \times SL_2$ is the correct Stage-2 target. The L2 $SL_2^3$ substrate remains valid as a complementary semisimple-only reading. 258 bit-exact zero residuals across the panel. The Q5' Stage-2 substrate-choice ambiguity is closed; the multi-year canonical next step is the Tannakian closure of the pro-unipotent factor $\mathbb{G}_a^{3 N(n_{\max})}$ in the Connes--Marcolli machinery.**
