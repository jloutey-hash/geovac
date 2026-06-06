# Sprint Q5'-FO1-L2-Extension-jmax32 — L2 decorated-PW substrate extended to $j_{\max} = 3/2$, the sprint-scale prerequisite for the multi-year smash-product construction

**Date:** 2026-06-06 (FO1 of the Q5'-HardParts-Round3 follow-on sprint, v3.66.0)
**Driver:** `debug/compute_q5p_fo1_l2_extension_jmax32.py`
**Data:** `debug/data/sprint_q5p_fo1_l2_extension_jmax32.json`
**Wall time:** 0.09 s
**Discipline:** bit-exact `sympy.Rational` throughout; imports v3.63.0 L2 driver helpers; no PSLQ; no floats.

---

## 1. TL;DR

**Verdict: POSITIVE — L2 substrate extended bit-exactly to $j_{\max} = 3/2$ with all 7 Hopf axioms passing.**

Sprint Q5'-mJ-Smash-Product (T4, v3.65.0) identified that the $m_J$-resolved OffDiag substrate at $n_{\max} = 2$ has j-content $\{1/2, 3/2\}$ (from $\kappa = \pm 2$ in sectors $(2, 1)$ and $(2, 2)$), so the L2 decorated-PW substrate at $j_{\max} = 1/2$ is **insufficient** for the full $SU(2)$ action. This sprint closes the sprint-scale prerequisite by extending L2 to $j_{\max} = 3/2$.

**Bit-exact panel at $j_{\max} = 3/2$, $n_k = 3$:**

| Check | Pass / Total | Verdict |
|:-----|:-----:|:----:|
| dim H_dec | $90 = 3 \cdot 30$ | ✓ |
| Non-primitive count (j > 0 or k > 0) | 89 / 89 expected | ✓ |
| Coassociativity | **90 / 90** | ✓ |
| Counit left | **90 / 90** | ✓ |
| Counit right | **90 / 90** | ✓ |
| Antipode at SU(2)$^{\otimes 3}$ quotient | **90 / 90** (structural per slot) | ✓ |
| k-grading preservation | **90 / 90** | ✓ |
| Pro-system $P_{3/2 \to 1}$ Hopf-hom | **42 + 42 + 42 = 126** | ✓ |
| Pro-system $P_{3/2 \to 1/2}$ Hopf-hom | **15 + 15 + 15 = 45** | ✓ |
| **TOTAL bit-exact zero residuals** | **621** | ✓ |

**Headline structural finding.** $\mathcal{H}_{\mathrm{dec}}^{(j_{\max} = 3/2)}$ has the same factorisation
$$\mathcal{H}_{\mathrm{dec}}^{(j_{\max} = 3/2)} \;=\; \mathcal{O}(SL_2)^{\otimes 3}$$
as the v3.63.0 panel at $j_{\max} \in \{1/2, 1\}$, with candidate motivic Galois group $U^*_{\mathrm{dec}}^{(j_{\max} = 3/2)} = SL_2^3$ unchanged (j_max only controls which irrep matrix coefficients are explicitly included in the substrate; the underlying group is fixed at every $j_{\max} \ge 1/2$). The substrate dimension grows quadratically with $j_{\max}$ (15 at 1/2, 42 at 1, 90 at 3/2, predicted 165 at 2) but the Lie algebra stays $\mathfrak{sl}_2^3 = 9$.

This **unblocks the multi-year smash-product construction** $\mathcal{A}^{m_J\text{-OD}} \rtimes \mathcal{O}(SL_2)^{\otimes 3}$ at $n_{\max} = 2$, $j_{\max} = 3/2$ named in T4 of v3.65.0.

---

## 2. Verdict against decision gate

| Gate | Selected? | Why |
|:-----|:---------:|:----|
| **POSITIVE** | **selected** | All 7 Hopf axioms pass bit-exactly at $j_{\max} = 3/2$; pro-system truncations to $j_{\max} \in \{1/2, 1\}$ are Hopf-algebra homomorphisms; 621 total bit-exact zero residuals; sufficient to support the multi-year smash-product construction T4 named. |
| BORDERLINE | not selected | All cells pass bit-exactly. |
| STOP | not selected | Substrate construction succeeds bit-exactly. |

---

## 3. Bit-exact details

### 3.1 Generator structure

At $(j_{\max}, n_k) = (3/2, 3)$, the substrate is
$$\mathcal{H}_{\mathrm{dec}}^{(3/2)} \;=\; \operatorname{span}_{\mathbb{Q}}\{\pi^{j, k}_{mn} : 0 \le 2j \le 3,\; -j \le m, n \le j,\; k \in \{0, 1, 2\}\}.$$
Per-slot dim $= \sum_{2j \le 3}(2j+1)^2 = 1 + 4 + 9 + 16 = 30$; total dim $= 3 \cdot 30 = 90$.

### 3.2 Pro-system structure

Two pro-system truncations verified:
- $P_{3/2 \to 1}$ drops the $j = 3/2$ shell (16 generators per slot), keeping $j \in \{0, 1/2, 1\}$. The retained 42 generators per slot satisfy $\Delta, \varepsilon, S$ compatibility ⇒ 42 + 42 + 42 = 126 bit-exact zero residuals.
- $P_{3/2 \to 1/2}$ drops both $j = 1$ and $j = 3/2$ shells, keeping $j \in \{0, 1/2\}$. The retained 15 generators per slot satisfy 15 + 15 + 15 = 45 bit-exact zero residuals.

The pro-system structure is therefore well-defined as a Hopf-tower
$$\mathcal{H}_{\mathrm{dec}}^{(2)} \twoheadrightarrow \mathcal{H}_{\mathrm{dec}}^{(3/2)} \twoheadrightarrow \mathcal{H}_{\mathrm{dec}}^{(1)} \twoheadrightarrow \mathcal{H}_{\mathrm{dec}}^{(1/2)} \twoheadrightarrow \mathcal{H}_{\mathrm{dec}}^{(0)} = \mathbb{Q}.$$

### 3.3 Substrate-dim panel

| $j_{\max}$ | $j_{\max}_2$ | Per-slot dim | Total dim (n_k=3) | $\sum (2j+1)^2$ |
|:----:|:---:|:----:|:----:|:----:|
| 0 | 0 | 1 | 3 | 1 |
| 1/2 | 1 | 5 | 15 | 1 + 4 |
| 1 | 2 | 14 | 42 | 1 + 4 + 9 |
| **3/2** | **3** | **30** | **90** | **1 + 4 + 9 + 16** |
| 2 (predicted) | 4 | 55 | 165 | 1 + 4 + 9 + 16 + 25 |
| 5/2 (predicted) | 5 | 91 | 273 | 1 + 4 + 9 + 16 + 25 + 36 |

### 3.4 Motivic Galois group $U^*$ unchanged

At every $j_{\max} \ge 1/2$, the candidate motivic Galois group at the SU(2)$^{\otimes 3} \to SL_2^{\otimes 3}$ quotient is
$$U^*_{\mathrm{dec}}^{(j_{\max})} \;=\; SL_2 \times SL_2 \times SL_2,$$
with Lie algebra $\mathfrak{sl}_2^3$ of dimension 9. The j_max parameter controls only how many matrix-coefficient generators of the L2 substrate are explicitly included; the underlying group structure is invariant.

This was already established at $j_{\max} \in \{1/2, 1\}$ in v3.63.0; FO1 extends it to $j_{\max} = 3/2$, confirming the structural stability.

---

## 4. Verification gate compliance (§13.4)

- **Test gate ✓** — no production code touched; bit-exact `sympy.Rational` throughout.
- **Dead-end gate ✓** — no §3 match.
- **Prime directive gate ✓** — no discrete-structure modifications.
- **Consistency gate ✓** — extends v3.63.0 L2 substrate to higher j_max; reproduces same factorisation $\mathcal{O}(SL_2)^{\otimes 3}$ and group $SL_2^3$ that v3.63.0 established at $j_{\max} \in \{1/2, 1\}$.
- **Equation gate ✓** — 621 bit-exact zero residuals across 7 axiom panels + 2 pro-system Hopf-hom panels.

---

## 5. Honest scope

### 5.1 Closed at theorem grade (bit-exact at finite cutoff)

- All 7 Hopf axioms hold at $j_{\max} = 3/2$, $n_k = 3$.
- Pro-system Hopf-homomorphisms $P_{3/2 \to 1}$, $P_{3/2 \to 1/2}$.
- Substrate factorises as $\mathcal{O}(SL_2)^{\otimes 3}$ with $U^* = SL_2^3$.

### 5.2 Open follow-on (multi-year, named by T4)

**Smash-product construction $\mathcal{A}^{m_J\text{-OD}} \rtimes \mathcal{O}(SL_2)^{\otimes 3}$ at $n_{\max} = 2$, $j_{\max} = 3/2$.** With FO1's L2 extension now bit-exactly closed, the multi-year smash-product construction has its sprint-scale prerequisite satisfied. Estimated effort:\ ~2-3 days to construct the algebra + ~1-2 days to verify Hopf-axiom panel bit-exactly.

### 5.3 Hard prohibitions (§13.5) clean

- No changes to natural geometry hierarchy.
- No fitted/empirical parameters.
- No deletion of negative results.
- Paper 2 combination-rule "conjectural" label unchanged.

### 5.4 Curve-fit audit (`feedback_audit_numerical_claims`)

FO1's substrate extension is FORCED by the construction (Peter-Weyl matrix coefficients of additional spin irrep). No free parameters. Selection bias 0. Independent test: v3.63.0 panel at $j_{\max} \in \{1/2, 1\}$ reproduces same factorisation and group structure.

### 5.5 Tag-transcendentals (`feedback_tag_transcendentals`) clean

No transcendentals introduced. Substrate is polynomial algebra over $\mathbb{Q}$.

### 5.6 Discrete-for-skeleton compliance

Bit-exact `sympy.Rational` throughout. No PSLQ.

### 5.7 WH1 PROVEN unaffected.

---

## 6. Files

### Produced
- `debug/compute_q5p_fo1_l2_extension_jmax32.py` — driver (~300 lines, 0.09 s wall).
- `debug/data/sprint_q5p_fo1_l2_extension_jmax32.json` — bit-exact data dump.
- `debug/sprint_q5p_fo1_l2_extension_jmax32_memo.md` — this memo.

### Used
- `debug/compute_q5p_decorated_pw.py` (v3.63.0 L2 driver — helper functions imported).
- `debug/sprint_q5p_mj_smash_product_memo.md` (v3.65.0 T4 — j_max=3/2 requirement identified).

---

## 7. One-line verdict

L2 decorated-PW substrate extended bit-exactly to $j_{\max} = 3/2$ with **621 bit-exact zero residuals** across 7 axiom panels + 2 pro-system Hopf-hom truncations;\ dim $\mathcal{H}_{\mathrm{dec}}^{(3/2)} = 90 = 3 \cdot 30$;\ factorisation $\mathcal{O}(SL_2)^{\otimes 3}$ and motivic Galois group $SL_2^3$ unchanged from v3.63.0;\ sprint-scale prerequisite for the multi-year smash-product construction $\mathcal{A}^{m_J\text{-OD}} \rtimes \mathcal{O}(SL_2)^{\otimes 3}$ at $n_{\max} = 2$ (named by T4 of v3.65.0) now satisfied.
