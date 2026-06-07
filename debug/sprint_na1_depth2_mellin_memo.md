# Sprint NA-1 depth-2 Mellin test

**Date:** 2026-06-06 (evening, follow-up to Hain-Brown PSLQ negative)
**Sprint:** NA-1 depth-2 Mellin test — Reading A (primitive / abelianisation) vs Reading B (shuffle / free non-abelian) disambiguation; opportunistic Hain–Brown probe at depth 2
**Driver:** `debug/sprint_na1_depth2_mellin_compute.py`
**Data:** `debug/data/na1_depth2_mellin_results.json`
**Wall time:** 38 s (closed-form M2 via Paper 28 Thm 1 + mpmath alternating-series acceleration; PSLQ at 50/100/200 dps).
**Discipline:** `sympy` for symbolic closed-form sanity; `mpmath.nsum(method='alternating')` for slow alternating series (J_eff(s_tot=3) is conditionally convergent at the order ∼ 1/n); PSLQ at three precisions with cross-precision agreement filter; algebraic identities verified termwise on the closed-form Dirac S³ spectrum.

---

## 1. TL;DR

**Verdict against the original decision gate:** **NEITHER Reading A nor Reading B** is supported on the joint depth-2 Mellin of an M2 × M3 pair built on the Camporesi–Higuchi Dirac S³ substrate. A **third structural reading (Reading C — diagonal-collapse)** is bit-exactly forced instead:

> **Theorem (NA-1).** Let $D$ be the Camporesi–Higuchi Dirac on $S^3$ and $\gamma_P : |n,\ldots\rangle \mapsto (-1)^n |n,\ldots\rangle$ the vertex-parity grading on the eigenshells. For every integer $s_{\mathrm{tot}} \ge 3$,
> $$J(s_1, s_2) := \frac{1}{\Gamma(s_1)\Gamma(s_2)} \int_0^\infty\!\!\!\int_0^\infty t_1^{s_1-1} t_2^{s_2-1}\, \mathrm{Tr}\big(D^2 e^{-t_1 D^2} \cdot \gamma_P D\, e^{-t_2 D^2}\big)\, dt_1 dt_2$$
> depends only on $s_{\mathrm{tot}} = s_1 + s_2$ and collapses bit-exactly to a depth-1 parity-graded Mellin moment at shifted argument:
> $$J(s_1, s_2) \;=\; J_{\mathrm{eff}}(s_{\mathrm{tot}}) \;=\; M_3^{\gamma_P}(s_{\mathrm{tot}} - 1),$$
> where $M_3^{\gamma_P}(s) := \sum_{n \ge 0}(-1)^n (n+1)(n+2)\big((2n+3)/2\big)^{1-2s}$.

The depth-2 joint Mellin is a depth-1 object in disguise on the CH diagonal substrate. The deconcatenation-pair vs primitive-product diagnostic is **structurally invisible** on this substrate because $D^2$, $\gamma_P$, $D$, and $e^{-t D^2}$ are all simultaneously diagonal in the CH eigenbasis — operator ordering is trivial, and the joint trace algebraically degenerates to a single power $\lambda_n^{3 - 2 s_{\mathrm{tot}}}$.

**Hain–Brown at depth 2:** **NEGATIVE.** No GeoVac depth-2 period identifies with $E_4(2i)$, $E_6(2i)$, $\Delta(2i)$, or any Eichler period of $\Delta$ across all three precisions. Yesterday's HB negative at depth 1 is structurally inherited at depth 2 on this substrate.

**Additional substantive finding.** PSLQ-identification of the depth-1 sector $M_3^{\gamma_P}(s)$ at integer $s$ yields the closed form
$$M_3^{\gamma_P}(s) \;=\; 2^{2s-3}\big(\beta(2s-1) - \beta(2s-3)\big),$$
with $\beta(2k+1) \in \pi^{2k+1}\mathbb{Q}$ (Euler-number closed form for Dirichlet $\beta$ at odd integers). The "$M_3$" label on the CH parity-grading is therefore **pure-Tate (level 1)**, not level-4 cyclotomic. Paper 28's level-4 $M_3$ (Hurwitz at quarter-integer shift) is a **different operator structure** (vertex-restricted, not parity-graded-diagonal). The CH parity grading and the QED vertex parity are not the same M3 mechanism; the natural CH parity grading produces pure-Tate output.

**Decision-gate outcome (refined).** The original A/B/INCONCLUSIVE gate underspecified the geometry. The honest verdict is:

- **NEITHER on the diagonal CH substrate** — the substrate is *structurally insensitive* to the A vs B distinction because of operator simultaneous diagonalisation.
- **The test is non-vacuous only on a non-diagonal substrate** (e.g. the off-diagonal CH used in WH1's R3.5 / Sprint TS work, or a vertex-restricted QED graph where γ_P does NOT commute with the kernel).
- **Reading A vs Reading B is the correct dichotomy** but **requires a non-diagonal substrate to test**; running it on CH-diagonal is the equivalent of asking whether two scalars commute.

**Strategic implication.** The NA-1 question is correctly identified by yesterday's structural sprint (`debug/sprint_q5p_na1_non_abelian_probe_memo.md`) and the strategic synthesis (`debug/strategic_synthesis_2026_06_06_memo.md` §2 NA-1, §6 Recommendation A). The **right next sprint** is to repeat NA-1 on a **non-diagonal substrate** — either the off-diagonal CH operator system of WH1 R3.5 (where chirality-flipping E1 entries break simultaneous diagonalisation) or the Paper 28 QED vertex graph (where γ_P is realised by combinatorial vertex parity, not by the diagonal $(-1)^n$). This is a **2–4 week sprint**, comparable to the original NA-1 depth-2 ambition.

---

## 2. Setup

### 2.1 The joint depth-2 Mellin

Per Paper 32 §VIII + Paper 18 §III.7, the master Mellin engine identifies three sub-mechanisms via the index $k$ in
$$\mathcal{M}_k[f](s) := \frac{1}{\Gamma(s)}\int_0^\infty t^{s-1}\,\mathrm{Tr}(D^k\, f(D^2))\, dt,$$
with $k = 0$ (M1, Hopf base measure), $k = 1$ (M3, vertex parity / Hurwitz), $k = 2$ (M2, Seeley–DeWitt heat kernel).

The natural **depth-2 generalisation** considered by NA-1 is the joint moment over an M2 slot × M3 slot:
$$J(s_1, s_2) := \mathcal{M}_{k_1,\,k_2}\big[D^{k_1} e^{-t_1 D^2} \cdot \gamma_P\, D^{k_2} e^{-t_2 D^2}\big]$$
with $(k_1, k_2) = (2, 1)$. The vertex-parity grading $\gamma_P$ on the CH eigenshells is taken as the diagonal operator $\gamma_P |n, \chi, m_l\rangle = (-1)^n |n, \chi, m_l\rangle$.

### 2.2 Bit-exact spectrum

Camporesi–Higuchi Dirac on $S^3$:
- eigenvalues: $|\lambda_n| = (2n + 3)/2$, $n = 0, 1, 2, \ldots$
- degeneracies: $g_n = (n+1)(n+2)$ (sum over both chiralities + $m_l$ multiplicity)
- parity sign: $s_n = (-1)^n$

Operator traces:
$$\mathrm{Tr}(D^k e^{-t D^2}) = \sum_n g_n\, \lambda_n^k\, e^{-t \lambda_n^2}, \qquad \mathrm{Tr}(\gamma_P D^k e^{-t D^2}) = \sum_n s_n\, g_n\, \lambda_n^k\, e^{-t \lambda_n^2}.$$

### 2.3 The collapse

Because $D^2$, $D$, $\gamma_P$, and $e^{-t D^2}$ are simultaneously diagonal in the CH eigenbasis (where $D|n\rangle = \lambda_n |n\rangle$ and $\gamma_P|n\rangle = (-1)^n |n\rangle$), the joint operator trace algebraically degenerates:
$$\mathrm{Tr}\big(D^2 e^{-t_1 D^2}\, \gamma_P\, D\, e^{-t_2 D^2}\big) = \sum_n (-1)^n g_n\, \lambda_n^3\, e^{-(t_1+t_2)\lambda_n^2}.$$

Mellin-double-transforming:
$$J(s_1, s_2) = \sum_n (-1)^n g_n\, \lambda_n^3 \cdot \lambda_n^{-2 s_1} \cdot \lambda_n^{-2 s_2} = \sum_n (-1)^n g_n\, \lambda_n^{3 - 2 s_{\mathrm{tot}}},$$
where $s_{\mathrm{tot}} := s_1 + s_2$. This is a function of $s_{\mathrm{tot}}$ ONLY — the depth-2 label is decorative.

### 2.4 The depth-1 collapse identity

Comparing to the depth-1 M3-like parity-graded Mellin
$$M_3^{\gamma_P}(s) := \sum_n (-1)^n g_n\, \lambda_n^{1 - 2s},$$
the exponents agree iff $3 - 2 s_{\mathrm{tot}} = 1 - 2(s_{\mathrm{tot}} - 1)$ — which holds for ALL $s_{\mathrm{tot}}$. Therefore
$$\boxed{\quad J(s_1, s_2) \;=\; M_3^{\gamma_P}(s_{\mathrm{tot}} - 1) \quad}$$
**as a termwise algebraic identity** of the closed-form spectral sums.

---

## 3. Numerical verification

### 3.1 Cross-precision panel (50, 100, 200 dps)

| $s_{\mathrm{tot}}$ | $J_{\mathrm{eff}}(s_{\mathrm{tot}})$ | $M_3^{\gamma_P}(s_{\mathrm{tot}} - 1)$ | $|J - M_3|$ |
|:---:|:---|:---|:---|
| 3 | 0.367095965723842141735948 | 0.367095965723842141735948 | 0.0 (bit-exact at 200 dps) |
| 4 | 0.217693454541749468181476 | 0.217693454541749468181476 | 0.0 |
| 5 | 0.108693754030459055680869 | 0.108693754030459055680869 | 0.0 |
| 6 | 0.050582565975063081601577 | 0.050582565975063081601577 | 0.0 |
| 7 | 0.022881682741047986085672 | 0.022881682741047986085672 | 0.0 |

(Identity proved termwise on the closed-form series; numerical verification at $200$ dps confirms bit-exact agreement.)

### 3.2 PSLQ panel at 50, 100, 200 dps

For each $s_{\mathrm{tot}} \in \{3, \ldots, 7\}$, PSLQ on $J_{\mathrm{eff}}(s_{\mathrm{tot}})$ against a basis $\{M_2(s_1)M_3(s_2)$ for splits $s_1 + s_2 = s_{\mathrm{tot}}\}$ $\cup\, \{M_2(s),\, M_3(s)\}$ $\cup\,\{\pi^k\}_{k=1}^{10}$ $\cup\, \{\zeta(3), \zeta(5), G, \beta(4)\}$ $\cup\, \{E_4(2i), E_6(2i), \Delta(2i), \mathrm{Eichler}\ P_{0,1,2}\}$ at coefficient ceiling $10^6$, `maxsteps = 2000`. All three precisions return the same signature, agreeing across the cross-precision filter:

| $s_{\mathrm{tot}}$ | Cross-prec agreement | HB modular content | PSLQ relation (at 200 dps) |
|:---:|:---:|:---:|:---|
| 3 | **agree (3/3)** | **NONE** | $J_{\mathrm{eff}}(3) - M_3(2) = 0$ |
| 4 | **agree (3/3)** | **NONE** | $J_{\mathrm{eff}}(4) - M_3(3) = 0$ |
| 5 | **agree (3/3)** | **NONE** | $J_{\mathrm{eff}}(5) - M_3(4) = 0$ |
| 6 | **agree (3/3)** | **NONE** | $J_{\mathrm{eff}}(6) - M_3(5) = 0$ |
| 7 | **agree (3/3)** | **NONE** | $J_{\mathrm{eff}}(7) - M_3(6) = 0$ |

Every PSLQ-identification at every precision picked out the same single-basis-element relation $J_{\mathrm{eff}}(s_{\mathrm{tot}}) = M_3(s_{\mathrm{tot}} - 1)$ exactly. The collapse identity is what PSLQ finds, not the primitive product or the shuffle pair.

### 3.3 Factorisation diff panel (200 dps)

For each split $(s_1, s_2)$, the differences $J_{\mathrm{eff}} - M_2(s_1)\cdot M_3(s_2)$ (primitive) and $J_{\mathrm{eff}} - (M_2(s_1)M_3(s_2) + M_3(s_1)M_2(s_2))$ (shuffle) are all $O(1)$, not zero:

| $(s_1, s_2)$ | $J_{\mathrm{eff}} - \mathrm{primitive}$ | $J_{\mathrm{eff}} - \mathrm{shuffle}$ |
|:---:|:---:|:---:|
| $(2, 1)$ | $-0.321$ | (s2 < 2, undefined) |
| $(2, 2)$ | $-0.426$ | $-1.07$ |
| $(3, 1)$ | $+0.0514$ | (s2 < 2, undefined) |
| $(3, 2)$ | $-0.0467$ | $-0.428$ |
| $(4, 2)$ | $-0.0101$ | $-0.201$ |
| $(3, 3)$ | $-0.0416$ | $-0.134$ |
| $(4, 3)$ | $-0.0131$ | $-0.0591$ |
| $(5, 2)$ | $-7.75$ | $-7.84$ |

Conclusion: **zero primitive hits, zero shuffle hits, zero half-shuffle hits** — neither Reading A nor Reading B is empirically supported on the joint depth-2 Mellin under the CH-diagonal substrate.

### 3.4 Sub-finding: pure-Tate closed form for $M_3^{\gamma_P}(s)$

Algebraic derivation (verified bit-exact at 60 dps):
$$M_3^{\gamma_P}(s) = \sum_n (-1)^n [(n+3/2)^2 - 1/4]\,(n+3/2)^{1-2s} = \eta_H(2s-3, 3/2) - \tfrac{1}{4}\eta_H(2s-1, 3/2),$$
where $\eta_H(s, a) := \sum_n (-1)^n (n+a)^{-s}$. With the identity $\eta_H(s, 3/2) = 2^s (1 - \beta(s))$ (substituting $u = n + 3/2$ and reorganising the alternating Hurwitz at half-integer shift through the standard Dirichlet-$\beta$ identity for odd-integer shifts), one obtains the closed form
$$\boxed{\quad M_3^{\gamma_P}(s) = 2^{2s-3}\big(\beta(2s-1) - \beta(2s-3)\big). \quad}$$

Bit-exact verification at $s \in \{2, 3, 4\}$:
- $M_3^{\gamma_P}(2) = 2(\beta(3) - \beta(1)) = 2(\pi^3/32) - 2(\pi/4) = \pi^3/16 - \pi/2$ ✓
- $M_3^{\gamma_P}(3) = 8(\beta(5) - \beta(3)) = 8(5\pi^5/1536) - 8(\pi^3/32) = 5\pi^5/192 - \pi^3/4$ ✓
- $M_3^{\gamma_P}(4) = 32(\beta(7) - \beta(5)) = 32(61\pi^7/184320) - 32(5\pi^5/1536) = 61\pi^7/5760 - 5\pi^5/48$ ✓

Since $\beta(2k+1) = \frac{(-1)^k E_{2k}}{2(2k)!}\,(2\pi)^{2k+1} \in \pi^{2k+1}\mathbb{Q}$ (Euler-number closed form), every $M_3^{\gamma_P}(s)$ at integer $s$ lies in $\pi\cdot\mathbb{Q} \oplus \pi^3\cdot\mathbb{Q} \oplus \pi^5\cdot\mathbb{Q} \oplus \cdots$ — strictly **pure-Tate at odd weights**. This is a different period ring than Paper 28's level-4 cyclotomic $M_3$ sector (which uses $\zeta(s, 1/4) \pm \zeta(s, 3/4)$ and lives in $\beta(\mathrm{even})\cdot\mathbb{Q}$, the Catalan / β(4) ring).

**The "M3 sector" splits.** Two structurally distinct $M_3$ realisations exist in the GeoVac framework:

1. **Pure-Tate $M_3$ — diagonal parity grading on CH.** $M_3^{\gamma_P}(s) = 2^{2s-3}(\beta(2s-1) - \beta(2s-3))$, lives in $\bigoplus_k \pi^{2k+1}\mathbb{Q}$. Always pure-Tate, never modular, never level-4 cyclotomic.
2. **Level-4 cyclotomic $M_3$ — Hurwitz at quarter-integer shift (Paper 28 Theorem 3, vertex-parity discriminant $D_{\mathrm{even}}(s) - D_{\mathrm{odd}}(s) = 2^{s-1}(\beta(s) - \beta(s-2))$).** Lives in the Catalan / $\beta(4)$ extension, i.e. $\mathrm{MT}(\mathbb{Z}[i, 1/2])$.

Paper 55 §5 implicitly assumed all $M_3$ goes to category (2). The NA-1 diagonal substrate makes category (1) explicit and shows the CH parity grading **cannot reach** the level-4 cyclotomic content — that sits one step away on a different operator structure.

---

## 4. Interpretation against the decision gate

The original NA-1 strategic-synthesis decision gate (`debug/sprint_q5p_na1_non_abelian_probe_memo.md` §3 Part D; `debug/strategic_synthesis_2026_06_06_memo.md` §6 Recommendation A) was:

- **READING A WINS** if symmetric factorisation in pure MZV $\Rightarrow$ GeoVac IS the abelianization of cosmic Galois.
- **READING B WINS** if asymmetric (deconcatenation-pair) factorisation $\Rightarrow$ substrate needs shuffle Hopf enrichment.
- **INCONCLUSIVE** if joint Mellin not cleanly defined without depth-2 generalisation of case-exhaustion.
- **HAIN–BROWN at depth 2** if modular content surfaces in $J_{\mathrm{eff}}$.

This sprint adds a fourth alternative the gate did not anticipate:

> **READING C — diagonal collapse.** The joint depth-2 Mellin **algebraically degenerates** to a depth-1 object via simultaneous diagonalisation of $D^2, D, \gamma_P, e^{-t D^2}$ on the CH eigenbasis. The substrate is structurally **insensitive to the A vs B distinction** — both Readings predict depth-2 product structure that the collapse identity rules out term-by-term, regardless of which depth-2 product Hopf-coproduct one favours.

The honest reading is **Reading C is forced bit-exactly on the CH-diagonal substrate**, and **distinguishing A from B requires a different substrate**.

### 4.1 What this implies for the cosmic-Galois target

Three things, ordered from least to most surprising:

1. **The pure-Tate result confirms Paper 55 §6**: the *natural* joint M2 × M3 on the CH substrate factors trivially. The level-4 cyclotomic content does NOT enter via this depth-2 construction — it would only enter via the Paper 28 vertex-restricted M3 mechanism, which is operator-different.

2. **The collapse identity sharpens the master Mellin engine partition.** Paper 18 §III.7 states $k \in \{0, 1, 2\}$ classifies M1/M3/M2. NA-1 strengthens this: on the diagonal substrate, $k_1 + k_2$ at depth 2 maps to $k = k_1 + k_2 - 1$ at depth 1 via the algebraic collapse $D^{k_1} \gamma_P D^{k_2} = \gamma_P D^{k_1 + k_2}$ on the eigenbasis (since $\gamma_P$ commutes with $D$ on simultaneous eigenstates). **Depth-2 in the master Mellin engine is decorative on a diagonal substrate; the genuine depth-2 content lives off-diagonal.**

3. **The A vs B test is well-posed but mis-located on this substrate.** The substrate that respects the Hopf-coproduct primitivity-vs-shuffle distinction must have a non-trivial commutator structure between the M2 and M3 operators. On CH-diagonal, $[D^2, \gamma_P D] = 0$ identically (both diagonal). On the off-diagonal CH (WH1 R3.5), $\gamma_P$ no longer commutes with $D$ because the chirality-flipping E1 entries connect even-$n$ and odd-$n$ sectors — a genuine commutator emerges and the depth-2 joint trace ceases to collapse.

### 4.2 What does NOT change

- WH1 PROVEN unchanged.
- WH6 (BC-RH wall) unchanged.
- Paper 55 / 56 finite-cutoff Tannakian closure unchanged.
- The strategic-synthesis Recommendation B (explicit $U^*_{GV} \hookrightarrow \mathcal{U}_4$ injection, 2–3 weeks) is **strengthened**: NA-1's bit-exact collapse to depth-1 M3 + pure-Tate closed form of $M_3^{\gamma_P}$ confirms the GeoVac substrate's period content sits in $\mathrm{MT}(\mathbb{Z}[i, 1/2])$ at depth $\le 1$ as Paper 55 §5 classifies; depth-2 does not add new content via the CH-diagonal substrate.

### 4.3 The five-sentence diagnostic

> The NA-1 depth-2 Mellin test is the right test, run on the wrong substrate. On the CH-diagonal substrate, $D^2$ and $\gamma_P D$ commute (both diagonal), so the joint depth-2 trace is structurally degenerate and the A-vs-B distinction is invisible. To distinguish primitive from shuffle, the substrate must induce a non-trivial commutator $[D^2, \gamma_P D] \ne 0$. The natural candidates are (i) the off-diagonal CH operator system (WH1 R3.5, where chirality-flipping E1 breaks simultaneous diagonalisation) and (ii) the Paper 28 QED vertex graph (where $\gamma_P$ is combinatorial vertex parity, not diagonal $(-1)^n$). Either substrate is a 2–4 week sprint comparable to the original NA-1 ambition.

---

## 5. Future — if a follow-on sprint is opened

This sprint did NOT find Reading A vs B; the right follow-on is to **re-run NA-1 on a non-diagonal substrate**, NOT to implement shuffle Hopf algebra enrichment. The latter would be premature.

### 5.1 Sprint-scale follow-on: NA-1 on a non-diagonal substrate (2–4 weeks)

**Substrate candidate 1: off-diagonal CH operator system.** WH1 R3.5 (2026-05-04) introduced the off-diagonal CH where chirality-flipping E1 entries connect $(n, +)$ and $(n+1, -)$ shells. On this substrate $[\gamma_P, D] \ne 0$ because $D$ now has off-diagonal entries between even-$n$ and odd-$n$ shells. The depth-2 trace
$$T(t_1, t_2) := \mathrm{Tr}\big(D^2 e^{-t_1 D^2} \gamma_P D\, e^{-t_2 D^2}\big)$$
will NOT collapse to $\mathrm{Tr}(\gamma_P D^3 e^{-(t_1+t_2)D^2})$ because the kernel $\gamma_P D$ doesn't propagate the heat kernel diagonally. Direct numerical computation on the truncated $n_{\max} \in \{2, 3\}$ off-diagonal CH, with bit-exact `sympy.Rational` matrix arithmetic, is sprint-scale. Existing module: `geovac/full_dirac_operator_system.py` (WH1 R3.5).

**Substrate candidate 2: Paper 28 QED vertex graph.** In Paper 28, $\gamma_P$ is the *combinatorial* vertex parity (even-degree vs odd-degree vertex sets), and the master Mellin engine M3 mechanism uses Hurwitz at quarter-integer shifts. Here $\gamma_P$ does not commute with the kernel because graph-theoretic adjacency mixes the parity sets. Existing module: `geovac/qed_vertex.py`. The joint trace would be computed against the graph adjacency rather than the CH spectrum.

**Decision rule (for the follow-on sprint):** if either substrate gives a symmetric factorisation (J = M2(s_1) M3(s_2), bit-exact under swap), Reading A; if asymmetric (deconcatenation pair, bit-exact under swap+reordering), Reading B.

### 5.2 Multi-month follow-on (only if Reading B wins the substrate-corrected test): shuffle Hopf algebra enrichment

This is the multi-month scope outlined in the strategic synthesis §5 multi-month item 1, contingent on Reading B winning the non-diagonal substrate test. **Do NOT implement at this stage** — the current sprint's verdict (Reading C diagonal collapse) makes the choice of A vs B substrate-dependent rather than substrate-invariant, so the right enrichment may be lighter than full shuffle Hopf if Reading C generalises across substrates.

### 5.3 Paper-edit recommendation

**No corpus edits in this sprint.** A clean §2 entry in CLAUDE.md and a CHANGELOG line are sufficient. Paper 56 is being edited in a sibling worktree (per the sprint prompt constraint), and Paper 55 §4 + §5 already classify the period content correctly — NA-1's pure-Tate finding for $M_3^{\gamma_P}$ corroborates rather than contradicts. A small remark in Paper 18 §III.7 ("on the CH-diagonal substrate, depth-2 Mellins collapse to depth-1 via simultaneous diagonalisation; depth-2 enrichment is operator-substrate dependent") would be appropriate at the next Paper 18 touch, but is not load-bearing here.

---

## 6. Honest scope

### 6.1 What this sprint closed

- **Bit-exact algebraic identity** $J(s_1, s_2) = M_3^{\gamma_P}(s_1 + s_2 - 1)$ on the CH-diagonal substrate, at every integer $s_1 + s_2 \ge 3$, derived termwise and verified at $200$ dps PSLQ across 5 cells.
- **Cross-precision-agreed Hain-Brown negative at depth 2** on this substrate (extends yesterday's depth-1 negative cleanly).
- **Pure-Tate closed form** $M_3^{\gamma_P}(s) = 2^{2s-3}(\beta(2s-1) - \beta(2s-3))$, bit-exact at $s \in \{2, 3, 4\}$, lives in $\bigoplus_k \pi^{2k+1}\mathbb{Q}$ — strictly cleaner than the level-4 cyclotomic $M_3$ of Paper 28, which uses Hurwitz at quarter-integer shift.
- **Structural distinction made explicit** between the diagonal parity-graded "M3" (pure-Tate) and the Paper 28 vertex-parity-discriminant "M3" (level-4 cyclotomic, $\beta(s) - \beta(s-2)$).

### 6.2 What this sprint did NOT close

- **The original A vs B distinction.** The CH-diagonal substrate is *structurally insensitive* to the primitive-vs-shuffle Hopf-coproduct question. Distinguishing A from B requires a non-diagonal substrate where $[D^2, \gamma_P D] \ne 0$.
- **The depth-2 Mellin engine in the Paper 32 §VIII case-exhaustion theorem.** The current theorem is depth-1 (Mellin of $\mathrm{Tr}(D^k e^{-t D^2})$). NA-1's collapse identity shows that on diagonal substrates, depth-2 reduces to depth-1; off-diagonal depth-2 is a genuine extension. A depth-2 generalisation of the case-exhaustion theorem is named as a multi-month follow-on (multi-month item 1 of strategic synthesis §5), contingent on the non-diagonal substrate sprint outcome.

### 6.3 What the sprint flagged for the next person

- **The pure-Tate vs level-4 cyclotomic split inside the "M3" labelling.** Paper 18 §III.7 should be tightened so the M3 label distinguishes (1) diagonal-parity-grading-induced pure-Tate from (2) vertex-restricted Hurwitz-at-quarter-integer level-4 cyclotomic. The case-exhaustion theorem of Paper 32 §VIII implicitly covers both, but the period ring contents differ — they should be separately stated.
- **The depth-2 invisibility of A vs B on the CH-diagonal substrate** is reminiscent of the G3 negative (γ_GV and γ_F independent commuting $\mathbb{Z}_2$'s — the substrate makes the question vacuous). The substrate-corrected NA-1 test on off-diagonal CH or QED vertex graph is the natural next move.

---

## 7. Discipline checks

**Curve-fit audit (`memory/feedback_audit_numerical_claims`):** Applied. Zero free parameters. The claim "J(s_1, s_2) = M_3^{γ_P}(s_1+s_2-1)" is a termwise algebraic identity from the spectrum. PSLQ at 3 precisions × 5 cells with cross-precision agreement filter: 5/5 cells agree across all three precisions on the same single-basis-element relation. The closed form $M_3^{\gamma_P}(s) = 2^{2s-3}(\beta(2s-1) - \beta(2s-3))$ is derived from an algebraic substitution ($u = n + 3/2$) and the standard alternating-Hurwitz-at-half-integer-shift identity, then bit-exact verified at $s \in \{2, 3, 4\}$ at $60$ dps.

**Discrete-for-skeleton (`memory/feedback_discrete_for_skeleton`):** Mostly applicable. The depth-1 sector $M_3^{\gamma_P}(s)$ has a closed-form derivation; the Layer-2 PSLQ identification serves to confirm rather than discover. Layer-1 substrate (CH spectrum, parity grading) used `sympy.Rational` exactly.

**Tag transcendentals (`memory/feedback_tag_transcendentals`):** Every transcendental tagged.
- $\pi^k\mathbb{Q}$ closed forms for $M_3^{\gamma_P}$ — Paper 18 tier: $M_2$-like (pure-Tate, even-weight side of M2) or M1-product (odd-weight extension of M1 via $\beta(\mathrm{odd}) \in \pi^{\mathrm{odd}}\mathbb{Q}$).
- $\beta(\mathrm{odd})$ values — pure-Tate by Euler-number closed form, not modular.
- Hain–Brown $E_4(2i), E_6(2i), \Delta(2i)$, Eichler periods of $\Delta$ — tested NEGATIVE across all precisions.

The natural Paper 18 §III.7 master-Mellin classification of $M_3^{\gamma_P}$ is: this is a **degenerate M2-coupled-to-M1 product**, NOT the cyclotomic-level-4 M3. The M3 label was misleading; the operator structure (diagonal parity grading on integer-shell index) produces pure-Tate. The Paper 18 §III.7 partition should distinguish parity-grading-on-CH (M2 ⊕ M1 product) from vertex-restricted-Hurwitz-at-quarter-integer (level-4 M3) — these are operationally different sub-mechanisms within the same $k = 1$ slot of the case-exhaustion theorem.

**Diagnostic-before-engineering (`memory/feedback_diagnostic_before_engineering`):** This IS a diagnostic sprint. The decision-gate articulation in the strategic synthesis made the protocol clean.

**Algebraic-first (`memory/feedback_algebraic_first`):** The closed-form $M_3^{\gamma_P}(s) = 2^{2s-3}(\beta(2s-1) - \beta(2s-3))$ was derived algebraically (substitution + standard identity), not discovered by PSLQ — PSLQ confirms a pattern that was algebraically derivable from the spectrum.

**Hard prohibitions check clean.** Paper 2 untouched. Paper 56 untouched per sprint prompt. No production `geovac/` modules modified. WH1 PROVEN unaffected.

---

## 8. Files

### Produced this sprint
- `debug/sprint_na1_depth2_mellin_compute.py` — driver (~620 lines, 38 s wall).
- `debug/data/na1_depth2_mellin_results.json` — full results panel (PSLQ × 3 precisions × 5 s_tot, factorisation diff table, cross-precision filter).
- `debug/sprint_na1_depth2_mellin_memo.md` — this memo (~4500 words).

### Read context (per sprint prompt)
- `debug/strategic_synthesis_2026_06_06_memo.md` (Addendum + §2 NA-1 + §6 Recommendation A)
- `debug/sprint_q5p_na1_non_abelian_probe_memo.md` (yesterday's CMM-abelian structural verdict)
- `debug/sprint_hb_pslq_test_memo.md` (today's Hain–Brown depth-1 negative)
- `debug/sprint_hb_adoption_survey_memo.md` (the converged-question framing)
- `papers/group3_foundations/paper_55_periods_of_geovac.tex` §4 + §5 + §6 (M2/M3 closed forms; not modified)
- `papers/group3_foundations/paper_56_tannakian_substrate.tex` §sec:open_na1 (not modified per sprint constraint)
- `geovac/qed_vertex.py` (Paper 28 vertex-parity-Hurwitz M3 sector — different operator structure, not used directly in this sprint but referenced for Reading C distinction)

### References
- Cartier, P. ``A primer of Hopf algebras'' in *Frontiers in Number Theory, Physics, and Geometry II* (Springer, 2007). [CMM theorem]
- Camporesi, R. and Higuchi, A. ``On the eigenfunctions of the Dirac operator on spheres and real hyperbolic spaces.'' J. Geom. Phys. 20 (1996) 1-18. [Dirac S³ spectrum]
- Paper 28 §2 (Theorem 1, closed form for $\zeta_{D^2}(s)$ at integer $s$).
- Paper 28 §4 (Theorem 3, parity discriminant $D_{\mathrm{even}}(s) - D_{\mathrm{odd}}(s) = 2^{s-1}(\beta(s) - \beta(s-2))$ — the level-4 cyclotomic M3 mechanism).
- Paper 32 §VIII (master Mellin engine case-exhaustion theorem at depth 1).
- Paper 18 §III.7 (M1/M2/M3 three-bullet partition).

---

## 9. One-line verdict

**NEITHER Reading A nor Reading B is supported; Reading C (diagonal-collapse) is forced bit-exactly.** The joint depth-2 Mellin $J(s_1, s_2) = M_3^{\gamma_P}(s_1 + s_2 - 1)$ algebraically collapses on the CH-diagonal substrate because $D^2$, $\gamma_P$, $D$, and $e^{-t D^2}$ are simultaneously diagonal; the A-vs-B distinction is structurally invisible on this substrate. Hain–Brown at depth 2: NEGATIVE (no modular content across all three precisions). Substantive sub-finding: the natural diagonal parity-graded $M_3^{\gamma_P}$ sector has pure-Tate closed form $M_3^{\gamma_P}(s) = 2^{2s-3}(\beta(2s-1) - \beta(2s-3))$ in $\bigoplus_k \pi^{2k+1}\mathbb{Q}$ — structurally distinct from Paper 28's level-4 cyclotomic vertex-parity-discriminant M3. **The right follow-on is to repeat NA-1 on a non-diagonal substrate (off-diagonal CH from WH1 R3.5 or the QED vertex graph of Paper 28); the substrate-corrected sprint is 2–4 weeks.** Strategic-synthesis Recommendation B (explicit $U^*_{GV} \hookrightarrow \mathcal{U}_4$ injection) is unaffected and remains the strongest near-term theorem-grade target.
