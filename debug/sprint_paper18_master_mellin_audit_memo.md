# Sprint Paper 18 §III.7 master Mellin engine audit (NA-1 operator-vs-slot refinement)

**Date:** 2026-06-06
**Sprint:** Paper 18 §III.7 audit and refinement to incorporate NA-1 finding (operator-vs-slot distinction).
**Verdict:** **POSITIVE — descriptive refinement applied.** Paper 18 §III.7 master Mellin engine description now reflects the operator-vs-slot distinction surfaced by Sprint NA-1 (Reading C / `rem:gammaP_vs_vertex_restricted` of Paper 55). Three-pass clean compile, 27 pages (was 26). No new undefined references introduced.

---

## 1. The audit question

The NA-1 finding (2026-06-06; `debug/sprint_na1_depth2_mellin_memo.md`) established that the **same parity-grading operator $\gamma_P$** on the natural Camporesi–Higuchi diagonal substrate produces two structurally distinct period-ring outputs:

- At Mellin slot exponent $r = 2s - 1$ (odd integer): $M_3^{\gamma_P}(s) = 2^{2s-3}(\beta(2s-1) - \beta(2s-3))$, **pure-Tate** in $\bigoplus_k \pi^{2k+1}\cdot\mathbb{Q}$.
- At integer $s$ via $D_{\mathrm{even}}(s) - D_{\mathrm{odd}}(s) = 2^{s-1}(\beta(s) - \beta(s-2))$: **level-4 cyclotomic** in $\mathcal{MT}(\mathbb{Z}[i, 1/2])$ at level 4.

Both probes act on the same $\gamma_P$; the period content (pure-Tate vs level-4 cyclotomic) is selected by the **Mellin slot exponent**, not by a different operator. Paper 55 captured this in §subsec:m3_diagonal_collapse + §subsec:m3_galois_descent with new Theorems and the explicit Remark `rem:gammaP_vs_vertex_restricted`.

The audit question: does Paper 18 §III.7's authoritative description of the master Mellin engine reflect this distinction? **Answer at audit start: no.** Paper 18 §III.7 stated the M1/M2/M3 partition via $k \in \{0, 1, 2\}$ as "operator order selects sub-mechanism" but did not separate (i) the sub-mechanism index $k$ that fixes which member of the operator family $\mathrm{Tr}(D^k \cdot e^{-tD^2})$ is being read, from (ii) the slot exponent $s$ that selects which periods within the sub-mechanism's period ring are accessed.

## 2. What was refined in Paper 18 §III.7

One new paragraph inserted between the "M3 sectional convention refinement (Sprint L2-D, May 2026)" paragraph and the "Compact / non-compact partition of the $4/\pi$ rate (Sprint M-Z, May 2026)" paragraph. Title: **"Operator-vs.-slot refinement of the M3 sub-mechanism (Sprint NA-1, June 2026)."** Content:

1. **Statement of the distinction.** The Mellin-slot index $k \in \{0, 1, 2\}$ classifies the sub-mechanism (M1 / M3 / M2). It does NOT by itself fix the period ring of the output. The Mellin slot exponent $s$ — the argument at which the moment is evaluated — selects which periods inside the sub-mechanism's natural period ring are accessed.

2. **The M3 instance from NA-1.** Definition of $M_3^{\gamma_P}(r) := \sum_n (-1)^n g_n |\lambda_n|^{-r}$ as the parity-graded Mellin function of the Camporesi–Higuchi Dirac spectrum, with explicit pure-Tate slot ($r = 2s - 1$ odd integer; closed form $2^{2s-3}(\beta(2s-1) - \beta(2s-3))$) and level-4 cyclotomic slot (integer $s$ on $D_{\mathrm{even}} - D_{\mathrm{odd}}$; Paper 28 Theorem 3). Both probes act on the same $\gamma_P$; period content selected by slot exponent.

3. **Scope clarification.** This refinement does NOT change (a) the case-exhaustion statement (Paper 32 §VIII), (b) the M3 sectional convention refinement above it, or (c) the period-ring stratification $M_3 \subset \mathcal{MT}(\mathbb{Z}[i, 1/2])$ at level $\le 4$. It sharpens the description of the M3 row: the sub-mechanism index $k = 1$ commits to the $\eta$-invariant Mellin family $\mathcal{M}[\mathrm{Tr}(D \cdot e^{-tD^2})]$, but the period-ring location within $M_3$ at each evaluation is determined by the slot exponent.

4. **Cross-references** to Paper 55 §subsec:m3_diagonal_collapse and §subsec:m3_galois_descent for the explicit demonstration including the NA-1 diagonal-collapse identity $J(s_1, s_2) = M_3^{\gamma_P}(s_1 + s_2 - 1)$.

5. **Extension to M2 and M1.** Identical operator-vs-slot reading is recorded for M2 (Seeley–DeWitt at $k = 2$): operator $D^2$ is fixed, but the Mellin slot exponent $s$ selects which Seeley–DeWitt coefficient $a_k(D^2)$ is read off via $\mathrm{Res}_{s = (d-2k)/2}\,\zeta_{D^2}(s)$. The pure-Tate M2 closure $M_2 \subset \bigoplus_k \pi^{2k}\cdot\mathbb{Q}$ describes the period ring; the slot exponent selects the graded piece. **M1 is degenerate** in this sense: the trivial heat kernel collapses to the Hopf-base measure $\mathrm{Vol}(S^d)$ and the period ring $\mathbb{Q}[\pi, \pi^{-1}]$ has no internal Mellin-slot stratification to select.

## 3. What was preserved

The M1/M2/M3 partition statement itself was preserved verbatim — operator order $k \in \{0, 1, 2\}$ remains the classifying index for sub-mechanism. The audit was a **descriptive refinement**, not a structural change to the partition.

- The master Mellin engine equation $\pi$-source = $\mathcal{M}[\mathrm{Tr}(D^k \cdot e^{-tD^2})](s)$ with $k \in \{0, 1, 2\}$ is unchanged.
- The case-exhaustion theorem (Paper 32 §VIII) is unchanged.
- The period-ring stratification (eq:m_engine_period_stratification in Paper 18; eq:ring_nesting in Paper 55) is unchanged.
- The Mellin domain partition (M1 ↔ propinquity rates; M2 ↔ heat-kernel/spectral-action convergence; M3 ↔ vertex-restricted parity-character spectral sums) is unchanged.
- The M3 sectional convention refinement paragraph (Sprint L2-D) is unchanged.
- The M1 period-ring closure paragraph (Sprint M1 Pure-Tate) is unchanged.

## 4. Bibitem and reference changes

- Added new bibitem `loutey_paper55` ("Periods of GeoVac: Cyclotomic Mixed-Tate Classification of the Master Mellin Engine") to Paper 18 §bibliography, immediately after `loutey_paper31`.
- All existing references (`loutey_paper28`, `loutey_paper32`, `glanois2015`, `deligne2010`, `fathizadeh_marcolli2016`) already present.
- No new undefined `\ref` or `\cite` introduced. The three Paper 55 internal labels (`thm:m3_gammaP_closed_form`, `prop:vertex_parity_descent`, `thm:na1_diagonal_collapse`) are referenced descriptively (e.g., "Paper 55 §M3-diagonal-collapse Theorem on $M_3^{\gamma_P}$ closed form") rather than via `\ref`, since `\ref` to another document's label does not resolve.

## 5. Equation verification (CLAUDE.md §13.4a)

All closed forms cited in the new Paper 18 §III.7 paragraph cross-reference to existing Paper 55 Theorems with their own verification:

| Equation in new Paper 18 paragraph | Paper 55 source | Verification status |
|:-----------------------------------|:----------------|:-------------------|
| $M_3^{\gamma_P}(s) = 2^{2s-3}(\beta(2s-1) - \beta(2s-3))$ | Theorem `thm:m3_gammaP_closed_form` (lines 1107–1130 of paper_55) | Bit-exact PSLQ at $s \in \{2, 3, 4, 5, 6\}$ at 50, 100, 200 dps in `debug/sprint_na1_depth2_mellin_compute.py`; algebraic proof in proof sketch |
| $D_{\mathrm{even}}(s) - D_{\mathrm{odd}}(s) = 2^{s-1}(\beta(s) - \beta(s-2))$ | Paper 28 Theorem 3 | Already verified in Paper 28 (PSLQ at 80 digits per CLAUDE.md §6 Paper 28 description) |
| $J(s_1, s_2) = M_3^{\gamma_P}(s_1 + s_2 - 1)$ | Theorem `thm:na1_diagonal_collapse` (lines 1056–1097 of paper_55) | Bit-exact at $s_{\mathrm{tot}} \in \{3, 4, 5, 6, 7\}$ at 200 dps; data in `debug/data/na1_depth2_mellin_results.json` |
| $\bigoplus_k \pi^{2k+1}\cdot\mathbb{Q}$ pure-Tate placement of $M_3^{\gamma_P}$ at integer $s$ | Euler-number closed form $\beta(2k+1) \in \pi^{2k+1}\cdot\mathbb{Q}$ | Classical; in Glanois 2015 |
| $\mathcal{MT}(\mathbb{Z}[i, 1/2])$ at level 4 placement of $D_{\mathrm{even}} - D_{\mathrm{odd}}$ | Paper 55 Proposition `prop:vertex_parity_descent` | Deligne 2010 + Glanois 2015 + Paper 28 closed forms |

No new equations introduced that require new tests. The refinement is descriptive: every quantitative claim is cross-referenced to a Paper 55 or Paper 28 equation that already has its own verification.

## 6. Three-pass compile

Baseline (pre-edit): 26 pages, three pre-existing undefined refs (`paper7`, `thm:T9`, `sec:tensor_product_verification`) that predate this audit.

Post-edit: **27 pages**, same three pre-existing undefined refs, no new ones. Three-pass clean (pdflatex exit 0 on pass 1, 2, 3). Confirmed by `git stash` baseline comparison: the pre-existing warnings appear at identical line numbers in both baseline and post-edit logs, modulo the one-page line-number shift from the new paragraph insertion.

## 7. Honest scope

This is a **descriptive refinement** of Paper 18 §III.7, not a new structural finding. The structural finding (operator-vs-slot distinction) was captured in Paper 55 §subsec:m3_diagonal_collapse + §subsec:m3_galois_descent in today's NA-1 sprint. Paper 18 §III.7's role as the authoritative description of the master Mellin engine made it the natural place for the distinction to be reflected at the engine-description level.

The audit did NOT uncover any structural inconsistency in the master Mellin engine framework itself. No new Theorem in Paper 18 was needed. No other paper required edits.

## 8. Audit gate verdict

POSITIVE per the decision gate: Paper 18 §III.7 description updated to reflect operator-vs-slot distinction; consistent with Paper 55 §subsec:m3_diagonal_collapse + §subsec:m3_galois_descent; Paper 18 compiles three-pass clean.

## 9. Named open question

One small open item flagged in the new paragraph: the operator-vs-slot reading was demonstrated explicitly only for M3 (closed-form pure-Tate vs level-4 cyclotomic). The text states that an "identical" reading is available for M2 (which Seeley–DeWitt coefficient is read off by the residue at $s = (d-2k)/2$), but this M2 instance of the operator-vs-slot reading is informal rather than supported by a sprint-grade closed-form analog of NA-1. A natural follow-on would be to demonstrate an analogous "diagonal-collapse identity" for M2 (e.g., joint Mellin of $D^2 \cdot e^{-t_1 D^2}$ × $D^2 \cdot e^{-t_2 D^2}$ on the CH diagonal substrate) and identify which $\pi^{2k}$ slot is selected at each $(s_1, s_2)$. Not load-bearing for the audit; flagged here in case a future sprint wants to test it.

## 10. Files modified / created

**Modified:**
- `papers/group3_foundations/paper_18_exchange_constants.tex` — new paragraph (≈65 lines of LaTeX) inserted between the M3 sectional convention paragraph and the Compact / non-compact $4/\pi$ partition paragraph in §III.7. One new bibitem `loutey_paper55` added.

**Created (this memo):**
- `debug/sprint_paper18_master_mellin_audit_memo.md` (this file)

**Not modified:** Paper 55, Paper 56, Paper 28, Paper 32, any other paper. CHANGELOG.md and CLAUDE.md unchanged per sprint scope.

## 11. Cross-references

- Sprint NA-1 source memo: `debug/sprint_na1_depth2_mellin_memo.md` (Reading C / diagonal-collapse identity)
- Paper 55 §subsec:m3_diagonal_collapse (Theorem `thm:na1_diagonal_collapse`, Theorem `thm:m3_gammaP_closed_form`, Remark `rem:gammaP_vs_vertex_restricted`)
- Paper 55 §subsec:m3_galois_descent (Proposition `prop:vertex_parity_descent`)
- Paper 28 Theorem 3 ($\chi_{-4}$ identity for $D_{\mathrm{even}} - D_{\mathrm{odd}}$)
- Paper 32 §VIII case-exhaustion theorem (the master Mellin engine sits inside this theorem)
- Paper 18 §III.7 master-mechanism reading (Track TS-E3, May 2026)
- Paper 18 §III.7 mechanism-as-domain sharpening (Sprint MR-A/B/C, May 2026)
- Paper 18 §III.7 M3 sectional convention refinement (Sprint L2-D, May 2026)
- Paper 18 §III.7 M1 period-ring closure (Sprint M1 Pure-Tate, June 2026)
- `memory/mellin_taxonomy_engine.md`
- `memory/mellin_engine_domain_partition.md`

---

**Audit summary in one sentence.** Paper 18 §III.7 now records that the sub-mechanism index $k \in \{0, 1, 2\}$ in the master Mellin engine fixes which operator family $\mathrm{Tr}(D^k \cdot e^{-tD^2})$ is read, while the Mellin slot exponent $s$ selects which periods within that sub-mechanism's period ring are accessed — with the M3 instance worked out explicitly (pure-Tate at odd-integer slot exponent vs level-4 cyclotomic at integer-$s$ slot exponent, both from the same parity grading $\gamma_P$, per Paper 55 §subsec:m3_diagonal_collapse and §subsec:m3_galois_descent).
