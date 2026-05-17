# Sprint L2 Paper-Edits Round 1 — Summary Memo

**Sprint:** Paper-edits application round 1 for Sprints L2-B + L2-C + L2-D
**Date:** 2026-05-16
**Verdict:** **CLOSED.** All seven planned edits applied; all four modified papers compile to PDF in a three-pass workflow with no undefined references introduced by these edits. The pre-existing Paper 32 compile-blocker (undefined `\Tcal`/`\sthree`/`\Op`/`\nmaxa`/`\nmaxb` macros in pre-existing §VIII.D content) was resolved via a small surgical fix to the preamble (`\providecommand` definitions), which is documented as a pre-existing-bug fix rather than as new content.

**Companion deliverables:**

- `papers/synthesis/paper_32_spectral_triple.tex` (+~280 lines new content; ~10 lines preamble fix)
- `papers/observations/paper_34_projection_taxonomy.tex` (+~35 lines; +1 row status update)
- `papers/core/paper_18_exchange_constants.tex` (+~25 lines)
- `papers/standalone/paper_42_modular_hamiltonian_four_witness.tex` (+~40 lines)
- `CLAUDE.md` (+3 sprint bullets in §2, ~1 sentence appended to §1.7 WH1 entry)

**Source memos and data consumed:**

- `debug/l2_b_krein_construction_memo.md` (Krein space, Peskin-Schroeder chiral basis, axiom table)
- `debug/l2_c_lorentzian_dirac_memo.md` (Lorentzian Dirac via vdD Prop 4.1, i^t sign, chirality finding)
- `debug/l2_d_connes_axiom_audit_31_memo.md` (BBB Table 1 signs at (4, 6), six-axiom panel, M3 convention-dependence)
- `debug/data/l2_b_krein_construction.json`, `l2_c_lorentzian_dirac.json`, `l2_d_connes_axiom_audit_31.json`

---

## §1. Edits applied (per-paper)

### Edit 1: Paper 32 §IV — Connes axiom audit table extended to (m, n) = (4, 6) Lorentzian columns

**File:** `papers/synthesis/paper_32_spectral_triple.tex` lines ~1080–1180 (around the existing `\subsection{Summary}` block, before `\begin{theorem}[Rigorous identification]`)

**Changes:**
- Extended the caption of `tab:axiom_audit` to clarify that it covers signature $(3, 0)$ and to cross-reference the new Lorentzian-extension table.
- Added `\paragraph{Lorentzian-extension audit at (m, n) = (4, 6) ...}` introducing the L2-B/C/D scope and citing `\cite{bizi_brouder_besnard2018}`.
- Added new `tab:axiom_audit_lorentzian` with 8 rows comparing KO-dim 3 (Riemannian) vs (m, n) = (4, 6) Lorentzian:
  - $J^2$, $JD$, $J\chi$, $J\eta$, $\eta\chi$ pass bit-exact with BBB-predicted signs at $(m, n) = (4, 6)$;
  - $\chi D$ fails (structural finding);
  - Order-zero / order-one finite-resolution residuals match the Riemannian-side audit.
- Added `\paragraph{Structural finding on $\chi D$.}` (~50 lines) explaining the mechanism (chirality-as-$\gamma^5$ + chirality-diagonal $D_{\mathrm{GV}}$), naming the three resolutions R1/R2/R3, and recommending R1 + R2 in parallel for L2-E.

**Diff summary:** +~110 lines net. Honest framing per CLAUDE.md §13.5: the $\{\chi, D\} \neq 0$ result is documented as a structural feature of GeoVac's chirality-diagonal $D_{\mathrm{GV}}$, not a bug; the resolution path mirrors the existing Riemannian R3.5/Paper 42 truthful-vs-offdiag trade.

### Edit 2: Paper 32 §VIII.E — Three new subsections L2-B, L2-C, L2-D

**File:** `papers/synthesis/paper_32_spectral_triple.tex` lines ~3088–3340 (new subsections inserted after the existing `sec:lorentzian_audit_paper32` and before `sec:modular_hamiltonian_l1`)

**Three new subsections added:**

1. **§VIII.E.B — Krein space construction (Sprint L2-B), label `sec:l2b_krein`** (~125 lines). Defines $\mathcal{K}_{n_{\max}, N_t}$ with equation `eq:l2b_krein_space`; states the Peskin-Schroeder chiral basis convention (West-coast metric, $\gamma^5$ diagonal, `FullDiracLabel.chirality` $\leftrightarrow$ $\gamma^5$ eigenvalue); defines $J$ with equation `eq:l2b_J_def`; lists the four Krein axioms with bit-exact residuals; documents the load-bearing Riemannian-limit check (bit-identical at $N_t = 1$ at every $n_{\max} \in \{1, 2, 3\}$); cites the implementation module and tests.

2. **§VIII.E.C — Lorentzian Dirac operator $D_L$ on $S^3 \times \mathbb{R}$ (Sprint L2-C), label `sec:l2c_lorentzian_dirac`** (~110 lines). vdD Prop 4.1 lift with equations `eq:l2c_dirac_decomposition` and `eq:l2c_lorentzian_dirac`; centered FD + Dirichlet zero BC structurally forced for anti-Hermitian $\partial_t$ (equation `eq:l2c_partial_t`); $i^t = +i$ derived from $D_L^\times = D_L$ requirement; bit-exact Krein-self-adjointness + Riemannian-limit recovery + real-spectrum-on-$\mathcal{K}^+$ at 9 panel cells; structural finding on $\{\gamma^5, D_L\}$ flagged and forwarded to L2-D.

3. **§VIII.E.D — Connes axiom audit at $(m, n) = (4, 6)$ (Sprint L2-D), label `sec:l2d_axiom_audit_31`** (~145 lines). BBB Table 1 signs at $(4, 6)$ verified directly from `\cite{bizi_brouder_besnard2018}` v2; 4-spinor charge conjugation $U_4 = i\gamma^2$; lift to Krein space with equation `eq:l2d_J_lorentzian`; six-axiom panel with the four BBB-predicted-sign axioms bit-exact + axiom (v) structural finding + finite-resolution (vi)–(vii); the three resolutions R1/R2/R3 are catalogued; M3 trivialization L0 prediction documented as CONVENTION-DEPENDENT (falsified under Paper 28 $n_{\mathrm{Fock}}$-parity reading, trivially confirmed under chirality-pairing); BBB sign-flip $\{J, \gamma^5\} = 0$ identified as spectral-triple-axiom statement, not vertex-parity-sum statement.

**Bibitem added:**

- `\bibitem{vandungen2016}` — K. van~den~Dungen, "Krein spectral triples and the fermionic action," Math. Phys. Anal. Geom. 19, 4 (2016), arXiv:1505.01939. (BBB 2018 was already present in the bibliography.)

**Diff summary:** +~280 lines net + 1 new bibitem.

### Edit 3: Paper 34 §V.E — Lorentzian transfer audit refinement

**File:** `papers/observations/paper_34_projection_taxonomy.tex` lines ~5722, ~5755 (existing §V.E `sec:lorentzian_transfer_audit`)

**Changes:**

1. **Row 27 (Wick rotation) mechanism column updated** from "THE BRIDGE; Sprint L1 deepens to operator-system level." to "THE BRIDGE; Sprint L2-B/C verified Krein construction at finite cutoff (bit-exact Riemannian-limit recovery)." This is informational; the bucket assignment (W.M.F. / High confidence) is unchanged.

2. **New `\paragraph{Sprint L2-D refinement of the M3 trivialization prediction (2026-05-16)}` appended** (~25 lines) below the existing "Two predicted trivializations at $(3, 1)$" paragraph. Documents Sprint L2-D's verdict: M3 trivialization is convention-dependent; falsified under Paper 28's $n_{\mathrm{Fock}}$-parity convention (the convention that accesses Catalan $G$ and $\beta(4)$); trivially confirmed under chirality-pairing convention but tautologically. Structural refinement: BBB sign-flip $\{J, \gamma^5\} = 0$ is a spectral-triple-axiom statement, not a vertex-parity-sum statement. Cross-references Paper 32 §VIII.D and Paper 18 §III.7.

3. **Cross-reference fix**: replaced `\ref{sec:l2d_axiom_audit_31}` (which would not resolve cross-paper) with the literal "§VIII.D" string.

**Diff summary:** +~35 lines net + 1 row status update.

### Edit 4: Paper 18 §III.7 — Master Mellin engine sectional-convention sharpening

**File:** `papers/core/paper_18_exchange_constants.tex` lines ~957–981 (inserted between "Mechanism-as-domain sharpening" paragraph and "Discrete $c_1$ verification" paragraph)

**Changes:** Added `\paragraph{M3 sectional convention refinement (Sprint L2-D, May 2026)}` (~25 lines). Documents that M3's mechanism is sectional in the chosen parity convention, not in spacetime signature: the L0 prediction holds trivially under chirality-pairing but falsifies under Paper 28's $n_{\mathrm{Fock}}$-parity reading. Names the closed-form identity $D_{\mathrm{even}}(s) - D_{\mathrm{odd}}(s) = 2^{s-1}(\beta(s) - \beta(s-2))$ as the M3 access point that survives. Clarifies that the BBB sign-flip $\{J, \gamma^5\} = 0$ is a spectral-triple-axiom statement, not a vertex-parity-sum statement, and does not change the master Mellin engine domain partition itself.

Cross-references Paper 32 §VIII.D and Paper 28 Theorem 3 (via existing `\cite{loutey_paper32}` and `\cite{loutey_paper28}` bibitems).

**Diff summary:** +~25 lines.

### Edit 5: Paper 42 §10 — Open Question O1 progress update

**File:** `papers/standalone/paper_42_modular_hamiltonian_four_witness.tex` lines ~1568–1620 (extending the existing `(O1) Lorentzian extension at signature (3, 1)` block)

**Changes:** Added an `\emph{Progress update (Sprints L2-B/C/D, 2026-05-16)}` block (~40 lines). Documents that the Hilbert-space and Lorentzian-Dirac levels of the $(3, 1)$ extension have landed:
- L2-B Krein space (`geovac/krein_space_construction.py`) with all axioms bit-exact and load-bearing Riemannian-limit recovery passing;
- L2-C Lorentzian Dirac (`geovac/lorentzian_dirac.py`) via vdD Prop 4.1 with Krein-self-adjointness + Riemannian-limit recovery bit-exact;
- L2-D Connes axiom audit (`geovac/connes_axiom_audit_31.py`) with four BBB-predicted axioms bit-exact and the fifth ($\{\gamma^5, D\} = 0$) flagged as the load-bearing scope finding requiring R1+R2 resolution.

Notes that WH1 PROVEN (Paper 38) is NOT re-opened — the Krein construction is structurally additive (a new (3, 1) layer on top of the Riemannian foundation), not a re-test of the Riemannian foundation. O1 status remains OPEN; the operator-system-level closure is the named L2-E target.

Cross-references `\cite{paper32}` (Paper 32 §IV / §VIII.D) and `\cite{paper38}` (already in the bibliography).

**Diff summary:** +~40 lines.

### Edit 6: CLAUDE.md §2 — Three sprint bullets for L2-B/C/D

**File:** `CLAUDE.md` lines ~345 (inserted above the existing Sprint L1 bullet)

**Changes:** Three new entries added in chronological order (L2-B, L2-C, L2-D, all dated 2026-05-16), each ~250–450 words, matching the style of existing §2 entries. Each bullet includes:
- Verdict (CLOSED / CLOSED-WITH-STRUCTURAL-FINDING);
- Headline numerical result (residual at bit-exact, panel size);
- Key structural finding;
- Files (memo, module, tests, data);
- Cross-references to what was unblocked / what's next.

The L2-D bullet documents the BBB Table 1 sign verification, the four BBB-predicted-sign axioms passing bit-exact, the structural finding on $\{\chi, D_L\}$ failing on truthful $D_{\mathrm{GV}}$, the three resolutions R1/R2/R3 with recommended R1+R2 path, and the M3 trivialization convention-dependence finding (FALSIFIED under $n_{\mathrm{Fock}}$-parity, trivially-confirmed under chirality-pairing).

**Diff summary:** +3 sprint bullets, ~150-300 words each.

### Edit 7: CLAUDE.md §1.7 WH1 entry — L2-B/C closure sentence

**File:** `CLAUDE.md` line 64 (appended to the end of the existing WH1 entry)

**Changes:** Added one sentence at the end of the WH1 PROVEN entry noting that Sprint L2-B (Krein space) and L2-C (Lorentzian Dirac) both passed bit-identical Riemannian-limit checks at finite $n_{\max}$, confirming the Camporesi-Higuchi spinor structure is structurally compatible with the Cl(3, 1) Krein-space embedding. Explicitly notes that WH1 PROVEN is not re-opened — the Lorentzian extension is structurally additive (new (3, 1) layer on top of the proven Riemannian foundation), not a re-test.

**Diff summary:** +1 sentence (~70 words).

---

## §2. Compilation status

All four papers were compiled with `pdflatex -interaction=nonstopmode -halt-on-error` in a three-pass workflow.

### Paper 32 (`paper_32_spectral_triple.tex`): **SUCCESS** (after pre-existing macro fix)

- **48 pages**, 637,866 bytes.
- **Pre-existing compile blocker resolved.** The file had a pre-existing pre-this-sprint usage of undefined macros `\Tcal`, `\Op`, `\sthree`, `\nmaxa`, `\nmaxb` in §VIII.D (frontier-of-field paragraphs, lines ~1357–1377). These macros were referenced but never defined in the preamble. I confirmed via `git stash` that the file did not compile clean before my edits applied. I added a minimal `\providecommand` block to the preamble:
  ```latex
  \providecommand{\Tcal}{\mathcal{T}}
  \providecommand{\Op}{\mathcal{O}}
  \providecommand{\sthree}{S^3}
  \providecommand{\nmaxa}{n_{\max, a}}
  \providecommand{\nmaxb}{n_{\max, b}}
  ```
  Definitions are conservative and consistent with the usage contexts (calligraphic / subscripted-mathematical-italic). The preamble comment annotates this as restoring missing definitions from pre-existing content.
- **Remaining warnings**: pre-existing — one undefined citation (`suijlekom_book2015` typo for `vansuijlekom_book2015`); two undefined references (`sec:tensor_product_verification`, `thm:eta_trivialization`). None are introduced by my edits.
- **None of my new content uses any of the previously undefined macros** (verified by grep on my insertion range).

### Paper 34 (`paper_34_projection_taxonomy.tex`): **SUCCESS**

- **91 pages**, 1,108,451 bytes.
- **Remaining warnings**: pre-existing — three undefined references (`sec:matches`, `sec:curvature_coefficients`, `tab:catalog_off`). All flagged in prior commits; not introduced by my edits.
- **First-attempt issue fixed in-flight**: my initial edit used `\ref{sec:l2d_axiom_audit_31}` which does not resolve cross-paper; replaced with the literal "§VIII.D" string.

### Paper 18 (`paper_18_exchange_constants.tex`): **SUCCESS**

- **(REVTeX style)** 658,710 bytes. Page count not extracted from this log but the PDF generated cleanly.
- **Remaining warnings**: standard hyperref Unicode warnings (pre-existing, unrelated to my edits).
- BibTeX run between passes 1 and 2 to resolve natbib citations; this is the standard REVTeX workflow.

### Paper 42 (`paper_42_modular_hamiltonian_four_witness.tex`): **SUCCESS**

- **23 pages**, 594,985 bytes.
- **Zero undefined references after three passes.** All citations resolve.
- All warnings are standard layout warnings (`\hbox`, font shape) — no semantic warnings.

---

## §3. Honest-discipline preservation

Per CLAUDE.md §13.5, the following discipline was preserved through all edits:

- **WH1 PROVEN is NOT re-opened.** Every L2-B/C/D edit explicitly notes that the Lorentzian extension is structurally additive (a new $(3, 1)$ layer on top of the proven Riemannian foundation), not a re-test of WH1's Riemannian foundation. The WH1 §1.7 entry's "Status: PROVEN" classification is untouched.

- **The $\{\gamma^5, D_L\} \neq 0$ finding is documented as a structural feature of GeoVac, not as a bug.** Every appearance (Paper 32 §IV, §VIII.E.D; CLAUDE.md L2-C and L2-D bullets) frames it as a consequence of two locked Riemannian-side conventions (chirality-as-$\gamma^5$ from L2-B, chirality-diagonal $D_{\mathrm{GV}}$ from Paper 32 §III) plus BBB's universal axiom — a mutually-inconsistent triple where the resolution is to pick R1+R2 in parallel, mirroring the existing Riemannian R3.5/Paper 42 trade.

- **M3 trivialization L0 prediction is documented as convention-dependent.** Honest reporting: the prediction is falsified under one natural convention ($n_{\mathrm{Fock}}$-parity, the one that accesses Catalan $G$) and trivially confirmed under another (chirality-pairing). The structural refinement is named without hiding the negative.

- **No fitted parameters added; no new physics introduced; no "conjectural" labels removed.** The combination rule $K = \pi(B + F - \Delta)$ in Paper 2 is not touched by any of these edits (no cross-references to it).

- **The pre-existing Paper 32 compile-blocker fix is documented in-place** with a preamble comment noting the date and sprint origin of the fix, and labeled as "restoring missing definitions" rather than as new content. The fix is conservative (uses `\providecommand` so any later top-level macro definitions take precedence).

- **Bibliography integrity preserved.** The new `vandungen2016` bibitem in Paper 32 follows the existing entry format; `\cite{bizi_brouder_besnard2018}` was already in Papers 32, 34, and 42 and not re-added; `\cite{loutey_paper32}` / `\cite{loutey_paper28}` already in Paper 18 bibliography.

---

## §4. Cross-paper consistency check

The new content was checked for consistency across the four papers:

- **Cl(3, 1) sign conventions:** consistent across Paper 32 §IV table, Paper 32 §VIII.E.B/C/D, Paper 42 §10 O1 update. All four BBB-predicted-sign axioms at $(4, 6)$ are stated identically: $\varepsilon = +1$, $\varepsilon'' = -1$, $\kappa = -1$, $\kappa'' = +1$, giving $J^2 = +I$, $\{J, \chi\} = 0$, $\{J, \eta\} = 0$, $\{\eta, \chi\} = 0$.

- **Riemannian-limit bit-exactness language:** consistent across all locations. Phrased as "bit-identical at $N_t = 1$" or "Frobenius residual = 0.0 in float64" depending on context. The mechanism (Kronecker product of real $\{0, 1\}$-valued permutation with identity) is explained once in Paper 32 §VIII.E.B and referenced from the other locations.

- **Chirality finding phrasing:** the structural finding $\{\gamma^5, D_L\} \neq 0$ has a single canonical phrasing — "GeoVac's chirality-diagonal $D_{\mathrm{GV}}$ commutes with $\gamma^5$ in the Peskin-Schroeder chiral basis" — used in Paper 32 §IV, Paper 32 §VIII.E.D, Paper 42 §10 O1, CLAUDE.md L2-C/L2-D bullets. The three resolutions R1/R2/R3 are named identically everywhere they appear.

- **M3 trivialization convention dependence:** documented identically in Paper 32 §VIII.E.D, Paper 34 §V.E, Paper 18 §III.7. Convention names are consistent: "$n_{\mathrm{Fock}}$-parity convention" / "Paper 28 §QED-vertex convention" vs "chirality-pairing convention." The L0 prediction is consistently labeled as falsified under the first and trivially-confirmed under the second.

- **Cross-references between papers:** Paper 18 cites Paper 32 §VIII.D for L2-D; Paper 34 cites Paper 32 §VIII.D for L2-D; Paper 42 §10 cites Paper 32 §IV / §VIII.D and Paper 38. All cross-references use existing bibitems.

---

## §5. Items NOT changed

Per CLAUDE.md §13.5 access control and §13.8 paper edit policy:

- **No production code modified.** Only `.tex` files and `CLAUDE.md`.
- **Paper 2 (`papers/observations/paper_2_alpha.tex`) not touched.** No cross-references to Paper 2 introduced by any L2-B/C/D edit.
- **Paper 38 not modified.** Cross-referenced from Paper 42 §10 update only, using the existing `\cite{paper38}` bibitem.
- **CLAUDE.md §1.5 / §1.6 / §4 / §5 / §13 / §14 not touched.** All edits sit within PM-editable sections per §13.5.
- **Paper 32 §IV `tab:axiom_audit` Riemannian-side rows not modified.** Only the caption was extended and a new table `tab:axiom_audit_lorentzian` added.

---

## §6. Forward links

The edits applied here document the L2-B/C/D landing and set up Sprint L2-E (the Krein-level Paper 42 redo). Sprint L2-E should:

- Build the Lorentzian-side modular Hamiltonian $K_L = -\log \Delta_L$ on the Krein-positive cone $\mathcal{K}^+$.
- Test the R1+R2 trade in parallel (truthful $D_{\mathrm{GV}}$ for axioms requiring $JD = +DJ$, offdiag CH for axioms requiring $\{\chi, D\} = 0$).
- Verify the Lorentzian analog of Paper 42's $\sigma_{2\pi}(O) = O$ at the operator-system level on the Krein space.
- Document whether the M3 trivialization actually changes any measurable spectral output (currently the L2-D analysis shows it does not, under the only convention that accesses Paper 28's M3 content).

Sprint L2-G synthesis (if/when it happens) can compose the L2-A through L2-E memos into a standalone Lorentzian-extension paper analogous to Paper 42's role for the Riemannian closure.

---

End of memo.
