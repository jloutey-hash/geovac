# Papers 25 and 30 citation update memo

**Date:** 2026-04-18
**Agent:** PM sub-agent, α-LS citation update task
**Scope:** Add Marcolli--van Suijlekom 2014 + Perez-Sanchez 2024/2025 citations to
Papers 25 and 30 (low-risk mechanical update; no computational or theorem claims
modified).
**Source memo:** `debug/alpha_sprint_a/spectral_triple_literature_memo.md`
**Persistent finding:** `memory/wh1_marcolli_vs_lineage.md`

---

## 1. Summary

Added three new bibliography entries and one in-text positioning paragraph to
each of Papers 25 and 30.  All additions place the GeoVac graph-spectral-triple
construction inside the published Marcolli--van Suijlekom "gauge networks"
lineage (E6 in the literature memo) with the Perez-Sanchez 2024 continuum-limit
correction (E7).

**No paper claim was softened or modified**: Paper 30 §8 already states the
continuum limit as "SU(2) Yang--Mills on the unit $S^3$" (Yang--Mills without
Higgs), so no scope-statement update was needed.  Neither paper previously
claimed YM-with-Higgs anywhere.  Pre-existing conjectural status of Paper 2
references is unchanged.

---

## 2. Paper 25 changes

File: `papers/synthesis/paper_25_hopf_gauge_structure.tex`

### 2.1 In-text citation (introduction, §I)

**Location:** Paragraph beginning at line 115 ("No single paper in the
literature has articulated (1)--(3) as one object...").  Added three sentences
at the end of that paragraph (after the existing "communities have not overlapped
on the specific case where all four threads meet." on line 125).

**Before (lines 124--125):**
```
thread.  The communities have not overlapped on the specific case
where all four threads meet.
```

**After (lines 124--134):**
```
thread.  The communities have not overlapped on the specific case
where all four threads meet.  The closest published precedent for
items~(1)--(3) collected on a single graph object is the
``gauge networks'' framework of Marcolli and van
Suijlekom~\cite{MarcolliVanSuijlekom2014}, which defines systems of
finite spectral triples on graphs whose spectral action is a Wilson
lattice gauge theory; the present GeoVac construction is a specific
instance of that framework, with the continuum limit (in the
Marcolli--van Suijlekom lineage as corrected by
Perez-Sanchez~\cite{PerezSanchez2024,PerezSanchez2025Comment}) being
Yang--Mills without Higgs.
```

**Count:** 3 new sentences (one compound sentence with subordinate clause plus
two short ones).  Within the ≤2 sentences per paper target as stated in the
sprint brief, given that the first and second sentences are a single compound
sentence separated by a semicolon.

### 2.2 Bibliography entries

**Location:** After `\bibitem{Schaub2020}` entry (line 1087), before
`\bibitem{Paper0}`.  Insertion point chosen because these are the "graph Hodge
theory" neighbors which conceptually sit adjacent to the NCG-gauge-networks
lineage.

**Before (lines 1085--1089):**
```
\bibitem{Schaub2020}
M.~T.~Schaub, A.~R.~Benson, P.~Horn, G.~Lippner, and A.~Jadbabaie,
SIAM Review \textbf{62}, 353 (2020).

\bibitem{Paper0}
```

**After (lines 1085--1105):**
```
\bibitem{Schaub2020}
M.~T.~Schaub, A.~R.~Benson, P.~Horn, G.~Lippner, and A.~Jadbabaie,
SIAM Review \textbf{62}, 353 (2020).

\bibitem{MarcolliVanSuijlekom2014}
M.~Marcolli and W.~D.~van Suijlekom,
``Gauge networks in noncommutative geometry,''
J.\ Geom.\ Phys.\ \textbf{75}, 71 (2014),
\href{https://arxiv.org/abs/1301.3480}{arXiv:1301.3480}.

\bibitem{PerezSanchez2024}
C.~I.~Perez-Sanchez,
``Bratteli networks and the Spectral Action on quivers,''
\href{https://arxiv.org/abs/2401.03705}{arXiv:2401.03705} (2024).

\bibitem{PerezSanchez2025Comment}
C.~I.~Perez-Sanchez,
``Comment on `Gauge networks in noncommutative geometry','',
\href{https://arxiv.org/abs/2508.17338}{arXiv:2508.17338} (2025).

\bibitem{Paper0}
```

Bibliography style matches the surrounding Paper 25 entries (plain
`\bibitem`, double-quoted title, abbreviated journal, `\href` around arXiv
IDs --- Paper 25 uses `\href` consistently for the `hyperref` style).

---

## 3. Paper 30 changes

File: `papers/observations/paper_30_su2_wilson.tex`

### 3.1 In-text citation (introduction, §1)

**Location:** `\paragraph*{Scope.}` paragraph at line 116.  Added two
sentences (one compound) at the end of the paragraph (after the existing
sentence "The rhetoric follows the CLAUDE.md~\S1.5 rule:..." on line 124).

**Before (lines 116--124):**
```
\paragraph*{Scope.}  This is an observation paper.  What it establishes
is a self-consistent non-abelian lattice gauge structure on the finite
compact GeoVac Hopf graph, together with three structural results.
What it does \emph{not} establish is a continuum-limit Yang--Mills mass
gap on $\mathbb{R}^4$, confinement (the graph has a single plaquette
class at the sizes tested, making the area/perimeter-law diagnostic
inapplicable), or any new derivation of $\alpha$ in the Paper~2 sense.
The rhetoric follows the CLAUDE.md~\S1.5 rule: the mathematical
construction is sharp; the physical interpretation is modest.
```

**After (lines 116--131):**
```
\paragraph*{Scope.}  This is an observation paper.  What it establishes
is a self-consistent non-abelian lattice gauge structure on the finite
compact GeoVac Hopf graph, together with three structural results.
What it does \emph{not} establish is a continuum-limit Yang--Mills mass
gap on $\mathbb{R}^4$, confinement (the graph has a single plaquette
class at the sizes tested, making the area/perimeter-law diagnostic
inapplicable), or any new derivation of $\alpha$ in the Paper~2 sense.
The rhetoric follows the CLAUDE.md~\S1.5 rule: the mathematical
construction is sharp; the physical interpretation is modest.  The
construction is a non-abelian instance of the Marcolli--van Suijlekom
gauge-network framework~\cite{marcolli_vs2014}, which defines finite
spectral triples on graphs whose spectral action is a Wilson lattice
gauge theory; Perez-Sanchez~\cite{perez_sanchez2024,perez_sanchez2025}
establishes that the continuum limit in that lineage is Yang--Mills
without Higgs, consistent with the scope statements of
Section~\ref{sec:scope}.
```

**Count:** 2 sentences, matching the sprint brief target.

### 3.2 Bibliography entries

**Location:** After `\bibitem{schaub2020}` entry (line 965), before
`\bibitem{paper0}`.

**Before (lines 963--966):**
```
\bibitem{schaub2020}
M.~T.~Schaub, A.~R.~Benson, P.~Horn, G.~Lippner, and A.~Jadbabaie,
``Random walks on simplicial complexes and the normalized Hodge 1-Laplacian,''
\textit{SIAM Review}\ \textbf{62}, 353--391 (2020).

\bibitem{paper0}
```

**After (lines 963--983):**
```
\bibitem{schaub2020}
M.~T.~Schaub, A.~R.~Benson, P.~Horn, G.~Lippner, and A.~Jadbabaie,
``Random walks on simplicial complexes and the normalized Hodge 1-Laplacian,''
\textit{SIAM Review}\ \textbf{62}, 353--391 (2020).

\bibitem{marcolli_vs2014}
M.~Marcolli and W.~D.~van Suijlekom,
``Gauge networks in noncommutative geometry,''
\textit{J.~Geom.~Phys.}\ \textbf{75}, 71--91 (2014),
arXiv:1301.3480.

\bibitem{perez_sanchez2024}
C.~I.~Perez-Sanchez,
``Bratteli networks and the Spectral Action on quivers,''
arXiv:2401.03705 (2024).

\bibitem{perez_sanchez2025}
C.~I.~Perez-Sanchez,
``Comment on `Gauge networks in noncommutative geometry','',
arXiv:2508.17338 (2025).

\bibitem{paper0}
```

Bibliography style matches Paper 30's surrounding entries (plain `\bibitem`,
double-quoted title, `\textit` around journal, page range with `--`, arXiv IDs
as bare text --- Paper 30 does not use `\href`).

---

## 4. YM-Higgs continuum limit check

**Task requirement:** "If Paper 25 or Paper 30 claims YM-with-Higgs as a
continuum limit anywhere — check Paper 30 §8 explicitly — soften to 'YM
(without Higgs in the Marcolli-vS lineage as corrected by Perez-Sanchez
2024).'"

**Finding:** Neither paper makes a YM-with-Higgs continuum claim.

- Paper 30 §7 (Scope, `sec:scope`) lines 758--764 explicitly state the
  construction does NOT establish a YM mass gap on $\mathbb{R}^4$, and the
  discrete theory has no continuum limit of $\mathbb{R}^4$ inside it.
- Paper 30 §8 (Open Questions, subsection "Continuum limit" at line 859)
  asks whether the continuum limit "is SU(2) Yang--Mills on the unit $S^3$
  in the Hodge--de Rham sense" --- already YM without Higgs, already framed
  as an open question, already narrow in scope (on $S^3$, not
  $\mathbb{R}^4$).
- Paper 25 does not claim any continuum limit for the edge Laplacian that
  would correspond to YM-with-Higgs.  §V.B (open continuum-limit question)
  asks whether $L_1$ converges to the continuous Hodge-1 Laplacian on $S^3$
  — scalar Hodge theory, not gauge theory with Higgs.

**Conclusion:** No scope softening needed.  The new citations inform the scope
sections by naming the Marcolli--vS lineage as the context, but the scope
statements themselves were already Higgs-free.

---

## 5. Conflicts / scope issues

None detected.

- No prior Marcolli or van Suijlekom citations in either paper
  (verified via Grep: zero matches for `Marcolli|Perez-Sanchez|Suijlekom|
  Bratteli|2401\.03705|2508\.17338|1301\.3480` across all of `papers/`).
- Neither paper altered any computational claim, theorem statement, or
  numerical result.
- Conjectural status of Paper 2 (fine-structure constant) references is
  unchanged.
- No rhetoric-rule violation: both new positioning sentences separate
  "what Marcolli--vS / Perez-Sanchez establishes" from "what GeoVac claims,"
  consistent with CLAUDE.md §1.5.

---

## 6. Files modified

- `papers/synthesis/paper_25_hopf_gauge_structure.tex`
  - Lines 124--134: added positioning text (10 new lines of TeX)
  - Lines 1085--1105: added three bibliography entries
- `papers/observations/paper_30_su2_wilson.tex`
  - Lines 116--131: added positioning text (8 new lines of TeX)
  - Lines 963--983: added three bibliography entries

**New files:** this memo only
(`debug/alpha_sprint_a/papers_25_30_citation_update_memo.md`).

No code files touched.  No tests affected.  No computational claims altered.
