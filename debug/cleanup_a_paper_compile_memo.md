# Cleanup Track A — Paper compile fixes

**Date:** 2026-05-23
**Scope:** Fix three pre-existing LaTeX compile failures in Papers 17, 32, 34
**Verdict:** **ALL-FIXED-3-PASS-CLEAN**

## §1. Per-paper diagnosis

### Paper 17 (`papers/group2_quantum_chemistry/paper_17_composed_geometries.tex`)

**Error (line 11):**
```
! Undefined control sequence.
\@AF@join ...xpandafter \@title@join@ \@title@aux
l.11 \thanks{Paper 17 in the Geometric Vacuum series}
```

**Root cause:** revtex4-2 requires the title block (`\title`, `\thanks`,
`\author`, `\affiliation`, `\date`) to appear **inside** the `document`
environment. The original Paper 17 layout placed the title block in the
**preamble** (before `\begin{document}`), then opened the document
right before the abstract:

```latex
\title{...}
\thanks{Paper 17 in the Geometric Vacuum series}   % <-- in preamble, fails
\author{...}
\affiliation{...}
\date{...}
\begin{document}                                    % <-- too late
\begin{abstract}
```

`\thanks` is a `revtex` macro that requires the title-block apparatus
to be already active; in the preamble, `\thanks` is just a generic
LaTeX footnote command and tries to chain through revtex's
`\@AF@join` machinery, which is only bound after `\begin{document}`.
The companion Papers 18 and 19 (same author block style) had this
structure correctly already (`\begin{document}` first, then `\title`),
which is why they compile cleanly with the apparently-identical
`\thanks{...}` line.

**Pre-existing:** yes. The `\begin{document}` mis-placement predates
all 2026-05 edits; nothing in today's modular-arc paper-edit sprint
moved this line.

### Paper 32 (`papers/group1_operator_algebras/paper_32_spectral_triple.tex`)

**Error (line 1362, AND line 3794 found by the post-fix scan):**
```
! Missing $ inserted.
l.1362 .../track\_ni\_spectral\_triple\_zenodo.md}
```

**Root cause:** In each affected `\texttt{...}` argument, the underscore
escapes were applied to ALL underscores except the FIRST. The original
strings read

```latex
\texttt{papers/group4_quantum_computing/track\_ni\_spectral\_triple\_zenodo.md}
\texttt{papers/group1_operator_algebras/paper\_43\_lorentzian\_extension\_outline.md}
```

The unescaped underscores in `group4_quantum_computing` and
`group1_operator_algebras` cause `_` to fall through to math-mode
subscript syntax in horizontal mode, which raises "Missing `$`
inserted." Once one underscore is unescaped, every subsequent character
is interpreted in math mode and the rest of the `\texttt{}` content
breaks.

**Pre-existing:** yes. Both lines were authored before today and
nothing in the modular-arc paper-edit sprint modified the relevant
file regions. Confirmed via `git log` (last touch was v2.51.0 paper
folder reorganization, which copied content verbatim).

### Paper 34 (`papers/group6_precision_observations/paper_34_projection_taxonomy.tex`)

**Error (line 2446):** Identical failure mode to Paper 32:
```latex
\texttt{papers/group1_operator_algebras/paper\_43\_lorentzian\_extension\_outline.md}
```
The first underscore in `group1_operator_algebras` is unescaped.

**Pre-existing:** yes. Same provenance pattern as the Paper 32
instances.

## §2. Per-paper fix applied

### Paper 17 — diff summary

Moved `\begin{document}` from line 21 to line 9 (right after
`\newtheorem` declarations), so it now precedes the title block.
Removed the duplicate `\begin{document}` line that originally sat
before the abstract.

```latex
%% ========================================================================
\begin{document}                                       % <-- new position

\title{Composed Natural Geometries ...}
\thanks{Paper 17 in the Geometric Vacuum series}       % <-- now inside document
\author{J.~Loutey}
\affiliation{Independent Researcher, Kent, Washington}
\date{March 21, 2026}

%% ========================================================================
%% ABSTRACT
%% ========================================================================
\begin{abstract}                                       % <-- \begin{document} removed from here
```

This brings Paper 17's preamble layout in line with Papers 18 and 19.
No substantive content touched.

### Papers 32 and 34 — diff summary

Substituted `\texttt{path/with_underscores}` → `\verb|path/with_underscores|`
in three locations:

- Paper 32 line 1362: `track_ni_spectral_triple_zenodo.md` path
- Paper 32 line 3794: `paper_43_lorentzian_extension_outline.md` path
- Paper 34 line 2446: `paper_43_lorentzian_extension_outline.md` path

The `\verb|...|` form is the canonical LaTeX idiom for typesetting
verbatim content with metacharacters; it suppresses underscore
interpretation entirely without per-character escaping. It renders
identically to `\texttt{}` for these paths. The choice of pipe
delimiter `|` is safe because the contents never contain `|`.

The two existing in-line literals on Paper 34 lines 2441/2442/2443
already used `\verb|...|` correctly — the line 2446 instance was the
odd one out using `\texttt{}` with hand-escaped underscores; my fix
makes it consistent with its neighbors.

## §3. Compile verification

Three-pass `pdflatex -interaction=nonstopmode` results (post-fix):

| Paper | Errors pass 3 | Pages | PDF size (B) | Verdict |
|:------|:--------------|:------|:-------------|:--------|
| 17 | **0** | 13 | 489 273 | CLEAN |
| 32 | **0** | 53 | 664 836 | CLEAN |
| 34 | **0** | 112 | 1 300 927 | CLEAN |

Pass-to-pass byte sizes converged identically between passes 2 and 3
for Papers 32 and 34 (cross-references resolved); Paper 17 also
converged (same byte size).

Helper diagnostic added: `debug/scan_unescaped_underscores.py` scans
papers for unescaped underscores inside `\texttt{}`. After the three
fixes, it returns 0 hits across all six group folders + synthesis. This
script is a one-off cleanup tool; it can stay in `debug/` for future
audits.

## §4. Status of all group papers (regression check)

### Explicitly verified clean (this sprint)

| Paper | Folder | Pass-3 errors | PDF status |
|:------|:-------|:--------------|:-----------|
| Paper 17 | group2 | 0 | ✓ produced |
| Paper 18 | group3 | 0 | ✓ produced |
| Paper 19 | group2 | 0 | ✓ produced |
| Paper 32 | group1 | 0 | ✓ produced |
| Paper 34 | group6 | 0 | ✓ produced |
| Paper 39 | group1 | 0 | ✓ produced |
| Paper 43 | group1 | 0 | ✓ produced |
| Paper 44 | group1 | 0 | ✓ produced |
| Paper 45 | group1 | 0 | ✓ produced |
| Paper 46 | group1 | 0 | ✓ produced |
| Paper 47 | group1 | 0 | ✓ produced |

All Papers 18 / 19 (regression target, explicitly named in prompt) and
the recent math.OA papers I personally added a pass over compile
cleanly and produce PDFs.

### Out-of-scope pre-existing issues found (flagged but NOT fixed)

These were already broken before this sprint and are unrelated to the
three target failures. Per prompt instructions, I document them and
take no action.

| Paper | Issue | Pre-existing? |
|:------|:------|:--------------|
| Paper 38 | `pdfTeX error (font expansion): auto expansion is only possible with scalable fonts` | Yes — known MiKTeX environmental issue, CLAUDE.md §6 Paper 42 entry already notes "microtype disabled to avoid MiKTeX font-expansion environmental issue (same fix as Papers 38/39/40)." Paper 38 still has microtype enabled. Local MiKTeX install lacks scalable-font config. |
| Paper 40 | Same `pdfTeX font expansion` error | Same root cause as Paper 38. |
| Paper 42 | 6 math-mode errors (`Missing $`, `Extra }`) emitted around p19. PDF still produced (686 KB) so errors are non-fatal in the user's environment, but flagged here. | Pre-existing — none of today's edits touched the relevant region. |

Recommended follow-on (out of scope for Track A): a small "microtype
audit + Paper 42 math-mode fix" sprint that mirrors what was
already done for Papers 39, 42, 45, 46, 47 (which compile clean).
Estimate: ~30 min.

### Not verified (out of scope)

I did not compile Papers in `papers/archive/` or `papers/synthesis/`,
or papers from groups 4 and 5 (no prompt-mandated check). The
unescaped-underscore scanner confirmed zero issues of *that specific
class* across all six group folders + synthesis; other classes of
issue may exist.

## §5. Pre-existing issues NOT in scope

Two issue classes surfaced during the verification sweep that are out
of scope for Track A but worth recording:

1. **Microtype / pdfTeX scalable-fonts (Papers 38, 40).** These two
   papers fail with `pdfTeX error (font expansion): auto expansion is
   only possible with scalable fonts`. The fix is to either disable
   `microtype` (as already done in Papers 39, 42, 45, 46, 47 per
   CLAUDE.md §6) or install scalable-font support in the local MiKTeX
   tree. The CLAUDE.md note for Paper 42 (§6 entry) explicitly says
   "(same fix as Papers 38/39/40)" — that fix appears to have NOT been
   applied to Papers 38 and 40 themselves. Mechanical fix.

2. **Paper 42 math-mode escapes around page 19.** Six `! Missing $`
   and `! Extra }, or forgotten $` errors appear at compile time, but
   PDF still produces. Pre-existing; not introduced by today's sprint.
   Likely a `\verb` interaction with a math expression or an
   unescaped `<` / `>` in math content. Requires a focused diagnostic
   read of the affected region.

Neither is a Track A target.

### Git status

`git status --short` shows the expected target files modified
(`paper_17_composed_geometries.tex`, `paper_32_spectral_triple.tex`,
`paper_34_projection_taxonomy.tex`) plus generated PDFs and the
helper scanner script. Other modified files in the working tree
(Papers 18, 19, balanced_coupled.py, multi_zeta_orbitals.py, etc.,
plus many untracked debug/ files) are PRE-EXISTING from prior
Sprint $\alpha$-1/2/3 and Sprint M-Y/M-Z/M-H1 work and are not part
of this sprint's intended modifications.

## Verdict

**ALL-FIXED-3-PASS-CLEAN**

All three target papers (17, 32, 34) now compile three-pass clean
with 0 errors, produce valid PDFs at expected page counts, and have
zero unescaped-underscore-in-`\texttt{}` regressions per the
diagnostic scanner. No substantive content was touched; only the
syntactic compile blockers were fixed.

Regression targets (Papers 18, 19) verified still clean.
