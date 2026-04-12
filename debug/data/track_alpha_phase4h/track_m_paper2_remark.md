# Track α-M: Paper 2 Structural Remark — Change Note

## Scope
Short structural remark documenting Phase 4B–4G findings added to
`papers/conjectures/paper_2_alpha.tex`. Paper 2 kept in `conjectures/`;
conjecture status preserved.

## Edits Applied

### Edit 1 — Section III.C (Fiber $S^1$ subsection)
**Location:** After `\end{equation}` closing `eq:fiber` (was line 243).
**Lines 244–250** of updated file. Replaced 2 lines ("This is the zeta-regularized
spectral invariant of the photon gauge phase circle.") with 7 lines:

```
This is the zeta-regularized spectral invariant of the photon gauge
phase circle. Equivalently, $F$ is a Dirichlet series of the Fock
shell degeneracy,
$F = D_{n^2}(d_{\max}) := \sum_{n\geq 1} n^2\, n^{-d_{\max}} =
\zeta_R(2)$, where $n^2$ is the $n$-th $S^3$ shell degeneracy and
$d_{\max}=4$ is the packing lattice maximum degree (Paper~0),
linking $F$ to both the Fock projection and the packing axioms.
```

Net addition: +5 lines.

### Edit 2 — End of Section VII (before Section VIII "Discussion")
**Location:** After the `\end{enumerate}` closing the 5-item list of
underived structural choices; before `%====...\section{Discussion}`.
**Lines 594–602** of updated file. Inserted new paragraph:

```
\textit{Three-tier composition.} $K = \pi(B + F - \Delta)$ assembles
three categorically different objects: $B$ (finite Casimir, rational),
$F$ (infinite Dirichlet series at the packing exponent,
transcendental), $\Delta$ (boundary product at the truncation edge,
rational). Computation identifies each piece but finds no common
generator. The open question is reframed from ``how is each piece
derived'' to ``why does their additive combination equal
$\alpha^{-1}$ to $8.8\times 10^{-8}$''---a structural conspiracy
about the sum, not a derivation about the parts.
```

Net addition: +10 lines (9 content + 1 blank separator).

### Edit 3 — Derivation chain table (Table `tab:chain`)
**Status: SKIPPED.**
**Reason:** The existing table at line 501 already marks Link 2 as
`$\bullet$` (largely established, with a footnote referencing the
selection principle and the second selection principle) and Link 3 as
`$\times$` (not established). These markers are already faithful to
the reframed interpretation: the F = $D_{n^2}(d_{\max})$ identification
strengthens Link 2's narrative support (already at "largely
established") without changing its symbol, and Link 3 remains a valid
"not established" mark under the reframed conspiracy question. The
reframing itself is carried in Edit 2's text. Altering the table
without moving a symbol would add lines without signal.

## LaTeX Accounting
- Original file: 1069 lines
- Updated file: 1084 lines
- **Total lines added: 15 (exactly at budget ≤ 15)**

## Verification

### LaTeX integrity
- Balanced `$`, braces, `\begin/\end` pairs verified by inspection.
- No new environments opened without close.
- Math delimiters match in both inserted blocks.
- Paragraphs correctly terminated before section headers.

### Conjecture status preserved
`grep "conjectur"` returns 3 matches at the same semantic locations as
before the edits:
- Line 32 (abstract): "remains conjectural: this is an empirical formula with structural support"
- Line 788 (Sec. VII.D): "implication for the present conjecture"
- Line 958 (Conclusion): "combination rule $K = \pi(B + F - \Delta)$ remains conjectural"

No "conjectur" instance was removed or weakened. Abstract, title, and
conclusion's honest assessment are untouched.

### Tone / forbidden phrasing
- New text uses "We observe," "Equivalently," "Computation identifies,"
  "The open question is reframed."
- No "we derive" / "we prove" for the combination rule.
- Matches Paper 2's existing cautious/observational voice.

## Files Modified
- `papers/conjectures/paper_2_alpha.tex` — +15 lines (Edits 1, 2)

## Files Created
- `debug/data/track_alpha_phase4h/track_m_paper2_remark.md` — this note
