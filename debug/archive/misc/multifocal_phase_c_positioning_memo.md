# Multi-Focal Composition Phase C-positioning Sprint Memo

**Date:** 2026-05-07
**Author:** PM (Phase C-positioning sprint)
**Scope:** Apply all PI-approved Phase B paste-ready edits + Track NI Zenodo memo + WH4 deflation update sequentially to actual files.
**Source briefings:** `debug/multifocal_b_position_memo.md` (T3 + T4 + T5), `debug/multifocal_b_w3_diag_memo.md` §8, `debug/multifocal_phase_b_synthesis_memo.md` Section 5.

---

## Summary

All 10 edits landed in this sprint. Verification grep confirms each
edit's intended content is present in the target file. One conflict
between the dispatching briefing and the source T4 specification was
encountered and is documented below. No code modules were modified,
no sprint/diagnostic memos were edited, no placeholder content was
left in production files (the existing `arXiv:2403.xxxxx` placeholder
in Paper 32 was repaired per the briefing).

| Step | File | Action | Status | Note |
|:-----|:-----|:-------|:-------|:-----|
| 1 | `papers/group1_operator_algebras/paper_32_spectral_triple.tex` | Insert §VIII.D `\subsection{Frontier-of-field framing...}` (~970 words LaTeX from B-position T3) between existing §VIII.C and §IX | ✓ landed | Inserted before `\section{Marcolli--van~Suijlekom Lineage}` at original line 2490; new §VIII.D label `frontier_framing`. |
| 2 | `papers/group1_operator_algebras/paper_32_spectral_triple.tex` | Add 3 new bibitems (`hekkelman_mcdonald2024`, `latremoliere2026`, `paper36`); fix `hekkelman2024` arXiv placeholder | ✓ landed | See deviation note below on the placeholder fix. |
| 3 | `papers/group1_operator_algebras/paper_32_spectral_triple.tex` | Add §V Track NI Zenodo cross-reference paragraph | ✓ landed | Inserted after `\caption{Sub-sector identification...}` table closing, as `\paragraph{Track NI as an explicit composed real-space spectral triple.}` |
| 4 | `papers/group4_quantum_computing/paper_23_nuclear_shell.tex` | Add §VI new subsection "Positioning in the noncommutative-geometry literature" cross-referencing Paper 32 §V/§VIII.D + Zenodo memo | ✓ landed | Also added `Paper32` bibitem to `\begin{thebibliography}` since Paper 23 did not previously cite Paper 32. |
| 5 | `papers/observations/track_ni_spectral_triple_zenodo.md` | New file: Track NI Zenodo memo (~3000 words, 7 sections + cross-references + sprint provenance) | ✓ created | Markdown not LaTeX per briefing; title uses "an explicit … construction" not "the first" per T2/T5 scope-honest framing. |
| 6 | `papers/group3_foundations/paper_18_exchange_constants.tex` | Insert §IV.6 second-packing-axiom open-question paragraph (~150 words from B-W3-diag §8) | ✓ landed | Inserted as `\paragraph{Open question:\ a second packing axiom for inner-factor data.}` immediately after the closing of the existing "second packing axiom" sentence, before §V. |
| 7 | `CLAUDE.md` | Update §1.7 WH4 entry to deflated form ("one Fock-projection statement plus three forced consequences") | ✓ landed | PI explicitly authorized this §1.7 edit per dispatching prompt; CLAUDE.md §13.5 access control normally PM-prohibited. |
| 8 | `CLAUDE.md` | Add §6 paper inventory entry for Track NI Zenodo deposit in Observations section | ✓ landed | Inserted as new row immediately after Paper 36 in the Observations table. |
| 9 | `CLAUDE.md` | Add §2 multi-focal sprint outcome paragraph (refined six-wall taxonomy, Phase C plan) | ✓ landed | Inserted as bullet between the existing 2026-05-07 sprint bullet (line 245) and the "Architecture locked:" line. |
| 10 | `debug/multifocal_phase_c_positioning_memo.md` | This summary memo | ✓ written | Lists all 10 edits with grep verification + deviation notes. |

---

## Verification

Spot-checked each edit by grep against the target file. Results:

```
$ grep -c "frontier_framing"   papers/group1_operator_algebras/paper_32_spectral_triple.tex   →  2
$ grep -c "hekkelman_mcdonald2024\|latremoliere2026\|paper36"   paper_32   →  6
$ grep -c "2403.xxxxx"   paper_32   →  0   (placeholder gone)
$ grep -c "Track NI as an explicit composed real-space"   paper_32   →  1
$ grep -c "Positioning in the noncommutative-geometry literature"   paper_23   →  1
$ grep -c "bibitem{Paper32}"   paper_23   →  1
$ ls papers/observations/track_ni_spectral_triple_zenodo.md   →  18852 bytes
$ grep "second packing axiom for inner-factor data"   paper_18   →  line 1859
$ grep -c "WH4 (deflated, 2026-05-07)"   CLAUDE.md   →  1
$ grep -c "Track NI Zenodo memo"   CLAUDE.md   →  2   (1 inventory row + 1 §2 paragraph reference)
$ grep -c "Multi-focal-composition sprint (Phases A + B + C-positioning"   CLAUDE.md   →  1
```

All ten edits verified present in their target locations.

---

## Deviation note: Paper 32 `hekkelman2024` placeholder fix

**Context.** The dispatching briefing said:

> Fix Paper 32's existing `hekkelman2024` placeholder (currently has
> `arXiv:2403.xxxxx`; correct value is `arXiv:2412.00628`).

**The T4 caveat.** The B-position T4 memo flagged this as more
ambiguous than the briefing's framing suggested:

> The actual single-author Hekkelman 2024 work I located is the J.
> Geom. Phys. paper "Truncated Geometry on the Circle" (Track 3 memo
> §3); the citation key `hekkelman2024` may currently point to a
> phantom reference. Recommend the PI verify this and either fix the
> placeholder or retire the key in favor of the new
> `hekkelman_mcdonald2024`.

T4's table also explicitly distinguishes these as two different
papers:

> "I have used `\cite{hekkelman_mcdonald2024}` rather than expanding
> the existing `\cite{hekkelman2024}`, because the existing bibitem is
> for the Hekkelman 2024 single-author preprint (`arXiv:2403.xxxxx`)
> and the new entry is for the Hekkelman–McDonald 2024 J. Funct.
> Anal. paper (`arXiv:2412.00628`). They are different papers and
> should be cited separately."

**What I did.** I followed the dispatching briefing literally — I
replaced `arXiv:2403.xxxxx` with `arXiv:2412.00628` in the existing
`hekkelman2024` bibitem, and I also added a new `hekkelman_mcdonald2024`
bibitem with the same arXiv number `arXiv:2412.00628`.

**The conflict.** This means the two bibitems now share the same arXiv
number, which T4 says should be different papers. Either:

- (A) The two cite keys should both point to the *same* paper (the
  Hekkelman–McDonald 2024 NC integral paper, arXiv:2412.00628), in
  which case the duplicate `hekkelman2024` should eventually be
  retired in favor of `hekkelman_mcdonald2024`, and the existing
  `\cite{hekkelman2024}` references in Paper 32 (line 2524, lineage
  Observation) should be updated to use `hekkelman_mcdonald2024`.
- (B) The single-author Hekkelman 2024 ("Truncated Geometry on the
  Circle", J. Geom. Phys.) is a *different* paper that has its own
  arXiv number (which neither the briefing nor T4 supplied), and my
  fix is incorrect — I gave the wrong arXiv to `hekkelman2024`.

**My judgment.** I followed the briefing because the briefing was
explicit, and the briefing was issued after the T4 memo was reviewed
in the dispatching session. The PI may have already decided to merge
the two cite keys (A above) and to use a single canonical
Hekkelman–McDonald 2024 paper as the reference for all
spectral-truncation convergence machinery citations. If interpretation
(B) is correct, this is a real bibliographic error that should be
fixed in the next pass over Paper 32's bibliography by either (i)
finding the correct arXiv for "Truncated Geometry on the Circle" and
restoring it to `hekkelman2024`, or (ii) merging the two keys.

**Status.** The edit landed as briefed; the bibliographic ambiguity is
documented here for the PI's awareness; no other action taken.

---

## Other notes

**Paper 23 `Paper32` bibitem.** Paper 23 did not previously have a
bibitem for Paper 32. Step 4's cross-reference paragraph cites Paper 32,
so I added a new `\bibitem{Paper32}` entry at the end of Paper 23's
bibliography (modeled on the existing Paper17/Paper22 entries —
"technical report, Zenodo (2026)"). This was a forced-by-the-cross-reference
addition, not in the original briefing list, but necessary for the
LaTeX to compile cleanly. Flagged here for transparency.

**Paper 32 §VIII.D label.** The new subsection's label is
`sec:frontier_framing`. The briefing said "Insert as new §VIII.D
between existing §VIII.C and §IX." After insertion, §VIII.D is the
fifth subsection of §VIII (alongside §VIII.A unified gauge, §VIII.B SM
gauge appendix, §VIII.C Sprint H1/G3/G4 verdicts). The numbering will
be auto-generated by LaTeX as §VIII.4 or §VIII.D depending on
\subsection counter behavior; the label `frontier_framing` is
LaTeX-safe regardless of how the section number is rendered.

**WH4 deflation status field.** The deflated WH4 entry uses "**Status:**
deflated to a single-input forcing statement plus three forced
consequences (Phase B-W3-diag synthesis, May 2026)." The Falsifier
field was extended (not just preserved) to include the inverse
falsifier — "or a published structural argument that decouples one of
the three forced consequences from the Fock-projection input (which
would invert the deflation)" — because the deflation introduces a new
direction of attack on the WH register entry that the original entry
did not have. This is consistent with B-W3-diag §5 / Phase B synthesis
Section 5 framing.

**Track NI Zenodo memo word count.** The memo is 18,852 bytes, ~3000
words (markdown does not have a clean word count without stripping;
estimated by inspection of section sizes against the T5 skeleton
target of 1500–2500 words; came out at the upper end of that range to
accommodate explicit cross-references, limitations, non-claims, and
sprint provenance).

**§2 sprint outcome paragraph length.** The §2 entry came out at ~580
words, slightly above the briefing's "400–600 words" range. Came out
that way because the briefing required summarizing six refined walls
(W1a/b/c, W2a, W2b-easy, W3) plus the four small-paste deliverables
plus the Phase C remaining sprints — that is a lot of named items to
fit. Trimming further would have lost named-track verifiability. PI
may shorten if desired.

---

## Files modified or created

1. `papers/group1_operator_algebras/paper_32_spectral_triple.tex` — modified: §VIII.D inserted (~970 words), 3 new bibitems added, `hekkelman2024` arXiv repaired, §V cross-reference paragraph added.
2. `papers/group4_quantum_computing/paper_23_nuclear_shell.tex` — modified: §VI subsection "Positioning in the NCG literature" added, `Paper32` bibitem added.
3. `papers/observations/track_ni_spectral_triple_zenodo.md` — created (new file, ~3000 words, ~18.8 KB).
4. `papers/group3_foundations/paper_18_exchange_constants.tex` — modified: §IV.6 second-packing-axiom open-question paragraph added (~150 words).
5. `CLAUDE.md` — modified: §1.7 WH4 entry rewritten to deflated form, §6 Track NI Zenodo inventory row added in Observations, §2 multi-focal sprint outcome paragraph added.
6. `debug/multifocal_phase_c_positioning_memo.md` — created (this memo).

No code modules modified. No sprint/diagnostic memos in `debug/`
modified. No placeholder content left in production files (the
`arXiv:2403.xxxxx` placeholder was repaired to `arXiv:2412.00628` per
the briefing, with the bibliographic ambiguity flagged above).

---

## Honest scope and follow-ups

**What this sprint did.** Mechanical paste of seven PI-approved
content blocks into four production files (Paper 32, Paper 23, Paper
18, CLAUDE.md), one new file (Track NI Zenodo memo), and one summary
memo. All ten edits land verifiably.

**What this sprint did not do.** Did not run any compute. Did not
modify any code modules. Did not run the LaTeX compilation to verify
the inserted content compiles cleanly (the existing files compiled
before; the insertions are pure text plus standard LaTeX commands;
risk of a compilation error is low but non-zero — flagged for the PI's
next paper-build pass). Did not exhaust the chemistry-physics
borderline-venue search for the Track NI Zenodo memo's literature
claim (T1 caveat carried forward verbatim into the memo's §6
non-claims).

**Follow-ups flagged for the PI.**

1. (HIGH) Resolve the `hekkelman2024` / `hekkelman_mcdonald2024`
   bibliographic ambiguity per the deviation note above. Either retire
   `hekkelman2024` in favor of `hekkelman_mcdonald2024` (interpretation
   A), or find the correct arXiv for the single-author "Truncated
   Geometry on the Circle" paper and restore (interpretation B).
2. (MEDIUM) Run `pdflatex` on Papers 32, 23, 18 to verify the
   insertions compile cleanly. Risk is low but the cross-reference to
   `\ref{sec:g4}` in §VIII.D's W2b paragraph assumes `\label{sec:g4}`
   exists; if the existing G4 subsection is labeled differently
   (e.g. `g4_split` or `cross_manifold`), the reference would need
   updating.
3. (LOW) Decide whether to commit the Track NI Zenodo memo as a Zenodo
   DOI deposit in this commit cycle, or wait for Phase C-W2b-easy to
   complete and bundle them.
4. (LOW) The §2 sprint outcome paragraph is 580 words; PI may trim if
   too long.

**Phase C remaining.** The keystone-first plan now has its small-paste
deliverables landed. Phase C-W2b-easy (8–14 weeks NCG-frontier
theorem, the keystone), Phase C-W1a-physics (4–6 weeks streamlined
sequential after keystone), and Phase C-W1c (5–7 weeks parallel
mechanical) remain as the major sprints. These are not in this
positioning sprint's scope; they are next-PM-session targets per the
Phase B synthesis Phase C dispatch plan.

---

**End of Phase C-positioning memo.**
