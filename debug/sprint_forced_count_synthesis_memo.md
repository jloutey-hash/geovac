# Sprint Forced-Count Synthesis memo (2026-06-03)

## TL;DR

Paper 32 §VIII gets a new Theorem block (`thm:forced_count` +
`cor:bertrand_analogy` + `rem:forced_count_with_seam`) that
crystallizes a result already distributed across eight sub-sprint
paragraphs in the same section. **No new mathematics.** The
deliverable is to name what was already proved as the Forced-Count
Theorem and frame it as the inner-factor analogue of Bertrand.

## Sprint context

The user's "deep move" scoping conversation (this session) identified
three candidate deep moves on the Yukawa frontier:

1. **Forced-count theorem** — name the count $\dim \mathcal{M}(D_F)$ as
   a structural consequence of forced inputs (sprint-scale, low risk,
   modest novelty).
2. **Force $N_{\mathrm{gen}}$ via Hopf-tower-to-representation
   extension** — multi-month, real risk of falsification, paper-grade
   prize.
3. **Force Yukawa values** — theorem-blocked at the AC level (Door 4
   wall); multi-year, no current shortcut.

The PI picked sequence "1 first, then 2." This memo covers (1).

On reading Paper 32 §VIII before drafting, I found that **all the
ingredients were already there** — the Sprint H1 paragraph states the
$1024 \to 512 \to 128$ chain explicitly; the Door 4b/4d/4e/4f
paragraphs prove the algebra forcing; the Direction 2 paragraph
(2026-06-03, today) establishes the packing-unreachable residue; and
the Yukawa-PSLQ memo (also today) supplies empirical confirmation of
the Door 4 wall. **What was missing was a single Theorem-block
crystallization** with a name, an explicit hypothesis list, the
reduction chain as an equation, and the Bertrand analogy stated as a
corollary.

## What was added

Insertion point:\ Paper 32 §VIII, after the Direction 2 paragraph
(closing at `\texttt{debug/seam\_packing\_scoping\_memo.md}`) and
before the `\paragraph{Files.}` block.

Three new environments:

### `thm:forced_count` (Forced $D_F$ moduli dimension)
Hypotheses:
1. Bertrand × complex-Hopf-tower $n \le 3$ + Upgrade B (Door 4e
   sphere-Lie-group axiom) forces
   $\mathcal{A}_F = \mathbb{C} \oplus \mathbb{H} \oplus M_3(\mathbb{C})$.
2. Canonical CCM SM representation on $\mathcal{H}_F$ with
   $N_{\mathrm{gen}}$ external input.
3. Standard real-spectral-triple axioms (Hermiticity, KO-dim-6
   chirality, $J$-reality, order-one).

Conclusion: $\dim_{\mathbb{R}} \mathcal{M}(D_F) = 128$ at
$N_{\mathrm{gen}} = 1$, via chain $1024 \to 512 \to 128 \to 128$;
$N_{\mathrm{gen}}^2$ scaling for full mixing; 8 per-gen for diagonal
slice.

Proof sketch:\ direct counting at each step, with the Hermiticity,
chirality, $J$-reality, and order-one reductions explicitly named.
Verified bit-exactly at $n_{\max} \in \{2, 3\}$ in
`geovac/standard_model_triple.py` (45 tests passing) — i.e., this is
not new computation, it cites the existing H1 module.

### `cor:bertrand_analogy` (Bertrand analogy)
States the structural identity:\ Bertrand forces the FORM ($-k/r$ or
$kr^2$) without forcing the COUPLING $k$; Theorem~\ref{thm:forced_count}
forces the DIMENSION of $\mathcal{M}(D_F)$ without forcing the
$D_F$-entries.  Same Bertrand seam.

### `rem:forced_count_with_seam` (Negative companions)
Pairs the positive theorem with two negatives:
- **Door 4 seam theorem + Yukawa-PSLQ empirical confirmation:** values
  inside $\mathcal{M}(D_F)$ live in a Dirichlet ring disjoint from
  every forced outer ring (theorem); $162$-cell PSLQ sweep at
  $M \le 10^3$ returned zero hits (empirical, this same session).
- **Direction 2 packing-reach theorem:** $N_{\mathrm{gen}}$ and inner
  KO-dim are packing-unreachable in the current Paper 0 construction.

Net characterization of the full forced/free seam at the $D_F$ level:
- FORCED: the moduli space $\mathcal{M}(D_F)$ and its dimension.
- FREE: (i) point in moduli space (Yukawa values), (ii) multiplicity
  ($N_{\mathrm{gen}}$), (iii) real-structure signature (inner KO-dim).

## What this is not

This is not a new theorem — every ingredient is in the H1, Door 4, and
Direction 2 memos / Paper 32 paragraphs. It is a **crystallization**.

This is not a route to forcing Yukawa values. The values remain on the
free side of the Door 4 seam.  Read 2 (Hopf-tower-to-representation
extension to force $N_{\mathrm{gen}}$) is the next deep move; this
sprint is the language-setup for it.

## Compile-check

Paper 32 compiles clean (64 pages, no LaTeX errors). The new
environments use `\begin{theorem}` / `\begin{corollary}` /
`\begin{remark}` which are pre-defined in the Paper 32 preamble.

## Cross-references

**Strengthens / consolidates:**
- `debug/h1_higgs_inner_fluctuation_memo.md` — H1 sprint that produced
  the $1024 \to 128$ chain.
- `debug/sprint_door4_series_closure_memo.md` — Door 4 series that
  forced $\mathcal{A}_F$.
- `debug/seam_packing_scoping_memo.md` — Direction 2 scoping that named
  $N_{\mathrm{gen}}$ / KO-dim as packing-unreachable.
- `debug/sprint_yukawa_pslq_memo.md` — empirical companion (this
  session).

**Updates:**
- Paper 32 §VIII gets the three new environments (theorem, corollary,
  remark) between the Direction 2 paragraph and the Files paragraph.
- CLAUDE.md §2 gets a one-liner.
- No memory entry (the sprint outcome is in §2 + the theorem
  reference; no cross-session fact not derivable from those).

**Does NOT change:**
- The H1 Yukawa non-selection theorem (still the load-bearing forcing
  result; the new Theorem just consolidates and Bertrand-frames it).
- The Door 4 seam theorem (still the value-forcing wall).
- The Direction 2 NO-GO on packing (still the $N_{\mathrm{gen}}$ /
  KO-dim wall).
- WH1 PROVEN status (this is downstream of Paper 38, not affected).

## Open follow-on

**Read 2 (forced $N_{\mathrm{gen}}$, the candidate Hopf-tower-to-
representation extension):**  scoping target. The Theorem above
explicitly imports $N_{\mathrm{gen}}$ as external; the Direction 2
paragraph names it as packing-unreachable. A "sibling axiom" — new
primitive emitting fiber multiplicities — is the only structurally-
honest door. Multi-month, paper-grade.  Not pursued this session.

## 6. Honest scope

- **Theorem-grade content (given its hypotheses)**: `thm:forced_count`
  is at theorem-grade conditional on three inputs being themselves
  proven: (1) Bertrand × Hopf-tower + Upgrade B forcing of
  $\mathcal{A}_F$ (proven by Door 4b/4d/4e/4f at PARTIAL-DOOR FINAL,
  with Upgrade B's sphere-Lie-group axiom adopted as foundational);
  (2) canonical CCM SM rep with $N_{\mathrm{gen}}$ as external input
  (Direction 2 names this as packing-unreachable, i.e. genuinely
  external); (3) standard real-spectral-triple axioms. The proof is
  direct counting + bit-exact verification at $n_{\max} \in \{2, 3\}$ in
  the existing `geovac/standard_model_triple.py` module (45 tests
  passing, established in Sprint H1). The Theorem itself is a synthesis
  / naming move: no new computation, no new mathematics; the content
  was distributed across H1 + Door 4 + Direction 2 prose.
- **Interpretive content**: `cor:bertrand_analogy` is structural
  framing, not a separate theorem. The Bertrand analogy is real
  (both force the form/dimension without forcing the coupling/values)
  but the Corollary doesn't compute anything beyond the framing.
- **Cross-reference content**: `rem:forced_count_with_seam` cites three
  results from elsewhere (Door 4 wall; Direction 2 NO-GO; Yukawa-PSLQ
  empirical companion); it doesn't reprove them. The Yukawa-PSLQ
  empirical companion is what's new in the Remark — the structural
  Door 4 wall now has an empirical confirmation layer at 162 cells.
- **Named open follow-ons (deferred)**: Read 2 (Hopf-tower-to-rep
  shortcut for forcing $N_{\mathrm{gen}}$) — scoped this session at
  NO-GO, see `debug/sprint_read2_n_gen_scoping_memo.md`. Read 3 (force
  Yukawa values) — theorem-blocked at the AC level by Door 4 (i) wall,
  closed multi-year deferred. Sibling-axiom direction (new Paper 0
  primitive emitting fiber multiplicities) — multi-year, no current
  handle.

## Files

- `papers/group1_operator_algebras/paper_32_spectral_triple.tex` —
  three new environments added after the Direction 2 paragraph.
- `debug/sprint_forced_count_synthesis_memo.md` — this file.
- `CLAUDE.md` §2 — one-liner.
