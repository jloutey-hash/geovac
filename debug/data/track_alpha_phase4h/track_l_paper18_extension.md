# Track alpha-L: Paper 18 Three-Tier Extension — Change Note

**Phase:** 4H alpha documentation sprint
**Date:** 2026-04-10
**Target file:** `papers/core/paper_18_exchange_constants.tex`
**Lines added:** ~168 lines (new file length ~1702 lines, up from 1534)
**Insertion point:** After line 830 (end of second Observation in
Section V, "Connection to the Fine Structure Constant"), before the
section break that introduces Section VI ("Observable Classification
by Transcendental Content"). The new content now occupies lines
832–998 of the revised file.

## Subsection title

`\subsection{Three-tier structure of the combination rule}`
(`\label{sec:alpha_three_tier}`)

This is Section V's first subsection; the preceding prose, the
\begin{conjecture}, and the two \begin{observation} blocks are
untouched. The subsection is formatted with four
`\paragraph{...}` blocks plus a "Reframed open question"
paragraph, matching Section VI's mixed subsection/paragraph style.

## Summary of additions

1. **\paragraph{1. A $\kappa$--$B$ identity from the Fock weight.}**
   States the six exact rational matrix elements
   $\langle n,l | (p^2 + p_0^2)^{-2} | n,l \rangle_{p_0=1}$ for the
   Paper 2 shell set and the weighted-sum identity
   $\sum (2l+1) l(l+1) \langle n,l|w|n,l\rangle_{p_0=1} = 63/4 = 6 B |\kappa|$.
   Labelled as Eq. \ref{eq:kappa_B_identity}. First exact rational
   link between $\kappa = -1/16$ and $B = 42$.

2. **\paragraph{2. A Dirichlet identity for $F$.}** Defines
   $D_{n^2}(s) = \sum_n n^2 \cdot n^{-s} = \zeta_R(s - 2)$
   (Eq. \ref{eq:Dn2_def}) and gives
   $D_{n^2}(d_{\max}) = \zeta_R(2) = \pi^2/6 = F$
   (Eq. \ref{eq:F_dirichlet}). Three independent natural readings
   of $s = 4$ (d_max from Paper 0, ambient dim of R^4, $2 N_{init}$).

3. **\paragraph{3. Three-tier classification.}** Small tabular
   (4 rows, 4 cols: Component / Value / Type / Origin) summarizing
   the categorical differences between B (finite sum), F (infinite
   Dirichlet), and Delta (boundary product). Concludes that K
   assembles objects from three different tiers rather than
   specializing a single generator.

4. **\paragraph{4. Negative results on candidate mechanisms.}**
   Seven-item enumeration of eliminated mechanisms: Hopf-twist,
   higher S^3 Casimirs, Hopf graph quotient, discrete S^1 fiber
   zeta, continuous S^1 fiber via Mellin/heat-kernel (with Jacobi
   inversion explanation), S^5 spectral geometry, 3,228-candidate
   ζ-scan for Delta. Followed by the structural reason
   (round-sphere LB operators never produce π² + rational
   additively) and the Paper 24 Bargmann–Segal consistency note.

5. **\paragraph{Reframed open question.}** Restates the open
   question as why the additive combination $B + F - \Delta$,
   multiplied by $\omega_2 = \pi$, equals $\alpha^{-1}$ to
   $8.8 \times 10^{-8}$—i.e., a conspiracy question about three
   categorically different objects, not a derivation question
   about any single component. Explicit final sentence:
   "Paper~2's identification remains conjectural."

## Complications / notes

- **Bibliography:** There is no `loutey_paper24` entry in any
  Paper 18 bib file. The initial draft contained
  `\cite{loutey_paper24}`; I replaced it with a bare `Paper~24`
  reference (no cite) to avoid a `?` in the compiled output.
  If Paper 24 later gets a bib entry, this can be upgraded.
- **Section V had no \subsection before.** Section V was flat
  prose ending with two observations. Adding the first
  \subsection creates a slight visual asymmetry but matches
  Section VI's (and the rest of the paper's) standard format.
  The overall section structure is unchanged; only the new
  subsection is appended.
- **Equation labels:** Three new equation labels added
  (`eq:fock_weight_values`, `eq:kappa_B_identity`, `eq:Dn2_def`,
  `eq:F_dirichlet`). None collide with existing labels in the
  paper (verified by grep).
- **Environment balance:** 4 equation pairs, 1 center pair,
  1 tabular pair inserted, all properly matched.
- **Conjecture status:** The original \begin{conjecture} at line
  791 and the two \begin{observation} blocks at lines 803 and 819
  are untouched. The new subsection explicitly says "Paper 2's
  identification remains conjectural" at its final sentence and
  hedges throughout with "we observe," "to our knowledge,"
  "appears to occupy," "partially answered," etc.
- **Tone:** Matches Paper 18's existing observational voice.
  Avoids any assertive claim about Paper 2's cubic identity.
  All data is credited to "an ongoing computational investigation"
  with "exact-rational symbolic computation" and
  "high-precision arithmetic"—no specific track IDs are cited
  in the published text.

## Sources

- Phase 4B α-C: `debug/data/track_alpha_phase4b/track_c_analysis.md`
  (Fock weight matrix elements, 63/4 identity, κ·B link)
- Phase 4F α-J: `debug/data/track_alpha_phase4f/track_j_analysis.md`
  (D_{n^2}(d_max) = ζ(2) = F identity)
- Phase 4G α-K: `debug/data/track_alpha_phase4g/track_k_analysis.md`
  (Δ irreducibility, 3,228 candidate ζ-scan, Laurent expansion)
- Phase 4D α-H: `debug/data/track_alpha_phase4d/track_h_analysis.md`
  (Mellin / heat-kernel π-linearity mechanism)
