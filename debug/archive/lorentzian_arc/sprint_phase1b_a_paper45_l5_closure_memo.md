# Sprint Phase 1B-A — Paper 45 L5 K⁺-preservation bookkeeping closure (memo)

**Sprint:** Phase 1B-A (apply L5 K⁺-preservation bookkeeping closure to Paper 45 per the 1A diagnostic).
**Date:** 2026-05-24.
**Predecessors:**
- `debug/l3b_2d_propinquity_assembly_memo.md` (sister L5 closure under operator-norm seminorm on the strong-form natural substrate; confirms the Latrémolière-2017 §5 bookkeeping shape transports cleanly under the operator-norm Lipschitz seminorm and is purely a transcription job).
- `debug/l3a_1_lorentzian_operator_system_memo.md` (the L3a-1 trivial-multiplier K⁺-preservation result that underwrites the natural-substrate K⁺-preservation prerequisite).
- Paper 44 (Krein operator-system substrate of Paper 45).

**Status:** Three edits applied in-place to `papers/group1_operator_algebras/paper_45_lorentzian_propinquity.tex`; three-pass clean compile (20 pages, zero new errors, zero new undefined references). No new mathematics introduced.

---

## §1. Edits applied

### Edit 1: New Remark after `rem:Kplus_trivial` (in §3, L1$'$ section)

Added `\begin{remark}[Validity of the $\Kplus$-restricted tunneling pair]\label{rem:Kplus_tunnel_valid}` immediately after the existing `rem:Kplus_trivial`. The new remark states explicitly that every element of the standard Latrémolière tunneling pair $(\Bjoint, \Pjoint)$ commutes with $\JL$ at the operator level — citing Lemma L4(e) (joint Berezin Krein-positivity preservation, factor-wise via $[\gamma^0, M^{\spat}] = 0$ on chirality-doubled spatial block + $[I_{\Nt}, M^{\temp}] = 0$ on temporal slot) and the Stinespring orthogonal-projection structure of $\Pjoint$. Concludes that the $\Kplus$-restricted pair $(\Pplus\,\Bjoint(\cdot)\,\Pplus,\, \Pplus\,\Pjoint(\cdot)\,\Pplus)$ is a valid UCP tunneling pair between the $\Kplus$-restricted triples per `def:weak_form_propinquity`. Cross-references Paper 44 §3–4 for the L1$'$ substrate-level trivial-multiplier K⁺-preservation result.

### Edit 2: Expanded proof of `prop:reach_height` (in §5, L5 propinquity assembly section)

Replaced the prior four-paragraph proof with a structured five-paragraph version:
- New header sentence: `\emph{Bookkeeping under the $\Kplus$-restricted Lipschitz seminorm.}` (initially attempted `\paragraph{...}` but it triggered an amsthm proof-environment `trivlist` error; replaced with `\emph{...}` which is the canonical alternative).
- Each constituent paragraph now explicitly cites the analytical lemma: `$\mathrm{reach}_{B}$ via L4(c)`, `$\mathrm{reach}_{P}$ via L2 + L4(b)`, `$\mathrm{height}_{B}$ via L3 + L4(d) + Stein–Weiss factor-wise`, `$\mathrm{height}_{P} = 0$ (projection)`.
- Connects each constituent to the corresponding Latrémolière 2017 §3.4 / §3.5 definition (reach_B = approximate-identity deficit, reach_P = partial-inverse-duality deficit, height_B = Lipschitz-distortion bound).
- The reach_B paragraph cites the new `rem:Kplus_tunnel_valid` for the K⁺-restriction preservation.

### Edit 3: Expanded proof of `thm:main` (Main Theorem)

Replaced the prior compact 4-line proof with a 7-line version that:
- Opens with explicit citation of Latrémolière 2018/2023 §5 Theorem 5.5 (weak-form metric-spectral-triple propinquity bound for direct UCP tunneling pair) as the framework reference.
- Explicitly states that the propinquity $\Lambda$ is the infimum over UCP tunneling pairs of the four-term max (`eq:propinquity_bound`).
- Substitutes the four constituent bounds from `prop:reach_height` explicitly in the chain of inequalities.
- Cites `rem:Kplus_tunnel_valid` for the well-definedness of the $\Kplus$-restricted UCP tunneling pair.
- Records the four-term-max $\to$ dominated-term reduction with explicit justification: since $\Cthreejoint \le 1$ asymptotically (Lemma L3), $\mathrm{height}_{B} \le \mathrm{reach}_{B}$, so the max is $\gammajoint$ at qualitative-rate level. The $\Cthreejoint \cdot \gammajoint$ form is retained to record the finite-cutoff envelope-aware constant from Remark `envelope_v2`.

### Edit 4 (additional, per directive): §1.4 status update on L5 bookkeeping closure

Added a new sentence at the end of the §1.4 preamble (before the G1 entry) noting that the Latrémolière 2017 §5 $\Kplus$-preservation bookkeeping step is now closed in §5 via explicit operator-level $[\JL, \cdot] = 0$ transcription (`rem:Kplus_tunnel_valid` + expanded `prop:reach_height` proof), inheriting the trivial-multiplier $\Kplus$-preservation result of Paper 44 §3–4. Notes that prior to this transcription the step was treated as straightforward given the L3a-1 substrate result; the explicit assembly in §5 closes it formally and is the substantive bookkeeping content of the present paper. The G1 entry itself (strong-form Lorentzian propinquity without K⁺ restriction) is structurally distinct and remains open (closed downstream by Paper 46, per the existing bibliography).

The directive's "(was: 'straightforward given L3a-1'; should now read: 'closed in §5 via explicit transcription')" phrasing was implemented at the preamble level rather than mutating the G1 item text itself, because G1 refers to the strong-form question (different from the L5 K⁺-preservation bookkeeping gap) and conflating them would lose information. The new preamble sentence makes the closure visible to a reader of §1.4 without altering the meaning of G1.

---

## §2. Latrémolière 2017 citation

The directive specified "the canonical one (J. Math. Pures Appl. 2017 or arXiv preprint version) consistently with how other Paper 45 references treat it." Paper 45's existing bibliography entry `latremoliere_metric_st_2017` is the Adv. Math. 2023 published version with arXiv preprint November 2018 — the canonical reference used throughout Paper 45 for the metric-spectral-triple weak-form propinquity construction. All three edits cite this entry consistently with the rest of the paper. No bibliography modifications needed.

The directive's parenthetical "Latrémolière 2017 §5" is preserved in the prose throughout the edits (the §5 is the Latrémolière 2017/2023 paper's §5 — the metric-spectral-triple weak-form propinquity bound and its proof shape).

---

## §3. Compile verification

**Three-pass pdflatex compile:**
- Pass 1: `pdflatex -interaction=nonstopmode paper_45_lorentzian_propinquity.tex` — completes with the initial cross-reference warning (`Label(s) may have changed. Rerun to get cross-references right.`) which is normal first-pass behavior.
- Pass 2: re-run resolves cross-references.
- Pass 3: stable, zero `Rerun` warnings, zero undefined references.

**Final output:** `paper_45_lorentzian_propinquity.pdf`, **20 pages**, 605,348 bytes.

**Page count change:** 19 → 20 pages (one additional page from the new remark + expanded proofs).

**Pre-existing warnings (unchanged from prior state, environmental):**
- 2 `Package hyperref Warning: Token not allowed in a PDF string (Unicode)` — from Latrémolière's é character in bookmark generation (purely cosmetic, doesn't affect the PDF body).
- ~6 `Overfull \hbox` warnings — pre-existing line-breaking issues unrelated to this sprint, including paths in §appendix B (`debug/data/l3b_2_sub_sprint_D.json`) and a few text overruns; all under 90pt and not affecting correctness.
- 1 `Underfull \hbox (badness 10000)` in the van den Dungen bibliography entry — pre-existing.
- The pre-existing `microtype` disable comment at line 42 (MiKTeX font-expansion environmental issue) is retained.

**No new warnings introduced by these edits.**

---

## §4. Anything unexpected

**One LaTeX issue encountered and resolved:** the initial Edit 2 used `\paragraph{Bookkeeping under the $\Kplus$-restricted Lipschitz seminorm.}` as the header sentence, which triggered

```
! LaTeX Error: Something's wrong--perhaps a missing \item.
l.1203 \end{proof}
```

The cause is amsthm's `proof` environment uses an internal `trivlist` structure, and `\paragraph` (which itself opens a structural group) interacts badly with the trivlist's `\item` expectation at `\end{proof}`. Replaced `\paragraph{...}` with `\emph{...}` (the canonical alternative for in-proof structural emphasis); compile cleared immediately. The substantive content of the bookkeeping header was preserved.

No other surprises. The sister L3b-2d memo's three-step bookkeeping pattern (reach_B / reach_P / height_B / height_P with explicit lemma citations) transports verbatim to the K⁺-weak-form setting under the K⁺-restricted Lipschitz seminorm; the only difference is that Paper 45 cites L4 properties on the K⁺-restricted operator system (where the prerequisite K⁺-preservation has to be made explicit), while L3b-2d cites them on the natural-substrate operator system (where K⁺-preservation is trivial because $\{J, M\} = 0$ never engages on chirality-symmetric multipliers).

---

## §5. Status

**COMPLETE.** Three edits applied in-place (with one additional small §1.4 preamble update per directive constraint); three-pass clean compile at 20 pages with zero new errors and zero new undefined references. The Latrémolière 2017 §5 K⁺-preservation bookkeeping is now closed in Paper 45 §5 via explicit operator-level transcription, inheriting the L3a-1 trivial-multiplier K⁺-preservation result of Paper 44 §3–4. No new mathematics introduced.
