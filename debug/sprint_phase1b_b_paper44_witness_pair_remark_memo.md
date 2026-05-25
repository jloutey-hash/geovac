# Sprint Phase 1B-b — Paper 44 witness-pair remark per Connes–vS §IV

**Date:** 2026-05-24
**Sprint:** Phase 1B-b (1C diagnostic follow-on)
**Status:** COMPLETE
**Deliverable:** New Remark `rem:cvs_sharpening` added to Paper 44 §5.3 immediately after `rem:comparison`, sharpening the witness-pair structural reading per the Riemannian Connes–vS §IV C*-envelope vs operator-system distinction (Paper 32 §III `rem:operator_system`).

---

## §1. What was added

Single new `\begin{remark}[Connes--vS \S IV sharpening:\ envelope-relative informative content]\label{rem:cvs_sharpening}` block (~330 words) inserted between `rem:comparison` and `\section{Trivial Krein-positive restriction}` in `papers/group1_operator_algebras/paper_44_lorentzian_operator_system.tex`.

No other content modified. Pre-existing Remarks `rem:op_system`, `rem:both_real` (the §5.3 subsection-as-mechanism-explanation), `rem:comparison`, `rem:state_level` untouched.

## §2. Structural reading established

The Remark articulates a four-step argument:

1. **Recall the Riemannian template.** Connes–vS §IV reads the operator-system data as the surviving algebraic content vs the "non-informative full matrix algebra" $C^*$-envelope. On the Toeplitz $C(S^1)^{(n)}$ and on Paper 32's scalar Fock-projected $S^3$ Hilbert space, the witness pair (closure failure) and the propagation-2 result are two facets of one envelope: $\dim(\Op_{n_{\max}}) < \dim(\mathrm{env}) = N^2$ strictly, but $\dim(\Op^2) = N^2$.

2. **Achievable-envelope reading transports verbatim.** Under $\mathcal{V}_\achievable = \mathcal{V}_\Weyl \otimes \mathcal{D}_{N_t}$, the Lorentzian witness pair exhibits closure failure at 38.12% (n_max=2, every $N_t$ tested) and 35.89% (n_max=3) per the L3a-1 JSON data, while $\dim((\Op^L)^2) = \dim_\achievable$ exactly — Paper 32's Connes–vS reading transports cleanly.

3. **Full-envelope reading exhibits a second, deeper structural obstruction.** Under the full $\Bcal(\Krein)$ envelope, the witness-pair residual is supplemented by a structural ceiling: chirality-block-diagonal closure + commutative temporal subalgebra together cannot generate chirality-flipping or non-diagonal-temporal directions of $\Bcal(\Krein)$ at any finite $k$ (Theorem `thm:prop_envelope_dependent` eq.~\eqref{eq:prop_full_inf}). This is structural, not a closure failure of the truncation.

4. **Envelope-relative informative content.** On the achievable envelope, Paper 44 preserves Paper 32's structural reading exactly. On the full envelope, the operator system is "even more sub-informative" because it hits a Lorentzian-geometric ceiling. The achievable envelope is therefore named as the natural target for the Connes–vS propagation question on the Lorentzian extension; the full-envelope $\prop = \infty$ is recorded as a structural ceiling, not a failure of the achievable-envelope reading.

## §3. Cross-references resolved

Cross-references established:
- `Paper~32~\S~III` `rem:operator_system` (`paper32` citekey) — opening recall of the Riemannian Connes–vS reading.
- `Paper~32 Definition~3.2` — scalar Fock-projected Hilbert space cited.
- `Theorem~\ref{thm:witness_pair}` — the witness-pair theorem this Remark refers to (residual values).
- `Theorem~\ref{thm:prop_envelope_dependent}` — both `eq:prop_ach` and `eq:prop_full_inf` referenced.
- `debug/data/l3a_1_lorentzian_operator_system.json` — numerical residual values cited (38.12% and 35.89% values).
- `\cite{connes_vs2021}` and `\cite[\S\,3.1]{connes_vs2021}` — both already in Paper 44 bib.
- New label `cvs_sharpening` registered in `paper_44_lorentzian_operator_system.aux` (verified one occurrence).

## §4. Verification

Three-pass `pdflatex` compile:
- PASS 1 exit 0
- PASS 2 exit 0
- PASS 3 exit 0
- Output: 18 pages (unchanged from pre-edit page count; the new Remark fits within the existing §5 layout)
- Zero undefined references
- Zero new LaTeX warnings
- Pre-existing hyperref unicode "Token not allowed in PDF string" warnings preserved (cosmetic, environmental, not from this edit)

## §5. What this Remark does NOT do

- Does NOT introduce a new theorem, proposition, definition, or computation.
- Does NOT rewrite any existing Remark or section.
- Does NOT modify the Connes axiom audit, the envelope dichotomy theorem, or the witness pair theorem.
- Does NOT change Paper 44's page count or section structure.
- Does NOT introduce any new bibliography entries.

The Remark is a pure framing-and-cross-reference enhancement that makes the structural parallel with Paper 32 §III explicit and names the achievable envelope as the natural Connes–vS target on the Lorentzian extension.

## §6. Final status

**COMPLETE.** Paper 44 §5 now contains the sharpened witness-pair-via-Connes–vS-§IV reading required by the 1C diagnostic. The achievable vs full envelope dichotomy is now framed inside the Connes–vS informative-vs-non-informative reading, with both the Riemannian-template parallel and the genuinely-new Lorentzian sharpening made explicit.

No unexpected issues. Paper 44 remains 18 pages, three-pass clean.
