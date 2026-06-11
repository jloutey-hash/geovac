# Sprint Q1'-Phase-2.D — Paper 49 drafting

**Date:** 2026-05-25 (Phase-2.D drafting; immediate follow-on to Phase-2.B.3 numerical panel verification 2026-05-25).

**Sprint position:** Q1'-Phase-2.D of the staged sprint structure (Phase-1 → Phase-2.A → {Phase-2.B.1, Phase-2.B.2} → Phase-2.B.3 → **Phase-2.D**). Final stage of the Q1' arc. Drafts Paper 49 as the 11th math.OA standalone in the GeoVac series.

**Predecessors (load-bearing):**
- `debug/sprint_q1prime_light_diagnostic_memo.md` — Q1'-Light diagnostic
- `debug/sprint_q1prime_phase1_case_a_stepping_stone_memo.md` — Phase-1 Case A stepping stone
- `debug/sprint_q1prime_phase2a_oslpls_category_memo.md` — Phase-2.A OSLPLS category
- `debug/sprint_q1prime_phase2b1_wflip_morphisms_memo.md` — Phase-2.B.1 W^flip on morphisms
- `debug/sprint_q1prime_phase2b2_bridge_theorem_memo.md` — Phase-2.B.2 Bridge Theorem closure
- `debug/sprint_q1prime_phase2b3_panel_verification_memo.md` — Phase-2.B.3 numerical panel
- `papers/group1_operator_algebras/paper_48_krein_ms_bridge.tex` — Paper 48 (direct predecessor; conventions, macros, structure template)
- `papers/group1_operator_algebras/paper_46_strong_form_lorentzian_propinquity.tex` — Paper 46 Appendix B enlarged substrate

**Status:** PAPER DRAFTED. Single multi-section LaTeX file at arXiv-ready status pending PI metadata sign-off.

---

## Headline verdict

**DONE arXiv-ready.** Paper 49 drafted as a single multi-section math.OA standalone integrating the six Q1' phase memos into a unified 32-page manuscript. Three-pass clean LaTeX compilation, zero undefined references, zero errors. 57 bibliography entries integrating Paper 48 lineage + Connes 1973 + Bratteli-Robinson + Uhlmann + Lindblad + Wilde + the GeoVac math.OA series (Papers 24, 29, 32, 38, 39, 40, 42, 43, 44, 45, 46, 47, 48).

The 11th math.OA standalone in the GeoVac series is now available at:
- **Source:** `papers/group1_operator_algebras/paper_49_oslpls_strong_form_bridge.tex`
- **Output:** `papers/group1_operator_algebras/paper_49_oslpls_strong_form_bridge.pdf`

---

## Drafting log

### Section structure (final)

§1 Introduction (~9 pp): Q1' question, three structural moves, twin paradox as quantum information headline, honest scope, roadmap, related work.

§2 Setup (~5 pp): Paper 48 K+-weak-form bridge recap, Paper 46 Appendix B enlarged substrate, Q1'-Light bimodal-Abelianness finding, Mondino-Sämann definitions, Connes-Rovelli + Connes 1973 cocycle Radon-Nikodym derivative (C1)-(C4).

§3 The OSLPLS category (~6 pp): Theorem 3.1 (axiom transport, 8 MS axioms with classification), Definition 3.5 OSLPLS object 5-tuple (A1)-(A5), Definition 3.7 OSLPLS morphism UCP map with (M1)-(M5), Lemma 3.8 modular intertwining automaticity, Theorem 3.10 commutative MS as faithful sub-category via embedding functor ι.

§4 The Krein-side bridge candidate W^flip (~5 pp): Phase-1 chirality-graded stepping stone recap, Definition 4.2 W^flip on objects, Theorem 4.4 Krein-side bridge candidate verifies OSLPLS axioms, Definition 4.5 source morphism class, Definition 4.6 W^flip on morphisms via Connes-vS contravariant UCP convention, Theorem 4.7 functoriality, Theorem 4.8 commutative-case restriction agreement with Paper 48 W functor.

§5 Off-orbit structural activation at full M≠0 (~2 pp): why Phase-1 (M=0) fails to activate Decomposition O Case (iii), why Phase-2 (full M≠0) activates it, three-orbit triples and the off-orbit case, twin-paradox preview remark.

§6 The Connes-Rovelli thermal-time stack across distinct KMS states (~8 pp, **substantive heart**): Lemma 6.1 orbit-KMS lemma, §6.2 Connes cocycle Radon-Nikodym intertwiner, Theorem 6.3 triple-intersection cocycle identity (TICI), §6.4 cocycle entropy production and strict super-additivity (Theorem 6.5), §6.5 physical reading as twin paradox as quantum information, §6.6 honest scope of the thermal-time stack construction with substrate-composition acknowledgment.

§7 The aggregate Bridge Theorem 6.4'-Q1' (~4 pp): Theorem 7.1 B1' structural correspondence, Theorem 7.3 B2' off-orbit super-additivity, Theorem 7.4 B3' pre-compactness inheritance, Theorem 7.5 B4' convergence transport, Theorem 7.6 aggregate Bridge Theorem.

§8 Numerical verification panel (~3 pp): Table 8.1 Lambda joint propinquity bit-exact at $\{(2,3),(3,5),(4,7)\}$, Table 8.2 Uhlmann monotonicity deficits (substantively positive 66-81 nats), Table 8.3 Riemannian-limit recovery bit-exact at $\Nt=1$, §8.4 propagation number $\le 4$ verified at smallest cell. Honest scope.

§9 Synthetic-Lorentzian implications (~3 pp): strict super-additivity as quantum information-theoretic twin paradox (substantive novelty); OSLPLS as candidate non-commutative Lorentzian pre-length space (with the caveat that OSLPLS contains MS rather than directly generalizes MS); positioning vs Mondino-Sämann-Sormani-Cavalletti synthetic-Lorentzian-GH community.

§10 Open questions (~4 pp): Q2' non-commutative MS extension (multi-month NCG-research target); closed-form cocycle entropy production deficit; higher-genus orbit topology / TICI extension; cross-KMS thermal time across modular intertwiners; G3 cross-manifold (blocked at NCG-framework level by Paper 24 §V Coulomb/HO asymmetry); G4 inner-factor calibration data (W3 question of Paper 32; orthogonal to OSLPLS bridge).

§11 Conclusion (~2 pp): summary of three structural moves; headline result of Bridge Theorem 6.4'-Q1' + twin-paradox-as-quantum-information; honest scope + place in math.OA series.

### LaTeX/macro design choices

**Macros locked.** Reused Paper 48's macro set (`\nmax, \Nt, \sthree, \Krein, \Acal, \Mcal, \Ucal, \DGV, \DL, \JL, \BWvac, \Kalpha, \Wfun, \KreinMM, \LorPLG`, etc.) verbatim for full Krein-side and operator-system continuity. Added new macros for OSLPLS-specific structure:

- `\OSLPLS` for the target category $\mathbf{OSLPLS}_{\mathrm{cov}}$
- `\KreinMMflip` for source category at full M≠0
- `\KreinMMflipMzero` for source category at M=0
- `\LorPLGZtwo` for Phase-1 chirality-graded LPLS target
- `\ellOS` for $\ell^{\mathrm{OS}}$, `\tauOS` for $\tau^{\mathrm{OS}}$
- `\Mflipfull` for $\mathcal{M}^{L,\mathrm{flip}}_{\mathrm{full}}$
- `\MflipMzero` for $\mathcal{M}^{L,\mathrm{flip}}_{M=0}$
- `\Wflip` for the strong-form bridge functor
- `\WflipMzero` for the Phase-1 chirality-graded bridge functor
- `\Mflip{N}{L}{M}` for chirality-flip generator, `\Mspat{N}{L}{M}` for spatial multiplier
- `\osXkrein` for the source Krein PPQMS object
- `\Xfrak` for the OSLPLS abstract object
- `\Morb{n}`, `\Oorb{n}`, `\omegaOrb{n}` for orbit-indexed sub-operator-systems / orbits / effective KMS states
- `\cocyc{i}{j}` for the Connes cocycle Radon-Nikodym derivative
- `\DelS{i}{j}` for the cocycle entropy production

**Preamble conventions.** Matched Paper 48 verbatim: `\documentclass[11pt,a4paper]{article}`, UTF-8 unicode chars for Greek letters (sigma, pi, alpha, beta, etc.), `microtype` disabled (MiKTeX environmental issue), `\usepackage[numbers,sort&compress]{natbib}`, amsthm theorem environments.

### Bibliography

57 bibitems integrating:
- **5 published references new to Paper 49** (not in Paper 48): Bratteli-Robinson 1981 (load-bearing for Theorem 5.3.10 KMS structure), Connes 1973 (load-bearing for cocycle Radon-Nikodym derivative), Uhlmann 1977 (load-bearing for relative-entropy monotonicity), Lindblad 1975 (load-bearing for data-processing inequality), Wilde 2017 (modern textbook account); Tomita 1967, Takesaki 1970 (modular theory foundations); Latrémolière 2026 spectral-C1 + Martinetti 2026 adjacent NCG-causal (Q1' concurrent-work re-check).
- **Paper 48 lineage published references** retained verbatim (Allen-Burtscher, Bisognano-Wichmann, Bizi-Brouder-Besnard, Camporesi-Higuchi, Chamseddine-Connes, Che-Perales-Sormani, Connes 1994, Connes-Rovelli, Connes-vS, Farsi-Latrémolière, Franco-Eckstein, Hartle-Hawking, Hekkelman-McDonald, Ketterer, Kubota, Kunzinger-Sämann, Latrémolière lineage 5 entries, Leimbach-vS, Marcolli-vS, Minguzzi-Suhr, Mondino-Ryborz-Sämann, Mondino-Sämann pointed, Müller, Nieuviarts, Sakovich-Sormani, Sewell, Sormani-Vega, Strohmaier, Unruh, van den Dungen).
- **Internal GeoVac preprints**: Papers 24, 29, 32, 38, 39, 40, 42, 43, 44, 45, 46, 47, 48 (13 entries).

### Compile issues caught and fixed

1. **Double-subscript on `\Mflipfull_{1}`**: pdflatex fatal error. The macro `\Mflipfull` expands to `\mathcal{M}^{L,\mathrm{flip}}_{\mathrm{full}}` which already has a subscript. Adding `_{1}` after produces a double subscript. **Fixed** by spelling out `\Mcal^{L,\mathrm{flip}}_{\mathrm{full},1}` in the source morphism class definition (§4.5 Def 4.5 (S2)).

2. **Double-subscript on `\BWvac_{2}`**: similar issue (`\BWvac` expands to `\omega_{W}^{L}` which has a subscript `W`). **Fixed** by spelling out `\omega_{W,2}^{L}` in §4.5 Def 4.5 (S3).

3. **No `\end{X>` typos found** (the Paper 48 cleanup learning checked first). Clean.

After both fixes, three-pass clean LaTeX compile with zero errors and zero undefined references.

### Final metrics

| Metric | Value |
|:-------|:------|
| Source file | `papers/group1_operator_algebras/paper_49_oslpls_strong_form_bridge.tex` |
| Source lines | 2,771 |
| Source words | 12,539 |
| PDF output | `papers/group1_operator_algebras/paper_49_oslpls_strong_form_bridge.pdf` |
| PDF size | 766 KB |
| Page count | **32 pages** |
| Bibliography entries | **57 bibitems** (44 published + 13 internal GeoVac preprints) |
| Compile passes | 3 (cross-references resolved on pass 2) |
| Undefined references | 0 |
| Compile errors | 0 |
| Environmental warnings | newunicodechar redefining (harmless), float specifier 'h' → 'ht' (harmless) |

---

## Structural surprises during drafting

### None of substantive significance

The drafting was load-bearing-input-driven. The six Q1' phase memos collectively totaling ~40,000 words provided the complete substrate for Paper 49 at theorem-grade rigor; the drafting work was synthesis + reformatting + integration with the math.OA-paper structural convention rather than producing new content.

**Minor drafting decisions:**
1. **Twin-paradox-as-quantum-information framing in §1.3 and §9.1**: the strict super-additivity from the data-processing inequality is the headline novelty of Paper 49 beyond Paper 48. It's also a substantive structural observation —- to our knowledge, the first explicit identification of quantum relative-entropy monotonicity as the load-bearing ingredient for the reverse triangle inequality of synthetic Lorentzian geometry. This framing is given prominent placement in the abstract, §1.3 introduction, §6.5 thermal-time stack section, and §9.1 synthetic-Lorentzian implications.
2. **Honest scope on OSLPLS ⊃ MS vs OSLPLS = NC-MS** (§9.2 + §10.1): OSLPLS is the category that CONTAINS MS as a faithful sub-category via the embedding functor ι; it is NOT a direct synthetic extension of the MS concept of pre-length space. The genuine non-commutative MS pre-length space concept (Q2') remains the multi-month NCG-research frontier. This distinction is honored carefully in §9.2 (where OSLPLS is called a "candidate" non-commutative Lorentzian pre-length space) and §10.1 (where Q2' is named as the open question).
3. **Closed-form cocycle entropy production deficit open question** (§10.2): Theorem 6.5 establishes strict super-additivity qualitatively via the data-processing inequality but does not give a closed-form expression for the deficit $\Delta S^{\sigma_1 \to \sigma_2} + \Delta S^{\sigma_2 \to \sigma_3} - \Delta S^{\sigma_1 \to \sigma_3}$. The numerical panel verifies positivity but not closed form. This is named as the principal open question for Paper 49 follow-on work.

### Compression-pattern observation

Paper 49 drafting compressed to a single session at sprint cadence. The compression is consistent with the math.OA-paper substrate-inheritance pattern (Papers 38, 39, 40, 42, 43, 44, 45, 46, 47, 48 all compressed to single drafting sessions following theorem-grade memos). The Q1' arc as a whole — from Q1'-Light diagnostic (2026-05-24 morning) through Phase-2.D Paper 49 drafting (2026-05-25) — closed in approximately 24 hours of compute time across 7 sprint sub-memos and 1 paper draft. This is the structural-skeleton-scope compression pattern that the GeoVac math.OA series has exhibited throughout.

The 11th math.OA standalone caps a sustained ~3-week (May 2026) arc producing 4 papers (48, then 46/47 closures, then 49 as Q1' close-out). This is the natural rate at which the math.OA frontier has been advancing.

---

## Recommended next steps

### Priority 1: PI metadata sign-off

Paper 49 is arXiv-ready pending PI sign-off on metadata:
- Primary classification: math.OA
- Secondary classifications: math-ph, gr-qc
- arXiv abstract (use Paper 49 LaTeX abstract section verbatim)
- Submission date target: at PI's discretion
- License: as per Paper 48 convention

### Priority 2: arXiv submission

After PI sign-off, submit to arXiv via standard procedure. The Zenodo deposit can follow on the same cycle as previous math.OA standalones (Paper 38 released 2026-05-07).

### Priority 3 (optional): CLAUDE.md and MEMORY.md updates

CLAUDE.md §6 paper inventory should add the Paper 49 entry with title, status (drafted, arXiv-ready pending PI sign-off), and key result summary. Similar update to MEMORY.md cross-reference if appropriate.

The §2 Current Development Frontier may want a brief one-line entry per the §13.11 content discipline rule:
```
**Sprint Q1' arc closed (2026-05-25):** Paper 49 drafted as 11th math.OA standalone; strict-strong-form Krein-MS bridge via OSLPLS category + Connes-Rovelli thermal-time stack across distinct KMS states + twin-paradox-as-quantum-information. See debug/sprint_q1prime_phase2d_paper49_drafting_memo.md.
```

### Priority 4: Q2' multi-month NCG-research follow-on

The genuine non-commutative Mondino-Sämann pre-length space concept (Q2', §10.1 of Paper 49) remains the multi-month NCG-research frontier. This is the substantive direction if the PI elects to pursue the synthetic extension beyond OSLPLS. Estimated 6-12 months at sprint cadence; comparable to Latrémolière propinquity development.

Alternative direction: focus on closing the cocycle entropy production deficit closed form (Paper 49 §10.2 open question), which would give quantitative numerical predictions on the panel cells and enable direct comparison with the classical twin-paradox proper-time deficit formula. Sprint-scale 2-4 weeks at standard cadence.

---

## Cross-references

- **Paper 49**: `papers/group1_operator_algebras/paper_49_oslpls_strong_form_bridge.tex`
- **Q1' arc memos** (6): `debug/sprint_q1prime_*_memo.md`
- **Paper 48** (direct predecessor): `papers/group1_operator_algebras/paper_48_krein_ms_bridge.tex`
- **Paper 46 Appendix B** (enlarged substrate): `papers/group1_operator_algebras/paper_46_strong_form_lorentzian_propinquity.tex`
- **GeoVac math.OA series** (siblings): Papers 29, 32, 38, 39, 40, 42, 43, 44, 45, 46, 47, 48 (12 prior; Paper 49 is the 13th preprint but 11th math.OA standalone — Papers 29 and 32 sit at other tier/sub-classification levels in CLAUDE.md §6)

**End of memo.**
