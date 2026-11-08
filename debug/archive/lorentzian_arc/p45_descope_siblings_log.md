# P45 descope — sibling-paper erratum/descope edits log

**Date:** 2026-06-09. **Task:** apply erratum/descope notices to the five papers downstream of the Paper 45 descope (PI-approved Option C) + caveat remark in Paper 32. Canonical record: `debug/sprint_p45_hardening_phase1_memo.md`; agent memos `debug/p45_adversarial_L2_memo.md`, `debug/p45_adversarial_L5_memo.md`; verification artifact `debug/p45_kplus_seminorm_check.py`.

All papers under `papers/group1_operator_algebras/`. Papers 45 and 38 NOT touched (PM-owned). CLAUDE.md / MEMORY / CHANGELOG not touched. No splinter files.

---

## Paper 40 — `paper_40_unified_propinquity_convergence.tex`

L2-correction treatment (NOT the full template — no Lorentzian content).

| # | Location | Edit |
|---|----------|------|
| 1 | Weight-determination rationale before Lemma `lem:L2` (formerly "choice is determined by the cb-norm requirement") | Rewritten: weight equidistributes Plancherel mass; v1 cb-norm rationale flagged as erroneous (any normalised non-negative central kernel gives a UCP smoothing map, cb-norm 1); properties consumed by the proof — (a) and (d) — unaffected. |
| 2 | Lemma `lem:L2`(b) | Relabelled "Plancherel symbol" → "Plancherel mass distribution"; corrected statement distinguishing mass distribution dim V_π/Z from the true convolution-multiplier symbol m_K (= 1 at trivial rep, decreasing, doubled CG window). |
| 3 | Lemma `lem:L2`(c) | Replaced "cb-norm bound" with "Smoothing map is UCP": cb-norm exactly 1; sup_π dim V_π/Z ≤ 2/(Λ+1) relabelled as the mass-distribution maximum, "not the cb-norm of any map used in the proofs". |
| 4 | Proof parts (b), (c) of `lem:L2` | Rewritten in place: (b) Parseval mass-distribution identity + true symbol with SU(2) counterexample values (1, √2/3, 2/9) at n_max=2; (c) elementary UCP argument; Bożejko–Fendler 1991 noted in place as a discrete-groups theorem that does not apply (bibitem retained, non-orphaned); mass-distribution maximum derivation retained as a statement about the mass distribution. |
| 5 | New `\begin{remark}[Erratum, 2026-06-09]\label{rem:erratum_L2}` after the proof of `lem:L2` | States the corrections, that γ_Λ(G) closed form + universal 4/π asymptote are kernel moments and unaffected, that reach_P is reduced to a named gap (antiderivative-transference, Leimbach–vS; scalar sector corroborated at state-space level by GE–vS), and that the §L5 tunneling-pair vocabulary follows Paper 38, whose framework attribution carries a same-date erratum (device = van Suijlekom state-space GH framework). |
| 6 | Definition `def:berezin_general` | Removed the v1 Peter–Weyl weighted-sum display `eq:B_def_peter_weyl` (the two forms disagree at f = 1); convolution form declared the definition; multiplier-symbol weights on doubled CG window noted. No other `\ref` to the removed equation existed. |
| 7 | Proof of `lem:L5`, dual-reach step | Replaced the invalid "via the cb-norm bound of L2(c)" argument with the named-gap statement (forward cb-norm = 1 controls no partial inverse; Leimbach–vS-style transference needed; GE–vS state-space corroboration for the scalar sector); concluding combination flagged "modulo the reach_P named gap". |
| 8 | GE–vS positioning + bibitem | **Deviation: already present in v1.** Paper 40 already cites `gaudillot_estrada_vs2025` (abstract, §"Concurrent state-space result", bibitem with IMRN 2025 + arXiv:2310.14733) with exactly the requested positioning (state-space GH for all compact metric groups, no Dirac, no rate; Paper 40 contributes explicit rate + universal 4/π). No new bibitem needed; erratum remark and reach_P fix add two further cites of the existing key. |

**Compile:** three-pass clean, 26 pages, 0 errors, no undefined references/citations, no multiply-defined labels.

---

## Paper 46 — `paper_46_strong_form_lorentzian_propinquity.tex`

Full template.

| # | Location | Edit |
|---|----------|------|
| 1 | Abstract, first sentence | One-sentence qualification appended: main strong-form theorems descoped pending repair; the natural-substrate seminorm is degenerate; bound is an evaluation of a closed-form rate formula, not a quantum-metric distance. |
| 2 | New `\subsection*{Erratum and descope notice (2026-06-09)}` after the keywords block, before §1 | Full template: three upstream failures (K⁺ seminorm ≡ 0 with the {J_L, D_GV⊗I} = 0 mechanism, bit-exact at (2,3),(3,5); fabricated Latrémolière Thm 5.5/Def 3.4/3.5 — device is vS state-space GH; 2/(n_max+1) is mass-distribution max, smoothing map UCP cb-norm 1). Descope: Theorem `thm:main` + transports; Λ^strong = Λ^P45 "free upgrade" (`rem:free_upgrade`) explained as two evaluations of the same rate formula on a seminorm with kernel far exceeding the scalars. Surviving: Lemma `lem:temporal_invisibility` as bit-exact algebra reinterpreted as the *diagnosis* of the degeneracy; Appendix `app:enlarged` enlarged-substrate construction as algebra, flagged as the natural starting point for the repair ({J,M} = 0 generators have non-zero commutators). Reopened as named research target. |

**Compile:** three-pass clean, 25 pages, 0 errors, no undefined references/citations.

---

## Paper 47 — `paper_47_two_rate_hybrid_convergence.tex`

Full template.

| # | Location | Edit |
|---|----------|------|
| 1 | Abstract, inner-arrow sentence | Qualified: inner arrow + §`sec:g2_metric` G2-metric closure descoped pending upstream repair; outer arrow explicitly "unaffected by the descope". |
| 2 | New `\subsection*{Erratum and descope notice (2026-06-09)}` after the abstract, before `\tableofcontents` | Full template. Descoped: inner arrow (`thm:inner`, consumes Paper 45's Λ quantity); §7/§8 G2-metric closure (`thm:g2_metric`, Q1, Q2 — same degenerate quantity, additionally quoting the v1 constant 2/(n_max+1) as a rate). Surviving stated explicitly: outer arrow `thm:outer` (norm-resolvent, exponential tail O(e^{-|Im z| T/2}), operator-level, no Lipschitz seminorm) and three-carrier identification `thm:three_carriers` at the norm-resolvent level. |
| 3 | §`sec:g2_metric` head | Bracketed italic descope note (degenerate seminorm; temporal-Lipschitz-invisibility = diagnosis not feature; v1 constant quoted as rate; section retained as record, metric-level claims withdrawn). |
| 4 | §`sec:open` preamble + Q1 + Q2 | Preamble notes the hypertopology upgrade is descoped and Q1/Q2 revert to open; Q1 and Q2 status labels changed from "CLOSED" to "previously claimed CLOSED …; DESCOPED 2026-06-09, reopened". |
| 5 | **Pre-existing compile defects repaired (not introduced by this task; present in HEAD)** | (a) dangling cross-document `\ref{thm:g3_closure}` (a Paper 32 label) replaced by literal "Theorem 7.4" with the label name moved into the existing footnote; (b) missing `paper23` bibitem added (internal house style), fixing the orphan `\cite{paper23}` in Q3. |

**Compile:** three-pass clean, 19 pages, 0 errors, no undefined references/citations.

---

## Paper 48 — `paper_48_krein_ms_bridge.tex`

Full template.

| # | Location | Edit |
|---|----------|------|
| 1 | Abstract, first sentence | One-sentence qualification appended: domain category instantiated on the Paper 45 substrate whose Lipschitz seminorm is degenerate; metric-level theorems descoped pending repair; bridge survives as categorical design. |
| 2 | New `\subsection*{Erratum and descope notice (2026-06-09)}` after the abstract, before `\tableofcontents` | Full template. Descoped: domain category KreinMetaMet_pp as instantiated is degenerate as a metric object (instantiated objects not quantum (proper) metric spaces in the cited frameworks); T3 (`thm:bit_exact_panel`) transports formula values, not metric data; T6 (`thm:g2_wedge_closure`) descoped. Surviving: categorical design (Wick-rotation functor W, Connes–Rovelli thermal-time reading of the F2 mismatch) as a proposal conditional on the option-B repair (Toeplitz-compression temporal multipliers). Reopened as named research target. |
| 3 | §`sec:t3_bit_exact_panel` head | Bracketed italic descope note (formula values, not metric data). |
| 4 | §`sec:t6_g2_wedge_closure` head | Bracketed italic descope note (metric-level claim withdrawn pending repair; Paper 47 norm-resolvent ingredients unaffected). |

**Compile:** three-pass clean, 30 pages, 0 errors, no undefined references/citations.

---

## Paper 49 — `paper_49_oslpls_strong_form_bridge.tex`

Full template.

| # | Location | Edit |
|---|----------|------|
| 1 | Abstract, Main-result sentence | Qualification inserted: properties quantifying over the upstream Paper 45/46 metric quantity descoped pending repair; thermal-time-stack content state-level and survives. |
| 2 | New `\subsection*{Erratum and descope notice (2026-06-09)}` after the abstract, before `\tableofcontents` | Full template. Descoped: bit-exact Λ panel inheritance (§`sec:lambda_panel`) and the K⁺ framing inherited from Papers 45/48 (formula evaluations, not metric measurements); bridge-theorem properties of `thm:bridge_strong` quantifying over the degenerate quantity — in particular B4′ convergence transport (consumes Paper 46 strong-form bound) and the metric-ball reading of B3′. Surviving: cocycle entropy-production (max-divergence) deficit computations of §`sec:uhlmann_panel` as state-level numerics (modular theory + Datta chain inequality, no Lipschitz seminorm); TICI cocycle algebra (`thm:tici`) and strict super-additivity (`thm:strict_super_additivity`) as statements about KMS states; OSLPLS categorical design. Repair-relevant flag: enlarged substrate ({J,M} = 0 chirality-flipping generators) is where seminorms are genuinely non-zero. Reopened as named research target. |
| 3 | §`sec:lambda_panel` head | Bracketed italic descope note (bit-exact agreement = same closed-form formula on both sides). |
| 4 | **Pre-existing compile defects repaired (present in HEAD)** | (a) `\ref{sec:honest_scope}` in the abstract → existing label `sec:honest_scope_intro`; (b) `\ref{sec:strict_strong_form_bridge}` (no such label) → existing label `sec:bridge_theorem` (the aggregate Bridge Theorem section, matching the sentence's intent). |

**Deviation note:** instructions said "the Uhlmann relative-entropy deficit computations" survive; the paper's current revision computes these as **Datta max-divergence** cocycle deficits (the post-M1-fix surrogate, Thm `thm:strict_super_additivity` / §`sec:uhlmann_panel`, same 67/69/81-nat values). Erratum wording uses the paper's current terminology ("cocycle entropy-production (max-divergence) deficit … Datta chain inequality").

**Compile:** three-pass clean, 35 pages, 0 errors, no undefined references/citations.

---

## Paper 32 — `paper_32_spectral_triple.tex`

Caveat-remark treatment (NOT the full template).

| # | Location | Edit |
|---|----------|------|
| 1 | New `\begin{remark}[Status caveat, 2026-06-09]\label{rem:gh_status_caveat}` immediately after the proof sketch of Theorem `thm:gh_convergence` (§VIII), before `rem:gh_limit_id` | Four items per instructions: (i) L2 cb-norm conflation (mass distribution vs multiplier symbol; smoothing map UCP cb-norm 1; rate content — γ closed forms + 4/π — unaffected); (ii) reach_P dual estimate named gap (antiderivative-transference, Leimbach–vS Adv. Math. 439 (2024); scalar sectors on compact metric groups covered by GE–vS, cited); (iii) distance actually proved is van Suijlekom's state-space GH distance — "Latrémolière propinquity" label withdrawn (tunnels from quantum isometries, not UCP pairs; downstream Paper 45 v1 theorem/definition numbers nonexistent in the cited source); (iv) kernel condition holds for engineered offdiag CH Dirac (1/14 at n_max=2, 1/55 at n_max=3, verified) but fails for the truthful chirality-diagonal CH Dirac used in the theorem statement (10/14, 26/55) — truthful↔offdiag bridge named gap. Statement under repair; spatial content expected to survive. |
| 2 | `rem:wh1_proven` (WH1 keystone closure remark) | One sentence appended pointing to `rem:gh_status_caveat` (closure claim qualified; proof chain under repair; spatial content expected to survive in the vS framework). |
| 3 | Bibliography | New bibitem `gaudillot_estrada_vs2025` added (IMRN 2025, rnaf197; arXiv:2310.14733), Paper 32 plain bibitem style, after `vansuijlekom2024_ksystems`. |
| 4 | **Pre-existing compile defects repaired (9 dangling keys, all present in HEAD; none introduced by this task)** | (a) `\cite{suijlekom_book2015}` → existing `vansuijlekom_book2015` (typo'd key); (b) `\cite{loutey_paper28}` ×2 → existing `paper28`; (c) missing internal `paper51` bibitem added (orphan `\cite{paper51}` in the Mixed-Tate paragraph); (d) `\S\ref{sec:case_exhaustion}` → existing label `sec:pi_source_theorem` (the section housing `thm:pi_source_case_exhaustion`); (e) `\ref{sec:tensor_product_verification}` (cross-doc Paper 35 label) → literal "(Paper 35, tensor-product verification section)"; (f) `\ref{thm:eta_trivialization}` ×2 (theorem lives in Paper 18 §IV.6, no such label in Paper 32) → literal "the η-trivialization theorem (Paper 18 §IV.6)"; (g) `\ref{sec:w1e_refinement_f4_f6}`, `\ref{sec:w1e_schmidt_core_correlation}` (cross-doc Paper 19 labels) and `\ref{sec:conv_w1e_cross_domain_wall}` (cross-doc Paper 34 label) → `\S\texttt{...}` literal form, matching the house style already used two lines below for Paper 34's `sec:conv_w1e_f4_f6_refinement`. |
| 5 | Proof sketch of `thm:gh_convergence`, L2/L5 bullets | **Deliberately left as-is** (they repeat the v1 Plancherel-symbol / cb-norm / Bożejko–Fendler language): the instruction for Paper 32 is a caveat remark, not in-place rewriting; the caveat explicitly covers those items. |

**Compile:** three-pass clean, 84 pages, 0 errors, no undefined references/citations (all 9 pre-existing dangling keys resolved), no multiply-defined labels.

---

## Cross-cutting notes

- Erratum text in every paper sticks to what the memo establishes: three verified failures quoted with mechanism, bit-exact verification cells (2,3)/(3,5), vS-framework identification, mass-distribution/cb-norm correction, named-gap language for reach_P, kernel-condition numbers 1/14, 1/55 vs 10/14, 26/55, GE–vS scope (state-space GH, compact metric groups).
- No new orphan `\cite` keys anywhere; the only bibitem additions are `gaudillot_estrada_vs2025` (Paper 32), `paper23` (Paper 47), `paper51` (Paper 32), all immediately cited by pre-existing or new text.
- Log-check method note: per-key grep of wrapped LaTeX warnings can miss long keys (caught `sec:w1e_schmidt_core_correlation` this way); final verification used the wrap-immune aggregate lines ("There were undefined references/citations", "multiply defined") on all six logs — all zero.
- Papers 45 and 38 untouched per task boundary. `paper_29_ramanujan_hopf.tex` working-tree modifications pre-date this task (visible in starting git status) and were not touched.
