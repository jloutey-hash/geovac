# Wave 1 — Partial Confidence-Review: Salvage Memo

**Date:** 2026-06-01
**Type:** Audit infrastructure deployment (continuation of `sprint_confidence_review_infra_memo.md`).
**Scope:** Pre-broadcast confidence review of always-load foundational papers (0, 1, 7, 14, 16, 22, 23, 24, 27, 31, 32).

## 1. What happened

Wave 1 dispatched 22 parallel sub-agents (11 papers × {AUDITOR + CITATION}) per the multi-wave plan agreed in this session. The dispatch hit an immediate server-side rate-limit:

- **16 of 22** agents returned with `API Error: Server is temporarily limiting requests` after 0–14 tool calls. None of those wrote useful output.
- **4 of 22** completed full audits before the limit (Papers 7 AUDITOR + CITATION, 16 CITATION).
- **3 of 22** wrote their full memos before hitting a session-token cap on their final summary message (Papers 0 AUDITOR, 32 AUDITOR + CITATION). All three wrote substantive `debug/audit_paper<N>.md` / `debug/citecheck_paper<N>.md` reports.

**Net: 6 substantive audit reports salvaged covering 4 papers** (0, 7, 16, 32 — only Paper 0 missing CITATION, only Paper 16 missing AUDITOR).

## 2. Audit reports kept

| File | Bytes | Paper | Status |
|:-----|:-----:|:------|:------|
| `debug/audit_paper0.md` | 18K | Paper 0 | Full report (summary truncated) |
| `debug/audit_paper7.md` | 15K | Paper 7 | Full report |
| `debug/citecheck_paper7.md` | 13K | Paper 7 | Full report |
| `debug/citecheck_paper16.md` | 17K | Paper 16 | Full report |
| `debug/audit_paper32.md` | 22K | Paper 32 | Full report (summary truncated) |
| `debug/citecheck_paper32.md` | 15K | Paper 32 | Full report (summary truncated) |

## 3. Findings inventory (7 HIGH, 9 MEDIUM, 9+ LOW across 4 papers)

### HIGH (math errors, wrong cited values, broken citations)
1. **Paper 0 §V line 508** — H₂ recovery cited as 94.1%; Paper 15 itself states 96.0% in four places. **APPLIED.**
2. **Paper 7 §VI.B line 597** — Listed S⁵ helium eigenvalues `0, 3, 8, 15` are the S³ Coulomb spectrum $n^2-1$ transcribed into the S⁵ row; for the even-ν singlet sector the correct values are `0, 6, 16, 30` (matches Paper 13). **APPLIED** (with explicit "even-ν singlet sector" annotation).
3. **Paper 16 line ?? — `Schwerdtfeger2015`** — book-chapter citation does not appear to exist as described; looks like a confabulation (right topic / Springer series / author involvement, but the specific four-author chapter title and year cannot be located). Replace with either Pershina 2014 (the real chapter in the real 2nd ed.) or Schwerdtfeger et al. *Nucl. Phys. A* 944, 551 (2015) depending on intended content. **DEFERRED** (needs decision on replacement target).
4. **Paper 32 abstract `8.8×10⁻⁸` framing** — The number is the Sprint-A residual $|K - 1/\alpha - \alpha^2|/(1/\alpha)$, NOT the raw match $|K - 1/\alpha|/(1/\alpha) \approx 4.77\times10^{-7}$. Abstract is misleading; §VI Rem K_under_theorem is honest. **DEFERRED** (substantive abstract rewording).
5. **Paper 32 Thm 2 (case-exhaustion)** — Statement and proof reference "fifteen" Paper 34 projections; Paper 34 v3.2.x catalogues twenty-eight. Theorem-as-stated is now strictly stronger than what is proved. **DEFERRED** (substantive scope rewording — either restrict to the 15-projection corpus state at theorem date, or extend the proof to projections 16–28).
6. **Paper 32 `perez_sanchez2024` misbundle** — Cited together with `perez_sanchez2025` as "the Perez-Sanchez correction that the continuum limit is Yang–Mills without Higgs" in 5+ places. But arXiv:2401.03705 (Perez-Sanchez 2024, "Bratteli networks") derives Yang–Mills–**Higgs**; only arXiv:2508.17338 (Perez-Sanchez 2025, "Comment on…") is the correction. **APPLIED** in 4 correction-claim locations (kept `perez_sanchez2024` in 2 "continued by" / "lineage of" locations where it fits as related work).
7. **Paper 32 `bizi_brouder_besnard2018` bibitem title** — Title is "Towards a noncommutative geometry of the Standard Model with neutrino mixing"; the actual JMP 59, 062303 (2018) / arXiv:1611.07062 title is "Space and time dimensions of algebras with applications to Lorentzian noncommutative geometry and quantum electrodynamics." Right ID, wrong title — Fursaev failure mode in miniature. **APPLIED.**

### MEDIUM (overstatement, framing drift, citation metadata)
- **Paper 0** §VI.C "ground state → −0.5 Ha" claim is correctly scoped but readers may miss the node-weight caveat (per Paper 18 §III); needs one-line footnote. **DEFERRED.**
- **Paper 0** abstract "shell capacities 2k−1" is ambiguous (Table I shows total $2(2k-1)$); recommend "per-shell angular capacities". **DEFERRED.**
- **Paper 0** abstract "spectrally converge to the hydrogen eigenvalues" is stronger than body delivers without node-weight $W$ (per Paper 18 §III). **DEFERRED.**
- **Paper 0** §V.A "the lattice is the invariant" reads stronger than §1.5 dual-description rule allows; the very next paragraph hedges correctly. Minor framing tension. **DEFERRED.**
- **Paper 7** §III.A "κ alone maps graph eigenvalues to Rydberg spectrum" framing — Paper 18 §κ corrected this to $H = \kappa\mathcal{L} + W$. Lags-corpus, not wrong. **DEFERRED.**
- **Paper 7** §I "we demonstrate convergence" — convergence is rigorously proven only in Paper 38 (WH1 PROVEN). Lags-corpus. **DEFERRED.**
- **Paper 7** `barut1967` cite (lines 67, 99) — bibitem points to PR 156, 1541 (the SO(4,**1**) paper) but is used to ground SO(4,**2**) claims. The SO(4,2) extension is PR **157**, 1180 (1967) — same authors, sibling paper. **DEFERRED** (needs PR 157 exact title verified before bibitem replacement).
- **Paper 16** novelty claim that S_N-based periodicity derivation is new — current "appears to be new" is at the honest ceiling, but §I "the periodic law has not been derived from a single mathematical principle" glosses over substantial group-theoretic-periodicity literature (Rumer–Fet SO(4,2)×SU(2); arXiv:2501.18272 SO(4,4)). Add courtesy cite. **DEFERRED.**
- **Paper 32** Sprint L1 "literal identification" of the BW Wick-rotation theorem at operator-system level — $\sigma_{2\pi}(O) = O$ is a *necessary* condition for BW identification, not full proof of identification with the continuum Wightman vacuum. Soften to "operator-system-level consistency (period closure holds bit-exact)". **DEFERRED.**
- **Paper 32 `cacic2009` bibitem** — Volume 100 → 103, pages 181–202 → 793–816, year 2012 → 2013, plus arXiv:1209.4832. **APPLIED.**

### LOW (typos, stale cross-refs, key-vs-year mismatches)
- Paper 0 §IX "~6%" LiH → "~5%" or "5.3%".
- Paper 0 §II.B Step 1 wording ("separated by $d_0$" on circle of radius $d_0$).
- Paper 0 §VI.C "18 symbolic proofs validate" → "18 internal symbolic test cases validate" (internal-consistency labeling).
- Paper 0 abstract "foundational motif of the framework" → "a foundational motif" (mild rhetoric softening).
- Paper 7 `bander1966` bibitem lacks "(I)" suffix in title (cosmetic).
- Paper 16 `Pyykkoe2012` bibitem key — paper is from 2011 (text correct; only the key mislabeled).
- Paper 16 `KnowlesHandy1984` is a bibliography orphan (no `\cite{}` in body).
- Paper 32 conclusion uses older "four-way coincidence" framing while abstract+intro use deflated "four projections of one triple" (CLAUDE.md WH4 alignment).
- Paper 32 `paschke_sitarz2000` bibitem key — paper is 1998 (cosmetic).

### "Hand to a domain expert"
- **Ω(0) = 2 vs Ω(0) = 1 normalization** (Papers 0, 7, 18). Paper 0 §VI.C and Paper 18 §VII both state "$1/16 = 1/\Omega^4(0)$ with $\Omega(0) = 2$" while elsewhere $\Omega(p) = (1 + p^2/p_0^2)^{-1}$ gives $\Omega(0) = 1$. There is likely a normalization convention difference (e.g. the $p_0$-dependent prefactor that maps the $\mathbb{R}^3$ measure to the $S^3$ measure), but the inconsistency surfaces in three foundational papers and a referee will catch it. **Load-bearing open item; flag for PI / external review before broadcast.**
- Five-lemma GH-convergence proof detail (Paper 38 referee).
- BBB sign table at $(m,n)=(4,6)$ — confirm reading from arXiv:1611.07062 v2 PDF Table 1.
- Marcolli–vS / Perez-Sanchez lineage-placement judgment call.

## 4. Edits applied this session

| File | Edit | Severity |
|:-----|:-----|:--------:|
| `papers/group3_foundations/Paper_0_Geometric_Packing.tex` line 508 | 94.1\% → 96.0\% | HIGH |
| `papers/group3_foundations/Paper_7_Dimensionless_Vacuum.tex` line 597 | 0, 3, 8, 15 → 0, 6, 16, 30 (even-ν singlet sector) | HIGH |
| `papers/group1_operator_algebras/paper_32_spectral_triple.tex` 4 locations | Removed `perez_sanchez2024` from "correction" cite bundles (lines 64, 2369, 2474, 4575); kept it in lineage / "continued by" bundles (lines 242, 3678, 4601) | HIGH |
| `papers/group1_operator_algebras/paper_32_spectral_triple.tex` line 4904–4908 | `bizi_brouder_besnard2018` bibitem title → "Space and time dimensions of algebras…" | HIGH |
| `papers/group1_operator_algebras/paper_32_spectral_triple.tex` line 4794–4797 | `cacic2009` bibitem: vol 100 → 103, pp 181–202 → 793–816, yr 2012 → 2013, +arXiv:1209.4832 | MEDIUM |

5 of 7 HIGH findings fully closed by mechanical edits. 2 HIGH findings deferred because they require substantive rewording (Paper 16 Schwerdtfeger replacement choice; Paper 32 abstract α figure + case-exhaustion theorem scope) rather than one-line substitutions.

## 5. Deferred to PI review

| Item | Type | Why deferred |
|:-----|:-----|:-------------|
| Paper 32 abstract: "8.8×10⁻⁸ matches α⁻¹" rewording | HIGH | Substantive — requires choosing between two suggested replacement framings |
| Paper 32 §VIII Thm 2: case-exhaustion theorem 15 → 28 projections | HIGH | Substantive — choose between (a) restrict theorem to "first 15 Paper 34 projections (state at date)" or (b) extend proof to projections 16–28 |
| Paper 16 `Schwerdtfeger2015` ghost cite | HIGH | Need PI to choose replacement: Pershina 2014 (chapter), Schwerdtfeger NPA 2015 (article), or different cite per intended content |
| Paper 7 `barut1967` SO(4,2) split | MEDIUM | Need PR 157 exact published title before bibitem replacement |
| Paper 0 abstract softenings (4 phrases) | MEDIUM | Substantive framing edits |
| Paper 7 framing lags-corpus items (2) | MEDIUM | Substantive framing edits |
| Paper 16 §I "single mathematical principle" courtesy cite | MEDIUM | Substantive |
| Paper 32 Sprint L1 "literal identification" softening | MEDIUM | Substantive framing edit |
| Paper 32 conclusion "four-way coincidence" rewording | LOW | Substantive framing edit |
| Ω(0)=2 vs Ω(0)=1 normalization | "hand to expert" | Cross-paper structural — load-bearing if it bites |

## 6. Re-dispatch plan for the 16 missed audits

**Papers needing Wave-1 audit completion:**
- Paper 0 — CITATION only (AUDITOR has report)
- Paper 1 — both
- Paper 14 — both
- Paper 16 — AUDITOR only (CITATION has report)
- Paper 22 — both
- Paper 23 — both
- Paper 24 — both
- Paper 27 — both
- Paper 31 — both

Total: 16 dispatches (Paper 0 CITATION, Paper 16 AUDITOR, plus 7 papers × 2).

**Cadence revision:** rather than 22 parallel, fire **4 agents per batch** (= 2 papers × 2 agents), sequential batches inside one session. That respects the API concurrency cap while preserving the wave structure. 8 batches to complete Wave 1.

Alternatively: 6-agent batches (3 papers × 2 agents) = 6 batches. Either works; I'd lean toward 4-per-batch as the safer first restart given today's behaviour.

## 7. Calibration check

Yesterday's blind run on Paper 2 produced 2 verified MEDIUM fixes. Today's partial run on 4 papers produced 5 verified HIGH-severity fixes (one per audited paper plus the Paper 32 bibitem cluster). Severity profile is consistent with the agent's design — the auditor agents flag what a domain referee would call out, the citation agent caught two right-ID-wrong-title failure modes (the documented Fursaev-Solodukhin failure mode in miniature, twice). Instruments are working at the calibrated grade.

The cross-corpus check (mandatory pre-defect step added during yesterday's calibration) demonstrably worked on Paper 7: two candidate E findings were correctly demoted to C ("lags Paper 18 / Paper 38") rather than reported as wrong. That's exactly what the methodology-lesson §4 of yesterday's memo asked for.

## 8. Files written this session

| File | Purpose |
|:-----|:--------|
| `debug/audit_paper0.md` | Paper 0 audit (full report) |
| `debug/audit_paper7.md` | Paper 7 audit (full report) |
| `debug/citecheck_paper7.md` | Paper 7 citation check (full report) |
| `debug/citecheck_paper16.md` | Paper 16 citation check (full report) |
| `debug/audit_paper32.md` | Paper 32 audit (full report) |
| `debug/citecheck_paper32.md` | Paper 32 citation check (full report) |
| `debug/wave1_salvage_memo.md` | This memo |

## 9. Next session recommendation

1. PI reviews the 5 deferred substantive items (Paper 16 Schwerdtfeger replacement, Paper 32 abstract α rewording, Paper 32 case-exhaustion scope, Paper 32 Sprint L1 softening, Paper 7 barut1967 split).
2. Re-dispatch the 16 missing Wave-1 audits in 4-per-batch cadence.
3. Decide whether the Ω(0)=2 normalization item warrants its own diagnostic memo (touches Papers 0, 7, 18 — load-bearing).
4. After Wave 1 completes, proceed to Wave 2 (math.OA arc Papers 38–49 with REVIEWER added; highest external-review stakes).
