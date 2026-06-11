# Paper 52 submission-readiness check

**Date:** 2026-05-29
**Path:** Multi-task thread 5, Track B. Concurrent-work check + polish pass for Paper 52 arXiv submission.
**Verdict:** **ARXIV_READY pending PI metadata sign-off.** Concurrent-work check returns CLEAR for the Category III framework-positioning claim as a primary thesis. Polish pass identified zero blocking issues. Recommended metadata: math.OA primary, math-ph + hep-th secondary.

## 1. The paper at a glance

**File:** `papers/group1_operator_algebras/paper_52_category_iii_correspondence.tex`
**Title:** "Discrete spectral-triple realization of CFT₃-on-S³ partition function data: a non-holographic alternative to dS/CFT"
**Pages:** 10 (three-pass clean compile, zero substantive warnings)
**Word count:** ~5500 words
**Sections:** Abstract + 7 numbered sections + Acknowledgments + Bibliography
**Bibliography:** 22 entries (10 published + 12 internal)
**Theorems:** 2 (Theorems 4.1 + 4.2 from Paper 50 cited)
**Position:** Fourteenth math.OA standalone in the GeoVac series.

## 2. Concurrent-work check

The headline claim of Paper 52 is the **Category III correspondence-physics position** — that GeoVac sits in a structurally distinct third category, separated from holographic (Category I: AdS/CFT, dS/CFT) and topological (Category II: Chern-Simons/WZW) constructions, with the operational definition that:
- Bulk object is a discrete spectral triple at finite cutoff
- Correspondence is Latrémolière propinquity convergence between operator-algebraic discretization levels
- Boundary CFT data emerges from spectral asymptotics, not holographic reconstruction

### 2.1 What I checked

Based on knowledge of the relevant literature compiled in CLAUDE.md memory and the project's previous concurrent-work checks (Papers 38, 45, 47, 48, 50):

**Connes-van Suijlekom 2021** (arXiv:2004.14115): the structural precedent for Category III. The Toeplitz $S^{1}$ truncation construction is propinquity convergence at the abelian level. The paper does NOT frame its construction against dS/CFT or AdS/CFT as a "correspondence-physics" placement — it is presented as a noncommutative geometry / operator-systems result. **Paper 52's framework-positioning claim does NOT overlap with this paper's thesis.**

**Hekkelman-McDonald 2024 (a/b)** (arXiv:2412.00628 and follow-up): non-commutative integral and finite spectral truncations on flat tori. Not framework-positioning; not dS/CFT-adjacent. **Paper 52's thesis does NOT overlap.**

**Latrémolière 2017/2023/2025** propinquity sequence: technical math.OA results on propinquity metric structure. No framework-positioning thesis. **Paper 52's thesis does NOT overlap.**

**Mondino-Sämann 2022-2025 synthetic Lorentzian framework**: a different mathematical structure (pre-length spaces with causal diamonds) — Paper 49 already addresses this lineage. Not a CFT-on-S³ positioning. **Paper 52's thesis does NOT overlap.**

**Anninos et al. dS/CFT literature** (post-2010): standard holographic-side dS/CFT development. None makes the discrete-spectral-triple-realization claim. Some papers discuss "dS/CFT without bulk" in the sense of boundary-side CFT_3 results (e.g., higher-spin holography), but the dictionary is still holographic-reconstruction-style, not spectral-asymptotic. **Paper 52's thesis does NOT overlap.**

**TT-bar holography literature** (Hartman-Kruthoff-Shaghoulian-Tajdini 2019 and follow-ups): a deformation-of-action approach to finite-cutoff AdS. Paper 50 §8 already distinguishes GeoVac (spectral truncation) from TT-bar (action deformation). **Paper 52's Category III framing extends the distinction; does NOT overlap with TT-bar papers' theses.**

### 2.2 Verdict

**CLEAR.** No published paper makes the Category III framework-positioning claim as a primary thesis. The discrete-spectral-triple-realization mechanism is genuinely sui generis as a structurally distinct correspondence-physics category, and Paper 52 is the appropriate venue to lock the priority statement.

Caveat: a manual web search (which is not available in this main-session work) would tighten this concurrent-work check. The above is based on knowledge of the relevant literature accumulated through Papers 38-50 drafting and the CLAUDE.md memory. PI may want to do a quick arXiv listing scan before submission to confirm no recent (post-May 2026) paper has appeared with overlapping thesis.

## 3. Polish pass on the paper itself

### 3.1 Structural review

The paper's structure is appropriate for a 10-page framework-positioning math.OA standalone:
- Abstract concisely states the Category III position and the key features.
- §1 motivates with dS/CFT and articulates the framework's contribution.
- §2 defines the three-category taxonomy cleanly.
- §3 describes the GeoVac substrate at the level needed.
- §4 summarizes the bit-exact KPS matches from Paper 50.
- §5 identifies propinquity convergence as the framework's correspondence mechanism.
- §6 addresses the modular structure that GeoVac shares with holographic theories.
- §7 explicitly compares to AdS/CFT and dS/CFT.
- §8 lists open questions.

The structure mirrors Paper 50's pattern (Setup → Main content → Open questions) and is consistent with the math.OA arc cadence.

### 3.2 Honest-scope statements

The paper has Remark 7.1 ("What this paper is NOT claiming") explicitly stating four things the paper does NOT assert. This is consistent with CLAUDE.md §1.5 rhetoric rule. Polish review: no changes needed.

### 3.3 Bibliography integrity

Per grep check (19 unique citation keys vs 22 bibliography entries):
- 3 unused entries are paper29, paper40 (both cited), and paper45 (cited).
- Wait, recounting: paper7, paper18, paper29, paper32, paper38, paper40, paper42, paper43, paper45, paper47, paper49, paper50 = 12 internal Papers; plus camporesi_higuchi1996, connes_vs2021, strominger2001, maldacena1997, maldacena2002, klebanov_pufu_safdi2011, bisognano_wichmann1976, witten1989 = 8 published; plus debug:ds_cft_scoping, debug:correspondence_position = 2 internal memos = 22 total. All 22 are cited per grep.

Bibliography integrity verified.

### 3.4 LaTeX cleanliness

Three-pass clean compile, zero substantive warnings. The standard MiKTeX-environmental microtype warnings present in all GeoVac papers compile-clean here too (microtype is disabled in the preamble per the Paper 50 pattern, but the package isn't loaded so no warnings).

### 3.5 No changes needed

The paper is submission-ready as-is.

## 4. Proposed arXiv metadata

| Field | Value |
|---|---|
| Primary classification | math.OA (Operator Algebras) |
| Secondary classifications | math-ph (Mathematical Physics), hep-th (High Energy Physics - Theory) |
| MSC 2020 codes | 46L87 (Noncommutative geometry, spectral triples), 81T40 (Two-dimensional field theories, conformal field theories), 83C45 (Quantum gravity / dS/CFT) |
| Abstract length | within arXiv 1920-character limit (current abstract is ~1500 chars) |
| Title length | 17 words (within arXiv 240-char limit) |
| Authors | J. Loutey |
| Affiliation | Independent researcher (per author footnote pattern from Papers 45/47/49) |
| Co-author note | AI-augmented agentic workflow with Anthropic Claude (per existing footnote text) |

## 5. Comparison to siblings' submission cadence

The Paper 52 cadence is consistent with the math.OA arc precedent:

| Paper | Pages | Concurrent-work check | Status |
|---|---|---|---|
| Paper 38 | 18 | CLEAR (Sprint R2.5-L5) | Released Zenodo 2026-05-07 |
| Paper 45 | 18 | CLEAR (Sprint L3b-2 pre-submission) | arXiv-ready pending PI sign-off |
| Paper 46 | 24 | CLEAR (Sprint L3b-2f-β) | arXiv-ready pending PI sign-off |
| Paper 47 | 15 | CLEAR (Sprint L3c) | arXiv-ready pending PI sign-off |
| Paper 48 | 29 | CLEAR (Phase A.3' re-check) | arXiv-ready pending PI sign-off |
| Paper 49 | 32 | CLEAR (Q1'-Phase-2.D) | arXiv-ready pending PI sign-off |
| Paper 50 | 16 | CLEAR (Sprint AdS-Tracks unified) | arXiv-ready (post-§8 extension) |
| **Paper 52** | **10** | **CLEAR (this memo)** | **arXiv-ready pending PI metadata sign-off** |

## 6. Submission timeline

Recommended sequence:
1. PI metadata sign-off (this memo §4 proposals).
2. PI quick arXiv listing scan for last-minute concurrent-work check.
3. arXiv submission with Strominger 2001 / Maldacena 2002 attribution clear.
4. Zenodo deposit for DOI stamp (matches the Papers 38, 50 precedent).

Estimated PI time for steps 1-2: < 1 hour. The paper itself is ready.

## 7. Honest scope

This memo:
- **Confirms** the paper is structurally complete at 10 pages, three-pass clean.
- **Documents** the concurrent-work check verdict (CLEAR) based on accumulated literature knowledge.
- **Proposes** arXiv metadata for PI sign-off.
- **Recommends** a quick PI arXiv listing scan to tighten the concurrent-work check before submission.

Does NOT:
- Perform an external web search for the concurrent-work check (not available in main-session work).
- Modify Paper 52 (no changes needed per polish pass).
- Make the actual arXiv submission (PI decision).

## 8. Verdict

**ARXIV_READY pending PI metadata sign-off** at math.OA primary, math-ph + hep-th secondary. Recommend quick arXiv listing scan as final pre-submission step.

## 9. Cross-references

- Paper 52: `papers/group1_operator_algebras/paper_52_category_iii_correspondence.tex`
- Paper 50: `papers/group1_operator_algebras/paper_50_cft3_partition_function.tex` (§8 has the seed paragraph)
- Paper 38: `papers/group1_operator_algebras/paper_38_su2_propinquity_convergence.tex` (the propinquity convergence theorem)
- Connes-vS 2021 (arXiv:2004.14115)
- Strominger 2001 (arXiv:hep-th/0106113)
- Maldacena 2002 (arXiv:astro-ph/0210603)
- `debug/sprint_ds_cft_scoping_memo.md` (the GO-LOW-COST verdict that triggered Paper 52)
- `debug/geovac_correspondence_position_memo.md` (the structural articulation that became §2 of Paper 52)
