# Corpus-Wide Citation Audit — Master Memo

**Date:** 2026-06-02
**Scope:** Pre-broadcast bibliography survey across the GeoVac paper corpus.
**Companions:** `debug/cite_sweep_A.md` (math.OA Lorentzian — 10 papers, 371 bibitems), `debug/cite_sweep_B.md` (QED + precision — 11 papers, 190 bibitems), `debug/cite_sweep_C.md` (chemistry + foundations + synthesis — 13 papers, 190 bibitems), `debug/wave1_salvage_memo.md` + `debug/wave1_pi_recommendations.md` (Wave 1 results), and `debug/review_paper{0,1,7,14,16,22,23,24,27,31,32,38,39,40,42,43}.md` (per-paper review reports).

## 1. Executive summary

The dominant pre-broadcast risk across the corpus is **bibliographic**, not mathematical. Out of 31 papers checked for citations across Waves 1+2 + the three citation sweeps, the internal mathematics generally survives (Pass A audits found very few E-grade errors). What fails is the bibliography:

- **The Fursaev–Solodukhin failure mode** (real arXiv ID / DOI pointing to a different paper than the bibitem claims) has occurred in at least 9 distinct cite-keys across the corpus. The original incident logged in CLAUDE.md §3 was not isolated; it was the visible tip of an LLM-drafting metadata-fabrication pattern.
- **Cross-paper key collisions** occur (same cite-key resolves to different papers in different bibliographies). Most acute: `hekkelman_mcdonald2024` resolves to the real arXiv:2412.00628 in some papers and a ghost arXiv:2403.18619 in others.
- **Cluster concentration**: the damage is heavily concentrated in the math.OA Lorentzian arc (Papers 44–49 + Group 1 synthesis) — that subfamily was drafted in a single sprint window with shared bad-metadata DNA. Other clusters (Papers 28/33/34, Paper 23 precision-AMO, Paper 14 Trenev, Paper 16 Schwerdtfeger) have isolated misattributions specific to their drafting.
- **One structural issue:** Paper 20's bibliography file `paper_20_refs.bib` is missing from the repository — the paper does not compile as-is.

The mathematics is sound; the corpus's bibliographic infrastructure needs systematic remediation before broadcast.

## 2. Verdict matrix (post-fixes-as-of-this-memo)

| Verdict | Papers |
|:--------|:-------|
| **GREEN** (clean or near-clean) | 11, 12, 15, 17, 19, 29, 30, 35, 41, 50, 51, 52, 53, 54, FCI-Atoms |
| **YELLOW** (minor metadata, no broadcast blockers) | 0, 1, 7, 13, 18, 24, 25, 27, 32, 36, 39, 40, 42, FCI-Molecules, Paper 26, Group 3 synthesis |
| **YELLOW-RED → YELLOW** (post-fix) | 16 (math fix applied), 31 (5 fixes applied), 33 (chamseddine fix applied), 34 (chamseddine + Nieuviarts fixes applied) |
| **RED — autonomous-fix-pending** | 14, 22, 23, 28 (Sweep B may have been wrong — already correct?), 38 (multiple fixes), 44, 45, 46, 47, 48, 49, Group 1 synthesis, **Paper 20** |

Approximately **15 GREEN, 17 YELLOW, 14 RED** post-current-fixes. Pre-fix the RED count was ~22; some closed to YELLOW via the autonomous fixes this session.

## 3. Failure-mode taxonomy

### 3.1 Fursaev–Solodukhin (real ID, wrong paper)

| Cite-key | Wrong ID/title/etc. | Correct | Where found | Status |
|:---------|:--------------------|:--------|:------------|:-------|
| `hekkelman_mcdonald2024` (variant) | arXiv:2403.18619 + "Spectral truncations of $T^d$..." | arXiv:2403.18619 = OpenMP/Floyd-Warshall CS paper. Intended T^d work is **Leimbach–vS arXiv:2302.07877** (already cited as `leimbach_vs2024`) | Papers 38, 39, 43, 44, 45, 46, Group 1 synth (~7 papers) | Open — needs per-paper triage |
| `hekkelman_mcdonald2024` (correct variant) | arXiv:2412.00628 "A noncommutative integral..." | This IS the real Hekkelman-McDonald work | Papers 32, 47, 48, 49 | OK as-is (but key collision with above) |
| `hekkelman2022` | arXiv:2206.13744 | 2206.13744 = Kerr-Melvin BH paper; real Hekkelman 2022 = arXiv:2111.13865 | Papers 38, 39, 44, Group 1 synth | Open |
| `ucp_maps_2024` | arXiv:2410.15454 attributed to Hekkelman-McDonald-vS | 2410.15454 = Bhattacharyya et al. | Paper 38 | Open (CITE-CANT-FIND for the intended work) |
| `zhu_casini2020` | arXiv:2003.00315 / "Zhu, Casini, Dalmonte, Hauke" | Actual authors at 2003.00315 = **Zhang, Calabrese, Dalmonte, Rajabpour** | Paper 42 | Open |
| `bizi_brouder_besnard2018` | "Towards a noncommutative geometry of the Standard Model with neutrino mixing" | Actual = "Space and time dimensions of algebras with applications to Lorentzian NCG and QED" (JMP 59, 062303) | Papers 31, 32 (✓ fixed), 44–46 | Partially fixed |
| `bizi_brouder_besnard2018` (different variant) | "Spectral action in Lorentzian signature" + invented CQG 35 (2018), 175004 | Same actual JMP 59, 062303 paper | Papers 47, 48, 49 | Open — different wrong title than the one already fixed |
| `chamseddine_connes2010` | "The spectral action principle in NCG / + superconnection" + Phys. Rev. D 83, 045001 (2011) | Real arXiv DOI of Phys. Rev. D 83, 045001 (2011) per agent = Czajka-Mrowczynski "Collective Excitations of Supersymmetric Plasma." Gold copy = Fortsch. Phys. 58, 553–600 (2010) | Papers 33 (✓ fixed), 34 (✓ fixed). Paper 28 was a Sweep-B false flag — already correct. Paper 32 has multi-source variant (separate triage) | Mostly closed |
| `mondino_samann2024` | arXiv:2209.14384 attributed to Mondino-Sämann | Real 2209.14384 = Minguzzi-Suhr | Papers 45, 46, Group 1 synth | Open |
| `che_perales_sormani2025` | "R. Che" + wrong title (claimed) + wrong journal (Adv. Math.) | Author = **M. Che** (Mauricio); title = "Gromov's compactness theorem for the intrinsic timed-Hausdorff distance"; journal = DGA 103 | Papers 47, 48, 49 | Open |
| `Pachucki2023` | PRL 130, 023004 = "Three-Photon Exchange" | Actual = "Recoil Corrections" by Pachucki–Yerokhin | Paper 23 | Open |
| `PachuckiYerokhin2010` | PRL 104, 070403 = deuterium HFS | Actual = "Fine Structure of Heliumlike Ions" | Paper 23 | Open |
| `Eides2024` | DOI 10.1016/j.physletb.2024.139049 = Eides hyperfine | Actual = black-hole physics by Hadi-Akbarieh | Paper 23 | Open |
| `trenev2025` | arXiv:2311.03719 = electronic-structure LiH/H₂O Pauli counts | Actual = vibrational acetylene-like polyynes; cited Pauli counts not in this paper | Paper 14 | Open (source unknown) |
| `Schwerdtfeger2015` | Four-author book chapter in *Chemistry of Superheavy Elements* 2nd ed. | Does not appear to exist; intended replacement is Pershina 2014 chapter or Schwerdtfeger NPA 944 (2015) | Paper 16 | Open |

### 3.2 Wrong metadata (right paper, wrong year/vol/pages/author)

| Cite-key | Wrong | Correct | Where found |
|:---------|:------|:--------|:------------|
| `latremoliere2018` | Trans. AMS vol. 370 (2018) | Vol. **368 (2016)** | Papers 44, 45, 46, 47, 48, 49 — **6-paper propagation** |
| `latremoliere_metric_st_2017` | Adv. Math. 415 (2023), 108876, 88pp | Adv. Math. **404 (2022), 108393, 56pp** | Papers 44, 45, 46, 47, 48, 49, Group 1 synth — **7-paper propagation** |
| `leimbach_vs2024` author | "Lukas Leimbach" or "L. Leimbach" | **Malte Leimbach** | Papers 44, 45, 46, Group 1 synth |
| `avery_wen_avery1986` | JMP 27 (1986), 3-author | **Wen & Avery (2 authors), JMP 26 (1985)** | Paper 44 |
| `A.~Nieuviarts` / `B.~Nieuviarts` bibitem author | Various wrong initials | **G. Nieuviarts** (Gaston) | Papers 32 (✓ fixed), 42 (✓ fixed), 43 (✓ fixed), 44 (✓ fixed), 31 (✓ fixed), 34 (✓ fixed in this turn). Wave 2 may surface more. |
| `cacic2009` | Lett. Math. Phys. 100, 181–202 (2012) | LMP **103, 793–816 (2013)**, arXiv:1209.4832 | Paper 32 (✓ fixed) |
| `Pyykkoe2012` | bibitem key vs paper year | Paper is 2011 | Paper 16 cosmetic |
| `FriarPayne2005` | PRA 014501 | **PRC 014002** | Paper 23 |
| `Strohmaier 2006` year-in-key | actual 2000 | — | Paper 42 cosmetic |
| `Peruzzo` article number | 4213 | 5213 | Paper 26 |
| `Eides 2001` | mixes 2007-book title with 2001-journal pagination | — | Paper 36 |

### 3.3 Misbundles (correct papers, wrong joint cite)

`perez_sanchez2024` (real: "Bratteli networks…") bundled with `perez_sanchez2025` (real: "Comment on… Yang–Mills without Higgs") as if both were the "without-Higgs correction." Only 2025 is the correction; 2024 derives Yang–Mills-with-Higgs. Found and fixed in correction-claim contexts across Papers 30, 31, 32, 38, Group 1 synthesis. Lineage/"continued-by" cites correctly retain both 2024 and 2025.

### 3.4 Structural defects

- **Paper 20:** bibliography file `paper_20_refs.bib` is missing from the repository — paper cannot LaTeX-compile. Highest single-paper severity issue in the corpus.
- **XXXXX placeholders shipped to PDF:** `mondino_samann2024/2025` in Paper 47, `minguzzi_suhr2024` and `sakovich_sormani2024` in Papers 48, 49. Unfilled template literals reached publication.

### 3.5 Synthesis-paper inheritance

The synthesis papers in `papers/synthesis/` inherit the bad bibitems of their source papers:
- Group 1 synthesis = Lorentzian-arc cluster's bad bibitems (latremoliere, leimbach, hekkelman_mcdonald, hekkelman, mondino_samann) + a unique `perez_sanchez2024` wrong title.
- Group 3 synthesis = milder; only `perez_sanchez2024` wrong title.

## 4. Fixes applied this session

35 total mechanical edits across this session (yesterday + today):

### Yesterday's salvage (8 fixes)
Paper 0 §V H₂ figure; Paper 7 §VI.B math; Paper 32 (×4) perez_sanchez correction-bundles, bizi_brouder title, cacic2009 metadata.

### Today, content level
Paper 16 Eq.(4) `2(N-2)²` → `2(N-2)(N-1)` + Table I 11-row recompute; Paper 0 abstract 3 softenings.

### Today, citation level
- Paper 31: perez_sanchez bibitem split, bizi_brouder title fix, Nieuviarts ×2, line 200 cite update.
- Paper 30, 38, synthesis Group 1: perez_sanchez correction-bundle cleanup.
- Paper 32 (1 more), 42 (3), 43 (3), 44 (2): A.~Nieuviarts → G.~Nieuviarts.
- Paper 33: chamseddine_connes2010 bibitem (this turn).
- Paper 34: chamseddine_connes2010 bibitem + B.~Nieuviarts → G.~Nieuviarts (this turn).

## 5. Remaining work — what's autonomously applicable

The following fixes are now sufficiently well-characterized that they could be applied autonomously in a single batch with verification reads:

### 5.1 Metadata fixes (mechanical, replace strings)

| Cite-key | Where | Edit |
|:---------|:------|:-----|
| `latremoliere2018` Trans. AMS | Papers 44, 45, 46, 47, 48, 49 | Vol. 370 (2018) → **368 (2016)** |
| `latremoliere_metric_st_2017` Adv. Math. | Papers 44, 45, 46, 47, 48, 49, Group 1 synth | 415 (2023), 108876, 88pp → **404 (2022), 108393, 56pp** |
| `leimbach_vs2024` author | Papers 44, 45, 46, Group 1 synth | "Lukas"/"L." → **Malte / M.** |
| `hekkelman2022` arXiv ID | Papers 44 (and possibly 38 via re-grep) | 2206.13744 → **2111.13865** |
| `avery_wen_avery1986` | Paper 44 | JMP 27 (1986), 3-author → **JMP 26 (1985), Wen & Avery, 2-author** |

### 5.2 Bibitem replacements (replace whole bibitem block)

| Cite-key | Where | Replacement |
|:---------|:------|:-----------|
| `bizi_brouder_besnard2018` (Lorentzian-arc variant with "Spectral action in Lorentzian signature" title) | Papers 47, 48, 49 | Use the gold copy from Paper 32 (already corrected) |
| `mondino_samann2024` (in Papers 45, 46 — arXiv:2209.14384 is actually Minguzzi-Suhr) | Papers 45, 46 | Find Mondino-Sämann's correct arXiv ID for the intended Lorentzian pre-length-space work |
| Group 1 + Group 3 synth `perez_sanchez2024` wrong title | Both synth papers | Replace with correct title "Bratteli networks and the Spectral Action on quivers" |

### 5.3 Things that need PI input first

- **Paper 14 `trenev2025`** — source of the actual Pauli baseline numbers needs to be identified.
- **Paper 16 `Schwerdtfeger2015`** ghost cite — needs replacement target choice.
- **Paper 22 1.44% vs 6.06%** — physics judgment on pair-diagonal-m sparsity convention.
- **Paper 23 three CITE-MISATTRIBUTED** (Pachucki2023, PachuckiYerokhin2010, Eides2024) — need correct cites identified.
- **Paper 18 §III vs §VIII Ω-normalization** — convention choice affects Papers 0, 7.
- **Paper 20 missing `.bib` file** — needs recovery / reconstruction.
- **Paper 38 math error** — `4/π = 2·Vol(S¹)/Vol(SU(2))` identity coefficient; pick correct derivation (coefficient should be 4 to match the M1 Vol(S²)/π² reading).
- **`hekkelman_mcdonald2024` cite-key collision** — per-paper triage on which use-context calls for the real `2412.00628` vs the intended `leimbach_vs2024`.
- **XXXXX placeholders** — find the intended cites.
- **`che_perales_sormani2025`, `ucp_maps_2024`** — full bibitem rewrites with correct authors/titles.

## 6. Recommended next steps

1. **Single mechanical-fix sprint** (autonomous, ~20 edits): apply §5.1 and §5.2 patterns above. Closes most of the Lorentzian-arc REDs to YELLOW or GREEN.
2. **PI-decision pass** on §5.3 items. Most are 1-line fixes once the right cite/value/convention is named.
3. **Paper 20 recovery**: locate the missing `.bib` file, or reconstruct from text-level cite usage.
4. **Verify the previously-applied corrections compile cleanly**: a single LaTeX-compile sanity pass over Papers 31, 32, 33, 34 (and any other this-session edits) before broadcast.
5. **Wave 2 Batch 2 (math.OA Papers 44–46) and Batch 3 (Papers 47–49)** can resume *after* the mechanical-fix sprint — those audits will be substantively about the content (not the bibliography we already mapped) at that point.

## 7. Calibration note

The Sweep B agent miscalibrated on three claims this session:
- Claimed Paper 28 had a misattributed `chamseddine_connes2010` bibitem — Paper 28 actually has the correct gold-copy bibitem already.
- Claimed Papers 2 and 25 cite `perez_sanchez2024` for "without-Higgs" — grep found no `perez_sanchez` cites in either paper.

These were the only material miscalibrations across ~750 bibitems sampled by the three sweeps. The verify-the-verifier discipline (cross-corpus + exact-text checks) caught them at the PM level before any edits propagated. Confidence in the audit instruments remains high; the discipline is what kept Sweep B's miscalibration from contaminating the corpus.
