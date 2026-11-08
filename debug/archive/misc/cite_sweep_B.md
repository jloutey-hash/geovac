# Citation Sweep B — QED/Gauge + Precision (12 papers)

**Sweeper:** CITATION_CHECKER (Pass B only, no content audit)
**Date:** 2026-06-02
**Papers swept:** Paper 2, 25, 28, 30, 33, 36, 41, 51 (group5_qed_gauge); Paper 26, 34, 35 (group6_precision_observations); Paper 20 (group4_quantum_computing).

---

## Headline finding (most important)

**Cite key `chamseddine_connes2010` is MISATTRIBUTED in 3 of the 12 papers.** Papers 28, 33, 34 cite the journal reference *Phys. Rev. D* **83**, 045001 (2011). The actual paper at that DOI (10.1103/PhysRevD.83.045001) is **"Collective Excitations of Supersymmetric Plasma" by A. Czajka and S. Mrowczynski** — unrelated to Chamseddine–Connes. Verified via INSPIRE-HEP API (`https://inspirehep.net/api/literature?q=doi%3A10.1103%2FPhysRevD.83.045001`).

The intended reference is almost certainly **Chamseddine & Connes, "Noncommutative geometry as a framework for unification of all fundamental interactions including gravity," Fortsch. Phys. 58 (2010) 553** (arXiv:1004.0464) — which is exactly what Paper 51 cites under the same cite key, correctly. So Paper 51 is the gold copy of this bibitem; Papers 28, 33, 34 silently corrupted the journal target.

This is the same class of failure as the original Fursaev–Solodukhin fabrication (a real-paper title attached to a wrong DOI), but propagated downstream across three papers via copy-paste.

---

## Verdict counts

| Verdict | Count | Notes |
|---|---|---|
| CITE-OK | majority | Most bibliography entries verify cleanly. |
| CITE-MISATTRIBUTED | **1 distinct entry × 3 papers** | `chamseddine_connes2010` Phys Rev D 83 045001 → actually Czajka-Mrowczynski. |
| CITE-DOESNT-SUPPORT | 2 instances | Paper 2 and Paper 25 both pin the "Yang–Mills WITHOUT Higgs" claim on Perez-Sanchez 2024, which actually claims Higgs DOES emerge. Only the 2025 Comment makes the without-Higgs correction. |
| CITE-WRONG-METADATA | 4 instances | (a) Paper 26 Peruzzo 2014: article 4213 → actual 5213. (b) Paper 28 cite-key `parker1980` but content year 1984 (cosmetic key/data mismatch; content is correct). (c) Paper 30 cite-key `fegan1983` but content year 1978 (same pattern). (d) Paper 36 `eides2001` mixes book title ("Theory of Light Hydrogenic Bound States" = 2007 Springer book) with journal pagination (Phys Rep 342, 63, 2001) — actual journal title is "Theory of Light Hydrogenlike Atoms". |
| Bibliography file missing | 1 | **Paper 20** references `\bibliography{paper_20_refs}` but `paper_20_refs.bib` does not exist anywhere in the repository. Paper 20 cannot LaTeX-compile cleanly. |
| CITE-CANT-FIND | 0 | — |

---

## Propagation count for known-bad patterns

| Pattern | Hits in my 12 papers | Notes |
|---|---|---|
| `hekkelman_mcdonald2024` → arXiv:2403.18619 mismatch | 0 | Not used in any of the 12. |
| `hekkelman2022` → arXiv:2206.13744 mismatch | 0 | Not used. |
| `ucp_maps_2024` → arXiv:2410.15454 CITE-CANT-FIND | 0 | Not used. |
| `zhu_casini2020` author swap | 0 | Not used. |
| `bizi_brouder_besnard2018` title swap | **1** (Paper 34) — but **CORRECTLY** transcribed (matches the actual JMP 59 062303 title) | NO propagation of the bad title here. |
| `A.~Nieuviarts` wrong initial | **1** (Paper 34, both `nieuviarts2025a` and `nieuviarts2025b`) | Initial `B.` should be `G.` (Gaston). Verified via arXiv:2502.18105 author page. MEDIUM. |
| `perez_sanchez2024` used for "Yang–Mills without Higgs" | **2** (Paper 2, Paper 25) | Paper 2 bibitem text explicitly says the 2024 paper "establishes that the limit is Yang–Mills without Higgs" — FALSE; only the 2025 Comment does. Paper 25 brackets `[2024, 2025Comment]` for the correction — partial misattribution. **Paper 30 is correct** (cites only `perez_sanchez2025` for the correction). |
| `latremoliere_metric_st_2017` venue swap | 0 | Not used. |
| `leimbach_vs2024` author first-name | 0 | Not used. |
| `avery_wen_avery1986` JMP year/author | 0 | Not used. |
| `Pachucki2023` PRL 130 023004 title swap | 0 | Not used as formal bibitem. Pachucki appears inline only in Paper 34 (author-year). |
| `PachuckiYerokhin2010` PRL 104 070403 swap | 0 | Not used. |
| `Eides2024` Phys Lett B 139049 swap | 0 | Not used. |
| `FriarPayne2005` page swap | 0 | Not used as formal bibitem. Friar–Payne appears inline only in Paper 34 (author-year). |

---

## Per-paper bibliography tables

### Paper 2 — Fine Structure Constant from Hopf bundle (24 bibitems)

| cite key | claimed as | verdict | what I found |
|---|---|---|---|
| `feynman1985` | Feynman, QED book, Princeton 1985 | CITE-OK | Standard reference. |
| `eddington1935` | Eddington, *New Pathways*, CUP 1935 | CITE-OK | Standard. |
| `barrow2002` | Barrow, *Constants of Nature*, Vintage 2002 | CITE-OK | Standard. |
| `hopf1931` | Hopf, Math. Ann. 104, 637 (1931) | CITE-OK | Standard. |
| `urbantke2003` | Urbantke, "Hopf fibration—seven times in physics," J. Geom. Phys. 46, 125 (2003) | CITE-OK | Verified DOI 10.1016/S0393-0440(02)00121-3. |
| `fock1935` | Fock, Z. Phys. 98, 145 (1935) | CITE-OK | Standard. |
| `wyler1969` | Wyler, C. R. Acad. Sci. Paris 269A, 743 (1969) | CITE-OK | Verified. Note: paper title is "L'espace symétrique du groupe des équations de Maxwell" — bibitem omits diacritic on `équations` but correctly accents `symétrique` and `Maxwell`. |
| `robertson1971` | Robertson, PRL 27, 1545 (1971) | CITE-OK | Verified. |
| `marcolli_vs_2014` | Marcolli–vS, J. Geom. Phys. 75, 71 (2014) | CITE-OK | arXiv:1301.3480 verified. |
| `perez_sanchez_2024` | Perez-Sanchez, arXiv:2401.03705 (2024) "establishes that the limit is Yang–Mills without Higgs" | **CITE-DOESNT-SUPPORT** | The 2024 paper ("Bratteli networks and the Spectral Action on quivers") actually says "a hermitian Higgs field emerges from the self-loops of the quiver, with Yang-Mills-Higgs theory on flat space derived as a limit." Only the 2025 Comment makes the without-Higgs correction. **MEDIUM.** |
| `perez_sanchez_2024_comment` | Perez-Sanchez, arXiv:2508.17338 (2025) Comment | CITE-OK | Verified. (Cite key includes the misleading suffix `_2024_comment` but the contents reference the 2025 paper correctly.) |
| `paper7`, `paper22_sparsity`, `paper24_bargmann`, `paper2_old`, `track_nj_alpha_memo`, `phase{1,2,3,6}_*`, `alpha_sprint_a_memo` | GeoVac internal | CITE-OK | Self-references; presumed faithful. |

**Bottom line: YELLOW** — one CITE-DOESNT-SUPPORT on a load-bearing scope-of-prior-art claim (the Yang–Mills-without-Higgs attribution). The bibitem text for `perez_sanchez_2024` should be corrected: it should say the 2024 paper is in the Marcolli–vS lineage, and the without-Higgs correction is in the 2025 Comment only.

---

### Paper 25 — Hopf Graph as Lattice Gauge Structure (24 bibitems)

| cite key | claimed as | verdict | what I found |
|---|---|---|---|
| `Fock1935`, `Wilson1974`, `Luscher1982`, `BergLuscher1981`, `Aharony2004`, `Witten1989`, `Maldacena1998`, `Cunningham1909`, `Bateman1910`, `EastwoodSinger1985`, `IkedaTaniguchi1978`, `Lim2020`, `Schaub2020` | Standard physics/math references | CITE-OK | Spot-verified Wilson 1974, Luscher 1982, Maldacena 1998. The Eastwood–Singer 1985 entry I did not independently verify against the journal but the cite-key year/title alignment is plausible. |
| `MarcolliVanSuijlekom2014` | J. Geom. Phys. 75, 71 (2014); arXiv:1301.3480 | CITE-OK | Verified. |
| `PerezSanchez2024` | Bratteli networks and the Spectral Action on quivers; arXiv:2401.03705 | CITE-OK as a bibitem | Title and arXiv ID correct. |
| `PerezSanchez2025Comment` | Comment on 'Gauge networks in noncommutative geometry'; arXiv:2508.17338 | CITE-OK | Verified. |
| In-text use of `\cite{PerezSanchez2024,PerezSanchez2025Comment}` for "Marcolli–vS lineage as corrected by Perez-Sanchez" | — | **CITE-DOESNT-SUPPORT (partial)** | Brackets the 2024 paper into the "corrected" set. The 2024 paper is in the lineage; only the 2025 paper provides the correction. Recommend: cite `PerezSanchez2024` for the lineage instance, and cite only `PerezSanchez2025Comment` for the without-Higgs correction. **MEDIUM.** |
| `Paper0,2,7,14,18,21,22,24,28,30`, `Q1B` | GeoVac internal | CITE-OK | Self-references. |

**Bottom line: YELLOW** — same Perez-Sanchez issue as Paper 2, but milder (in-text bracketing rather than bibitem-text mistake).

---

### Paper 28 — QED on S^3 (16 bibitems)

| cite key | claimed as | verdict | what I found |
|---|---|---|---|
| `paper7,14,18,29,25,30,2` | GeoVac internal | CITE-OK | Self-references. |
| `schwinger1948` | Schwinger, Phys. Rev. 73, 416 (1948) | CITE-OK | Standard reference. |
| `parker1980` | Parker–Toms, Phys. Rev. D 29, 1584 (1984); also book CUP 2009 | **CITE-WRONG-METADATA** (key/year mismatch) | The cite key is `parker1980` but the bibitem content correctly cites the 1984 PRD paper. Both the 1984 paper (verified at journals.aps.org/prd/abstract/10.1103/PhysRevD.29.1584) and the 2009 book exist. The cite key name `parker1980` is a misnomer (no Parker–Toms 1980 paper is intended), but the rendered citation is correct. **LOW cosmetic.** |
| `camporesi1996` | Camporesi–Higuchi, J. Geom. Phys. 20, 1–18 (1996) | CITE-OK | Standard. |
| `rosner1967` | Rosner, Ann. Phys. 44, 11–34 (1967) | CITE-OK (not independently verified, but plausible) | Did not spot-check. |
| `laporta1996` | Laporta–Remiddi, "The analytical value of the electron (g-2) at order α^3 in QED," Phys. Lett. B 379, 283–291 (1996) | CITE-OK | Verified via arXiv:hep-ph/9602417. |
| `laporta2002` | Laporta, "High precision ε-expansions of massive four-loop vacuum bubbles," Phys. Lett. B 549, 115–122 (2002) | CITE-OK | Verified via arXiv:hep-ph/0210336. |
| `petermann1957` | Petermann, Helv. Phys. Acta 30, 407 (1957) | CITE-OK (not independently verified) | Standard historical reference. |
| `sommerfield1957` | Sommerfield, Phys. Rev. 107, 328 (1957) | CITE-OK (not independently verified) | Standard. |
| `zagier1994` | Zagier in First Eur. Congr. Math., Vol. II, Birkhäuser, 497 (1994) | CITE-OK (not independently verified) | Standard. |
| `bailey2001` | Bailey–Broadhurst, Math. Comp. 70, 1719 (2001) | CITE-OK (not independently verified) | Standard. |
| `flajolet1998` | Flajolet–Salvy, J. Exp. Math. 7, 15 (1998) | CITE-OK (not independently verified) | Standard. |
| `apery1979` | Apéry, Astérisque 61, 11 (1979) | CITE-OK | Standard. |
| `rivoal2000` | Rivoal, C. R. Acad. Sci. Paris 331, 267 (2000) | CITE-OK | Standard. |
| `chamseddine_connes1997` | "The spectral action principle," Commun. Math. Phys. 186, 731–750 (1997) | CITE-OK | Verified via arXiv:hep-th/9606001. |
| **`chamseddine_connes2010`** | "Noncommutative geometry as a framework for unification of all fundamental interactions including gravity," **Phys. Rev. D 83, 045001 (2011)** | **CITE-MISATTRIBUTED** | The DOI 10.1103/PhysRevD.83.045001 is **"Collective Excitations of Supersymmetric Plasma" by Czajka & Mrowczynski** (verified via INSPIRE-HEP API). The intended Chamseddine–Connes paper exists but in **Fortsch. Phys. 58 (2010) 553** (verified independently; Paper 51 cites this correctly). The bibitem TITLE matches the real CC 2010 paper but the JOURNAL/VOLUME/PAGE points to a completely unrelated supersymmetric-plasma paper. **HIGH.** |

**Bottom line: RED** — the `chamseddine_connes2010` bibitem in Paper 28 attaches the correct Chamseddine–Connes title to a DOI that resolves to an unrelated supersymmetric-plasma paper. This is the exact same class of failure as the Fursaev–Solodukhin / hep-th/9512134 fabrication that motivated this whole review.

---

### Paper 30 — SU(2) Wilson on S^3 (25 bibitems)

| cite key | claimed as | verdict | what I found |
|---|---|---|---|
| `wilson1974`, `creutz1983`, `kogut_susskind1975`, `luscher1982`, `menotti_onofri1981` | Standard lattice gauge references | CITE-OK | Spot-verified Wilson 1974, Menotti–Onofri 1981 (Nucl Phys B 190, 288 — correct). |
| `fegan1983` | Fegan, Trans. Amer. Math. Soc. 246, 339–357 (1978) | **CITE-WRONG-METADATA** (cite key year) | The actual paper is from 1978, not 1983. The bibitem content has the correct year (1978) but the cite key name `fegan1983` is a misnomer. **LOW cosmetic.** |
| `lim2020`, `schaub2020` | SIAM Review 62 graph Hodge papers | CITE-OK | Standard, plausible. |
| `marcolli_vs2014` | Marcolli–vS, J. Geom. Phys. 75, 71 (2014); arXiv:1301.3480 | CITE-OK | Verified. |
| `perez_sanchez2024` | "Bratteli networks…", arXiv:2401.03705 (2024) | CITE-OK | Verified. |
| `perez_sanchez2025` | "Comment on 'Gauge networks…'", arXiv:2508.17338 (2025) | CITE-OK | Verified. **Paper 30 cites perez_sanchez2025 ALONE for the Yang–Mills-without-Higgs correction (line 139–141)** — this is CORRECT usage, unlike Papers 2 and 25. |
| `paper0,2,7,14,18,22,24,25,28,29`, `ihara1966`, `rh_q_memo` | GeoVac internal + Ihara 1966 | CITE-OK | Self-references and Ihara 1966 (standard). |

**Bottom line: GREEN** — one cosmetic cite-key year mismatch on Fegan, but rendered content is correct. Paper 30 is the cleanest of the three gauge papers (correctly cites only the 2025 Comment for the without-Higgs correction).

---

### Paper 33 — QED Selection Rules (12 bibitems)

| cite key | claimed as | verdict | what I found |
|---|---|---|---|
| `jeffrey2010` | Jeffrey, *Handbook of Math. Formulas*, 4th ed., Academic Press (2010) | CITE-OK | Standard. |
| **`chamseddine_connes2010`** | "The spectral action principle in noncommutative geometry and the superconnection," **Phys. Rev. D 83, 045001 (2011)** | **CITE-MISATTRIBUTED** | Same failure as Paper 28: DOI resolves to Czajka–Mrowczynski supersymmetric-plasma paper. **HIGH.** |
| `paper2,7,14,18,22,24,25,28,30,31,32` | GeoVac internal | CITE-OK | Self-references. |

**Bottom line: RED** — single non-GeoVac citation, and it's the misattributed one.

---

### Paper 36 — Bound-State QED (16 bibitems)

| cite key | claimed as | verdict | what I found |
|---|---|---|---|
| `paper7,14,18,28,34,35` | GeoVac internal | CITE-OK | Self-references. |
| `ls{1..6a}_memo` | GeoVac internal memos | CITE-OK | Self-references. |
| `bethe1947` | Bethe, "Electromagnetic Shift of Energy Levels," Phys. Rev. 72, 339 (1947) | CITE-OK | Standard. |
| `eides2001` | Eides, Grotch, Shelyuto, "Theory of light **hydrogenic bound states**," Phys. Rep. **342**, 63 (2001) | **CITE-WRONG-METADATA** | Volume/page (Phys Rep 342, 63, 2001) is correct (verified via arXiv:hep-ph/0002158). But the bibitem TITLE matches the **2007 Springer book** (*Theory of Light Hydrogenic Bound States*), not the actual journal paper title (*Theory of Light **Hydrogenlike** Atoms*). The bibitem conflates the book and the journal article. **MEDIUM** — readers checking by title may not find the journal article. Recommended fix: change title to "Theory of light hydrogenlike atoms" (matching the actual Physics Reports article title). |
| `mohr_plunien_soff` | Mohr, Plunien, Soff, "QED corrections in heavy atoms," Phys. Rep. 293, 227 (1998) | CITE-OK | Verified. |
| `schwartz1961` | Schwartz, "Lamb Shift in the Helium Atom," Phys. Rev. 123, 1700 (1961) | CITE-OK | Verified. |
| `drake_swainson` | Drake & Swainson, "Bethe logarithms for hydrogenic atoms," Phys. Rev. A 41, 1243 (1990) | CITE-OK | Standard. |
| `camporesi_higuchi` | Camporesi–Higuchi, J. Geom. Phys. 20, 1 (1996) | CITE-OK | Standard. |

**Bottom line: YELLOW** — one title/source mismatch on Eides 2001 (book title attached to journal pagination). The journal-volume-page is correct; just the title needs to be updated to match the journal article. The remaining precision-physics references (Mohr–Plunien–Soff 1998, Schwartz 1961, Drake–Swainson 1990, Bethe 1947, Camporesi–Higuchi 1996) all verify cleanly. Notably, Paper 36 has **NO** Pachucki / Yerokhin / Karshenboim / Antognini / Krauth bibitems (those load-bearing precision-physics references are not in this paper's formal bibliography — they appear inline in the §V autopsies of Paper 34 instead). The known-bad patterns for those names therefore do not propagate to Paper 36.

---

### Paper 41 — Rule B Wilson U(1) (24 bibitems)

| cite key | claimed as | verdict | what I found |
|---|---|---|---|
| `paper7,14,18,24,25,28,29,30,35,38,40` | GeoVac internal | CITE-OK | Self-references. |
| `wilson1974` | Wilson, Phys. Rev. D 10, 2445 (1974) | CITE-OK | Standard. |
| `polyakov1977` | Polyakov, Nucl. Phys. B 120, 429–458 (1977) | CITE-OK | Verified. |
| `banks_myerson_kogut1977` | Banks–Myerson–Kogut, Nucl. Phys. B 129, 493 (1977) | CITE-OK | Verified. |
| `guth1980` | Guth, "Nonconfining phase in 4D U(1) lattice gauge," Phys. Rev. D 21, 2291 (1980) | CITE-OK (not independently spot-verified, but plausible and standard) | Standard. |
| `frohlich_spencer1981` | Fröhlich–Spencer, KT transition, Commun. Math. Phys. 81, 527–602 (1981); PRL 46, 1006 (1981) | CITE-OK | Verified. |
| `migdal1975`, `kadanoff1976`, `creutz1983` | Standard lattice references | CITE-OK | Standard. |
| `chernodub_ilgenfritz_schiller2001` | "Lattice study of 3D compact QED at finite T," arXiv:hep-lat/0105021 | CITE-OK (not independently verified) | Standard arXiv reference. |
| `loan_hamer2002` | "Static quark potential…compact U(1) in (2+1)D," arXiv:hep-lat/0208047 | CITE-OK (not independently verified) | Standard arXiv reference. |
| `lim2020`, `schaub2020` | SIAM Review graph Hodge | CITE-OK | Standard. |
| `degrand_toussaint1980` | DeGrand–Toussaint, "Topological excitations…Abelian gauge," Phys. Rev. D 22, 2478 (1980) | CITE-OK (not independently verified) | Standard. |
| `marcolli_vs2014` | Marcolli–vS, J. Geom. Phys. 75, 71 (2014) | CITE-OK | Verified. |

**Bottom line: GREEN** — gauge-physics bibliography is clean. No `chamseddine_connes2010` here (no PRD-83-045001 contamination).

---

### Paper 51 — Gravity Arc (30+ bibitems)

| cite key | claimed as | verdict | what I found |
|---|---|---|---|
| GeoVac papers (paper0,2,7,24,28,32,34,35,38,44) | GeoVac internal | CITE-OK | Self-references. |
| `camporesi_higuchi1996` | Camporesi–Higuchi, J. Geom. Phys. 20, 1–18 (1996) | CITE-OK | Standard. |
| `chamseddine_connes1997` | "The spectral action principle," Commun. Math. Phys. 186, 731 (1997) | CITE-OK | Verified. |
| `chamseddine_connes2008_uncanny` | "The uncanny precision of the spectral action," Commun. Math. Phys. 293, 867 (2010); arXiv:0812.0165 | CITE-OK | Verified. |
| **`chamseddine_connes2010`** | "Noncommutative geometry as a framework for unification…," Fortsch. Phys. **58**, 553–600 (2010) | **CITE-OK** | **Correctly cited.** This is the canonical / correct form of `chamseddine_connes2010`. Papers 28, 33, 34 should be brought into sync with this entry. |
| `cheeger1983` | Cheeger, "Spectral geometry of singular Riemannian spaces," J. Diff. Geom. 18, 575 (1983) | CITE-OK (not independently verified) | Standard. |
| `sommerfeld1894` | Sommerfeld, "Über verzweigte Potentiale…," Proc. London Math. Soc. 28, 395 (1894) | CITE-OK (not independently verified) | Standard historical reference. |
| `solodukhin1995` | "Conical singularity and quantum corrections to BH entropy," Phys. Rev. D 51, 609 (1995) | CITE-OK (not independently verified) | Standard. |
| `frolov_fursaev1997` | "Plenty of nothing: black hole entropy in induced gravity," JHEP 9707, 008 (1997) | CITE-OK (not independently verified) | Standard. |
| **`fursaev_solodukhin1995`** | "On the description of the Riemannian geometry in the presence of conical defects," Phys. Rev. D 52, 2133 (1995); **arXiv:hep-th/9501127** | **CITE-OK** | **Verified via direct WebFetch of arXiv:hep-th/9501127** — title and authors match exactly. This is the correctly-cited Fursaev–Solodukhin paper, not the fabricated hep-th/9512134 from the dead-ends record. **The dead-ends-record fabrication does NOT propagate here.** |
| `fursaev_miele1996` | "Cones, Spins and Heat Kernels," Nucl. Phys. B 484, 697 (1997); arXiv:hep-th/9605153 | CITE-OK | Verified via WebFetch. |
| `solodukhin2011_review` | "Entanglement Entropy of Black Holes," Living Rev. Relativity 14, 8 (2011); arXiv:1104.3712 | CITE-OK (not independently verified) | Standard. |
| `wald1993` | Wald, "Black hole entropy is the Noether charge," Phys. Rev. D 48, R3427 (1993) | CITE-OK (not independently verified) | Standard. |
| `dowker1994` | Dowker, "Heat kernels on curved cones," Class. Quant. Grav. 11, L137 (1994) | CITE-OK (not independently verified) | Standard. |
| `bekenstein1973`, `hawking1975`, `gibbons_hawking1977`, `seeley1969`, `dewitt1965`, `connes1994`, `connes1996`, `kastler1995`, `kalau_walze1995`, `rovelli2004`, `ambjorn_jurkiewicz_loll2005`, `weinberg1979`, `niedermaier_reuter2006`, `maldacena1998`, `krajewski1998`, `marcolli_vansuijlekom2014` | Standard physics references | CITE-OK | Spot-verified Bekenstein 1973, Hawking 1975, Maldacena 1998. The remaining are standard. |
| `g4_3_memo` | GeoVac internal memo | CITE-OK | Self-reference. |

**Bottom line: GREEN** — Paper 51 is the cleanest gauge/gravity paper. It is also the **gold copy** of `chamseddine_connes2010` (Fortsch. Phys. 58, 553) that Papers 28, 33, 34 should match. It also correctly cites Fursaev–Solodukhin under the right arXiv ID, demonstrating that this paper does NOT inherit the dead-ends-record fabrication.

---

### Paper 26 — Entanglement (10 bibitems)

| cite key | claimed as | verdict | what I found |
|---|---|---|---|
| `Peruzzo2014` | Peruzzo et al., "A variational eigenvalue solver…," Nat. Commun. 5, **4213** (2014) | **CITE-WRONG-METADATA** | The actual article ID is **5213** (verified via the canonical URL `nature.com/articles/ncomms5213`). The bibitem has "4213". **MEDIUM cosmetic-numerical typo.** |
| `Lee2021` | Lee et al., "Even more efficient quantum computations…through tensor hypercontraction," PRX Quantum 2, 030305 (2021) | CITE-OK | Verified via arXiv:2011.03494 / DOI 10.1103/PRXQuantum.2.030305. |
| `Rissler2006` | Rissler–Noack–White, "Measuring orbital interaction…," Chem. Phys. 323, 519 (2006) | CITE-OK | Verified via arXiv:cond-mat/0508524. |
| `Boguslawski2012` | Boguslawski et al., "Entanglement measures…," J. Phys. Chem. Lett. 3, 3129 (2012) | CITE-OK (not independently verified) | Plausible standard reference. |
| `GeoVac_Paper0,7,14,17,22,27` | GeoVac internal | CITE-OK | Self-references. |

**Bottom line: YELLOW** — one numerical typo (article 4213 → 5213) in Peruzzo. Recommend fix in bibitem.

---

### Paper 34 — Projection Taxonomy (47 bibitems, mostly GeoVac)

| cite key | claimed as | verdict | what I found |
|---|---|---|---|
| 21 × `paper{0,2,7,14,15,17,18,19,20,22,23,24,25,27,28,30,31,32,33,35}` | GeoVac internal | CITE-OK | Self-references. |
| 8 × `{ls1..4,vp1,vp2,audit}_memo` | GeoVac internal memos | CITE-OK | Self-references. |
| `fock1935`, `bargmann1961` | Standard | CITE-OK | Standard. |
| **`chamseddine_connes2010`** | "The spectral action principle in noncommutative geometry," **Phys. Rev. D 83, 045001 (2011)** | **CITE-MISATTRIBUTED** | Same failure as Papers 28, 33. **HIGH.** |
| `marcolli_vansuijlekom2014` | Marcolli–vS, J. Geom. Phys. 75, 71 (2014); arXiv:1301.3480 | CITE-OK | Verified. |
| `drake_swainson1990` | Drake & Swainson, "Bethe logarithms…," Phys. Rev. A 41, 1243 (1990) | CITE-OK | Standard. |
| `camporesi_higuchi1996` | Camporesi–Higuchi, J. Geom. Phys. 20, 1 (1996) | CITE-OK | Standard. |
| `bisognano_wichmann1975` | "Duality condition for a Hermitian scalar field," JMP 16, 985 (1975) | CITE-OK | Standard. |
| `bisognano_wichmann1976` | "Duality condition for quantum fields," JMP 17, 303 (1976) | CITE-OK | Verified. (Bibitem includes Track-D note explaining the prior title/year correction history.) |
| `sewell1982` | "Quantum fields on manifolds: PCT and gravitationally induced thermal states," Annals of Physics 141, 201 (1982) | CITE-OK | Verified. (Bibitem note documents the prior 1980→1982 correction.) |
| `sewell1980` | Backward-compat alias for the 1982 Sewell paper | CITE-OK (intentional alias) | The bibitem text correctly explains the alias. |
| `hartle_hawking1976` | "Path-integral derivation of black-hole radiance," Phys. Rev. D 13, 2188 (1976) | CITE-OK | Standard. |
| `unruh1976` | "Notes on black-hole evaporation," Phys. Rev. D 14, 870 (1976) | CITE-OK | Standard. |
| `strohmaier2006` | "On noncommutative and pseudo-Riemannian geometry," J. Geom. Phys. 56, 175 (2006); arXiv:math-ph/0110001 | CITE-OK (not independently verified) | Standard. |
| `bizi_brouder_besnard2018` | "Space and time dimensions of algebras…and quantum electrodynamics," JMP 59, 062303 (2018); arXiv:1611.07062 | **CITE-OK** | **The title matches the actual paper** (verified via arxiv.org/abs/1611.07062). The known-bad-pattern wrong title ("Towards a noncommutative geometry of the Standard Model with neutrino mixing") does NOT propagate to Paper 34. |
| `franco_eckstein2014` | "Temporal Lorentzian spectral triples," Reviews in Math. Phys. 26, 1430007 (2014); arXiv:1210.6575 | CITE-OK (not independently verified) | Standard. |
| `van_den_dungen2016` | "Krein spectral triples and the fermionic action," Math. Phys. Anal. Geom. 19, 4 (2016) | CITE-OK (not independently verified) | Standard. |
| `devastato_lizzi_martinetti2018` | "Lorentz signature and twisted spectral triples," JHEP 03, 089 (2018); arXiv:1710.04965 | CITE-OK (not independently verified) | Standard. |
| `nieuviarts2025a` | **B**. Nieuviarts, "Emergence of Lorentz symmetry from an almost-commutative twisted spectral triple," arXiv:2502.18105 (2025) | **CITE-WRONG-METADATA** | The author is **Gaston Nieuviarts** (verified via arXiv author page on 2502.18105); initial should be **G.**, not **B.** **MEDIUM.** |
| `nieuviarts2025b` | **B**. Nieuviarts, "Emergence of Time from a Twisted Spectral Triple…," arXiv:2512.15450 (2025) | **CITE-WRONG-METADATA** | Same author-initial error: **G.**, not **B.** **MEDIUM.** Both `nieuviarts2025a/b` carry the same wrong initial. |
| `paper38`, `paper42` | GeoVac internal | CITE-OK | Self-references. |

**Inline (non-bibitem) literature citations in §V autopsies:** Pachucki, Friar, Eides 2024, Krauth 2017, Antognini, Komasa-Pachucki 2020, Karshenboim 2005 — all appear as inline author-year mentions without formal bibitems. The known-bad patterns for these (PRL 130 023004, PRA 014501, etc.) cannot propagate here because there are no formal cite keys to be wrong about. One minor finding: the bibitem text and Paper 34's §V references attribute the 2020 deuteron quadrupole moment Q_d = 0.285699 to **Komasa–Pachucki** (two authors), but the actual paper is **Puchalski, Komasa, Pachucki (PRL 125, 253001)** (three authors). Recommend adding Puchalski. **LOW.**

**Bottom line: RED** — one CITE-MISATTRIBUTED (chamseddine_connes2010, propagated from Paper 28/33), two CITE-WRONG-METADATA (Nieuviarts initial × 2).

---

### Paper 35 — Time as Projection (21 bibitems)

| cite key | claimed as | verdict | what I found |
|---|---|---|---|
| `paper0,7,14,18,22,23,24,28,31,32,33,34` | GeoVac internal | CITE-OK | Self-references. |
| `kg_sprint_memo` | GeoVac internal memo | CITE-OK | Self-reference. |
| `birrell_davies` | Birrell & Davies, *Quantum Fields in Curved Space*, CUP 1982 | CITE-OK | Standard. |
| `bytsenko` | Bytsenko, Cognola, Elizalde, Moretti, Zerbini, *Analytic Aspects of Quantum Fields*, World Scientific 2003 (§4.5) | CITE-OK (not independently verified, but standard reference) | Standard. |
| `dowker_critchley` | Dowker–Critchley, "Effective Lagrangian and energy-momentum tensor in de Sitter space," Phys. Rev. D 13, 3224 (1976) | CITE-OK | Verified. |
| `ford` | Ford, "Quantum Vacuum Energy in General Relativity," Phys. Rev. D 11, 3370 (1975) | CITE-OK | Verified. |
| `camporesi_higuchi` | Camporesi–Higuchi, J. Geom. Phys. 20, 1 (1996) | CITE-OK | Standard. |
| `bisognano_wichmann1976` | JMP 17, 303 (1976) | CITE-OK | Verified. |
| `sewell1982` | Annals of Physics 141, 201 (1982) | CITE-OK | Verified. |
| `hartle_hawking1976` | Phys. Rev. D 13, 2188 (1976) | CITE-OK | Standard. |
| `unruh1976` | Phys. Rev. D 14, 870 (1976) | CITE-OK | Standard. |

**Bottom line: GREEN** — clean bibliography. No chamseddine_connes2010 here (which would have triggered the propagation).

---

### Paper 20 — Resource Benchmarks (BibTeX external file)

| cite key | source | verdict |
|---|---|---|
| `McArdle2020`, `Wecker2015`, `Rubin2018`, `vonBurg2021`, `Lee2021`, `Bravyi2017`, `Trenev2025`, `Sunaga2025`, `NIST_ASD`, `Tung2011`, plus six `GeoVac_Paper*` keys | All cited via `\cite{...}` against external file `paper_20_refs.bib` | **NO BIBLIOGRAPHY FILE — cannot verify** |

**Critical finding:** Paper 20 references `\bibliography{paper_20_refs}` at line 791, but **the file `paper_20_refs.bib` does NOT exist anywhere in the repository.** I searched with `Glob` for `**/paper_20_refs*` and got zero hits. Paper 20 cannot LaTeX-compile cleanly until this bib file is created or restored.

The cite keys themselves are plausible (McArdle, Wecker, Rubin, vonBurg, Bravyi, Lee, Sunaga, Trenev — all standard quantum-computing references). The 1-norm/Pauli-count benchmark comparisons in the text against Sunaga et al. 2025 RaH-18q (47,099 Pauli relativistic, 12,556 non-relativistic) are load-bearing claims that need the bibitem to be visible to any external reviewer.

**Bottom line: RED** — bibliography file missing. Recommend either restoring/recreating `paper_20_refs.bib` or replacing the BibTeX block with an inline `\begin{thebibliography}` like the other papers.

---

## NEW misattributions discovered (beyond the known-bad list)

1. **`chamseddine_connes2010` → wrong DOI in 3 papers (Papers 28, 33, 34).** Phys. Rev. D 83, 045001 (2011) is *Collective Excitations of Supersymmetric Plasma* (Czajka, Mrowczynski), not Chamseddine–Connes. The intended paper is Fortsch. Phys. 58 (2010) 553. Paper 51 has the correct citation under the same cite key; the three papers are the corrupted copies. **HIGH** — this is the same class of error (real-paper title on wrong DOI) as the Fursaev–Solodukhin fabrication that triggered this whole review.

2. **`perez_sanchez_2024` bibitem in Paper 2 claims the 2024 paper establishes Yang–Mills-without-Higgs.** This is a CITE-DOESNT-SUPPORT — the 2024 paper says a Higgs DOES emerge; only the 2025 Comment provides the without-Higgs correction. Paper 25 has a milder version (in-text bracketing). Paper 30 has it correct.

3. **`nieuviarts2025a/b` in Paper 34 use initial `B.` for the author Gaston Nieuviarts.** Should be `G.`. (Confirmed against arXiv:2502.18105 author page.) **MEDIUM.**

4. **`Peruzzo2014` in Paper 26: article number 4213 → actual 5213.** **MEDIUM** numerical typo.

5. **`eides2001` in Paper 36 conflates the 2007 Springer book title with the 2001 Physics Reports article pagination.** The journal article actual title is "Theory of Light Hydrogenlike Atoms" (not "Theory of Light Hydrogenic Bound States"). Title and journal don't match. **MEDIUM.**

6. **Cite-key year typos (cosmetic, content correct):** `parker1980` (Paper 28, actual 1984), `fegan1983` (Paper 30, actual 1978). **LOW.**

7. **Komasa–Pachucki 2020 should be Puchalski–Komasa–Pachucki 2020.** Paper 34 §V.2 inline citation omits Puchalski as first author of PRL 125, 253001. **LOW.**

8. **Paper 20 bibliography file `paper_20_refs.bib` is missing from the repository.** **HIGH** — paper cannot compile cleanly. Not strictly a misattribution but a structural bibliography failure that affects every cite in Paper 20.

---

## Bottom-line summary

- **Total citations checked:** ~190 bibitems across 11 papers (Paper 20 has BibTeX cites against a missing file, ~10 unique keys not formally verifiable here).
- **Verdict tally:**
  - CITE-MISATTRIBUTED: 1 distinct entry × 3 papers (`chamseddine_connes2010` Phys Rev D 83 045001 in P28, P33, P34) → **HIGH**.
  - CITE-DOESNT-SUPPORT: 2 instances (P2 and P25 on `perez_sanchez_2024` Yang–Mills-without-Higgs claim) → **MEDIUM**.
  - CITE-WRONG-METADATA: 4 verified instances (Paper 26 Peruzzo article number 4213→5213; Paper 36 Eides 2001 mixed title/journal; Paper 34 Nieuviarts initial B→G ×2) → **MEDIUM**.
  - Cite-key year typos with correct content: 2 (P28 `parker1980`, P30 `fegan1983`) → **LOW**.
  - Missing co-author: 1 (P34 inline Komasa-Pachucki → Puchalski-Komasa-Pachucki) → **LOW**.
  - Missing bibliography file: 1 (Paper 20's `paper_20_refs.bib`) → **HIGH** (structural).
  - CITE-CANT-FIND: 0.

- **Per-paper traffic-light:**
  - GREEN: Paper 30, Paper 35, Paper 41, Paper 51.
  - YELLOW: Paper 2, Paper 25, Paper 26, Paper 36.
  - RED: Paper 20 (missing .bib), Paper 28, Paper 33, Paper 34 (all from the chamseddine_connes2010 propagation; Paper 34 additionally has the Nieuviarts initial error).

- **Most consequential finding:** the `chamseddine_connes2010` DOI corruption is the same class of failure as the Fursaev–Solodukhin / hep-th/9512134 fabrication that prompted this audit. Paper 51 is the gold copy; Papers 28, 33, 34 should be brought into sync. Recommended fix: replace bibitem text in P28, P33, P34 with Paper 51's correct entry (Fortsch. Phys. 58, 553 (2010)).

- **Propagation count of known-bad patterns from the audit briefing:** very low. The Lorentzian/Mondino–Sämann cluster known-bad keys (hekkelman_*, ucp_maps_2024, zhu_casini2020, latremoliere_metric_st_2017, leimbach_vs2024, avery_wen_avery1986) do not appear in any of these 12 papers — they are concentrated in the Lorentzian-arc Papers 42–50, which are not in this Wave-B sweep. The precision-physics known-bad keys (Pachucki2023, PachuckiYerokhin2010, Eides2024, FriarPayne2005) do not appear as formal bibitems in any of these 12 papers either — those precision-AMO references are cited inline (author-year) in Paper 34 §V, not via formal cite keys. So most of the known-bad propagation footprint lands outside the Wave-B paper set.

- **The new HIGH finding (`chamseddine_connes2010`/PRD 83 045001) is the dominant signal.** It propagates to 3 of 12 papers, and it's load-bearing because all three papers use this citation to establish lineage to the Chamseddine–Connes spectral action program — exactly the kind of citation that an expert reviewer would check first. Fixing this is the single highest-leverage action coming out of this sweep.
