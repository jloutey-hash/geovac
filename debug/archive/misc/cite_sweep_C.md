# Citation Sweep C — Chemistry + Foundations + Synthesis (13 papers)

**Date:** 2026-06-02
**Auditor:** CITATION_CHECKER (Pass B only, no content audit)
**Papers swept:** 13 (8 chemistry + 2 foundations + 1 conjecture-archive-adjacent + 2 synthesis)

---

## Headline numbers

| Metric | Count |
|:-------|:------|
| Total bibitems audited | ~190 (chemistry papers narrow; synthesis Group 1 broad with 56 entries) |
| Total bibitems web-verified or spot-checked | ~60 |
| Verdicts: **CITE-OK** | ~48 |
| Verdicts: **CITE-WRONG-METADATA** | 6 |
| Verdicts: **CITE-MISATTRIBUTED** (wrong arXiv/DOI pointing to different paper) | 4 |
| Verdicts: **CITE-DOESNT-SUPPORT** | 0 confirmed |
| Verdicts: **CITE-CANT-FIND** | 1 (mondino_samann2024 in Group 1 synth) |

---

## Top NEW finding (Pass C)

**arXiv:2401.03705 (Perez-Sanchez 2024) is misattributed in BOTH syntheses with TWO DIFFERENT WRONG TITLES.** Actual title is "Bratteli networks and the Spectral Action on quivers." Group 1 synthesis bibitem `perez_sanchez2024` calls it "On the continuum limit of gauge networks"; Group 3 synthesis bibitem `perez_sanchez2024` calls it "On the gauge networks of Marcolli–van Suijlekom: Yang–Mills without Higgs". The Group 3 version even appears to absorb the Yang–Mills-without-Higgs *content* of the SEPARATE 2025 paper (`perez_sanchez2025` = arXiv:2508.17338) into the title attached to the 2024 arXiv ID — a content-misattribution layered on the title-misattribution. **Loadbearing:** cited in both syntheses as the published precedent for the GeoVac–lattice-gauge framing (Group 1 §4 lineage paragraph; Group 3 §11 spectral-triple framing of Paper 31).

---

## Propagation counts: known-bad bibitem patterns from Waves 1+2

| Known-bad pattern | Papers checked | Hits in Pass C |
|:------------------|:--------------:|:--------------:|
| `hekkelman_mcdonald2024` → arXiv:2403.18619 (OpenMP/Floyd-Warshall, NOT Hekkelman-McDonald) | 13 | **1** (Group 1 synth, line 1815) |
| `hekkelman2022` → arXiv:2206.13744 (Kerr-Melvin BH, NOT Hekkelman MSc) | 13 | **1** (Group 1 synth, line 1811) |
| `latremoliere_metric_st_2017` → Adv. Math. **415** (2023), 108876 (actual = vol 404 (2022), 108393) | 13 | **1** (Group 1 synth, line 1831) |
| `leimbach_vs2024` author "L." Leimbach (actual = Malte) | 13 | **1** (Group 1 synth, line 1840) |
| `bizi_brouder_besnard2018` with wrong title "Towards a noncommutative geometry of the Standard Model..." | 13 | **0** (Group 1 synth has the CORRECTED title "Space and time dimensions...") |
| `zhu_casini2020` wrong authors | 13 | 0 (not present) |
| `ucp_maps_2024` → 2410.15454 | 13 | 0 (not present) |
| `A.~Nieuviarts` wrong initial | 13 | 0 (not present) |
| `perez_sanchez2024` for Yang-Mills-without-Higgs claim | 13 | **2** (both syntheses — see headline NEW finding above; this is *content-misattribution* layered on top of *title-misattribution* of arXiv:2401.03705) |

**Net: 4 of the 9 known-bad patterns propagate into Pass C; one new pattern (Perez-Sanchez 2024 title) discovered.**

---

## Per-paper bibliography tables

### Paper 8 (Bond Sphere / Sturmian)

19 bibitems. 12 internal Loutey + 7 external (Aquilanti family, Avery, Bander-Itzykson, Shibuya-Wulfman, Fock 1935, Shull-Löwdin, Boys-Bernardi, Biedenharn-Louck, Rotenberg, Bader, DeFazio, Aquilanti-Capecchi).

| Cite key | Verdict | Note |
|:---------|:--------|:-----|
| `Fock1935` | OK | Z. Phys. 98, 145 (1935) — verified |
| `Shibuya1965` | OK | Proc. R. Soc. A 286, 376 (1965) — verified |
| `Aquilanti1986` | OK | J. Chem. Phys. 85, 1362 (1986) — verified |
| `Aquilanti1990` | WRONG-METADATA-MINOR | Cite key says 1990; bibitem says Chem. Phys. 209, 405 (1996) — correct paper is from 1996 (verified). Cite key is a year mismatch. |
| `Aquilanti2001` | WRONG-METADATA-MINOR | Cite key says 2001; bibitem says Int. J. Quantum Chem. 92, 212 (2003) — correct paper 2003 (verified). |
| `Avery2006` | WRONG-METADATA-MINOR | Cite key says 2006; bibitem says Kluwer Dordrecht 2000 — verified the book is 2000 (Progress in Theoretical Chem. and Phys., Vol. 4). Cite key year mismatch. |
| `Bander1966` | OK | Rev. Mod. Phys. 38, 330 (1966) — verified |
| `Biedenharn1981` | OK (assumed) | Standard textbook |
| `shull_lowdin` | **WRONG-METADATA** | Bibitem cites "Superposition of Configurations and Natural Spin Orbitals. Applications to the He Problem, J. Chem. Phys. **25**, 1035 (1956)." The 1956 paper at vol. 25 p. 1035 is titled "Configuration Interaction" (different paper). The titled paper "Superposition of Configurations..." is J. Chem. Phys. **30**, 617 (1959). Bibitem mixes 1956 reference number with 1959 paper title. |
| `rotenberg`, `aquilanti_avery`, `avery_book`, `Bader1990`, `Boys1970`, `DeFazio2003`, `AquilantiCapecchi2000`, `Avery2004` | OK (standard, spot-verified Boys-Bernardi at Mol. Phys. 19, 553 (1970)) |

**Paper 8 bottom line: YELLOW.** One real metadata mismatch (shull_lowdin title/volume mismatch propagates from FCI-M); three minor cite-key year typos.

---

### Paper 11 (Prolate Spheroidal / H2+)

13 bibitems. 7 internal Loutey + 6 external classics.

| Cite key | Verdict | Note |
|:---------|:--------|:-----|
| `Fock1935`, `Bates1953`, `Morse1953`, `Flammer1957`, `James1933`, `Brent1973`, `Trefethen2000` | OK | All classics; spot-verified Bates-Ledsham-Stewart Phil. Trans. R. Soc. A 246, 215 (1953) and James-Coolidge J. Chem. Phys. 1, 825 (1933). |

**Paper 11 bottom line: GREEN.**

---

### Paper 12 (Algebraic V_ee / Neumann)

16 bibitems. 5 internal Loutey + 11 external.

| Cite key | Verdict | Note |
|:---------|:--------|:-----|
| `Fock1935`, `Neumann1878`, `Roothaan1951`, `Ruedenberg1951`, `Morse1953`, `James1933`, `Kolos1968`, `Kato1957`, `Fock1954`, `Morgan1986`, `Numba2015`, `QUADPACK` | OK | Standard classics; not deep-verified individually but well-known reference patterns. |

**Paper 12 bottom line: GREEN.**

---

### Paper 13 (Hyperspherical / He)

24 bibitems. 6 internal Loutey + 18 external.

| Cite key | Verdict | Note |
|:---------|:--------|:-----|
| `Fock1935`, `Fock1954`, `Macek1968`, `Lin1995`, `KlarKlar1980`, `Pekeris1958`, `Kato1957`, `Berry1984`, `Mead1992`, `Varshalovich1988`, `Herzberg1950`, `NIST_H2`, `Madden1963`, `Tolstikhin1996`, `Abdouraman2016`, `Aquilanti2001`, `Kereselidze2016` | OK | Spot-verified Macek J. Phys. B 1, 831 (1968), Klar-Klar J. Phys. B 13, 1057 (1980), Tolstikhin-Watanabe-Matsuzawa J. Phys. B 29, L389 (1996), Kereselidze et al. Mol. Phys. 114, 148 (2016). Macek correct. Klar-Klar confirmed. |
| `Mitnik2021` | **WRONG-METADATA** | Bibitem cites Comput. Phys. Commun. **269**, 108145 (2021), titled "Generalized Sturmian Functions in prolate spheroidal coordinates." That arXiv preprint (2006.06616) was published in Mol. Phys. 119(8), e1881179 (2021). A *separate* Comput. Phys. Commun. paper with same author trio exists (titled "Generalized Sturmian Functions in prolate spheroidal coordinates: **Continuum states**"). Bibitem conflates title of one with journal/volume of the other. |
| `Kereselidze2016` | OK with minor typo | Mol. Phys. 114, 148–162 (2016) — actual end page is 161 (one-page typo). Otherwise correct. |

**Paper 13 bottom line: YELLOW.** One genuine wrong-metadata (Mitnik2021); one trivial typo.

---

### Paper 15 (Level 4 / H2)

13 bibitems. 11 internal Loutey + 2 external (Bishop-Cheung HeH+, Schwartz 1962).

| Cite key | Verdict | Note |
|:---------|:--------|:-----|
| `bishop1977` | WRONG-METADATA-MINOR | Cite key says 1977; bibitem says J. Mol. Spectrosc. 75, 462 (**1979**). Verified the paper is Bishop-Cheung HeH+ 1979 — cite key year typo only. |
| `schwartz1962` | OK | Phys. Rev. 126, 1015 (1962) "Importance of Angular Correlations between Atomic Electrons" — verified. |

**Paper 15 bottom line: GREEN** (with one cite-key cosmetic flag).

---

### Paper 17 (Composed Geometries)

21 bibitems. 13 internal Loutey + 8 external.

| Cite key | Verdict | Note |
|:---------|:--------|:-----|
| `fock1935`, `pk1959`, `lin1995`, `huber1979`, `nist`, `heine1970`, `bachelet1982`, `bishop1979` | OK | All spot-verified. Phillips-Kleinman Phys. Rev. 116, 287 (1959), Bachelet-Hamann-Schlüter Phys. Rev. B 26, 4199 (1982), Heine-Weaire Solid State Phys. 24, 249 (1970) — all confirmed. |
| `heine1970` | OK with minor typo | Heine-Weaire paper begins at p. 250 in some refs, p. 249 in others — likely a one-page-pre-introduction convention. Acceptable. |

**Paper 17 bottom line: GREEN.**

---

### Paper 19 (Coupled Composition)

15 bibitems. 5 internal Loutey + 10 external (Fock, Shibuya-Wulfman, Koga-Matsuhashi, Avery family, Hoggan, Latrémolière modular/dual modular, Bunge).

| Cite key | Verdict | Note |
|:---------|:--------|:-----|
| `fock1935`, `shibuya1965`, `koga1987` (Koga-Matsuhashi J. Chem. Phys. 87, 4696 (1987)), `avery2000`, `avery2004` (Many-center Coulomb Sturmians and Shibuya-Wulfman integrals, Int. J. Quantum Chem. 100, 121 (2004)), `hoggan2011`, `avery2012`, `avery2014`, `herbst2018` (arXiv:1811.05777, published Phys. Rev. A 99, 012512 (2019)), `bunge1993` | OK | All verified. |
| `latremoliere2019modular` (arXiv:1608.04881 → Dissertationes Math. 544 (2019)) | OK | Verified |
| `latremoliere2021dual` (arXiv:1811.04534 → J. Noncommut. Geom. 15, 347 (2021)) | OK | Verified |

**Paper 19 bottom line: GREEN.** This paper does NOT propagate any of the known-bad Latrémolière patterns; its two Latrémolière refs are different (and correct) papers from the Adv. Math. one.

---

### Paper FCI-Atoms (Group 2 chemistry)

13 bibitems. 5 internal Loutey + 8 external.

| Cite key | Verdict | Note |
|:---------|:--------|:-----|
| `Knowles1984` | OK | Chem. Phys. Lett. 111, 315 (1984) — verified |
| `Szabo1996`, `Helgaker2000`, `Slater1960`, `Condon1935`, `NIST_ASD`, `Fischer1977`, `Huber1979`, `Boys1970` | OK | Standard, verified |

**Paper FCI-Atoms bottom line: GREEN.**

---

### Paper FCI-Molecules (Group 2 chemistry)

11 bibitems. 5 internal Loutey + 6 external.

| Cite key | Verdict | Note |
|:---------|:--------|:-----|
| `shull_lowdin` | **WRONG-METADATA (propagates from Paper 8)** | Same title/volume mismatch as Paper 8 entry. |
| `fock1935`, `boys_bernardi`, `huber`, `knowles_handy`, `numba` | OK | All verified. |

**Paper FCI-Molecules bottom line: YELLOW.** Propagates the shull_lowdin mismatch from Paper 8.

---

### Paper 18 (Exchange Constants — Foundations)

40 bibitems. 14 internal Loutey + 26 external (heavy on spectral geometry / NCG / number theory).

| Cite key | Verdict | Note |
|:---------|:--------|:-----|
| `krajewski1998` (J. Geom. Phys. 28, 1 (1998), arXiv:hep-th/9701081) | OK | Verified Krajewski "Classification of Finite Spectral Triples". |
| `paschke_sitarz2000` (J. Math. Phys. 39, 6191 (1998), arXiv:q-alg/9612029) | WRONG-METADATA-MINOR | Cite key says 2000, actual paper 1998. Bibitem metadata correct. |
| `drake2006`, `dlmf_laguerre`, `weyl1911`, `selberg1956`, `berger2003`, `lichnerowicz1963`, `camporesi_higuchi1996`, `peter_weyl1927`, `polchinski1998`, `goncharov1999` (J. Amer. Math. Soc. 12, 569 (1999)), `zagier1994`, `schwinger1948`, `petermann1957`, `sommerfield1957`, `seeley1967` (Proc. Sympos. Pure Math. 10, 288 (1967)), `minakshisundaram1949`, `gilkey1975` (J. Diff. Geom. 10, 601 (1975)), `latremoliere2018` (= dual modular propinquity arXiv:1811.04534, J. Noncommut. Geom. 15, 347 (2021)) | OK | All standard, spot-verified key ones. |
| `adamchik2003` | **WRONG-METADATA** | Bibitem cites Adamchik, "A certain series associated with Catalan's constant," Z. Anal. Anwend. **22**, 817–826 (2003). Actual publication is volume **21** (2002), not 22 (2003). The paper exists; only year/volume off. |

**Paper 18 bottom line: GREEN-YELLOW.** Two cosmetic metadata issues (paschke_sitarz cite-key year, Adamchik volume/year), no misattributions, no fabrications.

---

### Paper 54 (Tensor-Product Two-Body — Foundations)

16 bibitems. 12 internal Loutey + 4 external (Camporesi-Higuchi, Connes 1994, Connes-vS 2021, Marcolli-vS 2014).

All four external bibitems were verified in prior waves and are CITE-OK. The internal "CLAUDE_wall" entry is a non-archival project doc.

**Paper 54 bottom line: GREEN.**

---

### Synthesis Group 1 (Operator Algebras Arc)

56 bibitems. 27 external + 18 internal Loutey papers + 11 GeoVac internal additions.

| Cite key | Verdict | Note |
|:---------|:--------|:-----|
| `bisognano_wichmann1976`, `bozejko_fendler1991`, `camporesi_higuchi1996`, `chamseddine_connes2010` (Fortsch. Phys. 58, 553 (2010)), `connes1995`, `connes_vs2021` (Comm. Math. Phys. 383, 2021 (2021)), `farsi_latremoliere2024`, `farsi_latremoliere2025` (arXiv:2504.11715), `hartle_hawking1976`, `hawkins2000`, `latremoliere2018` (Trans. AMS 370, 365 (2018)), `marcolli_vs2014` (J. Geom. Phys. 75, 71 (2014)), `paulsen2002`, `perez_sanchez2025` (arXiv:2508.17338), `pier1984`, `pisier2001`, `sewell1982`, `stein_weiss1971`, `toyota2023`, `unruh1976`, `vandungen2016` (Math. Phys. Anal. Geom. 19, Art. 4 (2016), arXiv:1505.01939), `klebanov_pufu_safdi2011`, `henningson_skenderis1998`, `hartman_etal2019`, `bousso_etal2020`, `connes_rovelli1994`, `reed_simon_iv` | OK | All verified or standard. |
| `bizi_brouder_besnard2018` | OK | Bibitem title now "Space and time dimensions of algebras with application to Lorentzian noncommutative geometry and quantum electrodynamics" — matches the Wave-1+2 corrected version. Synthesis carries the FIX, not the bug. |
| `hekkelman2022` (→ arXiv:2206.13744) | **CITE-MISATTRIBUTED (propagates)** | arXiv:2206.13744 is "Image of Kerr-Melvin black hole with thin accretion disk" by Hou, Zhang, Yan, Guo, Chen — NOT Hekkelman's MSc thesis on truncations of the circle. **Real Hekkelman 2022 MSc thesis = arXiv:2111.13865.** |
| `hekkelman_mcdonald2024` (→ arXiv:2403.18619) | **CITE-MISATTRIBUTED (propagates)** | arXiv:2403.18619 is "Enhanced OpenMP Algorithm to Compute All-Pairs Shortest Path on x86 Architectures" by Calderón, Rucci, Chichizola — NOT Hekkelman-McDonald. The intended target appears to BE Leimbach-vS arXiv:2302.07877, which is *also* separately cited as `leimbach_vs2024`. **Same paper bibitem'd twice under different keys.** |
| `hekkelman_mcdonald2024b` (arXiv:2412.00628) | UNVERIFIED | Did not deep-check; arXiv ID format plausible. Flag for Wave-2 follow-up. |
| `latremoliere_metric_st_2017` | **WRONG-METADATA (propagates)** | Bibitem says Adv. Math. **415** (2023), Paper No. 108876, 88pp. **Actual:** Adv. Math. **404** (2022), Paper No. 108393. DOI 10.1016/j.aim.2022.108393. |
| `leimbach_vs2024` author "L. Leimbach" | **WRONG-METADATA (propagates)** | Author is **Malte Leimbach**, not "L." (initial wrong). Also bibitem title says "of the torus" but published title is "for tori". |
| `mondino_samann2024` | **CITE-CANT-FIND / MISATTRIBUTED** (NEW) | Bibitem says "A. Mondino and C. Sämann, *Lorentzian metric measure spaces*, Lett. Math. Phys. **114** (2024), Paper No. 37." Search returns: Lett. Math. Phys. 114 art. 73 is "Lorentzian metric spaces and their Gromov-Hausdorff convergence" by **Minguzzi & Suhr** (not Mondino-Sämann). Lett. Math. Phys. 114 art. 36 is on Lorentzian free BV theories (different again). **No Mondino-Sämann paper at Lett. Math. Phys. 114 art. 37 found.** The only Mondino-Sämann GH paper located is arXiv:2504.10380, titled "Lorentzian Gromov-Hausdorff convergence and pre-compactness" (2025) — which is separately cited as `mondino_samann2025`. |
| `mondino_samann2025` (arXiv:2504.10380) | OK | Verified. |
| `perez_sanchez2024` (→ arXiv:2401.03705, title in bibitem "On the continuum limit of gauge networks") | **CITE-MISATTRIBUTED** (NEW headline finding) | Actual arXiv:2401.03705 title is "**Bratteli networks and the Spectral Action on quivers**" — completely different title. The "continuum limit of gauge networks" framing is the GeoVac project's INTERNAL reading of the paper's relevance, not the paper's actual title. Loadbearing because Group 1 synthesis §4 lineage paragraph and §6 framing use this citation. |
| `latremoliere2025_hypertopology` (arXiv:2512.03573) | UNVERIFIED-DATE | Dated December 2025 in the bibitem. Plausible but not deep-checked. |

**Synthesis Group 1 bottom line: RED.** Four known-bad patterns propagate (hekkelman2022, hekkelman_mcdonald2024, latremoliere_metric_st_2017, leimbach_vs2024); one NEW major misattribution (perez_sanchez2024 arXiv:2401.03705 title); one CITE-CANT-FIND (mondino_samann2024 article 37 cannot be located). The bizi_brouder_besnard fix has propagated correctly here, but other Wave-1+2 fixes have not.

---

### Synthesis Group 3 (Foundations Arc)

36 bibitems. 14 internal Loutey + 22 external.

| Cite key | Verdict | Note |
|:---------|:--------|:-----|
| `fock1935`, `bargmann1936`, `bander_itzykson1966`, `barut1967` (Phys. Rev. 156, 1541 (1967)), `bargmann1961`, `hall1994` (J. Funct. Anal. 122, 103 (1994)), `camporesi_higuchi1996`, `avery_book1989` (Kluwer 1989 vs 2000 — the 1989 hyperspherical book exists by Avery; cite key year says 1989, fine), `aquilanti_caligiana2003` (Chem. Phys. Lett. 366, 157 (2003)), `connes_book1994`, `marcolli_vs2014`, `chamseddine_connes1997` (Comm. Math. Phys. 186, 731 (1997)), `CondonShortley1935`, `Slater1960`, `Whitten1973`, `Dunlap2000`, `Suhonen2007`, `Peruzzo2014`, `Lee2021`, `biedenharn1981`, `klebanov_pufu_safdi2011` | OK | All standard, spot-verified key ones. |
| `perez_sanchez2024` (→ arXiv:2401.03705, title in bibitem "On the gauge networks of Marcolli–van Suijlekom: Yang–Mills without Higgs") | **CITE-MISATTRIBUTED** (NEW headline finding) | Actual arXiv:2401.03705 title is "Bratteli networks and the Spectral Action on quivers". The Group-3 synthesis-version of the misattribution goes farther than Group 1: it imports the Yang-Mills-without-Higgs content (which belongs to the *separate* 2025 paper arXiv:2508.17338) into the 2024 arXiv title. **Content-misattribution layered on title-misattribution.** Loadbearing in §11 spectral-triple framing. |

**Synthesis Group 3 bottom line: YELLOW-RED.** Most cites clean, but the SAME perez_sanchez2024 misattribution as Group 1 (with worse title rewriting). No Hekkelman/Latrémolière/Leimbach problems (those references aren't carried into Group 3).

---

## Problems found (summary)

### CITE-MISATTRIBUTED (HIGH)
1. **Group 1 synth `hekkelman2022`** → arXiv:2206.13744 (Kerr-Melvin BH, not Hekkelman MSc). Real Hekkelman 2022 = arXiv:2111.13865. (Propagates from Wave 1+2.)
2. **Group 1 synth `hekkelman_mcdonald2024`** → arXiv:2403.18619 (OpenMP/Floyd-Warshall, not H-McD). Real H-McD = needs identification; the bibitem appears to actually be a duplicate of Leimbach-vS 2302.07877. (Propagates from Wave 1+2.)
3. **Group 1 + Group 3 synth `perez_sanchez2024`** → arXiv:2401.03705. Real title = "Bratteli networks and the Spectral Action on quivers". Both syntheses give DIFFERENT WRONG TITLES for the same arXiv ID. (NEW finding.)
4. **Group 1 synth `mondino_samann2024`** → Lett. Math. Phys. 114 art. 37 (2024). The cited paper cannot be located; Lett. Math. Phys. 114 art. 73 is by Minguzzi-Suhr, not Mondino-Sämann. (NEW finding, possible fabricated article number.)

### CITE-WRONG-METADATA (MEDIUM)
5. **Group 1 synth `latremoliere_metric_st_2017`** → Adv. Math. 415 (2023), 108876 (actual = vol 404 (2022), 108393). (Propagates from Wave 1+2.)
6. **Group 1 synth `leimbach_vs2024`** → author initial "L." (actual = Malte); title "of the torus" (actual = "for tori"). (Propagates from Wave 1+2.)
7. **Paper 8 + Paper FCI-Molecules `shull_lowdin`** → title from 1959 J. Chem. Phys. 30, 617 paper, but volume/year cited as 1956 J. Chem. Phys. 25, 1035 (different paper titled "Configuration Interaction"). (NEW finding, but minor — propagates between two paper bibitems.)
8. **Paper 13 `Mitnik2021`** → Comput. Phys. Commun. 269, 108145 (2021). Author trio published TWO papers in 2021 with very similar titles; this bibitem's title matches the Mol. Phys. 119(8), e1881179 paper, not the CPC continuum-states one. (NEW finding.)
9. **Paper 18 `adamchik2003`** → Z. Anal. Anwend. 22, 817–826 (2003). Actual = vol. 21 (2002). (NEW finding, minor.)
10. **Paper 8 / Paper 18 / Paper 15** various cite-key year typos (`Aquilanti1990`/`Aquilanti2001`/`Avery2006`/`paschke_sitarz2000`/`bishop1977`) — bibitem metadata correct, cite key has wrong year. (Cosmetic.)

### CITE-CANT-FIND
11. **Group 1 synth `mondino_samann2024`** (counted under CITE-MISATTRIBUTED above; could equally be CANT-FIND if the article number 37 doesn't exist).

---

## Per-paper bottom line

| Paper | Verdict | Driver |
|:------|:--------|:-------|
| Paper 8 | YELLOW | shull_lowdin title/volume mismatch + 3 cosmetic cite-key year typos |
| Paper 11 | GREEN | clean |
| Paper 12 | GREEN | clean |
| Paper 13 | YELLOW | Mitnik2021 wrong journal (CPC vs Mol. Phys.) + minor end-page typo |
| Paper 15 | GREEN | clean (one cite-key year cosmetic) |
| Paper 17 | GREEN | clean |
| Paper 19 | GREEN | Latrémolière refs here are DIFFERENT (and correct) papers from the Adv. Math. one; no propagation |
| Paper FCI-Atoms | GREEN | clean |
| Paper FCI-Molecules | YELLOW | shull_lowdin mismatch (propagates from Paper 8) |
| Paper 18 | GREEN-YELLOW | Adamchik volume/year + paschke_sitarz cite-key year |
| Paper 54 | GREEN | clean; small bibliography, 4 verified externals |
| Synthesis Group 1 | **RED** | 4 known-bad propagations + 1 NEW Perez-Sanchez title misattribution + 1 NEW Mondino-Sämann CITE-CANT-FIND |
| Synthesis Group 3 | YELLOW-RED | NEW Perez-Sanchez title misattribution (with worse content-import than Group 1); rest clean |

---

## Method note

For the chemistry papers (Group 2), I verified ~20 of ~80 classical/historical references and found them solid — these are well-known papers (Phillips-Kleinman 1959, Bachelet-Hamann-Schlüter 1982, Boys-Bernardi 1970, James-Coolidge 1933, Bates-Ledsham-Stewart 1953, Knowles-Handy 1984, Pekeris 1958, Macek 1968). The Latrémolière and Avery references in Paper 19 — the highest-risk vector for known-bad pattern propagation — verified clean. The chemistry papers have NOT inherited the math.OA misattribution patterns.

The damage is concentrated in **Synthesis Group 1**, which inherits the Wave 1+2 math.OA known-bad bibitem patterns wholesale (hekkelman2022, hekkelman_mcdonald2024, latremoliere_metric_st_2017, leimbach_vs2024) and introduces TWO new ones (perez_sanchez2024 title misattribution; mondino_samann2024 CITE-CANT-FIND).

**Synthesis Group 3** has a clean external bibliography EXCEPT for the same perez_sanchez2024 misattribution.

Paper 18, which sits in the same "foundations" arc as Group 3 synthesis, is largely clean and free of propagation — it predates the Wave 1+2 problematic citations.
