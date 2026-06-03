# Citation-metadata verification — five-question audit

Date: 2026-06-02
Scope: Settle five GeoVac bibliography metadata questions against published
sources, deciding whether earlier audit-agent corrections were RIGHT, WRONG, or
PARTIAL.

---

## Q1. Latrémolière, "The Gromov–Hausdorff propinquity for metric spectral triples"

**Settled metadata.**
- Author: Frédéric Latrémolière
- Title: "The Gromov–Hausdorff propinquity for metric Spectral Triples"
- Journal: Advances in Mathematics, **volume 404 (2022), Paper No. 108393, 56 pp.**
- Preprint: arXiv:1811.10843 (Nov 2018, multiple revisions through 2022)

**Verdict on prior audit's claim ("Adv. Math. 404 (2022), 108393, 56 pp"):
RIGHT.**

The current GeoVac bibitem ("Adv. Math. 415 (2023), Paper No. 108876, 88 pp")
is incorrect on **all four** numerical fields (volume, year, paper number, page
count). The audit-proposed correction is the actual published metadata.

**URLs used.**
- https://arxiv.org/abs/1811.10843 (arXiv page lists "Journal ref: Adv. Math.
  404 (2022), Paper No. 108393, 56 pp")
- https://www.sciencedirect.com/science/article/abs/pii/S0001870822002092
  (Advances in Mathematics article landing page, confirms volume/year)

---

## Q2. Latrémolière, "The Quantum Gromov–Hausdorff Propinquity"

**Settled metadata.**
- Author: Frédéric Latrémolière
- Title: "The Quantum Gromov–Hausdorff Propinquity"
- Journal: **Trans. Amer. Math. Soc. 368 (2016), no. 1, pp. 365–411.**
- Preprint: arXiv:1302.4058 (Feb 2013, last revision Nov 2013)
- Article DOI prefix: S0002-9947-2015-06334-X (AMS recordkeeping)

**Verdict on prior audit's claim ("Trans. AMS 368 (2016)"): RIGHT.**

The current GeoVac bibitem in Papers 47, 48, 49 ("Trans. Amer. Math. Soc. 370
(2018), 365–411") is wrong on both volume (370 → 368) and year (2018 → 2016).
Page range 365–411 is correct. The audit-proposed correction is the actual
published metadata.

**URLs used.**
- https://arxiv.org/abs/1302.4058
- https://www.ams.org/journals/tran/2016-368-01/S0002-9947-2015-06334-X/
  (AMS journal record for vol. 368 no. 1, 2016)

---

## Q3. Leimbach first name

**Settled metadata.**
- First name: **Malte** (not Lukas).
- Full citation: Malte Leimbach and Walter D. van Suijlekom,
  "Gromov–Hausdorff convergence of spectral truncations for tori",
  Adv. Math. 439 (2024), 109496. arXiv:2302.07877.

**Verdict on prior audit's claim ("Malte, not Lukas"): RIGHT.**

The arXiv abstract page and the Advances in Mathematics ScienceDirect landing
page both list the author as **Malte Leimbach**. The five GeoVac papers
(38, 39, 40, 42, 43) that thank "Lukas Leimbach" in acknowledgments are wrong
on the first name.

**URLs used.**
- https://arxiv.org/abs/2302.07877 (author list: "Malte Leimbach and Walter
  D. van Suijlekom")
- https://www.sciencedirect.com/science/article/pii/S0001870824000112
  (Advances in Mathematics, vol. 439, article 109496, 2024)
- https://arxiv.org/abs/2302.07877v2

---

## Q4. Hekkelman 2022

**Settled metadata.**

(a) arXiv:2206.13744 is **NOT a Hekkelman paper.** It is "Image of Kerr-Melvin
black hole with thin accretion disk" by Yehui Hou, Zhenyu Zhang, Haopeng Yan,
Minyong Guo, and Bin Chen — a gr-qc paper on black-hole imaging with no
relation to spectral truncations or noncommutative geometry.

(b) arXiv:2111.13865 IS a real Hekkelman work:
- Author: Eva-Maria Hekkelman
- Title: "Truncated Geometry on the Circle"
- Submitted Nov 27, 2021; revised March 1, 2022
- Subject classification: math.QA, math-ph, math.OA
- Content: Toeplitz-matrix pure state spaces converge in Gromov-Hausdorff
  sense to state space on C(S¹) under Connes-distance metrics. This is the
  content GeoVac Paper 44's `hekkelman2022` citation needs.

(c) Hekkelman's Radboud PhD thesis "Truncated geometries from spectral
triples" is a separate work (the survey/thesis-level treatment), not the
same as 2111.13865.

Other Hekkelman 2022+ papers from her author page (https://emhekkelman.nl/papers/):
- "An application of singular traces to crystals and percolation"
  (with Azamov, McDonald, Sukochev, Zanin), J. Geom. Phys. 2022 — different
  topic (crystal/percolation singular traces).
- arXiv:2404.16338 "Multiple operator integrals, pseudodifferential
  calculus, and asymptotic expansions" (2024).
- arXiv:2412.00628 "A noncommutative integral on spectrally truncated
  spectral triples, and a link with quantum ergodicity" (2024).

**Verdict on prior audit's claim ("2206.13744 is wrong, correct is
2111.13865"): RIGHT.**

The GeoVac Paper 44 `hekkelman2022` bibitem pointing to arXiv:2206.13744 is
a wrong arXiv ID (off by ~3000 IDs into an unrelated gr-qc paper). The
correct arXiv ID for the Hekkelman 2022 spectral-truncation work cited in
Paper 44's spectral-truncation lineage is **2111.13865**, "Truncated Geometry
on the Circle".

**URLs used.**
- https://arxiv.org/abs/2206.13744 (Kerr-Melvin BH paper, confirmed unrelated)
- https://arxiv.org/abs/2111.13865 (Hekkelman, "Truncated Geometry on the
  Circle")
- https://emhekkelman.nl/papers/ (Hekkelman author publications list)

---

## Q5. Avery-Wen-Avery citation

**Settled metadata.**
- Authors: **Zhen-Yi Wen and John Avery** (2 authors, not 3).
- Title: "Some properties of hyperspherical harmonics"
- Journal: **J. Math. Phys. 26 (1985), no. 3, pp. 396–403** (March 1, 1985 issue).
- AIP DOI: 10.1063/1.526621

**Verdict on prior audit's claim ("Wen & Avery, JMP 26 (1985), 2-author"):
RIGHT.**

The GeoVac Paper 44 `avery_wen_avery1986` bibitem is wrong on three counts:
- Author count (3 → 2): there is no "Avery-Wen-Avery" 3-author paper; the
  canonical hyperspherical-Sturmian paper from this group in this period is
  the 2-author Wen-Avery paper.
- Year (1986 → 1985).
- Volume (27 → 26).

Page range and journal name (JMP) are correct. The bibitem key
`avery_wen_avery1986` itself should be renamed (e.g. `wen_avery1985`) to
match the actual citation.

**URLs used.**
- https://pubs.aip.org/aip/jmp/article/26/3/396/376832/Some-properties-of-hyperspherical-harmonics
  (AIP Journal of Mathematical Physics article landing page; confirms
  authors, title, volume, issue, year, pages)
- https://www.semanticscholar.org/paper/Some-properties-of-hyperspherical-harmonics-Wen-Avery/c519a5a0df2c11275a641c52ad9bb49e514680f1
  (Semantic Scholar record corroborating the two-author Wen-Avery citation)

---

## Summary table

| Q | Current GeoVac bibitem | Verified truth | Verdict |
|:--|:-----------------------|:---------------|:--------|
| Q1 | Adv. Math. 415 (2023), 108876, 88pp | Adv. Math. 404 (2022), 108393, 56pp | RIGHT |
| Q2 | Trans. AMS 370 (2018), 365–411 | Trans. AMS 368 (2016), no. 1, 365–411 | RIGHT |
| Q3 | "Lukas Leimbach" | "Malte Leimbach" | RIGHT |
| Q4 | hekkelman2022 → arXiv:2206.13744 | arXiv:2111.13865 ("Truncated Geometry on the Circle") | RIGHT |
| Q5 | "Avery-Wen-Avery" 3-author, JMP 27 (1986) | Wen & Avery 2-author, JMP 26 (1985), 396–403 | RIGHT |

All five audit-proposed corrections are confirmed by published sources. The
existing bibitems in the GeoVac corpus are wrong as audited; the corrections
should be applied as proposed.
