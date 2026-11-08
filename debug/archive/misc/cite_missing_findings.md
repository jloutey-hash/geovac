# Cite-Missing Findings — 7 GeoVac Bibitem Resolutions

Date: 2026-06-02

Mission: locate intended citations for 7 demonstrably wrong bibitems. Each section reports the best identification, confidence level, URLs, and a recommended replacement bibitem text.

---

## F1. Paper 14 `trenev2025` — source of headline LiH/H2O/BeH2 Pauli counts

**Confidence: HIGH (mixed verdict)**

**Finding:**
- The cited paper (arXiv:2311.03719, Trenev–Ollitrault–Harwood–Gujarati–Raman–Mezzacapo–Mostame, *Quantum* 2025) DOES exist and IS by these authors.
- HOWEVER the paper's full title is "Refining resource estimation for the quantum computation of **vibrational** molecular spectra through Trotter error analysis" — the GeoVac bibitem dropped the word **vibrational**, which is load-bearing: the paper is exclusively about vibrational structure of acetylene-like polyynes.
- The cited paper does **NOT** contain Table 5 with the LiH/H2O/BeH2 STO-3G/cc-pVDZ/cc-pVTZ Pauli counts (276/5851/63519, 551/8921/107382). Verified by reading the arXiv HTML directly.
- The STO-3G values (LiH 276, H2O 551) ARE standard widely-reported numbers from OpenFermion/Qiskit-Nature without a single canonical primary source. The cc-pVDZ/cc-pVTZ values (5851, 8921, 63519, 107382) do not appear in any external paper found in extensive web searches.

**Verdict:** The cited paper is misattributed for these Pauli counts. The counts are **most likely a GeoVac internal computation** (e.g., from the project's own OpenFermion/Qiskit-Nature pipeline), with the Trenev paper added by mistake. The PM should either:
1. Replace `\cite{trenev2025}` with an internal-computation note (no external cite), citing OpenFermion/Qiskit-Nature tooling; OR
2. Cite Trenev 2025 only as a comparison METHOD reference (resource-estimation methodology, not source of the specific LiH/H2O counts).

**URLs:**
- https://arxiv.org/abs/2311.03719
- https://arxiv.org/html/2311.03719v2
- https://quantum-journal.org/papers/q-2025-02-11-1630/

**Recommended bibitem text (corrected title):**
```
\bibitem{trenev2025}
D.~Trenev, P.~J.~Ollitrault, S.~M.~Harwood, T.~P.~Gujarati,
S.~Raman, A.~Mezzacapo, and S.~Mostame,
``Refining resource estimation for the quantum computation of
vibrational molecular spectra through Trotter error analysis,''
Quantum \textbf{9}, 1630 (2025); arXiv:2311.03719.
```

**Additional recommendation:** Paper 14 text should be amended where it claims the LiH/H2O Pauli counts come from this paper (lines 23, 286, 296, 713, 721, 958, 972, 2199, 2330, 2370, 2404). Best option: add an internal-computation citation (e.g., `\cite{loutey_internal_pauli_2026}` pointing to a Zenodo memo) and reserve `trenev2025` for the Trotter-error-analysis methodology reference at the front of the paper.

---

## F2. Paper 16 `Schwerdtfeger2015` — replacement

**Confidence: HIGH**

**Finding:**
- The Schwerdtfeger–Pašteka–Punnett–Bowman 2015 *Nucl. Phys. A* **944**, 551–577 paper IS real and confirmed via ADS (2015NuPhA.944..551S). Title: "Relativistic and quantum electrodynamic effects in superheavy elements."
- The Pershina 2014 chapter "Theoretical chemistry of the heaviest elements" in Schädel–Shaughnessy 2nd ed. *Chemistry of Superheavy Elements* (Springer 2014, pp. 135–240) ALSO exists.

**Reading Paper 16 context** (`papers/group4_quantum_computing/paper_16_periodicity.tex`):

The `\cite{Schwerdtfeger2015}` usage in Paper 16 should be checked for whether the claim is about (a) **relativistic+QED effects on superheavy elements** (favors Schwerdtfeger 2015 *Nucl. Phys. A* — a focused review article) or (b) **chemistry/predicted-orbital-structure of superheavy elements** (favors Pershina 2014 — a book chapter with periodic-table predictions).

The original (wrong) bibitem was a four-author book chapter "Relativistic effects in superheavy elements" — phrasing suggesting the chapter form. The best single replacement is therefore likely **Pershina 2014** as a structural fit to the bibliographic shape.

**Recommended bibitem texts** (PM should pick based on the actual claim being supported):

Option A (most likely fit, replaces the "book chapter" shape):
```
\bibitem{Pershina2014}
V.~Pershina,
``Theoretical chemistry of the heaviest elements,''
in \emph{The Chemistry of Superheavy Elements}, 2nd ed.,
edited by M.~Sch\"adel and D.~Shaughnessy
(Springer, Berlin, 2014), pp.~135--240.
```

Option B (if the cite supports a relativistic+QED-effects claim):
```
\bibitem{Schwerdtfeger2015}
P.~Schwerdtfeger, L.~F.~Pa\v{s}teka, A.~Punnett, and P.~O.~Bowman,
``Relativistic and quantum electrodynamic effects in superheavy elements,''
Nucl.~Phys.~A \textbf{944}, 551--577 (2015);
doi:10.1016/j.nuclphysa.2015.02.005.
```

**URLs:**
- https://ui.adsabs.harvard.edu/abs/2015NuPhA.944..551S
- https://www.sciencedirect.com/science/article/abs/pii/S0375947415000366
- https://link.springer.com/book/10.1007/978-3-642-37466-1

---

## F3a. Paper 23 `Pachucki2023` — three-photon-exchange paper misattribution

**Confidence: HIGH**

**Finding:**
- The cited PRL 130, 023004 (2023) is actually **"Recoil Corrections to the Energy Levels of Hydrogenic Atoms"** by Pachucki–Yerokhin (verified at journals.aps.org issue listing).
- The actual **"Three-Photon Exchange Nuclear Structure Correction in Hydrogenic Systems"** paper by Pachucki–Patkóš–Yerokhin is **Phys. Rev. A 97, 062511 (2018)**, arXiv:1803.10313.
- There is no published "Three-Photon Exchange" Pachucki–Patkóš–Yerokhin paper in PRL.

**Verdict:** The claimed title is correct (paper exists, is by the named authors), but published in PRA in 2018, not PRL in 2023.

**Recommended bibitem text:**
```
\bibitem{PachuckiPatkosYerokhin2018}
K.~Pachucki, V.~Patk\'o\v{s}, and V.~A.~Yerokhin,
``Three-photon exchange nuclear structure correction in hydrogenic systems,''
Phys.~Rev.~A \textbf{97}, 062511 (2018); arXiv:1803.10313.
```

**Note:** Paper 23 uses `\cite{Pachucki2023}` in two places (line 807 for "leading-order term of the Pachucki--Patkóš--Yerokhin 2023 Foldy--Wouthuysen-reduced two-particle Hamiltonian"; line 1111 for "next-order sign disagreement against Pachucki 2023"). The PM should also check whether the actual supporting claim is the 2018 three-photon-exchange paper, OR a different 2023 paper (e.g., Patkóš–Yerokhin–Pachucki PRL 131, 183001 (2023) "Higher-Order QED Corrections to the Hyperfine Splitting in 3He" — would be a different physics target). Recommend renaming the cite key to `PachuckiPatkosYerokhin2018` so the year mismatch surfaces explicitly.

**URLs:**
- https://journals.aps.org/pra/abstract/10.1103/PhysRevA.97.062511
- https://arxiv.org/abs/1803.10313
- https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.130.023004 (the WRONG paper at the claimed DOI)

---

## F3b. Paper 23 `PachuckiYerokhin2010` — deuterium HFS paper misattribution

**Confidence: MEDIUM**

**Finding:**
- The cited PRL 104, 070403 (2010) is actually **"Fine Structure of Heliumlike Ions and Determination of the Fine Structure Constant"** by Pachucki–Yerokhin (verified at journals.aps.org issue listing).
- No Pachucki–Yerokhin "Theoretical hyperfine splitting in hydrogenic atoms with deuteron nucleus" PRL was found in 2010.
- The closest real paper is **Pachucki 2011** (single author): "Nuclear Structure Corrections in Muonic Deuterium," PRL 106, 193007 (2011). This is the standard reference for the deuteron polarizability +44 ppm correction to hyperfine splitting.

**Verdict:** The claimed title and bibitem do not correspond to any published paper. The most likely intended reference is Pachucki 2011 PRL 106, 193007.

**Recommended bibitem text:**
```
\bibitem{Pachucki2011}
K.~Pachucki,
``Nuclear structure corrections in muonic deuterium,''
Phys.~Rev.~Lett.~\textbf{106}, 193007 (2011);
arXiv:1102.3833; doi:10.1103/PhysRevLett.106.193007.
```

**Note:** Paper 23 uses this cite for the claim "deuteron polarizability ($+44$~ppm, Pachucki--Yerokhin~2010… QCD-internal NN dynamics)". If the +44 ppm number specifically comes from Pachucki 2011, then renaming to `Pachucki2011` with updated authors (single) is correct. The PM should verify the +44 ppm is in Pachucki 2011 and not another paper. Alternative candidate: Kalinowski–Pachucki–Yerokhin 2018 (arXiv:1810.06601) on muonic deuterium — but that's 2018, not 2010.

**URLs:**
- https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.104.070403 (the WRONG paper)
- https://doi.org/10.1103/PhysRevLett.106.193007 (likely intended)
- https://arxiv.org/abs/1810.06601 (alternative)

---

## F3c. Paper 23 `Eides2024` — hyperfine paper misattribution

**Confidence: MEDIUM (best candidate identified, but uncertain)**

**Finding:**
- The cited DOI 10.1016/j.physletb.2024.139049 is actually **"Dynamical non-locality in the near-horizon region of a black hole with quantum time"** by Hadi–Akbarieh (Phys. Lett. B, 2024). Confirmed via SCOAP3 backend.
- There is NO 2024 Eides paper on hydrogen hyperfine splitting in Phys. Lett. B with that DOI.
- **Best candidate replacement: Eides–Shelyuto 2023** "Two-loop corrections to Lamb shift and hyperfine splitting in hydrogen via multi-loop methods," **JHEP 07 (2023) 211**, arXiv:2306.13369. This is the most recent Eides paper on hydrogen hyperfine splitting found in the search.

**Verdict:** The claimed DOI is wrong. The most likely intended reference is Eides–Shelyuto JHEP 2023. The PM should also consider Eides–Shelyuto 2023 *JHEP* 12 (2023) 147 (Wichmann-Kroll potential) as a sibling candidate. Could-not-find a true 2024 Eides hyperfine paper.

**Recommended bibitem text:**
```
\bibitem{EidesShelyuto2023}
M.~I.~Eides and V.~A.~Shelyuto,
``Two-loop corrections to Lamb shift and hyperfine splitting in
hydrogen via multi-loop methods,''
J.~High Energy Phys. \textbf{07} (2023) 211;
arXiv:2306.13369; doi:10.1007/JHEP07(2023)211.
```

**Note:** Paper 23 uses `\cite{Eides2024}` alongside `\cite{EidesGrotchShelyuto2007}` to support the Zemach correction at the operator level (line 675) and as a "scalar substitution" reference (line 882). If the Zemach-correction claim specifically traces to the 2007 review only, the 2024 cite may be excisable; if a more recent Eides-side compilation is desired, JHEP 2023 is the best available replacement.

**URLs:**
- https://scoap3-prod-backend.s3.cern.ch/media/files/88778/10.1016/j.physletb.2024.139049.xml (the WRONG paper)
- https://link.springer.com/article/10.1007/jhep07(2023)211
- https://arxiv.org/abs/2306.13369

---

## F4. `ucp_maps_2024` — Hekkelman–McDonald–vS misattribution

**Confidence: HIGH (could-not-find with confidence)**

**Finding:**
- arXiv:2410.15454 "Gromov-Hausdorff convergence of metric spaces of UCP maps" is by **Tirthankar Bhattacharyya, Ritul Duhan, and Chandan Pradhan** (verified at arxiv.org/abs/2410.15454 v2 Feb 2026). NOT by Hekkelman, McDonald, or van Suijlekom.
- No Hekkelman–McDonald–vS paper on UCP maps was found. The closest real Hekkelman–McDonald paper is **arXiv:2412.00628 (Dec 2024)** "A noncommutative integral on spectrally truncated spectral triples, and a link with quantum ergodicity" by Hekkelman–McDonald (NOT vS). This is about a noncommutative integral, not UCP maps specifically.
- van Suijlekom's solo work arXiv:2409.02773 "A generalization of K-theory to operator systems" is related but not on UCP maps.
- Foundational work: Connes–vS arXiv:2004.14115 (PRX Quantum) "Spectral truncations in noncommutative geometry and operator systems."

**Verdict:** The cite key `ucp_maps_2024` with arXiv:2410.15454 should be corrected to Bhattacharyya–Duhan–Pradhan attribution (the paper exists; only the author list was wrong). If the GeoVac claim is genuinely about Hekkelman–McDonald–vS UCP work, **could-not-find such a paper** — recommend dropping the cite or replacing with Hekkelman–McDonald 2024 (arXiv:2412.00628) if the supporting claim aligns.

**Recommended bibitem text (corrected attribution for arXiv:2410.15454):**
```
\bibitem{bhattacharyya_duhan_pradhan2024}
T.~Bhattacharyya, R.~Duhan, and C.~Pradhan,
``Gromov--Hausdorff convergence of metric spaces of UCP maps,''
arXiv:2410.15454 (2024); v2 February 2025.
```

**Alternative (if the citing claim is about Hekkelman–McDonald spectrally-truncated framework):**
```
\bibitem{hekkelman_mcdonald2024}
E.~Hekkelman and E.~A.~McDonald,
``A noncommutative integral on spectrally truncated spectral triples,
and a link with quantum ergodicity,''
J.~Funct.~Anal. (2025); arXiv:2412.00628.
```

**URLs:**
- https://arxiv.org/abs/2410.15454
- https://arxiv.org/abs/2412.00628
- https://www.sciencedirect.com/science/article/pii/S0022123625003362

---

## F5. `mondino_samann2024`, `mondino_samann2025`, `minguzzi_suhr2024`, `sakovich_sormani2024` — XXXXX placeholders

**Confidence: HIGH**

### F5a. `mondino_samann2025` (Lorentzian GH convergence)

**Finding:** Mondino–Sämann 2025 "Lorentzian Gromov–Hausdorff convergence and pre-compactness" is **arXiv:2504.10380** (v1 April 14, 2025; v4 Dec 9, 2025).

**Recommended bibitem text:**
```
\bibitem{mondino_samann2025}
A.~Mondino and C.~S\"amann,
``Lorentzian Gromov--Hausdorff convergence and pre-compactness,''
arXiv:2504.10380 (2025).
```

### F5b. `mondino_samann2024`

**Finding:** **Could-not-find with confidence** a single canonical Mondino–Sämann 2024 paper on Lorentzian pre-length spaces. The closest candidates:
- Cavalletti–Mondino "Optimal transport in Lorentzian synthetic spaces…" (arXiv:2004.08934, 2020 — not 2024, not Sämann)
- Mondino review papers
- The 2025 paper above is most likely the intended Mondino–Sämann reference; the 2024 entry may be a duplicate or stale-named version that should be merged into `mondino_samann2025`.

**Recommendation:** Consolidate `mondino_samann2024` and `mondino_samann2025` into a single bibitem using arXiv:2504.10380, OR identify the specific 2024 Mondino–Sämann work being cited. Honest verdict: **could-not-find** a distinct 2024 paper.

### F5c. `minguzzi_suhr2024` (Lorentzian metric spaces & GH-convergence)

**Finding:** Minguzzi–Suhr 2024 "Lorentzian metric spaces and their Gromov–Hausdorff convergence" is **arXiv:2209.14384**, published in **Lett. Math. Phys. 114 (2024), article 73**. Verified.

**Recommended bibitem text:**
```
\bibitem{minguzzi_suhr2024}
E.~Minguzzi and S.~Suhr,
``Lorentzian metric spaces and their Gromov--Hausdorff convergence,''
Lett.~Math.~Phys. \textbf{114}, article 73 (2024);
arXiv:2209.14384.
```

### F5d. `sakovich_sormani2024` (Various distances between space-times)

**Finding:** Sakovich–Sormani 2024 "Introducing Various Notions of Distances between Space-Times" is **arXiv:2410.16800** (v1 Oct 22, 2024; v3 Oct 15, 2025).

**Recommended bibitem text:**
```
\bibitem{sakovich_sormani2024}
A.~Sakovich and C.~Sormani,
``Introducing various notions of distances between space-times,''
arXiv:2410.16800 (2024).
```

---

## F6. `mondino_samann2024` mis-arXiv-ID (Papers 45, 46)

**Confidence: HIGH**

**Finding:** arXiv:2209.14384 IS by Minguzzi–Suhr ("Lorentzian metric spaces and their Gromov–Hausdorff convergence," LMP 114, 73, 2024), NOT Mondino–Sämann. The misattribution is confirmed.

**Resolution:** In Papers 45 and 46, the bibitem currently labeled `mondino_samann2024` and pointing to arXiv:2209.14384 should be:
1. Renamed to `minguzzi_suhr2024` (matching the correct authors)
2. The arXiv ID 2209.14384 is correct for Minguzzi–Suhr LMP 2024.
3. Any text referring to "Mondino–Sämann" in connection with arXiv:2209.14384 should be changed to "Minguzzi–Suhr."

The actual Mondino–Sämann work is **arXiv:2504.10380** (per F5a). If both citations are needed, both bibitems should exist with their correct attributions.

**Recommended bibitem text:** Same as F5c above.

---

## F7. `che_perales_sormani2025`

**Confidence: HIGH**

**Finding:** The paper is **"Gromov's Compactness Theorem for the Intrinsic Timed-Hausdorff Distance"** by **Mauricio Che, Raquel Perales, and Christina Sormani**, arXiv:2510.13069 (v1 Oct 15, 2025; v3 Feb 26, 2026), published in **Differential Geometry and its Applications Vol. 103 (2026)**.

The first author's first name is **Mauricio** (not "R."). The journal is **Differential Geometry and its Applications**, NOT Adv. Math.

**Recommended bibitem text:**
```
\bibitem{che_perales_sormani2025}
M.~Che, R.~Perales, and C.~Sormani,
``Gromov's compactness theorem for the intrinsic timed-Hausdorff
distance,''
Differ.~Geom.~Appl. \textbf{103} (2026);
arXiv:2510.13069 (2025).
```

**URLs:**
- https://arxiv.org/abs/2510.13069
- https://www.sciencedirect.com/journal/differential-geometry-and-its-applications (Vol. 103)

---

## Summary Table

| ID | Confidence | Resolution |
|:---|:-----------|:-----------|
| F1 | HIGH (mixed) | Cited paper exists but doesn't contain claimed counts; title incorrect (missing "vibrational"); Pauli counts are likely GeoVac-internal. Recommend: fix title, add internal-computation cite for the counts. |
| F2 | HIGH | Both Pershina 2014 (book chapter) and Schwerdtfeger 2015 (Nucl. Phys. A 944, 551) exist. Pershina better matches the bibitem's "book chapter" shape. |
| F3a | HIGH | "Three-Photon Exchange" Pachucki–Patkóš–Yerokhin is **PRA 97, 062511 (2018)**, NOT PRL 130, 023004 (2023). |
| F3b | MEDIUM | PRL 104, 070403 is Pachucki–Yerokhin heliumlike fine structure, not deuterium HFS. Likely intended: **Pachucki PRL 106, 193007 (2011)** "Nuclear Structure Corrections in Muonic Deuterium." |
| F3c | MEDIUM | physletb.2024.139049 is a Hadi-Akbarieh black-hole paper. Best Eides replacement candidate: **Eides–Shelyuto JHEP 07 (2023) 211**, arXiv:2306.13369. |
| F4 | HIGH (could-not-find Hekkelman-McDonald-vS UCP work) | arXiv:2410.15454 IS by Bhattacharyya–Duhan–Pradhan. Closest real Hekkelman–McDonald: arXiv:2412.00628 (no vS, no UCP focus). |
| F5 | HIGH (3 of 4) | mondino_samann2025 = arXiv:2504.10380; minguzzi_suhr2024 = arXiv:2209.14384, LMP 114:73 (2024); sakovich_sormani2024 = arXiv:2410.16800; mondino_samann2024 could-not-find as a distinct work. |
| F6 | HIGH | arXiv:2209.14384 is Minguzzi–Suhr (verified), NOT Mondino–Sämann. Rename bibitem. |
| F7 | HIGH | M. (Mauricio) Che, R. Perales, C. Sormani — arXiv:2510.13069 — Differ. Geom. Appl. 103 (2026), NOT Adv. Math. |
