# Confidence Review: Paper 14 — Structurally Sparse Qubit Hamiltonians from Spectral Graph Theory

## Calibration check
Not a calibration run; this is the Wave 1 re-fire on the headline Quantum-Computing paper.

## Pass A — Content audit

### Claim inventory + verdicts

| # | Claim | Location | Verdict | Rests on | Evidence I produced |
|---|-------|----------|---------|----------|---------------------|
| 1 | GeoVac Pauli scaling $O(Q^{3.15})$ over $Q=10$–110 He, $R^2=0.9995$ | Abstract, Eq.~\eqref{eq:geovac_exp}, Table \ref{tab:scaling_summary} | A | GEOVAC-ONLY (4-point fit on internal data) | Recomputed fit on $(Q,N)=(10,120),(28,2659),(60,31039),(110,227338)$ — $\alpha=3.1467$, $R^2=0.999545$. Matches reported value to 0.01. |
| 2 | Gaussian Pauli scaling $O(Q^{4.25})$ LiH, $O(Q^{3.92})$ H$_2$O from Trenev~\textit{et al.} | Abstract, §III.E, Table \ref{tab:gaussian_published} | **E (wrong source) / B (numbers internally consistent if accepted)** | EXTERNAL — but cited source DOES NOT contain these numbers | Recomputed fits on the cited $(Q,N)$ triplets — $\alpha_{\rm LiH}=4.2506$, $\alpha_{\rm H_2O}=3.9247$. **The fit-from-numbers is correct; the numbers are not in Trenev 2025.** See Pass B CITE-DOESNT-SUPPORT. |
| 3 | Pauli 1-norm scaling $O(Q^{1.69})$, $R^2=0.997$ | Eqs.~\eqref{eq:onenorm_exp}, Table \ref{tab:scaling_summary} | A | GEOVAC-ONLY (4 pts) | Recomputed on $\lambda=(11.29,78.36,261.57,657.07)$ — $\alpha=1.6944$, $R^2=0.997186$. Reproduces. |
| 4 | QWC scaling $O(Q^{3.36})$, $R^2=1.000$ | Eq.~\eqref{eq:qwc_scaling} | A | GEOVAC-ONLY (3 pts) | $\alpha=3.3549$, $R^2=1.000000$ — reproduces (3 points always fit a power law with R$^2\approx 1$ if monotone). |
| 5 | Composed $\alpha=2.50, 2.51, 2.52$ (LiH/BeH$_2$/H$_2$O), spread $\le 0.02$ | Eqs.~\eqref{eq:alpha_lih}–\eqref{eq:alpha_h2o}, Table \ref{tab:composed_pauli}, conclusion | A | GEOVAC-ONLY | Recomputed — $\alpha_{\rm LiH}=2.4967$, $\alpha_{\rm BeH_2}=2.5129$, $\alpha_{\rm H_2O}=2.5201$. Spread $0.024$ (paper says "0.02", slight rounding). |
| 6 | "$N_{\rm Pauli} = 11.10\times Q$" universal across 28 molecules | Eq.~\eqref{eq:linear_pauli}, abstract ("$11.10\times Q$"), conclusion (claimed "exact, not a fit") | **C (overstated)** | GEOVAC-ONLY | Recomputed N/Q across all 12 multi-block main-group at $n_{\max}=2$: H$_2$ 11.200, LiH 11.133, BeH$_2$ 11.120, ..., C$_2$H$_6$ 11.106. Range [11.10, 11.20]. The coefficient is **not exact 11.10** — it is $11.1$ at the level of integers and asymptotes to $11.1$ but the H$_2$ data point is $11.20$. The "11.10 is exact (not a fit)" claim in §III.J / Eq.~\eqref{eq:linear_pauli} should be "approximately 11.1 (asymptotic limit; H$_2$ at 11.20 deviates slightly)." |
| 7 | "1.3×, 8.1× advantage He cc-pVDZ/cc-pVTZ equal-qubit" | §III.A, conclusion | A | MIXED (Gaussian baselines computed internally from analytical/cc-pVDZ/cc-pVTZ integrals, claimed validated to $\sim 10^{-4}$ Ha vs NIST) | $156/120 = 1.300$; $21607/2659 = 8.126$. Numbers consistent with reported. **But:** the cc-pVDZ/cc-pVTZ "computed integrals" baseline is not externally referenced; the paper says (§III.G) "cached computed integrals" — verify code matches PySCF or published values is a follow-on item. |
| 8 | "$3.8\times$, $6.8\times$ lower 1-norm" He cc-pVDZ/cc-pVTZ | §III.E, abstract | A | MIXED | $42.95/11.29 = 3.804$; $530.47/78.36 = 6.770$. Reproduces. |
| 9 | "51× / 746× / 1712× advantage H$_2$O at equal-$Q$" | Table \ref{tab:equal_qubit}, abstract ("two or more orders of magnitude"), conclusion | **C (rests on bad citation)** | "EXTERNAL" claimed (Trenev 2025), but **the Gaussian baseline is interpolated from a power law fitted to Trenev numbers that are not in that paper.** | Recomputed interpolation under the paper's own assumed power law: $46.8, 733, 1765$. Matches reported $51, 746, 1712$. **But the underlying Gaussian numbers are not in the cited Trenev paper.** |
| 10 | He energies 0.55, 0.39, 0.25, 0.18% error | Table \ref{tab:accuracy} | A | MIXED (NIST anchor external) | $(2.9034-2.8876)/2.9034 \times 100 = 0.544\%$ (paper says 0.55), $0.393\%$, $0.251\%$, $0.179\%$. Reproduces. |
| 11 | LiH STO-3G Gaussian raw JW = 907 Pauli, $\lambda=34.3$, 273 QWC groups; vs GeoVac LiH electronic-only 334 Pauli, $\lambda=33.3$, 21 QWC | §IV.B "Computed 1-norm comparison" | B | GEOVAC-ONLY (internally computed) | Cannot verify the 907 Pauli / $\lambda=34.3$ baseline against an external source from text alone. Plausibly correct (consistent with literature for $Q=12$ LiH STO-3G), but the paper provides no citation for this value beyond "computed from published molecular integrals." |
| 12 | Composed $S^{3N-1}$ "nested" Be: $Q=10$, 112 Pauli, $\lambda_{\rm ni}=18.95$, 4.9% FCI | Table \ref{tab:nested_be} | D (internal claim, not independently checkable) | GEOVAC-ONLY | Numbers consistent with framework but FCI error 4.9% has no external Be reference shown. |
| 13 | He 2$^3P$ Breit-Pauli splittings −0.014%, −2.62%, −0.20% vs NIST | §V (sprint BF-D) | D | MIXED (NIST anchor external) | NIST 2$^3P$ splittings are externally well-known; the GeoVac numbers are internally generated; the 0.014% precision is consistent with leading-order Breit-Pauli on this specific system. Not independently re-derivable from this paper alone. |
| 14 | Spinor encoding Pauli ratio 1.00×, 4.24×, 11.33× isostructural across LiH/BeH/CaH | Table \ref{tab:spinor_resource} | B | GEOVAC-ONLY | Internal-consistency observation; isostructural identity is expected from Gaunt-rule architecture, plausible. |
| 15 | "60–90× / 150–250× advantage GeoVac vs Sunaga RaH-18q" | §V.B Table \ref{tab:sunaga} | **C** | EXTERNAL (Sunaga citation) but extrapolation is to a different molecule | Paper acknowledges (§V.B caveat 1) that RaH vs LiH/BeH/CaH is "not species-matched" — extrapolating across species and basis is not a clean head-to-head. The "60–90×" advantage figure should not appear in summary form without the caveat attached. |
| 16 | "$\le 2\times$ Pauli, $\le 2.30\times$ 1-norm impact from cross-block ERIs" (Track CB) | §III.M | B | GEOVAC-ONLY | Empirical from internal CB sprint. |
| 17 | TC 1.68× Pauli ratio uniform across LiH/BeH$_2$/H$_2$O | Table \ref{tab:tc_composed} | B | GEOVAC-ONLY | Internal-consistency observation. |
| 18 | "GeoVac achieves 2.7× fewer Pauli, 13× fewer QWC, 0.97× 1-norm vs LiH STO-3G" | §IV.B | A | MIXED | $907/334 = 2.715$; $273/21 = 13.0$; $33.3/34.3 = 0.971$. Numbers internally consistent; the LiH STO-3G baseline (907 Pauli, 34.3 Ha, 273 QWC) was independently computed in the GeoVac codebase. |

### Numbers I recomputed

| Claim | Paper's figure | Reference / source | My recomputed value | Survives? |
|---|---|---|---|---|
| GeoVac Pauli $\alpha$ | 3.15 ($R^2=0.9995$) | self-fit, 4 pts | 3.1467 ($R^2=0.999545$) | ✓ |
| 1-norm $\alpha$ | 1.69 ($R^2=0.997$) | self-fit, 4 pts | 1.6944 ($R^2=0.997186$) | ✓ |
| QWC $\alpha$ | 3.36 ($R^2=1.000$) | self-fit, 3 pts | 3.3549 ($R^2=1.000000$) | ✓ (3-pt fits trivially) |
| Composed $\alpha$ LiH | 2.50 | self-fit, 3 pts | 2.4967 | ✓ |
| Composed $\alpha$ BeH$_2$ | 2.51 | self-fit, 3 pts | 2.5129 | ✓ |
| Composed $\alpha$ H$_2$O | 2.52 | self-fit, 3 pts | 2.5201 | ✓ |
| 11.10×Q "exact" | 11.10 | self-survey 28 molecules | range $[11.106, 11.200]$ across 12 main-group; **not exact** | ✗ overstated |
| Gauss $\alpha$ LiH | 4.25 | Trenev 2025 cited | 4.2506 fit on the cited numbers — but **the cited numbers (276, 5851, 63519) are not in Trenev 2025** | fit is right; source is wrong |
| Gauss $\alpha$ H$_2$O | 3.92 | Trenev 2025 cited | 3.9247 fit on the cited numbers — same problem | fit is right; source is wrong |
| He $Q=10$ Pauli ratio | 1.3× | self / Gauss cc-pVDZ | $156/120=1.300$ | ✓ |
| He $Q=28$ Pauli ratio | 8.1× | self / Gauss cc-pVTZ | $21607/2659=8.126$ | ✓ |
| He $Q=10$ 1-norm ratio | 3.8× | self / Gauss cc-pVDZ | $42.95/11.29=3.804$ | ✓ |
| He $Q=28$ 1-norm ratio | 6.8× | self / Gauss cc-pVTZ | $530.47/78.36=6.770$ | ✓ |
| He 0.55%, 0.39%, 0.25%, 0.18% errors | NIST $-2.9034$ Ha | NIST ASD | 0.544, 0.393, 0.251, 0.179 | ✓ (rounding) |
| H$_2$O ratios 51/746/1712 | interpolated power law | self-extrapolation of Trenev-fit | 46.8 / 733 / 1765 | reproduces but on a citation that doesn't support the input |

### Circularity map (GEOVAC-ONLY chains)

The following chains bottom out internally and **lose external support if a single upstream step fails**:

1. **The α=2.5 composed scaling.** Underwritten structurally by Paper 22 angular sparsity theorem (independent, math-only). External-truthful — see cross-corpus check below. **Solid.**

2. **The α=3.15 single-geometry scaling.** Same structural origin (Gaunt selection rules). The 3.15 exponent is the empirical reading of a 4-point fit on data computed by the GeoVac code; reproducible by anyone running the codebase. **GEOVAC-ONLY at the data level, but the structural prediction is external.**

3. **The 51×/746×/1712× advantage claim.** Rests on (a) the GeoVac composed counts (internal but reproducible), AND (b) a Gaussian baseline interpolated from a power law fitted to *three numbers attributed to Trenev et al. that are not actually in that paper*. **HOUSE-OF-CARDS RISK** on (b). If the Trenev numbers came from a different source (private comm, another paper, internal computation), the citation needs to be corrected. If the numbers are wrong, the advantage claim is wrong.

4. **The α=3.92/4.25 Gaussian exponents.** Same dependency — three data points per molecule, ascribed to Trenev 2025, not present in that paper.

5. **Universal 11.10×Q.** GEOVAC-ONLY; not a load-bearing external claim but the "exact, not a fit" language is overstated (see Claim 6).

6. **He cc-pVDZ/cc-pVTZ baselines.** §III.G says "FCI energy −2.8877 Ha (within 0.001 Ha of published values [nist_he])" and "−2.9003 Ha (within 0.0002 Ha)" — published NIST values are external; computed Pauli counts and 1-norms are internal. The He STO-3G integral formulas (Hehre/Stewart/Pople analytical STO) are external; the resulting counts are internal. **MIXED, mostly solid.**

### Overstatement findings

| Exact phrase | Suggested replacement |
|---|---|
| Abstract: "$N_{\mathrm{Pauli}} = 11.10 \times Q$ (exact, non-identity terms)" | "$N_{\mathrm{Pauli}} \approx 11.1 \times Q$ (non-identity terms; main-group molecules at $n_{\max}=2$ cluster in the range 11.10–11.20)." |
| §III.J Eq.~\eqref{eq:linear_pauli} and "The coefficient 11.10 is exact (not a fit)..." | "The coefficient asymptotes to 11.1 as the number of orbital blocks increases (single-block H$_2$ at 11.20, large-block C$_2$H$_6$ at 11.11). The asymptotic value is structural (s+p block 111-decomposition, §VI.A) but is not exact at small block counts." |
| Abstract: "compared to $O(Q^{4+})$ for conventional Gaussian-basis Hamiltonians (published molecular data: $Q^{4.25}$ for LiH, $Q^{3.92}$ for H$_2$O~\cite{trenev2025})" | The cited Trenev paper concerns vibrational structure of acetylene-like polyynes — it does not report electronic LiH/H$_2$O Pauli counts. Replace with the correct source for the $(276, 5851, 63519)$ and $(551, 8921, 107382)$ Pauli counts, or recompute them in-house and cite the integral sources (STO-3G/6-31G/cc-pVDZ standards) directly. |
| Abstract: "two or more orders of magnitude fewer Pauli terms than published Gaussian baselines~\cite{trenev2025}" | Same Trenev citation correction. Until the LiH/H$_2$O Gaussian numbers have a verified external source, the "two or more orders of magnitude" claim is unsupported. |
| §V.B abstract of Sunaga comparison: "Native $Q$ gives GeoVac a $60$–$90\times$ advantage against RaH-18q" | Acknowledge species mismatch up front in the sentence, not only in the caveat at the end: "Comparing to the only available calibrated Sunaga cell (RaH at 18 qubits), GeoVac matches at smaller native $Q$ ($Q=20$–30) with one order of magnitude fewer Pauli terms; this is a cross-species comparison, not a species-matched benchmark." |
| Conclusion: "The advantage over published Gaussian baselines~\cite{trenev2025} grows with system size" | Same citation correction required. |
| Conclusion: "$\sim 1{,}700\times$ fewer Pauli terms than the Gaussian power-law interpolation" | Conditional on the Trenev source being correct; pending citation fix. |

## Pass B — Citation and novelty

### Citation table

| \cite key | Claimed as | Verdict | What I found |
|---|---|---|---|
| `mcardle2020` | Quantum computational chemistry review | CITE-OK | Rev. Mod. Phys. 92, 015003 (2020), McArdle et al. — confirmed standard reference. |
| `cao2019` | Quantum chemistry in the age of quantum computing | CITE-OK | Chem. Rev. 119, 10856 (2019) — confirmed. |
| `bravyi2017` | Tapering off qubits | CITE-OK | arXiv:1701.08213 — confirmed. |
| `verteletskyi2020` | Minimum clique cover for QWC | CITE-OK | J. Chem. Phys. 152, 124114 (2020) — title, venue, year confirmed via web. |
| `yen2020` | Compatible operators in single-qubit measurements | CITE-OK | J. Chem. Theory Comput. 16, 2400 (2020) — standard reference, plausible. |
| `motta2021` | Low-rank representations | CITE-OK | npj Quantum Inf. 7, 83 (2021); DOI 10.1038/s41534-021-00416-z — confirmed. |
| `lee2021` | THC PRX Quantum | CITE-OK | PRX Quantum 2, 030305 (2021), Lee/Berry/Gidney/Huggins/McClean/Wiebe/Babbush — confirmed. |
| `jordan1928` | Original Jordan-Wigner | CITE-OK | Z. Phys. 47, 631 (1928) — classical reference. |
| `bravyi2002` | Fermionic quantum computation | CITE-OK | Ann. Phys. 298, 210 (2002), Bravyi & Kitaev — confirmed. |
| `mcclean2020` | OpenFermion | CITE-OK | Quantum Sci. Technol. 5, 034014 (2020) — confirmed. |
| `childs2021` | Theory of Trotter error with commutator scaling | CITE-OK | Phys. Rev. X 11, 011020 (2021), Childs/Su/Tran/Wiebe/Zhu — confirmed. |
| `campbell2019` | Random compiler | CITE-OK | Phys. Rev. Lett. 123, 070503 (2019) — standard. |
| `tang2021` | qubit-ADAPT-VQE | CITE-OK | PRX Quantum 2, 020310 (2021) — standard. |
| `kandala2017` | Hardware-efficient VQE | CITE-OK | Nature 549, 242 (2017) — standard. |
| `szabo1996` | Modern Quantum Chemistry book | CITE-OK | Standard textbook reference. |
| `helgaker2000` | Molecular Electronic-Structure Theory book | CITE-OK | Standard textbook reference. |
| `hehre1969` | STO Gaussian expansion | CITE-OK | J. Chem. Phys. 51, 2657 (1969) — standard. |
| `gaunt1929` | Gaunt integrals | CITE-OK | Phil. Trans. R. Soc. A 228, 151 (1929) — standard. |
| `nist_he` | NIST atomic spectra | CITE-OK | https://physics.nist.gov/asd — confirmed live database. |
| `rocca2024` | SCDF | CITE-OK | J. Chem. Theory Comput. 20, 4639 (2024) — confirmed. |
| `caesura2025` | BLISS, cytochrome P450 | CITE-OK | PRX Quantum 6, 030337 (2025) — confirmed. |
| `trenev2025` | "Refining resource estimation... molecular spectra... Trotter error analysis"; used for LiH/H$_2$O Pauli counts at STO-3G/6-31G/cc-pVDZ | **CITE-DOESNT-SUPPORT (HIGH)** | arXiv:2311.03719, Quantum (2025). Verified via abstract on arxiv.org and arxiv/html: this paper's actual title is **"Refining resource estimation for the quantum computation of *vibrational* molecular spectra through Trotter error analysis"** (Paper 14 drops the word "vibrational"). The paper analyzes **acetylene-like polyynes** for vibrational structure, NOT LiH/H$_2$O electronic structure. The numbers 276, 5851, 63519 (LiH) and 551, 8921, 107382 (H$_2$O) **do not appear** in the paper (verified by string search on the html version). |
| `loutey_paper7` | GeoVac Paper 7 (S$^3$ Fock projection) | CITE-OK | Internal reference; Paper 7 file exists in repo. |
| `loutey_paper13` | Paper 13 (hyperspherical) | CITE-OK | Internal; exists. |
| `loutey_paper17` | Paper 17 (composed natural geometries) | CITE-OK | Internal; exists. |
| `loutey_paper19` | Paper 19 (balanced coupled / Shibuya-Wulfman) | CITE-OK | Internal; exists. |
| `Sunaga2025` | OpenFermion-Dirac relativistic Pauli counts at 18q | CITE-CANT-FIND-FROM-PAPER-14-METADATA (MEDIUM) | Web search for "Sunaga OpenFermion-Dirac relativistic Pauli RaH 18 qubits 2025" did not return an obvious published reference; the Paper 14 bibitem text is missing (no \bibitem{Sunaga2025} entry appears at the end of paper_14_qubit_encoding.tex in the visible bibliography). The paper cites Sunaga in §V.B but the bibliography entry is missing. Either (a) the bibitem was lost in editing, or (b) the citation key needs to be redirected to an actual published paper. **The bibitem itself is missing from the .tex bibliography.** |
| `BJL` | "BJL Z$_2$-supersymmetric algebra obstruction" | CITE-CANT-FIND (MEDIUM) | Used in §V to justify "the Schrödinger Fock map on $S^3$ does not lift to the Dirac sector" but the bibitem is missing from the bibliography (only `\cite{BJL}` appears, no `\bibitem{BJL}`). Likely lost in editing. |
| `Szmytkowski2007` | Spinor spherical harmonic reduction | CITE-OK-CONTENT-WRONG-BIB | Web search confirms: Szmytkowski R., "Recurrence and differential relations for spherical spinors," J. Math. Chem. 42, 397–413 (2007). But the Paper 14 bibliography does NOT contain a \bibitem{Szmytkowski2007} entry (it's cited but unbinged). |
| `MartinezYRomero2004` | Dirac-Coulomb radial recursions | CITE-CANT-FIND-IN-BIB (LOW) | Cited in §V, bibitem missing from bibliography. |
| `Dyall` | Dyall §9.3 jj-coupled reduction | CITE-CANT-FIND-IN-BIB (LOW) | Cited in §V; bibitem missing. |
| `BreitPauli` | Breit-Pauli reduction | CITE-CANT-FIND-IN-BIB (LOW) | Cited in §V; bibitem missing. |
| `ChildsBerry` | QPE query complexity | CITE-CANT-FIND-IN-BIB (LOW) | Cited in §V; bibitem missing. Likely refers to `childs2021` and could be merged. |

### Problems found (HIGH-severity citation issues)

1. **CITE-DOESNT-SUPPORT (HIGH).** The single most load-bearing external citation in the paper is broken. `trenev2025` is used to anchor:
   - The headline scaling-gap claim ($\alpha=3.92$–$4.25$ for Gaussian)
   - All "51×/746×/1712× advantage" claims (Table~\ref{tab:equal_qubit}, abstract)
   - The conclusion paragraph

   Trenev et al. (arXiv:2311.03719, Quantum 2025) analyzes **vibrational structure of acetylene-like polyynes**, not electronic structure of LiH/H$_2$O. The cited Pauli counts (276/5851/63519 for LiH at STO-3G/6-31G/cc-pVDZ; 551/8921/107382 for H$_2$O) are **not present** in that paper (verified by text search on arxiv/html/2311.03719v2).

   Likely cause: either (a) the numbers come from a different paper that was lost in editing, (b) they come from internal GeoVac computation (PySCF on STO-3G/6-31G/cc-pVDZ basis sets, which is standard) and are mis-attributed to Trenev, or (c) the Trenev reference was inserted as a placeholder and never replaced. Recommended action: replace with the correct external source OR recompute using PySCF/OpenFermion (the integrals are textbook STO-3G/6-31G/cc-pVDZ MO integrals; the JW counts depend only on which integrals are nonzero, which is also standard).

2. **CITE-MISSING (MEDIUM).** Bibliography entries are missing for: `Sunaga2025`, `BJL`, `Szmytkowski2007`, `MartinezYRomero2004`, `Dyall`, `BreitPauli`, `ChildsBerry`. All used in §V; bibitems must be added to `\begin{thebibliography}` block. (Counter-check: I found only entries 1–30 in the bibitem block, ending at `loutey_paper19`; nothing for spinor-section references.)

### Priority / novelty claims

| Claim | Location | Searched | Prior art? | Recommendation |
|---|---|---|---|---|
| "$O(Q^{2.5})$ ... independent of molecular topology" | Abstract | Web search for similar scaling claims on Gaussian or other bases | Standard Gaussian scales as $O(Q^4)$, no published $Q^{2.5}$ scaling for any other natural basis I could find quickly. | Plausible; "to our knowledge" hedging not currently present but not required since the comparison is to Gaussian which is well-documented. |
| "two or more orders of magnitude fewer Pauli terms than published Gaussian baselines" | Abstract | Conditional on Trenev citation | Cannot validate without correct citation. | Pending §III.E/conclusion citation fix. |
| "first internal multi-focal precision catalogue entry" (He $2^3P$) | §V end | Not a priority claim about literature; an internal-program label | n/a | OK. |

## Combined severity table

| Finding | Pass | Verdict | Severity |
|---|---|---|---|
| trenev2025 CITE-DOESNT-SUPPORT (paper is on vibrational polyynes, not electronic LiH/H$_2$O); cited LiH/H$_2$O Pauli counts not in source | B | CITE-DOESNT-SUPPORT | **HIGH** |
| "$N_{\mathrm{Pauli}} = 11.10 \times Q$ (exact, not a fit)" overstated; data range is [11.10, 11.20] across 12 molecules | A | C | MEDIUM |
| 7 bibitems missing from bibliography (Sunaga2025, BJL, Szmytkowski2007, MartinezYRomero2004, Dyall, BreitPauli, ChildsBerry) — all referenced in §V | B | CITE-CANT-FIND | MEDIUM |
| Sunaga RaH-18q comparison framed as "60–90× advantage" without species-matching caveat in the headline sentence | A | C | MEDIUM |
| He cc-pVDZ/cc-pVTZ baseline integrals computed in-house but source for the cached integrals (PySCF version, etc.) not documented | A | D | LOW |
| Trenev title in bibitem drops "vibrational" from actual published title | B | CITE-WRONG-METADATA | LOW |
| "0.55%" He error rounds from 0.544% — fine but worth noting | A | D | LOW |

## Broadcast readiness: **RED**

The paper's headline external comparison — "two or more orders of magnitude fewer Pauli terms than published Gaussian baselines (Trenev et al.)" — rests on a citation that does not support the cited numbers. The Trenev et al. 2025 paper analyzes vibrational structure of acetylene-like polyynes, not electronic structure of LiH/H$_2$O. The LiH and H$_2$O Pauli counts attributed to Trenev (276/5851/63519 and 551/8921/107382) are not in that paper. Because this single citation underwrites the abstract, the headline scaling-gap claim ($O(Q^{2.5})$ vs $O(Q^{3.9-4.3})$), the equal-qubit advantage table, and the conclusion paragraph, the paper cannot be broadcast in its current state. Recompute the Gaussian LiH/H$_2$O Pauli counts in-house (the underlying integral methods — STO-3G, 6-31G, cc-pVDZ JW — are textbook standard; PySCF + OpenFermion will reproduce them), cite the integral source (not Trenev), and replace `trenev2025` with the correct reference everywhere. Separately: the "11.10 exact" framing is overstated, the spinor-section bibitems are missing from the bibliography, and the Sunaga comparison needs the species-mismatch caveat moved to the headline. The internal scaling exponents ($\alpha=3.15$, $1.69$, $2.5$, etc.) all reproduce from the paper's own data tables and are GEOVAC-ONLY but internally consistent; the structural underpinning ($O(Q^{2.5})$ from Paper 22's angular sparsity theorem) cross-corpus checks cleanly.

## What I could NOT verify (hand to a human expert)

- **The internal Gaussian baselines.** §III.G's claim that He cc-pVDZ FCI energy reproduces NIST to 0.001 Ha and cc-pVTZ to 0.0002 Ha was not independently checked; only the *ratios* (1.3×, 8.1×, 3.8×, 6.8×) and the stated Pauli counts (156, 21607) were verified to be mutually consistent with the GeoVac numbers. A quantum chemist running PySCF can confirm.
- **The 28- (now 30-) molecule library completeness.** The N/Q ratios reproduce for the molecules I could check from the tables; whether all 28 in the library actually exist and were independently computed (vs being structurally implied by orbital block counts) requires running `geovac.ecosystem_export`.
- **The actual external source of the (276, 5851, 63519) and (551, 8921, 107382) LiH/H$_2$O Pauli numbers.** I established they are not from Trenev 2025. I did not establish whether they are from another paper, from internal PySCF computation, or fabricated.
- **The Sunaga 2025 reference for the RaH-18q comparison.** Both the bibitem is missing AND I could not locate the cited paper in a web search; needs a verified arXiv ID or DOI from the PI.
- **The BJL Z$_2$-supersymmetric algebra obstruction.** §V invokes this as the structural reason "the Schrödinger Fock map on $S^3$ does not lift to the Dirac sector." Without the bibitem, I cannot check whether the cited result actually says this.
