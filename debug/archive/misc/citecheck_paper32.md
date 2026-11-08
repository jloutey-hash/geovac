# Citation & Novelty Check: Paper 32 — *The GeoVac Spectral Triple: A Synthesis Paper Locking the Framework into the Marcolli–van Suijlekom Gauge-Network Lineage*

**Verdict: YELLOW** — the lineage citations exist and the load-bearing arXiv IDs point to the right works, but two bibliography entries have wrong metadata that would be embarrassing in front of an expert (`paschke_sitarz2000` year tag, `cacic2009` volume/pages/year), and one entry has the title swapped to a completely different paper (`bizi_brouder_besnard2018`). The load-bearing "Yang–Mills without Higgs" lineage claim is correctly attributed to `perez_sanchez2025` (arXiv:2508.17338) but is **mis-bundled with `perez_sanchez2024` (arXiv:2401.03705)**, which actually concludes Yang–Mills–HIGGS. This bundling appears in five places in the paper. The α-priority claim is appropriately hedged.

## Citation table (Wave 1: math.OA load-bearing)

| \cite key | claimed as | verdict | what I found |
|:----------|:-----------|:--------|:-------------|
| `marcolli_vs2014` | Marcolli, van Suijlekom, "Gauge networks in noncommutative geometry," J. Geom. Phys. 75, 71–91 (2014), arXiv:1301.3480 | **CITE-OK** | arXiv:1301.3480 exists; authors, title, journal, year all correct; abstract confirms it constructs finite spectral triples on quiver vertices with edges carrying connection data and derives Wilson lattice gauge action via spectral action. **However**, its continuum limit is Yang–Mills–Higgs (not "without Higgs"). |
| `perez_sanchez2024` | C.I. Perez-Sanchez, "Bratteli networks and the Spectral Action on quivers," arXiv:2401.03705 (2024) | **CITE-DOESNT-SUPPORT** | The paper exists with that title/author/year. **But it does NOT correct Marcolli-vS to Yang–Mills without Higgs.** It explicitly acknowledges Marcolli-vS ("paved the way... in independent directions") and derives "Yang–Mills–Higgs theory on flat space as a smooth limit." Paper 32 bundles this with perez_sanchez2025 in five places as "the Perez-Sanchez correction that the continuum limit is Yang–Mills without Higgs" — that claim is owned by perez_sanchez2025 alone. |
| `perez_sanchez2025` | C.I. Perez-Sanchez, "Comment on 'Gauge networks in noncommutative geometry'," arXiv:2508.17338 (2025) | **CITE-OK** | arXiv:2508.17338 exists; abstract confirms it is explicitly a comment on Marcolli-vS 2014 (J. Geom. Phys. 75:71–91, 2014) and explicitly demonstrates "the continuum limit of this theory is the Yang–Mills action functional, without a Higgs scalar." Load-bearing claim is supported. |
| `connes_vs2021` | Connes, van Suijlekom, "Spectral truncations in noncommutative geometry and operator systems," Commun. Math. Phys. 383, 2021–2067 (2021) | **CITE-OK** | Confirmed: CMP 383, 2021–2067 (2021), DOI 10.1007/s00220-020-03825-x, arXiv:2004.14115. Title/authors/venue/year all correct. Defines propagation number for operator systems, spectral truncations of S¹ as Toeplitz example — matches usage in Paper 32 §III. |
| `hekkelman2022` | E. Hekkelman, "Truncated geometries from spectral triples," PhD thesis, Radboud (2022) | Not separately verified | Hekkelman thesis exists at Radboud per other search hits; not load-bearing enough to chase. |
| `hekkelman_mcdonald2024` | E.-M. Hekkelman, E.A. McDonald, "A noncommutative integral on spectrally truncated spectral triples, and a link with quantum ergodicity," J. Funct. Anal. 289 (2025), accepted; arXiv:2412.00628 | **CITE-OK** | Confirmed: arXiv:2412.00628, submitted Dec 2024, accepted to J. Funct. Anal. with DOI 10.1016/j.jfa.2025.111154. Title and authors match. |
| `connes1995` | A. Connes, "Noncommutative geometry and reality," J. Math. Phys. 36, 6194–6231 (1995) | **CITE-OK** | Confirmed via NASA ADS and Connes's own publication list: JMP 36, 6194 (1995), DOI 10.1063/1.531241. This is the KO-dim sign-table source for the audit at §IV. |
| `camporesi_higuchi1996` | R. Camporesi, A. Higuchi, "On the eigenfunctions of the Dirac operator on spheres and real hyperbolic spaces," J. Geom. Phys. 20, 1–18 (1996) | **CITE-OK** | Confirmed: J. Geom. Phys. 20, 1–18 (Sept 1996), DOI 10.1016/0393-0440(95)00042-9, arXiv:gr-qc/9505009. (Paper 32 omits the arXiv ID but the published details are right.) |
| `chamseddine_connes2010` | Chamseddine, Connes, three papers (spectral action principle CMP 1997; "Why the Standard Model" JGP 2008; spectral action + superconnection PRD 2011) | Plausible (not chased) | Multi-paper bibitem; primary work (Chamseddine–Connes spectral action) is well-known CMP 186, 731–750 (1997). Worth a closer pass in a follow-on review if any of the three is load-bearing for a specific claim. |
| `krajewski1998` | T. Krajewski, "Classification of finite spectral triples," J. Geom. Phys. 28, 1–30 (1998) | **CITE-OK** | Confirmed: J. Geom. Phys. 28(1–2), 1–30 (1998), arXiv:hep-th/9701081. |
| `paschke_sitarz2000` | M. Paschke, A. Sitarz, "Discrete spectral triples and their symmetries," J. Math. Phys. 39, 6191–6205 (1998) | **CITE-WRONG-METADATA** | The journal/volume/pages/title are correct (JMP 39, 6191–6205, 1998, arXiv:q-alg/9612029), **but the bibitem KEY is `paschke_sitarz2000`** while the actual paper is 1998. Cosmetic but inconsistent. |
| `cacic2009` | B. Ćaćić, "Real structures on almost-commutative spectral triples," Lett. Math. Phys. 100, 181–202 (2012) | **CITE-WRONG-METADATA** | Title and author confirmed. **But the published version is Lett. Math. Phys. 103, 793–816 (2013)**, DOI 10.1007/s11005-013-0616-7, arXiv:1209.4832. Paper 32's volume (100), pages (181–202), and year (2012) are all wrong. The bibitem key `cacic2009` also doesn't match either the arXiv year (2012) or the published year (2013). |
| `latremoliere2026` | F. Latrémolière, "Spectral continuity of almost commutative manifolds for the C¹ topology on Riemannian metrics," arXiv:2603.19128 (2026) | **CITE-OK** | Confirmed: arXiv:2603.19128, submitted March 19, 2026, by Frédéric Latrémolière, title verbatim. |
| `bizi_brouder_besnard2018` | N. Bizi, C. Brouder, F. Besnard, "Towards a noncommutative geometry of the Standard Model with neutrino mixing," J. Math. Phys. 59, 062303 (2018), arXiv:1611.07062 | **CITE-MISATTRIBUTED** | Authors, journal/volume/year, and arXiv ID are correct **but the title is wrong**. arXiv:1611.07062 = JMP 59, 062303 (2018) is *"Space and time dimensions of algebras with applications to Lorentzian noncommutative geometry and quantum electrodynamics."* The cited title "Towards a noncommutative geometry of the Standard Model with neutrino mixing" is the title of a *different* Bizi–Brouder–Besnard adjacent work (or an earlier preprint title; not the JMP-published title). |
| `vandungen2016` | K. van den Dungen, "Krein spectral triples and the fermionic action," Math. Phys. Anal. Geom. 19, 4 (2016), arXiv:1505.01939 | **CITE-OK** | Confirmed: arXiv:1505.01939, MPAG 19, 4 (2016). Used in §VIII Lorentzian sections for the Krein–spectral-triple substrate. |
| `nieuviarts2025a` | A. Nieuviarts, "Emergence of Lorentz symmetry from an almost-commutative twisted spectral triple," arXiv:2502.18105 (2025) | **CITE-OK** | Confirmed: arXiv:2502.18105, submitted Feb 25, 2025, by Gaston Nieuviarts. (Paper 32 abbreviates as "A. Nieuviarts" — initial-only abbreviation is conventional, but the given name is Gaston.) |
| `witten2018_rmp` | E. Witten, "APS Medal... entanglement properties of quantum field theory," Rev. Mod. Phys. 90, 045003 (2018), arXiv:1803.04993 | **CITE-OK** | Title in Paper 32 is the APS Medal preamble; published RMP title is *"Notes on Some Entanglement Properties of Quantum Field Theory"* — same paper, different framing. RMP vol/page/year and arXiv ID all confirmed. |

## Problems found

### CITE-DOESNT-SUPPORT — `perez_sanchez2024` paired with the "Yang–Mills without Higgs" correction
**Location:** five places in the paper (lines 64, 242, 2369, 2473, 4575).
**Exact phrase affected (line 64 / abstract):** *"...the Perez-Sanchez correction~\cite{perez_sanchez2024,perez_sanchez2025} that the continuum limit is Yang–Mills without Higgs..."*
**What I found:** arXiv:2401.03705 (Perez-Sanchez Jan 2024) is "Bratteli networks and the Spectral Action on quivers." Its abstract derives "Yang–Mills–Higgs theory on flat space as a smooth limit," not Yang–Mills without Higgs. It explicitly says it "paved the way for some of our results in independent directions" — it is not a correction of Marcolli-vS. The "without Higgs" content lives entirely in arXiv:2508.17338 (perez_sanchez2025), which is explicitly titled "Comment on 'Gauge networks in noncommutative geometry'" and whose abstract states *"the continuum limit of this theory is the Yang–Mills action functional, without a Higgs scalar."*
**Fix:** Remove `perez_sanchez2024` from the five `[perez_sanchez2024,perez_sanchez2025]` bundles wherever the load-bearing claim is "Yang–Mills without Higgs." `perez_sanchez2024` can still be cited as a related/independent quiver-spectral-action work (it is, per its own abstract); but it is not the "correction." This is the single highest-impact finding.

### CITE-MISATTRIBUTED — `bizi_brouder_besnard2018`
**Location:** bibliography line 4904–4908.
**What I found:** Authors / venue / year / arXiv ID are correct (Bizi, Brouder, Besnard / JMP 59, 062303 / 2018 / arXiv:1611.07062), but the cited title *"Towards a noncommutative geometry of the Standard Model with neutrino mixing"* is wrong — the correct title is *"Space and time dimensions of algebras with applications to Lorentzian noncommutative geometry and quantum electrodynamics."*
**Fix:** Replace the title string. (This is the Fursaev–Solodukhin failure mode in miniature — right ID, wrong title — caught here at the bibliography level rather than after a memo.)

### CITE-WRONG-METADATA — `cacic2009`
**Location:** bibliography line 4794–4797.
**What I found:** Title and author correct. Volume, pages, year all wrong: cited as "Lett. Math. Phys. 100, 181–202 (2012)"; actual = Lett. Math. Phys. **103**, **793–816** (**2013**), arXiv:1209.4832.
**Fix:** Replace the three numbers and reconsider the bibitem key (`cacic2009` does not correspond to either the arXiv year 2012 or the publication year 2013).

### CITE-WRONG-METADATA — `paschke_sitarz2000` (cosmetic)
**Location:** bibliography line 4789–4792.
**What I found:** Body of the entry is correct (JMP 39, 6191–6205, 1998). The bibitem KEY contains "2000" but the paper is from 1998.
**Fix:** Either rename the key to `paschke_sitarz1998` or leave the cosmetic mismatch and accept it as legacy.

## Priority / novelty claims

| claim (verbatim) | location | searched | prior art found? | recommendation |
|:-----------------|:---------|:---------|:-----------------|:---------------|
| "Every structural ingredient of the GeoVac construction has published precedent; what is not duplicated is the α prediction at 8.8×10⁻⁸ (Paper 2) from a discrete spectral construction." | lines 67–71 (abstract), 2367–2371, 3682–3684, 4603–4606 | "fine structure constant discrete spectral action noncommutative geometry 1/137"; specific Connes/Chamseddine/Marcolli-vS searches | Web does not surface a peer-reviewed NCG/discrete-spectral derivation of α at 8.8×10⁻⁸ precision. Several non-peer-reviewed Medium / ResearchGate proposals exist (Pajuhaan, Pol Alan, etc.) but none are in the math.OA/NCG literature. | The claim is **appropriately phrased**: "not duplicated" is the right hedge (not "first in the literature"). I cannot upgrade this to "definitely first" (the honest ceiling on priority — Sec. 2 of CITATION_CHECKER spec) but I can say that nothing I found falsifies it. Leave as-is. |

## Bottom line

**Bibliography verdict: YELLOW.**

What is right: the load-bearing lineage placement is real. arXiv:1301.3480 (Marcolli–vS 2014) genuinely is a graph-spectral-triple framework whose spectral action yields Wilson lattice gauge theory. arXiv:2508.17338 (Perez-Sanchez 2025) genuinely demonstrates the Yang–Mills-without-Higgs correction. Connes–vS 2021 (CMP 383) is correctly cited. Camporesi–Higuchi 1996 (JGP 20, 1–18) is correctly cited. Connes 1995 (JMP 36, 6194) is correctly cited. The α priority is hedged correctly.

What would embarrass us:

1. **`perez_sanchez2024` bundled with `perez_sanchez2025` as "the Perez-Sanchez correction" in five places** — perez_sanchez2024 derives Yang–Mills–Higgs, not Yang–Mills without Higgs. An NCG referee will spot this immediately. This is the must-fix.
2. **`bizi_brouder_besnard2018` has the title of a different paper** — the Fursaev failure mode. Must-fix.
3. **`cacic2009` has wrong volume (100 → 103), wrong pages (181–202 → 793–816), wrong year (2012 → 2013)** — checkable in one click via the LMP web page. Must-fix.
4. `paschke_sitarz2000` key vs 1998 paper — cosmetic, can stay.

Top finding (one sentence): **`perez_sanchez2024` (arXiv:2401.03705) does NOT correct Marcolli–vS to "Yang–Mills without Higgs" — that result lives in `perez_sanchez2025` (arXiv:2508.17338) alone, and bundling them in five places overstates one work and misattributes the load-bearing correction to two papers when only one supports it.**

## Sources verified

- [arXiv:1301.3480 — Marcolli & van Suijlekom, "Gauge networks in noncommutative geometry"](https://arxiv.org/abs/1301.3480)
- [arXiv:2401.03705 — Perez-Sanchez, "Bratteli networks and the Spectral Action on quivers" (2024)](https://arxiv.org/abs/2401.03705)
- [arXiv:2508.17338 — Perez-Sanchez, "Comment on 'Gauge networks in noncommutative geometry'" (2025)](https://arxiv.org/abs/2508.17338)
- [arXiv:2004.14115 — Connes & van Suijlekom, "Spectral truncations in noncommutative geometry and operator systems" (CMP 383, 2021–2067)](https://arxiv.org/abs/2004.14115)
- [arXiv:2412.00628 — Hekkelman & McDonald, J. Funct. Anal. 2025](https://arxiv.org/abs/2412.00628)
- [arXiv:gr-qc/9505009 — Camporesi & Higuchi, J. Geom. Phys. 20, 1–18 (1996)](https://arxiv.org/abs/gr-qc/9505009)
- [arXiv:hep-th/9701081 — Krajewski, "Classification of finite spectral triples"](https://arxiv.org/abs/hep-th/9701081)
- [arXiv:q-alg/9612029 — Paschke & Sitarz, "Discrete spectral triples and their symmetries"](https://arxiv.org/abs/q-alg/9612029)
- [arXiv:1209.4832 — Ćaćić, "Real structures on almost-commutative spectral triples" (LMP 103, 793–816, 2013)](https://arxiv.org/abs/1209.4832)
- [arXiv:1505.01939 — van den Dungen, "Krein spectral triples and the fermionic action"](https://arxiv.org/abs/1505.01939)
- [arXiv:1611.07062 — Bizi, Brouder, Besnard, "Space and time dimensions of algebras..." (JMP 59, 062303, 2018)](https://arxiv.org/abs/1611.07062)
- [arXiv:2502.18105 — Nieuviarts, "Emergence of Lorentz symmetry..."](https://arxiv.org/abs/2502.18105)
- [arXiv:1803.04993 — Witten, "Notes on Some Entanglement Properties of Quantum Field Theory" (RMP 90, 045003, 2018)](https://arxiv.org/abs/1803.04993)
- [arXiv:2603.19128 — Latrémolière, "Spectral continuity of almost commutative manifolds..."](https://arxiv.org/abs/2603.19128)
- [Connes "Noncommutative geometry and reality" J. Math. Phys. 36, 6194 (1995)](https://alainconnes.org/wp-content/uploads/reality.pdf)
