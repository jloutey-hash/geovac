# Citation & Novelty Check: Paper 16 — The Periodic Table Recast as $S_N$ Representation Theory on $S^{3N-1}$

**Date:** 2026-06-01
**Wave:** 1 (S_N rep-theory load-bearing focus)
**Checker:** CITATION_CHECKER agent

## Citation table

| \cite key | claimed as | verdict | what I found (URL / arXiv / DOI) |
|:----------|:-----------|:--------|:----------------------------------|
| `Mendeleev1869` | Periodic table established empirically in 1869 | **CITE-OK** | Mendeleev, "Sootnoshenie svoistv s atomnym vesom elementov", *Zh. Russ. Khim. Obshch.* vol. 1, 1869. Confirmed via History-of-Information records and Bookvica. The cited journal name and year are correct (the paper begins on p. 60; the periodic table appears on p. 70). |
| `Fock1935` | Hyperspherical separation of N-electron Schrödinger / SO(4) on hydrogen | **CITE-OK** | Fock, "Zur Theorie des Wasserstoffatoms," *Z. Phys.* **98**, 145 (1935). Confirmed via NASA-ADS, Springer Link, viaLibri. DOI 10.1007/BF01336904. Note: this paper is the SO(4) hydrogen paper, NOT the N-electron hyperspherical paper — it is correctly used in the paper as one of two parts of the cited pair. |
| `Fock1954` | Hyperspherical separation framework for N-body systems | **CITE-CANT-INDEPENDENTLY-VERIFY (metadata plausible)** | The 1954 Fock paper on N-body hyperspherical harmonics in *Izvestiya Akademii Nauk SSSR, Ser. Fiz.* vol. 18 is a well-known and frequently-cited reference in the few-body literature. Volume number 18 (1954) and page 161 are commonly cited in the field. zbMATH confirms the journal series exists. Could not fetch the article itself (Soviet-era journal, not indexed online). No reason to doubt; flagging as could-not-independently-verify (not a problem). |
| `Avery2000` | Hyperspherical formalism; specifically Eq. (18) $\nu_{\min} = N - \lambda_1$ | **CITE-OK (existence) / CITE-CANT-INDEPENDENTLY-VERIFY (specific formula)** | Book confirmed: J. Avery, *Hyperspherical Harmonics and Generalized Sturmians*, Kluwer (Progress in Theoretical Chemistry and Physics, vol. 4), 2000, ISBN 9780792360872. Confirmed via Springer, AbeBooks, ResearchGate. The book IS the standard reference for symmetry-adapted hyperspherical harmonics including antisymmetric (fermionic) symmetrization. Could not independently confirm from public sources that the specific formula $\nu_{\min} = N - \lambda_1$ appears verbatim there; this requires book access. Plausible. |
| `Weyl1946` | Schur–Weyl duality, textbook material | **CITE-OK (mild metadata note)** | H. Weyl, *The Classical Groups: Their Invariants and Representations*, Princeton (2nd ed., 1946). First edition was 1939; the cited 1946 edition is the standard revised edition. Confirmed via Wikipedia, Princeton Press, Internet Archive. |
| `FultonHarris1991` | Schur–Weyl duality (paired with Weyl1946) | **CITE-OK** | W. Fulton & J. Harris, *Representation Theory: A First Course*, Graduate Texts in Mathematics 129, Springer, 1991, ISBN 0-387-97495-4. Confirmed via Cambridge Mathematical Gazette review and Springer Link. |
| `loutey_paper0` | GeoVac internal | **N/A** (internal cross-reference; not web-checkable) |
| `loutey_paper1` | GeoVac internal | **N/A** |
| `loutey_paper7` | GeoVac internal | **N/A** |
| `loutey_paper13` | GeoVac internal | **N/A** |
| `KnowlesHandy1984` | Cited in bibliography but NO in-text `\cite{KnowlesHandy1984}` found in body | **CITE-OK (but orphan)** | P.J. Knowles & N.C. Handy, "A new determinant-based full configuration interaction method," *Chem. Phys. Lett.* **111**, 315 (1984). Confirmed via NASA-ADS, ScienceDirect, Semantic Scholar. Citation exists and metadata is correct. However, this bibitem appears to have no `\cite{}` in the body of the paper — it's an orphan reference. |
| `Dirac1928` | Dirac equation, prediction of $Z\alpha \to 1$ instability | **CITE-OK** | P.A.M. Dirac, "The quantum theory of the electron," *Proc. R. Soc. Lond. A* **117**, 610 (1928). DOI 10.1098/rspa.1928.0023. Confirmed via Royal Society, NASA-ADS, scirp.org. |
| `Greiner1985` | $Z \sim 137$ Dirac instability, $Z \sim 170$ pair-creation | **CITE-OK** | W. Greiner, B. Müller, J. Rafelski, *Quantum Electrodynamics of Strong Fields*, Springer, Berlin, 1985. Standard reference for the Dirac dive and supercritical QED. Confirmed via Springer, Amazon, Google Books. |
| `Pyykkoe2012` | Dirac–Fock calculations on superheavy elements; Z ≤ 172 periodic table extension | **CITE-WRONG-METADATA (bib key only)** | The actual paper is P. Pyykkö, *Phys. Chem. Chem. Phys.* **13**, 161 (2011), DOI 10.1039/C0CP01575J — confirmed published online 22 October 2010, journal year 2011. The **bibitem text** correctly says 2011. The **bib key** says `Pyykkoe2012` which is misleading (should be `Pyykkoe2011`). Cosmetic, not a citation error per se. |
| `Schwerdtfeger2015` | Relativistic effects in superheavy elements (book chapter) | **CITE-CANT-FIND / likely CITE-MISATTRIBUTED** | The cited chapter — Schwerdtfeger, L.F. Pashtedsky, A. Pernpointner, B. Fricke, "Relativistic effects in superheavy elements," in *The Chemistry of Superheavy Elements*, 2nd ed., Springer, Berlin, **2015** — could not be found. **Multiple problems**: (i) The 2nd edition of *The Chemistry of Superheavy Elements* is edited by **M. Schädel & D. Shaughnessy**, published **2014** (not 2015). (ii) The chapter on "Theoretical Chemistry of the Heaviest Elements" in that book is by **V. Pershina** (single author), not Schwerdtfeger et al. (iii) The author "L.F. Pashtedsky" is a misspelling of **L.F. Pašteka** (Schwerdtfeger's actual collaborator). (iv) A real 2015 Schwerdtfeger–Pašteka paper exists but it is in *Nuclear Physics A* **944**, 551 (2015), "Relativistic and quantum electrodynamic effects in superheavy elements" — a journal review article, not a book chapter, with author list **Schwerdtfeger, Pašteka, Punnett, Bowman** (not Pernpointner/Fricke). The cited entry appears to be an amalgamation. See "Problems found" below. |

## Problems found (CITE-MISATTRIBUTED / DOESNT-SUPPORT / CANT-FIND)

### Problem 1 — `Schwerdtfeger2015` is likely fabricated / amalgamated

**Citation as printed:**
> P. Schwerdtfeger, L.F. Pashtedsky, A. Pernpointner, and B. Fricke,
> "Relativistic effects in superheavy elements,"
> in *The Chemistry of Superheavy Elements*, 2nd ed.
> (Springer, Berlin, 2015).

**Used in body:** Section VII.C "Three limits on the periodic table" — supports the claim that "nuclear instability ($Z \sim 126$), Dirac instability point nucleus ($Z \sim 137$), Dirac instability finite nucleus ($Z \sim 170$)" set the upper limit on the periodic table.

**What is verifiably wrong:**
- The 2nd edition of *The Chemistry of Superheavy Elements* (Springer) is by **Schädel & Shaughnessy** (eds.), **2014**, not 2015 (DOI 10.1007/978-3-642-37466-1). Confirmed via Amazon, Springer Link, AbeBooks.
- The chapter on theoretical chemistry of superheavy elements in that book is **"Theoretical Chemistry of the Heaviest Elements" by V. Pershina** (chapter 3, DOI 10.1007/978-3-642-37466-1_3), not by Schwerdtfeger et al.
- "L.F. Pashtedsky" is not a known researcher; the actual name is **L.F. Pašteka** (a known Schwerdtfeger collaborator).
- The real 2015 Schwerdtfeger–Pašteka piece is a *Nuclear Physics A* article (vol. 944, p. 551), title "Relativistic and quantum electrodynamic effects in superheavy elements," authored by Schwerdtfeger, Pašteka, Punnett, Bowman — different venue, different title, different author list.

**Risk classification:** RED. This citation looks like an LLM/transcription confabulation: the right topic, the right Springer book series, the right Schwerdtfeger involvement, but the **specific publication described doesn't appear to exist**. This is the same failure mode as the Fursaev-Solodukhin / hep-th/9512134 incident documented in this project's dead-ends record (CLAUDE.md §3).

**Recommendation:** Replace with one of:
- V. Pershina, "Theoretical Chemistry of the Heaviest Elements," in Schädel & Shaughnessy (eds.), *The Chemistry of Superheavy Elements*, 2nd ed., Springer 2014, ch. 3, DOI 10.1007/978-3-642-37466-1_3.
- P. Schwerdtfeger, L.F. Pašteka, A. Punnett, P.O. Bowman, "Relativistic and quantum electrodynamic effects in superheavy elements," *Nucl. Phys. A* **944**, 551 (2015).

The Pyykkö 2011 reference already covers the periodic-table-extension content; the new reference is needed primarily for the $Z \sim 170$ pair-creation limit, which is well-established and could also be supported by Greiner-Müller-Rafelski 1985 (already cited).

### Problem 2 — `Pyykkoe2012` bib key mislabel (cosmetic)

The bib key suggests 2012; the cited journal year is 2011 (DOI 10.1039/C0CP01575J, online 22 Oct 2010, vol. 13 issue dated 2011). The bibitem text correctly shows 2011. Cosmetic only — anyone reading the printed bibliography sees the right year. Rename to `Pyykkoe2011` for clarity at next pass.

### Problem 3 — `KnowlesHandy1984` is an orphan citation

The bibitem appears in the bibliography but I find no `\cite{KnowlesHandy1984}` in the body of the .tex. This is likely vestigial from an earlier draft. Either reinsert the `\cite` somewhere it supports a claim (e.g., the FCI computational context for testing the present framework), or remove the bibitem.

## Priority / novelty claims

| claim (verbatim) | location | searched | prior art found? | recommendation |
|:-----------------|:---------|:---------|:-----------------|:---------------|
| "This specific identification appears to be new." (re: universal $\nu = N - 2$ for $S < N/2$) | §II.B item 1 | Searches for "S_N representation theory" + "hyperspherical" + "periodic table"; for "minimum K hyperspherical Pauli Young diagram"; for prior group-theoretic treatments of the periodic table | **No prior art found** that derives the universal $\nu = N-2$ result for spin-$S < N/2$ atoms with this $S_N$/Young-diagram language. The closest published work uses the **SO(4,2)×SU(2) / Rumer–Fet group** approach (Pyykkö "essay on periodic tables"; "Group Theoretical Description of the Periodic System" MDPI Symmetry 14:137; arXiv:2501.18272 "Periodic Table and SO(4,4)") — these are different mathematical structures, not the $S_N$ approach Paper 16 uses. | Paper's phrasing "appears to be new" is **appropriately hedged**. Per CITATION_CHECKER protocol, web searches cannot establish "first in literature" — only downgrade. The current wording already is at the safe ceiling ("appears to be new"); recommend KEEP. |
| "The four-type classification (A/B/C/D)... provides a group-theoretic reformulation of the empirical pattern" | §II.B item 2 | Same as above + searched for prior periodic table classifications using Young-diagram shape | No prior art found that uses the specific Young-diagram-shape classification to partition into noble/alkali/alkaline-earth types. There is rich prior art on group-theoretic *periodicity* (Rumer–Fet, SO(4,2)×SU(2)) but not on this specific reformulation. | Phrasing is appropriately framed as "reformulation" rather than discovery — KEEP. |
| "The per-pair angular momentum limit ($\mufree/\binom{N}{2} \to 4$)" | §II.B item 3 | Searched for "Pauli centrifugal" + "per pair" + "hyperspherical" | No prior art found. The result is an elementary consequence of the universal $\nu = N-2$ formula. | KEEP as stated. |
| "The Dirac instability as metric vs. topological singularity" | §II.B item 4 | Searched for prior work on this specific framing | The factual content (Dirac instability is a metric/relativistic effect on a fixed topology) is well-known; the **framing in terms of $S^{3N-1}$ topological smoothness** through $Z=137$ appears to be the new contribution. | The framing-level novelty is plausibly new; the *physical content* it organizes is standard. Phrasing is appropriately scoped. KEEP. |
| "the periodic law itself... has not been derived from a single mathematical principle" | §I, ¶2 | Searched for prior group-theoretic derivations of periodic law | Several prior frameworks (Rumer–Fet, SO(4,2)×SU(2)) claim "group-theoretic derivation of periodic law" — the paper should acknowledge these as related work. | YELLOW — the §I claim glosses over a substantial prior literature. Recommend adding a sentence to §II "Known results" or §VII.D acknowledging the SO(4,2)×SU(2) / Rumer–Fet lineage and stating how the $S_N$ approach is complementary (focuses on permutation symmetry; the Rumer–Fet group is a dynamical symmetry of a different kind). |

## Bottom line

**Verdict: YELLOW** — bibliography is mostly clean, but one likely-fabricated citation (`Schwerdtfeger2015`) is the kind of fault that has embarrassed the project before (Fursaev-Solodukhin / hep-th/9512134 dead-end). Two cosmetic issues (`Pyykkoe2012` bib key, `KnowlesHandy1984` orphan).

**Citations checked:** 15 (13 external + 4 GeoVac-internal not web-checkable, of which 4 are skipped).

**Problems by verdict:**
- CITE-MISATTRIBUTED / CITE-CANT-FIND: 1 (`Schwerdtfeger2015`) — load-bearing for §VII.C "Three limits" claim.
- CITE-WRONG-METADATA: 1 cosmetic (`Pyykkoe2012` bib key).
- Orphan citation: 1 (`KnowlesHandy1984` no in-text use).
- CITE-OK: 9.
- CITE-CANT-INDEPENDENTLY-VERIFY (not problems, just flagged): 2 (Fock1954, Avery2000 specific formula).

**Top finding:** `Schwerdtfeger2015` is the load-bearing risk. The cited four-author chapter with title "Relativistic effects in superheavy elements" in *The Chemistry of Superheavy Elements* 2nd ed. (Springer 2015) does not appear to exist. The Springer 2nd edition is dated 2014 (Schädel/Shaughnessy eds.), the relevant chapter is by Pershina alone, and the real 2015 Schwerdtfeger paper is in *Nucl. Phys. A* with a different title and a different author list ("Pashtedsky" is a misspelling of "Pašteka"; "Pernpointner" and "Fricke" are not on the 2015 NPA paper). Recommend replacement with either Pershina 2014 (book chapter) or Schwerdtfeger et al. 2015 *NPA* (journal review) depending on which content is intended to be supported.

**S_N representation theory (load-bearing focus, per Wave 1 instruction):**
The two textbook citations for $S_N$ rep theory (Weyl 1946, Fulton-Harris 1991) are correct and standard. The load-bearing formula $\nu_{\min} = N - \lambda_1$ (Eq. 18, attributed to Avery 2000) is plausible but could not be independently verified from public sources — flagged. The novelty claim that the **universal $\nu = N-2$ identification + four-type A/B/C/D classification** is new in the $S_N$ framework returned no prior art in my searches; current wording ("appears to be new") is appropriately at the honest ceiling. **Recommend adding** a one-sentence acknowledgment of the SO(4,2)×SU(2) / Rumer–Fet group-theoretic-periodicity literature in §II or §VII.D so the reader sees the present work as complementary, not as unaware of a substantial prior literature.

## Sources

- [Fock 1935, NASA-ADS](https://ui.adsabs.harvard.edu/abs/1935ZPhy...98..145F)
- [Avery 2000, Hyperspherical Harmonics and Generalized Sturmians (AbeBooks listing)](https://www.abebooks.com/9780792360872/Hyperspherical-Harmonics-Generalized-Sturmians-Progress-0792360877/plp)
- [Weyl 1946, Internet Archive](https://archive.org/details/classicalgroupst0000weyl)
- [Fulton-Harris 1991, Springer Link](https://link.springer.com/book/10.1007/978-1-4612-0979-9)
- [Knowles-Handy 1984, NASA-ADS](https://ui.adsabs.harvard.edu/abs/1984CPL...111..315K/abstract)
- [Dirac 1928, Royal Society](https://royalsocietypublishing.org/doi/10.1098/rspa.1928.0023)
- [Greiner-Müller-Rafelski 1985, Springer Link](https://link.springer.com/book/10.1007/978-3-642-82272-8)
- [Pyykkö 2011, RSC](https://pubs.rsc.org/en/content/articlehtml/2011/cp/c0cp01575j)
- [Schädel & Shaughnessy (eds.) 2014, *The Chemistry of Superheavy Elements* 2nd ed., Springer Link](https://link.springer.com/book/10.1007/978-3-642-37466-1)
- [Pershina 2014 chapter, Springer Link](https://link.springer.com/chapter/10.1007/978-3-642-37466-1_3)
- [Schwerdtfeger et al. 2015 NPA, NASA-ADS](https://ui.adsabs.harvard.edu/abs/2015NuPhA.944..551S)
- [Mendeleev 1869 periodic table, History of Information](https://www.historyofinformation.com/detail.php?id=2898)
- [SO(4,2)×SU(2) group-theoretic periodic system, MDPI Symmetry 14:137](https://www.mdpi.com/2073-8994/14/1/137)
- [Periodic Table and SO(4,4), arXiv:2501.18272](https://arxiv.org/pdf/2501.18272)
