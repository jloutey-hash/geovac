# Confidence Review: Paper 31 — The Universal/Coulomb Partition of the GeoVac Spectral Triple

## Pass A — Content audit

### Claim inventory + verdicts

| # | Claim | Location | Verdict | Rests on | Evidence I produced |
|---|---|---|---|---|---|
| 1 | The GeoVac framework decomposes cleanly into a universal sector (depending only on the algebra $\mathcal{A}$ and angular structure) and a potential-specific sector (depending on $D$) | Abstract, §II "The decomposition," §VI Conclusion | A (organizational reading, well-grounded) | EXTERNAL (Connes spectral-triple split + GeoVac papers) | The split is the standard $\mathcal{A}$/$D$ separation in spectral-triple theory; the paper itself states ("the decomposition is not novel as a piece of spectral triple formalism", §II ¶7) that the framing is GeoVac-organizational. Verdict A as a piece of organization. |
| 2 | Angular sparsity ERI density: 7.81% at $l_{\max}=1$, 2.76% at $l_{\max}=2$, 1.44% at $l_{\max}=3$, 0.90% at $l_{\max}=4$, 0.62% at $l_{\max}=5$, bit-identical across Coulomb, HO, Woods–Saxon, square well, and Yukawa at $l_{\max}=2$ | §III Eq. eri_density | B (internal cross-check) | GEOVAC-ONLY (depends entirely on Paper 22's numerics) | Cross-checked vs `papers/group3_foundations/paper_22_angular_sparsity.tex` lines 23–25, 255–258, 343–347: all five numbers match Paper 22 verbatim, including the five-potential bit-identical pattern. No external verification possible without re-implementing the Gaunt enumeration. |
| 3 | Fock projection rigidity theorem (Paper 23): $-Z/r$ is unique central potential with $\ell$-independent spectrum within each $n$ | §IV.A Theorem 1 | B (well-established for SO(4)/Coulomb in textbooks; specific GeoVac formulation rests on Paper 23) | MIXED (SO(4) symmetry of Coulomb is textbook; the "rigidity-as-uniqueness-of-Fock-projection" formulation is Paper 23 internal) | Bander–Itzykson 1966 verified online; SO(4) → $\ell$-degeneracy is a textbook result. The wrapping as a "rigidity theorem" of the Fock projection is GeoVac. |
| 4 | Bargmann–Segal rigidity theorem (Paper 24): 3D isotropic HO is unique central potential whose spectrum arises from Euler operator on $H^2(S^5)$ restricted to symmetric (N,0) irreps of SU(3) | §IV.B Theorem 2 | B (internal; rests on Paper 24) | GEOVAC-ONLY (the rigidity claim is Paper 24's contribution) | Jauch–Hill 1940 (SU(3) symmetry of 3D HO), Borel–Weil, and Bargmann–Segal–Hall are textbook inputs as the paper acknowledges. The composite rigidity wrapping is GeoVac-internal. |
| 5 | "Coulomb is honestly $\pi$-laden at every level above the bare graph; the HO is bit-exactly $\pi$-free" | §V; Table tab:pi_distribution | B (rests on Papers 7/18/24/28) | GEOVAC-ONLY | The $\pi$-free certificate for the Bargmann–Segal lattice is Paper 24; the transcendental taxonomy is Paper 18. Internal-consistency check, not independent verification. |
| 6 | The Coulomb/HO asymmetry has **three layers, not two** (spectrum-computing role of $L_0$, calibration $\pi$, Wilson physical content of non-abelian gauge) — §V title, Table tab:three_layers, §V.C "Why three and not two" | §V throughout | **C (overstated/stale: the corpus now has four layers)** | GEOVAC-ONLY | Paper 24 §V "asymmetry_layer4" (added Sprint L2-F.1, 2026-05-17, per CLAUDE.md §6 Paper 24 row); Paper 25 lines 956–963 also describe a four-layer reading (with Sprint ST-SU3 closure). Paper 31's three-layer framing is now superseded by the corpus's four-layer framing. |
| 7 | Sprint 5 result on $S^5$ gauge structure: U(1) Wilson–Hodge transfers verbatim; non-abelian SU(3) extension "is not natural" (transitions between $(N,0)$ and $(N+1,0)$ are Clebsch–Gordan intertwiners, not group elements) | §V.A Layer 3 | B (rests on Paper 25 Sprint 5 / Sprint ST-SU3) | GEOVAC-ONLY | Paper 25 confirms the qualitative finding; the published-elsewhere SU(3)-on-$S^5$ NO-GO is Sprint ST-SU3 / Paper 25 internal. |
| 8 | "Composed-block $O(Q^{2.5})$ Pauli scaling … 51× to 1712× advantage over Gaussian baselines across the GeoVac molecule library" | §VI.A bullet 2 | **C (overstated)** | GEOVAC-ONLY | Paper 14 says "two or more orders of magnitude" (line 46) and presents specific numbers (e.g., LiH 334 Pauli vs STO-3G 276 — line 911). The "51× to 1712×" exact ratios come from CLAUDE.md §1.5 (LiH/BeH2/H2O), not from Paper 14's body. Paper 31 attributes the precise figures to \cite{paper14}, which doesn't carry them. Replace with "two or more orders of magnitude" matching Paper 14's own framing, or cite Paper 20 / CHANGELOG. |
| 9 | Parametric/central sector: $D \mapsto D + c \cdot \mathbb{1}$ leaves every Connes axiom unchanged (bounded commutators, order-one, real structure $J$, $\gamma$-grading, KO-dimension) | §VI.C | A on commutator/order-one; D on the $J$/$\gamma$ claims (signature-dependent) | EXTERNAL (linear algebra) | Bounded commutator $[D+c\mathbb{1}, a] = [D, a]$ is trivial. The $\{J, D\} = 0$ vs $JD = +DJ$ structure depends on KO-dimension: in GeoVac's KO-dim 3 ($JD = +DJ$), an additive real shift preserves $J(D+c) = (D+c)J$ — fine. The paper's blanket "every Connes axiom unchanged" is correct in the GeoVac convention but would be sign-dependent for other KO-dimensions. Worth a half-sentence clarification rather than a defect call. |
| 10 | Sprint L0 signature partition: 17 transfers-freely / 4 Wick-map-free / 5 Euclidean-specific / 2 mixed of Paper 34's 28 projections | §VIII.A | B | GEOVAC-ONLY | Per `debug/lorentzian_l0_audit_memo.md` (cited) and Paper 34 §V.E. Internal classification — no external check possible. |
| 11 | "BBB sign table at $(m,n)=(3,1)$ … M3 sub-mechanism … trivializes on chirality-symmetric spectrum" | §VIII.D | D (predicted/unverified; the paper labels it a Sprint L2 prediction) | GEOVAC-ONLY | The paper is appropriately hedged ("Sprint L2 prediction worth verifying"). No claim to a proof. Verdict D = unverifiable here. |
| 12 | "Empirical verification of the partition at 28/28 projections" (signature partition) | §VIII.F | B | GEOVAC-ONLY | Cross-corpus confirmation via Paper 32 §VIII.E and Paper 42 §§5–8 (per CLAUDE.md §6). Internal-only. |
| 13 | Two-body verification: gauged tensor-product spectral action gives 100% $m$-conserving, 100% Gaunt-compatible, pure $k=0$ monopole at both $n_{\max}=2,3$; Pearson correlation 0.41–0.58 decreasing with $n_{\max}$ for the spectral-action two-body; Dirac resolvent 0.81/0.75 also decreasing | §IX.A,B | B | GEOVAC-ONLY | Per CLAUDE.md §2 sprint entries 2026-06-01 (resolvent two-body CLEAN NEGATIVE, sprint nuclear + tensor-product spectral action). Numbers consistent with the listed memos. Tested only inside the framework. |
| 14 | "Conversion between Gegenbauer triple integral on $S^3$ and Slater integral in flat coordinates is the Fock projection … the conformal factor $r_{\rm chordal} = 2R\sin(\chi_{\rm geodesic}/2R)$" | §IX.C | A (geometric) | EXTERNAL | Chordal–geodesic identity is the standard Fock chordal-distance result; correctly stated. |

### Numbers I recomputed

| Claim | Paper's figure | Independent reference | Recomputed value/error | Survives? |
|---|---|---|---|---|
| At $l_{\max}=3$ "the new system's two-body Hamiltonian is at least 98.56% sparse" (§VI.A) | 98.56% | Paper 31 §III: 1.44% nonzero | $100 - 1.44 = 98.56$ ✓ | Yes |
| At $l_{\max}=2$ 97.24% zero (implicitly in tab:eri_density) | 2.76% nonzero | Paper 22 line 343 ("97.24 & 0.00 & 2.76") | $100 - 2.76 = 97.24$ ✓ | Yes |
| Paper 22 ERI density numbers carry to Paper 31 §III | 7.81 / 2.76 / 1.44 / 0.90 / 0.62 | Paper 22 §III tables | All match bit-exact | Yes |
| "Five potentials Coulomb / HO / Woods–Saxon / square well / Yukawa at $l_{\max}=2$" | bit-identical zero patterns | Paper 22 lines 343–347 | Match | Yes |

### Circularity map

**Bottoms out EXTERNAL:**
- Connes spectral-triple formalism (Connes 1994 book, verified standard).
- Marcolli–vS gauge networks (arXiv:1301.3480, verified existence + content).
- Bander–Itzykson 1966 SO(4) of Coulomb (verified textbook).
- Lichnerowicz 1963 (textbook, verified).
- Jauch–Hill 1940 SU(3) of HO (textbook).
- Bargmann 1961 / Hall 1994 (textbook).
- Hartle–Hawking 1976 / Sewell 1982 / Bisognano–Wichmann 1976 / Unruh 1976 (all standard).
- BBB 2018 (arXiv:1611.07062, verified existence; SEE PASS B citation issues).
- Chordal-distance identity (§IX.C, textbook geometry).

**Bottoms out GEOVAC-ONLY (the chains most at risk):**
- The five-numbers ERI density and five-potential verification (Claim 2) rest entirely on Paper 22's Python implementation.
- The two rigidity theorems (Claims 3, 4) wrap textbook symmetry results in GeoVac-specific "rigidity as uniqueness of projection" framing; the wrapping is Paper 23/24's contribution.
- The three-layer (now four-layer) Coulomb/HO sharpening (Claim 6) rests on Paper 25 + CLAUDE.md §2 sprint records.
- The "51× to 1712×" Pauli ratio (Claim 8) traces to CLAUDE.md §1.5, not to Paper 14's body.
- The 100%/Gaunt/k=0 two-body diagnostics (Claim 13) rest on Sprint 2026-06-01 memos, none yet cross-checked outside GeoVac.

**MIXED:**
- Claim 3 (Fock rigidity): textbook SO(4) input + GeoVac-internal wrapping.
- Claim 4 (Bargmann–Segal rigidity): textbook SU(3) + Borel–Weil + Bargmann inputs + GeoVac-internal wrapping.
- Claim 9 (parametric central sector): trivial linear algebra + GeoVac convention.

### Overstatement findings

- **§V title "The Three-Layer Sharpening (Sprint 5)"**, and the entire §V framing, are STALE relative to the corpus. Per CLAUDE.md §6 Paper 24 row and Paper 25 lines 956–963, the asymmetry is now FOUR layers (the Sprint L2-F.1 modular-Hamiltonian addition added a fourth layer). Recommended replacement: §V title → "Sharpening: three layers in the Sprint 5 reading, four layers in the post-L2-F.1 reading"; §V.C "Why three and not two" → "Why three (post-Sprint 5) and four (post-Sprint L2-F.1)"; add one paragraph or footnote citing Paper 24 §V `subsec:asymmetry_layer4` and Paper 25's four-layer framing. Severity MEDIUM (load-bearing per the brief).
- **§VI.A bullet 2:** "51× to 1712× across the GeoVac molecule library" → "two or more orders of magnitude across the GeoVac molecule library, per Paper 14 §VI and Paper 20 §IV" (matching Paper 14's own quantitative wording). The exact ratios should be sourced from CLAUDE.md §1.5 or moved to a footnote with the specific molecules and references. Severity MEDIUM.
- **§I "almost verbatim"** for the lineage statement: appropriate caution would be "structurally a specialisation"; "verbatim" overstates the match. Paper 32 uses "structurally a specialization" (line 4571). Severity LOW.
- **Internal cross-reference defects**: lines 876–877 reference `sec:universal_sector` and `sec:potential_specific_sector`, but the actual labels are `sec:universal` and `sec:specific`. Will render as `??` in PDF. Severity LOW.

## Pass B — Citation and novelty

### Citation table

| \cite key | Claimed as | Verdict | What I found |
|---|---|---|---|
| `connes_book1994` | Spectral-triple data triple definition (§I, §II, §VI) | CITE-OK | Connes, *Noncommutative Geometry*, Academic Press 1994 — standard reference, correctly cited. |
| `chamseddine_connes1997` | Spectral action principle, Commun. Math. Phys. 186 731 (1997) | CITE-OK | Paper exists; standard reference; correctly cited. |
| `marcolli_vs2014` | Gauge networks, J. Geom. Phys. 75, 71 (2014); arXiv:1301.3480 | CITE-OK | Verified via arxiv.org/abs/1301.3480. Authors, journal, year, arXiv ID all correct. |
| `perez_sanchez2024` | "Yang–Mills without Higgs" continuum-limit correction, arXiv:2401.03705 (with "see also arXiv:2508.17338"); bibitem title "On the gauge networks of Marcolli–van Suijlekom: Yang–Mills without Higgs" | **CITE-MISATTRIBUTED + CITE-WRONG-METADATA + CITE-DOESNT-SUPPORT** | arXiv:2401.03705 actually titled "Bratteli networks and the Spectral Action on quivers" (verified). Its abstract derives Yang–Mills–**HIGGS** ("We show that a hermitian ('Higgs') matrix field emerges from the self-loops … derive the Yang–Mills–Higgs theory"). The "without Higgs" correction is solely arXiv:2508.17338 — a separate paper, "Comment on 'Gauge networks in noncommutative geometry'". Paper 31 bundles both arXiv IDs into one bibitem with the wrong title, mis-attributing the "without Higgs" content to the 2024 paper. **This is exactly the same misbundling the Paper 32 audit fixed yesterday**: Paper 32 has separate `perez_sanchez2024` and `perez_sanchez2025` bibitems with correct titles ("Bratteli networks…" and "Comment on…") and uses only `perez_sanchez2025` for the "without Higgs" correction. **Severity HIGH.** Fix: split into two bibitems matching Paper 32 (`perez_sanchez2024` = arXiv:2401.03705 = "Bratteli networks…"; `perez_sanchez2025` = arXiv:2508.17338 = "Comment on 'Gauge networks…'"). In §II line 200, change "Perez-Sanchez correction~\cite{perez_sanchez2024}" to "~\cite{perez_sanchez2025}". §I line 163 reads "Connes, Marcolli, and van Suijlekom~\cite{...,perez_sanchez2024}" — should add `perez_sanchez2025` for the without-Higgs angle, or be clear that 2024 is the more recent Bratteli-networks paper rather than the correction. |
| `lichnerowicz1963` | Lichnerowicz formula | CITE-OK | C. R. Acad. Sci. Paris vol. 257 (1963), standard reference. |
| `bander_itzykson1966` | Group theory and the hydrogen atom, Rev. Mod. Phys. 38, 330 and 346 (1966) | CITE-OK | Verified online: two-part paper at the cited pages, correct. |
| `fock1935` | Fock 1935 stereographic | CITE-OK | Standard Z. Phys. 98 145 (1935) reference. |
| `bargmann1961`, `hall1994` | Bargmann–Segal transform | CITE-OK | Standard references for the Bargmann–Segal transform; not separately verified but consistent with standard usage. |
| `jauchhill1940` | SU(3) symmetry of 3D HO, Phys. Rev. 57 641 (1940) | CITE-OK | Standard textbook reference. |
| `borelweil`, `Dyall`, `infeld_hull1951`, `cooper_khare1995` | Standard references | CITE-OK | All standard, consistent with usage. |
| `paper0`/`paper1`/`paper2_alpha`/`paper7`/`paper8_9`/`paper11`/`paper13`/`paper14`/`paper15`/`paper17`/`paper18`/`paper22`/`paper23`/`paper24`/`paper25`/`paper28`/`paper30`/`paper34`/`paper38` | Internal GeoVac papers | CITE-INTERNAL | Self-references; spot-checked Paper 22 numbers above (✓); spot-checked Paper 14 "51×–1712×" claim (does NOT carry the exact ratios — see overstatement above). |
| `curve_fit_audit` | Internal memo at `docs/curve_fit_audit_memo.md` | CITE-OK | Internal document, exists per CLAUDE.md §6. |
| `lorentzian_l0_memo` | Internal memo at `debug/lorentzian_l0_audit_memo.md` | CITE-OK | Internal document; bibitem path is correct. |
| `bizi_brouder_besnard2018` | "Towards a noncommutative geometry of the Standard Model with neutrino mixing," arXiv:1611.07062 (2018) | **CITE-WRONG-METADATA** | The arXiv ID is real and the year ranges of 2016–2018 are consistent (J. Math. Phys. 59:062303 (2018)), but the **actual title is "Space and time dimensions of algebras with applications to Lorentzian noncommutative geometry and quantum electrodynamics"** — verified online. Paper 31's title is a different paper / paraphrase. Severity MEDIUM. |
| `nieuviarts2025a` | "Emergence of Lorentz symmetry from an almost-commutative twisted spectral triple," arXiv:2502.18105, "A. Nieuviarts" | **CITE-WRONG-METADATA** (author initial) | Title is correct (verified via arxiv.org). **Author's first initial is G. (Gaston Nieuviarts), not A.** Severity LOW (cosmetic). |
| `nieuviarts2025b` | "Proceedings synthesis on twisted spectral triples and signature change," arXiv:2512.15450, "A. Nieuviarts" | **CITE-WRONG-METADATA** (author initial AND title) | arXiv:2512.15450 actually titled "Emergence of Time from a Twisted Spectral Triple in Almost-Commutative Geometry" (verified), submitted Dec 17, 2025. **Author's first initial is G. (Gaston), not A.** Paper 31's title "Proceedings synthesis on twisted spectral triples and signature change" is not the real title. Severity MEDIUM (paraphrased/incorrect title) + LOW (initial). |
| `hartle_hawking1976` | Phys. Rev. D 13, 2188 (1976) | CITE-OK | Verified. |
| `sewell1982` | Ann. Phys. 141, 201 (1982) | CITE-OK | Standard reference. |
| `bisognano_wichmann1976` | J. Math. Phys. 17, 303 (1976) | CITE-OK | Standard reference. |
| `unruh1976` | Phys. Rev. D 14, 870 (1976) | CITE-OK | Standard reference. |

### Problems found

- **HIGH — `perez_sanchez2024` misbundling.** Paper 31 has the EXACT same misbundling that Paper 32 was caught on yesterday: arXiv:2401.03705 derives Yang–Mills–**HIGGS** (verified abstract: "derive the Yang–Mills–Higgs theory on flat space as a smooth limit"); only arXiv:2508.17338 ("Comment on 'Gauge networks in noncommutative geometry'") gives the without-Higgs continuum limit ("the continuum limit of this theory is the Yang–Mills action functional, without a Higgs scalar"). Paper 31 bundles both arXiv IDs into one bibitem labeled `perez_sanchez2024`, gives the bibitem a synthesised title ("On the gauge networks of Marcolli–van Suijlekom: Yang–Mills without Higgs") that is not the actual title of either paper, and uses the citation in §II (line 200) and §I (line 163) to support the "without Higgs" claim — which is supported only by the 2025 paper. Paper 32 already has the correct split; Paper 31 needs the same fix.
- **MEDIUM — `bizi_brouder_besnard2018` wrong title.** The cited title "Towards a noncommutative geometry of the Standard Model with neutrino mixing" is not the title of arXiv:1611.07062, which is "Space and time dimensions of algebras with applications to Lorentzian noncommutative geometry and quantum electrodynamics" (J. Math. Phys. 59:062303 (2018)).
- **MEDIUM — `nieuviarts2025b` wrong title.** The cited title "Proceedings synthesis on twisted spectral triples and signature change" is not the title of arXiv:2512.15450, which is "Emergence of Time from a Twisted Spectral Triple in Almost-Commutative Geometry".
- **LOW — Nieuviarts author initial.** Both `nieuviarts2025a` and `nieuviarts2025b` give the author as "A. Nieuviarts"; the correct initial is G. (Gaston Nieuviarts).

### Priority / novelty claims

| Claim (verbatim) | Location | Searched | Prior art found? | Recommendation |
|---|---|---|---|---|
| "The decomposition is not novel as a piece of spectral triple formalism — it is the standard $\mathcal{A}$/$D$ split" | §II ¶7 | Connes 1994 + standard NCG | Yes — explicitly acknowledged | Verdict: appropriately scoped. No change. |
| "Every theorem cited is established in a previous paper of this series" (§I "Scope of contribution") | §I and §VI | Internal cross-check | Yes — the paper is honest about being framing | Verdict A. No change. |

The paper is appropriately humble about novelty — it explicitly says the partition framework is descriptive, not predictive, and disclaims new theoretical content. No false novelty claims to flag.

## Combined severity table

| Finding | Pass | Verdict | Severity |
|---|---|---|---|
| `perez_sanchez2024` misbundling — same defect Paper 32 was caught on | B | CITE-MISATTRIBUTED + CITE-DOESNT-SUPPORT | **HIGH** |
| §V three-layer framing stale vs corpus's four-layer (Paper 24 §V `subsec:asymmetry_layer4` + Paper 25 lines 956–963) | A | C (overstated via staleness; load-bearing per brief) | MEDIUM |
| §VI.A "51× to 1712× advantage" attributed to Paper 14 but Paper 14 says "two or more orders of magnitude" with different specifics | A/B | C / CITE-DOESNT-SUPPORT | MEDIUM |
| `bizi_brouder_besnard2018` wrong title | B | CITE-WRONG-METADATA | MEDIUM |
| `nieuviarts2025b` wrong title | B | CITE-WRONG-METADATA | MEDIUM |
| `nieuviarts2025a` and `nieuviarts2025b` wrong author initial (A. → G.) | B | CITE-WRONG-METADATA | LOW |
| Stale internal references `sec:universal_sector`, `sec:potential_specific_sector` (lines 876–877) | A | E (LaTeX, will render as ??) | LOW |
| §I "almost verbatim" overstates lineage match | A | C | LOW |

## Broadcast readiness: **YELLOW**

The paper is structurally sound — it explicitly frames itself as organizational, not as new theory, and that framing earns it broad clearance on the content side (most Pass A verdicts are A or B, not E). However, **broadcasting it in its current form would repeat the exact Perez-Sanchez 2024/2025 misbundling that Paper 32 was just fixed for** — the same single-bibitem-with-both-arXiv-IDs construction, the same false attribution of the "without Higgs" content to arXiv:2401.03705 (which actually derives Yang–Mills–Higgs in its abstract). This is HIGH severity because it's a citation-misattribution defect that an NCG reader will catch on first read (the Bratteli-networks paper is well known). The four-layer staleness in §V is MEDIUM because the brief explicitly flagged it as load-bearing, and the corpus already has the corrected version (Paper 24 §V, Paper 25 lines 956–963) — Paper 31 is the lagging instance. Two more wrong-title bibitems (BBB 2018, Nieuviarts 2025b) and two wrong-author-initial bibitems (Nieuviarts) are cosmetic but multiply with the Perez-Sanchez defect to suggest the bibliography was assembled from CLAUDE.md / memory without final arXiv verification. Fix the HIGH and MEDIUM items in one errata sprint and it goes GREEN.

## What I could NOT verify (hand to a human expert)

- The detailed content of the Marcolli–vS 2014 paper beyond its abstract (whether the inner-derivation / almost-commutative-extension framing of §II matches their specific construction in detail).
- The claim that Paper 22 verified bit-identical ERI patterns across five potentials is GeoVac-internal — no external comparison source exists.
- The two-body verification (§IX) Pearson correlations (0.41–0.58 for spectral action, 0.81/0.75 for Dirac resolvent) rest on Sprint 2026-06-01 memos; not independently reproduced.
- Whether the "M3 trivialization on chirality-symmetric spectrum" Sprint L2 prediction (§VIII.D) is correct — this is explicitly a prediction in the paper.
- Whether the BBB 2018 sign-table extension to KO-dim (3,1) GeoVac is consistent in detail (referenced as the Sprint L2 closure path).

---

## Summary line

- Paper title: **The Universal/Coulomb Partition of the GeoVac Spectral Triple**
- Verdict: **YELLOW**
- Pass A (A–E) counts: A=4 (Claims 1, 9-partial, 14, plus implicit framing), B=8 (Claims 2, 3, 4, 5, 7, 10, 12, 13), C=2 (Claims 6, 8), D=2 (Claims 9 sign-convention subpart, 11), E=0 (plus 2 LaTeX cross-reference errors at LOW severity)
- Pass B counts: CITE-OK ≈ 14, CITE-WRONG-METADATA = 3 (BBB title; Nieuviarts 2025a initial; Nieuviarts 2025b title+initial), CITE-MISATTRIBUTED + CITE-DOESNT-SUPPORT = 1 load-bearing (Perez-Sanchez 2024/2025 misbundle), CITE-CANT-FIND = 0
- HIGH = 1, MEDIUM = 4, LOW = 3
- **Top finding (one-liner): Paper 31 repeats the exact Perez-Sanchez 2024/2025 misbundling that Paper 32 was caught on yesterday — arXiv:2401.03705 derives Yang–Mills–HIGGS in its abstract, only arXiv:2508.17338 is the without-Higgs correction; Paper 31 needs the same two-bibitem split that was applied to Paper 32.**
