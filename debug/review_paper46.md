# Confidence Review: Paper 46 — Strong-form Lorentzian propinquity convergence on truncated SU(2)×U(1)_T Krein spectral triples

Reviewer: Confidence Reviewer (combined Pass A + Pass B).
Reviewed: 2026-06-02. Wave 2 Lorentzian-arc audit. Text-level only.
File: `papers/group1_operator_algebras/paper_46_strong_form_lorentzian_propinquity.tex` (2294 lines, 24 pages).

## Calibration check
Not a calibration run — no known-honest comparison answer staked. This is a substantive audit of a derivative-but-substantial paper in the math.OA Lorentzian arc.

## Pass A — Content audit

### Claim inventory + verdicts

| # | Claim (paraphrase) | Location | Verdict | Rests on | Evidence I produced |
|---|---|---|---|---|---|
| 1 | Main theorem: $\Lambda^{\rm strong} \le C_3^{\rm op}(n_{\max}) \cdot \gamma^{\rm joint} \to 0$ on the full Krein space under $L_{\rm op}$ | §6 Thm 4.1 (line 1105) | B | GEOVAC-ONLY (Papers 38/44/45 + five-lemma transport) | Five-lemma proof is internally consistent; rests entirely on Paper 45 cb-norm/Berezin and Paper 38 spatial L3 — neither externally re-verified here. |
| 2 | Closed form $C_3^{\rm op}(n_{\max}) = \sqrt{1 - 1/n_{\max}} \to 1^-$ | Eq (3.7) line 776; Lemma 4.1 line 766 | **A** | EXTERNAL (elementary calculus) | I recomputed: $\sup_{2\le N\le 2n_{\max}-1}\sqrt{(N-1)/(N+1)}$ saturates at $N=2n_{\max}-1$, giving $\sqrt{(2n_{\max}-2)/(2n_{\max})} = \sqrt{1-1/n_{\max}}$ exactly. Verified at $n_{\max}\in\{2,3,4,5,10,20\}$ to 6 digits. |
| 3 | Temporal-Lipschitz-invisibility identity (F1): $[D_L, a_s\otimes a_t] = i[D_{GV}, a_s]\otimes a_t$ bit-exact | Lemma 3.2 (line 591) | B | GEOVAC-ONLY (Paper 44 Prop 5.1 + chirality conventions) | The two algebraic ingredients ($[\gamma^0, M^{\rm spat}]=0$ on chirality-doubled lifts; $[\partial_t, M^{\rm temp}]=0$ on diagonal symbols) are claimed as Paper 44 Prop 5.1. I verified the second is trivial (commuting diagonals). The first depends on the specific chirality-doubled doubling convention. Internally consistent; no independent check. |
| 4 | Numerical panel: $\Lambda^{\rm strong}(2,3)=2.0746$, $(3,5)=1.6101$, $(4,7)=1.3223$ matches Paper 45 bit-exactly | Tab 1 line 1200 | B | GEOVAC-ONLY (driver `debug/l3b_2d_propinquity_assembly_compute.py`) | Not run by reviewer; values match Paper 45's published numbers verbatim. Bit-exactness is internal-consistency, not external truth. |
| 5 | Riemannian-limit recovery at $N_t=1$ is bit-exact (Frobenius residual 0.0 in float64) | Prop 6.1 line 1229 | B | GEOVAC-ONLY | Internal load-bearing falsifier; not independently re-run. Same pattern as Paper 45 Prop 6.2 — if it fails Paper 45 fails first. |
| 6 | Free-upgrade reading: same bound as Paper 45 on a strictly larger seminorm/state space | Rem 6.2 line 1143 | C → reword | GEOVAC-ONLY | Body shows $C_5^{\rm joint} = \max(1,1,C_3^{\rm op},0) = 1$. The "strictly larger" claim is true *as a statement-quantity* (different Lipschitz seminorm, different state space) but the *bound* is identical, so "no rate cost" is the right framing. Reads OK; mild risk of "free upgrade" sounding stronger than it is — see overstatement section. |
| 7 | $\gamma^{\rm joint}_{n_{\max},N_t,T} = O(\log n_{\max}/n_{\max} + T/N_t)$ asymptotic rate | §6 Thm 4.1 and L2/L4 | B | MIXED (Paper 38 Stein–Weiss closed-form + Paper 45 Fejér-on-circle) | Spatial rate $4/\pi$ is Paper 38 Lemma L2 (claimed verified there); temporal rate $O(T/N_t)$ is Fejér-on-$S^1$, classical (Katznelson). The MIXED label is because the spatial rate is GeoVac-internal even though the temporal half is textbook. |
| 8 | App A: Paper 45 envelope erratum (range $N\le 2n_{\max}-1$ not $N\le n_{\max}$) | App A lines 1428–1540 | C | EXTERNAL (algebra) + GEOVAC-ONLY (substrate claim) | The closed-form $\sqrt{1-1/n_{\max}}$ I verified externally. The substrate claim that natural-substrate multipliers reach $N=2n_{\max}-1$ rests on Paper 44 §3 (the Avery–Wen–Avery $\Delta n=\pm 1$ ladder). **STALENESS:** Paper 45 has already absorbed this erratum (`rem:envelope_v2` at line 842, plus erratum-pair eqs in Paper 45). Appendix A still reads "NOT applied to Paper~45 in this paper; we record it here for completeness" (line 1500–1501). Framing is now historically stale; the erratum *has been* applied. |
| 9 | App B: enlarged-substrate (strict-strong-form) bound $\Lambda^{\rm enlarged} = \Lambda^{\rm P45}$ bit-exact | Thm B.2 line 1854 | D | GEOVAC-ONLY | The proof is sketched, not fully derived. The crucial "gradient-norm absorption" step is given as a *sketch* (line 1763–1799) deferring the full derivation to "the $\beta$-L5 closure memo §6.4". A domain expert would want to see the full step. The numerical panel for the enlarged substrate is asserted, not displayed. The closure of G1' is **claimed** at theorem grade but is proof-sketched, not proof-complete in the paper itself. |
| 10 | "To our knowledge this is the first strong-form Lorentzian propinquity convergence theorem on truncated Krein spectral triples in the published math.OA literature" | Abstract line 198 | A (form) / D (substance) | EXTERNAL (no prior art found) | The hedge "To our knowledge" is correct per §13 honesty ceiling. I searched arXiv math.OA for "Lorentzian propinquity" and "Krein spectral triple propinquity" — no hits. Latrémolière 2025 (arXiv:2512.03573) is Riemannian non-compact, not Lorentzian/Krein. Cannot *confirm* novelty (per the agent rule); can confirm no public counterexample. |
| 11 | "$\prop_{\rm achievable}(\Op^L_{\rm enlarged}) = 1$" (enlarged substrate denser than natural) | App B eq (B.3) line 1713 | B | GEOVAC-ONLY | Cited to L3b-2a sub-sprint §4; not re-derived. The natural substrate has prop = 2 (Paper 44, Paper 32). The enlargement claim is plausible but not checked here. |
| 12 | $4/\pi$ rate constant = M1 Hopf-base measure signature; transports to joint setting | §6.3 line 1268 | B | GEOVAC-ONLY (Paper 38 §L2 + Paper 32 §VIII case-exhaustion) | This is a Paper 38/32 result imported; not re-derived. Per CLAUDE.md WH1 master Mellin engine, the M1 reading is plausible internal-consistency, not externally verified. |

### Numbers I recomputed

| Claim | Paper's figure | Independent reference | My recomputed value | Survives? |
|---|---|---|---|---|
| Closed form $\sqrt{1-1/n_{\max}}$ at $n_{\max}=2,3,4$ | 0.7071, 0.8165, 0.8660 | Elementary algebra | 0.707107, 0.816497, 0.866025 | YES — exact agreement to 6 digits. |
| $\sup_{N\le 2n_{\max}-1}(N-1)/\sqrt{N^2-1}$ saturates at $N=2n_{\max}-1$ | Implicit in eq 3.7 | Calculus (monotone increasing in $N$) | Verified monotone, saturates at $N=2n_{\max}-1$ | YES. |
| Tab 1 panel values (2.0746, 1.6101, 1.3223) | Paper 46 Tab 1 | Paper 45 Tab 1 | Verbatim match | INTERNALLY consistent (Paper 45 = Paper 46 same values). External recomputation would require running `debug/l3b_2d_propinquity_assembly_compute.py`. |
| Ratio $\Lambda(4,7)/\Lambda(2,3)=0.6374$ | line 1271 | $1.3223/2.0746$ | $0.6374$ | YES (computed). |

### Circularity map (GEOVAC-ONLY chains)

The paper's logical chain bottoms out almost entirely on prior GeoVac papers:

- **Paper 38 (Riemannian SU(2) propinquity)** supplies: (i) the spatial L2 cb-norm rate $4/\pi$; (ii) the spatial L3 per-harmonic Lichnerowicz bound $C_3^{(N)} = \sqrt{(N-1)/(N+1)}$; (iii) the spatial Berezin reconstruction. The closed form of $C_3^{\rm op}$ (verdict A) is the ONLY externally-verifiable substantive result in the paper. Every other content claim chains through Paper 38.
- **Paper 44 (Lorentzian operator-system substrate)** supplies: (i) the chirality-doubled scalar-multiplier substrate; (ii) Prop 5.1 ($[\gamma^0, M^{\rm spat}]=0$); (iii) propagation number $\prop = 2$. Note: "Paper 44 Prop 5.1" referenced at line 503 maps to what is actually `prop:krein_positivity_trivial` in Paper 44 (line 1011) — the numbering is approximate. Mostly a label fidelity issue (LOW).
- **Paper 45 (K⁺-weak-form predecessor)** supplies: (i) the L4(a–c,e) properties as seminorm-independent (Rem 4.5 line 1002); (ii) the joint cb-norm machinery (Lemma 4.2); (iii) the joint mass-concentration moment $\gamma^{\rm joint}$. The "free upgrade" reading is *defined* by equality to Paper 45's bound.
- **Paper 32 §VIII case-exhaustion theorem** supplies the M1 reading of $4/\pi$.
- **Paper 42, 43** supply: BW-α/BW-γ four-witness arc context (no direct theorem use).
- **Paper 47** is forward-referenced (G2 closure via norm-resolvent).

External chain endpoints exist for: Latrémolière 2018/2023 propinquity construction; Bożejko–Fendler central multipliers; van den Dungen 2016 Krein Dirac; BBB 2018 signature classification; Bhatia 1997 tensor-norm factorization. These are all standard.

**Headline circularity finding:** Paper 46 is structurally a transport theorem — its content is "Paper 45 closure under a different seminorm gives the same bound on a larger object." It does not (and does not claim to) introduce new external structure. The audit verdict on the strong-form reading therefore inherits the audit verdict on Paper 45 verbatim. If Paper 45's five-lemma transport is sound, Paper 46's is sound; if Paper 45 has a defect, Paper 46 inherits it.

### Overstatement findings

| Phrase | Suggested honest replacement | Reason |
|---|---|---|
| "We close … the strong-form leg" (abstract line 148) | "We prove … for the strong-form Lipschitz seminorm $L_{\rm op}$ on the natural chirality-doubled scalar-multiplier substrate" | "Close" overstates given Appendix B's G1' enlarged-substrate proof is sketch-grade and the strict-strong-form bound is asserted-with-sketch rather than fully derived. |
| "free upgrade phenomenon: a stronger theorem at no rate cost" (Rem 6.2 line 1164) | "the strong-form bounds a strictly larger Lipschitz seminorm on a strictly larger state space; the numerical bound equals the K⁺-weak-form bound" | "Free upgrade" is suggestive but a non-standard term that could mislead. The bound is the same; the *quantity bounded* is larger. Rephrase in literal terms once and then use the shorthand. |
| "strict-strong-form extension on the enlarged operator system" (App B title; closure assertion line 1957–1962) | "Sketch of strict-strong-form extension to the enlarged operator system" or move §B to a "follow-on" framing | App B's main theorem (Thm B.2) is proved via "Sketch" (line 1763) deferring the load-bearing absorption step to an internal memo. The claim "G1' is closed" (line 1957) overstates the rigor in this paper; it is closed *in the L3b-2f sprint memos*, not derived in the published text. |
| "Appendix A: … (NOT applied to Paper~45 in this paper; we record it here for completeness)" (line 1500) | Update to acknowledge Paper 45 has been amended | Stale statement: Paper 45 §sec contains `rem:envelope_v2` already absorbing the erratum. Either remove "NOT applied" or change to past tense. |
| "this is the first strong-form Lorentzian propinquity convergence theorem on truncated Krein spectral triples in the published math.OA literature" (abstract line 198) | Already hedged with "To our knowledge" — acceptable as written | The hedge is honored. No softening needed; cannot positively confirm absence (per §13 honesty ceiling). |

## Pass B — Citation and novelty

### Citation table

| \cite key | Bib entry says | Verdict | What I found |
|---|---|---|---|
| `bekka_harpe_valette2008` | Cambridge 2008 | CITE-OK | Standard monograph; not load-bearing here (cited only in Paper 45 context). |
| `bhatia1997` | Matrix Analysis, Springer 1997 | CITE-OK | Used for tensor-product operator-norm factorization (Cor 3.4 proof, Thm IV.2.6); claim supported. |
| `bisognano_wichmann1976` | J Math Phys 17 (1976), 303–321 | CITE-OK | Standard. |
| `bizi_brouder_besnard2018` | J Math Phys 59 (2018), 062303; arXiv:1611.07062 | CITE-OK | Verified arXiv:1611.07062 — Bizi, Brouder, Besnard, JMP 59, 062303 (2018). Title in bibitem says "application to Lorentzian noncommutative geometry and quantum electrodynamics" (singular) but the actual title says "applications" (plural). Cosmetic. |
| `bozejko_fendler1991` | Arch Math 57 (1991), 290–298 | CITE-OK (assumed) | Standard reference for Herz–Schur multipliers. Not web-verified individually; consistent with use. |
| `camporesi_higuchi1996` | J Geom Phys 20 (1996), 1–18; arXiv:gr-qc/9505009 | CITE-OK (assumed) | Used widely in GeoVac; consistent. |
| `connes_vs2021` | Comm Math Phys 383 (2021), 2021–2067; arXiv:2004.14115 | CITE-OK (assumed) | Standard. |
| `devastato_lizzi_martinetti2018` | J Math Phys 59 (2018), 112301 | CITE-OK (assumed) | Mentioned in concurrent-work footnote only. |
| `farsi_latremoliere2024` | preprint 2024 | CITE-INCOMPLETE | No arXiv ID given. Footnote-level only; low-stakes. |
| `farsi_latremoliere2025` | arXiv:2504.11715 | CITE-OK (assumed; Paper 47 cross-corpus check confirmed this entry) | Cited as Riemannian Dirac path; consistent. |
| `franco_eckstein2014` | Class Quantum Grav 30 (2013), 135007; arXiv:1212.5171 | CITE-WRONG-METADATA | Bibitem year is 2014 in key, 2013 in volume — pick one. Class. Quantum Grav. 30 was indeed 2013. LOW. |
| `geroch1967` | J Math Phys 8 (1967), 782–786 | CITE-OK | Standard. |
| `hartle_hawking1976` | Phys Rev D 13 (1976), 2188–2203 | CITE-OK | Standard. |
| `hekkelman_mcdonald2024`, `hekkelman_mcdonald2024b` | J Noncommut Geom / J Funct Anal "to appear"; arXiv:2403.18619 / arXiv:2412.00628 | CITE-OK (assumed) | Only cited in concurrent-work footnote. |
| `katznelson2004` | Harmonic analysis, CUP 2004 | CITE-OK | Standard. |
| `latremoliere_metric_st_2017` | Adv Math 404 (2022) Paper 108393; arXiv:1811.10843 | CITE-OK | Verified arXiv:1811.10843 — Latrémolière, "The Gromov-Hausdorff propinquity for metric Spectral Triples", Adv Math 404 (2022), Paper 108393. Bib-key name "...2017" is misleading (preprint was 2018, published 2022), but this is consistent with cross-corpus usage and was raised+resolved in CLAUDE.md prior fixes. **LOW (cosmetic key)**. |
| `latremoliere2018` | Trans AMS 368 (2016), 365–411 | CITE-WRONG-METADATA | Bib key says "2018" but the journal published it in 2016 (the *Gromov-Hausdorff propinquity*, not metric ST). Brief notes this was already fixed corpus-wide as an OK convention. **LOW (cosmetic key)**. |
| `latremoliere2025_hypertopology` | arXiv:2512.03573, Dec 2025 | CITE-OK | Verified — Frederic Latrémolière, "The quantum Gromov-Hausdorff Hypertopology on the class of pointed Proper Quantum Metric Spaces", arXiv:2512.03573, submitted Dec 3, 2025. |
| `leimbach_vs2024` | Adv Math 439 (2024), Paper 109496 | CITE-OK (assumed; brief noted prior corpus fix verified) | Concurrent-work footnote only. |
| **`mondino_samann2024`** | **Lett Math Phys 114 (2024), Paper 37; arXiv:2209.14384** | **CITE-MISATTRIBUTED** | **arXiv:2209.14384 is Minguzzi & Suhr, "Lorentzian metric spaces and their Gromov–Hausdorff convergence", Lett Math Phys 114 (2024) art. 73 — NOT Mondino-Sämann.** Same misattribution as Paper 45 (flagged in audit brief as open item). Mondino-Sämann do have a 2024 work, but it is `mondino_samann2025` (arXiv:2504.10380, "Lorentzian Gromov-Hausdorff convergence and pre-compactness", 2025) — and the citation `mondino_samann2024` in Paper 46 should either (a) be retitled to Minguzzi-Suhr with corrected authors, or (b) replaced with an actual Mondino-Sämann 2024 (e.g., the Kunzinger-Sämann 2018 Lorentzian length spaces paper, or one of the Mondino-Sämann ESI/preprint papers). **HIGH severity** if propagated to readership; **MEDIUM severity** in practice because Paper 46 only cites this in the concurrent-work footnote (line 381) and §7.7 Mondino-Sämann discussion (line 1393), not in any load-bearing proof step. |
| `mondino_samann2025` | arXiv:2504.10380, 2025 | CITE-OK | Verified — Mondino & Sämann, "Lorentzian Gromov-Hausdorff convergence and pre-compactness", arXiv:2504.10380, submitted April 14, 2025. |
| `che_perales_sormani2025` | arXiv:2510.13069, Oct 2025 | CITE-OK | Verified — Che, Perales, Sormani, "Gromov's Compactness Theorem for the Intrinsic Timed-Hausdorff Distance", arXiv:2510.13069, October 2025. Audit brief noted this is CORRECT in Paper 46; confirmed. |
| `nieuviarts2025a` | arXiv:2502.18105 v3, May 2025 | CITE-OK (assumed; cross-corpus checked) | Concurrent-work footnote only. |
| `nieuviarts2025b_proceedings` | arXiv:2512.15450, Dec 2025 rev May 2026 | CITE-OK (assumed; cross-corpus checked) | Concurrent-work footnote only. |
| `pier1984` | Amenable Locally Compact Groups, Wiley 1984 | CITE-OK | Standard; not load-bearing here. |
| `pisier2001` | Similarity Problems and CB Maps, LNM 1618, Springer 2nd ed 2001 | CITE-OK | Standard; used for cb-norm context. |
| `sewell1982` | Ann Phys (NY) 141 (1982), 201–224 | CITE-OK | Standard. |
| `strohmaier2006` | J Geom Phys 56 (2006), 175–195 | CITE-OK (assumed) | Standard pseudo-Riemannian geometry reference. |
| `takesaki1979` | Theory of Operator Algebras I, Springer 1979 | CITE-OK | Standard. |
| `toyota2023` | arXiv:2309.13469, 2023 | CITE-OK (assumed) | Concurrent-work footnote only. |
| `unruh1976` | Phys Rev D 14 (1976), 870–892 | CITE-OK | Standard. |
| `vandungen2016` | Math Phys Anal Geom 19 (2016), Art 4; arXiv:1505.01939 | CITE-OK | Verified arXiv:1505.01939 — van den Dungen, "Krein spectral triples and the fermionic action", MPAG 19 (2016). Cannot confirm Prop 4.1 in particular without full text access, but the paper does have a Proposition structure consistent with the Krein-Dirac claim (Paper 43, Paper 45 reuse same citation). |
| `de_groot2026_su11` | arXiv:2601.22171, Jan 2026 | CITE-OK | Verified — Jort de Groot, "Pseudo-Riemannian Spectral Triples for SU(1,1)", arXiv:2601.22171, January 2026. |
| `paper24`, `paper32`, `paper38`, `paper39`, `paper40_unified`, `paper42`, `paper43`, `paper44`, `paper45`, `paper47` | GeoVac internal preprints, DOI via Zenodo | CITE-OK (internal) | All present in `papers/group1_operator_algebras/` and `papers/group3_foundations/`. Internal references; treated as "Paper X says so" not as proof. |
| `memo_L3b_2a/b/c/d` | internal sprint memos | CITE-OK (internal) | Memo files exist in `debug/`. |

### Problems found

**CITE-MISATTRIBUTED (1)**:
- **`mondino_samann2024` (line 2114–2118 bibitem)** points to arXiv:2209.14384, which is in fact Minguzzi & Suhr (verified via arxiv.org). This is the same misattribution flagged in Paper 45's audit (open item per audit brief). Affected in-text uses:
  - Line 381 (concurrent-work footnote): `Mondino--S\"amann~\cite{mondino_samann2024,mondino_samann2025}` — concurrent-work paragraph names Mondino-Sämann twice; only the 2025 cite is correct.
  - In-text `\cite{mondino_samann2024}` reference count: **1 use** (line 381 footnote only). The Mondino-Sämann discussion in §7.7 (line 1393) cites only `mondino_samann2025`. So Paper 46 has lighter dependency on the broken cite than Paper 45 does.

**CITE-WRONG-METADATA (2)**:
- `franco_eckstein2014`: key year says 2014, journal year says 2013. Pick one.
- `bizi_brouder_besnard2018`: bibitem title says "application" singular, actual title says "applications" plural.

**CITE-INCOMPLETE (1)**:
- `farsi_latremoliere2024`: no arXiv ID nor volume given (says "preprint 2024"). Concurrent-work footnote only — low stakes.

### Priority / novelty claims

| Claim (verbatim) | Location | Searched | Prior art found? | Recommendation |
|---|---|---|---|---|
| "To our knowledge this is the first strong-form Lorentzian propinquity convergence theorem on truncated Krein spectral triples in the published math.OA literature" | Abstract line 198 | arXiv (math.OA, math-ph), keywords "Lorentzian propinquity" / "Krein spectral triple propinquity" / "Lorentzian Gromov-Hausdorff truncated Krein" | No counterexample found. Latrémolière 2025 (arXiv:2512.03573) is Riemannian non-compact. Mondino-Sämann 2025 is synthetic/pre-length-space, not operator-algebraic. | KEEP — hedge is appropriate per §13 honesty ceiling. Cannot confirm novelty; can confirm no public counterexample. |
| "this is the first strong-form Lorentzian propinquity convergence theorem … closing the named open problem in Paper~45 §1.4 G1" | Abstract line 198–201 | Same as above + Paper 45 G1 status | Paper 45's open-problem status was a self-declared internal label (G1); the "first … in published math.OA literature" claim is what's externally verifiable and is appropriately hedged. The G1-closure is internal-consistency. | KEEP, but note that "closes G1" is a GeoVac-internal claim, not an externally validated priority claim. |

## Combined severity table

| Finding | Pass | Verdict | Severity |
|---|---|---|---|
| `mondino_samann2024` is actually Minguzzi-Suhr (arXiv:2209.14384) | B | CITE-MISATTRIBUTED | **HIGH** (CLAUDE.md §3 Fursaev-Solodukhin pattern; though here used only in footnote) |
| App B "Theorem B.2 closes G1'" overstates: Thm B.2 is proved via "Sketch" deferring load-bearing step to internal memo | A | C / D | **MEDIUM** |
| App A framing stale: says Paper 45 not yet amended, but Paper 45 has `rem:envelope_v2` already | A | C | **MEDIUM** (factual misstatement about state of corpus) |
| "Free upgrade" phrasing is suggestive, may mislead | A | C | **LOW** |
| Paper 44 "Prop 5.1" reference label fidelity (numbering offset) | A | C | **LOW** |
| `franco_eckstein2014` year mismatch (key 2014 vs vol 2013) | B | CITE-WRONG-METADATA | **LOW** |
| `bizi_brouder_besnard2018` title typo "application" vs "applications" | B | CITE-WRONG-METADATA | **LOW** |
| `farsi_latremoliere2024` incomplete bibitem (no arXiv ID) | B | CITE-INCOMPLETE | **LOW** |
| `latremoliere_metric_st_2017` and `latremoliere2018` key-year cosmetic mismatch | B | LOW (per CLAUDE.md prior fix decision) | **LOW** |

HIGH: 1. MEDIUM: 2. LOW: 6.

A-count: 1. B-count: 7. C-count: 4 (one phrase overstatement + three framing/staleness). D-count: 1 (App B proof rigor + novelty caveat). E-count: 0 (no math errors found).

CITE-OK: ~25. CITE-MISATTRIBUTED: 1. CITE-WRONG-METADATA: 2. CITE-INCOMPLETE: 1. CITE-CANT-FIND: 0.

## Broadcast readiness: **YELLOW**

Paper 46 is in good mathematical shape. The closed-form constant $\sqrt{1-1/n_{\max}}$ is the only externally re-verifiable substantive content; I confirmed it to 6 digits and verified the supremum saturation algebra. The body of the paper executes a clean five-lemma transport from Paper 45 under the strong-form seminorm $L_{\rm op}$ with the temporal-Lipschitz-invisibility identity (Lemma 3.2) as the load-bearing analytical input; all of this is internally consistent and rests transparently on Papers 38, 44, 45. The Riemannian-limit recovery falsifier (Prop 6.1) is set up correctly. The "free upgrade" phenomenon is real-as-defined: same numerical bound on a strictly larger Lipschitz seminorm/state space.

The YELLOW (not GREEN) verdict has three drivers, in descending importance:

1. **The `mondino_samann2024` misattribution** (arXiv:2209.14384 is Minguzzi-Suhr, not Mondino-Sämann) is a true citation error of the same class as the Fursaev-Solodukhin failure logged in CLAUDE.md §3. It appears only in a concurrent-work footnote in Paper 46 (1 in-text use), so the load-bearing impact is small, but the *same defect* is present in Paper 45 and must be fixed corpus-wide before broadcast. Recommended fix: replace `mondino_samann2024` cite key entirely with a corrected `minguzzi_suhr2024` entry pointing to arXiv:2209.14384, and remove or relocate the cite in the concurrent-work paragraph.

2. **Appendix B's G1' closure is sketch-grade, not proof-complete.** Paper 46 announces in §1.2 (F3) and Appendix B that the strict-strong-form / enlarged-substrate problem (originally named in §7.2 as "the deepest open question on the Lorentzian side") is now CLOSED. But the L3 transport on the enlarged substrate (line 1763) is given as a *Sketch* with the load-bearing absorption step deferred to "the β-L5 closure memo §6.4". A domain reviewer (Latrémolière, vS, BBB-side referee) would likely flag this and ask for the explicit derivation. Recommended fix: either expand the sketch into a full proof in App B, or reframe App B as "Proof outline for strict-strong-form extension; full derivation deferred to internal sprint memos" and downgrade the "closure" language.

3. **Appendix A's "Paper 45 not yet amended" framing is stale.** Paper 45 (line 842 `rem:envelope_v2`) has already absorbed the envelope erratum. Recommended fix: update the parenthetical in Paper 46 line 1500 from "(NOT applied to Paper~45 in this paper; we record it here for completeness)" to "(applied to Paper~45 as `rem:envelope_v2` 2026-05-22; we record the derivation here for the present paper's L3 closure)".

No math errors found. No structural defect in the main theorem. Paper 46's transport theorem is internally sound on the natural substrate; broadcast risk is concentrated in (1) the citation defect (HIGH) and (2) the G1' overclaim (MEDIUM).

## What I could NOT verify (hand to a human expert)

- Whether van den Dungen 2016 Proposition 4.1 is *literally* the construction $D_L = i(\gamma^0 \otimes \partial_t + D_{GV} \otimes I_{N_t})$ as cited (line 158, line 453). Bib entry is correct, but full text of Prop 4.1 not visible to me. A domain expert with the paper open should re-check that the cited proposition is in fact what's being used.
- Whether the L2 cb-norm closed form `cbnorm{S_K^joint} = 2/(nmax+1)` from Paper 45 Lemma 4.2 (transported into Lemma 4.2 here) is itself correctly derived in Paper 45. This audit takes Paper 45's L2 as given.
- Whether the "Avery–Wen–Avery $\Delta n = \pm 1$ ladder on $n_{\max}$ shells produces shell-difference labels $N$ up to $2 n_{\max} - 1$" is correctly stated in Paper 44 §3. Plausible per CLAUDE.md substrate description; not re-derived here.
- Whether `debug/l3b_2d_propinquity_assembly_compute.py` produces the panel values 2.0746 / 1.6101 / 1.3223 as claimed. Driver path consistent with sub-sprint memo references; not run by reviewer.
- Whether the original L3b-2f $\beta$-L5 closure memo (referenced as "$\beta$-L5 closure memo §6.4" but no `\bibitem` entry) contains a fully rigorous derivation. The four sub-sprint memos for L3b-2a/b/c/d ARE in `\bibitem` entries; the $\beta$-L5 memo is referenced in text but lacks a `\bibitem` — this is a LOW-severity citation gap (the memo should be added to the bibliography if it is the load-bearing source for Theorem B.2's proof).

---

**Summary delivered to dispatcher.**

Title: Paper 46 — Strong-form Lorentzian propinquity convergence on truncated SU(2)×U(1)_T Krein spectral triples.
Verdict: **YELLOW**.
Pass A counts: **A=1, B=7, C=4, D=1, E=0.**
Pass B counts: **CITE-OK ~25, CITE-MISATTRIBUTED 1, CITE-WRONG-METADATA 2, CITE-INCOMPLETE 1, CITE-CANT-FIND 0.**
Severity totals: **HIGH=1, MEDIUM=2, LOW=6.**
Top finding: **`mondino_samann2024` bibitem points to arXiv:2209.14384, which is in fact Minguzzi & Suhr, "Lorentzian metric spaces and their Gromov–Hausdorff convergence" (Lett Math Phys 114, 2024, art. 73) — same CITE-MISATTRIBUTED defect as Paper 45, requiring corpus-wide fix before broadcast.**
