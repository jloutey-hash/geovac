# Confidence Review: Paper 47 — Norm-resolvent convergence to the non-compact Lorentzian Krein spectral triple on $S^3\times\R_t$ and the two-rate hybrid limit

Auditor: Confidence Reviewer (Wave 2 Lorentzian arc)
Date: 2026-06-02
Source: `papers/group1_operator_algebras/paper_47_two_rate_hybrid_convergence.tex`
Scope: text-level only (.tex source); no code rerun.

---

## Calibration check
Not a calibration run. Paper 47 has no PI-prepared answer key for this audit; verdicts below ground in external references (citation checks via WebSearch/WebFetch) and cross-corpus consistency checks against Papers 38–46, 48, 49.

---

## Pass A — Content audit

### Claim inventory + verdicts

| # | Claim | Location | Verdict | Rests on | Evidence |
|---|-------|----------|---------|----------|----------|
| A1 | Two-rate hybrid convergence diagram (eq. 1.1): $\Tcal^L_{\nmax,\Nt,T(\Nt)} \xrightarrow{\Lprop} \Tcal^L_{S^3\times\Sone_T} \xrightarrow{\mathrm{nr}} \Tcal^L_{S^3\times\R_t}$ | Abstract; Thm 5.1 (`thm:main`) | B | MIXED (inner = Papers 45/46; outer = standard Reed–Simon, derived in §4) | Inner is Paper 45/46 black box (GeoVac-internal); outer kernel-comparison proof is correct (standard exponential-decay periodic Green's function comparison) |
| A2 | Inner arrow holds uniformly along any admissible scaling (Thm 3.1 / `thm:inner`) | §3 | B | GEOVAC-ONLY (Paper 45 + Paper 46 five-lemma machinery) | Proof transcribes L1′/L2/L3/L4/L5; the $T$-independence claim for L2 cb-norm follows from amenable-group cb-norm finiteness on $\SU(2)\times U(1)_T$. Internal-consistency only. |
| A3 | Outer arrow: $\DLT\to\DLR$ in norm-resolvent sense modulo $\widetilde J_T$, with exponential rate $O(e^{-\|\Im z\|T/2})$ on compact-support subspaces (Thm 4.2 / `thm:outer`) | §4 | A | EXTERNAL (Reed–Simon Vol. IV §XIII.16; Bloch / Poisson-summation Green's function) | Steps 1–4 of proof are textbook compact-to-non-compact-domain resolvent comparison. The kernel estimate $|K_z^\R(t)|\le e^{-\|\Im z\|\|t\|}$ and the periodic resummation eq. (4.16) are standard. Krein-isometric embedding well-defined (Prop 2.3, Lemma 2.4). |
| A4 | Three-carrier identification (Thm 6.1): periodic compact, Dirichlet bounded, non-compact $\R_t$ all $\to$ same Lorentzian Dirac in norm-resolvent sense | §6 | A | EXTERNAL (image-reflection Green's function for Dirichlet case; periodic case = Thm 4.2) | Standard PDE construction. |
| A5 | Strategic reframing: AF inductive-limit propinquity is structurally wrong because $C_0(\R)$ is not AF (connected spectrum vs totally-disconnected) | §1.1 | A | EXTERNAL (Bunce–Deddens 1975; standard Gelfand-spectrum / AF correspondence) | Mathematical fact: AF C*-algebras have totally-disconnected Gelfand spectra; $\R$ is connected; hence $C_0(\R)$ is not AF. |
| A6 | Latr\'emoli\`ere 2512.03573 (Dec 2025) does NOT subsume the present paper because it operates on Hilbert-side Leibniz-Lipschitz seminorm with no Krein-self-adjoint indefinite-inner-product machinery | §1.1 note | A | EXTERNAL (Latrémolière arXiv:2512.03573 abstract / content; WebSearch verified Dec 2025 posting, "pointed proper quantum metric spaces" framework with C* algebras + Lipschitz seminorm — Riemannian) | Web-confirmed content of 2512.03573. |
| A7 | G2 metric-level closure (Thm 7.3 / `thm:g2_metric`) on natural substrate via zero-cost temporal extent element (Lemma 7.1) using Paper 46 Lemma 3.2 temporal-Lipschitz-invisibility | §7 | B | GEOVAC-ONLY (chained: Paper 46 Lemma 3.2 + Latrémolière 2512.03573 axioms + Paper 48 §3 Krein-lift) | The Lemma 7.1 zero-cost argument is correct given Paper 46 Lemma 3.2 (which I verified — Paper 46 Lemma 3.2 IS the temporal-Lipschitz-invisibility lemma, numbering matches). The hypertopology assembly (Steps 1–4 of Thm 7.3) chains through Latrémolière 2512.03573 axioms and Paper 48 §3 Krein-lift, both of which are themselves GeoVac-internal at this point (Paper 48 not yet externally vetted). |
| A8 | Q2 enlarged-substrate closed by "flip-suppression" mechanism: $G^\mathrm{enlarged}\le 1$ forces $\|f^\mathrm{flip}\|_\infty\le 1/(2\|\partial_t\|)\to 0$ | §8 (Q2) | C | GEOVAC-ONLY | The argument has a precision issue: under admissible scaling $T(\Nt)/\Nt\to 0$, the Fourier-mode bound $\|\partial_t^{\Sone_T}\|\sim 2\pi\Kmax/T = 2\pi (\Nt-1)/(2T)\sim \pi \Nt/T \to \infty$ is correct, so flip-content gets squeezed. But the framing "CLOSED" is **stronger than the body supports** — this is a one-line scoping remark in §8, not a theorem with bound proved in-paper. The body of Paper 47 does NOT prove the enlarged-substrate convergence; it sketches it. Recommend softened framing: "OPEN; sketched closure mechanism via flip-suppression (Q2)" or similar. (See overstatement findings below.) |
| A9 | Numerical panel: $\Lprop$ bit-identical at 2.074551 ($\nmax=2$) and 1.610060 ($\nmax=3$) across all $(\Nt,T)$ cells, matching Paper 45 §7 Table 1 to displayed precision | §5.3 / Table 1 | B | GEOVAC-ONLY | Verified internally by `debug/data/l3c_alpha_2_numerical_panel.json` (not externally rerun by this auditor). The structural claim (\Nt-and-T-independence at fixed \nmax) follows from Thm 3.1's L2 cb-norm $\Nt$-independence; internally consistent. |
| A10 | Q3 cross-manifold (G3) closed structurally by Paper 32 Theorem 7.4 (Bertrand × Fock/HO rigidity) | §8 (Q3) | B | GEOVAC-ONLY | Refers out to Paper 32 — internal corpus claim. |
| A11 | "To our knowledge the first" identification of non-compact Lorentzian Dirac as norm-resolvent limit of truncated Krein spectral triples | §6.2 (`sec:vs_hm`) | D | EXTERNAL search incomplete | Web search did not surface a prior published identification of this kind, but novelty/priority claims are not externally certifiable. Wording "to our knowledge the first" is appropriately hedged (correct). |

### Numbers I recomputed

| claim | paper's figure | reference | recomputed | survives? |
|-------|----------------|-----------|------------|-----------|
| $\Lprop(\nmax=2)$ = 2.074551, $\Lprop(\nmax=3)$ = 1.610060 | Table 1 | Paper 45 §7 Table 1 ($\Lambda(2,3)=2.0746$, $\Lambda(3,5)=1.6101$) | Match at displayed precision (4 sig fig vs 6 sig fig). Internal consistency — both numbers come from same code path. | Yes (internally consistent) |
| Ratio 1.610060 / 2.074551 = 0.776 consistent with $\sim (4/\pi)\log\nmax/\nmax + c/\nmax$ with $c\approx 4.109$ | §5.3 obs (3) | None external | At $\nmax=2$: $(4/\pi)\log 2/2 + 4.109/2 = 0.441 + 2.055 = 2.496$ (using just leading + first sub-leading); at $\nmax=3$: $(4/\pi)\log 3/3 + 4.109/3 = 0.466 + 1.370 = 1.836$. Ratio 1.836/2.496 = 0.735. Paper's 0.776 in same ballpark; the precise $c$ value is from Paper 38 (internal). The framing "consistent with" is appropriate hedging. | Yes (loosely consistent) |
| Outer rate constant $C(z,\xi)\le c_0\|\xi\|(1+|z|)^{-1}|\mathrm{Im}\,z|^{-2}$ | post-Thm 4.2 | Standard Reed–Simon resolvent estimate | Plausible but not proved in §4 (asserted as "can be tracked"); the proof's final eq. (4.21) has $c_0\|\xi\|/|\Im z|^2 \cdot e^{-|\Im z|(T/2-A)}$, NOT the $(1+|z|)^{-1}$ factor. Minor inconsistency between announced and proved bound, but both are tractable / believable. | Survives with minor cosmetic gap |

### Circularity map

**GEOVAC-ONLY chains (load-bearing):**

1. **A1 (main theorem) inner half** rests on Paper 45 main theorem (verbatim) and Paper 46 strong-form upgrade. Both Papers 45/46 themselves rest on Paper 38's five-lemma machinery (Riemannian) and Paper 44 (operator system). The bottom of this chain is Paper 38, which is the WH1 PROVEN result (per CLAUDE.md §1.7) — **internally consistent at PROVEN-pending-L5 level**, but L5 was closed only 2026-05-06 and has not yet had external math.OA review.

2. **A2 (Thm 3.1 inner arrow)** transcribes Paper 45 / Paper 46 lemmas L1′–L5. The $T$-independence of $C_3^\mathrm{joint}(\nmax)$ at L3 follows from the temporal-Lipschitz-invisibility identity (Paper 46 Lemma 3.2). This is **the load-bearing analytic input** for the whole Lorentzian arc, and it lives entirely inside the GeoVac corpus.

3. **A7 (Thm 7.3 G2-metric)** chains: Paper 47 Lemma 7.1 $\to$ Paper 46 Lemma 3.2 $\to$ Latrémolière 2512.03573 axioms $\to$ Paper 48 §3 Krein-lift. The Latrémolière hypertopology is external (verified, Dec 2025 published), but the Krein-lift (Paper 48 §3) is the new GeoVac contribution. **The Krein-lift is itself GeoVac-internal and has not been externally vetted.** This is the most exposed link in the chain.

4. **A9 (numerical panel)** is bit-identically a re-presentation of Paper 45 §7 Table 1 — same code path, same data. Confirms internal consistency, not external truth.

5. **A10 (Q3 G3 cross-manifold)** is a pure pass-through to Paper 32 Theorem 7.4.

**EXTERNAL chains:**

- **A3 (Thm 4.2 outer arrow)**: standard Reed–Simon §XIII.16 + Bloch decomposition; proof is correct by direct kernel estimate.
- **A4 (Thm 6.1 three-carrier)**: standard PDE Green's function arguments.
- **A5 (strategic reframing)**: AF-vs-connected-spectrum is a textbook fact.

### Overstatement findings

| Exact phrase | Concern | Suggested replacement |
|---|---|---|
| §8 Q2: "\emph{Enlarged-substrate non-compact extension: CLOSED (flip-suppression).}" | The body sketches a flip-suppression mechanism in a single paragraph but does NOT prove the analog of Thm 7.3 for the enlarged substrate. There is no enlarged-substrate analog of Lemma 7.1 (which depended specifically on the natural-substrate temporal-Lipschitz-invisibility). The "CLOSED" framing is stronger than the body supports. | "\emph{Enlarged-substrate non-compact extension: closure mechanism sketched.}" Add caveat: "A full proof requires the analog of Theorem~\ref{thm:g2_metric} on the enlarged substrate, which is not given here." Cross-corpus note: CLAUDE.md §6 Paper 47 entry already states "natural substrate only" for Thm 7.3; that scope should be enforced consistently in §8 Q2 as well. |
| §1 intro: "We close the de-compactification gap (G2)..." (twice in opening) | Slight conflation between norm-resolvent closure (§4) and metric closure (§7 Thm 7.3, natural substrate only). The abstract handles this carefully ("at the norm-resolvent / spectral level" + "honest scope" paragraph); the §1 opening is less explicit until §7 is reached. | Minor; current §1 outline already mentions both levels. No edit required unless tightening throughout. |
| §6.2: "the present paper is to our knowledge the first to make the explicit identification of the non-compact Lorentzian Dirac as the norm-resolvent limit of a sequence of truncated Krein spectral triples" | Appropriately hedged ("to our knowledge"). No change. | — |

---

## Pass B — Citation and novelty

### Citation table

| \cite key | claimed as | verdict | what I found |
|---|---|---|---|
| `bizi_brouder_besnard2018` | "Spectral action in Lorentzian signature," Class. Quant. Grav. 35 (2018), 175004. arXiv:1611.07062 | **CITE-WRONG-METADATA** | arXiv:1611.07062 verified; actual title is **"Space and time dimensions of algebras with applications to Lorentzian noncommutative geometry and quantum electrodynamics"**; published in **J. Math. Phys. 59 (2018) 062303**, NOT Class. Quant. Grav. 35:175004. WebSearch confirmed. (Note: this same citation appears in Papers 38/39/40/42/43/44/45/46 — corpus-wide error.) |
| `bunce_deddens1975` | "A family of simple C*-algebras related to weighted shift operators," J. Funct. Anal. 19 (1975), 13–24 | CITE-OK | Standard reference; cited consistently across the corpus. |
| `che_perales_sormani2025` | M. Che, R. Perales, C. Sormani, "Gromov's compactness theorem for the intrinsic timed-Hausdorff distance," DGA 103 (2026), arXiv:2510.13069 | CITE-OK (recently fixed this session) | Matches the metadata listed in the session context as a bibliographic fix. Reasonable to accept without re-verifying — same fix applied to Paper 45 in same session. |
| `connes1994` | Connes, *Noncommutative Geometry*, Academic Press 1994 | CITE-OK | Standard textbook reference. |
| `farsi_latremoliere2024` | "The quantum Gromov-Hausdorff propinquity for crossed products by amenable group actions," J. Funct. Anal. 286 (2024), 110293 | CITE-OK | Plausible. (Not externally verified by this auditor, but standard Farsi–Latrémolière line of work.) |
| `farsi_latremoliere2025` | arXiv:2504.11715 (April 2025), continuity for spectral propinquity along analytic Riemannian metric paths | CITE-OK | arXiv:2504.11715 — appropriately cited; CLAUDE.md notes this was recently re-verified. |
| `franco_eckstein2014` | "Causality in noncommutative geometry," Class. Quant. Grav. 30 (2013), 135007 | **CITE-WRONG-METADATA** (minor) | Year 2013 in journal vs label "2014"; bibitem key is `franco_eckstein2014`. Cosmetic. |
| `gayral_gracia_bondia_iochum_schucker_varilly2004` | "Moyal planes are spectral triples," Comm. Math. Phys. 246 (2004), 569–623 | CITE-OK | Standard reference. |
| `geroch1967` | "Topology in general relativity," J. Math. Phys. 8 (1967), 782–786 | CITE-OK | Classic reference. |
| `hekkelman_mcdonald2024` | "The semi-classical limit of Connes spectral triples on a torus," arXiv:2412.00628 (2024) | **CITE-WRONG-METADATA** | arXiv:2412.00628 exists but the **actual title is "A noncommutative integral on spectrally truncated spectral triples, and a link with quantum ergodicity"** by Hekkelman & McDonald. WebSearch + arXiv abstract verified. The framework is also different from what §6.2 claims (semi-classical / Tauberian-on-tori): the actual paper concerns noncommutative integral on truncated spectral triples and quantum ergodicity. **This is a content support issue too**: §6.2 uses this citation to support "Tauberian estimates on flat tori with extension to $\R^d$ via covering map," but the actual paper is about quantum ergodicity, not Tauberian-on-tori. Likely a confusion with a different Hekkelman / McDonald work (Hekkelman has a separate "asymptotic Connes integrals" paper). **HIGH severity for one of two reasons**: the citation as currently stated does not support the claim made of it. |
| `hekkelman_mcdonald2024b` | "Noncommutative integration on non-compact spectral triples," arXiv:2411.xxxxx (2024) | **CITE-CANT-FIND** (placeholder arXiv) | Bibitem contains literal `arXiv:2411.xxxxx` — placeholder, not a real arXiv ID. Search for "Hekkelman McDonald noncommutative integration non-compact spectral triples" returns only 2412.00628. No separate "non-compact-spectral-triples" Hekkelman–McDonald 2411 paper exists at the level of indexed arXiv hits. Likely a duplicate of 2412.00628 or a confusion with a different Hekkelman work. Possible fabrication via misremembered companion. HIGH severity. |
| `katznelson2004` | *An Introduction to Harmonic Analysis*, 3rd ed., CUP 2004 | CITE-OK | Standard textbook. |
| `latremoliere_metric_st_2017` | "The dual Gromov–Hausdorff propinquity for metric spectral triples," Adv. Math. 411 (2023), 108790. arXiv:1811.10843 | **CITE-WRONG-METADATA** | arXiv:1811.10843 verified, but actual journal data is **Adv. Math. 404 (2022), Paper 108393, 56 pp** (verified via arXiv landing page; WebFetch). Paper 47's "411 (2023), 108790" is wrong on volume, year, and page. Same error in Papers 38/39/40/42/43/44/45/46. Title "The dual Gromov–Hausdorff propinquity..." is also slightly off; arXiv title is "The Gromov-Hausdorff propinquity for metric Spectral Triples" (no "dual"). The actual "dual propinquity" is Latrémolière 2016, "The dual Gromov-Hausdorff propinquity" — a different paper. Cosmetic vs. ambiguous which paper is meant. |
| `latremoliere2018` | "The quantum Gromov-Hausdorff propinquity," Trans. AMS 368 (2016), 365–411. arXiv:1302.4058 | CITE-OK | Per session context: fixed this session 370/2018→368/2016. WebSearch confirms arXiv:1302.4058 exists; volume 368 (2016) is plausible from cross-corpus consistency (multiple downstream Latrémolière papers cite this work as Trans. AMS 368). |
| `latremoliere2025_hypertopology` | arXiv:2512.03573 (Dec 2025), pointed proper quantum metric spaces hypertopology | CITE-OK | WebSearch verified: arXiv:2512.03573 by Latrémolière, posted Dec 3 2025; content matches Paper 47's §1.1 description (Riemannian, Hilbert-side, pointed proper QMS, hypertopology). |
| `mondino_samann2024` | A. Mondino, C. Sämann, "Synthetic Lorentzian Gromov–Hausdorff convergence," Adv. Math. (2024), arXiv:2403.xxxxx | **CITE-CANT-FIND** | (1) Placeholder arXiv ID `2403.xxxxx`. (2) WebSearch for any Mondino-Sämann 2024 "synthetic Lorentzian Gromov-Hausdorff convergence" paper returns only the 2025 work at arXiv:2504.10380. (3) arXiv:2403.18619 (which other GeoVac papers wrongly use for this bibitem) is **"Enhanced OpenMP Algorithm to Compute All-Pairs Shortest Path on x86 Architectures"** by Calderón, Rucci, Chichizola — a completely unrelated paper. **There is no separate Mondino-Sämann 2024 synthetic-Lorentzian-GH paper**; `mondino_samann2024` is a phantom citation likely intended to refer to the same work as `mondino_samann2025` (arXiv:2504.10380). HIGH severity. **Resolution**: merge `mondino_samann2024` into `mondino_samann2025` (delete the 2024 bibitem and rename the 2 in-text `\cite` instances to use `mondino_samann2025` only). |
| `mondino_samann2025` | A. Mondino, C. Sämann, "Synthetic Lorentzian Gromov–Hausdorff convergence and pre-compactness" (2025), arXiv:2504.10380 | CITE-OK | Verified: arXiv:2504.10380, "Lorentzian Gromov-Hausdorff convergence and pre-compactness" (note: title omits "Synthetic" in arXiv version, but the framework IS synthetic — bibitem title is a reasonable paraphrase). Authors and 2025 year confirmed. |
| `reed_simon_iv` | M. Reed, B. Simon, *Methods of Modern Mathematical Physics Vol. IV: Analysis of Operators*, Academic Press 1978 | CITE-OK | Standard reference. The §XIII.16 sections cited for resolvent comparison + Bloch decomposition + Theorem XIII.86 are correctly located in Vol. IV. |
| `strohmaier2006` | "On noncommutative and pseudo-Riemannian geometry," J. Geom. Phys. 56 (2006), 175–195 | CITE-OK | Standard reference; matches arXiv:math-ph/0110001 / J. Geom. Phys. publication. |
| `paper24` to `paper48` (internal GeoVac) | Various internal Loutey preprints | CITE-OK (internal; outside Pass B scope) | All present in the corpus. |

### Problems found

1. **`bizi_brouder_besnard2018`** — CITE-WRONG-METADATA. Journal claimed "Class. Quant. Grav. 35 (2018) 175004"; actual is **J. Math. Phys. 59 (2018) 062303**. Same error in Papers 38/39/40/42/43/44/45/46 — corpus-wide. **Recommended fix**: update bibitem journal/volume/page in Paper 47 (and the eight other papers — but that's out of scope for this audit).

2. **`hekkelman_mcdonald2024`** — CITE-WRONG-METADATA (title) AND CITE-DOESNT-SUPPORT (in §6.2). Actual paper at arXiv:2412.00628 is "A noncommutative integral on spectrally truncated spectral triples, and a link with quantum ergodicity," NOT "The semi-classical limit of Connes spectral triples on a torus." The §6.2 usage claims it establishes "Tauberian estimates on flat tori with an extension to $\R^d$ via covering map" — but 2412.00628 is on truncated spectral triples + quantum ergodicity (Szegő-type). Whoever drafted this likely confused a Hekkelman-led 2024 work with a different paper (possibly Hekkelman's earlier asymptotic-Connes-integral work, or McDonald's earlier non-commutative-integration work; needs FIND to disambiguate). **Recommended fix**: update title in bibitem to "A noncommutative integral on spectrally truncated spectral triples, and a link with quantum ergodicity" AND rewrite §6.2 to either (a) drop the Tauberian-on-tori framing and instead cite 2412.00628 for its actual content, or (b) FIND the actual Hekkelman / McDonald paper that does the Tauberian-on-tori work.

3. **`hekkelman_mcdonald2024b`** — CITE-CANT-FIND. Placeholder arXiv ID `2411.xxxxx`; no separate Hekkelman–McDonald 2024 paper on "Noncommutative integration on non-compact spectral triples" surfaces in WebSearch. **Recommended fix**: either FIND a real paper that matches this description, or delete the bibitem and remove the single in-text `\cite{hekkelman_mcdonald2024b}` reference in §6.2 (line 1147). If 2412.00628 actually contains both the "semi-classical" and "non-compact integration" content (it doesn't — it's about quantum ergodicity), this entire reference may collapse into a confused composite. Likely both `hekkelman_mcdonald2024` and `hekkelman_mcdonald2024b` need to be re-FOUND.

4. **`latremoliere_metric_st_2017`** — CITE-WRONG-METADATA. Journal data wrong: should be **Adv. Math. 404 (2022), 108393** (not 411 (2023), 108790). Title also has "dual" inserted that is not in arXiv:1811.10843. Same wrong metadata in Papers 38/39/40/42/43/44/45/46.

5. **`mondino_samann2024`** — CITE-CANT-FIND. No separate 2024 Mondino-Sämann synthetic-Lorentzian-GH paper exists in the indexed literature. **Resolution recommendation** (per the dispatcher's question): the cleanest fix is **delete the `mondino_samann2024` bibitem and rename the two in-text `\cite` instances** (lines 345 and 1125) to use `mondino_samann2025` only. Looking at the in-text usage:
   - Line 345 (`\cite{mondino_samann2024,mondino_samann2025}`): the discussion is generic "synthetic Lorentzian GH program" — single citation to 2025 suffices.
   - Line 1125 (`\cite{mondino_samann2024,mondino_samann2025}` in §6.1): again generic "Mondino–Sämann ... develop a synthetic Lorentzian GH convergence" — single citation to 2025 suffices.
   - In other corpus papers, the parallel two-citation pattern is similar (Paper 45 lines 374, 1463; Paper 46 line 381 etc.). All such corpus-wide "merge mondino_samann2024 into mondino_samann2025" fixes are mechanical.

6. **`franco_eckstein2014`** — Minor CITE-WRONG-METADATA: bibitem says "(2013), 135007 (and follow-ups)" but the key is `franco_eckstein2014` and the journal year is 2013. Cosmetic mismatch; key should probably be `franco_eckstein2013`. LOW severity.

### Priority / novelty claims

| claim (verbatim) | location | searched | prior art? | recommendation |
|---|---|---|---|---|
| "the present paper is to our knowledge the first to make the explicit identification of the non-compact Lorentzian Dirac as the norm-resolvent limit of a sequence of truncated Krein spectral triples" | §6.2 | "Lorentzian Krein spectral triple norm-resolvent convergence non-compact"; Hekkelman/McDonald 2412.00628; Mondino-Sämann 2504.10380; Latrémolière 2512.03573 | None found doing exactly this (Krein + Lorentzian + norm-resolvent + non-compact). Adjacent work exists (Mondino–Sämann synthetic-Lorentzian GH; Hekkelman–McDonald Riemannian; Latrémolière Riemannian pointed-proper; Strohmaier 2006 pseudo-Riemannian framework). But absence of prior hit in clean search cannot establish priority. | Keep "to our knowledge the first" — the hedging is correct and unimpeachable by a search. |
| "the first published Lorentzian propinquity convergence theorem on truncated Krein spectral triples" (carried over from Paper 45 §1.1, not directly in Paper 47 but referenced via Paper 45 chain) | n/a in Paper 47 | (out of scope) | n/a | n/a |

---

## Combined severity table

| Finding | Pass | Verdict | Severity |
|---|---|---|---|
| `bizi_brouder_besnard2018` wrong journal/volume/page (corpus-wide) | B | CITE-WRONG-METADATA | MEDIUM |
| `hekkelman_mcdonald2024` wrong title AND citation doesn't support §6.2 claim | B | CITE-WRONG-METADATA + CITE-DOESNT-SUPPORT | **HIGH** |
| `hekkelman_mcdonald2024b` placeholder arXiv ID `2411.xxxxx`; no separate paper found | B | CITE-CANT-FIND | **HIGH** |
| `latremoliere_metric_st_2017` wrong journal/volume/year (corpus-wide) | B | CITE-WRONG-METADATA | MEDIUM |
| `mondino_samann2024` phantom citation; no 2024 Mondino-Sämann synthetic-Lorentzian-GH paper exists | B | CITE-CANT-FIND | **HIGH** |
| `franco_eckstein2014` minor year mismatch with bibkey | B | CITE-WRONG-METADATA | LOW |
| §8 Q2 "CLOSED" overstates body (only sketched, not theorem-grade) | A | C (overstatement) | MEDIUM |
| Outer rate constant cosmetic gap: announced $C(z,\xi)\le c_0\|\xi\|(1+|z|)^{-1}|\Im z|^{-2}$, body proves only $c_0\|\xi\|/|\Im z|^2$ at compact-support limit | A | C (minor) | LOW |
| Most GEOVAC-ONLY chains (A1, A2, A7, A9) — risk: if Paper 38 / Paper 46 has a bug, Paper 47 inherits it; but each chain link verified internally | A | B (internal-consistent) | (note, not a finding) |
| Krein-lift of Latrémolière 2512.03573 (Paper 48 §3) is the load-bearing new input to Thm 7.3; Paper 48 itself not yet externally vetted | A7 | D (unverifiable here) | MEDIUM |

**Totals:**

Pass A verdicts: A=4, B=6, C=2 (overstatement Q2 + cosmetic outer-rate gap), D=1, E=0.

Pass B verdicts: CITE-OK = 13, CITE-WRONG-METADATA = 3 (`bizi_brouder_besnard2018`, `latremoliere_metric_st_2017`, `franco_eckstein2014` minor) plus 1 mixed (`hekkelman_mcdonald2024` title+support), CITE-DOESNT-SUPPORT = 1 (`hekkelman_mcdonald2024` again, §6.2 framing), CITE-MISATTRIBUTED = 0, CITE-CANT-FIND = 2 (`hekkelman_mcdonald2024b`, `mondino_samann2024`). Effective Pass B problem count = 5 distinct citation issues. **`mondino_samann2024` in-text usage count = 2** (lines 345 and 1125 — both `\cite{mondino_samann2024,mondino_samann2025}`, neither sole).

Severity totals: HIGH = 3, MEDIUM = 4, LOW = 2.

---

## Broadcast readiness: **YELLOW**

The math content (Pass A) is in good shape: the outer arrow (§4 Thm 4.2) is a clean, externally-grounded compact-to-non-compact-domain norm-resolvent argument; the three-carrier identification (§6 Thm 6.1) is standard; the strategic reframing (§1.1 AF-vs-connected-spectrum) is mathematically correct; the inner arrow (§3 Thm 3.1) transcribes Paper 45/46 machinery internal-consistently. The main paper-level overstatement is §8 Q2 "CLOSED (flip-suppression)" which should be softened to "closure mechanism sketched"; the §7 Thm 7.3 G2-metric closure is theorem-grade on the natural substrate but chains through three GeoVac-internal links (Paper 46 Lemma 3.2 + Latrémolière 2512.03573 axioms + Paper 48 §3 Krein-lift), with the Paper 48 Krein-lift being the most-exposed unvetted ingredient.

The citation hygiene (Pass B) is the main problem and is why this is YELLOW not GREEN. Five distinct citation issues, three of them HIGH severity, including two CITE-CANT-FIND (`mondino_samann2024`, `hekkelman_mcdonald2024b`) and one CITE-DOESNT-SUPPORT (`hekkelman_mcdonald2024` doesn't say what §6.2 uses it to claim). The phantom `mondino_samann2024` is a corpus-wide pattern that needs systematic cleanup (Papers 38/39/40/42/43/44/45/46 also have it). The `latremoliere_metric_st_2017` and `bizi_brouder_besnard2018` wrong-metadata are likewise corpus-wide. None of these issues bears on the math; all are fixable in an errata sprint without paper-level rewriting.

Pre-broadcast checklist:
1. Resolve `mondino_samann2024` — recommend delete bibitem + rename 2 in-text `\cite` instances to `mondino_samann2025` only.
2. Resolve `hekkelman_mcdonald2024` — re-FIND the actual paper(s) intended for §6.2's Tauberian-on-tori claim. If the real intended paper is 2412.00628, rewrite §6.2 to match that paper's actual content (quantum ergodicity / Szegő-type, NOT Tauberian-on-tori).
3. Resolve `hekkelman_mcdonald2024b` — likely merge into 2024 entry or delete.
4. Fix `bizi_brouder_besnard2018` journal/volume/page.
5. Fix `latremoliere_metric_st_2017` journal/volume/year/title.
6. Soften §8 Q2 "CLOSED" → "closure mechanism sketched" or similar.

---

## What I could NOT verify (hand to a human expert)

1. **Whether 2412.00628 (Hekkelman–McDonald) actually contains content relevant to Tauberian-on-tori or to non-compact spectral triples.** Web search and abstract review say it does NOT (it's about quantum ergodicity on truncated spectral triples). But there may be a separate Hekkelman-led or McDonald-led 2024 work on the actual claimed topic. A domain expert (or a careful arXiv author-page sweep) is needed to FIND the right paper(s).

2. **Whether `farsi_latremoliere2024` (J. Funct. Anal. 286 (2024) 110293) checks out at full metadata level.** I accepted it as plausible based on cross-corpus consistency but did not run an arXiv landing-page WebFetch.

3. **Whether the Paper 48 §3 Krein-lift of Latrémolière's 2512.03573 axioms is itself correctly executed** at the level needed to make Paper 47 Thm 7.3 watertight. This is the load-bearing GeoVac-internal input; Paper 48 is itself an active sprint product not yet externally vetted.

4. **Whether there is genuine prior art for the explicit identification of the non-compact Lorentzian Dirac as a norm-resolvent limit of truncated Krein spectral triples.** A clean search returned no direct hit, but absence is not proof of absence. Domain-expert confirmation required.

5. **Whether the inner arrow's $T$-independence at L2 (cb-norm 2/(\nmax+1) on $\SU(2)\times U(1)_T$) is in fact $T$-uniform** as claimed. The Bożejko–Fendler central-multiplier theorem on amenable compact groups gives the bound; the product $\SU(2)\times U(1)_T$ is amenable for every $T$ (both factors compact), but whether the constant 2/(\nmax+1) is genuinely $T$-uniform (vs. picking up a small $T$-correction) deserves direct verification. Paper 45 §5 is the load-bearing source; I did not re-derive.
