# Phase 2 Connection-Finding — deep structural arcs
Date: 2026-06-03

## Methodology

Web search across arXiv, Google Scholar, ResearchGate, and Springer/Elsevier journal indexes.
Date window emphasis: 2020–2026, with 2024–2026 prioritized. Cross-checked author homepages
(van Suijlekom, Sämann, Connes) and followed up structural matches via WebFetch on arXiv
abstracts. Each theme was searched with 2–4 queries combining theme-specific keywords with
the names of likely interlocutors (Connes, Chamseddine, Marcolli, van Suijlekom, Hekkelman,
Latrémolière, Perez-Sanchez, Mondino, Sämann, Bochniak, Sitarz, Eckstein). No author was
contacted; this is a literature-only sweep.

## Per-theme findings

### Theme 1 — Master Mellin engine (Paper 18 §III.7)

- **Closest active programs** (2024–2026):
  1. **Marcolli & Pierpaoli, "Periods and motives in the spectral action of Robertson-Walker spacetimes"** (arXiv:1611.01815, building on Fathizadeh–Marcolli, arXiv:1407.5972 "Rationality of spectral action for Robertson-Walker metrics"). Heat-kernel Seeley–DeWitt coefficients of the Dirac-Laplacian on Robertson–Walker metrics are explicitly identified as **periods of mixed Tate motives** — algebraic differential forms integrated over semi-algebraic complements of unions of quadrics and hyperplanes. Continuing arc through Fathizadeh–Ghorbanpour–Marcolli "Bell polynomials and Brownian bridge in spectral gravity models on multifractal Robertson–Walker cosmologies" (arXiv:1811.02972), and Marcolli "Spectral action gravity and cosmological models" (CR Physique 2017 review).
  2. **Estrada–Fulling-style Mellin-transform expansion**: implicit in the 2024 *Complexity* paper on heat kernels for networks with long-range interactions (Kalala Mutombo, https://doi.org/10.1155/2024/6745905). Uses normalized Laplace and Mellin transforms of the k-path Laplacian — structurally adjacent to Paper 18 §III.7's M1/M2/M3 partition for graph operators.
  3. **Rejzner, "Renormalization and periods in perturbative AQFT"** (arXiv:1603.02748, in Springer 2020 *Progress in Mathematics* vol. 335). pAQFT β-function coefficients identified as Kontsevich–Zagier periods. Direct conceptual sibling of the master Mellin engine: identifies *which periods* appear as a function of the perturbative order, analogous to GeoVac's claim that k ∈ {0,1,2} indexes mechanism.
  4. **Connes–Marcolli, "Renormalization and motivic Galois theory"** (arXiv:math/0409306; book *Noncommutative Geometry, Quantum Fields and Motives* AMS Colloq. Pub. 55, 2008). Renormalization group as one-parameter subgroup of motivic Galois group. The most ambitious published "where do the transcendentals come from" program.
  5. **Bochniak–Dabrowski–Sitarz–Zalecki, "On geometric spectral functionals"** (May 2025 preprint). Investigates spectral functionals associated with Dirac and Laplace-type operators on manifolds — direct relevance to which functionals can produce M1 vs M2 content.

- **Connections worth building**:
  - **Mixed Tate motives ↔ M2**: Fathizadeh–Marcolli's "Seeley–DeWitt coefficients are mixed Tate periods" is the natural mathematical home for the M2 mechanism. Paper 18 §III.7 currently identifies M2 as "Seeley–DeWitt via Mellin of Tr(e^{-tD²}) at k=2"; the Fathizadeh–Marcolli result tells us *which* arithmetic class M2 lives in. **Could be cited as an independent characterization of M2's natural arithmetic home.**
  - **Periods ↔ master Mellin engine**: Rejzner's pAQFT periods program is the closest published analog of GeoVac's claim that *only* M1/M2/M3 transcendentals appear. Different setting (Minkowski pAQFT vs S³ spectral triple), but the structural question — "given a perturbative observable, which transcendentals can appear?" — is the same. **Could be cited as the closest existing operationalization of the question Paper 18 §III.7 answers.**
  - **Spectral action on fractal/multifractal substrates**: Fathizadeh–Marcolli's 2018 multifractal Robertson–Walker work introduces log-periodic corrections from poles of the fractal-string zeta. GeoVac's discrete spectral action on S³ is structurally analogous (compact substrate, discrete spectrum, similar Mellin expansion). Worth noting in Paper 18 as a finite-cutoff analog of their fractal-cutoff calculation.

- **Specific citations to add to Paper 18**:
  - Fathizadeh–Marcolli arXiv:1611.01815 (mixed Tate motives for SD coefficients) — should be added to §III.7 next to the M2 definition.
  - Rejzner arXiv:1603.02748 (Kontsevich–Zagier periods in pAQFT) — should be added to §III.7 master-Mellin remark as the closest published "which transcendentals appear" classification.
  - Connes–Marcolli arXiv:math/0409306 (motivic Galois theory of renormalization) — should be added as the broader framework into which the case-exhaustion theorem (Paper 32 §VIII) fits as a finite-cutoff slice.
  - Bochniak–Dabrowski–Sitarz–Zalecki "On geometric spectral functionals" (May 2025) — recent enough to be worth a single-line reference if the corpus is being audited.

- **Open question for the PI**: The Fathizadeh–Marcolli result says SD coefficients are mixed Tate periods on Robertson–Walker. **Are GeoVac's discrete SD coefficients on S³ also mixed Tate?** If yes, this is a non-trivial structural alignment with their program and a natural collaboration. If no, the discrete sector is *more* restricted than the continuum sector — a potentially novel claim worth highlighting.

### Theme 2 — Projection taxonomy (Papers 34, 35)

- **Closest active programs**:
  1. **Categorical Operational Physics** (Gogioso & Genovese, arXiv:1902.00343, and the broader Abramsky–Coecke CQM lineage). Most direct philosophical analog of "two-layer + projections" framing: a categorical separation of process structure (Layer 1 analog) from operational realization (Layer 2 analog). No published version that targets *transcendentals* as the discriminator the way Paper 35 does.
  2. **Döring–Isham topos-theoretic reformulation of QM** (broader program at https://golem.ph.utexas.edu/category/2008/04/algebraic_quantum_mechanics_an.html, Heunen–Landsman–Spitters lineage). Provides a topos-internal logic where physical content is *constructed by projection from a context*. Paper 34's "projections inject physics" framing maps cleanly onto this.
  3. **Bisognano–Wichmann literature 2024+**: Camassa et al. "The Bisognano-Wichmann property for non-unitary Wightman conformal field theories" (arXiv:2506.10625); Crisand "A geometric perspective on Algebraic Quantum Field Theory" (arXiv:2412.20410). These provide *theorem-grade* "where physics enters via signature/modular structure" classifications — direct analog of Paper 35's "π enters at temporal compactification."
  4. **Mondino–Sämann synthetic Lorentzian program** (Cavalletti–Mondino arXiv:2004.08934 in Cambridge J. Math. 2024; Manzano–Mosani–Sämann–Zoghlami arXiv:2512.05842 Dec 2025 on conformal Lorentzian length spaces; Sämann START project 2024–2029). A "where physics enters" program for Lorentzian curvature: separates the metric structure (Layer 1 analog) from the timelike Ricci-curvature realization (Layer 2 analog) via optimal-transport projections.

- **Connections worth building**:
  - **Bisognano–Wichmann ↔ Paper 35's π-from-temporal-compactification prediction**: BW theorems classify *where* modular structure forces π and 2π factors. Camassa et al. 2026 generalizes this to non-unitary CFTs. Paper 35's falsifiable prediction ("a GeoVac observable contains π iff its evaluation includes continuous integration over a temporal/spectral parameter") is a discrete-substrate analog of BW. **Could be cited as the natural framework into which Paper 35's prediction extends.**
  - **Mondino–Sämann ↔ Paper 34's projection axes**: Their synthetic Lorentzian program already has a Sämann START project funded 2024–2029 (5-year horizon). They are the natural mathematical partners for the Lorentzian arc (Papers 42–49). Paper 47–48 already cites them; the connection is established. Worth deepening for any Phase 5 / W-axis work.
  - **Topos / categorical QM ↔ Layer-1/Layer-2 split**: Döring–Isham's "physics is what's measurable in a context" is the closest philosophical anchor for Paper 35's two-layer framing. Paper 34 does not need to be reframed in topos language — but a single paragraph relating "projection" to "context-dependent valuation in a presheaf topos" would defuse the "no published framework does this" concern and give a citation home.

- **Specific citations to add to Paper 34/35**:
  - Camassa et al. arXiv:2506.10625 (BW for non-unitary Wightman CFTs) — single citation in Paper 35 §VI (falsifiable prediction) as the closest theorem-grade analog.
  - Crisand arXiv:2412.20410 (geometric perspective on AQFT) — citation in Paper 34 as a recent (Dec 2024) framework adjacent to the projection-taxonomy spirit.
  - Cavalletti–Mondino Cambridge J. Math. 2024 — already cited in Papers 47/48; worth a single line in Paper 35 as the Lorentzian-curvature analog of Paper 35's temporal-compactification mechanism.

- **Open question for the PI**: Paper 35's "π appears iff continuous temporal integration" is theorem-grade in GeoVac but is essentially *the* content of Bisognano–Wichmann translated to a discrete substrate. **Is Paper 35's prediction equivalent to BW for the GeoVac wedge KMS state?** If yes, then Paper 35 is providing a discrete-substrate proof of a known continuum theorem — substantively easier to publish and cite.

### Theme 3 — Spectral triple construction (Paper 32)

- **Closest active programs**:
  1. **Perez-Sanchez "Bratteli networks and the Spectral Action on quivers"** (arXiv:2401.03705, Jan 2024) and **"Comment on 'Gauge networks in noncommutative geometry'"** (arXiv:2508.17338, Aug 2025). Both critical 2024–2025 follow-ups to the Marcolli–vS 2014 paper. The Comment establishes that the continuum limit of Marcolli–vS gauge networks is **pure Yang–Mills (no Higgs)** — a refinement directly relevant to Paper 32 §VIII.B's gauge-network reading and Sprint H1's Higgs verdict. The Bratteli-networks paper introduces **prespectral triples** (discrete topological NC spaces hosting quiver representations) as the natural finite-dimensional substrate. GeoVac's S³ truncated triple sits in this lineage.
  2. **Hekkelman & McDonald, "A noncommutative integral on spectrally truncated spectral triples"** (arXiv:2412.00628, Dec 2024; J. Funct. Anal. accepted July 2025). Provides a quantum-ergodicity-compatible noncommutative integral for the Connes–vS truncated paradigm. Direct sibling to Paper 32's truncated-triple constructions; potentially load-bearing for any GeoVac "trace formula on truncated triple" claim.
  3. **Hekkelman PhD thesis, "Trace Formulas in Noncommutative Geometry"** (arXiv:2506.21950, June 2025). Systematizes heat-trace expansions, Dixmier traces, and density-of-states formulas across NCG. Bridges solid-state physics applications with operator-algebraic NCG.
  4. **van Suijlekom, "A generalization of K-theory to operator systems"** (arXiv:2409.02773, Sept 2024). Defines K₀ on operator systems, reducing to standard K-theory on C*-algebras. Direct relevance to Paper 32's operator-system-level constructions and Paper 44/45/46's Krein extensions.
  5. **Farsi–Latrémolière "Collapse in Noncommutative Geometry and Spectral Continuity"** (arXiv:2404.00240, Apr 2024). New necessary-and-sufficient condition for inductive-sequence convergence in propinquity. Direct relevance to Papers 38, 45, 46, 47, 48, 49 (the math.OA arc).
  6. **Farsi–Latrémolière–Packer "Convergence of inductive sequences of spectral triples for the spectral propinquity"** (Adv. Math. 437, 2024). Already cited by GeoVac Papers 47/48 in different forms.
  7. **Connes–Consani–Moscovici "Zeta zeros and prolate wave operators"** (arXiv:2310.18423, Ann. Funct. Anal. 2024) plus Connes & Consani Jan 2025 geometric class-field-theory extension. Different problem (Hilbert–Pólya conjecture) but same operator-algebraic toolkit; relevant for any future GeoVac RH-track resumption.

- **Connections worth building**:
  - **Perez-Sanchez Comment 2025 ↔ Paper 32 §VIII.B + Sprint H1**: The Comment's "Yang–Mills *without* Higgs in the continuum limit of Marcolli–vS gauge networks" is the *exact* structural result GeoVac Sprint H1 (POSITIVE-THIN) produced for the AC extension. This is currently noted in CLAUDE.md WH1 entry, but **the cross-reference is not in Paper 32 §VIII.B itself**. Should be added — this is the strongest external corroboration of the Sprint H1 verdict.
  - **Hekkelman–McDonald 2024 noncommutative integral ↔ Papers 38/45**: Their noncommutative integral on truncated triples is the natural Layer-2 observable companion to the GeoVac propinquity convergence. If GeoVac wants to make any "Connes-integral at finite cutoff" claim downstream, this is the published machinery to cite.
  - **Farsi–Latrémolière 2024 collapse paper ↔ Papers 38/40 GH-convergence**: New convergence criterion in propinquity. Should be checked against the five-lemma sequence — could provide an alternative proof path or a sharper bound at qualitative-rate level.
  - **van Suijlekom 2024 operator-system K-theory ↔ Paper 44 operator-system extension**: Direct sibling — same year, same author lineage, same operator-system setting. Should be cited in Paper 44's introduction as the K-theoretic complement to the Krein-substrate propinquity work.

- **Specific citations to add to Paper 32**:
  - Perez-Sanchez arXiv:2508.17338 (the *Comment*) — in §VIII.B alongside Marcolli–vS as the corrected continuum limit. Already in CLAUDE.md, missing from Paper 32 itself.
  - Perez-Sanchez arXiv:2401.03705 (Bratteli networks) — single citation in Paper 32 §III as a 2024 sibling construction.
  - Hekkelman–McDonald arXiv:2412.00628 — in Paper 32 §IV (Connes axiom audit) as the truncated-triple integration framework.
  - van Suijlekom arXiv:2409.02773 — in Paper 32 §VIII or Paper 44 §1 as the operator-system K-theory complement.

- **Open question for the PI**: Perez-Sanchez's 2024+2025 program is the most direct conversation partner for GeoVac's spectral-triple framing. Both arcs ended in similar "Yang–Mills no Higgs" verdicts. **A 1-page note co-citing Perez-Sanchez 2025 + Sprint H1 as independent confirmations of the Marcolli–vS-without-Higgs reading would be the natural cross-validation publication.** Email-scale outreach to Perez-Sanchez may be warranted.

### Theme 4 — Case-exhaustion theorem (Paper 32 §VIII)

- **Closest active programs**:
  1. **Marcolli–Fathizadeh mixed-Tate-motives program for SD coefficients** (same arXiv:1611.01815 as Theme 1). Provides a **theorem-grade arithmetic classification** of which periods can appear in heat-kernel expansions on Robertson–Walker spacetimes. Structurally identical question to the case-exhaustion theorem at the SD-coefficient level.
  2. **Brouder–Dang–Hélein wavefront-set program** (referenced via Rejzner pAQFT). Operationalist version of "which singularities can appear in perturbative observables."
  3. **Cacic "Moduli spaces of Dirac operators for finite spectral triples"** (arXiv:0902.2068). Foundational classification of which finite spectral triples can exist; structurally the discrete-substrate analog of the case-exhaustion question. No 2024–2026 explicit follow-up surfaced, but a literature audit on Cacic's recent work would close any gap.
  4. **Krajewski-diagram classification program** (Krajewski 1997 hep-th/9701081; Paschke–Sitarz 1996). Combinatorial enumeration of finite real spectral triples with prescribed KO-dimension. Discrete-substrate sibling of Paper 32 §VIII's exhaustion claim.
  5. **Camassa et al. BW property non-unitary CFTs** (arXiv:2506.10625, 2026) — theorem-grade classification of *when* the BW property holds. Conceptually adjacent: "which projections preserve which structural properties" is the kind of theorem-grade content the case-exhaustion theorem is.

- **Connections worth building**:
  - **Mixed-Tate-motives classification ↔ M2 case of the theorem**: Fathizadeh–Marcolli's result that SD coefficients live in the mixed-Tate-period sub-ring is *strictly the M2 case* of the case-exhaustion theorem. If GeoVac's discrete SD coefficients also live in the mixed-Tate ring, then Paper 32 §VIII's case-exhaustion theorem at k=2 follows from theirs. **This is a candidate published-proof transport for the M2 case.**
  - **Krajewski / Cacic finite-triple classification ↔ axiom-grade base**: The classification literature already provides theorem-grade enumeration of which finite spectral triples can host which gauge structures. Paper 32 §VIII implicitly relies on this; explicit citation would harden the foundations.
  - **BW + Camassa et al. ↔ classifications of projection-preservation**: The case-exhaustion theorem says "every π in any finite projection chain engages M1/M2/M3." This is structurally analogous to "every BW property in any AQFT net engages signature-modular-Wick-rotation." Worth flagging the analogy in Paper 32 §VIII Remark.

- **Specific citations to add to Paper 32 §VIII**:
  - Fathizadeh–Marcolli arXiv:1611.01815 — in the M2 case statement as the published transport theorem for mixed-Tate periods.
  - Krajewski hep-th/9701081 and Cacic arXiv:0902.2068 — in §VIII's introduction as the existing finite-spectral-triple classification literature.
  - Connes–Marcolli motivic Galois arXiv:math/0409306 — in §VIII closing remark as the broader framework into which the case-exhaustion theorem extends.

- **Open question for the PI**: **Are GeoVac's discrete SD coefficients on S³ explicitly mixed-Tate periods?** If yes, the M2 case of the case-exhaustion theorem inherits from Fathizadeh–Marcolli and becomes a corollary, not an independent claim. If no, the discrete sector is more restricted, which is publishable structural content in its own right. A short ~1-week sprint computing 3–4 SD coefficients in the mixed-Tate-period testable form would answer this definitively.

## Cross-cutting observations

1. **Hekkelman/van Suijlekom/Latrémolière is the natural collaboration cluster across themes 3 and 4.** Three 2024 papers (Hekkelman–McDonald noncommutative integral, van Suijlekom operator-system K-theory, Farsi–Latrémolière collapse) and one 2024 Adv. Math. paper (Farsi–Latrémolière–Packer inductive convergence) are all directly relevant to GeoVac's math.OA arc. The cluster is small (~4–6 active researchers), publishes regularly together, and is the natural peer group for the GeoVac papers 38–49. **Recommend at least one direct citation cycle.**

2. **Marcolli/Fathizadeh program is the natural conversation partner across themes 1 and 4.** Their published "SD coefficients are mixed-Tate periods on RW spacetimes" is the closest theorem-grade analog of Paper 18 §III.7 M2 + Paper 32 §VIII case-exhaustion theorem at k=2. They have a clear publication trajectory through 2018+ and continue active through Marcolli's PhD students. **Recommend a 1-week sprint to test whether GeoVac discrete SD coefficients on S³ live in the mixed-Tate ring.**

3. **Perez-Sanchez 2024–2025 is the natural conversation partner for the spectral-triple gauge-network arc (Paper 32 §VIII.B).** His Aug 2025 Comment establishes Yang–Mills-without-Higgs as the corrected continuum limit of Marcolli–vS — *the same verdict GeoVac Sprint H1 reached* on the discrete substrate. This convergence is currently invisible because Paper 32 doesn't cite him. **Single highest-value citation addition across the entire sweep.**

4. **Bisognano–Wichmann and Mondino–Sämann are the natural framing partners for Papers 34 and 35.** Both provide theorem-grade "where physics enters via signature/modular/causal structure" classifications. Paper 35's falsifiable π-prediction is a discrete-substrate BW statement; explicit framing would clarify the contribution.

5. **The "categorical / operationalist projection framework" question (Theme 2) returned weaker matches than expected.** No published 2024–2026 program explicitly does "two-layer skeleton-vs-calibration" the way Paper 34/35 do. This is mildly positive: Paper 35's two-layer framing appears genuinely novel as an *operational classification scheme*. (BW-style classifications exist, but as theorems within specific frameworks, not as the framework itself.)

## Highest-value follow-ups

1. **Add Perez-Sanchez arXiv:2508.17338 to Paper 32 §VIII.B.** Single load-bearing citation. The Yang–Mills-without-Higgs continuum-limit verdict published Aug 2025 corroborates Sprint H1's POSITIVE-THIN verdict on the discrete substrate. ~30 min of edit time. Highest cite-to-cost ratio in the sweep.

2. **1-week sprint: test whether GeoVac discrete SD coefficients on S³ are mixed-Tate periods.** Concrete falsifiable computation. If POSITIVE: Paper 32 §VIII M2 case inherits from Fathizadeh–Marcolli; Paper 18 §III.7 M2 mechanism gets a published arithmetic home; opens collaboration channel with Marcolli's group. If NEGATIVE: discrete sector is more arithmetically restricted than continuum — publishable structural content in its own right.

3. **Add Hekkelman–McDonald arXiv:2412.00628, van Suijlekom arXiv:2409.02773, Farsi–Latrémolière arXiv:2404.00240 to Papers 32/38/44/45/46.** Three citations across the math.OA arc, each ~10 min of edit time. Brings GeoVac into citation alignment with the natural peer cluster. Plausibly initiates reciprocal citation in their next round of papers.

## Limitations

- Search was English-language only; missed any French/German/Russian work in the area.
- WebSearch returns up to 10 results per query; this is enough for connection-finding but not for an exhaustive precedence audit. Any claim of "first in literature" should be checked separately.
- arXiv abstract-only WebFetch can miss technical content; the recommended citations were verified to abstract-grade only, not full-paper-grade.
- The "categorical / operationalist projection" theme (Theme 2) returned weaker matches than expected; this could be a genuine gap, or it could reflect that the framing question is unusual enough that searches don't cleanly hit it. A targeted manual search of Coecke / Heunen / Spitters recent papers may close any remaining gap.
- The case-exhaustion theorem (Theme 4) shares structure with both the mixed-Tate-period program and the Krajewski-diagram classification, but no published source uses the exact "case exhaustion across operator order × bundle type" framing GeoVac uses. The framing itself appears novel; the *content* sits inside known programs.
