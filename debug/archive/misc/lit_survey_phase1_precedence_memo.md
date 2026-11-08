# Phase 1 Precedence Audit — math.OA standalone novelty (Papers 45–50)
Date: 2026-06-03
Auditor: PM (Claude Opus 4.7) under dispatch from PI

## Methodology

- WebSearch (Anthropic web tool) sweeps over arXiv math.OA / math-ph / math.MG / math.DG / hep-th, 2024–2026 windows, with targeted 2018–2023 baselines for the Latrémolière / Connes–van Suijlekom / Mondino–Sämann / Eckstein–Franco / van den Dungen lineages.
- WebFetch on individual arXiv abstract pages for the highest-relevance hits (2512.03573 Latrémolière hypertopology; 2504.10380 Mondino–Sämann Lorentzian GH; 2506.10852 Braun–Sämann reconstruction theorem; 2412.00628 Hekkelman–McDonald; 2502.18105 / 2512.15450 Nieuviarts).
- Cross-checked existing CLAUDE.md scope-claim wording against current arXiv listings, with explicit attention to (i) the 2025 hypertopology and (ii) the 2024–2025 synthetic-Lorentzian-GH thread that has expanded since the May-2026 audits in CLAUDE.md.
- **Limitations**: web-tool access only; no direct arXiv-listing scrape; no DOI-resolution for journal-only references; cannot read full PDFs, only abstracts + introductions + section excerpts. A more exhaustive audit would (i) read the May–June 2026 listings on math.OA + math-ph day-by-day, (ii) check zbMATH / MathSciNet for journal-only items, (iii) ask Eva-Maria Hekkelman / Frédéric Latrémolière directly. Recommend re-running the audit one week before each arXiv submission.

---

## Per-paper findings

### Paper 45 — K⁺-restricted weak-form Lorentzian propinquity (Sprint L3b-2, arXiv-ready)

- **Closest existing work**:
  1. **van den Dungen 2016** ("Krein spectral triples and the fermionic action", *Math. Phys. Anal. Geom.* 19; arXiv:1505.01939) — defines Krein spectral triples; **no convergence content**.
  2. **Mondino–Sämann 2025** ("Lorentzian Gromov–Hausdorff convergence and pre-compactness", arXiv:2504.10380; v4 Dec 2025) — synthetic Lorentzian GH convergence + pre-compactness for ε-nets of causal diamonds, by Mondino & Sämann themselves. Pure differential geometry; **no operator algebras / spectral triples in scope**.
  3. **Braun–Sämann 2025** (arXiv:2506.10852, June 2025) — Gromov reconstruction theorem + measured GH convergence in Lorentzian geometry. Again synthetic, **no operator algebras**.
  4. **Nieuviarts 2025** (arXiv:2502.18105 v3 May 2025; proceeding arXiv:2512.15450 Dec 2025) — emergence of Lorentz / pseudo-Riemannian structure from twisted spectral triples on almost-commutative manifolds. Structural Wick-rotation morphism, **no propinquity convergence**.
- **Novelty verdict**: **NOVEL**. The "first published Lorentzian propinquity convergence theorem on truncated Krein spectral triples" wording in CLAUDE.md still stands as of June 3, 2026. The synthetic-Lorentzian-GH program is active (Mondino–Sämann themselves are publishing) but is categorically different (no Krein, no operator algebras). The almost-commutative-twisted-Wick-rotation program (Nieuviarts) is the closest spectral-triple-side concurrent thread but addresses *emergence*, not *convergence*.
- **Citations to add**: Mondino–Sämann arXiv:2504.10380 v4 (Dec 2025) should be cited in §1.4 alongside the existing Mondino–Sämann lineage as the most current synthetic-Lorentzian-GH convergence statement; explicitly note the operator-algebraic / synthetic categorical separation that Paper 48 then bridges. Braun–Sämann arXiv:2506.10852 (June 2025) is appropriate as a follow-up reference.
- **Notes**: Hekkelman–McDonald arXiv:2412.00628 (*J. Funct. Anal.* 2025) does NOT overlap — they treat noncommutative integration on truncated spectral triples for quantum ergodicity (Szegő limit), strictly Riemannian, no propinquity, no Krein. Confirm this in Paper 45 §1 with a citation-and-distinguish sentence.

### Paper 46 — Strong-form Lorentzian propinquity on natural chirality-doubled substrate

- **Closest existing work**: Same lineage as Paper 45. The strong-form-without-K⁺-restriction question with operator-norm Lipschitz seminorm $L_{op}(a) = \|[D_L, a]\|_{op}$ on a chirality-doubled scalar-multiplier substrate has no concurrent treatment in the literature surveyed. The closest Riemannian sibling is Latrémolière's own spectral propinquity (arXiv:1811.10843, 2018) and its 2024–2025 extensions (Farsi–Latrémolière 2504.11715, 2025); none of these address the Krein / Lorentzian case.
- **Novelty verdict**: **NOVEL**. The closed-form $C_3^{op}(n_{max}) = \sqrt{1 - 1/n_{max}}$ is original; the $\Lambda^{strong} = \Lambda^{P45}$ "free upgrade" is original; the chirality-doubled gradient-norm absorption in Appendix B is original.
- **Citations to add**: cite Farsi–Latrémolière arXiv:2504.11715 (April 2025, *Continuity for the spectral propinquity of the Dirac operators associated with an analytic path of Riemannian metrics*) as a sibling Riemannian-side continuity theorem that uses a Lichnerowicz-style commutator analysis comparable to Lemma 3.2's temporal-Lipschitz-invisibility. The structural analogy is worth flagging.
- **Notes**: Paper 46's strong-form result on the *natural* substrate is what is novel. The deeper question on the enlarged substrate with chirality-flipping generators (where $\{J, M\} = 0$) is genuinely open in the literature, not just GeoVac-internal; this should be explicitly stated.

### Paper 47 — Norm-resolvent convergence to non-compact Lorentzian Krein spectral triple

- **Closest existing work**:
  1. **Reed–Simon Vol. IV §XIII.16** (1978) — classical norm-resolvent + Bloch-decomposition machinery; cited per CLAUDE.md.
  2. **Behrndt–Holzmann–Stelzer 2025** (arXiv:2404.07784; PMC12287183 PMC 2025) — norm-resolvent convergence for Dirac operators with confining Lorentz-scalar δ-shell interactions. Mathematically adjacent (Dirac + norm-resolvent + Lorentz scalar) but addresses a completely different physical setting (PDE on $\mathbb{R}^3$ with localized potentials approximating boundary conditions), with no spectral-triple structure and no Krein-space content.
  3. **Latrémolière arXiv:2512.03573** (Dec 3, 2025, *Quantum Gromov–Hausdorff hypertopology on pointed proper quantum metric spaces*) — published Riemannian non-compact extension of the propinquity framework, at the *metric* level (not norm-resolvent / spectral level). The compact-to-noncompact direction is opposite to Paper 47's, and the underlying technology (inframetric on pointed proper QMS) is different.
- **Novelty verdict**: **NOVEL**. The two-rate hybrid limit composition (propinquity at coupled cells + norm-resolvent of the Lorentzian Dirac with exponential tail $O(e^{-\|\Im z\| T/2})$) is original; the three-carrier identification (periodic / Dirichlet / non-compact converge to same physical Lorentzian Dirac) is original.
- **Citations to add**: Behrndt–Holzmann–Stelzer arXiv:2404.07784 should be cited in §1 and §6 as the closest PDE-side concurrent thread (norm-resolvent for Lorentz-scalar Dirac), with explicit distinction that GeoVac's setting is operator-algebraic / spectral-triple-based, not local-PDE. This matters for venue framing — math-ph reviewers may know that thread.
- **Notes**: The strategic-reframing acknowledgment in Paper 47 §1.1 (post-Phase-1B-C) that the AF inductive-limit framing is structurally wrong and that norm-resolvent is the right tool is the kind of honest scope statement that strengthens a precedence claim, not weakens it. Keep that paragraph.

### Paper 48 — Krein-pointed proper QMS + Wick-rotation bridge to Mondino–Sämann

- **Closest existing work**:
  1. **Latrémolière arXiv:2512.03573** (Dec 2025) — Riemannian pointed-proper-QMS hypertopology. Paper 48's Krein-side construction is the formal Krein-lift of this framework. CLAUDE.md already acknowledges Latrémolière's paper as the Riemannian-side input.
  2. **Mondino–Sämann arXiv:2504.10380** + **Braun–Sämann arXiv:2506.10852** (April–June 2025) — synthetic Lorentzian GH convergence; the target category for Paper 48's bridge.
  3. **Cavalletti–Mondino 2024** (*J. Diff. Geom.*; arXiv:2004.08934 published version) — synthetic timelike Ricci curvature; foundational for the MS framework.
  4. **Connes–Rovelli 1994** (*Class. Quantum Grav.* 11) — thermal time hypothesis; foundational for Paper 48's bridge mechanism.
  5. **Eckstein–Franco 2015** (*Class. Quantum Grav.* 32; arXiv:1212.5171; also "Two roads to noncommutative causality" arXiv:1508.01917) — causality + Lorentzian distance in noncommutative geometry; the closest prior bridge attempt between operator algebras and Lorentzian structure.
- **Novelty verdict**: **NOVEL**. The Wick-rotation functor $W: \mathbf{KreinMetaMet}_{pp} \to \mathbf{LorPLG}_{cov}$, the bridge identification of $\ell^L$ as $\kappa_g \cdot \tau_{mod}^{\omega_W^L}$ on wedge boost orbits, the functorial (not isometric) reading via R2(a)+R2(c) composite, and the "first quantitative pLGH-convergence panel from operator-algebraic input" claim (Paper 45 panel transports verbatim) — none of these appear in the surveyed literature.
- **Citations to add**:
  - Eckstein–Franco arXiv:1212.5171 and arXiv:1508.01917 as the prior operator-algebraic Lorentzian-distance lineage (they did *not* bridge to Mondino–Sämann; the gap is real).
  - Beran–Sämann ("Hyperbolic angles in Lorentzian length spaces and timelike curvature bounds", 2024) for the synthetic-Lorentzian curvature side.
  - Verify the Latrémolière 2512.03573 citation appears in §3 (Krein-pointed proper QMS substrate) as the load-bearing Riemannian-side input.
- **Notes**: This is the most precedence-sensitive of the six papers because the synthetic-Lorentzian-GH program is genuinely active and Mondino–Sämann are publishing at sprint cadence. The bridge structure is what's GeoVac-original; the target category (MS pLGH) is not. **Strongly recommend** explicit one-paragraph "what's new vs MS / what we take from MS" in §1.

### Paper 49 — OSLPLS + twin paradox via Connes–Rovelli thermal time + Uhlmann

- **Closest existing work**:
  1. **Connes 1973** (*Ann. Sci. ENS*) — cocycle Radon–Nikodym derivative, source of the triple-intersection cocycle identity (TICI) construction.
  2. **Bratteli–Robinson Vol. II Thm 5.3.10** (1981) — modular flow on KMS states.
  3. **Uhlmann 1977** (*Comm. Math. Phys.* 54) — relative-entropy monotonicity.
  4. **Connes–Rovelli 1994** — thermal time hypothesis.
  5. **Chirco–Josset–Rovelli 2024** (arXiv:2604.08349; *Thermal time and irreversibility from non-commuting observables*) — closest concurrent work: relative-entropy quantification of irreversibility from modular structures, but addresses non-commuting observables in an accelerated frame, NOT a twin-paradox/cocycle stack across three KMS states. **Does not subsume Paper 49's strict super-additivity result.**
  6. **D'Auria–Liguori 2025** (arXiv:2502.12750, *Thermal time of noncommutative Minkowski spacetime*) — thermal time hypothesis applied to κ- and ρ-Minkowski. Different scope (single noncommutative spacetime, not multi-KMS-state cocycle).
  7. **Lemenant–Rovelli "Thermal time as unsharp observable"** (arXiv:2306.13774, *J. Math. Phys.* 65 (2024)) — observable-theoretic reformulation; orthogonal.
- **Novelty verdict**: **NOVEL**. The synthesis — Connes (cocycle) + Bratteli–Robinson (modular flow) + Connes–Rovelli (thermal time) + Uhlmann (rel-entropy monotonicity) → strict super-additivity of OSLPLS reverse triangle on off-orbit triples — does not appear in the literature surveyed. The twin-paradox-as-quantum-information framing is the kind of substantive new content the math.OA community will recognize as a contribution.
- **Citations to add**:
  - Chirco–Josset–Rovelli arXiv:2604.08349 as the closest concurrent operator-algebraic-relative-entropy-thermal-time thread; explicitly distinguish (their setting: non-commuting observables in single state vs Paper 49's setting: cocycle stack across three KMS states).
  - D'Auria–Liguori arXiv:2502.12750 as the noncommutative-Minkowski thermal-time follow-up to Connes–Rovelli; orthogonal application but worth flagging.
  - **MEMORY.md notes a fix**: per `m1_datta_max_divergence_replacement.md`, the Umegaki relative entropy chain inequality fails by ~5% on TICI; Datta 2009 D_max chain inequality holds. Verify Paper 49 cites **Datta 2009** (*IEEE Trans. Inf. Theory* 55, "Min- and max- relative entropies and a new entanglement monotone") instead of Uhlmann 1977 wherever the cocycle chain inequality is invoked. This is a load-bearing correction.
- **Notes**: This is the most physically-rich of the six and the most likely to attract attention. The Datta max-divergence replacement (per MEMORY) is a non-trivial correctness issue — re-check the proof of Thm 3.3-B2 reads D_max throughout where chain inequality is invoked.

### Paper 50 — F-theorem bit-exact match on truncated spectral triple

- **Closest existing work**:
  1. **Klebanov–Pufu–Safdi 2011** (*JHEP* 10:038; arXiv:1105.4598) — the F-theorem on $S^3$ for conformal scalar / Dirac; closed-form values $F_s = (\log 2)/8 - 3\zeta(3)/(16\pi^2)$ and $F_D = (\log 2)/4 + 3\zeta(3)/(8\pi^2)$ are theirs. **GeoVac's claim is a *match*, not a derivation of new values** — Paper 50's contribution is producing them from a truncated-spectral-triple spectral-zeta side.
  2. **Dowker** (multiple papers, 1990s–2010s) — spectral-zeta derivation of free energies on spheres; foundational for the $-(1/2)\zeta'(0)$ formula. Paper 50 must cite Dowker explicitly.
  3. **Pufu 2017** (*J. Phys. A* 50:443008, "F-theorem and F-maximization") — comprehensive review.
  4. **Casini–Huerta–Myers 2011** (*JHEP* 05:036) — entanglement-entropy route to F-theorem (independent derivation, gives same values).
  5. **Anninos–Anous–Denef** (2021–2024 review papers on sphere partition functions) — adjacent thread.
  6. **arXiv:2603.09799** ("The scheme independent 3-sphere free energy is not a monotone F-function", 2026) — current literature on scheme-dependence of F; should be cited.
  7. **Hekkelman–McDonald arXiv:2412.00628** (*J. Funct. Anal.* 2025) — closest "computation on a truncated spectral triple" concurrent thread, but they compute Szegő-limit noncommutative integrals, not F-theorem coefficients. **No overlap on the F-theorem direction.**
- **Novelty verdict**: **PARTIALLY NOVEL / OVERLAPS-DELIBERATELY**. The values themselves are not new (KPS 2011); the *route* (spectral-zeta on a truncated spectral triple), the *bit-exactness verification* (61+ digits + PSLQ + sympy symbolic), and the *orthogonal M2/M3 dual-basis decomposition* ($F_D + 2F_s = \log(2)/2$ isolates M2; $F_D - 2F_s = 3\zeta(3)/(4\pi^2)$ isolates M3) are GeoVac-original. The §7 extension to $S^5$ with new closed forms ($F^{S^5}$ for scalar via PSLQ $[32, -2, -2, 15]$; analytical Dirac $D'^{S^5}(0)$) appears genuinely new to the literature surveyed — search did not return concurrent $F^{S^5}$ closed forms.
- **Citations to add**:
  - Dowker (multiple papers) — must be cited for the spectral-zeta route on spheres.
  - Casini–Huerta–Myers 2011 (arXiv:1102.0440) — entanglement-entropy route, independent verification of KPS values.
  - arXiv:2603.09799 (2026 scheme-independent F-function paper) — current literature on F-theorem monotonicity.
  - Anninos–Anous–Denef sphere-partition-function reviews — for orientation in the broader community.
- **Notes**: The strongest claim is the *orthogonal M2/M3 decomposition* (the dual-basis projection isolating $\log 2$ vs $\zeta(3)$) — that's where Paper 50 contributes structural insight beyond KPS. Make sure §3.4–§3.5 emphasize this as the novel content; the bit-exact match is the *verification*, not the *claim*. Frame the paper this way in revision: "we verify the KPS values from a structurally different route, and the verification exposes an orthogonal M2/M3 decomposition that is invisible from the localization / entanglement-entropy routes."

---

## Cross-cutting observations

1. **Lorentzian propinquity (Papers 45/46/47) is GeoVac-only as of June 2026.** The synthetic-Lorentzian-GH program (Mondino–Sämann, Braun–Sämann, Cavalletti–Mondino) is actively publishing in 2025–2026 but is categorically separate (synthetic differential geometry, not operator algebras). The almost-commutative-twisted-Wick-rotation program (Nieuviarts) is the only spectral-triple-side concurrent thread but addresses *emergence* (Riemannian → Lorentzian via twist), not *convergence* (truncated → continuum). **The novelty claims for Papers 45/46/47 hold.**

2. **Paper 48 is the precedence-sensitive bridge.** The synthetic side (MS pLGH) is publicly published and developing fast. Paper 48's contribution is the *operator-algebraic input → synthetic output* direction. **Recommend writing §1 with explicit "what we take from MS / what we contribute to MS" framing before arXiv submission.**

3. **Paper 49 OSLPLS / twin paradox is the most physically novel.** The combination of Connes–Rovelli thermal time + Uhlmann (or Datta — see correction note) relative entropy on a cocycle stack across three KMS states is not in the surveyed literature. The closest concurrent thread is Chirco–Josset–Rovelli 2024 (arXiv:2604.08349), which is operator-algebraically adjacent but addresses a different question.

4. **Paper 50 is a verification-with-structural-insight, not a new F-theorem.** Frame it as such. The KPS 2011 values are well-known; the spectral-zeta route via $-(1/2)\zeta'(0)$ is standard (Dowker); the *truncated-spectral-triple* side is GeoVac. The M2/M3 orthogonal decomposition is the substantive contribution. Get the framing right in §1 to avoid reviewer pushback that "this is a check, not a result."

5. **The Latrémolière 2512.03573 hypertopology (Dec 2025) is the single most-relevant concurrent work across the whole arc.** It does NOT subsume any of Papers 45–49 (Riemannian, non-compact, metric-level), but it IS the load-bearing input for Paper 48's Krein-lift and an adjacent thread for Paper 47's non-compact direction. **Already acknowledged in CLAUDE.md / Paper 47 §1.1 / Paper 48 §3.** Re-verify the citation is in place in every paper that touches non-compact convergence.

6. **Mondino–Sämann arXiv:2504.10380 v4 (Dec 2025) and Braun–Sämann arXiv:2506.10852 (June 2025) are NEW since CLAUDE.md's May-2026 concurrent-work audit.** Both papers extend the synthetic Lorentzian GH program in 2025. Neither overlaps GeoVac's operator-algebraic claims, but both should be cited in Paper 48 (and Paper 45 §1.4 footnote) as current MS-program literature.

---

## Limitations

1. **No direct PDF access**: relied on WebSearch + WebFetch on abstract pages. A more thorough audit would download and search the full PDFs of the closest hits, particularly Latrémolière 2512.03573 (85 pages) and the Mondino–Sämann 2025 papers.
2. **Journal-only items**: did not check zbMATH / MathSciNet / Math. Reviews for items that may not be on arXiv (rare in this community, but possible).
3. **The Hekkelman/Hekkelman-McDonald thread** has more recent items on Hekkelman's homepage (https://www.emhekkelman.nl/) than the abstract page exposes. Recommend a targeted follow-up directly browsing her preprint list and Walter van Suijlekom's group page.
4. **The Eckstein–Franco / Bochniak / Polish school** on causality + noncommutative geometry has many published-but-not-on-arXiv items in *Class. Quantum Grav.* and *J. Phys. A*. Recommend a follow-up audit specifically against their JFA / CMP / PoS papers 2023–2026.
5. **The 2026 listings on math.OA** since v3.45 of CLAUDE.md was finalized (May 2026 audit window) have not been systematically swept day-by-day. Recommend a one-week-pre-submission re-audit for each of Papers 45–50 before arXiv upload.
6. **Datta max-divergence correction**: per MEMORY.md item `m1_datta_max_divergence_replacement.md`, Paper 49's Thm 3.3-B2 should be re-checked to ensure it uses Datta's $D_{max}$ chain rule (Datta 2009, IEEE Trans. Inf. Theory 55) rather than Umegaki's relative-entropy chain (which fails ~5% on TICI). This is a substantive correctness item, not just a citation update — recommend the PM verify directly before submission.

---

**Bottom line**: All six novelty claims survive the June 2026 audit. Papers 45–47 stand cleanly as first-in-literature in the Lorentzian-propinquity / Krein-spectral-triple-convergence space. Paper 48 needs framing care vs the active MS program. Paper 49 needs the Datta correction verified. Paper 50 needs framing as verification-with-structural-insight, not a new F-theorem. Recommended new citations: Mondino–Sämann 2504.10380 v4, Braun–Sämann 2506.10852, Behrndt–Holzmann–Stelzer 2404.07784, Chirco–Josset–Rovelli 2604.08349, Datta 2009 IEEE-IT-55, arXiv:2603.09799 scheme-independent F.
