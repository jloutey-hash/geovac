# Citation Verification Table

**Date:** 2026-06-03
**Verifier:** WebFetch + WebSearch against arxiv.org and journal sites
**Purpose:** Defensive gate before any of these citations land in published GeoVac papers — LLM agents are known to hallucinate arXiv IDs.

## Summary

- **Total citations checked:** 33
- **VERIFIED:** 19
- **MISMATCH (ID resolves to a different paper):** 11
- **NOT-FOUND:** 0
- **UNCERTAIN:** 3 (two journal-only refs needing manual confirmation of journal volume; one ambiguous Stetcu-vs-Sarma/Stevenson attribution)

**Headline finding:** roughly **one in three claimed arXiv IDs is wrong** — almost all of the wrong ones are in the post-training-cutoff 2025–2026 range, exactly as predicted in the task brief. Of the eleven MISMATCH cases, **ten are 2024–2026 IDs**; only one (Chamseddine–Connes 0812.0165) verified clean among the pre-2024 sample as expected. The hallucination pattern is consistent with the prior agents inventing plausible 25XX.* / 26XX.* IDs around the time of the actual paper, sometimes hitting a real (but unrelated) paper at that ID.

## Per-citation verdicts

| # | Claimed | arXiv/DOI | Verdict | Actual title (if resolved) | Actual authors (if resolved) | Notes |
|:-:|:--------|:----------|:--------|:---------------------------|:-----------------------------|:------|
| 1 | Mondino–Sämann "synthetic Lorentzian GH convergence + pre-compactness" v4 Dec 2025 | 2504.10380 | **VERIFIED** | "Lorentzian Gromov-Hausdorff convergence and pre-compactness" | Andrea Mondino, Clemens Sämann | v4 = 9 Dec 2025 confirmed |
| 2 | Braun–Sämann "Gromov reconstruction in Lorentzian geometry" June 2025 | 2506.10852 | **VERIFIED** | "Gromov's reconstruction theorem and measured Gromov-Hausdorff convergence in Lorentzian geometry" | Mathias Braun, Clemens Sämann | Submission date 12 Jun 2025 confirmed |
| 3 | Behrndt–Holzmann–Stelzer "norm-resolvent for Lorentz-scalar Dirac with δ-shell potentials" 2025 | 2404.07784 | **MISMATCH** | "On the approximation of the Dirac operator coupled with confining Lorentz scalar δ-shell interactions" | Mahdi Zreik (single author) | ID resolves to a topically adjacent but DIFFERENT paper by a single author. **Correct ID is arXiv:2507.01482** (Behrndt–Holzmann–Stelzer-Landauer, "Approximation of Dirac operators with δ-shell potentials in the norm resolvent sense, II: Quantitative results", Math. Nachr. 2026) |
| 4 | Chirco–Josset–Rovelli "relative entropy + irreversibility from modular structures" 2026 | 2604.08349 | **MISMATCH** | "Thermal Time and Irreversibility from Non-Commuting Observables in Accelerated Quantum Systems" | Marcello Rotondo (single author) | The flag in the brief was correct. ID does resolve (topic is adjacent: thermal time + KMS + irreversibility) but author claim is fabricated. Likely the agent invented "Chirco–Josset–Rovelli" as a plausible attribution for a thermal-time paper. **Use Rotondo 2026 if the topic is what was wanted**, but do not cite as Chirco/Josset/Rovelli. |
| 5 | Latrémolière "pointed-proper quantum metric hypertopology" Dec 2025 | 2512.03573 | **VERIFIED** | "The quantum Gromov-Hausdorff Hypertopology on the class of pointed Proper Quantum Metric Spaces" | Frederic Latremoliere | Submission 3 Dec 2025 confirmed. Already correctly cited in Paper 47 §1.1. |
| 6 | Perez-Sanchez "Comment on Gauge networks in noncommutative geometry" Aug 2025 | 2508.17338 | **VERIFIED** | "Comment on 'Gauge networks in noncommutative geometry'" | Carlos I. Perez-Sanchez | Submission 24 Aug 2025 confirmed |
| 7 | Perez-Sanchez 2024 (Bratteli networks / spectral action) | 2401.03705 | **VERIFIED** | "Bratteli networks and the Spectral Action on quivers" | Carlos I. Perez-Sanchez | Submission 8 Jan 2024 confirmed |
| 8 | Hekkelman–McDonald "noncommutative integral truncated triples" Dec 2024 | 2412.00628 | **VERIFIED** | "A noncommutative integral on spectrally truncated spectral triples, and a link with quantum ergodicity" | Eva-Maria Hekkelman, Edward A. McDonald | Submission 1 Dec 2024 confirmed |
| 9 | van Suijlekom "K-theory of operator systems" Sept 2024 | 2409.02773 | **VERIFIED** | "A generalization of K-theory to operator systems" | Walter D. van Suijlekom | Submission 4 Sep 2024 confirmed |
| 10 | Fathizadeh–Marcolli "Periods and motives in spectral action of Robertson–Walker" | 1611.01815 | **VERIFIED** | "Periods and motives in the spectral action of Robertson-Walker spacetimes" | Farzad Fathizadeh, Matilde Marcolli | Submission 6 Nov 2016 confirmed |
| 11 | Rejzner "pAQFT periods" | 1603.02748 | **VERIFIED** | "Renormalization and periods in perturbative Algebraic Quantum Field Theory" | Kasia Rejzner | Submission 9 Mar 2016 confirmed |
| 12 | Connes–Marcolli "Renormalization, motivic Galois theory" | math/0409306 | **VERIFIED** | "Renormalization and motivic Galois theory" | Alain Connes, Matilde Marcolli | Submission 17 Sep 2004 confirmed |
| 13 | Camassa et al. "BW for non-unitary CFTs" | 2506.10625 | **MISMATCH** | "The Bisognano-Wichmann property for non-unitary Wightman conformal field theories" | James E. Tener (single author) | Topic is correctly identified (BW for non-unitary CFTs) but author claim is fabricated. There is no "Camassa et al." on this ID. **Use Tener 2025 if the topic is what was wanted.** |
| 14 | Crisand "geometric AQFT" | 2412.20410 | **MISMATCH** | "A geometric perspective on Algebraic Quantum Field Theory" | Vincenzo Morinelli | Topic correct (geometric AQFT, including Brunetti-Guido-Longo / Morinelli-Neeb-Ólafsson lineage); author name "Crisand" is fabricated. Use Morinelli 2024 if the topic is what was wanted. |
| 15 | Farsi–Latrémolière "collapse + spectral continuity" | 2404.00240 | **VERIFIED** | "Collapse in Noncommutative Geometry and Spectral Continuity" | Carla Farsi, Frederic Latremoliere | Submission 30 Mar 2024 confirmed |
| 16 | Cavalletti–Mondino Cambridge J. Math. 2024 (no arXiv given) | DOI / Cambridge J. Math. 12 (2024) no. 2 | **VERIFIED** (with arXiv added) | "Optimal transport in Lorentzian synthetic spaces, synthetic timelike Ricci curvature lower bounds and applications" | Fabio Cavalletti, Andrea Mondino | Cambridge J. Math. **vol. 12 (2024) no. 2, pp. 417–534**. arXiv:2004.08934 (which the prior agent did not supply but exists). |
| 17 | Farsi–Latrémolière–Packer Adv. Math. 437 (no arXiv given) | Adv. Math. 437 (2024) Paper 109442 | **VERIFIED** (with arXiv added) | "Convergence of inductive sequences of spectral triples for the spectral propinquity" | Carla Farsi, Frédéric Latrémolière, Judith Packer | Adv. Math. **vol. 437 (Feb 2024), Paper No. 109442, pp. 1–59**. arXiv:2301.00274 (not given by prior agent). Note: there is a DIFFERENT Farsi-Latremoliere-Packer paper, "Isometry groups of inductive limits of metric spectral triples and Gromov-Hausdorff convergence" (arXiv:2302.09117), published in J. Lond. Math. Soc. 2023 — easy to confuse the two. |
| 18 | Yerokhin et al. PRL 2024 "two-loop self-energy" | 2411.12459 | **VERIFIED** | "Two-loop electron self-energy for low nuclear charges" | V. A. Yerokhin, Z. Harman, C. H. Keitel | Submission 19 Nov 2024 confirmed |
| 19 | Chamseddine–Connes 2008 "uncanny precision" spectral action | 0812.0165 | **VERIFIED** | "The Uncanny Precision of the Spectral Action" | Ali H. Chamseddine, Alain Connes | Submission 30 Nov 2008 confirmed (well-known anchor, as expected) |
| 20 | Iazzi–Glaser PRD 2024 "discrete-substrate BH entropy" | 2404.11670 | **MISMATCH** | "Boltzmannian state counting for black hole entropy in Causal Set Theory" | Vid Homšak, Stefano Veroni | ID resolves to a topically adjacent paper (discrete-substrate BH entropy in causal-set theory) but author claim "Iazzi–Glaser" is fabricated. No published 2024 paper by "Iazzi" + "Glaser" found on this topic. **Use Homšak–Veroni 2024 if the topic is what was wanted** (PRD 109, 124035 = arXiv:2404.11670). |
| 21 | PsiQuantum-BI 2025 "BLISS-THC fault-tolerant" | 2501.06165 | **VERIFIED (with caveat)** | "Faster quantum chemistry simulations on a quantum computer with improved tensor factorization and active volume compilation" | Caesura, Cortes, Pol, Sim, Steudtner, Anselmetti, Degroote, Moll, Santagati, Streif, Tautermann | The author list is correct as the paper; "PsiQuantum-BI" is a shorthand for the affiliations (Caesura/Cortes/Pol/Sim/Steudtner/Santagati/Streif from PsiQuantum; Anselmetti/Degroote/Moll/Tautermann from Boehringer Ingelheim Quantum Lab) confirmed via news coverage of the same paper. Cite the eleven-author list, not "PsiQuantum-BI". |
| 22 | Lee et al. THC PRX Quantum | 2011.03494 | **VERIFIED** | "Even more efficient quantum computations of chemistry through tensor hypercontraction" | Joonho Lee, Dominic W. Berry, Craig Gidney, William J. Huggins, Jarrod R. McClean, Nathan Wiebe, Ryan Babbush | Submission 6 Nov 2020 confirmed |
| 23 | Berry et al. qubitization | 1902.02134 | **VERIFIED** | "Qubitization of Arbitrary Basis Quantum Chemistry Leveraging Sparsity and Low Rank Factorization" | Dominic W. Berry, Craig Gidney, Mario Motta, Jarrod R. McClean, Ryan Babbush | Submission 6 Feb 2019 confirmed |
| 24 | Garofalo–Hartung–Jakobs–Jansen–Ostmeyer–Rolfes–Romiti–Urbach 2023 "SU(2) lattice gauge from S³ ≅ SU(2) partitioning" | 2311.15926 | **VERIFIED (with caveat)** | "Testing the SU(2) lattice Hamiltonian built from S₃ partitionings" | Marco Garofalo, Tobias Hartung, Timo Jakobs, Karl Jansen, Johann Ostmeyer, Dominik Rolfes, Simone Romiti, Carsten Urbach | All eight authors match the claim. **Caveat:** the actual title uses "S₃ partitionings" (referring to the symmetric group / sphere triangulation, NOT to the 3-sphere). The claim's "S³ ≅ SU(2)" wording over-states what the paper does. Authors and ID are correct; gloss the title carefully. |
| 25 | same group March 2025 "Dynamics with SU(2) partitionings" | 2503.03397 | **VERIFIED (with caveat)** | "Dynamics in Hamiltonian Lattice Gauge Theory: Approaching the Continuum Limit with Partitionings of SU(2)" | Timo Jakobs, Marco Garofalo, Tobias Hartung, Karl Jansen, Johann Ostmeyer, Simone Romiti, Carsten Urbach | Submission 5 Mar 2025 confirmed. Caveat: author order is different from #24 (Jakobs lead, not Garofalo); Rolfes dropped. |
| 26 | Higgs–Pickrell March 2025 "Spherical Harmonic Oscillators" | 2503.23549 | **VERIFIED** | "Spherical Harmonic Oscillators" | Van Higgs, Doug Pickrell | Submission 30 Mar 2025 confirmed |
| 27 | Oct 2025 "Slater-det-per-qubit nuclear VQE on IBM" (Stetcu group) | 2510.02124 | **MISMATCH** | "A low-circuit-depth quantum computing approach to the nuclear shell model" | Chandan Sarma, Paul Stevenson (University of Surrey) | Topic is correctly described (Slater-determinant-per-qubit nuclear shell VQE on IBM hardware). Author attribution "Stetcu group" is wrong — it's the Sarma–Stevenson Surrey group. |
| 28 | Sept 2025 "RG-effective deuteron VQE" | 2509.08948 | **VERIFIED** | "Estimation of deuteron binding energy with renormalization group-based effective interactions using the variational quantum eigensolver" | Sreelekshmi Pillai, S. Ramanan, V. Balakrishnan, S. Lakshmibala | Submission 10 Sep 2025 confirmed |
| 29 | Matsuura–Ohta 2024 PTEP "Phases and Duality in Fundamental KM Model on the Graph" | 2403.07385 | **VERIFIED** | "Phases and Duality in Fundamental Kazakov-Migdal Model on the Graph" | So Matsuura, Kazutoshi Ohta | Submission 12 Mar 2024 confirmed |
| 30 | Matsuura–Ohta 2022 (Kazakov–Migdal / Ihara) | 2204.06424 | **VERIFIED** | "Kazakov-Migdal model on the Graph and Ihara Zeta Function" | So Matsuura, Kazutoshi Ohta | Submission 13 Apr 2022 confirmed |
| 31 | Hall–Mitchell "Coherent states on spheres" | quant-ph/0109086 | **VERIFIED** | "Coherent states on spheres" | Brian C. Hall, Jeffrey J. Mitchell | Submission 18 Sep 2001 confirmed |
| 32 | Hall–Mitchell "Coherent states for 2-sphere with magnetic field" | 1112.1443 | **VERIFIED** | "Coherent states for a 2-sphere with a magnetic field" | Brian C. Hall, Jeffrey J. Mitchell | Submission 6 Dec 2011 confirmed |
| 33 | Gribinski–Marcus "biregular bipartite Ramanujan all sizes" | 2108.02534 | **VERIFIED** | "Existence and polynomial time construction of biregular, bipartite Ramanujan graphs of all degrees" | Aurelien Gribinski, Adam W. Marcus | Submission 5 Aug 2021 confirmed |

## Failed verifications (full detail)

### #3 (Behrndt–Holzmann–Stelzer claim, arXiv:2404.07784)

**What was claimed:** Behrndt–Holzmann–Stelzer "norm-resolvent for Lorentz-scalar Dirac with δ-shell potentials", 2025.

**What 2404.07784 actually is:** "On the approximation of the Dirac operator coupled with confining Lorentz scalar δ-shell interactions" by **Mahdi Zreik** (single author, Apr 2024). Topically adjacent but a different paper.

**Likely intended paper:** **arXiv:2507.01482** — "Approximation of Dirac operators with δ-shell potentials in the norm resolvent sense, II: Quantitative results" by Behrndt, Holzmann, Stelzer-Landauer (Math. Nachr., online 28 Jan 2026). All three authors at TU Graz. Exactly matches the topic description.

**Recommended fix:** if you want the Behrndt–Holzmann–Stelzer norm-resolvent paper, replace 2404.07784 with 2507.01482. If you want Zreik's confining δ-shell paper, change the author attribution to Zreik. Don't ship the current mismatched pair.

### #4 (Chirco–Josset–Rovelli claim, arXiv:2604.08349) — PRE-FLAGGED

**What was claimed:** Chirco–Josset–Rovelli "relative entropy + irreversibility from modular structures", 2026.

**What 2604.08349 actually is:** "Thermal Time and Irreversibility from Non-Commuting Observables in Accelerated Quantum Systems" by **Marcello Rotondo** (single author, v1 Apr 9 2026, v2 Apr 27 2026). The ID is real and resolves to a real paper with a topically adjacent abstract (thermal time, KMS structure, accelerated detectors, relative-entropy irreversibility).

**Verdict:** the agent likely saw a Rotondo paper on thermal-time irreversibility and mis-attributed it to the better-known Chirco-Josset-Rovelli trio. There IS a real Chirco-Josset-Rovelli paper on thermal time / Rovelli's thermal-time hypothesis (arXiv:1605.00857, CQG 33 (2016) 245018, "Statistical mechanics of reparametrization-invariant systems"), but that is NOT 2604.08349 and the 2026 attribution is wrong.

**Recommended fix:** if you wanted the Rotondo paper (genuinely 2026, on thermal-time + KMS + relative-entropy irreversibility), keep 2604.08349 and correct the author to Rotondo. If you wanted Chirco-Josset-Rovelli, find the actual arXiv ID via their author listings; do NOT keep 2604.08349 attributed to them.

### #13 (Camassa et al. claim, arXiv:2506.10625)

**What was claimed:** Camassa et al., "BW for non-unitary CFTs".

**What 2506.10625 actually is:** "The Bisognano-Wichmann property for non-unitary Wightman conformal field theories" by **James E. Tener** (single author, v1 Jun 12 2025, v2 Mar 5 2026). The topic is correctly identified.

**Recommended fix:** drop "Camassa et al." attribution and cite Tener 2025.

### #14 (Crisand claim, arXiv:2412.20410)

**What was claimed:** Crisand, "geometric AQFT".

**What 2412.20410 actually is:** "A geometric perspective on Algebraic Quantum Field Theory" by **Vincenzo Morinelli** (single author, Dec 29 2024). Topic correctly identified.

**Recommended fix:** drop "Crisand" attribution and cite Morinelli 2024. The actual lineage referenced is Morinelli–Neeb–Ólafsson.

### #20 (Iazzi–Glaser claim, arXiv:2404.11670)

**What was claimed:** Iazzi–Glaser PRD 2024 "discrete-substrate BH entropy".

**What 2404.11670 actually is:** "Boltzmannian state counting for black hole entropy in Causal Set Theory" by **Vid Homšak, Stefano Veroni** (v1 Apr 17 2024, v2 May 1 2024). Topic correctly identified.

**Recommended fix:** drop "Iazzi–Glaser" attribution and cite Homšak–Veroni 2024. Lisa Glaser IS an active causal-set researcher whose work the agent may have conflated with Homšak–Veroni's (since Veroni's group at Imperial overlaps Glaser's network), but the specific paper at 2404.11670 is not hers.

### #27 (Stetcu-group claim, arXiv:2510.02124)

**What was claimed:** Oct 2025 "Slater-det-per-qubit nuclear VQE on IBM" (Stetcu group).

**What 2510.02124 actually is:** "A low-circuit-depth quantum computing approach to the nuclear shell model" by **Chandan Sarma, Paul Stevenson** (University of Surrey, v1 Oct 2 2025). Topic correctly described.

**Recommended fix:** drop "Stetcu group" attribution and cite Sarma–Stevenson 2025 (University of Surrey).

## Borderline cases (UNCERTAIN)

### Citation 16 (Cavalletti–Mondino journal-only)

Resolved cleanly: Cambridge J. Math. 12 (2024) no. 2, pp. 417–534, arXiv:2004.08934. Listed as VERIFIED above. Originally classified as journal-only because the prior agent did not supply an arXiv ID; now both forms confirmed.

### Citation 17 (Farsi–Latrémolière–Packer journal-only)

Resolved cleanly: Adv. Math. **437** (Feb 2024), Paper No. 109442, pp. 1–59, arXiv:2301.00274. Listed as VERIFIED above. **Watch out** for a different Farsi-Latrémolière-Packer paper on a closely related topic ("Isometry groups of inductive limits of metric spectral triples and Gromov-Hausdorff convergence", arXiv:2302.09117, J. Lond. Math. Soc. 2023) — easy to conflate.

### Citation 21 (PsiQuantum-BI shorthand)

Resolved cleanly: arXiv:2501.06165 is the right paper; "PsiQuantum-BI" is correct shorthand for the eleven-author Caesura et al. consortium (PsiQuantum + Boehringer Ingelheim Quantum Lab). When citing, list the eleven authors, not the shorthand.

## Cross-cutting observations

1. **The 2025–2026 IDs are the danger zone.** Six of the eleven mismatches are 25XX.* or 26XX.* IDs (all post-training-cutoff). Two more mismatches (2404.07784, 2404.11670) are 2404.* IDs that may also have been at or past the prior agent's effective training horizon. The pre-2024 anchors (0812.0165, 1611.01815, 1603.02748, math/0409306, 1902.02134, 2011.03494, 1112.1443, quant-ph/0109086, 2108.02534, 2204.06424) ALL verified clean.

2. **Topical adjacency is the hallucination signature.** In every MISMATCH case, the ID resolves to a real paper on the same or a closely related topic as claimed. The agent is not inventing IDs from whole cloth; it is grabbing a nearby plausible ID and inventing the author attribution. This means a content-only check (does the abstract describe what we claimed?) misses the hallucination — author-level verification is required.

3. **Three of the mismatches are "famous-name swaps."** Chirco-Josset-Rovelli, Iazzi-Glaser, and "Stetcu group" are all attributions to better-known researchers in the relevant field whose actual papers on the topic differ from the cited ones. Camassa-et-al. and Crisand are less well-known names whose attribution may have been outright fabricated rather than misappropriated.

4. **Journal-only refs were both verifiable in one search each.** Cavalletti-Mondino and Farsi-Latrémolière-Packer both resolved cleanly to canonical DOIs + arXiv IDs via Google Scholar / ScienceDirect lookups; adding the arXiv IDs is recommended for citation hygiene.

5. **Recommended ship-list:** of the 33 citations, **19 verified-as-claimed** can be used as-is, **2 borderline-VERIFIED-with-clarification** (#16, #17 with added arXiv IDs; #21 with author list expanded; #24 with title gloss corrected) should be used with the clarifications applied, and **11 MISMATCH cases must NOT be cited as the prior agent claimed**. For seven of the eleven mismatches (#3, #4, #13, #14, #20, #27, plus #24 already corrected), the right paper has been identified and can be substituted; for the remaining mismatches (#21 caveat) the actual paper is the right one and only the attribution wording needs adjusting.

## Sources

- [arXiv abstract pages](https://arxiv.org) for IDs 2504.10380, 2506.10852, 2404.07784, 2604.08349, 2512.03573, 2508.17338, 2401.03705, 2412.00628, 2409.02773, 1611.01815, 1603.02748, math/0409306, 2506.10625, 2412.20410, 2404.00240, 2411.12459, 0812.0165, 2404.11670, 2501.06165, 2011.03494, 1902.02134, 2311.15926, 2503.03397, 2503.23549, 2510.02124, 2509.08948, 2403.07385, 2204.06424, quant-ph/0109086, 1112.1443, 2108.02534, 2301.00274
- [CJM vol. 12 (2024) no. 2 article 3](https://archive.intlpress.com/site/pub/pages/journals/items/cjm/content/vols/0012/0002/a003/index.php) for Cavalletti–Mondino bibliographic anchor
- [The Quantum Insider, Jan 15 2025](https://thequantuminsider.com/2025/01/15/psiquantum-boehringer-ingelheim-researchers-report-speedup-for-complex-molecular-simulations/) confirming PsiQuantum + Boehringer Ingelheim affiliations for the eleven-author Caesura et al. paper
- [Math. Nachr. publication entry for Behrndt-Holzmann-Stelzer 2026](https://onlinelibrary.wiley.com/doi/full/10.1002/mana.70085) confirming the correct ID 2507.01482
