# Phase 4 Long-Tail Connection Scout

Date: 2026-06-03

## Methodology

Four parallel web search sweeps (≈16 queries total) across 2015–2026, refined with
targeted follow-ups on specific authors/arXiv IDs once initial hits surfaced. For
each theme I aimed to identify (a) one or two specific recent papers that are
direct adjacencies (not just topical overlap), (b) the senior author or group
behind those papers as a potential conversation partner, and (c) what each
connection would actually buy GeoVac (citation hardening, response target,
collaboration angle, or a structural cross-check).

All arXiv IDs below were surfaced by search and confirmed via the search-result
snippets and abstracts. I did not download PDFs; claims about paper contents are
those reported in the abstracts/result snippets and should be re-verified before
being cited in a paper. I flag uncertainty where it exists.

The Marcolli–van Suijlekom 2014 + Pérez-Sánchez 2024/2025 lineage was already
known to CLAUDE.md and is cross-cut into Themes 3 and 4 below; the new findings
are how active and graph-shaped the post-2022 work is.

---

## Per-theme findings

### Theme 1 — Bargmann-Segal lattice (Paper 24)

**Adjacent active programs.**

- **Van Higgs & Doug Pickrell (Arizona, math).** "Spherical Harmonic Oscillators,"
  arXiv:2503.23549 (Mar 2025, rev Oct 2025). Constructs an SO(d)-invariant
  self-adjoint second-order operator on the d-sphere with the d-dim spherical
  analogue of the harmonic oscillator; ground state (1+r²)^{−ω}, compact resolvent
  (discrete spectrum), spectrum multiplicity-free for ω > 0. Explicitly discusses
  a "holomorphic Hermitian line bundle with canonical connection" — i.e. the
  Hardy / Segal–Bargmann side. Higgs is junior; Pickrell is the senior name.
- **Brian Hall (Notre Dame).** The canonical Segal–Bargmann-on-compact-Lie-group
  thread (1994; 2002 quant-ph/0109086 "Coherent states on spheres" with J.
  Mitchell; 2012 arXiv:1112.1443 "Coherent states for a 2-sphere with a magnetic
  field"). Hall's program is the literature root that Paper 24 should cite.
- **Quivers + AF-spectral-triple lineage** (Pérez-Sánchez 2024 + 2022 Lifting
  Bratteli Diagrams arXiv:2207.04466) — adjacent because Bratteli/AF is the
  natural inductive-limit setting in which one might lift Paper 24's
  finite-N_max Bargmann-Segal lattice toward a continuum spectral triple.
- **Coupled-SUSY Segal–Bargmann (Cameron L. Williams).** arXiv:2110.08995 (2021,
  pub. 2022 in CAOT). Direct generalization of SB to coupled-SUSY oscillators;
  niche but very close in shape to the Paper 24 setup.

**Specific papers to cite or build bridges to.**

- arXiv:2503.23549 (Higgs–Pickrell 2025): natural §III citation for Paper 24
  where the SO(d)-invariant HO operator and its ground state (1+r²)^{−ω} are
  discussed — Paper 24's Hardy-sector ground state is the same object in the S⁵
  case. Their "potentially relevant to chiral models in 2d QFT" line is a
  conversation hook: GeoVac is asking what the d ≥ 3 case is for, and chiral 2d
  is interesting but not where the framework lives.
- arXiv:quant-ph/0109086 (Hall–Mitchell, "Coherent states on spheres"): the
  canonical SB-on-spheres reference. Paper 24's bibliography should ground its
  Hardy-sector framing on this paper if it doesn't already.
- arXiv:1112.1443 (Hall–Mitchell, 2012): adds a magnetic field — the natural
  toy model for a connection / Wilson U(1) on the S^d HO Hardy sector, i.e. the
  HO analog of the Coulomb-side Paper 25/30 gauge story.
- arXiv:2110.08995 (Williams, coupled-SUSY SB): GeoVac's HO ladder + Bargmann
  lattice is in the same algebraic neighborhood; worth pinging.

**Potential collaborators / response targets.**

- **Doug Pickrell (Arizona)** is the natural conversation partner — his March
  2025 paper is the closest published thing to Paper 24's S^d HO Hardy
  construction, and he explicitly raised the d ≥ 3 question.
- **Brian Hall (Notre Dame)** would be the canonical math.OA / coherent-states
  reviewer for Paper 24 if PI wants to get expert eyes outside the GeoVac
  ecosystem before any Zenodo push beyond the current state. Hall's program is
  exactly the lineage Paper 24's HO rigidity theorem sits in.

**What each connection offers.**

- Pickrell: structural cross-check on the d = 5 spectrum + ground-state form,
  and a potential answer to "what does Paper 24's HO rigidity theorem look like
  to a working geometric-quantization person?"
- Hall: literature root + senior endorsement candidate; Paper 24's claim of
  π-freeness on the Hardy sector is a sharper statement than I see anywhere in
  Hall's chain and would interest him.
- Williams (coupled SUSY): low-cost ping; very direct algebraic adjacency.
- Bratteli/AF lineage (Pérez-Sánchez, the 2207.04466 thread): the natural route
  to lift Paper 24 to a continuum-limit spectral triple if and when WH1's
  HO-side analog is sprinted.

---

### Theme 2 — Nuclear shell model qubit Hamiltonians (Paper 23)

**Adjacent active programs.**

- **Slater-determinant-per-qubit encoding (Pittsburgh / IBM lineage)** —
  arXiv:2510.02124 (Oct 2025, pub Springer Discover Quantum Science 2026).
  Maps each Slater determinant to a qubit (rather than each single-particle
  state). Tested on ⁶,⁷,⁸,⁹Li, ¹⁸F, ²¹⁰Po, ²¹⁰Pb on IBM ibm_pittsburgh hardware
  with ZNE error mitigation; < 4% deviation from shell model. Directly
  comparable to Paper 23's Jordan–Wigner-per-orbital choice.
- **Gray-code encoding lineage** — arXiv:2504.11689 (April 2025, "Advancing
  quantum simulations of the nuclear shell model with Gray-code-based
  resource-efficient protocols"), with earlier roots in Romero et al. and
  Stetcu et al. ³⁸Ar, ⁶Li with VQE / qubit-ADAPT-VQE / VQD on noisy and noisy +
  mitigated hardware. This is the principal qubit-count-vs-circuit-depth
  alternative to Paper 23's encoding.
- **Comparison of VQE algorithms in p-shell nuclei** — arXiv:2507.13819 (Aug
  2025, pub Phys. Rev. C). UCC vs ADAPT on ⁶He → ¹⁰B, with the finding that
  Slater-determinant references are best for both and ADAPT wins near magic,
  UCC wins mid-shell. This is the natural benchmark frame for Paper 23 to
  position against.
- **Trapped-ion hardware extension** — arXiv:2509.20642 (Sept 2025, "Bridging
  Quantum Computing and Nuclear Structure: Atomic Nuclei on a Trapped-Ion
  Quantum Computer"). Adds trapped-ion as a second hardware platform alongside
  the IBM superconducting work.
- **Deuteron binding via RG-effective interactions + VQE** —
  arXiv:2509.08948 (Sept 2025). Same molecule as Paper 23's headline; very
  natural direct comparison target.
- **Trapped-ion two-nucleon single-step** — arXiv:2512.12798 (Dec 2025,
  "Single-step Quantum Simulation of Two Nucleons"). Two-nucleon at single
  step on trapped ions.
- **¹²C rotational band qubit mapping** — arXiv:2412.06979 (Dec 2024, "Efficacious
  qubit mappings for quantum simulations of the ¹²C rotational band").

**Specific papers to cite or build bridges to.**

- arXiv:2510.02124 (Slater-det-per-qubit, 2025): direct comparison target for
  Paper 23's "deuteron 16 qubits, 592 Pauli" headline. Their 22- and 29-qubit
  ²¹⁰Po/²¹⁰Pb on hardware is a useful upper-end anchor — Paper 23's "26
  qubits, 614 Pauli" composed nuclear-electronic POC is in the same ballpark.
- arXiv:2509.08948 (deuteron binding via RG, 2025): direct deuteron-binding
  competitor; Paper 23's Minnesota + Moshinsky–Talmi approach is structurally
  different and the comparison is informative for both sides.
- arXiv:2507.13819 (UCC vs ADAPT in p-shell, 2025): the natural frame inside
  which to position Paper 23's resource counts.
- arXiv:2401.03705 (Pérez-Sánchez Bratteli networks) + arXiv:2508.17338
  (Pérez-Sánchez comment on Marcolli–vS): cross-pollination — Paper 23's
  composed nuclear-electronic POC has the same shape as a quiver / Bratteli
  network with two factor algebras (a nuclear shell C* and an electronic
  C*-truncation joined by a hyperfine "edge"). This is a Theme-4 connection
  too.
- Moshinsky-bracket recent computational work: ScienceDirect Comp. Phys. Comm.
  paper 2024 ("Moshinsky brackets for a wide range of quantum numbers using
  generating functions") and arXiv:2104.12515 (2021). Both are relevant for
  Paper 23's appendix on Moshinsky–Talmi bracket computation.

**Potential collaborators / response targets.**

- **Ionel Stetcu (LANL / nuclear-VQE lineage)** is the senior name across
  Stetcu et al. 2022 + the post-2022 deuteron-VQE thread. Stetcu's group is
  the most natural reviewer / cite-back partner for Paper 23.
- **Alessandro Roggero (Trento)** is the other senior name (Quantum Simulation
  of Nuclear Inelastic Scattering 2006.01369, and adjacent). Both Stetcu and
  Roggero will recognize the Minnesota potential + Moshinsky framing
  immediately.
- **Pooja Siwach + Calvin Johnson** lineage (sd/pf-shell VQE, Discover Quantum
  Science 2026) for the Slater-det-per-qubit comparison.
- **Calvin Johnson (San Diego State)** — long-running nuclear shell model
  classical-side authority + lifelong shell-code maintainer (BIGSTICK); would
  be the natural classical-side reviewer for the Minnesota / Moshinsky check
  on Paper 23.

**What each connection offers.**

- Stetcu / Roggero: in-discipline conversation partners; their immediate
  question will be "what does the composed nuclear-electronic encoding buy?"
  Paper 23 §VI (26 q, 614 Pauli with hyperfine validation against the 21cm
  gap) is the answer, and it's an angle they don't have.
- Slater-det-per-qubit (2510.02124): direct comparison head-to-head against
  Paper 23's encoding; gives Paper 23 an external numerical anchor.
- 2509.08948 deuteron-binding RG: same molecule, different encoding;
  comparison sharpens both papers.

---

### Theme 3 — Hopf-graph Ramanujan (Paper 29)

**Adjacent active programs.**

- **Matsuura–Ohta Kazakov–Migdal-on-graphs lineage** — arXiv:2204.06424
  (JHEP09(2022)178), arXiv:2208.14032 (PTEP 2022, Wilson loops + Ihara),
  arXiv:2403.07385 (PTEP 2024 Aug, "Phases and Duality in the Fundamental
  Kazakov–Migdal Model on the Graph"). Already known to CLAUDE.md §RH Sprint
  2 from the 2204.06424 hit; the 2024 paper extends to the *fundamental*
  Kazakov–Migdal model and explicitly discusses duality, which is a structural
  near-neighbor of Paper 29's Hopf-U(1) block decomposition.
- **Marcus–Spielman–Srivastava + extensions** — original 2013 interlacing
  families I (Annals of Math 2015), and arXiv:2108.02534 (Gribinski–Marcus,
  2021) extending to biregular bipartite Ramanujan of every degree and every
  size. The MSS thread is canonical and Paper 29's §6.1 Alon–Boppana crossing
  result lives in this neighborhood.
- **Quantum-graph spectral zeta + GUE lineage** — see e.g.
  arXiv:1810.08664 (spectral statistics of quantum circulant graphs, intermediate
  statistics) and the older Bolte–Harrison thread on quantum-graph zeta functions
  (arXiv:0710.4063 and follow-ups; "Zeta functions of quantum graphs" PMID
  ResearchGate 45884110). Paper 29's RH-M finding (CV ≈ 0.35–0.40 GUE-like on
  D(s) zeros) is unusual in this literature and would be a positive citation
  target.
- **Yakaboylu 2024 (JPhysA)** half-line Hilbert–Pólya candidate, already in
  CLAUDE.md §Sprint 2 RH-I. Direct adjacency, not duplication.
- **Smilansky / Berkolaiko quantum-graph spectral zeta program** — pre-2020
  canonical; if Paper 29 is to be sent for math-reviewer eyes the natural senior
  candidates are in this lineage.
- **Ramanujan complexes for quantum LDPC codes** — arXiv:2004.07935 (Evra–
  Kaufman–Zémor, decodable QLDPC codes beyond √n via high-dim expanders) and
  arXiv:2504.15087 (explicit lossless vertex expanders, 2025). Tangentially
  adjacent: Paper 29 is single-graph Ramanujan, not the QLDPC code-distance
  question; this is "neighboring discipline, won't be a cite" rather than
  "build a bridge."

**Specific papers to cite or build bridges to.**

- arXiv:2204.06424 + arXiv:2208.14032 (Matsuura–Ohta KM-Ihara): already named
  in Paper 29 §6 / §RH-I; the 2024 follow-up arXiv:2403.07385 should be added
  to Paper 29's bibliography on the next pass — it sharpens the duality reading
  that Paper 29's Hopf-U(1) block decomposition shadow-matches.
- arXiv:2108.02534 (Gribinski–Marcus 2021, biregular Ramanujan of all sizes):
  Paper 29's RH-D (finite-size statement) lives in the same "Ramanujan past
  the asymptotic regime is hard" neighborhood. Citation hardening.
- arXiv:1810.08664 (Maioli–Smilansky quantum-circulant graph intermediate
  statistics): adjacent to Paper 29's RH-M GUE finding, and probably the
  natural comparison background for the spectral CV ≈ 0.35–0.40 reading.

**Potential collaborators / response targets.**

- **So Matsuura + Kazutoshi Ohta (Japan, the JHEP09(2022)178 + PTEP 2022
  + PTEP 2024 chain)** are the closest active program — explicitly graph-
  Ihara-Wilson-loop in QFT, with a 2024 paper that extends to duality. Most
  natural conversation partners for Paper 29 + the Paper 25/30/41 Wilson
  arc combined (cross-cuts with Theme 4).
- **Adam Marcus (Princeton / EPFL)** — Ramanujan-construction lineage senior;
  the natural place to send Paper 29 for "is the finite-size statement
  Alon–Boppana crossing at V ~ 30–60 new?" feedback.
- **Uri Smilansky / Gregory Berkolaiko (or Brian Winn)** — quantum-graph
  spectral-statistics program; the natural reviewer for the RH-M GUE finding
  if PI wants out-of-ecosystem eyes.

**What each connection offers.**

- Matsuura–Ohta: very close formal program; their 2024 fundamental KM result
  is exactly the kind of thing that would be enriched by Paper 25/30's
  observation that L₁ = B^T B is a discrete-Hodge kinetic term + their Hopf-U(1)
  block decomposition that Paper 29 §5.3 validates.
- Marcus: probably the only person who will care whether the finite-size
  crossing at V ~ 30–60 is publishable as a Ramanujan-graph-theory statement
  in its own right.
- Smilansky / Berkolaiko: quantum-graph framing is the natural disciplinary
  home for the RH-M GUE finding even though GeoVac's framing is spectral-triple
  rather than quantum-graph.

---

### Theme 4 — SU(2)/SU(3) Wilson on physical graphs (Papers 25, 30, 41)

**Adjacent active programs.**

- **Garofalo–Hartung–Jakobs–Jansen–Ostmeyer–Rolfes–Romiti–Urbach SU(2)-on-S³
  partitioning lineage** — arXiv:2311.15926 (PoS LATTICE 2023; Nov 2023) and
  the follow-up arXiv:2503.03397 (March 2025, "Dynamics in Hamiltonian Lattice
  Gauge Theory: Approaching the Continuum Limit with Partitionings of SU(2)").
  This is the most direct near-neighbor I found: they partition the S³ ≅ SU(2)
  group manifold into discrete vertices, the link operators are unitary and
  diagonal, the canonical momenta are finite-difference operators approximating
  Lie derivatives, and they implement the standard Wilson Hamiltonian. **This is
  structurally Paper 30 done from a different motivation (digitization for
  quantum computing rather than spectral-graph theory).** Same headline object,
  different community. Should be cited prominently in Paper 30 if it's not
  already.
- **Pérez-Sánchez Bratteli networks** — arXiv:2401.03705 (Jan 2024, rev Feb
  2025) — Wilsonian Yang–Mills lattice gauge as continuum limit of a quiver
  spectral triple, with a Higgs field emerging from self-loops. This is the
  direct line to Marcolli–vS that CLAUDE.md WH1 already calls out, but Paper
  30 / Paper 25 / Paper 41 should also cite it to anchor the
  graph-Wilson-gauge-as-spectral-triple side.
- **Pérez-Sánchez comment on Marcolli–vS** — arXiv:2508.17338 (Aug 2025) —
  clarifies that the continuum limit of MvS is Yang–Mills *without* Higgs.
  Paper 30's Marcolli–vS lineage citation should pair MvS 2014 with this
  comment.
- **Finite-group Cayley-graph gauge theory** — arXiv:2503.17301 (March 2025,
  "Finite group gauge theory on graphs and gravity-like modes"; ScienceDirect
  pub 2025). Generalizes Cayley-graph differential calculus gauge story on
  arbitrary finite group G with gauge-network structure; discusses
  Lorentzian-version well-behaved + bulk/boundary effects on finite chains.
  This is the structural neighbor on the *non-Lie-group* side and is interesting
  as a comparison case (GeoVac is on SU(2), 2503.17301 is on finite G).
- **Sourav Chatterjee Wilson-loop-expectations program** — arXiv:2001.05627
  (CMP 2020, "Wilson loop expectations in lattice gauge theories with finite
  gauge groups") and arXiv:2001.07453 (2020 Wilson loops in finite Abelian
  lattice gauge theories). The probability-theory side of finite-group Wilson
  loops; Paper 25/30/41 are constructive/algebraic but the Chatterjee program
  gives a sharp asymptotic/concentration-of-measure framing.
- **Orbifold-lattice approach** — arXiv:2401.05731 (JHEP05(2024)234, "Toward
  QCD on quantum computer: orbifold lattice approach"). Quantum-computing
  motivated, but uses orbifold structures on the lattice that are adjacent to
  Paper 25's S³ → S² Hopf quotient.

**Specific papers to cite or build bridges to.**

- arXiv:2311.15926 + arXiv:2503.03397 (Garofalo et al., S³ partitioning SU(2)
  Wilson) — **Paper 30 should cite both prominently in §1 (positioning) and
  §6 (comparison).** They are doing the same SU(2)-Hamiltonian-on-S³ Wilson
  construction with a completely different toolkit (sphere partitioning
  digitization for quantum computing) and getting compatible structure. This
  is the strongest single bridge candidate in this scout.
- arXiv:2401.03705 + arXiv:2508.17338 (Pérez-Sánchez Bratteli networks +
  comment on MvS) — Paper 25 and Paper 30 should both cite to anchor the
  Marcolli–vS-without-Higgs reading already named in CLAUDE.md WH1.
- arXiv:2503.17301 (finite-group gauge on graphs, 2025) — Paper 30 / Paper 41
  open-questions section; useful comparison case.
- arXiv:2001.05627 (Chatterjee Wilson loops finite gauge): probably a Paper
  29 + Paper 30 joint cite under "rigorous Wilson-loop framings."

**Potential collaborators / response targets.**

- **Carsten Urbach (Bonn) and Karl Jansen (DESY)** lead the
  Garofalo et al. SU(2)-on-S³ partitioning program. They are the closest
  external program working on what is essentially Paper 30's headline
  object, and the natural conversation partners. The mismatch in
  motivation (they want quantum-computer-implementable digitization; GeoVac
  wants the spectral-graph object) is the angle to pitch.
- **Carlos Pérez-Sánchez (Heidelberg)** is the right senior NCG / quiver-
  spectral-action partner. His 2024 + 2025 papers are the direct math-side
  upgrade path for the Paper 25 → Paper 30 → Paper 32 lineage.
- **Sourav Chatterjee (Stanford)** for the rigorous-finite-group-gauge-
  loops side.
- **Walter van Suijlekom (Nijmegen)** is the canonical senior NCG reviewer
  here and via Marcolli–vS the literature root. Already implicit in
  CLAUDE.md WH1; explicit ping for the Paper 30 release would make sense.

**What each connection offers.**

- Urbach / Jansen: another community doing the same construction. The
  comparison is information for both sides and is the most concrete
  collaboration angle in this scout.
- Pérez-Sánchez: math-side upgrade path — Paper 30 / 32 / 41 could be
  re-read through the Bratteli / Yang-Mills-without-Higgs lens, and this
  is a natural pre-Zenodo math review.
- van Suijlekom: senior-NCG endorsement; Paper 30 is structurally inside
  his program's territory.
- Chatterjee: probability-theory framing for the Wilson-loop asymptotics
  Paper 30 §6 reports.

---

## Cross-cutting observations

- **Marcolli–vS 2014 + Pérez-Sánchez 2024/2025 is the single piece of
  external math literature that touches three GeoVac themes at once** —
  Theme 1 (the Bratteli / AF inductive-limit lift for Paper 24's Hardy
  sector), Theme 2 (the quiver / spectral-action reading of the composed
  nuclear-electronic POC), and Theme 4 (the YM-without-Higgs continuum
  limit). This was already named in CLAUDE.md WH1 in the Riemannian
  context; the explicit reading is that it should be cited in Papers
  23, 24, 25, 30, 41 not just Papers 25, 30, 32.
- **The Matsuura–Ohta KM-Ihara line and the Garofalo et al. SU(2)-on-S³
  line are doing GeoVac-shaped work in two different communities**
  (math-physics QFT-on-graphs vs lattice-QCD-for-quantum-computers). Paper
  30 + Paper 29 sit at the intersection of both, and both communities
  are unaware of each other in the search results.
- **The "first GUE signature" finding in Paper 29 RH-M is genuinely
  unusual** in the quantum-graph + Ihara literature I surveyed; that
  literature is dominated by Berry–Tabor intermediate statistics on
  Ihara-side and bulk-randomness GUE only via Hashimoto-type
  non-backtracking matrices. CLAUDE.md §RH-M's caveat (zeros not on a
  single critical line) is well-taken, but the GUE-on-spectral-zeta
  result is more novel than the intermediate-statistics framing might
  suggest and is the right thing to lead with for Paper 29's external
  pitch.
- **Pickrell's "potentially relevant to chiral models in 2d QFT"
  remark on the d = 2 spherical HO is exactly the AdS_3/CFT_2-adjacent
  angle that Paper 50's S^5 extension already touches via the F-theorem.**
  There may be a Paper-24-meets-Paper-50 thread here, mediated by the
  Hardy-sector ground state.

---

## Top 5 specific actionable connections

Ranked by combined "specific to GeoVac" × "tractable next step":

1. **Cite arXiv:2311.15926 + arXiv:2503.03397 (Garofalo–Hartung–Jansen et
   al. SU(2)-on-S³ partitioning) in Paper 30 §1 and §6.** This is the
   single closest external near-neighbor in the entire scout — the same
   Wilson Hamiltonian on the same S³ ≅ SU(2) group manifold from a
   different (quantum-computing-digitization) motivation. Both
   communities can use the comparison. Concrete next step: Urbach
   (Bonn) / Jansen (DESY) email pointing at the Paper 30 Zenodo.

2. **Cite arXiv:2503.23549 (Higgs–Pickrell) in Paper 24 §III and ping
   Pickrell (Arizona).** The published-March-2025 d = SO(d) HO with the
   (1+r²)^{−ω} ground state and the holomorphic-line-bundle structure is
   the closest published-math analog of Paper 24's S^5 Hardy-sector
   construction. Pickrell's open d ≥ 3 question is exactly Paper 24's
   d = 5 case. Probably the highest-density single-paper match across
   the four themes.

3. **Add arXiv:2401.03705 (Pérez-Sánchez Bratteli networks) +
   arXiv:2508.17338 (Pérez-Sánchez comment on Marcolli–vS) to Papers 25,
   30, 41 bibliographies.** CLAUDE.md WH1 already names these; the next
   step is mechanical bibliography update across Papers 25/30/41/32 so
   that the Marcolli–vS-without-Higgs reading is explicit on the gauge
   side of the GeoVac corpus, not just in the WH1 register.

4. **Cite arXiv:2510.02124 (Slater-det-per-qubit, 2025) and
   arXiv:2509.08948 (RG-deuteron-VQE, 2025) in Paper 23 §VI.** Both are
   from 2025 and both are direct comparison-target competitors for
   Paper 23's deuteron Pauli count and the composed
   nuclear-electronic POC. The Paper 23 release is meaningfully
   strengthened by an explicit head-to-head table against these
   competitors. Stetcu (LANL) and/or Calvin Johnson (San Diego State)
   are the natural reviewers.

5. **Add arXiv:2403.07385 (Matsuura–Ohta 2024 PTEP "Phases and Duality")
   to Paper 29's bibliography and check whether the duality reading
   they expose maps onto the Hopf-U(1) block decomposition in Paper 29
   §5.3.** This is the natural follow-up to Sprint 2 RH-I's
   identification of arXiv:2204.06424. If the dualities match, it's a
   joint cite + potential conversation with Matsuura/Ohta (Keio
   University / Meiji Gakuin University).

---

## Limitations

- I did not download any of the cited PDFs. Section-level claims (e.g. "Pickrell
  explicitly discusses a holomorphic Hermitian line bundle") came from search
  snippets and abstracts; before citing in a GeoVac paper these need a one-pass
  read against the actual PDF.
- The senior-author attributions are inferred from author order, affiliation
  surface in search results, and a small amount of prior knowledge; PI may
  want to verify Stetcu / Roggero / Johnson / Urbach / Jansen / Pickrell /
  Hall / Matsuura / Ohta / Marcus / Smilansky / van Suijlekom / Pérez-Sánchez
  attributions against the actual papers before reaching out.
- I made no attempt to find pre-2015 canonical references that are likely
  already in the relevant papers' bibliographies (Brian Hall's foundational
  Segal–Bargmann work, original Marcolli–vS 2014, Connes–Chamseddine 1996/97,
  etc.) — task scope was post-2015.
- I did not search arXiv listings directly for May–June 2026, so the freshest
  preprints are likely under-represented.
- I did not check whether any of the citations recommended above are
  *already* in the relevant paper bibliographies. The recommendation to add
  them assumes they are not; PI should verify before mechanical bibliography
  update.
- Paper 24's specific S^5 Hardy claim is sharper than anything I see in the
  Hall / Pickrell literature, but I can't certify "no near-duplicate exists"
  — that would require a dedicated math.OA search beyond this scout's scope.
