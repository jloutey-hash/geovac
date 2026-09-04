# GeoVac Project Guidelines

## 1. Project Identity

**Name:** GeoVac (The Geometric Vacuum)
**Version:** v5.6.0 (September 3, 2026)
**Mission:** Spectral graph theory approach to computational quantum chemistry. The discrete graph Laplacian is a dimensionless, scale-invariant topology (unit S3) that is mathematically equivalent to the Schrodinger equation via Fock's 1935 conformal projection. This equivalence is exploited computationally to replace expensive continuous integration with O(N) sparse matrix eigenvalue problems.

**Mission statement (adopted 2026-08-29, PI direction):** GeoVac charts the forced/free boundary of quantum physics. For every structure -- quantum number, selection rule, degeneracy, sparsity pattern, convergence rate, physical constant -- the program renders one of three verdicts: **FORCED** (derived from the packing construction, exactly, with a frozen falsifier), **FREE** (an exchange constant, with its projection chain named and its minimal transcendental content classified), or **WALL** (a proven obstruction with the mechanism pinned). The deliverable is the atlas of that boundary. Under this statement the corpus is one program: the chemistry/QC arc surveys how far the forced side reaches computationally; the 40+ documented negatives are the boundary itself, measured; the precision program (SS1.8) is the atlas's experimental interface; the periods/transcendence work (Papers 18/34/54-59) is the coordinate system for the free side; and the QA apparatus is what makes the atlas trustworthy. *Scope note:* this is the internal research mission (the register of SS1.7); papers remain under the SS1.5 rhetoric rule -- "forced" is atlas vocabulary, not ontology language for publication.

**Authoritative source rule:** The papers in `papers/group1_operator_algebras/`, `papers/group2_quantum_chemistry/`, `papers/group3_foundations/`, `papers/group4_quantum_computing/`, `papers/group5_qed_gauge/`, `papers/group6_precision_observations/`, and `papers/synthesis/` are the authoritative source for all physics. If any documentation (README, CHANGELOG, code comments) conflicts with the papers, the papers win. Flag the conflict to the user rather than silently resolving it. (Papers were reorganized from the previous `core/`, `methods/`, `applications/`, `synthesis/`, `standalone/`, `observations/`, `conjectures/` layout into six audience-targeted groups on 2026-05-22.)

**Project context:** GeoVac is an independent research project with no institutional affiliation, developed using an AI-augmented agentic workflow. The principal investigator (Josh Loutey — `J.~Loutey` on the papers) provides scientific direction and quality control; implementation and documentation drafting are performed collaboratively with LLMs (Anthropic Claude). The primary dissemination channel is GitHub + Zenodo (DOI-stamped releases). The papers across the six audience-targeted group folders (see §6) are written to academic standards but are not submitted to traditional journals. The project's viability case rests on its corpus of verified structural results — discretization- and encoding-structure theorems, honest negative results, and the benchmarked computational artifact (pip-installable, 37-system library) that demonstrates them. The tool is a research instrument, not a production-chemistry replacement (see the benchmarking rule in §1.5 and `docs/claims_register.md`). Do not suggest formatting papers for specific journals or pursuing traditional peer review unless the user asks. (Viability-case sentence reworded 2026-06-10 per PI direction, accessibility plan Phase 3.)

---

## 1.5. Positioning & Framing

GeoVac is a discretization framework that exploits the natural geometry of separable quantum systems to produce sparse graph Hamiltonians. It is **not** a new foundation for quantum mechanics — it works because the graph Laplacian converges to the known continuous Laplace-Beltrami operators in the continuum limit (Paper 7).

**Rhetoric rule:** The GeoVac framework is demonstrably conformally equivalent to standard continuous quantum mechanics at every level where it has been tested (Paper 7 provides 18 symbolic proofs for the S³ case; subsequent papers verify operator convergence for each natural geometry). The papers should present this conformal equivalence as the primary result. The discrete graph topology and the continuous Laplace–Beltrami operators are dual descriptions connected by proven conformal maps — the mathematics supports both readings, and the interpretation of which is more fundamental is left as a choice for the reader. Avoid language that asserts ontological priority of either description. Lead with concrete computational results (sparsity, scaling, accuracy, structural insight) rather than interpretive claims about the nature of quantum mechanics.

**Lead with concrete advantages:** O(V) sparsity, angular momentum selection rules baked into the basis, zero-parameter construction from nuclear charges and geometry alone, efficient qubit encodings (Paper 14). These structural properties are the framework's actual selling points — not philosophical claims about the nature of quantum mechanics.

**The π-free graph principle:** The GeoVac graph is π-free: all eigenvalues, degeneracies, and coupling coefficients are integers or rationals. Transcendental numbers (π, exponential integrals, spectral zeta values) enter exclusively when projecting the graph onto continuous manifolds for comparison with experiment or for computational convenience. This observation (Paper 18) motivates the design principle: stay on the graph whenever possible, and when projection is necessary, identify the minimal transcendental content (the exchange constant) required. The exchange constant taxonomy (intrinsic, calibration, embedding, flow) classifies projections by what determines them and predicts the computational cost of each departure from the graph.

**Quantum simulation positioning:** The classical solver investigation (v2.0.6-23) characterized what the natural-geometry hierarchy can and cannot achieve for ground-state PES; the primary computational value proposition is quantum simulation. The composed architecture produces qubit Hamiltonians whose Pauli count is EXACTLY linear in Q across molecules at fixed basis — N_Pauli = 27.90 x Q universal (main-group; 30.03 d-block), vs O(Q^3.9-4.3) Gaussian system-size scaling — with 54x-317x fewer Pauli terms at equal qubits (H2O, n_max=1-3), structural Gaunt sparsity compatible with all downstream optimizations (tapering, grouping, tensor factorization), and a 37-system library (35 composed molecules + He + H2) across three periodic-table rows. The within-molecule basis exponent is 3.17 (universal across all systems, small-basis fit; local slope rises toward ~3.8 by n_max=4 — the linearity in Q, not the basis exponent, is the load-bearing advantage). Position for VQE/NISQ (Pauli count, QWC groups) rather than fault-tolerant QPE; composed LiH is near parity with minimal-basis STO-3G on raw counts (838 Pauli @ Q30 vs 907 raw; 3.0x MORE than the 276 reduced; 1-norm 1.007x, 34.5 vs 34.3 Ha) — the decisive wins are the exact Q-linearity and the 76x cc-pVDZ reduction. (All numbers under the exact global-M_L rule, corrected 2026-08-29 — the prior 11.10xQ / O(Q^2.5) / 51x-1712x / 334 / 190x line was measured under the retired pair-diagonal sign error; see debug/sprint_eri_evaluator_defects_memo.md.) Full track-by-track chronicle: `docs/development_frontier_archive.md` + CHANGELOG.

**General composed builder foundation:** Paper 16 (Chemical Periodicity as S_N Representation Theory) provides the group-theoretic foundation for the general composed builder's atomic classification. The structure types (A/B/C/D/E), the universal angular quantum number ν = N−2, and the recursive core-valence decomposition are all derived from S_N representation theory and implemented in `atomic_classifier.py`.

**PK limitation and research path:** The composed framework's PK pseudopotential is the accuracy bottleneck (5.3-26% R_eq error); six modification attempts failed (Section 3). The balanced coupled Hamiltonian (Papers 19/20) is the PK-free alternative: cross-center V_ne via multipole expansion with exact Gaunt termination; LiH 4e FCI reaches 0.20% energy at n_max=3 with structural R_eq drift (~8.8%) — energy converges excellently, geometry drifts; optimal for single-point quantum simulation at fixed geometries. Resource anchor: composed LiH 838 Pauli vs STO-3G 907 vs cc-pVDZ 63,519 (exact rule, 2026-08-29). Chronicle: `docs/development_frontier_archive.md`.

**Benchmarking rule:** When comparing to other methods, always use the strongest available baseline (cc-pVTZ or better for atoms, explicitly correlated methods for molecules), not just STO-3G. If the comparison is unfavorable, say so honestly and identify what the framework offers instead (sparsity, scaling, structural insight). Position the framework as a computationally principled alternative, not as a replacement for production quantum chemistry.

---

## 1.6. Project Phase

**Phase 1 (v0.9.x-v1.x): Foundation.** Graph Laplacian equivalence proof, natural geometry hierarchy, LCAO diagnostic arc, bond sphere theory, 18 symbolic proofs (Paper 7).

**Phase 2 (v2.0.0-v2.0.23): Classical solver investigation.** Systematic exploration of all solver architectures (adiabatic, coupled-channel, 2D variational) across all levels (2-5, 4N). 30+ completed tracks, 40+ documented negative results. Key outcomes: H2 at 96.0% D_e (Level 4), LiH R_eq at 5.3% (Level 5 composed), full N-electron equilibrium without PK (Level 4N). Investigation complete: all solver x PK x basis combinations exhausted, structural accuracy ceilings characterized, exchange constant taxonomy established (Paper 18).

**Phase 3 (v2.0.24-v2.7.0): Quantum simulation.** The composed architecture's structural sparsity (exact Q-linearity of the Pauli count across molecules at fixed basis, N_Pauli = 27.90 x Q; block-diagonal ERIs) is the framework's primary computational advantage. (Corrected 2026-08-29: the retired pair-diagonal rule gave O(Q^2.5) here and named the *within-molecule* basis exponent as the advantage; under the exact rule that exponent is 2.82 on the two-point fit and ~3.8 locally by n_max=4, and it is the across-molecule linearity -- not the basis exponent -- that is load-bearing.) Classical solver results serve as validation benchmarks for the quantum Hamiltonians. Next steps: quantum resource estimation, hardware-aware circuit compilation, experimental collaboration.

**Phase 4 (v2.7.0+): Nuclear extension and framework delineation.** Extension of the GeoVac framework to nuclear shell model Hamiltonians using the hyperspherical (HO + spin-orbit) basis, together with a precise delineation of which parts of the framework transfer to non-Coulomb systems. Key outcomes: (1) the potential-independent angular sparsity theorem (Paper 22) — ERI density depends only on l_max, not on V(r), verified at l_max=3 as 6.06% (global-M_L headline; 1.44% was the retired pair-diagonal rule, corrected 2026-08-29); (2) the deuteron and He-4 qubit Hamiltonians (Paper 23, Tracks NE/NF) using the Minnesota NN potential, Moshinsky-Talmi brackets, and a two-species tensor-product JW encoding; (3) the Fock projection rigidity theorem (Track NH) — the S^3 conformal projection is unique to -Z/r and does NOT transfer to HO or Woods-Saxon; (4) the composed nuclear-electronic deuterium PoC (Track NI) demonstrating the block architecture at 26 qubits with hyperfine validation; (5) the Bargmann-Segal lattice (Paper 24, Track NK) — discrete graph encoding of the 3D HO on the holomorphic sector of S^5, bit-exactly π-free in exact rational arithmetic at every N_max (verified at N_max=5: 56 nodes, 165 edges, zero irrationals); (6) the universal vs Coulomb-specific partition as the Phase 4 conceptual result, completed by the Coulomb/HO asymmetry analysis in Paper 24 (calibration π is structurally Coulomb-specific — tied to second-order Riemannian operators with nonlinear projections, has no analog for the HO's first-order complex-analytic projection). Track NJ memo definitively shelved the nuclear → alpha connection (originally "shelve," upgraded after Sprint 2/Paper 24 to "definitively shelved" based on the structural impossibility, not just empirical absence). HO rigidity theorem (Theorem 3 of Paper 24) is the structural dual of the Fock rigidity theorem.

---

## 1.7. Working Hypotheses (Internal Register)

This section is a **bold-claim register**, distinct from the rhetoric of the papers. Papers remain cautious under §1.5 (dual-description framing, no ontological priority); this register is what sub-agents and the PM may *reason from* during synthesis work. Nothing here appears in papers unless promoted after its falsifier clears. Full status chronicles (the sprint-by-sprint evidence trail, April–June 2026) live in `docs/wh_register_history.md`; this section holds only claim, falsifier, and current status.

**Governance:**
- PM may update a WH's "Status" line based on sprint evidence — REPLACE the line and move superseded text to `docs/wh_register_history.md`; never append (§13.11 rule 9).
- Adding, retiring, or promoting a WH to paper-level claim requires explicit PI direction.
- Retired WHs move to the history doc with rationale; never silently deleted.
- No WH is a license to bypass the rhetoric rule in papers or the verification gates in §13.4.

---

**WH1 — GeoVac is an almost-commutative spectral triple.** A = functions on the Fock-projected S³ graph; H = scalar/spinor state space; D = Camporesi–Higuchi Dirac; non-abelian gauge structure enters as inner derivations of the almost-commutative extension A ⊗ M_n(ℂ) (Marcolli–van Suijlekom lineage; Papers 25/30/32).
*Falsifier:* a GeoVac observable demonstrably inconsistent with any spectral-action expansion, or a violated structural axiom (order-one, reality).
*Status:* **PROVEN — unconditional (2026-06-10; constant re-priced 2026-09-03).** Paper 38: the discrete truncations converge to the round-S³ spectral triple in van Suijlekom's state-space GH distance at rate (c + o(1))·log n_max/n_max, on the truthful CH substrate (translation-seminorm metrization; frozen falsifier `tests/test_p38_action_seminorm.py`). The constant is c = 4/π in the rule Paper 40 declares (Cas(ad) = h∨), under which the SU(2) geodesic distance is the rotation angle — so P38 is P40's rank-1 case *in one consistently applied convention* (settled 2026-09-03; `tests/test_p38_metric_convention.py`). That rule is NOT the field-standard one: Kac's basic form (θ|θ) = 2 gives Cas(ad) = 2h∨, radius √2 and constant 2√2/π, and the unit sphere gives 2/π. So 4/π is convention-dependent, not canonical — an earlier version of this line said otherwise. On the unit S³ the same bound reads 2/π. The "SU(2) = twice the circle" reading is the metric scale. Lemma L5's height bound was REFUTED 2026-09-03 (FULL run #4: height_B ≡ 1 by a finite-band witness, so the bound fails for every n_max ≥ 6) — the unconditional theorem uses only reach-type estimates and is untouched, so the keystone stands and one proof path falls. The Lorentzian extension (Papers 45–49) is DESCOPED (P45 annihilation theorem); repair path = Toeplitz temporal compressions (see WH7).

**WH2 — Paper 18 is the Seeley-DeWitt + ζ-invariant decomposition of this spectral triple.** The transcendental taxonomy is the structured output of spectral-action geometry, organized by operator order × bundle type.
*Falsifier:* a transcendental in a GeoVac observable that cannot be placed in the grid.
*Status:* three of four axis-quadrants filled; mixed-Tate sharpening POSITIVE (2026-06-03) — M2 on S³ sits in the pure-Tate sub-ring ⊕_k π^{2k}·ℚ (Fathizadeh–Marcolli inherited). See `debug/sprint_mixed_tate_test_memo.md`.

**WH3 — The lattice exists a priori; match to physics is evidence, not derivation.** The packing construction is independent of known physics; its persistent match (Fock S³, nuclear magic numbers, Dirac fine structure, Pauli sparsity) is evidence physics is hosted by a discrete spectral triple, not that the lattice was reverse-engineered.
*Falsifier:* the packing construction fails to force the (n, l, m, s) structure or the n²−1 spectrum; or a match turns out to depend on a hidden physics-informed parameter.
*Status:* origin-story framing permits strong ontological claim internally; papers stay under §1.5.

**WH4 (deflated 2026-05-07) — The four-way S³ unity is one Fock-projection statement plus three forced consequences.** Bertrand + SO(4) force S³; the Hopf base, the CH spinor bundle, and SU(2) Wilson all follow from S³ = SU(2) parallelizability / maximal-torus structure. Does NOT extend to inner-factor selection (Yukawas, generations remain unforced).
*Falsifier:* a construction forcing one of the four roles onto a different manifold; or a published argument decoupling a consequence from the Fock input.
*Status:* deflated to a single-input forcing statement; outer structural unity essentially closed (Sprint TS-D + Paper 38).

**WH5 — α is a projection constant, not a derivable number.** K = π(B + F − Δ) composes three structurally independent spectral objects (finite Casimir trace; Fock Dirichlet ζ(2); Dirac boundary count 1/40); the right open question is why the sum equals α⁻¹, not how to derive each piece.
*Falsifier:* a spectral-triple construction deriving K as a single coefficient of a well-defined functional.
*Status:* TWELVE mechanisms eliminated (Phases 4B–4I + Sprint A + Sprint K-CC, including the T9 algebraic obstruction: no single CC heat-kernel expansion contains B, F, Δ as terms). Standing reading: three-regime projection coincidence. Paper 2 stays in Observations; combination rule labeled an Observation — not a conjecture or derivation (§13.5 hard prohibition; conjecture→observation downgrade 2026-06-14, PI direction — "conjecture" judged to carry unearned confidence that a derivation exists).

**WH6 — GeoVac's RH-adjacent object is the Dirac spectral zeta D(s), not classical ζ.** Internal GUE-like zero statistics (CV ≈ 0.35–0.40); the classical-RH bridge is closed by three independent walls (zeros not on one line; no spectral-triple-natural functional equation, 48 OoM; wrong Weyl class).
*Falsifier:* D(s) zeros on a single critical line at larger samples; or a natural functional equation closing the RH-O gap.
*Status:* paused (2026-04-18); if resumed, the target is D(s) and its spectral-action interpretation.

**WH7 — Time-discreteness is observer-compactification (registered 2026-06-10, PI direction).** The only temporal structure the framework can see metrically is compactified time, and compactified time is automatically discrete. Inputs: (i) Paper 35 — π enters exactly at temporal compactification (Matsubara 2πk/β); (ii) P45 annihilation theorem — ℝ-time is Lipschitz-invisible in the v1 architecture; (iii) Paper 47 — the three temporal carriers are spectrally indistinguishable. Reading: discreteness of time is supplied by the observer's compact integration window (KMS β = 2π, four-witness theorem, Connes–Rovelli thermal time) — time is the prototype free-side projection. Temporal restriction of the organizing observation below.
*Falsifier (primary):* the Toeplitz temporal-compression program — a metrically visible temporal algebra on a NON-compact carrier without compactification weakens WH7 to convention; a proven annihilation-type obstruction for the framework's whole non-compact temporal class forces it.
*Falsifier (secondary, from Paper 35):* a GeoVac observable containing π whose evaluation provably involves no temporal/spectral integration.
*Status:* REGISTERED (2026-06-10). Load-bearing leg is **(i)**: only compactification makes time discrete and injects π. The Toeplitz probes leave the primary falsifier leaning **weakens-to-convention** — time is metrically visible; the observer's window is what makes it discrete. The Lorentzian-propinquity chase is **structurally closed**: the truncated Bisognano–Wichmann boost is *compact* (integer modular spectrum), so every Lipschitz/boost seminorm is signature-blind and the Lorentzian leg is de-compactification = convention. Re-read 2026-08-21 as an **arithmetic (CM)** structure — the finite-cutoff modular circle is P56's ℚ(i) Hodge circle (μ₄, exact tests); the split-form "Wick = base change" mechanism was **REJECTED** (the scalar sector is compact but pure-Tate). Honest cap: input (iii) means discrete-vs-continuous may be empirically undecidable at every computable level; papers stay under §1.5. Full chronicle + superseded status text: `docs/wh_register_history.md`.

**WH8 — The Born measure is the exchange constant of the observation projection (registered 2026-07-01, PI direction).** The compact↔continuous junction accounts for the discreteness of outcome menus (records are bound/compact ⇒ spectra discrete) and the diagonality of effective observer states (compact-window KMS, four-witness theorem): the same projection that makes time discrete and injects 2π (WH7) injects |ψ|² at observation. The selection measure itself is skeleton-external, both directions now closed: *derivation*-direction closed structurally 2026-05-26 (Paper 34 §VIII Born entry — three routes all reduce to standard GNS/Gleason inheritance; sprint `ahha_born_rule_attempt`); *generation*-direction tested negative 2026-07-01 (the graph's counting measure = Born exactly on the flat locus and only there — exact-rational TV=1/3 on the generic 4-multiplet; TV=0.86 vs the He graph-native CI ground state, 90.8% on 1s²; the flat locus is dynamically unstable under the graph's own H, TV~t², exponent 2.00). Given the projection lattice the import is unique (Gleason, dim≥3) — the residual mystery is single outcomes, not the measure's form.
*Falsifier:* a skeleton-side rule (state-independent, built from the graph's discrete invariants — degeneracies, selection rules, rational couplings) reproducing Born statistics on non-flat physical states (breaks `tests/test_wh8_born_probe.py`); or a composition of more-primitive skeleton projections that carries the measure natively (would upgrade the import to derived and falsify the exchange-constant classification).
*Status:* REGISTERED (2026-07-01); Step-1 probe same day **NEGATIVE-as-expected** (`debug/sprint_wh8_born_probe_memo.md`; falsifier `tests/test_wh8_born_probe.py`; Paper 34 §VIII Born entry updated with the quantitative converse). The external-input Class-1 placement of the Born rule is upgraded classification→tested-negative. Honest cap: this fences, not solves, the definite-outcome problem — no interpretation on the market resolves single outcomes without added structure; papers stay under §1.5 (Paper 34 §III.28 + §VIII carry the tier-honest record).

---

**Organizing observation — discreteness is compactness (established).** The discrete, bit-exact skeleton is the closed/compact regime: compact groups have discrete spectra (Peter-Weyl), and this is the mathematical content of the S³ graph. The continuum is what remains when compactness is released. Calibration data is continuum (un-packed) data, which is why the skeleton is rational/forced and calibration is transcendental/free. Well-captured in Paper 18 §III (compactness thesis) and Paper 35 (temporal compactification injects π). The "second packing axiom" framing is retired; the 2026-05-30 confinement reframing is archived as an organizing reading (2026-05-31).

---
## 1.8. Multi-Observable Focal-Length Decomposition Program (RETIRED)

**Dormant.** Queued 2026-05-09 with five deliverables (Paper 34 §V.C.2–6);
**none was ever built** (verified 2026-09-01); successor = the v4.71.0
chemistry-error sprint. Verbatim: `docs/focal_length_program.md`. Reviving
it is a PI call.

---
## 2. Current Development Frontier

> Full sprint chronicles live in `CHANGELOG.md`. This section is a compact index. Sprint detail is in the memos linked below.

> **⚠ ORIENTATION — read before answering anything about polyatomics, accuracy, or which basis GeoVac uses.** Two facts a fresh session otherwise re-derives at the PI as a corrective lecture (PI direction 2026-08-14, after three derailed sessions): **(1) exact ≠ accurate.** The closed-form ERI arc (v4.77–79) bought *decidability*, not accuracy — the whole integral-exactness axis is worth ~0.003 Ha on H₂ and then flatlines, while basis size is worth 0.050 Ha and rising; v4.73.0 localized the chemistry defect to **100% max_n**, bit-invariant to angular and quadrature refinement. **(2) three bases are in play** — theory = Coulomb-Sturmian (shared p₀), `composed_qubit`/`two_center_eri` = **hydrogenic** (a = Z/n), `noci_engine` = STO shapes fitted by Gaussians (an *evaluator*, not a basis). Full map + the three-centre wall + "Paper 58 is not Avery and not Sturmian": `memory/polyatomic_state_of_play.md`.

> **Older entries (84 bullets, rounds 2–7) are in `docs/development_frontier_archive.md`.** §2 is the index; CHANGELOG.md is the chronicle.

- **Normalisation cross-check settled (2026-09-03, v5.5.0, PI direction):** Cas(ad) = h∨ on su(2) gives the rotation angle as the dual-Coxeter geodesic distance, so P38's moment is dual-Coxeter and 4/π is canonical; P38 restated on that sphere. The rule is the corpus's own, not the field-standard one (Kac: Cas(ad) = 2h∨), so 4/π is convention-dependent. The constant is the unit-sphere quotient Vol(S2)/Vol(S3) = 2/pi; M1's volume content stands, the "Hopf base" label is a misnomer. See carryforward Part K.
- **L5 crossing measured (2026-09-03, v5.4.4):** margin +0.053 at n_max = 6, −0.069 at 7 — the crossing is between 6 and 7; the papers now carry the measured table, not the extrapolated threshold.
- **L5 panel check corrected (2026-09-03, v5.4.3):** the guard added in v5.4.1 verified a small-cutoff coincidence — the panel height rises toward 1 while gamma falls, crossing near n_max = 6. Reframed as a measurement; panel-side quantity = open check. See carryforward J.4.
- **DELTA #3 (2026-09-03, v5.4.2):** DEFECTS — the v5.4.1 remediation was locus-by-locus, so the descoped readings survived at 12 loci; swept claim-wide. New: s/p splitting is a node-amplitude proxy (disconnected l-blocks), and the block spectrum is closed form with a proven O(n^-2) rate. See carryforward Part J.
- **/qa trunk FULL run #3 (2026-09-03, v5.4.1):** FAIL (55 MATERIAL, unseeded), remediated same day. Five re-pricings: graph→S³ convergence was κ·λ_max by construction; P38 constant 2/π on the unit S³ (4/π = rotation angle); Forced count 260 → 32 (representation bug); P32 Thm 1 scope; prop = 2 generic. See carryforward Part I.
- **/qa trunk FULL run #4 (2026-09-03, v5.6.0):** FAIL, remediated. Two gates examined NOTHING on trunk (C17 zero families, C21 zero annotations); s/p evidence leg = mod-3 artifact (withdrawn, but closed form gained); P38 L5 height bound FALSE (WH1 unaffected); three papers named operators they do not use. See trunk.carryforward.md Part M.
- **Trunk DELTA #4 (2026-09-03, v5.5.2):** DEFECTS. Two reviewers, same two LARGE: "4/pi is canonical" was my own over-claim (the corpus rule is twice Kac's), and the Hopf retraction missed ~10 loci. Three new guards could not fail. All fixed + C16-registered. See trunk.carryforward.md Part L.
- **v5.4.0 PI adjudications (2026-09-02):** P40 retitled (semisimple); finite triple prints sign triple (−,+,+), no KO label; Paper 7 convergence test added; P40 §L5 = state-space GH. See CHANGELOG v5.4.0.
- **QA gate: seeding opt-in (2026-09-02, PI direction):** blind calibration seeds are no longer default for DELTA/FULL runs; `/qa <target> seeded` invokes them. See .claude/commands/qa.md.
- **Trunk Part F + DELTA #1/#2 (2026-09-02, v5.3.1):** 30 rows closed; DELTA #1 9/9 DEFECTS (26), DELTA #2 8/9 DEFECTS (15), both remediated; PI items: P40 title, finite KO label. See trunk.carryforward.md Parts G/H.
- **Trunk run-#2 PI items (2026-09-02, v5.3.0):** circle Fejér 2/π (SU(2) = twice, Observation); P32 `prop:D_equiv` → Remark, KO-3 relabel; qa.md code-tier exception (gate change). See CHANGELOG v5.3.0.
- **/qa trunk FULL run #2 (2026-09-02, v5.2.7):** FAIL (19/21 seeds, 0/8 FP). P38's circle Fejér constant is 2/π (SU(2) 4/π stands); P32 `prop:D_equiv` unbacked; KO-dim label. Scope: `docs/qa/trunk.carryforward.md` Part F.
- **CLAUDE.md compaction (2026-09-01, v5.2.6):** 220 → 100 KB, verified lossless (27/27 conservation checks). §2 and §3 were 56% of the file; full text now in `docs/failed_approaches_ledger.md` + the frontier archive.
- **C11 could not fail (2026-09-01, v5.2.5):** proving a NEW gate criterion fires exposed that the existing one never could -- findings keyed papers-relative, predicate matches repo-relative, so a planted wrong title was detected, printed, and exited 0. Fixed + pinned. ...
- **Trunk C3 tier pass (2026-09-01, v5.2.5):** 257 inline tiers across 6 documents (0 before); run as an audit it found ~25 defects, one class dominating -- correct in the middle, overclaiming at both ends. +30 bibitems closed 30 unreachable prose ...
- **Re-cert sweep pre-flight (2026-09-01, v5.2.4):** 4 gates mis-scoped, C16 checked *nothing* for 5 of 11 targets; balanced λ was at LiH's bond length for every molecule; C21 widened corpus-wide -> 30 retired values fixed, all 11 targets PASS; Paper ...
- **Forced/free seam in the resource tables (2026-09-01, v5.2.4):** N_Pauli is geometry-invariant (FORCED, molecular analogue of Paper 22), N_QWC moves (greedy heuristic), λ moves (FREE). See tests/test_paper20_geometry_independence.py.
- **Test purpose resolved to zero unknowns (2026-08-31, v5.2.3):** PI rule -- no human sweep, the papers are the criterion. Docstring inference + falsifier/guard/negative tiers took unknown 14 -> 0, paper-backing 49% -> 55%; 3 superseded prolate tests archived. See debug/qa/test_suite_cost_memo.md.
- **C22 test-claim backing gate (2026-08-31, v5.2.2):** C13 run backwards -- every claim a test backs must still exist and not be retracted. 4 checks, all proven to fire; 3 debug/-importing paper tests baselined. See docs/qa/criteria.md C22.
- **Test-suite cost + purpose index (2026-08-31, v5.2.1):** full suite measured at 6.7 h (not 15 min); 27% of files over budget; xdist 2.77x adopted; reverse test->claim index built, 14 decay candidates. See debug/qa/test_suite_cost_memo.md.
- **Stage-4 cross-document ledger (2026-08-30, v5.2.0):** 13 ledger items closed; M9 was a 4-document pair-diagonal zombie; two 1-norm columns mixed conventions; angular gradient adds 78/78 L_z-violating entries. See debug/qa/stage4_ledger_notes.md.
- **ALL 11 QA TARGETS OWE RE-CERTIFICATION (measured 2026-08-31).** Every `docs/qa/*.done.md` says CERTIFIED; every one is stale (trunk, group1-6, synthesis, papers 58/59/60). Do NOT cite any of them as present-tense status. Compute, never remember: `python debug/qa/check_cert_staleness.py --detail`. Each record now carries a generated ...
- **Numeric registry + C21 gate (2026-08-30):** the corpus's numeric dependency graph is now a declared, maintained artifact (`debug/qa/numeric_registry.py`) with a gate that recomputes derivations and blocks retired values. Built after a census showed 352 numerals live at multiple loci / ...
- **group6 re-cert DELTA (2026-08-29, v5.1.11):** DEFECTS (11/11 seeds caught, 6/6 dims calibrated) -> same-day remediation; tautological D-pin rewritten; 13 one-locus-away survivors closed. FULL cert unlocked. See debug/qa/group6_full_run_2026_08_28_notes.md.
- **/qa group6 FULL cert #2 + remediation (2026-08-29, v5.1.10):** FAIL (9/9 dims calibrated) -> ledger fully remediated; Friar 3pi/8; Minnesota -0.55/+17.3; K-residual = P2 cubic root. See debug/qa/group6_full_run_2026_08_28_notes.md.
- **ERI evaluator defects -> exact rule everywhere (2026-08-29):** `q = mc - ma` sign error (= the whole content of "convention A"/CF-1) found in SEVEN modules; two evaluators bracketed truth (65 < 107 < 265), arbitrated by quadrature-from-definition. PI call: exact rule ...
- **UNFROZEN (2026-08-13, PI direction):** close-out freeze lifted; work/sparsity-boundary + tags v4.77-79 pushed. **Corrected 2026-08-26:** `main` is NOT at the v4.76.0 freeze point — it sits at `51f2eef` (2026-08-15), 58 commits past v4.76.0, and contains v4.77.0/v4.79.0. No mechanical freeze ...
- **Close-out freeze (2026-07-09, v4.76.0):** 59 per-paper Zenodo DOIs published; 443 dangling debug/ refs resurrected; dispositions PI-approved. See debug/sprint_closeout_distribution_memo.md + docs/project_closeout_plan.md.

**Best results by system type:**

| System | Result | Method | Paper |
|:-------|:-------|:-------|:-----:|
| He (atom) | 0.004% (cusp); 0.022% (raw) | 2D variational, cusp l_max=4 / raw l_max=7 | 13 |
| He (graph-native CI) | 0.19% | Zero-parameter, exact algebraic integrals, n_max=7 | 13 |
| H⁻ | Bound, over-binds 21% | Graph-native CI, Z_c≈1.84 boundary | 13 |
| PsH | 4.1% | Level 3, sign-flipped charge | 13 |
| H₂⁺ | 0.0002% | Spectral Laguerre | 11 |
| H₂ | 96.0% D_e | Level 4, l_max=6, 61 channels | 15 |
| LiH | R_eq 5.3% | Composed, l-dependent PK, l_max=2 | 17 |
| BeH₂ | R_eq 11.7% | Composed, full 1-RDM exchange | 17 |
| H₂O | R_eq 19.4% | Composed, 5-block, zero parameters | 17 |
| LiH (4N) | R_eq 63.5% | Full 4e mol-frame, PK-free | 17 |
| Composed Pauli | N = 27.90xQ exact | 54x-317x equal-qubit vs Gaussian, 35 molecules (exact rule 2026-08-29) | 14 |
| Atomic Pauli | ~Q^3.8 | 0.54x @ Q10 (unfavorable) to 1.5x @ Q28 vs cc-pVDZ/cc-pVTZ (exact rule) | 14 |
**Key structural results (details in papers and CHANGELOG.md):**
- **κ = −1/16 (v2.26.1; Observation, not derived — v4.13.0 QA):** the matching κ *coincides* with the geometric 1/16 of the Fock projection (a numerical observation; no derivation bridge — never "derived"). Not fitted. See `debug/probe_kappa_sprint_memo.md`, `memory/kappa_observation_not_derived.md`.
- **α structural decomposition (Phases 4B-4I, April 2026):** B=42 (Casimir), F=π²/6 (Fock Dirichlet at d_max), Δ=1/40 (Dirac degeneracy g_3). Three independent spectral homes; combination rule K=π(B+F−Δ) is a numerical observation, not derived (conjecture→observation downgrade 2026-06-14). 12 mechanisms eliminated. See Paper 2, CHANGELOG.
- **Spectral-action supertrace (v2.26.1):** SD cancellation theorem, Δ⁻¹=40 from Euler-Maclaurin, (−) sign = (-1)^F grading. Two-term exactness on S³ Dirac. See `debug/st_supertrace_sprint_memo.md`.
- **Nuclear systems (Paper 23):** Deuteron 16q/688 Pauli; He-4 16q/828 Pauli; composed nuclear-electronic deuterium 26q/710 Pauli (corrected 2026-08-22, v5.0.0; the retired 592/712/614 carried the moshinsky N_tot guard). Fock rigidity theorem (S³ unique to −Z/r).
- **Angular sparsity theorem (Paper 22):** ERI density depends only on l_max, not V(r). Universal across potentials.
- **Bargmann-Segal lattice (Paper 24):** HO on S⁵ Hardy space, π-free, HO rigidity theorem. Coulomb/HO asymmetry = 7 layers (7th: CM Hodge structure, v4.103.0; χ₋₄ periods transfer, ¼ of P28).
- **S⁵ gauge extension (v2.26.1):** U(1) Wilson transfers; SU(3) NOT natural; CP² quotient fails. See `debug/s5_gauge_structure_memo.md`.

**Classical solver status: INVESTIGATION COMPLETE (v2.0.24).** 30+ tracks exhausted all solver × PK × basis combinations. Structural ceilings characterized. Composed at l_max=2 is the production operating point. See CHANGELOG.md for full track history.

**Quantum computing status: ACTIVE FRONTIER (Paper 14).** Exact-rule numbers (2026-08-29): composed N_Pauli = 27.90×Q exactly linear across molecules (30.03 d-block), within-molecule exponent 3.17 universal (local slope → ~3.8), 54×–317× equal-qubit vs Gaussian (H₂O), ecosystem export (OpenFermion/Qiskit/PennyLane). Market test: composed LiH near parity with STO-3G raw (838 @ Q30 vs 907 raw; 3.0× more than reduced 276; 1-norm 1.007×, 34.5 vs 34.3 Ha) at more qubits; decisive wins are the exact Q-linearity + 76× vs cc-pVDZ + structural sparsity. Atomic path: exponent ~Q^3.8, and the equal-qubit advantage INVERTS at Q=10 (288 vs 156 cc-pVDZ = 0.54×) recovering to 1.5× at Q=28 — disclosed honestly per the benchmarking rule. (CF-1 is DISSOLVED: the pair-diagonal "convention" was a wrong-sign-q bug; every prior 11.10×Q / O(Q^2.5) / 51×–1712× / 334 / 190× / 2.51×-constant figure is retired and C17-registered.) 37-system `hamiltonian()` library (35 composed molecules + He + H2). See Paper 14, Paper 20, CHANGELOG Tracks AW-CA.

**Historical arc chronicles** (entropy/EP-2; cusp; balanced coupled; precision-atomic; RH sprint + S_min erratum; and the QED / spectral-triple-WH1 / Lorentzian / gauge / Dirac-on-S³ / chemistry / precision-multi-focal / gravity arcs) moved verbatim to `docs/development_frontier_archive.md` (2026-06-16 compaction). The §6 loading tiers + `papers/INDEX.md` are the current map; CHANGELOG.md is the canonical chronicle; the §2 one-liner index above + the best-results / key-structural tables are the live current-state quick reference.

**Scope:** See `SCOPE_BOUNDARY.md` for supported atoms/molecules (the 37-system `hamiltonian()` library; Z=1-56 via frozen cores).

**Architecture locked:** The LCAO/graph-concatenation approach (v0.9.x series) is superseded. All molecular work uses natural geometry (Papers 11, 13, 15, 17).

---
## 3. Approaches That Failed

Critical institutional memory. **Do not re-derive these dead ends.** Category table below, then a name+count index of the 93 per-sprint instances; the **full account of every row is `docs/failed_approaches_ledger.md`** (moved verbatim). Refinement layer: `docs/walls/register.md`. Sub-agents must check the ledger first.

|:---------|:-----:|:-----------|
| Category | Count | Key Lesson |
| LCAO / single-S³ molecular encoding | 3 | Graph Laplacian kinetic energy is R-independent; need natural geometry where separation occurs |
| PK modifications (projector, spectral, self-consistent, R-dependent, l-dependent on 2D) | 6 | PK provides coordinate-space exclusion that angular projectors categorically cannot replicate; PK is the irreducible cost of composed-geometry factorization |
| Cusp treatments (alpha-only, graph absorption, θ₁₂-adapted basis) | 3 | Cusp is 2D in (α, θ₁₂); 1/r₁₂ is an embedding exchange constant; Schwartz extrapolation is the correct approach |
| Inter-group antisymmetry (node exclusion, fiber bundle, spectral) | 3 | Antisymmetry requires shared coordinate system; composed geometry factorizes into incompatible coordinates |
| Full N-electron radial solvers (adiabatic, coupled-channel, 2D variational) | 3 | Adiabatic overcounts D_e; coupled-channel numerically unstable; 2D gives unbound D_e; angular basis is the bottleneck |
| Polyatomic coupling (Z_eff partition, classical repulsion, lone pair at Z_eff>4) | 3 | Orbital-level exchange coupling needed; lone pair Slater integrals unphysical at high Z_eff |
| Geometric elevation (blow-up, Lie algebra, S³×S³) | 1 | min/max boundary is physical; per-ρ diagonalization structurally irreducible |
| Diagnostic arcs (eigenchannel rotation, spheroidal compression, enhanced Z_eff, midpoint origin) | 4 | Various; see CHANGELOG.md |
| l_max convergence via 2D solver (variational_2d in composed pipeline) | 1 | l_max divergence is structural to PK/composed architecture, NOT from adiabatic approximation. |
| Sturmian CI (Coulomb, generalized) | 2 | Coulomb improves atomic FCI but Lowdin 1-norm inflated 2.8-4.5x; generalized loses within-config flexibility at larger basis; neither bypasses PK ceiling |
| TC Jastrow in adiabatic hyperspherical solver | 1 | TC replaces multiplicative V_ee (O(R) in angular eigenvalue) with first-derivative G_ang (O(1)); adiabatic separation requires V_ee as multiplicative operator; ~46% error; TC needs direct variational/FCI framework, not adiabatic. |
| TC in second quantization (composed qubit pipeline) | 1 | TC qubit Hamiltonians plateau at ~3.4% across n_max=2..5; standard FCI converges 5.3%→2.0%. |
| TC 2D variational cusp correction (basis doubling) | 3 | Three attempts failed (R sin(α) basis worsens, sin(α) negligible, [H,r₁₂]=0 identically). |
| TC angular gradient for l>0 orbitals | 1 | Adds 2.66× Pauli terms for 0.01 pp accuracy. |
| Coupled composition (cross-block ERIs replacing PK) | 1 | Cross-block ERIs add 2.56× Pauli, 2.30× 1-norm; 29% FCI error. |
| Single-center nested molecular LiH | 1 | R_eq 33.7% error; single Z can't represent both core and bond length scales. |
| Two-center charge-center nested LiH | 1 | 48.2% energy error; truncated basis cannot represent 1s² core off-center. |
| Heterogeneous nested (per-pair Z_eff) | 1 | Löwdin orthogonalization destroys Gaunt sparsity (1711 vs 120 Pauli, 14× inflation). |
| Fock energy-shell self-consistency for He | 1 | k²=−2E over-constrains 2-electron problem; SC 12.8%/5.1% vs variational 1.9%/1.6%. |
| SM-running origin for Δ = 1/40 (Paper 2 alpha) | 6 | Six tracks ruled out SM-running origin for Δ=1/40. |
| Schwartz-tail hot-node patch on He Z=2 graph-native CI (CUSP-2) | 1 | Schwartz tail correction worsens accuracy; 0.20% floor is small-Z graph-validity-boundary artifact (Z_c≈1.84), not cusp. |
| Energy graph for V_ee on S³ (Paper-12 analog search) | 1 | Pair-state graph dense (47%); cross-shell denominators don't close. |
| Darwin+MV for He/Li/Be 2p-doublet improvement | 1 | Both 2p states share l=1; Darwin=0 for l≥1; MV cancels in splitting. |
| Balanced + frozen-core PES for second-row molecules | 2 | NaH/MgH₂ overattract monotonically at n_max=2; frozen [Ne] hides core screening from cross-center V_ne. |
| Single-constant graph-to-continuum QED projection (C×F₂→α/(2π)) | 1 | C×F₂ grows with n_max, doesn't converge to α/(2π). |
| σ-vertex and direction-resolved vector QED on Fock graph (VQ-1..VQ-5) | 5 | Five tracks ruled out σ-vertex approaches to vector photon selection rules. |
| Co-exact mode q-labeling for SO(4) channel count / Ward / charge conjugation | 1 | Co-exact eigenvectors have right eigenvalues but wrong support (nearest-neighbor only); triangle inequality fails at high q. |
| Dirac-sector lift of Paper 2 α combination rule ingredients B, F, Δ | 3 | Three obstructions: B doesn't lift (zero mode forces (m−1)), F doesn't lift (Apéry Q-linear independence), Hopf-equivariant decomp doesn't produce B/F. |


### Per-sprint instances (index; full rows in the ledger)

| Approach (full row in the ledger) | N |
|:---|:---:|
| TC angular-gradient ERI assembly (2026-08-30) | 1 |
| Multi-focal spatial composition (Sprint HF, 2026-05-07) | 3 |
| Phillips-Kleinman cross-center barrier for second-row chemistry binding (2026-05-08) | 1 |
| Screened-Schrödinger valence basis (h1 diagonal only, Track 3, 2026-05-09) | 1 |
| W3 spectral-zeta calibration-data identification (Sprint W3, 2026-05-08) | 3 |
| Multi-focal Path C5 saturated basis (numerical-luck f closure, 2026-05-09) | 1 |
| Heuristic two-zeta screening for [Xe] core (CR67 + BBB93 Kr ratios, 2026-05-09) | 1 |
| All-positive-coefficient single-zeta hydrogenic basis lacks radial nodes (Track 2 diagnostic, 2026-05-09) | 1 |
| Multi-zeta physical Na valence basis substitution (Sprint α-Multi-zeta + α-PES, 2026-05-23) | 1 |
| Three-bucket M-Z partition (basis-closable cross-shift) FALSIFIED at NaH max_n=3 (2026-05-23) | 1 |
| Kernel-shape substitution as W1c-residual closure (Sprint F2, 2026-05-23) | 1 |
| Single-particle Pauli orthogonality (rank-1 PK on bonding orbital) as W1e closure (Sprint F4, 2026-05-23) | 1 |
| Mean-field core-bonding J-K (explicit-core Hartree) as W1e closure (Sprint F5, 2026-05-23) | 1 |
| Basis-level Schmidt orthogonalization of H 1s against Na [Ne] core as W1e closure (Day-1 diagnostic, 2026-05-23) | 1 |
| [Ne] core correlation as W1e closure mechanism (Day-1 literature estimate, 2026-05-23) | 1 |
| Basis enlargement to max_n=4 alone as W1e closure (Sprint F6, 2026-05-23) | 1 |
| UV refinement to close α > 1 spinor SC gap (G4-4c week 2, 2026-05-29) | 1 |
| Cutoff-function cure for IR over-count of $S_{\rm BH}$ at small Λ (G4-5c-IR-fix, 2026-05-29) | 1 |
| Attributing the G7/G4-2 Newton-constant factor of 2 to a scalar-vs-Dirac cone coefficient (2026-05-30) | 1 |
| Scalar (periodic BC) conical-defect SC extraction at sprint scale (G4-3c-proper / T1, G4-4e, 2026-05-28/29) | 1 |
| Naive K_var/K_disk → 4(l_max+1)(l_max+2) saturation prediction (G4-4b-c, 2026-05-29) | 1 |
| Small-t polynomial fit for Seeley-DeWitt $a_0$ extraction (G4-4d, 2026-05-29) | 1 |
| Naive $\phi(2)$ Mellin-moment prediction for $S_{\BH}$ cutoff dependence (G4-5d, 2026-05-29) | 1 |
| Bulk $\Lambda^4$ cosmological-constant extraction from 2D disk-Dirac alone (G4-5b, 2026-05-29) | 1 |
| v3.19.0 "Fursaev-Solodukhin spinor double-cover correction" mechanism attribution for Möbius α/(2α-1) (task #26, 2026-05-29) | 1 |
| Geometric-mean azimuthal as cheap UV cure for FD/Spec bracket (task #27, 2026-05-29) | 1 |
| Single-axis N_φ-sweep for G4-6a refined A-coefficient extraction (2026-05-29) | 1 |
| Isotropic multi-axis refinement for A-coefficient extraction (2026-05-31) | 1 |
| soft_IR_frac → 1/(2α) as the Möbius α>1 "mechanism" (Route C, demoted by Sprint GD-2 t-audit 2026-05-29) | 1 |
| Spurious "sign discrepancy" Track α'' thread 9 (2026-05-29, RESOLVED same day) | 1 |
| Möbius α/(2α−1) as a continuum α>1 closed form (B4, 2026-05-29, RETIRES the open thread) | 1 |
| Naive de-compactification of the finite-R Dirichlet disk (Paper 53 B3, 2026-05-29) | 1 |
| Species-II spatial fission-aperture discriminator (which-site entanglement frozen=non-binder / responsive=binder, 2026-05-30) | 1 |
| BW wedge entanglement entropy as Bekenstein-Hawking area law (BH-Phase0, 2026-05-31) | 1 |
| Furnstahl IR extrapolation for E1 polarizability (2026-06-01) | 3 |
| Sturmian basis for deuteron polarizability with Minnesota NN (2026-06-01) | 3 |
| LIT for deuteron E1 polarizability (2026-06-01) | 1 |
| NaH Z_orb scan as chemistry-wall diagnostic (2026-06-01) | 1 |
| Resolvent (D²)⁻¹ for two-body Coulomb interaction (2026-06-01) | 4 |
| J_GV²=−1 as a forcing handle for the inner ℍ factor (Door 4c, 2026-06-01) | 1 |
| Gauged tensor-product spectral action (full double-sum gauge field) for two-body Coulomb radial weights (2026-06-03) | 1 |
| Yukawa values in low-coefficient pure-Tate periods M1 ∪ M2 (Sprint Yukawa-PSLQ, 2026-06-03) | 1 |
| Hopf-tower-to-representation extension as shortcut for forcing N_gen = 3 (Sprint Read 2 scoping, 2026-06-03) | 1 |
| W1e chemistry corrections as outer-factor M1/M2/M3 periods (Sprint W1e period-class, 2026-06-04) | 1 |
| Propinquity-derived Trotter bound at production parameters (Sprint Trotter propinquity, 2026-06-04) | 1 |
| DMRG-on-FCIDUMP closes W1e on NaH (Sprint R3-B falsifier, 2026-06-07) | 1 |
| LiH composed qubit FCI binds (R3-A diagnostic, 2026-06-07) | 1 |
| Existing engineering kwargs close LiH over-binding (LiH kwarg sweep, 2026-06-07) | 1 |
| Explicit-core HF closes NaH binding (Sprint B.1, 2026-06-07) | 1 |
| W1e closure via off-diagonal cross-block h1 WITHOUT orthogonalization (SO(4)-breaking diagnostic, 2026-06-07) | 1 |
| Spectral action expansion has chemistry-side structural content at finite cutoff (Sprint spectral action expansion, 2026-06-07) | 1 |
| Camporesi-Higuchi κ-parity Z₂ as relativistic chemistry tapering stabilizer (Sprint CH κ-parity, 2026-06-07) | 1 |
| Direct-basis m_j-parity Z₂ as relativistic chemistry tapering stabilizer (Sprint m_j-parity direct, 2026-06-08) | 1 |
| Rotated-basis m_j-parity Z₂ as relativistic chemistry tapering stabilizer (Sprint m_j-parity rotated, 2026-06-08) | 1 |
| Spectral action $S(D)(R)$ of assembled MvS Dirac binds LiH (Sprint M-vS-2 + R-sweep, 2026-06-07) | 1 |
| Non-abelian M-vS gauge group reduces Pauli count beyond Z₂ Hopf-U(1) tapering (Sprint M-vS Gauge, 2026-06-07) | 2 |
| K⁺ compression of the Krein Dirac as a Lorentzian quantum-metric device (P45 descope, 2026-06-09) | 1 |
| Momentum-diagonal temporal multipliers as a metric temporal algebra (P45, 2026-06-09) | 1 |
| Grab-bag / single-weight PSLQ bases for weight-inhomogeneous graph-QED constants (2026-06-11) | 1 |
| Levin/EM series acceleration on log-modulated summands (2026-06-11) | 3 |
| Nested adaptive nsum cascade for depth-≥4 multiple-t-values (S^(4) stage-1, 2026-06-13) | 1 |
| `mpmath.sumem` (auto Euler-Maclaurin) for log-modulated trailing-t tails (S^(4) stage-1, 2026-06-13) | 1 |
| float64-contaminated terms inside a high-precision Levin sum (S^(4) stage-1, 2026-06-13) | 1 |
| Brute high-precision summation of b1=2 high-log trailing multiple-t-values (S^(4) stage-1, 2026-06-13) | 1 |
| sympy exact-Q rank for MZV-scale relation matrices (S^(4) stage-2, 2026-06-13) | 1 |
| SO(4) Wigner-D critical points predict LiH R_eq (harmonic phase-lock probe, 2026-03-11; surfaced 2026-06-16) | 1 |
| Paper 56 closed-immersion injectivity via per-sector η period map (2026-06-16) | 1 |
| Nested 2D adaptive quadrature for two-center overlap SUPPORT questions (Topos-3 attempt, 2026-07-05) | 1 |
| Graph counting/degeneracy measure as Born-rule generator (WH8 Step-1, 2026-07-01) | 1 |
| Overlap-slope law for the balanced-solver R_eq tilt (2026-07-06) | 1 |
| Löwdin retrofit of GeoVac integrals — the "sparsity-destroying option" (2026-07-07) | 1 |
| Non-orthogonal fermionic encoding to preserve Gaunt l-sparsity (2026-07-07) | 1 |
| "Exact up to one known seed" as the two-center ERI deliverable (2026-08-12) | 1 |
| Basis compactness as a quantum-resource advantage for Slater/Sturmian (QC-1, 2026-08-12) | 1 |
| Closing the *reducible* half of the three-centre ERI block first and deferring the hard half (Poly-2, 2026-08-14) | 1 |
| Dropping the three-centre ERI block to reach polyatomics with a two-centre engine (Poly-0, 2026-08-13) | 1 |
| Point-group symmetry (gerade lever) to rescue polyatomic ground-state metric conditioning (water C₂ᵥ, 2026-08-18) | 1 |
| Brute high-precision T2 by integrating the fibre over the family (2026-08-20) | 1 |
| ℚ(i) seam (P56 Dirac ↔ P59 T2) as a shared mechanism (2026-08-21) | 4 |
| Seam involution/fixed-point mechanism — τ=i pinned by a physical ℤ₂ (2026-08-21) | 1 |
| Low-rank/Woodbury flattening of the unique-center SW metric penalty (2026-08-21) | 1 |
| Analytic momentum-space factorization as a 1-norm lever ("analytic THC", 2026-08-21) | 1 |
| HO total-quanta conservation as an entanglement-rigidity mechanism (Papers 24 + 27, 2026-08-22) | 1 |
| Spectrum-aware QSVT interpolation as a scaling change for S^(−1/2) (2026-08-21) | 1 |
| TC three-body operator collapse via Gaunt/6j selection rules (2026-08-23) | 1 |
| Elliptic/CM-adapted basis for chemical accuracy (2026-08-23) | 3 |
| Symmetrising the TC operator to dodge the non-Hermitian quantum axis (2026-08-23) | 1 |
| Cheap + accurate electron-cusp on the Coulomb-Sturmian basis (2026-08-23) | 3 |
| Fixed / physics-determined GLOBAL γ for the cheap TC dressing (2026-08-24) | 1 |
| Position-dependent γ(r)∝n^{1/3} baked into the cheap non-Herm TC operator (2026-08-24) | 1 |
| Q̂₁₂ strong orthogonality as a quantum-encoding win for R12-CI (2026-08-26) | 1 |
| Per-pair γ / a second correlation length as the many-electron cusp fix (2026-08-26) | 1 |
| Hydrogenic per-n scaling (k_n = Z/n) as an L²-orthonormal Sturmian replacement (2026-08-26) | 1 |

---

## 3.5. Guardrail Papers

Certain papers contain proven theorems or exhaustive diagnostic arcs that constrain what approaches are viable. Sub-agents MUST load the relevant guardrail paper before beginning any investigation in the constrained domain. If a proposed investigation falls within a guardrail's scope, the sub-agent must:

1. Load the guardrail paper into context
2. Identify the specific theorem or negative result that applies
3. WARN THE USER: "Paper [X] documents a proven negative result for this class of approach: [theorem statement]. The proposed investigation [specific proposal] appears to fall within scope. Do you want to proceed anyway, or modify the approach?"
4. If the user proceeds: document in the investigation prompt that the guardrail was acknowledged and state why the new attempt differs from the documented negative result

This protocol exists because Track DF (6 sprints, April 2026) re-derived at significant cost two negative results that were already proven in Papers 8-9 and FCI-M. The re-derivation produced genuine new findings (6j sparsity theorem, compact encoding, three-layer structure), but Sprints 4 and 4B could have been avoided entirely if the guardrail papers had been loaded at the investigation design stage.

### Current Guardrail Papers

| Paper | Domain | Theorem/Result | Scope |
|:------|:-------|:---------------|:------|
| Paper 8-9 | Single-center molecular encoding | Sturmian Structural Theorem: in a shared-p₀ basis (H,S) share one SO(4) congruence → eigenvalues R-independent (no binding without R-dependent β_k(R)). **Scoped (v4.72.x, PI-reviewed):** the R-independence holds for the single-n D-matrix *approximation* (cross-block *imposed* = orthogonal D^(n); the backing test hard-codes it). The genuine cross-n Coulomb-Sturmian *overlap* couples n (computed: ⟨χ_1s\|χ_2s⟩=−0.47) and escapes the theorem — it binds (Avery SW closed forms; Herbst-Avery-Dreuw) — but hits the qubit-encoding l-mixing wall (see §3 non-orthogonal-encoding row). Paper 8 §Remark + `test_paper8_overlap_cross_n.py`. | ANY investigation proposing to encode a heteronuclear molecule in a single-center (single-Z, single-k, single-p₀) basis, including nested hyperspherical, bond sphere, Sturmian CI, or unified orbital approaches |
| FCI-M | Graph-concatenation molecular encoding | Graph Laplacian kinetic energy is R-independent in LCAO basis → monotonically attractive PES, no equilibrium | ANY investigation proposing to concatenate atom-centered graphs into a molecular graph without natural geometry coordinates |
| Track DF record | Nested hyperspherical molecular encoding | Three molecular variants tested (single-center, charge-center, heterogeneous Löwdin) — all NEGATIVE | ANY investigation proposing to place all molecular electrons in a single S^(3N-1) Hilbert space. Extends Paper 8-9's theorem to the hyperspherical setting with additional finding: Löwdin orthogonalization of mixed-exponent bases destroys Gaunt sparsity |

### Guardrail Triggers

When a user or investigation prompt proposes work in any of the following categories, the sub-agent MUST load the indicated guardrail paper BEFORE writing any implementation plan:

- **"single center" + "molecule"** → Load Papers 8-9
- **"unified basis" + "molecule"** → Load Papers 8-9
- **"shared exponent" OR "single k" OR "single p₀" + "molecule"** → Load Papers 8-9
- **"nested" + "molecule"** → Load Papers 8-9 + Track DF record
- **"graph" + "concatenat"** → Load FCI-M
- **"LCAO" + "graph"** → Load FCI-M
- **"Sturmian" + "molecule"** → Load Papers 8-9

After loading, the sub-agent must check whether the proposed approach falls within the guardrail's scope and warn the user if it does. The warning should be constructive, not blocking — it documents what DOESN'T work, which helps refine toward what MIGHT work.

---

## 4. The Dimensionless Vacuum Principle

The graph Laplacian is dimensionless and scale-invariant, topologically equivalent to the unit three-sphere S3. Physical energies emerge only through the energy-shell constraint p0^2 = -2E, which acts as a stereographic focal length. The 1/r Coulomb potential is not an input force law -- it is the coordinate distortion from stereographic projection (chordal distance identity). The universal kinetic scale kappa = -1/16 maps graph eigenvalues to the Rydberg spectrum. Eigenvalues of the Laplace-Beltrami operator on unit S3 are pure integers: lambda_n = -(n^2 - 1).

For molecules, the natural geometry shifts from S3 (atoms) to prolate spheroidal coordinates (Paper 11), hyperspherical coordinates (Paper 13), molecule-frame hyperspherical (Paper 15), or composed fiber bundles (Paper 17). The choice of geometry is determined by where separation of variables occurs. At Level 1, the separated coordinates are all angular (on S3), making the problem fully discrete. At Levels 2+, separation produces continuous radial-like coordinates alongside discrete angular channels — a consequence of SO(4) symmetry breaking by multi-center or multi-electron potentials (Papers 8-9, negative theorem).

**Prime directive (level-aware):** The framework preserves two distinct kinds of structure, and the rules differ by level:

*Angular/symmetry structure (discrete at all levels):* The quantum number labeling (n, l, m), the per-shell degeneracy 2l+1, the selection rules from Gaunt integrals, and the channel structure are combinatorial invariants derived from the packing construction (Paper 0). These must never be modified, approximated, or bypassed at any level. They are the framework's core identity.

*Radial/parametric structure (level-dependent):* At Level 1, Fock's SO(4) symmetry converts the radial coordinate into an angular coordinate on S3, making the entire problem discrete. The graph Laplacian is the exact object and must not be modified to artificially recover continuous differential terms (like 1/r or nabla^2). At Levels 2-4, the SO(4) symmetry is broken by multi-center or multi-electron physics, and the radial-like coordinates (internuclear distance R, hyperradius R_e, prolate spheroidal xi) currently require continuous numerical methods. However, this is not necessarily a fundamental feature of the physics — it may reflect algebraic structures that haven't been found yet. The project has already demonstrated this replacement in several cases: Paper 12's Neumann expansion replaced numerical quadrature of V_ee with algebraic recurrence relations; Gaunt integrals replace angular integration at all levels; the split-region Legendre expansion (Paper 15) terminates exactly via the 3j triangle inequality. In each case, what appeared to require continuous integration was replaced by algebraic evaluation once the right structure was identified.

The places where continuous numerical methods currently survive are: Z_eff screening quadrature in composed geometries, spline caching for the rho-collapse, finite-difference radial grids in hyperspherical solvers, and the adiabatic potential curves U(R). These are computational tools for evaluating coupling matrix elements between discrete channel structures. Standard numerical methods (grid refinement, basis convergence, quadrature improvement) are appropriate for improving accuracy in these areas. Finding algebraic replacements for these numerical steps — as Paper 12 did for V_ee — is an ongoing goal of the project.

The framework's long-term aspiration is: for every natural geometry, there exists an algebraic structure that computes all coupling matrix elements without spatial integration. This has been achieved for Level 1 (S3 graph eigenvalues), for angular couplings at all levels (Gaunt integrals), and for electron-electron repulsion in prolate spheroidal coordinates (Neumann expansion). It has not yet been achieved for the hyperradial coupling in Levels 3-4 or for Z_eff screening in Level 5. What must be preserved at every level is the discrete channel structure and selection rules, not the specific numerical method used to evaluate radial amplitudes within those channels.

**Algebraic-first, observation-aware (refinement after Papers 18 / 34, May 2026):** The aspiration above is the right discipline for Layer 1 (the bare graph) and for the discrete channel structure at every level. After Paper 34's two-layer framing, the project distinguishes algebraic content (Layer 1: π-free, integer / rational / algebraic-extension matrix elements) from observation-side content (Layer 2: where projections to physical observables introduce specific transcendentals — π via the Hopf measure, π^{2k} via the spectral action, ζ(2k) via even-zeta Dirichlet series, ζ(3) via half-integer Hurwitz, Catalan G via vertex parity, 2π via temporal compactification when an observer integrates over a finite time window). Irreducible transcendentals that survive algebraic decomposition are not failures of the algebraic-first discipline — they are the *content* of a specific Paper 34 projection and should be pinned to that projection rather than chased indefinitely. The diagnostic question for any quadrature wall is therefore two-headed: (i) is the wall a missing algebraic structure (decompose, per Paper 12 / cross-block V_ne / hypergeometric Slater), or (ii) is it the irreducible signature of an observation-side projection (catalogue against Paper 34 and stop)? Both answers are valid outcomes; conflating them was the error the May-2026 curve-fit audit flagged.

*Practical test:* If a proposed modification changes which quantum numbers label the states or which transitions are allowed, it violates the prime directive. If it changes how accurately the radial amplitude is computed within a given channel, it is a legitimate numerical improvement. If it adds a transcendental, the transcendental must be tagged to a Paper 34 projection — anonymous transcendentals are not allowed in production code or papers.

**Topological integrity tests:** The 18 symbolic proofs in `tests/test_fock_projection.py` and `tests/test_fock_laplacian.py` validate S3 conformal geometry. These must never be broken or bypassed. Run before any release.

---

## 5. Natural Geometry Hierarchy

The core organizational principle of the project. Each electron configuration has a natural coordinate system where the physics separates.

| Level | System | Natural Geometry | Best Result | Paper |
|:-----:|:-------|:-----------------|:------------|:-----:|
| 1 | H (1-center, 1e) | S3 (Fock) | lambda_max -> 2 d_max = 8; deficit 0.57% at n_max = 30 (a spectral *bound*, not an accuracy: E_0 = kappa*lambda_max by construction) | 7 |
| 2 | H2+ (2-center, 1e) | Prolate spheroid | 0.0002% (spectral) | 11 |
| 3 | He (1-center, 2e) | Hyperspherical | 0.004% (2D var, cusp l_max=4); 0.022% (raw l_max=7); 0.19% (graph-native CI n_max=7, 0 params, exact algebraic integrals) | 13 |

*Level 3 note:* The previous 0.05% result used the FD adiabatic solver (lucky error cancellation, non-variational). The adiabatic coupled-channel solver converges to a structural floor of 0.19-0.20% (v2.0.8). The 2D variational solver (Track DI, v2.6.0) breaks this floor: raw 0.022% at l_max=7 (tensor-product Laguerre × Gegenbauer basis, 8000 dim), cusp-corrected 0.004% at l_max=4. The 2D solver treats R and α simultaneously, capturing non-adiabatic R-α correlation that the adiabatic approximation misses. l_max convergence is monotonic; the per-channel angular basis (n_basis=40 Gegenbauer functions) is the convergence bottleneck, not partial-wave truncation.
| 4 | H2 (2-center, 2e) | Mol-frame hyperspherical | 96.0% D_e | 15 |
| 4N | LiH (2-center, 4e) | Full mol-frame hypersp. (SO(12)) | R_eq 63.5% (l_max=2, 2D variational; unbound D_e) | 17 |
| 5 | LiH (core+valence) | Composed (Level 3 + 4) | R_eq 5.3% | 17 |
| 5 | BeH₂ (polyatomic) | Composed (Level 3 + 4) + exchange | R_eq 11.7% | 17 |
| 5 | H₂O (triatomic) | Composed (Level 3 + 4) + lone pairs | R_eq 19.4% | 17 |
| Nested | Be (1-center, 4e) | S¹¹ (SO(12)), H-set | Q=10, 112 Pauli, 4.9% | Track DF |

*Level 4N note:* Level 4N is the exact N-electron generalization of Level 4's mol-frame hyperspherical coordinates (SO(3N) replacing SO(6), S_N antisymmetry replacing the gerade constraint). Level 5 (composed geometry) approximates Level 4N, trading exact inter-group antisymmetry for 144x angular compression via PK. The l_max=2 result (R_eq ≈ 1.1 bohr, 63.5% error) demonstrates that equilibrium exists without PK. Track AR (v2.0.23) confirmed that D_e overcounting was an adiabatic artifact: the 2D variational solver gives E_min = -7.79 Ha (variational bound respected, above exact -8.07) with D_e = -0.19 Ha (unbound), versus the adiabatic solver's D_e = +0.49 Ha (5.3x exact). The l_max=2 S₄ [2,2] angular basis is genuinely insufficient for LiH binding. Track AS (v2.0.23) confirmed composed encoding is categorically sparser than full N-electron encoding (334 vs 3,288 Pauli terms, 20x lower 1-norm). Composed geometry's 144x angular compression and structural sparsity are essential, not approximations of convenience.

The composed geometry (Level 5) is a fiber bundle: G_total = G_nuc semi-direct G_core(R) semi-direct G_val(R, core_state). Each electron group gets its own natural coordinate system, coupled via Z_eff screening and Phillips-Kleinman pseudopotential.

**Algebraic structure:** At every level, angular matrix elements are computed from quantum number labels and Wigner 3j symbols (via Gaunt integrals), with no spatial quadrature. The split-region Legendre expansion (Paper 15) terminates exactly via the 3j triangle inequality. At Level 2, the radial solver is fully algebraic for all m: σ states (m=0) use ordinary Laguerre three-term recurrence with zero numerical integration (v2.0.9); π/δ states (m≠0) use associated Laguerre basis L_n^{|m|}(x) with partial-fraction decomposition and Stieltjes integral recurrence, reducing non-algebraic content to a single transcendental seed e^a·E₁(a) (Track J, v2.0.10). At Level 3, the angular problem is fully algebraic, and the hyperradial overlap S and kinetic K matrices are algebraic via three-term Laguerre recurrence (Track H, v2.0.10: pentadiagonal M2 for S, tridiagonal derivative expansion for K, < 1e-14 relative error, 11× build speedup). The adiabatic eigenvalues μ(R) are proven transcendental (O(R) → O(R²) regime transition, v2.0.9 Track G) — the potential V_eff(R) must stay quadrature, point-by-point diagonalization is irreducible, though spectral radial solvers achieve 95-120× speedups. At Level 4, the spectral Laguerre basis achieves 16× dimension reduction for the hyperradial coordinate (Track I), and the spectral Jacobi basis achieves 20× dimension reduction for the angular sweep (Track K, 269× speedup, 1000→50 matrix dimension). The combined spectral solver reduces the angular sweep from 99% to ~50% of total cost. The angular eigenvalues μ(ρ) remain transcendental (computed by diagonalization at each ρ-point), but the per-point cost is now a 50×50 eigensolve rather than a 1000×1000 one. n_basis_radial=20 optimal (mild conditioning at n_basis≥25). Spatial quadrature enters for: Level 3 hyperradial potential matrix elements (V_eff transcendental), Level 4 angular eigenvalue sweeps (μ(ρ) transcendental, but spectral basis reduces cost 269×), Z_eff screening in composed geometries, and rho-collapse spline caching.

---

## 6. Paper Series

> Newcomer-facing status map with one-liners: `papers/INDEX.md`. Detailed per-paper notes (the long key-result descriptions, frozen at v3.110.0): `docs/paper_notes_archive.md`. Topic → paper lookup: `docs/topic_to_paper_lookup.md`. This section keeps only what agents need at dispatch time: loading tiers, the folder map, and live status flags.

### Loading tiers

**Always load** (framework identity): Paper 0 (packing axiom, K = −1/16) · 1 (spectral graph methods) · 7 (S³ proof, 18 symbolic proofs) · 14 (qubit encoding headline) · 16 (S_N periodicity) · 22 (angular sparsity theorem) · 23 (nuclear hub, Fock rigidity) · 24 (Bargmann-Segal S⁵, Coulomb/HO asymmetry) · 27 (entropy as projection) · 31 (universal/Coulomb partition) · 32 (the spectral triple; §VIII theorems).

**Load on topic** (full list and statuses in `papers/INDEX.md`): chemistry solvers → 8–9, 11, 12, 13, 15, 17, 19, FCI-A/M; QC resources → 20; QED/gauge/gravity → 2, 25, 28, 30, 33, 36, 41, 51; math.OA arc → 29, 38, 39, 40, 42–50, 52, 53; foundations/periods → 18, 54, 55, 56, 57; precision → 26, 34, 35.

**GUARDRAIL papers** — MUST load before any investigation in their domain (trigger words and protocol in §3.5): Papers 8–9 (single-center / unified-basis / shared-exponent molecular — Sturmian structural theorem), FCI-M (graph-concatenation molecular), Track DF record (nested hyperspherical).

### Folder organization (audience groups, reorganized 2026-05-22)

| Folder | Audience | Papers (.tex) |
|:-------|:---------|:-------------:|
| `papers/group1_operator_algebras/` | math.OA / NCG | 16 |
| `papers/group2_quantum_chemistry/` | quantum chemists | 11 |
| `papers/group3_foundations/` | mathematical physicists | 11 |
| `papers/group4_quantum_computing/` | QC / NISQ / VQE | 4 |
| `papers/group5_qed_gauge/` | HEP / gauge theory | 8 |
| `papers/group6_precision_observations/` | precision AMO | 4 |
| `papers/synthesis/` | cross-group narratives + field guide | 3 |
| `papers/archive/` | historical (3, 4, 5, 6, 10, 18v1, 21) | 7 |

### Live status flags every agent must know

- **Paper 38 UNCONDITIONAL** (2026-06-10, translation-seminorm metrization; falsifier `tests/test_p38_action_seminorm.py`). The WH1 keystone.
- **Paper 45 DESCOPED + partially rebuilt** (2026-06-09 K⁺ theorem withdrawn, falsifier `tests/test_p45_kplus_degeneracy.py`; 2026-06-10 product-carrier convergence restored in the action-seminorm framework, `prop:product_action_seminorm`, falsifier `tests/test_wh7_b1_joint.py` — signature-agnostic, NOT a Lorentzian claim). Do NOT cite pre-descope claims. **Paper 46 DESCOPED; Papers 47/48/49 PARTIAL** (in-paper Status notes; the norm-resolvent arrow and the TICI/cocycle algebra survive).
- **Paper 2 is an Observation**; the combination rule K = π(B + F − Δ) is labeled an Observation — never conjecture or derived (§13.5 hard prohibition; conjecture→observation downgrade 2026-06-14).
- **Paper 34** is the living projection catalogue (28 projections); **Paper 18 §III.7** is the master Mellin engine; tag every transcendental against both (memory rule).
- Papers are corrected **in place** (de-versioning directive 2026-06-10); git/Zenodo are the version record. No splinter files.

---

## 7. Code Architecture

> Full entry-point catalogue (~120 rows) and solver-method table: `docs/code_architecture.md` (live document — update there, not here). Most-used entry points:

| Task | Module · Entry point |
|:-----|:---------------------|
| Atomic lattice / Hamiltonian | `geovac/lattice.py` `GeometricLattice(max_n, nuclear_charge=Z)` · `geovac/atomic_solver.py` `AtomicSolver(max_n, Z)` |
| Multi-electron FCI | `geovac/lattice_index.py` `LatticeIndex(Z, n_electrons, max_n)`; direct CI at N_SD ≥ 5000 |
| Molecular spec + composed builder | `geovac/molecular_spec.py` `MolecularSpec` · `geovac/composed_qubit.py` `build_composed_hamiltonian(spec)` + `*_spec()` factories |
| Balanced coupled builder | `geovac/balanced_coupled.py` `build_balanced_hamiltonian(spec, nuclei)` |
| Ecosystem export | `geovac/ecosystem_export.py` `hamiltonian(name, tapered=None/'global'/'per_block'/'extended'/'full')` → `.to_qiskit()/.to_openfermion()/.to_pennylane()` |
| Z₂ tapering | `geovac/z2_tapering.py` `apply_hopf_tapering()` · `geovac/extended_tapering.py` |
| Frozen cores / cross-center V_ne | `geovac/neon_core.py` `FrozenCore(Z)` · `geovac/shibuya_wulfman.py` `compute_cross_center_vne()` |
| Slater integrals (exact) | `geovac/hypergeometric_slater.py` `compute_rk_float()` (threshold-dispatched n ≥ 5 → exact Fraction) |
| Operator system / Connes distance / GH | `geovac/operator_system.py` `TruncatedOperatorSystem(n_max)` · `geovac/connes_distance.py` · `geovac/gh_convergence.py` |
| Physical constants | no central module; `-1/16` may be used directly (§8); `ALPHA`/`C_LIGHT` live next to their modules |

---

## 8. Coding Standards

### Sparse vs Dense: Context-Dependent

- **Hamiltonian and CI matrices (N > 100):** Always `scipy.sparse` (csr_matrix, coo_matrix). Never densify.
- **Hot-loop lookup tables (ERI, h1 in direct CI):** Use dense NumPy when array fits in memory (n_spinorb <= ~300). `scipy.sparse._validate_indices` overhead (~24us/call) is prohibitive at 100K+ lookups.

Rule of thumb: sparse for the physics matrix (N_SD x N_SD), dense for orbital-index lookup tables (n_spinorb x n_spinorb or n_spatial^4).

### Type Hints Required

All function signatures must have type hints.

```python
def compute_ground_state(self, n_states: int = 1) -> Tuple[np.ndarray, np.ndarray]:
    ...
```

### Physical Constants

Import from `geovac.constants` or define at module top. No hardcoded magic numbers.
**Exception:** `-1/16` is the universal topological constant (can be used directly).

### Vectorization Over Loops

Avoid Python loops for graph operations; use NumPy masking/vectorization.

```python
mask = (n_values >= 1) & (l_values < n_values)
states_filtered = states[mask]
```

---

## 9. Workflow Protocols

### Theory Check Rule

Before implementing new physics:
1. Check `papers/` for the derivation
2. If code contradicts paper -> flag it and ask user
3. If changing physics in code -> prompt user to update papers

### Current-State Check

Before forming or reporting a verdict — **or proposing a new connection or direction from corpus state** (e.g. an /aha generative pass) — on any question, *especially* when resuming a thread from a `debug/` memo, verify the CURRENT state, not the snapshot. The papers (the section that *owns* the question) + CHANGELOG since the memo's date are canonical; `debug/` memos are dated snapshots the corpus moves past, and a memo's "open question / next step" may already be closed. First move when picking up a thread = read the owning paper section + post-memo CHANGELOG, THEN conclude *or propose*. (Added 2026-06-14: twice in one session a verdict was formed from an 8-day-stale synthesis memo — an already-proven theorem was reported as a future sprint, and an already-settled A-vs-B question was re-answered wrong; both answers were live in Paper 56. Widened 2026-07-08: generation is not exempt — an /aha pass pitched a connection as novel that the Marcolli–vS chemistry arc had already established a month earlier; the `/aha` skill quarantines that verification to its Phase B (B1). Standing rule: `memory/feedback_verify_current_state.md`.)

### Branch QA Review Protocol

The corpus is QA'd branch by branch. Cycle, in order: **(1)** synthesis update; **(2)** adversarial paper review (overclaim, §1.5 rhetoric, zombie citations, status drift); **(3)** re-run 1–2 if the audit forces material synthesis changes; **(4)** adversarial *code* review, one `code-reviewer` per paper — map each claim to its backing test, **RUN the tests**, and audit whether the test actually *proves* the claim; **(5)** disposition. Run in **dependency order — trunk roots (Papers 0, 1, 7, 32, 38) before the branches**, so a finding at a root re-prices everything above it.

**Claim → artifact rule (hard).** Every load-bearing paper claim maps to a backing test, recorded in `docs/claim_test_matrix.md`. A claim with no backing test is a **coverage gap** — logged, and raised to the PI if load-bearing — never a silent omission. New equations follow the `test_paper{N}_*` convention.

**Retraction → dependents rule (hard, added 2026-09-04, PI-directed).** When a claim is retracted, corrected, or re-tiered, its `check_retracted_terms.py` entry must declare `cited_by`: the documents whose *argument* rests on it — distinct from `files`, which is only where its wording might appear. Each dependent carries a review stamp or `None`; unstamped dependents fail the gate, and a new entry declaring nothing fails too (`cited_by: {}` if genuinely none). **Stamp outcomes, not intentions** — an edit is not a review. *Why, measured:* eight of the ten recurring defect classes in the v5.4.4..v5.7.3 arc had one shape — a claim owned by one document, cited by several, corrected in the owner, left standing in the citers (“monotone from below” was corrected in Paper 0 and left false in Paper 7 *in the same commit*). Patterns cannot catch this: a citer restates the claim in its own words, so four successive pattern rebuilds each missed loci and one sweep reported clean while five loci survived a spelling difference. A dependency list is enumerated once, from the argument, and cannot be defeated by spelling. *Not covered:* the mirror direction, where an owner **strengthens** a claim and the citers keep its weaker form — a blocklist has nothing to say about a claim that is merely too weak; that needs its own registry and is owed.

**Disposition.** PM fixes small issues directly (status drift, cross-ref hygiene, precision, missing caveats); PM raises to the PI any load-bearing claim with no/weak/false-positive backing, a test that proves less than the prose, a suspected bug in a keystone, or anything touching a hard prohibition.

Full protocol + the three QA principles (provenance visibility / fresh adversary / two-way verdict): `docs/branch_qa_protocol.md`. The gate is `/qa` (PI-invoked only) against `docs/qa/criteria.md`.
### Benchmarking Rule

After any modification to production code in `geovac/`:
1. Run `/regression` (default scope `touched`) — derives the test selection from `git diff` + import graph (consumer test files of every touched module) plus the topological-integrity baseline plus a small reproducible random sample. 30s–2min typical wall.
2. If the diff spans more than 2–3 modules, widen the `touched` selection — do **NOT** reach for `/regression full` at sprint close. **Measured 2026-08-31: the full scope is ~6.7 h serial, ~2.4 h with `-n auto`** (9,156 tests; 98 of 361 files over a 12 s budget; the cost is inherent to the Hamiltonian-building tests, not incidental). It is a **scheduled baseline**, not a close gate. Always pass `-n auto` when you do run it. See `debug/qa/test_suite_cost_memo.md`.
3. Verify the 18 symbolic S³ proofs pass (always included in `touched` and `topo` scopes).
4. Verify H2+ < 0.1% error (topological control).
5. Verify H2 Full CI < 1.0% error (accuracy control).
6. Report any speed regression > 10%.

The previous narrow 3-file allowlist (`tests/test_fock_projection.py`, `tests/test_fock_laplacian.py`, `tests/advanced_benchmarks.py`) silently let test rot accumulate when refactors cascaded into the dozens of test files that import from `composed_qubit`, `inter_fiber_coupling`, etc.  `/regression touched` removes the consumer-selection bottleneck by deriving it mechanically from the diff.  Use it as the standard discipline after any code edit.

### Clean Room Rule

- Generated plots -> `debug/plots/` or `papers/figures/`
- Generated data -> `debug/data/`
- Scripts -> `debug/`, `demo/`, `tests/`, or `benchmarks/` (never root)
- Documentation -> `docs/`

**Papers cite the permanent record, not transient `debug/` (policy, 2026-06-17 PI direction).** A paper must not cite a `debug/` sprint memo or script as backing. `debug/` is the transient clean-room dir (the bullets above) and is pruned over time, so a `\texttt{debug/...}` citation in a paper goes stale *by design* and leaves a dangling pointer (a reader following it hits a deleted file). The permanent homes are: the paper's own derivation, `CHANGELOG.md` (the sprint chronicle, §13.11), and `tests/` (the regression backing). The deterministic **C14** check (`debug/qa/check_file_refs.py`) gates the permanent code/artifact dirs (`geovac/`/`benchmarks/`/`demo/`) and reports `debug/` refs as **advisory** under this policy. *Standing debt (deferred sprint):* a corpus-wide sweep to neutralize the existing ~443 dangling `debug/` references (concentrated in papers 34/32/28; surfaced by the v4.21.0 C14 corpus sweep) is scheduled separately — the policy is adopted now, the ref-removal is its own task.

### Changelog Protocol

**Changelog granularity:** Each CHANGELOG entry should correspond to a releasable unit of work — a new feature, a completed diagnostic arc, a paper update, or a benchmark result. Do not create a new version entry for intermediate debugging steps, parameter sweeps, or exploratory runs. If a task requires multiple iterations to complete, document it as a single entry when the task concludes, noting the key findings and negative results. Intermediate data belongs in `debug/` with timestamped filenames, not in the CHANGELOG.

**Version numbering (revised 2026-08-22, PI direction).** **Default: bump the LAST number only (x.y.Z → x.y.Z+1).** The PM does not choose minor or major bumps; those are PI calls, made explicitly. The point of the change is to make version jumps *mean* something: when the second or first number moves, a reader should be able to infer that a significant corpus change happened, without having to read the CHANGELOG to find out. Under the old rule the PM bumped the minor version for every completed arc, which made minor bumps routine and therefore uninformative (v4.77 → v4.109 = 32 minor bumps).

So: patch (x.y.Z) is the standing default for everything — bug fixes, documentation, completed diagnostic arcs, paper updates, benchmark results, sprint closes. Minor (x.Y.0) and major (X.0.0) are reserved for the PI to call, and mark corpus-significant events: a retraction that moves published numbers, a change to the QA gate or the agent protocol, an architectural change, a reorganization of the paper series. If a sprint feels like it warrants more than a patch, say so in the session summary and let the PI decide — do not bump it unilaterally. Granularity is unchanged: a diagnostic arc that tests 10 hypotheses and finds 9 negative results is ONE version entry, not 10.

---

## 10. Validation Benchmarks

See `docs/validation_benchmarks.md` for the full benchmark table. Update that file (not this section) when adding new benchmarks.

<!-- Full table extracted to docs/validation_benchmarks.md on 2026-05-31 to reduce context-load cost. -->

---

## 11. Topic-to-Paper Lookup

See `docs/topic_to_paper_lookup.md` for the full table mapping topics to papers, sections, and loading tiers. Update that file (not this section) when adding new topic→paper mappings.

<!-- Full table extracted to docs/topic_to_paper_lookup.md on 2026-05-31 to reduce context-load cost. -->

---

## 12. Algebraic Registry

Tracks which matrix elements at each level are computed algebraically vs numerically. **Full registry** (Levels 2 / 3 / 4 / 4N / 5 + spin-ful Tier-2 tables): `docs/algebraic_registry.md` (live document — update statuses there, not here). Status vocabulary: **algebraic** (closed-form from quantum numbers) / **algebraic (implicit)** (defined by P = 0 with known coefficient ring; pointwise diagonalization is convenience, not necessity) / **algebraic-pending** (route identified, production still uses quadrature) / **numerical-required** (no known algebraic replacement). The §4 prime-directive test governs changes: anything touching quantum-number labels or selection rules is prohibited; improving radial-amplitude evaluation within a channel is legitimate.

---
## 13. Multi-Agent Protocol

The GeoVac project uses an AI-augmented agentic workflow with a formalized four-layer architecture for research direction, planning, and execution.

**13.1 Architecture (moved).** Four layers: research agents (`agents/*.md` — Leader / Explorer / Decomposer / Reviewer, PI-invoked), plan mode, the PM session, opt-in worker sub-agents. Detail + cost policy: `docs/multi_agent_protocol.md`.

### 13.2 PM Session Kickoff

Every PM session begins by reading CLAUDE.md and then executing the following:

1. Identify the current track(s) and relevant papers from the plan mode directive
2. Read those papers and any results from the previous session
3. Check the failed approaches summary (Section 3); if the current track touches PK, cusp, inter-group antisymmetry, or molecular encoding, read the full details in CHANGELOG.md before proceeding
4. Plan the session as main-session work (sequence of Read / Edit / Bash / etc.). Flag any tracks that are candidates for sub-agent dispatch (parallelizable, context-heavy) and ask the PI before launching
5. Identify which papers need updating based on the session's results (see 13.8)

**Sprint standard (2026-05-26 update):** The default workflow is **main-session work**. The PM reads, edits, computes, drafts memos, and commits directly. Sub-agent dispatch is opt-in:\ either the PI explicitly directs ("dispatch this in parallel," "run Explore for this lookup") or the PM identifies a case where dispatch genuinely beats sequential main-session work and asks first. When sub-agents are used, the standard prompt template (§13.3) applies.

Earlier convention was "one PM prompt per sprint, dispatching all tracks as parallel sub-agents." That convention is retired:\ it overspent context budget on tasks that could be done sequentially in main session at lower cost.

### 13.3 Sub-Agent Prompt Template (opt-in)

Used only when the PI directs sub-agent dispatch (see §13.2 update). Sub-agents read CLAUDE.md by default. Do not restate CLAUDE.md context inside the prompt. Keep prompts terse: target < 1500 words for diagnostics, < 2500 for implementation. The format below is the canonical shape; omit any line that's obvious from context.

```
TASK: [one sentence, the deliverable]
DECISION GATE: [what counts as GO / BORDERLINE / STOP, with thresholds]
FILES TO READ (beyond CLAUDE.md): [only files the agent wouldn't naturally find]
DO NOT MODIFY: [only if a non-obvious file is at risk]
OUTPUT:
  - [what to return in the response, brief]
  - [what files to write: debug/*.py drivers, debug/data/*.json, debug/*_memo.md]
  - [paper edits to apply, if any — apply directly per §13.8]
```

Constraints (failed approaches to avoid, structures to preserve, success criteria) are stated only when non-obvious from CLAUDE.md §3 and the paper guardrails. The agent has access to all of CLAUDE.md; trust it to find what it needs.

### 13.4 Verification Gates

Before the PM agent accepts a sub-agent result, it checks:

1. **Test gate:** Do all relevant tests pass? (Non-negotiable.)
2. **Dead-end gate:** Does the approach match any entry in the failed approaches table? If so, reject unless the sub-agent explicitly explains what is different this time.
3. **Prime directive gate:** Does the result modify any discrete structure — quantum number labeling, selection rules, channel structure, Gaunt integral coupling? If so, do NOT accept. Escalate to plan mode for human review.
4. **Consistency gate:** Does the result contradict any established result in the papers? If uncertain, flag in the session summary rather than accepting.
5. **Equation gate:** Does every equation in affected papers have a corresponding test in `tests/`? If a new equation was added to a paper in this session, was a verification test also added? If not, the session summary must flag this as an open item.

**13.4a Equation verification (moved).** **No equation enters a paper without a test that verifies it computationally.** Verification types, per-claim sufficiency bar, `test_paper{N}_*` convention: `docs/multi_agent_protocol.md`.

### 13.5 Hard Prohibitions

The following changes must NEVER be made by sub-agents or the PM agent:

- Any change to the natural geometry hierarchy (new levels, changed coordinate systems)
- Introduction of any fitted or empirical parameter
- Deletion or suppression of negative results from Section 3 or CHANGELOG.md
- Presenting the **combination rule K = π(B + F − Δ) (Paper 2)** as anything stronger than an **Observation**. The standing label is Observation — never "conjecture", "derived", "prediction", "theorem", or any stronger tier. The three ingredients B, F, Δ have independent derived spectral homes; their *combination* is a numerical coincidence with no first-principles derivation (12 single-mechanism derivations eliminated). (History: Conjectures → Core 2026-04-18 Sprint A; Core → Observations 2026-05-02 per the curve-fit audit memo `docs/curve_fit_audit_memo.md`; conjecture → observation label downgrade 2026-06-14 per PI direction — "conjecture" was judged to carry unearned confidence that a derivation exists. Applies at the combination-rule level regardless of Paper 2's folder.)

**CLAUDE.md access control:** The PM may edit CLAUDE.md for mechanical updates that keep documentation in sync with code and paper changes. The PM may NOT edit sections that define strategy, framing, or the PM's own operating rules.

| Section | PM may edit? | Examples of allowed edits |
|:--------|:-------------|:-------------------------|
| 1 (Project Identity) | Version number ONLY | Bump v2.0.25 → v2.0.26 |
| 1.5 (Positioning & Framing) | **NO** | — |
| 1.6 (Project Phase) | **NO** | — |
| 2 (Development Frontier) | Yes | Update best results, add/complete track summaries, update backlog |
| 3 (Failed Approaches) | Yes (append only) | Add new failed approach rows; never delete or modify existing entries |
| 4 (Dimensionless Vacuum) | **NO** | — |
| 5 (Natural Geometry Hierarchy) | **NO** | — |
| 6 (Paper Series) | Yes | Update file paths, loading guide descriptions, inventory tables, key results |
| 7-9 (Code/Coding/Workflow) | Yes | Add new entry points, update module paths |
| 10 (Validation Benchmarks) | Yes | Add new benchmark rows for new tests |
| 11 (Topic-to-Paper Lookup) | Yes | Add new topic → paper mappings |
| 12 (Algebraic Registry) | Yes | Update status (algebraic-pending → algebraic) when proven |
| 13 (Multi-Agent Protocol) | **NO** | — |
| 14 (Test Architecture) | **NO** | — |

### 13.6 Track Management

Active work is organized into tracks. Each track has:

- A name and one-sentence goal
- A list of relevant papers and code modules
- A current status (active / blocked / complete)
- A log of sub-agent dispatches and results (maintained by the PM in `debug/track_logs/`)

The PM agent maintains a brief track status file at `debug/track_logs/STATUS.md` that is updated at the end of each session. This file is read at the start of the next PM session to restore context.

### 13.7 General Guidance

**All changes are autonomous.** The test suite catches code errors. The PI catches framing errors in review. The PM's job is to make its best judgment, apply changes directly (including to papers), and produce a clear summary of what changed and what the results were.

The overriding principles are: **code changes are autonomous** (the test suite catches errors); **paper updates are autonomous** (adding results, tables, new-method subsections, correcting claims contradicted by new evidence, reframing based on findings). If a change is wrong, it will be caught in plan-mode review and corrected — this is faster than proposal cycles.

### 13.8 Paper Update Policy

Papers are the authoritative source for all physics (Section 1). Code that outpaces the papers creates documentation drift. PMs are expected to keep papers in sync with code results.

#### Paper edit policy

PMs may edit papers in any of the six group folders (`papers/group1_operator_algebras/`, `papers/group2_quantum_chemistry/`, `papers/group3_foundations/`, `papers/group4_quantum_computing/`, `papers/group5_qed_gauge/`, `papers/group6_precision_observations/`), the synthesis folder (`papers/synthesis/`), and the archive folder (`papers/archive/`) directly for ALL changes, including:

- Adding or updating benchmark tables and numerical results
- Adding new subsections documenting methods and results
- Correcting claims that are contradicted by new computational evidence
- Reframing results based on new findings (e.g., reclassifying transcendental → algebraic when proven)
- Updating abstracts and conclusions to reflect current best results

**Record cross-paper dependencies WHEN YOU WRITE THEM (hard, added 2026-09-04, PI direction).** When a paper claim you write or edit *rests on* a claim owned by another document, record that dependency in the same edit:

- if the owner claim already has a `check_retracted_terms.py` entry, add your document to its `cited_by` (see §9, retraction → dependents rule);
- if it does not — the common case, since most owner claims have never been retracted — note it in your claim's `docs/claim_test_matrix.md` row as `rests on: Paper N's <claim>`. Whoever later retracts that claim greps for it when building the entry.

*The operational test, because not every citation is a dependency:* **if the cited claim were withdrawn tomorrow, would this sentence have to change?** If yes, it is a dependency and gets recorded. If the citation is context, provenance, or courtesy, it is not.

*Why this sits at authoring time rather than retraction time:* at retraction time the dependency set has to be **reconstructed from memory**, and reconstruction is what failed — eight of the ten recurring defect classes in the v5.4.4..v5.7.3 arc were a claim corrected in its owner and left standing in its citers, including one corrected in Paper 0 and left false in Paper 7 *in the same commit*. When you are writing the dependent sentence you are looking straight at the thing you are relying on — that is the only moment the set is free to obtain. Patterns cannot substitute: a citing document restates the claim in its own words, which is why four successive pattern rebuilds each missed loci and one sweep reported clean while five loci survived a spelling difference.

**PMs must still NOT:**
- Introduce fitted or empirical parameters without PI direction
- Change the natural geometry hierarchy (new levels, changed coordinates)
- Delete or suppress negative results from Section 3
- Present the **combination rule K = π(B + F − Δ) in Paper 2** as anything stronger than an **Observation** (the standing label since the 2026-06-14 conjecture→observation downgrade; never "conjecture", "derived", or "theorem"). The prohibition is at the combination-rule level, not the paper-tier level: it applies regardless of whether Paper 2 sits in Conjectures, Core, or Observations. See §13.5.

**Splinter file prohibition:** PMs must edit papers in-place. Do NOT create separate .tex diff files, proposal files, or draft directories. Proposed changes go directly into the paper. If the change is wrong, `git revert` is cheaper than context-loading splinter files in plan-mode review.

### 13.9 Session Summary Format (MANDATORY)

Every PM session MUST end with a summary in the following format. This is the PI's primary review mechanism — without it, changes are invisible.

```
## Session Summary [date]

### Tracks
- Track XX: [status] — [one-line result]

### Results
[Tables of computed data, benchmark comparisons, scaling exponents,
convergence studies — whatever the track produced. Include numbers,
not just descriptions. The PI needs to see the data to evaluate
whether the result is correct.]

### Files Modified
- `path/to/file` — [one-line description of what changed]

### Files Created
- `path/to/file` — [one-line description]

### Decisions
- [Any changes to paper claims, framing, or theoretical arguments — briefly stated]
- [Any negative results or dead ends encountered]
```

The Results section goes before the file list because the data is what the PI reviews first. If a track produced no quantitative results (e.g., a pure documentation track), the Results section can be replaced with a brief description of what changed and why.

### 13.9a Two primitives for PI-desired behavior

Two distinct primitives for shaping how the PM behaves. They are NOT interchangeable, and using the wrong one is a frequent failure mode (caught and corrected 2026-05-26 in the v3.3.1 release):

**Standing rules — `memory/feedback_*.md`.** Behaviors the PM runs automatically whenever the trigger condition is met. The PI should NOT have to invoke these. Memory rules are loaded into every conversation and apply unconditionally. The body of each rule states (a) the trigger condition, (b) why the rule exists (incident or insight that produced it), and (c) how to apply it. Examples: `feedback_audit_numerical_claims.md` (any "X matches Y" claim triggers curve-fit audit); `feedback_diagnostic_before_engineering.md` (≥ 2 honest negatives on a wall triggers diagnostic-only sprint); `feedback_tag_transcendentals.md` (any transcendental appearance triggers Paper 18 + Mellin engine + Paper 34 classification).

**Slash commands — `.claude/commands/<name>.md`.** Actions the PI triggers at a specific moment chosen by the PI. Slash commands are NOT defaults — they are explicit invocations. Each command is a saved prompt that fires in current main-session context (no sub-agent dispatch, no CLAUDE.md re-load, no context multiplier).

**The distinction.** If the behavior should fire whenever a condition is met → memory rule. If the behavior fires at a specific PI-chosen moment → slash command. A behavior CAN be both:\ the memory rule is the primary mechanism (enforces the default), and the slash command is a **force-fire backup** for cases where the PM has failed to detect the trigger condition. **Putting a behavior behind a slash command alone makes it more optional, not less** — it shifts the responsibility for triggering to the PI, which is the opposite of what a discipline rule needs.

**13.9b Slash commands (moved).** Genuine triggers: `/aha`, `/sprint-close`, `/checkpoint`, `/qa`, `/walls`. Force-fire backups for standing memory rules: `/audit-claim`, `/diag`, `/transcendental-tag`. `/qa` and `/walls` are **PI-invoked only** — never self-triggered. Each command's own file is `.claude/commands/<name>.md`; the full table is in `docs/multi_agent_protocol.md`.

**13.10 Research-agent integration (moved).** Research agents **do not write code or modify files** — they propose, the PM executes, the tests verify. Invocation + handoff format: `docs/multi_agent_protocol.md`.

### 13.11 Content Discipline and Token Efficiency

CLAUDE.md is loaded into every PM session and every sub-agent dispatch. Its size is paid repeatedly and every edit invalidates the prompt cache. The rules below keep this cost bounded. They are hard rules, not preferences — they override the impulse to "be thorough" in CLAUDE.md, because thoroughness has homes other than CLAUDE.md.

**Where each kind of content lives:**

| Content | Home | Format |
|:--------|:-----|:-------|
| Sprint chronicle (what happened, in detail) | `CHANGELOG.md` | Full prose, no length limit |
| Sprint summary (CLAUDE.md §2 entry) | CLAUDE.md §2 | ≤ 30 words: name + date + verdict + memo path |
| Dead-end record, full | `docs/failed_approaches_ledger.md` | The unabridged row. **Append here first.** |
| Dead-end record, index (§3) | CLAUDE.md §3 | Category row, or name + count in the per-sprint index |
| Superseded §2 bullets | `docs/development_frontier_archive.md` | Verbatim, by compaction round |
| Superseded WH status text | `docs/wh_register_history.md` | Verbatim (§1.7 rule 9: replace, never append) |
| Sprint memo | `debug/*.md` | ONE canonical memo per sprint, ≤ 5000 words |
| Surviving structural findings | papers/group*/*.tex | The papers, not CLAUDE.md, are the permanent record; papers cite permanent records (CHANGELOG / the paper / tests), **never transient `debug/`** (§9 policy, C14) |
| Cross-session facts | `memory/*.md` | One-liner index entry in MEMORY.md, ≤ 200 chars |
| Active behavior rules | CLAUDE.md §13 | Short, dense, this section |

**Rules:**

1. **No synthesis memos.** One canonical memo per sprint. Cross-sprint synthesis lives in CHANGELOG.md or in a paper. The "comprehensive synthesis memo that supersedes earlier synthesis memos" pattern is forbidden.

2. **CLAUDE.md §2 entries are one-liners.** Format: `**Sprint NAME (YYYY-MM-DD):** Verdict in one sentence. See debug/MEMO.md.` Full sprint detail goes to CHANGELOG.md (`### Added` / `### Changed` per release entry). Existing multi-thousand-word §2 entries are technical debt to be compacted as touched.

3. **CLAUDE.md §3 dead-end rows are short.** Format: `Approach name (date) | count | One or two sentences: the lesson + path to memo.` Existing paragraph-length rows are technical debt.

4. **Memory files are for cross-session facts not derivable from CLAUDE.md or papers.** Do not auto-create memory files for sprint outcomes (sprint detail lives in CHANGELOG.md). MEMORY.md index entries strictly ≤ 200 chars — the system silently truncates beyond ~24KB and the warning has fired multiple times.

5. **Agent prompts: task + decision gate + specific files + output format.** Do not restate CLAUDE.md context inside agent prompts — sub-agents already load CLAUDE.md. Target: < 1500 words for diagnostic tasks, < 2500 words for implementation tasks. If the prompt is growing past 2500 words, the task is probably too big for one agent.

6. **Prefer Explore agent for read-only diagnostics.** Reserve general-purpose for tasks that require code modification or extensive computation. Explore has narrower context allocation.

7. **Don't dispatch sub-agents for tasks doable in the main session.** A single Read + Edit, a few-line script, or analysis under ~500 lines of context should be done by the PM directly. Sub-agents are for parallelizable work or context-heavy delegation.

8. **One canonical record per fact.** A fact appears in at most one of: CLAUDE.md, CHANGELOG.md, paper, memo, memory. Duplication is the failure mode that produced the current bloat.

9. **Status updates replace, never append (added 2026-06-10, PI-authorized).** When a WH status, a §6 status flag, or a paper-state description changes, REPLACE the existing text and move the superseded version to its history home (`docs/wh_register_history.md`, `docs/paper_notes_archive.md`, CHANGELOG.md). Chronicling-by-appending inside CLAUDE.md is the failure mode that regrew the file from 1,263 lines (2026-05-31 compaction) to 1,400 lines / 320 KB by 2026-06-10.

10. **Compaction is relocation, never deletion (2026-09-01).** Content moves *verbatim* to its archive; the section keeps a pointer and enough stub that the item stays discoverable by name. The 2026-09-01 round took the file 220 KB → 99.7 KB and was verified by a conservation check (all 123 §3 rows verbatim in the ledger, all 121 dead ends still named, all 100 §2 bullets kept or archived). **Run that check on any future compaction** — one that silently drops content is worse than a large file.

11. **Apply rules 2–3 when writing the entry, not later.** At that round, 85 of 100 §2 bullets were over the 30-word budget and 76 were dated to a *single month* — never an age problem, just a rule never enforced in bulk.

**Enforcement.** When the PM is tempted to write a long CLAUDE.md §2 bullet, a synthesis memo, or a verbose agent prompt — stop, move the content to its proper home, and write the short version in CLAUDE.md. Apply this to existing entries when touched, not as a one-time pass.

---

## 14. Test Architecture Policy

### Main pipeline tests (`tests/`)

The `tests/` directory contains ONLY tests that validate the main pipeline — the methods and results described in the active papers (Papers 0–18, FCI-A, FCI-M). Every test file in `tests/` should correspond to at least one paper result or a piece of infrastructure that the pipeline depends on.

**The rule:** If a paper documents a method as its primary result, that method's tests live in `tests/`. If a method is documented as failed (Section 3), superseded (architecture locked), or was a one-off diagnostic, its tests live in `tests/_archive/`.

### Archived tests (`tests/_archive/`)

Archived tests are organized into three tiers:

| Directory | Contents | When to use |
|-----------|----------|-------------|
| `tests/_archive/dead_ends/` | Tests for approaches documented as failed in Section 3 | Method confirmed dead; resolution documented |
| `tests/_archive/superseded/` | Tests for approaches replaced by a newer architecture | New method has its own tests in `tests/`; old method explicitly "architecture locked" |
| `tests/_archive/diagnostics/` | Tests from one-off investigation arcs | Diagnostic complete; resolution baked into main pipeline |

Archived tests are **not** collected by default (`--ignore=tests/_archive` in pytest config). They can be run explicitly via `pytest tests/_archive/`.

Archived tests have scientific value as institutional memory — they help prevent re-deriving failed approaches. **Do not delete archived tests.**

### Archived code (`geovac/_archive/`)

Follows the same three-tier structure (`dead_ends/`, `superseded/`, `auxiliary/`). Archived code modules are still importable from their `geovac._archive.*` paths. If a core module previously imported an archived module, a compatibility shim is left at the original location.

### When new work produces dead ends

When a track produces a negative result or a method is superseded:

1. Document the failure or supersession in CLAUDE.md Section 3 (failed approaches) or Section 2 (active frontier, marking the track COMPLETE)
2. Move the relevant tests to `tests/_archive/` in the same commit
3. Move the relevant code module to `geovac/_archive/` in the same commit
4. **Check for redirects first:** If any tests validate infrastructure that's still live (CI machinery, solver constructors, shared utilities), extract those tests and rewrite them to call the current architecture before archiving the rest
5. Update import paths in archived files

### Slow tests

Core tests that are computationally expensive are marked with `@pytest.mark.slow` and skipped by default. Run them with `pytest --slow`. A test qualifies as slow if it takes >10 seconds or performs heavy computation not needed for fast regression (e.g., cc-pVTZ integral engine, full VQE pipeline).

### Redirect-before-archive rule

Before archiving a test file, the PM must check: does this file test any infrastructure that the main pipeline still depends on? If yes, extract those tests into a core test file that calls the current API, then archive the rest. This prevents coverage gaps from forming when methods are superseded.

### Post-archive test commands

| What You Want | Command |
|---------------|---------|
| Fast regression (daily dev) | `pytest` |
| Full core + slow tests | `pytest --slow` |
| Archived tests only | `pytest tests/_archive/` |
| Everything | `pytest tests/ tests/_archive/ --slow` |
| Specific archive tier | `pytest tests/_archive/dead_ends/` |
| Topological proofs only | `pytest tests/test_fock_projection.py tests/test_fock_laplacian.py` |

---

## 15. Numeric Registry (the `\gvq` scheme)

**Load-bearing numbers in the papers are linked objects, not free text.**

- **Artifact:** `debug/qa/numeric_registry.py` — every load-bearing quantity
  with its canonical value, the **convention** it is stated in, its
  **provenance** (measured by what route / cited from what source), and, if
  derived, the **expression** rather than a stored number.
- **Citation form:** `\gvq{registry_key}{literal}` in the `.tex`. The macro
  renders the literal only, so the PDF stays self-contained and archivable
  and the `.tex` stays hand-editable; the key is a foreign key.
- **Gate:** C21 (`debug/qa/check_numeric_consistency.py`) recomputes every
  derivation, checks every annotation against the registry (accepting
  correct display rounding), blocks retired values, and reports unregistered
  multi-document numerals.
- **Query:** `--index <key>` names every locus citing a quantity, plus the
  derivations it feeds.

### Why it exists

A census of the group4 + group6 numeric surface found **352 numerals at more
than one locus, 201 across more than one document**. Successive QA delta
passes were each finding 15–20 genuine defects — about **5% of that coupled
surface per pass** — because manual review *samples* the dependency graph
rather than traversing it. The classes that kept surviving were relational,
not value-level: **twins** (the same quantity tabulated in two papers),
**derivations** (λ/Q, N², ratios, exponents fitted from printed columns), and
**conventions** (identity-in vs identity-out; a fit's point range).

### Operating rules

1. **A value moves → edit the registry FIRST, then run C21.** The gate names
   every stale locus, including in other papers. Do not hand-sweep; that is
   the process this replaces.
2. **Then re-read the prose at each locus C21 names.** This is the step the
   scheme exists to force. A sentence around a number is a claim keyed to
   that number's *magnitude* and it does not update itself — "far below" was
   true at 5.8× and false at 0.92×; "the d-block is the sparsest" reversed
   outright. Updating a numeral without re-reading its sentence is the
   defect, not the fix.
3. **Never register a value you have not measured or cited.** `provenance`
   is not decoration; a registry of guesses launders them into authority.
4. **Annotation must not change rendered text.** Wrapping a literal is a
   no-op by construction; verify it (strip all `\gvq` and diff) rather than
   assuming it.
5. **Do not resolve values at build time.** Deliberately rejected: it would
   update numerals while leaving the surrounding claims silently false, and
   it would break the self-containedness of DOI-stamped PDFs.

Full rationale: the `numeric_registry.py` docstring and `docs/qa/criteria.md`
(C21). Self-test: `tests/test_numeric_registry.py`.
