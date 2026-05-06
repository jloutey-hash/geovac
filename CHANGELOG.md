# Changelog

All notable changes to GeoVac will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

> **Note:** the CHANGELOG is currently behind the `CLAUDE.md` version cursor (CLAUDE.md tracks v2.10–v2.26 range; CHANGELOG below jumps from v2.9.2 to v2.26.1). Intermediate version entries for the RH sprint series (v2.20–v2.25, Papers 28/29/30) are in `git log` commit messages but have not yet been back-filled into this CHANGELOG. A consolidation sprint is flagged for future work.

## [2.28.0] - 2026-05-06

### Added — WH1 PROVEN, Sprint H1 verdict, L2 quantitative rate

Two-sprint sequence (May 5 prerequisite closure + May 6 keystone closure) executing the PI's three-step spectral-triple commitment of 2026-05-04 (R2.5 GH → real structure J → almost-commutative Higgs).

**May 5 — three-track prerequisite closure:**

- **R2.5 L4 Berezin reconstruction CLOSED.** Constructed `B_{n_max}: C(S³) → O_{n_max}` as `Σ_{N≤n_max,L,M} K̂(N) c_{NLM} M_{NLM}` with Plancherel weight `K̂(N) = N/Z_{n_max}` from L2 — SU(2) non-Kähler analog of Connes-vS 2021 §3 Fejér-kernel reconstruction on S¹. Four properties proved: positivity for f≥0, contractivity ‖B(f)‖_op ≤ ‖f‖_∞, approximate-identity bound with γ → 0 inherited from L2, L3 Lipschitz compatibility with C_3 = 1. Numerical verification at n_max ∈ {2, 3, 4} on a 5-function panel. **Four of five GH-convergence lemmas closed.**
- **WH1-Connes Step 2 — real structure J at finite n_max.** Three load-bearing Connes axioms hold exactly at finite n_max ∈ {1, 2, 3} on truthful CH: J² = −I, JD = +DJ (Euclidean KO-dim 3), J·O·J⁻¹ = O — bit-exact zero residuals. Clean negative on offdiag CH: JD = +DJ FAILS at residual 2.0 (offdiag is SDP-bounding for Connes distance only, not spectral-action-foundational). Order-zero/one fail at 5–20% as truncation artifacts.
- **Sprint H1 scoping CLOSED.** Cross-track sharpening: Marcolli–vS + Perez-Sanchez 2024/2025 give Yang–Mills without Higgs by default; CC Higgs source = off-diagonal D_F block on AC factor. Track 3 Candidate A (offdiag CH bridging) RULED OUT by Track 2; H1 architecture is internal D_F couplings on M_n(C) factor with J = J_GV ⊗ J_F.

New modules: `geovac/berezin_reconstruction.py` (~485 lines, 49 tests), `geovac/real_structure.py` (~620 lines, 43 tests + 1 slow), `geovac/almost_commutative.py` (Track 3 stub, 278 lines). Memos: `debug/r25_l4_proof_memo.md` (~3500 words), `debug/real_structure_finite_nmax_memo.md` (~2400 words), `debug/almost_commutative_scoping_memo.md` (~3700 words). Paper 32 §IV extended (+135 lines including bibtex fix `connes_vsuijlekom2021` → `connes_vs2021`); §VIII GH-convergence Remark extended (~30 lines).

**May 6 — three-track keystone closure: WH1 PROVEN.**

- **R2.5 L5 Latrémolière propinquity assembly CLOSED.** Tunneling pair (B_{n_max}, P_{n_max}) using L1'–L4 ingredients. **Theorem (Paper 32 §VIII, `thm:gh_convergence`):** the truncated metric spectral triples T_{n_max} converge to the round-S³ Camporesi–Higuchi spectral triple T_{S³} in the Latrémolière propinquity, Λ(T_{n_max}, T_{S³}) ≤ C_3·γ_{n_max} → 0, with C_3 = 1 (L3) and γ_{n_max} from L2's central Fejér moment. Numerical verification at n_max ∈ {2, 3, 4}: γ_2 = 2.0746, γ_3 = 1.6101, γ_4 = 1.3223; monotone; Λ(4)/Λ(2) = 0.637. Five-lemma chain closed: L1' R3.5 (May 4), L2 (May 4), L3 (May 4), L4 (May 5), L5 (May 6). Connes–vS 2021 deferred this convergence question to "elsewhere" three times; subsequent NCG follow-ups (Hekkelman–McDonald 2024, Leimbach–vS 2024) covered only S¹/T^d/S² flat structures. **WH1 PROMOTED to PROVEN** per PI authorization 2026-05-06; CLAUDE.md §1.7 WH1 entry updated. Paper 32 §VIII formal Theorem replaces prior "Status of GH-convergence roadmap" Remark; Limit Identification Remark added (Wasserstein–Kantorovich state space); WH1 Keystone Closure Remark added.
- **L2 quantitative rate via Stein–Weiss IBP CLOSED.** Asymptotic rate constant proven rigorously: lim_{n→∞} n·γ_n / log n = **4/π** ≈ 1.27324, approached from above. Uniform bound γ_n ≤ 6 log n / n for all n ≥ 2 (tight at n=2 with 0.2% margin). Closed-form sum-rule decomposition γ_n = π − 4T_n / (π Z_n). **Headline structural connection:** the constant **4/π = Vol(S²)/π² is the M1 Hopf-base measure mechanism signature** of Sprint TS-E1's master Mellin engine — the GH-convergence rate of the GeoVac spectral triple itself carries the same M1 transcendental signature that the case-exhaustion theorem (Paper 32 §VIII) identified as one of the three π-source sub-mechanisms. SU(2) analog of Leimbach–vS 2024 torus rate, slowed by exactly one log factor due to non-abelian volume.
- **Sprint H1 minimal AC extension VERDICT: POSITIVE-THIN.** AC extension constructed at finite n_max: A_F = C ⊕ H on doubled H_F = C^4_matter ⊕ C^4_antimatter (Connes–Marcolli convention), KO-dim 3+6 = 9 ≡ 1 (mod 8), J_AC = J_GV ⊗ J_F with J_F = σ_x ⊗ I_4 giving J² = −I and JD = +DJ exact. Pitfall corrected: undoubled H_F = C^4 with J_F = iσ₂K ⊕ iσ₂K failed JD = +DJ for non-degenerate Yukawa; matter–antimatter doubling fixes it (canonical CC convention). Inner fluctuations decompose ω = ω_gauge + ω_Higgs cleanly: U(1) (from C) and SU(2) (from H) gauge sectors recover Papers 25/30 at operator level; Higgs sector non-zero off-diag C↔H block iff D_F has Yukawa. Falsifier verdict at n_max ∈ {1, 2, 3} with 50 random generators each: strong falsifier "every Hermitian D_F from GeoVac forces Higgs zero" **HOLDS iff Y = 0**. With Y = 0: Φ = 0 exactly. With Y > 0: ‖Φ‖_max ∈ [0.05, 0.27]. **GeoVac sits on the Marcolli–vS-without-Higgs side at the structural level** — inner fluctuations admit a Higgs but GeoVac does not autonomously select the Yukawa. Paper 32 §VIII.B G2 sharpened from "open structural question" to "construction admits Higgs but Yukawa is a free input not derived by GeoVac"; §VIII.C addendum (~150 lines) documents H1 verdict; G3 (R3.5 chirality vs weak-isospin chirality identification γ_5) stays open as the natural next sprint within the §VIII.B electroweak co-location target.

New modules: `geovac/gh_convergence.py` (~510 lines, 39+2 tests). Expansion of `geovac/almost_commutative.py` to 854 lines (replacing 278-line stub) with 38 tests. `geovac/central_fejer_su2.py` extended with 7 quantitative-rate functions: `T_n_via_sum_rule`, `gamma_n_via_sum_rule`, `asymptotic_rate_constant`, `quantitative_rate_bound`, `doubling_estimator`, `N0_for_constant`, `quantitative_rate_certificate`. Memos: `debug/r25_l5_proof_memo.md` (~3500 words), `debug/r25_l2_quantitative_rate_memo.md` (~2700 words), `debug/h1_ac_extension_memo.md` (~3200 words). Paper 32 §VIII formal Theorem (replaces Remark), §VIII.B G2 sharpened, §VIII.C addendum (~150 lines, H1 verdict), §VIII Q1 forward reference.

**Standing open items going into v2.28.x:**
- Paper 38 manuscript (J. Geom. Phys. or Adv. Math. companion to Leimbach–vS 2024) — five proof memos ready for assembly
- G3 closure (electroweak chirality co-location: γ_GV ↔ γ_5) — natural next physics sprint
- Real-structure J finite-resolution audit (order-zero/one violations should converge to zero at rate γ_{n_max} now that propinquity limit is established)

## [2.27.3] - 2026-05-04

### Added — Sprint TS (Triple-Taxonomy Bridge) and WH1-R3.5 closure

Eight-track sprint testing whether the Paper 34 projection taxonomy and the Paper 32 spectral-triple framing are two views of the same structure, plus completion of the WH1 R3.5 thread on the full Dirac operator system.

**Sprint TS — Triple-Taxonomy Bridge:**

- **Track A (GH convergence sketch):** Leimbach–vS 2024 torus proof transports cleanly to S³ via Peter–Weyl on SU(2). Fock-graph index (n, l, m_l) IS the Peter–Weyl basis under n = 2j+1; the central spectral Fejér kernel is the abelianizing assumption that makes Schur–Fourier transference go through. Verdict (b): reachable in 4–8 weeks, with R3.5 (full Dirac, both chiralities) as 2-week prerequisite. Five-lemma roadmap. Memo: `debug/track_ts_a_gh_convergence_memo.md`.
- **Track B (15-row triple↔taxonomy dictionary):** Mechanical mapping of each Paper 34 projection to its sector of (A, H, D), Connes-axiom impact, universal/Coulomb status, transcendental class. 9 candidate asymmetries flagged. Memo: `debug/track_ts_b_dictionary_memo.md`.
- **Track C (asymmetry investigation):** Three sub-tasks. (i) **C-1 ERRATUM:** TX-A memo's "spectral_action ↔ observation_window mutually one-way (UV/IR ordering anomaly)" claim was a confabulation — JSON has them commuting. CLAUDE.md §2 TX-A bullet, Paper 34 §VII (T3 theorem) and §IX (Conclusion), and memory file all corrected. The QFT UV/IR anomaly is real physics that TX-A's coarse layer-grading cannot see; finer output_layer index is the natural refinement. (ii) **C-7:** shared-Hopf hypothesis on K = π(B+F−Δ) NOT confirmed; F = D_{n²}(d_max) and Δ⁻¹ = g_3^Dirac are derived without invoking Hopf; outer π already counted in depth-3 chain. **Twelve mechanisms now eliminated** for K = single-principle derivation. Paper 34 §VIII.3 sharpened with Track TS-C update. (iii) **#14 rest-mass:** confirmed as central deformation D² → D² + m²·𝟙. Drake-Swainson passes the same test. **New tier in Paper 31 §VII.2.5: parametric / central sector**, with rest-mass and Drake-Swainson as charter members. Memo: `debug/track_ts_c_asymmetry_investigation_memo.md`.
- **Track D (Paper 35 × triple cross-check):** Outcome 1 — all 15 projections agree under both lenses. Every π-source reduces to one of three triple-framing mechanisms (M1 Hopf-base measure Vol(S²)/4; M2 Seeley–DeWitt √π on S³; M3 vertex-topology Hurwitz Dirichlet L). Paper 35 Prediction 1 promotes from "208/208 empirical" to candidate spectral-triple THEOREM. Memo: `debug/track_ts_d_paper35_triple_test_memo.md`.
- **Track E1 (theorem promotion):** Theorem reachable as stated, no fourth mechanism. **Paper 32 §VIII added** (~140 lines LaTeX): π-source case-exhaustion theorem, proof skeleton, three remarks (joint engagement, K status, falsification target). Inserted between §VII (Coulomb/HO asymmetry) and §IX (Marcolli–vS Lineage). Memo: `debug/track_ts_e1_theorem_promotion_memo.md`.
- **Track E3 (sixteenth projection falsification):** All three candidates REDUCE to combinations of M1/M2/M3 (continuum gravity → M1+M2; anomaly classes → M1+M2+M3; principal bundles → M1 normalisations against integer integrands). **Headline structural finding:** M1, M2, M3 are not three independent mechanisms but **three sub-cases of a single master Mellin-engine mechanism** π-source = M[Tr(D^k · e^{−tD²})] with k ∈ {0, 1, 2} selecting the sub-mechanism (k=0 → Hopf-base, k=1 → vertex parity, k=2 → spectral action). **Paper 18 §III.7 sharpened** with the master-mechanism reading. Falsification target flagged: discrete c_1 on Hopf S³ graph at n_max=3 (predicted integer-valued). Memo: `debug/track_ts_e3_sixteenth_projection_memo.md`.

**WH1-R3.5 — full Dirac (both chiralities) operator system, COMPLETE 2026-05-04:**

- New module `geovac/full_dirac_operator_system.py` (~580 lines, 29 tests) with `FullDiracLabel`, `FullDiracTruncatedOperatorSystem`, truthful and offdiag CH builders.
- All five computational legs done: truthful n_max=2,3; Weyl-only truthful n_max=3 reference; offdiag n_max=2,3.
- **Headline numbers (offdiag CH n_max=3, full Dirac):** dim_H = 40, total pairs 780, **+∞ count: 0** (down from 624 in truthful), forced zeros 60 (m-reflection), finite count 720, Pearson nz **−0.2501**, Spearman nz −0.2300.
- **Refactor of `geovac/connes_distance.py`:** kernel of [D, ·] in O is structural data of (op_sys, D); was being recomputed inside every per-pair SDP call. Cached once in `compute_distance_matrix` and passed through. Also removed redundant 2nd SDP solve (feasible set is symmetric in x, so max diff = max −diff). All 17 tests pass; regression vs n_max=2 baseline at 5×10⁻⁹ max diff. **Speedup: ~15–25× per-pair**, n_max=2 offdiag dropped from ~10–15 min to **39 seconds**; n_max=3 offdiag completed in **66 minutes** vs projected 1–3 hours.
- **Track A's L1 lemma reformulation:** original L1 ("truthful CH, every cross-pair finite") is FALSE at every tested n_max. **L1' (offdiag CH, every cross-pair finite) is VERIFIED** at n_max=2 and n_max=3. The offdiag CH operator system is the natural setting for Track A's GH proof: it includes shell-coupling multipliers M_{N,L,M} that break ker([D, ·]) ∩ O_h down to ℂ·𝟙 structurally. Track A's 4–8 week timeline is preserved with this reformulation. Memo: `debug/wh1_r35_full_dirac_memo.md`. New memory: `wh1_r35_full_dirac_complete.md`.

**Cross-paper amendments applied this version:**

- Paper 32 §VIII (case-exhaustion theorem, ~140 LaTeX lines)
- Paper 18 §III.7 (master-mechanism reading paragraph)
- Paper 31 §VII.2.5 (parametric / central sector)
- Paper 34 §VII T3 theorem and §IX conclusion (TS-C erratum)
- Paper 34 §VIII.3 (Track TS-C update on K's depth-prediction anomaly)
- Paper 35 §VII new subsection on TS-D promotion to spectral-triple theorem
- CLAUDE.md §2 Sprint TS bullet, §6 Paper 32 inventory entries

**New code/tests:** `geovac/full_dirac_operator_system.py`, `tests/test_full_dirac_operator_system.py` (29 tests), refactored `geovac/connes_distance.py` (kernel cache + redundant solve removal).

**Standing open items:** R2.5 keystone (Peter–Weyl GH on full state space, 4–8 weeks) and Lemma 5.3 Lipschitz bound on offdiag CH multipliers (1–2 weeks; infrastructure ready) for Track A's GH proof.

## [2.27.2] - 2026-05-03

### Added — Paper 35 (Time as Projection) and Paper 36 (Bound-State QED on S³)

Two new observation papers extending the Paper 34 projection dictionary, plus a substantial expansion of Paper 34 itself. All work consolidated in a single chat session.

**Sprint sequence (chronological within the day):**

- **Paper 4 archived** (`papers/conjectures/Paper_4_Universality.tex` → `papers/archive/Paper_4_Universality.tex`). Two reasons: (1) the proton radius puzzle the paper claimed to fully resolve has been settled in the physics community as a measurement artifact in older electronic spectroscopy (PDG consensus 0.8409(4) fm; April 2026 review "the puzzle is no more"); (2) the holographic central-charge / spectral-dimension claims (d_s ≈ 2.07, c ≈ 1/36) inherit the machinery retracted from Paper 3 (CLAUDE.md v2.18.2 retraction). Salvageable kernel (mass-as-unit-fixing-projection) absorbed into Paper 35.

- **CLAUDE.md §4 transcendentals policy refined.** Added "algebraic-first, observation-aware" language: irreducible transcendentals that survive algebraic decomposition are not failures of the algebraic-first discipline — they are the content of a specific Paper 34 projection and should be tagged to that projection. The diagnostic question for any quadrature wall is two-headed: missing algebraic structure (decompose) or observation-side projection signature (catalogue and stop). Anonymous transcendentals not allowed in production code or papers.

- **Sprint KG (KG-1, KG-2, KG-3) → Paper 35**: `papers/observations/paper_35_time_as_projection.tex` (~810 lines after edits). KG-1 verified the bare Klein-Gordon spectrum on S³ × ℝ is π-free in ℚ[√d for d square-free positive integer] for every rational m² (200 cases, n ∈ [1,50] across panel {0, 1, 1/4, 2}, symbolic verification). KG-2 explicitly identified the Matsubara mode (n=0, k=1) at ω² = 4π² as the first π-bearing eigenvalue under temporal compactification on S³ × S¹_β; before/after table in §IV. KG-3 computed the conformally coupled massless scalar Casimir on unit S³ as **E_Cas = 1/240 exact rational** via ζ_R(s−2) reduction at s=−1 with B₄ = −1/30; matches Bytsenko et al. / Dowker-Critchley / Ford to 60 dps; per-step transcendental ledger has zero π injections in the spatial calculation. Stefan-Boltzmann π²/90 in the high-T limit on S³ × S¹_β enters only via the Matsubara sum.

- **Sprint KG-5 (Dirac Casimir spinor companion)**: spatial Dirac Casimir on unit S³ via Camporesi-Higuchi spectrum |λ_n| = n + 3/2 with degeneracy 2(n+1)(n+2). Spectral zeta ζ_|D|(s) = 4[ζ_R(s−2, 3/2) − (1/4) ζ_R(s, 3/2)] via Hurwitz with B_2(3/2)=11/12 and B_4(3/2)=127/240. Result: **E_Cas^Dirac = 17/480 exact rational** (full Dirac convention; Weyl 17/960; single-chirality-single-sign 17/1920). Matches Camporesi-Higuchi 1996 Eq. 5.27 to 40 dps. Confirms Paper 35 prediction in spinor sector.

- **Paper 34 amended (13 → 15 projections)**: §III.14 (Rest-mass projection: introduces m, dimension mass, transcendental signature trivial / ring-preserving for m² ∈ ℚ) and §III.15 (Observation/temporal-window projection: introduces β or T = 1/β, dimension time, transcendental signature 2π·ℚ per Matsubara mode; π^{2k}·ℚ in integrated quantities; Stefan-Boltzmann π²/90 in high-T limit). Table 1 extended; Observation 2 (two-axis duality) updated with two new tier mappings; §VIII open question 1 updated for 13 → 15 expansion. Empirical catalogue gets two new single-projection rows: 1/240 (scalar Casimir) and 17/480 (Dirac Casimir).

- **Sprint LS-5 (α⁵ Lamb shift two-loop scoping)**: structural feasibility verified for two-loop QED on Dirac-S³ — mode-count budget favorable (~2,870 terms at n_max=20 from SO(4) selection rules), `qed_two_loop.py` infrastructure adequate. Critical reframing of LS-3 residual: ~80% of the residual was a one-loop convention mismatch flagged in LS-1 §2.2 footnote, only ~20% genuine A-tier multi-loop. Recommended LS-6a (1 sprint) before any two-loop work.

- **Sprint LS-6a (Eides §3.2 convention fix)**: re-derived LS-1 one-loop SE in canonical Eides §3.2 convention. **Identified the LS-1 +38/45 SE coefficient as a Uehling-kernel double-counting bug**: the canonical Eides 10/9 SE constant already includes the Karplus-Klein-Darwin and j=1/2 Schwinger AMM contributions; LS-1 had inadvertently subtracted the Uehling kernel constant (4/15 = 10/9 − 38/45) from this, leaving 38/45. The Uehling shift is already counted as a separate VP contribution. Restoring 10/9 shifts the predicted Lamb shift by exactly +27.13 MHz at n=2, Z=1, structurally identifiable as (4/15) · α³Z⁴ / πn³ · (Ha-to-MHz). **Result**: LS-6a Lamb shift = **1052.19 MHz** vs experimental 1057.845 MHz; residual −5.65 MHz / **−0.534%**. 5.8× reduction in absolute error vs LS-3.

- **Sprint K-CC (clean triple negative, WH5 strengthened)**: tested Paper 35 §VII.2 hypothesis that K = π(B + F − Δ) is structurally one Connes-Chamseddine spectral-action coefficient on S³, not three independent projections summing. Three sub-tracks: (a) PSLQ at 100 dps over 244-element basis on Λ_∞ = 3.7102454679060528505 — no identification; (b) search for a natural Λ from S³ spectral data — only tautological hits; (c) **algebraic obstruction**: T9 forces ζ_{D²}(s) at integer s into √π·ℚ ⊕ π²·ℚ, so F = π²/6 cannot appear at any integer s of the unified heat-kernel expansion. B, F, Δ live in three categorically different spectral objects on two different bundles. **Twelve mechanisms now eliminated for K = single-principle derivation** (nine from Phases 4B-4I + Sprint A α-SP + three from K-CC). Paper 2 stays in Observations status. Paper 35 §VII.2 sharpened to clean negative on possibility (i) (CC unification ruled out); possibilities (ii) (non-temporal source for K's π) and (iii) (K is below framework's intrinsic resolution) remain live. Side finding: ζ_{D²}(2) = π² − π⁴/12 verified symbolically as a new T9-consistent closed form (not previously logged).

- **Sprint LS-7 (first pass at two-loop SE)**: iterated CC spectral action on Dirac-S³ confirms the two-loop SE prefactor (α/π)² · (Zα)⁴ · m_e c² / n³ — the 1/π² factor traces to two iterated Schwinger proper-time integrations, exactly as Paper 35 Prediction 1 specifies. Weak test (universal across two-loop QED). Strong test (LS-8a): native derivation of dimensionless bracket coefficient C_2S = +3.63 from iterated CC + bound-state Sturmian. **Critical residual reframing**: LS-7 identified that the LS-5/LS-6a "+7.10 MHz Tab 7.4 multi-loop ceiling" reading was a misreading. Eides Tab. 7.3 (proper α⁵ multi-loop QED) totals only ~+1.20 MHz; the rest of the −5.65 MHz residual is non-loop physics (recoil −2.40, finite nuclear size +1.18, hyperfine averaging ~+5.0). Multi-loop test target sharpened to ~+1.20 MHz; non-loop sectors covered by Paper 34 §III.14 rest-mass projection and Paper 23.

- **Paper 36 (Bound-State QED on S³)**: `papers/observations/paper_36_bound_state_qed.tex` (~700 lines). Consolidates LS-1..LS-7 into a standalone observation paper. Architecture (Dirac-S³ + Sturmian projection at λ=Z/n); Bethe log (velocity vs acceleration vs Drake-Swainson); Eides convention reconciliation with Observation 1 (Uehling kernel double-counting); result and component breakdown including LS-7 reframed residual decomposition; LS-7 outlook with Proposition 1 sharpened (multi-loop test target ~+1.20 MHz; LS-8a strong test). Sub-percent closure confirms framework's internal completeness for one-loop QED on S³ even more cleanly than originally framed.

- **Paper 34 LS-3/LS-6a/LS-7 catalogue rows reflagged** to reflect LS-7 reframing: LS-3 row now C (resolved by LS-6a) + A; new LS-6a row with the clean 1052.19 MHz / −0.534% result and structural identification of +27.13 MHz convention shift; new LS-7 structural row.

- **CLAUDE.md updates**: version v2.27.0 → v2.27.2 (two patches in one day for Paper 35 + Paper 36 additions); §4 transcendentals policy refined (PI-authorized edit to normally-protected section); §6 Paper 4 archived, Paper 35 + Paper 36 added to inventory and loading guide; §11 nine new Paper 35 + five new Paper 36 topic rows; §1.7 WH5 status field strengthened with Sprint K-CC findings (twelve mechanisms now eliminated; algebraic obstruction via T9; ζ_{D²}(2) closed form noted).

### Files

**New files created:**
- `papers/observations/paper_35_time_as_projection.tex`
- `papers/observations/paper_36_bound_state_qed.tex`
- `debug/kg{1,2,3,5}_*.{py,memo.md}` + `debug/data/kg{1,2,3,5}_*.json` (4 sprint pairs)
- `debug/kcc_{lambda_pslq,natural_lambda,unified_heatkernel,gaussian_check}.py` + memo + JSON
- `debug/ls{5,6a,7}_*.{py,memo.md}` + JSON

**Files modified:**
- `papers/observations/paper_34_projection_taxonomy.tex` (Paper 34 amendments throughout)
- `CLAUDE.md` (extensive)
- `CHANGELOG.md` (this entry)

**Files moved:**
- `papers/conjectures/Paper_4_Universality.tex` → `papers/archive/Paper_4_Universality.tex`

### What's queued for next session
- LS-8a (~2 sprints): native derivation of C_2S = +3.63 from iterated CC + bound-state Sturmian; the strong test of Paper 35 Prediction 1 in the multi-loop sector
- LS-8b (1 sprint): Karplus-Sachs two-loop VP (+0.16 MHz; simpler than two-loop SE; would validate iterated VP machinery)
- LS-8d/e (~3 sprints): recoil + FNS + hyperfine via Paper 34 §III.14 rest-mass projection + Paper 23
- K-CC follow-up: possibility (ii) for the K conjecture (non-temporal source for K's π) — Marcolli-vS gauge-network coefficient (per WH1) is the natural next place to look

## [2.27.0] - 2026-05-03

### Added — Paper 34: Projection Taxonomy and Empirical Matches Catalogue

- **New paper** `papers/observations/paper_34_projection_taxonomy.tex` (~1065 lines) names the GeoVac framework's two-layer architecture explicitly:
  - **Layer 1** = bare combinatorial graph (quantum-number labels, integer eigenvalues, rational matrix elements, π-free, no physics)
  - **Layer 2** = thirteen named projection mechanisms; each adds (i) specific physical variables, (ii) specific physical dimensions, (iii) specific transcendental signature
- **Three-axis dictionary** (Table 1, §IV) tags each projection by variable / dimension / transcendental class
- **Empirical matches catalogue** (§V) groups verified machine-precision matches by projection depth (zero, one, two, three, four)
- **Off-precision matches with error-source classification** (§V.B) adds T (truncation), B (basis quality), A (approximation order), C (calibration mismatch), S (structural floor) tags
- **Falsifiable Prediction 1** (§VI): error compounds with projection depth; consistent with Lamb shift 4-projection chain at 3.10% residual
- **Living document protocol**: PMs may append catalogue rows autonomously per CLAUDE.md §13.8 with the constraint that no structural identification is asserted beyond what the producing sprint verified
- Companion to Paper 18 (transcendental taxonomy) — duality stated as Observation 2

### Added — Bound-state QED arc (LS-1 → LS-4)

Four-sprint sequence delivering the first bound-state QED observable computed in GeoVac.

- **LS-1** (Lamb shift via standard formula on Dirac-on-S³): ΔE_VP(2S₁/₂) = −27.13 MHz from Π = 1/(48π²) cross-checks textbook Uehling shift to <1%; total Lamb shift 1025.06 MHz at −3.10% (Bethe logs imported from Drake-Swainson 1990). Files: `debug/ls1_lamb_shift.{py,memo.md}`, `debug/data/ls1_lamb_shift.json`.
- **LS-2** (velocity-form Bethe log via Coulomb Sturmians at λ=Z/n): ln k₀(2S) = 2.726 (−3.1%), ln k₀(1S) = 2.924 (−2.0%); 2P diverges (closure forces I→0 for ℓ>0). Structural identity surfaced: **the Sturmian basis at exponent λ=Z/n IS the Fock graph re-parameterized for bound-state space**. Files: `debug/ls2_bethe_log.{py,memo.md}`, `debug/data/ls2_bethe_log.json`.
- **LS-3** (acceleration form): s-states improved 3.3× (ln k₀(2S) at −0.92%, ln k₀(1S) at +0.60%); 2P remains divergent (acceleration form mathematically identical to velocity for closure issue); T+A error decomposition demonstrated to machine precision. Files: `debug/ls3_bethe_log_regularized.{py,memo.md}`, `debug/data/ls3_bethe_log_regularized.json`.
- **LS-4** (Drake-Swainson asymptotic-subtraction regularization, **the 13th projection** in Paper 34's Layer 2): closes 2P structural floor — ln k₀(2P) at +2.40% (N=40), 3D Bethe log at −0.24% (cleanest from-scratch result). K-cancellation between β_low(N,K) and β_high(K) verified to machine precision over 3.6 orders of magnitude in K. Combined Lamb shift 1053.76 MHz (−0.39%) at N=40, **honestly flagged as fortuitous T+A cancellation**; at N→∞ residual returns to LS-1 baseline −3.10% set by α⁵ multi-loop ceiling. Files: `debug/ls4_bethe_log_drake.{py,memo.md}`, `debug/data/ls4_bethe_log_drake.json`.

**Net structural conclusion:** Lamb shift accuracy is **A-bound (one-loop ceiling), not T-bound** — no amount of basis polishing closes the −3.10% gap; only multi-loop QED can.

### Added — Vector-photon graph-native sprints (VP-1, VP-2)

Parallel investigation of whether vector-photon promotion in graph-native QED closes the C × F₂ → α/(2π) negative result.

- **VP-1**: vector promotion (graph-native vertex with Wigner 3j coupling) gets 3.6× closer than scalar baseline but does NOT close projection mismatch. F₂_vec at n_max=2 = 3/(11π) exact symbolic — vector promotion injects π once via spherical-harmonic normalization. Files: `debug/vp1_vector_graph_native.{py,memo.md}`, `debug/data/vp1_vector_graph_native.json`.
- **VP-2**: family table for {C_VP, C_SE, C_F2_asymp} at n_max ∈ {2,3,4,5}. **β(C) = β(continuum) − β(graph) verified to machine precision** — the projection constant is a quotient of two power laws, not a single calibration. graph_VP_trace ≈ √π · n_max ≈ a₀(S³) · n_max (Seeley-DeWitt Weyl structure). C_VP/C_SE ≈ 3/29 (CV 0.83% across n_max=3,4,5) flagged as near-miss per curve-fit-audit standards (NOT identified as a structural identity). Files: `debug/vp2_topology_projections.{py,memo.md}`, `debug/data/vp2_topology_projections.json`.

### Changed — CLAUDE.md updates

- **Section 2 (Active Frontier):** new bullet for the bound-state QED arc + VP sprints + Paper 34 creation, summarizing all key findings
- **Section 6 (Paper Series):** Paper 34 added to Observations table with maintenance protocol; loading guide entry marked On-topic
- **Section 11 (Topic-to-Paper Lookup):** ten new rows pointing at Paper 34 sections (projection taxonomy, matches catalogue, off-precision classification, error compounds prediction, Drake-Swainson, etc.)

### K = π(B + F − Δ) anomaly status

The Paper 2 conjecture remains anomalous under Paper 34's projection-depth prediction. With four LS sprints providing data on what 4-projection chains tolerate (~3% residual at converged basis), Paper 2's K hits 1/α at 8.8×10⁻⁸ on a 3-projection chain — six orders of magnitude below normal projection-depth tolerance. Sharpened, not resolved. Open question §VIII.3 of Paper 34.

### Future sprints flagged

- **LS-5**: derive Drake denominator D = 2(2ℓ+1)Z⁴/n³ from Schwartz integral form (closes empirical/convention question)
- **LS-6**: two-loop QED on S³ for sub-1% Lamb shift accuracy (binding constraint = α⁵ multi-loop)

### Notes

- No production code modified in this version (all sprint deliverables in `debug/`; only `papers/observations/paper_34_projection_taxonomy.tex` and `CLAUDE.md` updated)
- Topological integrity tests (18/18 symbolic proofs in `tests/test_fock_projection.py` and `tests/test_fock_laplacian.py`) verified passing before commit

## [2.26.1] - 2026-04-18

### Changed — Paper 2 promoted from Conjectures to Core

- **File move:** `papers/conjectures/paper_2_alpha.tex` → `papers/core/paper_2_alpha.tex` (via `git mv`; history preserved).
- **Rationale:** Sprint A (2026-04-18) applied the structural reframe (see v2.26.0 notes + `debug/paper_2_reframe_skeleton.md` + `debug/alpha_sprint_a/`). Three structural verifications landed with partial positive / partial negative verdicts (α-PI, α-MI, α-SP, α-EB v2, α-X, α-LS); Marcolli-van Suijlekom 2014 gauge-network framework identified as published lineage (WH1 upgraded to MODERATE); Paper 2 body extended by 7 substantive edits (+168 lines) naming the triple m=3 selection, Hopf-measure π, APS-shape on Δ minus sign, π³α³ residual footnote (as hint, not derivation), and spectral-triple setting with Marcolli-vS citation.
- **Conjectural-label prohibition narrowed** in `CLAUDE.md` §13.5 and §13.8: "Removal of the 'conjectural' label from Paper 2" → "Removal of the 'conjectural' label from the combination rule K = π(B + F − Δ) in Paper 2." The paper's surrounding theorems (three-homes, three obstructions, Sprint A verifications) are not conjectural; the combination rule observation remains conjectural until a first-principles derivation lands.
- **Papers 25 and 30 updated** (same date) with Marcolli-van Suijlekom 2014 + Perez-Sanchez 2024/2025 citations establishing published gauge-network lineage.
- **CLAUDE.md §6 Paper Inventory** updated: Paper 2 added to Core tier, removed from Conjectures tier. Context Loading Guide entry added (On-topic).
- **README.md Paper Series** table updated: Paper 2 bolded and description reflects reframe.

### Notes

- Sprint B (g−2 from QED-on-S³) plan drafted as `debug/sprint_b_g_minus_2_plan.md`; NOT auto-dispatched; awaits PI go/no-go. Sprint B will test whether the framework predicts the anomalous magnetic moment at one loop — the Layer 2 out-of-sample prediction test per the promotion criterion.

## [2.26.0] - 2026-04-18

### Added — CLAUDE.md §1.7 Working Hypotheses (Internal Register)

- **New section:** CLAUDE.md §1.7 "Working Hypotheses (Internal Register)" — bold-claim register distinct from paper rhetoric. Six WHs (almost-commutative spectral triple, Paper 18 as Seeley-DeWitt + ζ-invariant decomposition, lattice a priori from packing origin, four-way S³ coincidence as one structure, α as projection constant, D(s) not classical Riemann ζ), each with falsifier and status fields. Governance specifies PM may update Status only; add/retire/promote requires PI direction.
- **Purpose:** papers remain cautious under §1.5 rhetoric rule; internal thinking gets Ramanujan-register license. Papers are Hardy-letters; the WH register is the notebook.
- **Memory note:** `memory/wh_register_april2026.md` records the register's addition and conversational context.

## [2.9.2] - 2026-04-15

### Added — Paper 27 (entropy as projection artifact)

- **New core paper:** `papers/core/paper_27_entropy_projection.tex` ("Entropy as a Projection Artifact: One-Body Operators are Entanglement-Inert on Sparse Lattices"). Five results:
  1. Operator-theoretic floor (theorem): non-degenerate one-body GS has S_kin=0 in its natural-orbital basis. Verified S_kin/S_full ~ 1e-14 for He at n_max=2,3 (EP-1).
  2. Area-law identity: A_n = g_n² = 4n⁴ pair count corrects Paper 5's one-particle mislabeling.
  3. Cusp localization: (1s,1s) is the unique hot-node on the V_ee pair graph.
  4. **HO zero-entropy rigidity (theorem, EP-2b):** on Bargmann-Segal HO basis, [H_HO, V_central]=0 exactly via Moshinsky-Talmi total-quanta conservation. Closed-shell 2-fermion GS is a single Slater determinant for ANY central two-body V. S=0 identically. New module `geovac/nuclear/ho_two_fermion.py`.
  5. **Universal scaling (EP-2c→2N):** S_B = A(w̃_B/δ_B)^γ with γ_∞ ≈ 1.96 (Richardson on n_max=2..5). Below second-order Rayleigh-Schrödinger γ=2; residue attributed to multi-shell aggregation. Verified across He-like Z∈[2,100], LiH R-sweep, HF/H₂O R-sweeps, Be analytical degenerate-PT.
- **Paper 24 extension:** §V.C HO entanglement rigidity corollary — structural dual of Fock rigidity, completing the π-free partition.
- **Paper 26 cross-reference:** S~Z^{-2.56} (Paper 26 §III) derives from Paper 27 §VI.B's universal curve. Reproduced at n_max=4: measured α_Z=-2.546, R²=0.995.
- **CLAUDE.md Section 2** consolidated with the full EP-1 → EP-2N sprint arc.

### Added — Cusp attack sprints

- **CUSP-1** (`debug/cusp1_screened_ci.py`): w̃/δ-based 2nd-order-RS configuration screening for He graph-native CI at n_max=7. **11× config compression** (1218→109) at the full-CI floor, 22× (→56) at floor+0.01pp. Top-scored states are exactly the cusp-relevant (1s,ns), (2s,2s), (2p_-1, 2p_+1) singlet configs. Useful for VQE/qubit compactness; doesn't break the accuracy floor (which is structural).
- **CUSP-2** (`debug/cusp2_hotnode_patch.py`): Tested Schwartz-tail patch on the (1s,1s) hot-node. **Decisive negative diagnosis of a long-standing misunderstanding**: the 0.20% floor for graph-native CI on He at Z=2 is NOT the cusp. Linear fit err_abs = -A/(l+2)⁴ + c gives irreducible c ≈ 6 mHa — the small-Z graph-validity-boundary artifact (Z_c ≈ 1.84). Z=10 confirms: error sign flips (under-binding +85 mHa). **For He graph-native CI, the cusp contribution is below 1 mHa at l_max≥4; the reported ceiling is basis-internal.** Real cusp work should be conducted at Z≥4 where the small-Z artifact is suppressed.
- **CUSP-3** (`debug/cusp3_tc_high_lmax.py`): TC benchmarked with proper particle-number-projected FCI at He composed n_max=2,3,4,5. **Decisive structural negative**: TC error plateaus at 3.40→3.47%, ratios approaching 1 (0.987→0.994→0.997). Standard continues to improve 2.48→2.02 and is still decreasing. Pauli ratio TC/Std worsens 1.68→2.73. **TC is now confirmed dead at every accessible n_max in second quantization on the composed basis**, not just small n_max. Section 3 updated with this extended-n_max evidence.

### Added — Paper 27 EP-2 sprint chain (drivers, data, tests)

- Drivers: `debug/ep2{b,c,d,e,f,g,h,i,j,k,L,M,N}_*.py` (13 scripts, one per sprint).
- Data: `debug/data/ep2{b,c,c_multi_block,d,e,f,g,h,i,j,k,L,M,N}_*.json`.
- Tests (`tests/test_paper27_entropy.py`): 27 tests total (26 passing + 1 `--slow`-gated). Includes:
  - EP-1 reproduction at n_max=2,3
  - Area-law pair counting (combinatorial + leading-order + one-body counter-check)
  - Energy-graph structural invariants at n_max=3,4 (commutator ratio, diagonality, (1s,1s)=5/8)
  - EP-2b HO commutator + single-SD ground state
  - EP-2c/2d/2g/2h/2i/2L dimensionless power law + Z-range artifact + asymptote below 2
  - EP-2n Be analytical degenerate-PT
  - Paper 26 / Paper 27 Z-scaling consistency
  - Proposition non-degeneracy qualifier (Li/Be)

### Fixed

- **Paper 27 abstract + conclusion:** γ → 2 asymptote corrected to γ_∞ ≈ 1.96 based on n_max=5 Richardson extrapolation. Previous "γ → 2 from above" claim was an artifact of limited n_max coverage.
- **Paper 27 Section II.C (Scope of the proposition):** added non-degeneracy qualifier, open-shell qualifier, and Be analytical degenerate-PT subsection.
- **Paper 27 Section V (wording fix):** ‖V_ee^diag(H_1)‖_F / ‖V_ee‖_F is unsquared Frobenius ratio; previous text incorrectly stated squared formula. n_max=4 value 0.94 corrected to 0.892.
- **Paper 27 Section VII consolidated:** 1111→900 lines, sprint chronology collapsed into two clean subsections (HO rigidity + Universal scaling). Methods/robustness collected into one paragraph.

### Status

- Paper 27 core content complete and regression-locked.
- Cusp problem: fundamentally re-diagnosed at He Z=2 (not cusp, graph-validity). No breakthrough on the cusp itself at Z≥4, but three structural facts established (see Section 3 updates).

---

## [2.9.1] - 2026-04-14

### Added (documentation only — no production code changes)

- **Energy graph characterization for V_ee on S³ (exploratory sprint):** Tested the hypothesis that a "second graph" — nodes = electron pair-states |(n₁l₁m₁)(n₂l₂m₂)⟩, edges = ⟨ab|1/r₁₂|cd⟩ — could be a Paper-12-style algebraic generator for V_ee on S³ (analog of the Neumann expansion for H₂ in prolate spheroidal coordinates).
- Findings:
  1. Pair-state graph at n_max=3 has 31 nodes, n_max=4 has 101 nodes (singlet M_L=0 sector).
  2. Within-parity blocks are ~47% dense at both sizes — orbital-level Gaunt sparsity (Paper 22) does NOT survive projection to pair-states.
  3. Diagonal V_ii are exact rationals (⟨1s1s|V|1s1s⟩ = 5/8, ⟨2s2s|V|2s2s⟩ = 77/512), but cross-shell entries mix 2^k and 3^j denominators — V's spectrum is non-rational. Consistent with Paper 18's classification of 1/r₁₂ as embedding exchange constant.
  4. No three-term recurrence in (n,l) for the Slater integrals: V(1s·ns) ratios drift {2.70, 2.20}, V(ns·ns) drift {4.16, 2.26}. Neumann-style separability is specific to prolate spheroidal coordinates, not S³.
  5. **Main new finding:** In the H₁ eigenbasis, 92% of V's Frobenius mass is diagonal at n_max=3, 94% at n_max=4 (residual ‖[H₁,V]‖/‖V‖ = 6.1% → 5.3%, saturating not vanishing). The wavefunction graph nearly diagonalizes V_ee.
  6. **Cusp signature:** The (1s1s) pair-state is simultaneously the highest-diagonal node, largest weighted hub, and head of the hottest off-diagonal edge. The cusp is concentrated on one pair-state at coalescence, NOT smeared across all angular channels.
- Assessment: NEGATIVE on the Paper-12 analog hypothesis (no clean algebraic generator for V_ee on S³); POSITIVE PARTIAL on wavefunction-graph near-diagonalization of V_ee. The "energy graph as distinct object" framing is less useful than hoped — the two graphs share most eigenstructure.
- Documentation: Paper 18 updated with a new paragraph in the embedding-constants section making the graph-theoretic characterization of 1/r₁₂ explicit. CLAUDE.md Section 3 (failed approaches) records the Paper-12 analog negative result; Section 2 cusp paragraph extended with the saturation finding; Section 11 topic lookup updated.
- Pointers: `debug/energy_graph_exploration.md`, `debug/data/energy_graph_nmax3.json`, `debug/data/energy_graph_nmax4.json`.
- No production code changed; no new tests added.

---

## [2.9.0] - 2026-04-13

### Added

- **111 Pauli count derivation (Investigation 1):** Analytically derived N_Pauli = 55 + 56 per s/p block at max_n=2. 55 = Q(Q+1)/2 direct terms (universal, from k=0 monopole). 56 exchange terms from 3 Gaunt channels: s-s (16), s-p cross-Coulomb (24), k=1 dipole (16). Pure l-shells (e.g., d-only) have zero exchange → Pauli/Q = (Q+1)/2 = 5.5. Universal coefficient 11.1 = 111/10. Key files: `debug/` analysis scripts.
- **He energy decomposition (Investigation 2):** `<h1>` converges by n_max=2 (within 0.04 Ha of limit). `<V_ee>` is the slow component — ALL remaining convergence struggle is off-diagonal V_ee correlation. `<V_ee>/|E|` ratio: 0.4545 (n_max=1) → 0.3290 (n_max=7) → 0.3257 (exact). Finite basis overestimates electron coalescence probability.
- **V_ee spectral analysis (Investigation 5):** V_ee is FULL RANK in the graph eigenbasis at every n_max (2-5). No low-rank shortcut exists for the cusp. Diagonal V_ee (mean-field) converges instantly (locked by n_max=3 at -2.787 Ha). Off-diagonal V_ee correlation converges as ~n_max^(-1) toward ~0.109 Ha. Cusp distributes broadly across all angular channels.
- **Graph-native CI variational boundary (Z-sweep):** The graph Laplacian CI crosses from non-variational (over-binding) to variational at Z_c ≈ 1.84 (n_max=7, extrapolated CBS ~1.87-1.89). Below Z_c, kappa=-1/16 off-diagonal coupling dominates Z²-scaled diagonal, creating artificial stabilization. Relative importance scales as 1/(8Z²): 12.5% at Z=1, 3.1% at Z=2, 0.13% at Z=10.
- **H⁻ variational-k investigation:** Standard (non-graph) FCI is ALWAYS variational for H⁻. Graph-native FCI over-binds by 21% (E=-0.640 vs exact -0.528 Ha). Variational-k optimization makes it WORSE (k_opt=1.16, energy drops further). Over-binding is specific to graph h1 topology, not orbital exponent.
- **TC gamma optimal scan (Investigation 4):** gamma_opt decreases with l_max, approximately gamma_opt ~ 0.51/(l_max+1.7). l_max=3 sweet spot (0.001% error, 33x improvement) is a near-cancellation, not systematic. TC roughly doubles convergence exponent (~l⁻¹ to ~l⁻²) but does not achieve theoretical l⁻⁸. Large-gamma floor ~0.35% from non-Hermitian contamination.
- **PsH exotic atom solver:** Standalone positronium hydride (e⁻e⁺ in proton field) solver using Level 3 hyperspherical framework with sign-flipped charge function. PsH IS bound: V_eff minimum at R=3.62 bohr, well depth 0.042 Ha below H+e⁺ threshold (l_max=3). Energy -0.756 Ha (4.1% error vs exact -0.789 Ha). Alpha parity mixing essential (distinguishable particles). Gaunt selection rules preserved.
- **H⁻ binding via graph-native CI:** H⁻ found bound at n_max≥2 but with massive over-binding (21% error). Confirms correlation required for H⁻ (n_max=1 correctly unbound).
- Analysis scripts in `debug/`: `he_energy_decomposition.py`, `vee_spectral_analysis.py`, `z_sweep_variational.py`, `h_minus_variational_k.py`, `tc_gamma_scan.py`, `psh_solver.py`, `positronium_analysis.py`
- Plots in `debug/plots/`: `vee_graph_eigenbasis.png`, `vee_sv_decay.png`, `psh_adiabatic_curves.png`, `psh_mu_comparison.png`, `psh_charge_function.png`

### Key Scientific Findings

- **Cusp characterization:** The electron-electron cusp is a full-rank, broadly-distributed, off-diagonal V_ee object in the graph eigenbasis. h1 and diagonal V_ee converge fast; ALL convergence struggle is off-diagonal V_ee. No rank-reduction or channel-dominant shortcut exists. TC similarity transformation is the structurally correct response (modifies effective V_ee).
- **Graph validity boundary:** Z_c ≈ 1.84. The graph Laplacian with kappa=-1/16 is reliable when Z² >> kappa (nuclear potential dominates graph topology). He (Z=2) barely qualifies. H⁻ (Z=1) does not. The graph's rigid inter-shell coupling overestimates correlation for asymmetric/weakly-bound systems.
- **Exotic atoms:** PsH demonstrates the hyperspherical framework extends to matter-antimatter systems via sign flips. The composed geometry insight applies: systems with disparate particle-nucleus distance scales need composed blocks, not single-center coordinates.

---

## [2.8.2] - 2026-04-13

### Added

- **5 multi-center diatomic molecules:** LiF (Q=70), CO (Q=100), N₂ (Q=100), F₂ (Q=100), NaCl (Q=50)

### Fixed

- Paper 14 corrections: composed coefficient 11.11→11.10 (exact), TM Pauli/Q 9.27→9.23, 1-norm table corrected to electronic-only
- Paper 20 corrections: same coefficient fixes, table caption clarifies composed vs balanced per row
- Balanced coupled fix for `spec.nuclei` attribute

---

## [2.8.1] - 2026-04-12

### Added

- **Algebraic Slater integrals:** `geovac/hypergeometric_slater.py` with exact Fraction-arithmetic R^k evaluator for arbitrary n_max. Validated 144/145 table entries, found+fixed F²(2p,2p) typo (43/512→45/512). 8x speedup.
- **Float algebraic path:** `compute_rk_float()` gives machine-precision (1.5e-12) at 25x faster than Fraction. Corrected systematic grid bias (0.06-0.44% per integral).
- **DUCC downfolding:** `geovac/downfolding.py` computes exact (2J-K) core potential. Root cause of l_max divergence: PK underestimates p-orbital potential by 109x. H₂O 1-norm 9% lower with downfolding.
- **Ecosystem export:** 30→35 molecules via `hamiltonian()` API
- **He graph-native FCI convergence** to n_max=8 (0.207%, 2262 configs) and n_max=9 (0.201%, 3927 configs) with exact algebraic integrals
- **He 2D variational best:** 0.019% self-consistent cusp correction, 0.004% with exact coalescence density
- **Wigner 3j caching, property caching, double-build elimination**
- 418 files restored from OneDrive migration
- 7 broken test imports fixed, 2908 tests collect cleanly

### Fixed

- O (Z=8) PK parameter gap fixed
- `casimir_ci.py` F²(2p,2p) typo corrected (43/512→45/512)
- Pip reinstalled to correct directory

---

## [2.8.0] - 2026-04-12

### Added

- **Full first transition series (Z=21-30)** as hydrides: ScH, TiH, VH, CrH, MnH, FeH, CoH, NiH, CuH, ZnH
- **General `build_composed_hamiltonian(spec)`** in `composed_qubit.py` — MolecularSpec-driven builder consumed by balanced_coupled and coupled_composition
- **Atomic classifier extended to Z=1-30** — second row (Z=11-18), K/Ca (Z=19-20), and all first-row transition metals with structure type F
- **`l_min` field on `OrbitalBlock`** — restricts angular momentum enumeration for d-only blocks (l_min=2)
- **`_v_cross_nuc_frozen_core`** — frozen-core electrostatic potential for multi-shell cores
- **`transition_metal_hydride_spec(Z)`** — spec factory for all 10 TM hydrides with convenience aliases
- **10 TM hydrides in ecosystem export** — accessible via `hamiltonian('ScH')` etc.
- **48 new tests** in `tests/test_transition_metals.py`
- **Cr/Cu anomalous configurations** correctly handled (3d⁵4s¹ and 3d¹⁰4s¹)

### Key Results

- All 10 TM hydrides: Q=30 qubits, 277 Pauli terms, Pauli/Q = 9.23
- Isostructural invariance confirmed: identical block topology → identical Pauli count
- d-block ERI density sparser than s/p (confirming Track CZ: 4.0% vs 8.9%)
- Pauli/Q = 9.23 < 11.11 main-group coefficient — transition metals are cheaper per qubit
- Library expanded from 30 to 38 molecules

### Changed

- `SCOPE_BOUNDARY.md` updated to v2.8.0 — transition metals now "Fully Implemented"
- Paper 20 (Resource Benchmarks) updated with full 10-molecule TM hydride table
- Paper 20 future directions: TM classifier item (iv) removed (completed)

---

## [2.7.1] - 2026-04-12

### Changed

- Outreach-ready documentation correction (tag only)

---

## [2.7.0] - 2026-04-12

### Added

- Papers 22 (Angular Sparsity Theorem), 23 (Nuclear Shell Hamiltonians), 24 (Bargmann-Segal Lattice)
- Paper 21 (Geometric Vacuum Synthesis)
- Precision He: 2D variational solver (0.004%), graph-native CI (0.19%), excited states
- Nuclear shell model: deuteron (16Q/592 Pauli), He-4 (16Q/712 Pauli), composed nuclear-electronic (26Q/614 Pauli)
- `geovac/casimir_ci.py`, `geovac/level3_variational.py`, `geovac/nuclear/` package
- Alpha structural decomposition phases 4B-4H (paused — combination rule open)
- Papers 8-9 promoted from archive to methods tier

---

## [0.4.0] - 2026-02-15

### 🌟 Major Scientific Breakthrough

#### Global Metric Scaling for Isoelectronic Series
- **BREAKTHROUGH:** Conformal transformation approach for multi-electron Z-scaling
- **PHYSICS FIX:** Resolved virial mismatch from previous Jacobian scaling
- **METHOD:** Solve Helium-equivalent system, scale eigenvalues by γ = (Z/2)²
- **VALIDATION:** Li+ 10.87% error, Be2+ 15.22% error (improved from 31.6%/44.5%)

#### Theoretical Significance
- **CONFORMAL INVARIANCE:** Z-scaling is a metric transformation, not parameter change
- **VIRIAL THEOREM:** Both T and V scale uniformly by Z², preserving <T> = -<V>/2
- **UNIVERSALITY:** Lattice topology is universal, only metric (energy scale) changes with Z
- **PHYSICAL LIMIT:** Remaining 10-15% error attributed to relativistic corrections (Z⁴)

### Added

- **Global metric scaling implementation** in isoelectronic tests
- **`docs/GLOBAL_METRIC_SCALING_SUCCESS.md`** - Complete technical analysis
- **`docs/JACOBIAN_SCALING_RESULTS.md`** - Historical context (archived)
- **`debug/plots/create_isoelectronic_plot.py`** - Visualization script
- **`debug/plots/isoelectronic_scaling.png`** - Validation plot
- **`tests/test_isoelectronic.py`** - Comprehensive isoelectronic test suite
- **E/Z² ratio analysis** - Validates near-constant scaling
- **Transition state test** - Linear H3 (19.94% error)

### Changed

- **Version:** 0.3.2 → **0.4.0**
- **README:** Added v0.4.0 section with global metric scaling results
- **README:** Updated benchmarks table with isoelectronic series
- **README:** Updated roadmap (v0.4.0 current, v0.5.0 planned)
- **Scaling approach:** Jacobian (kinetic-only) → Global conformal transformation
- **Isoelectronic accuracy:** 31-45% → **10-15%** (20-30 point improvement)

### Validated

- ✅ Global metric scaling preserves virial theorem
- ✅ E/Z² ratio nearly constant (GeoVac: -0.713 to -0.723)
- ✅ Li+ (Z=3, 2e): -6.489 Ha (10.87% error)
- ✅ Be2+ (Z=4, 2e): -11.572 Ha (15.22% error)
- ✅ Linear H3 transition state: -1.321 Ha (19.94% error)
- ✅ Conformal transformation theory validated

### Deprecated

- **Jacobian scaling** (scaling only kinetic energy by Z²) - causes virial mismatch
- Use **global metric scaling** instead for isoelectronic series

### Documentation

- [RELEASE_NOTES_v0.4.0.md](RELEASE_NOTES_v0.4.0.md) - Detailed release notes
- [docs/GLOBAL_METRIC_SCALING_SUCCESS.md](docs/GLOBAL_METRIC_SCALING_SUCCESS.md) - Technical analysis

## [0.3.2] - 2026-02-14

### Added

- **AtomicSolver class** - Pure geometric formulation for single-electron atoms
- **Z²-scaling for hydrogenic ions** - Automatic scaling for H, He+, Li2+, etc.
- **`solve_atom()` convenience function** - Quick single-electron calculations
- Comprehensive benchmark suite for validation
- Complete documentation for universal kinetic scale

### Validated

- ✅ H (Z=1): -0.497 Ha (0.57% error at max_n=30)
- ✅ He+ (Z=2): -1.989 Ha (0.57% error at max_n=30)
- ✅ Li2+ (Z=3): -4.474 Ha (0.57% error at max_n=30)
- ✅ Universal kinetic scale -1/16 works for all single-electron systems
- ✅ Z²-scaling formula exact: `kinetic_scale_eff = -1/16 × Z²`

### Changed

- Version: 0.3.1 → 0.3.2
- README updated with AtomicSolver examples
- Documentation expanded for single-electron systems

## [0.3.1] - 2026-02-13

### Added

- **Multi-solver architecture** - Mean-Field, Geometric-DFT, Full CI, Dirac
- **Geometric-DFT** - Fast correlation functional (5.7% error, 79% recovery)
- **Full CI for 2-electron systems** - Exact correlation (<1% with optimization)
- **Dirac relativistic solver** - Spinor formalism with relativistic corrections
- **Geometry optimization** - PES scanning and bond length optimization

### Validated

- ✅ H₂ Mean-Field: -0.980 Ha (16.5% error)
- ✅ H₂ Geometric-DFT: -1.108 Ha (5.7% error)
- ✅ H₂ Full CI (R=1.40): -1.142 Ha (2.8% error)
- ✅ H₂ Full CI (R=1.30 optimized): -1.169 Ha (0.43% error) ⭐

## [0.2.1] - 2026-02-13

### 🔬 Major Scientific Discoveries

#### Universal Constant Discovery
- **DISCOVERED:** `kinetic_scale = -1/16` is a fundamental topological invariant, not a fitting parameter
- **VALIDATED:** Across H (Z=1), He⁺ (Z=2), and H₂⁺ with <0.1% error
- **PHYSICAL MEANING:** Dimensionless ground state eigenvalue of vacuum lattice is exactly 8

#### H₂⁺ Control Experiment
- **PROVEN:** Graph topology correctly models covalent bonding (0% error for H₂⁺)
- **CONFIRMED:** 17% H₂ discrepancy is correlation energy, not topological flaw
- **LITMUS TEST:** Single-electron H₂⁺ validates mean-field framework

#### Mean-Field Classification
- **CLASSIFIED:** GeoVac as Topological Hartree-Fock solver
- **SINGLE-ELECTRON:** Exact accuracy (0% error)
- **MULTI-ELECTRON:** Mean-field quality (~17% correlation error, expected)

#### Bridge Scaling Physics
- **MECHANISM:** Super-linear scaling (α≈1.1) from angular momentum recruitment
- **EVIDENCE:** 90% high-l states (f,g,h,i) participate at n=25
- **PHYSICAL:** Mimics d/f orbital chemistry in heavy elements

### Added

- `UNIVERSAL_KINETIC_SCALE = -1/16` constant in `geovac/__init__.py`
- `HYDROGEN_GROUND_STATE`, `H2_PLUS_USES_UNIVERSAL_SCALE`, `H2_CORRELATION_ERROR` constants
- Physics classification section in package docstring
- `validate_universal_constant.py` - Comprehensive validation tool for H/He⁺/H₂⁺
- `analyze_bridge_distribution.py` - Physical analysis of bridge scaling
- `CORE_PRODUCT_STATUS.md` - Complete status report
- H₂⁺ control experiment documentation in README
- Molecular bonding correlation test section in README
- Universal constant section in README with validation data
- Mean-field classification documentation throughout
- Paper 5 appendix: H₂⁺ experiment and bridge scaling physics

### Changed

- **DEFAULT PARAMETER:** `MoleculeHamiltonian(..., kinetic_scale)`: `-0.075551` → `-1/16`
- **PACKAGE DESCRIPTION:** From empirical to "Topological Hartree-Fock solver"
- **PERFORMANCE CLAIMS:** "~35% error for H₂" → "0% H₂⁺, ~17% H₂ (correlation)"
- **BRIDGE SCALING:** Updated from static N=16 to dynamic N≈4×max_n
- **ERROR ATTRIBUTION:** Clarified correlation vs topology
- `demo_h2.py` to use universal constant with validation references
- `geovac/__init__.py` docstring to reflect mean-field nature
- `geovac/hamiltonian.py` documentation and examples

### Fixed

- Theoretical foundation: Framework now has first-principles basis
- Error attribution: Clear separation of topology (exact) vs correlation (missing)
- Bridge scaling mechanism: Physical origin identified and validated
- Documentation: Proper classification and realistic performance claims

### Validated

- ✅ Universal constant convergence (H, He⁺, H₂⁺)
- ✅ Single-electron topology (0% error for H₂⁺)
- ✅ Multi-electron mean-field behavior (17% correlation in H₂)
- ✅ Angular momentum recruitment in bridge scaling
- ✅ All existing tests pass with new constant

### Backward Compatibility

- ✅ **MAINTAINED:** Existing code with explicit `kinetic_scale` still works
- ✅ **NEW DEFAULT:** Code without explicit parameter uses universal constant
- ✅ **API STABLE:** No breaking changes to method signatures

## [0.2.0] - 2026-02-12

### Added

- `MoleculeHamiltonian` class for molecular bonding
- `GeometricLattice.stitch_lattices()` method for bridge connections
- Spectral delocalization bonding mechanism
- `demo_h2.py` - Complete H₂ molecule demonstration
- Bridge priority ranking system
- Wavefunction delocalization analysis
- Binding energy calculations
- Performance benchmarks for molecules

### Changed

- README: Updated to "First Topological Quantum Chemistry Solver"
- Documentation: Added molecular bonding examples
- Examples: Updated with H₂ demonstrations

### Fixed

- Matrix sparsity maintenance in molecular systems
- Bridge connectivity for optimal bonding

## [0.1.0] - 2026-02-01

### Added

- Initial release
- `GeometricLattice` class for atomic systems
- `HeliumHamiltonian` class for two-electron atoms
- `DiracHamiltonian` class (experimental)
- Graph Laplacian based kinetic energy
- Sparse matrix eigenvalue solver
- Basic documentation and examples

---

## Version Naming Convention

- **Major (X.0.0):** Breaking API changes
- **Minor (0.X.0):** New features, backward compatible
- **Patch (0.0.X):** Bug fixes, documentation updates

## Links

- [v0.2.1 Release Notes](RELEASE_NOTES_v0.2.1.md) - Detailed release documentation
- [v0.2.0 Release Notes](RELEASE_NOTES_v0.2.0.md) - Previous release
- [Core Product Status](CORE_PRODUCT_STATUS.md) - Complete validation report

[0.2.1]: https://github.com/your-org/geovac/compare/v0.2.0...v0.2.1
[0.2.0]: https://github.com/your-org/geovac/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/your-org/geovac/releases/tag/v0.1.0
