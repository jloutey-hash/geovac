# GeoVac Track Status

Last updated: 2026-05-08

---

## Track A: l_max Divergence in Composed Geometry

**Goal:** Determine whether the l_max divergence in LiH composed geometry originates in the Level 4 solver or the composition layer (PK pseudopotential), and fix it.

**Status:** COMPLETE (v2.0.6) — Root cause identified: intrinsic to adiabatic approximation with asymmetric charges. Fix identified: variational 2D solver (Paper 15 Sec VI.D). Blocked on Level 4 iteration time.

**Key findings:**
- Bare Level 4 drifts at +0.262 bohr/l_max with no Z_eff/PK/composition — divergence is intrinsic
- PK-induced symmetry breaking causes composed LiH drift (+0.303)
- Dead ends: algebraic PK projector (5.6x worse), enhanced Z_eff (destroys channel symmetry), single-channel DBOC (97% cancelled)

**Next step:** Reduce Level 4 iteration times for 2D solver integration into composition pipeline. Blocked until Levels 2-3 algebraicization is complete (Tracks C/D).

---

## Track B: Hyperradial Algebraicization (Level 3)

**Goal:** Replace the finite-difference angular solver at Level 3 with an algebraic Gegenbauer spectral basis, eliminating grid-resolution artifacts and enabling proper l_max convergence.

**Status:** COMPLETE (v2.0.6) — Algebraic angular solver implemented, coupled-channel integration working. Exact Q-matrix derived and implemented (Track D, 2026-03-28).

**Key results:**
- Algebraic solver: 10-30 dim matrix replaces 200-800 FD matrix
- Single-channel adiabatic: 0.16% at l_max=0 (removes lucky FD error cancellation from 0.05%)
- Coupled-channel: convergent 0.37% → 0.27% at l_max 1-3
- Exact Q-matrix (q_mode='exact'): 14-19% improvement over closure at l_max≥1
  - l_max=3: 0.219% (exact Q) vs 0.272% (closure) vs 0.651% (single-channel)
- l_max=0 overcorrection (1.1%) is structural — diagonal DBOC dominance, not fixable by Q improvement

**Backlog item (Q-matrix):** RESOLVED. Exact dP/dR derived algebraically from Hellmann-Feynman quantities. No finite differences needed. Verified against FD to <1e-9.

**Remaining path to sub-0.1%:** Increase n_channels and n_basis; or compute exact dP/dR for Q-matrix closure improvement (done), then increase coupled channels.

---

## Track C: Level 2 Algebraicization (Prolate Spheroidal)

**Goal:** Apply the Track B playbook to Paper 11's prolate spheroidal solver — replace FD radial grid with spectral basis.

**Status:** AUDIT COMPLETE (2026-03-28) — See `debug/level2_algebraic_audit.md`

**Relevant papers:** Paper 11, Mitnik et al. (2021, arXiv:2006.06616), Kereselidze et al. (2016)

**Relevant code:** `prolate_spheroidal_lattice.py`, `molecular_sturmian.py`

**Key findings:**
- Angular solver: 100% algebraic (spectral Legendre basis, 50x50 matrix)
- Radial solver: 100% FD (5000-point grid) — this is the target
- All radial matrix elements (kinetic, nuclear a*xi, confinement c²xi², centrifugal m²/(xi²-1)) have known algebraic replacements via Laguerre/Sturmian recurrence relations
- Literature provides two routes: (a) Mitnik et al. GSF with mapped Laguerre basis, (b) Kereselidze Heun equation Sturmians
- Expected impact: 5000-dim FD → ~15-dim spectral, accuracy from 0.70% to sub-0.01%, ~20x speedup
- CRITICAL: These are NOT the failed shared-p0 Sturmians from Papers 8-9. They operate within already-separated prolate spheroidal equations.

**Next step:** Implement spectral radial solver following Mitnik et al. mapped Laguerre approach. Estimated 2-3 sub-agent sessions.

---

## Track D: Level 3 Exact Q-Matrix

**Goal:** Replace Q ≈ P·P closure with exact Q = P·P + dP/dR to reduce coupled-channel overcorrection.

**Status:** COMPLETE (2026-03-28) — See `debug/level3_exact_q_results.md`

**Key results:**
- Algebraic dP/dR derived from Hellmann-Feynman quantities (no finite differences)
- Verified against central FD (delta=1e-7) to <1e-9 accuracy
- q_mode='exact' added to coupled-channel solver (additive, existing modes unchanged)
- 5 new tests added, all 15 coupled-channel tests pass
- l_max=0: marginal improvement (1.12% → 1.10%) — overcorrection is diagonal DBOC, not off-diagonal Q
- l_max≥1: 14-19% error reduction (significant)
- l_max=3: 0.219% (best result)
- Monotonic convergence preserved at all l_max

**Escalation for plan mode:** Consider updating CLAUDE.md Section 12 (Algebraic Registry) to add Q-matrix as algebraic (was algebraic-pending). Also update the backlog item text.

---

## Track Backlog

- **Qubit encoding publication readiness:** Paper 14 is the most publishable standalone result. Needs external reproduction.
- **H₂O accuracy:** Blocked on Level 4 angular basis at 6:1 charge asymmetry. May benefit from Track C results.
- **κ = -1/16 derivation:** Paper 2's second selection principle (2d_max = n²-1 at n=3) provides a structural route. Not currently active.
- **Public benchmarking:** Create a standalone reproduction script that an outsider can run to verify Paper 14's Pauli scaling claims against Gaussian baselines. High priority for credibility.
- **Level 2 spectral radial solver:** Ready for implementation (Track C audit complete). High priority — unblocks Level 4 improvements.
- **Level 3 n_channels convergence:** Test whether increasing n_channels from 3 to 5-10 with exact Q pushes He accuracy below 0.1%.

---

## Session Summary 2026-04-15 — Dirac-on-S³ Tier 1 Sprint

### Tracks
- D1: completed — dirac_s3.py module, 51/51 tests, π-free both sectors
- D2: completed — B does not lift; spectral gap obstruction
- D3: completed — F does not lift; Apéry Q-linear independence; ζ(3) new transcendental
- D4: completed — Hopf decomposition clean, no K-unification
- D5: completed — Paper 2 §IV rewrite, CLAUDE.md edits, Paper 18 ζ(3) entry (all applied)

### Results
Three-tier coincidence formally documented:
- B = 42 home = scalar Laplace-Beltrami (finite Casimir truncation with zero-mode factor (m-1))
- F = π²/6 home = scalar Fock-degeneracy Dirichlet series D_{n²}(d_max)
- Δ = 1/40 home = Dirac spectrum single-level g_3^Dirac(S³)
New transcendental: ζ(3) enters via Dirac Dirichlet series at s=4.
Hopf decomposition of Δ = 20+20 chirality split, alternating {4,5,4,5,4,5,4,5,4} over q ∈ {-2..+2}.

### Files Modified
- `papers/conjectures/paper_2_alpha.tex` — §IV rewrite applied
- `papers/core/paper_18_exchange_constants.tex` — new ζ(3) taxonomy entry
- `CLAUDE.md` — 5 mechanical edits (version bump, Phase 4I summary, §3 row, §10 benchmarks, §11 mappings)

### Files Created
- `geovac/dirac_s3.py`
- `tests/test_dirac_s3.py`
- `docs/dirac_s3_design_memo.md`
- `docs/dirac_s3_tier1_sprint_plan.md`
- `docs/dirac_s3_leader_prompt.md`
- `docs/dirac_s3_explorer_prompt.md`
- `docs/dirac_s3_verdict.md`
- `docs/paper2_section4_rewrite.tex` (proposal archive)
- `docs/claude_md_proposed_updates.md` (proposal archive)
- `debug/dirac_d{2,3,4}_*.py`, `debug/data/dirac_d{2,3,4}_*.json`, `debug/dirac_d{2,3,4}_memo.md`

### Decisions
- Common-generator theorem for α: dead — B, F, Δ have three categorically different spectral homes with proven non-unification obstructions.
- Tier 1b proof sprint: not opened — per sprint plan's "all negative" decision branch.
- Tier 2 (spin-ful composed qubit encoding): confirmed as next research phase.
- ζ(3) added to Paper 18 taxonomy as a new tier (odd-zeta, first-order operator origin).

---

## Session Summary 2026-04-15 — Dirac-on-S³ Tier 2 Sprint

### Tracks
- T0: completed — d_spinor(l_max) table (Paper 22 extension)
- T1: completed — spinor matrix elements in (κ, m_j) basis (117 tests pass)
- T2: completed — Breit-Pauli spin-orbit, H_SO = −Z⁴α²(κ+1)/[4n³l(l+½)(l+1)] (22 tests)
- T3: completed — spin-ful composed pipeline for LiH/BeH/CaH (13 regression tests + 164 pre-existing all pass)
- T4: completed — Sunaga 2025 market test (150–250× resource ratio at matched Q=18); fine-structure sign+OoM for He/Li/Be
- T5: completed — π-free spinor certificate + Paper 18 new subtier R_sp := ℚ(α²)[γ]/(γ²+(Zα)²−1) (25 tests)
- T6: completed — Paper 14 §V / Paper 20 Table / Paper 22 spinor section / Paper 18 subtier / CLAUDE.md v2.11.0 proposals drafted (drop-in files, NOT applied)

### Results
Three relativistic composed molecular Hamiltonians built end-to-end algebraically.
- Pauli ratio rel/scalar: 1.00×/2.42×/5.89× at n_max=1/2/3.
- 1-norm rel vs scalar: flat to slightly lower (QPE-favorable).
- d_spinor(l_max) ratio: 1/4 at l_max=0 → 0.92 at l_max=5.
- Sunaga RaH-18q: 47,099 Pauli vs GeoVac LiH 805 (0.017×) / CaH 534 (0.011×).
- Fine structure He/Li/Be: sign + OoM correct; absolute accuracy 50–200%.

### Files Modified (this sprint)
None yet — all Tier 2 proposals are drop-in files awaiting PI approval.

### Files Created
- geovac/dirac_matrix_elements.py (T1)
- geovac/spin_orbit.py (T2)
- geovac/composed_qubit_relativistic.py (T3)
- geovac/spinor_certificate.py (T5)
- tests/test_dirac_matrix_elements.py (T1)
- tests/test_spin_orbit.py (T2)
- tests/test_spin_ful_composed.py (T3)
- tests/test_spinor_certificate.py (T5)
- benchmarks/relativistic_comparison.py (T4)
- benchmarks/fine_structure_check.py (T4)
- debug/tier2_t0_spinor_density.py + data + memo (T0)
- docs/* (all Tier 2 memos and T6 proposals)
- debug/data/tier2_market/* (T4 data)

### Decisions
- Three-tier program of Tier 1 closed; Tier 2 engineering upgrade complete; Tier 3 (heavy atoms, γ corrections, QED extension) deferred.
- Paper 18 taxonomy now has executable per-tier certifiers (D1 scalar, T5 spinor); pattern to be extended for all future tiers.
- CLAUDE.md v2.11.0 bump awaits PI approval of the §1/§2/§3/§10/§11/§12 edits.

---

## Session Summary 2026-04-15 — Dirac-on-S³ Tier 3 Sprint

### Tracks
- T7: completed — gamma = sqrt(kappa^2-(Z*alpha)^2) radial corrections, exact for n_r=0, NR limit verified (31 tests)
- T8: completed — Darwin + mass-velocity alpha^4 ladder; Dirac formula exact for all (n,l,j) through n=4; honest negative on He/Li/Be splitting accuracy (43 tests)
- T9: completed — D^2 spectral zeta theorem: zeta_{D^2}(s) = 2^{2s-1}[lambda(2s-2)-lambda(2s)] = polynomial in pi^2 at every integer s; 4th cell of Paper 18 degenerate with calibration; operator order is the transcendental discriminant

### Results
- zeta_{D^2}(1) = -pi^2/4; zeta_{D^2}(2) = pi^2-pi^4/12; zeta_{D^2}(3) = pi^4(10-pi^2)/30; zeta_{D^2}(4) = pi^6(168-17pi^2)/1260
- No zeta(odd) content at any s — theorem, not numerical
- Darwin+MV don't improve 2p doublet splittings (shared l=1 -> Darwin=0, MV cancels)
- gamma radial <1/r> = Z(gamma*n_r+kappa^2)/(gamma*N_D^3), algebraic over Q(Z,alpha,gamma_kappa)

### Files Modified
- papers/core/paper_18_exchange_constants.tex — empty-cell closure paragraph + table row + bibitems added
- papers/core/paper_24_bargmann_segal.tex — spinor-Lichnerowicz corollary subsection + bibitem added
- papers/core/paper_14_qubit_encoding.tex — §V Dirac formula verification paragraph added
- CLAUDE.md — v2.12.0 + §2/§3/§10/§11/§12 edits
- debug/track_logs/STATUS.md — this session summary

### Files Created
- docs/tier3_verdict.md
- docs/paper18_empty_cell_proposal.tex
- docs/paper24_d2_corollary_proposal.tex
- docs/paper14_tier3_update_proposal.tex
- docs/claude_md_tier3_updates.md

### Decisions
- Paper 18 2x2 grid: 4th cell (2nd-order x spinor-bundle) is degenerate with calibration (pi^{even}). Grid effectively 3-tier: operator order is the primary discriminant.
- Fine-structure accuracy for multi-electron atoms requires Direction 3 (SS/SOO) — deferred.
- Full gamma radial corrections for n_r>=1 states require Kramers-Pasternak recursion — deferred.

---

## Session Summary 2026-05-08 — Three-Round Precision Catalogue Extension (post-Sprint MH)

### Tracks
- Round 1 Track 1 (Mu 1S-2S + Ps HFS): completed POSITIVE
- Round 1 Track 2 (chemistry re-test, LiH/NaH balanced coupled): completed; first-row no-op, second-row W1c-residual confirmed
- Round 1 Track 3 (CP² packing scoping): completed CLEAN NEGATIVE; rank-1 specificity
- Round 2 PK cross-center (NaH binding test): completed HONEST NEGATIVE
- Round 2 Mu 2S-2P Lamb: completed POSITIVE +0.013%
- Round 2 D 1S HFS: completed POSITIVE +40 ppm
- Round 3 Mu 1S HFS: completed POSITIVE-WITH-NUANCE +199 ppm (cleanest LS-8a isolation)
- Round 3 He 2³P fine structure: completed POSITIVE-PARTIAL (2 of 3 sub-percent)

### Results

**7 systems sub-percent on framework-native parts; 3+ axes covered:**

| System | Channel | Framework | Reference | Residual |
|:---|:---|:---|:---|:---|
| Mu 1S-2S | transition | rest-mass rescale H | Mu-MASS Crivelli 2018 | −0.11 ppm |
| Mu 2S-2P | Lamb shift | full Uehling + Drake-SW | Karshenboim 2005 | +0.013% |
| Mu 1S HFS | hyperfine | BF × (1+a_e)(1+a_μ) | Liu 1999 | +199 ppm |
| Ps 1S HFS | hyperfine | BF + Layer-2 annihilation | Ishida 2014 | +0.49% |
| D 1S HFS | hyperfine | BF strict | Wineland-Ramsey 1972 | +40 ppm |
| He 2³P₀-P₁ | fine struct | Breit-Pauli SO+SS+SOO | Pachucki/NIST | −0.014% |
| He 2³P₀-P₂ | fine struct | Breit-Pauli SO+SS+SOO | Pachucki/NIST | −0.201% |

**Coverage axes:** mass-hierarchy (H, μH, Mu, Ps); nuclear-spin (I=1/2, I=1); multi-focal-kind (external cross-register, internal multi-electron); QED-channel (Lamb shift, HFS, fine structure, transition).

**Key structural products:**
1. Internal multi-focal architecture (He 2³P) is angular-content-only — same mechanism as cross-register Roothaan termination
2. Empirical LS-8a wall calibration: ~200 ppm at α²(Zα) for HFS observables in clean isolation (Mu HFS)
3. CP² scoping: Paper 0 axioms rank-1 specific; W3 second-axiom question structurally distinct from "packing on different geometry"
4. Chemistry arc: PK cross-center 14.6% incremental on top of W1c's 17.5×; binding not recovered; named target = FrozenCore Z_eff(r) Schrödinger basis on Na valence

**Paper 34 row diff:** +6 §V machine-precision rows, +3 §V.B off-precision rows.

### Files Modified

- `papers/observations/paper_34_projection_taxonomy.tex` — 9 row additions across §V and §V.B
- `papers/observations/paper_36_bound_state_qed.tex` — Mu Lamb subsection + Mu triple closure subsection (~160 lines new)
- `papers/applications/paper_23_nuclear_shell.tex` — D HFS subsection (~80 lines) + 3 bibliography entries
- `papers/core/paper_14_qubit_encoding.tex` — He 2³P internal-multi-focal cross-reference
- `papers/core/paper_17_composed_geometries.tex` — §6.10 NEW "Architectural note: two solvers, two drift signatures"
- `papers/methods/paper_19_coupled_composition.tex` — NEW subsection "W1c-residual orthogonality wall"
- `geovac/balanced_coupled.py` — pk_cross_center kwarg added (bit-exact backward compat)
- `CLAUDE.md` §1 (version v2.32.1 → v2.33.0), §2 (8 sprint paragraphs), §3 (1 failed-approach row)
- `CHANGELOG.md` — v2.33.0 entry
- `debug/track_logs/STATUS.md` — this session summary

### Files Created

- `geovac/phillips_kleinman_cross_center.py` (PK cross-center module, ~310 lines)
- `tests/test_phillips_kleinman_cross_center.py` (21 tests, all pass)
- `debug/precision_catalogue_muonium.{py,memo.md,json}`
- `debug/precision_catalogue_positronium.{py,memo.md,json}`
- `debug/precision_catalogue_muonium_lamb.{py,memo.md,json}`
- `debug/precision_catalogue_deuterium_hfs.{py,memo.md,json}`
- `debug/precision_catalogue_muonium_hfs.{py,memo.md,json}`
- `debug/precision_catalogue_he_2_3p.{py,memo.md,json}`
- `debug/cp2_packing_scoping_memo.md` + `debug/data/cp2_packing_scoping_results.json`
- `debug/chemistry_solver_retest_lih*.{py,memo.md}` + `debug/chemistry_solver_retest_nah.{py,json}` + synthesis memo
- `debug/pk_cross_center_nah.{py,json}` + `debug/pk_cross_center_synthesis_memo.md`
- 6 new memory files (precision_catalogue_3rounds, internal_multifocal_angular_only, mu_hfs_ls8a_calibration, cp2_packing_rank1_specific, chemistry_arc_paused_w1c_residual, feedback_diagnostic_before_engineering)

### Decisions

- Multi-focal architecture verified across the full empirical-coverage matrix; precision catalogue arc has consistent momentum and is the durable evidence-producing direction.
- Chemistry arc paused per PI: a "more focused diagnostic arc may be needed; still feels like there's a missing something." Diagnostic-only sprint to be designed before another implementation sprint. Pattern: 3+ accumulated negatives in one direction → diagnostic before engineering (see new feedback memory).
- CP² and higher-rank packing closed cleanly as not packing-axiom directions; W3 second-packing-axiom question reframed as structurally distinct from "packing on different geometry" — about how parameters get values, not how a manifold gets cellular structure.
- Mu HFS sets the empirical LS-8a wall scale (~200 ppm at α²(Zα) for HFS observables in clean isolation) — Paper 35 Refined Prediction 1 quantitative anchor.

---

## Session Summary 2026-05-17 — Sprint L2-F.1 + Pythagorean Extension

### Tracks
- L2-F.1: completed POSITIVE — Probe Theorem 7.1 residual structure on Lorentzian wedge; bit-exact closed form r² = κ_g²·S(n)/(4π²) + D(n) via Pythagorean orthogonality
- Track 1 (parallel): completed UNIVERSAL — Six-witness HS-orthogonality probe; 18/18 cells bit-exact zero
- Track 2 (parallel): completed; SU(2) Wilson on S³ GO-WITH-PREREQS, SU(3) Bargmann on S⁵ NO-GO with four named obstructions
- Synthesis: completed — 5 paper edits across Papers 43/42/32/24
- Paper 24 cleanup: completed — preamble `\newtheorem` + `amsthm` + 2 cross-ref fixes (clean three-pass)
- Docs sync: completed — CLAUDE.md v2.45.0, MEMORY.md, CHANGELOG.md updated; 2 new memory files

### Results

**Pythagorean orthogonality (load-bearing structural finding).** On the hemispheric wedge of the Lorentzian Krein space at finite n_max, ⟨H_local, D_W^L⟩_HS = 0 bit-exact at every panel cell (max 8.9×10⁻¹⁶ across 18-cell six-witness panel; PSLQ-verified at 100 dps, ceiling 10⁶, n_max ∈ {1..6}). Closed form:

  r²(n; κ_g) = κ_g²·S(n)/(4π²) + D(n)
  S(n) = n(n+1)(n+2)(2n²+4n−1)/15  [cumulative ⟨m_j²⟩ wedge trace]
  D(n) = n(n+1)(n+2)(2n+1)(2n+3)/20  [cumulative ⟨|λ|²⟩ wedge trace]

The 1/π² is the master Mellin engine M1 Hopf-base-measure signature. Mechanism: H_local lives in the diagonal subspace of B(K_W) (function of m_j only), D_W^L lives in the off-diagonal subspace (Δn=±1, intertwines ±m_j chirality partners) — mutually orthogonal subspaces of operator space. Six-witness universality is a κ_g-linearity corollary of the BW result.

**Scoping verdicts.** SU(2) Wilson on S³ reachable with matter-coupling prereq (1-2 weeks of separate `geovac/su2_wilson_gauge.py` extension). SU(3) Wilson on S⁵ Bargmann NO-GO with four named structural obstructions (no spinor sector on (N,0) Hardy tower; no second-order/first-order distinction with separate Dirac; no half-integer wedge — m_l integer ⇒ period π not 2π; Coulomb/HO category mismatch resurfaces at modular-Hamiltonian level). Establishes **fourth layer** of Paper 24 §V Coulomb/HO asymmetry.

**Bit-exactness rule of thumb (PI-adopted heuristic).** New working principle: bit-exact closure = green light (skeleton operation); residuals = caution light (Layer 2 work); neither = drift detector needed. Used productively across all four tracks today; all stayed green-light.

### Files Modified
- `papers/standalone/paper_43_lorentzian_extension.tex` — §10.2 new subsec:pythagorean_orthogonality + Corollary; §11 O4 extended (~137 lines)
- `papers/standalone/paper_42_modular_hamiltonian_four_witness.tex` — §8 new rem:pythagorean_underlies_collapse + §10 O3 extended (~96 lines)
- `papers/synthesis/paper_32_spectral_triple.tex` — §VIII new rem:pythagorean_m1_closure (~104 lines); §VIII.C G4b paragraph revised
- `papers/core/paper_24_bargmann_segal.tex` — §V new subsec:asymmetry_layer4 (~122 lines); preamble: `+amsthm`, `+\newtheorem{theorem}`, `+\newtheorem{corollary}`; cross-ref fixes (sec:ho-rigidity-theorem → sec:rigidity ×3; sec:coulomb-ho-asymmetry → sec:asymmetry)
- `CLAUDE.md` — §1 v2.44.0 → v2.45.0; §1.7 WH1 sharpened (PROVEN maintained); §2 full sprint bullet; §6 4 paper entries × 2 tables = 8 edits; §11 3 new lookup rows
- `CHANGELOG.md` — v2.45.0 entry
- `MEMORY.md` — 2 one-line entries pointing at new memory files

### Files Created
- `memory/pythagorean_orthogonality.md` — project-type memory
- `memory/feedback_bit_exactness_rule.md` — feedback-type memory
- `debug/h_local_residual_pslq_compute.py`, `_closed_form_verify.py`, `_pslq_memo.md`
- `debug/data/h_local_residual_pslq_data.json`, `_closed_form.json`, `_final.json`
- `debug/six_witness_hs_orthogonality_compute.py`, `_memo.md`
- `debug/data/six_witness_hs_orthogonality.json`
- `debug/pythagorean_extension_scoping_memo.md`

### Decisions
- **Bit-exactness rule of thumb adopted** as ongoing PM heuristic. Documented in `memory/feedback_bit_exactness_rule.md`.
- **WH1 PROVEN status maintained**, not re-opened. L2-F.1 sharpens Theorem 7.1 with exact closed form; does not modify five-lemma GH-convergence proof.
- **SU(2) Wilson matter-coupling sprint named as next major direction** (3-5 weeks, gated on PI authorization). SU(3) NO-GO recorded as dictionary content.
- **Paper 32 §VIII.C G4b paragraph reframed** from "fourth layer" to "sibling of the four layers" since Paper 24 §V now formally records four layers explicitly.
- **Pre-existing Paper 24 bugs cleaned up** (same `\newtheorem` pattern as Paper 34 fix 2026-05-08; +amsthm for `\begin{proof}` envs; 2 underscore/hyphen cross-ref typos). Now compiles clean three-pass.
- **Pre-existing revtex4-2 `Note1` natbib quirk left as-is.** Benign warning, no PDF impact.
- **Strategic options on the table** for next session: (Option 2) formal proof of diagonal/off-diagonal subspace decomposition; (Option 3) SU(2) Wilson matter-coupling sprint, 3-5 weeks; (MEMORY.md compaction) separate sprint when PI has appetite.

### Honest Scope Flags
- Pythagorean orthogonality empirically PSLQ-verified at 100 dps + 18-cell panel; structural mechanism sketched in `debug/h_local_residual_pslq_memo.md`; formal proof at general (n_max, N_t) is named follow-on.
- MEMORY.md at 30.9 KB > 24.4 KB limit; pre-existed today.
- Paper 24 §V four-layer claim is structural, not numerical — no S⁵-side computation run; statement rests on the four named obstructions from the scoping memo.
