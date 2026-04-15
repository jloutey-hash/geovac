# GeoVac Track Status

Last updated: 2026-03-28

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
