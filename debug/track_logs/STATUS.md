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
