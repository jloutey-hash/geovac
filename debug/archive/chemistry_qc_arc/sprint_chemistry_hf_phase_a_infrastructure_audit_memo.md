# Sprint Chemistry HF — Phase A Infrastructure Audit

**Date:** 2026-05-25
**Sprint:** Pre-implementation audit for explicit-core 12-electron HF diagnostic on NaH (Reading A vs Reading B decision-gate)
**Mode:** Read-only audit; no file modifications.
**Predecessor:** `debug/sprint_chemistry_position_w1e_structural_memo.md` §5.1 (named the diagnostic).

---

## §0. Headline verdict

**Classification: (b) Significant infrastructure gaps; Phase B is 4.5 weeks critical path; sprint completes in 6–8 weeks total.**

The original Track 1 §5.1 estimate of "~2–4 weeks" was optimistic by 1.5–3×. The framework has strong foundations (two-electron prolate SCF, FrozenCore screening data, hypergeometric Slater integrals, multi-zeta orbital tabulations) but no N-electron HF solver, no core-orbital inversion mode, and no multi-center two-body integral assembly. Each of these is well-defined but non-trivial engineering work.

**Three sub-findings:**
1. **Multi-electron HF solver extension** (~1.5 weeks) — scaffold off `prolate_scf.py` two-electron grid SCF; need new Fock matrix logic, multi-electron ERI access, convergence hardening.
2. **Core orbital inversion** (~1 week) — straightforward wrapper around existing `neon_core.py` hydrogenic radials; low risk.
3. **Multi-center ERI assembly** (~1 week) — Slater integrals exist single-center via `hypergeometric_slater.py`; cross-center via prolate transformation valid for diatomics but needs validation.

---

## §1. Infrastructure inventory (per audit-question 1–7)

### 1.1 HF / SCF solvers
**EXISTS:** `geovac/prolate_scf.py` (~1005 lines) with `grid_scf()` (lines 577–769): full iterative SCF on 2D (ξ,η) finite-difference grid, Coulomb potential via azimuthal averaging (elliptic K, lines 301–393), eigenvalue solve via 2D Hamiltonian (lines 400–570), orbital damping (line 738). `geovac/prolate_heteronuclear_scf.py` (360 lines) with `heteronuclear_scf_energy()` (lines 36–128) for heteronuclear HF energy evaluation. Both: **two-electron closed-shell only.**

**MISSING:** Multi-electron HF solver for N > 2. No Fock matrix assembly for N-electron systems. No DIIS / level-shifting (only simple damping).

**Phase B work:** ~1.5 weeks. Risk: NaH HF convergence at stretched geometries historically requires DIIS or level-shifting — existing damping may be insufficient.

### 1.2 Frozen [Ne] core lifting
**EXISTS:** `geovac/neon_core.py` (FrozenCore class, lines 226–611) with Clementi-Raimondi Slater exponents (lines 48–58), hydrogenic `_hydrogenic_radial()` (lines 195–223), core density (lines 349–359), Z_eff(r) screening profile (lines 536–564), NIST/Clementi-Roetti core energies (lines 129–192).

**MISSING:** Inversion mode — no method exposing the 10 core orbital wavefunctions {1s, 2s, 2p×3} × 2 spins as explicit basis functions for downstream FCI. Atomic data exists; assembly wrapper does not.

**Phase B work:** ~1 week. Low risk.

### 1.3 Basis enumeration for 12 electrons at max_n=2
**EXISTS:** `geovac/molecular_spec.py` (OrbitalBlock, MolecularSpec dataclasses); `geovac/composed_qubit.py` (block-wise enumeration). Current NaH at max_n=2: 2-electron active in Q=20 (10 spatial × 2 spins).

**NEEDED:** Lift to 15 spatial orbitals (5 [Ne] core + 10 valence). Mechanical extension of existing block builder.

**Phase B work:** ~3 days. Low risk.

### 1.4 Two-electron integrals for the multi-center 12-electron problem
**EXISTS:** `geovac/hypergeometric_slater.py` (exact Fraction-arithmetic R^k integrals at arbitrary n_max, single-center hydrogenic). `geovac/prolate_scf.py` lines 202–269: azimuthal-averaged two-electron integrals via elliptic K kernel for prolate systems. `geovac/cross_block_h1.py` lines 1–74: one-body cross-center integrals, s-s pairs only.

**MISSING:** Multi-center two-body (ij|kl) integrals across the full 15-orbital basis (~50,625 unique integrals, sparse with ~0.1–1% non-negligible). Cross-center 2-body is flagged as named follow-on in `cross_block_h1.py` scope statement (line 22).

**Two approaches for Phase B:**
- **(A) Direct quadrature:** Extend cross-block to two-body via 3D Gauss-Legendre. Correct, ~2 weeks, slow.
- **(B) Analytical Slater + prolate transformation:** Use `hypergeometric_slater.py` for single-center; cross-center via prolate spheroidal coords + azimuthal averaging (valid for diatomics). ~1 week, fastest.

**Recommended:** Approach B + cache.

**Phase B work:** ~1 week. Medium risk: cross-center validation needs care (core 2p alignment with bond axis; Boys-integral analog for Slater orbitals not standard).

### 1.5 Multi-zeta basis support
**EXISTS:** `geovac/multi_zeta_orbitals.py` with Roothaan-HF orbital tabulations (Clementi-Roetti 1974, Bunge-Barrientos-Bunge 1993, Koga-Tatewaki 2000); Na valence physical-fit multi-zeta tabulated (lines 50–64); Xe-core multi-zeta loaded.

**STATUS:** Na valence multi-zeta ready to use. Core multi-zeta (Ar, Kr) not in scope yet — single-zeta CR67 core acceptable for first diagnostic.

**Phase B work:** 0 weeks additional.

### 1.6 Single-determinant HF precedent
**CLOSEST:** `prolate_heteronuclear_scf.py::heteronuclear_scf_energy()` (two-electron HF energy evaluation). `prolate_scf.py::eckart_scf_energy()` (Eckart variational SCF for H₂). `prolate_scf.py::relaxed_orbital_ci()` (CI with optimized orbitals).

**MISSING:** No N-electron HF anywhere in codebase. `casimir_ci.py`, `direct_ci.py`, `coupled_composition.py` are FCI, not HF.

### 1.7 Wall-time estimate for 12-electron HF at max_n=2
- Fock assembly per iteration: O(N_basis² × N_el × N_basis²) ≈ ~2,000 integral lookups; ~225 μs on modern CPU
- Per-iteration: Fock + 15×15 eigensolve + density update ≈ 1–2 ms
- SCF convergence: typically 20–50 iterations for NaH HF (with damping)
- Full SCF per geometry: ~50–200 ms
- 8-point PES: ~0.4–1.6 seconds wall time
- Plus ERI precomputation: ~0.1–1.0 seconds

**Actual compute is trivial.** Implementation/testing time dominates — Phase B 4.5 weeks, Phase C 1–1.5 weeks, Phase D 1 week.

---

## §2. Critical-path analysis: minimum Phase B

**Sequential dependency:**
1. **Core orbital inversion** (1 week, blocker for step 3) — extend `neon_core.py` with `get_core_orbitals_explicit()` returning {R_1s, R_2s, R_2p×3} on grid + orthonormalization.
2. **ERI assembly** (1 week, parallel to step 3) — Slater single-center via `hypergeometric_slater.py`; cross-center via prolate transformation; ERI dictionary caching.
3. **Multi-electron HF solver** (1.5 weeks, depends on 1+2) — scaffold on prolate grid infrastructure; Fock matrix assembly; damping + optional DIIS.
4. **Basis integration + FCI interface** (1 week, depends on step 3) — wire core+valence into MolecularSpec; validate orthonormality.
5. **Integration testing** (1 week, parallel to step 4) — unit tests; benchmark against Molpro/GAMESS HF reference on H₂, HeH⁺, NaH at equilibrium.

**Critical path: 4.5 weeks minimum.**

---

## §3. Revised sprint estimate

| Phase | Original | Revised | Why |
|:------|:--------:|:-------:|:----|
| A (audit) | — | done | This memo. |
| B (impl) | ~1 wk | **4.5 wk** | No N-electron HF; no core inversion; multi-center ERI assembly is real work, not glue. |
| C (hardening) | included | **1–1.5 wk** | NaH HF convergence at stretched R historically needs DIIS / level-shifting. |
| D (PES + analysis) | included | **1 wk** | 8-point sweep + reference comparison + verdict + paper edits. |
| Slack | — | **0.5 wk** | Integration testing, debugging. |
| **Total** | **2–4 wk** | **6–8 wk** | |

---

## §4. Risk catalogue

### 4.1 NaH HF convergence at stretched geometries — MEDIUM RISK
Standard issue (Szabo–Ostlund §3.5): HF near bonding-antibonding crossings exhibits oscillation. Existing damping (line 738 `prolate_scf.py`) is two-electron-only. Need DIIS or level-shifting. Mitigation: start at R=3.0 (equilibrium), march outward with adaptive damping; warm-start each geometry from previous converged density. **Cost:** +1 week Phase C.

### 4.2 Slater (STO) vs Gaussian (GTO) basis differences — LOW-MEDIUM RISK
Framework uses STOs; standard quantum chemistry uses GTOs. STO-HF is chemically sound but may converge slower or differ <1% from GTO-HF reference. **Cost:** Acceptable; document as known limitation.

### 4.3 Multi-center two-body integrals — prolate transformation validity — MEDIUM RISK
Prolate spheroidal coords are canonical for diatomics with bond on z-axis. Core 2p_x, 2p_y orbitals are off-axis — elliptic K azimuthal averaging may need explicit transformation. **Mitigation:** Validate cross-center integrals against analytical Boys-analog for Slaters; or fall back to 3D Gauss-Legendre quadrature for cross-center core integrals. **Cost:** +3–5 days validation in Phase B.

### 4.4 ERI dictionary size — LOW RISK
~50,000 unique integrals, ~0.5–1 GB memory for dense storage. Standard quantum chemistry caching technique.

### 4.5 Does framework HF give binding at all? — HIGH-IF-BROKEN, LOW-IF-NOT
Standard HF-on-NaH gives ~80% of experimental D_e (Clementi tables, modern Molpro/GAMESS). If framework HF doesn't reproduce this, indicates deeper framework bug. **Mitigation:** Unit-test framework HF against reference code on H₂ and HeH⁺ before NaH PES sweep. **Cost:** +2–3 days Phase B for validation.

---

## §5. Three options for the user

### Option A: Commit to the 6–8 week sprint
Full Phase B+C+D as described. Decision-gate verdict in 6–8 weeks. Multi-month CCSD-class follow-on (Reading A) or chemistry-domain architectural-scope paper (Reading B) becomes the next major arc after that.

**Cost:** 6–8 weeks of focused engineering. Real opportunity cost (math.OA polish, precision catalogue, Paper 49 follow-ons all paused).

### Option B: Rescope to a simpler diagnostic
Possible simpler diagnostics that might distinguish Reading A from B without building full 12-electron HF:

- **Literature HF reference comparison.** Look up published HF-on-NaH (Clementi tables, Molpro benchmark) and compare to framework's basis quality. If literature HF gives binding and framework basis is qualitatively comparable, infer Reading A. If framework basis is qualitatively different, more care needed. ~3 days. Weak verdict.

- **4-electron HF test on a simpler system.** Build minimal HF infrastructure for a 3–4 electron system (LiH at 4 electrons) where the framework's existing 2-electron HF scaffolding extends more directly. If framework 4-electron HF reproduces LiH binding, validates infrastructure path forward without committing to NaH. ~2 weeks. Validates path; doesn't directly settle W1e.

- **Variational HF energy bound from existing 2-electron FCI.** Use FCI energies to constrain HF energy via Pauli-pressure arguments. Sharper if a tight lower bound on HF energy can be proven. ~1 week diagnostic. Probably weak.

### Option C: Pause chemistry arc
Accept the W1e structural reading at current maturity (sprint-scale exhausted, multi-month tractable for components, three-class partition documented). Pivot to math.OA / precision catalogue / Paper 49 polish. Return to chemistry when the cost-benefit shifts.

**Cost:** None to the framework. The structural-skeleton-scope statement for chemistry remains "Reading A or B undetermined" — honestly recorded as open question.

---

## §6. Recommendation

**Honest assessment:** Option A (commit) is a real 6–8 week investment with uncertain strategic value. The W1e diagnostic settles a structural question about chemistry, but the chemistry arc has been "second-tier priority" relative to the math.OA arc for months. Running the full sprint is justified IF (i) the chemistry-domain structural-boundary verdict is load-bearing for the framework's positioning, (ii) Reading B (if confirmed) would be paper-worthy at synthesis level (it would — chemistry analog of H1 is genuinely new content), and (iii) opportunity cost is acceptable.

**Option B (rescope)** is the diagnostic-before-engineering-faithful choice. The literature-reference comparison (~3 days) is essentially free and might give a defensible Reading A inference without infrastructure investment. The 4-electron LiH test (~2 weeks) is a credible middle ground — builds half the infrastructure, validates the path, settles less of the question.

**Option C (pause)** is the rational default if the chemistry arc is not load-bearing for the framework's near-term direction. The W1e wall is documented; the three-class partition is recorded; the diagnostic is named and queued. Return-when-warranted is a legitimate choice.

**My recommendation: present all three options to the PI; this is a strategic-direction call, not a technical one.** The Phase A audit's value is the revised scope estimate (6–8 wk, not 2–4 wk); the cost-benefit shifts substantially at that scope.

---

**End of memo.**
