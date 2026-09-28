# Build plan — the π / M=16 marriage: reach the −8.047 window with a π-containing reference

**Written 2026-09-23 (v5.15.24), PI-directed, for execution by a FRESH context using subagents.**
Self-contained. Before starting, read (do NOT re-derive what they settle):
`memory/lih_marriage_state_of_play.md` (READ FIRST), `debug/sprint_lih_marriage_memo.md`
(Phases 0/0b/1/2/3 — what is built and validated), `debug/lih_marriage_build_plan.md` (the σ
marriage this extends), and CHANGELOG v5.15.24. This is the **near follow-on** the marriage flagged.

---

## THE GOAL
Run the already-built RI-free marriage machinery on a **π-containing M=16 reference** instead of the
σ-only M=12 one, and obtain a deterministic, RI-free, VMC-cross-checked LiH energy in the **−8.047
window**. This turns the σ marriage (−8.0372, a proof of concept) into the **actual deterministic
cross-check of the VMC-over-FCI −8.047** — the marriage's original stated goal.

## WHAT THIS WILL AND WILL NOT DO (honest cap, set before starting)
- The correlation ansatz stays **two-body pairwise** (F = Σ f(r_ij)). A π reference sits at the
  orbital ceiling ~−8.032; a two-body r₁₂ treatment recovers ~14–15 mHa on this system (measured:
  the σ marriage recovered 14.3 mHa; the VMC-Padé recovered ~15). So **expect E_R12 ≈ −8.046..−8.047**.
- This does **NOT** reach −8.062. That needs the nucleus-coupled r₁₂t² **3-body** term — a SEPARATE
  follow-on (its own plan), not this one. Do not sell this as −8.062.
- The real test is **agreement with the VMC-Padé −8.047 ± 0.001** (both are 2-body treatments on the
  same M=16 reference, so they must agree). Hitting a preconceived number is not the point; the
  deterministic-vs-stochastic agreement is.

## WHY THIS IS FEASIBLE (the hard analytic part is DONE)
The σ-only case was easy because every density was φ-independent (single angular mode m=0). π orbitals
carry circulation around the bond axis (m=±1; δ = ±2), so pair densities acquire e^{imφ} content and
need the **general-m** Coulomb/Neumann/triangle kernels. **Those general-m kernels are already built
and validated:** `neumann_exact.ExactNeumann` is generic in (l, m) and was checked to 1e-13 for
m=0..4 (Phase 0b); `triangle.coul_mode_potential(B2d, m, exact=True)`, `build_kernel_m`, `Kmu_modes`
are the general-m triangle machinery; the bridge reduction was validated on a p-like (π-flavored)
density in Phase 1 (the l=4 channel). So the deep analytics are done. What is new is bookkeeping and
compute scale, not physics.

## CURRENT STATE (verified 2026-09-23)
- σ marriage: E_R12 = −8.0372 (linexp −8.037236 / exp −8.037088), G2-cross-checked, on the 74-det
  M=12 σ reference (E_trunc −8.022918). This is the **Phase-B regression anchor**.
- Entry points that already exist: `debug/lih_marriage_phase2.py` `all_pieces(gem, ref)` +
  `assemble(pieces, E0_grid)` (the generalized enumerators, currently **m=0 only** — this is what
  Phase B extends); the g_Vee productionization (`get_tri_matrix` per-middle-vertex triangle cache,
  `coul_mode_stack` batched mode-m Neumann); `KG.USE_EXACT_NEUMANN=True`. `debug/lih_marriage_phase0.py`
  builds the coefficient tensor T + spectator traces for a reference (currently σ/m=0).
- `debug/lih_vmc.py` `build_lih_wavefunction(Jb=2, Lb=1, npi=1, Jpi=1, Lpi=0, alpha=1.0, core2=(4.5,1.6))`
  — the DEFAULT call already includes π (npi=1); the σ marriage used the npi=0 override. The VMC-over-FCI
  headline used the **complex-π M=16** basis (gate-6 passed there), so the reference side handles π and
  complex orbitals already. The linear Jastrow (1+c(F−F̄))Ψ_CI is in `lih_vmc.py` (Phase 3, additive).
- The state-of-play flags the M=16+π case as **~136 pair densities, T ≤ 3.4e8 before truncation,
  m-resolved kernels — LARGE, not blocked; truncate to leading dets.**

## THE BUILD — phased, one subagent per phase, one gate, one capture
Each phase: dispatch one agent with the phase text as its TASK, the named files, the gate as its
DECISION GATE; it returns a ≤400-word structured report + a `debug/*.py` driver + `debug/data/*.log`.
The PM integrates, updates this STATUS block, captures in Paper 12/19 + CHANGELOG, THEN dispatches the
next phase. Never hold more than one phase in the main context.

**Phase A — π reference + m-resolved tensor (the reference upgrade).**
Build the M=16 π-containing CI (the reference the VMC used; orbital ceiling ~−8.032) via
`build_lih_wavefunction` (π defaults), truncate to leading dets in the natural-orbital basis
(as Phase 0 did for σ), and build the coefficient tensor T with an **m-label on every pair density**
plus the m-resolved pair densities on the grid. Extend `debug/lih_marriage_phase0.py`'s tensor builder
to carry complex orbitals and m-labels. **Gate G0-π:** (i) spectator traces of T reproduce the 1-/2-RDMs
of the truncated vector to ~1e-14 (the Phase-0 bookkeeping check, now complex/m-resolved); (ii) grid
⟨T⟩, ⟨V_ne⟩, ⟨V_ee⟩ vs the engine ≤ 1 mHa, with ⟨V_ee⟩ now summing the m-resolved Neumann channels
(the general-m operator is validated, so a miss here is a bookkeeping bug). **STOP if G0-π fails** →
report which m-channel/pair; do not proceed to the energy.

**Phase B — extend the 7 enumerators to m-channels (the load-bearing regression).**
Generalize `all_pieces`'s enumerators (currently m=0) to sum over angular modes with the **L_z
selection rule** (total circulation conserved around the axis — this keeps the m-coupling SPARSE, it
is a help not a cost). **Gate G1-π (regression, mandatory):** feeding a **σ-only reference** (npi=0)
through the m-extended enumerators must reproduce the banked σ marriage **−8.0372** and the per-piece
σ values to <0.1 mHa (i.e. the m-extension collapses to the known case when π coefficients are zero).
This is the "did the extension break the proven case" gate — do NOT skip it or start it at a context
tail. If G1-π fails, the m-bookkeeping is wrong; report which piece.

**Phase C — the M=16 energy + VMC cross-check (g_Vee is the wall-clock wall).**
Run `all_pieces('linexp', ref_M16)` + `assemble(...)` on the exact path. **g_Vee is the cost** (136
m-resolved pairs, several m-channels/triple; the Kmu caching + `coul_mode_stack` pay off more here but
there are more distinct (vertex, m) combos). Budget for it: cache aggressively, and if it exceeds wall
budget, **truncate the det list harder** (trade a little accuracy) rather than abandon — report the
truncation. Then **Gate G2-π:** VMC of (1+c(F−F̄))Ψ_CI on the M=16 reference (standard local energy,
NOT the gradient form), E_VMC == E_R12 within ±1 mHa; and gate-6 (c=0) reproduces the M=16 E_trunc.

## DECISION GATE
- **GO:** E_R12 ∈ ~[−8.050, −8.043], > −8.0705, G0-π/G1-π/G2-π pass, **AND E_R12 agrees with the
  VMC-Padé −8.047 ± 0.001 within ~2 mHa** → the deterministic RI-free cross-check of the −8.047 number
  is achieved; capture as the LiH marriage's headline (Paper 19, and the §2/§5 tables can then carry
  the deterministic −8.047-class number alongside the VMC).
- **BORDERLINE:** G2-π agrees but E_R12 is a few mHa short of −8.047 → likely the det truncation;
  widen it, re-run Phase C only.
- **STOP:** G0-π or G1-π fails (bookkeeping), or g_Vee is genuinely intractable even after truncation
  → report the blocker; the machinery correctness is not in question (G1-π guards it), only the scale.

## ANCHORS
−8.0372 (σ regression — Phase B must reproduce) · −8.032 (M=16 orbital ceiling, the reference floor) ·
**−8.047 ± 0.001 (the VMC target this cross-checks — the real test)** · −8.0705 (exact, variational
falsifier, must stay above).

## TRAPS / DO NOT RE-DERIVE
- The general-m Coulomb/Neumann/triangle machinery is BUILT and validated (Phase 0b/1) — use it, do
  not rebuild. `exact=True` everywhere; never the legacy kink-crossing sum.
- **Phase B's σ-zeroing regression to −8.0372 is load-bearing** — it is what proves the m-extension
  didn't break the proven case. Do not skip it.
- Complex orbitals: the π basis is complex (e^{±imφ}); the tensor and enumerators must carry complex
  coefficients (the VMC already ran "complex-π M=16" and passed gate-6, so the reference side is fine).
- Do NOT double-count: the π reference already captures some sideways correlation, so r₁₂ recovers
  what's left. The **VMC-Padé −8.047 is the ground truth** for "2-body on M=16"; agreement with THAT,
  not with a preconceived −8.047, is the test.
- g_Vee wall-clock is the real risk — plan the compute, truncate dets if needed, don't let it silently
  blow up. Pin `*_NUM_THREADS=1` before numpy import if spawning MC workers (BLAS oversubscription bit
  Phase 1).
- Do NOT promise −8.062 (that's the 3-body follow-on). Do NOT change the §2/§5 best-results until GO.

## FOLLOW-ON (after this GO): the 3-body path to −8.062
Add the nucleus-coupled r₁₂t²-type geminal to F (the leaf dressings acquire a radial weight; the bridge
structure is unchanged) and see whether the same reductions carry it. That is the elegant, RI-free
alternative to a brute-force 3-body Jastrow, and the platform this build leaves behind.

## STATUS
2026-09-23: plan written.
2026-09-24: **Phases A+B DONE (module `debug/lih_pi_m16.py`); Phase C running.**
- **A+B pass.** M=16 π reference: E_full −8.029045, and its m-resolved coefficient tensor built.
  **G0-π PASS** (tensor traces vs RDMs 3.3e-16; grid ⟨T⟩/⟨V_ne⟩/⟨V_ee⟩ vs engine all 0.0000 mHa —
  the m-resolved π Coulomb validated to machine precision). Enumerators extended to m-channels via an
  L_z circulation mask (Σ mA=0 at each integrated vertex). **G1-π PASS** — a σ reference through the
  m-enumerators reproduces the Phase-2 σ marriage per piece to 0.000000 mHa and E_R12 −8.037236 by
  construction (the m-extension collapses exactly to the proven case).
- **Phase C, trunc=1e-2 (nact=64): gated GO.** E_R12(linexp) = **−8.042181** (E_trunc −8.025184),
  variational, beats the M=16 ceiling −8.029/−8.032, 4.8 mHa from VMC-Padé −8.047. **G2-π PASS** at
  high VMC statistics: E_VMC −8.042500 ± 0.000172 vs analytic −8.042181, d = 0.319 mHa (the first
  short VMC's 1.1 mHa was its own ±0.55 mHa noise, resolved by 10× samples; gate-6 +0.5σ). Result +
  preserved artifact `debug/data/lih_pi_m16_energy_trunc1e-2.npz`.
- **MEMORY WALL (the real constraint — RAM, not time).** g_Vee's kept-pair triple field is dense
  ~1.4·nact³·NG·8 bytes. Measured nact(trunc): 1e-2→64 (9 GB), 5e-3→81 (19 GB), 3e-3→100 (36 GB),
  1e-3(full)→196 (**267 GB, infeasible**). On 24.5 GB free the tightest feasible is **5e-3**. The full
  reference (E_full −8.029) is only ~3.9 mHa below the 1e-2 E_trunc, so **feasible M=16 caps ~−8.042
  to −8.043; −8.047 needs the full reference, which is memory-walled.** The clean path to −8.047 is the
  **O(nnz) rewrite of the kept-pair Coulomb** (exploit tensor sparsity, don't materialize nact³) — the
  same O(nnz) trick the σ enumerators already use; engineering, not compute-churn.
- **trunc=5e-3 attempt: KILLED by the OS memory-pressure guard, inside g_Vee at 544s** (the cheap
  pieces finished — E_trunc −8.026118, Fbar 3.498, σ² 0.193, h_Vne 0.319, g_Vne −3.907, h_T −0.388,
  g_T 1.861 — then g_Vee's dense nact³ field crossed the machine's headroom). This **empirically
  confirms the memory wall**: ~19 GB for nact=81 exceeds safe headroom on 31 GB total / ~24 GB free.
  Not restarted (per the guard). **So on this hardware trunc=1e-2 (−8.042181, gated) is the practical
  ceiling.** Tighter needs the O(nnz) g_Vee rewrite or a bigger machine.
- **DELIVERABLE of the π follow-on: E_R12 = −8.042181 (trunc 1e-2), G2-π PASS (0.32 mHa), variational,
  beats the M=16 ceiling.** The π step is confirmed (σ −8.037 → π −8.042), cross-validated. −8.047 is
  memory-walled short; the clean path is the O(nnz) kept-pair-Coulomb rewrite (a focused engineering
  sprint — the σ enumerators already contract O(nnz); only g_Vee's kept-pair Coulomb materializes the
  dense field), OR simply run the existing code at trunc≤3e-3 on a ≥64 GB machine.
- **Owed capture (PI decision point):** fold −8.042181 into Paper 19 [MARRIAGE] + a memo + CHANGELOG at
  the next /checkpoint; and decide whether to build the O(nnz) rewrite to reach −8.047.
- Probe `debug/lih_pi_m16_probe_nact.py` (the nact/RAM curve); driver `debug/lih_pi_m16_energy.py`
  (`--trunc`,`--geminal`); VMC `debug/lih_pi_m16_vmc.py` (`--nwalk/--nsweep/--nburn/--seed`).

2026-09-25: **g_Vee diagnostic DONE -> both algebraic levers STOP, sparsity is the lever.** (`debug/lih_gvee_diag.py`)
- Cost by topology (pi 1e-2): TRI 60% of wall (mode-m Neumann solves, 400 s) but 0.18 GB; **TREE4 owns the dense
  field (7.07 GB = 64^3*NG*8)**; PAIR 9%; BRIDGE 1.5%. |F2Vee|/|g| = 65, no class droppable (TRI ~241 mHa in E).
- Cheap pair-space route exact (8e-9 Ha in E) but only 9%/0.3% of cost -> STOP. gamma-series over exact r12^k pieces
  converges (18-35 terms, float64-safe, sidesteps spheroidal wave functions) but ~100x dearer -> STOP. Both logged in
  `docs/failed_approaches_ledger.md` + CLAUDE.md S3.
- Unique index triples in T: 1,156 at 1e-2, **56,644 at 1e-3 (vs 7.5M dense) -> ~4.3 GB vs 267 GB**, so the O(nnz)
  rewrite should fit trunc 1e-3 (nact 196, E_trunc -8.028633) in RAM.
- **O(nnz) rewrite DONE** (`debug/lih_gvee_sparse.py`, `all_pieces_sparse(gem, refm, T0, mem_cap_gb)`; driver
  `debug/lih_pi_m16_energy_sparse.py`). Symbolic replay of the dense vertex order, evaluation restricted to T's unique
  index projections, memoised compound rows, streaming Coulomb, triangle solves restricted to (cm,L,R) in T, randomness
  replayed in dense order. **R1 PASS** (pi 1e-2: E_R12 -8.042181269 == dense to 1.8e-15 Ha; peak 3.06 GB vs 9 GB;
  701 s vs 1350 s). **R2 PASS** (sigma 74-det -8.0372357, SVD-noise level vs banked). **5e-3 (previously OS-killed):
  E_R12 = -8.042694**, 3.41 GB, 962 s. Guard `tests/test_lih_gvee_sparse.py` (@slow, 4 passed ~12 min; tiny pi ref)
  + `debug/firetest_lih_gvee_sparse.py` 5/5 fired (the kept-pair circulation-selection plant is redundant -- T
  conserves circulation -- documented, cannot fire).
- **trunc 1e-3 DONE (sparse path): E_R12(linexp) = -8.043494** (E_trunc -8.028633, corr -14.86 mHa, c 0.20597,
  Fbar 3.498634), variational, 3.5 mHa from VMC-Pade -8.047. Peak RSS 8.76 GB, wall 8921 s (2.5 h). nact 196,
  T nnz 129,600. Artifact `debug/data/lih_pi_m16_energy_sparse_1e-03.npz` (copied to `lih_pi_m16_energy.npz` as the
  VMC input; the 1e-2 input is preserved as `lih_pi_m16_energy_trunc1e-2.npz`). G2-pi high-stats VMC running.
- **Convergence + reading:** 1e-2 -8.042181 (corr -17.0) -> 5e-3 -8.042694 (-16.6) -> 1e-3 -8.043494 (-14.9).
  The reference gained 3.4 mHa but E_R12 only 1.3 mHa: the richer CI already carries correlation the geminal
  supplied. **The reference is now exhausted as an explanation** (0.41 mHa from the full M=16 CI the VMC used), so the
  remaining ~3.1 mHa to -8.047 must sit in the CORRELATION ANSATZ: ours is linear in F, (1+c(F-Fbar))Psi_CI with one
  fixed r e^{-0.5r} geminal; the VMC used an exponential Pade Jastrow (all powers of u). HYPOTHESIS, not yet measured.
- **Next lever (cheap, existing machinery):** a multi-geminal linear expansion {Psi_CI, (F_g1-Fbar)Psi_CI, ...,
  (F_gk-Fbar)Psi_CI} over several ranges gamma -> a (1+k)x(1+k) generalized eigenproblem. Every matrix element keeps
  <= 3 f-edges (same topology the RI-free reducer already handles exactly); cost ~k^2/2 x one g_Vee on the cross terms.

2026-09-25 (later): **G2-pi on the trunc-1e-3 result.** 8 parallel seeds (nwalk 4000, nsweep 2500, nburn 1000),
inverse-variance combined: E_VMC(linear) = -8.043887 +/- 0.22 mHa vs analytic E_R12 -8.043494 -> d = -0.39 mHa
(1.5 sigma), **G2 PASS**. gate-6 (bare Psi_CI, c=0) = -8.027469 +/- 0.38 mHa vs E_trunc -8.028633 -> d = +1.16 mHa
(3.1 sigma) -- a TENSION (at 1e-2 the same check passed at +0.5 sigma with long burn-in). Leading suspect: the short
seeds used half the original burn-in on a 90-det Psi. Re-running 8 seeds with nburn 5000
(`debug/data/lih_pi_m16_vmc_1e-3_burn_s20*.log`); the original long run (nburn 2000, nsweep 20000) is still going.
Do NOT call the 1e-3 number fully gated until gate-6 is resolved.
- 2026-09-25: first burn-in re-run INVALID (my launch error: `lih_pi_m16_vmc.py` counts --nsweep as TOTAL sweeps
  including --nburn; nburn 5000 > nsweep 2500 collected zero samples and all 8 seeds crashed). Relaunched with
  nwalk 2000, nburn 4000, nsweep 6500 (2500 sampling sweeps). The earlier batch's effective burn-in was 1000 of 2500.
- 2026-09-26: **gate-6 tension RESOLVED as a burn-in artifact.** The original long run (nwalk 4000, nburn 2000,
  nsweep 20000) first pass: gate-6 E_VMC = -8.02833 +/- 0.00030 vs E_trunc -8.028633 -> +0.30 mHa (1.0 sigma) PASS.
  The short seeds (nburn 1000 of 2500) were biased +1.16 mHa. CONSEQUENCE: their G2 linear estimate (-0.39 mHa vs
  analytic) may carry burn-in bias too, so G2 is NOT final until the nburn-4000 seeds' linear pass reports.
- 2026-09-26: **G2-pi on trunc 1e-3 = PASS, FULLY GATED.** 8 seeds (nwalk 2000, nburn 4000, nsweep 6500):
  gate-6 -8.028898 +/- 0.40 mHa vs E_trunc -8.028633 (-0.66 sigma, PASS); linear E_VMC = -8.043631 +/- 0.25 mHa vs
  analytic E_R12 -8.043494 -> -0.14 mHa (-0.56 sigma, PASS). The earlier 3.1 sigma offset was burn-in. Next lever toward
  -8.047: the multi-geminal linear expansion (hypothesis: the gap is the linear-in-F single-geminal ansatz).
- 2026-09-26: **High-stats G2 on 1e-3 PASS** (nwalk 4000, 20k sweeps): gate-6 +1.0 sigma; linear -8.043600 +/- 0.19 mHa vs
  analytic -8.043494 (-0.105 mHa).
- 2026-09-26: **Multi-geminal hypothesis TESTED (VMC linear method, trunc-1e-2 reference): WRONG, better reason.** The gap is
  geminal SHAPE (Kato cusp), not count. linexp0.5 -8.04244; single exp(-2r) -8.04781 (direct VMC agrees -8.047773); Pade
  same ref -8.04729; k=2/3/4 -8.04901/-8.04913/-8.04918 (saturated); F^2 +0.013 mHa (worthless). Files
  `debug/lih_multigeminal_{vmc,analyze,coefs,direct}.py`, `debug/data/multigeminal_vmc.log`.
  NEXT (cheap, gated): analytic marriage at trunc 1e-2 with exp geminal, gamma in {1, 1.5, 2} (~12 min each on the sparse
  path) -- gate = match the already-measured VMC -8.04781 for gamma=2 on the same reference; check NXI/NETA convergence
  (short range). THEN decide whether the 2.5 h trunc-1e-3 run is needed. GAM is a module constant in
  geovac/lih_r12ci/kernels.py. Caveat: linexp c 0.8% (3.8 sigma) off analytic, energy fine -- component check owed.
- 2026-09-26: **Analytic exp-geminal at 1e-2 GO**: gamma 1.5 -8.048138 (best), 2.0 -8.047573 (0.2 mHa vs VMC twin, <=1 sigma;
  grid-converged). **g_Vee has a real analytic defect** (drives the linexp c offset; <=~0.2 mHa in E; suspect F2Vee
  raw moment; unfixed). Memo debug/sprint_lih_exp_geminal_memo.md. NEXT (PI): trunc-1e-3 exp gamma=1.5 run (2.5 h,
  9 GB) and/or the g_Vee term-group audit.
- 2026-09-26: g_Vee triangle-weight defect FIXED (debug/sprint_gvee_defect_audit_memo.md). **trunc 1e-3 exp gamma=1.5:
  E_R12 = -8.049045** (after re-assembly fix), analytic/variational, 1.7 mHa below VMC-over-FCI -8.047. No VMC twin at
  1e-3 yet. Candidate new corpus-best LiH energy (PI to decide on headline after a VMC check).
- 2026-09-27: VMC twin exp g1.5 1e-3 BORDERLINE: gate-6 -0.59 mHa (-1.5 sig); linear -8.049962 +/- 0.19 vs analytic
  -8.049045 (-0.92 mHa, -4.8 sig). Mostly tracks the bare-reference offset (burn-in-type). Resolution option: one long
  high-stat run (nwalk 4000, nburn 2000, nsweep 20000) as for linexp. Headline update in CLAUDE.md/Paper 19 held for PI.
