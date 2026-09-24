# Build plan — the "marriage": exact RI-free multibody r₁₂ machinery on a multi-determinant reference

**Written 2026-09-22 (v5.15.23), PI-directed, for execution by a FRESH context using subagents.**
Self-contained. Before starting: read `memory/lih_marriage_state_of_play.md` (READ FIRST) and
this file; then the owning paper section (Paper 12 `sec:r12`, the [AUDITED 2026-09-22] paragraph)
+ CHANGELOG after v5.15.23, in case anything moved. Do NOT re-derive the two-halves picture.

---

## THE GOAL
Run the analytic, resolution-of-identity-free explicit-r₁₂ engine (`geovac/lih_r12ci`) on a
**multi-determinant σ reference** instead of the 2-orbital `|1s_A² 1s_B²|` it is welded to, and
obtain a **deterministic, rigorous, variational** LiH energy from it. Two things are won at once:
1. a second, independent, error-bar-free rigorous LiH number cross-checking the VMC −8.047;
2. the first run of the bridge/triangle reductions on **real non-isotropic pair densities under a
   non-block-factorized reference** — the path the 2026-09-22 audit rates as the factorization's
   weakest link. **The marriage IS the test of the factorization's generality.** A GO upgrades the
   "bridge addition" reduction from validated-on-models to validated-on-real-densities; a FAIL at
   G-leaf is the more important finding (the soft spot becoming a wall).

## WHAT IT WILL AND WILL NOT DO (honest cap, set before starting)
The correlation ansatz stays pairwise, F = Σ f(r_ij). The "3- and 4-body" is the INTEGRALS that
products of pairwise geminals generate, handled exactly. So expect **E ≈ −8.045..−8.055, near the
VMC −8.047, NOT −8.062.** The −8.062 additive target needs a nucleus-coupled (r₁₂t²-type, e-e-n)
geminal term — that is the FOLLOW-ON, for which this build is the platform. Do not sell this build
as the −8.062 result.

## CURRENT STATE (verified 2026-09-22)
- Orbital ceilings: −8.032 (M>16) / −8.029 (M=16) / −8.030 (M=17). VMC 2-body Padé: −8.04731 ± 0.00073.
  Analytic r12-CI on the 2-MO reference: −7.942 (exp) / −7.917 (linexp). Exact −8.0705.
- Audit B verdict: NOT a wall. Every reduction takes a separable pair-density product as an argument
  (`_reduce` `gVee.py:106-148`; `contract_raw`; `triangle_raw` `triangle.py:147-159`;
  `neumann_potential`; IBP identities `hT.py:21-25`, `gT.py:7-11`). The reference lives ONLY in
  coefficient lists and closed-form 1s shortcuts (welding map: `energy.py:47-51,70-75,94-103,192-259`;
  `basis.py:37-59`; `kernels.py:43-44,59-62,85-92,120-136`; `hT.py:54-59,75-80`; `hVee.py:40-46,
  84-105,108-152`; `gVee.py:97-103,151-166,191-215,232`; `gVne.py:33-38,70-103`; `gT.py:41-44,72-86`;
  `triangle.py:39,203-219`; `assembly.py:79-94,167-201,204-250,274-281,297,303-305`).
- Audit A caveats that BIND this build: (i) termination bound is `L ≤ 2(l_bridge + l_leaf)`, NOT
  `2·l_bridge` (an L=4 channel is nonzero at 24σ for p-leaves; `r12ci_4e_wall_q2_reducibility.py:99`
  hard-codes the wrong one — do not copy it); (ii) the 1-D radial leaf dressing is the isotropic-leaf
  special case — general leaves need the two-centre Neumann potential (`psi_yuk` = Neumann − smooth,
  which carries a ~0.5% F̄-amplified isotropic bias, `gT.py:95-102`: RE-MEASURE it); (iii) the
  spectator-integration bookkeeping assumed `|Φ₀|² = D(12)D(34)` — with a CI, spectator electrons
  integrate to δ by orthonormality instead (3-index traces of the coefficient tensor).
- Tools that exist: `debug/lih_vmc.py` (real-space Ψ_CI evaluator `psi_ci_full`, analytic orbital
  values/gradients/Laplacians `orbital_vgl`, cached wavefunctions `debug/data/lih_vmc_wf_*.pkl`);
  `debug/fci_fast.py`; `debug/prolate_float_eri.py`; `tests/test_lih_r12ci.py` (6 @slow, the
  2-MO energies — the regression anchor); `tests/test_lih_vmc.py` (5 fast).

## THE BUILD — phased so each phase is ONE subagent task with ONE gate
Each phase: dispatch one agent with the phase text below as its TASK, the listed files, and the
gate as its DECISION GATE; it returns a ≤400-word structured report + a `debug/*.py` driver and a
`debug/data/*.log`. The PM integrates, updates this file's STATUS line, and captures in Paper 12 /
CHANGELOG **before** dispatching the next phase. Never hold more than one phase in the main context.

**Phase 0 — Reference + tensor (no r₁₂ yet).**
Reference = the σ-only core-enriched CI: `lih_core2exp_probe.run(Jb=2,Lb=1,npi=0,core2=(4.5,1.6))`
(E = −8.02335), eigenvector via `lih_vmc.build_lih_wavefunction` / `fci_ground_vector`; truncate
|c| > 1e-3 (~20–60 dets), renormalise; K ≈ 8–10 σ orbitals → K(K+1)/2 pair densities on the
`kernels.py` grid, plus the coefficient tensor T[A,B,C,D] over pair indices and its spectator
traces (3-index). **Gate G0:** grid ⟨T⟩, ⟨V_ne⟩, ⟨V_ee⟩ of the truncated vector vs the ERI-engine
values ≤ 1 mHa. This is ALSO the tight-core grid test (`kernels.py:43-44` was sized for ζ=2.7; the
ζ=4.5 core hit the v5.15.7 grid wall). **STOP if G0 fails** → a graded/atom-centred grid comes first.

**Phase 1 — G-leaf (the factorization's generality test; do this before any energy).**
Pick one REAL non-1s pair density from Phase 0 (core×bond or 2s×bond). Evaluate the 4-body bridge
⟨f₁₂f₃₄/r₁₃⟩ with that density on a leaf by the reduction (general bound 2(l_bridge+l_leaf); Neumann
dressing) and by brute-force MC. **Gate:** agreement ≤ 1e-3 relative (NOT the 5% of the old Q2 check)
AND the L-channel spectrum consistent with the general bound. Also re-measure the `psi_yuk` isotropic
bias on that density. **STOP if it fails** — report it as the factorization's soft spot becoming a
wall; that is a paper-grade negative result (register it in `docs/walls/register.md`).

**Phase 2 — Enumerators over T (the 7 generalisations).**
(1) CI vector → T + spectator traces (Phase 0); (2) F̄, σ² as a 0/2-edge no-operator enumerator
(replaces `assembly.py:188-201`); (3) `FkVne` over T (`:167-176`); (4) `FkSv2` over T with G_pq from
`lih_vmc.orbital_vgl` (`:178-185`, `hT.py:54-59`); (5) `eval_triple/pair_fc/coul_only` over T,
dropping the 2/4 pair-equivalence (`gVee.py:151-215,232`; `assembly.py:297`); (6) Cov[F,Y],
Cov[F,E] via the same 1-edge×kernel enumerator (replaces `cov_FA_FB`); (7) E0 grid-consistent from
the k=0 pieces (replaces `:305`). Rewrite as tensor contractions over pair indices — the current
per-term Python loop (216 × 10⁴ `_reduce` calls) would be hours. **Gate G1 (regression):** feed the
2-MO determinant through T → reproduce −7.917 (linexp) / −7.942 (exp) to < 0.1 mHa.

**Phase 3 — Energy + VMC cross-check.**
Assemble the 2×2 → E_R12. **Gate G2 (the VMC cross-check — this is the rigorous validation):**
Ψ_trial = (1 + c(F − F̄)) Ψ_CI is real-space-evaluable (`psi_ci_full` × a scalar); add a linear
"Jastrow" (1 + cF form; ∇ and ∇² trivial) to `lih_vmc.py`, run VMC at the 2×2's c, and require
E_VMC = E_R12 within ±1 mHa. Also require the variational bound E_R12 > −8.0705.

## DECISION GATE
- **GO:** E_R12 < −8.032 and > −8.0705, G0/G-leaf/G1/G2 all pass → the first deterministic rigorous
  LiH number from the RI-free engine, and the factorization validated on real densities. Capture in
  Paper 12 (upgrade the [AUDITED] tier paragraph) and Paper 19.
- **BORDERLINE:** G2 agrees but E_R12 ≥ −8.032 (the truncation/basis too small) → widen the det
  truncation, re-run Phase 3 only.
- **STOP:** G0 fails (grid) or G-leaf fails (factorization) → report the blocker as the result;
  do not report an ungated energy.

## ANCHORS
−8.032 (must beat) · −8.047 ± 0.001 (VMC; expect to land near it) · −7.917 (G1 regression) ·
−8.0705 (must stay above — frozen falsifier).

## TRAPS / DO NOT RE-DERIVE
- Do not re-assess "is it blocked" — audit B settled it (UNDONE, tractable). Build.
- Do not hard-code `L_term = 2*l_bridge` — that is the s-leaf case (audit A).
- Do not use the 1-D radial leaf dressing on non-isotropic leaves.
- Do not use the gradient-form kinetic estimator in the VMC cross-check (infinite node variance;
  v5.15.22) — the standard local energy with an analytic Jastrow Laplacian.
- Do not promise −8.062 from this build (pairwise ansatz caps ≈ −8.047).
- Do not start Phase 2 at the tail of a long context; Phases 0–1 are cheap, Phase 2 is the big one.

## FOLLOW-ON (after GO): the elegant path to −8.062
Add an r₁₂t²-type (nucleus-coupled) geminal term to F and see whether the same reductions carry it
(the leaf dressings acquire a radial weight; the bridge structure is unchanged). If they do, the
e-e-n correlation the additive-F12 rich basis captured is reachable analytically, RI-free — without
brute-force 3-body QMC. That is the forced/free-seam question worth the effort.

## STATUS
2026-09-22: plan written; audits A+B complete.
2026-09-22 Phase 0 DONE, verdict STOP-on-G0 (⟨V_ee⟩ only). Reference −8.023354 reproduced (M=Mk=12,
4356 dets); |c|>1e-3 in the canonical basis keeps 954 dets, so the reference was rotated to the
natural orbitals of the full CI (E unchanged to <1 µHa): 74 dets, norm 0.999931, E_trunc −8.022918
(+0.435 mHa). Tensor T (78 unordered pairs, 38,786 nonzeros, 22 active pairs) + spectator traces
== 1-RDM/2-RDM to 1.7e-16 (bookkeeping generalisation DONE). G0: ⟨T⟩, ⟨V_ne⟩ agree to 0.0000 mHa;
⟨V_ee⟩ +8.09 mHa, NOT the tight-core density resolution (STO controls exact to 1e-11 incl. ζ=4.5)
but the prolate-Neumann ordered kernel P_l(ξ<)Q_l(ξ>) integrated across its diagonal kink by
cumulative GL sums (Li-1s self-Coulomb +5e-3 rel; NXI 72→144→288: +8.09→+2.03→+0.51 mHa, O(NXI⁻²)).
Artifacts: `debug/lih_marriage_phase0.py`, `debug/data/lih_marriage_phase0_ref.npz` (NO basis),
`debug/data/lih_marriage_phase0.log`. Phase-1 leaf candidate: ρ_(NO0,NO1) core×bond (p-like).
2026-09-23 Phase 0b DONE, GO. `geovac/lih_r12ci/neumann_exact.py` (`ExactNeumann`: barycentric
interpolant of the η-moments through the 72 GL nodes; P-side one 60-pt Gauss, exact for the
≤deg-109 polynomial; Q-side geometric panels ratio 4 in ξ−1; m≥1 via the (ξ²−1)^{m/2} split;
operator (5,35,72,72), 0.29 s build, cached by `kernels.exact_neumann`), opt-in via
`kernels.USE_EXACT_NEUMANN` (read at call time) or `exact=True` on `hVee.neumann_potential`,
`triangle.coul_mode_potential`, `gVee.psi_coul/psi_yuk`. Legacy path bit-identical to HEAD
3d4a279 (27/27 arrays). Anchors: 5ζ/8 at ζ=1.0/1.6/2.6875/4.5 to ≤1.7e-13 (legacy +1.7e-3..
+8.5e-3); two-centre Hartree closed forms ≤6.5e-14; m=0..4 solid-harmonic closed forms ≤2.7e-13;
non-isotropic ρ_(NO0,NO1) mode potentials vs brute-force pointwise ≤1.2e-12. **G0 PASSES:** ⟨T⟩,
⟨V_ne⟩, ⟨V_ee⟩ all within 1e-6 mHa of the engine; total −8.02291830 = E_trunc. Interpolant error
2e-11 at ζ=4.5 (not the 1e-13 the brief guessed). Side change kept: `energy.py` Stage-1 f-integrals
lazy (PEP 562), import 60.5 s → 0.5 s, values bit-identical. G-e: 2-MO `energy()` on the exact
kernel = linexp −7.917199 (−0.38 mHa; dE −29.4 vs VMC −29.5), exp −7.943605 (−1.6 mHa; now 2.4 mHa
BELOW the same-geminal VMC 2×2 −7.94121 — the `gT.py:92-102` cancellation, psi_yuk I_Y[aa,aa] bias
6.6e-3 → 7.8e-7). Banked legacy numbers/tests untouched. Guard `tests/test_lih_r12ci_neumann_exact.py`
(7 passed, 0.9 s) + `debug/firetest_lih_r12ci_neumann_exact.py` (5/5 plants fired). Artifacts
`debug/lih_marriage_phase0b{,_headcheck}.py`, `debug/data/lih_marriage_phase0b{,_headcheck,_ge}.log`,
`debug/data/lih_marriage_phase0_ref_exact.npz` (USE THIS from here on).
**Consequence for G1 (Phase 2):** state the regression per kernel path — legacy path must give
−7.9168 / −7.9420, exact path −7.917199 / −7.943605, each to < 0.1 mHa.
2026-09-23 Phase 1 DONE, GO. **G-leaf PASSES on real non-isotropic leaves** (the whole point of the
build): the exact-RI-free bridge ⟨ρ₁ρ₂ρ₃ρ₄ f₁₂f₃₄/r₁₃⟩ reproduces brute-force MC under the 74-det
CI natural-orbital reference. C0 (isotropic control) PASS 3.3e-7; **C1 (signed p-like core×bond leaf,
bridge=bond) PASS** via the 6-D semi-brute route |reduced−MC|=1.07e-8 vs tol 1.73e-8 (the signed-leaf
12-D floors at σ/|I|=1.5e-3 and agrees at 0.3σ); C2 (directional bond leaf, bridge=core) PASS via 12-D
at 1.8e-4. C1 L-spectrum has a nonzero l=4 channel (0.0017) ⇒ the general bound 2(l_bridge+l_leaf) is
operative, NOT the s-leaf 2·l_bridge (audit A(i) confirmed). Coarse-grid f-kernel dressing errors:
integral-level few×1e-6 (Kf 3.7e-6, Kf2 3.9e-6, psi_yuk 1.6e-9 on the real leaf) — **adequate for
<1 mHa energies; the f-kernel does NOT need the exact-integration treatment the Coulomb kernel got**
(pointwise ~1e-3 is near-nodal only). Guard `tests/test_lih_r12ci_bridge.py` (5 fast, fire-tested 4/4:
isotropic-shortcut and L≤2 plants both fire) — the registered bridge coverage gap is CLOSED. Root-caused
the earlier no-traceback crashes to BLAS thread oversubscription (12 MC workers × multithreaded BLAS),
fixed by pinning *_NUM_THREADS=1. Artifacts `debug/lih_marriage_phase1.py`,
`debug/data/lih_marriage_phase1{,_c2}.log`, `..._phase1.npz`.
Captured: Paper 12 [VALIDATED 2026-09-23] note + [AUDITED] stale-clause fix; Paper 19 energy clause;
walls N=4 entry; claim_test_matrix (bridge row → BACKED, generality row → VALIDATED); memo Phase 1.
2026-09-23 Phase 2 DONE, G1 PASS (exact path). All 7 enumerators generalised to O(nnz) contractions
over the Phase-0 tensor T (`debug/lih_marriage_phase2.py`); on the exact kernel every piece matches the
live `assembly.energy()` oracle to <0.1 mHa (worst 0.015), E_R12 = −7.943536 (exp) / −7.917198 (linexp),
consistent with Phase 0b. Legacy-path pure-f pieces match exactly; the Coulomb/Yukawa offsets are the
legacy assembly's own known inconsistencies, gone on the exact path. **74-det energy NOT reached:**
g_Vee's 216-triple over 55 active pairs (vs 3) exceeds wall-clock on repeated Kmu mode-builds; all other
74-det pieces computed + self-consistent (E_trunc −8.0229 ✓). Machinery proven correct; only g_Vee perf
blocks. Owed minor: stale `assembly.ANALYTIC_REF` g_Vne/g_T/g_Vee (live oracle authoritative).
2026-09-23 Phase 3 DONE — **VERDICT GO. The marriage is complete.** g_Vee productionized (`get_tri_matrix`
per-middle-vertex cache + `coul_mode_stack` batched mode-m Neumann; exact re-check on 2-MO g_Vee <0.001
mHa; 55-pair g_Vee ~1040 s). **Marriage energy E_R12 = −8.0372 Ha** on the 74-det M=12 σ reference
(E_trunc −8.022918): linexp −8.037236 (c +0.20539), exp −8.037088 (c −0.22637), agree 0.15 mHa,
variational, ~5 mHa below the −8.032 ceiling. **G2 VMC cross-check PASS:** linexp E_VMC −8.03755±0.00041
(Δ 0.311 mHa), exp −8.03704±0.00003 (Δ 0.046 mHa); gate-6 c=0 −8.0216±0.0009 vs −8.0229 (+1.5σ). The
first deterministic RI-free variational LiH energy from the many-det explicit-r₁₂ engine, cross-validated
two ways. **Honest: −8.0372 is NOT −8.045..−8.055 and NOT corpus-best (−8.047 VMC stays best).** The
σ-only M=12 reference tops out ~−8.0234 (no π), ~9.4 mHa above the −8.032 ceiling; the ~14 mHa 2-body
correlation the marriage recovers matches the VMC-Padé ~15 mHa → the shortfall is the REFERENCE, not the
machinery. Artifacts `debug/lih_marriage_phase3.py`, `..._phase3.npz`; `lih_vmc.py` linear Jastrow (additive,
fast tests 17/17). Captured: Paper 19 [MARRIAGE 2026-09-23], Paper 12 validated note, matrix, memo Phase 3.
**FOLLOW-ONS (genuine, not cheap):** (a) π-containing M=16 reference → the −8.047 window (LARGE, m-resolved,
general-m triangle; widening M=12 σ can't help, tops out −8.0234); (b) r₁₂t² 3-body geminal → toward −8.062.
Owed at PI /checkpoint: CHANGELOG + version bump + commit; a LaTeX compile of Papers 12/19; the stale
`assembly.ANALYTIC_REF` g_Vne/g_T/g_Vee re-measure.
