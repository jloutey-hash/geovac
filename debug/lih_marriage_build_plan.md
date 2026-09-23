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
2026-09-22: plan written; audits A+B complete; nothing built. Next: Phase 0 (one agent).
