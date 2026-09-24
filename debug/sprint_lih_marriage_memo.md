# Sprint memo — the LiH "marriage" (exact RI-free r₁₂ machinery on a multi-determinant reference)

Canonical sprint memo (one per sprint, §13.11). Plan: `debug/lih_marriage_build_plan.md`.
READ-FIRST context: `memory/lih_marriage_state_of_play.md`. Started 2026-09-22 (v5.15.23 base).

## Phase 0 — reference + pair-index tensor + gate G0 (2026-09-22) — STOP on G0 (⟨V_ee⟩ only)

Driver `debug/lih_marriage_phase0.py`; artifact `debug/data/lih_marriage_phase0_ref.npz`
(84 keys, natural-orbital basis); log `debug/data/lih_marriage_phase0.log`.

- **Reference reproduced:** σ-only core-enriched CI, `build_lih_wavefunction(Jb=2,Lb=1,npi=0,
  core2=(4.5,1.6))`: E = −8.023354 (banked −8.02335), M = Mk = 12, 4356 dets, cond(S_o) 5e5,
  R = 3.015 in both the probe and `lih_r12ci/energy.py`.
- **Truncation needed the natural-orbital basis.** |c| > 1e-3 in the canonical (Löwdin-like) MO
  basis keeps 954 dets (no dominant determinant); rotating to the natural orbitals of the full CI
  and re-solving in the same 12-orbital space (E identical to < 1 µHa) gives a leading det
  |NO0² NO1²| with c = 0.987 and **74 dets at |c| > 1e-3**, norm 0.999931, E_trunc = −8.022918
  (engine h1/eri via RDMs == Rayleigh quotient to 2e-15), **E_trunc − E_full = +0.435 mHa**.
- **Tensor T over unordered pair indices** A = (p ≤ q), 78 pairs, ρ_A = φ_pφ_q with the symmetry
  factor in the coefficients; P = |ψ|²/4 normalised, electrons 1,2 α / 3,4 β. T is 6084×6084 with
  only 38,786 nonzeros (22 active α/β pairs; stored COO). Spectator traces T3aab/T3abb (78³),
  T2aa/T2ab/T2bb (78²), T1a/T1b. **Trace-vs-RDM max |Δ| = 1.7e-16**; independent 2×2-minor
  product-operator route max rel 2.4e-15; ∫P = 1 to 1e-12; singlet symmetry 3e-14. The
  spectator-integration bookkeeping the 2-MO code assumed away (|Φ₀|² = D(12)D(34)) is DONE
  for a general CI.
- **G0 table (grid / engine / Δ mHa):** ⟨T⟩ 8.024175 / 8.024175 / 0.0000 (Laplacian and gradient
  forms agree); ⟨V_ne⟩ −20.440798 / −20.440798 / 0.0000; **⟨V_ee⟩ 3.406766 / 3.398680 / +8.086**;
  total −8.014832 / −8.022918 / +8.086.
- **Diagnosis — NOT the tight-core grid wall the plan anticipated.** Per-NO normalisation errors
  ≤ 3.7e-6; STO controls (norm, ⟨1/r⟩, ⟨T⟩, Hartree-dressed self-Coulomb) exact to 1e-11..1e-13
  INCLUDING ζ = 4.5. The whole error is the prolate-Neumann ordered kernel P_l(ξ<)Q_l(ξ>)
  integrated across its diagonal kink by cumulative Gauss–Legendre sums on 72 ξ-nodes: Neumann
  self-Coulomb relative errors +5.0e-3 (ζ = 2.69), +8.5e-3 (4.5), +2.8e-3 (1.6), ~1e-3 (H). The
  (0,0|0,0) Li-1s core-NO term (occ 1.997) carries +8.01 of the +8.09 mHa. Scan on scratch grids:
  NXI 72→144→288 gives +8.09→+2.03→+0.51 mHa (clean O(NXI⁻²)); NETA and LMAX have zero effect;
  raising xi_max at fixed NXI worsens it.
- **Implication for the banked 2-MO energies:** the same kernel underlies `geovac/lih_r12ci`'s
  −7.942 / −7.917, so those sit on a few-mHa Coulomb quadrature error that partially cancels
  (the `gT.py:92-102` note already records a piece-level cancellation). To be quantified once the
  exact operator exists (Phase 0b).
- **Phase-1 leaf candidate:** ρ_(NO0,NO1) = core×bond transition density — largest weight in T
  (2.84); Legendre content about A: a0 = 0, a1 = 0.51, a2 = 0.23 (p-like ⇒ l_leaf ≥ 1, bound
  2(l_bridge + 1)); NO0 is isotropic 1s (a1 = 0.008). Alternatives ρ_(NO1,NO1) (a1 = 0.79),
  ρ_(NO0,NO2) (quadrupolar, a2 = 0.29).
- **Structural note for Phase 2:** T is 0.1% dense — contract pair-space kernels
  K[a,a'] = Σ_AB E[(a,a'),AB] W[A,B] on the 22 active pairs rather than T explicitly.
- Wall: 389 s first run (321 s float64 ERI, cached `debug/data/lih_marriage_phase0_prim.npz`),
  71 s rerun. Nothing in `geovac/lih_r12ci/` or `debug/lih_vmc.py` touched.

**PM decision (2026-09-22):** the fix is an exact ordered-integral Neumann operator on the
EXISTING 72×44 grid, added as a new code path (legacy path bit-identical so the 2-MO regression
`tests/test_lih_r12ci.py` is untouched), rather than NXI ≳ 300 (every NG×NG kernel would grow
~17×, > 1 GB each). Algebraic-first: the densities are entire in (ξ,η) (r_A = a(ξ+η) is linear),
so the only non-smoothness is the kernel's kink, and the partial integrals can be done exactly.

## Phase 0b — exact ordered-integral Neumann operator + G0 re-run (2026-09-23) — GO

(First attempt was killed by an API rate limit after writing the operator and the switch wiring;
the second attempt inventoried and completed it. Drivers `debug/lih_marriage_phase0b.py`,
`debug/lih_marriage_phase0b_headcheck.py`; logs `debug/data/lih_marriage_phase0b{,_headcheck,_ge}.log`.)

- **Design** (`geovac/lih_r12ci/neumann_exact.py`, class `ExactNeumann`): barycentric Lagrange
  interpolant of the raw η-moments g_l(ξ) through the 72 GL nodes; the P-side partial integral
  ∫_1^{ξ_i} P_l^m g_l by one 60-point Gauss rule (the integrand is a polynomial of degree ≤ 109 in
  the mapped variable — exact); the Q-side ∫_{ξ_i}^{ξmax} Q_l^m g_l on geometric panels (ratio 4 in
  ξ−1) toward the log-steep lower endpoint, 60 points; m ≥ 1 through the (ξ²−1)^{m/2} split.
  Operator A(l, m) of shape (5, 35, 72, 72), 7.3 MB, 11,340 sub-points, 0.29 s to build, cached by
  `kernels.exact_neumann(lmax, mmax)`. Opt-in: `kernels.USE_EXACT_NEUMANN` (module attribute read
  at call time) or `exact=True` on `hVee.neumann_potential`, `triangle.coul_mode_potential`,
  `gVee.psi_coul/psi_yuk` (cache keyed by the flag). Measured interpolant error: 1e-15 for ζ ≤ 2.69,
  2e-11 at ζ = 4.5 (c71/c0 = 2e-9) — the brief's "~1e-13" guess was wrong; integrals unaffected.
- **Legacy bit-identity:** 27/27 arrays max |Δ| = 0.0 against the committed HEAD 3d4a279 modules
  (extracted to a scratch package); the earlier pristine capture was itself HEAD-produced.
- **G-a, isotropic 1s self-Coulomb vs 5ζ/8** (exact rel | legacy rel): ζ = 1.0: 3.1e-14 | +1.68e-3;
  1.6: 5.8e-14 | +2.83e-3; 2.6875: 1.0e-13 | +4.96e-3; 4.5: 1.65e-13 | +8.53e-3.
- **G-b, two-centre Hartree closed forms:** (aa|bb) 6.5e-14, (aa|ab) 1.2e-14, (bb|ab) 5.0e-14
  (legacy 5.8e-5 / 2.8e-3 / 3.4e-4). Note `energy.V_aabb…` are the *geminal* (f) integrals and
  `VABAB_REF` is (ab|f|ab), not Coulomb anchors. (ab|ab) Coulomb: exact 0.007752164; legacy-72
  0.007771069; legacy on 400×160 0.007752813; importance-MC 0.007683 ± 4.5e-5 (−1.5σ).
- **G-c, general m:** (i) m = 0..4 solid-harmonic closed forms at both centres ≤ 2.7e-13 (legacy
  2e-3..7e-2). (ii) ρ_(NO0,NO1) × ρ_cyl^m mode potentials: pointwise brute-force (no interpolant)
  1.2e-12 / 1.1e-14 / 8e-14 for m = 0/1/2 (legacy-72: 1e-2 / 5.6e-2 / 1.8e-1); integrated
  Richardson R1(576,1152) 4.6e-8..6.3e-8 (R2 ~1e-10). The brief's pointwise Richardson was
  ill-posed (legacy doubling ratios 1.06–11.7, not 4) and was replaced.
- **G-d, G0 re-run (grid / engine / Δ mHa):** ⟨T⟩ 8.02417516 / 8.02417516 / +1e-6; ⟨V_ne⟩
  −20.44079796 / −20.44079796 / +1e-6; ⟨V_ee⟩ 3.39867962 / 3.39867962 / −1e-6 (legacy +8.0863);
  total −8.02291830 = E_trunc. **G0 PASSES.** Residual diffuse-NO self-Coulomb error −2e-6 is
  xi_max tail truncation, not the kernel.
- **G-e, the 2-MO engine on the exact kernel (324 s):** linexp −7.917199 (−0.378 mHa vs banked
  −7.9168; dE −29.4 vs VMC −29.5 mHa — improved); exp −7.943605 (−1.605 mHa vs banked −7.9420;
  now 2.4 mHa BELOW the same-geminal VMC 2×2 −7.94121). Exactly the ~−7.9435 that the
  `gT.py:92-102` note predicted: the psi_yuk isotropic bias (I_Y[aa,aa] 6.6e-3 → 7.8e-7) was one
  side of a documented piece-level cancellation (g_T −7.8 mHa, g_Vee −2.45 mHa). Banked numbers
  and `tests/test_lih_r12ci.py` untouched; the remaining piece-level errors are in the f-kernel
  dressings (Kf, Kf2) on the coarse grid — to be quantified in Phase 1.
- **Side change kept:** `energy.py` Stage-1 f-dressings/f-integrals built lazily (PEP 562
  `__getattr__`); consumers are two `__main__` blocks only; from-import verified; values
  bit-identical; `import kernels` 60.5 s → 0.5 s (57–61 s deferred to first access).
- **Guard:** `tests/test_lih_r12ci_neumann_exact.py` — 7 passed in 0.89 s (5ζ/8 at ζ = 4.5 and
  2.6875 to 1e-7; the legacy operator asserted > 1e-3 off — the wrong answer pinned; two-centre
  closed form to 1e-7). `debug/qa/fire_test.py` has no registry (convention = `debug/firetest_*.py`):
  `debug/firetest_lih_r12ci_neumann_exact.py` → 5/5 plants FIRED.
- **PM regression check (2026-09-23):** `tests/test_lih_r12ci.py --slow` on the modified engine
  → 6 passed in 335.9 s (`debug/data/regression_lih_r12ci_slow.log`); fast guards
  `test_lih_r12ci_neumann_exact.py` + `test_lih_vmc.py` → 12 passed in 3.9 s
  (`debug/data/marriage_fast_guards.log`). The legacy energy path is intact end-to-end.
- **For Phase 1/2:** use `debug/data/lih_marriage_phase0_ref_exact.npz`; m = 0 dressing
  `hVee.neumann_potential(rho2d, LMAX=34, exact=True)`; general-m
  `triangle.coul_mode_potential(B2d, m, exact=True)` (B must carry ρ_cyl^m); cost 35 matvecs of
  72×72 per potential (78 potentials + G0 in 13 s).

## Phase 1 — G-leaf (the factorization's generality gate) (2026-09-23) — GO (C0 + C1 PASS)

Driver `debug/lih_marriage_phase1.py`; logs `debug/data/lih_marriage_phase1.log` (C0+C1),
`..._c2.log` (C2 follow-on); artifact `..._phase1.npz`. **The bridge factorization HOLDS on a real
non-isotropic leaf under a non-block-factorized (74-det CI natural-orbital) reference.**

- **Retune applied (previous agent's diagnosis).** Gate rule = R1 (few×1e-7 on the C0 closed form,
  100× inside the 1e-3 gate); convergence check = a new coarse rule RC (few×1e-7, ~0.6× R1 cost) —
  R2 (~40× R1 on the 12-prim NO density) retired. Scratch grid 90×54. Cache + artifact banked after
  each dressing/config (nothing lost on interrupt). **Root-caused a hard crash (no traceback) in the
  C1 R1 dressing: BLAS thread oversubscription** — 12 spawn workers × multi-threaded BLAS on 16 cores;
  the coarser RC stayed below BLAS's size threshold and survived, R1's larger matmuls crossed it.
  Fixed by pinning `OMP/MKL/OPENBLAS/NUMEXPR_NUM_THREADS=1` before numpy imports (the driver now sets
  this itself); memory was never the issue (17.7 GB free). RC↔R1 dressing agree to 2.19e-6.
- **C0 (isotropic control, σ-gate config): PASS, carrier 12D.** (i)=2.85624479e-2, (ii)=2.85624570e-2,
  (i)−(ii)=−3.2e-7; 12-D MC (dip+quad, 1.44e8) 2.85621e-2±2.9e-6, |(i)−MC|=3.3e-7 vs tol 2.86e-5;
  6-D bridge 2.85661e-2±5.4e-6 (+0.7σ). Reproduces the σ-gate at 100× tighter tol.
- **C1 (THE gate — signed p-like leaf ρ_(NO0,NO1), bridge ρ_(NO1,NO1)): PASS, carrier 6D.**
  (i)=1.72777754e-5, (ii)=1.72778893e-5, (i)−(ii)=−6.6e-6 (production f-dressing error on I, negligible);
  RC↔R1 on I = 4.1e-6; scratch-grid (90×54) independence 2.2e-8; L-series 34→50 shift 0. **The 12-D
  MC on the signed leaf floors at σ_MC/|I|=1.5e-3 (dipole)/2.8e-3 (dip+quad) even with the exact-mean
  control variate, and both agree with (i) at ±0.3σ.** The **6-D semi-brute bridge** (trusted grid
  dressing interpolated onto MC'd bridge electrons — isolates the Coulomb/Neumann step) reaches
  σ_MC/|I|=5.0e-4 and gives **|(i)−MC|=1.07e-8 vs tol 1.73e-8 ⇒ PASS**. L-spectrum per-l fractions
  0.654/0.335/0.0088/0.00043/**0.00173**/… — nonzero l=4 (audit A(i): channels beyond 2·l_bridge);
  dressed-leaf a0..a4 = 0/0.663/0.070/0.011/0.002 (p-like, from leaf a0..a4=0/0.511/0.232/0.147/0.076).
- **Dressing-error table (production Kf/Kf2/psi_yuk vs exact, 72×44).** Isotropic 1s: Kf `int`
  ≤5e-8, Kf2 `int` ≤4e-7, psi_yuk `int` ≤3e-8. **Real leaf ρ_(NO0,NO1): Kf sup 1.8e-5 / int 3.7e-6;
  Kf2 sup 2.5e-5 / int 3.9e-6; psi_yuk (exact Neumann) sup 5.5e-7 / int 1.6e-9** (legacy
  hT.yukawa_pot_iso vs closed 8e-5..2e-4 — the isotropic bias the plan flagged). **Verdict for Phase
  2/3: the production f-kernel dressing is adequate for <1 mHa energies (integral-level error few×1e-6,
  ≪ gate); it does NOT need the exact-integration treatment the Coulomb kernel got.** The ptw ~1e-3 is
  confined to near-nodal (|Ψ|→0) points and does not enter the integral.
- **Guard + fire test WRITTEN and VALIDATED.** `tests/test_lih_r12ci_bridge.py` (5 tests, 0.91 s, fast):
  C0/C1 production-reduction regressions through the exact Neumann bridge, the C1 gate (reduction vs
  MC within the bar + recorded PASS), and two built-in fire tests (the 1-D isotropic shortcut on the
  signed leaf deviates >gate; the L≤2 truncation deviates >gate). `debug/firetest_lih_r12ci_bridge.py`
  → 4/4 plants FIRED (hVee dead switch, P/Q interchange, and each fire test made vacuous).
- **C2 (positive strongly-directional bond leaf ρ_(NO1,NO1), near-isotropic core bridge ρ_(NO0,NO0)):
  PASS, carrier 12D** (`..._c2.log`, merged into the artifact). (i)=7.59193253e-2, (ii)=7.59189931e-2,
  (i)−(ii)=+4.4e-6, RC↔R1=3.4e-5. Positive leaf ⇒ the 12-D full-brute MC reaches σ_MC/|I|=1.7e-4 and
  gives |(i)−MC12|=1.37e-5=+1.8e-4 rel ⟹ PASS. Leaf a0..a4=1/0.791/0.540/0.340/0.199 (strongly
  directional); L-spectrum 0.467/0.369/0.133/0.0267/0.00265/… (l up to ~4 significant — richer than
  C1, again beyond 2·l_bridge). The 6-D cross-check sits +6.08e-4 (−3.5σ at its 1.7e-4 bar) — a small
  barycentric grid-interpolation offset on the directional leaf, within the gate; the **unbiased 12-D
  is the trustworthy brute force here and confirms (i) at 1.8e-4**.
- **Net verdict: GO.** The exact-RI-free bridge factorization ⟨ρ₁ρ₂ρ₃ρ₄ f₁₂f₃₄/r₁₃⟩ reproduces
  brute-force MC on real, non-isotropic pair densities (signed core×bond and positive bond×core) under
  the 74-det CI natural-orbital reference — the audit's "weakest link" holds. Upgrades the N=4 bridge
  reduction from validated-on-models to validated-on-real-densities. Next: Phase 2 (enumerators over T).

## Phase 2 — the 7 enumerators generalised over T + gate G1 (2026-09-23) — G1 PASS (exact path)

Driver `debug/lih_marriage_phase2.py`; artifact/log `debug/data/lih_marriage_phase2.{npz,log}`;
74-det preview attempt `..._phase2_preview.{py,log}`. Nothing in `geovac/lih_r12ci/` or
`debug/lih_vmc.py` modified (legacy `assembly.energy()` kept as the G1 oracle).

- **All 7 enumerators generalised** to ⟨O⟩ = Σ_ABCD T[A,B,C,D]·contract(O; ρ_A,ρ_B,ρ_C,ρ_D),
  message-passed to pair-space objects (W[A,B]=∫∫ρ_A K ρ_B, 3-index junctions, kept-pair Coulomb)
  then contracted over T's 38,786 nonzeros — O(nnz), never the dense 55⁴. (1) T+traces from Phase 0
  (2-MO T == ¼ DP⊗DP to 0.0); (2) F̄/σ² = 6/36 f-edge diagrams; (3) h/g_Vne = V_ne one-body multiplier
  over the 4 electrons; (4) h/g_T = drift² on the same T with ρ→G_pq (reuse `hT.Gpq`) + geminal PartB/
  gT23; (5) g/h_Vee = kept-Coulomb-pair reducer + RI-free triangle; (6) Cov[F,Y]/Cov[F,E] = mixed
  2-edge enumerator; (7) E0 = grid-consistent ⟨T⟩+⟨V_ne⟩+⟨V_ee⟩+V_NN.
- **G1 = PASS both geminals on the EXACT path (the Phase-3 path).** Every piece within <0.1 mHa of the
  live path-consistent oracle — worst 0.0149 mHa (exp g_T), 0.0016 (linexp g_T); E0_grid = −7.887821
  (0.0015 mHa vs −7.887822). **E_R12 = −7.943536 (exp) / −7.917198 (linexp)** — 0.0024 / 0.0005 mHa
  vs the live oracle, consistent with Phase 0b's −7.943605 / −7.917199 (the 0.069 mHa exp gap is a
  secondary Phase0b-vs-Phase2 exact-path detail, not load-bearing).
- **Legacy path:** the pure-f pieces (σ², h_Vne, g_Vne, g_Vee) match the legacy oracle EXACTLY
  (≤1e-4 mHa) — the enumeration combinatorics are correct; the Coulomb/Yukawa pieces carry the
  documented legacy-Neumann + psi_yuk offsets (worst h_T 4.8/9.6 mHa), an internal inconsistency of
  the legacy assembly (cherry-picked closed forms + hardcoded fine-grid E0), removed by the exact path.
- **Minor finding (owed):** the hardcoded `assembly.ANALYTIC_REF` g_Vne/g_T/g_Vee are slightly stale;
  the live `energy()` oracle matches the generalised values exactly. Re-measure and update those
  reference entries (they only feed the `__main__` self-check, dev<0.5 mHa — not load-bearing).
- **74-det preview: NOT completed (ungated).** The non-g_Vee pieces on the real tensor are computed
  and self-consistent (F̄=1.783, σ²=0.163, ⟨V_ne⟩=−20.441, g_Vne=−3.373, h_T=0.681, g_T=1.562,
  E0_grid→E_trunc=−8.0229 ✓ Phase 0). The full E_R12 is blocked by **g_Vee's 216-triple over 55 active
  pairs** (~275 middle-vertex Kmu mode-builds + many Coulomb solves exceed wall-clock); the sparse/
  vectorized paths are in place and validated (2-MO g_Vee <0.001 mHa), only the mode-builds are the
  residual cost. Cost: G1 ~1140 s (dominated by 4× ~2-min live-oracle calls; the enumerators ~70 s each).
- **For Phase 3:** `all_pieces(gem, ref); assemble(pieces, pieces['E0_grid'])` with
  `KG.USE_EXACT_NEUMANN=True`. The one productionization owed = cache/share the Kmu mode-builds across
  triples so g_Vee over the 55-pair reference is tractable; then the marriage E_R12 computes and the
  G2 VMC cross-check (linear 1+cF Jastrow in `lih_vmc.py`) can run.

## Phase 3 — g_Vee productionization + the marriage energy + gate G2 (2026-09-23) — GO

Driver `debug/lih_marriage_phase3.py`; g_Vee caching added to `debug/lih_marriage_phase2.py`; a linear
Jastrow added to `debug/lih_vmc.py` (additive). Logs `debug/data/lih_marriage_phase3{,_m1,_m2}.log`;
artifact `..._phase3.npz`. No `geovac/lih_r12ci/*` production files changed. Fast guards 10/10.

- **g_Vee productionization (M1 enabler).** `get_tri_matrix` caches the per-middle-vertex triangle
  matrix M_cm[a,b]=Σ_m Σ_rr λ_rr⟨ρ_a v_rr, C_m(ρ_b v_rr)⟩ (shared across all triples through that
  vertex) + `coul_mode_stack` batches the exact mode-m Neumann l-loop over all rank×R-density rows;
  the Kmu builds were already shared via `_MODE_CACHE`. **Exactness re-check on the 2-MO reference:
  g_Vee reproduces the live oracle to <0.001 mHa and E_R12 reproduces the frozen −7.917198 / −7.943536
  exactly**, so the caching is a pure speed change. The 55-pair g_Vee is now tractable (~1040–1050 s),
  which was the sole wall-clock blocker.
- **M1 — the marriage energy (74-det M=12 σ reference, E0_grid = E_trunc = −8.022918):**
  **linexp E_R12 = −8.037236 (c = +0.20539, corr −14.32 mHa); exp E_R12 = −8.037088 (c = −0.22637,
  corr −14.17 mHa).** The two independent geminals agree to 0.15 mHa; both variational (> −8.0705) and
  both beat the M=16 orbital ceiling −8.032 by ~5 mHa.
- **M2 — gate G2 (VMC of the same trial (1+c(F−F̄))Ψ_CI, standard local energy, analytic Laplacian):**
  gate-6 (c=0) −8.0216 ± 0.0009 vs E_trunc −8.0229 (+1.5σ) PASS; **linexp E_VMC = −8.03755 ± 0.00041
  vs E_R12 −8.037236, |Δ| = 0.311 mHa PASS; exp E_VMC = −8.03704 ± 0.00003 vs −8.037088, |Δ| = 0.046
  mHa PASS.** The analytic RI-free assembly and the independent real-space VMC agree.
- **VERDICT: GO** (E_R12 < −8.032, > −8.0705, G2 passes). **The first deterministic, RI-free,
  variational LiH energy from the many-determinant explicit-r₁₂ engine = −8.0372 Ha, cross-validated
  two ways (two geminals + VMC).** This closes the marriage: the RI-free multibody machinery is proven
  correct on real densities (Phase 1) AND runs end-to-end to a real, cross-checked energy (Phase 3).
- **HONEST framing (do NOT over-claim).** −8.0372 is NOT the −8.045..−8.055 the honest cap predicted,
  and it is NOT the corpus-best LiH energy — VMC-over-FCI −8.047 stays best, and −8.062 stays the target.
  The gap is the REFERENCE, not the machinery: the M=12 σ-only 74-det reference tops out at E_trunc
  −8.0229 (the full M=12 σ CI is only −8.02335), ~9.4 mHa above the M=16 ceiling −8.032, because it has
  no π orbitals; the 2-body correlation the marriage captures (~14.2 mHa) MATCHES the VMC-Padé ~15 mHa.
  So the marriage recovers the right 2-body correlation from a deliberately small reference — a
  methodological milestone, not a new best number.
- **Follow-ons (both genuine, not cheap re-runs):** (a) to reach the −8.047 window, the marriage on the
  **π-containing M=16 reference** (136 pair densities, m-resolved kernels, general-m triangle — the
  LARGE case the state-of-play flagged; widening the M=12 σ truncation cannot help, it tops out at
  −8.0234); (b) toward −8.062, the nucleus-coupled **r₁₂t² 3-body geminal** on the same RI-free
  reductions (the plan's elegant follow-on). Cost this phase: M1 2382 s, M2 1925 s (single-thread).
