# Track log — prolate-native (H2-style) LiH: the accuracy route

**Opened:** 2026-09-20 (PI-directed). **PM holds this thread.**

## ===== RESUME HERE — HANDOFF for a fresh context (2026-09-22, v5.15.17) =====

**WHERE WE ARE.** Route C is COMPLETE on geometry. The from-scratch prolate all-electron
LiH **bond length is at experiment: R_eq = +0.2%, π-converged** (v5.15.16, tagged; best in
the corpus — beats composed 5.3%, balanced 8.8%, frozen-core +5.5%). The **energy** is
characterized but NOT closed: the radial-ladder follow-on (v5.15.17) reaches **chemical
accuracy on H2 (99.66%)** but LiH plateaus **~58 mHa above exact** (−8.012 vs −8.070),
capped by the **4-electron determinant wall + core-core correlation** (the single-orbital
analytic core → dense FCI captures ~zero core correlation). NOT radial incompleteness.

**THE NEXT TASK (PI-directed 2026-09-22): lean on the r₁₂ arc for the LiH ENERGY.**
Current-state probe (this session): `geovac/lih_r12ci.energy()` is a **2×2 PoC** —
reference Φ0 = |1s_A²1s_B²| (ionic single-ζ, **E0 = −7.888**, 182 mHa above exact) + ONE
James-Coolidge geminal, → E_R12 = −7.917 (linexp) / −7.942 (exp). **The ceiling is the
REFERENCE, not the geminal** — the geminal pulls a healthy −29 to −53 mHa of short-range
correlation RI-free (incl. the core pair), but sits on a terrible base, so its energy is
actually WORSE than Route C's −8.012. **The path:** F12 synthesis — the r₁₂ geminal
correction on a GOOD base (Route C's −8.012 reference, `debug/prolate_allelectron_c4.py`).

**HONEST SCOPING (so the fresh context doesn't underestimate it).**
- Substantial rebuild: the r₁₂ integrals (Fbar, h_*, g_* in `geovac/lih_r12ci/`) are hard-
  wired to the ionic 2×2 reference; applying the geminal on Route C's multi-orbital reference
  means re-deriving them there + handling F12 double-counting (the geminal must add only the
  correlation the orbital basis MISSES). Slow-mpf (dps=60; minutes–hours/point).
- Realistic ceiling: the geminal adds ~30–50 mHa of the missing (short-range/core-core)
  correlation → LiH plausibly ~−8.04/−8.06 (near-chemical); full <2 mHa on 4e is uncertain.
- **DIAGNOSTIC-BEFORE-ENGINEERING first** (§ memory rule): a design pass to pin the achievable
  ceiling before the full build. Do NOT launch the big mpf build blind.
  **DONE 2026-09-22** (`debug/lih_r12_ceiling_probe.py`, memo `debug/sprint_lih_r12_ceiling_diagnostic_memo.md`):
  the core-core deficit splits ~50/50 into a **cusp half** (geminal-reachable, PLATEAUS at
  ~50% no matter how many r₁₂ terms) and a **radial in-out half** (a 2nd core ORBITAL, NOT a
  geminal d.o.f.). So the LITERAL synthesis (geminal on the single-core-orbital −8.012 base)
  tops out at **~−8.042** (≈25 mHa above exact, ~2×, NOT chemical). **RE-SCOPED: do the cheap
  2nd-core-exponent experiment FIRST** (existing C4 `build_Xtab_s` mixed-exponent code, no new
  integrals) — it banks the radial half (~30 mHa → ~−8.042) AND resolves whether the determinant
  wall admits a 2nd core orbital, the question the whole strategy hinges on. Commit to the
  multi-session two-center 4e mpf geminal build ONLY for the cusp remainder, after step 1.
  **STEP-1 EXPERIMENT DONE 2026-09-22** (`debug/lih_core2exp_probe.py`, data
  `debug/data/lih_core2exp*.log`): adding a 2nd/3rd core STO exponent (tight, ζ≈4.5) to
  Route C's ladder recovers the radial half — **ΔE = −17 mHa (1 partner), −23.5 mHa (2
  partners)**, trending to ~30; robust across valence (moderate-valence confirm −20.1 mHa).
  **Determinant wall RESOLVED favorably:** core orbitals beat valence per unit budget
  (minimal-valence + 2 cores M=8 = −8.018 already BEATS full Route C M=16 = −8.012), so the
  radial half fits within M≤16, and the geminal (cusp half) is a low-rank {Φ0,FΦ0}
  correction that doesn't consume the budget. **BANKED (`bank` mode, M=16, full-π valence):
  3 cores [+4.5,+1.6] → E=−8.02905 (41 mHa from exact), beating Route C's −8.012 by 17 mHa
  at the SAME budget (ΔE −20.1 mHa at matched valence).** Roadmap: Route C −8.012 → core
  enrichment **−8.029** (dense M=16, existing code) → +geminal cusp half (~30 mHa) →
  **~−8.06 near-chemical** → break M=16 (float64/frugal solver) for the last radial+valence
  → chemical. Honest: the FIRST ~20 mHa radial is budget-cheap (done); the FULL radial half
  (~30) is dense-M=16-limited (needs M>16). Next (PI call): (2) F12 geminal build on the
  −8.029 base for the cusp remainder; and/or (1b) float64 productionization to break M=16 +
  make it a sweep. See `debug/sprint_lih_r12_ceiling_diagnostic_memo.md`.
  **PRODUCTIONIZATION DONE (v5.15.19, `debug/{prolate_float_eri,fci_fast}.py`):** sparse FCI
  (bit-exact, 280s→~10s at M=16, UNBLOCKS M>16) + float64 ERI (energy-exact ~5 µHa, ~2×) →
  combined 3× at M=16. **M>16 PUSH (`push` mode, `debug/data/lih_push_m16plus.log`):** the
  orbital lever is NEAR-EXHAUSTED — M=16 −8.029 → M=20 **−8.03170** (+4th core/bond-J3/π-J2 each
  only 1–3 mHa, diminishing). **Pure-orbital LiH ceiling ≈ −8.032 (38 mHa from exact) = new
  corpus best LiH energy** (vs Route C −8.012). The remaining ~38 mHa is the CUSP (orbitals
  can't reach it — ceiling-diagnostic prediction confirmed). **→ the r₁₂ geminal (item 1) is
  the ONLY path to chemical accuracy, target ~30 mHa cusp on the −8.032 base.**
  **ADDITIVE-F12 CAPSTONE (2026-09-22, `debug/lih_additive_f12.py`):** the geminal engine
  (`geovac/lih_r12ci`) is welded to a crude 2-MO ionic reference (E0=−7.888 hardcoded) → can't
  sit on −8.032 without a reference rewrite. So the STANDARD F12 move: add the r₁₂ CUSP
  correction (He-like 2e machinery) to the orbital ceiling. Core-core cusp = **30.1 mHa** (Li²⁺
  Z=3, solid/transferable; He control 26.5 matches known Hylleraas). **−8.032 + core cusp =
  −8.0618 (−8.2 mHa from exact = NEAR-CHEMICAL);** budget pins the valence cusp at ~8 mHa
  (free-H⁻ 15.7 overshoots → bonded valence ~2× smaller) → ~−8.070. **Additive ESTIMATE, not
  variational** — the rigorous number needs the reference-generalization rewrite (deferred).
  **The marriage lands near-chemical.** See `debug/sprint_lih_additive_f12_memo.md`.

**DO NOT RE-DERIVE (closed this session):** the geometry (Route C solved it, +0.2%);
core BREATHING (Phase 0, inert −0.06pp) and core dipole POLARIZATION (Phase 1, −0.07pp,
off by ~50×) — the drift is the frozen-core *approximation*, not core multipole rigidity;
the r₁₂ N=4 wall (already solved — the full 2-center 4e RI-free integrals are assembled;
the limit is the base, not a wall).

**KEY FILES.** Route C engines: `debug/prolate_atomcentered_core.py` (C1), `prolate_mixed_eri.py`
(C2, `build_Xtab_pair`), `prolate_allelectron_c4.py` (C4, `build_Xtab_s`, the −8.012 base +
the +0.2% geometry), `prolate_energy_ladder.py` (v5.15.17, the ladder + the float64
Laguerre×Legendre conditioning fix). r₁₂ arc: `geovac/lih_r12ci/` (the 2×2 PoC). This track
log = the full chronicle; CHANGELOG v5.15.9–17.

**STILL OWED (do at the r₁₂-synthesis close):** (1) `tests/` regression backing the +0.2%
geometry / H2+π gate (a fast reduced-basis proxy — the full runs are ~340 s/pt at dps=60);
(2) Paper 19 sharpen — geometry paper-grade + the mechanism (frozen-core *approximation*, not
rigidity) + energy determinant-walled with the r₁₂-frugal path named; cite the permanent record
(CHANGELOG/tests), NOT `debug/`.

**SPEED.** Everything is dps=60 mpf → slow. The float64 productionization (the re-basing in
`prolate_energy_ladder.py` is a start) turns minutes→seconds and is the enabler for any real
sweeping; scope it if the r₁₂-synthesis needs many points.

## ===== END HANDOFF =====

## Goal
Build LiH the way H2 is built — in the prolate-spheroidal natural geometry with a
**variational, contraction-capable basis** — so the bond orbital can tighten and
the R_eq drift closes. This is the named breach path for the CHEM-ACCURACY /
balanced-LiH wall proven in v5.15.2 (`debug/sprint_balanced_reqdrift_mechanism_memo.md`):
the drift is the finite-extent screening of a *fixed* Z_orb=1 orbital, inexpressible
in the parameter-free construction; a variational prolate basis removes that
constraint by construction.

Scope trick to stay in the proven regime: **freeze the Li 1s² core** (as balanced
already does) → a **2-valence-electron two-center problem** = the N=2 sweet spot
where explicit-r₁₂ works (H2; HeH⁺, v5.14.11). Avoids the N=4 explicit-r₁₂ wall.

## Status
- **PoC (HeH⁺ R_eq)**: heteronuclear 2e variational-prolate gives the right geometry?
  - (2,2) basis: **R_eq = 1.398 vs ref 1.463 bohr = −4.5%** (inward), vs the composed/
    balanced recipe's **+8.8% outward and WORSENING with basis**. Half the error,
    opposite sign, basis-limited character. `debug/data/heh_req_scan.json`.
  - Convergence check at (3,2): **R_eq = 1.399 (−4.4%), stable** vs (2,2)'s −4.5%;
    ENERGY deepens variationally (−2.96704 → −2.96819). So: energy converges, but
    the R_eq residual did NOT shrink with RADIAL (j) basis. Reading: the ~4% inward
    residual is angular/correlation-limited, not radial — and NOT the composed
    recipe's structural OUTWARD wall (which worsens with basis). At R≈R_e the α
    optimum pegged at the grid ceiling 2.0, so a small part of the −4.4% is α-grid
    under-optimization at large R (biases inward); untested whether wider α + higher
    l shrinks it.
  - **α is NOT the lever** (verified): wider α grid [1.4..3.4] → optimum stays low
    (1.4–1.9), R_eq unchanged (−4.7%). Robust to radial + α.
  - **ANGULAR (l) IS the lever — CONVERGENCE DEMONSTRATED (2026-09-20):** (2,3),
    l=2→3, drops R_eq **−4.5% → −0.5%** (1.4556 vs ref 1.4632) and energy
    **11.6 → 0.9 mHa** off the reference. The prolate recipe CONVERGES to the correct
    HeH⁺ geometry AND energy; the residual was heteronuclear *polarization* (angular),
    exactly as physics predicts. **PoC COMPLETE — decisively GO.** Design input for
    LiH: needs adequate ANGULAR (l≥3) basis, not just radial. Contrast: composed
    +8.8% OUTWARD & worsening (structural wall) vs prolate −0.5% INWARD & converging
    (basis rate). `debug/data/heh_req_scan.json` (last run = (2,3)).
- **Scoping**: GO. Pieces identified below; the sparsity objection that walled the
  cheap route does NOT apply (this is an accuracy CI, not a qubit encoding).
- **Increment 1 (core-Hartree screening): PROTOTYPED + VALIDATED** (2026-09-20,
  `debug/prolate_core_hartree.py`). Closed form `V_H(r_A)=(2/r_A)[1−(1+zc r_A)e^{−2zc r_A}]`
  (zc=Z−5/16=2.6875) matches a direct 3D integral to <5e-3 at r_A∈{0.3..4.0}; exact
  two-electron-charge limit (V_H(20)·20=2.0000); r_A→0 → 2zc=5.375. Prolate-basis
  matrix `<g_i|V_H|g_j>` builds via proven quadrature, diagonal positive (repulsive
  screening). The "one new piece" is LOW RISK, empirically confirmed.

## Machinery map (verified 2026-09-20)
- **Engine scaffold**: `debug/prolate_r12_mpf.py` — 2e, 2-center, heteronuclear,
  explicit-r₁₂, prolate (Hylleraas) basis. `assemble_hetero(basis, R, alpha, Z_A,
  Z_B)` (takes arbitrary Z_A → Z_A=3 for Li), `vne_hetero_mpf`, `solve_canonical`
  (float64; **== mpf solve**, verified −2.967037; the mpf *eigensolve* was the n³
  bottleneck — use float64 solve for scans, mpf assembly for the entries).
- **NOT reusable as-is**: `shibuya_wulfman.py` (cross-center V_ne) and
  `neon_core.py` (FrozenCore) are **hydrogenic-basis**, not prolate. The core-shield
  must be built in the prolate basis. (`neon_core` also only covers Z≥11 cores; Li's
  1s² core is currently just the constant E_core=−7.2799.)

## The one new physics piece: frozen Li 1s² core screening the valence
The valence electrons on the prolate grid feel, from the Li side:
1. **Bare −Z_A/r_A = −3/r_A** — already in `vne_hetero_mpf` (Z_A=3). ✓
2. **Core Hartree (repulsive, screens 1→net Li²⁺)** — CLOSED FORM one-body potential:
   `V_H(r_A) = (2/r_A)[1 − (1+ζ_c r_A) e^{−2ζ_c r_A}]`, ζ_c = He-like Li²⁺ 1s
   exponent (≈ Z−5/16 = 2.6875). Net long-range −3/r_A + V_H → −1/r_A (Li²⁺);
   short-range → −3/r_A (valence penetrates). r_A = (R/2)(ξ+η) → a prolate one-body
   integral of the **same class as V_ne** (function of ξ+η). LOW RISK.
3. **Core–valence exchange K** — non-local one-body from antisymmetry with the 1s²
   core; a 2e integral (core-valence-valence-core). MODERATE. Likely small (tight
   core × diffuse valence); test Hartree-only first, then add.
4. **Core–valence orthogonality** — Gram–Schmidt the prolate valence basis against
   the 1s_Li core on the grid, to stop variational collapse into the core. EASY.
   Sparsity is irrelevant here (accuracy CI), so explicit orthogonalization is fine
   — this is exactly what the cheap route could NOT afford.

## Increment 2-3 FIRST LOOK (2026-09-20, `debug/lih_frozen_core_first.py`)
Assembled the frozen-core LiH valence Hamiltonian (assemble_hetero Z_A=3,Z_B=1 →
T+V_ne+V_ee; + the 2e core-Hartree V_H, factorized at p=0 from increment 1's 1e
V_H⊗S). At R=3.015, (2,2), α=1.0:
- **Bookkeeping VALIDATED.** E_tot(R) = E_val + E_core(−7.2799) + 3/R − V_Hdens(R).
  Plugging the *expected* valence energy (−1.122) gives **exactly −8.070** (LiH ref)
  → the R-dependent assembly (V_ne −3/r_A−1/r_B, core-shield V_H, V_NN, core–H
  attraction −V_Hdens) is correct.
- **Projector CONFIRMED NECESSARY (empirically).** Measured E_val = −1.730, i.e.
  **0.6 Ha too low** → E_tot −8.678 BELOW exact −8.070 (non-variational). The
  valence collapses into the −3/r_A core well without a keep-out constraint (the
  core-shield V_H only screens the net charge, not the short-range well). So
  increment 2b (Huzinaga level-shift H→H+λ(P_1s(1)+P_1s(2)), P_1s=|1s_Li⟩⟨1s_Li|,
  needs ⟨g|1s_Li⟩ prolate overlaps) is the required NEXT step; an R_eq from the
  collapsed valence is meaningless and is NOT quoted.

## Increment 2b (Huzinaga projector) + increment 3 (first LiH R_eq) — DONE (2026-09-20)
- **Projector WORKS** (`debug/lih_projector_test.py`): H→H+λ(P_1s(1)+P_1s(2)),
  P_1s=|1s_Li⟩⟨1s_Li| (prolate σ overlaps by quadrature, factorized like V_H). At
  R=3.015, (2,2): E_val rises −1.730 (λ=0, collapsed) → **plateau −1.064 by λ≥100**,
  E_tot −8.012 (VARIATIONAL, above exact −8.070). Collapse removed; orthogonality
  enforced. λ=1000 used (deep in the plateau).
- **FIRST prolate LiH R_eq (`debug/lih_req_scan.py`, (2,2), λ=1000, Hartree-only,
  no exchange/r₁₂): R_eq = 3.065 bohr = +1.7%** vs exp 3.015, E_min −8.0146.
  **vs composed/balanced +8.8% (3.28) & worsening — a ~4-5× improvement at the
  crudest level.** Wider α [1.0..2.0] → +2.2% (E flat in α, −8.01462; a* 1.3→1.6),
  so **robust to α** (α not the lever, same as HeH⁺). Honest headline: **R_eq ≈ +2%
  (3.06-3.08 bohr), Hartree-only, (2,2), robust to α**; exact value coarse-grid-
  limited. Headroom (per HeH⁺): angular l=3 (the geometry lever), core-valence
  exchange, valence r₁₂. **Strategy validated END-TO-END on the real target.**
  NOT yet a paper claim — Hartree-only PoC; needs exchange + r₁₂ + converged basis.
  Files: `debug/lih_frozen_core_first.py`, `lih_projector_test.py`, `lih_req_scan.py`.

  **QUALIFIED 2026-09-20 — the +2% does NOT survive angular refinement, and the
  reason is a confound.** The (2,3) angular scan gives R_eq **+7.1%** (3.228 bohr),
  WORSE than l=2's +2.2% and approaching composed's +8.8% (energy improves, −8.041).
  So the +2% is a small-basis (l=2), p=0 snapshot, NOT a converged result. **BUT the
  LiH scans ran p=0 (NO r₁₂) while the HeH⁺ PoC that converged (−4.5%→−0.5%) had
  p={0,1} (r₁₂ ON)** — confirmed (`grep p_set`). LiH was forced to p=0 because the
  core-shield V_H and the projector are p=0-factorized. So the l=3 divergence is very
  likely the MISSING r₁₂ (a high angular basis over-polarizing with no correlation to
  absorb it), NOT a genuine prolate-LiH failure — but the frozen core (fixed/unpolar-
  ized, the v5.15.2 drift mechanism on the core) is an alternative suspect not yet
  ruled out. **INCONCLUSIVE until the r₁₂-coupled V_H + projector is built.** Do NOT
  quote +2% as the prolate-LiH result; it is unconverged. "Strategy validated
  end-to-end" holds for HeH⁺ (with r₁₂); for LiH it is OPEN pending r₁₂.

## Increment 4 (core-valence exchange) — DONE (2026-09-20, `debug/lih_exchange.py`)
K_ij = ⟨g_i 1s_A|1/r₁₂|1s_A g_j⟩ by prolate quadrature with the φ-integral done
analytically (elliptic K, captures the integrable 1/r₁₂ singularity); K2e = K1e⊗S +
S⊗K1e; core Fock = 2J−K = V_H − K. Result at R=3.015, (2,2): **exchange = −0.9 mHa**
(K1e symmetric, diag>0, correct sign, still variational −8.013). **Small, as scoped**
(valence ⊥ tight core → modest overlap). So the remaining 57 mHa gap is CORRELATION
(r₁₂) + BASIS (angular/radial), NOT exchange; exchange is negligible for R_eq.

## Increment 6 (r₁₂-coupled V_H + projector) — BUILT + VALIDATED; verdict scan RUNNING (2026-09-20)
`debug/lih_r12_coupled.py` rebuilds BOTH p=0-factorized add-ons as r₁₂-coupled 2e
prolate-quadrature matrices, so LiH can run with `p_set=(0,1)` (apples-to-apples
with the HeH⁺ PoC). Route = φ-integrated kernels `Φ_P[a,b]=∫∫dφ₁dφ₂ r₁₂^P`
(P=0→(2π)²; P=1→2π·(R/2)·4√(A+B)·E(m), elliptic E; P=2→(2π)²(R/2)²A) + single-φ
`J_p` for the projector's partial overlap; vectorized over COMBINED powers
(matrix element depends only on jᵢ+jⱼ etc.), so it is fast at l=3 / p={0,1}.
- **V_H (local, summed over electrons):** `VH_ij = Q1_P[(J,L),(K,M)] + Q2_P[…]`,
  `Q_P = (w·V_H·ξ^Jη^L·e^{-2αξ})·Φ_P·(w·ξ^Kη^M·e^{-2αξ})`.
- **Huzinaga projector (NON-local):** `P1_ij=(2π/N₁ₛ)Σ_b w_b f2ᵢf2ⱼ B̃ᵢB̃ⱼ`,
  `B̃ᵢ[b]=Σ_a w_a f1ᵢ[a]·1s[a]·J_{pᵢ}[a,b]` (a `W diag(w) Wᵀ` product), + P2 (elec-2).
- **VALIDATED (`python debug/lih_r12_coupled.py` / scratch val12):** on a p=0-only
  basis BOTH reduce to the existing references — V_H vs `_twoelec_VH` max rel
  **5.6e-4**, projector vs `lih_projector_test.build` **1.2e-4** (diagonals match to
  6 digits; residual = grid-limited near-zero odd-parity off-diagonals). p=0 reduction
  is the frozen falsifier for the formulas.
- **End-to-end single point (l=2, R=3.015, λ=1000, l_neu=20):** p=0 E_tot=**−8.01208**
  (reproduces the track-log baseline exactly), r₁₂-on E_tot=**−8.04149**; r₁₂ lowers
  E_val by **−29.4 mHa** (correlation binds, correct sign), stays variational
  (>−8.070), 162/162 kept. Pipeline consistent (mpf S/T/V_ne/V_ee + float64 V_H/proj).
- **VERDICT — FROZEN CORE IS THE CULPRIT (PI's bet RIGHT), 2026-09-20.**
  `debug/lih_r12_req_scan.py`, p_set=(0,1), λ=1000, l_neu=20, α∈{1.0,1.4},
  R∈{2.70,2.85,3.015,3.20,3.45}. Logs `debug/data/lih_r12_req_l{2,3}.log`:

  | | l=2 | l=3 | trend |
  |---|---|---|---|
  | LiH r₁₂ OFF (p=0) | +2.2% | +7.1% | diverges outward |
  | **LiH r₁₂ ON** | **+3.4%** (3.118 bohr) | **+5.5%** (3.180 bohr) | **still diverges outward** |
  | HeH⁺ r₁₂ ON | −4.5% | −0.5% | converges |

  WITH r₁₂ the R_eq drift STILL grows with angular basis and STILL points outward —
  the OPPOSITE of HeH⁺'s convergence at the identical (2,3) basis. r₁₂ absorbs part of
  the over-polarization (drift magnitude +7.1%→+5.5% at l=3, energy −8.041→−8.055,
  correlation ~14 mHa) but does NOT restore the HeH⁺ pattern. Energies variational
  throughout. **λ-plateau confirmed** (l=2, R=3.015, r₁₂: E_val flat −1.0965→−1.0953
  across λ=30→10⁴; λ=0 collapses to −1.888/E_tot −8.836) → the outward drift is NOT a
  projector artifact. **Solve health (l=3, R=3.20):** cond(S)=6.4e12 (≪ float64 1e16),
  float64 canonical sweep keeps 288/288 vectors (no subspace dropped), and a full mpf
  cross-solve AGREES to **0.000 mHa** (both E_val=−1.087800) → the float64 scan geometry
  is exact for this purpose. The ONE variable differing between converging HeH⁺ and drifting
  LiH at identical basis is the **frozen core** → clean controlled comparison. The
  fixed/unpolarized Z_orb=1 core is the v5.15.2 finite-extent-screening mechanism ON
  THE CORE; r₁₂ (a valence-valence correlation lever) cannot reach it. **The cure is an
  all-electron / variational-core LiH = the N=4 build** (hits the N=4 explicit-r₁₂ wall
  for full correlation; the geometry, not full correlation, is what N=4 must fix).
  Honest baseline: prolate frozen-core (+5.5% @ l=3) improves on the BALANCED recipe
  (+8.8%, Paper 19) but does NOT beat COMPOSED l-dep-PK (5.3% @ l_max=2, Paper 17) —
  a mechanism diagnostic, not a new best LiH R_eq.

## Increment 7 (all-electron / variational-core FCI) — ENGINE BUILT + VALIDATED; verdict scan RUNNING (2026-09-20)
The increment-6 verdict said the cure is a variational CORE = all 4 electrons active,
R-adaptive, no frozen core. Built as a multi-exponent prolate FCI reusing the
validated `geovac/prolate_scf.py` grid machinery (`debug/prolate_allelectron_fci.py`):
- **1-particle basis:** one-electron eigen-MOs from `get_orbital_on_grid` at SEVERAL
  (Z_eff_A, Z_eff_B) scales (spanning the tight Li core Z≈2.9 → diffuse valence Z≈0.7),
  Löwdin-orthogonalized. h1 across different-exponent orbitals via
  `T φ_q = ε_q φ_q + A(gen_q)φ_q` (kinetic → computable multiplicative attraction
  matrices `_attraction_matrix`). ERIs via `compute_vee_integral` (elliptic-K).
- **FCI:** spin-orbital Slater-Condon, Sz=0 block, V_ee PAIRWISE → **sidesteps the N=4
  explicit-r₁₂ 4-body wall** (that wall is Hylleraas-specific, not CI).
- **VALIDATED:** (a) HF-level single determinant reproduces `eckart_scf_energy` EXACTLY
  (−1.05264 at H₂ Z=1, ε/J₀₀ identical) → FCI + integral wiring correct. (b) The bare
  single-exponent MO basis is radially inflexible (H₂ stuck at −1.067); the
  MULTI-exponent basis converges: H₂ −1.067→**−1.134 (proper HF)**→−1.142→−1.143 as
  scales are added (residual ~30 mHa is σ→π angular correlation, omitted; σ sets the
  geometry). (c) LiH all-electron FCI converges variationally toward −8.07:
  −7.671→−7.727→−7.746 (Nsc=2/3/4, σ-only, M=8; basis-limited ~0.32 Ha, dominated by
  R-INDEPENDENT core-correlation, so R_eq is still the meaningful read-out).
- **VERDICT (σ-only): INCONCLUSIVE — the σ-only all-electron PES COLLAPSES INWARD, and
  the missing lever is ANGULAR (π), 2026-09-20.** `debug/lih_allelectron_req_scan.py`
  (log `debug/data/lih_allelectron_req.log`): no minimum in R∈[2.60,3.45]; E_tot is
  MONOTONE decreasing toward small R (Nsc4: −7.759@2.60 → −7.698@3.45). The electronic
  R-slope is too STEEP (+0.41 vs the +0.33 = Z_AZ_B/R² needed to balance V_NN) →
  inward collapse. This is the MIRROR of the frozen-core failure (which was too-weak
  slope → outward drift). Nsc3 and Nsc4 agree on the over-steep slope (basis-stable).
- **Diagnosis (`debug/data/lih...` diffuse test):** the diffuse-H⁻ (ionic Li⁺H⁻)
  hypothesis is directionally RIGHT but far too weak — adding diffuse scales
  (M=8→12) lowered E slightly more at large R, flattening the slope +0.0370→+0.0335
  (~10% of the +0.037 needed for a min near 3.0), while cond(S) blew up 1e4→3e6
  (σ linear dependence walls out further diffuse refinement). So the residual
  over-steepness is NOT a diffuse-basis deficiency; it is the missing **angular (π)
  polarization** — the same lever that was decisive for HeH⁺ (−4.5%→−0.5% at l=2→3)
  and the frozen-core study.
- **What this does and does NOT settle.** Engine BUILT + VALIDATED (HF == eckart
  exactly; multi-exp H₂ converges −1.067→−1.134→−1.143; LiH FCI variational). But the
  σ-only truncation introduces its OWN (inward) basis error, so this run neither
  confirms nor refutes "variational core cures the drift" — it is INCONCLUSIVE pending
  π. The increment-6 conclusion (the frozen core IS a real problem) stands; the cure's
  demonstration is not yet in hand.
- **π (generalized-m) ERIs BUILT + VALIDATED, and they do NOT fix it (PI-directed,
  2026-09-20).** Added the Fourier-resolved azimuthal Coulomb kernel
  `K_μ = 2π·F_|μ|(a,b)` (Cohl–Tohline toroidal expansion; `F_0=4K/√(a+b)`,
  `F_1=(4/b)[aK/√(a+b)−√(a+b)E]`, `F_2` by the toroidal recurrence) + the M_L selection
  `m_p+m_r=m_q+m_s` (`_azimuthal_kernels`, `vee_m`, `build_mo_integrals_full`).
  VALIDATED: μ=0 reproduces `compute_vee_integral` to 1e-10; μ=0,1,2 match direct
  numerical Δφ integration to 1e-13; H₂+π lowers E correctly (−1.133→−1.144 toward
  −1.1745). Guards `tests/test_prolate_allelectron_fci.py` (3 pass, fire-testable).
  **But π does NOT flatten the LiH inward slope** — σ-only +0.017 → +π(2) +0.019 →
  +π(4) +0.019 (marginally MORE inward). Angular is not the missing lever here.
- **THE REAL WALL (both diffuse-basis and π ruled out): basis convergence + conditioning.**
  The all-electron multi-exponent FCI is 0.4–0.45 Ha ABOVE the −8.07 reference at
  achievable M, and cond(S)→1e6 as scales are added (σ linear dependence). A method
  0.4 Ha from exact cannot locate a bond (D_e≈0.08 Ha), so the inward collapse is a
  R-dependent basis-incompleteness (BSSE-like) artifact, NOT physics — and it is not
  cured by the angular channel. **A clean all-electron LiH R_eq verdict is NOT reachable
  with this grid-orbital multi-exponent FCI at practical accuracy/conditioning.**
- **What IS established (honest scope):** (i) the engine + π machinery are validated,
  reusable deliverables; (ii) the increment-6 diagnosis stands (frozen core = the
  rigidity that caused the +5.5% OUTWARD drift); (iii) the two failures BRACKET the
  true answer with opposite sign (frozen-core rigid-core → outward; all-electron
  variational-core-but-basis-limited → inward), and standard all-electron QC is known
  to give good LiH R_eq (~3.0 bohr) in adequate bases — so the variational core removes
  the frozen-core rigidity, but this particular engine can't demonstrate it cleanly.
- **Routes to an actual verdict (PI call):** (a) multi-orbital prolate SCF (HF orbitals
  converge far faster per-function → likely a sensible R_eq; needs non-local exchange on
  the grid — substantial); (b) counterpoise-correct the FCI PES (addresses the BSSE
  artifact directly); (c) accept the bracketing + standard-QC argument as the soft
  verdict and stop.

## Review cycle (PI-directed: "does the cancellation idea spark?") — mechanism PINNED to a NUMERICAL grid wall (2026-09-20)
PI hunch: maybe H₂/HeH⁺ had similar error cancellations worth exploiting. Reviewed;
the hunch was productive — it pinned the mechanism (which an earlier guess mis-called
BSSE) and tightened the bracket, though it is not a clean fix.
- **The grid-FCI engine is SOUND.** H₂ R_eq in the SAME engine that fails for LiH gives
  a clean minimum EXACTLY at R=1.40 (−1.083/−1.124/**−1.133**/−1.121/−1.088 over
  1.0–2.1), despite ~0.04 Ha basis error. So the inward collapse is NOT a universal
  engine bias. `debug/data/h2_gridfci_req.log`.
- **The mechanism is the tight Li core's R-DEPENDENT GRID RESOLUTION, not BSSE.** An
  ISOLATED Li²⁺ 1s² (physically R-independent) swings **0.36 Ha** across R in this
  engine, and a ghost-H basis makes ZERO difference (`lih_core_bsse_probe.log`) → not
  basis borrowing. It is r_A=(R/2)(ξ+η): the uniform ξ-grid resolves the tight core
  better at small R. **Confirmed a grid artifact:** the core R-spread SHRINKS with
  refinement (0.098→0.062→0.049 Ha at N_grid=40→60→80) and E→ true −7.28
  (`lih_core_grid_refine.log`). The 0.36/0.05–0.1 Ha core artifact swamps the 0.088 Ha
  bond → the inward collapse. H₂ has no tight core → immune.
- **Fragment/cancellation correction (the PI's idea) OVER-corrects.** E_CP(R) =
  E_LiH(4e; Z=3,1) − E_core(Li²⁺ 2e; Z=3,0) on identical basis specs (identical
  orbitals/grid/ERIs; only V_ne differs) flips inward collapse → **+22% OUTWARD**
  (R_eq 3.69). The isolated core over-estimates the artifact present in the screened
  in-molecule core → not cleanly separable. `lih_fragment_corrected.log`.
- **NET (tighter bracket, mechanism pinned):** raw all-electron → inward collapse
  (min<2.6); fragment-corrected → +22%; frozen-core (incr 6) → +5.5%; truth 3.015 is
  inside. The clean number is blocked by a **numerical wall** (R-tied core grid
  resolution, converges too slowly on a uniform grid), NOT by physics. The physics
  conclusion stands (frozen core = rigidity culprit; variational core = cure, bracketed
  + standard-QC). **Real fix = a graded grid concentrating resolution at the nuclei**
  (or an analytic/exact-integral core), a well-defined numerical target — not a new
  method. Guards `tests/test_prolate_allelectron_fci.py` (π-ERI machinery).

### Graded-grid test cycle (PI-directed) — MARGINAL; the core R-dependence is co-dominated by BASIS completeness, not grid resolution (2026-09-20)
Built a graded prolate quadrature grid (`debug/prolate_graded_grid.py`,
`get_orbital_graded`): a COMPOSITE two-region Gauss rule with an R-adaptive fine core
panel [1, 1+C/(ζR)] following the tight-core decay e^{−ζR(ξ−1)}.
- **Real bug found + fixed (engine hardening, keeper):** `_attraction_matrix` evaluated
  Z/r numerically near the focus (ξ→1,η→−1, r→0) → −30 Ha garbage on any focus-sampling
  grid. The 1/r singularity CANCELS the Jacobian analytically: (Z_A/r_A+Z_B/r_B)·J =
  (R/2)²[Z_A(ξ−η)+Z_B(ξ+η)]. Rewritten to the analytic form → grid-robust. Identical on
  the stock grid (which never samples the focus), so all prior results stand; guards pass.
- **The grading itself is only MARGINAL.** Isolated Li²⁺ core spread over the bonding
  range [2.6,3.75]: stock 0.081 Ha → graded 0.069 Ha (~15%), still ≈ the 0.088 Ha bond.
  And brute uniform refinement (N_grid=80 → 0.049) beats the graded N=48 — so a simple
  composite grid is not even the best use of points.
- **Why grading can't crack it:** the core's R-dependence is ~HALF grid-integration
  (partly fixable by refinement/grading) and ~HALF **two-center-basis completeness for a
  tight atomic core** (the Coulomb-Sturmian-like orbitals at fixed ζ span the true 1s²
  with R-dependent fidelity — a BASIS issue no grid touches). Neither half alone is <<
  the bond. The clean LiH number stays blocked; a proper fix needs a better tight-core
  BASIS (atom-centered/analytic core), not just a grid — a bigger change than a grid tweak.
- **Verdict unchanged:** physics conclusion stands (frozen core = culprit; variational
  core = cure, bracketed + standard QC); the bespoke clean number is a well-characterized
  numerical wall with now TWO named halves (grid + tight-core basis completeness).

## Next increments (remaining) — reordered after the l=3 divergence
5. **DONE: angular l=2→3 (p=0) — DIVERGED (+7.1%), confounded by missing r₁₂.** Not
   the clean HeH⁺ convergence; see the QUALIFIED note above.
6. **DONE (machinery): r₁₂-coupled V_H + projector built + validated (see section
   above); the load-bearing R_eq verdict scan is RUNNING.** Answers: does prolate LiH
   converge (r₁₂ was the lever) or drift (frozen core is the culprit)?
7. If (6) still drifts → the frozen core is implicated (v5.15.2 mechanism on the
   core); the fix would be an all-electron (variational-core) LiH = the N=4 build.
   Force-decompose (v5.14.8-style) to attribute the drift to fixed-core vs valence.
8. Finer R-grid + total-energy accounting → converged prolate LiH R_eq for a paper claim.

## Build increments (each validated before the next)
1. **Core-Hartree one-body term** `V_H` in the prolate basis (quadrature first,
   Neumann/closed-form later). Validate: net (−3/r_A + V_H) matrix → screened Li²⁺;
   single-center matched-exponent diagonal vs the analytic two-1s Coulomb value.
2. **Core–valence orthogonalization** of the valence basis vs 1s_Li. Validate: no
   below-core collapse; valence 1s-like state removed.
3. **First LiH valence energy + R_eq** (Hartree + orthogonality, no exchange, no
   r₁₂). Validate: R_eq near experiment (1.595 Å = 3.015 bohr); compare to the
   balanced +8.8% drift — expect a much smaller, basis-convergent residual.
4. **Add core–valence exchange K**; measure its effect on E and R_eq.
5. **Add valence r₁₂** (the HeH⁺ engine already has it) → chemical-accuracy energy
   at the correct geometry.

## Reference values (verify before quoting)
- LiH experimental R_e = 1.595 Å = **3.015 bohr**; D_e ≈ 2.43 eV; total E ≈ −8.070 Ha.
- HeH⁺ R_e = 1.4632 bohr, E ≈ −2.97869 Ha (used for the PoC).

## The discriminating test (r₁₂), and the PI's hypothesis (2026-09-20)
**PI's bet: the FROZEN CORE is the culprit** (not the missing r₁₂). The r₁₂-coupled
build discriminates: run LiH WITH r₁₂ (apples-to-apples with the HeH⁺ PoC that
converged). If R_eq then converges toward the reference → r₁₂ was the issue (PI's
bet wrong). If it STILL drifts outward with angular basis → the frozen core is the
culprit (PI's bet right; the fixed/unpolarized core = the v5.15.2 drift mechanism on
the core → the fix is an all-electron/variational-core LiH, the N=4 build). Either
verdict is a clean atlas result (forced/free/WALL). Build route: r₁₂-coupled V_H +
projector by 2e quadrature (adapt `heh_probe.vne_hetero_quad`: swap V_ne for V_H /
the 1s-projector density; the engine already gives S/T/V_ne/V_ee at p={0,1}).

## Build-architecture wrinkle found while scoping increment 3 (2026-09-20)
The HeH⁺ engine (`assemble_hetero`) is a **2-electron Hylleraas** (explicitly
correlated) basis, NOT an orbital-product CI. Consequences for frozen-core LiH:
- **Core-Hartree V_H** (increment 1) adds fine — a one-body potential summed over
  both electrons, ⟨g_i|V_H(r1A)+V_H(r2A)|g_j⟩, same structure as the V_ne sum
  `vne_hetero` already builds. Needs a V_H analog of `vne_hetero_mpf`. Straightforward.
- **Core–valence orthogonality is the real wrinkle.** In a Hylleraas 2e basis there
  is no clean "valence orbital" to Gram–Schmidt against the 1s core. Without it the
  valence variationally COLLAPSES into the core (the net Li potential −3/r_A+V_H is
  still −∞ at r_A→0). Fix = a **Huzinaga level-shift projector** H→H+λ(P_1s(1)+P_1s(2)),
  P_1s=|1s_Li⟩⟨1s_Li|, needing ⟨g_i|1s_Li⟩ one-electron overlaps (computable). This is
  EXACT orthogonality (not the crude PK pseudopotential), and sparsity is irrelevant
  for an accuracy CI — so it is the *improvement* over PK, not a re-run of it.
- **Two build paths:**
  - **(A) frozen-core 2e Hylleraas** (reuse HeH⁺ engine): + V_H + Huzinaga projector
    + exchange. Stays 2-electron (sweet spot). The wrinkle is the projector. RECOMMENDED.
  - **(B) all-electron 4e prolate orbital-CI**: orthogonality automatic, but needs a
    NEW 4e prolate CI engine (doesn't exist) and hits the N=4 r₁₂ wall for full
    correlation. Bigger. Not recommended as the first build.
- **This is a multi-session build** (V_H-in-Hylleraas, projector, exchange, R-scan),
  each piece validated. NOT a one-turn job.

## Open risks
- The −4.5% HeH⁺ residual must be shown to CONVERGE with basis (gate running). If it
  stalls, the prolate route has its own (smaller) geometry limit — still better than
  composed, but not "solved".
- Core–valence exchange magnitude unknown; may or may not matter for R_eq.
- Guardrail note: two-center prolate is NOT single-center (Papers 8-9) nor
  concatenation (FCI-M) nor nested (Track DF) — none block this; the walls audit
  already scoped the composed/PK negatives as MIS-SCOPED for the prolate route.

## RELAXABLE-CORE ROUTE — Phase 0 (core BREATHING) = NO-GO, decisive (2026-09-21, PI-chosen)
PI chose the relaxable-analytic-core route: keep the clean frozen-core ANALYTIC
representation (no grid wall — the increment-7 blocker) but let the core RELAX.
Phase 0 = the cheapest relaxation: let the 1s² core BREATHE (one variational knob,
the exponent ζ_c), staying spherical. `debug/lih_r12_breathing.py` (l=2, λ=1000,
l_neu=20; assembles the zc-independent mpf S,H once per (R,α) and sweeps ζ_c cheaply).
- **Wiring validated:** fixed-ζ_c reproduces the increment-6 baseline BIT-FOR-BIT
  (E_tot=−8.04149 @ R=3.015; R_eq(frozen)=3.118=+3.42%, matches the +3.4% l=2 verdict
  log). So the R_eq *shift* is the breathing effect, not a code artifact — the grid/
  basis offset cancels in the fixed-vs-free comparison (identical assembly).
- **VERDICT: NO-GO — breathing is INERT.** ζ_c*(R) is flat ≈2.70 across R∈[2.70,3.45]
  (2.700→2.696, a whisker off the isolated 2.6875), worth 0.07–0.16 mHa; R_eq moves
  3.1181→3.1163 bohr = **−0.06 pp** (+3.42%→+3.36%). ζ_c* trends the *right* way
  (more contracted at small R) but the magnitude is ~2% of the drift. Radial core
  SIZE is decoupled from the bond (tight, stiff core: E_core curvature 2; valence
  lives far from it). Data: `debug/data/lih_breathing_l2.log`.
- **Reading — confirms the v5.15.2 mechanism's own prediction:** the core must change
  SHAPE (a bond-directed dipole = screening that FOLLOWS the bond), not SIZE.
  → **Phase 1 = core POLARIZATION** (anisotropic dipole screening + harmonic core-
  response cost, one dipole knob d, minimized at each R like ζ_c). New machinery:
  the anisotropic V_H^{dip}(r_A)·cos θ_A term is a NEW integral class (cos θ_A =
  z_A/r_A; z_A=(R/2)ξη, r_A=(R/2)(ξ+η)), NOT the existing (ξ+η) one-body class.

## RELAXABLE-CORE ROUTE — Phase 1 (core POLARIZATION) = WEAK/NO, decisive (2026-09-21)
Built the polarizable analytic core: φ_core = 1s(ζc) + λ·2p_z(ζc), the dipolar
1s·2p_z cross density and its closed-form l=1 Hartree potential
V_H^dip(r_A)=(4π/3)[r_A⁻²∫₀^{r_A}u r'³ + r_A∫_{r_A}^∞ u], u=4λ(ζc⁴/π)r e^{−2ζc r}.
`debug/prolate_core_dipole.py` — VALIDATED (exact dipole tail d=4/ζc to 6 digits;
large-r vs direct 1e-3; singular small-r reference CONVERGES to closed as n grows).
Model polarizability α_c=9/ζc⁴=0.1725 (parameter-free; true Li⁺=0.1925). Driver
`debug/lih_r12_polarization.py`: anisotropic screening built by REUSING the validated
`vh_coupled` with its pointwise potential swapped monopole→V_H^dip(r_A)·cosθ_A,
cosθ_A=(ξη+1)/(ξ+η); E_tot(R,d)=e_val[H+VH_mono+d·VHDIP+proj]+ZZ/R−V_H_mono(R)+E_core
+d²/(2α_c)−Z_H·d·V_H^dip(R); minimized over d at each R (d=0 ≡ frozen baseline, bit-
for-bit −8.04149). l=2, λ=1000, l_neu=20. Data `debug/data/lih_polarization_l2.log`.
- **VERDICT: WEAK/NO.** R_eq: frozen 3.1181(+3.42%) → pol@model-α 3.1162(+3.36%, −0.06pp)
  → pol@true-α 3.1160(+3.35%, −0.07pp). d*(R) FLAT ≈−0.013 across R∈[2.70,3.45] (so no
  PES tilt → no R_eq shift), gain ~0.5 mHa. **true-α ≈ model-α → NOT a model-fidelity
  limit.** To close +3.4% by polarization alone would need ~50× the physical core α
  (~9 a.u., a valence-scale value) — unphysical by a wide margin.
- **Physical reason (clean, worth keeping):** LiH is IONIC (Li⁺H⁻). The valence sits
  on H; its repulsion pushes the core *away* from H (∂e_val/∂d≈+0.18) *harder* than the
  bare H nucleus pulls it toward H (drive≈0.11) → net d*<0 (core polarizes AWAY from H⁻)
  and small. The ionicity that DEFINES LiH starves the core of a polarizing field, and
  the Li⁺ core is intrinsically stiff (α=0.19). Core dipole polarization cannot be the
  lever for ionic LiH.
- **JOINT PHASE-0+1 CONCLUSION — the relaxable-analytic-core route is CLOSED.** Neither
  natural analytic core-relaxation multipole (breathing −0.06pp; dipole polarization
  −0.06/−0.07pp) closes the drift; both are ~R-independent (offset, not tilt). The
  increment-6 "frozen core is the culprit" (by elimination) stands as an OBSERVATION,
  but the *cure hypothesis* — "let the analytic core relax" — is REFUTED: the residual
  is NOT core-multipole rigidity (off by ~50×). It lives in the frozen-core APPROXIMATION
  STRUCTURE (Huzinaga projector hardness / fixed spherical screening shape / absent
  core-valence correlation), which only a full all-electron treatment removes — = the
  increment-7 route, blocked by the tight-core numerical/basis wall. **§3 dead-end
  candidate; PI fork: bank the soft verdict (Paper 19 + this insight) vs. commit to the
  all-electron + atom-centered-tight-core build.** Deliverables kept: `prolate_core_dipole.py`
  (validated dipole-Hartree, reusable), `lih_r12_polarization.py`, `lih_r12_breathing.py`.

## ROUTE C (all-electron, atom-centered ANALYTIC tight core) — PI-chosen 2026-09-21
Goal: the clean all-electron LiH R_eq that increment 7 (all-electron on the grid)
could not reach — its tight Li core swung 0.36 Ha across R (R-tied grid resolution,
r_A=(R/2)(ξ+η)). Cure: compute every core-involving integral ANALYTICALLY (R-accurate
by construction) with the tight core atom-centered, all 4 e⁻ active (no frozen core,
no projector) so core-valence correlation + full core relaxation are IN. Stays in the
prolate natural geometry (NOT LCAO — guardrailed). Reuses the analytic engine
(`prolate_recondition.py` S/T/V_ne, `neumann_vee_general_m.py` general-m V_ee) and
increment-7's integral-source-AGNOSTIC FCI solver `fci_energy(h1,eri,M,nelec)`
(consumes spatial/orthonormal h1[M,M] + chemist eri[M,M,M,M]; add V_NN yourself).
Machinery survey: subagent 2026-09-21 (reusability map + the η-exponential obstacle).

**The obstacle (named by the survey): the η-EXPONENTIAL.** A Li-centered 1s STO
e^{−ζr_A}=e^{−ζ(R/2)ξ}·e^{−ζ(R/2)η} carries an η-exponential the analytic engine's
basis (ξ^j η^l e^{−αξ}, polynomial-η) lacks. Every core-involving integral reduces to
ONE new primitive M_η(q,β)=∫_{−1}^1 η^q e^{−βη}dη (β→0 recovers the existing polynomial
`_mom_eta`), because the tight-core Jacobian (ξ²−η²)/(ξ+η)=(ξ−η) removes the 1/r_A
denominator — no new integral *machinery*, one elementary moment (a by-parts recurrence).

### C1 (the GATE) — PASS, decisive (2026-09-21). `debug/prolate_atomcentered_core.py`.
Primitives: M_xi(n,c)=c^{−(n+1)}Γ(n+1,c) (upper-incomplete-gamma, arbitrary c) +
M_eta(q,β) (recurrence, self-checked vs quad to 6e-51). Isolated Li core (ζ=2.6875)
through the prolate moments at c=β=ζR:
- ⟨χ|χ⟩ = π/ζ³ = 0.1618469255 to **rel 1e-50**, R-SPREAD **2e-51** across R∈[2.70,3.45];
- ⟨χ|1/r_A|χ⟩/N = ζ = 2.68750000 to 1e-50; E_1s = −ζ²/2 = −3.61133 exactly, R-flat.
**The 0.36 Ha wall becomes ~1e-50 → the core-representation problem that sank increment
7 is SOLVED. C1 = GO.** Next: C2 (mixed core-valence ERIs, analytic Neumann) — the
harder piece: extend the V_ee η-side (Ytab) with M_η(q,β) AND the X-table to per-index
radial exponents (survey's secondary obstacle; currently single-α / two-block stitch).

### C2 diagnostic — grid ERIs are R-INACCURATE for the core; the analytic path IS required (2026-09-21). `debug/eri_core_grid_diagnostic.py`, data `debug/data/eri_core_grid_diagnostic.log`.
Tested the easy-hybrid hope (reuse grid ERIs for core classes). It fails. At the FCI's own
grid resolution (N_grid=44, ξ_max=15), with analytic-exact normalization (norm=1 to 1e-15,
so NOT a density/norm issue — purely the 1/r₁₂ kernel on a sharp density; the elliptic-K
singularity makes the sharp core HARDER to integrate than the one-body 1/r, whose focus
singularity cancelled the Jacobian → C1's 1e-50):
- (cc|cc): R-spread **26 mHa**, abs **+132 mHa** (exact 5ζc/8=1.6797) — ~30% of D_e, 16× chem-acc.
- (cc|vv),(cv|cv): R-spread ~4.2 mHa each, abs +21.5 mHa (one core density).
- (vv|vv) control: 2.86 mHa spread, +15 mHa (fastest-converging; may stay on grid).
R-inaccuracy scales with tight-core density content. Grid refinement is a slow wall (cc|cc:
26→21→5.3→1.8 mHa spread at N=44→48→96→160; 13× the FCI points still leaves 1.8/+9.3 mHa) —
same signature as increment-7's one-body core. **VERDICT: C2 must build the analytic
mixed core-involving ERIs. Lean: Legendre-η expansion of the core (e^{−αc η}=Σ(2l+1)(−1)^l
i_l(αc)P_l(η)) turns cores into η-polynomial ProductFns at exponent αc → dissolves the
η-exponential (reuses the validated V_ee η-side), leaving only the mixed-ξ-exponent extension
(the two-block combined-rate stitch, generalized to the full (pq|rs) tensor). (cc|cc) special-
cased to its closed form 5ζc/8. Deep Neumann-code map dispatched to a subagent.**

### C2 build plan (from the Neumann-code map, 2026-09-21). Route A confirmed.
The analytic V_ee (mpf ref `prolate_recondition.vee_mp`/`_build_Xtab_mp`; float
`neumann_vee_general_m.build_Xtab`) = Neumann sum over l of a SEPARABLE
X_l(ξ)·Y1(η)·Y2(η)·pref with the (ξ²−η²) Jacobian.
- **η-side is already GENERAL** (`_mom_eta` takes any η-polynomial) → the Legendre-η
  core (L≈24, de-risked to 1e-15) feeds it UNCHANGED. σ-only ⇒ pure Legendre moments
  ∫η^Q P_l(η) = `_cl_m0` (`debug/prolate_r12_mpf.py:67`). NO η surgery.
- **Single-exponent hard-wire is localized**: `_build_Xtab_mp:334` (c=2α, two_c=2c) +
  the P1↔P2 symmetrization (build_Xtab:442). All blocks (`_mono_moments`,`_B_table`,
  `_seed_B_closed`,`_L_moments`,`_corr`) already take arbitrary rate.
- **The ONE new primitive: `build_Xtab_pair(...,α1,α2)`** — X_l with per-electron rates
  c1=2α1, c2=2α2 (correction B-table at c1+c2, NO P1↔P2 symmetrization). Unavoidable
  either route: (cc|vv) is c1=2αc≠c2=2αv, which the two-block stitch CANNOT synthesize
  (it is single-c called at three values). Route B (η-exponential Y-table) would need
  this SAME extension PLUS new η code → Route A strictly less.
- Assembly (STEP 3): (pq|rs) drives build_Xtab_pair at α1=(αp+αq)/2, α2=(αr+αs)/2;
  only 3 rate-pairs for a 2-exp set {2αc,αc+αv,2αv}². (cc|cc)→closed form 5ζc/8. h1 +
  Löwdin reuse `assemble_hetero_2block` one-body; feed `fci_energy(h1,eri,M,4)`.
- **Validation gates (airtight)**: G1 reduction falsifier α1=α2 → bit-for-bit vs
  `_build_Xtab_mp` (<1e-25); G2 (cc|cc)=5ζc/8=1.67969; G3 (vv|vv)=5ζv/8; G4 (cc|vv)/
  (cv|cv) vs the closed-form 1s Hartree-potential integral; G5 R-independence (grid
  swung (cc|cc) 26 mHa → analytic target ~1e-30). build_Xtab_pair + σ-only assembler
  dispatched to a subagent (`debug/prolate_mixed_eri.py`), gates as hard acceptance.

### C2 CRUX DONE + PM-VERIFIED (2026-09-21). `debug/prolate_mixed_eri.py` (mp.dps=60, ~11s).
`build_Xtab_pair(...,α1,α2)` (mixed-ξ X-table, c1≠c2, correction B-table at c1+c2, no
P1↔P2 symmetry) + σ-only assembler `eri_sigma`. Route A (Legendre-η core, L=24). ALL 5
gates PASS, re-run by the PM (not just the subagent's word):
- G1 reduction falsifier: build_Xtab_pair|_{c1=c2} vs `_build_Xtab_mp` **0.00e+00** over
  825 entries (incl. s>0,m>0) — bit-for-bit.
- G2 (cc|cc)=5ζc/8=1.679687500 (Neumann path 8e-26); G3 (vv|vv)=5ζv/8 to 1e-46.
- G4 (cc|vv),(cv|cv) vs an INDEPENDENT radial-quad reference (no shared code) to ~3e-30;
  (cc|vv) genuinely exercises c1≠c2 (8.1 vs 2.4).
- G5 R-independence: R-SPREAD **0.00e+00** (relerr 5e-32→6e-28 growing with R as the
  forward Q_l recurrence sheds digits — fine at dps=60). **The 26 mHa (cc|cc) / ~4.2 mHa
  cross grid wall is annihilated to ~1e-28 — the two-body analog of C1's 1e-50.**
Convention: unit-normalized orbitals (correct for FCI; the (2/r) 2e-core sibling of
(cc|vv) would be 2×, printed for visibility). **C2 = GO.** Remaining to the R_eq answer:
C3 (wire C1 one-body + C2 σ ERIs → σ-only all-electron FCI = decisive intermediate vs
increment-7's σ-only inward collapse, now with R-accurate integrals) and then the general-m
π/δ core-valence ERIs (build_Xtab_pair's general path is G1-validated at c1=c2; needs the
c1≠c2 exercise for π) + full FCI = C4. σ-only won't be the final geometry (π is the lever,
per HeH⁺ −4.5%→−0.5% at l=2→3) but isolates "did beating the wall fix the collapse?".

### C3 DONE — the σ-only all-electron LiH BINDS; the integral wall WAS increment-7's collapse (2026-09-21). `debug/prolate_allelectron_analytic_fci.py`, data `debug/data/lih_analytic_sigma_req.log`.
Wired C1 one-body + C2 σ ERIs → increment-7's integral-agnostic `fci_energy(h1,eri,M,4)`.
The ONE new piece = a UNIFIED analytic σ one-body engine `one_body_sigma` over the SAME
Route-A `Orbital` objects the ERIs use (core = Legendre-η STO, valence = prolate ProductFn
/ centre-B STO `sto_orbital_B`), so one-body + ERI + FCI share the identical basis. S, V_ne
(heteronuclear −Z_A⟨1/r_A⟩−Z_B⟨1/r_B⟩), T all closed-form in the C1 ξ-moments A_k(c) + the
elementary η-moments. V_NN=3/R, nelec=4, Löwdin S^{−1/2}.
- **CONTROLS ALL PASS.** One-body vs `prolate_recondition` single-e `_ov/_kin/_vne`: S
  bit-exact, H 1.9e-16 (C-A wiring foundation). C-B isolated core (Route-A): norm=1,
  T=ζc²/2=3.6113, ⟨1/r_A⟩=ζc=2.6875, E(Z_A=ζc,0)=−3.6113, **R-spread 0.00e+00** — the
  increment-7 0.36-Ha core wall is gone. C-A H2 (2 valence e, no core, Z=1,1, R=1.40):
  variational, converges 64.9%→**91.9% D_e** (M=2→6), grid-FCI cross-check −1.067 (the
  documented single-exp plateau). C-C variational (all E_tot>−8.070) TRUE. C-D dissociation
  E(8.0)=−7.899 finite, above R_eq region.
- **VERDICT: BINDS, and R_eq IMPROVES with basis** (log `lih_analytic_sigma_req.log`):

  | basis | R_eq | drift vs 3.015 | E_min |
  |---|---|---|---|
  | minimal (M=5) | 3.141 | **+4.2%** | −7.9873 |
  | extended (M=8) | 3.071 | **+1.9%** | −8.0044 (@R=3.015) |

  A clear INTERIOR minimum (E decreases 2.40→~R_eq then rises to E(8.0)=−7.899, D_e(basis)
  ≈88 mHa vs exp ~92 mHa); NOT increment-7's monotone inward collapse.
- **THE KEY READOUT (plainly): beating the integral wall CHANGED the behaviour.** The
  R-accurate analytic core (C1/C2) removes increment-7's grid core artifact (the 0.36-Ha
  R-swing), and the σ-only all-electron LiH then BINDS instead of collapsing inward — so
  increment-7's inward collapse was the GRID CORE ARTIFACT, not the σ-only truncation.
  This is the clean confirmation of the increment-6/7 "variational core is the cure"
  hypothesis that the grid engine could not demonstrate.
- **Trend beats the frozen core AND composed.** R_eq +4.2%→+1.9% IMPROVES with basis — the
  OPPOSITE of frozen-core (+5.5% and WORSENING with l) and composed's structural wall.
  Extended σ-only (+1.9%) already beats composed l-dep-PK (5.3%, Paper 17) and balanced
  (8.8%, Paper 19), and this is BEFORE the π angular lever (HeH⁺: −4.5%→−0.5% at l=2→3),
  which should tighten it further. Energy also basis-converging (−7.987→−8.004 at R=3.015,
  toward −8.070).
- **Scope / honesty.** σ-only PoC (M≤8, single ProductFn/STO exponents; not spectroscopic —
  ~66 mHa above exact at M=8). Cost ~50s/pt (M=5), ~320s/pt (M=8) at dps=60; the many
  distinct orbital exponents make per-(c1,c2) X-table memoization ineffective here. Next =
  C4 (general-m π/δ core-valence ERIs via `build_Xtab_pair`'s c1≠c2 general path + full FCI)
  → the converged prolate all-electron LiH R_eq for a paper claim.
- **PM-VERIFIED (2026-09-22):** re-ran `validate` + `h2` + `lih` independently; every number
  above reproduced bit-for-bit (deterministic mpf). Controls: one-body S=0/H=1.9e-16, core
  R-spread 0.00e+00, H2 wiring 91.9% D_e variational. LiH raw minima: minimal at R=3.20
  (−7.98733), extended at R=3.015 (−8.00435) → fits +4.2%/+1.9% confirmed; BINDS, variational,
  dissociates (E(8.0)=−7.899). Result accepted as established.

### C4 DONE — general-m (π) ERIs built + validated; the CONVERGED prolate LiH R_eq = 3.02 bohr (+0.2%), π-converged (2026-09-22). `debug/prolate_allelectron_c4.py`, data `debug/data/lih_analytic_c4.log`.
The σ-only C2 assembler is generalized to μ≠0. **THE ONE NEW ANALYTIC PRIMITIVE:
`build_Xtab_s(m, s1, s2, ...)`** — generalizes C2's `build_Xtab_pair` from a single
(m,s) shared by both electrons to INDEPENDENT weights s1 (electron 1) / s2 (electron 2),
which is unavoidable for ERIs (e.g. the Coulomb (σσ|ππ) has s1=0, s2=1). At s1=s2 it
reduces to `build_Xtab_pair`→`_build_Xtab_mp` (hence G-REDUCE is bit-for-bit).
`eri_general` uses definite signed-m orbitals (like the grid `vee_m`): selection
m_p−m_q=m_s−m_r, Neumann order m=m_p−m_q, d^|m|P_l on both η-sides (weights s1,s2),
ξ-side from `build_Xtab_s` at the pair rates. NO cos-basis mult=2 (a definite-m orbital
gives a single Neumann term with the (2π)² φ factor = `eri_sigma`'s prefactor).
`one_body_general` adds the μ>0 gradient + μ² azimuthal kinetic (block-diagonal in
(μ, signed m)). Cores stay μ=0 Route-A Legendre STOs; π valence = bond-centred μ=1
ProductFns at m=±1. FCI reuses increment-7's `fci_energy(h1,eri,M,4)` unchanged.
- **ALL GATES PASS (`debug/data/lih_analytic_c4.log`).**
  - **one-body**: μ=1 monomial vs `pr._ov/_vne/_kin` max rel S=0/H=1.1e-15; μ=0 reduces
    to `one_body_sigma` bit-for-bit (0.0e+00).
  - **G-REDUCE**: `eri_general|_{μ=0}` == `eri_sigma` **0.0e+00** over 7 integral types
    incl. the s1≠s2 (cc|vv) and mixed bond/core/valence.
  - **G-PI-REF (independent)**: the grid Cohl–Tohline toroidal kernel
    (`prolate_allelectron_fci._azimuthal_kernels`/`vee_m`) CONVERGES to the analytic
    values as N_grid→∞ for μ=0,1,2 (μ=0 8.300→8.232→8.206→[Rich 8.190] vs analytic
    8.179; μ=1 3.820→3.676→3.622→[3.589] vs 3.566; μ=2 8.44→7.79→7.55→[7.40] vs 7.295,
    grid-limited δ). Same convergence rate across μ ⇒ the analytic values are correct
    and the gap is pure grid discretization.
  - **G-RIND**: atom-centred core–π ERIs (cc|π_A π_A), (π_A c|c π_A) R-spread **0.00e+00**
    (analytic core R-exact even through the μ=1 machinery). NB the FCI's own π are
    BOND-centred and correctly R-dependent, like the σ bond functions.
  - **G-H2PI (make-or-break)**: H2 σ-only 91.91% → σ+1π 97.10% → σ+2π 97.95% D_e (→98.5%
    at 3–4 π shells). π climbs decisively past the σ-only wall toward the known ~99% —
    the π assembly is validated end-to-end through the full FCI.
  - **G-VAR**: all LiH E_tot > −8.070 (variational throughout).
- **VERDICT — π (the angular lever) TIGHTENS LiH R_eq to experiment, and it CONVERGES:**

  | basis (extended σ + π) | R_eq | drift vs 3.015 | E_min |
  |---|---|---|---|
  | C3 σ-only | 3.071 | **+1.9%** | −8.0044 |
  | +1 π shell (m=±1) | 3.020 | **+0.2%** | −8.00794 |
  | +2 π shells | 3.021 | **+0.2%** | −8.01071 |

  Clean interior minimum at R≈3.015 (E rises both sides, dissociates: E(4.0)=−7.990,
  E(6.0)=−7.946). **R_eq is CONVERGED w.r.t. the π basis (1π=2π=+0.2%)** — the +1.9%→+0.2%
  jump is the π polarization, exactly the HeH⁺ precedent (−4.5%→−0.5% at l=2→3). This is
  the increment-6/7 "variational core is the cure" hypothesis DEMONSTRATED CLEANLY: the
  R-accurate analytic core (C1/C2) removes the grid artifact, and the angular (π) lever
  then converges the geometry to experiment.
- **Beats every prior LiH R_eq DECISIVELY:** composed l-dep-PK 5.3% (Paper 17), balanced
  8.8% (Paper 19), frozen-core prolate +5.5% (incr 6), C3 σ-only +1.9% — and even HeH⁺'s
  own −0.5%. The +5.5%-OUTWARD-and-worsening frozen-core wall is gone.
- **HONEST SCOPE — the GEOMETRY is paper-grade; the ENERGY is not (yet).** R_eq=+0.2%
  (π-converged) is a genuine near-spectroscopic geometry and is the deliverable Route C set
  out to get. But E_min is still ~60 mHa above exact (−8.011 vs −8.070) because each orbital
  carries a SINGLE exponent (no multi-exponent radial ladder / completeness) — so the energy
  and D_e remain basis-limited (D_e unreliable at large R; the separated-atom limit is poorly
  spanned by bond-centred functions). Geometry is the meaningful readout (as throughout this
  track); the energy needs a multi-exponent radial ladder (the prolate_recondition route) to
  reach chemical accuracy — a scoped future build, orthogonal to the now-solved angular lever.
- Cost ~340s/pt (M=10, dps=60) / ~375s/pt (M=12); X-tables memoized by (m,s1,s2,c1,c2).
  Deliverables: `debug/prolate_allelectron_c4.py` (validated, reusable), the log, this note.
- **PM-VERIFIED (2026-09-22):** ran the make-or-break **G-H2PI** myself (σ 91.91% → σ+1π
  97.10% → σ+2π 97.95% — bit-for-bit) AND independently re-ran the LiH σ+1π bond-range
  points (`debug/data/c4_pm_verify.log`): R=2.850 −8.00695, 3.015 **−8.00794**, 3.200
  −8.00683 → interior minimum at 3.015 confirmed (matches the agent log exactly,
  deterministic). **C4 result accepted as established: LiH R_eq = +0.2% (π-converged),
  geometry paper-grade; energy ~60 mHa above exact (single-exponent radial PoC).**

## ROUTE C — COMPLETE (2026-09-22). The from-scratch prolate all-electron LiH GEOMETRY is at experiment (+0.2%, π-converged); best in the corpus. Chain closed: breathing NO → polarization NO → the frozen-core APPROXIMATION is the culprit → R-accurate analytic core (C1 1e-50 / C2 1e-28) → binds (C3, +1.9%) → π converges to experiment (C4, +0.2%). Energy-to-chemical-accuracy (multi-exponent radial ladder) is the scoped orthogonal follow-on. Owed at sprint-close: CHANGELOG entry; Paper 19 sharpen (with tests/ backing, cite permanent record not debug/); MEMORY index.

## ROUTE C — ENERGY follow-on (RADIAL LADDER) — done, det-wall-bounded (2026-09-22). `debug/prolate_energy_ladder.py`, data `debug/data/prolate_energy_ladder.log`.
The C4 geometry is solved but E was ~60 mHa high (single-exponent radial incompleteness). Added a radial+azimuthal ladder ON TOP of the C4 engine (no new integral code — C4's `build_Xtab_s` already takes mixed exponents/weights; a ladder is just more OrbitalM's). **The one new enabling piece: a float64 Laguerre(ξ)×Legendre(η) RE-BASING** (block-diagonal per channel, from `prolate_recondition.laguerre_coeffs/legendre_coeffs`) — the same conditioning fix Paper 12/prolate_recondition use. Validated: reproduces the monomial energy bit-for-bit while cutting cond(S) 2.9e9→1.4e5.
- **G-REDUCE PASS** (size-1 ladder == C4 assemble_energy_m, 2.2e-15).
- **The conditioning wall is REAL and the re-basing BREAKS it (measured):** monomial float64 δ-radial J=4 REGRESSES to 98.93% (66/75 kept); re-based → 99.611% (75/75, monotone).
- **H2 (make-or-break control):** ladder climbs σ 92.41 → +π 99.22 → +δ 99.60 → +φ 99.63 → deep 99.664% (M=93, 0.59 mHa, all kept, monotone, variational, α≈1.2 optimal). **Reaches CHEMICAL ACCURACY**, ~3× better than C4's σ+π (97.1%). Saturation map: radial saturates fast, σ-only caps ~92.5%, η-angular (l) saturates ~L=2, the **azimuthal μ ladder is the lever** (diminishing returns). Did NOT hit the literal 99.8% gate: the last ~0.13pp is the partial-wave e-e cusp — prolate_recondition reaches 99.77–99.81% only via the far-more-compact **2-electron product CI** at (5,5,2) (identical span, but orbital-FCI's ndet=M² blows up if l→5 on all channels), and the cusp is fully closed only by explicit r₁₂ (the R12-CI arc, H2 99.97%). So the radial-ladder APPROACH/WIRING is confirmed sound; 99.8% is a representation-compactness/cusp limit, not a ladder failure.
- **LiH energy @ R=3.015 (the actual target):** radial ladder + π lowers E MONOTONICALLY & variationally −7.995 → −8.003 → −8.009 → **−8.012** (M=6→16), then **PLATEAUS ~58 mHa above exact (−8.070) at the 4-electron FCI determinant wall** (ndet=C(M,2)², caps M~16). Conditioning is NOT the LiH limit (cond ≤3.6e6, all kept) — the determinant count is. Residual is dominated by core-core correlation (the analytic Li 1s core is a SINGLE orbital → FCI captures zero core correlation, ~40 mHa) + valence completeness, both needing more orbitals than dense FCI ndet allows.
- **VERDICT:** the radial ladder is the right cure for the ENERGY and is validated to chemical accuracy on H2; for LiH the binding constraint is the 4e **determinant wall**, not the ladder or conditioning. Closing the full LiH gap needs a determinant-frugal solver (CASCI/DMRG/selected-CI over the ladder) OR explicit r₁₂ core-valence correlation (the R12-CI arc) — NOT more radial functions in dense FCI. C4 geometry (+0.2%) stands; G-GEOM re-fit is gated on first reaching the converged energy, which the det wall prevents here. Deliverables: `debug/prolate_energy_ladder.py` (reusable: builders + the re-basing solve + gates), the log.
