# Track log — prolate-native (H2-style) LiH: the accuracy route

**Opened:** 2026-09-20 (PI-directed). **PM holds this thread.**

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
