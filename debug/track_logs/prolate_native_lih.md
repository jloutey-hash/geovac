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

## Next increments (remaining) — reordered after the l=3 divergence
5. **DONE: angular l=2→3 (p=0) — DIVERGED (+7.1%), confounded by missing r₁₂.** Not
   the clean HeH⁺ convergence; see the QUALIFIED note above.
6. **CRITICAL NEXT: r₁₂-coupled V_H + projector.** Extend the core-shield V_H and the
   Huzinaga projector from p=0-factorized to the p={0,1} (r₁₂-coupled) matrix
   elements, so LiH can be run WITH r₁₂ like the HeH⁺ PoC that converged. This is the
   apples-to-apples test and the load-bearing open question: does prolate LiH converge
   (r₁₂ was the lever) or drift (frozen core is the culprit)? Real work, not free.
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
