# LiH R12 F12-synthesis — ceiling diagnostic (design pass, before the big build)

**Date:** 2026-09-22 (v5.15.17 HEAD). **PI-directed:** diagnostic-before-engineering,
"pin the achievable ceiling before the big build — don't launch it blind."
**Probe:** `debug/lih_r12_ceiling_probe.py` (float64, ~seconds, self-validating).
**Verdict:** the literal synthesis (geminal on Route C's −8.012) has a ceiling of
**~−8.042** (≈28 mHa above exact, a real ~2× gain, NOT chemical). Chemical accuracy
needs a **2nd core orbital** (radial half) FIRST — a cheap existing-code change — with
the geminal supplying only the cusp half. **Recommendation: do the cheap 2nd-core-exponent
experiment before committing to the multi-session two-center 4e mpf geminal rebuild.**

---

## STEP-1 EXPERIMENT — DONE, prediction CONFIRMED (2026-09-22, PI-approved)

`debug/lih_core2exp_probe.py` (adds a 2nd/3rd core STO at a different exponent to Route
C's ladder; no engine change — the `(cc|cc)=5ζ/8` shortcut correctly ignores
cross-exponent core pairs). All at R=3.015, dps=60.

**A 2nd core "shape" recovers the radial half the geminal cannot** (minimal valence
bond(1,0), baseline 1-core E=−7.99468):

| core partners | ΔE | note |
|:--|:--:|:--|
| +ζ₂=4.0 (tight) | −16.9 mHa | one partner ≈ half the radial half |
| +ζ₂=5.5 (tight) | −17.0 mHa | tight is the lever (in–out dodge) |
| +ζ₂=1.6 (loose) | −8.8 mHa | loose alone weaker |
| +[4.5, 1.6] | −23.5 mHa | two partners → toward the ~30 mHa radial half |
| +[4.5, 8.0] | −22.2 mHa | diminishing returns |

**Two decisive readouts:**
1. **The radial half is real and reachable by orbitals** — ~17 mHa (1 partner) → ~23 mHa
   (2 partners) → trending to the diagnostic's ~30 mHa. Exactly the half the geminal
   plateaus below.
2. **Core orbitals beat valence orbitals per unit of determinant budget** — minimal
   valence + 2 cores (M=8, E=−8.018) **already beats full Route C** (M=16, E=−8.012).
   So the determinant-wall worry is resolved *favorably*: the radial half costs only 2–3
   orbitals and is more efficient than the valence Route C was spending its budget on.
   The wall is NOT the blocker for this half.

**Moderate-valence confirm** (bond(2,1), no π): 1-core E=−8.00329 → +[4.5,1.6]
E=−8.02335, **ΔE=−20.1 mHa**. The core-enrichment gain is robust across valence contexts
(minimal −23.5, moderate −20.1 mHa). (This no-π config is not comparable in absolute terms
to Route C's −8.012, which carries π; the transferable quantity is ΔE.)

**Determinant budget — RESOLVED.** Full Route-C valence + 2–3 core exponents would be
M~18–19 (past the dense-FCI ceiling ~14400), but that is NOT required: (a) core orbitals
beat valence per unit budget (M=8 core-enriched −8.018 > M=16 valence-heavy −8.012), so a
core-enriched allocation reaches the radial half within M≤16; and (b) the geminal (cusp
half) is a low-rank {Φ0, FΦ0} correction, not more determinants — it does not consume the
budget. So the dense wall does not block the combined route.

## STEP-1 BANKED NUMBER (2026-09-22) — core enrichment at full-π valence, M≤16

`debug/lih_core2exp_probe.py bank`, data `debug/data/lih_core2exp_bank.log`. Valence =
bond(2,1)+1π(1,0) both runs; only core count differs.

| config | M | nd | E | from exact |
|:--|:--:|:--:|:--:|:--:|
| A) 1 core (matched baseline) | 14 | 8281 | −8.00898 | −61.0 mHa |
| **B) 3 cores [4.5, 1.6]** | 16 | 14400 | **−8.02905** | **−41.0 mHa** |

**ΔE(core enrichment) = −20.1 mHa** at matched valence. **B beats Route C's −8.012 by
17 mHa at the SAME M=16 budget** — pure orbital reallocation (2 valence → core), existing
code. Banked best pure-orbital LiH energy = **−8.029 (41 mHa from exact)**.

**Honest read on the determinant wall.** The FIRST ~20 mHa of the radial half is
budget-cheap (demonstrated). Reaching the FULL radial half (~30 mHa, isolated-core probe)
is wall-limited: within dense M=16 you can't hold 3–4 core orbitals AND full valence, so
~10 mHa of radial + valence completeness stays out of reach until M>16 (a float64
productionization to break the dense ceiling, or a determinant-frugal solver).

**Updated quantitative marriage (the roadmap):**
- Route C (1 core, M=16): −8.012 (58 mHa short)
- + core enrichment (dense M=16, DEMONSTRATED, existing code): **−8.029** (41 mHa short)
- + geminal (cusp half, ~30 mHa, the F12 build, low-rank on top): → **~−8.059 (near-chemical, ~11 mHa)**
- + break M=16 (float64 / frugal solver) for the last radial+valence: → **chemical (<2 mHa)**

Both halves needed; neither alone reaches chemical accuracy. The cheap half is now a
concrete banked number. **The F12 geminal build is justified for its correct role — the
~30 mHa cusp remainder on the −8.029 core-enriched base — and lands the total near-chemical
(~−8.06); the final <2 mHa additionally needs breaking the dense M=16 ceiling.**

**Next-step options (PI call):** (2) the F12 geminal build on the −8.029 core-enriched base
(cusp remainder → ~−8.06); and/or (1b) the float64 productionization of the Route-C engine
(the handoff's minutes→seconds enabler), which both lets core enrichment break M=16 for the
last radial/valence mHa AND makes the whole thing a sweep instead of a few slow points.

## The question
Route C is geometry-complete (R_eq +0.2%) but energy-walled at **E=−8.012** (58 mHa above
exact −8.070) by the 4e FCI **determinant wall**: the analytic Li 1s core is a SINGLE
orbital, so dense FCI captures ZERO core-core correlation (~40–57 mHa), and ndet=C(M,2)²
caps M~16. The PI-directed path is F12 synthesis — an r₁₂ geminal correction on Route C's
good base. Before the substantial mpf rebuild (the r₁₂ integrals are hard-wired to the
ionic 2×2 reference; re-deriving them on Route C's multi-orbital base + F12 double-counting
= multi-session), pin how much a geminal can actually buy.

## The probe
The dominant residual is the Li 1s² core-core pair. Isolated Li⁺ (He-like Z=3) is the clean,
transferable proxy (the core is barely perturbed by the distant H⁻; core correlation is
local). Because Route C's core holds ZERO correlation, the geminal's core contribution is
**double-counting-free** — a rare clean case. 2×2..(1+n) generalized eigenproblem
{Φ0, b_k Φ0}, Φ0=1s(ζ)² at **ζ=2.6875 fixed** (= Route C's frozen analytic core exponent,
NOT re-optimized — the realistic scenario), Hylleraas coordinates, gradient-form kinetic
energy. Validated: norm=1, ⟨1/r₁⟩=ζ, T00=ζ² exact; He {u,t²}=−2.892 reproduces the classic
3-term Hylleraas −2.9024; 5-term → 98% of He's 57 mHa. Numbers trustworthy.

## The finding — the core correlation is two roughly-equal, near-orthogonal halves
Li⁺ core deficit vs the single-ζ reference = **60.1 mHa** (≈13.7 mHa orbital-shape +
~43 mHa true correlation; the isolated-core analog of Route C's 58 mHa LiH gap).

| lever | basis | He rec | Li⁺ rec | reachable by |
|:--|:--|:--:|:--:|:--|
| r₁₂ cusp | {u} | 48% | 48% | geminal |
| r₁₂ cusp | {u, u²} | 52% | 52% | geminal — **PLATEAUS** (more u-terms, cond↑, no gain) |
| radial in–out | {t²} | 51% | 45% | a 2nd core **orbital** |
| radial in–out | {t², s} | 51% | 45% | orbital — plateaus |
| **both** | {u, t²} | 79% | 74% | geminal + 1 radial orbital (near-additive) |
| rich | {u,t²,s,u²,ut²} | 98% | 95% | small explicitly-correlated basis |

**Load-bearing:** a geminal ALONE — any number of r₁₂ terms — caps at **~50%** of the
core-core deficit (~30 of ~60 mHa), because it rides a *fixed single-ζ* radial core. The
other half is a **radial orbital degree of freedom** (t²=(r₁−r₂)² = one electron in, one
out) the geminal cannot express. This is exactly why **H2 reached 99.97%** (radial ladder
AND r₁₂) while **Be R12-CI reached 19%** (crude single geminal, neither) — the two corpus
anchors that bracketed the ceiling now have a mechanism.

## Ceiling of the literal synthesis
Route C −8.012 + geminal core (cusp half, ~30 mHa) + a few mHa short-range valence (partly
double-counted) ⟹ **~−8.042 to −8.048**, i.e. **≈22–28 mHa above exact**. Real (58→~25 mHa,
~2×) but **not chemical accuracy** — and less than the handoff's optimistic −8.06, because
that assumed the geminal reaches 30–50 mHa; it plateaus at ~30.

## The reframe (the actionable output)
The residual is dominated by the frozen single-orbital core, and it splits cusp/radial ~50/50:
- **radial half (~30 mHa)** = a 2nd tight core exponent (1s+1s′). Route C's C4 engine
  (`build_Xtab_s`) ALREADY takes mixed exponents/weights — **no new integral code**. One
  extra orbital, M~16→~17. Directly tests the track log's stated blocker (does the
  determinant wall actually block a 2nd *core* orbital, or was the ladder simply all-valence?).
- **cusp half (~30 mHa)** = the geminal. This is the piece that genuinely needs the
  expensive two-center 4e RI-free machinery — and it's the half where the corpus's RI-free
  reduction is "ahead of F12."

Both core pieces are double-counting-free (core corr = 0 in Route C). Double-counting risk
is confined to the valence.

## Recommendation (GO, re-scoped — cheap experiment BEFORE the big build)
1. **FIRST (cheap, ~1 session, existing code):** add a 2nd (and 3rd) tight core exponent to
   Route C's radial ladder; measure the core correlation recovered and whether the
   determinant wall admits it. Predicted ~+30 mHa → E≈−8.042 with NO geminal.
   - If it recovers the radial half → the geminal is then needed only for the final cusp
     half (~30 mHa), and {radial orbital + geminal} → ~95% → **~−8.068, chemical**. The
     expensive geminal build is justified *for the cusp remainder*, and you know exactly
     how much it must buy.
   - If the determinant wall blocks the 2nd core orbital → dense FCI is exhausted AND the
     geminal cannot reach the radial half either (it plateaus at cusp) → chemical accuracy
     needs a **determinant-frugal solver** (CASCI/DMRG/selected-CI over the ladder, the
     track log's own named alternative), and the geminal only ever buys ~30 mHa regardless.
     This reprioritizes toward the frugal solver.
2. **THEN (only if step 1 banks the radial half):** the two-center 4e mpf geminal build for
   the cusp remainder. The `geovac/lih_r12ci` machinery (4-body bridge + 3-body triangle,
   RI-free) is the substrate; re-derive on the enriched base + handle valence double-counting.

Either way, step 1 is the correct next move: it's cheap, reuses validated code, resolves the
determinant-wall question the whole strategy hinges on, and prevents launching a multi-session
mpf build that (on the current single-core-orbital base) tops out at ~half the residual.

## Honest scope / caveats
- Isolated-core proxy: the ~50/50 split and ~60 mHa magnitude are for the isolated Li⁺ core;
  embedding in LiH shifts by a few mHa (core correlation is local/transferable — safe).
- The valence (H⁻-side) short-range correlation is a separate, smaller piece not probed here;
  Route C's ladder already captures much of it.
- ζ fixed at 2.6875 is deliberate (it IS Route C's core); the split is a property of the
  fixed-core scenario, which is the one that governs the synthesis.
- Deliverables: `debug/lih_r12_ceiling_probe.py` (validated, reusable). No production code
  touched; no paper claim yet (this is a design pass).
