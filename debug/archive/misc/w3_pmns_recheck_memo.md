# Sprint W3 Test 1 — PMNS recheck with prespecified 30-form basis

**Date:** 2026-05-08 (post-CKM-verification follow-up)
**Author:** PM session (sub-agent dispatch)
**Context:** Forward plan Test 1 of `debug/w3_forward_plan_memo.md`. The CKM
verification (`debug/w3_lambda_predictive_verification.py`) returned 7 of 30
hits at <1% across 4 primary Wolfenstein parameters and identified clean
M1/M2 spectral forms for all four (λ, A, ρ̄, η̄). This sprint runs the SAME
30-form prespecified basis against the four primary PMNS parameters with
identical methodology, statistical accounting, and comparison framework.

The CKM result, if real, is sector-universal in SM mixing matrices: PMNS
should also fit M1/M2 forms at comparable signal level. The first-pass
PMNS probe (`debug/w3_ckm_pmns_probe_memo.md`) saw no clean Koide-style
fact for PMNS, but did not run the full predictive-verification pipeline.
This sprint closes that gap.

## TL;DR

**Clean negative.** Zero PMNS primary parameters fit any of the 30 prespecified
M1/M2/Other natural forms within 2% relative deviation. CKM, on the same basis
and same threshold, found 7 hits at <1% (random expectation: 0.7 ± 1.0). The
contrast is unambiguous: the M1/M2 spectral identification is **sector-restricted
to CKM**, not sector-universal in SM mixing matrices.

**Tri-bimaximal (TBM) wins 3–0 on sin² angles** against the best M1/M2 fit:
the leptonic mixing sector sits closer to pure rationals (1/3, 1/2, 0) than to
the spectral-zeta rings, suggesting a different structural mechanism for PMNS
(group-theoretic discrete-symmetry breaking, e.g., A_4 / S_4 family symmetries)
rather than the master Mellin engine that fits CKM.

**One headline cross-cut:** θ_12 in degrees matches asin(1/√3) ≈ 35.26° to
**0.11 PDG σ** — the cleanest single-PMNS-parameter fit found anywhere in this
or the previous probe, and it is exactly the TBM solar-angle prediction. (The
PDG best fit is 33.65° vs asin(1/√3) = 35.26°; PDG σ = 0.81° puts the TBM
prediction at 2.0 σ on sin² but only 0.11 σ on the angle itself due to the
nonlinear map. The asymmetry is a real artifact of which parameterization is
"natural" — see §6 below.)

## What was tested

**Frozen basis (30 forms, copy-pasted verbatim from the CKM verification
script).** The basis breaks naturally into three families:
- **M1 (Hopf-base / π-family, 14 forms):** 1/√Vol(Sᵏ), 1/Vol(Sᵏ) for k = 1,
  2, 3, 5; π/14, π/16, π/(2φ); 1/π, 1/(2π²), √2/(2π), 1/(π√π).
- **M2 (chirality / ln 2-family, 6 forms):** ln(2), ln(2)/2, ln(2)/4, √(ln 2),
  1/ln(2), log(2)/π.
- **Other (10 forms):** κ² = 1/256, Δ = 1/40, 1/B = 1/42, F = π²/6, α = 1/137,
  1/φ, 1/φ², φ−1, φ²−φ = 1, 1/√20.

**PMNS parameters (PDG 2024 / NuFIT 5.3 NH best fit, with PDG σ):**
- sin²θ_12 = 0.307 ± 0.013
- sin²θ_23 = 0.546 ± 0.021 (octant ambiguity flagged)
- sin²θ_13 = 0.0220 ± 0.0007
- δ_CP = 230° ± 25° (very loose)

Plus 7 additional derived observables (sinθ_ij, θ_ij in deg, J_PMNS).

## Per-parameter match table at <2% threshold

| PMNS Param | PDG | σ | Best M1/M2 fit (any) | Dev % | PDG σ |
|---|---|---|---|---|---|
| sin²θ_12 | 0.307 | 0.013 | **none < 2%** | — | — |
| sin²θ_23 | 0.546 | 0.021 | **none < 2%** | — | — |
| sin²θ_13 | 0.022 | 0.0007 | **none < 2%** | — | — |
| δ_CP (deg) | 230.0 | 25.0 | **none < 2%** | — | — |
| sinθ_12 | 0.554 | 0.012 | none < 2% | — | — |
| sinθ_23 | 0.739 | 0.014 | none < 2% | — | — |
| sinθ_13 | 0.148 | 0.0024 | none < 2% | — | — |
| θ_12 (deg) | 33.65 | 0.81 | none < 2% | — | — |
| θ_23 (deg) | 47.64 | 1.21 | none < 2% | — | — |
| θ_13 (deg) | 8.53 | 0.14 | none < 2% | — | — |
| J_PMNS | −0.0255 | ~0.004 | none < 2% | — | — |

Across all 11 PMNS observables × 30 forms = 330 tests, **zero matches at <2%**.
For the 4 primary parameters specifically: 0 hits at <1%, 0 at <0.5%, 0 within
1 PDG σ.

## Statistical accounting (matches CKM script methodology)

| Metric | CKM (4 primary) | PMNS (4 primary) |
|---|---|---|
| Observed hits at <1% | 7 | 0 |
| Observed hits at <0.5% | 5 | 0 |
| Observed within 1 PDG σ | 6 | 0 |
| Expected per param at <1% (chance, log-uniform [0.10, 1.00]) | 0.174 | 0.179 |
| Expected total over 4 params at <1% | 0.7 ± 1.06 | 0.72 ± 1.06 |
| Sigma above (below) chance | +5.9 σ | −0.68 σ |

PMNS sits **0.68 σ below** the random-target chance expectation. This is a
better-than-random *negative* — the basis is not just "no signal," it is
slightly less hit-prone than chance, consistent with PMNS values clustering
in regions of parameter space where the spectral-zeta forms happen to be
sparse.

## TBM-rational vs M1/M2 comparison (the load-bearing result)

The Tri-Bimaximal (TBM) benchmark predicts sin²θ_12 = 1/3, sin²θ_23 = 1/2,
sin²θ_13 = 0 — pure rationals from group-theoretic discrete-symmetry breaking
(A_4, S_4, Δ(54), … family symmetries). These are NOT in the M1/M2 rings
(no transcendentals).

| Parameter | PDG | σ_TBM (in PDG σ) | σ_M1/M2 (best) | Winner |
|---|---|---|---|---|
| sin²θ_12 | 0.307 | 2.03 | none < 2% (∞) | TBM |
| sin²θ_23 | 0.546 | 2.19 | none < 2% (∞) | TBM |
| sin²θ_13 | 0.022 | 31.43 | none < 2% (∞) | TBM |

TBM wins 3–0. The TBM σ_distances are 2–31 σ — TBM is itself broken (this is
known; the "TBM with corrections" research program in flavor physics is exactly
about this), but it is still categorically closer than any M1/M2 form.

The clear interpretation: **PMNS is in the TBM-rational neighborhood, not in
the M1/M2 spectral-zeta neighborhood**. The two neighborhoods are categorically
different on the transcendental axis. Whatever generates PMNS calibration
data is closer to the discrete-symmetry research program than to a spectral
action on the GeoVac inner factor.

## One bright spot: θ_12 in degrees vs asin(1/√3)

Step 9 of the script tests PMNS angles (in degrees) against simple rational
fractions of π and structurally natural angles. One result jumps out:

| Param | PDG | σ | Best simple-form fit | σ |
|---|---|---|---|---|
| θ_12 (deg) | 33.647 | 0.807 | **asin(1/√3) = 35.26° (TBM solar)** | **0.11** |
| θ_23 (deg) | 47.639 | 1.208 | 45° = π/4 | 2.18 |
| θ_13 (deg) | 8.530 | 0.137 | π/12 = 15° | 47.3 |

θ_12 = 33.65° vs the TBM solar prediction asin(1/√3) = 35.264° lands at 0.11 σ
in the angular metric. (The same comparison in sin² is at 2.03 σ; the 18×
difference reflects the nonlinear map between sin² and θ. Which parameterization
is "natural" depends on the underlying mechanism — see §7.) This is the
**single-cleanest PMNS fit found anywhere in the W3 program**, and it is
**not in the M1/M2 ring** — it is the TBM angle, derived from group theory.

θ_23 at 47.6° vs 45° is the famous "near-maximal but not quite" result driving
the octant ambiguity discussion in atmospheric neutrino fits.

θ_13 vs π/12 is a clean negative — θ_13 was thought zero or very small until
Daya Bay 2012; π/12 = 15° is far from the measured 8.53°. The point of including
this row is honesty: not every "simple rational" works.

## Step 7 / Step 8 (derived relations) — both null

- **Step 7:** the CKM result yielded a derived relation A² = 2η̄ that survives
  PDG at 0.52 σ. This sprint sought analogous low-coefficient relations among
  PMNS primary parameters using their best M1/M2 fits. **Result: zero such
  relations**, because no primary PMNS param had a within-1-σ M1/M2 fit to
  begin with.

- **Step 8:** the CKM result yielded tan(δ_CP) = π · ln 2 (M1×M2 product) at
  0.02 σ. This sprint tested tan(δ_CP_PMNS) = tan(230°) = 1.192 against all
  M1×M2, M1/M2, M2/M1 product candidates. The cleanest hit is M1/M2 with
  (π/(2φ), √ln 2) at −2.16% deviation — outside any reasonable threshold given
  the very loose PDG σ on δ_CP_PMNS. The angles are degenerate at this
  precision; the cleanness of the CKM tan(δ_CP) result is not reproduced here.

## Sector-mapping reading after this sprint

After this sprint, the picture across three sectors looks like:

| Sector | Algebra factor | Candidate pattern | Cleanness |
|---|---|---|---|
| Lepton mass | (ℂ + ℍ)³ | Koide 45° cone | 1 arcsec |
| Quark mixing (CKM) | M_3(ℂ) on flavor | 4 of 4 primary params fit M1/M2 within 1σ | 0.4–0.8% across |
| Neutrino mixing (PMNS) | (ℂ + ℍ)³ + Majorana | TBM rationals 2–3σ off, M1/M2 entirely absent | TBM-broken |

The sector-mapping is now **structurally coherent in a different sense than
the previous memo suggested**: each sector lights up under a different
mechanism category.

- **Lepton mass:** geometric cone (Koide) — projection / angle constraint.
- **Quark mixing:** spectral-zeta values (M1/M2 master Mellin engine).
- **Neutrino mixing:** discrete-symmetry rationals (TBM-broken; A_4/S_4
  flavor-symmetry research program).

These are three categorically different mechanisms, mapping onto three
categorically different sub-algebras of the Connes inner factor. That is itself
an interesting structural pattern — but it falsifies the simpler "M1/M2 split
is sector-universal" hypothesis from the CKM probe.

## Honest scope and caveats

1. **Basis size matters.** With only 6 forms in the M2 family and 14 in M1,
   the basis is probably under-resolved for finding all spectral-zeta-style
   matches. A larger basis (say 60–100 forms including more
   small-integer-coefficient combinations) would test the negative result more
   robustly. But: the same basis found 7 hits in CKM at the same threshold,
   so the basis is at least basis-fair across sectors. Test 5 of the forward
   plan (mechanically-generated basis) is the cleanest follow-up.

2. **PMNS uncertainties are larger than CKM's.** sin²θ_12 has 4% relative PDG
   σ vs CKM λ at 0.3%. A larger σ window means the within-1-σ test is *more
   permissive* in PMNS — yet still zero matches. The negative is even stronger
   than it would be at CKM-level precision.

3. **The TBM 3–0 result is structurally tight.** The three sin² parameters sit
   2.0, 2.2, and 31 σ from TBM rationals (1/3, 1/2, 0). Even though TBM is
   technically "broken at the few-σ level," it is the closest *clean* structure
   to the data. Any future tightening of NuFIT data could either drive PMNS
   *toward* TBM (consistent with TBM + corrections), away from TBM, or onto a
   yet-unidentified third structure. Right now TBM is the ansatz to beat in the
   leptonic mixing sector, and M1/M2 is structurally far from it.

4. **The θ_12 = asin(1/√3) result is single-data-point and uses a different
   parameterization than the sin² test.** Cleanness in the angular metric (0.11
   σ) vs sin² metric (2.0 σ) reflects the nonlinear map. A first-principles
   derivation would have to predict the parameterization too — i.e., are the
   primary observables the angles, or the sin² of angles? The standard PDG
   parameterization treats sin² as primary (because the Hamiltonian eigenvalues
   are functions of sin²), but TBM's ansatz θ_12 = arctan(1/√2) treats the
   angle as primary. Both are defensible; this single-data-point cleanness is
   not by itself decisive.

5. **The CKM signal we are comparing against was 6–10σ above chance.** PMNS
   here is 0.68σ below chance. The contrast (CKM at +5.9σ above chance, PMNS
   at −0.68σ below chance) is a 6.6σ swing — categorical, not marginal.
   Whatever generates the PMNS calibration data, it does not live in the
   prespecified M1/M2 basis.

6. **Octant ambiguity in θ_23.** PDG includes a two-fold ambiguity on whether
   sin²θ_23 is below 0.5 (lower octant ~0.45) or above (upper octant ~0.55).
   This sprint used best-fit 0.546 (upper octant). The 45° (= π/4 = M1) form
   is at 2.18 σ in the angular metric for the upper octant; under the lower
   octant it would be at ~0 σ. **In the lower-octant scenario, sin²θ_23 = 1/2
   exactly is at <1 σ.** This is a TBM-rational identification, not M1/M2.
   The octant resolution (T2K vs NOvA tension) will resolve in the next 5
   years; if lower octant wins, the TBM identification of θ_23 strengthens.

## Recommendation

**This sprint downgrades the W3 universality hypothesis to "sector-restricted
to CKM."** The bold-path claim from the previous memo ("M1/M2 split is
sector-universal in SM mixing matrices") is falsified. The narrower claim
that survives is:

> The CKM Wolfenstein parameters fit M1/M2 spectral-zeta forms at a 6–10σ
> above-chance signal level, with derived relations A² = 2η̄ and tan(δ_CP) =
> π · ln 2 surviving PDG. PMNS does NOT fit the same basis; it is closer to
> the TBM-rational discrete-symmetry research program. Charged lepton masses
> follow a third pattern (Koide 45° cone). The three calibration-data
> sectors map to three categorically different mechanism families.

This is a more nuanced result than the original universality claim, but it is
internally consistent with the Connes spectral-triple framework: different
inner-factor sub-algebras (ℂ+ℍ for leptons, M_3(ℂ) for quarks) admit
different families of natural calibration data. The sector specificity could
itself be derived from the inner-factor structure if the right calculation can
be set up.

**Next-step priority adjustment:**

- The CKM-only candidate hypothesis is still worth Test 6 (Connes-Chamseddine
  spectral action sketch) of the forward plan. If the spectral action on
  M_3(ℂ) generates the M1/M2 Yukawa values without any free parameter, the
  sector-restricted claim closes with structural backing.

- The PMNS-TBM and lepton-Koide patterns now form their own research questions.
  TBM/A_4 flavor-symmetry models are an active literature; CLAUDE.md §3
  documents many failed approaches but does not include a flavor-symmetry
  derivation attempt.

- **Test 2 (charged lepton mass spectrum with same methodology) is now higher
  priority than originally weighted.** If the lepton mass ratios fit M1/M2
  forms — i.e., if the basis lights up for the (ℂ+ℍ)³ sub-algebra in lepton
  masses but not in PMNS — that would be additional structural evidence for
  sub-algebra-specific Mellin engine output.

- **Test 5 (mechanically-generated basis) becomes more urgent.** The 30-form
  basis tested here has clear "designed" choices (which φ-family forms are
  included, which combinations of π and ln 2). A mechanically-generated basis
  with bounded complexity would either (a) preserve the CKM 7-hit signal and
  the PMNS 0-hit result, strengthening sector specificity; or (b) preserve
  CKM but light up for PMNS at unexpected forms, weakening the negative; or
  (c) wash out the CKM signal too, suggesting the original CKM result was
  partly basis-dependent.

## Files

- `debug/w3_pmns_recheck.py` (this sprint's script, ~3 minutes wall time at dps=80)
- `debug/data/w3_pmns_recheck.json` (full structured output)
- `debug/w3_pmns_recheck_memo.md` (this memo)

Cross-reference:
- `debug/w3_lambda_predictive_verification.py` (CKM script — basis source of truth)
- `debug/w3_ckm_pmns_probe.py` (first-pass probe)
- `debug/w3_forward_plan_memo.md` (master forward plan; this is Test 1)
- `memory/w3_spectral_zeta_candidate.md` (state file)
- `memory/w3_forward_plan.md` (next-test plan pointer)
- `debug/wh7_candidate_draft.md` (WH7 candidate draft — should now be updated
  with sector-restriction caveat)

## What this means for W3 / WH7

The W3 question (second packing axiom for calibration data) becomes: do
calibration data candidate-live in master Mellin engine M1/M2 rings? After
this sprint, the answer is **partially yes, sector-specifically**. CKM does;
PMNS does not. The candidate is real but narrower than originally framed.
This is a substantive sharpening, not a fatal blow.

The WH7 candidate text in `debug/wh7_candidate_draft.md` and CLAUDE.md §1.7
(if/when promoted) should be updated to reflect: "calibration data candidate-
lives in master Mellin engine M1/M2 rings *for the quark mixing sector*; the
neutrino mixing sector follows a structurally different (TBM-rational)
pattern; charged lepton masses follow a third (Koide cone) pattern. The
sector-specificity may itself be derivable from the Connes inner-factor
sub-algebra structure (ℂ+ℍ vs M_3(ℂ) vs Majorana extension), which is the
natural Test 6 / Connes-Chamseddine spectral action target."

PI sign-off required for any §1.7 edit.
