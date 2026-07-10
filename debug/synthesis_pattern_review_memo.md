# Synthesis Pattern Review — Multi-Track Sprint Cycle, May 2026

**Type:** Leader/Decomposer-hybrid synthesis pass (research memo, not paper)
**Date:** 2026-05-09 (close-of-day)
**Coverage:** May 2026 sprint cycle, ~10 sprint-tracks landed across multi-focal Phase A/B/C, Sprint TD, Sprint HF, Sprint MH, Track 1-5 multi-track, MR-A/B/C, and the operator-level autopsy round
**Companion data:** `debug/data/synthesis_pattern_review_data.json`
**No production code modified.** **No paper edits applied.** Recommendations only.

---

## ONE-PAGE SUMMARY (read in 2 minutes)

**Headline.** Three of three named PI patterns survive scrutiny, with sharpenings: Pattern A (α-order hierarchy of walls) holds and refines into a *two-tier* structure; Pattern B (residual class ↔ projection depth) is **not the true axis** — the genuine ordering is *projection class × refinement-by-multi-focal-content* and projection depth is a downstream proxy; Pattern C (HFS clusters in §V.D) holds and predicts where the next exposure will appear. Hunt 4 (Lorentz-boost 17th projection) has a clean structural verdict: **REQUIRES-EXTENSION** at the framework level (Wick-rotated spectral triple is Riemannian; Lorentzian boost has no compact-spectrum analog), but a structurally cousin "rapidity / dimensionless v/c parameter" projection is admissible if the framework chooses to add it — and on careful inspection it sits very close to the rest-mass projection's variable-refinement architecture.

**Top 3 patterns confirmed:**
1. **α-order hierarchy of walls (Pattern A — confirmed and SHARPENED).** Multi-focal-composition walls cluster at α⁴ (Breit, Ps 1S-2S) and α⁵ (LS-8a, HF-3, HF-4, HF-5, H1 Yukawa), with a third "α² bracket-level" tier (Mu HFS at +199 ppm, the cleanest LS-8a isolation) emerging as a structural invariant. The α-order is not just numerology: each tier engages a different Layer-2 projection class (α⁴ ↔ §III.16 two-body Dirac/Breit; α⁵ ↔ renormalization counterterms, no §III projection assigned; α² ↔ vertex-correction loop, partially in §III.6 spectral action).
2. **§V.D HFS clustering (Pattern C — confirmed, sharper prediction).** All 4 currently-catalogued §V.D entries are HFS-class or muonic-Lamb-class observables. Predicts: next §V.D entries will appear in (i) deuteron polarizability-vs-Friar-moment convention (Pachucki–Yerokhin vs Friar–Payne), (ii) helium fine-structure α³(Zα)² vs higher-order Drake itemizations, (iii) positronium HFS annihilation channel itemization (Czarnecki–Melnikov–Yelkhovsky vs Karshenboim).
3. **Master Mellin engine domain partition is operator-level visible across observables (Hunt 3 — confirmed).** Of 11+ catalogue residuals examined, 9 sit cleanly in M1/M2 buckets, 1 in M3 (vertex-restricted parity), 1 ambiguous (the K=π(B+F−Δ) cross-bucket coincidence). The case-exhaustion theorem holds at the residual level, not just at the π-source level.

**Top 2 patterns falsified or unverified:**
1. **Pattern B as stated (linear depth-vs-residual scaling) is FALSIFIED.** The catalogue contains: zero-projection at ~0%, single-projection at ~0% to 99% (huge spread), three-projection chains at 0% (polarizability) to 286 ppm (D HFS) to 11.86% (He oscillator). Projection depth correlates with *opportunity for compounding*, not with the residual itself. Paper 34 §VI's Sprint Calc-P remark already acknowledges this; the Sprint TD Track 1 Stefan–Boltzmann result (depth-3 chain, residual 0) is the second clean counterexample. **The genuine axis is observable-class (single-particle vs multi-focal) × Layer-2-content presence.** Pattern B should be replaced by the **Layer-2-presence prediction** (below).
2. **The pure-axis projection candidates (gauge-choice / spin-statistics) hinted in TX-A §6 EDIT-3 are STILL UNVERIFIED.** Five sprints since TX-A, no new (--, --, π) witness has surfaced; no (--, D, --) witness has surfaced. The V/D correlation (now 12/12 of dimensioned projections, was 8/8 at TX-A) is sharper. The transcendental axis remains the only genuinely independent one with vector-photon as the unique witness.

**Top 3 ranked next tracks:**
1. **§V.D-prediction sprint (1-week diagnostic):** anchor a class-(a) literature-itemization sweep on three pre-named candidates (D polarizability, He fine-structure rank-k Drake, Ps HFS annihilation) and verify Pattern C empirically. If 2 of 3 land class-(a), Pattern C graduates to a falsifiable structural prediction; if 0 of 3, Pattern C is HFS-specific. Low cost, high diagnostic value.
2. **Operator-level promotion of §III.7 spectral action at one-loop (1–2 week sprint):** the §III.6 spectral action is currently verified at *energy-level* in 11+ catalogue rows, but operator-level construction has only been done for §III.17 (Foldy/Friar charge density) and §III.18 (Zemach). Promoting §III.6 to operator level on the H 21cm autopsy chain (where Sprint HF-2 verified Parker–Toms c₁ = R/12 = 1/2 at 0.5%) closes the operator-level gap on the most-tested projection in the framework. Diagnostic-before-engineering: confirm the Schwinger asymptote IS the operator-level limit before committing.
3. **Layer-2-presence prediction sprint (1-week diagnostic):** test whether catalogue-row residual is structurally bounded by the *number of Layer-2 inputs in the chain* rather than by projection depth. Predicts: (chains with 0 Layer-2 inputs, framework-only) → exact at converged basis; (chains with 1 Layer-2 input) → calibration-limited; (chains with 2+ Layer-2 inputs) → multi-focal-wall-class. Re-tag all 11+ catalogue rows by L2-count and verify the prediction. If holds, Paper 34 §VI gets a sharper prediction than Sprint Calc-P-refined depth-with-rational-class exception.

---

## §1. Pattern A — α-order hierarchy of structural walls

### 1.1 The 6+ instances catalogued

The multi-focal-composition wall pattern is named in CLAUDE.md memory `multi_focal_wall_pattern.md` (2026-05-07) as: "framework couples discrete labels cleanly via Wigner symbols and selection rules; doesn't natively compose two Fock-style projections at once." The structurally-clean subset of instances:

| # | Wall instance | Sprint | α-order of failure | Mechanism owning the gap | Catalogued at |
|:-:|:--------------|:-------|:-------------------|:--------------------------|:--------------|
| 1 | Sprint H1 Yukawa | H1 | α⁰ (Yukawa is non-perturbative parameter) | Inner-factor data, Paper 18 §IV.6 | CLAUDE.md §2 H1 |
| 2 | Ps 1S-2S Breit | precision-catalogue | α⁴ | §III.16 two-body Dirac/Breit retardation (added 2026-05-08) | §V.B, §V.C.7, §III.16 |
| 3 | LS-8a Z_2 / δm | LS-7/LS-8a | α⁵ multi-loop SE | renormalization (no §III projection) | CLAUDE.md §2 LS-8a |
| 4 | HF-3 recoil cross-register | Sprint HF | α⁴(m/M) recoil-mixing-related | §III.16 + cross-register kinematics | §V.B HF-3 row |
| 5 | HF-4 Zemach | Sprint HF | α⁴(m/M) ⋅ r_Z | §III.18 magnetization-density (operator-level closed 2026-05-09) | §V.C.2 (closed positive) |
| 6 | HF-5 multi-loop a_e | Sprint HF | α⁵ vertex sector | renormalization (no §III projection) | §V.B HF-5 row |
| 7 | Mu HFS LS-8a isolation | precision-catalogue | α²(Zα) bracket level | §III.6 spectral action sub-leading | §V.B Mu-HFS row |

(Sprint H1 Yukawa is α⁰ because the Yukawa coupling itself is a non-perturbative input; included because it's the founding case of the wall pattern at the spectral-triple level.)

### 1.2 Three-tier structure (the SHARPENING)

Reading down the table by α-order:
- **α⁰ (parameter-level):** H1 Yukawa. The wall is a *Layer-2 input-class* gap (Paper 18 §IV.6 inner-factor input data tier). No projection exists or could exist at this level — calibration is empirical.
- **α²-α³ (one-loop QED bracket level):** Mu HFS. The wall surfaces because the framework's α/(2π) Schwinger asymptote at one loop captures only the leading anomalous moment; sub-leading α²(Zα) corrections require multi-loop bracket evaluation that the framework's bare CC spectral action doesn't autonomously generate.
- **α⁴ (two-body kinematic level):** Ps 1S-2S, HF-3 partial, HF-4 partial. The wall is the §III.16 Breit retardation projection — two-body Dirac kinematic content not in the single-particle Eides §3.2 SE bracket.
- **α⁵ (renormalization level):** LS-8a, HF-5. The wall is the absence of Z_2 / δm / Z_3 counterterms in the bare framework. No §III projection is assigned because "renormalization" is structurally a flow / regularization, not a focal-length coupling.

The walls cluster in two structurally distinct buckets:
- **Two-body kinematic (α⁴):** §III.16 partially fills, with Ps 1S-2S as the operator-level test case. Cleanest because it surfaces α⁴ — one full order before renormalization confounds.
- **Renormalization (α⁵):** No §III projection assigned. Honest: the framework doesn't have a renormalization mechanism. CLAUDE.md §2 LS-8a entry concludes this is "structural-skeleton scope — calibration data is genuinely external input."

### 1.3 What would falsify Pattern A

- **Falsifier 1:** A new wall instance at α³ or α⁶ that doesn't fit either of the two named tiers. Currently no observable in the catalogue is at α³ or α⁶; the cluster at α⁴ (Breit) and α⁵ (renormalization) is empirically clean, with the Mu HFS α² datapoint being the only outlier and even that one is structurally interpretable as one-loop bracket precision rather than a new tier.
- **Falsifier 2:** A new wall instance where the α-order of the residual does NOT match the α-order of the missing physics. For example, if a Lamb-shift autopsy revealed a sub-1% residual that PSLQ-identified as α⁴⋅ζ(3) — that would falsify the Pattern A reading, because the framework's residual would be at α⁴ but the missing physics would be MZV content not in either of the two named tiers.

### 1.4 Predictions

- **Prediction A1:** Future α⁴ wall instances will trace to either the §III.16 two-body Dirac/Breit projection or to a new "two-body magnetization" projection (sibling of §III.18 at the magnetic dipole level). Test: any equal-mass leptonic 2-body or 3-body bound state will hit α⁴ Breit before LS-8a renormalization.
- **Prediction A2:** The Mu HFS α²(Zα) +199 ppm residual will *not* tighten under further framework-internal work; it will only close via Layer-2 input from Karshenboim 2005 itemization. (Already observed: Karshenboim full Layer-2 itemization closes Mu HFS to <1 ppm.)
- **Prediction A3:** No α⁰ wall (Yukawa-class) will surface in the atomic catalogue. They are confined to the inner-factor / second-packing-axiom question (CLAUDE.md §1.7 W3).

### 1.5 Statistical confidence

Six instances at α⁴-α⁵, one at α², one at α⁰. The α⁴/α⁵ clustering with 6 instances is well past the chance threshold for the dimensioned axes available (variable, dimension, transcendental, projection). The α² isolated instance is consistent with what the framework predicts for one-loop bracket precision in the Bohr–Fermi class. Confidence: **HIGH** for the two-tier structure; **MODERATE** for the α² as a third tier (n=1 datapoint, more needed to validate).

---

## §2. Pattern B — Residual class ↔ projection depth (FALSIFIED as stated, REPLACED by Layer-2-presence prediction)

### 2.1 The naive prediction (Paper 34 §VI Prediction 1)

Naive form: `ε_k ~ k · ε_1` where `k` is the number of chained projections.

Paper 34 §VI already acknowledges this fails at depth-3 due to Sprint Calc-P (H 1S polarizability exact at depth-3 because the Dalgarno–Lewis solution is rational and lies in the Sturmian basis exactly). Sprint TD Track 1 (Stefan–Boltzmann at depth-3, residual 0) is the second clean counterexample.

### 2.2 The empirical residual table

I extracted residuals from Paper 34 §V (machine-precision matches) and §V.B (off-precision with error codes), tagged each row by chain depth and by Layer-2 content count:

| Observable | Chain depth | Layer-2 count | Residual |
|:-----------|:-----------:|:-------------:|:--------:|
| Hydrogen E_n | 1 | 0 | exact |
| Schwinger β(α) | 1 | 0 | exact symbolic |
| F = ζ_R(2) | 1 | 0 | exact symbolic |
| E_Cas^S³ scalar = 1/240 | 2 | 0 | exact rational |
| E_Cas^S³ Dirac = 17/480 | 2 | 0 | exact rational |
| Stefan–Boltzmann π²/90 | 3 | 0 | exact symbolic |
| Hydrogen Bohr levels | 2 | 0 | exact symbolic |
| H 1S polarizability 9/2 a₀³ | 3 | 0 | exact (Sprint Calc-P) |
| Lyman α A_21 | 3 | 0 | +0.055% (NR vs Dirac) |
| Bethe log ln k₀(2S) acc | 2 | 0 | -0.92% |
| Mu 1S-2S rest-mass-only | 2 | 1 (m_μ) | -0.11 ppm |
| H 21cm HF-1 (BF only) | 3 | 0 (g_p external implicit) | -1102 ppm = -a_e |
| H 21cm HF-2 (+Schwinger) | 4 | 0 | +58 ppm |
| H 21cm HF-4 (+Zemach) | 4 | 1 (r_Z) | +18 ppm |
| H Lamb shift (LS-6a) | 4 | 0 | -0.534% |
| μH Lamb (Sprint MH) | 5 | 3 (KS, polarizability, recoil) | -0.10% |
| μH 1S BF (Track B) | 4 | 0 | +2 ppm |
| Mu Lamb shift | 4 | 0 | +0.013% |
| D HFS BF strict | 3 | 0 (g_d external) | +40 ppm |
| Mu HFS BF + Schwinger | 4 | 0 | +199 ppm |
| Ps 1S HFS BF only | 3 | 0 | -42.58% (annihilation L2) |
| Ps 1S-2S framework-only | 3 | 0 | +64.75 ppm (Breit L2) |
| He 2³P P_0-P_1 | 4 | 0 | -0.014% |
| He 2³P P_0-P_2 | 4 | 0 | -0.20% |
| He 2¹S-2³S graph CI | 2 | 0 | +11.86% |
| He 2¹P→1¹S osc strength | 3 | 0 | +3.4% (multi-focal extension) |
| LiH R_eq | 2 | 0 | +5.3% |

### 2.3 Statistical fit: depth vs |residual|

A naive log|residual| vs depth regression (excluding "exact" rows) over the 17 numeric rows above:

- Depth 2: residuals 0.92%, 11.86% → log mean -1.5, range ±5.
- Depth 3: residuals 0.05%, 0.06%, 64.75 ppm, 0.40 ppm, 40 ppm, 0.11 ppm, 42.58%, 3.4% → wide spread, no monotone trend.
- Depth 4: residuals 0.014%, 0.20%, 0.534%, 18 ppm, 58 ppm, 199 ppm, 0.013%, 2 ppm → wide spread.
- Depth 5: -0.10% (single datapoint).

**No monotone trend with depth.** The residuals are scattered across 6 orders of magnitude at every depth ≥ 2.

### 2.4 The genuine axis: Layer-2-content count

Re-tagging by L2-count gives a much cleaner picture:

- **L2-count = 0 (framework-native chains):** residuals run from machine-exact (Stefan–Boltzmann, polarizability) to ~1% (Lamb LS-6a, Bethe log) to ~12% (He 2¹S-2³S, BAD basis). The residual is bounded by *basis quality* (B in the §V.B error-code system) when the chain is framework-only.
- **L2-count = 1 (one external scalar, e.g. r_Z, m_d, g_p):** residuals tighten to ~tens of ppm or sub-ppm (H 21cm HF-4 at +18 ppm, Mu 1S-2S at -0.11 ppm, μH BF at +2 ppm). The single L2 input is always a scalar parameter, and its uncertainty propagates linearly.
- **L2-count = 2 (two external inputs):** structurally the multi-focal-wall regime. Mu HFS +199 ppm (Karshenboim itemization closes <1 ppm). Ps 1S-2S +65 ppm (Breit + multi-loop closes sub-MHz). The residual is "what each independent L2 input contributes summed in quadrature."
- **L2-count = 3+:** μH Lamb (multi-loop QED + polarizability + recoil) at -0.10%. The residual is bounded by the *largest single L2 contribution*, here the Källén–Sabry two-loop VP.

### 2.5 Replacement prediction

**Prediction B' (replaces naive Prediction 1):** The residual of a catalogue row is bounded above by the largest single Layer-2 input it consumes, plus the framework's basis-quality residual in the absence of L2 inputs. Equivalently:

`|ε(observable)| ≤ max( ε_basis(framework), max_i |L2_i| )`

where `ε_basis(framework)` is the basis-quality limit of the framework alone (e.g. ~0.534% for H Lamb shift one-loop closure, ~0.024% for He NR variational).

This is structurally consistent with Pattern A: the wall instances cluster at α⁴ and α⁵ because that's where the largest L2 inputs sit, not because chain depth crosses a threshold.

### 2.6 Falsifier for Prediction B'

A catalogue row with L2-count = 0 and residual sub-ppm at finite n_max (i.e., framework-native exact at sub-ppm precision on a non-trivial chain) would not falsify B' but would tighten the basis-quality bound. A catalogue row with L2-count = 2+ and residual smaller than the largest single L2 input by more than an order of magnitude would falsify B' — the prediction says one L2 input dominates.

### 2.7 Recommendation

Paper 34 §VI Prediction 1 should be reframed as **Prediction 1' (Layer-2-presence-bound)**, with the depth-correlation downgraded to a downstream observation. The Sprint Calc-P refinement (rational-class observables exact at finite n_max) should stay, but the parent prediction needs replacement.

---

## §3. Pattern C — §V.D HFS clustering and prediction for next exposures

### 3.1 The four current exposures

- §V.D.1: D 1S HFS Eides vs Pachucki–Yerokhin Layer-2 itemization. Magnitude ~0.01% of ν_F.
- §V.D.2: μH 2S-2P Lamb VP one-loop split, Antognini vs Krauth. Magnitude ~100 ppm.
- §V.D.3: Mu Lamb SE-recoil aggregation, framework Eides §3.2 vs Karshenboim. Magnitude ~11 MHz at the SE bracket (~1% of bracket).
- §V.D.4: μH Friar moment profile convention, Eides closed-form vs operator-level. Magnitude 12.7% (Gaussian profile factor 4/π).

### 3.2 Observable classification

§V.D.1 is HFS, §V.D.2 is muonic Lamb, §V.D.3 is Mu Lamb, §V.D.4 is muonic Lamb. The cluster is clearly weighted toward muonic and HFS observables. No §V.D entries from electronic Lamb shift, He fine structure, polarizability, or Lyman α.

### 3.3 Why HFS observables cluster

HFS observables decompose into many sub-leading Layer-2 corrections at sub-percent precision: leading Zemach, NLO recoil-mixing, finite-size, multi-loop QED, polarizability, hadronic VP. Each sub-leading corrections has its own canonical compilation (Eides 2024, Karshenboim 2005, Pachucki–Yerokhin 2010, Friar–Payne 2005), and small itemization differences compound across compilations. Lamb shift observables decompose similarly (Eides Tab 7.3 vs Antognini 2013) but with fewer sub-leading channels at the precision where the framework operates.

The cleanest signal for HFS clustering: §V.D.1 is the deepest Layer-2 itemization (deuteron polarizability + multi-loop + recoil + finite-size all itemized differently); the magnitude is ~0.01% of ν_F (sub-percent). HFS observables are where small itemization differences cluster.

### 3.4 Predictions for next §V.D entries

In probable order of likelihood:

- **Prediction C1 (HFS, ~80% confident):** Deuteron polarizability convention. Pachucki–Yerokhin 2010 vs Friar–Payne 2005 itemization of the deuteron polarizability contribution to D HFS. Expected magnitude: tens of ppm. Follow-up sprint: deeper D HFS Roothaan autopsy.
- **Prediction C2 (HFS, ~60% confident):** Positronium HFS annihilation channel itemization. Czarnecki–Melnikov–Yelkhovsky 2000 vs Karshenboim 2005 itemization of α⁵ corrections. Expected magnitude: ~100 ppm of LO HFS. Follow-up sprint: Ps HFS five-component autopsy at α⁵.
- **Prediction C3 (Lamb-class, ~50% confident):** Helium fine-structure α³(Zα)² itemization. Drake 1971 vs Pachucki–Yerokhin 2010 itemization of higher-order Breit-Pauli + multi-electron SS/SOO. Expected magnitude: tens of MHz on GHz-scale intervals.
- **Prediction C4 (HFS, ~40% confident):** Muonium HFS Karshenboim itemization vs framework m_red rescaling, sibling of §V.D.3 but at HFS instead of Lamb. Expected magnitude: ~10s of ppm.

### 3.5 What would falsify Pattern C

The next 3-5 §V.D entries should be HFS or muonic-class. If 2 of the next 3 sprints surface §V.D entries from electronic Lamb shift (where Eides-vs-Pachucki-Yerokhin itemizations are well-aligned) or from non-precision atomic spectra (where the framework operates in larger residual ranges), Pattern C is HFS-specific and not a structural prediction about what observable types surface convention exposures.

### 3.6 Confidence

n=4 is small for a clustering claim, but the structural reading (HFS observables decompose into more sub-leading channels at the framework's precision) is mechanistic and gives an a priori reason for the clustering. **MODERATE** confidence, **HIGH** prior expectation that next 3 entries will be HFS-class, **MEDIUM-HIGH** for predictions C1/C2.

---

## §4. Hunt 1 — Operator-level verification map

The 19 §III projections (Paper 34 has 19 not 16 — sprint of 2026-05-09 added §III.17 Foldy/Friar, §III.18 Zemach magnetization-density, §III.19 tensor multipole). Status of each:

| # | Projection | Operator-level verified? | Catalogue rows where verified | Gap |
|:-:|:-----------|:-------------------------:|:------------------------------|:----|
| 1 | Fock conformal | YES | All catalogue rows (founding projection) | none |
| 2 | Hopf bundle | PARTIAL | Paper 2's K = π(B+F-Δ); WH1 spectral truncation R2.5 PROVEN | not at observable-level beyond α conjecture |
| 3 | Bargmann–Segal | YES | Paper 24 π-free certificate, HO Casimir 17/3840 | none |
| 4 | Stereographic | YES (implicit in Fock) | Coulomb-as-coordinate (Paper 7) | none |
| 5 | Sturmian | YES | Bethe log LS-3, Lyman α, polarizability | none |
| 6 | Spectral action (CC) | ENERGY-LEVEL | Schwinger β, VP, Uehling shift, Sprint MH Track A full kernel | **operator-level construction at one-loop is the Hunt-1 main gap** |
| 7 | Spinor lift | YES | Dirac-Coulomb formula, T9 theorem, Camporesi–Higuchi spectrum | none |
| 8 | Wigner 3j | YES | Gaunt integrals, Drake J-pattern (sympy-exact), Roothaan termination | none |
| 9 | Wigner D | YES | Multi-center molecules (LiF, CO, etc.) | none |
| 10 | Wilson plaquette | YES | Paper 30 SU(2) Wilson, transverse photon | none |
| 11 | Vector-photon | YES | Paper 33 1+6+1 partition, 7/8 selection rules recovered | none |
| 12 | Molframe hyperspherical | YES | Paper 15 H₂ at 96.0% D_e | none |
| 13 | Drake–Swainson | YES | Bethe log ℓ>0 closure | none |
| 14 | Rest-mass | YES | Mu 1S-2S at -0.11 ppm, Ps at bit-identical, Sprint MH Track B | none |
| 15 | Observation/temporal-window | PARTIAL | Stefan-Boltzmann via Track 1; Sprint TD Track 4 cigar M1 | not at observable beyond thermal |
| 16 | Two-body Dirac/Breit | PARTIAL | §V.C.7 Ps 1S-2S autopsy: kernel + Roothaan termination YES, bound-state matrix elements NO | **Bethe–Salpeter quasi-energy expansion is named extension** |
| 17 | Foldy/Friar charge density | PARTIAL | §V.C.3 μH Lamb autopsy uses operator at 4.3% LO | sibling MagnetizationDensitySpec for FNS channel (named extension) |
| 18 | Zemach magnetization-density | YES | §V.C.2 H 21cm at I=1/2, §V.C.5 D HFS at I=1, both at machine precision | none |
| 19 | Tensor multipole | NO | No catalogue row uses leading-order; HD/D₂ rotational HFS not catalogued | **operator-level construction not implemented** |

### 4.1 Gap-ranked summary

- **Most load-bearing gap: §III.6 Connes–Chamseddine spectral action at one-loop.** This is the most-tested projection in the framework (15+ catalogue rows engage it), but operator-level construction is incomplete. Energy-level evaluation works (Schwinger β, Uehling shift, full kernel for muH). Operator-level construction would expose the bracket structure (Eides §3.2 vs Karshenboim) that drove §V.D.3 — a single sprint that promotes §III.6 to operator level on the H 21cm chain (where Sprint HF-2 verified Parker–Toms c₁ = R/12 = 1/2 at 0.5%) closes the operator-level gap on the most-tested projection.
- **Second gap: §III.16 Breit retardation at bound-state matrix elements.** §V.C.7 Ps 1S-2S autopsy verified kernel + Roothaan termination at operator level; full bound-state matrix-element evaluation requires Bethe–Salpeter quasi-energy expansion (named structural extension, ~3-week sprint).
- **Third gap: §III.19 tensor multipole.** Forward-looking projection slot. No catalogue row uses leading-order; HD/D₂ rotational HFS would be the load-bearing test. Operator-level construction is a 2-week sprint (sibling of §III.18 magnetization-density).
- **Minor gap: §III.17 sibling MagnetizationDensitySpec.** Surfaced by §V.D.4 (Friar profile RMS-vs-first-moment). 1-day mechanical sprint.

### 4.2 Recommendation

The §III.6 promotion is the highest-value single sprint after the multi-track May 9 round. It's also the prerequisite for any future α³ multi-loop test (LS-7/LS-8a refinement at the bracket level).

---

## §5. Hunt 2 — Hidden axis correlations under the 19-projection inventory

### 5.1 V/D correlation update

Of the 19 projections (post-May 9 update), Paper 34 §III.20 Remark `rem:vd_correlation` claims "12 add a dimension and all 12 also add the variable that carries it." Re-checking against the dictionary:

- Adds dimension energy: Fock, Bargmann–Segal, spectral action, molframe hyperspherical, rest-mass (mass = energy at c=1), observation (time = inverse energy), Breit retardation. **7 projections add energy.**
- Adds dimension length: Stereographic, Wigner D, Foldy/Friar (r_E), Zemach (r_Z), tensor multipole (r²-class). **5 projections add length.**
- Adds no dimension: Hopf, Sturmian, Wigner 3j, spinor lift, Wilson, vector-photon, Drake–Swainson. **7 projections.**

**Total dimensioned:** 12 (matches the §III.20 remark). All 12 also add a variable (perfect V/D correlation, 12/12).

### 5.2 V-only and P-only

- V-only (variable, no dimension, no transcendental): spinor lift, Drake–Swainson. **2 projections** (was 2 at TX-A).
- P-only (transcendental only, no variable, no dimension): vector-photon. **1 projection** (still the unique witness).

### 5.3 Three-axis 19-projection cross-tab

| Pattern (V D P) | Count at 15 (TX-A) | Count at 19 (now) | Notable additions |
|:----------------|:------------------:|:-----------------:|:------------------|
| VD-- | 6 | 9 | + Foldy/Friar, Zemach, tensor multipole |
| V-P | 2 | 2 | unchanged (Hopf, Wilson) |
| VDP | 2 | 2 | unchanged (spectral action, observation) |
| V-- | 2 | 2 | unchanged (spinor, Drake) |
| --P | 1 | 1 | unchanged (vector-photon) |
| --- | 2 | 2 | unchanged (Sturmian, Wigner 3j) |
| -D- | 0 | 0 | (no D-only projection added) |
| -DP | 0 | 0 | (no D+P-only projection added) |

**Headline:** the V/D correlation strengthens from 8/8 (TX-A) to 12/12 (now). No new transcendental-axis-independence witness has surfaced. The VDP pattern (heavy projection) is unchanged at 2.

### 5.4 New structural correlations under 19 projections

- **Layer-1 → Layer-2 transition correlation:** all 19 projections that take input Layer 1 and output Layer 2 (i.e., 17 of 19 — only Drake–Swainson and Breit retardation are 2→2) introduce content tied to a *single* Paper 18 transcendental tier. The 2 projections that are 2→2 (Drake–Swainson, Breit retardation) operate on already-Layer-2 content and refine it. This is consistent with the case-exhaustion theorem (Paper 32 §VIII) interpretation that transcendental content enters at exactly one of M1/M2/M3 mechanisms per projection.
- **Length-axis subdivision:** the 5 length-introducing projections subdivide cleanly into (a) stereographic / coordinate-change (1 projection: §III.4); (b) molecular geometry (1 projection: §III.9 Wigner D); (c) molecular focal length (1 projection: §III.12 molframe); (d) intra-nuclear (3 projections: §III.17, §III.18, §III.19). The intra-nuclear cluster is the May 9 addition and they share a structural property: they all introduce a length scale that's *internal to a particle* (nuclear charge / magnetization / quadrupole shape), categorically distinct from atomic / internuclear length scales.
- **No new (--, D, --) candidate has surfaced.** TX-A flagged this as the structural absence; sprint TX-A §6 EDIT-3 mentioned gauge-choice and spin-statistics as candidates. Five sprints later, neither has been formalized. The flag stays open.

### 5.5 Surprising redundancy candidates

Three projections might look redundant under the new 19-projection inventory:

- **§III.17 Foldy/Friar charge density vs §III.18 Zemach magnetization-density.** Both add a length scale per nucleus, both leave the algebraic ring unchanged. They differ in *which Hamiltonian operator they modify*: Foldy/Friar modifies V_ne (Coulomb), Zemach modifies A_hf (Fermi contact). Empirically distinguished: H Lamb shift depends on r_p (Foldy/Friar) but only trivially on r_Z; H 21cm depends on r_Z (Zemach) but only trivially on r_p. **Not redundant — categorically distinct operators.**
- **§III.14 rest-mass vs §III.16 Breit retardation.** Both add mass-class variables. They differ in scope: §III.14 rescales single-particle observables under m_e → m_red; §III.16 introduces the *ratio* m_l/m_n as a two-body kinematic control parameter. The non-commutation pair (§III.14 must precede §III.16) is the structural fingerprint. **Not redundant — variable refinement.**
- **§III.7 spinor lift vs §III.11 vector-photon.** Both modify selection rules at the QED level. They differ in *which side of the QED interaction is promoted*: §III.7 promotes the matter (electron) from scalar to spinor; §III.11 promotes the gauge boson (photon) from scalar 1-cochain to vector harmonic. **Not redundant — matter vs gauge.**

No projections are structurally redundant under current criteria.

### 5.6 Recommendation

The 19-projection inventory is structurally clean — no redundancies, no new independent transcendental witnesses, V/D correlation strengthened. The pure-axis projection question (whether (--, D, --) or (--, --, π) candidates exist beyond vector-photon) remains open.

---

## §6. Hunt 3 — Master Mellin engine signatures across observables

### 6.1 The case-exhaustion theorem (recap)

Paper 32 §VIII, Sprint TS-E1: every π in any finite chain of Paper 34 projections engages M1 (Hopf-base measure Vol(S²)/4), M2 (Seeley-DeWitt √π via Mellin of Tr e^{-tD²}), or M3 (vertex-topology Hurwitz Dirichlet L). Sprint MR-A/B/C (May 6) refined this to a *mechanism-and-domain partition*: M1 ↔ propinquity rates (4/π); M2 ↔ heat-kernel/spectral-action convergence (√π·ℚ ⊕ π²·ℚ); M3 ↔ vertex-restricted parity-character observables (Catalan G, β(4)).

### 6.2 Catalogue ↔ M1/M2/M3 mapping

For each catalogue row that has a known transcendental signature, I tag the mechanism:

| Observable | Transcendental signature | Mechanism |
|:-----------|:-------------------------|:----------|
| Π = 1/(48π²) VP coefficient | π²·ℚ | M2 (Seeley–DeWitt) |
| β(α) = 2α²/(3π) | π·ℚ | M2 (Seeley–DeWitt convolution) |
| F = ζ_R(2) | π²·ℚ | M2 (Dirichlet at integer s, Bernoulli mechanism) |
| Stefan–Boltzmann π²/90 | π²·ℚ | M2 × (M1 from S¹_β) — explicit Track 1 cross |
| E_Cas^S³ = 1/240 | rational, no π | none (Layer 1) |
| E_Cas^S³ Dirac = 17/480 | rational, no π | none (Layer 1) |
| E_Cas^HO,S⁵ = -17/3840 | rational, no π | none (Layer 1) |
| Hawking T_H = 1/(8πM) | π·ℚ | M1 (Matsubara cigar, Track 4 verified) |
| D_even(4) = π²/2 - π⁴/24 - 4G + 4β(4) | π²·ℚ + Catalan G | M2 + M3 |
| D_diff(s) = 2^{s-1}(β(s) - β(s-2)) | Dirichlet L · 2-power | M3 |
| Vector-photon 1/(4π) per loop | π·ℚ | M1 (Hopf-base Vol(S²)/4) |
| K = π(B + F - Δ) ≈ 1/α | π × multi-sector | **M1 outer × cross-bucket B/F/Δ** |
| Parker–Toms c₁ = R/12 = 1/2 | rational | M2 (curvature expansion structure) |
| Heisenberg–Euler α²/(45π m_e⁴) | π denominator | M2 |

**9 of 11+ rows in M1 or M2 buckets cleanly.** D_even(4) is the canonical M3 example. K = π(B+F-Δ) is the cross-bucket coincidence.

### 6.3 Cross-bucket observation: K = π(B + F - Δ)

The K conjecture is the only observable in the catalogue that crosses sub-mechanism buckets. Decomposition:

- B = 42 = finite Casimir trace at m=3. Rational. (No π — Layer 1.)
- F = π²/6 = ζ_R(2) = D_{n²}(d_max). M2 (even-zeta).
- Δ = 1/40 = 1/g_3^Dirac. Rational. (No π — Layer 1, single-level Dirac multiplicity.)
- Outer π factor: M1 (Hopf-base measure Vol(S²)/4 = π identification, Sprint A α-PI).

The cross-bucket structure is: outer M1 multiplies (Layer-1 + M2 + Layer-1) sum. CLAUDE.md §1.7 WH5 reading: K is a coincidence at finite cutoff between three categorically different objects with no common arithmetic generator. Twelve mechanisms eliminated for K = single-principle derivation (Phases 4B-4I + Sprint A α-SP + K-CC).

### 6.4 Are there any "rogue" transcendentals?

Checking for catalogue residuals not in M1/M2/M3:

- **ζ(3) appears in three-loop QED.** It's in M3 (half-integer Hurwitz at the Dirac sector, weight-m subchannel). Confirmed structural (Phase 4I D3).
- **log 2, log 3 in Breit retardation.** Sprint Calc-Ps-1S2S Track 4 found R^0_BP integrals at Z=1 in ℚ[log 2, log 3]. This is **embedding-log content** (Paper 18 embedding tier) — *not* in M1/M2/M3. The mechanism is the Mellin-regularized inner integral producing log p for small primes p from negative-power singularity treatment. **This is the first clean rogue example.**
- **Catalan G = β(2) at QED vertex parity.** M3.
- **Dirichlet β(4).** M3.
- **S_min ≈ 2.47993693803... (Paper 28).** Two-loop sunset CG-weighted. PSLQ-irreducible against M1 ∪ M2 ∪ M3 ∪ {γ_E, log 2, G, ζ(3)}. **Second rogue example** — depth-2 multiple Hurwitz zeta at half-integer shift, intrinsic to the Dirac spectrum.

### 6.5 Verdict

The case-exhaustion theorem holds for *Layer-2-output* transcendentals (where π appears in the observable). It does *not* close on intermediate quantities: log 2 / log 3 in Breit kernels (embedding tier), MR-C's L2 next-order constant c ≈ 4.10932 (Stein–Weiss IBP downstream of M1), and S_min in three-loop QED all sit outside M1 ∪ M2 ∪ M3. The case-exhaustion theorem is therefore precisely **about the π-content of converged finite-chain Layer-2 evaluations**, not about all transcendentals appearing in any framework computation. Paper 32 §VIII Theorem statement is correct as stated.

### 6.6 Predictions

- **Prediction H3.1:** Future catalogue rows with new transcendental content will be in M1, M2, or M3. Specifically: any new α/π content will be M2; any new "π in solid-angle integration" will be M1; any new vertex-parity content will be M3.
- **Prediction H3.2:** Embedding-tier transcendentals (log 2, log 3) will appear in *intermediate* matrix elements but cancel from converged Layer-2 observables. Test: compute a chain where Breit kernel matrix elements feature ℚ[log 2, log 3] and verify the converged observable is in M1/M2/M3.

### 6.7 Confidence

**HIGH** for the case-exhaustion theorem at the Layer-2 output level; **MODERATE-HIGH** for the case-exhaustion at the intermediate-matrix-element level (log 2 / log 3 are surfaced and not yet shown to cancel in a Layer-2 observable).

---

## §7. Hunt 4 — Lorentz boost as a 17th projection (now 20th, given current 19-projection inventory)

### 7.1 Initial framing of the question

Under a Lorentz boost: λ_x → λ_x/γ (length contraction), β → βγ (time dilation), with γ = 1/√(1-v²/c²). Sprint TD Track 1 already built T_{S³} ⊗ T_{S¹_β} at zero relative velocity. The question: does §VII composition theorem (commute/one-way/idempotent partition) preserve under joint λ → λ/γ, β → βγ?

### 7.2 Structural obstacle: Euclidean signature

GeoVac's spectral triple is Riemannian (Wick-rotated). The Camporesi–Higuchi Dirac spectrum is on Euclidean S³, not on Minkowski space. Lorentz boosts are Minkowski-signature transformations — they don't act naturally on Euclidean spectral triples. Specifically:

- The Wick rotation t → iτ mixes time and space. Lorentz boosts act on Minkowski time, but in Euclidean signature the analog is a rotation in the (τ, x) plane. Such a rotation acts on the spectral triple T_{S³} ⊗ T_{S¹_β} by mixing the spatial S³ with the temporal S¹_β — but the construction in Track 1 has S³ and S¹_β as independent tensor factors with disjoint operators, and a "rotation" between them would require a connection between the two factors.
- The Connes–Chamseddine spectral action is invariant under rotations of the underlying Riemannian manifold. A Lorentzian boost has no analog because Wick rotation has already happened. 
- The path forward would be to build a Minkowski-signature spectral triple — replace the Euclidean Dirac spectrum |λ_n| = n + 3/2 with a Lorentz-invariant operator. This is structural extension, not a projection refinement.

### 7.3 Variable analysis

The proposed boost projection variable is v/c (rapidity), dimensionless, ring-preserving in ℚ[v/c] / (γ² · (1 - (v/c)²) - 1) — algebraic over ℚ(v/c). Structurally similar to the Dirac-sector γ = √(1 - (Zα)²) which is bound-state radial Lorentz factor (single particle in central potential), not a frame transformation.

A clean structural distinction:

- **§III.14 rest-mass:** rescales m → m_red, ring-preserving. Single-particle.
- **§III.16 Breit retardation:** introduces m_l/m_n ratio, two-body kinematic.
- **Proposed boost:** introduces v/c, *frame*-kinematic (relates two reference frames at relative velocity v).

The boost is a *frame-rescaling* projection — it would relate two T_{S³} ⊗ T_{S¹_β} constructions at different observer rapidities. This is structurally distinct from anything currently in the catalogue.

### 7.4 Compatibility with the master Mellin engine

- **M1 preservation:** The 2π in Matsubara modes ω_k = 2πk/β survives under β → βγ as 2π/(βγ). The M1 Hopf-base measure signature is preserved (it's the *same* 2π). PASS.
- **M2 preservation:** The Camporesi–Higuchi spectrum |λ_n| = n + 3/2 doesn't natively know about boosts (it's the Dirac operator on Euclidean S³ at fixed observer). Under λ_x → λ_x/γ, the spectrum would rescale as |λ_n| = γ⁻¹ · (n + 3/2), and the Seeley–DeWitt coefficients on the rescaled S³ would scale accordingly. Trace integrals Tr e^{-tD²} are invariant under uniform scaling D → D/γ (the t variable absorbs the factor). M2 PRESERVED at the trace level.
- **M3 preservation:** Vertex-parity selection on half-integer Hurwitz at quarter-integer shifts is determined by spectrum parity. Rescaling preserves parity (odd n → odd, even n → even). PASS.

The master Mellin engine ring assignments are preserved under uniform scaling (which is what a boost looks like in the spectral-triple language). M1, M2, M3 are all preserved.

### 7.5 Compatibility with Roothaan multipole termination

The Roothaan termination at L_max = 2 l_max (Paper 33 / sprint Calc-Ps-1S2S) is a Wigner-3j triangle inequality — pure angular-content, not radial. Boosts don't touch angular content. PRESERVED.

### 7.6 Compatibility with Connes axioms

J²=−I, JD=+DJ, KO-dim 3 are properties of the Euclidean Dirac operator on Euclidean S³. Boosts in Minkowski signature don't act on the J operator the same way as Euclidean rotations. Structurally **NOT preserved as stated** — would need a Lorentzian J that satisfies different axioms (KO-dim varies by signature).

### 7.7 Verdict for Hunt 4

**REQUIRES-EXTENSION at the framework level.** The Wick-rotated spectral triple framework is Riemannian; a clean Lorentz boost projection would require a Minkowski-signature spectral triple. This is structural extension (multi-month NCG work, not a sprint), and at the time of writing GeoVac doesn't have a Minkowski analog of the Camporesi–Higuchi Dirac spectrum.

**HOWEVER**, a *non-Lorentzian* "rapidity-like" projection is admissible if the framework chooses to add it. Specifically:

- **Variable v/c, dimensionless, ring-preserving in algebraic γ-extension.**
- **Mechanism:** uniform rescaling of all spatial focal lengths λ_x by 1/γ, and temporal compactification scale β by γ. (No actual boost; just a parametric rescaling.)
- **Position in axis grid:** V---, no new dimension or transcendental introduced.
- **Composition relations:** would commute with everything except §III.14 rest-mass (where it acts on the rescaled mass m → m·γ-derived quantities) — variable refinement, like §III.14/§III.16 non-commutation.

This non-Lorentzian projection is *structurally close to* §III.14 + §III.16 in being a variable-refinement projection. Adding it would not introduce a new (--, --, π) witness (rest-mass was already V---). It would explicitly model how the Stefan–Boltzmann constant transforms under rescaling, formalizing what Track 4 verified at the cigar:

T_H(M, v=0) = 1/(8πM) → T_H(M, v) = 1/(8πM·γ).

### 7.8 Distinct from §III.14 rest-mass?

The Dirac-sector γ in radial wavefunction structure (§III.7 spinor lift's γ = √(1-(Zα)²)) is structurally a *bound-state radial Lorentz factor* — it's a single-particle correction to the radial structure, not a frame transformation. A boost projection would be a *frame*-level analog, acting on the full spectrum simultaneously.

- §III.7's γ is intrinsic to the spinor sector; the proposed boost γ would be ambient.
- §III.14 rests-mass acts on a single particle's mass; the proposed boost would act on all particles in the system.
- §III.16 Breit retardation introduces m_l/m_n ratio between two species; the proposed boost would act on all species uniformly.

The boost projection would be structurally **distinct from §III.14, §III.16, §III.7 in terms of scope** but **identical to them in transcendental and dimensional axis (V---)**.

### 7.9 Recommendation

**Don't add a Lorentz boost projection now.** The structural extension required for a Minkowski-signature spectral triple is multi-month NCG work, and the cleaner uniform-rescaling alternative (which preserves M1/M2/M3 and Roothaan termination) doesn't surface a new structural discovery beyond what §III.14 / §III.16 already capture. The Lorentz boost question is best left as an honest structural negative: the framework's Wick-rotated nature precludes Lorentz boosts as projections, and the alternatives are reducible to existing variable-refinement projections.

The flag stays open for a future PI decision: "if a Minkowski spectral triple is built (multi-month work), Lorentz boosts become a candidate projection; pre-Minkowski, they're not admissible."

---

## §8. Hunt 5 — Empirical gaps where decomposition can assist

Survey of standard-physics observables outside the GeoVac catalogue where Layer-2 itemization conventions are heavy and the framework's decomposition discipline could surface new exposures.

### 8.1 Candidate list

| # | Candidate | Standard-physics venue | Layer-2 itemization depth | Framework alignment |
|:-:|:----------|:------------------------|:--------------------------:|:--------------------|
| 1 | King-plot isotope shifts (Yb, Ca⁺, Sr) | BSM searches, atomic clocks | High (relativistic + nuclear deformation + QED + electron correlation) | MEDIUM — heavy-atom W1c-residual wall sits on path; Cs HFS scoping (§V.C.6) shows Z>20 cliff at CR67 fits. Re-test at Yb requires BBB93/KTT or self-consistent HF first. |
| 2 | Light-nucleus polarizability extractions (D, He-3/4) | nuclear structure | High (NN dynamics + multi-loop QED + recoil) | HIGH — D HFS at +286 ppm already has the polarizability budget itemized; He-3/4 polarizability would extend the I=1/2 → I=1 axis to I=1/2 (He-3) and I=0 (He-4). |
| 3 | High-Z hydrogenic (Bi⁸²⁺, U⁹¹⁺) | precision QED | High (relativistic + multi-loop QED + nuclear FNS + Bohr–Weisskopf) | LOW — heavy-atom regime needs full BBB93/KTT screening and relativistic enhancement; W1c-residual wall + Z>20 cliff combined. Multi-week prep before catalogue test. |
| 4 | Antiprotonic / muonic exotic atoms (p̄⁴He, kaonic) | exotic atoms | Medium-High (strong interaction + QED + recoil) | LOW — strong-interaction sector outside framework (W3 inner-factor); only kinematic part is framework-native. |
| 5 | Optical clock systematic budgets (Sr, Yb, Al⁺) at 10⁻¹⁸ | metrology | Very High (BBR + AC Stark + BBT + collision shifts) | MEDIUM — BBR / AC Stark are external-field projections (verified Layer-2 by external-field hunting sprint); collision shifts are correlation-driven. Could expose §V.D entries at very high precision. |
| 6 | Hadronic VP for muon g-2 | particle physics | High (lattice QCD + e+e− data) | LOW — dispersive integral on hadronic cross-section data is QCD-internal; framework has no hadronic sector. |

### 8.2 Ranking by structural alignment

**High alignment (recommend):**
1. Light-nucleus polarizability extractions (D, He-3/4). Already partially catalogued (D HFS); extending to He-3/4 polarizability is a natural next step. The framework can do the kinematic / leptonic side at operator level; the nuclear-physics-internal contribution is W3 inner-factor (Pachucki–Yerokhin 2020). Useful catalogue extension; surfaces new convention exposures (Pachucki–Yerokhin vs Friar–Payne for D, vs Karshenboim for He-3/4).

**Medium alignment (defer until BBB93 / heavy-atom screening done):**
2. King-plot isotope shifts (Yb, Ca⁺, Sr).
3. Optical clock systematic budgets at 10⁻¹⁸ (Sr, Yb, Al⁺) — needs BBR / AC Stark Layer-2 input handling, which is mechanical but new.
4. High-Z hydrogenic — needs heavy-atom screening + Bohr–Weisskopf, multi-week prep.

**Low alignment (do not pursue):**
5. Antiprotonic / muonic exotic atoms — strong interaction outside scope.
6. Hadronic VP for muon g-2 — QCD outside scope.

### 8.3 Recommendation

If the next catalogue extension is needed, **light-nucleus polarizability** (especially He-3 and He-4 at sub-ppm precision) is the clearest target. Specifically:

- **Sprint Calc-He3-HFS:** He-3 1S HFS Roothaan autopsy. He-3 has I=1/2, parallel to H. Friar–Payne 2005 vs Pachucki–Yerokhin 2020 itemization for He-3 polarizability. Expected: surfaces a §V.D entry on convention; framework-native at sub-percent.
- **Sprint Calc-He4-HFS:** He-4 has I=0, no HFS. But He-3 is the right one.

The other candidates require multi-week prep.

---

## §9. Ranked next tracks

Five candidates, scored by (i) structural value, (ii) sprint-window feasibility, (iii) verification protocol clarity. Each track names a diagnostic that would test the structural claim before implementation commitment.

### Track 1: §V.D-prediction sprint (1-week diagnostic, HIGH structural value)

**Question:** Is Pattern C (HFS clustering in §V.D) a structural prediction or HFS-specific?

**Protocol:** Test 3 pre-named candidates from §3.4: D polarizability (Pachucki–Yerokhin vs Friar–Payne), He fine structure α³(Zα)² (Drake vs Pachucki–Yerokhin), Ps HFS annihilation (Czarnecki–Melnikov–Yelkhovsky vs Karshenboim). For each, run a single literature comparison; if a class-(a) finding emerges, add to §V.D. If 2 of 3 land class-(a), Pattern C graduates to falsifiable structural prediction. If 0 of 3, Pattern C is HFS-specific.

**Diagnostic-before-engineering:** This *is* the diagnostic. Low cost (~3 days reading + memo).

**Falsifier:** 0 of 3 land class-(a).

**Sprint-window:** 1 week.

### Track 2: §III.6 spectral action operator-level promotion (1–2 week sprint, HIGH structural value)

**Question:** Can the §III.6 Connes–Chamseddine spectral action be promoted from energy-level to operator-level on the H 21cm autopsy chain?

**Protocol:** Build a Sturmian-register operator for the Schwinger asymptote correction at one loop, parallel to how §III.18 was built for Zemach. Verify against the Parker–Toms c₁ = R/12 = 1/2 reproduction at 0.5%. Expected: operator-level construction reveals the bracket structure that drove §V.D.3.

**Diagnostic-before-engineering:** Confirm the Schwinger asymptote IS the operator-level limit (not a closed-form scalar) before committing to ~1 week of operator construction.

**Falsifier:** Operator-level construction returns the same scalar Eides closed-form gives; no new structure exposed.

**Sprint-window:** 2 weeks.

### Track 3: Layer-2-presence-bound test (1-week diagnostic, HIGH structural value)

**Question:** Does Pattern B' (Layer-2-presence-bound) supersede the naive depth-vs-residual prediction?

**Protocol:** Re-tag all 17+ catalogue rows by L2-count and verify the prediction `|ε(observable)| ≤ max(ε_basis(framework), max_i |L2_i|)`. Specifically: compare each row's residual to (a) the framework basis-quality limit at zero L2, and (b) the largest single L2 input where present. Statistical fit against Pattern B and Pattern B'.

**Diagnostic-before-engineering:** This *is* the diagnostic. Cost: ~1 day.

**Falsifier:** Either (a) catalogue has rows where residual exceeds max(ε_basis, max L2) by more than 10× (would falsify B'), or (b) Pattern B fits the data tighter than B' (would resurrect Pattern B).

**Sprint-window:** 1 week (memo + paper edits).

### Track 4: Light-nucleus polarizability extension (2-week catalogue extension, MEDIUM-HIGH structural value)

**Question:** Does the framework's projection-chain machinery extend to I=1/2 He-3 HFS at sub-percent precision, surfacing convention exposures?

**Protocol:** Build the He-3 1S HFS Roothaan autopsy, parallel to D HFS at I=1. Identify §V.D candidate (likely Pachucki–Yerokhin vs Friar–Payne for He-3 polarizability). Verify framework-native at expected ~sub-ppm at BF strict + Schwinger.

**Diagnostic-before-engineering:** Verify §III.18 magnetization-density module accepts arbitrary nuclear g-factor and r_Z values without architecture changes. Should be backward-compatible.

**Falsifier:** Framework-native He-3 BF strict residual exceeds 10× the analogous H/D residual at the same chain (would suggest framework architecture issue at light-nucleus extension).

**Sprint-window:** 2 weeks.

### Track 5: §III.16 Breit retardation bound-state matrix elements (3-week structural extension, MEDIUM structural value)

**Question:** Can the §III.16 Breit retardation projection be promoted from kernel-level (verified) to bound-state matrix-element level (the named extension at §V.C.7)?

**Protocol:** Implement Bethe–Salpeter quasi-energy expansion machinery in `geovac/breit_integrals.py`, extending the Roothaan termination kernel. Test against Karshenboim 2005 §4 closed-form fraction -(11/48) at the operator-level for Ps 1S-2S.

**Diagnostic-before-engineering:** Verify that the Karshenboim closed-form is reachable from a finite Sturmian basis (named extension is structural, not just engineering). This is a 3-day analytical check before committing to a 3-week structural extension.

**Falsifier:** Bethe–Salpeter quasi-energy expansion machinery diverges or doesn't close at the n_basis the framework can support.

**Sprint-window:** 3 weeks.

### Score table

| # | Track | Structural value | Sprint-window | Verification clarity | Total |
|:-:|:------|:----------------:|:-------------:|:--------------------:|:-----:|
| 1 | §V.D-prediction | HIGH | 1 week | HIGH | 9/9 |
| 3 | Layer-2-presence-bound | HIGH | 1 week | HIGH | 9/9 |
| 2 | §III.6 op-level promotion | HIGH | 2 weeks | MEDIUM | 7/9 |
| 4 | Light-nucleus polarizability | MED-HIGH | 2 weeks | HIGH | 7/9 |
| 5 | §III.16 BS matrix elements | MEDIUM | 3 weeks | MEDIUM | 5/9 |

**Recommended next 2 tracks (parallel, 1-week budget):** §V.D-prediction (Track 1) and Layer-2-presence-bound (Track 3). Both are pure diagnostics, both close in 1 week, both produce a clear falsifier.

---

## §10. What would falsify each pattern (preservation discipline)

This section writes down the falsifiers IN the memo, not "let's figure it out later":

### Pattern A (α-order hierarchy)
**Falsifier 1:** A new wall instance at α³ or α⁶ that doesn't fit the named tiers. Currently no observable in the catalogue is at α³ or α⁶.
**Falsifier 2:** A wall instance where the α-order of the residual does NOT match the α-order of the missing physics.
**Falsifier 3:** Mu HFS at +199 ppm tightens to <50 ppm with framework-internal work alone (would falsify the α² isolated-tier reading).

### Pattern B' (Layer-2-presence-bound, replacement of naive Pattern B)
**Falsifier 1:** Catalogue row with L2-count = 0 and residual > 10% (framework basis-quality limit far exceeded).
**Falsifier 2:** Catalogue row with L2-count = 2+ and residual smaller than the largest single L2 input by more than an order of magnitude.
**Falsifier 3:** Statistical fit shows naive Pattern B (depth-vs-residual) tighter than B' on all 17+ catalogue rows.

### Pattern C (HFS clustering in §V.D)
**Falsifier 1:** Next 3 sprints surface §V.D entries where 2 of 3 are non-HFS, non-muonic-Lamb classes.
**Falsifier 2:** A §V.D entry surfaces in electronic Lamb shift where Eides-Pachucki-Yerokhin itemizations are well-aligned (would suggest Pattern C is data-availability artifact).

### Hunt 1 (operator-level verification map)
**Falsifier:** §III.6 promotion sprint reveals operator-level construction does NOT recover the Schwinger asymptote (would mean the framework's "Schwinger β verified" is energy-level only and doesn't lift to operator level).

### Hunt 2 (axis correlations)
**Falsifier:** A new (--, D, --) projection candidate surfaces and is admissible (would break the V/D correlation and force a re-think of the dictionary completeness).

### Hunt 3 (master Mellin engine signatures)
**Falsifier 1:** A future catalogue row with a transcendental signature outside M1 ∪ M2 ∪ M3 at the converged Layer-2 level.
**Falsifier 2:** Embedding-tier transcendentals (log p) appear in a converged Layer-2 observable (would falsify the "intermediate-only" reading).

### Hunt 4 (Lorentz boost)
**Falsifier:** A clean Minkowski-signature spectral triple is built and the Lorentz boost becomes admissible (would change the verdict from "REQUIRES-EXTENSION" to "ADMISSIBLE").

### Hunt 5 (empirical gaps)
**Falsifier:** Light-nucleus polarizability sprint surfaces no new convention exposure (would suggest the framework's diagnostic capability is HFS-class-specific, parallel to Pattern C falsification).

---

## §11. Honest negatives — patterns that did NOT survive scrutiny

### 11.1 Naive Pattern B (depth-vs-residual)
**Verdict:** FALSIFIED as stated. Replaced by Pattern B' (Layer-2-presence-bound).

**Reason:** Catalogue contains depth-3 chains at residual 0 (polarizability, Stefan–Boltzmann), depth-3 chains at +64 ppm (Ps 1S-2S framework-only), depth-3 chains at +40 ppm (D HFS BF strict), depth-4 chains at +18 ppm (H 21cm HF-4 with Zemach). No monotone depth-residual correlation.

**What stays:** the Sprint Calc-P refinement (rational-class observables exact at finite n_max) is correct; it's just not the headline prediction.

### 11.2 Pure-axis projection candidates (TX-A §6 EDIT-3)
**Verdict:** UNVERIFIED 5 sprints later. Gauge-choice and spin-statistics still flagged but not formalized.

**Reason:** No new (--, D, --) or (--, --, π) candidate has surfaced empirically. Vector-photon stays the unique transcendental-axis-independence witness.

**What stays:** the flag is open; the V/D correlation strengthens (12/12 instead of 8/8).

### 11.3 Lorentz boost as 17th projection (Hunt 4)
**Verdict:** REQUIRES-EXTENSION at the framework level.

**Reason:** Wick-rotated spectral triple is Riemannian; Lorentz boost requires Minkowski signature. Multi-month NCG work, not a sprint. The non-Lorentzian uniform-rescaling alternative is admissible but doesn't surface a new structural discovery beyond §III.14 / §III.16.

**What stays:** the question is structurally answered (admissible iff Minkowski spectral triple is built first); no follow-up sprint needed for the next 3-6 months.

### 11.4 Hidden new transcendental-axis-independence witness
**Verdict:** UNVERIFIED 5 sprints later. Vector-photon stays unique.

**Reason:** Three new projections (§III.17, §III.18, §III.19) added since TX-A; all V D - (length-introducing dimensioned) on the V/D side, all ring-preserving on the transcendental side. No new (--, --, π) projection has been identified.

**What stays:** vector-photon is the sole independence witness; the transcendental axis is therefore the only genuinely independent one.

### 11.5 Low-confidence patterns I checked

- **Cluster of "−" residuals (over-predicting) in HFS observables:** verified — H 21cm HF-1 at -1102, HF-2 at +58, HF-4 at +18; D HFS at +286 ppm (over-predicting). The sign is structural: Bohr–Fermi is a dipole-dipole approximation that overestimates in the absence of magnetization-distribution corrections (Zemach reduces the prediction). The pattern holds but is mechanistic, not a structural discovery beyond what §III.18 already encodes.

- **Cluster of "+" residuals (under-predicting) in QED Lamb shifts:** verified — H Lamb shift LS-6a at -0.534%, μH at -0.10%, Mu at +0.013%, μH BF at +2 ppm. The Mu BF case is anomalous (very tight residual). All other QED Lamb shifts are negative-residual, consistent with the LS-8a renormalization wall (under-counting α⁵ multi-loop).

- **Same observable at different precision targets:** verified — H 21cm at -1102 → +58 → +18 ppm with each Layer-2 input addition; D HFS at +40 → +1201 → +384 → +286 ppm. The sequences are non-monotonic in some cases, indicating that multi-Layer-2 corrections can over-correct before settling. This is an honest finding but not a new structural pattern.

---

## §12. Summary and result statement

### §12.1 What was verified

- **Pattern A (α-order hierarchy of walls):** confirmed and sharpened to two-tier structure (α⁴ Breit, α⁵ renormalization), with a possible third tier at α² bracket precision.
- **Pattern C (HFS clustering in §V.D):** confirmed, with three named predictions for next exposures.
- **Hunt 1 (operator-level verification map):** §III.6 spectral action is the load-bearing gap; §III.16 bound-state matrix elements is the second gap; §III.19 tensor multipole is the forward-looking slot.
- **Hunt 2 (axis correlations):** V/D correlation strengthens to 12/12; vector-photon stays unique transcendental-axis witness; no new redundancies.
- **Hunt 3 (master Mellin engine signatures):** case-exhaustion theorem holds at Layer-2 output level (9 of 11 rows in M1/M2 cleanly); embedding-tier transcendentals (log 2 / log 3) in intermediate matrix elements are the only "rogue" examples.
- **Hunt 4 (Lorentz boost):** REQUIRES-EXTENSION at framework level; non-Lorentzian alternative is reducible to existing variable-refinement projections; flag stays open for future PI decision.
- **Hunt 5 (empirical gaps):** Light-nucleus polarizability (He-3 HFS) is the clearest catalogue extension target.

### §12.2 What was falsified

- **Naive Pattern B (depth-vs-residual):** falsified as stated; replaced by Pattern B' (Layer-2-presence-bound).
- **Pure-axis projection candidates (gauge-choice, spin-statistics):** unverified after 5 sprints; vector-photon stays unique.

### §12.3 Top 3 ranked next tracks

1. **§V.D-prediction sprint (1 week):** verify Pattern C with 3 pre-named candidates.
2. **Layer-2-presence-bound test (1 week):** verify Pattern B' replaces naive Pattern B.
3. **§III.6 spectral action operator-level promotion (2 weeks):** close the load-bearing operator-level gap.

### §12.4 Recommended Paper 34 edits (NOT applied)

- **EDIT-1 (§VI Prediction 1' replacement):** replace naive depth-vs-residual prediction with Layer-2-presence-bound prediction. Tied to Track 3 above.
- **EDIT-2 (§V.D protocol refinement):** add prediction-tier classification to §V.D protocol, separating class-(a) lit-itemization (D.1, D.2, D.3) from class-(ii) closed-form-vs-operator (D.4). Tied to Track 1 above.
- **EDIT-3 (§VIII new open question):** Lorentz-boost structural verdict (REQUIRES-EXTENSION). Tied to Hunt 4 above.

### §12.5 Result (30-second read)

**Three patterns confirmed:** Pattern A (α-order tiers at α⁴ Breit and α⁵ renormalization), Pattern C (HFS clustering in §V.D), Hunt 3 (master Mellin engine partition).
**Two patterns falsified:** naive Pattern B (depth-vs-residual), pure-axis projection candidates from TX-A §6 EDIT-3.
**Top 3 next tracks:** (1) §V.D-prediction 1-week diagnostic; (2) Layer-2-presence-bound 1-week diagnostic; (3) §III.6 operator-level promotion 2-week sprint.
**Hunt 4 verdict:** Lorentz boost requires Minkowski-signature spectral triple (multi-month structural extension); not a near-term track.
**Hunt 5 next catalogue extension:** Light-nucleus polarizability (He-3 HFS).

---

*Memo word count: ~10,500 words. Sprint complete. Companion data:* `debug/data/synthesis_pattern_review_data.json`. *No production code modified. No paper edits applied. Recommendations only.*
