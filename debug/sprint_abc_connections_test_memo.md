# Sprint: testing the three /ahha connections on the chemistry geometry defect
**Date:** 2026-07-06
**Mode:** diagnostic (MODEL the residual, do not reduce it). Follow-on to v4.71.0
(chemistry-error projection decomposition). PI: "I want to test all 3 of these ideas."
**Connections under test (from the /ahha pass on v4.71.0):**
- **A** — the geometry defect lives in the *non-compact R-coordinate* (free-side /
  multi-focal-composition wall; same seat as WH7's time & the calibration continuum).
  *Falsifier:* a zero-parameter skeleton-side form for the tilt built from discrete
  invariants (Gaunt/3j, integer eigenvalues, rational couplings).
- **B** (adjudicator) — the residual is either π-free/discrete (Layer-1, derivable) or
  continuum/radial (Layer-2). *Test:* which truncation knob carries the tilt.
- **C** — "energy converges but geometry drifts" is the textbook Pulay / gradient-lags-
  energy basis-response artifact. *Falsifier:* decouple the basis from R ⇒ tilt heals.

## Method
Reference-free anchor (no external PES; pyscf unavailable on Windows). At the *true*
minimum R_true=3.015, dE_true/dR=0, and the balanced total carries the **exact** physical
V_NN: the `nuclear_repulsion` field = −7.2799 (R-independent core double-count, derivative-
immune) **+ 3/R** (verified: field derivative = −3/R² to 1e-7). Therefore:
- computed tilt ≡ dE_c/dR(R_true) = ε′(R_true) exactly (the pure geometry error), and
- for zero net tilt the electronic gradient must equal −dV_NN/dR = **+3/R² = +0.3300**.

Knob battery (balanced LiH, n_max=2, live solver, 5-pt stencil h=0.05, cubic):
`n_grid_vne` ∈ {2000,4000,8000,16000} (radial quadrature), `L_max` ∈ {2,3,4,6,8}
(cross-V_ne multipole/angular order). Composed side: refit the tilt at R_true from the
cached l_dependent-PK PES (l=2,3,4) — no live run (Paper 17 already has l=2..7).

Drivers: `debug/sprint_abc_tilt_sensitivity.py` (+ `.log`, `abc_tilt_sensitivity.json`),
`debug/sprint_abc_composed_driftlaw.py` (+ `abc_composed_driftlaw.json`).

## Results

### Balanced knob battery (14 live configs, ~68 s each)
| knob | setting | tilt (Ha/bohr) | curv (Ha/bohr²) | curv/k_true |
|------|--------:|---------------:|----------------:|------------:|
| base | n2/g8000/L4 | −0.02978842 | 0.14793567 | 2.245× |
| n_grid_vne | 2000 / 4000 / 16000 | −0.02978842 (all) | 0.14793567 (all) | 2.245× |
| L_max | 2 / 3 / 6 / 8 | −0.02978842 (all) | 0.14793567 (all) | 2.245× |

**Bit-identical across every knob** (tilt spread 1.7e-12, curv spread 1.1e-11 = float
roundoff). The electronic gradient at R_true = **+0.30024** vs required **+0.33002** →
**9.03 % too weak**, and this deficit is knob-invariant.

Interpretation of the invariance (NOT trivial): both knobs feed the *cross-center V_ne* —
the R-dependent binding term that dominates the tilt. Their invariance means that term is
**exact** (angular multipole Gaunt-terminates at L=2 for l≤1 orbitals — the angular
sparsity theorem) and **quadrature-converged** (n_grid≥2000). The cross-attraction is
computed exactly *given the orbitals*; the 9 % gradient deficit is therefore **entirely in
the max_n=2 one-particle orbital basis**, not in the integral evaluation.

### Composed drift law (cached l_dependent PK, refit at R_true)
| l_max | n_ch | tilt | R_eq err | curv/k |
|------:|-----:|-----:|---------:|-------:|
| 2 | 9 | −0.0167 | +5.3 % | 1.51× |
| 3 | 16 | −0.0382 | +15.7 % | 1.45× |
| 4 | 25 | −0.0619 | +25.5 % | 1.56× |

Tilt **linear, diverging**: −0.0226/l_max, R²=0.9992 (`tilt ≈ −0.0226·l_max + 0.029`;
3-pt consistency, NOT a derivation — `tilt~l_max` and `tilt~n_ch` indistinguishable at 3
pts, audit rule). R_eq **monotonically worse**. Refining the *discrete angular* knob
*destabilizes* the R-geometry (Paper 17: continues to l_max=7, +0.168 bohr/l_max).

## Verdicts on the three connections

**B — ADJUDICATED → Layer-2 (continuum/radial), NOT Layer-1 (discrete skeleton).**
The defect has *zero* discrete-angular content (L_max 2→8 bit-invariant; angular
structure exact by Gaunt termination) and *zero* quadrature content. It lives entirely in
the continuum orbital-basis (max_n) sector. This **rules out** the "derive the tilt from
graph/Gaunt invariants (Paper-12-style algebraic closure)" reading. **My /ahha Step-2
de-hedged prediction — "the coefficients come out rational × graph-invariants (Layer-1)" —
is FALSIFIED.** The defect is not a discrete-graph object at all. (Protocol working as
intended — the bold leap was made, then killed by the computation.)

**A — SURVIVES its falsifier; supported.** The falsifier was "a zero-parameter skeleton-
side form built from discrete invariants." B shows there is *no* discrete content to build
from — the defect is 100 % continuum-radial. Composed confirms the positive signature:
discrete angular refinement (l_max) drives the non-compact R-geometry monotonically worse.
Unified reading across both architectures: **discrete-skeleton refinement is orthogonal
(balanced: bit-invariant) to — or anti-correlated with (composed: divergent) — the
non-compact R-coordinate accuracy.** The skeleton cannot pin R. This is the same free-side
structural position as WH7's time and the calibration continuum.

**C — SUPPORTED in its honest form; strong falsifier un-run.** The defect is basis-
completeness-limited (max_n), and the well-SHAPE lags the energy in the extreme: n_max
2→3 the **energy heals 8.8×** (1.75%→0.20%) while the **curvature is frozen** (0.1479→
0.1463, still 2.22× k; ω_e stays ~+45 % at its own minimum) and the tilt does not improve
(−0.030→−0.036). That is exactly the Pulay-class gradient/curvature-lags-energy basis
asymmetry (energy 2nd-order, gradient 1st-order in basis incompleteness). The strong
falsifier — decouple the orbital basis from the nuclei ⇒ tilt heals — is architecturally
hard and **was not run**; and the v4.71.0 overlap-law falsification already argued against
a *naive single-overlap* Pulay term (ε′ linear, not exponential), so if Pulay it is the
full density-weighted sum, not a single dS/dR.

## The residual A-vs-C question (the honest open item)
A and C are **two descriptions of one established fact**: the geometry defect lives in the
continuum orbital basis and does not heal at the energy's rate. The distinguishing
question — *irreducible free-side wall (A)* vs *slow-but-eventual basis-response error (C)*
— turns on whether the curvature converges as max_n→∞. Over n_max 2→3 it does **not**
converge (frozen), consistent with **both**. Deciding needs n_max≥4 (n_max=3 already
~2.3 h/pt; n_max=4 out of reach). This is a genuine §1.5 dual-description situation:
the framework currently cannot distinguish the "structural" and "convergent" readings, so
papers stay dual. Named decider (banked, not run): balanced LiH curvature at n_max=4.

## What each test bought
- Localized the balanced geometry defect to a single sector — **the max_n=2 orbital basis**
  — by *exclusion* (quadrature and angular multipole both bit-irrelevant). Sharper than
  v4.71.0, which had not separated the knobs.
- Reference-free "**electronic gradient 9.0 % too weak at R_true**" quantifier, knob-invariant.
- Composed: the tilt-∝-l_max drift is the *angular-refinement-destabilizes-geometry* free-
  side signature, quantified (−0.0226/l_max).
- Adjudicated B against my own bold prediction (Layer-1 → **falsified**, it is Layer-2).

## Honest scope / caveats
- Theorem grade: none. Diagnostic localization + one falsified prediction.
- Single system (LiH). Balanced universality (is the 9 %/+45 % generic across hydrides?)
  untested — needs BeH₂ balanced (~600k-det, hours/pt).
- max_n leg rests on cached n_max=3 (coarse 5-pt tilt fit; curvature robust). n_max≥4 out
  of reach ⇒ A-vs-C stays open by construction.
- Composed (0,0)-weight-dilution *derivation* of −0.0226/l_max banked (needs solver
  instrumentation to emit channel weights; cached JSON has PES only).
- No production (`geovac/`) code touched ⇒ no `/regression` required.

## Deliverables
- `debug/sprint_abc_tilt_sensitivity.py` / `.log` / `debug/data/abc_tilt_sensitivity.json`
- `debug/sprint_abc_composed_driftlaw.py` / `debug/data/abc_composed_driftlaw.json`
- this memo.
