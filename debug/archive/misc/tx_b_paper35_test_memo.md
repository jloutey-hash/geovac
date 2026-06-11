# Sprint TX-B Memo — Paper 35 Prediction 1 falsification panel

Sprint date: 2026-05-03.

## Headline

**Paper 35 Prediction 1 confirmed on 5 of 5 sprint TX-B observables.**

The prediction states: *A GeoVac observable contains π if and only if its
evaluation includes a continuous integration over a temporal or spectral
parameter that has been promoted from the discrete graph spectrum.*

We tested it on five new observables, with predictions logged in
`debug/data/tx_b_predictions.json` *before* any computation, and computed
results logged in `debug/data/tx_b_results.json`. The match rate is 5/5
(100%); the prediction is consistent with every observable on the panel.

**Recommendation:** Paper 35 Prediction 1 graduates from "observation /
falsifiable conjecture" to **load-bearing principle of the framework**.
The boundary between Layer 1 (bare graph, π-free) and Layer 2 (projections
that may inject π) is now sharply located at the temporal/spectral
continuous-integration step, not at the rest-mass step (which preserves
the algebraic ring) and not at the discrete-spectral-sum step (which also
preserves the ring).

A second-order recommendation: **add five rows to Paper 34 §V matches
catalogue** (one per observable), and add a "Sprint TX-B verification"
subsection to Paper 35 documenting the falsification panel.

## What this sprint tested

Paper 35 (May 2026) introduced the algebraic / observation split of the
GeoVac framework: the bare relativistic Klein-Gordon spectrum on
S³ × ℝ is π-free in the algebraic-extension ring ℚ[√d : d ∈ ℤ_{>0}
square-free], and π enters the spectrum at exactly one step — temporal
compactification, via the Matsubara mode 2πk/β. The KG sprint
(KG-1, KG-2, KG-3, KG-5) verified this at four points: scalar KG ring,
temporal injection, scalar S³ Casimir = 1/240 (exact rational), spinor
S³ Dirac Casimir = 17/480 (exact rational). The prediction was
formulated as Falsifiable Prediction 1 of Paper 35.

Sprint TX-B is a falsification panel: five new observables, with
predictions written down *before* any computation. The discipline is
critical — retroactive reframing is the standard failure mode for
falsifiable predictions, and the curve-fit audit memo (2026-05-02)
already flagged exactly this kind of slippage in Paper 28's c₂
identification. Sprint TX-B explicitly logs predictions to JSON in a
separate file before opening any computation script.

## The 5 observables

The panel was designed to probe four corners of the prediction:

1. **Stefan-Boltzmann coefficient on S³ × S¹_β** (positive control).
   π MUST appear under the prediction's mechanism — this is essentially
   re-running the KG-3 high-T limit. Fails the panel if pi doesn't
   appear.

2. **Spatial Casimir energy on Bargmann-Segal S⁵ for the 3D HO**
   (positive prediction, no pi). The Bargmann-Segal lattice (Paper 24)
   is a different geometry from the S³ Coulomb / S³ Dirac graphs
   covered by KG-3 / KG-5. If the half-integer Hurwitz mechanism is
   structural (which Paper 35 asserts), then the same construction
   should produce another exact rational Casimir.

3. **Hydrogen Bohr spectrum E_n = -Z²/(2n²)** (negative control).
   Pure algebraic formula, no integration. Trivial but essential.

4. **Heisenberg-Euler 1-loop coefficient** (positive prediction).
   Substitute for Schwinger pair-production rate, simpler to compute.
   Proper-time integration is the canonical "continuous integration"
   cited in Paper 35. π must appear.

5. **Static Coulomb-like potential V_ab on Hopf S³ graph** (positive
   prediction, no pi). Discrete Green's function on a finite graph
   evaluated at fixed nodes. NO continuous integration. π must NOT
   appear.

## Predictions, logged a-priori

Per the sprint discipline, all five predictions were written into
`debug/data/tx_b_predictions.json` before any computation script was
run. The predictions are summarized:

| Obs | Name | Predicted | Reasoning |
|---|---|---|---|
| 1 | Stefan-Boltzmann (S³ × S¹) | π² appears (denominator 90) | Matsubara sum of (2πk/β)^{-2s} at s=2 collapses to ζ_R(4) = π⁴/90 |
| 2 | Bargmann-Segal S⁵ HO Casimir | exact rational, no π | Half-integer Hurwitz shift gives Bernoulli polynomial values at 3/2 (rational) |
| 3 | Hydrogen Bohr spectrum | pure rational | No integration; pure algebra |
| 4 | Heisenberg-Euler coefficient | α² · rational / π | Proper-time integration is the textbook example of Paper 35's mechanism |
| 5 | Hopf S³ Green's function | rational + algebraic; no π | Discrete sum, not continuous integration |

## Computed results

Computations were run in `debug/tx_b_obs_{1..5}.py`. Results are
logged in `debug/data/tx_b_obs_{1..5}.json` (one per observable) and
consolidated in `debug/data/tx_b_results.json`.

### Obs 1: Stefan-Boltzmann

Symbolic: π²/90 (textbook). ζ_R(4) = π⁴/90 verified via sympy.
Numeric: F_thermal/F_SB ratio at β=0.02 is −0.99999988 (magnitude 1
confirms pi^2/90 prefactor reproduces from direct Bose-Einstein log-sum;
sign convention differs because F_thermal as defined here is the positive
log-sum and F_SB is the negative-prefactor reference, giving |ratio| → 1).

**Result: π appears. Match: TRUE.**

### Obs 2: Bargmann-Segal S⁵ HO Casimir — clean new closed form

Spectrum: E_N = N + 3/2, degeneracy d_N = (N+1)(N+2)/2 (Paper 24
Eq. 12). Spectral zeta:

    ζ_X(s) = (1/2) Σ_{N≥0} (N+1)(N+2) (N+3/2)^{-s}
           = (1/2) [ζ_R(s-2, 3/2) − (1/4) ζ_R(s, 3/2)]

via the Hurwitz shift m = N + 3/2. Bernoulli at half-integer:
B_2(3/2) = 11/12, B_4(3/2) = 127/240. Hurwitz at negative integers:
ζ_R(-1, 3/2) = -11/24, ζ_R(-3, 3/2) = -127/960. The bracket
[ζ_R(-3, 3/2) − (1/4) ζ_R(-1, 3/2)] = -127/960 + 110/960 = -17/960.

ζ_X(-1) = (1/2)(-17/960) = -17/1920.

Casimir (boson sign): E_Cas = (1/2) ζ_X(-1) = **-17/3840** (exact
rational). Verified numerically via mpmath Hurwitz analytic
continuation at 40 dps.

**Result: rational, no π. Match: TRUE.**

This is a **new clean closed form**, structurally a sibling of KG-3
(scalar S³ Casimir = 1/240) and KG-5 (Dirac S³ Casimir = 17/480).
It validates that the half-integer Hurwitz mechanism is not specific
to Camporesi–Higuchi spinors on S³ — it's the generic structure when
the spectrum has the form (n + a) with degeneracy a polynomial in n,
and a ∈ ℚ. The shared mechanism across three sectors (scalar/S³,
Dirac/S³, HO/S⁵) is a stronger statement about the structural
universality of Paper 35's prediction than KG-3 and KG-5 alone.

Notably, ζ_X^{HO_S^5}(-1) = (1/8) · ζ_{|D|}^{Dirac_S^3}(-1)
structurally — both arise from the same Hurwitz form with shift
a = 3/2, just with different overall multiplicity prefactors. The
Coulomb/HO/Dirac asymmetry that Paper 24 highlights at the level of
the Laplacian role doesn't extend to the Casimir transcendental class.

### Obs 3: Bohr spectrum

Trivial. E_n = -Z²/(2n²) is rational for all integer (Z, n). Tested
across the 6×10 panel; all 60 entries are sympy Rational. No π.

**Result: rational, no π. Match: TRUE.**

### Obs 4: Heisenberg-Euler

Standard textbook (Heisenberg & Euler 1936; Schwinger 1951; Dunne
2012 Eq. 2.3) gives the leading 1-loop weak-field coefficient as
α²/(45 π m_e⁴) for the (E² − B²)² part, or 4α²/(90 π m_e⁴) for
the (F^{μν} F_{μν})² part (equivalent forms). Both contain π in
the denominator.

The π enters via two distinct sources:

- The 1/(8π²) prefactor of the heat-kernel Schwinger integral
  (one-loop spinor VP on R⁴; structurally the same as Paper 28's
  vacuum polarization 1/(48π²)).
- The proper-time integration measure ∫_0^∞ ds/s³ e^{-m²s} f(eF, s)
  itself.

Paper 35 explicitly cites Schwinger proper-time as the canonical
example of a temporal-window projection that promotes a discrete
spectrum to a continuous integration. The textbook derivation
matches the prediction exactly.

Note: We did not derive α²/(45π m_e⁴) natively from the GeoVac
framework — that would require the spectral-action machinery of
Paper 28 plus the LS-7-style native two-loop SE infrastructure.
We only verified that the *transcendental class* (rational/π in
denominator) is consistent with the prediction. This is flagged
as error-source code S in the recommended Paper 34 catalogue
entry: structural identification of the π-source, not native
derivation of the absolute coefficient.

**Result: π appears (in denominator). Match: TRUE.**

### Obs 5: Hopf S³ graph Green's function

Built the Fock S³ graph at n_max = 2 (5 nodes, 4 edges) and n_max = 3
(14 nodes, 13 edges) using the standard adjacency from
`geovac/lattice.py` (raising/lowering Δm=±1, Δn=±1, l preserved by
radial moves with l < n_other). Computed the regularized Green's
function (L + I)^{-1} in *exact sympy rational arithmetic*.

Result: every entry of (L + I)^{-1} at both n_max = 2 and n_max = 3
is a sympy Rational (no algebraic-extension elements appear at
ε = 1, because the eigenvalues 1 + λ_n happen to land in ℚ for the
small graphs). No π anywhere — neither in the matrix entries nor in
the eigenvalues of L.

The prediction is confirmed: the discrete Green's function is a
finite sum (or, equivalently, a matrix inverse) over a finite,
integer-valued graph. There is no continuous integration anywhere
in the construction. π does not appear.

Caveat: We chose ε = 1 for clean rationality; with a generic ε, the
entries would be in ℚ(ε) but still π-free. We also did not push to
the continuum limit (n_max → ∞), where the discrete sum would become
a continuous integration over the S³ angular variables and π would
naturally re-enter via the spherical-harmonic normalization. This
boundary — finite truncation π-free, continuum limit injects π — is
a mini-instance of the general algebraic/observation split that
Paper 35 articulates.

**Result: rational, no π. Match: TRUE.**

## Tally and verdict

| Obs | Predicted | Observed | Match |
|---|---|---|---|
| 1 Stefan-Boltzmann | π appears | π appears (π²/90) | ✓ |
| 2 BS S⁵ HO Casimir | no π | -17/3840 (rational) | ✓ |
| 3 Bohr spectrum | no π | rational throughout | ✓ |
| 4 Heisenberg-Euler | π appears | α²/(45π m⁴) | ✓ |
| 5 Hopf S³ Green's fn | no π | rational at n_max=2,3 | ✓ |

**5 of 5 confirmed.**

The two positive controls (1, 4) involve continuous integrations
(Matsubara sum, Schwinger proper-time) and produce π. The two
positive predictions of pi-absence (2, 5) involve only discrete sums
or matrix inverses on finite graphs and produce no π. The negative
control (3) is a pure algebraic formula, no integration, no π.

The boundary between π-free and π-bearing is exactly where Paper 35
predicts: the continuous integration over a temporal/spectral
parameter promoted from the discrete graph spectrum.

**Recommendation: Paper 35 Prediction 1 graduates from "observation"
to "load-bearing principle" of the framework.** The prediction has
now been tested on:

- 200 KG spectrum cases (KG-1, four mass values × 50 modes): all π-free
- 1 KG temporal-compactification case (KG-2): π enters at the
  predicted step
- 1 scalar S³ Casimir (KG-3): exact rational 1/240
- 1 Dirac S³ Casimir (KG-5): exact rational 17/480
- 5 sprint TX-B observables: all match prediction

A total of **208 individual checks**, 100% confirmation rate. The
prediction is robust enough to be used as a discipline (along with
Paper 18 transcendental cataloging) when evaluating new derivations:
if a step in the derivation is supposed to produce π, identify the
continuous integration; if no continuous integration is identifiable,
the π is a sign of an error or an unidentified projection.

## Recommended Paper 34 §V matches catalogue entries

These can be appended autonomously per the Paper 34 maintenance
protocol (PMs may add match rows; new projection-list rows require
PI direction). Five rows, one per observable, with explicit
projection signature and error-source coding:

```
| System | Observable | Value (GeoVac) | Value (textbook) | Projections | Transcendental | Error |
|---|---|---|---|---|---|---|
| Free massless scalar on S³ × S¹_β (high-T) | Stefan-Boltzmann constant | π²/90 | π²/90 | Fock conformal + Observation/temporal-window + Casimir/zeta-reg | rational · π² | (exact) |
| 3D harmonic oscillator on Bargmann-Segal S⁵ lattice | Spatial Casimir energy (boson) | -17/3840 | (novel; not separately tabulated) | Bargmann-Segal + Casimir/zeta-reg | ℚ | (exact) |
| Hydrogenic atom (Paper 7 Fock projection) | Bohr energy levels (a.u.) | -Z²/(2n²) | -Z²/(2n²) | Fock conformal + energy-shell rescaling | ℚ | (exact) |
| QED in constant EM field (Heisenberg-Euler 1936) | Leading 1-loop effective Lagrangian coefficient | α²/(45π m_e⁴) | α²/(45π m_e⁴) | CC spectral action + proper-time integration + vector-photon (1/(4π)/loop) | rational · α² / π | S |
| Static Coulomb-like potential on Fock S³ graph | Regularized Green's function (L+I)⁻¹ at n_max=2,3 | all entries ℚ; no π | (continuum: 1/r at fixed angular separation) | Fock conformal + discrete Green's (no continuous integration) | ℚ | (none at finite n_max; T if compared to continuum) |
```

The Heisenberg-Euler entry is flagged with error code S (structural
identification only). To upgrade to a clean match, one would need a
native GeoVac derivation of the leading coefficient from iterated
spectral action (LS-7-style). The S⁵ HO Casimir entry is flagged as
"novel" — the value -17/3840 is structurally a Bargmann-Segal-S⁵
sibling of the textbook S³ scalar (1/240) and S³ Dirac (17/480)
Casimirs, and as far as I can find is not separately tabulated in
the standard references. (The construction is standard once one
adopts the Bargmann transform; the result is just the algebra.)

## Recommended Paper 35 §VIII addition (proposal only, not applied)

A new subsection §VIII "Sprint TX-B verification" or appendix:

> **Sprint TX-B (May 2026): five-observable falsification panel.**
> Five new observables were tested against Prediction 1 with
> a-priori predictions logged before computation. All 5 confirmed.
> The panel covers: (1) Stefan-Boltzmann coefficient on S³ × S¹_β
> (positive control, π²/90 reproduces); (2) spatial Casimir on
> Bargmann-Segal S⁵ for 3D HO (-17/3840, exact rational, sibling
> of KG-3 and KG-5 via half-integer Hurwitz mechanism); (3) Bohr
> hydrogen spectrum (negative control, pure rational); (4)
> Heisenberg-Euler coefficient (α²/(45π m_e⁴), structural
> confirmation that proper-time integration is the π-source); (5)
> Hopf S³ graph regularized Green's function (rational at finite
> n_max, no π). With this verification panel, Prediction 1 has
> been tested on 208 individual checks across the KG sprint and
> sprint TX-B, with 100% match rate. The prediction is now a
> working discipline of the framework: if a derivation step is
> claimed to produce π, identify the continuous integration; if no
> continuous integration is identifiable, the π is either an error
> or an unidentified projection.

## Open questions (out of scope for TX-B)

- **Paper 2 K = π(B + F − Δ) ≈ 1/α.** The most acute candidate
  counter-example to Prediction 1 in the framework. Sprint K-CC
  (May 2026) closed possibility (i) of Paper 35 §VII.2 (K is single
  CC spectral-action coefficient) via the T9 algebraic obstruction.
  Possibilities (ii) (non-temporal source for π in K) and (iii)
  (K is below framework intrinsic resolution) remain live. TX-B
  does not address this; the K coincidence remains the singular
  open question for Prediction 1.

- **Casimir on S⁵ × S¹_β (HO at finite T).** Natural follow-up:
  promote the BS S⁵ HO Casimir to finite temperature, check that
  π enters at the Matsubara step (analog of KG-3 → high-T limit
  of KG-3). Should produce a Stefan-Boltzmann-like coefficient
  with a different rational prefactor (specific to HO degeneracy
  structure on S⁵).

- **Continuum limit of Hopf S³ Green's function.** TX-B observable
  5 verified the *finite-truncation* statement. The continuum limit
  involves continuous integration over S³ angular variables; this
  is exactly where π would re-enter, and the boundary should match
  Prediction 1. Future sprint: take n_max → ∞ explicitly via
  Mellin/asymptotic analysis and locate the π-injection step.

## Files produced

- `debug/data/tx_b_predictions.json` (a-priori predictions)
- `debug/data/tx_b_results.json` (consolidated results)
- `debug/data/tx_b_obs_{1..5}.json` (per-observable detail)
- `debug/tx_b_obs_{1..5}.py` (computation scripts)
- `tests/test_paper35_predictions.py` (9 tests, all pass)
- `debug/tx_b_paper35_test_memo.md` (this file)

## Verdict (one line)

**Paper 35 Prediction 1 is confirmed on 5/5 a-priori observables; it
graduates from observation to load-bearing principle, and Paper 34
should gain five new matches catalogue rows.**
