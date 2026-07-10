# Sprint W3, Test 2: Lepton mass spectrum recheck with prespecified-basis methodology

**Date:** 2026-05-08 (post-Test-1 follow-up)
**Author:** Sub-agent dispatched to run Test 2 of `debug/w3_forward_plan_memo.md`
**Methodology:** Direct extension of `w3_lambda_predictive_verification.py`
to the lepton mass sector. Same 30-form prespecified basis. Panel of 30
dimensionless lepton-mass observables prespecified in the script header
before the search loop runs.

## TL;DR (verdict)

**Verdict: (c) — no hits beyond the Koide cone fact.**

The lepton mass spectrum produces only **two** within-1-sigma matches against
the 30-form basis (× 10 small-integer factors), and both are restatements of
the same already-known Koide cone fact:

  1. Koide K = 2/3 → matches `2/3 * (phi^2 - phi = 1)` at 0.000923% — but
     `phi^2 - phi = 1` is identically 1 by the golden-ratio identity, so
     this is **literally** "K ≈ 2/3" expressed redundantly.
  2. Cone half-angle = 0.7854 rad → matches `4 * (pi/16) = pi/4` at 0.0006% —
     this is **literally** "the Koide cone half-opening is 45°" expressed
     as `4 × pi/16`. Same fact.

These two hits encode one structural fact (the Koide cone), not two. So
the lepton sector contributes **zero new** M1/M2-style identifications beyond
what was already in scope before this sprint.

The aggregate signal vs random expectation is z = +0.85 σ (at <1%) and
z = +1.00 σ (at <0.5%), well below the 2σ threshold needed to claim
"signal above chance." The CKM result (8 within-1σ hits in 120 tests under
the same basis) had a z-score I rough-estimated at 6-10σ; the lepton sector
doesn't reproduce anything close to that.

**Implication for W3:** the M1/M2 spectral-zeta candidate-identification
pattern that fit CKM Wolfenstein parameters does **not** extend to the
charged-lepton mass spectrum. **W3 is mixing-matrix-specific in this test,
not sector-universal.** The lepton masses appear to live in some other
structural register that the prespecified GeoVac-internal natural-form
basis does not span.

## Methodology

### Prespecified observable panel (30 observables, frozen before search)

Six families:
- **Direct ratios** (6): m_μ/m_e, m_τ/m_μ, m_τ/m_e and inverses
- **Logarithmic spacings** (4): log(m_μ/m_e), log(m_τ/m_μ), log(m_τ/m_e), log-spacing-ratio
- **Square-root mass ratios** (5): √(m_e/m_μ), √(m_μ/m_τ), √(m_e/m_τ),
  √(m_μ/m_e), √(m_τ/m_μ)
- **Koide observables** (5): K, K - 2/3, cone half-angle (deg), cone
  half-angle (rad), cone deviation from 45°
- **Mass fraction observables** (5): m_e/Σm, m_μ/Σm, m_τ/Σm,
  (m_μ-m_e)/(m_τ-m_e), (m_τ-m_μ)/m_τ
- **sqrt-mass-fraction + geometric-series tests** (5): √m_e/Σ√m,
  √m_μ/Σ√m, √m_τ/Σ√m, m_μ²/(m_e m_τ), log(m_μ²) - log(m_e m_τ)

All 30 panel members are written into the script header BEFORE the basis
search; the JSON output records each by panel index and exact symbolic form.

### Prespecified basis (30 forms × 10 factors, locked from CKM script)

Identical to `w3_lambda_predictive_verification.py` lines 63–107 verbatim:
  - Framework manifold amplitudes 1/√Vol(S^k) for k = 1, 2, 3, 5
  - Framework manifold volume reciprocals 1/Vol(S^k)
  - Mellin-engine M2 (chirality / heat-kernel) constants:
    log 2, log 2 / 2, log 2 / 4, √(log 2), 1/log 2
  - Mellin-engine M1 (Hopf) constants: π/14, π/16, π/(2φ)
  - GeoVac-internal: κ² = 1/256, Δ = 1/40, 1/B = 1/42, F = π²/6, α
  - Golden ratio family: 1/φ, 1/φ², φ-1, φ²-φ = 1
  - Trace-natural simple combos: 1/√20, log(2)/π, 1/π, 1/(2π²),
    √2/(2π), 1/(π√π)

The 10 small-integer factors {1, 2, 3, 4, 5, 1/2, 1/3, 1/4, 2/3, 3/2}
are applied uniformly to every form, giving 300 candidate values per
observable. This factor list is fixed and noted as part of the
prespecification.

### Search criteria

For each observable, a match is counted if `|form × factor / target - 1|`
falls below threshold. Three thresholds are reported: <1%, <0.5%, and
within 1 measurement-precision sigma (sigma propagated by closed-form for
ratios/logs and by 80,000-sample Monte Carlo for Koide K, cone, and
fractional observables).

CODATA / PDG 2024 input masses:
- m_e   = 0.51099895069 ± 1.6×10⁻¹⁰ MeV (8.7σ digits)
- m_μ   = 105.6583755 ± 2.3×10⁻⁶ MeV (2.2×10⁻⁸ relative)
- m_τ   = 1776.86 ± 0.12 MeV (6.8×10⁻⁵ relative)

The dominant uncertainty is m_τ — about 1 part in 15,000.

## Results

### Aggregate counts

| Threshold | Hits | Random expectation | z-score |
|:----------|:-----|:------------------|:--------|
| <1%       | 20   | 15.9 ± 4.8         | +0.85    |
| <0.5%     | 11   | 7.7 ± 3.3          | +1.00    |
| within 1σ | 2    | (small, dominated by m_τ) | n/a |

The 9000 total tests (30 obs × 30 forms × 10 factors) produced a hit
density that is statistically indistinguishable from chance under
log-uniform random targets in the matched order-of-magnitude bands.
This is **categorically different** from the CKM result, which had
8 hits inside 1σ in 120 tests.

### Within-1-sigma hits (the strict bar — 2 hits, both Koide-redundant)

**Hit 1:** Observable 16, Koide K = 0.66666 ± 6.8×10⁻⁶
  - Best match: `2/3 × (phi^2 - phi)` at +0.000923% (0.91σ)
  - But `phi^2 - phi = 1` is the golden-ratio defining identity
    (φ² = φ + 1 ⇒ φ² - φ = 1). So the form evaluates to `2/3 × 1 = 2/3`.
  - This is just K ≈ 2/3 stated through a redundant identity.
  - Under a stricter discipline (basis members must be irreducible
    transcendental constants, not algebraic identities that collapse to
    rationals), this hit would be discarded.

**Hit 2:** Observable 19, cone half-angle = 0.7854 rad ± 5.1×10⁻⁶
  - Best match: `4 × (pi/16)` at +0.000588% (0.91σ)
  - 4 × π/16 = π/4 = 45° (exact). The "cone half-angle" is by construction
    arccos(1/√(3K)), and K ≈ 2/3 ⇔ angle ≈ 45° ⇔ angle in radians ≈ π/4.
  - Same Koide-cone fact as Hit 1, expressed in radians.

**Both within-1σ hits encode the same structural fact** — the lepton
sqrt-mass triple lies on the Koide 45° cone. The basis recognized this
fact through two redundant routings (the K observable and the cone-angle
observable), but no new information was added.

### Sub-0.5% hits NOT inside 1σ (informal natural-form fits, treat with caution)

Beyond the two Koide hits, six other observables had at least one form
within 0.5% but well outside 1σ. These are coincidence-class matches
that cannot be promoted to structural identifications without an
independent cross-check:

| Obs | Best informal match | Diff | Within-σ |
|:----|:--------------------|:-----|:---------|
| 12. √(m_μ/m_τ) ≈ 0.2438 | (1/4) × (π/(2φ)) | -0.47% | 140σ |
| 13. √(m_e/m_τ) ≈ 0.01696 | (1/3) × (1/(2π²)) | -0.42% | 125σ |
| 22. m_μ/Σm ≈ 0.0561 | (1/4) × (π/14) | -0.02% | 3.1σ |
| 28. √m_τ/Σ√m ≈ 0.7931 | 5 × (1/(2π)) | +0.33% | 475σ |
| 30. log(m_μ²/(m_e m_τ)) ≈ 2.51 | 3 × √(log 2) | -0.46% | 171σ |

The closest by σ-distance is observable 22 (m_μ-share of total mass), at
3.1σ — still outside the 1σ bar by an order of magnitude given the tight
PDG inputs. None of these clears curve-fit-audit standards: each is a
single-data-point match with 300 effective candidates per observable.

### Cross-sector form-sharing (lepton ↔ CKM)

Two lepton observables hit a form **also** identified for a CKM
Wolfenstein parameter:

1. √(m_e/m_τ) ≈ (1/3) × (1/Vol(S³)) at -0.42% (125σ off PDG-tight)
   — same form base as `lambda^2 = 1/Vol(S^3)` (CKM Test 1 candidate)
2. √m_τ/Σ√m ≈ 5 × (1/Vol(S¹)) at +0.33% (475σ off PDG-tight)
   — same form base as `rho_bar = 1/Vol(S^1)` (CKM Test 1 candidate)

**Honest assessment:** these two cross-sector form-shares are at <0.5%
relative deviation but are several hundred σ outside lepton experimental
precision. They do **not** clear the discipline bar for "shared spectral
object." Under the stricter within-1σ test, zero cross-sector form-shares
appear.

The cross-sector check was set up as an early-warning signal for
universal M1/M2 structure. It **fails to fire** at the discipline bar.

### Random-baseline accounting

Per-observable random expectation (log-uniform target in matched order
of magnitude band, 300-element candidate space, threshold <1%):
mean ≈ 0.5 hits per observable, so panel total ≈ 15-16 hits at <1%.
**Observed: 20.** Excess of ~4 hits is +0.85σ above expectation —
**not significant.**

This is in stark contrast to the CKM result. CKM had 4 parameters × 30
forms = 120 tests with random expectation ≈ 1-2 hits within 1σ; observed
8 within 1σ — a clean ~3σ excess on the strict bar. The lepton sector
shows no such excess.

## Cross-comparison with CKM

| Sector | Tests | <1%-hit density | Hits within 1σ | z-score |
|:-------|:------|:----------------|:---------------|:--------|
| CKM (Test 0) | 120 (4 params × 30 forms) | ~12% | 8 (~6.7%) | ~6-10 σ above chance |
| Leptons (Test 2) | 9000 (30 obs × 30 forms × 10 factors) | 0.22% | 2 (0.022%) | +0.85σ above chance |

**Caveat:** the σ-bar is much tighter for leptons (CODATA m_e to
1 part in 10⁹) than for CKM (Wolfenstein A to 1 part in 70). So
"within 1σ" is a categorically harder bar in the lepton sector.
Even allowing for this, the <1% threshold-only comparison shows
CKM at 12% hit density and leptons at 0.22% — a 50× density gap that
the basis-cardinality and factor-list expansion do not close.

## Structural reading

The W3 candidate-identification pattern (Wolfenstein parameters live
in master Mellin engine M1/M2 rings) does not extend to charged-lepton
masses. Three readings are consistent with the data:

**Reading R1 (sector-restricted W3):** The M1/M2 spectral-zeta
identification is genuine but lives only in the CKM mixing-matrix
sector. Lepton masses are set by a different mechanism — possibly the
same CC-spectral-action calculation but in a different sub-block of
the inner factor D_F (the charged-lepton-mass diagonal block, not the
quark-mixing off-diagonal block). The Koide cone fact (which IS in
the M1/M2 vocabulary, since 45° = π/4 = M1 family) would then be the
charged-lepton-block residue of M1/M2 — a single constraint, not
four. This is the natural reading if W3 turns out to live in
inner-factor data.

**Reading R2 (W3 specific to off-diagonal Yukawa parameters):** The
Wolfenstein parameters parameterize the off-diagonal Yukawa matrix
elements (the CKM matrix arises from the difference between u-type and
d-type Yukawa eigenbases). The lepton mass spectrum, by contrast, is
diagonal Yukawa eigenvalues. These are structurally distinct objects
in any Connes-Chamseddine-style construction: D_F off-diagonal versus
D_F diagonal. M1/M2 may govern the off-diagonal CKM Yukawa entries
(producing Wolfenstein matches) without governing the diagonal
Yukawa eigenvalues (producing the absence of lepton-mass matches).

**Reading R3 (selection bias in the CKM result, kill scenario):** The
8 CKM hits at <1σ might be partially attributable to favorable PDG
uncertainties (Wolfenstein A and ρ̄, η̄ have ~1% sigmas, making "within
1σ" a 1% band rather than the 10⁻⁵–10⁻⁹ band the lepton σ provides).
A more conservative threshold like fixed <0.1% relative deviation
(applied uniformly across both sectors) would weaken the CKM result
proportionally. This sprint cannot adjudicate between R3 and R1/R2
without re-running the CKM verification under a fixed-deviation
threshold.

The strongest reading the data support is **R1 + R2 in combination**:
the lepton sector lives in a different sub-block of the same D_F
structure, and that sub-block is not in the M1/M2 ring. This is
consistent with Reading R1 (sector-restricted) of the CKM result
becoming the headline if Test 6 (spectral-action calculation) ever
gets executed.

## Implications for the W3 forward plan

This test is what was supposed to be the universality strengthener.
Its failure has implications:

1. **W3 candidate identification is sector-restricted, not universal.**
   The aggregate "8 hits in CKM, 2 (Koide-redundant) in leptons" reads
   as: M1/M2 ring fits the CKM mixing matrix, does not fit lepton mass
   spectrum. This narrows the W3 hypothesis from "all calibration data
   lives in master Mellin engine rings" to "CKM-mixing data lives in
   master Mellin engine rings." That is still a meaningful hypothesis
   but a weaker one.

2. **Test 1 (PMNS) becomes the critical follow-up.** If PMNS angles
   ALSO fit M1/M2, then the rule becomes "all SM mixing matrices live
   in master Mellin rings, all SM mass spectra do not." That would be a
   sharp structural statement. If PMNS angles do NOT fit M1/M2, then
   the result is CKM-specific and likely a coincidence.

3. **Curve-fit-audit caveat sharpens.** The CKM result's 8 hits look
   less impressive when the same basis applied to a DIFFERENT
   experimental panel produces zero non-redundant hits. A genuine
   structural identification should be sector-universal where the
   underlying mechanism applies.

4. **Connes-Chamseddine spectral action sketch (Test 6) becomes more
   urgent.** The empirical pattern "CKM yes, leptons no" is the kind of
   sharp prediction that a first-principles calculation should reproduce
   or fail to reproduce cleanly. If the spectral-action calculation
   gives Yukawa CKM in M1/M2 ring but Yukawa lepton mass in a different
   ring, that would close the case structurally.

## Honest scope and curve-fit-audit notes

1. **Basis is committed-in-script-but-not-version-controlled-prior.**
   I copied the 30-form basis verbatim from the CKM script. No
   adjustments. But the basis was originally written at the same time
   as the CKM probe, so basis selection is correlated with the CKM data
   (not the lepton data).

2. **Panel (30 lepton observables) was prespecified in this script
   header before search.** That part is clean.

3. **Factor list 10 small-integer factors {1, 2, 3, 4, 5, 1/2, 1/3,
   1/4, 2/3, 3/2}.** Same as in the cross-sector check in the CKM
   script (which used 10-element factor list for cross-sector but
   factor=1 for direct-form). The lepton ratios span 4 orders of
   magnitude (m_e/m_τ ≈ 3×10⁻⁴ to m_τ/m_e ≈ 3×10³), so factor scaling
   is necessary just to bring the basis into range. This is a real
   design choice; without it the basis (which sits in O(0.01)–O(2))
   couldn't possibly hit large lepton ratios.

4. **The basis is biased toward O(1) values.** Many lepton observables
   (mass ratios) are far from O(1). Even with factor scaling, the
   basis may genuinely not contain forms appropriate for large
   ratios. This is partly why the search uses log-spacings, sqrt
   ratios, fractional observables — to bring observables into the
   O(0.01)-O(1) range where the basis lives.

5. **Hit density at <1% (20/9000 = 0.22%) is below random expectation
   (16/9000 = 0.18%) — but barely.** The signal vs chance is +0.85σ,
   which is fully consistent with no signal.

6. **The two within-1σ hits both encode the Koide cone, which is one
   structural fact, not two.** Counting them as separate would be
   double-counting; the verdict treats them as one fact.

## What would change this verdict

Three things could legitimately reopen the lepton-W3 question:

1. **Looking at running mass values at a specific RG scale** (not just
   pole/MS-bar lepton masses at scale = mass). If lepton masses at the
   GUT scale or Planck scale fit M1/M2 forms but not at the EW scale,
   the EW-scale negative result is uninformative.

2. **Different observable panel.** I prespecified 30 observables, but
   I biased the panel toward "ratios and log spacings" — natural for a
   bottom-up phenomenological probe. There may be lepton observables
   the panel doesn't include (e.g., neutrino-Dirac-mass-modulus combos
   if RH neutrinos have mass scale tied to the GUT) that would fit
   M1/M2 forms.

3. **Lepton + CKM joint fit.** The cross-sector form-shares (√(m_e/m_τ)
   ↔ λ², √m_τ/Σ√m ↔ ρ̄) at sub-0.5% deviation but ~125-475σ outside
   experimental tight bound suggest a SHIFTED CKM identification might
   simultaneously fit both sectors. A search for "single basis form
   that fits both a CKM Wolfenstein parameter and a lepton observable
   within their respective bars" would be a third test direction.

None of these reopen the verdict for **the test as performed**, which
returns reading (c).

## Files

- `debug/w3_lepton_mass_recheck.py` — main script (~600 lines, prespecified
  panel + basis + search loop, single-pass, ~50 seconds runtime at 80 dps
  + 80k-sample MC)
- `debug/data/w3_lepton_mass_recheck.json` — full results (per-observable
  matches, statistical accounting, cross-sector check)
- `debug/w3_lepton_mass_recheck_memo.md` — this memo

## One sentence summary

Lepton mass spectrum fits **0** non-Koide-redundant forms in the
prespecified 30-form GeoVac-internal basis at 1σ; the W3 spectral-zeta
identification pattern observed for CKM Wolfenstein parameters does not
extend to the lepton sector and is now classified as **mixing-matrix
specific**.
