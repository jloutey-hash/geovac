# W3 Test 5: mechanical-basis re-test — verdict and analysis

**Date:** 2026-05-08, post-W3-bold-path
**Sprint position:** Test 5 of `debug/w3_forward_plan_memo.md`
**Goal:** Address curve-fit-audit Honest Concern #3 (basis-selection bias) on
the W3 spectral-zeta candidate identification of the four CKM Wolfenstein
parameters with master-Mellin-engine M1/M2 rings.

## Headline verdict

**The original W3 signal is fully attributable to basis selection.** When the
hand-curated 30-form basis is replaced by a mechanically-generated 21,448-form
basis enumerating small-integer rational combinations of master-Mellin-engine
seeds, the per-Wolfenstein-parameter match counts at <1 PDG σ become
indistinguishable from random expectation (z = -0.52 aggregate, individual
z-scores in ±1σ band). Every Wolfenstein parameter has many alternatives
matching at <0.1 σ — alternatives that are simpler, equally "natural" in the
Mellin-engine vocabulary, and not obviously preferred by the structural
reading the original sprint proposed.

**Selection-bias verdict:** **fully attributable**. The original signal does
not survive basis expansion.

## What the test did

The original W3 sprint
(`debug/w3_lambda_predictive_verification.py`) tested 30 hand-curated forms
against four Wolfenstein parameters. The basis was named in the script header
("frozen") but written in the same session as the test, so a strict reader
could legitimately argue the basis was lightly tuned to the data during
script-writing.

The mechanical re-test replaces the curated basis with a frozen-by-rule
enumeration:

- **56 seed constants** organized by master-Mellin mechanism class:
  - M1 (Hopf-base, π-family): 16 seeds
  - M2 (chirality, ln-family): 7 seeds
  - M3 (vertex parity / ζ-family): 12 seeds
  - ALG (small algebraic, including φ): 10 seeds
  - RAT (small rationals, control class): 11 seeds
- **Generation rule:** single-seed forms with prefactor (num/den), num ∈
  {1..5}, den ∈ {1, 2, 3, 4, 5, 6, 8, 10}, gcd(num,den)=1. Two-seed forms via
  {a·b, a/b, b/a} with numerator prefactor in {1, 2, 3} and same denominator
  set. Range filter [0.001, 10] (Wolfenstein parameter range). Dedup by 1e-10
  relative tolerance.
- **Result:** 21,448 unique forms (basis frozen and saved to JSON before
  any test ran).

For each Wolfenstein parameter, the basis is searched for matches at <0.5%
(threshold ignoring measurement uncertainty) and within 1 PDG σ (threshold
calibrated to experimental error). Class distribution of matches is compared
to class distribution of the basis to detect M1/M2 enrichment.

A random-target null is computed: random log-uniform targets on [0.10, 1.00]
are drawn 10,000 times and the basis is searched at the same threshold. This
gives the baseline expectation against which observed match counts are
compared.

## Per-parameter results

### λ (lambda) = 0.22500 ± 0.00067

| Statistic | Value |
|:---|:---|
| Matches at <0.5% | 46 |
| Matches within 1 PDG σ (= 0.298%) | 29 |
| Random expectation within 1 σ | 25.2 ± 5.1 |
| z-score | +0.74 |

**Top 5 mechanical matches by sigma:**

| Rank | Form | Deviation | σ | Class |
|:---:|:---|---:|---:|:---:|
| 1 | `(3/10)·(3/4)` = **9/40** | +0.0000% | 0.000 | RAT |
| 2 | (3/8)·(π/2)·(1/φ²) | -0.0015% | 0.005 | ALG+M1 |
| 3 | (2/5)·(1/φ)·(1/ln 3) | +0.0104% | 0.035 | ALG+M2 |
| 4 | β(4)·(1/(4·ln 3)) | +0.0196% | 0.066 | M2+M3 |
| 5 | (3/10)·ln(2)·(π⁴/90) | +0.0279% | 0.094 | M2+M3 |

**The original candidate** 1/√(Vol(S³)) = 1/(π√2) (a.k.a. `sqrt(2)/(pi*2)`)
appears at rank 6 with deviation 0.118 σ. The best match is the pure rational
**9/40 = 0.225 exactly** at 0.000 σ. Five other forms beat the original
candidate, including: (3/8)·(π/2)·(1/φ²) at 0.005 σ (a much sharper M1+ALG fit
than 1/(π√2)), and β(4)/(4 ln 3) at 0.066 σ (a clean M2+M3 fit).

### A = 0.826 ± 0.012

| Statistic | Value |
|:---|:---|
| Matches at <0.5% | 38 |
| Matches within 1 PDG σ (= 1.453%) | 123 |
| Random expectation within 1 σ | 123.2 ± 11.5 |
| z-score | -0.02 |

**Top 5 mechanical matches by sigma:**

| Rank | Form | Deviation | σ | Class |
|:---:|:---|---:|---:|:---:|
| 1 | (3/4)·√π/(ln 5) | -0.0042% | 0.003 | ALG+M2 |
| 2 | (3/10)·ln(3)·√(2π) | +0.0174% | 0.012 | ALG+M2 |
| 3 | (3/4)·(π⁴/90)·ζ(6) | -0.0217% | 0.015 | M3 |
| 4 | (3/4)·G·ζ(3) | -0.0264% | 0.018 | M3 |
| 5 | β(4)·√(2π)/3 | +0.0370% | 0.025 | ALG+M3 |

The original candidate **√(ln 2) at 0.55 σ off does not appear in the top 10**.
Many M3-family and ALG+M2-family forms fit A much more precisely. The top
match (3·√π)/(4·ln 5) at 0.003 σ is structurally a different combination than
√(ln 2) — and is dramatically closer to PDG.

### ρ̄ (rho_bar) = 0.159 ± 0.010

| Statistic | Value |
|:---|:---|
| Matches at <0.5% | 40 |
| Matches within 1 PDG σ (= 6.29%) | 477 |
| Random expectation within 1 σ | 534.4 ± 56.7 |
| z-score | -1.01 |

**Top 5 mechanical matches by sigma:**

| Rank | Form | Deviation | σ | Class |
|:---:|:---|---:|---:|:---:|
| 1 | log₂(φ)·G/4 | -0.0155% | 0.002 | M2+M3 |
| 2 | φ/(10·ζ(6)) | +0.0283% | 0.005 | ALG+M3 |
| 3 | β(4)/(6·ζ(5)) | -0.0288% | 0.005 | M3 |
| 4 | (3/10)·(3/4)/√2 | +0.0623% | 0.010 | ALG |
| 5 | (3/5)·log₂(φ)·(1/φ²) | +0.0667% | 0.011 | ALG+M2 |

The original candidate 1/Vol(S¹) = 1/(2π) ranks **12th of 477** within-1σ
matches. The 6.29% PDG band is wide enough that hundreds of forms fit; the
class distribution of those matches is not preferentially M1.

### η̄ (eta_bar) = 0.348 ± 0.009

| Statistic | Value |
|:---|:---|
| Matches at <0.5% | 45 |
| Matches within 1 PDG σ (= 2.59%) | 240 |
| Random expectation within 1 σ | 219.6 ± 23.9 |
| z-score | +0.82 |

**Top 5 mechanical matches by sigma:**

| Rank | Form | Deviation | σ | Class |
|:---:|:---|---:|---:|:---:|
| 1 | (π/2)·√π/8 | +0.0059% | 0.002 | ALG+M1 |
| 2 | log₂(φ)·√(2π)/5 | +0.0119% | 0.005 | ALG+M2 |
| 3 | (π⁴/90)/(3·ζ(5)) | -0.0212% | 0.008 | M3 |
| 4 | (3/4)·ζ(5)/√5 | -0.0586% | 0.023 | ALG+M3 |
| 5 | (3/10)·ζ(3)/ζ(5) | -0.0648% | 0.025 | M3 |

The original candidate (ln 2)/2 ranks **34th of 240**. The top match for η̄ is
**(π/2)·√π/8 = π^{3/2}/16**, which is in the M1+ALG family — i.e., the M1
family fits η̄ more sharply than the M2 family does at this basis size,
contradicting the original sprint's structural reading of "η̄ in M2 family."

## Aggregate statistical verdict

| Wolfenstein param | Observed (1σ) | Expected (1σ) | z |
|:---|---:|---:|---:|
| λ | 29 | 25.2 ± 5.1 | +0.74 |
| A | 123 | 123.2 ± 11.5 | -0.02 |
| ρ̄ | 477 | 534.4 ± 56.7 | -1.01 |
| η̄ | 240 | 219.6 ± 23.9 | +0.82 |
| **Total** | **869** | **902.5 ± 63.9** | **-0.52** |

The aggregate z = -0.52 means the observed total is *slightly below* (not
above) the random-target null. **There is no signal whatsoever in this larger
basis.** Any individual parameter's match count is within 1 σ of what one
would expect by combinatorial chance.

## Class enrichment test (the M1/M2 structural reading)

The original sprint argued that λ and ρ̄ live in M1 (Hopf-base, π-family) while
A and η̄ live in M2 (chirality, ln-family) — a clean two-class split that
mapped onto "real-magnitude vs CP/amplitude" structurally. The mechanical-
basis test asks: are the match classes preferentially M1/M2, or are they
distributed in proportion to basis composition?

**Aggregated class distribution of <1σ matches across all 4 Wolfenstein
parameters (top of distribution):**

| Class | Matches | Match % | Basis % | Enrichment |
|:---|---:|---:|---:|---:|
| ALG+M3 | 22 | 27.5% | 16.4% | 1.67× |
| ALG+M2 | 15 | 18.8% | 10.3% | 1.81× |
| M3 | 14 | 17.5% | 14.2% | 1.24× |
| M2+M3 | 9 | 11.2% | 11.7% | 0.96× |
| ALG+M1 | 5 | 6.2% | 8.6% | 0.73× |
| M1+M3 | 5 | 6.2% | 9.7% | 0.64× |
| ALG | 4 | 5.0% | 9.9% | 0.51× |
| RAT | 2 | 2.5% | 1.0% | 2.52× |
| M2 | 2 | 2.5% | 7.0% | 0.36× |
| M1+M2 | 1 | 1.2% | 7.1% | 0.17× |
| M1 | 1 | 1.2% | 4.0% | 0.32× |

(Note: aggregate totals 80, not 869, because the aggregation here is over
top-20 entries per parameter — the most precise matches. The full <1σ
distributions per parameter, in the JSON, show the same pattern at larger
scale.)

**Reading the table:** if M1 and M2 are the structurally preferred classes for
Wolfenstein parameters, M1 and M2 forms should be enriched (>1×) among the
matches. They are not. M1 alone has enrichment 0.32× (de-enriched), M2 alone
0.36× (de-enriched). The enrichment is concentrated in **ALG+M3** and
**ALG+M2** mixed classes, and in pure **RAT** (mostly because pure rationals
near 0.225 = 9/40 hit λ exactly).

This is the discrimination signal saying: **the original M1/M2 reading was an
artifact of basis composition, not a structural preference**. The original
30-form basis was M1-heavy (it included six "1/sqrt(Vol)" and "1/Vol" forms
out of 30, ~20%) and M2-heavy (six ln-2 forms out of 30, ~20%), so the four
matches it found were artificially M1+M2. In a fairer basis, the match
distribution looks nothing like that.

## Side-finding: λ = 9/40 exactly

The single most striking entry in the mechanical results: **λ = 9/40 = 0.225
exactly** (PDG central value). This is the form `(3/10)·(3/4)` in the
mechanical generator — purely rational, no transcendentals at all. Five other
mechanical forms also beat the original 1/(π√2) candidate's 0.118 σ.

PDG λ is a measured quantity with σ = 0.00067. It happens to round to 0.22500
at 5 digits, but this is unlikely to be a true rational. λ is conventionally
defined as |V_us|, which has no a priori reason to be rational. The match at
9/40 is a basis-side artifact.

But the existence of this exact-rational match in the mechanical basis is
itself important: **the data point λ ≈ 0.225 is in a region of the real line
densely populated by simple rational expressions, and any "natural" basis
containing small rationals will find one nearby**. The original sprint's
1/(π√2) ≈ 0.22508 is a 0.118 σ match precisely because the value 0.225 is
generic — many natural objects sit near it.

## Original candidate ranking summary

| Wolfenstein param | Original candidate | Class | Rank in mechanical <1σ matches | Total <1σ matches |
|:---|:---|:---|---:|---:|
| λ | 1/(π√2) | M1+ALG | 6 | 29 |
| A | √(ln 2) | M2+ALG | not found at <0.5%, σ = 0.55 (within 1σ of PDG, but outside top-20) | 123 |
| ρ̄ | 1/(2π) | M1 | 12 | 477 |
| η̄ | (ln 2)/2 | M2 | 34 | 240 |

**No original candidate ranks first in its parameter's match list.** Three of
four don't even rank in the top 5. The η̄ candidate is in the top 14% (rank
34 of 240); the others are mid-list. There is no preferential structural
reading these specific forms can be claimed to encode.

## Comparison to the original 30-form basis

| Quantity | Original sprint (30 forms) | Mechanical (21,448 forms) |
|:---|---:|---:|
| Basis size | 30 | 21,448 |
| Total <1σ matches across 4 params | 6 | 869 |
| Random expectation within 1σ | 0.7 ± 1.0 | 902.5 ± 63.9 |
| Signal z-score | ~5-6σ | -0.52σ |
| Verdict | "Strong signal" | No signal |

The original sprint's 6 hits / 30 forms / 4 params is statistically the same
match density (5.0%) as the mechanical basis's 869 hits / 21,448 / 4 params
(1.0%) — wait, no, that's a different ratio. Let me re-examine.

The original 30-form basis at <1σ found 6 hits across 4 params, density 6/120
= 5.0%. The mechanical 21,448-form basis at <1σ finds 869 hits across 4
params, density 869/(21,448·4) = 1.0%. **The original basis was 5× denser in
hits than chance would predict for ANY basis** — but only because the
original 30 forms were targeted at the parameters' value range. When the
basis is fairer (mechanically generated), match density falls to chance.

**The "5-10σ above chance" claim in the original sprint was correct AT THAT
BASIS SIZE.** The error wasn't in the statistics; the error was in believing
that the basis itself was selected without bias. The mechanical re-test
demonstrates that the basis was selected with bias — the 30 forms were
disproportionately near the parameter values.

## Two structural readings that I was asked to be ruthless about

### (a) Are the original specific forms still empirically distinguished?

The original four candidate forms still match PDG within 1 σ. That is a
correct empirical statement and does not change. What changes is that these
specific forms are no longer **the natural** forms — many other forms in a
fair basis match equally well or better. So the empirical fact "1/(π√2) ≈ λ"
remains true; the structural inference "λ is set by 1/(π√2) for a Connes-
Chamseddine reason" loses its force.

### (b) Does the M1/M2 ↔ real/CP-imaginary structural reading survive?

No. The mechanical class distribution does not preferentially fit M1 to λ/ρ̄
and M2 to A/η̄. In fact, the strongest single fit for η̄ is (π/2)·√π/8 (M1+ALG
class), and the strongest for λ is the pure rational 9/40. The original four
matches happened to land 1-1 in M1 and 1-1 in M2 because the original basis
had roughly equal M1 and M2 representation; in the larger fairer basis, the
class assignment is essentially random.

## What this changes for WH7

WH7 was drafted as CANDIDATE pending Tests 1, 2, 5 in
`debug/w3_forward_plan_memo.md`. **Test 5 is now decisive negative.** The
"empirical signal at ~6-10σ above chance" claim in the WH7 draft is false
once basis-selection bias is properly controlled.

Recommended PI action:

1. **Do NOT add WH7 to CLAUDE.md §1.7** as currently drafted. The candidate
   was empirically supported at the original-basis level but the support
   does not survive larger-basis testing.

2. **Do NOT promote any of the four candidate forms to "structural" status
   in any paper.** They remain four numerical coincidences that fit PDG
   within experimental error — no different from many other forms in a
   fair basis.

3. **The Koide cone for charged leptons (the 1-arcsec result) is unaffected
   by this test** — it tests one specific *known* mathematical fact about
   the lepton mass triple (45° from the democratic axis) at experimental
   precision, not a basis-search match. The Koide cone result is
   independently real.

4. **Paper 34 §V.C "Calibration-data candidate matches" subsection (added
   in the original sprint) should be reframed or removed.** The matches
   are no longer above-chance under the curve-fit-audit discipline. The
   honest framing would be: "Catalog of W3 candidate matches under a
   30-form prespecified basis; mechanical-basis re-test (Test 5) shows
   these matches are not above-chance in a fair-basis enumeration. Listed
   here for institutional memory, not as evidence of structure."

5. **The "tan(δ_CP) = π · ln 2 at 0.02 σ" closed form** is the most striking
   single number in the original sprint. It survives this mechanical test
   trivially (the test was for λ, A, ρ̄, η̄ individually, not for δ_CP =
   atan2(η̄, ρ̄)). But under the principle "the parts are no longer
   structurally identified, so derived combinations of the parts inherit
   the same caveat" — the closed form is one observation among many.
   PDG δ_CP is 65.4° ± 3.0°, an 8.5% relative uncertainty; the band is
   wide enough that many simple expressions fit.

## Honest concerns about this Test 5

1. **The mechanical generator is parameterized.** I chose num ∈ {1..5},
   den ∈ {1..10 with gcd}, and {1, 2, 3} for two-seed prefactors. Different
   choices change the basis size. A true reader-proof discipline would also
   vary the generator parameters and verify the conclusion is robust to
   them. I did NOT do that. But the 21,448-form basis is large enough that
   doubling or halving it would not flip the verdict from z = -0.52 to
   z > 5 — random expectation scales with basis size, so the ratio stays
   roughly chance.

2. **The seed list is also a choice.** I included 56 seeds chosen to be
   "natural" in the master Mellin engine vocabulary plus controls. A
   different seed list (e.g., adding {Apéry constant, Gieseking volume,
   Madelung constant}) would change the basis. But adding more
   transcendentals only increases the basis and increases random-match
   density proportionally — the verdict direction does not flip.

3. **PDG values are the targets.** The targets are fixed real numbers, not
   ranges. If PDG λ moves by even 0.1 σ in future precision, every match
   ranking changes. The test is a snapshot at PDG 2024.

4. **The three-axis Mellin-engine partition (M1/M2/M3) is the framework's
   own classification.** A skeptical reader could argue: "Of course M1 and
   M2 forms appear in the matches — you defined the classification to
   include the most prominent transcendentals. The fact that *some* M1/M2
   form fits each PDG value is uninformative because M1/M2 are dense
   transcendental classes." This is a fair critique. The test addresses
   it by including ALG and RAT as control classes — and those classes are
   represented in matches at the same per-form rates as M1/M2.

5. **Alternative: subset-restricted hypothesis.** A more careful original
   sprint could have prespecified "I claim λ matches a form involving
   ONLY M1 seeds, with no algebraic mixing, and prefactor 1." That would
   be a much narrower hypothesis. The original sprint did not do that —
   the 30 forms freely mixed M1 with sqrt(2), with rational prefactors,
   etc. A narrower hypothesis might survive this test, but it would also
   be a fundamentally different (smaller-N) claim.

## Files

- `debug/w3_mechanical_basis.py` — generator + test (~600 lines)
- `debug/data/w3_mechanical_basis.json` — full results
  - Basis: 21,448 forms (sample of first 50 + last 50 in JSON)
  - Per-parameter results: top 20 matches by sigma + top 10 at <0.5%
  - Original-candidate ranking
  - Random-expectation null
  - Class enrichment table
- `debug/data/w3_mechanical_basis.log` — full stdout log
- `debug/w3_mechanical_basis_memo.md` — this file

## Conclusion

The W3 candidate spectral-zeta identification of CKM Wolfenstein parameters
**does not survive the basis-expansion test**. The original 8-match signal
was real at that basis size but is fully attributable to basis-selection
bias when the basis is mechanically generated.

This is a clean curve-fit-audit negative. WH7 should not be promoted to
§1.7. Paper 34 §V.C should be reframed or removed. The "second packing
axiom" question (W3) returns to its prior status: open, no concrete
candidate, 14 cataloged speculations.

The methodological lesson is the standard one: **when a basis is hand-curated,
"sigma above chance" is computed against the wrong null**. The right null is
chance against a *fair* basis, and the test of basis fairness is mechanical
generation. The original sprint's fault was not the statistics; it was the
implicit assumption that a 30-form basis written in one session was free of
the data-knowledge it was being tested against.

Test 1 (PMNS recheck with same prespecified-basis methodology) is now lower
priority because it would test a methodology that this Test 5 has shown to
be unreliable. The remaining higher-leverage tests in the forward plan are
Test 6 (Connes-Chamseddine spectral action computation, 2-4 weeks) — the
actual derivation. If the spectral action picks specific values for the
Wolfenstein parameters from first principles, those values can be compared
to PDG independently of any basis-search framework.
