# Sprint K-CC: Is K = π(B + F − Δ) one CC spectral-action coefficient?

**Date:** 2026-05-03
**Type:** Focused scoping sprint
**Verdict:** **CLEAN NEGATIVE on all three sub-tracks. WH5 confirmed.**

## What we tested

Paper 35's algebraic/observation split makes Paper 2's α conjecture
sharper. The explicit π in K = π(B + F − Δ) is presumably a
"temporal/observation projection π" (Paper 34 §III.15: 2π·ℚ per Matsubara
mode, π^(2k)·ℚ in integrated quantities). The Connes–Chamseddine spectral
action Tr f(D²/Λ²) is structurally a single such temporal projection — Λ
has dimension energy and the trace runs over the spectrum. So if K is one
CC-spectral-action coefficient (rather than three independent projections,
the WH5 reading), it should be ONE projection, not three, and Paper 34's
projection-depth model is not violated by the 8.8×10⁻⁸ residual.

This sprint asked: can K be derived as a single CC-projection statement on
unit S³? Three parallel sub-tracks tested this independently. Any clean
positive resolves the K anomaly; all three clean negatives confirm WH5
(K is irreducible projection coincidence below the framework's intrinsic
resolution).

The setting was fixed by Sprint A (April 2026), which established the
SD two-term exactness theorem on D² on unit S³:

   K_heat(t) = Tr exp(-tD²) = (√π/2) t^{-3/2} − (√π/4) t^{-1/2} + O(e^{-π²/t}),

so setting Tr exp(-D²/Λ²) = K/π with K = 1/α reduces to the depressed cubic

       2 Λ³ − Λ = 4 (K/π) / √π,                                          (*)

with Cardano root Λ_∞ ≈ 3.7102448927… (CODATA 2018 α^{-1} = 137.035999084).
Sprint A's verdict was that the match is tautological (IVT + monotonicity);
this sprint hardens that verdict by checking whether Λ_∞ has any independent
home, and whether B, F, Δ unify in one heat-kernel computation.

## (a) Λ_∞ number-theoretic analysis (PSLQ)

**Verdict: NEGATIVE. Λ_∞ has no integer-relation identification in
framework or spectral invariants.**

PSLQ at 100 dps (120 dps headroom) against five basis sets:

| Basis | Size | Identification (Λ-coeff ≠ 0)? |
|:------|:----:|:------------------------------|
| Framework invariants {B, F, Δ, κ, |λ_n|, g_n, π, √π}  | 25 | **None**. PSLQ finds tautologies (g_2 + g_4 = g_5 = 42 ✓; B + F − Δ − (B+F) + Δ = 0; π = (√π)²) but no relation with nonzero Λ coefficient. |
| Spectral invariants {Vol(S³), SD coefficients, Hurwitz ζ at half-integer args} | 27 | **None**. Same outcome — intra-basis identities only. |
| Bernoulli polynomials at half-integers | 18 | PSLQ failed (B_1(1/2) = 0 in basis). |
| Pure π powers {π^{-2}, …, π², √π powers} | 10 | **None** at 10⁻⁹⁰ tol; trivial π = (√π)². |
| Combined products B^p F^q Δ^r π^s √π^t, p,q,r,s,t ∈ {-1,0,1} | 244 | **None** (one trivial π/√π = √π identity). |

A refined script (`kcc_lambda_pslq_v2.py`) forced the Λ coefficient to be
nonzero. Results:

- **PSLQ on {Λ, Λ³, 1, K/(π√π)} ⇒ +1 · Λ − 2 · Λ³ + 4 · K/(π√π) = 0** at
  10⁻⁹⁰ tolerance (positive control: this is exactly the depressed cubic
  (*), as expected from Sprint A).
- **PSLQ on {Λ, Λ³, 1, B, F, Δ, B+F−Δ, π, π², √π, 1/π, 1/√π}** without K:
  no identification at any tolerance. Λ_∞ is **not** algebraically expressible
  in {B, F, Δ, π} alone — its dependence on K = 1/α is irreducible.
- **PSLQ on {Λ, π powers}** alone: no identification.

So Λ_∞ depends on K non-trivially: removing K from the basis kills the
identification. The cubic (*) is the only closed form Λ_∞ has.

## (b) Natural Λ from S³

**Verdict: NEGATIVE. No natural S³ selection criterion (independent of K)
lands within 10⁻³ of Λ_∞.**

Three categories of criterion:

**(A) Smooth-cutoff Casimir trace targets.** N_smooth(Λ) = Σ g_n^Dirac χ(|λ_n|/Λ)
with cutoff χ ∈ {sharp, Gaussian, exp, (1+x²)^{-1}, e^{-x²}(1+x²)} solved
for Λ at natural targets {B = 42, g_3 = 40, Σ g_k = 80, Vol(S³) = 2π², B+F}:

- target = K/π = 43.62 (Gaussian): hits Λ_∞ to ~1e-7 — **tautological** by
  Sprint A's SD identity. **The Gaussian-cutoff trace at Λ_∞ literally
  equals K/π by construction**: at Λ_∞ the SD asymptotic (√π/2)Λ³ − (√π/4)Λ
  agrees with the full sum to 10⁻⁵⁵ (verified in `kcc_gaussian_check.py`),
  so this "hit" is the depressed cubic restated. Not an independent criterion.
- target = B+F = 43.6449 (Gaussian, **no K input**): hits Λ_∞ at relative
  distance 1.86×10⁻⁴. Cleanest non-tautological near-miss.
- target = B = 42 (Gaussian): Λ = 3.6649, **1.22% off** Λ_∞.
- target = g_3 = 40 (Gaussian): Λ = 3.6073, **2.78% off** Λ_∞.
- target = Σ g_k = 80 (Gaussian): Λ = 4.523, **22% off**.
- target = Vol(S³) (Gaussian): Λ = 2.873, **23% off**.

The B+F near-miss at 1.86×10⁻⁴ is suggestive but does NOT meet the 10⁻³
threshold-with-margin standard the sprint plan defined. It is also still
tied to the same cubic structure (it's just choosing a slightly different
Λ along the same SD curve). And the residual 1.86×10⁻⁴ is **3 orders of
magnitude looser** than the 8.8×10⁻⁸ residual that K = π(B+F−Δ) achieves
against 1/α — so even if "Gaussian trace = B+F" were elevated to a derivation,
it would not explain the K precision.

**(B) Composite-eigenvalue candidates.** Λ = |λ_n|, √(|λ_n|² ± 1), midpoints
(|λ_n| + |λ_m|)/2, |λ_n| ± 1, |λ_n| ± 1/2: 37 candidates, **0 within 10⁻³**.
Best candidate is √(|λ_2|² + 1) = √(7² /4 + 1) = √(53/4) ≈ 3.640, **1.89% off**.

**(C) Critical-point search.** The smooth-cutoff heat trace
Σ g_n exp(-(|λ_n|/Λ)²) is monotone in 1/Λ², so it has no interior critical
point. The sharp-cutoff Casimir density is a step function with no smooth
critical point. **No natural functional has a stationary point at Λ_∞.**

## (c) Unified CC heat-kernel for B, F, Δ

**Verdict: NEGATIVE. B, F, Δ are NOT three terms of one CC heat-kernel
expansion.**

Closed forms verified symbolically via the half-integer Hurwitz formula
ζ(s, 3/2) = (2^s − 1) ζ(s) − 2^s:

   ζ_{D²}(s) = 2 ζ(2s−2, 3/2) − (1/2) ζ(2s, 3/2)

   ζ_{D²}(2) = π² − π⁴/12       ≈ 1.7521
   ζ_{D²}(3) = (computed numerically; mixed Hurwitz)
   ζ_{D²}(4) = π⁶ (168 − 17π²)/1260   ≈ 0.16536
   ζ_{D²}(5) ≈ 0.04246
   ζ_{D²}(6) ≈ 0.01193

**B = 42** is the Paper 7 §VI scalar Casimir trace truncated at m = 3:
B(m) = Σ_{n≤m, l<n} (2l+1) l(l+1). Sanity verified, B(3) = 42 exact. This is
a **finite truncation of a different operator** (scalar Laplace–Beltrami,
not Dirac D²), separate from any heat-kernel coefficient.

**F = π²/6** does NOT appear as ζ_{D²}(s) at any integer s. The closest is
ζ_{D²}(2) = π² − π⁴/12, which equals (6F) − F²·(6/π² · π²/12)... but is
**not** 1·F or any small-integer multiple of F. Also F does not appear in
the SD coefficients a_0 = a_1 = √π, a_2 = √π/8 (T9 theorem: SD coefficients
of D² on unit S³ live in √π · ℚ, no π² content).

**Δ⁻¹ = 40** is g_3^Dirac = 2(n+1)(n+2) at n=3 — the **single-level
degeneracy** at the cutoff edge. It is not a heat-kernel coefficient at all
(it appears as the Euler–Maclaurin upper-boundary term of the truncated
mode-count sum, per CLAUDE.md Phase 4H finding F2, but that is a discrete
boundary value, not an SD or ζ coefficient of the heat trace).

The three live in three different spectral objects:

| Object | Where it lives | Computation type |
|:-------|:---------------|:-----------------|
| B = 42 | Scalar Laplace–Beltrami on S³, finite truncation at m=3 | Finite Casimir sum |
| F = π²/6 | Scalar Fock-degeneracy Dirichlet at exponent d_max = 4 | Infinite Dirichlet ζ_R(2) |
| Δ⁻¹ = 40 | Dirac D² spectrum, single-level degeneracy at n=3 | Boundary mode count |

These are genuinely **three computations**, on **two different bundles**
(scalar/spinor), of **three different types** (finite trace / infinite series /
single-level value). No single CC heat-kernel expansion unifies them.

## Verdict on the K hypothesis

**Negative across all three sub-tracks. WH5 confirmed.**

K = π(B + F − Δ) is **not** a single CC spectral-action coefficient on unit
S³ at any natural Λ. Specifically:

- The Λ that solves Tr exp(-D²/Λ²) = K/π exists (Sprint A) but is not
  selected by any natural S³ criterion independent of K (sub-track b).
- Λ_∞ has no number-theoretic identification in {B, F, Δ, π, √π, |λ_n|, g_n}
  (sub-track a). Its only closed form is the depressed cubic (*) which
  takes K as input.
- B, F, Δ do not appear as three terms of one heat-kernel expansion
  (sub-track c). They are three independent spectral computations on two
  bundles with three categorically different object types.

This **strengthens** WH5 (the Working Hypothesis that K is a coincidence
of three independent projections at one finite cutoff, not derivable from
a single principle): the only "single principle" interpretation, the CC
spectral action with f = exp(-x²), is a tautological rewrite of K = K with
no autonomous selection of Λ.

The Gaussian-cutoff "B + F" near-miss at 1.86×10⁻⁴ is the cleanest
non-tautological hint, but its residual is 4 orders of magnitude looser
than the K residual itself, so it cannot explain the K precision either.

## Implications for WH1 and WH5

**WH5 (K is projection coincidence, not derivable from one principle):**
**Strengthened.** Three more mechanisms eliminated. The CC framework
provides the right functional form (SD two-term exact) but does not
autonomously select Λ. K stays the irreducible "why does the sum equal α⁻¹"
question, with each ingredient (B, F, Δ) having a clean spectral home but
the addition rule resisting unification.

**WH1 (GeoVac is an almost-commutative spectral triple):** **Unchanged.**
Sprint K-CC tests CC spectral-action *derivation* of K, not the spectral-triple
*structure*. WH1 remains MODERATE-STRONG (per Sprint A α-LS Marcolli–vS
literature grounding). The fact that K is not one CC coefficient doesn't
demote the spectral-triple framing; it just confirms that the K combination
rule is below the framework's intrinsic resolution. Marcolli–vS gauge
networks (the closest published precedent for our structure) also do not
predict α; that prediction remains GeoVac's open conjecture.

## Implications for Paper 2 / Paper 35

**Paper 2 (Observations status, per 2026-05-02 audit memo):** **No change
recommended.** Paper 2 already documents K as a numerical coincidence
without first-principles derivation. Sprint K-CC adds three more eliminated
mechanisms to the "ways K is not derived" record:

   (i) Λ_∞ has no PSLQ identification in framework/spectral invariants.
   (ii) No natural S³ selection criterion independent of K picks Λ_∞.
   (iii) B, F, Δ do not unify in one heat-kernel expansion.

Suggested addition to Paper 2's open-questions section: a short remark
that the CC framework provides the exact two-term form on unit S³
(Sprint A) but does not autonomously select Λ; the K = π(B+F−Δ) match
remains structurally a three-projection coincidence in the WH5 sense,
not a one-projection CC coefficient.

**Paper 35 §VII.2 (if it speculates K = single CC projection):**
**Should be sharpened to clean negative**, citing Sprint K-CC. Specifically,
the temporal-projection π in K is structurally consistent with a single
Matsubara/temporal projection (per Paper 34 §III.15), but the CC heat-kernel
realization fails: there is no Λ selection principle, and B, F, Δ live in
three different sectors. The Paper 34 projection-depth model
(error compounds with depth) is **not** rescued by collapsing K to one
projection — that collapse fails. K stays the §VIII.3 anomaly (3-projection
chain producing 6-orders-below-tolerance residual).

## Honest limits

1. **PSLQ basis completeness.** A negative PSLQ result rules out
   identifications in the chosen basis. We tested 25 framework invariants,
   27 spectral invariants, 244 small-exponent products, plus pure π powers.
   We did not test, e.g., Catalan's constant G, Apéry's ζ(3), products of
   half-integer Hurwitz ζ values, or Bernoulli polynomials evaluated at
   irrational arguments. A nontrivial Λ identification in some richer basis
   cannot be ruled out by this sprint. However, the structural negative on
   sub-track (c) (B, F, Δ in three different sectors) is independent of
   PSLQ and stands.

2. **CODATA precision in Λ_∞.** CODATA 2018 α^{-1} = 137.035999084(21) gives
   Λ_∞ to ~12 digits. Sprint A's quoted reference Λ = 3.7102454679060528505
   (20 digits) corresponds to a different (more precise or differently sourced)
   K value. The 7th-decimal disagreement between the two reference values is
   ~5.8×10⁻⁷, comparable to the CODATA uncertainty propagated through the
   cubic. This affects only the very last digits of Λ_∞ and does not change
   any PSLQ verdict at the 10⁻⁹⁰ tolerance.

3. **Cutoff function space.** Sub-track (b) tested 5 cutoff shapes
   {sharp, Gaussian, exp, (1+x²)^{-1}, e^{-x²}(1+x²)}. The CC formalism
   in principle allows any positive Schwartz function; we did not search
   the full space. However, the SD two-term identity is universal across
   smooth cutoffs (Sprint A), so the depressed-cubic structure is robust;
   the question is whether some special cutoff selects Λ_∞ from its target
   value — and the search found no such cutoff selecting Λ_∞ from any
   K-independent target.

4. **Other gauge-network frameworks.** Marcolli–vS gauge networks are the
   closest published precedent (Sprint A α-LS finding), but they do not
   predict α. Other frameworks (Connes–Chamseddine NCG Standard Model,
   Chamseddine–Mukhanov "geometry from quanta") use the spectral action to
   FIT the gauge-coupling unification scale; they too do not autonomously
   select Λ. So the negative here is consistent with the broader NCG
   literature, not a peculiarity of GeoVac's setup.

## Files

- `debug/kcc_lambda_pslq.py` + `debug/data/kcc_lambda_pslq.json` — sub-track (a)
- `debug/kcc_lambda_pslq_v2.py` + `debug/data/kcc_lambda_pslq_v2.json` — refined PSLQ with Λ-coefficient constraint
- `debug/kcc_natural_lambda.py` + `debug/data/kcc_natural_lambda.json` — sub-track (b)
- `debug/kcc_gaussian_check.py` + `debug/data/kcc_gaussian_check.json` — Gaussian-tautology verification
- `debug/kcc_unified_heatkernel.py` + `debug/data/kcc_unified_heatkernel.json` — sub-track (c)
- `debug/kcc_sprint_memo.md` — this file
