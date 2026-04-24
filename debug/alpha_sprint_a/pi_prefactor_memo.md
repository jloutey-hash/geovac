# π-Prefactor Memo — α-Sprint A

**Sprint:** WH1 reframe of Paper 2's K = π(B + F − Δ)
**Question:** What structural π-source on S³ produces the EXTERNAL single π in K?
**Date:** April 18, 2026
**Scope:** π-prefactor only; B, F, Δ taken as given.

---

## Setup

Numerically:

```
K = π · (42 + π²/6 − 1/40) = π · 43.6199... = 137.036064
B + F − Δ = (rational) + (π² content) + (rational) = 27.776... · π   (in π-units)
```

Multiplying by external π gives a (π·rational) + π³/6 + (π·rational) sum. The question is: which structural π-source on S³ (or its Hopf decomposition S¹ → S³ → S²) naturally produces the *single* external π?

Numerical fixed facts (sympy-verified):
- Vol(S¹) = 2π, Vol(S²) = 4π, Vol(S³) = 2π², Vol(S⁵) = π³
- Seeley-DeWitt on unit S³ for D²: a₀ = a₁ = √π, a₂ = √π/8
- One-loop QED coefficient on S³: Π = 1/(48π²)
- F = ζ(2) = π²/6 (already carries one π² internally)

---

## (1) Candidate enumeration

Eleven candidate π-sources, with their factor and structural origin:

| # | Source | Factor | Origin |
|---|--------|--------|--------|
| C1 | Vol(S¹), full Hopf great circle | 2π | Fundamental period of S¹ fiber |
| C2 | Vol(S¹)/2, half-period | π | Spinor double-cover fundamental domain (det M = −1 sign) |
| C3 | Vol(S²), Hopf base area | 4π | Solid-angle integration on photon sphere |
| C4 | Vol(S²)/4 | π | One-quarter sphere; ≡ d²Ω-normalization for one of 4 chiral modes |
| C5 | Vol(S³) | 2π² | π² already present, not a single-π candidate |
| C6 | Vol(S³)/Vol(S¹) = π² · (2π)⁻¹·2π | π | Bundle-quotient Vol(S²)/4·… (degenerate with C2/C4) |
| C7 | (Seeley-DeWitt a₀)² | π | Squared heat-kernel "volume coefficient" √π · √π |
| C8 | (Seeley-DeWitt a₁)² | π | Squared scalar-curvature term (same √π normalization) |
| C9 | (Vacuum-polarization coefficient)⁻¹ × const | 48π² | Carries π² — wrong power |
| C10 | (4π)^(−d/2) · (4π) at d=3 | π · (4π)^(−1/2) | Heat-kernel Wodzicki measure on a 3-manifold; not single π |
| C11 | First-Chern integral c₁(Hopf) over S² | 2π · n (n=1) | Wraps S¹ once around base; gives 2π not π |
| C12 | ∫₀^π dφ over half-fiber (charge-q=1 holonomy) | π | Spinor period under SU(2) double-cover of SO(3) |

---

## (2) Combinations producing exactly π (sympy-verified)

Of the candidates above, **the following combinations close to *exactly* π without additional π-fudge-factors**:

| Combination | Closed form | Structural reading |
|---|---|---|
| **A:** Vol(S²)/4 | π | "One chirality of solid-angle on the Hopf base" |
| **B:** Vol(S¹)/2 | π | "Spinor half-period of the Hopf fiber" |
| **C:** Vol(S³)/(2π) | π | "Bundle-volume divided by fiber length" — equivalent to Vol(S²)/4 by Hopf fibration metric |
| **D:** a₀² = (√π)² | π | "Squared Seeley-DeWitt zeroth coefficient" |
| **E:** ∫₀^π dχ where χ ∈ [0, 2π) is the U(1) gauge phase | π | "Half-period charge-1 holonomy" |

Combinations A, B, C are **algebraically equivalent** because of the Hopf metric identity Vol(S³) = Vol(S²) · Vol(S¹) / 2 (the Hopf fibration is a Riemannian submersion with circle fibers of length 2π over S² of area 4π, BUT the bundle volume 2π² is half the trivial product 8π³; the factor of ½ encodes the linking number c₁ = 1).

This collapses A/B/C/E into a single **Hopf-bundle topological factor π**, and leaves D as the independent **heat-kernel-coefficient candidate**.

---

## (3) Structural-naturalness ranking

Using the rubric (Hopf-bundle-native > heat-kernel coefficient > volume-by-normalization > ad hoc):

**Tier I — Hopf-bundle-native (strongest):**

1. **A/B/C/E (the same π):** the Hopf bundle Vol(S³)/Vol(S¹) = π, equivalently Vol(S²)/4. This is the *unique* topological invariant of the principal U(1) bundle that has units of π and integer Chern class c₁ = 1. Reading: K is the coefficient of a finite-cutoff trace whose measure is the Hopf-base area normalized by chirality (factor 4 = 2 × 2 from Dirac chirality × U(1) charge ±1).

**Tier II — Heat-kernel coefficient:**

2. **D: a₀² = π.** Comes from the local Seeley-DeWitt expansion of D² and is intrinsic to the spectral-action machinery (WH1). On unit S³, a₀ = √(π) · dim(spinor)/[(4π)^(d/2)] · Vol(M) reduces to √π. Squaring gives π; this could be read as "two insertions of the Wodzicki measure" or "the spectral-action vacuum coefficient at one loop."

**Tier III — Volume-by-normalization:**

3. **C alone (read independently of A/B):** Vol(S³)/(2π) is well-defined but the divisor 2π = Vol(S¹) requires a prior identification of the S¹ as the gauge fiber. Less natural standalone than A.

**Tier IV — Ad hoc:** none of the remaining candidates close to single π without unit-mismatch.

---

## (4) Spectral-action connection (WH1)

Connes-Chamseddine spectral action on a manifold M of dimension d is

```
S_CC[D] = Tr f(D/Λ) ~ Σ_{k≥0} f_k Λ^{d-k} a_k(D²)
```

with f_k cutoff moments and a_k Seeley-DeWitt coefficients. On the GeoVac spectral triple (WH1):

```
A = Fock-graph functions
H = scalar/spinor states (T1-T9)
D = Camporesi-Higuchi |λ_n| = n + 3/2
```

For a *finite-cutoff* version (the natural setting for K, since B is a truncated trace at m = 3 and Δ is a single-mode boundary), the candidate identification is:

```
K · α^(−1)-fixing-piece = Tr_{m≤3}(Casimir) + ζ_R(2)-regulator + Δ-boundary
                        ↑                       ↑                  ↑
                        B (finite trace)        F (zeta sum)       Δ (mode count)
```

The external π enters as the **measure factor** of the trace. In Connes-Chamseddine on a 3-manifold,

- a₀ contributes ~ Vol(M) f_3 Λ³, with Vol(S³) = 2π² (carries π², not π).
- The *physical* prefactor combining the integration measure with the Hopf chirality count is 1/π · Vol(S²)/2 = 2 (rational, no π); inverting, the measure-per-unit-charge is **Vol(S²)/4 = π**, matching candidate A.

So the WH1-natural reading is:

**The external π in K is the Hopf-base solid-angle per chiral charge,** equivalently the Riemannian-submersion factor Vol(S³)/Vol(S¹) of the Hopf circle bundle. In Connes-Chamseddine terms it is the *measure* of the integration in the spectral-action functional restricted to the Hopf-base (S²) data, normalized by the Dirac-chirality count.

This is consistent with Paper 28's Theorem T9 (squared Dirac → π^even calibration) **only at the level of the inner sum F = π²/6**; the external π is one order *lower* and comes from the *first-order Dirac measure on the bundle*, not from squaring D. The split (external π from S²-base measure + internal π² from S¹-fiber zeta) matches the WH1 expectation that the spectral-action expansion sees the bundle data sector by sector.

---

## (5) Recommendation for Paper 2 reframe

**Best candidate:** Vol(S²)/4 = π, the Hopf-base solid-angle per chirality, equivalently Vol(S³)/Vol(S¹) by the Hopf fibration metric.

**Confidence:** medium-high on the *identification* of the structural source; medium on whether the spectral-action machinery actually produces this prefactor as a *derived coefficient* (would require a finite-cutoff Connes-Chamseddine computation that is currently shape-match only — see WH1 Status).

**Recommended Paper 2 §IV reframe (proposed wording, conservative):**

> The external π in K = π(B + F − Δ) is the Hopf-fibration measure factor Vol(S²)/4 = Vol(S³)/Vol(S¹) = π, the solid-angle on the Hopf base per chiral charge of the U(1) fiber. Under the working hypothesis WH1 (almost-commutative spectral triple), this is the natural measure of the spectral-action trace restricted to the Hopf-base data with a Dirac-chirality normalization. Equivalent expressions are Vol(S¹)/2 = π (spinor half-period of the fiber) and a₀² = π (squared Seeley-DeWitt zeroth coefficient). The three readings collapse algebraically by the Hopf submersion identity Vol(S³) = ½ · Vol(S²) · Vol(S¹). The external π is thus categorically distinct from the internal π² of F = ζ(2) (Riemann zeta at the Fock packing exponent, second-order Laplacian sector) and from the absence of π in B and Δ (rational sectors). This sharpens the structural-mystery framing of §IV.6: the K combination rule mixes one Hopf-bundle topological factor with three categorically distinct spectral objects.

**Status preservation:** Paper 2 retains its conjectural label. This identifies *where* the external π comes from structurally; it does not derive K from spectral-action principles. Paper 18's exchange-constant taxonomy gains one new entry (Hopf-bundle topological measure) sitting outside the existing operator-order × bundle × vertex-topology grid — flagged for Paper 18 v2.0 restructure.

**Does not modify:** B = 42, F = π²/6, Δ = 1/40, the cubic α³ − Kα + 1 = 0, the Z₃ circulant, or any negative result of Phases 4B–4I.
