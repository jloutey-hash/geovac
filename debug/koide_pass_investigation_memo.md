# Koide cone focused investigation — memo

**Date:** 2026-05-18
**Author:** PM session, focused sub-agent dispatch
**Context:** Focused Koide-specific pass after Sprint W3 (2026-05-08) tested
general lepton mass spectrum and found Koide-cone-only signal. Framework
state now includes Papers 38–45 (Lorentzian propinquity, modular Hamiltonian,
inner-factor Mellin engine theorem 2026-05-07).
**Output files:**
- `debug/koide_pass_investigation.py` (analysis script, 100 dps mpmath)
- `debug/data/koide_pass_investigation.json` (structured results)

## Verdict (TL;DR)

**CLEAN NEGATIVE. Koide remains a surviving empirical anchor outside the
framework's autonomous scope.**

The Koide cone holds at sub-σ precision (K = 2/3 at −0.91σ, θ = 45° at
−0.95 arcsec) at PDG 2024 precision. No structural mapping to the framework
emerged from any of four angles tested. The "hits" found by the PSLQ search
are all tautological restatements of the Koide constraint itself, exactly
mirroring the Sprint W3 (2026-05-08) lepton mass recheck finding. The
framework's inner-factor Mellin engine theorem (CLAUDE.md §2, 2026-05-07)
predicts this outcome structurally: Yukawas live in the inner-factor
Dirichlet ring, categorically disjoint from the outer M1/M2/M3 mechanism,
and the framework does NOT autonomously select Yukawa values.

## Part 1: Numerical verification

| Quantity | Value | Target | Deviation |
|----------|-------|--------|-----------|
| K = ΣM / (Σ√M)² | 0.66666051148 | 2/3 = 0.66666666666 | −6.15 × 10⁻⁶ |
| K − 2/3 in σ_K (PDG, τ-dominated) | — | — | **−0.91σ** |
| cos²(θ) where θ = ∠(√M, 1) | 0.50000461643 | 1/2 | +4.62 × 10⁻⁶ |
| θ (degrees) | 44.99973550 | 45 | **−0.95 arcsec** |

Koide holds at sub-σ precision at PDG 2024 measurement uncertainty (τ
dominates: σ_τ = 0.12 MeV propagates to σ_K ≈ 6.77 × 10⁻⁶).

## Part 2: Focused PSLQ on Koide-specific quantities

### Basis design

7 targets (3 mass ratios + θ + cos²(θ) + 2 residuals) × 47 forms × 14
factors = 658 candidates per target. Basis is **hand-curated but structurally
motivated**:

- **M1 Hopf-base measure family** (14 forms): π/k for k ∈ {2,3,4,6,8,12,16},
  multiples 1/(2π²), 1/(4π), etc.
- **M2 Seeley-DeWitt family** (8 forms): √π, ζ(2), ζ(3), ζ(4), log 2, etc.
- **M3 Hurwitz / Dirichlet-L family** (4 forms): Catalan G, β(2), etc.
- **SU(3) Casimir rationals** (11 forms): 1/3, 2/3, 1/8, 1/9, 1/10, 27, etc.
- **Framework algebraics** (6 forms): κ² = 1/256, Δ = 1/40, 1/B = 1/42,
  F = π²/6, α
- **Golden ratio** (4 forms): φ, 1/φ, φ², φ−1 (selection-bias control)

Random-coincidence expectation at 1% tolerance: ~92 hits across 7 targets
(double-sided). Observed: **16 hits at 1%, 15 hits at 0.01%**.

### Headline finding: hits are tautological

All 15 hits at 0.01% tolerance fall into **two equivalence classes**, both
of which are literal restatements of the Koide constraint:

1. **θ_radians group (5 hits, all rel_err ≈ 5.9 × 10⁻⁶):**
   `4 · π/16`, `3 · π/12`, `2 · π/8`, `1.5 · π/6`, `1 · π/4`
   — all equal π/4. This is the Koide cone half-opening; the search basis
   contains π/16, π/12, π/8, π/6, π/4 so any rational multiple landing on
   π/4 trips the filter. **One fact, five expressions.**

2. **cos²(θ) group (5 of 7 hits, rel_err ≈ 9.2 × 10⁻⁶):**
   `1.5 · 1/3`, `0.75 · 2/3`, `1 · 1/2`, `3 · 1/6`, `4 · 1/8`
   — all equal 1/2. This is cos²(π/4) = 1/2; the SU(3) Casimir basis
   contains 1/3, 2/3, 1/2, 1/6, 1/8 so rational multiples landing on 1/2
   trip the filter. **One fact, five expressions.**

This is **exactly the same pattern** the Sprint W3 lepton mass recheck
(2026-05-08) found (memo `debug/w3_lepton_mass_recheck_memo.md`): the two
matches it logged were `2/3 · (φ² − φ = 1) = 2/3` and `4 · π/16 = π/4`,
both literally encoding K = 2/3 / θ = π/4 via different identities.

### Headline finding: zero hits on actual three-generation content

The three Koide-relevant mass ratios:
- r₁ = √(m_e/m_μ) = 0.06954…
- r₂ = √(m_e/m_τ) = 0.01696…
- r₃ = √(m_μ/m_τ) = 0.24385…

returned **0 hits at 0.1% and 0.01% tolerance**. The single 1% hit (r₂ vs
1/(6π²) at 0.42% rel_err) sits well inside the 92-hit random-coincidence
expectation and well outside PDG precision. **The actual generation content
of the lepton triple does not land on any framework-internal natural form.**

## Part 3: Structural angles

### (A) Master Mellin engine k ∈ {0,1,2} mapping

**No mapping.** The master Mellin engine produces a continuous family of
values 𝓜[Tr(D^k e^{−tD²})](s); selecting three discrete values
(m_e, m_μ, m_τ) would require three specific (k, t, s) selections with no
internal mechanism. The inner-factor Mellin engine theorem (CLAUDE.md §2,
2026-05-07) places Yukawas in the inner Dirichlet ring ℚ[y_i^{−2s}],
categorically disjoint from outer M1/M2/M3. Koide is an inner-factor
constraint and is not addressable by the master Mellin engine.

### (B) SU(3) 45° check

**No natural 45°.** SU(3) simple roots have **120°** between them (Cartan
matrix of A₂, verified: angle(α₁, α₂) = 120° to 10 dp). Fundamental weights
have 60° to each other. The natural angles in the SU(3) weight lattice are
60° and 120°, not 45°. The Koide cone is **not** an SU(3) representation-
theoretic fact.

### (C) Inner-factor Mellin engine theorem

**Tautological reformulation, not derivation.** Koide CAN be expressed as
a ratio of inner Dirichlet values at half-integer s:

- Σᵢ y_i^{−2s} at s = −1/2 gives Σᵢ y_i = m_e + m_μ + m_τ
- Σᵢ y_i^{−2s} at s = −1/4 gives Σᵢ √y_i = √m_e + √m_μ + √m_τ
- K = (sum at s=−1/2) / (sum at s=−1/4)²

This reformulation **expresses** Koide without deriving it: the K = 2/3
constraint becomes a polynomial relation among the three Yukawa values
that does not follow from any framework axiom. Half-integer s sits outside
the integer-s regime where the inner Dirichlet ring is rigorously
characterized; analytic continuation does not help because the values are
just the unconstrained sums Σ y_i and Σ √y_i. The inner-factor theorem
says the inner ring is **constrained** but **not determined** by the
framework — Koide's K = 2/3 is exactly the 1-parameter constraint the
framework cannot autonomously generate.

## Part 4: Overall verdict and falsification protocol

**Verdict: CLEAN NEGATIVE, as expected.**

Three convergent results:

1. **Numerical**: Koide holds at sub-σ precision, consistent with PDG 2024
   and unchanged from Sprint W3 / W3-LepTriple measurements.

2. **PSLQ**: 15/15 hits at 0.01% tolerance are tautological restatements of
   K = 2/3 / θ = π/4 (the Koide constraint expressed redundantly through
   the basis × factor expansion). Zero hits on actual three-generation
   content (r₁, r₂, r₃ sqrt-mass ratios). Aggregate signal indistinguishable
   from random expectation (16 observed vs 92 expected at 1%).

3. **Structural**: No k ∈ {0,1,2} mapping in master Mellin engine; no
   natural 45° in SU(3); inner-factor reformulation expresses but does
   not derive Koide.

The framework's structural-skeleton scope (CLAUDE.md §2 pattern
crystallization, 2026-05-07) is confirmed: GeoVac determines selection
rules, transcendental rings, and structural identities cleanly, but does
not autonomously generate inner-factor calibration data. Koide joins the
list of irreducible-but-natural empirical anchors that sit outside the
framework's autonomous scope, alongside K = π(B + F − Δ) (the α near-miss),
the L2 next-order constant c ≈ 4.10932…, S_full(GS) atomic correlation
entropy, and the CKM/PMNS mixing parameters.

**Falsification protocol** (when the verdict could change):
- A new structural axiom on the inner Dirichlet ring ℚ[y_i^{−2s}] that
  forces a non-trivial polynomial relation among the y_i would change
  the verdict. The 45°-cone constraint would have to follow as a corollary.
- A geometric construction that produces the (ℂ ⊕ ℍ)³ lepton sector with
  an internal 45°-rotation symmetry would be a candidate (the
  W3-LepTriple memo flagged this as the natural second-packing-axiom
  target). No concrete proposal has emerged in 7 sprints of attempts.
- A non-perturbative spectral-action computation on the Connes–Chamseddine
  almost-commutative geometry that derives a 45°-constraint at the
  generation level would also flip the verdict. Such a computation has
  not been attempted in the literature to our knowledge.

**Status**: WH7 is not promoted. Paper 2 (Observations) is unaffected.
Paper 34 §V.D "Literature convention exposures" is unaffected — Koide is
not a convention exposure but an external empirical fact the framework
acknowledges and does not address.

**Recommendation**: do not open another Koide-specific sprint until either
a concrete structural proposal emerges (e.g., a Lie-algebraic 45° construction
in (ℂ ⊕ ℍ)³ or a non-perturbative inner-factor axiom). The W3 closure
remains the correct standing position.

## Files

- `debug/koide_pass_investigation.py` (~280 lines, 100 dps mpmath)
- `debug/data/koide_pass_investigation.json` (structured results, all hits)
- `debug/koide_pass_investigation_memo.md` (this memo)

No production code modified. Paper edits not warranted (verdict is clean
negative, no new structural content).
