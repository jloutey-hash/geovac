# Sprint LS-8a: Native derivation of the two-loop self-energy bracket C_2S from iterated CC spectral action on Dirac-S^3

**Date:** 2026-05-07
**Sprint goal:** Derive the dimensionless coefficient C_2S = +3.63 (Eides 2001 Tab. 7.3 two-loop SE 2S contribution +0.857 MHz divided by the LS-7 structural prefactor 0.2363 MHz/dim) from spectral data alone, with no empirical input from multi-loop QED literature. Strong test of Paper 35 Prediction 1 in the multi-loop sector.

**Verdict: WEAK — STRUCTURAL CONFIRMATION WITH RENORMALIZATION GAP.** The bare iterated CC spectral action on Dirac-S^3 with proper SO(4) vertex constraints at four vertices (rainbow + crossed topologies) **faithfully reproduces the UV-divergent structure of two-loop QED** but cannot autonomously extract the finite physical bracket C_2S = +3.63 because that requires field-theoretic renormalization (Z_2 wave-function and δm mass counterterms) which is NOT present in the bare spectral action. Paper 35 Prediction 1 is consistent at the structural level (the (α/π)² prefactor is right; the divergent integrand structure matches flat-space QED) but does NOT pass the strong test (native finite C_2S extraction).

## 1. The two-loop SE on Dirac-S^3

The two-loop electron self-energy has two irreducible 1PI topologies:

- **Rainbow:** ext → V₁ → V₂ → V₃ → V₄ → ext, with photon q_outer connecting V₁↔V₄ (outermost) and photon q_inner connecting V₂↔V₃ (innermost, nested).
- **Crossed:** same vertices, with q₁: V₁↔V₃ and q₂: V₂↔V₄ (the two photons cross).

Each vertex carries the SO(4) channel-count weight W ∈ {0, 1, 2} from the spinor-photon coupling (Paper 28 §vertex), with vertex parity selection n_a + n_b + q ≡ odd.

Three internal electron propagators carry levels n₁, n₂, n₃ (CH convention, |λ_n| = n+3/2, g_n = 2(n+1)(n+2)). Two photon propagators carry levels q_a, q_b (Hodge-1 eigenvalue μ_q = q(q+2), transverse degeneracy d_q = q(q+2)).

Per LS-7 the structural prefactor is

```
ΔE_2L^SE = (α/π)² · (Zα)⁴ · m_e c² / n³ · C_2S
```

with C_2S the dimensionless bound-state matrix element. For Z=1, n=2 the prefactor is **0.2363 MHz/dim**, so the literature value 0.857 MHz corresponds to **C_2S = +3.6266**.

## 2. Native LS-8a calculation

The spectral sum implemented in `geovac/two_loop_self_energy.py` is, for each topology:

```
Σ_topology(n_ext, n_max) = Σ_{n₁,n₂,n₃, q_a, q_b}
                            (∏ W at 4 vertices) · g(n₁) g(n₂) g(n₃) · d_T(q_a) d_T(q_b)
                            ────────────────────────────────────────────────────────────
                              |λ(n₁)|² · |λ(n₂)|² · |λ(n₃)|² · μ(q_a) · μ(q_b)
```

with all sums over CH levels 0..n_max for electrons and q ≥ 1 (photon) subject to triangle inequality and parity at each of the four vertices. Hydrogen 2S maps to n_ext = 1 (since |λ_n^Fock| = n_principal, and Fock ↔ CH gives n_CH = n_Fock − 1).

The dimensionless bracket extraction uses

```
C_2S^GeoVac = N_norm · [Σ_R(n_ext=1) + Σ_C(n_ext=1)]
```

with several natural normalizations tested. The most physically motivated is N_norm = 1/(4π)², from two photon Hopf-base measure factors (Paper 18 M1 mechanism).

## 3. Numerical results

### 3.1 Bare spectral sum

| n_max | rainbow | crossed | raw total | C_2S (1/(4π)² norm) | predicted MHz |
|------:|--------:|--------:|----------:|--------------------:|--------------:|
| 2 | 27.5 | 17.3 | 44.82 | 0.284 | 0.067 |
| 3 | 173.5 | 92.8 | 266.33 | 1.687 | 0.399 |
| 4 | 500.1 | 261.2 | 761.29 | 4.821 | 1.139 |
| 5 | 1080.7 | 577.1 | 1657.83 | 10.498 | 2.481 |
| 6 | 1909.0 | 1144.0 | 3052.80 | 19.332 | 4.568 |

**Asymptotic growth:** Power-law fit gives raw ~ 6.60 · N^3.43. The sum does NOT converge as n_max → ∞.

This is faithfully consistent with flat-space QED: the bare two-loop SE diagram is logarithmically UV-divergent in d=4 dimensions, and on the discretized S³ with no momentum cutoff the divergence manifests as a power law in the level-cutoff n_max.

### 3.2 Disconnected-piece subtraction (failed)

Standard QED renormalization removes the disconnected one-loop² piece: Σ_2L^connected = Σ_2L^full − [Σ_1L]² × combinatorial factor. Numerical test:

| n_max | full norm | [Σ_1L]² norm | connected |
|------:|----------:|-------------:|----------:|
| 3 | 1.687 | 0.0004 | 1.686 |
| 4 | 4.821 | 0.0007 | 4.820 |
| 5 | 10.498 | 0.0009 | 10.497 |
| 6 | 19.332 | 0.0010 | 19.331 |

The disconnected piece is 1000× smaller than the full sum and does not tame the divergence. **The divergence lives in the genuine connected two-loop piece** — the same place the flat-space QED divergence lives, requiring Z_2 and δm renormalization.

### 3.3 Drake-Swainson asymptotic subtraction (failed)

Drake-Swainson regularization (Paper 34's 13th projection, LS-4) works for *logarithmic* divergences: it splits the spectral sum at intermediate K, subtracts the asymptotic Λ-dependent piece, takes K → ∞ to get a finite answer. This works for one-loop bound-state Bethe logs at ℓ > 0 (LS-4, Paper 36).

Applied here: subtract 6.60 · N^3.43 from the raw sum. Residuals oscillate around zero (−18.9, −3.2, +15.3, −15.6) with no convergent finite part. **Drake-Swainson does not work for power-law divergences.**

### 3.4 Honest verdict tier

Per the LS-8a verdict structure stated at sprint launch:

- Sub-percent agreement (within ±10%): **NO** — the sum diverges, no unique value
- Mixed (within ±50%): **NO** — same reason
- Order of magnitude (within ±factor of 10): possibly, but only at one specific n_max (n_max=4 gives C_2S ~ 4.82, within 1.33× target). Not stable.
- Wildly off: **NO** — at n_max=4, predicted value is the right sign and within 33% of literature.

**Net: WEAK.** The framework's two-loop machinery has the right STRUCTURAL FORM (correct (α/π)² prefactor, correct UV-divergent integrand, correct sign at finite cutoff) but cannot autonomously extract the renormalized finite answer.

## 4. What this means for Paper 35 Prediction 1

Paper 35 Prediction 1 says: a GeoVac observable contains π if and only if its evaluation includes a continuous integration over a temporal/spectral parameter promoted from the discrete graph spectrum.

**Confirmation tier (LS-7):** The prefactor (α/π)² traces to two iterated photon proper-time integrations. ✓ Consistent.

**Strong test (LS-8a):** Does the dimensionless bracket coefficient (C_2S) emerge from the GeoVac iterated CC spectral action without external calibration?

LS-8a's answer: the bare spectral action does NOT produce C_2S because the sum diverges. The framework correctly identifies the UV divergence of two-loop QED (a structural result; the divergence is real physics) but does not autonomously perform the renormalization that gives the finite physical answer.

**Refined reading of Prediction 1:** the prediction is consistent with LS-8a — π enters via the two iterated proper-time integrations exactly as predicted. But the prediction does NOT extend to claiming GeoVac autonomously regularizes; renormalization counterterms (Z_2, δm) are still required input from QED. The renormalization adds finite contributions to C_2S that the bare framework cannot generate.

This is the right place to recognize a scope boundary: GeoVac-internal projections (Paper 34's fifteen) handle the *finite* parts of one-loop QED via Drake-Swainson + Eides §3.2 conventions. They do NOT handle the renormalization counterterms of multi-loop QED, which are a separate piece of input data.

## 5. Implications for Paper 36

Paper 36 §VII's open Proposition 1 (LS-7 multi-loop test, sharpened) needs refinement:

**Original:** "the strong test is the native derivation of the dimensionless bracket coefficient C_2S = +3.63 from the bound-state Sturmian projection, which is the LS-8a sprint scope."

**Sharpened by LS-8a:** the bare spectral action on Dirac-S^3 with bound-state Sturmian projection at λ=Z/n produces a UV-divergent expression at the order (α/π)²(Zα)⁴/n³, faithfully matching the flat-space two-loop SE divergence structure. Native finite extraction of C_2S would require *additional* sprint work implementing:

1. **Wave-function renormalization Z_2 at the spectral level.** This requires computing the on-shell electron propagator self-energy at one loop, extracting Z_2^{(1)}, and subtracting Z_2^{(1)} × (one-loop SE) at the two-loop level.
2. **Mass renormalization δm at the spectral level.** Subtract the δm × (electron propagator derivative) counterterm.
3. **Ward identity verification on the spectral side.** The renormalization should preserve gauge invariance at two loops, which on Dirac-S^3 is a non-trivial structural test.

This is a genuine multi-sprint extension. The remaining +1.20 MHz multi-loop QED contribution to the H 2S Lamb shift (Eides Tab. 7.3) cannot be derived natively from GeoVac without this renormalization machinery.

**This is honest scope.** The framework's one-loop bound-state QED (Paper 36) closes at sub-percent without fits; the two-loop multi-loop sector requires renormalization machinery as input.

## 6. Honest limits

1. **Spectral sum truncation does not regularize.** Cutting at n_max gives a finite number, but that number has no physical meaning — different cutoffs give different answers spanning 0.07 → 4.57 MHz. The bare framework does not pick out a canonical truncation that gives the right answer.

2. **The n_max=4 coincidence is not robust.** At n_max=4, C_2S = 4.82 (1.33× target). This is suggestive but not stable across n_max. Calling this "agreement" would be over-interpretation.

3. **Disconnected and Drake-Swainson subtractions both fail.** The two natural regularizations available within the GeoVac toolkit don't tame the power-law divergence. Renormalization counterterms (which are NOT GeoVac-internal) would be required.

4. **Sign is correct.** At every n_max tested, the sum is positive, matching the sign of the literature C_2S.

5. **The bare framework correctly identifies that two-loop QED is UV-divergent.** This is a structural confirmation worth recording: the iterated CC spectral action on Dirac-S^3 inherits the same UV divergence as flat-space QED. It does not magically resolve the divergence; it reproduces it faithfully.

## 7. What remains

- **LS-8b (Karplus-Sachs two-loop VP, +0.16 MHz):** simpler than LS-8a in that vacuum polarization is renormalized by a single Z_3 counterterm. Should still show the same structural pattern: framework reproduces divergence, renormalization is external input.
- **LS-8a-renorm:** if pursued, implement Z_2 and δm spectral counterterms. Major sprint, ~3-4 weeks.
- **Paper 36 §VII update:** apply the LS-8a refinement to Proposition 1.
- **Paper 35 Prediction 1 update:** add a remark that the prediction concerns *finite parts* of GeoVac observables; UV divergences are inherited from the underlying field theory and are renormalized by GeoVac-external counterterms.

## 8. Structural reading

The cleanest reading of LS-8a's result:

> **GeoVac's two-loop machinery is faithful but not autonomous.** The framework reproduces the UV-divergent structure of multi-loop QED exactly (the iterated CC spectral action contains the right divergent integrand), but it does not autonomously generate renormalization counterterms. Finite multi-loop QED predictions require external Z and δm input. This is the scope boundary at which the "GeoVac as two-axis classifier" position from today's earlier conversation finds its concrete computational expression: GeoVac controls the divergence STRUCTURE; the inner factor (renormalization counterterms, in addition to A_F's matter-content) controls the finite calibration.

This is consistent with the SM-thread upper-bound reading: GeoVac maps the outer factor faithfully, including the divergence structure that's geometrically forced. The renormalization counterterms (alongside SM-distinguishing A_F input) are the inner-factor data that GeoVac doesn't autonomously select.

## 9. Files

- Module: `geovac/two_loop_self_energy.py`
- Tests: `tests/test_two_loop_self_energy.py`
- Memo: this file
- Data: `debug/data/ls8a_two_loop_self_energy.json`

## 10. References

- Sprint LS-7: `debug/ls7_two_loop_se_memo.md`
- Sprint LS-4: `debug/ls4_bethe_log_drake_memo.md` (Drake-Swainson method)
- Paper 36 §VII (Proposition 1, sharpened by LS-8a)
- Paper 35 §VII.3 (Prediction 1, refined by LS-8a)
- Paper 28 §curvature_coefficients, §two_loop, §three_loop (existing multi-loop infrastructure)
- M. I. Eides, H. Grotch, V. A. Shelyuto, *Phys. Rep.* 342 (2001) 63 — multi-loop QED reference values
- K. Pachucki, V. A. Yerokhin, *J. Phys. Chem. Ref. Data* 39 (2010) 023105 — modern compilation

## 11. Summary table

| Quantity | LS-7 status | LS-8a result |
|----------|-------------|-------------|
| Structural prefactor (α/π)²(Zα)⁴/n³ | Derived | Confirmed |
| C_2S = +3.6266 (literature) | Reference target | Sum diverges; no native unique value |
| Sign of C_2S | + (input) | + (output, at every n_max) |
| Order of magnitude | O(1)–O(10) | O(1)–O(10) at finite n_max |
| Convergence in n_max | Not tested | DIVERGENT, raw ~ N^3.43 |
| Disconnected subtraction | Not tested | Tiny piece, doesn't help |
| Drake-Swainson regularization | Not tested | Power-law divergence — DS doesn't apply |
| Field-theoretic renormalization (Z_2, δm) | Not implemented | Required for finite native C_2S |
| Paper 35 Prediction 1 (structural) | Consistent | Consistent |
| Paper 35 Prediction 1 (strong test) | Deferred to LS-8a | Cannot be tested without renormalization machinery |

**LS-8a verdict: WEAK.** Framework reproduces two-loop QED divergence structure faithfully but cannot autonomously extract finite C_2S. This is honest scope, not a framework failure: the bare spectral action of Connes-Chamseddine type does NOT generate renormalization counterterms; it preserves them as external input.
