# Sprint A — Paper 24 π³ cross-check for α-EB v2

**Scope.** Test whether π³ = Vol(S⁵) is a structurally natural prefactor
in Paper 24's Bargmann-Segal / S⁵ spectral-action machinery, or whether
the α-EB v2 coincidence `R_predict = K − 1/α − α² ≈ π³α³` (0.25% at 80 dps)
is numerical only.

**One-line verdict.** **PARTIAL (negative-leaning).** π³ = Vol(S⁵)
appears in **every** Seeley-DeWitt coefficient on round S⁵ as the
integration-measure factor, so the π³ prefactor itself is
structurally unsurprising once an S⁵ spectral-action object is
admitted. **But no S⁵ spectral-action coefficient produces a
contribution at the right order (α³ against leading 1/α) to
structurally derive the α³ multiplier.** The α-EB v2 finding should
**not** motivate a Paper 2 §IV.G integration beyond a structural-hint
footnote.

---

## 1. Paper 24 π-free vs heat-kernel π reconciliation

Paper 24 Theorem 1 states that the Bargmann-Segal graph is
bit-exactly π-free in exact rational arithmetic at every finite
N_max: every diagonal entry, every adjacency weight, and the full
spectrum ℏω(N + 3/2) are rational. Paper 24 §III repeats this:
"the single π in the entire construction appears in the Gaussian
normalization π⁻³ of the continuous Bargmann measure, which is
never used by the discrete graph."

The π³ candidate here **does not violate this claim.** Paper 24's
π-free certificate is scoped to *three* objects:

1. the node diagonal (HO spectrum),
2. the adjacency edge weights (Wigner 3j × radial matrix elements),
3. the spectrum of D − A and its Hodge-decomposition partners.

All three are exact rationals. Vol(S⁵) = π³ is a **different**
object: it is the integral of the volume form over the *continuous*
S⁵ manifold. This appears in heat-kernel traces Tr exp(−tP), in
Seeley-DeWitt coefficients a_k, and in the Connes-Chamseddine
spectral-action asymptotic expansion. These continuum objects
co-exist with the π-free discrete graph without contradiction,
because they are integrated invariants of the **ambient** manifold,
not invariants of the discrete graph itself.

**Reconciliation rule:** the Paper 24 π-free claim applies to the
discrete spectrum and graph. The π³ candidate lives in the
integration-measure of the ambient continuous manifold. Both can
hold simultaneously, and any structural reading must be stated in
the continuum spectral-action register, not the discrete
graph register. This distinction is respected throughout below.

---

## 2. Seeley-DeWitt a₀, a₂, a₄ on round S⁵

Exact sympy computation (script:
`debug/alpha_sprint_a/compute_s5_sd_coeffs.py`, data:
`debug/data/alpha_sprint_a_s5_sd_coeffs.json`). All values on
**unit** round S⁵ (R = 1). Conventions match
`geovac/qed_vacuum_polarization.py` on S³ (dim_S = 4 spinor,
Lichnerowicz endomorphism E = R_scalar/4 for the squared Dirac
operator; 4π prefactor kept outside a_k).

**Volume and curvature invariants.** Vol(S⁵_{R=1}) = 2π³ / Γ(3) = π³
(exact). On constant-curvature S⁵: R_scalar = 20, |Ric|² = 80,
|Riem|² = 40.

**Scalar Laplacian.**
- a₀(scalar) = π³
- a₂(scalar) = (10/3) π³ = (1/6)·R_scalar·Vol = (20/6)·π³
- a₄(scalar) = (16/3) π³ = (1/360)·(5·400 − 2·80 + 2·40)·π³

**Squared Dirac D² with Lichnerowicz E = R_scalar/4.**
- a₀(Dirac) = 4 π³ (dim_S = 4 × Vol)
- a₂(Dirac) = −(20/3) π³ = 4·(R_scalar/6 − R_scalar/4)·Vol
  (the Lichnerowicz shift makes this negative)
- a₄(Dirac) = (14/3) π³ = 4·(1/360)·[5·R² − 2·|Ric|² + 2·|Riem|²
  − 30·R·E + 60·E²]·Vol, with R·E simplification

**π-content per coefficient.** Every single coefficient factorizes
as π³ × (rational). No independent transcendental is injected at
any order. This is the Paper 24 analog of the GeoVac S³ module
finding (`qed_vacuum_polarization.py`) that a₀,a₁,a₂ on unit S³
are √π, √π, √π/8 respectively — pure rational × π^{d/2} with d
fixing the pi-power.

**Key observation.** On a 5-manifold d = 5 is odd, so
(4π)^{−d/2} = (4π)^{−5/2} carries an explicit √π factor outside
a_k, and every a_k carries π³ from Vol(S⁵) inside. **π³ = Vol(S⁵)
is therefore the natural π-prefactor of every scalar or spinor SD
coefficient on S⁵.** This is a clean, shape-level positive.

---

## 3. Does π³ enter at order α³ structurally?

On a 5-manifold the Connes-Chamseddine asymptotic expansion of the
spectral action is
```
Tr f(D/Λ) ~ f_5 Λ⁵ a_0 + f_3 Λ³ a_2 + f_1 Λ a_4 + O(Λ⁻¹).
```
All powers of Λ are odd (d = 5 has no log term; that would require
d even). Each SD coefficient a_k is a rational multiple of π³
= Vol(S⁵).

For π³ to enter at order α³ against the leading 1/α of K, one
would need:

**Requirement (a).** A coupling Λ ↔ 1/α somewhere. In standard CC
on a 4-manifold, Λ is a fixed UV cutoff (unification scale); the
gauge coupling g² ≈ α is determined by
`g² = (12 / 4!) · a_0 / a_4` (Chamseddine-Connes 1996 Eq. 3.18
on d = 4), which expresses 1/α in terms of SD coefficients, NOT α³.

**Requirement (b).** A natural mechanism for a specific SD
coefficient to couple at the 3rd power of α. On a 5-manifold the
only α's in the game enter through the gauge sector via the
minimally coupled Dirac operator D + iA; these give at most α²
(through the F_μν F^μν piece in a_4 at d = 4, or equivalents on
d = 5). α³ contributions in CC are loop-level (two-loop QED),
not tree-level spectral-action terms.

**Neither is available in Paper 24.** The Bargmann-Segal
construction is scalar and first-order; it does not include a
coupled Dirac operator on S⁵ (and indeed the Paper 24 asymmetry
theorem says the natural S⁵ operator is first-order and
complex-analytic, not second-order Riemannian — so even the
squared-Dirac pipeline above is not the native S⁵ operator
of Paper 24, it is the S³-module analog extended here only as a
sanity reference).

**Numerical sanity.** R_predict = 1.208e−5 = π³ × α³ × 1.0025.
α³ ≈ 3.9e−7; π³ ≈ 31.006; π³·α³ ≈ 1.205e−5. Rational multipliers
on SD coefficients are O(1) to O(10) (10/3, 16/3, 20/3, 14/3).
Obtaining α³ itself would require cubing the quantity that
produces α (which in CC comes from a single a_0/a_4 ratio), not a
new SD coefficient. There is no shape-level reason for a₀·α, a₂·α,
or a₄·α to pick up α³ from S⁵ structure alone.

---

## 4. Verdict

**PARTIAL (negative-leaning).**

- The **π³ prefactor** in R_predict has a clean structural home on
  S⁵: every Seeley-DeWitt coefficient on round S⁵ is π³ × rational.
  In this limited sense, π³ = Vol(S⁵) is a *shape-level* fit.
- The **α³ multiplier** has **no** structural home on S⁵. Standard
  CC does not produce α³ coefficients in the spectral-action
  expansion at tree level on 5-manifolds; Paper 24 provides no
  additional mechanism.
- The 0.25% gap (R_predict / (π³ α³) = 1.0025) is currently
  unexplained. The next-order coefficient C ≈ 0.344 (memo v2) is
  not uniquely determined by a single CODATA data point.
- **Non-falsifying numerical coincidence.** The find is consistent
  with WH5 (α is a projection constant, not derivable) and with
  WH1's Sprint A findings (shape-match only at the APS-flavor
  level, no explicit derivation). It does not *refute* any Paper
  24 or Paper 2 claim.

---

## 5. Recommendations

**For Paper 2 §IV.G.** **Defer the α-EB v2 integration.**
Specifically:

- Do NOT add "R_predict = π³ α³" as a substantive finding in §IV.G
  or the Open Questions register, because the α³ multiplier is not
  derived.
- DO record the observation as a footnote or a tightly scoped
  remark in Paper 2 §IV.G open questions, with the caveat *"π³
  equals Vol(S⁵), but no S⁵ spectral-action coefficient produces
  an α³ contribution; the numerical coincidence at 0.25% therefore
  remains structurally unexplained."*
- Keep Paper 2 conjectural (no change to that status).
- Do not list this under the Phase 4B–4I positive structural
  identifications (κ↔B, F = D_{n²}(d_max), Δ⁻¹ = g_3^Dirac, ζ(3)
  in Dirac). Those are closed algebraic identities; this is a
  numerical match at the 3rd-decimal level.

**For Paper 24.** **No changes.** The π-free certificate (Theorem 1)
is intact; the continuum Vol(S⁵) = π³ does not enter any discrete
Bargmann-Segal object. The S⁵ SD coefficients computed here are
standard textbook values on constant-curvature S⁵ (Vassilevich,
Branson-Gilkey), not new Paper 24 material.

**For WH5 (§1.7 register).** Update the WH5 status line to note
that the α-EB v2 π³ observation **survives partial-shape cross-check
on S⁵ SD coefficients but does not upgrade to derivation status**:
π³ is the Vol(S⁵) integration-measure factor on every S⁵ SD
coefficient (shape consistent), but α³ is not produced by any
Paper-24 / CC spectral-action coefficient at the matching order
(derivation absent). Net: the post-cubic residual remains a
numerical coincidence with structural-hint status, not a
derivation. This is consistent with WH5's core thesis.

**For next sprint target.** If pursued further, probe **loop-level**
CC contributions, not tree-level SD coefficients on S⁵: does any
two-loop CC contribution produce a relative α³·Vol(S⁵) shift to 1/α?
This is outside Paper 24 and outside the current Phase 4B–4I scope.

---

**Files.**
- `debug/alpha_sprint_a/compute_s5_sd_coeffs.py` (this script)
- `debug/data/alpha_sprint_a_s5_sd_coeffs.json` (numerical record)
- `debug/alpha_sprint_a/paper24_pi_cubed_crosscheck_memo.md` (this memo)

**Cross-reference.**
- α-EB v2 memo: `debug/alpha_sprint_a/error_budget_memo_v2.md`
- Paper 24: `papers/core/paper_24_bargmann_segal.tex`
- S³ SD module: `geovac/qed_vacuum_polarization.py`
- Paper 2 §IV: `papers/conjectures/paper_2_alpha.tex`, lines 283–666
- CLAUDE.md §1.7 WH1/WH5 (Sprint A register entries)
