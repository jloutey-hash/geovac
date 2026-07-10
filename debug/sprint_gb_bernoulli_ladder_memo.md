# Sprint GB — the Bernoulli ladder and the functional-equation reading of the two-layer split

**Date:** 2026-05-29. **Type:** diagnostic-only (no production code, no paper edits applied — paper capture recommended, see §5). **Verdict:** STRUCTURAL-CORRESPONDENCE (not generic; not a derived theorem).

## 1. The question

Motivated by the -1/12 conversation: the gravity sector keeps producing the same small rationals (conical tip 1/12, replica-entropy derivative 1/6, per-t UV target 1/24, scalar Casimir 1/240), and the α-conjecture carries F = π²/6 = ζ(2). Are these shadows of a single Bernoulli ladder B_{2n}, glued by the ζ functional equation — i.e. is the framework's **two-layer split (Paper 34/35)** the reflection ζ(s) ↔ ζ(1-s)?

**Null (generic):** every ζ(2n) and ζ(1-2n) carries a Bernoulli by the standard formulas; finding B₂ in the lowest gravity coefficient and lowest α-coefficient is *forced* by lowest-order, carries no GeoVac information.

## 2. Exact results (sympy, skeleton quantities)

The ladder, two windows per rung:

| rung | M2 / transcendental window | skeleton / rational window |
|:----:|:--------------------------:|:--------------------------:|
| B₂ = 1/6   | ζ(2) = π²/6   | ζ(-1) = -1/12 |
| B₄ = -1/30  | ζ(4) = π⁴/90  | ζ(-3) = 1/120 |
| B₆ = 1/42   | ζ(6) = π⁶/945 | ζ(-5) = -1/252 |

GeoVac coefficients placed on the ladder (all exact, FIT = sympy-zero residual):

| coefficient | source | ladder placement |
|:------------|:-------|:-----------------|
| conical tip 1/12 | G4-2 | -ζ(-1) = **B₂/2** (skeleton) |
| replica entropy derivative 1/6 | G4-4f | -2ζ(-1) = **B₂** (skeleton) |
| per-t UV target 1/24 | v3.20.0 | -ζ(-1)/2 = **B₂/4** = string c/24 (skeleton) |
| scalar Casimir S³ 1/240 | Paper 35 | +ζ(-3)/2 = **B₄ rung** (skeleton) |
| α-conj F = π²/6 | Paper 2 | π²·B₂ = ζ(2) = **B₂, M2 window** |

**Controls behave:**
- κ = -1/16 (Fock Jacobian 1/Ω⁴) — no clean ζ form. Correctly off-ladder (different mechanism).
- Dirac Casimir 17/480 — residual 1/32 from the scalar B₄ rung; needs the half-integer Hurwitz shift (ζ(-3,1/2) = -7/960). Correctly off the *scalar* rung; lives in the half-integer sector = **scalar/spinor = M2/M3 split**.

## 3. The non-trivial internal check

The B₂ rung is not three free coincidences. The replica-entropy derivative is *forced* by the tip:
d/dα [(1/12)(1/α - α)] |_{α=1} = -1/6 = -B₂ = 2ζ(-1).
The entropy (B₂) is the α-derivative of the tip coefficient (B₂/2). Exact. This is a structural consistency check, not a free fit.

## 4. Audit (per feedback_audit_numerical_claims)

- **Free params:** zero; all coefficients independently pre-derived; fits are exact equalities.
- **Selection bias:** controlled. κ off-ladder (right reason); Dirac off scalar-rung (right reason). Not "everything fits."
- **Alternatives:** generic null defeated by (a) B₄ rung *independently* populated by a *different* observable (Casimir, not horizon); (b) replica derivative internally forced.
- **Robustness:** exact, not numerical.
- **Independent test:** B₄ rung is the independent test of the B₂ observation — lands. Spinor off-rung is a second independent prediction — lands.

**Verdict: NOT generic.** But it is a **structural correspondence / reading**, at the same epistemic level as the Paper 32 §VIII case-exhaustion theorem's "structural correspondence" — NOT a forcing proof.

## 5. The new content

We discovered the two-layer split (Layer 1 π-free / Layer 2 π-bearing) empirically via the bit-exact mapping; we separately have the master Mellin engine M2 ring. This sprint **unifies them under the ζ functional equation**:

> **Layer 2 observables carry ζ(2n) = rational·π^{2n} (M2 transcendental ring). Layer 1 skeleton coefficients carry ζ(1-2n) = rational (π-free residue). The functional equation is the projection map. The scalar/spinor distinction is integer/half-integer Hurwitz shift = M2/M3.** The master Mellin engine M2 = 𝓜[Tr(D²e^{-tD²})] is the Mellin object whose positive-even-integer values are the Layer-2 transcendentals and whose analytic continuation to negative-odd integers gives the Layer-1 gravity coefficients.

This is the **rational-residue extension of the case-exhaustion theorem**: the theorem already says all π is 𝓜[Tr(D^k e^{-tD²})]; this says the rational coefficients are the negative-integer continuation of the *same* Mellin object.

The -1/12 in the black-hole conical tip is ζ(-1), the B₂ skeleton residue — genuinely the analytic continuation of 1+2+3+…, sitting in the Layer-1 window of the same B₂ whose Layer-2 window is the F = π²/6 of the α conjecture.

## 6. Honest scope / what this is NOT

- A *reading*, not a theorem. "Two-layer split IS the functional equation" is a correspondence supported by the ladder, not a derivation.
- Only B₂ and B₄ rungs are populated by current data; B₆ rung is empty (no GeoVac observable identified at ζ(6)/ζ(-5) yet — a falsifiable target: if the framework forces a ζ(-5)=-1/252 skeleton coefficient somewhere, it would be the third rung).
- The spinor sector (17/480) needs the M3 half-integer machinery; the ladder as stated is the scalar (M2) ladder. A parallel half-integer Hurwitz ladder for the spinor sector is the natural follow-on.

## 7. Recommended paper capture (NOT applied — flagged for PI)

Natural home: **Paper 32 §VIII** (case-exhaustion theorem) as a Remark — "rational-residue extension via the functional equation" — OR **Paper 35** (two-layer split) as a subsection. Recommend Paper 32 §VIII because the case-exhaustion theorem is the precise object this extends. Interpretive/framing content (§1.5-sensitive: it's a "reading"), so flagged rather than auto-applied.

## Files
- `debug/sprint_gb_bernoulli_ladder.py` (first pass, sign typo + crash)
- `debug/sprint_gb_ladder_final.py` (corrected, clean)
