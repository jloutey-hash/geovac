# NA-1 Reading A vs B — period-value irreducibility is NOT a discriminator

**Date:** 2026-06-14. **Status:** clarifying negative. **Reading A stands** (Sprint JLO-Depth2, 2026-06-07, already in Paper 56 §sec:open_na1).

## What was attempted
A fresh attempt to decide Reading A (abelianization) vs Reading B (free
non-abelian / shuffle) for GeoVac's cosmic-Galois substrate, using the
**period-value content** of the k=3 self-energy substrate (S^(3)) rather than
the JLO cyclic cocycle. Drivers: `debug/sprint_na1_s3_bracket_coproduct.py`,
`debug/sprint_na1_adversarial_check.py`. Data in `debug/data/na1_*.json`.

The argument: GeoVac's S^(3) closed form contains the irreducible depth-2 MZV
ζ(5,3) (gated 1.15e-198, frozen `tests/test_s3_w10_identification.py`).
ζ(5,3) is the period dual to the Lie bracket [f₃,f₅]; the abelianization kills
brackets; so (claimed) GeoVac carries the bracket → NOT-A.

## Why it FAILS (the invalid inference)
**Irreducibility ≠ bracket.** PSLQ confirms ζ(5,3) is irreducible to *products*
of single zetas (verified to dps 120, maxcoeff 10⁹; literature-solid). But
"irreducible / not reducible to depth-1" does NOT imply "non-abelian / bracket":

- By **Cartier–Milnor–Moore**, a cocommutative (abelian-dual) Hopf algebra has
  genuinely new **primitive** generators at *every* degree. A new irreducible
  generator at weight 8 is consistent with an **abelian** substrate — it is a
  fresh primitive atom, not a bracket.
- The discriminator is the **coproduct** (primitive Δ′=0 vs deconcatenation
  Δ′=ζ(5)⊗ζ(3)), NOT period-value irreducibility. Period-value irreducibility
  separates "generator" from "product"; it does **not** separate "primitive
  generator" (A) from "bracket" (B).
- Sprint JLO-Depth2 §4.4 **explicitly warned about this exact error**: depth-2
  content not reducible to depth-1 "looks like a Reading B signal at face value
  … [but] is a structural fact about commutative algebras." The "no weight-8
  Mellin slot ⇒ can't be a fresh atom" closure in the adversarial driver is
  WRONG: loop order indexes new primitives; loop-order-3 gives a new primitive
  whose value is ζ(5,3), substrate stays abelian.

## The correct measurement (already done)
**JLO-Depth2 (2026-06-07)** measured the *coproduct* directly via the depth-2
JLO cyclic cocycle on the CH spectral triple: bit-exactly **S₃-symmetric /
cocommutative** (40/40 zero swap-asymmetric residuals, 6/6 S₃-orbit identical,
n_max=2 and 3), strictly stronger than commutativity-of-𝒜 forces. → **Reading
A**. GeoVac's substrate is abelian; the injection U*_GV ↪ 𝒰₄^ab is the
structurally correct, complete statement; shuffle-Hopf T(V) enrichment is NOT
empirically demanded.

## Reconciliation (both facts are true, no contradiction)
GeoVac's *observables* land on bracket-VALUED motivic periods (ζ(5,3)) — that
is real. GeoVac's *substrate* is abelian — also real. The period value is an
element of the TARGET ring 𝒢₄ (which has brackets); the substrate coproduct on
the generator producing it is primitive. Producing a bracket-valued number ≠
carrying a bracket in the substrate. This is exactly the injection's content.

## If the question is ever reopened
The correct next probe (per JLO-Depth2 §6.1) is the **depth-3 JLO cocycle**
(S₄-cocommutativity) on **S^(4)** — a coproduct-level test, NOT a period-value
test, and NOT the shuffle enrichment build. Period-value irreducibility is
retired as a discriminator.

## Discipline
Curve-fit audit clean (PSLQ relations exact / null, zero free params). The
ERROR was interpretive, not numerical: the computations (ζ(5,3) irreducible)
are correct; the inference (⇒ Reading B) is invalid. Caught by the
authoritative-source / consistency gate before any paper edit. No corpus edit
made; Paper 56 §sec:open_na1 Reading A is correct and unchanged.
