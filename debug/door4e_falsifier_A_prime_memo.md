# Door 4e — Falsifier A' on DAS: does the AC tensor product force ℍ?

*Tests the literal Falsifier A' flagged in Door 4d: do any natural AC
tensor-product compatibility conditions distinguish ℍ from M₂(ℂ) at the
n=2 rung? If yes, DAS upgrades from PARTIAL-DOOR to FULL-DOOR. If no
(expected), DAS is CONFIRMED as a genuine new principle requiring a
structural axiom upgrade (paper-level), not derivable from existing CCM.*
*Structure-only / compatibility probe. NO Yukawa value selected
(guardrail).*
Started + closed 2026-06-02. Driver: `debug/door4e_falsifier_A_prime.py`.
Data: `debug/data/door4e_falsifier_A_prime.json`.

---

## TL;DR verdict — **PARTIAL-DOOR (CONFIRMED)**

The literal Falsifier A' is **NEGATIVE on all three candidate
compatibility conditions** — exactly as Door 4d flagged would likely
happen. The standard CCM AC tensor-product construction does **not**
provide any morphism that forces the n=2 inner algebra to be ℍ rather
than M₂(ℂ).

The substantive content of Door 4e is **not** the negative outcome
(which was expected) but the **input-comparison** it enables:

> **Standard CCM ℍ-selection** requires THREE axioms (second-order
> condition + complex chiral fermion rep + 2N² = 32 dimension count).
> **DAS ℍ-selection** requires ONE axiom (the rung-n sphere's
> division-algebra structure is inherited by the inner algebra).

Both are imports relative to the bare Connes axioms. **DAS is leaner.**
The path to FULL-DOOR is to **elevate DAS's single axiom — the
sphere-Lie-group correspondence — to a foundational construction
principle** ("the inner algebra at rung n is the *-algebra whose unit
group IS the rung-n sphere as a Lie group, when one exists; else
fallback M_n(ℂ)"). This is a paper-level structural-axiom upgrade,
not a theorem about existing CCM machinery.

---

## What was tested

Three candidate compatibility conditions between the inner algebra and
the AC tensor-product structure, each natural and each independently
plausible as a forcing handle:

### C1. SU(2)-equivariance under Ad action

SU(2) acts on M₂(ℂ) by Ad (conjugation): g · m · g⁻¹. A real
*-subalgebra A ⊂ M₂(ℂ) is **SU(2)-Ad-invariant** if it is closed under
this conjugation.

**Result (bit-exact, 64 SU(2) × ℍ samples):** ℍ is SU(2)-Ad-invariant
(0 violations), and M₂(ℂ) is trivially SU(2)-Ad-invariant (the whole
space). Both are SU(2)-Ad-invariant real *-subalgebras of M₂(ℂ).

**Distinguishing power:** **NONE**. The MINIMAL non-trivial SU(2)-Ad-invariant
real *-subalgebra of M₂(ℂ) containing I is **ℍ**; the MAXIMAL is **M₂(ℂ)**.
Selecting between them requires a MINIMALITY criterion that is itself
extra input.

### C2. Hopf U(1)-equivariance (scalar phase)

The Hopf bundle S³ → S² has structure group U(1) acting on S³ ⊂ C² by
right scalar phase: (z₁, z₂) → (z₁ e^{iθ}, z₂ e^{iθ}) = e^{iθ}·v.

**Result (bit-exact, 64 random M ∈ M₂(ℂ) samples):** every 2×2 complex
matrix commutes with the scalar phase e^{iθ}I bit-exactly (max
commutator norm = 0). Scalar phases are central; the U(1) action is
trivially compatible with **every** matrix action.

**Distinguishing power:** **NONE**. The condition is vacuous.

### C3. Principal-bundle compatibility (Hopf U(1) fiber)

The Hopf bundle S³ → S² is a principal U(1)-bundle. An inner algebra A
acting on S³ ⊂ C² by left matrix multiplication is **principal-bundle
compatible** if its action commutes with the right U(1) fiber action.

**Result (bit-exact, 32 random samples):** ℍ-action and M₂(ℂ)-action
both satisfy h · (q · e^{iθ}) = (h · q) · e^{iθ} at machine precision
(ℍ: max error 1.88e−15; M₂(ℂ): max error 1.34e−15). Both are
principal-bundle compatible by associativity of matrix multiplication
and centrality of scalar phases.

**Distinguishing power:** **NONE**.

### Net Falsifier A' verdict

| Compatibility condition | ℍ passes? | M₂(ℂ) passes? | Distinguishes? |
|:-----------------------|:--------:|:-------------:|:--------------:|
| C1: SU(2)-Ad-invariance | ✓ (min) | ✓ (max) | **NO** |
| C2: Hopf U(1)-equivariance | ✓ | ✓ (trivially) | **NO** |
| C3: Principal-bundle compatibility | ✓ | ✓ | **NO** |

**The literal Falsifier A' is NEGATIVE.** No natural AC tensor-product
compatibility condition forces ℍ over M₂(ℂ). The standard CCM
construction is genuinely too permissive at the n=2 rung.

---

## Why this is the SUBSTANTIVE content of Door 4e

The negative outcome above was the expected one (Door 4d flagged it).
What's substantive is what the negative outcome tells us about the
**comparison between CCM's and DAS's ℍ-selection imports**.

### Standard CCM (Chamseddine–Connes 2008) ℍ-selection

Three independent axioms beyond the bare Connes spectral-triple data:

1. **Second-order condition:** [[D, a], [J b J⁻¹, c]] = 0 for all
   a, b, c ∈ A_F, with **nonzero** D_F (Yukawa-carrying finite Dirac).
   Together with the action of A_F on the chiral fermion sector, this
   constrains the L-block algebra to be quaternionic.
2. **Complex chiral fermion representation:** the L-doublet must be
   paired with the anti-L-doublet via charge conjugation (i.e., the
   2-dim complex rep is realized as ℍ acting on its own underlying C²,
   not as two independent 2-dim complex reps).
3. **Dimension count 2N² = 32:** the total Hilbert-space dimension fixes
   the rank of the finite triple, sub-constraining A_F to lie inside
   M₂(ℍ) ⊕ M₄(ℂ).

The combination selects ℂ ⊕ ℍ ⊕ M₃(ℂ).

### DAS (this work, Door 4d) ℍ-selection

One axiom:

> **DAS:** the inner algebra at rung n is the natural associative real
> normed division algebra whose unit-norm sphere IS the rung-n Hopf
> bundle total space S^(2n−1). When no such division algebra exists
> (Hurwitz's theorem at n=3: no associative real normed division
> algebra has dim 6), fall back to the matrix algebra M_n(ℂ) reproducing
> SU(n).

The combination selects ℂ ⊕ ℍ ⊕ M₃(ℂ).

### Why this comparison matters

Both derivations land on the SM finite algebra. Neither is derivable
from the other; they are **alternative ℍ-selectors**, not redundant.

- **CCM is more standard** (integrates with existing NCG machinery; the
  three axioms are all well-established Connes-Marcolli content).
- **DAS is leaner** (one geometric criterion, reading the rung-n sphere's
  Lie-group structure directly off the Bertrand × Hopf-tower output
  that GeoVac already extracts for the gauge group).

By **input minimality**, DAS wins. By **integration with established
literature**, CCM wins. **Neither beats the other across both axes
simultaneously.**

The path to FULL-DOOR for DAS is therefore not a theorem about existing
CCM machinery (Falsifier A' shows that's vacuous) but a paper-level
**structural-axiom upgrade**: adopting DAS's single criterion as a
natural construction principle to replace CCM's three axioms.

---

## The structural-upgrade path for FULL-DOOR

Three candidate axioms could elevate DAS from "consistent extension" to
"natural construction principle." Each agrees with Door 4b at n=1, n=3
and closes the fork at n=2 by selecting ℍ.

### Upgrade A — Minimal-AC axiom

> The inner algebra at rung n is the **minimal** real *-subalgebra of
> M_n(ℂ) that (a) contains the identity, (b) is closed under Ad SU(n),
> and (c) is non-abelian (to support a non-trivial gauge group).

- Selects at n=2: ℍ (the unique 4-real-dim minimal non-abelian
  SU(2)-Ad-invariant *-subalgebra of M₂(ℂ) containing I; the only
  candidates are ℂ·I — trivial, no gauge structure — and ℍ itself,
  with M₂(ℂ) being non-minimal).
- Naturalness: **strong**. Minimality is a standard guiding principle
  in NCG. The principle "use the minimal algebra reproducing the gauge
  group" is a SHARPER version of CCM's "maximal sub-algebra satisfying
  order-one."

### Upgrade B — Sphere-Lie-group axiom (**recommended**)

> When the rung-n Hopf-bundle total sphere S^(2n−1) has a Lie-group
> structure (i.e., for n=1, 2: from the ℂ and ℍ division algebras), the
> inner algebra at rung n is the *-algebra whose group of unitaries IS
> this Lie group. When S^(2n−1) has no Lie-group structure (n=3: S⁵ is
> not a Lie group; n=4: S⁷ is a Lie group only non-associatively via
> octonions), fall back to the minimal matrix algebra M_n(ℂ) reproducing
> SU(n).

- Selects at n=2: ℍ (because S³ = Sp(1) = U(ℍ) is the unique Lie-group
  structure on S³).
- Naturalness: **very strong**. Reads the algebraic content directly off
  the geometric content GeoVac already extracts via Bertrand × Hopf-tower.
  No new input data — just reads the Hopf-rung sphere FULLY (Lie
  structure + topology) where the existing construction reads it
  PARTIALLY (gauge group only).

### Upgrade C — Division-ring axiom

> The inner algebra at rung n is, when possible, a **division ring**
> (every nonzero element invertible). When no division ring of the
> right Wedderburn type exists, fall back to M_n(ℂ).

- Selects at n=2: ℍ (the unique 4-real-dim division ring of
  complex-realizable type).
- Naturalness: **moderate**. "Prefer division rings" is a common
  preference but less structurally tied to the Hopf-rung geometry than
  upgrades A or B.

### Recommendation

**Upgrade B (sphere-Lie-group axiom).** It reads the Hopf-rung sphere's
FULL geometric content rather than just the gauge group, uses NO input
beyond what GeoVac already extracts, and is the leanest extension of
the existing Bertrand × Hopf-tower argument. Upgrades A and C are good
alternative framings but invoke axioms (minimality, division-ring
preference) less directly tied to the Hopf tower.

---

## Inner-algebra status after Door 4e

```
FORCED (Bertrand x Hopf-tower, same principle as gauge content):
  • factor count = 3
  • Wedderburn / matrix-algebra type of each factor
  • n=1 factor = C
  • n=3 factor = M_3(C)

PARTIAL-DOOR (CONFIRMED, Door 4e):
  • n=2 factor = H (vs M_2(C))
      DAS handle (Door 4d): geometric criterion S^(2n-1) = unit norm
        sphere of n-th associative real normed division algebra
      Falsifier A' verdict (Door 4e): the LITERAL test on standard AC
        compatibility conditions (C1: SU(2)-Ad-invariance; C2: Hopf U(1)
        equivariance; C3: principal-bundle compatibility) returns
        NEGATIVE on all three -- no condition forces H over M_2(C) via
        existing CCM tensor-product machinery.
      Substantive content:
        CCM imports H via 3 axioms; DAS imports H via 1 axiom.
        Both are imports; DAS is leaner; neither is derivable from
        the other.
      FULL-DOOR upgrade path (paper-level, sprint-scale):
        Adopt Upgrade B (sphere-Lie-group axiom) as foundational:
        "the inner algebra at rung n is the *-algebra whose unit group
        IS the rung-n sphere as a Lie group, when one exists; else
        fallback M_n(C)." This reads the Hopf-rung sphere FULLY where
        the existing construction reads it PARTIALLY.

FREE (no construction reaches them):
  • generation count N_gen   (invisible to inner automorphisms)
  • inner KO-dimension       (imposed to hit SM total 1 mod 8)
  • Yukawa eigenvalues       (seam theorem; out of scope)
```

Net: **DAS is CONFIRMED as a genuine alternative ℍ-selector, leaner
than CCM's standard three-axiom path.** The PARTIAL-DOOR status of Door
4d is preserved; the next substantive step is a paper-level argument
that adopting Upgrade B is justified relative to (a) CCM's
three-axiom path and (b) the alternative upgrades A and C.

---

## The next sharpest falsifier

**Falsifier A'' (paper-level, sprint-to-month):** Establish that
adopting Upgrade B (sphere-Lie-group axiom) **does not contradict** any
existing GeoVac result, and check whether it **introduces any new
downstream constraints** (e.g., on the inner KO-dim, on the chirality
grading γ_F, on the Higgs sector) that the existing framework either
supports or contradicts. Three outcomes:

1. **Upgrade B introduces a new constraint that the framework
   supports** → DAS gains independent corroboration; promote to working
   axiom alongside (or replacing) CCM's three.
2. **Upgrade B introduces a new constraint that the framework
   contradicts** → Upgrade B is RULED OUT; the fork goes back to
   "admitted-not-forced via combined-J sign table (Door 4c) + leaner
   ℍ-import via DAS (Door 4d) but no FULL-DOOR closure."
3. **Upgrade B introduces no new constraints** (purely a
   reformulation of CCM's content) → DAS is "consistent with CCM but
   with leaner axiom count," PARTIAL-DOOR maintained, no further
   forcing extracted.

Outcome 1 would be the deepest result; outcomes 2 and 3 are also
informative.

**Falsifier B (deferred):** Force the generation count N_gen and the
inner KO-dimension from a packing principle. Same as Door 4b's
Falsifier B; no known handle; deferred.

---

## Honest scope — what was NOT shown

- The literal Falsifier A' was tested on the THREE candidate
  compatibility conditions C1, C2, C3 listed above. A FOURTH compatibility
  condition (e.g., a HIGHER-order Connes axiom not tested here — the
  second-order condition with nonzero D_F is the CCM workhorse and is
  partially tested in Door 4c) might still force ℍ via existing CCM
  machinery. The recorded NEGATIVE for the three tested conditions does
  not preclude a positive for an as-yet-untested fourth — but the
  natural candidates are exhausted.
- No claim is made that adopting Upgrade B is the **right** thing to do
  for GeoVac. The recommendation is structural (Upgrade B is the leanest
  and most geometrically natural); whether to adopt it is a paper-level
  judgment that belongs to the PI and the broader NCG community.
- No Yukawa value, generation count, or KO-dim is selected. The probe
  stays strictly on the algebra-structure side at the n=2 rung. H1 / W3 /
  Koide negatives (CLAUDE.md §3) untouched. The seam theorem (Paper 32
  §VIII Door 4) separating algebra factors from Yukawa values is
  unaffected.

---

## Files / sources

- Driver: `debug/door4e_falsifier_A_prime.py`
- Data: `debug/data/door4e_falsifier_A_prime.json`
- Prior memos:
  - `debug/door4_gauge_yukawa_boundary_memo.md` (seam theorem)
  - `debug/door4b_inner_algebra_forcing_memo.md` (factor count + n=1/3 forced)
  - `debug/door4c_j_signtable_audit_memo.md` (combined-J handle CLOSED NEGATIVE)
  - `debug/door4d_division_algebra_sphere_memo.md` (DAS introduced, PARTIAL-DOOR)
- Paper 32 §VIII.B/C: `papers/group1_operator_algebras/paper_32_spectral_triple.tex`
- Code:
  - `geovac/almost_commutative.py` (`ElectroweakFiniteTriple`, `quaternion_to_matrix`)
  - `geovac/standard_model_triple.py`
  - `geovac/real_structure.py`
- Memory: `bertrand_sm_gauge_truncation.md`,
  `wh_register_april2026.md` (WH4 deflated, "single-input forcing
  statement"; DAS adds inner-factor content)
- Literature:
  - Hurwitz, A. (1898). Classification of normed division algebras.
  - Chamseddine, A. H., & Connes, A. (2008). *Why the Standard Model.*
    J. Geom. Phys. 58, 38. (The CCM ℍ-selector via three-axiom path
    that DAS provides an alternative to.)
  - Connes, A. (1995); van Suijlekom, W. D. (2015). (KO-dim sign tables;
    Connes axiom framework.)
  - Adams, J. F. (1958), "On the non-existence of elements of Hopf
    invariant one." (S^(2n-1) is a Lie group only for n = 1, 2, 4,
    and only for n = 1, 2 with associative multiplication — the deep
    topological / division-algebra fact underwriting Upgrade B.)

---

## One-line summary for forcing-catalogue update

> **Door 4e (2026-06-02): PARTIAL-DOOR (CONFIRMED).** Literal Falsifier
> A' returned NEGATIVE on three candidate AC tensor-product
> compatibility conditions (SU(2)-Ad-invariance, Hopf U(1)-equivariance,
> principal-bundle compatibility) — none distinguish ℍ from M₂(ℂ) at
> n=2. Standard CCM AC is genuinely too permissive. **Substantive
> content:** CCM ℍ-selection requires 3 axioms (second-order condition
> + complex chiral fermion rep + dimension count); DAS ℍ-selection
> requires 1 axiom (sphere = unit-norm subset of division algebra). Both
> are imports relative to bare Connes axioms; DAS is leaner. The
> FULL-DOOR upgrade path is to adopt **Upgrade B (sphere-Lie-group
> axiom)** as foundational: "inner algebra at rung n = *-algebra whose
> unit group IS rung-n sphere as a Lie group, when one exists; else
> fallback M_n(ℂ)." Paper-level structural-axiom upgrade,
> sprint-reachable. Falsifier A'' (does Upgrade B introduce new
> downstream constraints?) is the next thread.
