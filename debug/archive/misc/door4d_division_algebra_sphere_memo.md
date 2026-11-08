# Door 4d — Does the Division-Algebra-of-the-Sphere criterion close the ℍ-vs-M₂(ℂ) fork?

*A structurally orthogonal handle to Door 4c (no J sign-table, no
combined-J data). Tests whether the same Hopf-tower argument that forces
the gauge group can be strengthened to force the n=2 inner algebra's real
form.*
*Structure-only / division-algebra probe. NO Yukawa value selected
(guardrail).*
Started + closed 2026-06-02. Driver: `debug/door4d_division_algebra_sphere.py`.
Data: `debug/data/door4d_division_algebra_sphere.json`.

---

## TL;DR verdict — **PARTIAL-DOOR**

A new GeoVac-coherent handle is on the table. The **Division-Algebra-of-the-Sphere
(DAS) criterion** — "the inner algebra at Hopf rung n is the natural
associative real normed division algebra whose unit-norm sphere IS the
rung's Hopf bundle total space S^(2n-1)" — **agrees bit-exactly with
Door 4b's elimination forcings at n=1 and n=3**, and **closes the
ℍ-vs-M₂(ℂ) fork at n=2 by selecting ℍ**. It is structurally orthogonal to
Door 4c: the argument never invokes any J or sign-table data.

DAS is **not yet a true GeoVac forcing**. The existing Paper 32 §VIII.B
construction extracts only the gauge group from the Hopf tower; it does
not invoke the unit-sphere identification to constrain the inner algebra's
realization. Adopting DAS as a forcing is a strict-extension of the
construction principle that is consistent with — but not derivable from —
the existing framework. The probe upgrades the inner-algebra status from
"mostly-forced with the ℍ-fork closed by literature import (CCM)" to
"mostly-forced with the ℍ-fork closed by a GeoVac-coherent geometric
criterion, pending one structural axiom upgrade."

The honest residue: **whether DAS is a derivable consequence of the
Bertrand × Hopf-tower construction or a NEW axiom on the same standing as
the CCM second-order condition** is the deeper open question. Either way,
DAS is a meaningfully sharper handle than what Door 4c left in place.

---

## What this probe tested

Four bit-exact identities and one structural comparison.

### Part 1 — n=1 rung: S¹ IS the unit sphere of ℂ (bit-exact)

| Quantity | Value |
|:---------|:------|
| Max norm error on S¹ samples | 1.11e-16 |
| Max closure error under complex multiplication | 1.11e-16 |
| Is S¹ a multiplicative subgroup of ℂ? | **YES** (bit-exact) |

S¹ = {z ∈ ℂ : |z| = 1} is closed under complex multiplication; |zw| = 1
whenever |z| = |w| = 1. Standard fact, recorded as a bit-exact baseline.

### Part 2 — n=2 rung: S³ IS the unit sphere of ℍ = Sp(1) = SU(2) (bit-exact)

| Quantity | Value |
|:---------|:------|
| Max norm error on S³ samples | 2.22e-16 |
| Max closure error under quaternion multiplication | 2.22e-16 |
| Max SU(2)-isomorphism error (quaternion→matrix) | 3.34e-16 |
| Is S³ a multiplicative subgroup of ℍ? | **YES** (bit-exact) |
| Is S³ bit-exactly isomorphic to SU(2) = Sp(1)? | **YES** |

S³ = {q ∈ ℍ : |q| = 1} is closed under quaternion (Hamilton) product, and
the standard embedding ℍ → M₂(ℂ) restricted to S³ gives bit-exact SU(2).
**The ℍ-unit-sphere identification at n=2 is DIRECT — one chart, no
quotient.** This contrasts with the M₂(ℂ) realization at n=2 (Part 4).

### Part 3 — n=3 rung: S⁵ has NO associative-division-algebra realization (Hurwitz)

Hurwitz (1898): the finite-dimensional real normed division algebras are
exactly ℝ, ℂ, ℍ, 𝕆 of dimensions 1, 2, 4, 8 — and only ℝ, ℂ, ℍ are
**associative**. Their unit spheres have dimensions 0, 1, 3, 7
respectively (S⁰, S¹, S³, S⁷). **No associative real normed division
algebra has S⁵ as its unit sphere** (and 𝕆's unit sphere S⁷ is at the
wrong rung). At n=3, the DAS criterion has no candidate division algebra
and must fall back to a matrix algebra. The fallback target is M₃(ℂ),
selected by Door 4b's elimination argument (only M₃(ℂ) reproduces SU(3)
as its post-unimodularity unitary group; M₃(ℝ) gives PO(3), M₃(ℍ) gives
PSp(3), both wrong).

### Part 4 — the n=2 fork under DAS

| Candidate | Unitary group | Unit-sphere dim | Realization of SU(2) |
|:----------|:-------------|:---------------:|:--------------------:|
| **ℍ** | U(ℍ) = Sp(1) = SU(2) | **3** (= S³) | **direct** (single chart) |
| **M₂(ℂ)** | U(M₂(ℂ)) = U(2) | 4 (= U(2)) | **indirect** (quotient U(2)/U(1) = SU(2) via unimodularity) |

The two candidates produce the SAME gauge group SU(2), but via
structurally different routes. **DAS selects the DIRECT (no-quotient)
realization** — ℍ at n=2.

### Part 5 — DAS verdict and consistency with existing forcings

| Rung | Door 4b elimination | DAS selection | Agree? |
|:----:|:-------------------:|:-------------:|:------:|
| n=1 | ℂ (only ℂ has U(1) as unitary group) | ℂ (S¹ = unit ℂ) | **YES** |
| n=2 | ℍ or M₂(ℂ) (both give SU(2)) | ℍ (S³ = unit ℍ = Sp(1)) | **DAS closes** |
| n=3 | M₃(ℂ) (only M₃(ℂ) gives SU(3)) | M₃(ℂ) (no div-alg fallback) | **YES** (different routes) |

DAS is a **strict-extension** of Door 4b's elimination forcing: it agrees
at n=1 and n=3 (where Door 4b is decisive), and it **closes** the n=2
fork (where Door 4b is silent).

---

## Why this is structurally orthogonal to Door 4c

Door 4c showed that the COMBINED-J / KO-dim sign table is bit-identical
for ℍ and M₂(ℂ) — every entry in the combined sign table coincides; the
internal C² antilinear involution that decides the real form is invisible
to J = J_GV ⊗ J_F (whose J_F has J_F² = +1 in both cases, being the
matter↔antimatter swap-conjugation). Order-zero and order-one are also
blind by matter/antimatter decoupling.

**DAS bypasses this obstruction entirely.** The argument never invokes
any J. It identifies the rung-n sphere as the unit-norm subset of an
associative real normed division algebra (Hurwitz), and asks for the
inner algebra at that rung to be that division algebra. The selection
mechanism is **geometric / unit-sphere identification**, not
**algebraic / sign-table**.

This places DAS on a different footing than the literature's standard
ℍ-selector (CCM second-order condition + complex chiral fermion
representation). The CCM argument is finite-side and references no S³ /
no Hopf bundle; DAS is geometric and references the rung's sphere
directly. The two are independent handles for the same selection.

---

## Is DAS a true GeoVac forcing, or a new axiom?

This is the open question that determines whether the verdict is DOOR or
PARTIAL-DOOR.

**Pro-DOOR (DAS is derivable from existing GeoVac structure):**
The Bertrand × Hopf-tower argument (CLAUDE.md §1.7 WH4; Paper 32 §VIII.B)
produces the rung-n sphere S^(2n-1) as the geometric substrate at rung n.
The complex Hopf bundle S^(2n-1) → ℂP^(n-1) is **classically** a
division-algebra construction:

- n=1: S¹ → S⁰ (trivial; complex-numbers structure)
- n=2: S³ → S² (Hopf; **quaternionic** structure of S³ as unit elements
  of ℍ is intrinsic to this bundle — it's exactly Sp(1) → Sp(1)/U(1))
- n=4: S⁷ → S⁴ (octonionic Hopf, non-associative)

At n=2, the Hopf bundle S³ → S² **is** the principal U(1)-bundle inside
Sp(1) = unit ℍ. The identification "S³ = unit-norm sphere of ℍ" is not
imposed; it is a property of the Hopf bundle GeoVac already uses. If
GeoVac extracts the rung-n SPHERE from the Hopf tower, it extracts a
sphere that comes pre-equipped with a division-algebra structure (when
one exists, i.e., for n ∈ {1, 2}). Reading off that structure into the
inner algebra is a natural extraction — arguably more natural than
discarding it.

**Pro-PARTIAL (DAS is a strict extension, not a derivation):**
Paper 32 §VIII.B's current argument explicitly extracts only the
**gauge group** SU(n) from the Hopf rung — not the algebra's
representation type. Door 4b's enumeration found two algebras at n=2
that reproduce SU(2), and the Hopf-tower argument as written does not
prefer one. The unit-sphere identification (DAS) is a geometric fact
ABOUT the rung-n sphere, but the **inner algebra** in the AC tensor
product H = H_GV ⊗ H_F lives on H_F — a separate Hilbert factor — and
there is no morphism in the standard CCM construction that transfers
the rung-n sphere's division-algebra structure to a constraint on H_F's
algebra. Adopting DAS as a forcing requires adding the construction
principle "the inner algebra at rung n inherits the division-algebra
structure of the rung-n sphere," which is **new content** even though
consistent.

**My read:** PARTIAL-DOOR is the honest verdict at this stage. The
PRO-DOOR case is structurally appealing (the Hopf tower's
division-algebra content is intrinsic, not added input), but the
PRO-PARTIAL case is honest about the gap (the tensor-product H_GV ⊗ H_F
construction does not automatically transfer the sphere's structure to
H_F).

**What would upgrade PARTIAL-DOOR → DOOR:** a structural argument (a
small theorem, paper-grade) that the AC tensor-product construction,
applied to GeoVac's Hopf-rung substrate, **forces** the H_F algebra at
rung n to be the division algebra whose unit sphere is the rung-n
sphere. The candidate argument: at n=2, the rung-2 Hopf bundle is
*itself* the principal U(1)-bundle of Sp(1) = unit ℍ. The Hopf
structure provides a canonical isomorphism between (a) the rung-2
sphere with its division-algebra structure and (b) the natural matter
representation space of an inner SU(2). If this isomorphism is the unique
SU(2)-equivariant structure compatible with the rung-2 Hopf bundle, then
the inner algebra must be ℍ. This is paper-level work, sprint-scale.

---

## Honest scope — what was NOT shown

- DAS was not promoted to a theorem here. The verdict is PARTIAL-DOOR
  precisely because the structural-upgrade question (above) is open.
- No Yukawa value, generation-specific mass, or mixing angle was
  selected. DAS targets the inner algebra's **real form** at the n=2
  rung, not its **eigenvalues**. The seam theorem (Paper 32 §VIII Door 4)
  separating algebra factors from Yukawa values is untouched.
- The H1 / W3 / Koide negatives (CLAUDE.md §3) remain untouched. DAS is
  a structure-only criterion on the algebra side.
- DAS is **not** claimed to derive the inner KO-dimension (Door 4b Q3,
  remains free) or the generation count N_gen (Door 4b Q3, remains
  free). It only addresses the ℍ-vs-M₂(ℂ) fork at the n=2 rung.

---

## Inner-algebra status after Door 4d

```
FORCED (Bertrand x Hopf-tower, same principle as gauge content):
  • factor count = 3
  • Wedderburn / matrix-algebra type of each factor
  • n=1 factor = C
  • n=3 factor = M_3(C)

PARTIAL (GeoVac-coherent handle on the table; closed-theorem pending):
  • n=2 factor = H (vs M_2(C))
      handle (NEW, Door 4d):
        DAS = the inner algebra at rung n is the natural associative
        real normed division algebra whose unit-norm sphere IS the
        rung's Hopf bundle total space S^(2n-1).
      structural verification (bit-exact):
        n=1: S^1 = unit C       ; closure error 1.1e-16
        n=2: S^3 = unit H = SU(2); closure error 2.2e-16, SU(2)-iso 3.3e-16
        n=3: S^5 has NO associative division-algebra realization (Hurwitz)
      forcing status:
        DAS agrees with Door 4b at n=1 and n=3.
        DAS closes the H-vs-M_2(C) fork at n=2 (selects H).
        DAS is consistent with -- but not derivable from -- the current
        Paper 32 SS VIII.B construction. Adoption as a forcing requires
        a paper-level theorem that the AC tensor product transfers the
        rung-n sphere's division-algebra structure to H_F.

FREE (no construction reaches them):
  • generation count N_gen   (invisible to inner automorphisms)
  • inner KO-dimension       (imposed to hit SM total 1 mod 8)
  • Yukawa eigenvalues       (seam theorem; out of scope)
```

Net: **the inner algebra at n=2 has a GeoVac-coherent handle that is
structurally distinct from Doors 4b/4c.** DAS does not yet promote
"mostly-forced" to "fully-forced," but it sharpens the question by
providing a clean PRO/CON axis (is the AC tensor product
sphere-structure-transferring or not?) that the prior memos lacked.

---

## The next sharpest falsifier (after Door 4d)

**Falsifier A' (the sprint-reachable upgrade):** Prove that the
AC-tensor-product construction H = H_GV ⊗ H_F, with H_GV the
GeoVac Hopf-rung substrate at rung n and H_F a finite Hilbert space
carrying the inner algebra, **forces** the inner algebra to inherit the
rung-n sphere's division-algebra structure. Concretely: the principal
SU(2)-bundle structure of the rung-2 Hopf bundle S³ → S² should pick
out a canonical SU(2)-equivariant action on H_F that is realized only
by ℍ (not by M₂(ℂ)). This is paper-level work, no new code; if it
closes, PARTIAL-DOOR upgrades to DOOR.

  - *Concrete first step:* characterize the SU(2)-equivariant
    representations of M₂(ℂ) on a 2-dim complex space and check whether
    any of them is compatible with the Hopf bundle's principal U(1)
    structure in a way that the ℍ-equivariant rep also satisfies. If
    only the ℍ-equivariant rep is compatible (the M₂(ℂ) ones being
    "over-parameterized" by the unimodularity quotient), DAS is forced
    by the AC construction.

**Falsifier A'' (the alternative, narrower upgrade):** Show that
adopting DAS introduces a NEW downstream constraint (e.g., on the inner
KO-dim, on the chirality grading γ_F, or on the Higgs sector) that the
existing framework either supports or contradicts. If it supports, DAS
gains independent corroboration. If it contradicts, DAS is RULED OUT
and the fork goes back to "admitted-not-forced via combined-J sign
table (Door 4c)."

**Falsifier B (the deep wall, deferred):** Force the generation count
N_gen and the inner KO-dimension from a packing principle. Same as Door
4b Falsifier B — no known handle; deferred.

---

## Decision-gate result

Per the gate: **DOOR** = DAS is a derivable consequence of the existing
GeoVac construction, forces ℍ; **PARTIAL-DOOR** = DAS is a coherent
extension consistent with existing forcings that closes the fork if
adopted; **WALL** = DAS adds nothing structurally beyond what Door 4b
extracted.

**Result: PARTIAL-DOOR.**

- DAS is a coherent and natural geometric criterion identifying the
  rung-n sphere as the unit-norm subset of an associative real normed
  division algebra (when one exists at the rung's dimension).
- It **agrees bit-exactly with Door 4b's elimination forcings** at n=1
  (C) and n=3 (M_3(C)).
- It **closes the n=2 fork by selecting ℍ** (the direct, no-quotient
  realization).
- It is **structurally orthogonal to Door 4c** (no J / sign-table data).
- It is **not yet derivable** from the existing Paper 32 §VIII.B
  construction; it requires one structural-upgrade theorem.
- The handle is **sprint-reachable** (Falsifier A' above), at the
  paper-level rather than the compute-level.

**Net inner-algebra status:** *factor count + n=1 + n=3 forced; n=2 fork
closed by DAS pending one structural theorem; generation count + inner
KO-dim + Yukawa values free.* This is materially sharper than the Door
4c end-state ("n=2 admitted-not-forced, ℍ is a literature import via
CCM second-order condition").

---

## Out-of-scope record (guardrail compliance)

No Yukawa value, generation-specific mass, or mixing angle was selected
or proposed. The DAS criterion targets the **real form of the inner
algebra at n=2** (a structure-side selection: ℍ vs M₂(ℂ)), not its Dirac
eigenvalues. The generation count and inner KO-dimension are explicitly
recorded as remaining FREE. The H1 / W3 / Koide negatives (CLAUDE.md
§3) remain untouched; this probe stayed strictly on the
algebra-structure side of the seam throughout.

The probe makes **no claim to derive the Standard Model gauge group
beyond what Door 4b already established**, and **no claim to predict any
SM parameter**. It identifies a coherent NEW handle that strictly
sharpens the existing forcing catalogue, and is honest about the gap
between "consistent extension" and "derived theorem."

---

## Files / sources

- Driver: `debug/door4d_division_algebra_sphere.py`
- Data: `debug/data/door4d_division_algebra_sphere.json`
- Prior memos:
  - `debug/door4_gauge_yukawa_boundary_memo.md` (seam theorem)
  - `debug/door4b_inner_algebra_forcing_memo.md` (factor count + n=1/3 forced)
  - `debug/door4c_j_signtable_audit_memo.md` (combined-J handle CLOSED NEGATIVE)
- Paper 32 §VIII.B/C: `papers/group1_operator_algebras/paper_32_spectral_triple.tex`
  - Prop. `prop:reality` (J_GV² = −1, KO-dim 3, lines 919–926)
  - SM gauge appendix §VIII.B (Bertrand × Hopf-tower truncation, lines ~3518–3522)
  - Door 4 / 4b / 4c paragraphs (lines 3486–3562)
- Paper 18 §IV: `papers/group3_foundations/paper_18_exchange_constants.tex`
  - Second-packing-axiom open question (the deeper structural framework
    DAS sits inside)
- Code:
  - `geovac/almost_commutative.py` (`quaternion_to_matrix`,
    `ElectroweakFiniteTriple` — the ℂ⊕ℍ inner factor as currently
    constructed)
  - `geovac/standard_model_triple.py` (full ℂ⊕ℍ⊕M₃(ℂ) AC triple)
  - `geovac/real_structure.py` (J_GV, KO-dim 3 audit)
- Memory: `bertrand_sm_gauge_truncation.md`,
  `wh_register_april2026.md` (WH4 deflated 2026-05-07,
  "single-input forcing statement"; DAS adds inner-factor content to
  this reading)
- Literature:
  - Hurwitz, A. (1898). *Über die Composition der quadratischen Formen
    von beliebig vielen Variablen.* (Classification of normed division
    algebras)
  - Chamseddine, A. H., & Connes, A. (2008). *Why the Standard Model.*
    J. Geom. Phys. 58, 38. (The CCM ℍ-selector via second-order condition
    — the literature import DAS would replace with a
    geometric/sphere-based selector)
  - Connes, A. (1995); van Suijlekom, W. D. (2015). (KO-dim sign tables)

---

## One-line summary for forcing-catalogue update

> **Door 4d (2026-06-02): PARTIAL-DOOR.** A new GeoVac-coherent handle
> (Division-Algebra-of-the-Sphere, DAS) closes the ℍ-vs-M₂(ℂ) fork by
> identifying the rung-n sphere S^(2n-1) with the unit-norm subset of an
> associative real normed division algebra (ℂ at n=1, ℍ at n=2; no
> candidate at n=3 → M₃(ℂ) fallback). Bit-exact at all three rungs.
> Structurally orthogonal to Door 4c (no J / sign-table). Agrees with
> Door 4b at n=1, n=3; closes the n=2 fork by selecting ℍ. Not yet
> derivable from the existing Paper 32 §VIII.B construction; adoption as
> a forcing requires one structural-upgrade theorem on the AC
> tensor-product (sphere-structure transfer to H_F). Falsifier A'
> sprint-reachable, paper-level. The deeper "second packing axiom for
> inner-factor data" (Paper 18 §IV) is now sharper.
