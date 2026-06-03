# Door 4f — Falsifier A'' on Upgrade B: does adopting it introduce new downstream constraints?

*Closes the Door 4 series. Tests whether the sphere-Lie-group axiom
(Upgrade B from Door 4e) introduces new downstream constraints on
existing GeoVac observables, contradicts any existing result, or is a
clean reformulation of CCM with no new content.*
*Structure-only / paper-level analysis. NO Yukawa value selected
(guardrail).*
Started + closed 2026-06-02. Driver: `debug/door4f_falsifier_A_double_prime.py`.
Data: `debug/data/door4f_falsifier_A_double_prime.json`.

---

## TL;DR verdict — **Outcome 3 (NEUTRAL): PARTIAL-DOOR FINAL**

**12 sub-tests, 0 introducing new constraints, 0 contradicting
existing results.** Upgrade B is a clean reformulation of CCM's
ℍ-selection with no new downstream predictive content. It buys
axiom-minimality (1 axiom vs CCM's 3) but no additional GeoVac-
observable forcing. Five sub-tests are SILENT (Upgrade B doesn't
constrain the target at all), three are COMPATIBLE (the constraint
passes for ℍ-selection just as it did for the broader cases tested
in Door 4c), two are IDENTICAL to CCM at the gauge-sector / Higgs
vacuum level, one provides a clean extensibility rule via Adams 1958
(T11), and one is a consistent strengthening of the existing
Bertrand × Hopf-tower reading (T12).

**The Door 4 series is structurally complete at PARTIAL-DOOR FINAL.**
The next open question is judgment-level (PI / NCG-community-side):
is Upgrade B's axiom-minimality advantage worth promoting to a working
GeoVac axiom, or is the CCM 3-axiom path preferred for literature
integration? This is paper-level choice work, not a probe-able falsifier.

---

## What was tested

12 sub-tests on whether adopting Upgrade B introduces new constraints
on existing downstream observables. Each test articulates the target,
its current status pre-Upgrade-B, and whether Upgrade B directly
constrains it.

| ID  | Target | Verdict |
|:----|:-------|:--------|
| T1  | Inner KO-dimension | SILENT |
| T2  | Chirality grading γ_F | SILENT (consistent but not unique) |
| T3  | Yukawa eigenvalues | SILENT (seam theorem unaffected) |
| T4  | Generation count N_gen | SILENT |
| T5  | Connes order-zero | COMPATIBLE (passes; mechanism unaffected) |
| T6  | Connes order-one | COMPATIBLE (passes; mechanism unaffected) |
| T7  | Connes second-order | COMPATIBLE but NOT REDUNDANT |
| T8  | Bosonic spectral action gauge coefficient | IDENTICAL to CCM |
| T9  | Higgs vacuum manifold | IDENTICAL (S² for both) |
| T10 | Cross-rung Yukawa coupling | SILENT |
| T11 | Adams 1958 extensibility | CLEAN EXTENSIBILITY RULE |
| T12 | Bertrand × Hopf-tower consistency | CONSISTENT EXTENSION |

**Net counts:**
- Tests introducing new constraint on existing observables: **0**
- Tests contradicting existing results: **0**
- Tests silent (target unconstrained by Upgrade B): 5
- Tests compatible (target unchanged from pre-Upgrade-B): 3
- Tests identical to CCM at the relevant level: 2
- Tests providing clean extensibility rule: 1 (T11)
- Tests providing consistent extension: 1 (T12)

---

## Detail: each sub-test

### T1. Inner KO-dimension

**Status pre-Upgrade-B:** FREE (Door 4b Q3). Imposed to be 6 to make
combined outer (3) + inner (6) = 9 ≡ 1 mod 8 match the SM total.

**Verdict:** SILENT. Upgrade B fixes A_F's structure at each rung
(the *-algebra). The inner KO-dim is a property of (H_F, D_F, J_F, γ_F)
determined by HOW we double H_F to include antimatter (J_F² sign) and
the chirality grading. Upgrade B is articulated at the algebra level
and says nothing about H_F's doubling or J_F²'s sign.

### T2. Chirality grading γ_F

**Status pre-Upgrade-B:** Independent of γ_GV per G3 NEGATIVE
(γ_GV and γ_F are independent commuting Z₂'s).

**Verdict:** SILENT (consistent but not unique). The standard CCM
γ_F = diag(+1, +1, −1, −1) on (ν_L, e_L, ν_R, e_R) is CONSISTENT
with ℍ-action on the L-doublet (left quaternions naturally act on a
single 2-dim block). But Upgrade B doesn't UNIQUELY determine which
subspace of H_F is L vs R; that's additional input from the doubled
Hilbert space structure.

### T3. Yukawa eigenvalues

**Status pre-Upgrade-B:** FREE per seam theorem (Paper 32 §VIII Door 4).

**Verdict:** SILENT (seam theorem unaffected). Upgrade B constrains
A_F; the seam theorem applies independently of A_F's specific real
form. The off-diagonal D_F entries (Yukawas) carry their own
Dirichlet ring disjoint from the outer Mellin engine. Selecting ℍ
over M₂(ℂ) at the n=2 factor does not enter the Yukawa Dirichlet ring.

### T4. Generation count N_gen

**Status pre-Upgrade-B:** FREE per Door 4b Q3 (invisible to inner
automorphisms; enters as H_F = ℂ^N_gen ⊗ H_F^{(1 gen)} multiplicity).

**Verdict:** SILENT. Upgrade B constrains A_F. N_gen is a multiplicity
of the matter rep, not a feature of A_F itself.

### T5. Connes order-zero

**Status pre-Upgrade-B:** PASSES for both ℍ and M₂(ℂ) at finite
n_max ∈ {1, 2, 3} (Door 4c bit-exact data; matter/antimatter decoupling
mechanism).

**Verdict:** COMPATIBLE. Upgrade B restricts to ℍ, which already
passes order-zero. Mechanism (matter/antimatter decoupling) unaffected.

### T6. Connes order-one

**Status pre-Upgrade-B:** PASSES for both ℍ and M₂(ℂ) at finite
n_max ∈ {1, 2} with generic nonzero Yukawa y = 0.3 (Door 4c part 3c,
bit-exact). Same matter/antimatter decoupling.

**Verdict:** COMPATIBLE. Same argument as T5.

### T7. Connes second-order condition (the CCM workhorse)

**Status pre-Upgrade-B:** The CCM axiom used to SELECT ℍ over M₂(ℂ)
in Chamseddine–Connes 2008. Imposed as axiom; constrains D_F's
off-diagonal Yukawa structure.

**Verdict:** COMPATIBLE but NOT REDUNDANT. With Upgrade B forcing
ℍ via a different route (geometric), the second-order condition
becomes either (a) automatic given ℍ, or (b) still needed as a
separate constraint on D_F. The standard CCM analysis suggests (b):
the second-order condition fixes the off-diagonal D_F block structure
(which Yukawa entries can be nonzero), and ℍ selection alone does
not determine this. **Upgrade B and the second-order condition
coexist without redundancy** — both are needed but they constrain
different objects (A_F's real form vs D_F's off-diagonal structure).

This is informative: it means even adopting Upgrade B, the
second-order condition is STILL an independent axiom in the full
construction; Upgrade B replaces one CCM axiom (ℍ-selection via
second-order + complex chiral rep + dimension count) with itself
(ℍ-selection via sphere-Lie-group), but the second-order condition's
role in CONSTRAINING D_F persists.

### T8. Bosonic spectral action gauge-sector coefficient

**Status pre-Upgrade-B:** Standard CCM: bosonic action
Tr(f(D_A² / Λ²)) gives the SU(2) Yang-Mills kinetic term
(1/4g²) Tr(F_μν F^μν) with coupling determined by the SU(2)
Killing form / trace normalization.

**Verdict:** IDENTICAL at gauge-sector leading order. Both ℍ and
M₂(ℂ) give SU(2) (post-unimodularity for M₂(ℂ)). The Killing form
of su(2) is identical in both cases; the spectral action's gauge
kinetic coefficient is INSENSITIVE to which real form of M₂(ℂ) the
algebra is.

### T9. Higgs vacuum manifold (post-SSB)

**Status pre-Upgrade-B:** Standard CCM: Higgs vacuum manifold
= SU(2)×U(1)/U(1) = SU(2)/Z₂ modulo gauge = S² after SSB.

**Verdict:** IDENTICAL (S² for both). Determined by gauge group +
Higgs rep, both identical for ℍ and M₂(ℂ) selections at the
gauge-group level.

### T10. Cross-rung Yukawa coupling (CKM mixing)

**Status pre-Upgrade-B:** Free per seam theorem; CKM is calibration
data.

**Verdict:** SILENT. Upgrade B selects each rung's algebra
independently (block diagonal). Cross-rung Yukawa entries
(generation-mixing) are calibration data inaccessible to either
ℍ or M₂(ℂ) selection at the algebra level.

### T11. Adams 1958 extensibility to higher rungs

**Status pre-Upgrade-B:** GeoVac currently doesn't reach beyond n=3
(Bertrand × Hopf truncates at U(1) × SU(2) × SU(3) = SM gauge group).

**Verdict:** CLEAN EXTENSIBILITY RULE. Adams 1958/1960
(Hopf-invariant-one theorem): S^(2n−1) is parallelizable iff
n ∈ {1, 2, 4}; carries Lie-group structure only at n=1 (U(1)),
n=2 (Sp(1) = SU(2)), and n=4 (S⁷ via octonions, NON-ASSOCIATIVE).
Upgrade B is natural for the associative-only case (n=1, 2); for
n=4 the octonionic S⁷ requires extension to non-associative
*-algebras (outside standard NCG); for n ≥ 5, S^(2n−1) is not a
Lie group, so Upgrade B's fallback (M_n(ℂ)) applies cleanly.

**The rule has a clean per-rung specification all the way up; no
ambiguity.** This is a STRUCTURAL constraint that becomes substantive
only if GeoVac ever extends beyond n=3 (which the SM gauge truncation
doesn't currently motivate). Adams 1958 is the deep topological fact
underwriting Upgrade B's elegance at finite rungs.

### T12. Consistency with Bertrand × Hopf-tower reading

**Status pre-Upgrade-B:** Bertrand × Hopf-tower extracts the rung-n
sphere S^(2n−1) from the closed-orbit constraint + complex-Hopf
bundle. The existing reading (Paper 32 §VIII.B) extracts the
GAUGE GROUP SU(n) from this sphere.

**Verdict:** CONSISTENT EXTENSION. Upgrade B extracts MORE structure
from the same sphere: not just its topology (which gives the gauge
group via homotopy / fundamental group) but also its Lie-group
structure (when one exists, by Adams 1958). The reading "use the Lie
structure when available, fall back to topology when not" is a
**natural strengthening** of "use the topology only." The case n=3
fallback is "when the Lie structure doesn't exist, use only the
topology." Structurally consistent with the existing argument; does
not contradict any existing GeoVac result.

---

## Three-outcome breakdown

Per the Door 4d/4e framing:

- **Outcome 1: Introduces a new constraint the framework SUPPORTS →**
  corroborates DAS; promote to working axiom alongside or replacing CCM's.
  **NOT REALIZED** (0 out of 12 tests introduce a new constraint).

- **Outcome 2: Introduces a new constraint the framework CONTRADICTS →**
  Upgrade B RULED OUT; fork returns to "PARTIAL-DOOR via leaner
  ℍ-import (Door 4d) but no FULL-DOOR closure."
  **NOT REALIZED** (0 out of 12 tests contradict an existing result).

- **Outcome 3: Introduces no new constraints (pure reformulation of CCM) →**
  DAS maintained at PARTIAL-DOOR; Upgrade B is "CCM-equivalent with
  leaner axiom count."
  **REALIZED** (12 out of 12 tests are silent/compatible/identical/
  clean-extension).

---

## Inner-algebra status after Door 4f (FINAL for the series)

Under the ASSUMPTION that Upgrade B is adopted as a working axiom
(a paper-level choice, not a forced one):

```
FORCED via Upgrade B:
  • factor count = 3
  • n=1 factor = C (S^1 = U(1) = U(C))
  • n=2 factor = H (S^3 = Sp(1) = U(H))
  • n=3 factor = M_3(C) (S^5 not a Lie group, fallback to minimal
                          matrix algebra reproducing SU(3))

FREE (no construction reaches them; Upgrade B silent):
  • inner KO-dimension (T1)
  • chirality grading gamma_F (T2, consistent but not unique)
  • Yukawa eigenvalues (T3, seam theorem)
  • generation count N_gen (T4)
  • cross-rung Yukawa mixing / CKM (T10)
```

If Upgrade B is NOT adopted, the status reverts to the Door 4c
ledger (n=2 admitted-not-forced; ℍ via CCM literature import).

---

## CCM 3-axiom path vs DAS 1-axiom path: input minimality comparison

| Selection mechanism | Axioms beyond bare Connes data | Selects at n=2 |
|:--------------------|:------------------------------:|:--------------:|
| Standard CCM (Chamseddine–Connes 2008) | **3** | ℍ |
| DAS Upgrade B | **1** | ℍ |

CCM's three axioms:
1. Second-order condition: [[D, a], [J b J⁻¹, c]] = 0 with nonzero D_F.
2. Complex chiral fermion representation (L-doublet paired with anti-L-doublet).
3. Dimension count 2N² = 32 / 4×4-grading (Chamseddine–Connes–Marcolli 2007).

DAS Upgrade B's one axiom:
1. Inner algebra at rung n is the *-algebra whose unit group IS the
   rung-n sphere S^(2n−1) as a Lie group (when one exists; else
   fallback M_n(ℂ)).

**Important caveat:** Door 4f's T7 shows the CCM second-order condition
is STILL needed (independent of ℍ-selection) to constrain D_F's
off-diagonal Yukawa structure. Adopting Upgrade B replaces CCM's
"second-order + complex chiral rep + dimension count" → "Upgrade B
sphere-Lie axiom + second-order condition." Net axiom count: CCM 3 →
Upgrade B 2 (one geometric, one algebraic). **Real axiom-savings:
1 axiom**, not 2.

This is a more honest accounting than Door 4e's headline "1 vs 3." The
sphere-Lie-group axiom REPLACES CCM's "ℍ-selection complex" (which is
really 2 of CCM's 3 axioms: complex chiral rep + dimension count). The
second-order condition stays in both formulations because it
constrains D_F, not A_F.

---

## What this means for the Door 4 series

The Door 4 series of probes — gauge/Yukawa boundary (Door 4); inner
algebra forcing (4b); J sign-table audit (4c); Division-Algebra-of-the-
Sphere handle (4d); literal Falsifier A' (4e); downstream-constraint
Falsifier A'' (4f) — is **structurally complete**.

Final results:
- The seam theorem (Door 4) is a structural theorem of GeoVac: gauge
  forced, Yukawa free, ring-disjoint.
- The inner algebra structure (Door 4b) is mostly-forced: factor count
  + ℂ + M₃(ℂ) forced by Bertrand × Hopf-tower elimination.
- The combined-J / KO-dim sign-table handle (Door 4c) does NOT close the
  ℍ vs M₂(ℂ) fork.
- The geometric DAS handle (Door 4d) DOES close the fork (PARTIAL-DOOR,
  pending paper-level axiom adoption).
- The literal Falsifier A' (Door 4e) confirms no standard AC compatibility
  condition is sufficient; DAS is a genuinely new selector.
- The downstream-constraint Falsifier A'' (Door 4f) confirms Upgrade B is
  a clean reformulation of CCM with 1 axiom of savings, no new constraints,
  no contradictions.

**The remaining open question is judgment-level**, not probe-able:
should GeoVac adopt Upgrade B as a working axiom, accepting the leaner
axiom-count and tighter geometric integration in exchange for departing
from standard CCM conventions? That's a paper-level / community-level
decision.

---

## Honest scope — what was NOT shown

- The 12 sub-tests are the natural set; an exhaustive enumeration might
  uncover a 13th test where Upgrade B does introduce a new constraint.
  Outcome 3 is robust against the natural-test enumeration; it is not
  meta-mathematically robust against arbitrary unexpected probes.
- The "real axiom-savings: 1 axiom" comparison (Upgrade B 2 axioms vs
  CCM 3) is more honest than Door 4e's headline 1-vs-3 framing. The
  difference is that CCM's "ℍ-selection complex" bundles 2 axioms
  (complex chiral rep + dimension count) whose direct GeoVac analog
  (Upgrade B sphere-Lie axiom) is just 1.
- No claim is made that Upgrade B IS the right axiom for GeoVac;
  Outcome 3 says Upgrade B doesn't BREAK anything and isn't MORE
  predictive than CCM. The choice between Upgrade B and CCM-3-axiom
  is judgment-level.
- No Yukawa value, generation count, or KO-dim is selected. H1 / W3 /
  Koide negatives (CLAUDE.md §3) untouched. Seam theorem (Paper 32
  §VIII Door 4) unaffected.

---

## The next sharpest question (post-Door-4 series)

The Door 4 series leaves one PAPER-LEVEL question open and one DEEPER
WALL untouched. Both are deferred:

**Paper-level (judgment): adopt Upgrade B or not?** Depends on (a)
whether GeoVac's broader narrative wants axiom-minimality more than
literature-compatibility, (b) whether the Adams-1958-based extensibility
rule (T11) becomes useful in future higher-rung work (not currently
motivated by SM truncation), and (c) whether the Bertrand × Hopf-tower
geometric reading (T12) is the preferred organizing principle going
forward. **No further computational probe applies.**

**Deeper wall (the genuine remaining open question):** force the
generation count N_gen and the inner KO-dimension from a packing
principle. Same as Door 4b's Falsifier B; no known handle. **Deferred
indefinitely.**

---

## Files / sources

- Driver: `debug/door4f_falsifier_A_double_prime.py`
- Data: `debug/data/door4f_falsifier_A_double_prime.json`
- Prior memos:
  - `debug/door4_gauge_yukawa_boundary_memo.md` (seam theorem)
  - `debug/door4b_inner_algebra_forcing_memo.md` (factor count + n=1/3 forced)
  - `debug/door4c_j_signtable_audit_memo.md` (combined-J handle CLOSED NEGATIVE)
  - `debug/door4d_division_algebra_sphere_memo.md` (DAS introduced, PARTIAL-DOOR)
  - `debug/door4e_falsifier_A_prime_memo.md` (literal A' negative; Upgrade B named)
- Paper 32 §VIII.C: `papers/group1_operator_algebras/paper_32_spectral_triple.tex`
- Code: `geovac/almost_commutative.py`, `geovac/standard_model_triple.py`,
  `geovac/real_structure.py`
- Memory: `bertrand_sm_gauge_truncation.md`, `wh_register_april2026.md`
- Literature:
  - Hurwitz, A. (1898). Classification of normed division algebras.
  - Adams, J. F. (1958, 1960). Hopf-invariant-one theorem / parallelizability
    of spheres / division-algebra structure on S^(2n-1).
  - Chamseddine, A. H., & Connes, A. (2008). *Why the Standard Model.*
    J. Geom. Phys. 58, 38.
  - Chamseddine, A. H., Connes, A., & Marcolli, M. (2007).
    Adv. Theor. Math. Phys. 11, 991. (Dimension count 2N² = 32 /
    4×4-grading argument.)
  - Connes, A. (1995); van Suijlekom, W. D. (2015). (KO-dim sign
    tables; Connes axiom framework.)

---

## One-line summary for forcing-catalogue update

> **Door 4f (2026-06-02): PARTIAL-DOOR FINAL (Outcome 3 NEUTRAL).** 12
> sub-tests on whether adopting Upgrade B (sphere-Lie-group axiom)
> introduces new downstream constraints: 0 new constraints, 0
> contradictions, 5 silent, 3 compatible, 2 identical to CCM, 1 clean
> extensibility rule (T11 via Adams 1958), 1 consistent strengthening
> of Bertrand × Hopf-tower reading (T12). Upgrade B is a clean
> reformulation of CCM's ℍ-selection with **1 axiom of net savings**
> (after honest accounting: Upgrade B 2 axioms total — Upgrade B + CCM
> second-order condition for D_F structure — vs CCM 3 axioms). No new
> predictive content for existing GeoVac observables. **Door 4 series
> structurally complete at PARTIAL-DOOR FINAL.** Remaining open
> questions are judgment-level (adopt Upgrade B vs CCM-3) and the
> deeper wall (force N_gen / inner KO-dim from a packing principle —
> deferred indefinitely; no known handle).
