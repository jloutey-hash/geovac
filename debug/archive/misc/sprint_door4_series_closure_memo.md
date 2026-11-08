# Sprint Door 4 series closure (4d → 4e → 4f) — 2026-06-02

*Closes the inner-algebra ℍ-vs-M₂(ℂ) fork at PARTIAL-DOOR FINAL via a
geometric handle (Door 4d), confirms via the literal AC-compatibility
test (Door 4e), and audits the downstream-constraint balance (Door 4f).
Net axiom savings vs standard CCM: one. Door 4 series structurally
complete.*

Three sub-sprint memos (canonical per-sub-sprint detail):
- `debug/door4d_division_algebra_sphere_memo.md`
- `debug/door4e_falsifier_A_prime_memo.md`
- `debug/door4f_falsifier_A_double_prime_memo.md`

Three sub-sprint drivers + data:
- `debug/door4{d,e,f}_*.py`
- `debug/data/door4{d,e,f}_*.json`

---

## 1. The arc in one paragraph

Door 4c (yesterday) closed NEGATIVE on the combined-J / KO-dim sign-table
handle for the ℍ-vs-M₂(ℂ) fork at the n=2 inner-factor rung. The fork
stayed open as ADMITTED-NOT-FORCED, with ℍ being a literature import
from Chamseddine–Connes 2008. Today's three sub-sprints introduce a
structurally orthogonal handle (Door 4d, the Division-Algebra-of-the-
Sphere criterion, DAS), test it against the standard AC tensor-product
compatibility conditions (Door 4e, literal Falsifier A'), and audit it
against the downstream-observable catalogue (Door 4f, Falsifier A'').
The arc lands at **PARTIAL-DOOR FINAL**: DAS via Upgrade B (sphere-Lie-
group axiom) is a clean reformulation of CCM ℍ-selection with one axiom
of net savings, no new predictive content, no contradictions. The
Door 4 series is structurally complete.

---

## 2. Sub-sprint 1: Door 4d — DAS introduced (PARTIAL-DOOR)

### Question
Is there a GeoVac-coherent handle on the ℍ-vs-M₂(ℂ) fork that is
structurally orthogonal to Door 4c's combined-J sign table?

### Answer
Yes — the Division-Algebra-of-the-Sphere criterion: "the inner algebra
at Hopf rung n is the natural associative real normed division algebra
whose unit-norm sphere IS the rung's Hopf bundle total space S^(2n-1)."

### Bit-exact verification (driver `debug/door4d_division_algebra_sphere.py`)
| Rung | Sphere | DAS algebra | Closure error | Verdict |
|:----:|:------:|:-----------:|:-------------:|:--------|
| n=1 | S¹ | ℂ | 1.11e−16 | S¹ = unit ℂ ✓ |
| n=2 | S³ | ℍ | 2.22e−16; SU(2)-iso 3.34e−16 | S³ = unit ℍ = Sp(1) = SU(2), direct ✓ |
| n=3 | S⁵ | (none) | — | Hurwitz: no candidate at dim 6 → fallback M₃(ℂ) |

### Verdict
**PARTIAL-DOOR.** DAS agrees with Door 4b's elimination forcings at
n=1 (ℂ) and n=3 (M₃(ℂ)) and closes the n=2 fork by selecting ℍ —
the direct (no-quotient) realization. Structurally orthogonal to
Door 4c (no J / sign-table anywhere). Not yet derivable from the
existing Paper 32 §VIII.B construction, which extracts only the gauge
group from the Hopf tower; adoption as a forcing requires the
sprint-reachable Falsifier A' theorem (does the AC tensor product
transfer the rung-n sphere's division-algebra structure to H_F?).

---

## 3. Sub-sprint 2: Door 4e — Falsifier A' literal test (PARTIAL-DOOR CONFIRMED)

### Question
Does the standard CCM AC tensor-product construction provide ANY
morphism that forces ℍ over M₂(ℂ) at n=2? Three natural compatibility
conditions tested.

### Results (driver `debug/door4e_falsifier_A_prime.py`)
| Condition | ℍ passes? | M₂(ℂ) passes? | Distinguishes? |
|:----------|:---------:|:-------------:|:--------------:|
| C1: SU(2)-equivariance under Ad action | ✓ (minimal) | ✓ (maximal) | **NO** |
| C2: Hopf U(1)-equivariance (scalar phase) | ✓ (trivially) | ✓ (vacuously) | **NO** |
| C3: Principal-bundle compatibility | ✓ (1.88e−15) | ✓ (1.34e−15) | **NO** |

### Substantive finding (not the negative outcome)
The literal Falsifier A' returning negative was *expected* (Door 4d
flagged it). The substantive content is the **import comparison**:
standard CCM uses three axioms beyond bare Connes data to select ℍ
(second-order condition + complex chiral fermion rep + 2N²=32 dimension
count); DAS uses one axiom (sphere-Lie-group). DAS is leaner.

### Verdict
**PARTIAL-DOOR CONFIRMED.** No AC compatibility condition fires. The
FULL-DOOR upgrade path is named: **Upgrade B (sphere-Lie-group axiom)**
— "the inner algebra at rung n is the *-algebra whose unit group IS
the rung-n sphere S^(2n−1) as a Lie group, when one exists; else
fallback M_n(ℂ)." Adopting Upgrade B is paper-level structural work,
not a theorem about existing CCM machinery.

---

## 4. Sub-sprint 3: Door 4f — Falsifier A'' downstream-constraint audit (PARTIAL-DOOR FINAL)

### Question
Does adopting Upgrade B introduce NEW downstream constraints on
existing GeoVac observables? Three possible outcomes:
- Outcome 1: New constraint the framework supports → corroborates DAS.
- Outcome 2: New constraint the framework contradicts → Upgrade B ruled out.
- Outcome 3: No new constraints → clean CCM-equivalent reformulation.

### Results (driver `debug/door4f_falsifier_A_double_prime.py`)
12 sub-tests, distributed across:

| Category | Count | Sub-tests |
|:---------|:-----:|:----------|
| SILENT (target unconstrained by Upgrade B) | 5 | T1 inner KO-dim; T2 γ_F; T3 Yukawa; T4 N_gen; T10 cross-rung CKM |
| COMPATIBLE (pre-existing constraint still passes for ℍ) | 3 | T5 order-zero; T6 order-one; T7 second-order |
| IDENTICAL to CCM at observable level | 2 | T8 spectral action gauge coefficient; T9 Higgs vacuum manifold |
| CLEAN EXTENSIBILITY RULE | 1 | T11 Adams 1958 |
| CONSISTENT STRENGTHENING of existing reading | 1 | T12 Bertrand × Hopf-tower |

**Net: 0 new constraints, 0 contradictions.** Outcome 3 realized cleanly.

### The substantive correction (T7 finding)
The 1-vs-3 axiom framing from Door 4e **over-counts the savings**. The
CCM second-order condition stays in BOTH formulations because it
constrains D_F's off-diagonal Yukawa structure, not A_F. So Upgrade B
replaces 2 of CCM's 3 ℍ-selection axioms (complex chiral rep + 2N²=32
dimension count) — NOT the second-order condition itself. **Net axiom
savings: 1, not 2.** This correction is now in Paper 32 §VIII.C as the
"Honest axiom-savings audit" paragraph.

### Verdict
**PARTIAL-DOOR FINAL.** Door 4 series structurally complete.

---

## 5. Aggregate forcing-ledger update

```
FORCED (Bertrand × Hopf-tower, unchanged from Door 4b):
  • factor count = 3
  • n=1 factor = C
  • n=3 factor = M_3(C)

CLOSEABLE via Upgrade B (per Door 4d/4e/4f, pending paper-level adoption):
  • n=2 factor = H (sphere-Lie-group axiom, 1 net axiom of savings vs CCM)

FREE (no construction reaches them; out of scope for the Door 4 series):
  • inner KO-dimension (Door 4f T1 confirms silent)
  • chirality grading γ_F (Door 4f T2 confirms silent; consistent but not unique)
  • Yukawa eigenvalues (seam theorem, Door 4f T3)
  • generation count N_gen (Door 4b Q3, Door 4f T4 confirms silent)
  • cross-rung Yukawa mixing / CKM (Door 4f T10)
```

---

## 6. Honest scope

### Closed at theorem grade
**Nothing.** The Door 4 series is structural-analysis work, not
theorem derivation. The forcings Door 4b establishes (n=1 = ℂ; n=3 =
M₃(ℂ); factor count = 3) remain theorem-grade per Bertrand × Hopf-tower
elimination, but those were closed before today's sprint.

### Structural / paper-level findings
- **DAS criterion (Door 4d)**: a structurally orthogonal handle on the
  ℍ-vs-M₂(ℂ) fork. Geometric / sphere-based; no J or sign-table data
  invoked.
- **Upgrade B (Door 4e)**: the sphere-Lie-group axiom as the FULL-DOOR
  upgrade path; "inner algebra at rung n = *-algebra whose unit group
  IS the rung-n sphere as a Lie group, when one exists; else fallback
  M_n(ℂ)." Underwritten by Adams 1958 / Hurwitz 1898.
- **Honest axiom-savings audit (Door 4f T7)**: Upgrade B replaces 2 of
  CCM's 3 ℍ-selection axioms but the CCM second-order condition stays
  in both formulations. Net savings: 1.

### Numerical observations (bit-exact spot-checks)
- S¹ = unit ℂ (closure 1.11e−16)
- S³ = unit ℍ = Sp(1) = SU(2) (closure 2.22e−16; SU(2)-iso 3.34e−16)
- S⁵ has no associative-division-algebra realization (Hurwitz, structural)
- AC compatibility conditions C1/C2/C3 all blind to ℍ vs M₂(ℂ)
  (machine precision)

These are *confirmations of expected structural facts*, not surprises.
No numerical-coincidence claim is made; nothing requires the
[[feedback_audit_numerical_claims]] audit pattern.

### Transcendentals
**None introduced.** The Door 4 series is purely algebraic / structural.
No π, ζ, G, or log appears anywhere in the bit-exact data, the
structural arguments, or the Paper 32 paragraph edits. The
[[feedback_tag_transcendentals]] discipline doesn't fire because no
transcendental appears.

### Named open follow-ons
1. **Paper-level adoption choice (judgment-level).** Should GeoVac
   adopt Upgrade B as a working axiom (axiom-minimality) or stay with
   CCM's three-axiom path (literature compatibility)? Not a probe-able
   falsifier; PI / community call.
2. **Deep wall (deferred indefinitely).** Force the generation count
   N_gen and the inner KO-dimension from a packing principle. Same as
   Door 4b's Falsifier B; no known handle.
3. **Optional retroactive softening of Door 4e memo's "1 vs 3" headline.**
   Not done today; the corrected accounting now lives in Paper 32 §VIII.C
   and the Door 4f memo. Door 4e memo is left as the historical record.

### What was NOT shown
- The 12 Falsifier A'' sub-tests are the natural set; an exhaustive
  enumeration might find a 13th where Upgrade B does introduce a new
  constraint. Outcome 3 is robust against the natural-test enumeration,
  not against arbitrary unexpected probes.
- No claim that Upgrade B IS the right axiom for GeoVac; Outcome 3
  says it doesn't break anything and isn't more predictive than CCM.
  The choice is judgment-level.
- No Yukawa, generation, KO-dim value selected anywhere. The H1 / W3 /
  Koide / probe-of-the-Higgs-Yukawa negatives in CLAUDE.md §3 are
  untouched. The seam theorem (Paper 32 §VIII Door 4) is untouched.

---

## 7. Files modified or created

### Created
- `debug/door4d_division_algebra_sphere.py`, `.json`, `_memo.md`
- `debug/door4e_falsifier_A_prime.py`, `.json`, `_memo.md`
- `debug/door4f_falsifier_A_double_prime.py`, `.json`, `_memo.md`
- This umbrella memo: `debug/sprint_door4_series_closure_memo.md`

### Modified
- `docs/forcing_catalogue.md` — three new entries under Graduation tests
  (Door 4d/4e/4f) updating the Door 4 status from "ℍ via literature
  import" to "ℍ closeable via Upgrade B at net 1-axiom savings."
- `CLAUDE.md` §2 — three one-liner entries (one per sub-sprint), to be
  consolidated into one umbrella entry per `/sprint-close` step 3.
- `CLAUDE.md` §1 version cursor: to be bumped from v3.43.0 to v3.44.0 per
  step 2 of `/sprint-close`.
- `papers/group1_operator_algebras/paper_32_spectral_triple.tex` §VIII.C
  — three new paragraphs after the existing Door 4c paragraph: "A
  geometric handle on the n=2 fork (Door 4d)", "Falsifier A' on DAS
  (Door 4e)", "Closure of the Door 4 series (Door 4f)" + the load-
  bearing "Honest axiom-savings audit" paragraph.

### Not modified
- No `geovac/` code touched. Topological integrity gate (§9 benchmarking
  rule) does not need to re-run; no risk of regression to the 18/18
  symbolic proofs or to existing accuracy benchmarks.
- No tests added or modified. (Per §13.4a, the Paper 32 paragraph edits
  add no new equations — they are descriptive paragraphs about the
  forcing-catalogue findings, not numerical claims requiring tests.)
- §3 dead-end table not modified (PARTIAL-DOOR is forward-progress, not
  a ruled-out approach). Hard prohibitions list (§13.5) untouched —
  Paper 2 combination rule conjectural label preserved; no fitted
  parameters introduced; natural geometry hierarchy unchanged.

---

## 8. The one-line summary for the forcing catalogue
(matches the entry now in `docs/forcing_catalogue.md` under
Graduation tests)

> **Door 4 series (4d/4e/4f), 2026-06-02:** PARTIAL-DOOR FINAL. The
> ℍ-vs-M₂(ℂ) fork is closeable via the **sphere-Lie-group axiom
> (Upgrade B)**: "the inner algebra at rung n is the *-algebra whose
> unit group IS the rung-n sphere as a Lie group, when one exists;
> else fallback M_n(ℂ)." Net axiom savings vs CCM: **1** (Upgrade B
> replaces 2 of CCM's 3 ℍ-selection axioms; the CCM second-order
> condition stays in both formulations to constrain D_F's off-diagonal
> Yukawa structure). 12-sub-test downstream-constraint audit returns
> Outcome 3 NEUTRAL (0 new constraints, 0 contradictions). Door 4
> series structurally complete; remaining open: paper-level adoption
> choice + the deep wall on N_gen / inner-KO-dim from a packing
> principle.
