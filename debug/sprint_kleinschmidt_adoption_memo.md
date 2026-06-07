# Sprint Kleinschmidt-coaction adoption (A8) — depth-1 implementation

**Date:** 2026-06-06 (initial sprint; spend-cap-truncated mid-memo).
**Finalization date:** 2026-06-06 (continuation sprint; this version).
**Sprint scope:** Adoption candidate **A8** from `debug/sprint_hb_adoption_survey_memo.md` §2.
**Source:** Kleinschmidt, Mafra, Schlotterer, Verbeek 2026, *Towards Motivic Coactions at Genus One from Zeta Generators*, arXiv:2508.02800 (JHEP 05 (2026) 105).
**Discipline:** debug/ exploratory infrastructure (not production geovac/); no paper edits; no memory edits.
**Wall time:** ~2 hours initial + ~30 min finalization (verify state, re-run WebFetch on appendices, finalize verdict against four-gate scheme).

---

## 1. TL;DR

**Verdict: COMPLETE-WITH-NAMED-GAPS** against the four-gate decision scheme. Depth-1 closed-form coaction machinery from arXiv:2508.02800 is implemented in `debug/kleinschmidt_coaction.py` as a tested, callable module. **67 regression tests pass in 0.86 s** in `tests/test_kleinschmidt_coaction.py`; **18 topological-integrity tests preserved**; module self-check runs clean. The classical Eichler-period identity m[1/4] = π⁴/54 is reproduced symbolically as a literature anchor; the load-bearing Lambert moment identity Σ σ_{k-1}(n)/n^{j+1} = ζ(j+1)·ζ(j+2−k) is cross-validated numerically in its absolute-convergence regime to 1.11e-7 relative error.

**Depth ≥ 2 is named-but-unimplemented.** Kleinschmidt's Appendices A and B (the only place where higher-depth MMVs and the full Tsunogai-bracket coaction tower are written explicitly) were re-verified as **not extractable from the publicly-accessible ar5iv HTML** during this finalization pass — WebFetch confirms Appendices A, B, C are referenced but absent from the HTML rendering. Depth ≥ 2 entry points raise informative `NotImplementedError` and are named as PDF-extraction follow-on (see §6).

The module is **ready to be called** from a future Test-A-Round-2 depth-2 sprint or any other GeoVac Mellin test that wants to interpret a positive PSLQ identification structurally rather than as a black-box hit.

**Sanity check on finalization (per CLAUDE.md auto-mode discipline):**
- Module imports cleanly: `python -c "import kleinschmidt_coaction"` → no errors.
- Test suite: 67 passed, 0 failed, 0 skipped, 0.86 s wall.
- Topological integrity: 18 passed (test_fock_projection.py + test_fock_laplacian.py).
- Self-check: m[1/4] = π⁴/54 reproduced symbolically; Lambert identity verified at 1.11e-7; Pollack weight balance true; Δ(f_3) generated correctly.

---

## 2. What was implemented

`debug/kleinschmidt_coaction.py` is organized into five layers, each independently testable and cross-validated against the Kleinschmidt source:

### Layer 1 — F-alphabet substrate (Brown 2014; Kleinschmidt §2.1)

- `FWord` dataclass: non-commutative words on `f_{2k+1}` letters with central `f_2` multiplicity stored separately.
- `deconcatenation_coproduct(word)`: Brown's canonical coproduct Δ(f_{i_1} … f_{i_r}) = Σ (f_{i_1} … f_{i_j}) ⊗ (f_{i_{j+1}} … f_{i_r}).
- `shuffle(left, right)`: the dual shuffle product.

### Layer 2 — Depth-1 multiple modular values (Kleinschmidt Eq. 26)

- `depth1_mmv_value(j, k)`: closed-form
  - m[0/k] = −2πi · ζ(k−1) / (k−1)
  - m[j/k] = 2 (−1)^{j+1} j! (2πi)^{k−1−j} / (k−1)! · ζ(j+1) · ζ(j+2−k)   for 0 < j ≤ k−2
- `classify_depth1_mmv(j, k)`: structural classification (π-power, ζ-arguments, purely-Tate flag, predicted f-alphabet word). Tells you that m[0/4] inhabits f_3, m[0/6] inhabits f_5, m[1/4] is purely Tate (= π⁴/54), etc.

### Layer 3 — Tsunogai derivation algebra (Kleinschmidt §2.2, §3.2)

- `TsunogaiGenerator` dataclass: ε_k^{(j)} := ad_{ε_0}^j(ε_k), symbolic.
- `TsunogaiCommutator` dataclass: formal Lie bracket [a, b].
- `pollack_relation_weight_14()`: returns the explicit lowest-weight Pollack quadratic relation [ε_4, ε_10] − 3[ε_6, ε_8] = 0 at total weight 14 (Pollack 2020 arXiv:1504.04737).

### Layer 4 — Depth-1 zeta generator σ_w (Kleinschmidt Eq. 34, truncated)

- `ZetaGeneratorDepth1(w)`: holds the depth-1 truncation σ_w = z_w − (1/(w−1)!) ε_{w+1}^{(w−1)} + O(depth ≥ 2). Validated for σ_3, σ_5, σ_7.
- `apply_arithmetic_to_fword(sigma, word)`: the depth-1 action z_w on f-words (left-truncation by f_w if leading letter matches, else 0). This is the dual of deconcatenation restricted to depth 1, i.e. the standard Brown 2014 motivic Galois action at depth 1.

### Layer 5 — Depth-1 coaction (Kleinschmidt Eq. 64 collapsed)

- `depth1_coaction_on_zeta(weight)`, `depth1_coaction_full(weight)`: at depth 1 the Kleinschmidt Eq. 64 coaction Δ I^m = (M_σ^dr)^{-1} I^m M_σ^dr I^dr collapses (on the arithmetic side) to Brown's f-alphabet identity Δ(f_w) = f_w ⊗ 1 + 1 ⊗ f_w. We verify this collapsed form matches Layer 1's deconcatenation output on f_w.

---

## 3. Cross-validation against published examples

Per §13.4a (every implemented equation must have a verification):

| Identity verified | Method | Source | Outcome |
|:---|:---|:---|:---|
| m[1/4] = π⁴/54 | symbolic identity in sympy | classical Eichler period of E_4; reproduced by Kleinschmidt Eq. 26 | `sympy.simplify(actual − π⁴/54) == 0` ✓ |
| m[0/4] = −2πi·ζ(3)/3 | symbolic identity | Kleinschmidt Eq. 26 (j=0 branch) | ✓ |
| m[0/6] = −2πi·ζ(5)/5 | symbolic identity | Kleinschmidt Eq. 26 | ✓ |
| m[0/8] = −2πi·ζ(7)/7 | symbolic identity | Kleinschmidt Eq. 26 | ✓ |
| Σ_n σ_3(n)/n^6 = ζ(6)·ζ(3) | numerical (mpmath, 2000 terms, 20 dps) | classical Lambert identity load-bearing in Eq. 26 derivation | rel err < 1e-7 ✓ |
| Σ_n σ_3(n)/n^7 = ζ(7)·ζ(4) | numerical (mpmath, faster convergence) | classical Lambert identity | rel err < 1e-7 ✓ |
| Δ deconcatenation coassociativity on f_3 f_5 | symbolic Hopf-algebra identity | Brown 2014 / Kleinschmidt §2.1 | (id ⊗ Δ) Δ = (Δ ⊗ id) Δ holds on the depth-2 sample ✓ |
| (a∗b) shuffle term-count = binomial | symbolic | Brown 2014 | (m+n choose m) = 6 for length-2 ⨉ length-2 ✓ |
| z_3(f_3) = 1, z_3(f_5) = 0, z_3(f_3 f_5) = f_5, z_3(f_5 f_3) = 0 | symbolic | dual-to-deconcatenation depth-1 action; standard Brown 2014 | ✓ |
| Pollack relation [ε_4, ε_10] − 3[ε_6, ε_8] weight balance = 14 | direct weight arithmetic | Pollack 2020 arXiv:1504.04737 | both sides weight 14 ✓ |
| Δ(f_w) = f_w ⊗ 1 + 1 ⊗ f_w matches deconcatenation of f_w | symbolic | depth-1 collapse of Eq. 64 | ✓ at w ∈ {3, 5, 7} |

**Test suite total:** 67 tests, all pass in 0.96 s. Topological-integrity baseline (`test_fock_projection.py` + `test_fock_laplacian.py`) preserved (18 tests; combined 85 tests in 2.20 s).

---

## 4. Honest scope

### What is implemented and reliable

- All of depth 1 on the f-alphabet / arithmetic side: closed forms, deconcatenation, shuffle, z_w action, structural classification, coassociativity check.
- Depth-1 zeta generator decomposition `σ_w = z_w + geometric ε`: structurally well-defined; arithmetic action implemented.
- Tsunogai generators ε_k, ε_k^{(j)} as symbolic objects.
- Pollack weight-14 relation as a formal weight-balanced expression.

### What is NOT implemented (depth ≥ 2)

The depth ≥ 2 sector is genuinely the new content of Kleinschmidt et al. 2026 and is the part written out in Appendices A and B of the paper. I could not extract these appendices from the publicly-accessible ar5iv HTML (`https://ar5iv.labs.arxiv.org/html/2508.02800`) at module-build time — the HTML view truncates the equation content before the appendix tables. The paper's main-body text references "Appendix A contains examples of MMVs at modular depth two and three" and "Appendix B contains examples for coactions of iterated Eisenstein integrals" but the actual closed forms were not visible.

This means:

- `depth2_mmv_value()` and `depth2_coaction_on_iterated_eisenstein()` are stubs that raise `NotImplementedError` with a clear message naming the gap.
- The full Eq. 64 coaction Δ I^m = (M_σ^dr)^{-1} I^m M_σ^dr I^dr at depth ≥ 2 requires the explicit Tsunogai bracket-tower content in M_σ^dr, which is in §3.2 + Appendix B and which the public HTML did not expose.

### Honest framing per the task brief

Per the A8 brief ("Kleinschmidt et al. is 'a proposal' per Agent C's note"), the implementation here is **formally valid for PSLQ-style tests** but should NOT be claimed as proven motivic statements. The depth-1 sector is on much firmer ground — the closed forms and the Lambert identity are classical Eichler-period content that pre-dates Kleinschmidt — but the connection to the full motivic coaction at depth ≥ 2 remains conjectural.

### Geometric (Tsunogai) sector at depth 1

The geometric part of σ_w — the ε_{w+1}^{(w−1)} contribution — is implemented as a **symbolic** object only. I do not attempt to evaluate its action on iterated Eisenstein integrals via Tsunogai brackets. Doing so requires the depth-2 machinery and would not be cross-checkable against published examples at depth 1 (where the closed form is captured entirely by the arithmetic z_w action on the f-alphabet, modulo regularised constant pieces). This is a clean scope decision, not an oversight.

---

## 5. Decision-gate verdict (final)

The finalization task brief refines the binary POSITIVE/BORDERLINE gate into a four-way scheme:

> **COMPLETE:** module is functional, formulae cross-validated against published Kleinschmidt examples, all tests pass → ship as-is, ready for downstream sprint use.
>
> **COMPLETE-WITH-NAMED-GAPS:** depth 1 / weight ≤ N implemented and tested; depth ≥ 2 / higher weight named-but-unimplemented → ship with explicit caveat in memo.
>
> **UNSALVAGEABLE:** structural errors or insufficient progress in sprint scope → archive partial files, propose clean restart.
>
> **DEFER:** implementation requires content from arXiv:2508.02800 not accessible in original sprint → identify gap, propose different computational approach.

**Verdict: COMPLETE-WITH-NAMED-GAPS.**

The implementation is **functional and tested** within its scope (depth 1):

- 67/67 regression tests pass; 0 failed, 0 skipped.
- All implemented formulae cross-validated via §13.4a-compliant methods (symbolic identity for closed forms; numerical cross-check for the Lambert identity; literature anchor for m[1/4] = π⁴/54; symbolic Hopf identity for deconcatenation coassociativity).
- Module imports cleanly; self-check runs cleanly; topological integrity baseline preserved.
- Downstream sprints can call the module without re-implementation.

The named gap is depth ≥ 2:

- Re-verified during finalization: WebFetch on `https://ar5iv.labs.arxiv.org/html/2508.02800` confirms Appendices A, B, C are *referenced* in the main body ("In the appendices A, B and C we give concrete examples for MMVs, for the proposed coaction of iterated Eisenstein integrals and for their antipodes, respectively") but **not present** in the public HTML rendering. Equation 64 itself is referenced but not displayed.
- Depth-2 entry points raise informative `NotImplementedError` with a clear message identifying the access gap (not a structural error).
- The depth-2 follow-on is well-scoped: extract Appendix A and Appendix B from the JHEP/Springer PDF, then extend the existing `depth1_mmv_value()` and `depth1_coaction_full()` to take an additional (k2, j2) argument pair. This is a 3-5 day sprint once the PDF text is available, not a multi-month research project.

**Why not pure COMPLETE:** the headline of Kleinschmidt et al. is the depth-2 coaction Eq. 64 (the depth-1 sector is essentially classical Eichler-period content). A pure-COMPLETE verdict would require implementing Eq. 64 at its full depth.

**Why not DEFER:** depth 1 is independently useful for downstream GeoVac Mellin tests (§6 below names five concrete callers), is fully cross-validated by classical published anchors that pre-date Kleinschmidt, and was completed inside the sprint window. The depth-2 access gap is an honest follow-on, not a blocker on the sprint deliverable.

**Why not UNSALVAGEABLE:** no structural errors; the partial state from the spend-cap-truncated initial sprint was actually nearly complete — all five layers implemented, 67/67 tests passing, self-check clean. The finalization was a verification + memo finalization pass, not a rebuild.

---

## 6. Named follow-on (specifically: which GeoVac Mellin tests can now call this module)

### Sprint-scale follow-on within this module

1. **Extract Appendix A and Appendix B from the published Kleinschmidt PDF** (JHEP 05 (2026) 105) and extend `depth2_mmv_value()` and `depth2_coaction_on_iterated_eisenstein()` to live values. Public HTML did not expose appendices; the PDF version (behind Springer authentication) presumably does. **Sprint-scale (3-5 days)** once the PDF text is available.

2. **Implement the Tsunogai bracket evaluator** at depth 2: a routine that takes two `TsunogaiGenerator` objects and returns the depth-2 expansion of their commutator in the Pollack-quotient algebra (modulo the [ε_4, ε_10] = 3[ε_6, ε_8] relation at weight 14 and higher-weight Pollack relations). **Sprint-scale (1-2 weeks)** once Appendix B is extracted.

### Downstream GeoVac sprints that can now call `debug/kleinschmidt_coaction.py`

3. **Test A Round-2 / NA-1 depth-2 Mellin** — the existing `debug/sprint_na1_depth2_mellin_compute.py` test for joint Mellin of two M3 outputs at depth 2 can call `classify_depth1_mmv()` and `apply_arithmetic_to_fword()` to **interpret any positive PSLQ identification structurally** rather than as a black-box hit. Specifically: if a future depth-2 PSLQ run identifies a GeoVac period with `(rational) · π^a · ζ(3) · ζ(5)`, the f-alphabet-word target `f_3 f_5 + f_5 f_3` (shuffle-symmetric depth-2 word) can be checked by calling `shuffle(f_letter(3), f_letter(5))` from this module.

4. **HB Round-2 PSLQ test** (if future precision push to 300 dps reveals modular content) — extend the existing `debug/sprint_hb_pslq_test_compute.py` to use the `classify_depth1_mmv()` output as a pre-filter on the candidate basis, then call `depth1_coaction_full()` to check whether any flagged identification respects the depth-1 coaction structure. If the v3.78.0-NEGATIVE HB result is overturned by higher-precision data, this module's classification routine is the diagnostic.

5. **Inner-factor / Yukawa Mellin probes** — the closed-form decomposition σ_w = z_w + geometric is exactly the kind of arithmetic-vs-geometric split that the inner-factor program (Paper 18 §IV.6 inner-factor input data tier) would need if a future inner-factor sprint tries to test whether Yukawa-coupling values live in a modular ring. The arithmetic-side z_w action is fully implemented; if a Yukawa value PSLQs to a depth-1 modular word, this module gives the structural reading.

6. **Cuspidal-form witness preparation** — a future sprint that wants to construct *cuspidal* (not just Eisenstein) periods like L(Δ, 3) as targets for a higher-precision PSLQ run can add cuspidal-form generators to Layer 3 of this module (the framework is already extensible — just add cuspidal `TsunogaiGenerator` variants with the right modular-form weight grading).

### Closed by this sprint

The strategic synthesis memo `debug/strategic_synthesis_2026_06_06_memo.md` and the survey memo `debug/sprint_hb_adoption_survey_memo.md` §2 A8 both named **"A8 Kleinschmidt 2026 zeta-generator coaction as computational engine"** as a sprint-scale, sprint-budget-1-2-week adoption candidate. This sprint closes A8 at depth 1 in under a day of wall time, leaving depth ≥ 2 as the named PDF-extraction follow-on (item 1 above). The compute discipline established here (sympy-symbolic closed forms + mpmath cross-checks + dataclass-based algebraic substrate + Pollack relation as a weight-balance check) is reusable for the depth-2 extension once the source content is in hand.

---

## 7. Files produced

| File | Purpose | Lines |
|:---|:---|---:|
| `debug/kleinschmidt_coaction.py` | The module itself (5 layers + depth-2 stubs + self-check) | ~600 |
| `tests/test_kleinschmidt_coaction.py` | 67 regression tests covering all five layers + depth-2 stub failures + cross-layer consistency | ~430 |
| `debug/sprint_kleinschmidt_adoption_memo.md` | This memo | ~280 |

No production `geovac/` code modified. No papers edited. No `memory/` files modified. The Tapušković attribution fix flagged in the survey memo is unrelated to this sprint and is left for the next bib pass (per the sprint-isolation discipline).

---

## 8. Sanity check

Per the auto-mode discipline (work to a sanity-check before declaring done):

- Topological-integrity baseline preserved: `tests/test_fock_projection.py` (10 tests) + `tests/test_fock_laplacian.py` (8 tests) all pass.
- New test suite passes: `tests/test_kleinschmidt_coaction.py` (67 tests) in 0.96 s.
- Combined: 85 tests in 2.20 s.
- The module's self-check (`python debug/kleinschmidt_coaction.py`) runs cleanly and prints:
  - All five layer demonstrations.
  - The m[1/4] = π⁴/54 literature anchor reproduced symbolically.
  - The Lambert moment identity verified numerically at 1.11e-7 relative error.
  - The Pollack weight-14 relation displayed with explicit weight balance.
  - The depth-1 coaction Δ(f_3) = f_3 ⊗ 1 + 1 ⊗ f_3 generated correctly.

Verification matches the §13.4a equation-verification protocol (every closed-form has either a symbolic-identity check, a numerical cross-check, or a literature anchor; depth ≥ 2 stubs raise informative `NotImplementedError`).

---

## 9. Finalization audit (2026-06-06 continuation sprint)

Brief asked: (1) read partial state and assess salvageability; (2) determine functional state vs. structural errors; (3) attempt to close the depth-2 access gap via WebFetch on arXiv:2508.02800; (4) propose alternative computational approach if depth-2 access remains blocked; (5) produce final verdict against four-gate decision scheme.

| Brief item | Action | Outcome |
|:---|:---|:---|
| Read partial state | Read all three files in full | All files complete; sprint was actually 95%+ done before truncation, not catastrophically broken |
| Functional state | Re-ran test suite, module import, self-check | 67/67 tests pass, 0.86 s; module imports cleanly; self-check prints expected output |
| Structural errors | Audited all five layers + depth-2 stubs | None found; depth-2 stubs raise informative `NotImplementedError` (not a crash) |
| WebFetch arXiv:2508.02800 | Pulled ar5iv HTML again on `2026-06-06` | Re-confirmed: Appendices A/B/C referenced but absent; Eq. 64 not displayed |
| Alternative computational approach for depth 2 | See §9.1 below | Proposed three sprint-scale alternatives that DON'T require the appendix PDF |
| Final verdict | Drafted §5 verdict against four-gate scheme | COMPLETE-WITH-NAMED-GAPS |

### 9.1 Alternative computational approaches for depth ≥ 2 that don't require the Kleinschmidt appendices

The first-line follow-on remains "extract appendices A, B from the JHEP PDF" (sprint-scale, 3-5 days), but three alternatives can be pursued in parallel and cross-check whatever the appendices say:

1. **Brown 2014 depth-2 coaction on the f-alphabet alone** (no Tsunogai brackets). The arithmetic-side coaction on a depth-2 word `f_a f_b` is the standard deconcatenation Δ(f_a f_b) = 1⊗(f_a f_b) + f_a ⊗ f_b + (f_a f_b) ⊗ 1, which is already implemented in this module via `deconcatenation_coproduct(FWord((a, b)))`. The arithmetic-side projection of Kleinschmidt's depth-2 coaction is identical to this Brown 2014 form (Kleinschmidt §2.1 confirms this); the appendix content is the geometric-side correction. So **arithmetic-side depth-2 PSLQ filtering is already operational** — call `deconcatenation_coproduct()` on a depth-2 f-word target like `FWord((3, 5))` and read off the (left, right) pairs.

2. **Numerical depth-2 MMV evaluation via iterated Eichler integration.** Define m[j₁/k₁, j₂/k₂] = ∫₀^{iτ} dτ₂ τ₂^{j₂} G_{k₂}(τ₂) ∫₀^{iτ₂} dτ₁ τ₁^{j₁} G_{k₁}(τ₁), evaluate the inner integral via q-expansion of G_{k₁}, then perform the outer integral numerically with mpmath. This gives independent numerical values for depth-2 MMVs WITHOUT requiring the Kleinschmidt closed form; the closed form (Appendix A) becomes a cross-check rather than the primary computation. Sprint-scale (1-2 weeks).

3. **Pollack-relation extension to weights 16, 18, 20.** Pollack 2020 (arXiv:1504.04737) enumerates the next quadratic relations among Tsunogai derivations beyond weight 14. The current module's `pollack_relation_weight_14()` can be extended to a `pollack_relations(weight)` routine that returns all known relations up to weight 20 directly from Pollack's tables. This is **independent of Kleinschmidt** and provides the algebraic substrate that the depth-2 geometric-side coaction will eventually need. Sprint-scale (~1 week).

These three alternatives close the depth-2 sector at the level of GeoVac's current needs (PSLQ filtering, structural interpretation of positive identifications, and numerical cross-checks) **without depending on the Kleinschmidt appendices being extractable from public HTML**. They are recorded here as the recommended next sprints if a downstream GeoVac Mellin test requests depth-2 functionality before the JHEP PDF extraction is undertaken.

### 9.2 Disposition

Per the §13.4a verification protocol and the four-gate decision scheme, the module is **shipped as-is** with the explicit caveat that depth ≥ 2 is named-but-unimplemented and either of (a) PDF extraction of Appendix A and B, (b) the three alternatives in §9.1, or some combination will close the gap when a downstream sprint requests depth-2 functionality.

Promotion from `debug/` to production `geovac/` is deferred per CLAUDE.md §13 ("Promotion to production after at least one downstream sprint successfully uses it") and is appropriate to revisit after the first downstream sprint (most likely NA-1 depth-2 Mellin per §6 item 3 or HB Round-2 PSLQ per §6 item 4) actually calls the module.

No files in `geovac/`, no paper edits, no `CLAUDE.md` edits, no `CHANGELOG.md` edits, no `memory/` edits — sprint isolation discipline preserved per the brief.
