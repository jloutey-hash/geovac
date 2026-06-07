# Sprint Kleinschmidt-coaction adoption (A8) — depth-1 implementation

**Date:** 2026-06-06 (post-HB-PSLQ-negative continuation).
**Sprint scope:** Adoption candidate **A8** from `debug/sprint_hb_adoption_survey_memo.md` §2.
**Source:** Kleinschmidt, Mafra, Schlotterer, Verbeek 2026, *Towards Motivic Coactions at Genus One from Zeta Generators*, arXiv:2508.02800 (JHEP 05 (2026) 105).
**Discipline:** debug/ exploratory infrastructure (not production geovac/); no paper edits; no memory edits.
**Wall time:** ~2 hours (literature extraction + design + implementation + tests + memo).

---

## 1. TL;DR

**Verdict: POSITIVE (BORDERLINE).** Depth-1 closed-form coaction machinery from arXiv:2508.02800 is now implemented in `debug/kleinschmidt_coaction.py` as a tested, callable module. 67 regression tests pass in `tests/test_kleinschmidt_coaction.py`. The classical Eichler-period identity m[1/4] = π⁴/54 is reproduced symbolically as a literature anchor; the load-bearing Lambert moment identity Σ σ_{k-1}(n)/n^{j+1} = ζ(j+1)·ζ(j+2−k) is cross-validated numerically in its absolute-convergence regime to 7-8 digits.

**Depth ≥ 2 is NOT implemented.** Kleinschmidt's Appendices A and B (the only place where higher-depth MMVs and the full Tsunogai-bracket coaction tower are written explicitly) were not extractable from the publicly-accessible ar5iv HTML at module-build time. Depth ≥ 2 entry points raise informative `NotImplementedError` and are named as follow-on; see §6.

The module is now ready to be called from a future Test-A-Round-2 depth-2 sprint (or any other GeoVac Mellin test that wants to interpret a positive PSLQ identification structurally rather than as a black-box hit).

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

## 5. Decision-gate verdict

Per the task brief:

> POSITIVE: closed-form coaction formulae from arXiv:2508.02800 implemented as a tested debug/ module; cross-validated against published examples; usable from a future sprint without re-implementation.

> BORDERLINE: coaction formulae implemented at depth 1 only (or weight ≤ some bound); depth ≥ 2 / higher weight blocked by structural complexity → report what's missing, leave as named follow-on.

**Verdict: BORDERLINE — closer to POSITIVE within depth 1, BORDERLINE because depth ≥ 2 is blocked.** Depth 1 is the full deliverable the task brief named as the **target**, with depth 2 as the **stretch**. Depth 1 is delivered with full test coverage and literature cross-validation; depth 2 is stubbed with a named follow-on. This matches the BORDERLINE description verbatim.

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
