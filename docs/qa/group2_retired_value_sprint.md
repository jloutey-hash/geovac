# group2 — retired-value remediation sprint (scope)

> # COMPLETE — 2026-09-01. All 24 remediated; `/qa group2` unblocked.
>
> `check_numeric_consistency.py --gate group2` → **PASS**, and all eleven
> targets now pass. Escape / internal-title / test-ref / file-ref gates green
> on group2; all four papers compile with zero undefined references.
>
> **The material finding: Conjecture 1 no longer holds as stated.** Paper 19
> framed it as "≲1,500 Pauli terms, roughly 4.5× the composed count of 334" —
> so its *content* is a ratio and its absolute figure was derived from a count
> that has since been retired. Under the exact rule the **ratio leg survives**
> (2,726 / 838 = 3.25×, comfortably under 4.5×) while **both absolute legs
> fail**: 2,726 exceeds the ~1,500 the ratio implied, and exceeds the STO-3G
> count of 907. The "below Level 4N 3,288" leg survives, but by 17% rather
> than 3.7×. Three places called it "confirmed"; all three now carry the split
> verdict, and the conjecture is restated in its rule-independent ratio form.
> This is exactly what a numeral swap would have destroyed.
>
> Also fixed en route, all **pre-existing**: eight unresolvable
> cross-document `\ref`s in Paper 17 (labels living in Papers 19 and 34 —
> LaTeX renders them `??`; replaced with the target section titles), three
> dangling `\ref`s and two bad `\cite`s in Paper 13, and the `2659` assertion
> in `tests/test_sturmian_qubit.py` that skipped behind `@pytest.mark.slow`.
>
> ---
>
> **PROGRESS LOG — 4 of 24 (Papers 8 and 13), then completed below.**
>
> **Paper 8 (3 loci) — DONE.** Both of its claims were re-checked, not just
> the numerals, using `debug/qa/remeasure_paper8_bu2.py`:
> - ERI counts `65 / 1492` → **`107 / 3800`**, and the *basis-independence*
>   claim survives **exactly** — Sturmian and standard both measured, both
>   give the same counts.
> - the Track BU-2 comparison re-measured end to end: `112 vs 120` →
>   **`280 vs 288`** (the Pauli advantage shrinks, ~7% → ~3%);
>   `31.4 vs 11.3` → **`26.3 vs 11.2`** (2.8× → 2.4×);
>   `2627 vs 2659` → **`14,047 vs 14,079`**;
>   `349 vs 78` → **`242 vs 74`** (4.5× → 3.3×).
>   The Löwdin-inflation conclusion stands, with smaller factors.
> - `tests/test_sturmian_qubit.py` still asserted the retired `2659` **behind
>   `@pytest.mark.slow`**, so it skipped on every default run while the
>   n_max=2 case beside it had been corrected. Updated to `14079` and
>   un-marked (the build is ~2 s). Same "guard that never runs" pattern as the
>   balanced-λ test.
>
> **Paper 13 (1 locus) — DONE**, plus five **pre-existing** defects found on
> the way: three dangling `\ref`s (`sec:liouville`, `sec:dboc` — both
> subsections exist but carried **no label**, so the labels were added rather
> than the references repointed; `sec:coupled_channel` — no such section,
> repointed to `sec:radial`) and two `\cite{paper15}` with no bibitem (wrong
> key + missing entry). Paper 13 now compiles with zero undefined references.
> Its scaling sentence also carried a **second** retired number the list did
> not flag — `O(Q^{4.60})` for Gaussians, outside Paper 14's measured
> `Q^{3.9}`–`Q^{4.3}` band — so the comparison was re-framed: the advantage is
> **0.1–0.5 in exponent**, not the far larger margin the retired pair implied.
>
> **Spilled outside group2: BeH₂'s balanced λ was also at LiH's geometry.**
> The 2026-08-31 fix corrected the second/third-row table; BeH₂ lives in a
> different table and was missed. Corrected to its own R = 2.502 bohr
> (`328.110 → 323.427` incl-identity, `289.581 → 286.893` non-identity;
> Pauli 8,868 and QWC 1,220 unchanged), across **Papers 14 and 20** (four
> loci, including a `0.88× → 0.86×` ratio and a `12% → 14%`) and the registry,
> with the wrong-geometry values retired. H₂O and LiH are unaffected — H₂O's
> nuclei come from a fixed function that never took R, and LiH's own R *is*
> 3.015.
>
> **Remaining: 20 loci — Paper 17 (4) and Paper 19 (16).** Values already
> measured and ready to apply:
>
> | locus | retired | corrected |
> |:------|:--------|:----------|
> | P17 L711, P19 L864/874 BeH₂ 1-norms | 304.7 vs 354.9 | **323.4 vs 374.9** (both incl-identity) |
> | P17 L763 LiH balanced regression | 878 | **2,726** |
> | P17 L1670 Löwdin inflation | 1,711 vs 120 | **4,743 vs 288** (16.5×, was 14.3×) |
> | P17 L1740 composed vs full N-electron | 334 vs 3,288 | **838 vs 3,288** — the 3,288 reproduces *exactly* (full 4e, l_max=2, Q=10) and is **not** retired; only the composed side was. The "20× lower 1-norm" also survives (738.5/38.3 ≈ 19.3×). |
> | P19 ×10 LiH balanced Pauli | 878 | **2,726** |
> | P19 ×5 LiH composed Pauli | 334 | **837** non-identity / **838** with |
> | P19 L865 H₂O composed Pauli | 778 | **1,953** non-identity / **1,954** with |
>
> Paper 19's resource table (lines ~863–865) is the one piece still needing a
> full re-measurement rather than substitution, per §4(b) below.

**Scoped 2026-08-31** from the Phase-0 C21 widening
(`debug/qa/phase0_gate_audit_notes.md` §F5). This is a **prerequisite** for
`/qa group2`: sending a panel at a paper set that still asserts the retired
pair-diagonal convention wastes the panel on defects a gate already names.

Not started. No group2 paper has been edited.

---

## 1. What this is

The 2026-08-29 exact-rule ERI correction re-priced the composed/balanced
resource numbers corpus-wide. It swept the loci C21 was watching — which at
the time meant **group3, group4 and group6 only**. group2 was certified
2026-06-28, two months before the correction, and was never in a numeric
gate's scope. It therefore still carries the pair-diagonal vintage.

Widening C21 to all eleven targets surfaced **24 live retired values** in four
group2 papers. Every one was spot-checked in context: they are present-tense
assertions, not disclosed history.

## 2. The 24, by paper

| paper | count |
|:------|------:|
| `paper_19_coupled_composition.tex` | **16** |
| `paper_17_composed_geometries.tex` | 4 |
| `Paper_8_Bond_Sphere_Sturmian.tex` | 3 |
| `paper_13_hyperspherical.tex` | 1 |

## 3. The 24, by retired quantity

| n | quantity | retired → canonical |
|--:|:---------|:--------------------|
| 10 | `lih_balanced_pauli` | 878 → **2,726** |
| 5 | `lih_composed_pauli` | 334 → **837** non-identity (**838** identity-in) |
| 3 | `beh2_composed_pk_lambda` | 354.9 → **374.866** |
| 2 | `he_n2_pauli` | 120 → **287** |
| 1 | `he_n3_pauli` | 2,659 → **14,078** |
| 1 | `block_sp_n2_eri` | 65 → **107** |
| 1 | `exp_pauli_4pt` | 3.15 → **3.773** |
| 1 | `h2o_composed_pauli` | 778 → **1,953** |

Regenerate the live list at any time with:

```
python debug/qa/check_numeric_consistency.py --gate group2
```

## 4. Why this is not a find-and-replace

Three reasons, in increasing order of how much work they add.

**(a) The prose around each number is a claim keyed to its magnitude.**
CLAUDE.md §15 rule 2. Paper 19 line 151 currently reads

> "854 Pauli terms (2.56× the composed value of 334)"

Swapping 334 → 837 leaves the ratio 2.56× wrong and the sentence
self-contradicting. Every locus needs its sentence re-derived, not its numeral
replaced.

**(b) Paper 19 carries a whole resource table.** Lines ~863–865 tabulate
LiH / BeH₂ / H₂O composed-vs-balanced Pauli counts, ratios and 1-norms — every
cell pair-diagonal vintage, and the ratio column is derived from the two count
columns. The table must be re-measured from the live builders, not patched
cell by cell. This is the same shape as the Stage-4 ledger work, and the
`debug/qa/` drivers from that pass are reusable.

**(c) Conventions.** `lih_composed_pauli` is 837 non-identity / 838
identity-inclusive. Which one a locus wants depends on what that table's other
columns do — the identity-in/identity-out mix is precisely the defect class
C21 check E exists for. Decide per table, not per number.

**(d) Two loci are outside the retired-value list but in the same family.**
The Phase-0 remediation of the *other* four targets found a `190×` cc-pVDZ
ratio that had **no registry entry at all** (canonical 76×). It is now
registered, so C21 will catch any group2 instance — but the lesson stands:
re-read each passage for retired content beyond the flagged numeral.

## 5. Suggested shape

1. Re-measure Paper 19's resource table from the live builders; record the run
   in a memo alongside the Stage-4 ledger notes.
2. Work paper by paper, largest first (19 → 17 → 8 → 13), re-reading each
   sentence and fixing derived ratios alongside their inputs.
3. Register anything found that has no entry (the 190× lesson), proving
   two-way discrimination per the C21 maintenance rule.
4. Re-run `check_numeric_consistency.py --gate group2` to zero, plus C19
   (escape corruption) and C10 (compiles) on the four papers.
5. Widen `tests/test_numeric_registry.py::test_gate_passes_on_current_corpus`
   to all eleven targets once group2 is clean — it currently parameterises
   only over C21's legacy three-target map, so none of this is
   regression-protected.

## 6. Then, and only then

`/qa group2` — as a DELTA run scoped to the diff, per
`docs/qa/recert_sweep_plan.md`. The plan's Phase 4 placement for group2 is
**superseded**: group2 is a Layer-2 target, and this sprint is what moves it
back to reviewable.
