# The balanced λ column is computed at LiH's bond length — for every molecule

**Found 2026-08-31**, while verifying that the Stage-4 ledger's "balanced QWC
column" open item was genuinely closed. It was — but the verification walked
into a larger defect one column over.

**Status.** *Code FIXED, papers NOT YET UPDATED.*

PI direction 2026-08-31: choose the option that reinforces the mission
(charting the forced/free seam). That is **per-molecule experimental geometry**
— a λ projected through a fictitious shared R describes no physical system, so
it can carry no falsifier: it is neither forced nor a *named* free-side
projection, which is the one thing the three-verdict rule forbids.

Done:
- `build_balanced_hamiltonian` now takes `R=None` and resolves the geometry
  from the spec; the old 3.015 default is gone and the fallback path raises
  rather than substituting silently.
- The corrected λ column is measured for all 12 second/third-row molecules
  (§2b below).
- The seam itself is pinned by
  `tests/test_paper20_geometry_independence.py` (see the corrected blast-radius
  table below).

**LANDED 2026-08-31.** Both tables now carry the corrected columns and an
explicit `R (bohr)` column, so the free-side input is named rather than
inherited; both captions state the convention and disclose the retired cells.
`tests/test_paper20_balanced_lambda.py` is re-synced (all twelve rows at their
own geometry, QWC added to each pinned row, and the two cheapest rows promoted
out of `--slow`). The registry was **inverted**: it had been populated from the
2026-08-30 pass, so it held the LiH-geometry values as canonical and the
correct ones as retired — C21 flagged the corrected papers on the first run
after the fix. The other ten rows of the column, which had no registry keys at
all, are now registered.

Verification: C21 PASS on every target except group2 (its 24 are a separate,
already-scoped sprint); Papers 14 and 20 compile with zero undefined
references; escape gate clean; gate-coverage matrix PASS; 17/17 registry tests;
9 passed / 11 slow-skipped across the two geometry test files, and 13/13 with
`--slow`.

---

## 1. The mechanism

`build_balanced_hamiltonian(spec, R=3.015, ...)` — `R` defaults to **3.015
bohr, LiH's bond length**, and the function does **not** read the geometry from
the spec it was handed. `spec.R` is consulted only to *correct V_NN* from the
spec's own R to the requested R (the 2026-06-07 W1e-Projection-Audit fix), so
the nuclei end up placed at the passed R regardless of which molecule it is.

Called as `build_balanced_hamiltonian(nah_spec())`, it therefore computes
**NaH at LiH's bond length**. The docstring is honest about the default —
"Used for LiH backward compatibility"; "If None, defaults to LiH geometry" —
but nothing warns a caller that the default silently applies to every other
molecule too.

## 2. The measurement

Every λ printed in Paper 14 `tab:second_row` (and the matching column of
Paper 20 `tab:molecules`) reproduces at **R = 3.015**, and **none** reproduces
at the molecule's own bond length:

| molecule | its own R | paper λ_ni | at R = 3.015 | at its own R | matches |
|:---------|----------:|-----------:|-------------:|-------------:|:--------|
| NaH  | 3.566 |  20.6 |  20.6073 |  19.5821 | R = 3.015 |
| MgH₂ | 3.261 | 110.5 | 110.5274 | 111.8003 | R = 3.015 |
| HCl  | 2.409 | 869.7 | 869.7200 | 866.1325 | R = 3.015 |
| H₂S  | 2.534 | 879.6 | 879.6158 | 873.6323 | R = 3.015 |
| PH₃  | 2.680 | 895.3 | 895.2688 | 889.3493 | R = 3.015 |
| SiH₄ | 2.800 | 914.1 | 914.1448 | 909.1540 | R = 3.015 |

No molecule's own R is 3.015. The error is 0.6%–5.2%, systematic, and in both
directions (NaH and HCl are over-stated, MgH₂ under-stated).

**Blast radius, corrected 2026-08-31 after widening the sample.** An earlier
version of this memo said "the Pauli *and QWC* columns are geometry-independent,
so the blast radius is exactly one column." That was measured on NaH and HCl
only, and it is **wrong for QWC**. Widening to the rest of the library:

| column | behaviour under a geometry change | seam verdict |
|:-------|:----------------------------------|:-------------|
| N_Pauli | **bit-identical** at every molecule tested | **FORCED** — fixed by the angular selection rules |
| N_QWC | **moves**: H2S 1264↔1282, PH3 1427↔1477, SiH4 1566↔1551, MgH₂ 903↔877 | **heuristic** — greedy grouping, algorithm-path dependent |
| λ | moves, up to 12.7% (KH) | **FREE** — projection through the bond length |

So the blast radius is **two** columns, not one, and the two move for different
reasons. NaH / LiH / HCl happen to be QWC-insensitive, which is why three
molecules looked like an invariance. Three agreeing points are not an
invariance — the same lesson as the arXiv magnitude claim earlier the same day.

The three-way split is the mission-relevant result: one table, three columns,
three different places on the forced/free seam. Pinned in
`tests/test_paper20_geometry_independence.py`, including the QWC leg asserted
in the *negative* so the false invariance cannot be re-derived.

## 2b. The corrected column (measured at each molecule's own R)

`python debug/qa/remeasure_balanced_lambda.py`

| molecule | R (bohr) | published λ (at 3.015) | corrected λ | Δ% |
|:---------|---------:|-----------------------:|------------:|-----:|
| NaH  | 3.566 |  20.6 |  19.6 | −4.94 |
| MgH₂ | 3.261 | 110.5 | 111.8 | +1.18 |
| HCl  | 2.409 | 869.7 | 866.1 | −0.41 |
| H₂S  | 2.534 | 879.6 | 873.6 | −0.68 |
| PH₃  | 2.680 | 895.3 | 889.3 | −0.66 |
| SiH₄ | 2.800 | 914.1 | 909.2 | −0.54 |
| KH   | 4.243 |  32.2 |  28.1 | **−12.66** |
| CaH₂ | 3.807 | 128.0 | 123.8 | −3.30 |
| HBr  | 2.670 | 877.0 | 875.3 | −0.19 |
| H₂Se | 2.760 | 883.4 | 883.4 | 0.00 |
| AsH₃ | 2.820 | 885.2 | 885.2 | +0.01 |
| GeH₄ | 2.870 | 877.3 | 874.8 | −0.28 |

KH is the worst at −12.7%, unsurprisingly: its bond length (4.243) is furthest
from the 3.015 that was substituted. The two rows that barely move (H₂Se,
AsH₃) do so by coincidence, not by insensitivity — their own R happens to sit
where λ is locally flat.

## 3. Why nothing caught it — the guard existed and was invisible

`tests/test_paper20_balanced_lambda.py` was written 2026-07-02 for exactly this
column, and it holds the **physically correct** values: it passes NaH's real
R = 3.566 and HCl's real R = 2.409 explicitly, pinning λ = 19.6 and 866.1.

Then:

1. The test is marked `@pytest.mark.slow`, so **12 of its 13 cases are skipped
   by default** (`pytest tests/test_paper20_balanced_lambda.py` → "1 passed,
   12 skipped").
2. The 2026-08-30 re-measurement pass called the builder **without** `R`, got
   the LiH-geometry values, and wrote them into both papers — recording it in
   the captions as *closing a stale cell*: "NaH $19.6 \to 20.6$ and HCl
   $866.1 \to 869.7$".
3. The test still passes, because it uses its own explicit R.

So the paper and its backing test now disagree about what the column means, and
**both are internally consistent**, so no gate fires. The correction moved the
paper *away* from its test, not toward it.

The four rows without an explicit R in the test (MgH₂, H₂S, PH₃, SiH₄) inherit
the same default, so the test agrees with the paper there — at the wrong
geometry, in both places.

## 4. What the PI has to decide

Not "is this wrong" — the geometry is wrong either way — but **what the column
is for**:

- **(a) Per-molecule geometry.** Each row at its own equilibrium bond length.
  Physically meaningful; changes all six published λ values; makes the existing
  test's NaH/HCl values canonical and requires the other four rows re-measured
  and the test's implicit-default rows given explicit R.
- **(b) Fixed comparison geometry.** All rows at one R deliberately, so the
  column isolates electronic-structure differences from bond-length
  differences. Defensible for a *resource* table — but then it must say so in
  the caption, the R must be chosen on purpose rather than inherited from a
  backward-compatibility default, and 3.015 (LiH's) is an odd choice for a
  second-row table.

Either way the builder's signature should stop defaulting to a specific
molecule's geometry — take `R=None` and either read `spec.R` or raise.

## 5. Recommended, whichever branch is chosen

1. Make `R=None` the default and resolve from `spec.R`, so no caller can
   silently compute one molecule at another's geometry.
2. Drop `@pytest.mark.slow` from the cheap rows — NaH builds in 0.3 s and MgH₂
   in 4 s. A guard that never runs is not a guard; these two would have caught
   the drift the day it landed.
3. Re-measure the column under the chosen convention and state the convention
   in both captions.
4. Register the λ values so C21 guards them (they are currently unregistered
   multi-document numerals — this column is precisely the "twin across two
   papers" class the registry exists for).

## 6. Reproduce

```
python debug/qa/verify_balanced_qwc.py          # all six rows, default R
python -m pytest tests/test_paper20_balanced_lambda.py --slow -q
```

`debug/qa/verify_balanced_qwc.py` also documents the identity-term convention
trap in the returned dict (`N_pauli_composed` excludes the identity,
`N_pauli` includes it, `one_norm` includes it while the papers print λ
without it) — three different conventions in one return value.
