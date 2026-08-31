# Stage-4 cross-document ledger (2026-08-30)

**Verdict: DEFECTS found and remediated.** This is the cross-document pass
that the delta-1..7 cycle kept deferring — the class of defect that lives
*between* two documents and is therefore invisible to a reviewer scoped to
either one.

The pass was run against the ledger accumulated across delta-4..7 plus the
newly-built C21 numeric gate. Every value below is **measured**, not
inferred from a sibling cell; where a claim could not be re-derived it was
withdrawn rather than re-priced.

---

## 1. What the ledger held, and what happened to each item

| id | claim | verdict |
|:--|:--|:--|
| M3 | P20 keeps a scalar-vs-rel λ comparator built from a removed column and retired values | **WITHDRAWN** |
| M4 | `tab:sunaga` projects "17–32×" from two retired inputs, in *both* papers | **WITHDRAWN** (column struck) |
| M6 | P14 Conclusion prints 6.8× where 7.1× is measured; the discussion sentence self-contradicts | **FIXED** (3 loci) |
| M7 | P20 prose "CaH sits at 4.74×" vs its own table's 4.783× | **FIXED** |
| M8 | P14 "NaH and KH both yield 223/239" vs table 559/575 | **FIXED + measured** |
| M9 | P26 quotes the pair-diagonal 2.76% as "the angular ERI density" | **FIXED — and it was a 4-document cluster** |
| M10 | my own garbled fit-range sentence at P14:448 | **FIXED** |
| M11 | P20 asserts a basis-axis inference it retracts 50 lines later | **FIXED** |
| NIT 1 | 54× vs 51× equal-qubit low end across 3 docs | **already consistent** (no live 51×) |
| NIT 2 | balanced 1-norm 1,511 vs 1,440 | **FIXED — both were wrong** (1,612.0) |
| NIT 3 | LiH composed λ 27.7 vs 34.0, undeclared | **FIXED — a whole-column convention mix** |
| NIT 4 | 2,434,441 identity convention in `tab:pauli` footnote | **FIXED** |
| NIT 7 | dangling `sec:block_pauli_count` | **FALSE POSITIVE** — label exists, C10 shows 0 undefined refs |

## 2. The two findings worth more than the ledger that produced them

### (a) M9 was not a number, it was a four-document zombie

P26 quoting `2.76%` looked like a one-cell citation slip. Tracing it found
the *cause*: until the 2026-08-29 wrong-sign-`q` fix, the production
enumerator `angular_zero_count` silently returned the **pair-diagonal**
density `D_pd` while enforcing the global-`M_L` rule textually — the buggy
coefficient zeroed every `m`-changing multipole, so the global constraint
was redundant. The corpus wrote that artifact down as a *fact about the
pipeline*:

- `P22 tab:sparsity` — printed 97.24 / 0.00 / 2.76 for five potentials.
  Re-running the **same** production routine today gives 91.48 / 0.00 / 8.52.
  The table was a fossil of the bug.
- `P22` Theorem-3 note + `group3` synthesis caption — "the density the
  **composed pipeline realizes**". The fix made that false.
- `P26` — quoted `D_pd` as *the* density, overstating symmetry-enforced
  sparsity by **3.1×**.
- `tests/test_paper22_density.py` — a `NOTE (FLAG, do not "fix" here)`
  pointing at a test name that no longer exists and calling CF-1 open.

A "do not fix" flag aimed at a deleted test is worse than no flag: it tells
the next reader the discrepancy is known and deliberate.

Regenerating the table also surfaced something the fossil concealed —
**Yukawa binds only 17 of the 18 states** at the tested screening
parameters, so its row was never a matched-orbital-count comparison. The
invariance test already groups by orbital count, so this is a presentation
gap, not a new defect; it now carries a footnote instead of a silent `18`.

Both retired phrasings are now in the **C16 registry** with two-way
discrimination proven against the actual pre-fix text recovered from git
(fires on all three retired loci, silent on all three corrected).

### (b) The 1-norm columns were mixing two correct conventions

Not one stale number — a *systematic* interleave, in two papers:

| | Pauli | λ |
|:--|:--|:--|
| `P14 tab:balanced` | composed non-identity, balanced **identity-included** | composed **identity-included**, balanced non-identity |
| `P20 tab:resources` | both identity-included | composed identity-included, balanced **non-identity** |

`P14 tab:balanced`'s caption *declared* "Pauli counts exclude the identity
term" while its balanced column did the opposite. Every individual number
was correct; only the pairing was wrong. This is precisely the class C21's
check E exists for, and it could not see these cells because they carried no
`\gvq` annotation — now they do.

Ground truth, measured 2026-08-30 from the live builders:

| system | comp all / ni | bal all / ni | λ comp all / ni | λ bal all / ni | comp+PK (all) |
|:--|--:|--:|--:|--:|--:|
| LiH | 838 / 837 | 2726 / 2725 | 34.036 / 27.715 | 78.156 / 75.231 | 38.552 |
| BeH₂ | 1396 / 1395 | 8868 / 8867 | 68.812 / 54.286 | 328.110 / 289.581 | 374.866 |
| H₂O | 1954 / 1953 | 19742 / 19741 | 372.464 / 186.761 | 1612.010 / 1439.751 | 28061.879 |

## 3. C21 earned its build again

`tab:composed_onenorm` carried a caption saying its λ column "awaits a
convention-matched re-measurement" — while the prose beneath spent those
cells as live facts, including a headline **78× PK-partitioning reduction**
and a `λ = 355 at Q = 50` that Paper 20 separately identifies as a
*deprecated legacy builder value*. C21 flagged three cells; reading the
region showed the whole table was pre-exact-rule.

A standing hedge nobody discharges is how a wrong number keeps getting
quoted. So it was measured, not hedged again:

| | n=1 | n=2 | n=3 / n=4 |
|:--|--:|--:|--:|
| H₂ λ | — | 8.307 | 46.498 / 163.585 |
| LiH λ | 15.347 | 38.552 | 225.235 |
| BeH₂ λ | 211.199 | 374.866 | 801.693 |
| H₂O λ | 18914.898 | 28061.879 | — |

Trotter columns were **regenerated, not re-derived**: the printed table
already obeyed `r = λt/√(2ε)` at every cell to the digit (LiH n=2:
37.33/√0.002 = 834.7 → 835), so the same map carries the new λ.

Consequences: H₂O partitioning **78× → 75.3×**; LiH PK share
**10.9% → 11.7%**; λ_elec **33/66/359 → 34.0/68.8/372.5**. Two claims
**survived** re-measurement unchanged and are worth naming, because a pass
that only ever cuts is miscalibrated: the H₂O PK share (98.7%) and the
`~2,387 Ha` PK diagonal — measured `2386.769`, exactly four times, on the
O-side valence blocks.

## 4. Gate work

- **C21 scope `group3` added.** The eri-density family is owned there, and
  `--gate group3` did not report "unscoped" — it **crashed** with a
  `TypeError` from `_files(None)`. An unknown gate now fails loudly, so a
  typo'd scope reads as a missing scope rather than a broken checker.
- **`density` kind added to `_KINDS`.** Angular ERI density carries a
  *selection-rule* convention (global-`M_L` vs pair-diagonal), not an
  identity convention — a genuinely new axis for check E.
- **`test_every_convention_string_parses` caught three of my own entries.**
  That guard exists because `convention` is free text: an unrecognised
  phrasing drops silently out of check E, so the entry *looks* guarded and
  is not. Worth recording that it fired on the person who wrote it.
- **C10 caught a cross-document `\ref`** I introduced (`tab:resources` lives
  in P20, not P14). It only catches this because it was hardened to scan for
  undefined references rather than trust `pdflatex`'s exit code.

## 5. Process note

The bash-heredoc backslash-halving rule fired **twice more** this pass, on
edits I judged small enough to shortcut. It is not a size-dependent hazard.

## 6. Green state

| gate | group3 | group4 | group6 |
|:--|:--|:--|:--|
| C10 compiles / 0 undefined refs | PASS (12 docs) | PASS (5) | PASS (5) |
| C11 internal titles | PASS | PASS | PASS |
| C13 test refs | PASS | PASS | PASS |
| C14 file refs | PASS | PASS | PASS |
| C15 inline arXiv | PASS | PASS | PASS |
| C16 retracted terms | PASS | PASS | PASS |
| C17 headline numbers | PASS | PASS | PASS |
| C18 duration language | PASS | PASS | PASS |
| C19 eaten escapes | PASS | PASS | PASS |
| C20 inline attributions | PASS | PASS | PASS |
| C21 numeric consistency | PASS | PASS (83 annotations) | PASS |

Tests: 103 passed / 14 skipped across the registry, headline-number,
spinor-ordering, scaling, P22-density, nuclear-sparsity and topological sets.

## 6b. Two more clusters, found by tightening the gate as a diagnostic

Running the discrimination harness on the newly-retired values showed four
of nine never firing *even on the pre-fix text* — the exact "an entry that
fires on nothing is worse than none" failure. Diagnosing it found two
things.

### `tab:tc_composed` was pre-exact-rule behind a blanket exemption

C21 exempted on a ±8-line window that did double duty: finding
require-context **and** accepting disclosure markers. One retirement note
five lines below `tab:tc_composed` therefore exempted the whole table, so
its retired counts `334 / 556 / 778` — `334` is named as retired in
CLAUDE.md itself — were reported clean. **A disclosure about one number is
not a disclosure about its neighbours.**

Re-measured (electronic-only, PK partitioned — the caption's own basis):

| | Q | std | TC | ratio | λ std → TC |
|:--|--:|--:|--:|--:|--:|
| LiH | 30 | 837 | 1,353 | 1.616 | 27.715 → 35.561 (1.283×) |
| BeH₂ | 50 | 1,395 | 2,255 | 1.616 | 54.286 → 70.169 (1.293×) |
| H₂O | 70 | 1,953 | 3,157 | 1.616 | 186.761 → 234.093 (1.253×) |

The **uniformity claim survives** — identical to three decimals — with the
constant moving `1.68× → 1.616×`. The 1-norm overheads do not: `9–16%`
with H₂O at exactly `1.00×` becomes a tight `25–29%` band, so that
coincidence and the PK-dominance explanation built on it are gone. The
vintage note's *scope caveat* (fixed-basis ⇒ says nothing about the
basis-axis exponent) was **kept**, not retired along with its numbers.

### The angular gradient's added support is entirely unphysical

Paper 14 quoted the angular-gradient multipliers (`2.66×`, `2.31×`,
`2.67×`, total `4.49×`) as plain measurements. The operator producing them
is defective: **every entry it adds violates L_z** — measured `78/78` for
composed LiH.

Verified rather than inherited from CLAUDE.md, and the convention control
was load-bearing: the composed ERI is stored in **chemist** order, so the
rule is `m_a + m_c == m_b + m_d`. The physicist reading reports **126
violations on the plain composed tensor** — a confidently wrong answer. The
control (zero violations on plain composed *and* radial-only TC) is what
makes the 78 mean anything.

The multipliers are **withdrawn, not re-measured** — re-measuring a broken
operator yields a number you must disclaim in the same sentence. The
paper's *verdict* (angular gradient = negative for quantum efficiency)
stands and is now over-determined. Backed at the level the claim is made:
`tests/test_tc_angular.py::test_composed_angular_adds_only_lz_violating_entries`
(the pre-existing companion pins the block level, 26 entries, different
index convention — a gap where a claim looked backed and was not).

### Gate change: the exemption window is now separate, and tighter

Swept 8/4/3/2/1 in both directions before changing anything:

| window | fires on pre-fix 66.0 | false positives on corrected corpus |
|--:|:--|--:|
| ±8 | no | 0 |
| ±4 | **yes** | 0 |
| ±3 | **yes** | 0 |
| ±2 | yes | 1 (`
ightarrow` disclosure) |
| ±1 | yes | 6 |

`EXEMPT_WINDOW = 3` is the widest setting that gains sensitivity at zero
cost. `WINDOW = 8` is **kept** for require-context, where a caption
legitimately sits many rows above its cell — that was why 8 was chosen and
it is still right for that job. `
ightarrow` joined `	o` in EXEMPT so
the tighter setting does not depend on which arrow macro an author picked.
Pinned by `test_exemption_window_is_tighter_than_the_context_window`, which
asserts both directions on synthetic text.

## 7. Still open

- **101 unregistered multi-document numerals** in group4, 24 in group3
  (C21 check D, advisory). Visible maintenance debt, not defects — the
  registry now names what is guarded, so what is *not* guarded is countable
  for the first time.
- The `tab:molecules` / `tab:second_row` **balanced QWC columns** print
  `---` under a "pending re-measurement" note. They are not pending: the
  builder returns `n_qwc_balanced` directly, and it reproduces
  `tab:resources`'s first-row values (629 / 1220 / 1494) exactly. Being
  measured for the full library.

## 8. Honest scope

**Closed at measurement grade** (re-run from the live builders, values in
this memo): the composed/balanced Pauli and 1-norm family for LiH/BeH2/H2O;
`tab:composed_onenorm` across H2/LiH/BeH2/H2O; the TC composed table;
`tab:sparsity` for five potentials; the balanced QWC column for all twelve
second/third-row molecules; NaH = KH and the other five isostructural
pairs; the H2O PK diagonal (2386.769 x4).

**Closed at proof grade**: the angular gradient adds only L_z-violating
support (78/78 composed LiH, with a zero-violation control on two clean
tensors); the C21 exemption-window split (two-way discrimination on
synthetic text); the two new C16 entries (two-way discrimination against
pre-fix text recovered from git).

**Withdrawn, not re-priced** -- no valid input survives to rebuild them:
the `tab:sunaga` matched-Q projection and its 17-32x claim (both papers);
Paper 20's scalar-vs-relativistic lambda comparator; Paper 14's
angular-gradient resource multipliers.

**Explicitly NOT established.** The TC/standard Pauli ratio is uniform
*across molecules at fixed basis* (n_max = 2, 1.616 to three decimals).
That is NOT evidence the TC overhead leaves the basis-axis exponent alone
-- the original vintage note flagged exactly this single-point-fallacy
risk, and re-measuring the numbers does not discharge the caveat.  It is
retained verbatim in the paper.

**Named open follow-ons:**

1. **142 unregistered multi-document numerals** corpus-wide (101 group4,
   24 group3, advisory).  Not defects -- but the registry now makes the
   unguarded surface countable, which it was not before.
2. **The angular-gradient assembly bug itself is unfixed.**  This pass
   disclosed and fenced it; nobody has re-derived the assembly.  Deliberately
   not filtered: the L_z violations are the only detector, and masking them
   would hide whatever the same mis-assignment does to entries that land in
   a *valid* M_L sector.
3. **C21 scopes cover group3/4/6 only.**  Every other group's numerals are
   unguarded by check C.  The M9 trace reached group3 only because a defect
   pointed there; there is no reason to think the other groups are cleaner.
4. **Repo health WARNs stand**: CLAUDE.md 218 KB (budget 150), debug/
   top-level 1368 files (budget 600).  Neither blocks; both are growing.
