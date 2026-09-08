# /qa 58/59/60 FULL + group1 DELTA — findings and remediation

**Invoked:** PI, 2026-09-07 (`/qa full each of the papers we edited recently`);
scope selected by the PI as **58/59/60 FULL + group1 DELTA**.
**Mode:** unseeded (default since 2026-09-02). Standing calibration record cited
per `.claude/commands/qa.md`: trunk FULL #2 19/21 seeds 0/8 FP; DELTA #1 9/9,
0/8 FP; DELTA #2 8/9, 0/8 FP.

---

## Verdict

**FAIL, remediated.** One LARGE (Paper 60, at two loci), plus nine instrument
defects and a set of SMALLs. Every finding below was verified by the PM against
primary text or code before acceptance, and every instrument fix was proven to
discriminate in **both** directions before being recorded green.

The single most useful thing this run produced is not a paper fix — it is that
**five gates were reporting PASS while examining nothing**, and one gate was
mis-scoped in a way that made three of my own earlier "12/12 PASS" reports
meaningless. Those are recorded first.

---

## 1. The LARGE: Paper 60's `eq:sublinear`

### What the paper said

> `‖M‖₁ ∼ K^0.84` … "grows **sublinearly**" … "The mechanism is the one Avery
> notes for Coulomb potentials: the **pure-number matrix elements** decrease
> with quantum number …, so the sum grows more slowly than the matrix
> dimension."

### What the measurement says

Re-measured with the paper's own `gen_configs` / `solve`, pushed past its
largest fitted point (K = 164) to K = 340:

| n | K | ‖M‖₁ | local slope | diag share |
|---|---|---|---|---|
| 7 | 74 | 85.359 | — | 64.6 % |
| 8 | 100 | 109.869 | 0.838 | 62.4 % |
| 9 | 130 | 136.975 | 0.840 | 60.5 % |
| 10 | **164** | 166.539 | 0.841 | 58.8 % | ← the paper's last point
| 11 | 202 | 198.796 | 0.850 | 57.2 % |
| 12 | 244 | 234.237 | 0.868 | 55.5 % |
| 13 | 290 | 272.775 | 0.882 | 53.9 % |
| 14 | 340 | 315.069 | **0.906** | 52.2 % |

```
global fit, window K ≤ 164   : K^0.8399   ← reproduces the published 0.84 exactly
global fit, extended K ≤ 340 : K^0.8541
nuclear diagonal T⁰ = Z·ΣR_ν : K^0.704    (local slope FALLING, 0.708 → 0.697)
off-diagonal T′              : K^1.049    (SUPERlinear)
```

**Two defects, not one.**

1. **The window is unstated.** 0.84 is a fit over K = 74–164 and the local slope
   rises monotonically outside it. "Grows sublinearly" as a regime is not
   established; the exponent trends toward 1.
2. **The mechanism is backwards.** The pure-number block `T′` — the part the
   paper credits, and the part that is attractive for a quantum device because
   it is generated on-chip from quantum-number labels — is the **superlinear**
   one. The sublinearity is carried entirely by the nuclear diagonal
   `T⁰ = Z R_ν`, whose own slope is flat to falling. The total exponent drifts
   up precisely because `T′`'s share of the norm grows.

**The paper already contained the correct statement.** §7 (molecular) says the
advantage "rides on the *clean diagonal* `T⁰ = Z R_ν`" and "does not survive the
loss of the diagonal `T⁰`". So §4 contradicted §7 inside one document, and §7
was right. The fix aligns §4 to §7 rather than inventing anything.

**What survives, and it is most of the claim:** ‖M‖₁ does grow more slowly than
the matrix dimension over every basis that can be computed, and the contrast
with the L² Löwdin superlinear inflation is real. The encoding claim stands;
its *asymptotic reading* and its *mechanism* were wrong.

### Remediation

- Paper 60: abstract, intro, `eq:sublinear`, the mechanism paragraph, §5 echo
  and the conclusion restated; new `eq:sublinear_split` records the T⁰/T′ split
  with its measurement.
- `papers/synthesis/group2_quantum_chemistry_synthesis.tex`: the Paper-60 block
  carried the bare "sublinear 1-norm" with none of the qualification (this was
  the second locus of the LARGE). Rewritten. Two SMALLs closed in the same
  block: "the one genuine lever" → the paper names **three** (gerade sector,
  large separation, one-electron-only metric), and the isoenergetic posing is
  **Avery's** and was credited to nobody.
- Numeric registry: four entries (`p60_onenorm_exponent`,
  `p60_onenorm_exponent_sonly`, `p60_T0_exponent`, `p60_Tprime_exponent`), so
  the **window and the mechanism are part of the registered object** and C21
  rule 2 forces a prose re-read if any of them moves. Paper 60 annotated with
  five `\gvq` citations; no-op verified by strip-and-diff.
- C16 entry `p60-sublinear-as-regime` with `cited_by`, proven to fire on the
  five retired phrasings and stay silent on the five corrected ones.
- `docs/claim_test_matrix.md`: row re-tiered; **`eq:sublinear_split` row added**;
  three separate uses of "asymptote" for the windowed 0.84 corrected.
- Backing test `test_paper60_sublinearity_is_carried_by_the_nuclear_diagonal`,
  written as a separate activity from the fix and **fire-tested in both
  directions** (perturbing the diagonal's K-scaling → FIRED; making T′ sublinear
  → FIRED).

---

## 2. Instrument defects — five gates were examining nothing

These matter more than any single paper finding, because each one made a PASS
that carried no information.

| # | Gate | Defect | Status |
|---|---|---|---|
| G1 | **C16** | Reported PASS while printing "this run checks NOTHING" for paper_59/60; my `grep '^RESULT'` filtered the warning out. | Fixed — `_grounded` self-audit ported from C17; errors out when no selected entry declares a locus in scope |
| G2 | **C16** | `cited_by` inert; two L5 entries missing most of their loci; the braced `\mathrm{height}_{P}` spelling defeated the pattern | Fixed + widened |
| G3 | **C9** | paper_59/60 done-records recorded C9 as "N/A"/"OUT OF SCOPE" on a stale premise | Rewritten as GATING |
| G4 | **C21** | Examined **0 annotations** on papers 59, 60, 61 | **Closed** — 8 / 1 / 2 annotations now checked. The salience report could not help: it lists only *multi-document* numerals, and its four "measurement-shaped" candidates on Paper 59 were bibliography volume/page numbers. The real quantities had to be read out of the papers |
| G5 | **C17** | 2 of 7 Paper-60 families grounded; **0 families** on paper_61 (the `_grounded` self-audit correctly *errored* rather than printing a hollow PASS) | **Closed** — family `p61-t2-superseded-19th-digit` added, guarding the superseded 19th digit of the pre-certification T2 anchor; discrimination-proven both ways. 4 / 1 families now grounded on 60 / 61 |
| D1 | done-records | 58/59/60 used `--gate group2`, so the gates ran against the wrong file set | Fixed → `--gate paper_NN` |
| D3 | **C11** | 27.6 % blind spot — keyed bibitem prefixes and `\emph{Paper N: Title}` forms unmatched | Fixed; **16 title defects** then found and fixed corpus-wide |
| D4 | **fire_test** | Reported a *skipped* test as "DID NOT FIRE / FAIL" | Fixed — `--slow`, plus NOT RUN → INCONCLUSIVE |
| M3 | **C22** | Check D matched only `import debug.x`; the `sys.path.insert(REPO/"debug")` + bare-import idiom was invisible — and that idiom was the corpus's **only** real instance | Fixed (`DEBUG_SYSPATH`), discrimination-tested, baselined |

Also: **Paper 61 was an orphan** — covered by no `qa_scopes` entry, so no gate
had ever examined it. Added, with a standing assertion in `qa_scopes.selftest()`
that computes `set(files) - covered` so a future orphan fails the selftest
instead of going unnoticed.

**My own error, recorded:** I reported "12/12 gates PASS across all seven
scopes" for scopes that *excluded the papers under certification*. Caught when
C10's group2 file list turned out not to contain 58/59/60. Corrected and re-run
against the correct single-paper scopes.

---

## 3. Paper 58 — code findings

- **MATERIAL-1 (blind decider).** `test_paper58_aabb_decided_census_is_195_of_195`
  asserts `(n_nonzero, n_gaunt, n_accidental) == (195, 0, 0)`. **Fire-tested:**
  hard-wiring `_decide_zero → return False` (the decider never decides a zero)
  leaves the test **green**, because 195/0/0 is exactly what a dead decider
  produces. The headline's content is that the decider *looked*; the tuple alone
  cannot distinguish looking from being asleep.
  *Fixed* by `test_paper58_census_deciders_are_alive`, pinning the decider in
  both directions on a real census expression (fire-tested: dead → FIRED,
  always-true → FIRED).
  *Measured while writing it:* **no quartet of `_ORBS` is Gaunt-forbidden at all
  (0 of 625)**, so `n_gaunt == 0` is **forced by the basis, not discovered**.
  Recorded as its own test rather than left to read as a measurement.
- **MATERIAL-3 (C22 blindness).** Above. The dependency is baselined rather than
  removed, on the grounds that the same H₂ number has permanent backing in
  `test_paper58_qfd.py` from `geovac/qfd_assemble.py`. Porting
  `step1_native_molecule` into `geovac/` is the standing repair.
- **MATERIAL-2 (dps-dependent QFD guard).** Examined; the assertion is
  meaningful (both values are created at `dps+25` and compared at a global
  `dps = 50`), but the guard's strength rides on a **global mutable**
  `mp.mp.dps`, and the docstring's "84-digit certified value" is checked to 38
  digits against a 60-digit constant. **Owed**, not yet fixed.

---

## 4. group1 DELTA

- One live C16 locus at `paper_46:1124` (`\mathrm{height}_{P} = 0`). My earlier
  withdrawal note sat 8 lines away and C16's window is ±5, so the gate was
  right to keep calling it live. The proposition's caption now leads with
  **REFUTED** inside the window.
- **A fix of mine broke the build:** the same scope-note edit had silently
  consumed `\end{proposition}`, so `paper_46` did not compile. Restored.
  (Second time this pass — the Paper 61 bibitem fix used macros undefined in
  that file, making C11 green while C10 failed.)
- L5-withdrawal remediation reached its citers: Papers 40, 46, 52 and the group1
  synthesis.
- **Owed:** the remaining group1 SMALLs (S1–S4, S6–S8), and **Paper 53's
  height-leg**, which appears to rest on the same withdrawn L5 assembly —
  whether that is a descope or a genuine repair is a **PI adjudication**.

---

## 5. The pattern this run repeats

Every substantive *mathematical* claim examined here survived independent
re-derivation. What failed, again, was the *apparatus*: five gates examining
nothing, one guard that could not fail, three of my own PASS reports that were
scoped away from the thing being certified, and two of my own fixes that broke
what they touched. This is the same asymmetry CLAUDE.md §9 already records from
the v5.4.4..v5.7.3 arc, and it is the argument for the guard-writing rule and
the gate self-audit rule rather than for more review effort.

One new instance worth naming: **C21 failed on me within a minute of my making
it non-vacuous.** I annotated the extended-range `0.854` with the windowed key
(0.84), and the gate I had just un-vacuumed caught it. That is what a working
gate looks like.


---

## 6. Follow-on pass, same day

### Paper 61: a definition of done, and its group

Paper 61 was split out of Paper 59 on 2026-09-06 and belonged to **no scope at
all** until this arc. Giving it a single-paper scope satisfied the corpus-wide
orphan check — while `/qa group3` still walked straight past it.

It is a **group3** paper: `papers/group3_foundations/`, the periods / Tannakian
arc, siblings 55–57. Now in the `group3` scope and in `docs/qa/group3.done.md`.
`docs/qa/paper_61.done.md` written: scope, **gating** C9, dimensions, ten
enumerated C8 headlines with tiers and backing, seven ranked watch-notes.

**New group-membership assertion** in `qa_scopes.selftest()` — a paper in
`papers/groupN_*/` must be in the `groupN` scope unless declared in
`GROUP_SCOPE_EXEMPT`. **Fire-tested:** removing 61 again makes the selftest FAIL
*while the orphan check still reports clean*. That is precisely why the narrower
assertion was needed, and it is the general lesson: **"is it in some scope?" and
"is it in its own group's scope?" are different questions**, and only the second
one matches how runs are actually dispatched.

### Paper 53 — the PI-flagged height leg

**It is not Paper 38's L5 transported.** L5's height was a
reconstruction-*defect* quantity, refuted by a finite-band `f` with `B(f) = 0`.
That same `f` gives `‖∇Bf‖ = 0` here, which **satisfies** non-expansiveness.
Different quantity, different failure mode — the witness does not transfer.

**But the leg has its own defect.** Paper 53 claimed the plane Berezin is
gradient-non-expansive and listed that among the *established* load-bearing
ingredients. `B` is a radial convolution, so `∇(Bf) = B(∇f)` and the sharp
constant is the kernel's L¹ norm:

| s | 0.6 | 0.75 | 1 | 1.5 | 2 | 3 | 5 | 8 |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| ‖K‖₁ | 5.72 | 3.23 | **2.01** | 1.42 | 1.23 | 1.09 | 1.02 | 1.004 |

`> 1` at every finite order — the kernel `J_{s+1}(r)/r^{s+1}` oscillates, so it
is negative somewhere — approaching 1 only as `s → ∞`. And it is
**Λ-independent** by scaling, so the reported "ratio rising to 1 as Λ→∞" cannot
have been the operator's gradient gain at all.

Corrected (`rem:height_constant`), backed by
`tests/test_paper53_height_constant.py` (3 legs, fire-tested both ways).
**[OPEN — PI]** whether a constant height preserves `Λ_prop ≤ Cγ_Λ → 0`. I did
not decide it: that is a question about Latrémolière's bookkeeping, not a
numerical one.

### group1 residue (carryforward U.2) closed

P40's cross-manifold future work aimed at Paper 39's **withdrawn** Pythagorean
route → redirected to the lifted-state route. `paper_43_..._outline.md` still
asserted the descoped "literal identification" at **both** loci (U.2 listed
one). P51 was never affected — it is group5, mis-filed in the U.2 enumeration.

**The mirror case:** the synthesis presented the proved `k=2` case as carrying
the `√k` triangle constant, while Paper 39 says that constant "belongs to the
*abandoned* (B,P)-pair route and is not a rate constant of this theorem". Owner
strengthened, summary kept the weaker form — the direction C16 and `cited_by`
structurally cannot catch. Registered for this occurrence; the class stays open.

### Corrections to my own work in this pass

Recorded because they are all one failure mode — **searching by filename when
the key is something else**:

- I said Paper 61 had **no** claim-matrix rows. It has **ten**; I grepped the
  filename, not the row key. (The same wrong-key mistake that made C11 unable
  to fail.)
- I logged a coverage gap for the disc−8 CM leg that
  `test_routeC_momentum.py::test_cosmic_galois_cm_periods_are_gamma_values`
  already backs to 1e-50 — again, I searched only `test_paper59_*`.
- The 59→61 split itself made this mistake: nine loci in two files kept the old
  attribution because the sweep keyed on `test_paper59_*` and
  `test_routeC_momentum.py` is not named that.
- My first `k=2` C16 pattern missed `$k=2$`; its replacement then
  **false-positived**, because C16 also scans a newline-joined copy where
  `[^
]{0,60}` stops bounding anything. Use `[^.
]`.
- My Paper 53 fix was locus-by-locus; the C16 entry I had written minutes
  earlier caught the two loci I left.

Two gates caught their own author within minutes of being written or widened
(C21 on the annotation I mis-keyed, C16 on the Paper 53 loci I missed). That is
the strongest evidence in this run that the instruments are working.

### A gate that read a non-certification as a certification

`check_cert_staleness.py` took the **max date anywhere in a record** as its
certification date. Consequences, all found by writing records that state their
status honestly:

- "STATUS: RE-RUN 2026-09-07 — NOT re-certified" made 58/59/60 read as
  **`current`** — the healthiest verdict in the table.
- "STATUS: NEVER CERTIFIED" made Paper 61 read as `current` too; the decline
  regex knew "NOT certified" and not "NEVER CERTIFIED" — the third
  spelling-defeats-the-pattern miss of the day.
- Three group dates were **overstated**: group3 2026-08-29 (real: 08-24),
  group4 08-30 (08-24), group6 08-30 (08-24), each picked up from
  "Post-certification touch" notes rather than a certifying line.
- group6's STATUS says **NOT CERTIFIED — superseded**, and the table said OWED.

Fixed in three steps, two of which were wrong: same-line-CERTIFIED-only lost
trunk's date outright (its status line is "FROZEN — certified PASS" with the
date elsewhere); date-comparison alone still mislabelled group6. The landing
rule is a hybrid — prefer the newest **CERTIFIED-asserting** dated line, fall
back to max-date only when a record has none *and report that fallback*, and
let the record's own **STATUS line** decide declines. Verified row by row
against all twelve records. `--selftest` added (there was none) pinning each of
the four failures, and fire-tested: reverting to max-date reproduces
`max_date_overstates: cert_date '2026-08-29', want '2026-08-24'`.

**The general lesson:** this is the audit built to catch stale certifications,
and it could be fooled by a date. An instrument one level up from the thing it
checks is not automatically more trustworthy than what it checks.
