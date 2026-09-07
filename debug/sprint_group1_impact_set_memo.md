# Sprint memo — the group1 impact-set cluster

**Date:** 2026-09-07 · **Version:** v5.10.8 · **Owed since:** v5.10.1 (2026-09-05)
**Canonical memo for this sprint** (§13.11 rule 1).

---

## 1. What this was

The v5.10.1 claim-impact reviewer found, by *reading the impact set rather than
the diff*, that two corrections had never reached their citers:

- the **L5 withdrawal** (2026-09-03, FULL #4 — Paper 38's height leg refuted), and
- the **Lorentzian descope** (2026-06 — Paper 45's K⁺ annihilation theorem).

It was logged as "almost all group1 (non-trunk) — a focused group1 remediation"
and deferred. This is that pass.

**Current-state check first (§9):** four versions had shipped since the finding
was logged, so the loci were re-verified rather than assumed still open. All were
still open.

---

## 2. The framing result: C16 said `clean`, twice, and was wrong both times

Running C16 on group1 returned **PASS** — which is exactly what v5.10.1
predicted, since a citer restates a claim in its own words. But two of the misses
were not paraphrase at all. They were **pattern-vs-typesetting**, the class
CLAUDE.md §9 already records (`$s/p$-lift` vs `s/p-lift`):

### Hole 1 — the braced subscript
`l5-height-bound-achieved` matched `height_P`, `height P`, `height}P`, but not
`\mathrm{height}_{P}`. The alternation was `height\}?[_ ]?P`: it allowed a brace
*before* the underscore but not after it. The corpus writes the braced form.

Consequence: the group1 synthesis carried
`contributes $\mathrm{height}_{P} = 0$` — **the refuted claim, stated as live, in
the document that summarises the paper it was refuted in** — and the gate called
the entry clean for four days. Fixed with `\{?`. Verified against all four
spellings; the gate named the locus on the next run.

### Hole 2 — a file list built from the wrong end
`hopf-base-label-for-4-over-pi` declared seven `files`. All seven were documents
where the label had been *fixed* in the v5.5.x pass. **None was a document that
merely cites the constant.** P42, P43 and P20 were therefore never scanned, and
all three carried the retired reading.

Its trailing alternative also required `is the` / `as the` before the label,
missing the parenthetical apposition — "asymptotic rate $4/\pi$ (the Hopf-base
measure …)" — which is how P20 and P42 actually write it.

*The generalisable lesson:* a retraction entry's `files` list should be built
from the **citation graph of the claim**, not from the edit history of the fix.
The fix touches owners; the claim travels to citers. This is the same asymmetry
`cited_by` exists to capture, appearing one level down in `files`.

**Gate-first order was followed** (§9): register/repair the pattern → prove it
discriminates → let the gate enumerate → fix every named locus → re-run until
dry. Never fix-then-register.

**Score:** 5 loci found by the gate *after* the two repairs; 4 more found only by
reading. Patterns and reading each caught what the other could not.

---

## 3. Content remediated (9 loci, 7 documents)

| Document | Defect | Severity |
|:--|:--|:--|
| group1 synthesis §L5 | withdrawn (B,P) assembly presented as "completing the proof"; false `height_P = 0` | **LARGE** |
| P47 proof | L5's reach/**height** bookkeeping asserted to transport verbatim — inside a proof | **LARGE** |
| P42 abstract | continuum Lorentzian propinquity "the named open frontier" (superseded June 2026) | MATERIAL |
| P43 abstract | same | MATERIAL |
| P44 §recommendation | recommends the compactification that is now known to be the obstruction | MATERIAL |
| P40 intro ×3, P42 crossref, P44 list, P47 header, field guide | five-lemma echoes presenting L5 as live | SMALL |
| P42 ×2, P20 ×1 | 4/π attributed to a Hopf base-to-total ratio / SU(2)/U(1) Haar quotient | MATERIAL |

### The synthesis rewrite

Its §L5 said L5 "completes the proof" while §tensor_extension, two sections
later, called that same (B,P) route *abandoned* — the document contradicted
itself. Replaced with **"What carries the theorem, and the withdrawn Lemma L5"**:

1. the live argument — UCP compression S(f) = P M_f P, dual map υ, four facts,
   S∘υ = the conjugation average Φ against the Fejér kernel; **no Berezin map, no
   partial inverse, hence no reach and no height** — which is *why* the
   refutation leaves it untouched;
2. L5's refutation — finite-band witness gives an f with B(f) = 0 whose height
   quantity is exactly 1, so height_B ≡ 1 and the bound fails for every n_max ≥ 6;
3. the net effect — **one proof path fell and the keystone stands; WH1 unaffected.**

### P47 — flagged, not repaired

Its inner-arrow proof claimed L5's height bookkeeping transports verbatim. The
honest repair is *not* a change of constants: the cell-wise argument has to be
rebuilt on the lifted-state route, and that has not been done. Recorded as a
**second, independent** reason to read `thm:inner` as descoped, alongside the
degenerate-seminorm reason already in its statement. A proof not re-derived is
not a proof to assert.

### The Lorentzian statement, sharpened

P42/P43/P44 now carry the mechanism rather than a bare "descoped": these papers'
own result (i) — an **integer** modular spectrum — is what makes the
finite-cutoff modular flow a compact β = 2π circle rather than a non-compact
boost, hence signature-blind seminorms. P44 additionally records the general
form:

> Signature is a non-compact phenomenon. Quantum-metric convergence is a
> finite-dimensional-approximant phenomenon. Finite-dimensionality forces
> discrete spectrum, so *every* finite truncation compactifies the modular flow.

This also explains P43's own "signature-independent residual" finding, which had
been reported as a curiosity.

### What was deliberately *not* changed

"Hopf-base measure" remains the legitimate historical name for the k = 0 Mellin
slot and for Paper 25's Vol(S²)/4 = π — 163 uses corpus-wide, declared historical
in Paper 18. Only reading **4/π itself** as a base-to-total ratio was retired. The
gate entry's own note already drew this line correctly; the fix respects it.

---

## 4. Verification

- C16 + C17 PASS on trunk, group1, group2, group3, group4; `--dependencies` PASS.
- 8 further deterministic gates PASS.
- Six papers + synthesis + field guide: 3-pass compile, **0 errors, 0 undefined**.
- 61 tests pass, 1 skipped — including the WH1 keystone
  (`test_p38_action_seminorm`, `test_p38_metric_convention`) and the 18 symbolic
  S³ proofs.
- Three macro slips of my own (`\KS`, `\Kplus`, `\SU` in documents that do not
  define them) were caught by compiling and fixed.
- One missing `\bibitem{paper45}` in P42, introduced by my own new citation,
  caught by the compile and added.

**Pre-existing, not mine, not fixed:** Paper 20 carries 51 undefined references.
Verified against a stashed baseline — identical count without this pass's edits.
Out of scope (group4). Logged as owed.

---

## 5. Owed

1. **P47's inner-arrow rebuild** on the lifted-state route — the only way to
   restore that item; currently flagged open.
2. **Paper 20's 51 undefined references** (pre-existing, group4).
3. A **`files`-list audit for every C16 entry**, on the lesson above: entries
   whose file list was built from the fix rather than from the citation graph
   will have the same blind spot. This pass fixed one entry; the class is
   unaudited.
4. Paper 61 cert + a clean re-delta remain the standing precondition for the next
   certifying FULL run (carried from v5.10.6).
