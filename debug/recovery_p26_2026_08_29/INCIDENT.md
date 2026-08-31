# INCIDENT: Paper 26 source truncated and over-reverted (2026-08-29)

**Status: content preserved, source reconstruction OWED.**

## What happened

Two errors, compounding.

1. **Truncation.** A cleanup line in my edit script read
   `io.open(P, "w").write(io.open(P).read())`. Python evaluates the object
   expression before the arguments, so the open-for-write **truncated the
   file before the read** — writing back an empty string.
   `paper_26_entanglement.tex` went to 0 bytes.

2. **Over-reversion.** I recovered with
   `git checkout papers/.../paper_26_entanglement.tex`, which restores the
   last **committed** state. The working tree also held **uncommitted edits
   from earlier in the same session**, so the checkout discarded more than
   the truncation had. The correct first move would have been to look for a
   non-git copy (the built PDF) *before* touching git.

## What was lost from the working tree

Relative to the pre-truncation source, the restored HEAD version (709 lines)
is missing:

- **The §II.C rewrite (a PI-delegated LARGE from earlier this session).** HEAD
  still carries the **false `1/(8Z^2)` "matches" claim**; the corrected text
  (measured `Z^-0.85` decay of the off-diagonal energy fraction, plus the
  explicit negative and the measure characterisation) is gone. *This is the
  most important loss — an uncorrected false claim is now back in the source.*
- The ceiling-attainment tier language and the `1.837/1.785/1.644` hub values.
- All of 2026-08-29's ERI-arc edits: §III `42.4% -> 17.1%` and the counts,
  the basis-size trend reversal, the molecular cells, the ceiling retraction,
  the core-closure scoping.

## What survives (authoritative)

- `paper_26_PRE_TRUNCATION.pdf` — the built PDF from the pre-truncation
  source (8 pages, compiled 21:49, errors=0). **This is the complete record
  of the lost state**, missing only the last 3 edits of `p26_final.py`.
- `paper_26_PRE_TRUNCATION_text.txt` — `pdftotext` extraction of the above.
- The edit scripts (`fix_p26_*.py`, `p26_final.py`) — every 2026-08-29 change
  as exact old/new string pairs. Their *new* strings are the target text;
  their *old* anchors no longer match HEAD, which is precisely the diagnostic
  for what the earlier-session remediation had changed.

## Reconstruction plan

1. Re-author §II.C in HEAD from the PDF text (the `Z^-0.85` result and the
   explicit negative) — **do this first**, it removes a live false claim.
2. Re-apply the 2026-08-29 chain in order: `fix_p26_paper` → `_molecular` →
   `_ceiling` → `_icv_theorem` → `p26_final` (paper part), re-anchoring each
   to the reconstructed text.
3. Apply the post-order-fix final values (already measured, in `p26_final.py`
   docstring): hub `0.148/1.063/0.311`, sampled max `0.69/1.34/1.32`,
   degeneracies `4/7/6`, exact core closure at N/O/F.
4. Compile three-pass; run `tests/test_paper26_entanglement.py` (already
   pinned to the final values — **the tests are intact and are the spec**).
5. Re-run C10/C16/C17/C19 on group6.

## Note on the tests

`tests/test_paper26_entanglement.py` was **not** affected — all of today's
pins landed before the truncation and are correct. The suite therefore
encodes the target numbers and is the reconstruction's acceptance test.

## Process lessons

- **Never `open(path,"w")` on a file you are also reading in the same
  expression.** Read fully, close, then write.
- **On data loss, inventory non-git artifacts (built PDFs, caches, logs)
  BEFORE running any git restore** — `git checkout` is not a neutral undo
  when the working tree holds uncommitted work.
- The standing bit-exactness/verification discipline did its job here: the
  loss was detected within one command because the compile check reported
  `errors=2` where it had been clean.

---

## Update (same session): the loss is DEEPER than first scoped

Replaying the 2026-08-29 chain against HEAD failed on the very first anchor.
HEAD reads `42.4% -> 99.2%` for the sparsity step, whereas the pre-truncation
source read `42.4% -> 100%`. So HEAD predates **the cert-2/delta saturation
correction as well** ("saturation is complete at 625/625, not 620/625"), not
just the SS II.C rewrite. The uncommitted working tree held remediation from
*several* prior sessions.

**Consequence:** anchor-by-anchor replay is the wrong method — it would
produce a hybrid that is neither the old nor the new paper. The
reconstruction must be a single deliberate pass with the preserved PDF as
the source of truth and the (intact) test suite as acceptance.

**Current file state: HYBRID.** It is the committed base plus one recovered
subsection (SS II.C, which removed a live false claim and is independently
verified against the PDF text). A `% SOURCE STATE WARNING` block is prepended
to the .tex so the state cannot be mistaken for current.

**Do not compile, cite, or QA this file until the reconstruction lands.**

---

## RESOLVED (same session)

Reconstruction complete, done as a single deliberate pass against the
preserved PDF rather than anchor-by-anchor replay.

**Restored:**
- SS II.C — the PI-delegated LARGE: the false `1/(8Z^2)` "matches" claim
  replaced by the measured `Z^-0.85` decay (12x-79x above the naive
  estimate), with the explicit negative and the "no quantitative bridge
  from kappa" statement.
- SS III — full rebuild: 107/625 = 17.1%, the graded theta-profile
  (0.427/0.702/0.830/1.000) replacing the retired step function, the
  generator-dependence and cutoff-robustness paragraphs (smallest nonzero
  0.0172), n_max=4 at 57,700 (7.12%) with the canonical-unique 15,293
  (7.07%), the basis-size trend REVERSAL (~M^-0.49), and both dated
  correction records.
- SS IV — the n_max provenance note (values had come from n_max=3/4 runs
  while the text said n_max=2).
- SS V — I_cv table extended per element (Li 0.228, Be 3.7e-3, B 1.7e-3,
  C 6.6e-4, N--F < 1e-14 exact), the "constrains the basis rather than
  justifying the factorization" scoping, molecular S_core 0.008 /
  S_bond 0.330 / ratio ~40.
- The N--F hub paragraph with the degeneracy caveat, today's final values
  (0.148/1.063/0.311, maxima 0.69/1.34/1.32, degeneracies 4/7/6) and the
  **ceiling retraction** (ln 8 / ln 16 "attained exactly" was an artifact
  of the diagonal-only entropy routine).
- Abstract, provenance-tier paragraph and conclusions synced; all retired
  `O(Q^2.5)` mentions replaced with the exact-linearity statement.

**Verification:** compiles three-pass clean (7 pp, 0 errors, 0 undefined
refs); `tests/test_paper26_entanglement.py` +
`test_paper26_molecular_entanglement.py` + `test_paper27_entropy_locus.py`
= 12 passed, 1 skipped; C10/C11/C16/C17/C18/C19 all PASS on group6.
Content coverage cross-checked by `pdftotext` against the preserved target.

**Note on values:** the preserved PDF is the *intermediate* state
(post-ERI-fix, pre-Condon-Shortley-order-fix). Where the two disagree the
**tests** are the spec, since they carry the final post-order-fix
measurements. This affects the I_cv threshold (exact closure moves from
Z>=5 to Z>=7), the degeneracies (4/7/6) and the hub values.

**Standing rule written:** `memory/feedback_never_truncate_on_read.md`.
