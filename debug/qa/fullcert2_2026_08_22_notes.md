# /qa FULL certifying run #2 — 2026-08-22

**Target:** papers 58, 59, 60, 19 (group2) + 56, 24 (group3) + the group2 and
group3 syntheses.
**Shape:** FULL — all dimensions, whole-paper enumeration, completeness-critic.
**Seed key:** `debug/qa/fullcert2_seed_key.json` — 13 seeds, 6 controls.
**Criteria:** `docs/qa/criteria.md` C1–C19, frozen before dispatch. Nothing
relaxed; the C19 gate and the registry-discrimination rule added since the
previous run both tighten.

## Precondition note (stated, not waived)

The run-shapes rule makes a **CLEAN-DELTA** the precondition for a full
certifying run. The preceding delta returned **DEFECTS** — eleven genuine
incompletions plus three defective gates/registry entries — and all were
remediated afterwards, so the current diff has **never been seen by a
calibrated panel**. This run was fired at PI direction without an intervening
delta cycle. The mitigation applied: every reviewer prompt names the recently
rewritten loci explicitly and instructs the reviewer to weight them, so the
delta scope is folded into the full run rather than assumed clean.

## Step 1 — deterministic layer (whole target, before dispatch)

| Gate | Result |
|:--|:--|
| C5/C12 K-label | PASS |
| C10 compiles + undefined refs (clean-start) | PASS — 8/8 documents |
| C11 internal titles | PASS |
| C13 paper↔test refs | PASS |
| C14 paper↔file refs | PASS |
| C15 inline arXiv IDs | PASS |
| C16 retracted terms / zombies | PASS |
| C17 headline-number registry | PASS (selftest 6 passed) |
| C18 duration language | PASS |
| C19 eaten-escape corruption | PASS (selftest 10 pos / 5 neg, 0 failures) |

Zero build strays. Paper 56's three warnings are font-substitution only.

## Step 3 — seeds (worktree only)

13 seeds, **classes deliberately varied** from the previous FULL and DELTA runs
so nothing is pattern-learnable.

| Dimension | Seeds | Classes |
|:--|:--:|:--|
| code | 6 | structural identity → tautology; self-referential comparison; gutted magnitude pin; **gutted tripwire**; exactness downgraded to float; tautological test |
| claims | 2 | status overstatement of a paper frontier ([OPEN] → given); **hard-prohibition touch (K relabelled "derived")** |
| citation | 3 | arXiv ID resolving to a different paper ×2; nonexistent theorem number |
| synthesis | 2 | open question asserted settled; **refuted keystone re-asserted (zombie)** |

Two are deliberately cross-cutting: the K-prohibition seed should be caught by
the claims reviewer *and* by the C5/C12 deterministic screen, and the gutted
tripwire tests whether a reviewer notices a guard that reports green while
being unable to fire — the exact defect found in a registry entry last cycle.

6 controls, all drawn from corrections made in this session (the proved
endomorphism closed form, the maxima-not-range restatement, the two convention
percentages, the softened external attributions, the census provenance with its
n_max tag, and the corrected quadrature paragraph). Several are *softened* or
*narrowed* forms — flagging them as overclaims would be a false positive, which
is the specificity axis this run most needs to measure.

## Step 4 — panel

11 reviewers, all Opus, all path-pinned to
`C:/Users/jlout/Desktop/geovac-qa-fc2` and forbidden the real corpus:

- **code** (1 per paper): 58, 59, 60, 19, 56, 24 — each instructed to
  independently re-derive the published constants rather than trust the tests'
  arithmetic, and to verify regression guards can actually fire
- **claims**: A = 58/59/60, B = 19/56/24 (B carries the enumerate-every-K mandate)
- **citation**: A = 58/59/60/19, B = 56/24 + both syntheses — both carrying an
  explicit mandate to digit-audit the classic-era entries, the class that
  defeated a previous pass
- **synthesis**: both documents, C9-currency exercised mechanically

Completeness-critic dispatched after the panel returns.

---

*(Calibration scorecard, findings and verdict appended below.)*

---

## The run's most consequential finding: an undocumented truncation in production code

**Found by** the Paper 24 code reviewer; **verified independently** by the PM
(direct computation + derivation); **fixed at PI direction** after a review of
the introducing commit's context.

### The defect

`geovac/nuclear/moshinsky.py::lab_to_relative_matrix_element` carried

```python
if N_bra != N_ket:
    return 0.0
```

zeroing every two-body matrix element between states of different total HO
quantum number. Paper 24 justified it as *"any central V(r_12) preserves
N_tot"*. That is false.

### Why it is false

The Moshinsky **bracket** conserves N — that is a property of the HO coordinate
transformation. The **potential matrix element** does not. A central V(r_rel)
is diagonal in the CM quantum numbers and in l_rel but COUPLES different
relative-n. Writing the element out, the bra bracket forces
2n + l + 2N_cm + Lam = N_bra and the ket bracket forces
2n' + l + 2N_cm + Lam = N_ket, hence n - n' = (N_bra - N_ket)/2. The discarded
elements are exactly those with n != n', and they are large:

| element | value (Minnesota S=0, b=1) |
|:--|--:|
| <n_rel=0,l=0 \| V \| n_rel=0,l=0>  (kept) | -0.81 MeV |
| <n_rel=0,l=0 \| V \| n_rel=1,l=0>  (zeroed) | **+17.21 MeV** |
| <n_rel=0,l=0 \| V \| n_rel=2,l=0>  (zeroed) | **+17.82 MeV** |

Through the production API the same shows up as lab-frame couplings of
**+7.06** and **+8.15 MeV** against a **-0.55 MeV** diagonal.

### Why it was fixed rather than escalated further

§9 makes a suspected keystone bug a PI-raise, and it was raised. The PI directed
a fix *after reviewing the introducing commit's context*, which is what settled
it. That review found:

- introduced undocumented in `8d692a0` (v2.7.0, 2026-04-12);
- **Paper 23 — which depends on this module — never states an N_tot
  restriction**, only "Moshinsky--Talmi brackets";
- `docs/nuclear_electronic_embedding_spec.md` from the same commit has no such
  language.

So it was never a declared model convention. The only prose justifying it is
the false claim in Paper 24.

### The fix

Bra and ket are now decomposed at their own N and matched on CM quantum
numbers — which is what V's structure actually requires. The genuine parity
selection rule (N_bra - N_ket even) is kept as a short-circuit. A
`conserve_N=True` flag reproduces the old behaviour so published numbers stay
checkable, mirroring the `faithful` flag that preserves the pre-fix
double-excitation phase convention.

### Measured impact

**Paper 24 — entanglement-rigidity corollary REFUTED.** All five backing tests
fail against the corrected physics:

| N_max | E0 (MeV) | S | ‖[H_HO,V]‖/‖H_HO‖ |
|--:|--:|--:|--:|
| 2 | 21.6538 | **0.118** | **0.74** |
| 3 | 21.6279 | **0.123** | **0.63** |
| 4 | 21.5442 | **0.136** | **0.67** |

against published S = 0 exactly and commutator < 1e-15. Note the fifth failing
test, `E0_basis_independent`: that independence was **itself the bug's
signature** — the ground state was frozen in a single N block, so E0 was
trivially identical across basis sizes. It now decreases with basis size, which
is correct variational behaviour.

**Paper 23 — resource counts increase ~16%; conclusions survive.**
268 passed, 3 failed; only the Pauli counts moved.

| system | published | corrected | delta |
|:--|--:|--:|--:|
| deuteron | 592 | **688** | +16.2% |
| He-4 (no Coulomb) | 712 | **828** | +16.3% |
| He-4 (with Coulomb) | 712 | **828** | +16.3% |

Qubit counts unchanged (16Q each); magic gaps and spin-orbit unaffected. The
increase is bounded by the N_max-truncated model space, which is why it is
+16% rather than orders of magnitude.

### A memory was carrying the false claim

`memory/ep2b_ho_zero_entropy.md` recorded it in full — and carried its own
tell: it asserted N_tot conservation for "ANY central V" and then, two
sentences later, said Coulomb (also central) breaks it. The two statements
contradict each other, and that contradiction was the signature of the error
sitting in the index of every session. Rewritten as a corrected record.

### Method lesson

Five tests passed for four months against a hard zero produced by restricting
the evaluation object rather than by the physics holding — the TC-qubit-space
failure mode (CLAUDE.md §3) recurring in a new place. **A bit-exact zero on a
quantity that should merely be small is a reason to look harder, not a reason
to relax.**

### Owed (PI decisions, not PM)

- Paper 24: rewrite the entanglement-rigidity corollary and its "structural
  dual of Paper 27's nonzero Coulomb scaling" framing.
- Paper 23: update 592/712 to 688/828 and re-check any downstream sparsity
  narrative.
- Paper 27 §VII.A: the Prediction-1 verification cited EP-2b; re-check.
- The five Paper 24 tests and three Paper 23 resource tests now encode the old
  numbers and must be updated to the corrected values (or pinned under
  `conserve_N=True` with the convention stated).

---

# VERDICT — FULL certifying run #2

**FAIL.**

The panel was calibrated; the target was not clean. Both halves of that
sentence are load-bearing, and the three-way vocabulary is what lets them be
stated separately: this is *"target not done"*, not *"reviewer not
trustworthy."*

## Calibration scorecard

| | result |
|---|---|
| Sensitivity | **13 / 13 planted defects caught** |
| Specificity | **0 false positives on 6 known-good controls** |
| Dimensions exercised | 4 of 4 (code, claims, citations, synthesis) |
| Dimensions uncalibrated | none |

One dimension required recovery. The **citation reviewer for group3
de-calibrated on the Sonnet tier**, then again on the first Opus re-dispatch —
the second time by grading a theorem number GROUNDED while noting *in the same
row* that it had not refetched the entry. Re-dispatched a third time with the
explicit **CONFIRMED / WRONG / UNVERIFIABLE** vocabulary and a rule that
nothing unfetched may be upgraded on plausibility, it caught the planted
number immediately *and* separately CONFIRMED the paper's genuine citation of
the same work. That it discriminated rather than blanket-flagging is what
makes the recovery count.

Because every gating dimension ended calibrated, the verdict is a real
measurement of the target rather than INCONCLUSIVE.

## Why FAIL

Verified MATERIAL defects at seed-free locations, across all four dimensions.
Three are worth naming as classes rather than instances.

**1. A false published result, held up by production code (code dimension).**
Paper 24's two-fermion entanglement-rigidity corollary — and, as this session
traced, Paper 27's EP-2b "rigidity theorem", one of *its* two headline results
at INTERNAL THEOREM tier — asserted that the closed-shell HO ground-state
entropy is identically zero for any central two-body potential. Both the claim
and its stated mechanism are false. The Moshinsky–Talmi **bracket** conserves
total oscillator quanta; the potential **matrix element** does not. The
published zero came from an undocumented `if N_bra != N_ket: return 0.0` in
`geovac/nuclear/moshinsky.py`, which discarded couplings (+17.2 MeV) roughly
twenty times the diagonal terms it kept. Five tests had passed against it for
four months.

**2. Gates that fire correctly and then scope their verdict away.** The
highest-value finding of the run came from the **completeness-critic**, not
from any calibrated dimension: C5/C12 *detected* a planted K-rule violation in
Paper 24 and *discarded its own verdict*, because the gate hard-coded a
six-document TRUNK set and downgraded everything else to advisory — while
§13.5 makes that prohibition corpus-wide. C10 failed the same way in a
different key: it certified by exit code, and `pdflatex -halt-on-error` exits
0 on undefined references, so "C10 green" had been reported indefinitely for a
class it never measured. Twelve documents carried dangling refs the moment a
real check was written.

**3. Remediation introducing defects.** Consistent with the delta-run
experience: several findings sat in passages the previous cycle had just
rewritten, including two TAB-corrupted LaTeX macros that compile silently and
render as literal text — the class C19 now exists to catch.

## Standing caveat on the verdict's reach

C3 requires the provenance tier inline. Two papers in this target carry **zero**
inline tier tags and a third exactly one, so the claims dimension had almost no
C3 surface to audit there. That absence read as clean during the run. It is
not: **a criterion the target does not exercise is unmeasured, not passed.**
Whatever the next run concludes about C3 on this target should be read against
that gap until the tags exist.

## Disposition

Remediation is complete for everything the PM may decide alone, and is recorded
above and in the CHANGELOG:

- `moshinsky.py` N_tot guard removed; `conserve_N` flag added to reproduce the
  retired behaviour on demand.
- Paper 24 corollary and Paper 27 §sec:pred1 rewritten as retractions carrying
  the corrected measurements; abstracts, provenance tiers, section lead-ins and
  conclusions of both brought into line.
- Paper 23 and the group4 synthesis resource counts updated (592→688,
  712→828, 1-norms 342.2→383.7, 466.9→511.8, 462.4→507.2). **The structural
  claim survives exactly**: 828/688 = +20.3%, identical to the retired
  712/592 = +20.3%, so "≈20% more terms for a 12.25× larger Hilbert space"
  is unchanged.
- Eight tests rewritten against corrected physics; a tripwire
  (`test_paper24_n_tot_guard_stays_removed`) now fails if the guard is
  silently reinstated or its default flipped.
- C16 and C17 registries extended per their HARD RULES; both entries verified
  to fire (each caught live zombies before remediation) and to exempt correctly
  after.
- Blast radius measured: 361 passed / 2 failed across the nuclear + paper
  suites, the two failures being exactly the Paper 27 tests since rewritten.
  No other consumer regressed.

**One correction made during remediation, worth recording because it is the
same failure mode as the defect being remediated.** The entropy figures first
written into the retraction (0.118 / 0.123 / 0.136) were wrong —
mis-transcribed from an exploratory run. The reproducible values on the two
papers' *bit-identical* 1-RDM builders are **0.0671 / 0.0716 / 0.0833** nats
(normalized trace-one spatial 1-RDM, natural log). They were caught by checking
the number rather than the test, since the assertion band as first written
passed on both the right and the wrong value. Corrected at all six loci.

## What blocks PASS

Nothing further is owed by the PM. A PASS requires a **clean delta run over
the remediated surface** — which, per the run-shape rule, is the precondition
for firing the next full certifying pass, and per the hard rule added this
session, must re-test *these specific defects* rather than trust this entry.
Two items are genuinely PI-scoped and are not remediated here:

- the corpus-wide C5/C12 scope decision (the gate now covers all `.tex`; the
  question of whether to sweep the advisory AUDIT list is a PI call);
- the inline-tier coverage gap behind the C3 caveat above.
