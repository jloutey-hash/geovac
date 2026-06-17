# Paper 57 P5 "packing-reachability" validation — diagnostic memo (2026-06-16)

## Verdict

**INHERENTLY-A-RESTATEMENT.** P5's headline 98.3% accuracy cannot be made
non-circular. The `packing_reachable` column in `debug/principle_hunt_audit.py`
is a one-to-one relabel of the F/C `status` column it then "predicts," and the
only criterion that could make a witness-existence check discriminate forced
from free — reading whether the cited witness is a *forcing derivation* or a
*non-selection theorem* — IS the F/C determination. There is no symmetric,
F/C-blind witness-existence check that recovers the 98.3%.

Recommendation: **validate before reducing is the wrong frame for P5 — there is
nothing to validate.** Paper 57 should present the 98.3% as an
**internal-consistency check**, not a discovery, with an explicit
circularity caveat. The genuinely validatable structural content lives in P1
(multi-focal depth), P2 (period class), and the two-family decomposition — NOT
in P5. Those should carry the discriminator narrative; P5 should be demoted to
"the consolidated restatement of the F/C boundary," which is its honest role.

---

## What was checked (artifacts built)

Driver work was done against the existing `debug/principle_hunt_audit.py`
catalogue (60 entries, mirrors `docs/forced_free_seam.md`). No new production
code; three diagnostic probes run inline.

### Probe 1 — is `packing_reachable` independent of `status`? (cross-tab)

Cross-tabulating the two hand-assigned columns over all 60 entries:

| status | packing_reachable | count |
|:------:|:-----------------:|:-----:|
| F | yes | 35 |
| A | conditional | 1 |
| C | no | 23 |
| C | conditional | 1 (I3) |

The map `status → packing_reachable` is **perfectly deterministic** with exactly
one deliberate exception (I3, tagged `conditional` rather than `no`). And I3 is
*precisely* the single P5 misclassification — it is `conditional` so that it
graduates to F if the open Hopf-base identification lands. So P5 reproduces the
F/C label on 59/60 entries by construction, and disagrees on the 1 entry that
was hand-perturbed off the diagonal. The 98.3% measures the consistency of two
columns that are one-to-one by design. **This is the tautology, made numeric.**

### Probe 2 — can a symmetric (polarity-blind) witness-existence check discriminate?

The genuine-validation path requires determining `packing_reachable` by an
*independent* criterion: "does a cited theorem / bit-exact test actually exist
for this entry?", judged the same way for every entry.

Tested directly: **every one of the 60 entries carries a cited witness in the
corpus.** The forced entries cite forcing derivations (Paper 32 §VIII Forced-Count,
Paper 18 master Mellin engine, Paper 38 GH convergence, Paper 51 gravity closed
forms, …). The calibration entries cite *non-selection theorems* that exist in
the corpus just as solidly:

- F2 (N_gen): `thm:n_gen_non_selection` — 12 refs in Paper 32
- F3 (KO-dim): `thm:ko_dim_non_selection` — 5 refs in Paper 32
- E6 (K-rule): `thm:no_single_mechanism_K` — present in Paper 32
- D5/D6 (cutoff): `thm:cutoff_function_external` — present in Paper 32
- F1 (Yukawa): Sprint H1 non-selection theorem + 162-cell PSLQ memo

A **polarity-blind** witness-existence check ("is there *a* cited theorem/test?")
therefore returns TRUE for all 60 entries → predicts F for all → accuracy
36/60 = **60.0%** (the base rate of forced entries). It cannot discriminate at
all.

### Probe 3 — what makes the check discriminate, and is that thing independent?

The only thing separating a "yes" from a "no" is the **polarity** of the
witness: a forcing derivation vs. a non-selection theorem. Reading that polarity
is the F/C judgment itself (the catalogue's own `Status conventions` in
`docs/forced_free_seam.md` define F = "witness: a theorem or bit-exact
verification [that forces it]" and C = "witness: a non-selection theorem or
documented negative search"). The Sprint C2 memo states this outright
(`debug/sprint_c2_principle_hunt_memo.md` line 84): the tag "was determined by
examining the corpus for witness derivations (forced entries) and explicit no-go
theorems / non-selection theorems (calibration entries)" — i.e. the criterion
applied is *already conditioned on which side the entry is known to be on*.

**Conclusion: there is no F/C-blind operationalization of `packing_reachable`.**
Either it's polarity-blind (60% — useless) or polarity-aware (98.3% — but
polarity-aware = F/C-aware = circular). P5 cannot occupy a middle ground.

---

## Assessment of the paper's own three "reasons P5 is non-trivial" (§6.2)

The paper (§6.1–§6.2 in the .tex; reproduced in the C2 memo) offers three
defenses. Each is addressed:

1. **"The tagging procedure works on new entries."** True but does not rescue
   the *reported number*. Classifying a genuinely new entry still requires
   reading the polarity of whatever witness the corpus offers — i.e. deciding
   whether it's forced or free. That is a re-derivation of the F/C label, not an
   independent test of it. A pre-registered blind-prediction protocol (classify
   N new physical inputs as F/C *before* checking the corpus verdict, then
   score) WOULD be a genuine test — but (a) it tests the *analyst's ability to
   apply the F/C definition*, not a structural predictor P5 distinct from F/C,
   and (b) it has not been run, so it cannot back the 98.3% headline. The 98.3%
   is retrospective on a fully-labeled catalogue.

2. **"Predictive content (I3 graduates if Hopf-base lands; sibling axiom
   re-classifies F2)."** These are predictions about *re-labeling under changed
   inputs*, which any consistent definition has. "If the Higgs direction is
   identified with the Hopf base, it becomes forced" is "if we find a forcing
   derivation, it's forced" — definitional, not P5-specific.

3. **"Failure-mode decomposition matches independent prior work
   (multi-focal-wall pattern + η-trivialization)."** This is the ONE genuinely
   non-circular thread — but the independent content lives in the **two-family
   decomposition and the MF / period columns**, not in P5. Those columns are
   tagged separately from `packing_reachable` and DO discriminate (see below).
   P5 merely inherits the split; it is not what carries it.

---

## The genuinely validatable content (what to lead with instead)

The MF (multi-focal depth) and P (period class) columns ARE tagged independently
of `packing_reachable`, and they carry real structural signal scored against the
F/C tag:

- **P1 (MF=1 ⇒ forced): 86.7%.** Catches family 1 cleanly.
- **P2 (period ≠ outside ⇒ forced): 86.7%.** Catches family 1 cleanly.
- **P4 (compactness inheritance): 96.7%.**

And the two-family decomposition has cleanly different *independent-column*
signatures on the 24 calibration entries:

| Family | n | MF | period=outside |
|:-------|:-:|:--:|:--------------:|
| multi-focal composition | 17 | all MF=2 | 15/17 |
| inner-factor input data | 7 | all MF=1 | 2/7 |

This split — uniform MF=2 / period-outside for family 1 vs. uniform MF=1 /
period-not-outside for family 2 — is a real, testable pattern (and it matches
the independently-named multi-focal-wall pattern and η-trivialization tier).
**This is the discovery.** P5 is the restatement that sits on top of it.

---

## Recommended wording for Paper 57

Two concrete edits (NOT applied — for PI decision):

### (a) Recast §5.5 "P5" and Table 1 — present P5 as the consolidated F/C
boundary, not a fifth discriminator at 98.3%.

Replace the "98.3% accuracy" framing with:

> **P5 (packing-reachability) is the consolidated statement of the forced/free
> boundary itself, not an independent discriminator of it.** An entry is tagged
> packing-reachable exactly when the corpus carries a forcing derivation for it,
> and packing-unreachable exactly when the corpus carries a non-selection
> theorem or documented negative search. Because "forcing derivation exists" is
> the definition of Forced (F) and "non-selection theorem / negative search"
> is the definition of Calibration (C), the packing-reachable tag is — by
> construction — the F/C label under a different name. The 98.3% agreement with
> the F/C status (59/60, the lone exception I3 deliberately tagged
> *conditional* to encode an open identification) is therefore an
> **internal-consistency check**: it verifies that the catalogue applied its own
> F/C conventions consistently across 60 entries. It is not a discovered
> structural predictor, and we do not present it as one. The genuine
> discriminator content is carried by the multi-focal-depth (P1) and
> period-class (P2) axes and by the two-family decomposition
> (Observation~\ref{obs:two_families}), all three of which are tagged
> independently of the F/C status and which P5 organizes rather than predicts.

In Table 1, either drop the P5 row or annotate it: "P5 — consolidated F/C
boundary (definitionally the status column); 98.3% = internal-consistency check,
not independent accuracy."

### (b) Rewrite §6.1 "Is P5 the principle, or tautological?" to state the
resolved verdict rather than leaving it open.

Replace the "three reasons to think P5 is non-trivial" with the honest finding:

> A direct check settles this. We attempted to determine packing-reachability
> by an F/C-blind criterion — "does the corpus carry a cited theorem or
> bit-exact test for this entry?" — applied identically to every entry. It does
> not discriminate: all 60 entries carry such a witness (the forced entries via
> forcing derivations, the calibration entries via the non-selection theorems of
> Paper~32 §VIII and the negative-search memos), so an F/C-blind
> witness-existence check predicts "forced" for all 60 and scores at the 60%
> base rate. The 98.3% is recovered only by reading each witness's *polarity*
> (forcing vs. non-selection), and that polarity reading is the F/C judgment
> itself. We therefore conclude that P5 **cannot be made non-circular**: it is
> the forced/free boundary restated, and its accuracy against the F/C label is
> an internal-consistency check. This does not weaken the catalogue — the
> structural discovery is the two-family decomposition and its independent
> MF / period signatures (§5.1–§5.4), which P5 consolidates. The open question
> is correspondingly sharpened: not "is P5 tautological" (it is), but "does the
> packing-reachable boundary admit a meta-theorem" — i.e. a structural
> characterization that would let one decide forcedness *without* first locating
> the forcing-or-non-selection witness. That meta-theorem is the multi-year
> NCG-research target of §6.4.

Also update the abstract: change "The fifth candidate … hits 98.3% accuracy" to
"The fifth candidate (packing-reachability) consolidates the forced/free
boundary; its 98.3% agreement with the F/C status is an internal-consistency
check confirming the catalogue applies its conventions uniformly, not an
independent discriminator (the independent signal is carried by multi-focal
depth, period class, and the two-family decomposition)."

---

## Honest scope of this memo

- The verdict is structural, not numerical-precision-dependent: it rests on the
  cross-tab (Probe 1) and the polarity argument (Probes 2–3), both of which are
  exact statements about the catalogue, not estimates.
- The two-family decomposition's independent content (P1/P2/MF/period) is
  genuine and is NOT being demoted — only P5 is.
- The cited non-selection theorems (`thm:n_gen_non_selection`, etc.) genuinely
  exist in Paper 32 (verified by grep); their existence is in fact what *defeats*
  the symmetric witness-existence check, since they make the calibration entries
  witness-backed too.
- This memo recommends wording only; per the task it does not edit any paper.

## Files
- `debug/paper57_p5_validation_memo.md` — this memo
- (probes were run inline against `debug/principle_hunt_audit.py`; no new driver
  needed — the existing catalogue already exposes both the `status` and
  `packing_reachable` columns whose 1-to-1 relation is the whole finding)
