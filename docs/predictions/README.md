# Pre-registered predictions

One file per prediction set, named `YYYY-MM-DD_<topic>.md`. Each file is committed on its own,
before any computation that scores it. After scoring, a disposition is appended; the predictions
above it are never edited.

## Required sections

1. **Header.** Registered date, what prompted it, what it extends, commit status.
2. **The question.** What is released or tested, in a few sentences.
3. **What I know before computing (disclosed).** Everything the author already knows or has
   derived by hand, and where it came from. Nothing computed.
4. **Model.** Exactly what will be computed, with conventions.
5. **Predictions.** Each with an ID, a statement, and a falsifier: the outcome that refutes it.
6. **Novelty status, per prediction (added 2026-10-19, PI-directed).** One of:
   - **KNOWN**: the outcome is in the literature. Cite the primary source and exact location.
     A hit is a check of the method, not a finding.
   - **PROBABLY KNOWN**: a strong lead that was not confirmed. Say what was searched.
   - **NOT FOUND**: searched, nothing found. Say where and how thoroughly. Only a hit on a
     NOT FOUND prediction can count as a finding, and only after the search is repeated at
     scoring time.
   Set this before computing, so the hit rate on unknown outcomes cannot be adjusted afterwards.
7. **Controls.** At least one positive control with a known answer, and a negative control
   where one makes sense.
8. **Disposition** (appended after scoring). Each prediction HELD / FAILED / PARTLY, with the
   evidence, and its final class: **FINDING** (held, NOT FOUND), **CHECK** (held, KNOWN),
   **MISS** (failed), or **UNCLEAR**.

## Why the novelty field exists

A correct prediction of a known result checks the bookkeeping. A correct prediction of an
unknown result is a finding. Without the field the two look the same in the record, and the
program's hit rate on new outcomes cannot be measured.
