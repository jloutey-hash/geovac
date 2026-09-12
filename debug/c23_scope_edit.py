"""C23 scope change (2026-09-12, PI-approved): run near authorship, not only at
certification.  Driver: C23 run #1's two best catches were its two NEWEST claims
(Prop D written 2026-09-11, the flat-limit degeneracy written 2026-09-12), both
caught within 48 hours of authorship; the five older ones had survived three
DELTA runs and a FULL run."""
from pathlib import Path

P = Path("docs/qa/criteria.md")
s = P.read_text(encoding="utf-8")

OLD = ("**Scope.** Runs on FULL certification runs, not on every DELTA — the cost is a "
       "literature scan per keystone claim. A DELTA run inherits the last FULL run's C23 "
       "verdicts unless the claim itself changed.")

NEW = ("**Scope (revised 2026-09-12 after run #1, PI direction).** C23 has two triggers, and the "
       "second is the important one.\n\n"
       "1. **At FULL certification**, over the target's `[SYMBOLIC]`-tier claims. A DELTA run "
       "inherits the last FULL run's verdicts unless the claim itself changed.\n"
       "2. **At authorship**, on any *new* claim matching a priority signature above — a clean "
       "closed-form constant, an external field entered sideways, or a derivation under a page. "
       "One scan, at the moment the claim is written.\n\n"
       "*Why trigger 2 exists, measured.* Run #1 audited eight claims and returned six PRIOR ART. "
       "Its two most valuable catches were its two **newest** claims — a \"Proposition\" written "
       "the day before the run, which was the Löwdin symmetry-preservation property known since "
       "Slater–Koster (1954), and a degeneracy structure written the same morning, which was the "
       "*flat limit* of the radial-basis-function literature. Both were caught within 48 hours of "
       "being written. The other four had survived three DELTA runs and a FULL run, which is the "
       "cost of certification-only scoping: a rediscovery sits in a paper, gets cited by its own "
       "corpus, accretes dependents, and the eventual re-attribution becomes a sweep rather than "
       "an edit. Catching it at authorship costs one scan and no dependents.\n\n"
       "*Cost control.* Trigger 2 is per-claim, not per-paper, and only for claims matching a "
       "priority signature — not every measured number. If a sprint produces no `[SYMBOLIC]`-tier "
       "claim, it fires not at all.")

assert s.count(OLD) == 1, "C23 scope paragraph not found"
P.write_text(s.replace(OLD, NEW), encoding="utf-8")
print("criteria.md: C23 scope revised to two triggers")
