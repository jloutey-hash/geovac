"""KNOWN-GAP 1/5, step 6 -- the three "Proposition D" loci the fixed pattern found.

Sequence worth recording, because it is the gate-first discipline paying out
exactly as qa.md predicts:

  1. wrote the entry with pattern `Proposition~?D` (LaTeX tilde only);
  2. the discrimination proof caught that it misses the PLAIN-TEXT form;
  3. widening it to `Proposition[~\\s]?D` immediately surfaced THREE live loci
     the narrow pattern had walked straight past.

"A locus list from a grep is exhaustive; a locus list from memory is a sample."

All three carry the retired LABEL, not the retired mathematics -- the result is
true and its backing is sound.  What was demoted (C23 run #1) is its status as
a Proposition of ours: it is Loewdin symmetry preservation specialised to the
l grading, known since Slater-Koster (1954).  The paper already dropped the
label (0 occurrences); the claim matrix and the backing test had not.

Fixed by REPLACING the label with the attributed description rather than
exempting it -- an exemption would keep the retired framing alive under a flag.
The test FUNCTION name is left alone: it is an identifier the claim matrix and
C13/C22 resolve against, and it does not match the pattern.

Idempotent.
"""
from __future__ import annotations

import sys

EDITS = [
    ("docs/claim_test_matrix.md",
     "| 60 | §molecular [SYMBOLIC] — Proposition D: a block-diagonal congruence",
     "| 60 | §molecular [SYMBOLIC + PRIOR ART] — the block-diagonal congruence "
     "result (**re-attributed 2026-09-12, C23 run #1: this is Löwdin "
     "symmetry-preservation specialised to the `l` grading, known since "
     "Slater-Koster 1954; what the paper claims is the `l`-vs-`m` application, "
     "not a new proposition**): a block-diagonal congruence"),

    ("tests/test_paper60_kms_attribution.py",
     "  6. Proposition D: a block-diagonal congruence cannot orthogonalize a metric\n"
     "     that is not itself block diagonal -- the l-selection loss is independent\n"
     "     of conditioning.",
     "  6. The block-diagonal congruence result: a block-diagonal congruence cannot\n"
     "     orthogonalize a metric that is not itself block diagonal -- the\n"
     "     l-selection loss is independent of conditioning.  PRIOR ART (C23 run #1,\n"
     "     2026-09-12): this is Loewdin symmetry preservation specialised to the l\n"
     "     grading, known since Slater-Koster (1954).  What is claimed here is the\n"
     "     l-versus-m application, not a proposition of ours."),

    ("tests/test_paper60_kms_attribution.py",
     "# ------------------------------------------------------------- Proposition D",
     "# ------------------------- block-diagonal congruence (Loewdin/Slater-Koster)"),
]


def main() -> int:
    applied = 0
    for path, old, new in EDITS:
        with open(path, encoding="utf-8") as fh:
            t = fh.read()
        if old not in t:
            print(f"  skip {path} (anchor absent -- already applied?)")
            continue
        if t.count(old) != 1:
            print(f"  MISS {path}: count={t.count(old)}")
            return 2
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(t.replace(old, new))
        applied += 1
        print(f"  ok   {path}")
    print(f"applied {applied}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
