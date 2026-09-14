"""DELTA #2 -- correct the chronicle's own 1.7x mesh-disagreement claim.

Both v5.11.7's correction block and the v5.11.8 entry state that the uniform
and graded meshes differ by ~1.7x on the converged c=3 value and that the third
digit is therefore undetermined.  Measured to convergence on both meshes, that
is false: they agree to 0.43%, the converged value is 2.15e-07, and the 1.29e-07
figure was a reference resolved 4x more coarsely than the quantity it was
referencing.

This is the FOURTH mechanism offered for one number, and the chronicle must not
carry the third and fourth as if they were settled.

Idempotent.
"""
from __future__ import annotations

import sys

CL = "CHANGELOG.md"

CORRECTION = (
    "**The conclusion survives and the value IS determined.** App. A's `~1e-7` is "
    "vindicated at **2.15e-07**. Both meshes converge and they AGREE: uniform reaches "
    "`2.14665e-07` at 1.5M points (step 1.010), graded `r=Rt^2` reaches `2.15582e-07` "
    "at 240k (step 1.010), the graded absolute route gives `2.16319e-07`, and the two "
    "meshes' c=5 ABSOLUTE values agree to `7.5e-10`. *(Corrected 2026-09-13, the fourth "
    "reading of this one number: an earlier version of this paragraph said the meshes "
    "\"differ by ~1.7x, so no converged value should be quoted\". That came from "
    "referencing c=3 against a c=12 box at the SAME point count -- 4x coarser spacing "
    "than the quantity being resolved -- and the resulting ratio is not even monotone "
    "in npts: `2.86e-07 / 1.29e-07 / 2.02e-07` at 250k / 600k / 1.5M. The ladder was "
    "stopped at 600k while still climbing at 1.067x, and \"still moving\" was read as "
    "\"cannot be determined\".)* The guard built on the plateau pinned `(100000, "
    "250000)`, **the only adjacent pair in the ladder under a 10% window** (0.54%, "
    "23.15%, 6.69%), and the reviewer fired it by planting MORE resolution -- a change "
    "with no physics in it, verbatim the defect `delta_fix_05` was written to remove "
    "from the guard before it. It is now replaced by an order-plus-non-collapse leg at "
    "the resolutions it can afford, plus a graded-mesh leg that reaches and pins "
    "2.15e-07 at 240k -- converged, and ~6x cheaper than the uniform route. "
    "**Four mechanisms have now been offered for this number and only the paper's own "
    "claim survived all four.**"
)

OLD_A = (
    "**The conclusion survives; the mechanism and the third digit do not.** App. A's "
    "`~1e-7` is vindicated as an ORDER by every route -- uniform mesh `1.29e-7`, graded "
    "mesh `2.16e-7`, in-repo relative `1.6-2.1e-7` -- and the revert was right. But the "
    "two discretizations differ by ~1.7x, so no converged value should be quoted and "
    "v5.11.7 quoted one. The guard built on the false mechanism pinned `(100000, "
    "250000)`, **the only adjacent pair in the ladder under a 10% window** (measured: "
    "0.54%, 23.15%, 6.69%), and the reviewer fired it by planting MORE resolution -- a "
    "change with no physics in it. That is verbatim the defect `delta_fix_05` was "
    "written to remove from the guard before it. Replaced by one asserting the ORDER "
    "and the NON-collapse, which is what is robust, naming both wrong answers."
)

OLD_B = (
    "**What survives is the conclusion, and it is the part that mattered:** App. A's "
    "`~1e-7` is vindicated as an ORDER by every route (uniform mesh 1.29e-7, graded "
    "mesh 2.16e-7, in-repo relative 1.6-2.1e-7), and the 2026-09-13 revert was right. "
    "**What falls is the third digit and the mechanism** -- the two discretizations "
    "differ by ~1.7x, so no converged value should be quoted, and this entry quoted "
    "one."
)

NEW_B = (
    "**What survives is the conclusion, and it is the part that mattered:** App. A's "
    "`~1e-7` is vindicated, at a converged **2.15e-07** -- uniform `2.14665e-07` at "
    "1.5M, graded `2.15582e-07` at 240k, agreeing to 0.43% -- and the 2026-09-13 revert "
    "was right. *(Corrected 2026-09-13: an earlier version of this line said the two "
    "discretizations differ by ~1.7x so no converged value should be quoted. They do "
    "not differ; the 1.29e-07 was a c=12 reference at the same point count as the c=3 "
    "box it referenced, hence 4x coarser, and its ratio is not monotone in npts.)*"
)


def main() -> int:
    with open(CL, encoding="utf-8") as fh:
        t = fh.read()
    applied, skipped, missed = [], [], []
    for name, marker, old, new in [
        ("v5.11.8-block", "the fourth\nreading of this one number" if False else
         "the fourth reading of this one number", OLD_A, CORRECTION),
        ("v5.11.7-block", "They do not differ; the 1.29e-07 was a c=12 reference", OLD_B, NEW_B),
    ]:
        if marker in t:
            skipped.append(name)
            continue
        if t.count(old) != 1:
            missed.append((name, t.count(old)))
            continue
        t = t.replace(old, new)
        applied.append(name)
    with open(CL, "w", encoding="utf-8") as fh:
        fh.write(t)
    for n in applied:
        print(f"  ok    {n}")
    for n in skipped:
        print(f"  skip  {n} (already applied)")
    for n, c in missed:
        print(f"  MISS  {n}: anchor count={c}")
    print(f"applied {len(applied)}, skipped {len(skipped)}, MISSED {len(missed)}")
    return 3 if missed else 0


if __name__ == "__main__":
    sys.exit(main())
