"""DELTA #2 -- correct the two false statements DELTA #1 left in the chronicle.

The CHANGELOG is the corpus's canonical chronicle, so a false MEASUREMENT
sitting in it is exactly the zombie class the C16 machinery exists to stop.
Two, both from the v5.11.7 entry:

(1) "It PLATEAUS."  It does not.  A reviewer sampled npts=400000 -- a point
    nobody had taken -- and the column climbs past the flat spot; an
    independent PM ladder reproduced it digit for digit:

        100000  1.55917e-07    --
        250000  1.56764e-07  1.005x   <- the "plateau"
        400000  1.93056e-07  1.232x
        600000  2.05980e-07  1.067x

    The entry's OWN table already printed 2.0598e-07 at 600k beneath the word
    PLATEAU, as a "0.76x step", unreconciled.  The 1.005x step is a crossing:
    the quantity is relative to a c=5 box whose own truncation (7.69e-08
    absolute, against c=3's 1.29e-07) is comparable and partially cancels
    there.  The CONCLUSION the entry drew -- App. A's ~1e-7 vindicated, the
    revert correct -- survives, because every discretization agrees on the
    ORDER.  The third digit does not: uniform mesh 1.29e-7, graded mesh
    2.16e-7, a factor 1.7 apart.

(2) "The per-state threshold sits 4.6x below the true separation."  At the
    test's ACTUAL configuration the margin is 3.5x (1.764% against a 0.5%
    threshold); 4.6x is from a different, K=100 configuration.

Idempotent.
"""
from __future__ import annotations

import sys

CL = "CHANGELOG.md"

EDITS = [
    ("plateau", "CORRECTED 2026-09-13 (DELTA #2): IT DOES NOT PLATEAU",
     'It PLATEAUS. The paper\'s figure is right, the original "VINDICATED" wording was '
     'right, and the correction was the error. Restored with the full ladder, which is '
     'stronger evidence than the original had. `test_c3_sits_below_the_quadrature_floor` '
     '-- a guard that PASSED while naming a false mechanism -- replaced by one asserting '
     'convergence.',
     '~~It PLATEAUS.~~ **CORRECTED 2026-09-13 (DELTA #2): IT DOES NOT PLATEAU.** A '
     'reviewer sampled `npts=400000`, a point nobody had taken, and the column climbs '
     'straight past the flat spot -- `1.55917e-07 / 1.56764e-07 / 1.93056e-07 / '
     '2.05980e-07` at 100k/250k/400k/600k -- and an independent PM ladder reproduced '
     'that digit for digit. **The table above already said so and was not read:** its '
     'last row rises 31% over the row marked `0.99x`, printed as a "0.76x step" beneath '
     'the word PLATEAU. The `1.005x` step is a CROSSING, not a plateau -- this quantity '
     'is relative to a `c=5` box whose own truncation (7.69e-08 absolute, against '
     '`c=3`\'s 1.29e-07) is comparable and partially cancels there. **What survives is '
     'the conclusion, and it is the part that mattered:** App. A\'s `~1e-7` is '
     'vindicated as an ORDER by every route (uniform mesh 1.29e-7, graded mesh 2.16e-7, '
     'in-repo relative 1.6-2.1e-7), and the 2026-09-13 revert was right. **What falls is '
     'the third digit and the mechanism** -- the two discretizations differ by ~1.7x, so '
     'no converged value should be quoted, and this entry quoted one. The replacement '
     'guard asserted a plateau over the only adjacent pair in the ladder that satisfies '
     'it, and a reviewer fired it by planting MORE resolution; it is now replaced again '
     'by one asserting the order and the non-collapse. **Third grid-convergence trap of '
     'this arc, second to reach the record, and this one was read too SHORT in the '
     'opposite direction from the first.**'),

    ("margin", "3.5x at the test's own configuration",
     'The per-state threshold sits 4.6x below the true separation and infinitely above '
     'an exactly-zero noise floor.',
     'The per-state threshold sits 3.5x below the true separation (1.764% against a 0.5% '
     'threshold) and infinitely above an exactly-zero noise floor. *(Corrected '
     '2026-09-13: this said 4.6x, which is the margin at a different, K=100 '
     'configuration, not at the one the test runs -- 3.5x at the test\'s own '
     'configuration.)*'),
]


def main() -> int:
    with open(CL, encoding="utf-8") as fh:
        t = fh.read()
    applied, skipped, missed = [], [], []
    for name, marker, old, new in EDITS:
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
