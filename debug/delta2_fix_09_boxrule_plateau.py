"""DELTA #2 / M1 + M2 -- retire the "plateau", keep the vindication.

WHAT I MEASURED MYSELF (debug/p60_c3_extended_ladder.py, uniform mesh,
n_max=8, l_max=1), after the code reviewer sampled npts=400000 -- a point
nobody had taken:

    npts       c=3 rel        step
    100000   1.55917e-07       --
    250000   1.56764e-07     1.005x
    400000   1.93056e-07     1.232x     <- breaks the "plateau"
    600000   2.05980e-07     1.067x

The 0.99x step at 100k->250k is a CROSSING, not a plateau.  The in-repo
quantity is measured against a c=5 box at the SAME npts, and that reference
carries its own truncation -- measured absolutely against a c=12 box, c=5 is
7.69e-08 while c=3 is 1.29e-07, i.e. COMPARABLE.  Near 250k the two partially
cancel and the difference stalls; past it the cancellation unwinds and the
column climbs.

SO THE PLATEAU IS FALSE, AND IT WAS MY CLAIM.  Worse, my own published table
printed 2.0598e-07 at 600k against 1.5676e-07 at 250k -- a 31% RISE -- as a
"0.76x step" directly beneath the word PLATEAU, and I did not reconcile it.
Third time in this arc that a grid ladder was read too short, and the second
that reached the record.

WHAT SURVIVES, AND IT IS THE PART THAT MATTERED.  The paper's App. A says c=3
gives ~1e-7.  Every route agrees on the ORDER: uniform mesh vs a c=12
reference gives 1.29e-07, the reviewer's graded mesh gives 2.16e-07, the
in-repo relative quantity gives 1.6-2.1e-07 across resolutions.  The appendix
figure is VINDICATED and the 2026-09-13 revert was right.  What is NOT
established is the third digit: the two discretizations disagree by ~1.7x, so
no single converged value should be quoted, and this file previously quoted
one.

THE GUARD.  `test_c3_converges_and_vindicates_the_appendix_figure` asserts a
plateau over the hard-coded pair (100000, 250000) -- the ONLY adjacent pair in
the ladder that satisfies it (measured: 0.54%, 23.15%, 6.69%).  The reviewer
fire-tested it by planting (250000, 600000) -- more resolution, no physics --
and it FIRED.  That is verbatim the defect delta_fix_05 was written to remove
from the PREVIOUS guard.  Replaced by one that asserts the ORDER and the
NON-collapse, which is what is robust, and names both wrong answers.

Idempotent.
"""
from __future__ import annotations

import sys

T = "tests/test_paper60_box_rule.py"
MTX = "docs/claim_test_matrix.md"
CL = "CHANGELOG.md"

DOC_OLD = """and had to be redone at higher resolution.  The full ladder, because three
points were not enough to read it:

    npts      c=1         c=2         c=3        step(c=3)
    12000   3.6130e-02  6.5629e-04  2.5709e-05      --
    40000   3.6166e-02  6.8727e-04  2.1107e-06    12.2x
   100000   3.6169e-02  6.8984e-04  1.5592e-07    13.5x
   250000   3.6170e-02  6.9025e-04  1.5676e-07     0.99x   <- PLATEAU

All three columns converge.  c=3 settles at ~1.6e-7, so it is a genuine
box-truncation quantity and **the paper's appendix figure of ~1e-7 is
vindicated**.

Read the c=3 column to 100k only and it looks like an unconverged quantity
falling 12-13x per step; the plateau appears at the next point.  Both a
reviewer and this author drew the wrong conclusion from the three-point version
before extending it."""

DOC_NEW = """and had to be redone at higher resolution -- twice, because the second reading
was also too short.  The full ladder:

    npts      c=1         c=2         c=3        step(c=3)
    12000   3.6130e-02  6.5629e-04  2.5709e-05      --
    40000   3.6166e-02  6.8727e-04  2.1107e-06    0.082x
   100000   3.6169e-02  6.8984e-04  1.5592e-07    0.074x
   250000   3.6170e-02  6.9025e-04  1.5676e-07    1.005x   <- looks flat
   400000   3.6170e-02  6.9030e-04  1.9306e-07    1.232x   <- it is NOT
   600000   3.6170e-02  6.9031e-04  2.0598e-07    1.067x

c=1 and c=2 plateau.  **c=3 does not.**  The 1.005x step is a CROSSING, not a
plateau:  this column is measured against a c=5 box at the SAME npts, and that
reference carries its own truncation -- absolutely, against a c=12 box, c=5 is
7.69e-08 where c=3 is 1.29e-07, i.e. comparable.  Near 250k they partially
cancel; past it the cancellation unwinds and the column climbs.

**What is established is the ORDER, and that is what the appendix claims.**
Uniform mesh against a c=12 reference: 1.29e-07.  A graded mesh: 2.16e-07.  The
in-repo relative quantity: 1.6-2.1e-07.  So App. A's ~1e-7 is VINDICATED, and
the third digit is NOT determined -- the two discretizations disagree by ~1.7x.
Do not quote a converged value here; an earlier version of this docstring quoted
1.6e-7 and was wrong.

Read to 100k this column looks like a vanishing grid artifact; read to 250k it
looks converged.  Both readings were made, in that order, and both were wrong."""

TEST_OLD = '''@pytest.mark.slow
def test_c3_converges_and_vindicates_the_appendix_figure():
    """The c=3 column is a real box-truncation quantity, and it converges.

    WRONG ANSWER REJECTED: "c=3 is under the quadrature floor, so its value is
    a grid artifact and the appendix's ~1e-7 is only an upper bound."  That is
    what the ladder looks like if you stop at 100k, and both a reviewer and the
    author of this file concluded it before extending the ladder.  It is false:
    the column PLATEAUS at the next point (1.5592e-07 -> 1.5676e-07, a step of
    0.99x), so ~1.6e-7 is the converged truncation error and the appendix
    figure is right.

    Asserting convergence rather than a value keeps this honest at whatever
    resolution it runs: a genuine artifact could not hold still.
    """
    nmax, box = 8, 3.0 * 8 ** 2
    global NPTS
    saved = NPTS
    try:
        vals = []
        for npts in (100000, 250000):
            NPTS = npts
            vals.append(_rel_err(nmax, box))
    finally:
        NPTS = saved
    lo, hi = vals
    assert abs(hi - lo) / lo < 0.10, (
        f"c=3 must PLATEAU, not keep falling: {lo:.4e} -> {hi:.4e}")
    assert hi < 1e-6, f"the converged c=3 error should be ~1e-7, got {hi:.4e}"'''

TEST_NEW = '''@pytest.mark.slow
def test_c3_is_a_genuine_truncation_of_order_1e_minus_7():
    """The c=3 box error is real and of order 1e-7 -- not a vanishing artifact.

    WRONG ANSWER REJECTED (1): "c=3 sits under the quadrature floor, so it
    falls without limit with resolution and the appendix's ~1e-7 is only an
    upper bound."  Measured, it stops falling and then RISES:  1.559e-7 at
    100k, 1.568e-7 at 250k, 1.931e-7 at 400k, 2.060e-7 at 600k.  A grid
    artifact does not rise.

    WRONG ANSWER REJECTED (2): "it plateaus at 1.5676e-07, so that is the
    converged value."  That is what THIS FILE asserted on 2026-09-13 and it is
    false -- the 250k->400k step moves 23%.  The flat-looking step at
    100k->250k is a crossing: the quantity is relative to a c=5 box whose own
    truncation (7.7e-8 absolute) is comparable to c=3's, so the two partially
    cancel there.  The previous guard pinned exactly that pair -- the only
    adjacent pair in the ladder under a 10% window -- and a reviewer fired it
    by planting MORE resolution, a change with no physics in it.

    So this asserts what is robust across discretizations, which is also what
    App. A actually claims: the ORDER, and the non-collapse.  It deliberately
    asserts no plateau and no third digit, because the uniform and graded
    meshes disagree by ~1.7x on the converged value (1.29e-7 vs 2.16e-7).
    """
    nmax, box = 8, 3.0 * 8 ** 2
    global NPTS
    saved = NPTS
    try:
        vals = []
        for npts in (100000, 600000):
            NPTS = npts
            vals.append(_rel_err(nmax, box))
    finally:
        NPTS = saved
    lo, hi = vals
    for v, n in zip(vals, (100000, 600000)):
        assert 1e-8 < v < 1e-6, (
            f"c=3 must be a genuine ~1e-7 truncation at npts={n}, got {v:.4e} "
            f"-- below 1e-8 would mean a vanishing grid artifact, above 1e-6 "
            f"would contradict App. A")
    assert hi > 0.9 * lo, (
        f"c=3 must NOT collapse with resolution (that is the artifact reading "
        f"this rejects): {lo:.4e} at 100k -> {hi:.4e} at 600k")'''

MTX_OLD = ("**Measurement note:** the first draft's table was GRID-limited and its "
           "guards failed; redone at converged resolution the paper's own ~1e-7 figure "
           "is **VINDICATED** — c=3 reads 2.6e-5 / 2.1e-6 / 1.6e-7 / 1.57e-7 at "
           "12k / 40k / 100k / 250k points, i.e. it falls 12-13x per step and then "
           "PLATEAUS (0.99x), so it is a genuine box-truncation quantity. "
           "**Recorded because it cost two wrong conclusions:** read to 100k only, the "
           "column looks unconverged, and on 2026-09-13 both a code reviewer and the "
           "PM concluded from the three-point version that c=3 was a quadrature floor "
           "and re-tiered this row to \"consistent but not confirmed\". Extending the "
           "ladder one point refuted that; the reviewer withdrew its own finding and "
           "the row is restored. Guard thresholds were left resolution-robust (they "
           "hold at 12k, 40k and 100k, not only above the hard-coded NPTS), which is "
           "the one improvement that survived.")

MTX_NEW = ("**Measurement note (THIRD reading, 2026-09-13):** App. A's ~1e-7 for c=3 is "
           "**VINDICATED as an ORDER**, and the third digit is **not determined**. "
           "Ladder: 2.57e-5 / 2.11e-6 / 1.559e-7 / 1.568e-7 / 1.931e-7 / 2.060e-7 at "
           "12k / 40k / 100k / 250k / 400k / 600k. **c=1 and c=2 plateau; c=3 does "
           "not** — the flat-looking 1.005x step at 250k is a CROSSING, because this "
           "quantity is relative to a c=5 box whose own truncation (7.69e-8 absolute, "
           "vs c=3's 1.29e-7) is comparable and partially cancels there. Absolute "
           "against a c=12 reference: uniform mesh 1.29e-7, graded mesh 2.16e-7 — the "
           "two discretizations differ by ~1.7x, so no converged value is quoted. "
           "**Recorded because this ladder was misread twice, in opposite directions, "
           "by two different parties:** first a code reviewer and the PM read it to "
           "100k and called c=3 a quadrature floor (re-tiering a correct claim); then "
           "the PM read it to 250k and called it a plateau, publishing a table whose "
           "own last row rose 31% beneath the word PLATEAU. A reviewer sampling 400k — "
           "a point nobody had taken — broke it, and an independent PM ladder "
           "reproduced that digit for digit. The guard now asserts the order and the "
           "non-collapse, never a plateau at a hard-coded pair.")


def main() -> int:
    loaded: dict[str, str] = {}
    applied, missed = [], []
    for path, name, marker, old, new in [
        (T, "docstring", "it is NOT", DOC_OLD, DOC_NEW),
        (T, "guard", "test_c3_is_a_genuine_truncation", TEST_OLD, TEST_NEW),
        (MTX, "matrix-note", "THIRD reading, 2026-09-13", MTX_OLD, MTX_NEW),
    ]:
        if path not in loaded:
            with open(path, encoding="utf-8") as fh:
                loaded[path] = fh.read()
        t = loaded[path]
        if marker in t:
            print(f"  skip  {name} (already applied)")
            continue
        if t.count(old) != 1:
            missed.append((name, t.count(old)))
            continue
        loaded[path] = t.replace(old, new)
        applied.append(name)
    for path, t in loaded.items():
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(t)
    for n in applied:
        print(f"  ok    {n}")
    for n, c in missed:
        print(f"  MISS  {n}: anchor count={c}")
    print(f"applied {len(applied)}, MISSED {len(missed)}")
    return 3 if missed else 0


if __name__ == "__main__":
    sys.exit(main())
