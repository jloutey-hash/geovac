"""REVERT delta_fix_05's M3 "correction" -- it corrected a CORRECT claim.

Sequence, recorded because the lesson is mine:

  1. The code reviewer filed M3: the c=3 box-rule column "never plateaus", so
     the claim-matrix's "redone at converged resolution ... VINDICATED" was
     false and the figure was only a bound.
  2. I verified its THREE tabulated points (12k/40k/100k), saw the 12-13x
     per-step fall, agreed, and rewrote the claim matrix and the test docstring
     to a weaker "consistent but not confirmed" reading.
  3. The reviewer then extended its own ladder to 250k and 600k, found a
     PLATEAU, and withdrew M3 itself.
  4. I measured it independently: c=3 = 1.5592e-07 at 100k and 1.5676e-07 at
     250k -- a step of 0.99x.  FLAT.

So the column is a genuine, convergent box-truncation quantity, the paper's
appendix figure of ~1e-7 is accurate, and the original "VINDICATED" wording was
right.  My error was accepting a measurement-based finding on the reviewer's
sample instead of extending the measurement myself -- the exact discipline I
had been applying everywhere else today, and the third grid-convergence trap of
this arc.

What is RESTORED: the vindication, now stated with the full ladder, which is
stronger evidence than the original had.

What is KEPT from delta_fix_05: the resolution-robust thresholds.  Those remain
a real improvement (they hold at 12k, 40k and 100k rather than only above the
hard-coded NPTS) and the reviewer's demoted finding N6 stands as a NIT.

What is REMOVED: `test_c3_sits_below_the_quadrature_floor`.  It passes, but its
docstring asserts a FALSE mechanism -- that c=3 is under the grid floor.  A
guard that passes while naming the wrong reason is worse than none.  Replaced
by one that asserts the true and stronger fact: the column CONVERGES.

Idempotent.
"""
from __future__ import annotations

import sys

T = "tests/test_paper60_box_rule.py"
MTX = "docs/claim_test_matrix.md"

DOC_OLD = """and had to be redone at higher resolution.  Read the table carefully, because
the obvious conclusion from it is wrong:

  * c=1 and c=2 PLATEAU across resolutions (third digit stable), so those
    columns measure box truncation;
  * c=3 does NOT plateau -- it falls 12.2x then 13.5x per step.  At c=3 the box
    error is already BELOW the quadrature floor, so that column measures the
    grid mismatch between the two discretizations, not truncation.

So the paper's appendix figure of ~1e-7 at c=3 is CONSISTENT with this data and
is not contradicted -- but it is NOT confirmed by it, and an earlier version of
this docstring said the figure was "vindicated at converged resolution".  It is
a bound: the c=3 box error is below 1.6e-7 and below the floor at every grid
tested."""

DOC_NEW = """and had to be redone at higher resolution.  The full ladder, because three
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

TEST_OLD = '''@pytest.mark.slow
def test_c3_sits_below_the_quadrature_floor():
    """Why the c=3 number cannot be quoted as a box-truncation error.

    WRONG ANSWER REJECTED: "the c=3 column measures box truncation, so its
    value is the rule's residual at c=3."  It does not.  c=1 and c=2 are
    resolution-STABLE while c=3 FALLS with resolution -- the signature of a
    quantity already under the grid floor.  That asymmetry is the evidence,
    and asserting it is what keeps the paper's ~1e-7 an upper bound rather
    than a confirmed value.
    """
    nmax, box = 8, lambda c: c * nmax ** 2
    global NPTS
    saved = NPTS
    try:
        vals = {}
        for npts in (12000, 40000):
            NPTS = npts
            vals[npts] = [_rel_err(nmax, box(c)) for c in (1.0, 2.0, 3.0)]
    finally:
        NPTS = saved
    lo, hi = vals[12000], vals[40000]
    for i, name in ((0, "c=1"), (1, "c=2")):
        assert abs(hi[i] - lo[i]) / lo[i] < 0.10, (
            f"{name} must be resolution-stable: {lo[i]:.4e} -> {hi[i]:.4e}")
    assert lo[2] > 5.0 * hi[2], (
        f"c=3 must FALL with resolution (it is under the floor): "
        f"{lo[2]:.4e} -> {hi[2]:.4e}")'''

TEST_NEW = '''@pytest.mark.slow
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

MTX_OLD = ("**Measurement note (corrected 2026-09-13):** the c=1 and c=2 columns "
           "PLATEAU across resolutions and measure box truncation; the c=3 column does "
           "NOT (2.6e-5 at 12k, 2.1e-6 at 40k, 1.6e-7 at 100k -- falling 12-13x per "
           "step), so at c=3 the box error is already below the quadrature floor and "
           "that column measures grid mismatch. The paper's ~1e-7 is therefore "
           "CONSISTENT and not contradicted, but **not confirmed** -- an earlier "
           "version of this row said \"VINDICATED at converged resolution\", which "
           "overstated it. The guards now assert only resolution-robust legs, plus the "
           "stable-vs-falling asymmetry that is the actual evidence.")

MTX_NEW = ("**Measurement note:** the first draft's table was GRID-limited and its "
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


def main() -> int:
    with open(T, encoding="utf-8") as fh:
        t = fh.read()
    with open(MTX, encoding="utf-8") as fh:
        m = fh.read()
    if "vindicates_the_appendix_figure" in t and "VINDICATED**" in m:
        print("ALREADY APPLIED")
        return 1
    for nm, s, hay in (("doc", DOC_OLD, t), ("test", TEST_OLD, t),
                       ("matrix", MTX_OLD, m)):
        if hay.count(s) != 1:
            print(f"  {nm} anchor count={hay.count(s)}; ABORT")
            return 2
    t = t.replace(DOC_OLD, DOC_NEW).replace(TEST_OLD, TEST_NEW)
    m = m.replace(MTX_OLD, MTX_NEW)
    with open(T, "w", encoding="utf-8") as fh:
        fh.write(t)
    with open(MTX, "w", encoding="utf-8") as fh:
        fh.write(m)
    print("applied: vindication restored with the full ladder; "
          "floor-guard replaced by a convergence guard")
    return 0


if __name__ == "__main__":
    sys.exit(main())
