"""DELTA remediation -- the three code findings, all on my own box-rule work.

M2 + M3 are one error with two faces, and it is the SAME error I was caught by
twice yesterday: concluding from a measurement without checking whether the
measurement had converged.

MEASURED (n_max=8, rel. err of ||T'||_1^off vs a c=5 reference at matched npts):

    npts      c=1         c=2         c=3        e2/e3
    12000   3.6130e-02  6.5629e-04  2.5709e-05     25.5
    40000   3.6166e-02  6.8727e-04  2.1107e-06    325.6
   100000   3.6169e-02  6.8984e-04  1.5592e-07   4424.3

The c=1 and c=2 columns PLATEAU (third digit stable).  The c=3 column does NOT:
it falls 12.2x then 13.5x per resolution step, with no sign of settling.

What that means, stated correctly: at c=3 the box-truncation error is BELOW the
quadrature floor at every resolution tested, so what the c=3 column measures is
the grid mismatch between the box=192 and box=320 discretizations, not box
truncation.  The paper's appendix figure of ~1e-7 is therefore CONSISTENT with
the data and NOT CONTRADICTED -- but it is not confirmed either, and I wrote
"redone at converged resolution the paper's own ~1e-7 figure is VINDICATED" in
both the claim matrix and the test docstring.  That is an overclaim.  Corrected
to a bound.

M2: two of the three assertions (`e2 > 100*e3`, `e3 < 1e-5`) are FALSE at 12k
and pass only above the hard-coded NPTS = 40000 -- i.e. they pin a resolution
choice rather than the rule.  Replaced with legs that hold at 12k, 40k and 100k,
plus the one that actually carries the content: c=1 and c=2 are stable across
resolutions while c=3 falls, which IS the evidence that c=3 sits under the floor.

M1: `tests/test_paper60_full_shell_family.py` is cited by ZERO claim-matrix rows.
The FULL run's finding was "no test, no registry key and no C17 family"; the
remediation closed the test leg and left the registration leg open, so deleting
the file would leave every gate green.  C22 checks rows->tests and never
tests->rows, so this backing was invisible.

Idempotent.
"""
from __future__ import annotations

import sys

T = "tests/test_paper60_box_rule.py"
MTX = "docs/claim_test_matrix.md"

OLD_DOC = """and had to be redone at higher resolution -- at n_max=8 the c=3 column reads
2.6e-5 at 12k points and 1.6e-7 at 100k.  This vindicates the paper: its
appendix says c=3 gives ~1e-7, and at converged resolution it does."""

NEW_DOC = """and had to be redone at higher resolution.  Read the table carefully, because
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

OLD_TEST = '''@pytest.mark.slow
def test_the_coefficient_ordering_pins_c_equals_three():
    """(b) c = 1 does not suffice; the rule's coefficient cannot drift down.

    Each increment of c must buy at least two orders of magnitude, so
    'c = 1 is close enough' is excluded by the ordering rather than by a
    tolerance someone can loosen.
    """
    nmax = 8
    e1, e2, e3 = (_rel_err(nmax, c * nmax ** 2) for c in (1.0, 2.0, 3.0))
    # measured at this resolution: 3.6e-2, 6.9e-4, 2.1e-6 -- so c=1 -> c=2 buys
    # ~52x and c=2 -> c=3 buys ~330x.  Thresholds sit well under both, so they
    # exclude "c=1 is close enough" without pinning tighter than the quantity's
    # own spread across resolutions.
    assert e1 > 20.0 * e2, f"c=2 must be far better than c=1: {e1:.2e}, {e2:.2e}"
    assert e2 > 100.0 * e3, f"c=3 must be far better than c=2: {e2:.2e}, {e3:.2e}"
    assert e3 < 1e-5, f"c=3 should reach 1e-6 or below, got {e3:.2e}"'''

NEW_TEST = '''@pytest.mark.slow
def test_the_coefficient_ordering_pins_c_equals_three():
    """(b) c = 1 does not suffice, asserted RESOLUTION-ROBUSTLY.

    An earlier version asserted `e2 > 100*e3` and `e3 < 1e-5`.  Both are FALSE
    at 12k points (25.5 and 2.6e-5) and pass only above the hard-coded NPTS, so
    they pinned a resolution choice rather than the rule -- planting
    NPTS 40000 -> 12000, a change with no physics in it, made them fire.
    Every leg below holds at 12k, 40k AND 100k.
    """
    nmax = 8
    e1, e2, e3 = (_rel_err(nmax, c * nmax ** 2) for c in (1.0, 2.0, 3.0))
    assert e1 > 20.0 * e2, f"c=2 must be far better than c=1: {e1:.2e}, {e2:.2e}"
    assert e2 > 5.0 * e3, f"c=3 must be better than c=2: {e2:.2e}, {e3:.2e}"
    assert e3 < 1e-4, f"c=3 must be well under the c=2 level, got {e3:.2e}"


@pytest.mark.slow
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

MTX_OLD = ("**Measurement note:** the "
           "first draft's table was GRID-limited (c=3 read 2.6e-5 at 12k points, 1.6e-7 "
           "at 100k) and its guards failed; redone at converged resolution the paper's "
           "own ~1e-7 figure is VINDICATED.")
MTX_NEW = ("**Measurement note (corrected 2026-09-13):** the c=1 and c=2 columns "
           "PLATEAU across resolutions and measure box truncation; the c=3 column does "
           "NOT (2.6e-5 at 12k, 2.1e-6 at 40k, 1.6e-7 at 100k -- falling 12-13x per "
           "step), so at c=3 the box error is already below the quadrature floor and "
           "that column measures grid mismatch. The paper's ~1e-7 is therefore "
           "CONSISTENT and not contradicted, but **not confirmed** -- an earlier "
           "version of this row said \"VINDICATED at converged resolution\", which "
           "overstated it. The guards now assert only resolution-robust legs, plus the "
           "stable-vs-falling asymmetry that is the actual evidence.")

MTX_ROW = (
    "| 60 | sec:atomic (abstract L35) — **the full-shell growth rule**, which the "
    "paper calls \"the substantive finding\" of its atomic section: the sublinearity "
    "belongs to the BASIS-GROWTH RULE, not the isoenergetic construction. Grown as "
    "full hydrogenic shells (`l_max = n_max - 1`) the same construction gives a "
    "SUPERLINEAR `||T'||_1 ~ K^1.07` while the total exponent rises 0.849 -> 0.911 "
    "and never reaches 1; `0.867` is the K=35->56 rung, NOT an interior value of the "
    "K=56–220 window (first interior slope 0.879, global fit 0.893) | "
    "`tests/test_paper60_full_shell_family.py```::test_full_shell_tprime_block_is_"
    "superlinear`` + ``::test_full_shell_total_rises_but_never_reaches_one`` + "
    "``::test_full_shell_endpoints_name_the_right_window`` (all slow) | tracked "
    "`geovac/sturmian_{secular,variational}.py` | **NEW 2026-09-13** | BACKED-SOUND. "
    "**Registration gap closed:** the FULL run found this claim had \"no test, no "
    "registry key and no C17 family\"; the 2026-09-12 remediation added the test and "
    "left it cited by NO row, so deleting the file would have left every gate green "
    "(C22 checks rows->tests, never tests->rows). Fire-tested: switching to the "
    "fixed-`l_max` family FIRES the superlinear guard; claiming 0.867 sits inside the "
    "K=56–220 window FIRES the endpoint guard. **Guard weaker than the prose (owed):** "
    "the superlinear leg asserts `fit > 1.02` where the measured value is 1.0727, and "
    "the 0.893 global fit is asserted nowhere. Independently reproduced two cutoffs "
    "past the window (K=286 -> 0.9209, K=364 -> 0.9304, increments decaying "
    "monotonically), so the claim holds beyond its stated range. rests on: "
    "eq:sublinear (the fixed-`l_max` family this one is the counter-case to) |\n")
MTX_ANCHOR = "| 60 | App. A `eq:boxrule` —"


def main() -> int:
    with open(T, encoding="utf-8") as fh:
        t = fh.read()
    with open(MTX, encoding="utf-8") as fh:
        m = fh.read()
    n = 0
    if "below the quadrature floor" in t:
        print("  skip test (already applied)")
    else:
        for nm, s in (("doc", OLD_DOC), ("test", OLD_TEST)):
            if t.count(s) != 1:
                print(f"  {nm} anchor count={t.count(s)}; ABORT")
                return 2
        t = t.replace(OLD_DOC, NEW_DOC).replace(OLD_TEST, NEW_TEST)
        with open(T, "w", encoding="utf-8") as fh:
            fh.write(t)
        n += 1
        print("  ok   box-rule docstring + resolution-robust guards")
    if "corrected 2026-09-13" in m:
        print("  skip matrix note (already applied)")
    elif m.count(MTX_OLD) == 1:
        m = m.replace(MTX_OLD, MTX_NEW)
        n += 1
        print("  ok   matrix: VINDICATED -> bound")
    else:
        print(f"  MISS matrix note: count={m.count(MTX_OLD)}")
        return 3
    if "the full-shell growth rule" in m:
        print("  skip full-shell row (already applied)")
    elif m.count(MTX_ANCHOR) == 1:
        m = m.replace(MTX_ANCHOR, MTX_ROW + MTX_ANCHOR)
        n += 1
        print("  ok   matrix: full-shell row added")
    else:
        print(f"  MISS full-shell anchor: count={m.count(MTX_ANCHOR)}")
        return 4
    with open(MTX, "w", encoding="utf-8") as fh:
        fh.write(m)
    print(f"applied {n}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
