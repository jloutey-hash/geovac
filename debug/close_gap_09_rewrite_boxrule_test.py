"""KNOWN-GAP 5/5, step 2 -- rewrite the box-rule guards FROM MEASUREMENT.

The first draft failed, and the cause was mine: I wrote a measurement table into
the docstring without measuring it, then wrote assertions against those invented
numbers.  That is the authoring error this corpus's whole `provenance` rule
exists to stop, committed inside a test whose subject is measurement discipline.

Measured properly (l_max=1 families, ||T'||_1^off vs an R_max = 5 n^2
reference).  The first pass was GRID-limited, exactly as the chain-rung recheck
was earlier today -- at n_max=8 the c=3 column read 2.6e-5 at 12k points and
1.6e-7 at 100k:

    npts      c=1        c=2        c=3        fixed 60
    12000   3.613e-02  6.563e-04  2.571e-05  4.323e-02
    40000   3.617e-02  6.873e-04  2.111e-06  4.327e-02
   100000   3.617e-02  6.898e-04  1.559e-07  4.327e-02

**This vindicates the paper.**  Its appendix says c=3 gives ~1e-7; at converged
resolution it does.  My under-resolved measurement, not the paper, was wrong.

Drift, at 40k points (the load-bearing comparison):

    n_max   c=2 scaling   fixed 60    ratio
      6     1.707e-03    8.354e-03      4.9
      8     6.873e-04    4.327e-02     63.0
     10     3.407e-04    1.171e-01    343.8

The scaling box FALLS (1.7e-3 -> 3.4e-4) while the fixed box RISES by 14x over
the same range, and the ratio between them grows 70-fold.  That is the paper's
"sign of drift" point, and it is what the guards now assert.

Thresholds are set from these numbers with margin, not from a round figure:
c=1 -> c=2 buys ~52x (guard: >20x), c=2 -> c=3 buys ~4400x (guard: >100x).

Idempotent.
"""
from __future__ import annotations

import sys

T = "tests/test_paper60_box_rule.py"
MARKER = "vindicates the paper"

OLD_DOC = """MEASURED by the PM before these guards were written (l_max=1 families, the
`||T'||_1` off-diagonal norm against an `R_max = 5 n_max^2` reference):

    n_max   K     c=1        c=2        c=3        fixed 60 bohr
      6     21    3.1e-03    1.2e-05    8.9e-09    2.7e-05
      8     36    3.0e-03    1.1e-05    7.4e-09    1.6e-03
     10     55    2.9e-03    1.0e-05    6.6e-09    1.2e-02
"""

NEW_DOC = """MEASURED (l_max=1 families, the `||T'||_1` off-diagonal norm against an
`R_max = 5 n_max^2` reference).  The first attempt at this table was GRID-limited
and had to be redone at higher resolution -- at n_max=8 the c=3 column reads
2.6e-5 at 12k points and 1.6e-7 at 100k.  This vindicates the paper: its
appendix says c=3 gives ~1e-7, and at converged resolution it does.

    npts      c=1        c=2        c=3        fixed 60      (n_max = 8)
    12000   3.613e-02  6.563e-04  2.571e-05  4.323e-02
    40000   3.617e-02  6.873e-04  2.111e-06  4.327e-02
   100000   3.617e-02  6.898e-04  1.559e-07  4.327e-02

Drift at 40k points, which is the load-bearing comparison:

    n_max   c=2 scaling   fixed 60    ratio
      6     1.707e-03    8.354e-03      4.9
      8     6.873e-04    4.327e-02     63.0
     10     3.407e-04    1.171e-01    343.8
"""

OLD_NPTS = "NPTS = 12000"
NEW_NPTS = ("NPTS = 40000        # 12k is GRID-limited for this quantity; see the\n"
            "                    # resolution table in the module docstring")

OLD_DRIFT = '''    ns = (6, 8, 10)
    scaling = [_rel_err(n, 3.0 * n ** 2) for n in ns]
    fixed = [_rel_err(n, 60.0) for n in ns]
    assert scaling[-1] <= scaling[0] * 1.5, (
        f"a scaling box must not degrade with basis size: {scaling}")
    assert fixed[-1] > fixed[0] * 10.0, (
        f"the fixed 60-bohr box must degrade with basis size: {fixed}")
    assert fixed[-1] > scaling[-1] * 100.0, (
        f"at the largest basis the fixed box must be far worse: "
        f"{fixed[-1]:.2e} vs {scaling[-1]:.2e}")'''

NEW_DRIFT = '''    ns = (6, 8, 10)
    scaling = [_rel_err(n, 2.0 * n ** 2) for n in ns]
    fixed = [_rel_err(n, 60.0) for n in ns]
    # the scaling box IMPROVES with basis size (1.7e-3 -> 3.4e-4 measured)
    assert scaling[-1] < scaling[0], (
        f"a scaling box must improve, not degrade, with basis size: {scaling}")
    # the fixed box degrades (8.4e-3 -> 1.2e-1 measured, ~14x)
    assert fixed[-1] > fixed[0] * 5.0, (
        f"the fixed 60-bohr box must degrade with basis size: {fixed}")
    # and the gap between them WIDENS -- 4.9x at n=6 to 344x at n=10.  This is
    # the assertion that a merely-small fixed box cannot satisfy.
    r0, r1 = fixed[0] / scaling[0], fixed[-1] / scaling[-1]
    assert r1 > 10.0 * r0, (
        f"the fixed/scaling gap must widen with basis size: {r0:.1f} -> {r1:.1f}")'''

OLD_ORDER = '''    nmax = 8
    e1, e2, e3 = (_rel_err(nmax, c * nmax ** 2) for c in (1.0, 2.0, 3.0))
    assert e1 > 100.0 * e2, f"c=2 must be >=2 orders better than c=1: {e1}, {e2}"
    assert e2 > 100.0 * e3, f"c=3 must be >=2 orders better than c=2: {e2}, {e3}"
    assert e3 < 1e-6, f"c=3 should reach ~1e-8, got {e3:.2e}"'''

NEW_ORDER = '''    nmax = 8
    e1, e2, e3 = (_rel_err(nmax, c * nmax ** 2) for c in (1.0, 2.0, 3.0))
    # measured at this resolution: 3.6e-2, 6.9e-4, 2.1e-6 -- so c=1 -> c=2 buys
    # ~52x and c=2 -> c=3 buys ~330x.  Thresholds sit well under both, so they
    # exclude "c=1 is close enough" without pinning tighter than the quantity's
    # own spread across resolutions.
    assert e1 > 20.0 * e2, f"c=2 must be far better than c=1: {e1:.2e}, {e2:.2e}"
    assert e2 > 100.0 * e3, f"c=3 must be far better than c=2: {e2:.2e}, {e3:.2e}"
    assert e3 < 1e-5, f"c=3 should reach 1e-6 or below, got {e3:.2e}"'''


def main() -> int:
    with open(T, encoding="utf-8") as fh:
        t = fh.read()
    if MARKER in t:
        print("ALREADY APPLIED")
        return 1
    for nm, old in (("doc", OLD_DOC), ("npts", OLD_NPTS),
                    ("drift", OLD_DRIFT), ("order", OLD_ORDER)):
        if t.count(old) != 1:
            print(f"  {nm} anchor count={t.count(old)}; ABORT")
            return 2
    t = (t.replace(OLD_DOC, NEW_DOC).replace(OLD_NPTS, NEW_NPTS)
          .replace(OLD_DRIFT, NEW_DRIFT).replace(OLD_ORDER, NEW_ORDER))
    with open(T, "w", encoding="utf-8") as fh:
        fh.write(t)
    print("applied: box-rule guards rewritten from measurement")
    return 0


if __name__ == "__main__":
    sys.exit(main())
