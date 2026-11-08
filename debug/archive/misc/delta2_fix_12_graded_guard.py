"""DELTA #2 -- the guard leg that actually reaches the converged c=3 value.

Written as a separate pass from the docstring correction it completes, and
reviewed by asking what wrong answer it rejects rather than whether it passes.

WRONG ANSWER REJECTED, named: "the converged c=3 value cannot be pinned,
because the uniform and graded meshes disagree by ~1.7x (1.29e-7 vs 2.16e-7)."
That is what this repo asserted for part of 2026-09-13.  They agree to 0.43%.

The existing `test_c3_is_a_genuine_truncation_of_order_1e_minus_7` runs at 100k
and 600k on the UNIFORM mesh, where the column is still climbing, so it can
only assert the order.  The graded mesh `r = R t^2` resolves the near-origin
region where the integrand lives and converges at 240k -- about 6x cheaper than
the ~1.5M the uniform mesh needs -- so it CAN assert the value.  That is the
gap this closes.

Idempotent.
"""
from __future__ import annotations

import sys

T = "tests/test_paper60_box_rule.py"

NEW_TEST = '''

@pytest.mark.slow
def test_c3_converged_value_on_the_graded_mesh():
    """The converged c=3 box truncation is 2.15e-07, and the meshes agree.

    WRONG ANSWER REJECTED: "the converged value cannot be pinned, because the
    uniform and graded meshes disagree by ~1.7x (1.29e-07 vs 2.16e-07)."  This
    repo asserted that for part of 2026-09-13 and it is false -- they agree to
    0.43%.  The 1.29e-07 came from referencing c=3 against a c=12 box at the
    SAME point count, i.e. with 4x coarser spacing than the quantity it was
    resolving; that ratio is not even monotone in npts (2.86e-07 / 1.29e-07 /
    2.02e-07 at 250k / 600k / 1.5M).

    The graded mesh r = R t^2 resolves the near-origin region where the
    integrand lives, so it converges at 240k points where the uniform mesh
    needs ~1.5M.  Measured on the graded mesh: 2.045e-07 / 2.134e-07 /
    2.15582e-07 at 60k / 120k / 240k, step 1.010.  The uniform ladder reaches
    2.14665e-07 at 1.5M.  Cross-check: the c=5 ABSOLUTE value reads
    17.2714838375 graded against 17.2714838504 uniform, agreeing to 7.5e-10.

    Asserts BOTH convergence (the last step must be small -- an unconverged
    quantity cannot hold still) and the value, because both are established.
    The 3% tolerance admits the measured mesh-to-mesh spread with ~7x margin.
    """
    nmax = 8
    saved = (SS.r, SS.dr, SS.r2, SS.R_MAX, SS.N_GRID,
             SS._ctrap_fwd, SS._ctrap_rev)
    try:
        def rel(npts: int) -> float:
            vals = []
            for c in (3.0, 5.0):
                SV.set_grid(c * nmax ** 2, npts, kind="grade")
                cfgs = SS.build_configs(SV.family(nmax, 1))
                M = SS.build_M(cfgs, Z=Z_HE)
                D = np.diag([Z_HE * x.Rnu for x in cfgs])
                Tp = M - D
                vals.append(float(np.abs(Tp).sum()
                                  - np.abs(np.diag(Tp)).sum()))
            return abs(vals[0] - vals[1]) / abs(vals[1])

        mid, top = rel(120_000), rel(240_000)
    finally:
        (SS.r, SS.dr, SS.r2, SS.R_MAX, SS.N_GRID,
         SS._ctrap_fwd, SS._ctrap_rev) = saved
        SS._GAUNT_CACHE.clear()
        SS.reset_caches()

    assert top / mid < 1.05, (
        f"the graded mesh must have CONVERGED by 240k, not still be climbing: "
        f"{mid:.5e} -> {top:.5e} is a step of {top / mid:.3f}x")
    assert abs(top - 2.155e-07) / 2.155e-07 < 0.03, (
        f"converged c=3 truncation should be 2.15e-07 (uniform agrees at "
        f"2.14665e-07); got {top:.5e}")
'''


def main() -> int:
    with open(T, encoding="utf-8") as fh:
        t = fh.read()
    if "test_c3_converged_value_on_the_graded_mesh" in t and "def test_c3_converged_value" in t:
        print("already applied")
        return 0
    anchor = "\n\n@pytest.mark.slow\ndef test_c3_is_a_genuine_truncation_of_order_1e_minus_7():"
    if t.count(anchor) != 1:
        print(f"  MISS anchor count={t.count(anchor)}")
        return 3
    with open(T, "w", encoding="utf-8") as fh:
        fh.write(t.replace(anchor, NEW_TEST + anchor, 1))
    print("  ok    test_c3_converged_value_on_the_graded_mesh added")
    return 0


if __name__ == "__main__":
    sys.exit(main())
