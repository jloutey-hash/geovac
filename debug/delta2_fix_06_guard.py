"""DELTA #2 -- the guard for the CORRECTED water-control claim.

Written as a SEPARATE pass from the fix (CLAUDE.md Sec.9 guard-writing rule),
and reviewed by asking what wrong answer it would accept rather than whether it
passes.

THE WRONG ANSWER IT REJECTS, named:  "the uniform band control halves the
growth exponent, 1.96 -> 0.98, so banding alone is most of the lever and the
rotation is a refinement."  That is what paper 60 asserted until 2026-09-13.
It is false twice over -- the uniform control COMMUTES with the rotation, so it
cannot speak to the rotation at all, and its exponent is 1.950 against the raw
1.967, i.e. banding alone does essentially nothing.

The existing `test_water_needs_the_null_direction_rotation` asserts the
selective/unrotated blow-up (>4.5x per doubling, >20x raw) and remains the
load-bearing discrimination.  What it does NOT assert, and what the paper now
prints, is the UNIFORM column's exponent -- the number that was wrong.  A guard
that does not pin the number that was wrong is not a guard for this defect.

Idempotent.
"""
from __future__ import annotations

import sys

T = "tests/test_paper60_preconditioner.py"

NEW_TEST = '''

def test_uniform_banding_alone_is_not_the_lever():
    """The exponent the paper got wrong, pinned.

    WRONG ANSWER REJECTED: "uniform banding halves the growth exponent
    (1.96 -> 0.98), so it is most of the lever."  Paper 60 printed exactly that
    until 2026-09-13; /qa paper_60 DELTA #2 measured it.  The uniform column
    runs N^1.950 against the raw N^1.967 -- the two are within 1% of each other
    and banding alone buys essentially nothing.  The 0.98 was arithmetic on a
    two-point range read per-doubling instead of per-decade.

    This is the number, not the mechanism: the mechanism (that the uniform
    control commutes with the rotation and is therefore blind to it) is pinned
    by `test_uniform_band_preconditioner_commutes_with_the_rotation` above, and
    the discrimination is carried by
    `test_water_needs_the_null_direction_rotation` below.
    """
    ns = (12, 24, 48, 96)
    raw, uni = [], []
    for n in ns:
        A = _water_A1(n)
        P = np.zeros_like(A)
        P[:n, :n] = tridiag(n)
        P[n:, n:] = tridiag(n)
        P_is = inv_sqrt(P)
        raw.append(np.linalg.cond(A))
        uni.append(np.linalg.cond(P_is @ A @ P_is))

    ln = np.log(np.array(ns, float))
    e_raw = float(np.polyfit(ln, np.log(np.array(raw)), 1)[0])
    e_uni = float(np.polyfit(ln, np.log(np.array(uni)), 1)[0])

    assert 1.85 < e_raw < 2.05, f"raw exponent should be ~1.97, got {e_raw:.3f}"
    assert 1.85 < e_uni < 2.05, (
        f"uniform-banding exponent should be ~1.95 -- NOT halved to ~0.98, which "
        f"is what the paper claimed before 2026-09-13 -- got {e_uni:.3f}")
    assert abs(e_uni - e_raw) < 0.10, (
        f"uniform banding must NOT materially change the exponent: raw "
        f"{e_raw:.3f} vs uniform {e_uni:.3f}")
'''


def main() -> int:
    with open(T, encoding="utf-8") as fh:
        t = fh.read()
    if "test_uniform_banding_alone_is_not_the_lever" in t:
        print("already applied")
        return 0
    anchor = "\n\ndef test_water_needs_the_null_direction_rotation():"
    if t.count(anchor) != 1:
        print(f"  MISS anchor count={t.count(anchor)}")
        return 3
    with open(T, "w", encoding="utf-8") as fh:
        fh.write(t.replace(anchor, NEW_TEST + anchor, 1))
    print("  ok    test_uniform_banding_alone_is_not_the_lever added")
    return 0


if __name__ == "__main__":
    sys.exit(main())
