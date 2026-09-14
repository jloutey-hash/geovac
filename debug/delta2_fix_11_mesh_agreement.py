"""DELTA #2 / M3 -- the FOURTH wrong mechanism on one claim, and it is mine.

The c=3 box-truncation number has now carried four different explanations, each
offered as a correction of the last, and the first three were wrong:

  1. original          "VINDICATED at converged resolution"      (3 points)
  2. delta_fix_05      "c=3 sits below the quadrature floor"     (3 points, opposite)
  3. delta_fix_06      "it PLATEAUS at 1.5676e-07"               (4 points)
  4. delta2_fix_09     "the meshes disagree ~1.7x, so the third
                        digit is not determined"                 (FALSE -- they agree)

The paper's App. A claim ("c=3 gives ~1e-7") was correct every single time.
Every failure was in the EVIDENCE offered for it.

WHY (4) WAS WRONG, measured here and confirmed independently.  My "absolute"
route set `ref = tprime_off(12 * n^2, 600000)` -- a 768-bohr box at the SAME
point count as the 192-bohr box it was the reference for, so its spacing was
4x coarser than the quantity it was supposed to resolve.  Run up the ladder it
is not even monotone:

    rel(c3 vs uniform c12, both at npts)
      250k  2.86187e-07
      600k  1.29081e-07   <- the 1.29e-7 I recorded
      1.5M  2.02361e-07

That is a reference carrying ~9.2e-08 of its own grid error against a 2.16e-07
signal, not a mesh disagreement.

WHAT IS ACTUALLY TRUE.  Both meshes converge, to the same value:

    uniform  vs c=5:  1.559 / 1.568 / 1.931 / 2.060 / 2.126 / 2.14665 e-07
                      at 100k / 250k / 400k / 600k / 1M / 1.5M   (step 1.010)
    graded   vs c=5:  2.045 / 2.134 / 2.15582 e-07
                      at 60k / 120k / 240k                       (step 1.010)
    graded, absolute against a graded c=12 reference:  2.16319e-07

    cross-mesh check: the c=5 ABSOLUTE value reads 17.2714838375 (graded 240k)
    against 17.2714838504 (uniform 1.5M) -- agreement 7.5e-10.

So the converged c=3 truncation is 2.15e-07 (spread 0.4% across meshes, 0.8%
including the absolute route), the meshes AGREE, and the third digit IS
determined.  My uniform ladder stopped at 600k while still climbing at 1.067x
and I read "still moving" as "cannot be determined".

The guard's ASSERTIONS stay as they are: it runs at 100k and 600k, where the
value is genuinely still climbing, so order-plus-non-collapse is the right
thing to assert THERE.  What changes is every statement of mechanism.

Idempotent.
"""
from __future__ import annotations

import sys

T = "tests/test_paper60_box_rule.py"
MTX = "docs/claim_test_matrix.md"
CL = "CHANGELOG.md"

DOC_OLD = """**What is established is the ORDER, and that is what the appendix claims.**
Uniform mesh against a c=12 reference: 1.29e-07.  A graded mesh: 2.16e-07.  The
in-repo relative quantity: 1.6-2.1e-07.  So App. A's ~1e-7 is VINDICATED, and
the third digit is NOT determined -- the two discretizations disagree by ~1.7x.
Do not quote a converged value here; an earlier version of this docstring quoted
1.6e-7 and was wrong.

Read to 100k this column looks like a vanishing grid artifact; read to 250k it
looks converged.  Both readings were made, in that order, and both were wrong."""

DOC_NEW = """**The converged value is 2.15e-07, and App. A's ~1e-7 is VINDICATED.**  Both
meshes converge and they AGREE:

    uniform vs c=5, 1M / 1.5M       :  2.126e-07 / 2.14665e-07  (step 1.010)
    graded  vs c=5, 120k / 240k     :  2.134e-07 / 2.15582e-07  (step 1.010)
    graded, absolute vs a graded c=12:  2.16319e-07
    cross-check: the c=5 ABSOLUTE value is 17.2714838375 (graded 240k) against
    17.2714838504 (uniform 1.5M) -- the two meshes agree to 7.5e-10.

The graded mesh `r = R t^2` resolves the near-origin region where the integrand
lives, so it converges at ~6x fewer points; the uniform mesh needs ~1.5M.

**Four mechanisms have now been offered for this one number and the first three
were wrong** -- "converged" (3 points), "below the quadrature floor" (3 points,
opposite direction), "PLATEAUS at 1.5676e-07" (4 points), and "the meshes
disagree by 1.7x so the value is undetermined" (they agree to 0.4%).  The last
came from referencing c=3 against a c=12 box at the SAME point count, i.e. with
4x coarser spacing than the thing it was resolving; that ratio is not even
monotone in npts (2.86e-07 / 1.29e-07 / 2.02e-07 at 250k / 600k / 1.5M).  The
paper's claim was right every time; only the evidence kept failing."""

TEST_OLD = '''    So this asserts what is robust across discretizations, which is also what
    App. A actually claims: the ORDER, and the non-collapse.  It deliberately
    asserts no plateau and no third digit, because the uniform and graded
    meshes disagree by ~1.7x on the converged value (1.29e-7 vs 2.16e-7).
    """'''

TEST_NEW = '''    The converged value IS known -- 2.15e-07, with the uniform and graded
    meshes agreeing to 0.4% (2.14665e-07 at 1.5M points against 2.15582e-07 at
    240k graded).  An earlier version of this docstring said they disagreed by
    1.7x; that came from referencing c=3 against a c=12 box at the same point
    count, hence 4x coarser spacing than the quantity it was resolving.

    This test nonetheless asserts only the ORDER and the NON-COLLAPSE, because
    it runs at 100k and 600k where the column is still genuinely climbing
    (1.559e-07 -> 2.060e-07, step 1.067x).  Pinning 2.15e-07 at these
    resolutions would be pinning a number the test cannot reach.  See
    `test_c3_converged_value_on_the_graded_mesh` for the leg that does reach it.
    """'''

MTX_OLD = ("**Measurement note (THIRD reading, 2026-09-13):** App. A's ~1e-7 for c=3 is "
           "**VINDICATED as an ORDER**, and the third digit is **not determined**. "
           "Ladder: 2.57e-5 / 2.11e-6 / 1.559e-7 / 1.568e-7 / 1.931e-7 / 2.060e-7 at "
           "12k / 40k / 100k / 250k / 400k / 600k. **c=1 and c=2 plateau; c=3 does "
           "not** — the flat-looking 1.005x step at 250k is a CROSSING, because this "
           "quantity is relative to a c=5 box whose own truncation (7.69e-8 absolute, "
           "vs c=3's 1.29e-7) is comparable and partially cancels there. Absolute "
           "against a c=12 reference: uniform mesh 1.29e-7, graded mesh 2.16e-7 — the "
           "two discretizations differ by ~1.7x, so no converged value is quoted. ")

MTX_NEW = ("**Measurement note (FOURTH reading, 2026-09-13 — and the first three were "
           "wrong):** App. A's ~1e-7 for c=3 is **VINDICATED**, and the converged value "
           "is **2.15e-7**. Uniform ladder: 2.57e-5 / 2.11e-6 / 1.559e-7 / 1.568e-7 / "
           "1.931e-7 / 2.060e-7 / 2.126e-7 / 2.14665e-7 at 12k…1.5M (step 1.010). "
           "Graded mesh `r=Rt²`: 2.045e-7 / 2.134e-7 / 2.15582e-7 at 60k/120k/240k. "
           "**The two meshes AGREE to 0.43%**, and their c=5 absolute values agree to "
           "7.5e-10. The flat-looking 1.005x step at 250k is a CROSSING, not a plateau. "
           "**Four mechanisms were offered for this one number and only the paper's own "
           "claim survived all four:** \"converged\" (3 points), \"below the quadrature "
           "floor\" (3 points, opposite direction), \"PLATEAUS at 1.5676e-7\" (4 points), "
           "and \"the meshes disagree 1.7x so the value is undetermined\" — the last from "
           "referencing c=3 against a c=12 box at the SAME point count, i.e. 4x coarser "
           "spacing than the quantity being resolved, a ratio that is not even monotone "
           "(2.86e-7 / 1.29e-7 / 2.02e-7 at 250k/600k/1.5M). ")


def main() -> int:
    loaded: dict[str, str] = {}
    applied, skipped, missed = [], [], []
    for path, name, marker, old, new in [
        (T, "module-docstring", "The converged value is 2.15e-07", DOC_OLD, DOC_NEW),
        (T, "guard-docstring", "The converged value IS known", TEST_OLD, TEST_NEW),
        (MTX, "matrix-note", "FOURTH reading", MTX_OLD, MTX_NEW),
    ]:
        if path not in loaded:
            with open(path, encoding="utf-8") as fh:
                loaded[path] = fh.read()
        t = loaded[path]
        if marker in t:
            skipped.append(name)
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
    for n in skipped:
        print(f"  skip  {n} (already applied)")
    for n, c in missed:
        print(f"  MISS  {n}: anchor count={c}")
    print(f"applied {len(applied)}, skipped {len(skipped)}, MISSED {len(missed)}")
    return 3 if missed else 0


if __name__ == "__main__":
    sys.exit(main())
