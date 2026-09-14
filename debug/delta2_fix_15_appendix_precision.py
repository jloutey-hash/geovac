"""DELTA #2 -- fix the upgrade I just wrote, which over-stated its own precision.

`delta2_fix_14_final.py` upgraded App. A from "c=3 gives ~1e-7" to
"2.16e-7 ... cross-checked on two independent quadrature rules that agree to
0.08%".  The value is right; the AGREEMENT FIGURE is attached to the wrong
quantity.  Two quantities are in play:

  * ABSOLUTE box truncation, against a well-resolved c=12 reference:
        uniform @1.5M  2.16145e-07
        graded  @240k  2.16319e-07      -> agree to 0.080%
  * RELATIVE to a c=5 box -- which is what App. A's sentence actually states,
    since the whole paragraph is "measured relative error against a
    R_max = 5 n_max^2 reference":
        uniform @1.5M  2.14665e-07
        graded  @240k  2.15582e-07      -> agree to 0.43%

So 0.08% belongs to the absolute form, not to the form the sentence quotes, and
on the relative form the two rules differ in the THIRD digit (2.15 vs 2.156).
Quoting "2.16e-7 ... agree to 0.08%" for the relative quantity claims a
precision the relative measurement does not have.

Restated to name both forms and give each its own spread.  This is the fifth
correction to this one number in a day, and the second where the defect was
introduced by the correction of the previous one -- which is the argument for
stating the quantity a number refers to, every time, rather than the number
alone.

Idempotent.
"""
from __future__ import annotations

import sys

P = "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex"

OLD = (
    "$c=3$ gives $2.16\\times10^{-7}$ --- converged, and cross-checked on two\n"
    "independent quadrature rules (a uniform mesh at $1.5\\times10^{6}$ points and a\n"
    "graded $r=R t^{2}$ mesh at $2.4\\times10^{5}$) that agree to $0.08\\%$;\\ an\n"
    "earlier version of this sentence gave only the order."
)

NEW = (
    "$c=3$ gives $2.15$--$2.16\\times10^{-7}$, converged.  Two independent quadrature\n"
    "rules were run to convergence:\\ a uniform mesh at $1.5\\times10^{6}$ points and a\n"
    "graded $r=Rt^{2}$ mesh at $2.4\\times10^{5}$.  On the relative form quoted here\n"
    "they read $2.147$ and $2.156\\times10^{-7}$ (a spread of $0.43\\%$);\\ on the\n"
    "\\emph{absolute} box truncation, measured against a well-resolved\n"
    "$c=12$ reference, they read $2.161$ and $2.163\\times10^{-7}$ ($0.08\\%$).  An\n"
    "earlier version of this sentence gave only the order of magnitude."
)


def main() -> int:
    with open(P, encoding="utf-8") as fh:
        t = fh.read()
    if "on the relative form quoted here" in t:
        print("already applied")
        return 0
    n = t.count(OLD)
    if n != 1:
        print(f"  MISS: anchor count={n}")
        return 3
    with open(P, "w", encoding="utf-8") as fh:
        fh.write(t.replace(OLD, NEW))
    print("  ok    App. A now names both forms and gives each its own spread")
    return 0


if __name__ == "__main__":
    sys.exit(main())
