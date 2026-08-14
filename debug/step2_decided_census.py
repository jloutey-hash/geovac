"""Step 2: a DECIDED census of the (AA|BB) block -- zeros proven, not thresholded.

Paper 58's `g` row is COUNTED: it says which entries the selection rules PERMIT,
never whether each permitted entry is actually nonzero. Its own text is explicit
that "the $g$ figures are ... claims about how many entries the selection rules
permit, not verifications that each permitted entry is nonzero", and the
McMurchie-Davidson corroboration decides zeros by a 1e-10 float threshold rather
than deciding them.

That gap is now closable for the (AA|BB) block. Increment 1c gives the class in
closed form, and every pi cancels, leaving

    (ab|cd) = A_0(R) + sum_j A_j(R) e^{-lambda_j R},   A_j rational, lambda_j rational

so for algebraic R the family {1, e^{-lambda_j R}} is linearly independent over
the algebraics (Lindemann-Weierstrass) and the entry vanishes IFF every A_j does.
That is a decision, not a measurement.

THREE OUTCOMES, and the third is the interesting one:

  nonzero            decided nonzero
  zero, by Gaunt     the crude count permits it, but the actual Gaunt
                     coefficient vanishes -- a symmetry zero the counting missed
  zero, ACCIDENTAL   Gaunt coefficients nonzero, yet the radial combination
                     cancels. A genuine accidental zero, invisible to counting
                     and indistinguishable from a small number to thresholding.

Census configuration is Paper 58's own: Z_A = 3, Z_B = 1, n_max = 2, R = 3 bohr.

Run from repo root:  python debug/step2_decided_census.py
"""

from __future__ import annotations

import sys
import time
from fractions import Fraction
from itertools import product
from pathlib import Path

import sympy as sp

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from geovac.two_center_eri import (  # noqa: E402
    R_s, aabb_closed_form, multipole_decomposition,
)

Z_A, Z_B, R_CEN = Fraction(3), Fraction(1), sp.Integer(3)
ORBS = [(1, 0, 0), (2, 0, 0), (2, 1, -1), (2, 1, 0), (2, 1, 1)]


def side_permitted(a, b, Mv):
    """The crude per-side Gaunt feasibility the counting census uses."""
    l1, l2 = a[1], b[1]
    return any((l1 + l2 + L) % 2 == 0 and abs(Mv) <= L
               for L in range(abs(l1 - l2), l1 + l2 + 1))


def gaunt_all_zero(Za, a, b, Zb, c, d):
    """True if every surviving (L_A,M_A)x(L_B,M_B) pair has a vanishing product."""
    tA = multipole_decomposition(Za, *a, Za, *b)
    tB = multipole_decomposition(Zb, *c, Zb, *d)
    for LA, MA, gA, _r, _bb in tA:
        for LB, MB, gB, _r2, _b2 in tB:
            if MA + MB == 0 and sp.simplify(gA * gB) != 0:
                return False
    return True


def decide_zero(expr, Rval):
    """Lindemann decision: group by exponential rate, require every A_j(R) = 0."""
    groups: dict = {}
    for term in sp.Add.make_args(sp.expand(expr)):
        rate, coeff = sp.Integer(0), sp.Integer(1)
        for f in sp.Mul.make_args(term):
            if isinstance(f, sp.exp):
                arg = sp.expand(f.args[0])
                d = -sp.diff(arg, R_s)
                rate += d
                coeff *= sp.exp(sp.expand(arg + d * R_s))
            else:
                coeff *= f
        key = sp.nsimplify(rate)
        groups[key] = groups.get(key, 0) + coeff
    return all(sp.simplify(c.subs(R_s, Rval)) == 0 for c in groups.values())


def main() -> None:
    print("Step 2 -- DECIDED census of the (AA|BB) block\n")
    print(f"Paper 58 config: Z_A = {Z_A}, Z_B = {Z_B}, n_max = 2, R = {R_CEN} bohr")
    print("zeros DECIDED by Lindemann separation, not thresholded\n")

    permitted = []
    for p, q, r, s in product(range(len(ORBS)), repeat=4):
        a, b, c, d = (ORBS[p], ORBS[q], ORBS[r], ORBS[s])
        MA, MB = b[2] - a[2], d[2] - c[2]
        if MA + MB != 0:
            continue
        if side_permitted(a, b, MA) and side_permitted(c, d, MB):
            permitted.append((a, b, c, d))
    print(f"  permitted by the counting rules : {len(permitted)} of 625")

    t0 = time.time()
    n_nonzero = n_gaunt = n_accidental = 0
    accidents = []
    seen: dict = {}
    for a, b, c, d in permitted:
        key = (a, b, c, d)
        if key in seen:
            continue
        if gaunt_all_zero(Z_A, a, b, Z_B, c, d):
            n_gaunt += 1
            seen[key] = "gaunt"
            continue
        e = aabb_closed_form(Z_A, a, b, Z_B, c, d)
        if decide_zero(e, R_CEN):
            n_accidental += 1
            accidents.append(key)
            seen[key] = "accidental"
        else:
            n_nonzero += 1
            seen[key] = "nonzero"

    n = len(seen)
    print(f"  decided ({time.time()-t0:.0f}s):")
    print(f"    nonzero                        : {n_nonzero:4d}  ({100*n_nonzero/n:5.1f}%)")
    print(f"    zero, Gaunt coefficient        : {n_gaunt:4d}  ({100*n_gaunt/n:5.1f}%)"
          f"   <- symmetry zeros the counting missed")
    print(f"    zero, ACCIDENTAL               : {n_accidental:4d}  ({100*n_accidental/n:5.1f}%)"
          f"   <- invisible to counting AND to thresholding")
    if accidents:
        print("\n  accidental zeros:")
        for a, b, c, d in accidents[:12]:
            print(f"    ({a}{b}|{c}{d})")
        if len(accidents) > 12:
            print(f"    ... and {len(accidents)-12} more")

    genuine = n_nonzero
    print(f"\n  So of {n} entries the counting calls permitted, {genuine} are")
    print(f"  genuinely nonzero: the counted density overstates the true density")
    print(f"  by a factor {n/genuine:.3f} on this block." if genuine else "")


if __name__ == "__main__":
    main()
