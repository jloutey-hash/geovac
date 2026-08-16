"""Bet 2, extension: does the (AA|BB) "no accidental zeros" survive n_max=3?

n_max=3 brings 3s/3p/3d, i.e. RADIAL NODES -- the first place a genuine
accidental cancellation could occur, and the regime Paper 58 leaves "counted
only". Full decision is intractable (~10^5 permitted entries, symbolic), so we
SAMPLE the permitted entries that contain at least one n=3 orbital (the new
ones), seeded for reproducibility, and decide them by the same Lindemann
separation.

Run:  python debug/bet2_nmax3_sample.py --n 120
"""

from __future__ import annotations

import argparse
import random
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

# n_max = 3 : 1 + 4 + 9 = 14 orbitals per center
ORBS3 = [(1, 0, 0),
         (2, 0, 0), (2, 1, -1), (2, 1, 0), (2, 1, 1),
         (3, 0, 0), (3, 1, -1), (3, 1, 0), (3, 1, 1),
         (3, 2, -2), (3, 2, -1), (3, 2, 0), (3, 2, 1), (3, 2, 2)]

CONFIGS = [("census-like 3:1", 3, 1), ("equal 2:2", 2, 2)]
R_LIST = [sp.Integer(3), sp.Rational(5, 2)]


def side_permitted(a, b, Mv):
    l1, l2 = a[1], b[1]
    return any((l1 + l2 + L) % 2 == 0 and abs(Mv) <= L
               for L in range(abs(l1 - l2), l1 + l2 + 1))


def gaunt_all_zero(Za, a, b, Zb, c, d):
    tA = multipole_decomposition(Za, *a, Za, *b)
    tB = multipole_decomposition(Zb, *c, Zb, *d)
    for LA, MA, gA, _r, _bb in tA:
        for LB, MB, gB, _r2, _b2 in tB:
            if MA + MB == 0 and sp.simplify(gA * gB) != 0:
                return False
    return True


def decide_zero(expr, Rval):
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


def permitted_with_n3(orbs):
    out = []
    for p, q, r, s in product(range(len(orbs)), repeat=4):
        a, b, c, d = orbs[p], orbs[q], orbs[r], orbs[s]
        if max(a[0], b[0], c[0], d[0]) < 3:      # only the genuinely-new ones
            continue
        MA, MB = b[2] - a[2], d[2] - c[2]
        if MA + MB != 0:
            continue
        if side_permitted(a, b, MA) and side_permitted(c, d, MB):
            out.append((a, b, c, d))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=120, help="sample size per config")
    ap.add_argument("--seed", type=int, default=58)
    args = ap.parse_args()

    pool = permitted_with_n3(ORBS3)
    print(f"Bet 2 n_max=3 extension: permitted entries containing an n=3 orbital"
          f" = {len(pool)}")
    rng = random.Random(args.seed)
    sample = rng.sample(pool, min(args.n, len(pool)))
    print(f"sampling {len(sample)} (seed {args.seed})\n")

    for label, ZA, ZB in CONFIGS:
        print(f"=== {label}  (Z_A={ZA}, Z_B={ZB}) ===")
        pc = decide_zero(aabb_closed_form(Fraction(ZA), (2, 1, 1), (2, 1, -1),
                                          Fraction(ZB), (2, 1, 0), (2, 1, 0)),
                         sp.Integer(3))
        print(f"  positive control decided zero: {pc}")
        t0 = time.time()
        gz = 0
        sym = []
        for a, b, c, d in sample:
            if gaunt_all_zero(Fraction(ZA), a, b, Fraction(ZB), c, d):
                gz += 1
                continue
            sym.append(((a, b, c, d), aabb_closed_form(
                Fraction(ZA), a, b, Fraction(ZB), c, d)))
        print(f"  gaunt-zero(angular)={gz}  radially-live={len(sym)}  "
              f"(build {time.time()-t0:.0f}s)")
        for Rval in R_LIST:
            acc = []
            for key, e in sym:
                if decide_zero(e, Rval):
                    acc.append(key)
            flag = f"  <-- {len(acc)} ACCIDENTAL" if acc else "  clean"
            print(f"    R={str(Rval):>4}: nonzero={len(sym)-len(acc):3d} "
                  f"accidental={len(acc):2d}{flag}")
            for key in acc[:8]:
                print(f"        {key}")
        print()


if __name__ == "__main__":
    main()
