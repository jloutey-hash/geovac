"""Bet 2: is the (AA|BB) "no accidental zeros" result special to the census config?

Paper 58 / Step 2 DECIDED the (AA|BB) block on ONE configuration
(Z_A=3, Z_B=1, n_max=2, R=3): all 195 permitted entries genuinely nonzero,
0 accidental zeros, 0 missed symmetry zeros -> counted density == true density.

This sweep asks whether that is general or config-specific. Two facts, read off
the machinery, focus the experiment:

  * The PERMITTED set is charge- and R-independent (pure angular Gaunt/m
    feasibility). So it is 195 for every (Z_A, Z_B) at n_max=2.
  * The "Gaunt symmetry zero" count (n_gaunt) is angular too -> charge-independent.

So the ONLY quantity that can move across configs is the count of ACCIDENTAL
radial zeros (Gaunt nonzero, closed-form radial sum cancels). The sharpest places
one could hide:
  - equal charges Z_A = Z_B (degenerate cross-center exponential rates; and
    N2/F2 in Paper 58's own swap table are homonuclear);
  - rate-coincidence charges (4:2 puts a 2s_A exponent = 1s_B exponent);
  - special rational R (a coefficient polynomial hitting a root).

Efficiency: aabb_closed_form(..., no R) returns the symbolic form in R_s, so we
build each entry ONCE per charge pair and decide it at every R cheaply.

Positive control included: an M-violating quartet is truly zero, and the decider
must flag it -- so "0 accidental" means "none found", not "decider blind".

Run from repo root:  python debug/bet2_aabb_accidental_zero_sweep.py
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

ORBS_NMAX2 = [(1, 0, 0), (2, 0, 0), (2, 1, -1), (2, 1, 0), (2, 1, 1)]

# (label, Z_A, Z_B) -- control first, then the sharp probes
CONFIGS = [
    ("census 3:1 (control)", 3, 1),
    ("equal   1:1",          1, 1),
    ("equal   2:2",          2, 2),
    ("hetero  2:1",          2, 1),
    ("ratecoin 4:2",         4, 2),
]

R_LIST = [sp.Integer(1), sp.Integer(2), sp.Integer(3), sp.Integer(4),
          sp.Integer(5), sp.Rational(5, 2), sp.Rational(7, 3)]


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
    """Lindemann decision: group by exponential rate, require every A_j(R)=0."""
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


def permitted_entries(orbs):
    out = []
    for p, q, r, s in product(range(len(orbs)), repeat=4):
        a, b, c, d = orbs[p], orbs[q], orbs[r], orbs[s]
        MA, MB = b[2] - a[2], d[2] - c[2]
        if MA + MB != 0:
            continue
        if side_permitted(a, b, MA) and side_permitted(c, d, MB):
            out.append((a, b, c, d))
    return out


def positive_control(Z_A, Z_B):
    """An M-violating quartet is identically zero; the decider must catch it."""
    e = aabb_closed_form(Fraction(Z_A), (2, 1, 1), (2, 1, -1),
                         Fraction(Z_B), (2, 1, 0), (2, 1, 0))
    return decide_zero(e, sp.Integer(3))


def run_config(label, Z_A, Z_B, orbs, entries, limit=None):
    print(f"\n=== {label}   (Z_A={Z_A}, Z_B={Z_B}, n_max=2) ===")
    ok = positive_control(Z_A, Z_B)
    print(f"  positive control (M-violating quartet decided zero): {ok}")
    if not ok:
        print("  !! decider FAILED its positive control -- results below are suspect")

    todo = entries if limit is None else entries[:limit]
    t0 = time.time()

    # classify each entry once: gaunt-zero, or keep its symbolic closed form
    gaunt_zeros = 0
    sym = []  # (key, expr) for the radially-live entries
    for a, b, c, d in todo:
        if gaunt_all_zero(Fraction(Z_A), a, b, Fraction(Z_B), c, d):
            gaunt_zeros += 1
            continue
        e = aabb_closed_form(Fraction(Z_A), a, b, Fraction(Z_B), c, d)
        sym.append(((a, b, c, d), e))
    build_s = time.time() - t0
    print(f"  permitted={len(todo)}  gaunt-zero(angular)={gaunt_zeros}  "
          f"radially-live={len(sym)}  (build {build_s:.0f}s)")

    # decide accidental zeros at every R, reusing the symbolic forms
    accidental_any = set()
    print(f"  {'R':>6} | {'nonzero':>8} {'accidental':>11}")
    for Rval in R_LIST:
        nz = acc = 0
        acc_here = []
        for key, e in sym:
            if decide_zero(e, Rval):
                acc += 1
                acc_here.append(key)
                accidental_any.add((key, str(Rval)))
            else:
                nz += 1
        flag = "  <-- ACCIDENTAL ZERO(S)" if acc else ""
        print(f"  {str(Rval):>6} | {nz:8d} {acc:11d}{flag}")
        for key in acc_here[:6]:
            a, b, c, d = key
            print(f"           ({a}{b}|{c}{d})")

    return {"label": label, "gaunt": gaunt_zeros, "live": len(sym),
            "accidental": accidental_any}


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None,
                    help="cap entries per config (smoke test)")
    ap.add_argument("--only", type=int, default=None,
                    help="run only CONFIGS[:only]")
    args = ap.parse_args()

    orbs = ORBS_NMAX2
    entries = permitted_entries(orbs)
    print(f"Bet 2 -- (AA|BB) accidental-zero sweep")
    print(f"permitted entries at n_max=2 (charge/R-independent): {len(entries)}")

    cfgs = CONFIGS if args.only is None else CONFIGS[:args.only]
    results = []
    for label, ZA, ZB in cfgs:
        results.append(run_config(label, ZA, ZB, orbs, entries, limit=args.limit))

    print("\n================ SUMMARY ================")
    for r in results:
        n_acc = len(r["accidental"])
        verdict = "CLEAN (no accidental zeros)" if n_acc == 0 else \
            f"{n_acc} accidental (entry,R) hits"
        print(f"  {r['label']:24s}  live={r['live']:3d}  gaunt={r['gaunt']:2d}  "
              f"-> {verdict}")
        for (key, Rs) in list(r["accidental"])[:8]:
            print(f"        {key} at R={Rs}")


if __name__ == "__main__":
    main()
