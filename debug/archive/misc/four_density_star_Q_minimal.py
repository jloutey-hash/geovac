"""FD3 step Q (minimality probe, mod p): S's module maps onto F_4's piece (B5 S != 0, R4' irreducible),
so a right factor of W9 shorter than W9 that annihilates S must be X4 o V with V a proper right factor
of B5 and X4 ~ F_4.  For V in {1, theta-3, L1, R4} test whether F_4 is a sub of W9/V (step M's
extended box, both parities); for each hit the order-(4 + ord V) right factor exists and is then
tested on S (step P machinery).  Also: algebraic solutions of B5 and R4 (other order-1 subs)."""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction
from pathlib import Path

import sympy as sp

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import four_density_star_F_hyperexp as F  # noqa: E402
import four_density_star_J_gauge as J  # noqa: E402
import four_density_star_ops as O  # noqa: E402
import four_density_star_P_score as PS  # noqa: E402

u = O.u
t = sp.Symbol("t")
OUT = HERE / "data" / "four_density_star_Q.json"


def main() -> None:
    rec = {}
    t0 = time.time()
    ops = PS.build()
    NT = 1100
    R = J.refs(NT)
    f4 = R["F4"][0]
    extra = sp.expand((5 * t**2 + 4) ** 3 * (5 * t**4 + 2 * t**2 + 2) ** 3)
    for nm in ["R4", "L1", "TH3", "ONE"]:
        if nm == "ONE":
            quo = ops["W9"]
        else:
            quo, rem = O.rdiv(ops["W9"], ops[nm])
            assert all(c == 0 for c in rem)
        Qq, _ = O.to_Q(quo)
        res = [J.test(Qq, f4, 4, 74, 4, 3, 4, a=10, NT=NT, par=par, extra=extra) for par in (0, 1)]
        rec[f"F4 sub of W9/{nm}"] = res
        print(f"F4 sub of W9/{nm} (order {len(quo)-1}): {res}  ({time.time()-t0:.0f}s)", flush=True)
        OUT.write_text(json.dumps(rec, indent=1))
    for nm in ["R4", "B5"]:
        Qq, _ = O.to_Q(ops[nm])
        mu = max(r for r in sp.roots(sp.Poly(Qq[-1], u).as_expr(), u) if r.is_rational)
        alg = F.algebraic_solutions(Qq, [0, -1, -2, -3], [0, Fraction(-1, 2), Fraction(-3, 2)], [0, -1, -2], mu, N=500, verbose=False)
        rec[f"algebraic solutions of {nm}"] = alg
        print(f"algebraic solutions of {nm}: {alg}", flush=True)
    OUT.write_text(json.dumps(rec, indent=1, default=str))


if __name__ == "__main__":
    main()
