"""FD1-FD3 step K: which reference families sit where.  L = M15 o R4 (step H, exact), R4 = R3 o (theta-3).
(1) R3 vs F_3 (sub and quotient);  (2) M15 = L/R4: algebraic solutions; sub tests for F_4, F_3, Legendre.
mod p (step J machinery)."""
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
import four_density_star_I_w9 as I  # noqa: E402
import four_density_star_J_gauge as J  # noqa: E402
import four_density_star_ops as O  # noqa: E402

u = O.u
OUT = HERE / "data" / "four_density_star_K.json"


def main() -> None:
    rec = {}
    t0 = time.time()
    Q = [sp.sympify(q) for q in json.loads((HERE / "data" / "four_density_star_A.json").read_text())["L"]["Q"]]
    R4Q = I.primitive([sp.sympify(q) for q in json.loads((HERE / "data" / "four_density_star_H.json").read_text())["R4_Q"]])
    L = O.from_Q(Q)
    R4 = O.from_Q(R4Q)
    M15, rem = O.rdiv(L, R4)
    assert all(c == 0 for c in rem)
    M15Q, _ = O.to_Q(M15)
    M15Q = I.primitive(M15Q)
    rec["M15_Q"] = [str(q) for q in M15Q]
    R3, _ = O.rdiv(R4, O.from_Q([u - 3]))
    R3Q = I.primitive(O.to_Q(R3)[0])
    rec["R3_Q"] = [str(q) for q in R3Q]
    print(f"M15 order {len(M15)-1}, {len(M15Q)} t-terms; t^0: {sp.factor(M15Q[0])}; top: {sp.factor(M15Q[-1])}  ({time.time()-t0:.0f}s)", flush=True)
    NT = 700
    R = J.refs(NT)
    # (1) R3 vs F3
    for side, QQ, f in (("R", R3Q, R["F3"][0]), ("L", J.adjoint(R3Q), J.adjoint_ref("F3", NT))):
        res = J.scan(QQ, f, 3, 16, maxpow=2, NT=NT)
        rec[f"R3|F3|{side}"] = res
        print(f"R3 vs F3 ({side}): {res}  ({time.time()-t0:.0f}s)", flush=True)
    # (2) M15
    mu = max(r for r in sp.roots(sp.Poly(M15Q[-1], u).as_expr(), u) if r.is_rational)
    alg = F.algebraic_solutions(M15Q, [0, -1, -2, -3], [0, Fraction(-1, 2), Fraction(-3, 2)], [0, -1, -2], mu, N=500, verbose=False)
    rec["M15_algebraic"] = alg
    print(f"M15 algebraic solutions: {alg}  ({time.time()-t0:.0f}s)", flush=True)
    for name in ("F4", "F3", "Leg(t^2)", "Leg(t^2/4)", "Leg(t^2/9)"):
        f, r, _ = R[name]
        res = J.scan(M15Q, f, r, 20, maxpow=2, NT=NT)
        rec[f"M15|{name}|R"] = res
        print(f"M15 sub {name}: {res}  ({time.time()-t0:.0f}s)", flush=True)
        OUT.write_text(json.dumps(rec, indent=1, default=str))


if __name__ == "__main__":
    main()
