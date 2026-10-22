"""FD1-FD3 step H: the exponent-2 log-free solution y2 of L at t = 0 has minimal annihilator order 4
(step E, mod p, degree <= 24).  Make it exact: guess R4 = sum_{j<=4, i<=12} a_{ji} t^(2i) theta^j with
R4 y2 = 0 from the exact Frobenius series (nullspace over QQ), check L = M o R4 exactly (zero
remainder), and compare R4 with F_4's operator Q4 = 9u^2(u-2)^2 - 2u^2(5u^2+1) t^2 + (u+1)^4 t^4."""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction
from pathlib import Path

import sympy as sp
from sympy.polys.matrices import DomainMatrix

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import four_density_star_G_peel1 as G  # noqa: E402
import four_density_star_ops as O  # noqa: E402

u = O.u
t = sp.Symbol("t")
OUT = HERE / "data" / "four_density_star_H.json"


def guess(ser, order, I, nterms):
    """operators sum_{j<=order, i<=I} a_{ji} t^(2i) theta^j killing ser (dict exp -> Fraction)."""
    lam = min(ser)
    cols, labels = [], []
    for j in range(order + 1):
        for i in range(I + 1):
            col = {}
            for e, c in ser.items():
                col[e + 2 * i] = c * Fraction(e) ** j
            cols.append(col)
            labels.append((j, i))
    exps = [lam + 2 * n for n in range(nterms)]
    rows = [[sp.Rational(col.get(e, Fraction(0)).numerator, col.get(e, Fraction(0)).denominator) for col in cols] for e in exps]
    M = DomainMatrix.from_list_sympy(len(rows), len(cols), rows).convert_to(sp.QQ)
    ns = M.nullspace().to_Matrix()
    out = []
    for r in range(ns.shape[0]):
        Q = [0] * (I + 1)
        for v, (j, i) in zip(ns.row(r), labels):
            Q[i] += v * u**j
        out.append([sp.expand(q) for q in Q])
    return out


def main() -> None:
    rec = {}
    Q = [sp.sympify(q) for q in json.loads((HERE / "data" / "four_density_star_A.json").read_text())["L"]["Q"]]
    t0 = time.time()
    NT = 2 * 150
    y2 = G.frob_exact(Q, 2, NT)
    print(f"y2: {len(y2)} terms ({time.time()-t0:.0f}s); first: {[str(y2[e]) for e in sorted(y2)[:4]]}", flush=True)
    for order in (3, 4):
        sols = guess(y2, order, 12, 140)
        print(f"order {order}, t-degree <= 24: nullity {len(sols)}  ({time.time()-t0:.0f}s)", flush=True)
        if sols and order == 4:
            break
    # take the solution of minimal t-degree: reduce by gcd and pick the lowest-degree combination
    R4Q = sols[0]
    while R4Q and R4Q[-1] == 0:
        R4Q.pop()
    g = sp.gcd_list([sp.Poly(q, u).content() for q in R4Q if q != 0] or [1])
    print("R4 (Q-form):", [sp.factor(q) for q in R4Q], flush=True)
    rec["R4_Q"] = [str(q) for q in R4Q]
    rec["nullity_order4_deg24"] = len(sols)
    # exact right division
    L = O.from_Q(Q)
    R4 = O.from_Q(R4Q)
    M, rem = O.rdiv(L, R4)
    rec["L_div_R4_rem_zero"] = all(c == 0 for c in rem)
    print("L = M o R4 exact:", rec["L_div_R4_rem_zero"], f"({time.time()-t0:.0f}s)", flush=True)
    MQ, _ = O.to_Q(M)
    rec["M_Q"] = [str(q) for q in MQ]
    rec["M_order"] = len(M) - 1
    # compare with F_4's operator directly
    Q4 = [9 * u**2 * (u - 2) ** 2, -2 * u**2 * (5 * u**2 + 1), (u + 1) ** 4]
    rec["R4_equals_Q4_up_to_scalar"] = len(R4Q) == 3 and all(sp.expand(a * Q4[0].coeff(u, 4) - b * R4Q[0].coeff(u, 4)) == 0 for a, b in zip(R4Q, Q4))
    print("R4 == Q4 up to scalar:", rec["R4_equals_Q4_up_to_scalar"], flush=True)
    OUT.write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
