"""FD1-FD3 step I: the order-9 right factor W9 of L (the minimal annihilator of the log-free
exponent-5 Frobenius solution y5 of L at t = 0, free slots 0; order 9 found mod p in step E).
(a) minimal t-degree of an order-9 annihilator of y5, mod p;
(b) exact W9 from the exact series (nullspace over QQ, smallest degree);
(c) L = N o W9 exactly; W9 = X5 o R4 exactly (R4 from step H, normalized); R4 = R3 o (theta - 3).
Saves Q-forms to four_density_star_I.json."""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction
from pathlib import Path

import sympy as sp

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import density_loop_C_factor as C  # noqa: E402
import four_density_star_F_hyperexp as F  # noqa: E402
import four_density_star_G_peel1 as G  # noqa: E402
import four_density_star_H_order4 as H  # noqa: E402
import four_density_star_ops as O  # noqa: E402

u = O.u
t = sp.Symbol("t")
P = C.P
OUT = HERE / "data" / "four_density_star_I.json"


def nullity_mod(ser, order, I, nterms):
    lam = min(ser)
    cols = []
    for j in range(order + 1):
        for i in range(I + 1):
            cols.append({e + 2 * i: c * pow(e % P, j, P) % P for e, c in ser.items()})
    exps = [lam + 2 * n for n in range(nterms)]
    rows = [[col.get(e, 0) for col in cols] for e in exps]
    return len(cols) - C.rank_mod(rows)


def apply_Q(Qf, ser):
    out = {}
    for i, q in enumerate(Qf):
        pq = sp.Poly(q, u)
        for e, c in ser.items():
            v = Fraction(str(pq.eval(e)))
            if v:
                out[e + 2 * i] = out.get(e + 2 * i, 0) + v * c
    top = max(ser) + 1
    return {e: c for e, c in out.items() if c and e < top}


def primitive(Q):
    cs = []
    for q in Q:
        cs += [sp.Rational(c) for c in sp.Poly(q, u).coeffs()]
    den = sp.ilcm(*[c.q for c in cs])
    num = sp.igcd(*[(c * den).p for c in cs])
    return [sp.expand(q * den / num) for q in Q]


def main() -> None:
    rec = {}
    t0 = time.time()
    Q = [sp.sympify(q) for q in json.loads((HERE / "data" / "four_density_star_A.json").read_text())["L"]["Q"]]
    L = O.from_Q(Q)
    # R4 normalized
    R4Q = primitive([sp.sympify(q) for q in json.loads((HERE / "data" / "four_density_star_H.json").read_text())["R4_Q"]])
    print("R4:", [sp.factor(q) for q in R4Q], flush=True)
    rec["R4_Q"] = [str(q) for q in R4Q]
    R4 = O.from_Q(R4Q)
    R3, rem = O.rdiv(R4, O.from_Q([u - 3]))
    assert all(c == 0 for c in rem)
    R3Q, _ = O.to_Q(R3)
    R3Q = primitive(R3Q)
    print("R3 = R4/(theta-3):", [sp.factor(q) for q in R3Q], flush=True)
    rec["R3_Q"] = [str(q) for q in R3Q]
    # (a) z = R4 y5 (y5: exact log-free exponent-5 Frobenius solution of L, free slots 0)
    y5 = G.frob_exact(Q, 5, 2 * 260)
    z = apply_Q(R4Q, y5)
    zp = {e: c.numerator % P * pow(c.denominator % P, P - 2, P) % P for e, c in z.items()}
    print(f"z = R4 y5: lowest exponent {min(z)}  ({time.time()-t0:.0f}s)", flush=True)
    best = None
    for order in (1, 2, 3, 4, 5):
        for I in range(2, 16):
            n = nullity_mod(zp, order, I, 230)
            if n:
                best = (order, I)
                break
        print(f"  z: order {order}: {'nullity %d at t-degree %d' % (n, 2*I) if best else 'none up to t-degree 30'}", flush=True)
        if best:
            break
    rec["X5_order_degree"] = best
    # (b) exact X5 from z
    order, I = best
    sols = H.guess(z, order, I, 200)
    print(f"exact nullity order {order}, degree {2*I}: {len(sols)}  ({time.time()-t0:.0f}s)", flush=True)
    X5Q = primitive(sols[0])
    while X5Q[-1] == 0:
        X5Q.pop()
    X5 = O.from_Q(X5Q)
    W9 = O.compose(X5, R4)
    W9Q, _ = O.to_Q(W9)
    W9Q = primitive(W9Q)
    rec["W9_Q"] = [str(q) for q in W9Q]
    N10, rem = O.rdiv(L, W9)
    rec["L_div_W9_zero"] = all(c == 0 for c in rem)
    print("L = N o (X5 o R4) exact:", rec["L_div_W9_zero"], f"({time.time()-t0:.0f}s)", flush=True)
    X5Q, _ = O.to_Q(X5)
    X5Q = primitive(X5Q)
    N10Q, _ = O.to_Q(N10)
    rec["X5_Q"] = [str(q) for q in X5Q]
    rec["N10_Q"] = [str(q) for q in primitive(N10Q)]
    print("X5:", [sp.factor(q) for q in X5Q], flush=True)
    lead = lambda QQ, k: sp.factor(sum(sp.Poly(q, u).coeff_monomial(u**k) * t ** (2 * i) for i, q in enumerate(QQ)))
    print("X5 order", len(X5) - 1, "leading:", lead(X5Q, len(X5) - 1), flush=True)
    print("W9 order", len(W9) - 1, "leading:", lead(W9Q, len(W9) - 1), flush=True)
    OUT.write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
