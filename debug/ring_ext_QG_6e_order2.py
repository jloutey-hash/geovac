"""QG step 6e (gamma = 3/2): the order-2 right factor of L9 through the exponent-0 Frobenius solution
y0, exactly: smallest-degree operator sum_{i<=d} t^(2i) P_i(theta), deg P_i <= 2, killing the exact
rational series of y0; then its singular points and local exponents, and an exact right division
L9 = M7 o R2."""
from __future__ import annotations
import json, sys
from fractions import Fraction
from pathlib import Path
import sympy as sp
from sympy.polys.matrices import DomainMatrix
sys.path.insert(0, str(Path(__file__).parent))
import ring_ext_QG_ops as O
from ring_ext_QG_ops import u, z, tt as t

DATA = Path(__file__).parent / "data"


def main():
    rec = json.loads((DATA / "ring_ext_QG_5c_ylog_g1.5_r9_d25.json").read_text())
    Q = [sp.Poly(sp.sympify(q, locals={"u": u}), u) for q in rec["Q"]]
    ev = lambda i, x: Fraction(int(Q[i].eval(x)))
    N = 140
    c = [Fraction(1)]
    for n in range(1, N):
        s = sum(ev(i, 2 * (n - i)) * c[n - i] for i in range(1, len(Q)) if n - i >= 0)
        c.append(-s / ev(0, 2 * n))
    for d in range(1, 30):
        cols = []
        for i in range(d + 1):
            for j in range(3):
                cols.append([(Fraction(2 * (n - i)) ** j * c[n - i]) if n - i >= 0 else Fraction(0) for n in range(N)])
        rows = [[sp.Rational(col[n].numerator, col[n].denominator) for col in cols] for n in range(N)]
        M = DomainMatrix.from_list_sympy(N, len(cols), rows).convert_to(sp.QQ)
        ns = M.nullspace().to_Matrix()
        if ns.shape[0]:
            v = list(ns.row(0))
            QR = [sp.expand(sum(v[3 * i + j] * u**j for j in range(3))) for i in range(d + 1)]
            den = sp.lcm([sp.denom(cf) for q in QR for cf in sp.Poly(q, u).coeffs()])
            QR = [sp.expand(q * den) for q in QR]
            g = sp.gcd_list([cf for q in QR for cf in sp.Poly(q, u).coeffs()])
            QR = [sp.expand(q / g) for q in QR]
            print("order-2 annihilator of y0: t^2-degree", d, "nullity", ns.shape[0], flush=True)
            for i, q in enumerate(QR):
                print("  t^%d:" % (2 * i), sp.factor(q), flush=True)
            A = O.theta_to_D(QR, z, 2)
            print("  leading:", sp.factor(A[-1]), flush=True)
            for t0 in (0, 2, 3, sp.Rational(5, 2), sp.Rational(1, 2)):
                print("  exponents at", t0, O.indicial(A, t0), flush=True)
            print("  at inf", O.indicial_inf(A), flush=True)
            L9c = O.Q_to_coeffs([sp.sympify(q, locals={"u": u}) for q in rec["Q"]])
            R2c = O.Q_to_coeffs(QR)
            Mq, rem = O.rdiv(L9c, R2c)
            print("  L9 = M7 o R2 exactly:", all(sp.simplify(x) == 0 for x in rem), flush=True)
            (DATA / "ring_ext_QG_6e_order2_g1.5.json").write_text(json.dumps(
                {"R2_Q": [str(q) for q in QR], "remainder_zero": all(sp.simplify(x) == 0 for x in rem),
                 "M7": [str(x) for x in Mq]}, indent=1))
            return


if __name__ == "__main__":
    main()
