"""QG step 6f: exact right factor through one Frobenius solution.  Given an operator (theta-form,
t^(2i) Q_i), its (optionally adjoint) Frobenius solution at exponent lam, and an order r: find the
smallest t^2-degree d (mod p first), then the exact operator by rational linear algebra on the exact
series, check the exact right division, and report singular points / exponents.
usage: python ring_ext_QG_6f_rightfactor.py <op-json> <key> <adj|-> <lam> <r> <out-tag>"""
from __future__ import annotations
import json, sys, time
from fractions import Fraction
from pathlib import Path
import sympy as sp
from sympy.polys.matrices import DomainMatrix
sys.path.insert(0, str(Path(__file__).parent))
import ring_ext_QG_ops as O
import density_loop_C_factor as C
from ring_ext_QG_ops import u, z, tt as t

DATA = Path(__file__).parent / "data"
P = C.P


def frob_exact(Q, lam, N):
    Qp = [sp.Poly(q, u) for q in Q]
    ev = lambda i, x: Fraction(str(Qp[i].eval(x)))
    c = [Fraction(1)]
    for n in range(1, N):
        s = sum(ev(i, lam + 2 * (n - i)) * c[n - i] for i in range(1, len(Q)) if n - i >= 0)
        c.append(-s / ev(0, lam + 2 * n))
    return c


def cols_for(c, lam, r, d, N, modp):
    cols = []
    for i in range(d + 1):
        for j in range(r + 1):
            col = []
            for n in range(N):
                m = n - i
                v = Fraction(lam + 2 * m) ** j * c[m] if m >= 0 else Fraction(0)
                col.append(v.numerator % P * pow(v.denominator % P, P - 2, P) % P if modp else v)
            cols.append(col)
    return cols


def main():
    path, key, adj, lam, r, tag = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]), sys.argv[6]
    Q = [sp.sympify(q, locals={"u": u}) for q in json.loads(Path(path).read_text())[key]]
    if adj == "adj":
        Q = [sp.expand(q.subs(u, -u - 2 * i)) for i, q in enumerate(Q)]
    N = 420
    t0 = time.time()
    c = frob_exact(Q, lam, N)
    print("series", f"{time.time()-t0:.0f}s", flush=True)
    dmin = None
    for d in range(1, 80):
        U = (r + 1) * (d + 1)
        if U + 40 > N:
            break
        cols = cols_for(c, lam, r, d, N, True)
        rows = [[col[n] for col in cols] for n in range(N)]
        if U - C.rank_mod(rows) > 0:
            dmin = d
            break
    print("mod-p minimal t^2-degree:", dmin, flush=True)
    if dmin is None:
        return
    cols = cols_for(c, lam, r, dmin, N, False)
    rows = [[sp.Rational(col[n].numerator, col[n].denominator) for col in cols] for n in range(N)]
    M = DomainMatrix.from_list_sympy(N, len(cols), rows).convert_to(sp.QQ)
    ns = M.nullspace().to_Matrix()
    print("exact nullity", ns.shape[0], f"{time.time()-t0:.0f}s", flush=True)
    v = list(ns.row(0))
    QR = [sp.expand(sum(v[(r + 1) * i + j] * u**j for j in range(r + 1))) for i in range(dmin + 1)]
    den = sp.lcm([sp.denom(cf) for q in QR for cf in sp.Poly(q, u).coeffs()])
    QR = [sp.expand(q * den) for q in QR]
    g = sp.gcd_list([cf for q in QR for cf in sp.Poly(q, u).coeffs()])
    QR = [sp.expand(q / g) for q in QR]
    A = O.theta_to_D(QR, z, 2)
    print("leading:", sp.factor(A[-1]), flush=True)
    for t0_ in (0, 2, 3, sp.Rational(5, 2), sp.Rational(1, 2)):
        try:
            print("  exponents at", t0_, O.indicial(A, t0_), flush=True)
        except Exception as e:
            print("  exponents at", t0_, "err", e)
    print("  at inf", O.indicial_inf(A), flush=True)
    Lc, Rc = O.Q_to_coeffs(Q), O.Q_to_coeffs(QR)
    Mq, rem = O.rdiv(Lc, Rc)
    ok = all(sp.simplify(x) == 0 for x in rem)
    print("exact right division remainder zero:", ok, flush=True)
    (DATA / f"ring_ext_QG_6f_{tag}.json").write_text(json.dumps(
        {"Q": [str(q) for q in QR], "remainder_zero": ok, "adj": adj, "lam": lam, "r": r, "d": dmin,
         "quotient": [str(x) for x in Mq]}, indent=1))


if __name__ == "__main__":
    main()
