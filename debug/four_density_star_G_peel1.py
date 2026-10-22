"""FD1-FD3 step G: peel the order-one right factors of L exactly.
Loop: algebraic solutions of the current quotient (mod p, step F), then the exact solution by exact
Frobenius series + nullspace over QQ, then theta - theta(y)/y and exact right division (remainder
must be zero).  Saves the factors and the quotient (Q-form) to four_density_star_G.json."""
from __future__ import annotations

import itertools
import json
import sys
import time
from fractions import Fraction
from pathlib import Path

import sympy as sp
from sympy.polys.matrices import DomainMatrix

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import four_density_star_F_hyperexp as F  # noqa: E402
import four_density_star_ops as O  # noqa: E402

u = O.u
t = sp.Symbol("t")
OUT = HERE / "data" / "four_density_star_G.json"


def frob_exact(Q, lam, NT):
    Qp = [sp.Poly(q, u) for q in Q]
    ev = lambda i, x: Fraction(str(Qp[i].eval(sp.Rational(x)))) if True else None
    c = [Fraction(1)]
    for n in range(1, NT // 2):
        s_ = sum(ev(i, lam + 2 * (n - i)) * c[n - i] for i in range(1, len(Q)) if n - i >= 0)
        dd = ev(0, lam + 2 * n)
        if dd == 0:
            if s_ != 0:
                return None
            c.append(Fraction(0))
        else:
            c.append(-s_ / dd)
    return {lam + 2 * n: cn for n, cn in enumerate(c) if cn}


def bser(alpha, scale, N):
    out, c = {}, Fraction(1)
    for n in range(N // 2 + 1):
        if c:
            out[2 * n] = c
        c = c * (Fraction(alpha) - n) / (n + 1) * (-Fraction(scale))
    return out


def mulx(a, b, top):
    out = {}
    for e1, c1 in a.items():
        for e2, c2 in b.items():
            if e1 + e2 < top:
                out[e1 + e2] = out.get(e1 + e2, 0) + c1 * c2
    return out


def exact_solution(Q, al, be, ga, dmax, N=160):
    roots = sorted(r for r in set(sp.roots(sp.Poly(Q[0], u).as_expr(), u)) if r.is_rational)
    sols = []
    for lam in roots:
        if lam > dmax:
            continue
        s_ = frob_exact(Q, int(lam), N)
        if s_ is not None:
            sols.append((lam, s_))
    w = mulx(mulx(bser(-Fraction(al), 1, N), bser(-Fraction(be), Fraction(1, 4), N), N + 4), bser(-Fraction(ga), Fraction(1, 9), N), N + 4)
    lo = min(l for l, _ in sols)
    top = N + int(lo) - 4
    cols = [mulx(s_, w, top) for _, s_ in sols]
    rows = [[sp.Rational(col.get(e, 0).numerator, col.get(e, 0).denominator) for col in cols]
            for e in range(int(sp.floor(dmax)) + 1, top)]
    M = DomainMatrix.from_list_sympy(len(rows), len(cols), rows).convert_to(sp.QQ)
    ns = M.nullspace().to_Matrix()
    assert ns.shape[0] == 1, ns.shape
    vec = list(ns.row(0))
    poly = 0
    for v, col in zip(vec, cols):
        for e, c in col.items():
            if e <= dmax:
                poly += v * sp.Rational(c.numerator, c.denominator) * t**e
    y = sp.factor(poly * (1 - t**2) ** sp.Rational(al) * (1 - t**2 / 4) ** sp.Rational(be) * (1 - t**2 / 9) ** sp.Rational(ga))
    return y


def main() -> None:
    rec = json.loads(OUT.read_text()) if OUT.exists() else {}
    Q = [sp.sympify(q) for q in json.loads((HERE / "data" / "four_density_star_A.json").read_text())["L"]["Q"]]
    cur = O.from_Q(Q)
    factors = []
    for step in range(6):
        t0 = time.time()
        Qc, Dc = O.to_Q(cur)
        mu = max(r for r in sp.roots(sp.Poly(Qc[-1], u).as_expr(), u) if r.is_rational)
        _, found = F.algebraic_solutions(Qc, [0, -1, -2, -3], [0, Fraction(-1, 2), Fraction(-3, 2)], [0, -1, -2], mu, N=500, verbose=False)
        print(f"step {step}: order {len(cur)-1}; algebraic-solution boxes: {found}  ({time.time()-t0:.0f}s)", flush=True)
        if not found:
            break
        al, be, ga, nul, dmax = found[0]
        y = exact_solution(Qc, Fraction(al), Fraction(be), Fraction(ga), sp.Rational(dmax))
        r = sp.factor(sp.cancel(t * sp.diff(y, t) / y))
        print(f"   y = {y};  theta y / y = {r}", flush=True)
        R = O.from_theta_t([-r, sp.Integer(1)])
        qq, rem = O.rdiv(cur, R)
        assert all(c == 0 for c in rem), "nonzero remainder"
        factors.append({"y": str(y), "r": str(r)})
        cur = qq
    Qc, Dc = O.to_Q(cur)
    rec["order1_factors_bottom_up"] = factors
    rec["quotient_Q"] = [str(q) for q in Qc]
    rec["quotient_order"] = len(cur) - 1
    print("quotient order", len(cur) - 1, " leading:", sp.factor(sum(sp.Poly(q, u).coeff_monomial(u ** (len(cur) - 1)) * t ** (2 * i) for i, q in enumerate(Qc))), flush=True)
    OUT.write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
