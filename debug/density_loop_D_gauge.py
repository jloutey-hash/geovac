"""Loop test PT2, step D: which known families the triangle and the chain control contain.

Series in t, mod p = 2^61 - 1.  A D-module relation q(t) y = sum_{j<r} p_j(t) theta^j f between a
solution y of one operator and a solution f of a reference operator (polynomials q, p_j of degree
<= D, q != 0) shows that y lies in the Q(t)-span of f and its derivatives.

References: f2 = 2F1(1/2,1/2;1;t^2) (F_2's Legendre equation, the two-density elliptic rung),
            f3 = 3F2(1/2,1/2,1/2;1,1;t^2/4) (F_3's operator 4 th^3 - t^2 (th+1)^3, the chain's K3 rung).
Triangle: the order-1 right factor of L_G (from its regular solution at t = 0) is found exactly;
y = L1 h, h the t^3 Frobenius solution of L_G, is a solution of the order-2 left factor M2.
Chain: every non-resonant Frobenius solution of its order-7 operator; minimal annihilator order;
relation to f3 and f2.
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).parent))
import density_loop_C_factor as C  # noqa: E402

P = C.P
u = sp.Symbol("u")
OUT = Path(__file__).parent / "data" / "density_loop_D.json"
NT = 300   # number of t-exponent slots kept


def inv(a):
    return pow(a % P, P - 2, P)


def fr(x):
    x = Fraction(x)
    return x.numerator % P * inv(x.denominator) % P


def frob_t(Q, lam_t, NT=NT):
    """Solution t^lam (sum c_n t^{2n}) of sum_i t^{2i} Q_i(theta_t); returned as dict t-exp -> coeff."""
    Qp = [sp.Poly(q, u) for q in Q]

    def ev(i, x):
        val = Fraction(0)
        for (m,), c in zip(Qp[i].monoms(), Qp[i].coeffs()):
            val += Fraction(int(sp.Rational(c).p), int(sp.Rational(c).q)) * Fraction(x) ** m
        return fr(val)
    c = [1]
    for n in range(1, NT // 2):
        s_ = sum(ev(i, lam_t + 2 * (n - i)) * c[n - i] for i in range(1, len(Q)) if n - i >= 0) % P
        d = ev(0, lam_t + 2 * n)
        if d == 0:
            raise ZeroDivisionError(f"resonance at n={n}")
        c.append((-s_) * inv(d) % P)
    return {lam_t + 2 * n: cn for n, cn in enumerate(c)}


def hyper_t(a_list, b_list, NT=NT, scale=Fraction(1)):
    """pFq(a; b; scale * t^2) as a t-series."""
    out, term = {}, Fraction(1)
    for n in range(NT // 2):
        out[2 * n] = fr(term)
        num = Fraction(1)
        for a in a_list:
            num *= (Fraction(a) + n)
        den = Fraction(n + 1)
        for b in b_list:
            den *= (Fraction(b) + n)
        term = term * num / den * scale
    return out


def theta(ser, j=1):
    return {e: c * pow(e % P, j, P) % P for e, c in ser.items()}


def shift(ser, d):
    return {e + d: c for e, c in ser.items()}


def relation_count(y, fs, D, lo=0, hi=None):
    """Dimension of {(q, p_j)} with q y = sum_j p_j f_j, deg <= D, restricted to t-exponents < hi.
    Returns (dimension, dimension with q forced to 0)."""
    exps = set()
    for s_ in [y] + fs:
        exps |= set(s_)
    m0 = min(exps)
    hi = hi if hi is not None else m0 + NT - 2 * D - 4
    cols = []
    for d in range(D + 1):
        cols.append(shift(y, d))
    for f in fs:
        for d in range(D + 1):
            cols.append({e: (-c) % P for e, c in shift(f, d).items()})
    rows_idx = list(range(m0, hi))
    M = [[col.get(e, 0) for col in cols] for e in rows_idx]
    rank = C.rank_mod(M)
    Mq0 = [row[D + 1:] for row in M]
    rank0 = C.rank_mod(Mq0)
    return len(cols) - rank, (len(cols) - (D + 1)) - rank0


def minimal_order(Q, ser, maxorder, D):
    """Smallest r <= maxorder such that ser has a nonzero annihilator sum_{j<=r,d<=D} a t^d theta^j."""
    for r in range(1, maxorder + 1):
        fs = [theta(ser, j) for j in range(r + 1)]
        cols = [shift(f, d) for f in fs for d in range(D + 1)]
        exps = sorted(set().union(*cols))
        rows = [[col.get(e, 0) for col in cols] for e in exps[: len(cols) + 40] if e < min(ser) + NT - 2 * D - 4]
        if len(cols) - C.rank_mod(rows) > 0:
            return r
    return None


def main() -> None:
    rec = {}
    B = json.loads((Path(__file__).parent / "data" / "density_loop_B.json").read_text())
    QG = [sp.sympify(q) for q in B["QG"]]
    QC = [sp.sympify(q) for q in B["QC"]]
    f2 = hyper_t([Fraction(1, 2)] * 2, [1])
    f3 = hyper_t([Fraction(1, 2)] * 3, [1, 1], scale=Fraction(1, 4))
    F2s = [f2, theta(f2)]
    F3s = [f3, theta(f3), theta(f3, 2)]
    # sanity: f3 vs f2 should be unrelated at order < 3 (control: f3 related to itself)
    print("control f3~f3:", relation_count(f3, F3s, 4), " f3~f2:", relation_count(f3, F2s, 12), flush=True)

    # ---- triangle
    g0 = frob_t(QG, 0)
    h = frob_t(QG, 3)
    print("triangle: minimal order of g0:", minimal_order(QG, g0, 3, 6), " of h:", minimal_order(QG, h, 3, 12), flush=True)
    # exact order-1 factor for g0: rational g0'/g0 -- find (a0 + a1 t^2) theta g0 + (b0 + b1 t^2) g0 = 0
    from itertools import product as _p
    # solve over Q with exact Fraction series (first 12 terms)
    t = sp.Symbol("t")
    Qp = [sp.Poly(q, u) for q in QG]
    cex = [Fraction(1)]
    for n in range(1, 12):
        s_ = sum(Fraction(str(Qp[i].as_expr().subs(u, 2 * (n - i)))) * cex[n - i] for i in range(1, 3) if n - i >= 0)
        cex.append(-s_ / Fraction(str(Qp[0].as_expr().subs(u, 2 * n))))
    g0ex = sum(sp.Rational(c.numerator, c.denominator) * t ** (2 * n) for n, c in enumerate(cex))
    a0, a1, b0, b1 = sp.symbols("a0 a1 b0 b1")
    expr = sp.expand((a0 + a1 * t**2) * t * sp.diff(g0ex, t) + (b0 + b1 * t**2) * g0ex)
    eqs = [expr.coeff(t, k) for k in range(0, 22, 2)]
    Mx = sp.Matrix([[sp.diff(e, v) for v in (a0, a1, b0, b1)] for e in eqs])
    ns = Mx.nullspace()
    print("order-1 factor for g0: nullspace dim", len(ns), [list(v) for v in ns], flush=True)
    rec["g0_order1"] = [str(list(v)) for v in ns]
    sol = ns
    if len(ns) == 1:
        vv = ns[0]
        A = sp.expand(vv[0] + vv[1] * t**2)
        Bb = sp.expand(vv[2] + vv[3] * t**2)
        gclosed = sp.simplify(sp.exp(-sp.integrate(sp.simplify(Bb / A) / t, t)))
        print("   g0 closed form:", sp.factor(gclosed), " check:", sp.simplify(sp.series(gclosed, t, 0, 10).removeO() - sp.series(g0ex, t, 0, 10).removeO()), flush=True)
        rec["g0_closed"] = str(sp.factor(gclosed))
        # y = L1 h = A theta h + Bb h, as t-series mod p
        Ap = sp.Poly(sp.expand(A), t)
        Bp = sp.Poly(sp.expand(Bb), t)
        y = {}
        th_h = theta(h)
        for (k,), cf in zip(Ap.monoms(), Ap.coeffs()):
            for e, c in th_h.items():
                y[e + k] = (y.get(e + k, 0) + fr(Fraction(str(cf))) * c) % P
        for (k,), cf in zip(Bp.monoms(), Bp.coeffs()):
            for e, c in h.items():
                y[e + k] = (y.get(e + k, 0) + fr(Fraction(str(cf))) * c) % P
        y = {e: c for e, c in y.items() if c}
        print("   minimal order of y = L1 h:", minimal_order(None, y, 2, 12), flush=True)
        for D in (6, 12, 20):
            r2, r2q0 = relation_count(y, F2s, D)
            r3, r3q0 = relation_count(y, F3s, D)
            print(f"   y vs f2 (deg {D}): relations {r2} (with q=0: {r2q0});  y vs f3: {r3} (q=0: {r3q0})", flush=True)
            rec[f"tri_y_vs_f2_D{D}"] = [r2, r2q0]
            rec[f"tri_y_vs_f3_D{D}"] = [r3, r3q0]

    # ---- chain control
    Q0 = sp.Poly(QC[0], u)
    roots = sorted(set(int(r) for r in sp.roots(Q0.as_expr(), u)))
    print("chain exponents at t=0:", sp.factor(QC[0]), flush=True)
    for lam in roots:
        try:
            yc = frob_t(QC, lam)
        except ZeroDivisionError as e:
            print(f"  lam {lam}: resonant ({e})", flush=True)
            continue
        mo = minimal_order(QC, yc, 7, 10)
        r3, r3q0 = relation_count(yc, F3s, 16)
        r2, r2q0 = relation_count(yc, F2s, 16)
        print(f"  lam {lam}: minimal order {mo}; vs f3 {r3} (q=0 {r3q0}); vs f2 {r2} (q=0 {r2q0})", flush=True)
        rec[f"chain_lam{lam}"] = {"minimal_order": mo, "vs_f3": [r3, r3q0], "vs_f2": [r2, r2q0]}
    OUT.write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
