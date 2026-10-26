"""Step B of the minimal-order test: the constants c_i = B(Y, g_i), B the bilinear concomitant of L15.

L15 = sum_k c_k(x) theta^k (theta = theta_t, coefficients left), L15* = sum_k (-theta)^k o c_k (G.adjoint).  Lagrange:
g L(y) - L*(g) y = theta B(y, g), B(y, g) = sum_k sum_{j<k} (-1)^j theta^j(c_k g) theta^(k-1-j) y.  For y = Y (L15 Y = 0)
and g = g_i (L15* g_i = 0, pf_29_cosocle) B is constant on each interval free of singular points.  The order-one left
factor q_c of L15 with Sol(q_c*) = <g_1 + c g_2> has right cofactor M_c killing Y iff c_1 + c c_2 = 0; over Q(x) such a
factor exists iff c_1 / c_2 is rational.
g_i = x^(-a) P_i(x) s(x) / den(x), s the solution of F_2* in closed form (sympy integrate of its log-derivative, real
branch on the interval; a common constant factor cancels from c_1 / c_2).  Taylor series at x0 = t^2 to order 16 in
mpmath; theta^k Y from the 60-digit grid (pf_grid_d60_h0.03125.json, ring_ext_QG_8_onY.thetaY).
Control: g_1 with P_1 perturbed in one coefficient (not a solution of L15*): B must then vary with t.
usage: python pf_30_concomitant.py <t list>"""
from __future__ import annotations
import json, sys
from fractions import Fraction
from pathlib import Path
import mpmath as mp
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import ring_ext_QG_8_onY as E

DATA = Path(__file__).parent / "data"
NS = 17
X = sp.Symbol("x")


def ser_poly(coeffs, x0):
    """Taylor coefficients at x0 (a Fraction) of a polynomial (ascending integer coefficients), first NS, computed
    exactly (the coefficients are large and alternate) and only then rounded"""
    c = [int(v) for v in coeffs]
    num, dq = x0.numerator, x0.denominator
    out = []
    for _ in range(NS):
        if not c:
            out.append(mp.mpf(0)); continue
        acc, rest = Fraction(0), []
        for v in reversed(c):                       # synthetic division by (x - x0): value and quotient
            acc = acc * x0 + v
            rest.append(acc)
        out.append(mp.mpf(rest[-1].numerator) / rest[-1].denominator)
        c = rest[:-1][::-1]
    return out


def smul(a, b):
    return [mp.fsum(a[i] * b[n - i] for i in range(n + 1)) for n in range(NS)]


def spow(a, r):
    """(a_0 + a_1 d + ...)^r, a_0 > 0 or r integer"""
    a0 = a[0]
    out = [a0 ** r if r == int(r) else (abs(a0) ** r)]
    for n in range(1, NS):
        out.append(mp.fsum((r * k - (n - k)) * a[k] * out[n - k] for k in range(1, n + 1)) / (n * a0))
    return out


def theta(s, x0):
    d = [(n + 1) * s[n + 1] for n in range(NS - 1)] + [mp.mpf(0)]
    return [2 * (x0 * d[n] + (d[n - 1] if n else 0)) for n in range(NS)]


def s_closed(S):
    c0 = sp.Poly([int(v) for v in reversed(S[0].coeffs())], X).as_expr()
    c1 = sp.Poly([int(v) for v in reversed(S[1].coeffs())], X).as_expr()
    I = sp.integrate(sp.apart(sp.cancel(-c0 / (2 * X * c1)), X), X)
    fac = []
    for term in sp.Add.make_args(sp.expand(I)):
        cf, rest = term.as_coeff_Mul()
        lg = [a for a in sp.Mul.make_args(rest) if isinstance(a, sp.log)]
        if len(lg) != 1:
            raise ValueError(f"unexpected term {term}")
        other = sp.Mul(*[a for a in sp.Mul.make_args(rest) if not isinstance(a, sp.log)])
        fac.append((sp.Poly(lg[0].args[0], X).clear_denoms()[1], Fraction(int(sp.Rational(cf * other).p), int(sp.Rational(cf * other).q))))
    return fac


def g_series(P, den, a, fac, x0q, x0):
    out = ser_poly(P, x0q)
    out = smul(out, spow(ser_poly(den, x0q), -1))
    out = smul(out, spow(ser_poly([0, 1], x0q), -a))
    for poly, r in fac:
        out = smul(out, spow(ser_poly([int(v) for v in reversed(poly.all_coeffs())], x0q), mp.mpf(r.numerator) / r.denominator))
    return out


def concomitant(L, g, thY, x0q, x0):
    tot = []
    for k in range(1, len(L)):
        h = smul(ser_poly([int(v) for v in L[k].coeffs()], x0q), g)
        thj = h
        for j in range(k):
            tot.append((-1) ** j * thj[0] * thY[k - 1 - j])
            thj = theta(thj, x0)
    return mp.fsum(tot), max(abs(v) for v in tot)


def main():
    ts = sys.argv[1].split(",")
    rec = json.loads((DATA / "pf_grid_d60_h0.03125.json").read_text())
    mp.mp.dps = 90
    L15 = G.parse_op(json.loads((DATA / "pf_19_quotient_Q.json").read_text())["L15_lam_1/2"])
    F2 = G.parse_op(json.loads((DATA / "pf_pieces_tetra_exact.json").read_text())["factor_2"])
    cs = json.loads((DATA / "pf_29_cosocle.json").read_text())
    a = cs["box"][0]
    den = [int(v) for v in cs["den_rest"]]
    Ps = [[int(v) for v in r["P"]] for r in cs["rows"]]
    Pbad = list(Ps[0]); Pbad[len(Pbad) // 2] += 10 ** 6
    fac = s_closed(G.adjoint(F2))
    print("s = prod f^r:", [(str(f.as_expr()), str(r)) for f, r in fac], flush=True)
    out = {"t": {}}
    for t in ts:
        mp.mp.dps = 90
        x0q = Fraction(t) ** 2
        x0 = mp.mpf(x0q.numerator) / x0q.denominator
        thY = E.thetaY(rec, t, len(L15))
        row = {}
        for nm, P in (("g1", Ps[0]), ("g2", Ps[1]), ("g1 perturbed", Pbad)):
            g = g_series(P, den, a, fac, x0q, x0)
            B, big = concomitant(L15, g, thY, x0q, x0)
            row[nm] = {"B": mp.nstr(B, 40), "largest_term": mp.nstr(big, 5)}
        row["c1/c2"] = mp.nstr(mp.mpf(row["g1"]["B"]) / mp.mpf(row["g2"]["B"]), 40)
        out["t"][t] = row
        print(t, row, flush=True)
    (DATA / f"pf_30_concomitant_{ts[0]}.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
