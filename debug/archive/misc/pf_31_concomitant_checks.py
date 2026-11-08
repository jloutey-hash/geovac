"""Checks on pf_30_concomitant.
(1) Code check at full precision: s_1, the closed-form solution of the order-one RIGHT factor F_1 of L15 (exact
    division), is an exact solution of L15 in a class other than the 9/4 class of g_1, g_2, so B(s_1, g_i) = 0
    exactly; computed with theta^k s_1 from its exact series it must vanish to the working precision relative to the
    largest term, at every t.
(2) Precision scaling on Y: B(Y, g_i) with theta^k Y from the h grid and from the 2h grid (every second node).
    A real value is unchanged; a value at the noise floor (largest term x quadrature error) grows with the error.
usage: python pf_31_concomitant_checks.py <t list>"""
from __future__ import annotations
import json, sys
from fractions import Fraction
from pathlib import Path
import mpmath as mp
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import ring_ext_QG_8_onY as E
import pf_30_concomitant as C

DATA = Path(__file__).parent / "data"


def theta_values(ser, x0, J):
    out, cur = [], ser
    for _ in range(J + 1):
        out.append(cur[0]); cur = C.theta(cur, x0)
    return out


def main():
    ts = sys.argv[1].split(",")
    rec = json.loads((DATA / "pf_grid_d60_h0.03125.json").read_text())
    mp.mp.dps = 90
    L15 = G.parse_op(json.loads((DATA / "pf_19_quotient_Q.json").read_text())["L15_lam_1/2"])
    pc = json.loads((DATA / "pf_pieces_tetra_exact.json").read_text())
    F1, F2 = G.parse_op(pc["factor_1"]), G.parse_op(pc["factor_2"])
    cs = json.loads((DATA / "pf_29_cosocle.json").read_text())
    a, den = cs["box"][0], [int(v) for v in cs["den_rest"]]
    Ps = [[int(v) for v in r["P"]] for r in cs["rows"]]
    facg = C.s_closed(G.adjoint(F2))
    fac1 = C.s_closed(F1)
    out = {}
    for t in ts:
        x0q = Fraction(t) ** 2
        x0 = mp.mpf(x0q.numerator) / x0q.denominator
        s1 = [mp.mpf(1)] + [mp.mpf(0)] * (C.NS - 1)
        for poly, r in fac1:
            s1 = C.smul(s1, C.spow(C.ser_poly([int(v) for v in reversed(poly.all_coeffs())], x0q),
                                   mp.mpf(r.numerator) / r.denominator))
        th_s1 = theta_values(s1, x0, len(L15) - 1)
        thY = E.thetaY(rec, t, len(L15))
        thY2 = E.thetaY(rec, t, len(L15), every=2)
        row = {}
        for nm, P in (("g1", Ps[0]), ("g2", Ps[1])):
            g = C.g_series(P, den, a, facg, x0q, x0)
            b1, m1 = C.concomitant(L15, g, th_s1, x0q, x0)
            bh, mh = C.concomitant(L15, g, thY, x0q, x0)
            b2, m2 = C.concomitant(L15, g, thY2, x0q, x0)
            row[nm] = {"B(s1)/largest": mp.nstr(abs(b1) / m1, 3),
                       "B(Y,h)/largest": mp.nstr(abs(bh) / mh, 3), "B(Y,2h)/largest": mp.nstr(abs(b2) / m2, 3),
                       "B(Y,h)": mp.nstr(bh, 8), "B(Y,2h)": mp.nstr(b2, 8)}
        row["Y quad h vs 2h"] = mp.nstr(abs(thY2[0] / thY[0] - 1), 3)
        out[t] = row
        print(t, row, flush=True)
    (DATA / f"pf_31_concomitant_checks_{ts[0]}.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
