"""Data for tests/test_paper59_offslice_y_operator.py (written as a separate step from the PF results).

Stores the exact order-19 tetrahedron operator for S = Y/t^3 in x = t^2 (D-form, sum_k polys[k](x) (d/dx)^k S = 0,
integer coefficients, from data/ib_tetra_pf_op_Q.json) and theta_t^k Y, k = 0..19, at t = 0.8, 1.7, 2.6, 3.8 from
the 60-digit grid (data/pf_grid_d60_h0.03125.json), 60 significant digits.

Before writing, checks that the D-form converts to the stored theta-form Q (sum_i t^(2i) Q_i(theta_t) Y) up to a
common x-power and scalar, so the test may rebuild the theta form from the D form alone.
usage: python build_offslice_y_guard_data_2026_11_01.py"""
from __future__ import annotations

import json
import sys
from fractions import Fraction as Fr
from pathlib import Path

import mpmath as mp
import sympy as sp

sys.path.insert(0, str(Path(__file__).parent))
import ring_ext_QG_8_onY as E  # noqa: E402

DATA = Path(__file__).parent / "data"
OUT = Path(__file__).parent.parent / "tests" / "data" / "paper59_offslice_y_2026_11_01.json"
U = sp.Symbol("u")


def d_to_theta(polys):
    """x^19 * x^(3/2) * sum_k p_k(x) D^k (x^(-3/2) Y) = sum_k p_k(x) x^(19-k) prod_{j<k} (u/2 - 3/2 - j) Y."""
    R = len(polys) - 1
    out = {}
    for k, p in enumerate(polys):
        P = sp.Integer(1)
        for j in range(k):
            P *= (U / 2 - sp.Rational(3, 2) - j)
        P = sp.Poly(sp.expand(P), U)
        for i, c in enumerate(p):
            if c:
                e = i + R - k
                out[e] = out.get(e, sp.Poly(0, U)) + P * c
    return out


def main():
    d = json.loads((DATA / "ib_tetra_pf_op_Q.json").read_text())
    polys = [[int(v) for v in p] for p in d["polys"]]
    conv = d_to_theta(polys)
    Qs = [sp.Poly(sp.sympify(q), U) for q in d["Q"]]
    lo = min(conv)
    keys = sorted(conv)
    assert len(keys) == len(Qs), (len(keys), len(Qs))
    # common scalar
    ratio = None
    for i, e in enumerate(keys):
        a, b = conv[e], Qs[i]
        if b.is_zero:
            assert a.is_zero
            continue
        r = Fr(int(sp.Poly(a, U).LC().p), int(sp.Poly(a, U).LC().q)) / Fr(int(b.LC()))
        if ratio is None:
            ratio = r
        assert a == b * sp.Rational(ratio.numerator, ratio.denominator), e
    print("D-form -> theta form matches stored Q; x shift", lo, "scalar", ratio)

    grid = json.loads((DATA / "pf_grid_d60_h0.03125.json").read_text())
    rec = {"polys": [[str(v) for v in p] for p in polys],
           "cuts_x": ["1/4", "1", "9/4", "4", "9", "49/4", "16", "81/4"],
           "ordinary_x": ["3", "25/4", "25"],
           "thetaY": {}}
    mp.mp.dps = 68
    for t in ("0.8", "1.7", "2.6", "3.8"):
        th = E.thetaY(grid, t, 19)
        rec["thetaY"][t] = [mp.nstr(v, 60) for v in th[:20]]
        print("t", t, "Y", mp.nstr(th[0], 25), flush=True)
    OUT.write_text(json.dumps(rec))
    print(OUT, OUT.stat().st_size)


if __name__ == "__main__":
    main()
