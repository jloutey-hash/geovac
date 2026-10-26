"""Data for tests/test_paper59_y_minimal_order.py (written as a separate step from the results of pf_27-pf_34).
Operators as ascending integer coefficient lists of c_k(x) (L = sum_k c_k(x) theta_t^k, x = t^2): M13 (pf_34), R_S
(pf_21), R_V (pf_10_rvmap), F_1, F_2, F_5 (pf_pieces_tetra_exact), z_1, z_2 (pf_12_hom.legendre_source); theta^k Y,
k = 0..15, at t = 1.7 and 2.6 from the 60-digit grid (h), 60 significant digits.
usage: python pf_35_guard_data.py"""
from __future__ import annotations
import json, sys
from pathlib import Path
import mpmath as mp
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import pf_lib  # noqa: F401
import pf_12_hom as P
import ring_ext_QG_8_onY as E

DATA = Path(__file__).parent / "data"
OUT = Path(__file__).parent.parent / "tests" / "data" / "paper59_y_2026_10_26.json"


def ints(L):
    return [[str(int(v)) for v in c.coeffs()] for c in L]


def main():
    pc = json.loads((DATA / "pf_pieces_tetra_exact.json").read_text())
    rec = {
        "M13": ints(G.parse_op(json.loads((DATA / "pf_34_m13_exact.json").read_text())["M13"])),
        "R_S": ints(G.parse_op(json.loads((DATA / "pf_21_pieces.json").read_text())["R_S"]["op"])),
        "R_V": ints(__import__("pf_10_rvmap").RV),
        "F_1": ints(G.parse_op(pc["factor_1"])), "F_2": ints(G.parse_op(pc["factor_2"])),
        "F_5": ints(G.parse_op(pc["factor_5"])),
        "z1": ints(P.legendre_source("z1")), "z2": ints(P.legendre_source("z2")),
        "cuts_x": ["1/4", "1", "9/4", "4", "9", "49/4", "16", "81/4"],
        "thetaY": {},
    }
    grid = json.loads((DATA / "pf_grid_d60_h0.03125.json").read_text())
    mp.mp.dps = 68
    for t in ("1.7", "2.6"):
        th = E.thetaY(grid, t, 15)
        rec["thetaY"][t] = [mp.nstr(v, 60) for v in th[:16]]
    OUT.write_text(json.dumps(rec))
    print(OUT, OUT.stat().st_size)


if __name__ == "__main__":
    main()
