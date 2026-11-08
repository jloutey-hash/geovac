"""PF minimality: Y's Taylor data (of S = Y/t^3 in x) at the monodromy base point B = 0.12 + 0.12i, carried from the
grid at t_s through the upper half-plane by the stepping of pf_11_monodromy.py (init copied from ib_tetra_pf_sheets.py).
The span of the monodromy orbit of this vector is the solution space of Y's minimal annihilator.
usage: python pf_11c_y_at_base.py <grid json> <t_s> <dps> <out tag>"""
from __future__ import annotations
import json, sys, time
from pathlib import Path
import mpmath as mp
sys.path.insert(0, str(Path(__file__).parent))
import ring_ext_QG_8_onY as E
import pf_11_monodromy as P
from sympy.functions.combinatorial.numbers import stirling as _st


def init(grid, tstr, R):
    rec = json.loads(Path(grid).read_text())
    t = mp.mpf(tstr)
    x = t ** 2
    th = E.thetaY(rec, tstr, R - 1)
    thx = []
    for j in range(R):
        cf = [mp.binomial(j, i) * mp.mpf(-3) ** (j - i) / mp.mpf(2) ** j for i in range(j + 1)]
        thx.append(x ** mp.mpf(-1.5) * mp.fsum(cf[i] * th[i] for i in range(j + 1)))
    a = [mp.mpc(mp.fsum(mp.mpf(int(_st(k, j, kind=1, signed=True))) * thx[j] for j in range(k + 1))
                / x ** k / mp.factorial(k)) for k in range(R)]
    return a, x


def main():
    grid, ts, dps, tag = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
    P.setup(dps)
    t0 = time.time()
    R = P.R
    a, xs = init(grid, ts, R)
    XB, H = P.ST["XB"], P.ST["H"]
    a = P.path(a, [mp.mpc(xs), mp.mpc(xs, H), mp.mpc(XB, H)])
    rec = {"grid": grid, "t_s": ts, "dps": dps, "base": "0.12+0.12i", "Y_at_B": P.enc(a), "seconds": time.time() - t0}
    (P.DATA / f"pf_11c_y_at_base_{tag}.json").write_text(json.dumps(rec))
    print(f"saved ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
