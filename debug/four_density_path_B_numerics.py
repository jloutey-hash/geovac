"""FP step B: normalization N = 16/pi^4 and the momentum-space control (registered: P at one t agrees with an
independent momentum-space evaluation to 1e-8 or better); grid convergence of P.  Also the control that
the same position code with H -> chi (leaf d removed: H_0 = int min(r,r') r' K0(r') dr' = chi) gives the chain
C(t) = (8/pi^3) t^3 int chi^2 K0(t x) dx at its registered value (density_loop_B_ode.C_pos)."""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import mpmath as mp

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import four_density_path_num as PN  # noqa: E402

OUT = HERE / "data" / "four_density_path_B.json"


def main() -> None:
    rec = json.loads(OUT.read_text()) if OUT.exists() else {}
    t0 = time.time()
    # grid convergence at 40 digits
    G1 = PN.Grid(40, 32)
    G2 = PN.Grid(40, 64)
    for tv in ["0.7", "1.6"]:
        p1 = G1.thetas(tv, 0)[0]
        p2 = G2.thetas(tv, 0)[0]
        rec[f"P_t{tv}"] = {"hinv32": mp.nstr(p1, 35), "hinv64": mp.nstr(p2, 35), "rel": mp.nstr(abs(p1 / p2 - 1), 3)}
        print(f"P({tv}) hinv32 {mp.nstr(p1, 30)}  hinv64 {mp.nstr(p2, 30)}  rel {mp.nstr(abs(p1/p2-1), 3)}", flush=True)
    rec["Hinf"] = mp.nstr(G2.Hinf, 35)
    # H(x) check: H'' = -K0 chi at a point, by finite differences of the quadrature definition
    mp.mp.dps = 30
    xx = mp.mpf(2)
    F = lambda r: mp.besselk(0, r) * PN.chi(r)
    Hq = lambda x: mp.quad(lambda r: r * F(r), [0, x]) + x * mp.quad(F, [x, mp.inf])
    d2 = mp.diff(Hq, xx, 2)
    rec["H_second_derivative_check"] = mp.nstr(abs(d2 + F(xx)), 3)
    print("H'' + K0 chi at x=2:", rec["H_second_derivative_check"], flush=True)
    # leaf removed: P with H -> chi, (8/pi^3) normalization -> chain
    import density_loop_B_ode as B
    mp.mp.dps = 40
    c_grid = G2.thetas("0.7", 0, which="chi")[0] * (mp.pi / 2)    # N=16/pi^4 -> 8/pi^3: factor pi/2
    mp.mp.dps = 25
    c_ref = B.C_pos("0.7")
    rec["chain_leaf_removed"] = {"grid": mp.nstr(c_grid, 22), "C_pos": mp.nstr(c_ref, 22), "rel": mp.nstr(abs(c_grid / c_ref - 1), 3)}
    print("leaf removed (H -> chi):", rec["chain_leaf_removed"], flush=True)
    OUT.write_text(json.dumps(rec, indent=1))
    # momentum control at 20 digits
    mp.mp.dps = 20
    for tv in ["0.7", "1.6"]:
        t1 = time.time()
        pm = PN.P_mom(tv)
        pp = mp.mpf(rec[f"P_t{tv}"]["hinv64"])
        rec[f"mom_t{tv}"] = {"P_mom": mp.nstr(pm, 20), "P_pos": mp.nstr(pp, 20), "rel": mp.nstr(abs(pm / pp - 1), 3)}
        print(f"momentum t={tv}:", rec[f"mom_t{tv}"], f"({time.time()-t1:.0f}s)", flush=True)
        OUT.write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
