"""FP step L: local exponents (indicial polynomials) of R3 (step K), M8 and W12 at t = 1, 2, 3 and infinity."""
from __future__ import annotations
import json, sys
from pathlib import Path
import sympy as sp
HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import four_density_star_E_exps as E  # noqa: E402
u, t, rho = E.u, E.t, E.rho
OUT = HERE / "data" / "four_density_path_L.json"


def exps(Q, name, pts=(1, 2, 3)):
    Q = [sp.sympify(q).subs(sp.Symbol("u"), u) for q in Q]
    Lc = E.theta_coeffs(Q)
    out = {"t=0": str(sp.factor(Q[0].subs(u, rho))), "inf": str(sp.factor(Q[-1].subs(u, rho)))}
    for tt in pts:
        ip, lo = E.indicial(Lc, sp.Integer(tt))
        out[f"t={tt}"] = str(ip)
    print(name, out, flush=True)
    return out


def main():
    rec = {}
    K = json.loads((HERE / "data" / "four_density_path_K.json").read_text())
    rec["R3"] = exps(K["R3_Q"], "R3")
    OUT.write_text(json.dumps(rec, indent=1))
    if "all" in sys.argv:
        rec["M8"] = exps(K["M8_Q"], "M8")
        OUT.write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
