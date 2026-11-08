"""FD2/FD3 step P: the right factors of L on S, by quadrature (theta_t^j S cached by step D at 45 digits).
  B5 = (theta - r) o X1 o E2 o (theta - 3)   (below F_4's piece; FD2: B5 S != 0 at two t)
  W9 = R4' o B5                              (FD3 candidate; registered residual <= 1e-25 at four t,
                                              perturbation of one coefficient of R4' by 1/200 > 1e-5)
also the lower right factors theta - 3, L1 = E2 o (theta-3), R4 = X1 o L1, and L itself."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import mpmath as mp
import sympy as sp

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import four_density_star_O_eval as EV  # noqa: E402
import four_density_star_ops as O  # noqa: E402
import four_density_star_I_w9 as I  # noqa: E402

u = O.u
t = sp.Symbol("t")
OUT = HERE / "data" / "four_density_star_P.json"


def build():
    Q = [sp.sympify(q) for q in json.loads((HERE / "data" / "four_density_star_A.json").read_text())["L"]["Q"]]
    L = O.from_Q(Q)
    TH3 = O.from_Q([u - 3])
    E2 = O.from_Q([(u - 3) ** 2, -u * (u - 2)])
    L1 = O.compose(E2, TH3)
    R4 = O.from_Q(I.primitive([sp.sympify(q) for q in json.loads((HERE / "data" / "four_density_star_H.json").read_text())["R4_Q"]]))
    Ld = json.loads((HERE / "data" / "four_density_star_L.json").read_text())
    r = sp.sympify(Ld["M15_peel"]["r"]).subs(sp.Symbol("t"), t)
    THr = O.from_theta_t([-r, sp.Integer(1)])
    B5 = O.compose(THr, R4)
    Nd = json.loads((HERE / "data" / "four_density_star_N.json").read_text())
    R4p = O.from_Q([sp.sympify(q) for q in Nd["R4p_Q"]])
    W9 = O.compose(R4p, B5)
    return dict(L=L, TH3=TH3, L1=L1, R4=R4, B5=B5, R4p=R4p, W9=W9)


def perturbed(R4p, k=None, factor=sp.Rational(201, 200)):
    """change one coefficient: the theta^k coefficient of R4' (function-left, monic top) scaled."""
    R = [c / R4p[-1] for c in R4p]
    k = 2 if k is None else k
    R[k] = R[k] * factor
    return R


def main() -> None:
    ops = build()
    rec = {}
    _, rem = O.rdiv(ops["L"], ops["W9"])
    rec["L_div_W9_zero"] = all(c == 0 for c in rem)
    print("L = N o W9 exact:", rec["L_div_W9_zero"], flush=True)
    d = json.loads(EV.DFILE.read_text())
    ts = [k.split("@")[0] for k in d if k.endswith("@45")]
    Wp = O.compose(perturbed(ops["R4p"]), ops["B5"])
    names = ["TH3", "L1", "R4", "B5", "W9", "L"]
    for tv in ts:
        row = {}
        for nm in names:
            tot, big, rel = EV.evaluate(ops[nm], tv, 45)
            row[nm] = {"value": mp.nstr(tot, 12), "largest": mp.nstr(big, 6), "rel": mp.nstr(rel, 4)}
        tot, big, rel = EV.evaluate(Wp, tv, 45)
        row["W9_perturbed"] = {"value": mp.nstr(tot, 12), "largest": mp.nstr(big, 6), "rel": mp.nstr(rel, 4)}
        rec[tv] = row
        print(f"t={tv}: " + "  ".join(f"{k}: {v['rel']}" for k, v in row.items()), flush=True)
    OUT.write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
