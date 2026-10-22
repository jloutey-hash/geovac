"""FD3 step S: the order-8 right factor W8 = R4'' o R4 (step R; F_4's piece directly on R4, skipping the
order-1 piece theta - r) on S, and its perturbation control; exact checks L = N o W8 and W9 = Z o W8."""
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
import four_density_star_P_score as PS  # noqa: E402

OUT = HERE / "data" / "four_density_star_S.json"


def main() -> None:
    ops = PS.build()
    Rd = json.loads((HERE / "data" / "four_density_star_R.json").read_text())
    R4pp = O.from_Q([sp.sympify(q) for q in Rd["R4p_Q"]])
    W8 = O.compose(R4pp, ops["R4"])
    rec = {"map_box": Rd["box"], "M_A_mod_Q4_zero": Rd["M14_A_mod_Q4_zero"], "div_zero": Rd["M14_div_R4p_zero"]}
    _, rem = O.rdiv(ops["L"], W8)
    rec["L_div_W8_zero"] = all(c == 0 for c in rem)
    Z, rem2 = O.rdiv(ops["W9"], W8)
    rec["W9_div_W8_zero"] = all(c == 0 for c in rem2)
    print("L = N o W8:", rec["L_div_W8_zero"], "  W9 = Z o W8:", rec["W9_div_W8_zero"], " ord Z", len(Z) - 1, flush=True)
    Wp = O.compose(PS.perturbed(R4pp), ops["R4"])
    d = json.loads(EV.DFILE.read_text())
    for tv in [k.split("@")[0] for k in d if k.endswith("@45")]:
        a = EV.evaluate(W8, tv, 45)
        b = EV.evaluate(Wp, tv, 45)
        rec[tv] = {"W8_rel": mp.nstr(a[2], 4), "W8_perturbed_rel": mp.nstr(b[2], 4), "largest": mp.nstr(a[1], 6)}
        print(f"t={tv}: W8 S rel {mp.nstr(a[2], 4)}   perturbed {mp.nstr(b[2], 4)}", flush=True)
    OUT.write_text(json.dumps(rec, indent=1, default=str))


if __name__ == "__main__":
    main()
