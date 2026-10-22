"""FP step S: R3 = L1 o R2 exactly (R2 = annihilator of R3's exponent-3 log-free solution, order 2, t-degree 6),
R2's irreducibility (exponents at 0 are 1 (log forced) and 3; a hyperexponential solution would be log-free at 0,
hence y3, whose annihilator has order 2), and a map from Legendre's operator (F_2) into R2."""
from __future__ import annotations
import json, sys, time
from pathlib import Path
import sympy as sp
HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import four_density_star_G_peel1 as Gp  # noqa: E402
import four_density_star_H_order4 as H  # noqa: E402
import four_density_star_I_w9 as I9  # noqa: E402
import four_density_star_ops as O  # noqa: E402
import four_density_path_gauge as G  # noqa: E402
u, t = O.u, sp.Symbol("t")
OUT = HERE / "data" / "four_density_path_S.json"


def main():
    rec = {}
    K = json.loads((HERE / "data" / "four_density_path_K.json").read_text())
    R3Q = [sp.sympify(q).subs(sp.Symbol("u"), u) for q in K["R3_Q"]]
    y3 = Gp.frob_exact(R3Q, 3, 2 * 60)
    sols = H.guess(y3, 2, 3, 50)
    print("R2 nullity (order 2, t-degree 6):", len(sols), flush=True)
    R2Q = I9.primitive(sols[0])
    while R2Q[-1] == 0:
        R2Q.pop()
    print("R2:", [sp.factor(q) for q in R2Q], flush=True)
    R3 = O.from_Q(R3Q)
    R2 = O.from_Q(R2Q)
    L1, rem = O.rdiv(R3, R2)
    rec["R3_div_R2_zero"] = all(c == 0 for c in rem)
    rec["R2_Q"] = [str(q) for q in R2Q]
    rec["L1_delta_w"] = [str(c.as_expr()) for c in L1]
    print("R3 = L1 o R2:", rec["R3_div_R2_zero"], " L1 =", [sp.factor(c.as_expr()) for c in L1], flush=True)
    Leg = [u**2, -(u + 1) ** 2]
    for par, (lo, hi) in {"even": (-12, 20), "odd": (-11, 21)}.items():
        n = G.nullity(R2Q, Leg, 0, 2, {"a": 4, "b": 4}, lo, hi, NT=300)
        print(f"Legendre -> R2 {par}: nullity {n}", flush=True)
        rec[f"Leg_to_R2_{par}"] = n
    OUT.write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
