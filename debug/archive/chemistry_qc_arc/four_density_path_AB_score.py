"""FP scoring: residuals on P (cached theta_t^j P at 32 and 60 digits) of
  W12 = Y4 o Z2 o Z2' o L1 o R2 o (theta - 2)   (FP3 candidate: <= 1e-25 at 32 digits at four t)
  W12p = Y4p o R8, Y4p = Y4 (monic) with its delta^2 coefficient times 201/200 (perturbation fixed before evaluation,
         same choice as the star record: theta^2 coefficient of the F_4 piece)          (> 1e-5 required)
  R8 = Z2 o Z2' o L1 o R2 o (theta - 2), the right factor below the F_4 piece          (FP2: > 1e-5 at two t)
  and the lower right factors theta - 2, R2 o (theta-2), R3 o (theta-2), Z2' o R3 o (theta-2)."""
import json, sys
from pathlib import Path
import mpmath as mp
import sympy as sp
HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import four_density_star_ops as O  # noqa: E402
import four_density_path_M_eval as ME  # noqa: E402
u = O.u
OUT = HERE / "data" / "four_density_path_AB.json"


def Qop(lst):
    return O.from_Q([sp.sympify(q).subs(sp.Symbol("u"), u) for q in lst])


def main():
    V = json.loads((HERE / "data" / "four_density_path_V.json").read_text())
    K = json.loads((HERE / "data" / "four_density_path_K.json").read_text())
    S = json.loads((HERE / "data" / "four_density_path_S.json").read_text())
    X = json.loads((HERE / "data" / "four_density_path_X.json").read_text())
    D = json.loads((HERE / "data" / "four_density_path_D.json").read_text())
    W12 = Qop(D["W12"]["Q"])
    Y4, R8 = Qop(V["Y4_Q"]), Qop(V["R8_Q"])
    TH2 = O.from_Q([u - 2])
    R2 = Qop(S["R2_Q"])
    R3 = Qop(K["R3_Q"])
    Z2p = Qop(X["Z2p_Q"])
    ops = {"theta-2": TH2, "R2o(theta-2)": O.compose(R2, TH2), "R3o(theta-2)": O.compose(R3, TH2),
           "Z2'oR3o(theta-2)": O.compose(Z2p, O.compose(R3, TH2)), "R8": R8, "W12": W12}
    # identity check: W12 = Y4 o R8 (monic comparison)
    comp = O.compose(Y4, R8)
    ok = all(a * W12[-1] == b * comp[-1] for a, b in zip(comp, W12))
    Ym = [c / Y4[-1] for c in Y4]
    Ym[2] = Ym[2] * sp.Rational(201, 200)
    ops["W12_perturbed"] = O.compose(Ym, R8)
    rec = {"W12_eq_Y4oR8": ok}
    for dps in (32, 60):
        for tv in ME.avail(dps):
            row = {}
            for nm, op in ops.items():
                tot, big, rel = ME.evaluate(op, tv, dps)
                row[nm] = mp.nstr(rel, 4)
            rec[f"{tv}@{dps}"] = row
            print(f"t={tv} @{dps}:", row, flush=True)
    OUT.write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
