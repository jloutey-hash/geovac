"""P3e-3 step I: does G(t;1,1) need N2?  L_G = N2 R3 R1 (theta+2) exactly (steps B, F, H).
Evaluate W G for W = R3 o R1 o (theta+2) (order 5) by quadrature with theta_t^j applied under the
integral (step A's evaluator).  W G = 0 at every t means N2 is not needed: G is annihilated by K3
(R3, equivalent to F_3's operator) above two elementary pieces.
Control: the same with one coefficient of R3 changed by 1/198 must not vanish.  Also: L_G G at
the same precision (the operator itself), and V(t;1,1)-free sanity: (R1 (theta+2)) G != 0.
  python debug/p3e2w_slice_I_needs.py <dps>"""
from __future__ import annotations
import json, sys, time
from pathlib import Path
import mpmath as mp
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import density_loop_D10_k3_exact_check as D10  # noqa: E402
import density_loop_D11_chain_minimal as D11  # noqa: E402
import p3e2w_slice_A_operator as A  # noqa: E402

t, u = D10.t, D10.u
HERE = Path(__file__).parent


def main():
    dps = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    H = json.loads((HERE / "data" / "p3e2w_slice_H.json").read_text())
    R3 = [sp.sympify(x) for x in H["R3"]]
    R1 = [-(1 + 7 * t**2), 1 - t**2]                   # (1 - t^2) theta - (1 + 7 t^2)
    E2 = D10.compose(R1, [sp.Integer(2), sp.Integer(1)])  # R1 o (theta + 2)
    W = D10.compose(R3, E2)
    QW, den = D11.to_Q(W)
    R3p = [R3[0], R3[1], sp.cancel(R3[2] * sp.Rational(199, 198)), R3[3]]
    QWp, _ = D11.to_Q(D10.compose(R3p, E2))
    QE2, _ = D11.to_Q(E2)
    QG = [sp.sympify(q) for q in json.loads((HERE / "data" / "p3e2w_slice_A.json").read_text())["QG"]]
    print("W order", len(W) - 1, " cleared denominator", sp.factor(den), flush=True)
    mp.mp.dps = dps
    rec = {"dps": dps, "QW": [str(q) for q in QW], "rows": {}}
    for ts in ["0.5", "0.9", "1.3", "2.5", "3.5"]:
        t0 = time.time()
        g = A.LGnum([sp.Integer(1)], ts, None)
        wg = A.LGnum(QW, ts, None)
        cp = A.LGnum(QWp, ts, None)
        e2 = A.LGnum(QE2, ts, None)
        lg = A.LGnum(QG, ts, None)
        rec["rows"][ts] = {"G": mp.nstr(g, 20), "WG": mp.nstr(wg, 5), "control": mp.nstr(cp, 5),
                           "E2G": mp.nstr(e2, 8), "LG": mp.nstr(lg, 5)}
        print(f"t {ts}: G = {mp.nstr(g, 15)}  (den W) G = {mp.nstr(wg, 5)}  control = {mp.nstr(cp, 5)}  "
              f"R1(th+2)G = {mp.nstr(e2, 6)}  L_G G = {mp.nstr(lg, 5)}  ({time.time()-t0:.0f}s)", flush=True)
    (HERE / "data" / f"p3e2w_slice_I_{dps}.json").write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
