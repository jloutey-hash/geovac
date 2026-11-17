"""Scoring driver for record P63L1 (Paper 63, tree fibre at l = 1): the eps^-3/2 branch at t = -2.

Object (debug/p63l1_engine.py): chain of three, slice (t;1,1), s-type root of scale t, two p_z leaves
rho_p = d/dz [K0(r)/(2 pi^2)], Coulomb kernels.  Fibre = t^3 int K0(t r) H(r) dr, H = (2/(3pi)) r^2 g^2,
g = (2/pi)[K1 - chi/r^2] = -1/r^2 + E.  At t = -(2 - eps) + i0 the branch is pi (2-eps)^3 S(eps) with
      S(eps) = int_1^oo I0((2 - eps) x) H_sing(x) dx,   H_sing = (2/(3pi)) x^2 E(x)^2   (e^-2x piece of H).
Closed form (debug/p63l1_controls.py, C3, from the tail E ~ sqrt(2/pi) r^-1/2 e^-r):
      eps^(3/2) S(eps) -> A = 1/(3 pi^2),  eps^(3/2) S = A (1 + kappa eps) + O(eps^(3/2)).

Computes, for eps = 1e-2 .. 1e-4 (9 log-spaced points):
  S at two precisions: float64 (log-grid trapezoid, scaled Bessels, SA's scheme) and mpmath (dps 30,
  composite Gauss-Legendre); their relative difference;
  the fitted log-slope d log S / d log eps (all points, and the last 5), local slopes;
  a Richardson-type extrapolation of eps^(3/2) S in u = sqrt(eps): least squares on {1, u^2, u^3, u^4}
  (the exponent pattern of the expansion: no u^1 term) over all 9 points and over the smallest 6;
  the ratio of the extrapolated amplitude to the closed form A, and of the fitted u^2 coefficient to A kappa.

Usage: python debug/p63l1_run.py [--object l1|s] [--out PATH]
  --object s is the toy (the s-type chain of three: A = 1/pi^2, SA's n = 3 branch times (2/pi)^3).
Default output: debug/data/p63l1_run.json
"""
import argparse
import json
import sys
import time
from pathlib import Path

import mpmath as mp
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import p63l1_engine as E  # noqa: E402

CLOSED = {
    # exact expressions derived in p63l1_controls.py C3 (sympy), from the leaf-potential tails only
    "l1": {"A": 1 / (3 * np.pi**2), "A_expr": "1/(3 pi^2)", "kappa": 47 / 8},
    "s": {"A": 1 / np.pi**2, "A_expr": "1/pi^2", "kappa": 31 / 8},
}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--object", choices=["l1", "s"], default="l1")
    ap.add_argument("--out", default=str(HERE / "data" / "p63l1_run.json"))
    ap.add_argument("--dps", type=int, default=30)
    ap.add_argument("--ngl", type=int, default=40)
    args = ap.parse_args()

    hf = E.h_l1_f if args.object == "l1" else E.h_s_f
    hm = E.h_l1_mp if args.object == "l1" else E.h_s_mp
    A, kappa = CLOSED[args.object]["A"], CLOSED[args.object]["kappa"]
    mp.mp.dps = args.dps
    eps_list = list(np.logspace(-2, -4, 9))
    rows = []
    t0 = time.time()
    for eps in eps_list:
        sf = E.branch_float(eps, hf)
        sm = E.branch_mp(mp.mpf(eps), hm, ngl=args.ngl)
        sm_f = float(sm)
        rows.append({"eps": eps, "S_float64": sf, "S_mp": mp.nstr(sm, 25), "rel_float_vs_mp": (sf - sm_f) / sm_f,
                     "eps32_S_mp": eps**1.5 * sm_f, "ratio_eps32_S_to_A": eps**1.5 * sm_f / A})
        print(f"eps={eps:.3e}  S_mp={mp.nstr(sm, 15)}  float rel {(sf-sm_f)/sm_f:.1e}  "
              f"eps^1.5 S / A = {eps**1.5*sm_f/A:.6f}  [{time.time()-t0:.0f}s]", flush=True)

    e = np.array(eps_list)
    S = np.array([float(mp.mpf(r["S_mp"])) for r in rows])
    le, lS = np.log(e), np.log(S)
    slope_all = float(np.polyfit(le, lS, 1)[0])
    slope_last5 = float(np.polyfit(le[-5:], lS[-5:], 1)[0])
    local = [float((lS[i + 1] - lS[i]) / (le[i + 1] - le[i])) for i in range(len(e) - 1)]

    u = np.sqrt(e)
    y = e**1.5 * S

    def fit(sel):
        M = np.stack([u[sel] ** 0, u[sel] ** 2, u[sel] ** 3, u[sel] ** 4], 1)
        c, *_ = np.linalg.lstsq(M, y[sel], rcond=None)
        return c

    c_all = fit(slice(None))
    c_small = fit(slice(3, None))
    out = {
        "object": args.object,
        "closed_form": {"A": A, "A_expr": CLOSED[args.object]["A_expr"], "kappa": kappa},
        "precision": {"mp_dps": args.dps, "gauss_legendre_nodes": args.ngl, "float_npts": 400001},
        "rows": rows,
        "log_slope_fit_all": slope_all,
        "log_slope_fit_last5": slope_last5,
        "local_log_slopes": local,
        "richardson_all9": {"A_fit": float(c_all[0]), "ratio_to_closed": float(c_all[0] / A),
                            "u2_coeff_over_A_kappa": float(c_all[1] / (A * kappa))},
        "richardson_smallest6": {"A_fit": float(c_small[0]), "ratio_to_closed": float(c_small[0] / A),
                                 "u2_coeff_over_A_kappa": float(c_small[1] / (A * kappa))},
        "max_rel_float_vs_mp": float(max(abs(r["rel_float_vs_mp"]) for r in rows)),
        "runtime_s": time.time() - t0,
    }
    print(f"log-slope all {slope_all:.6f}, last5 {slope_last5:.6f}; local {np.round(local, 5)}")
    print(f"Richardson A (all 9) = {c_all[0]:.10g}, ratio {c_all[0]/A:.8f}; (smallest 6) ratio {c_small[0]/A:.8f}")
    print(f"u^2 coefficient / (A kappa): all {c_all[1]/(A*kappa):.5f}, smallest 6 {c_small[1]/(A*kappa):.5f}")
    Path(args.out).write_text(json.dumps(out, indent=1))
    print("wrote", args.out)


if __name__ == "__main__":
    main()
