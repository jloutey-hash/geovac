"""PF1-PF4: the Mellin image of the kite's operator (kite-first route, record QG step 2a).  From the exact kite operator
sum_i z^i P_i(theta_z) annihilating phi(z), ring_ext_QG_ops.mellin_Y_from_phi gives an exact operator M_Y annihilating
Y(t) = C t^3 int z^(1/2) (z + t^2)^(-3/2) phi(z) dz, hence Y_log.  Checks: (i) M_Y on the float Y_log series (exact
recurrence residual on all coefficients); (ii) the guessed minimal annihilator L of Y_log right-divides M_Y exactly
(pseudo-division over Q), i.e. L is a right factor of the kite-first upper bound.

Object: the exponential-geminal correlated triangle Y(t) at the fingerprint fibre with vertex 1 = density rho_i of
scale a = 1, vertex 2 = density rho_j of scale b = 2, vertex 3 = apex density rho_k of scale t; geminal
exp(-gamma_1 r_13) with gamma_1 = 3/2, geminal exp(-gamma_2 r_23) with gamma_2 = 3, Coulomb 1/r_12 (QG conventions).
usage: python pf_mellin.py <kite op json (Q in theta_z)> <series json> <N> <dps> <out tag> [<Y_log op json (Q)>]"""
from __future__ import annotations
import json, sys, time
from pathlib import Path
sys.set_int_max_str_digits(10**7)
import sympy as sp
from flint import arb, ctx
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import ring_ext_QG_ops as O
import ring_ext_QG_guess as Gs

DATA = Path(__file__).parent / "data"
u = O.u


def main():
    kf, sf, N, dps, tag = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
    yf = sys.argv[6] if len(sys.argv) > 6 else None
    t0 = time.time()
    Qphi = [sp.sympify(q, locals={"u": u}) for q in json.loads(Path(kf).read_text())["Q"]]
    QY = O.mellin_Y_from_phi(Qphi)
    LY = G.from_Q(QY)
    print(f"M_Y: order {len(LY)-1}, x-degree {G.op_degree(LY)}  ({time.time()-t0:.0f}s)", flush=True)
    rec = {"kite_op": kf, "M_Y_order": len(LY) - 1, "M_Y_x_degree": G.op_degree(LY), "M_Y_Q": [str(q) for q in QY]}
    # (i) on the float Y_log series
    ctx.dps = dps
    phi = [arb(s) for s in json.loads(Path(sf).read_text())["phi_n"]][:N]
    c, b = [], arb(1)
    for n in range(N):
        c.append(b * phi[n])
        b = b * (arb(-3) / 2 - n) / (n + 1)
    P = {}
    for i, q in enumerate(QY):
        for (j,), v in sp.Poly(q, u).terms():
            P[(i, j)] = int(v)
    ex = Gs.exact_residuals(P, c, 3, 2, N)
    rec["M_Y_on_Ylog_series_residual"] = ex.mid().str(3, radius=False)
    print("M_Y on the Y_log float series (all %d coefficients): %s" % (N, rec["M_Y_on_Ylog_series_residual"]), flush=True)
    if yf:
        L = G.from_Q([sp.sympify(q, locals={"u": u}) for q in json.loads(Path(yf).read_text())["Q"]])
        m, Qt, rem = G.pseudo_rdiv(LY, L)
        ok = G.is_zero_op(rem)
        rec["L_right_divides_M_Y_exactly"] = bool(ok)
        rec["quotient_order"] = len(Qt) - 1 if ok else None
        print(f"the guessed Y_log operator (order {len(L)-1}) right-divides M_Y exactly: {ok}; quotient order {len(Qt)-1}  ({time.time()-t0:.0f}s)", flush=True)
    (DATA / f"pf_mellin_{tag}.json").write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
