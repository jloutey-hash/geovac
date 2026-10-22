"""FP step X: X4 = Z2 o Z2' exactly.  N2 = annihilator of the exponent -1 log-free solution of the adjoint X4*
(order 2), X4* = M o N2 exactly; Z2 = N2* (adjoint, a left factor of X4); Z2' = X4 left-divided by Z2 (exact).
Adjoint w.r.t. dw/w in delta form (four_density_star_ops.adjoint)."""
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
u, t = O.u, sp.Symbol("t")
OUT = HERE / "data" / "four_density_path_X.json"


def main():
    rec = {}
    t0 = time.time()
    V = json.loads((HERE / "data" / "four_density_path_V.json").read_text())
    X4Q = [sp.sympify(q).subs(sp.Symbol("u"), u) for q in V["X4_Q"]]
    X4sQ = [sp.expand(q.subs(u, -u - 2 * i)) for i, q in enumerate(X4Q)]
    z = Gp.frob_exact(X4sQ, -1, 2 * 220)
    sols = H.guess(z, 2, 38, 200)
    print("N2 nullity (order 2, t-degree 76):", len(sols), f"({time.time()-t0:.0f}s)", flush=True)
    N2Q = I9.primitive(sols[0])
    while N2Q[-1] == 0:
        N2Q.pop()
    X4 = O.from_Q(X4Q)
    X4s = O.adjoint(X4)
    N2 = O.from_Q(N2Q)
    # N2 is in Q-form for the operator in theta_t; the adjoint in delta form must match X4s (same convention)
    Mq, rem = O.rdiv(X4s, N2)
    rec["X4adj_div_N2_zero"] = all(c == 0 for c in rem)
    print("X4* = M o N2 exact:", rec["X4adj_div_N2_zero"], f"({time.time()-t0:.0f}s)", flush=True)
    Z2 = O.adjoint(N2)
    Z2p, rem = O.ldiv(X4, Z2)
    rec["X4_eq_Z2_o_Z2p"] = all(c == 0 for c in rem)
    chk = O.compose(Z2, Z2p)
    rec["recompose_equals_X4"] = all(a * X4[-1] == b * chk[-1] for a, b in zip(chk, X4)) if len(chk) == len(X4) else False
    print("X4 = Z2 o Z2' (left division remainder zero):", rec["X4_eq_Z2_o_Z2p"], " recompose:", rec["recompose_equals_X4"], flush=True)
    rec["Z2_Q"] = [str(q) for q in O.to_Q(Z2)[0]]
    rec["Z2p_Q"] = [str(q) for q in O.to_Q(Z2p)[0]]
    for nm in ["Z2_Q", "Z2p_Q"]:
        Q = [sp.sympify(q) for q in rec[nm]]
        print(nm, "t-terms", len(Q), " t^0:", sp.factor(Q[0]), " top:", sp.factor(Q[-1]), flush=True)
        print("   lead:", sp.factor(sum(sp.Poly(q, u).coeff_monomial(u**2) * t**(2*i) for i, q in enumerate(Q))), flush=True)
    OUT.write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
