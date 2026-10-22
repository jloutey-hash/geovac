"""FP step AA: the composite exact map from F_4's operator Q4 (eq:four_density_op) into Y4 (the order-4 left factor
of W12): A = (C o G) mod Q4 (right-division remainder), C from step U (R8 P = C P_h4), G = t^-3 (theta-4)^2 (theta-2)^2 theta^2
(P_h4 = (16/pi^4) G f4, f4 = int K0(t x) K0(x)^3 dx).  Exact checks: Y4 o A = 0 mod Q4; image annihilator = Y4."""
import json, sys, time
from pathlib import Path
import sympy as sp
HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "tests"))
import test_paper59_generic_slice as GS  # noqa: E402
import four_density_star_ops as O  # noqa: E402
import four_density_path_gauge as G  # noqa: E402
t, u = GS.t, GS.u
OUT = HERE / "data" / "four_density_path_AA.json"


def main():
    t0 = time.time()
    V = json.loads((HERE / "data" / "four_density_path_V.json").read_text())
    CQ = [sp.sympify(q, locals={"u": u}) for q in V["C_Q"]]
    while CQ and CQ[-1] == 0:
        CQ.pop()
    Cth = GS._coeff_list(CQ)
    Gth = [sp.expand(c * t**-3) for c in GS._coeff_list([sp.expand((u - 4)**2 * (u - 2)**2 * u**2)])]
    Q4 = [9 * u**2 * (u - 2) ** 2, -2 * u**2 * (5 * u**2 + 1), (u + 1) ** 4]
    Q4th = GS._coeff_list(Q4)
    CG = GS._compose(Cth, Gth)
    _, A = GS._rdiv(CG, Q4th)
    A = [sp.factor(sp.cancel(a)) for a in A]
    print("A (theta-form, order", len(A) - 1, "):", f"({time.time()-t0:.0f}s)", flush=True)
    for a in A:
        print("   ", a, flush=True)
    Y4Q = [sp.sympify(q).subs(sp.Symbol("u"), G.u) for q in V["Y4_Q"]]
    Q4g = [q.subs(u, G.u) for q in Q4]
    Ag = [a.subs(t, G.t) for a in A]
    ok, R, odd = G.exact_checks(Y4Q, Q4g, Ag)
    Y4 = O.from_Q(Y4Q)
    monic = [c / Y4[-1] for c in Y4]
    same = R is not None and all(x == y for x, y in zip(R, monic))
    print("Y4 o A = 0 mod Q4:", ok, " injective, image annihilator = Y4:", same, " odd:", odd, f"({time.time()-t0:.0f}s)", flush=True)
    OUT.write_text(json.dumps({"A_theta_form": [str(a) for a in A], "Y4oA_mod_Q4_zero": ok, "image_ann_equals_Y4": same, "odd": odd}, indent=1))


if __name__ == "__main__":
    main()
