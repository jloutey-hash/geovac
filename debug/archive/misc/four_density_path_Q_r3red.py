"""FP step Q: is R3 (step K) reducible?  Order-1 right factors = algebraic solutions of R3; order-1 left factors =
algebraic solutions of the adjoint R3* (Q*_i(u) = Q_i(-u - 2i)).  Box of four_density_star_F_hyperexp with the
exponents allowed by R3's local data (t = 1: 0, 1, -3; t = 2: 0, 1, -3/2) and their adjoint counterparts.
Also: the Mellin operators of the odd-weight piece t^3 int K0(t x) x^3 K0(x)^2 dx and of F_3 (weight x^0)."""
from __future__ import annotations
import json, sys
from fractions import Fraction
from pathlib import Path
import sympy as sp
HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import four_density_star_F_hyperexp as F  # noqa: E402
import density_loop_B_ode as Bm  # noqa: E402
u, z = F.u, Bm.z


def main():
    K = json.loads((HERE / "data" / "four_density_path_K.json").read_text())
    R3 = [sp.sympify(q).subs(sp.Symbol("u"), u) for q in K["R3_Q"]]
    R3s = [sp.expand(q.subs(u, -u - 2 * i)) for i, q in enumerate(R3)]
    for name, Q in [("R3", R3), ("R3*", R3s)]:
        top = [r for r in sp.roots(sp.Poly(Q[-1], u).as_expr(), u) if r.is_rational]
        mu = max(top)
        print(name, "exponents t=0:", sp.factor(Q[0]), " top:", sp.factor(Q[-1]), flush=True)
        e1 = [0, -1, -2, -3, 1, 2, 3]
        e2 = [0, Fraction(-1, 2), Fraction(-3, 2), Fraction(1, 2), Fraction(3, 2), -1, 1]
        sols, found = F.algebraic_solutions(Q, e1, e2, [0], mu + 12, N=300, verbose=False)
        print(f"  {name}: log-free exponents {sols}; boxes with algebraic solutions: {found}", flush=True)
    vK2 = Bm.sym_power([z, 0], 2, 2)
    Qodd = Bm.mellin_t_ode(vK2, 1, 3, lambda sg: sg**2)
    QF3 = Bm.mellin_t_ode(vK2, 1, 0, lambda sg: sg**2)
    print("odd-weight piece t^3 int K0(tx) x^3 K0^2:", [sp.factor(q) for q in Qodd], flush=True)
    print("F_3:", [sp.factor(q) for q in QF3], flush=True)


if __name__ == "__main__":
    main()
