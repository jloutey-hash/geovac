"""FP step H: the operator of the far-field piece S_a(t) = t^3 int K0(t x) chi(x) dx (H -> H(oo); Mellin route,
alpha = -2, beta = 3, v = chi's operator) and whether it is a right factor of W12 (step D) / W11 = W12/(theta-2)."""
from __future__ import annotations
import json, sys
from pathlib import Path
import sympy as sp
HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import density_loop_B_ode as B  # noqa: E402
import four_density_star_ops as O  # noqa: E402
import four_density_path_F_struct as F  # noqa: E402
u, z = O.u, B.z


def main():
    v = [sp.Integer(0), z - 3, sp.Integer(-4) * 0 + 4]  # placeholder
    v = [0, (z - 3) * 0, 0]
    # chi's operator as v-list: theta^3 - 4 theta^2 - (z-3) theta = ((theta-2)^2 - 1 - z) theta
    v = [sp.Integer(0), -(z - 3), sp.Integer(-4), sp.Integer(1)]
    Qa = B.mellin_t_ode(v, -2, 3, lambda sg: sg**2)
    print("S_a operator:", [sp.factor(q) for q in Qa], flush=True)
    W = O.from_Q(F.load())
    Ra = O.from_Q(Qa)
    _, rem = O.rdiv(W, Ra)
    print("S_a operator right factor of W12:", all(c == 0 for c in rem), flush=True)


if __name__ == "__main__":
    main()
