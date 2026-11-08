"""FD1-FD3 step L: exact pieces found so far, and the next peel.
  L1 = (theta-3)^3 - t^2 theta(theta-3)(theta-2) = E2 o (theta-3)   (the one-leaf moment's operator,
       mellin_t_ode on chi with weight x^-1), a right factor of L (checked here again);
  E2 = (theta-3)^2 - t^2 theta(theta-2): the map f -> t^3 (1 + theta) f sends every solution of
       Legendre's equation theta^2 - t^2 (theta+1)^2 (F_2's operator) to a solution of E2 (exact);
  R4 = X1 o L1 (step H; exact), X1 of order one;
  M15 = L / R4; its algebraic solution (step K, box al=be=ga=0, degree <= 7) made exact and peeled:
  M15 = M14 o (theta - r).
Saves Q-forms."""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction
from pathlib import Path

import sympy as sp

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import four_density_star_G_peel1 as G  # noqa: E402
import four_density_star_I_w9 as I  # noqa: E402
import four_density_star_ops as O  # noqa: E402

u = O.u
t = sp.Symbol("t")
OUT = HERE / "data" / "four_density_star_L.json"


def main() -> None:
    rec = {}
    t0 = time.time()
    Q = [sp.sympify(q) for q in json.loads((HERE / "data" / "four_density_star_A.json").read_text())["L"]["Q"]]
    L = O.from_Q(Q)
    L1 = O.from_Q([(u - 3) ** 3, -u * (u - 3) * (u - 2)])
    E2 = O.from_Q([(u - 3) ** 2, -u * (u - 2)])
    TH3 = O.from_Q([u - 3])
    rec["L1_is_E2_TH3"] = O.compose(E2, TH3) == L1
    _, rem = O.rdiv(L, L1)
    rec["L1_right_factor_of_L"] = all(c == 0 for c in rem)
    # E2 vs Legendre: the map A = t^3 (1 + theta_t) = w^(3/2) (1 + 2 delta); conjugate by w^(3/2):
    # E2' = w^(-3/2) E2 w^(3/2) acts as E2 with delta -> delta + 3/2; check E2' o (1 + 2 delta) = 0 mod Leg
    Leg = O.from_Q([u**2, -(u + 1) ** 2])
    E2c = [O.K(0)] * len(E2)
    sh = sp.Rational(3, 2)
    for k, c in enumerate(E2):          # c (delta + 3/2)^k
        for j in range(k + 1):
            E2c[j] += c * sp.binomial(k, j) * sh ** (k - j)
    A = [O.K(1), O.K(2)]
    comp = O.compose(E2c, A)
    _, rem = O.rdiv(comp, Leg)
    rec["E2_Legendre_map_exact"] = all(c == 0 for c in rem)
    print("L1 = E2 o (theta-3):", rec["L1_is_E2_TH3"], " L1 right factor of L:", rec["L1_right_factor_of_L"],
          " t^3(1+theta) maps Sol(Legendre) into Sol(E2):", rec["E2_Legendre_map_exact"], flush=True)
    # R4 = X1 o L1
    R4Q = I.primitive([sp.sympify(q) for q in json.loads((HERE / "data" / "four_density_star_H.json").read_text())["R4_Q"]])
    R4 = O.from_Q(R4Q)
    X1, rem = O.rdiv(R4, L1)
    rec["R4_is_X1_L1"] = all(c == 0 for c in rem)
    rec["X1"] = [str(c.as_expr()) for c in X1]
    print("R4 = X1 o L1:", rec["R4_is_X1_L1"], " X1 =", [sp.factor(c.as_expr()) for c in X1], flush=True)
    M15, rem = O.rdiv(L, R4)
    assert all(c == 0 for c in rem)
    M15Q, _ = O.to_Q(M15)
    y = G.exact_solution(M15Q, Fraction(0), Fraction(0), Fraction(0), sp.Integer(7), N=200)
    r = sp.factor(sp.cancel(t * sp.diff(y, t) / y))
    print(f"M15 algebraic solution y = {y}; theta y / y = {r}  ({time.time()-t0:.0f}s)", flush=True)
    R = O.from_theta_t([-r, sp.Integer(1)])
    M14, rem = O.rdiv(M15, R)
    rec["M15_peel"] = {"y": str(y), "r": str(r), "zero_remainder": all(c == 0 for c in rem)}
    print("M15 = M14 o (theta - r) exact:", rec["M15_peel"]["zero_remainder"], flush=True)
    M14Q, D = O.to_Q(M14)
    M14Q = I.primitive(M14Q)
    rec["M14_Q"] = [str(q) for q in M14Q]
    rec["M14_cleared_den"] = str(sp.factor(D))
    print(f"M14: order {len(M14)-1}, {len(M14Q)} t-terms, cleared denominator {sp.factor(D)}", flush=True)
    OUT.write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
