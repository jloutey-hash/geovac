"""CF record (d68e6fe5), after scoring CF1-CF2.  Not registered; diagnostics.
 (a) the t-twisted L_Leg equals the EC source operator exactly; the CF2 map equals the EC map up to scale.
 (b) why CF1 failed: L_Leg on the cosh part C(t) = int K0(tx) K0(x) cosh(g x) dx (t > g - 1) and on the
     sinh part S(t); H = C - S.
 (c) the quotient Q3 with L5 = Q3 o L_Leg (route 2), and the values phi = L_Leg H.
Usage: python debug/c2b_feyn_C_diag.py"""
import json
import sys
from pathlib import Path

import mpmath as mp
import sympy as sp

sys.path.insert(0, str(Path(__file__).parent))
from c2b_feyn_A import leg_op, mellin_route, theta_to_d, compose_D, t  # noqa: E402

x = sp.Symbol("x")
DATA = Path(__file__).parent / "data"
ROOT = Path(__file__).parent.parent


def part_vals(tv, g, kind):
    w = {"cosh": lambda y: mp.cosh(g * y), "sinh": lambda y: mp.sinh(g * y), "exp": lambda y: mp.exp(-g * y)}[kind]
    f0 = lambda y: mp.besselk(0, tv * y) * mp.besselk(0, y) * w(y)
    f1 = lambda y: -y * mp.besselk(1, tv * y) * mp.besselk(0, y) * w(y)
    f2 = lambda y: y**2 * (mp.besselk(0, tv * y) + mp.besselk(1, tv * y) / (tv * y)) * mp.besselk(0, y) * w(y)
    return [mp.quad(f, [0, 1, 4, 16, 64, mp.inf]) for f in (f0, f1, f2)]


def apply(cs, vals, tv):
    terms = [sp.lambdify(t, c, "mpmath")(tv) * v for c, v in zip(cs, vals)]
    return sum(terms), max(abs(v) for v in terms)


def quotient(L, R):
    L = [sp.cancel(v) for v in L]
    m = len(R) - 1
    Q = {}
    while len(L) - 1 >= m and any(v != 0 for v in L):
        k = len(L) - 1
        c = sp.cancel(L[-1] / R[-1])
        Q[k - m] = c
        CR = compose_D(R, k - m)
        L = [sp.cancel(L[i] - c * (CR[i] if i < len(CR) else 0)) for i in range(len(L))]
        while L and L[-1] == 0:
            L.pop()
    return [Q.get(i, 0) for i in range(max(Q) + 1)], L


if __name__ == "__main__":
    out = {}
    B = json.load(open(DATA / "c2b_feyn_B_map.json"))
    src = [sp.sympify(v) for v in B["source"]]
    ec = [-352 * x**2 + 380 * x + 125, 256 * x**4 + 640 * x**3 - 2480 * x**2,
          256 * x**5 - 1344 * x**4 - 1680 * x**3 + 500 * x**2]
    ratio = [sp.cancel(a / b) for a, b in zip(src, ec)]
    out["twisted_LLeg_equals_EC_source"] = len(set(ratio)) == 1
    print("twisted L_Leg = EC source up to scale:", out["twisted_LLeg_equals_EC_source"], ratio)
    T = json.load(open(ROOT / "tests" / "data" / "paper59_c2b_2026_10_22.json"))
    run = B["runs"][-1]
    n0c, n0e = sp.sympify(run["N0"]), sp.sympify(T["N0"].replace("^", "**"))
    n1c, n1e = sp.sympify(run["N1"]), sp.sympify(T["N1"].replace("^", "**"))
    r0, r1 = sp.cancel(n0c / n0e), sp.cancel(n1c / n1e)
    Ec, Ee = sp.sympify(run["E"]), sp.sympify(T["E"].replace("^", "**"))
    out["map_equals_EC_map"] = r0 == r1 and sp.cancel(Ec / Ee).is_number
    print("CF2 map = EC map up to scale:", out["map_equals_EC_map"], r0, r1, sp.cancel(Ec / Ee))

    mp.mp.dps = 40
    g = mp.mpf(3) / 2
    L = leg_op(sp.Rational(3, 2))
    diag = {}
    for tv in [mp.mpf("0.8"), mp.mpf("1.7"), mp.mpf("3.2")]:
        row = {}
        for kind in ("cosh", "sinh", "exp"):
            r, big = apply(L, part_vals(tv, g, kind), tv)
            row[kind] = float(abs(r) / big)
        diag[str(tv)] = row
        print("t =", tv, "relative L_Leg residual:", row, flush=True)
    out["parts"] = diag

    L5 = theta_to_d(mellin_route(sp.Rational(3, 2)))
    Q, rem = quotient(L5, L)
    out["Q3"] = [str(sp.factor(c)) for c in Q]
    out["Q3_remainder_zero"] = len(rem) == 0
    print("Q3 (L5 = Q3 o L_Leg):", [sp.factor(c) for c in Q])
    json.dump(out, open(DATA / "c2b_feyn_C_diag.json", "w"), indent=1)
