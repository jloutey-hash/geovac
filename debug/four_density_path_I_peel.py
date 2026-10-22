"""FP step I: W11 = W12/(theta - 2) exactly, then peel order-one right factors (algebraic solutions, mod p
search in the box (1-t^2)^-al (1-t^2/4)^-be (1-t^2/9)^-ga of four_density_star_F_hyperexp, made exact by
four_density_star_G_peel1.exact_solution, exact right division)."""
from __future__ import annotations
import json, sys, time
from fractions import Fraction
from pathlib import Path
import sympy as sp
HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import four_density_star_F_hyperexp as F  # noqa: E402
import four_density_star_G_peel1 as G  # noqa: E402
import four_density_star_ops as O  # noqa: E402
import four_density_path_F_struct as FS  # noqa: E402
u, t = O.u, sp.Symbol("t")
OUT = HERE / "data" / "four_density_path_I.json"


def main():
    rec = {}
    W = O.from_Q(FS.load())
    cur, rem = O.rdiv(W, O.from_Q([u - 2]))
    assert all(c == 0 for c in rem)
    factors = []
    for step in range(6):
        t0 = time.time()
        Qc, Dc = O.to_Q(cur)
        roots_top = [r for r in sp.roots(sp.Poly(Qc[-1], u).as_expr(), u) if r.is_rational]
        mu = max(roots_top)
        exps0 = sorted(set(r for r in sp.roots(sp.Poly(Qc[0], u).as_expr(), u)))
        print(f"step {step}: order {len(cur)-1}, exponents at 0 {exps0}, max exponent at oo {mu}", flush=True)
        _, found = F.algebraic_solutions(Qc, [0, -1, -2, -3, -4], [0, Fraction(-1, 2), -1, Fraction(-3, 2), -2], [0, -1, -2], mu, N=500, verbose=False)
        print(f"   boxes with algebraic solutions: {found} ({time.time()-t0:.0f}s)", flush=True)
        if not found:
            break
        al, be, ga, nul, dmax = found[0]
        y = G.exact_solution(Qc, Fraction(al), Fraction(be), Fraction(ga), sp.Rational(dmax))
        r = sp.factor(sp.cancel(t * sp.diff(y, t) / y))
        print(f"   y = {y}; theta y / y = {r}", flush=True)
        R = O.from_theta_t([-r, sp.Integer(1)])
        qq, rem = O.rdiv(cur, R)
        assert all(c == 0 for c in rem), "nonzero remainder"
        factors.append({"y": str(y), "r": str(r)})
        cur = qq
    Qc, Dc = O.to_Q(cur)
    rec["order1_factors_bottom_up_after_theta-2"] = factors
    rec["quotient_Q"] = [str(q) for q in Qc]
    rec["quotient_order"] = len(cur) - 1
    print("quotient order", len(cur) - 1, flush=True)
    OUT.write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
