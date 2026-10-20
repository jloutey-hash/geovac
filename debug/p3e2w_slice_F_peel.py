"""P3e-3 step F: peel R1 = (theta-1) - t^2 (theta+7) (solution t/(1-t^2)^4, step E) off M6:
M6 = M5 o R1 exactly (density_loop_D7.rdiv).  Then, on M5: mod-p factor tests against f3, f2,
f2(t^2/4) (sub R / quotient L, rational gauge, D5 machinery), local exponents, algebraic
solutions of M5 and of its adjoint (= order-1 left factors)."""
from __future__ import annotations
import json, sys
from fractions import Fraction
from pathlib import Path
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import density_loop_D_gauge as D  # noqa: E402
import density_loop_D4_factors as F4  # noqa: E402
import density_loop_D5_rational as F5  # noqa: E402
import density_loop_D7_m5 as F7  # noqa: E402
import p3e2w_slice_D_exponents as DX  # noqa: E402
import p3e2w_slice_E_hyperexp as EX  # noqa: E402

u, rho = sp.symbols("u rho")
HERE = Path(__file__).parent


def exps(Q):
    Lc = DX.theta_coeffs(Q)
    e = {"0": sp.factor(Q[0].subs(u, rho)), "inf": sp.factor(Q[-1].subs(u, rho))}
    for t0 in (1, 2):
        e[str(t0)] = DX.indicial_at(Lc, t0)[0]
    return e


def roots_of(expr):
    return [r for r in sp.roots(sp.Poly(expr, rho)).keys()]


def factor_tests(name, Q, rec):
    NT = F4.NT
    refs = [("f3", D.hyper_t([Fraction(1, 2)] * 3, [1, 1], NT=NT, scale=Fraction(1, 4)), 3),
            ("f2", D.hyper_t([Fraction(1, 2)] * 2, [1], NT=NT), 2),
            ("f2(t^2/4)", D.hyper_t([Fraction(1, 2)] * 2, [1], NT=NT, scale=Fraction(1, 4)), 2)]
    for ref, f, r in refs:
        for side, QQ in [("R", Q), ("L", F4.adjoint(Q))]:
            hit, where = 0, None
            for b in range(4):
                for c in range(4):
                    hit = F5.test(QQ, f, r, 10, b, c)
                    if hit:
                        where = (b, c)
                        break
                if hit:
                    break
            rec[f"{name}|{ref}|{side}"] = [hit, where]
            print(f"  {name} {ref:10s} ({side}): {'FOUND nullity %d at (t^2-1)^%d (t^2-4)^%d' % (hit, *where) if hit else 'none'}", flush=True)


def main():
    M6 = [sp.sympify(q) for q in json.loads((HERE / "data" / "p3e2w_slice_B.json").read_text())["M6"]]
    M5 = F7.rdiv(M6, [u - 1, -(u + 7)])
    print("M5 = M6 / R1:", flush=True)
    for i, q in enumerate(M5):
        print(f"  t^{2*i}:", sp.factor(q), flush=True)
    rec = {"M5": [str(sp.expand(q)) for q in M5]}
    e = exps(M5)
    print("  exponents:", {k: str(v) for k, v in e.items()}, flush=True)
    rec["M5_exps"] = {k: str(v) for k, v in e.items()}
    factor_tests("M5", M5, rec)
    # algebraic solutions of M5 (right factors) and of M5* (left factors)
    for nm, Q in [("M5", M5), ("M5*", F4.adjoint(M5))]:
        ee = exps(Q)
        e1 = [r for r in roots_of(ee["1"]) if r.is_rational]
        e2 = [r for r in roots_of(ee["2"]) if r.is_rational]
        ei = [r for r in roots_of(ee["inf"]) if r.is_rational]
        print(f"  {nm} algebraic solutions (exps t=1 {e1}, t=2 {e2}, inf {ei}):", flush=True)
        rec[f"{nm}_alg"] = EX.algebraic_solutions(Q, e1, e2, ei)
    (HERE / "data" / "p3e2w_slice_F.json").write_text(json.dumps(rec, indent=1, default=str))


if __name__ == "__main__":
    main()
