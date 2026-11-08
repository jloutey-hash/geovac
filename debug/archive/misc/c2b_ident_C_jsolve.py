"""EC3 stage C (docs/predictions/2026-10-22_c2b_identification.md): every rational J of degree 6 with
poles only at x = 0, 1/4, 25/4, infinity.

J = P^3 / (lam D),  J - 1728 = Q^2 / (lam D),  i.e.  P^3 - Q^2 = 1728 lam D,
D = x^n (4x - 1)^a (4x - 25)^b,  n + a + b + c = 6, all widths >= 1 (c = pole order at infinity).
P monic of degree 2 fixes the scaling; then Q is monic (Q -> -Q gives the same J).

Elimination, exact: the x^5, x^4, x^3 coefficients are linear in q2, q1, q0; the x^2, x^1, x^0
coefficients give e2, e1, e0 in (p0, p1, lam); r12 = Res_lam(e2, e1), r10 = Res_lam(e2, e0); their
common factor (4 p0 - p1^2)^2 is the degenerate square P (lam = 0) and is removed; rr = Res_p1(r12, r10)
is univariate in p0. Roots: numerical at 80 digits from the exact polynomials, so the list is complete;
each solution is checked on e2, e1, e0, identified, and re-checked exactly when rational.

Usage: python debug/c2b_ident_C_jsolve.py [planted]"""
import itertools
import json
import sys
from pathlib import Path

import mpmath as mp
import sympy as sp

x = sp.Symbol("x")
p0, p1, q0, q1, q2, lam = sp.symbols("p0 p1 q0 q1 q2 lam")
OUT = Path(__file__).parent / "data"
mp.mp.dps = 80
TINY = mp.mpf(10) ** -30


def D_of(pts, widths):
    D = sp.Integer(1)
    for s, w in zip(pts, widths):
        D *= x**w if s == 0 else (sp.denom(s) * x - sp.numer(s)) ** w
    return sp.expand(D)


def roots_mp(coeffs):
    cs = list(coeffs)
    while cs and abs(cs[0]) == 0:
        cs.pop(0)
    if len(cs) < 2:
        return []
    return mp.polyroots(cs, maxsteps=4000, extraprec=4000)


def coeffs_mp(expr, var, env):
    out = []
    syms = list(env.keys())
    for c in sp.Poly(expr, var).all_coeffs():
        out.append(mp.mpc(sp.lambdify(syms, c, "mpmath")(*env.values())))
    return out


def safe_identify(v):
    if abs(v) < mp.mpf(10) ** -40:
        return "0"
    try:
        return mp.identify(v)
    except ValueError:
        return None


def solve_D(D):
    P = x**2 + p1 * x + p0
    Q = x**3 + q2 * x**2 + q1 * x + q0
    E = sp.Poly(sp.expand(P**3 - Q**2 - 1728 * lam * D), x)
    co = {k: E.coeff_monomial(x**k) for k in range(7)}
    assert sp.expand(co[6]) == 0
    s2 = sp.solve(co[5], q2)[0]
    s1 = sp.solve(co[4].subs(q2, s2), q1)[0]
    s0 = sp.solve(co[3].subs({q2: s2, q1: s1}), q0)[0]
    sub = {q2: s2, q1: s1, q0: s0}
    e = [sp.expand(co[k].subs(sub)) for k in (2, 1, 0)]
    r12 = sp.Poly(sp.resultant(e[0], e[1], lam), p0, p1)
    r10 = sp.Poly(sp.resultant(e[0], e[2], lam), p0, p1)
    g = sp.gcd(r12, r10)
    note = None
    if g.total_degree() > 0:
        fg = sp.factor(g.as_expr())
        assert sp.simplify(fg / (4 * p0 - p1**2) ** 2).is_number, fg
        note = str(fg)
        r12, r10 = sp.quo(r12, g), sp.quo(r10, g)
    rr = sp.resultant(r12.as_expr(), r10.as_expr(), p1)
    assert rr != 0
    univ = sp.factor_list(rr)[1]
    fe = [sp.lambdify((p0, p1, lam), ei, "mpmath") for ei in e]
    fs = {str(k): sp.lambdify((p0, p1, lam), v, "mpmath") for k, v in sub.items()}
    sols, keys = [], []
    for f, _ in univ:
        if sp.degree(f, p0) < 1:
            continue
        for v0 in roots_mp([mp.mpf(int(c)) for c in sp.Poly(f, p0).all_coeffs()]):
            for v1 in roots_mp(coeffs_mp(r12.as_expr(), p1, {p0: v0})):
                for vl in roots_mp(coeffs_mp(e[0], lam, {p0: v0, p1: v1})):
                    if abs(vl) < TINY:
                        continue
                    if max(abs(fi(v0, v1, vl)) for fi in fe) > mp.mpf(10) ** -40 * (1 + abs(vl)) ** 3:
                        continue
                    if any(abs(v0 - a) < TINY and abs(v1 - b) < TINY and abs(vl - c) < TINY for a, b, c in keys):
                        continue
                    keys.append((v0, v1, vl))
                    vals = {"p0": v0, "p1": v1, "lam": vl}
                    vals.update({k: f2(v0, v1, vl) for k, f2 in fs.items()})
                    real = all(abs(v.imag) < mp.mpf(10) ** -40 for v in vals.values())
                    ident = {k: (safe_identify(v.real) if real else None) for k, v in vals.items()}
                    sols.append({"values": {k: mp.nstr(v, 45) for k, v in vals.items()},
                                 "identified": ident, "real": real, "factor_degree": int(sp.degree(f, p0))})
    return sols, note, [int(sp.degree(f, p0)) for f, _ in univ]


def exact_check(D, ident):
    try:
        v = {k: sp.Rational(sp.nsimplify(ident[k], rational=True)) for k in ("p0", "p1", "lam", "q2", "q1", "q0")}
    except (TypeError, ValueError, KeyError, AttributeError, sp.SympifyError):
        return None
    P = x**2 + v["p1"] * x + v["p0"]
    Q = x**3 + v["q2"] * x**2 + v["q1"] * x + v["q0"]
    if sp.expand(P**3 - Q**2 - 1728 * v["lam"] * D) != 0:
        return None
    if sp.degree(sp.gcd(P, Q), x) > 0:
        return None
    return {"P": str(P), "Q": str(Q), "lam": str(v["lam"]), "J": str(sp.factor(P**3 / (v["lam"] * D)))}


def main(planted=False):
    if planted:
        Pt = x**2 + 3 * x - 1
        Qt = sp.expand(sp.series(sp.sqrt(1 + 3 / x - 1 / x**2) ** 3, x, sp.oo, 4).removeO() * x**3)
        R = sp.Poly(sp.expand(Pt**3 - Qt**2), x)
        lam0 = sp.LC(R) / 1728
        D = sp.expand(R.as_expr() / (1728 * lam0))
        sols, note, degs = solve_D(D)
        ex = [exact_check(D, s["identified"]) for s in sols if s["real"]]
        ok = any(e and sp.expand(sp.sympify(e["P"]) - Pt) == 0 for e in ex)
        print("planted P =", Pt, " D =", sp.factor(D), " recovered exactly:", ok, " solutions:", len(sols),
              " removed:", note)
        return
    pts = [sp.Integer(0), sp.Rational(1, 4), sp.Rational(25, 4)]
    allres = []
    for n, a, b in itertools.product(range(1, 5), repeat=3):
        c = 6 - n - a - b
        if c < 1:
            continue
        D = D_of(pts, (n, a, b))
        sols, note, degs = solve_D(D)
        exact = [exact_check(D, s["identified"]) if s["real"] else None for s in sols]
        print(f"widths (0, 1/4, 25/4, inf) = ({n},{a},{b},{c}): {len(sols)} solution(s); p0-factor degrees {degs}",
              flush=True)
        for s, e in zip(sols, exact):
            print("   real:", s["real"], "| p0 =", s["values"]["p0"][:28], "| lam =", s["values"]["lam"][:28],
                  "| exact J:", (e["J"][:160] if e else None), flush=True)
        allres.append({"widths": [n, a, b, c], "common_factor_removed": note, "p0_factor_degrees": degs,
                       "solutions": sols, "exact": exact})
    json.dump(allres, open(OUT / "c2b_ident_C_jsolve.json", "w"), indent=1, default=str)


if __name__ == "__main__":
    main(planted=len(sys.argv) > 1 and sys.argv[1] == "planted")
