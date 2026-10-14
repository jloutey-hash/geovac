"""Three-electron chain integral: fingerprint check (pre-registered:
docs/predictions/2026-10-13_three_electron_k3.md).

Chain integral I = <rho_a(1) rho_b(2) rho_c(3) / (r12 r13)> = int rho_a V_b V_c.
Momentum form: I = (1/4 pi^4) int d3p d3q rho_b(p) rho_c(q) rho_a(-p-q) / (p^2 q^2).

Checks
  F1  angle integral of the Yukawa kernel e^{-D Delta}/Delta over cos(p,q) is
      elementary for every D (no log, no new root)                    [sympy]
  F2  reduced 2-D fibre == direct 3-D (p, q, x) fibre at D = 0        [numeric]
  F3  six branch lines: no four concurrent; triple-point locus         [sympy]
  F4  p-fixed slice is genus 1 with p-dependent j-invariant            [sympy]
      generic line slice has 6 distinct branch points (genus 2)       [numeric]
  F5  one-centre degeneration: momentum route (reduced 2-D) ==
      exact position-space closed form, which contains only logs      [sympy+numeric]
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.integrate import dblquad, tplquad

OUT: dict = {}


# ---------------------------------------------------------------- F1
def f1_angle_integral() -> None:
    c, D, p, q, x = sp.symbols("c D p q x", positive=True)
    u = p**2 + q**2 + 2 * p * q * x
    Dl = sp.Symbol("Delta", positive=True)
    # du = 2 p q dx ;  u = (Delta^2 - 1)/c  =>  du = 2 Delta dDelta / c
    # int dx e^{-D Delta}/Delta = (1/(2pq)) int (2 Delta/c) e^{-D Delta}/Delta dDelta
    anti = sp.integrate(sp.exp(-D * Dl) / Dl * 2 * Dl / c, Dl) / (2 * p * q)
    anti = sp.simplify(anti)
    # verify by differentiating back in x
    expr = anti.subs(Dl, sp.sqrt(c * u + 1))
    target = sp.exp(-D * sp.sqrt(c * u + 1)) / sp.sqrt(c * u + 1)
    resid = sp.simplify(sp.diff(expr, x) - target)
    # D = 0 Slater fingerprint density 1/Delta^5 (the D=0 value of G, Paper 59)
    anti5 = sp.integrate(Dl**-5 * 2 * Dl / c, Dl) / (2 * p * q)
    expr5 = anti5.subs(Dl, sp.sqrt(c * u + 1))
    resid5 = sp.simplify(sp.diff(expr5, x) - (c * u + 1) ** sp.Rational(-5, 2))
    has_log = anti.has(sp.log) or anti5.has(sp.log)
    print("F1 antiderivative (Yukawa):", anti, " residual:", resid)
    print("F1 antiderivative (1/Delta^5):", sp.simplify(anti5), " residual:", resid5)
    OUT["F1"] = {"yukawa_antiderivative": str(anti), "residual": str(resid),
                 "d0_antiderivative": str(sp.simplify(anti5)),
                 "residual_d0": str(resid5), "log_present": bool(has_log)}


# ---------------------------------------------------------------- F2
def rho(k: float, c: float) -> float:
    return (c * k * k + 1.0) ** -2.5


def f2_reduced_vs_direct(ca: float, cb: float, cc: float) -> None:
    # fibre: int p^2 dp q^2 dq dx  rho_b(p) rho_c(q) rho_a(|p+q|) / (p^2 q^2)
    #      = int dp dq dx rho_b rho_c rho_a
    def direct_integrand(x, q, p):
        u = p * p + q * q + 2 * p * q * x
        return rho(p, cb) * rho(q, cc) * (ca * u + 1.0) ** -2.5

    # map p, q in [0, inf) via p = t/(1-t)
    def direct_mapped(x, tq, tp):
        p = tp / (1 - tp); q = tq / (1 - tq)
        return direct_integrand(x, q, p) / (1 - tp) ** 2 / (1 - tq) ** 2

    direct, derr = tplquad(direct_mapped, 0, 1, 0, 1, -1, 1,
                           epsabs=1e-12, epsrel=1e-10)

    def reduced_mapped(tq, tp):
        p = tp / (1 - tp); q = tq / (1 - tq)
        if p * q < 1e-14:
            # limit: int dx rho_a -> 2 rho_a(max(p,q))
            ang = 2.0 * (ca * max(p, q) ** 2 + 1.0) ** -2.5
        else:
            Dm = np.sqrt(ca * (p - q) ** 2 + 1.0)
            Dp = np.sqrt(ca * (p + q) ** 2 + 1.0)
            ang = (Dm ** -3 - Dp ** -3) / (3.0 * ca * p * q)
        return rho(p, cb) * rho(q, cc) * ang / (1 - tp) ** 2 / (1 - tq) ** 2

    reduced, rerr = dblquad(reduced_mapped, 0, 1, 0, 1, epsabs=1e-13, epsrel=1e-12)
    print(f"F2 c=({ca},{cb},{cc})  direct={direct:.12f}  reduced={reduced:.12f}"
          f"  diff={abs(direct - reduced):.2e}")
    OUT.setdefault("F2", []).append({"c": [ca, cb, cc], "direct": direct,
                                     "reduced": reduced,
                                     "abs_diff": abs(direct - reduced)})


# ---------------------------------------------------------------- F3
def f3_line_configuration() -> None:
    a, b, c = sp.symbols("sa sb sc", positive=True)  # square roots of c_a, c_b, c_c
    I = sp.I
    lines = {  # coefficient vectors (p, q, w) of the six lines
        "b+": (b, 0, I), "b-": (b, 0, -I),
        "c+": (0, c, I), "c-": (0, c, -I),
        "a+": (a, a, I), "a-": (a, a, -I),     # the (p+q) term
    }
    triples = {}
    for trio in itertools.combinations(lines, 3):
        det = sp.factor(sp.Matrix([lines[k] for k in trio]).det())
        triples["".join(trio)] = str(det)
    # four concurrent lines need every sub-triple determinant to vanish
    quad_possible = []
    for quad in itertools.combinations(lines, 4):
        dets = [sp.Matrix([lines[k] for k in t]).det()
                for t in itertools.combinations(quad, 3)]
        sol = sp.solve([sp.expand(d) for d in dets], [a], dict=True)
        sol = [s for s in sol if s[a].is_positive is not False and s[a] != 0]
        if sol:
            quad_possible.append(("".join(quad), str(sol)))
    print("F3 triple determinants (nonzero unless on a locus):")
    for k, v in triples.items():
        print("   ", k, v)
    print("F3 quadruple-concurrency solutions with positive c:", quad_possible)
    OUT["F3"] = {"triple_dets": triples, "quadruple_concurrent": quad_possible}


# ---------------------------------------------------------------- F4
def quartic_j(coeffs):
    a0, a1, a2, a3, a4 = coeffs
    I = 12 * a0 * a4 - 3 * a1 * a3 + a2**2
    J = 72 * a0 * a2 * a4 + 9 * a1 * a2 * a3 - 27 * a0 * a3**2 - 27 * a4 * a1**2 - 2 * a2**3
    return sp.simplify(1728 * 4 * I**3 / (4 * I**3 - J**2))


def f4_slices(ca: float, cb: float, cc: float) -> None:
    p, q = sp.symbols("p q")
    ca_, cc_ = sp.Rational(ca).limit_denominator(1000), sp.Rational(cc).limit_denominator(1000)
    quart = sp.Poly(sp.expand((cc_ * q**2 + 1) * (ca_ * (p + q) ** 2 + 1)), q)
    co = [quart.coeff_monomial(q**k) for k in (4, 3, 2, 1, 0)]
    j = quartic_j(co)
    dj = sp.simplify(sp.diff(j, p))
    jvals = {str(pv): float(j.subs(p, pv)) for pv in (0, sp.Rational(1, 2), 1, 2, 5)}
    print("F4 j(p) of the p-fixed slice at sample p:", jvals)
    print("F4 dj/dp identically zero?", dj == 0)
    # generic line slice q = lam p + mu: sextic in p, count distinct roots
    rng = np.random.default_rng(1)
    distinct = []
    for _ in range(5):
        lam, mu = rng.normal(), rng.normal()
        poly = np.polymul(np.polymul([cb, 0, 1], [cc * lam**2, 2 * cc * lam * mu, cc * mu**2 + 1]),
                          [ca * (1 + lam) ** 2, 2 * ca * (1 + lam) * mu, ca * mu**2 + 1])
        r = np.roots(poly)
        dmin = min(abs(r[i] - r[k]) for i in range(6) for k in range(i + 1, 6))
        distinct.append({"lam": lam, "mu": mu, "n_roots": len(r), "min_root_sep": float(dmin)})
    print("F4 generic line slices (6 distinct roots => genus 2):",
          [(d["n_roots"], round(d["min_root_sep"], 4)) for d in distinct])
    OUT["F4"] = {"j_samples": jvals, "dj_dp_zero": bool(dj == 0),
                 "line_slices": distinct}


# ---------------------------------------------------------------- F5
def f5_one_centre(za: int, zb: int, zc: int) -> None:
    r = sp.symbols("r", positive=True)

    def dens(z):  # normalized 1s^2 density, total charge 1
        return z**3 / sp.pi * sp.exp(-2 * z * r)

    def pot(z):
        return 1 / r - sp.exp(-2 * z * r) * (z + 1 / r)

    exact = sp.simplify(sp.integrate(4 * sp.pi * r**2 * dens(za) * pot(zb) * pot(zc),
                                     (r, 0, sp.oo)))
    exact_f = float(exact)

    def rt(k, z):  # Fourier transform of the density
        return 16.0 * z**4 / (k * k + 4.0 * z * z) ** 2

    def red(tq, tp):
        p = tp / (1 - tp); q = tq / (1 - tq)
        if p * q < 1e-14:
            ang = 2.0 * rt(max(p, q), za)
        else:
            ang = 8.0 * za**4 / (p * q) * (1.0 / ((p - q) ** 2 + 4 * za * za)
                                           - 1.0 / ((p + q) ** 2 + 4 * za * za))
        return rt(p, zb) * rt(q, zc) * ang / (1 - tp) ** 2 / (1 - tq) ** 2

    val, _ = dblquad(red, 0, 1, 0, 1, epsabs=1e-13, epsrel=1e-12)
    mom = 2.0 / np.pi**2 * val
    atoms = sorted({type(f).__name__ for f in sp.preorder_traversal(exact)
                    if isinstance(f, sp.Function)})
    print(f"F5 zeta=({za},{zb},{zc}) exact = {exact}  = {exact_f:.12f}")
    print(f"F5 momentum route = {mom:.12f}  diff = {abs(mom - exact_f):.2e};"
          f" functions in closed form: {atoms}")
    OUT.setdefault("F5", []).append({"zeta": [za, zb, zc], "exact": str(exact),
                                     "exact_float": exact_f, "momentum": mom,
                                     "abs_diff": abs(mom - exact_f),
                                     "functions": atoms})


def main() -> None:
    f1_angle_integral()
    for cs in [(0.21, 0.09, 0.16), (0.05, 0.24, 0.13), (0.25, 0.25, 0.01)]:
        f2_reduced_vs_direct(*cs)
    f3_line_configuration()
    f4_slices(0.21, 0.09, 0.16)
    for zs in [(1, 1, 1), (1, 2, 3), (2, 1, 3)]:
        f5_one_centre(*zs)
    Path("debug/data").mkdir(parents=True, exist_ok=True)
    Path("debug/data/three_electron_k3_check.json").write_text(json.dumps(OUT, indent=2, default=str))


if __name__ == "__main__":
    main()
