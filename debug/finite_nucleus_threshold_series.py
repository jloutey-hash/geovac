"""Derive the s-wave defect of a uniformly charged nucleus as a series in x = ZR (owed item from
docs/predictions/2026-10-18_finite_nucleus.md: the measured -80/63 and the 2/5 on x^2 ln x).

Route: Seaton's limit (confirmed to 1e-15) lets the defect be computed at E = 0. Inside, with s = r/R,
u'' = -x (3 - s^2) u, so u is a power series in x. Outside, the energy-normalized zero-energy pair is
f = sqrt(2r) J_1(sqrt(8r)), g = -sqrt(2r) Y_1(sqrt(8r)) (Wronskian 2/pi, checked below), Z = 1.
Matching u = f cos(Delta) + g sin(Delta) at r = R = x gives tan(Delta) = (x f' - lam f)/(lam g - x g'),
lam = u_s(1)/u(1), and the defect is mu = Delta/pi.
"""
from __future__ import annotations

import sympy as sp

x, s, Lg = sp.symbols("x s Lg")
ORDER = 5  # keep x^5 in mu


def interior_lambda(order: int):
    us = [s]
    for j in range(1, order + 1):
        rhs = -(3 - s**2) * us[j - 1]
        uj = sp.integrate(sp.integrate(rhs, (s, 0, s)), (s, 0, s))
        us.append(sp.expand(uj))
    u = sum(x**j * us[j] for j in range(order + 1))
    lam = sp.diff(u, s).subs(s, 1) / u.subs(s, 1)
    return sp.series(lam, x, 0, order + 1).removeO()


def exterior(K: int = 8):
    """f, g as series in x with log(2x) kept as a symbol Lg; x d/dx applied before substituting."""
    gam = sp.EulerGamma
    f = sum((-1) ** k * (2 * x) ** (k + 1) / (sp.factorial(k) * sp.factorial(k + 1)) for k in range(K))
    # Y_1(z) = (2/pi) ln(z/2) J_1(z) - 2/(pi z) - (1/pi) sum_k (psi(k+1) + psi(k+2)) (-1)^k (z/2)^(2k+1)/(k!(k+1)!)
    # with z/2 = sqrt(2x), ln(z/2) = log(2x)/2; g = -sqrt(2x) Y_1 written as a series in x
    psi = lambda m: sp.harmonic(m - 1) - gam
    sumterm = sum((psi(k + 1) + psi(k + 2)) * (-1) ** k * (2 * x) ** (k + 1) / (sp.factorial(k) * sp.factorial(k + 1))
                  for k in range(K))
    g = -(sp.log(2 * x) / sp.pi * f - 1 / sp.pi - sumterm / sp.pi)
    xf = x * sp.diff(f, x)
    xg = x * sp.diff(g, x)
    sub = lambda e: sp.expand(e.subs(sp.log(2 * x), Lg))
    return sub(f), sub(g), sub(xf), sub(xg)


def main() -> None:
    lam = interior_lambda(ORDER)
    f, g, xf, xg = exterior()
    t = (xf - lam * f) / (lam * g - xg)
    tser = sp.series(t, x, 0, ORDER + 1).removeO()
    mu = sp.expand(sp.series(sp.atan(tser) / sp.pi, x, 0, ORDER + 1).removeO())
    lead = -sp.Rational(2, 5) * x**2
    rel = sp.expand(sp.series(mu / lead, x, 0, ORDER - 1).removeO())
    rel = sp.collect(sp.expand(rel.subs(Lg, sp.log(2) + sp.Symbol("lnx"))), x)
    print("mu =", sp.collect(mu, x))
    print("mu / (-(2/5) x^2) =", rel)
    # Wronskian check of the exterior pair at a sample point (numerically)
    import mpmath as mp
    mp.mp.dps = 30
    r0 = mp.mpf("0.37")
    fN = lambda r: mp.sqrt(2 * r) * mp.besselj(1, mp.sqrt(8 * r))
    gN = lambda r: -mp.sqrt(2 * r) * mp.bessely(1, mp.sqrt(8 * r))
    W = fN(r0) * mp.diff(gN, r0) - mp.diff(fN, r0) * gN(r0)
    print("Wronskian * pi / 2 =", mp.nstr(W * mp.pi / 2, 20))
    # series vs mpmath for g at a small x
    gS = sp.N(g.subs(Lg, sp.log(2 * sp.Rational(1, 100))).subs(x, sp.Rational(1, 100)), 25)
    print("g series vs bessely at x = 0.01:", gS, mp.nstr(gN(mp.mpf(1) / 100), 25))


if __name__ == "__main__":
    main()
