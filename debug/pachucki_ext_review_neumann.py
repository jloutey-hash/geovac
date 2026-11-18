"""Reviewer's independent two-centre integrals (own code; no geovac import, no Pachucki closed form).

Prolate spheroidal coordinates xi = (rA+rB)/R, eta = (rA-rB)/R, Neumann expansion of 1/r12 (m = 0 part only, the
densities are axially symmetric):
    1/r12 -> (2/R) sum_l (2l+1) P_l(xi<) Q_l(xi>) P_l(eta1) P_l(eta2),   d^3r = (R/2)^3 (xi^2-eta^2) dxi deta dphi.
For a density rho = W(xi, eta) e^{-al xi - be eta}, W polynomial (al = R(a+b)/2, be = R(a-b)/2):
    I[W] = (2/R)(R/2)^6 (2pi)^2 sum_l (2l+1) * 2 * int_1^oo Q_l g_l e^{-al x} [ int_1^x P_l g_l e^{-al y} dy ] dx,
    g_l(xi) = int_{-1}^{1} W(xi, eta) P_l(eta) e^{-be eta} deta      (W already includes the (xi^2-eta^2) Jacobian).
The inner integral is done in closed form (polynomial x exponential), the eta moments and the outer integral by
mpmath quadrature.  Cases: K for 1s, 2psig, 3dsig, ppsig, pppi0 (normalised STO products, A at eta = -1 side:
rA = R(xi+eta)/2, rB = R(xi-eta)/2, z_A = (rA^2 + R^2 - rB^2)/(2R) = R(1 + xi eta)/2, z_B = z_A - R), and Pachucki's
master integral f0 = R * int int e^{-a r1A - b r1B - a r2A - b r2B}/(r1A r1B r2A r2B r12) d1 d2/(4pi)^2.
Usage: python debug/pachucki_ext_review_neumann.py CASE a b R [dps] [lmax]
"""
import sys
from fractions import Fraction

import mpmath as mp
import sympy as sp


def neumann(case, a, b, R, dps=40, lmax=26, verbose=False):
    mp.mp.dps = dps + 25
    a, b, R = sp.Rational(str(a)), sp.Rational(str(b)), sp.Rational(str(R))
    xi, eta = sp.symbols("xi eta")
    rA, rB = R * (xi + eta) / 2, R * (xi - eta) / 2
    zA = R * (1 + xi * eta) / 2          # A at origin, B at z = R:  z = (rA^2 + R^2 - rB^2)/(2R)
    zB = zA - R
    jac = xi ** 2 - eta ** 2
    pi2 = sp.pi ** 2
    if case == "1s":
        W, N2 = jac, a ** 3 * b ** 3 / pi2
    elif case == "2psig":
        W, N2 = jac * zA, a ** 5 * b ** 3 / pi2
    elif case == "3dsig":
        W, N2 = jac * (3 * zA ** 2 - rA ** 2), a ** 7 * b ** 3 / (18 * pi2)
    elif case == "ppsig":
        W, N2 = jac * zA * zB, a ** 5 * b ** 5 / pi2
    elif case == "pppi0":
        # conj(2p+1_A) 2p+1_B = (1/2) N_A N_B (x^2 + y^2) e^{..},  x^2+y^2 = rA^2 - zA^2
        W, N2 = jac * (rA ** 2 - zA ** 2), a ** 5 * b ** 5 / (4 * pi2)
    elif case == "f0":
        # 1/(rA rB) = (2/R)^2 / jac  cancels the Jacobian; prefactor handled below
        W, N2 = sp.Integer(1), None
    else:
        raise ValueError(case)
    Wp = sp.Poly(sp.expand(W), xi, eta)
    al = mp.mpf(R * (a + b) / 2)
    be = mp.mpf(R * (a - b) / 2)
    Rm = mp.mpf(R)
    kmax = max(m[1] for m in Wp.monoms())
    total = mp.mpf(0)
    per = []
    for l in range(lmax + 1):
        Pl = [mp.mpf(c) for c in reversed(sp.Poly(sp.legendre(l, xi), xi).all_coeffs())]   # ascending
        mom = [mp.quad(lambda e, k=k: e ** k * mp.legendre(l, e) * mp.exp(-be * e), mp.linspace(-1, 1, 5))
               for k in range(kmax + 1)]
        g = {}
        for (j, k), c in zip(Wp.monoms(), Wp.coeffs()):
            g[j] = g.get(j, mp.mpf(0)) + mp.mpf(c) * mom[k]
        gl = [g.get(j, mp.mpf(0)) for j in range(max(g) + 1)]
        # q = P_l * g_l (ascending coefficients)
        q = [mp.mpf(0)] * (len(Pl) + len(gl) - 1)
        for i, x in enumerate(Pl):
            for j, y in enumerate(gl):
                q[i + j] += x * y
        # h = sum_j q^{(j)} / al^{j+1}:  d/dx[-h e^{-al x}] = q e^{-al x}
        h = [mp.mpf(0)] * len(q)
        d = list(q)
        j = 0
        while any(d):
            for i, x in enumerate(d):
                h[i] += x / al ** (j + 1)
            d = [i * d[i] for i in range(1, len(d))]
            j += 1
        h1 = mp.polyval(h[::-1], 1) * mp.exp(-al)

        def outer(x):
            A = h1 - mp.polyval(h[::-1], x) * mp.exp(-al * x)
            return mp.legenq(l, 0, x, type=3).real * mp.polyval(gl[::-1], x) * mp.exp(-al * x) * A

        val = mp.quad(outer, [1, mp.mpf("1.25"), 2, 4, 8, 16, 40, mp.inf])
        term = (2 * l + 1) * 2 * val
        per.append(term)
        total += term
        if verbose:
            print(l, mp.nstr(term, 8), flush=True)
        if l > kmax + 2 and abs(term) < mp.mpf(10) ** (-(dps + 8)) * abs(total):
            break
    if case == "f0":
        # f0 = R/(4pi)^2 * (2/R)(2pi)^2 (R/2)^6 (2/R)^4 * S = R^2/8 * S ... (R/2)^6 (2/R)^4 = R^2/4; times (2/R) R /4
        pref = Rm * (2 / Rm) * (Rm / 2) ** 6 * (2 / Rm) ** 4 / 4
    else:
        pref = mp.mpf(sp.N(N2 * pi2, dps + 30)) / mp.pi ** 2 * (2 / Rm) * (Rm / 2) ** 6 * (2 * mp.pi) ** 2
    return pref * total, per


if __name__ == "__main__":
    case, a, b, R = sys.argv[1:5]
    dps = int(sys.argv[5]) if len(sys.argv) > 5 else 40
    lmax = int(sys.argv[6]) if len(sys.argv) > 6 else 26
    v, per = neumann(case, Fraction(a), Fraction(b), Fraction(R), dps, lmax)
    print(case, a, b, R, mp.nstr(v, dps), "terms", len(per), "last", mp.nstr(per[-1], 3))
