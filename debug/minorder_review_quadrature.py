"""Reviewer's from-scratch quadrature of the hybrid (2p0 2p0 | 1s 1s_B), Z_A = 3, Z_B = 6 (and controls).

No geovac integral code is used for the reference value: the potential of the one-centre density 2p0^2 is written
from the Laplace expansion (L = 0, 2) with incomplete gamma functions, and integrated against 1s_A 1s_B in prolate
spheroidal coordinates (mpmath Gauss-Legendre, tanh-sinh fallback).  The closed form is then read from
geovac.two_center_eri.hybrid_closed_form and compared.

Orbitals (hydrogenic, normalised): 2p0 = Z^{5/2}/(4 sqrt(2 pi)) r cos(theta) e^{-Z r/2};  1s = sqrt(Z^3/pi) e^{-Z r}.
Writes debug/data/minorder_review_quadrature.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction as F
from pathlib import Path

import mpmath as mp
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

mp.mp.dps = 18


def potential_2p0sq(Z, r, ct):
    """Potential at (r, cos theta) of the density [2p0_Z]^2 centred at the origin."""
    ap = mp.mpf(Z)                       # pair rate
    N2 = mp.mpf(Z) ** 5 / (32 * mp.pi)   # density = N2 r^2 e^{-ap r} cos^2 = N2 r^2 e^{-ap r} (1/3 + 2/3 P2)
    tot = mp.mpf(0)
    for L, cL in ((0, mp.mpf(1) / 3), (2, mp.mpf(2) / 3)):
        PL = mp.mpf(1) if L == 0 else (3 * ct * ct - 1) / 2
        inner = mp.gammainc(L + 5, 0, ap * r) / ap ** (L + 5)          # int_0^r s^{L+4} e^{-ap s} ds
        outer = mp.gammainc(4 - L, ap * r) / ap ** (4 - L)             # int_r^oo s^{3-L} e^{-ap s} ds
        tot += 4 * mp.pi / (2 * L + 1) * cL * PL * N2 * (inner / r ** (L + 1) + outer * r ** L)
    return tot


def hybrid(ZA, ZB, Rv):
    Rv = mp.mpf(Rv)
    nc = mp.sqrt(mp.mpf(ZA) ** 3 / mp.pi)
    nd = mp.sqrt(mp.mpf(ZB) ** 3 / mp.pi)

    def integrand(xi, eta):
        rA = Rv * (xi + eta) / 2
        rB = Rv * (xi - eta) / 2
        if rA == 0:
            return mp.mpf(0)
        ct = (1 + xi * eta) / (xi + eta)
        return (xi * xi - eta * eta) * potential_2p0sq(ZA, rA, ct) * mp.exp(-ZA * rA - ZB * rB)

    # xi = 1 + t, split the t range; the integrand decays like e^{-(ZA+ZB) R t / 2}
    scale = 2 / ((ZA + ZB) * Rv)
    cuts = [0, scale, 4 * scale, 12 * scale, 40 * scale, 120 * scale]
    val = mp.quad(lambda t, eta: integrand(1 + t, eta), cuts, [-1, -0.5, 0, 0.5, 1], maxdegree=6)
    return 2 * mp.pi * (Rv / 2) ** 3 * nc * nd * val


def main():
    from geovac.two_center_eri import R_s as R, hybrid_closed_form
    S1, P0 = (1, 0, 0), (2, 1, 0)
    out = {}
    for (ZA, ZB) in ((3, 6), (3, 4), (3, 3)):
        f = hybrid_closed_form(F(ZA), P0, P0, S1, F(ZB), S1)
        for Rv in ("1.4", "0.5"):
            q = hybrid(ZA, ZB, mp.mpf(Rv))
            cf = mp.mpf(str(sp.N(f.subs(R, sp.Rational(Rv)), 30)))
            out[f"Z=({ZA},{ZB}) R={Rv}"] = {"quadrature": mp.nstr(q, 20), "closed_form": mp.nstr(cf, 20),
                                           "rel_diff": mp.nstr(abs(q - cf) / abs(cf), 3)}
            print(f"Z=({ZA},{ZB}) R={Rv}: quad={mp.nstr(q, 20)} closed={mp.nstr(cf, 20)} "
                  f"rel={mp.nstr(abs(q - cf) / abs(cf), 3)}", flush=True)
            (ROOT / "debug" / "data" / "minorder_review_quadrature.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
