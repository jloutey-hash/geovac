"""Increment 3a: the (xi, eta) expansion of a general two-center orbital product.

FOUNDATION for the exchange class. Phase 0-e argued -- on paper -- that the whole
exchange integrand is a POLYNOMIAL in (xi, eta) times separable exponentials, and
that this is what makes the class tractable at all. It only ever machine-checked
that for 1s orbitals, where the claim is trivial. Everything else in increment 3
sits on top of it, so it gets verified for general (n, l, m) first.

THE CLAIM, precisely. With A at the origin, B at R zhat,

    r_A = R(xi+eta)/2      z_A  = R(1+xi eta)/2
    r_B = R(xi-eta)/2      rho^2 = (R^2/4)(xi^2-1)(1-eta^2)

then for chi_a^A = R_{n_a l_a}(r_A) Y_{l_a m_a}(Om_A),

    conj(chi_a^A) chi_b^B
        = P(xi, eta) * [(xi^2-1)(1-eta^2)]^{(|m_a|+|m_b|)/2}
          * exp(-p xi - q eta) * e^{i(m_b - m_a) phi}

with P a POLYNOMIAL, p = (alpha+beta)R/2 and q = (alpha-beta)R/2.

Two supporting facts, both used and both checked here:

  * R_nl(r)/r^l is a polynomial (R_nl starts at r^l), and r^l Y_lm(Om) is a solid
    harmonic -- homogeneous polynomial in (x, y, z). So the only non-polynomial
    piece is rho^{|m|}, which is what the half-power above collects.
  * |m_a| + |m_b| + |sigma| is ALWAYS even, with sigma = m_a - m_b forced by the
    azimuthal integral. So once the kernel's own (1-eta^2)^{|sigma|/2} is folded
    in, the half-power becomes an integer and nothing fractional survives.

    (same sign: |m_a|+|m_b|+||m_a|-|m_b|| = 2 max;
     opposite:  |m_a|+|m_b|+(|m_a|+|m_b|) = 2(|m_a|+|m_b|).)

Run from repo root:  python debug/inc3a_spheroidal_expansion.py
"""

from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path

import numpy as np
import sympy as sp

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from geovac.two_center_eri import (  # noqa: E402
    _legendre_rational_part, _sph_norm, radial_norm, radial_poly,
)

xi, eta, Rsym = sp.symbols("xi eta R", positive=True)
Z1, Z2, Z3 = Fraction(1), Fraction(2), Fraction(3)


def solid_harmonic_poly(l: int, m: int, z, rho2):
    """r^l P_l^m(cos th), with the rho^{|m|} factor STRIPPED, as a polynomial.

    P_l^m(w) = part(w) (1-w^2)^{|m|/2}, so with w = z/r,

        r^l P_l^m = [r^{l-|m|} part(z/r)] * rho^{|m|}

    and r^{l-|m|} part(z/r) is a polynomial in (z, r^2) because `part` has degree
    l-|m| and parity (-1)^{l-|m|} -- only even powers of r survive, and
    r^2 = rho^2 + z^2.
    """
    w = sp.Symbol("w")
    part = _legendre_rational_part(l, m, w)
    mm = abs(m)
    out = sp.Integer(0)
    poly = sp.Poly(sp.expand(part), w)
    for (j,), c in poly.terms():
        e = l - mm - j                       # power of r left over
        assert e >= 0 and e % 2 == 0, f"odd leftover power r^{e} at l={l},m={m}"
        out += c * z ** j * (rho2 + z ** 2) ** (e // 2)
    return sp.expand(out)


def two_center_product(ZA, oa, ZB, ob, R=Rsym):
    """conj(chi_a^A) chi_b^B -> (P, half_power, p, q).

    The product equals  P(xi,eta) * [(xi^2-1)(1-eta^2)]^{half_power}
                        * exp(-p xi - q eta) * e^{i(m_b-m_a) phi}.
    """
    (na, la, ma), (nb, lb, mb) = oa, ob
    ca, alpha = radial_poly(ZA, na, la)
    cb, beta = radial_poly(ZB, nb, lb)
    Na, Nb = radial_norm(ZA, na, la), radial_norm(ZB, nb, lb)

    r_A = R * (xi + eta) / 2
    r_B = R * (xi - eta) / 2
    z_A = R * (1 + xi * eta) / 2
    z_B = z_A - R
    rho2 = R ** 2 / 4 * (xi ** 2 - 1) * (1 - eta ** 2)

    # R_nl(r)/r^l is a polynomial in r; r^2 is polynomial in (xi,eta), and the
    # radial polynomials have every power >= l, so no odd r survives either.
    def radial_over_rl(coeffs, N, l, r, r2):
        out = sp.Integer(0)
        for k, c in coeffs.items():
            e = k - l
            assert e >= 0, f"radial power {k} below l={l}"
            out += N * c * (r ** e if e % 2 else r2 ** (e // 2))
        return sp.expand(out)

    radA = radial_over_rl(ca, Na, la, r_A, sp.expand(r_A ** 2))
    radB = radial_over_rl(cb, Nb, lb, r_B, sp.expand(r_B ** 2))

    # conj(Y_l m) = (-1)^m Y_{l,-m}: same modulus structure, so the polynomial
    # part is the m -> m one and the phase rides on e^{i(m_b-m_a)phi}.
    angA = _sph_norm(la, ma) * solid_harmonic_poly(la, ma, z_A, rho2)
    angB = _sph_norm(lb, mb) * solid_harmonic_poly(lb, mb, z_B, rho2)

    half = sp.Rational(abs(ma) + abs(mb), 2)
    scale = (R / 2) ** (abs(ma) + abs(mb))      # rho^{|m|} = (R/2)^{|m|} [...]^{|m|/2}
    P = sp.expand(radA * radB * angA * angB * scale)
    p = (alpha + beta) * R / 2
    q = (alpha - beta) * R / 2
    return P, half, p, q


# ------------------------------------------------------------------ validation

def _chi_numeric(Z, n, l, m, r, th, ph):
    c, a = radial_poly(Z, n, l)
    N = radial_norm(Z, n, l)
    rad = sum(float(N * cc) * r ** k for k, cc in c.items()) * np.exp(-float(a) * r)
    return rad * complex(sp.Ynm(l, m, th, ph).expand(func=True).evalf())


def main() -> None:
    print("Increment 3a -- (xi, eta) expansion of a two-center orbital product\n")
    Rv = 2.5

    print("V1  is the constructed P actually a POLYNOMIAL in (xi, eta)?")
    cases = [((1, 0, 0), (1, 0, 0)), ((2, 1, 0), (1, 0, 0)),
             ((2, 1, 1), (2, 1, 1)), ((2, 1, 1), (2, 1, -1)),
             ((3, 2, 2), (2, 1, 1)), ((3, 2, 0), (3, 2, 0)),
             ((3, 1, -1), (3, 2, 1))]
    for oa, ob in cases:
        P, half, _p, _q = two_center_product(Z3, oa, Z1, ob, sp.Rational(str(Rv)))
        is_poly = sp.Poly(P, xi, eta) is not None and P.is_polynomial(xi, eta)
        sigma = oa[2] - ob[2]
        tot = abs(oa[2]) + abs(ob[2]) + abs(sigma)
        print(f"    {oa} x {ob}: polynomial={is_poly}  deg={sp.Poly(P,xi,eta).total_degree()}"
              f"  half={half}  |m_a|+|m_b|+|sigma|={tot} "
              f"{'even OK' if tot % 2 == 0 else 'ODD -- CLAIM FAILS'}")

    print("\nV2  pointwise: does the expansion reproduce conj(chi_a) chi_b?")
    worst = 0.0
    for oa, ob in cases:
        P, half, p, q = two_center_product(Z3, oa, Z1, ob, sp.Rational(str(Rv)))
        Pf = sp.lambdify((xi, eta), P, "numpy")
        pf, qf = float(p), float(q)
        for xv, ev, phv in ((1.7, 0.3, 0.7), (2.4, -0.6, 2.1), (1.15, 0.85, 4.4)):
            rA = Rv * (xv + ev) / 2
            rB = Rv * (xv - ev) / 2
            zA = Rv * (1 + xv * ev) / 2
            rho = Rv / 2 * np.sqrt((xv ** 2 - 1) * (1 - ev ** 2))
            thA = np.arccos(np.clip(zA / rA, -1, 1))
            thB = np.arccos(np.clip((zA - Rv) / rB, -1, 1))
            direct = np.conj(_chi_numeric(Z3, *oa, rA, thA, phv)) \
                * _chi_numeric(Z1, *ob, rB, thB, phv)
            built = (Pf(xv, ev)
                     * ((xv ** 2 - 1) * (1 - ev ** 2)) ** float(half)
                     * np.exp(-pf * xv - qf * ev)
                     * np.exp(1j * (ob[2] - oa[2]) * phv))
            worst = max(worst, abs(direct - built))
        print(f"    {oa} x {ob}: worst so far {worst:.2e}")
    print(f"\n    worst |direct - built| = {worst:.2e}  "
          f"{'OK' if worst < 1e-10 else 'FAIL'}")


if __name__ == "__main__":
    main()
