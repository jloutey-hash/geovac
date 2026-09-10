"""General-l (m=0) two-centre one-body integrals, closed form.

WHY THIS IS SMALL.  The s-only closed forms in geovac/qfd_core.py stalled the
molecular Sturmian build at an 11.7 mHa l-truncation floor.  Extending to
general (l,m) in closed form looked like a derivation sprint -- but the floor is
entirely SIGMA character (H2+'s ground state is 1 sigma_g, so only m=0
contributes), and for m=0 the solid harmonic is a POLYNOMIAL in the prolate
coordinates:

    r_C^l P_l(cos theta_C)  is a polynomial in (xi, eta), degree l in each
    (verified symbolically for l = 0..4)

and since xi = (r_A + r_B)/R, eta = (r_A - r_B)/R, it is a polynomial in
(r_A, r_B) -- exactly what qfd_core.I2c already integrates.  So the extension
reuses the tested Mulliken auxiliary machinery monomial by monomial;  nothing
new is derived, only multiplied.

Normalisation composes just as simply.  For m=0 the angular factor is
sqrt((2l+1)/4pi) P_l, so a two-centre element picks up
sqrt((2l+1)(2l'+1)) / (4pi) -- and at l = l' = 0 that is the 1/(4pi) qfd_core
already carries.  So l=0 MUST reproduce qfd_core exactly, which is the test.
"""
import os, sys
from fractions import Fraction
import sympy as sp
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from geovac.qfd_core import I2c
from geovac.two_center_eri import radial_poly, radial_norm

_XI, _ETA = sp.symbols("xi eta")


def solid_harm_rArB(l: int, center: str, R) -> dict:
    """r_C^l P_l(cos theta_C) as {(i, j): coeff} over r_A^i r_B^j."""
    if center == "A":
        expr = (R * (_XI + _ETA) / 2) ** l * sp.legendre(l, (_XI * _ETA + 1) / (_XI + _ETA))
    else:
        expr = (R * (_XI - _ETA) / 2) ** l * sp.legendre(l, (_XI * _ETA - 1) / (_XI - _ETA))
    expr = sp.expand(sp.cancel(sp.together(expr)))
    rA, rB = sp.symbols("rA rB")
    expr = sp.expand(expr.subs({_XI: (rA + rB) / R, _ETA: (rA - rB) / R}))
    p = sp.Poly(expr, rA, rB)
    return {(i, j): sp.nsimplify(c) for (i, j), c in zip(p.monoms(), p.coeffs())}


def orbital_rArB(center: str, Z, n: int, l: int, R) -> tuple:
    """Spatial polynomial {(i,j): coeff} and decay rate for one orbital.

    R_nl(r) = sum_k c_k r^k e^{-a r} has lowest power r^l, so factoring r^l out
    leaves non-negative powers;  multiply by the solid-harmonic polynomial.
    """
    Zf = Z if isinstance(Z, Fraction) else Fraction(Z)
    c, a = radial_poly(Zf, n, l)
    N = radial_norm(Zf, n, l)
    sh = solid_harm_rArB(l, center, R)
    out = {}
    for k, v in c.items():
        assert k >= l, "radial power below r^l"
        for (i, j), w in sh.items():
            key = (i + (k - l), j) if center == "A" else (i, j + (k - l))
            out[key] = out.get(key, 0) + N * v * w
    return out, sp.Rational(Zf.numerator, Zf.denominator * n)


def inv_r_sigma(oi, oj, which: str, R):
    """<i| 1/r_C |j>, both orbitals m=0, on OPPOSITE centres.  oi=(center,Z,n,l)."""
    assert oi[0] != oj[0], "this routine is the genuinely two-centre case"
    pA, aA = orbital_rArB(*(oi if oi[0] == "A" else oj), R=R)
    pB, aB = orbital_rArB(*(oj if oi[0] == "A" else oi), R=R)
    ang = sp.sqrt((2 * oi[3] + 1) * (2 * oj[3] + 1))
    shift = (-1, 0) if which == "A" else (0, -1)
    tot = sp.Integer(0)
    for (i1, j1), v1 in pA.items():
        for (i2, j2), v2 in pB.items():
            i, j = i1 + i2 + shift[0], j1 + j2 + shift[1]
            if i < -1 or j < -1:
                continue
            tot += v1 * v2 * I2c(i, j, aA, aB, R)
    return sp.expand(ang * tot / (4 * sp.pi))


if __name__ == "__main__":
    from geovac import qfd_core as Q
    R = sp.Integer(2)
    print("VALIDATION: l=0 must reproduce qfd_core._inv_r EXACTLY")
    ok = True
    for (nA, nB) in ((1, 1), (1, 2), (2, 3)):
        for which in ("A", "B"):
            mine = inv_r_sigma(("A", Fraction(1), nA, 0), ("B", Fraction(1), nB, 0), which, R)
            theirs = Q._inv_r(("A", Fraction(1), nA), ("B", Fraction(1), nB), which, R)
            d = sp.simplify(mine - theirs)
            ok &= (d == 0)
            print("   n=(%d,%d) C=%s : difference = %s" % (nA, nB, which, d))
    print("   l=0 reproduction EXACT:", ok)
