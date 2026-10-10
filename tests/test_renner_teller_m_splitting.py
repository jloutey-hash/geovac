"""Bending breaks the azimuthal label m: exact angular ladder (rung 1 of the
m-seam probe, PI-directed 2026-10-10).

Model: one centre with a single l-shell, two equal point charges at unit
directions (s, 0, +c) and (s, 0, -c), c = sqrt(1 - s^2), s = sin(a) where a is
half the bend away from linear (s = 0 is the linear molecule, axis z; the
charges bend toward +x).  The shell operator is

    V = sum_k w_k sum_{charges} P_k(n . rhat),

with the multipole weights w_k left SYMBOLIC: they stand for the radial factors
(for a real atom w_k = -Z <r_<^k / r_>^{k+1}>), which carry the transcendental
content.  Everything here is exact rational arithmetic on the angular part.

What is pinned (each test names the wrong answer it rejects):
  * linear geometry: each +-m pair exactly degenerate and uncoupled;
  * |m| = 1: eigenvalue splitting (in-plane - out-of-plane) = (6/5) w_2 s^2,
    exactly -- no other state in the p shell mixes (the xz terms of the two
    charges cancel);
  * |m| = L: the DIRECT (diagonal) splitting from multipole k is zero for
    k < 2L and proportional to s^{2L} at k = 2L (L = 1, 2, 3) -- the
    Renner-Teller ladder s^{2|m|};
  * |m| = 2: x^2-y^2 also mixes with the m = 0 orbital z^2 through an
    O(s^2) off-diagonal, so the eigenvalue splitting at s^4 is
    direct (5/6) w_4 PLUS a relay term -9 (2w_2 - 5w_4)^2 / (14 (12w_2 + 5w_4)):
    the power is forced, the coefficient is a rational function of the
    radial weights.

Prior art: Renner, Z. Phys. 92, 172 (1934); Walsh, J. Chem. Soc. 2260 (1953)
(the in-plane partner is stabilised on bending).  Context in the corpus:
Paper 32 rem:multicenter_composition (linear = reducible sigma/pi, any bend
welds one irreducible block).  Not yet backing a paper equation.
"""
from __future__ import annotations

from typing import Dict, List

import numpy as np
import pytest
import sympy as sp
from sympy import Rational, gamma

x, y, z = sp.symbols('x y z', real=True)
s, c = sp.symbols('s c', real=True)
w2, w4, w6 = sp.symbols('w2 w4 w6', nonzero=True)


def sphere_int(poly: sp.Expr) -> sp.Expr:
    """Exact integral of a polynomial in (x, y, z) over the unit sphere."""
    total = sp.Integer(0)
    for (a, b, d), coef in sp.Poly(sp.expand(poly), x, y, z).terms():
        if a % 2 or b % 2 or d % 2:
            continue
        total += coef * 2 * gamma(Rational(a + 1, 2)) * gamma(Rational(b + 1, 2)) \
            * gamma(Rational(d + 1, 2)) / gamma(Rational(a + b + d + 3, 2))
    return sp.simplify(total)


def bent_multipole(k: int) -> sp.Expr:
    """sum over the two bent charges of P_k(n . rhat)."""
    return sp.expand(sum(sp.legendre(k, s * x + sg * c * z) for sg in (1, -1)))


def eliminate_c(expr: sp.Expr) -> sp.Expr:
    """Rewrite in s alone using c^2 = 1 - s^2; odd powers of c must cancel."""
    poly = sp.Poly(sp.expand(expr), c)
    out = sp.Integer(0)
    for (p,), coef in poly.terms():
        assert p % 2 == 0, "odd power of c survived: the two charges are not mirror images"
        out += coef * (1 - s**2) ** (p // 2)
    return sp.expand(out)


def shell_matrix(basis: Dict[str, sp.Expr], weights: Dict[int, sp.Expr]) -> sp.Matrix:
    """Hamiltonian in the L2-normalised (orthogonal) real basis, function of s."""
    V = sum(w * bent_multipole(k) for k, w in weights.items())
    names = list(basis)
    norms = {n: sphere_int(basis[n] ** 2) for n in names}
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            assert sphere_int(basis[a] * basis[b]) == 0, "basis not orthogonal"
    H = sp.zeros(len(names), len(names))
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            H[i, j] = eliminate_c(sphere_int(basis[a] * basis[b] * V)
                                  / sp.sqrt(norms[a] * norms[b]))
    return H


P_SHELL = {'x': x, 'y': y, 'z': z}
D_SHELL = {'x2-y2': x**2 - y**2, 'xy': x * y, 'z2': 2 * z**2 - x**2 - y**2,
           'xz': x * z, 'yz': y * z}
# real |m| = L pair in the l = L shell: Re / Im of (x + i y)^L (in-plane first)
PAIRS = {L: (sp.expand(sp.re(sp.expand((x + sp.I * y) ** L))),
             sp.expand(sp.im(sp.expand((x + sp.I * y) ** L)))) for L in (1, 2, 3)}


def direct_splitting(L: int, k: int) -> sp.Expr:
    """Diagonal (in-plane minus out-of-plane) from multipole k alone."""
    pin, pout = PAIRS[L]
    V = bent_multipole(k)
    e_in = sphere_int(pin * pin * V) / sphere_int(pin * pin)
    e_out = sphere_int(pout * pout * V) / sphere_int(pout * pout)
    return sp.expand(eliminate_c(e_in - e_out))


def test_integrator_matches_independent_quadrature() -> None:
    """Exact monomial integrator == Gauss-Legendre x trapezoid quadrature
    (an independent route) on bent matrix elements at s = 0.37.
    Rejects: a wrong closed-form sphere integral."""
    sv = 0.37
    cv = np.sqrt(1 - sv**2)
    tg, wg = np.polynomial.legendre.leggauss(40)
    phi = np.linspace(0, 2 * np.pi, 64, endpoint=False)
    T, PH = np.meshgrid(tg, phi, indexing='ij')
    X, Y, Z = np.sqrt(1 - T**2) * np.cos(PH), np.sqrt(1 - T**2) * np.sin(PH), T
    W = np.outer(wg, np.full(phi.size, 2 * np.pi / phi.size))
    for f in (x**2 - y**2, x * y, x**3 - 3 * x * y**2):
        for k in (2, 4):
            integrand = f * f * bent_multipole(k)
            exact = float(sphere_int(integrand).subs({s: sv, c: cv}))
            fn = sp.lambdify((x, y, z, s, c), integrand, 'numpy')
            quad = float(np.sum(W * fn(X, Y, Z, sv, cv)))
            assert abs(exact - quad) < 1e-12 * max(1.0, abs(exact))


def test_linear_geometry_pairs_degenerate_and_uncoupled() -> None:
    """At s = 0 every +-m pair is exactly degenerate and uncoupled.
    Rejects: a symmetry leak that splits m at the linear geometry."""
    weights = {2: w2, 4: w4}
    for basis, pairs in ((P_SHELL, [('x', 'y')]),
                         (D_SHELL, [('x2-y2', 'xy'), ('xz', 'yz')])):
        H = shell_matrix(basis, weights).subs(s, 0)
        names = list(basis)
        for a, b in pairs:
            i, j = names.index(a), names.index(b)
            assert sp.simplify(H[i, i] - H[j, j]) == 0
            assert sp.simplify(H[i, j]) == 0


def test_m1_splitting_is_six_fifths_s_squared_exactly() -> None:
    """p shell: H is diagonal for every bend, and E_x - E_y = (6/5) w2 s^2.
    With attractive charges (w2 < 0) the in-plane p_x drops (Walsh).
    Rejects: splitting linear in the bend; a wrong coefficient; sigma-pi
    mixing inside the p shell."""
    H = shell_matrix(P_SHELL, {2: w2})
    assert H == sp.diag(*[H[i, i] for i in range(3)]), "p shell mixes under the symmetric bend"
    assert sp.expand(H[0, 0] - H[1, 1] - Rational(6, 5) * w2 * s**2) == 0


@pytest.mark.parametrize("L", [1, 2, 3])
def test_direct_splitting_opens_at_s_to_the_2L(L: int) -> None:
    """Direct coupling of +L and -L needs multipole 2L: zero for every
    k < 2L, and a single s^{2L} term at k = 2L.
    Rejects: an |m| = L pair breaking before order s^{2L}."""
    for k in range(0, 2 * L, 2):
        assert direct_splitting(L, k) == 0
    top = sp.Poly(direct_splitting(L, 2 * L), s)
    assert top.monoms() == [(2 * L,)]
    assert top.coeffs()[0] > 0  # in-plane partner on the charge side


def test_direct_coefficients_are_rational() -> None:
    """The angular prices: 6/5 (|m| = 1, k = 2) and 5/6 (|m| = 2, k = 4).
    Rejects: any non-rational or mis-normalised angular coefficient."""
    assert direct_splitting(1, 2) == Rational(6, 5) * s**2
    assert direct_splitting(2, 4) == Rational(5, 6) * s**4


def test_m2_relay_through_sigma_sets_the_s4_coefficient() -> None:
    """d shell: x^2-y^2 couples to z^2 at O(s^2), so the eigenvalue
    splitting of the |m| = 2 pair has NO s^2 term and an s^4 coefficient
    5/6 w4 - 9 (2w2 - 5w4)^2 / (14 (12w2 + 5w4)).
    Rejects: the m = 2 pair opening at s^2; a coefficient that ignores the
    relay through the m = 0 orbital."""
    H = shell_matrix(D_SHELL, {2: w2, 4: w4})
    names = list(D_SHELL)
    ia, ixy, iz = names.index('x2-y2'), names.index('xy'), names.index('z2')
    # the coupling exists and is O(s^2)
    off = sp.Poly(H[ia, iz], s)
    assert min(m[0] for m in off.monoms()) == 2
    # xy is alone in its symmetry class; x2-y2/z2 form a 2x2 block
    for j in range(5):
        if j != ixy:
            assert H[ixy, j] == 0
    # x2-y2 branch of the 2x2 block, lam = H_aa + delta, from the exact
    # characteristic equation delta = H_az^2 / (H_aa - H_zz + delta),
    # solved as a power series in s by fixed-point iteration.
    assert sp.simplify(H[ia, ia].subs(s, 0) - H[iz, iz].subs(s, 0)) != 0
    delta = sp.Integer(0)
    for _ in range(3):
        delta = sp.series(H[ia, iz] ** 2 / (H[ia, ia] - H[iz, iz] + delta),
                          s, 0, 6).removeO()
    split = sp.series(H[ia, ia] + delta - H[ixy, ixy], s, 0, 6).removeO()
    split = sp.expand(sp.simplify(split))
    assert sp.simplify(split.coeff(s, 1)) == 0
    assert sp.simplify(split.coeff(s, 2)) == 0
    assert sp.simplify(split.coeff(s, 3)) == 0
    expected = Rational(5, 6) * w4 - 9 * (2 * w2 - 5 * w4) ** 2 / (14 * (12 * w2 + 5 * w4))
    assert sp.simplify(split.coeff(s, 4) - expected) == 0
    # independent numeric route: diagonalise at small s and divide by s^4
    vals = {w2: -1.0, w4: -0.3}
    Hn = sp.lambdify(s, H.subs(vals), 'numpy')
    for sv in (2e-2, 1e-2):
        ev = np.linalg.eigvalsh(np.array(Hn(sv), dtype=float))
        e_a = min(ev, key=lambda e: abs(e - float(H[ia, ia].subs(vals).subs(s, sv))))
        e_xy = float(H[ixy, ixy].subs(vals).subs(s, sv))
        ratio = (e_a - e_xy) / sv**4
        assert abs(ratio - float(expected.subs(vals))) < 5e-3
