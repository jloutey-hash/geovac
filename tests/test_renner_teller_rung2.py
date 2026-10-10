"""Bending breaks m with real nuclei: one electron, three centres (rung 2 of
the m-seam probe, PI-directed 2026-10-10; rung 1 is
tests/test_renner_teller_m_splitting.py).

Model: one electron; a central nucleus Z_c at the origin; two protons at
distance R in directions (s, 0, +c), (s, 0, -c), c = sqrt(1 - s^2), s = sin of
half the bend (s = 0 is linear along z).  Basis: hydrogenic orbitals of Z_c
(n <= n_max, l <= l_max) with real harmonic angular parts.  The protons'
attraction is the multipole expansion

    -1/|r - R n_i| = -sum_k (r_<^k / r_>^{k+1}) P_k(n_i . rhat),

which TERMINATES inside the basis (Gaunt: k <= l + l'), so every matrix
element is exact: angular part a polynomial in s with rational (times
square-root normalisation) coefficients, radial part a closed form in
Q(R)[e^{-alpha R}], alpha = Z_c (1/n + 1/n').  No quadrature.

Guardrail (CLAUDE.md S3.5, Papers 8-9), acknowledged by the PI 2026-10-10:
this is a single-centre basis for a molecule.  It is outside the Sturmian
structural theorem, which concerns the shared-p0 single-n D-matrix form
(R-independent spectrum); here the cross-centre attraction is explicitly
R- and angle-dependent, and the question is not binding or R_eq but the
bending law of the +-m splittings at FIXED R, a property of the exact matrix.

Pinned (each test names the wrong answer it rejects):
  * the bend is exact in symmetry: (y-parity, z-parity) blocks never mix and
    each +-m pair is degenerate at s = 0;
  * the DIRECT part of the |m| = 1 splitting is rung 1's 6/5 times the exact
    radial weight -I_2(2p,2p;R), whose closed form carries only e^{-Z_c R}
    (no pi) and tends to the point-charge value -<r^2>/R^3;
  * the full |m| = 1 coefficient (second-order, exact matrix) also has a
    RELAY through the sigma states via the dipole term (k = 1), which rung 1
    could not have; it equals the extrapolated exact eigenvalue splitting;
  * the relay's size AND sign are set by where the sigma levels sit (FREE):
    with pure hydrogenic levels the accidental 2s/2p degeneracy makes it
    dominate at large R and reverse the ordering; a quantum-defect-like
    shift of the s levels tames it.  So the Walsh ordering is not forced
    by this model (which has no hydrogen-centred functions);
  * the |m| = 2 pair has zero s^2 coefficient and opens at s^4.
"""
from __future__ import annotations

from functools import lru_cache
from math import factorial
from typing import Dict, List, Tuple

import mpmath as mp
import pytest
import sympy as sp
from sympy import Rational, gamma
from sympy.physics.hydrogen import R_nl

x, y, z = sp.symbols('x y z', real=True)
s, c = sp.symbols('s c', real=True)
r, Rs = sp.symbols('r R', positive=True)

ANG = {  # real harmonic polynomials, l <= 2
    (0, 's'): sp.Integer(1),
    (1, 'x'): x, (1, 'y'): y, (1, 'z'): z,
    (2, 'x2-y2'): x**2 - y**2, (2, 'xy'): x * y, (2, 'xz'): x * z,
    (2, 'yz'): y * z, (2, 'z2'): 2 * z**2 - x**2 - y**2,
}


def solid_harmonic(l: int, m: int, kind: str) -> sp.Expr:
    """Real solid harmonic r^l P_l^m(z/r) {cos, sin}(m phi) as a polynomial:
    (x + i y)^m times r^{l-m} (d^m P_l / dt^m)(z / r)."""
    t = sp.Symbol('t')
    q = sp.Poly(sp.diff(sp.legendre(l, t), t, m), t)
    rad = sp.expand(sum(coef * z**j * (x**2 + y**2 + z**2) ** ((l - m - j) // 2)
                        for (j,), coef in q.terms()))
    w = sp.expand((x + sp.I * y) ** m)
    return sp.expand((sp.re(w) if kind == 'c' else sp.im(w)) * rad)


for _l in (3, 4):  # generated labels l{l}c{m} / l{l}s{m}; l <= 2 keep their names
    ANG[(_l, 'l%dc0' % _l)] = solid_harmonic(_l, 0, 'c')
    for _m in range(1, _l + 1):
        ANG[(_l, 'l%dc%d' % (_l, _m))] = solid_harmonic(_l, _m, 'c')
        ANG[(_l, 'l%ds%d' % (_l, _m))] = solid_harmonic(_l, _m, 's')


def sphere_int(poly: sp.Expr) -> sp.Expr:
    """Exact integral of a polynomial in (x, y, z) over the unit sphere."""
    total = sp.Integer(0)
    for (a, b, d), coef in sp.Poly(sp.expand(poly), x, y, z).terms():
        if a % 2 or b % 2 or d % 2:
            continue
        total += coef * 2 * gamma(Rational(a + 1, 2)) * gamma(Rational(b + 1, 2)) \
            * gamma(Rational(d + 1, 2)) / gamma(Rational(a + b + d + 3, 2))
    return sp.simplify(total)


def eliminate_c(expr: sp.Expr) -> sp.Expr:
    """Rewrite in s alone (c^2 = 1 - s^2); odd powers of c must cancel."""
    out = sp.Integer(0)
    for (p,), coef in sp.Poly(sp.expand(expr), c).terms():
        assert p % 2 == 0, "odd power of c: the two protons are not mirror images"
        out += coef * (1 - s**2) ** (p // 2)
    return sp.expand(out)


@lru_cache(maxsize=None)
def angular(a: Tuple[int, str], b: Tuple[int, str], k: int) -> sp.Expr:
    """sum over both protons of <Y_a| P_k(n_i . rhat) |Y_b>, normalised, poly in s."""
    op = sum(sp.legendre(k, s * x + sg * c * z) for sg in (1, -1))
    num = sphere_int(ANG[a] * ANG[b] * sp.expand(op))
    if num == 0:
        return sp.Integer(0)
    norm = sp.sqrt(sphere_int(ANG[a] ** 2) * sphere_int(ANG[b] ** 2))
    return eliminate_c(sp.simplify(num / norm))


def _lower_inc(j: int, alpha: sp.Expr) -> sp.Expr:
    """int_0^R r^j e^{-alpha r} dr, exact."""
    return factorial(j) / alpha ** (j + 1) * (
        1 - sp.exp(-alpha * Rs) * sum((alpha * Rs) ** i / factorial(i) for i in range(j + 1)))


def _upper_inc(j: int, alpha: sp.Expr) -> sp.Expr:
    """int_R^oo r^j e^{-alpha r} dr, exact (j >= 0)."""
    return factorial(j) / alpha ** (j + 1) * sp.exp(-alpha * Rs) * sum(
        (alpha * Rs) ** i / factorial(i) for i in range(j + 1))


@lru_cache(maxsize=None)
def radial(na: int, la: int, nb: int, lb: int, k: int, Zc: int) -> sp.Expr:
    """I_k = int R_a R_b r_<^k / r_>^{k+1} r^2 dr, closed form in R."""
    alpha = Rational(Zc, na) + Rational(Zc, nb)
    prod = sp.expand(R_nl(na, la, r, Zc) * R_nl(nb, lb, r, Zc) * sp.exp(alpha * r))
    poly = sp.Poly(sp.simplify(prod), r)
    total = sp.Integer(0)
    for (p,), coef in poly.terms():
        total += coef * (_lower_inc(p + 2 + k, alpha) / Rs ** (k + 1))
        assert p + 1 - k >= 0, "Gaunt termination should keep r-powers nonnegative"
        total += coef * Rs ** k * _upper_inc(p + 1 - k, alpha)
    return sp.simplify(total)


def basis(n_max: int, l_max: int = 2) -> List[Tuple[int, int, str]]:
    return [(n, l, lab) for n in range(1, n_max + 1)
            for (l, lab) in ANG if l < n and l <= l_max]


def parity(lab: str) -> Tuple[int, int]:
    """(y-parity, z-parity) of a real harmonic: exact symmetries of the bend."""
    p = ANG[next(key for key in ANG if key[1] == lab)]
    return (1 if sp.expand(p.subs(y, -y) - p) == 0 else -1,
            1 if sp.expand(p.subs(z, -z) - p) == 0 else -1)


def hamiltonian_series(n_max: int, Zc: int = 1, l_max: int = 2) -> Tuple[List, List[sp.Matrix]]:
    """[H0, H1, H2, H3, H4] with H(s) = sum_j s^j H_j, entries exact in R."""
    B = basis(n_max, l_max)
    N = len(B)
    H = [sp.zeros(N, N) for _ in range(5)]
    for i, (na, la, a) in enumerate(B):
        H[0][i, i] += -Rational(Zc**2, 2 * na**2)
        for j, (nb, lb, b) in enumerate(B):
            for k in range(abs(la - lb), la + lb + 1):
                ang = angular((la, a), (lb, b), k)
                if ang == 0:
                    continue
                rad = radial(na, la, nb, lb, k, Zc)
                for (p,), coef in sp.Poly(ang, s).terms():
                    if p <= 4:
                        H[p][i, j] += -rad * coef
    return B, H


def numeric(Hs: List[sp.Matrix], Rv) -> List[mp.matrix]:
    return [mp.matrix([[mp.mpf(sp.N(e.subs(Rs, Rv), 50)) for e in row] for row in M.tolist()])
            for M in Hs]


def block_eigs(B, Hn: List[mp.matrix], sv, par) -> List:
    idx = [i for i, b in enumerate(B) if parity(b[2]) == par]
    M = mp.matrix(len(idx), len(idx))
    for a, i in enumerate(idx):
        for b, j in enumerate(idx):
            M[a, b] = sum(Hn[p][i, j] * mp.mpf(sv) ** p for p in range(5))
    E, V = mp.eigsy(M)
    return [E[t] for t in range(len(idx))], V, idx


def pair_state(B, Hn, par, labels) -> Tuple[mp.mpf, mp.matrix]:
    """Lowest linear-geometry eigenvector of `par` living on `labels` only."""
    E, V, idx = block_eigs(B, Hn, 0, par)
    best = None
    for t in range(len(E)):
        w = sum(V[a, t] ** 2 for a, i in enumerate(idx) if B[i][2] in labels)
        if w > 1 - mp.mpf(10) ** -30 and (best is None or E[t] < best[0]):
            vec = mp.matrix(len(B), 1)
            for a, i in enumerate(idx):
                vec[i] = V[a, t]
            best = (E[t], vec)
    assert best is not None, "no pure state on %s" % (labels,)
    return best


def second_order(B, Hn, par, labels) -> Tuple[mp.mpf, mp.mpf]:
    """(direct, relay) s^2 coefficients for the state on `labels` in block `par`."""
    E0, psi = pair_state(B, Hn, par, labels)
    direct = (psi.T * Hn[2] * psi)[0]
    E, V, idx = block_eigs(B, Hn, 0, par)
    relay = mp.mpf(0)
    for t in range(len(E)):
        phi = mp.matrix(len(B), 1)
        for a, i in enumerate(idx):
            phi[i] = V[a, t]
        cpl = (phi.T * Hn[1] * psi)[0]
        if abs(E[t] - E0) < mp.mpf(10) ** -25:
            assert abs(cpl) < mp.mpf(10) ** -25
            continue
        relay += cpl ** 2 / (E0 - E[t])
    return direct, relay


def tracked_level(B, Hn, sv, par, E0) -> mp.mpf:
    E, _, _ = block_eigs(B, Hn, sv, par)
    return min(E, key=lambda e: abs(e - E0))


IN, OUT = (1, 1), (-1, 1)   # (y, z) parities of pi_x / pi_y and of x2-y2 / xy
R_TEST = 2


@pytest.fixture(scope='module')
def model3():
    mp.mp.dps = 50
    B, H = hamiltonian_series(3)
    return B, H, numeric(H, R_TEST)


def test_radial_closed_form_vs_quadrature() -> None:
    """Closed-form I_k (incomplete-gamma route) == direct mpmath quadrature.
    Rejects: a wrong radial weight."""
    mp.mp.dps = 30
    for (na, la, nb, lb, k) in [(2, 1, 2, 1, 2), (2, 0, 2, 1, 1), (3, 2, 3, 0, 2), (3, 1, 2, 1, 0)]:
        closed = mp.mpf(sp.N(radial(na, la, nb, lb, k, 1).subs(Rs, R_TEST), 30))
        f = sp.lambdify(r, R_nl(na, la, r, 1) * R_nl(nb, lb, r, 1) * r**2, 'mpmath')
        quad = mp.quad(lambda t: f(t) * t**k / mp.mpf(R_TEST) ** (k + 1), [0, R_TEST]) \
            + mp.quad(lambda t: f(t) * mp.mpf(R_TEST) ** k / t ** (k + 1), [R_TEST, mp.inf])
        assert abs(closed - quad) < mp.mpf(10) ** -20


def test_symmetry_blocks_and_linear_degeneracy(model3) -> None:
    """Every H_j respects (y, z) parity; at s = 0 pi_x/pi_y and the delta pair
    are degenerate to 1e-40.  Rejects: a symmetry leak in the cross-centre term."""
    B, H, Hn = model3
    for M in H:
        for i, bi in enumerate(B):
            for j, bj in enumerate(B):
                if parity(bi[2]) != parity(bj[2]):
                    assert M[i, j] == 0
    for lin, lout in ((('x',), ('y',)), (('x2-y2',), ('xy',))):
        e_in, _ = pair_state(B, Hn, IN, lin)
        e_out, _ = pair_state(B, Hn, OUT, lout)
        assert abs(e_in - e_out) < mp.mpf(10) ** -40


def test_direct_m1_part_is_rung1_price_times_exact_weight() -> None:
    """n_max = 2: the s^2 diagonal difference of 2p_x vs 2p_y is EXACTLY
    -(6/5) I_2(2p,2p;R); I_2 carries only e^{-Z_c R} (no pi) and
    R^3 I_2 -> <r^2>_{2p} = 30/Z_c^2.  Rejects: a direct term that does not
    reduce to rung 1; a pi or stray exponential in the price."""
    B, H = hamiltonian_series(2)
    labs = [b[2] for b in B]
    ix = [i for i, b in enumerate(B) if b[:2] == (2, 1) and b[2] == 'x'][0]
    iy = [i for i, b in enumerate(B) if b[:2] == (2, 1) and b[2] == 'y'][0]
    I2 = radial(2, 1, 2, 1, 2, 1)
    assert sp.simplify(H[2][ix, ix] - H[2][iy, iy] + Rational(6, 5) * I2) == 0
    assert not I2.has(sp.pi)
    # exponential content: exactly A(R) + B(R) e^{-R}, A and B rational in R, B != 0
    E = sp.Symbol('E', positive=True)
    I2e = sp.expand(sp.expand(I2).subs(sp.exp(Rs), 1 / E).subs(sp.exp(-Rs), E))
    I2e = sp.expand(sp.powsimp(I2e))
    assert not I2e.has(sp.exp)
    poly_E = sp.Poly(sp.together(I2e) * sp.denom(sp.together(I2e)), E)
    assert poly_E.degree() == 1 and set(poly_E.free_symbols) <= {Rs, E}
    assert sp.simplify(I2e.coeff(E, 1)) != 0 and sp.simplify(I2e.coeff(E, 0)) != 0
    assert sp.limit(I2 * Rs**3, Rs, sp.oo) == 30
    assert 'x' in labs and 'y' in labs


def test_m1_coefficient_has_a_relay_and_matches_exact_eigenvalues(model3) -> None:
    """n_max = 3, R = 2: a2(pi) = direct + relay from second-order theory on
    the exact matrix equals the Richardson-extrapolated exact eigenvalue
    splitting / s^2, and the relay is nonzero.
    Rejects: a coefficient that is rung 1's direct term alone; a splitting
    linear in the bend."""
    B, H, Hn = model3
    d_in, r_in = second_order(B, Hn, IN, ('x',))
    d_out, r_out = second_order(B, Hn, OUT, ('y',))
    a2 = (d_in + r_in) - (d_out + r_out)
    relay = r_in - r_out
    assert abs(relay) > mp.mpf(10) ** -3 * abs(a2)
    e0x, _ = pair_state(B, Hn, IN, ('x',))
    e0y, _ = pair_state(B, Hn, OUT, ('y',))
    def split(sv):
        return tracked_level(B, Hn, sv, IN, e0x) - tracked_level(B, Hn, sv, OUT, e0y)
    h = mp.mpf('1e-4')
    q1, q2 = split(h) / h**2, split(h / 2) / (h / 2) ** 2
    extrap = (4 * q2 - q1) / 3
    assert abs(extrap - a2) < mp.mpf(10) ** -12 * max(1, abs(a2))
    assert abs(split(h) / h) < mp.mpf(10) ** -3  # no linear term


def test_relay_is_set_by_level_structure_not_angles(model3) -> None:
    """The relay is FREE: at R = 16 the pure-hydrogenic 2s/2p degeneracy
    (accidental SO(4)) makes it dominate the direct term (~1/R vs ~1/R^3)
    and reverse the ordering; lifting the s levels by 0.1 Ha (a quantum-
    defect-like shift) cuts it more than tenfold and restores in-plane
    lowering.  The direct term is untouched by the shift.
    Rejects: a relay that is a fixed angular price; an ordering that is
    forced independent of the sigma level positions."""
    B, H, _ = model3
    rt_numeric = numeric(H, 16)
    def coeffs(shift):
        M = [m.copy() for m in rt_numeric]
        for i, b in enumerate(B):
            if b[1] == 0:
                M[0][i, i] += mp.mpf(shift)
        d_in, r_in = second_order(B, M, IN, ('x',))
        d_out, r_out = second_order(B, M, OUT, ('y',))
        return d_in - d_out, r_in - r_out
    d0, rel0 = coeffs(0)
    d1, rel1 = coeffs('-0.1')
    assert abs(d0 - d1) < mp.mpf(10) ** -40
    assert rel0 > 5 * abs(d0) and d0 + rel0 > 0
    assert abs(rel1) < abs(rel0) / 10 and d1 + rel1 < 0


def test_m2_pair_opens_at_s4(model3) -> None:
    """n_max = 3: the 3d delta pair has a2 = 0 (to 1e-35) and split/s^4 is
    finite and nonzero.  Rejects: the |m| = 2 pair breaking at s^2."""
    B, H, Hn = model3
    d_in, r_in = second_order(B, Hn, IN, ('x2-y2',))
    d_out, r_out = second_order(B, Hn, OUT, ('xy',))
    assert abs((d_in + r_in) - (d_out + r_out)) < mp.mpf(10) ** -35
    e0a, _ = pair_state(B, Hn, IN, ('x2-y2',))
    e0b, _ = pair_state(B, Hn, OUT, ('xy',))
    def split(sv):
        return tracked_level(B, Hn, sv, IN, e0a) - tracked_level(B, Hn, sv, OUT, e0b)
    q1, q2 = split(mp.mpf('1e-3')) / mp.mpf('1e-3') ** 4, split(mp.mpf('5e-4')) / mp.mpf('5e-4') ** 4
    assert abs(q1) > mp.mpf(10) ** -6
    assert abs(q1 - q2) < mp.mpf(10) ** -3 * abs(q1)
