"""Paper 24 sec:ks_frame and prop thm:asymmetry: the Kustaanheimo-Stiefel lift.

Backs: (1) the fibre-neutral (Q = 0) sector of the four-dimensional oscillator holds the hydrogen label
set, n^2 states at level N = 2(n - 1) and none at odd N, with l = 0..n-1 once each; (2) the SO(4)
Casimir 2(J+^2 + J-^2) = (n^2 - 1) on integer lattice data; (3) the projection is linear, Z = n p0
with p0 = 2 omega; (4) the lift's one pi, the fibre volume pi/4; (5) the spectral sums of both columns
carry pi^2 at even argument; plus the heat-trace remark of Paper 32 sec:coulomb_ho.
Record: docs/predictions/2026-10-19_ks_lift.md.

Independent of debug/ks_lift_door.py. That script pushed oscillator states DOWN to R^3; here the
Coulomb Sturmians are pulled UP to R^4 and the oscillator equation is checked on them. The Casimir
uses integer matrices with 2 J_z, not Fractions; the fibre volume uses Gaussian moments in R^4
against a spherical integral in R^3; the ladders use mpmath nsum.
Oscillator H = p_u^2/8 + 2 omega^2 |u|^2; z1 = u1 + i u2, z2 = u3 + i u4;
x1 + i x2 = 2 z1 conj(z2), x3 = |z1|^2 - |z2|^2, r = |z|^2.
"""
from __future__ import annotations

import itertools

import mpmath as mp
import numpy as np
import pytest
import sympy as sp


# ---------------------------------------------------------------- lattice side
def q_zero_basis(N: int) -> list[tuple[int, int, int, int]]:
    """Monomials A+^a A-^b B+^c B-^d of degree N with Q = (a - b) + (c - d) = 0."""
    return [m for m in itertools.product(range(N + 1), repeat=4)
            if sum(m) == N and (m[0] - m[1]) + (m[2] - m[3]) == 0]


def test_q_zero_sector_is_the_hydrogen_label_set():
    for N in range(0, 22):
        dim = len(q_zero_basis(N))
        if N % 2:
            assert dim == 0, N
        else:
            n = N // 2 + 1
            assert dim == n * n, (N, dim)


def test_l_content_once_each():
    # diagonal L_z = Jz+ + Jz-, with 2 Jz+ = a - c and 2 Jz- = -(b - d): multiplicity of m is n - |m|
    for n in range(1, 11):
        basis = q_zero_basis(2 * (n - 1))
        twice_m = [(a - c) - (b - d) for a, b, c, d in basis]
        for m in range(-(n - 1), n):
            assert twice_m.count(2 * m) == n - abs(m), (n, m)


def _op(basis, terms):
    """Integer matrix of sum(coef * X * d/dY) on the basis; terms = [(coef, X, Y)], X, Y in 0..3."""
    idx = {m: i for i, m in enumerate(basis)}
    M = np.zeros((len(basis), len(basis)), dtype=object)
    for j, m in enumerate(basis):
        for coef, X, Y in terms:
            if m[Y] == 0:
                continue
            new = list(m)
            new[Y] -= 1
            new[X] += 1
            M[idx[tuple(new)], j] += coef * m[Y]
    return M


def _four_casimir(basis, plus: bool):
    """4 J^2 = (2 Jz)^2 + 2 (J+ J- + J- J+), integer matrices."""
    A, Bv = (0, 2) if plus else (1, 3)          # (A+, B+) or (A-, B-)
    sign = 1 if plus else -1                    # SU(2)- acts as the conjugate doublet
    two_jz = _op(basis, [(sign, A, A), (-sign, Bv, Bv)])
    jp = _op(basis, [(sign, A, Bv)]) if plus else _op(basis, [(-1, Bv, A)])
    jm = _op(basis, [(sign, Bv, A)]) if plus else _op(basis, [(-1, A, Bv)])
    return two_jz.dot(two_jz) + 2 * (jp.dot(jm) + jm.dot(jp))


@pytest.mark.parametrize("n", [1, 2, 3, 4, 5])
def test_casimir_is_n_squared_minus_one(n):
    basis = q_zero_basis(2 * (n - 1))
    cp = _four_casimir(basis, True)
    cm = _four_casimir(basis, False)
    assert (cp == cm).all()                     # J+^2 = J-^2 on the Q = 0 level
    target = 2 * (n * n - 1) * np.eye(len(basis), dtype=int)
    assert (cp + cm == target).all()            # 2 (J+^2 + J-^2) = (n^2 - 1) 1, times 4/2


# ---------------------------------------------------------------- linear projection
u1, u2, u3, u4, w = sp.symbols("u1 u2 u3 u4 omega", positive=True)
U = (u1, u2, u3, u4)
Z1 = u1 + sp.I * u2
Z2 = u3 + sp.I * u4
XP = sp.expand(2 * Z1 * sp.conjugate(Z2))       # x1 + i x2
X3 = sp.expand(u1**2 + u2**2 - u3**2 - u4**2)
RR = sp.expand(u1**2 + u2**2 + u3**2 + u4**2)

# (n, l, harmonic polynomial r^l Y_lm in x)
STATES = [(1, 0, sp.Integer(1)), (2, 0, sp.Integer(1)), (2, 1, X3), (2, 1, XP),
          (3, 0, sp.Integer(1)), (3, 1, X3), (3, 2, X3 * XP), (3, 2, XP**2)]


def _osc_eigenvalue(n, l, harm, p0, nlevel=None):
    """(H F)/F for the Coulomb Sturmian (n, l) at exponent p0, pulled up to R^4."""
    poly = sp.expand(harm * sp.assoc_laguerre(n - l - 1, 2 * l + 1, 2 * p0 * RR))
    G = sp.exp(-p0 * RR)
    F = poly * G
    HF = -sum(sp.diff(F, v, 2) for v in U) / 8 + 2 * w**2 * RR * F
    return sp.simplify(sp.expand(sp.simplify(HF / G)) / poly)


@pytest.mark.parametrize("n,l,harm", STATES)
def test_sturmian_pulled_up_has_oscillator_eigenvalue_n_p0(n, l, harm):
    p0 = 2 * w
    E = _osc_eigenvalue(n, l, harm, p0)
    assert sp.simplify(E - n * p0) == 0, (n, l, E)   # E_osc = omega (N + 2) = Z = n p0


@pytest.mark.parametrize("n,l,harm", [(2, 0, sp.Integer(1)), (2, 1, X3), (3, 2, XP**2)])
def test_wrong_exponent_is_not_an_eigenstate(n, l, harm):
    E = _osc_eigenvalue(n, l, harm, 3 * w)
    assert sp.simplify(sp.diff(E, u1)) != 0 or sp.simplify(sp.diff(E, u3)) != 0


# ---------------------------------------------------------------- the price: pi / 4
def _gauss4(poly, c):
    """Exact integral of poly(u) exp(-c |u|^2) over R^4 by one-dimensional Gaussian moments."""
    total = sp.Integer(0)
    for mon, coef in sp.Poly(sp.expand(poly), *U).terms():
        term = coef
        for k in mon:
            if k % 2:
                term = 0
                break
            term *= sp.gamma(sp.Rational(k + 1, 2)) / c**sp.Rational(k + 1, 2)
        total += term
    return sp.simplify(total)


def _sph3(xpoly, c):
    """Integral of xpoly(x) exp(-c r) / r over R^3."""
    r, th, ph = sp.symbols("r theta phi", positive=True)
    x1s, x2s, x3s = sp.symbols("x1 x2 x3")
    sub = {x1s: r * sp.sin(th) * sp.cos(ph), x2s: r * sp.sin(th) * sp.sin(ph), x3s: r * sp.cos(th)}
    f = xpoly(x1s, x2s, x3s).subs(sub) * sp.exp(-c * r) / r * r**2 * sp.sin(th)
    return sp.simplify(sp.integrate(f, (ph, 0, 2 * sp.pi), (th, 0, sp.pi), (r, 0, sp.oo)))


FIBRE_CASES = [
    (lambda a, b, c: sp.Integer(1), sp.Integer(1)),
    (lambda a, b, c: c**2, X3**2),
    (lambda a, b, c: a**2 + b**2, sp.expand(XP * sp.conjugate(XP))),
    (lambda a, b, c: c**2 * sp.sqrt(a**2 + b**2 + c**2), X3**2 * RR),
]


@pytest.mark.parametrize("case", range(len(FIBRE_CASES)))
def test_r4_inner_product_is_pi_over_4_times_weighted_r3(case):
    c = sp.Rational(3, 2)
    xpoly, upoly = FIBRE_CASES[case]
    lhs = _gauss4(upoly, c)
    rhs = _sph3(xpoly, c)
    assert sp.simplify(lhs - sp.pi / 4 * rhs) == 0, (case, lhs, rhs)


# ---------------------------------------------------------------- spectral sums
def test_spectral_ladders_carry_pi_squared_on_both_sides():
    mp.mp.dps = 30
    pi = mp.pi
    cas2 = mp.nsum(lambda n: n**2 / (n**2 - 1)**2, [2, mp.inf])
    cas4 = mp.nsum(lambda n: n**2 / (n**2 - 1)**4, [2, mp.inf])
    euler4 = mp.nsum(lambda N: (N + 1) * (N + 2) / 2 / N**4, [1, mp.inf])
    energy4 = mp.nsum(lambda N: (N + 1) * (N + 2) / 2 / (N + mp.mpf(3) / 2)**4, [0, mp.inf])
    tol = mp.mpf(10)**-25
    assert abs(cas2 - (pi**2 / 12 + mp.mpf(1) / 16)) < tol
    assert abs(cas4 - (mp.mpf(5) / 256 - pi**2 / 96 + pi**4 / 720)) < tol
    assert abs(euler4 - (mp.zeta(2) + 3 * mp.zeta(3) + 2 * mp.zeta(4)) / 2) < tol
    assert abs(energy4 - (pi**2 / 4 - pi**4 / 48)) < tol


def test_ho_heat_trace_rational_for_the_operator_pi_for_its_square():
    # Paper 32 sec:coulomb_ho: sum g_N e^{-tN} = (1 - e^{-t})^{-3} with rational Laurent coefficients
    q, t = sp.symbols("q t")
    ser = sp.series(1 / (1 - q)**3, q, 0, 12).removeO()
    for N in range(12):
        assert ser.coeff(q, N) == sp.Rational((N + 1) * (N + 2), 2)
    lau = sp.series((1 - sp.exp(-t))**-3, t, 0, 4).removeO()
    coeffs = [lau.coeff(t, k) for k in range(-3, 4)]
    assert all(cf.is_rational for cf in coeffs) and coeffs[0] == 1
    # the functional of the square: sum g_N e^{-t (N + 3/2)^2} ~ sqrt(pi) / 8 * t^(-3/2)
    tt = 1e-6
    k = np.arange(0, 20000) + 1.5
    S = np.sum((k * k - 0.25) / 2 * np.exp(-tt * k * k))
    assert abs(S * tt**1.5 / (np.sqrt(np.pi) / 8) - 1) < 1e-4
