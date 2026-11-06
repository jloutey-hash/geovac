"""Paper 63 sec:dictionary -- every exponential-orbital integral is a derivative of a 1s Yukawa integral.

Guards for Lemma lem:orbital, Proposition prop:rung (eq:euler_weights) and the MEASURED checks that
follow it.  Each test names the wrong answer it rejects.

  euler identities      Delta^{-(2m+1)} = (1 + 2/(2m-1) c d/dc) Delta^{-(2m-1)}, m = 1..4, and
                        L^2 Delta^{-j} = (Delta^{-(j-2)} - Delta^{-j})/c, exactly; the coefficient 1
                        in place of 2 fails
  l = 1 inversion       Y_10(r-A) e^{-zeta s} = -zeta^{-1} [T_2 + zeta^{-1} T_1],
                        T_m = Y_10(-grad_A) (-d/dzeta)^m g_zeta, exactly in coordinates
  l = 2 inversion       Y_20(r-A) e^{-zeta s} = zeta^{-2} U_3 + 3 zeta^{-3} U_2 + 3 zeta^{-4} U_1, exactly
  rational transform    the radial Hankel transform of r^{n-1} e^{-zeta r} Y_l0 is p^l times a rational
                        function of p^2, equal to (-d/dzeta)^{n-l} [2^l l! p^l/(p^2+zeta^2)^{l+1}],
                        exactly, for (n,l) up to (4,3)
  two densities         int_0^inf Delta_1^{-3} Delta_2^{-3} dk and int Delta_1^{-5} Delta_2^{-5} dk equal
                        the Euler operators of eq:euler_weights applied to K(m)/sqrt(c_max), to 1e-25,
                        and the operator images contain only K and E (the Legendre rung); the operator
                        (1 + c d/dc) fails
  three densities       the two-dimensional fingerprint integral with Delta^{-3} (and Delta^{-5}) weights
                        equals the tree identity with the Euler-operator images of the 1s weight,
                        (1/pi) int prod_v w_v(x) dx with w = 2 x c^{-1} K_1(x/sqrt c)
                        (and (4/3) x c^{-1} K_1 + (2/3) x^2 c^{-3/2} K_0), at three scale triples; the
                        1s weight 2 c^{-1/2} K_0 in place of w fails
"""
from __future__ import annotations

import mpmath as mp
import numpy as np
import pytest
import sympy as sp
from scipy import integrate, special

# ----------------------------------------------------------------------------- eq:euler_weights
c, L = sp.symbols("c L", positive=True)
Delta = sp.sqrt(c * L**2 + 1)


def _euler(expr, coeff):
    return sp.simplify(expr + coeff * c * sp.diff(expr, c))


@pytest.mark.parametrize("m", [1, 2, 3, 4])
def test_euler_identity_raises_the_power(m):
    lhs = _euler(Delta ** (-(2 * m - 1)), sp.Rational(2, 2 * m - 1))
    assert sp.simplify(lhs - Delta ** (-(2 * m + 1))) == 0
    # the wrong coefficient (1 instead of 2/(2m-1)) is not the identity
    wrong = _euler(Delta ** (-(2 * m - 1)), 1)
    assert sp.simplify(wrong - Delta ** (-(2 * m + 1))) != 0


@pytest.mark.parametrize("j", [3, 5, 7])
def test_momentum_square_lowers_the_power(j):
    assert sp.simplify(L**2 * Delta ** (-j) - (Delta ** (-(j - 2)) - Delta ** (-j)) / c) == 0


# ----------------------------------------------------------------------------- Lemma lem:orbital, position space
x, y, z, Ax, Ay, Az = sp.symbols("x y z A_x A_y A_z", real=True)
zeta = sp.symbols("zeta", positive=True)
s = sp.sqrt((x - Ax) ** 2 + (y - Ay) ** 2 + (z - Az) ** 2)
g = sp.exp(-zeta * s) / s


def _dzeta(f, m):
    for _ in range(m):
        f = -sp.diff(f, zeta)
    return f


def test_l1_solid_harmonic_times_exponential_is_a_derivative_of_the_yukawa():
    # Y_10(r - A) = z - A_z ; Y_10(-grad_A) = -d/dA_z
    T = {m: -sp.diff(_dzeta(g, m), Az) for m in (1, 2)}
    target = (z - Az) * sp.exp(-zeta * s)
    combo = -(T[2] + T[1] / zeta) / zeta
    assert sp.simplify(combo - target) == 0
    # the sign-flipped combination is not it
    assert sp.simplify((T[2] - T[1] / zeta) / zeta - target) != 0


POINTS = [{x: sp.Rational(3, 2), y: -1, z: sp.Rational(1, 3), Ax: sp.Rational(1, 4), Ay: 2, Az: -sp.Rational(5, 2)},
          {x: 2, y: sp.Rational(1, 2), z: -3, Ax: -1, Ay: sp.Rational(1, 5), Az: 1},
          {x: -sp.Rational(1, 2), y: sp.Rational(7, 3), z: 1, Ax: 1, Ay: -1, Az: sp.Rational(1, 2)}]


def test_l2_solid_harmonic_times_exponential_is_a_derivative_of_the_yukawa():
    # Y_20(r - A) = 3 (z - A_z)^2 - s^2 ;  Y_20(-grad_A) = 2 d^2/dA_z^2 - d^2/dA_x^2 - d^2/dA_y^2
    # (checked exactly at rational coordinate points with zeta symbolic: the generic simplification
    #  of the second derivatives is minutes long and adds nothing)
    def Y20_op(f):
        return 2 * sp.diff(f, Az, 2) - sp.diff(f, Ax, 2) - sp.diff(f, Ay, 2)

    U = {m: Y20_op(_dzeta(g, m)) for m in (1, 2, 3)}
    target = (3 * (z - Az) ** 2 - s**2) * sp.exp(-zeta * s)
    combo = U[3] / zeta**2 + 3 * U[2] / zeta**3 + 3 * U[1] / zeta**4
    for pt in POINTS:
        assert sp.simplify((combo - target).subs(pt)) == 0
    assert sp.simplify((U[3] / zeta**2 - target).subs(POINTS[0])) != 0


# ----------------------------------------------------------------------------- Lemma lem:orbital, momentum space
r, p = sp.symbols("r p", positive=True)


def _hankel(n, l):
    """sqrt-free radial transform  int_0^inf r^2 j_l(p r) r^(n-1) e^(-zeta r) dr, done termwise in sin/cos."""
    jl = sp.expand_func(sp.jn(l, p * r))
    integrand = sp.expand(r ** (n + 1) * sp.exp(-zeta * r) * jl)
    return sp.simplify(sp.integrate(integrand, (r, 0, sp.oo), conds="none"))


@pytest.mark.parametrize("n,l", [(1, 0), (2, 0), (2, 1), (3, 1), (3, 2), (4, 3)])
def test_transform_of_slater_orbital_is_p_to_the_l_times_rational(n, l):
    direct = _hankel(n, l)
    # seed = transform of r^l e^{-zeta r} Y_l0 / r  (the l-th Yukawa); the orbital r^{n-1} e^{-zeta r} Y_l0
    # is its (n - l)-th zeta-derivative, since s^j e^{-zeta s} = (-d/dzeta)^{j+1} g with j = n - 1 - l
    seed = 2**l * sp.factorial(l) * p**l / (p**2 + zeta**2) ** (l + 1)
    predicted = _dzeta(seed, n - l)
    assert sp.simplify(direct - predicted) == 0
    ratio = sp.simplify(direct / p**l)
    assert ratio.is_rational_function(p), ratio


# ----------------------------------------------------------------------------- two densities: the Legendre rung
def _two_density_lhs(c1, c2, power, dps=40):
    with mp.workdps(dps):
        c1, c2 = mp.mpf(c1.p) / c1.q, mp.mpf(c2.p) / c2.q      # exact: mp.mpf(Rational) would pass through a double
        f = lambda k: (c1 * k * k + 1) ** (-mp.mpf(power) / 2) * (c2 * k * k + 1) ** (-mp.mpf(power) / 2)
        return mp.quad(f, [0, 1, 10, mp.inf])


C1, C2 = sp.symbols("c_1 c_2", positive=True)
PERIOD = sp.elliptic_k(1 - C2 / C1) / sp.sqrt(C1)      # c_1 > c_2 branch of eq:period


def _two_density_rhs(c1, c2, power, dps=40):
    """Operators of eq:euler_weights applied symbolically to K(m)/sqrt(c_max)."""
    assert c1 > c2
    expr = PERIOD
    ops = {3: [2], 5: [2, sp.Rational(2, 3)]}[power]
    for coeff in ops:
        expr = expr + coeff * C1 * sp.diff(expr, C1)
        expr = expr + coeff * C2 * sp.diff(expr, C2)
    funcs = {type(a) for a in expr.atoms(sp.Function)}
    assert funcs <= {sp.elliptic_k, sp.elliptic_e}, funcs       # the Legendre rung
    with mp.workdps(dps):
        return expr.subs({C1: sp.Rational(c1), C2: sp.Rational(c2)}).evalf(dps)


@pytest.mark.parametrize("c1,c2", [("23/100", "11/100"), ("1/4", "1/20"), ("1/5", "7/100")])
@pytest.mark.parametrize("power", [3, 5])
def test_two_density_weights_are_euler_images_of_the_period(c1, c2, power):
    with mp.workdps(40):
        lhs = _two_density_lhs(sp.Rational(c1), sp.Rational(c2), power)
        rhs = mp.mpf(str(_two_density_rhs(sp.Rational(c1), sp.Rational(c2), power)))
        assert abs(lhs - rhs) < mp.mpf(10) ** (-25) * abs(lhs), (lhs, rhs)


def test_two_density_wrong_operator_fails():
    c1, c2 = sp.Rational(23, 100), sp.Rational(11, 100)
    lhs = _two_density_lhs(c1, c2, 3)
    expr = PERIOD
    expr = expr + C1 * sp.diff(expr, C1)          # (1 + c d/dc), the wrong coefficient
    expr = expr + C2 * sp.diff(expr, C2)
    wrong = mp.mpf(str(expr.subs({C1: c1, C2: c2}).evalf(40)))
    assert abs(lhs - wrong) > mp.mpf("1e-3") * abs(lhs)


# ----------------------------------------------------------------------------- three densities: the tree identity with Euler images
def _w(power, cv, xv):
    """Position-space image of Delta^{-power} with scale c: the Euler operators applied to 2 c^{-1/2} K_0."""
    u = xv / np.sqrt(cv)
    if power == 1:
        return 2.0 * cv ** (-0.5) * special.k0(u)
    if power == 3:
        return 2.0 * xv / cv * special.k1(u)
    if power == 5:
        return (4.0 / 3.0) * xv / cv * special.k1(u) + (2.0 / 3.0) * xv**2 * cv ** (-1.5) * special.k0(u)
    raise ValueError(power)


def _tree_rhs(ca, cb, cc, power):
    f = lambda xv: _w(power, ca, xv) * _w(power, cb, xv) * _w(power, cc, xv) / np.pi
    val, err = integrate.quad(f, 0.0, np.inf, limit=400, epsabs=0, epsrel=1e-13)
    return val


def _chain_lhs(ca, cb, cc, power):
    """int_{R^2} dp dq Delta_b(p)^-pow Delta_c(q)^-pow Delta_a(p+q)^-pow (the chain of three, 1-D momenta)."""
    e = -power / 2.0
    f = lambda q, pp: (cb * pp * pp + 1) ** e * (cc * q * q + 1) ** e * (ca * (pp + q) ** 2 + 1) ** e
    # the integrand falls like |k|^(-3 power); the box [-R, R]^2 leaves a tail below 1e-14 relative
    R = 400.0
    val, err = integrate.dblquad(f, -R, R, -R, R, epsabs=0, epsrel=1e-11)
    return val


TRIPLES = [(0.21, 0.08, 0.13), (0.05, 0.25, 0.17), (0.11, 0.11, 0.24)]


@pytest.mark.parametrize("ca,cb,cc", TRIPLES)
@pytest.mark.parametrize("power", [3, 5])
def test_three_density_fibre_matches_tree_identity_with_euler_images(ca, cb, cc, power):
    lhs = _chain_lhs(ca, cb, cc, power)
    rhs = _tree_rhs(ca, cb, cc, power)
    assert abs(lhs - rhs) < 2e-9 * abs(lhs), (lhs, rhs)


def test_three_density_wrong_weight_fails():
    ca, cb, cc = TRIPLES[0]
    lhs = _chain_lhs(ca, cb, cc, 3)
    f = lambda xv: _w(1, ca, xv) * _w(1, cb, xv) * _w(1, cc, xv) / np.pi    # the 1s weight for a Delta^{-3} fibre
    wrong, _ = integrate.quad(f, 0.0, np.inf, limit=400)
    assert abs(lhs - wrong) > 1e-2 * abs(lhs)
