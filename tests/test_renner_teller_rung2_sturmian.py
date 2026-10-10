"""Bending breaks m in the Coulomb-Sturmian basis (rung 2S of the m-seam
probe, PI-directed 2026-10-10).  Same physics as
tests/test_renner_teller_rung2.py (one electron, central Z_c, two protons at
distance R in directions (s, 0, +-c)), but posed in Avery's way: one shared
exponent p0 for every basis function, and the Sturmian secular equation in
place of a Galerkin energy.

Derivation (scaled rho = p0 r, S = p0 R, E = -p0^2/2).  With Coulomb
Sturmians (-1/2 lap + 1/2) S_n = (n/rho) S_n normalised <S_i|1/rho|S_j> =
delta_ij / n_i, inserting psi = sum c_j S_j into the Schroedinger equation
and projecting gives the STANDARD Hermitian eigenproblem

    T(S) c = p0 c,   T(S) = diag(Z_c / n) + W(S),
    W_ij = <S_i| sum_protons 1/|rho - S n_m| |S_j>,

so the eigenvalues are p0 at fixed S (Shibuya-Wulfman's fixed-S posing, one
centre); fixed R is the self-consistent p0 = lambda(p0 R).  The multipole
expansion terminates by Gaunt; the radial products are poly(rho) e^{-2 rho},
so every matrix element is a closed form in Q(S) + Q(S) e^{-2S}: ONE
exponential, against the several e^{-Z(1/n+1/n')R} of the hydrogenic basis.

Guardrail (CLAUDE.md S3.5, Papers 8-9), acknowledged by the PI 2026-10-10:
this IS a shared-p0 single-centre Sturmian basis for a molecule.  The
structural theorem's R-independence is for the single-n D-matrix form; here
the cross-centre term is the exact multipole W(S) with cross-n couplings, so
the spectrum depends on S (pinned below), and the question is the bending law
at fixed R, not binding.
"""
from __future__ import annotations

import importlib.util
import os
from functools import lru_cache
from math import factorial
from typing import List, Tuple

import mpmath as mp
import pytest
import sympy as sp
from sympy import Rational

_here = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    '_rt2', os.path.join(_here, 'test_renner_teller_rung2.py'))
rt2 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(rt2)

rho = sp.Symbol('rho', positive=True)
Ss = rt2.Rs            # the boundary symbol; here it means S = p0 R
s = rt2.s
IN, OUT = rt2.IN, rt2.OUT


@lru_cache(maxsize=None)
def sturmian(n: int, l: int) -> sp.Expr:
    """Coulomb Sturmian in rho, normalised so int S^2 (1/rho) rho^2 drho = 1/n."""
    f = rho ** l * sp.exp(-rho) * sp.assoc_laguerre(n - l - 1, 2 * l + 1, 2 * rho)
    nrm = sp.integrate(sp.expand(f ** 2 * rho), (rho, 0, sp.oo))
    return sp.simplify(f / sp.sqrt(nrm * n))


@lru_cache(maxsize=None)
def radial_s(na: int, la: int, nb: int, lb: int, k: int) -> sp.Expr:
    """int S_a S_b rho_<^k / rho_>^{k+1} rho^2 drho with boundary S, closed form."""
    prod = sp.expand(sturmian(na, la) * sturmian(nb, lb) * sp.exp(2 * rho))
    total = sp.Integer(0)
    for (p,), coef in sp.Poly(sp.simplify(prod), rho).terms():
        total += coef * rt2._lower_inc(p + 2 + k, sp.Integer(2)) / Ss ** (k + 1)
        assert p + 1 - k >= 0
        total += coef * Ss ** k * rt2._upper_inc(p + 1 - k, sp.Integer(2))
    return sp.simplify(total)


def t_series(n_max: int, Zc: int = 1, l_max: int = 2) -> Tuple[List, List[sp.Matrix]]:
    """[T0..T4], T(S, s) = sum_j s^j T_j(S); eigenvalues are p0."""
    B = rt2.basis(n_max, l_max)
    N = len(B)
    T = [sp.zeros(N, N) for _ in range(5)]
    for i, (na, la, a) in enumerate(B):
        T[0][i, i] += Rational(Zc, na)
        for j, (nb, lb, b) in enumerate(B):
            for k in range(abs(la - lb), la + lb + 1):
                ang = rt2.angular((la, a), (lb, b), k)
                if ang == 0:
                    continue
                rad = radial_s(na, la, nb, lb, k)
                for (p,), coef in sp.Poly(ang, s).terms():
                    if p <= 4:
                        T[p][i, j] += rad * coef
    return B, T


def neg_numeric(T: List[sp.Matrix], Sv) -> List[mp.matrix]:
    """-T at fixed S (so rung-2's 'lowest' = largest p0 = most bound)."""
    return [-m for m in rt2.numeric(T, Sv)]


def fixed_r_p0(B, T_fns, R, sv, par, labels, p_guess) -> mp.mpf:
    """Self-consistent p0 = lambda(p0 R) on the branch that is pure on
    `labels` at s = 0, tracked by continuity in s."""
    def lam(p, target):
        Hn = [mp.matrix(f(p * R)) for f in T_fns]
        E, _, _ = rt2.block_eigs(B, [-h for h in Hn], sv, par)
        return -min(E, key=lambda e: abs(-e - target))
    def lam0(p):
        Hn = [-mp.matrix(f(p * R)) for f in T_fns]
        e, _ = rt2.pair_state(B, Hn, par, labels)
        return -e
    p0 = mp.findroot(lambda p: lam0(p) - p, p_guess)
    if sv == 0:
        return p0
    return mp.findroot(lambda p: lam(p, p0) - p, p0)


R_TEST = 2


@pytest.fixture(scope='module')
def model3():
    mp.mp.dps = 50
    B, T = t_series(3)
    return B, T


def test_sturmian_radials_vs_quadrature() -> None:
    """Closed forms == mpmath quadrature at S = 1.3; normalisation
    <S|1/rho|S> = 1/n.  Rejects: a wrong Sturmian radial integral."""
    mp.mp.dps = 30
    Sv = mp.mpf('1.3')
    for (na, la, nb, lb, k) in [(2, 1, 2, 1, 2), (2, 0, 3, 1, 1), (3, 2, 2, 0, 2), (3, 1, 1, 0, 1)]:
        closed = mp.mpf(sp.N(radial_s(na, la, nb, lb, k).subs(Ss, sp.Rational(13, 10)), 30))
        f = sp.lambdify(rho, sturmian(na, la) * sturmian(nb, lb) * rho**2, 'mpmath')
        quad = mp.quad(lambda t: f(t) * t**k / Sv ** (k + 1), [0, Sv]) \
            + mp.quad(lambda t: f(t) * Sv**k / t ** (k + 1), [Sv, mp.inf])
        assert abs(closed - quad) < mp.mpf(10) ** -20
    for n, l in [(1, 0), (2, 1), (3, 2), (3, 0)]:
        assert sp.simplify(sp.integrate(sturmian(n, l) ** 2 * rho, (rho, 0, sp.oo)) - Rational(1, n)) == 0


def test_single_exponential_price(model3) -> None:
    """Every T entry is A(S) + B(S) e^{-2S} with A, B rational in S: one
    transcendental for the whole basis, no pi.  Rejects: a stray exponential
    or a pi surviving the angular normalisation."""
    B, T = model3
    E = sp.Symbol('E', positive=True)
    for M in T:
        for e in set(M):
            if e == 0:
                continue
            assert not e.has(sp.pi)
            ee = sp.expand(sp.powsimp(sp.expand(e).subs(sp.exp(2 * Ss), 1 / E).subs(sp.exp(-2 * Ss), E)))
            assert not ee.has(sp.exp)
            assert sp.Poly(sp.together(ee) * sp.denom(sp.together(ee)), E).degree() <= 1


def test_guardrail_escape_spectrum_depends_on_S(model3) -> None:
    """The exact multipole W(S) makes the linear-geometry p0 spectrum move
    with S (Paper 8's R-independence is for the single-n D-matrix form).
    Rejects: an S-independent spectrum, i.e. the theorem's setting."""
    B, T = model3
    e1, _ = rt2.pair_state(B, neg_numeric(T, 1), IN, ('x',))
    e2, _ = rt2.pair_state(B, neg_numeric(T, 3), IN, ('x',))
    assert abs(e1 - e2) > mp.mpf('1e-3')


def test_symmetry_and_linear_degeneracy(model3) -> None:
    """(y, z)-parity blocks exact; pi and delta pairs degenerate at s = 0.
    Rejects: a symmetry leak in the Sturmian cross-centre term."""
    B, T = model3
    for M in T:
        for i, bi in enumerate(B):
            for j, bj in enumerate(B):
                if rt2.parity(bi[2]) != rt2.parity(bj[2]):
                    assert M[i, j] == 0
    Tn = neg_numeric(T, 1)
    for lin, lout in ((('x',), ('y',)), (('x2-y2',), ('xy',))):
        a, _ = rt2.pair_state(B, Tn, IN, lin)
        b, _ = rt2.pair_state(B, Tn, OUT, lout)
        assert abs(a - b) < mp.mpf(10) ** -40


def test_direct_m1_part_is_rung1_price(model3) -> None:
    """2p_x vs 2p_y s^2 diagonal of T = +(6/5) radial_s(2p,2p,k=2;S),
    the same angular price as rungs 1 and 2.  Rejects: a direct term that
    does not reduce to rung 1."""
    B, T = model3
    ix = B.index((2, 1, 'x'))
    iy = B.index((2, 1, 'y'))
    assert sp.simplify(T[2][ix, ix] - T[2][iy, iy] - Rational(6, 5) * radial_s(2, 1, 2, 1, 2)) == 0


def test_powers_at_fixed_S(model3) -> None:
    """Fixed S = 1: the p0 splitting of pi opens at s^2 (second-order theory
    == extrapolated exact eigenvalues, relay nonzero) and of delta at s^4
    (s^2 coefficient zero to 1e-35).  Rejects: pi linear in the bend; delta
    breaking at s^2; a pi coefficient without the relay."""
    B, T = model3
    Tn = neg_numeric(T, 1)
    d_in, r_in = rt2.second_order(B, Tn, IN, ('x',))
    d_out, r_out = rt2.second_order(B, Tn, OUT, ('y',))
    a2 = (d_in + r_in) - (d_out + r_out)
    assert abs(r_in - r_out) > mp.mpf('1e-3') * abs(a2)
    e0x, _ = rt2.pair_state(B, Tn, IN, ('x',))
    e0y, _ = rt2.pair_state(B, Tn, OUT, ('y',))
    def split(sv):
        return rt2.tracked_level(B, Tn, sv, IN, e0x) - rt2.tracked_level(B, Tn, sv, OUT, e0y)
    h = mp.mpf('1e-4')
    q1, q2 = split(h) / h**2, split(h / 2) / (h / 2) ** 2
    assert abs((4 * q2 - q1) / 3 - a2) < mp.mpf(10) ** -12 * max(1, abs(a2))
    dd = [rt2.second_order(B, Tn, par, lab) for par, lab in ((IN, ('x2-y2',)), (OUT, ('xy',)))]
    assert abs(sum(dd[0]) - sum(dd[1])) < mp.mpf(10) ** -35


def test_fixed_R_energy_split_matches_hydrogenic_power_law(model3) -> None:
    """Fixed R = 2 (self-consistent p0): E_pi_x - E_pi_y over s^2 is finite
    and stable between s = 2e-3 and 1e-3 (s^2 law at fixed R), and the
    branch at s = 0 is a genuine fixed point.  Rejects: a power law that
    holds only at fixed S."""
    B, T = model3
    T_fns = [sp.lambdify(Ss, M, 'mpmath') for M in T]
    p0 = fixed_r_p0(B, T_fns, R_TEST, 0, IN, ('x',), mp.mpf('0.6'))
    Tn = neg_numeric(T, p0 * R_TEST)
    e, _ = rt2.pair_state(B, Tn, IN, ('x',))
    assert abs(-e - p0) < mp.mpf(10) ** -30
    def dE(sv):
        px = fixed_r_p0(B, T_fns, R_TEST, sv, IN, ('x',), p0)
        py = fixed_r_p0(B, T_fns, R_TEST, sv, OUT, ('y',), p0)
        return -(px**2 - py**2) / 2
    h = [mp.mpf('2e-3'), mp.mpf('1e-3'), mp.mpf('5e-4')]
    q = [dE(v) / v**2 for v in h]
    rich1, rich2 = (4 * q[1] - q[0]) / 3, (4 * q[2] - q[1]) / 3
    assert abs(rich1) > mp.mpf('1e-4')
    assert abs(rich1 - rich2) < mp.mpf('1e-6') * abs(rich1)


@pytest.mark.slow
def test_hydrogenic_bound_set_is_incomplete_here() -> None:
    """Referee: a VARIATIONAL Galerkin in the fixed-k Sturmian p_x functions
    (H = k^2 (I - O/2) - k T(kR), H c = E O c; an upper bound) at n_max = 5,
    k = 1.2 lies within 5e-4 of the Sturmian secular pi energy, while the
    hydrogenic bound-state basis of rung 2 (n_max = 5) is > 0.1 Ha above it:
    the diffuse bound orbitals of Z_c = 1 cannot represent the compact pi
    state, so rung 2's R = 2 coefficients are not physical numbers.
    Rejects: treating the hydrogenic bound set as converged; a secular
    solution that is not a stationary point of the variational family."""
    mp.mp.dps = 40
    B, T = t_series(5, l_max=1)
    idx = [i for i, b in enumerate(B) if b[2] == 'x']
    O = mp.matrix(len(idx), len(idx))
    for a, i in enumerate(idx):
        for b, j in enumerate(idx):
            O[a, b] = mp.mpf(sp.N(sp.integrate(sturmian(B[i][0], 1) * sturmian(B[j][0], 1)
                                               * rho**2, (rho, 0, sp.oo)), 40))
    k = mp.mpf('1.2')
    Tk = mp.matrix(sp.lambdify(Ss, T[0], 'mpmath')(k * R_TEST))
    H = k**2 * (mp.eye(len(idx)) - O / 2) - k * mp.matrix(
        [[Tk[i, j] for j in idx] for i in idx])
    Li = mp.inverse(mp.cholesky(O))
    Ev = mp.eigsy(Li * H * Li.T)[0]
    e_gal = min(Ev[t] for t in range(len(idx)))
    T_fns = [sp.lambdify(Ss, M, 'mpmath') for M in T]
    p0 = fixed_r_p0(B, T_fns, R_TEST, 0, IN, ('x',), mp.mpf('1.2'))
    e_sec = -p0**2 / 2
    assert abs(e_sec - e_gal) < mp.mpf('5e-4')
    Bh, Hh = rt2.hamiltonian_series(5)
    e_hyd, _ = rt2.pair_state(Bh, rt2.numeric(Hh, R_TEST), IN, ('x',))
    assert e_hyd > e_gal + mp.mpf('0.1')
