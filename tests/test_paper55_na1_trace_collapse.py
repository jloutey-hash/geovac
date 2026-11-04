"""Paper 55 Sec. "Diagonal-substrate collapse of the depth-2 joint Mellin".

Guards (group3 guard pass, 2026-11-03) for three theorems that were backed
only by debug drivers whose "bit-exact" checks compared a sum with itself:

* ``thm:na1_trace_functional_collapse`` (eq:na1_J_general_closed) --
  for ANY self-adjoint gamma (not commuting with D) the double Mellin
  transform of Tr(D^2 e^{-t1 D^2} gamma D e^{-t2 D^2}) depends only on
  s_tot = s1 + s2 and equals sum_i gamma~_ii lambda_i^{3 - 2 s_tot}.
  Independent route: 2-D numerical quadrature of the operator trace built
  from matrix exponentials (no eigendecomposition inside the integrand).
  Wrong answers excluded: (a) off-diagonal entries of gamma~ contribute;
  (b) J depends on s1, s2 separately.

* ``thm:na1_diagonal_collapse`` (eq:na1_collapse) -- the specialisation to
  the Camporesi-Higuchi spectrum with vertex parity.  See the docstring of
  ``test_ch_specialisation_needs_one_chirality``: with BOTH chiralities
  (+-lambda_n, as the paper's state labels |n, chi, m_l> indicate) the
  trace vanishes identically, and the identity J = M3^{gammaP}(s_tot - 1)
  holds per chirality (i.e. for |D|).  The test pins both facts.

* ``thm:m3_gammaP_closed_form`` (eq:m3_gammaP_closed) --
  M3^{gammaP}(s) = 2^{2s-3} (beta(2s-1) - beta(2s-3)), s = 2..6.
  Independent route: the alternating series summed directly at 90 dps
  versus the closed form evaluated through the Euler-number formula
  beta(2k+1) = (-1)^k E_{2k} pi^{2k+1} / (4^{k+1} (2k)!) (no Hurwitz zeta,
  no shared code with the summation side), agreement < 1e-70.
"""
from __future__ import annotations

import mpmath as mp
import numpy as np
import pytest
import sympy as sp
from scipy.linalg import expm
from scipy.special import gamma as Gamma


# ---------------------------------------------------------------------------
# thm:na1_trace_functional_collapse
# ---------------------------------------------------------------------------

def _random_pair(seed: int = 7):
    rng = np.random.default_rng(seed)
    n = 4
    # D with spectrum of both signs, |lambda| in [0.9, 1.6]
    lam = np.array([0.9, -1.1, 1.3, -1.6])
    Q, _ = np.linalg.qr(rng.normal(size=(n, n)))
    D = Q @ np.diag(lam) @ Q.T
    G = rng.normal(size=(n, n))
    gamma = (G + G.T) / 2  # generic self-adjoint, does NOT commute with D
    return D, gamma


_U_MAX = 9.0     # t = u^2 up to 81: e^{-81 * 0.81} is far below tolerance
_N_NODES = 300
_X, _W = np.polynomial.legendre.leggauss(_N_NODES)
_U = (_X + 1) * _U_MAX / 2
_WU = _W * _U_MAX / 2


def _J_quadrature(D: np.ndarray, gamma: np.ndarray, s1: float, s2: float) -> float:
    """Double Mellin transform of the operator trace, by 2-D quadrature.

    Substitution t = u^2 (dt = 2u du) on each axis, tensor Gauss-Legendre in
    (u1, u2); the heat kernels are matrix exponentials (scipy expm), so no
    eigendecomposition of D enters the integrand.
    """
    D2 = D @ D
    T = _U ** 2
    A = np.array([D2 @ expm(-t * D2) for t in T])          # (N, n, n)
    B = np.array([D @ expm(-t * D2) for t in T])
    w1 = _WU * 2 * _U * T ** (s1 - 1)
    w2 = _WU * 2 * _U * T ** (s2 - 1)
    Aw = np.einsum("a,aij->ij", w1, A)                      # sum over t1
    Bw = np.einsum("b,bij->ij", w2, B)                      # sum over t2
    val = np.trace(Aw @ gamma @ Bw)
    return float(val) / (Gamma(s1) * Gamma(s2))


def _J_closed(D: np.ndarray, gamma: np.ndarray, s_tot: float) -> float:
    lam, U = np.linalg.eigh(D)
    gt = U.T @ gamma @ U
    return float(np.sum(np.diag(gt) * np.sign(lam) * np.abs(lam) ** (3 - 2 * s_tot)))


def test_trace_functional_collapse_by_quadrature():
    D, gamma = _random_pair()
    assert np.linalg.norm(D @ gamma - gamma @ D) > 0.5  # genuinely non-commuting
    for (s1, s2) in [(1.0, 2.0), (2.0, 1.0), (1.5, 1.5)]:
        J = _J_quadrature(D, gamma, s1, s2)
        assert abs(J - _J_closed(D, gamma, s1 + s2)) < 1e-7, (s1, s2, J)


def test_offdiagonal_entries_would_be_detected():
    """Discrimination control: the alternative where off-diagonal entries of
    gamma~ also contribute differs from the quadrature by O(1)."""
    D, gamma = _random_pair()
    lam, U = np.linalg.eigh(D)
    gt = U.T @ gamma @ U
    wrong = float(np.sum(gt.sum(axis=1) * np.sign(lam) * np.abs(lam) ** (3 - 6)))
    J = _J_quadrature(D, gamma, 1.0, 2.0)
    assert abs(J - wrong) > 1e-2


def test_J_not_a_function_of_s1_alone():
    """J(1,2) = J(2,1) although s1 changes; and J(1,2) != J(1,3)."""
    D, gamma = _random_pair()
    a = _J_quadrature(D, gamma, 1.0, 2.0)
    b = _J_quadrature(D, gamma, 2.0, 1.0)
    c = _J_quadrature(D, gamma, 1.0, 3.0)
    assert abs(a - b) < 1e-7
    assert abs(a - c) > 1e-3


# ---------------------------------------------------------------------------
# thm:na1_diagonal_collapse -- Camporesi-Higuchi specialisation
# ---------------------------------------------------------------------------

def _ch_diag(n_cut: int, both_chiralities: bool):
    lam, par = [], []
    for n in range(n_cut + 1):
        g = (n + 1) * (n + 2)
        for sign in ((+1, -1) if both_chiralities else (+1,)):
            lam += [sign * (2 * n + 3) / 2] * g
            par += [(-1) ** n] * g
    return np.array(lam), np.array(par, dtype=float)


def _J_quad_diag(lam: np.ndarray, par: np.ndarray, s1: float, s2: float) -> float:
    """Quadrature of the trace on a diagonal substrate."""
    l2 = lam ** 2
    T = _U ** 2
    w1 = _WU * 2 * _U * T ** (s1 - 1)
    w2 = _WU * 2 * _U * T ** (s2 - 1)
    A = np.einsum("a,ai->i", w1, l2[None, :] * np.exp(-np.outer(T, l2)))
    B = np.einsum("b,bi->i", w2, lam[None, :] * np.exp(-np.outer(T, l2)))
    return float(np.sum(A * par * B)) / (Gamma(s1) * Gamma(s2))


def test_ch_specialisation_needs_one_chirality():
    """With both chiralities (+-lambda_n) the trace D^3 e^{-tD^2} gamma_P is odd
    in D and J vanishes identically; on one chirality (|D|) J equals the
    truncated M3^{gammaP}(s_tot - 1).  Paper prose writes D (both
    chiralities, states |n, chi, m_l>); the identity holds for |D|."""
    n_cut = 2
    lam_b, par_b = _ch_diag(n_cut, both_chiralities=True)
    assert abs(_J_quad_diag(lam_b, par_b, 1.5, 1.5)) < 1e-10
    lam_p, par_p = _ch_diag(n_cut, both_chiralities=False)
    for (s1, s2) in [(1.5, 1.5), (1.0, 2.0)]:
        stot = s1 + s2
        m3_trunc = sum((-1) ** n * (n + 1) * (n + 2) * ((2 * n + 3) / 2) ** (1 - 2 * (stot - 1))
                       for n in range(n_cut + 1))
        assert abs(_J_quad_diag(lam_p, par_p, s1, s2) - m3_trunc) < 1e-9


# ---------------------------------------------------------------------------
# thm:m3_gammaP_closed_form
# ---------------------------------------------------------------------------

def _beta_odd(k: int) -> mp.mpf:
    """beta(2k+1) from the Euler-number closed form (independent of zeta)."""
    E = sp.euler(2 * k)
    return mp.mpf(int((-1) ** k * E)) * mp.pi ** (2 * k + 1) / (
        mp.mpf(4) ** (k + 1) * mp.factorial(2 * k))


def _closed_form(s: int) -> mp.mpf:
    # beta(2s-1) = beta(2(s-1)+1), beta(2s-3) = beta(2(s-2)+1)
    return mp.mpf(2) ** (2 * s - 3) * (_beta_odd(s - 1) - _beta_odd(s - 2))


def _direct_sum(s: int) -> mp.mpf:
    def a(n):
        n = int(n)
        return (n + 1) * (n + 2) * (mp.mpf(2 * n + 3) / 2) ** (1 - 2 * s)
    return mp.nsum(lambda n: (-1) ** int(n) * a(n), [0, mp.inf], method="alternating")


@pytest.mark.parametrize("s", [2, 3, 4, 5, 6])
def test_m3_gammaP_closed_form_high_precision(s):
    mp.mp.dps = 90
    direct = _direct_sum(s)
    closed = _closed_form(s)
    assert abs(direct - closed) < mp.mpf(10) ** -70, (s, mp.nstr(direct - closed, 5))


def test_m3_gammaP_closed_form_rejects_neighbours():
    """The direct sum must not also match a neighbouring wrong closed form
    (sign of the second beta flipped; exponent of 2 shifted)."""
    mp.mp.dps = 50
    for s in (2, 3, 4):
        d = _direct_sum(s)
        plus = mp.mpf(2) ** (2 * s - 3) * (_beta_odd(s - 1) + _beta_odd(s - 2))
        shifted = mp.mpf(2) ** (2 * s - 2) * (_beta_odd(s - 1) - _beta_odd(s - 2))
        assert abs(d - plus) > 1e-3 and abs(d - shifted) > 1e-3
