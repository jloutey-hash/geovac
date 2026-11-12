"""Paper 38, Eq. (korovkin_rate) and the measured best constant of the bound d_GH <= gamma(a).

Claim (Theorem thm:main_unconditional, "Other amplitudes"): with theta = chi/2 and k = 2j+1 = 1..n,
gamma(a) = (4/pi) int_0^pi theta |sum_k a_k sin k theta|^2 dtheta / sum_k a_k^2   (dual-Coxeter),
and the sine taper a_k = sin(pi k/(n+1)) gives gamma <= 2 pi sin(pi/(2(n+1))) <= pi^2/(n+1) for every n.
Measured: n * min_a gamma(a) -> 2 lambda_2, lambda_2 the second Dirichlet eigenvalue of (-Delta)^(1/2) on (-1,1).

Wrong answers rejected (each named at its test):
- the reduction's normalisation (a factor-2 slip between dual-Coxeter and unit-S^3, or the j-window n vs (n+1)/2),
  caught against the repo's own gamma_rate on the Cesaro taper and an independent quadrature;
- the claim that any amplitude gives 1/n (the Fejer amplitude must break pi^2/(n+1) at large n);
- the unit-S^3 constant pi^2/(2(n+1)) quoted in the dual-Coxeter normalisation;
- lambda_1 (the even sector) or the asymptotic 7 pi/8 in place of lambda_2.
"""
import numpy as np
import pytest
from scipy.integrate import quad
from scipy.linalg import eigh


def form_matrix(n: int) -> np.ndarray:
    k = np.arange(1, n + 1, dtype=float)[:, None]
    l = np.arange(1, n + 1, dtype=float)[None, :]
    d = k - l
    odd = np.abs(np.rint(d)) % 2 == 1
    with np.errstate(divide="ignore"):
        T = np.where(odd, 1.0 / d**2 - 1.0 / (k + l) ** 2, 0.0)
    return np.pi * np.eye(n) - (4.0 / np.pi) * T


def gamma(a: np.ndarray) -> float:
    a = np.asarray(a, float)
    return float(a @ form_matrix(len(a)) @ a / (a @ a))


def gamma_quad(a: np.ndarray) -> float:
    """Independent route: integrate the definition directly."""
    k = np.arange(1, len(a) + 1)
    p2 = lambda t: np.dot(a, np.sin(k * t)) ** 2
    num = quad(lambda t: t * p2(t), 0, np.pi, limit=400)[0]
    return 4.0 / np.pi * num / float(np.dot(a, a))


def sine_taper(n: int) -> np.ndarray:
    return np.sin(np.pi * np.arange(1, n + 1) / (n + 1))


def test_reduction_matches_repo_gamma_rate_on_cesaro():
    """Rejects a normalisation or window slip in the reduction: the repo's gamma_rate (character sum,
    mpmath quadrature on [0, 2 pi]) must equal the matrix form on the Cesaro amplitude."""
    from geovac.central_fejer_su2 import gamma_rate
    for n in (5, 12):
        k = np.arange(1, n + 1)
        j = (k - 1) / 2
        a = np.sqrt(2 * j + 1) * (1 - 2 * j / (n + 1))
        assert abs(gamma(a) - float(gamma_rate(n, prec=30, use_cesaro=True))) < 1e-10


def test_matrix_matches_direct_quadrature():
    rng = np.random.default_rng(0)
    for n in (3, 8, 15):
        a = rng.standard_normal(n)
        assert abs(gamma(a) - gamma_quad(a)) < 1e-8


@pytest.mark.parametrize("n", list(range(1, 41)) + [100, 400])
def test_sine_taper_obeys_eq_korovkin_rate(n):
    g = gamma(sine_taper(n))
    assert g <= 2 * np.pi * np.sin(np.pi / (2 * (n + 1))) + 1e-12
    assert 2 * np.pi * np.sin(np.pi / (2 * (n + 1))) <= np.pi**2 / (n + 1) + 1e-15


def test_bound_is_not_the_unit_sphere_constant():
    """Rejects quoting pi^2/(2(n+1)) in the dual-Coxeter normalisation: the taper exceeds it."""
    n = 400
    assert gamma(sine_taper(n)) > np.pi**2 / (2 * (n + 1))


def test_fejer_amplitude_breaks_the_1_over_n_bound():
    """Rejects 'any amplitude gives pi^2/(n+1)': the Fejer amplitude of Def. central_fejer,
    a_j = sqrt(2j+1), i.e. a_k = sqrt(k) (gamma ~ (4/pi) log n / n), exceeds it at n = 400
    (n gamma = 11.72 vs 9.85), so the taper is doing the work."""
    n = 400
    assert gamma(np.sqrt(np.arange(1, n + 1))) > np.pi**2 / (n + 1)


def _lambda2_galerkin(M: int = 160) -> float:
    """(-Delta)^(1/2) on (-1,1): phi_m = sqrt(1-x^2) U_m, (-Delta)^(1/2) phi_m = (m+1) U_m; odd sector m odd."""
    ms = np.arange(1, M, 2)
    Q = 2 * M + 10
    th = np.arange(1, Q + 1) * np.pi / (Q + 1)
    w = np.pi / (Q + 1) * np.sin(th) ** 3          # Gauss-Chebyshev-2 weight times sqrt(1-x^2)
    U = np.sin(np.outer(ms + 1, th)) / np.sin(th)
    mass = (U * w) @ U.T
    stiff = np.diag((ms + 1) * np.pi / 2)
    return float(eigh(stiff, mass, eigvals_only=True)[0])


def test_best_constant_is_lambda2_of_the_cauchy_interval_operator():
    """Measured identification. Rejects lambda_1 = 1.1578 (even sector) and the asymptotic 7 pi/8 = 2.7489."""
    lam2 = _lambda2_galerkin()
    assert abs(lam2 - 2.7547547) < 2e-7            # Kwasnicki 2012 Table 2 brackets [2.75473, 2.75483]
    v = []
    for n in (400, 800):
        lmin = eigh(form_matrix(n), eigvals_only=True, subset_by_index=[0, 0])[0]
        v.append(n * lmin / 2)
    richardson = 2 * v[1] - v[0]
    assert v[0] < v[1] < lam2                       # rises towards it from below
    assert abs(richardson - lam2) < 2e-6
