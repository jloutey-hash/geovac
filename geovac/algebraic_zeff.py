"""
Algebraic Z_eff evaluation via Laguerre spectral expansion.

Replaces numerical quadrature (cumulative_trapezoid + CubicSpline) in the
Z_eff screening function with closed-form evaluation using the regularized
lower incomplete gamma function.

The core electron density n(r) is projected onto a Laguerre spectral basis:

    n(r) = (4*beta^3) * sum_k c_k * r^2 * L_k(2*beta*r) * exp(-2*beta*r)

where L_k is the Laguerre polynomial, beta is a decay parameter matched to
the density's exponential fall-off, and c_k are expansion coefficients
determined by least-squares projection from the grid density.

The cumulative electron count N_core(r) = integral_0^r n(r') dr' then becomes:

    N_core(r) = C_inf - exp(-X) * P(X)

where X = 2*beta*r, C_inf is a rational constant (the total electron count
in the infinite limit), and P(X) is a polynomial with rational coefficients
derived from Laguerre polynomial coefficients and factorials.

Transcendental content: a single exponential e^{-2*beta*r} (the Stieltjes seed
analog for the screening problem). All other operations are rational arithmetic
on Laguerre polynomial coefficients and factorials.

References:
  - Paper 13, Section XII (algebraic structure)
  - Paper 17, Section III (Z_eff screening usage in composed geometries)
  - Track H (algebraic Laguerre matrix elements pattern)
"""

import numpy as np
from scipy.special import eval_laguerre, factorial
from typing import Tuple, Optional
from math import comb as _comb


def _laguerre_coefficients(k: int) -> np.ndarray:
    """Return the monomial coefficients of L_k(x) = sum_j a_j x^j.

    L_k(x) = sum_{j=0}^k (-1)^j C(k,j) x^j / j!

    Parameters
    ----------
    k : int
        Laguerre polynomial degree.

    Returns
    -------
    coeffs : ndarray of shape (k+1,)
        coeffs[j] is the coefficient of x^j.
    """
    coeffs = np.zeros(k + 1)
    for j in range(k + 1):
        coeffs[j] = (-1)**j * _comb(k, j) / factorial(j, exact=True)
    return coeffs


def _lower_incomplete_gamma_int(a: int, x: np.ndarray) -> np.ndarray:
    """Compute the lower incomplete gamma function gamma(a, x) for integer a >= 1.

    gamma(a, x) = (a-1)! * [1 - e^{-x} * sum_{m=0}^{a-1} x^m / m!]

    This is exact (no numerical integration) and vectorized over x.

    Parameters
    ----------
    a : int
        Must be a positive integer.
    x : ndarray
        Evaluation points (must be >= 0).

    Returns
    -------
    result : ndarray, same shape as x
        gamma(a, x) at each point.
    """
    assert a >= 1, f"a must be >= 1, got {a}"
    fact_a_minus_1 = float(factorial(a - 1, exact=True))

    # Compute partial sum: sum_{m=0}^{a-1} x^m / m!
    partial_sum = np.ones_like(x, dtype=float)  # m=0 term
    term = np.ones_like(x, dtype=float)
    for m in range(1, a):
        term = term * x / m
        partial_sum += term

    exp_neg_x = np.exp(-x)
    result = fact_a_minus_1 * (1.0 - exp_neg_x * partial_sum)
    result = np.where(x <= 0, 0.0, result)
    return result


class LaguerreZeffExpansion:
    """Precomputed Laguerre expansion for fast algebraic Z_eff evaluation.

    After fitting, stores a polynomial P(X) and constant C_inf such that:
        N_core(r) = C_inf - exp(-X) * P(X)
    where X = 2*beta*r.

    This representation requires only one exponential evaluation plus a
    polynomial evaluation per call, making it O(n_basis) per r-point
    with small constant factor (no loops over incomplete gamma functions).
    """

    def __init__(self, coeffs: np.ndarray, beta: float, n_electrons: float = 2.0) -> None:
        self.coeffs = coeffs
        self.beta = beta
        self.n_electrons = n_electrons
        self.n_basis = len(coeffs)

        # Precompute the polynomial-exponential representation:
        # N_core(r) = C_inf - exp(-X) * P(X) where X = 2*beta*r
        #
        # From the Laguerre expansion:
        # N_core = sum_k c_k * (1/2) * sum_j a_{k,j} * gamma(j+3, X)
        # gamma(j+3, X) = (j+2)! * [1 - e^{-X} * sum_{m=0}^{j+2} X^m/m!]
        #
        # So: N_core = sum_j d_j * (j+2)! - e^{-X} * sum_j d_j * (j+2)! * sum_m X^m/m!
        # where d_j = (1/2) * sum_k c_k * a_{k,j}

        # Step 1: Compute combined power coefficients d_j
        k_max = self.n_basis - 1
        j_max = k_max  # Highest power from L_{k_max}
        d = np.zeros(j_max + 1)
        for k in range(self.n_basis):
            if abs(coeffs[k]) < 1e-30:
                continue
            lag_c = _laguerre_coefficients(k)
            d[:k + 1] += 0.5 * coeffs[k] * lag_c

        # Step 2: C_inf = sum_j d_j * (j+2)!
        self._C_inf = 0.0
        for j in range(j_max + 1):
            self._C_inf += d[j] * float(factorial(j + 2, exact=True))

        # Step 3: Polynomial P(X) = sum_m p_m * X^m
        # where p_m = sum_{j: j+2 >= m} d_j * (j+2)! / m!
        # The maximum power in P is j_max + 2
        m_max = j_max + 2
        self._poly_coeffs = np.zeros(m_max + 1)
        for j in range(j_max + 1):
            if abs(d[j]) < 1e-30:
                continue
            f_j2 = float(factorial(j + 2, exact=True))
            for m in range(j + 3):  # m = 0, 1, ..., j+2
                f_m = float(factorial(m, exact=True))
                self._poly_coeffs[m] += d[j] * f_j2 / f_m

    def n_core(self, r: np.ndarray) -> np.ndarray:
        """Evaluate N_core(r) using fast polynomial form with stability guard.

        Primary path: C_inf - exp(-X)*P(X) via Horner's method (fast).
        Fallback: per-term incomplete gamma summation when polynomial form
        produces out-of-range values (catastrophic cancellation detection).

        Parameters
        ----------
        r : ndarray
            Radial distances (bohr).

        Returns
        -------
        N_core : ndarray
            Cumulative core electron count, clamped to [0, n_electrons].
        """
        r_arr = np.atleast_1d(np.asarray(r, dtype=float))
        X = 2.0 * self.beta * r_arr

        # Fast path: polynomial-exponential form
        poly = np.zeros_like(X)
        for m in range(len(self._poly_coeffs) - 1, -1, -1):
            poly = poly * X + self._poly_coeffs[m]
        result = self._C_inf - np.exp(-X) * poly
        result = np.where(r_arr <= 0, 0.0, result)

        # Stability guard: detect catastrophic cancellation
        # Three failure modes:
        # 1. Value outside [-0.5, n_electrons+0.5] (gross error)
        # 2. Value is NaN
        # 3. Non-monotonicity: N_core should be monotonically non-decreasing;
        #    any decrease indicates catastrophic cancellation in the
        #    polynomial-exponential subtraction. This catches the Be (Z=4)
        #    case where C_inf ≈ 2 and exp(-X)*P(X) ≈ 2 at intermediate X.
        bad_mask = (result < -0.5) | (result > self.n_electrons + 0.5) | np.isnan(result)

        # Non-monotonicity detection: check for decreasing N_core values
        # along increasing r (only meaningful when r is sorted)
        sort_idx = np.argsort(r_arr)
        sorted_result = result[sort_idx]
        running_max = np.maximum.accumulate(sorted_result)
        # A point is bad if it drops below the running max by more than
        # a small tolerance (allow for numerical noise)
        mono_bad_sorted = sorted_result < (running_max - 0.01)
        # Map back to original order
        mono_bad = np.zeros_like(bad_mask)
        mono_bad[sort_idx] = mono_bad_sorted
        bad_mask = bad_mask | mono_bad
        if np.any(bad_mask):
            bad_idx = np.where(bad_mask)[0]
            X_bad = X[bad_idx]
            stable = np.zeros(len(bad_idx), dtype=float)
            for k in range(self.n_basis):
                if abs(self.coeffs[k]) < 1e-30:
                    continue
                lag_c = _laguerre_coefficients(k)
                for j in range(k + 1):
                    if abs(lag_c[j]) < 1e-30:
                        continue
                    a = j + 3
                    weight = 0.5 * self.coeffs[k] * lag_c[j]
                    stable += weight * _lower_incomplete_gamma_int(a, X_bad)
            result[bad_idx] = stable

        return np.clip(result, 0.0, self.n_electrons)

    def z_eff(self, r: np.ndarray, Z: float) -> np.ndarray:
        """Evaluate Z_eff(r) = Z - N_core(r).

        Parameters
        ----------
        r : ndarray
            Radial distances (bohr).
        Z : float
            Nuclear charge.

        Returns
        -------
        z_eff : ndarray
            Effective nuclear charge.
        """
        r_arr = np.atleast_1d(np.asarray(r, dtype=float))
        result = Z - self.n_core(r_arr)
        result = np.where(r_arr <= 0, Z, result)
        return result

    def density(self, r: np.ndarray) -> np.ndarray:
        """Evaluate n(r) from the Laguerre expansion.

        Parameters
        ----------
        r : ndarray
            Radial distances (bohr).

        Returns
        -------
        n_r : ndarray
            Radial number density.
        """
        r_arr = np.atleast_1d(np.asarray(r, dtype=float))
        x = 2.0 * self.beta * r_arr
        envelope = r_arr**2 * np.exp(-2.0 * self.beta * r_arr) * (4.0 * self.beta**3)

        # Evaluate Laguerre sum vectorized
        laguerre_sum = np.zeros_like(r_arr, dtype=float)
        for k in range(self.n_basis):
            if abs(self.coeffs[k]) < 1e-30:
                continue
            laguerre_sum += self.coeffs[k] * eval_laguerre(k, x)

        return np.maximum(envelope * laguerre_sum, 0.0)


def fit_density_laguerre(
    r_grid: np.ndarray,
    n_r: np.ndarray,
    n_basis: int = 30,
    beta: Optional[float] = None,
) -> Tuple[np.ndarray, float]:
    """Project grid density onto Laguerre spectral basis.

    Given n(r) on a grid, compute expansion coefficients c_k such that:

        n(r) approx= (4*beta^3) * sum_k c_k * r^2 * L_k(2*beta*r) * exp(-2*beta*r)

    Parameters
    ----------
    r_grid : ndarray of shape (N,)
        Radial grid points (bohr).
    n_r : ndarray of shape (N,)
        Radial number density on the grid.
    n_basis : int
        Number of Laguerre basis functions.
    beta : float, optional
        Exponential decay parameter. If None, estimated from the density.

    Returns
    -------
    coeffs : ndarray of shape (n_basis,)
        Expansion coefficients c_k.
    beta : float
        The decay parameter used.
    """
    if beta is None:
        total = np.trapezoid(n_r, r_grid)
        if total > 0:
            r_mean = np.trapezoid(r_grid * n_r, r_grid) / total
            beta = 1.5 / r_mean if r_mean > 0.01 else 2.0
        else:
            beta = 2.0

    x = 2.0 * beta * r_grid
    envelope = r_grid**2 * np.exp(-2.0 * beta * r_grid) * (4.0 * beta**3)

    # Build design matrix: A[i, k] = envelope_i * L_k(x_i)
    A = np.zeros((len(r_grid), n_basis))
    for k in range(n_basis):
        A[:, k] = envelope * eval_laguerre(k, x)

    # Ridge regression for stability
    AtA = A.T @ A
    Atb = A.T @ n_r
    ridge = 1e-12 * np.trace(AtA) / n_basis
    coeffs = np.linalg.solve(AtA + ridge * np.eye(n_basis), Atb)

    return coeffs, beta


def n_core_algebraic(
    r: np.ndarray,
    coeffs: np.ndarray,
    beta: float,
    n_electrons: float = 2.0,
) -> np.ndarray:
    """Evaluate N_core(r) algebraically (convenience wrapper).

    Creates a LaguerreZeffExpansion and evaluates. For repeated calls,
    use LaguerreZeffExpansion directly.
    """
    expansion = LaguerreZeffExpansion(coeffs, beta, n_electrons)
    return expansion.n_core(r)


def z_eff_algebraic(
    r: np.ndarray,
    Z: float,
    coeffs: np.ndarray,
    beta: float,
    n_electrons: float = 2.0,
) -> np.ndarray:
    """Evaluate Z_eff(r) algebraically (convenience wrapper)."""
    expansion = LaguerreZeffExpansion(coeffs, beta, n_electrons)
    return expansion.z_eff(r, Z)


def density_algebraic(
    r: np.ndarray,
    coeffs: np.ndarray,
    beta: float,
) -> np.ndarray:
    """Evaluate n(r) from the Laguerre expansion (convenience wrapper)."""
    expansion = LaguerreZeffExpansion(coeffs, beta)
    return expansion.density(r)
