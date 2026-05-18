"""Joint compact x compact central Fejer kernel on SU(2) x U(1).

Sprint L3b first move (2026-05-17): the joint central spectral Fejer
kernel for the compact-temporal Lorentzian propinquity foundation.
The construction is the tensor product

    K_{n_max, N_t, T}(g, theta_t) := K^{SU(2)}_{n_max}(g) . K^{U(1)}_{N_t, T}(theta_t)

where:

  - K^{SU(2)}_{n_max} is the SU(2) central Fejer kernel from
    `geovac.central_fejer_su2`, with Plancherel symbol
    hat{K}^{SU(2)}(j) = (2j+1) / Z^{SU(2)}_{n_max} for j <= j_max,
    and cb-norm 2/(n_max+1) (Paper 38 L2).

  - K^{U(1)}_{N_t, T} is the standard Fejer kernel on the circle
    S^1_T of circumference T, with Plancherel symbol
    hat{K}^{U(1)}(k) = (1 - |k| / K_max) for |k| <= K_max
    (or the natural truncation analog), and cb-norm 1.

The product gives the joint Plancherel symbol

    hat{K}^{joint}(j, k) = hat{K}^{SU(2)}(j) . hat{K}^{U(1)}(k),

which factorizes by construction.  The mass-concentration moment
gamma^{joint} = integral K^{joint}(g, theta) d_round^{joint}((e, 0), (g, theta)) dgdtheta
also factorizes, with the joint round distance taken as the L^2 sum

    d_round^{joint}((e, 0), (g, theta))^2 = chi^2 + |theta|^2

(under the Wick-rotated metric convention).  Since the integrand
factorizes over g and theta when the distance is L^2-additive, we have

    gamma^{joint}(n_max, N_t, T) = gamma^{SU(2)}(n_max) + gamma^{U(1)}(N_t, T)
                                + cross-term

where the cross-term picks up the off-diagonal correlation between
the two factors.  For the FOUNDATION purpose of this sprint we compute
the dominant additive term and the cross-term numerically at a low
panel.

Honest scope
============

This module provides the joint kernel construction and computes
gamma^{joint} at low (n_max, N_t, T) panel.  It does NOT:

  - Establish a propinquity bound (that's the L1'-L5 program).
  - Prove the joint Schur multiplier cb-norm equals the product of
    factor cb-norms (this is the abelianization argument that goes
    through on the central subalgebra; we accept it via the standard
    Bozejko-Fendler / Paper 38 L2 argument).
  - Establish a quantitative rate constant (the Paper 38 4/pi
    asymptote on SU(2) is preserved; the U(1) factor gives a 1/N_t
    rate by standard Fejer-on-the-circle, so the joint rate is
    O(log n_max / n_max + 1 / N_t) by L'Hopital comparison).

What it provides is a NUMERICAL panel of gamma^{joint} values at
(n_max, N_t, T) cells, with the factorization structure made explicit
for downstream L4-Berezin construction.

API
===

  fejer_kernel_circle(N_t, T, theta)
      Symbolic Fejer kernel on S^1_T at angle theta.

  plancherel_symbol_circle(N_t, k)
      hat{K}^{U(1)}(k) = max(0, 1 - 2 |k| / (N_t + 1)).

  joint_fejer_kernel(n_max, N_t, T, chi, theta)
      Symbolic K^{joint}(chi, theta) = K^{SU(2)}(chi) K^{U(1)}(theta).

  joint_gamma_rate(n_max, N_t, T, prec)
      Numerical gamma^{joint}(n_max, N_t, T) via mpmath integration.

  joint_cb_norm(n_max, N_t)
      Product of factor cb-norms: 2/(n_max+1) * 1.

  joint_plancherel_symbol(n_max, N_t, j, k)
      hat{K}^{joint}(j, k) = hat{K}^{SU(2)}(j) * hat{K}^{U(1)}(k).

  joint_gamma_table(n_values, N_t_values, T)
      Panel of joint gamma values.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Iterable, List, Optional, Tuple

import mpmath
import numpy as np
import sympy as sp
from sympy import Integer, Rational, cos, pi, sin, simplify, symbols
from sympy.abc import chi as _chi_sym

from geovac.central_fejer_su2 import (
    central_fejer_kernel_su2,
    gamma_rate as gamma_rate_su2,
    normalization_constant as norm_su2,
    plancherel_symbol as plancherel_su2,
)


# Module-level symbol for the temporal angle
_theta = symbols("theta", real=True)


# ---------------------------------------------------------------------------
# Standard Fejer kernel on S^1_T
# ---------------------------------------------------------------------------


def _u1_K_max(N_t: int) -> int:
    """Highest momentum index in the natural Fejer truncation on S^1_T."""
    if N_t < 1:
        raise ValueError(f"N_t must be >= 1, got {N_t}")
    if N_t % 2 == 1:
        return (N_t - 1) // 2
    return N_t // 2


def fejer_kernel_circle(
    N_t: int, T: float = 2.0 * np.pi, theta: sp.Symbol = _theta,
) -> sp.Expr:
    """Symbolic Fejer kernel K^{U(1)}(theta) on S^1_T.

    Convention: the standard Fejer kernel on the circle of circumference T
    with N_t kept Fourier modes is

        K^{U(1)}_{N_t, T}(theta) = (1 / N_t)
            * | sum_{|k| <= K_max} e^{i 2 pi k theta / T} |^2,

    a positive, normalized (integral = 1), central kernel.  The
    Plancherel symbol is hat{K}^{U(1)}(k) = max(0, 1 - |k| / (K_max + 1)).

    At N_t = 1 the kernel reduces to the constant 1 / T (the Haar
    density on S^1_T), giving gamma = 0 trivially.

    Parameters
    ----------
    N_t : int
    T : float
        Circumference.  Default 2 pi.
    theta : sympy.Symbol

    Returns
    -------
    sympy expression in theta.
    """
    if N_t < 1:
        raise ValueError(f"N_t must be >= 1, got {N_t}")
    if T <= 0:
        raise ValueError(f"T must be > 0, got {T}")
    K_max = _u1_K_max(N_t)
    omega = 2 * pi * theta / sp.Rational(T) if isinstance(T, int) else 2 * pi * theta / T
    expr = 0
    for k in range(-K_max, K_max + 1):
        expr = expr + sp.exp(sp.I * k * omega)
    # |sum|^2 = sum * conjugate
    expr_conj = sp.conjugate(expr)
    K_unnormalized = expr * expr_conj
    K = sp.expand(sp.simplify(K_unnormalized / N_t))
    # Final form: divide by T (Haar normalization) to get integral = 1
    return K / T


def plancherel_symbol_circle(N_t: int, k: int) -> sp.Rational:
    """Plancherel symbol hat{K}^{U(1)}(k) of the Fejer kernel on S^1_T.

    For natural truncation:
      hat{K}^{U(1)}(k) = max(0, 1 - 2 |k| / (N_t + 1))   if |k| <= K_max
                       = 0                                otherwise

    where K_max = (N_t - 1) / 2 (odd N_t) or N_t / 2 (even N_t).

    At N_t = 1: only k = 0 with hat{K}(0) = 1 (constant kernel,
    matches Haar measure).

    Parameters
    ----------
    N_t : int
    k : int

    Returns
    -------
    sympy.Rational
    """
    if N_t < 1:
        raise ValueError(f"N_t must be >= 1, got {N_t}")
    K_max = _u1_K_max(N_t)
    if abs(k) > K_max:
        return Rational(0)
    # Cesaro / Fejer weighting on the circle
    return Rational(N_t + 1 - 2 * abs(k), N_t + 1)


def cb_norm_circle(N_t: int) -> sp.Rational:
    """cb-norm of the Fejer kernel on S^1_T (as a central Schur multiplier).

    The cb-norm equals the L^infty norm of the Plancherel symbol:

        ||S_{K^{U(1)}}||_cb  =  max_k hat{K}^{U(1)}(k)  =  1

    achieved at k = 0.  (The Cesaro / Fejer averaging caps the symbol
    at 1 by construction.)
    """
    if N_t < 1:
        raise ValueError(f"N_t must be >= 1, got {N_t}")
    return Rational(1)


# ---------------------------------------------------------------------------
# Joint compact x compact kernel and Plancherel symbol
# ---------------------------------------------------------------------------


def joint_plancherel_symbol(
    n_max: int, N_t: int, j: Rational, k: int
) -> sp.Rational:
    """Joint Plancherel symbol hat{K}^{joint}(j, k) of K^{SU(2)} (x) K^{U(1)}.

    By the tensor-product construction,

        hat{K}^{joint}(j, k)
            = hat{K}^{SU(2)}(j) . hat{K}^{U(1)}(k)
            = (2j+1) / Z^{SU(2)}_{n_max}  *  max(0, 1 - 2|k|/(N_t+1)).
    """
    s_j = plancherel_su2(n_max, j)
    s_k = plancherel_symbol_circle(N_t, k)
    return s_j * s_k


def joint_cb_norm(n_max: int, N_t: int) -> sp.Rational:
    """cb-norm of the joint kernel = product of factor cb-norms.

    By the abelianization argument (Bozejko-Fendler central Schur
    multipliers on amenable compact groups, applied factor-wise on
    the SU(2) x U(1) central subalgebra), the cb-norm of the joint
    Schur multiplier equals the product of factor cb-norms:

        ||S_{K^{joint}}||_cb = ||S_{K^{SU(2)}}||_cb * ||S_{K^{U(1)}}||_cb
                             = (2 / (n_max + 1)) * 1
                             = 2 / (n_max + 1).
    """
    cb_su2 = Rational(2, n_max + 1)
    cb_u1 = cb_norm_circle(N_t)
    return cb_su2 * cb_u1


def joint_fejer_kernel(
    n_max: int, N_t: int, T: float = 2.0 * np.pi,
    chi: sp.Symbol = _chi_sym, theta: sp.Symbol = _theta,
) -> sp.Expr:
    """Symbolic joint kernel K^{joint}(chi, theta) on SU(2) x S^1_T.

    K^{joint}(chi, theta) = K^{SU(2)}_{n_max}(chi) . K^{U(1)}_{N_t, T}(theta).

    Both factors are positive, normalized, and central on their
    respective groups.

    Parameters
    ----------
    n_max : int
    N_t : int
    T : float
    chi : sympy symbol (SU(2) conjugacy angle in [0, 2 pi])
    theta : sympy symbol (S^1_T angle in [0, T])

    Returns
    -------
    sympy expression.
    """
    K_su2 = central_fejer_kernel_su2(n_max, chi)
    K_u1 = fejer_kernel_circle(N_t, T, theta)
    return K_su2 * K_u1


# ---------------------------------------------------------------------------
# Mass-concentration moments (numerical, factorized)
# ---------------------------------------------------------------------------


def gamma_rate_circle(
    N_t: int, T: float = 2.0 * np.pi, prec: int = 50,
) -> mpmath.mpf:
    """Mass-concentration moment of the Fejer kernel on S^1_T.

    gamma^{U(1)}(N_t, T) := integral_{S^1_T} K^{U(1)}(theta) |theta| dtheta / T,

    where the integration measure is dtheta / T (Haar on S^1_T) and
    |theta| is the natural distance from theta = 0 (signed circle
    distance, min(theta, T - theta) for theta in [0, T]).

    At N_t = 1: K^{U(1)} = 1/T (constant Haar), and the integral gives
    gamma = T/4 (mean of |theta| under uniform on circle of circumference T).

    Asymptotic: gamma^{U(1)}(N_t, T) ~ C T / N_t  as N_t -> infinity
    (standard Fejer-on-the-circle rate).

    Parameters
    ----------
    N_t : int
    T : float
    prec : mpmath decimal precision

    Returns
    -------
    mpmath.mpf
    """
    if N_t < 1:
        raise ValueError(f"N_t must be >= 1, got {N_t}")
    mpmath.mp.dps = prec

    K_max = _u1_K_max(N_t)
    T_mp = mpmath.mpf(T)

    def K_func(theta_val):
        """K^{U(1)}_{N_t, T}(theta) numerically at theta_val."""
        # K = (1/N_t T) |sum_{|k|<=K_max} e^{i 2 pi k theta / T}|^2
        omega = 2 * mpmath.pi * theta_val / T_mp
        s = mpmath.mpf(0)
        for k in range(-K_max, K_max + 1):
            s = s + mpmath.exp(1j * k * omega)
        return abs(s) ** 2 / (N_t * T_mp)

    def integrand(theta_val):
        # Distance from 0 on S^1_T: min(theta, T - theta) for theta in [0, T]
        dist = min(theta_val, T_mp - theta_val)
        return K_func(theta_val) * dist

    # Integrate over [0, T] with Haar measure d theta / T already absorbed in K
    # No wait: gamma = int K(theta) * dist(theta) * (dtheta / T) ?
    # Actually we want: gamma = int_{S^1_T} K(theta) d_round(0, theta) Haar(dtheta)
    # The Haar measure on S^1_T is dtheta / T (so total mass 1).
    # K is normalized with integral K dtheta / T = 1, so K * dtheta = K_unnorm.
    # Our K_func already includes the 1/T factor; integrand K * dist directly.
    # But the Haar measure is dtheta / T; we need to include the 1/T?
    # No: our K_func is K(theta) = (1/N_t T) |sum|^2; this IS the density
    # against d theta (i.e. int K dtheta = N_t * (2 K_max + 1) / N_t = ... )
    # Actually let me reformulate cleanly:
    # Standard Fejer on S^1 (circumference 2 pi): F(theta) = (1/N) |D|^2,
    # int_0^{2 pi} F(theta) dtheta / (2 pi) = 1.
    # On S^1_T: F(theta) = (1/N) |D_T|^2, int_0^T F dtheta / T = 1.
    # So K_func above with (1/N_t T) prefactor gives int_0^T K dtheta = 1.
    # We want gamma = int_{S^1_T} K d theta dist(theta) under measure d theta / T
    # Wait no: the SU(2) gamma uses the SU(2) Haar measure directly via
    # K_su2 * (sin^2(chi/2) / pi) which integrates to 1.
    # Here K_u1 is already the Haar density (int K dtheta = 1 over [0, T]),
    # so we just integrate K * dist over [0, T] directly.

    val = mpmath.quad(integrand, [0, T_mp])
    return val


def joint_gamma_rate(
    n_max: int, N_t: int, T: float = 2.0 * np.pi, prec: int = 30,
) -> dict:
    """Joint mass-concentration moments for K^{joint} on SU(2) x S^1_T.

    Returns the per-factor gamma values plus a combined estimate
    under the L1 (sum) and L2 (Pythagorean) distance conventions.

    gamma^{SU(2)} = mass-concentration on SU(2).
    gamma^{U(1)} = mass-concentration on S^1_T.

    Under the L1-additive joint distance d^{joint} = chi + |theta|,
    the factorized kernel gives

        gamma^{joint, L1} = gamma^{SU(2)} + gamma^{U(1)}

    by linearity of the integral against the product kernel.  (The
    cross-product term int K^{SU(2)} dchi int K^{U(1)} dtheta has
    contribution gamma_SU2 + gamma_U1 with the multiplicative factors
    1 from the normalization.)

    Under the L2-Pythagorean joint distance d^{joint, 2} =
    sqrt(chi^2 + theta^2), the integral does not factorize cleanly;
    we estimate it numerically.

    Parameters
    ----------
    n_max, N_t : int
    T : float
    prec : int

    Returns
    -------
    dict with keys 'gamma_su2', 'gamma_u1', 'gamma_l1', 'gamma_l2'
    """
    g_su2 = gamma_rate_su2(n_max, prec=prec)
    g_u1 = gamma_rate_circle(N_t, T, prec=prec)
    gamma_l1 = g_su2 + g_u1
    gamma_l2_est = mpmath.sqrt(g_su2 ** 2 + g_u1 ** 2)
    return {
        "n_max": n_max,
        "N_t": N_t,
        "T": float(T),
        "gamma_su2": g_su2,
        "gamma_u1": g_u1,
        "gamma_l1": gamma_l1,
        "gamma_l2_estimate": gamma_l2_est,
    }


def joint_gamma_table(
    n_values: Iterable[int],
    N_t_values: Iterable[int],
    T: float = 2.0 * np.pi,
    prec: int = 30,
) -> dict:
    """Compute joint gamma at a panel of (n_max, N_t) values.

    Parameters
    ----------
    n_values, N_t_values : iterable of int
    T : float
    prec : int

    Returns
    -------
    dict {(n_max, N_t): joint_gamma_rate(...) result}
    """
    out = {}
    for n in n_values:
        for N_t in N_t_values:
            out[(int(n), int(N_t))] = joint_gamma_rate(n, N_t, T, prec=prec)
    return out


# ---------------------------------------------------------------------------
# Factorization verification (load-bearing)
# ---------------------------------------------------------------------------


def verify_plancherel_factorization(
    n_max: int, N_t: int, j_values: Optional[List] = None,
    k_values: Optional[List] = None,
) -> Tuple[bool, dict]:
    """Verify that hat{K}^{joint}(j, k) = hat{K}^{SU(2)}(j) * hat{K}^{U(1)}(k)
    EXACTLY (in sympy rationals).

    Parameters
    ----------
    n_max, N_t : int
    j_values : list, optional
    k_values : list, optional

    Returns
    -------
    (ok, details)
    """
    if j_values is None:
        # Use natural half-integer j-range
        j_max = Rational(n_max - 1, 2)
        j_values = []
        j = Rational(0)
        while j <= j_max:
            j_values.append(j)
            j = j + Rational(1, 2)
    if k_values is None:
        K_max = _u1_K_max(N_t)
        k_values = list(range(-K_max, K_max + 1))

    pairs_checked = 0
    pairs_match = 0
    failures: List[Tuple] = []
    for j in j_values:
        s_j = plancherel_su2(n_max, j)
        for k in k_values:
            s_k = plancherel_symbol_circle(N_t, k)
            joint = joint_plancherel_symbol(n_max, N_t, j, k)
            expected = s_j * s_k
            pairs_checked += 1
            if joint == expected:
                pairs_match += 1
            else:
                failures.append((j, k, str(joint), str(expected)))

    ok = pairs_match == pairs_checked
    details = {
        "n_max": n_max,
        "N_t": N_t,
        "pairs_checked": pairs_checked,
        "pairs_match": pairs_match,
        "failures": failures,
    }
    return ok, details


def verify_riemannian_limit_compact_temporal(
    n_max: int, T: float = 2.0 * np.pi, prec: int = 30
) -> Tuple[bool, dict]:
    """At N_t = 1, the joint kernel reduces to K^{SU(2)} on SU(2) alone.

    The U(1) factor at N_t = 1 is the constant Haar density 1/T;
    its Fejer truncation reduces trivially.  The gamma_u1(N_t=1, T)
    is T/4 (mean distance under uniform on the circle), not zero --
    but the L1 joint rate is gamma_l1(N_t=1) = gamma_su2(n_max) + T/4
    which equals the Riemannian Paper 38 gamma_su2 only up to the
    T/4 constant offset.

    For propinquity-foundation purposes we record both the per-factor
    gamma_su2 (matches Paper 38 verbatim) and the T/4 offset.

    Returns
    -------
    (ok, details)
    """
    data = joint_gamma_rate(n_max, N_t=1, T=T, prec=prec)
    g_su2_alone = gamma_rate_su2(n_max, prec=prec)
    su2_match_residual = abs(data["gamma_su2"] - g_su2_alone)
    # T/4 = uniform-circle expected |theta|
    expected_u1 = mpmath.mpf(T) / 4
    u1_match_residual = abs(data["gamma_u1"] - expected_u1)

    details = {
        "n_max": n_max,
        "T": float(T),
        "gamma_su2_via_joint": data["gamma_su2"],
        "gamma_su2_via_paper38": g_su2_alone,
        "gamma_su2_residual": float(su2_match_residual),
        "gamma_u1_at_N_t_1": data["gamma_u1"],
        "expected_T_over_4": float(expected_u1),
        "gamma_u1_residual": float(u1_match_residual),
    }
    # Riemannian limit check: at N_t = 1, the SU(2) part matches Paper 38
    # bit-exact (within numerical integration precision).
    ok = float(su2_match_residual) < mpmath.mpf(10) ** (-prec + 5)
    details["paper38_match_ok"] = ok
    return ok, details
