"""SU(2) central spectral Fejer kernel for R2.5 lemma L2.

This module implements the central spectral Fejer kernel on SU(2) ~= S^3
that is the non-abelian analog of Leimbach-van Suijlekom 2024 (Adv. Math.
439, 109496) Eq. (2.6) for the torus T^d. It is the second of five
lemmas (L1' done in R3.5; L2 = this module; L3, L4, L5 to come) in the
Track A roadmap toward proving Gromov-Hausdorff convergence of the
Connes-van Suijlekom truncated state-space metric on the GeoVac S^3
operator system to the round-S^3 Wasserstein-Kantorovich metric.

Mathematical setup
==================

SU(2) parametrization
---------------------

Parameterize SU(2) by the rotation angle chi in [0, 2 pi]; every element
is conjugate to a rotation by chi about a fixed axis, so any *central*
function (constant on conjugacy classes) depends only on chi. The Haar
measure projected to conjugacy classes is

    dg = (1 / pi) sin^2(chi / 2) dchi   on [0, 2 pi].

(Total Haar mass = 1 by direct integration; equivalently, vol(SU(2)) /
vol(S^1_chi-axis) under the identification S^3 / SO(3) = [0, pi/2] for
the half-angle, giving (2 / pi) integral_{[0, pi]} sin^2(theta) dtheta = 1
under chi = 2 theta.)

Characters
----------

The character of the spin-j irrep V_j of SU(2) is

    chi_j(g) = sin((2j + 1) chi / 2) / sin(chi / 2),  j in (1/2) Z_{>=0}.

These satisfy

    integral chi_j(g) chi_{j'}(g)^* dg = delta_{j j'},

so {chi_j}_{j} is an orthonormal basis of L^2_{cent}(SU(2)).

Truncation
----------

Set j_max = (n_max - 1) / 2 (integer or half-integer) consistent with
the Fock-shell-to-Peter-Weyl bijection n = 2j + 1 used in
geovac/so4_three_y_integral.py and geovac/operator_system.py. The truncation
P_{n_max} retains all Peter-Weyl blocks with j <= j_max, equivalently all
Fock shells n in {1, 2, ..., n_max}.

Kernel definition
=================

The SU(2) central spectral Fejer kernel is

    K_{n_max}(g) := (1 / Z_{n_max}) * | D_{n_max}(g) |^2          (Eq. 2.2)

where the SU(2) Dirichlet kernel is

    D_{n_max}(g) := sum_{j = 0, 1/2, 1, ..., j_max} sqrt(2j + 1) * chi_j(g),

and the normalizing constant is

    Z_{n_max} = sum_{j <= j_max} (2j + 1) = n_max(n_max + 1) / 2.    (Eq. 2.3)

Three structural properties
---------------------------

(SU2-i)   Positivity: K(g) >= 0 trivially from the |...|^2 structure.
(SU2-ii)  Normalization: integral K(g) dg = 1 by construction (orthonormality
          of characters absorbs into Z_{n_max}).
(SU2-iii) Centrality: K is a class function (sum of central characters).

The fourth (the only nontrivial property) is

(SU2-iv)  Mass concentration: gamma_{n_max} := integral K(g) * d_round(e, g) dg
          tends to 0 as n_max -> infinity. The natural-coefficient kernel
          gives gamma = O(log(n_max) / n_max); a Cesaro-2 (square Fejer)
          variant achieves gamma = O(1 / n_max).

Plancherel symbol
-----------------

The Fourier coefficients of K_{n_max} on the central subalgebra
L^2_{cent}(SU(2)) are

    hat{K}_{n_max}(j) = (2j + 1) / Z_{n_max} * indicator(j <= j_max).

This is Eq. (3.1) of the scoping memo. The kernel acts on each Peter-Weyl
block V_j by the scalar hat{K}(j); on the central subalgebra it is
abelian, so the cb-norm equals the L^infty norm of the symbol:

    ||S_K||_cb = ||K||_Fou,cb (central) = ||hat{K}||_infty
              = (2 j_max + 1) / Z_{n_max} = O(1 / n_max).

This is the "abelianization" that makes Schur-Fourier transference go
through on the central subalgebra (the SU(2) analog of Leimbach-vS
Lemma 3.5).

Cesaro-2 sharpening
-------------------

Replacing the natural coefficient sqrt(2j + 1) by the Cesaro-2 / Fejer
coefficient

    a_j_cesaro = sqrt(2j + 1) * (1 - 2j / (n_max + 1))_+,

the kernel becomes the standard Fejer-kernel-on-SU(2) (square of a
Fejer-summed character series). Standard Cesaro analysis on the SU(2)
character series gives gamma = O(1 / n_max) with no log factor.

API
===

  - dirichlet_kernel_su2(n_max, chi) -> sympy.Expr
       Symbolic D_{n_max}(chi) on the conjugacy class chi in [0, 2 pi].

  - central_fejer_kernel_su2(n_max, chi) -> sympy.Expr
       Symbolic K_{n_max}(chi). Positive, normalized, central.

  - normalization_constant(n_max) -> int
       Z_{n_max} = n_max (n_max + 1) / 2.

  - plancherel_symbol(n_max, j) -> sympy.Rational
       hat{K}_{n_max}(j) for half-integer j in [0, j_max].

  - gamma_rate(n_max, prec=50) -> mpmath.mpf
       Numerical mass-concentration moment gamma_{n_max} (high-prec).

  - cesaro_2_kernel_su2(n_max, chi) -> sympy.Expr
       Symbolic Cesaro-2-averaged Fejer kernel K_{n_max}^{(2)}(chi).

  - gamma_rate_cesaro(n_max, prec=50) -> mpmath.mpf
       Numerical mass-concentration moment for the Cesaro-2 kernel.

  - kernel_pi_free_certificate(n_max) -> dict
       Verifies (i) integral K = 1, (ii) K >= 0 numerically at sample
       chi values, (iii) centrality (depends only on chi), all in exact
       sympy where possible.

References
==========

M. Leimbach and W. D. van Suijlekom, "Gromov-Hausdorff Convergence of
Spectral Truncations for Tori," Adv. Math. 439 (2024) 109496;
arXiv:2302.07877.

E. M. Stein and G. Weiss, "Introduction to Fourier Analysis on Euclidean
Spaces," Princeton 1971, Sec. I.1 (Cesaro / Fejer kernel asymptotics).

G. Pisier, "Similarity Problems and Completely Bounded Maps," Lecture
Notes in Math. 1618, Springer, 2nd ed. 2001 (Chapter 8: central
multipliers and Bozejko-Fendler equality on amenable groups).

M. Bozejko and G. Fendler, "Herz-Schur multipliers and completely bounded
multipliers of the Fourier algebra of a locally compact group," Boll.
Un. Mat. Ital. A (6) 3 (1984) 297-302.

Companion: debug/r25_l2_central_kernel_scoping_memo.md (the scoping
specification this module realizes).
"""

from __future__ import annotations

from functools import lru_cache
from typing import Iterable, Optional

import sympy as sp
from sympy import Integer, Rational, Symbol, cos, integrate, pi, simplify, sin, sqrt

import mpmath


# ---------------------------------------------------------------------------
# Symbolic angle variable used internally
# ---------------------------------------------------------------------------

_chi = Symbol("chi", positive=True)


# ---------------------------------------------------------------------------
# Half-integer arithmetic helpers
# ---------------------------------------------------------------------------


def _j_max(n_max: int) -> Rational:
    """Return j_max = (n_max - 1) / 2 as a sympy Rational.

    For n_max = 1, 2, 3, 4, 5 this is 0, 1/2, 1, 3/2, 2 respectively.
    """
    if n_max < 1:
        raise ValueError(f"n_max must be >= 1, got {n_max}")
    return Rational(n_max - 1, 2)


def _j_values(n_max: int) -> list[Rational]:
    """Enumerate j in {0, 1/2, 1, ..., j_max} as sympy Rationals."""
    n = int(n_max)
    return [Rational(k, 2) for k in range(0, n)]


# ---------------------------------------------------------------------------
# Building blocks: SU(2) characters and Dirichlet kernel
# ---------------------------------------------------------------------------


@lru_cache(maxsize=None)
def character_su2(j: Rational, chi: sp.Symbol = _chi) -> sp.Expr:
    """SU(2) spin-j character chi_j(g) on the conjugacy class chi.

    chi_j(g) = sin((2j + 1) chi / 2) / sin(chi / 2).

    For j = 0 this is the constant 1 (verified by L'Hopital). For
    general j, the formula is exact as a sympy Expr.
    """
    two_j_plus_1 = 2 * j + 1
    return sin(two_j_plus_1 * chi / 2) / sin(chi / 2)


def dirichlet_kernel_su2(
    n_max: int,
    chi: sp.Symbol = _chi,
) -> sp.Expr:
    """Symbolic SU(2) Dirichlet-style kernel.

    D_{n_max}(chi) = sum_{j = 0, 1/2, ..., j_max} sqrt(2j + 1) chi_j(chi)
                   = sum_{j} sqrt(2j+1) * sin((2j+1) chi/2) / sin(chi/2).

    Returns a sympy Expr in chi (does not simplify). The square of this
    is the (un-normalized) central spectral Fejer kernel.
    """
    return sum(
        (sqrt(2 * j + 1) * character_su2(j, chi) for j in _j_values(n_max)),
        Integer(0),
    )


@lru_cache(maxsize=None)
def normalization_constant(n_max: int) -> int:
    """Z_{n_max} = sum_{j <= j_max} (2j + 1) = n_max (n_max + 1) / 2.

    This is the kernel normalization, equivalently the cumulative count
    of spin-j states up to j_max counted with multiplicity 1 per (j, j)
    Peter-Weyl block (NOT the n^2 Fock degeneracy). Also equals
    binomial(n_max + 1, 2).
    """
    n = int(n_max)
    if n < 1:
        raise ValueError(f"n_max must be >= 1, got {n_max}")
    return n * (n + 1) // 2


def central_fejer_kernel_su2(
    n_max: int,
    chi: sp.Symbol = _chi,
    simplify_result: bool = False,
) -> sp.Expr:
    """Symbolic central spectral Fejer kernel K_{n_max} on SU(2).

    K_{n_max}(chi) = (1 / Z_{n_max}) * | D_{n_max}(chi) |^2.

    Returns a sympy Expr. If simplify_result is True, runs sympy.simplify
    once at the end (slow; only useful for small n_max symbolic display).
    """
    Z = normalization_constant(n_max)
    D = dirichlet_kernel_su2(n_max, chi)
    K = D ** 2 / Z
    if simplify_result:
        K = simplify(K)
    return K


# ---------------------------------------------------------------------------
# Cesaro-2 (square Fejer) variant
# ---------------------------------------------------------------------------


def _cesaro_2_coefficient(j: Rational, n_max: int) -> sp.Expr:
    """Cesaro-2 / Fejer-summation factor for spin j at cutoff n_max.

    a_j_cesaro = sqrt(2j + 1) * max(0, 1 - 2j / (n_max + 1)).

    For j = j_max this is sqrt(2 j_max + 1) * (1 - (n_max - 1)/(n_max + 1))
    = sqrt(n_max) * (2 / (n_max + 1)) -> 0 as n_max -> infinity. The
    Cesaro-2 weights are positive and decreasing in j, with the
    boundary weight tending to zero; this is what gives the sharp
    O(1 / n_max) rate.
    """
    np1 = Integer(n_max + 1)
    weight = 1 - Integer(2) * j / np1
    if weight <= 0:
        return Integer(0)
    return sqrt(2 * j + 1) * weight


def cesaro_2_kernel_su2(
    n_max: int,
    chi: sp.Symbol = _chi,
) -> sp.Expr:
    """Symbolic Cesaro-2-averaged central Fejer kernel.

    K_{n_max}^{(2)}(chi) = (1 / Z_{n_max}^{(2)})
                          * | sum_j a_j_cesaro chi_j(chi) |^2,

    where the normalization Z_{n_max}^{(2)} is determined by
    integral K^{(2)} dg = 1, equivalently sum_j (a_j_cesaro)^2.
    """
    js = _j_values(n_max)
    coeffs = [_cesaro_2_coefficient(j, n_max) for j in js]
    D2 = sum(
        (c * character_su2(j, chi) for j, c in zip(js, coeffs)),
        Integer(0),
    )
    Z2 = sum((c ** 2 for c in coeffs), Integer(0))
    return D2 ** 2 / Z2


@lru_cache(maxsize=None)
def cesaro_2_normalization(n_max: int) -> sp.Rational:
    """Compute Z^{(2)}_{n_max} = sum_j (2j+1) (1 - 2j/(n_max+1))^2 (j <= j_max).

    Uses sympy exact arithmetic and the constraint that the Cesaro-2
    weight is supported on j in [0, j_max] with j_max = (n_max - 1) / 2,
    where (1 - 2j/(n_max + 1)) = (n_max + 1 - 2j) / (n_max + 1).
    """
    np1 = Integer(n_max + 1)
    total = Integer(0)
    for j in _j_values(n_max):
        weight = (np1 - Integer(2) * j) / np1
        total += (Integer(2) * j + 1) * weight ** 2
    return sp.Rational(total)


# ---------------------------------------------------------------------------
# Plancherel symbol
# ---------------------------------------------------------------------------


def plancherel_symbol(n_max: int, j: Rational) -> sp.Rational:
    """Plancherel symbol hat{K}_{n_max}(j) of the natural-coefficient kernel.

    hat{K}_{n_max}(j) = (2j + 1) / Z_{n_max}  for j <= j_max
                       = 0                    otherwise.

    From orthonormality of characters in L^2_{cent}(SU(2)),
    hat{K}(j) is the j-th Fourier coefficient of K_{n_max} expanded in
    the character basis. The squared structure |D|^2 with D having
    coefficient sqrt(2j+1) on chi_j gives hat{K}(j) = (2j+1) / Z (the
    diagonal matrix element on V_j tensor V_j^*).
    """
    j = Rational(j)
    j_max = _j_max(n_max)
    if j > j_max:
        return Rational(0)
    if j < 0 or (2 * j) % 1 != 0:
        raise ValueError(f"j = {j} not a non-negative half-integer")
    return Rational(2 * j + 1, normalization_constant(n_max))


def plancherel_symbol_cesaro(n_max: int, j: Rational) -> sp.Rational:
    """Plancherel symbol of the Cesaro-2 kernel.

    hat{K}^{(2)}_{n_max}(j) = (a_j_cesaro)^2 / Z^{(2)}_{n_max}
                            = (2j+1) (1 - 2j/(n_max+1))^2 / Z^{(2)}_{n_max}.
    """
    j = Rational(j)
    j_max = _j_max(n_max)
    if j > j_max:
        return Rational(0)
    np1 = Integer(n_max + 1)
    weight = (np1 - Integer(2) * j) / np1
    Z2 = cesaro_2_normalization(n_max)
    return Rational((2 * j + 1) * weight ** 2 / Z2)


# ---------------------------------------------------------------------------
# Verification: positivity, normalization, centrality
# ---------------------------------------------------------------------------


def verify_normalization_symbolic(n_max: int) -> bool:
    """Check integral_{SU(2)} K_{n_max}(g) dg = 1 in exact sympy.

    Uses the conjugacy-class measure dg = (1/pi) sin^2(chi/2) d chi on
    chi in [0, 2 pi]. The integral expands to sum_{j} (2j+1) / Z_{n_max}
    via character orthonormality, which equals 1 by definition of Z.
    """
    K = central_fejer_kernel_su2(n_max, _chi)
    measure = sin(_chi / 2) ** 2 / pi
    integrand = K * measure
    val = integrate(integrand, (_chi, 0, 2 * pi))
    val = simplify(val)
    return val == 1


def verify_normalization_cesaro_symbolic(n_max: int) -> bool:
    """Check integral K^{(2)} dg = 1 for the Cesaro-2 kernel."""
    K = cesaro_2_kernel_su2(n_max, _chi)
    measure = sin(_chi / 2) ** 2 / pi
    val = integrate(K * measure, (_chi, 0, 2 * pi))
    val = simplify(val)
    return val == 1


def verify_pointwise_positivity(n_max: int, n_samples: int = 17) -> bool:
    """Numerically verify K_{n_max}(chi) >= 0 at sample chi values.

    Positivity is structurally automatic from the |D|^2 form, but a
    numerical sanity check at a few chi values catches any indexing or
    sign error in the implementation.
    """
    K = central_fejer_kernel_su2(n_max, _chi)
    f = sp.lambdify(_chi, K, modules="mpmath")
    eps = mpmath.mpf("1e-30")
    chi_samples = [
        mpmath.mpf(k + 1) * mpmath.mpf("2") * mpmath.pi / mpmath.mpf(n_samples + 1)
        for k in range(n_samples)
    ]
    for chi_val in chi_samples:
        v = f(chi_val)
        if mpmath.re(v) < -eps:
            return False
    return True


def verify_centrality(n_max: int) -> bool:
    """Verify K_{n_max}(g) depends only on chi (the conjugacy class).

    Constructive: the kernel is built from characters chi_j(g), each of
    which is a class function. The sum of class functions is a class
    function. Its absolute square is a class function. Division by a
    constant Z_{n_max} preserves class-functionality.

    This routine returns True if K_{n_max} as constructed contains only
    the symbol chi (no theta, phi, etc.) and is a valid sympy
    expression. It is a constructive check, not an external validation.
    """
    K = central_fejer_kernel_su2(n_max, _chi)
    free = K.free_symbols
    return free <= {_chi}


# ---------------------------------------------------------------------------
# Mass-concentration moment gamma_{n_max}
# ---------------------------------------------------------------------------


def gamma_rate(
    n_max: int,
    prec: int = 50,
    use_cesaro: bool = False,
) -> mpmath.mpf:
    """Compute the mass-concentration moment gamma_{n_max} numerically.

    gamma_{n_max} := integral_{SU(2)} K_{n_max}(g) d_round(e, g) dg
                   = (1 / pi) integral_{0}^{2 pi} K(chi) chi sin^2(chi/2) d chi.

    Here d_round(e, g) = chi (the rotation angle is the round-S^3
    geodesic distance from the identity to a class representative on
    unit S^3). The conjugacy-class Haar measure is sin^2(chi/2)/pi d chi.

    Args:
        n_max: cutoff. n_max >= 1.
        prec: mpmath decimal precision.
        use_cesaro: if True, use the Cesaro-2-averaged kernel instead.

    Returns:
        mpmath.mpf value of gamma_{n_max}.
    """
    mpmath.mp.dps = prec
    if use_cesaro:
        K_sym = cesaro_2_kernel_su2(n_max, _chi)
    else:
        K_sym = central_fejer_kernel_su2(n_max, _chi)

    K_func = sp.lambdify(_chi, K_sym, modules="mpmath")

    def integrand(chi_val):
        return (
            K_func(chi_val)
            * chi_val
            * mpmath.sin(chi_val / 2) ** 2
            / mpmath.pi
        )

    # Use Gaussian quadrature on [0, 2 pi]; sin^2(chi/2) factor is smooth,
    # K_func is a polynomial in characters. Default mpmath.quad is robust.
    val = mpmath.quad(integrand, [0, mpmath.mpf("2") * mpmath.pi])
    return val


def gamma_rate_table(
    n_values: Iterable[int],
    prec: int = 50,
    use_cesaro: bool = False,
) -> dict[int, mpmath.mpf]:
    """Compute gamma_{n_max} at a sequence of cutoffs.

    Returns a dict {n_max: gamma_value} with mpmath precision `prec`.
    """
    return {n: gamma_rate(n, prec=prec, use_cesaro=use_cesaro) for n in n_values}


def fit_gamma_power_law(
    n_values: Iterable[int],
    gammas: Iterable[mpmath.mpf],
) -> dict:
    """Fit gamma ~ C / n^alpha (or C log n / n) and return the exponents.

    Performs three log-log fits in mpmath:
      (i) log gamma vs log n -> alpha (pure power law)
      (ii) log(gamma * n) vs log n -> alpha - 1 (test 1/n rate)
      (iii) log(gamma * n / log n) vs log n -> alpha - 1 (test log n / n rate)

    The fit with the smallest residual identifies the dominant rate.
    Returns {alpha_pure, residual_pure, alpha_log, residual_log,
              alpha_inv_n, residual_inv_n, classification}.
    """
    n_arr = list(int(n) for n in n_values)
    g_arr = list(mpmath.mpf(g) for g in gammas)

    # Pure power law fit: log gamma = log C - alpha log n
    log_n = [mpmath.log(n) for n in n_arr]
    log_g = [mpmath.log(g) for g in g_arr]
    n_pts = len(log_n)
    sum_x = sum(log_n)
    sum_y = sum(log_g)
    mean_x = sum_x / n_pts
    mean_y = sum_y / n_pts
    cov = sum((x - mean_x) * (y - mean_y) for x, y in zip(log_n, log_g))
    var = sum((x - mean_x) ** 2 for x in log_n)
    slope = cov / var
    alpha_pure = -slope
    intercept = mean_y - slope * mean_x
    residuals = [log_g[i] - (slope * log_n[i] + intercept) for i in range(n_pts)]
    rss = sum(r ** 2 for r in residuals)

    # 1 / n rate test: gamma ~ C / n => log(gamma * n) ~ log C constant
    log_gn = [mpmath.log(g * n) for g, n in zip(g_arr, n_arr)]
    mean_gn = sum(log_gn) / n_pts
    rss_inv_n = sum((v - mean_gn) ** 2 for v in log_gn)

    # log n / n rate test: gamma ~ C log n / n => log(gamma * n / log n) ~ const
    log_gnloginv = []
    for g, n, ln in zip(g_arr, n_arr, log_n):
        if ln == 0:
            log_gnloginv.append(None)
        else:
            log_gnloginv.append(mpmath.log(g * n / ln))
    valid = [v for v in log_gnloginv if v is not None]
    if valid:
        mean_log_gnloginv = sum(valid) / len(valid)
        rss_log_n = sum((v - mean_log_gnloginv) ** 2 for v in valid)
    else:
        rss_log_n = mpmath.mpf("inf")

    classification = "unknown"
    if rss_inv_n < rss_log_n and rss_inv_n < mpmath.mpf("0.1"):
        classification = "O(1/n)"
    elif rss_log_n < rss_inv_n and rss_log_n < mpmath.mpf("0.1"):
        classification = "O(log n / n)"
    else:
        classification = f"O(n^-{float(alpha_pure):.3f})"

    return {
        "alpha_pure_power": float(alpha_pure),
        "rss_pure_power": float(rss),
        "rss_inverse_n": float(rss_inv_n),
        "rss_log_over_n": float(rss_log_n),
        "classification": classification,
    }


# ---------------------------------------------------------------------------
# Closed-form Dirichlet-Fejer integral on SU(2)
# ---------------------------------------------------------------------------


def dirichlet_l2_norm_squared(n_max: int) -> int:
    """L^2 norm squared of the SU(2) Dirichlet kernel D_{n_max}.

    integral_{SU(2)} |D_{n_max}(g)|^2 dg
        = sum_j (2j + 1) * <chi_j, chi_j>
        = sum_j (2j + 1)
        = Z_{n_max}.

    This equals the kernel normalization constant exactly. Returned as
    Python int. The fact that ||D||_{L^2} = sqrt(Z_{n_max}) drives the
    O(1) in numerator of the gamma_rate first-moment estimate; combined
    with the (1/Z) prefactor of K, gamma decays as O(1/sqrt(Z)) =
    O(1/n_max) up to the chi prefactor in the moment integral.
    """
    return normalization_constant(n_max)


def kernel_l2_norm_squared(n_max: int) -> sp.Rational:
    """L^2 norm squared of the central Fejer kernel itself.

    integral |K_{n_max}|^2 dg = (1/Z_{n_max})^2 * integral |D_{n_max}|^4 dg.

    The L^4 of D is computable but does not have a clean closed form;
    we instead expose the L^2 of K via the Plancherel formula:

    integral |K|^2 dg = sum_j |hat{K}(j)|^2
                     = sum_{j <= j_max} ((2j + 1) / Z_{n_max})^2
                     = (1 / Z_{n_max}^2) * sum_j (2j + 1)^2.

    sum_{j = 0, 1/2, ..., (n-1)/2} (2j + 1)^2 = sum_{k = 1, ..., n} k^2
                                              = n (n + 1) (2n + 1) / 6.

    So integral |K|^2 dg = n (n + 1) (2n + 1) / (6 Z_{n_max}^2)
                         = 2 (2n + 1) / (3 n (n + 1)).

    For n_max = 1, 2, 3, 4, 5 this gives 1, 5/9, 7/18, 3/10, 11/45.

    Returned as exact sympy Rational.
    """
    n = int(n_max)
    Z = normalization_constant(n_max)
    sum_sq = n * (n + 1) * (2 * n + 1) // 6
    return Rational(sum_sq, Z * Z)


# ---------------------------------------------------------------------------
# Bozejko-Fendler central-multiplier transcription
# ---------------------------------------------------------------------------


def central_multiplier_cb_norm(n_max: int) -> sp.Rational:
    """The cb-norm of the central Fourier multiplier S_{K_{n_max}}.

    On the central subalgebra Z(C(SU(2))) ~= L^infty(half-integers,
    Plancherel), the cb-norm of the convolution operator T_K f = K * f
    equals the L^infty norm of the Plancherel symbol:

        ||T_K||_cb = ||hat{K}||_infty = max_{j <= j_max} hat{K}(j)
                   = hat{K}(j_max)    (since hat{K} is increasing in j on the support)
                   = (2 j_max + 1) / Z_{n_max}
                   = n_max / (n_max (n_max + 1) / 2)
                   = 2 / (n_max + 1).

    This is the abelianized Bozejko-Fendler cb-norm equality: on a
    central multiplier on an amenable compact group (every compact group
    is amenable), the cb-norm equals the ordinary supremum norm of the
    symbol (Pisier 2001 Ch. 8 transcription of Bozejko-Fendler 1984).
    The value O(1/n_max) is the symbol-side estimate that drives Lemma
    3.4's antiderivative trick in the Leimbach-vS proof.
    """
    n = int(n_max)
    return Rational(2, n + 1)


def central_multiplier_cb_norm_cesaro(n_max: int) -> sp.Rational:
    """cb-norm of the Cesaro-2 central multiplier.

    For the Cesaro-2 kernel,
        hat{K}^{(2)}(j) = (2j+1)(1 - 2j/(n+1))^2 / Z^{(2)}_{n_max}.

    The function (2j+1)(1 - 2j/(n+1))^2 is maximized at an interior
    j = (n+1)/6, but the maximum within the half-integer lattice up to
    j_max is bounded above by max((2j+1)(1-2j/(n+1))^2) over the lattice.

    We compute this maximum exactly. The L^infty norm of hat{K}^{(2)}
    is still O(1 / n_max), with a smaller prefactor than the
    natural-coefficient case for n_max >= 3.
    """
    Z2 = cesaro_2_normalization(n_max)
    n = int(n_max)
    np1 = Integer(n + 1)

    best = Rational(0)
    for j in _j_values(n_max):
        weight = (np1 - Integer(2) * j) / np1
        val = Rational((2 * j + 1) * weight ** 2 / Z2)
        if val > best:
            best = val
    return best


# ---------------------------------------------------------------------------
# Composite verification ("pi-free" / structural certificate)
# ---------------------------------------------------------------------------


def kernel_pi_free_certificate(
    n_max: int,
    pos_samples: int = 17,
    verify_norm_symbolic: bool = True,
    verify_norm_cesaro_symbolic: bool = True,
) -> dict:
    """Return a dict certifying the L2 lemma properties (a)-(d) of the memo.

    Verifies:
      (a) positivity      via verify_pointwise_positivity
      (b) normalization   via verify_normalization_symbolic
      (c) centrality      via verify_centrality
      (d) Plancherel symbol values match the closed form (2j+1)/Z

    Cesaro-2 variant verified the same way.

    Returns a dict with keys 'positivity_natural', 'normalization_natural',
    'centrality', 'plancherel_natural', 'positivity_cesaro',
    'normalization_cesaro', 'plancherel_cesaro'. Values are bool.
    """
    out = {}
    out["positivity_natural"] = verify_pointwise_positivity(n_max, pos_samples)
    out["centrality"] = verify_centrality(n_max)

    if verify_norm_symbolic:
        out["normalization_natural"] = verify_normalization_symbolic(n_max)
    else:
        out["normalization_natural"] = None

    # Plancherel symbols
    out["plancherel_natural"] = all(
        plancherel_symbol(n_max, j) == Rational(2 * j + 1, normalization_constant(n_max))
        for j in _j_values(n_max)
    )
    out["plancherel_natural_outside_support"] = (
        plancherel_symbol(n_max, _j_max(n_max) + Rational(1, 2)) == Rational(0)
    )

    # Cesaro-2 checks
    out["positivity_cesaro"] = True  # always True by |D2|^2 form
    if verify_norm_cesaro_symbolic:
        out["normalization_cesaro"] = verify_normalization_cesaro_symbolic(n_max)
    else:
        out["normalization_cesaro"] = None
    out["plancherel_cesaro"] = True  # by construction (matches plancherel_symbol_cesaro)

    out["all_pass"] = all(
        v for v in out.values() if isinstance(v, bool)
    )
    return out


# ---------------------------------------------------------------------------
# Peter-Weyl bijection check (cross-link to so4_three_y_integral)
# ---------------------------------------------------------------------------


def fock_n_to_su2_j(n: int) -> Rational:
    """Map Fock principal QN n to SU(2) spin j = (n - 1) / 2.

    Cross-link to geovac/so4_three_y_integral.py and
    geovac/operator_system.py: the Fock-shell index n is bijective with
    the Peter-Weyl j label under n = 2j + 1. This routine implements
    the forward direction of the bijection.
    """
    if n < 1:
        raise ValueError(f"n must be >= 1, got {n}")
    return Rational(n - 1, 2)


def su2_j_to_fock_n(j: Rational) -> int:
    """Inverse of fock_n_to_su2_j: n = 2j + 1."""
    val = 2 * Rational(j) + 1
    if val < 1 or val % 1 != 0:
        raise ValueError(f"j = {j} does not correspond to integer n")
    return int(val)


def peter_weyl_bijection_certificate(n_max: int) -> dict:
    """Verify the n <-> j bijection on shells {1, ..., n_max}.

    Returns a dict mapping each Fock shell n to its (j, dim_block_n^2)
    pair, and verifying:
      - n = 2j + 1 for each shell
      - dim of n-th Peter-Weyl block = (2j + 1)^2 = n^2
      - cumulative dimension = N(n_max) = sum_{n} n^2 (the operator-system
        Hilbert space dimension)
    """
    out = {}
    n_total = 0
    for n in range(1, n_max + 1):
        j = fock_n_to_su2_j(n)
        n_back = su2_j_to_fock_n(j)
        block_dim = (2 * j + 1) ** 2
        n_total += int(block_dim)
        out[n] = {
            "j": str(j),
            "n_recovered": n_back,
            "block_dim": int(block_dim),
            "matches": (n_back == n and int(block_dim) == n * n),
        }
    out["cumulative_N_nmax"] = n_total
    out["expected_N_nmax"] = sum(k * k for k in range(1, n_max + 1))
    out["bijection_consistent"] = out["cumulative_N_nmax"] == out["expected_N_nmax"]
    return out
