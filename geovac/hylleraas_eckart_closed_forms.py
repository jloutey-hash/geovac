"""
Hylleraas-Eckart double-alpha master integral closed forms
==========================================================

Closed-form symbolic + cached numerical evaluator for the Hylleraas-Eckart
master integral

    I_HE^cosh(L, M, N; alpha, B)
        = integral_0^infty ds exp(-2 alpha s) s^L
          * integral_0^s du u^{N+1}
          * integral_{-u}^u dt t^{2M} (s^2 - t^2) cosh(B t)

This is the Hylleraas-Eckart extension of Bethe-Salpeter §32. The cosh(B t)
factor with B = beta_p +/- beta_q encodes the inner/outer screening
asymmetry of the Eckart 1933 trial function for He excited states.

CONVENTION (singlet vs triplet)
-------------------------------

The Eckart basis functions, in the (s, t, u) coordinate system with
s = r_1 + r_2, t = r_1 - r_2, u = r_12, are:

    phi_p^S(s, t, u) = exp(-alpha s) * cosh(beta_p t) * s^l_p * t^(2 m_p) * u^n_p  (singlet)
    phi_p^T(s, t, u) = exp(-alpha s) * sinh(beta_p t) * s^l_p * t^(2 m_p) * u^n_p  (triplet)

Cross-products that arise in matrix elements <phi_p | O | phi_q> reduce to
the SAME cosh master integral, just with different sign combinations:

    phi_p^S phi_q^S = exp(-2 alpha s) * (1/2)[cosh(B_+ t) + cosh(B_- t)] * s^L t^(2M) u^N
    phi_p^T phi_q^T = exp(-2 alpha s) * (1/2)[cosh(B_- t) - cosh(B_+ t)] * s^L t^(2M) u^N

where B_+ = beta_p + beta_q, B_- = beta_p - beta_q, L = l_p + l_q,
M = m_p + m_q, N = n_p + n_q. So one master integral I_HE^cosh suffices
for both sectors. A separate "sinh master integral" is NOT needed:
sinh(B t) * t^{2M} is odd in t and integrates to zero over [-u, u];
sinh(B t) * t^{2M+1} is even (used for kinetic-energy mixed terms only)
but the scoping plan (Track 2) handles kinetic T via quadrature, not
closed forms, so we never invoke that integral analytically.

The triplet basis-set convention (sinh * t^{2m}) means phi^T identically
vanishes at beta = 0 (sinh(0) = 0). The β -> 0 limit recovers the
single-alpha triplet (which uses t^{2m+1} with no sinh) by extracting an
overall β factor:

    sinh(beta t) * t^{2m} = beta * t^{2m+1} + O(beta^3)

so the variational triplet energy converges smoothly to the single-alpha
triplet energy in the beta -> 0 limit, just not bit-identically the
basis functions themselves. This is unlike the singlet sector where the
beta -> 0 reduction is bit-identical.

CONVERGENCE
-----------

The s-integration requires 2*alpha > |B|. For He 2^1S0 with alpha~1.5,
beta~0.5: B_+ = 1.0, B_- = 0.0, comfortably below 2*alpha = 3.0. For
He 2^3S1 with alpha~0.485, beta~1.485 (Bethe-Salpeter §32 Table 13),
single-beta cross-products give B_+ = 2*1.485 = 2.97 ~ 2*alpha = 0.97 *
3.06 = within tolerance only if alpha is bumped up; this is one of the
parameter-regime questions Track 2 must validate.

CLOSED FORM (for L = M = N = 0)
-------------------------------

    I_HE^cosh(0, 0, 0; alpha, B) = 64 / (4 alpha^2 - B^2)^3

B -> 0: 64 / (4 alpha^2)^3 = 1 / alpha^6, matching the existing
single-alpha master I(0,0,0; alpha) = 1/alpha^6 exactly.

For (L, M, N) > (0, 0, 0), the closed form is a rational function of
(alpha, B) with polynomial numerator and denominator of the form
(4 alpha^2 - B^2)^k * (smaller polynomial factors). All coefficients
are integers. By even-power symmetry (cosh is even in B, the integrand
is even in B since (s^2-t^2) and t^{2M} are even in t), the closed form
is a rational function of (alpha^2, B^2).

DISK CACHE
----------

Sympy symbolic evaluation is slow (~5-10 s per call at total degree
L+2M+N ~ 4, scaling roughly as (total_degree+5)^2.5). For production use,
this module precomputes closed forms once and caches integer-coefficient
polynomial representations to disk. Subsequent calls do polynomial
arithmetic in (alpha^2, B^2) at machine speed.

The cache is keyed by a schema version + (L, M, N). On schema change
(e.g., bug fix to closed-form derivation), bump SCHEMA_VERSION and the
cache rebuilds from scratch.

Architecture
------------

* I_HE_cosh_symbolic(L, M, N) -> sympy Expr  : SLOW symbolic closure.
* I_HE_cosh_polynomial(L, M, N) -> Polynomial : extract integer polynomial form.
* I_HE_cosh_numeric(L, M, N, alpha, B) -> float : fast evaluator (uses cache).
* precompute_table(L_max, M_max, N_max) -> Dict : one-shot precomputation.
* load_table() / save_table() : disk persistence.

References
----------
* Eckart, C. (1930). Phys. Rev. 36, 878. (Asymmetric variational trial.)
* Bethe, H. A. & Salpeter, E. E. (1957). Quantum Mechanics of One- and
  Two-Electron Atoms. Section 32. (Hylleraas master integral; Eckart
  trial functions and their He excited-state accuracies in Table 13.)
* Scoping memo: debug/hylleraas_eckart_scoping_memo.md (2026-05-18).

Author: GeoVac Track 1 Hylleraas-Eckart implementation.
"""

from __future__ import annotations

import json
import os
import pickle
import time
from dataclasses import dataclass
from fractions import Fraction
from typing import Dict, Optional, Tuple

import numpy as np
import sympy as sp


# ---------------------------------------------------------------------------
# Schema versioning
# ---------------------------------------------------------------------------

SCHEMA_VERSION = "v1.0"
"""Schema version for the cached closed-form table. Bump on any change
to the closed-form derivation, the polynomial encoding, or the
(alpha, B) variable convention. Cache is invalidated on mismatch."""

DEFAULT_CACHE_PATH = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "_cache",
    "hylleraas_eckart_closed_forms.pkl",
)


# ---------------------------------------------------------------------------
# Symbolic master integral
# ---------------------------------------------------------------------------

# Reserved sympy symbols (module-level so they round-trip through pickling).
# Both are declared positive: alpha by physics (decay rate); B by sympy
# convergence-tractability (with B real-but-not-necessarily-positive, the
# inner integrals return Piecewise expressions that block the outer
# s-integration). The closed form is even in B by t-parity of the
# integrand, so the result extends to any real B by substitution.
_ALPHA = sp.Symbol("alpha", positive=True, real=True)
_B = sp.Symbol("B", positive=True, real=True)


def I_HE_cosh_symbolic(L: int, M: int, N: int) -> sp.Expr:
    """Symbolic Hylleraas-Eckart cosh master integral.

    I_HE^cosh(L, M, N; alpha, B)
        = integral_0^infty ds exp(-2 alpha s) s^L
          * integral_0^s du u^{N+1}
          * integral_{-u}^u dt t^{2M} (s^2 - t^2) cosh(B t)

    Closes in elementary form (rational in (alpha, B)) under the
    convergence condition 2*alpha > B (B > 0 by convention here).

    Parameters
    ----------
    L, M, N : int
        Non-negative integers giving powers of s, t (squared), u in the
        integrand. L = l_p + l_q, etc. from the basis function combination.

    Returns
    -------
    sympy.Expr
        Rational function of _ALPHA and _B. Even in B by t-parity of
        the integrand.

    Notes
    -----
    SLOW: ~1-10 s for total degree L+2M+N <= 4, scaling roughly as
    (total_degree+5)^2.5. Use I_HE_cosh_numeric for production calls.

    Implementation: sympy returns a Piecewise (convergent branch +
    fallback). We extract the convergent branch (B < 2*alpha) and
    simplify it.
    """
    if L < 0 or M < 0 or N < 0:
        raise ValueError(f"L, M, N must be >= 0, got ({L}, {M}, {N})")

    s, u, t = sp.symbols("s u t", positive=True, real=True)
    # Inner-to-outer integration: t in [-u, u], u in [0, s], s in [0, inf].
    integrand = (
        sp.exp(-2 * _ALPHA * s)
        * s**L
        * t**(2 * M)
        * u**(N + 1)  # u^N from basis * u from Jacobian
        * (s**2 - t**2)
        * sp.cosh(_B * t)
    )
    r_t = sp.integrate(integrand, (t, -u, u))
    r_u = sp.integrate(r_t, (u, 0, s))
    r_s = sp.integrate(r_u, (s, 0, sp.oo))
    expr = _extract_convergent_branch(r_s)
    return sp.simplify(expr)


def _extract_convergent_branch(expr: sp.Expr) -> sp.Expr:
    """If expr is a Piecewise, pick the convergent branch (B < 2*alpha).

    Sympy returns the closed-form integral as a Piecewise with the
    convergent branch as the first piece (matching B < 2*alpha) and an
    unevaluated Integral fallback for B >= 2*alpha. We need the closed
    form, so we always take the first piece's expression.

    If expr is not a Piecewise, return as-is.
    """
    if isinstance(expr, sp.Piecewise):
        # First piece is the convergent closed form.
        return expr.args[0].expr
    return expr


def I_HE_cosh_b_zero_check(L: int, M: int, N: int) -> sp.Expr:
    """Symbolic difference I_HE^cosh(L,M,N; alpha, 0) - I_single(L,M,N; alpha).

    Returns simplified expression which should be exactly zero.
    """
    he = I_HE_cosh_symbolic(L, M, N)
    he_b0 = sp.simplify(he.subs(_B, 0))
    single = _single_alpha_master_symbolic(L, M, N)
    return sp.simplify(he_b0 - single)


def _single_alpha_master_symbolic(L: int, M: int, N: int) -> sp.Expr:
    """Existing single-alpha master I(L, M, N; alpha) in symbolic form.

    I(L, M, N; alpha) = 4 (N + 4M + 6) (L + N + 2M + 5)!
                       -----------------------------------
                       (2M+1)(2M+3)(N+2M+3)(N+2M+5)
                       * (2 alpha)^(-(L + N + 2M + 6))

    See geovac/hylleraas_r12.py::_hylleraas_master_int_factorial.
    """
    num = 4 * (N + 4 * M + 6) * sp.factorial(L + N + 2 * M + 5)
    den = (2 * M + 1) * (2 * M + 3) * (N + 2 * M + 3) * (N + 2 * M + 5)
    power = L + N + 2 * M + 6
    return sp.Rational(num, den) / (2 * _ALPHA) ** power


# ---------------------------------------------------------------------------
# Polynomial extraction (alpha^2, B^2)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ClosedForm:
    """Compact integer-rational representation of a master integral.

    The general closed form factors as

        I(L, M, N; alpha, B)
            = B^{b_factor} * numerator(alpha, B^2) / denominator(alpha, B^2)

    where numerator and denominator are polynomials in (a := alpha,
    b := B^2) with integer-rational coefficients, and `b_factor` is an
    explicit integer power of B that is factored out (typically 0 for
    cosh masters, 1 for sinh masters).

    For cosh-type masters (I, J, K, master_C_gen), the integrand is
    even in t and the integral is even in B, so b_factor = 0 and only
    even B-powers appear inside the polynomials.

    For sinh-type masters (master_S_gen), the integrand is also even in
    t (sinh*odd-t = even), but the integral is ODD in B (sinh changes
    sign under B -> -B). We factor out one explicit B to keep the inner
    polynomials in (alpha, B^2).

    The polynomials are stored as nested coefficient lists:

        numerator[i][j]   = coefficient of alpha^i * (B^2)^j in numerator
        denominator[i][j] = coefficient of alpha^i * (B^2)^j in denominator

    Both are Fractions for exact arithmetic.
    """

    L: int
    M: int
    N: int
    num_coeffs: Tuple[Tuple[Fraction, ...], ...]
    den_coeffs: Tuple[Tuple[Fraction, ...], ...]
    b_factor: int = 0  # extra explicit B^{b_factor} multiplier; 0 (cosh) or 1 (sinh)

    @property
    def num_degree_alpha(self) -> int:
        return len(self.num_coeffs) - 1

    @property
    def num_degree_b2(self) -> int:
        return max((len(row) for row in self.num_coeffs), default=0) - 1

    @property
    def den_degree_alpha(self) -> int:
        return len(self.den_coeffs) - 1

    @property
    def den_degree_b2(self) -> int:
        return max((len(row) for row in self.den_coeffs), default=0) - 1

    def evaluate(self, alpha: float, B: float) -> float:
        """Evaluate the closed form at numerical (alpha, B).

        Returns B^{b_factor} * num(alpha, B^2) / den(alpha, B^2).

        Convergence requires 2*alpha > |B|; we do NOT check that here
        (the polynomial evaluation is well-defined everywhere except at
        zeros of the denominator).
        """
        a = float(alpha)
        b = float(B) ** 2
        num = _eval_poly2d_float(self.num_coeffs, a, b)
        den = _eval_poly2d_float(self.den_coeffs, a, b)
        if den == 0.0:
            raise ZeroDivisionError(
                f"Master integral denominator vanishes at alpha={alpha}, B={B} "
                f"for (L={self.L}, M={self.M}, N={self.N}). Likely on a "
                f"branch where 2 alpha = |B| (convergence boundary)."
            )
        result = num / den
        if self.b_factor != 0:
            result *= float(B) ** self.b_factor
        return result

    def evaluate_exact(self, alpha: Fraction, B: Fraction) -> Fraction:
        """Evaluate at rational (alpha, B). Useful for tests and exact regression."""
        a = alpha
        b = B * B
        num = _eval_poly2d_exact(self.num_coeffs, a, b)
        den = _eval_poly2d_exact(self.den_coeffs, a, b)
        if den == 0:
            raise ZeroDivisionError(
                f"Master integral denominator vanishes at alpha={alpha}, B={B}"
            )
        result = num / den
        if self.b_factor != 0:
            result *= B ** self.b_factor
        return result


def _eval_poly2d_float(coeffs: Tuple[Tuple[Fraction, ...], ...],
                      a: float, b: float) -> float:
    """Horner-style evaluation of nested polynomial coeffs[i][j] * a^i b^j."""
    result = 0.0
    a_power = 1.0
    for row in coeffs:
        row_val = 0.0
        b_power = 1.0
        for c in row:
            row_val += float(c) * b_power
            b_power *= b
        result += a_power * row_val
        a_power *= a
    return result


def _eval_poly2d_exact(coeffs: Tuple[Tuple[Fraction, ...], ...],
                      a: Fraction, b: Fraction) -> Fraction:
    """Exact-rational polynomial evaluation."""
    result = Fraction(0)
    a_power = Fraction(1)
    for row in coeffs:
        row_val = Fraction(0)
        b_power = Fraction(1)
        for c in row:
            row_val += c * b_power
            b_power *= b
        result += a_power * row_val
        a_power *= a
    return result


def _sympy_to_closed_form(L: int, M: int, N: int,
                          expr: sp.Expr) -> ClosedForm:
    """Convert a sympy expression in (_ALPHA, _B) to a ClosedForm.

    Strategy:
      1. Factor as numerator/denominator over (_ALPHA, _B).
      2. Cast each as a sympy Poly in (_ALPHA, _B).
      3. Verify each is even in _B (only even-degree B monomials).
      4. Collapse alpha -> a := alpha (kept as is) and B^{2j} -> b^j (b := B^2).

    The result is two polynomials in (a, b) where a := alpha^2, b := B^2
    only AFTER we verify the alpha part also has consistent even-power
    structure (in our case, alpha appears to all powers, both even and
    odd, because of (2 alpha)^k factors; so we just use a := alpha rather
    than alpha^2).

    Compromise: store as polynomials in (a, b) = (alpha, B^2). Numerator
    and denominator are then integer-rational polynomials.
    """
    # Factor as numerator / denominator. We use sp.cancel which puts the
    # expression over a single common denominator AND cancels common factors,
    # ensuring the resulting numerator/denominator are polynomials in
    # (alpha, B). Without this, individual sub-fractions can have odd-B
    # factors in their denominators that only cancel when fully combined
    # (B-parity is a property of the sum, not of each rational sub-piece).
    expr_canceled = sp.cancel(expr)
    num_sym, den_sym = sp.fraction(expr_canceled)
    num_sym = sp.expand(num_sym)
    den_sym = sp.expand(den_sym)

    # Each must be a polynomial in (_ALPHA, _B) with even B-degree only.
    num_poly_ab = sp.Poly(num_sym, _ALPHA, _B)
    den_poly_ab = sp.Poly(den_sym, _ALPHA, _B)

    # Check evenness in B.
    for poly, name in [(num_poly_ab, "numerator"), (den_poly_ab, "denominator")]:
        for monom in poly.monoms():
            if monom[1] % 2 != 0:
                raise ValueError(
                    f"Closed form {name} for (L={L}, M={M}, N={N}) has "
                    f"odd power of B (monomial alpha^{monom[0]} * B^{monom[1]}); "
                    f"this violates the t-parity invariant of the master integral."
                )

    # Collapse B^2 -> b. Now polynomials are in (alpha, b) with all integer
    # exponents on b.
    num_coeffs = _poly_to_nested_alpha_b2(num_poly_ab)
    den_coeffs = _poly_to_nested_alpha_b2(den_poly_ab)

    return ClosedForm(L=L, M=M, N=N,
                      num_coeffs=num_coeffs,
                      den_coeffs=den_coeffs)


def _poly_to_nested_alpha_b2(poly: sp.Poly
                             ) -> Tuple[Tuple[Fraction, ...], ...]:
    """Convert sympy Poly in (alpha, B) (even in B) to nested Fraction
    tuple coeffs[i][j] = coefficient of alpha^i * (B^2)^j.

    Only even B-monomials contribute; the parity check is upstream.
    """
    # Build alpha-row dict, then materialize as nested tuple.
    coeff_dict: Dict[Tuple[int, int], Fraction] = {}
    for monom, c in poly.terms():
        i_alpha, j_B = monom
        if j_B % 2 != 0:
            # Defensive; should have been caught upstream.
            raise ValueError("odd B-monomial leaked through parity check")
        j_b2 = j_B // 2
        coeff_dict[(i_alpha, j_b2)] = _sympy_to_fraction(c)

    if not coeff_dict:
        return ((Fraction(0),),)

    max_i = max(k[0] for k in coeff_dict.keys())
    max_j = max(k[1] for k in coeff_dict.keys())

    rows = []
    for i in range(max_i + 1):
        row = []
        for j in range(max_j + 1):
            row.append(coeff_dict.get((i, j), Fraction(0)))
        # Trim trailing zeros (compact).
        while len(row) > 1 and row[-1] == 0:
            row.pop()
        rows.append(tuple(row))

    # Trim trailing zero rows.
    while len(rows) > 1 and all(c == 0 for c in rows[-1]):
        rows.pop()
    return tuple(rows)


def _sympy_to_fraction(c: sp.Expr) -> Fraction:
    """Convert a sympy rational to a Python Fraction."""
    r = sp.Rational(c)
    return Fraction(int(r.p), int(r.q))


# ---------------------------------------------------------------------------
# Cache (in-memory + disk)
# ---------------------------------------------------------------------------

# In-memory cache for the current Python process.
_CACHE_TABLE: Dict[Tuple[int, int, int], ClosedForm] = {}
_CACHE_LOADED_FROM_DISK = False


def I_HE_cosh_polynomial(L: int, M: int, N: int) -> ClosedForm:
    """Closed-form polynomial representation of I_HE^cosh(L, M, N) (overlap).

    Production path: routes through the algebraic recurrence engine
    (geovac/hylleraas_eckart_recurrence.py) on cache miss. The sympy
    oracle (I_HE_cosh_symbolic) is retained for testing/verification only,
    not in the production path.

    Disk cache (when present) and in-memory cache use the same recurrence-
    engine encoding, so all paths agree.

    Parameters
    ----------
    L, M, N : int
        Non-negative integers.

    Returns
    -------
    ClosedForm
        Integer-rational polynomial representation of the closed form.
    """
    key = (L, M, N)
    if key in _CACHE_TABLE:
        return _CACHE_TABLE[key]

    # Production path: recurrence engine (microseconds per cell).
    from geovac.hylleraas_eckart_recurrence import I_HE_cosh_polynomial_recurrence
    cf = I_HE_cosh_polynomial_recurrence(L, M, N)
    _CACHE_TABLE[key] = cf
    return cf


# In-memory caches for J and K master integrals. Not yet persisted to
# disk (matrix assembly populates them on demand at ~ms per call; if
# this proves too slow we'll add a separate disk-cache file).
_CACHE_TABLE_J: Dict[Tuple[int, int, int], ClosedForm] = {}
_CACHE_TABLE_K: Dict[Tuple[int, int, int], ClosedForm] = {}


def J_HE_cosh_polynomial(L: int, M: int, N: int) -> ClosedForm:
    """Closed-form polynomial representation of J_HE^cosh(L, M, N).

    Nuclear-attraction master integral (V_ne); integrand has no
    (s^2 - t^2) factor (cancelled by 1/r_1 + 1/r_2 = 4s/(s^2-t^2)).

    See geovac/hylleraas_eckart_recurrence.py::J_HE_cosh_polynomial_recurrence.
    """
    key = (L, M, N)
    if key in _CACHE_TABLE_J:
        return _CACHE_TABLE_J[key]
    from geovac.hylleraas_eckart_recurrence import J_HE_cosh_polynomial_recurrence
    cf = J_HE_cosh_polynomial_recurrence(L, M, N)
    _CACHE_TABLE_J[key] = cf
    return cf


def K_HE_cosh_polynomial(L: int, M: int, N: int) -> ClosedForm:
    """Closed-form polynomial representation of K_HE^cosh(L, M, N).

    Electron-electron repulsion master integral (V_ee); u-power is N
    (not N+1) because 1/r_12 = 1/u cancels the Jacobian u.

    See geovac/hylleraas_eckart_recurrence.py::K_HE_cosh_polynomial_recurrence.
    """
    key = (L, M, N)
    if key in _CACHE_TABLE_K:
        return _CACHE_TABLE_K[key]
    from geovac.hylleraas_eckart_recurrence import K_HE_cosh_polynomial_recurrence
    cf = K_HE_cosh_polynomial_recurrence(L, M, N)
    _CACHE_TABLE_K[key] = cf
    return cf


def I_HE_cosh_numeric(L: int, M: int, N: int,
                      alpha: float, B: float) -> float:
    """Fast numerical evaluator for I_HE^cosh(L, M, N; alpha, B).

    Uses the cached polynomial closed form. ~10 us per call after the
    cache is warm, vs ~1-10 s for the sympy-symbolic path.

    Parameters
    ----------
    L, M, N : int
        Non-negative integers.
    alpha : float
        Decay parameter (positive).
    B : float
        Asymmetry parameter (any sign; result is even in B).

    Returns
    -------
    float
        Value of the master integral.

    Raises
    ------
    ZeroDivisionError
        If 2*alpha = |B| (convergence boundary).
    ValueError
        If alpha <= 0.
    """
    if alpha <= 0:
        raise ValueError(f"alpha must be positive, got {alpha}")
    cf = I_HE_cosh_polynomial(L, M, N)
    return cf.evaluate(alpha, B)


# ---------------------------------------------------------------------------
# Precomputation driver
# ---------------------------------------------------------------------------

def precompute_table(max_total_degree: int = 8,
                     verbose: bool = False
                     ) -> Dict[Tuple[int, int, int], ClosedForm]:
    """Precompute closed forms for all (L, M, N) with L + 2M + N <= max_total_degree.

    Parameters
    ----------
    max_total_degree : int
        Maximum value of L + 2M + N. Default 8 (covers omega <= 4 basis pairs).
        For ω <= 8 production target, use 20.
    verbose : bool
        Print per-cell progress.

    Returns
    -------
    Dict[(L, M, N), ClosedForm]
        Computed closed forms. Also populates the in-memory _CACHE_TABLE.
    """
    table: Dict[Tuple[int, int, int], ClosedForm] = {}
    n_total = 0
    for L in range(max_total_degree + 1):
        for M in range((max_total_degree - L) // 2 + 1):
            for N in range(max_total_degree - L - 2 * M + 1):
                n_total += 1

    if verbose:
        print(f"[precompute_table] computing {n_total} cells "
              f"(max_total_degree={max_total_degree}) ...")

    t_start = time.time()
    count = 0
    for L in range(max_total_degree + 1):
        for M in range((max_total_degree - L) // 2 + 1):
            for N in range(max_total_degree - L - 2 * M + 1):
                key = (L, M, N)
                if key in _CACHE_TABLE:
                    table[key] = _CACHE_TABLE[key]
                    count += 1
                    continue

                t0 = time.time()
                cf = I_HE_cosh_polynomial(L, M, N)
                dt = time.time() - t0
                table[key] = cf
                count += 1

                if verbose:
                    print(f"  [{count}/{n_total}] (L={L}, M={M}, N={N}): "
                          f"deg_alpha={cf.num_degree_alpha}+{cf.den_degree_alpha}, "
                          f"deg_B^2={cf.num_degree_b2}+{cf.den_degree_b2}, "
                          f"{dt:.2f} s")

    elapsed = time.time() - t_start
    if verbose:
        print(f"[precompute_table] total {elapsed:.1f} s for {n_total} cells.")
    return table


def save_table(table: Dict[Tuple[int, int, int], ClosedForm],
               path: Optional[str] = None) -> str:
    """Save a closed-form table to disk.

    Parameters
    ----------
    table : dict
        Table from precompute_table.
    path : str, optional
        Output path. Default: DEFAULT_CACHE_PATH.

    Returns
    -------
    str
        Path the table was saved to.
    """
    if path is None:
        path = DEFAULT_CACHE_PATH

    os.makedirs(os.path.dirname(path), exist_ok=True)
    payload = {
        "schema_version": SCHEMA_VERSION,
        "table": table,
    }
    with open(path, "wb") as f:
        pickle.dump(payload, f)
    return path


def load_table(path: Optional[str] = None,
               populate_cache: bool = True
               ) -> Optional[Dict[Tuple[int, int, int], ClosedForm]]:
    """Load a closed-form table from disk.

    Parameters
    ----------
    path : str, optional
        Input path. Default: DEFAULT_CACHE_PATH.
    populate_cache : bool
        If True, also populate the in-memory _CACHE_TABLE.

    Returns
    -------
    dict or None
        Loaded table, or None if the file doesn't exist or has a stale
        schema version.
    """
    global _CACHE_LOADED_FROM_DISK

    if path is None:
        path = DEFAULT_CACHE_PATH
    if not os.path.exists(path):
        return None

    with open(path, "rb") as f:
        payload = pickle.load(f)

    if payload.get("schema_version") != SCHEMA_VERSION:
        # Stale cache; do not load.
        return None

    table = payload["table"]
    if populate_cache:
        _CACHE_TABLE.update(table)
        _CACHE_LOADED_FROM_DISK = True
    return table


# ---------------------------------------------------------------------------
# Module-level lazy load
# ---------------------------------------------------------------------------

def _try_load_default_cache() -> None:
    """Attempt to load the default disk cache on first import.

    Silent no-op if the cache doesn't exist or has a stale schema.
    """
    try:
        load_table()
    except Exception:
        # Never let a malformed cache file break import.
        pass


_try_load_default_cache()
