"""Category 3 -- one-centre Slater radial repulsion integrals as EXACT RATIONALS.

These are the cleanest entries in the table: they carry no digit count at all,
because they are not decimal approximations of anything.  The Slater integral

    R^k(n1 l1, n3 l3; n2 l2, n4 l4)
        = int int R_{n1l1}(r1) R_{n3l3}(r1) (r_<^k / r_>^{k+1})
                  R_{n2l2}(r2) R_{n4l4}(r2) r1^2 r2^2 dr1 dr2

has a polynomial-times-exponential integrand with rational coefficients at unit
orbital exponent, so the whole double integral collapses to a rational number.
``geovac.hypergeometric_slater.compute_rk_algebraic`` returns it as a Python
``Fraction``: infinite precision, no rounding, no convergence question.

Two independent verifications are recorded per entry:

  1. ROUTE CHECK.  A nested adaptive quadrature of the DEFINING double integral,
     built here from ``scipy.special.eval_genlaguerre`` and the textbook
     normalisation -- it never forms a Laguerre product expansion or a T-kernel,
     so it shares no algebra with the exact route.  Agreement is reported as an
     absolute difference (typically ~1e-16, the float64 floor).
  2. FLOAT-PATH CHECK.  ``compute_rk_float`` for max(n) <= 4, which is a
     genuinely separate pure-float implementation.  For max(n) >= 5 that function
     DISPATCHES to the exact Fraction path and casts, so it is NOT an independent
     check there; the entry says so explicitly rather than quietly counting it.

The Z-scaling law is R^k(Z) = Z * R^k(Z=1); every value below is at Z = 1.
"""
from __future__ import annotations

import warnings
from fractions import Fraction
from typing import Any, Dict, List, Tuple

import numpy as np
from scipy import integrate, special

from geovac.hypergeometric_slater import compute_rk_algebraic, compute_rk_float

from ._common import entry, frac_str

#: (n1,l1,n3,l3,n2,l2,n4,l4,k) quartets, chosen to span the float path
#: (max n <= 4) and the exact-Fraction dispatch path (max n >= 5), and to cover
#: k = 0 through k = 4.
QUARTETS: List[Tuple[Tuple[int, ...], str]] = [
    ((1, 0, 1, 0, 1, 0, 1, 0, 0), "R^0(1s,1s;1s,1s) -- the helium direct integral"),
    ((2, 0, 2, 0, 2, 0, 2, 0, 0), "R^0(2s,2s;2s,2s)"),
    ((2, 1, 2, 1, 2, 1, 2, 1, 0), "R^0(2p,2p;2p,2p)"),
    ((2, 1, 2, 1, 2, 1, 2, 1, 2), "R^2(2p,2p;2p,2p)"),
    ((1, 0, 2, 0, 1, 0, 2, 0, 0), "R^0(1s,2s;1s,2s)"),
    ((3, 2, 3, 2, 3, 2, 3, 2, 4), "R^4(3d,3d;3d,3d)"),
    ((2, 0, 3, 1, 2, 0, 3, 1, 1), "R^1(2s,3p;2s,3p)"),
    ((5, 0, 5, 0, 5, 0, 5, 0, 0), "R^0(5s,5s;5s,5s) -- exact-Fraction dispatch path"),
    ((6, 1, 6, 1, 6, 1, 6, 1, 2), "R^2(6p,6p;6p,6p) -- exact-Fraction dispatch path"),
]

_FLOAT_PATH_INDEPENDENT_MAX_N = 4


def _R_nl(n: int, l: int, r):
    """Hydrogenic radial function at unit orbital exponent, textbook normalisation."""
    norm = np.sqrt((2.0 / n) ** 3 * special.factorial(n - l - 1)
                   / (2 * n * special.factorial(n + l)))
    x = 2.0 * r / n
    return norm * x ** l * np.exp(-r / n) * special.eval_genlaguerre(n - l - 1,
                                                                     2 * l + 1, x)


def _rk_quadrature(n1, l1, n3, l3, n2, l2, n4, l4, k) -> float:
    """Nested adaptive quadrature of the DEFINING Slater double integral."""
    def inner(r1: float) -> float:
        lo, _ = integrate.quad(
            lambda r2: _R_nl(n2, l2, r2) * _R_nl(n4, l4, r2) * r2 ** (k + 2),
            0.0, r1, limit=300, epsabs=1e-15, epsrel=1e-13)
        hi, _ = integrate.quad(
            lambda r2: _R_nl(n2, l2, r2) * _R_nl(n4, l4, r2) * r2 ** (1 - k),
            r1, np.inf, limit=300, epsabs=1e-15, epsrel=1e-13)
        return (_R_nl(n1, l1, r1) * _R_nl(n3, l3, r1) * r1 ** 2
                * (lo / r1 ** (k + 1) + hi * r1 ** k))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        val, _ = integrate.quad(inner, 0.0, np.inf, limit=300,
                                epsabs=1e-14, epsrel=1e-12)
    return val


def build(mode: str = "full") -> List[Dict[str, Any]]:
    fast = (mode == "fast")
    rows: List[Dict[str, Any]] = []

    for q, label in QUARTETS:
        n1, l1, n3, l3, n2, l2, n4, l4, k = q
        exact = compute_rk_algebraic(*q)
        assert isinstance(exact, Fraction)
        max_n = max(n1, n3, n2, n4)

        flt = compute_rk_float(*q)
        float_delta = abs(float(exact) - flt)
        if max_n <= _FLOAT_PATH_INDEPENDENT_MAX_N:
            float_note = (f"independent pure-float implementation "
                          f"(compute_rk_float) gives {flt!r}, differing from the "
                          f"exact rational by {float_delta:.2e} -- float64 "
                          f"round-off only")
        else:
            float_note = (f"compute_rk_float dispatches to the exact Fraction "
                          f"path at max(n) = {max_n} and casts, so it is NOT an "
                          f"independent check here (it returns {flt!r}, delta "
                          f"{float_delta:.2e} = the float cast)")

        if fast:
            quad_note = "(independent quadrature skipped in fast mode)"
            quad_delta = None
        else:
            ref = _rk_quadrature(*q)
            quad_delta = abs(float(exact) - ref)
            quad_note = (f"independent nested quadrature of the defining double "
                         f"integral gives {ref:.16g}, differing from the exact "
                         f"rational by {quad_delta:.2e} (float64 floor)")

        rows.append(entry(
            f"slater.R{k}.{n1}{l1}{n3}{l3}_{n2}{l2}{n4}{l4}",
            "slater_rational",
            label,
            frac_str(exact),
            "exact",
            "Exact rational arithmetic, "
            "geovac.hypergeometric_slater.compute_rk_algebraic: the associated "
            "Laguerre polynomials are expanded with exact Fraction coefficients, "
            "the pair products are collected as polynomial-times-exponential "
            "terms, and the ordered double integral is evaluated term by term as "
            "incomplete Gamma functions at rational arguments.  No floating point "
            "anywhere in the route; the result is a Python Fraction.",
            f"EXACT -- infinitely many correct digits, because the value is a "
            f"rational number and not an approximation.  Decimal form "
            f"{float(exact)!r}.  Two verifications: (i) {quad_note}; "
            f"(ii) {float_note}.  Backing test: "
            f"tests/test_hypergeometric_slater.py.",
            value_kind="rational",
            decimal_value=repr(float(exact)),
            transcendence_class="{} (rational -- no transcendental content)",
            quadrature_abs_agreement=quad_delta,
            float_path_abs_delta=float_delta,
            float_path_independent=bool(max_n <= _FLOAT_PATH_INDEPENDENT_MAX_N),
            backing_test="tests/test_hypergeometric_slater.py",
            geovac_entry_point="geovac.hypergeometric_slater.compute_rk_algebraic",
        ))

    return rows
