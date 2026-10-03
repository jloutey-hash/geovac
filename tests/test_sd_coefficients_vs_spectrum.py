"""Independent-route guard: production Seeley-DeWitt coefficients vs the
exact Camporesi--Higuchi spectrum on the unit S^3.

Backs Paper 28 eq:seeley_dewitt / eq:sd_two_term and the correction of
2026-10-02 (/qa trunk DELTA #26): the production function
`geovac.qed_vacuum_polarization.seeley_dewitt_coefficients_s3` had the
Lichnerowicz endomorphism with the wrong sign, dropped it from a_1 and
omitted the spin-curvature term, returning the withdrawn values
[retracted 2026-10-02: sd-coeff-sign-error] a_0 = a_1 = sqrt(pi),
a_2 = sqrt(pi)/8 [retracted 2026-10-02: sd-coeff-sign-error].  The existing tests in test_qed_vacuum_polarization.py
re-derived those values from the same curvature polynomial, so they could
not see the error.

This guard uses a route that shares no code with the curvature
polynomials: the heat trace summed directly over the Dirac spectrum,
eigenvalues +-(n + 3/2) with multiplicity (n+1)(n+2) per sign per
2-component spinor (doubled for the 4-component convention the production
code uses).  On S^3 the asymptotic expansion terminates (Jacobi theta
inversion), so

    sum_n 4 (n+1)(n+2) exp(-t (n+3/2)^2)
        = a_0 t^{-3/2} + a_1 t^{-1/2} + a_2 t^{1/2} + O(exp(-pi^2/t)).

What wrong answer each test rejects:
  * test_trace_matches_production_coefficients: any production a_0, a_1,
    a_2 that does not reproduce the summed trace -- in particular the
    retired (sqrt(pi), sqrt(pi), sqrt(pi)/8), whose residual at t = 0.02 is
    ~ 1.5 sqrt(pi) t^{-1/2} ~ 19.
  * test_coefficients_extracted_from_spectrum: extracts a_0, a_1, a_2 from
    the summed trace alone (three small-t points, linear solve) and compares
    with production; a sign error in a_1 or a nonzero a_2 fails it.
"""

from __future__ import annotations

import mpmath as mp
import pytest
import sympy as sp

from geovac.qed_vacuum_polarization import seeley_dewitt_coefficients_s3


@pytest.fixture(autouse=True)
def _restore_mp_dps():
    """These tests raise mpmath's GLOBAL precision; restore it so other
    tests in the same xdist worker are not affected (a PSLQ test failed
    under contamination, 2026-10-03)."""
    old = mp.mp.dps
    yield
    mp.mp.dps = old


def _spectral_trace(t: mp.mpf) -> mp.mpf:
    """Tr exp(-t D^2) on unit S^3, 4-component spinors, summed directly."""
    return mp.nsum(lambda n: 4 * (n + 1) * (n + 2)
                   * mp.e ** (-t * (n + mp.mpf(3) / 2) ** 2), [0, mp.inf])


def _production() -> tuple:
    c = seeley_dewitt_coefficients_s3()
    return tuple(mp.mpf(str(sp.N(c[k], 40))) for k in ("a0", "a1", "a2"))


def test_trace_matches_production_coefficients() -> None:
    mp.mp.dps = 40
    a0, a1, a2 = _production()
    for t in (mp.mpf("0.01"), mp.mpf("0.02"), mp.mpf("0.05")):
        lhs = _spectral_trace(t)
        rhs = a0 * t ** mp.mpf(-1.5) + a1 * t ** mp.mpf(-0.5) + a2 * t ** mp.mpf(0.5)
        assert abs(lhs - rhs) < mp.mpf("1e-25"), (t, lhs, rhs)


def test_coefficients_extracted_from_spectrum() -> None:
    mp.mp.dps = 40
    ts = [mp.mpf("0.01"), mp.mpf("0.02"), mp.mpf("0.04")]
    A = mp.matrix([[t ** mp.mpf(-1.5), t ** mp.mpf(-0.5), t ** mp.mpf(0.5)] for t in ts])
    b = mp.matrix([_spectral_trace(t) for t in ts])
    x = mp.lu_solve(A, b)
    a0, a1, a2 = _production()
    assert abs(x[0] - a0) < mp.mpf("1e-20")
    assert abs(x[1] - a1) < mp.mpf("1e-20")
    assert abs(x[2] - a2) < mp.mpf("1e-20")
    # and the spectrum's own values, independent of production:
    assert abs(x[0] - mp.sqrt(mp.pi)) < mp.mpf("1e-20")
    assert abs(x[1] + mp.sqrt(mp.pi) / 2) < mp.mpf("1e-20")
    assert abs(x[2]) < mp.mpf("1e-20")
