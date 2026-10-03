"""Paper 28 thm:parity scope note / Paper 24 corollary correction (2026-10-03,
/qa trunk DELTA #26 self-check): the integer-s pi^even property holds for
Delta_LB and D^2 (spectra n^2 - 1 and (n+3/2)^2, rational partial-fraction
poles), NOT for every second-order
operator.  Counterexample: the spinor connection Laplacian
nabla*nabla = D^2 - 3/2 on unit S^3 (Lichnerowicz, R/4 = 3/2), with

    zeta_{nabla*nabla}(2)
        = (1/12) (14 pi tan(pi sqrt6 / 2) / sqrt6 + 5 pi^2 sec^2(pi sqrt6 / 2)),

which is outside Q[pi^2].

What wrong answer each test rejects:
  * test_closed_form_matches_spectral_sum: a wrong closed form (any other
    rational coefficients, or a pi^even polynomial) -- the direct spectral
    sum pins it to 30 digits.
  * test_not_pi_even: the retired claim that the value lies in
    Q + Q pi^2 + Q pi^4 -- PSLQ over that basis must find no relation with
    small coefficients (the tan/sec relation is pinned by the first test).
"""

from __future__ import annotations

import mpmath as mp
import pytest


@pytest.fixture(autouse=True)
def _restore_mp_dps():
    """These tests raise mpmath's GLOBAL precision; restore it so other
    tests in the same xdist worker are not affected (a PSLQ test failed
    under contamination, 2026-10-03)."""
    old = mp.mp.dps
    yield
    mp.mp.dps = old


def _zeta_conn_lap(s: int) -> mp.mpf:
    # CH eigenvalues of D^2: (n+3/2)^2, multiplicity 2(n+1)(n+2) (2-spinor);
    # nabla*nabla = D^2 - 3/2.
    return mp.nsum(lambda n: 2 * (n + 1) * (n + 2)
                   / ((n + mp.mpf(3) / 2) ** 2 - mp.mpf(3) / 2) ** s, [0, mp.inf])


def _closed_form() -> mp.mpf:
    r6 = mp.sqrt(6)
    x = mp.pi * r6 / 2
    return (14 * mp.pi * mp.tan(x) / r6 + 5 * mp.pi ** 2 / mp.cos(x) ** 2) / 12


def test_closed_form_matches_spectral_sum() -> None:
    mp.mp.dps = 40
    assert abs(_zeta_conn_lap(2) - _closed_form()) < mp.mpf("1e-30")


def test_not_pi_even() -> None:
    mp.mp.dps = 40
    z = _zeta_conn_lap(2)
    rel = mp.pslq([z, 1, mp.pi ** 2, mp.pi ** 4], maxcoeff=10 ** 6, maxsteps=10 ** 6)
    assert rel is None
