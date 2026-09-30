"""Paper 32, two corrections of /qa trunk FULL #10 (2026-09-30).

(1) The alpha-decomposition item "F is an infinite scalar spectral zeta":
    F = pi^2/6 is zeta_{1 - Delta_{S^3}}(2) = sum_{n>=1} n^2 * (n^2)^{-2}
    (eigenvalue n^2 of 1 - Delta on the degree-(n-1) harmonics, Fock
    degeneracy n^2).  It is NOT zeta_{-Delta_{S^3}}(2) =
    sum_{n>=2} n^2 / (n^2 - 1)^2 ~= 0.88497, which the paper printed
    before 2026-09-30.

    WRONG ANSWER REJECTED: "F = zeta_{Delta_{S^3}}(2)" -- the Delta sum is
    pinned to 0.884967... and asserted != pi^2/6.

(2) The pi-source case-exhaustion theorem's converse (iii) => (i) is
    WITHDRAWN (DELTA #21, 2026-09-30; an intermediate FULL #10 re-scope to
    'M1 or M2 only' is also retired).  Witnesses: M1's rational core B = 42
    (proven); M2's spinor zeta at half-integer s and M3 at s = 2 (analytic
    continuation) and s = 4 (convergent) -- no KNOWN pi.  The PSLQ helper
    below is EVIDENCE only and its scope is narrow: it tests q1*pi^k + q2
    for integer k in [-4, 4]; it is blind to sqrt(pi), polynomials in pi,
    and products of pi with other constants.

    WRONG ANSWER REJECTED: an M3 output at s = 2 that is a rational
    multiple of a power of pi (tested against pi^k, k = -4..4, with small
    rational coefficients via PSLQ at 50 digits).
"""
from __future__ import annotations

import mpmath as mp
import pytest


@pytest.fixture(autouse=True)
def _hp():
    old = mp.mp.dps
    mp.mp.dps = 50
    yield
    mp.mp.dps = old


def test_f_is_zeta_of_one_minus_laplacian():
    # 1 - Delta has eigenvalue n^2, degeneracy n^2, on degree-(n-1) harmonics
    zeta_1mD = mp.nsum(lambda n: n**2 * (n**2) ** (-2), [1, mp.inf])
    assert abs(zeta_1mD - mp.pi**2 / 6) < mp.mpf(10) ** -40


def test_f_is_not_zeta_of_laplacian():
    # -Delta has eigenvalue n^2 - 1 (zero mode n = 1 excluded)
    zeta_mD = mp.nsum(lambda n: n**2 / (n**2 - 1) ** 2, [2, mp.inf])
    assert abs(zeta_mD - mp.mpf("0.88496703342411321824")) < mp.mpf(10) ** -18
    assert abs(zeta_mD - mp.pi**2 / 6) > mp.mpf("0.7")


def _beta(s):
    """Dirichlet beta via Hurwitz zeta; beta(0) = 1/2."""
    return (mp.zeta(s, mp.mpf(1) / 4) - mp.zeta(s, mp.mpf(3) / 4)) / mp.mpf(4) ** s


def test_m3_output_at_s2_carries_no_pi():
    # PSLQ evidence, not proof (Catalan's G is not known to be irrational).
    assert abs(_beta(0) - mp.mpf(1) / 2) < mp.mpf(10) ** -40
    assert abs(_beta(2) - mp.catalan) < mp.mpf(10) ** -40
    s = 2
    m3 = mp.mpf(2) ** (s - 1) * (_beta(s) - _beta(s - 2))
    assert abs(m3 - 2 * (mp.catalan - mp.mpf(1) / 2)) < mp.mpf(10) ** -40
    # No small rational relation m3 = q * pi^k for k = -4..4 (with 1 as a
    # further basis element, i.e. m3 is not q1*pi^k + q2).
    for k in range(-4, 5):
        if k == 0:
            continue
        rel = mp.pslq([m3, mp.pi**k, mp.mpf(1)], maxcoeff=10**6, maxsteps=10**5)
        assert rel is None or rel[0] == 0, (k, rel)


def _no_small_pi_relation(x):
    """PSLQ EVIDENCE (not proof) that x is not q1*pi^k + q2 for k = -4..4
    with coefficients up to 1e6."""
    for k in range(-4, 5):
        if k == 0:
            continue
        rel = mp.pslq([x, mp.pi**k, mp.mpf(1)], maxcoeff=10**6, maxsteps=10**5)
        if rel is not None and rel[0] != 0:
            return False
    return True


def test_converse_withdrawn_m3_convergent_witness_s4():
    """DELTA #21 (2026-09-30): the converse (iii)=>(i) is WITHDRAWN.  M3's
    closed form 2^{s-1}(beta(s) - beta(s-2)) at the CONVERGENT point s = 4
    (the s = 2 value needs analytic continuation) is 8(beta(4) - G): checked
    against the direct Dirichlet sums, and PSLQ finds no small pi relation.
    WRONG ANSWER REJECTED: an M3 output at s = 4 that is q*pi^k + r."""
    direct_b4 = mp.nsum(lambda n: (-1) ** n / (2 * n + 1) ** 4, [0, mp.inf])
    assert abs(direct_b4 - _beta(4)) < mp.mpf(10) ** -40
    m3 = 8 * (_beta(4) - mp.catalan)
    assert abs(m3 - mp.mpf(2) ** 3 * (_beta(4) - _beta(2))) < mp.mpf(10) ** -40
    assert _no_small_pi_relation(m3)


def test_converse_withdrawn_m2_half_integer_witness():
    """DELTA #21: M2's spectral zeta (Paper 28 T9) at half-integer s = 5/2 is
    14 zeta(3) - (31/2) zeta(5) -- no pi -- so engaging M2 does not force pi.
    Computed from the Dirac spectrum |lambda| = n + 3/2, deg 2(n+1)(n+2).
    WRONG ANSWER REJECTED: an M2 spectral-zeta value at half-integer s with a
    small pi relation (i.e. 'M2 always injects pi')."""
    s = mp.mpf(5) / 2
    z = mp.nsum(lambda n: 2 * (n + 1) * (n + 2) * (n + mp.mpf(3) / 2) ** (-2 * s), [0, mp.inf])
    assert abs(z - (14 * mp.zeta(3) - mp.mpf(31) / 2 * mp.zeta(5))) < mp.mpf(10) ** -30
    assert _no_small_pi_relation(z)
