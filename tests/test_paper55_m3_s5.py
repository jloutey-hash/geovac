"""Backing test for Paper 55 S^5 M3 results: thm:m3_s5_parity, thm:m3_s5_chi4,
cor:m3_s5_cyclotomic.

Promoted from debug/sprint_a6_m3_s5_derivation.py so the S^5 specialisation of
the M3 cyclotomic-mixed-Tate classification carries tracked, regression-
protected backing instead of an in-proof debug/ citation (Clean-Room policy,
group3 re-cert 2026-08-24).

Genuine content (two independent routes, not self-comparison):
  * The S^5 Dirac spectral zeta D(s) = sum_n g_n |lambda_n|^{-s} with
    g(m)=(m^2-9/4)(m^2-1/4)/3, m=n+5/2, is computed BOTH by a direct
    high-precision spectral sum AND by the Hurwitz closed form; the two must
    agree.  A wrong degeneracy g(m) or a wrong closed form fails.
  * The parity discriminant D_even - D_odd carries the Dirichlet-beta
    L-value chi_{-4} = beta(s) (level-4 content), which does NOT reduce to the
    pure-pi (level <=2) ring -- the cor:m3_s5_cyclotomic marker.
"""
from __future__ import annotations

import mpmath as mp
import pytest

DPS = 60


def g_of_m(m):
    """S^5 Dirac shell degeneracy as a polynomial in m=n+5/2 (a6 derivation)."""
    return (m**2 - mp.mpf(9) / 4) * (m**2 - mp.mpf(1) / 4) / 3


def D_s5_direct(s, N=400000):
    """Direct spectral sum sum_{n>=0} g(n+5/2) (n+5/2)^{-s}, Euler-Maclaurin tail."""
    mp.mp.dps = DPS
    s = mp.mpf(s)
    # partial sum + analytic tail via mpmath.sumem (Euler-Maclaurin) is unstable
    # for polynomial*power; use nsum which handles smooth g(m) m^{-s} cleanly.
    return mp.nsum(lambda n: g_of_m(n + mp.mpf(5) / 2) * (n + mp.mpf(5) / 2) ** (-s),
                   [0, mp.inf])


def Z_full_closed(s):
    """Z(s)=sum_{m in Z+1/2, m>=5/2} m^{-s} = (2^s-1) zeta_R(s) - 2^s - (2/3)^s."""
    return (mp.mpf(2) ** s - 1) * mp.zeta(s) - mp.mpf(2) ** s - (mp.mpf(2) / 3) ** s


def D_s5_closed(s):
    """Hurwitz closed form: (1/3)Z(s-4) - (5/6)Z(s-2) + (3/16)Z(s)."""
    s = mp.mpf(s)
    return (mp.mpf(1) / 3 * Z_full_closed(s - 4)
            - mp.mpf(5) / 6 * Z_full_closed(s - 2)
            + mp.mpf(3) / 16 * Z_full_closed(s))


def D_diff_closed(s):
    """Closed form for D_even - D_odd via diff_Z(s) = 2^s beta(s) - 2^s + (2/3)^s."""
    s = mp.mpf(s)

    def diff_Z(sv):
        return mp.mpf(2) ** sv * _beta(sv) - mp.mpf(2) ** sv + (mp.mpf(2) / 3) ** sv

    return (mp.mpf(1) / 3 * diff_Z(s - 4) - mp.mpf(5) / 6 * diff_Z(s - 2)
            + mp.mpf(3) / 16 * diff_Z(s))


def D_diff_direct(s, N=400000):
    """Direct (D_even - D_odd): n even (m=5/2,9/2,...) minus n odd (m=7/2,11/2,...)."""
    mp.mp.dps = DPS
    s = mp.mpf(s)
    even = mp.nsum(lambda k: g_of_m(2 * k + mp.mpf(5) / 2) * (2 * k + mp.mpf(5) / 2) ** (-s),
                   [0, mp.inf])
    odd = mp.nsum(lambda k: g_of_m(2 * k + 1 + mp.mpf(5) / 2) * (2 * k + 1 + mp.mpf(5) / 2) ** (-s),
                  [0, mp.inf])
    return even - odd


def _beta(s):
    """Dirichlet beta L(s, chi_{-4}) = 4^{-s}[zeta(s,1/4) - zeta(s,3/4)]."""
    return mp.mpf(4) ** (-s) * (mp.zeta(s, mp.mpf(1) / 4) - mp.zeta(s, mp.mpf(3) / 4))


@pytest.mark.parametrize("s", [7, 8, 9, 10])
def test_direct_sum_equals_closed_form(s):
    """S^5 Dirac zeta: direct spectral sum == Hurwitz closed form (two routes)."""
    mp.mp.dps = DPS
    direct = D_s5_direct(s)
    closed = D_s5_closed(s)
    assert abs(direct - closed) < mp.mpf(10) ** (-(DPS - 20)), \
        f"s={s}: direct {direct} != closed {closed}"


@pytest.mark.parametrize("s", [7, 8, 9, 10])
def test_parity_discriminant_carries_chi4_beta(s):
    """D_even - D_odd == closed form built from Dirichlet beta (chi_{-4} marker)."""
    mp.mp.dps = DPS
    direct = D_diff_direct(s)
    closed = D_diff_closed(s)
    assert abs(direct - closed) < mp.mpf(10) ** (-(DPS - 20)), \
        f"s={s}: parity direct {direct} != beta-closed {closed}"


@pytest.mark.parametrize("s", [8, 10])
def test_beta_is_genuinely_chi4_L_value(s):
    """cor:m3_s5_cyclotomic marker: the transcendental in the parity discriminant
    is L(s, chi_{-4}) = beta(s) -- verified by matching the Hurwitz form used in
    the closed form above to the defining alternating odd-integer sum
    sum_{j>=0} (-1)^j (2j+1)^{-s} (two independent routes to the chi_{-4} L-value)."""
    mp.mp.dps = DPS
    hurwitz = _beta(mp.mpf(s))
    alt_sum = mp.nsum(lambda j: mp.mpf(-1) ** j / (2 * j + 1) ** mp.mpf(s), [0, mp.inf])
    assert abs(hurwitz - alt_sum) < mp.mpf(10) ** (-(DPS - 20)), \
        f"beta({s}) Hurwitz form != chi_{{-4}} alternating sum ({hurwitz} vs {alt_sum})"


def test_level4_content_is_real_and_nonzero():
    """The parity split is genuine: D_even - D_odd is nonzero and equals the
    chi_{-4}-bearing closed form (i.e. the level-4 content does not cancel)."""
    mp.mp.dps = DPS
    d = D_diff_direct(8)
    assert abs(d) > mp.mpf(10) ** (-6), "parity discriminant vanished (split not real)"
    assert abs(d - D_diff_closed(8)) < mp.mpf(10) ** (-(DPS - 20))
    # positive control: beta at ODD argument IS a pi-power (beta(3)=pi^3/32),
    # confirming the chi_{-4} L-function has the expected functional structure.
    assert abs(_beta(3) - mp.pi ** 3 / 32) < mp.mpf(10) ** (-(DPS - 20))


def test_nontautology_guard_wrong_degeneracy_fails():
    """A perturbed degeneracy g(m)+1 must break the direct==closed agreement."""
    mp.mp.dps = DPS
    s = mp.mpf(9)
    wrong = mp.nsum(lambda n: (g_of_m(n + mp.mpf(5) / 2) + 1) * (n + mp.mpf(5) / 2) ** (-s),
                    [0, mp.inf])
    assert abs(wrong - D_s5_closed(9)) > mp.mpf(10) ** (-6), \
        "wrong degeneracy still matched the closed form -- test is not discriminating"
