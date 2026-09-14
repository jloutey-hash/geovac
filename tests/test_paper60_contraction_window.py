"""Paper 60: the contraction window behind the two-center conditioning law.

Backs the two paragraphs added after `eq:chirp_decay` (CHANGELOG v5.11.4):

  (i)  the pi^2 of `eq:sigma_law` is TRUNCATION-side, not Bessel content --
       the same constant appears in a band-limited concentration problem with
       no Bessel function anywhere in it, and the SW near-null direction is
       exactly the vector that attains it;
  (ii) the law is independent of any smooth positive radial weight on the Fock
       sphere, WITH the vanishing-weight control that makes (ii) a measurement
       rather than an insensitivity.

Written as a separate pass from the edits it protects (CLAUDE.md Sec. 9).  Each
guard names, in its docstring, the specific WRONG ANSWER it exists to reject;
all are fire-tested via debug/qa/fire_test.py.
"""
from __future__ import annotations

import numpy as np
import pytest

from geovac.sturmian_sigma_law import (
    COLLAPSE_CONSTANT,
    generalized_sigma_max,
    sw_cross_block,
    theta2_band_matrix,
    weighted_blocks,
)

M_FAST = 120001
KR = 2.0


def _richardson(vals: list[float]) -> float:
    """First-order Richardson in 1/n on a doubling ladder."""
    return 2.0 * vals[-1] - vals[-2]


def _band_min(n: int) -> float:
    return float(np.min(np.linalg.eigvalsh(theta2_band_matrix(n))))


def _theta_moment2(v: np.ndarray, n: int, M: int = 200001) -> float:
    """<theta^2> of f(theta) = sum_a v_a sin(a chi), theta = pi - chi."""
    chi = np.linspace(1e-9, np.pi, M)
    a = np.arange(1, n + 1)
    f = v @ np.sin(np.outer(a, chi))
    w = f ** 2
    return float(np.trapezoid(w * (np.pi - chi) ** 2, chi) / np.trapezoid(w, chi))


# --------------------------------------------------------------- anchors
def test_theta2_closed_form_matches_quadrature():
    """ANCHOR. Rejects: a mistyped closed form that nothing else would catch.

    Every downstream claim reads the closed form, so if its Hankel term or its
    diagonal were wrong the whole file would be measuring a different matrix
    and agreeing with itself.
    """
    n, M = 60, 200001
    chi = np.linspace(0.0, np.pi, M)
    wq = np.full(M, chi[1] - chi[0])
    wq[0] *= 0.5
    wq[-1] *= 0.5
    a = np.arange(1, n + 1)
    smat = np.sin(np.outer(a, chi))
    quad = (2.0 / np.pi) * (smat * (wq * (np.pi - chi) ** 2)) @ smat.T
    assert np.max(np.abs(theta2_band_matrix(n) - quad)) < 1e-10


def test_unit_weight_reproduces_the_tracked_sw_block():
    """ANCHOR. Rejects: the weighted pipeline being a separate universe.

    Without this, the weight-independence agreement below could hold inside a
    self-consistent but wrong construction that never touches the SW metric.
    """
    n = 40
    A, B = weighted_blocks(KR, n, weight=None, M=200001)
    assert np.max(np.abs(A - np.eye(n))) < 1e-10
    assert np.max(np.abs(B - sw_cross_block(KR, n, M=200001))) < 1e-9


# ------------------------------------------- (i) the pi^2 is truncation-side
def test_band_minimum_is_pi_squared_with_no_bessel_present():
    """Rejects: "the pi^2 of eq:sigma_law is Bessel/j0 content".

    This computation contains no Bessel function of any kind -- only the
    geometric matrix <theta^2> -- and still returns pi^2.  If the constant were
    carried by j0 it could not appear here.
    """
    ns = [80, 160, 320, 640]
    vals = [_band_min(n) * n * n for n in ns]
    assert vals == sorted(vals), "n^2 lam_min must rise monotonically to pi^2"
    assert all(v < np.pi ** 2 for v in vals), "must approach pi^2 from below"
    assert abs(_richardson(vals) - np.pi ** 2) < 2e-3


def test_the_constant_tracks_the_zero_order_not_the_value_pi_squared():
    """Rejects: a guard that returns pi^2 by construction rather than measuring.

    Raising the symbol's zero from quadratic to quartic must move the constant
    well away from pi^2; if it does not, the test above proves nothing about
    where pi^2 comes from.
    """
    n, M = 200, 200001
    chi = np.linspace(0.0, np.pi, M)
    wq = np.full(M, chi[1] - chi[0])
    wq[0] *= 0.5
    wq[-1] *= 0.5
    a = np.arange(1, n + 1)
    smat = np.sin(np.outer(a, chi))
    quartic = (2.0 / np.pi) * (smat * (wq * (np.pi - chi) ** 4)) @ smat.T
    lam4 = float(np.min(np.linalg.eigvalsh(quartic)))
    # quartic zero => lam_min ~ c_2/n^4, so n^2 lam_min must COLLAPSE, not sit at pi^2
    assert lam4 * n ** 2 < 0.1 * np.pi ** 2


def test_near_null_direction_attains_the_band_minimum():
    """Rejects: "the degeneracy is a diffuse property spread over the sphere".

    The SW near-null direction must be the OPTIMAL concentrator at the p = 0
    pole -- n*rms(theta) -> pi -- and a lower singular direction must be
    materially worse.  Asserting only the first would accept a construction in
    which every direction looks the same.
    """
    ns, tops = [40, 80, 160], []
    second = None
    for n in ns:
        U, _, _ = np.linalg.svd(sw_cross_block(KR, n, M=M_FAST))
        tops.append(np.sqrt(_theta_moment2(U[:, 0], n)) * n)
        if n == ns[-1]:
            second = np.sqrt(_theta_moment2(U[:, 1], n)) * n
    assert tops == sorted(tops), "n*rms must rise toward pi"
    assert abs(_richardson(tops) - np.pi) < 0.02
    assert second > 1.5 * tops[-1], (
        f"a lower singular direction must be far less concentrated "
        f"(got {second:.3f} vs top {tops[-1]:.3f})")


def test_factorisation_holds_on_the_near_null_direction():
    """Rejects: "1-sigma_max and <theta^2> are two unrelated numbers".

    The paragraph asserts 1-sigma_max = (kR)^2 <theta^2>/24 ON that direction,
    which is what licenses reading the conditioning law as a window statement.
    """
    for kr in (0.5, 2.0, 4.0):
        n = 160
        U, sig, _ = np.linalg.svd(sw_cross_block(kr, n, M=M_FAST))
        lhs = 1.0 - float(sig[0])
        rhs = (kr ** 2 / 24.0) * _theta_moment2(U[:, 0], n)
        assert abs(lhs / rhs - 1.0) < 0.02, f"kR={kr}: {lhs:.6e} vs {rhs:.6e}"


def test_minimiser_is_the_dirichlet_mode_with_the_antipodal_parity():
    """Rejects: "the minimiser is the Dirichlet ground state", stated bare.

    Paper 60 says the band minimiser IS that mode, which is the sentence that
    makes the Bessel-free measurement a change of representation rather than an
    independent route.  It is true only with the alternation: the basis index
    is chi while the Dirichlet mode is natural in theta, and
    sin(a chi) = (-1)^(a+1) sin(a theta).

    Both halves are asserted, because the second is what carries the claim --
    WITHOUT the parity factor the two vectors are exactly ORTHOGONAL, so a
    guard checking only "correlates with sin(pi a/(n+1))" would fail, and one
    checking only the alternating form would not show the factor is load
    bearing.
    """
    for n in (80, 320):
        v = np.linalg.eigh(theta2_band_matrix(n))[1][:, 0]
        v = v / np.linalg.norm(v)
        a = np.arange(1, n + 1)
        plain = np.sin(np.pi * a / (n + 1))
        plain = plain / np.linalg.norm(plain)
        alt = ((-1.0) ** (a + 1)) * plain
        assert abs(float(v @ alt)) > 1 - 1e-4, (
            f"n={n}: minimiser must be the ALTERNATING Dirichlet mode")
        assert abs(float(v @ plain)) < 1e-3, (
            f"n={n}: the unalternated mode must be orthogonal — if it is not, "
            f"the antipodal parity factor is not doing anything")


def test_near_null_direction_is_that_same_mode():
    """Rejects: "the band minimiser and the SW near-null vector are unrelated".

    The paragraph's argument needs the cross block's near-null direction to BE
    the concentrator, not merely to share its constant.
    """
    n = 160
    U, _, _ = np.linalg.svd(sw_cross_block(KR, n, M=M_FAST))
    v = U[:, 0] / np.linalg.norm(U[:, 0])
    a = np.arange(1, n + 1)
    alt = ((-1.0) ** (a + 1)) * np.sin(np.pi * a / (n + 1))
    alt = alt / np.linalg.norm(alt)
    assert abs(float(v @ alt)) > 1 - 1e-4


# ---------------------------------------- (ii) weight-independence + control
SMOOTH = {
    "W=1": None,
    "W=1+0.8cos": lambda c: 1.0 + 0.8 * np.cos(c),
    "W=2+sin": lambda c: 2.0 + np.sin(c),
    "W=exp(-chi)": lambda c: np.exp(-c),
}


def _collapse(weight, n: int, kr: float = KR) -> float:
    A, B = weighted_blocks(kr, n, weight=weight, M=M_FAST)
    return (1.0 - generalized_sigma_max(A, B)) * (n / kr) ** 2


@pytest.mark.parametrize("label", list(SMOOTH))
def test_smooth_weights_preserve_the_collapse_constant(label):
    """The hopeful half. Rejects: a weight-dependent law.

    On its own this guard is NOT sufficient -- an insensitive pipeline would
    also pass it.  The control below is what makes it evidence.
    """
    got = _collapse(SMOOTH[label], 160)
    assert abs(got - COLLAPSE_CONSTANT) < 0.012, f"{label}: {got:.5f}"


def test_vanishing_weight_control_moves_the_constant():
    """THE LOAD-BEARING HALF. Rejects: "the agreement above is an insensitivity".

    A weight that vanishes AT the degeneracy (W = 1+cos(chi) -> 0 at chi = pi)
    must move the constant far outside the band the smooth weights occupy.  If
    this passes with the control agreeing, the measurement means nothing.
    """
    got = _collapse(lambda c: 1.0 + np.cos(c), 160)
    assert got > 2.0 * COLLAPSE_CONSTANT, (
        f"control must break the collapse, got {got:.5f} against "
        f"pi^2/24 = {COLLAPSE_CONSTANT:.5f}")


def test_weight_independence_holds_at_the_exponent_too():
    """Rejects: "only the constant survives; the exponent drifts with W".

    The paragraph claims BOTH. Fitted exponent must sit near -2 for every
    weight, control included.
    """
    for label, W in list(SMOOTH.items()) + [
            ("control", lambda c: 1.0 + np.cos(c))]:
        ns = (40, 80, 160)
        d = [(1.0 - generalized_sigma_max(*weighted_blocks(KR, n, W, M_FAST)))
             for n in ns]
        expo = float(np.polyfit(np.log(np.array(ns, float)), np.log(d), 1)[0])
        assert abs(expo + 2.0) < 0.06, f"{label}: exponent {expo:+.3f}"
