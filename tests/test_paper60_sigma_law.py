"""Paper 60 -- the principal-angle (sigma-spectrum) law for the SW metric.

Backs the derived conditioning law (CHANGELOG v4.103.0; Paper 60 sec:molecular):

  1. spec(S) = {1 +/- sigma_k} and cond(S) = (1+sigma_max)/(1-sigma_max) exactly
     (SW intra block = identity), and the v4.73.0 composition-wall commutator is
     the companion functional max_k sigma_k sqrt(1-sigma_k^2) of the SAME spectrum.
  2. 1 - sigma_max = (s^2/24) pi^2 / n^2 asymptotically: the collapse variable
     x = n/s satisfies (1-sigma_max) x^2 -> pi^2/24, so the conditioning
     exponent is exactly 2; the fitted N^1.85 window slope is pre-asymptotic.
  3. cond(gerade) -> 2/(1 + min j0) = 2.555041..., R- and N-independent.
  4. The commutator norm saturates at 1/2 (some sigma_k near 1/sqrt2) once the
     spectrum densifies -- it carries no conditioning-severity information.
"""
import numpy as np
import pytest

from geovac.sturmian_sigma_law import (
    COLLAPSE_CONSTANT,
    assemble_two_center,
    commutator_direct,
    commutator_from_sigma,
    cond_from_sigma,
    gerade_constant,
    sigma_spectrum,
    sw_cross_block,
)


@pytest.mark.parametrize("s,nmax", [(1.4, 4), (2.0, 6), (4.0, 10)])
def test_spectrum_and_cond_identity(s, nmax):
    C = sw_cross_block(s, nmax, M=400001)
    S = assemble_two_center(C)
    sig = sigma_spectrum(C)
    ev = np.sort(np.linalg.eigvalsh(S))
    pred = np.sort(np.concatenate([1.0 - sig, 1.0 + sig]))
    assert np.max(np.abs(ev - pred)) < 1e-12
    cond = np.linalg.cond(S)
    assert abs(cond - cond_from_sigma(sig)) / cond < 1e-10


@pytest.mark.parametrize("s,nmax", [(1.4, 10), (2.0, 6)])
def test_commutator_identity(s, nmax):
    C = sw_cross_block(s, nmax, M=400001)
    sig = sigma_spectrum(C)
    assert abs(commutator_direct(C) - commutator_from_sigma(sig)) < 1e-10


def test_sigma_law_collapse():
    """(1 - sigma_max) (n/s)^2 -> pi^2/24 monotonically from below."""
    s = 2.0
    ratios = []
    for n in (40, 80, 120):
        sig = sigma_spectrum(sw_cross_block(s, n, M=150001))
        ratios.append((1.0 - sig.max()) * (n / s) ** 2 / COLLAPSE_CONSTANT)
    assert ratios[0] > 0.94 and ratios[1] > 0.97 and ratios[2] > 0.985
    assert ratios[0] < ratios[1] < ratios[2] < 1.0


def test_window_slope_is_preasymptotic():
    """The Paper-60 fit window (N = 4..24, s = 1.4) reads a slope well below 2,
    while the local slope at n = 40 -> 80 is already near 2."""
    s = 1.4
    ns = np.arange(2, 13)
    y = []
    for n in ns:
        sig = sigma_spectrum(sw_cross_block(s, int(n), M=200001))
        y.append((1.0 + sig.max()) / (1.0 - sig.max()))
    beta = np.polyfit(np.log(2.0 * ns), np.log(y), 1)[0]
    assert 1.70 < beta < 1.95     # the published 1.85 is a window reading
    sig40 = sigma_spectrum(sw_cross_block(s, 40, M=150001))
    sig80 = sigma_spectrum(sw_cross_block(s, 80, M=150001))
    local = np.log((1 - sig40.max()) / (1 - sig80.max())) / np.log(2.0)
    assert local > 1.93


def test_gerade_constant():
    """cond(I + C) -> 2/(1 + min j0) = 2.555041..., independent of s."""
    assert abs(gerade_constant() - 2.555041) < 1e-5
    for s in (1.4, 2.0):
        ev = np.linalg.eigvalsh(sw_cross_block(s, 120, M=150001))
        cond_g = (1.0 + ev.max()) / (1.0 + ev.min())
        assert abs(cond_g - 2.555041) < 0.01


def test_commutator_saturates():
    """Once the spectrum densifies the commutator pins at its 1/2 ceiling
    while cond diverges -- the norm is not a severity metric."""
    C = sw_cross_block(1.4, 20, M=200001)
    sig = sigma_spectrum(C)
    assert cond_from_sigma(sig) > 200.0
    assert abs(commutator_from_sigma(sig) - 0.5) < 0.02


def test_sw_matches_paper60_measured_matrix():
    """The tie to physics: the sine-basis symbol construction reproduces the
    paper's earlier triple-cross-validated measured cond(S) values (basis 6:
    30.1 / 15.4 / 4.4 at R = 1.4 / 2.0 / 4.0)."""
    for s_red, target in ((1.4, 30.1), (2.0, 15.4), (4.0, 4.4)):
        sig = sigma_spectrum(sw_cross_block(s_red, 3, M=300001))
        cond = cond_from_sigma(sig)
        assert abs(cond - target) / target < 0.02
