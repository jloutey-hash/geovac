"""Backing for Paper 60 sec:molecular's [MEASURED] "one direction" paragraph
(added 2026-09-12).

Why this test exists.  A frames-theoretic reading of the two-center degeneracy
was adopted on 2026-09-11 and is WRONG for this basis: it held that
completeness of the ONE-CENTRE set forces lam_min -> 0 (if g lies in the closed
span of {f_i} then lam_min(G_N) <= dist(g, V_N)^2 -> 0).  That hypothesis is
measurable here, because the intra-centre block of the SW metric is exactly the
identity, so the A-set is orthonormal in that metric and Bessel's inequality
gives the deficit directly:

    eps_N^2 = 1 - sum_{i<=N} <chi^A_i, chi^B_1>^2  >= 0.

It does not go to zero -- it plateaus around 0.38 / 0.70 / 0.91 for kR = 1/2/4.
So the one-centre set is far from complete in the molecular metric, the
hypothesis fails, and the degeneracy is NOT "one span already contains the
other".  What sigma_max -> 1 says is narrower: some COMBINATION is captured.
The degeneracy is one direction, not a diffuse property of the basis.

The claim is therefore a PAIR, and is tested as a pair: one quantity collapses
(1 - sigma_max ~ N^-2) while the other stays put (eps_N^2 plateaus).  A guard
asserting only one half could be satisfied by a bug that made everything
collapse, or everything plateau.
"""
import numpy as np
import pytest

from geovac.sturmian_sigma_law import sw_cross_block

M_QUAD = 200_001
KR_VALUES = (1.0, 2.0, 4.0)
# measured plateau of eps_N^2, flat over N = 16..256
EXPECTED_PLATEAU = {1.0: 0.380, 2.0: 0.696, 4.0: 0.907}


def bessel_deficit(s, N):
    """1 - sum_i <chi^A_i, chi^B_1>^2 in the SW metric (A-set orthonormal there)."""
    C = sw_cross_block(s, N, M=M_QUAD)
    return float(1.0 - np.sum(C[:, 0] ** 2))


@pytest.mark.parametrize("kR", KR_VALUES)
def test_one_centre_set_is_not_complete_in_the_molecular_metric(kR):
    """eps_N^2 plateaus strictly above zero -- the Prop-A hypothesis FAILS here.

    WRONG ANSWER REJECTED: eps_N -> 0, i.e. the reading that the one-centre span
    already contains a displaced Sturmian.  Asserted three ways a decaying
    sequence cannot satisfy: the value at N=256 must exceed 0.3, the ratio
    between N=32 and N=256 must be within 1% of unity (a plateau, not a decay),
    and a fitted log-log slope must be flatter than 0.01 in magnitude.
    """
    vals = {N: bessel_deficit(kR, N) for N in (32, 64, 128, 256)}
    for N, v in vals.items():
        assert 0.0 < v < 1.0, f"deficit out of Bessel range at N={N}: {v}"

    assert vals[256] > 0.3, f"deficit collapsed to {vals[256]:.4f}"
    assert abs(vals[256] / vals[32] - 1.0) < 0.01, f"not a plateau: {vals}"
    slope = np.polyfit(np.log(list(vals)), np.log(list(vals.values())), 1)[0]
    assert abs(slope) < 0.01, f"deficit is decaying, slope {slope:+.4f}"
    assert abs(vals[256] - EXPECTED_PLATEAU[kR]) < 0.005, (
        f"plateau {vals[256]:.4f} != recorded {EXPECTED_PLATEAU[kR]}")


def test_the_degeneracy_is_one_direction_not_the_whole_basis():
    """The PAIRED claim: 1 - sigma_max collapses while eps_N^2 does not.

    This is the content of the paragraph and the reason the rank-(M-1) rotation
    works.  Both halves are asserted together.

    WRONG ANSWER REJECTED: either collapse of the pair.  If the basis really did
    become globally redundant, eps_N^2 would fall too (first assertion); if
    there were no near-dependence at all, 1 - sigma_max would not fall like
    N^-2 (second).  A single-sided guard admits one of those; this one admits
    neither.
    """
    kR = 2.0
    Ns = (32, 64, 128, 256)
    deficits, gaps = [], []
    for N in Ns:
        C = sw_cross_block(kR, N, M=M_QUAD)
        deficits.append(float(1.0 - np.sum(C[:, 0] ** 2)))
        gaps.append(float(1.0 - np.linalg.svd(C, compute_uv=False)[0]))

    # half 1: the deficit does NOT collapse
    assert min(deficits) > 0.3, f"deficits collapsed: {deficits}"
    assert max(deficits) - min(deficits) < 0.01, f"deficits not flat: {deficits}"

    # half 2: the spectral gap DOES collapse, quadratically
    for a, b in zip(gaps, gaps[1:]):
        assert 3.0 < a / b < 5.0, f"gap not falling as N^-2: {gaps}"
    assert gaps[-1] < gaps[0] / 50, f"gap barely moved: {gaps}"

    # and therefore: orders of magnitude apart at the largest N
    assert deficits[-1] / gaps[-1] > 1e3, (
        f"deficit {deficits[-1]:.3e} is not vastly larger than the gap {gaps[-1]:.3e}")
