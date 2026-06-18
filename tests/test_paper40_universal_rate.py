"""Genuine backing for Paper 40's universal 4/pi rate constant.

WHAT IS ROBUST (and tested here):
  * rank-1 (SU(2)): the rate constant is 4/pi exactly -- the doubling
    estimator of the closed-form sum-rule converges to 4/pi (production
    code geovac.central_fejer_su2, also the Paper 38 keystone).
  * rank-2 machinery correctness: the generic compact-Lie Weyl-integration
    used for SU(3)/Sp(2)/G2 integrates Haar to 1.0.
  * the A-over-B discrimination: at G2 (the cleanest discriminator, |W|=12),
    the extracted rate constant sits decisively closer to Reading A
    (c = 4/pi ~ 1.273) than to Reading B (c = 2|W|/pi^r = 24/pi^2 ~ 2.432).

WHAT IS NOT a robust theorem (honest scope; calibrated in Paper 40 prose):
  the *exact* value c(G) = 4/pi at rank >= 2.  The leading-constant
  extraction is fit-sensitive -- 2-param vs 3-param Stein-Weiss fits and the
  min_L cut scatter the raw constant (Sp(2) 3-param across cuts: 3.12, 1.67,
  0.73, 0.13, -0.64).  The clean table values (SU(3) 1.243, Sp(2) 1.087, G2
  1.177) come from a specific 2-param fit + generic->canonical rescaling.
  The full rank-uniform analytical proof is a named gap; see Paper 40
  Theorem (universality) and debug/qa/_resurrected/ for the recovered
  computation + provenance.
"""
import math
import os
import sys

import pytest


# ---------------------------------------------------------------------------
# rank-1 anchor: SU(2) rate constant is 4/pi exactly (production code)
# ---------------------------------------------------------------------------

def test_su2_rate_constant_is_4_over_pi():
    """SU(2) doubling estimator a_n -> 4/pi (the rigorous rank-1 case)."""
    from geovac.central_fejer_su2 import doubling_estimator
    import mpmath as mp
    target = 4 / mp.pi
    a_small = doubling_estimator(64)
    a_large = doubling_estimator(512)
    # converges to 4/pi from above, monotone improving
    assert a_large < a_small
    assert abs(a_large - target) < abs(a_small - target)
    assert abs(a_large - target) < 0.05  # within 5% by n=512 (slow log convergence)


# ---------------------------------------------------------------------------
# rank-2: recovered genuine Weyl-integration machinery
# ---------------------------------------------------------------------------

_RES = os.path.join(os.path.dirname(__file__), "..", "debug", "qa", "_resurrected")
_HAVE_RES = os.path.isfile(os.path.join(_RES, "sp2_g2_rate_constant.py"))
_skip = pytest.mark.skipif(not _HAVE_RES,
                           reason="recovered rank-2 rate drivers not present")


@pytest.fixture(scope="module")
def _rank2():
    sys.path.insert(0, os.path.abspath(_RES))
    from dirac_triangle_extended_verify import build_A, build_C2, build_G2
    from sp2_g2_rate_constant import run_panel, verify_haar_normalization
    return dict(build_A=build_A, build_C2=build_C2, build_G2=build_G2,
                run_panel=run_panel, verify_haar=verify_haar_normalization)


@_skip
@pytest.mark.slow
def test_rank2_weyl_integration_haar_normalized(_rank2):
    """The rank-2 Weyl integration is correct: int_G 1 dg = 1 for SU(3),
    Sp(2), G2 (a real correctness check of the machinery the universality
    numbers ride on)."""
    for build in (_rank2["build_A"], _rank2["build_C2"], _rank2["build_G2"]):
        g = build(2) if build is _rank2["build_A"] else build()
        haar, _ = _rank2["verify_haar"](g, n_quad=80)
        assert abs(haar - 1.0) < 1e-6


@_skip
@pytest.mark.slow
def test_g2_rate_reading_A_over_B(_rank2):
    """G2 is the cleanest A-vs-B discriminator (|W|=12).  The extracted rate
    constant (canonical normalisation c_can = c_generic/sqrt(6)) is monotone
    decreasing toward Reading A (4/pi ~ 1.273) and sits decisively below
    Reading B (24/pi^2 ~ 2.432), i.e. closer to A than B.

    This tests the robust universality content (A over B), NOT the exact value
    4/pi (which is fit-sensitive -- see module docstring)."""
    g2 = _rank2["build_G2"]()
    rows, _ = _rank2["run_panel"](g2, "G2", [48, 84, 144, 240, 420, 600])
    sqrt6 = math.sqrt(6.0)
    c_can = [r["L*gamma/log(L)"] / sqrt6 for r in rows]
    A = 4 / math.pi          # 1.273
    B = 24 / math.pi ** 2    # 2.432
    # rate vanishes (gamma decreasing)
    gammas = [r["gamma"] for r in rows]
    assert gammas[-1] < gammas[0]
    # tail constant decreasing toward A, and decisively below B
    assert c_can[-1] < c_can[2]
    assert c_can[-1] < 2.0 < B
    # closer to A than to B
    assert abs(c_can[-1] - A) < abs(c_can[-1] - B)
