"""Genuine backing for Paper 40's universal 4/pi rate constant.

WHAT IS ROBUST (and tested here, DEFAULT-collected):
  * rank-1 (SU(2)): the rate constant is 4/pi exactly -- the doubling
    estimator of the closed-form sum-rule converges to 4/pi (production
    code geovac.central_fejer_su2, also the Paper 38 keystone).
  * general-G C_3 <= 1 (the L3 / Dirac-triangle keystone): for SU(3), Sp(2)
    and G2, |D(lambda) - D(lambda')| <= sqrt(C(sigma)) for every sigma in
    lambda (x) lambda'^* across a reduced dominant-weight panel (fail_count = 0,
    sup ratio < 1).  This is the general-rank content the universality rides on;
    it was previously verified only by an UNCALLED driver on a pruning-scheduled
    path (fixed v4.44.0 -- /qa group1 batch3 Flavor-B backfill).
  * rank-2 machinery correctness: the generic compact-Lie Weyl-integration
    used for SU(3)/Sp(2)/G2 integrates Haar to 1.0.
  * the A-over-B discrimination (G2, |W|=12): the extracted rate constant sits
    decisively closer to Reading A (c = 4/pi ~ 1.273) than Reading B
    (c = 2|W|/pi^r = 24/pi^2 ~ 2.432).  [@slow, full panel]

WHAT IS NOT a robust theorem (honest scope; calibrated in Paper 40 prose):
  the *exact* value c(G) = 4/pi at rank >= 2.  The leading-constant extraction
  is fit-sensitive -- 2-param vs 3-param Stein-Weiss fits and the min_L cut
  scatter the raw constant (Sp(2) 3-param across cuts: 3.12, 1.67, 0.73, 0.13,
  -0.64).  The clean table values (SU(3) 1.243, Sp(2) 1.087, G2 1.177) come from
  a specific 2-param fit + generic->canonical rescaling; a small live panel gives
  a larger raw c_can ~ 1.1-1.8.  The full rank-uniform analytical proof is a
  named gap; see Paper 40 Theorem (universality).

Backing drivers (permanent): tests/rank2_rate_support/ (moved out of the prunable
debug/qa/_resurrected/ on 2026-06-23, v4.44.0).
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
# rank-2 backing drivers (permanent test-support package)
# ---------------------------------------------------------------------------

_SUPPORT = os.path.join(os.path.dirname(__file__), "rank2_rate_support")
if _SUPPORT not in sys.path:
    sys.path.insert(0, os.path.abspath(_SUPPORT))


@pytest.fixture(scope="module")
def _rank2():
    # dirac-triangle (C_3 / L3) machinery + the rate (Weyl-integration) machinery
    from dirac_triangle_extended_verify import (
        build_A, build_C2, build_G2, panel_dominant_weights,
        run_panel as dirac_run_panel,
    )
    from sp2_g2_rate_constant import (
        run_panel as rate_run_panel, verify_haar_normalization,
    )
    return dict(build_A=build_A, build_C2=build_C2, build_G2=build_G2,
                panel_dominant_weights=panel_dominant_weights,
                dirac_run_panel=dirac_run_panel,
                rate_run_panel=rate_run_panel, verify_haar=verify_haar_normalization)


def _groups(r):
    """(label, algebra, rank) for the three rank-2 test groups."""
    return [("SU(3)", r["build_A"](2), 2),
            ("Sp(2)", r["build_C2"](), 2),
            ("G2", r["build_G2"](), 2)]


# ---------------------------------------------------------------------------
# general-G C_3 <= 1 (the L3 keystone) -- DEFAULT, fast reduced panel.
# Wires the previously-uncalled verify_dirac_triangle / run_panel.
# ---------------------------------------------------------------------------

def test_general_G_dirac_triangle_C3_leq_1(_rank2):
    """For SU(3), Sp(2), G2 the Dirac triangle inequality holds on a reduced
    dominant-weight panel: every ordered pair (lambda, lambda') satisfies
    |D(lambda) - D(lambda')| <= sqrt(C(sigma)) for all sigma in the tensor
    decomposition (fail_count = 0), with sup ratio strictly < 1.  This is the
    general-rank C_3 <= 1 / L3 content -- the keystone the universality rides on.
    """
    for label, la, _rank in _groups(_rank2):
        weights = _rank2["panel_dominant_weights"](la.rank, 2)  # sum a_i <= 2
        res = _rank2["dirac_run_panel"](la, weights, label, verbose_violations=False)
        assert res["fail_count"] == 0, (
            f"{label}: Dirac-triangle C_3<=1 violated in {res['fail_count']}/"
            f"{res['total_pairs']} pairs; max_ratio={res['max_ratio']:.4f}, "
            f"first violations={res['violations'][:3]}"
        )
        # non-vacuous: the panel actually exercises the bound (sup ratio in (0, 1])
        assert 0.0 < res["max_ratio"] <= 1.0 + 1e-9, (
            f"{label}: max_ratio={res['max_ratio']:.6f} not in (0, 1] "
            f"(vacuous panel or genuine C_3>1)"
        )
        assert res["total_pairs"] >= 9  # reduced panel still has real content


# ---------------------------------------------------------------------------
# rank-2 Weyl-integration correctness -- DEFAULT, fast (n_quad = 40).
# ---------------------------------------------------------------------------

def test_rank2_weyl_integration_haar_normalized(_rank2):
    """The rank-2 Weyl integration is correct: int_G 1 dg = 1 for SU(3), Sp(2),
    G2 (a real correctness check of the machinery the universality numbers ride
    on).  Fast at n_quad=40; the heavier n_quad=80 version is the slow test."""
    for label, la, _rank in _groups(_rank2):
        haar, _ = _rank2["verify_haar"](la, n_quad=40)
        assert abs(float(haar) - 1.0) < 1e-6, f"{label}: Haar = {float(haar)}"


# ---------------------------------------------------------------------------
# SLOW: comprehensive full-panel checks (the 977/5641-cell-class versions).
# ---------------------------------------------------------------------------

@pytest.mark.slow
def test_general_G_dirac_triangle_full_panel(_rank2):
    """Comprehensive general-G C_3<=1 on the VALIDATED group-specific
    Casimir-appropriate panels (SU(3) p+q<=5, Sp(2) a+b<=3, G2 a+b<=2 -- the
    bounds the original driver validated ALL-PASS; G2 irreps grow fast, so its
    panel is the smallest).  fail_count = 0, sup ratio < 1 on each.

    HONEST SCOPE: the all-sigma triangle is PANEL-BOUNDED.  Beyond these Casimir
    bounds it fails for extreme weight pairs (e.g. G2 (1,0) vs (0,4): a
    small-Casimir sigma gives ratio ~2.4 > 1) while the PRV / max-Casimir sigma
    still dominates.  Paper 40's C_3 = 1 claim is the *asymptotic* PRV-summand
    bound (the existence of a dominating sigma), of which these panels are
    empirical corroboration -- NOT a uniform all-weights all-sigma theorem.
    """
    panels = [("SU(3) p+q<=5", _rank2["build_A"](2), 2, 5),
              ("Sp(2) a+b<=3", _rank2["build_C2"](), 2, 3),
              ("G2 a+b<=2", _rank2["build_G2"](), 2, 2)]
    for label, la, rank, tot in panels:
        weights = _rank2["panel_dominant_weights"](rank, tot)
        res = _rank2["dirac_run_panel"](la, weights, label, verbose_violations=False)
        assert res["fail_count"] == 0, (
            f"{label}: {res['fail_count']}/{res['total_pairs']} fail; "
            f"max_ratio={res['max_ratio']:.4f}"
        )
        assert 0.0 < res["max_ratio"] <= 1.0 + 1e-9
        assert res["total_pairs"] >= 36


@pytest.mark.slow
def test_rank2_weyl_integration_haar_n_quad_80(_rank2):
    """Haar normalization at the heavier n_quad=80 quadrature."""
    for label, la, _rank in _groups(_rank2):
        haar, _ = _rank2["verify_haar"](la, n_quad=80)
        assert abs(float(haar) - 1.0) < 1e-6, f"{label}: Haar = {float(haar)}"


@pytest.mark.slow
def test_g2_rate_reading_A_over_B(_rank2):
    """G2 is the cleanest A-vs-B discriminator (|W|=12).  The extracted rate
    constant (canonical normalisation c_can = c_generic/sqrt(6)) is monotone
    decreasing toward Reading A (4/pi ~ 1.273) and sits decisively below
    Reading B (24/pi^2 ~ 2.432), i.e. closer to A than B.

    This tests the robust universality content (A over B), NOT the exact value
    4/pi (which is fit-sensitive -- see module docstring)."""
    g2 = _rank2["build_G2"]()
    rows, _ = _rank2["rate_run_panel"](g2, "G2", [48, 84, 144, 240, 420, 600])
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
