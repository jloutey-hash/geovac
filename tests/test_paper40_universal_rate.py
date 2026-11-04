"""Backing for Paper 40's leading Fejer-moment constant (a kernel constant, not a convergence rate).

WHAT IS TESTED (DEFAULT-collected):
  * SU(2), Paper 40 Lemma L2(d): gamma_Lambda(SU(2)) = (2/pi) log Lambda / Lambda
    + O(1/Lambda) in Paper 40's own variable Lambda (Cas(j) = j(j+1) <= Lambda^2,
    so Lambda = sqrt(j(j+1)) ~ n/2 for the shell n = 2j + 1).  This is the same
    function as Paper 38's (4/pi) log n / n in the shell index n; both are pinned
    from production code (geovac.central_fejer_su2), and the per-Lambda slope is
    measured directly against Lambda, so the retracted reading "4/pi per Lambda"
    fails.
  * the general-G Dirac-triangle inequality (DT): for SU(3), Sp(2) and G2,
    |D(lambda) - D(lambda')| <= sqrt(C(sigma)) for every sigma in
    lambda (x) lambda'^* across a reduced dominant-weight panel (fail_count = 0,
    sup ratio < 1).  (DT) is a representation-theoretic inequality; it is NOT
    a check of the Lipschitz constant C_3 (= 1 by compression, Paper 40 Lemma
    L3) and no statement about the leading constant rides on it.
  * rank-2 machinery correctness: the generic compact-Lie Weyl-integration
    used for SU(3)/Sp(2)/G2 integrates Haar to 1.0.

WHAT IS NOT TESTED, and why: whether the leading constant is the same for every
compact simple G is [OPEN] (Paper 40, Remark thm:universal_constant).  The
rank-2 extraction pipeline divides by sqrt(C_gen(ad)/h^vee) while leaving the
distances unscaled; on SU(2) that convention gives (2/pi)/sqrt2 ~ 0.450, not the
anchor 2/pi, so its rank-2 outputs (SU(3) 1.243, Sp(2) 1.087, G2 1.177) are in a
mismatched convention and discriminate neither for nor against a
group-independent constant.  The former G2 "Reading A over Reading B" assertion
was withdrawn with that verdict (2026-11-03); no rank-2 constant is asserted here
until a consistent-convention recomputation exists.

Backing drivers (permanent): tests/rank2_rate_support/.
"""
import math
import os
import sys

import pytest


# ---------------------------------------------------------------------------
# rank-1 anchor: SU(2), production code
# ---------------------------------------------------------------------------

def _lam(n: int) -> float:
    """Paper 40's Lambda for the SU(2) shell n = 2j + 1: Lambda^2 = j(j + 1)."""
    return math.sqrt((n * n - 1) / 4.0)


def test_su2_constant_is_2_over_pi_per_lambda():
    """Lemma L2(d): the slope of Lambda * gamma against log Lambda is 2/pi.

    Rejects the retracted "4/pi per Lambda" (off by a factor 2) and any
    normalisation error above ~1.5% (e.g. 3.9/pi-type scale slips)."""
    from geovac.central_fejer_su2 import gamma_n_via_sum_rule
    n1, n2 = 512, 1024
    g1, g2 = float(gamma_n_via_sum_rule(n1)), float(gamma_n_via_sum_rule(n2))
    slope = (_lam(n2) * g2 - _lam(n1) * g1) / (math.log(_lam(n2)) - math.log(_lam(n1)))
    assert abs(slope / (2 / math.pi) - 1) < 0.015, f"per-Lambda slope {slope:.5f}"


def test_su2_constant_is_4_over_pi_per_shell_index():
    """Paper 38's form of the same function: the doubling estimator in the
    shell index n converges to 4/pi (monotone, within 1.5% at n = 512)."""
    from geovac.central_fejer_su2 import doubling_estimator
    a_small = float(doubling_estimator(64))
    a_large = float(doubling_estimator(512))
    target = 4 / math.pi
    assert abs(a_large - target) < abs(a_small - target)
    assert abs(a_large / target - 1) < 0.015, f"doubling estimator {a_large:.5f}"


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
# general-G Dirac-triangle inequality (DT) -- DEFAULT, fast reduced panel.
# (DT) is a representation-theoretic inequality on Kostant eigenvalues and
# Casimirs; it is NOT a check of the Lipschitz constant C_3, which is 1 by
# compression (Paper 40 Lemma L3; tests/test_paper40_l3_compression.py).  The
# bridge from (DT) to C_3 was withdrawn 2026-10-02.
# Wires the previously-uncalled verify_dirac_triangle / run_panel.
# ---------------------------------------------------------------------------

def test_general_G_dirac_triangle_inequality_panel(_rank2):
    """For SU(3), Sp(2), G2 the Dirac triangle inequality holds on a reduced
    dominant-weight panel: every ordered pair (lambda, lambda') satisfies
    |D(lambda) - D(lambda')| <= sqrt(C(sigma)) for all sigma in the tensor
    decomposition (fail_count = 0), with sup ratio strictly < 1.  A statement
    about eigenvalues, not about the commutator norm (the bridge is withdrawn).
    """
    for label, la, _rank in _groups(_rank2):
        weights = _rank2["panel_dominant_weights"](la.rank, 2)  # sum a_i <= 2
        res = _rank2["dirac_run_panel"](la, weights, label, verbose_violations=False)
        assert res["fail_count"] == 0, (
            f"{label}: Dirac-triangle inequality (DT) violated in {res['fail_count']}/"
            f"{res['total_pairs']} pairs; max_ratio={res['max_ratio']:.4f}, "
            f"first violations={res['violations'][:3]}"
        )
        # non-vacuous: the panel actually exercises the bound (sup ratio in (0, 1])
        assert 0.0 < res["max_ratio"] <= 1.0 + 1e-9, (
            f"{label}: max_ratio={res['max_ratio']:.6f} not in (0, 1] "
            f"(vacuous panel or a (DT) violation)"
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
    """Comprehensive general-G Dirac-triangle inequality (DT; not a C_3 check) on the VALIDATED group-specific
    Casimir-appropriate panels (SU(3) p+q<=5, Sp(2) a+b<=3, G2 a+b<=2 -- the
    bounds the original driver validated ALL-PASS; G2 irreps grow fast, so its
    panel is the smallest).  fail_count = 0, sup ratio < 1 on each.

    HONEST SCOPE (corrected 2026-09-15): the all-sigma triangle is
    PANEL-BOUNDED by the decomposition driver's validated range.  The panels
    below (dimension-clean) show fail_count == 0.  An EARLIER version of this
    note cited G2 (1,0) vs (0,4) "ratio ~2.4 > 1" as a real DT counterexample
    beyond the panel;  that is WITHDRAWN -- the 2.4 is an artifact of
    dirac_triangle_extended_verify.py's tensor_product, which is not
    dimension-conserving outside its validated panel (it reports a
    Schur-impossible trivial summand for that pair).  The dimension-correct DT
    value there is < 1 (holds).  What the code establishes is therefore:
    all-sigma DT holds on every dimension-clean panel, plus the asymptotic
    PRV-summand bound (a dominating sigma).  The paper's Cor L3_closure claims
    the interior closure at all ranks via an ANALYTICAL Steinberg argument that
    the code does not verify;  that tier is a primary-math question, not settled
    here.
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
