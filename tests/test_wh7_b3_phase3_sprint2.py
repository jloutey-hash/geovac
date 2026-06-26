# -*- coding: utf-8 -*-
"""Frozen falsifier for B3 Phase 3 Sprint 2: evenness mechanism, cost-functional
comparison, and the wedge substrate (2026-06-10).

Pins:
  (1) the EVENNESS MECHANISM -- for kicks commuting with the boost ([G,K] = 0,
      zero K-weight transfer), the kicked legs satisfy c12(eps) = c23(-eps)
      exactly (bit-level), forcing quadratic leading order; for weight-
      transferring kicks the identity breaks at O(eps) and a signed linear
      term appears (this REPLACES the Sprint-1 "top-eigenspace weight
      selection rule" suspicion and EXPLAINS the Sprint-1 sign caveat);
  (2) the per-leg slopes do NOT individually vanish for the commuting classes
      -- the cancellation is between legs, not within them;
  (3) Umegaki/Jeffreys chain inequality fails GENERICALLY on modular-orbit
      triples (Sprint-1's 0/96 was reference-state luck); D_max/Bures/trace
      hold by theorem (0 violations);
  (4) no cost functional gives state-independent positive excess on the
      ensemble (the distributional closure of the Sprint-1 caveat);
  (5) wedge substrate (HemisphericWedge, unfolded K_W, odd positive integer
      spectrum): anchors bit-level, orbit period = pi (U_pi = -I), D_max chain
      0 violations, the commuting/non-commuting dichotomy transfers, and the
      operational interval functional (flow parameter recovered from state
      pairs by trace-distance matching) is well-defined and additive.

Driver: tests/wh7_support/wh7_b3_phase3_sprint2.py (deterministic seeds; pinned by this
reference).
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent / "wh7_support"))
import wh7_b3_phase3_sprint2 as p32  # noqa: E402

COMMUTING = ("(1.0,0.0) spacelike", "(2.0,0.0) spacelike")
NONCOMMUTING = ("(0.5,0.5) spacelike", "(2.0,1.0) spacelike",
                "(1.0,1.0) null", "(1.5,1.5) timelike", "(2.0,2.0) timelike")


@pytest.fixture(scope="module")
def res():
    return p32.run()


# ---------- T3: the evenness mechanism ----------------------------------------
def test_base_legs_equal_by_flow_translation(res):
    assert res["T3_mechanism"]["base_leg_equality"] < 1e-10


def test_evenness_identity_commuting_classes(res):
    tab = res["T3_mechanism"]["classes"]
    for name in COMMUTING:
        e = tab[name]
        assert e["K_commutator_norm"] < 1e-12, name
        assert e["leg_evenness_residual"] < 1e-10, name
        assert e["kind"] == "quadratic", name
        assert abs(e["D_plus"]) < 1e-2 and abs(e["D_minus"]) < 1e-2, name


def test_evenness_broken_noncommuting_classes(res):
    tab = res["T3_mechanism"]["classes"]
    for name in NONCOMMUTING:
        e = tab[name]
        assert e["K_commutator_norm"] > 0.5, name
        assert e["leg_evenness_residual"] > 1e-3, name
        assert e["kind"] == "smooth_linear", name
        assert abs(e["D_plus"]) > 1e-2, name


def test_commutator_norm_equals_weight_transfer(res):
    """[K, G] norm = 2|m'| exactly: the weight-transfer witness."""
    tab = res["T3_mechanism"]["classes"]
    expected = {"(0.5,0.5) spacelike": 1.0, "(2.0,1.0) spacelike": 2.0,
                "(1.0,1.0) null": 2.0, "(1.5,1.5) timelike": 3.0,
                "(2.0,2.0) timelike": 4.0}
    for name, val in expected.items():
        assert abs(tab[name]["K_commutator_norm"] - val) < 1e-10, name


def test_per_leg_slopes_cancel_not_vanish(res):
    """The selection mechanism is leg cancellation, NOT per-leg vanishing:
    the rho-slot first-order term is nonzero for the commuting classes even
    though the total derivative vanishes."""
    tab = res["T3_mechanism"]["classes"]
    for name in COMMUTING:
        e = tab[name]
        assert abs(e["analytic_rho_slot_slope"]) > 1e-2, name
        assert abs(e["analytic_rho_slot_slope"]
                   - e["numerical_rho_slot_slope"]) < 1e-2, name


# ---------- T2: cost-functional chain comparison -------------------------------
def test_chain_theorem_functionals_hold(res):
    for k in ("dmax", "bures", "trace"):
        assert res["T2_chain_panel"][k]["chain_violations"] == 0, k
        assert res["T2_chain_panel"][k]["n_cells"] == 96, k


def test_umegaki_chain_fails_generically(res):
    """Sprint-1 saw 0/96 Umegaki violations on its panel; at this reference
    state the failure is generic (negative baseline deficit) -- the Paper 49
    Datta-replacement lesson reproduced on modular-orbit triples."""
    assert res["T2_chain_panel"]["umegaki"]["chain_violations"] > 50
    assert res["T2_chain_panel"]["umegaki"]["baseline_deficit"] < 0
    assert res["T2_chain_panel"]["jeffreys"]["chain_violations"] > 50


# ---------- T1: ensemble sign statistics ---------------------------------------
def test_no_state_independent_positive_excess(res):
    for k, sm in res["T1_ensemble"].items():
        assert sm["frac_positive_overall"] < 0.999, k


def test_bures_positive_on_low_transfer_classes(res):
    """Observed (panel-wide, not a theorem): the Bures excess is positive on
    every ensemble cell for the three low-weight-transfer classes."""
    pc = res["T1_ensemble"]["bures"]["per_class"]
    for name in ("(0.5,0.5) spacelike", "(1.0,0.0) spacelike",
                 "(2.0,0.0) spacelike"):
        assert pc[name]["frac_positive"] > 0.99, name


def test_dmax_least_sign_stable(res):
    """D_max (max-eigenvalue non-smoothness) is the least sign-stable cost
    on the ensemble, while its chain inequality is the theorem anchor."""
    sm = res["T1_ensemble"]
    for k in ("umegaki", "bures", "trace", "jeffreys"):
        assert (sm["dmax"]["frac_positive_overall"]
                < sm[k]["frac_positive_overall"]), k


# ---------- T4: wedge substrate -------------------------------------------------
def test_wedge_anchors_bit_level(res):
    w = res["T4_wedge"]
    assert w["dim_W"] * 2 == w["dim_full"]
    assert w["kms_residual"] < 1e-12
    assert w["flow_invariance"] < 1e-12
    assert w["orbit_period_pi_residual"] < 1e-12
    # period pi: Sprint-3 correction — this is the spinor (all-odd two_m_j)
    # spin-statistics grading, holding upstairs too, NOT a folding effect
    assert w["U_pi_is_minus_identity"] < 1e-12
    assert w["KW_spectrum_min"] >= 1.0           # genuine Boltzmann state


def test_wedge_chain_universal(res):
    assert res["T4_wedge"]["chain"]["dmax_chain_violations"] == 0


def test_wedge_transfer_dichotomy(res):
    """The commuting/non-commuting dichotomy transfers to the proper
    substrate: Dw = 0 kicks are even (bit-level) with quadratic scaling;
    Dw >= 2 kicks break evenness at O(eps) with a nonzero linear term."""
    tab = res["T4_wedge"]["transfer_table"]
    for r in tab["0"]:
        assert r["K_commutator_norm"] < 1e-12
        assert r["leg_evenness_residual"] < 1e-10
        assert abs(r["D_plus"]) < 1e-2
        assert 1.6 < r["fit_power"] < 2.4
    for dw in ("2", "4"):
        for r in tab[dw]:
            assert r["K_commutator_norm"] > 1.0
            assert r["leg_evenness_residual"] > 1e-8
            assert abs(r["D_plus"]) > 1e-3


def test_wedge_interval_functional(res):
    """ell = flow parameter is an operational STATE-PAIR functional on
    orbits: recovered by trace-distance matching (well-defined mod pi)
    and additive on ordered triples."""
    w = res["T4_wedge"]
    assert w["interval_recovery_max_err"] < 5e-3
    assert w["interval_additivity_max_err"] < 5e-3
    assert w["min_orbit_distance"] > 0.01
