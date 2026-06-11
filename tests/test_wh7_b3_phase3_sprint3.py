# -*- coding: utf-8 -*-
"""Frozen falsifier for B3 Phase 3 Sprint 3: cone-graded admissibility, cutoff
stability, Bures refutation (2026-06-10).

Pins:
  (1) the wedge folding REORGANIZES the Phase-2 causal classes: the (2,1)
      spacelike class is annihilated (folds to zero), and the (2,2) timelike
      top class becomes flow-commuting on the wedge (m' -> -m' identification
      turns maximal weight transfer into zero absolute-weight transfer) --
      admissibility is the flow-commutation grading, NOT the causal grading;
  (2) penalty magnitudes organize by K-weight transfer (Spearman r > 0.9)
      at least as well as by the causal ratio q_F;
  (3) the band-limited interval-penalty structure is bit-exactly
      cutoff-independent across n_max = 2..5 (Gibbs normalization cancels in
      the D_max ratio; band-limited conjugations leave the bulk untouched),
      while naive trace-norm orbit scale decays exactly as const/Z (the
      finite-cutoff shadow of type III);
  (4) PERIOD ATTRIBUTION CORRECTED: U_pi = -1 already on the FULL Dirac
      space at every n_max (all two_m_j odd) -- the period-pi is the spinor
      spin-statistics grading, not a wedge-folding effect (corrects the
      Sprint-2 attribution);
  (5) Bures excess positivity on flow-commuting kicks is REFUTED (574/2400
      adversarial cells negative; the Sprint-2 class generators themselves go
      negative on the widened panel) -- NO tested functional has a
      sign-definite even-sector penalty; only the chain inequality is
      universal.

Driver: debug/wh7_b3_phase3_sprint3.py (deterministic seeds; pinned by this
reference).
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "debug"))
import wh7_b3_phase3_sprint3 as p33  # noqa: E402


@pytest.fixture(scope="module")
def res():
    return p33.run()


# ---------- T1: folding reorganizes the cone classes ---------------------------
def test_wedge_construction(res):
    t1 = res["T1_cone_admissibility"]
    assert t1["R_involution_residual"] < 1e-12
    assert t1["dim_wedge"] == 9
    assert t1["weights"] == [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 2.0]


def test_spacelike_21_class_annihilated(res):
    fi = res["T1_cone_admissibility"]["fold_info"]["(2.0,1.0) spacelike"]
    assert fi["fold_norm"] < 1e-12


def test_timelike_top_becomes_admissible_on_wedge(res):
    """Upstairs (2,2) has maximal transfer 4 (linear penalty, median |D+| =
    0.26); on the wedge it commutes with K_W and is even at bit level."""
    t1 = res["T1_cone_admissibility"]
    fi = t1["fold_info"]["(2.0,2.0) timelike"]
    assert fi["fold_norm"] > 0.9
    assert fi["K_commutator_wedge"] < 1e-12
    ew = t1["ensemble_wedge"]["(2.0,2.0) timelike"]
    assert ew["median_abs_Dp"] < 1e-2
    assert ew["median_even"] < 1e-12
    assert t1["ensemble_full"]["(2.0,2.0) timelike"]["median_abs_Dp"] > 0.1


def test_m0_classes_stay_admissible_on_wedge(res):
    t1 = res["T1_cone_admissibility"]
    for name in ("(1.0,0.0) spacelike", "(2.0,0.0) spacelike"):
        assert t1["fold_info"][name]["K_commutator_wedge"] < 1e-12, name
        assert t1["ensemble_wedge"][name]["median_even"] < 1e-12, name


def test_transfer_classes_stay_linear_on_wedge(res):
    t1 = res["T1_cone_admissibility"]
    for name in ("(0.5,0.5) spacelike", "(1.5,1.5) timelike",
                 "(1.0,1.0) null"):
        assert t1["fold_info"][name]["K_commutator_wedge"] > 0.5, name
        assert t1["ensemble_wedge"][name]["median_abs_Dp"] > 1e-2, name


def test_transfer_organizes_at_least_as_well_as_cone(res):
    org = res["T1_cone_admissibility"]["organization_full"]
    assert org["transfer"]["spearman_r"] > 0.9
    assert org["transfer"]["spearman_r"] >= org["qF"]["spearman_r"] - 1e-9


# ---------- T2: bit-exact band stability ----------------------------------------
def test_period_pi_holds_upstairs_at_every_cutoff(res):
    """Correction of the Sprint-2 attribution: the full-Dirac U_pi = -1 at
    every n_max (odd two_m_j spectrum); period-pi is NOT a folding effect."""
    for n, r in res["T2_cutoff_stability"]["per_n"].items():
        assert r["U_pi_plus_identity_full"] < 1e-12, n


def test_band_limited_penalties_cutoff_independent(res):
    per_n = res["T2_cutoff_stability"]["per_n"]
    for n in ("2", "3", "4"):
        r = per_n[n]
        assert r["delta_E02_dw0_vs_ref"] < 1e-9, n
        assert r["delta_E02_dw2_vs_ref"] < 1e-9, n
        assert r["delta_Dp_dw2_vs_ref"] < 1e-7, n


def test_evenness_cutoff_uniform(res):
    for n, r in res["T2_cutoff_stability"]["per_n"].items():
        assert r["evenness_dw0"] < 1e-12, n


def test_trace_orbit_scale_decays_as_inverse_Z(res):
    """orbit_scale * Z is constant across cutoffs (observed 0.1805...):
    naive trace objects vanish in the limit; only band-relative (D_max)
    objects survive -- the finite-cutoff shadow of type III."""
    per_n = res["T2_cutoff_stability"]["per_n"]
    prods = [per_n[n]["orbit_scale_trace"] * per_n[n]["Z"]
             for n in ("2", "3", "4", "5")]
    assert max(prods) - min(prods) < 1e-6
    # and the scale itself decays
    scales = [per_n[n]["orbit_scale_trace"] for n in ("2", "3", "4", "5")]
    assert scales == sorted(scales, reverse=True)


def test_dmax_interval_recovery_cutoff_stable(res):
    for n, r in res["T2_cutoff_stability"]["per_n"].items():
        assert r["interval_recovery_err_dmax"] < 1e-3, n


# ---------- T3: Bures positivity refuted -----------------------------------------
def test_bures_positivity_refuted(res):
    t3 = res["T3_bures_adversarial"]
    assert t3["n_cells"] == 2400
    assert t3["max_K_commutator"] < 1e-12        # kicks genuinely admissible
    assert t3["bures_negative"] > 100            # observed 574
    assert t3["worst_bures_excess"] < -1e-3


def test_sprint2_class_generators_go_negative_on_wider_panel(res):
    """The Sprint-2 'Bures positive on low-transfer classes' observation was
    a narrow-panel artifact (theta <= 0.6, t_total = 1.0, eps = 0.2): the
    same generators go negative once theta and t_total widen."""
    t3 = res["T3_bures_adversarial"]
    assert t3["class_gen_negative"] > 10         # observed 134/600


def test_no_functional_sign_definite_even_sector(res):
    """D_max shows the same indefiniteness on the admissible sector."""
    assert res["T3_bures_adversarial"]["dmax_negative"] > 100
