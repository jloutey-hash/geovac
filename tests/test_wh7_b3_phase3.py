# -*- coding: utf-8 -*-
"""Frozen falsifier for B3 Phase 3 Sprint 1: state-level cost structure on modular
orbits (2026-06-10, v3.115.0).

Pins: KMS/flow/period anchors; orbit injectivity; the D_max chain inequality
(0 violations on the 96-cell panel — the universal twin-direction fact); the
primary-configuration excess positivity + bimodal scaling split; AND the honest
caveat (excess sign is reference-state dependent, documented by T5).

Driver: debug/wh7_b3_phase3_state_intervals.py (deterministic seeds; pinned by
this reference).
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "debug"))
import wh7_b3_phase3_state_intervals as p3  # noqa: E402


@pytest.fixture(scope="module")
def res():
    return p3.run()


def test_kms_flow_period_anchors(res):
    assert res["T0"]["kms_residual"] < 1e-12
    assert res["T0"]["flow_invariance"] < 1e-12
    assert res["T0"]["orbit_period_2pi_residual"] < 1e-12


def test_orbit_injective(res):
    assert res["T1"]["min_orbit_distance"] > 0.05


def test_dmax_chain_inequality_universal(res):
    assert res["T2"]["n_cells"] == 96
    assert res["T2"]["dmax_chain_violations"] == 0


def test_primary_config_excess_positive(res):
    for name, entry in res["T3_T4_class_table"].items():
        for tag in ("raw", "perp"):
            for row in entry[tag]["rows"]:
                assert row["excess_dmax"] > -1e-9, (name, tag, row)


def test_bimodal_scaling_split_primary_config(res):
    tab = res["T3_T4_class_table"]
    for name in ("(0.5,0.5) spacelike", "(1.0,0.0) spacelike", "(2.0,0.0) spacelike"):
        assert tab[name]["raw"]["fit_power"] > 1.6, name
    for name in ("(2.0,1.0) spacelike", "(1.0,1.0) null",
                 "(1.5,1.5) timelike", "(2.0,2.0) timelike"):
        assert 0.9 < tab[name]["raw"]["fit_power"] < 1.5, name


def test_tangent_confound_ruled_out(res):
    """Projection of the flow-tangent component changes the powers by < 0.1 —
    the bimodal split is NOT a reparametrization artifact."""
    for name, entry in res["T3_T4_class_table"].items():
        assert abs(entry["raw"]["fit_power"] - entry["perp"]["fit_power"]) < 0.1


def test_state_dependence_caveat_documented(res):
    """T5: at least one (seed, class) cell has nonpositive excess (p ~ 0 marker).
    Freezing the CAVEAT: excess-over-baseline positivity is reference-state
    dependent; only the absolute chain inequality is universal."""
    flat = [p for ps in res["T5_robustness"].values() for p in ps]
    assert any(p < 0.5 for p in flat)
    assert any(p > 1.6 for p in flat)
