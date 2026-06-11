# -*- coding: utf-8 -*-
"""Frozen falsifier for B3 Phase 1: boost/modular-flow seminorm structure
(2026-06-10, v3.113.0). Pins the Paper 45 Q1 boost-candidate probe results:
system rank 55; sigma_{2pi} closure + half-period band-parity grading bit-exact;
boost-alone kernel dim 9 (boost-invariant multipliers); frame kernel dim 1;
causal-classifier sign pattern with the b=1 top weight on the cone.

Driver: debug/wh7_b3_boost_seminorm_probe.py (pinned by this reference).
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "debug"))
import wh7_b3_boost_seminorm_probe as b3  # noqa: E402


@pytest.fixture(scope="module")
def res():
    return b3.run()


def test_system_rank_55(res):
    assert res["T0_rank"]["rank"] == 55


def test_sigma_2pi_closure_bit_exact(res):
    assert res["T1_flow"]["sigma_2pi_residual"] < 1e-12


def test_half_period_band_parity_grading(res):
    assert res["T1_flow"]["band_parity_grading_residual"] < 1e-12


def test_boost_alone_kernel_dim_9(res):
    assert res["T2_boost_kernel"]["dim_ker"] == 9
    assert res["T2_boost_kernel"]["system_invariance_residual"] < 1e-10


def test_frame_kernel_dim_1(res):
    assert res["T3_frame_kernel"]["dim_ker"] == 1


def test_causal_classifier_sign_pattern(res):
    rows = {(r["b"], r["abs_mp"]): r for r in res["T4_causal_classes"]
            if r["norm"] == "op"}
    # purely spacelike classes (q < 0 strictly)
    for key in [(0.5, 0.5), (1.0, 0.0), (1.5, 0.5), (2.0, 0.0), (2.0, 1.0)]:
        assert rows[key]["q_max"] < -1e-6, key
    # timelike top weights from b = 3/2 up
    for key in [(1.5, 1.5), (2.0, 2.0)]:
        assert rows[key]["q_min"] > 1e-6, key
    # b = 1 top weight straddles the cone with its minimum ON the null ray
    r = rows[(1.0, 1.0)]
    assert abs(r["q_min"]) < 1e-6 and r["q_max"] > 0.1
