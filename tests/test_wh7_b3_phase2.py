# -*- coding: utf-8 -*-
"""Frozen falsifier for B3 Phase 2: exact cone ratios + causal-form signature
(2026-06-10, v3.114.0).

Pins: (a) the exact rational HS cone table q_F = (2m'^2 - b(b+1))/(b(b+1))
(sympy CG arithmetic — compression invisible to the HS cone); (b) the
class-diagonal definite-type inertia of the causal form (graded cone, NOT a
Lorentzian signature); (c) the op-norm null ray at C^1_{+-1,0}; (d) the
NEGATIVE: rate-level reverse triangle fails on the timelike sector.

Driver: tests/wh7_support/wh7_b3_phase2_cone_structure.py (pinned by this reference).
"""
import sys
from pathlib import Path

import pytest
from sympy import Rational, S

sys.path.insert(0, str(Path(__file__).resolve().parent / "wh7_support"))
import wh7_b3_phase2_cone_structure as p2  # noqa: E402


def test_exact_hs_cone_ratios():
    half = S(1) / 2
    assert p2.q_exact(half, half, half) == Rational(-1, 3)
    assert p2.q_exact(S(1), S(1), S(1)) == 0              # null locus exact
    assert p2.q_exact(S(1), S(0), S(1)) == Rational(-1)
    assert p2.q_exact(S(3) / 2, S(3) / 2, S(3) / 2) == Rational(1, 5)
    assert p2.q_exact(S(3) / 2, half, S(3) / 2) == Rational(-13, 15)
    assert p2.q_exact(S(2), S(2), S(2)) == Rational(1, 3)
    assert p2.q_exact(S(2), S(1), S(2)) == Rational(-2, 3)


def test_hs_cone_equals_symbol_classifier():
    """q_F = (2 m'^2 - b(b+1)) / (b(b+1)) — compression-invisible."""
    for b, mp in [(S(1) / 2, S(1) / 2), (S(1), S(1)), (S(3) / 2, S(1) / 2),
                  (S(2), S(2)), (S(2), S(1))]:
        ideal = Rational(2 * mp ** 2 - b * (b + 1), b * (b + 1))
        assert p2.q_exact(b, mp, b) == ideal, (b, mp)


@pytest.fixture(scope="module")
def resB():
    return p2.part_B()


def test_global_inertia(resB):
    assert tuple(resB["inertia_global"]) == (18, 7, 30)


def test_timelike_sectors_positive_definite(resB):
    assert tuple(resB["per_class"]["(1.5, 1.5)"]["inertia"]) == (8, 0, 0)
    assert tuple(resB["per_class"]["(2.0, 2.0)"]["inertia"]) == (10, 0, 0)


def test_b1_top_class_identically_null(resB):
    assert tuple(resB["per_class"]["(1.0, 1.0)"]["inertia"]) == (0, 6, 0)


def test_op_norm_null_ray_at_m0(resB):
    assert abs(resB["b1_top_null_residual"][0]) < 1e-12


def test_reverse_triangle_fails_at_rate_level(resB):
    """The frozen NEGATIVE: timelike sector is positive-definite, so tau is
    subadditive — no operator-rate twin paradox."""
    t = resB["tau_superadd"]
    assert t["trials"] > 50
    assert t["violations"] == t["trials"]
