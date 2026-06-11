# -*- coding: utf-8 -*-
"""Frozen falsifier for B1: joint S^3 x S^1_T action-seminorm convergence machinery
(2026-06-10, v3.112.0). Pins Paper 45 Proposition `prop:product_action_seminorm`.

Imports the driver debug/wh7_b1_joint_product_gh.py (pinned at debug/ top level by
this reference). Module import builds the SU(2) quadrature grid (~seconds).
"""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "debug"))
import wh7_b1_joint_product_gh as b1  # noqa: E402

K = 2
N_T = 2 * K + 1


def test_pw_orthonormality():
    assert np.max(np.abs(b1.GRAM - np.eye(b1.NW))) < 1e-12


def test_pure_temporal_exact_in_joint():
    Sq = b1.shift_matrix(K, 1)
    L = b1.L_joint(np.kron(np.eye(b1.NW), Sq), K)
    assert abs(L - 2 * np.pi / b1.T_TIME) < 1e-10


def test_pure_spatial_matches_factor_seminorm():
    A = b1.compress_spatial(b1.band_vals(0.5, 0, 0))
    assert abs(b1.L_joint(np.kron(A, np.eye(N_T)), K) - b1.L_spatial(A)) < 1e-10


def test_joint_kernel_condition():
    const = b1.compress_spatial(np.full(b1.NG, 1.7))
    assert b1.L_joint(np.kron(const, np.eye(N_T)), K) < 1e-12
    A = b1.compress_spatial(b1.band_vals(1.0, 1, 1))
    B = b1.shift_matrix(K, 1) + b1.shift_matrix(K, -1)
    assert b1.L_joint(np.kron(A, B), K) > 0.1


def test_leibniz_envelope():
    A = b1.compress_spatial(b1.band_vals(0.5, 0, 0))
    B = b1.shift_matrix(K, 1) + b1.shift_matrix(K, -1)
    Ls, Lt = b1.L_spatial(A), b1.L_temporal(B, K)
    nA, nB = np.linalg.norm(A, 2), np.linalg.norm(B, 2)
    L = b1.L_joint(np.kron(A, B), K)
    assert max(Ls * nB, nA * Lt) - 1e-9 <= L <= Ls * nB + nA * Lt + 1e-9


def test_riemannian_limit_nt1_bit_exact():
    for (j, a, bb) in [(0.5, 0, 0), (1.0, 1, 1)]:
        A = b1.compress_spatial(b1.band_vals(j, a, bb))
        assert abs(b1.L_joint(np.kron(A, np.eye(1)), 0) - b1.L_spatial(A)) < 1e-13
