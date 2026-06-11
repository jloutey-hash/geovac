# -*- coding: utf-8 -*-
"""Frozen falsifier for the WH7 Step-1 Toeplitz temporal probe (2026-06-10, v3.111.0).

Pins the four load-bearing facts of debug/wh7_toeplitz_temporal_probe.py:
  1. L(S_q) = 2*pi*q/T exactly (translation seminorm = continuum Lipschitz on modes)
  2. P45 momentum-diagonal architecture is invisible (L = 0) -- the control
  3. kernel condition: L = 0 iff f constant on the band-limited algebra
  4. fixed physical frequency survives de-compactification (L = omega, T-independent)

If any of these regress, the WH7 register entry (CLAUDE.md S1.7) and the Lorentzian
rebuild path are both affected.
"""
import numpy as np
import pytest


def shift_matrix(K, q):
    N = 2 * K + 1
    S = np.zeros((N, N), dtype=complex)
    for k in range(-K, K + 1):
        if -K <= k + q <= K:
            S[k + q + K, k + K] = 1.0
    return S


def seminorm(F, K, T, n_grid=300):
    w = 2 * np.pi * np.arange(-K, K + 1) / T
    Dt = np.diag(w)
    best = np.linalg.norm(Dt @ F - F @ Dt, 2)
    for s in np.linspace(T / (2 * n_grid), T / 2, n_grid):
        U = np.exp(-1j * w * s)
        Fs = (U[:, None] * F) * np.conj(U)[None, :]
        best = max(best, np.linalg.norm(Fs - F, 2) / s)
    return float(best)


@pytest.mark.parametrize("T", [1.0, 2 * np.pi])
@pytest.mark.parametrize("q", [1, 2, 3])
def test_single_mode_seminorm_exact(T, q):
    K = 8
    L = seminorm(shift_matrix(K, q), K, T)
    assert abs(L - 2 * np.pi * q / T) < 1e-10


def test_p45_architecture_invisible():
    K = 8
    rng = np.random.default_rng(0)
    g = np.diag(rng.normal(size=2 * K + 1))  # non-constant g(D_t), momentum-diagonal
    assert seminorm(g, K, T=2 * np.pi) < 1e-12


def test_kernel_condition():
    K, T = 6, 2 * np.pi
    assert seminorm(2.5 * shift_matrix(K, 0), K, T) < 1e-12          # constants -> 0
    F = shift_matrix(K, 1) + shift_matrix(K, -1)                     # cos mode -> > 0
    assert seminorm(F, K, T) > 0.5


def test_lipschitz_domination():
    K, T = 16, 2 * np.pi
    F = 0.5 * (shift_matrix(K, 1) + shift_matrix(K, -1))             # f = cos t, Lip = 1
    L = seminorm(F, K, T)
    assert L <= 1.0 + 1e-9
    assert L > 0.9                                                    # window ratio near 1


@pytest.mark.parametrize("T", [1.0, 4.0, 16.0])
def test_fixed_frequency_decompactification(T):
    q = int(round(T))                                                 # omega = 2*pi fixed
    K = 2 * q + 2
    L = seminorm(shift_matrix(K, q), K, T)
    assert abs(L - 2 * np.pi) < 1e-10
