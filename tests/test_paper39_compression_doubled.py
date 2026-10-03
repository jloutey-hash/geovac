"""Paper 39 Lemma L3-T (2026-10-02): the joint Lipschitz comparison constant is
C_3^(2) = 1 relative to the classical product Lipschitz norm, by compression on
the Clifford-doubled module D = D_a (x) s1 (x) 1 + 1 (x) s2 (x) D_b.

The proof has three inputs: (a) the true-window projection P_a (x) 1 (x) P_b
commutes with D (each P is a spectral projection of its factor's Dirac);
(b) the continuum symbol c_a(d_a F) s1 + s2 c_b(d_b F) squares to
-(|d_a F|^2 + |d_b F|^2) because s1 and s2 ANTICOMMUTE, so its norm is
||grad F||_inf; (c) compression contracts operator norms.

The S^3 inputs (spectral truncation, Clifford symbol) are Paper 38's; the
product-specific step is (b).  These tests check it on the flat torus, where
D_j = -i d/dtheta_j on Fourier modes |k| <= K and the truncation IS a spectral
projection of D_j.  What they reject:
  * "the joint constant is the triangle bound sqrt(2) relative to L(f(x)g)",
    i.e. that the classical-norm constant cannot be 1: for F = sin x sin y,
    ||grad F||_inf = 1 while the triangle route allows ||grad f|| ||g|| +
    ||f|| ||grad g|| = 2;  the compressed commutator stays <= 1 and tends to 1
    (test_simple_tensor_constant_one);
  * "the anticommuting Clifford pair is decoration": with the same factors but
    commuting couplings (no s1/s2), F = cos x + cos y drives the commutator to
    |sin x| + |sin y| -> 2 > sqrt(2) = ||grad F||_inf, so the bound FAILS
    (test_anticommutation_is_load_bearing);
  * "the bound is loose": the doubled commutator tends to ||grad F||_inf
    (sharpness, both tests).
"""
from __future__ import annotations

import numpy as np

S1 = np.array([[0, 1], [1, 0]], dtype=complex)
S2 = np.array([[0, -1j], [1j, 0]], dtype=complex)
I2 = np.eye(2, dtype=complex)


def _factor(K: int):
    ks = np.arange(-K, K + 1)
    n = len(ks)
    D = np.diag(ks.astype(complex))
    cos = np.zeros((n, n), dtype=complex)
    sin = np.zeros((n, n), dtype=complex)
    for i in range(n - 1):  # cos = (e^{i th} + e^{-i th})/2, sin = (e^{i th} - e^{-i th})/(2i)
        cos[i, i + 1] = cos[i + 1, i] = 0.5
        sin[i + 1, i] = 1 / (2j)
        sin[i, i + 1] = -1 / (2j)
    return D, cos, sin, np.eye(n, dtype=complex)


def _comm_norm(A: np.ndarray, M: np.ndarray) -> float:
    return float(np.linalg.norm(A @ M - M @ A, 2))


def _kron3(a, b, c):
    return np.kron(np.kron(a, b), c)


def test_simple_tensor_constant_one() -> None:
    """F = sin x sin y: ||grad F||_inf = 1; the triangle route would allow 2."""
    prev = 0.0
    for K in (4, 8, 16):
        D, _cos, sin, I = _factor(K)
        Dd = _kron3(D, S1, I) + _kron3(I, S2, D)
        M = _kron3(sin, I2, sin)
        c = _comm_norm(Dd, M)
        assert c <= 1.0 + 1e-12, (K, c)
        assert c >= prev - 1e-12
        prev = c
    assert prev > 0.99, prev


def test_anticommutation_is_load_bearing() -> None:
    """F = cos x + cos y: ||grad F||_inf = sqrt(2).  Doubled: <= sqrt(2), -> sqrt(2).
    Commuting couplings (no Clifford pair): exceeds sqrt(2), -> 2."""
    target = np.sqrt(2.0)
    prev = 0.0
    for K in (4, 8, 16):
        D, cos, _sin, I = _factor(K)
        Dd = _kron3(D, S1, I) + _kron3(I, S2, D)
        Dn = _kron3(D, I2, I) + _kron3(I, I2, D)
        M = _kron3(cos, I2, I) + _kron3(I, I2, cos)
        cd = _comm_norm(Dd, M)
        cn = _comm_norm(Dn, M)
        assert cd <= target + 1e-12, (K, cd)
        assert cd >= prev - 1e-12
        prev = cd
        if K >= 8:
            assert cn > target + 0.4, (K, cn)
    assert prev > 0.99 * target, prev


def test_truncation_commutes_with_doubled_dirac() -> None:
    """Input (a): the product truncation commutes with the doubled Dirac.  Embed
    the K-window in a 2K-window and check P D = D P with P the window projector."""
    K, Kbig = 3, 6
    D, _c, _s, I = _factor(Kbig)
    Dd = _kron3(D, S1, I) + _kron3(I, S2, D)
    ks = np.arange(-Kbig, Kbig + 1)
    p = np.diag((np.abs(ks) <= K).astype(complex))
    P = _kron3(p, I2, p)
    assert np.linalg.norm(P @ Dd - Dd @ P) < 1e-12
