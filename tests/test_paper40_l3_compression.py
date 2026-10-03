"""Paper 40 Lemma L3 (2026-10-02): the Lipschitz comparison holds with C_3(G) = 1
at every rank, by compression.

The proof has two inputs: (a) the truncation projection P commutes with the
Dirac operator (left-equivariance of the Kostant operator; P is a left-isotypic
projection), and (b) compression contracts operator norms, so
    ||[P D P, P M P]|| = ||P [D, M] P|| <= ||[D, M]||.
Sharpness: P -> 1 strongly, so the compressed norm tends to ||[D, M]||.

What these tests reject:
  * "the bound holds for ANY projection" -- a projection NOT commuting with D
    can make ||[PDP, PMP]|| exceed ||[D, M]|| (test_hypothesis_is_load_bearing);
  * "the bound survives without [D,P] = 0" (rejected by exhibiting bound
    failure; the identity [D_P, PMP] = P[D,M]P itself is not tested
    separately, and no non-abelian group instance is computed);
  * "the constant can be lowered below 1 uniformly in the cutoff" -- on the
    circle (the simplest compact Lie group, D = -i d/dtheta, f = cos) the
    compressed norm climbs to ||f'||_inf = 1 (test_sharpness_on_circle).
"""
from __future__ import annotations

import numpy as np


def _herm(rng, n):
    a = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
    return (a + a.conj().T) / 2


def _spectral_projection(D, keep):
    w, v = np.linalg.eigh(D)
    idx = np.argsort(np.abs(w))[:keep]
    return v[:, idx] @ v[:, idx].conj().T


def _norm(a):
    return np.linalg.norm(a, 2)


def test_compression_contracts_when_projection_commutes():
    rng = np.random.default_rng(40)
    for trial in range(50):
        n = 24
        D, M = _herm(rng, n), _herm(rng, n)
        P = _spectral_projection(D, keep=rng.integers(3, n))
        assert _norm(D @ P - P @ D) < 1e-10
        DP, MP = P @ D @ P, P @ M @ P
        lhs = DP @ MP - MP @ DP
        # identity used in the proof
        assert _norm(lhs - P @ (D @ M - M @ D) @ P) < 1e-9
        # the bound
        assert _norm(lhs) <= _norm(D @ M - M @ D) + 1e-9


def test_hypothesis_is_load_bearing():
    """Without [D, P] = 0 the bound can fail: the commuting hypothesis is what the
    proof uses (supplied in Paper 40 by left-equivariance of the Kostant operator)."""
    rng = np.random.default_rng(7)
    # Witness: M commutes with D, so [D, M] = 0, but the compressions of two commuting
    # operators by a projection that commutes with neither need not commute.
    D = np.diag([1.0, 2.0, 3.0])
    M = np.diag([5.0, 1.0, 7.0])
    assert _norm(D @ M - M @ D) == 0.0
    q, _ = np.linalg.qr(rng.standard_normal((3, 2)))
    P = q @ q.T                                 # generic rank-2 projection
    assert _norm(D @ P - P @ D) > 1e-3
    DP, MP = P @ D @ P, P @ M @ P
    assert _norm(DP @ MP - MP @ DP) > 1e-3, "non-commuting P should break the bound here"


def test_sharpness_on_circle():
    """Circle: D = -i d/dtheta on Fourier modes |k| <= K, f = cos(theta), so
    [D, M_f] = -i f' and ||f'||_inf = 1.  The truncation P_K commutes with D; the
    compressed commutator norm stays <= 1 and tends to 1 as K grows."""
    prev = 0.0
    for K in (4, 16, 64, 256):
        ks = np.arange(-K, K + 1)
        D = np.diag(ks.astype(float))
        n = len(ks)
        M = np.zeros((n, n))
        for i in range(n - 1):                  # cos = (e^{i th} + e^{-i th}) / 2
            M[i, i + 1] = M[i + 1, i] = 0.5
        c = _norm(D @ M - M @ D)
        assert c <= 1.0 + 1e-12
        assert c >= prev - 1e-12
        prev = c
    assert prev > 0.999, prev
