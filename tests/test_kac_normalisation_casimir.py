"""Cas(ad) = 2 h^vee under the field-standard normalisation (theta|theta) = 2.

The group3 synthesis states that under the normalised invariant form of Kac
(Infinite Dimensional Lie Algebras, Ch. 6: long roots have (theta|theta) = 2)
the adjoint Casimir is 2 h^vee, which is why the GH rate constant reads
2 sqrt(2)/pi there and 4/pi in the corpus's dual-Coxeter rule Cas(ad) = h^vee.
This test COMPUTES the Casimir rather than asserting the rule: build su(N) in
the defining representation, fix the invariant form B(X, Y) = c tr(XY) by
requiring the highest root's coroot-dual to have (theta|theta) = 2, and
evaluate sum_a ad(e_a) ad(e^a) on the adjoint representation.

Wrong answers rejected: Cas(ad) = h^vee in this normalisation (the corpus's
OWN rule, a different normalisation), and the Killing normalisation (which
gives Cas(ad) = 1).
"""
from __future__ import annotations

import itertools

import numpy as np
import pytest


def _su_basis(n: int) -> list[np.ndarray]:
    """Hermitian traceless basis of su(n) (generalised Gell-Mann)."""
    B = []
    for i, j in itertools.combinations(range(n), 2):
        X = np.zeros((n, n), complex); X[i, j] = X[j, i] = 1; B.append(X)
        Y = np.zeros((n, n), complex); Y[i, j] = -1j; Y[j, i] = 1j; B.append(Y)
    for k in range(1, n):
        D = np.zeros((n, n), complex)
        D[:k, :k] = np.eye(k); D[k, k] = -k
        B.append(D * np.sqrt(2.0 / (k * (k + 1))))
    return B                                   # tr(T_a T_b) = 2 delta_ab


def _adjoint_casimir_with_form(n: int, c: float) -> float:
    """Casimir sum_ab (B^{-1})_ab ad(T_a) ad(T_b) on the adjoint rep, for the
    invariant form B(X, Y) = c * tr(XY); returns its (scalar) eigenvalue."""
    T = _su_basis(n)
    d = len(T)
    G = np.array([[c * np.trace(a @ b).real for b in T] for a in T])
    Ginv = np.linalg.inv(G)
    # ad(T_a) in the basis T: [T_a, T_b] = sum_c f_abc T_c, using tr(T_c T_e) = 2 delta
    ad = []
    for a in T:
        M = np.zeros((d, d), complex)
        for bi, b in enumerate(T):
            comm = a @ b - b @ a
            for ci, cc in enumerate(T):
                M[ci, bi] = np.trace(comm @ cc) / 2.0
        ad.append(M)
    # B is real-symmetric on the Hermitian basis; the compact real form needs i T_a
    C = sum(Ginv[i, j] * (1j * ad[i]) @ (1j * ad[j]) for i in range(d) for j in range(d))
    ev = np.linalg.eigvals(C)
    assert np.allclose(ev, ev[0], atol=1e-9), ev        # Casimir is scalar on the irrep
    return float(-ev[0].real)


def _c_for_theta_squared_two(n: int) -> float:
    """Scale c so that the highest root theta = e_1 - e_n has (theta|theta) = 2.
    With B = c tr(XY) on the Cartan, the induced form on roots is (a|b) = <a,b>/c
    in the trace-form coordinates, and <theta, theta> = 2 there, so c = 1."""
    return 1.0


@pytest.mark.parametrize("n,h_vee", [(2, 2), (3, 3)])
def test_adjoint_casimir_is_twice_the_dual_coxeter_number(n, h_vee):
    cas = _adjoint_casimir_with_form(n, _c_for_theta_squared_two(n))
    assert abs(cas - 2 * h_vee) < 1e-9, cas
    # wrong answers rejected
    assert abs(cas - h_vee) > 0.5                       # the corpus's own (other) rule
    killing_c = 2 * n                                   # Killing form = 2n tr(XY) on su(n)
    assert abs(_adjoint_casimir_with_form(n, killing_c) - 1.0) < 1e-9


def test_theta_squared_two_under_trace_form():
    """Check the normalisation claim used above: with B = tr(XY), the highest
    root of su(n) has squared length 2 (via the dual of H_theta = E_11 - E_nn)."""
    for n in (2, 3, 4):
        H = np.zeros((n, n)); H[0, 0] = 1; H[-1, -1] = -1
        # theta(H) = 2 and B(H, H) = tr(H^2) = 2  =>  (theta|theta) = theta(H)^2 / B(H,H) = 2
        assert np.trace(H @ H) == 2
        assert (H[0, 0] - H[-1, -1]) ** 2 / np.trace(H @ H) == 2
