"""Closed-form Tomita modular Hamiltonian (geovac.modular_hamiltonian_lorentzian).

Since 2026-10-03 the Krein-side Tomita structure builds K_TT = -log Delta in closed
form, K_TT = (log rho)^T (x) I - I (x) log rho, instead of diagonalising the
(dim_W^2 x dim_W^2) Delta (which ran out of memory at n_max = 3, N_t >= 11).

Rejected wrong answers: a sign error in either Kronecker term, or a transposed
factor, which agree with the eigendecomposition route only when rho is diagonal
(the production case) -- hence the NON-diagonal, well-conditioned rho here.
"""
from __future__ import annotations

import numpy as np

from geovac.modular_hamiltonian_lorentzian import for_bisognano_wichmann_lorentzian


def test_closed_form_k_tt_matches_eigendecomposition_on_nondiagonal_rho() -> None:
    lmh = for_bisognano_wichmann_lorentzian(n_max=2, N_t=1, T_max=1.0)
    d = lmh.dim_W_L
    rng = np.random.default_rng(7)
    A = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    H = (A + A.conj().T) / 400           # rho well conditioned: eigh of Delta is accurate
    ts = lmh.tomita_structure(H_local=H)
    w, V = np.linalg.eigh(ts.Delta_matrix)
    K_ref = (V * -np.log(w)) @ V.conj().T
    assert np.linalg.norm(K_ref) > 0.5
    assert np.linalg.norm(ts.K_TT - K_ref) < 1e-12


def test_delta_is_conjugation_by_rho() -> None:
    """Delta(vec X) = vec(rho X rho^{-1}) for a non-diagonal rho."""
    lmh = for_bisognano_wichmann_lorentzian(n_max=2, N_t=1, T_max=1.0)
    d = lmh.dim_W_L
    rng = np.random.default_rng(3)
    A = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    ts = lmh.tomita_structure(H_local=(A + A.conj().T) / 400)
    X = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    vecX = X.T.reshape(-1)
    lhs = (ts.Delta_matrix @ vecX).reshape(d, d).T
    rhs = ts.rho @ X @ np.linalg.inv(ts.rho)
    assert np.linalg.norm(lhs - rhs) < 1e-10 * np.linalg.norm(rhs)
