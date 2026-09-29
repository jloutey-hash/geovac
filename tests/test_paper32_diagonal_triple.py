"""Paper 32: the diagonal representative (A_GV, H_GV, D_GV) is a SPECTRAL, not
a metric, object, and its finite truncation is EVEN.

Claims guarded (papers/group1_operator_algebras/paper_32_spectral_triple.tex):

  (A) Proposition prop:bounded_commutators + the [SYMBOLIC PROOF] paragraph
      after it: D_GV and pi(f) are both diagonal in |n_D, kappa, m_j>, so
      [D_GV, pi(f)] = 0 identically for every f in A_GV (no Lipschitz
      geometry, no inner fluctuations).  The metric role is carried by the
      NON-diagonal multipliers of the truncated operator system O_{n_max},
      where the commutators are nonzero.
      WRONG ANSWER REJECTED: a nonzero commutator on the diagonal
      representative (e.g. a non-diagonal D_GV, or a non-diagonal pi(f)).

  (B) Proposition prop:odd(b): gamma_GV = sigma_x (x) 1 (exchanging the two
      chirality copies) is a self-adjoint involution that commutes with
      pi(A_GV) and anticommutes with D_GV = sigma_z (x) (n + 1/2).
      WRONG ANSWER REJECTED: "the finite truncation is odd" -- the natural
      diagonal candidate sigma_z (x) 1 COMMUTES with D_GV instead, so a
      grading built that way fails the anticommutation.

pi(f) is realised as the paper states it (Eq. eq:representation): diagonal
in the Dirac basis, the value at each basis vector depending only on the
Fock node above it and not on the chirality label.  Any such diagonal
pullback is covered: f is drawn at random on the node set.

Fire-test plants (debug/qa/fire_test.py):
  * --plant-in geovac/full_dirac_operator_system.py: D_GV builder returns a
    tridiagonal (non-diagonal) matrix -> test_diagonal_commutator_vanishes FIRES;
  * in this file: pi(f) given a nearest-neighbour off-diagonal term -> FIRES;
  * --plant-in geovac/chirality_grading.py: sigma_x -> sigma_z -> the
    anticommutation test FIRES.
"""
from __future__ import annotations

import numpy as np
import pytest

from geovac.chirality_grading import build_gamma_GV
from geovac.full_dirac_operator_system import (
    FullDiracTruncatedOperatorSystem,
    camporesi_higuchi_full_dirac_matrix,
    full_dirac_basis,
)


def _fock_node(b) -> tuple:
    """Fock node (n, l, m_l) above a full-Dirac label: n, l are carried
    directly; m_l is a deterministic bridge from m_j (m_j - 1/2 clipped into
    [-l, l]).  Only 'diagonal and chirality-independent' is load-bearing."""
    m_l = (b.two_m_j - 1) // 2
    m_l = max(-b.l, min(b.l, m_l))
    return (b.n_fock, b.l, m_l)


def _pi(f: dict, basis) -> np.ndarray:
    """pi(f): the diagonal pullback of a function f on Fock nodes."""
    diag = np.array([f[_fock_node(b)] for b in basis], dtype=np.complex128)
    return np.diag(diag)


def _random_f(n_max: int, rng) -> dict:
    nodes = [(n, l, m) for n in range(1, n_max + 1)
             for l in range(n) for m in range(-l, l + 1)]
    return {v: complex(rng.standard_normal(), rng.standard_normal()) for v in nodes}


@pytest.mark.parametrize("n_max", [2, 3])
def test_diagonal_commutator_vanishes(n_max):
    """(A) [D_GV, pi(f)] = 0 identically, for random f on the node set.

    n_max >= 2 only: at n_max = 1 there is a single Fock node, so pi(f) is a
    scalar multiple of 1 and the commutator vanishes for ANY D (vacuous; see
    test_n_max_1_case_is_vacuous).  Here pi(f) is required to be genuinely
    non-scalar, so a non-diagonal D_GV cannot pass."""
    basis = full_dirac_basis(n_max)
    D = camporesi_higuchi_full_dirac_matrix(basis)
    rng = np.random.default_rng(n_max)
    for _ in range(5):
        P = _pi(_random_f(n_max, rng), basis)
        d = np.diag(P)
        assert np.ptp(d.real) + np.ptp(d.imag) > 1e-3, "pi(f) is scalar: vacuous"
        comm = D @ P - P @ D
        assert np.max(np.abs(comm)) == 0.0, (
            f"n_max={n_max}: [D_GV, pi(f)] != 0 (max {np.max(np.abs(comm)):.3e})")


def test_n_max_1_case_is_vacuous():
    """LABELLED VACUOUS, kept only to document it: at n_max = 1 the node set is
    the single node (1, 0, 0), so pi(f) = f(1,0,0) * 1 and [D, pi(f)] = 0 holds
    for EVERY D.  This asserts the vacuity (it carries no evidence for (A));
    the load-bearing cases are n_max >= 2 above."""
    basis = full_dirac_basis(1)
    P = _pi(_random_f(1, np.random.default_rng(1)), basis)
    assert np.allclose(P, P[0, 0] * np.eye(len(basis)))


@pytest.mark.parametrize("n_max", [2, 3])
def test_operator_system_commutators_are_nonzero(n_max):
    """(A, contrast) the metric role lives on the non-diagonal multipliers of
    O_{n_max}: the same D_GV has nonzero commutators there, so the vanishing
    above is a property of the diagonal representative, not of D_GV."""
    op = FullDiracTruncatedOperatorSystem(n_max)
    D = camporesi_higuchi_full_dirac_matrix(op.basis)
    norms = [np.linalg.norm(D @ M - M @ D, 2) for M in op.multiplier_matrices]
    assert max(norms) > 0.1, f"n_max={n_max}: all [D, M] ~ 0 (max {max(norms):.2e})"


@pytest.mark.parametrize("n_max", [1, 2, 3])
def test_gamma_is_selfadjoint_involution(n_max):
    """(B) gamma_GV^2 = 1 and gamma_GV = gamma_GV^dagger."""
    g = build_gamma_GV(n_max).gamma_matrix()
    assert np.max(np.abs(g @ g - np.eye(g.shape[0]))) == 0.0
    assert np.max(np.abs(g - g.conj().T)) == 0.0


@pytest.mark.parametrize("n_max", [1, 2, 3])
def test_gamma_commutes_with_diagonal_algebra(n_max):
    """(B) [gamma_GV, pi(f)] = 0 for f on the node set (pi(f) is chirality-blind)."""
    basis = full_dirac_basis(n_max)
    g = build_gamma_GV(n_max).gamma_matrix()
    rng = np.random.default_rng(10 + n_max)
    for _ in range(5):
        P = _pi(_random_f(n_max, rng), basis)
        assert np.max(np.abs(g @ P - P @ g)) == 0.0


@pytest.mark.parametrize("n_max", [1, 2, 3])
def test_gamma_anticommutes_with_D_GV(n_max):
    """(B) {gamma_GV, D_GV} = 0 bit-exactly -> the finite truncation is EVEN.
    Contrast: sigma_z (x) 1 commutes with D_GV, so it is not a grading."""
    basis = full_dirac_basis(n_max)
    D = camporesi_higuchi_full_dirac_matrix(basis)
    g = build_gamma_GV(n_max).gamma_matrix()
    assert np.max(np.abs(g @ D + D @ g)) == 0.0
    # nondegenerate: gamma is not zero on either chirality block
    assert np.linalg.norm(g @ D) > 0
    gz = build_gamma_GV(n_max, convention="sigma_z").gamma_matrix()
    assert np.max(np.abs(gz @ D + D @ gz)) > 1.0


@pytest.mark.parametrize("n_max", [1, 2, 3])
def test_sigma_z_candidate_commutes_with_D_GV(n_max):
    """(B, the paper's contrast statement) sigma_z (x) 1 COMMUTES with D_GV
    bit-exactly -- the reason the natural diagonal candidate is not a grading.
    (The anticommutator > 1 above does not imply this; this pins it.)"""
    basis = full_dirac_basis(n_max)
    D = camporesi_higuchi_full_dirac_matrix(basis)
    gz = build_gamma_GV(n_max, convention="sigma_z").gamma_matrix()
    assert np.allclose(gz @ gz, np.eye(gz.shape[0]))
    assert not np.allclose(gz, np.eye(gz.shape[0]))  # not trivially the identity
    assert np.max(np.abs(gz @ D - D @ gz)) == 0.0
