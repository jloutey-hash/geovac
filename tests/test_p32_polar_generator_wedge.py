"""Paper 32 / Paper 42 Def 5.1: the polar rotation generator does NOT preserve the wedge.

Backs the 2026-10-03 correction (/qa trunk FULL #12): K = 2J_z anticommutes with
the wedge reflection R_polar (m_j -> -m_j), so its compression P_W K P_W to the
"hemispheric" wedge is the zero operator, and the wedge operator K_alpha^W used
for the wedge KMS state is a diagonal relabelling (two_m_j on the m_j > 0 states),
not a restriction of K.

Rejected wrong answer: "K = 2J_z preserves the wedge and K_alpha^W is its
restriction" -- then P_W K P_W would be nonzero and carry the spectrum of
K_alpha^W.
"""
from __future__ import annotations

import numpy as np
import pytest

from geovac.modular_hamiltonian import for_bisognano_wichmann


@pytest.mark.parametrize("n_max", [2, 3, 4])
def test_polar_generator_anticommutes_with_wedge_reflection(n_max: int) -> None:
    mh = for_bisognano_wichmann(n_max=n_max)
    K = mh.K_geometric
    R = mh.wedge.reflection_matrix(mh.basis)
    assert np.allclose(R @ R, np.eye(R.shape[0]))
    assert np.linalg.norm(K) > 1.0
    assert np.linalg.norm(K @ R + R @ K) == pytest.approx(0.0, abs=1e-12)


@pytest.mark.parametrize("n_max", [2, 3, 4])
def test_wedge_compression_of_polar_generator_is_zero(n_max: int) -> None:
    mh = for_bisognano_wichmann(n_max=n_max)
    K = mh.K_geometric
    P = mh.wedge.projection_matrix(mh.basis)
    assert np.trace(P).real == pytest.approx(K.shape[0] / 2)
    assert np.linalg.norm(P @ K @ P) == pytest.approx(0.0, abs=1e-12)


@pytest.mark.parametrize("n_max", [2, 3])
def test_wedge_operator_is_a_relabelling_not_the_compression(n_max: int) -> None:
    """K_alpha^W has the nonzero odd spectrum two_m_j > 0, which the (zero)
    compression of K cannot carry."""
    mh = for_bisognano_wichmann(n_max=n_max)
    KW = mh.restrict_K_alpha_to_wedge()
    ev = np.sort(np.linalg.eigvalsh(KW).real)
    assert ev.min() >= 1.0 - 1e-12
    assert np.allclose(ev, np.round(ev)) and np.all(np.round(ev) % 2 == 1)
    P = mh.wedge.projection_matrix(mh.basis)
    compressed_ev = np.linalg.eigvalsh(P @ mh.K_geometric @ P).real
    assert np.max(np.abs(compressed_ev)) == pytest.approx(0.0, abs=1e-12)
