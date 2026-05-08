"""Tests for geovac.phillips_kleinman_cross_center.

Sprint: post-Track-2 chemistry-solver re-test, named engineering closure
for the W1b-residual orthogonality wall (Phase C-W1c memo §6).
"""

from __future__ import annotations

import numpy as np
import pytest

from geovac.phillips_kleinman_cross_center import (
    _core_orbitals_for_Z,
    _hydrogenic_radial_at_points,
    _phi_integral_factor,
    _real_sh_normalization,
    compute_pk_cross_center_barrier,
    cross_center_overlap,
)


# ---------------------------------------------------------------------------
# Real SH normalization sanity
# ---------------------------------------------------------------------------


def test_real_sh_norm_l0_m0():
    val = _real_sh_normalization(0, 0)
    assert val == pytest.approx(np.sqrt(1.0 / (4 * np.pi)), rel=1e-12)


def test_real_sh_norm_l1_m0():
    val = _real_sh_normalization(1, 0)
    assert val == pytest.approx(np.sqrt(3.0 / (4 * np.pi)), rel=1e-12)


def test_real_sh_norm_l1_mpm1():
    val_pos = _real_sh_normalization(1, 1)
    val_neg = _real_sh_normalization(1, -1)
    expected = np.sqrt(2.0 * 3.0 / (4 * np.pi) * 1.0 / 2.0)  # (l-|m|)!/(l+|m|)! = 0!/2! = 1/2
    assert val_pos == pytest.approx(expected, rel=1e-12)
    assert val_neg == pytest.approx(expected, rel=1e-12)


def test_phi_integral_factor():
    assert _phi_integral_factor(0, 0) == pytest.approx(2 * np.pi)
    assert _phi_integral_factor(1, 1) == pytest.approx(np.pi)
    assert _phi_integral_factor(-1, -1) == pytest.approx(np.pi)
    assert _phi_integral_factor(1, -1) == pytest.approx(0.0)
    assert _phi_integral_factor(0, 1) == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# Hydrogenic radial wavefunction normalization
# ---------------------------------------------------------------------------


def test_hydrogenic_radial_normalization_1s():
    """Integral |R_{1s}|^2 r^2 dr = 1 for analytical normalization."""
    r = np.linspace(0.001, 30.0, 10000)
    R = _hydrogenic_radial_at_points(1.0, 1, 0, r)
    norm = np.trapezoid(R ** 2 * r ** 2, r)
    assert norm == pytest.approx(1.0, abs=1e-3)


def test_hydrogenic_radial_normalization_2p():
    r = np.linspace(0.001, 30.0, 10000)
    R = _hydrogenic_radial_at_points(1.0, 2, 1, r)
    norm = np.trapezoid(R ** 2 * r ** 2, r)
    assert norm == pytest.approx(1.0, abs=1e-3)


# ---------------------------------------------------------------------------
# Cross-center overlap sanity
# ---------------------------------------------------------------------------


def test_cross_center_overlap_self_at_zero():
    """Overlap of an orbital with itself at R~=0 should be ~1."""
    S = cross_center_overlap(
        Z_orb_p=1.0, n_p=1, l_p=0, m_p=0,
        Z_orb_c=1.0, n_c=1, l_c=0, m_c=0,
        R_AB=0.001, n_grid_r=4000,
    )
    assert S == pytest.approx(1.0, abs=1e-3)


def test_cross_center_overlap_decays_at_large_R():
    """Overlap of two well-separated 1s orbitals should be small."""
    S = cross_center_overlap(
        Z_orb_p=1.0, n_p=1, l_p=0, m_p=0,
        Z_orb_c=1.0, n_c=1, l_c=0, m_c=0,
        R_AB=10.0,
    )
    assert abs(S) < 0.01


def test_cross_center_overlap_azimuthal_selection():
    """Different real-SH m must give zero overlap."""
    S = cross_center_overlap(
        Z_orb_p=1.0, n_p=2, l_p=1, m_p=1,
        Z_orb_c=10.0, n_c=2, l_c=1, m_c=0,
        R_AB=3.0,
    )
    assert S == 0.0


def test_cross_center_overlap_decays_with_Z():
    """Compact orbitals (high Z) overlap less at fixed R."""
    S_low = cross_center_overlap(
        Z_orb_p=1.0, n_p=1, l_p=0, m_p=0,
        Z_orb_c=1.0, n_c=1, l_c=0, m_c=0,
        R_AB=3.0,
    )
    S_high = cross_center_overlap(
        Z_orb_p=1.0, n_p=1, l_p=0, m_p=0,
        Z_orb_c=10.0, n_c=1, l_c=0, m_c=0,
        R_AB=3.0,
    )
    assert abs(S_low) > abs(S_high)


# ---------------------------------------------------------------------------
# Core orbital data
# ---------------------------------------------------------------------------


def test_core_orbitals_for_Na():
    """[Ne] core has 5 orbitals: 1s, 2s, 2p_{-1,0,+1}."""
    orbs = _core_orbitals_for_Z(11, 'Ne')
    assert len(orbs) == 5
    nlms = [(o['n'], o['l'], o['m']) for o in orbs]
    assert (1, 0, 0) in nlms
    assert (2, 0, 0) in nlms
    assert (2, 1, -1) in nlms
    assert (2, 1, 0) in nlms
    assert (2, 1, 1) in nlms


def test_core_orbital_energies_negative():
    """Bound core orbitals have negative HF eigenvalues."""
    orbs = _core_orbitals_for_Z(11, 'Ne')
    for o in orbs:
        assert o['energy'] < 0


def test_core_orbital_unsupported_returns_empty():
    """Unsupported core type returns empty list (PK reverts to zero)."""
    orbs = _core_orbitals_for_Z(19, 'Ar')
    assert orbs == []


def test_core_orbitals_first_row_returns_empty():
    """First-row Z (no frozen core) returns empty."""
    orbs = _core_orbitals_for_Z(3, 'Ne')  # Lithium with wrong core_type
    # _core_orbitals_for_Z requires explicit core_type but Li (Z=3) isn't
    # in the [Ne] table, so should return empty
    assert orbs == []


# ---------------------------------------------------------------------------
# PK barrier matrix
# ---------------------------------------------------------------------------


def test_pk_barrier_first_row_is_zero():
    """First-row Z_nuc has no frozen core: PK barrier is zeros."""
    H_states = [(1, 0, 0), (2, 0, 0), (2, 1, -1), (2, 1, 0), (2, 1, 1)]
    delta = compute_pk_cross_center_barrier(
        Z_orb=1.0, valence_states=H_states,
        Z_nuc=3.0,  # Li (no frozen core, Z<11)
        R_AB=3.015,
    )
    assert np.allclose(delta, 0.0, atol=1e-15)


def test_pk_barrier_symmetric():
    """PK barrier matrix should be symmetric (S real -> S diag(w) S^T sym)."""
    H_states = [(1, 0, 0), (2, 0, 0), (2, 1, -1), (2, 1, 0), (2, 1, 1)]
    delta = compute_pk_cross_center_barrier(
        Z_orb=1.0, valence_states=H_states,
        Z_nuc=11.0, R_AB=3.5,
    )
    assert np.allclose(delta, delta.T, atol=1e-12)


def test_pk_barrier_repulsive_at_E_v_zero():
    """With E_v=0, PK barrier is purely repulsive (positive eigenvalues)."""
    H_states = [(1, 0, 0), (2, 0, 0), (2, 1, -1), (2, 1, 0), (2, 1, 1)]
    delta = compute_pk_cross_center_barrier(
        Z_orb=1.0, valence_states=H_states,
        Z_nuc=11.0, R_AB=3.5, E_valence_ref=0.0,
    )
    eigvals = np.linalg.eigvalsh(delta)
    # All eigenvalues should be >= 0 (Delta H = S diag(|E_c|) S^T is PSD)
    assert np.all(eigvals >= -1e-12)


def test_pk_barrier_decays_with_R():
    """Barrier should decrease as R increases (overlaps shrink)."""
    H_states = [(1, 0, 0)]
    d_short = compute_pk_cross_center_barrier(
        Z_orb=1.0, valence_states=H_states,
        Z_nuc=11.0, R_AB=2.5,
    )
    d_long = compute_pk_cross_center_barrier(
        Z_orb=1.0, valence_states=H_states,
        Z_nuc=11.0, R_AB=5.0,
    )
    assert d_short[0, 0] > d_long[0, 0] > 0


def test_pk_barrier_direction_z_axis():
    """Direction (0,0,1) should match z-axis (no rotation)."""
    H_states = [(1, 0, 0), (2, 0, 0), (2, 1, -1), (2, 1, 0), (2, 1, 1)]
    d_no_dir = compute_pk_cross_center_barrier(
        Z_orb=1.0, valence_states=H_states,
        Z_nuc=11.0, R_AB=3.0,
    )
    d_with_dir = compute_pk_cross_center_barrier(
        Z_orb=1.0, valence_states=H_states,
        Z_nuc=11.0, R_AB=3.0,
        direction=(0.0, 0.0, 1.0),
    )
    assert np.allclose(d_no_dir, d_with_dir, atol=1e-10)


def test_pk_barrier_NaH_magnitude():
    """Sanity: NaH PK barrier on H 1s at R=3.5 is in mHa range."""
    H_states = [(1, 0, 0)]
    delta = compute_pk_cross_center_barrier(
        Z_orb=1.0, valence_states=H_states,
        Z_nuc=11.0, R_AB=3.5,
    )
    # H 1s self-barrier: should be sub-Ha (small but nonzero)
    val = delta[0, 0]
    assert 1e-4 < val < 1e-1


def test_pk_barrier_returns_correct_shape():
    """Return matrix shape matches valence_states length."""
    states = [(1, 0, 0), (2, 0, 0), (2, 1, 0)]
    delta = compute_pk_cross_center_barrier(
        Z_orb=1.0, valence_states=states,
        Z_nuc=11.0, R_AB=3.0,
    )
    assert delta.shape == (3, 3)
