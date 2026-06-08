"""Regression tests for `geovac.symmetry_adapted_basis`.

Verifies:
  - Combined rotation U_total = U_swap @ U_hopf is orthogonal to < 1e-12
  - U_hopf and U_swap commute (structural)
  - Rotated h1 is block-diagonal across Z_2^3 sectors (max off-sector < 1e-10)
  - Rotated ERI is block-structured (4-index sector products preserved)
  - Hidden-Z_2 candidates pass commutator audit against full H
  - Fully-symmetric ground state preserved under extended+hidden tapering
"""

from __future__ import annotations

import numpy as np
import pytest


@pytest.fixture(scope='module')
def lih_spec_default():
    from geovac.molecular_spec import lih_spec
    return lih_spec()


@pytest.fixture(scope='module')
def lih_spec_max_n1():
    """Tiny LiH for exact-spectrum comparison."""
    from geovac.molecular_spec import lih_spec
    return lih_spec(max_n=1)


@pytest.fixture(scope='module')
def beh2_spec_default():
    from geovac.molecular_spec import beh2_spec
    return beh2_spec()


@pytest.fixture(scope='module')
def beh2_nuclei():
    from geovac.balanced_coupled import _get_nuclei_for_beh2
    return _get_nuclei_for_beh2(2.5)


@pytest.fixture(scope='module')
def h2o_spec_default():
    from geovac.molecular_spec import h2o_spec
    return h2o_spec()


@pytest.fixture(scope='module')
def h2o_nuclei():
    from geovac.balanced_coupled import _get_nuclei_for_h2o
    return _get_nuclei_for_h2o()


# ---------------------------------------------------------------------------
# Rotation construction
# ---------------------------------------------------------------------------

def test_combined_rotation_orthogonal(lih_spec_default):
    from geovac.symmetry_adapted_basis import build_symmetry_adapted_rotation

    U, _, _, _, ot, diag = build_symmetry_adapted_rotation(lih_spec_default)
    M = len(ot)
    err = float(np.max(np.abs(U @ U.T - np.eye(M))))
    assert err < 1e-12, f"orthogonality residual {err:.2e}"
    assert diag['orthogonality_err'] < 1e-12


def test_hopf_and_swap_commute_no_nuclei(lih_spec_default):
    """Without nuclei → U_swap = I, so they trivially commute."""
    from geovac.symmetry_adapted_basis import build_symmetry_adapted_rotation

    _, _, _, _, _, diag = build_symmetry_adapted_rotation(lih_spec_default)
    assert diag['commute_hopf_swap'] < 1e-12


def test_hopf_and_swap_commute_with_swap(beh2_spec_default, beh2_nuclei):
    """With BeH2 H1↔H2 swap → rotations act on distinct labels; must commute."""
    from geovac.symmetry_adapted_basis import build_symmetry_adapted_rotation

    _, _, _, _, _, diag = build_symmetry_adapted_rotation(
        beh2_spec_default, beh2_nuclei,
    )
    assert diag['commute_hopf_swap'] < 1e-12, \
        f"Hopf and swap don't commute: {diag['commute_hopf_swap']:.2e}"


def test_sector_labels_have_three_axes(lih_spec_default):
    from geovac.symmetry_adapted_basis import build_symmetry_adapted_rotation

    _, labels, _, _, _, _ = build_symmetry_adapted_rotation(lih_spec_default)
    assert labels.shape[1] == 3
    assert set(labels.ravel().tolist()) <= {-1, 1}


# ---------------------------------------------------------------------------
# Block-diagonality
# ---------------------------------------------------------------------------

def test_rotated_h1_block_diagonal_lih(lih_spec_default):
    from geovac.symmetry_adapted_basis import decompose_hamiltonian_by_sector

    sectors = decompose_hamiltonian_by_sector(
        lih_spec_default, builder='composed',
    )
    meta = sectors.pop(('__meta__',))
    assert meta['max_offsector_h1'] < 1e-10, (
        f"h1 max off-sector = {meta['max_offsector_h1']:.2e}"
    )


def test_rotated_eri_block_diagonal_lih(lih_spec_default):
    from geovac.symmetry_adapted_basis import decompose_hamiltonian_by_sector

    sectors = decompose_hamiltonian_by_sector(
        lih_spec_default, builder='composed',
    )
    meta = sectors.pop(('__meta__',))
    # ERI off-sector entries must vanish (the rotation respects Z_2^3)
    assert meta['max_offsector_eri'] < 1e-10, (
        f"eri max off-sector = {meta['max_offsector_eri']:.2e}"
    )


def test_rotated_h1_block_diagonal_beh2(beh2_spec_default, beh2_nuclei):
    from geovac.symmetry_adapted_basis import decompose_hamiltonian_by_sector

    sectors = decompose_hamiltonian_by_sector(
        beh2_spec_default, nuclei=beh2_nuclei, builder='composed',
    )
    meta = sectors.pop(('__meta__',))
    assert meta['max_offsector_h1'] < 1e-10, (
        f"BeH2 h1 max off-sector = {meta['max_offsector_h1']:.2e}"
    )
    assert meta['max_offsector_eri'] < 1e-10


def test_sector_count_lih(lih_spec_default):
    """LiH default expected sectors: 3 (no swap, only Hopf × ell)."""
    from geovac.symmetry_adapted_basis import decompose_hamiltonian_by_sector

    sectors = decompose_hamiltonian_by_sector(
        lih_spec_default, builder='composed',
    )
    sectors.pop(('__meta__',))
    # Z_2^3 has 8 cells; only those with nonzero orbital count are returned.
    # LiH: 3 sub-blocks × per-(n,l,m) = 15 orbitals split into
    # (+1,+1,+1), (-1,-1,+1), (+1,-1,+1)
    assert len(sectors) == 3


def test_sector_dims_sum_to_M(lih_spec_default):
    """Sum of sector dims should equal M (orbital count)."""
    from geovac.z2_tapering import _enumerate_orbitals
    from geovac.symmetry_adapted_basis import decompose_hamiltonian_by_sector

    M = len(_enumerate_orbitals(lih_spec_default))
    sectors = decompose_hamiltonian_by_sector(
        lih_spec_default, builder='composed',
    )
    sectors.pop(('__meta__',))
    total = sum(s['dim'] for s in sectors.values())
    assert total == M


# ---------------------------------------------------------------------------
# Hidden-Z_2 detection + audit
# ---------------------------------------------------------------------------

def test_hidden_z2_candidates_audit_pass_lih(lih_spec_default):
    from geovac.symmetry_adapted_basis import find_hidden_z2_in_sector

    candidates = find_hidden_z2_in_sector(
        lih_spec_default, builder='composed',
    )
    # Every emitted candidate should pass the commutator audit when
    # constructed by the per-component Z-string method (the structure of
    # the composed builder guarantees per-sub-block conservation).
    assert len(candidates) > 0
    for c in candidates:
        assert c['is_valid'], (
            f"candidate {c['sector']}, c={c['component_index']} "
            f"failed audit: residual {c['audit_residual']:.2e}"
        )


def test_hidden_z2_saves_qubits_lih(lih_spec_default):
    """LiH (composed, no swap): expect strictly positive hidden-Z_2 savings."""
    from geovac.symmetry_adapted_basis import (
        extended_plus_hidden_tapered_from_spec,
    )
    from geovac.extended_tapering import extended_tapered_from_spec

    ext = extended_tapered_from_spec(
        lih_spec_default, use_hopf=True, use_ell_parity=True,
        use_atom_swap=False,
    )
    eph = extended_plus_hidden_tapered_from_spec(
        lih_spec_default, use_hopf=True, use_ell_parity=True,
        use_atom_swap=False,
    )
    save = ext['Q_tapered'] - eph['Q_tapered_hidden']
    # LiH expected: ext Q=22, eph Q=20 → 2 hidden saved
    assert save >= 1, (
        f"LiH expected at least 1 hidden saving; got {save}"
    )


# ---------------------------------------------------------------------------
# Spectrum preservation (small-Q only)
# ---------------------------------------------------------------------------

def test_ground_state_preserved_lih_max_n1(lih_spec_max_n1):
    """LiH max_n=1 is small enough to fully diagonalize. The fully-symmetric
    ground state must match between extended and extended+hidden tapering
    to machine precision."""
    from openfermion.linalg import get_sparse_operator
    from geovac.extended_tapering import extended_tapered_from_spec
    from geovac.symmetry_adapted_basis import (
        extended_plus_hidden_tapered_from_spec,
    )

    ext = extended_tapered_from_spec(
        lih_spec_max_n1, use_hopf=True, use_ell_parity=True,
    )
    eph = extended_plus_hidden_tapered_from_spec(
        lih_spec_max_n1, use_hopf=True, use_ell_parity=True,
    )

    if ext['Q_tapered'] > 16:
        pytest.skip("too large for full-diag verification")

    def lowest(qop, n):
        sparse = get_sparse_operator(qop, n_qubits=n).todense()
        eigs = np.linalg.eigvalsh(np.asarray(sparse))
        return float(np.min(eigs))

    E_ext = lowest(ext['qubit_op_tapered'], ext['Q_tapered'])
    E_eph = lowest(
        eph['qubit_op_tapered_plus_hidden'],
        eph['Q_tapered_hidden'],
    )
    # The ground state must be in the surviving sector of the tighter tapering
    assert abs(E_ext - E_eph) < 1e-10, (
        f"GS shifted: ext={E_ext:.10f}, eph={E_eph:.10f}, "
        f"delta={abs(E_ext-E_eph):.2e}"
    )


# ---------------------------------------------------------------------------
# Verification panel coverage
# ---------------------------------------------------------------------------

def test_panel_no_offsector_leakage():
    """Across the representative panel, block-diagonality holds at the same
    machine-precision threshold."""
    from geovac.molecular_spec import (
        lih_spec, beh2_spec, h2o_spec, nh3_spec, ch4_spec, hf_spec,
    )
    from geovac.balanced_coupled import (
        _get_nuclei_for_beh2, _get_nuclei_for_h2o,
    )
    from geovac.symmetry_adapted_basis import decompose_hamiltonian_by_sector

    panel = [
        ('LiH', lih_spec(), None),
        ('HF', hf_spec(), None),
        ('BeH2', beh2_spec(), _get_nuclei_for_beh2(2.5)),
        ('H2O', h2o_spec(), _get_nuclei_for_h2o()),
        ('NH3', nh3_spec(), None),
        ('CH4', ch4_spec(), None),
    ]
    for name, spec, nuc in panel:
        sectors = decompose_hamiltonian_by_sector(
            spec, nuclei=nuc, builder='composed',
        )
        meta = sectors.pop(('__meta__',))
        assert meta['max_offsector_h1'] < 1e-10, (
            f"{name} h1 off-sector = {meta['max_offsector_h1']:.2e}"
        )
        assert meta['max_offsector_eri'] < 1e-10, (
            f"{name} eri off-sector = {meta['max_offsector_eri']:.2e}"
        )
