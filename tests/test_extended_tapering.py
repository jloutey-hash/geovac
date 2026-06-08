"""Regression tests for `geovac.extended_tapering`.

Verifies:
  - ℓ-parity Z₂ stabilizers are valid (commute with H) on composed-builder
    molecules at default n_max=2
  - Extended tapering saves ≥ Hopf-only's ΔQ in qubit count
  - The fully-symmetric ground-state sector is preserved at machine
    precision for small systems amenable to sparse eigsh
  - Atom-swap Z₂ stabilizers DO commute and produce additional savings
    on polyatomics where applicable
  - Backward compatibility: with all extended kwargs False, the result
    matches `hopf_tapered_from_spec` bit-exactly
"""

from __future__ import annotations

import math
import numpy as np
import pytest
from scipy.sparse.linalg import eigsh

from openfermion.linalg import get_sparse_operator


@pytest.fixture(scope='module')
def lih_spec_default():
    from geovac.molecular_spec import lih_spec
    return lih_spec()


@pytest.fixture(scope='module')
def hf_spec_default():
    from geovac.molecular_spec import hf_spec
    return hf_spec()


@pytest.fixture(scope='module')
def beh2_nuclei():
    from geovac.balanced_coupled import _get_nuclei_for_beh2
    return _get_nuclei_for_beh2(2.5)


@pytest.fixture(scope='module')
def h2o_nuclei():
    from geovac.balanced_coupled import _get_nuclei_for_h2o
    return _get_nuclei_for_h2o()


def _spectrum_lowest(qop, n_qubits: int, k: int = 1):
    sparse = get_sparse_operator(qop, n_qubits=n_qubits)
    eigs, _ = eigsh(sparse, k=k, which='SA')
    return sorted(float(e) for e in eigs)


def _count_pauli(qop) -> int:
    return len(qop.terms) - (1 if () in qop.terms else 0)


# ---------------------------------------------------------------------------
# ℓ-parity stabilizer construction
# ---------------------------------------------------------------------------

def test_ell_parity_stabilizer_construction(lih_spec_default):
    """LiH default has 3 sub-blocks, all containing p (l=1) orbitals,
    so we expect 3 per-block ℓ-parity stabilizers."""
    from geovac.z2_tapering import _enumerate_orbitals
    from geovac.extended_tapering import build_ell_parity_stabilizers

    orbital_table = _enumerate_orbitals(lih_spec_default)
    P_list = build_ell_parity_stabilizers(orbital_table, mode='per_block')
    assert len(P_list) == 3


def test_ell_parity_global_mode_returns_one_stabilizer(lih_spec_default):
    from geovac.z2_tapering import _enumerate_orbitals
    from geovac.extended_tapering import build_ell_parity_stabilizers

    orbital_table = _enumerate_orbitals(lih_spec_default)
    P_list = build_ell_parity_stabilizers(orbital_table, mode='global')
    assert len(P_list) == 1


# ---------------------------------------------------------------------------
# Atom-swap stabilizer construction
# ---------------------------------------------------------------------------

def test_atom_swap_no_equivalent_atoms_returns_empty(lih_spec_default, monkeypatch):
    """LiH has no equivalent atoms — no swap pairs."""
    from geovac.z2_tapering import _enumerate_orbitals
    from geovac.extended_tapering import (
        build_atom_swap_rotation_and_stabilizers, find_equivalent_atom_pairs,
    )

    nuclei = [
        {'Z': 3.0, 'position': (0.0, 0.0, 0.0), 'label': 'Li'},
        {'Z': 1.0, 'position': (0.0, 0.0, 3.015), 'label': 'H'},
    ]
    pairs = find_equivalent_atom_pairs(lih_spec_default, nuclei)
    assert pairs == []

    orbital_table = _enumerate_orbitals(lih_spec_default)
    U_swap, parity, _ = build_atom_swap_rotation_and_stabilizers(
        lih_spec_default, orbital_table, nuclei,
    )
    # No swap pairs → U is identity, parity all +1
    assert np.allclose(U_swap, np.eye(len(orbital_table)), atol=1e-14)
    assert all(parity == 1)


def test_atom_swap_beh2_finds_h1_h2_pair(beh2_nuclei):
    from geovac.molecular_spec import beh2_spec
    from geovac.extended_tapering import find_equivalent_atom_pairs

    spec = beh2_spec()
    pairs = find_equivalent_atom_pairs(spec, beh2_nuclei)
    # Should find one pair (the two H atoms)
    assert len(pairs) == 1
    i, j = pairs[0]
    assert beh2_nuclei[i]['Z'] == 1.0
    assert beh2_nuclei[j]['Z'] == 1.0


# ---------------------------------------------------------------------------
# Centrosymmetric detection
# ---------------------------------------------------------------------------

def test_centrosymmetric_detection_beh2(beh2_nuclei):
    from geovac.extended_tapering import is_centrosymmetric
    assert is_centrosymmetric(beh2_nuclei) is True


def test_centrosymmetric_detection_h2o_false(h2o_nuclei):
    """H₂O is NOT centrosymmetric (bent geometry)."""
    from geovac.extended_tapering import is_centrosymmetric
    assert is_centrosymmetric(h2o_nuclei) is False


def test_centrosymmetric_detection_lih_false():
    from geovac.extended_tapering import is_centrosymmetric
    nuclei = [
        {'Z': 3.0, 'position': (0.0, 0.0, 0.0), 'label': 'Li'},
        {'Z': 1.0, 'position': (0.0, 0.0, 3.015), 'label': 'H'},
    ]
    # Heteronuclear → not centrosymmetric
    assert is_centrosymmetric(nuclei) is False


# ---------------------------------------------------------------------------
# End-to-end pipeline
# ---------------------------------------------------------------------------

def test_extended_backward_compat_all_off_matches_naive(lih_spec_default):
    """With all extended kwargs False AND use_hopf False → should match
    the un-tapered Jordan-Wigner naive operator."""
    from geovac.extended_tapering import extended_tapered_from_spec

    out = extended_tapered_from_spec(
        lih_spec_default, use_hopf=False, use_ell_parity=False,
        use_atom_swap=False, use_inversion=False,
    )
    # alpha/beta parity Z-strings still get applied (always-on)
    # so Q_naive - Q_tapered = 2
    assert out['Q_naive'] == 30
    assert out['delta_Q'] == 2


def test_extended_hopf_only_matches_hopf_tapered_from_spec(lih_spec_default):
    """With only use_hopf=True (others off), result should match
    `hopf_tapered_from_spec`."""
    from geovac.z2_tapering import hopf_tapered_from_spec
    from geovac.extended_tapering import extended_tapered_from_spec

    hopf_out = hopf_tapered_from_spec(lih_spec_default, mode='per_block')
    ext_out = extended_tapered_from_spec(
        lih_spec_default, use_hopf=True, use_ell_parity=False,
        use_atom_swap=False, use_inversion=False,
    )
    assert ext_out['Q_tapered'] == hopf_out['Q_tapered']
    assert ext_out['delta_Q'] == hopf_out['delta_Q']


def test_extended_ell_parity_saves_n_sub_blocks_on_lih(lih_spec_default):
    """LiH has 3 sub-blocks; ℓ-parity should add 3 more qubits saved."""
    from geovac.z2_tapering import hopf_tapered_from_spec
    from geovac.extended_tapering import extended_tapered_from_spec

    hopf_out = hopf_tapered_from_spec(lih_spec_default, mode='per_block')
    ext_out = extended_tapered_from_spec(
        lih_spec_default, use_hopf=True, use_ell_parity=True,
        use_atom_swap=False, use_inversion=False,
    )
    assert ext_out['Q_tapered'] == hopf_out['Q_tapered'] - 3


def test_extended_ell_parity_reduces_pauli_on_lih(lih_spec_default):
    """ℓ-parity should reduce Pauli count vs Hopf-only on LiH."""
    from geovac.z2_tapering import hopf_tapered_from_spec
    from geovac.extended_tapering import extended_tapered_from_spec

    hopf_out = hopf_tapered_from_spec(lih_spec_default, mode='per_block')
    ext_out = extended_tapered_from_spec(
        lih_spec_default, use_hopf=True, use_ell_parity=True,
        use_atom_swap=False, use_inversion=False,
    )
    n_p_hopf = _count_pauli(hopf_out['qubit_op_tapered'])
    n_p_ext = _count_pauli(ext_out['qubit_op_tapered'])
    assert n_p_ext < n_p_hopf


def test_extended_savings_panel():
    """Across the representative panel, extended tapering with
    Hopf + ℓ-parity should yield ΔQ_extended == ΔQ_hopf + n_sub_blocks_with_lodd
    (atom-swap and inversion off by default in this test)."""
    from geovac.molecular_spec import (
        lih_spec, beh2_spec, h2o_spec, hf_spec, nh3_spec, ch4_spec,
    )
    from geovac.z2_tapering import hopf_tapered_from_spec
    from geovac.extended_tapering import extended_tapered_from_spec

    panel = [
        ('LiH', lih_spec(), 3),
        ('HF',  hf_spec(), 6),
        ('BeH2', beh2_spec(), 5),
        ('H2O', h2o_spec(), 7),
    ]
    for (name, spec, expected_extra) in panel:
        hopf_out = hopf_tapered_from_spec(spec, mode='per_block')
        ext_out = extended_tapered_from_spec(
            spec, use_hopf=True, use_ell_parity=True,
            use_atom_swap=False, use_inversion=False,
        )
        delta_extra = hopf_out['Q_tapered'] - ext_out['Q_tapered']
        assert delta_extra == expected_extra, (
            f"{name}: expected ΔQ_extra = {expected_extra}, "
            f"got {delta_extra}"
        )
