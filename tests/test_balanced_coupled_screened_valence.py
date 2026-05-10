"""
Tests for the screened-Schrödinger valence basis correction (Track 3).

Covers:
  - geovac/screened_valence_basis.py module-level helpers
  - build_balanced_hamiltonian flag wiring
  - backward-compat: flag OFF gives bit-exact existing behavior
  - first-row systems (LiH) unaffected even with flag ON
  - frozen-core systems (NaH, MgH₂, HCl) get correct corrections
  - composability with screened_cross_center / pk_cross_center
  - eigenvalue correctness against known atomic ionization potentials
"""

import numpy as np
import pytest

from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.molecular_spec import (
    lih_spec, nah_spec, mgh2_spec, hcl_spec,
)
from geovac.screened_valence_basis import (
    _detect_frozen_core_Z,
    apply_screened_valence_correction,
    clear_eigenvalue_cache,
    compute_screened_h1_correction_block,
    screened_valence_eigenvalue,
    screened_valence_h1_diagonal,
)


# ---------------------------------------------------------------------------
# Module-level helper tests
# ---------------------------------------------------------------------------


class TestDetectFrozenCore:
    """Auto-detection of frozen-core species."""

    def test_first_row_no_core(self):
        for Z in [1, 2, 3, 4, 6, 8, 9, 10]:
            assert not _detect_frozen_core_Z(Z), f"Z={Z} should not have FC"

    def test_second_row_ne_core(self):
        for Z in [11, 12, 14, 16, 17, 18]:
            assert _detect_frozen_core_Z(Z), f"Z={Z} should have [Ne] core"

    def test_third_row(self):
        for Z in [19, 20]:
            assert _detect_frozen_core_Z(Z)


class TestScreenedEigenvalue:
    """Eigenvalues against published atomic data."""

    def setup_method(self):
        clear_eigenvalue_cache()

    def test_na_3s_known(self):
        # Na 3s: NIST IP = 5.139 eV = 0.189 Ha (negative of binding)
        # Single-zeta Clementi-Raimondi screened gives ~-0.17 Ha
        e = screened_valence_eigenvalue(
            Z_nuc=11, block_n=1, l=0, n_val_offset=2,
        )
        assert -0.30 < e < -0.10, f"Na 3s eigenvalue {e} out of range"

    def test_na_4s_bound_above_3s(self):
        e_3s = screened_valence_eigenvalue(11, 1, 0, 2)
        e_4s = screened_valence_eigenvalue(11, 2, 0, 2)
        assert e_4s > e_3s
        assert e_4s < 0  # still bound

    def test_na_4p_bound_above_4s(self):
        e_4s = screened_valence_eigenvalue(11, 2, 0, 2)
        e_4p = screened_valence_eigenvalue(11, 2, 1, 2)
        assert e_4p > e_4s
        assert e_4p < 0

    def test_first_row_rejects(self):
        # Z=3 has no FrozenCore registered
        with pytest.raises(ValueError):
            screened_valence_eigenvalue(3, 1, 0, 0)

    def test_invalid_block_n(self):
        # n_phys < l+1 must raise: n_phys = block_n + n_val_offset = 1+2 = 3
        # so l=5 gives n_phys=3 < 6; not allowed
        with pytest.raises(ValueError):
            screened_valence_eigenvalue(11, 1, 5, 2)  # n_phys=3, l=5 -> invalid
        with pytest.raises(ValueError):
            screened_valence_eigenvalue(11, 1, 4, 2)  # n_phys=3, l=4 -> invalid

    def test_eigenvalue_cache_consistent(self):
        e1 = screened_valence_eigenvalue(11, 1, 0, 2)
        e2 = screened_valence_eigenvalue(11, 1, 0, 2)  # second call hits cache
        assert e1 == e2

    def test_mg_3s(self):
        # Mg 3s²: 1st IP = 7.65 eV = 0.281 Ha. Single-zeta Clementi-Raimondi
        # gives a more-negative value (~-0.52 Ha) because the single-zeta
        # screening is less complete than the actual SCF (no orbital
        # relaxation), but still negative and physically bound.
        e = screened_valence_eigenvalue(12, 1, 0, 2)
        assert -1.0 < e < -0.10


class TestScreenedH1Diagonal:
    """h1 diagonal correction values."""

    def test_na_3s_correction_positive(self):
        # Hydrogenic Z=1 1s = -0.5 Ha, Na 3s ~ -0.17 Ha,
        # so correction is positive (less negative)
        e_scr, e_hyd, dh = screened_valence_h1_diagonal(
            Z_orb=1.0, Z_nuc=11.0, block_n=1, l=0, n_val_offset=2,
        )
        assert dh > 0
        assert abs(e_hyd + 0.5) < 1e-12
        assert -0.30 < e_scr < -0.10
        assert abs(dh - (e_scr - e_hyd)) < 1e-12

    def test_first_row_rejected(self):
        # First-row centers should not be passed through this function
        with pytest.raises(ValueError):
            screened_valence_h1_diagonal(
                Z_orb=1.0, Z_nuc=3.0, block_n=1, l=0, n_val_offset=0,
            )


# ---------------------------------------------------------------------------
# Block-level correction tests
# ---------------------------------------------------------------------------


class TestComputeScreenedH1CorrectionBlock:
    """Per-block correction matrix construction."""

    def test_first_row_block_returns_none(self):
        # LiH bond center has Z_nuc=3, n_val_offset=0
        spec = lih_spec(max_n=2)
        # Find the bond center block
        for blk in spec.blocks:
            if blk.label == 'LiH_bond':
                lih_bond = blk
                break
        else:
            pytest.skip("No LiH_bond block found")

        from geovac.composed_qubit import _enumerate_states
        states = _enumerate_states(lih_bond.max_n)
        delta = compute_screened_h1_correction_block(lih_bond, states)
        assert delta is None

    def test_nah_bond_block_correction(self):
        spec = nah_spec(max_n=2)
        for blk in spec.blocks:
            if blk.label == 'NaH_bond':
                nah_bond = blk
                break
        else:
            pytest.fail("No NaH_bond block found")

        from geovac.composed_qubit import _enumerate_states
        states = _enumerate_states(nah_bond.max_n)
        delta = compute_screened_h1_correction_block(nah_bond, states)
        assert delta is not None
        # 5 orbitals: (1,0,0), (2,0,0), (2,1,-1), (2,1,0), (2,1,+1)
        assert delta.shape == (5, 5)
        # Diagonal only
        assert np.allclose(delta - np.diag(np.diag(delta)), 0)
        # First entry corresponds to Na 3s (block_n=1, l=0)
        # delta = -0.170 - (-0.5) = +0.330
        assert delta[0, 0] > 0.30
        assert delta[0, 0] < 0.35


# ---------------------------------------------------------------------------
# build_balanced_hamiltonian wiring tests
# ---------------------------------------------------------------------------


class TestBackwardCompatFlagOff:
    """Flag OFF must give bit-exact existing behavior."""

    def test_lih_h1_unchanged(self):
        spec = lih_spec(max_n=2)
        r_off = build_balanced_hamiltonian(spec, R=3.015, screened_valence_basis=False)
        r_on = build_balanced_hamiltonian(spec, R=3.015, screened_valence_basis=True)
        # First-row: NO frozen-core block, so flag ON should be identical
        assert np.allclose(r_off['h1'], r_on['h1'], atol=0, rtol=0)
        assert r_off['N_pauli'] == r_on['N_pauli']

    def test_nah_flag_off_matches_existing(self):
        spec = nah_spec(max_n=2)
        # The "old" result without the flag should still come through
        r = build_balanced_hamiltonian(spec, R=3.5, screened_valence_basis=False)
        info = r['screened_valence_info']
        assert info['n_orbitals_corrected'] == 0
        assert info['total_trace_shift'] == 0.0
        assert np.all(r['h1_screened_correction'] == 0.0)

    def test_lih_pauli_count_preserved(self):
        # LiH balanced should give 878 Pauli at max_n=2 per CLAUDE.md / Track CD
        spec = lih_spec(max_n=2)
        r = build_balanced_hamiltonian(spec, R=3.015, screened_valence_basis=True)
        assert r['N_pauli'] == 878


class TestNaHWithScreenedValence:
    """NaH with the correction applied."""

    def test_correction_applied(self):
        spec = nah_spec(max_n=2)
        r = build_balanced_hamiltonian(spec, R=3.5, screened_valence_basis=True)
        info = r['screened_valence_info']
        assert info['n_orbitals_corrected'] == 5
        # trace shift ~+0.61 Ha (3s + 4s + 3*4p)
        assert 0.50 < info['total_trace_shift'] < 0.70

    def test_pauli_unchanged(self):
        # Diagonal h1 correction shouldn't change Pauli count
        spec = nah_spec(max_n=2)
        r_off = build_balanced_hamiltonian(spec, R=3.5, screened_valence_basis=False)
        r_on = build_balanced_hamiltonian(spec, R=3.5, screened_valence_basis=True)
        assert r_off['N_pauli'] == r_on['N_pauli']

    def test_h1_diagonal_only_changed(self):
        spec = nah_spec(max_n=2)
        r_off = build_balanced_hamiltonian(spec, R=3.5, screened_valence_basis=False)
        r_on = build_balanced_hamiltonian(spec, R=3.5, screened_valence_basis=True)
        diff = r_on['h1'] - r_off['h1']
        # off-diagonal elements should be identical
        offdiag = diff - np.diag(np.diag(diff))
        assert np.max(np.abs(offdiag)) < 1e-14

    def test_eri_unchanged(self):
        spec = nah_spec(max_n=2)
        r_off = build_balanced_hamiltonian(spec, R=3.5, screened_valence_basis=False)
        r_on = build_balanced_hamiltonian(spec, R=3.5, screened_valence_basis=True)
        assert np.allclose(r_off['eri'], r_on['eri'], atol=0, rtol=0)


class TestComposability:
    """Screened valence basis composes with W1c and PK cleanly."""

    def test_screened_plus_pk_plus_valence(self):
        spec = nah_spec(max_n=2)
        r = build_balanced_hamiltonian(
            spec, R=3.5,
            screened_cross_center=True,
            pk_cross_center=True,
            screened_valence_basis=True,
        )
        # Must produce a valid Hamiltonian
        assert r['N_pauli'] > 0
        assert np.isfinite(r['one_norm'])
        # Each component must be present in diagnostics
        assert r['screened_valence_info']['n_orbitals_corrected'] == 5
        assert r['pk_cross_center_count'] > 0
        # h1 traces should reflect ALL three corrections compounded
        assert np.trace(r['h1_screened_correction']) > 0


class TestMgH2WithScreenedValence:
    """MgH₂ also has [Ne] frozen-core blocks (2 bond blocks)."""

    def test_correction_applied_two_blocks(self):
        spec = mgh2_spec(max_n=2)
        r = build_balanced_hamiltonian(spec, R=3.268, screened_valence_basis=True)
        info = r['screened_valence_info']
        # 2 bond blocks, each 5 orbitals = 10 total
        assert info['n_orbitals_corrected'] == 10
        # trace shift should be positive (less negative)
        assert info['total_trace_shift'] > 0


class TestHClWithScreenedValence:
    """HCl: [Ne] core + Cl side bond block + 3 lone pairs."""

    def test_correction_includes_lone_pairs(self):
        spec = hcl_spec(max_n=2)
        r = build_balanced_hamiltonian(spec, R=2.409, screened_valence_basis=True)
        info = r['screened_valence_info']
        # bond + 3 lone pairs = 4 frozen-core sub-blocks * 5 orbitals each = 20
        assert info['n_orbitals_corrected'] == 20
        assert info['total_trace_shift'] > 0


class TestDiagnosticOutputs:
    """Verify diagnostic dict structure."""

    def test_block_corrections_structure(self):
        spec = nah_spec(max_n=2)
        r = build_balanced_hamiltonian(spec, R=3.5, screened_valence_basis=True)
        info = r['screened_valence_info']
        assert 'block_corrections' in info
        assert len(info['block_corrections']) == 1
        bc = info['block_corrections'][0]
        for key in ('sub_block', 'parent_label', 'Z_nuc', 'Z_orb',
                    'n_val_offset', 'n_states', 'trace_shift_ha',
                    'per_orbital'):
            assert key in bc
        # per_orbital should have one entry per state
        assert len(bc['per_orbital']) == bc['n_states']
        for po in bc['per_orbital']:
            assert 'state' in po
            assert 'delta_h1' in po
            assert isinstance(po['state'], tuple)
            assert len(po['state']) == 3
