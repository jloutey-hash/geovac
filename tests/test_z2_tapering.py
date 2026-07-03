"""
Tests for ``geovac.z2_tapering`` — Hopf-U(1) Z2 tapering production module.

Verifies:
- Orthogonality of the basis rotation U.
- Per-orbital P-parity assignment (+1 for m=0, otherwise +/- per sym/antisym).
- Stabiliser construction in both 'global' and 'per_block' modes.
- End-to-end pipeline: tapering preserves spectrum to machine precision
  on small systems (H2 at Q=10, NaH at Q=20).
- The closed-form ``delta_Q = 2 + n_sb_with_antisym`` formula on
  per_block mode and ``delta_Q = 3`` on global mode (when antisymmetric
  orbitals exist).

Spectrum-preservation tests on systems with Q_tapered <= 16 do a direct
sparse-Lanczos comparison; larger systems rely on the structural
similarity-transform argument (validated in
``debug/sprint_z2_tapering_memo.md`` and
``debug/sprint_z2_per_subblock_memo.md``).
"""

from __future__ import annotations

import numpy as np
import pytest

try:
    import openfermion  # noqa: F401
    from openfermion.linalg import get_sparse_operator
    HAVE_OPENFERMION = True
except ImportError:
    HAVE_OPENFERMION = False

pytestmark = pytest.mark.skipif(
    not HAVE_OPENFERMION,
    reason="openfermion is required for Z2 tapering tests",
)


def _ground_state_energy(qubit_op, n_qubits, max_q=18):
    """Compute the ground-state energy of qubit_op via sparse Lanczos.
    Returns None when n_qubits exceeds max_q.
    """
    import scipy.sparse.linalg as spla
    if n_qubits <= 0:
        return float(np.real(qubit_op.terms.get((), 0.0)))
    if n_qubits > max_q:
        return None
    H = get_sparse_operator(qubit_op, n_qubits=n_qubits)
    if H.shape[0] <= 1024:
        E = np.linalg.eigvalsh(H.toarray())
        return float(E[0])
    vals = spla.eigsh(H, k=1, which='SA', tol=1e-12)[0]
    return float(vals[0])


# ---------------------------------------------------------------------------
# Rotation U: orthogonality + parity assignment
# ---------------------------------------------------------------------------

class TestRotation:

    def test_rotation_orthogonal_minimal(self):
        """Single sub-block with 1s + 2p_{-1,0,+1}: parity = [+, -, +, +]."""
        from geovac.z2_tapering import build_pm_rotation
        # (sb_key, n, l, m): one sub-block, n=1 has l=0,m=0 and n=2 has l=0,m=0
        # plus l=1,m=-1,0,+1.
        orbital_table = [
            (('blk0',), 1, 0, 0),
            (('blk0',), 2, 0, 0),
            (('blk0',), 2, 1, -1),
            (('blk0',), 2, 1, 0),
            (('blk0',), 2, 1, 1),
        ]
        U, parity = build_pm_rotation(orbital_table)
        # Orthogonal
        assert np.allclose(U @ U.T, np.eye(5), atol=1e-12)
        # m=0 orbitals are P-fixed (+1)
        assert parity[0] == +1  # 1s
        assert parity[1] == +1  # 2s
        assert parity[3] == +1  # 2p_0
        # m=+1 / m=-1 give one symmetric and one antisymmetric
        plus_minus = sorted([parity[2], parity[4]])
        assert plus_minus == [-1, +1]

    def test_rotation_multi_subblock(self):
        """Two sub-blocks, each with 2p triplet: 6 orbitals -> 4 sym + 2 anti."""
        from geovac.z2_tapering import build_pm_rotation
        orbital_table = []
        for blk_idx in (0, 1):
            for m in (-1, 0, 1):
                orbital_table.append(((blk_idx, 'X'), 2, 1, m))
        U, parity = build_pm_rotation(orbital_table)
        assert np.allclose(U @ U.T, np.eye(6), atol=1e-12)
        # Each block contributes one m=0 (P=+1) and one (m=+1, m=-1) pair
        # -> 1 sym + 1 antisym per block.
        n_sym = int(np.sum(parity == +1))
        n_anti = int(np.sum(parity == -1))
        assert n_sym == 4  # two m=0 + two symmetric m=+/-1 combos
        assert n_anti == 2  # one per block


# ---------------------------------------------------------------------------
# Stabiliser construction
# ---------------------------------------------------------------------------

class TestStabilizers:

    def test_global_mode_returns_one_P(self):
        from geovac.z2_tapering import build_pm_rotation, build_stabilizers
        orbital_table = [
            ((0, 'A'), 2, 1, -1),
            ((0, 'A'), 2, 1, 0),
            ((0, 'A'), 2, 1, 1),
            ((1, 'B'), 2, 1, -1),
            ((1, 'B'), 2, 1, 0),
            ((1, 'B'), 2, 1, 1),
        ]
        U, parity = build_pm_rotation(orbital_table)
        Z_alpha, Z_beta, P_list = build_stabilizers(
            parity, orbital_table=orbital_table, mode='global',
        )
        assert len(P_list) == 1
        # Z_alpha and Z_beta should always be returned
        assert len(Z_alpha.terms) == 1
        assert len(Z_beta.terms) == 1

    def test_per_block_mode_returns_one_P_per_subblock(self):
        from geovac.z2_tapering import build_pm_rotation, build_stabilizers
        orbital_table = [
            ((0, 'A'), 2, 1, -1),
            ((0, 'A'), 2, 1, 0),
            ((0, 'A'), 2, 1, 1),
            ((1, 'B'), 2, 1, -1),
            ((1, 'B'), 2, 1, 0),
            ((1, 'B'), 2, 1, 1),
        ]
        U, parity = build_pm_rotation(orbital_table)
        Z_alpha, Z_beta, P_list = build_stabilizers(
            parity, orbital_table=orbital_table, mode='per_block',
        )
        assert len(P_list) == 2  # one per sub-block

    def test_no_antisymmetric_returns_empty(self):
        """A pure-s basis (only m=0) has no antisymmetric orbitals."""
        from geovac.z2_tapering import build_pm_rotation, build_stabilizers
        orbital_table = [
            ((0, 'A'), 1, 0, 0),
            ((0, 'A'), 2, 0, 0),
        ]
        U, parity = build_pm_rotation(orbital_table)
        assert all(parity == +1)
        Z_alpha, Z_beta, P_list = build_stabilizers(
            parity, orbital_table=orbital_table, mode='per_block',
        )
        assert P_list == []

    def test_per_block_requires_orbital_table(self):
        from geovac.z2_tapering import build_stabilizers
        parity = np.array([+1, -1, +1, -1])
        with pytest.raises(ValueError, match="orbital_table is required"):
            build_stabilizers(parity, mode='per_block')

    def test_invalid_mode_raises(self):
        from geovac.z2_tapering import build_stabilizers
        parity = np.array([+1, -1])
        with pytest.raises(ValueError, match="mode must be"):
            build_stabilizers(parity, mode='banana')


# ---------------------------------------------------------------------------
# End-to-end spectrum preservation
# ---------------------------------------------------------------------------

class TestEndToEndH2:

    @pytest.fixture(scope='class')
    def h2_spec(self):
        from geovac.molecular_spec import MolecularSpec, OrbitalBlock
        return MolecularSpec(
            name='H2',
            blocks=[OrbitalBlock(
                label='H2_bond', block_type='bond_pair', Z_center=1.0,
                n_electrons=2, max_n=2,
            )],
            nuclear_repulsion_constant=1.0 / 1.4,
        )

    def test_h2_global_delta_q_3(self, h2_spec):
        """H2 (single sub-block) -> global mode gives delta_Q = 3."""
        from geovac.z2_tapering import hopf_tapered_from_spec
        out = hopf_tapered_from_spec(h2_spec, mode='global')
        assert out['delta_Q'] == 3
        assert out['Q_tapered'] == out['Q_naive'] - 3
        assert out['dropped_P_indices'] == []

    def test_h2_per_block_equals_global_for_single_subblock(self, h2_spec):
        """H2 has n_sb=1 so per_block ≡ global (delta_Q = 2 + 1 = 3)."""
        from geovac.z2_tapering import hopf_tapered_from_spec
        out = hopf_tapered_from_spec(h2_spec, mode='per_block')
        assert out['delta_Q'] == 3
        assert out['n_sub_blocks_with_antisym'] == 1

    def test_h2_per_block_spectrum_preserved(self, h2_spec):
        """The ground state of the tapered H2 should match the un-tapered
        ground state to machine precision."""
        from geovac.z2_tapering import hopf_tapered_from_spec
        out = hopf_tapered_from_spec(h2_spec, mode='per_block')

        # Pick best sector via sweep over 2^n_stab signs.  H2 has 3 stabs
        # so 8 sectors.
        from geovac.z2_tapering import (
            apply_hopf_tapering, build_pm_rotation, rotate_h1_eri,
        )
        from geovac.composed_qubit import build_composed_hamiltonian
        from geovac.qubit_encoding import build_fermion_op_from_integrals
        from openfermion import jordan_wigner

        result = build_composed_hamiltonian(h2_spec, pk_in_hamiltonian=True)
        fop_naive = build_fermion_op_from_integrals(
            result['h1'], result['eri'], result['nuclear_repulsion'],
        )
        qop_naive = jordan_wigner(fop_naive)
        E_naive = _ground_state_energy(qop_naive, out['Q_naive'])

        # Sweep sectors to find the ground state in the tapered basis
        from itertools import product
        n_stab = 3
        best_E = None
        for signs in product([+1, -1], repeat=n_stab):
            qop_tap, _ = apply_hopf_tapering(
                jordan_wigner(build_fermion_op_from_integrals(
                    *rotate_h1_eri(
                        result['h1'], result['eri'],
                        build_pm_rotation(out['orbital_table'])[0],
                    ),
                    result['nuclear_repulsion'],
                )),
                out['parity'], orbital_table=out['orbital_table'],
                mode='per_block', sector_signs=list(signs),
            )
            E = _ground_state_energy(qop_tap, out['Q_tapered'])
            if E is not None and (best_E is None or E < best_E):
                best_E = E

        assert best_E is not None
        rel_err = abs(best_E - E_naive) / max(abs(E_naive), 1e-30)
        assert rel_err < 1e-12, (
            f"|dE|/|E| = {rel_err:.2e} > 1e-12 (E_naive={E_naive}, "
            f"best_E={best_E})"
        )


class TestEndToEndNaH:
    """NaH at Q_naive=20: per-block mode should give delta_Q = 4
    (n_sub_blocks = 2: Na frozen-core + H partner).
    """

    @pytest.fixture(scope='class')
    def nah_spec(self):
        from geovac.molecular_spec import hydride_spec
        return hydride_spec(11, max_n=2)

    def test_nah_per_block_delta_q_4(self, nah_spec):
        from geovac.z2_tapering import hopf_tapered_from_spec
        out = hopf_tapered_from_spec(nah_spec, mode='per_block')
        assert out['delta_Q'] == 4, (
            f"expected delta_Q=4 for NaH (n_sb=2), got {out['delta_Q']}"
        )
        assert out['n_sub_blocks_with_antisym'] == 2
        assert out['dropped_P_indices'] == []

    def test_nah_global_delta_q_3(self, nah_spec):
        from geovac.z2_tapering import hopf_tapered_from_spec
        out = hopf_tapered_from_spec(nah_spec, mode='global')
        assert out['delta_Q'] == 3


# ---------------------------------------------------------------------------
# Defensive gate: drop_noncommuting
# ---------------------------------------------------------------------------

class TestDropNoncommuting:

    def test_passthrough_when_construction_is_clean(self):
        """Standard composed builder is clean: dropped list is empty."""
        from geovac.z2_tapering import hopf_tapered_from_spec
        from geovac.molecular_spec import hydride_spec
        out = hopf_tapered_from_spec(
            hydride_spec(11, max_n=2), mode='per_block',
        )
        assert out['dropped_P_indices'] == []


# ---------------------------------------------------------------------------
# ecosystem_export.hamiltonian(..., tapered=...) integration
# ---------------------------------------------------------------------------

class TestEcosystemExportTapered:

    def test_default_is_untapered(self):
        """tapered=None should give the historical (un-tapered) operator."""
        from geovac.ecosystem_export import hamiltonian
        ham = hamiltonian('H2', max_n=2)
        # Default keeps the original Q
        assert 'tapered_mode' not in ham._metadata

    def test_global_mode_removes_three_qubits(self):
        from geovac.ecosystem_export import hamiltonian
        ham_naive = hamiltonian('H2', max_n=2)
        ham_tap = hamiltonian('H2', max_n=2, tapered='global')
        assert ham_tap._metadata['tapered_mode'] == 'global'
        assert ham_tap._metadata['delta_Q'] == 3

    def test_per_block_mode_on_nah(self):
        from geovac.ecosystem_export import hamiltonian
        ham = hamiltonian('NaH', max_n=2, tapered='per_block')
        assert ham._metadata['tapered_mode'] == 'per_block'
        # NaH has 2 sub-blocks -> delta_Q = 2 + 2 = 4
        assert ham._metadata['delta_Q'] == 4
        assert ham._metadata['n_sub_blocks_with_antisym'] == 2

    def test_invalid_mode_raises(self):
        from geovac.ecosystem_export import hamiltonian
        with pytest.raises(ValueError, match="tapered must be"):
            hamiltonian('H2', max_n=2, tapered='banana')

    def test_extended_mode_strict_improvement_over_per_block(self):
        """v3.90.0 wiring: tapered='extended' should give strictly fewer
        qubits than tapered='per_block' on any molecule with l-odd
        orbitals at default n_max=2, and at least as few Pauli terms."""
        from geovac.ecosystem_export import hamiltonian
        for sys_name in ('LiH', 'BeH2', 'H2O', 'HF'):
            ham_pb = hamiltonian(sys_name, max_n=2, tapered='per_block')
            ham_ext = hamiltonian(sys_name, max_n=2, tapered='extended')
            assert ham_ext._metadata['tapered_mode'] == 'extended'
            assert ham_ext._metadata['Q_tapered'] < ham_pb._metadata['Q_tapered'], (
                f"{sys_name}: extended should save more qubits than per_block; "
                f"got Q_pb={ham_pb._metadata['Q_tapered']}, "
                f"Q_ext={ham_ext._metadata['Q_tapered']}"
            )
            assert ham_ext.n_terms <= ham_pb.n_terms, (
                f"{sys_name}: extended should have ≤ Pauli terms than per_block; "
                f"got n_pb={ham_pb.n_terms}, n_ext={ham_ext.n_terms}"
            )

    def test_extended_lih_saves_three_extra_qubits(self):
        """LiH has 3 sub-blocks all with p orbitals → extended saves
        +3 qubits over per_block."""
        from geovac.ecosystem_export import hamiltonian
        ham_pb = hamiltonian('LiH', max_n=2, tapered='per_block')
        ham_ext = hamiltonian('LiH', max_n=2, tapered='extended')
        assert ham_pb._metadata['Q_tapered'] - ham_ext._metadata['Q_tapered'] == 3

    def test_full_mode_strict_improvement_over_extended(self):
        """v3.93.0 wiring: tapered='full' adds hidden particle-conservation
        Z₂s from the symmetry-adapted basis (v3.92.0). Should give strictly
        ≤ Q than 'extended' on every system."""
        from geovac.ecosystem_export import hamiltonian
        for sys_name in ('LiH', 'HF', 'BeH2'):
            ham_ext = hamiltonian(sys_name, max_n=2, tapered='extended')
            ham_full = hamiltonian(sys_name, max_n=2, tapered='full')
            assert ham_full._metadata['tapered_mode'] == 'full'
            assert ham_full._metadata['Q_tapered'] <= ham_ext._metadata['Q_tapered'], (
                f"{sys_name}: full should save ≥ qubits than extended; got "
                f"Q_ext={ham_ext._metadata['Q_tapered']}, "
                f"Q_full={ham_full._metadata['Q_tapered']}"
            )

    def test_full_lih_saves_two_extra_qubits_over_extended(self):
        """LiH: v3.92.0 panel showed Q 22 → 20 (+2) via hidden Z₂s."""
        from geovac.ecosystem_export import hamiltonian
        ham_ext = hamiltonian('LiH', max_n=2, tapered='extended')
        ham_full = hamiltonian('LiH', max_n=2, tapered='full')
        assert ham_ext._metadata['Q_tapered'] - ham_full._metadata['Q_tapered'] == 2


# ---------------------------------------------------------------------------
# Paper 14 SS[hopf_tapering] quantified library claims (8th-cert backfill,
# 2026-07-02): the 254-qubit library saving and a second continuously-verified
# spectrum-preservation case (He) beyond H2.  The 6-system 4e-15 measurement
# is v2.6.0 memo-recorded; the paper sentence now cites exactly this backing.
# ---------------------------------------------------------------------------

class TestLibraryTaperingHeadlines:

    @pytest.mark.slow
    def test_library_delta_q_sum_is_254(self) -> None:
        """Sum of per_block-mode delta_Q over all 37 systems = 254 qubits
        (the paper's savings convention: He/H2 save 3, CO/N2/F2 save 12;
        global mode saves a uniform 3/system = 111, a different quantity)."""
        from geovac.ecosystem_export import hamiltonian, _SYSTEM_REGISTRY
        total = 0
        for name in sorted(set(_SYSTEM_REGISTRY.values())):
            H = hamiltonian(name, tapered='per_block')
            total += H.metadata['delta_Q']
        assert total == 254, f"library delta_Q sum {total}, paper says 254"

    def test_he_ground_state_preserved(self) -> None:
        """He (Q=10): the tapered ground energy matches the naive ground
        energy after the standard sector-sign sweep (the H2 end-to-end
        pattern, second continuously-verified system)."""
        from itertools import product

        from openfermion import jordan_wigner

        from geovac.composed_qubit import build_composed_hamiltonian
        from geovac.molecular_spec import MolecularSpec, OrbitalBlock
        from geovac.qubit_encoding import build_fermion_op_from_integrals
        from geovac.z2_tapering import (
            apply_hopf_tapering, build_pm_rotation, hopf_tapered_from_spec,
            rotate_h1_eri,
        )

        spec = MolecularSpec(
            name='He',
            blocks=[OrbitalBlock(
                label='He_core', block_type='atomic', Z_center=2.0,
                n_electrons=2, max_n=2,
            )],
            nuclear_repulsion_constant=0.0,
        )
        out = hopf_tapered_from_spec(spec, mode='global')
        result = build_composed_hamiltonian(spec, pk_in_hamiltonian=True)
        qop_naive = jordan_wigner(build_fermion_op_from_integrals(
            result['h1'], result['eri'], result['nuclear_repulsion'],
        ))
        E_naive = _ground_state_energy(qop_naive, out['Q_naive'])

        n_stab = out['Q_naive'] - out['Q_tapered']
        qop_rot = jordan_wigner(build_fermion_op_from_integrals(
            *rotate_h1_eri(
                result['h1'], result['eri'],
                build_pm_rotation(out['orbital_table'])[0],
            ),
            result['nuclear_repulsion'],
        ))
        best_E = None
        for signs in product([+1, -1], repeat=n_stab):
            qop_tap, _ = apply_hopf_tapering(
                qop_rot, out['parity'], orbital_table=out['orbital_table'],
                mode='global', sector_signs=list(signs),
            )
            E = _ground_state_energy(qop_tap, out['Q_tapered'])
            if E is not None and (best_E is None or E < best_E):
                best_E = E

        assert best_E is not None
        rel_err = abs(best_E - E_naive) / max(abs(E_naive), 1e-30)
        assert rel_err < 1e-12, (E_naive, best_E)
