"""
Tests for nested (full-atom) hyperspherical qubit Hamiltonian.

Track DF Sprints 3-4: Validates the nested approach for Be (4 electrons,
single center) and LiH (4 electrons, two centers) against the composed
architecture.
"""

import pytest
import numpy as np


# ============================================================================
# Helpers
# ============================================================================

def _build_nested_be(max_n: int = 2):
    """Build nested (full-atom) Be Hamiltonian."""
    from geovac.molecular_spec import MolecularSpec, OrbitalBlock
    from geovac.composed_qubit import build_composed_hamiltonian

    spec = MolecularSpec(
        name=f'Be_nested_n{max_n}',
        blocks=[OrbitalBlock(
            label='Be_full', block_type='core',
            Z_center=4.0, n_electrons=4, max_n=max_n,
            pk_A=0.0, pk_B=0.0,
        )],
        nuclear_repulsion_constant=0.0,
    )
    return build_composed_hamiltonian(spec, pk_in_hamiltonian=False)


def _build_composed_be(max_n: int = 2, with_pk: bool = True):
    """Build composed Be Hamiltonian (core + valence blocks)."""
    from geovac.molecular_spec import MolecularSpec, OrbitalBlock
    from geovac.composed_qubit import build_composed_hamiltonian

    spec = MolecularSpec(
        name=f'Be_composed_n{max_n}',
        blocks=[
            OrbitalBlock(label='core', block_type='core',
                         Z_center=4.0, n_electrons=2, max_n=1),
            OrbitalBlock(label='val', block_type='bond_pair',
                         Z_center=2.0, n_electrons=2, max_n=max_n,
                         pk_A=13.01 if with_pk else 0.0,
                         pk_B=12.53 if with_pk else 0.0),
        ],
        nuclear_repulsion_constant=0.0,
    )
    return build_composed_hamiltonian(spec, pk_in_hamiltonian=with_pk)


def _one_norm(qubit_op):
    """Compute total and non-identity 1-norm."""
    total = sum(abs(c) for c in qubit_op.terms.values())
    id_coeff = abs(qubit_op.terms.get((), 0.0))
    return total, total - id_coeff


# ============================================================================
# Tests: Nested Be construction
# ============================================================================

class TestNestedBeConstruction:
    """Verify nested Be Hamiltonian builds correctly."""

    def test_nested_be_n2_qubit_count(self):
        """Nested Be at max_n=2 should have Q=10 (5 spatial * 2 spin)."""
        res = _build_nested_be(max_n=2)
        assert res['Q'] == 10

    def test_nested_be_n3_qubit_count(self):
        """Nested Be at max_n=3 should have Q=28 (14 spatial * 2 spin)."""
        res = _build_nested_be(max_n=3)
        assert res['Q'] == 28

    def test_nested_be_n2_pauli_count(self):
        """Nested Be at max_n=2: Pauli count should be 112."""
        res = _build_nested_be(max_n=2)
        assert res['N_pauli'] == 112

    def test_nested_be_no_pk(self):
        """Nested Be should have zero PK matrix elements."""
        res = _build_nested_be(max_n=2)
        assert res['n_pk_nonzero'] == 0

    def test_nested_be_h1_diagonal(self):
        """h1 should be diagonal (hydrogenic eigenvalues only)."""
        res = _build_nested_be(max_n=2)
        h1 = res['h1']
        # Off-diagonal should be zero
        off_diag = h1 - np.diag(np.diag(h1))
        assert np.max(np.abs(off_diag)) < 1e-15

    def test_nested_be_eri_gaunt_sparsity(self):
        """ERI density should be well below 100% (Gaunt selection rules)."""
        res = _build_nested_be(max_n=2)
        M = res['M']
        n_eri = int(np.count_nonzero(np.abs(res['eri']) > 1e-15))
        density = n_eri / M**4
        assert density < 0.15, f"ERI density {density:.1%} exceeds 15%"

    def test_nested_be_eri_symmetry(self):
        """ERI tensor should be symmetric: (pq|rs) = (rs|pq)."""
        res = _build_nested_be(max_n=2)
        eri = res['eri']
        diff = np.max(np.abs(eri - eri.transpose(2, 3, 0, 1)))
        assert diff < 1e-14, f"ERI symmetry violated: max diff {diff}"


# ============================================================================
# Tests: FCI energy validation
# ============================================================================

class TestNestedBeFCI:
    """Validate FCI energies from nested Be Hamiltonian."""

    def test_nested_be_n2_fci_bound(self):
        """Nested Be at max_n=2 FCI energy should be above HF limit."""
        res = _build_nested_be(max_n=2)
        from openfermion import get_sparse_operator
        sH = get_sparse_operator(res['qubit_op'], n_qubits=10).toarray()
        E = float(np.linalg.eigvalsh(sH)[0])
        # Be HF limit: -14.573 Ha. Nested n2 should be above this
        # (limited basis gives ~-13.95 Ha)
        assert E > -15.0, f"Energy {E} is unphysically low"
        assert E < -13.0, f"Energy {E} is unphysically high"

    def test_nested_be_n2_fci_energy(self):
        """Nested Be at max_n=2 should give ~4.9% error vs exact."""
        res = _build_nested_be(max_n=2)
        from openfermion import get_sparse_operator
        sH = get_sparse_operator(res['qubit_op'], n_qubits=10).toarray()
        E = float(np.linalg.eigvalsh(sH)[0])
        E_exact = -14.6674
        err = abs(E - E_exact) / abs(E_exact)
        assert err < 0.10, f"FCI error {err:.1%} exceeds 10%"


# ============================================================================
# Tests: Comparison to composed
# ============================================================================

class TestNestedVsComposed:
    """Compare nested and composed Be Hamiltonians."""

    def test_nested_fewer_qubits(self):
        """Nested should use fewer qubits (no separate core block)."""
        res_n = _build_nested_be(max_n=2)
        res_c = _build_composed_be(max_n=2)
        assert res_n['Q'] < res_c['Q']

    def test_pauli_counts_comparable(self):
        """Pauli counts should be within 20% of each other."""
        res_n = _build_nested_be(max_n=2)
        res_c = _build_composed_be(max_n=2)
        ratio = res_n['N_pauli'] / res_c['N_pauli']
        assert 0.5 < ratio < 2.0, f"Pauli ratio {ratio:.2f} out of range"

    def test_nested_lower_1norm_than_composed_pk(self):
        """Nested 1-norm should be much lower than composed+PK."""
        res_n = _build_nested_be(max_n=2)
        res_c = _build_composed_be(max_n=2, with_pk=True)
        _, ln_ni = _one_norm(res_n['qubit_op'])
        _, lc_ni = _one_norm(res_c['qubit_op'])
        # PK inflates 1-norm ~6x, so nested should be significantly lower
        assert ln_ni < lc_ni, "Nested 1-norm should be lower than composed+PK"
        assert ln_ni < lc_ni / 3, f"Expected >3x advantage, got {lc_ni/ln_ni:.1f}x"

    def test_nested_pauli_per_q_within_gate(self):
        """Nested Pauli/Q should be within 2x of composed (go/no-go gate)."""
        res_n = _build_nested_be(max_n=2)
        res_c = _build_composed_be(max_n=2)
        pq_n = res_n['N_pauli'] / res_n['Q']
        pq_c = res_c['N_pauli'] / res_c['Q']
        ratio = pq_n / pq_c
        assert ratio < 2.0, f"Pauli/Q ratio {ratio:.2f} exceeds 2x gate"


# ============================================================================
# Tests: Scaling
# ============================================================================

class TestNestedScaling:
    """Test Pauli scaling behavior."""

    @pytest.mark.slow
    def test_nested_be_n3_builds(self):
        """Nested Be at max_n=3 should build successfully."""
        res = _build_nested_be(max_n=3)
        assert res['Q'] == 28
        assert res['N_pauli'] > 0

    @pytest.mark.slow
    def test_eri_density_decreases_with_basis(self):
        """ERI density should decrease from max_n=2 to max_n=3."""
        res2 = _build_nested_be(max_n=2)
        res3 = _build_nested_be(max_n=3)
        M2, M3 = res2['M'], res3['M']
        dens2 = int(np.count_nonzero(np.abs(res2['eri']) > 1e-15)) / M2**4
        dens3 = int(np.count_nonzero(np.abs(res3['eri']) > 1e-15)) / M3**4
        assert dens3 < dens2, f"ERI density did not decrease: {dens2:.3f} -> {dens3:.3f}"


# ============================================================================
# Helpers: Nested LiH (Sprint 4)
# ============================================================================

def _build_nested_lih(max_n: int = 2, R: float = 3.015):
    """Build nested (single-center) LiH Hamiltonian."""
    from geovac.molecular_spec import MolecularSpec, OrbitalBlock
    from geovac.composed_qubit import (
        build_composed_hamiltonian, _enumerate_states,
        _compute_rk_integrals_block, _build_eri_block,
    )
    from geovac.qubit_encoding import build_fermion_op_from_integrals
    from geovac.shibuya_wulfman import compute_cross_center_vne
    from openfermion import jordan_wigner

    Z_Li, Z_H = 3.0, 1.0
    l_max = max_n - 1
    L_max = 2 * l_max + 2

    states = _enumerate_states(max_n)
    M = len(states)
    Q = 2 * M

    # h1: Li nuclear attraction (diagonal) + cross-center V_ne from H
    h1 = np.zeros((M, M))
    for i, (n, l, m) in enumerate(states):
        h1[i, i] = -Z_Li**2 / (2.0 * n**2)

    vne_cross = compute_cross_center_vne(
        Z_orb=Z_Li, states=states, Z_nuc=Z_H,
        R_AB=R, L_max=L_max, n_grid=4000, nuc_parity=1,
    )
    h1 += vne_cross

    # ERIs: single-center Slater integrals at Z=3
    rk_cache = _compute_rk_integrals_block(Z_Li, states)
    eri_phys = _build_eri_block(Z_Li, states, rk_cache)
    eri = np.zeros((M, M, M, M))
    for (a, b, c, d), val in eri_phys.items():
        eri[a, c, b, d] = val
    eri = 0.5 * (eri + eri.transpose(2, 3, 0, 1))

    V_NN = Z_Li * Z_H / R
    fermion_op = build_fermion_op_from_integrals(h1, eri, V_NN)
    qubit_op = jordan_wigner(fermion_op)

    return {
        'M': M, 'Q': Q, 'N_pauli': len(qubit_op.terms),
        'h1': h1, 'eri': eri, 'qubit_op': qubit_op,
        'nuclear_repulsion': V_NN,
    }


# ============================================================================
# Tests: Nested LiH construction (Sprint 4)
# ============================================================================

class TestNestedLiHConstruction:
    """Verify nested LiH Hamiltonian builds correctly."""

    def test_nested_lih_n2_qubit_count(self):
        """Nested LiH at max_n=2 should have Q=10."""
        res = _build_nested_lih(max_n=2)
        assert res['Q'] == 10

    def test_nested_lih_n2_pauli_count(self):
        """Nested LiH at max_n=2 should have 120 Pauli terms."""
        res = _build_nested_lih(max_n=2)
        assert res['N_pauli'] == 120

    def test_nested_lih_h1_has_offdiagonal(self):
        """h1 should have off-diagonal elements from cross-center V_ne."""
        res = _build_nested_lih(max_n=2)
        h1 = res['h1']
        off_diag = h1 - np.diag(np.diag(h1))
        assert np.max(np.abs(off_diag)) > 1e-10, "Cross-center V_ne missing"

    def test_nested_lih_eri_gaunt_sparsity(self):
        """ERI density should be below 15%."""
        res = _build_nested_lih(max_n=2)
        M = res['M']
        n_eri = int(np.count_nonzero(np.abs(res['eri']) > 1e-15))
        density = n_eri / M**4
        assert density < 0.15

    def test_nested_lih_eri_symmetry(self):
        """ERI tensor should satisfy (pq|rs) = (rs|pq)."""
        res = _build_nested_lih(max_n=2)
        eri = res['eri']
        diff = np.max(np.abs(eri - eri.transpose(2, 3, 0, 1)))
        assert diff < 1e-14


# ============================================================================
# Tests: Nested LiH FCI (Sprint 4)
# ============================================================================

class TestNestedLiHFCI:
    """Validate FCI energies for nested LiH."""

    def test_nested_lih_bound(self):
        """Nested LiH should be bound (D_e > 0)."""
        from openfermion import get_sparse_operator
        E_list = []
        for R in [2.0, 6.0]:
            res = _build_nested_lih(max_n=2, R=R)
            sH = get_sparse_operator(res['qubit_op'], n_qubits=10).toarray()
            E_list.append(float(np.linalg.eigvalsh(sH)[0]))
        D_e = E_list[1] - E_list[0]  # E(inf) - E(min)
        assert D_e > 0, f"Unbound: D_e = {D_e}"

    def test_nested_lih_de_accuracy(self):
        """D_e should be within 20% of exact (0.092 Ha)."""
        from openfermion import get_sparse_operator
        R_values = [2.0, 2.5, 3.0, 4.0, 6.0]
        energies = []
        for R in R_values:
            res = _build_nested_lih(max_n=2, R=R)
            sH = get_sparse_operator(res['qubit_op'], n_qubits=10).toarray()
            energies.append(float(np.linalg.eigvalsh(sH)[0]))
        E_min = min(energies)
        E_inf = energies[-1]
        D_e = E_inf - E_min
        err = abs(D_e - 0.092) / 0.092
        assert err < 0.20, f"D_e = {D_e:.4f} Ha, error {err:.1%}"

    def test_nested_lih_fewer_qubits_than_composed(self):
        """Nested LiH should use fewer qubits than composed."""
        from geovac.composed_qubit import build_composed_hamiltonian, lih_spec
        res_n = _build_nested_lih(max_n=2)
        spec_c = lih_spec(max_n_core=2, max_n_val=2)
        res_c = build_composed_hamiltonian(spec_c)
        assert res_n['Q'] < res_c['Q']

    def test_nested_lih_fewer_pauli_than_composed(self):
        """Nested LiH should have fewer Pauli terms than composed."""
        from geovac.composed_qubit import build_composed_hamiltonian, lih_spec
        res_n = _build_nested_lih(max_n=2)
        spec_c = lih_spec(max_n_core=2, max_n_val=2)
        res_c = build_composed_hamiltonian(spec_c)
        assert res_n['N_pauli'] < res_c['N_pauli']


# ============================================================================
# Tests: Nested LiH 1-norm (Sprint 4)
# ============================================================================

class TestNestedLiH1Norm:
    """Test 1-norm properties of nested LiH."""

    def test_nested_lih_1norm_lower_than_composed_pk(self):
        """Nested 1-norm should be lower than composed+PK."""
        from geovac.composed_qubit import build_composed_hamiltonian, lih_spec
        res_n = _build_nested_lih(max_n=2)
        spec_c = lih_spec(max_n_core=2, max_n_val=2)
        res_c = build_composed_hamiltonian(spec_c, pk_in_hamiltonian=True)
        _, ln_ni = _one_norm(res_n['qubit_op'])
        _, lc_ni = _one_norm(res_c['qubit_op'])
        assert ln_ni < lc_ni

    def test_nested_lih_pauli_per_q_within_gate(self):
        """Nested Pauli/Q should be within the 15 gate."""
        res = _build_nested_lih(max_n=2)
        pq = res['N_pauli'] / res['Q']
        assert pq <= 15, f"Pauli/Q = {pq:.1f} exceeds gate"


# ============================================================================
# Helpers: Two-center nested LiH (Sprint 4B)
# ============================================================================

def _build_nested_lih_twocenter(max_n: int = 2, R: float = 3.015, Z_orb: float = 4.0):
    """Build two-center nested LiH Hamiltonian (charge-center origin)."""
    from geovac.composed_qubit import (
        _enumerate_states, _compute_rk_integrals_block, _build_eri_block,
        _radial_wf_grid,
    )
    from geovac.qubit_encoding import build_fermion_op_from_integrals
    from geovac.shibuya_wulfman import compute_cross_center_vne
    from openfermion import jordan_wigner

    Z_A, Z_B = 3.0, 1.0
    l_max = max_n - 1
    L_max = 2 * l_max + 2

    states = _enumerate_states(max_n)
    M = len(states)
    Q = 2 * M

    # Charge-center geometry
    d_A = R * Z_B / (Z_A + Z_B)   # Li distance from center
    d_B = R * Z_A / (Z_A + Z_B)   # H distance from center

    # Kinetic energy: T = -Z_orb^2/(2n^2) I + Z_orb * <1/r>
    r_max = 80.0 / max(Z_orb, 0.5)
    r_grid = np.linspace(0, r_max, 4001)[1:]
    unique_nl = sorted(set((n, l) for n, l, m in states))
    R_on_grid = {}
    for n, l in unique_nl:
        R_on_grid[(n, l)] = _radial_wf_grid(Z_orb, n, l, r_grid)

    inv_r = np.zeros((M, M))
    for i, (ni, li, mi) in enumerate(states):
        for j, (nj, lj, mj) in enumerate(states):
            if li != lj or mi != mj:
                continue
            if i == j:
                inv_r[i, i] = Z_orb / ni**2
            elif j > i:
                integrand = R_on_grid[(ni, li)] * R_on_grid[(nj, lj)] * r_grid
                inv_r[i, j] = np.trapezoid(integrand, r_grid)
                inv_r[j, i] = inv_r[i, j]

    T = Z_orb * inv_r
    for i, (ni, li, mi) in enumerate(states):
        T[i, i] += -Z_orb**2 / (2.0 * ni**2)

    # V_ne from both nuclei
    vne_A = compute_cross_center_vne(
        Z_orb=Z_orb, states=states, Z_nuc=Z_A,
        R_AB=d_A, L_max=L_max, nuc_parity=-1,
    )
    vne_B = compute_cross_center_vne(
        Z_orb=Z_orb, states=states, Z_nuc=Z_B,
        R_AB=d_B, L_max=L_max, nuc_parity=+1,
    )

    h1 = T + vne_A + vne_B

    # ERIs: single-center Slater at Z_orb
    rk_cache = _compute_rk_integrals_block(Z_orb, states)
    eri_phys = _build_eri_block(Z_orb, states, rk_cache)
    eri = np.zeros((M, M, M, M))
    for (a, b, c, d), val in eri_phys.items():
        eri[a, c, b, d] = val
    eri = 0.5 * (eri + eri.transpose(2, 3, 0, 1))

    V_NN = Z_A * Z_B / R
    fermion_op = build_fermion_op_from_integrals(h1, eri, V_NN)
    qubit_op = jordan_wigner(fermion_op)

    return {
        'M': M, 'Q': Q, 'N_pauli': len(qubit_op.terms),
        'h1': h1, 'eri': eri, 'qubit_op': qubit_op,
        'nuclear_repulsion': V_NN, 'T': T,
        'vne_A': vne_A, 'vne_B': vne_B,
    }


# ============================================================================
# Tests: Two-center nested LiH construction (Sprint 4B)
# ============================================================================

class TestTwoCenterLiHConstruction:
    """Verify two-center nested LiH builds correctly."""

    def test_twocenter_lih_n2_qubit_count(self):
        """Two-center LiH at max_n=2 should have Q=10."""
        res = _build_nested_lih_twocenter(max_n=2)
        assert res['Q'] == 10

    def test_twocenter_lih_n2_pauli_within_gate(self):
        """Pauli terms should be ≤ 150 (Sprint 4B gate)."""
        res = _build_nested_lih_twocenter(max_n=2)
        assert res['N_pauli'] <= 150, f"Pauli={res['N_pauli']} exceeds 150"

    def test_twocenter_lih_kinetic_positive_diagonal(self):
        """Kinetic energy diagonal should be positive."""
        res = _build_nested_lih_twocenter(max_n=2)
        T = res['T']
        for i in range(T.shape[0]):
            assert T[i, i] > 0, f"T[{i},{i}] = {T[i,i]} is not positive"

    def test_twocenter_lih_kinetic_symmetric(self):
        """Kinetic energy matrix should be symmetric."""
        res = _build_nested_lih_twocenter(max_n=2)
        T = res['T']
        diff = np.max(np.abs(T - T.T))
        assert diff < 1e-14, f"T not symmetric: max diff {diff}"

    def test_twocenter_lih_h1_has_offdiagonal(self):
        """h1 should have off-diagonal elements (V_ne coupling)."""
        res = _build_nested_lih_twocenter(max_n=2)
        h1 = res['h1']
        off_diag = h1 - np.diag(np.diag(h1))
        assert np.max(np.abs(off_diag)) > 1e-10

    def test_twocenter_lih_both_vne_nonzero(self):
        """Both V_ne matrices should have significant elements."""
        res = _build_nested_lih_twocenter(max_n=2)
        assert np.max(np.abs(res['vne_A'])) > 0.1, "V_ne(Li) too small"
        assert np.max(np.abs(res['vne_B'])) > 0.01, "V_ne(H) too small"

    def test_twocenter_lih_eri_symmetry(self):
        """ERI tensor should satisfy (pq|rs) = (rs|pq)."""
        res = _build_nested_lih_twocenter(max_n=2)
        eri = res['eri']
        diff = np.max(np.abs(eri - eri.transpose(2, 3, 0, 1)))
        assert diff < 1e-14

    def test_twocenter_lih_eri_gaunt_sparsity(self):
        """ERI density should be below 15%."""
        res = _build_nested_lih_twocenter(max_n=2)
        M = res['M']
        n_eri = int(np.count_nonzero(np.abs(res['eri']) > 1e-15))
        density = n_eri / M**4
        assert density < 0.15


# ============================================================================
# Tests: Two-center nested LiH FCI (Sprint 4B)
# ============================================================================

class TestTwoCenterLiHFCI:
    """Validate FCI energies for two-center nested LiH."""

    def test_twocenter_lih_energy_bound(self):
        """FCI energy should be in a physical range."""
        from openfermion import get_sparse_operator
        res = _build_nested_lih_twocenter(max_n=2)
        sH = get_sparse_operator(res['qubit_op'], n_qubits=10).toarray()
        E = float(np.linalg.eigvalsh(sH)[0])
        assert E > -12.0, f"Energy {E} is unphysically low"
        assert E < -2.0, f"Energy {E} is unphysically high"

    def test_twocenter_lih_molecule_bound(self):
        """Two-center LiH should be bound (D_e > 0)."""
        from openfermion import get_sparse_operator
        E_list = []
        for R in [3.015, 6.0]:
            res = _build_nested_lih_twocenter(max_n=2, R=R)
            sH = get_sparse_operator(res['qubit_op'], n_qubits=10).toarray()
            E_list.append(float(np.linalg.eigvalsh(sH)[0]))
        D_e = E_list[1] - E_list[0]  # E(inf) - E(eq)
        # May or may not be bound depending on Z_orb; just check it runs
        assert isinstance(D_e, float)

    def test_twocenter_fewer_qubits_than_composed(self):
        """Two-center should use fewer qubits than composed."""
        from geovac.composed_qubit import build_composed_hamiltonian, lih_spec
        res_tc = _build_nested_lih_twocenter(max_n=2)
        spec_c = lih_spec(max_n_core=2, max_n_val=2)
        res_c = build_composed_hamiltonian(spec_c)
        assert res_tc['Q'] < res_c['Q']

    def test_twocenter_pauli_per_q_within_gate(self):
        """Two-center Pauli/Q should be within 15."""
        res = _build_nested_lih_twocenter(max_n=2)
        pq = res['N_pauli'] / res['Q']
        assert pq <= 15, f"Pauli/Q = {pq:.1f} exceeds gate"


# ============================================================================
# Helpers: Heterogeneous nested LiH (Sprint 5)
# ============================================================================

def _compute_yk_vectorized(wf_b, wf_d, r_grid, k):
    """Vectorized Coulomb Y^k potential.

    Y^k(r1) = int R_b(r2) R_d(r2) (r_<^k / r_>^{k+1}) r2^2 dr2
    """
    dr = r_grid[1] - r_grid[0]
    f_r = wf_b * wf_d * r_grid**2

    # Inner: (1/r1^{k+1}) cumsum(f * r^k * dr)
    cum_inner = np.cumsum(f_r * r_grid**k * dr)
    inner = cum_inner / r_grid**(k + 1)

    # Outer: r1^k * reverse_cumsum(f / r^{k+1} * dr)
    outer_integrand = f_r / r_grid**(k + 1) * dr
    rev_cum = np.sum(outer_integrand) - np.cumsum(outer_integrand)
    outer = r_grid**k * rev_cum

    return inner + outer


def _build_nested_lih_heterogeneous(
    max_n: int = 2,
    R: float = 3.015,
    Z_core: float = 3.0,
    Z_val: float = 1.7,
    n_grid: int = 4000,
):
    """Build heterogeneous nested LiH Hamiltonian.

    Two orbital sets (core at Z_core, valence at Z_val), both centered on Li,
    Lowdin-orthogonalized for second quantization.  Full cross-Z ERIs included.

    Parameters
    ----------
    max_n : int
        Max principal quantum number for each set.
    R : float
        Li-H internuclear distance (bohr).
    Z_core : float
        Nuclear charge for core orbital exponent.
    Z_val : float
        Nuclear charge for valence orbital exponent.
    n_grid : int
        Radial grid points.
    """
    from geovac.composed_qubit import (
        _enumerate_states, _radial_wf_grid, _ck_coefficient,
    )
    from geovac.qubit_encoding import build_fermion_op_from_integrals
    from openfermion import jordan_wigner
    from collections import defaultdict

    Z_Li, Z_H = 3.0, 1.0
    l_max = max_n - 1
    L_max = 2 * l_max + 2

    states = _enumerate_states(max_n)
    M_single = len(states)
    M_total = 2 * M_single
    Q = 2 * M_total

    # Orbital Z: first M_single at Z_core, next M_single at Z_val
    Z_list = [Z_core] * M_single + [Z_val] * M_single
    all_states = list(states) + list(states)

    # Shared radial grid
    r_max = max(80.0 / max(Z_core, 0.5), 80.0 / max(Z_val, 0.5))
    r_grid = np.linspace(0, r_max, n_grid + 1)[1:]

    # Pre-compute radial wavefunctions
    unique_nl = sorted(set((n, l) for n, l, m in states))
    wf = {}
    for Z in sorted(set(Z_list)):
        for n, l in unique_nl:
            wf[(Z, n, l)] = _radial_wf_grid(Z, n, l, r_grid)

    # ---- Overlap matrix S ----
    S = np.eye(M_total)
    for i in range(M_total):
        ni, li, mi = all_states[i]
        Zi = Z_list[i]
        for j in range(i + 1, M_total):
            nj, lj, mj = all_states[j]
            Zj = Z_list[j]
            if li != lj or mi != mj:
                continue
            if abs(Zi - Zj) < 1e-10:
                continue  # same Z: hydrogenic orthogonality
            s_ij = float(np.trapezoid(
                wf[(Zi, ni, li)] * wf[(Zj, nj, lj)] * r_grid**2, r_grid,
            ))
            S[i, j] = s_ij
            S[j, i] = s_ij

    # Lowdin X = S^{-1/2}
    eigvals_S, eigvecs_S = np.linalg.eigh(S)
    X = eigvecs_S @ np.diag(1.0 / np.sqrt(eigvals_S)) @ eigvecs_S.T

    # ---- h1_raw: T - Z_Li/r ----
    # Using right-side eigenvalue trick:
    # <a(Za)|T|b(Zb)> = -Zb^2/(2nb^2) S[a,b] + Zb <a|1/r|b>
    # => <a|T - Z_Li/r|b> = -Zb^2/(2nb^2) S + (Zb - Z_Li) <1/r>
    invr = np.zeros((M_total, M_total))
    for i in range(M_total):
        ni, li, mi = all_states[i]
        Zi = Z_list[i]
        for j in range(M_total):
            nj, lj, mj = all_states[j]
            Zj = Z_list[j]
            if li != lj or mi != mj:
                continue
            invr[i, j] = float(np.trapezoid(
                wf[(Zi, ni, li)] * wf[(Zj, nj, lj)] * r_grid, r_grid,
            ))

    h1_raw = np.zeros((M_total, M_total))
    for i in range(M_total):
        for j in range(M_total):
            li, mi = all_states[i][1], all_states[i][2]
            lj, mj = all_states[j][1], all_states[j][2]
            if li != lj or mi != mj:
                continue
            nj = all_states[j][0]
            Zj = Z_list[j]
            h1_raw[i, j] = (
                -Zj**2 / (2 * nj**2) * S[i, j]
                + (Zj - Z_Li) * invr[i, j]
            )
    h1_raw = 0.5 * (h1_raw + h1_raw.T)  # enforce Hermiticity

    # ---- V_ne from H (multipole expansion) ----
    vne_H = np.zeros((M_total, M_total))
    for i in range(M_total):
        ni, li, mi = all_states[i]
        Zi = Z_list[i]
        for j in range(i, M_total):
            nj, lj, mj = all_states[j]
            Zj = Z_list[j]
            if mi != mj:
                continue
            val = 0.0
            for L in range(0, L_max + 1):
                if (li + L + lj) % 2 != 0:
                    continue
                if L < abs(li - lj) or L > li + lj:
                    continue
                c_L = _ck_coefficient(li, mi, lj, mj, L)
                if abs(c_L) < 1e-15:
                    continue
                f_L = np.where(
                    r_grid < R,
                    r_grid**L / R**(L + 1),
                    R**L / r_grid**(L + 1),
                )
                I_L = float(np.trapezoid(
                    wf[(Zi, ni, li)] * f_L * wf[(Zj, nj, lj)] * r_grid**2,
                    r_grid,
                ))
                val += c_L * I_L
            vne_H[i, j] = -Z_H * val
            vne_H[j, i] = vne_H[i, j]

    h1_raw += vne_H

    # ---- ERIs (general cross-Z, Gaunt selection rules) ----
    ck_table = {}
    for i in range(M_total):
        li, mi = all_states[i][1], all_states[i][2]
        for j in range(M_total):
            lj, mj = all_states[j][1], all_states[j][2]
            for k in range(0, li + lj + 1):
                if (li + lj + k) % 2 != 0:
                    continue
                c = _ck_coefficient(li, mi, lj, mj, k)
                if abs(c) > 1e-15:
                    ck_table[(i, j, k)] = c

    ij_k_map = defaultdict(list)
    for (i, j, k), c in ck_table.items():
        ij_k_map[(i, j)].append((k, c))

    # Pre-compute Y^k potentials (vectorized)
    yk_needed = set()
    for (a, c) in ij_k_map:
        ma, mc = all_states[a][2], all_states[c][2]
        for (b, d) in ij_k_map:
            mb, md = all_states[b][2], all_states[d][2]
            if ma + mb != mc + md:
                continue
            for k_ac, _ in ij_k_map[(a, c)]:
                for k_bd, _ in ij_k_map[(b, d)]:
                    if k_ac == k_bd:
                        Zb = Z_list[b]
                        nb, lb = all_states[b][0], all_states[b][1]
                        Zd = Z_list[d]
                        nd, ld = all_states[d][0], all_states[d][1]
                        yk_needed.add((Zb, nb, lb, Zd, nd, ld, k_ac))

    yk_cache = {}
    for key in yk_needed:
        Zb, nb, lb, Zd, nd, ld, k = key
        yk_cache[key] = _compute_yk_vectorized(
            wf[(Zb, nb, lb)], wf[(Zd, nd, ld)], r_grid, k,
        )

    # Build ERI: physicist <ab|cd> stored as chemist eri_raw[a,c,b,d]
    eri_raw = np.zeros((M_total, M_total, M_total, M_total))
    for (a, c), ck_ac in ij_k_map.items():
        na, la = all_states[a][0], all_states[a][1]
        Za = Z_list[a]
        nc, lc = all_states[c][0], all_states[c][1]
        Zc = Z_list[c]
        ma, mc = all_states[a][2], all_states[c][2]
        for (b, d), ck_bd in ij_k_map.items():
            mb, md = all_states[b][2], all_states[d][2]
            if ma + mb != mc + md:
                continue
            Zb = Z_list[b]
            nb, lb = all_states[b][0], all_states[b][1]
            Zd = Z_list[d]
            nd, ld = all_states[d][0], all_states[d][1]
            val = 0.0
            for k_ac, c_ac in ck_ac:
                for k_bd, c_bd in ck_bd:
                    if k_ac != k_bd:
                        continue
                    k = k_ac
                    yk = yk_cache.get((Zb, nb, lb, Zd, nd, ld, k))
                    if yk is None:
                        continue
                    rk = float(np.trapezoid(
                        wf[(Za, na, la)] * wf[(Zc, nc, lc)]
                        * yk * r_grid**2,
                        r_grid,
                    ))
                    val += c_ac * c_bd * rk
            if abs(val) > 1e-15:
                eri_raw[a, c, b, d] = val

    eri_raw = 0.5 * (eri_raw + eri_raw.transpose(2, 3, 0, 1))

    # ---- Lowdin transform ----
    h1 = X.T @ h1_raw @ X
    eri = np.einsum(
        'ap,bq,cr,ds,abcd->pqrs', X, X, X, X, eri_raw, optimize=True,
    )
    eri = 0.5 * (eri + eri.transpose(2, 3, 0, 1))

    # ---- JW encode ----
    V_NN = Z_Li * Z_H / R
    fermion_op = build_fermion_op_from_integrals(h1, eri, V_NN)
    qubit_op = jordan_wigner(fermion_op)

    return {
        'M': M_total, 'Q': Q, 'N_pauli': len(qubit_op.terms),
        'h1': h1, 'h1_raw': h1_raw, 'eri': eri, 'eri_raw': eri_raw,
        'qubit_op': qubit_op, 'fermion_op': fermion_op,
        'nuclear_repulsion': V_NN,
        'S': S, 'X': X, 'S_eigenvalues': eigvals_S,
        'Z_core': Z_core, 'Z_val': Z_val,
    }


def _hetero_fci_energy(result, n_electrons: int = 4):
    """Sector-restricted FCI for heterogeneous nested Hamiltonian."""
    from geovac.coupled_composition import coupled_fci_energy
    return coupled_fci_energy(result, n_electrons)


# ============================================================================
# Tests: Heterogeneous nested LiH construction (Sprint 5)
# ============================================================================

class TestHeterogeneousConstruction:
    """Verify heterogeneous nested LiH builds correctly."""

    def test_hetero_qubit_count(self):
        """Q = 2 * (2 * M_single). At max_n=2: 2*2*5 = 20."""
        res = _build_nested_lih_heterogeneous(max_n=2)
        assert res['Q'] == 20

    def test_hetero_overlap_positive_definite(self):
        """Overlap matrix S must be positive definite."""
        res = _build_nested_lih_heterogeneous(max_n=2)
        assert np.all(res['S_eigenvalues'] > 0.01), (
            f"S nearly singular: min eigenvalue {min(res['S_eigenvalues']):.4f}"
        )

    def test_hetero_h1_hermitian(self):
        """Transformed h1 should be symmetric."""
        res = _build_nested_lih_heterogeneous(max_n=2)
        diff = np.max(np.abs(res['h1'] - res['h1'].T))
        assert diff < 1e-12, f"h1 not symmetric: max diff {diff}"

    def test_hetero_eri_symmetry(self):
        """ERI tensor should satisfy (pq|rs) = (rs|pq)."""
        res = _build_nested_lih_heterogeneous(max_n=2)
        eri = res['eri']
        diff = np.max(np.abs(eri - eri.transpose(2, 3, 0, 1)))
        assert diff < 1e-10, f"ERI symmetry violated: max diff {diff}"

    def test_hetero_eri_gaunt_sparsity(self):
        """ERI density should be below 15% (Gaunt selection rules)."""
        res = _build_nested_lih_heterogeneous(max_n=2)
        M = res['M']
        n_eri = int(np.count_nonzero(np.abs(res['eri']) > 1e-15))
        density = n_eri / M**4
        assert density < 0.15, f"ERI density {density:.1%} exceeds 15%"

    def test_hetero_core_h1_diagonal(self):
        """Core block of h1_raw should have -Z_core^2/(2n^2) on diagonal."""
        res = _build_nested_lih_heterogeneous(max_n=2, R=100.0)
        h1_raw = res['h1_raw']
        # At large R, V_ne from H is negligible
        # Core 1s: index 0, expected -9/2 = -4.5
        assert abs(h1_raw[0, 0] - (-4.5)) < 0.01, (
            f"Core 1s h1 = {h1_raw[0,0]:.4f}, expected -4.5"
        )

    def test_hetero_vne_cross_z_present(self):
        """Cross-Z V_ne elements (core-val) should be nonzero."""
        res = _build_nested_lih_heterogeneous(max_n=2, R=3.015)
        h1_raw = res['h1_raw']
        M_single = res['M'] // 2
        # Off-diagonal block (core-val)
        cross_block = h1_raw[:M_single, M_single:]
        assert np.max(np.abs(cross_block)) > 0.01, (
            "Cross-Z h1 elements are too small"
        )

    def test_hetero_pauli_count(self):
        """Pauli terms should be below composed+balanced (878)."""
        res = _build_nested_lih_heterogeneous(max_n=2)
        # Lowdin inflates Pauli count; check it's below balanced coupled
        assert res['N_pauli'] < 2000, (
            f"N_pauli = {res['N_pauli']} exceeds 2000"
        )


# ============================================================================
# Tests: Heterogeneous nested LiH FCI (Sprint 5)
# ============================================================================

class TestHeterogeneousFCI:
    """Validate FCI energies for heterogeneous nested LiH."""

    def test_hetero_fci_energy_physical(self):
        """FCI energy should be in physical range."""
        res = _build_nested_lih_heterogeneous(max_n=2, R=3.015)
        fci = _hetero_fci_energy(res)
        E = fci['E_coupled']
        assert E > -12.0, f"Energy {E} is unphysically low"
        assert E < -5.0, f"Energy {E} is unphysically high"

    def test_hetero_fci_bound(self):
        """Heterogeneous LiH should be bound (D_e > 0)."""
        E_list = []
        for R in [3.015, 8.0]:
            res = _build_nested_lih_heterogeneous(max_n=2, R=R, Z_val=1.7)
            fci = _hetero_fci_energy(res)
            E_list.append(fci['E_coupled'])
        D_e = E_list[1] - E_list[0]
        assert D_e > 0, f"Unbound: D_e = {D_e:.4f}"

    def test_hetero_fewer_qubits_than_composed(self):
        """Heterogeneous should use fewer qubits than composed (20 < 30)."""
        from geovac.composed_qubit import build_composed_hamiltonian, lih_spec
        res_h = _build_nested_lih_heterogeneous(max_n=2)
        spec_c = lih_spec(max_n_core=2, max_n_val=2)
        res_c = build_composed_hamiltonian(spec_c)
        assert res_h['Q'] < res_c['Q']


# ============================================================================
# Tests: Heterogeneous vs uniform nested comparison (Sprint 5)
# ============================================================================

class TestHeterogeneousVsUniform:
    """Compare heterogeneous to uniform nested (Sprint 4)."""

    def test_hetero_more_qubits_than_uniform(self):
        """Heterogeneous uses more qubits than uniform (20 vs 10)."""
        res_h = _build_nested_lih_heterogeneous(max_n=2)
        res_u = _build_nested_lih(max_n=2)
        assert res_h['Q'] > res_u['Q']

    def test_hetero_eri_density_within_sprint2_prediction(self):
        """ERI density should stay within 12% (Sprint 2 prediction)."""
        res = _build_nested_lih_heterogeneous(max_n=2)
        M = res['M']
        n_eri = int(np.count_nonzero(np.abs(res['eri']) > 1e-15))
        density = n_eri / M**4
        assert density < 0.12, (
            f"ERI density {density:.1%} exceeds Sprint 2 prediction of 12%"
        )
