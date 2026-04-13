"""
Tests for new molecule qubit Hamiltonians: HF, NH3, CH4 (Track BJ).

Validates that the general composed builder produces correct qubit
Hamiltonians for hydrogen fluoride, ammonia, and methane using
ab initio PK parameters from Track BI.
"""

import numpy as np
import pytest

from geovac.composed_qubit import build_composed_hamiltonian
from geovac.molecular_spec import hf_spec, nh3_spec, ch4_spec

# Legacy per-molecule builders no longer exist; tests use spec-driven path
build_composed_hf = None
build_composed_nh3 = None
build_composed_ch4 = None


# ---------------------------------------------------------------------------
# Expected orbital counts at max_n=2
# ---------------------------------------------------------------------------
# _enumerate_states(2) gives: 1s(1,0,0), 2s(2,0,0), 2p(2,1,-1),
#   (2,1,0), (2,1,1) = 5 spatial orbitals
#
# HF: 1 core(5) + 1 bond(5+5) + 3 lone(5 each) = 5 + 10 + 15 = 30 spatial
# NH3: 1 core(5) + 3 bonds(5+5 each) + 1 lone(5) = 5 + 30 + 5 = 40 spatial
# CH4: 1 core(5) + 4 bonds(5+5 each) = 5 + 40 = 45 spatial

M_HF_EXPECTED = 30
M_NH3_EXPECTED = 40
M_CH4_EXPECTED = 45

Q_HF_EXPECTED = 2 * M_HF_EXPECTED    # 60
Q_NH3_EXPECTED = 2 * M_NH3_EXPECTED   # 80
Q_CH4_EXPECTED = 2 * M_CH4_EXPECTED   # 90


# ---------------------------------------------------------------------------
# HF tests
# ---------------------------------------------------------------------------

class TestHF:
    """Validate HF (hydrogen fluoride) qubit Hamiltonian."""

    @pytest.fixture(scope='class')
    def result(self):
        return build_composed_hf(max_n_core=2, max_n_val=2, verbose=False)

    @pytest.fixture(scope='class')
    def result_via_general(self):
        spec = hf_spec(max_n_core=2, max_n_val=2)
        return build_composed_hamiltonian(spec, pk_in_hamiltonian=True, verbose=False)

    def test_builds_without_error(self, result):
        assert result is not None

    def test_spatial_orbitals(self, result):
        assert result['M'] == M_HF_EXPECTED

    def test_qubit_count(self, result):
        assert result['Q'] == Q_HF_EXPECTED

    def test_pauli_count_positive(self, result):
        assert result['N_pauli'] > 0
        assert isinstance(result['N_pauli'], int)

    def test_h1_shape(self, result):
        M = result['M']
        assert result['h1'].shape == (M, M)

    def test_eri_shape(self, result):
        M = result['M']
        assert result['eri'].shape == (M, M, M, M)

    def test_nuclear_repulsion_finite(self, result):
        nuc = result['nuclear_repulsion']
        assert np.isfinite(nuc)

    def test_h1_pk_shape(self, result):
        M = result['M']
        assert result['h1_pk'].shape == (M, M)

    def test_block_count(self, result):
        # 5 blocks: core + bond + 3 lone pairs
        assert len(result['blocks']) == 5

    def test_spec_name(self, result):
        assert result['spec_name'] == 'HF'

    def test_cross_block_eri_zero(self, result):
        """Cross-block ERIs must be exactly zero."""
        eri = result['eri']
        blocks = result['blocks']
        M = result['M']
        # Build block ranges
        ranges = []
        for b in blocks:
            start = b['center_offset']
            end = start + b['center_M']
            if b['partner_M'] > 0:
                end = b['partner_offset'] + b['partner_M']
            ranges.append((start, end))
        # Check that cross-block ERI entries are zero
        for i, (s1, e1) in enumerate(ranges):
            for j, (s2, e2) in enumerate(ranges):
                if i == j:
                    continue
                cross = eri[s1:e1, s1:e1, s2:e2, s2:e2]
                assert np.allclose(cross, 0.0), (
                    f"Cross-block ERI nonzero between blocks {i} and {j}"
                )

    def test_wrapper_matches_general(self, result, result_via_general):
        """Wrapper function matches direct general builder call."""
        assert result['Q'] == result_via_general['Q']
        assert result['N_pauli'] == result_via_general['N_pauli']
        np.testing.assert_allclose(result['h1'], result_via_general['h1'], atol=1e-14)

    def test_no_pk_option(self):
        """Builder runs with include_pk=False."""
        r = build_composed_hf(include_pk=False, verbose=False)
        assert r['N_pauli'] > 0
        assert np.allclose(r['h1_pk'], 0.0)


# ---------------------------------------------------------------------------
# NH3 tests
# ---------------------------------------------------------------------------

class TestNH3:
    """Validate NH3 (ammonia) qubit Hamiltonian."""

    @pytest.fixture(scope='class')
    def result(self):
        return build_composed_nh3(max_n_core=2, max_n_val=2, verbose=False)

    @pytest.fixture(scope='class')
    def result_via_general(self):
        spec = nh3_spec(max_n_core=2, max_n_val=2)
        return build_composed_hamiltonian(spec, pk_in_hamiltonian=True, verbose=False)

    def test_builds_without_error(self, result):
        assert result is not None

    def test_spatial_orbitals(self, result):
        assert result['M'] == M_NH3_EXPECTED

    def test_qubit_count(self, result):
        assert result['Q'] == Q_NH3_EXPECTED

    def test_pauli_count_positive(self, result):
        assert result['N_pauli'] > 0
        assert isinstance(result['N_pauli'], int)

    def test_h1_shape(self, result):
        M = result['M']
        assert result['h1'].shape == (M, M)

    def test_eri_shape(self, result):
        M = result['M']
        assert result['eri'].shape == (M, M, M, M)

    def test_nuclear_repulsion_finite(self, result):
        nuc = result['nuclear_repulsion']
        assert np.isfinite(nuc)

    def test_h1_pk_shape(self, result):
        M = result['M']
        assert result['h1_pk'].shape == (M, M)

    def test_block_count(self, result):
        # 5 blocks: core + 3 bonds + 1 lone pair
        assert len(result['blocks']) == 5

    def test_spec_name(self, result):
        assert result['spec_name'] == 'NH3'

    def test_cross_block_eri_zero(self, result):
        """Cross-block ERIs must be exactly zero."""
        eri = result['eri']
        blocks = result['blocks']
        ranges = []
        for b in blocks:
            start = b['center_offset']
            end = start + b['center_M']
            if b['partner_M'] > 0:
                end = b['partner_offset'] + b['partner_M']
            ranges.append((start, end))
        for i, (s1, e1) in enumerate(ranges):
            for j, (s2, e2) in enumerate(ranges):
                if i == j:
                    continue
                cross = eri[s1:e1, s1:e1, s2:e2, s2:e2]
                assert np.allclose(cross, 0.0), (
                    f"Cross-block ERI nonzero between blocks {i} and {j}"
                )

    def test_wrapper_matches_general(self, result, result_via_general):
        assert result['Q'] == result_via_general['Q']
        assert result['N_pauli'] == result_via_general['N_pauli']
        np.testing.assert_allclose(result['h1'], result_via_general['h1'], atol=1e-14)

    def test_no_pk_option(self):
        r = build_composed_nh3(include_pk=False, verbose=False)
        assert r['N_pauli'] > 0
        assert np.allclose(r['h1_pk'], 0.0)


# ---------------------------------------------------------------------------
# CH4 tests
# ---------------------------------------------------------------------------

class TestCH4:
    """Validate CH4 (methane) qubit Hamiltonian."""

    @pytest.fixture(scope='class')
    def result(self):
        return build_composed_ch4(max_n_core=2, max_n_val=2, verbose=False)

    @pytest.fixture(scope='class')
    def result_via_general(self):
        spec = ch4_spec(max_n_core=2, max_n_val=2)
        return build_composed_hamiltonian(spec, pk_in_hamiltonian=True, verbose=False)

    def test_builds_without_error(self, result):
        assert result is not None

    def test_spatial_orbitals(self, result):
        assert result['M'] == M_CH4_EXPECTED

    def test_qubit_count(self, result):
        assert result['Q'] == Q_CH4_EXPECTED

    def test_pauli_count_positive(self, result):
        assert result['N_pauli'] > 0
        assert isinstance(result['N_pauli'], int)

    def test_h1_shape(self, result):
        M = result['M']
        assert result['h1'].shape == (M, M)

    def test_eri_shape(self, result):
        M = result['M']
        assert result['eri'].shape == (M, M, M, M)

    def test_nuclear_repulsion_finite(self, result):
        nuc = result['nuclear_repulsion']
        assert np.isfinite(nuc)

    def test_h1_pk_shape(self, result):
        M = result['M']
        assert result['h1_pk'].shape == (M, M)

    def test_block_count(self, result):
        # 5 blocks: core + 4 bonds
        assert len(result['blocks']) == 5

    def test_spec_name(self, result):
        assert result['spec_name'] == 'CH4'

    def test_cross_block_eri_zero(self, result):
        """Cross-block ERIs must be exactly zero."""
        eri = result['eri']
        blocks = result['blocks']
        ranges = []
        for b in blocks:
            start = b['center_offset']
            end = start + b['center_M']
            if b['partner_M'] > 0:
                end = b['partner_offset'] + b['partner_M']
            ranges.append((start, end))
        for i, (s1, e1) in enumerate(ranges):
            for j, (s2, e2) in enumerate(ranges):
                if i == j:
                    continue
                cross = eri[s1:e1, s1:e1, s2:e2, s2:e2]
                assert np.allclose(cross, 0.0), (
                    f"Cross-block ERI nonzero between blocks {i} and {j}"
                )

    def test_wrapper_matches_general(self, result, result_via_general):
        assert result['Q'] == result_via_general['Q']
        assert result['N_pauli'] == result_via_general['N_pauli']
        np.testing.assert_allclose(result['h1'], result_via_general['h1'], atol=1e-14)

    def test_no_pk_option(self):
        r = build_composed_ch4(include_pk=False, verbose=False)
        assert r['N_pauli'] > 0
        assert np.allclose(r['h1_pk'], 0.0)


# ---------------------------------------------------------------------------
# Electron count consistency
# ---------------------------------------------------------------------------

class TestElectronCounts:
    """Verify electron counts are consistent across all new molecules."""

    def test_hf_electron_count(self):
        spec = hf_spec()
        total_e = sum(b.n_electrons for b in spec.blocks)
        assert total_e == 10  # F(9) + H(1)

    def test_nh3_electron_count(self):
        spec = nh3_spec()
        total_e = sum(b.n_electrons for b in spec.blocks)
        assert total_e == 10  # N(7) + 3*H(1)

    def test_ch4_electron_count(self):
        spec = ch4_spec()
        total_e = sum(b.n_electrons for b in spec.blocks)
        assert total_e == 10  # C(6) + 4*H(1)


# ---------------------------------------------------------------------------
# Scaling consistency (all 10-electron molecules should have similar Q^~2.5)
# ---------------------------------------------------------------------------

class TestScalingConsistency:
    """Verify Pauli counts are consistent with Q^2.5 scaling regime."""

    @pytest.fixture(scope='class')
    def all_results(self):
        """Build all three molecules and collect (Q, N_pauli)."""
        results = {}
        for name, builder in [
            ('HF', lambda: build_composed_hf(verbose=False)),
            ('NH3', lambda: build_composed_nh3(verbose=False)),
            ('CH4', lambda: build_composed_ch4(verbose=False)),
        ]:
            r = builder()
            one_norm = sum(abs(v) for v in r['qubit_op'].terms.values())
            results[name] = {
                'Q': r['Q'], 'N_pauli': r['N_pauli'],
                'one_norm': one_norm, 'wall_time': r['wall_time_s'],
                'nuc_repul': r['nuclear_repulsion'],
                'n_blocks': len(r['blocks']),
                'M': r['M'],
            }
        return results

    def test_pauli_counts_reasonable(self, all_results):
        """Pauli counts should be in a reasonable range for 10-electron molecules."""
        for name, data in all_results.items():
            Q = data['Q']
            N = data['N_pauli']
            # At Q^2.5 scaling, Q=60 -> ~27,885; Q=80 -> ~57,243; Q=90 -> ~76,837
            # Allow wide range since block structure affects exact count
            assert N > 100, f"{name}: N_pauli={N} too small"
            assert N < 200000, f"{name}: N_pauli={N} too large"

    def test_print_benchmark_table(self, all_results):
        """Print benchmark summary (always passes, for data collection)."""
        print("\n" + "=" * 80)
        print("TRACK BJ BENCHMARK TABLE — New Molecules (max_n=2)")
        print("=" * 80)
        print(f"{'System':<8} {'Blocks':>6} {'M':>5} {'Q':>5} "
              f"{'N_pauli':>10} {'1-norm (Ha)':>12} {'V_nuc (Ha)':>12} "
              f"{'Time (s)':>10}")
        print("-" * 80)
        for name in ['HF', 'NH3', 'CH4']:
            d = all_results[name]
            print(f"{name:<8} {d['n_blocks']:>6} {d['M']:>5} {d['Q']:>5} "
                  f"{d['N_pauli']:>10,} {d['one_norm']:>12.2f} "
                  f"{d['nuc_repul']:>12.4f} {d['wall_time']:>10.1f}")
        print("=" * 80)
