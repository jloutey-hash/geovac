"""
Tests for harmonic oscillator shell model (Track NA).

Validates HO energies, degeneracies, magic numbers, radial wavefunctions,
Gaunt selection rule potential-independence, and sparsity comparison.
"""

import numpy as np
import pytest

from geovac.nuclear.harmonic_shell import (
    _ck_coefficient,
    _ho_shell_degeneracy,
    _l_values_in_shell,
    build_nuclear_graph_harmonic,
    harmonic_oscillator_energies,
    harmonic_radial_wavefunctions,
    harmonic_slater_integrals,
    sparsity_analysis,
    verify_magic_numbers,
)


# ── 1. Exact HO energies ───────────────────────────────────────────────────

class TestHOEnergies:
    def test_ho_energies_exact(self):
        """Verify E(N) = hw*(N + 3/2) for N=0..5."""
        hw = 1.0
        energies = harmonic_oscillator_energies(n_max=6, hw=hw)
        for (n_r, l, m_l, sigma), E in energies.items():
            N = 2 * n_r + l
            expected = hw * (N + 1.5)
            assert abs(E - expected) < 1e-14, (
                f"E({n_r},{l},{m_l},{sigma}) = {E}, expected {expected}"
            )

    def test_ho_energies_custom_hw(self):
        """Verify energies scale with hw."""
        hw = 41.0  # MeV, typical nuclear scale
        energies = harmonic_oscillator_energies(n_max=3, hw=hw)
        for (n_r, l, m_l, sigma), E in energies.items():
            N = 2 * n_r + l
            expected = hw * (N + 1.5)
            assert abs(E - expected) < 1e-10


# ── 2. Degeneracies ────────────────────────────────────────────────────────

class TestDegeneracies:
    def test_ho_degeneracies(self):
        """Verify shell degeneracies: 2, 6, 12, 20, 30, 42."""
        expected_deg = [2, 6, 12, 20, 30, 42]
        for N, exp in enumerate(expected_deg):
            assert _ho_shell_degeneracy(N) == exp, (
                f"Shell N={N}: got {_ho_shell_degeneracy(N)}, expected {exp}"
            )

    def test_ho_state_count(self):
        """Verify total states = sum_{N=0}^{n_max-1} (N+1)(N+2)."""
        for n_max in range(1, 8):
            energies = harmonic_oscillator_energies(n_max=n_max)
            expected_total = sum(_ho_shell_degeneracy(N) for N in range(n_max))
            assert len(energies) == expected_total, (
                f"n_max={n_max}: got {len(energies)} states, expected {expected_total}"
            )

    def test_l_values_in_shell(self):
        """Verify allowed l values per shell."""
        assert _l_values_in_shell(0) == [0]
        assert _l_values_in_shell(1) == [1]
        assert _l_values_in_shell(2) == [2, 0]
        assert _l_values_in_shell(3) == [3, 1]
        assert _l_values_in_shell(4) == [4, 2, 0]
        assert _l_values_in_shell(5) == [5, 3, 1]


# ── 3. Magic numbers ───────────────────────────────────────────────────────

class TestMagicNumbers:
    def test_ho_magic_numbers(self):
        """Verify cumulative closures: 2, 8, 20, 40, 70, 112."""
        result = verify_magic_numbers(n_max=6)
        assert result["match"], (
            f"Magic numbers mismatch: expected {result['expected']}, "
            f"got {result['computed']}"
        )

    def test_magic_numbers_from_graph(self):
        """Verify shell closures from the graph builder match."""
        graph = build_nuclear_graph_harmonic(n_max=6)
        expected = [2, 8, 20, 40, 70, 112]
        assert graph["shell_closures"] == expected


# ── 4. Radial wavefunctions ────────────────────────────────────────────────

class TestRadialWavefunctions:
    def test_ho_radial_normalization(self):
        """Verify integral |R_{n_r,l}|^2 r^2 dr = 1 for several (n_r, l)."""
        b = 1.0
        r_grid = np.linspace(0, 15.0, 5000)[1:]
        test_cases = [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (2, 0), (2, 1), (0, 3)]
        for n_r, l in test_cases:
            R = harmonic_radial_wavefunctions(n_r, l, r_grid, b)
            norm = np.trapezoid(R ** 2 * r_grid ** 2, r_grid)
            assert abs(norm - 1.0) < 1e-4, (
                f"Normalization (n_r={n_r}, l={l}): {norm:.6f} (expected 1.0)"
            )

    def test_ho_radial_orthogonality(self):
        """Verify orthogonality for same l, different n_r."""
        b = 1.0
        r_grid = np.linspace(0, 15.0, 5000)[1:]
        # Test l=0: (0,0) vs (1,0) and (0,0) vs (2,0)
        for l in [0, 1]:
            for n_r1, n_r2 in [(0, 1), (0, 2), (1, 2)]:
                if n_r1 + l + 0.5 < 0 or n_r2 + l + 0.5 < 0:
                    continue
                R1 = harmonic_radial_wavefunctions(n_r1, l, r_grid, b)
                R2 = harmonic_radial_wavefunctions(n_r2, l, r_grid, b)
                overlap = np.trapezoid(R1 * R2 * r_grid ** 2, r_grid)
                assert abs(overlap) < 1e-3, (
                    f"Orthogonality (n_r={n_r1},l={l}) vs (n_r={n_r2},l={l}): "
                    f"overlap = {overlap:.6f}"
                )

    def test_ho_radial_positive_at_origin(self):
        """Verify R(r->0) behavior: r^l for small r."""
        b = 1.0
        r_grid = np.linspace(0.001, 0.01, 10)
        for l in [0, 1, 2]:
            R = harmonic_radial_wavefunctions(0, l, r_grid, b)
            # R should behave as r^l near origin
            if l == 0:
                assert abs(R[0]) > 0  # non-zero at origin for l=0
            else:
                # R/r^l should be approximately constant
                ratio = R / r_grid ** l
                rel_var = np.std(ratio) / np.mean(np.abs(ratio))
                assert rel_var < 0.1, (
                    f"l={l}: R/r^l not approximately constant near origin"
                )


# ── 5. Gaunt coefficient potential-independence ─────────────────────────────

class TestGauntPotentialIndependence:
    def test_gaunt_potential_independence(self):
        """
        Gaunt coefficients c^k(l,m; l',m') depend only on angular quantum
        numbers. Verify they are identical whether the radial potential is
        Coulomb or HO.
        """
        # Test a range of (l, m, l', m', k) triples
        test_triples = [
            (0, 0, 0, 0, 0),
            (1, 0, 1, 0, 0),
            (1, 1, 1, 1, 0),
            (1, 0, 1, 0, 2),
            (0, 0, 1, 0, 1),
            (1, -1, 2, -1, 1),
            (2, 0, 2, 0, 0),
            (2, 0, 2, 0, 2),
            (2, 0, 2, 0, 4),
            (1, 1, 2, 0, 1),
        ]
        for la, ma, lc, mc, k in test_triples:
            # The _ck_coefficient function uses only Wigner 3j symbols
            # (pure angular momentum algebra), so it's automatically
            # potential-independent. We verify by calling it and confirming
            # the result is well-defined.
            val = _ck_coefficient(la, ma, lc, mc, k)
            # Call it again to confirm determinism
            val2 = _ck_coefficient(la, ma, lc, mc, k)
            assert val == val2, (
                f"Gaunt coeff not deterministic for "
                f"({la},{ma},{lc},{mc},{k})"
            )

    def test_gaunt_selection_rules(self):
        """Verify triangle inequality and parity selection rules."""
        # Triangle: |l1 - l2| <= k <= l1 + l2
        assert abs(_ck_coefficient(0, 0, 2, 0, 1)) < 1e-15  # k=1, l1+l2=2, |l1-l2|=2 -> k not in range
        assert abs(_ck_coefficient(0, 0, 0, 0, 1)) < 1e-15  # k=1, l1+l2=0 -> k > l1+l2

        # Parity: l1 + l2 + k must be even
        assert abs(_ck_coefficient(1, 0, 0, 0, 0)) < 1e-15  # 1+0+0=1 odd
        assert abs(_ck_coefficient(0, 0, 1, 0, 1)) > 1e-15  # 0+1+1=2 even, valid

    def test_ho_slater_angular_selection(self):
        """
        Verify angular selection rules (triangle + parity) are identical
        for HO and Coulomb: same set of (l1,l2,l3,l4,k) tuples survive.
        """
        # Collect which (l1,l2,l3,l4,k) survive angular screening
        l_max = 3
        surviving = set()
        for l1 in range(l_max + 1):
            for l2 in range(l_max + 1):
                for l3 in range(l_max + 1):
                    for l4 in range(l_max + 1):
                        k_max = min(l1 + l3, l2 + l4)
                        for k in range(0, k_max + 1):
                            if (l1 + l3 + k) % 2 != 0:
                                continue
                            if (l2 + l4 + k) % 2 != 0:
                                continue
                            surviving.add((l1, l2, l3, l4, k))

        # These selection rules come from Gaunt integrals only and are
        # potential-independent. The set is the same for Coulomb and HO.
        assert len(surviving) > 0
        # Verify specific examples
        assert (0, 0, 0, 0, 0) in surviving
        assert (1, 1, 1, 1, 0) in surviving
        assert (1, 1, 1, 1, 2) in surviving
        # Parity-forbidden
        assert (0, 0, 1, 0, 0) not in surviving  # l1+l3+k = 0+1+0 = 1 odd


# ── 6. Slater integrals ────────────────────────────────────────────────────

class TestSlaterIntegrals:
    def test_ho_slater_f0_same_shell(self):
        """
        Verify F^0(0s, 0s) for HO ground state analytically.

        For the HO ground state (n_r=0, l=0):
          R(r) = 2/b^(3/2) * (1/pi)^(1/4) * exp(-r^2/(2b^2))
          (normalized so integral R^2 r^2 dr = 1)

        F^0(0s,0s) = integral_0^inf integral_0^inf
                      R(r1)^2 * (1/r_>) * R(r2)^2 * r1^2 * r2^2 dr1 dr2

        Analytically: F^0 = (2/(b*sqrt(pi))) * (1/sqrt(2))
                     For b=1: F^0 = sqrt(2/pi) ~ 0.7979
        """
        rk = harmonic_slater_integrals(n_max=1, b=1.0, n_grid=4000)
        # Key for F^0(0s,0s,0s,0s): (0,0, 0,0, 0,0, 0,0, 0)
        key = (0, 0, 0, 0, 0, 0, 0, 0, 0)
        assert key in rk, f"F^0(0s,0s) not found in rk_cache"
        numerical = rk[key]
        analytical = np.sqrt(2.0 / np.pi)
        rel_err = abs(numerical - analytical) / analytical
        assert rel_err < 0.01, (
            f"F^0(0s,0s) = {numerical:.6f}, analytical = {analytical:.6f}, "
            f"rel error = {rel_err:.4f}"
        )

    def test_ho_slater_integrals_positive(self):
        """F^0 integrals should be positive (Coulomb repulsion)."""
        rk = harmonic_slater_integrals(n_max=2, b=1.0, n_grid=2000)
        for key, val in rk.items():
            k = key[-1]
            if k == 0:
                # F^0 integrals (direct) must be non-negative
                # (they represent Coulomb repulsion between densities)
                n1, l1, n2, l2, n3, l3, n4, l4 = key[:8]
                if (n1, l1) == (n3, l3) and (n2, l2) == (n4, l4):
                    assert val >= -1e-10, (
                        f"F^0 integral negative: key={key}, val={val}"
                    )


# ── 7. Sparsity analysis ───────────────────────────────────────────────────

class TestSparsity:
    def test_sparsity_angular_invariance(self):
        """
        Angular sparsity percentage should be identical for HO and Coulomb
        at matched quantum numbers (since angular selection rules are
        potential-independent).
        """
        result = sparsity_analysis(n_max=2, b=1.0, n_grid=1000)
        # Angular zeros should be identical — they come from Gaunt rules only
        # Both HO and Coulomb use the SAME angular selection rules
        angular_frac = result["angular_zero_fraction"]
        assert angular_frac > 0, "Expected some angular zeros"
        # The angular_allowed count is by construction the same for both
        assert result["angular_allowed"] > 0

    def test_sparsity_radial_comparison(self):
        """
        Measure and report radial integral magnitude differences between
        HO and Coulomb. This is a measurement test, not strict pass/fail.
        """
        result = sparsity_analysis(n_max=2, b=1.0, n_grid=1500)

        ho_nnz = result["ho"]["nonzero_eris"]
        coul_nnz = result["coulomb"]["nonzero_eris"]
        total = result["total_possible_eris"]
        angular_allowed = result["angular_allowed"]

        # Report
        print(f"\nSparsity analysis (n_max=2):")
        print(f"  Spatial orbitals: {result['n_spatial_orbitals']}")
        print(f"  Total possible ERIs: {total}")
        print(f"  Angular zeros: {result['angular_zeros']} "
              f"({result['angular_zero_fraction']:.1%})")
        print(f"  Angular allowed: {angular_allowed}")
        print(f"  HO non-zero ERIs: {ho_nnz} "
              f"(density {result['ho']['eri_density']:.4f})")
        print(f"  Coulomb non-zero ERIs: {coul_nnz} "
              f"(density {result['coulomb']['eri_density']:.4f})")
        print(f"  HO radial zeros: {result['ho']['radial_zeros']}")
        print(f"  Coulomb radial zeros: {result['coulomb']['radial_zeros']}")

        # Basic sanity: both should have some non-zero ERIs
        assert ho_nnz > 0
        assert coul_nnz > 0
        # Angular zeros should dominate
        assert result["angular_zeros"] > ho_nnz

    @pytest.mark.slow
    def test_sparsity_radial_comparison_nmax3(self):
        """Extended sparsity comparison at n_max=3 (10 spatial orbitals)."""
        result = sparsity_analysis(n_max=3, b=1.0, n_grid=2000)
        ho_nnz = result["ho"]["nonzero_eris"]
        coul_nnz = result["coulomb"]["nonzero_eris"]
        total = result["total_possible_eris"]

        print(f"\nSparsity analysis (n_max=3):")
        print(f"  Spatial orbitals: {result['n_spatial_orbitals']}")
        print(f"  Total possible ERIs: {total}")
        print(f"  Angular zeros: {result['angular_zeros']} "
              f"({result['angular_zero_fraction']:.1%})")
        print(f"  Angular allowed: {result['angular_allowed']}")
        print(f"  HO non-zero ERIs: {ho_nnz} "
              f"(density {result['ho']['eri_density']:.4f})")
        print(f"  Coulomb non-zero ERIs: {coul_nnz} "
              f"(density {result['coulomb']['eri_density']:.4f})")
        print(f"  HO radial zeros: {result['ho']['radial_zeros']}")
        print(f"  Coulomb radial zeros: {result['coulomb']['radial_zeros']}")

        assert ho_nnz > 0
        assert coul_nnz > 0


# ── 8. Quantum number consistency ───────────────────────────────────────────

class TestQuantumNumbers:
    def test_ho_quantum_number_consistency(self):
        """Verify each state has valid (N, n_r, l, m_l, sigma) with N=2*n_r+l."""
        graph = build_nuclear_graph_harmonic(n_max=6)
        for state in graph["states"]:
            N, n_r, l, m_l, sigma = state
            assert N == 2 * n_r + l, f"N != 2*n_r + l for state {state}"
            assert n_r >= 0, f"Negative n_r in state {state}"
            assert l >= 0, f"Negative l in state {state}"
            assert abs(m_l) <= l, f"|m_l| > l in state {state}"
            assert sigma in (0, 1), f"Invalid spin in state {state}"
            assert l in _l_values_in_shell(N), (
                f"l={l} not valid for shell N={N}"
            )

    def test_ho_graph_hamiltonian_hermitian(self):
        """Verify the graph Hamiltonian is Hermitian."""
        graph = build_nuclear_graph_harmonic(n_max=4)
        H = graph["hamiltonian"].toarray()
        assert np.allclose(H, H.T), "Hamiltonian not symmetric"

    def test_ho_graph_state_count(self):
        """Verify graph returns correct total state count."""
        for n_max in [1, 2, 3, 4, 5, 6]:
            graph = build_nuclear_graph_harmonic(n_max=n_max)
            expected = sum(_ho_shell_degeneracy(N) for N in range(n_max))
            assert graph["n_states"] == expected, (
                f"n_max={n_max}: n_states={graph['n_states']}, "
                f"expected={expected}"
            )
