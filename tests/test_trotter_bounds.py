"""
Tests for geovac.trotter_bounds — Trotter error bounds and cost metrics.

Validates:
  1. 1-norm correctness for hand-built operators
  2. Trotter step formula gives correct ceiling
  3. Monotonicity: more steps -> smaller error bound
  4. H2 STO-3G vs GeoVac He comparison with printed table
  5. Scaling: lambda/Q for GeoVac He at nmax=2 vs nmax=3

Author: GeoVac Development Team
Date: March 2026
"""

import math
import warnings

import numpy as np
import pytest

from openfermion import QubitOperator

from geovac.trotter_bounds import (
    pauli_1norm,
    max_coefficient,
    coefficient_statistics,
    first_order_trotter_bound,
    trotter_steps_required,
    TrotterAnalysis,
    analyze_trotter_cost,
    compare_trotter_cost,
    _pauli_string_to_binary,
    _symplectic_inner_product_matrix,
    pauli_commutator_bound,
    analyze_commutator_cost,
)
from geovac.gaussian_reference import h2_sto3g, build_qubit_hamiltonian
from geovac.lattice_index import LatticeIndex
from geovac.qubit_encoding import JordanWignerEncoder

warnings.filterwarnings('ignore', category=UserWarning)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def h2_sto3g_qop():
    """H2 STO-3G QubitOperator."""
    _, qop, _ = build_qubit_hamiltonian(h2_sto3g())
    return qop


@pytest.fixture(scope="module")
def geovac_he_nmax2():
    """GeoVac He at nmax=2."""
    return LatticeIndex(
        n_electrons=2, max_n=2, nuclear_charge=2,
        vee_method='slater_full', h1_method='hybrid', fci_method='matrix',
    )


@pytest.fixture(scope="module")
def geovac_he_nmax3():
    """GeoVac He at nmax=3."""
    return LatticeIndex(
        n_electrons=2, max_n=3, nuclear_charge=2,
        vee_method='slater_full', h1_method='hybrid', fci_method='matrix',
    )


# ---------------------------------------------------------------------------
# 1-norm correctness
# ---------------------------------------------------------------------------

class TestPauli1Norm:
    """Tests for pauli_1norm and related coefficient functions."""

    def test_known_coefficients(self) -> None:
        """Hand-built operator: 0.5 X0 + 0.3 Z1 - 0.2 X0Z1."""
        op = (
            QubitOperator(((0, 'X'),), 0.5)
            + QubitOperator(((1, 'Z'),), 0.3)
            + QubitOperator(((0, 'X'), (1, 'Z')), -0.2)
        )
        assert pauli_1norm(op) == pytest.approx(1.0, abs=1e-12)

    def test_identity_only(self) -> None:
        """Pure identity operator: H = 3.14 * I."""
        op = QubitOperator((), 3.14)
        assert pauli_1norm(op) == pytest.approx(3.14, abs=1e-12)

    def test_empty_operator(self) -> None:
        """Empty operator has zero 1-norm."""
        op = QubitOperator()
        assert pauli_1norm(op) == 0.0

    def test_max_coefficient(self) -> None:
        """max_coefficient picks the largest |c_i|."""
        op = (
            QubitOperator(((0, 'X'),), 0.5)
            + QubitOperator(((1, 'Z'),), -0.8)
            + QubitOperator(((0, 'Y'),), 0.3)
        )
        assert max_coefficient(op) == pytest.approx(0.8, abs=1e-12)

    def test_max_coefficient_empty(self) -> None:
        """Empty operator returns 0."""
        assert max_coefficient(QubitOperator()) == 0.0

    def test_coefficient_statistics(self) -> None:
        """Statistics dict has correct keys and values."""
        op = (
            QubitOperator(((0, 'X'),), 0.4)
            + QubitOperator(((1, 'Z'),), 0.6)
        )
        stats = coefficient_statistics(op)
        assert stats['one_norm'] == pytest.approx(1.0, abs=1e-12)
        assert stats['max_coeff'] == pytest.approx(0.6, abs=1e-12)
        assert stats['median_coeff'] == pytest.approx(0.5, abs=1e-12)
        assert stats['mean_coeff'] == pytest.approx(0.5, abs=1e-12)
        assert stats['n_terms'] == 2
        assert stats['n_qubits'] == 2


# ---------------------------------------------------------------------------
# Trotter step formula
# ---------------------------------------------------------------------------

class TestTrotterFormulas:
    """Tests for Trotter error bound and step count formulas."""

    def test_trotter_bound_known_value(self) -> None:
        """
        For lambda=2.0, t=1.0, r=10:
        error <= (1/10)^2 * 4 / 2 = 0.02
        """
        # Build operator with 1-norm = 2.0
        op = (
            QubitOperator(((0, 'X'),), 1.0)
            + QubitOperator(((1, 'Z'),), 1.0)
        )
        assert pauli_1norm(op) == pytest.approx(2.0, abs=1e-12)

        bound = first_order_trotter_bound(op, time=1.0, n_steps=10)
        expected = (1.0 / 10)**2 * 2.0**2 / 2.0  # = 0.02
        assert bound == pytest.approx(expected, abs=1e-12)

    def test_trotter_steps_known_value(self) -> None:
        """
        For lambda=2.0, t=1.0, eps=0.02:
        r >= 1.0 * 2.0 / sqrt(2 * 0.02) = 2.0 / 0.2 = 10.0
        """
        op = (
            QubitOperator(((0, 'X'),), 1.0)
            + QubitOperator(((1, 'Z'),), 1.0)
        )
        r = trotter_steps_required(op, time=1.0, epsilon=0.02)
        assert r == 10

    def test_trotter_steps_ceiling(self) -> None:
        """
        For lambda=2.0, t=1.0, eps=0.019:
        r >= 2.0 / sqrt(0.038) = 2.0 / 0.19494... = 10.259...
        ceil -> 11
        """
        op = (
            QubitOperator(((0, 'X'),), 1.0)
            + QubitOperator(((1, 'Z'),), 1.0)
        )
        r = trotter_steps_required(op, time=1.0, epsilon=0.019)
        assert r == 11

    def test_monotonicity(self) -> None:
        """More Trotter steps -> smaller error bound."""
        op = (
            QubitOperator(((0, 'X'),), 1.5)
            + QubitOperator(((1, 'Z'),), 0.7)
            + QubitOperator(((0, 'Y'), (1, 'X')), -0.3)
        )
        errors = [
            first_order_trotter_bound(op, time=1.0, n_steps=r)
            for r in [1, 2, 5, 10, 50, 100]
        ]
        for i in range(len(errors) - 1):
            assert errors[i] > errors[i + 1], (
                f"Error should decrease: r={[1,2,5,10,50,100][i]} "
                f"gave {errors[i]}, r={[1,2,5,10,50,100][i+1]} gave {errors[i+1]}"
            )

    def test_invalid_n_steps(self) -> None:
        """n_steps < 1 raises ValueError."""
        op = QubitOperator(((0, 'X'),), 1.0)
        with pytest.raises(ValueError, match="n_steps"):
            first_order_trotter_bound(op, time=1.0, n_steps=0)

    def test_invalid_epsilon(self) -> None:
        """epsilon <= 0 raises ValueError."""
        op = QubitOperator(((0, 'X'),), 1.0)
        with pytest.raises(ValueError, match="epsilon"):
            trotter_steps_required(op, time=1.0, epsilon=0.0)


# ---------------------------------------------------------------------------
# TrotterAnalysis
# ---------------------------------------------------------------------------

class TestTrotterAnalysis:
    """Tests for analyze_trotter_cost."""

    def test_h2_sto3g(self, h2_sto3g_qop) -> None:
        """Analysis has consistent fields for H2 STO-3G."""
        analysis = analyze_trotter_cost(h2_sto3g_qop)
        assert analysis.n_terms == 15
        assert analysis.n_qubits == 4
        assert analysis.one_norm > 0
        assert analysis.max_coefficient > 0
        assert analysis.max_coefficient <= analysis.one_norm
        assert analysis.trotter_steps_eps3 > 0
        assert analysis.trotter_steps_eps6 > analysis.trotter_steps_eps3

    def test_summary_string(self, h2_sto3g_qop) -> None:
        """Summary contains key metrics."""
        analysis = analyze_trotter_cost(h2_sto3g_qop)
        s = analysis.summary()
        assert "lambda=" in s
        assert "r(1e-3)=" in s
        assert "r(1e-6)=" in s


# ---------------------------------------------------------------------------
# GeoVac vs Gaussian comparison
# ---------------------------------------------------------------------------

class TestComparison:
    """Trotter cost comparison: GeoVac He vs H2 STO-3G."""

    def test_comparison_table(
        self,
        h2_sto3g_qop,
        geovac_he_nmax2: LatticeIndex,
    ) -> None:
        """Print comparison table and verify basic ordering."""
        enc = JordanWignerEncoder(geovac_he_nmax2)
        geovac_qop = enc.build_qubit_operator()

        table = compare_trotter_cost(
            geovac_qop, h2_sto3g_qop,
            labels=('GeoVac He n2', 'H2 STO-3G'),
        )
        print(table)

        a_geo = analyze_trotter_cost(geovac_qop)
        a_h2 = analyze_trotter_cost(h2_sto3g_qop)

        # Both should have positive 1-norms
        assert a_geo.one_norm > 0
        assert a_h2.one_norm > 0


# ---------------------------------------------------------------------------
# Scaling test: lambda/Q
# ---------------------------------------------------------------------------

class TestScaling:
    """Verify 1-norm scaling with basis size."""

    def test_lambda_per_qubit(
        self,
        geovac_he_nmax2: LatticeIndex,
        geovac_he_nmax3: LatticeIndex,
    ) -> None:
        """
        lambda/Q should grow sub-linearly with Q if sparsity helps.

        For a dense Hamiltonian, lambda ~ Q^2 so lambda/Q ~ Q.
        For a sparse Hamiltonian, lambda/Q should grow more slowly.
        """
        enc2 = JordanWignerEncoder(geovac_he_nmax2)
        enc3 = JordanWignerEncoder(geovac_he_nmax3)
        qop2 = enc2.build_qubit_operator()
        qop3 = enc3.build_qubit_operator()

        a2 = analyze_trotter_cost(qop2)
        a3 = analyze_trotter_cost(qop3)

        lpq2 = a2.one_norm / a2.n_qubits
        lpq3 = a3.one_norm / a3.n_qubits
        q_ratio = a3.n_qubits / a2.n_qubits

        print(f"\n  nmax=2: Q={a2.n_qubits}, lambda={a2.one_norm:.4f}, "
              f"lambda/Q={lpq2:.4f}")
        print(f"  nmax=3: Q={a3.n_qubits}, lambda={a3.one_norm:.4f}, "
              f"lambda/Q={lpq3:.4f}")
        print(f"  Q ratio: {q_ratio:.2f}, lambda/Q ratio: {lpq3/lpq2:.2f}")

        # lambda/Q should grow, but sub-linearly relative to Q growth.
        # If lambda ~ Q^alpha with alpha < 2, then lambda/Q ~ Q^(alpha-1)
        # grows slower than Q itself.  We check that lambda/Q doesn't
        # grow faster than Q (which would indicate alpha > 2, i.e. worse
        # than dense scaling).
        assert lpq3 / lpq2 < q_ratio, (
            f"lambda/Q grew faster than Q itself: "
            f"{lpq3/lpq2:.2f}x vs Q ratio {q_ratio:.2f}x"
        )


# ---------------------------------------------------------------------------
# Commutator bound tests
# ---------------------------------------------------------------------------

class TestBinaryConversion:
    """Tests for Pauli string -> binary vector conversion."""

    def test_binary_conversion_identity(self) -> None:
        """Identity term () converts to all-zeros."""
        x, z = _pauli_string_to_binary((), n_qubits=4)
        assert np.all(x == 0)
        assert np.all(z == 0)
        assert x.shape == (4,)

    def test_binary_conversion_xyz(self) -> None:
        """X, Y, Z on specific qubits convert correctly."""
        # X0 -> x[0]=1, z[0]=0
        x, z = _pauli_string_to_binary(((0, 'X'),), n_qubits=3)
        assert x[0] == 1 and z[0] == 0
        assert x[1] == 0 and z[1] == 0

        # Y1 -> x[1]=1, z[1]=1
        x, z = _pauli_string_to_binary(((1, 'Y'),), n_qubits=3)
        assert x[1] == 1 and z[1] == 1

        # Z2 -> x[2]=0, z[2]=1
        x, z = _pauli_string_to_binary(((2, 'Z'),), n_qubits=3)
        assert x[2] == 0 and z[2] == 1

        # X0 Z2 -> x=[1,0,0], z=[0,0,1]
        x, z = _pauli_string_to_binary(((0, 'X'), (2, 'Z')), n_qubits=3)
        assert list(x) == [1, 0, 0]
        assert list(z) == [0, 0, 1]


class TestSymplecticProduct:
    """Tests for symplectic inner product (commutativity check)."""

    def test_symplectic_commuting_pair(self) -> None:
        """Z0 and Z1 act on different qubits -> commute (s=0)."""
        import numpy as np
        X = np.array([[0, 0], [0, 0]], dtype=np.uint8)  # Z0, Z1
        Z = np.array([[1, 0], [0, 1]], dtype=np.uint8)
        S = _symplectic_inner_product_matrix(X, Z)
        assert S[0, 1] == 0
        assert S[1, 0] == 0

    def test_symplectic_anticommuting_pair(self) -> None:
        """X0 and Z0 anticommute (s=1)."""
        import numpy as np
        X = np.array([[1], [0]], dtype=np.uint8)  # X0, Z0
        Z = np.array([[0], [1]], dtype=np.uint8)
        S = _symplectic_inner_product_matrix(X, Z)
        assert S[0, 1] == 1
        assert S[1, 0] == 1


class TestCommutatorBound:
    """Tests for the commutator-based Trotter error bound."""

    def test_comm_bound_two_anticommuting_terms(self) -> None:
        """
        H = a*X0 + b*Z0: X and Z anticommute.
        comm_sum = |a|*|b| (one anticommuting pair).
        comm_bound = t^2 * |a|*|b|.
        onenorm_bound = t^2 * (|a|+|b|)^2 / 2.
        """
        a, b = 0.5, 0.3
        op = QubitOperator(((0, 'X'),), a) + QubitOperator(((0, 'Z'),), b)
        result = pauli_commutator_bound(op, time=1.0, epsilon=1e-3)

        assert result['n_anticommuting_pairs'] == 1
        assert result['n_total_pairs'] == 1
        assert result['anticommuting_fraction'] == pytest.approx(1.0)
        assert result['comm_bound'] == pytest.approx(a * b, rel=1e-10)
        assert result['onenorm_bound'] == pytest.approx((a + b)**2 / 2, rel=1e-10)
        assert result['tightening_ratio'] < 1.0

    def test_comm_bound_commuting_terms(self) -> None:
        """
        H = a*Z0 + b*Z1: Z0 and Z1 commute (different qubits).
        comm_bound = 0 (no anticommuting pairs).
        """
        op = QubitOperator(((0, 'Z'),), 0.5) + QubitOperator(((1, 'Z'),), 0.3)
        result = pauli_commutator_bound(op, time=1.0, epsilon=1e-3)

        assert result['n_anticommuting_pairs'] == 0
        assert result['comm_bound'] == 0.0
        assert result['comm_trotter_steps'] == 0

    def test_comm_bound_tighter_than_onenorm(
        self,
        geovac_he_nmax2: LatticeIndex,
    ) -> None:
        """For GeoVac He nmax=2, comm_bound < onenorm_bound."""
        enc = JordanWignerEncoder(geovac_he_nmax2)
        qop = enc.build_qubit_operator()
        result = pauli_commutator_bound(qop, time=1.0, epsilon=1e-3)

        assert result['comm_bound'] < result['onenorm_bound'], (
            f"Commutator bound {result['comm_bound']:.6f} should be tighter "
            f"than 1-norm bound {result['onenorm_bound']:.6f}"
        )
        assert result['tightening_ratio'] < 1.0
        print(f"\n  He nmax=2: tightening ratio = {result['tightening_ratio']:.4f}, "
              f"anticomm fraction = {result['anticommuting_fraction']:.4f}")

    def test_anticommuting_fraction_decreases(
        self,
        geovac_he_nmax2: LatticeIndex,
        geovac_he_nmax3: LatticeIndex,
    ) -> None:
        """
        Anticommuting fraction should decrease with system size.

        As the system grows, selection rules make more pairs commute
        (they act on disjoint qubit subsets).
        """
        enc2 = JordanWignerEncoder(geovac_he_nmax2)
        enc3 = JordanWignerEncoder(geovac_he_nmax3)
        qop2 = enc2.build_qubit_operator()
        qop3 = enc3.build_qubit_operator()

        r2 = pauli_commutator_bound(qop2, time=1.0, epsilon=1e-3)
        r3 = pauli_commutator_bound(qop3, time=1.0, epsilon=1e-3)

        print(f"\n  nmax=2: anticomm fraction = {r2['anticommuting_fraction']:.4f}")
        print(f"  nmax=3: anticomm fraction = {r3['anticommuting_fraction']:.4f}")

        assert r3['anticommuting_fraction'] < r2['anticommuting_fraction'], (
            f"Anticommuting fraction should decrease: nmax=3 "
            f"({r3['anticommuting_fraction']:.4f}) >= nmax=2 "
            f"({r2['anticommuting_fraction']:.4f})"
        )
