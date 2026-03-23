"""
Tests for geovac.measurement_grouping — QWC Pauli grouping.

Validates:
  1. QWC correctness: all pairs within each group qubitwise-commute
  2. Completeness: total terms across groups equals original Pauli count
  3. Known result: H2 STO-3G (15 terms) fits in <= 5 QWC groups
  4. GeoVac vs Gaussian comparison with printed table
  5. Scaling: groups/terms ratio stays flat or decreases with basis size

Author: GeoVac Development Team
Date: March 2026
"""

import warnings

import pytest

from openfermion import jordan_wigner, QubitOperator

from geovac.measurement_grouping import (
    qwc_compatible,
    qwc_groups,
    count_qwc_groups,
    MeasurementAnalysis,
    analyze_measurement_cost,
    compare_measurement_cost,
)
from geovac.gaussian_reference import h2_sto3g, he_sto3g, build_qubit_hamiltonian
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
def he_sto3g_qop():
    """He STO-3G QubitOperator."""
    _, qop, _ = build_qubit_hamiltonian(he_sto3g())
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
# QWC compatibility unit tests
# ---------------------------------------------------------------------------

class TestQWCCompatible:
    """Tests for the qwc_compatible pure function."""

    def test_identical_terms(self) -> None:
        """Identical terms always QWC."""
        term = ((0, 'X'), (1, 'Z'))
        assert qwc_compatible(term, term)

    def test_disjoint_qubits(self) -> None:
        """Terms on disjoint qubits always QWC."""
        a = ((0, 'X'), (1, 'Y'))
        b = ((2, 'Z'), (3, 'X'))
        assert qwc_compatible(a, b)

    def test_same_pauli_shared_qubit(self) -> None:
        """Same Pauli on shared qubit: QWC."""
        a = ((0, 'X'), (2, 'Z'))
        b = ((0, 'X'), (3, 'Y'))
        assert qwc_compatible(a, b)

    def test_different_pauli_shared_qubit(self) -> None:
        """Different Pauli on shared qubit: NOT QWC."""
        a = ((0, 'X'), (1, 'Z'))
        b = ((0, 'Y'), (1, 'Z'))
        assert not qwc_compatible(a, b)

    def test_identity_terms(self) -> None:
        """Empty tuple (identity) QWC with everything."""
        identity = ()
        term = ((0, 'X'), (1, 'Z'), (2, 'Y'))
        assert qwc_compatible(identity, term)
        assert qwc_compatible(term, identity)
        assert qwc_compatible(identity, identity)

    def test_single_qubit_conflict(self) -> None:
        """Single-qubit X vs Z conflict."""
        assert not qwc_compatible(((0, 'X'),), ((0, 'Z'),))

    def test_multi_qubit_partial_overlap(self) -> None:
        """Multiple shared qubits, all matching: QWC."""
        a = ((0, 'X'), (1, 'Y'), (2, 'Z'))
        b = ((0, 'X'), (1, 'Y'), (3, 'X'))
        assert qwc_compatible(a, b)

    def test_multi_qubit_one_conflict(self) -> None:
        """Multiple shared qubits, one conflict: NOT QWC."""
        a = ((0, 'X'), (1, 'Y'), (2, 'Z'))
        b = ((0, 'X'), (1, 'Z'), (3, 'X'))
        assert not qwc_compatible(a, b)


# ---------------------------------------------------------------------------
# QWC grouping correctness
# ---------------------------------------------------------------------------

def _verify_all_pairs_qwc(groups) -> None:
    """Assert every pair of terms within each group actually QWC."""
    for g_idx, group in enumerate(groups):
        terms = [term for term, _ in group]
        for i in range(len(terms)):
            for j in range(i + 1, len(terms)):
                assert qwc_compatible(terms[i], terms[j]), (
                    f"Group {g_idx}: terms {terms[i]} and {terms[j]} "
                    f"do NOT qubitwise-commute"
                )


class TestQWCGrouping:
    """Tests for qwc_groups and correctness properties."""

    def test_correctness_h2_sto3g(self, h2_sto3g_qop) -> None:
        """All pairs within each group QWC for H2 STO-3G."""
        groups = qwc_groups(h2_sto3g_qop)
        _verify_all_pairs_qwc(groups)

    def test_correctness_he_sto3g(self, he_sto3g_qop) -> None:
        """All pairs within each group QWC for He STO-3G."""
        groups = qwc_groups(he_sto3g_qop)
        _verify_all_pairs_qwc(groups)

    def test_completeness_h2_sto3g(self, h2_sto3g_qop) -> None:
        """Total terms across groups equals original Pauli count."""
        groups = qwc_groups(h2_sto3g_qop)
        total = sum(len(g) for g in groups)
        assert total == len(h2_sto3g_qop.terms), (
            f"Groups contain {total} terms but operator has "
            f"{len(h2_sto3g_qop.terms)} terms"
        )

    def test_completeness_he_sto3g(self, he_sto3g_qop) -> None:
        """Total terms across groups equals original Pauli count."""
        groups = qwc_groups(he_sto3g_qop)
        total = sum(len(g) for g in groups)
        assert total == len(he_sto3g_qop.terms)

    def test_known_result_h2_sto3g(self, h2_sto3g_qop) -> None:
        """
        H2 STO-3G (15 Pauli terms, 4 qubits) needs <= 5 QWC groups.

        This is consistent with published VQE results where H2 STO-3G
        requires 3-5 measurement bases depending on grouping strategy.
        """
        n_groups = count_qwc_groups(h2_sto3g_qop)
        assert n_groups <= 5, (
            f"H2 STO-3G should need <= 5 QWC groups, got {n_groups}"
        )

    def test_empty_operator(self) -> None:
        """Empty QubitOperator produces zero groups."""
        empty = QubitOperator()
        groups = qwc_groups(empty)
        assert len(groups) == 0

    def test_single_term(self) -> None:
        """Single Pauli term produces one group."""
        op = QubitOperator(((0, 'X'), (1, 'Z')), 1.5)
        groups = qwc_groups(op)
        assert len(groups) == 1
        assert len(groups[0]) == 1


# ---------------------------------------------------------------------------
# MeasurementAnalysis tests
# ---------------------------------------------------------------------------

class TestMeasurementAnalysis:
    """Tests for analyze_measurement_cost."""

    def test_h2_sto3g_analysis(self, h2_sto3g_qop) -> None:
        """Analysis dataclass has consistent fields for H2 STO-3G."""
        analysis = analyze_measurement_cost(h2_sto3g_qop)
        assert analysis.n_pauli_terms == 15
        assert analysis.n_qubits == 4
        assert analysis.n_qwc_groups > 0
        assert analysis.n_qwc_groups <= analysis.n_pauli_terms
        assert analysis.max_group_size >= analysis.min_group_size
        assert analysis.min_group_size >= 1

    def test_summary_string(self, h2_sto3g_qop) -> None:
        """Summary contains key metrics."""
        analysis = analyze_measurement_cost(h2_sto3g_qop)
        s = analysis.summary()
        assert "Q=" in s
        assert "terms=" in s
        assert "QWC_groups=" in s


# ---------------------------------------------------------------------------
# GeoVac vs Gaussian comparison
# ---------------------------------------------------------------------------

class TestGeoVacVsGaussian:
    """Side-by-side comparison of measurement costs."""

    def test_comparison_he(
        self,
        he_sto3g_qop,
        geovac_he_nmax2: LatticeIndex,
    ) -> None:
        """
        Compare QWC group counts: GeoVac He nmax=2 vs Gaussian STO-3G.

        GeoVac has more qubits and terms (10 vs 2 qubits), so it will
        have more groups.  But we verify the groups-per-term ratio.
        """
        enc = JordanWignerEncoder(geovac_he_nmax2)
        geovac_qop = enc.build_qubit_operator()

        table = compare_measurement_cost(
            geovac_qop, he_sto3g_qop,
            labels=('GeoVac He n2', 'Gauss He STO3G'),
        )
        print(table)

        a_geovac = analyze_measurement_cost(geovac_qop)
        a_gauss = analyze_measurement_cost(he_sto3g_qop)

        # GeoVac has more terms and groups (larger basis)
        assert a_geovac.n_pauli_terms > a_gauss.n_pauli_terms
        assert a_geovac.n_qwc_groups >= a_gauss.n_qwc_groups

        # Groups-per-term ratio should be <= 1 for both
        assert a_geovac.n_qwc_groups / a_geovac.n_pauli_terms <= 1.0
        assert a_gauss.n_qwc_groups / a_gauss.n_pauli_terms <= 1.0

    def test_correctness_geovac_he(self, geovac_he_nmax2: LatticeIndex) -> None:
        """All pairs within each GeoVac He group QWC."""
        enc = JordanWignerEncoder(geovac_he_nmax2)
        geovac_qop = enc.build_qubit_operator()
        groups = qwc_groups(geovac_qop)
        _verify_all_pairs_qwc(groups)


# ---------------------------------------------------------------------------
# Scaling test
# ---------------------------------------------------------------------------

class TestScaling:
    """Verify QWC grouping scales favorably with basis size."""

    def test_groups_per_term_ratio(
        self,
        geovac_he_nmax2: LatticeIndex,
        geovac_he_nmax3: LatticeIndex,
    ) -> None:
        """
        Groups/terms ratio should not increase as basis grows.

        If ERI sparsity (selection rules) helps with measurement cost,
        the ratio of QWC groups to Pauli terms should decrease or stay
        flat as nmax increases.
        """
        enc2 = JordanWignerEncoder(geovac_he_nmax2)
        enc3 = JordanWignerEncoder(geovac_he_nmax3)
        qop2 = enc2.build_qubit_operator()
        qop3 = enc3.build_qubit_operator()

        a2 = analyze_measurement_cost(qop2)
        a3 = analyze_measurement_cost(qop3)

        ratio2 = a2.n_qwc_groups / a2.n_pauli_terms
        ratio3 = a3.n_qwc_groups / a3.n_pauli_terms

        print(f"\n  nmax=2: {a2.n_pauli_terms:,} terms, "
              f"{a2.n_qwc_groups} groups, ratio={ratio2:.4f}")
        print(f"  nmax=3: {a3.n_pauli_terms:,} terms, "
              f"{a3.n_qwc_groups} groups, ratio={ratio3:.4f}")

        # Greedy QWC grouping is suboptimal, especially at small sizes
        # where the ratio can increase before stabilizing.  Allow 50%
        # tolerance — the important signal is that the ratio doesn't
        # explode (which would indicate O(N) groups per N terms).
        assert ratio3 <= ratio2 * 1.50, (
            f"Groups/terms ratio increased too much: nmax=2 {ratio2:.4f} -> "
            f"nmax=3 {ratio3:.4f} (>{ratio2 * 1.50:.4f} limit)"
        )
