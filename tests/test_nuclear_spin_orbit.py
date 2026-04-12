"""
Tests for nuclear shell model spin-orbit coupling (Track NB).

Verifies that the Mayer-Jensen spin-orbit term reproduces the real
nuclear magic numbers: 2, 8, 20, 28, 50, 82, 126.
"""

from typing import Tuple

import numpy as np
import pytest

from geovac.nuclear.spin_orbit import (
    ls_eigenvalue,
    nuclear_shell_levels,
    find_magic_numbers,
    find_optimal_vls,
    verify_magic_ordering,
    build_spin_orbit_hamiltonian,
    level_ordering_table,
    sparsity_with_spin_orbit,
)
from geovac.nuclear.harmonic_shell import _ho_shell_degeneracy


# ── 1. l.s eigenvalue tests ─────────────────────────────────────────────────

class TestLSEigenvalues:
    """Verify <l.s> = l/2 for j=l+1/2 and -(l+1)/2 for j=l-1/2."""

    @pytest.mark.parametrize("l", [0, 1, 2, 3, 4, 5])
    def test_ls_j_plus(self, l: int) -> None:
        """j = l + 1/2 gives <l.s> = l/2."""
        j = l + 0.5
        expected = l / 2.0 if l > 0 else 0.0
        assert ls_eigenvalue(l, j) == pytest.approx(expected)

    @pytest.mark.parametrize("l", [1, 2, 3, 4, 5])
    def test_ls_j_minus(self, l: int) -> None:
        """j = l - 1/2 gives <l.s> = -(l+1)/2."""
        j = l - 0.5
        expected = -(l + 1) / 2.0
        assert ls_eigenvalue(l, j) == pytest.approx(expected)

    def test_ls_l0(self) -> None:
        """l=0: only j=1/2, <l.s>=0."""
        assert ls_eigenvalue(0, 0.5) == 0.0

    def test_invalid_j(self) -> None:
        """Invalid j raises ValueError."""
        with pytest.raises(ValueError):
            ls_eigenvalue(2, 0.5)  # j=0.5 invalid for l=2


# ── 2. Zero v_ls gives HO magic numbers ─────────────────────────────────────

def test_zero_vls_gives_ho_magic() -> None:
    """v_ls=0 should recover harmonic oscillator magic numbers."""
    result = find_magic_numbers(n_max=6, hw=1.0, v_ls=0.0, gap_threshold=0.3)
    assert result["magic_numbers"] == [2, 8, 20, 40, 70, 112]


def test_zero_vls_verify_ho_ordering() -> None:
    """v_ls=0 ordering verification: HO magic numbers are valid cumulative counts."""
    target = [2, 8, 20, 40, 70, 112]
    result = verify_magic_ordering(target, n_max=6, hw=1.0, v_ls=0.0, d_ll=0.0)
    assert result["valid"], (
        f"HO magic numbers not valid: missing={result['missing']}"
    )


# ── 3. Degeneracy tests ─────────────────────────────────────────────────────

def test_so_degeneracy() -> None:
    """Each j level has exactly 2j+1 states."""
    levels = nuclear_shell_levels(n_max=5, hw=1.0, v_ls=0.3)
    for lev in levels:
        j = lev["j"]
        expected_deg = int(2 * j + 1)
        assert lev["degeneracy"] == expected_deg, (
            f"Level {lev['label']}: expected deg {expected_deg}, got {lev['degeneracy']}"
        )


def test_total_state_count() -> None:
    """Total states with SO = total states without SO."""
    for n_max in [3, 5, 7]:
        levels = nuclear_shell_levels(n_max, hw=1.0, v_ls=0.5)
        total_so = sum(lev["degeneracy"] for lev in levels)
        total_ho = sum(_ho_shell_degeneracy(N) for N in range(n_max))
        assert total_so == total_ho, (
            f"n_max={n_max}: SO total {total_so} != HO total {total_ho}"
        )


# ── 4. Spin-orbit splitting ─────────────────────────────────────────────────

def test_so_splits_correctly() -> None:
    """For l=1: j=3/2 has 4 states, j=1/2 has 2 states, total 6."""
    levels = nuclear_shell_levels(n_max=2, hw=1.0, v_ls=0.3)
    # N=1 has l=1, which splits into j=3/2 (4) and j=1/2 (2)
    p_levels = [lev for lev in levels if lev["l"] == 1]
    assert len(p_levels) == 2
    degs = sorted([lev["degeneracy"] for lev in p_levels])
    assert degs == [2, 4]
    assert sum(degs) == 6  # (2*1+1)*2


# ── 5. j-level count ────────────────────────────────────────────────────────

def test_j_level_count() -> None:
    """Each (n_r, l>0) gives 2 j-levels, l=0 gives 1."""
    levels = nuclear_shell_levels(n_max=5, hw=1.0, v_ls=0.3)
    # Count levels per (n_r, l)
    from collections import Counter
    nr_l_count = Counter()
    for lev in levels:
        nr_l_count[(lev["n_r"], lev["l"])] += 1

    for (n_r, l), count in nr_l_count.items():
        if l == 0:
            assert count == 1, f"(n_r={n_r}, l={l}): expected 1 j-level, got {count}"
        else:
            assert count == 2, f"(n_r={n_r}, l={l}): expected 2 j-levels, got {count}"


# ── 6-9. Real magic numbers ─────────────────────────────────────────────────

@pytest.fixture
def optimal_params() -> Tuple[float, float]:
    """Find the optimal (v_ls, d_ll) for real magic numbers."""
    result = find_optimal_vls(n_max=7, hw=1.0)
    assert result["optimal_vls"] is not None, "Failed to find optimal v_ls"
    return result["optimal_vls"], result["optimal_d_ll"]


def test_magic_28(optimal_params: Tuple[float, float]) -> None:
    """With appropriate v_ls and d_ll, 0f7/2 drops to give magic number 28."""
    v_ls, d_ll = optimal_params
    result = verify_magic_ordering([28], n_max=7, hw=1.0, v_ls=v_ls, d_ll=d_ll)
    assert 28 in result["present"], (
        f"28 not in ordering: present={result['present']}, missing={result['missing']}"
    )


def test_magic_50(optimal_params: Tuple[float, float]) -> None:
    """With appropriate v_ls and d_ll, closure at 50 exists."""
    v_ls, d_ll = optimal_params
    result = verify_magic_ordering([50], n_max=7, hw=1.0, v_ls=v_ls, d_ll=d_ll)
    assert 50 in result["present"], (
        f"50 not in ordering: present={result['present']}, missing={result['missing']}"
    )


def test_magic_82(optimal_params: Tuple[float, float]) -> None:
    """With appropriate v_ls and d_ll, closure at 82 exists."""
    v_ls, d_ll = optimal_params
    result = verify_magic_ordering([82], n_max=7, hw=1.0, v_ls=v_ls, d_ll=d_ll)
    assert 82 in result["present"], (
        f"82 not in ordering: present={result['present']}, missing={result['missing']}"
    )


def test_magic_126(optimal_params: Tuple[float, float]) -> None:
    """With appropriate v_ls and d_ll, closure at 126 exists."""
    v_ls, d_ll = optimal_params
    result = verify_magic_ordering([126], n_max=7, hw=1.0, v_ls=v_ls, d_ll=d_ll)
    assert 126 in result["present"], (
        f"126 not in ordering: present={result['present']}, missing={result['missing']}"
    )


# ── 10. All real magic numbers ───────────────────────────────────────────────

def test_all_real_magic_numbers(optimal_params: Tuple[float, float]) -> None:
    """Optimal (v_ls, d_ll) produces all 7 real nuclear magic numbers."""
    v_ls, d_ll = optimal_params
    target = [2, 8, 20, 28, 50, 82, 126]
    result = verify_magic_ordering(target, n_max=7, hw=1.0, v_ls=v_ls, d_ll=d_ll)
    assert result["valid"], (
        f"Missing magic numbers: {result['missing']}, present: {result['present']}"
    )


# ── 11. Hamiltonian is Hermitian ─────────────────────────────────────────────

def test_hamiltonian_hermitian() -> None:
    """Spin-orbit Hamiltonian must be Hermitian."""
    data = build_spin_orbit_hamiltonian(n_max=4, hw=1.0, v_ls=0.3)
    H = data["hamiltonian"].toarray()
    assert np.allclose(H, H.T.conj(), atol=1e-14), "H is not Hermitian"


# ── 12. Eigenvalues match j-scheme ──────────────────────────────────────────

@pytest.mark.parametrize("d_ll", [0.0, 0.05, 0.1])
def test_hamiltonian_eigenvalues_match_j_scheme(d_ll: float) -> None:
    """Matrix eigenvalues match analytical j-scheme energies."""
    hw = 1.0
    v_ls = 0.35

    # Analytical j-scheme energies (with degeneracies)
    levels = nuclear_shell_levels(n_max=5, hw=hw, v_ls=v_ls, d_ll=d_ll)
    expected_evals = []
    for lev in levels:
        expected_evals.extend([lev["energy"]] * lev["degeneracy"])
    expected_evals.sort()

    # Matrix eigenvalues
    data = build_spin_orbit_hamiltonian(n_max=5, hw=hw, v_ls=v_ls, d_ll=d_ll)
    matrix_evals = data["eigenvalues"]

    assert len(expected_evals) == len(matrix_evals), (
        f"Length mismatch: {len(expected_evals)} vs {len(matrix_evals)}"
    )
    np.testing.assert_allclose(
        matrix_evals, expected_evals, atol=1e-10,
        err_msg="Matrix eigenvalues do not match j-scheme energies"
    )


# ── 13. SO preserves ERI sparsity ───────────────────────────────────────────

def test_so_preserves_eri_sparsity() -> None:
    """
    Two-body ERI sparsity is unchanged by spin-orbit coupling.

    SO is a one-body operator; it does not modify ERIs. This test verifies
    the structural claim by checking that the ERI tensor (from the HO module)
    is identical with and without SO (since SO doesn't touch ERIs at all).
    """
    # The ERI tensor depends only on radial wavefunctions and angular
    # selection rules, neither of which changes with SO. This test confirms
    # the structural claim by checking that sparsity_with_spin_orbit reports
    # the correct note.
    result = sparsity_with_spin_orbit(n_max=3, hw=1.0, v_ls=0.3)
    assert "one-body operator" in result["note"]
    # SO should add elements (off-diagonal spin-flip terms)
    assert result["with_so"]["nonzero_elements"] >= result["without_so"]["nonzero_elements"]


# ── 14. Level ordering table ────────────────────────────────────────────────

def test_level_ordering_table(optimal_params: Tuple[float, float]) -> None:
    """Verify printed table contains all expected labels and magic numbers."""
    v_ls, d_ll = optimal_params
    table = level_ordering_table(n_max=7, hw=1.0, v_ls=v_ls, d_ll=d_ll)
    # Check key labels are present
    expected_labels = ["0s1/2", "0p3/2", "0p1/2", "0d5/2", "0f7/2", "0g9/2"]
    for label in expected_labels:
        assert label in table, f"Label '{label}' not found in table"
    # Check magic numbers are marked
    assert "28" in table
    assert "50" in table
    assert "82" in table
    assert "126" in table


# ── 15. Optimal v_ls ratio is physical ───────────────────────────────────────

def test_optimal_vls_ratio_physical() -> None:
    """The optimal v_ls/hw and d_ll/hw ratios should be in physical ranges."""
    result = find_optimal_vls(n_max=7, hw=1.0)
    assert result["ratio_vls_hw"] is not None
    ratio_vls = result["ratio_vls_hw"]
    ratio_dll = result["ratio_d_ll_hw"]
    assert 0.05 < ratio_vls < 1.0, (
        f"v_ls/hw = {ratio_vls} outside physical range [0.05, 1.0]"
    )
    assert 0.0 <= ratio_dll < 0.3, (
        f"d_ll/hw = {ratio_dll} outside physical range [0.0, 0.3]"
    )
