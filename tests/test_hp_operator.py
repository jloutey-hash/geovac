"""Tests for Track RH-N Hilbert-Polya operator module.

Verifies:
(1) Each construction produces a Hermitian matrix.
(2) Each construction reproduces the input spectrum (to numerical precision
    for diagonal/tridiagonal/companion; to the Szegő asymptotic for toeplitz).
(3) The diagonal construction gives exactly diag(sorted eigenvalues) up to sort.
(4) analyze_hp_structure runs on all four and reports sensible defaults.
(5) compare_to_dirac runs without error and produces a consistent Pearson r.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

from geovac.hp_operator import (
    build_hp_operator_from_eigenvalues,
    analyze_hp_structure,
    compare_to_dirac,
    construction_summary,
    load_rhm_zeros,
)


# -- Small synthetic spectrum ------------------------------------------

TEST_EIGS = np.array(
    [1.0, 2.5, 3.7, 5.1, 6.9, 8.2, 10.0, 12.4, 15.0, 17.8], dtype=float
)


def test_diagonal_exact_spectrum():
    """Diagonal construction recovers eigenvalues exactly."""
    H = build_hp_operator_from_eigenvalues(TEST_EIGS, construction="diagonal")
    eigs = np.sort(np.linalg.eigvalsh(H))
    assert np.allclose(eigs, np.sort(TEST_EIGS), atol=1e-12)


def test_tridiagonal_exact_spectrum():
    """Tridiagonal Jacobi construction recovers eigenvalues to numerical precision."""
    H = build_hp_operator_from_eigenvalues(
        TEST_EIGS, construction="tridiagonal", seed=0
    )
    assert H.shape == (len(TEST_EIGS), len(TEST_EIGS))
    eigs = np.sort(np.linalg.eigvalsh(H))
    assert np.allclose(eigs, np.sort(TEST_EIGS), atol=1e-8)


def test_companion_exact_spectrum():
    """Unitary-conjugated diagonal construction recovers eigenvalues."""
    H = build_hp_operator_from_eigenvalues(
        TEST_EIGS, construction="companion", seed=1
    )
    eigs = np.sort(np.linalg.eigvalsh(H))
    assert np.allclose(eigs, np.sort(TEST_EIGS), atol=1e-10)


def test_hermitian_all_constructions():
    """All four constructions produce (numerically) symmetric matrices."""
    for name in ("diagonal", "tridiagonal", "toeplitz", "companion"):
        H = build_hp_operator_from_eigenvalues(TEST_EIGS, construction=name)
        asym = np.linalg.norm(H - H.T)
        frob = np.linalg.norm(H)
        assert asym < 1e-10 * max(frob, 1.0), (
            f"construction {name}: asymmetry {asym} vs frob {frob}"
        )


def test_diagonal_sort_matches_eigenvalues_sort():
    """Diagonal construction's diagonal IS the input order (unsorted)."""
    H = build_hp_operator_from_eigenvalues(TEST_EIGS, construction="diagonal")
    assert np.allclose(np.diag(H), TEST_EIGS)
    # and off-diagonals are exactly zero
    assert np.max(np.abs(H - np.diag(np.diag(H)))) == 0.0


def test_toeplitz_is_toeplitz():
    """The toeplitz construction actually gives a Toeplitz matrix."""
    H = build_hp_operator_from_eigenvalues(TEST_EIGS, construction="toeplitz")
    n = H.shape[0]
    # Check H[j,k] depends only on |j-k|
    for d in range(n):
        vals = []
        for j in range(n):
            for k in range(n):
                if abs(j - k) == d:
                    vals.append(H[j, k])
        vals = np.array(vals)
        assert np.max(np.abs(vals - vals[0])) < 1e-10, (
            f"Toeplitz matrix not Toeplitz at distance {d}"
        )


def test_analyze_hp_structure_runs():
    """Structure analysis produces a dict with the expected keys."""
    H = build_hp_operator_from_eigenvalues(TEST_EIGS, construction="tridiagonal")
    rep = analyze_hp_structure(H)
    for key in (
        "trace",
        "frobenius",
        "density",
        "offdiag_density",
        "min_eig",
        "max_eig",
        "is_symmetric",
    ):
        assert key in rep, f"missing key {key}"
    assert rep["is_symmetric"]


def test_analyze_with_labels_block_structure():
    """Providing labels gives block statistics; synthetic block-diagonal H
    should have block_diagonal_fraction == 1.0."""
    # Block-diag with two 2x2 blocks
    H = np.zeros((4, 4), dtype=float)
    H[:2, :2] = np.array([[1.0, 0.3], [0.3, 2.0]])
    H[2:, 2:] = np.array([[3.0, -0.1], [-0.1, 4.0]])
    labels = [(0,), (0,), (1,), (1,)]
    rep = analyze_hp_structure(H, labels=labels)
    assert rep["block_stats"]["block_diagonal_fraction"] == pytest.approx(1.0, abs=1e-12)
    # Group sizes
    assert rep["block_stats"]["group_sizes"][0] == 2
    assert rep["block_stats"]["group_sizes"][1] == 2


def test_compare_to_dirac_runs():
    """Dirac comparison returns numerical fields on a small H."""
    H = build_hp_operator_from_eigenvalues(TEST_EIGS, construction="diagonal")
    rep = compare_to_dirac(H, n_max=5, rescale=True)
    for key in (
        "dirac_eigs_used_first_n",
        "H_eigs_sorted",
        "H_eigs_rescaled_to_dirac_range",
        "residual_frobenius",
        "pearson_r_sorted_eigs",
    ):
        assert key in rep
    # After rescaling the min and max should match
    H_rescaled = rep["H_eigs_rescaled_to_dirac_range"]
    dirac = rep["dirac_eigs_used_first_n"]
    assert abs(H_rescaled[0] - dirac[0]) < 1e-10
    assert abs(H_rescaled[-1] - dirac[-1]) < 1e-10
    # Pearson r for a strictly increasing sequence should be positive
    assert rep["pearson_r_sorted_eigs"] > 0.5


def test_load_rhm_zeros_roundtrip(tmp_path):
    """load_rhm_zeros reads a fixture dict and returns sorted floats."""
    import json

    data = {
        "zeros": {
            "D_full": [[2.5, 30.0], [3.0, 20.0], [2.9, 40.0]],
        }
    }
    path = tmp_path / "zeros.json"
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f)
    ims = load_rhm_zeros(str(path), which="D_full")
    assert list(ims) == [20.0, 30.0, 40.0]


def test_construction_summary_full_roundtrip():
    """construction_summary produces a dict with all four constructions."""
    summary = construction_summary(TEST_EIGS, seed=0)
    for name in ("diagonal", "tridiagonal", "toeplitz", "companion"):
        assert name in summary
        entry = summary[name]
        # each has structure + dirac_comparison + spectrum_max_abs_err
        assert "spectrum_max_abs_err" in entry
        assert "structure" in entry
        assert "dirac_comparison" in entry
    # exact constructions reproduce the spectrum tightly
    assert summary["diagonal"]["spectrum_max_abs_err"] < 1e-12
    assert summary["tridiagonal"]["spectrum_max_abs_err"] < 1e-8
    assert summary["companion"]["spectrum_max_abs_err"] < 1e-10


def test_input_validation_bad_construction():
    """Unknown construction name raises ValueError."""
    with pytest.raises(ValueError):
        build_hp_operator_from_eigenvalues(TEST_EIGS, construction="blah")


def test_input_validation_bad_shape():
    """Non-1-D eigenvalue input raises."""
    with pytest.raises(ValueError):
        build_hp_operator_from_eigenvalues(
            np.array([[1.0, 2.0]]), construction="diagonal"
        )
