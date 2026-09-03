"""Paper 32, Remark ``rem:D_GV_no_graph_form`` -- module-contents pin.

The 2026-09-02 trunk certification run withdrew Paper 32's ``graph form''
Dirac operator (``def:D_GV_graph`` + ``prop:D_equiv``): the proposition was
tautological and the module it cited holds no hopping operator.  The Remark
that replaced it makes four concrete statements about what the code
actually contains.  This file pins each of them so the Remark cannot drift
back into the withdrawn wording:

  (1) ``geovac.dirac_matrix_elements`` exposes closed-form angular
      (Szmytkowski) and radial (hydrogenic) matrix elements in the
      (kappa, m_j) basis -- and NO hopping / adjacency / edge-set operator.
  (2) The D_GV used for the axiom checks is the diagonal
      ``camporesi_higuchi_full_dirac_matrix``: every eigenvalue is a
      half-integer, each with a positive-integer multiplicity, and the
      spectrum is symmetric (chirality flips the sign).
  (3) The real structure J_GV is ``geovac.real_structure.build_J_full_dirac``
      and satisfies the KO-3 signs on THAT operator: J^2 = -1, JD = +DJ.
  (4) ``geovac.dirac_lattice.DiracLattice`` is a different object: its
      adjacency is an E1-dipole selection-rule graph, not D_GV -- the
      adjacency has off-diagonal support, D_GV has none.

Tier: guard (pins a withdrawn claim; no new physics).
"""

from __future__ import annotations

import inspect
from fractions import Fraction

import numpy as np
import pytest

import geovac.dirac_matrix_elements as dme
from geovac.dirac_lattice import DiracLattice
from geovac.full_dirac_operator_system import (
    camporesi_higuchi_full_dirac_matrix,
    full_dirac_basis,
)
from geovac.real_structure import build_J_full_dirac


# ---------------------------------------------------------------------------
# (1) dirac_matrix_elements: matrix elements, not a hopping operator
# ---------------------------------------------------------------------------


def test_dirac_matrix_elements_exposes_closed_form_matrix_elements():
    for name in ("angular_matrix_sigma", "angular_matrix_L",
                 "angular_matrix_r_hat", "radial_matrix_element"):
        assert callable(getattr(dme, name)), name


def test_dirac_matrix_elements_has_no_hopping_operator():
    """The withdrawn definition said the module documents 'the edge set and
    weights' of a hopping operator.  It does not: no public callable builds
    an adjacency, edge set, hopping matrix, or graph Dirac."""
    forbidden = ("hopping", "adjacency", "edge", "graph_dirac", "build_d")
    public = [n for n, obj in vars(dme).items()
              if not n.startswith("_") and callable(obj)
              and getattr(obj, "__module__", None) == dme.__name__]
    offenders = [n for n in public if any(f in n.lower() for f in forbidden)]
    assert offenders == [], offenders
    # And the module does not import scipy.sparse (no sparse adjacency).
    src = inspect.getsource(dme)
    assert "scipy.sparse" not in src


# ---------------------------------------------------------------------------
# (2) the diagonal D_GV: half-integer spectrum, integer multiplicities
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max", [1, 2, 3, 4])
def test_D_GV_is_diagonal_half_integer_with_integer_multiplicities(n_max):
    basis = full_dirac_basis(n_max)
    D = camporesi_higuchi_full_dirac_matrix(basis)
    # diagonal
    assert np.count_nonzero(D - np.diag(np.diag(D))) == 0
    diag = np.real(np.diag(D))
    assert np.allclose(np.imag(np.diag(D)), 0.0)
    # every eigenvalue is a half-integer
    for lam in diag:
        assert Fraction(float(lam)).limit_denominator(2).denominator == 2, lam
    # multiplicities are positive integers and symmetric under sign flip
    vals, counts = np.unique(diag, return_counts=True)
    assert all(c >= 1 for c in counts)
    mult = dict(zip(vals.tolist(), counts.tolist()))
    for v, c in mult.items():
        assert mult[-v] == c, (v, c, mult.get(-v))
    # |lambda| = n_fock + 1/2 with multiplicity n_fock (n_fock + 1) per
    # chirality (= sum_{l<n} (2l+2)); both chiralities together give the
    # Camporesi-Higuchi degeneracy 2 (n_D+1)(n_D+2) with n_D = n_fock - 1.
    for n in range(1, n_max + 1):
        assert mult[n + 0.5] == n * (n + 1), (n, mult[n + 0.5])


# ---------------------------------------------------------------------------
# (3) J_GV on that operator: KO-3 signs
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max", [1, 2, 3])
def test_J_full_dirac_KO3_signs_on_the_diagonal_D(n_max):
    basis = full_dirac_basis(n_max)
    D = camporesi_higuchi_full_dirac_matrix(basis)
    J = build_J_full_dirac(n_max)
    U = J.U
    assert U.shape == D.shape
    # J = U K, so J^2 = U conj(U) and J D J^{-1} = U conj(D) U^dagger.
    JJ = U @ np.conj(U)
    assert np.allclose(JJ, -np.eye(len(basis)))
    JDJinv = U @ np.conj(D) @ U.conj().T
    assert np.allclose(JDJinv, D)


# ---------------------------------------------------------------------------
# (4) DiracLattice is a different operator
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max", [2, 3])
def test_dirac_lattice_adjacency_is_not_D_GV(n_max):
    lat = DiracLattice(n_max, mode="s3")
    A = lat.adjacency.toarray()
    # E1-dipole adjacency has off-diagonal support and zero diagonal ...
    assert np.count_nonzero(A) > 0
    assert np.count_nonzero(np.diag(A)) == 0
    # ... whereas the CH eigenvalues it also carries are purely diagonal data.
    lam = lat.dirac_eigenvalues
    assert lam.shape == (lat.num_states,)
    assert not np.allclose(A, np.diag(lam))
