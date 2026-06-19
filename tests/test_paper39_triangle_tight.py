"""Paper 39 backing test: the joint Leibniz bound is the TIGHT TRIANGLE bound,
not the (false) Pythagorean identity.

Corrects the withdrawn keystone (2026-06-18, /qa group1 Bite B-3). Paper 39's
Lemma L3-T once claimed the operator-norm Pythagorean identity
    ||[D_ab, M_f (x) M_g]||^2 = ||[D_a,M_f]||^2 ||M_g||^2 + ||M_f||^2 ||[D_b,M_g]||^2
on the Connes-Marcolli graded module, giving C_3^(2) < 1. That identity is FALSE
(graded anticommutation does not imply operator-norm orthogonality). This test
builds the joint commutator from the REAL Camporesi-Higuchi Dirac matrices and the
graded multipliers and confirms, over every single-harmonic tensor pair:

  (1) the Pythagorean ratio ||A+B||^2 / (||A||^2 + ||B||^2) reaches its MAXIMUM 2
      (the would-be identity is maximally violated -- the two Leibniz terms add
      constructively at a pair of top-shell raising harmonics, ||A||=||B||);
  (2) the triangle ratio ||A+B|| / (||A|| + ||B||) reaches its MAXIMUM 1
      (the triangle bound is TIGHT / achieved).

Hence the correct constant is C_3^(2) <= sqrt(2) (-> sqrt(2)), not < 1. The
convergence theorem is unaffected (sqrt(2) is a finite constant). See Paper 39
Remark "Withdrawn Pythagorean refinement" and geovac.gh_convergence_tensor
.c3_full_triangle_bound.

A = [D_a, M_f] (x) M_g ;  B = (gamma_a M_f) (x) [D_b, M_g] ;  A+B = [D_ab, M_f(x)M_g].
"""

from __future__ import annotations

import numpy as np
import pytest

from geovac.full_dirac_operator_system import (
    full_dirac_basis,
    camporesi_higuchi_full_dirac_matrix,
    FullDiracTruncatedOperatorSystem,
)


def _opn(X: np.ndarray) -> float:
    return float(np.linalg.norm(X, 2)) if X.size else 0.0


def _build(n_max: int):
    basis = full_dirac_basis(n_max)
    D = camporesi_higuchi_full_dirac_matrix(basis)
    gamma = np.diag([float(b.chirality) for b in basis]).astype(complex)
    Ms = [np.asarray(M, dtype=complex)
          for M in FullDiracTruncatedOperatorSystem(n_max).multiplier_matrices]
    return D, gamma, Ms


@pytest.mark.slow
@pytest.mark.parametrize("n_max_a,n_max_b", [(2, 2), (2, 3)])
def test_joint_leibniz_triangle_tight_pythagorean_false(n_max_a, n_max_b) -> None:
    Da, ga, Msa = _build(n_max_a)
    Db, _gb, Msb = _build(n_max_b)
    Ib = np.eye(Db.shape[0], dtype=complex)
    Dab = np.kron(Da, Ib) + np.kron(ga, Db)

    max_pyth_ratio = 0.0  # ||A+B||^2 / (||A||^2 + ||B||^2)
    max_tri_ratio = 0.0   # ||A+B||   / (||A|| + ||B||)
    n_pairs = 0
    for Mf in Msa:
        Ca = Da @ Mf - Mf @ Da
        if _opn(Ca) < 1e-9 or _opn(Mf) < 1e-9:
            continue
        for Mg in Msb:
            Cb = Db @ Mg - Mg @ Db
            if _opn(Cb) < 1e-9 or _opn(Mg) < 1e-9:
                continue
            Mfg = np.kron(Mf, Mg)
            lhs = _opn(Dab @ Mfg - Mfg @ Dab)        # ||A+B||
            nA = _opn(Ca) * _opn(Mg)
            nB = _opn(Mf) * _opn(Cb)
            max_pyth_ratio = max(max_pyth_ratio, lhs**2 / (nA**2 + nB**2))
            max_tri_ratio = max(max_tri_ratio, lhs / (nA + nB))
            n_pairs += 1

    assert n_pairs > 0
    # (2) triangle bound is an upper bound and is TIGHT (achieved = 1)
    assert max_tri_ratio <= 1.0 + 1e-9, (
        f"triangle bound must hold; got max ratio {max_tri_ratio}"
    )
    assert max_tri_ratio == pytest.approx(1.0, abs=1e-6), (
        f"triangle bound must be TIGHT (=1); got {max_tri_ratio}"
    )
    # (1) Pythagorean identity is FALSE -- maximally violated (ratio reaches 2)
    assert max_pyth_ratio == pytest.approx(2.0, abs=1e-6), (
        f"Pythagorean ratio should reach 2 (identity maximally violated); "
        f"got {max_pyth_ratio}"
    )


def test_c3_constant_is_triangle_not_pythagorean() -> None:
    """The closed-form C_3 bound -> sqrt(2) (triangle), exceeds 1 beyond (3,3)."""
    from geovac.gh_convergence_tensor import c3_full_triangle_bound as c3
    # c3(n_max,n_max) sups over N in [2, n_max+1]; corner (n_max+1,n_max+1):
    #   sqrt( (2 n_max)^2 / (2 (n_max+1)^2 - 2) ) = sqrt( 2 n_max / (n_max+2) ).
    assert c3(1, 1) == pytest.approx(np.sqrt(2 / 3), abs=1e-9)   # corner (2,2): sqrt(4/6)
    assert c3(2, 2) == pytest.approx(1.0, abs=1e-9)              # corner (3,3): sqrt(16/16)
    assert c3(3, 3) == pytest.approx(np.sqrt(6 / 5), abs=1e-9)   # corner (4,4): sqrt(36/30)
    assert c3(4, 4) > 1.0                                         # exceeds 1
    assert c3(50, 50) < np.sqrt(2)                               # approaches but < sqrt(2)
    assert c3(200, 200) > 1.39                                   # -> sqrt(2)
