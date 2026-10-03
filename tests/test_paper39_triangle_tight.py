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

Hence the triangle route cannot do better than sqrt(2), and the Pythagorean
'< 1' is false.  The constant itself is C_3^(2) = 1 by compression on the
Clifford-doubled module (Paper 39 Lemma L3-T, 2026-10-02;
tests/test_paper39_compression_doubled.py).  See Paper 39 Remark "Withdrawn
Pythagorean refinement".

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


def test_c3_record_formula_withdrawn() -> None:
    """The per-irrep closed form (sqrt(2n/(n+2)) on the diagonal) rests on Paper 38's
    FALSE per-harmonic bound ||[D, M_NLM]|| <= N-1 and is withdrawn as a value of
    C_3^(2) (Paper 39, 2026-10-02); the constant is 1 by compression.  Rejects the
    record function being presented as a live Lipschitz constant: its docstring
    must open with the WITHDRAWN marker, and the module's theorem constant must not
    be built from it."""
    from geovac.gh_convergence_tensor import (
        c3_full_triangle_bound as c3, compute_tensor_propinquity_bound)
    doc = (c3.__doc__ or "").lstrip()
    assert doc.startswith("[WITHDRAWN"), doc[:80]
    b = compute_tensor_propinquity_bound(3, 4, gamma_prec=15)
    gmax = max(b.gamma_a, b.gamma_b)
    assert abs(b.propinquity_bound_theorem - 2.0 * gmax) < 1e-9 * gmax
    assert abs(b.propinquity_bound_theorem - 2.0 * c3(3, 4) * gmax) > 1e-3
