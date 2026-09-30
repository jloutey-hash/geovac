"""Paper 32 axiom audit, item 5 + "Reading of the order-condition failures" +
"Net status of the real structure": the order-zero residual GROWS with the
cutoff; it does not decay.

Paper 32 (Sec. axiom audit, item 5; [MEASURED], 2026-09-28 re-measurement,
Weyl sector) prints:

  * failing order-zero pairs [a, J b J^{-1}] != 0:
        124/196 (n_max=2), 2594/3025 (n_max=3), 17970/19600 (n_max=4)
    i.e. failing fraction 63% -> 86% -> 92%;
  * max-abs residual 5.5e-2 -> 7.9e-2 -> 9.6e-2 at n_max = 2, 3, 4;
  * on a FIXED set of generators (the degree-1 harmonics, N = 2) the
    operator-norm residual runs 0.051, 0.068, 0.076, 0.081 at n_max = 2..5,
    with ratio to the growing generator norm 0.225 -> 0.208, still
    drifting; whether it tends to a nonzero limit is NOT established.

WRONG ANSWER REJECTED: "the residuals decay as n_max -> infinity" (the
finite-resolution-artifact reading an earlier version of the paper inferred
from GH convergence, now refuted by the measurement).  A residual sequence
that is flat or decreasing, or a fixed-generator residual heading to 0,
fails these tests.

Fire-test plants (debug/qa/fire_test.py, --plant-in geovac/real_structure.py):
  * verify_order_zero residual divided by the Hilbert-space dimension
    (a decaying residual)  -> FIRES;
  * apply_to_operator result divided by dim (J-conjugate shrinks with the
    cutoff, so every commutator decays) -> FIRES.
"""
from __future__ import annotations

import numpy as np
import pytest

from geovac.real_structure import build_J_weyl, verify_order_zero
from geovac.spinor_operator_system import SpinorTruncatedOperatorSystem


# Paper 32 item 5 printed values (Weyl sector)
PRINTED = {
    2: (124, 196, 5.5e-2),
    3: (2594, 3025, 7.9e-2),
    4: (17970, 19600, 9.6e-2),
}
FIXED_GEN_PRINTED = {2: 0.051, 3: 0.068, 4: 0.076, 5: 0.081}


def _order_zero(n_max: int):
    op = SpinorTruncatedOperatorSystem(n_max)
    J = build_J_weyl(n_max)
    ok, failures = verify_order_zero(J, op.multiplier_matrices, tol=1e-9)
    n_total = len(op.multiplier_matrices) ** 2
    max_res = max((r for _, _, r in failures), default=0.0)
    return len(failures), n_total, max_res


def _fixed_generator_residual(n_max: int):
    """Max operator-norm ||[a, J b J^{-1}]|| over the degree-1 (N = 2)
    multipliers a, b, and the max operator norm of those generators."""
    op = SpinorTruncatedOperatorSystem(n_max)
    J = build_J_weyl(n_max)
    idx = [i for i, lab in enumerate(op.multiplier_labels) if lab[0] == 2]
    assert len(idx) == 4, f"expected the 4 degree-1 harmonics, got {len(idx)}"
    worst = 0.0
    for i in idx:
        a = op.multiplier_matrices[i]
        for j in idx:
            b = J.apply_to_operator(op.multiplier_matrices[j])
            worst = max(worst, float(np.linalg.norm(a @ b - b @ a, 2)))
    gen_norm = max(float(np.linalg.norm(op.multiplier_matrices[i], 2)) for i in idx)
    return worst, gen_norm


def test_order_zero_failure_counts_and_residuals_grow():
    """Counts and max residuals match the printed values and INCREASE with n_max."""
    rows = {n: _order_zero(n) for n in (2, 3, 4)}
    for n, (nf, nt, mr) in rows.items():
        pf, pt, pmr = PRINTED[n]
        assert nt == pt, f"n_max={n}: {nt} pairs, paper prints {pt}"
        assert nf == pf, f"n_max={n}: {nf} failing, paper prints {pf}"
        assert abs(mr - pmr) < 5e-4, f"n_max={n}: max residual {mr:.4e} vs printed {pmr}"
    res = [rows[n][2] for n in (2, 3, 4)]
    frac = [rows[n][0] / rows[n][1] for n in (2, 3, 4)]
    assert res[0] < res[1] < res[2], f"residuals do not grow: {res}"
    assert frac[0] < frac[1] < frac[2], f"failing fraction does not grow: {frac}"


def test_fixed_generator_residual_does_not_decay():
    """Degree-1 generators, n_max = 2, 3, 4: residual 0.051 -> 0.068 -> 0.076,
    increasing, and a fixed nonzero fraction (~0.21) of the generator norm."""
    vals = {}
    for n in (2, 3, 4):
        r, g = _fixed_generator_residual(n)
        vals[n] = (r, g)
        assert abs(r - FIXED_GEN_PRINTED[n]) < 1e-3, (
            f"n_max={n}: fixed-generator residual {r:.4f} vs printed {FIXED_GEN_PRINTED[n]}")
    rs = [vals[n][0] for n in (2, 3, 4)]
    assert rs[0] < rs[1] < rs[2], f"fixed-generator residual does not grow: {rs}"
    # not decaying to zero: bounded below by the n_max=2 value at every cutoff
    assert min(rs) >= 0.05
    # ratio to the (growing) generator norm stays near 0.21 through n_max=4
    # (still drifting; the limit is not established), not heading to 0
    ratio4 = vals[4][0] / vals[4][1]
    assert 0.19 < ratio4 < 0.23, f"residual/generator-norm at n_max=4 = {ratio4:.3f}"


@pytest.mark.slow
def test_fixed_generator_residual_n_max_5():
    """The printed n_max = 5 endpoint 0.081 (~17 s): still increasing, ratio ~0.21."""
    r4, _ = _fixed_generator_residual(4)
    r5, g5 = _fixed_generator_residual(5)
    assert abs(r5 - FIXED_GEN_PRINTED[5]) < 1e-3, f"n_max=5 residual {r5:.4f}"
    assert r5 > r4
    assert 0.19 < r5 / g5 < 0.23


# Paper 32 audit item 5 (model order one, full-Dirac, truthful CH):
# 43/196 failing at n_max = 2 (max residual 1.0e-1), 1337/3025 at n_max = 3
# (2.0e-1).  Paper 32 states both counts are regression-pinned here and that
# order one GROWS over n_max = 2, 3 (FULL #10, 2026-09-30: these were
# printed but never asserted, and n_max = 3 was never run).
ORDER_ONE_PRINTED = {2: (43, 196, 1.0e-1), 3: (1337, 3025, 2.0e-1)}


def _order_one(n_max: int):
    from geovac.full_dirac_operator_system import (
        FullDiracTruncatedOperatorSystem, camporesi_higuchi_full_dirac_matrix)
    from geovac.real_structure import build_J_full_dirac, verify_order_one
    op = FullDiracTruncatedOperatorSystem(n_max)
    D = camporesi_higuchi_full_dirac_matrix(op.basis)
    J = build_J_full_dirac(n_max)
    ok, failures = verify_order_one(J, op.multiplier_matrices, D, tol=1e-9)
    n_total = len(op.multiplier_matrices) ** 2
    max_res = max((r for _, _, r in failures), default=0.0)
    return len(failures), n_total, max_res


def test_order_one_counts_pinned_and_grow():
    """Model order-one counts match the printed 43/196 and 1337/3025, and both
    the failing fraction and the max residual INCREASE from n_max = 2 to 3.
    WRONG ANSWER REJECTED: 'order one fails only at the smallest cutoff and
    recovers' (a decaying count or residual)."""
    rows = {n: _order_one(n) for n in (2, 3)}
    for n, (nf, nt, mr) in rows.items():
        pf, pt, pmr = ORDER_ONE_PRINTED[n]
        assert nt == pt, f"n_max={n}: {nt} pairs, paper prints {pt}"
        assert nf == pf, f"n_max={n}: {nf} failing, paper prints {pf}"
        assert abs(mr - pmr) < 5e-3, f"n_max={n}: max residual {mr:.4e} vs printed {pmr}"
    assert rows[2][0] / rows[2][1] < rows[3][0] / rows[3][1]
    assert rows[2][2] < rows[3][2]
