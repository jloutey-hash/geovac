"""Operator-system backing for Paper 46 Lemma L3 (C_3^op = C_3 = 1, Paper 38 L3).

Replaces the withdrawn pure-sympy `test_paper46_C3op_closed_form.py`
(archived under tests/_archive/dead_ends/). The earlier draft asserted a
per-harmonic constant C_3^(N) = sqrt((N-1)/(N+1)) (operator-norm normalization
||[D_GV, M]|| <= C_3^(N) ||M||) with envelope sup C_3^op = sqrt(1 - 1/n_max),
attributed to Paper 38 L3. That is FALSE on the operator system (the /qa group1
cert5 diagnostic, debug/cert5_p46_c3op_diagnostic.py). This test freezes the
measured facts that motivated the withdrawal and verifies the surviving content:

  (a) the natural-substrate multipliers span shell-difference labels
      N <= 2*n_max - 1 (the Avery-Wen-Avery envelope RANGE refinement -- TRUE);
  (b) the per-harmonic operator-norm ratios ||[D_GV, M]|| / ||M|| do NOT realize
      sqrt((N-1)/(N+1)) -- they exceed it, and their supremum GROWS with n_max
      rather than approaching sqrt(1 - 1/n_max);
  (c) the envelope-max monopole Y^(3)_{2*n_max-1,0,0} COMMUTES with D_GV
      (||[D_GV, M]|| = 0 bit-exact) -- it is the loosest, not the tightest, generator,
      refuting the "asymptotically tight on the envelope-max harmonic" claim.

The actual L3 comparison constant is C_3 = 1 under the gradient/translation
seminorm (Paper 38 Lemma L3: ||[D_CH, M_f]|| = ||grad f||_inf on the continuum),
backed by tests/test_p38_action_seminorm.py. See Paper 46 Remark rem:env_tightness
and Appendix app:correction.
"""
from __future__ import annotations

import numpy as np
import pytest

from geovac.full_dirac_operator_system import (
    FullDiracTruncatedOperatorSystem,
    camporesi_higuchi_full_dirac_matrix,
)


def _opnorm(M: np.ndarray) -> float:
    return float(np.linalg.norm(M, ord=2))


def _comm(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    return A @ B - B @ A


def _ratios_by_N(n_max: int) -> dict:
    """Map degree N -> list of (L, M, ||[D,M]||, ||M||, ratio) over the substrate."""
    ops = FullDiracTruncatedOperatorSystem(n_max)
    D = camporesi_higuchi_full_dirac_matrix(ops.basis)
    by_N: dict = {}
    for (N, L, M), mat in zip(ops.multiplier_labels, ops.multiplier_matrices):
        lop = _opnorm(_comm(D, mat))
        mn = _opnorm(mat)
        by_N.setdefault(N, []).append((L, M, lop, mn, lop / mn if mn > 1e-14 else 0.0))
    return by_N


@pytest.mark.parametrize("n_max", [3, 4, 5])
def test_substrate_spans_awa_envelope(n_max: int) -> None:
    """(a) RANGE refinement stands: substrate reaches N = 2*n_max - 1."""
    by_N = _ratios_by_N(n_max)
    assert (2 * n_max - 1) in by_N, f"envelope-max N={2*n_max-1} absent from substrate"
    assert max(by_N) == 2 * n_max - 1


@pytest.mark.parametrize("n_max", [3, 4, 5])
def test_per_harmonic_opnorm_ratio_does_not_realize_sqrt(n_max: int) -> None:
    """(b) The withdrawn per-harmonic constant sqrt((N-1)/(N+1)) is violated."""
    by_N = _ratios_by_N(n_max)
    exceed = 0
    for N in by_N:
        if N < 2:
            continue
        c3N = np.sqrt((N - 1) / (N + 1))
        if max(r[4] for r in by_N[N]) > c3N + 1e-9:
            exceed += 1
    assert exceed >= 1, "op-norm ratios should exceed the withdrawn sqrt((N-1)/(N+1))"


@pytest.mark.parametrize("n_max", [3, 4, 5])
def test_sup_opnorm_ratio_nowhere_near_sqrt_envelope(n_max: int) -> None:
    """(b') sup op-norm ratio is nowhere near the withdrawn sqrt(1-1/n_max)."""
    by_N = _ratios_by_N(n_max)
    sup = max(r[4] for v in by_N.values() for r in v)
    assert sup > np.sqrt(1 - 1 / n_max) + 0.5


def test_sup_ratio_grows_with_n_max() -> None:
    """(b'') sup op-norm ratio GROWS with n_max (does not approach 1)."""
    sups = []
    for n_max in (3, 4, 5):
        by_N = _ratios_by_N(n_max)
        sups.append(max(r[4] for v in by_N.values() for r in v))
    assert sups[0] < sups[1] < sups[2], f"sup should grow, got {sups}"


@pytest.mark.parametrize("n_max", [3, 4, 5])
def test_envelope_max_monopole_commutes(n_max: int) -> None:
    """(c) envelope-max monopole commutes with D_GV (loosest generator)."""
    by_N = _ratios_by_N(n_max)
    Nenv = 2 * n_max - 1
    mono = [r for r in by_N[Nenv] if r[0] == 0 and r[1] == 0]
    assert mono, f"no L=0,M=0 multiplier at envelope-max N={Nenv}"
    lop = mono[0][2]
    assert lop < 1e-9, (
        f"envelope-max monopole N={Nenv} should commute with D_GV "
        f"(||[D,M]||=0), got {lop}"
    )
