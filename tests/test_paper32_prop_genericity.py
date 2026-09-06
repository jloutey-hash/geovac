"""Backing test for Paper 32 rem: after prop:propagation_number:
the propagation number prop = 2 is the GENERIC value for a *-closed unital
operator system, not a structural specificity of the GeoVac S^3 system.

This backs the 2026-09-03 DOWNGRADE of cor:structural_specificity (prop=2 was
retitled from "structural specificity" to "generic value").  FULL trunk cert
#7 (2026-09-06) flagged it as a load-bearing claim with no committed test
(the claims/code reviewers confirmed it TRUE by hand: 8/8 at dim 14 in M5,
3/3 at dim 55 in M14).  This makes that a regression.

The propagation number is computed INLINE by iterated-product span rank
(prop = smallest k with dim(V^k) = N^2), independent of
geovac.operator_system.propagation_number -- so this is a genuine cross-route
check of the CLAIM, not a re-run of the production routine.

Discrimination is built in via controls that must NOT return 2:
  - the full matrix algebra M_N (dim N^2) has prop = 1;
  - the diagonal subalgebra (commutative) never saturates, prop = -1.
A test that returned 2 for those would be vacuous; these assertions prove the
inline prop genuinely varies with the subspace.
"""

from __future__ import annotations

import numpy as np
import pytest

RNG = np.random.default_rng(20260906)


def _herm(N: int) -> np.ndarray:
    a = RNG.standard_normal((N, N)) + 1j * RNG.standard_normal((N, N))
    return a + a.conj().T


def _rank(mats: list[np.ndarray], tol: float = 1e-9) -> int:
    """dim of the complex span of a list of N x N matrices."""
    M = np.stack([m.reshape(-1) for m in mats])
    s = np.linalg.svd(M, compute_uv=False)
    return int((s > tol * max(1.0, s[0])).sum())


def _random_star_closed_unital(N: int, d: int) -> list[np.ndarray]:
    """A *-closed unital operator system of complex dimension d in M_N.

    The complex span of Hermitian matrices is automatically *-closed
    ((sum c_i H_i)^* = sum conj(c_i) H_i lies in the same complex span), and
    including the identity makes it unital.  We add random Hermitian
    generators until the span has dimension exactly d.
    """
    assert 1 <= d <= N * N
    basis = [np.eye(N, dtype=complex)]
    guard = 0
    while _rank(basis) < d:
        cand = basis + [_herm(N)]
        if _rank(cand) > _rank(basis):
            basis = cand
        guard += 1
        if guard > 20 * d:
            raise RuntimeError("could not reach dimension d")
    assert _rank(basis) == d
    return basis


def _reduce(mats: list[np.ndarray], N: int, tol: float = 1e-9) -> list[np.ndarray]:
    """Reduce a set of N x N matrices to an orthonormal basis of their span,
    so iterated products cannot blow up combinatorially."""
    M = np.stack([m.reshape(-1) for m in mats])
    _, s, Vh = np.linalg.svd(M, full_matrices=False)
    r = int((s > tol * max(1.0, s[0])).sum())
    return [Vh[t].reshape(N, N) for t in range(r)]


def _prop(basis: list[np.ndarray], N: int, max_k: int = 8) -> int:
    """Smallest k with dim(V^k) = N^2; -1 if not reached within max_k.

    The span is reduced to an independent basis each iteration (<= N^2 mats),
    so a non-terminating target cannot cause combinatorial blow-up.
    """
    target = N * N
    span = _reduce(list(basis), N)
    for k in range(1, max_k + 1):
        if len(span) == target:
            return k
        span = _reduce([a @ b for a in span for b in basis], N)
    return -1


def test_prop_is_generically_two_M5_dim14():
    """dim-14 operator system in M5: prop = 2 across seeds (claims-reviewer 8/8)."""
    N, d = 5, 14
    got = [_prop(_random_star_closed_unital(N, d), N) for _ in range(8)]
    assert all(p == 2 for p in got), got


def test_prop_is_generically_two_M14_dim55():
    """dim-55 operator system in M14: prop = 2 across seeds (claims-reviewer 3/3)."""
    N, d = 14, 55
    got = [_prop(_random_star_closed_unital(N, d), N) for _ in range(3)]
    assert all(p == 2 for p in got), got


def test_controls_do_not_return_two():
    """The discrimination controls: prop varies with the subspace, so the
    ==2 assertions above are not vacuous."""
    N = 5
    full = [m for m in np.eye(N * N)]  # placeholder, replaced below
    # full matrix algebra: a basis of matrix units -> dim N^2 -> prop 1
    units = []
    for i in range(N):
        for j in range(N):
            e = np.zeros((N, N), dtype=complex)
            e[i, j] = 1.0
            units.append(e)
    assert _rank(units) == N * N
    assert _prop(units, N) == 1
    # diagonal subalgebra: commutative, never saturates -> prop -1
    diag = [np.diag((np.arange(N) == k).astype(complex)) for k in range(N)]
    assert _rank(diag) == N
    assert _prop(diag, N) == -1
