"""Paper 60, sec:molecular, "[DERIVED + MEASURED] On the sigma basis no
preconditioner of that kind exists".

Backs the claims that, on the full sigma basis (n <= N, l < n, m = 0 per centre):

* the cross block C of the Shibuya--Wulfman metric fills its singular spectrum
  toward 1 as the basis grows (each centre's set is a unitary translate of the
  other's, so the limit [[I, U], [U*, I]] has spectrum {0, 2});
* the softness is between different l -- every fixed-l diagonal block stays
  well conditioned while the full metric does not;
* no centre-block-diagonal SPD preconditioner improves cond(S) (the
  (u, +-v) lemma);
* canonical orthogonalisation at threshold 1e-6 keeps the H2+ accuracy, while a
  fixed coarse threshold (1e-2) stalls it.

Wrong answers each test rejects are named in its docstring.
"""
from __future__ import annotations

import numpy as np
import pytest

import geovac.sturmian_sigma_basis as sb

H2P_EXACT_R2 = -1.102634214   # H2+ electronic energy at R = 2 bohr


def _cross_sv(N: int, kR: float) -> np.ndarray:
    S, _, labels = sb.build(N, 1.0, kR)
    m = len(labels)
    return np.linalg.svd(S[:m, m:], compute_uv=False)


def test_cross_block_spectrum_fills_toward_one() -> None:
    """Rejects "the sigma-basis near-dependence is a fixed number of directions"
    (the s-sector one-direction structure): the fraction of cross-block singular
    values above 0.9 must grow with N and pass 0.4 by 21 functions per centre."""
    fr = []
    for N in range(2, 7):
        sv = _cross_sv(N, 3.0)
        fr.append(float(np.mean(sv > 0.9)))
    assert fr[0] == 0.0
    assert all(b >= a for a, b in zip(fr, fr[1:])), fr
    assert fr[-1] > 0.4, fr


def test_softness_is_between_l_not_within() -> None:
    """Rejects "a per-l block treatment reaches the soft directions": in the
    gerade block at kR = 1, the worst fixed-l diagonal block grows less than
    10x from N = 3 to N = 5 (and stays < 60), while the full gerade block grows
    more than 1e6x over the same ladder."""
    def worst_and_full(N: int) -> tuple:
        S, _, labels = sb.build(N, 1.0, 1.0)
        g, _ = sb.gerade_ungerade(S, labels)
        ls = [l for (_, l) in labels]
        worst = 0.0
        for l in set(ls):
            idx = [i for i, x in enumerate(ls) if x == l]
            ev = np.linalg.eigvalsh(g[np.ix_(idx, idx)])
            worst = max(worst, ev[-1] / ev[0])
        ev = np.linalg.eigvalsh(g)
        return worst, ev[-1] / ev[0]
    w3, f3 = worst_and_full(3)
    w5, f5 = worst_and_full(5)
    assert w5 < 60.0 and w5 / w3 < 10.0, (w3, w5)
    assert f5 / f3 > 1e6, (f3, f5)


def test_centre_block_diagonal_preconditioner_cannot_help() -> None:
    """Rejects "a per-centre (block-Jacobi-type) preconditioner bounds the sigma
    metric": for random SPD P = diag(P_A, P_B), cond(P^-1 S) >= cond(S)."""
    S, _, labels = sb.build(4, 1.0, 3.0)
    m = len(labels)
    ev = np.linalg.eigvalsh(S)
    c0 = ev[-1] / ev[0]
    rng = np.random.default_rng(7)
    for _ in range(20):
        blocks = []
        for _c in range(2):
            A = rng.standard_normal((m, m))
            blocks.append(A @ A.T + 0.1 * np.eye(m))
        P = np.zeros_like(S)
        P[:m, :m], P[m:, m:] = blocks
        L = np.linalg.cholesky(P)
        Li = np.linalg.inv(L)
        evp = np.linalg.eigvalsh(Li @ S @ Li.T)
        assert evp[-1] / evp[0] >= c0 * (1 - 1e-9)


def _energy_at_tol(N: int, tol: float, monkeypatch) -> float:
    monkeypatch.setattr(sb, "CANON_TOL", tol)
    return sb.galerkin_energy(N, 2.0)


def test_canonical_orthogonalisation_keeps_accuracy(monkeypatch) -> None:
    """Rejects "dropping the soft directions costs the sigma-basis accuracy":
    at threshold 1e-6, N = 5 (15 functions per centre), the H2+ error stays
    below 1e-5 Ha, equal to the undropped value."""
    e6 = _energy_at_tol(5, 1e-6, monkeypatch) - H2P_EXACT_R2
    e12 = _energy_at_tol(5, 1e-12, monkeypatch) - H2P_EXACT_R2
    assert 0 < e6 < 1e-5, e6
    assert abs(e6 - e12) < 1e-7


def test_fixed_coarse_threshold_stalls(monkeypatch) -> None:
    """Rejects "any fixed threshold is safe": at threshold 1e-2 the N = 5 error
    is more than 10x the undropped error."""
    e2 = _energy_at_tol(5, 1e-2, monkeypatch) - H2P_EXACT_R2
    e12 = _energy_at_tol(5, 1e-12, monkeypatch) - H2P_EXACT_R2
    assert e2 > 10 * e12, (e2, e12)


def test_metric_spectrum_is_one_plus_minus_cross_singular_values() -> None:
    """Backs the structural step behind the {0, 2} limit: with identity intra-centre
    blocks, spec(S) = {1 +- sigma_k(C)} exactly, so S's spectrum is symmetric about 1
    and its softness is the cross block's near-unit singular values. Rejects a metric
    whose intra blocks are not the identity (e.g. full-V0 weighting), for which the
    pairing fails."""
    S, _, labels = sb.build(4, 1.0, 3.0)
    m = len(labels)
    sv = np.linalg.svd(S[:m, m:], compute_uv=False)
    pred = np.sort(np.concatenate([1 - sv, 1 + sv]))
    assert np.max(np.abs(np.linalg.eigvalsh(S) - pred)) < 1e-9
