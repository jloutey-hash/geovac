"""Paper 38, Proposition prop:gauge_window: gauge fields (inner fluctuations) on the true CH window.

D_w = D_CH (x) 1 + c(omega) on L^2(S) (x) C^k, truncated A_n = P D_w P (P = window |lambda| <= n + 1/2).
Construction reused from test_p38_true_ch_window (Clebsch-Gordan multipliers, D = 2 sigma.J_R + 3/2).

Checked:
(a) pure gauge u = pi_2 (spin J = 1/2): levels m <= n - 2 are eigenvalues of A_n with EXACTLY their
    multiplicity k (m+1)(m+2); the top level m = n - 1 is not reproduced;
(b) a position-dependent field omega = sum_a f_a theta^a (x) X_a: truncated eigenvalues converge to D_w's,
    within the proposition's factorial bound, and the gauge transform by pi_2 converges too;
(c) kernel of T -> [A_n, T] on O_n (x) M_2: scalars (generic), diagonal (abelian), scalars (pure gauge).

Wrong answers rejected (each named at its test):
- "truncation is gauge covariant" (the top level of the pure-gauge truncation would match);
- "the pure-gauge levels are reproduced only approximately" / "with extra copies" (exact count asserted);
- "a non-pure gauge field is reproduced exactly" (errors are nonzero but inside the factorial bound);
- "an abelian field metrizes the matrix factor" (kernel 2, not 1);
- "the truncated kernel equals the continuum kernel" (pure gauge: 1 on the window, 4 in the continuum).
"""
from fractions import Fraction
from functools import lru_cache
import itertools
import math

import numpy as np
import pytest

import test_p38_true_ch_window as W

H = Fraction(1, 2)
K = 2
SIG = [np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]], complex),
       np.diag([1, -1]).astype(complex)]


def _herm(Y):
    return 0.5 * (Y + Y.conj().T)


@lru_cache(maxsize=None)
def _space(nb):
    labels, idx, D, ev, Uev = W._build(nb)
    N = len(labels)
    sig = []
    for a in range(3):
        S = np.zeros((N, N), complex)
        for (j, m, mp, s), i in idx.items():
            for t in (0, 1):
                S[idx[(j, m, mp, t)], i] += SIG[a][t, s]
        sig.append(S)
    Mf = [_herm(W._mult_op(H, A, B, labels, idx)) for (A, B) in ((H, H), (H, -H), (-H, H))]
    Ms = [H, -H]
    U = np.zeros((K * N, K * N), complex)
    for a in range(2):
        for b in range(2):
            E = np.zeros((2, 2))
            E[a, b] = 1
            U += np.kron(W._mult_op(H, Ms[a], Ms[b], labels, idx), E)
    low = np.kron(Uev[:, np.abs(ev) <= 1.5 + 1e-9], np.eye(K))
    U /= np.sqrt(np.real(np.trace(low.conj().T @ U @ U.conj().T @ low)) / low.shape[1])
    return labels, idx, D, ev, Uev, sig, Mf, U


def _X(kind):
    rng = np.random.default_rng(5)
    if kind == "generic":
        return [_herm(rng.standard_normal((K, K)) + 1j * rng.standard_normal((K, K))) for _ in range(3)]
    return [np.diag(rng.standard_normal(K)).astype(complex) for _ in range(3)]


def _ops(nb, kind="generic"):
    _, _, D, ev, Uev, sig, Mf, U = _space(nb)
    Dk = np.kron(D, np.eye(K))
    X = _X(kind)
    Dw = Dk + _herm(sum(np.kron(Mf[a] @ sig[a], X[a]) for a in range(3)))
    return Dk, Dw, U, ev, Uev


def _trunc(Op, ev, Uev, n):
    V = np.kron(Uev[:, np.abs(ev) <= n + 0.5 + 1e-9], np.eye(K))
    return V, V.conj().T @ Op @ V


@pytest.mark.parametrize("n", [2, 3, 4])
def test_pure_gauge_levels_exact_below_the_boundary_layer(n):
    """(a): levels m <= n - 2 with exactly multiplicity k (m+1)(m+2). Rejects approximate reproduction and
    extra copies."""
    Dk, _, U, ev, Uev = _ops(n + 3)
    _, A = _trunc(U @ Dk @ U.conj().T, ev, Uev, n)
    lam = np.linalg.eigvalsh(_herm(A))
    for m in range(n - 1):
        for sgn in (1, -1):
            cnt = int(np.sum(np.abs(lam - sgn * (m + 1.5)) < 1e-9))
            assert cnt == K * (m + 1) * (m + 2)


@pytest.mark.parametrize("n", [1, 2, 3, 4])
def test_pure_gauge_top_level_is_broken(n):
    """(a), sharpness: the top level m = n - 1 is NOT reproduced. Rejects 'truncation is gauge covariant'."""
    Dk, _, U, ev, Uev = _ops(n + 3)
    _, A = _trunc(U @ Dk @ U.conj().T, ev, Uev, n)
    lam = np.linalg.eigvalsh(_herm(A))
    m = n - 1
    full = K * (m + 1) * (m + 2)
    worst = max(np.sort(np.abs(lam - s * (m + 1.5)))[full - 1] for s in (1, -1))
    assert 0.3 < worst <= 3.0 + 1e-9  # measured 1/3 .. 3 (paper: errors 1/3 to 3)


def _band(Op, ev, Uev):
    """Largest change of |lambda| that Op makes between D eigenvectors (within the build)."""
    V = np.kron(Uev, np.eye(K))
    lamk = np.repeat(ev, K)
    M = np.abs(V.conj().T @ Op @ V)
    i, j = np.nonzero(M > 1e-10)
    return float(np.max(np.abs(np.abs(lamk[i]) - np.abs(lamk[j]))))


def _errors(Dw, ev, Uev, ns):
    _, Aref = _trunc(Dw, ev, Uev, 7)
    ref = np.linalg.eigvalsh(_herm(Aref))
    ref = ref[np.abs(ref) < 2.6]
    out = {}
    for n in ns:
        _, A = _trunc(Dw, ev, Uev, n)
        lam = np.linalg.eigvalsh(_herm(A))
        out[n] = [(mu, float(np.min(np.abs(lam - mu)))) for mu in ref]
    return out


def test_band_limited_field_converges():
    """(b), measured: the unit-strength field's worst error on |mu| < 2.6 falls at every step n = 3..6, to
    below 1e-4 (reference: the n = 7 window)."""
    _, Dw, _, ev, Uev = _ops(9)
    errs = [max(e for _, e in v) for v in _errors(Dw, ev, Uev, range(3, 7)).values()]
    assert errs[-1] < 1e-4
    assert all(errs[i + 1] < errs[i] for i in range(len(errs) - 1))


def test_factorial_bound_holds_where_it_is_not_vacuous():
    """(b), the bound: on a field scaled by 0.3 (beta / b < 1, so eps_q < beta and the bound bites), every
    eigenvalue with |mu| < 2.6 is within eps_q / (1 - eps_q / beta) of spec A_n, n = 3..6; at least 20 checks
    are non-vacuous. The bound is far from tight (measured errors sit below (eps_q / beta)^8), so stronger
    bounds also pass; the wrong answer it rejects is 'the non-pure field is reproduced exactly' (fire-tested by
    planting eps ~ 0). At unit strength beta / b ~ 1.25 and the bound is vacuous for n <= 6."""
    Dk, Dw1, _, ev, Uev = _ops(9)
    Dw = Dk + 0.3 * (Dw1 - Dk)
    B = Dw - Dk
    beta = float(np.linalg.norm(B, 2))
    b = _band(B, ev, Uev)
    assert b == pytest.approx(2.0)
    assert beta < b
    checked = 0
    for n, rows in _errors(Dw, ev, Uev, range(3, 7)).items():
        for mu, e in rows:
            q = math.floor((n + 0.5 - abs(mu)) / b)
            eps = beta * (beta / b) ** q / math.factorial(q)
            if eps < beta:
                assert e <= eps / (1 - eps / beta) + 1e-12, (n, mu, e, eps)
                checked += 1
    assert checked >= 20


def test_gauge_transformed_field_converges_to_the_same_spectrum():
    """(b), covariance in the limit: A_n(omega^u) -> spec D_w, differing from A_n(omega) at finite n."""
    NB = 9
    _, Dw, U, ev, Uev = _ops(NB)
    Dwu = U @ Dw @ U.conj().T
    _, Aref = _trunc(Dw, ev, Uev, 7)
    ref = np.linalg.eigvalsh(_herm(Aref))
    ref = ref[np.abs(ref) < 2.6]
    errs = {}
    for n in (4, 6):
        _, A = _trunc(Dwu, ev, Uev, n)
        lam = np.linalg.eigvalsh(_herm(A))
        errs[n] = max(float(np.min(np.abs(lam - mu))) for mu in ref)
        _, A0 = _trunc(Dw, ev, Uev, n)
        lam0 = np.linalg.eigvalsh(_herm(A0))
        errs[(n, 0)] = max(float(np.min(np.abs(lam0 - mu))) for mu in ref)
    assert errs[6] < 1e-3
    assert errs[4] > 2 * errs[(4, 0)]  # finite-n non-covariance: the transform lags


def _op_basis(n, nb):
    labels, idx, D, ev, Uev = _space(nb)[:5]
    Vs = Uev[:, np.abs(ev) <= n + 0.5 + 1e-9]
    vecs = []
    for Nb in range(1, 2 * n + 3):
        J = Fraction(Nb - 1, 2)
        Ms = [J - q for q in range(Nb)]
        for M, Mp in itertools.product(Ms, Ms):
            vecs.append((Vs.conj().T @ W._mult_op(J, M, Mp, labels, idx) @ Vs).reshape(-1))
    vecs = np.array(vecs)
    r = int(np.linalg.matrix_rank(vecs, tol=1e-9))
    _, _, Vh = np.linalg.svd(vecs, full_matrices=False)
    d = Vs.shape[1]
    out = []
    for B in Vh[:r]:
        for p in range(K * K):
            e = np.zeros((K, K), complex)
            e[p // K, p % K] = 1
            out.append(np.kron(B.reshape(d, d), e))
    return out


def _kernel(basis, G, with_gap=False):
    M = np.array([(G @ B - B @ G).reshape(-1) for B in basis]).T
    sv = np.linalg.svd(M, compute_uv=False)
    kd = len(basis) - int((sv > 1e-8 * sv[0]).sum())
    return (kd, float(np.min(sv[sv > 1e-8 * sv[0]]))) if with_gap else kd


@pytest.mark.parametrize("n", [1, 2, 3])
def test_kernels_generic_abelian_pure(n):
    """(c): generic -> 1, abelian -> 2 (diagonal), pure gauge -> 1 (continuum: 4). Rejects 'an abelian field
    metrizes' and 'the truncated kernel equals the continuum kernel'."""
    nb = n + 2
    basis = _op_basis(n, nb)
    Dk, Dg, U, ev, Uev = _ops(nb, "generic")
    _, Da, _, _, _ = _ops(nb, "abelian")
    got = {}
    for name, Op in (("generic", Dg), ("abelian", Da), ("pure", U @ Dk @ U.conj().T)):
        _, A = _trunc(Op, ev, Uev, n)
        got[name] = _kernel(basis, A)
    assert got == {"generic": 1, "abelian": 2, "pure": 1}
    _, A = _trunc(U @ Dk @ U.conj().T, ev, Uev, n)
    gap = _kernel(basis, A, with_gap=True)[1]
    assert gap == pytest.approx({1: 2.0, 2: 4.0 / 3.0, 3: 5 ** 0.5 / 2}[n], abs=1e-9)


@pytest.mark.parametrize("n", [1, 2])
def test_commuting_constants_lie_in_the_kernel(n):
    """(c), proved half: 1 (x) b with b commuting with the field's values (b diagonal, abelian field)."""
    nb = n + 2
    _, Da, _, ev, Uev = _ops(nb, "abelian")
    V, A = _trunc(Da, ev, Uev, n)
    d = V.shape[1] // K
    T = np.kron(np.eye(d), np.diag([1.0, -2.0]))
    assert np.abs(A @ T - T @ A).max() < 1e-10
