"""Paper 38, Proposition prop:matrix_extension: the true CH window tensored with M_k.

Claims checked on the true Camporesi-Higuchi window (construction reused from test_p38_true_ch_window):
(a) the seminorm of G = SU(2)_L x SU(2)_R acting by Ad W (x) Ad pi_k(h) has kernel exactly the scalars on
    O_n (x) M_k (pi_k irreducible);
(c) the Dirac route max(||[D (x) 1, .]||, ||[1 (x) D_F, .]||) with M_k on C^k has kernel of dimension >= k;
and the structural input of (b): the window projection P_n and D commute with every generator of G.

Wrong answers rejected (each named at its test):
- "Paper 38's left-translation seminorm already metrizes the matrix factor" (its kernel is 1 (x) M_k);
- "SU(2)_R need not rotate the matrices" (right action without Ad pi_k: kernel 1 (x) M_k again);
- "a reducible matrix rotation suffices" (pi = pi_2 (+) pi_1 on M_3: kernel > 1);
- "the finite Dirac operator metrizes M_k" (kernel k, not 1);
- "the right action is plain right translation, spinor slot untouched" (then P_n and D are not invariant).
"""
from fractions import Fraction
import itertools

import numpy as np
import pytest

import test_p38_true_ch_window as W

SIG = [np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]], complex),
       np.diag([1, -1]).astype(complex)]


def _generators(n_max, spin_slot=True):
    labels, idx, D, ev, U = W._build(n_max)
    N = len(labels)
    Lg = [np.zeros((N, N), complex) for _ in range(3)]
    Rg = [np.zeros((N, N), complex) for _ in range(3)]
    for j in sorted({lab[0] for lab in labels}):
        d = int(2 * j) + 1
        Jm = W._spin_mats(j)
        ms = [j - k for k in range(d)]
        ids = [idx[(j, m, mp, s)] for m in ms for mp in ms for s in (0, 1)]
        I_d, I_2 = np.eye(d), np.eye(2)
        for a in range(3):
            Lg[a][np.ix_(ids, ids)] = np.kron(np.kron(Jm[a], I_d), I_2)
            Rg[a][np.ix_(ids, ids)] = np.kron(np.kron(I_d, Jm[a]), I_2)
            if spin_slot:
                Rg[a][np.ix_(ids, ids)] += np.kron(np.kron(I_d, I_d), SIG[a] / 2)
    V = U[:, np.abs(ev) <= n_max + 0.5 + 1e-9]
    return Lg, Rg, V, D


def _op_basis(n_max):
    _, vecs, Dc, dimH, _ = W._window(n_max)
    dim_O = int(np.linalg.matrix_rank(vecs, tol=1e-9))
    _, _, Vh = np.linalg.svd(vecs, full_matrices=False)
    return [b.reshape(dimH, dimH) for b in Vh[:dim_O]], Dc, dimH


def _tensor_basis(basisO, k):
    out = []
    for B in basisO:
        for a, b in itertools.product(range(k), range(k)):
            e = np.zeros((k, k), complex)
            e[a, b] = 1
            out.append(np.kron(B, e))
    return out


def _kernel_dim(basis, gens):
    M = np.array([np.concatenate([(G @ B - B @ G).reshape(-1) for G in gens]) for B in basis]).T
    return len(basis) - int(np.linalg.matrix_rank(M, tol=1e-8))


def _irrep(k):
    return W._spin_mats(Fraction(k - 1, 2))


def _setup(n_max, k, rep=None):
    Lg, Rg, V, _ = _generators(n_max)
    basisO, Dc, dimH = _op_basis(n_max)
    Lc = [V.conj().T @ G @ V for G in Lg]
    Rc = [V.conj().T @ G @ V for G in Rg]
    rep = _irrep(k) if rep is None else rep
    left = [np.kron(G, np.eye(k)) for G in Lc]
    right = [np.kron(G, np.eye(k)) + np.kron(np.eye(dimH), rep[a]) for a, G in enumerate(Rc)]
    return _tensor_basis(basisO, k), left, right, Dc, dimH


CASES = [(1, 2), (1, 3), (2, 2), (2, 3)]


@pytest.mark.parametrize("n_max,k", CASES)
def test_left_translations_alone_leave_the_matrix_factor_unmetrized(n_max, k):
    """Rejects 'Paper 38's seminorm already metrizes M_k': kernel is 1 (x) M_k, dim k^2."""
    basis, left, _, _, _ = _setup(n_max, k)
    assert _kernel_dim(basis, left) == k * k


@pytest.mark.parametrize("n_max,k", CASES)
def test_combined_action_kernel_is_the_scalars(n_max, k):
    """Part (a): SU(2)_L x SU(2)_R with Ad pi_k on the matrices is ergodic on O_n (x) M_k."""
    basis, left, right, _, _ = _setup(n_max, k)
    assert _kernel_dim(basis, left + right) == 1


@pytest.mark.parametrize("n_max", [1, 2])
def test_right_action_without_matrix_rotation_fails(n_max):
    """Rejects 'SU(2)_R need not rotate the matrices': kernel stays 1 (x) M_k."""
    k = 2
    basis, left, right0, _, _ = _setup(n_max, k, rep=[np.zeros((k, k), complex)] * 3)
    assert _kernel_dim(basis, left + right0) == k * k


@pytest.mark.parametrize("n_max", [1, 2])
def test_reducible_matrix_rotation_fails(n_max):
    """Rejects 'any rotation of the matrices will do': pi_2 (+) pi_1 on M_3 leaves a kernel > 1."""
    j2 = _irrep(2)
    red = [np.block([[j2[a], np.zeros((2, 1))], [np.zeros((1, 2)), np.zeros((1, 1))]]) for a in range(3)]
    basis, left, right, _, _ = _setup(n_max, 3, rep=red)
    assert _kernel_dim(basis, left + right) > 1


@pytest.mark.parametrize("n_max,k", CASES)
def test_finite_dirac_operator_cannot_metrize_the_factor(n_max, k):
    """Part (c). Rejects 'D_CH (x) 1 and 1 (x) D_F metrize O_n (x) M_k': kernel dimension k for generic D_F."""
    basis, _, _, Dc, dimH = _setup(n_max, k)
    rng = np.random.default_rng(7 + k)
    X = rng.standard_normal((k, k)) + 1j * rng.standard_normal((k, k))
    DF = X + X.conj().T
    kd = _kernel_dim(basis, [np.kron(Dc, np.eye(k)), np.kron(np.eye(dimH), DF)])
    assert kd == k


@pytest.mark.parametrize("n_max", [1, 2, 3])
def test_window_and_dirac_are_invariant_under_the_spinor_isometries(n_max):
    """Input of (b): [P_n, W] = [D, W] = 0 for left and right(+spin) generators. Rejects the right action
    as plain right translation with the spinor slot untouched (then D is not invariant)."""
    Lg, Rg, V, D = _generators(n_max)
    P = V @ V.conj().T
    for G in Lg + Rg:
        assert np.abs(P @ G - G @ P).max() < 1e-9
        assert np.abs(D @ G - G @ D).max() < 1e-9
    _, Rbad, _, _ = _generators(n_max, spin_slot=False)
    assert max(np.abs(D @ G - G @ D).max() for G in Rbad) > 0.1


@pytest.mark.parametrize("n_max", [1, 2])
def test_spinor_isometries_preserve_the_operator_system(n_max):
    """Spin-blindness: [G, O_n] lies in O_n for every left and right(+spin) generator, so tau preserves
    O_n (x) M_k. Rejects a 'right action' that leaves the operator system (fire-tested with a random one)."""
    Lg, Rg, V, _ = _generators(n_max)
    basisO, _, _ = _op_basis(n_max)
    gens = [V.conj().T @ G @ V for G in Lg + Rg]
    span = np.array([B.reshape(-1) for B in basisO])
    dim_O = int(np.linalg.matrix_rank(span, tol=1e-9))
    for G in gens:
        comm = np.array([(G @ B - B @ G).reshape(-1) for B in basisO])
        assert int(np.linalg.matrix_rank(np.vstack([span, comm]), tol=1e-8)) == dim_O


@pytest.mark.parametrize("n_max", [1, 2])
def test_degenerate_finite_dirac_kernel_is_sum_of_squared_multiplicities(n_max):
    """Part (c) with D_F = diag(1, 1, 2) on C^3: kernel dim = 2^2 + 1^2 = 5, not k = 3 and not the
    number of distinct eigenvalues (2). Rejects the 'spectral projections' count."""
    k = 3
    basis, _, _, Dc, dimH = _setup(n_max, k)
    DF = np.diag([1.0, 1.0, 2.0]).astype(complex)
    kd = _kernel_dim(basis, [np.kron(Dc, np.eye(k)), np.kron(np.eye(dimH), DF)])
    assert kd == 5
