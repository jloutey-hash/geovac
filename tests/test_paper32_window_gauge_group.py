"""Paper 32, scope paragraph after thm:h1_yukawa_non_selection: Dong's gauge group on the true window.

Record: docs/predictions/2026-10-31_window_gauge_group.md (WG, GO).  Written as a separate step after the
result (CLAUDE.md Sec. 9, guard-writing rule).  No debug/ import: the Gram construction is restated here, and
the operator system is the tracked Paper 38 guard's Route A (tests/test_p38_true_ch_window.py).

Dong (SIGMA 18 (2022) 060): G(E) = {U unitary : U E U^* = E}.  The Lie algebra of its identity component is
g = {X anti-Hermitian : [X, E] in E}, the null space of the PSD form G = sum_k L_k^dag (1 - P_E) L_k,
L_k X = [X, e_k], over an orthonormal basis e_k of E.
"""
import importlib.util
import itertools
import os
from fractions import Fraction

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location('p38_for_wg', os.path.join(HERE, 'test_p38_true_ch_window.py'))
p38 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(p38)


def onb(mats, tol=1e-9):
    V = np.array([m.reshape(-1) for m in mats])
    _, s, Vh = np.linalg.svd(V, full_matrices=False)
    N = mats[0].shape[0]
    return [Vh[i].reshape(N, N) for i in range(int(np.sum(s > tol * s[0])))]


def gauge_algebra_dim(E):
    N = E[0].shape[0]
    I = np.eye(N)
    G = np.kron(I, sum(e.conj() @ e.T for e in E)) + np.kron(sum(e.conj().T @ e for e in E), I)
    for e in E:
        G -= np.kron(e, e.conj()) + np.kron(e.conj().T, e.T)
    R = np.array([(e @ b.conj().T - b.conj().T @ e).T.reshape(-1) for e in E for b in E])
    G -= R.conj().T @ R
    w = np.linalg.eigvalsh(0.5 * (G + G.conj().T))
    top = w[-1]
    null = w[w < 1e-9 * top]
    rest = w[w >= 1e-9 * top]
    return len(null), rest.min() / top


def toeplitz(n):
    mats = []
    for k in range(-(n - 1), n):
        T = np.zeros((n, n), complex)
        for j in range(n):
            if 0 <= j + k < n:
                T[j + k, j] = 1.0
        mats.append(T)
    return onb(mats)


def window(n_max):
    _, vecs, Dc, dimH, _ = p38._window(n_max)
    return onb([v.reshape(dimH, dimH) for v in vecs]), Dc


def test_dongs_circle_is_reproduced_with_and_without_complex_conjugation():
    """Control: Dong's Cor. 4.3 gives an identity component U(1) x U(1) for Toep_n at every n, so dim g = 2;
    also after conjugating by a random unitary (complex matrices).

    Wrong answers rejected: a Gram form that drops the in-E projection and so counts only the commutant
    (dim 1).  (Conjugating the basis vectors is NOT a wrong answer for these systems: Toep_n and the true
    window's O are both closed under entrywise conjugation, measured to 5e-15, so it changes nothing.)
    """
    rng = np.random.default_rng(7)
    for n in (3, 5, 6):
        d, gap = gauge_algebra_dim(toeplitz(n))
        assert d == 2 and gap > 1e-2, (n, d, gap)
        W, _ = np.linalg.qr(rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)))
        d, gap = gauge_algebra_dim(onb([W @ e @ W.conj().T for e in toeplitz(n)]))
        assert d == 2 and gap > 1e-2, (n, d, gap)


def test_window_gauge_algebra_is_a_phase_plus_the_isometries():
    """WG1/WG2: on the true window dim g = 7 at n_max = 1, 2, 3 with a clean gap, and the seven known
    generators (identity, left translations, right translations with spin) lie in g, so
    g = u(1) + su(2)_L + su(2)_R.

    Wrong answers rejected: a surviving local gauge direction (dim > 7), a missing isometry (dim < 7), a
    right translation without the spin part (it does not commute with D, leaves the window, and is not in g).
    """
    for n_max in (1, 2, 3):
        E, Dc = window(n_max)
        d, gap = gauge_algebra_dim(E)
        assert d == 7 and gap > 0.05, (n_max, d, gap)
    labels, idx, D, ev, U = p38._build(2)
    keep = np.abs(ev) <= 2.5 + 1e-9
    V = U[:, keep]
    E, _ = window(2)
    B = np.array([e.reshape(-1) for e in E])
    n = len(labels)
    for a in range(3):
        for with_spin in (True, False):
            R = np.zeros((n, n), complex)
            for j in sorted({lab[0] for lab in labels}):
                dd = int(2 * j) + 1
                Jm = p38._spin_mats(j)[a]
                ms = [j - k for k in range(dd)]
                ids = [idx[(j, m, mp, s)] for m in ms for mp in ms for s in (0, 1)]
                sig = (p38.SX, p38.SY, p38.SZ)[a] / 2 if with_spin else 0 * p38.SX
                R[np.ix_(ids, ids)] = np.kron(np.eye(dd), np.kron(Jm, np.eye(2)) + np.kron(np.eye(dd), sig))
            commutes = bool(np.linalg.norm(R @ D - D @ R) < 1e-12)
            assert commutes is with_spin, (a, with_spin)
            if with_spin:
                X = V.conj().T @ R @ V
                for e in E:
                    c = (X @ e - e @ X).reshape(-1)
                    assert np.linalg.norm(c - B.T @ (B.conj() @ c)) < 1e-12


def test_inversion_is_implemented_by_a_unitary_that_anticommutes_with_D():
    """WG3: the inversion g -> g^{-1}, D^J_{MM'} -> (-1)^{M'-M} D^J_{-M',-M}, is implemented on the window by a
    unitary unique up to phase, which anticommutes with D (orientation reversal); the identity map gives
    the scalars only (control).

    Wrong answers rejected: an inversion that no unitary implements (null space 0), a reflection that
    commutes with D, and a solution space larger than one (an extra symmetry).
    """
    n_max = 2
    labels, idx, D, ev, Uev = p38._build(n_max)
    keep = np.abs(ev) <= n_max + 0.5 + 1e-9
    V = Uev[:, keep]
    Dc = V.conj().T @ D @ V
    keys = []
    for N in range(1, 2 * n_max + 1):
        J = Fraction(N - 1, 2)
        Ms = [J - k for k in range(N)]
        keys += [(J, M, Mp) for M, Mp in itertools.product(Ms, Ms)]
    comp = {k: V.conj().T @ p38._mult_op(*k, labels, idx) @ V for k in keys}
    N = V.shape[1]
    I = np.eye(N)

    def solve(iota):
        G = np.zeros((N * N, N * N), complex)
        for k, e in comp.items():
            sgn, k2 = iota(*k)
            L = np.kron(I, e.T) - np.kron(sgn * comp[k2], I)
            G += L.conj().T @ L
        w, v = np.linalg.eigh(0.5 * (G + G.conj().T))
        null = int(np.sum(w < 1e-9 * w[-1]))
        Us = v[:, 0].reshape(N, N)
        return null, Us / np.sqrt(np.trace(Us.conj().T @ Us).real / N)

    null, Uid = solve(lambda J, M, Mp: (1, (J, M, Mp)))
    assert null == 1 and np.linalg.norm(Uid @ Dc - Dc @ Uid) < 1e-10
    null, Ui = solve(lambda J, M, Mp: ((-1) ** int(Mp - M), (J, -Mp, -M)))
    assert null == 1, null
    assert np.linalg.norm(Ui.conj().T @ Ui - I) < 1e-10
    assert np.linalg.norm(Ui @ Dc + Dc @ Ui) < 1e-10
    assert np.linalg.norm(Ui @ Dc - Dc @ Ui) > 1.0


def _preserves(U, E):
    B = np.array([e.reshape(-1) for e in E])
    worst = 0.0
    for e in E:
        c = (U @ e @ U.conj().T).reshape(-1)
        worst = max(worst, np.linalg.norm(c - B.T @ (B.conj() @ c)) / np.linalg.norm(c))
    return worst


def test_full_group_inputs_commutant_scalars_and_centre_but_not_sign_D():
    """The inputs of the full-group reduction (record WG follow-up, debug/wg_full_group_memo.md, not
    pre-registered): the commutant of O is the scalars (so O generates M_N, its own C*-envelope), and among unitary
    functions of D the level-alternating one (a centre element of Spin(4)) preserves O while sign(D) does not.

    Wrong answers rejected: a commutant larger than the scalars (the envelope argument would fail), sign(D) taken
    as a gauge symmetry, and a test that cannot see a symmetry (the centre element must pass).
    """
    for n_max in (2, 3):
        E, Dc = window(n_max)
        N = Dc.shape[0]
        I = np.eye(N)
        G = sum((np.kron(I, e.T) - np.kron(e, I)).conj().T @ (np.kron(I, e.T) - np.kron(e, I)) for e in E)
        w = np.linalg.eigvalsh(0.5 * (G + G.conj().T))
        assert int(np.sum(w < 1e-9 * w[-1])) == 1
        ev, V = np.linalg.eigh(0.5 * (Dc + Dc.conj().T))
        vals = sorted(set(np.round(ev, 9)))
        P = {v: V[:, np.abs(ev - v) < 1e-6] @ V[:, np.abs(ev - v) < 1e-6].conj().T for v in vals}
        lev = sum(((-1) ** i) * P[v] for i, v in enumerate(vals))
        sgn = sum((1.0 if v > 0 else -1.0) * P[v] for v in vals)
        assert _preserves(lev, E) < 1e-10, _preserves(lev, E)
        assert _preserves(sgn, E) > 0.5, _preserves(sgn, E)
