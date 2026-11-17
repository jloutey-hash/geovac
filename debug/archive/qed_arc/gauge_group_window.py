"""Dong's gauge group of the true Camporesi-Higuchi window operator system (thread 1, 2026-10-31).

Dong (SIGMA 18 (2022) 060) defines G(E) = {U unitary in C*(E) : U E U^* = E}.  For the true
window, C*(O) = M_N (Verhoeven 2026, Thm 4.10).  The Lie algebra of the identity component is
    g(E) = {X anti-Hermitian : [X, e] in E for every e in E},
a LINEAR condition.  Since E is *-closed, the complex solution space of [X, E] in E has complex
dimension equal to the real dimension of g(E).  It is the null space of the PSD Gram form
    G = sum_k L_k^dag (1 - P_E) L_k,     L_k : X -> [X, e_k],
with {e_k} an orthonormal basis of E and P_E the orthogonal projector onto E in M_N.

    python debug/gauge_group_window.py control      # Toeplitz system, Dong's case (dim g = 2)
    python debug/gauge_group_window.py inclusion    # the 7 isometry/phase generators lie in g
    python debug/gauge_group_window.py score        # dim g on the true window (after registration)

O^CH is built by the tracked Paper 38 guard's Route A (tests/test_p38_true_ch_window.py).
"""
import importlib.util
import itertools
import json
import os
import sys
import time
from fractions import Fraction

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(HERE, 'data', 'gauge_group_window.jsonl')
_spec = importlib.util.spec_from_file_location('p38w', os.path.join(ROOT, 'tests', 'test_p38_true_ch_window.py'))
p38 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(p38)


def orthonormal_basis(mats, tol=1e-9):
    """Orthonormal (Frobenius) basis of span(mats), as N x N matrices."""
    V = np.array([m.reshape(-1) for m in mats])
    _, s, Vh = np.linalg.svd(V, full_matrices=False)
    r = int(np.sum(s > tol * s[0]))
    N = mats[0].shape[0]
    return [Vh[i].reshape(N, N) for i in range(r)]


def gram(E):
    """G with null space = {X : [X, E] in E}; row-major vec, X -> vec(X)."""
    N = E[0].shape[0]
    I = np.eye(N)
    S_a = sum(e.conj() @ e.T for e in E)          # sum conj(e) e^T
    S_b = sum(e.conj().T @ e for e in E)          # sum e^dag e
    G = np.kron(I, S_a) + np.kron(S_b, I)
    for e in E:
        G -= np.kron(e, e.conj()) + np.kron(e.conj().T, e.T)
    # subtract the in-E part: rows r_{kl} with <b_l, [X, e_k]> = Tr(C X), C = e_k b_l^dag - b_l^dag e_k
    rows = []
    for e in E:
        for b in E:
            bd = b.conj().T
            C = e @ bd - bd @ e
            rows.append(C.T.reshape(-1))
            if len(rows) >= 4096:
                R = np.array(rows)
                G -= R.conj().T @ R
                rows = []
    if rows:
        R = np.array(rows)
        G -= R.conj().T @ R
    return 0.5 * (G + G.conj().T)


def null_dim(G, rel=1e-9):
    w = np.linalg.eigvalsh(G)
    top = max(abs(w[-1]), 1.0)
    gap_lo = w[w < rel * top]
    gap_hi = w[w >= rel * top]
    return len(gap_lo), (float(gap_lo.max()) if len(gap_lo) else None), (float(gap_hi.min()) if len(gap_hi) else None)


def residual(X, E):
    """|| (1 - P_E) [X, e] || summed over the basis, relative to ||[X, e]||."""
    B = np.array([e.reshape(-1) for e in E])
    num = den = 0.0
    for e in E:
        c = (X @ e - e @ X).reshape(-1)
        proj = B.T @ (B.conj() @ c)
        num += np.linalg.norm(c - proj) ** 2
        den += np.linalg.norm(c) ** 2
    if den < 1e-24:                  # X commutes with all of E (the identity): trivially in g
        return 0.0, 0.0
    return float(np.sqrt(num / den)), float(np.sqrt(den))


def toeplitz_system(n):
    """Toep_n: compressions of e^{ik theta}, |k| <= n-1, to span{e^{ij theta}, j < n}."""
    mats = []
    for k in range(-(n - 1), n):
        T = np.zeros((n, n), complex)
        for j in range(n):
            if 0 <= j + k < n:
                T[j + k, j] = 1.0
        mats.append(T)
    return orthonormal_basis(mats)


def window_system(n_max):
    _, vecs, Dc, dimH, ev = p38._window(n_max)
    mats = [v.reshape(dimH, dimH) for v in vecs]
    return orthonormal_basis(mats), Dc, dimH


def isometry_generators(n_max):
    """Identity, left translations (spin-j on m), right translations with spin (J_R on m' + sigma/2),
    each restricted to the window (they commute with D, so they preserve it)."""
    labels, idx, D, ev, U = p38._build(n_max)
    keep = np.abs(ev) <= n_max + 0.5 + 1e-9
    V = U[:, keep]
    n = len(labels)
    gens = {'identity': np.eye(n, dtype=complex)}
    for a, name in enumerate('xyz'):
        L = np.zeros((n, n), complex)
        R = np.zeros((n, n), complex)
        for j in sorted({lab[0] for lab in labels}):
            d = int(2 * j) + 1
            Jm = p38._spin_mats(j)[a]
            ms = [j - k for k in range(d)]
            ids = [idx[(j, m, mp, s)] for m in ms for mp in ms for s in (0, 1)]
            sig = (p38.SX, p38.SY, p38.SZ)[a] / 2
            L[np.ix_(ids, ids)] = np.kron(Jm, np.eye(2 * d))
            R[np.ix_(ids, ids)] = np.kron(np.eye(d), np.kron(Jm, np.eye(2)) + np.kron(np.eye(d), sig))
        gens['L' + name] = L
        gens['R' + name] = R
    comm_D = {k: float(np.linalg.norm(g @ D - D @ g)) for k, g in gens.items()}
    return {k: 1j * (V.conj().T @ g @ V) for k, g in gens.items()}, comm_D


def write(row):
    print(json.dumps(row), flush=True)
    with open(OUT, 'a') as f:
        f.write(json.dumps(row) + '\n')


def control():
    for n in range(2, 8):
        E = toeplitz_system(n)
        t = time.time()
        d, lo, hi = null_dim(gram(E))
        write({'case': 'toeplitz', 'n': n, 'dim_E': len(E), 'dim_g': d, 'null_max': lo, 'gap_min': hi,
               'dong_expects': 2, 't': round(time.time() - t, 2)})


def control_complex():
    """Same Toeplitz systems conjugated by a random unitary: complex matrices, dim g must stay 2.
    A second control of the Gram form on complex matrices.  It is NOT a test of basis conjugation: these
    systems (and the true window's O) are closed under entrywise conjugation, so that change is invisible
    here (record WG, erratum)."""
    rng = np.random.default_rng(20261031)
    for n in range(2, 8):
        Z = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
        W, _ = np.linalg.qr(Z)
        E = orthonormal_basis([W @ e @ W.conj().T for e in toeplitz_system(n)])
        d, lo, hi = null_dim(gram(E))
        write({'case': 'toeplitz_conjugated', 'n': n, 'dim_E': len(E), 'dim_g': d, 'null_max': lo,
               'gap_min': hi, 'dong_expects': 2})


def inclusion():
    for n_max in (1, 2, 3):
        E, Dc, dimH = window_system(n_max)
        gens, cD = isometry_generators(n_max)
        res = {k: residual(X, E)[0] for k, X in gens.items()}
        span = int(np.linalg.matrix_rank(np.array([X.reshape(-1) for X in gens.values()]), tol=1e-9))
        write({'case': 'inclusion', 'n_max': n_max, 'dimH': dimH, 'dim_E': len(E),
               'max_rel_residual': max(res.values()), 'commutes_with_D': max(cD.values()),
               'rank_of_7_generators': span})


def score(levels):
    for n_max in levels:
        t = time.time()
        E, Dc, dimH = window_system(n_max)
        d, lo, hi = null_dim(gram(E))
        gens, _ = isometry_generators(n_max)
        Gs = np.array([X.reshape(-1) for X in gens.values()])
        write({'case': 'window', 'n_max': n_max, 'dimH': dimH, 'dim_E': len(E), 'dim_g': d,
               'null_max': lo, 'gap_min': hi, 'rank_of_7_generators': int(np.linalg.matrix_rank(Gs, tol=1e-9)),
               't': round(time.time() - t, 1)})


def reflection(n_max):
    """WG3 (reported): is there a unitary U on the window with U A U^* = iota(A) for the
    inversion g -> g^{-1}?  D^J_{MM'}(g^{-1}) = conj D^J_{M'M}(g) = (-1)^{M'-M} D^J_{-M',-M}(g)
    (Condon-Shortley).  Linear in U: the null space of U -> U e - iota(e) U over the basis.
    Candidate sign conventions are all tried and reported, with the identity map as a control
    (it must give the scalars: null dim 1)."""
    labels, idx, D, ev, Uev = p38._build(n_max)
    keep = np.abs(ev) <= n_max + 0.5 + 1e-9
    V = Uev[:, keep]
    pairs = []
    for N in range(1, 2 * n_max + 1):
        J = Fraction(N - 1, 2)
        Ms = [J - k for k in range(N)]
        for M, Mp in itertools.product(Ms, Ms):
            pairs.append((J, M, Mp))
    comp = {k: V.conj().T @ p38._mult_op(*k, labels, idx) @ V for k in pairs}
    cands = {
        'identity (control)': lambda J, M, Mp: (1, (J, M, Mp)),
        'inversion, CS sign': lambda J, M, Mp: ((-1) ** int(Mp - M), (J, -Mp, -M)),
        'inversion, no sign': lambda J, M, Mp: (1, (J, -Mp, -M)),
        'index swap': lambda J, M, Mp: (1, (J, Mp, M)),
    }
    out = {}
    N = V.shape[1]
    I = np.eye(N)
    for name, f in cands.items():
        G = np.zeros((N * N, N * N), complex)
        for k, e in comp.items():
            sgn, k2 = f(*k)
            ie = sgn * comp[k2]
            Lk = np.kron(I, e.T) - np.kron(ie, I)      # vec_row(U e - ie U)
            G += Lk.conj().T @ Lk
        wv, vv = np.linalg.eigh(0.5 * (G + G.conj().T))
        top = wv[-1]
        null = int(np.sum(wv < 1e-9 * top))
        unitary_dev = None
        if null:
            Usol = vv[:, 0].reshape(N, N)
            Usol = Usol / np.sqrt(np.trace(Usol.conj().T @ Usol).real / N)
            unitary_dev = float(np.linalg.norm(Usol.conj().T @ Usol - I))
            anti = float(np.linalg.norm(Usol @ p38._window(n_max)[2] + p38._window(n_max)[2] @ Usol))
            comm = float(np.linalg.norm(Usol @ p38._window(n_max)[2] - p38._window(n_max)[2] @ Usol))
        out[name] = {'null_dim': null, 'unitary_dev': unitary_dev,
                     'D_anticommutator': anti if null else None, 'D_commutator': comm if null else None}
    write({'case': 'reflection', 'n_max': n_max, **out})


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'control':
        control()
    elif cmd == 'control_complex':
        control_complex()
    elif cmd == 'inclusion':
        inclusion()
    elif cmd == 'reflection':
        for n in [int(x) for x in sys.argv[2:]] or [1, 2]:
            reflection(n)
    elif cmd == 'score':
        score([int(x) for x in sys.argv[2:]] or [1, 2, 3])
