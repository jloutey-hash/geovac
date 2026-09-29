"""Paper 38 true-window freeze (owed by the 2026-09-28 /qa trunk FULL #9).

Port of debug/p38_true_ch_window.py (an independent build: no geovac import).
Paper 38 (papers/group1_operator_algebras/paper_38_su2_propinquity_convergence.tex):
  * lem:band_injectivity proof + Sec. "Substrate: the true window": on the
    true Camporesi-Higuchi window (complete Dirac eigenspaces, |lambda| <=
    n_max + 1/2) the compression of the band-N multipliers has rank N^2 for
    every N <= 2 n_max and rank 0 above; dim O_{n_max} = 5, 30, 91 at
    n_max = 1, 2, 3.
  * rem:dirac_degeneracy: the kernel of a -> [D, a] on O_{n_max} is exactly
    the scalars (dim 1, and 1 in O) at n_max = 1, 2, 3.
  * the Dirac spectrum is +-(n + 3/2), n >= 0, multiplicity (n+1)(n+2) each sign.

WRONG ANSWERS REJECTED:
  * the doubled-model envelope N <= 2 n_max - 1 (the top band N = 2 n_max
    missing), which is what the frozen test_p38_action_seminorm.py sees;
  * a Dirac-commutator kernel larger than the scalars (the doubled model's
    10/14 at n_max = 2).

Realization: L^2(SU(2)) (x) C^2 with basis sqrt(2j+1) D^j_{m m'} (x) e_s;
D = 2 sigma . J_R + 3/2 (J_R = spin-j matrices on m', identity on m);
multiplication by the spin-J matrix coefficient D^J_{M M'} acts by
Clebsch-Gordan coefficients on (m, m') (identity on spin).

Fire-test plant: the compression T = V^dag A V has its sign-mixing blocks
(between positive and negative Dirac levels) zeroed -- the doubled model's
defect -> the band-rank and kernel tests FIRE.

SECOND, INDEPENDENT CONSTRUCTION (Route B, /qa trunk DELTA #14 guard pass).
Paper 38 states that two independent constructions agree and are frozen here;
Route A above is the Clebsch-Gordan build.  Route B (ported from the /qa
reviewer's probe) uses NO Clebsch-Gordan code: the matrix coefficients
sqrt(2j+1) D^j_{ab}(g) are SAMPLED at 600 random SU(2) points (Wigner-D via
the matrix exponential of its own spin matrices), and each multiplication
operator is obtained by least-squares expansion of the sampled pointwise
product in the sampled basis.  It carries its own spin matrices, its own
Dirac assembly, and its own window cut, so a defect in one route leaves the
other intact.  Asserted for n_max = 1, 2 by both routes: dim O = 5, 30; band
ranks N^2 for N <= 2 n_max and 0 above; kernel of [D, .] on O = the scalars;
and the two routes agree band by band.
Fire-test plant (Route B only): its window cut n_max + 1/2 -> n_max - 1/2
-> the Route-B and cross-route tests FIRE while Route A stays green.
"""
from __future__ import annotations

import itertools
from fractions import Fraction
from functools import lru_cache
from typing import Dict, List, Tuple

import numpy as np
import pytest
from sympy import Rational
from sympy.physics.wigner import clebsch_gordan

SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]], complex)
SZ = np.array([[1, 0], [0, -1]], complex)

DIM_O_PRINTED = {1: 5, 2: 30, 3: 91}


@lru_cache(maxsize=None)
def _cg(j1: Fraction, m1: Fraction, j2: Fraction, m2: Fraction, J: Fraction, M: Fraction) -> float:
    if m1 + m2 != M:
        return 0.0
    R = lambda x: Rational(x.numerator, x.denominator)
    return float(clebsch_gordan(R(j1), R(j2), R(J), R(m1), R(m2), R(M)))


def _spin_mats(j: Fraction):
    ms = [j - k for k in range(int(2 * j) + 1)]
    d = len(ms)
    jz = np.diag([float(m) for m in ms]).astype(complex)
    jp = np.zeros((d, d), complex)
    for a, m in enumerate(ms):
        if a > 0:
            jp[a - 1, a] = np.sqrt(float(j * (j + 1) - m * (m + 1)))
    return (jp + jp.conj().T) / 2, (jp - jp.conj().T) / (2j), jz


@lru_cache(maxsize=None)
def _build(n_max: int):
    js = [Fraction(k, 2) for k in range(0, n_max + 1)]
    labels: List[Tuple[Fraction, Fraction, Fraction, int]] = []
    for j in js:
        ms = [j - k for k in range(int(2 * j) + 1)]
        for m, mp, s in itertools.product(ms, ms, (0, 1)):
            labels.append((j, m, mp, s))
    idx = {lab: i for i, lab in enumerate(labels)}
    D = np.zeros((len(labels), len(labels)), complex)
    for j in js:
        d = int(2 * j) + 1
        jx, jy, jz = _spin_mats(j)
        blk_r = 2 * (np.kron(jx, SX) + np.kron(jy, SY) + np.kron(jz, SZ)) + 1.5 * np.eye(2 * d)
        blk = np.kron(np.eye(d), blk_r)
        ms = [j - k for k in range(d)]
        ids = [idx[(j, m, mp, s)] for m in ms for mp in ms for s in (0, 1)]
        D[np.ix_(ids, ids)] = blk
    ev, U = np.linalg.eigh(D)
    return labels, idx, D, ev, U


def _mult_op(J: Fraction, M: Fraction, Mp: Fraction, labels, idx) -> np.ndarray:
    A = np.zeros((len(labels), len(labels)), complex)
    for (j1, m1, m1p, s) in labels:
        i1 = idx[(j1, m1, m1p, s)]
        for (j2, m2, m2p, s2) in labels:
            if s2 != s or m2 != m1 + M or m2p != m1p + Mp:
                continue
            if not (abs(j1 - J) <= j2 <= j1 + J):
                continue
            c = _cg(j1, m1, J, M, j2, m2) * _cg(j1, m1p, J, Mp, j2, m2p)
            if c != 0.0:
                A[idx[(j2, m2, m2p, s2)], i1] = np.sqrt(float((2 * j1 + 1) / (2 * j2 + 1))) * c
    return A


@lru_cache(maxsize=None)
def _window(n_max: int):
    """Band ranks, compressed operators, and the compressed Dirac on the window."""
    labels, idx, D, ev, U = _build(n_max)
    keep = np.abs(ev) <= n_max + 0.5 + 1e-9
    V = U[:, keep]
    ev_kept = ev[keep]
    Dc = V.conj().T @ D @ V
    ranks: Dict[int, int] = {}
    vecs_all = []
    for N in range(1, 2 * n_max + 3):
        J = Fraction(N - 1, 2)
        Ms = [J - k for k in range(N)]
        vecs = []
        for M, Mp in itertools.product(Ms, Ms):
            A = _mult_op(J, M, Mp, labels, idx)
            T = V.conj().T @ A @ V
            vecs.append(T.reshape(-1))
        ranks[N] = int(np.linalg.matrix_rank(np.array(vecs), tol=1e-9))
        vecs_all.extend(vecs)
    return ranks, np.array(vecs_all), Dc, V.shape[1], ev_kept


@pytest.mark.parametrize("n_max", [1, 2, 3])
def test_dirac_spectrum_and_multiplicities(n_max):
    """Spectrum +-(n + 3/2), multiplicity (n+1)(n+2) per sign, on the window."""
    _, _, _, ev, _ = _build(n_max)
    keep = np.abs(ev) <= n_max + 0.5 + 1e-9
    kept = np.round(ev[keep], 9)
    for n in range(n_max):
        for sgn in (+1, -1):
            mult = int(np.sum(np.abs(kept - sgn * (n + 1.5)) < 1e-7))
            assert mult == (n + 1) * (n + 2), f"n_max={n_max}: {sgn*(n+1.5)} x{mult}"
    assert len(kept) == 2 * sum((n + 1) * (n + 2) for n in range(n_max))
    # nothing on the window outside the CH ladder
    ladder = {sgn * (n + 1.5) for n in range(n_max) for sgn in (1, -1)}
    assert all(min(abs(x - y) for y in ladder) < 1e-7 for x in kept)


@pytest.mark.parametrize("n_max", [1, 2, 3])
def test_band_ranks_full_up_to_2nmax_zero_above(n_max):
    """rank N^2 for N <= 2 n_max (INCLUDING the top band N = 2 n_max, which the
    doubled model misses), 0 for N = 2 n_max + 1, 2 n_max + 2; dim O = 5, 30, 91."""
    ranks, vecs, _, _, _ = _window(n_max)
    for N, r in ranks.items():
        expect = N * N if N <= 2 * n_max else 0
        assert r == expect, f"n_max={n_max}, band N={N}: rank {r}, expected {expect}"
    assert ranks[2 * n_max] == (2 * n_max) ** 2 > 0  # the top band is present
    dim_O = int(np.linalg.matrix_rank(vecs, tol=1e-9))
    assert dim_O == DIM_O_PRINTED[n_max] == sum(N * N for N in range(1, 2 * n_max + 1))


@pytest.mark.parametrize("n_max", [1, 2, 3])
def test_dirac_commutator_kernel_is_the_scalars(n_max):
    """ker(a -> [D, a]) on O is exactly C.1: kernel dim 1 and 1 in O."""
    _, vecs, Dc, dimH, _ = _window(n_max)
    dim_O = int(np.linalg.matrix_rank(vecs, tol=1e-9))
    _, _, Vh = np.linalg.svd(vecs, full_matrices=False)
    basis = Vh[:dim_O]
    comm = np.array([(Dc @ b.reshape(dimH, dimH) - b.reshape(dimH, dimH) @ Dc).reshape(-1)
                     for b in basis])
    rc = int(np.linalg.matrix_rank(comm, tol=1e-9))
    assert dim_O - rc == 1, f"n_max={n_max}: kernel dim {dim_O - rc} (scalars = 1)"
    one = np.eye(dimH).reshape(-1)
    assert int(np.linalg.matrix_rank(np.vstack([basis, one]), tol=1e-9)) == dim_O, "1 not in O"


# ---------------------------------------------------------------------------
# Route B: sampled-function / Wigner-D least-squares build (no Clebsch-Gordan)
# ---------------------------------------------------------------------------
_B_PAULI = (np.array([[0, 1], [1, 0]], complex),
            np.array([[0, -1j], [1j, 0]], complex),
            np.diag([1, -1]).astype(complex))
_B_SAMPLES = 600


def _b_spin(j2: int):
    """Spin-(j2/2) matrices (J_x, J_y, J_z), built independently of Route A."""
    j = j2 / 2
    ms = [j - k for k in range(j2 + 1)]
    jz = np.diag(ms).astype(complex)
    jp = np.zeros((j2 + 1, j2 + 1), complex)
    for a, m in enumerate(ms):
        if a > 0:
            jp[a - 1, a] = np.sqrt(j * (j + 1) - m * (m + 1))
    return [(jp + jp.conj().T) / 2, (jp - jp.conj().T) / (2j), jz]


def _b_samples():
    rng = np.random.default_rng(7)
    out = []
    for _ in range(_B_SAMPLES):
        q = rng.standard_normal(4)
        q /= np.linalg.norm(q)
        out.append(np.array([[q[0] + 1j * q[3], q[2] + 1j * q[1]],
                             [-q[2] + 1j * q[1], q[0] - 1j * q[3]]]))
    return out


def _b_wigner(j2: int, g: np.ndarray) -> np.ndarray:
    """D^{j2/2}(g): write g = exp(-i H), H = sum_a c_a sigma_a / 2, and
    exponentiate the same rotation in the spin-(j2/2) representation."""
    from scipy.linalg import expm
    w, V = np.linalg.eig(g)
    H = V @ np.diag(-np.angle(w)) @ np.linalg.inv(V)
    c = [np.trace(H @ s_).real for s_ in _B_PAULI]
    J = _b_spin(j2)
    return expm(-1j * sum(ci * Ji for ci, Ji in zip(c, J)))


@lru_cache(maxsize=None)
def _route_b(n_max: int):
    gs = _b_samples()
    # sampled basis up to j2 = 3 n_max + 1 (enough for every product below)
    cols, labs = [], []
    Dcache = {}
    for j2 in range(3 * n_max + 2):
        Ds = np.array([_b_wigner(j2, g) for g in gs])
        Dcache[j2] = Ds
        for a in range(j2 + 1):
            for b in range(j2 + 1):
                cols.append(np.sqrt(j2 + 1) * Ds[:, a, b])
                labs.append((j2, a, b))
    F = np.array(cols).T
    H_idx = [i for i, l in enumerate(labs) if l[0] <= n_max]
    Hlabs = [labs[i] for i in H_idx]
    dimf = len(Hlabs)
    full = [(l, s_) for l in Hlabs for s_ in (0, 1)]
    pos = {f: i for i, f in enumerate(full)}
    Dop = np.zeros((len(full), len(full)), complex)
    for (j2, a, b), s1 in full:
        J = _b_spin(j2)
        for (j2p, a2, b2), s2 in full:
            if j2p != j2 or a2 != a:
                continue
            val = sum(2 * J[k][b2, b] * _B_PAULI[k][s2, s1] for k in range(3))
            if b2 == b and s2 == s1:
                val += 1.5
            Dop[pos[((j2p, a2, b2), s2)], pos[((j2, a, b), s1)]] = val
    ev_b, U_b = np.linalg.eigh(Dop)
    keep_b = np.abs(ev_b) <= n_max + 0.5 + 1e-9
    Vb = U_b[:, keep_b]
    Dc_b = Vb.conj().T @ Dop @ Vb
    Fp = np.linalg.pinv(F)
    ranks: Dict[int, int] = {}
    allv = []
    for N in range(1, 2 * n_max + 3):
        Ds = Dcache[N - 1]
        vecs = []
        for M, Mp in itertools.product(range(N), range(N)):
            fvals = Ds[:, M, Mp]
            A = np.zeros((dimf, dimf), complex)
            for c, i in enumerate(H_idx):
                coef = Fp @ (fvals * F[:, i])
                A[:, c] = coef[H_idx]
            T = Vb.conj().T @ np.kron(A, np.eye(2)) @ Vb
            vecs.append(T.reshape(-1))
        ranks[N] = int(np.linalg.matrix_rank(np.array(vecs), tol=1e-7))
        allv.extend(vecs)
    return ranks, np.array(allv), Dc_b, Vb.shape[1]


@pytest.mark.parametrize("n_max", [1, 2])
def test_route_b_band_ranks_dim_O_and_kernel(n_max):
    """Route B (sampled Wigner-D, no CG): band ranks N^2 for N <= 2 n_max (top
    band included), 0 above; dim O = 5, 30; kernel of [D, .] on O = C.1."""
    ranks, allv, Dc, dh = _route_b(n_max)
    for N, r in ranks.items():
        expect = N * N if N <= 2 * n_max else 0
        assert r == expect, f"Route B n_max={n_max}, band N={N}: rank {r}, expected {expect}"
    dim_O = int(np.linalg.matrix_rank(allv, tol=1e-7))
    assert dim_O == DIM_O_PRINTED[n_max], f"Route B n_max={n_max}: dim O {dim_O}"
    _, _, Vh = np.linalg.svd(allv, full_matrices=False)
    B = Vh[:dim_O]
    comm = np.array([(Dc @ b.reshape(dh, dh) - b.reshape(dh, dh) @ Dc).reshape(-1) for b in B])
    rc = int(np.linalg.matrix_rank(comm, tol=1e-7))
    assert dim_O - rc == 1, f"Route B n_max={n_max}: kernel dim {dim_O - rc}"
    one = np.eye(dh).reshape(-1)
    assert int(np.linalg.matrix_rank(np.vstack([B, one]), tol=1e-7)) == dim_O, "1 not in O (B)"


@pytest.mark.parametrize("n_max", [1, 2])
def test_routes_a_and_b_agree(n_max):
    """The two independent constructions agree: same window dimension, same
    rank in every band, same dim O."""
    ranks_a, vecs_a, _, dimH_a, _ = _window(n_max)
    ranks_b, allv_b, _, dimH_b = _route_b(n_max)
    assert dimH_a == dimH_b, f"n_max={n_max}: window dims {dimH_a} vs {dimH_b}"
    assert ranks_a == ranks_b, f"n_max={n_max}: band ranks A {ranks_a} vs B {ranks_b}"
    assert (int(np.linalg.matrix_rank(vecs_a, tol=1e-9))
            == int(np.linalg.matrix_rank(allv_b, tol=1e-7)))
