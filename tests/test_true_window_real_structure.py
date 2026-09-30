"""Real structure J and order-zero / order-one on the TRUE Camporesi-Higuchi window.

Backs Paper 32 (axiom audit: item 5, "Net status of the real structure",
Table tab:axiom_audit) and Paper 38 (Sec. "Substrate: the true window"):
the first measurement of the real-structure axioms on the object the
convergence theorem is about -- the complete Dirac eigenspaces
|lambda| <= n_max + 1/2 of the round-S^3 spinor Dirac, O = P C(S^3) P --
rather than on the doubled Fock / Weyl MODEL objects.

Realization (tests/test_p38_true_ch_window.py): L^2(SU(2)) (x) C^2, basis
sqrt(2j+1) D^j_{m m'} (x) e_s, D = 2 sigma.J_R + 3/2.
Closed-form J = U K,  U e_{j,m,m',s} = (-1)^{m-m'} (i sigma_y)_{s's} e_{j,-m,-m',s'}
(pointwise complex conjugation of the spinor, twisted by i sigma_y).

Measured facts pinned here (driver: debug/true_window_real_structure.py):
  * sign triple (epsilon, epsilon') = (-1, +1) on the window at n_max = 1..3;
  * the Dirac + isometry-equivariant antiunitaries form a (2 n_max)-dim family
    of block phases, every member has J^2 = -1, and NONE anticommutes with D;
  * J O J^{-1} = O exactly, and J b^* J^{-1} = b (the compressed functions);
  * order zero FAILS already at n_max = 1 (12/25 pairs), residual grows:
    max 1/2, 5/9, 11/16 (n_max = 1, 2, 3); failing fraction 48% -> 64% -> 71%;
  * order one FAILS: max-abs 3/4, 5/3, 21/8; op-norm 0.75, 1.863, 2.902;
  * degree-1 generators: order-zero op-norm 1/2, 5/9, 11/16, 19/25 (n_max=1..4),
    order-one 3/4, 7/9, 13/16, 21/25 -- increasing toward the sup-norm scale 1;
  * mechanism: every residual is a Toeplitz boundary term -- it vanishes on all
    Dirac levels n <= top - (max(N_a, N_b) - 1);
  * an independent sampled construction (Route B: no Clebsch-Gordan, J from
    pointwise conjugation of sampled functions) reproduces every per-pair norm
    to 1e-10, and independently gives commutant = scalars with alg(O) = B(H)
    (dim 16, 256 at n_max = 1, 2 by word closure);
  * a third construction (Route C: polynomials on R^4, exact sphere moments,
    no Wigner-D / Clebsch-Gordan / sampling) gives the same commutant, grading
    space, J^2 = -1, JD = +DJ and J O J^{-1} = O at n_max = 1, 2;
  * J relations (J^2 = -1, JD = +DJ, window and O preserved) pinned at
    n_max = 1..4 (dim O = 5, 30, 91, 204).
Not pinned here: alg(O) = B(H) at n_max = 3 (dim 1600) -- the word closure at
d = 40 is memory-heavy; n_max = 3 oddness rests on Route A's commutant.
"""
from __future__ import annotations

import itertools
from fractions import Fraction
from functools import lru_cache

import numpy as np
import pytest

from test_p38_true_ch_window import SX, SY, SZ, _build, _mult_op, _spin_mats

ISY = np.array([[0, 1], [-1, 0]], complex)
TOL = 1e-10


def _U(labels, idx):
    U = np.zeros((len(labels), len(labels)), complex)
    for (j, m, mp, s) in labels:
        for s2 in (0, 1):
            if ISY[s2, s] != 0:
                U[idx[(j, -m, -mp, s2)], idx[(j, m, mp, s)]] = (-1) ** int(m - mp) * ISY[s2, s]
    return U


@lru_cache(maxsize=None)
def _setup(n_max: int, bands: int | None = None):
    labels, idx, D, ev, Ue = _build(n_max)
    keep = np.abs(ev) <= n_max + 0.5 + 1e-9
    V = Ue[:, keep]
    lev = np.rint(np.abs(ev[keep]) - 1.5).astype(int)
    Dc = V.conj().T @ D @ V
    U = _U(labels, idx)
    Jw = V.conj().T @ U @ V.conj()
    labs, amb = [], []
    Nrange = range(1, 2 * n_max + 1) if bands is None else [bands]
    for N in Nrange:
        J = Fraction(N - 1, 2)
        Ms = [J - k for k in range(N)]
        for M, Mp in itertools.product(Ms, Ms):
            labs.append((N, M, Mp))
            amb.append(_mult_op(J, M, Mp, labels, idx))
    amb = np.array(amb)
    Ts = V.conj().T @ amb @ V
    Bs = Jw @ np.transpose(Ts, (0, 2, 1)) @ Jw.conj().T  # J b^* J^{-1}
    return dict(labels=labels, idx=idx, D=D, V=V, lev=lev, Dc=Dc, U=U, Jw=Jw,
                labs=labs, amb=amb, Ts=Ts, Bs=Bs)


def _stats(n_max: int, order: int):
    S = _setup(n_max)
    Ts, Bs, V, Dc = S["Ts"], S["Bs"], S["V"], S["Dc"]
    Xs = Ts if order == 0 else Dc @ Ts - Ts @ Dc
    n = len(Ts)
    maxabs = np.zeros((n, n))
    opn = np.zeros((n, n))
    for i in range(n):
        C = Xs[i] @ Bs - Bs @ Xs[i]
        maxabs[i] = np.abs(V @ C @ V.conj().T).reshape(n, -1).max(axis=1)
        opn[i] = np.linalg.norm(C, ord=2, axis=(1, 2))
    return maxabs, opn


# ---------------------------------------------------------------- sign triple
@pytest.mark.parametrize("n_max", [1, 2, 3, 4])
def test_sign_triple_is_KO3(n_max):
    """(epsilon, epsilon') = (-1, +1) on the true window.
    Rejects: J^2 = +1 (a real, not quaternionic, structure), J D = -D J
    (KO-1/5/7-type sign), and a J that leaks out of the window."""
    S = _setup(n_max)
    U, D, Jw, Dc, V = S["U"], S["D"], S["Jw"], S["Dc"], S["V"]
    I = np.eye(len(U))
    assert np.abs(U @ U.conj() + I).max() < 1e-13           # J^2 = -1
    assert np.abs(U @ U.conj() - I).max() > 1.0              # not +1
    assert np.abs(U @ D.conj() - D @ U).max() < 1e-13        # JD = +DJ
    assert np.abs(U @ D.conj() + D @ U).max() > 1.0          # not -DJ
    W = U @ V.conj()
    assert np.abs(V @ (V.conj().T @ W) - W).max() < 1e-12    # window preserved
    d = Jw.shape[0]
    assert np.abs(Jw @ Jw.conj().T - np.eye(d)).max() < 1e-12
    assert np.abs(Jw @ Jw.conj() + np.eye(d)).max() < 1e-12
    assert np.abs(Jw @ Dc.conj() - Dc @ Jw).max() < 1e-12


@pytest.mark.parametrize("n_max", [1, 2])
def test_equivariant_J_family_forces_signs(n_max):
    """Solving U conj(X) = s X U for X = D (s = +-1) and the isometry generators
    J_L^a, J_R^a + sigma^a/2 (s = -1): the epsilon' = +1 family is one phase per
    irreducible block (dim 2 n_max), contains the closed form, every member has
    J^2 = -1; the epsilon' = -1 family is EMPTY.
    Rejects: 'the sign epsilon is a phase convention' and 'an equivariant J
    anticommuting with D exists'."""
    S = _setup(n_max)
    labels, idx, V, Dc = S["labels"], S["idx"], S["V"], S["Dc"]
    d = V.shape[1]
    n = len(labels)
    gens = []
    for a in range(3):
        JL = np.zeros((n, n), complex)
        KR = np.zeros((n, n), complex)
        for j in sorted({l[0] for l in labels}):
            dj = int(2 * j) + 1
            Sm = _spin_mats(j)
            ms = [j - k for k in range(dj)]
            ids = [idx[(j, m, mp, s)] for m in ms for mp in ms for s in (0, 1)]
            JL[np.ix_(ids, ids)] = np.kron(np.kron(Sm[a], np.eye(dj)), np.eye(2))
            KR[np.ix_(ids, ids)] = (np.kron(np.kron(np.eye(dj), Sm[a]), np.eye(2))
                                    + np.kron(np.eye(dj * dj), (SX, SY, SZ)[a]) / 2)
        gens += [V.conj().T @ JL @ V, V.conj().T @ KR @ V]
    I = np.eye(d)
    dims = {}
    for sD in (1, -1):
        M = np.vstack([np.kron(X.conj().T, I) - s * np.kron(I, X)
                       for X, s in [(Dc, sD)] + [(G, -1) for G in gens]])
        _, sv, Vh = np.linalg.svd(M, full_matrices=False)
        null = int((sv < 1e-9 * sv.max()).sum())
        dims[sD] = null
        if sD == 1:
            sols = Vh[-null:].conj()
            v = S["Jw"].reshape(-1, order="F")
            assert np.abs(sols.T @ (sols.conj() @ v) - v).max() < 1e-10
            rng = np.random.default_rng(3)
            Ug = ((rng.standard_normal(null) + 1j * rng.standard_normal(null)) @ sols)
            u, _, vh = np.linalg.svd(Ug.reshape(d, d, order="F"))
            Uu = u @ vh
            assert np.abs(Uu @ Uu.conj() + I).max() < 1e-10
    assert dims == {1: 2 * n_max, -1: 0}, dims


def test_ambient_conjugation_of_multipliers():
    """U conj(A_{J M M'}) U^dag = (-1)^{M-M'} A_{J,-M,-M'} (conj D^J_{MM'} =
    (-1)^{M-M'} D^J_{-M,-M'}), n_max = 2, exactly.  Rejects the phase-less
    map e_{m m'} -> e_{-m,-m'} (which fails for M - M' odd)."""
    S = _setup(2)
    U, amb, labs = S["U"], S["amb"], S["labs"]
    li = {l: i for i, l in enumerate(labs)}
    odd_seen = False
    for (N, M, Mp), A in zip(labs, amb):
        B = (-1) ** int(M - Mp) * amb[li[(N, -M, -Mp)]]
        assert np.abs(U @ A.conj() @ U.conj().T - B).max() < 1e-13
        odd_seen |= int(M - Mp) % 2 == 1 and np.abs(B).max() > 0
    assert odd_seen


# ---------------------------------------------------------------- J O J^{-1}
@pytest.mark.parametrize("n_max", [1, 2, 3, 4])
def test_J_preserves_operator_system(n_max):
    """J O J^{-1} = O (every J T J^{-1} lies in span O), and J b^* J^{-1} = b.
    Rejects a J that maps a compressed function out of the operator system."""
    S = _setup(n_max)
    Ts, Jw = S["Ts"], S["Jw"]
    flat = Ts.reshape(len(Ts), -1)
    r = int(np.linalg.matrix_rank(flat, tol=1e-9))
    assert r == {1: 5, 2: 30, 3: 91, 4: 204}[n_max]
    _, _, Vh = np.linalg.svd(flat, full_matrices=False)
    Q = Vh[:r]
    JTJ = (Jw @ Ts.conj() @ Jw.conj().T).reshape(len(Ts), -1)
    assert np.abs(JTJ @ Q.conj().T @ Q - JTJ).max() < 1e-12
    assert np.abs(S["Bs"] - Ts).max() < 1e-12


# ---------------------------------------------------------------- order zero / one
ORDER0 = {1: (12, 25, 1 / 2), 2: (572, 900, 5 / 9), 3: (5842, 8281, 11 / 16)}
ORDER1 = {1: (4, 25, 3 / 4, 0.75), 2: (513, 900, 5 / 3, 1.8633899812),
          3: (5856, 8281, 21 / 8, 2.9020466916)}


def test_order_zero_fails_and_grows():
    """Order zero [a, J b^* J^{-1}] on the true window: failing pairs 12/25,
    572/900, 5842/8281, max residual (ambient basis = op norm) 1/2, 5/9, 11/16.
    Rejects: 'order zero holds on the true window' (0 failing), 'fails only
    from n_max = 2 like the model' (n_max = 1 fails here), and 'decays'."""
    fr, mx = [], []
    for n, (nf, nt, m) in ORDER0.items():
        ma, op = _stats(n, 0)
        assert ma.size == nt
        assert int((ma > TOL).sum()) == nf, (n, int((ma > TOL).sum()))
        assert abs(ma.max() - m) < 1e-12 and abs(op.max() - m) < 1e-12, (n, ma.max(), op.max())
        fr.append(nf / nt)
        mx.append(m)
    assert fr[0] < fr[1] < fr[2] and mx[0] < mx[1] < mx[2]


def test_order_one_fails_and_grows():
    """Order one [[D, a], J b^* J^{-1}]: failing 4/25, 513/900, 5856/8281;
    max-abs 3/4, 5/3, 21/8; op-norm 0.75, 1.8634, 2.9020.  Rejects 'order one
    holds on the true window' and 'decays'."""
    ops = []
    for n, (nf, nt, m, o) in ORDER1.items():
        ma, op = _stats(n, 1)
        assert int((ma > TOL).sum()) == nf, (n, int((ma > TOL).sum()))
        assert abs(ma.max() - m) < 1e-12, (n, ma.max())
        assert abs(op.max() - o) < 1e-8, (n, op.max())
        ops.append(op.max())
    assert ops[0] < ops[1] < ops[2]


FIXED0 = {1: 1 / 2, 2: 5 / 9, 3: 11 / 16, 4: 19 / 25}
FIXED1 = {1: 3 / 4, 2: 7 / 9, 3: 13 / 16, 4: 21 / 25}


@pytest.mark.parametrize("n_max", [1, 2, 3, 4])
def test_degree1_generator_residuals(n_max):
    """Fixed generators (the 4 degree-1 matrix coefficients, N = 2): order-zero
    op-norm 1/2, 5/9, 11/16, 19/25 and order-one 3/4, 7/9, 13/16, 21/25 at
    n_max = 1..4 -- increasing, not decaying.  Rejects the 'finite-resolution
    artifact that vanishes as n_max grows' reading on fixed generators."""
    S = _setup(n_max, 2)
    Ts, Bs, Dc = S["Ts"], S["Bs"], S["Dc"]
    w0 = max(np.linalg.norm(a @ b - b @ a, 2) for a in Ts for b in Bs)
    DT = [Dc @ a - a @ Dc for a in Ts]
    w1 = max(np.linalg.norm(x @ b - b @ x, 2) for x in DT for b in Bs)
    assert abs(w0 - FIXED0[n_max]) < 1e-12, w0
    assert abs(w1 - FIXED1[n_max]) < 1e-12, w1
    if n_max > 1:
        assert FIXED0[n_max] > FIXED0[n_max - 1] and FIXED1[n_max] > FIXED1[n_max - 1]


@pytest.mark.parametrize("n_max", [2, 3])
def test_residual_is_a_boundary_term(n_max):
    """Mechanism: [PfP, PgP] = P g (1-P) f P - P f (1-P) g P (and likewise at
    order one, [D, f] = c(df) commuting with g pointwise), and a band-N
    multiplier moves the Dirac level by <= N - 1, so every residual vanishes on
    all levels n <= top - (max(N_a, N_b) - 1).  The check is non-vacuous (the
    inner region is non-empty for some pair).  Rejects a residual living in the
    interior of the window (a bulk failure of commutativity)."""
    S = _setup(n_max)
    Ts, Bs, Dc, lev = S["Ts"], S["Bs"], S["Dc"], S["lev"]
    Ns = [l[0] for l in S["labs"]]
    top = n_max - 1
    checked = 0
    for order in (0, 1):
        Xs = Ts if order == 0 else Dc @ Ts - Ts @ Dc
        for i in range(len(Ts)):
            C = Xs[i] @ Bs - Bs @ Xs[i]
            for j in range(len(Ts)):
                inner = lev <= top - (max(Ns[i], Ns[j]) - 1)
                if inner.any():
                    checked += 1
                    assert np.abs(C[j][inner, :]).max() < 1e-12
                    assert np.abs(C[j][:, inner]).max() < 1e-12
    assert checked > 0


# ---------------------------------------------------------------- Route B
@lru_cache(maxsize=None)
def _route_b_build(n_max: int, n_samples: int = 900):
    """Sampled Wigner-D build: no Clebsch-Gordan code, and J's function part is
    the least-squares expansion of the pointwise-conjugated sampled basis (no
    (-1)^{m-m'} formula).  Returns the window Dirac Dc, J (Jw) and the
    compressed multipliers Ts, all in Route B's own eigenbasis."""
    from scipy.linalg import expm
    P = (SX, SY, SZ)

    def spin(j2):
        j = j2 / 2
        ms = [j - k for k in range(j2 + 1)]
        jp = np.zeros((j2 + 1, j2 + 1), complex)
        for a, m in enumerate(ms):
            if a > 0:
                jp[a - 1, a] = np.sqrt(j * (j + 1) - m * (m + 1))
        return [(jp + jp.conj().T) / 2, (jp - jp.conj().T) / (2j), np.diag(ms).astype(complex)]

    def wig(j2, g):
        w, Vv = np.linalg.eig(g)
        H = Vv @ np.diag(-np.angle(w)) @ np.linalg.inv(Vv)
        c = [np.trace(H @ s_).real for s_ in P]
        return expm(-1j * sum(ci * Ji for ci, Ji in zip(c, spin(j2))))

    rng = np.random.default_rng(5)
    gs = []
    for _ in range(n_samples):
        q = rng.standard_normal(4)
        q /= np.linalg.norm(q)
        gs.append(np.array([[q[0] + 1j * q[3], q[2] + 1j * q[1]],
                            [-q[2] + 1j * q[1], q[0] - 1j * q[3]]]))
    cols, labs, Dc_ = [], [], {}
    for j2 in range(3 * n_max + 2):
        Ds = np.array([wig(j2, g) for g in gs])
        Dc_[j2] = Ds
        for a in range(j2 + 1):
            for b in range(j2 + 1):
                cols.append(np.sqrt(j2 + 1) * Ds[:, a, b])
                labs.append((j2, a, b))
    F = np.array(cols).T
    Fp = np.linalg.pinv(F)
    H = [i for i, l in enumerate(labs) if l[0] <= n_max]
    Hl = [labs[i] for i in H]
    nf = len(Hl)
    full = [(l, s) for l in Hl for s in (0, 1)]
    pos = {f: i for i, f in enumerate(full)}
    Dop = np.zeros((2 * nf, 2 * nf), complex)
    for (j2, a, b), s1 in full:
        J = spin(j2)
        for (j2p, a2, b2), s2 in full:
            if j2p == j2 and a2 == a:
                v = sum(2 * J[k][b2, b] * P[k][s2, s1] for k in range(3))
                v += 1.5 if (b2 == b and s2 == s1) else 0
                Dop[pos[((j2p, a2, b2), s2)], pos[((j2, a, b), s1)]] = v
    ev, Ue = np.linalg.eigh(Dop)
    V = Ue[:, np.abs(ev) <= n_max + 0.5 + 1e-9]
    Dc = V.conj().T @ Dop @ V
    Uf = np.array([(Fp @ np.conj(F[:, i]))[H] for i in H]).T
    Jw = V.conj().T @ np.kron(Uf, ISY) @ V.conj()
    Ts = []
    for N in range(1, 2 * n_max + 1):
        for M, Mp in itertools.product(range(N), range(N)):
            fv = Dc_[N - 1][:, M, Mp]
            A = np.array([(Fp @ (fv * F[:, i]))[H] for i in H]).T
            Ts.append(V.conj().T @ np.kron(A, np.eye(2)) @ V)
    Ts = np.array(Ts)
    return dict(Dc=Dc, Jw=Jw, Ts=Ts)


def _route_b_norms(n_max: int):
    """Per-pair order-0/1 op norms from the Route B build, and the J^2 + 1 error."""
    R = _route_b_build(n_max)
    Dc, Jw, Ts = R["Dc"], R["Jw"], R["Ts"]
    Bs = Jw @ np.transpose(Ts, (0, 2, 1)) @ Jw.conj().T
    out = []
    for order in (0, 1):
        Xs = Ts if order == 0 else Dc @ Ts - Ts @ Dc
        out.append(np.array([np.linalg.norm(Xs[i] @ Bs - Bs @ Xs[i], ord=2, axis=(1, 2))
                             for i in range(len(Ts))]))
    return out, float(np.abs(Jw @ Jw.conj() + np.eye(len(Jw))).max())


@pytest.mark.parametrize("n_max", [1, 2])
def test_route_b_reproduces_per_pair_norms(n_max):
    """Independent sampled route: J^2 = -1, and every per-pair order-zero and
    order-one op norm equals Route A's to 1e-10 (same pair ordering).
    Rejects an agreement that holds only for the maxima, and a Route-A J that
    differs from pointwise spinor conjugation."""
    (o0b, o1b), j2err = _route_b_norms(n_max)
    assert j2err < 1e-9
    _, o0a = _stats(n_max, 0)
    _, o1a = _stats(n_max, 1)
    assert np.abs(o0a - o0b).max() < 1e-10
    assert np.abs(o1a - o1b).max() < 1e-10


def _commutant_and_grading_dims(Ts, Dc):
    """(dim of {X : [X, a] = 0 for all a in Ts}, dim of that space further
    restricted to X D = -D X).  Null spaces by SVD, relative cut 1e-9."""
    d = Dc.shape[0]
    I = np.eye(d)
    rows = [np.kron(a, I) - np.kron(I, a.T) for a in Ts]
    sv = np.linalg.svd(np.vstack(rows), compute_uv=False)
    svg = np.linalg.svd(np.vstack(rows + [np.kron(Dc, I) + np.kron(I, Dc.T)]),
                        compute_uv=False)
    return int((sv < 1e-9 * sv.max()).sum()), int((svg < 1e-9 * svg.max()).sum())


def _generated_algebra_rank(Ts):
    """Dimension of the unital algebra generated by O: span{1} + O, then
    repeatedly left-multiplied by an orthonormal basis of O until the span
    stops growing (words of every length).  Returns (dimension, word length
    at which it stabilised)."""
    d = Ts.shape[1]

    def orth(M):
        _, sv, Vh = np.linalg.svd(M, full_matrices=False)
        return Vh[: int((sv > 1e-9 * sv.max()).sum())]

    G = orth(Ts.reshape(len(Ts), -1)).reshape(-1, d, d)
    Q = orth(np.vstack([np.eye(d).reshape(1, -1), G.reshape(len(G), -1)]))
    for length in range(2, 4 * d):
        prods = np.einsum("aij,bjk->abik", G, Q.reshape(-1, d, d)).reshape(-1, d * d)
        Qn = orth(np.vstack([Q, prods]))
        if len(Qn) == len(Q):
            return len(Q), length - 1
        Q = Qn
    return len(Q), None


@pytest.mark.parametrize("n_max", [1, 2, 3])
def test_true_window_is_odd_no_grading(n_max):
    """Paper 32 prop:odd(c), Paper 38 substrate note: the true window is ODD.
    The commutant of O^CH (null space of X -> [X, a] over a basis of O) is the
    scalars, and no nonzero X both commutes with O and anticommutes with D.
    Rejects: 'the finite truncation is even' carried over from the
    chirality-doubled model (sigma_x (x) 1 would be a nonzero solution)."""
    S = _setup(n_max)
    comm, grad = _commutant_and_grading_dims(S["Ts"], S["Dc"])
    assert comm == 1          # commutant = C * 1
    assert grad == 0          # no grading


@pytest.mark.parametrize("n_max", [1, 2])
def test_route_b_commutant_and_algebra(n_max):
    """Paper 32 prop:odd(c), SECOND route: the sampled-function construction
    (Route B -- no Clebsch-Gordan code) gives commutant = scalars, grading
    space {0}, and the unital algebra it generates is all of B(H) (dimension
    d^2 = 16, 256 by explicit word closure; the span stops growing after
    words of length 4 at both n_max, asserted).  Rejects a commutant result
    that holds only for the Clebsch-Gordan build.

    NEGATIVE CONTROL -- ASSERTED (promoted 2026-09-29).  Zeroing Route B's
    blocks that mix positive and negative Dirac levels (the doubled-model
    defect; Dc is diagonal in Route B's eigenbasis, so the mask is the sign
    pattern of its diagonal) must NOT pass the same checks:  n_max = 2 gives
    commutant 2, grading 0, algebra dim 128 < d^2 = 256 (this row carries
    the control).  At n_max = 1 every multiplier mixes the sign sectors, so the
    masked system collapses to the scalars (algebra dim 1); only that is
    asserted -- the commutant/grading counts there are SVD-noise on a
    (near-)zero matrix and are NOT pinned (the earlier (5, 8, 1) pin froze an
    artifact: grading <= commutant in exact arithmetic).  This proves the
    commutant and algebra-rank checks can return something other than 1 / d^2
    on this build, so the true-window pass is not vacuous."""
    R = _route_b_build(n_max)
    Ts, Dc = R["Ts"], R["Dc"]
    d = Dc.shape[0]
    assert d == {1: 4, 2: 16}[n_max]
    comm, grad = _commutant_and_grading_dims(Ts, Dc)
    assert comm == 1 and grad == 0, (comm, grad)
    dim_alg, word_len = _generated_algebra_rank(Ts)
    assert dim_alg == d * d, dim_alg
    assert word_len == 4, word_len
    # negative control: the doubled model (sign-mixing blocks zeroed)
    assert np.abs(Dc - np.diag(np.diag(Dc))).max() < 1e-12
    sgn = np.sign(np.diag(Dc).real)
    Tz = Ts * (sgn[:, None] == sgn[None, :])
    comm_z, grad_z = _commutant_and_grading_dims(Tz, Dc)
    dim_z, _ = _generated_algebra_rank(Tz)
    if n_max == 1:
        assert dim_z == 1, dim_z
    else:
        assert (comm_z, grad_z, dim_z) == (2, 0, 128), (comm_z, grad_z, dim_z)


# ---------------------------------------------------------------- Route C
# Third, fully independent build (ported from the /qa DELTA #15 reviewer's
# independent construction): polynomials of degree <= n_max on S^3 in R^4 with
# exact sphere moments; spin fields from quaternion RIGHT multiplication, the
# i/2 normalisation fixed numerically by [J1, J2] = i J3; D = 2 sigma.J + 3/2;
# J = (i sigma_y) K with K = complex conjugation of coefficients in a REAL
# orthonormal basis.  No Wigner-D, no Clebsch-Gordan, no sampling.
_R4 = (np.array([[0, -1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, -1, 0]], float),
       np.array([[0, 0, -1, 0], [0, 0, 0, -1], [1, 0, 0, 0], [0, 1, 0, 0]], float),
       np.array([[0, 0, 0, -1], [0, 0, 1, 0], [0, -1, 0, 0], [1, 0, 0, 0]], float))


def _monos(deg: int):
    out = []
    for k in range(deg + 1):
        for c in itertools.combinations_with_replacement(range(4), k):
            e = [0] * 4
            for i in c:
                e[i] += 1
            out.append(tuple(e))
    return out


@lru_cache(maxsize=None)
def _sphere_moment(a: tuple) -> float:
    """Normalised <x^a> on the unit S^3 (Dirichlet-type Gamma formula)."""
    from math import exp, lgamma
    if any(x % 2 for x in a):
        return 0.0
    b = [(x + 1) / 2 for x in a]
    return exp(sum(lgamma(x) for x in b) - lgamma(sum(b)) - (4 * lgamma(0.5) - lgamma(2)))


def _madd(a, b):
    return tuple(x + y for x, y in zip(a, b))


@lru_cache(maxsize=None)
def _route_c_build(n_max: int):
    ms = _monos(n_max)
    mi = {m: i for i, m in enumerate(ms)}
    n = len(ms)
    G = np.array([[_sphere_moment(_madd(a, b)) for b in ms] for a in ms])
    w, Ug = np.linalg.eigh(G)
    keep = w > 1e-10 * w.max()
    B = Ug[:, keep] / np.sqrt(w[keep])                 # real orthonormal basis
    dimf = B.shape[1]
    assert dimf == sum((k + 1) ** 2 for k in range(n_max + 1)), dimf
    Xs = []
    for R in _R4:
        M = np.zeros((n, n))
        for j, m in enumerate(ms):
            for p in range(4):
                if m[p] == 0:
                    continue
                for q in range(4):
                    if R[p, q] == 0:
                        continue
                    e = list(m)
                    e[p] -= 1
                    e[q] += 1
                    M[mi[tuple(e)], j] += m[p] * R[p, q]
        Xs.append(B.T @ G @ M @ B)
    for c in (0.5j, -0.5j):
        Jg = [c * X for X in Xs]
        if np.abs(Jg[0] @ Jg[1] - Jg[1] @ Jg[0] - 1j * Jg[2]).max() < 1e-9:
            break
    else:
        raise AssertionError("no normalisation gives [J1, J2] = i J3")
    D = 2 * sum(np.kron(Ja, s_) for Ja, s_ in zip(Jg, (SX, SY, SZ))) + 1.5 * np.eye(2 * dimf)
    ev, Ue = np.linalg.eigh(D)
    V = Ue[:, np.abs(ev) <= n_max + 0.5 + 1e-9]
    Dc = V.conj().T @ D @ V
    Jw = V.conj().T @ np.kron(np.eye(dimf), ISY) @ V.conj()
    Tp = []
    for e in _monos(2 * n_max + 1):
        Mm = np.array([[_sphere_moment(_madd(_madd(a, b), e)) for b in ms] for a in ms])
        Tp.append(V.conj().T @ np.kron(B.T @ Mm @ B, np.eye(2)) @ V)
    Tp = np.array(Tp)
    return dict(Dc=Dc, Jw=Jw, Ts=Tp)


@pytest.mark.parametrize("n_max", [1, 2])
def test_route_c_polynomial_window(n_max):
    """Paper 32 prop:odd(c) and audit item 6, THIRD route (polynomials on R^4):
    window spectrum +-(n + 3/2) with multiplicity (n+1)(n+2) each, window dim
    4 / 16, dim O = 5 / 30, commutant = scalars, grading space {0}, J^2 = -1,
    J D = +D J (and NOT -D J), J O J^{-1} = O.  Rejects the commutant result
    being an artifact of Wigner-D / Clebsch-Gordan machinery shared by
    Routes A and B.  (Measured, not asserted here: at n_max = 2, zeroing this
    route's sign-mixing blocks gives commutant 2; at n_max = 1 it gives
    commutant 5 / algebra dim 1.)"""
    C = _route_c_build(n_max)
    Dc, Jw, Ts = C["Dc"], C["Jw"], C["Ts"]
    d = Dc.shape[0]
    assert d == {1: 4, 2: 16}[n_max]
    expect = []
    for k in range(n_max):
        mult = (k + 1) * (k + 2)
        expect += [k + 1.5] * mult + [-(k + 1.5)] * mult
    w = np.sort(np.linalg.eigvalsh(Dc))
    assert np.abs(w - np.sort(expect)).max() < 1e-9, w
    flat = Ts.reshape(len(Ts), -1)
    r = int(np.linalg.matrix_rank(flat, tol=1e-9))
    assert r == {1: 5, 2: 30}[n_max]
    comm, grad = _commutant_and_grading_dims(Ts, Dc)
    assert comm == 1 and grad == 0, (comm, grad)
    I = np.eye(d)
    assert np.abs(Jw @ Jw.conj() + I).max() < 1e-12                 # J^2 = -1
    assert np.abs(Jw @ Dc.conj() - Dc @ Jw).max() < 1e-12            # JD = +DJ
    assert np.abs(Jw @ Dc.conj() + Dc @ Jw).max() > 1.0              # not -DJ
    _, _, Vh = np.linalg.svd(flat, full_matrices=False)
    Q = Vh[:r]
    JTJ = (Jw @ Ts.conj() @ Jw.conj().T).reshape(len(Ts), -1)
    assert np.abs(JTJ @ Q.conj().T @ Q - JTJ).max() < 1e-10