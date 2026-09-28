"""Paper 22, Sec. "Selection-rule sparsity and the Monte Carlo sign problem".

Backs the [MEASURED] claims of that subsection on He's PHYSICAL Hamiltonian in
the k=Z hydrogenic basis (h1 = diag(-Z^2/2n^2) exactly -- NO (1/16)A graph
hop -- plus exact Slater V_ee), singlet, M_L = 0:

  * the Fock-graph adjacency is bipartite under (n+m) parity, so the gauge
    (-1)^(n+m) makes it stoquastic (Delta = 0);
  * frozen cycle rank / frustrated fundamental cycles / stoquastic gap Delta
    of the ground-state connected sector at n_max = 2, 3 (Table tab:sign);
  * destroying the Gaunt zeros by a random orbital rotation raises the cycle
    rank and Delta by > 10x (n_max = 3, fixed seed);
  * the native signs are NOT better than random signs on the same graph.

Wrong answers these reject (see debug/firetest_paper22_sign_holonomy.py):
no gauge fix, one V_ee sign flipped, the retracted hop Hamiltonian swapped
in, an identity "rotation".  Driver + log: debug/sign_holonomy_probe.py,
debug/data/sign_holonomy.log.
"""
from collections import deque

import numpy as np
import pytest

Z = 2
TOL = 1e-12


# ----------------------------------------------------------------- builders
def _orbitals(nmax):
    from geovac.lattice import GeometricLattice
    lat = GeometricLattice(max_n=nmax)
    return lat, list(lat.states)


def _he_tensors(nmax):
    """Physical one-body (diagonal) and exact V_ee tensor <ab|cd>."""
    from geovac.casimir_ci import two_electron_integral
    lat, orbs = _orbitals(nmax)
    ns = len(orbs)
    h = np.diag([-Z * Z / (2.0 * n * n) for (n, l, m) in orbs])
    V = np.zeros((ns, ns, ns, ns))
    for a in range(ns):
        for b in range(ns):
            for c in range(ns):
                for d in range(ns):
                    if orbs[a][2] + orbs[b][2] != orbs[c][2] + orbs[d][2]:
                        continue
                    V[a, b, c, d] = two_electron_integral(
                        *orbs[a], *orbs[b], *orbs[c], *orbs[d], float(Z))
    return h, V, orbs


def _singlet(h, V, orbs, ml0=True):
    ns = h.shape[0]
    pairs = [(a, b) for a in range(ns) for b in range(a, ns)
             if (not ml0) or orbs[a][2] + orbs[b][2] == 0]
    I = np.eye(ns)
    Hp = np.kron(h, I) + np.kron(I, h) + V.reshape(ns * ns, ns * ns)
    P = np.zeros((ns * ns, len(pairs)))
    for k, (a, b) in enumerate(pairs):
        P[a * ns + b, k] += 1
        P[b * ns + a, k] += 1
        P[:, k] /= np.linalg.norm(P[:, k])
    H = P.T @ Hp @ P
    H = 0.5 * (H + H.T)
    return H


def _rotate(h, V, U):
    return U.T @ h @ U, np.einsum("abcd,ai,bj,ck,dl->ijkl", V, U, U, U, U,
                                  optimize=True)


def _random_orth(n, rng):
    Q, R = np.linalg.qr(rng.standard_normal((n, n)))
    return Q * np.sign(np.diag(R))


# ----------------------------------------------------------------- signed graph
def _edges(H):
    r, c = np.triu_indices(H.shape[0], 1)
    w = H[r, c]
    m = np.abs(w) > TOL
    return r[m], c[m], w[m]


def _bfs_gauge(n, r, c, w):
    adj = [[] for _ in range(n)]
    for a, b, x in zip(r, c, w):
        adj[a].append((b, x))
        adj[b].append((a, x))
    s = np.zeros(n)
    for root in range(n):
        if s[root] != 0:
            continue
        s[root] = 1
        q = deque([root])
        while q:
            i = q.popleft()
            for j, x in adj[i]:
                if s[j] == 0:
                    s[j] = -np.sign(x) * s[i]
                    q.append(j)
    return s


def _n_components(n, r, c):
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for a, b in zip(r, c):
        parent[find(a)] = find(b)
    return len({find(i) for i in range(n)})


def _stoq(H):
    D = np.diag(np.diag(H))
    return D - np.abs(H - D)


def _e0(H):
    return float(np.linalg.eigvalsh(H)[0])


def _gs_sector(H):
    n = H.shape[0]
    r, c, _ = _edges(H)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for a, b in zip(r, c):
        parent[find(a)] = find(b)
    v = np.linalg.eigh(H)[1][:, 0]
    root = find(int(np.argmax(np.abs(v))))
    idx = [i for i in range(n) if find(i) == root]
    return H[np.ix_(idx, idx)]


def _stats(H):
    n = H.shape[0]
    r, c, w = _edges(H)
    beta = len(w) - n + _n_components(n, r, c)
    s = _bfs_gauge(n, r, c, w)
    frust = int(np.sum(s[r] * s[c] * w > 0))
    return dict(dim=n, edges=len(w), npos=int(np.sum(w > 0)), beta=beta,
                frust=frust, E0=_e0(H), delta=_e0(H) - _e0(_stoq(H)))


_CACHE = {}


def _phys(nmax):
    if nmax not in _CACHE:
        h, V, orbs = _he_tensors(nmax)
        _CACHE[nmax] = (h, V, orbs, _singlet(h, V, orbs))
    return _CACHE[nmax]


# ----------------------------------------------------------------- tests
@pytest.mark.parametrize("nmax", [3, 4])
def test_fock_adjacency_gauge_free(nmax):
    """(-1)^(n+m) maps A to -A: bipartite, stoquastic after gauge, Delta = 0.
    Control: the ungauged A has positive off-diagonals (the gauge does work)."""
    lat, orbs = _orbitals(nmax)
    A = lat.adjacency.toarray()
    assert np.max(A - np.diag(np.diag(A))) > 0          # not stoquastic as is
    G = np.diag([(-1.0) ** (n + m) for (n, l, m) in orbs])
    gA = G @ A @ G
    assert np.max(gA - np.diag(np.diag(gA))) <= 0.0     # stoquastic after gauge
    assert np.allclose(gA, -A)                           # bipartite: GAG = -A
    assert abs(_e0(A) - _e0(_stoq(A))) < 1e-12           # Delta = 0


@pytest.mark.parametrize("nmax", [2, 3])
def test_physical_hamiltonian_is_builder_with_diagonal_h1(nmax, monkeypatch):
    """Our H equals build_graph_native_fci with h1 replaced by its diagonal,
    and the physical one-body operator has no off-diagonals."""
    import geovac.casimir_ci as C
    h, V, orbs, H = _phys(nmax)
    assert np.count_nonzero(h - np.diag(np.diag(h))) == 0
    orig = C._build_graph_h1

    def diag_only(Zz, nm):
        h1, ob = orig(Zz, nm)
        return np.diag(np.diag(h1)), ob
    monkeypatch.setattr(C, "_build_graph_h1", diag_only)
    Hb = C.build_graph_native_fci(Z=Z, n_max=nmax)
    assert Hb.shape == H.shape
    assert np.max(np.abs(Hb - H)) < 1e-12


# frozen from debug/data/sign_holonomy.log (2026-09-25)
FROZEN = {
    2: dict(dim=5, edges=10, npos=6, beta=6, frust=4,
            E0=-2.8334051759, delta=9.131622926e-04),
    3: dict(dim=19, edges=153, npos=85, beta=135, frust=82,
            E0=-2.8429083710, delta=2.290287137e-02),
}


@pytest.mark.parametrize("nmax", [2, 3])
def test_frozen_gs_sector_frustration(nmax):
    g = _stats(_gs_sector(_phys(nmax)[3]))
    f = FROZEN[nmax]
    for k in ("dim", "edges", "npos", "beta", "frust"):
        assert g[k] == f[k], (k, g[k], f[k])
    assert g["E0"] == pytest.approx(f["E0"], abs=1e-9)
    assert g["delta"] == pytest.approx(f["delta"], rel=1e-7)
    assert g["frust"] > 0 and g["delta"] > 1e-4        # genuinely frustrated


def test_gaunt_zeros_contain_sign_problem_n3():
    """Same-dimension random orbital rotation (Gaunt zeros destroyed, E0
    unchanged) raises cycle rank and Delta by > 10x in the full singlet space."""
    h, V, orbs, _ = _phys(3)
    Hn = _singlet(h, V, orbs, ml0=False)
    U = _random_orth(len(orbs), np.random.default_rng(1000))
    hr, Vr = _rotate(h, V, U)
    Hr = _singlet(hr, Vr, orbs, ml0=False)
    sn, sr = _stats(Hn), _stats(Hr)
    assert sn["dim"] == sr["dim"] == 105
    assert sr["E0"] == pytest.approx(sn["E0"], abs=1e-10)   # pure basis change
    assert sn["beta"] == 404 and sr["beta"] == 5356
    assert sr["beta"] > 10 * sn["beta"]
    assert sr["delta"] > 10 * sn["delta"]


def test_native_signs_not_better_than_random_n3():
    """Native Delta exceeds the mean Delta of random signs on the same graph."""
    H = _gs_sector(_phys(3)[3])
    r, c, w = _edges(H)
    D = np.diag(np.diag(H))
    e_st = _e0(_stoq(H))
    rng = np.random.default_rng(20260925)
    dr = []
    for _ in range(5):
        M = np.zeros_like(H)
        M[r, c] = np.abs(w) * rng.choice([-1.0, 1.0], size=len(w))
        dr.append(_e0(D + M + M.T) - e_st)
    native = _e0(H) - e_st
    assert native > 1.5 * float(np.mean(dr))
