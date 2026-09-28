"""Sign problem as Z2 holonomy -- re-measured on He's PHYSICAL Hamiltonian (2026-09-25).

Physical He Hamiltonian in the k=Z hydrogenic basis: h1 = diag(-Z^2/2n^2) exactly
(no (1/16)A graph hop), V_ee = exact Slater integrals (geovac.casimir_ci
.two_electron_integral).  The only off-diagonals are V_ee matrix elements.

Measures (per n_max = 2,3,4, singlet):
  - fraction of positive off-diagonals
  - cycle rank beta = E - V + C of the off-diagonal graph
  - frustrated fundamental cycles after a BFS spanning-tree +/-1 gauge
    (= non-tree edges left positive after the gauge; 0 <=> balanced <=> curable)
  - stoquastic gap Delta = E0(H) - E0(H_st), H_st = diag(H) - |offdiag(H)|
Both on the ground-state connected component (GS) and on the full space.
Native basis vs same-dimension random orthogonal orbital rotation (Gaunt lost),
radial-only rotation (Gaunt kept), random-sign control on the native graph.
Also: one-electron Fock adjacency bipartiteness; composed LiH re-verification.
"""
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"

import sys
import time
import warnings
from collections import deque

import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import eigsh

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
warnings.simplefilter("ignore")

TOL = 1e-12
Z = 2


# ------------------------------------------------------------------ tensors
def he_tensors(nmax, hop=False):
    from geovac.lattice import GeometricLattice
    from geovac.casimir_ci import two_electron_integral
    lat = GeometricLattice(max_n=nmax)
    orbs = list(lat.states)
    ns = len(orbs)
    h = np.diag([-Z * Z / (2.0 * n * n) for (n, l, m) in orbs])
    if hop:
        h = h + lat.adjacency.toarray() / 16.0
    V = np.zeros((ns, ns, ns, ns))
    for a in range(ns):
        for b in range(ns):
            for c in range(ns):
                for d in range(ns):
                    if orbs[a][2] + orbs[b][2] != orbs[c][2] + orbs[d][2]:
                        continue
                    V[a, b, c, d] = two_electron_integral(*orbs[a], *orbs[b],
                                                          *orbs[c], *orbs[d],
                                                          float(Z))
    return h, V, orbs


def singlet_H(h, V, keep=None):
    """Singlet pair basis (a<=b), physics-convention <ab|cd>."""
    ns = h.shape[0]
    pairs = [(a, b) for a in range(ns) for b in range(a, ns)]
    if keep is not None:
        pairs = [p for p in pairs if keep(*p)]
    I = np.eye(ns)
    Hp = (np.kron(h, I) + np.kron(I, h) + V.reshape(ns * ns, ns * ns))
    P = np.zeros((ns * ns, len(pairs)))
    for k, (a, b) in enumerate(pairs):
        P[a * ns + b, k] += 1
        P[b * ns + a, k] += 1
        P[:, k] /= np.linalg.norm(P[:, k])
    H = P.T @ Hp @ P
    return 0.5 * (H + H.T), pairs


def rotate(h, V, U):
    return U.T @ h @ U, np.einsum("abcd,ai,bj,ck,dl->ijkl", V, U, U, U, U,
                                  optimize=True)


def random_orth(n, rng):
    Q, R = np.linalg.qr(rng.standard_normal((n, n)))
    return Q * np.sign(np.diag(R))


# ------------------------------------------------------------------ signed-graph tools
def e0(H):
    H = np.asarray(H.toarray() if sp.issparse(H) else H)
    return float(np.linalg.eigvalsh(H)[0])


def edges_of(H):
    H = np.asarray(H)
    r, c = np.triu_indices(H.shape[0], 1)
    w = H[r, c]
    m = np.abs(w) > TOL
    return r[m], c[m], w[m]


def bfs_gauge(n, r, c, w):
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


def stoq(H):
    H = np.asarray(H)
    D = np.diag(np.diag(H))
    return D - np.abs(H - D)


def graph_stats(H, n_rand=5, rng=None):
    H = np.asarray(H)
    n = H.shape[0]
    r, c, w = edges_of(H)
    E = len(w)
    A = sp.csr_matrix((np.ones(E), (r, c)), shape=(n, n))
    ncomp, _ = connected_components(A + A.T, directed=False)
    beta = E - n + ncomp
    s = bfs_gauge(n, r, c, w)
    frust = int(np.sum(s[r] * s[c] * w > 0))
    E0 = e0(H)
    Est = e0(stoq(H))
    out = dict(dim=n, edges=E, comps=ncomp, beta=beta,
               posfrac=float(np.mean(w > 0)) if E else 0.0,
               frust=frust, E0=E0, delta=E0 - Est)
    if rng is not None and n_rand:
        d = np.diag(np.diag(H))
        dr = []
        for _ in range(n_rand):
            sg = rng.choice([-1.0, 1.0], size=E)
            M = np.zeros_like(H)
            M[r, c] = np.abs(w) * sg
            dr.append(e0(d + M + M.T) - Est)
        out["delta_rand"] = float(np.mean(dr))
        out["delta_rand_all"] = dr
    return out


def gs_component(H):
    H = np.asarray(H)
    A = sp.csr_matrix(np.abs(H) > TOL)
    _, lab = connected_components(A, directed=False)
    w, v = np.linalg.eigh(H)
    k = lab[np.argmax(np.abs(v[:, 0]))]
    idx = np.where(lab == k)[0]
    return H[np.ix_(idx, idx)], idx


def report(tag, H, rng, n_rand=5):
    Hg, _ = gs_component(H)
    g = graph_stats(Hg, n_rand=n_rand, rng=rng)
    f = graph_stats(H, n_rand=0)
    print(f"[{tag}] FULL dim={f['dim']} comps={f['comps']} edges={f['edges']} "
          f"beta={f['beta']} frust={f['frust']} Delta={f['delta']:.4e}")
    rs = f"  Delta_randsign(avg{n_rand})={g['delta_rand']:.4e}" if n_rand else ""
    print(f"      GS   dim={g['dim']} edges={g['edges']} beta={g['beta']} "
          f"posfrac={g['posfrac']:.3f} frust={g['frust']} E0={g['E0']:.8f} "
          f"Delta={g['delta']:.4e}{rs}")
    sys.stdout.flush()
    return f, g


# ------------------------------------------------------------------ sections
def one_body():
    from geovac.lattice import GeometricLattice
    print("\n=== one-electron Fock graph ===")
    for nmax in (3, 5):
        lat = GeometricLattice(max_n=nmax)
        A = lat.adjacency.toarray()
        par = np.array([(n + m) % 2 for (n, l, m) in lat.states])
        rr, cc = np.nonzero(A)
        bip = bool(np.all(par[rr] != par[cc]))
        G = np.diag([(-1.0) ** (n + m) for (n, l, m) in lat.states])
        gA = G @ A @ G
        maxoff = float(np.max(gA - np.diag(np.diag(gA))))
        dA = e0(A) - e0(stoq(A))
        hphys = np.diag([-Z * Z / (2.0 * n * n) for (n, l, m) in lat.states])
        print(f"n_max={nmax}: nodes={A.shape[0]} bipartite(n+m)={bip} "
              f"max offdiag of G A G={maxoff:+.1f} (<=0 => stoquastic) "
              f"Delta(bare A)={dA:.1e}; physical h1 offdiag nnz="
              f"{int(np.count_nonzero(hphys - np.diag(np.diag(hphys))))}")


def he(nmax, n_seeds=5):
    print(f"\n=== He singlet, n_max={nmax} ===")
    t0 = time.time()
    h, V, orbs = he_tensors(nmax)
    ns = len(orbs)
    mv = np.array([o[2] for o in orbs])
    # cross-check vs the production builder with h1 made diagonal (physical)
    import geovac.casimir_ci as C
    orig = C._build_graph_h1

    def diag_only(Zz, nm):
        h1, ob = orig(Zz, nm)
        return np.diag(np.diag(h1)), ob
    C._build_graph_h1 = diag_only
    try:
        Hb = C.build_graph_native_fci(Z=Z, n_max=nmax)
    finally:
        C._build_graph_h1 = orig
    H0, _ = singlet_H(h, V, keep=lambda a, b: mv[a] + mv[b] == 0)
    print(f"tensors {time.time()-t0:.1f}s; ns={ns}; M_L=0 dim={H0.shape[0]}; "
          f"max|mine - build_graph_native_fci(diag h1)|={np.max(np.abs(H0-Hb)):.1e}; "
          f"h offdiag nnz={int(np.count_nonzero(h-np.diag(np.diag(h))))}")
    Hhop = C.build_graph_native_fci(Z=Z, n_max=nmax)
    rng = np.random.default_rng(20260925 + nmax)
    res = {}
    res["phys_ML0"] = report(f"PHYSICAL M_L=0 n_max={nmax}", H0, rng)
    res["hop_ML0"] = report(f"HOP (retracted) M_L=0 n_max={nmax}", Hhop, rng)
    Hn, _ = singlet_H(h, V)
    res["native"] = report(f"PHYSICAL native all-M_L n_max={nmax}", Hn, rng)
    # radial rotation (Gaunt kept)
    U = np.zeros((ns, ns))
    groups = {}
    for i, (n, l, m) in enumerate(orbs):
        groups.setdefault((l, m), []).append(i)
    ru = {}
    for (l, m), ids in groups.items():
        if l not in ru:
            ru[l] = random_orth(len(ids), rng)
        U[np.ix_(ids, ids)] = ru[l]
    hr, Vr = rotate(h, V, U)
    Hr, _ = singlet_H(hr, Vr)
    res["radial"] = report(f"radial-rotated (Gaunt kept) n_max={nmax}", Hr, rng, n_rand=0)
    rot = []
    for sd in range(n_seeds):
        rs = np.random.default_rng(1000 + sd)
        hf, Vf = rotate(h, V, random_orth(ns, rs))
        Hf, _ = singlet_H(hf, Vf)
        rot.append(report(f"fully-rotated seed={1000+sd} n_max={nmax}", Hf, rng, n_rand=0))
    fn, gn = res["native"]
    bf = np.array([r[0]["beta"] for r in rot])
    df = np.array([r[0]["delta"] for r in rot])
    print(f"SUMMARY n_max={nmax}: full-space beta native={fn['beta']} rotated="
          f"{bf.min()}..{bf.max()} ratio={bf.mean()/max(fn['beta'],1):.1f}x; "
          f"full-space Delta native={fn['delta']:.4e} rotated={df.min():.4e}..{df.max():.4e} "
          f"ratio={df.mean()/fn['delta']:.1f}x")
    g0 = res["phys_ML0"][1]
    print(f"SUMMARY n_max={nmax}: GS-sector Delta native={g0['delta']:.4e} vs random-sign "
          f"{g0['delta_rand']:.4e} (seeds {np.round(g0['delta_rand_all'],5).tolist()})")
    return res


def lih():
    print("\n=== composed LiH re-verification (build_composed_hamiltonian, continuum integrals) ===")
    from geovac.molecular_spec import lih_spec
    from geovac.composed_qubit import build_composed_hamiltonian
    import openfermion as of
    out = build_composed_hamiltonian(lih_spec())
    Q = out["Q"]
    H = of.get_number_preserving_sparse_operator(out["fermion_op"], Q, 4,
                                                 spin_preserving=False)
    H = sp.csr_matrix(H.real)
    ncomp, lab = connected_components(abs(H) > 0, directed=False)
    w, v = eigsh(H, k=1, which="SA", tol=1e-10)
    k = lab[np.argmax(np.abs(v[:, 0]))]
    idx = np.where(lab == k)[0]
    Hg = H[idx][:, idx].toarray()
    g = graph_stats(Hg, n_rand=3, rng=np.random.default_rng(7))
    print(f"Q={Q} N_pauli={out['N_pauli']} N=4 dim={H.shape[0]} nnz={H.nnz} components={ncomp}; "
          f"GS comp dim={g['dim']} edges={g['edges']} beta={g['beta']} frust={g['frust']} "
          f"E0={g['E0']:.8f} (global E0 {w[0]:.8f}) Delta={g['delta']:.4e} "
          f"Delta_rand={g['delta_rand']:.4e}")


if __name__ == "__main__":
    T0 = time.time()
    one_body()
    for nm in (2, 3, 4):
        he(nm)
    try:
        lih()
    except Exception as exc:
        print(f"LiH FAILED: {type(exc).__name__}: {exc}")
    print(f"\nTOTAL {time.time()-T0:.1f}s")
