"""MC-discrete scoping probe (2026-09-25). Scratch diagnostic, not production.

P1: where does the sign problem enter GeoVac's discrete Hamiltonians?
    Signed-graph (frustration) analysis of the off-diagonal graph:
      - raw positive-off-diagonal fraction (convention-dependent)
      - cycle rank beta = E - V + C
      - diagonal +/-1 gauge by BFS spanning tree -> frustrated non-tree edges
        (balanced  <=>  0  <=>  sign problem curable by a diagonal gauge)
      - greedy local search upper bound on the frustration index
      - stoquastic gap Delta = E0(H) - E0(H_st), H_st = diag(H) - |offdiag(H)|
        (gauge-invariant; = asymptotic decay rate of the average sign in a
        projector / FCIQMC walk without annihilation)
    Systems:
      (a) one-electron Fock graph, hybrid h1 (exact diag + kappa*(-A))
      (b) He 2e singlet, hybrid h1 + exact Slater V_ee (LatticeIndex slater_full)
          native basis / radial-rotated (Gaunt kept) / fully rotated (Gaunt lost)
          + random-sign control on the native graph
      (c) composed LiH (Q=30) fixed-N=4 sector of the JW qubit Hamiltonian
          (= determinant space)
P2: Fock-graph heat kernel closed form vs expm; positivity of the one-body
    reference propagator after the bipartite gauge.
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
from scipy.linalg import expm, eigh

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
warnings.simplefilter("ignore")

RNG = np.random.default_rng(20260925)


def e0(H):
    H = sp.csr_matrix(H)
    n = H.shape[0]
    if n <= 1500:
        return float(np.linalg.eigvalsh(H.toarray())[0])
    return float(eigsh(H, k=1, which="SA", tol=1e-10)[0][0])


def gs(H):
    H = sp.csr_matrix(H)
    n = H.shape[0]
    if n <= 1500:
        w, v = np.linalg.eigh(H.toarray())
        return w[0], v[:, 0]
    w, v = eigsh(H, k=1, which="SA", tol=1e-10)
    return w[0], v[:, 0]


def edges_of(H, tol):
    C = sp.triu(sp.csr_matrix(H), k=1).tocoo()
    m = np.abs(C.data) > tol
    return C.row[m], C.col[m], C.data[m]


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
    return s, adj


def n_frust(s, r, c, w):
    return int(np.sum(s[r] * s[c] * w > 0))


def greedy(s, adj, r, c, w, sweeps=50):
    s = s.copy()
    for _ in range(sweeps):
        improved = False
        for i in range(len(s)):
            bad = sum(1 for j, x in adj[i] if s[i] * s[j] * x > 0)
            if bad * 2 > len(adj[i]):
                s[i] = -s[i]
                improved = True
        if not improved:
            break
    return s


def gs_component(H):
    H = sp.csr_matrix(H)
    ncomp, lab = connected_components(abs(H) > 0, directed=False)
    E, v = gs(H)
    k = lab[np.argmax(np.abs(v))]
    idx = np.where(lab == k)[0]
    return H[idx][:, idx], ncomp


def analyze(name, H, tol=1e-12, n_rand=5, do_comp=True):
    H = sp.csr_matrix(H)
    H = 0.5 * (H + H.T)
    full_n = H.shape[0]
    ncomp_full = None
    if do_comp:
        H, ncomp_full = gs_component(H)
    n = H.shape[0]
    r, c, w = edges_of(H, tol)
    E = len(w)
    ncomp, _ = connected_components(abs(H) > tol, directed=False)
    beta = E - n + ncomp
    posfrac = float(np.mean(w > 0)) if E else 0.0
    s, adj = bfs_gauge(n, r, c, w)
    f_bfs = n_frust(s, r, c, w)
    s_g = greedy(s, adj, r, c, w)
    f_greedy = n_frust(s_g, r, c, w)
    E0, v0 = gs(H)
    sg = np.where(v0 >= 0, 1.0, -1.0)
    f_gs = n_frust(sg, r, c, w)
    wf_gs = float(np.sum(np.abs(w[sg[r] * sg[c] * w > 0])) / max(np.sum(np.abs(w)), 1e-300))
    d = H.diagonal()
    Hst = sp.diags(d) + sp.csr_matrix((-np.abs(w), (r, c)), shape=(n, n))
    Hst = Hst + sp.triu(Hst, 1).T
    Hst = sp.csr_matrix(Hst)
    # fix: we added the upper triangle twice to itself? rebuild cleanly
    U = sp.csr_matrix((-np.abs(w), (r, c)), shape=(n, n))
    Hst = sp.diags(d) + U + U.T
    E0st = e0(Hst)
    delta = E0 - E0st
    drand = []
    for _ in range(n_rand):
        sg_r = RNG.choice([-1.0, 1.0], size=E)
        Ur = sp.csr_matrix((np.abs(w) * sg_r, (r, c)), shape=(n, n))
        drand.append(e0(sp.diags(d) + Ur + Ur.T) - E0st)
    out = dict(name=name, full_dim=full_n, n_sectors=ncomp_full, dim=n, E=E,
               beta=beta, posfrac=posfrac, f_bfs=f_bfs, f_greedy=f_greedy,
               f_gs=f_gs, wfrac_gs=wf_gs, E0=E0, E0st=E0st, delta=delta,
               delta_rand=float(np.mean(drand)) if drand else float("nan"))
    print(f"[{name}] full_dim={full_n} sectors={ncomp_full} GS-sector dim={n} "
          f"edges={E} beta={beta} raw_pos_frac={posfrac:.3f}")
    print(f"    frustrated edges: BFS-gauge={f_bfs}  greedy<= {f_greedy}  "
          f"GS-sign-gauge={f_gs} (|H| weight frac {wf_gs:.3e})")
    print(f"    E0={E0:.8f}  E0(stoq)={E0st:.8f}  Delta={delta:.3e}  "
          f"Delta_randsign(avg{n_rand})={out['delta_rand']:.3e}")
    sys.stdout.flush()
    return out


# ---------------------------------------------------------------- (a)
def fock_one_body():
    from geovac.lattice import GeometricLattice
    print("\n=== (a) one-electron Fock graph ===")
    for nmax in (3, 5):
        lat = GeometricLattice(max_n=nmax, nuclear_charge=1)
        A = lat.adjacency
        # bipartite check with parity (n+m)
        par = np.array([(n + m) % 2 for (n, l, m) in lat.states])
        C = A.tocoo()
        bip = bool(np.all(par[C.row] != par[C.col]))
        # hybrid h1 = diag(-Z^2/2n^2) + kappa*(-A), kappa=-1/16 -> +A/16
        d = np.array([-0.5 / n ** 2 for (n, l, m) in lat.states])
        H = sp.diags(d) + (-1.0 / 16.0) * (-A)
        print(f"n_max={nmax}: nodes={A.shape[0]}  parity-(n+m) bipartite={bip}")
        analyze(f"Fock h1 hybrid n_max={nmax}", H, n_rand=3, do_comp=False)


# ---------------------------------------------------------------- (b)
def he_tensors(nmax):
    from geovac.lattice_index import LatticeIndex
    idx = LatticeIndex(n_electrons=2, max_n=nmax, nuclear_charge=2,
                       vee_method="slater_full", h1_method="hybrid")
    h = idx._H1_spatial.toarray()
    ns = h.shape[0]
    V = np.zeros((ns, ns, ns, ns))
    for (a, b, c, dd), x in idx._eri.items():
        V[a, b, c, dd] = x
    return idx, h, V, idx.lattice.states


def singlet_H(h, V):
    ns = h.shape[0]
    I = np.eye(ns)
    Hp = np.kron(h, I) + np.kron(I, h) + V.reshape(ns * ns, ns * ns)
    cols = []
    for a in range(ns):
        for b in range(a, ns):
            v = np.zeros(ns * ns)
            v[a * ns + b] += 1
            v[b * ns + a] += 1
            cols.append(v / np.linalg.norm(v))
    P = np.array(cols).T
    return P.T @ Hp @ P


def rotate(h, V, U):
    h2 = U.T @ h @ U
    V2 = np.einsum("abcd,ai,bj,ck,dl->ijkl", V, U, U, U, U, optimize=True)
    return h2, V2


def random_orth(n):
    Q, R = np.linalg.qr(RNG.standard_normal((n, n)))
    return Q * np.sign(np.diag(R))


def he_block(nmax):
    print(f"\n=== (b) He 2e singlet, n_max={nmax} ===")
    idx, h, V, states = he_tensors(nmax)
    ns = h.shape[0]
    print(f"spatial orbitals={ns}; max|V - V^T(pair)|="
          f"{np.max(np.abs(V - V.transpose(2, 3, 0, 1))):.1e}; "
          f"complex? {np.iscomplexobj(V)}")
    Hn = singlet_H(h, V)
    res = [analyze(f"He native n_max={nmax}", Hn)]
    # exact-h1 variant (h1 diagonal): isolate V_ee as sole off-diagonal source
    hd = np.diag(np.diag(h))
    res.append(analyze(f"He exact-h1 (V_ee only) n_max={nmax}", singlet_H(hd, V)))
    # radial rotation: mix n within each (l,m) -> Gaunt selection rules kept
    U = np.zeros((ns, ns))
    groups = {}
    for i, (n, l, m) in enumerate(states):
        groups.setdefault((l, m), []).append(i)
    radial_U = {}
    for (l, m), ids in groups.items():
        if l not in radial_U:
            radial_U[l] = random_orth(len(ids))
        U[np.ix_(ids, ids)] = radial_U[l]
    h2, V2 = rotate(h, V, U)
    res.append(analyze(f"He radial-rotated (Gaunt kept) n_max={nmax}", singlet_H(h2, V2)))
    Uf = random_orth(ns)
    h3, V3 = rotate(h, V, Uf)
    res.append(analyze(f"He fully-rotated (Gaunt lost) n_max={nmax}", singlet_H(h3, V3)))
    # FCIQMC space: LatticeIndex spin-orbital determinants
    Hsd = idx.assemble_hamiltonian()
    res.append(analyze(f"He determinant space (FCIQMC) n_max={nmax}", Hsd))
    return res


# ---------------------------------------------------------------- (c)
def lih_composed():
    print("\n=== (c) composed LiH, JW qubit Hamiltonian, N=4 sector ===")
    from geovac.molecular_spec import lih_spec
    from geovac.composed_qubit import build_composed_hamiltonian
    import openfermion as of
    t0 = time.time()
    out = build_composed_hamiltonian(lih_spec())
    Q, fop = out["Q"], out["fermion_op"]
    print(f"Q={Q} N_pauli={out['N_pauli']} build {time.time()-t0:.1f}s")
    t0 = time.time()
    H = of.get_number_preserving_sparse_operator(fop, Q, 4, spin_preserving=False)
    H = sp.csr_matrix(H.real)
    print(f"N=4 sector dim={H.shape[0]} nnz={H.nnz} ({time.time()-t0:.1f}s)")
    res = [analyze("LiH composed Q=30 N=4 (global-min sector)", H, n_rand=3)]
    # --- physical sector: per-spec-block electron counts, Sz = 0 ---
    from openfermion.linalg.sparse_tools import _iterate_basis_
    ref = np.array([i < 4 for i in range(Q)])
    basis = np.array(list(_iterate_basis_(ref, 4, False)), dtype=int)
    spec = lih_spec()
    sub = out["blocks"]
    sizes = [b["n_orbitals"] for b in sub]
    # group sub-blocks into spec blocks
    groups, k = [], 0
    for blk in spec.blocks:
        nsub = 2 if getattr(blk, "has_h_partner", False) else 1
        groups.append((blk.label, blk.n_electrons, list(range(k, k + nsub))))
        k += nsub
    offs = np.cumsum([0] + sizes)
    ncomp, lab = connected_components(abs(H) > 0, directed=False)
    for conv in ("interleaved", "blocked"):
        M = Q // 2
        if conv == "interleaved":
            spat = np.arange(Q) // 2
            spin = np.arange(Q) % 2
        else:
            spat = np.arange(Q) % M
            spin = np.arange(Q) // M
        occ_blocks = []
        for (lbl, ne, subs) in groups:
            orb = np.zeros(M, bool)
            for s_ in subs:
                orb[offs[s_]:offs[s_ + 1]] = True
            occ_blocks.append(basis[:, orb[spat]].sum(1))
        sz = basis[:, spin == 0].sum(1) - basis[:, spin == 1].sum(1)
        labels = np.stack(occ_blocks + [sz], 1)
        # consistency: label constant inside each connected component?
        bad = 0
        for c_ in range(ncomp):
            ids = np.where(lab == c_)[0]
            if len(ids) > 1 and np.any(labels[ids] != labels[ids[0]]):
                bad += 1
        target = np.array([g[1] for g in groups] + [0])
        phys = np.where(np.all(labels == target, 1))[0]
        print(f"spin convention '{conv}': components with non-constant block labels={bad}; "
              f"physical-sector dets={len(phys)}; groups={[(g[0], g[1]) for g in groups]}")
        if bad == 0 and len(phys):
            Hp = H[phys][:, phys]
            nc, lb = connected_components(abs(Hp) > 0, directed=False)
            print(f"  physical sector splits into {nc} components (further M_L/per-block Sz)")
            res.append(analyze(f"LiH composed PHYSICAL sector ({conv})", Hp, n_rand=3))
            break
    return res


def ml_truncation_check():
    """Hybrid h1 carries L+- (m -> m+-1) graph edges, which do not conserve M_L.
    Compare E0 of the He singlet with all M_L vs restricted to M_L = 0 configs
    (the build_graph_native_fci enumeration)."""
    print("\n=== (b2) M_L conservation of the hybrid h1 ===")
    for nmax in (2, 3, 4):
        idx, h, V, states = he_tensors(nmax)
        ns = h.shape[0]
        mvals = np.array([m for (n, l, m) in states])
        C = sp.coo_matrix(h - np.diag(np.diag(h)))
        brk = int(np.sum(mvals[C.row] != mvals[C.col]))
        Hs = singlet_H(h, V)
        pairs = [(a, b) for a in range(ns) for b in range(a, ns)]
        ml = np.array([mvals[a] + mvals[b] for a, b in pairs])
        keep = np.where(ml == 0)[0]
        Eall = np.linalg.eigvalsh(Hs)[0]
        E0ml = np.linalg.eigvalsh(Hs[np.ix_(keep, keep)])[0]
        hc = h.copy()
        for i in range(ns):
            for j in range(ns):
                if mvals[i] != mvals[j]:
                    hc[i, j] = 0.0
        Econs = np.linalg.eigvalsh(singlet_H(hc, V))[0]
        print(f"n_max={nmax}: h1 off-diag entries with m!=m' = {brk}; "
              f"E0(all M_L)={Eall:.6f}  E0(M_L=0 restricted)={E0ml:.6f}  "
              f"E0(m-hops removed)={Econs:.6f}")


# ---------------------------------------------------------------- P2
def heat_kernel():
    print("\n=== P2: Fock-graph heat kernel closed form ===")
    from geovac.lattice import GeometricLattice
    from geovac.lattice_spectrum import path_eigenvectors, path_spectrum
    nmax = 5
    lat = GeometricLattice(max_n=nmax, nuclear_charge=1)
    A = lat.adjacency.toarray()
    L = np.diag(A.sum(1)) - A
    N = len(lat.states)
    for t in (0.1, 1.0, 10.0):
        K_num = expm(-t * L)
        K_cf = np.zeros((N, N))
        for l in range(nmax):
            ids = [i for i, (n, ll, m) in enumerate(lat.states) if ll == l]
            ns_, ms_ = nmax - l, 2 * l + 1

            def pk(mm):
                if mm == 1:
                    return np.ones((1, 1))
                V = path_eigenvectors(mm)
                lam = path_spectrum(mm)
                return V @ np.diag(np.exp(-t * lam)) @ V.T
            Kl = np.kron(pk(ns_), pk(ms_))
            # order ids by (n, m) to match kron (n outer, m inner)
            key = sorted(ids, key=lambda i: (lat.states[i][0], lat.states[i][2]))
            K_cf[np.ix_(key, key)] = Kl
        err = np.max(np.abs(K_num - K_cf))
        rows = np.max(np.abs(K_num.sum(1) - 1))
        print(f"t={t}: max|expm - closed form|={err:.2e}; row-sum-1 dev={rows:.1e}; "
              f"min entry={K_num.min():.2e} (Markov kernel: >=0)")
    # one-body reference propagator for the hybrid h1 after the (-1)^(n+m) gauge
    d = np.array([-0.5 / n ** 2 for (n, l, m) in lat.states])
    H1 = np.diag(d) + A / 16.0
    g = np.array([(-1.0) ** (n + m) for (n, l, m) in lat.states])
    G = np.diag(g)
    P = expm(-1.0 * (G @ H1 @ G))
    print(f"gauged hybrid h1: max offdiag={np.max(G@H1@G - np.diag(d)):.3e} "
          f"(<=0 stoquastic); e^(-tau h1) min entry={P.min():.3e} (>=0 => sign-free one-body walk)")
    # ungauged propagator has negative entries
    P0 = expm(-1.0 * H1)
    print(f"ungauged e^(-tau h1) negative-entry fraction={np.mean(P0 < -1e-15):.3f}")


if __name__ == "__main__":
    T0 = time.time()
    fock_one_body()
    heat_kernel()
    for nmax in (2, 3, 4):
        he_block(nmax)
    ml_truncation_check()
    try:
        lih_composed()
    except Exception as exc:  # report, do not hide
        print(f"LiH section FAILED: {type(exc).__name__}: {exc}")
    print(f"\nTOTAL {time.time()-T0:.1f}s")
