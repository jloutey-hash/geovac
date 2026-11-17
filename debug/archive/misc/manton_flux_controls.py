"""Non-scored controls and applicability checks for the DRAFT pre-registration MF-Hopf ("Manton flux test",
debug/manton_flux_prereg_draft.md).

Modes:
  python debug/manton_flux_controls.py manton      # C1: Manton's 24-cell U(1) Hopf bundle (CMP 113, 341 (1987))
  python debug/manton_flux_controls.py trivial     # C2: product bundle octahedron x Z4, same machinery
  python debug/manton_flux_controls.py applicable  # A1-A4: does the GeoVac label graph admit Manton's structure?
  python debug/manton_flux_controls.py all

Nothing here computes a flux, holonomy or Chern sum on the GeoVac graph.
"""
from __future__ import annotations

import itertools
import os
import sys
from typing import Dict, List, Tuple

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

Z4 = [1, 1j, -1, -1j]
PHASE = {0: 0.0, 1: np.pi / 2, 2: np.pi, 3: -np.pi / 2}  # Manton's branch: 1, i, -1, -i -> 0, pi/2, pi, -pi/2


def z4_index(g: complex) -> int:
    for k, h in enumerate(Z4):
        if abs(g - h) < 1e-9:
            return k
    raise ValueError(g)


# ------------------------------------------------------------------ generic discrete-bundle machinery
class DiscreteBundle:
    """Total space: points P (any labels), a symmetric neighbour relation, a free Z4 action act(p, k)."""

    def __init__(self, points: List, nbr: Dict[int, set], act, fibre_of: List[int], base_xyz: np.ndarray):
        self.P, self.nbr, self.act, self.fib, self.xyz = points, nbr, act, fibre_of, base_xyz
        nf = max(fibre_of) + 1
        self.fibres = [[i for i in range(len(points)) if fibre_of[i] == f] for f in range(nf)]
        self.rep = [f[0] for f in self.fibres]
        links = set()
        for i in range(len(points)):
            for j in nbr[i]:
                if fibre_of[i] != fibre_of[j]:
                    links.add(tuple(sorted((fibre_of[i], fibre_of[j]))))
        self.links = sorted(links)
        adj = {f: set() for f in range(nf)}
        for a, b in self.links:
            adj[a].add(b)
            adj[b].add(a)
        tri = []
        for a, b, c in itertools.combinations(range(nf), 3):
            if b in adj[a] and c in adj[a] and c in adj[b]:
                # orient outward: (xb - xa) x (xc - xa) . xa > 0
                if np.dot(np.cross(base_xyz[b] - base_xyz[a], base_xyz[c] - base_xyz[a]), base_xyz[a]) < 0:
                    b, c = c, b
                tri.append((a, b, c))
        self.tri = tri
        # lift choices on link (a, b), from the representative of a
        self.choices = {(a, b): sorted(j for j in nbr[self.rep[a]] if fibre_of[j] == b) for (a, b) in self.links}

    def group_elt(self, p: int, q: int) -> int:
        """k with q = act(p, k) (p, q on one fibre)."""
        for k in range(4):
            if self.act(p, k) == q:
                return k
        raise ValueError("not on one fibre")

    def holonomy(self, conn, loop: Tuple[int, ...]) -> int:
        p0 = self.rep[loop[0]]
        p = p0
        for b in list(loop[1:]) + [loop[0]]:
            p = self.step(conn, p, b)
        return self.group_elt(p0, p)

    def step(self, conn, p: int, b: int) -> int:
        a = self.fib[p]
        k = self.group_elt(self.rep[a], p)
        if (a, b) in conn:
            return self.act(self.choices[(a, b)][conn[(a, b)]], k)
        # reverse: find the point q on fibre b whose chosen lift b->a... equivalently the lift a->b reversed:
        # the chosen lift of (b, a) maps rep[b] -> s on fibre a; s = act(rep[a], k_s); so from rep[a] the reversed
        # lift goes to act(rep[b], -k_s), and from act(rep[a], k) to act(rep[b], k - k_s).
        s = self.choices[(b, a)][conn[(b, a)]]
        ks = self.group_elt(self.rep[a], s)
        return self.act(self.rep[b], (k - ks) % 4)

    def connections(self):
        keys = self.links
        for bits in itertools.product(range(2), repeat=len(keys)):
            yield dict(zip(keys, bits))

    def global_sections(self) -> int:
        nf = len(self.fibres)
        count = 0
        for sel in itertools.product(*self.fibres):
            ok = all((sel[b] in self.nbr[sel[a]]) for a, b in self.links)
            count += ok
        return count


def report_bundle(B: DiscreteBundle, name: str, orient_sign: int = 1) -> Dict:
    two = all(len(v) == 2 for v in B.choices.values())
    chern, flat, hol0 = {}, 0, {}
    for conn in B.connections():
        hs = [B.holonomy(conn, t) for t in B.tri]
        if all(h == 0 for h in hs):
            flat += 1
        if orient_sign == 1:
            c = sum(PHASE[h] for h in hs) / (2 * np.pi)
        elif orient_sign == -1:   # inward faces, Manton's outward branch reused naively on inverse holonomies
            c = sum(PHASE[(-h) % 4] for h in hs) / (2 * np.pi)
        else:                     # inward faces with the mirrored branch (phase(h^-1) = -phase(h))
            c = sum(-PHASE[h] for h in hs) / (2 * np.pi)
        c = round(c, 9)
        chern[c] = chern.get(c, 0) + 1
        hol0[hs[0]] = hol0.get(hs[0], 0) + 1
    out = dict(points=len(B.P), fibres=len(B.fibres), links=len(B.links), triangles=len(B.tri),
               two_lifts_per_link=two, connections=sum(chern.values()), flat=flat,
               chern_hist=chern, first_triangle_holonomy_hist={str(Z4[k]): v for k, v in sorted(hol0.items())},
               global_sections=B.global_sections())
    print(f"[{name}]")
    for k, v in out.items():
        print(f"  {k}: {v}")
    return out


# ------------------------------------------------------------------ C1: Manton's 24-cell
def manton_bundle() -> DiscreteBundle:
    pts = []
    for i in range(4):
        for s in (1, -1):
            v = np.zeros(4)
            v[i] = s
            pts.append(v)
    for signs in itertools.product((1, -1), repeat=4):
        pts.append(0.5 * np.array(signs, float))
    zs = [(p[0] + 1j * p[1], p[2] + 1j * p[3]) for p in pts]
    n = len(zs)
    assert n == 24 and all(abs(abs(a) ** 2 + abs(b) ** 2 - 1) < 1e-12 for a, b in zs)

    def find(z) -> int:
        for k, w in enumerate(zs):
            if abs(w[0] - z[0]) < 1e-9 and abs(w[1] - z[1]) < 1e-9:
                return k
        raise ValueError

    actt = {(p, k): find((zs[p][0] * Z4[k], zs[p][1] * Z4[k])) for p in range(n) for k in range(4)}
    act = lambda p, k: actt[(p, k % 4)]  # noqa: E731
    assert all(actt[(p, k)] != p for p in range(n) for k in (1, 2, 3)), "Z4 not free"
    nbr = {i: {j for j in range(n) if j != i and abs(np.linalg.norm(pts[i] - pts[j]) - 1) < 1e-9} for i in range(n)}
    assert all(len(v) == 8 for v in nbr.values())
    assert all(actt[(j, k)] in nbr[actt[(i, k)]] for i in range(n) for j in nbr[i] for k in range(4)), "links not equivariant"
    # fibres and base points Z = z1/z2 on the Riemann sphere
    fib, seen = [-1] * n, 0
    for p in range(n):
        if fib[p] < 0:
            for k in range(4):
                fib[actt[(p, k)]] = seen
            seen += 1
    xyz = np.zeros((seen, 3))
    for p in range(n):
        z1, z2 = zs[p]
        # Hopf map to S^2: (2 Re z1 conj z2, 2 Im z1 conj z2, |z1|^2 - |z2|^2)
        w = z1 * np.conj(z2)
        xyz[fib[p]] = [2 * w.real, 2 * w.imag, abs(z1) ** 2 - abs(z2) ** 2]
    return DiscreteBundle(zs, nbr, act, fib, xyz)


def trivial_bundle() -> DiscreteBundle:
    """Octahedron x Z4, link (v, g) ~ (v', g) and (v, g) ~ (v', g i) for each base edge v < v' (oriented)."""
    xyz = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]], float)
    edges = [(a, b) for a, b in itertools.combinations(range(6), 2) if abs(np.dot(xyz[a], xyz[b])) < 0.5]
    pts = [(v, k) for v in range(6) for k in range(4)]
    idx = {p: i for i, p in enumerate(pts)}
    nbr = {i: set() for i in range(len(pts))}
    for a, b in edges:
        for k in range(4):
            for dk in (0, 1):
                i, j = idx[(a, k)], idx[(b, (k + dk) % 4)]
                nbr[i].add(j)
                nbr[j].add(i)
    act = lambda p, k: idx[(pts[p][0], (pts[p][1] + k) % 4)]  # noqa: E731
    fib = [p[0] for p in pts]
    return DiscreteBundle(pts, nbr, act, fib, xyz)


# ------------------------------------------------------------------ A1-A4: applicability on the GeoVac graph
def applicability(nmax: int) -> Dict:
    from geovac.lattice import GeometricLattice
    lat = GeometricLattice(max_n=nmax)
    S = list(lat.states)
    A = (np.asarray(lat.adjacency.toarray()) != 0).astype(int)
    N0 = len(S)
    E = [(i, j) for i in range(N0) for j in range(i + 1, N0) if A[i, j]]
    # A1: square 2-cells (all 4-cycles; the graph is triangle-free) and the homology of the square complex
    tri = sum(1 for i, j, k in itertools.combinations(range(N0), 3) if A[i, j] and A[j, k] and A[i, k])
    sq = set()
    for a, b, c, d in itertools.permutations(range(N0), 4):
        if A[a, b] and A[b, c] and A[c, d] and A[d, a] and not A[a, c] and not A[b, d]:
            cyc = (a, b, c, d)
            r = min(range(4), key=lambda t: cyc[t])
            cyc = cyc[r:] + cyc[:r]
            if cyc[1] > cyc[3]:
                cyc = (cyc[0], cyc[3], cyc[2], cyc[1])
            sq.add(cyc)
    sq = sorted(sq)
    eidx = {e: k for k, e in enumerate(E)}
    d1 = np.zeros((N0, len(E)), int)
    for k, (i, j) in enumerate(E):
        d1[i, k], d1[j, k] = -1, 1
    d2 = np.zeros((len(E), len(sq)), int)
    for f, cyc in enumerate(sq):
        for t in range(4):
            u, v = cyc[t], cyc[(t + 1) % 4]
            d2[eidx[(min(u, v), max(u, v))], f] = 1 if u < v else -1
    assert not (d1 @ d2).any()
    r1 = np.linalg.matrix_rank(d1)
    r2 = np.linalg.matrix_rank(d2) if sq else 0
    b0 = N0 - r1
    b1 = len(E) - r1 - r2
    b2 = len(sq) - r2
    # A2: n-preserving graph automorphisms and whether any acts without fixed nodes
    import networkx as nx
    from networkx.algorithms.isomorphism import GraphMatcher
    G = nx.Graph()
    for i, s in enumerate(S):
        G.add_node(i, n=s[0])
    G.add_edges_from(E)
    auts = list(GraphMatcher(G, G, node_match=lambda x, y: x["n"] == y["n"]).isomorphisms_iter())
    free = [a for a in auts if all(a[i] != i for i in range(N0))]
    # A3: does the Hopf-fibre U(1), exp(i phi (L_z + A_z)), or exp(i phi L_z), permute the labels?
    import importlib.util
    spec = importlib.util.spec_from_file_location("so4_guard", os.path.join(ROOT, "tests", "test_atlas_so4_commutant.py"))
    so4 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(so4)
    if nmax != so4.NMAX:
        so4.NMAX = nmax
        so4.STATES = [(n, l, m) for n in range(1, nmax + 1) for l in range(n) for m in range(-l, l + 1)]
        so4.IDX = {st: i for i, st in enumerate(so4.STATES)}
        so4.DIM = len(so4.STATES)
    assert so4.STATES == S
    L, Av = so4.generators()
    phi = 0.7

    def monomial_cols(G_: np.ndarray) -> int:
        w, V = np.linalg.eigh(G_)
        U = V @ np.diag(np.exp(1j * phi * w)) @ V.conj().T
        return int(max((np.abs(U[:, c]) > 1e-10).sum() for c in range(U.shape[1])))
    out = dict(nmax=nmax, nodes=N0, edges=len(E), triangles=tri, squares=len(sq),
               b0=b0, b1_graph=len(E) - r1,
               b1_square_complex=b1, b2_square_complex=b2,
               n_preserving_automorphisms=len(auts), fixed_point_free=len(free),
               fibre_LzAz_max_support_per_column=monomial_cols(L[2] + Av[2]),
               Lz_max_support_per_column=monomial_cols(L[2]))
    print(f"[applicability n_max={nmax}]")
    for k, v in out.items():
        print(f"  {k}: {v}")
    return out


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    if mode in ("manton", "all"):
        Bm = manton_bundle()
        report_bundle(Bm, "C1 Manton 24-cell, outward orientation")
        report_bundle(Bm, "C1' Manton 24-cell, inward faces, outward branch reused", orient_sign=-1)
        report_bundle(Bm, "C1'' Manton 24-cell, inward faces, mirrored branch", orient_sign=0)
    if mode in ("trivial", "all"):
        Bt = trivial_bundle()
        report_bundle(Bt, "C2 trivial product bundle octahedron x Z4 (Manton's branch)")
        # C2b: additive flux relative to the global section g = 0 (exists here, not on Manton's bundle)
        tot = set()
        for conn in Bt.connections():
            a = {}
            for (u, v), bit in conn.items():
                dk = Bt.P[Bt.choices[(u, v)][bit]][1] - Bt.P[Bt.rep[u]][1]
                a[(u, v)], a[(v, u)] = dk % 4, -(dk % 4)
            tot.add(sum(a[(x, y)] + a[(y, z)] + a[(z, x)] for x, y, z in Bt.tri))
        print(f"[C2b additive flux sum over 8 faces, in quarter-turns, all 4096 connections] values: {sorted(tot)}")
    if mode in ("applicable", "all"):
        for nm in (3, 4):
            applicability(nm)
