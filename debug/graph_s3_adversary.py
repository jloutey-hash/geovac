"""Adversarial checks on the graph->S^3 refutation (2026-09-25).

(a) components/zero modes, binary + topological weights
(b) fock_graph_hodge docstring claim: L0/L1 eigenvalues integers = n^2-1?
(c) topological weights docstring: 'D-A eigenvalues scaling as -1/n^2'?
(d) He graph-native CI: effect of the graph off-diagonals (on vs zeroed)
(e) label-lattice SO(4) Casimir: diagonal n^2-1 (trivial check)
"""
import os
for v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[v] = "1"
import sys
import warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import numpy as np
from scipy.sparse.csgraph import connected_components
from geovac.lattice import GeometricLattice
import geovac.casimir_ci as cc

print("== (a)/(c) components and spectra")
for tw in (False, True):
    for nm in (3, 6):
        lat = GeometricLattice(max_n=nm, topological_weights=tw)
        A = lat.adjacency.toarray()
        L = np.diag(A.sum(1)) - A
        ncomp, _ = connected_components(lat.adjacency, directed=False)
        ev = np.linalg.eigvalsh(L)
        nz = int(np.sum(np.abs(ev) < 1e-10))
        print(f"tw={tw} n_max={nm} V={len(ev)} comps={ncomp} zero-modes={nz} "
              f"integer-eigs={int(np.sum(np.abs(ev-np.round(ev))<1e-9))}/{len(ev)} "
              f"top5={np.round(np.sort(ev)[-5:],4)}")

print("== (b) hodge module")
try:
    import geovac.fock_graph_hodge as fh
    names = [n for n in dir(fh) if not n.startswith("_")]
    print("public:", names)
except Exception as e:  # noqa
    print("import fail", e)

print("== (e) SO(4) Casimir on label lattice: L^2 + A^2 per shell = n^2-1 (diag)")
for nm in (4,):
    lat = GeometricLattice(max_n=nm)
    print(sorted(set((n, n * n - 1) for n, l, m in lat.states)))

print("== (d) He graph-native CI, graph off-diagonals ON vs ZEROED")
E_exact = -2.903724377
orig = cc._build_graph_h1
def diag_only(Z, n_max):
    h, orb = orig(Z, n_max)
    return np.diag(np.diag(h)), orb
for nm in (3, 4, 5):
    cc._build_graph_h1 = orig
    H = cc.build_graph_native_fci(2, nm)
    e_on = np.linalg.eigvalsh(H)[0]
    cc._build_graph_h1 = diag_only
    H0 = cc.build_graph_native_fci(2, nm)
    e_off = np.linalg.eigvalsh(H0)[0]
    cc._build_graph_h1 = orig
    print(f"n_max={nm} dim={H.shape[0]} E_on={e_on:.6f} ({100*(e_on-E_exact)/abs(E_exact):+.4f}%) "
          f"E_diag={e_off:.6f} ({100*(e_off-E_exact)/abs(E_exact):+.4f}%)")
