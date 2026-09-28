"""Independent adversarial verification of the heat-kernel-swap W1-W4 legs
and their scope across GeoVac graph constructions (2026-09-25).

Written from scratch; does NOT import debug/heat_kernel_swap_probe.py.
Routes:
  (a) graph-algorithmic: read the production adjacency, scipy csgraph
      connected components, dense per-component diagonalisation;
  (b) closed form: each l-block is P_{n_max-l} x P_{2l+1} (free-end paths),
      spectrum {(2-2cos(j pi/N1)) + (2-2cos(k pi/N2))}.
Route (b) is cross-checked against route (a) where (a) is affordable.
"""
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"

import sys
import time
from typing import Dict, List, Tuple

import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components
from scipy.linalg import expm
from scipy.optimize import minimize_scalar

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from geovac.lattice import GeometricLattice  # noqa: E402


def lap(A: sp.spmatrix) -> np.ndarray:
    A = sp.csr_matrix(A)
    d = np.asarray(A.sum(axis=1)).ravel()
    return (sp.diags(d) - A).toarray()


def graph_report(A: sp.spmatrix, lvals: np.ndarray, name: str,
                 dense_ok: bool = True) -> Dict:
    A = sp.csr_matrix(A)
    coo = sp.triu(A, k=1).tocoo()
    cross_l = int(np.sum(lvals[coo.row] != lvals[coo.col]))
    ncomp, _ = connected_components(A, directed=False)
    kerdim = None
    if dense_ok:
        ev = np.linalg.eigvalsh(lap(A))
        kerdim = int(np.sum(np.abs(ev) < 1e-9))
    print(f"  {name:48s} V={A.shape[0]:6d} E={coo.nnz:7d} cross-l={cross_l:6d} "
          f"components={ncomp:4d} dimker={kerdim}")
    return dict(cross_l=cross_l, ncomp=ncomp, kerdim=kerdim)


def closed_form_spectrum(n_max: int) -> np.ndarray:
    out = []
    for l in range(n_max):
        N1, N2 = n_max - l, 2 * l + 1
        a = 2 - 2 * np.cos(np.arange(N1) * np.pi / N1)
        b = 2 - 2 * np.cos(np.arange(N2) * np.pi / N2)
        out.append((a[:, None] + b[None, :]).ravel())
    return np.sort(np.concatenate(out))


def count_top(n_max: int, thresh: float) -> int:
    """#{eigs with 8 - lam <= thresh}, closed form, O(N log N) per block."""
    tot = 0
    for l in range(n_max):
        N1, N2 = n_max - l, 2 * l + 1
        x = np.sort(2 + 2 * np.cos(np.arange(N1) * np.pi / N1))  # 4 - a
        y = np.sort(2 + 2 * np.cos(np.arange(N2) * np.pi / N2))
        tot += int(np.searchsorted(y, thresh - x, side="right").sum())
    return tot


def main() -> None:
    t0 = time.time()
    print("=" * 78)
    print("W1: components / cross-l edges / dim ker  (production GeometricLattice)")
    print("=" * 78)
    for nm in (2, 3, 5, 8, 12):
        for tw in (False, True):
            lat = GeometricLattice(max_n=nm, topological_weights=tw)
            lv = np.array([s[1] for s in lat.states])
            r = graph_report(lat.adjacency, lv,
                             f"GeometricLattice n={nm} topo_w={tw}")
            assert r["cross_l"] == 0 and r["ncomp"] == nm and r["kerdim"] == nm

    # Block = grid P_{n-l} x P_{2l+1}: compare each component's spectrum
    # (route a) with the closed form (route b).
    print("\nBlock-structure check (per-component dense eig vs closed form):")
    for nm in (5, 10, 20, 30):
        lat = GeometricLattice(max_n=nm)
        A = lat.adjacency.tocsr()
        ncomp, lab = connected_components(A, directed=False)
        allev = []
        for c in range(ncomp):
            idx = np.where(lab == c)[0]
            allev.append(np.linalg.eigvalsh(lap(A[idx][:, idx])))
        allev = np.sort(np.concatenate(allev))
        cf = closed_form_spectrum(nm)
        dev = float(np.max(np.abs(allev - cf)))
        n_gt2 = int(np.sum(allev > 2 + 1e-12))
        print(f"  n_max={nm:3d}  V={len(cf):6d}  max|direct-closed|={dev:.2e}  "
              f"#lam>2 (direct)={n_gt2}  lam_max={allev[-1]:.10f}")
    # Full-matrix dense check at small n_max (no component split at all)
    for nm in (5, 12):
        lat = GeometricLattice(max_n=nm)
        ev = np.linalg.eigvalsh(lap(lat.adjacency))
        print(f"  full dense n_max={nm}: max|eig-closed|="
              f"{np.max(np.abs(ev - closed_form_spectrum(nm))):.2e}")

    print("\n" + "=" * 78)
    print("W2: top-edge accumulation and the kappa=-1/16 hydrogen count")
    print("=" * 78)
    f = lambda nu: np.pi ** 2 * ((1 - nu) ** -2 + (2 * nu) ** -2)
    res = minimize_scalar(f, bounds=(1e-3, 1 - 1e-3), method="bounded",
                          options={"xatol": 1e-12})
    C = float(res.fun)
    print(f"  C = min_nu pi^2[(1-nu)^-2+(2nu)^-2] = {C:.6f} at nu*={res.x:.6f} "
          f"(1/(1+2^(2/3)) = {1/(1+2**(2/3)):.6f})")
    for nm in (5, 30, 100):
        ev = closed_form_spectrum(nm)
        print(f"  n_max={nm}: #lam>2 strict = {int(np.sum(ev > 2 + 1e-9))}; "
              f"#lam>=2 = {int(np.sum(ev > 2 - 1e-9))}; "
              f"#lam==2 (ties) = {int(np.sum(np.abs(ev - 2) < 1e-9))}")
    print("  n_max   n^2(8-lam_max)   #(n^2(8-lam)<=1.5C)/n   #(<=2C)/n   top5 rescaled")
    for nm in (50, 100, 200, 400, 800):
        c15 = count_top(nm, 1.5 * C / nm ** 2)
        c2 = count_top(nm, 2.0 * C / nm ** 2)
        # top five rescaled eigenvalues
        vals = []
        for l in range(nm):
            N1, N2 = nm - l, 2 * l + 1
            x = 2 + 2 * np.cos(np.arange(N1 - 3, N1) * np.pi / N1)
            y = 2 + 2 * np.cos(np.arange(max(N2 - 3, 0), N2) * np.pi / N2)
            vals.extend((x[:, None] + y[None, :]).ravel().tolist())
        top5 = np.sort(np.array(vals))[:5] * nm ** 2
        print(f"  {nm:5d}   {top5[0]:12.6f}   {c15/nm:10.4f}   {c2/nm:10.4f}   "
              + " ".join(f"{v:.3f}" for v in top5))

    # hydrogen via the production solver
    from geovac.atomic_solver import AtomicSolver
    for nm in (5, 10, 30):
        s = AtomicSolver(max_n=nm, Z=1)
        E, _ = s.compute_ground_state(n_states=3)
        print(f"  AtomicSolver n_max={nm}: lowest 3 energies = "
              + ", ".join(f"{e:.6f}" for e in np.sort(E))
              + "   (exact: -0.5, -0.125 x4)")

    print("\n" + "=" * 78)
    print("W3: spectral dimension d_s(t) = -2 dlnTr/dlnt (closed-form spectrum)")
    print("=" * 78)
    for nm in (50, 200):
        ev = closed_form_spectrum(nm)
        row = []
        for t in (1.0, 3.0, 10.0, 30.0, 100.0, 1000.0):
            Z = np.sum(np.exp(-t * ev))
            dZ = -np.sum(ev * np.exp(-t * ev))
            row.append(f"t={t:g}:{-2 * t * dZ / Z:.3f}")
        print(f"  n_max={nm}: " + "  ".join(row))
    # S^3 reference
    n = np.arange(1, 200001, dtype=float)
    row = []
    for tau in (1e-4, 1e-3, 1e-2):
        w = n ** 2 * np.exp(-tau * (n ** 2 - 1))
        row.append(f"tau={tau:g}:{2 * tau * np.sum(w * (n**2 - 1)) / np.sum(w):.3f}")
    print("  S^3 (sum n^2 e^{-tau(n^2-1)}): " + "  ".join(row))

    print("\n" + "=" * 78)
    print("W4: heat kernel between blocks; same-block pair fraction")
    print("=" * 78)
    lat = GeometricLattice(max_n=6)
    K = expm(-1.0 * lap(lat.adjacency))
    lv = np.array([s[1] for s in lat.states])
    cross = lv[:, None] != lv[None, :]
    print(f"  n_max=6, t=1: max |K_ij| over cross-l pairs = {np.max(np.abs(K[cross])):.2e}; "
          f"min K_ij over same-l = {np.min(K[~cross]):.2e}")
    for nm in (10, 100, 1000, 10000):
        l = np.arange(nm)
        Vl = (nm - l) * (2 * l + 1.0)
        frac = np.sum(Vl ** 2) / np.sum(Vl) ** 2
        print(f"  n_max={nm:6d}: n_max * same-block fraction = {nm * frac:.5f} (claim 6/5=1.2)")

    print("\n" + "=" * 78)
    print("SCOPE: other constructions")
    print("=" * 78)
    # LatticeIndex one-body operators
    try:
        import warnings
        warnings.simplefilter("ignore")
        from geovac.lattice_index import LatticeIndex
        for meth, vee in (("graph", "chordal"), ("hybrid", "slater"),
                          ("exact", "slater_full")):
            li = LatticeIndex(n_electrons=2, max_n=4, nuclear_charge=2,
                              vee_method=vee, h1_method=meth)
            H1 = sp.csr_matrix(li._H1_spatial)
            lv = np.array([s[1] for s in li.lattice.states])
            mv = np.array([s[2] for s in li.lattice.states])
            off = sp.triu(H1, k=1).tocoo()
            off_nz = np.abs(off.data) > 1e-14
            cl = int(np.sum(lv[off.row][off_nz] != lv[off.col][off_nz]))
            dm = int(np.sum(mv[off.row][off_nz] != mv[off.col][off_nz]))
            ncomp, _ = connected_components((abs(H1) > 1e-14).astype(float),
                                            directed=False)
            print(f"  LatticeIndex h1={meth:6s} n_max=4: offdiag={int(off_nz.sum())} "
                  f"cross-l={cl} m-changing={dm} components={ncomp}")
    except Exception as exc:  # pragma: no cover
        print(f"  LatticeIndex check failed: {exc!r}")

    # Native Dirac graphs
    from geovac.dirac_lattice import DiracLattice
    from geovac.ihara_zeta_dirac import build_dirac_s3_graph
    for nm in (2, 3, 5, 8):
        for mode in ("atomic", "s3"):
            dl = DiracLattice(n_max=nm, mode=mode)
            lv = np.array([lab.l for lab in dl.labels])
            graph_report(dl.adjacency, lv, f"DiracLattice(RuleB) n={nm} mode={mode}")
        for rule in ("A", "B"):
            A, labels, _, _ = build_dirac_s3_graph(nm, rule)
            lv = np.array([lab.l for lab in labels])
            graph_report(sp.csr_matrix(A.astype(float)), lv,
                         f"ihara build_dirac_s3_graph rule={rule} n={nm}")
    # Rule-B Dirac graph Laplacian spectrum: any S^3 structure?
    dl = DiracLattice(n_max=5, mode="s3")
    ev = np.linalg.eigvalsh(lap(dl.adjacency))
    print(f"  DiracLattice n=5 s3: L spectrum min/2nd/max = {ev[0]:.3e}/{ev[1]:.4f}/{ev[-1]:.4f}; "
          f"#integer eigs = {int(np.sum(np.abs(ev - np.round(ev)) < 1e-9))}/{len(ev)}")

    # Bargmann-Segal S^5 graph
    from geovac.nuclear.bargmann_graph import build_bargmann_graph
    for Nm in (2, 3, 5):
        g = build_bargmann_graph(Nm)
        Ad = g.adjacency_dense()
        lv = np.array([nd[1] for nd in g.nodes]) if hasattr(g, "nodes") else None
        if lv is None:
            from geovac.nuclear.bargmann_graph import enumerate_nodes
            lv = np.array([nd[1] for nd in enumerate_nodes(Nm)])
        graph_report(sp.csr_matrix(Ad), lv, f"BargmannGraph N_max={Nm}")

    # FockGraphHodge docstring claim: graph L0 spectrum = n^2-1 ?
    ev = np.linalg.eigvalsh(lap(GeometricLattice(max_n=3).adjacency))
    print(f"  GeometricLattice n=3 L spectrum: {np.round(ev, 6).tolist()}")
    print(f"  (continuum n^2-1 with multiplicity n^2 would be 0 x1, 3 x4, 8 x9)")

    print(f"\nwall time {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
