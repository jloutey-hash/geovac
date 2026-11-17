"""Forced checks and non-scored controls for the draft pre-registration HG
(debug/hopf_grouping_prereg_draft.md): Paper 25's base graph under the rotation grouping it uses and under a
Hopf (chiral) grouping of each shell.

Modes:
  python debug/hopf_grouping_controls.py forced    # SO(4) identities, fixed by symmetry (run)
  python debug/hopf_grouping_controls.py controls  # rotation recipe reproduces Paper 25 (run)
  python debug/hopf_grouping_controls.py score     # the Hopf base graph: SCORED, NOT RUN before registration

Generators: reused from tests/test_atlas_so4_commutant.py (hydrogenic (n, l, m), n <= 3, Condon-Shortley; L from
the ladder algebra, A from Pauli's shell identity). Graph: geovac.lattice.GeometricLattice(max_n=3), whose node
order equals the test's STATES order (checked below).
"""
from __future__ import annotations

import os
import sys
from fractions import Fraction
from typing import Dict, List, Tuple

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
import sympy as sp
from sympy.physics.quantum.cg import CG

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from geovac.lattice import GeometricLattice  # noqa: E402
import importlib.util  # noqa: E402

_spec = importlib.util.spec_from_file_location("so4_guard", os.path.join(ROOT, "tests", "test_atlas_so4_commutant.py"))
_so4 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_so4)
# Shell cutoff: HG_NMAX (default 3, the guard's own).  The guard's generators() reads its module globals, so
# rebinding them here builds the same SO(4) generators on more shells (added 2026-10-31 for the n_max = 4 scope).
_n = int(os.environ.get("HG_NMAX", "3"))
if _n != _so4.NMAX:
    _so4.NMAX = _n
    _so4.STATES = [(n, l, m) for n in range(1, _n + 1) for l in range(n) for m in range(-l, l + 1)]
    _so4.IDX = {st: i for i, st in enumerate(_so4.STATES)}
    _so4.DIM = len(_so4.STATES)
STATES, IDX, DIM, NMAX, generators, _comm = (_so4.STATES, _so4.IDX, _so4.DIM, _so4.NMAX, _so4.generators,
                                             _so4._comm)

Group = Tuple[int, sp.Rational]  # (n, label)


def graph() -> np.ndarray:
    lat = GeometricLattice(max_n=NMAX)
    assert list(lat.states) == STATES, "graph node order differs from the generator basis"
    return np.asarray(lat.adjacency.toarray(), dtype=float)


# ---------------------------------------------------------------- groupings as membership matrices
def rotation_M() -> Tuple[List[Group], sp.Matrix]:
    groups = sorted({(n, sp.Integer(l)) for (n, l, m) in STATES})
    M = sp.zeros(DIM, len(groups))
    for (n, l, m) in STATES:
        M[IDX[(n, l, m)], groups.index((n, sp.Integer(l)))] = 1
    return groups, M


def hopf_M_exact() -> Tuple[List[Group], sp.Matrix]:
    """Primary recipe: group (n, mu) = eigenspace of J^-_z = (L_z - A_z)/2 with eigenvalue mu inside shell n.
    Occupation of node (n, l, m) in group (n, mu) = |<j, m - mu; j, mu | l, m>|^2, j = (n - 1)/2 (exact)."""
    groups: List[Group] = []
    for n in range(1, NMAX + 1):
        j = sp.Rational(n - 1, 2)
        groups += [(n, -j + k) for k in range(n)]
    M = sp.zeros(DIM, len(groups))
    for (n, l, m) in STATES:
        j = sp.Rational(n - 1, 2)
        for g, (ng, mu) in enumerate(groups):
            if ng != n:
                continue
            mp = sp.Integer(m) - mu
            if abs(mp) > j:
                continue
            M[IDX[(n, l, m)], g] = sp.nsimplify(CG(j, mp, j, mu, l, m).doit() ** 2)
    return groups, M


def chiral_projectors(sign: int) -> Dict[Group, np.ndarray]:
    """Numerical route: eigenspace projectors of J^{sign}_z = (L_z + sign*A_z)/2 within each shell."""
    L, A = generators()
    Jz = (L[2] + sign * A[2]) / 2
    out: Dict[Group, np.ndarray] = {}
    for n in range(1, NMAX + 1):
        ids = [IDX[s] for s in STATES if s[0] == n]
        w, v = np.linalg.eigh(Jz[np.ix_(ids, ids)])
        for mu2 in sorted(set(int(round(2 * x)) for x in w)):
            sel = np.abs(2 * w - mu2) < 1e-8
            P = np.zeros((DIM, DIM), complex)
            P[np.ix_(ids, ids)] = v[:, sel] @ v[:, sel].conj().T
            out[(n, sp.Rational(mu2, 2))] = P
    return out


# ---------------------------------------------------------------- the two edge rules
def base_occupation(M: sp.Matrix, A: np.ndarray) -> sp.Matrix:
    """Primary edge rule: W_full = M^T A M; base adjacency = off-diagonal part."""
    W = M.T * sp.Matrix(A.astype(int)) * M
    for i in range(W.shape[0]):
        W[i, i] = 0
    return W


def base_frobenius(P: Dict[Group, np.ndarray], groups: List[Group], A: np.ndarray) -> np.ndarray:
    """Alternative edge rule: W[g, g'] = ||P_g A P_g'||_F^2, g != g'."""
    k = len(groups)
    W = np.zeros((k, k))
    for a, g in enumerate(groups):
        for b, h in enumerate(groups):
            if a != b:
                W[a, b] = float(np.linalg.norm(P[g] @ A @ P[h]) ** 2)
    return W


def laplacian(W: sp.Matrix) -> sp.Matrix:
    k = W.shape[0]
    Lp = -W
    for i in range(k):
        Lp[i, i] = sum(W[i, j] for j in range(k))
    return Lp


def coordinate_projectors(groups: List[Group], M: sp.Matrix) -> Dict[Group, np.ndarray]:
    out = {}
    for g, lab in enumerate(groups):
        P = np.zeros((DIM, DIM), complex)
        for u in range(DIM):
            if M[u, g] == 1:
                P[u, u] = 1.0
        out[lab] = P
    return out


# ---------------------------------------------------------------- modes
def forced() -> None:
    L, A = generators()
    eps = lambda i, j, k: np.sign((j - i) * (k - i) * (k - j))  # noqa: E731
    err = 0.0
    for i in range(3):
        for j in range(3):
            err = max(err, np.abs(_comm(L[i], L[j]) - 1j * sum(eps(i, j, k) * L[k] for k in range(3))).max(),
                      np.abs(_comm(L[i], A[j]) - 1j * sum(eps(i, j, k) * A[k] for k in range(3))).max(),
                      np.abs(_comm(A[i], A[j]) - 1j * sum(eps(i, j, k) * L[k] for k in range(3))).max())
    print(f"F0 so(4) closure, max error {err:.1e}")
    Jp = [(L[a] + A[a]) / 2 for a in range(3)]
    Jm = [(L[a] - A[a]) / 2 for a in range(3)]
    L2 = sum(G @ G for G in L)
    Jp2 = sum(G @ G for G in Jp)
    Jm2 = sum(G @ G for G in Jm)
    tot = {"L2": 0.0, "Jp2": 0.0, "Jm2": 0.0}
    for n in range(1, NMAX + 1):
        ids = [IDX[s] for s in STATES if s[0] == n]
        tL, tp, tm = (float(np.real(np.trace(X[np.ix_(ids, ids)]))) for X in (L2, Jp2, Jm2))
        tot["L2"] += tL
        tot["Jp2"] += tp
        tot["Jm2"] += tm
        cas_dev = max(np.abs(Jp2[np.ix_(ids, ids)] - (n * n - 1) / 4 * np.identity(n * n)).max(),
                      np.abs(Jm2[np.ix_(ids, ids)] - (n * n - 1) / 4 * np.identity(n * n)).max())
        print(f"F1 shell n={n}: Tr L^2 = {tL:.12f}, Tr J+^2 = {tp:.12f}, Tr J-^2 = {tm:.12f}, "
              f"ratio J+/L = {tp / tL if tL else float('nan'):.12f}; J+-^2 - (n^2-1)/4 max {cas_dev:.1e}")
    print(f"F1 totals n<={NMAX}: B = Tr L^2 = {tot['L2']:.10f}; Tr J+^2 = {tot['Jp2']:.10f}; Tr J-^2 = {tot['Jm2']:.10f}")
    nn = sp.Symbol("n", positive=True, integer=True)
    ll = sp.Symbol("l", integer=True, nonnegative=True)
    shell_B = sp.summation((2 * ll + 1) * ll * (ll + 1), (ll, 0, nn - 1))
    print(f"F1 symbolic: sum_l (2l+1) l(l+1) - 2 * n^2 (n^2-1)/4 = {sp.simplify(shell_B - nn**2 * (nn**2 - 1) / 2)}")

    Pm = chiral_projectors(-1)
    Pp = chiral_projectors(+1)
    groups, M = hopf_M_exact()
    d1 = d2 = d3 = d4 = 0.0
    for g, lab in enumerate(groups):
        n, mu = lab
        occ_num = np.real(np.diag(Pm[lab]))
        occ_ex = np.array([float(M[u, g]) for u in range(DIM)])
        d1 = max(d1, np.abs(occ_num - occ_ex).max())
        d2 = max(d2, np.abs(np.real(np.diag(Pp[lab])) - occ_num).max())
        # per-group L^2 trace is uniform: n (n^2 - 1)/2
        d3 = max(d3, abs(float(np.real(np.trace(Pm[lab] @ L2))) - n * (n * n - 1) / 2))
        # mu -> -mu with m -> -m
        mir = groups.index((n, -mu))
        for (nn_, l, m) in STATES:
            if nn_ == n:
                d4 = max(d4, abs(float(M[IDX[(n, l, m)], g] - M[IDX[(n, l, -m)], mir])))
        assert abs(np.trace(Pm[lab]).real - n) < 1e-10
    print(f"F2 exact CG^2 occupations vs numerical J-_z eigenprojector diagonals: max dev {d1:.1e}")
    print(f"F3 m+ grouping vs m- grouping, diagonal occupations: max dev {d2:.1e}  ((-1)^l maps one to the other)")
    print(f"F4 per-group Tr(P_g L^2) - n(n^2-1)/2: max dev {d3:.1e}  (uniform within a shell)")
    print(f"F5 mu -> -mu, m -> -m occupation symmetry: max dev {d4:.1e}")
    rows = [sum(M[u, g] for g in range(len(groups))) for u in range(DIM)]
    cols = [sum(M[u, g] for u in range(DIM)) for g in range(len(groups))]
    print(f"F6 membership rows sum to {set(rows)}; group weights {[str(c) for c in cols]} "
          f"(rotation: {[2 * int(l) + 1 for (_, l) in rotation_M()[0]]})")
    sn = [(n, [str(M[IDX[(n, 0, 0)], g]) for g, (ng, _) in enumerate(groups) if ng == n]) for n in range(1, NMAX + 1)]
    print(f"F7 s-node occupation per shell-n group: {sn}  (expected 1/n each)")
    for n in range(1, NMAX + 1):
        assert sum(c for c, (ng, _) in zip(cols, groups) if ng == n) == n * n
    print("F8 per-shell group weights sum to n^2 under both groupings, so F = sum n^2/n^4 is grouping-blind")
    Ag = graph()
    print(f"F9 1^T A 1 = {int(Ag.sum())} = 2 N_1; equals sum of W_full = M^T A M under ANY grouping with unit rows")


def controls() -> None:
    Ag = graph()
    groups, M = rotation_M()
    W = base_occupation(M, Ag)
    Lp = laplacian(W)
    x = sp.Symbol("x")
    print(f"C1 rotation sectors {[(n, int(l)) for (n, l) in groups]}, weights {[2 * int(l) + 1 for (_, l) in groups]}")
    nz = [((groups[i][0], int(groups[i][1])), (groups[j][0], int(groups[j][1])), W[i, j])
          for i in range(len(groups)) for j in range(i + 1, len(groups)) if W[i, j] != 0]
    print(f"C2 occupation rule, rotation: links {nz}")
    print(f"C3 occupation rule, rotation: charpoly {sp.factor(Lp.charpoly(x).as_expr())}, eigenvalues {Lp.eigenvals()}")
    P = coordinate_projectors(groups, M)
    Wf = base_frobenius(P, groups, Ag)
    print(f"C4 Frobenius rule, rotation: max |W_frob - W_occ| = "
          f"{np.abs(Wf - np.array(W.tolist(), dtype=float)).max():.1e}")
    # basis-independence control for the Frobenius rule: rotate each sector's basis by a random unitary
    rng = np.random.default_rng(7)
    U = np.zeros((DIM, DIM), complex)
    for lab in groups:
        ids = [u for u in range(DIM) if P[lab][u, u] == 1]
        Z = rng.normal(size=(len(ids), len(ids))) + 1j * rng.normal(size=(len(ids), len(ids)))
        Q, _ = np.linalg.qr(Z)
        U[np.ix_(ids, ids)] = Q
    Prot = {lab: U @ P[lab] @ U.conj().T for lab in groups}
    print(f"C5 Frobenius rule after a random unitary inside each sector (same subspaces): max change "
          f"{np.abs(base_frobenius(Prot, groups, Ag) - Wf).max():.1e}")
    # zero-mode counter on the rotation base
    ev = np.linalg.eigvalsh(np.array(Lp.tolist(), dtype=float))
    print(f"C6 zero-mode counter on rotation base: {int(np.sum(np.abs(ev) < 1e-9))} (Paper 25: 3)")


def score() -> None:  # NOT RUN before registration
    Ag = graph()
    groups, M = hopf_M_exact()
    W = base_occupation(M, Ag)
    Lp = laplacian(W)
    x = sp.Symbol("x")
    cp = Lp.charpoly(x).as_expr()
    print("HG groups", groups)
    print("HG W (exact)", W)
    print("HG charpoly", sp.factor(cp))
    ev = np.linalg.eigvalsh(np.array(Lp.tolist(), dtype=float))
    print("HG eigenvalues", ev, "zero modes", int(np.sum(np.abs(ev) < 1e-9)))
    inter = sum(W[i, j] for i in range(len(groups)) for j in range(i + 1, len(groups)))
    print("HG inter-group weight", inter, "intra (fibre) weight", Fraction(int(Ag.sum()) // 2) - Fraction(str(inter)))
    Wf = base_frobenius(chiral_projectors(-1), groups, Ag)
    print("HG-alt Frobenius W", np.round(Wf, 12))
    print("HG-alt eigenvalues", np.linalg.eigvalsh(np.diag(Wf.sum(1)) - Wf))


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "forced"
    {"forced": forced, "controls": controls, "score": score}[mode]()
