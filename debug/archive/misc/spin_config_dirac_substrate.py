"""
Spin-configuration / symmetry decomposition on the DIRAC substrate.
===================================================================

Sequel to debug/spin_config_symmetry_decomposition.py (scalar Fock graph).
Question: when spin is coupled in (the Camporesi-Higuchi Dirac-S^3 graph,
DiracLattice, E1-dipole Rule B), do the symmetric spin structures survive,
and does the genuine-spin (chirality) reflection still act?

Two reflections tested:
  sigma_mj : (n, kappa, m_j) -> (n, kappa, -m_j)   [MAGNETIC; relativistic
                                                    cousin of m -> -m]
  sigma_chi: chirality flip chi -> -chi (kappa -> partner)  [SPIN reflection;
                                                    the one the kappa/m_j-parity
                                                    sprints found resists]

Flat Z2 spin structures preserving coupling = H^1(G; F2) = ker(B mod 2),
dim beta_1.  A reflection preserves a spin structure iff it fixes it in the
cycle space; #symmetric = 2^dim(fixed subspace).
"""
from __future__ import annotations

import numpy as np
from geovac.dirac_lattice import DiracLattice


# --------------------------- F2 linear algebra ----------------------------- #
def f2_rref(M):
    R = (np.array(M, dtype=np.int64) % 2).copy()
    rows, cols = R.shape
    pivots, r = [], 0
    for c in range(cols):
        piv = next((rr for rr in range(r, rows) if R[rr, c]), None)
        if piv is None:
            continue
        R[[r, piv]] = R[[piv, r]]
        for rr in range(rows):
            if rr != r and R[rr, c]:
                R[rr] = (R[rr] + R[r]) % 2
        pivots.append(c); r += 1
        if r == rows:
            break
    return R, pivots


def f2_rank(M):
    return 0 if M.size == 0 else len(f2_rref(M)[1])


def f2_nullspace(M):
    M = np.array(M, dtype=np.int64) % 2
    rows, cols = M.shape
    R, piv = f2_rref(M)
    piv_set = set(piv)
    free = [c for c in range(cols) if c not in piv_set]
    basis = []
    for f in free:
        x = np.zeros(cols, dtype=np.int64); x[f] = 1
        for i, pc in enumerate(piv):
            x[pc] = R[i, f] % 2
        basis.append(x % 2)
    return (np.array(basis, dtype=np.int64).T if basis
            else np.zeros((cols, 0), dtype=np.int64))


def f2_independent_rows(C):
    """Return r row indices of C (E x r, rank r) forming an invertible r x r block."""
    Ct = (np.array(C, dtype=np.int64) % 2).T
    _, piv = f2_rref(Ct)        # pivot cols of C^T = independent rows of C
    return piv


def f2_inverse(M):
    M = np.array(M, dtype=np.int64) % 2
    n = M.shape[0]
    aug = np.concatenate([M, np.eye(n, dtype=np.int64)], axis=1)
    R, piv = f2_rref(aug)
    return R[:, n:] % 2


def f2_matmul(A, B):
    return (np.array(A, dtype=np.int64) @ np.array(B, dtype=np.int64)) % 2


# ------------------------------ analysis ----------------------------------- #
def analyze(n_max, mode='atomic'):
    dl = DiracLattice(n_max, mode=mode)
    labels = [(lab.n_fock, lab.kappa, lab.two_m_j) for lab in dl.labels]
    V = len(labels)
    idx = {lab: i for i, lab in enumerate(labels)}
    A = dl.adjacency.toarray()
    edges = [(i, j) for i in range(V) for j in range(i + 1, V) if A[i, j]]
    E = len(edges)
    edge_idx = {e: k for k, e in enumerate(edges)}

    # components / betti
    # beta0 via connected components of A
    seen = set(); comps = 0
    for s in range(V):
        if s in seen:
            continue
        comps += 1; stack = [s]
        while stack:
            u = stack.pop()
            if u in seen:
                continue
            seen.add(u)
            stack.extend(int(v) for v in np.nonzero(A[u])[0] if v not in seen)
    b0 = comps
    b1 = E - V + b0

    # cycle space = ker(B mod 2)
    B = np.zeros((V, E), dtype=np.int64)
    for k, (i, j) in enumerate(edges):
        B[i, k] = 1; B[j, k] = 1
    C = f2_nullspace(B)          # E x beta1
    r = C.shape[1]
    assert r == b1, (r, b1)

    # precompute for fast action matrices
    if r > 0:
        rows = f2_independent_rows(C)
        Cinv = f2_inverse(C[rows, :])
    else:
        rows, Cinv = [], None

    def action(P):
        """F2 action matrix on cycle space for node-permutation P (assumed automorphism)."""
        q = np.empty(E, dtype=np.int64)
        for k, (i, j) in enumerate(edges):
            a, c = int(P[i]), int(P[j])
            e = (a, c) if a < c else (c, a)
            q[k] = edge_idx[e]
        QC = np.zeros((E, r), dtype=np.int64)
        QC[q] = C                      # (QC)[q[e]] = C[e]
        return f2_matmul(Cinv, QC[rows, :])

    def is_auto(P):
        for (i, j) in edges:
            a, c = int(P[i]), int(P[j])
            e = (a, c) if a < c else (c, a)
            if e not in edge_idx:
                return False
        return True

    # sigma_mj : (n,kappa,2mj) -> (n,kappa,-2mj)
    P_mj = np.array([idx[(n, k, -tm)] for (n, k, tm) in labels], dtype=np.int64)

    res = {"n_max": n_max, "mode": mode, "V": V, "E": E,
           "beta0": b0, "beta1": b1, "configs": 2 ** b1}

    auto_mj = is_auto(P_mj)
    rec = {"is_automorphism": auto_mj, "fixed_point_free_on_nodes":
           all(P_mj[i] != i for i in range(V))}
    if auto_mj and r > 0:
        M = action(P_mj)
        dfix = r - f2_rank((M + np.eye(r, dtype=np.int64)) % 2)
        rec.update(symmetric_configs=2 ** dfix, fixed_dim=int(dfix))
    res["sigma_mj"] = rec

    # chirality structure: try to build chi-flip as a node permutation
    # chi = +1 for kappa<0 (j=l+1/2), chi=-1 for kappa>0 (j=l-1/2)
    n_pos = sum(1 for (_, k, _) in labels if k < 0)   # chi=+1
    n_neg = sum(1 for (_, k, _) in labels if k > 0)   # chi=-1
    res["chirality_counts"] = {"chi_plus(kappa<0)": n_pos, "chi_minus(kappa>0)": n_neg,
                               "balanced": n_pos == n_neg}
    # per-(n,l) block sizes to expose the 2l+2 vs 2l mismatch
    from collections import defaultdict
    blk = defaultdict(lambda: {"neg": 0, "pos": 0})
    for (n, k, tm) in labels:
        l = abs(k) - 1 if k < 0 else k
        blk[(n, l)]["neg" if k < 0 else "pos"] += 1
    res["nl_block_sizes"] = {f"n={n},l={l}": dict(v) for (n, l), v in sorted(blk.items())}
    return res


def pretty(res):
    print(f"\n{'='*72}")
    print(f"  DIRAC  n_max={res['n_max']} ({res['mode']})   V={res['V']} E={res['E']}"
          f"  beta0={res['beta0']} beta1={res['beta1']}")
    print(f"  coupling-preserving Z2 spin structures = 2^beta1 = {res['configs']}")
    print(f"{'-'*72}")
    s = res["sigma_mj"]
    print(f"  sigma_mj  (m_j -> -m_j, MAGNETIC reflection)")
    print(f"      automorphism            : {s['is_automorphism']}")
    print(f"      fixed-point-free nodes  : {s['fixed_point_free_on_nodes']}  (Kramers)")
    if "symmetric_configs" in s:
        print(f"      symmetric spin structs  : {s['symmetric_configs']}  "
              f"(fixed dim {s['fixed_dim']})")
    cc = res["chirality_counts"]
    print(f"  sigma_chi (chi -> -chi, SPIN reflection)")
    print(f"      chi=+1 (kappa<0) nodes  : {cc['chi_plus(kappa<0)']}")
    print(f"      chi=-1 (kappa>0) nodes  : {cc['chi_minus(kappa>0)']}")
    print(f"      balanced?               : {cc['balanced']}  "
          f"-> {'reflection POSSIBLE' if cc['balanced'] else 'NO node bijection: SPIN reflection OBSTRUCTED'}")
    print(f"      per-(n,l) block sizes (neg=chi+, pos=chi-):")
    for key, v in res["nl_block_sizes"].items():
        print(f"        {key:10s}  chi+={v['neg']:2d}  chi-={v['pos']:2d}"
              f"   {'<-- mismatch '+str(v['neg'])+' vs '+str(v['pos']) if v['neg']!=v['pos'] else ''}")


if __name__ == "__main__":
    for nm in (2, 3, 4):
        pretty(analyze(nm, mode='atomic'))
