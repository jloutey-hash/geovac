"""
Spin-configuration / symmetry-decomposition probe on the Layer-1 Fock graph.
============================================================================

Question (Josh, 2026-06-08): on the bare Layer-1 substrate, which arrangements
preserve coupling (adjacency), and how do the graph symmetries decompose across
them?  Submission: different projections can carry different spin configurations.

Formalization.
  - "Arrangements preserving coupling" at the Z2 level = flat Z2 connections on
    the graph = the cycle space H^1(G; F2) = ker(B mod 2), dim = beta_1.
    (B = signed incidence matrix.)  There are 2^beta_1 of them.
  - A graph symmetry sigma (node permutation preserving adjacency) acts on this
    cycle space as an F2-linear map M_sigma.
  - Configs that ALSO preserve the symmetry sigma = fixed subspace ker(M_sigma - I),
    dim d_fix -> 2^d_fix symmetric configs.  #orbits under <sigma> = (2^b1 + 2^d_fix)/2.

Symmetries tested:
  sigma_m : (n,l,m) -> (n,l,-m)              [Hopf Z2 / m-reflection; PHYSICAL]
  sigma_n : (n,l,m) -> (n_max+l+1-n, l, m)   [radial inversion; per-l grid
                                              reflection; HIDDEN, not physical]

Everything is exact F2 arithmetic.
"""
from __future__ import annotations

import numpy as np

from geovac.fock_graph_hodge import FockGraphHodge


# --------------------------------------------------------------------------- #
# F2 linear algebra (numpy int arrays mod 2)
# --------------------------------------------------------------------------- #
def f2_rref(M):
    """Reduced row echelon form over F2; returns (R, pivot_cols)."""
    R = (np.array(M, dtype=np.int64) % 2).copy()
    rows, cols = R.shape
    pivots = []
    r = 0
    for c in range(cols):
        piv = None
        for rr in range(r, rows):
            if R[rr, c]:
                piv = rr
                break
        if piv is None:
            continue
        R[[r, piv]] = R[[piv, r]]
        for rr in range(rows):
            if rr != r and R[rr, c]:
                R[rr] = (R[rr] + R[r]) % 2
        pivots.append(c)
        r += 1
        if r == rows:
            break
    return R, pivots


def f2_rank(M):
    if M.size == 0:
        return 0
    _, piv = f2_rref(M)
    return len(piv)


def f2_nullspace(M):
    """Basis of {x : M x = 0 mod 2} as columns of an array (cols x dim)."""
    M = np.array(M, dtype=np.int64) % 2
    rows, cols = M.shape
    R, piv = f2_rref(M)
    piv_set = set(piv)
    free = [c for c in range(cols) if c not in piv_set]
    basis = []
    for f in free:
        x = np.zeros(cols, dtype=np.int64)
        x[f] = 1
        for i, pc in enumerate(piv):
            x[pc] = R[i, f] % 2
        basis.append(x % 2)
    if not basis:
        return np.zeros((cols, 0), dtype=np.int64)
    return np.array(basis, dtype=np.int64).T


def f2_solve_in_basis(C, b):
    """Solve C a = b over F2 where C (E x r) has full column rank. Returns a (r,)."""
    E, r = C.shape
    aug = np.concatenate([C % 2, (b % 2).reshape(-1, 1)], axis=1)
    R, piv = f2_rref(aug)
    a = np.zeros(r, dtype=np.int64)
    for i, pc in enumerate(piv):
        if pc < r:
            a[pc] = R[i, -1] % 2
    return a % 2


# --------------------------------------------------------------------------- #
# Build graph, cycle space, symmetry actions
# --------------------------------------------------------------------------- #
def analyze(n_max: int):
    fg = FockGraphHodge(n_max)
    states = list(fg.states)              # [(n,l,m), ...]
    edges = list(fg.edges)                # [(i,j), ...]  i<j
    V, E = fg.n_nodes, fg.n_edges
    b0, b1 = fg.betti_0, fg.betti_1

    idx = {s: k for k, s in enumerate(states)}
    edge_idx = {e: k for k, e in enumerate(edges)}

    # Signed incidence B (V x E); cycle space = ker(B mod 2)
    B = np.zeros((V, E), dtype=np.int64)
    for k, (i, j) in enumerate(edges):
        B[i, k] = 1
        B[j, k] = 1  # mod 2, sign irrelevant
    cycles = f2_nullspace(B)              # E x beta1, columns are cycle vectors
    beta1_check = cycles.shape[1]

    # --- node permutations ---
    def perm_from_map(fn):
        P = np.empty(V, dtype=np.int64)
        for s, k in idx.items():
            P[k] = idx[fn(s)]
        return P

    sigma_m = perm_from_map(lambda s: (s[0], s[1], -s[2]))
    sigma_n = perm_from_map(lambda s: (n_max + s[1] + 1 - s[0], s[1], s[2]))

    def is_automorphism(P):
        for (i, j) in edges:
            a, c = int(P[i]), int(P[j])
            e = (a, c) if a < c else (c, a)
            if e not in edge_idx:
                return False
        return True

    def edge_perm(P):
        """Induced permutation on edge indices (assumes P is an automorphism)."""
        q = np.empty(E, dtype=np.int64)
        for k, (i, j) in enumerate(edges):
            a, c = int(P[i]), int(P[j])
            e = (a, c) if a < c else (c, a)
            q[k] = edge_idx[e]
        return q

    def action_on_cycles(P):
        """F2 matrix M (beta1 x beta1): sigma maps cycle basis col k to C @ M[:,k]."""
        if beta1_check == 0:
            return np.zeros((0, 0), dtype=np.int64)
        q = edge_perm(P)
        M = np.zeros((beta1_check, beta1_check), dtype=np.int64)
        for k in range(beta1_check):
            c = cycles[:, k]
            cq = np.zeros(E, dtype=np.int64)
            cq[q] = c                      # permute edges: (Qc)[q[e]] = c[e]
            a = f2_solve_in_basis(cycles, cq)
            M[:, k] = a
        return M

    out = {
        "n_max": n_max, "V": V, "E": E, "beta0": b0, "beta1": b1,
        "beta1_cyclecheck": beta1_check,
        "coupling_preserving_configs": 2 ** b1,
        "symmetries": {},
    }

    for name, P in [("sigma_m (m->-m, Hopf, physical)", sigma_m),
                    ("sigma_n (radial inversion, hidden)", sigma_n)]:
        auto = is_automorphism(P)
        rec = {"is_automorphism": auto}
        if auto and beta1_check > 0:
            M = action_on_cycles(P)
            d_fix = beta1_check - f2_rank((M + np.eye(beta1_check, dtype=np.int64)) % 2)
            n_sym = 2 ** d_fix
            n_orbits = (2 ** b1 + n_sym) // 2
            rec.update({
                "fixed_subspace_dim": int(d_fix),
                "symmetry_preserving_configs": int(n_sym),
                "broken_symmetry_configs": int(2 ** b1 - n_sym),
                "num_orbits": int(n_orbits),
                "action_matrix_F2": M.tolist(),
            })
        elif auto:
            rec.update({"fixed_subspace_dim": 0, "symmetry_preserving_configs": 1,
                        "broken_symmetry_configs": 0, "num_orbits": 1})
        out["symmetries"][name] = rec

    return out, cycles, states, edges


def pretty(out):
    print(f"\n{'='*70}")
    print(f"  n_max = {out['n_max']}   V={out['V']}  E={out['E']}  "
          f"beta0={out['beta0']}  beta1={out['beta1']}")
    print(f"  coupling-preserving Z2 spin configs = 2^beta1 = "
          f"{out['coupling_preserving_configs']}")
    print(f"{'-'*70}")
    for name, rec in out["symmetries"].items():
        print(f"  {name}")
        if not rec["is_automorphism"]:
            print(f"      NOT a graph automorphism")
            continue
        print(f"      automorphism: YES")
        print(f"      fixed-subspace dim       : {rec['fixed_subspace_dim']}")
        print(f"      symmetry-preserving cfgs : {rec['symmetry_preserving_configs']}")
        print(f"      symmetry-breaking  cfgs  : {rec['broken_symmetry_configs']}")
        print(f"      # orbits of <sigma>      : {rec['num_orbits']}")


if __name__ == "__main__":
    for nm in (2, 3, 4):
        out, cycles, states, edges = analyze(nm)
        pretty(out)

    # Explicit enumeration of the 4 sectors at n_max=3 with Hopf action
    print(f"\n{'='*70}")
    print("  n_max=3 explicit: 4 sectors, plaquette occupation under Hopf m->-m")
    print(f"{'-'*70}")
    out, cycles, states, edges = analyze(3)
    b1 = out["beta1"]
    # cycle basis columns; identify which (n,l,m) plaquette each is by its edges
    for k in range(cycles.shape[1]):
        es = [edges[e] for e in range(len(edges)) if cycles[e, k]]
        nodes = sorted(set([states[i] for e in es for i in e]))
        print(f"   cycle {k}: plaquette nodes {nodes}")
    Mm = np.array(out["symmetries"]["sigma_m (m->-m, Hopf, physical)"]["action_matrix_F2"])
    print(f"\n   Hopf action on cycle basis (F2):\n{Mm}")
    print("\n   sectors (c0,c1):  00=trivial  11=both-twisted  10/01=single-plaquette")
    for c in [(0, 0), (1, 0), (0, 1), (1, 1)]:
        v = np.array(c) % 2
        img = (Mm @ v) % 2
        fixed = np.array_equal(img, v)
        print(f"     {c} -> {tuple(int(x) for x in img)}   "
              f"{'SYMMETRIC (preserves Hopf)' if fixed else 'broken -> orbit partner'}")
