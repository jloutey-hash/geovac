"""Review check of Lemma 1 on trees the memo could not finish: least-order a with a H pure, modulo a prime.

Reviewer's own implementation of the module model (rules read from the memo text):
  theta z = 2z; theta K0 = D; theta D = z K0; theta P_y = P'_y; theta P'_y = -P'_y - z K0 prod_{c child of y} P_c.
Finds the least theta-order N (rank at a random point), then the least z-degree d with a polynomial solution
a = sum_{j<=N, k<=d} c_jk z^k theta^j, a H in the pure top, by a linear solve over F_p (random projection of the
non-pure coordinates, evaluation at random points).  Reports the weighted symbol at (1, -t) mod p: its value at
t = -(n-1), and whether it splits into roots among 0, +-1, ..., +-(n-2).
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
from flint import nmod_mat, nmod_poly

P = 2147483629
ROOT = Path(__file__).resolve().parent.parent
T0 = time.time()
rng = np.random.default_rng(20261118)


def log(*a) -> None:
    print(f"[{time.time() - T0:7.1f}s]", *a, flush=True)


def make_tree(n: int, edges: List[Tuple[int, int]], root: int):
    adj = {v: [] for v in range(n)}
    for a, b in edges:
        adj[a].append(b)
        adj[b].append(a)
    children = {}

    def rec(v, p):
        children[v] = [w for w in adj[v] if w != p]
        for w in children[v]:
            rec(w, v)
    rec(root, -1)
    return children


def padd(d: Dict, key, poly: List[int]) -> None:
    cur = d.get(key)
    if cur is None:
        d[key] = list(poly)
        return
    if len(cur) < len(poly):
        cur.extend([0] * (len(poly) - len(cur)))
    for i, c in enumerate(poly):
        cur[i] = (cur[i] + c) % P


def theta_vec(children, vec: Dict) -> Dict:
    out: Dict = {}
    for (S, kept, i), c in vec.items():
        k = len(S)
        padd(out, (S, kept, i), [(2 * m * cm) % P for m, cm in enumerate(c)])
        if k - i > 0:
            padd(out, (S, kept, i + 1), [((k - i) * cm) % P for cm in c])
        if i > 0:
            padd(out, (S, kept, i - 1), [0] + [(i * cm) % P for cm in c])
        for idx, (yv, st) in enumerate(kept):
            rest = kept[:idx] + kept[idx + 1:]
            if st == 0:
                padd(out, (S, tuple(sorted(rest + ((yv, 1),))), i), c)
            else:
                padd(out, (S, kept, i), [(-cm) % P for cm in c])
                newkept = tuple(sorted(rest + tuple((w, 0) for w in children[yv])))
                padd(out, (S | frozenset([yv]), newkept, i), [0] + [(-cm) % P for cm in c])
    res = {}
    for k_, v in out.items():
        while v and v[-1] == 0:
            v.pop()
        if v:
            res[k_] = v
    return res


def horner(poly: List[int], zs: np.ndarray) -> np.ndarray:
    acc = np.zeros_like(zs)
    for c in reversed(poly):
        acc = (acc * zs + c) % P
    return acc


def analyse(name: str, n: int, edges, root: int, max_order: int = 80) -> Dict:
    children = make_tree(n, edges, root)
    v = {(frozenset(), tuple(sorted((w, 0) for w in children[root])), 0): [0, 1]}
    vs = [v]
    pure = lambda key: len(key[0]) == n - 1  # noqa: E731
    z0 = np.array([int(rng.integers(2, P - 1))], dtype=np.int64)
    N = None
    for j in range(1, max_order + 1):
        vs.append(theta_vec(children, vs[-1]))
        cols = sorted({k_ for u in vs for k_ in u if not pure(k_)}, key=str)
        M = nmod_mat(len(vs), len(cols), [int(horner(u.get(k_, [0]), z0)[0]) for u in vs for k_ in cols], P)
        if M.rank() < len(vs):
            N = j
            break
    if N is None:
        return {"label": name, "status": "no relation up to order %d" % max_order}
    ncols = len(cols)
    log(f"{name}: least theta-order N = {N}; non-pure coordinates {ncols}")
    r = N + 3
    maxdeg = max(len(pl) for u in vs for pl in u.values())
    Rm = rng.integers(0, 1 << 20, size=(ncols, r), dtype=np.int64)

    def solve(d: int):
        nunk = (N + 1) * (d + 1)
        m = max(nunk // r + 6, d + maxdeg + 6)
        zs = rng.integers(2, P - 1, size=m, dtype=np.int64)
        V = np.zeros((m, N + 1, ncols), dtype=np.int64)
        for j, u in enumerate(vs):
            for ci, k_ in enumerate(cols):
                if k_ in u:
                    V[:, j, ci] = horner(u[k_], zs)
        PV = np.matmul(V, Rm) % P                               # (m, N+1, r); entries < 2^31 * 2^20 * ncols
        zp = np.ones((m, d + 1), dtype=np.int64)
        for k in range(1, d + 1):
            zp[:, k] = (zp[:, k - 1] * zs) % P
        # rows (i, rr); columns (j, k)
        E = (PV.transpose(0, 2, 1)[:, :, :, None] * zp[:, None, None, :]) % P   # (m, r, N+1, d+1)
        E = E.reshape(m * r, nunk)
        Mx = nmod_mat(m * r, nunk, E.ravel().tolist(), P)
        X, nullity = Mx.nullspace()
        return X, nullity, nunk

    d_try = 40
    while True:
        X, nullity, nunk = solve(d_try)
        log(f"{name}: z-degree bound {d_try}: nullity {nullity}")
        if nullity >= 1:
            break
        d_try = int(d_try * 1.6) + 5
        if (N + 1) * (d_try + 1) > 9000:
            return {"label": name, "theta_order": N, "status": "z-degree above cap (%d)" % d_try}
    d_min = d_try - (nullity - 1)
    X, nullity, nunk = solve(d_min)
    log(f"{name}: least z-degree {d_min}: nullity {nullity}")
    assert nullity == 1, nullity
    a = {}
    for j in range(N + 1):
        for k in range(d_min + 1):
            c = int(X[j * (d_min + 1) + k, 0])
            if c:
                a[(k, j)] = c
    # verify a H pure on the FULL (unprojected) coordinates, exactly mod p as polynomials
    tot: Dict = {}
    for (k, j), c in a.items():
        for key, poly in vs[j].items():
            padd(tot, key, [0] * k + [(c * cm) % P for cm in poly])
    nonpure_left = sum(1 for key, poly in tot.items() if not pure(key) and any(poly))
    pure_nonzero = sum(1 for key, poly in tot.items() if pure(key) and any(poly))
    w = max(2 * k + j for (k, j) in a)
    lead = [0] * (N + 1)
    for (k, j), c in a.items():
        if 2 * k + j == w:
            lead[j] = (lead[j] + c * (-1) ** j) % P
    lp = nmod_poly(lead, P)
    p = -(n - 1)
    val_p = int(lp(p % P))
    deg = lp.degree()
    rem = lp
    roots = {}
    for rt in range(-(n + 1), n + 2):
        lin = nmod_poly([(-rt) % P, 1], P)
        mult = 0
        while rem.degree() > 0:
            q, rr = divmod(rem, lin)
            if rr.degree() >= 0 and not rr.is_zero():
                break
            rem = q
            mult += 1
        if mult:
            roots[rt] = mult
    res = {"label": name, "n": n, "root": root, "edges": edges, "prime": P, "theta_order": N, "z_degree": d_min,
           "weight": w, "symbol_degree_in_t": deg, "symbol_at_p_mod_prime": val_p, "nonzero_at_p": bool(val_p != 0),
           "roots_found": {str(k_): v_ for k_, v_ in roots.items()}, "splits_completely": bool(rem.degree() == 0),
           "max_abs_root": max((abs(k_) for k_ in roots), default=0),
           "roots_within_n_minus_2": bool(all(abs(k_) <= n - 2 for k_ in roots)),
           "nonpure_coordinates_left_in_aH": nonpure_left, "pure_coordinates_nonzero_in_aH": pure_nonzero}
    log(f"{name}: weight {w}; symbol degree {deg}; roots {roots}; splits {res['splits_completely']}; "
        f"nonzero at t={p}: {res['nonzero_at_p']}; a H non-pure left {nonpure_left}, pure nonzero {pure_nonzero}")
    return res


TREES = {
    "star3_centre": (3, [(0, 1), (0, 2)], 0),                      # control: memo -4 t^3 (t-1)^2 (t+1)^2, weight 13
    "path4_end": (4, [(0, 1), (1, 2), (2, 3)], 0),                 # control: memo weight 32
    "star4": (4, [(0, 1), (0, 2), (0, 3)], 0),                     # control: memo weight 56
    "path5_end": (5, [(0, 1), (1, 2), (2, 3), (3, 4)], 0),         # control: memo weight 100
    "fork5_root0": (5, [(0, 1), (1, 2), (1, 3), (3, 4)], 0),       # branched, not in the memo's exact table
    "fork5_root4": (5, [(0, 1), (1, 2), (1, 3), (3, 4)], 4),
    "fork5_root2": (5, [(0, 1), (1, 2), (1, 3), (3, 4)], 2),
    "path6_end": (6, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5)], 0),
    "broom6": (6, [(0, 1), (1, 2), (2, 3), (3, 4), (3, 5)], 0),    # branched, n = 6 (even, logarithmic case)
    "path5_centre": (5, [(0, 1), (1, 2), (0, 3), (3, 4)], 0),
    "fork5_root_branch": (5, [(0, 1), (1, 2), (1, 3), (3, 4)], 1),
}

if __name__ == "__main__":
    which = sys.argv[1:] or list(TREES)
    outp = ROOT / "debug" / "data" / "tree_attain_review_modp.json"
    data = json.loads(outp.read_text()) if outp.exists() else {}
    for name in which:
        n, edges, root = TREES[name]
        t1 = time.time()
        try:
            res = analyse(name, n, edges, root)
        except MemoryError as e:
            res = {"label": name, "status": "MemoryError %s" % e}
        res["seconds"] = round(time.time() - t1, 1)
        data[name] = res
        outp.write_text(json.dumps(data, indent=1))
    log("done")
