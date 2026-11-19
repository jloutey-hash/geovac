"""tree_attain_l, step (i): for Coulomb trees with l > 0 densities, the least-order a = sum c_jk z^k theta^j with
a G in the pure top (n - 1 Bessel factors, no potential left), modulo a prime; its weighted symbol at (1, -t).

G = x^2 <P_0(-d) prod Phi_w> is the single function whose transform is the fibre (root derivative moved by parts).
Model and rules: debug/tree_attain_l_model.py.  Method as in debug/tree_attain_review_modp.py (rank at a random
point for the order, then a polynomial solve at random evaluation points with a random projection of rows).
Usage: python debug/tree_attain_l_algebra.py [tree names]
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path
from typing import Dict

import numpy as np
from flint import nmod_mat, nmod_poly

sys.path.insert(0, str(Path(__file__).resolve().parent))
import tree_attain_l_model as M  # noqa: E402

P = 2147483629
ROOT = Path(__file__).resolve().parent.parent
T0 = time.time()
rng = np.random.default_rng(20261119)


def log(*a) -> None:
    print(f"[{time.time() - T0:7.1f}s]", *a, flush=True)


def horner(poly, zs: np.ndarray) -> np.ndarray:
    acc = np.zeros_like(zs)
    for c in reversed(poly):
        acc = (acc * zs + int(c)) % P
    return acc


def analyse(name: str, n: int, edges, root: int, dens: Dict[int, str], max_order: int = 80) -> Dict:
    tree = M.LTree(n, edges, root, {v: M.Pdict(dens[v]) for v in range(n)})
    mdl = M.Model(tree, P)
    G = mdl.G()
    if not G:
        return {"label": name, "status": "G = 0 identically in the model"}
    vs = [G]
    pure = mdl.pure
    z0 = np.array([int(rng.integers(2, P - 1))], dtype=np.int64)
    N = None
    for j in range(1, max_order + 1):
        vs.append(mdl.theta(vs[-1]))
        cols = sorted({k_ for u in vs for k_ in u if not pure(k_)}, key=str)
        Mx = nmod_mat(len(vs), len(cols), [int(horner(u.get(k_, [0]), z0)[0]) for u in vs for k_ in cols], P)
        if Mx.rank() < len(vs):
            N = j
            break
        if time.time() - T0 > 1100:
            return {"label": name, "status": "time cap at order %d" % j}
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
        PV = np.zeros((m, N + 1, r), dtype=np.int64)
        for ci in range(ncols):                       # chunked to avoid int64 overflow
            PV = (PV + V[:, :, ci:ci + 1] * Rm[ci][None, None, :]) % P
        zp = np.ones((m, d + 1), dtype=np.int64)
        for k in range(1, d + 1):
            zp[:, k] = (zp[:, k - 1] * zs) % P
        E = (PV.transpose(0, 2, 1)[:, :, :, None] * zp[:, None, None, :]) % P
        E = E.reshape(m * r, nunk)
        Mx = nmod_mat(m * r, nunk, E.ravel().tolist(), P)
        X, nullity = Mx.nullspace()
        return X, nullity, nunk

    d_try = 12
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
    # a G on the FULL coordinates, exactly mod p
    tot: Dict = {}
    for (k, j), c in a.items():
        for key, poly in vs[j].items():
            mdl.padd(tot, key, poly, k, c)
    tot = mdl.trim(tot)
    nonpure_left = sum(1 for key in tot if not pure(key))
    pure_keys = [key for key in tot if pure(key)]
    pure_parities = sorted({key[3] for key in pure_keys})
    min_rpow = min((M.Model.rpow_key(key[3], next(m for m, c in enumerate(tot[key]) if c)) for key in pure_keys),
                   default=None)
    w = max(2 * k + j for (k, j) in a)
    lead = [0] * (N + 1)
    for (k, j), c in a.items():
        if 2 * k + j == w:
            lead[j] = (lead[j] + c * (-1) ** j) % P
    lp = nmod_poly(lead, P)
    p = -(n - 1)
    val_p = int(lp(p % P))
    val_mp = int(lp((n - 1) % P))
    rem = lp
    roots = {}
    for rt in range(-(n + 2), n + 3):
        lin = nmod_poly([(-rt) % P, 1], P)
        mult = 0
        while rem.degree() > 0:
            q, rr = divmod(rem, lin)
            if not rr.is_zero():
                break
            rem = q
            mult += 1
        if mult:
            roots[rt] = mult
    res = {"label": name, "n": n, "root": root, "edges": edges, "densities": dens, "prime": P,
           "G_coordinates": len(G), "theta_order": N, "z_degree": d_min, "weight": w,
           "symbol_degree_in_t": lp.degree(), "nonzero_at_t=-(n-1)": bool(val_p != 0),
           "nonzero_at_t=+(n-1)": bool(val_mp != 0),
           "roots_found": {str(k_): v_ for k_, v_ in roots.items()}, "splits_completely": bool(rem.degree() == 0),
           "max_abs_root": max((abs(k_) for k_ in roots), default=0),
           "roots_within_n_minus_2": bool(all(abs(k_) <= n - 2 for k_ in roots)),
           "nonpure_coordinates_left_in_aG": nonpure_left, "pure_coordinates_nonzero_in_aG": len(pure_keys),
           "pure_r_parities": pure_parities, "pure_min_r_power": min_rpow}
    log(f"{name}: weight {w}; roots {roots}; splits {res['splits_completely']}; nonzero at t={p}: "
        f"{res['nonzero_at_t=-(n-1)']}; aG non-pure left {nonpure_left}, pure nonzero {len(pure_keys)}, "
        f"parities {pure_parities}, min r-power {min_rpow}")
    return res


TREES = {
    # n = 3, the three shapes on record (root = 0)
    "A_pz_s_pz_rootmid": (3, [(0, 1), (0, 2)], 0, {0: "s", 1: "pz", 2: "pz"}),
    "B1_pz_pz_d_rootend": (3, [(0, 1), (1, 2)], 0, {0: "h2", 1: "pz", 2: "pz"}),
    "B2_pz_pz_cancel_rootend": (3, [(0, 1), (1, 2)], 0, {0: "cancel", 1: "pz", 2: "pz"}),
    # n = 3, an l = 1 root and an l = 2 leaf
    "C_pzroot_mid_pz_s": (3, [(0, 1), (0, 2)], 0, {0: "pz", 1: "pz", 2: "s"}),
    "D_d_leaf_d_leaf_sroot": (3, [(0, 1), (0, 2)], 0, {0: "s", 1: "h2", 2: "h2"}),
    "E_s_pz_pz_path_sroot_end": (3, [(0, 1), (1, 2)], 0, {0: "s", 1: "pz", 2: "pz"}),
    # controls (s-type, against the memo's table): star 3 weight 13, path 3 weight 8
    "ctrl_star3": (3, [(0, 1), (0, 2)], 0, {0: "s", 1: "s", 2: "s"}),
    "ctrl_path3": (3, [(0, 1), (1, 2)], 0, {0: "s", 1: "s", 2: "s"}),
    # n = 4
    "F_star4_sroot_pz_pz_s": (4, [(0, 1), (0, 2), (0, 3)], 0, {0: "s", 1: "pz", 2: "pz", 3: "s"}),
    "G_star4_pzroot_pz_s_s": (4, [(0, 1), (0, 2), (0, 3)], 0, {0: "pz", 1: "pz", 2: "s", 3: "s"}),
    "H_path4_sroot_end_internal_pz": (4, [(0, 1), (1, 2), (2, 3)], 0, {0: "s", 1: "pz", 2: "pz", 3: "s"}),
    "I_star4_sroot_d_d_s": (4, [(0, 1), (0, 2), (0, 3)], 0, {0: "s", 1: "h2", 2: "h2", 3: "s"}),
    "J_path4_droot_end": (4, [(0, 1), (1, 2), (2, 3)], 0, {0: "h2", 1: "pz", 2: "pz", 3: "s"}),
    "K_path4_dz2leaf_lap": (4, [(0, 1), (1, 2), (0, 3)], 0, {0: "lap", 1: "s", 2: "dz2", 3: "s"}),
    # n = 5
    "P_path5_droot_end": (5, [(0, 1), (1, 2), (2, 3), (3, 4)], 0, {0: "h2", 1: "pz", 2: "s", 3: "pz", 4: "s"}),
}

if __name__ == "__main__":
    which = sys.argv[1:] or list(TREES)
    outp = ROOT / "debug" / "data" / "tree_attain_l_algebra.json"
    data = json.loads(outp.read_text()) if outp.exists() else {}
    for name in which:
        n, edges, root, dens = TREES[name]
        t1 = time.time()
        try:
            res = analyse(name, n, edges, root, dens)
        except MemoryError as e:
            res = {"label": name, "status": "MemoryError %s" % e}
        res["seconds"] = round(time.time() - t1, 1)
        data[name] = res
        outp.write_text(json.dumps(data, indent=1))
    log("done")
