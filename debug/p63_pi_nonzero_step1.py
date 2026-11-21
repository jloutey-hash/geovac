"""Step (i) for non-homogeneous P_v: least-order a with a G pure, mod p; symbol roots; purity; parity; grading.
Uses tests/_tree_attain_l.py's model unchanged, with extra density entries added at runtime."""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction
from pathlib import Path

import numpy as np
from flint import nmod_mat, nmod_poly

ROOT = Path(r"C:\Users\jlout\Desktop\Project_Geometric")
sys.path.insert(0, str(ROOT / "tests"))
import _tree_attain_l as TL  # noqa: E402

TL.DENS.update({
    "one_p_dz": {(0, 0): Fraction(1), (1, 0): Fraction(1)},            # 1 + d_z (degrees 0, 1)
    "one_p_dz2": {(0, 0): Fraction(1), (2, 0): Fraction(1)},           # 1 + d_z^2 (degrees 0, 2)
    "mix012": {(0, 0): Fraction(1), (1, 0): Fraction(1), (2, 0): Fraction(1)},   # 1 + d_z + d_z^2
    "pz_p_lap": {(1, 0): Fraction(1), (0, 1): Fraction(1, 2)},         # d_z + lap/2 (degrees 1, 2)
})
P = TL.PRIME
rng = np.random.default_rng(20261121)
T0 = time.time()


def horner(poly, zs):
    acc = np.zeros_like(zs)
    for c in reversed(poly):
        acc = (acc * zs + int(c)) % P
    return acc


def analyse(name, n, edges, dens, max_order=60):
    tree = TL.LTree(n, edges, 0, dens)
    mdl = TL.Model(tree, P)
    G = mdl.G()
    grading = sorted({p for v in TL.theta_tower(mdl, G, 4) for p in TL.graded_degree_parities(v)})
    vs = [G]
    z0 = np.array([int(rng.integers(2, P - 1))], dtype=np.int64)
    N = None
    for j in range(1, max_order + 1):
        vs.append(mdl.theta(vs[-1]))
        cols = sorted({k for u in vs for k in u if not mdl.pure(k)}, key=str)
        Mx = nmod_mat(len(vs), len(cols), [int(horner(u.get(k, [0]), z0)[0]) for u in vs for k in cols], P)
        if Mx.rank() < len(vs):
            N = j
            break
    assert N is not None
    ncols = len(cols)
    r = N + 3
    maxdeg = max(len(pl) for u in vs for pl in u.values())
    Rm = rng.integers(0, 1 << 20, size=(ncols, r), dtype=np.int64)

    def solve(d):
        nunk = (N + 1) * (d + 1)
        m = max(nunk // r + 6, d + maxdeg + 6)
        zs = rng.integers(2, P - 1, size=m, dtype=np.int64)
        V = np.zeros((m, N + 1, ncols), dtype=np.int64)
        for j, u in enumerate(vs):
            for ci, k in enumerate(cols):
                if k in u:
                    V[:, j, ci] = horner(u[k], zs)
        PV = np.zeros((m, N + 1, r), dtype=np.int64)
        for ci in range(ncols):
            PV = (PV + V[:, :, ci:ci + 1] * Rm[ci][None, None, :]) % P
        zp = np.ones((m, d + 1), dtype=np.int64)
        for k in range(1, d + 1):
            zp[:, k] = (zp[:, k - 1] * zs) % P
        E = (PV.transpose(0, 2, 1)[:, :, :, None] * zp[:, None, None, :]) % P
        E = E.reshape(m * r, nunk)
        X, nul = nmod_mat(m * r, nunk, E.ravel().tolist(), P).nullspace()
        return X, nul

    d = 8
    while True:
        X, nul = solve(d)
        if nul >= 1:
            break
        d = int(d * 1.6) + 5
    dmin = d - (nul - 1)
    X, nul = solve(dmin)
    assert nul == 1, nul
    a = {}
    for j in range(N + 1):
        for k in range(dmin + 1):
            c = int(X[j * (dmin + 1) + k, 0])
            if c:
                a[(k, j)] = c
    aG = TL.apply_operator(mdl, a, G)
    nonpure, pure = TL.purity_report(mdl, aG)
    par = TL.pure_r_parities(mdl, aG)
    w = TL.weight(a)
    lead = [int(c) % P for c in TL.symbol_at_1_minus_t(a)]
    lp = nmod_poly(lead, P)
    rem, roots = lp, {}
    for rt in range(-(n + 2), n + 3):
        lin = nmod_poly([(-rt) % P, 1], P)
        m_ = 0
        while rem.degree() > 0:
            q, rr = divmod(rem, lin)
            if not rr.is_zero():
                break
            rem, m_ = q, m_ + 1
        if m_:
            roots[rt] = m_
    res = {"name": name, "n": n, "dens": dens, "Pi": str(TL.Pi(tree)), "principal_avg": str(TL.principal_average(tree)),
           "theta_order": N, "z_degree": dmin, "weight": w, "roots": roots, "splits": rem.degree() == 0,
           "nonzero_at_-(n-1)": int(lp((-(n - 1)) % P)) != 0, "aG_nonpure": nonpure, "aG_pure_coords": pure,
           "pure_parities": par, "grading_parities_theta_tower": grading, "secs": round(time.time() - T0, 1)}
    print(json.dumps(res, default=str), flush=True)
    return res


TREES = {
    "N1_sroot_star_1pdz_1pdz": (3, [(0, 1), (0, 2)], {0: "s", 1: "one_p_dz", 2: "one_p_dz"}),
    "N2_1pdzroot_star_1pdz_s": (3, [(0, 1), (0, 2)], {0: "one_p_dz", 1: "one_p_dz", 2: "s"}),
    "N3_1pdz2root_star_s_s": (3, [(0, 1), (0, 2)], {0: "one_p_dz2", 1: "s", 2: "s"}),
    "N4_sroot_path_1pdz_1pdz": (3, [(0, 1), (1, 2)], {0: "s", 1: "one_p_dz", 2: "one_p_dz"}),
    "N6_mix012root_path_pzlap_pz": (3, [(0, 1), (1, 2)], {0: "mix012", 1: "pz_p_lap", 2: "pz"}),
    "N5_path4_1pdzroot_1pdz_pz_s": (4, [(0, 1), (1, 2), (2, 3)], {0: "one_p_dz", 1: "one_p_dz", 2: "pz", 3: "s"}),
}

if __name__ == "__main__":
    out = {}
    for nm in (sys.argv[1:] or list(TREES)):
        n, e, d = TREES[nm]
        out[nm] = analyse(nm, n, e, d)
    Path(sys.argv[0]).with_suffix(".json").write_text(json.dumps(out, indent=1, default=str))
