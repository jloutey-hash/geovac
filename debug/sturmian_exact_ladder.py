"""Grid-free K-ladder for Paper 60's 1-norm exponents.

Builds ONE exact secular matrix at the largest basis (lmax=3, n<=N_MAX) and then
reads every smaller rung off as a SUBMATRIX, because the Goscinskian bases are
nested: config (l, na, nb) with na, nb <= N is in every basis with N' >= N.
That makes the whole ladder cost one O(K^2) assembly instead of eight.

Reported per rung:
    K                configuration count
    ||M||_1          entrywise 1-norm of M = diag(Z R_nu) + T'     (paper eq:sublinear)
    ||M||_1^diag     sum_i |M_ii|
    ||M||_1^off      sum_{i!=j} |M_ij|  ==  off-diagonal of T'     (paper ||T'||_1^off)
    ||T0||_1         Z sum_nu R_nu                                  (paper eq:sublinear_split)
    ||T'||_1         sum_ij |T'_ij|  INCLUDING the diagonal         (the actual T' block)
    ||T'||_1^diag    sum_i |T'_ii|
"""
from __future__ import annotations

import json
import sys
import time

import numpy as np

import debug.sturmian_exact_slater as ex

Z = 2.0


def build_full(n_max: int, lmax: int = 3, verbose: bool = True, dps: int = 60):
    ex.set_dps(dps)
    ex.assert_precision(50)
    tuples = ex.gen_configs(lmax, {l: n_max for l in range(lmax + 1)})
    cfgs = ex.build_exact_configs(tuples)
    K = len(cfgs)
    if verbose:
        print("full basis: lmax=%d n_max=%d  K=%d" % (lmax, n_max, K), flush=True)
    Tp = np.zeros((K, K))
    t0 = time.time()
    for i in range(K):
        ci = cfgs[i]
        for j in range(i, K):
            cj = cfgs[j]
            g = ci.norm * cj.norm * ex.exact_repulsion_terms(ci.terms, cj.terms)
            Tp[i, j] = Tp[j, i] = -g
        if verbose and (i % 10 == 0 or i == K - 1):
            print("  row %4d/%d  t=%.1fs  Rk-cache=%d"
                  % (i, K, time.time() - t0, len(ex._RK_CACHE)), flush=True)
    Rnu = np.array([c.Rnu for c in cfgs])
    return tuples, Rnu, Tp


def rung(tuples, Rnu, Tp, n_max: int, lmax: int = 3):
    """Sub-select the rung with all n <= n_max and report its norms."""
    idx = [i for i, (l, na, nb) in enumerate(tuples) if nb <= n_max and l <= lmax]
    t = Tp[np.ix_(idx, idx)]
    rn = Rnu[idx]
    K = len(idx)
    M = t.copy()
    M[np.diag_indices(K)] += Z * rn
    d = np.abs(np.diag(M)).sum()
    tot = np.abs(M).sum()
    return dict(
        K=K,
        M1=float(tot),
        M1_diag=float(d),
        M1_off=float(tot - d),
        T0=float(Z * rn.sum()),
        Tp1=float(np.abs(t).sum()),
        Tp1_diag=float(np.abs(np.diag(t)).sum()),
        E0=float(-np.sort(np.linalg.eigvalsh(M))[-1] ** 2 / 2),
    )


def fit(xs, ys):
    return float(np.polyfit(np.log(xs), np.log(ys), 1)[0])


def main(n_max_top: int = 10, lmax: int = 3, out: str = None, dps: int = 60):
    tuples, Rnu, Tp = build_full(n_max_top, lmax, dps=dps)
    rows = [rung(tuples, Rnu, Tp, n, lmax) for n in range(4, n_max_top + 1)]
    rows = [r for r in rows if r["K"] > 0]
    print()
    hdr = ("  N    K      ||M||_1     diag        off        ||T0||_1     "
           "||T'||_1    ||T'||diag    E0")
    print(hdr)
    for n, r in zip(range(4, n_max_top + 1), rows):
        print("%3d %4d  %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f  %.6f"
              % (n, r["K"], r["M1"], r["M1_diag"], r["M1_off"], r["T0"],
                 r["Tp1"], r["Tp1_diag"], r["E0"]))
    if out:
        with open(out, "w") as fh:
            json.dump(rows, fh, indent=1)
        print("\nwrote", out)
    return rows


if __name__ == "__main__":
    nm = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    o = sys.argv[2] if len(sys.argv) > 2 else None
    dd = int(sys.argv[3]) if len(sys.argv) > 3 else 60
    main(nm, 3, o, dd)
