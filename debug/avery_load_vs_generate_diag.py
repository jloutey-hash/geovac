"""Load-vs-generate diagnostic for Avery's on-device angular-generation idea (atoms).

Question: for a single-centre atomic basis (all n <= N, l < n, all m), is it
cheaper to LOAD every distinct two-electron integral into a qubitized block
encoding, or to load only (i) an angular table that depends on l_max alone and
(ii) the m-independent radial table, and combine them on the device?

Scheme LOAD  : sparse coherent alias sampling over the L distinct nonzero
               integrals (Berry et al., Quantum 3, 208 (2019)); one LCU term
               per integral, lambda_V = sum |(ab|cd)|.
Scheme FACT  : LCU term per (k, ab, cd) multipole contribution, prepared in two
               stages: stage 1 over (k, l's, m's) with weight
               |angular| * sum_n |R^k|; stage 2 conditioned on (k, l's) over the
               n-quartet with weight |R^k| / sum_n |R^k|.  This is exact (the
               product of the two stage weights is the term weight) and needs
               no on-device 3j arithmetic.  Cost: two QROAMs, of sizes T1
               (angular, l_max only) and T2 (radial).  Price: lambda_F =
               sum_k sum |term_k| >= lambda_V (multipole terms no longer cancel).

QROAM Toffoli cost (clean ancilla): min over k=2^j of ceil(L/k) + b(k-1).
Word size b = keep bits (20) + ceil(log2 table) alt-index bits, for BOTH
schemes; orbital-index words are NOT charged to LOAD (conservative toward LOAD).
Symmetry-unique counts = ordered nonzero count / 4 (4-fold symmetry of complex-
orbital integrals; slightly undercounts LOAD's table, again conservative).

Orbitals: hydrogenic, Z=1 (all integrals scale by Z; ratios are Z-free).
Angular: exact Gaunt integrals (sympy).  Radial: geovac.hypergeometric_slater.
"""
from __future__ import annotations

import json
import math
import sys
import time
from functools import lru_cache
from typing import Dict, List, Tuple

import numpy as np
from sympy.physics.wigner import gaunt

sys.path.insert(0, ".")
from geovac.hypergeometric_slater import get_rk_float  # noqa: E402

KEEP_BITS = 20


@lru_cache(maxsize=None)
def G(l1: int, m1: int, l2: int, m2: int, l3: int, m3: int) -> float:
    """Gaunt integral  int Y_l1m1 Y_l2m2 Y_l3m3 dOmega  (complex, CS phase)."""
    return float(gaunt(l1, l2, l3, m1, m2, m3))


def A1(la: int, ma: int, lb: int, mb: int, k: int, q: int) -> float:
    """int Y*_{la ma} Y*_{kq} Y_{lb mb}."""
    return (-1) ** ((ma + q) & 1) * G(la, -ma, k, -q, lb, mb)


def A2(lc: int, mc: int, ld: int, md: int, k: int, q: int) -> float:
    """int Y*_{lc mc} Y_{kq} Y_{ld md}."""
    return (-1) ** (mc & 1) * G(lc, -mc, k, q, ld, md)


def orbitals(N: int) -> List[Tuple[int, int, int]]:
    return [(n, l, m) for n in range(1, N + 1) for l in range(n) for m in range(-l, l + 1)]


def eri(o: List[Tuple[int, int, int]], a: int, b: int, c: int, d: int) -> float:
    """(ab|cd) = int phi_a* phi_b (1/r12) phi_c* phi_d, direct (validation path)."""
    (na, la, ma), (nb, lb, mb), (nc, lc, mc), (nd, ld, md) = o[a], o[b], o[c], o[d]
    q = mb - ma
    if mc - md != q:
        return 0.0
    tot = 0.0
    for k in range(abs(la - lb), la + lb + 1):
        if (la + lb + k) % 2 or k < abs(lc - ld) or k > lc + ld or (lc + ld + k) % 2:
            continue
        ang = 4 * math.pi / (2 * k + 1) * A1(la, ma, lb, mb, k, q) * A2(lc, mc, ld, md, k, q)
        if ang == 0.0:
            continue
        tot += ang * get_rk_float(na, la, nb, lb, nc, lc, nd, ld, k)
    return tot


def qroam_toffoli(L: float, b: int) -> int:
    if L <= 1:
        return 0
    best = None
    k = 1
    while k <= 4 * L:
        c = math.ceil(L / k) + b * (k - 1)
        best = c if best is None else min(best, c)
        k *= 2
    return best


def word(L: float) -> int:
    return KEEP_BITS + max(1, math.ceil(math.log2(max(L, 2))))


def validate() -> Dict[str, float]:
    o = orbitals(2)
    idx = {x: i for i, x in enumerate(o)}
    s1, s2, p0 = idx[(1, 0, 0)], idx[(2, 0, 0)], idx[(2, 1, 0)]
    checks = {
        "J(1s,1s) = 5/8": (eri(o, s1, s1, s1, s1), 5 / 8),
        "J(1s,2s) = 17/81": (eri(o, s1, s1, s2, s2), 17 / 81),
        "K(1s,2s) = 16/729": (eri(o, s1, s2, s2, s1), 16 / 729),
        "(2p0 2p0|2p0 2p0) = 93/512 + (4/25)45/512": (eri(o, p0, p0, p0, p0), 93 / 512 + 4 / 25 * 45 / 512),
    }
    out = {}
    for name, (got, want) in checks.items():
        err = abs(got - want)
        print(f"  check {name:45s} got {got:.12f} want {want:.12f} err {err:.1e}")
        assert err < 1e-10, name
        out[name] = err
    # symmetry spot check
    rng = np.random.default_rng(0)
    for _ in range(200):
        a, b, c, d = rng.integers(0, len(o), 4)
        v = eri(o, a, b, c, d)
        assert abs(v - eri(o, c, d, a, b)) < 1e-12
        assert abs(v - eri(o, b, a, d, c)) < 1e-12
    print("  symmetry (ab|cd)=(cd|ab)=(ba|dc): 200 random quartets OK")
    return out


def run(N: int) -> Dict[str, float]:
    o = orbitals(N)
    M = len(o)
    shells = sorted({(n, l) for n, l, _ in o})
    lmax = N - 1
    # pair lists keyed by q: (orbital pair index (a,b), shell pair, l pair, m pair)
    pairs_by_q: Dict[int, List[Tuple[int, int]]] = {}
    for a in range(M):
        for b in range(M):
            q = o[b][2] - o[a][2]
            pairs_by_q.setdefault(q, []).append((a, b))

    lam_V = 0.0
    lam_F = 0.0
    nnz_V = 0
    nnz_terms = 0
    ang_keys = set()      # (k, la, ma, lb, mb, lc, mc, ld, md) nonzero
    rad_keys = set()      # (k, na, la, nb, lb, nc, lc, nd, ld) used

    for q, plist in pairs_by_q.items():
        P = len(plist)
        V = np.zeros((P, P))
        for k in range(0, 2 * lmax + 1):
            a1 = np.array([A1(o[a][1], o[a][2], o[b][1], o[b][2], k, q)
                           if (o[a][1] + o[b][1] + k) % 2 == 0 and abs(o[a][1] - o[b][1]) <= k <= o[a][1] + o[b][1]
                           else 0.0 for a, b in plist])
            # electron-2 factor for pair (c,d) uses q = mc - md, i.e. pair stored as (d,c)?  No:
            # plist holds (a,b) with q = mb - ma.  For electron 2 we need pairs (c,d) with
            # mc - md = q, i.e. (d,c) in plist.  Build A2 on reversed pairs.
            a2 = np.array([A2(o[d][1], o[d][2], o[c][1], o[c][2], k, q)
                           if (o[c][1] + o[d][1] + k) % 2 == 0 and abs(o[c][1] - o[d][1]) <= k <= o[c][1] + o[d][1]
                           else 0.0 for c, d in plist])
            # here row index -> pair (a,b) for electron 1; column index j -> plist[j]=(c',d')
            # interpreted as electron-2 pair (c,d) = (d', c').
            r1 = np.nonzero(a1)[0]
            r2 = np.nonzero(a2)[0]
            if len(r1) == 0 or len(r2) == 0:
                continue
            Rk = np.zeros((len(r1), len(r2)))
            for ii, i in enumerate(r1):
                a, b = plist[i]
                for jj, j in enumerate(r2):
                    dp, cp = plist[j]
                    c, d = cp, dp
                    Rk[ii, jj] = get_rk_float(o[a][0], o[a][1], o[b][0], o[b][1],
                                              o[c][0], o[c][1], o[d][0], o[d][1], k)
                    rad_keys.add((k, o[a][0], o[a][1], o[b][0], o[b][1],
                                  o[c][0], o[c][1], o[d][0], o[d][1]))
                    ang_keys.add((k, o[a][1], o[a][2], o[b][1], o[b][2],
                                  o[c][1], o[c][2], o[d][1], o[d][2]))
            term = (4 * math.pi / (2 * k + 1)) * np.outer(a1[r1], a2[r2]) * Rk
            lam_F += float(np.abs(term).sum())
            nnz_terms += int(np.count_nonzero(np.abs(term) > 1e-14))
            V[np.ix_(r1, r2)] += term
        lam_V += float(np.abs(V).sum())
        nnz_V += int(np.count_nonzero(np.abs(V) > 1e-14))

    L = nnz_V / 4
    T1 = len(ang_keys) / 4
    T2 = len(rad_keys) / 4
    c_load = qroam_toffoli(L, word(L))
    c_fact = qroam_toffoli(T1, word(T1)) + qroam_toffoli(T2, word(T2))
    res = dict(N=N, M=M, Q=2 * M, shells=len(shells), lmax=lmax,
               L_load=L, T1_angular=T1, T2_radial=T2,
               table_ratio=L / (T1 + T2),
               toffoli_load=c_load, toffoli_fact=c_fact,
               toffoli_ratio=c_load / c_fact,
               lambda_V=lam_V, lambda_F=lam_F, lambda_ratio=lam_F / lam_V,
               net_ratio_load_over_fact=(lam_V * c_load) / (lam_F * c_fact))
    return res


if __name__ == "__main__":
    print("Validation (hydrogenic Z=1, textbook Slater integrals):")
    val = validate()
    rows = []
    Nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    for N in range(2, Nmax + 1):
        t0 = time.time()
        r = run(N)
        r["seconds"] = round(time.time() - t0, 1)
        rows.append(r)
        print(f"N={N} M={r['M']:3d} L={r['L_load']:12.0f} T1={r['T1_angular']:9.0f} T2={r['T2_radial']:9.0f} "
              f"tab={r['table_ratio']:7.1f}x  Toff load={r['toffoli_load']:7d} fact={r['toffoli_fact']:6d} "
              f"({r['toffoli_ratio']:5.2f}x)  lamF/lamV={r['lambda_ratio']:5.2f}  "
              f"NET={r['net_ratio_load_over_fact']:5.2f}x  [{r['seconds']}s]", flush=True)
    with open("debug/data/avery_load_vs_generate_diag.json", "w") as f:
        json.dump(dict(validation=val, rows=rows), f, indent=1)
