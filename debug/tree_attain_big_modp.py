"""Least-order a (Lemma 1 of debug/tree_attain_memo.md) for larger s-type trees, modulo p = 2147483629.

Same module model as debug/tree_attain_review_modp.py and tests/_tree_attain.py:
  theta z = 2z; theta K0 = D; theta D = z K0; theta P_y = P'_y; theta P'_y = -P'_y - z K0 prod_{c child of y} P_c,
coordinates (S, kept, i); H = z * (empty S, all root children in state 0, i = 0); pure = |S| = n - 1.

Speed-ups over the reviewer's script (which solved one large nullspace for the polynomial a directly):
  1. Orbit reduction.  H is invariant under the automorphisms of the rooted tree and theta is equivariant, so
     A H lives in the invariant part.  Coordinates are orbit sums e_O (orbit = decorated rooted tree up to
     isomorphism, plus i).  theta e_O = sum_O' (alpha + beta z) e_O' with integer alpha, beta.  a H is pure iff
     all non-pure orbit coefficients vanish (orbits have disjoint supports).
  2. The least theta-order N is the first j with theta^0 H .. theta^j H dependent mod pure over C(z) (rank at
     random points).  The relation over C(z) is one-dimensional; its primitive polynomial representative is the
     least-z-degree a (any polynomial solution is a C[z]-multiple of it), i.e. the reviewer's a.  We get it by
     solving at m random points (c_N = 1), rational reconstruction of each c_j(z), clearing denominators and
     dividing by the content gcd.  m doubles until the EXACT check (a H pure, polynomial arithmetic mod p) passes.
Outputs theta-order, z-degree, weight, the symbol sigma(a)(1, -t) mod p with its roots, value at t = -(n-1),
and the purity counts.
"""
from __future__ import annotations

import json
import os
import sys
import time
from fractions import Fraction
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
from flint import nmod_mat, nmod_poly

P = int(os.environ.get("TREE_BIG_PRIME", "2147483629"))
ROOT = Path(__file__).resolve().parent.parent
T0 = time.time()


def log(*a) -> None:
    print(f"[{time.time() - T0:8.1f}s]", *a, flush=True)


# ----------------------------------------------------------------------------------------------- tree + keys
def make_children(n: int, edges, root: int) -> Dict[int, List[int]]:
    adj = {v: [] for v in range(n)}
    for a, b in edges:
        adj[a].append(b)
        adj[b].append(a)
    ch: Dict[int, List[int]] = {}

    def rec(v, p):
        ch[v] = sorted(w for w in adj[v] if w != p)
        for w in ch[v]:
            rec(w, v)
    rec(root, -1)
    assert len(ch) == n
    return ch


def theta_key(ch, key):
    """theta(key) = sum (alpha + beta z) key' ; returns list of (key', alpha, beta).  z-derivative part separate."""
    S, kept, i = key
    k = len(S)
    out = []
    if k - i > 0:
        out.append(((S, kept, i + 1), k - i, 0))
    if i > 0:
        out.append(((S, kept, i - 1), 0, i))
    for idx, (y, st) in enumerate(kept):
        rest = kept[:idx] + kept[idx + 1:]
        if st == 0:
            out.append(((S, tuple(sorted(rest + ((y, 1),))), i), 1, 0))
        else:
            out.append(((S, kept, i), -1, 0))
            newkept = tuple(sorted(rest + tuple((w, 0) for w in ch[y])))
            out.append(((S | frozenset([y]), newkept, i), 0, -1))
    return out


def canon(ch, root, key) -> Tuple[str, int]:
    S, kept, i = key
    kd = dict(kept)

    def enc(v):
        if v == root:
            lab = "R"
        elif v in S:
            lab = "S"
        elif v in kd:
            lab = "k%d" % kd[v]
        else:
            lab = "x"
        return lab + "(" + ",".join(sorted(enc(c) for c in ch[v])) + ")"
    return enc(root), i


def build_orbit_model(n, edges, root, use_orbits=True):
    ch = make_children(n, edges, root)
    key0 = (frozenset(), tuple(sorted((w, 0) for w in ch[root])), 0)
    seen = {key0}
    stack = [key0]
    th = {}
    while stack:
        kk = stack.pop()
        th[kk] = theta_key(ch, kk)
        for k2, _, _ in th[kk]:
            if k2 not in seen:
                seen.add(k2)
                stack.append(k2)
    keys = sorted(seen, key=lambda q: (len(q[0]), str(q)))
    if use_orbits:
        lab = {kk: canon(ch, root, kk) for kk in keys}
    else:
        lab = {kk: kk for kk in keys}
    orbs = sorted(set(lab.values()), key=str)
    oid = {o: j for j, o in enumerate(orbs)}
    size = np.zeros(len(orbs), dtype=np.int64)
    pure = np.zeros(len(orbs), dtype=bool)
    for kk in keys:
        size[oid[lab[kk]]] += 1
        pure[oid[lab[kk]]] = len(kk[0]) == n - 1
    R = len(orbs)
    Ma = [[Fraction(0)] * R for _ in range(R)]
    Mb = [[Fraction(0)] * R for _ in range(R)]
    for kk in keys:
        src = oid[lab[kk]]
        for k2, al, be in th[kk]:
            dst = oid[lab[k2]]
            Ma[dst][src] += Fraction(al, int(size[dst]))
            Mb[dst][src] += Fraction(be, int(size[dst]))
    A = np.zeros((R, R), dtype=np.int64)
    B = np.zeros((R, R), dtype=np.int64)
    for a_ in range(R):
        for b_ in range(R):
            assert Ma[a_][b_].denominator == 1 and Mb[a_][b_].denominator == 1, "orbit matrix not integral"
            A[a_, b_] = int(Ma[a_][b_])
            B[a_, b_] = int(Mb[a_][b_])
    h0 = oid[lab[key0]]
    assert size[h0] == 1
    return dict(R=R, A=A, B=B, pure=pure, h0=h0, nkeys=len(keys), size=size)


def theta_step(M, F: np.ndarray) -> np.ndarray:
    """F: (R, L) coefficients mod P of sum_O f_O(z) e_O.  Returns theta F, shape (R, L+1)."""
    R, L = F.shape
    G = np.zeros((R, L + 1), dtype=np.int64)
    ks = (2 * np.arange(L, dtype=np.int64)) % P
    G[:, :L] = (F * ks[None, :]) % P
    G[:, :L] = (G[:, :L] + (M["A"] @ F) % P) % P
    G[:, 1:] = (G[:, 1:] + (M["B"] @ F) % P) % P
    return G


# ----------------------------------------------------------------------------------------------- mod-p linear algebra
def matmul_mod(X: np.ndarray, Y: np.ndarray) -> np.ndarray:
    """(X @ Y) mod P exactly, X, Y int64 in [0, P), via 16-bit limbs and float64 BLAS (exact below 2^53)."""
    assert X.shape[1] <= (1 << 20)
    x1, x0 = (X >> 16).astype(np.float64), (X & 0xFFFF).astype(np.float64)
    y1, y0 = (Y >> 16).astype(np.float64), (Y & 0xFFFF).astype(np.float64)
    hh = (x1 @ y1).astype(np.int64) % P
    mid = ((x1 @ y0).astype(np.int64) + (x0 @ y1).astype(np.int64)) % P
    ll = (x0 @ y0).astype(np.int64) % P
    s16 = (1 << 16) % P
    s32 = (1 << 32) % P
    return (hh * s32 % P + mid * s16 % P + ll) % P


def powers(zs: np.ndarray, L: int) -> np.ndarray:
    V = np.ones((L, len(zs)), dtype=np.int64)
    for k in range(1, L):
        V[k] = (V[k - 1] * zs) % P
    return V


def evals(Fs: List[np.ndarray], rows: np.ndarray, zs: np.ndarray) -> np.ndarray:
    """values[j, r, i] = (theta^j H)[rows[r]](zs[i])."""
    L = max(F.shape[1] for F in Fs)
    V = powers(zs, L)
    stack = np.zeros((len(Fs) * len(rows), L), dtype=np.int64)
    for j, F in enumerate(Fs):
        stack[j * len(rows):(j + 1) * len(rows), :F.shape[1]] = F[rows]
    out = matmul_mod(stack, V)
    return out.reshape(len(Fs), len(rows), len(zs))


def interpolate(zs: np.ndarray, Y: np.ndarray, Mz: nmod_poly, chunk: int = 512) -> np.ndarray:
    """Coefficients (m, N) of the polynomials of degree < m with values Y[i, :] at zs[i]."""
    m = len(zs)
    c = np.array([int(v) for v in Mz.coeffs()], dtype=np.int64)  # length m+1, monic
    out = np.zeros((m, Y.shape[1]), dtype=np.int64)
    for s0 in range(0, m, chunk):
        a = zs[s0:s0 + chunk]
        Q = np.zeros((m, len(a)), dtype=np.int64)  # Q[k, i] = coeff k of M(z)/(z - a_i)
        Q[m - 1] = c[m]
        for k in range(m - 1, 0, -1):
            Q[k - 1] = (c[k] + a * Q[k]) % P
        dv = np.zeros(len(a), dtype=np.int64)  # M'(a_i) = q_i(a_i), Horner
        for k in range(m - 1, -1, -1):
            dv = (dv * a + Q[k]) % P
        winv = np.array([pow(int(v), P - 2, P) for v in dv], dtype=np.int64)
        Yw = (Y[s0:s0 + chunk] * winv[:, None]) % P
        out = (out + matmul_mod(Q, Yw)) % P
    return out


def ratrecon(f: nmod_poly, Mz: nmod_poly, m: int):
    r0, r1 = Mz, f
    t0, t1 = nmod_poly([0], P), nmod_poly([1], P)
    while r1.degree() >= m // 2:
        q, r = divmod(r0, r1)
        r0, r1 = r1, r
        t0, t1 = t1, t0 - q * t1
    return r1, t1  # numerator, denominator


def least_order(name, n, edges, root, use_orbits=True, rng_seed=20261119, m0=None, m_cap=40000, return_a=False):
    rng = np.random.default_rng(rng_seed)
    t_start = time.time()
    M = build_orbit_model(n, edges, root, use_orbits)
    R = M["R"]
    nonpure = np.where(~M["pure"])[0]
    purei = np.where(M["pure"])[0]
    log(f"{name}: n={n} root={root}: {M['nkeys']} keys -> {R} coordinates ({len(nonpure)} non-pure, "
        f"{len(purei)} pure); max orbit size {int(M['size'].max())}")
    F = np.zeros((R, 2), dtype=np.int64)
    F[M["h0"], 1] = 1
    Fs = [F]
    zt = rng.integers(2, P - 1, size=2, dtype=np.int64)
    N = None
    for j in range(1, len(nonpure) + 2):
        Fs.append(theta_step(M, Fs[-1]))
        E = evals(Fs, nonpure, zt)
        rk = max(nmod_mat(E[:, :, i].T.tolist(), P).rank() for i in range(2))
        if rk < j + 1:
            N = j
            break
    assert N is not None
    t_order = time.time() - t_start
    log(f"{name}: least theta-order N = {N} (bound {len(nonpure)}); {t_order:.1f}s")
    # row selection: N rows independent on the first N columns at a random point
    zr = rng.integers(2, P - 1, size=1, dtype=np.int64)
    E = evals(Fs, nonpure, zr)[:, :, 0].T  # (Rnp, N+1)
    Mt = nmod_mat(E[:, :N].T.tolist(), P)  # (N, Rnp): pick pivot columns
    rr = Mt.rref()[0]
    piv = []
    rl = rr.tolist()
    for row in rl:
        for c, v in enumerate(row):
            if int(v) != 0:
                piv.append(c)
                break
    piv = piv[:N]
    assert len(piv) == N
    rows = nonpure[np.array(piv)]
    m = m0 or max(64, 8 * N)
    while True:
        t1 = time.time()
        zs = np.unique(rng.integers(2, P - 1, size=m + 10, dtype=np.int64))[:m]
        C = np.zeros((m, N), dtype=np.int64)
        cb = max(16, int(2e8 // ((N + 1) * N * 8)))
        for s0 in range(0, m, cb):
            Ev = evals(Fs, rows, zs[s0:s0 + cb])  # (N+1, N, chunk)
            for ii in range(Ev.shape[2]):
                A_ = nmod_mat(Ev[:N, :, ii].T.tolist(), P)
                b_ = nmod_mat([[(-int(v)) % P] for v in Ev[N, :, ii]], P)
                x = A_.solve(b_)
                C[s0 + ii] = [int(v) for v in x.entries()]
        t_solve = time.time() - t1
        # interpolate each c_j through the m points (Lagrange, chunked synthetic division + BLAS)
        Mz = nmod_poly([1], P)
        for z in zs:
            Mz = Mz * nmod_poly([(-int(z)) % P, 1], P)
        cl = interpolate(zs, C, Mz)  # (m, N) ndarray, coefficient k of c_j at [k, j]
        nums, dens = [], []
        ok = True
        for j in range(N):
            f = nmod_poly(cl[:, j].tolist(), P)
            num, den = ratrecon(f, Mz, m)
            if den.is_zero() or num.degree() + den.degree() > m - 20:
                ok = False
                break
            nums.append(num)
            dens.append(den)
        t_rec = time.time() - t1 - t_solve
        if ok:
            Lc = nmod_poly([1], P)
            for den in dens:
                g = Lc.gcd(den)
                Lc = Lc * (den // g)
            a_polys = [nums[j] * (Lc // dens[j]) for j in range(N)] + [Lc]
            g = a_polys[0]
            for q in a_polys[1:]:
                g = g.gcd(q)
            a_polys = [q // g for q in a_polys]
            d = max(q.degree() for q in a_polys)
            # exact purity check on all coordinates
            Lmax = max(Fx.shape[1] for Fx in Fs) + d + 2
            tot = [nmod_poly([0], P) for _ in range(R)]
            for j in range(N + 1):
                aj = a_polys[j]
                if aj.is_zero():
                    continue
                Fj = Fs[j]
                for o in range(R):
                    row = Fj[o]
                    if row.any():
                        tot[o] = tot[o] + aj * nmod_poly([int(v) for v in row], P)
            nonpure_left = sum(1 for o in nonpure if not tot[o].is_zero())
            pure_nz = sum(1 for o in purei if not tot[o].is_zero())
            log(f"{name}: m={m}: solves {t_solve:.1f}s, reconstruction {t_rec:.1f}s, z-degree {d}, "
                f"non-pure left {nonpure_left}, pure nonzero {pure_nz}")
            if nonpure_left == 0:
                break
        else:
            log(f"{name}: m={m}: reconstruction degree too high, doubling")
        m = int(m * 1.6)
        if m > m_cap:
            return {"label": name, "n": n, "root": root, "theta_order": N, "status": f"m above cap {m_cap}"}
    w = max(2 * q.degree() + j for j, q in enumerate(a_polys) if not q.is_zero())
    lead = [0] * (N + 1)
    for j, q in enumerate(a_polys):
        if (w - j) % 2 == 0 and not q.is_zero():
            k = (w - j) // 2
            if k <= q.degree():
                lead[j] = (int(q.coeffs()[k]) * (-1) ** j) % P
    lp = nmod_poly(lead, P)
    pt = -(n - 1)
    val = int(lp(pt % P))
    rem = lp
    roots = {}
    for rt in range(-(n + 1), n + 2):
        lin = nmod_poly([(-rt) % P, 1], P)
        mult = 0
        while rem.degree() > 0:
            q, r = divmod(rem, lin)
            if not r.is_zero():
                break
            rem, mult = q, mult + 1
        if mult:
            roots[rt] = mult
    # any other roots mod p?  (factor what is left)
    other = [] if rem.degree() <= 0 else [str(f) for f, e in rem.factor()[1]]
    res = {"label": name, "n": n, "root": root, "edges": edges, "prime": P, "method": "least-order (C(z) relation, "
           "rational reconstruction, exact mod-p purity check)", "orbit_reduced": use_orbits,
           "coordinates": R, "nonpure_coordinates": int(len(nonpure)), "theta_order": N, "z_degree": d,
           "weight": w, "symbol_degree_in_t": lp.degree(), "symbol_at_-(n-1)_mod_p": val,
           "nonzero_at_-(n-1)": bool(val != 0), "roots": {str(k): v for k, v in roots.items()},
           "splits_completely_over_listed_integers": bool(rem.degree() <= 0), "other_factors": other,
           "max_abs_root": max((abs(k) for k in roots), default=0),
           "roots_within_n_minus_2": bool(all(abs(k) <= n - 2 for k in roots)),
           "nonpure_left_in_aH": nonpure_left, "pure_coordinates_nonzero_in_aH": pure_nz,
           "points_used": m, "seconds": round(time.time() - t_start, 1)}
    res["verdict"] = "CONFIRMED" if (res["roots_within_n_minus_2"] and res["splits_completely_over_listed_integers"]
                                     and res["nonzero_at_-(n-1)"] and nonpure_left == 0 and pure_nz > 0) else "FAILED"
    log(f"{name}: N={N} d={d} w={w} roots {roots} other {other} nonzero at {pt}: {val != 0}; "
        f"aH non-pure {nonpure_left}, pure nonzero {pure_nz} -> {res['verdict']} ({res['seconds']}s)")
    if return_a:
        return res, a_polys, Fs, M
    return res


TREES = {
    # controls (reviewer mod-p / memo exact values)
    "star3": (3, [(0, 1), (0, 2)], 0),
    "path4_end": (4, [(0, 1), (1, 2), (2, 3)], 0),
    "star4": (4, [(0, 1), (0, 2), (0, 3)], 0),
    "path4_inner": (4, [(0, 1), (1, 2), (2, 3)], 1),
    "path5_centre": (5, [(0, 1), (1, 2), (0, 3), (3, 4)], 0),
    "fork5_root0": (5, [(0, 1), (1, 2), (1, 3), (3, 4)], 0),
    "path6_end": (6, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5)], 0),
    # targets
    "fork5_root_branch": (5, [(0, 1), (1, 2), (1, 3), (3, 4)], 1),
    "star6": (6, [(0, 1), (0, 2), (0, 3), (0, 4), (0, 5)], 0),
    "spider6_root0": (6, [(0, 1), (1, 2), (1, 3), (3, 4), (3, 5)], 0),
    "broom6": (6, [(0, 1), (1, 2), (2, 3), (3, 4), (3, 5)], 0),
    "star7": (7, [(0, k) for k in range(1, 7)], 0),
    "path7_end": (7, [(k, k + 1) for k in range(6)], 0),
    "path7_centre": (7, [(3, 2), (2, 1), (1, 0), (3, 4), (4, 5), (5, 6)], 3),
    "star8": (8, [(0, k) for k in range(1, 8)], 0),
    "star8_leaf": (8, [(0, k) for k in range(1, 8)], 1),
    "caterpillar8": (8, [(0, 1), (1, 2), (1, 3), (3, 4), (3, 5), (5, 6), (5, 7)], 0),
    "spider8_centre": (8, [(0, 1), (1, 2), (0, 3), (3, 4), (0, 5), (5, 6), (6, 7)], 0),
    "path8_end": (8, [(k, k + 1) for k in range(7)], 0),
}

if __name__ == "__main__":
    args = sys.argv[1:]
    use_orbits = "--full" not in args
    names = [a for a in args if not a.startswith("--")] or list(TREES)
    tag = ("" if use_orbits else "_full") + ("" if P == 2147483629 else "_p%d" % P)
    for name in names:
        n, edges, root = TREES[name]
        outp = ROOT / "debug" / "data" / f"tree_attain_big_{name}{tag}.json"
        try:
            res = least_order(name + tag, n, edges, root, use_orbits)
        except MemoryError as e:
            res = {"label": name, "status": "MemoryError %s" % e}
        outp.write_text(json.dumps(res, indent=1))
    log("done")
