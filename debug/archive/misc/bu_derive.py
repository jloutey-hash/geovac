"""Record BU derivation mod p: from the first-order system U_{L+1} = A(L) U_L, F_L = c(L) . U_L (debug/bu_system.py,
derived from the Bessel rules and IBP; no data), derive
  (1) the recurrence in L by the cyclic-vector method: rows r_k(L) = c(L+k) A(L+k-1) ... A(L); the first k = m with
      r_0..r_m dependent over F_p(L) gives sum_k p_k(L) F_{L+k} = 0 for EVERY solution U of the system;
  (2) the bond-angle operator: constants a_ij with D = sum a_ij x^j T^i such that the coefficient-space image of D,
      applied to the row-vector sequence n -> r_{n-lam}(lam), vanishes identically in lam (i.e. D annihilates
      c . U for every solution U).  Both are linear algebra on the derived symbolic identities, evaluated at random
      points lam of F_p (Schwartz-Zippel); no coefficient data of the target is used.
"""
import json
import random
import sys
import time

import flint
import sympy as sp

PRIMES = [(1 << 61) - 1, 4611686018427387847]


BOUNDARY = True         # impose the identity at integer centres l < R+Q from L = 0 (see boundary_rows)
VERSION = 'v1'          # 'v1': 18-dim system with two masters per piece; 'v2': 14-dim rank-one system


def load_system():
    cache = f'data/bu_system_{VERSION}.json'
    try:
        d = json.load(open(cache))
        return d['N'], d['A'], d['c']
    except (OSError, ValueError, KeyError):
        pass
    import bu_system as S
    A, c = S.build() if VERSION == 'v1' else S.build1()
    Nn = A.shape[0]
    L = S.L

    def nd(e):
        num, den = sp.fraction(sp.together(e))
        pn, pd = sp.Poly(num, L), sp.Poly(den, L)
        # integer coefficients, low degree first
        cn = [sp.Rational(v) for v in reversed(pn.all_coeffs())]
        cd = [sp.Rational(v) for v in reversed(pd.all_coeffs())]
        return [[str(v) for v in cn], [str(v) for v in cd]]
    Aj = {f'{i},{j}': nd(A[i, j]) for i in range(Nn) for j in range(Nn) if A[i, j] != 0}
    cj = {str(j): nd(c[0, j]) for j in range(Nn) if c[0, j] != 0}
    idx = S.IDX if VERSION == 'v1' else S.IDX1
    json.dump(dict(N=Nn, A=Aj, c=cj, idx={f'{k[0]}:{k[1]}': v for k, v in idx.items()}), open(cache, 'w'), indent=0)
    return Nn, Aj, cj


class Sys:
    def __init__(self, p):
        self.p = p
        N, Aj, cj = load_system()
        self.N = N

        def poly(lst):
            vals = [sp.Rational(v) for v in lst]
            return flint.nmod_poly([int(flint.nmod(v.p, p) / flint.nmod(v.q, p)) for v in vals], p)
        self.A = {tuple(map(int, k.split(','))): (poly(v[0]), poly(v[1])) for k, v in Aj.items()}
        self.c = {int(k): (poly(v[0]), poly(v[1])) for k, v in cj.items()}

    def Amat(self, lam):
        lam = flint.nmod(int(lam), self.p)
        M = flint.nmod_mat(self.N, self.N, self.p)
        for (i, j), (n, d) in self.A.items():
            M[i, j] = n(lam) / d(lam)
        return M

    def cvec(self, lam):
        lam = flint.nmod(int(lam), self.p)
        v = flint.nmod_mat(1, self.N, self.p)
        for j, (n, d) in self.c.items():
            v[0, j] = n(lam) / d(lam)
        return v

    def rows(self, lam, kmax):
        """r_k(lam), k = 0..kmax, as 1 x N nmod_mats:  r_k = c(lam+k) A(lam+k-1) ... A(lam)."""
        P = flint.nmod_mat(self.N, self.N, [1 if i == j else 0 for i in range(self.N) for j in range(self.N)], self.p)
        out = []
        for k in range(kmax + 1):
            out.append(self.cvec(lam + k) * P)
            P = self.Amat(lam + k) * P
        return out


def stack(rows, p):
    N = rows[0].ncols()
    return flint.nmod_mat(len(rows), N, [int(r[0, j]) for r in rows for j in range(N)], p)


def recurrence_order(S, trials=3):
    rnd = random.Random(1)
    orders = []
    for _ in range(trials):
        lam = rnd.randrange(10 ** 6, S.p - 10 ** 6)
        rows = S.rows(lam, S.N + 1)
        for m in range(1, S.N + 2):
            if stack(rows[:m + 1], S.p).rank() <= m:
                orders.append(m)
                break
    assert len(set(orders)) == 1, orders
    return orders[0]


def normalized_coeffs(S, lam, m):
    """p_k(lam) / p_m(lam), k = 0..m, from the 1-dim left kernel of [r_0; ...; r_m]."""
    M = stack(S.rows(lam, m), S.p)
    Mt = M.transpose()                      # N x (m+1); right kernel = left kernel of M
    X, nul = Mt.nullspace()
    assert nul == 1, nul
    v = [X[i, 0] for i in range(m + 1)]
    return [x / v[m] for x in v]


def ratrecon(vals, pts, p):
    """Rational function r/s through (pts, vals) with deg r + deg s < len(pts), by extended Euclid; s monic."""
    P = flint.nmod_poly.interpolate(pts, vals) if hasattr(flint.nmod_poly, 'interpolate') else None
    if P is None:
        # Lagrange / Newton interpolation
        P = flint.nmod_poly([0], p)
        Mx = flint.nmod_poly([1], p)
        for x, yv in zip(pts, vals):
            xv = flint.nmod(x, p)
            cur = P(xv)
            corr = (flint.nmod(int(yv), p) - cur) / Mx(xv)
            P = P + Mx * flint.nmod_poly([int(corr)], p)
            Mx = Mx * flint.nmod_poly([-x % p, 1], p)
    else:
        Mx = flint.nmod_poly([1], p)
        for x in pts:
            Mx = Mx * flint.nmod_poly([-x % p, 1], p)
    n = len(pts)
    r0, r1 = Mx, P
    t0, t1 = flint.nmod_poly([0], p), flint.nmod_poly([1], p)
    best = None
    while r1.degree() >= 0:
        if r1.degree() + t1.degree() < n and (best is None or r1.degree() + t1.degree() < best[2]):
            best = (r1, t1, r1.degree() + t1.degree())
        q, rr = divmod(r0, r1)
        r0, r1 = r1, rr
        t0, t1 = t1, t0 - q * t1
    r, s, _ = best
    lc = s.coeffs()[-1]
    return r * flint.nmod_poly([int(1 / flint.nmod(int(lc), p))], p), s * flint.nmod_poly([int(1 / flint.nmod(int(lc), p))], p)


def derive_recurrence(S, m, npts):
    rnd = random.Random(2)
    pts = [rnd.randrange(10 ** 6, S.p - 10 ** 6) for _ in range(npts)]
    vals = [normalized_coeffs(S, lam, m) for lam in pts]
    fr = [ratrecon([v[k] for v in vals], pts, S.p) for k in range(m + 1)]
    # common denominator -> polynomial coefficients
    den = flint.nmod_poly([1], S.p)
    for _, s in fr:
        den = den * s // den.gcd(s)
    polys = [r * (den // s) for r, s in fr]
    g = polys[0]
    for q in polys[1:]:
        g = g.gcd(q)
    polys = [q // g for q in polys]
    return polys


def check_recurrence(polys, seq, p):
    m = len(polys) - 1
    bad = []
    for l in range(len(seq) - m):
        lv = flint.nmod(l, p)
        s = sum((polys[k](lv) * flint.nmod(seq[l + k], p) for k in range(m + 1)), flint.nmod(0, p))
        if int(s) != 0:
            bad.append(l)
    return bad


# ---------------- angle operator from the system ----------------

def Xv(seq, n0, p):
    """local multiplication by x on a list of row vectors indexed n = n0 + t (n0 a field element)."""
    out = []
    for t in range(1, len(seq) - 1):
        n = flint.nmod(n0 + t, p)
        a = n / (2 * n - 1)
        b = (n + 1) / (2 * n + 3)
        out.append(seq[t - 1] * a + seq[t + 1] * b)
    return out


def Tv(seq, n0, p):
    out = []
    for t in range(1, len(seq) - 1):
        n = flint.nmod(n0 + t, p)
        a = -(n - 1) * n / (2 * n - 1)
        b = (n + 1) * (n + 2) / (2 * n + 3)
        out.append(seq[t - 1] * a + seq[t + 1] * b)
    return out


def op_rows(S, lam, R, Q):
    """For each (i, j): the row vector (x^j T^i applied to n -> r_{n-lam}(lam)) at the centre n = lam + R + Q.
    Computed as e_c^T X^j T^i (left row vectors on the window n = lam + t, t = 0..2K), then times the window."""
    K = R + Q
    rows = S.rows(lam, 2 * K)
    return op_rows_core(stack(rows, S.p), lam, R, Q, S.p, S.N)


def boundary_rows(S, l, R, Q, base=None):
    """Integer centre l >= 0: the window n = l-K..l+K with rows r_n(0) = c(n) A(n-1)...A(0) for n >= 0 and zero
    for n < 0 (the coefficient-space X, T never couple n < 0 to n >= 0: their n=0 -> -1 and -1 -> 0 coefficients
    vanish), i.e. the identity for every solution started at L = 0, with no passage through negative L."""
    K = R + Q
    if base is None:
        base = S.rows(0, l + K)
    p = S.p
    W = 2 * K + 1
    ent = []
    for t in range(W):
        n = l - K + t
        ent.extend([int(base[n][0, j]) for j in range(S.N)] if n >= 0 else [0] * S.N)
    Rm = flint.nmod_mat(W, S.N, ent, p)
    return op_rows_core(Rm, l - K, R, Q, p, S.N)


def op_rows_core(Rm, n0, R, Q, p, N):
    K = R + Q
    W = 2 * K + 1
    nn = [flint.nmod(n0 + t, p) for t in range(W)]
    # (X f)_n = a_n f_{n-1} + b_n f_{n+1};  (T f)_n = c_n f_{n-1} + d_n f_{n+1}
    xa = [n / (2 * n - 1) for n in nn]
    xb = [(n + 1) / (2 * n + 3) for n in nn]
    tc = [-(n - 1) * n / (2 * n - 1) for n in nn]
    td = [(n + 1) * (n + 2) / (2 * n + 3) for n in nn]
    zero = flint.nmod(0, p)

    def left(v, a, b):
        """row vector v (indexed by window position) times the band operator: (v Op)_m = v_{m+1} a_{m+1} + v_{m-1} b_{m-1}."""
        out = [zero] * W
        for t, vt in enumerate(v):
            if int(vt) == 0:
                continue
            if t >= 1:
                out[t - 1] += vt * a[t]
            if t + 1 < W:
                out[t + 1] += vt * b[t]
        return out
    keys = [(i, j) for i in range(R + 1) for j in range(Q + 1)]
    vecs = {}
    v = [zero] * W
    v[K] = flint.nmod(1, p)
    xs = [v]
    for j in range(Q):
        xs.append(left(xs[-1], xa, xb))
    for j in range(Q + 1):
        u = xs[j]
        for i in range(R + 1):
            vecs[(i, j)] = u
            if i < R:
                u = left(u, tc, td)
    Lm = flint.nmod_mat(len(keys), W, [int(x) for k in keys for x in vecs[k]], p)
    P = Lm * Rm                               # len(keys) x N
    return {k: [P[a, b] for b in range(N)] for a, k in enumerate(keys)}


def angle_operator(S, R, Q, npts, seed=3):
    rnd = random.Random(seed)
    keys = [(i, j) for i in range(R + 1) for j in range(Q + 1)]
    entries = []
    nrow = 0
    for _ in range(npts):
        lam = rnd.randrange(10 ** 6, S.p - 10 ** 6)
        o = op_rows(S, lam, R, Q)
        for comp in range(S.N):
            entries.extend(int(o[k][comp]) for k in keys)
            nrow += 1
    if BOUNDARY:
        K = R + Q
        base = S.rows(0, 2 * K)
        for l in range(K):
            o = boundary_rows(S, l, R, Q, base)
            for comp in range(S.N):
                entries.extend(int(o[k][comp]) for k in keys)
                nrow += 1
    M = flint.nmod_mat(nrow, len(keys), entries, S.p)
    X, nul = M.nullspace()
    vecs = [[X[a, b] for a in range(len(keys))] for b in range(nul)]
    return keys, nul, vecs


if __name__ == '__main__':
    mode = sys.argv[1]
    p = int(sys.argv[2]) if len(sys.argv) > 2 else PRIMES[0]
    if len(sys.argv) > 1 and '--v2' in sys.argv:
        VERSION = 'v2'
        sys.argv.remove('--v2')
    t0 = time.time()
    S = Sys(p)
    RANK = recurrence_order(S)
    if mode == 'order':
        print('system dim', S.N, 'recurrence order', recurrence_order(S), round(time.time() - t0, 1), 's')
    elif mode == 'rec':
        m = RANK
        npts = int(sys.argv[3])
        polys = derive_recurrence(S, m, npts)
        print('order', m, 'degrees', [q.degree() for q in polys], round(time.time() - t0, 1), 's', flush=True)
        json.dump(dict(p=p, order=m, polys=[[int(c) for c in q.coeffs()] for q in polys]),
                  open(f'data/bu_recurrence_{VERSION}_{p}.json', 'w'))
    elif mode == 'angle':
        R, Q = int(sys.argv[3]), int(sys.argv[4])
        npts = int(sys.argv[5]) if len(sys.argv) > 5 else ((R + 1) * (Q + 1)) // RANK + 15
        keys, nul, vecs = angle_operator(S, R, Q, npts)
        print('R', R, 'Q', Q, 'nullity', nul, round(time.time() - t0, 1), 's', flush=True)
        if nul:
            json.dump(dict(p=p, R=R, Q=Q, keys=keys, vecs=[[int(v) for v in vec] for vec in vecs]),
                      open(f'data/bu_angle_{VERSION}_{p}_{R}_{Q}.json', 'w'))


def nullity_auto(S, R, Q, margin=15, seed=3):
    npts = ((R + 1) * (Q + 1)) // RANK + margin
    return angle_operator(S, R, Q, npts, seed)


if __name__ == '__main__' and sys.argv[1] == 'scan':
    out = {}
    for R in range(int(sys.argv[3]), int(sys.argv[4]) + 1):
        lo, hi = -1, None
        Q = 8
        while True:                                  # doubling to find a Q with nullity >= 1
            _, nul, _ = nullity_auto(S, R, Q)
            print('R', R, 'Q', Q, 'nullity', nul, round(time.time() - t0), flush=True)
            if nul:
                hi = Q
                break
            lo = Q
            Q *= 2
            if Q > 400:
                break
        if hi is None:
            out[R] = None
            continue
        while hi - lo > 1:
            mid = (lo + hi) // 2
            _, nul, _ = nullity_auto(S, R, mid)
            print('R', R, 'Q', mid, 'nullity', nul, round(time.time() - t0), flush=True)
            if nul:
                hi = mid
            else:
                lo = mid
        _, nul, _ = nullity_auto(S, R, hi)
        out[R] = (hi, nul)
        print('R', R, 'minQ', hi, 'nullity', nul, flush=True)
    json.dump(out, open(f'data/bu_angle_scan_{VERSION}_{p}.json', 'w'))
