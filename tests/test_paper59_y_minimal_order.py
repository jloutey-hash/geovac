"""Paper 59 sec:density_ring: the minimal operator of the off-slice correlated triangle Y, and its K3 piece.

Records docs/predictions/2026-10-23_price_formula_slots.md (addenda 6, 7), 2026-10-25_geometric_route_y.md and
2026-10-25_tetrahedron_route_y.md.  Written as a separate step from the results; no debug/ import.  Data
tests/data/paper59_y_2026_10_26.json (operators as ascending integer coefficients of c_k(x), L = sum_k c_k theta^k,
theta = t d/dt, x = t^2; theta^k Y at two points from the 60-digit grid).

Five guards, each with the wrong answer it rejects.
(1) M13 (order 13, exact over Q) annihilates Y at t = 1.7 and 2.6, below 1e-40 of its largest x-group.  Rejected:
    M13 with its dominant coefficient at that t times 1 + 1/200.
(2) Y needs the K3 piece: R_S, the right factor of L15 without it, leaves more than 1e-6.  Rejected: "Y lies in
    Sol(R_S)" (the same measure gives below 1e-40 for M13 in (1)).
(3) slow.  The K3 piece R_V is M13's only irreducible left factor: Hom(S*, M13*) mod p (box a 20, b 12, numerator
    degree 500, M13's apparent factor squared) is 1 for R_V and 0 for z_1, z_2, F_1, F_2, F_5.  Rejected: an
    order-one left factor (which would allow an operator of order 12): the planted F_5 o M13 gives 1.
(4) R_V is irreducible and orthogonal: no admissible order-one shape for R_V or R_V*; the one admissible shape in the
    exterior square holds no element; one quadratic relation among its solutions, no linear one among its
    Wronskians.  Rejected: "reducible" (the same code finds both Wronskians of a direct sum of two Legendre
    operators) and "symplectic" (it finds the Wronskian relation of AESZ #1).
(5) Jordan blocks of R_V: 2 + 2 at x = 0, 3 + 1 at infinity, from the log filtration mod p.  Rejected: the other
    pairing (the same code returns 3 + 1 for theta_x^3 (theta_x - 1) and 2 + 2 for theta_x^2 (theta_x - 1)^2).
Not re-run here: the sunrise relation and L15 (about 0.5 MB), the exact maps, the concomitant tests.
"""
from __future__ import annotations

import itertools
import json
from fractions import Fraction as Fr
from math import comb
from pathlib import Path

import mpmath as mp
import numpy as np
import pytest
import sympy as sp
from flint import fmpz_poly, nmod_mat, nmod_poly

DATA = json.loads((Path(__file__).parent / "data" / "paper59_y_2026_10_26.json").read_text())
P1, P2 = 1073741789, 1073740783
CUTS = [Fr(s) for s in DATA["cuts_x"]]
U = sp.Symbol("u")
X = sp.Symbol("x")


def op(key):
    return [fmpz_poly([int(v) for v in c]) for c in DATA[key]]


def from_q(qs):
    """sum_i x^i Q_i(u), u = theta_t -> list of c_k(x)"""
    r = max(sp.Poly(q, U).degree() for q in qs)
    c = [[0] * len(qs) for _ in range(r + 1)]
    for i, q in enumerate(qs):
        for (j,), v in sp.Poly(q, U).terms():
            c[j][i] = int(v)
    return [fmpz_poly(row) for row in c]


# ---------------------------------------------------------------- generic operator algebra (theta = 2 x d/dx)
def dth(p):
    return p.derivative() * p.__class__([0, 2], *([p.modulus()] if hasattr(p, "modulus") else []))


def zero_like(p):
    return p.__class__([], *([p.modulus()] if hasattr(p, "modulus") else []))


def trim(L):
    while len(L) > 1 and L[-1].is_zero():
        L = L[:-1]
    return L


def add(A, B):
    n = max(len(A), len(B))
    z = zero_like(A[0])
    return trim([(A[k] if k < len(A) else z) + (B[k] if k < len(B) else z) for k in range(n)])


def theta_left(B):
    """theta o B"""
    out = [dth(c) for c in B] + [zero_like(B[0])]
    for k, c in enumerate(B):
        out[k + 1] = out[k + 1] + c
    return out


def compose(A, B):
    res = [zero_like(A[0])]
    cur = B
    for k, a in enumerate(A):
        if k:
            cur = theta_left(cur)
        res = add(res, [a * c for c in cur])
    return res


def adjoint(L):
    res = [zero_like(L[0])]
    for k, c in enumerate(L):
        cur = [c]
        for _ in range(k):
            cur = theta_left(cur)
        res = add(res, [(-1) ** k * x for x in cur])
    return res


def modp(L, p):
    return [nmod_poly([int(v) % p for v in c.coeffs()], p) for c in L]


def qform(L):
    """rows i: [coefficient of u^k in Q_i] as ints"""
    D = max(c.degree() for c in L)
    return [[int(c.coeffs()[i]) if i <= c.degree() and c.degree() >= 0 else 0 for c in L] for i in range(D + 1)]


def strip_x(L):
    """divide out the common power of x (left multiplication: the solutions are unchanged)"""
    m = min(next(i for i, v in enumerate(c.coeffs()) if int(v)) for c in L if c.degree() >= 0)
    return [fmpz_poly([int(v) for v in c.coeffs()][m:]) if c.degree() >= 0 else c for c in L]


def t0_exponents(L):
    L = strip_x(L)
    q0 = sum(int(c.coeffs()[0]) * U ** k for k, c in enumerate(L) if c.degree() >= 0)
    out = []
    for r, m in sp.roots(sp.Poly(q0, U)).items():
        out += [Fr(int(r.p), int(r.q))] * m
    return sorted(out)


def at_infinity(L):
    D = max(c.degree() for c in L)
    return [fmpz_poly(([int(v) for v in c.coeffs()] + [0] * (D + 1 - len(c.coeffs())))[::-1]) * (-1) ** k
            for k, c in enumerate(L)]


# ---------------------------------------------------------------- (1), (2): operators on Y
def residual(L, t, perturb=False):
    mp.mp.dps = 80
    th = [mp.mpf(v) for v in DATA["thetaY"][t]]
    t2 = mp.mpf(t) ** 2
    D = max(c.degree() for c in L)
    mono = {}
    for k, c in enumerate(L):
        for i, v in enumerate(c.coeffs()):
            if int(v):
                mono[(k, i)] = mp.mpf(int(v)) * t2 ** i * th[k]
    if perturb:
        kmax = max(mono, key=lambda key: abs(mono[key]))
        mono[kmax] *= mp.mpf(201) / 200
    groups = [mp.fsum(v for (k, i), v in mono.items() if i == j) for j in range(D + 1)]
    return abs(mp.fsum(groups)) / max(abs(g) for g in groups)


@pytest.mark.parametrize("t", ["1.7", "2.6"])
def test_m13_annihilates_y(t):
    M13 = op("M13")
    assert len(M13) - 1 == 13
    assert residual(M13, t) < mp.mpf(10) ** -40
    assert residual(M13, t, perturb=True) > mp.mpf(10) ** -8


@pytest.mark.parametrize("t", ["1.7", "2.6"])
def test_y_needs_the_k3_piece(t):
    RS = op("R_S")
    assert len(RS) - 1 == 11
    assert residual(RS, t) > mp.mpf(10) ** -6


# ---------------------------------------------------------------- (3): left factors of M13
def strip_known(lead, p):
    e = lead
    facs = [nmod_poly([0, 1], p)] + [nmod_poly([(-c.numerator) % p, c.denominator % p], p) for c in CUTS]
    for f in facs:
        while (e % f).is_zero():
            e = e // f
    return e


def source_extra(S):
    cuts = [fmpz_poly([-c.numerator, c.denominator]) for c in CUTS]
    e = fmpz_poly([1])
    for f, _ in S[-1].factor()[1]:
        if f.degree() < 1 or f == fmpz_poly([0, 1]) or f == -fmpz_poly([0, 1]):
            continue
        if any(f == c or f == -c for c in cuts):
            continue
        e = e * f
    return e


def logfree_series(S, p, N):
    Q = qform(strip_x(S))
    for e in sorted(set(t0_exponents(S))):
        ex = e / 2
        ep = ex.numerator * pow(ex.denominator, -1, p) % p
        ev = lambda row, s: sum(c * pow(s, k, p) for k, c in enumerate(row)) % p
        c, ok = [1], True
        for n in range(1, N):
            acc = sum(ev(Q[j], 2 * (ep + n - j)) * c[n - j] for j in range(1, min(n, len(Q) - 1) + 1)) % p
            d = ev(Q[0], 2 * (ep + n))
            if d == 0:
                if acc:
                    ok = False
                    break
                c.append(0)
                continue
            c.append(-acc * pow(d, -1, p) % p)
        if ok:
            return ex, c
    raise RuntimeError("no log-free exponent")


def inv_series(d, N, p):
    return [int(v) for v in nmod_poly(d, p).inverse_series_trunc(N).coeffs()] + [0] * N


def hom_nullity(S, F, nu, p, den_extra, a=20, b=12, M=500):
    """dim of maps t^nu x^-a den^-1 sum_{j,m} p_jm x^m theta^j from Sol(S) into Sol(F) (F mod p)"""
    rS = len(S) - 1
    Ucols = rS * (M + 1)
    D = max(c.degree() for c in F)
    N = Ucols + D + 200
    ex, f = logfree_series(S, p, N)
    ET = 2 * ex.numerator * pow(ex.denominator, -1, p) % p
    den = nmod_poly([1], p)
    for c in CUTS:
        den = den * nmod_poly([(-c.numerator) % p, c.denominator % p], p) ** b
    den = den * nmod_poly([int(v) % p for v in source_extra(S).coeffs()], p) ** b * den_extra ** 2
    w = inv_series([int(v) for v in den.coeffs()], N, p)[:N]
    cols = np.zeros((N, Ucols), dtype=np.int64)
    nn = np.arange(N, dtype=np.int64)
    lam = (ET + 2 * nn) % p
    fa = np.array(f + [0] * (N - len(f)), dtype=np.int64) % p
    k = 0
    for j in range(rS):
        tj = fa.copy()
        for _ in range(j):
            tj = (tj * lam) % p
        g = [int(v) for v in (nmod_poly([int(v) for v in tj], p) * nmod_poly(w, p)).coeffs()][:N]
        gv = np.array(g + [0] * (N - len(g)), dtype=np.int64)
        for m in range(M + 1):
            cols[m:, k] = gv[: N - m]
            k += 1
    e2 = ((nu - 2 * a) % p + ET) % p
    Q = [[int(v) % p for v in row] for row in qform([fmpz_poly([int(v) for v in c.coeffs()]) for c in F])]
    Nout = N - len(Q)
    uu = (e2 + 2 * nn) % p
    V = np.zeros((len(Q), N), dtype=np.int64)
    for i, row in enumerate(Q):
        acc = np.zeros(N, dtype=np.int64)
        for c in reversed(row):
            acc = (acc * uu + c) % p
        V[i] = acc
    W = np.zeros((Nout, N), dtype=np.int64)
    for i in range(len(Q)):
        idx = np.arange(0, Nout - i)
        W[idx + i, idx] = V[i, idx]
    Pm = nmod_mat(Nout, N, [int(v) for v in W.ravel()], p) * nmod_mat(N, Ucols, [int(v) for v in cols.ravel()], p)
    return Ucols - Pm.rank()


@pytest.mark.slow
def test_k3_piece_is_the_only_left_factor_of_m13():
    p = P1
    M13 = modp(op("M13"), p)
    ext = strip_known(M13[-1], p)
    target = adjoint(M13)
    got = {}
    for key in ("R_V", "z1", "z2", "F_1", "F_2", "F_5"):
        Sa = adjoint(op(key))
        got[key] = [hom_nullity(Sa, target, nu, p, ext) for nu in (0, 1)]
    F5 = op("F_5")
    planted = adjoint(compose(modp(F5, p), M13))
    got["planted"] = [hom_nullity(adjoint(F5), planted, nu, p, ext) for nu in (0, 1)]
    assert got["planted"] == [1, 0]
    assert got["R_V"] == [0, 1]
    for key in ("z1", "z2", "F_1", "F_2", "F_5"):
        assert got[key] == [0, 0], (key, got[key])


# ---------------------------------------------------------------- (4): irreducible and orthogonal
NT, NTV = 100, 90


def fm(q, p):
    q = Fr(q)
    return q.numerator % p * pow(q.denominator % p, -1, p) % p


def taylor_basis(L, x0, p):
    """solutions at a regular point x0 as power series in h = x - x0, mod p (exact recursion in d/dh form)"""
    h = sp.Symbol("h")
    r = len(L) - 1
    dcoef = [sp.Integer(0)] * (r + 1)            # L = sum_j a_j(h) D^j
    for k, c in enumerate(L):
        poly = sum(int(v) * (x0 + h) ** i for i, v in enumerate(c.coeffs()))
        th = {0: sp.Integer(1)}                  # theta^k = sum_j b_j(h) D^j, theta = 2 (x0 + h) D
        for _ in range(k):
            nxt = {}
            for j, bj in th.items():
                nxt[j + 1] = nxt.get(j + 1, 0) + 2 * (x0 + h) * bj
                nxt[j] = nxt.get(j, 0) + 2 * (x0 + h) * sp.diff(bj, h)
            th = nxt
        for j, bj in th.items():
            dcoef[j] += poly * bj
    A = [[fm(v, p) for v in reversed(sp.Poly(sp.expand(a), h).all_coeffs())] if a != 0 else [0] for a in dcoef]
    sols = []
    for j0 in range(r):
        y = [0] * (NT + r)
        y[j0] = 1
        for n in range(NT + r - r):            # coefficient of h^n in L y determines y[n + r]
            acc = 0
            for j, aj in enumerate(A):
                for i, ai in enumerate(aj):
                    m = n - i + j
                    if 0 <= m < len(y) and ai and not (j == r and i == 0):
                        ff = 1
                        for q in range(j):
                            ff = ff * (m - q) % p
                        acc = (acc + ai * ff * y[m]) % p
            lead = A[r][0] * 1
            ff = 1
            for q in range(r):
                ff = ff * (n + r - q) % p
            if n + r < len(y):
                y[n + r] = (-acc) * pow(lead * ff % p, -1, p) % p
        sols.append(y[:NT])
    return sols


def mul(a, b, p):
    return ([int(v) for v in (nmod_poly(a, p) * nmod_poly(b, p)).coeffs()] + [0] * NT)[:NT]


def theta_s(s, x0, p):
    d = [(n + 1) * s[n + 1] % p for n in range(NT - 1)] + [0]
    x0p = fm(x0, p)
    return [(2 * (x0p * d[n] + (d[n - 1] if n else 0))) % p for n in range(NT)]


def rank(vecs, p):
    return nmod_mat(len(vecs), NTV, [v % p for vec in vecs for v in vec[:NTV]], p).rank() if vecs else 0


def binom_series(c, E, p):
    out, coef = [], 1
    Ep, cp = fm(E, p), fm(c, p)
    for n in range(NT):
        out.append(coef * pow(cp, n, p) % p)
        coef = coef * ((Ep - n) % p) % p * pow(n + 1, -1, p) % p
    return out


def class_extreme(vals, low=True):
    out = {}
    for e in vals:
        key = e - (e.numerator // e.denominator)
        if key not in out or (e < out[key] if low else e > out[key]):
            out[key] = e
    return out


def exponents(L):
    """x-exponents at 0, at infinity (growth) and at each finite root of the leading coefficient"""
    ex = {"0": [e / 2 for e in t0_exponents(L)], "inf": [e / 2 for e in t0_exponents(at_infinity(L))]}
    ex["inf"] = [-e for e in ex["inf"]]
    lead = sp.Poly([int(v) for v in reversed(L[-1].coeffs())], X)
    for rt in sp.roots(lead):
        if rt == 0:
            continue
        a = Fr(int(rt.p), int(rt.q))
        h = sp.Symbol("h")
        dco = {}
        for k, c in enumerate(L):           # indicial at x = a in d/dx form: theta = 2 x d/dx
            poly = sum(int(v) * (a + h) ** i for i, v in enumerate(c.coeffs()))
            th = {0: sp.Integer(1)}
            for _ in range(k):
                nxt = {}
                for j, bj in th.items():
                    nxt[j + 1] = nxt.get(j + 1, 0) + 2 * (a + h) * bj
                    nxt[j] = nxt.get(j, 0) + 2 * (a + h) * sp.diff(bj, h)
                th = nxt
            for j, bj in th.items():
                dco[j] = dco.get(j, 0) + sp.expand(poly * bj)
        ords = {j: sp.Poly(v, h).monoms()[-1][0] for j, v in dco.items() if v != 0}
        m = min(o - j for j, o in ords.items())
        s = sp.Symbol("s")
        ind = sum(sp.Poly(dco[j], h).coeff_monomial(h ** (m + j)) * sp.ff(s, j) for j in dco if dco[j] != 0 and m + j >= 0)
        rts = sp.roots(sp.Poly(sp.expand(ind), s))
        ex[a] = sorted(sum([[Fr(int(r.p), int(r.q))] * mm for r, mm in rts.items()], []))
    return ex


def shapes(ex, kind):
    def comb_(vals):
        if kind == "one":
            return list(vals)
        if kind == "sym":
            return [vals[i] + vals[j] for i in range(len(vals)) for j in range(i, len(vals))]
        return [vals[i] + vals[j] for i in range(len(vals)) for j in range(i + 1, len(vals))]
    f0s = sorted(class_extreme(comb_(ex["0"])).values())
    mus = class_extreme(comb_(ex["inf"]), low=False)
    fin = [a for a in ex if a not in ("0", "inf")]
    per = [sorted(class_extreme([v - (1 if kind == "wedge" else 0) for v in comb_(ex[a])]).values()) for a in fin]
    out = []
    for f0 in f0s:
        for cmb in itertools.product(*per):
            s = f0 + sum(cmb, Fr(0))
            for mu in mus.values():
                N = mu - s
                if N >= 0 and N.denominator == 1:
                    out.append((f0, dict(zip(fin, cmb)), int(N)))
    return out


def hyp_count(funcs, shape, x0, p):
    F0, E, N = shape
    pre = binom_series(Fr(1) / x0, F0, p)
    for a, Ea in E.items():
        pre = mul(pre, binom_series(Fr(1) / (x0 - a), Ea, p), p)
    cols = [([0] * k + pre)[:NT] for k in range(N + 1)]
    return (N + 1) - (rank(cols + funcs, p) - rank(funcs, p))


def squares(Y, x0, p):
    n = len(Y)
    W = [[(a - b) % p for a, b in zip(mul(Y[i], theta_s(Y[j], x0, p), p), mul(theta_s(Y[i], x0, p), Y[j], p))]
         for i in range(n) for j in range(i + 1, n)]
    S = [mul(Y[i], Y[j], p) for i in range(n) for j in range(i, n)]
    return W, S


def form_report(Y, ex, x0, p):
    W, S = squares(Y, x0, p)
    rep = {"wedge dep": len(W) - rank(W, p), "sym dep": len(S) - rank(S, p)}
    for kind, funcs in (("one", Y), ("wedge", W)):
        sh = shapes(ex, kind)
        rep[f"{kind} shapes"] = len(sh)
        rep[f"{kind} hyp"] = sum(hyp_count(funcs, s_, x0, p) for s_ in sh)
    return rep


def legendre_pair():
    return from_q([U ** 2, -(U + 1) ** 2]), from_q([9 * U ** 2, -(U + 1) ** 2])


def test_k3_piece_is_irreducible_and_orthogonal():
    p, x0 = P1, Fr(2)
    RV = op("R_V")
    rv = form_report(taylor_basis(RV, x0, p), exponents(RV), x0, p)
    rva = form_report(taylor_basis(adjoint(RV), x0, p), exponents(adjoint(RV)), x0, p)
    assert rv["one hyp"] == 0 and rva["one hyp"] == 0            # no order-one or order-three right factor
    assert rv["wedge shapes"] >= 1 and rv["wedge hyp"] == 0       # no order-two factor
    assert rv["sym dep"] == 1 and rv["wedge dep"] == 0            # symmetric form, no alternating one
    # the same code on a direct sum of Legendre operators finds both Wronskians (a planted order-two factor)
    A, B = legendre_pair()
    Ys = taylor_basis(A, x0, p) + taylor_basis(B, x0, p)
    h = Fr(-1, 2)
    exs = {"0": [Fr(0)] * 4, "inf": [h] * 4, Fr(1): [Fr(0), Fr(0), Fr(0), Fr(1)], Fr(9): [Fr(0), Fr(0), Fr(0), Fr(1)]}
    assert form_report(Ys, exs, x0, p)["wedge hyp"] == 2
    # and on AESZ #1 the Wronskian relation of a symplectic operator
    Q5 = from_q([U ** 4, sp.expand(-5 * (5 * U + 2) * (5 * U + 4) * (5 * U + 6) * (5 * U + 8))])
    q5 = form_report(taylor_basis(Q5, x0, p), exponents(Q5), x0, p)
    assert q5["wedge dep"] == 1 and q5["sym dep"] == 0


# ---------------------------------------------------------------- (5): Jordan blocks
def jordan_blocks(L, where, p=P1):
    Lw = L if where == "0" else at_infinity(L)
    Q = qform(Lw)
    rts = t0_exponents(Lw)
    r = len(L) - 1
    cls = {}
    for s in rts:
        cls.setdefault(s - 2 * (s.numerator // (2 * s.denominator)), []).append(s)
    out = {}
    for key, lst in cls.items():
        s0 = min(lst)
        N = int((max(lst) - s0) / 2) + 25
        d = []
        for K in range(r):
            n = K + 1
            rows = []
            for ne in range(N + 1):
                block = [[0] * ((N + 1) * n) for _ in range(n)]
                for j in range(min(ne, len(Q) - 1) + 1):
                    a = (fm(s0, p) + 2 * (ne - j)) % p
                    Amat = [[(a if i == m else ((i + 1) if m == i + 1 else 0)) % p for m in range(n)] for i in range(n)]
                    R = [[0] * n for _ in range(n)]
                    for c in reversed(Q[j]):
                        R = [[sum(R[i][m] * Amat[m][jj] for m in range(n)) % p for jj in range(n)] for i in range(n)]
                        for i in range(n):
                            R[i][i] = (R[i][i] + c) % p
                    for i in range(n):
                        for m in range(n):
                            block[i][(ne - j) * n + m] = (block[i][(ne - j) * n + m] + R[i][m]) % p
                rows += block
            Mx = nmod_mat(len(rows), (N + 1) * n, [v for row in rows for v in row], p)
            d.append((N + 1) * n - Mx.rank())
        blocks = []
        for K in range(r):
            more = d[K] - (d[K - 1] if K else 0)
            nxt = (d[K + 1] - d[K]) if K + 1 < r else 0
            blocks += [K + 1] * (more - nxt)
        out[key] = sorted(blocks, reverse=True)
    return out


def test_k3_piece_jordan_blocks():
    T = sp.Symbol("T")
    c22 = from_q([sp.expand((T ** 2 * (T - 1) ** 2).subs(T, U / 2) * 16)])
    c31 = from_q([sp.expand((T ** 3 * (T - 1)).subs(T, U / 2) * 16)])
    assert list(jordan_blocks(c22, "0").values()) == [[2, 2]]
    assert list(jordan_blocks(c31, "0").values()) == [[3, 1]]
    RV = op("R_V")
    assert list(jordan_blocks(RV, "0").values()) == [[2, 2]]
    assert list(jordan_blocks(RV, "inf").values()) == [[3, 1]]
