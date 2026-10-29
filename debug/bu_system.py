"""Record BU: the first-order system in L for case N1, built from the IBP reductions of debug/bu_ibp.py and the
three-term rules of i_L, k_L.  No coefficient data enters the construction.

Coordinates U_L (18, all formal: transcendental constants travel inside the coordinates, so every relation has
coefficients in Q(L)):
    C1 = D2^{-L} mK (x) Phi1        Phi1 = (int_0^D1 u_L, int_0^D1 u_{L+1}),  family I            [4]
    C2 = D2^{-L} mI (x) Phi2        Phi2 = (int_D1^D2 u_L, int_D1^D2 u_{L+1}), family K          [4]
    C3 = D2^{L}  mI  phi3           phi3 = int_D2^oo u_{L+1}, family H                            [2]
    E1 = e^{beta D1} (D1/D2)^L mK (x) mI     (boundary values at r = D1)                          [4]
    E2 = e^{beta D2} mI (x) nK                (boundary values at r = D2)                          [4]
mK = (k_L(y), k_{L+1}(y)), mI = (i_L(y), i_{L+1}(y)), nK = (k_L(w), k_{L+1}(w)); y = zH D1, w = zH D2.
U_{L+1} = A(L) U_L and F_L / (2 pi N) = c(L) . U_L.
"""
import sympy as sp

import bu_ibp as ib

L = ib.L
D1, D2 = sp.Rational(9, 5), sp.Rational(13, 5)
zX, zH = sp.Rational(13, 10), sp.Integer(1)
y, w = zH * D1, zH * D2
assert ib.beta == -zX and ib.z == zH

# index maps
IDX = {}
for blk, n in [('C1', 4), ('C2', 4), ('C3', 2), ('E1', 4), ('E2', 4)]:
    for t in range(n):
        IDX[(blk, t)] = len(IDX)
N = len(IDX)


def pair(a, b):
    """flattened index of a 2x2 tensor coordinate (a = first factor index, b = second)."""
    return 2 * a + b


def bessel_shift(eps, t):
    """2x2 shift of (Y_L, Y_{L+1}) at argument t: Y_{L+2} = Y_L - eps (2L+3)/t Y_{L+1} (eps = +1 for i, -1 for k)."""
    return sp.Matrix([[0, 1], [1, -eps * (2 * L + 3) / t]])


def build():
    A = sp.zeros(N, N)
    c = sp.zeros(1, N)
    AK, AI, AN = bessel_shift(-1, y), bessel_shift(1, y), bessel_shift(-1, w)
    rel = {f: ib.relations(f) for f in ('I', 'K', 'H')}

    def lin(expr):
        return ib.coeffs(expr)

    # ---- E1, E2 blocks
    for a in range(2):
        for b in range(2):
            for a2 in range(2):
                for b2 in range(2):
                    A[IDX[('E1', pair(a, b))], IDX[('E1', pair(a2, b2))]] += (D1 / D2) * AK[a, a2] * AI[b, b2]
                    A[IDX[('E2', pair(a, b))], IDX[('E2', pair(a2, b2))]] += AI[a, a2] * AN[b, b2]

    # ---- C1: D2^{-L} mK (x) Phi1, family I on [0, D1]; B(k,j) = D1^k u_j(D1) (the r = 0 end vanishes for k >= 0)
    # u_j(L)(D1) = e^{beta D1} D1^{L+2+j} i_{L+j}(y)
    for s, expr in enumerate(rel['I']['next']):
        for nm, cf in lin(expr).items():
            kind, k, j = ib.parse(nm)
            for a in range(2):
                for a2 in range(2):
                    if kind == 'J':
                        m = ib.MASTERS['I'].index((k, j))
                        A[IDX[('C1', pair(a, s))], IDX[('C1', pair(a2, m))]] += cf * AK[a, a2] / D2
                    else:
                        assert k >= 0
                        A[IDX[('C1', pair(a, s))], IDX[('E1', pair(a2, j))]] += cf * AK[a, a2] * D1 ** (2 + j + k) / D2
    # ---- C2: D2^{-L} mI (x) Phi2, family K on [D1, D2]
    # B(k,j) = D2^k e^{beta D2} D2^{L+2+j} k_{L+j}(w) - D1^k e^{beta D1} D1^{L+2+j} k_{L+j}(y)
    for s, expr in enumerate(rel['K']['next']):
        for nm, cf in lin(expr).items():
            kind, k, j = ib.parse(nm)
            for a in range(2):
                for a2 in range(2):
                    if kind == 'J':
                        m = ib.MASTERS['K'].index((k, j))
                        A[IDX[('C2', pair(a, s))], IDX[('C2', pair(a2, m))]] += cf * AI[a, a2] / D2
                    else:
                        A[IDX[('C2', pair(a, s))], IDX[('E2', pair(a2, j))]] += cf * AI[a, a2] * D2 ** (k + 1 + j)
                        A[IDX[('C2', pair(a, s))], IDX[('E1', pair(j, a2))]] += -cf * AI[a, a2] * D1 ** (k + 2 + j) / D2
    # ---- C3: D2^L mI phi3, family H on [D2, oo); phi3 = J(0,1) (first master); B(k,j) = -D2^k e^{beta D2}
    # D2^{1-L-j} k_{L+j}(w)
    nxt = lin(rel['H']['next'][0])
    assert set(n for n in nxt if n.startswith('J')) == {'J_0_1'}
    for nm, cf in nxt.items():
        kind, k, j = ib.parse(nm)
        for a in range(2):
            for a2 in range(2):
                if kind == 'J':
                    A[IDX[('C3', a)], IDX[('C3', a2)]] += cf * AI[a, a2] * D2
                else:
                    A[IDX[('C3', a)], IDX[('E2', pair(a2, j))]] += -cf * AI[a, a2] * D2 ** (k + 2 - j)

    # ---- output
    c[IDX[('C1', pair(0, 0))]] += -2 / D2 * (2 * L + 1)
    c[IDX[('C1', pair(1, 0))]] += -2 / D2 * (-y)
    c[IDX[('C1', pair(0, 1))]] += -2 / D2 * zH
    c[IDX[('C2', pair(0, 0))]] += -2 / D2 * (2 * L + 1)
    c[IDX[('C2', pair(1, 0))]] += -2 / D2 * y
    c[IDX[('C2', pair(0, 1))]] += -2 / D2 * (-zH)
    h_int, r2h = (lin(e) for e in rel['H']['integrand'])
    for weight_i, lc in [((2 * L + 1, y), h_int), ((-zH, 0), r2h)]:
        for nm, cf in lc.items():
            kind, k, j = ib.parse(nm)
            for a in range(2):
                wa = weight_i[a]
                if wa == 0:
                    continue
                if kind == 'J':
                    assert (k, j) == (0, 1)
                    c[IDX[('C3', a)]] += -2 * wa * cf
                else:
                    c[IDX[('E2', pair(a, j))]] += -2 * wa * cf * (-D2 ** (k + 1 - j))
    A = A.applyfunc(lambda e: sp.factor(sp.together(e)))
    c = c.applyfunc(lambda e: sp.factor(sp.together(e)))
    return A, c


def numeric_U(Lv, dps=50):
    """U_L by mpmath quadrature and Bessel values (real numbers; for checking only)."""
    import mpmath as mp
    mp.mp.dps = dps
    bt = -mp.mpf(13) / 10
    d1, d2, yy, ww = mp.mpf(9) / 5, mp.mpf(13) / 5, mp.mpf(9) / 5, mp.mpf(13) / 5

    def iL(n, t):
        return mp.sqrt(mp.pi / (2 * t)) * mp.besseli(n + mp.mpf(1) / 2, t)

    def kL(n, t):
        return mp.sqrt(2 / (mp.pi * t)) * mp.besselk(n + mp.mpf(1) / 2, t)

    uI = lambda n, x: mp.e ** (bt * x) * x ** (n + 2) * iL(n, x)
    uK = lambda n, x: mp.e ** (bt * x) * x ** (n + 2) * kL(n, x)
    uH = lambda n, x: mp.e ** (bt * x) * x ** (1 - n) * kL(n, x)
    mK = [kL(Lv, yy), kL(Lv + 1, yy)]
    mI = [iL(Lv, yy), iL(Lv + 1, yy)]
    nK = [kL(Lv, ww), kL(Lv + 1, ww)]
    P1 = [mp.quad(lambda x: uI(Lv + j, x), [0, d1]) for j in (0, 1)]
    P2 = [mp.quad(lambda x: uK(Lv + j, x), [d1, d2]) for j in (0, 1)]
    p3 = mp.quad(lambda x: uH(Lv + 1, x), [d2, mp.inf])
    U = [mp.mpf(0)] * N
    for a in range(2):
        for b in range(2):
            U[IDX[('C1', pair(a, b))]] = d2 ** (-Lv) * mK[a] * P1[b]
            U[IDX[('C2', pair(a, b))]] = d2 ** (-Lv) * mI[a] * P2[b]
            U[IDX[('E1', pair(a, b))]] = mp.e ** (bt * d1) * (d1 / d2) ** Lv * mK[a] * mI[b]
            U[IDX[('E2', pair(a, b))]] = mp.e ** (bt * d2) * mI[a] * nK[b]
        U[IDX[('C3', a)]] = d2 ** Lv * mI[a] * p3
    return U


def numeric_F(Lv, dps=50):
    """F_L / (2 pi N) directly by quadrature of the defining radial integral (checking only)."""
    import mpmath as mp
    mp.mp.dps = dps
    bt = -mp.mpf(13) / 10
    d1, d2 = mp.mpf(9) / 5, mp.mpf(13) / 5

    def iL(n, t):
        return mp.sqrt(mp.pi / (2 * t)) * mp.besseli(n + mp.mpf(1) / 2, t)

    def kL(n, t):
        return mp.sqrt(2 / (mp.pi * t)) * mp.besselk(n + mp.mpf(1) / 2, t)

    def G(x):
        a, b = (x, d1) if x < d1 else (d1, x)
        return -2 * ((2 * Lv + 1) * iL(Lv, a) * kL(Lv, b) + a * iL(Lv + 1, a) * kL(Lv, b) - b * iL(Lv, a) * kL(Lv + 1, b))

    def wl(x):
        return x ** Lv / d2 ** (Lv + 1) if x < d2 else d2 ** Lv / x ** (Lv + 1)

    return mp.quad(lambda x: x ** 2 * mp.e ** (bt * x) * G(x) * wl(x), [0, d1, d2, mp.inf])


if __name__ == '__main__':
    import mpmath as mp
    A, c = build()
    worst_sys, worst_out = 0, 0
    for Lv in range(0, 6):
        U0, U1 = numeric_U(Lv), numeric_U(Lv + 1)
        An = [[mp.mpf(sp.Rational(A[i, j].subs(L, Lv)).p) / sp.Rational(A[i, j].subs(L, Lv)).q for j in range(N)]
              for i in range(N)]
        for i in range(N):
            v = sum(An[i][j] * U0[j] for j in range(N))
            worst_sys = max(worst_sys, abs(v - U1[i]) / max(abs(U1[i]), mp.mpf(10) ** -40))
        cv = sum(mp.mpf(sp.Rational(c[0, j].subs(L, Lv)).p) / sp.Rational(c[0, j].subs(L, Lv)).q * U0[j]
                 for j in range(N))
        F = numeric_F(Lv)
        worst_out = max(worst_out, abs(cv - F) / abs(F))
        print(Lv, mp.nstr(cv, 20), mp.nstr(F, 20), flush=True)
    print('system check worst rel', mp.nstr(worst_sys, 3), ' output check worst rel', mp.nstr(worst_out, 3))


# ---------------------------------------------------------------------------------------------------------------
# Rank-one version (bu_ibp.relations1): one master per piece, phi_P = int u_{L+1}, first-order recurrence.
# Coordinates (14): C1 = D2^{-L} mK phi1 [2], C2 = D2^{-L} mI phi2 [2], C3 = D2^{L} mI phi3 [2], E1 [4], E2 [4].
IDX1 = {}
for _blk, _n in [('C1', 2), ('C2', 2), ('C3', 2), ('E1', 4), ('E2', 4)]:
    for _t in range(_n):
        IDX1[(_blk, _t)] = len(IDX1)
N1 = len(IDX1)


def _bd_targets(piece, k, j, a2):
    """g^L m[a2] r0^k u_j(L)(r0), summed over the piece's ends with signs, as [(coef, coordinate index)]."""
    if piece == 'C1':      # family I on [0, D1]; the r = 0 end vanishes (u_j ~ r^{2L+2+j}, k >= -1)
        assert k >= -1
        return [(D1 ** (2 + j + k), IDX1[('E1', pair(a2, j))])]
    if piece == 'C2':      # family K on [D1, D2]
        return [(D2 ** (2 + j + k), IDX1[('E2', pair(a2, j))]), (-D1 ** (2 + j + k), IDX1[('E1', pair(j, a2))])]
    if piece == 'C3':      # family H on [D2, oo)
        return [(-D2 ** (1 - j + k), IDX1[('E2', pair(a2, j))])]
    raise KeyError(piece)


def build1():
    A = sp.zeros(N1, N1)
    c = sp.zeros(1, N1)
    AK, AI, AN = bessel_shift(-1, y), bessel_shift(1, y), bessel_shift(-1, w)
    for a in range(2):
        for b in range(2):
            for a2 in range(2):
                for b2 in range(2):
                    A[IDX1[('E1', pair(a, b))], IDX1[('E1', pair(a2, b2))]] += (D1 / D2) * AK[a, a2] * AI[b, b2]
                    A[IDX1[('E2', pair(a, b))], IDX1[('E2', pair(a2, b2))]] += AI[a, a2] * AN[b, b2]
    pieces = {'C1': ('I', AK, 1 / D2, -2 / D2, [(2 * L + 1, -y), (zH, 0)]),
              'C2': ('K', AI, 1 / D2, -2 / D2, [(2 * L + 1, y), (-zH, 0)]),
              'C3': ('H', AI, D2, sp.Integer(-2), [(2 * L + 1, y), (-zH, 0)])}
    for piece, (fam, AM, g, pref, weights) in pieces.items():
        rel = ib.relations1(fam)
        for nm, cf in ib.coeffs(rel['next'][0]).items():
            kind, k, j = ib.parse(nm)
            for a in range(2):
                for a2 in range(2):
                    if kind == 'J':
                        assert (k, j) == (0, 1)
                        A[IDX1[(piece, a)], IDX1[(piece, a2)]] += g * cf * AM[a, a2]
                    else:
                        for co, t in _bd_targets(piece, k, j, a2):
                            A[IDX1[(piece, a)], t] += g * cf * AM[a, a2] * co
        for wts, expr in zip(weights, rel['integrand']):
            for nm, cf in ib.coeffs(expr).items():
                kind, k, j = ib.parse(nm)
                for a in range(2):
                    if wts[a] == 0:
                        continue
                    if kind == 'J':
                        assert (k, j) == (0, 1)
                        c[IDX1[(piece, a)]] += pref * wts[a] * cf
                    else:
                        for co, t in _bd_targets(piece, k, j, a):
                            c[t] += pref * wts[a] * cf * co
    A = A.applyfunc(lambda e: sp.factor(sp.together(e)))
    c = c.applyfunc(lambda e: sp.factor(sp.together(e)))
    return A, c


def numeric_U1(Lv, dps=50):
    import mpmath as mp
    mp.mp.dps = dps
    U = numeric_U(Lv, dps)
    bt = -mp.mpf(13) / 10
    d1, d2, yy = mp.mpf(9) / 5, mp.mpf(13) / 5, mp.mpf(9) / 5

    def iL(n, t):
        return mp.sqrt(mp.pi / (2 * t)) * mp.besseli(n + mp.mpf(1) / 2, t)

    def kL(n, t):
        return mp.sqrt(2 / (mp.pi * t)) * mp.besselk(n + mp.mpf(1) / 2, t)
    mK = [kL(Lv, yy), kL(Lv + 1, yy)]
    mI = [iL(Lv, yy), iL(Lv + 1, yy)]
    p1 = mp.quad(lambda x: mp.e ** (bt * x) * x ** (Lv + 3) * iL(Lv + 1, x), [0, d1])
    p2 = mp.quad(lambda x: mp.e ** (bt * x) * x ** (Lv + 3) * kL(Lv + 1, x), [d1, d2])
    p3 = mp.quad(lambda x: mp.e ** (bt * x) * x ** (-Lv) * kL(Lv + 1, x), [d2, mp.inf])
    V = [mp.mpf(0)] * N1
    for a in range(2):
        V[IDX1[('C1', a)]] = d2 ** (-Lv) * mK[a] * p1
        V[IDX1[('C2', a)]] = d2 ** (-Lv) * mI[a] * p2
        V[IDX1[('C3', a)]] = d2 ** Lv * mI[a] * p3
    for t in range(4):
        V[IDX1[('E1', t)]] = U[IDX[('E1', t)]]
        V[IDX1[('E2', t)]] = U[IDX[('E2', t)]]
    return V


def check1():
    import mpmath as mp
    A, c = build1()
    ws, wo = 0, 0
    for Lv in range(0, 6):
        U0, U1 = numeric_U1(Lv), numeric_U1(Lv + 1)

        def q(e):
            v = sp.Rational(e.subs(L, Lv))
            return mp.mpf(v.p) / v.q
        for i in range(N1):
            v = sum(q(A[i, j]) * U0[j] for j in range(N1))
            ws = max(ws, abs(v - U1[i]) / abs(U1[i]))
        cv = sum(q(c[0, j]) * U0[j] for j in range(N1))
        wo = max(wo, abs(cv - numeric_F(Lv)) / abs(cv))
    return ws, wo
