"""Record B2 (docs/predictions/2026-10-28_two_electron_bottom_up.md): the two-electron three-centre integral (XY|XZ),
case W1, as a first-order system in L built from the label rules.  No coefficient data enters this file.

Integrand.  With g_j(r) = r^2 e^{-zX r} G_L^{(D_j)}(r) (debug/ao_exact.py),
    F = int int g1 g2 r<^L / r>^(L+1) = I_C + I_T,
    I_C = int_0^oo dt a1(t) Phi2(t),  I_T = int_0^oo dt a2(t) Phi1(t),  a_j = g_j t^{-L-1},  Phi_j(t) = int_0^t g_j s^L.
On a piece g_j s^L = -2 (M0 u_L + zH M1 u_{L+1}) and g_j t^{-L-1} = -2 (M0 v_L + zH M1 t^2 v_{L+1}), with
    u_L = e^{beta r} r^{L+2} Y_L(z r)   (module U[Y]),   v_L = e^{beta r} r^{1-L} Y_L(z r)   (module V[Y]),
Y = i below D_j (multipliers M0 = (2L+1) k_L(y_j) - y_j k_{L+1}(y_j), M1 = k_L(y_j)) and Y = k above
(M0 = (2L+1) i_L(y_j) + y_j i_{L+1}(y_j), M1 = -i_L(y_j)), y_j = zH D_j, beta = -zX, z = zH.
Phi_j(t) = int_0^t (inner form) for t < D_j and Phi_j(t) = Phi_j(oo) - int_t^oo (outer form) for t > D_j, so the
nested integrals are integrals of the module N[VY, UY, side] = span{v_c psi, v_c u_d}, psi = int_0^t u_{L+1}
(side 0) or int_t^oo u_{L+1} (side oo).

Machinery (all exact mod p at a point L = lam of F_p; nothing is fitted):
  - each module: basis, d/dr matrix and L-shift matrix with Laurent-in-r entries;
  - IBP:  int d/dr(r^k e_i) = [r^k e_i]  reduces integrals to masters + boundary symbols B(k, i) (row echelon);
  - atoms: Bessel pairs (Y_L, Y_{L+1})(z t0) at t0 in {D1, D2}; integral masters (module, interval, master);
  - coordinates: monomials e^{tag} D1^{a L} D2^{b L} prod(atoms); the shift of a monomial is the product of the
    atoms' shifts (integral atoms: reduced master shift plus boundary values at the interval ends); closure from the
    output expression gives U_{L+1} = A(L) U_L, F_L = c(L) . U_L.
"""
from fractions import Fraction

import flint

D1, D2 = Fraction(17, 10), Fraction(5, 2)
ZX, ZH = Fraction(7, 5), Fraction(11, 10)
BETA, Z = -ZX, ZH
INF = None


class Ctx:
    """Field + caches for one evaluation point lam."""

    def __init__(self, p, lam):
        self.p = p
        self.L = flint.nmod(int(lam), p)
        self.red = {}

    def F(self, q):
        q = Fraction(q)
        return flint.nmod(q.numerator, self.p) / flint.nmod(q.denominator, self.p)


# ----------------------------------------------------------------------------------------------------------- modules
# entries: {i: [(j, power of r, fn(ctx) -> nmod)]}

def _eps(Y):
    return 1 if Y == 'i' else -1


def mod_U(Y):
    e = _eps(Y)
    D = {0: [(0, 0, lambda c: c.F(BETA)), (0, -1, lambda c: 2 * c.L + 2), (1, -1, lambda c: c.F(e * Z))],
         1: [(1, 0, lambda c: c.F(BETA)), (1, -1, lambda c: c.F(1)), (0, 1, lambda c: c.F(e * Z))]}
    S = {0: [(1, 0, lambda c: c.F(1))],
         1: [(0, 2, lambda c: c.F(1)), (1, 0, lambda c: -c.F(e) * (2 * c.L + 3) / c.F(Z))]}
    return dict(name=f'U{Y}', basis=['u0', 'u1'], D=D, S=S, kmin=-1, kmax=5, prefer=[(0, 1), (0, 0)])


def mod_V(Y):
    e = _eps(Y)
    D = {0: [(0, 0, lambda c: c.F(BETA)), (0, -1, lambda c: c.F(1)), (1, 1, lambda c: c.F(e * Z))],
         1: [(1, 0, lambda c: c.F(BETA)), (1, -1, lambda c: -(2 * c.L + 2)), (0, -1, lambda c: c.F(e * Z))]}
    S = {0: [(1, 0, lambda c: c.F(1))],
         1: [(0, -2, lambda c: c.F(1)), (1, -2, lambda c: -c.F(e) * (2 * c.L + 3) / c.F(Z))]}
    return dict(name=f'V{Y}', basis=['v0', 'v1'], D=D, S=S, kmin=-5, kmax=5, prefer=[(0, 1), (0, 0)])


MODS = {}


def get_mod(name):
    if name in MODS:
        return MODS[name]
    if name[0] == 'U':
        m = mod_U(name[1])
    elif name[0] == 'V':
        m = mod_V(name[1])
    else:                                   # N<VY><UY><side>, side '0' or 'oo'
        m = mod_N(name[1], name[2], name[3:])
    MODS[name] = m
    return m


def u_shift_rank_one(ctx, UY):
    """psi(L+1) = rho psi(L) + sum gamma_kj [r^k u_j]: the reduced shift of the U master J(0,1)."""
    res = reduce_shift(ctx, f'U{UY}', (0, 1))
    masters, bds = res
    assert set(masters) <= {(0, 1)}, masters
    return masters.get((0, 1), ctx.F(0)), bds


def mod_N(VY, UY, side):
    V, U = get_mod(f'V{VY}'), get_mod(f'U{UY}')
    sig = 1 if side == '0' else -1
    basis = ['n0', 'n1', 'p00', 'p01', 'p10', 'p11']
    ix = {b: t for t, b in enumerate(basis)}
    D = {t: [] for t in range(6)}
    for c in (0, 1):
        for (c2, pw, fn) in V['D'][c]:
            D[ix[f'n{c}']].append((ix[f'n{c2}'], pw, fn))
        D[ix[f'n{c}']].append((ix[f'p{c}1'], 0, (lambda c_: (lambda ctx: ctx.F(sig)))(c)))
        for d in (0, 1):
            for (c2, pw, fn) in V['D'][c]:
                D[ix[f'p{c}{d}']].append((ix[f'p{c2}{d}'], pw, fn))
            for (d2, pw, fn) in U['D'][d]:
                D[ix[f'p{c}{d}']].append((ix[f'p{c}{d2}'], pw, fn))
    S = {t: [] for t in range(6)}
    for c in (0, 1):
        for (c2, pw, fn) in V['S'][c]:
            # v_c(L+1) psi(L+1) = sum Sv (rho n_c2 + sig' sum gamma r^k p_c2j)
            def f_rho(ctx, fn=fn):
                return fn(ctx) * u_shift_rank_one(ctx, UY)[0]
            S[ix[f'n{c}']].append((ix[f'n{c2}'], pw, f_rho))
            for k in range(-6, 8):
                for j in (0, 1):
                    def f_g(ctx, fn=fn, k=k, j=j):
                        g = u_shift_rank_one(ctx, UY)[1].get((k, j))
                        return ctx.F(0) if g is None else fn(ctx) * g * ctx.F(sig)
                    S[ix[f'n{c}']].append((ix[f'p{c2}{j}'], pw + k, f_g))
        for d in (0, 1):
            for (c2, pw, fn) in V['S'][c]:
                for (d2, pw2, fn2) in U['S'][d]:
                    S[ix[f'p{c}{d}']].append((ix[f'p{c2}{d2}'], pw + pw2, (lambda a, b: (lambda ctx: a(ctx) * b(ctx)))(fn, fn2)))
    # side 0 is integrated from r = 0, where [r^k e] vanishes for k > -(2L+3): kmin = -6 is valid for L >= 2
    # (the derived recurrence is therefore claimed for L >= 2 and checked on the data from L = 0, see b2_check)
    kmin = -6
    return dict(name=f'N{VY}{UY}{side}', basis=basis, D=D, S=S, kmin=kmin, kmax=6,
                prefer=[(0, t) for t in range(6)])


# ------------------------------------------------------------------------------------------------------------- IBP

def ibp(ctx, name):
    """Row-reduced IBP system of a module at ctx.L: returns (masters, solve) with solve((k, i)) ->
    ({master (k, i): coef}, {(k, i): coef of B(k, i)})."""
    key = ('ibp', name)
    if key in ctx.red:
        return ctx.red[key]
    m = get_mod(name)
    n = len(m['basis'])
    ks = range(m['kmin'], m['kmax'] + 1)
    eqs = []
    unk = set()
    for k in ks:
        for i in range(n):
            row = {}
            if k != 0:
                row[(k - 1, i)] = row.get((k - 1, i), ctx.F(0)) + ctx.F(k)
            for (j, pw, fn) in m['D'][i]:
                v = fn(ctx)
                if int(v) != 0:
                    row[(k + pw, j)] = row.get((k + pw, j), ctx.F(0)) + v
            eqs.append((row, (k, i)))
            unk |= set(row)
    pref = [t for t in m['prefer'] if t in unk]
    others = sorted(unk - set(pref), key=lambda t: (-abs(t[0]), t[0], t[1]))
    cols = others + pref[::-1]
    bcols = [b for _, b in eqs]
    ci = {c: t for t, c in enumerate(cols)}
    bi = {b: len(cols) + t for t, b in enumerate(bcols)}
    ncol = len(cols) + len(bcols)
    ent = [0] * (len(eqs) * ncol)
    for r, (row, b) in enumerate(eqs):
        for c, v in row.items():
            ent[r * ncol + ci[c]] = int(v)
        ent[r * ncol + bi[b]] = int(-ctx.F(1))
    M = flint.nmod_mat(len(eqs), ncol, ent, ctx.p)
    R, rank = M.rref()
    pivots = {}
    for r in range(rank):
        for c in range(ncol):
            if int(R[r, c]) != 0:
                pivots[c] = r
                break
    free = [c for c in range(len(cols)) if c not in pivots]
    masters = [cols[c] for c in free]

    def solve(t):
        if t not in ci:
            raise KeyError(('outside window', name, t))
        c = ci[t]
        if c not in pivots:
            return {t: ctx.F(1)}, {}
        r = pivots[c]
        ms, bs = {}, {}
        for c2 in free:
            v = R[r, c2]
            if int(v) != 0:
                ms[cols[c2]] = -v
        for b, c2 in bi.items():
            v = R[r, c2]
            if int(v) != 0:
                bs[b] = -v
        return ms, bs
    ctx.red[key] = (masters, solve)
    return masters, solve


EDGE = 2   # a master within EDGE of the window edge is treated as a window artifact


def reduce(ctx, name, t):
    masters, solve = ibp(ctx, name)
    ms, bs = solve(t)
    m = get_mod(name)
    for (k, i) in ms:
        if k <= m['kmin'] - 1 or k >= m['kmax'] - 1:
            raise RuntimeError(f'{name}: target {t} depends on edge master {(k, i)}')
    return ms, bs


def reduce_shift(ctx, name, master):
    """int r^k e_i(L+1) for master (k, i), reduced at L."""
    key = ('shift', name, master)
    if key in ctx.red:
        return ctx.red[key]
    m = get_mod(name)
    k0, i0 = master
    ms_tot, bs_tot = {}, {}
    for (j, pw, fn) in m['S'][i0]:
        v = fn(ctx)
        if int(v) == 0:
            continue
        ms, bs = reduce(ctx, name, (k0 + pw, j))
        for a, b in ms.items():
            ms_tot[a] = ms_tot.get(a, ctx.F(0)) + v * b
        for a, b in bs.items():
            bs_tot[a] = bs_tot.get(a, ctx.F(0)) + v * b
    out = ({a: b for a, b in ms_tot.items() if int(b)}, {a: b for a, b in bs_tot.items() if int(b)})
    ctx.red[key] = out
    return out


# ------------------------------------------------------------------------------------------------- monomial algebra
# monomial: (tag, g1, g2, factors) ; factor = (atom, comp); atom = ('P', Y, t0) or ('J', modname, A, B, master)

ONE = (Fraction(0), 0, 0, ())


def mmul(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2], tuple(sorted(a[3] + b[3], key=repr)))


def emul(x, y):
    out = {}
    for a, u in x.items():
        for b, v in y.items():
            k = mmul(a, b)
            out[k] = out.get(k, 0) + u * v
    return out


def eadd(x, y, s=1):
    out = dict(x)
    for k, v in y.items():
        out[k] = out.get(k, 0) + s * v
    return out


def escale(x, s):
    return {k: v * s for k, v in x.items()}


def geom(t0, n):
    """D^{n L} for t0 in {D1, D2}."""
    return (Fraction(0), n if t0 == D1 else 0, n if t0 == D2 else 0, ())


def point_value(ctx, name, i, t0, k):
    """t0^k e_i(L)(t0) as an expression; t0 in {D1, D2}; 0 and infinity give 0 (windows chosen so)."""
    if t0 is INF or t0 == 0:
        return {}
    F = ctx.F
    m = get_mod(name)
    b = m['basis'][i]
    if name[0] == 'U':
        d = int(b[1])
        mono = mmul((BETA * t0, 0, 0, (((('P', name[1], Z * t0)), d),)), geom(t0, 1))
        return {mono: F(t0) ** (2 + d + k) if 2 + d + k >= 0 else F(1) / F(t0) ** (-(2 + d + k))}
    if name[0] == 'V':
        c = int(b[1])
        mono = mmul((BETA * t0, 0, 0, ((('P', name[1], Z * t0), c),)), geom(t0, -1))
        return {mono: F(t0 ** (1 - c + k))}
    VY, UY, side = name[1], name[2], name[3:]
    if b[0] == 'n':
        c = int(b[1])
        v = point_value(ctx, f'V{VY}', c, t0, k)
        psi_atom = ('J', f'U{UY}', 0, t0, (0, 1)) if side == '0' else ('J', f'U{UY}', t0, INF, (0, 1))
        return emul(v, {(Fraction(0), 0, 0, ((psi_atom, 0),)): F(1)})
    c, d = int(b[1]), int(b[2])
    return emul(point_value(ctx, f'V{VY}', c, t0, k), point_value(ctx, f'U{UY}', d, t0, 0))


def integral_expr(ctx, name, A, B, t):
    """int_A^B r^k e_i for t = (k, i), as an expression in integral atoms and boundary monomials."""
    ms, bs = reduce(ctx, name, t)
    out = {}
    for mk, v in ms.items():
        out = eadd(out, {(Fraction(0), 0, 0, ((('J', name, A, B, mk), 0),)): v})
    for (k, i), v in bs.items():
        out = eadd(out, escale(point_value(ctx, name, i, B, k), v))
        out = eadd(out, escale(point_value(ctx, name, i, A, k), v), -1)
    return out


def atom_shift(ctx, atom):
    """Shift L -> L+1 of an atom component, as expression at L."""
    key = ('atom', atom)
    if key in ctx.red:
        return ctx.red[key]
    (a, comp) = atom
    F = ctx.F
    if a[0] == 'P':
        _, Y, x = a
        if comp == 0:
            out = {(Fraction(0), 0, 0, ((a, 1),)): F(1)}
        else:
            out = {(Fraction(0), 0, 0, ((a, 0),)): F(1),
                   (Fraction(0), 0, 0, ((a, 1),)): -F(_eps(Y)) * (2 * ctx.L + 3) / F(x)}
    else:
        _, name, A, B, mk = a
        ms, bs = reduce_shift(ctx, name, mk)
        out = {}
        for mk2, v in ms.items():
            out = eadd(out, {(Fraction(0), 0, 0, ((('J', name, A, B, mk2), 0),)): v})
        for (k, i), v in bs.items():
            out = eadd(out, escale(point_value(ctx, name, i, B, k), v))
            out = eadd(out, escale(point_value(ctx, name, i, A, k), v), -1)
    ctx.red[key] = out
    return out


def mono_shift(ctx, mono):
    tag, g1, g2, facs = mono
    F = ctx.F
    s = F(D1) ** g1 if g1 >= 0 else F(1) / F(D1) ** (-g1)
    s = s * (F(D2) ** g2 if g2 >= 0 else F(1) / F(D2) ** (-g2))
    out = {(tag, g1, g2, ()): s}
    for f in facs:
        out = emul(out, atom_shift(ctx, f))
    return {k: v for k, v in out.items() if int(v) != 0}


# ------------------------------------------------------------------------------------------------------- the output

def P(Y, x, c):
    return (Fraction(0), 0, 0, ((('P', Y, x), c),))


def multipliers(ctx, j, inner):
    """(M0, M1) expressions for electron j (centre D_j), inner (t < D_j) or outer form."""
    Dj = D1 if j == 1 else D2
    y = ZH * Dj
    F = ctx.F
    if inner:
        return ({P('k', y, 0): 2 * ctx.L + 1, P('k', y, 1): -F(y)}, {P('k', y, 0): F(1)})
    return ({P('i', y, 0): 2 * ctx.L + 1, P('i', y, 1): F(y)}, {P('i', y, 0): -F(1)})


def tvar(name, k, i):
    """placeholder for the t-dependent factor r^k e_i of module `name` (resolved by integration)."""
    return (Fraction(0), 0, 0, ((('T', name, k, i), 0),))


def integrate(ctx, expr, A, B):
    """Replace every t-dependent placeholder ('T', name, k, i) by int_A^B r^k e_i."""
    out = {}
    for mono, v in expr.items():
        tag, g1, g2, facs = mono
        tf = [f for f in facs if f[0][0] == 'T']
        rest = (tag, g1, g2, tuple(f for f in facs if f[0][0] != 'T'))
        assert len(tf) == 1, mono
        _, name, k, i = tf[0][0]
        out = eadd(out, emul({rest: v}, integral_expr(ctx, name, A, B, (k, i))))
    return out


def phi_total(ctx, j):
    """Phi_j(oo) = int_0^oo g_j s^L."""
    Dj = D1 if j == 1 else D2
    out = {}
    for inner, name, A, B in [(True, 'Ui', 0, Dj), (False, 'Uk', Dj, INF)]:
        M0, M1 = multipliers(ctx, j, inner)
        e = eadd(escale(emul(M0, integral_expr(ctx, name, A, B, (0, 0))), -2),
                 escale(emul(M1, integral_expr(ctx, name, A, B, (0, 1))), -2 * ctx.F(ZH)))
        out = eadd(out, e)
    return out


def nested_term(ctx, ja, jp, A, B):
    """int_A^B a_ja(t) Phi_jp(t) dt for a piece [A, B] of t."""
    F = ctx.F
    Dja = D1 if ja == 1 else D2
    Djp = D1 if jp == 1 else D2
    mid_below_a = (B is not INF and B <= Dja)
    mid_below_p = (B is not INF and B <= Djp)
    VY = 'i' if mid_below_a else 'k'
    Ma0, Ma1 = multipliers(ctx, ja, mid_below_a)
    # a(t) = -2 (Ma0 v0 + zH Ma1 t^2 v1): list of (multiplier expr, r power, c)
    a_terms = [(escale(Ma0, -2), 0, 0), (escale(Ma1, -2 * F(ZH)), 2, 1)]
    out = {}
    if mid_below_p:
        UY, side = 'i', '0'
        Mp0, Mp1 = multipliers(ctx, jp, True)
        sgn_b = 1
    else:
        UY, side = 'k', 'oo'
        Mp0, Mp1 = multipliers(ctx, jp, False)
        sgn_b = -1
    # Phi(t) inner: -2 (Mp0 Psi_L + zH Mp1 psi);  outer: Phi(oo) + 2 (Mp0 Psi_L + zH Mp1 psi)
    # Psi_L = int u_L over [0,t] (or [t,oo)) = rho' psi + sum gamma' [s^k u_d]_{0}^{t} (or _{t}^{oo})
    Uname = f'U{UY}'
    ms, bs = reduce(ctx, Uname, (0, 0))
    assert set(ms) <= {(0, 1)}, ms
    rho = ms.get((0, 1), F(0))
    pre = -2 if side == '0' else 2
    Nname = f'N{VY}{UY}{side}'
    # psi-coefficient: pre*(Mp0 rho + zH Mp1);  t^k u_d coefficient: pre*Mp0*gamma*sgn_b
    psi_coef = escale(eadd(escale(Mp0, rho), escale(Mp1, F(ZH))), F(pre))
    ud_terms = [(escale(Mp0, F(pre) * v * F(sgn_b)), k, d) for (k, d), v in bs.items()]
    expr = {}
    for (Ma, pa, c) in a_terms:
        expr = eadd(expr, emul(emul(Ma, psi_coef), {tvar(Nname, pa, c): F(1)}))
        for (Mu, k, d) in ud_terms:
            idx = 2 + 2 * c + d          # p_cd
            expr = eadd(expr, emul(emul(Ma, Mu), {tvar(Nname, pa + k, idx): F(1)}))
    out = integrate(ctx, expr, A, B)
    if not mid_below_p:                  # Phi(oo) * int_A^B a
        aint = {}
        for (Ma, pa, c) in a_terms:
            aint = eadd(aint, emul(Ma, {tvar(f'V{VY}', pa, c): F(1)}))
        out = eadd(out, emul(phi_total(ctx, jp), integrate(ctx, aint, A, B)))
    return out


def output(ctx):
    """F_L = I_C + I_T as an expression in coordinates (monomials)."""
    pieces = [(0, D1), (D1, D2), (D2, INF)]
    out = {}
    for ja, jp in [(1, 2), (2, 1)]:
        for A, B in pieces:
            out = eadd(out, nested_term(ctx, ja, jp, A, B))
    return {k: v for k, v in out.items() if int(v) != 0}


def closure(ctx, start):
    """All monomials reachable from `start` under the shift; returns (list of monomials, shifts)."""
    seen = list(start)
    idx = {m: t for t, m in enumerate(seen)}
    shifts = {}
    t = 0
    while t < len(seen):
        m = seen[t]
        s = mono_shift(ctx, m)
        shifts[m] = s
        for m2 in s:
            if m2 not in idx:
                idx[m2] = len(seen)
                seen.append(m2)
        t += 1
    return seen, shifts
