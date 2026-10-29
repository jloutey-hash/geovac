"""Record BU (docs/predictions/2026-10-28_bottom_up_angle_operator.md): integration-by-parts reduction of the radial
families of case N1, derived symbolically (sympy, L symbolic). No coefficient data is used anywhere in this file.

Three rank-2 families u_L(r) (e^{beta r} times a modified spherical Bessel function of z r), beta = -zX, z = zH:
    I : u_L = e^{beta r} r^{L+2} i_L(z r)     (piece r in [0, D1])
    K : u_L = e^{beta r} r^{L+2} k_L(z r)     (piece r in [D1, D2])
    H : u_L = e^{beta r} r^{1-L} k_L(z r)     (piece r in [D2, oo))
Rules used (checked numerically in `check()`): i_{L+1} = i_{L-1} - (2L+1)/t i_L, k_{L+1} = k_{L-1} + (2L+1)/t k_L,
i_L' = i_{L+1} + L/t i_L, k_L' = -k_{L+1} + L/t k_L, i_{L+1}' = i_L - (L+2)/t i_{L+1}, k_{L+1}' = -k_L - (L+2)/t k_{L+1}.
For each family d/dr (u0, u1) = M(r) (u0, u1) with u0 = u_L, u1 = u_{L+1}, and u_{L+2} = s0 u0 + s1 u1.

IBP identities  int d/dr(r^k u_j) dr = [r^k u_j]  hold for the formal integrals of debug/na_exact.py (its monomial rule
b J_n + n J_{n-1} = [r^n e^{br}] is the fundamental theorem on e^{br} r^n, in every exponential part separately), so
they reduce J(k, j) = int r^k u_j to two masters plus boundary symbols B(k, j) = [r^k u_j(L)] (upper minus lower end).
Masters: I, K: J(0,0) = int u_L, J(0,1) = int u_{L+1};  H: J(0,1), J(-2,1)  (for H, J(0,0) is NOT independent of J(0,1)
modulo boundaries: the reduction gives int u_L = (200 L / 69) int u_{L+1} + boundary terms).
"""
import sympy as sp

L, r = sp.symbols('L r')
beta = sp.Rational(-13, 10)
z = sp.Integer(1)

# family: (M as {(j, l): [(coef, power of r)]}, shift s0, s1 as [(coef, power)])
FAMILIES = {
    'I': ({(0, 0): [(beta, 0), (2 * L + 2, -1)], (0, 1): [(z, -1)],
           (1, 1): [(beta, 0), (1, -1)], (1, 0): [(z, 1)]},
          [(1, 2)], [(-(2 * L + 3) / z, 0)]),
    'K': ({(0, 0): [(beta, 0), (2 * L + 2, -1)], (0, 1): [(-z, -1)],
           (1, 1): [(beta, 0), (1, -1)], (1, 0): [(-z, 1)]},
          [(1, 2)], [((2 * L + 3) / z, 0)]),
    'H': ({(0, 0): [(beta, 0), (1, -1)], (0, 1): [(-z, 1)],
           (1, 1): [(beta, 0), (-(2 * L + 2), -1)], (1, 0): [(-z, -1)]},
          [(1, -2)], [((2 * L + 3) / z, -2)]),
}
MASTERS = {'I': [(0, 0), (0, 1)], 'K': [(0, 0), (0, 1)], 'H': [(0, 1), (-2, 1)]}
# integrals the piece integrands need: I, K: int u_L, int u_{L+1};  H: int u_L, int r^2 u_{L+1}
INTEGRAND = {'I': [(0, 0), (0, 1)], 'K': [(0, 0), (0, 1)], 'H': [(0, 0), (2, 1)]}


def Jsym(k, j):
    return sp.Symbol(f'J_{k}_{j}'.replace('-', 'm'))


def Bsym(k, j):
    return sp.Symbol(f'B_{k}_{j}'.replace('-', 'm'))


def parse(name):
    kind, k, j = name.split('_')
    return kind, int(k.replace('m', '-')), int(j)


def ibp_equations(name, kmin, kmax):
    M, _, _ = FAMILIES[name]
    eqs = []
    for k in range(kmin, kmax + 1):
        for j in (0, 1):
            lhs = k * Jsym(k - 1, j) if k != 0 else sp.Integer(0)
            for l in (0, 1):
                for c, e in M.get((j, l), []):
                    lhs += c * Jsym(k + e, l)
            eqs.append(sp.expand(lhs - Bsym(k, j)))
    return eqs


def reduce_targets(name, targets):
    """Each target (k, j) as a combination of the family's masters and boundary symbols (window search)."""
    masters = {Jsym(*m) for m in MASTERS[name]}
    ks = [t[0] for t in targets] + [m[0] for m in MASTERS[name]]
    for width in range(0, 6):
        for lo in range(0, width + 1):
            kmin, kmax = min(ks) - lo, max(ks) + 1 + (width - lo)
            eqs = ibp_equations(name, kmin, kmax)
            unk = sorted({s for e in eqs for s in e.free_symbols if s.name.startswith('J_')} - masters,
                         key=lambda s: s.name)
            sols = list(sp.linsolve(eqs, unk))
            if len(sols) != 1:
                continue
            sol = dict(zip(unk, sols[0]))
            out, ok = {}, True
            for t in targets:
                if Jsym(*t) in masters:
                    out[t] = Jsym(*t)
                    continue
                if Jsym(*t) not in sol:
                    ok = False
                    break
                expr = sp.expand(sol[Jsym(*t)])
                if {s for s in expr.free_symbols if s.name.startswith('J_')} - masters:
                    ok = False
                    break
                out[t] = expr
            if ok:
                return out
    raise RuntimeError(f'no reduction for {name} {targets}')


def relations(name):
    """'next': master_i(L+1) for each master, and 'integrand': each INTEGRAND target, in masters(L) and B(L)."""
    _, s0, s1 = FAMILIES[name]
    nxt_lin = []                 # master (k, j) at L+1 as a linear combination of J(.,.) at L
    for k, j in MASTERS[name]:
        if j == 0:
            nxt_lin.append([(sp.Integer(1), (k, 1))])
        else:
            nxt_lin.append([(c, (k + e, 0)) for c, e in s0] + [(c, (k + e, 1)) for c, e in s1])
    targets = sorted({t for lin in nxt_lin for _, t in lin} | set(INTEGRAND[name]))
    red = reduce_targets(name, targets)
    nxt = [sp.expand(sum(c * red[t] for c, t in lin)) for lin in nxt_lin]
    integ = [red[t] for t in INTEGRAND[name]]
    return dict(next=nxt, integrand=integ)


def coeffs(expr):
    """{symbol name: rational function of L} for an expression linear in J_ and B_ symbols."""
    expr = sp.expand(expr)
    syms = sorted({s for s in expr.free_symbols if s.name[:2] in ('J_', 'B_')}, key=lambda s: s.name)
    out = {}
    for s in syms:
        c = sp.factor(sp.together(expr.coeff(s)))
        if c != 0:
            out[s.name] = c
    rest = sp.simplify(expr - sum(sp.Symbol(n) * c for n, c in out.items()))
    assert rest == 0, rest
    return out


def check(dps=40):
    """Numerical (mpmath) check of the derivative systems, the shift rules and every reduced relation, as identities
    between definite integrals on a test interval [0.6, 2.9]."""
    import mpmath as mp
    mp.mp.dps = dps
    bt, zz = mp.mpf(-13) / 10, mp.mpf(1)

    def iL(n, t):
        return mp.sqrt(mp.pi / (2 * t)) * mp.besseli(n + mp.mpf(1) / 2, t)

    def kL(n, t):
        return mp.sqrt(2 / (mp.pi * t)) * mp.besselk(n + mp.mpf(1) / 2, t)

    fam = {'I': lambda n, x: mp.e ** (bt * x) * x ** (n + 2) * iL(n, zz * x),
           'K': lambda n, x: mp.e ** (bt * x) * x ** (n + 2) * kL(n, zz * x),
           'H': lambda n, x: mp.e ** (bt * x) * x ** (1 - n) * kL(n, zz * x)}

    def num(c, Lv):
        q = sp.Rational(sp.sympify(c).subs(L, Lv))
        return mp.mpf(q.p) / q.q

    worst = 0
    for name, u in fam.items():
        M, s0, s1 = FAMILIES[name]
        x0 = mp.mpf('1.7')
        for Lv in (0, 3, 7):
            for j in (0, 1):
                d = mp.diff(lambda x: u(Lv + j, x), x0)
                rhs = sum(num(c, Lv) * x0 ** e * u(Lv + l, x0) for l in (0, 1) for c, e in M.get((j, l), []))
                worst = max(worst, abs(d - rhs) / abs(d))
            sh = sum(num(c, Lv) * x0 ** e * u(Lv, x0) for c, e in s0) + \
                sum(num(c, Lv) * x0 ** e * u(Lv + 1, x0) for c, e in s1)
            worst = max(worst, abs(sh - u(Lv + 2, x0)) / abs(u(Lv + 2, x0)))
        rel = relations(name)
        A, B = mp.mpf('0.6'), mp.mpf('2.9')
        for Lv in (0, 2, 5):
            def val(nm, Lq=Lv):
                kind, k, j = parse(nm)
                if kind == 'J':
                    return mp.quad(lambda x: x ** k * u(Lq + j, x), [A, B])
                return B ** k * u(Lq + j, B) - A ** k * u(Lq + j, A)
            checks = []
            for (k, j), e in zip(MASTERS[name], rel['next']):
                checks.append((e, mp.quad(lambda x: x ** k * u(Lv + 1 + j, x), [A, B])))
            for (k, j), e in zip(INTEGRAND[name], rel['integrand']):
                checks.append((e, mp.quad(lambda x: x ** k * u(Lv + j, x), [A, B])))
            for e, target in checks:
                v = sum(num(c, Lv) * val(n) for n, c in coeffs(e).items())
                worst = max(worst, abs(v - target) / abs(target))
    return worst


if __name__ == '__main__':
    for name in FAMILIES:
        rel = relations(name)
        for key in ('next', 'integrand'):
            for e in rel[key]:
                print(name, key, coeffs(e))
    print('numerical check, worst relative error:', check())


# ---------------------------------------------------------------------------------------------------------------
# Rank-one reduction.  Allowing IBP identities with negative powers k (for the piece r in [0, D1] the r = 0 end of
# [r^k u_j] still vanishes for the k that occur: u_j ~ r^{2L+2+j}), every family has the relation
#     int u_L = rho(L) int u_{L+1} + boundary terms,
# so the single master phi = J(0,1) = int u_{L+1} obeys a FIRST-order inhomogeneous recurrence.  (The second
# cohomology class, J(-2,1), is never needed.)  MASTERS1 = the master set used for the reduction; only J(0,1) may
# survive in the relations that are used.

MASTERS1 = {'I': [(0, 1), (-2, 1)], 'K': [(0, 1), (-2, 1)], 'H': [(0, 1), (-2, 1)]}


def relations1(name):
    """'next': phi(L+1) = int u_{L+2} and 'integrand': the INTEGRAND targets, in J(0,1) and B(L) only."""
    global MASTERS
    saved = MASTERS
    MASTERS = dict(saved)
    MASTERS[name] = MASTERS1[name]
    try:
        _, s0, s1 = FAMILIES[name]
        lin = [(c, (e, 0)) for c, e in s0] + [(c, (e, 1)) for c, e in s1]
        targets = sorted({t for _, t in lin} | set(INTEGRAND[name]))
        red = reduce_targets(name, targets)
    finally:
        MASTERS = saved
    def canon(e):
        return sum(sp.Symbol(n) * c for n, c in coeffs(e).items())
    nxt = canon(sum(c * red[t] for c, t in lin))
    integ = [canon(red[t]) for t in INTEGRAND[name]]
    for e in [nxt] + integ:
        assert Jsym(-2, 1) not in e.free_symbols, (name, e)
    return dict(next=[nxt], integrand=integ)


def check1(dps=40):
    """Numerical check of relations1 on [0.6, 2.9] and, for family I, on [0, 1.8] (the r = 0 end)."""
    import mpmath as mp
    mp.mp.dps = dps
    bt = mp.mpf(-13) / 10

    def iL(n, t):
        return mp.sqrt(mp.pi / (2 * t)) * mp.besseli(n + mp.mpf(1) / 2, t)

    def kL(n, t):
        return mp.sqrt(2 / (mp.pi * t)) * mp.besselk(n + mp.mpf(1) / 2, t)

    fam = {'I': lambda n, x: mp.e ** (bt * x) * x ** (n + 2) * iL(n, x),
           'K': lambda n, x: mp.e ** (bt * x) * x ** (n + 2) * kL(n, x),
           'H': lambda n, x: mp.e ** (bt * x) * x ** (1 - n) * kL(n, x)}

    def num(c, Lv):
        q = sp.Rational(sp.sympify(c).subs(L, Lv))
        return mp.mpf(q.p) / q.q
    worst = 0
    for name, u in fam.items():
        rel = relations1(name)
        ivs = [(mp.mpf('0.6'), mp.mpf('2.9'))] + ([(mp.mpf(0), mp.mpf('1.8'))] if name == 'I' else []) + \
            ([(mp.mpf('2.6'), mp.inf)] if name == 'H' else [])
        for A, B in ivs:
            for Lv in (0, 1, 2, 5):
                def val(nm, Lq=Lv):
                    kind, k, j = parse(nm)
                    if kind == 'J':
                        return mp.quad(lambda x: x ** k * u(Lq + j, x), [A, B])
                    up = 0 if B == mp.inf else B ** k * u(Lq + j, B)
                    lo = 0 if A == 0 else A ** k * u(Lq + j, A)
                    return up - lo
                checks = [(rel['next'][0], mp.quad(lambda x: u(Lv + 2, x), [A, B]))]
                for (k, j), e in zip(INTEGRAND[name], rel['integrand']):
                    checks.append((e, mp.quad(lambda x: x ** k * u(Lv + j, x), [A, B])))
                for e, target in checks:
                    v = sum(num(c, Lv) * val(n) for n, c in coeffs(e).items())
                    worst = max(worst, abs(v - target) / abs(target))
    return worst
