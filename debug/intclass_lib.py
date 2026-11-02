"""Record IC (docs/predictions/2026-11-01_integrated_class.md): is a linear differential operator TRIANGULAR over
Qbar(w), i.e. are all its composition factors of order one?  Everything mod one prime p (p < 2^31, so numpy int64
products are exact before reduction).

Operator = list c[0..R] of flint.nmod_poly, L = sum_k c_k(w) (d/dw)^k.

Instrument ("order-one peel"): find a hyperexponential solution y = e^{lam w} * prod f_i^{e_i} * N(w) (lam from the
slope-one edge at infinity, e_i from the indicial roots at the finite singular factors f_i of the leading
coefficient, N a polynomial with the degree bound from the indicial polynomial at infinity), right-divide
L = Q o (d - y'/y) exactly (the remainder is checked to be zero), repeat on Q.  If the operator reaches order 0, every
composition factor has order one (TRIANGULAR).  If an operator of order >= 2 has no hyperexponential solution, its
minimal submodule has dimension >= 2, i.e. an irreducible composition factor of order >= 2 exists (STUCK).  For the
stuck residual we report its local data, whether its adjoint has an order-one right factor (= an order-one left
factor of the residual), repeated indicial roots (a repeated exponent forces a logarithm), and the exponential
types at infinity.
"""
from __future__ import annotations

import time
from fractions import Fraction
from math import comb

import flint
import numpy as np

SMALL_INT = 10 ** 6


# ---------------------------------------------------------------------------------------------- polynomial helpers
def zp(p):
    return flint.nmod_poly([], p)


def cp(v, p):
    return flint.nmod_poly([int(v) % p], p)


def signed(r, p):
    r = int(r) % p
    return r - p if r > p // 2 else r


def ratrec(a, p, bound=1000):
    """small rational with a = n/d mod p, |n|, d <= bound, else None"""
    a = int(a) % p
    r0, r1, s0, s1 = p, a, 0, 1
    while r1 > bound:
        q = r0 // r1
        r0, r1, s0, s1 = r1, r0 - q * r1, s1, s0 - q * s1
    if s1 == 0 or abs(s1) > bound:
        return None
    fr = Fraction(r1, s1)
    if (fr.numerator - a * fr.denominator) % p:
        return None
    return fr


def to_fp(fr, p):
    fr = Fraction(fr)
    return fr.numerator % p * pow(fr.denominator % p, p - 2, p) % p


def content_reduce(c):
    g = None
    for ck in c:
        if ck.degree() >= 0:
            g = ck if g is None else g.gcd(ck)
            if g.degree() == 0:
                break
    if g is None or g.degree() <= 0:
        lc = c[-1].leading_coefficient()
        return [ck * (1 / lc) for ck in c]
    out = [divmod(ck, g)[0] for ck in c]
    lc = out[-1].leading_coefficient()
    return [ck * (1 / lc) for ck in out]


def ff_poly(k, p):
    """falling factorial n(n-1)...(n-k+1) as nmod_poly in n"""
    f = flint.nmod_poly([1], p)
    for m in range(k):
        f *= flint.nmod_poly([(-m) % p, 1], p)
    return f


def ff_val(n, k, p):
    v = 1
    for m in range(k):
        v = v * ((n - m) % p) % p
    return v


# ---------------------------------------------------------------------------------------------- operator algebra
def from_theta(cth, p, s=1):
    """sum_k c_k(x) theta^k with theta = s x d/dx  ->  D-form.  theta^k = s^k sum_j S2(k,j) x^j D^j."""
    R = len(cth) - 1
    S = [[0] * (R + 1) for _ in range(R + 1)]
    S[0][0] = 1
    for k in range(1, R + 1):
        for j in range(1, k + 1):
            S[k][j] = j * S[k - 1][j] + S[k - 1][j - 1]
    out = [zp(p) for _ in range(R + 1)]
    for k, ck in enumerate(cth):
        for j in range(k + 1):
            if S[k][j]:
                xj = flint.nmod_poly([0] * j + [pow(s, k, p) * S[k][j] % p], p)
                out[j] += ck * xj
    return content_reduce(out)


def pullback_quadratic(c, a0, a2, p):
    """x = a0 + a2 w^2 (a0, a2 in F_p).  d/dx = s^{-1} d/dw, s = dx/dw = 2 a2 w.
    d_x^k = s^{-(2k-1)} sum_m E_{k,m} d_w^m;  E_{k+1} = -(2k-1) s' E_k + s (E_k' + shift E_k)."""
    R = len(c) - 1
    xw = flint.nmod_poly([a0 % p, 0, a2 % p], p)
    s = flint.nmod_poly([0, 2 * a2 % p], p)
    sd = s.derivative()
    E = {1: {1: flint.nmod_poly([1], p)}}
    for k in range(1, R):
        new = {}
        for m, e in E[k].items():
            new[m] = new.get(m, zp(p)) + cp(-(2 * k - 1), p) * sd * e + s * e.derivative()
            new[m + 1] = new.get(m + 1, zp(p)) + s * e
        E[k + 1] = new
    out = [zp(p) for _ in range(R + 1)]
    for k, ck in enumerate(c):
        ckw = ck(xw) if ck.degree() >= 0 else zp(p)
        if k == 0:
            out[0] += ckw * s ** (2 * R - 1)
            continue
        for m, e in E[k].items():
            out[m] += ckw * e * s ** (2 * R - 1 - (2 * k - 1))
    return content_reduce(out)


def shift_lambda(c, lam, p):
    """e^{-lam w} L e^{lam w}:  d -> d + lam"""
    R = len(c) - 1
    out = [zp(p) for _ in range(R + 1)]
    for k, ck in enumerate(c):
        lp = 1
        for j in range(k, -1, -1):
            out[j] += ck * (comb(k, j) * lp % p)
            lp = lp * lam % p
    return out


def gauge(c, f, e, p):
    """operator for z with y = f^e z (e in F_p):  L(f^e z) = f^{e-R} sum_j z^(j) b_j,
    b_j = sum_{k>=j} c_k binom(k,j) P_{k-j} f^{R-k+j},  (f^e)^(i) = f^{e-i} P_i, P_{i+1} = P_i' f + (e-i) f' P_i."""
    R = len(c) - 1
    fd = f.derivative()
    Pl = [flint.nmod_poly([1], p)]
    for i in range(R):
        Pl.append(Pl[-1].derivative() * f + cp(e - i, p) * fd * Pl[-1])
    fpow = [flint.nmod_poly([1], p)]
    for i in range(R):
        fpow.append(fpow[-1] * f)
    out = [zp(p) for _ in range(R + 1)]
    for k, ck in enumerate(c):
        if ck.degree() < 0:
            continue
        for j in range(k + 1):
            out[j] += ck * (comb(k, j) % p) * Pl[k - j] * fpow[R - k + j]
    return content_reduce(out)


def adjoint(c, p):
    """L* = sum_k (-1)^k d^k c_k = sum_j [sum_{k>=j} (-1)^k binom(k,j) c_k^(k-j)] d^j"""
    R = len(c) - 1
    out = [zp(p) for _ in range(R + 1)]
    for k, ck in enumerate(c):
        d = ck
        for i in range(k + 1):          # i = k - j
            j = k - i
            out[j] += d * ((-1) ** k * comb(k, j) % p)
            d = d.derivative()
    return content_reduce(out)


def apply_op(c, y, p):
    """L(y) for a polynomial y"""
    out = zp(p)
    d = y
    for ck in c:
        out += ck * d
        d = d.derivative()
    return out


def right_divide_order1(c, A, B, p):
    """L = Q o (d - A/B).  Returns (Q as polynomial operator, remainder_is_zero)."""
    R = len(c) - 1
    Bd = B.derivative()
    U = [A]
    for i in range(R):
        U.append(U[-1].derivative() * B - cp(i + 1, p) * U[-1] * Bd)
    Bp = [flint.nmod_poly([1], p)]
    for i in range(R + 1):
        Bp.append(Bp[-1] * B)
    Q = [zp(p) for _ in range(R)]       # q_k = Q[k] / B^{R-1-k}
    Q[R - 1] = c[R]
    for m in range(R - 1, 0, -1):       # Q_{m-1} = a_m B^{R-m} + sum_{k>=m} binom(k,k-m) Q_k U_{k-m}
        acc = c[m] * Bp[R - m]
        for k in range(m, R):
            acc += Q[k] * U[k - m] * (comb(k, k - m) % p)
        Q[m - 1] = acc
    rem = c[0] * Bp[R]
    for k in range(R):
        rem += Q[k] * U[k]
    newc = [Q[k] * Bp[k] for k in range(R)]
    return content_reduce(newc), rem.degree() < 0


# ---------------------------------------------------------------------------------------------- local data
def valuation(poly, f):
    m, q = 0, poly
    if q.degree() < 0:
        return 10 ** 9, q
    while True:
        d, r = divmod(q, f)
        if r.degree() >= 0:
            return m, q
        q = d
        m += 1


def indicial_at_factor(c, f, p):
    """Indicial data at the roots of an F_p-irreducible factor f of the leading coefficient.
    Returns dict(regular, mu, deg_chi, roots=[(r_int_mod_p, mult)])  (roots in F_p only; deg_chi - sum mult are the
    exponents not in F_p, necessarily irrational)."""
    R = len(c) - 1
    fd = f.derivative()
    vals = [valuation(ck, f) for ck in c]
    ords = [v[0] for v in vals]
    mu = min(ords[k] - k for k in range(R + 1))
    regular = (ords[R] - R) == mu
    lead = {}
    for k in range(R + 1):
        if ords[k] - k == mu:
            g = vals[k][1] % f
            lead[k] = (g * (fd ** ords[k] % f)) % f
    kmax = max(lead)
    # chi(n) = sum_k lead_k(t) ff_k(n);  h_j(n) = coefficient of t^j
    dfa = max(f.degree(), 1)
    H = [[0] * (kmax + 1) for _ in range(dfa)]
    for k, lk in lead.items():
        ffk = ff_poly(k, p).coeffs()
        lc_ = lk.coeffs()
        for j, tv in enumerate(lc_):
            if int(tv) == 0:
                continue
            for i, nv in enumerate(ffk):
                H[j][i] = (H[j][i] + int(tv) * int(nv)) % p
    G = None
    for row in H:
        h = flint.nmod_poly(row, p)
        if h.degree() >= 0:
            G = h if G is None else G.gcd(h)
    roots = [(int(r), int(m)) for r, m in G.roots()] if G is not None and G.degree() > 0 else []
    return dict(regular=regular, mu=mu, deg_chi=kmax, roots=roots)


def classify_root(r, p):
    s = signed(r, p)
    if abs(s) <= SMALL_INT:
        return ('int', s)
    q = ratrec(r, p)
    if q is not None:
        return ('rat', q)
    return ('other', r)


def finite_points(c, p):
    """every F_p-irreducible factor of lc: indicial data, classified roots"""
    out = []
    for f, mult in c[-1].factor()[1]:
        d = indicial_at_factor(c, f, p)
        d['factor'] = f
        d['fdeg'] = f.degree()
        d['croots'] = [(classify_root(r, p), m) for r, m in d['roots']]
        ints = sorted(v for (k, v), m in d['croots'] if k == 'int' for _ in range(m))
        d['ints'] = ints
        d['apparent_type'] = d['regular'] and d['deg_chi'] == len(c) - 1 and len(ints) == len(c) - 1 \
            and len(set(ints)) == len(ints) and min(ints) >= 0
        d['repeated'] = [cr for cr, m in d['croots'] if m >= 2]
        d['unaccounted'] = d['deg_chi'] - sum(m for _, m in d['roots'])
        out.append(d)
    return out


def infinity_data(c, p):
    """degrees, rank flag, slope-one edge polynomial and its F_p roots"""
    R = len(c) - 1
    deg = [ck.degree() for ck in c]
    dmax = max(deg)
    rank_le_1 = deg[R] == dmax
    E = flint.nmod_poly([int(c[k].leading_coefficient()) if deg[k] == dmax else 0 for k in range(R + 1)], p)
    roots = [(int(r), int(m)) for r, m in E.roots()] if E.degree() > 0 else []
    unsplit = E.degree() - sum(m for _, m in roots)
    return dict(deg=deg, rank_le_1=rank_le_1, edge_degree=E.degree(), lam_roots=roots, lam_unsplit=unsplit)


# ---------------------------------------------------------------------------------------------- solutions
def poly_solutions(b, p, Nmax=6000):
    """polynomial solutions of sum_j b_j(w) d^j; returns (list of nmod_poly basis, info)"""
    R = len(b) - 1
    deg = [bj.degree() for bj in b]
    m = max(deg[j] - j for j in range(R + 1) if deg[j] >= 0)
    chi = zp(p)
    for j in range(R + 1):
        if deg[j] >= 0 and deg[j] - j == m:
            chi += ff_poly(j, p) * int(b[j].leading_coefficient())
    cand = [signed(r, p) for r, _ in chi.roots()] if chi.degree() > 0 else []
    cand = [n for n in cand if 0 <= n <= SMALL_INT]
    if not cand:
        return [], dict(N0=None)
    N0 = max(cand)
    if N0 > Nmax:
        return None, dict(N0=N0, too_big=True)
    nrow = N0 + m + 1
    ncol = N0 + 1
    M = np.zeros((max(nrow, 1), ncol), dtype=np.int64)
    bv = [np.array([int(v) for v in bj.coeffs()], dtype=np.int64) if deg[j] >= 0 else None for j, bj in enumerate(b)]
    for j in range(R + 1):
        if bv[j] is None:
            continue
        L = len(bv[j])
        for n in range(j, ncol):
            fv = ff_val(n, j, p)
            if fv == 0:
                continue
            r0 = n - j
            M[r0:r0 + L, n] = (M[r0:r0 + L, n] + (bv[j] * fv) % p) % p
    Mf = flint.nmod_mat(M.tolist(), p)
    X, nul = Mf.nullspace()
    sols = []
    for col in range(nul):
        sols.append(flint.nmod_poly([int(X[i, col]) for i in range(ncol)], p))
    return sols, dict(N0=N0, nrow=nrow, ncol=ncol)


def find_hyperexp(c, p, Nmax=6000, max_combos=64, verbose=False):
    """one hyperexponential solution of L over Qbar(w) (searched mod p), or None.
    Returns (dict(lam, gauges=[(f, e)], N), diag)."""
    R = len(c) - 1
    fps = finite_points(c, p)
    inf = infinity_data(c, p)
    diag = dict(order=R, inf=inf, points=[dict(fdeg=d['fdeg'], regular=d['regular'], deg_chi=d['deg_chi'],
                                               croots=[(k, str(v), m) for (k, v), m in d['croots']],
                                               unaccounted=d['unaccounted'], apparent_type=d['apparent_type'])
                                          for d in fps], flags=[])
    if not inf['rank_le_1']:
        diag['flags'].append('rank>1 at infinity')
    if inf['lam_unsplit']:
        diag['flags'].append(f"lambda edge has {inf['lam_unsplit']} roots outside F_p")
    # exponent classes mod Z per factor: rational non-integer classes
    per = []
    for d in fps:
        classes = {}
        for (k, v), m in d['croots']:
            if k == 'int':
                classes.setdefault(Fraction(0), []).append(v)
            elif k == 'rat':
                cl = v - (v.numerator // v.denominator)
                classes.setdefault(cl, []).append(v)
        if any(k == 'other' for (k, v), m in d['croots']):
            diag['flags'].append(f"exponent not small rational at a degree-{d['fdeg']} factor")
        per.append((d, classes))
    nonint = [(i, sorted(cl for cl in classes if cl != 0)) for i, (d, classes) in enumerate(per)]
    nonint = [(i, cls) for i, cls in nonint if cls]
    combos = [[]]
    for i, cls in nonint:
        combos = [cb + [(i, cl)] for cb in combos for cl in [Fraction(0)] + cls]
        if len(combos) > max_combos:
            diag['flags'].append('too many exponent combinations; truncated')
            combos = combos[:max_combos]
            break
    lams = [r for r, _ in inf['lam_roots']]
    if 0 not in lams:
        lams = [0] + lams
    for cb in combos:
        choice = dict(cb)
        gauges = []
        for i, (d, classes) in enumerate(per):
            cl = choice.get(i, Fraction(0))
            vals = classes.get(cl, [])
            if not vals:
                continue
            # exponents in this class: e = cl + integer; pole bound -> smallest
            mn = min(vals)
            e = mn if cl != 0 else min(0, mn)
            if e == 0:
                continue
            gauges.append((d['factor'], e))
        cg = c
        for f, e in gauges:
            cg = gauge(cg, f, to_fp(e, p), p)
        for lam in lams:
            cl_ = shift_lambda(cg, lam, p) if lam else cg
            sols, info = poly_solutions(cl_, p, Nmax)
            if sols is None:
                diag['flags'].append(f"degree bound {info['N0']} > Nmax")
                continue
            if sols:
                return dict(lam=lam, gauges=gauges, N=sols[0], N0=info['N0'], nsol=len(sols)), diag
    return None, diag


def log_derivative(sol, p):
    """u = y'/y for y = e^{lam w} prod f^e N, as (A, B) with gcd 1"""
    N = sol['N']
    B = N
    A = N.derivative() + N * sol['lam']
    for f, e in sol['gauges']:
        fd = f.derivative()
        A = A * f + B * fd * to_fp(e, p)
        B = B * f
    g = A.gcd(B)
    if g.degree() > 0:
        A, B = divmod(A, g)[0], divmod(B, g)[0]
    return A, B


def peel(c, p, Nmax=6000, max_steps=None, verbose=True, t_budget=7200):
    """order-one peel.  Returns dict(verdict, steps, residual_diag, adjoint_has_order1, residual)."""
    t0 = time.time()
    steps = []
    cur = content_reduce(c)
    while len(cur) - 1 >= 2:
        if time.time() - t0 > t_budget:
            return dict(verdict='BUDGET', steps=steps, residual=cur)
        sol, diag = find_hyperexp(cur, p, Nmax)
        if sol is None:
            res = dict(verdict='STUCK', steps=steps, residual_order=len(cur) - 1, residual_diag=diag,
                       residual=cur)
            adj = adjoint(cur, p)
            sol2, d2 = find_hyperexp(adj, p, Nmax)
            res['adjoint_has_order1'] = sol2 is not None
            res['adjoint_flags'] = d2['flags']
            res['time'] = time.time() - t0
            return res
        A, B = log_derivative(sol, p)
        new, ok = right_divide_order1(cur, A, B, p)
        if not ok:
            raise RuntimeError('right division left a remainder')
        steps.append(dict(order=len(cur) - 1, lam=signed(sol['lam'], p),
                          gauges=[(g.degree(), str(e)) for g, e in sol['gauges']], Ndeg=sol['N'].degree(),
                          N0=sol['N0'], nsol=sol['nsol'], flags=diag['flags'],
                          newdeg=max(x.degree() for x in new), t=round(time.time() - t0, 1)))
        if verbose:
            print('   peel', steps[-1], flush=True)
        cur = new
        if max_steps is not None and len(steps) >= max_steps:
            return dict(verdict='MAXSTEPS', steps=steps, residual=cur)
    # an order-one remainder is itself an order-one composition factor
    return dict(verdict='TRIANGULAR', steps=steps, time=time.time() - t0, last_order=len(cur) - 1)


# ---------------------------------------------------------------------------------------------- classification
def classify_residual(res, p):
    """Pre-registered reading of a STUCK residual (record IC):
      clean        no flag at the stuck step (lambda edge split over F_p, exponents small rationals, degree bounds
                   within Nmax) -- only a clean STUCK counts for the gate
      lam          the exponential types at infinity (distinct lambda with multiplicity)
      logs         rational points (degree-1 factors, integer exponents, not apparent-type) where a logarithm is
                   forced (log-free formal solutions fewer than exponents)
      type         'fuchsian-elliptic-type'  order 2, one exponential type only (regular-type infinity up to a
                                              twist e^{lam w}), a forced log at some finite point (irreducible
                                              order 2 + unipotent local monodromy => Galois group contains SL2)
                   'confluent-genus1-type'   >= 4 distinct lambda
                   'confluent-genus0-type'   order 2 with exactly 2 distinct lambda (Bessel-like)
                   'fuchsian-higher'         one exponential type, order > 2
                   'other'                   anything else
    """
    d = res['residual_diag']
    cur = res['residual']
    lam = []
    for r, m in d['inf']['lam_roots']:
        q = classify_root(r, p)
        lam.append((str(q[1]), m))
    nd = len(lam)
    order = res['residual_order']
    logs = []
    for pt in finite_points(cur, p):
        if pt['fdeg'] == 1 and not pt['apparent_type'] and pt['croots']                 and all(k == 'int' for (k, v), m in pt['croots']):
            x0 = (-int(pt['factor'].coeffs()[0]) * pow(int(pt['factor'].coeffs()[1]), p - 2, p)) % p
            ints = pt['ints']
            if log_free_count(cur, x0, p, ints) < len(ints):
                logs.append(signed(x0, p))
    single = nd == 1
    if order == 2 and single and logs:
        typ = 'fuchsian-elliptic-type'
    elif nd >= 4:
        typ = 'confluent-genus1-type'
    elif order == 2 and nd == 2:
        typ = 'confluent-genus0-type'
    elif single and order > 2:
        typ = 'fuchsian-higher'
    else:
        typ = 'other'
    genuine = [dict(fdeg=pt['fdeg'], croots=pt['croots'], regular=pt['regular']) for pt in d['points']
               if not pt['apparent_type']]
    return dict(clean=not d['flags'], flags=d['flags'], order=order, lam=lam, logs=logs, type=typ,
                adjoint_has_order1=res.get('adjoint_has_order1'), genuine_points=genuine)


# ---------------------------------------------------------------------------------------------- minimal operator
def trim(c):
    c = list(c)
    while len(c) > 1 and c[-1].degree() < 0:
        c.pop()
    return c


def d_power_compose(c, m, p):
    """d^m o L = sum_k sum_i binom(m,i) c_k^(i) d^(k+m-i)"""
    R = len(c) - 1
    out = [zp(p) for _ in range(R + m + 1)]
    for k, ck in enumerate(c):
        d = ck
        for i in range(m + 1):
            out[k + m - i] += d * (comb(m, i) % p)
            d = d.derivative()
    return out


def right_rem(L1, L2, p):
    """pseudo-remainder of L1 by L2 (left multiples only), order < ord L2"""
    L1 = trim(L1)
    b = len(L2) - 1
    while len(L1) - 1 >= b and not (len(L1) == 1 and L1[0].degree() < 0):
        a = len(L1) - 1
        T = d_power_compose(L2, a - b, p)
        l1, lt = L1[-1], T[-1]
        g = l1.gcd(lt)
        f1, ft = divmod(lt, g)[0], divmod(l1, g)[0]
        L1 = trim([f1 * x - ft * y for x, y in zip(L1, T)])
        if len(L1) == 1 and L1[0].degree() < 0:
            return L1
        L1 = content_reduce(L1)
    return L1


def gcrd(L1, L2, p):
    A, B = trim(L1), trim(L2)
    if len(A) < len(B):
        A, B = B, A
    while not (len(B) == 1 and B[0].degree() < 0):
        Rm = right_rem(A, B, p)
        A, B = B, Rm
        if len(A) == 1:                  # order 0: gcrd trivial
            break
    return content_reduce(A)


# ---------------------------------------------------------------------------------------------- forced logarithms
def shift_op(c, x0, p):
    """L in the local coordinate w - x0: coefficient lists (ascending) of c_k(x0 + t)"""
    out = []
    for ck in c:
        out.append(ck(flint.nmod_poly([x0 % p, 1], p)) if ck.degree() >= 0 else zp(p))
    return out


def log_free_count(c, x0, p, ints, extra=40):
    """number of log-free formal solutions sum_{n>=e_min} c_n t^n at t = w - x0 (regular point of the indicial
    data with integer exponents ints).  A log is forced iff this is < len(ints) (counted with multiplicity)."""
    cs = shift_op(c, x0, p)
    R = len(cs) - 1
    vs = [ck.coeffs() for ck in cs]
    # L(t^n) = sum_k sum_m a_{k,m} n^(k) t^(n-k+m); delta = min_k (val_k - k)
    vals = []
    for v in vs:
        nz = [i for i, a in enumerate(v) if int(a)]
        vals.append(nz[0] if nz else 10 ** 9)
    delta = min(vals[k] - k for k in range(R + 1))

    def chi(j, n):     # coefficient of t^(n + delta + j) in L(t^n)
        s = 0
        for k in range(R + 1):
            m = k + delta + j
            if 0 <= m < len(vs[k]):
                s = (s + int(vs[k][m]) * ff_val(n, k, p)) % p
        return s
    emin, emax = min(ints), max(ints)
    roots = sorted(set(ints))
    dim = len(roots)
    coef = {}                           # n -> vector over the free parameters
    cons = []
    for n in range(emin, emax + 1):
        rhs = [0] * dim
        for j in range(1, n - emin + 1):
            cj = chi(j, n - j)
            if cj and (n - j) in coef:
                rhs = [(r - cj * v) % p for r, v in zip(rhs, coef[n - j])]
        c0 = chi(0, n)
        if c0 == 0:
            if n in roots:
                cons.append(rhs[:])      # sum chi_j c_{n-j} must vanish
                v = [0] * dim
                v[roots.index(n)] = 1
                coef[n] = v
            else:
                coef[n] = [0] * dim
        else:
            inv = pow(c0, p - 2, p)
            coef[n] = [r * inv % p for r in rhs]
    if cons:
        M = flint.nmod_mat(cons, p)
        rk = M.rank()
    else:
        rk = 0
    return dim - rk


def minimal_operator(seq, p, Rmax=24, margin=30, extra_orders=2):
    """Legendre-coefficient sequence -> its minimal bond-angle operator (D-form in x).  The AO3 search
    (na_modguess) returns the smallest order FINDABLE within the rows, which can exceed the minimal order (an
    order-degree trade-off: a left multiple M o L_min).  The minimal operator is the gcrd of independent
    annihilators: every null vector at (R0, Q0) and the minimal-degree ones at R0 + 1 .. R0 + extra_orders."""
    import na_modguess as G
    import na_modanalyze as A
    B = G.Basis([flint.nmod(int(v), p) for v in seq], p)
    found = []
    R0 = None
    for R in range(1, Rmax + 1):
        Q0, nul = G.min_degree_bisect(B, R, margin, len(seq))
        if Q0 is not None:
            R0 = R
            found.append((R, Q0, nul))
            break
    if R0 is None:
        return None
    ops = []
    for R in range(R0, R0 + extra_orders + 1):
        Q0, nul = G.min_degree_bisect(B, R, margin, len(seq))
        if Q0 is None:
            continue
        if R > R0:
            found.append((R, Q0, nul))
        nv, keys, (N, nrow) = G.nullity(B, R, Q0, margin)
        for col in range(nv):
            ops.append(A.standard_form(keys, [N[i, col] for i in range(len(keys))], p))
    g = ops[0]
    for o in ops[1:]:
        g = gcrd(g, o, p)
    return dict(found=found, c=content_reduce(g), order=len(g) - 1, degree=max(x.degree() for x in g))
