"""QFD assembly: closed-form S / h / g tensors -> Loewdin -> FCI, at any precision.

The production path here calls NO quadrature routine.  Every matrix element is a
sympy expression from `qfd_core` (one-electron, one-center ERI, exchange) or from
`geovac.two_center_eri` ((AA|BB), hybrid); the only numerics are (a) evaluating
those expressions at a requested `dps` and (b) the linear algebra (Loewdin
S^{-1/2} + FCI eigenvalue) which is done in mpmath at the same `dps`.

Orbital spec: (center, Z_orbital, n) with center in {"A", "B"}.  Z_orbital is the
exponent parameter (rate a = Z/n), NOT necessarily the nuclear charge; the
nuclear charges are supplied separately.
"""
from __future__ import annotations

from itertools import combinations

import sympy as sp
from mpmath import mp

from geovac import qfd_core as Q


# ------------------------------------------------------------------ one-electron

def build_S_h(orbs, ZA, ZB, R):
    """Exact S and h_core matrices (lists of lists of sympy expressions)."""
    n = len(orbs)
    S = [[None] * n for _ in range(n)]
    h = [[None] * n for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            S[i][j] = S[j][i] = Q.overlap(orbs[i], orbs[j], R)
            h[i][j] = h[j][i] = Q.h_core(orbs[i], orbs[j], ZA, ZB, R)
    return S, h


def h_core_check(orbs, ZA, ZB, R):
    """Second, independent closed-form route to h (hydrogenic eigen-trick).

        h_ij = E_j S_ij + (zeta_j - Z_own) <i|1/r_own|j> - Z_other <i|1/r_other|j>

    with E_j = -zeta_j^2/(2 n_j^2).  Exact for any zeta because the correction
    term restores the difference between the orbital's own screening charge and
    the physical nuclear charge.
    """
    n = len(orbs)
    out = [[None] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            oj = orbs[j]
            zj = sp.nsimplify(oj[1])
            Ej = -zj ** 2 / (2 * oj[2] ** 2)
            own = oj[0]
            other = "B" if own == "A" else "A"
            Zown = ZA if own == "A" else ZB
            Zoth = ZA if other == "A" else ZB
            out[i][j] = sp.expand(
                Ej * Q.overlap(orbs[i], oj, R)
                + (zj - Zown) * Q._inv_r(orbs[i], oj, own, R)
                - Zoth * Q._inv_r(orbs[i], oj, other, R))
    return out


# ------------------------------------------------------------------ two-electron

def classify(quartet):
    """Return ('one-center'|'aabb'|'hybrid'|'exchange', reordered indices)."""
    ca, cb, cc, cd = (o[0] for o in quartet)
    if ca == cb == cc == cd:
        return "one-center", quartet
    rho1_one = (ca == cb)
    rho2_one = (cc == cd)
    if rho1_one and rho2_one:
        return "aabb", quartet
    if rho1_one or rho2_one:
        return "hybrid", quartet
    return "exchange", quartet


def eri_closed_form(quartet, R, tau_max: int = 10, exchange_hp_dps=None):
    """(ab|cd) in chemist notation, as a sympy expression, plus its class tag.

    `exchange_hp_dps` routes the exchange class through `qfd_core.exchange_hp`,
    which evaluates every closed-form factor at that precision and accumulates
    the (non-terminating, heteronuclear) tau series numerically.  The returned
    value is then an mpmath mpf rather than a sympy expression, together with the
    per-tau list so the truncation tail can be audited.
    """
    a, b, c, d = quartet
    kind, _ = classify(quartet)
    if kind == "one-center":
        return Q.one_center_eri(a, b, c, d), "one-center"
    if kind == "aabb":
        assert a[1] == b[1] and c[1] == d[1], "aabb needs one Z per side"
        return (Q.aabb_closed_form(a[1], (a[2], 0, 0), (b[2], 0, 0),
                                   c[1], (c[2], 0, 0), (d[2], 0, 0), R),
                "(AA|BB)")
    if kind == "hybrid":
        if a[0] == b[0]:                     # rho1 is the one-center pair
            P, Qb, U, V = a, b, c, d
        else:                                # (ab|cd) = (cd|ab)
            P, Qb, U, V = c, d, a, b
        if U[0] != P[0]:                     # rho2 real => swap its members
            U, V = V, U
        assert U[0] == P[0] and V[0] != P[0]
        assert P[1] == Qb[1] == U[1], "hybrid needs one Z on the shared center"
        return (Q.hybrid_closed_form(P[1], (P[2], 0, 0), (Qb[2], 0, 0),
                                     (U[2], 0, 0), V[1], (V[2], 0, 0), R),
                "hybrid")
    # exchange: put the A-center member of each density first
    if a[0] != "A":
        a, b = b, a
    if c[0] != "A":
        c, d = d, c
    assert a[0] == c[0] == "A" and b[0] == d[0] == "B"
    assert a[1] == c[1] and b[1] == d[1], "exchange needs one Z per center"
    tm = tau_max((a, b, c, d)) if callable(tau_max) else tau_max
    if exchange_hp_dps is not None:
        tot, per = Q.exchange_hp(a[1], (a[2], 0, 0), (b[2], 0, 0),
                                 b[1], (c[2], 0, 0), (d[2], 0, 0),
                                 R, tau_max=tm, dps=exchange_hp_dps)
        EXCHANGE_TAILS[(a[2], b[2], c[2], d[2])] = per
        return tot, "exchange"
    return (Q.exchange_closed_form(a[1], (a[2], 0, 0), (b[2], 0, 0),
                                   b[1], (c[2], 0, 0), (d[2], 0, 0),
                                   R, tau_max=tm),
            "exchange")


EXCHANGE_TAILS: dict = {}


def build_g(orbs, R, tau_max: int = 10, verbose: bool = False,
            exchange_hp_dps=None):
    """Exact chemist-notation g tensor as a dict {(p,q,r,s): (expr, tag)}."""
    n = len(orbs)
    out: dict = {}
    canon: dict = {}
    for p in range(n):
        for q in range(n):
            for r in range(n):
                for s in range(n):
                    key = _canon_key(p, q, r, s)
                    if key not in canon:
                        expr, tag = eri_closed_form(
                            (orbs[p], orbs[q], orbs[r], orbs[s]), R, tau_max,
                            exchange_hp_dps=exchange_hp_dps)
                        canon[key] = (expr, tag)
                        if verbose:
                            print(f"    ({p}{q}|{r}{s}) [{tag}]")
                    out[(p, q, r, s)] = canon[key]
    return out, canon


def _canon_key(p, q, r, s):
    """8-fold permutational canonical key for real orbitals."""
    cands = [(p, q, r, s), (q, p, r, s), (p, q, s, r), (q, p, s, r),
             (r, s, p, q), (s, r, p, q), (r, s, q, p), (s, r, q, p)]
    return min(cands)


# --------------------------------------------------------------- numerics (mpmath)

def _mpm(expr, dps):
    if isinstance(expr, mp.mpf):
        return +expr
    return mp.mpf(str(sp.N(expr, dps + 10)))


def evaluate(S, h, gmap, n, dps):
    with mp.workdps(dps + 15):
        Sm = mp.matrix(n, n)
        hm = mp.matrix(n, n)
        for i in range(n):
            for j in range(n):
                Sm[i, j] = _mpm(S[i][j], dps)
                hm[i, j] = _mpm(h[i][j], dps)
        gcache: dict = {}
        gm = {}
        for k, (expr, _tag) in gmap.items():
            key = id(expr)
            if key not in gcache:
                gcache[key] = _mpm(expr, dps)
            gm[k] = gcache[key]
    return Sm, hm, gm


def lowdin(Sm, n):
    """S^{-1/2}, symmetric, in mpmath."""
    E, V = mp.eigsy(Sm)
    X = mp.matrix(n, n)
    for i in range(n):
        for j in range(n):
            acc = mp.mpf(0)
            for k in range(n):
                acc += V[i, k] * V[j, k] / mp.sqrt(E[k])
            X[i, j] = acc
    return X


def transform(X, hm, gm, n):
    ht = mp.matrix(n, n)
    for i in range(n):
        for j in range(n):
            acc = mp.mpf(0)
            for p in range(n):
                for q in range(n):
                    acc += X[p, i] * X[q, j] * hm[p, q]
            ht[i, j] = acc
    # g transform, done in stages to stay O(n^5)
    g1 = {}
    for i in range(n):
        for q in range(n):
            for r in range(n):
                for s in range(n):
                    acc = mp.mpf(0)
                    for p in range(n):
                        acc += X[p, i] * gm[(p, q, r, s)]
                    g1[(i, q, r, s)] = acc
    g2 = {}
    for i in range(n):
        for j in range(n):
            for r in range(n):
                for s in range(n):
                    acc = mp.mpf(0)
                    for q in range(n):
                        acc += X[q, j] * g1[(i, q, r, s)]
                    g2[(i, j, r, s)] = acc
    g3 = {}
    for i in range(n):
        for j in range(n):
            for k in range(n):
                for s in range(n):
                    acc = mp.mpf(0)
                    for r in range(n):
                        acc += X[r, k] * g2[(i, j, r, s)]
                    g3[(i, j, k, s)] = acc
    gt = {}
    for i in range(n):
        for j in range(n):
            for k in range(n):
                for m in range(n):
                    acc = mp.mpf(0)
                    for s in range(n):
                        acc += X[s, m] * g3[(i, j, k, s)]
                    gt[(i, j, k, m)] = acc
    return ht, gt


def fci_ground(ht, gt, n_sp, n_elec):
    """Dense FCI ground-state energy in mpmath (Slater-Condon, chemist g)."""
    sos = [(p, sg) for p in range(n_sp) for sg in (0, 1)]
    dets = [tuple(sorted(c)) for c in combinations(range(2 * n_sp), n_elec)]
    dim = len(dets)

    def spat(so):
        return sos[so][0]

    def spin(so):
        return sos[so][1]

    def anti(pq, rs):
        p, q = pq
        r, s_ = rs
        val = mp.mpf(0)
        if spin(p) == spin(r) and spin(q) == spin(s_):
            val += gt[(spat(p), spat(r), spat(q), spat(s_))]
        if spin(p) == spin(s_) and spin(q) == spin(r):
            val -= gt[(spat(p), spat(s_), spat(q), spat(r))]
        return val

    H = mp.matrix(dim, dim)
    for I, di in enumerate(dets):
        occ = set(di)
        e = mp.mpf(0)
        for p in di:
            e += ht[spat(p), spat(p)]
        for ai in range(len(di)):
            for bi in range(ai + 1, len(di)):
                e += anti((di[ai], di[bi]), (di[ai], di[bi]))
        H[I, I] = e
        for J in range(I + 1, dim):
            dj = dets[J]
            occj = set(dj)
            diff_i = sorted(occ - occj)
            diff_j = sorted(occj - occ)
            nd = len(diff_i)
            if nd > 2:
                continue
            if nd == 1:
                p, q = diff_i[0], diff_j[0]
                if spin(p) != spin(q):
                    continue
                phase = (-1) ** (di.index(p) + dj.index(q))
                val = ht[spat(p), spat(q)]
                for r in di:
                    if r == p:
                        continue
                    val += anti((p, r), (q, r))
                H[I, J] = H[J, I] = phase * val
            else:
                p1, p2 = diff_i
                q1, q2 = diff_j
                phase = (-1) ** (di.index(p1) + di.index(p2)
                                 + dj.index(q1) + dj.index(q2))
                H[I, J] = H[J, I] = phase * anti((p1, p2), (q1, q2))
    ev = mp.eigsy(H, eigvals_only=True)
    return min(ev)


def total_energy(S, h, gmap, orbs, n_elec, ZA, ZB, R, dps):
    """FCI + V_NN at precision `dps`, everything in mpmath."""
    n = len(orbs)
    with mp.workdps(dps + 25):
        Sm, hm, gm = evaluate(S, h, gmap, n, dps + 15)
        X = lowdin(Sm, n)
        ht, gt = transform(X, hm, gm, n)
        e = fci_ground(ht, gt, n, n_elec)
        vnn = mp.mpf(str(sp.N(sp.nsimplify(ZA) * sp.nsimplify(ZB)
                              / sp.nsimplify(R), dps + 25)))
        tot = e + vnn
    return tot, e, vnn


# ------------------------------------------------- the closed-form H2 PES
# In the minimal 1s/1s basis (zeta=1) the singlet-Sigma_g+ CI space is
# {|g^2>, |u^2>}, so the FCI ground state is the quadratic-formula root of a
# 2x2 with closed-form entries: E(R) is ONE symbolic expression with atoms
# {exp, E_1, log, EulerGamma}.  Certified constants (50 digits, Newton on the
# exact derivative; tests/test_paper58_qfd.py):
#   R_eq = 1.667999966972748724927046410307264513349804147483
#   D_e  = 0.11865036209829531185349776433388467169190968328627
#   k    = 0.25470393430796998115272793572435567267825178990899

def h2_closed_form_E(R):
    """E(R) for minimal-basis H2 (zeta=1) as a single sympy expression."""
    oA, oB = ("A", 1, 1), ("B", 1, 1)
    S = Q.overlap(oA, oB, R)
    hAA = Q.h_core(oA, oA, 1, 1, R)
    hAB = Q.h_core(oA, oB, 1, 1, R)
    g, _ = build_g([oA, oB], R, tau_max=2)
    G = {k: v[0] for k, v in g.items()}
    g0, gAB, gx = G[(0, 0, 0, 0)], G[(0, 0, 1, 1)], G[(0, 1, 0, 1)]
    gh = G.get((0, 0, 0, 1), G.get((0, 1, 1, 1)))
    hgg = (hAA + hAB) / (1 + S)
    huu = (hAA - hAB) / (1 - S)
    Jgg = (g0 + gAB + 2 * gx + 4 * gh) / (2 * (1 + S) ** 2)
    Juu = (g0 + gAB + 2 * gx - 4 * gh) / (2 * (1 - S) ** 2)
    Kgu = (g0 - gAB) / (2 * (1 - S ** 2))
    H11 = 2 * hgg + Jgg
    H22 = 2 * huu + Juu
    return (H11 + H22) / 2 - sp.sqrt(((H11 - H22) / 2) ** 2 + Kgu ** 2) + 1 / R


def h2_pes_certify(dps_hi: int = 60, dps_lo: int = 45):
    R = sp.Symbol("R", positive=True)
    E = h2_closed_form_E(R)
    dE, d2E = sp.diff(E, R), sp.diff(E, R, 2)
    fE = sp.lambdify(R, E, modules="mpmath")
    f = sp.lambdify(R, dE, modules="mpmath")
    f2 = sp.lambdify(R, d2E, modules="mpmath")

    def newton(dps, tol):
        with mp.workdps(dps):
            x = mp.mpf("1.6")
            for _ in range(80):
                step = f(x) / f2(x)
                x -= step
                if abs(step) < tol:
                    break
            return x

    req_hi = newton(dps_hi, mp.mpf(10) ** (-dps_hi + 5))
    req_lo = newton(dps_lo, mp.mpf(10) ** (-dps_lo + 5))
    with mp.workdps(dps_hi):
        out = {
            "Req": mp.nstr(req_hi, 50),
            "Req_crosscheck_diff": mp.nstr(abs(req_hi - req_lo), 4),
            "E_Req": mp.nstr(fE(req_hi), 50),
            "k": mp.nstr(f2(req_hi), 50),
            "De": mp.nstr(mp.mpf(-1) - fE(req_hi), 50),
            "dissociation_E60_plus_1": mp.nstr(fE(mp.mpf(60)) + 1, 8),
        }
    return E, out


