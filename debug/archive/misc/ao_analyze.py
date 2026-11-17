"""Analyse an operator found by debug/ao_guess.py (record AO, docs/predictions/2026-10-28_angle_operator.md).

Input: the null vector a_ij of  sum a_ij x^j T^i,  T = (1 - x^2) d/dx.  Output: the standard form
sum_k c_k(x) d^k/dx^k (polynomial coefficients), its finite singular points (roots of c_R), the local exponents at each
(indicial polynomial; 'irregular' if the Fuchs condition fails), and the Fuchs test at infinity.
"""
import mpmath as mp


def padd(p, q):
    n = max(len(p), len(q))
    return [(p[i] if i < len(p) else 0) + (q[i] if i < len(q) else 0) for i in range(n)]


def pmul(p, q):
    out = [mp.mpf(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i + j] += a * b
    return out


def pder(p):
    return [i * p[i] for i in range(1, len(p))] or [mp.mpf(0)]


ONE_MINUS_X2 = [mp.mpf(1), mp.mpf(0), mp.mpf(-1)]


def T_powers(R):
    """T^i as {k: poly} for i = 0..R."""
    ops = [{0: [mp.mpf(1)]}]
    for i in range(1, R + 1):
        prev, new = ops[-1], {}
        for k, b in prev.items():
            new[k] = padd(new.get(k, [mp.mpf(0)]), pmul(ONE_MINUS_X2, pder(b)))
            new[k + 1] = padd(new.get(k + 1, [mp.mpf(0)]), pmul(ONE_MINUS_X2, b))
        ops.append(new)
    return ops


def standard_form(keys, v):
    R = max(i for i, _ in keys)
    ops = T_powers(R)
    c = {}
    for (i, j), a in zip(keys, v):
        xj = [mp.mpf(0)] * j + [a]
        for k, b in ops[i].items():
            c[k] = padd(c.get(k, [mp.mpf(0)]), pmul(xj, b))
    return [c.get(k, [mp.mpf(0)]) for k in range(R + 1)]


def trim(p, tol):
    s = max(abs(a) for a in p)
    p = list(p)
    while len(p) > 1 and abs(p[-1]) < tol * s:
        p.pop()
    return p


def taylor_at(p, x0):
    """Coefficients of p in powers of (x - x0)."""
    n = len(p)
    out = []
    q = list(p)
    for m in range(n):
        out.append(mp.polyval(q[::-1], x0) / mp.factorial(m))
        q = pder(q)
    return out


def order_at(tc, tol, scale):
    for m, a in enumerate(tc):
        if abs(a) > tol * scale:
            return m
    return len(tc)


def exponents_at(c, x0, tol, mult):
    """mult = multiplicity of x0 as a root of c_R (from the root cluster); orders of the other c_k by tolerance."""
    R = len(c) - 1
    tcs = [taylor_at(ck, x0) for ck in c]
    scale = max(max(abs(a) for a in t) for t in tcs)
    ords = [order_at(t, tol, scale) for t in tcs]
    ords[R] = mult
    delta = ords[R] - R
    if any(ords[k] - k < delta for k in range(R)):
        return 'irregular', ords
    # indicial polynomial  sum_k lead_k e(e-1)...(e-k+1)
    poly = [mp.mpf(0)]
    for k in range(R + 1):
        lead = tcs[k][k + delta] if 0 <= k + delta < len(tcs[k]) else 0
        fall = [mp.mpf(1)]
        for m in range(k):
            fall = pmul(fall, [mp.mpf(-m), mp.mpf(1)])
        poly = padd(poly, [lead * f for f in fall])
    poly = trim(poly, tol)
    roots = mp.polyroots(poly[::-1], maxsteps=200, extraprec=200) if len(poly) > 1 else []
    return sorted(roots, key=lambda z: mp.re(z)), ords


def fuchs_at_infinity(c, tol):
    degs = [len(trim(ck, tol)) - 1 for ck in c]
    R = len(c) - 1
    return all(degs[k] - k <= degs[R] - R for k in range(R + 1)), degs


def analyse(keys, v, tol=None):
    tol = tol or mp.mpf(10) ** (-mp.mp.dps // 3)
    c = standard_form(keys, v)
    c = [trim(ck, tol) for ck in c]
    lead = c[-1]
    roots = mp.polyroots(lead[::-1], maxsteps=400, extraprec=400)
    # cluster near-equal roots
    clusters = []
    for z in roots:
        for cl in clusters:
            if abs(z - cl[0]) < mp.mpf(10) ** (-mp.mp.dps // 8):
                cl.append(z)
                break
        else:
            clusters.append([z])
    out = []
    for cl in clusters:
        z = sum(cl) / len(cl)
        ex, ords = exponents_at(c, z, tol, len(cl))
        out.append(dict(x=z, mult=len(cl), exponents=ex, orders=ords))
    reg_inf, degs = fuchs_at_infinity(c, tol)
    return dict(R=len(c) - 1, degrees=degs, points=out, regular_at_infinity=reg_inf)
