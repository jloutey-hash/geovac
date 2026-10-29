"""Exact (mod p) analysis of an operator found by debug/na_modguess.py (record NA).
Null vector a_ij of sum a_ij x^j T^i  ->  standard form sum_k c_k(x) d^k/dx^k (nmod_poly coefficients);
at a rational point x0: orders of vanishing of each c_k (exact), the Fuchs condition, the indicial polynomial and
which candidate exponents are its roots;  at infinity: the Fuchs condition on the degrees.
"""
import flint


def T_powers(R, p):
    one_minus_x2 = flint.nmod_poly([1, 0, -1], p)
    ops = [{0: flint.nmod_poly([1], p)}]
    for _ in range(R):
        prev, new = ops[-1], {}
        for k, b in prev.items():
            new[k] = new.get(k, flint.nmod_poly([], p)) + one_minus_x2 * b.derivative()
            new[k + 1] = new.get(k + 1, flint.nmod_poly([], p)) + one_minus_x2 * b
        ops.append(new)
    return ops


def standard_form(keys, v, p):
    R = max(i for i, _ in keys)
    ops = T_powers(R, p)
    c = [flint.nmod_poly([], p) for _ in range(R + 1)]
    for (i, j), a in zip(keys, v):
        xj = flint.nmod_poly([0] * j + [int(a)], p)
        for k, b in ops[i].items():
            c[k] += xj * b
    return c


def shift(poly, x0, p):
    """Coefficients of poly(x0 + t) in t."""
    out = []
    q = poly
    fact = flint.nmod(1, p)
    m = 0
    while True:
        out.append(q(x0) / fact)
        if q.degree() <= 0:
            break
        q = q.derivative()
        m += 1
        fact *= m
    return out


def local(c, x0, p, candidates):
    """At x0 (nmod): orders, regular?, indicial polynomial values at each candidate exponent (zero = root)."""
    R = len(c) - 1
    tcs = [shift(ck, x0, p) for ck in c]
    def order(t):
        for m, a in enumerate(t):
            if int(a) != 0:
                return m
        return 10 ** 9
    ords = [order(t) for t in tcs]
    delta = ords[R] - R
    regular = all(ords[k] - k >= delta for k in range(R))
    if not regular:
        return dict(orders=ords, regular=False)
    def indicial(e):
        s = flint.nmod(0, p)
        for k in range(R + 1):
            lead = tcs[k][k + delta] if 0 <= k + delta < len(tcs[k]) else flint.nmod(0, p)
            f = flint.nmod(1, p)
            for m in range(k):
                f *= (e - m)
            s += lead * f
        return s
    roots = {name: int(indicial(e)) == 0 for name, e in candidates.items()}
    # indicial polynomial as nmod_poly in e
    poly = flint.nmod_poly([], p)
    for k in range(R + 1):
        lead = tcs[k][k + delta] if 0 <= k + delta < len(tcs[k]) else flint.nmod(0, p)
        f = flint.nmod_poly([1], p)
        for m in range(k):
            f *= flint.nmod_poly([-m, 1], p)
        poly += f * flint.nmod_poly([int(lead)], p)
    return dict(orders=ords, regular=True, delta=delta, candidate_roots=roots, indicial=poly)


def fuchs_at_infinity(c):
    R = len(c) - 1
    degs = [ck.degree() for ck in c]
    return all(degs[k] - k <= degs[R] - R for k in range(R + 1)), degs


def valuation(poly, f):
    m = 0
    q = poly
    if q.degree() < 0:
        return 10 ** 9, q
    while True:
        d, r = divmod(q, f)
        if r.degree() >= 0:
            return m, q
        q = d
        m += 1


def all_points(c, p, Emax=60):
    """For every irreducible factor f of the leading coefficient: multiplicity, regular?, and the non-negative integer
    exponents e <= Emax that are roots of the indicial polynomial (computed in F_p[t]/f; leading Taylor coefficient of
    c_k at a root theta is g_k(theta) f'(theta)^m_k with c_k = f^m_k g_k).  'apparent-type' = regular with R distinct
    integer exponents in [0, Emax]."""
    R = len(c) - 1
    out = []
    for f, mult in c[R].factor()[1]:
        fp = f.derivative()
        vals = [valuation(ck, f) for ck in c]
        ords = [v[0] for v in vals]
        delta = ords[R] - R
        regular = all(ords[k] - k >= delta for k in range(R))
        rec = dict(factor=f, degree=f.degree(), multiplicity=mult, regular=regular)
        if regular:
            lead = []
            for k in range(R + 1):
                if ords[k] == k + delta:
                    g = vals[k][1] % f
                    lead.append((g * pow(fp, ords[k], f) if False else (g * (fp ** ords[k] % f)) % f))
                else:
                    lead.append(flint.nmod_poly([], p))
            ints = []
            for e in range(0, Emax + 1):
                s = flint.nmod_poly([], p)
                for k in range(R + 1):
                    fall = 1
                    for m in range(k):
                        fall *= (e - m)
                    s += lead[k] * flint.nmod(fall, p)
                if (s % f).degree() < 0:
                    ints.append(e)
            rec['integer_exponents'] = ints
            rec['apparent_type'] = len(ints) == R
        out.append(rec)
    return out
