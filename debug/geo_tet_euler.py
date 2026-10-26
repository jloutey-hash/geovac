"""Blind geometric route, part 2: Euler characteristics of the sector cuts of the tetrahedron.

For each sector (the set S of integer lines kept: Coulomb C, geminals g1, g2), the cut
variety of the 3-loop Baikov representation is the complement

    U_S = { (A, B, T, d_e for e not in S) :  A B T f != 0 },
    f = det(2 Gram)  with the lines of S on shell, the others free (ISPs),

and |chi(U_S)| counts the sector's master integrals for non-integer density exponents
(Lee-Pomeransky / Bitoun et al. criterion).  chi is read from point counts: #U_S(F_p)
is fitted by an integer polynomial at primes where it is polynomial, and evaluated at 1.
The sunrise (no integer line) is also counted in its 2-loop form, C^3 minus
{A B T lambda = 0}, lambda the Kallen function of (d_a, d_b, d_t).
(a,b,g1,g2) = (1,2,3/2,3), x = 7/3.  Geometry only: no operator is computed.
"""
import json
import sys
from fractions import Fraction as Fr
import numpy as np

a2, b2, g1s, g2s = Fr(1), Fr(4), Fr(9, 4), Fr(9)
X = Fr(7, 3)


def primes(lo, hi):
    return [n for n in range(lo, hi) if n > 1 and all(n % d for d in range(2, int(n**0.5) + 1))]


def red(q, p):
    q = Fr(q)
    return (q.numerator * pow(q.denominator, -1, p)) % p


def det2G(dC, da, db, dt, dg1, dg2, p):
    g11, g12, g13 = (2 * dC) % p, (dC + da - dg1) % p, (dC + db - dg2) % p
    g22, g23, g33 = (2 * da) % p, (da + db - dt) % p, (2 * db) % p
    return (g11 * (g22 * g33 - g23 * g23) - g12 * (g12 * g33 - g23 * g13)
            + g13 * (g12 * g23 - g22 * g13)) % p


def count_sector(S, p):
    """#{A B T f != 0} with lines in S on shell and the other integer lines free."""
    xa, A2, B2 = red(X, p), red(a2, p), red(b2, p)
    on = {'C': 0, 'g1': (-red(g1s, p)) % p, 'g2': (-red(g2s, p)) % p}
    free = [e for e in ('C', 'g1', 'g2') if e not in S]
    r = np.arange(p, dtype=np.int64)
    total = 0
    # outer loops: A and the free integer lines; inner (B, T) vectorized
    import itertools
    for Av in range(1, p):
        da = (Av - A2) % p
        for fv in itertools.product(range(p), repeat=len(free)):
            vals = dict(on)
            vals.update(dict(zip(free, fv)))
            Bv = r[:, None]
            Tv = r[None, :]
            db = (Bv - B2) % p
            dt = (Tv - xa) % p
            f = det2G(vals['C'], da, db, dt, vals['g1'], vals['g2'], p)
            ok = (Bv != 0) & (Tv != 0) & (f != 0)
            total += int(ok.sum())
    return total


def count_sunrise_2loop(p):
    xa, A2, B2 = red(X, p), red(a2, p), red(b2, p)
    r = np.arange(p, dtype=np.int64)
    total = 0
    for Av in range(1, p):
        da = (Av - A2) % p
        Bv = r[:, None]
        Tv = r[None, :]
        db = (Bv - B2) % p
        dt = (Tv - xa) % p
        lam = (da * da + db * db + dt * dt - 2 * da * db - 2 * db * dt - 2 * dt * da) % p
        ok = (Bv != 0) & (Tv != 0) & (lam != 0)
        total += int(ok.sum())
    return total


def fit(pts, deg):
    ps = np.array([q for q, _ in pts], dtype=float)
    Ns = np.array([n for _, n in pts], dtype=float)
    coef = np.polyfit(ps, Ns, deg)
    ic = [int(round(c)) for c in coef]
    exact = [q for q, n in pts if int(np.polyval(ic, q)) == n]
    return ic, len(exact), len(pts), int(np.polyval(ic, 1))


jobs = {
    'top {g1,g2,C}': (('C', 'g1', 'g2'), 3, primes(11, 120)),
    '{g1,g2}': (('g1', 'g2'), 4, primes(11, 60)),
    '{g1,C}': (('C', 'g1'), 4, primes(11, 60)),
    '{g2,C}': (('C', 'g2'), 4, primes(11, 60)),
    '{g1}': (('g1',), 5, primes(11, 48)),
}
which = sys.argv[1:] or list(jobs) + ['sunrise_2loop']
out = {}
bad = 2 * 3 * 5 * 7
for name in which:
    if name == 'sunrise_2loop':
        pts = [(p, count_sunrise_2loop(p)) for p in primes(11, 150) if bad % p]
        deg = 3
    else:
        S, deg, ps = jobs[name]
        pts = [(p, count_sector(S, p)) for p in ps if bad % p]
    # split by p mod 4 and mod 3 to expose Artin dependence
    res = {}
    for label, cond in [('all', lambda q: True), ('p=1 mod 4', lambda q: q % 4 == 1),
                        ('p=3 mod 4', lambda q: q % 4 == 3)]:
        sub = [(q, n) for q, n in pts if cond(q)]
        if len(sub) > deg + 1:
            ic, ne, nt, chi = fit(sub, deg)
            res[label] = {'coeffs': ic, 'exact_at': f'{ne}/{nt}', 'chi=N(1)': chi}
    out[name] = {'counts': pts, 'fits': res}
    print(name, json.dumps(res), flush=True)
with open('debug/data/geo_tet_euler_%s.json' % ('all' if not sys.argv[1:] else '_'.join(
        w.replace('{', '').replace('}', '').replace(',', '').replace(' ', '') for w in which)), 'w') as fh:
    json.dump(out, fh, indent=1)
