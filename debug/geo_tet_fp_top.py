"""Blind geometric route, part 2: point counts on the top-sector cut of the tetrahedron.

Baikov form of the 3-loop vacuum tetrahedron in D = 3 (exponent of the Gram/Cayley-Menger
factor = (D-L-1)/2 = -1/2).  The top-sector maximal cut puts the three integer lines on
shell (Coulomb d_C = 0, geminals d_g1 = -g1^2, d_g2 = -g2^2) and keeps the three density
lines with their half-integer powers.  Its cohomology is that of the double cover

    V_x :  y^2 = A B T f(A,B,T),   f = 4 det G_cut,
    A = d_a + a^2, B = d_b + b^2, T = d_t + x   (the three Fock factors),

over (A,B,T) in A^3.  (a,b,g1,g2) = (1,2,3/2,3).  We compute

  S(p,x) = sum_{A,B,T in F_p} chi(A B T f)        (chi = Legendre symbol),

and compare with the closed form derived from the linearity of f in T (geo_tet_symbolic S2):

  S = -chi(-1) + p * sum_{mu} ( -chi(-1) - chi(-mu alpha beta) ),
  mu over the F_p-roots of g1^2 mu^2 - (g1^2+g2^2-x) mu + g2^2, alpha = g1^2-a^2, beta = g2^2-b^2,

which contains only finite-order characters times powers of p (Artin-Tate): no trace of
Frobenius of any curve or surface.  Control: the same cut with a MASSIVE line opposite the
apex (d_C = -m^2) makes f quadratic in T and the fibre an elliptic curve.

Also: the Euler characteristic of the cut complement {A B T f != 0} (number of top-sector
masters for generic exponents), from point counts at split primes.
No differential operator for Y or for any sector is computed.
"""
import json
from fractions import Fraction as Fr
import numpy as np

a2, b2, g1s, g2s = Fr(1), Fr(4), Fr(9, 4), Fr(9)


def primes(lo, hi):
    out = []
    for n in range(lo, hi):
        if n > 1 and all(n % d for d in range(2, int(n**0.5) + 1)):
            out.append(n)
    return out


def red(q, p):
    q = Fr(q)
    return (q.numerator * pow(q.denominator, -1, p)) % p


def legendre_table(p):
    t = np.zeros(p, dtype=np.int64)
    sq = (np.arange(1, p, dtype=np.int64) ** 2) % p
    t[np.unique(sq)] = 1
    t[(t == 0) & (np.arange(p) != 0)] = -1
    return t


def topcut_sum(p, x, mC2=Fr(0), count=False):
    chi = legendre_table(p)
    A = np.arange(p, dtype=np.int64)
    G1, G2 = red(g1s, p), red(g2s, p)
    xa = red(x, p)
    P = (A - red(a2, p) + G1) % p           # P = d_a + g1^2
    Q = (A - red(b2, p) + G2) % p           # Q = d_b + g2^2 (indexed by B)
    S = 0
    N = 0
    m = red(mC2, p)
    for iA in range(p):
        Pv = int(P[iA])
        # f as a function of (B, T): 4 det G with d_C = -mC2, d_g = -g^2
        Bv = A[:, None]
        Tv = A[None, :]
        Qv = Q[:, None]
        dt = (Tv - xa) % p
        da = (iA - red(a2, p)) % p
        db = (Bv - red(b2, p)) % p
        if m == 0:
            f = (G2 * Pv * Pv + G1 * Qv * Qv - (G1 + G2 + dt) * Pv * Qv) % p
        else:
            # general 4 det G with d_C = -m, d_g1 = -G1, d_g2 = -G2
            dC = (-m) % p
            dg1, dg2 = (-G1) % p, (-G2) % p
            u = dC
            uv1 = (dC + da - dg1) % p
            uv2 = (dC + db - dg2) % p
            v1v2 = (da + db - dt) % p
            # 8 det G = 8[u(d_a d_b - v1v2^2/4) - (uv1/2)(uv1 d_b/2 - v1v2 uv2/4) + (uv2/2)(uv1 v1v2/4 - d_a uv2/2)]
            #  -> 4 det G * 2 computed with the half-integers cleared: use G2 = 4*det(2*Gram)/... simpler:
            # det of the integer matrix 2*Gram = 8 det G
            g11, g12, g13 = (2 * u) % p, uv1, uv2
            g22, g23, g33 = (2 * da) % p, v1v2, (2 * db) % p
            f = (g11 * (g22 * g33 - g23 * g23) - g12 * (g12 * g33 - g23 * g13)
                 + g13 * (g12 * g23 - g22 * g13)) % p
        val = (iA * Bv % p) * Tv % p * f % p
        S += int(chi[val].sum())
        if count:
            N += int((val != 0).sum())
    return S, N


def closed_form(p, x):
    chi = legendre_table(p)
    G1, G2 = red(g1s, p), red(g2s, p)
    al, be = red(g1s - a2, p), red(g2s - b2, p)
    xa = red(x, p)
    c1 = (G1 + G2 - xa) % p
    roots = [mu for mu in range(p) if (G1 * mu * mu - c1 * mu + G2) % p == 0]
    S = -int(chi[p - 1])
    for mu in roots:
        S += p * (-int(chi[p - 1]) - int(chi[(-mu * al * be) % p]))
    return S, len(roots)


out = {'rows': [], 'control_massive': [], 'euler': []}
xs = [Fr(7, 3), Fr(5), Fr(13, 2), Fr(31, 5), Fr(2, 7), Fr(-3, 4)]
bad_den = 2 * 3 * 5 * 7
match = 0
total = 0
for p in primes(11, 160):
    if bad_den % p == 0:
        continue
    for x in xs:
        xa = red(x, p)
        if xa in (0, red(Fr(9, 4), p), red(Fr(81, 4), p)):
            continue
        S, _ = topcut_sum(p, x)
        C, nr = closed_form(p, x)
        total += 1
        match += int(S == C)
        out['rows'].append({'p': p, 'x': str(x), 'S': S, 'closed': C, 'n_mu_roots': nr,
                            'S_mod_p_centered': ((S + p // 2) % p) - p // 2})
out['closed_form_matches'] = f'{match}/{total}'

# control: massive line opposite the apex.  Mod p both are constant; the difference is at
# the next p-adic level: q = (S - r)/p with r the centred residue.  Artin-Tate data keep q a
# bounded integer combination of characters; a curve's Frobenius trace makes it wander ~ sqrt p.
for p in primes(11, 110):
    if bad_den % p == 0:
        continue
    row = {'p': p, 'massive': [], 'massless': []}
    for x in xs:
        for key, m2 in (('massive', Fr(5, 3)), ('massless', Fr(0))):
            S, _ = topcut_sum(p, x, mC2=m2)
            r = ((S + p // 2) % p) - p // 2
            q1 = (S - r) // p
            r1 = ((q1 + p // 2) % p) - p // 2
            row[key].append({'x': str(x), 'S': S, 'r': r, 'q_centred_mod_p': r1})
    out['control_massive'].append(row)

# Euler characteristic of {A B T f != 0} at primes where the two mu-lines split
for x in [Fr(7, 3), Fr(5)]:
    pts = []
    for p in primes(11, 120):
        if bad_den % p == 0:
            continue
        xa = red(x, p)
        if xa in (0, red(Fr(9, 4), p), red(Fr(81, 4), p)):
            continue
        _, nr = closed_form(p, x)
        if nr != 2:
            continue
        _, N = topcut_sum(p, x, count=True)
        pts.append((p, N))
    ps = np.array([q for q, _ in pts], dtype=float)
    Ns = np.array([n for _, n in pts], dtype=float)
    coef = np.polyfit(ps, Ns, 3)
    resid = float(np.max(np.abs(np.polyval(coef, ps) - Ns)))
    icoef = [int(round(c)) for c in coef]
    exact = all(np.polyval(icoef, q) == n for q, n in pts)
    out['euler'].append({'x': str(x), 'n_primes': len(pts), 'poly_coeffs_p3_to_p0': icoef,
                         'fit_resid': resid, 'integer_poly_exact_at_all_primes': exact,
                         'chi = N(1)': int(np.polyval(icoef, 1))})

with open('debug/data/geo_tet_fp_top.json', 'w') as fh:
    json.dump(out, fh, indent=1)
print('closed-form matches', out['closed_form_matches'])
for r in out['rows'][:12]:
    print(r)
print('control (S mod p, centred), massive vs massless line opposite the apex:')
for r in out['control_massive']:
    print(r['p'], 'massive q:', [e['q_centred_mod_p'] for e in r['massive']],
          '| massless q:', [e['q_centred_mod_p'] for e in r['massless']])
print('euler', out['euler'])
