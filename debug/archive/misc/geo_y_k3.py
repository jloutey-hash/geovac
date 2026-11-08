"""Geometric route to Y off the slice: the three-Fock (K3) piece.

The cut through the three densities at the origin (scales a=1, b=2 and the apex t)
carries three Fock dispersion factors.  By Paper 59 eq:tree their holomorphic period is
the three-scale Bessel moment V3(1,2,t) = int K0(x)K0(2x)K0(tx) dx, the period of the
double plane  y^2 = (p^2+a^2)(q^2+b^2)((p+q)^2+t^2)  (six lines, a K3 surface).

This script (no operator for Y or its kite is computed anywhere):
  1. puts Bailey's moduli (BBBG eq. 33) in the variable x = t^2 and shows
     lambda_pm = (x-3 +- sqrt((x-1)(x-9)))/(2x), conjugate over the conic w^2=(x-1)(x-9);
  2. rationally parametrizes the conic, x = (9-u^2)/(1-u^2), and finds
     lambda_+ = (u+1)/(3-u), lambda_- = lambda_+(-u): each curve is the Legendre family
     pulled back by a Moebius map, so the K3 piece is the tensor induction of the
     Legendre local system from the u-line to the x-line;
  3. checks Bailey against quadrature in all three sorting regimes;
  4. checks numerically that the swap-odd product combination branches like sqrt
     at t=1 and t=3 while V3 is analytic there (monodromy = swap);
  5. counts points over F_p: the character sum of the double plane, mod p, equals
     +- a_p(E_+) a_p(E_-), and a_p(E_+)^2 != a_p(E_-)^2 at some p (no isogeny, even up
     to twist), at several rational u.
"""
import json
import sympy as sp
import mpmath as mp

out = {}

# ---------------------------------------------------------------- 1. Bailey moduli in x
x, u = sp.symbols('x u')
# regime t > 3 (c = t largest, a=1, b=2): P = k+^2 + k-^2, Q = k+ k-
P1 = (x + 1 - 4) / x
Q1sq = 1 / x
disc1 = sp.factor(sp.together(P1**2 - 4 * Q1sq))
# regime 1 < t < 2 (c = 2, a = 1, b = t)
P2 = (4 + 1 - x) / 4
Q2sq = sp.Rational(1, 4)
disc2 = sp.factor(sp.together(P2**2 - 4 * Q2sq))
out['disc_regime_t_gt_3'] = str(disc1)
out['disc_regime_1_t_2'] = str(disc2)

# ---------------------------------------------------------------- 2. conic parametrization
xu = (9 - u**2) / (1 - u**2)
wu = 8 * u / (1 - u**2)
assert sp.simplify(wu**2 - (xu - 1) * (xu - 9)) == 0
lam_p = sp.factor(sp.simplify((xu - 3 + wu) / (2 * xu)))
lam_m = sp.factor(sp.simplify((xu - 3 - wu) / (2 * xu)))
out['lambda_plus_u'] = str(lam_p)
out['lambda_minus_u'] = str(lam_m)
assert sp.simplify(lam_m - lam_p.subs(u, -u)) == 0


def jleg(l):
    return sp.factor(256 * (l**2 - l + 1)**3 / (l**2 * (l - 1)**2))


jp, jm = jleg(lam_p), jleg(lam_m)
# symmetric functions of j_pm are rational in x: verify by expressing in x
s_sum = sp.factor(sp.simplify(jp + jm))
s_prod = sp.factor(sp.simplify(jp * jm))
out['j_plus_u'] = str(jp)
# singular fibres of E_+ on the u-line
out['E_plus_bad_u'] = {'lambda=0': str(sp.solve(sp.numer(lam_p), u)),
                       'lambda=1': str(sp.solve(sp.numer(sp.together(lam_p - 1)), u)),
                       'lambda=inf': str(sp.solve(sp.denom(lam_p), u))}
out['x_at_u'] = {str(v): str(sp.limit(xu, u, v)) for v in [-3, -1, 0, 1, 3, sp.oo]}

# Legendre-data of each conjugate at the two preimages of x=0 (u=+-3) and x=inf (u=+-1)
loc = {}
for uv in [3, -3, 1, -1, 0, sp.oo]:
    lp = sp.limit(lam_p, u, uv)
    lm = sp.limit(lam_m, u, uv)
    loc[str(uv)] = (str(lp), str(lm))
out['lambda_pm_at_preimages'] = loc

# ---------------------------------------------------------------- 3. Bailey vs quadrature
mp.mp.dps = 30


def V3_quad(a, b, c):
    g = lambda s: mp.besselk(0, a * s) * mp.besselk(0, b * s) * mp.besselk(0, c * s)
    return mp.quad(g, [0, mp.mpf(1) / 16, mp.mpf(1) / 2, 2, 8, mp.inf])


def V3_bailey(a, b, c):
    a, b, c = sorted([mp.mpf(a), mp.mpf(b), mp.mpf(c)])
    s1 = mp.sqrt((c + a)**2 - b**2)
    s2 = mp.sqrt((c - a)**2 - b**2 + 0j)
    kp, km = (s1 + s2) / (2 * c), (s1 - s2) / (2 * c)
    K = lambda k: mp.ellipk(k * k)
    Kc = lambda k: mp.ellipk(1 - k * k)
    return mp.pi / (2 * c) * (K(km) * Kc(kp) + K(kp) * Kc(km))


chk = {}
for t in ['0.5', '1.5', '2.5', '4']:
    q = V3_quad(1, 2, mp.mpf(t))
    bl = V3_bailey(1, 2, mp.mpf(t))
    chk[t] = float(abs(q - bl) / abs(q))
out['bailey_rel_err'] = chk

# ---------------------------------------------------------------- 4. swap monodromy at t=3, t=1
# In regime t>3 write the four products; the swap-odd one is D = K(k-)K(k+') - K(k+)K(k-').
# Approach t -> 3+ : if the monodromy is the swap, D ~ sqrt(t-3) * analytic, and the swap-even
# combinations are analytic.  Fit D(3+e)/sqrt(e) at shrinking e: it should tend to a constant.


def prods(t, regime):
    t = mp.mpf(t)
    if regime == 'hi':
        a, b, c = mp.mpf(1), mp.mpf(2), t
    else:
        a, b, c = sorted([mp.mpf(1), t, mp.mpf(2)])
    s1 = mp.sqrt((c + a)**2 - b**2)
    s2 = mp.sqrt((c - a)**2 - b**2 + 0j)
    kp, km = (s1 + s2) / (2 * c), (s1 - s2) / (2 * c)
    K = lambda k: mp.ellipk(k * k)
    Kc = lambda k: mp.ellipk(1 - k * k)
    even = K(km) * Kc(kp) + K(kp) * Kc(km)
    odd = K(km) * Kc(kp) - K(kp) * Kc(km)
    return even, odd


sw = {}
for e in ['1e-4', '1e-6', '1e-8']:
    ev, od = prods(3 + mp.mpf(e), 'hi')
    sw['t=3+' + e] = {'odd/sqrt(e)': mp.nstr(od / mp.sqrt(mp.mpf(e)), 12),
                      'even': mp.nstr(ev, 15)}
for e in ['1e-4', '1e-6', '1e-8']:
    # regime t<1: sorted (t,1,2); t = 1 - e ; s2 = sqrt((2-t)^2-1) -> vanishes at t=1
    ev, od = prods(1 - mp.mpf(e), 'lo')
    sw['t=1-' + e] = {'odd/sqrt(e)': mp.nstr(od / mp.sqrt(mp.mpf(e)), 12),
                      'even': mp.nstr(ev, 15)}
out['swap_check'] = sw

# ---------------------------------------------------------------- 5. point counts over F_p


def chi(v, p):
    v %= p
    if v == 0:
        return 0
    return 1 if pow(v, (p - 1) // 2, p) == 1 else -1


def ap_legendre(lam_num, lam_den, p):
    # E: y^2 = X(X-1)(X-lam);  a_p = -sum chi
    if lam_den % p == 0:
        return None
    lam = lam_num * pow(lam_den, p - 2, p) % p
    if lam in (0, 1):
        return None
    return -sum(chi(X * (X - 1) * (X - lam), p) for X in range(p))


def k3_charsum(xn, xd, p):
    # S = sum_{P,Q} chi((P^2+1)(Q^2+4)((P+Q)^2+x)), x = xn/xd
    if xd % p == 0:
        return None
    xv = xn * pow(xd, p - 2, p) % p
    A = [chi(P * P + 1, p) for P in range(p)]
    B = [chi(Q * Q + 4, p) for Q in range(p)]
    C = [chi(R * R + xv, p) for R in range(p)]
    s = 0
    for P in range(p):
        if A[P] == 0:
            continue
        aP = A[P]
        for Q in range(p):
            s += aP * B[Q] * C[(P + Q) % p]
    return s


primes = [p for p in sp.primerange(11, 140)]
pc = {}
for uv in [sp.Rational(1, 2), sp.Rational(5, 7), sp.Rational(7, 3), sp.Rational(2, 5)]:
    xv = sp.nsimplify(xu.subs(u, uv))
    lp = sp.nsimplify(lam_p.subs(u, uv))
    lm = sp.nsimplify(lam_m.subs(u, uv))
    rows = []
    for p in primes:
        if any(sp.fraction(v)[1] % p == 0 or sp.fraction(v)[0] % p == 0 for v in [xv, xv - 1, xv - 9, lp, lm, lp - 1, lm - 1]):
            continue
        apl = ap_legendre(int(sp.fraction(lp)[0]), int(sp.fraction(lp)[1]), p)
        aml = ap_legendre(int(sp.fraction(lm)[0]), int(sp.fraction(lm)[1]), p)
        S = k3_charsum(int(sp.fraction(xv)[0]), int(sp.fraction(xv)[1]), p)
        if apl is None or aml is None or S is None:
            continue
        prod = apl * aml
        rows.append({'p': p, 'S': S, 'S_mod_p': S % p, 'ap+': apl, 'ap-': aml,
                     'match_plus': (S - 2 - prod) % p == 0, 'match_minus': (S - 2 + prod) % p == 0})
    pc[str(uv)] = {'x': str(xv), 'lam+': str(lp), 'lam-': str(lm), 'rows': rows,
                   'all_match_pm': all(r['match_plus'] or r['match_minus'] for r in rows),
                   'n_rows': len(rows),
                   'n_trivial': sum(1 for r in rows if (r['ap+'] * r['ap-']) % r['p'] == 0),
                   'ap2_differ_at': [r['p'] for r in rows if r['ap+']**2 != r['ap-']**2]}
out['point_counts'] = pc

with open('debug/data/geo_y_k3.json', 'w') as f:
    json.dump(out, f, indent=1, default=str)

for k, v in out.items():
    if k == 'point_counts':
        for uu, d in v.items():
            print('u =', uu, 'x =', d['x'], 'lam+ =', d['lam+'], 'lam- =', d['lam-'],
                  '| rows', d['n_rows'], 'all S = +-a+a- mod p:', d['all_match_pm'],
                  '| trivial rows', d['n_trivial'], '| a+^2 != a-^2 at', d['ap2_differ_at'][:6])
    else:
        print(k, ':', v)
