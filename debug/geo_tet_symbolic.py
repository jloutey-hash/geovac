"""Blind geometric route, part 2: the whole tetrahedron, including the top sector.

Exact symbolic facts (sympy).  No differential operator for Y, its kite or any subsector
is computed anywhere.  Notation: Baikov / dual-tetrahedron squared lengths d_e of the six
lines (Euclidean momenta squared); a line is on shell at d_e = -m_e^2.

Dual tetrahedron P, R, Q1, Q2 with edges
   C = P-R (Coulomb, massless), a = P-Q1, b = P-Q2, t = Q1-Q2 (apex density),
   g1 = R-Q1, g2 = R-Q2.
Opposite pairs: (C, t), (a, g2), (b, g1) -- the same pairing as the position-space
tetrahedron O,1,2,3.

  S1  Gram (= CM/288 up to a constant) of the three loop momenta; d_e^2 coefficient is
      -d_opp(e)/4: the apex line's square appears only through the Coulomb line.
  S2  Top-sector maximal cut (C, g1, g2 on shell): 4 det G is the binary quadratic
      g2^2 P^2 + g1^2 Q^2 - (g1^2+g2^2+d_t) P Q, P = d_a+g1^2, Q = d_b+g2^2, LINEAR in d_t.
  S3  Leading Landau singularity (all six on shell): det G = 0 is linear in x with root
      x_LS = -(g1^2 b^2 - g2^2 a^2)(g1^2-a^2-g2^2+b^2)/((g1^2-a^2)(g2^2-b^2)); zero at
      (a,b,g1,g2) = (1,2,3/2,3) because g1 b = g2 a.
  S4  The D = 3 one-loop ring (g1, g2, 0) in the Davydychev-Delbourgo/solid-angle form:
      D^2 + det M = P'Q'(w^2+(g1+g2)^2)/4 (log divisor = three planes), det M the same
      binary quadratic as S2; on the Gram boundary (collinear momenta) det M is a perfect
      square and D +- i sqrt(det M) split into linear forms.
  S5  Restriction of the six-line Fock double plane to the strata of S4: genus of each.
  S6  The apex-line integral of the top cut is exact (a total derivative) on the cut.
"""
import json
import sympy as sp

out = {}
dC, da, db, dt, dg1, dg2 = sp.symbols('d_C d_a d_b d_t d_g1 d_g2')
a, b, g1, g2, x, mC = sp.symbols('a b g1 g2 x m_C', positive=True)

# ---------------------------------------------------------------- S1 Gram of the dual tetrahedron
uu = dC
v1v1, v2v2 = da, db
uv1 = (dC + da - dg1) / 2
uv2 = (dC + db - dg2) / 2
v1v2 = (da + db - dt) / 2
G = sp.Matrix([[uu, uv1, uv2], [uv1, v1v1, v1v2], [uv2, v1v2, v2v2]])
detG = sp.expand(G.det())
opp = {dC: dt, dt: dC, da: dg2, dg2: da, db: dg1, dg1: db}
s1 = {}
for e, o in opp.items():
    c2 = sp.expand(sp.Poly(detG, e).coeff_monomial(e**2))
    s1[str(e)] = str(c2)
    assert sp.simplify(c2 + o / 4) == 0
out['S1_coeff_of_square'] = s1

# ---------------------------------------------------------------- S2 top cut
P, Q, T = sp.symbols('P Q T')
cut = {dC: 0, dg1: -g1**2, dg2: -g2**2}
detG_cut = sp.expand(detG.subs(cut))
sub = {da: P - g1**2, db: Q - g2**2}
f2 = sp.expand(4 * detG_cut.subs(sub))
target = sp.expand(g2**2 * P**2 + g1**2 * Q**2 - (g1**2 + g2**2 + dt) * P * Q)
assert sp.simplify(f2 - target) == 0
out['S2_4detG_topcut'] = str(sp.factor(f2))
out['S2_degree_in_dt'] = sp.Poly(f2, dt).degree()
# slopes of the two lines Q = mu P at d_t = -x (apex Fock factor T = d_t + x = 0)
mu = sp.symbols('mu')
quad_mu = sp.expand(f2.subs(Q, mu * P).subs(dt, -x) / P**2)
disc_mu = sp.factor(sp.discriminant(quad_mu, mu))
out['S2_mu_quadratic_at_T0'] = str(quad_mu)
out['S2_mu_discriminant'] = str(disc_mu)
# with massive Coulomb the cut is quadratic in d_t (control)
detG_cut_m = sp.expand(detG.subs({dC: -mC**2, dg1: -g1**2, dg2: -g2**2}))
out['S2_massiveC_degree_in_dt'] = sp.Poly(detG_cut_m, dt).degree()
out['S2_massiveC_dt2_coeff'] = str(sp.Poly(detG_cut_m, dt).coeff_monomial(dt**2))

# ---------------------------------------------------------------- S3 leading singularity
onshell = {dC: 0, dg1: -g1**2, dg2: -g2**2, da: -a**2, db: -b**2, dt: -x}
ls = sp.factor(sp.expand(4 * detG.subs(onshell)))
out['S3_4detG_all_onshell'] = str(ls)
xls = sp.solve(sp.Eq(ls, 0), x)
xls_f = [sp.factor(r) for r in xls]
out['S3_x_LS'] = [str(r) for r in xls_f]
pt = {a: 1, b: 2, g1: sp.Rational(3, 2), g2: 3}
out['S3_x_LS_at_point'] = [str(sp.nsimplify(r.subs(pt))) for r in xls_f]
out['S3_4detG_at_point'] = str(sp.factor(ls.subs(pt)))
# generic control point
pt2 = {a: 1, b: 2, g1: sp.Rational(3, 2), g2: sp.Rational(5, 2)}
out['S3_x_LS_generic_control'] = [str(sp.nsimplify(r.subs(pt2))) for r in xls_f]

# ---------------------------------------------------------------- S4 one-loop ring, solid angle form
p2, q2, w2 = sp.symbols('p2 q2 w2')   # squared momenta at electrons 1, 2, 3
m1, m2, m3 = g1, g2, 0                # lines: g1 (1-3), g2 (2-3), C (1-2)
# vertex between g1 and C carries p, between g2 and C carries q, between g1 and g2 carries w
M = sp.Matrix([[m1**2, (m1**2 + m2**2 + w2) / 2, (m1**2 + m3**2 + p2) / 2],
               [(m1**2 + m2**2 + w2) / 2, m2**2, (m2**2 + m3**2 + q2) / 2],
               [(m1**2 + m3**2 + p2) / 2, (m2**2 + m3**2 + q2) / 2, m3**2]])
detM = sp.expand(M.det())
M12, M13, M23 = M[0, 1], M[0, 2], M[1, 2]
Dlin = sp.expand(m1 * m2 * m3 + M12 * m3 + M13 * m2 + M23 * m1)
Pp, Qp = p2 + g1**2, q2 + g2**2
assert sp.simplify(Dlin - (g2 * Pp + g1 * Qp) / 2) == 0
s4 = sp.factor(sp.expand(Dlin**2 + detM))
assert sp.simplify(s4 - Pp * Qp * (w2 + (g1 + g2)**2) / 4) == 0
out['S4_D'] = str(sp.factor(Dlin))
out['S4_D2_plus_detM'] = str(s4)
assert sp.simplify(4 * detM + target.subs({P: Pp, Q: Qp, dt: w2})) == 0
out['S4_detM_equals_minus_topcut_quadric'] = True
# collinear (Gram boundary) w = p + q in 1-D signed momenta
pp, qq = sp.symbols('p q')
coll = {p2: pp**2, q2: qq**2, w2: (pp + qq)**2}
detM_c = sp.factor(sp.expand(4 * detM.subs(coll)))
out['S4_4detM_collinear'] = str(detM_c)
sq = sp.sqrt(detM_c)
out['S4_sqrt_4detM_collinear'] = str(sp.factor(sp.simplify(sq)))
Dc = sp.expand(Dlin.subs(coll))
I = sp.I
root = pp * qq * (pp + qq) + g2**2 * pp + g1**2 * qq
assert sp.expand(root**2 - detM_c) == 0
plus = sp.factor(sp.expand(2 * Dc + I * root))
minus = sp.factor(sp.expand(2 * Dc - I * root))
out['S4_2D_plus_i_sqrt'] = str(plus)
out['S4_2D_minus_i_sqrt'] = str(minus)

# ---------------------------------------------------------------- S5 strata of the six-line double plane
X = sp.symbols('X')


def genus_of(expr, var):
    """genus of y^2 = expr(var) (expr polynomial): from the number of odd-multiplicity roots."""
    poly = sp.Poly(sp.expand(expr), var)
    fac = sp.factor_list(poly.as_expr(), extension=sp.I)
    odd = 0
    for f, m in fac[1]:
        if m % 2 == 1:
            odd += sp.Poly(f, var).degree()
    deg = poly.degree()
    nb = odd + (1 if (odd % 2 == 1) else 0)   # point at infinity branches if odd count
    return max((nb - 2) // 2, 0), odd, deg


six = (pp**2 + a**2) * (qq**2 + b**2) * ((pp + qq)**2 + X)
num = {a: 1, b: 2, g1: sp.Rational(3, 2), g2: 3, X: sp.Rational(7, 3)}
strata = {
    'p = i g1 (geminal-1 pole)': six.subs(pp, I * g1),
    'q = i g2 (geminal-2 pole)': six.subs(qq, I * g2),
    'p+q = i(g1+g2) (geminal pair)': six.subs(qq, I * (g1 + g2) - pp),
    'p = -i g1 with q = i g2 (point: pseudo-threshold w = i(g2-g1))': six.subs({pp: -I * g1, qq: I * g2}),
}
s5 = {}
for name, e in strata.items():
    var = qq if name.startswith('p = i g1') else pp
    if name.startswith('p = -i g1 with'):
        s5[name] = str(sp.factor(e.subs({a: 1, b: 2, g1: sp.Rational(3, 2), g2: 3}).subs(X, x)))
        continue
    gen, odd, deg = genus_of(e.subs(num), var)
    s5[name] = {'genus': gen, 'odd_roots': odd, 'degree': deg}
out['S5_strata_genus'] = s5

# ---------------------------------------------------------------- S6 exactness of the apex integral on the cut
cc, beta = sp.symbols('c beta', positive=True)
F = -2 * sp.sqrt(cc - beta * T) / (cc * sp.sqrt(T))
dF = sp.simplify(sp.diff(F, T) - T**sp.Rational(-3, 2) * (cc - beta * T)**sp.Rational(-1, 2))
out['S6_dF_minus_integrand'] = str(dF)
assert dF == 0
# the other basis form T^(-1/2)(c - beta T)^(-1/2) dT: its closed period is c-independent
Tm = sp.symbols('Tm', positive=True)
per = sp.integrate(T**sp.Rational(-1, 2) * (cc - beta * T)**sp.Rational(-1, 2), (T, 0, cc / beta))
out['S6_closed_period_of_T^-1/2'] = str(sp.simplify(per))

with open('debug/data/geo_tet_symbolic.json', 'w') as fh:
    json.dump(out, fh, indent=1)
for k, v in out.items():
    print(k, ':', v)
