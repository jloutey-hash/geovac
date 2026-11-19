"""Independent review recomputation (not the memo's model code).

Tree B1: root h_2 = d_z^2 - (d_x^2+d_y^2)/2, middle p_z, leaf p_z (chain). Hat normalisation:
rho = P(d) K0(r), channel-L potential solves (th-L)(th+L+1) Phi = -r^2 W_L, th = r d/dr.

Hand derivation (reviewer):
  leaf: d_z K0 = -K1 P_1 = r^-1 D P_1, D = th K0.      phi1: (th-1)(th+2) phi1 = -x D
  middle: W = r^-1 D phi1 P_1^2 -> channel 2: (2/3) r^-1 D phi1.
  Phi2: (th-2)(th+3) Phi2 = -(2/3) x D phi1.   Let u = (th+3) Phi2, so (th-2) u = -(2/3) x D phi1.
  root (relocated): <h_2(d)(f P_2)> = (1/5) r^-2 (th+1)(th+3) f   [ladder: d_z(fP_2)|_1 = (2/5) r^-1 (th+3) f,
                     d_z(g P_1)|_0 = (1/3) r^-1 (th+2) g; lap part averages to 0]
  G = x^2 <h_2 Psi> = (1/5)(th+1) u = (1/5)(3u - (2/3) x D phi1).

Module coordinates over Q(x): [u, K0 phi1, D phi1, K0 thphi1, D thphi1, K0^2, K0 D, D^2].
Find the least theta-order a (coefficients in Q(x)) with a G having zero non-pure coordinates; clear denominators;
read the top-weight symbol and phi(a)'s leading coefficient at t; compare with the memo's -4t(t-1)^2(t+1)^2.
"""
import json
import sympy as sp

x, t = sp.symbols('x t')
N = 8  # coordinates


def theta_vec(v):
    """theta acting on a coordinate vector with coefficients in Q(x)."""
    out = [sp.Integer(0)] * N
    # product rule part: theta on the coefficient
    for i in range(N):
        out[i] += x * sp.diff(v[i], x)
    c = v
    # u: th u = 2u - (2/3) x D phi1
    out[0] += 2 * c[0]
    out[2] += -sp.Rational(2, 3) * x * c[0]
    # K0 phi1 -> D phi1 + K0 thphi1
    out[2] += c[1]
    out[3] += c[1]
    # D phi1 -> x^2 K0 phi1 + D thphi1
    out[1] += x**2 * c[2]
    out[4] += c[2]
    # K0 thphi1 -> D thphi1 + K0(-thphi1 + 2 phi1 - x D)
    out[4] += c[3]
    out[3] += -c[3]
    out[1] += 2 * c[3]
    out[6] += -x * c[3]
    # D thphi1 -> x^2 K0 thphi1 + D(-thphi1 + 2 phi1 - x D)
    out[3] += x**2 * c[4]
    out[4] += -c[4]
    out[2] += 2 * c[4]
    out[7] += -x * c[4]
    # pure: K0^2 -> 2 K0 D ; K0 D -> D^2 + x^2 K0^2 ; D^2 -> 2 x^2 K0 D
    out[6] += 2 * c[5]
    out[7] += c[6]
    out[5] += x**2 * c[6]
    out[6] += 2 * x**2 * c[7]
    return [sp.together(sp.expand(e)) for e in out]


G = [sp.Rational(3, 5), 0, -sp.Rational(2, 15) * x, 0, 0, 0, 0, 0]
G = [sp.Integer(e) if e == 0 else e for e in G]
vecs = [G]
res = {}
for order in range(1, 9):
    vecs.append(theta_vec(vecs[-1]))
    M = sp.Matrix([[vecs[j][i] for j in range(order + 1)] for i in range(5)])  # non-pure rows
    ns = M.nullspace()
    if ns:
        c = ns[0]
        den = sp.lcm([sp.fraction(sp.together(e))[1] for e in c])
        c = [sp.factor(sp.cancel(e * den)) for e in c]
        g = sp.gcd_list([sp.expand(e) for e in c if e != 0])
        c = [sp.cancel(e / g) for e in c]
        res['theta_order'] = order
        print('least theta order', order)
        # a = sum c_j(x) theta^j ; check evenness in x
        even = all(sp.simplify(e.subs(x, -x) - e) == 0 for e in c)
        odd = all(sp.simplify(e.subs(x, -x) + e) == 0 for e in c)
        print('coefficients even in x:', even, ' odd:', odd)
        res['even'] = bool(even)
        # weights: z = x^2 weight 2, theta weight 1
        polys = [sp.Poly(sp.expand(e), x) for e in c]
        w = max(p.degree() + j for j, p in enumerate(polys) if not p.is_zero)
        print('weight', w, ' z-degree', max(p.degree() for p in polys) // 2)
        res['weight'] = int(w)
        res['z_degree'] = int(max(p.degree() for p in polys) // 2)
        # phi: z^m th^i -> t^-2m (-1)^i th_t^(2m+i) at top order
        lead = sum(p.coeff_monomial(x**(w - j)) * t**(-(w - j)) * (-1)**j
                   for j, p in enumerate(polys) if w - j >= 0)
        lead = sp.factor(sp.simplify(lead * t**w))
        print('t^w * leading coeff of phi(a):', lead)
        res['symbol_times_tw'] = str(lead)
        # pure part of aG
        aG = [sp.Integer(0)] * N
        for j, e in enumerate(c):
            for i in range(N):
                aG[i] += e * vecs[j][i]
        aG = [sp.factor(sp.cancel(e)) for e in aG]
        print('aG non-pure:', aG[:5])
        print('aG pure:', aG[5:])
        res['aG_pure'] = [str(e) for e in aG[5:]]
        res['aG_nonpure_zero'] = all(e == 0 for e in aG[:5])
        res['a_coeffs'] = [str(e) for e in c]
        break

json.dump(res, open('debug/data/tree_attain_l_review_symbol.json', 'w'), indent=1)
