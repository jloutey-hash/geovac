"""Geometric route to Y off the slice: the cut polygons and the pieces they carry.

No operator for Y or its kite is computed.  Inputs are the six lines of the tetrahedron
(position-space vertices O, 1, 2, 3):
    O-1 density a=1 (index 3/2, Fock factor), O-2 density b=2, O-3 apex density t,
    1-3 geminal g1=3/2 (index 2, rational), 2-3 geminal g2=3, 1-2 Coulomb (index 1, massless).

1. Cut polygons.  A singular value of t is a closed, flat polygon whose sides are the
   apex line and the other lines of a cut (cocycle) of the tetrahedron, sides of length
   = mass (dual-tetrahedron picture: momentum conservation at a vertex closes a triangle
   of masses; a cut closes a polygon).  Enumerate all cuts that contain O-3, all sign
   patterns, and the 'adjacent' triangles that are not cuts.
2. Threshold exponent per cut: beta = ((k-1)(D-1)-2)/2 + sum_i (1-nu_i), D=3;
   integer -> unipotent (log), half-integer -> eigenvalue -1 (Picard-Lefschetz with
   the square-root twist of the Fock lines).
3. Fock count per cut -> rung of the piece the cut carries (density-topology rule of
   Paper 59 extended): one Fock factor -> genus 0, two -> elliptic curve, three -> K3.
4. The two elliptic pieces: the BBBG quartic y^2=(k^2+m^2)((k-i g)^2+t^2) (two Fock
   factors shifted by the geminal pole) -- its j-invariant, equal to Legendre's with
   z=(t^2+m^2-g^2)/(2mt); singular fibres over x=t^2.
5. The three-Fock piece for general a,b: discriminant of Bailey's moduli.
6. Slice reduction a=b=1, g1=g2=g: which pieces coincide.
"""
import json
import itertools
import sympy as sp

R = sp.Rational
a, b, t = R(1), R(2), sp.Symbol('t', positive=True)
g1, g2 = R(3, 2), R(3)
x = sp.Symbol('x')

lines = {'O1': ('dens', a, R(3, 2)), 'O2': ('dens', b, R(3, 2)), 'O3': ('dens', 't', R(3, 2)),
         '13': ('gem', g1, 2), '23': ('gem', g2, 2), '12': ('coul', 0, 1)}
verts = ['O', '1', '2', '3']
out = {}

# ---------------------------------------------------------------- 1. cuts containing O-3
cuts = []
for r in range(1, 4):
    for side in itertools.combinations(verts, r):
        S = set(side)
        if 'O' not in S:
            continue
        crossing = [e for e in lines if (e[0] in S) != (e[1] in S)]
        if 'O3' in crossing:
            cuts.append((tuple(sorted(S)), crossing))
cut_table = []
sing_x = set()
for S, crossing in cuts:
    others = [e for e in crossing if e != 'O3']
    masses = [lines[e][1] for e in others]
    vals = set()
    for signs in itertools.product([1, -1], repeat=len(masses)):
        v = abs(sum(s * m for s, m in zip(signs, masses)))
        if v != 0:
            vals.add(v)
    k = len(crossing)
    nus = [lines[e][2] for e in crossing]
    beta = R((k - 1) * 2 - 2, 2) + sum(1 - nu for nu in nus)
    nfock = sum(1 for e in crossing if lines[e][0] == 'dens')
    cut_table.append({'side_with_O': S, 'lines': crossing, 't_values': sorted(vals),
                      'x_values': sorted(v**2 for v in vals), 'beta': str(beta),
                      'beta_type': 'log (unipotent)' if beta.q == 1 else 'eigenvalue -1 (sqrt)',
                      'fock_factors': nfock,
                      'rung': {1: 'genus 0 (algebraic)', 2: 'elliptic curve', 3: 'K3 surface'}[nfock]})
    sing_x |= {v**2 for v in vals}
out['cuts'] = [{k_: (str(v_) if not isinstance(v_, (list, tuple)) else [str(w) for w in v_])
                for k_, v_ in c.items()} for c in cut_table]
out['singular_x_predicted'] = [str(v) for v in sorted(sing_x)]
given = [R(1, 4), R(1), R(9, 4), R(4), R(9), R(49, 4), R(16), R(81, 4)]
out['matches_given_eight'] = sorted(sing_x) == sorted(given)
# adjacent (non-cut) triangles: apex with a density and the geminal sharing its vertex
adj = {'(t,a,g1)': [a + g1, abs(a - g1)], '(t,b,g2)': [b + g2, abs(b - g2)]}
out['adjacent_not_cuts_x'] = {k_: [str(v**2) for v in vs] for k_, vs in adj.items()}

# ---------------------------------------------------------------- 4. elliptic pieces


def quartic_j(coeffs):
    """j-invariant of y^2 = quartic(k) via the invariants I, J of a binary quartic."""
    a4, a3, a2, a1, a0 = coeffs
    I = 12 * a4 * a0 - 3 * a3 * a1 + a2**2
    J = 72 * a4 * a2 * a0 + 9 * a3 * a2 * a1 - 27 * a4 * a1**2 - 27 * a0 * a3**2 - 2 * a2**3
    return sp.simplify(1728 * 4 * I**3 / (4 * I**3 - J**2))


k = sp.Symbol('k')
tt = sp.Symbol('tt', positive=True)
ell = {}
for name, m, g in [('E_A (cut O2|13: density a=1, geminal g2=3)', a, g2),
                   ('E_B (cut O1|23: density b=2, geminal g1=3/2)', b, g1)]:
    Qk = sp.expand((k**2 + m**2) * ((k - sp.I * g)**2 + tt**2))
    co = sp.Poly(Qk, k).all_coeffs()
    jq = sp.simplify(quartic_j(co))
    z = (tt**2 + m**2 - g**2) / (2 * m * tt)
    lam = (1 - z) / 2
    jl = sp.simplify(256 * (lam**2 - lam + 1)**3 / (lam**2 * (lam - 1)**2))
    same = sp.simplify(jq - jl) == 0
    jx = sp.factor(sp.simplify(jl.subs(tt, sp.sqrt(x))))
    num, den = sp.fraction(jx)
    poles = sp.roots(sp.Poly(den, x))
    deg_inf = sp.degree(num, x) - sp.degree(den, x)
    ell[name] = {'j_quartic_equals_j_legendre': bool(same), 'j(x)': str(jx),
                 'poles_of_j_in_x_with_order': {str(r_): int(mu) for r_, mu in poles.items()},
                 'order_of_pole_at_infinity': int(deg_inf),
                 'lambda(t)': str(sp.factor(sp.simplify(lam)))}
out['elliptic_pieces'] = ell

# ---------------------------------------------------------------- 5. three-Fock piece, general a,b
A, B = sp.symbols('A B', positive=True)
P = (x + A**2 - B**2) / x  # k+^2 + k-^2 in the regime t largest (c=t)
Q2 = A**2 / x
out['three_fock_discriminant_general'] = str(sp.factor(sp.together(P**2 - 4 * Q2)))
out['three_fock_lambda_sum_general'] = str(sp.simplify(P))

# ---------------------------------------------------------------- 6. slice reduction
gg = sp.Symbol('g', positive=True)
slice_vals = {'density pair (a+b, |a-b|)': [2, 0], 'geminal pair (g1+g2, |g1-g2|)': [2 * gg, 0],
              'a with g2': [1 + gg, sp.Abs(1 - gg)], 'b with g1': [1 + gg, sp.Abs(1 - gg)]}
out['slice_singular_t'] = {k_: [str(v) for v in vs] for k_, vs in slice_vals.items()}
out['slice_three_fock_lambda_sum'] = str(sp.simplify(P.subs({A: 1, B: 1})))  # =1: E- = E+ twisted
out['slice_note'] = ('On the slice lambda_+ + lambda_- = 1, so E_- is E_+ up to a twist and '
                     'E_+ (x) E_- = Sym^2 + wedge^2 (ranks 3 + 1); E_A and E_B coincide.')

with open('debug/data/geo_y_pieces.json', 'w') as f:
    json.dump(out, f, indent=1, default=str)
for k_, v in out.items():
    if k_ == 'cuts':
        for c in v:
            print(c)
    else:
        print(k_, ':', v)
