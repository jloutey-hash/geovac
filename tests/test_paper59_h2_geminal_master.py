"""Guard for Paper 59 sec:f12, the two-centre Slater-geminal master integral (record
docs/predictions/2026-10-27_h2_geminal_elliptic.md, HG1, HG2, HG4).

The master f(r) = r int int e^{-u3 r1A - u2 r1B - w2 r2A - w3 r2B - w1 r12} / (16 pi^2 r1A r1B r2A r2B r12)
obeys (Lesiuk and Moszynski, PRA 86, 052513 (2012), Eq. 54)
    w1^2 r f'''' + 2 w1^2 f''' + s1 r f'' + s1 f' + s2 r f = 0,
with s1, s2 transcribed from their Eqs. 32-33.  Claims guarded:
  * HG2: Delta(u1) = w1^2 u1^4 + s1 u1^2 + s2 is 4 x the D = 3 Gram determinant of the tetrahedron's loop
    momenta with every line on shell.  This rederives s1 and s2 independently of the transcription.
  * HG1: under r = D / sqrt(T1), rho = T1 / (T1 - T2), with T1, T2 the roots of w1^2 T^2 + s1 T + s2, the operator
    is a constant multiple of Paper 59 eq:pf.
  * The Bessel limit: at w1 = 0 the operator is second order, r f'' + f' - p^2 r f, p^2 = -s2/s1.
  * HG4: s2 = 0 (so rho = 1, the equal-scale collapse) on the James-Coolidge, Kolos-Wolniewicz, (AA|AA) and (AB|AB)
    patterns; s2 != 0 (genus one) on (AA|BB), (AA|AB) and a generic unequal point.
"""
import sympy as sp

u1, u2, u3, w1, w2, w3, r, x, y, z = sp.symbols('u1 u2 u3 w1 w2 w3 r x y z')

S1 = (-u2**2*w1**2 - u2**2*w2**2 + u2**2*w3**2 - u3**2*w1**2 + u3**2*w2**2 - u3**2*w3**2
      + w1**4 - w1**2*w2**2 - w1**2*w3**2)
S2 = (u2**4*w2**2 + u2**2*u3**2*w1**2 - u2**2*u3**2*w2**2 - u2**2*u3**2*w3**2 - u2**2*w1**2*w2**2
      + u2**2*w2**4 - u2**2*w2**2*w3**2 + u3**4*w3**2 - u3**2*w1**2*w3**2 - u3**2*w2**2*w3**2
      + u3**2*w3**4 + w1**2*w2**2*w3**2)


def gram_on_shell():
    """Gram determinant of k1, k2, k3 with k_i^2 = -u_i^2 and (k_i - k_j)^2 = -m_ij^2 (LM Eq. 13 labels)."""
    sq = {1: -u1**2, 2: -u2**2, 3: -u3**2}
    m = {frozenset((3, 2)): w1, frozenset((1, 3)): w2, frozenset((2, 1)): w3}

    def kk(i, j):
        if i == j:
            return sq[i]
        return (sq[i] + sq[j] + m[frozenset((i, j))]**2) / 2
    return sp.expand(sp.Matrix(3, 3, lambda i, j: kk(i + 1, j + 1)).det())


def lm_operator(f, a, b, c):
    return a*r*sp.diff(f, r, 4) + 2*a*sp.diff(f, r, 3) + b*r*sp.diff(f, r, 2) + b*sp.diff(f, r) + c*r*f


def test_delta_is_four_times_the_on_shell_gram_determinant():
    delta = sp.expand(w1**2*u1**4 + S1*u1**2 + S2)
    assert sp.expand(4*gram_on_shell() - delta) == 0


def test_lm_operator_is_eq_pf_rescaled():
    a, T1, T2, D = sp.symbols('a T1 T2 D', positive=True)
    L = sp.Function('L')
    rho = T1 / (T1 - T2)
    f = L(sp.sqrt(T1)*r)
    lhs = lm_operator(f, a, -a*(T1 + T2), a*T1*T2).subs(r, D/sp.sqrt(T1)).doit()
    pf = (rho*D*sp.diff(L(D), D, 4) + 2*rho*sp.diff(L(D), D, 3) + (1 - 2*rho)*D*sp.diff(L(D), D, 2)
          + (1 - 2*rho)*sp.diff(L(D), D) - (1 - rho)*D*L(D))
    d = sp.symbols('d0:5')

    def coeffs(e):
        e = sp.expand(e.doit())
        for k in range(4, 0, -1):
            e = e.subs(sp.Derivative(L(D), (D, k)), d[k])
        e = sp.expand(e.subs(L(D), d[0]))
        assert not e.has(L), "terms not in the derivatives of L(D): wrong rescaling"
        return [sp.simplify(e.coeff(dk)) for dk in d]
    cl, cr = coeffs(lhs), coeffs(pf)
    const = sp.simplify(cl[4] / cr[4])
    assert const != 0
    assert const.free_symbols <= {a, T1, T2}
    assert all(sp.simplify(cl[k] - const*cr[k]) == 0 for k in range(5))
    # the (s1, s2) of LM are the Vieta coefficients of the roots T1, T2
    T = sp.Symbol('T')
    roots = sp.solve(w1**2*T**2 + S1*T + S2, T)
    assert sp.simplify(sp.expand(-w1**2*(roots[0] + roots[1]) - S1)) == 0
    assert sp.simplify(sp.expand(w1**2*roots[0]*roots[1] - S2)) == 0


def test_w1_zero_is_the_bessel_equation():
    F = sp.Function('F')
    a0, b0 = sp.expand(S1.subs(w1, 0)), sp.expand(S2.subs(w1, 0))
    expr = lm_operator(F(r), 0, a0, b0)
    target = a0*(r*sp.diff(F(r), r, 2) + sp.diff(F(r), r) + (b0/a0)*r*F(r))
    assert sp.simplify(expr - target) == 0
    assert a0 != 0


def _s2(sub):
    return sp.factor(sp.expand(S2.subs(sub)))


def test_symmetric_patterns_collapse_to_rho_one():
    for sub in ({u2: x, u3: x, w2: x, w3: x},           # James-Coolidge
                {u3: x, w2: x, u2: y, w3: y},           # Kolos-Wolniewicz / symmetric exchange
                {u3: 2*z, w2: 2*z, u2: 0, w3: 0},       # (AA|AA)
                {u2: z, u3: z, w2: z, w3: z}):          # (AB|AB)
        assert _s2(sub) == 0, sub
    # s2 = 0 puts a root at T = 0, so rho = T1/(T1 - 0) = 1
    T = sp.Symbol('T')
    roots = sp.solve(w1**2*T**2 + S1.subs({u2: x, u3: x, w2: x, w3: x})*T, T)
    nonzero = [t for t in roots if sp.simplify(t) != 0]
    assert len(nonzero) == 1 and sp.simplify(nonzero[0] / (nonzero[0] - 0)) == 1


def test_lopsided_patterns_are_genus_one():
    assert sp.expand(_s2({u3: x, w3: x, u2: y, w2: y}) - (x**2 - y**2)**2*(2*x**2 + 2*y**2 - w1**2)) == 0
    assert sp.expand(_s2({u3: 2*z, w3: 2*z, u2: 0, w2: 0}) - 16*z**4*(8*z**2 - w1**2)) == 0
    assert sp.expand(_s2({u3: 2*z, u2: 0, w2: z, w3: z}) - z**4*(16*z**2 - 3*w1**2)) == 0
    generic = {u2: sp.Rational(1, 3), u3: sp.Rational(5, 2), w2: sp.Rational(7, 5), w3: sp.Rational(2, 7),
               w1: sp.Rational(3, 4)}
    assert S2.subs(generic) != 0
    assert (S1**2 - 4*w1**2*S2).subs(generic) != 0
