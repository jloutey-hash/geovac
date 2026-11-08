"""HG0, HG1, HG4 of docs/predictions/2026-10-27_h2_geminal_elliptic.md (symbolic part).

s1, s2 are transcribed from Lesiuk-Moszynski PRA 86, 052513 (2012), Eqs. 32-33; HG2 rederives Delta separately.
"""
import sympy as sp
import mpmath as mp

u2, u3, w2, w3, w1, r, D, lam, rho, t1, t2 = sp.symbols('u2 u3 w2 w3 w1 r D lam rho t1 t2')

s1 = (-u2**2*w1**2 - u2**2*w2**2 + u2**2*w3**2 - u3**2*w1**2 + u3**2*w2**2 - u3**2*w3**2
      + w1**4 - w1**2*w2**2 - w1**2*w3**2)
s2 = (u2**4*w2**2 + u2**2*u3**2*w1**2 - u2**2*u3**2*w2**2 - u2**2*u3**2*w3**2 - u2**2*w1**2*w2**2
      + u2**2*w2**4 - u2**2*w2**2*w3**2 + u3**4*w3**2 - u3**2*w1**2*w3**2 - u3**2*w2**2*w3**2
      + u3**2*w3**4 + w1**2*w2**2*w3**2)
# alternative printed form of s1 (LM Eq. 152) as a transcription cross-check
s1_alt = w1**4 - w1**2*(u2**2 + u3**2 + w2**2 + w3**2) + (u2**2 - u3**2)*(w3**2 - w2**2)
print("s1 transcription self-consistent (Eq.32 vs Eq.152):", sp.expand(s1 - s1_alt) == 0)


def op_LM(f, a, b, c, x):
    return a*x*sp.diff(f, x, 4) + 2*a*sp.diff(f, x, 3) + b*x*sp.diff(f, x, 2) + b*sp.diff(f, x) + c*x*f


def op_pf(L, rh, x):
    return (rh*x*sp.diff(L, x, 4) + 2*rh*sp.diff(L, x, 3) + (1 - 2*rh)*x*sp.diff(L, x, 2)
            + (1 - 2*rh)*sp.diff(L, x) - (1 - rh)*x*L)


# ---------------- HG0 ----------------
print("\n== HG0 ==")
# (a) w1 = 0 limit
a0 = sp.expand(s1.subs(w1, 0)); b0 = sp.expand(s2.subs(w1, 0))
F = sp.Function('F')
expr = op_LM(F(r), 0, a0, b0, r)
p2 = -b0/a0
target = a0*(r*sp.diff(F(r), r, 2) + sp.diff(F(r), r) - p2*r*F(r))
print("HG0(a) w1=0 reduces to r f'' + f' - p^2 r f with p^2=-s2/s1:", sp.simplify(expr - target) == 0)
# (b) eq:laplace annihilated by eq:pf
mp.mp.dps = 40


def Lap(Dv, rh, k):
    # k-th derivative of L(D) = int_1^inf e^{-Dx} dx / sqrt((x^2-1)(rho x^2+1-rho))
    return mp.quad(lambda x: (-x)**k*mp.e**(-Dv*x)/mp.sqrt((x**2 - 1)*(rh*x**2 + 1 - rh)), [1, 2, mp.inf])


rh = mp.mpf(1)/3
for Dv in (1, 2):
    Ls = [Lap(Dv, rh, k) for k in range(5)]
    res = rh*Dv*Ls[4] + 2*rh*Ls[3] + (1 - 2*rh)*Dv*Ls[2] + (1 - 2*rh)*Ls[1] - (1 - rh)*Dv*Ls[0]
    print(f"HG0(b) rho=1/3 D={Dv}: residual {mp.nstr(abs(res), 3)} (scale {mp.nstr(abs(Ls[0]), 3)})")
# (c) special cases
c1 = s2.subs({u2: 1, u3: 1, w2: 1, w3: 1, w1: sp.Rational(3, 2)})
c2 = s2.subs({u3: 1, w2: 1, u2: 2, w3: 2, w1: 1})
print("HG0(c) s2 at JC point:", c1, "; at symmetric-exchange point:", c2)

# ---------------- HG1 ----------------
print("\n== HG1 ==")
# Generic identity: a r f'''' + 2a f''' + b r f'' + b f' + c r f with a s^4 + b s^2 + c = a (s^2 - T1)(s^2 - T2).
a, T1, T2 = sp.symbols('a T1 T2', nonzero=True)
b = -a*(T1 + T2); c = a*T1*T2
Lf = sp.Function('L')
t1s = sp.sqrt(T1)
f_of_r = Lf(t1s*r)                      # f(r) = L(D), D = t1 r
lhs = op_LM(f_of_r, a, b, c, r)
rh_map = T1/(T1 - T2)
# express in D
Dsym = sp.Symbol('Dsym', positive=True)
lhs_D = lhs.subs(r, Dsym/t1s).doit()
lhs_D = sp.simplify(lhs_D.replace(lambda e: isinstance(e, sp.Subs), lambda e: e.doit()))
rhsD = op_pf(Lf(Dsym), rh_map, Dsym)
# compare coefficient by coefficient in the derivatives of L at Dsym
derivs = [Lf(Dsym)] + [sp.Derivative(Lf(Dsym), (Dsym, k)) for k in range(1, 5)]
dum = sp.symbols('d0:5')


def coeffs(e):
    e = sp.expand(e.doit())
    for k in range(4, -1, -1):
        e = e.subs(derivs[k], dum[k])
    return [sp.simplify(sp.expand(e).coeff(dd)) for dd in dum]


cl = coeffs(lhs_D); cr = coeffs(rhsD)
ratio = sp.simplify(cl[4]/cr[4])
ok = all(sp.simplify(cl[k] - ratio*cr[k]) == 0 for k in range(5))
print("HG1 generic identity holds:", ok, "; constant =", ratio)
# check that LM's (s1, s2) fit the form b = -a(T1+T2), c = a T1 T2 with a = w1^2 (Vieta) - trivially true for roots
# of w1^2 t^4 + s1 t^2 + s2; verify numerically at a random rational point
pt = {u2: sp.Rational(1), u3: sp.Rational(2), w2: sp.Rational(1), w3: sp.Rational(2), w1: sp.Rational(1)}
S1 = s1.subs(pt); S2 = s2.subs(pt); A = w1.subs(pt)**2
roots = sp.solve(A*sp.Symbol('T')**2 + S1*sp.Symbol('T') + S2, sp.Symbol('T'))
print("LM Table I point: s1 =", S1, " s2 =", S2, " disc =", S1**2 - 4*A*S2, " roots t^2 =", roots)

# ---------------- HG4 ----------------
print("\n== HG4 ==")
x, y, z = sp.symbols('x y z', positive=True)
disc = sp.expand(s1**2 - 4*w1**2*s2)


def classify(name, sub):
    S2v = sp.factor(sp.expand(s2.subs(sub)))
    Dv = sp.factor(sp.expand(disc.subs(sub)))
    print(f"{name}: s2 = {S2v};  s1^2-4w1^2 s2 = {Dv}")


classify("(a) JC u2=u3=w2=w3=x", {u2: x, u3: x, w2: x, w3: x})
classify("(b) KW u3=w2=x,u2=w3=y", {u3: x, w2: x, u2: y, w3: y})
classify("(c) u3=w3=x,u2=w2=y", {u3: x, w3: x, u2: y, w2: y})
classify("(d1) (AA|AA) u3=w2=2z,u2=w3=0", {u3: 2*z, w2: 2*z, u2: 0, w3: 0})
classify("(d2) (AB|AB) all z", {u2: z, u3: z, w2: z, w3: z})
classify("(d3) (AA|BB) u3=w3=2z,u2=w2=0", {u3: 2*z, w3: 2*z, u2: 0, w2: 0})
classify("(d4) (AA|AB) u3=2z,u2=0,w2=w3=z", {u3: 2*z, u2: 0, w2: z, w3: z})
pt_e = {u2: sp.Rational(1, 3), u3: sp.Rational(5, 2), w2: sp.Rational(7, 5), w3: sp.Rational(2, 7), w1: sp.Rational(3, 4)}
print("(e) generic unequal point: s2 =", s2.subs(pt_e), " disc =", disc.subs(pt_e))

# ---------------- HG5 (reported) ----------------
print("\n== HG5 (reported) ==")
for name, sub in [("(AA|BB)", {u3: 2*z, w3: 2*z, u2: 0, w2: 0}), ("(AA|AB)", {u3: 2*z, u2: 0, w2: z, w3: z})]:
    for wv in (sp.Rational(1, 2), 1, sp.Rational(3, 2)):
        sb = dict(sub); sb = {k: v.subs(z, 1) if hasattr(v, 'subs') else v for k, v in sb.items()}; sb[w1] = wv
        S1 = sp.nsimplify(s1.subs(sb)); S2 = sp.nsimplify(s2.subs(sb)); A = sp.nsimplify(wv**2)
        T = sp.Symbol('T'); rts = sp.solve(A*T**2 + S1*T + S2, T)
        rhos = [sp.nsimplify(sp.simplify(p/(p - q))) for p, q in [(rts[0], rts[1]), (rts[1], rts[0])]] if len(rts) == 2 else 'double root'
        print(f"{name} zeta=1 w1={wv}: s1={S1} s2={S2} t^2 roots={rts} rho={rhos}")
