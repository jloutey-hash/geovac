"""Phase-1: numbers implied by the hand model registered in docs/predictions/2026-10-27_kite_ninth_point_origin.md.
Model (by hand, not yet checked by computer): all five lines on shell, momenta coplanar with w, Coulomb momentum q
with q^2 = -eps. With K = k.(w-k), c1 = k.q = -(u+eps)/2, c2 = (w-k).q = (v+eps)/2, u = A-G1, v = B-G2:
   C(K, eps) = 4 eps K^2 - 2 (u+eps)(v+eps) K + B (u+eps)^2 + A (v+eps)^2 - 4 eps A B = 0,   z = 2K - A - B.
This script only evaluates the model. It runs no reduction, no Groebner basis, no critical-point count."""
import sympy as sp

A, B, G1, G2, K, e, z = sp.symbols("A B G1 G2 K eps z")
u, v = A - G1, B - G2
C = 4 * e * K**2 - 2 * (u + e) * (v + e) * K + B * (u + e)**2 + A * (v + e)**2 - 4 * e * A * B
# eps = 0: linear in K
K0 = sp.solve(C.subs(e, 0), K)[0]
z0 = sp.factor(2 * K0 - A - B)
print("model z9 =", z0)
z9_rec = (A * G2 - B * G1) * (A - B - G1 + G2) / ((A - G1) * (B - G2))
print("equals the measured formula:", sp.simplify(z0 - z9_rec) == 0)
# curve in (z, eps): substitute K = (z + A + B)/2
Cz = sp.expand(C.subs(K, (z + A + B) / 2))
print("curve C(z, eps) =", sp.collect(sp.factor_terms(Cz), [e, z]))
# slope dz/deps at (z9, 0)
dzde = sp.factor(-sp.diff(Cz, e) / sp.diff(Cz, z)).subs(e, 0)
P2 = {A: 1, B: sp.Rational(9, 4), G1: sp.Rational(49, 9), G2: sp.Rational(64, 9)}
z9 = z0.subs(P2)
kappa = sp.nsimplify(sp.simplify(dzde.subs(P2).subs(z, z9)))
print("P2: z9 =", z9, " kappa = dz/deps =", kappa, "=", sp.N(kappa, 12))
Cz2 = sp.expand(Cz.subs(P2))
print("P2 curve:", sp.factor(Cz2))
for ev in [sp.Rational(1, 400), sp.Rational(1, 1600), sp.Rational(1, 6400)]:
    roots = sp.Poly(Cz2.subs(e, ev), z).all_roots()
    print("eps =", ev, "roots z:", [sp.N(r, 15) for r in roots],
          " near-branch (z - z9)/eps:", [sp.N((r - z9) / ev, 12) for r in roots if abs(sp.N(r - z9)) < 0.05])
# rational points on the conic through (z9, 0): line z = z9 + sigma*eps
sg = sp.Symbol("sigma")
lin = sp.factor(Cz2.subs(z, z9 + sg * e))
print("on z = z9 + sigma eps:", lin)
sol = [s for s in sp.solve(sp.cancel(lin / e), e)]
print("eps(sigma) =", sol)
# threshold controls: z = -(c +- m)^2, c in {a+g2, a-g2} = {11/3, -5/3}
for m in [sp.Rational(1, 20), sp.Rational(1, 40), sp.Rational(1, 80)]:
    print("m =", m, " control z:", [-(c + sg2 * m)**2 for c in (sp.Rational(11, 3), sp.Rational(-5, 3)) for sg2 in (1, -1)])
