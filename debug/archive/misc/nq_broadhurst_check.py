"""Compare Broadhurst arXiv:2212.01962 Eqs. 8-9 (anomalous thresholds s+-, Minkowski s) with the record's curve
C(z, eps) = 0 (Euclidean z = -s), mapping his (m1..m5)^2 = (A, B, eps, G1, G2). Credit check, not a scored item."""
import sympy as sp
A, B, G1, G2, e, z = sp.symbols("A B G1 G2 eps z")
lam = lambda a, b, c: a*a + b*b + c*c - 2*a*b - 2*b*c - 2*c*a
u, v = A - G1, B - G2
K = (z + A + B) / 2
C = 4*e*K**2 - 2*(u + e)*(v + e)*K + B*(u + e)**2 + A*(v + e)**2 - 4*e*A*B
alpha = (A - G1)*(B - G2)/e - e
for vals in [dict(A=1, B=sp.Rational(9, 4), G1=sp.Rational(49, 9), G2=sp.Rational(64, 9), eps=sp.Rational(1, 400)),
             dict(A=3, B=sp.Rational(7, 5), G1=sp.Rational(2, 9), G2=11, eps=sp.Rational(5, 7))]:
    sb = {A: vals["A"], B: vals["B"], G1: vals["G1"], G2: vals["G2"], e: vals["eps"]}
    D = sp.sqrt(lam(sb[A], sb[e], sb[G1]) * lam(sb[B], sb[e], sb[G2]))
    mid = (sb[A] + sb[B] - 2*sb[e] + sb[G1] + sb[G2] - alpha.subs(sb)) / 2
    for sgn in (1, -1):
        s = mid + sgn * D / (2*sb[e])
        print(vals, "s%s" % ("+" if sgn > 0 else "-"), sp.N(s, 15), " C(z=-s) =", sp.N(C.subs(sb).subs(z, -s), 15))
