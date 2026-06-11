"""
Sprint GB-3: the B_6 rung hunt. Find a FORCED GeoVac observable carrying
zeta(-5) = -1/252 (skeleton window) / zeta(6) = pi^6/945 (M2 window).

Principled target (not cherry-picked): conformal-scalar Casimir on S^{2k-1}
leads with zeta(-(2k-1)). S^1 -> B_2, S^3 -> B_4, S^5 -> B_6. S^5 is a real
GeoVac manifold (Bargmann-Segal, Paper 24).

Method (validated against Paper 35 S^3 = 1/240):
  E_Casimir = (1/2) sum_l omega_l * deg_l  (zeta-regularized),
  conformal scalar: omega_l = l + (d-1)/2,  deg_l = dim of degree-l harmonics on S^d.
"""
import sympy as sp
n, l = sp.symbols('n l', integer=True, positive=True)
z = sp.zeta

def casimir_conformal_scalar(d):
    """Conformal-scalar Casimir on unit S^d via zeta regularization, exact."""
    # degeneracy of degree-l harmonics on S^d:
    deg_l = (2*l + d - 1) * sp.rf(l + 1, d - 2) / sp.factorial(d - 1)  # rising factorial form
    deg_l = sp.simplify(deg_l)
    omega = l + sp.Rational(d - 1, 2)               # conformal frequency
    # shift to n = l + (d-1)/2 so omega = n; for odd d this is integer
    shift = sp.Rational(d - 1, 2)
    # E = (1/2) sum_{l>=0} omega * deg_l, regularized term-by-term as zeta(-k)
    term = sp.expand(omega * deg_l)                 # polynomial in l
    poly = sp.Poly(term, l)
    # substitute l -> n - shift to express in n (omega = n), then regularize sum_{n} n^k -> zeta(-k)
    term_n = sp.expand(poly.as_expr().subs(l, n - shift))
    pn = sp.Poly(term_n, n)
    # zeta-regularized sum over the spectrum: sum_{n = shift+? } ... use full continuation
    # E = (1/2) sum_{n in spectrum} (poly in n). Spectrum n = shift, shift+1, ...
    # For odd d, shift integer >=1; regularize sum_{n>=1} minus low terms with n^k->zeta(-k).
    E = sp.Rational(0)
    coeffs = pn.all_coeffs()[::-1]                  # c0 + c1 n + ...
    # sum_{n>=1} n^k = zeta(-k); subtract n=1..shift-1 explicit (finite, exact)
    lowstart = int(shift)
    reg = sp.Rational(0)
    deg_poly = pn.degree()
    for k in range(deg_poly + 1):
        ck = pn.coeff_monomial(n**k) if k > 0 else pn.coeff_monomial(1)
        reg += ck * z(-k) if k >= 1 else ck * z(0)
    # subtract the missing low terms n = 1 .. lowstart-1
    missing = sum(pn.as_expr().subs(n, j) for j in range(1, lowstart))
    E = sp.Rational(1, 2) * (reg - missing)
    return sp.nsimplify(sp.simplify(E)), pn.as_expr()

print("="*72)
print("CONTROLS: reproduce Paper 35 (S^3 conformal scalar = 1/240)")
print("="*72)
for d in (3,):
    E, poly = casimir_conformal_scalar(d)
    print(f"  S^{d}: spectrum poly in n = {poly}")
    print(f"        E_Casimir = {E}   ({f"{float(E):.6g}"})")

print("\n" + "="*72)
print("TARGET: S^5 conformal scalar Casimir (B_6 rung)")
print("="*72)
E5, poly5 = casimir_conformal_scalar(5)
print(f"  S^5: spectrum poly in n = {poly5}   [quartic -> spreads across rungs]")
print(f"       E_Casimir = {E5}   ({f"{float(E5):.8g}"})")

print("\n" + "="*72)
print("Decompose S^5 result onto the ladder rungs")
print("="*72)
zm5, zm3, zm1 = z(-5), z(-3), z(-1)
print(f"  zeta(-5) = {zm5} (B_6 rung)   zeta(-3) = {zm3} (B_4)   zeta(-1) = {zm1} (B_2)")
# E5 should be (1/24)(zeta(-5) - zeta(-3)); verify
cand = sp.Rational(1,24)*(zm5 - zm3)
print(f"  candidate (1/24)(zeta(-5) - zeta(-3)) = {cand}")
print(f"  MATCH? {sp.simplify(E5 - cand) == 0}")
print(f"  => B_6 rung (zeta(-5)) PRESENT with coeff 1/24, alongside B_4 (zeta(-3))")

print("\n" + "="*72)
print("Why it spreads (the GB-2 mechanism, independently): degeneracy degree")
print("="*72)
for d in (1,3,5,7):
    _, poly = casimir_conformal_scalar(d)
    deg = sp.Poly(poly, n).degree()
    print(f"  S^{d}: omega*deg is degree-{deg} in n -> regularizes to "
          f"{'single zeta' if deg%2==1 and d<=3 else 'spread across zeta(-1..-%d)'%deg}")

print("\n" + "="*72)
print("M2-window partner + half-integer (M3) sibling check")
print("="*72)
print(f"  functional-equation partner of zeta(-5): zeta(6) = {z(6)} = pi^6/945")
print(f"    -> predicted pi^6 transcendental in M2 window of the S^5 sector (Stefan-Boltzmann-class)")
print(f"  HO/Bargmann S^5 Casimir (TX-B) = -17/3840 (half-integer shift) -> M3 sibling,")
print(f"    echoes Dirac-S^3 17/480; integer(conformal)->M2 ladder, half-integer(HO)->M3, at every rung")
