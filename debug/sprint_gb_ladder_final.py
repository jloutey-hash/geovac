"""Sprint GB final: the Bernoulli ladder, corrected signs + clean verdict."""
import sympy as sp
pi = sp.pi
z = sp.zeta
def line(a, b): print(f"  {a:<52} {b}")

print("="*80)
print("THE LADDER: rung B_{2n}, two windows glued by the functional equation")
print("="*80)
print(f"{'rung':<6}{'M2 window  zeta(2n)':<26}{'skeleton window zeta(1-2n)':<28}")
for n in (1,2,3):
    line(f"B_{2*n}={sp.bernoulli(2*n)}", f"zeta({2*n})={z(2*n)}    zeta({1-2*n})={z(1-2*n)}")

print("\n" + "="*80)
print("GeoVac coefficients placed on the ladder (exact, signs fixed)")
print("="*80)
checks = [
    ("conical tip 1/12 (G4-2)",            sp.Rational(1,12),  -z(-1),       "-zeta(-1) = B_2/2"),
    ("replica entropy deriv 1/6 (G4-4f)",  sp.Rational(1,6),   -2*z(-1),     "-2 zeta(-1) = B_2"),
    ("per-t UV 1/24 (v3.20.0)",            sp.Rational(1,24),  -z(-1)/2,     "-zeta(-1)/2 = B_2/4  [= string c/24]"),
    ("scalar Casimir S^3 1/240 (Pap.35)",  sp.Rational(1,240),  z(-3)/2,     "+zeta(-3)/2  [B_4 rung, skeleton]"),
    ("alpha-conj F = pi^2/6 (Pap.2)",      z(2),                pi**2*sp.bernoulli(2), "pi^2 B_2 = zeta(2) [B_2, M2 window]"),
]
allfit = True
for name, c, cand, desc in checks:
    ok = sp.simplify(c - cand) == 0
    allfit &= ok
    line(f"{name}", f"{'FIT ' if ok else 'MISS'}  {desc}")

print("\n  CONTROLS (should NOT sit on the ladder):")
kappa = sp.Rational(-1,16)
dirac_cas = sp.Rational(17,480)
line("kappa = -1/16 (Fock Jacobian 1/Omega^4)",
     "no clean zeta form" if all(sp.simplify(kappa-x)!=0 for x in [-z(-1),-z(-1)/2,z(-3)/2,-z(-3)/2,-2*z(-1)]) else "FIT?!")
# Dirac Casimir: integer vs half-integer
line("Dirac Casimir 17/480 vs scalar rung zeta(-3)/2",
     f"residual {sp.nsimplify(dirac_cas - z(-3)/2)}  -> NOT on scalar rung")
line("  Hurwitz zeta(-3,1/2) (half-integer shift)", z(-3, sp.Rational(1,2)))
line("  => 17/480 lives in half-integer sector = scalar/spinor = M2/M3 split", "")

print("\n" + "="*80)
print("INTERNAL CONSISTENCY (the non-trivial check, not a free coincidence)")
print("="*80)
a = sp.symbols('alpha', positive=True)
tip = sp.Rational(1,12)*(1/a - a)
d = sp.diff(tip,a).subs(a,1)
line("tip coeff = B_2/2; d/dalpha|_1 (entropy) = ", f"{d} = -B_2 = 2 zeta(-1)")
line("  replica derivative FORCED by tip (not independent):", sp.simplify(d - 2*z(-1))==0)

print("\n" + "="*80)
print(f"VERDICT INPUTS: B_2 rung 3-fold populated & internally forced; "
      f"B_4 rung independently populated (scalar Casimir);")
print(f"  spinor (17/480) correctly off-rung (half-integer); control (kappa) correctly off-rung.")
print(f"  all primary fits exact: {allfit}")
