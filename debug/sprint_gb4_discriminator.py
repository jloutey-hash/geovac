"""
Sprint GB-4: the discriminator. Does a NON-Casimir GeoVac observable land
on B_4 / B_6 -- and specifically, is the POSITIVE-EVEN (M2) window of each
rung lit, by different physics than the Casimir (which lit the negative
window)? Plus a real negative control: the combinatorial core must be OFF.
"""
import sympy as sp
x = sp.symbols('x', positive=True)
z = sp.zeta

print("="*74)
print("1. THERMAL window: Stefan-Boltzmann / blackbody carries zeta(d+1)")
print("   (Bose-Einstein integral I_d = int_0^inf x^d/(e^x-1) dx = Gamma(d+1)zeta(d+1))")
print("="*74)
for d, spacetime in [(3,"S^3 x S^1  (4D)"), (5,"S^5 x S^1  (6D)")]:
    I = sp.gamma(d+1)*z(d+1)          # textbook Bose-Einstein integral
    I_closed = sp.nsimplify(sp.simplify(I))
    zval = z(d+1)
    print(f"  {spacetime}:  I_{d} = Gamma({d+1})zeta({d+1}) = {I_closed}")
    print(f"      rung content: zeta({d+1}) = {zval}  (B_{d+1} M2 window)")
# cross-check 4D against Sprint TD Track 1's -pi^2/90 T^4 structure
print(f"\n  4D cross-check: zeta(4) = {z(4)} = pi^4/90  (matches Sprint TD Track 1 SB)")

print("\n" + "="*74)
print("2. THE TWO-WINDOW TABLE: each rung lit from BOTH sides by DIFFERENT physics")
print("="*74)
rungs = [
    ("B_2", z(-1), "S^1 Casimir / conical tip (vacuum)",  z(2),  "S^1 x S^1 thermal"),
    ("B_4", z(-3), "S^3 Casimir 1/240 (vacuum)",          z(4),  "S^3 x S^1 Stefan-Boltzmann (thermal)"),
    ("B_6", z(-5), "S^5 Casimir -31/60480 (vacuum)",      z(6),  "S^5 x S^1 Stefan-Boltzmann (thermal)"),
]
print(f"  {'rung':<5}{'skeleton zeta(1-2k)':<22}{'M2 zeta(2k)':<16}{'two windows, two observables'}")
for name, zneg, csrc, zpos, tsrc in rungs:
    print(f"  {name:<5}{str(zneg)+'  ('+csrc.split('(')[0].strip()+')':<22}")
    print(f"       skeleton: {str(zneg):<10} <- {csrc}")
    print(f"       M2:       {str(zpos):<10} <- {tsrc}")

print("\n" + "="*74)
print("3. The two windows ARE the functional equation = temperature inversion")
print("="*74)
print("  Casimir (vacuum, negative-integer zeta)  <-- T<->1/T -->  thermal (positive-even zeta)")
print("  The Casimir<->blackbody duality on S^d x S^1 IS the zeta functional equation")
print("  (temperature-inversion symmetry). Same Bernoulli rung, both windows, distinct physics.")

print("\n" + "="*74)
print("4. NEGATIVE CONTROL: the combinatorial core must be OFF the ladder")
print("="*74)
core = {"B (alpha Casimir trace) = 42": sp.Integer(42),
        "Delta (alpha) = 1/40": sp.Rational(1,40),
        "kappa (Fock Jacobian) = -1/16": sp.Rational(-1,16)}
rung_vals = [z(-1),z(-3),z(-5),z(2),z(4),z(6),
             z(-1)/2,-z(-1),z(-3)/2,2*z(-1),-z(-1)/2]   # all simple rung normalizations seen
for name, val in core.items():
    on = any(sp.simplify(val - r)==0 for r in rung_vals)
    print(f"  {name:<34} on a rung? {on}")
print("  (Sprint TD Track 5 also: GeoVac correlation entropy S_full(GS) PSLQ-null off the engine.)")
print("  => the ladder is the SPECTRAL-GEOMETRY sector's functional-equation structure,")
print("     NOT a universal 'everything in GeoVac carries a zeta'.")
