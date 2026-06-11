"""GB-5 (grounding): is the multi-electron hyperspherical sphere the SAME
sphere that carries the gravity ladder rungs? Check the dimension sequences
and the shared harmonic-degeneracy polynomial."""
import sympy as sp
l, N, k = sp.symbols('l N k', integer=True, nonnegative=True)

print("="*70)
print("1. The TWO dimensional ladders")
print("="*70)
print("  Multi-electron hyperspherical (Avery):  N electrons -> S^(3N-1)")
for Ne, atom in [(1,"H "),(2,"He"),(3,"Li"),(4,"Be")]:
    print(f"     {atom} ({Ne}e): S^{3*Ne-1}")
print("  Gravity conformal-Casimir ladder:        rung B_2k -> S^(2k-1)")
for kk in (1,2,3,4):
    print(f"     B_{2*kk}: S^{2*kk-1}")
print("  --> they COINCIDE only at S^5 (He, 2 electrons  <->  B_6 rung).")
print("      Li (S^8) and Be (S^11) are EVEN/odd off the gravity odd-sphere ladder.")

print("\n" + "="*70)
print("2. The SHARED object: degree-l scalar harmonic dimension on S^5")
print("="*70)
# dim of degree-l harmonics on S^d:
def harm_dim(d, lv):
    return sp.simplify((2*lv + d - 1) * sp.rf(lv + 1, d - 2) / sp.factorial(d - 1))
d5 = sp.simplify(harm_dim(5, l))
print(f"  He hyperspherical (SO(6) on S^5) scalar degeneracy:  {sp.factor(d5)}")
print(f"  GB-3 gravity Casimir used the SAME polynomial:        (l+1)(l+2)^2(l+3)/12")
print(f"  match? {sp.simplify(d5 - (l+1)*(l+2)**2*(l+3)/12) == 0}")
print(f"  degree in l: {sp.Poly(d5, l).degree()}  (quartic -> this is the multiplicity")
print(f"  that (a) counts He correlation channels AND (b) spread the Casimir across")
print(f"  zeta(-5),zeta(-3) and obstructed the functional equation in GB-2/GB-4)")

print("\n" + "="*70)
print("3. The founding motif (Fock 1935) the two share")
print("="*70)
print("  Fock:   hydrogen momentum space -> S^3, makes hidden SO(4) manifest")
print("  Avery:  N-electron config space -> S^(3N-1), makes correlation calculable")
print("  Bargmann: 3D HO phase space    -> S^5 (Hardy), makes spectrum pi-free")
print("  Gravity: round S^(2k-1) Casimir -> Bernoulli rungs")
print("  Common move: ELEVATE DIMENSION so a hard structure becomes a sphere harmonic.")
