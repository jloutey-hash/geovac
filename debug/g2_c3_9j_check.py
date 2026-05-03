"""
g2_c3_9j_check.py -- Verify if the 9j irrational is real or a convention artifact.

Standard 9j symbols are RATIONAL (sums of products of 6j symbols, which are
sums of products of factorials). Check if sympy's wigner_9j has √ factors.
"""

from sympy.physics.wigner import wigner_6j, wigner_9j, clebsch_gordan
from sympy import Rational as R, sqrt, simplify, S


# The 9j symbol {a,b,c; d,e,f; g,h,j} = Σ_k (-1)^{2k}(2k+1) × 6j×6j×6j
# Only k=1/2 contributes for our case

a, b, c = R(1), R(1,2), R(1,2)
d, e, f = R(1), R(0), R(1)
g, h, j = R(1), R(1,2), R(1,2)

print("=== 9j symbol verification ===")
print(f"9j({a},{b},{c}; {d},{e},{f}; {g},{h},{j})")

# Direct from sympy
w9j = wigner_9j(a, b, c, d, e, f, g, h, j)
print(f"  sympy result: {w9j}")
print(f"  float: {float(w9j.evalf(30)):.15e}")

# Manual computation using the standard 9j decomposition:
# {a b c}
# {d e f}  = Σ_k (-1)^(2k)(2k+1) {a b c}{d e f}{g h j}
# {g h j}                           {k j f}{k b h}{k d a}
#
# where each 6j symbol {j1 j2 j3; j4 j5 j6} has triads
# (j1,j2,j3), (j1,j5,j6), (j4,j2,j6), (j4,j5,j3)

print("\n=== Manual 6j decomposition ===")
# k ranges from max(|a-j|,|d-h|,|b-f|) to min(a+j, d+h, b+f)
# a=1,j=1/2 -> |1/2| to 3/2; d=1,h=1/2 -> 1/2 to 3/2; b=1/2,f=1 -> 1/2 to 3/2
# So k = 1/2, 3/2

for k_2 in [1, 3]:
    k = R(k_2, 2)
    phase = (-1)**(2*k_2) * (2*k + 1)
    print(f"\n  k={k}: phase_factor = (-1)^{2*k_2} * {2*k+1} = {phase}")

    # 6j{a,b,c; k,j,f} = {1, 1/2, 1/2; k, 1/2, 1}
    try:
        s6j_1 = wigner_6j(a, b, c, k, j, f)
        print(f"    6j{{1,1/2,1/2; {k},1/2,1}} = {s6j_1}  rational? {s6j_1.is_Rational}")
    except Exception as ex:
        print(f"    6j_1 error: {ex}")
        s6j_1 = S(0)

    # 6j{d,e,f; k,b,h} = {1, 0, 1; k, 1/2, 1/2}
    try:
        s6j_2 = wigner_6j(d, e, f, k, b, h)
        print(f"    6j{{1,0,1; {k},1/2,1/2}} = {s6j_2}  rational? {s6j_2.is_Rational}")
    except Exception as ex:
        print(f"    6j_2 error: {ex}")
        s6j_2 = S(0)

    # 6j{g,h,j; k,d,a} = {1, 1/2, 1/2; k, 1, 1}
    try:
        s6j_3 = wigner_6j(g, h, j, k, d, a)
        print(f"    6j{{1,1/2,1/2; {k},1,1}} = {s6j_3}  rational? {s6j_3.is_Rational}")
    except Exception as ex:
        print(f"    6j_3 error: {ex}")
        s6j_3 = S(0)

    product = phase * s6j_1 * s6j_2 * s6j_3
    print(f"    contribution = {simplify(product)}  rational? {simplify(product).is_Rational}")

# Try alternative decomposition formula from Varshalovich
# 9j{j1 j2 j3; j4 j5 j6; j7 j8 j9}
# = Σ_k (-1)^(2k)(2k+1) {j1 j4 j7}{j2 j5 j8}{j3 j6 j9}
#                          {j8 k j6}{j4 k j6... no, let me use the standard one

# Actually the standard decomposition (Varshalovich Ch. 10 Eq. 1):
# {j11 j12 j13}
# {j21 j22 j23} = Σ_x (-1)^(2x)(2x+1)
# {j31 j32 j33}
#   × {j11 j12 j13} × {j21 j22 j23} × {j31 j32 j33}
#     {j32 x  j23}    {j13 x  j31}    {j21 x  j12}

print("\n=== Varshalovich decomposition ===")
j11, j12, j13 = a, b, c  # 1, 1/2, 1/2
j21, j22, j23 = d, e, f  # 1, 0, 1
j31, j32, j33 = g, h, j  # 1, 1/2, 1/2

# x range: governed by triangle conditions in the three 6j symbols
# {j11,j12,j13; j32,x,j23} -> triads include (j32,x,j23)=(1/2,x,1): x in {1/2, 3/2}
# Also (j11,x,j23)=(1,x,1): x in {0,1,2}; (j12,j32,x)=(1/2,1/2,x): x in {0,1}
# Intersection: x in {1/2, 3/2} ∩ {0,1,2} ∩ {0,1} -> need half-int AND int, impossible
# Wait, the x must be same in all three. Let me recheck triangle rules.

# Actually x ranges: from each 6j we get constraints
# 6j{j11,j12,j13; j32,x,j23}: triads (j11,j12,j13)✓, (j11,x,j23)=(1,x,1), (j32,j12,j23)=(1/2,1/2,1)✓, (j32,x,j13)=(1/2,x,1/2)
#   From (1,x,1): |0|<=x<=2; from (1/2,x,1/2): |0|<=x<=1
#   So x ∈ {0, 1} (integer) or x ∈ {1/2} (half-int)?  No - (1,x,1) requires x integer or half-int with 1+x+1 integer -> x integer.
#   And (1/2,x,1/2) requires 1/2+x+1/2 integer -> x integer. So x ∈ {0, 1}.

# 6j{j21,j22,j23; j13,x,j31}: (1,0,1; 1/2,x,1)
#   Triads: (1,0,1)✓, (1,x,1): x∈{0,1,2}, (1/2,0,1)=FAIL unless... wait (j13,j22,j31)=(1/2,0,1): need |1/2-0|<=1<=1/2+0 -> 1<=1/2, FAIL
#   So this 6j = 0 always! Because triad (j13,j22,j31)=(1/2,0,1) violates triangle: |1/2-0|=1/2, max=1/2, but j31=1 > 1/2.

# Hmm, that would make the 9j=0, but sympy says sqrt(6)/18. Let me re-check which decomposition formula to use.

# The CORRECT formula (Varshalovich 10.2.1):
# {j1 j2 j3}
# {j4 j5 j6}  = Σ_g (-1)^(2g)(2g+1) {j1 j2 j3}{j4 j5 j6}{j7 j8 j9}
# {j7 j8 j9}                            {j6 j9 g}{j2 g j8}{j4 g j1... no

# Let me just use the simplest correct formula. The 9j in terms of 6j is:
# 9j(j1..j9) = Σ_g (-1)^(2g)(2g+1) {j1,j4,j7}{j2,j5,j8}{j3,j6,j9}
#                                     {j8,j9,g}{j4,g,j6}{g,j1,j2}...

# This is getting confusing with conventions. Let me just verify numerically.
print("  (Skipping manual Varshalovich — verifying numerically instead)")
print(f"\n  sympy 9j = {w9j}")
print(f"  sympy 9j is_Rational = {w9j.is_Rational}")
print(f"  sympy 9j**2 = {simplify(w9j**2)}")
print(f"  sympy 9j**2 is_Rational = {simplify(w9j**2).is_Rational}")

# KEY TEST: is 9j² rational? If yes, the √ comes from the 9j convention
# Standard 9j symbols ARE sums of products of 6j symbols (all rational).
# But the UNITARY 9j (used in some conventions) includes a sqrt(2j+1) normalization.
# Check if sympy uses the standard or unitary convention.
print(f"\n  9j² = {simplify(w9j**2)} = {float(simplify(w9j**2).evalf())}")
print(f"  1/54 = {float(R(1,54))} ... 9j² = 1/54? {simplify(w9j**2 - R(1,54)) == 0}")
print(f"  6/324 = 1/54, so 9j = sqrt(6)/18 means 9j² = 6/324 = 1/54")
print(f"  Check: {R(6, 324)} = {R(6, 324) == R(1, 54)}")

# Compare
diff = simplify(w9j - product)
print(f"\n  sympy - manual = {diff}")

# Now check ALL relevant 9j symbols
print("\n=== All 9j symbols for the vertex ===")
cases = [
    # (jgL, jgR, jg_tot, j_int)
    (1, 0, 1, R(1,2)),
    (1, 0, 1, R(3,2)),
    (0, 1, 1, R(1,2)),
    (0, 1, 1, R(3,2)),
]

for jgL, jgR, jg, j_int in cases:
    w = wigner_9j(R(1), R(1,2), R(1,2),
                   R(jgL), R(jgR), R(jg),
                   R(1), R(1,2), j_int)
    print(f"  9j(1,1/2,1/2; {jgL},{jgR},{jg}; 1,1/2,{j_int}) = {w}  rational? {w.is_Rational}")

# Key test: does the vertex REDUCED matrix element involve only 6j symbols?
# The vertex amplitude is:
# Σ_mL CG(jL,mL;jR,mR|j,m) × CG(jL',mL';jR',mR'|j',m') × CG(jL,mL;jgL,mgL|jL',mL') × CG(jR,mR;jgR,mgR|jR',mR')
# = sqrt((2j+1)(2j'+1)) × CG(j,m;jg_tot,mg|j',m') × 9j{jL,jR,j; jgL,jgR,jg; jL',jR',j'}
# where the CG couples j and jg_tot to j'

# BUT: the photon is NOT coupled to a definite jg_tot in our diagram!
# The photon has (jgL, mgL, jgR, mgR) separately, not a total (jg, mg).
# So we need to couple the photon first:
# V(j,m; j',m'; jgL,mgL,jgR,mgR) = Σ_{jg,mg} CG(jgL,mgL;jgR,mgR|jg,mg) × V_reduced(j,j';jg)

print("\n=== Factored vertex: couple photon first ===")
# For n_ext=n_int=1, the vertex v1 involves:
# V = Σ_{jg} CG(jgL,mgL;jgR,mgR|jg,mg) × CG(j_ext,mj_ext;jg,mg|j_int,mj_int)
#     × sqrt((2j_ext+1)(2j_int+1)) × 9j{jEL,jER,j_ext; jgL,jgR,jg; jIL,jIR,j_int}

# The 9j is the reduced matrix element. The CGs handle the m-dependence.
# Test: compute V for specific m values both ways

jEL, jER = R(1), R(1,2)
jIL, jIR = R(1), R(1,2)
j_ext = R(1,2)

# Pick mj_ext=+1/2, mj_int=+1/2
mj_ext_2 = +1
mj_int_2 = +1

# Channel (jgL=1, jgR=0):
jgL, jgR = R(1), R(0)
mgL_2 = 0  # try mgL=0, mgR=0
mgR_2 = 0

# Direct computation
from debug.g2_c3_racah import vertex_amp_exact
v_direct = vertex_amp_exact(2, 1, 1, mj_ext_2,
                             2, 1, 1, mj_int_2,
                             2, 0, mgL_2, mgR_2)
print(f"  Direct vertex (mgL=0, mgR=0): {v_direct}")

# Factored: couple photon (jgL=1, jgR=0) to total jg
# CG(1,0; 0,0 | jg, mg) = δ_{jg,1} δ_{mg,0}
# So jg=1, mg=0 only.
jg = R(1)
mg_2 = 0

cg_photon = clebsch_gordan(R(1), R(0), R(1), R(0), R(0), R(0))
cg_ext = clebsch_gordan(j_ext, jg, R(1,2), R(mj_ext_2, 2), R(mg_2, 2), R(mj_int_2, 2))
w9 = wigner_9j(jEL, jER, j_ext, jgL, jgR, jg, jIL, jIR, R(1,2))
dim_factor = sqrt((2*j_ext + 1) * (2*R(1,2) + 1))

v_factored = cg_photon * cg_ext * dim_factor * w9
print(f"  CG_photon = {cg_photon}")
print(f"  CG_ext = {cg_ext}")
print(f"  9j = {w9}")
print(f"  dim = {dim_factor}")
print(f"  Factored: {simplify(v_factored)}")
print(f"  Direct: {v_direct}")
print(f"  Ratio factored/direct: {simplify(v_factored/v_direct) if v_direct != 0 else 'N/A'}")
