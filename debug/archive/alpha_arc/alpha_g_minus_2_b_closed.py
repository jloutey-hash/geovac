"""
Attempt closed form for B(n_ext, n_int=1).

KEY SIMPLIFICATION: For (n_ext, n_int=1), vertex selection requires
n_ext + 1 + q odd, so q must have the same parity as n_ext.
The constraint |n_ext - 1| <= q <= n_ext + 1 gives q = n_ext as the
ONLY allowed value (since q_lo = n_ext-1 and q_hi = n_ext+1, and
(n_ext)+(1)+(n_ext) = 2*n_ext+1 is always odd).

Wait, let me verify more carefully:
  q_lo = max(1, |n_ext - 1|) = n_ext - 1 for n_ext >= 2, = 1 for n_ext = 1
  q_hi = n_ext + 1
  vertex_allowed: n_ext + 1 + q odd
    q = n_ext - 1: n_ext + 1 + n_ext - 1 = 2*n_ext -> EVEN -> NOT allowed
    q = n_ext:     n_ext + 1 + n_ext = 2*n_ext + 1 -> ODD -> ALLOWED
    q = n_ext + 1: n_ext + 1 + n_ext + 1 = 2*n_ext + 2 -> EVEN -> NOT allowed

So indeed q = n_ext is the ONLY allowed loop-photon mode.
mu_q = n_ext * (n_ext + 2)

This means B(n_ext, n_int=1) has the form:
  B = [sum of CG products] / (lam_int^4 * mu_q)
  = [CG sum] / ((5/2)^4 * n_ext*(n_ext+2))
  = [CG sum] * 16 / (625 * n_ext * (n_ext + 2))

The CG sum involves the tree-level vertex V(n_ext -> n_int=1; q=n_ext)
and the probe vertex V(n_int=1 -> n_int=1; q=1).

Let's compute the CG sums symbolically and look for patterns.
"""

from sympy import Rational, sqrt, simplify, S, cancel, factor
from sympy.physics.wigner import clebsch_gordan
import json

# Load exact results
with open('debug/data/alpha_g_minus_2_exact_forms.json') as f:
    data = json.load(f)
results = data['results']


def half_ints(j):
    result = []
    m = -j
    while m <= j + Rational(1, 100):
        result.append(m)
        m += 1
    return result


# The V_magnetic closed form
def V_magnetic_closed(n):
    return Rational(2, 3) * (sqrt(Rational(n+3, n+1)) - sqrt(Rational(n, n+2)))


# The propagator denominator
def propagator_denom(n_ext):
    """lam_int^4 * mu_q where n_int=1, q=n_ext"""
    lam_int = Rational(5, 2)
    mu_q = Rational(n_ext * (n_ext + 2))
    return lam_int**4 * mu_q  # = 625/16 * n_ext*(n_ext+2)


# Let's compute B * propagator_denom to get just the CG sum part
print("=== B * propagator_denom = pure CG sum ===")
for r in results:
    n = r['n_ext']
    denom = float(propagator_denom(n))
    b = r['B_nint1_float']
    cg_sum = b * denom
    print(f"  n_ext={n}: B={b:.6e}, denom={denom:.2f}, CG_sum = {cg_sum:.10f}")

# Let's also compute B * denom * V_mag to get lam^2*F2 * denom * V_mag^2
# Actually, the simplest object is:
# lam^2 * F2 * denom / V_mag
# = lam^2 * B / V_mag * denom / V_mag
# = lam^2 * (CG_sum / denom) / V_mag * denom / V_mag
# Hmm, that's circular.

# Better: let's look at CG_sum = B * lam_int^4 * mu_q directly
# and check if CG_sum / V_mag^2 is clean

print("\n=== CG_sum / V_magnetic^2 ===")
for r in results:
    n = r['n_ext']
    denom = float(propagator_denom(n))
    b = r['B_nint1_float']
    v_mag = float(V_magnetic_closed(n))
    cg_sum = b * denom
    ratio = cg_sum / v_mag**2
    print(f"  n_ext={n}: CG_sum/V_mag^2 = {ratio:.10f}")

# And CG_sum / V_mag
print("\n=== CG_sum / V_magnetic ===")
for r in results:
    n = r['n_ext']
    denom = float(propagator_denom(n))
    b = r['B_nint1_float']
    v_mag = float(V_magnetic_closed(n))
    cg_sum = b * denom
    ratio = cg_sum / v_mag
    print(f"  n_ext={n}: CG_sum/V_mag = {ratio:.10f}")

# And CG_sum * (n+1)*(n+2)
print("\n=== CG_sum * (n+1)*(n+2) ===")
for r in results:
    n = r['n_ext']
    denom = float(propagator_denom(n))
    b = r['B_nint1_float']
    cg_sum = b * denom
    ratio = cg_sum * (n+1) * (n+2)
    print(f"  n_ext={n}: CG_sum*(n+1)*(n+2) = {ratio:.10f}")

# CG_sum * n*(n+2) -- this removes the propagator mu_q
print("\n=== CG_sum * n*(n+2) = B * lam_int^4 * mu_q^2 ===")
for r in results:
    n = r['n_ext']
    denom = float(propagator_denom(n))
    b = r['B_nint1_float']
    cg_sum = b * denom
    ratio = cg_sum * n * (n+2)
    print(f"  n_ext={n}: CG_sum*n*(n+2) = {ratio:.10f}")

# Let's be more systematic. The exact B expressions have sqrt factors.
# B = CG_sum / (625/16 * n*(n+2))
# So CG_sum = B * 625/16 * n*(n+2)

print("\n\n=== Exact CG_sum expressions ===")
for r in results:
    n = r['n_ext']
    # Parse the exact B
    b_str = r['B_nint1_exact']
    from sympy import sympify
    b_exact = sympify(b_str)
    denom_exact = Rational(625, 16) * n * (n + 2)
    cg_exact = simplify(b_exact * denom_exact)
    print(f"  n_ext={n}: CG_sum = {cg_exact}")
    print(f"          = {float(cg_exact):.12f}")

# Now let's look at CG_sum * 9/4 / [(n+3)/(n+1) - n/(n+2)]
# = CG_sum * 9/4 / V_mag_sq_factor
# where V_mag = (2/3) * sqrt[X], X = (n+3)/(n+1) - n/(n+2) +cross term
# Actually V_mag^2 = (4/9)[(n+3)/(n+1) + n/(n+2) - 2*sqrt(n(n+3)/((n+1)(n+2)))]

# Let's simplify: look at B * (625/16) * n*(n+2) * (3/2) / V_mag
# = CG_sum * (3/2) / V_mag = F2 * (625/16) * n*(n+2) * (3/2) / 1
# Hmm, let me think again.

# F2 = B / V_mag, lam^2 * F2 -> c0
# So B ~ c0 * V_mag / lam^2 ~ c0 * (4/(3n)) / n^2 = 4*c0/(3*n^3)
# And CG_sum = B * (625/16) * n*(n+2) ~ (4*c0/(3*n^3)) * (625/16) * n^2
# = 4*c0*625 / (3*16*n) = 625*c0/(12*n) for large n

# Check: CG_sum * n ->?
print("\n=== CG_sum * n (should approach constant) ===")
for r in results:
    n = r['n_ext']
    denom = float(propagator_denom(n))
    b = r['B_nint1_float']
    cg_sum = b * denom
    ratio = cg_sum * n
    print(f"  n_ext={n}: CG_sum*n = {ratio:.10f}")

# CG_sum * n * 12/625 ->?
print("\n=== CG_sum * n * 12/625 (should approach c0) ===")
for r in results:
    n = r['n_ext']
    denom = float(propagator_denom(n))
    b = r['B_nint1_float']
    cg_sum = b * denom
    ratio = cg_sum * n * 12 / 625
    print(f"  n_ext={n}: CG_sum*n*12/625 = {ratio:.10f}")

# 625*c0/12 = ?
c0 = 0.0054313585
print(f"\n  625*c0/12 = {625*c0/12:.10f}")

# Now the KEY question: what is the exact CG sum structure?
# For n_ext=n, n_int=1, q=n: the loop photon has q=n
# jgL for loop = (n+1)/2, jgR = (n-1)/2 (or vice versa)
# Since q=n, the photon channels are:
#   (jgL, jgR) = ((n+1)/2, (n-1)/2) and ((n-1)/2, (n+1)/2)

# For vertex 1: n_ext -> n_int=1 via q=n
#   electron: (jL_ext, jR_ext) = ((n+1)/2, n/2), j_ext = 1/2
#   internal: (jL_int, jR_int) = (1, 1/2), j_int = 1/2
#   Photon: (jgL, jgR)
#   Triangle: |jL_ext - jgL| <= jL_int <= jL_ext + jgL
#             |(n+1)/2 - jgL| <= 1 <= (n+1)/2 + jgL
#   For jgL = (n+1)/2: |0| <= 1 <= n+1 -> YES (for n >= 1)
#   For jgL = (n-1)/2: |1| <= 1 <= n -> YES (1 <= 1 for n >= 2, |1| = 1 <= 1 for n=2)
#   And: |jR_ext - jgR| <= jR_int <= jR_ext + jgR
#        |n/2 - jgR| <= 1/2 <= n/2 + jgR
#   For jgR = (n-1)/2: |n/2 - (n-1)/2| = |1/2| = 1/2 <= 1/2 -> YES (tight!)
#   For jgR = (n+1)/2: |n/2 - (n+1)/2| = |(-1/2)| = 1/2 <= 1/2 -> YES (tight!)

# So both photon channels are allowed for vertex 1.
# The triangle inequality is tight in jR, meaning the vertex is maximally constrained.

print("\n\n=== Vertex 1 channel analysis: n_ext -> n_int=1 via q=n_ext ===")
for n in range(1, 8):
    jEL = Rational(n+1, 2)
    jER = Rational(n, 2)
    jIL = Rational(1)  # n_int=1
    jIR = Rational(1, 2)
    q = n

    print(f"\n  n_ext={n}: jE=({jEL},{jER}), jI=(1,1/2), q={q}")

    for jgL, jgR in [(Rational(q+1, 2), Rational(q-1, 2)),
                      (Rational(q-1, 2), Rational(q+1, 2))]:
        if jgL < 0 or jgR < 0:
            continue
        # Triangle checks
        check_L = (abs(jEL - jgL) <= jIL <= jEL + jgL)
        check_R = (abs(jER - jgR) <= jIR <= jER + jgR)
        print(f"    Photon ({jgL},{jgR}): L_tri={check_L} [{abs(jEL-jgL)}<=1<={jEL+jgL}], "
              f"R_tri={check_R} [{abs(jER-jgR)}<=1/2<={jER+jgR}]")

# The probe vertex: n_int=1 -> n_int=1 via q=1
# This is just the tree-level probe on the internal line
# V_probe_magnetic(n_int=1) = V_magnetic_closed(1) = 2/3*(sqrt(4/2) - sqrt(1/3))
# = 2/3*(sqrt(2) - 1/sqrt(3))
print(f"\n  Probe vertex V_mag(n_int=1) = {float(V_magnetic_closed(1)):.10f}")

# Since there's only one internal j value (j_int = 1/2 since n_int=1 -> jL=1, jR=1/2)
# and the probe vertex V(+1/2) - V(-1/2) is the same V_magnetic we already have,
# the B expression factorizes partially.

# Actually the full one-loop has structure:
# B = sum_q sum_{channels} sum_{m_int, m_int'} V1 * probe * V3 / (lam_int^4 * mu_q)
# Since q = n_ext (unique), and the probe is on n_int=1:
# B = (1/propagator) * sum_{channels, m} V1(n_ext->1; q=n_ext) *
#     probe(1->1; q=1) * V3(1->n_ext; q=n_ext)

# The probe vertex on n_int=1 is just V_magnetic(1) summed over probe polarizations
# But it's not just a number -- it depends on (j_int, m_int, m_int')

# Let me compute the "dressed probe" explicitly
print("\n\n=== Probe on n_int=1 (dressed) ===")
# j_int ranges from |jIL - jIR| to jIL + jIR = |1 - 1/2| to 3/2 = 1/2 to 3/2
# But j_ext = 1/2 (external electron), and our B computation
# sums over j_int and m_int as well

# For the g-2 extraction, we need:
# B(mj_ext=+1/2) - B(mj_ext=-1/2)
# The probe contributes V(j_int, m_int -> j_int, m_int'; probe)
# which is a matrix in (m_int, m_int') space

# Since n_int=1: j_int_min = |1 - 1/2| = 1/2, j_int_max = 1 + 1/2 = 3/2
# j_int can be 1/2 or 3/2

# For j_int = 1/2: this IS the j=1/2 sector we already computed
# For j_int = 3/2: this is a DIFFERENT sector

# B sums over BOTH j_int values!

print("  n_int=1: jIL=1, jIR=1/2")
print("  j_int values: 1/2, 3/2")
print("  The probe V(j_int=1/2) is our V_magnetic(1)")
print("  The probe V(j_int=3/2) is a separate contribution")

# This is getting complex. Let me just compute B_exact * propagator_denom
# symbolically and try to factor it.

print("\n\n=== Factoring exact CG_sum expressions ===")
for r in results:
    n = r['n_ext']
    from sympy import sympify
    b_exact = sympify(r['B_nint1_exact'])
    denom_exact = Rational(625, 16) * n * (n + 2)
    cg_exact = simplify(b_exact * denom_exact)

    # Try to express in terms of V_magnetic(n)
    v_mag = V_magnetic_closed(n)
    v_mag_1 = V_magnetic_closed(1)

    # Check CG_sum / (V_mag * V_mag_1)
    ratio = simplify(cg_exact / (v_mag * v_mag_1))
    print(f"  n_ext={n}: CG_sum / (V_mag(n)*V_mag(1)) = {ratio} = {float(ratio):.10f}")

# Check CG_sum / V_mag
print("\n=== CG_sum / V_mag(n_ext) ===")
for r in results:
    n = r['n_ext']
    b_exact = sympify(r['B_nint1_exact'])
    denom_exact = Rational(625, 16) * n * (n + 2)
    cg_exact = simplify(b_exact * denom_exact)
    v_mag = V_magnetic_closed(n)
    ratio = simplify(cg_exact / v_mag)
    print(f"  n_ext={n}: CG_sum / V_mag = {ratio} = {float(ratio):.10f}")

# CG_sum / (V(1,0)_n * V(1,0)_1)
print("\n=== CG_sum / (V(1,0)_n * V(1,0)_1) ===")
for r in results:
    n = r['n_ext']
    b_exact = sympify(r['B_nint1_exact'])
    denom_exact = Rational(625, 16) * n * (n + 2)
    cg_exact = simplify(b_exact * denom_exact)
    v10_n = Rational(2, 3) * sqrt(Rational(n+3, n+1))
    v10_1 = Rational(2, 3) * sqrt(Rational(4, 2))  # = 2*sqrt(2)/3
    ratio = simplify(cg_exact / (v10_n * v10_1))
    print(f"  n_ext={n}: CG_sum / (V10_n * V10_1) = {ratio} = {float(ratio):.10f}")

# Let's try: does CG_sum = a_n * V(1,0)_n + b_n * V(0,1)_n for some rational a_n, b_n?
print("\n=== CG_sum as linear combination of V(1,0) and V(0,1) ===")
for r in results:
    n = r['n_ext']
    b_exact = sympify(r['B_nint1_exact'])
    denom_exact = Rational(625, 16) * n * (n + 2)
    cg_exact = simplify(b_exact * denom_exact)

    v10 = Rational(2, 3) * sqrt(Rational(n+3, n+1))
    v01 = -Rational(2, 3) * sqrt(Rational(n, n+2))

    # Solve: CG = a * v10 + b * v01
    # This has two unknowns with one equation (unless CG has two sqrt terms)
    # Actually CG_exact has multiple sqrt terms typically
    print(f"  n_ext={n}: CG_sum = {cg_exact}")
    print(f"             V(1,0) = {v10} = {float(v10):.10f}")
    print(f"             V(0,1) = {v01} = {float(v01):.10f}")

    # Check if CG_sum has the same sqrt arguments as V_mag
    # V_mag has sqrt((n+3)/(n+1)) and sqrt(n/(n+2))
    # CG_sum from B has various sqrt terms...

print("\n\n=== DONE ===")
