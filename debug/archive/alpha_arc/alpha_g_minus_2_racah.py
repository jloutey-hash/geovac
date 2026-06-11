"""
Racah/6j reduction of the one-loop vertex correction B(n_ext, n_int=1).

The vertex amplitude V(s -> t; photon) involves 4 CG coefficients:
  V = sum_{mL1, mR1} C(jL_s, jR_s, j_s; mL1, m_s-mL1, m_s)
                     * C(jL_t, jR_t, j_t; mL1+mgL, m_s-mL1+mgR, m_t)
                     * C(jL_s, jgL, jL_t; mL1, mgL, mL1+mgL)
                     * C(jR_s, jgR, jR_t; m_s-mL1, mgR, m_s-mL1+mgR)

The one-loop B sums products of such vertices over internal m quantum numbers
and photon polarizations. This has the structure of a recoupling, expressible
via 6j or 9j symbols.

KEY: For n_int=1, the internal electron has (jL_int, jR_int) = (1, 1/2).
The loop photon with q = n_ext has (jgL, jgR) channels.
The probe photon has q = 1.

Let me use the Wigner-Eckart theorem to extract the REDUCED matrix elements
and then express B in terms of 6j symbols.

Actually, simpler: each vertex amplitude is a product of TWO recoupling
coefficients (one for L sector, one for R sector), each of which is a
CG coefficient. The sum over magnetic quantum numbers of a product of
CG coefficients gives 6j symbols.

The vertex V(jL_s, jR_s, j_s, m_s -> jL_t, jR_t, j_t, m_t; jgL, jgR, mgL, mgR)
factorizes as:
  V = sum_{mL1} C(jL_s, jR_s, j_s)_{mL1, m_s-mL1, m_s}
             * C(jL_t, jR_t, j_t)_{mL1+mgL, m_t-mL1-mgL, m_t}
             * C(jL_s, jgL, jL_t)_{mL1, mgL, mL1+mgL}
             * C(jR_s, jgR, jR_t)_{m_s-mL1, mgR, m_t-mL1-mgL}

The key observation: for a DIAGONAL (self-energy type) vertex where source=target,
the sum over all photon polarizations (mgL, mgR) of |V|^2 reduces to a product
of 6j symbols.

For the magnetic moment, we need V_up - V_dn (not |V|^2), so it's the linear
vertex amplitude, not the squared one.

Let me think about this differently. The FULL one-loop diagram has:
  B(m_ext) = sum_{q, channels} sum_{m_int, m_int'}
             V_1(ext, m_ext -> int, m_int; loop_photon)
             * V_probe(int, m_int -> int, m_int'; probe q=1)
             * V_2(int, m_int' -> ext, m_ext; loop_photon)
             / (lambda_int^4 * mu_q)

V_1 and V_2 share the SAME loop photon, so there's also a sum over
photon polarizations (sum_{mgL, mgR}).

For n_int=1, q=n_ext (unique), the structure is:
  B(m_ext) = sum over all internal m's and photon m's of:
    [CG product vertex1] * [CG product probe] * [CG product vertex2]
    / ((5/2)^4 * n*(n+2))

Since vertex1 goes (n_ext -> 1) and vertex2 goes (1 -> n_ext) with
the SAME photon quantum numbers, vertex2 is the "time-reverse" of vertex1.

For a Hermitian theory, V_2(int -> ext; q) = V_1(ext -> int; q)^*
but since all CG coefficients are real, V_2 = V_1.

Actually, the one-loop self-energy (without probe) would give:
  sum_{mgL, mgR} V(ext->int) * V(int->ext) = sum |V|^2
which is manifestly positive.

The probe vertex breaks this structure by inserting a different operator
between vertex1 and vertex2.

Let me just compute more data points (up to n=15) using the existing
CG-based code and look for patterns in the exact expressions.
"""

from sympy import (Rational, sqrt, simplify, S, sympify, nsimplify,
                   cancel, factor, collect, radsimp, together, oo)
from sympy.physics.wigner import clebsch_gordan, wigner_6j, racah
import json
import time


def half_ints(j):
    result = []
    m = -j
    while m <= j + Rational(1, 100):
        result.append(m)
        m += 1
    return result


def V_magnetic_closed(n):
    return Rational(2, 3) * (sqrt(Rational(n+3, n+1)) - sqrt(Rational(n, n+2)))


def V10(n):
    return Rational(2, 3) * sqrt(Rational(n+3, n+1))


def V01(n):
    return -Rational(2, 3) * sqrt(Rational(n, n+2))


def vertex_amp_pol(j_sL, j_sR, j_s, mj_s,
                   j_tL, j_tR, j_t, mj_t,
                   jgL, jgR, mgL, mgR):
    """Vertex amplitude for specific polarizations."""
    total = S.Zero
    for mL1 in half_ints(j_sL):
        mR1 = mj_s - mL1
        if abs(mR1) > j_sR:
            continue
        mL2 = mL1 + mgL
        if abs(mL2) > j_tL:
            continue
        mR2 = mR1 + mgR
        if abs(mR2) > j_tR:
            continue
        if mL2 + mR2 != mj_t:
            continue
        c1 = clebsch_gordan(j_sL, j_sR, j_s, mL1, mR1, mj_s)
        if c1 == 0:
            continue
        c2 = clebsch_gordan(j_tL, j_tR, j_t, mL2, mR2, mj_t)
        if c2 == 0:
            continue
        c3 = clebsch_gordan(j_sL, jgL, j_tL, mL1, mgL, mL2)
        c4 = clebsch_gordan(j_sR, jgR, j_tR, mR1, mgR, mR2)
        total += c1 * c2 * c3 * c4
    return total


def vertex_allowed(n1, n2, q):
    if q < 1 or q < abs(n1 - n2) or q > n1 + n2:
        return False
    return (n1 + n2 + q) % 2 == 1


def compute_B_nint1(n_ext, j_ext=Rational(1, 2)):
    """Compute B(n_ext, n_int=1) exactly using CG coefficients.

    B = sum over photon channels and m's of:
      V1(ext,m_ext -> int,m_int; loop) * V_probe(int,m_int -> int,m_int') * V2(int,m_int' -> ext,m_ext; loop)
      / (lam_int^4 * mu_q)

    For n_int=1: jIL=1, jIR=1/2, lam_int=5/2, q=n_ext, mu_q=n_ext*(n_ext+2)
    """
    n_int = 1
    jEL = Rational(n_ext + 1, 2)
    jER = Rational(n_ext, 2)
    jIL = Rational(1)
    jIR = Rational(1, 2)

    lam_int = Rational(5, 2)
    q = n_ext
    mu_q = Rational(q * (q + 2))
    prop = lam_int**4 * mu_q

    # Loop photon channels for (n_ext -> 1; q=n_ext)
    loop_channels = []
    for jgL, jgR in [(Rational(q+1, 2), Rational(q-1, 2)),
                      (Rational(q-1, 2), Rational(q+1, 2))]:
        if jgL < 0 or jgR < 0:
            continue
        # Check triangle for vertex1: jEL+jgL->jIL and jER+jgR->jIR
        if not (abs(jEL - jgL) <= jIL <= jEL + jgL):
            continue
        if not (abs(jER - jgR) <= jIR <= jER + jgR):
            continue
        loop_channels.append((jgL, jgR))

    # Probe photon channels for (1 -> 1; q=1)
    probe_channels = []
    for jpL, jpR in [(Rational(1), Rational(0)),
                      (Rational(0), Rational(1))]:
        if not (abs(jIL - jpL) <= jIL <= jIL + jpL):
            continue
        if not (abs(jIR - jpR) <= jIR <= jIR + jpR):
            continue
        probe_channels.append((jpL, jpR))

    # j_int ranges for the internal electron on n_int=1 shell
    # j_int from |jIL - jIR| to jIL + jIR = 1/2 to 3/2
    j_int_vals = [Rational(1, 2), Rational(3, 2)]

    B_up = S.Zero
    B_dn = S.Zero

    for jgL, jgR in loop_channels:
        for mgL in half_ints(jgL):
            for mgR in half_ints(jgR):
                for jpL, jpR in probe_channels:
                    for mpL in half_ints(jpL):
                        for mpR in half_ints(jpR):
                            for j_int in j_int_vals:
                                for m_int in half_ints(j_int):
                                    for m_int_p in half_ints(j_int):
                                        # Vertex 1: (ext, +1/2) -> (int, m_int)
                                        v1_up = vertex_amp_pol(
                                            jEL, jER, j_ext, Rational(1, 2),
                                            jIL, jIR, j_int, m_int,
                                            jgL, jgR, mgL, mgR)
                                        if v1_up == 0:
                                            continue

                                        # Probe: (int, m_int) -> (int, m_int')
                                        v_probe = vertex_amp_pol(
                                            jIL, jIR, j_int, m_int,
                                            jIL, jIR, j_int, m_int_p,
                                            jpL, jpR, mpL, mpR)
                                        if v_probe == 0:
                                            continue

                                        # Vertex 2: (int, m_int') -> (ext, +1/2)
                                        v2_up = vertex_amp_pol(
                                            jIL, jIR, j_int, m_int_p,
                                            jEL, jER, j_ext, Rational(1, 2),
                                            jgL, jgR, mgL, mgR)

                                        B_up += v1_up * v_probe * v2_up

                                        # Same for ext m = -1/2
                                        v1_dn = vertex_amp_pol(
                                            jEL, jER, j_ext, Rational(-1, 2),
                                            jIL, jIR, j_int, m_int,
                                            jgL, jgR, mgL, mgR)
                                        if v1_dn == 0:
                                            continue

                                        v2_dn = vertex_amp_pol(
                                            jIL, jIR, j_int, m_int_p,
                                            jEL, jER, j_ext, Rational(-1, 2),
                                            jgL, jgR, mgL, mgR)

                                        B_dn += v1_dn * v_probe * v2_dn

    B_mag = simplify((B_up - B_dn) / prop)
    return B_mag


# Compute B for n_ext = 1..10
print("=== Computing B(n_ext, n_int=1) for n_ext = 1..10 ===")
print("(This uses exact sympy CG coefficients)")
print()

all_data = []
for n in range(1, 11):
    t0 = time.time()
    B = compute_B_nint1(n)
    dt = time.time() - t0

    V_mag = V_magnetic_closed(n)
    F2 = simplify(B / V_mag)
    lam = Rational(2*n + 3, 2)
    lam2_F2 = simplify(lam**2 * F2)

    # CG_sum = B * propagator
    prop = Rational(625, 16) * n * (n + 2)
    cg_sum = simplify(B * prop)

    # CG_sum / V_mag^2
    cg_over_v2 = simplify(cg_sum / V_mag**2)

    # CG_sum * (n+1)*(n+2)
    cg_prod = simplify(cg_sum * (n+1) * (n+2))

    print(f"n_ext={n:2d}: B = {float(B):.10e}, V_mag = {float(V_mag):.10f}, "
          f"lam^2*F2 = {float(lam2_F2):.10f} ({dt:.1f}s)")
    print(f"         CG_sum = {float(cg_sum):.12f}, CG/V^2 = {float(cg_over_v2):.10f}, "
          f"CG*(n+1)(n+2) = {float(cg_prod):.10f}")

    all_data.append({
        'n': n,
        'B': float(B),
        'B_exact': str(B),
        'V_mag': float(V_mag),
        'F2': float(F2),
        'lam2_F2': float(lam2_F2),
        'cg_sum': float(cg_sum),
        'cg_sum_exact': str(cg_sum),
        'cg_over_v2': float(cg_over_v2),
        'cg_prod': float(cg_prod),
    })

# Analysis with 10 data points
print("\n\n=== CG_sum / V_mag^2: linearity analysis (10 points) ===")
for d in all_data:
    print(f"  n={d['n']:2d}: {d['cg_over_v2']:.12f}")

print("\n  Differences:")
for i in range(1, len(all_data)):
    diff = all_data[i]['cg_over_v2'] - all_data[i-1]['cg_over_v2']
    print(f"  n={all_data[i]['n']}: diff = {diff:.12f}")

# Richardson extrapolation on the differences
diffs = [all_data[i]['cg_over_v2'] - all_data[i-1]['cg_over_v2'] for i in range(1, len(all_data))]
print(f"\n  Last diff = {diffs[-1]:.12f}")
print(f"  1/6.3 = {1/6.3:.12f}")
print(f"  Trying to identify limiting slope...")
# For large n, CG/V^2 ~ a + b*n + c/n + ...
# So diff ~ b + c*(1/n - 1/(n-1)) + ... = b - c/n^2 + ...
# Richardson: b = diff_n + n^2 * (diff_n - diff_{n-1}) / (2n - 1)
for i in range(2, len(diffs)):
    n = i + 2
    rich = diffs[i] + n**2 * (diffs[i] - diffs[i-1]) / (2*n - 1)
    print(f"  n={n}: diff={diffs[i]:.12f}, Richardson_b = {rich:.12f}")


print("\n\n=== CG_sum*(n+1)*(n+2): linearity analysis (10 points) ===")
for d in all_data:
    print(f"  n={d['n']:2d}: {d['cg_prod']:.12f}")

print("\n  Differences:")
prods = [d['cg_prod'] for d in all_data]
for i in range(1, len(prods)):
    diff = prods[i] - prods[i-1]
    print(f"  n={all_data[i]['n']}: diff = {diff:.12f}")

prod_diffs = [prods[i] - prods[i-1] for i in range(1, len(prods))]
# Richardson on the product differences
print("\n  Richardson extrapolation on slope of CG*(n+1)(n+2):")
for i in range(2, len(prod_diffs)):
    n = i + 2
    rich = prod_diffs[i] + n**2 * (prod_diffs[i] - prod_diffs[i-1]) / (2*n - 1)
    print(f"  n={n}: diff={prod_diffs[i]:.12f}, Richardson = {rich:.12f}")


print("\n\n=== lam^2 * F2 convergence (10 points) ===")
for d in all_data:
    print(f"  n={d['n']:2d}: lam^2*F2 = {d['lam2_F2']:.12f}")

# Richardson on lam^2*F2
l2f2 = [d['lam2_F2'] for d in all_data]
print("\n  Neville table for lam^2*F2:")
# Build Richardson/Neville table
table = list(l2f2)
for k in range(1, len(table)):
    new_table = []
    for i in range(len(table) - 1):
        n_i = i + 1 + k - 1
        n_ik = i + 1
        # Standard Richardson: T_{k+1} = (n_i^2 * T_k[i+1] - n_ik^2 * T_k[i]) / (n_i^2 - n_ik^2)
        # But this isn't quite right for our sequence. Let me use Aitken/Shanks.
        pass
    break

# Simple Richardson assuming error ~ 1/n:
# T_extrap = (n2*T(n2) - n1*T(n1)) / (n2 - n1) with n1, n2 = consecutive
print("\n  Richardson (error ~ 1/n):")
for i in range(1, len(l2f2)):
    n1 = i
    n2 = i + 1
    extrap = (n2 * l2f2[i] - n1 * l2f2[i-1]) / (n2 - n1)
    print(f"  n=({n1},{n2}): extrap = {extrap:.12f}")

# Richardson assuming error ~ 1/n^2:
print("\n  Richardson (error ~ 1/n^2):")
for i in range(1, len(l2f2)):
    n1 = i
    n2 = i + 1
    extrap = (n2**2 * l2f2[i] - n1**2 * l2f2[i-1]) / (n2**2 - n1**2)
    print(f"  n=({n1},{n2}): extrap = {extrap:.12f}")


# Now try the EXACT ratio CG_sum / V(1,0) for n=6 (simplest case)
# and look for 6j structure
print("\n\n=== 6j symbol test ===")
# For the loop vertex (n_ext -> 1; q=n_ext), the CG sum has the structure
# of a recoupling coefficient. Let's see if we can express it via 6j.

# The tree-level probe on n_ext gives V(1,0) and V(0,1).
# These involve 4 CG products summed over mL.
# The sum over mL of C(jL,jR,j;mL,mR,m) * C(jL,jgL,jL';mL,mgL,mL')
# is a 6j-like object.

# Specifically, the vertex amplitude is:
# V = sum_mL C(jL_s,jR_s,j;mL,mR,m) * C(jL_s,jgL,jL_t;mL,mgL,mL+mgL) *
#     C(jR_s,jgR,jR_t;mR,mgR,mR+mgR) * C(jL_t,jR_t,j;mL+mgL,mR+mgR,m+mgL+mgR)

# For the PROBE (n -> n, q=1, same j_ext):
# This is a self-coupling, and we can use the Wigner-Eckart theorem:
# The sum over probe polarizations mgL, mgR of the probe vertex
# gives a tensor operator acting on the external state.

# For the magnetic moment, we need the rank-1 (vector) component
# proportional to m_j (the Zeeman-like term).

# Let me try a direct computation: for each n_ext, compute
# V_mag^2 * (something) and see if it matches a 6j expression.

# The Wigner-Eckart theorem says:
# <j m | T^k_q | j m'> = (-1)^{j-m} * 3j(j,k,j; -m,q,m') * <j || T^k || j>
# For k=1, q=0: <j,1/2|T^1_0|j,1/2> - <j,-1/2|T^1_0|j,-1/2>
# = (-1)^{j-1/2} * 3j(j,1,j;-1/2,0,1/2) * <j||T^1||j>
#   - (-1)^{j+1/2} * 3j(j,1,j;1/2,0,-1/2) * <j||T^1||j>

# For j=1/2:
# 3j(1/2,1,1/2; -1/2,0,1/2) = (-1)^{1/2-1/2+0} / sqrt(3) * ...
# Actually let me just compute the 3j symbols
from sympy.physics.wigner import wigner_3j
val = wigner_3j(Rational(1,2), 1, Rational(1,2), Rational(-1,2), 0, Rational(1,2))
print(f"  3j(1/2,1,1/2; -1/2,0,1/2) = {val} = {float(val):.10f}")
val2 = wigner_3j(Rational(1,2), 1, Rational(1,2), Rational(1,2), 0, Rational(-1,2))
print(f"  3j(1/2,1,1/2; 1/2,0,-1/2) = {val2} = {float(val2):.10f}")

# For j=1/2: V_up - V_dn = 2 * reduced_ME * 3j(1/2,1,1/2;-1/2,0,1/2)
# (up to phase)

# The key observation: V_magnetic = V_up - V_dn is proportional to
# the REDUCED matrix element of the rank-1 tensor part of the vertex.
# For the tree-level vertex, this gives a 6j symbol.

# Let me compute the tree-level V_mag using 6j and check:
print("\n=== Tree-level V_mag from 6j ===")
for n in range(1, 8):
    jL = Rational(n+1, 2)
    jR = Rational(n, 2)
    j_ext = Rational(1, 2)

    # For probe q=1, the photon channels are (1,0) and (0,1)
    v_mag_total = S.Zero

    for jpL, jpR in [(Rational(1), Rational(0)), (Rational(0), Rational(1))]:
        if jpL < 0 or jpR < 0:
            continue
        if not (abs(jL - jpL) <= jL <= jL + jpL):
            continue
        if not (abs(jR - jpR) <= jR <= jR + jpR):
            continue

        # The reduced matrix element for the rank-1 part involves a 6j symbol
        # <(jL,jR)j || T^1(jpL,jpR) || (jL,jR)j>
        # For a tensor product coupling:
        # T^k = [T^{jpL} x T^{jpR}]^k
        # The reduced ME involves:
        # sum over j' of 6j{jL,jR,j; jL,jR,j'} * ...

        # Actually, the direct CG computation is already fast.
        # Let me just verify the closed forms match.
        pass

    v_formula = V_magnetic_closed(n)
    print(f"  n={n}: V_mag(formula) = {float(v_formula):.10f}")

# Try: CG_sum = V_mag(n) * R(n) where R(n) is rational * n * stuff
# From the data: CG/V_mag approaches ~0.2 slowly
# CG * (n+1)(n+2) = p + q*n + r/n where q ~0.280

# If CG*(n+1)(n+2) -> p + q*n for large n, then:
# CG -> (p + q*n) / ((n+1)(n+2)) ~ q/n for large n
# F2 = CG / (prop * V_mag) = CG * 16 / (625 * n * (n+2) * V_mag)
# ~ (q/n) * 16 / (625 * n^2 * 4/(3n))
# = q * 16 * 3n / (625 * n^3 * 4)
# = q * 12 / (625 * n)   ... wait, this gives F2 ~ 1/n, so lam^2*F2 ~ n
# That's wrong. Let me redo:

# Actually lam = n + 3/2
# F2 = B / V_mag = CG_sum / (prop * V_mag)
# prop = 625/16 * n*(n+2)
# V_mag ~ 4/(3n)
# CG ~ q/(n) [from CG*(n+1)(n+2) ~ q*n, so CG ~ q*n/((n+1)(n+2)) ~ q/n]
# Wait: CG*(n+1)(n+2) ~ q*n means CG ~ q*n/n^2 = q/n

# F2 = CG / (625/16 * n*(n+2) * 4/(3n))
# = (q/n) / (625/16 * n * n * 4/(3n))
# = (q/n) / (625*4*n / (16*3))
# = (q/n) * 48 / (2500*n)
# = 48*q / (2500 * n^2)

# lam^2 * F2 ~ n^2 * 48*q / (2500 * n^2) = 48*q / 2500

# So c0 = 48*q / 2500 = 12*q/625

print(f"\n\n=== Prediction: c0 = 12*q/625 ===")
# Need q = lim_{n->inf} [CG*(n+1)*(n+2) - CG*(n)*(n+1)]
# From Richardson on prod_diffs

# Use the 10-point data for better Richardson
# The prod_diffs converge from below toward a limit
# Let's do careful Richardson

print("All computed lam^2*F2 values:")
for d in all_data:
    print(f"  n={d['n']:2d}: lam^2*F2 = {d['lam2_F2']:.12f}")

# Shanks transform
def shanks(seq):
    """Shanks transformation to accelerate convergence."""
    n = len(seq)
    if n < 3:
        return seq[-1]
    results = []
    for i in range(n - 2):
        denom = (seq[i+2] - seq[i+1]) - (seq[i+1] - seq[i])
        if abs(denom) < 1e-30:
            results.append(seq[i+2])
        else:
            results.append(seq[i+2] - (seq[i+2] - seq[i+1])**2 / denom)
    return results

print("\n  Shanks transform on lam^2*F2:")
s1 = shanks(l2f2)
for i, v in enumerate(s1):
    print(f"  S1[{i}] = {v:.12f}")

if len(s1) >= 3:
    s2 = shanks(s1)
    print("  Shanks^2:")
    for i, v in enumerate(s2):
        print(f"  S2[{i}] = {v:.12f}")

    if len(s2) >= 3:
        s3 = shanks(s2)
        print("  Shanks^3:")
        for i, v in enumerate(s3):
            print(f"  S3[{i}] = {v:.12f}")


# Also try to identify the constant with more precision
print("\n\n=== Identifying c0 ===")
# Best estimate from Shanks
if len(s3) > 0:
    c0_est = s3[-1]
elif len(s2) > 0:
    c0_est = s2[-1]
else:
    c0_est = s1[-1]
print(f"  Best c0 estimate = {c0_est:.15f}")
print(f"  1/c0 = {1/c0_est:.10f}")

# Test candidates
import math
candidates = [
    ("1/180", 1/180),
    ("1/184", 1/184),
    ("1/185", 1/185),
    ("pi/576", math.pi/576),
    ("pi/578", math.pi/578),
    ("pi/580", math.pi/580),
    ("1/(12*pi)", 1/(12*math.pi)),
    ("1/(36*pi)", 1/(36*math.pi)),
    ("4/(9*125)", 4/1125),
    ("16/2943", 16/2943),
    ("1/sqrt(33880)", 1/math.sqrt(33880)),
    ("8/(1125*sqrt(3))", 8/(1125*math.sqrt(3))),
    ("1/(25*sqrt(30))", 1/(25*math.sqrt(30))),
]
for name, val in candidates:
    rel_err = (val - c0_est) / c0_est
    print(f"  {name:20s} = {val:.12f}, rel err = {rel_err:.6e}")
