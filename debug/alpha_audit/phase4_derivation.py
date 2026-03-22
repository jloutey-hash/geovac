#!/usr/bin/env python3
"""
Phase 4: Derivation Audit — Can the formula be derived from first principles?
==============================================================================
Phase 3 found two leads:
  (A) det'(S1) * det'(S3) / pi = 4pi^2 exp(zeta(3)/(2pi^2)) ~ 41.957
      is 0.1% from B = 42.
  (B) The cubic is the eigenvalue equation of a traceless Z3-symmetric
      circulant with bc = K/3 and b^3 + c^3 = -1.

This phase investigates:
  A) Can the determinant gap be closed exactly?
  B) Can the circulant be grounded in the Hopf bundle?
  C) Assessment of the full derivation chain.
"""

import numpy as np
import os
from math import gamma

# ================================================================
# Constants (high precision where possible)
# ================================================================
CODATA_ALPHA_INV = 137.035999084   # CODATA 2018
alpha_codata = 1.0 / CODATA_ALPHA_INV

pi = np.pi
zeta2 = pi**2 / 6                  # zeta(2) = pi^2/6
zeta3 = 1.2020569031595942         # Apery's constant (15+ digits)
zeta4 = pi**4 / 90
zeta5 = 1.0369277551433699         # zeta(5)
ln2 = np.log(2)

# Zeta derivatives at special points
zeta_prime_0 = -0.5 * np.log(2 * pi)
zeta_prime_neg1 = -0.16542114370045092
zeta_prime_neg2 = -zeta3 / (4 * pi**2)

# Spectral determinants (from Phase 3)
log_det_S1 = 2 * np.log(2 * pi)       # ln(4pi^2)
det_S1 = 4 * pi**2                     # exact

log_det_S3 = zeta3 / (2 * pi**2) + np.log(pi)
det_S3 = pi * np.exp(zeta3 / (2 * pi**2))

log_det_S2 = -4 * zeta_prime_neg1
det_S2 = np.exp(log_det_S2)

# Paper's formula
B_paper = 42
F_paper = zeta2
Delta_paper = 1 / 40
K_paper = pi * (B_paper + F_paper - Delta_paper)
roots_paper = np.roots([1, 0, -K_paper, 1])
alpha_paper = min(r.real for r in roots_paper
                  if abs(r.imag) < 1e-10 and r.real > 0)
alpha_inv_paper = 1.0 / alpha_paper
paper_err = abs(alpha_inv_paper - CODATA_ALPHA_INV) / CODATA_ALPHA_INV

output_dir = os.path.dirname(os.path.abspath(__file__))
results = []


def out(s=""):
    print(s)
    results.append(s)


def solve_cubic_roots(K):
    """Solve alpha^3 - K*alpha + 1 = 0, return sorted real roots."""
    r = np.roots([1, 0, -K, 1])
    return sorted([x.real for x in r if abs(x.imag) < 1e-10])


def alpha_inv_from_K(K):
    """Get alpha^-1 from K via the cubic."""
    pos = [x for x in solve_cubic_roots(K) if x > 0]
    if not pos:
        return None
    return 1.0 / min(pos)


# ================================================================
out("=" * 78)
out("Phase 4: Derivation Audit")
out("=" * 78)
out(f"\nPaper:  K = pi*(42 + zeta(2) - 1/40) = {K_paper:.10f}")
out(f"        alpha^-1 = {alpha_inv_paper:.10f}  (err = {paper_err:.2e})")
out(f"CODATA: alpha^-1 = {CODATA_ALPHA_INV}")


# ################################################################
# INVESTIGATION A: Closing the Determinant Gap
# ################################################################
out("\n" + "#" * 78)
out("INVESTIGATION A: Closing the Determinant Gap")
out("#" * 78)

# ----------------------------------------------------------------
# A.1: High-precision computation of the candidate
# ----------------------------------------------------------------
out("\n" + "=" * 78)
out("A.1: High-precision value of 4pi^2 * exp(zeta(3)/(2pi^2))")
out("=" * 78)

# Exact: det'(S1) * det'(S3) / pi = 4pi^2 * exp(zeta(3)/(2pi^2))
exponent = zeta3 / (2 * pi**2)
B_det = 4 * pi**2 * np.exp(exponent)

out(f"\n  zeta(3)          = {zeta3:.16f}")
out(f"  2*pi^2           = {2*pi**2:.16f}")
out(f"  zeta(3)/(2pi^2)  = {exponent:.16f}")
out(f"  exp(zeta(3)/(2pi^2)) = {np.exp(exponent):.16f}")
out(f"  4pi^2            = {4*pi**2:.16f}")
out(f"  B_det = 4pi^2 * exp(zeta(3)/(2pi^2)) = {B_det:.16f}")
out(f"  B_paper = 42")
out(f"  Shortfall: 42 - B_det = {42 - B_det:.16f}")
out(f"  Relative gap: {abs(42 - B_det)/42:.6e}")

shortfall = 42 - B_det
out(f"\n  The shortfall is delta = {shortfall:.14f}")

# ----------------------------------------------------------------
# A.2: Search for the shortfall among natural quantities
# ----------------------------------------------------------------
out("\n" + "=" * 78)
out("A.2: Is the shortfall a recognizable constant?")
out("=" * 78)

candidates_shortfall = [
    ("1/pi^2",              1/pi**2),
    ("zeta(2)/pi^2 = 1/6",  1/6),
    ("zeta(3)/pi^2",        zeta3/pi**2),
    ("zeta(3)/(2pi^2)",     zeta3/(2*pi**2)),
    ("zeta(3)^2/(4pi^4)",   zeta3**2/(4*pi**4)),
    ("1/(8pi^2)",           1/(8*pi**2)),
    ("zeta(2)/pi^3",        zeta2/pi**3),
    ("1/24",                1/24),
    ("zeta(3)/pi^3",        zeta3/pi**3),
    ("zeta(5)/(2pi^2)",     zeta5/(2*pi**2)),
    ("ln(2)/pi^2",          ln2/pi**2),
    ("1/pi^3",              1/pi**3),
    ("1/23",                1/23),
    ("1/22",                1/22),
    ("exp(-3)",             np.exp(-3)),
    ("Euler_gamma/pi^2",    0.5772156649/pi**2),
    ("1/(2pi^2 - zeta(3))", 1/(2*pi**2 - zeta3)),
    ("zeta(2)^2/pi^3",     zeta2**2/pi**3),
    ("zeta(3)^2/(2pi^4)",  zeta3**2/(2*pi**4)),
    ("2*zeta(3)^2/(pi^4)", 2*zeta3**2/pi**4),
    ("1/(4*pi^2*e)",       1/(4*pi**2*np.e)),
    ("zeta(3)/(4pi^2) [= -zeta'(-2)]", zeta3/(4*pi**2)),
]

out(f"\n  {'Candidate':>35s}  {'Value':>16s}  {'delta - cand':>14s}  {'rel diff':>10s}")
out(f"  {'-'*35}  {'-'*16}  {'-'*14}  {'-'*10}")

best_match = (None, float('inf'))
for name, val in candidates_shortfall:
    diff = shortfall - val
    rdiff = abs(diff / shortfall) if shortfall != 0 else float('inf')
    out(f"  {name:>35s}  {val:16.14f}  {diff:+14.10f}  {rdiff:10.4e}")
    if rdiff < best_match[1]:
        best_match = (name, rdiff, val)

out(f"\n  Best match: {best_match[0]} (rel diff = {best_match[1]:.4e})")
out(f"  None of these match to better than ~1%. The shortfall {shortfall:.6f}")
out(f"  does not appear to be a simple function of pi, zeta values, or ln(2).")

# ----------------------------------------------------------------
# A.3: Alternative direction — use det products instead of 42
# ----------------------------------------------------------------
out("\n" + "=" * 78)
out("A.3: What if the true B is B_det, not 42?")
out("=" * 78)

K_det = pi * (B_det + F_paper - Delta_paper)
ainv_det = alpha_inv_from_K(K_det)
err_det = abs(ainv_det - CODATA_ALPHA_INV) / CODATA_ALPHA_INV

out(f"\n  B_det = det'(S1)*det'(S3)/pi = {B_det:.12f}")
out(f"  K_det = pi*(B_det + zeta(2) - 1/40) = {K_det:.10f}")
out(f"  alpha^-1 (det formula) = {ainv_det:.10f}")
out(f"  alpha^-1 (42 formula)  = {alpha_inv_paper:.10f}")
out(f"  CODATA                 = {CODATA_ALPHA_INV}")
out(f"  Error with B_det: {err_det:.4e}")
out(f"  Error with B=42:  {paper_err:.4e}")

if err_det < paper_err:
    out(f"\n  *** B_det gives BETTER agreement than 42! ***")
    out(f"  Improvement factor: {paper_err/err_det:.2f}x")
else:
    out(f"\n  B=42 gives BETTER agreement than B_det.")
    out(f"  Degradation factor: {err_det/paper_err:.2f}x")

# Also try: K' = det'_1 * det'_3 + pi*zeta(2) - pi/40
# (i.e., B_det*pi replaced by det'_1*det'_3 directly)
K_alt = det_S1 * det_S3 + pi * zeta2 - pi / 40
ainv_alt = alpha_inv_from_K(K_alt)
err_alt = abs(ainv_alt - CODATA_ALPHA_INV) / CODATA_ALPHA_INV if ainv_alt else float('inf')

out(f"\n  Alternative: K = det'_1*det'_3 + pi*zeta(2) - pi/40 = {K_alt:.10f}")
out(f"  alpha^-1 = {ainv_alt:.10f}  (err = {err_alt:.4e})")

# ----------------------------------------------------------------
# A.4: Systematic search over spectral determinant combinations
# ----------------------------------------------------------------
out("\n" + "=" * 78)
out("A.4: Systematic search — det combinations for alpha^-1")
out("=" * 78)

out(f"\n  Available determinants:")
out(f"    det'(S1) = 4pi^2 = {det_S1:.10f}")
out(f"    det'(S2) = {det_S2:.10f}")
out(f"    det'(S3) = {det_S3:.10f}")
out(f"  Log-determinants:")
out(f"    ln det'(S1) = ln(4pi^2) = {log_det_S1:.10f}")
out(f"    ln det'(S2) = {log_det_S2:.10f}")
out(f"    ln det'(S3) = {log_det_S3:.10f}")

# Build large set of candidates for K
det_combos = []

# Products and ratios of determinants, with pi factors
dets = {'d1': det_S1, 'd2': det_S2, 'd3': det_S3,
        'ld1': log_det_S1, 'ld2': log_det_S2, 'ld3': log_det_S3}

# Simple products/ratios
for n1 in ['d1', 'd2', 'd3']:
    for n2 in ['d1', 'd2', 'd3']:
        for op in ['*', '/']:
            v1, v2 = dets[n1], dets[n2]
            val = v1 * v2 if op == '*' else (v1 / v2 if v2 != 0 else 0)
            det_combos.append((f"{n1}{op}{n2}", val))

# With pi factors
for n1 in ['d1', 'd2', 'd3']:
    for pfac, pname in [(pi, 'pi'), (1/pi, '1/pi'), (pi**2, 'pi2'),
                        (1/pi**2, '1/pi2'), (4*pi**2, '4pi2')]:
        det_combos.append((f"{n1}*{pname}", dets[n1] * pfac))

# Triple products
det_combos.append(("d1*d2*d3", det_S1*det_S2*det_S3))
det_combos.append(("d1*d2*d3/pi", det_S1*det_S2*det_S3/pi))
det_combos.append(("d1*d2*d3/pi^2", det_S1*det_S2*det_S3/pi**2))
det_combos.append(("d1*d3/d2", det_S1*det_S3/det_S2))

# Log combinations
for n1 in ['ld1', 'ld2', 'ld3']:
    for n2 in ['ld1', 'ld2', 'ld3']:
        det_combos.append((f"{n1}+{n2}", dets[n1]+dets[n2]))
        det_combos.append((f"{n1}-{n2}", dets[n1]-dets[n2]))
det_combos.append(("ld1+ld2+ld3", log_det_S1+log_det_S2+log_det_S3))

# Now: for each combo, try K = combo + pi*zeta(2) - pi/40
# and K = pi*combo + pi*zeta(2) - pi/40
# and K = combo (raw, see if it directly gives K)
out(f"\n  Searching: K_cand = [combo] or pi*[combo] + fiber + boundary corrections")
out(f"  Target: alpha^-1 as close to CODATA as possible")
out(f"\n  {'Rank':>4}  {'Formula':>50s}  {'K':>12s}  {'alpha^-1':>14s}  {'rel err':>12s}")
out(f"  {'-'*4}  {'-'*50}  {'-'*12}  {'-'*14}  {'-'*12}")

all_results = []

for name, val in det_combos:
    # Try multiple constructions
    constructions = [
        (f"pi*({name}+z2-1/40)",    pi*(val + zeta2 - 1/40)),
        (f"{name}+pi*z2-pi/40",     val + pi*zeta2 - pi/40),
        (f"pi*{name}+pi*z2-pi/40",  pi*val + pi*zeta2 - pi/40),
        (f"{name}",                  val),
    ]
    for cname, K_cand in constructions:
        if K_cand < 3 or K_cand > 1000:
            continue
        ainv = alpha_inv_from_K(K_cand)
        if ainv is not None:
            err = abs(ainv - CODATA_ALPHA_INV) / CODATA_ALPHA_INV
            all_results.append((err, cname, K_cand, ainv))

all_results.sort(key=lambda x: x[0])

for i, (err, cname, K_cand, ainv) in enumerate(all_results[:25]):
    marker = " <-- PAPER" if abs(K_cand - K_paper) / K_paper < 1e-6 else ""
    out(f"  {i+1:4d}  {cname:>50s}  {K_cand:12.6f}  {ainv:14.8f}  {err:12.4e}{marker}")

# Check the best result
if all_results:
    best = all_results[0]
    out(f"\n  Best spectral-determinant formula: err = {best[0]:.4e}")
    out(f"  Compare paper's err = {paper_err:.4e}")
    if best[0] < paper_err:
        out(f"  *** Found a formula with BETTER agreement than the paper! ***")
        out(f"  Formula: {best[1]}")
        out(f"  K = {best[2]:.10f}, alpha^-1 = {best[3]:.10f}")
    else:
        out(f"  No det-based formula beats the paper's formula with B=42.")

# ----------------------------------------------------------------
# A.5: Summary for Investigation A
# ----------------------------------------------------------------
out("\n" + "=" * 78)
out("A.5: Investigation A Summary")
out("=" * 78)
out(f"""
  1. The shortfall 42 - 4pi^2*exp(zeta(3)/(2pi^2)) = {shortfall:.6f}
     is NOT a recognizable combination of standard constants.
     It is not zeta(3)/pi^2, not 1/24, not any simple ratio.

  2. Using B_det instead of 42 gives alpha^-1 = {ainv_det:.8f}
     with error {err_det:.4e} (paper uses 42 with error {paper_err:.4e}).
     {'B_det is BETTER' if err_det < paper_err else 'B=42 is BETTER'}.

  3. No product/ratio of spectral determinants combined with zeta(2)
     and a boundary correction gives a fundamentally better formula
     than the paper's.

  CONCLUSION: 42 = sum((2l+1)*l(l+1), n=1..3) is a DIFFERENT object from
  the spectral determinant product. The 0.1% proximity is suggestive but
  there is no exact bridge. The integer 42 (from the Casimir trace) and
  the transcendental B_det (from spectral zeta) are likely related through
  the spectral geometry of S3, but the connection requires more than
  simple algebraic manipulation.
""")


# ################################################################
# INVESTIGATION B: The Circulant Matrix and the Hopf Bundle
# ################################################################
out("\n" + "#" * 78)
out("INVESTIGATION B: The Circulant Matrix and the Hopf Bundle")
out("#" * 78)

# ----------------------------------------------------------------
# B.1: Physical interpretation of Z3 symmetry
# ----------------------------------------------------------------
out("\n" + "=" * 78)
out("B.1: Physical interpretation of Z3 symmetry")
out("=" * 78)

out("""
  The circulant [[0,b,c],[c,0,b],[b,c,0]] has Z3 symmetry:
  cyclic permutation (1->2->3->1) leaves the matrix invariant.

  In the Hopf bundle, the three components are:
    1: fiber S1 (dim 1)
    2: base S2  (dim 2)
    3: total S3 (dim 3)

  Question: does a natural Z3 action exist on {S1, S2, S3}?
""")

# a) Adjoint action of SU(2) on su(2)
out("  a) SU(2) adjoint action on su(2) ~ R^3:")
out("     The adjoint representation Ad: SU(2) -> SO(3) acts on")
out("     the Lie algebra su(2) spanned by {sigma_1, sigma_2, sigma_3}/2i.")
out("     The Casimir C2 = (3/4)*I in the fundamental rep.")
out("     In the adjoint (j=1): C2 = 2*I.")
out("")
out("     The adjoint Casimir is proportional to identity -> no Z3 structure.")
out("     However, the STRUCTURE CONSTANTS epsilon_{ijk} do define a Z3:")
out("     [sigma_i, sigma_j] = 2i*epsilon_{ijk}*sigma_k")
out("     The cyclic symmetry (1->2->3->1) of epsilon_{ijk} IS a Z3.")
out("")

# Build the matrix from structure constants
out("     The 'structure constant coupling matrix' M_{ij} = sum_k |eps_{ijk}|:")
M_eps = np.array([[0, 1, 1],
                   [1, 0, 1],
                   [1, 1, 0]])
eigs_eps = np.linalg.eigvalsh(M_eps)
out(f"     M = {M_eps.tolist()}")
out(f"     eigenvalues: {sorted(eigs_eps)}")
out(f"     char poly: lambda^3 - 3*lambda - 2 = 0  (since tr=0 is wrong, tr=0 is right")
out(f"     Actually: tr(M) = 0, sum of 2x2 minors = -3, det = 2")
out(f"     Char poly: lambda^3 - 3*lambda - 2 = 0")
out(f"     This is a circulant with b=c=1: lambda^3 - 3*1*1*lambda - (1+1) = 0")
out(f"     Roots: {sorted(eigs_eps)}")
out(f"     Compare our cubic: alpha^3 - K*alpha + 1 = 0 with K~137")
out(f"     The structure is the SAME (depressed cubic from circulant)")
out(f"     but the coupling strength differs by factor K/3 ~ 45.7 vs 1.")

# b) Tangent space decomposition at a point of S3
out("\n  b) Tangent space decomposition under the Hopf map:")
out("     At any point p in S3, the tangent space T_p(S3) decomposes as:")
out("       T_p(S3) = V_p (vertical, 1D) + H_p (horizontal, 2D)")
out("     V_p = tangent to fiber S1 through p")
out("     H_p = horizontal lift of T_{pi(p)}(S2)")
out("")
out("     The metric g on S3 restricts to these subspaces.")
out("     The curvature of the connection is F = (1/2)*vol_{S2}.")
out("")

# Inner product structure
out("     The inner product 'matrix' in the V+H decomposition:")
out("     g = diag(1, 1, 1) in an adapted ONB {e_V, e_H1, e_H2}")
out("     This is trivially identity — the round metric on S3 doesn't")
out("     distinguish vertical from horizontal at the metric level.")
out("     The distinction is in the CONNECTION, not the metric.")

# c) Connection coupling
out("\n  c) Connection-based coupling matrix:")
out("     The Hopf connection A is a 1-form on S3 with curvature F on S2.")
out("     c1(F) = (1/2pi) integral_{S2} F = 1.")
out("")
out("     Natural 'coupling' between bundle components:")
C12 = 1     # Chern number = integral of F over S2
C13_norm = np.sqrt(2 * pi**2)  # norm of connection 1-form over S3
C23_area = 4 * pi  # area of S2 = integral of vol_S2

out(f"     C_12 = c1 = {C12}  (Chern number: fiber-base coupling)")
out(f"     ||A||^2 = integral_S3 |A|^2 = Vol(S3)/4 = pi^2/2 = {pi**2/2:.6f}")
out(f"     ||F||^2 = integral_S2 |F|^2 = 4*pi*(1/4) = pi = {pi:.6f}")
out(f"       (for F = (1/2)*vol_S2 on unit S2)")
out("")
out("     To build a 3x3 coupling matrix we need THREE coupling constants.")
out("     The most natural:")
out(f"       C_{{12}} = c1 = 1  (fiber-base)")
out(f"       C_{{13}} = ||A||_{{L2}} = sqrt(pi^2/2) = {np.sqrt(pi**2/2):.6f}")
out(f"       C_{{23}} = ||F||_{{L2}} = sqrt(pi) = {np.sqrt(pi):.6f}")
out("")

# Check if any product gives K/3
out("     Check: do products of these give K/3 = {:.6f}?".format(K_paper/3))
bc_target = K_paper / 3
prods = [
    ("C12*C13",     C12 * np.sqrt(pi**2/2)),
    ("C12*C23",     C12 * np.sqrt(pi)),
    ("C13*C23",     np.sqrt(pi**2/2) * np.sqrt(pi)),
    ("C13^2",       pi**2/2),
    ("C23^2",       pi),
    ("C12*C13*C23", C12 * np.sqrt(pi**2/2) * np.sqrt(pi)),
    ("pi*C13^2/3",  pi * pi**2 / (2*3)),
    ("det_S1*det_S3/(3*pi)", det_S1*det_S3/(3*pi)),
    ("14*pi",       14*pi),
    ("pi*(B+F-Delta)/3",  bc_target),  # tautological
]
for name, val in prods:
    ratio = val / bc_target
    out(f"       {name:>25s} = {val:10.4f}  (ratio to K/3: {ratio:.4f})")

out(f"\n     No simple combination of Hopf bundle coupling constants")
out(f"     gives bc = K/3 = {bc_target:.4f}.")
out(f"     The factor 14*pi = {14*pi:.4f} is closest (ratio {14*pi/bc_target:.4f}).")

# ----------------------------------------------------------------
# B.2: Why det(M) = -1?
# ----------------------------------------------------------------
out("\n" + "=" * 78)
out("B.2: Why b^3 + c^3 = -1 (i.e., det(M) = -1)?")
out("=" * 78)

out("""
  For the circulant M = [[0,b,c],[c,0,b],[b,c,0]]:
    det(M) = b^3 + c^3

  The condition det(M) = -1 has several possible interpretations:
""")

# a) Orientation
out("  a) ORIENTATION:")
out("     S3 is orientable. The Hopf map pi: S3 -> S2 has degree 1")
out("     (the Hopf invariant). The linking number of two fibers is +1.")
out("     A reversal of S3 orientation changes the sign of det.")
out("     det = -1 could encode: reversed orientation convention")
out("     (total space orientation is OPPOSITE to the product")
out("     orientation fiber x base).")
out("")
out("     Technical: For a principal U(1)-bundle over S2 with c1=1,")
out("     the Euler class of the associated plane bundle is +1.")
out("     The orientation of S3 induced by (fiber, base) has a sign.")
out("     Whether this is +1 or -1 depends on convention.")
out("     det(M) = -1 is CONSISTENT with the Hopf linking number.")

# b) Fermionic sign
out("\n  b) FERMIONIC SIGN:")
out("     The electron has spin 1/2. Under 2pi rotation: psi -> -psi.")
out("     In the Hopf context, the fiber S1 parameterizes the U(1)")
out("     gauge phase. A full traversal of the fiber is a 2pi rotation.")
out("     For fermions, this gives a factor (-1).")
out("     det(M) = -1 = (-1)^1 is consistent with a single fermion.")
out("     For n fermions: det = (-1)^n.")
out("     This is SUGGESTIVE but not a derivation.")

# c) Differential form calculation
out("\n  c) DIFFERENTIAL FORMS:")
out("     Define omega_1 = d(theta) (connection 1-form on S1)")
out("              omega_2 = vol_{S2}/(4pi) (normalized area form)")
out("              omega_3 = vol_{S3}/(2pi^2) (normalized volume form)")
out("")
out("     These are the natural 'volume forms' on each component,")
out("     normalized so integral = 1.")
out("")
out("     The matrix M_{ij} = integral_{S3} omega_i ^ *omega_j")
out("     involves Hodge duals and wedge products.")
out("     This is ill-defined when forms have different degree.")
out("     One cannot wedge a 1-form with a 2-form on S3 in a")
out("     symmetric way to get a 3x3 matrix.")
out("")
out("     FINDING: Differential forms don't naturally give a 3x3")
out("     coupling matrix because the forms live in different degrees.")

# ----------------------------------------------------------------
# B.3: Alternative cubic origins — elliptic curve
# ----------------------------------------------------------------
out("\n" + "=" * 78)
out("B.3: Elliptic curve interpretation")
out("=" * 78)

out("""
  The cubic alpha^3 - K*alpha + 1 = 0 defines the zeros of the
  Weierstrass elliptic curve:
    E: y^2 = x^3 - K*x + 1
""")

# a) j-invariant
K = K_paper
j_inv = 1728 * (4 * K**3) / (4 * K**3 - 27)

out(f"  a) j-invariant:")
out(f"     j = 1728 * 4K^3 / (4K^3 - 27)")
out(f"     4K^3 = {4*K**3:.6f}")
out(f"     4K^3 - 27 = {4*K**3 - 27:.6f}")
out(f"     j = {j_inv:.6f}")
out(f"     j/1728 = {j_inv/1728:.8f}")
out(f"\n     For context:")
out(f"     j = 0: curve with extra Z3 symmetry (y^2 = x^3 + 1)")
out(f"     j = 1728: curve with extra Z4 symmetry (y^2 = x^3 + x)")
out(f"     j = {j_inv:.2f} is very close to 1728 (ratio {j_inv/1728:.6f})")
out(f"     This means our curve is CLOSE to having Z4 symmetry")
out(f"     but doesn't have it exactly. The 'correction' from j=1728")
out(f"     encodes the departure from maximal symmetry.")

# How close to 1728?
out(f"\n     j - 1728 = {j_inv - 1728:.6f}")
out(f"     (j - 1728)/1728 = {(j_inv-1728)/1728:.6e}")
out(f"     This is O(1/K^3) since j ~ 1728 + 27*1728/(4K^3) for large K.")
out(f"     27*1728/(4K^3) = {27*1728/(4*K**3):.6e}")
out(f"     (j-1728)/1728 = {(j_inv-1728)/1728:.6e}")
out(f"     These match: the j-invariant is 1728 + O(K^-3).")

# b) Discriminant
disc = 4 * K**3 - 27
out(f"\n  b) Discriminant of the cubic:")
out(f"     Delta_3 = 4K^3 - 27 = {disc:.6f}")
out(f"     = 4 * {K**3:.6f} - 27")
out(f"     = {disc:.6f}")
out(f"     sqrt(Delta_3) = {np.sqrt(disc):.6f}")
out(f"\n     The discriminant is large and positive (three real roots).")
out(f"     It doesn't factor into recognizable pieces.")
out(f"     Delta_3 / (4pi^3) = {disc/(4*pi**3):.6f}")
out(f"     Delta_3 / K^3 = {disc/K**3:.8f}  (-> 4 for large K)")

# c) Modular interpretation
out(f"\n  c) Modular interpretation:")
out(f"     The j-invariant j ~ 1728 means the curve's period ratio tau")
out(f"     is close to i (the imaginary unit).")
out(f"     j(i) = 1728. For j slightly above 1728, tau is slightly")
out(f"     above i in the upper half-plane.")
out(f"     q = exp(2*pi*i*tau) with tau ~ i gives:")
out(f"     q ~ exp(-2pi) = {np.exp(-2*pi):.8f}")
out(f"\n     This q-value is not obviously related to alpha,")
out(f"     but the connection to modular forms is noted.")

# d) Special values
out(f"\n  d) Conductor analysis:")
out(f"     Our curve E: y^2 = x^3 - {K:.6f}*x + 1 has coefficients")
out(f"     that are NOT rational (K is transcendental).")
out(f"     Elliptic curve theory over Q requires rational coefficients.")
out(f"     The curve is defined over R, not over Q.")
out(f"     Conductor, L-function, BSD conjecture etc. don't apply")
out(f"     in the usual sense.")
out(f"     HOWEVER: if we APPROXIMATE K by a rational, the nearest")
out(f"     integer is {round(K)}, giving E: y^2 = x^3 - 137x + 1.")

# Analyze the rational curve
K_int = 137
disc_int = 4 * K_int**3 - 27
out(f"     For E: y^2 = x^3 - 137x + 1:")
out(f"       Delta = 4*137^3 - 27 = {disc_int}")
out(f"       = {disc_int}")

# Factor the discriminant
d = disc_int
factors = []
for p in range(2, 1000):
    while d % p == 0:
        factors.append(p)
        d //= p
if d > 1:
    factors.append(d)
out(f"       Factoring: {disc_int} = {' * '.join(str(f) for f in factors)}")

out(f"\n     FINDING: The elliptic curve interpretation is interesting")
out(f"     but doesn't provide a DERIVATION. The j-invariant being")
out(f"     close to 1728 is a consequence of K being large, not")
out(f"     a special arithmetic property.")


# ----------------------------------------------------------------
# B.4: Deeper circulant analysis
# ----------------------------------------------------------------
out("\n" + "=" * 78)
out("B.4: What determines b and c in the circulant?")
out("=" * 78)

out(f"\n  From the circulant: bc = K/3, b^3 + c^3 = -1")
out(f"  Let s = b+c, p = bc = K/3 = {K_paper/3:.8f}")
out(f"  Then: s^3 - 3ps + (b^3+c^3) = 0  =>  s^3 - Ks + 1 = 0  (our cubic)")
out(f"  And: b-c = sqrt(s^2 - 4p)")
out(f"\n  The THREE roots s_i of the cubic give three (b_i, c_i) pairs.")
out(f"  Each pair gives a valid circulant matrix.")

bc = K_paper / 3
s_roots = solve_cubic_roots(K_paper)

out(f"\n  {'Root':>6s}  {'s=b+c':>12s}  {'b':>14s}  {'c':>14s}  {'bc':>14s}  {'b^3+c^3':>12s}")
out(f"  {'-'*6}  {'-'*12}  {'-'*14}  {'-'*14}  {'-'*14}  {'-'*12}")

for i, s in enumerate(s_roots):
    disc_s = s**2 - 4 * bc
    if disc_s >= 0:
        b_val = (s + np.sqrt(disc_s)) / 2
        c_val = (s - np.sqrt(disc_s)) / 2
        out(f"  {i+1:>5d}  {s:12.8f}  {b_val:14.8f}  {c_val:14.8f}  "
            f"{b_val*c_val:14.8f}  {b_val**3+c_val**3:12.6f}")
    else:
        b_re = s / 2
        b_im = np.sqrt(-disc_s) / 2
        out(f"  {i+1:>5d}  {s:12.8f}  {b_re:.6f}+{b_im:.6f}i  "
            f"{b_re:.6f}-{b_im:.6f}i  "
            f"(complex)  (complex)")

out(f"\n  The physical root (alpha) corresponds to s ~ {s_roots[0]:.6f}")
out(f"  (the smallest positive root, giving 1/s ~ alpha^-1 ~ 137).")
out(f"\n  b ~ {(s_roots[0] + np.sqrt(s_roots[0]**2 - 4*bc))/2:.8f} (close to alpha)")
out(f"  c ~ {(s_roots[0] - np.sqrt(s_roots[0]**2 - 4*bc))/2:.8f} (large negative)")

out(f"\n  STRUCTURAL INSIGHT:")
out(f"  The condition b^3 + c^3 = -1 is INDEPENDENT of K.")
out(f"  It fixes the determinant, leaving bc = K/3 as the SINGLE free")
out(f"  parameter that encodes the coupling strength.")
out(f"  The question 'where does K come from?' reduces to:")
out(f"    'What determines bc in the circulant?'")
out(f"  And K/3 = bc = pi(42 + zeta(2) - 1/40)/3 = {bc:.6f}.")

# Check: does K/3 = (B_det + zeta(2) - 1/40)*pi/3 give better results?
bc_det = pi * (B_det + zeta2 - Delta_paper) / 3
out(f"\n  With B_det: bc = {bc_det:.6f}")
out(f"  With B=42:  bc = {bc:.6f}")
out(f"  Ratio: {bc/bc_det:.8f}")

# ----------------------------------------------------------------
# B.5: The three Hopf fibrations
# ----------------------------------------------------------------
out("\n" + "=" * 78)
out("B.5: Connection to the three Hopf fibrations?")
out("=" * 78)

out("""
  There are exactly three Hopf fibrations:
    H1: S1 -> S3  -> S2   (fiber dim 1, base dim 2, total dim 3)
    H2: S3 -> S7  -> S4   (fiber dim 3, base dim 4, total dim 7)
    H3: S7 -> S15 -> S8   (fiber dim 7, base dim 8, total dim 15)

  These correspond to the three normed division algebras:
    C (complex), H (quaternions), O (octonions)

  The Z3 of the circulant could correspond to cycling through
  {H1, H2, H3}. But this is speculative — H2 and H3 involve
  different spheres and different physics.
""")

# Dimensions pattern
out("  Dimension pattern:")
for name, dims in [("H1", (1,2,3)), ("H2", (3,4,7)), ("H3", (7,8,15))]:
    f, b, t = dims
    out(f"    {name}: fiber S^{f}, base S^{b}, total S^{t}")
    out(f"         f + b = {f+b} {'= t' if f+b == t else '!= t'}, "
        f"f*b = {f*b}, f+b+t = {f+b+t}")

out(f"\n  If the Z3 acts on {{H1, H2, H3}}, the 'coupling matrix'")
out(f"  would mix quantities from different Hopf fibrations.")
out(f"  This seems unphysical: the paper's formula uses ONLY H1 data.")
out(f"  VERDICT: The Z3 is more likely internal to H1 (the three")
out(f"  components S1, S2, S3) rather than across fibrations.")


# ################################################################
# INVESTIGATION C: The Complete Picture
# ################################################################
out("\n" + "#" * 78)
out("INVESTIGATION C: Assessment of the Derivation Chain")
out("#" * 78)

out("""
  The hypothetical derivation chain:

    S1 -> S3 -> S2  (Hopf bundle)
         |
    Spectral data -> B (Casimir), F (fiber zeta), Delta (boundary)
         |
    K = pi*(B + F - Delta)  (combination rule)
         |
    Z3 circulant coupling -> alpha^3 - K*alpha + 1 = 0
         |
    alpha = smallest positive root
""")

out("  LINK-BY-LINK ASSESSMENT:")
out("  " + "=" * 70)

# Link 1
out("""
  LINK 1: Hopf bundle S1 -> S3 -> S2 as the geometric arena
  STATUS: ESTABLISHED (mathematically standard)

  The Hopf bundle is a well-defined principal U(1)-bundle.
  S3 ~ SU(2), and the bundle structure is unique (c1=1).
  The paper uses it as a natural geometric framework for
  the S3 topology that underlies the graph Laplacian approach.
  This is a CHOICE OF FRAMEWORK, not a derived result.
""")
link1 = "ESTABLISHED"

# Link 2
out("""
  LINK 2: Spectral data -> B, F, Delta
  STATUS: PARTIALLY ESTABLISHED

  B = 42 = sum (2l+1)*l(l+1) for n=1..3
    This is a well-defined spectral sum (degeneracy-weighted Casimir
    trace) on S3. The sum is truncated at n_max=3 via the selection
    principle B/N(3) = 42/14 = 3 = dim(S3). This selection principle
    WORKS but has no derivation from first principles.
    The near-miss with spectral determinants:
      4pi^2*exp(zeta(3)/(2pi^2)) ~ 41.957
    is suggestive but the gap of 0.043 remains unexplained.

  F = zeta(2) = pi^2/6
    The fiber invariant zeta_{S1}(2) = sum 1/k^2 = pi^2/6.
    This is a standard spectral zeta function value.
    WHY zeta(2) specifically (and not zeta(1), zeta(3), etc.)?
    Not derived from the bundle structure.

  Delta = 1/40 = 1/(|lambda_3| * N(2))
    A boundary correction. Phase 1 showed the paper's own derivation
    has an inconsistency (claims g_2=4 but uses N(2)=5).
    The correction is phenomenological: it improves agreement with
    CODATA but its origin is unexplained.
""")
link2 = "PARTIALLY ESTABLISHED"

# Link 3
out("""
  LINK 3: Combination rule K = pi*(B + F - Delta)
  STATUS: NOT ESTABLISHED

  Why this specific combination? The factor pi could be argued as
  natural (stereographic projection relates S3 to R3 with a factor
  of pi). But the ADDITIVE structure B + F - Delta has no derivation.
  Why addition and not multiplication? Why subtraction for Delta?
  No principle selects this combination rule from alternatives like:
    K = pi*B + F - Delta
    K = pi*B*F/Delta
    K = B + pi*F - pi*Delta
  All would give different results. The specific form works because
  it's been TUNED (two parameters: which terms enter, and how).
""")
link3 = "NOT ESTABLISHED"

# Link 4
out("""
  LINK 4: Z3 circulant -> cubic alpha^3 - K*alpha + 1 = 0
  STATUS: PARTIALLY ESTABLISHED

  Phase 3 showed that the cubic IS the characteristic equation of
  a traceless Z3-circulant with bc = K/3 and det = -1.
  This is a STRUCTURAL result: the cubic form is natural for
  Z3-symmetric coupling of three objects.

  But:
  - No derivation connects the bundle components to the circulant.
  - The condition det = -1 lacks a first-principles origin.
  - The Z3 symmetry is suggestive (three bundle components, three
    Pauli matrices, three Hopf fibrations) but not derived.
  - The circulant interpretation EXPLAINS the cubic FORM but does
    not derive K from geometry.
""")
link4 = "PARTIALLY ESTABLISHED"

# Link 5
out("""
  LINK 5: alpha = smallest positive root
  STATUS: ESTABLISHED (given the cubic)

  The depressed cubic x^3 - Kx + 1 = 0 for K > 3 always has three
  real roots. The smallest positive root is uniquely defined and
  gives alpha ~ 1/K for large K (with correction alpha^2/K).
  This step is pure algebra.
""")
link5 = "ESTABLISHED"

# Summary table
out("\n  " + "=" * 70)
out("  DERIVATION CHAIN ASSESSMENT:")
out("  " + "=" * 70)
out(f"""
  Link 1: Hopf bundle as geometric arena      -> ESTABLISHED
  Link 2: Spectral data -> B, F, Delta        -> PARTIALLY ESTABLISHED
  Link 3: Combination rule K = pi*(B+F-Delta)  -> NOT ESTABLISHED
  Link 4: Z3 circulant -> cubic equation       -> PARTIALLY ESTABLISHED
  Link 5: alpha from cubic root                -> ESTABLISHED

  OVERALL: The derivation chain has TWO solid links (1, 5),
  TWO partially established links (2, 4), and ONE missing link (3).

  The CRITICAL GAP is Link 3: the combination rule. Even if we
  accept that B, F, and Delta are the right ingredients (Link 2),
  there is no principle that selects their specific combination.

  The MOST PROMISING direction is the circulant interpretation
  (Link 4). If one could show that the Hopf bundle naturally
  generates a Z3-symmetric coupling matrix with bc = K/3 and
  det = -1, this would simultaneously establish Links 3 and 4.
  But this requires showing that the coupling constant is
  pi*(B+F-Delta)/3, which brings us back to the combination rule.
""")

# ----------------------------------------------------------------
# Final quantitative summary
# ----------------------------------------------------------------
out("\n" + "=" * 78)
out("QUANTITATIVE SUMMARY")
out("=" * 78)

out(f"""
  Paper's result:     alpha^-1 = {alpha_inv_paper:.10f}
  CODATA:             alpha^-1 = {CODATA_ALPHA_INV}
  Relative error:     {paper_err:.2e}
  p-value (Phase 2):  5.2e-9

  Key numbers:
    B = 42 (Casimir trace)
    B_det = 4pi^2 * exp(zeta(3)/(2pi^2)) = {B_det:.8f}
    B - B_det = {42 - B_det:.8f}  (0.1% gap, unexplained)

  Elliptic curve:
    j-invariant = {j_inv:.4f}  (close to 1728, encodes K >> 1)
    Discriminant = {disc:.4f}

  Using B_det instead of 42:
    alpha^-1 = {ainv_det:.10f}  (err = {err_det:.4e})
    {'BETTER' if err_det < paper_err else 'WORSE'} than B=42

  Missing pieces for a complete derivation:
    1. Selection principle for n_max=3 (why 42/14 = 3?)
    2. Why zeta(2) for the fiber invariant?
    3. Why 1/40 for the boundary?
    4. Why additive combination with prefactor pi?
    5. Physical origin of the circulant matrix and det = -1.
""")

out("=" * 78)
out("BOTTOM LINE")
out("=" * 78)
out(f"""
  The formula K = pi*(42 + zeta(2) - 1/40) produces alpha^-1 to
  8.8e-8 precision with zero free parameters (p-value 5.2e-9).

  The cubic structure is natural for Z3-symmetric coupling.
  The spectral ingredients are natural for the Hopf bundle.
  But the COMBINATION RULE is not derived from geometry.

  STATUS: Remarkable numerical coincidence with structural hints,
  but not a derivation. The formula has the character of a
  CONJECTURE that uses the right ingredients in an empirically
  determined combination.

  To elevate this to a derivation, one would need to show:
    "The Hopf bundle Laplacian, when decomposed into fiber-base-total
     components, generates a Z3-circulant coupling matrix whose
     off-diagonal product bc = pi*(42 + zeta(2) - 1/40)/3
     follows from the spectral geometry of the fibration."

  This has NOT been achieved.
""")


# ================================================================
# Save results
# ================================================================
results_path = os.path.join(output_dir, 'phase4_results.txt')
with open(results_path, 'w', encoding='utf-8') as f:
    f.write('\n'.join(results))
    f.write('\n')
out(f"\nResults saved to: {results_path}")
