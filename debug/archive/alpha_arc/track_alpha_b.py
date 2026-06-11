"""
Track alpha-B: Packing pi as Hopf exchange constant.

Subtask 1: Compare sigma_0 = pi d_0^2 / 2 (Paper 0) against the
           Weyl exchange constant for S^2 (Paper 18).

Subtask 2: Investigate K/pi = 43.6199... as a representation-theoretic
           trace on S^3 (Peter-Weyl / Casimir).

Uses mpmath at 50 dps for near-miss scanning.
"""

from mpmath import mp, mpf, pi, zeta, sqrt, fabs, nstr

mp.dps = 50

# -------------------------------------------------------------
# Target value
# -------------------------------------------------------------
B_val = mpf(42)
F_val = pi**2 / 6
Delta_val = mpf(1) / 40

K = pi * (B_val + F_val - Delta_val)
K_over_pi = B_val + F_val - Delta_val

print("K            =", nstr(K, 20))
print("K/pi         =", nstr(K_over_pi, 20))
print("B            =", nstr(B_val, 10))
print("F = pi^2/6   =", nstr(F_val, 20))
print("Delta = 1/40 =", nstr(Delta_val, 10))
print()

# -------------------------------------------------------------
# Subtask 1: sigma_0 vs Weyl S^2
# -------------------------------------------------------------
# Paper 0 Step 1:
#   two points on circle r_1 = d_0, separated by d_0
#   pi d_0^2 shared area  -> sigma_0 = pi d_0^2 / 2  per state
#
# Paper 18 Weyl/Weyl-constant catalog for unit spheres:
#   S^d  eigenvalue counting N(Lambda) ~ (omega_d/(2pi)^d) vol(S^d) Lambda^(d/2)
#   For S^2:  omega_2 = pi (area of unit disk)
#             (2 pi)^2 = 4 pi^2
#             Weyl const = omega_2 / (2 pi)^2 = 1/(4 pi)
#   vol(S^2) = 4 pi
#   So  N(Lambda) ~ (1/(4 pi)) * 4 pi * Lambda = Lambda   (dimensionless)
#
# Therefore the combination that appears when converting counts
# to areas on the S^2 base is:  area-per-state = (area-of-S^2)/N
#
# In Paper 0 the packing plane is R^2 (not S^2); the two-point axis
# tangent to the shell. The "fundamental area per state" sigma_0 has
# dimension [length^2]. For a unit-radius S^2 the analogous per-state
# area is:
#     area-per-state (S^2) = (4 pi)/N    [dimensionless "area" on S^2]
#
# On Paper 0's Euclidean packing plane, for TWO initial states packed on
# the shell area pi d_0^2 (the disk of radius d_0), the area per state is
#     sigma_0 = pi d_0^2 / 2
# i.e. the constant pi/2 times d_0^2.

# Claim to test: is the prefactor pi/2 equal to a clean rational or
# pi-times-rational multiple of the Weyl S^2 constant?
weyl_S2 = mpf(1) / (4*pi)  # omega_2 / (2 pi)^2

ratio_a = (pi/2) / weyl_S2                     # expect 2 pi^2 if they match up to rational
ratio_b = (pi/2) / (4*pi)                      # /area(S^2)-based weyl
ratio_c = (pi/2) / pi                          # 1/2 (trivial sanity)

print("Subtask 1: sigma_0 prefactor relations")
print("  sigma_0 / d_0^2     = pi/2          =", nstr(pi/2, 20))
print("  Weyl S^2 const      = 1/(4 pi)      =", nstr(weyl_S2, 20))
print("  ratio (pi/2)/Weyl   =", nstr(ratio_a, 20), "= 2 pi^2?",
      nstr(2*pi**2, 20))
print("  ratio (pi/2)/vol(S2)=", nstr(ratio_b, 20))
print()

# -------------------------------------------------------------
# Subtask 2: Peter-Weyl Casimir trace on S^3
# -------------------------------------------------------------
# SO(4) = (SU(2) x SU(2))/Z_2. The principal shell n in Paper 2's
# packing labels the irrep (j,j) with j=(n-1)/2; dim V_n = n^2.
#
# Casimir of SU(2) irrep j: j(j+1). Two copies: c_2(V_n) = 2 j(j+1)
# with j = (n-1)/2:
#   c_2(n) = 2 * ((n-1)/2) * ((n+1)/2) = (n^2 - 1)/2.
#
# This matches Paper 2 Eq (14): C_4(n) = (n^2-1)/2 and equals the
# shell-averaged SO(3) Casimir.  Note also the graph eigenvalue
# lambda_n = -(n^2-1), so |lambda_n| = 2 c_2(n).
#
# Truncated traces through n_max = 3:
def casimir4(n):  # SO(4) Casimir eigenvalue on shell n
    return (mpf(n)**2 - 1) / 2

def dimV(n):       # shell degeneracy = n^2
    return mpf(n)**2

n_max = 3

T_c2  = sum(dimV(n)*casimir4(n)          for n in range(1, n_max+1))
T_c22 = sum(dimV(n)*casimir4(n)**2       for n in range(1, n_max+1))
T_lam = sum(dimV(n)*(mpf(n)**2 - 1)      for n in range(1, n_max+1))
T_l3  = sum(dimV(n)*casimir4(n)**3       for n in range(1, n_max+1))

# Paper 2's B = 42 decomposition:
# B = sum_{n=1}^{3} sum_{l=0}^{n-1} (2l+1) l(l+1)
# Inner sum closed form: n^2 (n^2 - 1)/2 = n^2 * c_2(n). Thus
#    B = sum_{n=1}^{3} dim(V_n) * c_2(n) = T_c2.
print("Subtask 2: Peter-Weyl traces on S^3 (truncated n=1..3)")
print("  sum n^2 * c_2(n)        =", nstr(T_c2, 20), "  <-- this is B")
print("  sum n^2 * c_2(n)^2      =", nstr(T_c22, 20))
print("  sum n^2 * |lambda_n|    =", nstr(T_lam, 20))
print("  sum n^2 * c_2(n)^3      =", nstr(T_l3, 20))
print()

# Confirm B = 42 via Paper 2's identity
B_from_casimir = T_c2
print("B (from Casimir trace) =", nstr(B_from_casimir, 10), " (expected 42)")
print()

# -------------------------------------------------------------
# Near-miss scan for K/pi
# -------------------------------------------------------------
target = K_over_pi
print("Searching for K/pi =", nstr(target, 20))

# Candidates: traces plus small corrections
candidates = []

# integer base 42
candidates.append(("42", mpf(42)))
candidates.append(("T_c2 (=B)", T_c2))
candidates.append(("T_c22", T_c22))
candidates.append(("T_lam", T_lam))

# 42 + rational corrections
for num in range(1, 4):
    for den in [6, 8, 10, 20, 40, 60]:
        candidates.append((f"42 + {num}/{den}", mpf(42) + mpf(num)/den))
        candidates.append((f"42 - {num}/{den}", mpf(42) - mpf(num)/den))

# 42 + zeta(2) - rational
for den in [8, 10, 20, 30, 40, 50, 60]:
    candidates.append((f"42 + pi^2/6 - 1/{den}", mpf(42) + F_val - mpf(1)/den))

# Pure zeta combinations
candidates.append(("42 + zeta(2)", mpf(42) + zeta(2)))
candidates.append(("42 + zeta(2) - zeta(3)/8", mpf(42) + zeta(2) - zeta(3)/8))

# 42 + F - 1/(|lambda_3| * N(2))
candidates.append(("42 + F - 1/(8*5)", mpf(42) + F_val - mpf(1)/40))

# Pure traces
# c_2 trace of full shell (not degeneracy weighted)
T_plain_c2 = sum(casimir4(n) for n in range(1, n_max+1))
candidates.append(("sum c_2 only (no deg)", T_plain_c2))
# degeneracy alone
T_deg = sum(dimV(n) for n in range(1, n_max+1))
candidates.append(("sum n^2", T_deg))

# Report sorted by |value - target|
print()
print("Near-miss table (sorted by |value - target|):")
print(f"  {'candidate':40s}  {'value':25s}  {'|diff|':15s}  {'rel':15s}")
rows = []
for name, val in candidates:
    diff = fabs(val - target)
    rel = diff/fabs(target)
    rows.append((rel, name, val, diff))
rows.sort()
for rel, name, val, diff in rows[:15]:
    print(f"  {name:40s}  {nstr(val, 15):25s}  "
          f"{nstr(diff, 8):15s}  {nstr(rel, 6):15s}")
print()

# -------------------------------------------------------------
# Is 42 = sum n^2 (n^2 - 1) / 2 for n=1..3? (B identity)
# -------------------------------------------------------------
b_terms = []
for n in range(1, 4):
    t = mpf(n)**2 * (mpf(n)**2 - 1) / 2
    b_terms.append(int(t))
print(f"B terms n=1..3:  {b_terms}  sum = {sum(b_terms)}  (expected 42)")
# 1*0/2 + 4*3/2 + 9*8/2 = 0 + 6 + 36 = 42 ✓

# -------------------------------------------------------------
# Delta = 1/40 as Casimir correction
# -------------------------------------------------------------
# 1/40 = 1/(8 * 5) = 1/(|lambda_3| * N(2))
# Also: 1/40 = 1/(c_2^max * N(2)*2) etc. Let's enumerate rational
# corrections that give 1/40 from Casimir data:
print()
print("Delta = 1/40 candidate decompositions:")
decomps = [
    ("1/(|lambda_3| * N(2))",  mpf(1)/(8*5)),
    ("1/(c_2(3)^2)",           mpf(1)/casimir4(3)**2),   # 1/16
    ("1/(c_2(3) * c_2(2))",    mpf(1)/(casimir4(3)*casimir4(2))),
    ("1/(dim V_3 * c_2(3))",   mpf(1)/(dimV(3)*casimir4(3))),
    ("1/(2 * dim V_3 * N(2))", mpf(1)/(2*dimV(3)*5)),
    ("1/(T_c2 - 2)",           mpf(1)/(T_c2 - 2)),
]
for name, val in decomps:
    match = "MATCH" if fabs(val - Delta_val) < mpf("1e-12") else "-"
    print(f"  {name:30s} = {nstr(val, 15):20s} {match}")

# -------------------------------------------------------------
# Serialise
# -------------------------------------------------------------
import json
out = {
    "inputs": {
        "B": str(int(B_val)),
        "F_symbolic": "pi^2/6",
        "F_value": nstr(F_val, 25),
        "Delta": "1/40",
        "K_over_pi": nstr(K_over_pi, 25),
        "K": nstr(K, 25),
    },
    "subtask1_packing_vs_weyl": {
        "paper0_sigma0_prefactor": "pi/2",
        "weyl_S2_constant": nstr(weyl_S2, 25),
        "vol_S2": nstr(4*pi, 25),
        "ratio_sigma0_prefactor_to_weyl": nstr(ratio_a, 25),
        "ratio_equals_2_pi_squared": bool(fabs(ratio_a - 2*pi**2) < mpf("1e-30")),
        "verdict": (
            "sigma_0 and the S^2 Weyl constant are BOTH pi-type geometric "
            "exchanges, but they are not the same object. sigma_0 sits on "
            "a EUCLIDEAN packing plane (R^2) and has units [length^2], "
            "whereas the Weyl constant 1/(4 pi) is a dimensionless count "
            "density on the unit S^2 base of the Hopf bundle. The outer pi "
            "in K is a pi^1 of type 'area-per-state on a 2-manifold', which "
            "matches the class but NOT the specific Weyl constant of S^2."
        ),
    },
    "subtask2_casimir_trace": {
        "T_c2": nstr(T_c2, 25),
        "T_c2_squared": nstr(T_c22, 25),
        "T_lambda": nstr(T_lam, 25),
        "T_c2_cubed": nstr(T_l3, 25),
        "B_identity": "B = sum_{n=1}^{3} n^2 * (n^2-1)/2 = 42",
        "B_from_casimir_trace": int(T_c2),
        "closest_nearmiss": [
            {
                "name": name,
                "value": nstr(val, 25),
                "abs_diff": nstr(diff, 15),
                "rel_diff": nstr(rel, 15),
            }
            for rel, name, val, diff in rows[:5]
        ],
    },
    "delta_decomposition_candidates": [
        {"form": name, "value": nstr(val, 25),
         "matches_1_over_40": bool(fabs(val - Delta_val) < mpf("1e-12"))}
        for name, val in decomps
    ],
}
import os
outdir = "c:/Users/jlout/OneDrive/Desktop/Project_Geometric/debug/data/track_alpha_phase4b"
os.makedirs(outdir, exist_ok=True)
with open(os.path.join(outdir, "track_b_packing_pi.json"), "w") as f:
    json.dump(out, f, indent=2)

print()
print("wrote", os.path.join(outdir, "track_b_packing_pi.json"))
