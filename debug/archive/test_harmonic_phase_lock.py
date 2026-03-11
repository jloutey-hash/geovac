"""
Harmonic Phase Locking Test — Pure theoretical diagnostic.

Hypothesis: The critical points (maxima, minima, inflections) of the SO(4)
Wigner D-matrix elements D^n_{(l',m'),(l,m)}(γ) map to physically meaningful
bond lengths R* via the gamma->R inverse stereographic map.

If any critical point maps to R ≈ 3.015 bohr (LiH R_eq) under a physically
motivated p₀, that is a closed-form geometric prediction of R_eq from SO(4)
representation theory alone.

Date: 2026-03-11
Status: Diagnostic (no production code changes)
"""

import sys
import os
import numpy as np

# Add project root to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from geovac.wigner_so4 import wigner_D_so4, bond_angle


# ============================================================================
#  Configuration
# ============================================================================

N_POINTS = 10_000
GAMMA_MIN = 0.01
GAMMA_MAX = np.pi - 0.01

R_EQ = 3.015  # LiH experimental equilibrium distance (bohr)

P0_VALUES = {
    0.65: "Z_eff/n = 1.3/2 for Li 2s",
    0.92: "Li 2s: sqrt(2*0.42)",
    1.00: "baseline atomic (H 1s)",
    1.92: "sqrt(2*|E_Li_2s_exact|)",
    2.00: "converged dual-p0 (v0.9.34)",
}

# D-matrix elements to probe: (label, n, lp, mp, l, m)
ELEMENTS = [
    ("D1_00_00",  1, 0, 0, 0, 0),   # n=1: 1s→1s
    ("D2_00_00",  2, 0, 0, 0, 0),   # n=2: 2s→2s
    ("D2_10_10",  2, 1, 0, 1, 0),   # n=2: 2p0→2p0
    ("D2_00_10",  2, 0, 0, 1, 0),   # n=2: 2s→2p0 (off-diagonal)
    ("D3_00_00",  3, 0, 0, 0, 0),   # n=3: 3s→3s
]

NEAR_HIT = (2.8, 3.2)
STRONG_HIT = (2.9, 3.1)


# ============================================================================
#  Helper: γ → R inverse map
# ============================================================================

def gamma_to_R(gamma: float, p0: float) -> float:
    """Inverse stereographic map: gamma, p0 -> R = 1/p_R where p_R = p₀|tan(γ/2)|."""
    tan_half = np.abs(np.tan(gamma / 2.0))
    if tan_half < 1e-15:
        return np.inf
    p_R = p0 * tan_half
    return 1.0 / p_R


# ============================================================================
#  Step 1: Dense angular sweep — compute D-matrix elements
# ============================================================================

print("=" * 80)
print("HARMONIC PHASE LOCKING TEST")
print("=" * 80)
print(f"\nGrid: {N_POINTS} points, gamma in [{GAMMA_MIN:.4f}, {GAMMA_MAX:.4f}]")
print(f"R_eq(LiH) = {R_EQ} bohr\n")

gamma_grid = np.linspace(GAMMA_MIN, GAMMA_MAX, N_POINTS)

# Compute D-matrix element values across the grid
D_values = {}
for label, n, lp, mp, l, m in ELEMENTS:
    print(f"  Computing {label} (n={n}, l'={lp}, m'={mp}, l={l}, m={m}) ...")
    vals = np.array([wigner_D_so4(n, lp, mp, l, m, g) for g in gamma_grid])
    D_values[label] = vals
    print(f"    range: [{vals.min():.6f}, {vals.max():.6f}]")

print()


# ============================================================================
#  Step 2: Find critical points
# ============================================================================

def find_critical_points(gamma_grid: np.ndarray, values: np.ndarray):
    """Find maxima, minima, and steepest inflection points."""
    dv = np.gradient(values, gamma_grid)
    d2v = np.gradient(dv, gamma_grid)

    # Zero-crossings of dv (maxima and minima)
    sign_changes = np.where(np.diff(np.sign(dv)))[0]

    maxima = []
    minima = []
    for idx in sign_changes:
        # Interpolate zero-crossing
        g0, g1 = gamma_grid[idx], gamma_grid[idx + 1]
        dv0, dv1 = dv[idx], dv[idx + 1]
        if abs(dv1 - dv0) < 1e-30:
            continue
        g_star = g0 - dv0 * (g1 - g0) / (dv1 - dv0)
        d2v_star = d2v[idx]  # approximate

        if d2v_star < 0:
            v_star = np.interp(g_star, gamma_grid, values)
            maxima.append((g_star, v_star))
        elif d2v_star > 0:
            v_star = np.interp(g_star, gamma_grid, values)
            minima.append((g_star, v_star))

    # Zero-crossings of d2v (inflection points) — find steepest
    sign_changes_d2 = np.where(np.diff(np.sign(d2v)))[0]
    inflections = []
    for idx in sign_changes_d2:
        g0, g1 = gamma_grid[idx], gamma_grid[idx + 1]
        d2v0, d2v1 = d2v[idx], d2v[idx + 1]
        if abs(d2v1 - d2v0) < 1e-30:
            continue
        g_star = g0 - d2v0 * (g1 - g0) / (d2v1 - d2v0)
        dv_star = np.interp(g_star, gamma_grid, dv)
        inflections.append((g_star, abs(dv_star)))

    # Sort inflections by steepness
    inflections.sort(key=lambda x: -x[1])

    return maxima, minima, inflections


critical_points = {}
for label, n, lp, mp, l, m in ELEMENTS:
    vals = D_values[label]
    val_range = vals.max() - vals.min()
    if val_range < 1e-6:
        print(f"  {label}: CONSTANT (range={val_range:.2e}), skipping critical points")
        critical_points[label] = {'maxima': [], 'minima': [], 'inflections': [],
                                  'constant': True}
        continue
    maxima, minima, inflections = find_critical_points(gamma_grid, vals)
    critical_points[label] = {
        'maxima': maxima,
        'minima': minima,
        'inflections': inflections,
        'constant': False,
    }
    print(f"  {label}: {len(maxima)} max, {len(minima)} min, {len(inflections)} inflection")


# ============================================================================
#  Step 3: Map to physical R
# ============================================================================

p0_list = sorted(P0_VALUES.keys())


def format_R(R_val):
    """Format R value, handling inf."""
    if R_val > 1000:
        return "   inf  "
    return f"{R_val:8.3f}"


def classify_hit(R_val):
    """Classify as strong hit, near hit, or miss."""
    if STRONG_HIT[0] <= R_val <= STRONG_HIT[1]:
        return "STRONG"
    elif NEAR_HIT[0] <= R_val <= NEAR_HIT[1]:
        return "NEAR"
    return ""


# Collect all results for output
output_lines = []


def out(line=""):
    print(line)
    output_lines.append(line)


out("\n" + "=" * 100)
out("D-matrix critical points and R* mapping")
out("=" * 100)

header = f"{'Element':<16} {'Type':<12} {'gamma*':>8}  {'D-value':>10}"
for p0 in p0_list:
    header += f"  {'p0='+str(p0):>10}"
out(header)
out("-" * len(header))

near_hits = []
strong_hits = []

for label, n, lp, mp, l, m in ELEMENTS:
    cp = critical_points[label]
    if cp.get('constant', False):
        out(f"{label:<16} CONSTANT (trivial, no critical points)")
        out("")
        continue

    for cp_type, points in [('maximum', cp['maxima']),
                            ('minimum', cp['minima']),
                            ('inflection', cp['inflections'][:3])]:  # top 3 inflections
        for i, pt in enumerate(points):
            if cp_type == 'inflection':
                g_star, steepness = pt
                v_star = np.interp(g_star, gamma_grid, D_values[label])
            else:
                g_star, v_star = pt

            line = f"{label:<16} {cp_type:<12} {g_star:8.4f}  {v_star:10.6f}"

            for p0 in p0_list:
                R_star = gamma_to_R(g_star, p0)
                line += f"  {format_R(R_star)}"

                hit = classify_hit(R_star)
                if hit == "STRONG":
                    strong_hits.append((label, cp_type, g_star, p0, R_star))
                elif hit == "NEAR":
                    near_hits.append((label, cp_type, g_star, p0, R_star))

            out(line)

    out("")


# ============================================================================
#  Step 4: D-matrix values at R_eq = 3.015 bohr
# ============================================================================

out("\n" + "=" * 100)
out(f"D-matrix values at R_eq = {R_EQ} bohr")
out("=" * 100)

header2 = f"{'Element':<20}"
for p0 in p0_list:
    header2 += f"  {'p0='+str(p0):>12}"
out(header2)

sub_header = f"{'':20}"
for p0 in p0_list:
    gamma_eq = bond_angle(R_EQ, p0)
    sub_header += f"  {'g='+f'{gamma_eq:.4f}':>12}"
out(sub_header)
out("-" * len(header2))

for label, n, lp, mp, l, m in ELEMENTS:
    line = f"{label:<20}"
    for p0 in p0_list:
        gamma_eq = bond_angle(R_EQ, p0)
        D_val = wigner_D_so4(n, lp, mp, l, m, gamma_eq)
        line += f"  {D_val:12.6f}"
    out(line)

out("")

# Also report the derivative at R_eq to see if we're near a critical point
out(f"{'dD/dg at R_eq':<20}")
out("-" * len(header2))
for label, n, lp, mp, l, m in ELEMENTS:
    line = f"{label:<20}"
    for p0 in p0_list:
        gamma_eq = bond_angle(R_EQ, p0)
        # Numerical derivative
        dg = 0.001
        D_plus = wigner_D_so4(n, lp, mp, l, m, gamma_eq + dg)
        D_minus = wigner_D_so4(n, lp, mp, l, m, gamma_eq - dg)
        dD = (D_plus - D_minus) / (2 * dg)
        line += f"  {dD:12.6f}"
    out(line)

out("")


# ============================================================================
#  Summary: Hits
# ============================================================================

out("\n" + "=" * 100)
out("RESULTS SUMMARY")
out("=" * 100)

if strong_hits:
    out(f"\nStrong hits (R* in [{STRONG_HIT[0]}, {STRONG_HIT[1]}]):")
    for label, cp_type, g_star, p0, R_star in strong_hits:
        out(f"  {label:16s} {cp_type:12s}  g*={g_star:.4f}  p0={p0}  "
            f"R*={R_star:.4f} bohr  ({P0_VALUES[p0]})")
else:
    out(f"\nStrong hits (R* in [{STRONG_HIT[0]}, {STRONG_HIT[1]}]): NONE")

if near_hits:
    out(f"\nNear hits (R* in [{NEAR_HIT[0]}, {NEAR_HIT[1]}]):")
    for label, cp_type, g_star, p0, R_star in near_hits:
        out(f"  {label:16s} {cp_type:12s}  g*={g_star:.4f}  p0={p0}  "
            f"R*={R_star:.4f} bohr  ({P0_VALUES[p0]})")
else:
    out(f"\nNear hits (R* in [{NEAR_HIT[0]}, {NEAR_HIT[1]}]): NONE")


# ============================================================================
#  p₀ values and their physical motivation
# ============================================================================

out("\n" + "-" * 60)
out("p0 values tested:")
for p0 in p0_list:
    out(f"  p0 = {p0:<5.2f}  {P0_VALUES[p0]}")

out("\n" + "=" * 100)


# ============================================================================
#  Save to file
# ============================================================================

os.makedirs(os.path.join(os.path.dirname(__file__), 'data'), exist_ok=True)
outpath = os.path.join(os.path.dirname(__file__), 'data',
                       'harmonic_phase_lock_analysis.txt')
with open(outpath, 'w') as f:
    f.write('\n'.join(output_lines))
    f.write('\n')

print(f"\nOutput saved to: {outpath}")
