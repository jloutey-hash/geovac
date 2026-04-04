"""Track BI: Level 3 hyperspherical validation at multiple Z values + ab initio PK parameters."""

import sys
import json
import time
import traceback

sys.path.insert(0, r"c:\Users\jlout\OneDrive\Desktop\Project_Geometric")

import numpy as np

# Approximate exact non-relativistic energies for He-like ions (Drake/NIST)
EXACT_ENERGIES = {
    2: -2.9037,
    3: -7.2799,
    4: -13.6556,
    5: -22.0310,
    6: -32.4062,
    7: -44.7814,
    9: -75.5316,
}

Z_VALUES = [2, 3, 4, 5, 6, 7, 9]

results = {}

# ========== Part 1: Level 3 solver validation ==========
print("=" * 70)
print("PART 1: Level 3 Hyperspherical Solver Validation")
print("=" * 70)

from geovac.hyperspherical_radial import solve_helium

for Z in Z_VALUES:
    print(f"\n--- Z={Z} ---")
    t0 = time.time()
    try:
        # For higher Z, the wavefunction is more compact. Scale R_max down.
        R_max = max(10.0, 30.0 / (Z / 2))
        result = solve_helium(Z=Z, l_max=2, verbose=True, R_max=R_max)
        E = result['energy']
        dt = time.time() - t0
        E_exact = EXACT_ENERGIES[Z]
        err_pct = abs((E - E_exact) / E_exact) * 100
        print(f"  E_computed = {E:.6f} Ha")
        print(f"  E_exact    = {E_exact:.6f} Ha")
        print(f"  Error      = {err_pct:.4f}%")
        print(f"  Time       = {dt:.1f}s")
        results[f"Z{Z}"] = {
            "E_computed": float(E),
            "E_exact": float(E_exact),
            "error_pct": float(err_pct),
            "time_s": float(dt),
            "R_max": R_max,
            "status": "OK",
        }
    except Exception as e:
        dt = time.time() - t0
        print(f"  FAILED: {e}")
        traceback.print_exc()
        results[f"Z{Z}"] = {
            "status": "FAILED",
            "error_msg": str(e),
            "time_s": float(dt),
        }

# ========== Part 2: Ab initio PK parameters ==========
print("\n" + "=" * 70)
print("PART 2: Ab Initio PK Parameters")
print("=" * 70)

from geovac.core_screening import CoreScreening
from geovac.ab_initio_pk import AbInitioPK

PK_Z_VALUES = [3, 4, 5, 6, 7, 9]

for Z in PK_Z_VALUES:
    print(f"\n--- Z={Z} PK ---")
    key = f"Z{Z}"
    t0 = time.time()
    try:
        R_max = max(10.0, 30.0 / (Z / 2))
        cs = CoreScreening(Z=Z, l_max=2, R_max=R_max)
        cs.solve(verbose=True)

        pk = AbInitioPK(cs, n_core=2)
        dt = time.time() - t0

        print(f"  A = {pk.A:.4f} Ha*bohr^2")
        print(f"  B = {pk.B:.4f} 1/bohr^2")
        print(f"  r_core = {pk.r_core:.4f} bohr")
        print(f"  E_core/e = {pk.E_core_per_electron:.4f} Ha")
        print(f"  E_val_est = {pk.E_val_est:.4f} Ha")
        print(f"  Time = {dt:.1f}s")

        # Also get all B candidates
        b_cands = {}
        for method, data in pk.B_candidates.items():
            b_cands[method] = {"r": float(data["r"]), "B": float(data["B"])}

        if key not in results:
            results[key] = {}
        results[key].update({
            "A": float(pk.A),
            "B": float(pk.B),
            "r_core": float(pk.r_core),
            "E_core_per_electron": float(pk.E_core_per_electron),
            "E_val_est": float(pk.E_val_est),
            "B_candidates": b_cands,
            "pk_status": "OK",
        })

    except Exception as e:
        dt = time.time() - t0
        print(f"  PK FAILED: {e}")
        traceback.print_exc()
        if key not in results:
            results[key] = {}
        results[key].update({
            "pk_status": "FAILED",
            "pk_error_msg": str(e),
        })

# ========== Part 3: Z² scaling comparison ==========
print("\n" + "=" * 70)
print("PART 3: Z² Scaling Comparison")
print("=" * 70)

# Anchor to Li (Z=3) if available
if "Z3" in results and "A" in results["Z3"]:
    A_ref = results["Z3"]["A"]
    B_ref = results["Z3"]["B"]
    Z_ref = 3
    print(f"\nAnchor: Z={Z_ref}, A={A_ref:.4f}, B={B_ref:.4f}")

    for Z in PK_Z_VALUES:
        key = f"Z{Z}"
        if key in results and "A" in results[key]:
            A_z2 = A_ref * (Z / Z_ref) ** 2
            B_z2 = B_ref * (Z / Z_ref) ** 2
            A_actual = results[key]["A"]
            B_actual = results[key]["B"]
            A_err = abs((A_actual - A_z2) / A_actual) * 100 if A_actual != 0 else float('inf')
            B_err = abs((B_actual - B_z2) / B_actual) * 100 if B_actual != 0 else float('inf')
            print(f"  Z={Z}: A_ab_initio={A_actual:.4f}, A_z2={A_z2:.4f} ({A_err:.1f}% err)")
            print(f"         B_ab_initio={B_actual:.4f}, B_z2={B_z2:.4f} ({B_err:.1f}% err)")
            results[key]["A_z2_from_Li"] = float(A_z2)
            results[key]["B_z2_from_Li"] = float(B_z2)
            results[key]["A_z2_err_pct"] = float(A_err)
            results[key]["B_z2_err_pct"] = float(B_err)

# Save JSON
outpath = r"c:\Users\jlout\OneDrive\Desktop\Project_Geometric\debug\track_bi\pk_parameters_first_row.json"
with open(outpath, 'w') as f:
    json.dump(results, f, indent=2)
print(f"\nResults saved to {outpath}")
print("\nDone.")
