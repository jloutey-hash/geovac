"""Extended convergence + diagnostics for the He 2^1P -> 1^1S oscillator strength.

Run n_max=6,7 and decompose the oscillator strength into the omega and
|<z>|^2 pieces to see where the error is coming from.

Also: compute the SINGLE-CONFIGURATION baseline (no CI) to see if the
problem is multi-electron correlation or single-particle.

Key hypothesis from v1: f converges slowly. We need to know whether:
 (a) The transition energy omega converges with n_max (it should, like
     ground-state energy CI does)
 (b) The dipole matrix element converges (it should if Sturmian closure works)
"""

import sys
import time
import json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
from calc_track_he_oscillator_v1 import (
    compute_oscillator_strength,
    F_DRAKE_2P_1S,
)

PROJECT_ROOT = Path(__file__).resolve().parent.parent

# Single-config sanity at Z=2 from sympy:
# E(1s)^2 HF = 2 * (-Z^2/2) + F0(1s,1s) = -4 + 5/8 = -27/8 = -3.375
# E(1s,2p) HF = -Z^2/2(1/1 + 1/4) + V_ee for (1s)(2p) singlet
# But our basis is k_orb=Z=2 so 1s/2s/2p use the same exponent k=2
# This gives a CI calculation, not pure HF.

print("Extended convergence study (Path A, k_orb=Z=2):")
print("=" * 80)
print(f"{'n_max':>5} {'configs':>10} {'omega(Ha)':>10} {'<z>':>10} {'<z>^2':>10} "
      f"{'f_length':>10} {'err%':>9} {'wall(s)':>9}")
print("-" * 80)

E_1S_drake = -2.903724377
E_2P_drake = -2.123843087   # NIST 2^1P at -2.12384
# ω_drake = 0.77988 Ha
omega_drake = E_2P_drake - E_1S_drake

results = []
for n_max in [2, 3, 4, 5, 6, 7]:
    t0 = time.time()
    try:
        r = compute_oscillator_strength(n_max=n_max, Z=2, k_orb=2.0, verbose=False)
        elapsed = time.time() - t0
        r["wall_seconds"] = elapsed
        err = (r["f_length"] - F_DRAKE_2P_1S) / F_DRAKE_2P_1S * 100.0
        omega_err = (r["omega_Ha"] - omega_drake) / omega_drake * 100.0
        E_1S_err = abs(r["E_1S_Ha"] - E_1S_drake) / abs(E_1S_drake) * 100.0
        E_2P_err = abs(r["E_2P_Ha"] - E_2P_drake) / abs(E_2P_drake) * 100.0
        n_configs_total = r["n_configs_1S"] + r["n_configs_2P"]
        print(f"{n_max:>5} {n_configs_total:>10} {r['omega_Ha']:>10.4f} "
              f"{r['dipole_z_au']:>10.4f} {r['dipole_squared_au']:>10.4f} "
              f"{r['f_length']:>10.4f} {err:>9.2f} {elapsed:>9.1f}")
        results.append({
            **r,
            "error_pct": err,
            "omega_error_pct": omega_err,
            "E_1S_error_pct": E_1S_err,
            "E_2P_error_pct": E_2P_err,
        })
    except Exception as e:
        print(f"{n_max:>5} FAILED: {type(e).__name__}: {e}")

print("-" * 80)
print(f"{'EXACT':>5} {'':>10} {omega_drake:>10.4f} {'':>10} {'':>10} "
      f"{F_DRAKE_2P_1S:>10.4f} {0.00:>9.2f}")
print()
print(f"E(1^1S) Drake exact = {E_1S_drake:.6f} Ha")
print(f"E(2^1P) Drake exact = {E_2P_drake:.6f} Ha")
print(f"omega exact = {omega_drake:.6f} Ha = {omega_drake*27.21138:.4f} eV")

# Save to JSON
out = {
    "results": results,
    "drake_reference": {
        "E_1S_Ha": E_1S_drake,
        "E_2P_Ha": E_2P_drake,
        "omega_Ha": omega_drake,
        "f_2P_1S": F_DRAKE_2P_1S,
    },
}
out_file = PROJECT_ROOT / "debug" / "data" / "he_oscillator_extended.json"
with open(out_file, "w") as fp:
    json.dump(out, fp, indent=2, default=float)
print(f"\nSaved: {out_file}")
