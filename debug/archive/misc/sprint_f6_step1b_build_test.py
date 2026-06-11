"""Sprint F6 Step 1b — empirical build-time test at NaH max_n=4.

Confirm actual build time matches Step 1 estimate before committing to full
FCI. Builds ONE Hamiltonian at R=3.5 with the full F3 stack (W1c + multi-zeta
+ cross-block h1) and reports timing breakdown.

Decision:
  - If build < 30 min: confident GO for Step 2 at R=3.5 and R=10
  - If build 30-90 min: GO with single-point only at R=3.5 for the binding test
  - If build > 90 min: STOP, reconsider compute budget
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.stdout.reconfigure(line_buffering=True)

from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.molecular_spec import nah_spec


def main():
    out_path = Path("debug/data/sprint_f6_step1b_build_test.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("Sprint F6 Step 1b - empirical build test at NaH max_n=4")
    print("=" * 70)

    spec = nah_spec(max_n=4)
    R = 3.5

    print(f"\nBuilding NaH max_n=4 at R={R} with full F3 stack...")
    print(f"  screened_cross_center=True, multi_zeta_basis=True,")
    print(f"  cross_block_h1=True")
    print(f"  (M=60, Q=120, expected fci_dim=3600 for 2e singlet)")

    t0 = time.perf_counter()
    try:
        ham = build_balanced_hamiltonian(
            spec, R=R, verbose=True,
            screened_cross_center=True,
            multi_zeta_basis=True,
            cross_block_h1=True,
        )
        elapsed = time.perf_counter() - t0
        success = True
        M = ham["M"]
        h1_shape = ham["h1"].shape
        eri_shape = ham["eri"].shape
        nuc_rep = ham["nuclear_repulsion"]
        xblock_info = ham.get("cross_block_h1_info", {})
        h1_norm = float(np.linalg.norm(ham["h1"]))
        eri_norm = float(np.linalg.norm(ham["eri"]))
        max_h1 = float(np.max(np.abs(ham["h1"])))
        max_eri = float(np.max(np.abs(ham["eri"])))
        error = None
    except (NotImplementedError, Exception) as e:
        elapsed = time.perf_counter() - t0
        success = False
        M = None
        h1_shape = None
        eri_shape = None
        nuc_rep = None
        xblock_info = None
        h1_norm = None
        eri_norm = None
        max_h1 = None
        max_eri = None
        error = f"{type(e).__name__}: {e}"
        print(f"\nERROR during build: {error}")

    print(f"\nBuild complete (elapsed: {elapsed:.1f}s = {elapsed/60:.2f}min)")
    if success:
        print(f"  M = {M}")
        print(f"  h1 shape: {h1_shape}")
        print(f"  eri shape: {eri_shape}")
        print(f"  nuclear repulsion: {nuc_rep:.4f} Ha")
        print(f"  ||h1||_F = {h1_norm:.4f}, max|h1| = {max_h1:.4f}")
        print(f"  ||eri||_F = {eri_norm:.4f}, max|eri| = {max_eri:.4f}")
        print(f"  cross-block h1 info: {xblock_info}")
    else:
        print(f"  BUILD FAILED: {error}")

    if elapsed < 30 * 60:
        gate = "GO_FULL_2POINT"
        reason = "Build < 30min, safe for both R=3.5 and R=10"
    elif elapsed < 90 * 60:
        gate = "GO_SINGLE_POINT"
        reason = "Build 30-90min, recommend single R=3.5 first"
    else:
        gate = "STOP_OVERBUDGET"
        reason = "Build > 90min, reconsider compute budget"

    print(f"\n{'='*70}")
    print(f"GATE: {gate} - {reason}")
    print(f"{'='*70}")

    out = {
        "sprint": "F6 Step 1b - empirical build test",
        "date": "2026-05-23",
        "R": R,
        "success": success,
        "build_time_s": elapsed,
        "build_time_min": elapsed / 60,
        "M": M,
        "h1_shape": list(h1_shape) if h1_shape else None,
        "eri_shape": list(eri_shape) if eri_shape else None,
        "nuclear_repulsion": nuc_rep,
        "h1_frobenius_norm": h1_norm,
        "eri_frobenius_norm": eri_norm,
        "max_abs_h1": max_h1,
        "max_abs_eri": max_eri,
        "cross_block_h1_info": xblock_info,
        "error": error,
        "gate": gate,
        "gate_reason": reason,
    }
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"\nWrote: {out_path}")


if __name__ == "__main__":
    main()
