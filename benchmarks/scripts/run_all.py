#!/usr/bin/env python3
"""
GeoVac 3-Electron Extension — Run All Deliverables
====================================================

Entry point that executes all three deliverables for the 3-electron
sparse dynamics paper and writes figures to ./figures/.

Usage:
    python run_all.py                  # full run
    python run_all.py --quick          # reduced steps for testing
    python run_all.py --nmax 5         # override max_n for scaling benchmark
    python run_all.py --gaps-only      # run only measurement gap diagnostics

Deliverable 1 — LatticeIndex (3-electron FCI for Lithium)
  Builds the topology index for Li (Z=3, n_electrons=3, max_n=4),
  assembles the sparse FCI Hamiltonian, and reports the ground state.

Deliverable 2 — Benchmarks
  Test A: Unitarity — 10,000 CN steps on Li single-particle lattice.
  Test B: Scaling   — nmax sweep, two-panel log-log figure.

Deliverable 3 — Langevin AIMD
  LiH at T=300K for 600 steps. Reports energy conservation tolerance.

Gap Diagnostics (--gaps-only or appended to full run):
  Gap 1: NVE control — VelocityVerlet on LiH without thermostat (600 steps)
  Gap 2: Li energy audit — GeoVac vs exact at max_n=4, 6, 8
  Gap 3: Two-panel scaling figure — sparsity mask audit (NNZ vs max_n)

All known systematic errors are logged as warnings, not suppressed:
  - V_ee Mulliken ~12% overestimation (larger for same-orbital pairs)
  - ~17% mean-field correlation error in AIMD energy surface

Figures written to ./figures/:
  unitarity_norm_vs_step.png   (Test A)
  scaling_log_log.png          (Test B — primary publishable figure)
  aimd_lih_langevin.png        (Deliverable 3)

Author: GeoVac Development Team
Date: February 2026
"""

import argparse
import os
import sys
import warnings
from typing import Optional

import numpy as np

# Ensure geovac is importable from repo root
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def run_deliverable_1(max_n: int = 4) -> None:
    """
    Deliverable 1: Build LatticeIndex for Li (3-electron FCI) and report
    ground state.
    """
    from geovac.lattice_index import LatticeIndex

    print("\n" + "=" * 65)
    print("DELIVERABLE 1 — LatticeIndex: Li 3-Electron FCI")
    print("=" * 65)
    print(f"\n  Z=3, n_electrons=3, max_n={max_n}")

    # Warnings are intentional — do not suppress
    idx = LatticeIndex(n_electrons=3, max_n=max_n, nuclear_charge=3)

    stats = idx.assembly_stats()
    print(f"\n  Topology index built:")
    print(f"    Spatial states:    {stats['n_spatial']}")
    print(f"    Spin-orbitals:     {stats['n_sp']}")
    print(f"    Slater dets:       {stats['n_sd']:,}")
    print(f"    H1 nnz:            {stats['h1_nnz']:,}")
    print(f"    H1 off-diag edges: {stats['h1_offdiag_edges']:,}")

    print(f"\n  Conformal weights (first 5 spatial states):")
    for i, (n, l, m) in enumerate(idx.lattice.states[:5]):
        print(f"    ({n},{l},{m}): Omega = {idx.conformal_weights[i]:.6f}")

    # Assemble Hamiltonian and solve
    print(f"\n  Assembling FCI Hamiltonian...")
    H = idx.assemble_hamiltonian()
    print(f"  H shape: {H.shape}, nnz: {H.nnz:,}, "
          f"sparsity: {1 - H.nnz / H.shape[0]**2:.4f}")

    print(f"\n  Solving for ground state (sparse eigsh)...")
    eigvals, eigvecs = idx.compute_ground_state(n_states=3)

    print(f"\n  Three lowest eigenvalues:")
    for i, E in enumerate(eigvals):
        print(f"    E_{i} = {E:.6f} Ha")

    print(f"\n  Reference (exact Li FCI): E_0 = -7.478 Ha")
    print(f"  Mulliken V_ee overestimation shifts E_0 significantly above exact.")
    print(f"  This is the known systematic error documented in LatticeIndex.")

    # Report dominant SD
    psi0 = eigvecs[:, 0]
    dominant_idx = int(np.argmax(np.abs(psi0)))
    dominant_sd = idx.sd_basis[dominant_idx]
    dominant_states = [idx.sp_states[j] for j in dominant_sd]
    print(f"\n  Dominant SD in ground state (weight {psi0[dominant_idx]**2:.3f}):")
    for s in dominant_states:
        print(f"    ({s[0]},{s[1]},{s[2]},{'up' if s[3]==0 else 'dn'})")


def run_deliverable_2(
    quick: bool = False,
    nmax_max: int = 8,
    output_dir: str = "./figures",
) -> dict:
    """Deliverable 2: Run both benchmark tests. Returns dict with unitarity result."""
    from geovac.benchmark import run_unitarity_test, run_scaling_benchmark

    results = {}

    # Test A — Unitarity
    n_steps = 1_000 if quick else 10_000
    res_a = run_unitarity_test(
        max_n=5, n_steps=n_steps, dt=0.01, output_dir=output_dir
    )
    results["unitarity_max_dev"] = res_a.get("max_deviation")
    if res_a["passed"]:
        print(f"\n  [PASS] Unitarity: max |<psi|psi>-1| = {res_a['max_deviation']:.3e} < 1e-12")
    else:
        print(f"\n  [FAIL] Unitarity: max |<psi|psi>-1| = {res_a['max_deviation']:.3e}")

    # Test B — Scaling
    nmax_range = list(range(3, min(nmax_max, 9) + 1))
    res_b = run_scaling_benchmark(nmax_range=nmax_range, output_dir=output_dir)
    results["scaling"] = res_b
    if res_b["alpha"] is not None:
        print(f"\n  [{'PASS' if res_b['passed'] else 'FAIL'}] "
              f"Scaling: alpha = {res_b['alpha']:.3f} "
              f"({'< 2.0 sub-quadratic' if res_b['passed'] else '>= 2.0'})")

    return results


def run_deliverable_3(
    quick: bool = False,
    output_dir: str = "./figures",
) -> None:
    """Deliverable 3: Langevin AIMD for LiH at T=300K."""
    from geovac.aimd import run_lih_aimd

    n_steps = 60 if quick else 600
    result = run_lih_aimd(
        max_n=3,
        n_steps=n_steps,
        T_kelvin=300.0,
        gamma=0.01,
        dt=1.0,
        R_init=3.0,
        output_dir=output_dir,
        seed=42,
    )

    print(f"\n  AIMD energy drift: {result['energy_drift_pct']:.4f}% "
          f"(std/|mean| over {n_steps} steps)")
    print(f"  Note: Langevin NVT — fluctuations reflect thermal bath coupling,")
    print(f"        not numerical integration error.")


def _print_pub_summary(
    gap_results: dict,
    unitarity_max_dev: Optional[float] = None,
) -> None:
    """Print the fixed-format publication readiness summary block."""
    print("\n\n" + "=" * 38)
    print("=== PUBLICATION READINESS SUMMARY ===")
    print("=" * 38)

    # [NVE CONTROL]
    nve = gap_results.get("nve", {})
    if isinstance(nve, dict) and "error" in nve:
        print(f"[NVE CONTROL]     ERROR: {nve['error']}")
    elif isinstance(nve, dict) and "drift_pct" in nve:
        print(f"[NVE CONTROL]     energy drift: {nve['drift_pct']:.6f}%")
    else:
        print("[NVE CONTROL]     not run")

    # [LI ENERGY]
    audit = gap_results.get("audit")
    if isinstance(audit, dict) and "error" in audit:
        print(f"[LI ENERGY]       ERROR: {audit['error']}")
    elif isinstance(audit, list):
        valid_rows = [r for r in audit if not r.get("skipped", True)]
        if valid_rows:
            best = min(valid_rows, key=lambda r: r["pct_err"])
            print(f"[LI ENERGY]       error at max_n={best['max_n']}: "
                  f"{best['pct_err']:.2f}% vs exact")
        else:
            print("[LI ENERGY]       all max_n runs skipped (resource limit)")
    else:
        print("[LI ENERGY]       not run")

    # [SPARSITY]
    scaling = gap_results.get("scaling", {})
    if isinstance(scaling, dict) and "error" in scaling:
        print(f"[SPARSITY]        ERROR: {scaling['error']}")
    elif isinstance(scaling, dict):
        beta = scaling.get("nnz_beta")
        if beta is not None:
            tag = "  [WARNING: super-quadratic]" if beta > 2.0 else ""
            print(f"[SPARSITY]        NNZ growth exponent: {beta:.3f}{tag}")
        else:
            print("[SPARSITY]        NNZ growth fit not available (insufficient data)")
    else:
        print("[SPARSITY]        not run")

    # [UNITARITY]
    if unitarity_max_dev is not None:
        print(f"[UNITARITY]       max |<psi|psi>-1|: {unitarity_max_dev:.3e}")
    else:
        print("[UNITARITY]       (run deliverable 2 for unitarity test)")

    print("=" * 38)


def run_gaps(
    output_dir: str = "./figures",
    nmax_max: int = 8,
    unitarity_max_dev: Optional[float] = None,
) -> dict:
    """
    Run all three measurement gap diagnostics and print publication readiness summary.

    Gap 1: NVE control — VelocityVerlet on LiH without thermostat.
    Gap 2: Li energy audit — GeoVac vs exact at max_n=4, 6, 8.
    Gap 3: Two-panel scaling figure — NNZ growth exponent (sparsity audit).
    """
    from geovac.aimd import run_li_nve
    from geovac.benchmark import run_li_energy_audit, run_scaling_benchmark

    gap_results: dict = {}

    # --- Gap 1: NVE control ---
    print("\n\n[GAP 1 — NVE Control Run]\n")
    try:
        nve_res = run_li_nve(
            max_n=3, n_steps=600, dt=1.0,
            R_init=3.0, v_init=0.0, output_dir=output_dir,
        )
        gap_results["nve"] = nve_res
    except Exception as exc:
        print(f"\n  ERROR in Gap 1 (NVE): {exc}")
        gap_results["nve"] = {"error": str(exc)}

    # --- Gap 2: Li energy audit ---
    print("\n\n[GAP 2 — Li Energy Audit]\n")
    try:
        audit_rows = run_li_energy_audit(
            max_n_list=[4, 6, 8], output_dir=output_dir
        )
        gap_results["audit"] = audit_rows
    except Exception as exc:
        print(f"\n  ERROR in Gap 2 (energy audit): {exc}")
        gap_results["audit"] = {"error": str(exc)}

    # --- Gap 3: Two-panel scaling figure ---
    print("\n\n[GAP 3 — Sparsity / Scaling Figure]\n")
    try:
        nmax_range = list(range(3, min(nmax_max, 9) + 1))
        scale_res = run_scaling_benchmark(nmax_range=nmax_range, output_dir=output_dir)
        gap_results["scaling"] = scale_res
    except Exception as exc:
        print(f"\n  ERROR in Gap 3 (scaling): {exc}")
        gap_results["scaling"] = {"error": str(exc)}

    _print_pub_summary(gap_results, unitarity_max_dev)
    return gap_results


def main() -> None:
    parser = argparse.ArgumentParser(
        description="GeoVac 3-Electron Extension: Run All Deliverables"
    )
    parser.add_argument(
        "--quick", action="store_true",
        help="Reduced steps for testing (1000 CN steps, 60 AIMD steps)"
    )
    parser.add_argument(
        "--nmax", type=int, default=8,
        help="Max nmax for scaling benchmark (default 8)"
    )
    parser.add_argument(
        "--output-dir", type=str, default="./figures",
        help="Output directory for figures (default ./figures)"
    )
    parser.add_argument(
        "--skip", nargs="*", choices=["1", "2", "3"],
        help="Skip deliverables, e.g. --skip 1 3"
    )
    parser.add_argument(
        "--gaps-only", action="store_true",
        help="Run only the three gap diagnostics (NVE, Li energy audit, scaling) and exit"
    )
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # --- Gaps-only mode ---
    if args.gaps_only:
        print("=" * 65)
        print("GeoVac: Gap Diagnostics Only")
        print("=" * 65)
        print(f"  Output:  {os.path.abspath(args.output_dir)}")
        run_gaps(
            output_dir=args.output_dir,
            nmax_max=args.nmax,
            unitarity_max_dev=None,
        )
        return

    # --- Full run ---
    skip = set(args.skip or [])

    print("=" * 65)
    print("GeoVac: 3-Electron Sparse Dynamics Extension")
    print("=" * 65)
    print(f"  Output:  {os.path.abspath(args.output_dir)}")
    print(f"  Mode:    {'QUICK' if args.quick else 'FULL'}")
    if skip:
        print(f"  Skipping deliverables: {skip}")

    results_summary = {}
    unitarity_max_dev: Optional[float] = None

    if "1" not in skip:
        print("\n\n[DELIVERABLE 1 — LatticeIndex]\n")
        try:
            run_deliverable_1(max_n=4)
            results_summary["deliverable_1"] = "OK"
        except Exception as exc:
            print(f"\n  ERROR in Deliverable 1: {exc}")
            results_summary["deliverable_1"] = f"ERROR: {exc}"

    if "2" not in skip:
        print("\n\n[DELIVERABLE 2 — Benchmarks]\n")
        try:
            d2_res = run_deliverable_2(
                quick=args.quick,
                nmax_max=args.nmax,
                output_dir=args.output_dir,
            )
            unitarity_max_dev = d2_res.get("unitarity_max_dev")
            results_summary["deliverable_2"] = "OK"
        except Exception as exc:
            print(f"\n  ERROR in Deliverable 2: {exc}")
            results_summary["deliverable_2"] = f"ERROR: {exc}"

    if "3" not in skip:
        print("\n\n[DELIVERABLE 3 — AIMD]\n")
        try:
            run_deliverable_3(quick=args.quick, output_dir=args.output_dir)
            results_summary["deliverable_3"] = "OK"
        except Exception as exc:
            print(f"\n  ERROR in Deliverable 3: {exc}")
            results_summary["deliverable_3"] = f"ERROR: {exc}"

    # --- Gap diagnostics (appended to full run) ---
    print("\n\n[GAP DIAGNOSTICS]\n")
    try:
        run_gaps(
            output_dir=args.output_dir,
            nmax_max=args.nmax,
            unitarity_max_dev=unitarity_max_dev,
        )
        results_summary["gap_diagnostics"] = "OK"
    except Exception as exc:
        print(f"\n  ERROR in gap diagnostics: {exc}")
        results_summary["gap_diagnostics"] = f"ERROR: {exc}"

    # --- Final summary ---
    print("\n\n" + "=" * 65)
    print("SUMMARY")
    print("=" * 65)
    for key, status in results_summary.items():
        icon = "OK" if status == "OK" else "!!"
        print(f"  [{icon}] {key}: {status}")
    print(f"\n  Figures written to: {os.path.abspath(args.output_dir)}/")
    print("=" * 65)


if __name__ == "__main__":
    main()
