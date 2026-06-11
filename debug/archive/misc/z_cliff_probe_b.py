"""
Z>20 cliff diagnostic — Probe (b): FD vs analytical hydrogenic solver consistency.

Hypothesis (b): The dense-uniform-FD path in _solve_screened_radial_log
(used by screened_psi_origin_squared) may diverge or converge poorly for
s-wave at l=0, especially at heavy-atom Z. The closeout sprint flagged
the original log-mesh having artifacts; the dense uniform grid was the
sprint's replacement. Test whether THAT replacement is faithful for the
KNOWN-answer hydrogenic case at varied Z.

For hydrogenic |psi_{ns}(0)|^2 = Z^3 / (pi n^3) is the analytical result.
Test at Z=1, 2, 6, 10, 30, 55, 80 for n=1, 2, 3, 6.

If FD converges cleanly at all Z (matching analytical), the closeout
sprint's solver is faithful and FD is NOT the dominant cliff cause.
If it breaks at high Z, that's a real engineering fix path.

Read-only: uses geovac.neon_core._solve_screened_radial_log directly.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from geovac.neon_core import _solve_screened_radial_log


def hydrogenic_psi0_sq_analytical(Z: float, n: int) -> float:
    """|psi_{ns}(0)|^2 = Z^3 / (pi n^3) (atomic units)."""
    return Z**3 / (np.pi * n**3)


def hydrogenic_E_n_analytical(Z: float, n: int) -> float:
    """E_n = -Z^2 / (2 n^2) (atomic units)."""
    return -Z**2 / (2.0 * n**2)


def bare_z_eff_callable(Z: float):
    """Constant Z_eff(r) = Z (no screening, hydrogenic)."""
    def _z(r):
        return np.full_like(np.atleast_1d(np.asarray(r, dtype=float)),
                            float(Z))
    return _z


def main():
    test_cases = [
        # (Z, n, label)
        (1.0, 1, "H 1s"),
        (1.0, 2, "H 2s"),
        (1.0, 6, "H 6s"),
        (2.0, 1, "He+ 1s"),
        (2.0, 6, "He+ 6s"),
        (6.0, 6, "C5+ 6s"),
        (10.0, 6, "Ne9+ 6s"),
        (30.0, 6, "Zn29+ 6s"),
        (55.0, 6, "Cs54+ 6s"),    # the case that matters for Cs HFS
        (80.0, 6, "Hg79+ 6s"),    # extreme heavy
    ]

    n_grid_panel = [50_000, 100_000, 200_000, 400_000]

    # Auto-scale r_max for the requested state
    def auto_r_max(Z, n):
        # Hydrogenic <r> ~ 1.5 n^2 / Z; r_max should be ~5-10x <r>
        rmax = max(20.0, 15.0 * n * n / Z)
        return float(rmax)

    rows = []

    print("=" * 110)
    print("Z>20 CLIFF DIAGNOSTIC — PROBE (b): FD solver hydrogenic faithfulness")
    print("=" * 110)
    print()
    print(f"{'system':12s}  {'n_grid':>8s}  {'r_max':>7s}  "
          f"{'E_FD':>14s}  {'E_anal':>14s}  {'rel_err_E':>11s}  "
          f"{'psi0_FD':>14s}  {'psi0_anal':>14s}  {'rel_err_psi0':>13s}  "
          f"{'wall_s':>8s}")
    print("-" * 110)

    for Z, n, label in test_cases:
        psi0_anal = hydrogenic_psi0_sq_analytical(Z, n)
        E_anal = hydrogenic_E_n_analytical(Z, n)
        rmax = auto_r_max(Z, n)

        for n_grid in n_grid_panel:
            t0 = time.perf_counter()
            try:
                E_fd, _u, _r, R0 = _solve_screened_radial_log(
                    Z=int(Z) if Z == int(Z) else 1, l=0, n_target=n,
                    z_eff_callable=bare_z_eff_callable(Z),
                    Z_origin=Z,
                    n_grid=n_grid,
                    r_max=rmax,
                )
                dt = time.perf_counter() - t0
                psi0_fd = R0**2 / (4.0 * np.pi)
            except Exception as e:
                dt = time.perf_counter() - t0
                E_fd = float('nan')
                psi0_fd = float('nan')
                err_msg = str(e)[:30]
                print(f"  {label:12s}  {n_grid:8d}  {rmax:7.2f}  "
                      f"  ERROR: {err_msg}")
                continue

            rel_err_E = 100.0 * (E_fd - E_anal) / E_anal
            rel_err_psi0 = 100.0 * (psi0_fd - psi0_anal) / psi0_anal

            print(f"{label:12s}  {n_grid:8d}  {rmax:7.2f}  "
                  f"{E_fd:14.6f}  {E_anal:14.6f}  {rel_err_E:11.3f}  "
                  f"{psi0_fd:14.4e}  {psi0_anal:14.4e}  "
                  f"{rel_err_psi0:13.2f}  {dt:8.2f}")

            rows.append({
                "Z": Z, "n": n, "label": label,
                "n_grid": n_grid, "r_max": rmax,
                "E_FD": float(E_fd), "E_analytical": E_anal,
                "rel_err_E_pct": float(rel_err_E),
                "psi0_FD": float(psi0_fd), "psi0_analytical": float(psi0_anal),
                "rel_err_psi0_pct": float(rel_err_psi0),
                "wall_s": float(dt),
            })

        print()

    # Convergence diagnosis: at the finest grid, what is the error?
    print("Finest-grid (n_grid=400k) summary:")
    print(f"{'system':12s}  {'rel_err_E%':>11s}  {'rel_err_psi0%':>14s}  {'verdict':30s}")
    print("-" * 70)

    final_rows = [r for r in rows if r["n_grid"] == 400_000]
    for r in final_rows:
        e_err = r["rel_err_E_pct"]
        p_err = r["rel_err_psi0_pct"]
        if abs(e_err) < 1.0 and abs(p_err) < 5.0:
            verdict = "FAITHFUL (sub-percent)"
        elif abs(e_err) < 5.0 and abs(p_err) < 20.0:
            verdict = "ACCEPTABLE"
        else:
            verdict = "FD BREAKS at this Z"
        print(f"{r['label']:12s}  {e_err:+11.3f}  {p_err:+14.2f}  {verdict:30s}")

    print()

    # Convergence rate from n_grid=50k -> 400k
    print("Convergence rate per system (psi0 error change with grid):")
    print(f"{'system':12s}  {'err@50k':>12s}  {'err@100k':>12s}  "
          f"{'err@200k':>12s}  {'err@400k':>12s}  {'monotone?':>9s}")
    print("-" * 80)

    for Z, n, label in test_cases:
        sys_rows = [r for r in rows if r["Z"] == Z and r["n"] == n]
        if len(sys_rows) < len(n_grid_panel):
            continue
        errs = [r["rel_err_psi0_pct"] for r in sorted(sys_rows, key=lambda r: r["n_grid"])]
        # Monotone toward zero?
        monotone = all(abs(errs[i+1]) <= abs(errs[i]) + 0.1 for i in range(len(errs)-1))
        print(f"{label:12s}  " + "  ".join(f"{e:+12.2f}" for e in errs)
              + f"  {'YES' if monotone else 'NO':>9s}")

    # Summary verdict
    cs_rows = [r for r in final_rows if abs(r["Z"] - 55) < 0.1]
    cs_at_55 = cs_rows[0] if cs_rows else None
    h_at_1 = next((r for r in final_rows if r["Z"] == 1.0 and r["n"] == 1), None)

    print()
    print("Z-scaling analysis:")
    if h_at_1 and cs_at_55:
        print(f"  H 1s  at n_grid=400k: rel_err(E)={h_at_1['rel_err_E_pct']:+.3f}%, "
              f"rel_err(psi0)={h_at_1['rel_err_psi0_pct']:+.2f}%")
        print(f"  Cs54+ 6s at n_grid=400k: rel_err(E)={cs_at_55['rel_err_E_pct']:+.3f}%, "
              f"rel_err(psi0)={cs_at_55['rel_err_psi0_pct']:+.2f}%")
        if abs(h_at_1["rel_err_psi0_pct"]) < 5.0 and abs(cs_at_55["rel_err_psi0_pct"]) < 5.0:
            verdict = ("FD SOLVER IS FAITHFUL across Z range tested. "
                       "FD convergence is NOT the dominant cliff mechanism.")
        elif abs(cs_at_55["rel_err_psi0_pct"]) > 3.0 * abs(h_at_1["rel_err_psi0_pct"]):
            verdict = ("FD SOLVER PARTIAL FAILURE at high Z. "
                       "FD error grows ~Z; needs fix for Cs.")
        else:
            verdict = ("FD SOLVER ACCEPTABLE; mild Z dependence in error.")
    else:
        verdict = "Cs/H test points missing"

    print()
    print(f"VERDICT (Probe b): {verdict}")

    out = {
        "probe": "b",
        "hypothesis": "FD solver in _solve_screened_radial_log breaks at heavy-atom Z",
        "rows": rows,
        "verdict": verdict,
    }
    out_path = Path(ROOT) / "debug" / "data" / "z_cliff_probe_b.json"
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    main()
