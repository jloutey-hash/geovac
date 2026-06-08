"""Verification panel for the S^2 sector restriction sprint.

For each molecule in the panel, build the untapered GeoVac composed
Hamiltonian, infer ``M``, and report:

  - full ``S_z = 0`` FCI dimension
  - singlet (``S = 0``) sector dimension
  - ratio (VQE speedup factor)
  - commutator residual ``||[S^2, H]||_1`` (should be ~0 to machine
    precision on the untapered builder).

Run from the worktree root:

    python debug/sprint_s2_sector_panel.py
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List

from geovac.ecosystem_export import hamiltonian
from geovac.spin_sector import (
    singlet_sector_dim,
    sz_zero_sector_dim,
    verify_s2_commutes,
)


PANEL: List[Dict[str, Any]] = [
    # (canonical_name, n_electrons)
    {"system": "H2", "n_electrons": 2},
    {"system": "LiH", "n_electrons": 4},
    {"system": "BeH2", "n_electrons": 6},
    {"system": "H2O", "n_electrons": 10},
    {"system": "NH3", "n_electrons": 10},
    {"system": "CH4", "n_electrons": 10},
    {"system": "HF", "n_electrons": 10},
]


def _qubit_count(op) -> int:
    max_idx = -1
    for term in op.terms.keys():
        for q, _p in term:
            if q > max_idx:
                max_idx = q
    return max_idx + 1


def run_panel(max_n: int = 2) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for entry in PANEL:
        sys_name = entry["system"]
        n_el = entry["n_electrons"]
        try:
            ham = hamiltonian(sys_name, max_n=max_n)
        except Exception as exc:
            rows.append({
                "system": sys_name,
                "n_electrons": n_el,
                "status": f"BUILD_ERROR: {type(exc).__name__}: {exc}",
            })
            continue

        op = ham.to_openfermion()
        Q = _qubit_count(op)
        M = Q // 2

        d_sz0 = sz_zero_sector_dim(M, n_el)
        d_sing = singlet_sector_dim(M, n_el)
        ratio = d_sz0 / d_sing if d_sing > 0 else float("inf")

        commutes, residual = verify_s2_commutes(op, M, atol=1e-10)

        rows.append({
            "system": sys_name,
            "n_electrons": n_el,
            "Q": Q,
            "M": M,
            "dim_Sz0": d_sz0,
            "dim_singlet": d_sing,
            "speedup_factor": round(ratio, 3),
            "commutator_residual": float(residual),
            "commutes": bool(commutes),
            "status": "ok",
        })

    return rows


def main() -> None:
    rows = run_panel(max_n=2)

    print()
    print("S^2 sector verification panel (max_n=2, untapered builder)")
    print("-" * 96)
    hdr = (
        f"{'System':<8} {'Q':>3} {'M':>3} {'N_e':>4} "
        f"{'dim(S_z=0)':>12} {'dim(S=0)':>12} {'speedup':>8} "
        f"{'commutes':>9} {'residual':>12}"
    )
    print(hdr)
    print("-" * 96)
    for r in rows:
        if r.get("status") != "ok":
            print(f"{r['system']:<8} {r['status']}")
            continue
        print(
            f"{r['system']:<8} {r['Q']:>3} {r['M']:>3} {r['n_electrons']:>4} "
            f"{r['dim_Sz0']:>12d} {r['dim_singlet']:>12d} "
            f"{r['speedup_factor']:>8.3f} "
            f"{str(r['commutes']):>9} {r['commutator_residual']:>12.3e}"
        )
    print("-" * 96)

    # JSON dump for downstream consumers (memo, paper edits).
    here = Path(__file__).parent if "__file__" in globals() else Path("debug")
    out = here / "data" / "sprint_s2_sector_panel.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rows, indent=2))
    print(f"\nWrote JSON panel to {out}")


if __name__ == "__main__":
    main()
