"""Compute the Dirac-S^3 Ihara zeta data for Track RH-C.

Produces factored closed forms, Ramanujan verdicts, and JSON-serializable
payloads for debug/data/ihara_zeta_dirac_s3.json.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import sympy as sp

sys.stdout.reconfigure(encoding="utf-8")

from geovac.ihara_zeta_dirac import (
    dirac_s3_ihara_zeta,
    describe_adjacency_rule,
    build_dirac_s3_graph,
)


def _json_safe(v):
    """Convert a single value to a JSON-serializable primitive."""
    if isinstance(v, (bool, int, str)) or v is None:
        return v
    if isinstance(v, float):
        return v if np.isfinite(v) else None
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (np.floating,)):
        return float(v) if np.isfinite(v) else None
    if isinstance(v, np.ndarray):
        return [_json_safe(x) for x in v.tolist()]
    if isinstance(v, complex):
        return {"re": float(v.real), "im": float(v.imag)}
    if isinstance(v, (np.complex128, np.complex64)):
        return {"re": float(v.real), "im": float(v.imag)}
    if isinstance(v, list):
        return [_json_safe(x) for x in v]
    if isinstance(v, tuple):
        return [_json_safe(x) for x in v]
    if isinstance(v, dict):
        return {str(k): _json_safe(x) for k, x in v.items()}
    if isinstance(v, sp.Basic):
        return str(v)
    return str(v)


all_rows = []

print("=" * 78)
print("Track RH-C: Dirac-S^3 Ihara zeta — full scan")
print("=" * 78)
print()
print("Rule A:", describe_adjacency_rule("A"))
print()
print("Rule B:", describe_adjacency_rule("B"))
print()
print("-" * 78)

for n_max in [1, 2, 3]:
    for rule in ["A", "B"]:
        try:
            r = dirac_s3_ihara_zeta(n_max=n_max, adjacency_rule=rule, factor=True)
        except Exception as e:
            print(f"n_max={n_max} rule={rule}: ERROR {type(e).__name__}: {e}")
            all_rows.append({
                "n_max": n_max,
                "adjacency_rule": rule,
                "error": f"{type(e).__name__}: {e}",
            })
            continue

        print(f"\nn_max={n_max} rule={rule}")
        print(
            f"  V={r['V']} E={r['E']} c={r['c']} beta1={r['r_betti1']} "
            f"q_max={r['q_max']} max_deg={max(r['degree_sequence']) if r['degree_sequence'] else 0}"
        )
        print(
            f"  rho(T)={r['spectral_radius']:.4f} "
            f"max|mu_nt|={r['max_abs_nontrivial']:.4f} "
            f"sqrt(q_max)={r['sqrt_q_max']:.4f} "
            f"deviation={r['ramanujan_deviation']:+.4f}"
        )
        print(f"  Ramanujan: {r['ramanujan_verdict']}")
        if r["zeta_inverse_factored"] is not None:
            factored = r["zeta_inverse_factored"]
            print(f"  zeta_G(s)^(-1) = {factored}")
        print(f"  per-kappa: {r['per_kappa_sizes']}")

        # Data payload for JSON
        zeros_arr = r["zeros"] if len(r["zeros"]) else np.array([], dtype=complex)
        row = {
            "n_max": n_max,
            "adjacency_rule": rule,
            "rule_description": describe_adjacency_rule(rule),
            "V": r["V"],
            "E": r["E"],
            "c": r["c"],
            "beta1": r["r_betti1"],
            "degree_sequence": r["degree_sequence"],
            "q_max": r["q_max"],
            "sqrt_q_max": r["sqrt_q_max"],
            "spectral_radius_rho": r["spectral_radius"],
            "max_abs_nontrivial": r["max_abs_nontrivial"],
            "ramanujan_verdict": r["ramanujan_verdict"],
            "ramanujan_deviation": r["ramanujan_deviation"],
            "zeta_inv_expanded": str(r["zeta_inverse_expanded"]),
            "zeta_inv_factored": str(r["zeta_inverse_factored"]),
            "zeros": [_json_safe(z) for z in zeros_arr.tolist()],
            "per_kappa_sizes": {str(k): v for k, v in r["per_kappa_sizes"].items()},
            "labels_triple": r["labels_triple"],
        }
        all_rows.append(row)

out_path = Path("debug/data/ihara_zeta_dirac_s3.json")
out_path.parent.mkdir(exist_ok=True, parents=True)
with out_path.open("w", encoding="utf-8") as fp:
    json.dump(
        {
            "generated_by": "debug/compute_ihara_zeta_dirac.py",
            "sprint": "hey-buddy-we-need-crystalline-sprout (Track RH-C)",
            "date": "2026-04-17",
            "rows": all_rows,
        },
        fp,
        indent=2,
        default=_json_safe,
    )
print("\n" + "=" * 78)
print(f"Wrote {out_path}")
