"""Yukawa-period PSLQ at MS-bar unification scale (~2e16 GeV).

Second-scale follow-on to debug/sprint_yukawa_pslq.py. Uses
SM-RG-evolved Yukawa values at the conventional GUT scale.

Reference values approximated from SM one-loop running with M_Z
inputs from Antusch-Maurer 2013 (arXiv:1306.6879) Table 1 and
Xing-Zhang-Zhou 2008 evolution. Precision is degraded by the running
itself; ~3-4 sig digits achievable for most flavours, ~2 digits for
light quarks (mu_s mass uncertainty dominant).

The framework's structural prediction (inner-factor restricted to
M1 union M2 by eta-trivialization) is scale-INDEPENDENT, so this is
a sanity check rather than an independent observable.
"""
from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

PROJ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJ))

import mpmath
from mpmath import mp, mpf, log, pi

from debug.sprint_yukawa_pslq import (
    m1_m2_basis, pslq_search, precision_verdict, TRANSFORMS,
)


# Yukawa values at scale mu = 2 x 10^16 GeV under SM RG running.
# Approximate values from Buttazzo et al. 2013 (arXiv:1307.3536) Fig. 1
# and Xing-Zhang-Zhou 2008 evolution.  Precision degraded by ~1-2 digits
# from M_Z due to RG-running uncertainty.

YUKAWA_INPUTS_UNIF = [
    # (name, y_central_str, y_uncertainty_str, sig_digits)
    ("e",   "2.794e-6",  "0.003e-6",   4),
    ("mu",  "5.900e-4",  "0.006e-4",   4),
    ("tau", "1.003e-2",  "0.012e-2",   4),
    ("u",   "2.99e-6",   "0.59e-6",    2),
    ("d",   "6.43e-6",   "0.59e-6",    2),
    ("s",   "1.28e-4",   "0.14e-4",    2),
    ("c",   "1.47e-3",   "0.07e-3",    3),
    ("b",   "6.71e-3",   "0.20e-3",    3),
    ("t",   "0.494",     "0.012",      3),
]


def compute_yukawa_values_unif(dps=60):
    mp.dps = dps
    out = {}
    for name, y_str, sig_str, digits in YUKAWA_INPUTS_UNIF:
        out[name] = {
            "y": mpf(y_str),
            "y_uncertainty": mpf(sig_str),
            "sig_digits": digits,
            "y_central_str": y_str,
            "y_uncertainty_str": sig_str,
        }
    return out


def run_sweep_one_transform(transform_label, transform_fn, transform_desc,
                            yukawa_data, basis, max_coeffs, dps=120):
    mp.dps = dps
    results = {
        "transform": transform_label,
        "transform_description": transform_desc,
        "basis_labels": [lbl for lbl, _ in basis],
        "max_coeffs_tested": list(max_coeffs),
        "per_fermion": {},
    }
    any_hit = False
    for name, info in yukawa_data.items():
        y = info["y"]
        try:
            t_val = transform_fn(y)
        except Exception as e:
            results["per_fermion"][name] = {"error": str(e)}
            continue
        cells = pslq_search(t_val, basis, max_coeffs=max_coeffs, dps=dps)
        for c in cells:
            c["precision_verdict"] = precision_verdict(info["sig_digits"], c["max_coeff"], len(basis))
        honest_hit = any(c.get("found") and c["precision_verdict"]["honest"] for c in cells)
        if honest_hit:
            any_hit = True
        results["per_fermion"][name] = {
            "y_value": str(y),
            "transform_value": str(t_val),
            "sig_digits": info["sig_digits"],
            "cells": cells,
            "any_honest_hit": honest_hit,
        }
    results["any_honest_hit_across_fermions"] = any_hit
    return results


def main(dps=120):
    mp.dps = dps
    timestamp = datetime.now(timezone.utc).isoformat()
    yukawa = compute_yukawa_values_unif(dps=dps)
    basis = m1_m2_basis(dps=dps)
    max_coeffs = (10, 100, 1000)

    print(f"Yukawa-period PSLQ at MS-bar unification scale ({timestamp})")
    print(f"basis size = {len(basis)} (M1 union M2 pure-Tate)")
    print(f"max_coeffs = {max_coeffs}")
    print()

    all_results = {
        "timestamp": timestamp,
        "scale": "MS-bar mu = 2e16 GeV (GUT scale)",
        "rg_source": "Buttazzo+2013 / Xing-Zhang-Zhou 2008 SM running",
        "basis": "M1 union M2 = {1, pi, pi^2, pi^4, pi^6, pi^8, 1/pi, 1/pi^2, 1/pi^4}",
        "transforms": [t[0] for t in TRANSFORMS],
        "max_coeffs": list(max_coeffs),
        "yukawa_inputs": YUKAWA_INPUTS_UNIF,
        "yukawa_values": {k: {
            "y": str(v["y"]),
            "y_uncertainty": str(v["y_uncertainty"]),
            "sig_digits": v["sig_digits"],
            "y_central_str": v["y_central_str"],
            "y_uncertainty_str": v["y_uncertainty_str"],
        } for k, v in yukawa.items()},
        "per_transform": {},
        "halt_history": [],
    }

    for tlabel, tfn, tdesc in TRANSFORMS:
        print(f"--- Transform: {tlabel} ({tdesc}) ---")
        result = run_sweep_one_transform(
            tlabel, tfn, tdesc, yukawa, basis, max_coeffs, dps=dps,
        )
        all_results["per_transform"][tlabel] = result
        any_hit = result["any_honest_hit_across_fermions"]
        for fname, fdata in result["per_fermion"].items():
            cells = fdata.get("cells", [])
            for c in cells:
                hit = "HIT" if c.get("found") else "no hit"
                honest = "(honest)" if c["precision_verdict"]["honest"] else "(precision-shy)"
                print(f"  {fname:4s} M={c['max_coeff']:5d}  {hit:7s} {honest:16s}"
                      f"  sig_dig={fdata['sig_digits']} need={c['precision_verdict']['digits_needed']}")
            print()
        if any_hit:
            all_results["halt_history"].append({
                "transform": tlabel,
                "reason": "honest_hit_found",
                "next_step": "curve_fit_audit",
            })
            print(f">>> HALT: honest hit in transform {tlabel}.")
            break
        else:
            all_results["halt_history"].append({
                "transform": tlabel,
                "reason": "clean_null",
                "next_step": "proceed",
            })
            print(f"<<< CLEAN NULL across all fermions for transform {tlabel}.")

    out_path = PROJ / "debug" / "data" / "sprint_yukawa_pslq_unif.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\nWrote results to {out_path}")
    return all_results


if __name__ == "__main__":
    main(dps=120)
